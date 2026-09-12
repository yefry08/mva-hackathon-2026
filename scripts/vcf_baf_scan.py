#!/usr/bin/env python3
"""Barrido de BAF por ventanas para alteraciones cromosomicas en mosaico.

Reemplaza el metodo de lobulos de vcf_aneuploidy.py, que estaba mal: partir en
0.5 una distribucion unimodal siempre devuelve dos lobulos, asi que ese script
reportaba fracciones celulares identicas en los 24 cromosomas. Eso era el
metodo, no el paciente.

La prueba correcta. Para un heterocigoto disomico con profundidad DP, el conteo
del alelo alterno es binomial(DP, 0.5), asi que la desviacion esperada de BAF
respecto de 0.5 es aproximadamente 0.3989/sqrt(DP). Un evento en mosaico corre
el BAF de todos los hets de la region, de modo que la desviacion observada supera
sistematicamente a la esperada. La razon observado/esperado es adimensional y
comparable entre ventanas con distinta profundidad, que es lo que el enfoque
anterior no lograba.

Un evento real tiene que cumplir las dos cosas: exceso de dispersion del BAF Y
desviacion de cobertura en la misma direccion (ganancia arriba, perdida abajo).
Solo exceso de BAF con cobertura plana tambien lo produce una disomia
uniparental en mosaico, que es un hallazgo distinto y se marca aparte.

Uso: python scripts/vcf_baf_scan.py [ruta.vcf.gz]
"""
from __future__ import annotations

import gzip
import json
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
def _default_vcf():
    """Descubre el VCF en data/ sin escribir el nombre de la muestra en el codigo."""
    d = ROOT / "data"
    hits = sorted(p for p in d.glob("*.vcf.gz")) if d.exists() else []
    return hits[0] if hits else d / "input.vcf.gz"


DEFAULT_VCF = _default_vcf()
OUT = ROOT / "work" / "baf_scan.json"

WINDOW = 10_000_000      # 10 Mb: resolucion util para eventos en mosaico
MIN_DP = 20
MIN_HET = 150            # menos hets que esto y la ventana no dice nada
EXCESS_FLAG = 1.25       # 25% mas dispersion que la mediana autosomica
COV_FLAG = 0.08          # 8% de desviacion de cobertura
MIN_RUN = 3              # un evento real ocupa varias ventanas seguidas, no una
AUTOSOME = re.compile(r"^(chr)?([1-9]|1\d|2[0-2])$")
CANON = re.compile(r"^(chr)?([1-9]|1\d|2[0-2]|X|Y)$")
# E|Z| de una normal estandar es sqrt(2/pi) = 0.7979. Ojo: 0.3989 es la densidad
# en cero, no la esperanza del valor absoluto; confundirlas duplica el exceso
# aparente en todas las ventanas por igual, que es como se detecto el error.
E_ABS_DEV = math.sqrt(2.0 / math.pi)


def main() -> int:
    vcf = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_VCF
    if not vcf.exists():
        print("no existe %s" % vcf)
        return 1

    win = defaultdict(lambda: {"het": 0, "obs": 0.0, "exp": 0.0, "dp": 0, "n": 0})

    with gzip.open(vcf, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            chrom, pos, ref, alt, filt, fmt, sample = f[0], f[1], f[3], f[4], f[6], f[8], f[9]
            if not CANON.match(chrom) or filt not in ("PASS", "."):
                continue
            if len(ref) != 1 or len(alt) != 1 or alt == ".":
                continue
            kv = dict(zip(fmt.split(":"), sample.split(":")))
            gt = kv.get("GT", "./.").replace("|", "/")
            parts = kv.get("AD", "").split(",")
            if len(parts) < 2 or not parts[0].isdigit() or not parts[1].isdigit():
                continue
            r, v = int(parts[0]), int(parts[1])
            dp = r + v
            if dp < MIN_DP:
                continue
            key = (chrom, int(pos) // WINDOW)
            w = win[key]
            w["dp"] += dp
            w["n"] += 1
            a = gt.split("/")
            if len(a) == 2 and "." not in a and a[0] != a[1]:
                w["het"] += 1
                w["obs"] += abs(v / dp - 0.5)
                w["exp"] += E_ABS_DEV * 0.5 / math.sqrt(dp)

    # Linea de base de cobertura: mediana de las ventanas autosomicas.
    auto_dp = sorted(w["dp"] / w["n"] for k, w in win.items()
                     if AUTOSOME.match(k[0]) and w["n"] >= 500)
    base_dp = auto_dp[len(auto_dp) // 2] if auto_dp else float("nan")

    rows, flagged = [], []
    for (chrom, idx), w in sorted(win.items(),
                                  key=lambda kv: ((0, int(kv[0][0].replace("chr", "")))
                                                  if kv[0][0].replace("chr", "").isdigit()
                                                  else (1, kv[0][0]), kv[0][1])):
        if w["het"] < MIN_HET or not w["exp"]:
            continue
        excess = (w["obs"] / w["het"]) / (w["exp"] / w["het"])
        cov = (w["dp"] / w["n"]) / base_dp
        row = {"chrom": chrom, "start": idx * WINDOW, "het": w["het"],
               "baf_excess": excess, "cov_ratio": cov}
        rows.append(row)

    auto_rows = [r for r in rows if AUTOSOME.match(r["chrom"])]
    # Normalizamos contra la mediana autosomica: el sesgo de referencia y los
    # errores de mapeo suben el exceso por igual en todo el genoma, y lo que
    # interesa es la desviacion local respecto de ese piso.
    med_excess = sorted(r["baf_excess"] for r in auto_rows)
    med_excess = med_excess[len(med_excess) // 2] if med_excess else 1.0
    for r in rows:
        r["baf_excess_norm"] = r["baf_excess"] / med_excess
        r["call"] = None
        if AUTOSOME.match(r["chrom"]) and r["baf_excess_norm"] >= EXCESS_FLAG:
            if r["cov_ratio"] >= 1 + COV_FLAG:
                r["call"] = "ganancia en mosaico"
            elif r["cov_ratio"] <= 1 - COV_FLAG:
                r["call"] = "perdida en mosaico"
            else:
                r["call"] = "exceso de BAF con cobertura plana (UPD en mosaico?)"

    # Un evento en mosaico ocupa un segmento cromosomico, no una ventana suelta.
    # Los picos de una sola ventana viven en centromeros y brazos acrocentricos,
    # donde el mapeo se rompe. Exigimos una corrida de MIN_RUN ventanas seguidas
    # con la misma lectura.
    by_chrom = defaultdict(list)
    for r in rows:
        by_chrom[r["chrom"]].append(r)
    for chrom, rs in by_chrom.items():
        rs.sort(key=lambda r: r["start"])
        run = []
        for r in rs + [None]:
            same = (r is not None and r["call"] is not None and run
                    and r["call"] == run[-1]["call"]
                    and r["start"] == run[-1]["start"] + WINDOW)
            if r is not None and r["call"] is not None and (not run or same):
                run.append(r)
                continue
            if len(run) >= MIN_RUN:
                flagged.append({"chrom": chrom, "start": run[0]["start"],
                                "end": run[-1]["start"] + WINDOW,
                                "windows": len(run),
                                "het": sum(x["het"] for x in run),
                                "baf_excess": sum(x["baf_excess_norm"] for x in run) / len(run),
                                "cov_ratio": sum(x["cov_ratio"] for x in run) / len(run),
                                "call": run[0]["call"]})
            run = [r] if (r is not None and r["call"] is not None) else []
    ex = sorted(r["baf_excess"] for r in auto_rows)
    print("ventanas de %d Mb con >=%d hets: %d autosomicas" % (WINDOW // 10**6, MIN_HET, len(auto_rows)))
    print("cobertura de base: %.1fx" % base_dp)
    print("exceso de BAF, distribucion autosomica:")
    for label, q in (("p05", 0.05), ("mediana", 0.5), ("p95", 0.95), ("max", 1.0)):
        if ex:
            print("  %-8s %.3f" % (label, ex[min(len(ex) - 1, int(q * (len(ex) - 1)))]))
    print("\nreferencia: 1.00 es exactamente lo que predice la binomial. Por sesgo de")
    print("referencia y errores de mapeo lo tipico queda algo por encima de 1.")

    singles = [r for r in rows if r.get("call") and AUTOSOME.match(r["chrom"])]
    print("\nventanas sueltas con senal: %d (la mayoria son centromeros y brazos"
          " acrocentricos)" % len(singles))

    if flagged:
        print("\nsegmentos de >=%d ventanas seguidas (%d):" % (MIN_RUN, len(flagged)))
        print("%-6s %12s %12s %5s %8s %9s %8s  %s" %
              ("chrom", "inicio", "fin", "vent", "hets", "excesoBAF", "cobert", "lectura"))
        for r in flagged:
            print("%-6s %12d %12d %5d %8d %9.2f %8.3f  %s" %
                  (r["chrom"], r["start"], r["end"], r["windows"], r["het"],
                   r["baf_excess"], r["cov_ratio"], r["call"]))
    else:
        print("\nningun segmento de >=%d ventanas consecutivas supera los umbrales." % MIN_RUN)
        print("Las %d ventanas sueltas de arriba son picos aislados; los mas altos caen en" % len(singles))
        print("brazos acrocentricos, centromeros y bloques de heterocromatina, que es donde")
        print("el mapeo se rompe. Un evento en mosaico real abarca un segmento cromosomico.")
        print("\nCon cobertura de %.0fx y ~%d hets por ventana, este barrido detecta eventos" % (base_dp, MIN_HET))
        print("por encima de una fraccion celular de aproximadamente 10-15 por ciento. Un")
        print("mosaicismo por debajo de eso, o presente en otro tejido, no aparece aqui:")
        print("es un limite de deteccion, no una ausencia de aneuploidia.")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({"window": WINDOW, "baseline_dp": base_dp,
                               "flagged": flagged, "windows": rows}, indent=2),
                   encoding="utf-8")
    print("\nagregados en %s" % OUT.relative_to(ROOT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
