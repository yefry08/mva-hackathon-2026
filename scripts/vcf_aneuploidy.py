#!/usr/bin/env python3
"""SUPERSEDIDO por scripts/vcf_baf_scan.py. No usar sus columnas de lobulos.

El metodo de lobulos de aqui esta mal: parte la distribucion de BAF en 0.5 y toma
la mediana de cada mitad, lo que en una distribucion unimodal siempre devuelve dos
"lobulos" y, por lo tanto, una fraccion celular aparente. Corrido sobre este caso
devolvio 0.22/0.18 identicos en los 24 cromosomas, que es la firma de un artefacto
de metodo. Las columnas de profundidad de este script si sirven y se conservan.

Perfil de aneuploidia en mosaico desde el VCF entregado. Tarea G22 del plan.

Solo agregados por cromosoma: mediana de profundidad, forma de la distribucion
del balance alelico (BAF) de los heterocigotos, y la fraccion celular implicada
si la desviacion viene de una ganancia o una perdida en mosaico.

Modelo. En un cromosoma disomico los hets se agrupan en BAF 0.5. Con una trisomia
presente en una fraccion f de las celulas, los hets se parten en dos lobulos en
(1+f)/(2+f) y 1/(2+f), y la profundidad sube por (2+f)/2. Con una monosomia en
fraccion f, los lobulos van a 1/(2-f) y (1-f)/(2-f), y la profundidad baja.
De la separacion entre lobulos se despeja f, y la profundidad dice si fue
ganancia o perdida. Que ambas coincidan es lo que separa la senal del artefacto.

Filtros: solo PASS, solo SNV bialelicos, DP>=20. Se excluyen contigs no canonicos,
que es donde viven los artefactos de mapeo.

Uso: python scripts/vcf_aneuploidy.py [ruta.vcf.gz]
"""
from __future__ import annotations

import gzip
import json
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
OUT = ROOT / "work" / "aneuploidy.json"

MIN_DP = 20
MAX_DP = 300
NBINS = 50
CANON = re.compile(r"^(chr)?([1-9]|1\d|2[0-2]|X|Y)$")
AUTOSOME = re.compile(r"^(chr)?([1-9]|1\d|2[0-2])$")


def median_from_hist(hist):
    total = sum(hist)
    if not total:
        return float("nan")
    half, acc = total / 2.0, 0
    for i, n in enumerate(hist):
        acc += n
        if acc >= half:
            return float(i)
    return float(len(hist) - 1)


def lobes(ab_hist):
    """Mediana de cada lado de 0.5. Con lobulos simetricos, su separacion da f."""
    left = [(i + 0.5) / NBINS for i in range(NBINS // 2)]
    right = [(i + 0.5) / NBINS for i in range(NBINS // 2, NBINS)]
    lw = ab_hist[: NBINS // 2]
    rw = ab_hist[NBINS // 2:]

    def wmedian(vals, weights):
        total = sum(weights)
        if not total:
            return float("nan")
        half, acc = total / 2.0, 0
        for v, w in zip(vals, weights):
            acc += w
            if acc >= half:
                return v
        return vals[-1]

    return wmedian(left, lw), wmedian(right, rw)


def frac_from_lobes(lo, hi, mode):
    """Despeja la fraccion celular f a partir del lobulo superior."""
    if hi != hi or hi <= 0.5:
        return float("nan")
    if mode == "gain":          # hi = (1+f)/(2+f)
        denom = 1.0 - hi
        if denom <= 0:
            return float("nan")
        f = (2.0 * hi - 1.0) / denom
    else:                        # perdida: hi = 1/(2-f)
        f = 2.0 - 1.0 / hi
    return f if 0.0 <= f <= 1.0 else float("nan")


def main() -> int:
    vcf = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_VCF
    if not vcf.exists():
        print("no existe %s" % vcf)
        return 1

    dp_hist = defaultdict(lambda: [0] * (MAX_DP + 1))
    ab_hist = defaultdict(lambda: [0] * NBINS)
    n_het = defaultdict(int)
    n_var = defaultdict(int)
    n_pgt = 0
    n_pid_groups = set()

    with gzip.open(vcf, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            chrom, ref, alt, filt, fmt, sample = f[0], f[3], f[4], f[6], f[8], f[9]
            if not CANON.match(chrom):
                continue
            if filt not in ("PASS", "."):
                continue
            if len(ref) != 1 or len(alt) != 1 or alt == ".":
                continue

            keys = fmt.split(":")
            vals = sample.split(":")
            kv = dict(zip(keys, vals))
            gt = kv.get("GT", "./.")
            if "PGT" in kv and kv["PGT"] not in (".", ""):
                n_pgt += 1
                pid = kv.get("PID", "")
                if pid and pid != ".":
                    n_pid_groups.add((chrom, pid))

            a = gt.replace("|", "/").split("/")
            if len(a) != 2 or "." in a:
                continue
            n_var[chrom] += 1

            ad = kv.get("AD", "")
            parts = ad.split(",")
            if len(parts) < 2 or not parts[0].isdigit() or not parts[1].isdigit():
                continue
            r, v = int(parts[0]), int(parts[1])
            dp = r + v
            if dp < MIN_DP:
                continue
            dp_hist[chrom][min(dp, MAX_DP)] += 1
            if a[0] != a[1]:
                n_het[chrom] += 1
                ab_hist[chrom][min(NBINS - 1, int((v / dp) * NBINS))] += 1

    chroms = sorted(dp_hist, key=lambda c: (0, int(c.replace("chr", "")))
                    if c.replace("chr", "").isdigit() else (1, c))
    auto = [c for c in chroms if AUTOSOME.match(c)]
    auto_dp = sorted(median_from_hist(dp_hist[c]) for c in auto)
    baseline = auto_dp[len(auto_dp) // 2] if auto_dp else float("nan")

    print("profundidad autosomica de referencia: %.1fx" % baseline)
    print("variantes con fase PGT: %d en %d bloques PID\n" % (n_pgt, len(n_pid_groups)))
    print("%-6s %9s %8s %7s %7s %7s %8s %8s  %s" %
          ("chrom", "SNV PASS", "medDP", "razon", "lob.inf", "lob.sup", "f ganan", "f perd", "lectura"))

    rows = {}
    for c in chroms:
        med = median_from_hist(dp_hist[c])
        ratio = med / baseline if baseline == baseline and baseline else float("nan")
        lo, hi = lobes(ab_hist[c])
        fg = frac_from_lobes(lo, hi, "gain")
        fl = frac_from_lobes(lo, hi, "loss")
        sep = hi - lo

        if AUTOSOME.match(c):
            if ratio > 1.10 and sep > 0.10:
                verdict = "GANANCIA en mosaico (cobertura y BAF coinciden)"
            elif ratio < 0.90 and sep > 0.10:
                verdict = "PERDIDA en mosaico (cobertura y BAF coinciden)"
            elif sep > 0.10:
                verdict = "BAF ancho sin cambio de cobertura: revisar artefacto"
            elif ratio > 1.10 or ratio < 0.90:
                verdict = "cobertura desviada sin BAF partido: probable mapeo"
            else:
                verdict = "disomico"
        else:
            verdict = "sexual, se interpreta aparte"

        print("%-6s %9d %8.1f %7.2f %7.3f %7.3f %8s %8s  %s" %
              (c, n_var[c], med, ratio, lo, hi,
               ("%.2f" % fg) if fg == fg else "-",
               ("%.2f" % fl) if fl == fl else "-", verdict))
        rows[c] = {"snv_pass": n_var[c], "het": n_het[c], "median_dp": med,
                   "dp_ratio": ratio, "lobe_low": lo, "lobe_high": hi,
                   "f_gain": fg, "f_loss": fl, "verdict": verdict,
                   "ab_hist": ab_hist[c]}

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({"baseline_dp": baseline, "pgt_variants": n_pgt,
                               "pid_blocks": len(n_pid_groups), "per_chrom": rows},
                              indent=2), encoding="utf-8")
    print("\nagregados en %s" % OUT.relative_to(ROOT))
    print("Un cromosoma solo cuenta como aneuploide si cobertura y BAF apuntan al mismo lado.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
