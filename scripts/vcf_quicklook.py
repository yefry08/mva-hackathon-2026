#!/usr/bin/env python3
"""Vistazo agregado al VCF entregado. Tarea G09 + G22 del plan.

Produce SOLO agregados: metadatos del header, conteos por cromosoma y la
distribucion del balance alelico de los heterocigotos. Ningun genotipo individual
sale de aqui, asi que la salida es segura para mirar y para el writeup.

El balance alelico por cromosoma es el test rapido de aneuploidia en mosaico: en
un cromosoma disomico los hets se agrupan en 0.5; una ganancia o perdida en una
fraccion de las celulas corre la mediana y parte la distribucion en dos lobulos.

Uso: python scripts/vcf_quicklook.py [ruta.vcf.gz]
"""
from __future__ import annotations

import gzip
import json
import re
import statistics
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
OUT = ROOT / "work" / "quicklook.json"

MIN_DP = 20          # profundidad minima para que el balance alelico signifique algo
NBINS = 20
def _mask_pattern():
    """Los identificadores a enmascarar se leen de .claude/sample_ids.txt, que
    esta en .gitignore. Escribirlos aqui seria meter en el repositorio publico
    justo lo que este regex existe para ocultar."""
    ids = ROOT / ".claude" / "sample_ids.txt"
    pats = []
    if ids.exists():
        for line in ids.read_text(encoding="utf-8").splitlines():
            line = line.split("#", 1)[0].strip()
            if line:
                pats.append(re.escape(line))
    pats.append(r"EX\d{6,}")
    return re.compile("|".join(pats))


MASK = _mask_pattern()


def mask(s: str) -> str:
    return MASK.sub("<MUESTRA>", s)


def main() -> int:
    vcf = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_VCF
    if not vcf.exists():
        print("no existe %s" % vcf)
        return 1

    meta = {"reference": None, "contigs": 0, "contig_style": None, "source": [],
            "samples": 0, "format_keys": set(), "info_keys": set()}
    per_chrom = defaultdict(lambda: {
        "variants": 0, "pass": 0, "het": 0, "hom_alt": 0, "snv": 0, "indel": 0,
        "ab": [], "dp_sum": 0, "dp_n": 0, "phased": 0, "pgt": 0,
    })

    with gzip.open(vcf, "rt", encoding="utf-8", errors="replace") as fh:
        gt_idx = ad_idx = dp_idx = None
        for line in fh:
            if line.startswith("##"):
                if line.startswith("##reference"):
                    meta["reference"] = line.strip().split("=", 1)[1]
                elif line.startswith("##contig"):
                    meta["contigs"] += 1
                    if meta["contig_style"] is None:
                        m = re.search(r"ID=([^,>]+)", line)
                        if m:
                            meta["contig_style"] = "chr" if m.group(1).startswith("chr") else "sin chr"
                elif line.startswith(("##source", "##GATKCommandLine", "##DRAGEN")):
                    meta["source"].append(line.strip()[:160])
                elif line.startswith("##FORMAT"):
                    m = re.search(r"ID=([^,]+)", line)
                    if m:
                        meta["format_keys"].add(m.group(1))
                elif line.startswith("##INFO"):
                    m = re.search(r"ID=([^,]+)", line)
                    if m:
                        meta["info_keys"].add(m.group(1))
                continue
            if line.startswith("#CHROM"):
                meta["samples"] = max(0, len(line.rstrip("\n").split("\t")) - 9)
                continue

            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            chrom, _, _, ref, alt, _, filt, _, fmt, sample = f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7], f[8], f[9]
            c = per_chrom[chrom]
            c["variants"] += 1
            if filt in ("PASS", "."):
                c["pass"] += 1
            if len(ref) == 1 and len(alt) == 1:
                c["snv"] += 1
            else:
                c["indel"] += 1

            keys = fmt.split(":")
            try:
                gt_idx = keys.index("GT")
            except ValueError:
                continue
            vals = sample.split(":")
            if gt_idx >= len(vals):
                continue
            gt = vals[gt_idx]
            if "|" in gt:
                c["phased"] += 1
            if "PGT" in keys:
                c["pgt"] += 1

            a = gt.replace("|", "/").split("/")
            if len(a) != 2 or "." in a:
                continue
            if a[0] != a[1]:
                c["het"] += 1
            elif a[0] != "0":
                c["hom_alt"] += 1
                continue
            else:
                continue

            if "AD" in keys:
                ad_idx = keys.index("AD")
                if ad_idx < len(vals):
                    parts = vals[ad_idx].split(",")
                    if len(parts) >= 2 and parts[0].isdigit() and parts[1].isdigit():
                        r, v = int(parts[0]), int(parts[1])
                        dp = r + v
                        if dp >= MIN_DP:
                            c["ab"].append(v / dp)
                            c["dp_sum"] += dp
                            c["dp_n"] += 1

    def sortkey(ch):
        base = ch.replace("chr", "")
        return (0, int(base)) if base.isdigit() else (1, base)

    main_chroms = [ch for ch in per_chrom
                   if re.fullmatch(r"(chr)?([1-9]|1\d|2[0-2]|X|Y|MT?)", ch)]
    main_chroms.sort(key=sortkey)

    print("== header ==")
    print("referencia      : %s" % mask(str(meta["reference"])))
    print("contigs         : %d (%s)" % (meta["contigs"], meta["contig_style"]))
    print("muestras        : %d" % meta["samples"])
    print("FORMAT          : %s" % ",".join(sorted(meta["format_keys"])))
    print("fase PGT/PID    : %s" % ("si" if {"PGT", "PID"} <= meta["format_keys"] else "no"))
    for s in meta["source"][:3]:
        print("origen          : %s" % mask(s))

    print("\n== por cromosoma (hets con DP>=%d) ==" % MIN_DP)
    print("%-6s %10s %9s %8s %8s %8s %7s" %
          ("chrom", "variantes", "hets", "medAB", "fuera", "DPmed", "fase"))
    summary = {}
    for ch in main_chroms:
        c = per_chrom[ch]
        ab = c["ab"]
        med = statistics.median(ab) if ab else float("nan")
        outside = (sum(1 for x in ab if x < 0.40 or x > 0.60) / len(ab)) if ab else float("nan")
        dpm = (c["dp_sum"] / c["dp_n"]) if c["dp_n"] else float("nan")
        print("%-6s %10d %9d %8.3f %7.1f%% %8.1f %6.1f%%" %
              (ch, c["variants"], c["het"], med, outside * 100, dpm,
               100.0 * c["phased"] / max(1, c["variants"])))
        summary[ch] = {"variants": c["variants"], "het": c["het"], "hom_alt": c["hom_alt"],
                       "snv": c["snv"], "indel": c["indel"], "median_ab": med,
                       "frac_ab_outside": outside, "mean_dp": dpm,
                       "phased_frac": c["phased"] / max(1, c["variants"]),
                       "ab_hist": _hist(ab)}

    total = sum(c["variants"] for c in per_chrom.values())
    print("\ntotal de variantes: %d en %d contigs" % (total, len(per_chrom)))
    print("\nLectura: en un cromosoma disomico la mediana de AB queda cerca de 0.500 y la")
    print("fraccion fuera de [0.40,0.60] es baja. Una desviacion sostenida en un cromosoma")
    print("entero es la firma de aneuploidia en mosaico, y es evidencia ortogonal fuerte.")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    meta["format_keys"] = sorted(meta["format_keys"])
    meta["info_keys"] = sorted(meta["info_keys"])
    OUT.write_text(json.dumps({"meta": meta, "per_chrom": summary}, indent=2), encoding="utf-8")
    print("\nagregados en %s" % OUT.relative_to(ROOT))
    return 0


def _hist(ab):
    h = [0] * NBINS
    for x in ab:
        i = min(NBINS - 1, int(x * NBINS))
        h[i] += 1
    return h


if __name__ == "__main__":
    sys.exit(main())
