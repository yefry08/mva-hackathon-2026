#!/usr/bin/env python3
"""Barrido dirigido de los genes del checkpoint mitotico. Fase 7b del plan.

Recorre el VCF una vez y reparte las variantes entre los genes de
ref/sac_genes.txt, usando las coordenadas GRCh38 del transcrito MANE Select.

Lo variante por variante se queda en disco (work/gene_sweep.tsv). A la salida
estandar solo va el resumen por gen: simbolo, conteos, profundidad. Los nombres
de gen y los rankings de genes son justamente lo que las reglas del hackathon
dejan conservar y publicar; las coordenadas del nino, no.

MVA es autosomica recesiva, asi que lo que interesa por gen es: un homocigoto
alterno, o dos o mas heterocigotos que puedan estar en trans. La frecuencia
poblacional todavia no entra: sin gnomAD local, un conteo alto de heterocigotos
son polimorfismos comunes, no candidatos.

Uso: python scripts/gene_sweep.py [ruta.vcf.gz]
"""
from __future__ import annotations

import gzip
import sys
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / os.environ.get("MVA_WORK", "work")
def _default_vcf():
    """Descubre el VCF en data/ sin escribir el nombre de la muestra en el codigo."""
    d = ROOT / "data"
    hits = sorted(p for p in d.glob("*.vcf.gz")) if d.exists() else []
    return hits[0] if hits else d / "input.vcf.gz"


DEFAULT_VCF = _default_vcf()
MANE = ROOT / "ref" / "MANE.summary.txt.gz"
GENES = ROOT / "ref" / "sac_genes.txt"
OUT_TSV = WORK / "gene_sweep.tsv"

PAD = 2000   # margen para captar sitios de splicing y UTR proximos


def refseq_to_chrom(acc: str) -> str:
    """NC_000019.10 -> 19. El VCF usa contigs sin prefijo."""
    if not acc.startswith("NC_"):
        return ""
    n = int(acc.split(".")[0].split("_")[1])
    if 1 <= n <= 22:
        return str(n)
    return {23: "X", 24: "Y", 12920: "MT"}.get(n, "")


def load_genes():
    wanted = []
    for line in GENES.read_text(encoding="utf-8").splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            wanted.append(line)
    spans = {}
    with gzip.open(MANE, "rt", encoding="utf-8", errors="replace") as fh:
        cols = next(fh).rstrip("\n").split("\t")
        i_sym, i_chr = cols.index("symbol"), cols.index("GRCh38_chr")
        i_start, i_end = cols.index("chr_start"), cols.index("chr_end")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            sym = f[i_sym]
            if sym not in wanted:
                continue
            chrom = refseq_to_chrom(f[i_chr])
            if not chrom:
                continue
            spans[sym] = (chrom, int(f[i_start]) - PAD, int(f[i_end]) + PAD)
    faltan = [g for g in wanted if g not in spans]
    if faltan:
        print("sin coordenadas MANE: %s" % ", ".join(faltan), file=sys.stderr)
    return wanted, spans


def main() -> int:
    vcf = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_VCF
    wanted, spans = load_genes()

    by_chrom = {}
    for sym, (chrom, start, end) in spans.items():
        by_chrom.setdefault(chrom, []).append((start, end, sym))
    for chrom in by_chrom:
        by_chrom[chrom].sort()

    stats = {sym: {"pass": 0, "filtered": 0, "hom_alt": 0, "het": 0, "hom_ref": 0,
                   "indel": 0, "phased": 0, "dp_sum": 0, "dp_n": 0, "pid": set()}
             for sym in spans}

    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    with gzip.open(vcf, "rt", encoding="utf-8", errors="replace") as fh, \
            OUT_TSV.open("w", encoding="utf-8") as out:
        out.write("gene\tchrom\tpos\tref\talt\tfilter\tgt\tad\tdp\tgq\tpgt\tpid\n")
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            chrom = f[0][3:] if f[0].startswith("chr") else f[0]
            spans_here = by_chrom.get(chrom)
            if not spans_here:
                continue
            pos = int(f[1])
            hit = None
            for start, end, sym in spans_here:
                if start <= pos <= end:
                    hit = sym
                    break
                if start > pos:
                    break
            if hit is None:
                continue

            ref, alt, filt, fmt, sample = f[3], f[4], f[6], f[8], f[9]
            kv = dict(zip(fmt.split(":"), sample.split(":")))
            gt = kv.get("GT", "./.")
            a = gt.replace("|", "/").split("/")
            s = stats[hit]
            if filt in ("PASS", "."):
                s["pass"] += 1
            else:
                s["filtered"] += 1
            if len(ref) != 1 or len(alt) != 1:
                s["indel"] += 1
            if len(a) == 2 and "." not in a:
                if a[0] != a[1]:
                    s["het"] += 1
                elif a[0] != "0":
                    s["hom_alt"] += 1
                else:
                    s["hom_ref"] += 1
            if "|" in gt or kv.get("PGT", ".") not in (".", ""):
                s["phased"] += 1
            if kv.get("PID", ".") not in (".", ""):
                s["pid"].add(kv["PID"])
            dp = kv.get("DP", "")
            if dp.isdigit():
                s["dp_sum"] += int(dp)
                s["dp_n"] += 1

            if filt in ("PASS", ".") and len(a) == 2 and "." not in a and a != ["0", "0"]:
                out.write("\t".join([hit, chrom, str(pos), ref, alt, filt, gt,
                                     kv.get("AD", ""), kv.get("DP", ""),
                                     kv.get("GQ", ""), kv.get("PGT", ""),
                                     kv.get("PID", "")]) + "\n")
                written += 1

    print("%-10s %5s %8s %7s %7s %7s %7s %8s %7s" %
          ("gen", "chr", "span kb", "PASS", "homALT", "het", "indel", "bloqPID", "DPmed"))
    for sym in wanted:
        if sym not in stats:
            print("%-10s  sin coordenadas MANE" % sym)
            continue
        chrom, start, end = spans[sym]
        s = stats[sym]
        dp = (s["dp_sum"] / s["dp_n"]) if s["dp_n"] else float("nan")
        print("%-10s %5s %8.1f %7d %7d %7d %7d %8d %7.1f" %
              (sym, chrom, (end - start) / 1000.0, s["pass"], s["hom_alt"],
               s["het"], s["indel"], len(s["pid"]), dp))

    tot_hom = sum(s["hom_alt"] for s in stats.values())
    print("\n%d variantes no de referencia escritas en %s" % (written, OUT_TSV.relative_to(ROOT)))
    print("%d genotipos homocigotos alternos en el panel." % tot_hom)
    print("\nOjo con la lectura: sin frecuencias poblacionales, la mayoria de esos")
    print("homocigotos y heterocigotos son polimorfismos comunes. El paso que convierte")
    print("esta tabla en candidatos es anotar gnomAD y quedarse con AF < 0.001.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
