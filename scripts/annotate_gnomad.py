#!/usr/bin/env python3
"""Anota las variantes del barrido con frecuencias de gnomAD v4.1. Fase 6 del plan.

Usa el cliente tabix remoto: por cada gen del panel pide a gnomAD el rango de
bytes que cubre el gen, en exomas y en genomas. Lo que sale de esta maquina es la
coordenada de un gen del panel, nunca una posicion del paciente.

Exomas y genomas se consultan por separado a proposito: los exomas tienen mucho
mas poder estadistico en la region codificante (AN del orden de 1.5e5 alelos),
pero no cubren intrones ni UTR. Una variante ausente de exomas puede ser
simplemente intronica, no rara. Los genomas cubren todo con menos poder. La
regla que aplicamos: manda la fuente con mayor AN en ese sitio.

Salida a stdout: solo resumen por gen. El detalle por variante queda en
work/gene_sweep.annot.tsv.

Uso: python scripts/annotate_gnomad.py [--max-af 0.001]
"""
from __future__ import annotations

import argparse
import re
import sys
from collections import defaultdict
import os
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from tabix_remote import RemoteTabix  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / os.environ.get("MVA_WORK", "work")
SWEEP = WORK / "gene_sweep.tsv"
OUT = WORK / "gene_sweep.annot.tsv"

BASE = "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf"
EXOMES = BASE + "/exomes/gnomad.exomes.v4.1.sites.chr%s.vcf.bgz"
GENOMES = BASE + "/genomes/gnomad.genomes.v4.1.sites.chr%s.vcf.bgz"

AF_RE = re.compile(r"(?:^|;)AF=([^;]+)")
AN_RE = re.compile(r"(?:^|;)AN=([^;]+)")
FAF_RE = re.compile(r"(?:^|;)fafmax_faf95_max=([^;]+)")
NHOM_RE = re.compile(r"(?:^|;)nhomalt=([^;]+)")


def parse_info(info: str):
    def grab(rx, cast=float):
        m = rx.search(info)
        if not m or m.group(1) in (".", ""):
            return None
        try:
            return cast(m.group(1))
        except ValueError:
            return None
    return {"af": grab(AF_RE), "an": grab(AN_RE, int),
            "faf": grab(FAF_RE), "nhomalt": grab(NHOM_RE, int)}


def fetch_gene(chrom: str, start: int, end: int):
    """Devuelve {(pos, ref, alt): registro} combinando exomas y genomas."""
    table = {}
    for url, source in ((EXOMES % chrom, "exome"), (GENOMES % chrom, "genome")):
        try:
            for line in RemoteTabix(url).query("chr" + chrom, start, end):
                f = line.split("\t")
                if len(f) < 8:
                    continue
                pos, ref, alt, info = int(f[1]), f[3], f[4], f[7]
                rec = parse_info(info)
                rec["source"] = source
                key = (pos, ref, alt)
                prev = table.get(key)
                # Manda la fuente con mayor numero de alelos observados.
                if prev is None or (rec["an"] or 0) > (prev["an"] or 0):
                    table[key] = rec
        except Exception as e:                      # noqa: BLE001
            print("  aviso: %s chr%s:%d-%d -> %s" % (source, chrom, start, end, e),
                  file=sys.stderr)
    return table


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-af", type=float, default=0.001)
    args = ap.parse_args()

    if not SWEEP.exists():
        print("falta %s: corre antes scripts/gene_sweep.py" % SWEEP)
        return 1

    rows = []
    with SWEEP.open(encoding="utf-8") as fh:
        header = next(fh).rstrip("\n").split("\t")
        for line in fh:
            rows.append(dict(zip(header, line.rstrip("\n").split("\t"))))

    spans = defaultdict(lambda: [None, 10**12, 0])
    for r in rows:
        s = spans[r["gene"]]
        s[0] = r["chrom"]
        s[1] = min(s[1], int(r["pos"]))
        s[2] = max(s[2], int(r["pos"]))

    per_gene = {}
    for gene, (chrom, lo, hi) in sorted(spans.items()):
        print("consultando gnomAD para %s (chr%s:%d-%d)" % (gene, chrom, lo, hi),
              file=sys.stderr)
        per_gene[gene] = fetch_gene(chrom, lo, hi)

    summary = defaultdict(lambda: {"rare_hom": 0, "rare_het": 0, "novel_hom": 0,
                                   "novel_het": 0, "common": 0, "total": 0})
    with OUT.open("w", encoding="utf-8") as out:
        out.write("\t".join(header + ["gnomad_af", "gnomad_an", "gnomad_faf95",
                                      "gnomad_nhomalt", "gnomad_source", "clase"]) + "\n")
        for r in rows:
            key = (int(r["pos"]), r["ref"], r["alt"])
            rec = per_gene.get(r["gene"], {}).get(key)
            gt = r["gt"].replace("|", "/").split("/")
            hom = len(gt) == 2 and gt[0] == gt[1] and gt[0] != "0"
            s = summary[r["gene"]]
            s["total"] += 1
            if rec is None:
                clase = "ausente_en_gnomAD"
                s["novel_hom" if hom else "novel_het"] += 1
            elif rec["af"] is None or rec["af"] <= args.max_af:
                clase = "rara"
                s["rare_hom" if hom else "rare_het"] += 1
            else:
                clase = "comun"
                s["common"] += 1
            out.write("\t".join([r.get(c, "") for c in header] + [
                "" if rec is None or rec["af"] is None else "%.6g" % rec["af"],
                "" if rec is None or rec["an"] is None else str(rec["an"]),
                "" if rec is None or rec["faf"] is None else "%.6g" % rec["faf"],
                "" if rec is None or rec["nhomalt"] is None else str(rec["nhomalt"]),
                "" if rec is None else rec["source"], clase]) + "\n")

    print("\n%-10s %7s %7s %8s %8s %9s %8s  %s" %
          ("gen", "total", "comun", "raraHOM", "raraHET", "ausenteHOM", "ausHET", "lectura recesiva"))
    for gene in sorted(summary):
        s = summary[gene]
        hom = s["rare_hom"] + s["novel_hom"]
        het = s["rare_het"] + s["novel_het"]
        if hom:
            lectura = "HOMOCIGOTO raro: candidato directo"
        elif het >= 2:
            lectura = "2+ raras en het: posible compuesto, hay que fasear"
        elif het == 1:
            lectura = "una sola rara en het: incompleto para recesivo"
        else:
            lectura = "nada raro"
        print("%-10s %7d %7d %8d %8d %9d %8d  %s" %
              (gene, s["total"], s["common"], s["rare_hom"], s["rare_het"],
               s["novel_hom"], s["novel_het"], lectura))

    print("\ndetalle por variante en %s" % OUT.relative_to(ROOT))
    print("Umbral AF <= %g. 'ausente_en_gnomAD' no es sinonimo de raro: en region" % args.max_af)
    print("no codificante puede ser simple falta de cobertura del recurso.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
