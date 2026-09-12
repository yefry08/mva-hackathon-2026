#!/usr/bin/env python3
"""Candidatos recesivos por gen: consecuencia local y fase. Fase 7 del plan.

Toma las variantes raras o ausentes de gnomAD del barrido y responde las dos
preguntas que deciden si un gen sigue vivo bajo modelo recesivo:

1. Donde cae la variante. Se clasifica contra la estructura exonica del
   transcrito MANE Select, en local, con el GFF ya descargado. Nada de mandar
   coordenadas del paciente a un anotador remoto.

2. Si dos variantes del mismo gen estan en trans. GATK deja la fase en PGT/PID:
   PID identifica el bloque y PGT dice el haplotipo. Dos variantes del mismo
   bloque con PGT opuesto (0|1 frente a 1|0) estan en trans, que es lo que exige
   un compuesto heterocigoto. Mismo PGT es cis, y eso mata al candidato: las dos
   copias irian al mismo alelo y el otro quedaria intacto.

A stdout solo va el resumen por gen. El detalle por variante queda en
work/candidates.tsv.

Uso: python scripts/candidates.py
"""
from __future__ import annotations

import gzip
from collections import defaultdict
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / os.environ.get("MVA_WORK", "work")
ANNOT = WORK / "gene_sweep.annot.tsv"
GFF = ROOT / "ref" / "MANE.gff.gz"
GENES = ROOT / "ref" / "sac_genes.txt"
OUT = WORK / "candidates.tsv"

SPLICE_CORE = 2      # +-2 del borde exonico: sitio canonico
SPLICE_REGION = 8    # +-8: region de splicing


def load_genes():
    out = []
    for line in GENES.read_text(encoding="utf-8").splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            out.append(line)
    return set(out)


def load_structure(panel):
    """Por gen: lista de exones y de CDS del transcrito MANE Select."""
    exons, cds = defaultdict(list), defaultdict(list)
    with gzip.open(GFF, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("exon", "CDS"):
                continue
            attrs = f[8]
            gene = ""
            for part in attrs.split(";"):
                if part.startswith("gene="):
                    gene = part[5:]
                    break
            if gene not in panel:
                continue
            rec = (int(f[3]), int(f[4]))
            (exons if f[2] == "exon" else cds)[gene].append(rec)
    for d in (exons, cds):
        for g in d:
            d[g].sort()
    return exons, cds


def classify(pos, exons, cds):
    in_cds = any(a <= pos <= b for a, b in cds)
    in_exon = any(a <= pos <= b for a, b in exons)
    if in_cds:
        return "codificante"
    if in_exon:
        return "exonico UTR"
    best = None
    for a, b in exons:
        for edge in (a, b):
            d = abs(pos - edge)
            best = d if best is None else min(best, d)
    if best is None:
        return "fuera del transcrito"
    if best <= SPLICE_CORE:
        return "sitio de splicing canonico"
    if best <= SPLICE_REGION:
        return "region de splicing"
    return "intronico (%d nt del exon mas cercano)" % best


def phase_verdict(variants):
    """variants: lista de dicts con pgt y pid. Devuelve (veredicto, detalle)."""
    blocks = defaultdict(list)
    sin_fase = 0
    for v in variants:
        if v.get("pid") and v.get("pgt"):
            blocks[v["pid"]].append(v["pgt"])
        else:
            sin_fase += 1
    if not blocks:
        return "sin fase", "ninguna variante rara cae en un bloque PID"
    for pid, pgts in blocks.items():
        if len(pgts) >= 2 and len(set(pgts)) >= 2:
            return "EN TRANS", "dos variantes en el mismo bloque con haplotipo opuesto"
    for pid, pgts in blocks.items():
        if len(pgts) >= 2:
            return "en cis", "dos variantes en el mismo bloque, mismo haplotipo"
    return "bloques distintos", ("%d variantes en %d bloques PID distintos y %d sin "
                                 "fase: la fase no se resuelve por lecturas"
                                 % (len(variants), len(blocks), sin_fase))


def main() -> int:
    if not ANNOT.exists():
        print("falta %s: corre antes scripts/annotate_gnomad.py" % ANNOT)
        return 1
    panel = load_genes()
    exons, cds = load_structure(panel)

    rows = []
    with ANNOT.open(encoding="utf-8") as fh:
        header = next(fh).rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(header, line.rstrip("\n").split("\t")))
            if r.get("clase") in ("rara", "ausente_en_gnomAD"):
                rows.append(r)

    by_gene = defaultdict(list)
    for r in rows:
        gene = r["gene"]
        pos = int(r["pos"])
        r["region"] = classify(pos, exons.get(gene, []), cds.get(gene, []))
        r["tipo"] = "SNV" if len(r["ref"]) == 1 and len(r["alt"]) == 1 else "indel"
        by_gene[gene].append(r)

    with OUT.open("w", encoding="utf-8") as out:
        out.write("\t".join(header + ["region", "tipo"]) + "\n")
        for gene in sorted(by_gene):
            for r in by_gene[gene]:
                out.write("\t".join([r.get(c, "") for c in header] +
                                    [r["region"], r["tipo"]]) + "\n")

    print("%-10s %4s %5s %6s %8s %10s  %s" %
          ("gen", "n", "cod", "splice", "UTR/int", "indel", "fase"))
    for gene in sorted(by_gene, key=lambda g: -len(by_gene[g])):
        vs = by_gene[gene]
        cod = sum(1 for v in vs if v["region"] == "codificante")
        spl = sum(1 for v in vs if "splicing" in v["region"])
        otro = len(vs) - cod - spl
        ind = sum(1 for v in vs if v["tipo"] == "indel")
        verdict, _ = phase_verdict(vs) if len(vs) >= 2 else ("-", "")
        print("%-10s %4d %5d %6d %8d %10d  %s" % (gene, len(vs), cod, spl, otro, ind, verdict))

    print("\n== genes que siguen vivos bajo modelo recesivo ==")
    vivos = 0
    for gene in sorted(by_gene):
        vs = by_gene[gene]
        if len(vs) < 2:
            continue
        verdict, detalle = phase_verdict(vs)
        impacto = [v for v in vs if v["region"] == "codificante" or "splicing" in v["region"]]
        if verdict == "en cis":
            estado = "DESCARTADO: las dos en el mismo haplotipo"
        elif not impacto:
            estado = "debil: ninguna cae en region codificante ni de splicing"
        elif verdict == "EN TRANS":
            estado = "FUERTE: compuesto en trans confirmado por lecturas"
        else:
            estado = "abierto: fase no resuelta, hace falta otra evidencia"
        vivos += 1
        print("  %-10s %d raras, %d con impacto potencial | %s | %s"
              % (gene, len(vs), len(impacto), verdict, estado))
        print("             %s" % detalle)
    if not vivos:
        print("  ninguno")

    print("\ndetalle por variante en %s" % OUT.relative_to(ROOT))
    print("Recordatorio: un bloque PID solo cubre variantes lo bastante cercanas como")
    print("para caer en la misma lectura. Dos variantes separadas por mas que el inserto")
    print("no se pueden fasear asi, y eso no es evidencia de cis ni de trans.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
