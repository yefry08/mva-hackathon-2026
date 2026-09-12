#!/usr/bin/env python3
"""Cruza los candidatos contra ClinVar. Fase 6 del plan.

ClinVar se descarga entero (190 MB) y se consulta en local: asi ninguna posicion
del paciente sale de la maquina, ni siquiera como rango de bytes.

Dos cruces distintos, que responden preguntas distintas:

1. Coincidencia exacta de variante. Si una de nuestras candidatas ya esta
   clasificada, eso es evidencia directa y ademas es lo que sostiene el criterio
   PS1 o PM5 de ACMG.
2. Contexto del gen. Cuantas variantes patogenicas tiene ese gen en ClinVar y con
   que enfermedad se asocian. Un gen con patogenicas conocidas para el fenotipo
   correcto pesa distinto de uno sin ninguna.

A stdout solo va nivel gen y clasificacion. El detalle queda en
work/candidates.clinvar.tsv.

Uso: python scripts/clinvar_annot.py
"""
from __future__ import annotations

import gzip
import re
from collections import Counter, defaultdict
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / os.environ.get("MVA_WORK", "work")
CAND = WORK / "candidates.tsv"
CLINVAR = ROOT / "ref" / "clinvar.vcf.gz"
GENES = ROOT / "ref" / "sac_genes.txt"
OUT = WORK / "candidates.clinvar.tsv"

GENEINFO = re.compile(r"(?:^|;)GENEINFO=([^;:]+)")
CLNSIG = re.compile(r"(?:^|;)CLNSIG=([^;]+)")
CLNDN = re.compile(r"(?:^|;)CLNDN=([^;]+)")
CLNREVSTAT = re.compile(r"(?:^|;)CLNREVSTAT=([^;]+)")

PATOGENICO = ("Pathogenic", "Likely_pathogenic")


def load_panel():
    out = []
    for line in GENES.read_text(encoding="utf-8").splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            out.append(line)
    return set(out)


def main() -> int:
    if not CAND.exists():
        print("falta %s: corre antes scripts/candidates.py" % CAND)
        return 1
    panel = load_panel()

    rows, header = [], None
    with CAND.open(encoding="utf-8") as fh:
        header = next(fh).rstrip("\n").split("\t")
        for line in fh:
            rows.append(dict(zip(header, line.rstrip("\n").split("\t"))))
    wanted = {(r["chrom"], int(r["pos"])): r for r in rows}

    exact = {}                       # (chrom,pos,ref,alt) -> registro ClinVar
    same_pos = defaultdict(list)     # otra variante en la misma posicion (PM5)
    gene_stats = defaultdict(Counter)
    gene_dx = defaultdict(Counter)

    with gzip.open(CLINVAR, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            info = f[7]
            m = GENEINFO.search(info)
            gene = m.group(1) if m else ""
            sig = (CLNSIG.search(info).group(1) if CLNSIG.search(info) else "")

            if gene in panel:
                gene_stats[gene][sig.split("|")[0]] += 1
                if sig.startswith(PATOGENICO):
                    dn = CLNDN.search(info)
                    if dn:
                        for d in dn.group(1).split("|"):
                            gene_dx[gene][d] += 1

            key = (f[0], int(f[1]))
            if key in wanted:
                rec = {"sig": sig, "rev": (CLNREVSTAT.search(info).group(1)
                                           if CLNREVSTAT.search(info) else ""),
                       "dn": (CLNDN.search(info).group(1) if CLNDN.search(info) else ""),
                       "ref": f[3], "alt": f[4]}
                r = wanted[key]
                if f[3] == r["ref"] and f[4] == r["alt"]:
                    exact[(r["chrom"], int(r["pos"]), r["ref"], r["alt"])] = rec
                else:
                    same_pos[(r["chrom"], int(r["pos"]))].append(rec)

    with OUT.open("w", encoding="utf-8") as out:
        out.write("\t".join(header + ["clinvar_sig", "clinvar_rev", "clinvar_dn",
                                      "clinvar_misma_pos"]) + "\n")
        for r in rows:
            k = (r["chrom"], int(r["pos"]), r["ref"], r["alt"])
            rec = exact.get(k)
            otras = same_pos.get((r["chrom"], int(r["pos"])), [])
            out.write("\t".join([r.get(c, "") for c in header] + [
                rec["sig"] if rec else "", rec["rev"] if rec else "",
                rec["dn"] if rec else "",
                ";".join("%s>%s:%s" % (o["ref"], o["alt"], o["sig"]) for o in otras)]) + "\n")

    print("== coincidencias exactas de nuestras candidatas en ClinVar ==")
    if exact:
        for (chrom, pos, ref, alt), rec in exact.items():
            gene = next(r["gene"] for r in rows
                        if r["chrom"] == chrom and int(r["pos"]) == pos)
            print("  %-10s %-28s revision: %s" % (gene, rec["sig"], rec["rev"]))
            print("             enfermedad: %s" % rec["dn"][:110])
    else:
        print("  ninguna. Las candidatas no estan clasificadas en ClinVar.")

    print("\n== otra variante clasificada en la misma posicion (soporta PM5) ==")
    if same_pos:
        for (chrom, pos), otras in same_pos.items():
            gene = next(r["gene"] for r in rows
                        if r["chrom"] == chrom and int(r["pos"]) == pos)
            sigs = ", ".join(sorted({o["sig"].split("|")[0] for o in otras}))
            print("  %-10s %d variantes en el mismo sitio: %s" % (gene, len(otras), sigs))
    else:
        print("  ninguna")

    print("\n== contexto de gen: patogenicas conocidas en ClinVar ==")
    print("%-10s %6s %8s  %s" % ("gen", "total", "patog", "enfermedad dominante"))
    genes_con_cand = sorted({r["gene"] for r in rows})
    for gene in genes_con_cand:
        c = gene_stats.get(gene, Counter())
        total = sum(c.values())
        patog = sum(n for s, n in c.items() if s.startswith(PATOGENICO))
        dx = gene_dx.get(gene, Counter()).most_common(1)
        dx = dx[0][0][:58] if dx else "-"
        print("%-10s %6d %8d  %s" % (gene, total, patog, dx))

    print("\ndetalle en %s" % OUT.relative_to(ROOT))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
