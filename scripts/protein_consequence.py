#!/usr/bin/env python3
"""Consecuencia proteica de las variantes candidatas, calculada en local.

Todo se resuelve con recursos publicos ya descargados: la estructura exonica del
transcrito MANE Select (GFF) y la secuencia del transcrito (FASTA). Ninguna
posicion del paciente se le manda a un anotador remoto; la unica descarga fue de
archivos de referencia completos, iguales para todo el mundo.

El mapeo es el de siempre: posicion genomica -> coordenada en el transcrito
sumando exones en el sentido de la hebra -> desplazamiento dentro del CDS ->
codon. Antes de traducir nada se comprueba que la base de referencia del VCF
coincida con la del transcrito en esa posicion. Si no coincide, el mapeo esta mal
y el resultado no vale: se reporta el desajuste en vez de inventar una proteina.

Uso: python scripts/protein_consequence.py
"""
from __future__ import annotations

import gzip
from collections import defaultdict
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / os.environ.get("MVA_WORK", "work")
CAND = WORK / "candidates.clinvar.tsv"
GFF = ROOT / "ref" / "MANE.gff.gz"
RNA = ROOT / "ref" / "MANE.rna.fna.gz"
SUMMARY = ROOT / "ref" / "MANE.summary.txt.gz"
OUT = WORK / "candidates.protein.tsv"

CODONS = {}
BASES = "TCAG"
AAS = ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")
for i, b1 in enumerate(BASES):
    for j, b2 in enumerate(BASES):
        for k, b3 in enumerate(BASES):
            CODONS[b1 + b2 + b3] = AAS[i * 16 + j * 4 + k]

THREE = {"A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln",
         "E": "Glu", "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys",
         "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp",
         "Y": "Tyr", "V": "Val", "*": "Ter"}
COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(s):
    return s.translate(COMP)[::-1]


def load_tx_for_genes(genes):
    tx = {}
    with gzip.open(SUMMARY, "rt", encoding="utf-8", errors="replace") as fh:
        cols = next(fh).rstrip("\n").split("\t")
        i_sym, i_nuc = cols.index("symbol"), cols.index("RefSeq_nuc")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[i_sym] in genes:
                tx[f[i_sym]] = f[i_nuc]
    return tx


def load_structure(accs):
    exons, cds, strand = defaultdict(list), defaultdict(list), {}
    with gzip.open(GFF, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("exon", "CDS"):
                continue
            parent = ""
            for part in f[8].split(";"):
                if part.startswith("Parent=rna-"):
                    parent = part[len("Parent=rna-"):]
                    break
            if parent not in accs:
                continue
            (exons if f[2] == "exon" else cds)[parent].append((int(f[3]), int(f[4])))
            strand[parent] = f[6]
    for d in (exons, cds):
        for a in d:
            d[a].sort()
    return exons, cds, strand


def load_sequences(accs):
    seqs, cur, buf = {}, None, []
    with gzip.open(RNA, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur in accs:
                    seqs[cur] = "".join(buf)
                cur, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
        if cur in accs:
            seqs[cur] = "".join(buf)
    return seqs


def genomic_to_tx(pos, exons, strand):
    """Desplazamiento 0-based dentro del transcrito, o None si cae en intron."""
    if strand == "+":
        off = 0
        for a, b in exons:
            if a <= pos <= b:
                return off + (pos - a)
            off += b - a + 1
    else:
        off = 0
        for a, b in reversed(exons):
            if a <= pos <= b:
                return off + (b - pos)
            off += b - a + 1
    return None


def main() -> int:
    rows, header = [], None
    with CAND.open(encoding="utf-8") as fh:
        header = next(fh).rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(header, line.rstrip("\n").split("\t")))
            if r["region"] in ("codificante", "exonico UTR") or "splicing" in r["region"]:
                rows.append(r)

    genes = {r["gene"] for r in rows}
    tx = load_tx_for_genes(genes)
    accs = set(tx.values())
    exons, cds, strand = load_structure(accs)
    seqs = load_sequences(accs)

    print("%-9s %-6s %-8s %-22s %s" % ("gen", "tipo", "clase", "consecuencia", "nota"))
    out_rows = []
    for r in sorted(rows, key=lambda r: (r["gene"], int(r["pos"]))):
        gene, pos, ref, alt = r["gene"], int(r["pos"]), r["ref"], r["alt"]
        acc = tx.get(gene)
        nota, conseq = "", ""
        if not acc or acc not in seqs:
            print("%-9s %-6s %-8s %-22s %s" % (gene, r["tipo"], r["clase"][:8],
                                               "sin transcrito", acc or "-"))
            continue
        st = strand[acc]
        seq = seqs[acc]
        tx_pos = genomic_to_tx(pos, exons[acc], st)
        if tx_pos is None:
            conseq, nota = "intronica", "no cae en exon del MANE Select"
        else:
            base = seq[tx_pos].upper()
            esperado = ref.upper() if st == "+" else revcomp(ref.upper())
            if base != esperado:
                conseq = "MAPEO INCONSISTENTE"
                nota = "el transcrito trae %s donde el VCF dice %s" % (base, esperado)
            else:
                cds_list = cds.get(acc, [])
                if not cds_list:
                    conseq, nota = "no codificante", "el transcrito no tiene CDS"
                else:
                    cds_start_g = cds_list[0][0] if st == "+" else cds_list[-1][1]
                    tx_cds0 = genomic_to_tx(cds_start_g, exons[acc], st)
                    cds_off = tx_pos - tx_cds0
                    if cds_off < 0:
                        conseq, nota = "UTR 5 prima", ""
                    else:
                        codon_i, codon_p = divmod(cds_off, 3)
                        s = tx_cds0 + codon_i * 3
                        codon = seq[s:s + 3].upper()
                        if len(codon) < 3:
                            conseq, nota = "fuera del marco", "codon incompleto"
                        else:
                            alt_b = alt.upper() if st == "+" else revcomp(alt.upper())
                            mut = list(codon)
                            mut[codon_p] = alt_b
                            mut = "".join(mut)
                            aa_ref, aa_alt = CODONS.get(codon, "?"), CODONS.get(mut, "?")
                            aa_pos = codon_i + 1
                            if aa_ref == aa_alt:
                                conseq = "sinonima"
                                nota = "p.%s%d=" % (THREE.get(aa_ref, aa_ref), aa_pos)
                            elif aa_alt == "*":
                                conseq = "NONSENSE (stop ganado)"
                                nota = "p.%s%dTer" % (THREE.get(aa_ref, aa_ref), aa_pos)
                            elif aa_ref == "*":
                                conseq = "perdida de stop"
                                nota = "p.Ter%d%s" % (aa_pos, THREE.get(aa_alt, aa_alt))
                            else:
                                conseq = "missense"
                                nota = "p.%s%d%s" % (THREE.get(aa_ref, aa_ref), aa_pos,
                                                     THREE.get(aa_alt, aa_alt))
                            if aa_pos == 1:
                                nota += " (codon de inicio)"
        print("%-9s %-6s %-8s %-22s %s" % (gene, r["tipo"], r["clase"][:8], conseq,
                                           (nota + "  " + (r.get("clinvar_sig") or ""))[:60]))
        r["consecuencia"], r["proteina"] = conseq, nota
        out_rows.append(r)

    with OUT.open("w", encoding="utf-8") as out:
        out.write("\t".join(header + ["consecuencia", "proteina"]) + "\n")
        for r in out_rows:
            out.write("\t".join([r.get(c, "") for c in header + ["consecuencia", "proteina"]]) + "\n")
    print("\ndetalle en %s" % OUT.relative_to(ROOT))
    print("Este calculo cubre SNV en el MANE Select. Indels, isoformas alternativas y")
    print("efectos de splicing profundos necesitan otra herramienta y quedan pendientes.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
