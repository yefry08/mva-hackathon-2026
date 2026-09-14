#!/usr/bin/env python3
"""Busqueda genomica ciega de compuestos truncante + missense. Sin panel, sin ClinVar.

Es el hueco que quedaba declarado. El cruce contra ClinVar solo ve lo clasificado;
el barrido de LoF bialelica solo ve truncante con truncante. Ninguno de los dos
cubre un compuesto de alelo nulo mas missense en un gen sin entradas en ClinVar,
que es precisamente la arquitectura del candidato de BUB1B.

Este script no sabe nada de MVA ni del checkpoint mitotico. Recorre el genoma
entero, traduce en local cada variante codificante, y busca genes que carguen al
mismo tiempo una truncante rara y un missense raro, o dos truncantes raras, o una
truncante rara en homocigosis. Si de ahi sale BUB1B, lo hace sin ninguna pista.

Filtro de frecuencia contra gnomAD por rangos de bytes, solo para los genes que ya
tienen al menos una truncante, que es la condicion necesaria de todo el patron.

Uso: python scripts/compound_genomewide.py
"""
from __future__ import annotations

import bisect
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from lof_genomewide import (CODONS, WORK, default_vcf, genomic_to_tx,  # noqa: E402
                            load_sequences, load_structures, revcomp, SPLICE_CORE)
from tabix_remote import query_con_reintentos  # noqa: E402

ROOT = HERE.parent
OUT = WORK / "compound_genomewide.tsv"
MAX_AF = 0.001
AF_RE = re.compile(r"(?:^|;)AF=([^;]+)")
EXOMES = ("https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/"
          "vcf/exomes/gnomad.exomes.v4.1.sites.chr%s.vcf.bgz")
LOF = ("nonsense", "frameshift")


def classify(pos, ref, alt, acc, exons, cds, meta, seqs):
    gene, _, strand = meta[acc]
    dl = len(alt) - len(ref)
    if dl != 0:
        return "frameshift" if dl % 3 else "inframe"
    if len(ref) != 1:
        return None
    seq = seqs.get(acc)
    if not seq:
        return None
    tx = genomic_to_tx(pos, exons[acc], strand)
    if tx is None:
        return None
    esperado = ref.upper() if strand == "+" else revcomp(ref.upper())
    if seq[tx].upper() != esperado:
        return None                     # desajuste de referencia: no se inventa nada
    blocks = cds.get(acc)
    cds_g = blocks[0][0] if strand == "+" else blocks[-1][1]
    tx0 = genomic_to_tx(cds_g, exons[acc], strand)
    off = tx - tx0
    if off < 0:
        return None
    ci, cp = divmod(off, 3)
    codon = seq[tx0 + ci * 3: tx0 + ci * 3 + 3].upper()
    if len(codon) != 3:
        return None
    ab = alt.upper() if strand == "+" else revcomp(alt.upper())
    mut = codon[:cp] + ab + codon[cp + 1:]
    a, b = CODONS.get(codon), CODONS.get(mut)
    if a is None or b is None or a == "*":
        return None
    if b == "*":
        return "nonsense"
    return "missense" if a != b else "sinonima"


def main() -> int:
    exons, cds, meta = load_structures()
    seqs = load_sequences(set(meta))
    idx = defaultdict(list)
    for acc, blocks in cds.items():
        _, chrom, _ = meta[acc]
        if chrom:
            for a, b in blocks:
                idx[chrom].append((a - SPLICE_CORE, b + SPLICE_CORE, acc))
    for c in idx:
        idx[c].sort()
    starts = {c: [x[0] for x in v] for c, v in idx.items()}
    print("transcritos: %d" % len(seqs), flush=True)

    per_gene = defaultdict(list)
    with gzip.open(default_vcf(), "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10 or f[6] not in ("PASS", "."):
                continue
            chrom = f[0][3:] if f[0].startswith("chr") else f[0]
            if chrom not in idx:
                continue
            pos, ref, alt = int(f[1]), f[3], f[4]
            if alt == "." or "," in alt:
                continue
            kv = dict(zip(f[8].split(":"), f[9].split(":")))
            gt = kv.get("GT", "./.").replace("|", "/").split("/")
            if len(gt) != 2 or "." in gt or gt == ["0", "0"]:
                continue
            j = bisect.bisect_right(starts[chrom], pos) - 1
            hit = None
            lim = j - 40
            while j >= 0 and j > lim:
                a, b, acc = idx[chrom][j]
                if a <= pos <= b:
                    hit = acc
                    break
                j -= 1
            if hit is None:
                continue
            conseq = classify(pos, ref, alt, hit, exons, cds, meta, seqs)
            if conseq in ("nonsense", "frameshift", "missense"):
                zyg = "hom" if gt[0] == gt[1] else "het"
                per_gene[meta[hit][0]].append((conseq, chrom, pos, ref, alt, zyg))

    con_lof = {g: v for g, v in per_gene.items() if any(x[0] in LOF for x in v)}
    print("genes con variantes codificantes no sinonimas: %d" % len(per_gene))
    print("de ellos con al menos una truncante: %d" % len(con_lof), flush=True)
    print("consultando gnomAD para esos genes...", flush=True)

    ranking, sin_datos = [], []
    for g, v in sorted(con_lof.items()):
        chrom = v[0][1]
        lo, hi = min(x[2] for x in v), max(x[2] for x in v)
        af = {}
        try:
            lineas = query_con_reintentos(EXOMES % chrom, "chr" + chrom, lo, hi)
        except RuntimeError as e:
            # Antes el gen se saltaba en silencio: un fallo de red podia sacar del
            # ranking al candidato verdadero sin que nadie se enterara.
            print("  SIN DATOS %s: %s" % (g, e), file=sys.stderr)
            sin_datos.append(g)
            continue
        for line in lineas:
            f = line.split("\t")
            if len(f) >= 8:
                m = AF_RE.search(f[7])
                if m and m.group(1) not in (".", ""):
                    af[(int(f[1]), f[3], f[4])] = float(m.group(1))
        raras = [x for x in v if af.get((x[2], x[3], x[4]), 0.0) <= MAX_AF]
        lof_hom = [x for x in raras if x[0] in LOF and x[5] == "hom"]
        lof_het = [x for x in raras if x[0] in LOF and x[5] == "het"]
        mis_het = [x for x in raras if x[0] == "missense" and x[5] == "het"]

        if lof_hom:
            patron = "truncante homocigota"
        elif len(lof_het) >= 2:
            patron = "dos truncantes"
        elif lof_het and mis_het:
            patron = "truncante + missense"
        else:
            continue
        ranking.append((g, patron, len(lof_hom), len(lof_het), len(mis_het)))

    orden = {"truncante homocigota": 0, "dos truncantes": 1, "truncante + missense": 2}
    ranking.sort(key=lambda r: (orden[r[1]], -(r[2] + r[3] + r[4]), r[0]))

    WORK.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", encoding="utf-8") as out:
        out.write("gene\tpatron\tlof_hom\tlof_het\tmissense_het\n")
        for r in ranking:
            out.write("%s\t%s\t%d\t%d\t%d\n" % r)

    print("\ngenes con truncante que gnomAD no pudo anotar: %d%s"
          % (len(sin_datos), (" -> " + ", ".join(sin_datos)) if sin_datos else ""))
    if sin_datos:
        print("RESULTADO INCOMPLETO: esos genes no se evaluaron. Repetir antes de interpretar.")
    print("genes con patron recesivo compatible, todo raro (AF <= %g): %d\n" % (MAX_AF, len(ranking)))
    print("%-14s %-22s %7s %7s %9s" % ("gen", "patron", "LoF hom", "LoF het", "miss het"))
    for r in ranking:
        marca = "   <- candidato del panel" if r[0] == "BUB1B" else ""
        print("%-14s %-22s %7d %7d %9d%s" % (r + (marca,)))

    tm = [r for r in ranking if r[1] == "truncante + missense"]
    print("\ncon el patron exacto truncante + missense: %d genes" % len(tm))
    if any(r[0] == "BUB1B" for r in ranking):
        print("BUB1B aparece sin panel, sin ClinVar y sin ninguna pista sobre MVA.")
    else:
        print("BUB1B NO aparece. Revisar: o el filtro lo pierde, o el hallazgo del panel")
        print("dependia de algo que este barrido no reproduce.")
    print("\ndetalle en %s" % OUT.relative_to(ROOT))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
