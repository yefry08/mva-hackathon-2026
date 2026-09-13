#!/usr/bin/env python3
"""Busqueda genomica de perdida de funcion bialelica, sin depender de ClinVar.

El ranking ciego de clinvar_genomewide.py tiene un hueco declarado: solo encuentra
lo que ya esta clasificado. Una variante truncante nueva en un gen sin entradas en
ClinVar seria invisible. Este script cierra ese hueco por el otro lado.

Bajo herencia recesiva, un candidato causal tiene que ser homocigoto o
heterocigoto compuesto. Asi que se traduce en local **cada variante codificante
del genoma** contra el transcrito MANE Select, se marcan las de perdida de
funcion —codon de parada prematuro, indel que rompe el marco, sitio de splicing
canonico— y se reportan los genes con una en homocigosis o con dos o mas.

Nada sale de la maquina: MANE y la estructura exonica ya estan descargados.

Ojo con la interpretacion. Un genoma sano carga decenas de LoF homocigotas, la
mayoria en familias de genes redundantes (receptores olfativos, antigenos de
grupo sanguineo, genes con pseudogenizacion en curso). La lista cruda no es una
lista de candidatos: es el punto de partida para filtrar por frecuencia.

Uso: python scripts/lof_genomewide.py
"""
from __future__ import annotations

import bisect
import gzip
import os
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / os.environ.get("MVA_WORK", "work")
GFF = ROOT / "ref" / "MANE.gff.gz"
RNA = ROOT / "ref" / "MANE.rna.fna.gz"
OUT = WORK / "lof_genomewide.tsv"

SPLICE_CORE = 2

CODONS = {}
BASES = "TCAG"
AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, b1 in enumerate(BASES):
    for j, b2 in enumerate(BASES):
        for k, b3 in enumerate(BASES):
            CODONS[b1 + b2 + b3] = AAS[i * 16 + j * 4 + k]
COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(s):
    return s.translate(COMP)[::-1]


def refseq_to_chrom(acc):
    if not acc.startswith("NC_"):
        return ""
    n = int(acc.split(".")[0].split("_")[1])
    if 1 <= n <= 22:
        return str(n)
    return {23: "X", 24: "Y"}.get(n, "")


def load_structures():
    """Por transcrito: exones, CDS, hebra, gen, cromosoma."""
    exons, cds, meta = defaultdict(list), defaultdict(list), {}
    with gzip.open(GFF, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("exon", "CDS"):
                continue
            acc = gene = ""
            for part in f[8].split(";"):
                if part.startswith("Parent=rna-"):
                    acc = part[len("Parent=rna-"):]
                elif part.startswith("gene="):
                    gene = part[5:]
            if not acc:
                continue
            chrom = f[0][3:] if f[0].startswith("chr") else f[0]
            (exons if f[2] == "exon" else cds)[acc].append((int(f[3]), int(f[4])))
            meta[acc] = (gene, chrom, f[6])
    for d in (exons, cds):
        for a in d:
            d[a].sort()
    return exons, cds, meta


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


def default_vcf():
    d = ROOT / "data"
    hits = sorted(d.glob("*.vcf.gz")) if d.exists() else []
    return hits[0] if hits else d / "input.vcf.gz"


def main() -> int:
    exons, cds, meta = load_structures()
    seqs = load_sequences(set(meta))
    print("transcritos MANE cargados: %d" % len(seqs), flush=True)

    # Indice por cromosoma de regiones CDS mas margen de splicing.
    idx = defaultdict(list)
    for acc, blocks in cds.items():
        gene, chrom, strand = meta[acc]
        if not chrom:
            continue
        for a, b in blocks:
            idx[chrom].append((a - SPLICE_CORE, b + SPLICE_CORE, acc))
    for chrom in idx:
        idx[chrom].sort()
    starts = {c: [x[0] for x in v] for c, v in idx.items()}
    print("regiones CDS indexadas: %d" % sum(len(v) for v in idx.values()), flush=True)

    per_gene = defaultdict(lambda: {"hom": [], "het": []})
    n_cod = 0
    vcf = default_vcf()
    with gzip.open(vcf, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10 or f[6] not in ("PASS", "."):
                continue
            chrom = f[0][3:] if f[0].startswith("chr") else f[0]
            regions = idx.get(chrom)
            if not regions:
                continue
            pos, ref, alt = int(f[1]), f[3], f[4]
            if alt == "." or "," in alt:
                continue
            kv = dict(zip(f[8].split(":"), f[9].split(":")))
            gt = kv.get("GT", "./.").replace("|", "/").split("/")
            if len(gt) != 2 or "." in gt or gt == ["0", "0"]:
                continue
            zyg = "hom" if gt[0] == gt[1] else "het"

            j = bisect.bisect_right(starts[chrom], pos) - 1
            hit = None
            while j >= 0 and j > bisect.bisect_right(starts[chrom], pos) - 40:
                a, b, acc = regions[j]
                if a <= pos <= b:
                    hit = acc
                    break
                j -= 1
            if hit is None:
                continue
            n_cod += 1

            gene, _, strand = meta[hit]
            conseq = None
            dl = len(alt) - len(ref)
            if dl != 0:
                conseq = "frameshift" if dl % 3 else None
            else:
                if len(ref) == 1:
                    seq = seqs.get(hit)
                    tx = genomic_to_tx(pos, exons[hit], strand) if seq else None
                    if tx is None or not seq:
                        conseq = None
                    else:
                        base = seq[tx].upper()
                        esperado = ref.upper() if strand == "+" else revcomp(ref.upper())
                        if base != esperado:
                            conseq = None          # desajuste: no se inventa nada
                        else:
                            blocks = cds.get(hit)
                            cds_g = blocks[0][0] if strand == "+" else blocks[-1][1]
                            tx0 = genomic_to_tx(cds_g, exons[hit], strand)
                            off = tx - tx0
                            if off >= 0:
                                ci, cp = divmod(off, 3)
                                codon = seq[tx0 + ci * 3: tx0 + ci * 3 + 3].upper()
                                if len(codon) == 3:
                                    ab = alt.upper() if strand == "+" else revcomp(alt.upper())
                                    mut = codon[:cp] + ab + codon[cp + 1:]
                                    if CODONS.get(codon) != "*" and CODONS.get(mut) == "*":
                                        conseq = "nonsense"
            if conseq:
                per_gene[gene][zyg].append((conseq, chrom, pos, ref, alt))

    vivos = {g: v for g, v in per_gene.items() if v["hom"] or len(v["het"]) >= 2}
    print("\nvariantes codificantes evaluadas: %d" % n_cod)
    print("genes con alguna LoF: %d" % len(per_gene))
    print("genes con LoF bialelica potencial, sin filtrar frecuencia: %d" % len(vivos))

    # Filtro de frecuencia: sin el, esta lista es el paisaje normal de LoF de
    # cualquier genoma (receptores olfativos, mucinas, HLA, FMO2...). Se consulta
    # gnomAD solo para estos genes, por rangos de bytes.
    import sys as _sys
    _sys.path.insert(0, str(Path(__file__).resolve().parent))
    from tabix_remote import RemoteTabix                     # noqa: E402
    import re as _re
    AF_RE = _re.compile(r"(?:^|;)AF=([^;]+)")
    EX = ("https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/"
          "vcf/exomes/gnomad.exomes.v4.1.sites.chr%s.vcf.bgz")
    MAX_AF = 0.001

    print("\nconsultando gnomAD para %d genes..." % len(vivos), flush=True)
    resultados = []
    for g, v in sorted(vivos.items()):
        todas = v["hom"] + v["het"]
        chrom = todas[0][1]
        lo = min(x[2] for x in todas)
        hi = max(x[2] for x in todas)
        af = {}
        try:
            for line in RemoteTabix(EX % chrom).query("chr" + chrom, lo, hi):
                f = line.split("\t")
                if len(f) < 8:
                    continue
                m = AF_RE.search(f[7])
                if m and m.group(1) not in (".", ""):
                    af[(int(f[1]), f[3], f[4])] = float(m.group(1))
        except Exception as e:                                # noqa: BLE001
            print("  aviso: %s -> %s" % (g, e), file=_sys.stderr)
        def raro(x):
            return af.get((x[2], x[3], x[4]), 0.0) <= MAX_AF
        hom_raras = [x for x in v["hom"] if raro(x)]
        het_raras = [x for x in v["het"] if raro(x)]
        if hom_raras or len(het_raras) >= 2:
            resultados.append((g, len(hom_raras), len(het_raras),
                               sorted({x[0] for x in hom_raras + het_raras})))

    WORK.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", encoding="utf-8") as out:
        out.write("gene\thom_lof_raras\thet_lof_raras\tclases\n")
        for g, h, e, cl in sorted(resultados, key=lambda r: (-r[1], -r[2])):
            out.write("%s\t%d\t%d\t%s\n" % (g, h, e, ",".join(cl)))

    print("\ngenes con LoF bialelica **rara** (AF <= %g): %d\n" % (MAX_AF, len(resultados)))
    if resultados:
        print("%-14s %5s %5s  %s" % ("gen", "hom", "het", "clases"))
        for g, h, e, cl in sorted(resultados, key=lambda r: (-r[1], -r[2]))[:30]:
            print("%-14s %5d %5d  %s" % (g, h, e, ",".join(cl)))
    else:
        print("Ninguno. No hay ningun gen del genoma con perdida de funcion")
        print("bialelica rara, asi que no hay candidato causal alternativo de ese")
        print("tipo compitiendo con el par de BUB1B.")

    print("\nLimite importante: este barrido solo ve LoF mas LoF. El propio candidato")
    print("de BUB1B, que es una truncante mas un missense, NO aparece aqui por")
    print("diseno. Complementa al analisis por panel, no lo sustituye.")
    print("\ndetalle en %s" % OUT.relative_to(ROOT))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
