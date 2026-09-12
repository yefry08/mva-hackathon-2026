#!/usr/bin/env python3
"""Cruce genomico completo contra ClinVar, en local. Ranking ciego y hallazgos secundarios.

El barrido por panel mira 32 genes, asi que por construccion no puede encontrar
nada fuera de ellos. Esto recorre el VCF entero y lo cruza contra todas las
variantes patogenicas o probablemente patogenicas de ClinVar, sin prior de genes.
Sirve para dos cosas:

1. Ranking ciego. Si la mejor evidencia del genoma completo sigue siendo la misma
   que encontro el panel, el resultado no depende del prior, y eso es justo lo que
   hay que demostrar en el writeup.
2. Hallazgos secundarios. El hackathon los admite y no penalizan el puntaje
   automatico, pero solo cuentan los que esten bien justificados.

Todo local: ClinVar ya esta descargado y el VCF nunca sale de la maquina.

Uso: python scripts/clinvar_genomewide.py
"""
from __future__ import annotations

import gzip
import re
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
def _default_vcf():
    """Descubre el VCF en data/ sin escribir el nombre de la muestra en el codigo."""
    d = ROOT / "data"
    hits = sorted(p for p in d.glob("*.vcf.gz")) if d.exists() else []
    return hits[0] if hits else d / "input.vcf.gz"


VCF = _default_vcf()
CLINVAR = ROOT / "ref" / "clinvar.vcf.gz"
OUT = ROOT / "work" / "clinvar_genomewide.tsv"

GENEINFO = re.compile(r"(?:^|;)GENEINFO=([^;:]+)")
CLNSIG = re.compile(r"(?:^|;)CLNSIG=([^;]+)")
CLNDN = re.compile(r"(?:^|;)CLNDN=([^;]+)")
CLNREV = re.compile(r"(?:^|;)CLNREVSTAT=([^;]+)")

# Solo clasificaciones que aguantan peso. Conflicting y VUS quedan fuera.
BUENAS = ("Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic")
# Revisiones debiles: una sola submision sin criterios no sostiene un hallazgo.
REV_DEBIL = ("no_assertion", "no_classification", "no_interpretation")


def load_clinvar():
    """{(chrom,pos,ref,alt): (gene, sig, rev, enfermedad)} solo P/LP."""
    table = {}
    with gzip.open(CLINVAR, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            info = f[7]
            m = CLNSIG.search(info)
            if not m:
                continue
            sig = m.group(1)
            if not sig.startswith(BUENAS):
                continue
            gene = GENEINFO.search(info)
            dn = CLNDN.search(info)
            rev = CLNREV.search(info)
            table[(f[0], int(f[1]), f[3], f[4])] = (
                gene.group(1) if gene else "",
                sig,
                rev.group(1) if rev else "",
                dn.group(1) if dn else "")
    return table


def main() -> int:
    print("cargando ClinVar...", flush=True)
    cv = load_clinvar()
    print("%d variantes P/LP en ClinVar" % len(cv), flush=True)

    hits = []
    with gzip.open(VCF, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            key = (f[0], int(f[1]), f[3], f[4])
            rec = cv.get(key)
            if rec is None:
                continue
            kv = dict(zip(f[8].split(":"), f[9].split(":")))
            gt = kv.get("GT", "./.").replace("|", "/")
            a = gt.split("/")
            if len(a) != 2 or "." in a or a == ["0", "0"]:
                continue
            zyg = "hom" if a[0] == a[1] else "het"
            hits.append({"chrom": f[0], "pos": f[1], "ref": f[3], "alt": f[4],
                         "filter": f[6], "gt": gt, "zyg": zyg,
                         "dp": kv.get("DP", ""), "gq": kv.get("GQ", ""),
                         "gene": rec[0], "sig": rec[1], "rev": rec[2], "dn": rec[3]})

    with OUT.open("w", encoding="utf-8") as out:
        cols = ["chrom", "pos", "ref", "alt", "filter", "gt", "zyg", "dp", "gq",
                "gene", "sig", "rev", "dn"]
        out.write("\t".join(cols) + "\n")
        for h in hits:
            out.write("\t".join(str(h[c]) for c in cols) + "\n")

    fuertes = [h for h in hits if h["filter"] in ("PASS", ".")
               and not h["rev"].startswith(REV_DEBIL)
               and (h["gq"].isdigit() and int(h["gq"]) >= 50)]

    print("\n%d genotipos no de referencia coinciden con una P/LP de ClinVar" % len(hits))
    print("%d de ellos pasan filtro, GQ>=50 y revision con criterios\n" % len(fuertes))

    por_gen = defaultdict(list)
    for h in fuertes:
        por_gen[h["gene"]].append(h)

    print("== genes con evidencia P/LP, ordenados por cigosidad y numero ==")
    print("%-12s %4s %4s %5s  %s" % ("gen", "hom", "het", "DPmed", "enfermedad"))
    def orden(g):
        hs = por_gen[g]
        return (-sum(1 for h in hs if h["zyg"] == "hom"), -len(hs), g)
    for gene in sorted(por_gen, key=orden)[:40]:
        hs = por_gen[gene]
        hom = sum(1 for h in hs if h["zyg"] == "hom")
        het = len(hs) - hom
        dps = [int(h["dp"]) for h in hs if h["dp"].isdigit()]
        dn = Counter(h["dn"].split("|")[0] for h in hs).most_common(1)[0][0]
        print("%-12s %4d %4d %5.0f  %s" % (gene, hom, het,
                                           sum(dps) / len(dps) if dps else 0, dn[:52]))

    print("\n== lectura ==")
    print("Un heterocigoto P/LP suelto en un gen recesivo es estado de portador, no")
    print("diagnostico: la poblacion general carga varios. Lo que importa bajo modelo")
    print("recesivo es un homocigoto, o un heterocigoto acompanado de una segunda")
    print("variante rara en el mismo gen, que es lo que el barrido del panel ya evaluo.")
    print("\ndetalle en %s" % OUT.relative_to(ROOT))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
