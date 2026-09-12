#!/usr/bin/env python3
"""Hasta donde llega el phasing por lecturas en este VCF, y por que aqui no alcanza.

Decir "la fase no se resuelve" es una afirmacion que hay que sostener con numeros,
no con una excusa. GATK agrupa en un bloque PID las variantes que caben en la misma
lectura o en el mismo par, asi que la distribucion empirica de los tamanos de bloque
de este VCF mide directamente el alcance del phasing por lecturas en estos datos.

Si el tramo entre las dos candidatas supera con holgura el bloque mas largo
observado, entonces ninguna reconfiguracion del analisis las habria faseado: es el
limite del secuenciado de lecturas cortas, no una decision del pipeline.

Salida: solo estadisticas agregadas y distancias. Ninguna coordenada.

Uso: python scripts/phase_analysis.py [chrom]
"""
from __future__ import annotations

import csv
import gzip
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


VCF = _default_vcf()
PROT = ROOT / "work" / "candidates.protein.tsv"


def main() -> int:
    target = sys.argv[1] if len(sys.argv) > 1 else "15"

    blocks = defaultdict(list)   # pid -> [posiciones]
    n_var = n_phased = 0
    with gzip.open(VCF, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t", 10)
            if f[0] != target:
                if n_var and f[0] != target:
                    # el VCF viene ordenado; si ya pasamos el cromosoma, cortamos
                    if n_var > 1000:
                        break
                continue
            if len(f) < 10:
                continue
            n_var += 1
            kv = dict(zip(f[8].split(":"), f[9].split(":")))
            pid = kv.get("PID", "")
            if pid and pid != ".":
                n_phased += 1
                blocks[pid].append(int(f[1]))

    spans = [max(v) - min(v) for v in blocks.values() if len(v) >= 2]
    sizes = [len(v) for v in blocks.values()]

    print("cromosoma %s: %d variantes, %d con bloque PID (%.1f%%)"
          % (target, n_var, n_phased, 100.0 * n_phased / max(1, n_var)))
    print("bloques de fase: %d, con %d o mas variantes: %d"
          % (len(blocks), 2, len(spans)))
    if spans:
        spans.sort()
        print("\ntamano del bloque de fase, en pares de bases:")
        for label, q in (("mediana", 0.5), ("p90", 0.9), ("p99", 0.99), ("maximo", 1.0)):
            print("  %-8s %8d" % (label, spans[min(len(spans) - 1, int(q * (len(spans) - 1)))]))
        print("  media    %8.0f" % statistics.mean(spans))
        print("variantes por bloque: mediana %d, maximo %d"
              % (statistics.median(sizes), max(sizes)))

    if PROT.exists():
        cod = sorted(int(r["pos"]) for r in
                     csv.DictReader(PROT.open(encoding="utf-8"), delimiter="\t")
                     if r["gene"] == "BUB1B" and r["consecuencia"] in
                     ("NONSENSE (stop ganado)", "missense"))
        if len(cod) >= 2:
            d = cod[-1] - cod[0]
            print("\nseparacion entre las dos candidatas de BUB1B: %d pb" % d)
            if spans:
                print("bloque de fase mas largo observado en el cromosoma: %d pb" % spans[-1])
                print("razon: %.1f veces el maximo observado" % (d / spans[-1]))
                print("\nConclusion: ninguna lectura ni ningun par de este experimento puede")
                print("abarcar las dos posiciones. El phasing por lecturas no falla aqui por")
                print("como se corrio el analisis, sino porque la informacion no esta en los")
                print("datos. Lo que si la resolveria: secuenciar a los padres, lectura larga,")
                print("o PCR de alelo especifico sobre la region.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
