#!/usr/bin/env python3
"""Cribado de sitios de splicing crípticos, con contexto y con nulo empírico.

La primera version de este script contaba como hallazgo cualquier dinucleotido
creado de una lista de cinco. Eso da positivo en mas de la mitad de las variantes
al azar, porque un cambio de base toca dos dinucleotidos y cinco de dieciseis
estan en la lista. Un test que se dispara siempre no informa nada, asi que se
reemplazo por este.

Que mide ahora. Un sitio de splicing no es un dinucleotido suelto, es un motivo
con contexto:

  donador   MAG | GTRAGT   el GT canonico mas el consenso a ambos lados
  aceptor   tracto de polipirimidinas ... AG | G   el AG precedido de pirimidinas

Asi que solo cuenta como sitio criptico plausible un GT o un AG creado que ademas
encaje razonablemente con su consenso. Y para saber si ese "razonablemente" vale
algo, el script estima el nulo: muestrea posiciones al azar del mismo intron,
aplica el mismo tipo de sustitucion, y mide cada cuanto aparece un sitio igual de
plausible por puro azar. Si la variante real no destaca contra ese fondo, no hay
senal.

Esto sigue sin ser SpliceAI. Es un cribado con su tasa de falsos positivos medida,
que es distinto de un cribado sin calibrar.

Uso: python scripts/splice_screen.py
"""
from __future__ import annotations

import csv
import gzip
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FASTA = ROOT / "ref" / "chr15.fa.gz"
CAND = ROOT / "work" / "candidates.clinvar.tsv"

FLANK = 30
N_NULL = 2000          # posiciones al azar para estimar el fondo
DONOR_MIN = 6          # de 8 posiciones de consenso
ACCEPTOR_PY_MIN = 0.65  # fraccion de pirimidinas en las 20 nt previas

# Consenso del donador: posiciones -3..+6 respecto del corte, con el GT fijo.
# MAG|GTRAGT -> se puntuan las posiciones degeneradas.
DONOR_CONSENSUS = [("M", "AC"), ("A", "A"), ("G", "G"),
                   ("G", "G"), ("T", "T"),
                   ("R", "AG"), ("A", "A"), ("G", "G"), ("T", "T")]


def load_seq():
    parts = []
    with gzip.open(FASTA, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.startswith(">"):
                parts.append(line.strip())
    return "".join(parts).upper()


def donor_score(win, i):
    """win[i:i+2] debe ser GT. Puntua el consenso MAG|GTRAGT alrededor."""
    start = i - 3
    if start < 0 or start + 9 > len(win):
        return -1
    s = win[start:start + 9]
    return sum(1 for base, (_, allowed) in zip(s, DONOR_CONSENSUS) if base in allowed)


def acceptor_ok(win, i):
    """win[i:i+2] debe ser AG. Exige tracto de polipirimidinas por delante."""
    up = win[max(0, i - 20):i]
    if len(up) < 12:
        return False, 0.0
    py = sum(1 for b in up if b in "CT") / len(up)
    return py >= ACCEPTOR_PY_MIN, py


def new_sites(ref_win, alt_win):
    """Sitios plausibles presentes en alt y ausentes en ref."""
    out = []
    for i in range(len(alt_win) - 1):
        d = alt_win[i:i + 2]
        if ref_win[i:i + 2] == d:
            continue
        if d == "GT":
            sc = donor_score(alt_win, i)
            if sc >= DONOR_MIN:
                out.append(("donador", i, "consenso %d/9" % sc))
        elif d == "AG":
            ok, py = acceptor_ok(alt_win, i)
            if ok:
                out.append(("aceptor", i, "polipirimidinas %.0f%%" % (py * 100)))
    return out


def null_rate(seq, center, alt_base, rng):
    """Cada cuanto una sustitucion al azar del entorno crea un sitio plausible."""
    hits = 0
    for _ in range(N_NULL):
        p = center + rng.randint(-5000, 5000)
        if p - FLANK < 0 or p + FLANK + 1 > len(seq):
            continue
        ref_win = seq[p - FLANK:p + FLANK + 1]
        if "N" in ref_win or ref_win[FLANK] == alt_base:
            continue
        alt_win = ref_win[:FLANK] + alt_base + ref_win[FLANK + 1:]
        if new_sites(ref_win, alt_win):
            hits += 1
    return 100.0 * hits / N_NULL


def main() -> int:
    rows = [r for r in csv.DictReader(CAND.open(encoding="utf-8"), delimiter="\t")
            if r["gene"] == "BUB1B"]
    if not rows:
        print("no hay candidatas de BUB1B")
        return 1

    print("cargando la referencia de chr15...", flush=True)
    seq = load_seq()
    rng = random.Random(20260912)
    print("%d pb cargadas\n" % len(seq))

    for r in sorted(rows, key=lambda r: r["region"]):
        pos, ref, alt = int(r["pos"]), r["ref"].upper(), r["alt"].upper()
        i = pos - 1
        if seq[i] != ref:
            print("%-30s DESAJUSTE con la referencia" % r["region"][:30])
            continue
        ref_win = seq[i - FLANK:i + FLANK + 1]
        alt_win = ref_win[:FLANK] + alt + ref_win[FLANK + 1:]
        sitios = new_sites(ref_win, alt_win)
        fondo = null_rate(seq, i, alt, rng)

        print("== %s" % r["region"])
        print("   %s>%s, clase %s" % (ref, alt, r["clase"]))
        if sitios:
            for tipo, j, det in sitios:
                print("   sitio criptico de %s a %+d pb (%s)" % (tipo, j - FLANK, det))
            print("   fondo: %.1f%% de las sustituciones %s al azar del entorno crean uno igual"
                  % (fondo, alt))
            print("   lectura: %s" % ("destaca sobre el fondo" if fondo < 2 else
                                      "dentro de lo esperable por azar"))
        else:
            print("   no crea ningun sitio plausible")
            print("   fondo: %.1f%% de las sustituciones al azar si lo crean" % fondo)
        print()

    print("Este cribado mide dinucleotido mas consenso, con su tasa de fondo. No")
    print("puntua fuerza relativa de sitios competidores ni elementos reguladores,")
    print("y no sustituye a SpliceAI ni a Pangolin, que siguen pendientes.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
