"""Tests del formateador de submission. Datos sinteticos: ninguna coordenada real.

Cubre lo que puede costar una submission: la conversion de contigs, los limites
del formato, el manejo del par heterocigoto compuesto, y la replica del scoring
del hackathon (rank points y F-max).
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from make_submission import build, fmax, norm_contig, rank_points, score, validate  # noqa: E402

FAILS = []


def check(cond, msg):
    if not cond:
        FAILS.append(msg)


# --- conversion de contigs -------------------------------------------------
check(norm_contig("15") == "chr15", "15 debe volverse chr15")
check(norm_contig("chr15") == "chr15", "chr15 se mantiene")
check(norm_contig("X") == "chrX", "X debe volverse chrX")
check(norm_contig("MT") == "chrM", "MT debe volverse chrM")
try:
    norm_contig("GL000220.1")
    FAILS.append("un contig no canonico debe reventar")
except ValueError:
    pass

# --- tabla de puntos -------------------------------------------------------
check(rank_points(1) == 100, "rank 1 son 100")
check(rank_points(2) == 50 and rank_points(3) == 50, "ranks 2-3 son 50")
check(rank_points(4) == 25 and rank_points(5) == 25, "ranks 4-5 son 25")
check(rank_points(6) == 10 and rank_points(10) == 10, "ranks 6-10 son 10")
check(rank_points(11) == 0, "mas de 10 no puntua")
check(rank_points(1, partial=True) == 50, "coincidencia parcial vale la mitad")

# --- construccion y validacion --------------------------------------------
cands = [
    {"chrom": "9", "pos": "1000000", "ref": "A", "alt": "G", "epcr": "0.9"},
    {"chrom": "9", "pos": "1000500", "ref": "C", "alt": "T", "epcr": "0.4"},
]
rows, problems = build(cands, "SYN01")
check(not [p for p in problems if not p.startswith("AVISO")], "candidatos limpios no dan problemas: %s" % problems)
check(rows[0]["epcr"] == 0.9, "las filas se ordenan por epcr descendente")
check(rows[0]["chrom_1"] == "chr9", "el contig se normaliza al construir")
check(rows[0]["finding_type"] == "primary", "finding_type por defecto es primary")

bad = [{"chrom": "1", "pos": "100", "ref": "A", "alt": "G", "epcr": "0"}]
_, problems = build(bad, "SYN01")
check(any("fuera de (0,1]" in p for p in problems), "epcr 0 debe rechazarse")

too_many = [{"chrom": "1", "pos": str(1000 + i), "ref": "A", "alt": "G",
             "epcr": str(0.9 - i / 100)} for i in range(11)]
_, problems = build(too_many, "SYN01")
check(any("maximo es 10" in p for p in problems), "11 filas deben rechazarse")

ties = [{"chrom": "1", "pos": "1000", "ref": "A", "alt": "G", "epcr": "0.5"},
        {"chrom": "2", "pos": "2000", "ref": "C", "alt": "T", "epcr": "0.5"}]
_, problems = build(ties, "SYN01")
check(any(p.startswith("AVISO") and "repetido" in p for p in problems),
      "los epcr repetidos deben avisar")

half = [{"chrom_1": "1", "pos_1": "1000", "ref_1": "A", "alt_1": "G",
         "chrom_2": "1", "epcr": "0.8"}]
_, problems = build(half, "SYN01")
check(any("incompleta" in p for p in problems), "media segunda variante debe rechazarse")

# --- scoring: compuesto completo vs parcial --------------------------------
truth = [("chr15", 40160000, "A", "G"), ("chr15", 40170000, "C", "T")]

full = [{"chrom_1": "chr15", "pos_1": 40160000, "ref_1": "A", "alt_1": "G",
         "chrom_2": "chr15", "pos_2": 40170000, "ref_2": "C", "alt_2": "T", "epcr": 0.95}]
s = score(full, truth)
check(s["rank_points"] == 100, "compuesto completo en rank 1 son 100 puntos")
check(abs(s["fmax"] - 1.0) < 1e-9, "compuesto completo da F-max 1.0")

partial = [{"chrom_1": "chr15", "pos_1": 40160000, "ref_1": "A", "alt_1": "G",
            "chrom_2": "", "epcr": 0.95}]
s = score(partial, truth)
check(s["rank_points"] == 50, "una sola variante del par vale la mitad")
check(abs(s["fmax"] - (2 * 1.0 * 0.5) / 1.5) < 1e-9, "F-max de recuperar una de dos es 2/3")

# --- la estrategia de llenar las 10 filas no daña si el orden es correcto ---
padded = [dict(full[0])] + [
    {"chrom_1": "chr%d" % (i + 1), "pos_1": 1000000 + i, "ref_1": "A", "alt_1": "G",
     "chrom_2": "", "epcr": 0.5 - i / 100} for i in range(9)]
s_padded = score(padded, truth)
check(s_padded["rank_points"] == 100, "el relleno por debajo no baja los rank points")
check(abs(s_padded["fmax"] - 1.0) < 1e-9, "el relleno por debajo no baja F-max")

# --- pero una equivocada por encima sí duele -------------------------------
wrong_first = [
    {"chrom_1": "chr1", "pos_1": 999999, "ref_1": "A", "alt_1": "G", "chrom_2": "", "epcr": 0.99},
    dict(full[0], epcr=0.90),
]
s_wrong = score(wrong_first, truth)
check(s_wrong["rank_points"] == 50, "la verdadera en rank 2 son 50 puntos")
check(s_wrong["fmax"] < 1.0, "una equivocada por encima baja F-max")

print("tests de submission: %d fallos" % len(FAILS))
for f in FAILS:
    print("  !", f)
sys.exit(1 if FAILS else 0)
