#!/usr/bin/env python3
"""Genera y valida el CSV de submission del Track 1.

El VCF entregado usa contigs sin prefijo (15) y el formulario pide chr15. Esa
conversion es el error tonto que cuesta una submission entera, asi que vive aqui,
con test.

Incluye una replica del scoring del hackathon (rank points y F-max, tal como los
calcula evaluation.py del Space) para poder calibrar el orden y la separacion de
los EPCR contra un set sintetico antes de gastar un envio.

Uso:
    python scripts/make_submission.py candidatos.tsv salida.csv --proband ID
    python scripts/make_submission.py --validate salida.csv
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

COLUMNS = ["proband_id", "chrom_1", "pos_1", "ref_1", "alt_1",
           "chrom_2", "pos_2", "ref_2", "alt_2", "epcr", "finding_type", "notes"]
MAX_ROWS = 10
BASES = re.compile(r"^[ACGTN]+$", re.IGNORECASE)
RANK_POINTS = [(1, 100), (3, 50), (5, 25), (10, 10)]


def norm_contig(c: str) -> str:
    """15 -> chr15, MT -> chrM. El formulario pide el estilo con prefijo."""
    c = str(c).strip()
    if not c:
        return ""
    c = c[3:] if c.lower().startswith("chr") else c
    c = c.upper() if c.upper() in ("X", "Y", "M", "MT") else c
    if c in ("MT", "M"):
        c = "M"
    if not re.fullmatch(r"([1-9]|1\d|2[0-2]|X|Y|M)", c):
        raise ValueError("contig no canonico: %r" % c)
    return "chr" + c


def rank_points(rank: int, partial: bool = False) -> int:
    pts = 0
    for limit, value in RANK_POINTS:
        if rank <= limit:
            pts = value
            break
    return pts // 2 if partial else pts


def fmax(rows, truth):
    """F-max a nivel de variante individual, barriendo los umbrales de EPCR."""
    truth = set(truth)
    best = (0.0, None)
    for t in sorted({r["epcr"] for r in rows}, reverse=True):
        called = set()
        for r in rows:
            if r["epcr"] >= t:
                called.add((r["chrom_1"], int(r["pos_1"]), r["ref_1"], r["alt_1"]))
                if r.get("chrom_2"):
                    called.add((r["chrom_2"], int(r["pos_2"]), r["ref_2"], r["alt_2"]))
        tp = len(called & truth)
        fp = len(called - truth)
        fn = len(truth - called)
        if tp == 0:
            continue
        prec = tp / (tp + fp)
        rec = tp / (tp + fn)
        f = 2 * prec * rec / (prec + rec)
        if f > best[0]:
            best = (f, t)
    return best


def score(rows, truth_rows):
    """truth_rows: lista de tuplas de variantes; un compuesto son dos tuplas."""
    rows = sorted(rows, key=lambda r: -r["epcr"])
    truth = {tuple(v) for v in truth_rows}
    points = 0
    for i, r in enumerate(rows, start=1):
        called = {(r["chrom_1"], int(r["pos_1"]), r["ref_1"], r["alt_1"])}
        if r.get("chrom_2"):
            called.add((r["chrom_2"], int(r["pos_2"]), r["ref_2"], r["alt_2"]))
        if called == truth:
            points = max(points, rank_points(i))
        elif called & truth:
            points = max(points, rank_points(i, partial=True))
    f, thr = fmax(rows, truth)
    return {"rank_points": points, "fmax": f, "fmax_threshold": thr}


def read_candidates(path: Path):
    with path.open(newline="", encoding="utf-8") as fh:
        sniff = csv.Sniffer().sniff(fh.read(4096), delimiters="\t,")
        fh.seek(0)
        return list(csv.DictReader(fh, dialect=sniff))


def build(candidates, proband: str):
    rows, problems = [], []
    for i, c in enumerate(candidates, start=1):
        try:
            row = {
                "proband_id": proband,
                "chrom_1": norm_contig(c.get("chrom_1") or c.get("chrom") or ""),
                "pos_1": int(str(c.get("pos_1") or c.get("pos")).strip()),
                "ref_1": (c.get("ref_1") or c.get("ref") or "").strip().upper(),
                "alt_1": (c.get("alt_1") or c.get("alt") or "").strip().upper(),
                "chrom_2": norm_contig(c["chrom_2"]) if c.get("chrom_2") else "",
                "pos_2": int(c["pos_2"]) if c.get("pos_2") else "",
                "ref_2": (c.get("ref_2") or "").strip().upper(),
                "alt_2": (c.get("alt_2") or "").strip().upper(),
                "epcr": float(c["epcr"]),
                "finding_type": (c.get("finding_type") or "primary").strip().lower(),
                "notes": (c.get("notes") or "").strip(),
            }
        except (ValueError, KeyError, TypeError) as e:
            problems.append("fila %d: %s" % (i, e))
            continue
        rows.append(row)
    rows.sort(key=lambda r: -r["epcr"])
    problems.extend(validate(rows))
    return rows, problems


def validate(rows):
    problems = []
    if len(rows) > MAX_ROWS:
        problems.append("%d filas: el maximo es %d" % (len(rows), MAX_ROWS))
    seen_epcr = {}
    for i, r in enumerate(rows, start=1):
        if not r["chrom_1"] or not r["pos_1"]:
            problems.append("fila %d sin variante primaria" % i)
        if not (0 < r["epcr"] <= 1):
            problems.append("fila %d: epcr %s fuera de (0,1]" % (i, r["epcr"]))
        if r["finding_type"] not in ("primary", "secondary"):
            problems.append("fila %d: finding_type invalido %r" % (i, r["finding_type"]))
        for tag in ("1", "2"):
            ref, alt = r["ref_" + tag], r["alt_" + tag]
            if ref and not BASES.match(ref):
                problems.append("fila %d: ref_%s invalido %r" % (i, tag, ref))
            if alt and not BASES.match(alt):
                problems.append("fila %d: alt_%s invalido %r" % (i, tag, alt))
        half = [bool(r["chrom_2"]), bool(r["pos_2"]), bool(r["ref_2"]), bool(r["alt_2"])]
        if any(half) and not all(half):
            problems.append("fila %d: segunda variante incompleta" % i)
        seen_epcr.setdefault(r["epcr"], []).append(i)
    for value, idxs in seen_epcr.items():
        if len(idxs) > 1:
            problems.append("AVISO: epcr %s repetido en las filas %s; el desempate es "
                            "por orden de envio y F-max pierde resolucion" % (value, idxs))
    return problems


def write(rows, out: Path):
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("input")
    ap.add_argument("output", nargs="?")
    ap.add_argument("--proband", default="")
    ap.add_argument("--validate", action="store_true")
    args = ap.parse_args()

    if args.validate:
        rows = []
        with open(args.input, newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                r["epcr"] = float(r["epcr"])
                rows.append(r)
        problems = validate(rows)
        print("%d filas" % len(rows))
        for p in problems:
            print("  !", p)
        print("OK" if not problems else "revisa lo de arriba")
        return 1 if any(not p.startswith("AVISO") for p in problems) else 0

    if not args.output:
        print("falta la ruta de salida")
        return 1
    rows, problems = build(read_candidates(Path(args.input)), args.proband)
    for p in problems:
        print("  !", p)
    blocking = [p for p in problems if not p.startswith("AVISO")]
    if blocking:
        print("no se escribe nada hasta que eso este limpio")
        return 1
    write(rows, Path(args.output))
    print("escritas %d filas en %s" % (len(rows), args.output))
    print("recordatorio: ninguna submission se envia sin confirmacion humana")
    return 0


if __name__ == "__main__":
    sys.exit(main())
