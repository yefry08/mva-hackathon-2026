#!/usr/bin/env python3
"""Purga de datos del paciente, exigida por las reglas del hackathon.

Las reglas obligan a borrar, dentro de los 30 dias del cierre (2026-10-24), todo
lo que lleve el genoma del nino: VCF/BAM/CRAM/FASTQ y cualquier copia o recorte,
tablas intermedias con genotipos a escala genomica, caches, y los prompts o logs
guardados en nuestros sistemas que contengan bloques de datos de variantes.

Lo que sobrevive: la lista entregada, los HPO, rankings de genes y vias, el
mecanismo, los candidatos de farmaco, el codigo y el reporte.

Uso:
    python scripts/purge.py            # simulacro: dice que borraria
    python scripts/purge.py --apply    # borra de verdad
    python scripts/purge.py --verify   # revisa que lo que queda no lleve genotipos
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HOME = Path.home()

# Directorios que se borran enteros.
PURGE_DIRS = [
    ROOT / "data",
    ROOT / "work",
    ROOT / "logs",
    ROOT / ".snakemake",
    # Transcripts locales de Claude Code de este proyecto: son "logs guardados en
    # nuestros sistemas con bloques de datos de variantes" segun la regla.
    HOME / ".claude" / "projects" / "C--mva",
]

# Extensiones que se borran esten donde esten.
PURGE_GLOBS = [
    "**/*.vcf", "**/*.vcf.gz", "**/*.vcf.gz.tbi", "**/*.bcf",
    "**/*.bam", "**/*.bai", "**/*.cram", "**/*.crai",
    "**/*.fastq", "**/*.fastq.gz", "**/*.fq.gz",
    "**/*.g.vcf.gz", "**/annotated.parquet",
]

# No se tocan.
KEEP = [
    ROOT / "submissions",
    ROOT / "reports",
    ROOT / "track2",
    ROOT / "pipelines",
    ROOT / "scripts",
    ROOT / ".claude" / "skills",
]

# Firma de contenido a escala de genotipo, para --verify.
GENOTYPE_SIGNS = [
    re.compile(r"\b(?:chr)?(?:[1-9]|1\d|2[0-2]|X|Y|MT)[\s\t]+\d{4,9}[\s\t]+\S+[\s\t]+[ACGTN]+[\s\t]+[ACGTN,]+"),
    re.compile(r"\b\d[/|]\d:\d+,\d+:\d+\b"),
]
VERIFY_EXT = {".tsv", ".csv", ".txt", ".json", ".jsonl", ".md", ".parquet"}
# Un pu\u00f1ado de variantes nombradas es un hallazgo; una tabla genomica es el dataset.
VERIFY_LINE_THRESHOLD = 25


def human(n: int) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return "%.1f %s" % (n, unit)
        n /= 1024.0
    return str(n)


def dir_size(p: Path) -> int:
    total = 0
    for f in p.rglob("*"):
        if f.is_file():
            try:
                total += f.stat().st_size
            except OSError:
                pass
    return total


def collect():
    targets = []
    for d in PURGE_DIRS:
        if d.exists():
            targets.append((d, dir_size(d)))
    seen = {t[0] for t in targets}
    for pattern in PURGE_GLOBS:
        for f in ROOT.glob(pattern):
            if not f.is_file():
                continue
            if any(str(f).startswith(str(d)) for d in seen):
                continue
            targets.append((f, f.stat().st_size))
    return targets


def verify() -> int:
    """Revisa que ningun artefacto conservado lleve genotipos a escala genomica."""
    offenders = []
    for keep in KEEP:
        if not keep.exists():
            continue
        for f in keep.rglob("*"):
            if not f.is_file() or f.suffix.lower() not in VERIFY_EXT:
                continue
            try:
                text = f.read_text(encoding="utf-8", errors="ignore")
            except OSError:
                continue
            hits = sum(1 for line in text.splitlines()
                       if any(s.search(line) for s in GENOTYPE_SIGNS))
            if hits >= VERIFY_LINE_THRESHOLD:
                offenders.append((f, hits))
    if offenders:
        print("ARTEFACTOS CONSERVADOS CON ASPECTO DE DATASET:")
        for f, hits in offenders:
            print("  %s (%d lineas con genotipos)" % (f.relative_to(ROOT), hits))
        print("\nRevisalos a mano: un pu\u00f1ado de variantes es un hallazgo, una tabla no.")
        return 1
    print("verify OK: nada conservado parece una tabla de genotipos.")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true", help="borra de verdad")
    ap.add_argument("--verify", action="store_true", help="solo revisa lo conservado")
    args = ap.parse_args()

    if args.verify:
        return verify()

    targets = collect()
    if not targets:
        print("nada que purgar.")
        return 0

    total = sum(s for _, s in targets)
    print("%s a borrar en %d objetivos:\n" % (human(total), len(targets)))
    for p, s in targets:
        try:
            shown = p.relative_to(ROOT)
        except ValueError:
            shown = p
        print("  %-10s %s" % (human(s), shown))

    if not args.apply:
        print("\nsimulacro. Para borrar de verdad: python scripts/purge.py --apply")
        return 0

    for p, _ in targets:
        if p.is_dir():
            shutil.rmtree(p, ignore_errors=True)
        else:
            try:
                p.unlink()
            except OSError as e:
                print("  no se pudo borrar %s: %s" % (p, e))
    print("\npurga aplicada. Ahora:")
    print("  1. python scripts/purge.py --verify")
    print("  2. revisa copias fuera de este arbol (nube, otras maquinas, papelera)")
    print("  3. envia la atestacion a RarediseaserealkidMVAhackathon2026@synapse.org")
    return 0


if __name__ == "__main__":
    sys.exit(main())
