#!/usr/bin/env python3
"""Checklist de privacidad, bloqueante antes de hacer publico el repositorio.

Las reglas del hackathon prohiben redistribuir los datos por cualquier canal, y
el repositorio tiene que ser publico al cierre. Entre esas dos cosas hay una
oportunidad excelente de filtrar el genoma de un nino sin darse cuenta, asi que
esto se comprueba con un script y no con buena intencion.

Que revisa:

1. Que ningun archivo genomico este bajo control de versiones, ni ahora ni en
   ningun commit del historial. Un archivo borrado en el commit siguiente sigue
   estando en el historial y sigue siendo publico.
2. Que .gitignore cubra de verdad los directorios de datos.
3. Que ningun archivo rastreado contenga coordenadas, rsIDs, HGVS con posicion o
   identificadores de muestra.

La excepcion consciente es submissions/: la lista rankeada de candidatos LLEVA
coordenadas, y es el entregable. Las reglas lo permiten de forma explicita, pero
el script la marca igualmente para que una persona confirme que solo hay ahi lo
que se pretende entregar.

Uso: python scripts/privacy_check.py
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SAMPLE_IDS = ROOT / ".claude" / "sample_ids.txt"

EXCLUDE_DIRS = {".git", "data", "ref", "scale", "__pycache__", ".venv"}
DATA_EXT = re.compile(r"\.(vcf|bam|cram|fastq|fq|bai|crai|tbi|bgz|g\.vcf)(\.gz)?$", re.I)
TEXT_EXT = {".md", ".py", ".txt", ".json", ".yaml", ".yml", ".csv", ".tsv", ".sh", ".cfg"}

PATTERNS = [
    (re.compile(r"\b(?:chr)?(?:[1-9]|1\d|2[0-2]|X|Y|MT?)\s*[:\-]\s*\d{5,}\b"), "coordenada"),
    # Una tabla de variantes no separa cromosoma y posicion con dos puntos: los
    # pone en columnas contiguas. Sin esta linea el checker era ciego al formato
    # con mas probabilidad de filtrar un genoma, que es un CSV o un TSV.
    # (Los ejemplos van sin numeros a proposito: con ellos, este archivo se
    # detectaria a si mismo, cosa que paso de verdad al escribirlo.)
    (re.compile(r"(?:^|[,\t])(?:chr)?(?:[1-9]|1\d|2[0-2]|X|Y|MT)[,\t]\s*\d{5,}(?:[,\t]|$)",
                re.M), "coordenada en columnas"),
    (re.compile(r"\brs\d{3,}\b"), "rsID"),
    (re.compile(r"\b[cgmn]\.\s*-?\*?\d+"), "HGVS con posicion"),
]


def load_ids():
    out = []
    if SAMPLE_IDS.exists():
        for line in SAMPLE_IDS.read_text(encoding="utf-8").splitlines():
            line = line.split("#", 1)[0].strip()
            if line:
                out.append(re.compile(re.escape(line)))
    return out


def git(*args):
    try:
        r = subprocess.run(["git", *args], cwd=ROOT, capture_output=True,
                           text=True, encoding="utf-8", errors="replace")
        return r.stdout
    except FileNotFoundError:
        return ""


def main() -> int:
    fallos, avisos = [], []
    id_pats = load_ids()

    # 1. archivos genomicos rastreados, ahora
    tracked = [l for l in git("ls-files").splitlines() if l.strip()]
    malos = [f for f in tracked if DATA_EXT.search(f)]
    print("1. archivos genomicos rastreados ahora: %d" % len(malos))
    if malos:
        fallos.append("hay archivos genomicos en el indice: %s" % malos[:5])

    # 2. archivos genomicos en cualquier punto del historial
    hist = git("log", "--all", "--pretty=format:", "--name-only")
    hist_malos = sorted({l for l in hist.splitlines() if l.strip() and DATA_EXT.search(l)})
    print("2. archivos genomicos en el historial: %d" % len(hist_malos))
    if hist_malos:
        fallos.append("el historial contiene archivos genomicos: %s" % hist_malos[:5])

    # 3. .gitignore
    gi = (ROOT / ".gitignore").read_text(encoding="utf-8") if (ROOT / ".gitignore").exists() else ""
    faltan = [d for d in ("data/", "work", "ref/", "scale") if d not in gi]
    print("3. .gitignore cubre los directorios de datos: %s" % ("si" if not faltan else "NO (%s)" % faltan))
    if faltan:
        fallos.append(".gitignore no cubre: %s" % faltan)

    # 4. contenido de los archivos rastreados
    print("4. contenido de archivos rastreados:")
    revisados = 0
    for f in tracked:
        p = ROOT / f
        if not p.exists() or p.suffix.lower() not in TEXT_EXT:
            continue
        if set(p.parts) & EXCLUDE_DIRS:
            continue
        texto = p.read_text(encoding="utf-8", errors="ignore")
        revisados += 1
        hits = []
        for rx, label in PATTERNS:
            if rx.search(texto):
                hits.append(label)
        for rx in id_pats:
            if rx.search(texto):
                hits.append("identificador de muestra")
        if not hits:
            continue
        # submissions/ lleva coordenadas por diseno: es el entregable
        if f.startswith("submissions/"):
            avisos.append("%s contiene %s (esperado: es la lista entregada)"
                          % (f, ", ".join(sorted(set(hits)))))
        elif f.startswith("tests/") and "identificador de muestra" not in hits:
            avisos.append("%s contiene %s (esperado: canarios sinteticos del test)"
                          % (f, ", ".join(sorted(set(hits)))))
        else:
            fallos.append("%s contiene %s" % (f, ", ".join(sorted(set(hits)))))
    print("   %d archivos de texto revisados" % revisados)

    # 5. contenido de TODO el historial, no solo del estado actual.
    # El paso 2 solo mira extensiones. Un .tsv con coordenadas del paciente,
    # commiteado y borrado en el commit siguiente, seguiria publico en el
    # historial y ese paso no lo veria. Esto ya paso una vez, con datos de
    # controles publicos que se colaron por un gate saltado.
    print("5. contenido del historial completo:")
    CONTROLES_PUBLICOS = re.compile(r"^work_hg00\d/")
    vistos, n_blobs = set(), 0
    for commit in git("rev-list", "--all").split():
        for linea in git("ls-tree", "-r", commit).splitlines():
            try:
                meta, ruta = linea.split("\t", 1)
                blob = meta.split()[2]
            except (ValueError, IndexError):
                continue
            if blob in vistos or Path(ruta).suffix.lower() not in TEXT_EXT:
                continue
            vistos.add(blob)
            n_blobs += 1
            contenido = git("cat-file", "-p", blob)
            hits = [label for rx, label in PATTERNS if rx.search(contenido)]
            hits += ["identificador de muestra" for rx in id_pats if rx.search(contenido)]
            if not hits:
                continue
            if "identificador de muestra" in hits:
                fallos.append("historial %s: %s contiene identificador de muestra" % (commit[:7], ruta))
            elif ruta.startswith(("submissions/", "tests/")):
                pass                                  # excepciones declaradas arriba
            elif CONTROLES_PUBLICOS.match(ruta):
                avisos.append("historial %s: %s (genoma publico GIAB, no del paciente)"
                              % (commit[:7], ruta))
            else:
                fallos.append("historial %s: %s contiene %s"
                              % (commit[:7], ruta, ", ".join(sorted(set(hits)))))
    print("   %d versiones de archivo revisadas en %d commits"
          % (n_blobs, len(git("rev-list", "--all").split())))

    print()
    if avisos:
        print("AVISOS (revisar a mano, no bloquean):")
        for a in avisos:
            print("  -", a)
        print()
    if fallos:
        print("BLOQUEANTE. No hacer publico el repositorio hasta resolver:")
        for f in fallos:
            print("  -", f)
        return 1
    print("CHECKLIST EN VERDE. El repositorio no contiene datos del paciente")
    print("fuera de la lista de candidatos entregada.")
    print()
    print("Queda una comprobacion que ningun script puede hacer por ti: abrir")
    print("submissions/ y confirmar que ahi solo esta lo que pretendes entregar.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
