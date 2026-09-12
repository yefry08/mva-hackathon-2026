#!/usr/bin/env python3
"""Egress guard: PreToolUse hook.

Bloquea cualquier llamada de herramienta que sacaria datos a nivel de paciente de
esta maquina. Regla 3 de CLAUDE.md. Codigo de salida 2 = bloquear; stderr se le
muestra al agente.

Alcance: WebFetch, WebSearch, toda herramienta mcp__*, y comandos Bash que invocan
un binario de red. Las lecturas y escrituras locales no se tocan: esto es sobre
egreso, no sobre si un agente puede mirar los datos.
"""
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

# Identificadores especificos del caso, uno por linea, comentarios con #.
SAMPLE_ID_FILE = ROOT / ".claude" / "sample_ids.txt"

PATTERNS = [
    (r"\b(?:chr)?(?:[1-9]|1\d|2[0-2]|X|Y|MT?)\s*[:\-]\s*\d{5,}\b", "coordenada genomica"),
    (r"\brs\d{3,}\b", "rsID"),
    (r"\b[cgmn]\.\s*-?\*?\d+", "HGVS con posicion"),
    (r"\bp\.\s*(?:[A-Z][a-z]{2}|[A-Z])\s*\d+", "HGVS de proteina"),
    (r"\b(?:NM_|NC_|ENST)\d+(?:\.\d+)?\s*:\s*[cgnp]\.", "HGVS con transcrito"),
    (r"\.(?:vcf|bam|cram|fastq|fq|bai|crai|tbi)(?:\.gz)?\b", "nombre de archivo genomico"),
    (r"(?:^|[\s/\"'=@])(?:data|work)/", "ruta de datos del paciente"),
    (r"(?:^|[\s/\"'=@])(?:data|work)\\", "ruta de datos del paciente"),
    (r"\b(?:chr)?(?:[1-9]|1\d|2[0-2]|X|Y|MT)[\s\t]+\d{4,9}[\s\t]+\S+[\s\t]+[ACGTN]+[\s\t]+[ACGTN,]+",
     "registro VCF (CHROM POS ID REF ALT)"),
    (r"\bGT[:=]\s*\d[/|]\d\b", "campo de genotipo"),
    (r"\b\d[/|]\d:\d+,\d+:\d+\b", "registro de genotipo VCF"),
]

NET_BINARIES = re.compile(
    r"\b(curl|wget|nc|ncat|telnet|ssh|scp|sftp|rsync|aws|gcloud|az|gh|hf|"
    r"huggingface-cli|http|httpie)\b"
)
GIT_PUSH = re.compile(r"\bgit\s+(push|remote\s+add)\b")


def load_sample_ids():
    pats = []
    if SAMPLE_ID_FILE.exists():
        for line in SAMPLE_ID_FILE.read_text(encoding="utf-8").splitlines():
            line = line.split("#", 1)[0].strip()
            if line:
                pats.append((re.escape(line), "identificador de muestra"))
    return pats


def is_egress(tool_name, payload):
    if tool_name in ("WebFetch", "WebSearch"):
        return True
    if tool_name.startswith("mcp__"):
        return True
    if tool_name in ("Bash", "PowerShell"):
        return bool(NET_BINARIES.search(payload) or GIT_PUSH.search(payload))
    return False


def main():
    try:
        event = json.load(sys.stdin)
    except Exception:
        return 0  # no rompemos la corrida por un payload raro
    tool_name = event.get("tool_name", "")
    tool_input = event.get("tool_input", {})
    payload = json.dumps(tool_input, ensure_ascii=False)

    if not is_egress(tool_name, payload):
        return 0

    hits = []
    for pattern, label in PATTERNS + load_sample_ids():
        m = re.search(pattern, payload, re.IGNORECASE)
        if m:
            hits.append("%s (%r)" % (label, m.group(0)[:40]))

    if hits:
        sys.stderr.write(
            "EGRESS-GUARD: llamada bloqueada a %s.\n"
            "Detectado: %s\n"
            "Regla 3 de CLAUDE.md: hacia servicios externos solo va nivel gen, via y "
            "mecanismo. Reformula sin coordenadas ni identificadores, o reporta "
            "status blocked. No rodees este hook.\n" % (tool_name, "; ".join(hits))
        )
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
