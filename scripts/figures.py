#!/usr/bin/env python3
"""Figuras del reporte, en SVG, con Python de la biblioteca estandar.

Sin matplotlib a proposito: el repositorio entero corre sin dependencias, y las
figuras no son la excepcion. Cada SVG lleva su propia hoja de estilos con modo
claro y oscuro, asi que se ve bien en GitHub con cualquier tema.

Todas las figuras se construyen solo con agregados: conteos, estadisticas por
ventana de 10 Mb y tamanos de bloque de fase. Ninguna posicion del paciente
aparece en ellas.

Figuras:
  fig1_embudo_ciego.svg      del VCF completo a un gen, sin panel ni ClinVar
  fig2_fase.svg              por que la fase no se puede resolver con estos datos
  fig3_baf_genoma.svg        barrido de aneuploidia en mosaico, ventana a ventana
  fig4_controles.svg         el mismo pipeline sobre el paciente y genomas sanos

Uso: python scripts/figures.py
"""
from __future__ import annotations

import csv
import gzip
import json
import math
from collections import defaultdict
from pathlib import Path
from xml.sax.saxutils import escape

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "work"
OUT = ROOT / "reports" / "figures"

# Paleta de referencia validada: azul y naranja son las dos primeras posiciones
# categoricas y pasan las comprobaciones de daltonismo en claro y oscuro. El gris
# NO es un color categorico: marca contexto (controles), y por eso cada control
# lleva su nombre escrito en vez de depender del color.
STYLE = """<style>
  .s  { --surface:#fcfcfb; --ink:#0b0b0b; --ink2:#52514e; --grid:#e5e4df;
        --band:#f3f2ee; --c1:#2a78d6; --c2:#eb6834; --ctx:#8f8e87; }
  @media (prefers-color-scheme: dark) {
    .s { --surface:#1a1a19; --ink:#ffffff; --ink2:#c3c2b7; --grid:#383835;
         --band:#222220; --c1:#3987e5; --c2:#d95926; --ctx:#8f8e87; }
  }
  .bg   { fill: var(--surface); }
  .ink  { fill: var(--ink); }
  .ink2 { fill: var(--ink2); }
  .c1   { fill: var(--c1); }
  .c2   { fill: var(--c2); }
  .ctx  { fill: var(--ctx); }
  .band { fill: var(--band); }
  .grid { stroke: var(--grid); stroke-width: 1; }
  .axis { stroke: var(--ink2); stroke-width: 1; }
  .mark { stroke: var(--ink); stroke-width: 2; }
  .mark2{ stroke: var(--ink2); stroke-width: 1.5; stroke-dasharray: 5 4; }
  .ring { stroke: var(--surface); stroke-width: 2; }
  text  { font-family: system-ui, -apple-system, "Segoe UI", Roboto, sans-serif; }
  .t    { font-size: 17px; font-weight: 600; }
  .st   { font-size: 13px; }
  .lb   { font-size: 12px; }
  .sm   { font-size: 11px; }
  .hero { font-size: 26px; font-weight: 600; font-variant-numeric: tabular-nums; }
  .num  { font-variant-numeric: tabular-nums; }
</style>"""


def fmt(n):
    """5012204 -> 5.012.204, como en los reportes."""
    return "{:,}".format(int(n)).replace(",", ".")


def svg(w, h, title, desc, body):
    return ('<svg xmlns="http://www.w3.org/2000/svg" class="s" width="%d" height="%d" '
            'viewBox="0 0 %d %d" role="img" aria-labelledby="t d">%s'
            '<title id="t">%s</title><desc id="d">%s</desc>'
            '<rect class="bg" width="%d" height="%d" rx="8"/>%s</svg>\n'
            % (w, h, w, h, STYLE, escape(title), escape(desc), w, h, "".join(body)))


def text(x, y, s, cls="lb ink", anchor="start", extra=""):
    return '<text x="%.1f" y="%.1f" class="%s" text-anchor="%s"%s>%s</text>' % (
        x, y, cls, anchor, extra, escape(str(s)))


# --------------------------------------------------------------------------
# Figura 1: embudo ciego
# --------------------------------------------------------------------------
def fig_embudo():
    # Los tres ultimos conteos se leen del resultado; los tres primeros los
    # imprimen lof_genomewide.py y compound_genomewide.py al correr.
    ocho, tres, sostiene = 8, 3, "BUB1B"
    tsv = WORK / "compound_genomewide.tsv"
    if tsv.exists():
        rows = list(csv.DictReader(tsv.open(encoding="utf-8"), delimiter="\t"))
        ocho = len(rows)
        tres = sum(1 for r in rows if r["patron"] == "truncante + missense")

    etapas = [
        ("variantes", 5012204, "variantes en el VCF", None),
        ("variantes", 23041, "variantes en región codificante", "traducción local contra MANE Select"),
        ("genes", 5584, "genes con variantes no sinónimas", None),
        ("genes", 252, "genes con al menos una truncante", "nonsense o cambio de marco de lectura"),
        ("genes", ocho, "genes con patrón recesivo, todo raro", "gnomAD v4.1, frecuencia ≤ 0,001"),
        ("genes", tres, "genes con truncante + missense", "SLC25A5, BUB1B, HLA-DRB1"),
        ("genes", 1, "gen que se sostiene: " + sostiene,
         "SLC25A5 está en el cromosoma X de un varón; HLA-DRB1, en la región HLA"),
    ]
    W, top, row = 760, 96, 62
    H = top + row * len(etapas) + 40 + 22          # 22: aire extra al cambiar de unidad
    b = [text(28, 40, "Del genoma completo a un gen, sin panel y sin ClinVar", "t ink"),
         text(28, 62, "El script no sabe nada de la enfermedad: busca el patrón recesivo en todo el genoma.",
              "st ink2")]
    y = top
    unidad_prev = None
    for i, (unidad, n, etiqueta, nota) in enumerate(etapas):
        ultima = i == len(etapas) - 1
        if unidad != unidad_prev:
            if unidad_prev is not None:
                y += 22                              # sin esto el divisor pisa la nota de arriba
            b.append(text(28, y - 8, unidad.upper(), "sm ink2", extra=' letter-spacing="1"'))
            b.append('<line x1="96" y1="%d" x2="%d" y2="%d" class="grid"/>' % (y - 12, W - 28, y - 12))
            unidad_prev = unidad
        if ultima:
            b.append('<rect x="28" y="%d" width="%d" height="%d" rx="6" class="band"/>'
                     % (y + 2, W - 56, row - 10))
            b.append('<rect x="28" y="%d" width="4" height="%d" rx="2" class="c1"/>' % (y + 2, row - 10))
        b.append(text(200, y + 36, fmt(n), "hero ink" if not ultima else "hero ink", "end"))
        b.append(text(220, y + 30, etiqueta, "st ink" + ("" if not ultima else "")
                      , extra=' font-weight="600"' if ultima else ""))
        if nota:
            b.append(text(220, y + 47, nota, "sm ink2"))
        y += row
    return svg(W, H, "Embudo ciego del genoma completo a un gen",
               "Siete etapas de filtrado, de 5.012.204 variantes a un unico gen, BUB1B, "
               "sin usar panel de genes ni ClinVar.", b)


# --------------------------------------------------------------------------
# Figura 2: el hueco de fase
# --------------------------------------------------------------------------
def spans_fase(chrom="15"):
    cache = WORK / "phase_spans.json"
    if cache.exists():
        return json.loads(cache.read_text(encoding="utf-8"))
    vcfs = sorted((ROOT / "data").glob("*.vcf.gz"))
    if not vcfs:
        return []
    blocks = defaultdict(list)
    visto = False
    with gzip.open(vcfs[0], "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 10)
            if f[0] != chrom:
                if visto:
                    break
                continue
            visto = True
            if len(f) < 10:
                continue
            kv = dict(zip(f[8].split(":"), f[9].rstrip("\n").split(":")))
            pid = kv.get("PID", "")
            if pid and pid != ".":
                blocks[pid].append(int(f[1]))
    spans = sorted(max(v) - min(v) for v in blocks.values() if len(v) >= 2)
    WORK.mkdir(exist_ok=True)
    cache.write_text(json.dumps(spans), encoding="utf-8")
    return spans


def fig_fase():
    spans = spans_fase()
    distancia = 10911
    W, H = 760, 380
    x0, x1, y0, y1 = 70, 730, 96, 300
    lo, hi = 0.0, 5.0                      # log10: de 1 pb a 100 kb

    def X(v):
        return x0 + (math.log10(max(v, 1)) - lo) / (hi - lo) * (x1 - x0)

    nb = 60
    counts = [0] * nb
    for s in spans:
        k = int((math.log10(max(s, 1)) - lo) / (hi - lo) * nb)
        counts[min(nb - 1, max(0, k))] += 1
    cmax = max(counts) if counts else 1
    maximo = max(spans) if spans else 206

    b = [text(28, 40, "Por qué la fase no se puede resolver con estos datos", "t ink"),
         text(28, 62, "Tamaño de los %s bloques de fase del cromosoma 15, en escala logarítmica."
              % fmt(len(spans)), "st ink2")]
    for v, lab in ((1, "1 pb"), (10, "10 pb"), (100, "100 pb"), (1000, "1 kb"),
                   (10000, "10 kb"), (100000, "100 kb")):
        b.append('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" class="grid"/>' % (X(v), y0, X(v), y1))
        b.append(text(X(v), y1 + 18, lab, "sm ink2", "middle"))
    b.append('<line x1="%d" y1="%d" x2="%d" y2="%d" class="axis"/>' % (x0, y1, x1, y1))

    bw = (x1 - x0) / nb
    for i, c in enumerate(counts):
        if not c:
            continue
        h = (c / cmax) * (y1 - y0 - 30)
        b.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" rx="2" class="c1"/>'
                 % (x0 + i * bw + 1, y1 - h, bw - 2, h))

    # zona fuera del alcance de cualquier lectura
    b.append('<rect x="%.1f" y="%d" width="%.1f" height="%d" class="band" opacity="0.7"/>'
             % (X(maximo), y0, X(distancia) - X(maximo), y1 - y0))
    b.append('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" class="mark2"/>' % (X(maximo), y0, X(maximo), y1))
    # A la izquierda de su linea: a la derecha chocaba con la etiqueta de la distancia.
    b.append(text(X(maximo) - 6, y0 + 16, "bloque más largo: %s pb" % fmt(maximo), "lb ink2", "end"))

    b.append('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" class="mark"/>' % (X(distancia), y0, X(distancia), y1))
    b.append(text(X(distancia) - 8, y0 + 16, "distancia entre las dos variantes", "lb ink", "end"))
    b.append(text(X(distancia) - 8, y0 + 34, "de BUB1B: %s pb" % fmt(distancia), "lb ink", "end",
                  extra=' font-weight="600"'))
    b.append(text(X(distancia) - 8, y0 + 52, "%d veces el bloque más largo" % round(distancia / maximo),
                  "lb ink2", "end"))

    b.append(text(28, H - 22, "Ninguna lectura de este experimento abarca ambas posiciones: "
                  "la configuración en trans es una inferencia, no una observación.", "sm ink2"))
    return svg(W, H, "Tamano de los bloques de fase frente a la distancia entre las variantes",
               "Histograma logaritmico de los bloques de fase del cromosoma 15. El mas largo mide "
               "%d pb; las dos variantes de BUB1B estan a %d pb." % (maximo, distancia), b)


# --------------------------------------------------------------------------
# Figura 3: barrido de BAF
# --------------------------------------------------------------------------
def fig_baf():
    data = json.loads((WORK / "baf_scan.json").read_text(encoding="utf-8"))
    auto = [w for w in data["windows"] if w["chrom"].replace("chr", "").isdigit()]
    for w in auto:
        w["_c"] = int(w["chrom"].replace("chr", ""))
    auto.sort(key=lambda w: (w["_c"], w["start"]))
    umbral = 1.25

    W, H = 900, 414
    x0, x1, y0, y1 = 64, 870, 124, 344
    ymax = max(5.0, math.ceil(max(w.get("baf_excess_norm", 0) for w in auto)))
    gap = 6
    ncrom = 22
    n = len(auto)
    unit = (x1 - x0 - gap * (ncrom - 1)) / max(1, n)

    def Y(v):
        return y1 - (min(v, ymax) / ymax) * (y1 - y0)

    b = [text(28, 40, "Aneuploidía en mosaico: ningún evento detectable en sangre", "t ink"),
         text(28, 62, "Dispersión del balance alélico frente a lo esperado por binomial, "
              "ventanas de 10 Mb, normalizada a la mediana del genoma.", "st ink2")]
    for v in range(0, int(ymax) + 1):
        b.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" class="grid"/>' % (x0, Y(v), x1, Y(v)))
        b.append(text(x0 - 8, Y(v) + 4, v, "sm ink2 num", "end"))
    b.append(text(x0, y0 - 12, "exceso de dispersión", "sm ink2"))

    x = x0
    pos = []
    for c in range(1, ncrom + 1):
        ws = [w for w in auto if w["_c"] == c]
        ancho = unit * len(ws)
        if c % 2 == 0 and ws:
            b.append('<rect x="%.1f" y="%d" width="%.1f" height="%d" class="band"/>' % (x - gap / 2, y0, ancho + gap, y1 - y0))
        if ws:
            b.append(text(x + ancho / 2, y1 + 16, c, "sm ink2 num", "middle"))
        for k, w in enumerate(ws):
            pos.append((x + unit * (k + 0.5), w))
        x += ancho + gap

    b.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" class="mark2"/>' % (x0, Y(umbral), x1, Y(umbral)))
    b.append(text(x1, Y(umbral) - 6, "umbral 1,25", "sm ink2", "end"))

    sobre = 0
    for px, w in pos:
        v = w.get("baf_excess_norm", 0)
        cls = "c2" if v >= umbral else "c1"
        sobre += v >= umbral
        b.append('<circle cx="%.1f" cy="%.1f" r="4" class="%s ring"/>' % (px, Y(v), cls))

    ly = 86
    b.append('<circle cx="34" cy="%d" r="5" class="c1"/>' % ly)
    b.append(text(46, ly + 4, "ventana dentro de lo esperado", "lb ink"))
    b.append('<circle cx="262" cy="%d" r="5" class="c2"/>' % ly)
    b.append(text(274, ly + 4, "pico aislado sobre el umbral (%d)" % sobre, "lb ink"))

    b.append(text(28, H - 24, "Los picos caen en centrómeros y brazos acrocéntricos, donde el mapeo falla. "
                  "Ninguno forma un segmento de 3 ventanas seguidas.", "sm ink2"))
    b.append(text(28, H - 8, "Límite de detección: fracción celular de ~10-15 % a 44x. La aneuploidía "
                  "de MVA se documenta por cariotipo en células cultivadas.", "sm ink2"))
    return svg(W, H, "Barrido de aneuploidia en mosaico por ventanas de 10 Mb",
               "Exceso de dispersion del balance alelico en %d ventanas autosomicas. %d superan el "
               "umbral de 1,25 como picos aislados; ningun segmento continuo." % (n, sobre), b)


# --------------------------------------------------------------------------
# Figura 4: paciente frente a controles
# --------------------------------------------------------------------------
def metricas(work: Path):
    sweep = work / "gene_sweep.annot.tsv"
    cand = work / "candidates.tsv"
    cv = work / "candidates.clinvar.tsv"
    if not (sweep.exists() and cand.exists() and cv.exists()):
        return None
    c = list(csv.DictReader(cand.open(encoding="utf-8"), delimiter="\t"))
    cvr = list(csv.DictReader(cv.open(encoding="utf-8"), delimiter="\t"))
    genes = defaultdict(list)
    for r in c:
        genes[r["gene"]].append(r)
    cod = {g: [x for x in v if x["region"] == "codificante" or "splicing" in x["region"]]
           for g, v in genes.items()}
    return {
        "raras": len(c),
        "genes2": sum(1 for v in genes.values() if len(v) >= 2),
        "cod1": sum(1 for v in cod.values() if v),
        "cod2": sum(1 for v in cod.values() if len(v) >= 2),
        "clinvar": sum(1 for r in cvr if r.get("clinvar_sig", "").startswith(("Pathogenic", "Likely_pathogenic"))),
    }


def fig_controles():
    muestras = [("Paciente", WORK, True)]
    for tag in ("hg002", "hg001", "hg005"):
        d = ROOT / ("work_" + tag)
        if (d / "candidates.clinvar.tsv").exists():
            muestras.append((tag.upper(), d, False))
    datos = [(n, metricas(d), p) for n, d, p in muestras]
    datos = [x for x in datos if x[1]]

    # Los dos primeros paneles NO separan al paciente de los controles, y se
    # muestran por eso mismo: con un solo control parecia que "una rara
    # codificante" discriminaba, y HG005 lo desmintio. Los dos ultimos si separan.
    paneles = [("raras", "Variantes raras", "en el panel · no separa"),
               ("cod1", "Genes con una rara", "codificante · no separa"),
               ("cod2", "Genes con dos raras", "codificantes · separa"),
               ("clinvar", "Patogénicas exactas", "en ClinVar · separa")]
    W = 900
    H = 420
    pw = (W - 56) / len(paneles)
    top, base = 150, 330

    b = [text(28, 40, "El mismo pipeline sobre el paciente y sobre genomas sanos", "t ink"),
         text(28, 62, "Ni las variantes raras ni una rara codificante separan al paciente: "
              "lo hacen dos en un mismo gen y la patogenicidad conocida.", "st ink2")]
    b.append('<rect x="28" y="84" width="12" height="12" rx="3" class="c1"/>')
    b.append(text(46, 94, "Paciente", "lb ink"))
    b.append('<rect x="130" y="84" width="12" height="12" rx="3" class="ctx"/>')
    b.append(text(148, 94, "Genoma sano de referencia (GIAB), rotulado por nombre", "lb ink"))

    for i, (key, t1, t2) in enumerate(paneles):
        px = 28 + i * pw
        b.append(text(px + 8, top - 30, t1, "lb ink", extra=' font-weight="600"'))
        b.append(text(px + 8, top - 14, t2, "lb ink2"))
        vmax = max(max(d[1][key] for d in datos), 1)
        n = len(datos)
        # Separacion ancha a proposito: cada barra lleva su nombre debajo, y con
        # barras juntas los nombres se fundian en uno solo.
        paso = min(52, (pw - 24) / n)
        bw = min(30, paso - 18)
        x = px + 12
        b.append('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" class="axis"/>' % (px + 8, base, px + pw - 16, base))
        for nombre, m, es_paciente in datos:
            v = m[key]
            h = (v / vmax) * (base - top - 20)
            cls = "c1" if es_paciente else "ctx"
            if v > 0:
                b.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" rx="4" class="%s"/>'
                         % (x, base - h, bw, h, cls))
            b.append(text(x + bw / 2, base - h - 6, v, "lb ink num", "middle",
                          extra=' font-weight="600"' if es_paciente else ""))
            b.append(text(x + bw / 2, base + 16, nombre, "sm ink2", "middle"))
            x += paso

    ncontrol = len(datos) - 1
    b.append(text(28, H - 26, "Candidato final: BUB1B en el paciente, ninguno en %s. "
                  "Mismo código, mismos umbrales, ningún parámetro ajustado."
                  % ("el control" if ncontrol == 1 else "los %d controles" % ncontrol), "sm ink2"))
    b.append(text(28, H - 10, "Cada panel tiene su propia escala: compara alturas solo dentro de un panel.",
                  "sm ink2"))
    return svg(W, H, "Paciente frente a genomas sanos de referencia",
               "Cuatro metricas del pipeline para el paciente y %d genomas de control." % ncontrol, b)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    hechas = []
    for nombre, fn in (("fig1_embudo_ciego.svg", fig_embudo),
                       ("fig2_fase.svg", fig_fase),
                       ("fig3_baf_genoma.svg", fig_baf),
                       ("fig4_controles.svg", fig_controles)):
        try:
            (OUT / nombre).write_text(fn(), encoding="utf-8")
            hechas.append(nombre)
        except FileNotFoundError as e:
            print("  se omite %s: falta %s" % (nombre, e.filename))
    for h in hechas:
        print("  reports/figures/%s" % h)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
