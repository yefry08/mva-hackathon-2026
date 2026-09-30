#!/usr/bin/env python3
"""Paisaje farmacologico y de esencialidad de los genes del checkpoint mitotico.

Todo a nivel de gen, con fuentes publicas y gratuitas, sin tocar datos del
paciente: lo unico que sale de la maquina son simbolos e identificadores Ensembl
de los 32 genes del panel.

Dos preguntas, y cada una alimenta una parte distinta del reporte:

1. Esencialidad (DepMap, via Open Targets). Si BUB1B es esencial para dividirse,
   la perdida bialelica completa no seria viable, y eso explica con datos por que
   el segundo alelo del paciente tiene que conservar algo de funcion: un missense
   hipomorfico, no un segundo nulo. DepMap viene de lineas tumorales, asi que la
   lectura es "la perdida no se tolera en celulas que proliferan", no mas.

2. Paisaje de farmacos (Open Targets y DGIdb). Que compuestos existen contra estos
   genes, en que fase, y en que direccion actuan. Un paciente hipomorfico tiene
   poca funcion de checkpoint; un inhibidor del checkpoint empuja hacia el mismo
   lado que la enfermedad. Si el paisaje es casi todo inhibidores, el filtro de
   realismo pediatrico del reporte deja de ser un argumento y pasa a ser un dato.

Salida: reports/drug_landscape.tsv (datos publicos a nivel de gen).

Uso: python scripts/drug_landscape.py
"""
from __future__ import annotations

import gzip
import json
import statistics
import sys
import time
import urllib.request
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GENES = ROOT / "ref" / "sac_genes.txt"
MANE = ROOT / "ref" / "MANE.summary.txt.gz"
OUT = ROOT / "reports" / "drug_landscape.tsv"

OT = "https://api.platform.opentargets.org/api/v4/graphql"
DGI = "https://dgidb.org/api/graphql"
ESENCIAL = -1.0          # umbral habitual de Chronos para "dependencia"


def post(url, query, intentos=4):
    data = json.dumps({"query": query}).encode()
    ultimo = None
    for n in range(intentos):
        try:
            req = urllib.request.Request(url, data=data,
                                         headers={"Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=90) as r:
                out = json.load(r)
            if "errors" in out:
                raise RuntimeError(out["errors"][0]["message"][:200])
            return out["data"]
        except Exception as e:                        # noqa: BLE001
            ultimo = e
            time.sleep(3 * 2 ** n)
    # Nunca un resultado vacio silencioso: si la fuente no responde, se sabe.
    raise RuntimeError("%s no respondio tras %d intentos: %s" % (url, intentos, ultimo))


def panel():
    genes = [l.split("#")[0].strip() for l in GENES.read_text(encoding="utf-8").splitlines()
             if l.split("#")[0].strip()]
    ens = {}
    with gzip.open(MANE, "rt", encoding="utf-8") as fh:
        cols = next(fh).rstrip("\n").split("\t")
        i_s, i_e = cols.index("symbol"), cols.index("Ensembl_Gene")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[i_s] in genes:
                ens[f[i_s]] = f[i_e].split(".")[0]
    return genes, ens


def open_targets(genes, ens):
    partes = []
    for i, g in enumerate(genes):
        if g in ens:
            partes.append('g%d: target(ensemblId:"%s") { approvedSymbol isEssential '
                          'tractability { modality label value } '
                          'depMapEssentiality { screens { geneEffect } } '
                          'drugAndClinicalCandidates { count rows { maxClinicalStage drug { name } } } }'
                          % (i, ens[g]))
    data = post(OT, "{ " + " ".join(partes) + " }")
    res = {}
    for v in data.values():
        if not v:
            continue
        effs = [s["geneEffect"] for t in (v.get("depMapEssentiality") or [])
                for s in t["screens"] if s["geneEffect"] is not None]
        dc = v.get("drugAndClinicalCandidates") or {}
        rows = dc.get("rows") or []
        res[v["approvedSymbol"]] = {
            "esencial": v.get("isEssential"),
            "n_lineas": len(effs),
            "efecto_mediana": statistics.median(effs) if effs else None,
            "frac_dependiente": (sum(e < ESENCIAL for e in effs) / len(effs)) if effs else None,
            "sm_tratable": any(x["value"] and x["modality"] == "SM" for x in v["tractability"] or []),
            "n_candidatos": dc.get("count") or 0,
            "fases": Counter(r["maxClinicalStage"] for r in rows),
            "farmacos": sorted({r["drug"]["name"] for r in rows if r.get("drug")}),
        }
    return res


def dgidb(genes):
    q = ('{ genes(names: [%s]) { nodes { name interactions { drug { name approved } '
         'interactionTypes { type directionality } } } } }'
         % ", ".join('"%s"' % g for g in genes))
    res = {}
    for n in post(DGI, q)["genes"]["nodes"]:
        tipos, direc, aprobados = Counter(), Counter(), set()
        for it in n["interactions"]:
            for t in it["interactionTypes"]:
                tipos[t["type"]] += 1
                direc[t.get("directionality") or "sin dato"] += 1
            if it["drug"]["approved"]:
                aprobados.add(it["drug"]["name"])
        res[n["name"]] = {"n_interacciones": len(n["interactions"]), "tipos": tipos,
                          "direccion": direc, "aprobados": sorted(aprobados)}
    return res


def main() -> int:
    genes, ens = panel()
    print("consultando Open Targets y DGIdb para %d genes (solo simbolos e IDs)..." % len(genes),
          flush=True)
    try:
        ot = open_targets(genes, ens)
        dg = dgidb(genes)
    except RuntimeError as e:
        print("ERROR: %s" % e, file=sys.stderr)
        return 2

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", encoding="utf-8") as out:
        out.write("gene\tesencial\tlineas_depmap\tefecto_mediana\tfrac_dependiente\t"
                  "sm_tratable\tcandidatos_clinicos\tfarmacos_clinicos\tinteracciones_dgidb\t"
                  "tipos_dgidb\taprobados_dgidb\n")
        for g in genes:
            o, d = ot.get(g, {}), dg.get(g, {})
            out.write("\t".join(str(x) for x in [
                g, o.get("esencial"), o.get("n_lineas"),
                "" if o.get("efecto_mediana") is None else "%.2f" % o["efecto_mediana"],
                "" if o.get("frac_dependiente") is None else "%.2f" % o["frac_dependiente"],
                o.get("sm_tratable"), o.get("n_candidatos", 0), ";".join(o.get("farmacos", [])),
                d.get("n_interacciones", 0),
                ";".join("%s:%d" % kv for kv in (d.get("tipos") or {}).items()),
                ";".join(d.get("aprobados", []))]) + "\n")

    print("\n%-9s %9s %8s %8s %11s %10s %9s" %
          ("gen", "esencial", "efecto", "depend.", "candidatos", "DGIdb", "aprobados"))
    for g in genes:
        o, d = ot.get(g, {}), dg.get(g, {})
        print("%-9s %9s %8s %7s%% %11d %10d %9d" % (
            g, "si" if o.get("esencial") else "no",
            "-" if o.get("efecto_mediana") is None else "%.2f" % o["efecto_mediana"],
            "-" if o.get("frac_dependiente") is None else "%.0f" % (100 * o["frac_dependiente"]),
            o.get("n_candidatos", 0), d.get("n_interacciones", 0), len(d.get("aprobados", []))))

    tipos = Counter()
    for d in dg.values():
        tipos.update(d["tipos"])
    aprob = sorted({a for d in dg.values() for a in d["aprobados"]})
    esenciales = [g for g in genes if ot.get(g, {}).get("esencial")]
    print("\n== lectura ==")
    b = ot.get("BUB1B", {})
    if b.get("efecto_mediana") is not None:
        print("BUB1B: esencial=%s, efecto mediano %.2f, dependiente en %.0f%% de %d lineas."
              % (b["esencial"], b["efecto_mediana"], 100 * b["frac_dependiente"], b["n_lineas"]))
    print("genes esenciales en el panel: %d de %d" % (len(esenciales), len(genes)))
    print("tipos de interaccion en DGIdb, todo el panel: %s" % dict(tipos.most_common()))
    print("farmacos aprobados que interactuan con algun gen del panel: %d" % len(aprob))
    if aprob:
        print("  %s" % ", ".join(aprob[:25]))
    print("\ndetalle en %s" % OUT.relative_to(ROOT))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
