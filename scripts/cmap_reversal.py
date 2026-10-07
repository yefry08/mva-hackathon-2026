#!/usr/bin/env python3
"""Firma transcripcional de la perdida de BUB1B y farmacos que la revierten.

Todo a nivel de gen y de compuesto, con fuentes publicas y gratuitas (LINCS
L1000 via Enrichr y L1000CDS2, Ma'ayan Lab; ChEMBL para fase clinica). Ningun
dato del paciente entra ni sale: la consulta es "que pasa en una celula sin
BUB1B", no "que tiene este nino".

Tres preguntas, cada una con su control:

1. Que hace la perdida de BUB1B al transcriptoma. Firma consenso de knockout
   CRISPR de LINCS L1000, enriquecida contra MSigDB Hallmark. Control: la misma
   prueba sobre las ~5.200 firmas de knockout de la libreria, para saber si lo
   que sale es propio de BUB1B o le pasa a cualquier gen esencial que se apaga.

2. Que compuestos revierten esa firma (L1000CDS2, modo reversion). Control: la
   misma consulta con firmas de knockout tomadas al azar. Un compuesto que
   revierte tambien la perdida de genes que nada tienen que ver con el
   checkpoint revierte "celula enferma", no "celula sin BUB1B".

3. Si los que quedan estan aprobados (ChEMBL, fase maxima).

Salidas (datos publicos):
  reports/cmap_signature.tsv   enriquecimiento Hallmark y su percentil nulo
  reports/cmap_reversal.tsv    compuestos, especificidad y fase clinica

Uso: python scripts/cmap_reversal.py
"""
from __future__ import annotations

import json
import math
import random
import sys
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "ref" / "cmap_cache"
OUT_SIG = ROOT / "reports" / "cmap_signature.tsv"
OUT_REV = ROOT / "reports" / "cmap_reversal.tsv"

ENRICHR = "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName="
CDS2 = "https://maayanlab.cloud/L1000CDS2/query"
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data/molecule.json?"

DIANA = "BUB1B"
# Genes del mismo modulo (cinetocoro y checkpoint) presentes en la libreria.
VIA = ["BUB1", "TTK", "CENPE", "AURKB", "PLK1", "ZWINT", "NDC80"]
N_NULO = 40
SEMILLA = 2026
TOP = 50


def _get(url: str, intentos: int = 4, timeout: int = 180) -> bytes:
    for i in range(intentos):
        try:
            return urllib.request.urlopen(url, timeout=timeout).read()
        except Exception as e:  # red: reintentar, nunca devolver vacio
            if i == intentos - 1:
                raise RuntimeError(f"fallo definitivo: {url}: {e}")
            time.sleep(3 * (i + 1))


def libreria(nombre: str) -> dict[str, list[str]]:
    CACHE.mkdir(parents=True, exist_ok=True)
    f = CACHE / f"{nombre}.txt"
    if not f.exists():
        f.write_bytes(_get(ENRICHR + nombre))
    sets = {}
    for linea in f.read_text().splitlines():
        p = linea.split("\t")
        genes = [g.split(",")[0].strip().upper() for g in p[2:] if g.strip()]
        if genes:
            sets[p[0]] = genes
    if not sets:
        raise RuntimeError(f"libreria vacia: {nombre}")
    return sets


# --- 1. enriquecimiento hipergeometrico, local ------------------------------

def _lchoose(n: int, k: int) -> float:
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hipergeom_sf(k: int, N: int, K: int, n: int) -> float:
    """P(X >= k) con X ~ Hipergeometrica(N poblacion, K exitos, n extraccion)."""
    if k <= 0:
        return 1.0
    tope = min(K, n)
    if k > tope:
        return 0.0
    # Primer termino en log; los siguientes por recurrencia, que decaen rapido.
    t = math.exp(_lchoose(K, k) + _lchoose(N - K, n - k) - _lchoose(N, n))
    total = t
    for i in range(k, tope):
        t *= (K - i) * (n - i) / ((i + 1) * (N - K - n + i + 1))
        total += t
        if t < total * 1e-12:
            break
    return min(1.0, total)


def enriquecer(genes: set[str], hallmarks: dict[str, set[str]], fondo: set[str]):
    g = genes & fondo
    out = {}
    for h, hs in hallmarks.items():
        k = len(g & hs)
        out[h] = (k, hipergeom_sf(k, len(fondo), len(hs), len(g)))
    return out


# --- 2. reversion en L1000CDS2 ----------------------------------------------

def revertir(nombre: str, up: list[str], dn: list[str]) -> list[dict]:
    CACHE.mkdir(parents=True, exist_ok=True)
    f = CACHE / f"cds2_{nombre.replace(' ', '_')}.json"
    if f.exists():
        return json.loads(f.read_text())
    body = {"data": {"upGenes": up, "dnGenes": dn},
            "config": {"aggravate": False, "searchMethod": "geneSet",
                       "share": False, "combination": False,
                       "db-version": "latest"},
            "metadata": []}
    for i in range(4):
        try:
            req = urllib.request.Request(
                CDS2, data=json.dumps(body).encode(),
                headers={"Content-Type": "application/json"})
            meta = json.loads(urllib.request.urlopen(req, timeout=240).read())["topMeta"]
            break
        except Exception as e:
            if i == 3:
                raise RuntimeError(f"L1000CDS2 fallo para {nombre}: {e}")
            time.sleep(5 * (i + 1))
    if not meta:
        raise RuntimeError(f"L1000CDS2 devolvio vacio para {nombre}")
    f.write_text(json.dumps(meta))
    return meta


def compuesto(m: dict) -> str:
    d = str(m.get("pert_desc", "")).strip()
    return (d if d and d != "-666" else str(m.get("pert_id", "?"))).lower()


# --- 3. fase clinica en ChEMBL ----------------------------------------------

def fase_chembl(nombre: str) -> str:
    if nombre.startswith("brd-"):
        return "sin_nombre"
    for campo in ("pref_name__iexact", "molecule_synonyms__molecule_synonym__iexact"):
        q = urllib.parse.urlencode({campo: nombre, "format": "json", "limit": 5})
        d = json.loads(_get(CHEMBL + q))
        fases = [float(m["max_phase"]) for m in d.get("molecules", [])
                 if m.get("max_phase") not in (None, "")]
        if fases:
            return str(int(max(fases)))
        if d.get("molecules"):
            return "0"
    return "no_en_chembl"


def main() -> int:
    ko = libreria("LINCS_L1000_CRISPR_KO_Consensus_Sigs")
    hm = {h: set(g) for h, g in libreria("MSigDB_Hallmark_2020").items()}
    firmas = defaultdict(dict)
    for k, genes in ko.items():
        gen, sentido = k.rsplit(" ", 1)
        firmas[gen][sentido.lower()] = genes
    firmas = {g: s for g, s in firmas.items() if "up" in s and "down" in s}
    if DIANA not in firmas:
        raise RuntimeError("BUB1B no esta en la libreria de knockouts")
    fondo = set().union(*(set(s["up"]) | set(s["down"]) for s in firmas.values()))
    hm = {h: s & fondo for h, s in hm.items() if len(s & fondo) >= 15}
    print(f"firmas de knockout: {len(firmas)}  fondo: {len(fondo)} genes  "
          f"hallmarks: {len(hm)}")

    # 1. Firma de BUB1B frente a todas las demas
    filas_sig = []
    for sentido in ("up", "down"):
        diana = enriquecer(set(firmas[DIANA][sentido]), hm, fondo)
        nulo = defaultdict(list)
        for gen, s in firmas.items():
            if gen == DIANA:
                continue
            for h, (_, p) in enriquecer(set(s[sentido]), hm, fondo).items():
                nulo[h].append(p)
        for h, (k, p) in diana.items():
            pct = sum(1 for q in nulo[h] if q <= p) / len(nulo[h])
            via = {g: enriquecer(set(firmas[g][sentido]), {h: hm[h]}, fondo)[h][1]
                   for g in VIA if g in firmas}
            filas_sig.append((sentido, h, k, p, pct, via))
    filas_sig.sort(key=lambda r: r[3])
    OUT_SIG.parent.mkdir(exist_ok=True)
    with OUT_SIG.open("w", encoding="utf-8") as f:
        f.write("sentido\thallmark\tgenes_solapados\tp_hipergeom\t"
                "fraccion_knockouts_igual_o_mas_enriquecidos\t"
                + "\t".join(f"p_{g}" for g in VIA) + "\n")
        for s, h, k, p, pct, via in filas_sig:
            f.write(f"{s}\t{h}\t{k}\t{p:.3g}\t{pct:.4f}\t"
                    + "\t".join(f"{via.get(g, float('nan')):.3g}" for g in VIA) + "\n")
    print("\nFirma de BUB1B, hallmarks con p < 1e-3:")
    for s, h, k, p, pct, via in filas_sig:
        if p < 1e-3:
            nvia = sum(1 for q in via.values() if q < 1e-3)
            print(f"  {s:4s} {h:40s} k={k:3d} p={p:.1e}  "
                  f"knockouts iguales o mas: {pct:.1%}  via con p<1e-3: {nvia}/{len(via)}")

    # 2. Reversion, con nulo de knockouts al azar
    rng = random.Random(SEMILLA)
    candidatos_nulo = sorted(g for g in firmas if g != DIANA and g not in VIA)
    nulos = rng.sample(candidatos_nulo, N_NULO)
    hits = {}
    for gen in [DIANA] + [g for g in VIA if g in firmas] + nulos:
        meta = revertir(gen, firmas[gen]["up"], firmas[gen]["down"])[:TOP]
        hits[gen] = meta
        print(f"  reversion {gen:10s} {len(meta)} firmas", file=sys.stderr)

    def presencia(genes):
        c = Counter()
        for g in genes:
            c.update({compuesto(m) for m in hits[g]})
        return c

    c_nulo = presencia(nulos)
    c_via = presencia([g for g in VIA if g in hits])
    n_via = sum(1 for g in VIA if g in hits)
    diana = defaultdict(lambda: {"n": 0, "mejor": 0.0, "celulas": set()})
    for m in hits[DIANA]:
        d = diana[compuesto(m)]
        d["n"] += 1
        d["mejor"] = max(d["mejor"], float(m.get("score", 0)))
        d["celulas"].add(str(m.get("cell_id")))

    # Replica independiente: knockdown por shRNA (otra tecnica, otro experimento).
    # Solo cuentan las lineas donde el propio BUB1B aparece entre los genes que
    # bajan: en las demas el knockdown no es verificable y la firma no es "perdida
    # de BUB1B". Consenso: genes que cambian en todas las lineas verificadas.
    lib_up = libreria("L1000_Kinase_and_GPCR_Perturbations_up")
    lib_dn = libreria("L1000_Kinase_and_GPCR_Perturbations_down")
    todas = sorted(k for k in lib_dn if k.startswith(f"{DIANA} knockdown"))
    verificadas = [k for k in todas if DIANA in lib_dn[k]]
    if len(verificadas) < 2:
        raise RuntimeError(f"knockdown verificado en {len(verificadas)} lineas")
    rep = {}
    for sentido, lib in (("up", lib_up), ("down", lib_dn)):
        c = Counter(x for k in verificadas for x in set(lib[k]))
        rep[sentido] = sorted(x for x, n in c.items() if n == len(verificadas))
    print(f"\nreplica shRNA: knockdown verificado en {len(verificadas)} de "
          f"{len(todas)} lineas ({', '.join(k.rsplit(' ', 1)[1] for k in verificadas)}); "
          f"consenso up={len(rep['up'])} down={len(rep['down'])}")
    for h in ("Interferon Alpha Response", "Interferon Gamma Response",
              "p53 Pathway", "Apoptosis"):
        k, p = enriquecer(set(rep["up"]), {h: hm[h]}, fondo)[h]
        print(f"  shRNA up {h:30s} k={k:3d} p={p:.1e}")
    for h in ("G2-M Checkpoint", "E2F Targets"):
        k, p = enriquecer(set(rep["down"]), {h: hm[h]}, fondo)[h]
        print(f"  shRNA down {h:28s} k={k:3d} p={p:.1e}")
    meta_rep = revertir(f"{DIANA}_shRNA_verificado", rep["up"], rep["down"])[:TOP]
    c_rep = Counter(compuesto(m) for m in meta_rep)

    # Coincidencia CRISPR-shRNA contra el azar: cuantas firmas de farmaco
    # comparte la lista shRNA con la de cualquier otro knockout de la corrida.
    sig = lambda meta: {m.get("sig_id") for m in meta}
    compartidas = sig(meta_rep) & sig(hits[DIANA])
    otros = [g for g in hits if g != DIANA]
    nulo_comp = [len(sig(meta_rep) & sig(hits[g])) for g in otros]
    p_comp = sum(1 for n in nulo_comp if n >= len(compartidas)) / len(otros)
    detalle = []
    for m in meta_rep:
        if m.get("sig_id") in compartidas:
            ov = m.get("overlap", {})
            detalle.append({"compuesto": compuesto(m), "sig_id": m.get("sig_id"),
                            "placa": str(m.get("sig_id")).split(":")[0],
                            "genes_que_sostienen_el_match_shrna":
                                sorted(ov.get("up/dn", []) + ov.get("dn/up", []))})
    print(f"  firmas compartidas CRISPR-shRNA: {len(compartidas)}; knockouts al "
          f"azar con igual o mas: {p_comp:.0%}")
    for d in detalle:
        print(f"    {d['compuesto']:26s} placa {d['placa']}  match shRNA por "
              f"{len(d['genes_que_sostienen_el_match_shrna'])} genes: "
              f"{d['genes_que_sostienen_el_match_shrna']}")

    # Lo que habria propuesto un pipeline ingenuo: consenso de las 8 lineas sin
    # comprobar que el knockdown funciono (BUB1B solo baja en 2 de ellas).
    ing = {}
    for sentido, lib in (("up", lib_up), ("down", lib_dn)):
        c = Counter(x for k in todas for x in set(lib[k]))
        ing[sentido] = sorted(x for x, n in c.items() if n >= 3)
    meta_ing = revertir(f"{DIANA}_shRNA", ing["up"], ing["down"])[:TOP]
    c_ing = Counter(compuesto(m) for m in meta_ing)
    ingenuo = [{"compuesto": n, "firmas_top50": k, "fase_chembl": fase_chembl(n)}
               for n, k in c_ing.most_common(10)]
    print(f"  pipeline ingenuo (8 lineas sin verificar): "
          f"{[(d['compuesto'], d['fase_chembl']) for d in ingenuo]}")

    filas = []
    for nombre, d in diana.items():
        filas.append({
            "compuesto": nombre, "firmas_top50": d["n"], "mejor_score": d["mejor"],
            "lineas": ",".join(sorted(d["celulas"])),
            "frac_nulo": c_nulo[nombre] / N_NULO,
            "frac_via": c_via[nombre] / n_via if n_via else float("nan"),
            "firmas_replica_shrna": c_rep[nombre],
            "fase_chembl": fase_chembl(nombre),
        })
    filas.sort(key=lambda r: (r["frac_nulo"], -r["firmas_top50"]))
    with OUT_REV.open("w", encoding="utf-8") as f:
        cols = list(filas[0])
        f.write("\t".join(cols) + "\n")
        for r in filas:
            f.write("\t".join(f"{r[c]:.3f}" if isinstance(r[c], float) else str(r[c])
                              for c in cols) + "\n")

    print(f"\nReversion de BUB1B: {len(filas)} compuestos distintos en el top {TOP}")
    print(f"{'compuesto':28s} {'firmas':>6s} {'nulo':>6s} {'via':>6s} "
          f"{'shRNA':>6s} {'fase':>12s}")
    for r in filas:
        print(f"{r['compuesto'][:28]:28s} {r['firmas_top50']:6d} "
              f"{r['frac_nulo']:6.0%} {r['frac_via']:6.0%} "
              f"{r['firmas_replica_shrna']:6d} {r['fase_chembl']:>12s}")
    especificos = [r for r in filas if r["frac_nulo"] <= 0.05]
    aprobados = [r for r in especificos if r["fase_chembl"] == "4"]
    replicados = [r for r in especificos if r["firmas_replica_shrna"] > 0]
    print(f"\nespecificos (en <=5% de knockouts al azar): {len(especificos)}; "
          f"aprobados (fase 4): {[r['compuesto'] for r in aprobados]}; "
          f"replican con shRNA: {[r['compuesto'] for r in replicados]}")
    resumen = {
        "firmas_knockout": len(firmas), "fondo_genes": len(fondo),
        "knockouts_nulo_reversion": N_NULO, "semilla": SEMILLA, "top": TOP,
        "shrna_lineas_total": len(todas),
        "shrna_lineas_knockdown_verificado": [k.rsplit(" ", 1)[1] for k in verificadas],
        "compuestos_reversion_bub1b": len(filas),
        "especificos_nulo_le_5pct": [r["compuesto"] for r in especificos],
        "firmas_compartidas_crispr_shrna": detalle,
        "p_firmas_compartidas_vs_knockouts": p_comp,
        "pipeline_ingenuo_top10": ingenuo,
    }
    (ROOT / "reports" / "cmap_summary.json").write_text(
        json.dumps(resumen, indent=2, ensure_ascii=False), encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
