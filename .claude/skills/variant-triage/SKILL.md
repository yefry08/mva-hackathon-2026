---
name: variant-triage
description: Reglas ACMG/AMP y criterios de filtrado especificos de este caso. Usar en cualquier tarea de clasificacion, filtrado o priorizacion de variantes.
---

# Triaje de variantes

## Modelo de herencia
MVA es **autosomica recesiva**. Prioriza:
1. Homocigotos con AF gnomAD < 0.001
2. Heterocigotos compuestos (dos variantes raras en el mismo gen)
3. Hemicigotos en regiones con deleccion

**Si encuentras un par heterocigoto compuesto, SIEMPRE envia ambas variantes.**
El scoring del hackathon da medio credito por identificar solo una de las dos.

## Umbrales
- AF gnomAD global < 0.001; popmax < 0.005
- CADD PHRED > 20 como senal, no como corte duro
- SpliceAI delta > 0.2 merece revision aunque sea sinonima o intronica
- AlphaMissense "likely pathogenic" como soporte, nunca como unico criterio

## No descartes
- Variantes non-coding en UTR y regiones regulatorias
- Variantes sinonimas cerca de sitios de splicing
- Genes sin asociacion previa a enfermedad (candidatos a gen nuevo)

## Doble ranking obligatorio
Produce SIEMPRE dos rankings y guardalos por separado:
- `ranking_ciego.tsv` — sin prior de genes MVA
- `ranking_prior.tsv` — con prior sobre BUB1B/CEP57/TRIP13 y red SAC

Documenta el delta. Un pipeline que solo encuentra la respuesta cuando se la das
no demuestra nada, y el panel lo va a notar.

## ACMG/AMP
Clasifica las top 50 documentando cada criterio aplicado (PVS1, PS1-4, PM1-6, PP1-5,
BA1, BS1-4, BP1-7). Sin la justificacion criterio por criterio, la clasificacion no cuenta.
