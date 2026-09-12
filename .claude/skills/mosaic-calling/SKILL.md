---
name: mosaic-calling
description: Protocolo para llamado de variantes en mosaico y perfilado de aneuploidia variegada. Usar en cualquier tarea que involucre VAF, MuTect2 tumor-only, MosaicHunter, DeepMosaic, B-allele frequency o cobertura por cromosoma.
---

# Llamado en mosaico

MVA es mosaico por definicion. Un pipeline germinal estandar pierde la senal causal.

## Regla central
VAF intermedios (~5-35%) son **senal**, no ruido. Nunca los filtres por defecto.

## Callers (corre los tres, reporta union e interseccion)
- **MuTect2 tumor-only** con panel of normals publico (evita artefactos recurrentes)
- **MosaicHunter** modelo single-sample
- **DeepMosaic** sobre el BAM final

## Filtros de artefacto obligatorios
- Sesgo de hebra (strand bias): descarta si todas las lecturas alternas van en una hebra
- Posicion en la lectura: descarta si el alelo alterno se concentra en los extremos
- Homopolimeros y repeticiones de baja complejidad
- Regiones de mapeo ambiguo (mappability < 1)
- Cobertura minima 30x en el sitio

## Validacion
Inspeccion visual en IGV de las 20 mejores. Screenshot guardado en pipelines/04_mosaic/igv/.
Una variante en mosaico sin inspeccion visual no entra a la lista final.

## Aneuploidia variegada
Perfila por cromosoma: cobertura normalizada + B-allele frequency. Un cromosoma con
ganancia/perdida en mosaico muestra desviacion de cobertura Y separacion de BAF del 0.5.
Estima la fraccion celular afectada. Esto es **evidencia ortogonal fuerte** del fenotipo
y casi nadie lo va a reportar: uselo en el writeup.
