# Clasificación ACMG/AMP — par candidato en BUB1B

Criterios de Richards et al. 2015 (PMID 25741868), con la revisión de PVS1 de
Abou Tayoun et al. 2018 (PMID 30192042). Transcrito de referencia NM_001211.6
(MANE Select), GRCh38. Proteína O60566, 1050 aminoácidos.

Se documenta cada criterio aplicado **y cada criterio que se decidió no aplicar**,
con el motivo. Lo segundo es lo que suele faltar y es donde se cuela el sesgo:
aplicar un criterio dudoso a favor sube la clasificación sin que nadie lo note.

---

## Variante 1 — p.Leu737Ter (nonsense)

| Dato | Valor |
|---|---|
| Consecuencia | codón de parada prematuro, exón 17 de 23 |
| Distancia a la última unión exón-exón | 748 nt corriente arriba |
| Predicción de NMD | sujeto a degradación (regla de los 50 nt) |
| gnomAD v4.1 | AF 7.87e-05, AN 1.461.846, **0 homocigotos** |
| ClinVar | Pathogenic/Likely_pathogenic, múltiples submitters, sin conflictos, *Mosaic variegated aneuploidy syndrome 1* |
| Calidad | DP 46, GQ 99, PASS |

### Criterios aplicados

**PVS1 (Very Strong).** Variante nula en un gen donde la pérdida de función
bialélica es el mecanismo establecido de enfermedad: Hanks et al. 2004 describió
MVA1 por mutaciones bialélicas en BUB1B (PMID 15475955). El árbol de decisión de
PVS1 exige además que el transcrito truncado no escape a NMD, y aquí no escapa:
el codón de parada queda 748 nt antes de la última unión exón-exón, muy por
encima del umbral de 50 nt. Se elimina también la salida del árbol por "exón
prescindible": el truncamiento en el residuo 737 de 1050 elimina el dominio
766-1050 completo.

**PM2_Supporting.** AF 7.87e-05 sin homocigotos en 1.46 millones de alelos. No es
ausente, pero para una enfermedad recesiva ultra-rara una frecuencia de portador
de ese orden es compatible. Se aplica en nivel *supporting*, siguiendo la
recomendación del SVI de no usar PM2 en nivel moderado como criterio de
frecuencia.

### Criterios que NO se aplican, y por qué

- **PM3** (detectada en trans con una variante patogénica, para recesivos): *no se
  aplica*. La fase no está resuelta. Este es el criterio que más pesaría y es el
  que los datos no permiten.
- **PS3** (estudios funcionales): no hay ensayo funcional de esta variante.
- **PP5** (fuentes reputadas la reportan patogénica): en desuso por recomendación
  del SVI. La entrada de ClinVar se cita como contexto, no como criterio.
- **PS4** (prevalencia en afectados): requiere conteo de probandos en literatura,
  que no se ha hecho.

### Clasificación

**PVS1 + PM2_Supporting → Patogénica / Probablemente patogénica.** Consistente con
la clasificación vigente en ClinVar por múltiples submitters sin conflictos.

---

## Variante 2 — p.Asn1002Lys (missense)

| Dato | Valor |
|---|---|
| Consecuencia | missense, exón 23 de 23 (último exón) |
| Dominio | dentro del dominio 766-1050 de O60566 |
| gnomAD v4.1 | AF 6.84e-07 (1 alelo en 1.46 millones), **0 homocigotos** |
| ClinVar | sin entrada |
| Calidad | DP 28, GQ 99, PASS |

### Criterios aplicados

**PM2_Supporting.** Un solo alelo en 1.46 millones, sin homocigotos. Prácticamente
privada.

### Criterios que NO se aplican, y por qué

- **PM1** (dominio funcional crítico bien establecido, sin variación benigna): *no
  se aplica*, y conviene explicar por qué, porque la tentación es fuerte. UniProt
  anota el tramo 766-1050 como "Protein kinase", y la variante cae dentro. Pero
  ese dominio es un **pseudoquinasa**: no tiene actividad catalítica, y su papel
  es de andamiaje, reclutando PP2A-B56 al cinetocoro (PMID 33207204, PMID
  35525552). Llamarlo "dominio quinasa" y deducir pérdida de catálisis sería un
  error de razonamiento. Además, PM1 exige ausencia de variación benigna en el
  dominio, y gnomAD observa 1120 missense en el gen contra 1295 esperadas: no hay
  depleción que sostenga el criterio.
- **PP2** (missense como mecanismo en un gen con poca variación missense benigna):
  *no se aplica*. La restricción de missense de gnomAD para BUB1B es z = 1.78, por
  debajo de cualquier umbral razonable.
- **PP3** (evidencia computacional convergente): *no se aplica*. No se pudo obtener
  REVEL, AlphaMissense ni CADD sin enviar la posición del paciente a un servicio
  externo, y el índice remoto de AlphaMissense no está disponible en el bucket
  público. Queda declarado como pendiente, no como evidencia ausente.
- **PM3**: *no se aplica*, misma razón que en la variante 1.

### Clasificación

**PM2_Supporting como único criterio → Variante de significado incierto (VUS).**

Esto es incómodo y es correcto. El segundo alelo de un compuesto heterocigoto
suele ser el que no está clasificado; decir que es patogénica porque encaja en la
hipótesis sería razonamiento circular.

---

## Interpretación conjunta

El caso clínico es MVA. El paciente porta un alelo nulo de BUB1B clasificado como
patogénico para MVA1 y una segunda variante missense prácticamente privada en el
mismo gen. Bajo herencia autosómica recesiva, esa es la arquitectura esperada, y
coincide con lo descrito desde el trabajo original de Hanks (PMID 15475955): los
pacientes rara vez son nulos homocigotos.

**Lo que sostiene la hipótesis:** el gen correcto para el fenotipo, un alelo
patogénico conocido con NMD predicho, un segundo alelo raro en el mismo gen,
ausencia de homocigotos en gnomAD para ambos, y calidad alta en las dos llamadas.

**Lo que no está demostrado:** que las dos variantes estén en trans. Los bloques
de fase de este VCF llegan a 206 pb como máximo y las variantes están a 10.911 pb.
Si estuvieran en cis, el paciente tendría un alelo intacto y la hipótesis se cae.

**Qué lo resolvería, en orden de coste:** genotipado dirigido de los padres; PCR
de alelo específico o secuenciación de lectura larga sobre la región; ensayo
funcional del missense en el contexto del reclutamiento de PP2A-B56.

**Cómo cambiaría la clasificación.** Confirmar la fase en trans habilita PM3 sobre
la variante 2. Con PM3 en nivel moderado, PM2_Supporting + PM3 sigue siendo VUS;
hace falta además PP3 o un ensayo funcional para llegar a probablemente
patogénica. Es decir: incluso resolviendo la fase, el segundo alelo necesita
evidencia funcional. Conviene decirlo antes de que lo diga el panel.

---

## Referencias

Todas verificadas contra PubMed en 2026-09-12.

- PMID 25741868 — Richards et al., *Genet Med* 2015. Estándares ACMG/AMP.
- PMID 30192042 — Abou Tayoun et al., *Hum Mutat* 2018. Interpretación de PVS1.
- PMID 15475955 — Hanks et al., *Nat Genet* 2004. Aneuploidía constitucional y
  predisposición tumoral por mutaciones bialélicas en BUB1B.
- PMID 33207204 — *Cell Rep* 2020. El dominio pseudoquinasa de BUBR1 promueve el
  reclutamiento de PP2A-B56 al cinetocoro.
- PMID 35525552 — *Methods Enzymol* 2022. Consideraciones sobre la fosforilación
  del pseudoquinasa BUBR1.
