# Mecanismo: de la variante al fenotipo

Cadena causal para el par candidato en BUB1B. Cada flecha lleva una cita
verificada contra PubMed el 2026-09-12. Donde la evidencia es de modelo tumoral y
no de enfermedad constitucional, se dice.

---

## La cadena, flecha por flecha

### 1. BUBR1 es un componente central del checkpoint de ensamblaje del huso

BUB1B codifica BUBR1, que forma parte del complejo de checkpoint mitótico (MCC),
el freno que impide la separación de cromátidas hasta que todos los cinetocoros
están correctamente unidos al huso. El MCC se cataliza en el propio cinetocoro
(PMID 36518060), y BUB1 y BUBR1 actúan además como andamios de la corona fibrosa
del cinetocoro (PMID 40938979).

### 2. El dominio 766-1050 no es una quinasa: es un andamio para PP2A-B56

Anotado en UniProt como "Protein kinase", ese dominio es en realidad un
pseudoquinasa sin actividad catalítica. Su función demostrada es reclutar la
fosfatasa PP2A-B56 al cinetocoro, lo que estabiliza las uniones
cinetocoro-microtúbulo y permite el silenciamiento ordenado del checkpoint
(PMID 33207204, PMID 35525552).

**Esto es lo que hace relevante a p.Asn1002Lys:** cae dentro de ese dominio. La
hipótesis mecanística no es pérdida de catálisis —no hay catálisis que perder—
sino interferencia con el reclutamiento de PP2A-B56. Y esa hipótesis es
falsable, que es lo que la hace útil.

### 3. La dosis de BUBR1 determina la tasa de missegregación

No es un interruptor, es un reóstato. Niveles reducidos de BUB1B se asocian a
aneuploidía (PMID 18699967), y en la dirección contraria, aumentar BubR1 protege
frente a aneuploidía y cáncer y extiende la vida sana en ratón (PMID 23242215).

Esto encaja con la arquitectura genética observada: un alelo nulo con NMD
predicho más un segundo alelo que probablemente reduce función sin abolirla. La
pérdida bialélica completa no es viable; lo que se ve en pacientes es dosis
residual baja.

### 4. Mutaciones bialélicas en BUB1B causan MVA1

Es el trabajo fundacional: aneuploidía constitucional y predisposición tumoral
por mutaciones bialélicas en BUB1B (PMID 15475955). El patrón de alelo truncante
más alelo hipomórfico está documentado en siete familias con separación
prematura de cromátidas y checkpoint defectuoso (PMID 16411201).

#### Por qué el segundo alelo tiene que conservar función: la esencialidad lo predice

Hay tres genes MVA con arquitectura descrita, y no se parecen entre sí. Cruzándolos
con la esencialidad medida en DepMap (1.258 líneas celulares, a través de Open
Targets; `scripts/drug_landscape.py`) aparece un patrón:

| Gen | Líneas que dependen de él | Arquitectura descrita en pacientes |
|---|---|---|
| CEP57 (MVA2) | **0%** · efecto mediano −0,08 | homocigoto truncante (PMID 30035751) |
| TRIP13 (MVA3) | **1%** · efecto mediano −0,24 | pérdida de función bialélica (PMID 28553959) |
| **BUB1B** (MVA1) | **77%** · efecto mediano −1,21 | truncante más alelo hipomórfico (PMID 16411201) |

Donde el gen no es esencial, el doble nulo se tolera y es lo que se ve. Donde sí lo
es, la célula no puede dividirse sin él, y lo que se describe es un alelo nulo
acompañado de otro que conserva algo de función. gnomAD va en el mismo sentido:
ninguna de las dos variantes del paciente aparece en homocigosis en 1,46 millones
de alelos.

**El paciente encaja en la predicción para BUB1B:** p.Leu737Ter es nulo, con NMD
verificado, y p.Asn1002Lys es un missense, el tipo de cambio que típicamente reduce
función sin abolirla. Un segundo alelo truncante habría sido, bajo este modelo,
más difícil de explicar que el que tiene.

Límites de este razonamiento, que conviene decir antes que el panel: son tres
genes, no una ley; DepMap mide dependencia en líneas tumorales que proliferan, no
en un embrión; y no aporta ningún criterio ACMG. Es coherencia mecanística, no
evidencia de patogenicidad del missense.

### 5. Checkpoint débil → missegregación → aneuploidía en mosaico

Si el freno mitótico es insuficiente, los errores de segregación se acumulan
célula a célula durante el desarrollo. Como cada error es independiente, el
resultado no es una trisomía uniforme sino un mosaico de cromosomas distintos en
proporciones distintas: literalmente aneuploidía variegada.

### 6. Aneuploidía en mosaico → fenotipo del desarrollo

El cuadro clínico de MVA incluye microcefalia y retraso del crecimiento, con
heterogeneidad importante entre pacientes (PMID 18548531). La microcefalia no es
obligatoria para el diagnóstico (PMID 16059936), lo que importa a la hora de
emparejar fenotipos: exigirla como criterio pierde pacientes.

### 7. Inestabilidad cromosómica → predisposición tumoral

MVA cursa con riesgo tumoral en la infancia: rabdomiosarcoma embrionario
(PMID 9916837) y tumor de Wilms, descritos en el mismo paciente de forma
secuencial (PMID 42595739).

### 8. La célula aneuploide vive bajo estrés, y ese estrés es explotable

El desbalance estequiométrico de proteínas impone estrés proteotóxico
(PMID 27308438), y la disrupción de proteostasis por aneuploidía deteriora la
función mitocondrial (PMID 40527892). Los micronúcleos que genera la
missegregación activan la vía cGAS-STING (PMID 29342134), que en cánceres
cromosómicamente inestables sostiene una señalización inflamatoria dependiente de
IL-6 (PMID 35705809).

#### Lo que muestra el transcriptoma de una célula sin BUB1B

Esa rama inflamatoria se puede mirar con datos públicos. La firma consenso de
knockout CRISPR de BUB1B en LINCS L1000 (247 genes que suben y 244 que bajan),
enriquecida contra MSigDB Hallmark y comparada con las otras 5.211 firmas de
knockout de la misma librería (`scripts/cmap_reversal.py`, `cmap_signature.tsv`):

| Programa | Sentido | p | Knockouts con igual o más enriquecimiento |
|---|---|---|---|
| Respuesta a interferón gamma | sube | 4e-18 | 2,5% |
| Respuesta a interferón alfa | sube | 4e-17 | 3,1% |
| TNF-alfa vía NF-kB | sube | 3e-17 | 2,3% |
| Apoptosis | sube | 6e-13 | 0,1% |
| Vía de p53 | sube | 3e-09 | 0,9% |
| Checkpoint G2-M | baja | 1e-12 | 0,8% |

La respuesta a interferón no es lo que le pasa a cualquier gen esencial cuando se
apaga: solo el 3% de los knockouts la muestran con esa fuerza. Y separa al módulo de
señalización del checkpoint del resto de la mitosis. La comparten BUB1 (p 1e-12),
ZWINT (1e-12) y, más débil, TTK (0,005); no la muestran AURKB, PLK1, CENPE ni
NDC80. Es lo que predice la cadena: perder la señal del checkpoint produce
missegregación y micronúcleos, mientras que perder AURKB o PLK1 detiene la mitosis.

**Lo que no sostiene, dicho antes de que lo diga el panel:**

- **No se replica con shRNA.** En el knockdown por shRNA, en las 2 líneas donde se
  verificó que BUB1B baja, la respuesta a interferón no aparece (0 genes). Hay dos
  lecturas, y estos datos no permiten elegir: que la respuesta depende de la
  técnica, o que depende de la dosis. Un knockout pierde la proteína entera; un
  knockdown deja proteína residual, que es lo más parecido al paciente. Si es la
  segunda, la rama inflamatoria pesa menos en este niño de lo que la literatura
  tumoral sugiere.
- **El interferón no prueba cGAS-STING.** Es compatible, no específico. La prueba
  es el experimento de la tabla de abajo.
- **Son líneas tumorales.** La misma reserva del resto de este paso.

**Este es el punto de entrada para reposicionamiento:** no corregir la variante,
sino actuar sobre las consecuencias celulares que el mecanismo predice.

---

## Dónde la cadena es más débil

Un reporte que solo enumera evidencia a favor no es un reporte, es un alegato.

1. **El paso 2 es una hipótesis, no una observación.** Que el dominio recluta
   PP2A-B56 está demostrado; que *esta variante concreta* lo altera, no. Es VUS
   por ACMG y así se declara.
2. **Los pasos 8 vienen de biología tumoral.** Estrés proteotóxico, cGAS-STING e
   IL-6 están caracterizados en líneas y tumores con CIN, no en un niño con
   aneuploidía constitucional. Extrapolar de una célula tumoral seleccionada por
   tolerar aneuploidía a un tejido del desarrollo es un salto, y lo nombramos.
3. **La fase no está resuelta.** Todo el argumento asume configuración en trans.
4. **El WGS de sangre no muestra la aneuploidía.** Nuestro barrido de BAF no
   detecta eventos por encima de ~10-15% de fracción celular. La aneuploidía de
   MVA se documenta por cariotipo en células cultivadas; la ausencia en sangre no
   contradice el diagnóstico, pero tampoco lo apoya de forma independiente.

---

## Qué predice el mecanismo, en cosas medibles

Un mecanismo que no predice nada medible no sirve para diseñar experimentos.

| Predicción | Cómo se mide |
|---|---|
| Separación prematura de cromátidas en linfocitos | cariotipo con conteo de PCS |
| Aneuploidía variegada dependiente de tejido | cariotipo en fibroblasto contra linfocito |
| Checkpoint débil ante estrés del huso | índice mitótico tras nocodazol en células del paciente |
| Micronúcleos aumentados | conteo de micronúcleos por inmunofluorescencia |
| Reclutamiento de PP2A-B56 alterado | inmunoprecipitación del cinetocoro en célula que expresa el missense |
| Firma inflamatoria por cGAS-STING | expresión de genes de interferón tipo I e IL-6 |
| Firma de pérdida de BUB1B, y en qué grado | RNA-seq de fibroblastos del paciente puntuado contra los 247/244 genes de `cmap_signature.tsv` |

Las filas de PP2A-B56 y cGAS-STING son las que discriminan la hipótesis específica
de p.Asn1002Lys frente a la genérica de "menos BUBR1". La última fila mide dosis. Si
las células del paciente reproducen la firma del knockout, con su respuesta a
interferón, la rama inflamatoria cuenta en este niño. Si se parecen al knockdown,
sin interferón, cuenta poco, y el candidato senolítico gana peso relativo frente a
cualquier idea antiinflamatoria.

---

## Referencias

Verificadas contra PubMed, 2026-09-12.

| PMID | Qué aporta |
|---|---|
| 36518060 | Formación del MCC catalizada en el cinetocoro |
| 40938979 | BUB1 y BUBR1 como andamios de la corona fibrosa |
| 33207204 | El dominio pseudoquinasa de BUBR1 recluta PP2A-B56 |
| 35525552 | Fosforilación del pseudoquinasa BUBR1 |
| 18699967 | Nivel reducido de BUB1B asociado a aneuploidía |
| 23242215 | Más BubR1 protege de aneuploidía y cáncer |
| 15475955 | Mutaciones bialélicas en BUB1B causan MVA1 |
| 16411201 | Siete familias, alelo monoalélico y checkpoint defectuoso |
| 23812934 | BubR1 bajo: la célula no sostiene la parada ante venenos del huso (verificada 2026-10-07) |
| 30035751 | Homocigoto truncante en CEP57 en MVA2 (verificada 2026-09-12) |
| 28553959 | Pérdida de función bialélica en TRIP13, MVA3 y Wilms (verificada 2026-09-30) |
| 18548531 | Heterogeneidad clínica y genética de MVA |
| 16059936 | La microcefalia no es obligatoria en MVA |
| 9916837 | MVA con rabdomiosarcoma embrionario |
| 42595739 | Wilms y rabdomiosarcoma secuenciales en MVA |
| 27308438 | Aneuploidía y estrés proteotóxico |
| 40527892 | Proteostasis y función mitocondrial en aneuploidía |
| 29342134 | CIN y respuesta a DNA citosólico |
| 35705809 | cGAS-STING e IL-6 en cánceres con CIN |
