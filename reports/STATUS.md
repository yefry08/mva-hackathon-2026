# STATUS

**Repositorio público:** https://github.com/yefry08/mva-hackathon-2026
(commit inicial `daed07a`, 45 archivos, 0 archivos genómicos verificado contra
la API de GitHub tras publicar).

Actualizado: 2026-09-12. Lo reescribe el agente de convergencia cada 2h una vez que
arranque la corrida; por ahora lo mantiene el orquestador a mano.

## Dónde estamos

Fase 0 (gates e infraestructura) en curso. Todavía no hay corrida multiagente: falta
la API key comercial para que los workers puedan leer datos a nivel de variante.

## Hecho

- **Reglas leídas y documentadas** en `RULES-FINDINGS.md`. Lo decisivo: los LLM de
  terceros están permitidos bajo la prueba de procesador (sin entrenamiento sobre
  inputs ni outputs, retención limitada), hay que declarar proveedor y plan en el
  writeup, y **no se pueden calificar outputs**.
- **Workspace en `C:\mva`**, fuera de OneDrive. Git inicializado, canary del
  `.gitignore` pasando.
- **Egress-guard** (`.claude/hooks/egress_guard.py`): 30/30 canarios bloqueados,
  10/10 consultas a nivel de gen permitidas (`tests/test_egress_guard.py`).
- **Purga** (`scripts/purge.py`): simulacro, `--apply` y `--verify`. Cubre `data/`,
  `work/`, `logs/` y los transcripts locales de Claude Code de este proyecto.
- **Formateador de submission** (`scripts/make_submission.py`): conversión de
  contigs, validación del formato y réplica del scoring del hackathon.
  `tests/test_submission.py` en verde (12 casos).
- **Datos**: descargados VCF (315 MB), índice y el `.docx` de fenotipo. Los 80 GB de
  FASTQ no se bajaron: sin BAM/CRAM en el dataset, realinear es un proyecto aparte.

## Hallazgos sobre el caso (todo agregado, ningún genotipo individual)

| Qué | Valor |
|---|---|
| Variantes | 5.01 M, muestra única, llamadas con GATK |
| Profundidad | ~44x, uniforme entre autosomas (razón 1.00–1.07) |
| Build | GRCh38, 2580 contigs, **sin prefijo `chr`** |
| Fase | **494.524 variantes con PGT en 254.655 bloques PID** |
| Sexo genético | Masculino (X con 8.5k hets de 117k; Y presente) |

**Aneuploidía a escala de cromosoma completo: negativa** con este filtrado. La
profundidad elevada de chr21/22/16 del primer vistazo desapareció al exigir PASS,
bialélico y DP≥20: era ruido de regiones repetitivas.

**Aneuploidía segmentaria en mosaico: negativa** en el barrido por ventanas de
10 Mb (`vcf_baf_scan.py`). El exceso de dispersión del BAF respecto de la binomial
tiene mediana 1.087 y p95 1.91 entre las 293 ventanas autosómicas. Las 34 ventanas
que pasan el umbral son picos aislados y caen todas en brazos acrocéntricos (13p,
14p, 15p, 21p, 22p), centrómeros (2, 10, 16, 17, 20) y bloques de heterocromatina
(9q12, 1q12, 1q21). Ninguna forma un segmento de 3 ventanas seguidas.

**Límite de detección, que es el dato para el writeup:** a 44x y ~150 hets por
ventana, este barrido ve eventos por encima de una fracción celular de 10-15%.
Un mosaicismo por debajo de eso, o presente en otro tejido, no aparecería. La
aneuploidía variegada de este niño se documenta por cariotipo en células
cultivadas; que el WGS de sangre no la vea no la contradice, y decir esto con el
límite cuantificado es más fuerte que no reportarlo.

**Corrección de método registrada.** El primer test de BAF (`vcf_aneuploidy.py`)
partía la distribución en 0.5 y llamaba lóbulos a las medianas de cada mitad. Eso
devuelve dos lóbulos incluso en una distribución unimodal, y produjo fracciones
celulares idénticas (0.22/0.18) en los 24 cromosomas: artefacto del método, no del
paciente. Sustituido por `vcf_baf_scan.py`, que compara la dispersión observada del
BAF contra lo que predice la binomial a esa profundidad, por ventanas de 10 Mb, y
solo llama evento cuando BAF y cobertura apuntan al mismo lado.

## Barrido de los 32 genes del checkpoint

`scripts/gene_sweep.py` reparte las variantes del VCF entre los genes del panel
usando las coordenadas del transcrito MANE Select (v1.5). 2369 variantes no de
referencia, 590 genotipos homocigotos alternos. El detalle por variante vive en
`work/gene_sweep.tsv` y no sale de ahí.

Dos huecos de cobertura que importan, porque un hueco esconde la segunda variante
de un compuesto:

| Gen | DP media | Contexto |
|---|---|---|
| TRIP13 (MVA3) | 20.7x | menos de la mitad de los 44x del genoma |
| BUB1B (MVA1) | 35.9x | por debajo de la media |
| PLK4 | 29.0x | una sola variante en todo el gen |

**Límite del método, para el writeup:** este barrido solo ve lo que el VCF
reporta. Una deleción homocigota o una región sin cobertura no produce líneas y
por lo tanto es invisible aquí. La evidencia negativa por gen ("no hay variantes
en X") no es concluyente sin datos de cobertura, y el dataset no trae BAM/CRAM ni
gVCF. Decirlo es parte del rigor; callarlo sería afirmar más de lo que se sabe.

## Anotación de frecuencias sin descargar gnomAD

`scripts/tabix_remote.py` implementa un cliente tabix por HTTP en Python puro:
descarga el índice `.tbi` (51 KB), calcula los rangos de bytes y pide solo los
bloques BGZF de los genes del panel. Probado contra gnomAD v4.1: 2781 registros
para una ventana de 20 kb de BUB1B, sobre un archivo de 7.4 GB que nunca se baja.

Esto elimina tres dependencias del plan original: no hace falta Docker, ni AWS, ni
descargar recursos de anotación. Y lo que sale de la máquina es la coordenada de
un gen del panel, nunca una posición del paciente.

## Candidato de Track 1

**BUB1B, heterocigoto compuesto.**

| Alelo | Consecuencia | Evidencia |
|---|---|---|
| 1 | p.Leu737Ter (nonsense) | ClinVar Pathogenic/Likely_pathogenic, múltiples submitters sin conflictos, anotada a *Mosaic variegated aneuploidy syndrome 1*. DP 46, GQ 99 |
| 2 | p.Asn1002Lys (missense) | rara en gnomAD, ausente de ClinVar. DP 28, GQ 99 |

Arquitectura consistente con lo publicado en MVA1: un alelo truncante más uno
missense. Las consecuencias proteicas se calcularon en local contra el MANE
Select, y el control de mapeo (base de referencia del VCF contra la del
transcrito) pasó en las tres variantes del gen.

**Ranking ciego, sin prior de genes.** El cruce genómico completo contra ClinVar
encuentra 7 genotipos P/LP, ninguno homocigoto: BUB1B, FLG, GNRHR, HK1, LZTR1,
PRSS1, RBM8A. BUB1B es el único que además tiene una segunda variante rara
codificante en el mismo gen. El hallazgo no depende del prior de MVA.

*Límite de ese ranking:* se apoya en ClinVar, así que solo encuentra lo ya
clasificado. Ese hueco se cerró con dos barridos que no usan ClinVar ni panel.

**Ranking ciego definitivo** (`compound_genomewide.py`): sobre todo el genoma, sin
panel ni ClinVar, 3 genes cargan una truncante rara más un missense raro:
SLC25A5, BUB1B y HLA-DRB1. SLC25A5 está en el cromosoma X y el paciente es varón,
así que un compuesto es imposible; HLA-DRB1 está en la región HLA. **BUB1B es el
único que se sostiene**, y es también el único con el conteo que predice la
biología: una variante por alelo. Hueco que sigue declarado: un compuesto de dos
missense en un gen sin entradas en ClinVar.

**Lo que no está resuelto: la fase.** Cuantificado, no supuesto. Los bloques de
fase PID de este VCF en chr15 tienen mediana de 8 pb, p99 de 82 pb y **máximo de
206 pb** sobre 8696 bloques. Las dos candidatas están a **10.911 pb: 53 veces el
bloque más largo observado**. Ninguna lectura ni ningún par de este experimento
puede abarcar ambas posiciones, así que el phasing por lecturas no falla por cómo
se corrió el análisis: la información no está en los datos. Lo resolverían
secuenciar a los padres, lectura larga, o PCR de alelo específico. Se entrega el
par completo, que es lo que el scoring premia, y el writeup lo dice así.

**Cribado de splicing: negativo en las tres variantes de BUB1B.** Ninguna crea un
sitio críptico plausible, entendiendo por plausible un GT con consenso donador de
al menos 6/9 o un AG con tracto de polipirimidinas de al menos 65%. El cribado
lleva su propio nulo: muestreando 2000 sustituciones al azar del entorno, entre
6.5% y 7.9% crean un sitio así por azar. La intrónica profunda queda con prior
bajo y se mantiene en la lista solo como fila de respaldo.

*Segunda corrección de método registrada.* La primera versión de este cribado
contaba cualquier dinucleótido creado de una lista de cinco, y daba positivo en
las tres variantes. Con un cambio de base tocando dos dinucleótidos y cinco de
dieciséis en la lista, ese test se dispara en más de la mitad de las variantes al
azar: no discriminaba nada. Se reemplazó por dinucleótido más consenso, con tasa
de fondo medida.

## Submission preparada (NO enviada)

`submissions/track1_v1.csv`, 5 filas, valida contra el formato oficial:

| # | EPCR | tipo | qué |
|---|---|---|---|
| 1 | 0.95 | primary | par compuesto BUB1B |
| 2 | 0.08 | primary | SGO1 missense rara |
| 3 | 0.05 | primary | BUB1B intrónica profunda, ausente de gnomAD |
| 4 | 0.02 | secondary | PRSS1, alelo P/LP conocido, gen de acción dominante |
| 5 | 0.01 | secondary | LZTR1, alelo P/LP conocido, gen de acción dominante |

`proband_id` = `PROBAND01`: es el valor de la plantilla oficial y la clave del
gold standard. Si el scoring empareja por clave, el nombre de muestra del VCF
daría cero sin mensaje de error; `PROBAND01` gana en los dos escenarios.

## Track 2

**Mecanismo** (`reports/MECHANISM_BUB1B.md`): cadena de ocho eslabones desde la
variante hasta el fenotipo, con 16 citas verificadas contra PubMed, más una
sección explícita de dónde la cadena es más débil. El giro conceptual: el dominio
766-1050 que UniProt llama "Protein kinase" es un pseudoquinasa sin actividad
catalítica cuyo papel es reclutar PP2A-B56 (PMID 33207204). La hipótesis para
p.Asn1002Lys no es pérdida de catálisis sino interferencia con ese andamiaje, que
es falsable en el banco.

**Reposicionamiento** (`reports/REPURPOSING.md`): tres candidatos con ficha de
evidencia a favor, evidencia en contra, seguridad pediátrica y experimento
falsable. Confianzas 0.35, 0.25 y 0.15, bajas a propósito porque ninguno tiene un
solo dato en MVA.

1. Eje senolítico — el más fuerte, porque la intervención se probó en el ratón
   progeroide **BubR1** (PMID 22048312), el modelo del gen correcto. Debilidad:
   fue una herramienta genética inducible, no un fármaco, y en un organismo
   envejecido, no en desarrollo.
2. Sirolimus — recambio proteico y tolerancia a aneuploidía. Contraargumento
   específico: inhibir mTOR frena el crecimiento, y el retraso del crecimiento ya
   es parte del fenotipo.
3. N-acetilcisteína — carga oxidativa. Vínculo indirecto, confianza baja.

**Descartados con razón escrita:** inhibidores de Aurora y de MPS1/TTK, y
antimitóticos en general. Aparecen en toda búsqueda de aneuploidía porque se usan
*contra* tumores con CIN.

**Precedente:** no existe ningún intento previo de reposicionamiento en MVA. Lo
único con respaldo hoy son los protocolos de vigilancia oncológica (PMID
39264246), y el reporte lo dice antes de proponer nada.

**Pendiente y declarado:** Robin, LINCS/CMap y DepMap no se corrieron por falta de
credenciales.

## Barrido genomico de LoF bialelica (sin ClinVar)

`lof_genomewide.py` traduce en local las 23.041 variantes codificantes del genoma
contra su transcrito MANE y marca las truncantes. 252 genes con alguna LoF, 104
con LoF bialelica potencial, y tras filtrar por frecuencia en gnomAD quedan **5**:
AGAP3, SERPINA1, ADAMTS1, HLA-DQA1 y POU6F2. Ninguno compite con BUB1B: cuatro
caen en regiones donde el alineamiento de lecturas cortas falla (paralogos, HLA)
o presentan conteos implausibles de truncantes raras, y POU6F2 se anota solo
porque tiene mutaciones germinales descritas en Wilms (PMID 15459955), que forma
parte del espectro tumoral de MVA.

Limite declarado: este barrido solo ve LoF mas LoF, asi que el propio candidato de
BUB1B (truncante mas missense) no aparece en el. Complementa al analisis por
panel; no lo sustituye.

## Bloqueos

1. **API key comercial de Anthropic** en `C:\mva\.env`. Sin eso, ningún agente lee
   datos a nivel de variante y la fase de priorización no arranca.
2. FutureHouse/Edison sin credencial ni términos verificados: Track 2 se queda en
   literatura a nivel de gen hasta entonces.

## Presupuesto

- Submissions: Track 1 **0/5 usadas** (`config.py` dice 6, la pestaña dice 5).
  Track 2 0/3.
- LLM: sin corrida multiagente todavía.
- Cómputo: 0 dólares. Todo se ha hecho en este laptop con Python puro.

## Lo próximo

1. Barrido de BAF por ventanas (corriendo).
2. Reparar `dispatch.py` y regenerar la cola con el plan actualizado.
3. Con la key: fast-path recesivo sobre el VCF y barrido de los 32 genes del SAC.
