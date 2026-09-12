# Track 1 — Methods

Rare Disease, Real Kid: The MVA Hackathon 2026.

## Resumen

El candidato es un **heterocigoto compuesto en BUB1B (MVA1)**:

- **p.Leu737Ter**, nonsense en el exón 17 de 23, con degradación por NMD predicha.
  Clasificada Pathogenic/Likely_pathogenic en ClinVar por múltiples submitters sin
  conflictos, asociada a *Mosaic variegated aneuploidy syndrome 1*.
- **p.Asn1002Lys**, missense en el dominio 766-1050, presente en 1 de 1.461.878
  alelos de gnomAD v4.1 y ausente de ClinVar.

Las dos se entregan como par en una sola fila, porque el scoring da crédito
completo solo si aparecen ambas.

## Datos y referencia

- VCF entregado, 5.012.204 variantes, muestra única, llamadas con GATK
  (VariantFiltration en el header), profundidad mediana 44x, sexo genético
  masculino.
- Build GRCh38 con 2580 contigs **sin prefijo `chr`**, del analysis set con
  decoys hs38d1. El formulario de submission pide `chr15`, así que la conversión
  de nomenclatura es un paso explícito con test unitario.
- No se usaron los 80 GB de FASTQ. El dataset no incluye BAM ni CRAM, de modo que
  realinear habría sido reprocesar desde cero sin beneficio para la pregunta:
  el VCF entregado ya contiene la respuesta.
- Recursos de anotación: MANE Select v1.5 (estructura exónica y secuencia de
  transcrito), gnomAD v4.1 (exomas y genomas), ClinVar, UniProt O60566, referencia
  de chr15 de UCSC.

## Pipeline

Todo en Python de la biblioteca estándar. Sin contenedores, sin clúster, sin GPU.

| Paso | Script | Qué hace |
|---|---|---|
| Perfil del VCF | `vcf_quicklook.py` | metadatos, conteos y balance alélico por cromosoma |
| Aneuploidía | `vcf_baf_scan.py` | dispersión del BAF contra lo esperado por binomial, ventanas de 10 Mb |
| Barrido dirigido | `gene_sweep.py` | reparte variantes entre 32 genes del checkpoint según MANE |
| Frecuencias | `annotate_gnomad.py` | gnomAD por rango de bytes remoto, sin descargar los archivos |
| Consecuencia y fase | `candidates.py` | exón/intrón/splicing contra MANE, y fase por bloques PID |
| ClinVar | `clinvar_annot.py`, `clinvar_genomewide.py` | cruce dirigido y cruce genómico ciego |
| Proteína | `protein_consequence.py` | traducción local con control de base de referencia |
| Splicing | `splice_screen.py` | sitios crípticos con consenso y tasa de fondo medida |
| Fase | `phase_analysis.py` | alcance empírico del phasing por lecturas |
| Submission | `make_submission.py` | formato, conversión de contigs, validación y réplica del scoring |

### Una decisión de arquitectura que resultó importante

gnomAD publica cada cromosoma como un archivo de varios GB junto a un índice
tabix de unos 50 KB, y el servidor acepta rangos de bytes. En vez de descargar
7.4 GB por cromosoma, se implementó un cliente tabix por HTTP en Python puro
(`tabix_remote.py`): baja el índice, calcula los bloques BGZF que cubren cada gen
y pide solo esos. Se traen unos pocos MB.

Además de ahorrar ancho de banda, esto resuelve un problema de datos: la consulta
que sale de la máquina es el rango de un gen del panel, no una posición del
paciente.

## Priorización

Modelo autosómico recesivo. Se busca, por gen, un homocigoto raro o dos o más
heterocigotas raras que puedan estar en trans.

- **Cero homocigotos raros** en el panel de 32 genes.
- 14 genes con dos o más variantes raras.
- **Un solo gen con dos variantes raras codificantes: BUB1B.** El único otro gen
  con alguna rara codificante es SGO1, con una.

## Ranking ciego frente a ranking con prior

Requisito autoimpuesto: el pipeline tiene que poder llegar al gen sin que se le
diga dónde mirar.

El cruce genómico completo contra ClinVar, sin prior de genes, devuelve **7
genotipos P/LP en todo el genoma, ninguno homocigoto**: BUB1B, FLG, GNRHR, HK1,
LZTR1, PRSS1 y RBM8A. BUB1B es el único que además tiene una segunda variante
rara codificante en el mismo gen.

*Límite:* este ranking ciego se apoya en ClinVar, así que solo encuentra lo ya
clasificado. Una variante truncante nueva en un gen sin entradas no aparecería. Un
ranking ciego completo exige anotación de consecuencia y frecuencia a escala
genómica, que no se corrió.

## Clasificación ACMG/AMP

Detalle criterio por criterio en `ACMG_BUB1B.md`, incluidos los criterios que se
decidieron **no** aplicar y por qué.

- p.Leu737Ter: PVS1 (NMD verificado contra la estructura exónica, no asumido) +
  PM2_Supporting → **Patogénica**.
- p.Asn1002Lys: PM2_Supporting como único criterio → **VUS**.

No se aplicaron PM1, PP2 ni PP3, con razón escrita en cada caso. Conviene destacar
uno: UniProt anota el dominio 766-1050 como "Protein kinase", pero es un
**pseudoquinasa** sin actividad catalítica cuyo papel es reclutar PP2A-B56
(PMID 33207204). Deducir pérdida de catálisis habría sido un error de
razonamiento.

**Incluso resolviendo la fase, el segundo alelo sigue siendo VUS.** PM2 más PM3 no
alcanza para probablemente patogénica; hace falta evidencia funcional.

## Lo que no se pudo resolver

**Fase.** Los bloques de fase PID de este VCF en chr15 tienen mediana de 8 pb, p99
de 82 pb y máximo de 206 pb sobre 8696 bloques. Las dos candidatas están a
**10.911 pb, 53 veces el bloque más largo observado**. Ninguna lectura ni ningún
par de este experimento puede abarcar ambas posiciones: la configuración en trans
es una inferencia, no una observación. Lo resolverían genotipar a los padres,
lectura larga o PCR de alelo específico.

**Aneuploidía.** El barrido de BAF no detecta eventos en mosaico ni a escala de
cromosoma ni segmentaria a 10 Mb. A 44x y con ~150 heterocigotos por ventana, el
método ve fracciones celulares por encima de 10-15%. La aneuploidía variegada se
documenta por cariotipo en células cultivadas; su ausencia en el WGS de sangre no
contradice el diagnóstico, y el límite de detección se reporta en vez de omitirse.

**Splicing.** El cribado propio de sitios crípticos es negativo en las tres
variantes de BUB1B, con tasa de fondo medida entre 6.5% y 7.9%. No sustituye a
SpliceAI ni a Pangolin, que quedan pendientes.

**Cobertura.** El barrido solo ve lo que el VCF reporta. Una deleción homocigota o
una región sin cobertura no produce líneas y es invisible. TRIP13, uno de los
genes MVA, aparece con profundidad media de 20.7x frente a 44x del genoma, así
que "no hay variantes en TRIP13" no es concluyente.

## Dos errores de método, encontrados y corregidos

Se documentan porque el proceso importa tanto como el resultado.

1. **Test de lóbulos de BAF.** La primera versión partía la distribución en 0.5 y
   llamaba lóbulos a las medianas de cada mitad. Eso devuelve dos lóbulos incluso
   en una distribución unimodal, y produjo fracciones celulares idénticas
   (0.22 y 0.18) en los 24 cromosomas. La firma de un artefacto de método es
   precisamente esa uniformidad. Sustituido por comparación contra la dispersión
   esperada por binomial, con exigencia de segmentos consecutivos.
2. **Cribado de splicing.** La primera versión contaba cualquier dinucleótido
   creado de una lista de cinco, y daba positivo en las tres variantes. Como un
   cambio de base toca dos dinucleótidos y cinco de dieciséis están en la lista,
   ese test se dispara en más de la mitad de las variantes al azar. Sustituido por
   dinucleótido más consenso, con nulo empírico de 2000 sustituciones.

También se corrigió un error en el propio diseño del panel: "RZZ" no es un gen
sino un complejo, y sus componentes son ZW10, ZWILCH y KNTC1. Y ClinVar muestra
**MAD1L1 como *Mosaic variegated aneuploidy syndrome 7***, de modo que la lista
canónica de genes MVA (BUB1B, CEP57, TRIP13) está desactualizada.

## Estrategia de submission

El formato admite 10 filas y F-max toma el máximo sobre los umbrales de EPCR. De
ahí se sigue que **rellenar filas por debajo no baja el puntaje mientras el orden
sea correcto**, y que lo único que hunde el resultado es una variante equivocada
por encima de la verdadera. El trabajo está en el orden y en la separación de los
EPCR, no en el volumen. La réplica del scoring incluida en `make_submission.py`
verifica ambas cosas antes de escribir el archivo.

Se entregan 5 filas: el par de BUB1B con EPCR 0.95, dos candidatos de respaldo
claramente por debajo, y dos hallazgos secundarios en genes de acción dominante
(PRSS1 y LZTR1), marcados como `secondary` y anotados para revisión clínica.

`proband_id` se fija en `PROBAND01`: es el valor de la plantilla oficial y la
clave del gold standard. Si el scoring empareja por clave, usar el nombre de
muestra del VCF daría cero sin ningún mensaje de error.

## Reproducibilidad

```bash
python scripts/gene_sweep.py data/<vcf>            # barrido por panel
python scripts/annotate_gnomad.py                  # frecuencias, remoto por rangos
python scripts/candidates.py                       # consecuencia y fase
python scripts/clinvar_annot.py                    # cruce dirigido
python scripts/clinvar_genomewide.py               # ranking ciego
python scripts/protein_consequence.py              # traducción local
python scripts/make_submission.py in.tsv out.csv --proband PROBAND01
python tests/test_submission.py                    # 12 casos, incluida la réplica del scoring
python tests/test_egress_guard.py                  # 30 canarios, 10 consultas legítimas
```

Requisitos: Python 3.12+, biblioteca estándar. Corre en un portátil de 4 núcleos y
8 GB de RAM. Coste de cómputo: 0 dólares.

## Divulgación de uso de IA

El análisis lo ejecutaron scripts deterministas de Python. Un modelo de lenguaje
(Claude, Anthropic) escribió esos scripts, interpretó sus salidas y redactó los
informes.

**Qué recibió el modelo:** símbolos de genes, consecuencias a nivel de proteína,
métricas de calidad, frecuencias poblacionales, clasificaciones de ClinVar y
estadísticas agregadas por cromosoma.

**Qué no recibió:** coordenadas genómicas del paciente, tablas de genotipos,
contenido del VCF, ni identificadores de muestra. Los scripts leyeron esos datos
en local y escribieron los resultados a disco; a la salida estándar, y por tanto
al contexto del modelo, solo llegaron agregados y hallazgos a nivel de gen o de
proteína, que es lo que las reglas del hackathon permiten conservar y publicar.

Un hook `PreToolUse` (`egress_guard.py`) bloquea de forma automática cualquier
llamada a un servicio externo cuyo contenido incluya coordenadas, rsIDs, HGVS con
posición, identificadores de muestra o rutas de datos del paciente. Su prueba de
red-team bloquea 30 de 30 canarios y deja pasar 10 de 10 consultas legítimas a
nivel de gen.

## Manejo y borrado de datos

Los datos viven en un único árbol purgable, fuera de cualquier carpeta
sincronizada a la nube. `scripts/purge.py` borra el VCF y sus derivados, los
índices, las tablas de genotipos y **los transcripts locales de la sesión de
Claude Code**, que contienen bloques de datos y entran en la lista de borrado de
las reglas. Modo `--verify` comprueba que ningún artefacto conservado tenga
aspecto de tabla de genotipos. La atestación al cierre la envía una persona.
