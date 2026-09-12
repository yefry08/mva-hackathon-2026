# PLAN — cola de tareas, presupuesto y supuestos frágiles

Estado: **propuesta, pendiente de visto bueno**. Nada de esto está en `queue.jsonl` todavía.

Leyenda: `C` = job de cómputo (Snakemake, sin LLM) · `A` = agente `claude -p` · `H` = gate humano.
Peso de agente: `l` liviano · `m` medio · `h` pesado (ver §3).

**Total: 330 tareas = 98 C + 226 A + 6 H.** De las 226 A, 40 son slots dinámicos (Z) que el
agente de convergencia instancia cuando existen los candidatos. No hay tareas de relleno: las 24
tareas "una por cromosoma" del scaffold anterior se colapsaron en una sola (K09).

---

## Fase 0 — Gates e infraestructura (22)

| ID | Tipo | Tarea |
|---|---|---|
| G01 | H | Leer reglas/DUA del hackathon: ¿permite procesar T1 con Claude (API Anthropic) y en AWS? Bloquea toda tarea A que lea T1. |
| G02 | H | Visto bueno de este plan. |
| G03 | C | Provisionar EC2 Linux (EBS cifrado ≥1.5 TB, sin puertos de entrada, acceso por SSM). |
| G04 | C | Repo fuera de OneDrive; `git init`; canary test del `.gitignore` (`data/canary.vcf.gz` invisible a git). |
| G05 | A·h | Hook `PreToolUse` egress-guard: bloquea `chr:pos`, rsID, HGVS con posición, IDs de muestra, rutas `data/`/`work/` en MCP/WebFetch/WebSearch/curl. |
| G06 | A·m | Red-team del egress-guard: 30 canarios deben bloquearse (30/30) y 10 queries a nivel gen deben pasar. |
| G07 | A·m | Leer las páginas Submit Track 1/2 → `submission-formatter/spec.md` + validador de formato. |
| G08 | C | Inventario de `data/`: md5, tipo y tamaño (ningún LLM lee contenido). |
| G09 | C | Header del VCF: build, contigs, nº de muestras (¿trío?), caller de origen. |
| G10 | C | Headers FASTQ/CRAM: plataforma, longitud de lectura, lanes, read groups. |
| G11 | A·m | Descripción clínica → términos HPO (local, sin APIs). |
| G12 | A·m | Segundo agente independiente sobre HPO; diff y reconciliación. |
| G13 | C | GRCh38 no-alt: descarga, `faidx`, dict, índice bwa-mem2. |
| G14 | C | Pull de contenedores pinneados (DeepVariant, GATK4, Manta, GRIDSS, MosaicHunter, DeepMosaic, MosaicForecast, MoChA, VEP, Exomiser, LIRICAL) → `versions.lock`. |
| G15 | A·m | Esqueleto Snakemake con las reglas C de este plan + perfil AWS. |
| G16 | A·m | Smoke test de futurehouse-mcp (1 query Crow a nivel gen); medir costo y latencia. |
| G17 | A·m | Instalar Robin, LLM vía LiteLLM, smoke run "Mosaic variegated aneuploidy syndrome"; medir costo. |
| G18 | A·m | Smoke test del MCP PubMed/Europe PMC: 10 PMIDs reales deben resolver y 3 inventados deben fallar. |
| G19 | A·h | Reparar `dispatch.py`: DAG real, pool continuo, `--allowedTools`, validación del contrato, tracker de costo, tipos C/A/H. |
| G20 | C | Dry-run del DAG: cero ciclos, cero huérfanos, orden topológico. |
| G21 | A·m | Fast-path: filtro recesivo sobre el VCF entregado → baseline de candidatos en la hora 2. |
| G22 | C | Balance alélico por cromosoma + ROH desde el VCF entregado (test rápido de aneuploidía y consanguinidad). |

## Fase 1 — Recursos de anotación (12, todas C)

R01 caché VEP GRCh38 + LOFTEE + UTRannotator · R02 dbNSFP · R03 CADD v1.7 (SNV+indel) ·
R04 SpliceAI precomputado + modelo para variantes nuevas · R05 AlphaMissense hg38 ·
R06 gnomAD v4.1 sites (exomas+genomas) + constraint · R07 ClinVar (release fechado) ·
R08 gnomAD-SV + DGV + datos AnnotSV · R09 HPO + Orphanet (+ OMIM genemap2 si hay licencia) ·
R10 datos Exomiser + LIRICAL · R11 GTEx mediana TPM ·
R12 PoN público GATK + blacklist ENCODE + mappability + estratificaciones GIAB.

## Fase 1 — QC (10)

| ID | Tipo | Tarea |
|---|---|---|
| Q01 | C | FastQC por FASTQ |
| Q02 | C | MultiQC consolidado |
| Q03 | C | mosdepth: media, uniformidad, % ≥20x |
| Q04 | C | VerifyBamID2: contaminación |
| Q05 | C | Sexo genético por cobertura X/Y vs metadata |
| Q06 | C | somalier: nuestras llamadas vs VCF entregado (sample swap) |
| Q07 | C | Duplicados, insert size, complejidad de librería |
| Q08 | C | Ancestría (PCA vs 1000G) → población gnomAD de referencia |
| Q09 | C | Cobertura por exón en los 32 genes del barrido S |
| Q10 | A·m | Interpretar QC y dar go/no-go para alinear |

## Fase 2 — Alineamiento (9)

| ID | Tipo | Tarea |
|---|---|---|
| L01 | C | Benchmark: 1M pares con bwa-mem2 → extrapolar tiempo y costo (va primero) |
| L02 | C | bwa-mem2 por read group |
| L03 | C | Merge + MarkDuplicates |
| L04 | C | BQSR |
| L05 | C | CRAM + índices; borrar BAM intermedios |
| L06 | C | Métricas (CollectWgsMetrics) |
| L07 | C | DRAGMAP como control ortogonal (solo si sobra presupuesto) |
| L08 | A·m | Comparar BWA vs DRAGMAP |
| L09 | A·m | Informe de alineamiento para el methods |

## Fase 3 — Llamado germinal (10)

| ID | Tipo | Tarea |
|---|---|---|
| V01 | C | DeepVariant WGS (VCF + gVCF) |
| V02 | C | GATK HaplotypeCaller (scatter interno) |
| V03 | C | Consenso DV∩HC + set de discordantes |
| V04 | A·m | Revisión de discordantes en genes candidatos |
| V05 | C | Concordancia vs VCF entregado |
| V06 | C | WhatsHap: phasing basado en lecturas |
| V07 | C | ROH (bcftools roh) + estimación de consanguinidad |
| V08 | C | Scan de UPD por ROH de cromosoma completo |
| V09 | C | mtDNA (Mutect2 modo mitocondria) |
| V10 | C | ExpansionHunter (expansiones de repeticiones) |

## Fase 4 — Mosaico (14)

| ID | Tipo | Tarea |
|---|---|---|
| M01 | C | Mutect2 tumor-only + PoN + recurso germinal gnomAD |
| M02 | C | MosaicHunter single-sample |
| M03 | C | DeepMosaic |
| M04 | C | MosaicForecast (cuarto caller, validado en single-sample) |
| M05 | C | Filtros de artefacto (hebra, posición en lectura, homopolímeros, mappability, blacklist) |
| M06 | C | Histograma de VAF global y por cromosoma |
| M07 | C | Consenso ≥2 callers |
| M08 | C | Curva de poder de detección por VAF según la cobertura |
| M09 | C | Set VAF 5–35 % en genes SAC |
| M10 | C | Snapshots IGV (igv-reports) del top 20 |
| M11 | A·m | Revisión visual de los snapshots (T1, local) |
| M12 | C | Tasa de falsos positivos: mismo pipeline sobre GIAB HG002 |
| M13 | A·m | ¿Hay segundo golpe somático o reversión? |
| M14 | A·m | Informe de mosaico |

## Fase 5 — CNV / SV / aneuploidía (13)

| ID | Tipo | Tarea |
|---|---|---|
| K01 | C | Manta |
| K02 | C | GRIDSS |
| K03 | C | CNVnator |
| K04 | C | CNVkit modo WGS (cn.mops pide cohorte; se documenta el cambio) |
| K05 | C | Consenso SV + AnnotSV |
| K06 | C | MoChA: alteraciones cromosómicas en mosaico por BAF + LRR |
| K07 | C | Cobertura normalizada por cromosoma y brazo (100 kb, corregida por GC) |
| K08 | C | BAF por cromosoma (hets de alta calidad) |
| K09 | A·m | Fracción celular aneuploide estimada por cromosoma (los 24 en una sola tarea) |
| K10 | A·m | Figura + tabla genome-wide |
| K11 | A·m | Contrastar el patrón con la literatura de MVA (tejido, % de células) |
| K12 | A·m | SV/CNV en genes SAC (deleción + missense = compuesto) |
| K13 | A·m | Informe CNV/aneuploidía |

## Fase 6 — Anotación (12, todas C)

N01 VEP+LOFTEE · N02 dbNSFP · N03 CADD · N04 SpliceAI precomputado · N05 SpliceAI de novo
para variantes no precomputadas · N06 AlphaMissense · N07 gnomAD · N08 ClinVar ·
N09 constraint + GTEx · N10 UTR/regulatorio (UTRannotator, cCREs ENCODE) · N11 anotación SV ·
N12 merge → `work/annotated.parquet`.

## Fase 7 — Priorización y ACMG (20)

| ID | Tipo | Tarea |
|---|---|---|
| P01 | C | Exomiser (AR+AD, sin prior de genes) |
| P02 | C | LIRICAL |
| P03 | C | Phen2Gene (local) |
| P04 | C | Exomiser sin fenotipo (control: solo variante) |
| P05 | H | AMELIE: es servicio web. **Por defecto se omite**; solo con aprobación y enviando únicamente HPO + genes |
| P06 | A·m | Consenso de rankings (RRF) + concordancia |
| P07 | C | Filtro recesivo: hom + comphet, AF<0.001, popmax<0.005 |
| P08 | C | Detección de comphet + fase (WhatsHap) |
| P09 | A·m | Pares sin fase: argumentar cis/trans |
| P10–P14 | A·h | ACMG/AMP del top 50 en 5 lotes de 10, criterio por criterio |
| P15 | A·m | Ranking CIEGO |
| P16 | A·m | Ranking CON PRIOR (MVA + SAC) |
| P17 | A·m | Delta ciego vs prior, explicado |
| P18 | A·m | Candidatos a gen nuevo (fuera del prior, score alto) |
| P19 | A·m | Hallazgos secundarios (lista ACMG SF vigente) |
| P20 | A·m | Informe de priorización |

## Fase 7b — Barrido por gen (32, A·l)

Cada tarea: todas las variantes del gen (codificantes, splicing, UTR, SV/CNV), VAF, cobertura por
exón, huecos <20x y una **declaración explícita de evidencia negativa** ("sin variantes; exones
1–23 ≥25x" también es un resultado).

S01 BUB1B · S02 CEP57 · S03 TRIP13 · S04 MAD2L1BP · S05 BUB1 · S06 BUB3 · S07 MAD1L1 ·
S08 MAD2L1 · S09 TTK · S10 CDC20 · S11 KNL1 · S12 ZW10 · S13 KNTC1 · S14 ZWILCH · S15 SPDL1 ·
S16 NDC80 · S17 NUF2 · S18 SPC24 · S19 SPC25 · S20 CENPE · S21 CENPF · S22 AURKB · S23 INCENP ·
S24 CDCA8 · S25 BIRC5 · S26 PLK1 · S27 ESPL1 · S28 PTTG1 · S29 SGO1 · S30 CENPA · S31 PLK4 ·
S32 CHMP4C

(Corrige el scaffold anterior: "RZZ" es un complejo, no un gen → KNTC1/ZW10/ZWILCH. Se añade
MAD2L1BP, que es el compañero de TRIP13.)

## Fase 8 — Set sintético, calibración y ranking (16)

| ID | Tipo | Tarea |
|---|---|---|
| X01 | C | Descargar GIAB HG002 ~30x (público) |
| X02 | A·m | Diseñar el set sintético: variantes MVA de ClinVar + 20 controles AR de otros genes, en hom/comphet |
| X03 | C | BAMSurgeon: spike-in germinal (VAF 50/100) |
| X04 | C | BAMSurgeon: spike-in mosaico (VAF 10/20/35) |
| X05 | C | Spike-in de deleción + missense (compuesto SV/SNV) |
| X06 | C | Pipeline completo sobre el sintético |
| X07 | A·m | Recall/precisión por clase y por VAF |
| X08 | A·h | Score compuesto (fenotipo + ACMG + VAF + constraint + aneuploidía), pesos ajustados en sintético |
| X09 | A·h | Calibración isotónica → probabilidades |
| X10 | A·m | Barrido del tamaño de la lista vs F-max esperado |
| X11 | A·m | Chequeo de comphet: las dos variantes presentes |
| X12 | A·m | Ablaciones: sin mosaico / sin SV / sin prior |
| X13 | A·m | submission-formatter: generar T1 #1 (sanity check) |
| X14 | H | **Aprobar el envío T1 #1** |
| X15 | A·m | Methods writeup T1 |
| X16 | A·h | Reproducibilidad: `make all` desde cero en una instancia limpia |

## Escalabilidad (6)

| ID | Tipo | Tarea |
|---|---|---|
| E01 | C | 30 phenopackets AR (GA4GH phenopacket-store) + VCFs de fondo 1000G |
| E02 | C | Spike-in de la variante de cada phenopacket en su fondo |
| E03 | C | Pipeline en modo VCF sobre los 30 casos |
| E04 | A·m | Recall top-1/5/10, tiempo y costo por caso |
| E05 | A·m | Caso "no diagnosticado" simulado: gen sin asociación OMIM |
| E06 | A·m | Sección de escalabilidad del reporte |

## Track 2 — Literatura (30, todo T2)

**Crow ×16 (A·l):** función de BUB1B en el SAC · CEP57 en el centrosoma · TRIP13 e inactivación
de MAD2 · genotipo-fenotipo en MVA · espectro tumoral en MVA · modelos de deficiencia de BUB1B ·
estrés proteotóxico y metabólico de la aneuploidía · estrés oxidativo en aneuploidía ·
senescencia por CIN · cGAS-STING y micronúcleos · letalidad sintética con CIN · autofagia y
desbalance estequiométrico · fármacos aprobados que modulan el checkpoint · seguridad pediátrica
de antimitóticos · firmas de aneuploidía en LINCS · HSF1/chaperonas como buffer de aneuploidía.

**Falcon ×8 (A·m):** revisión profunda de las primeras 8 consultas.

**Owl ×6 (A·l):** ¿reposicionamiento previo en MVA? · ¿rescate de BUB1B con moléculas pequeñas? ·
¿ensayos en síndromes de aneuploidía? · ¿LINCS en enfermedades de CIN? · ¿letalidad sintética en
enfermedad germinal no oncológica? · ¿MVA en organoides/iPSC?

## Track 2 — Mecanismo (8)

| ID | Tipo | Tarea |
|---|---|---|
| Y01 | A·h | Cadena causal completa: función → efecto de la variante → SAC → missegregación → aneuploidía en mosaico → CIN → fenotipo del desarrollo y predisposición tumoral. Una cita por flecha |
| Y02 | A·h | Efecto estructural de la variante (AlphaFold/PDB local; hacia fuera solo "missense en dominio X") |
| Y03 | A·m | Si hay splicing: isoforma predicha + propuesta de ensayo de minigen |
| Y04 | A·m | Red STRING/Reactome + nodos druggables |
| Y05 | A·m | Biomarcadores medibles (micronúcleos, PCS, índice mitótico, SA-β-gal) |
| Y06 | A·m | Evidencia en contra del mecanismo (obligatoria) |
| Y07 | A·m | Inventario de modelos celulares/animales disponibles para testear |
| Y08 | A·m | Diagrama del mecanismo |

## Track 2 — Reposicionamiento (12)

| ID | Tipo | Tarea |
|---|---|---|
| D01 | A·m | Robin, corrida por defecto (input: nombre de la enfermedad) |
| D02 | A·h | Robin con `prompts.py` adaptado al mecanismo (solo T2) |
| D03 | C | Open Targets: gen + red |
| D04 | C | DGIdb |
| D05 | C | ChEMBL |
| D06 | C | DrugBank si hay licencia; si no, DrugCentral |
| D07 | C | LINCS L1000: firma de knockdown del gen → reversores aprobados |
| D08 | C | DepMap: dependencias en líneas con alta CIN |
| D09 | C | TCGA/GDC: firmas de aneuploidía como contexto (no como cohorte comparable) |
| D10 | A·m | Integrar → 20 candidatos sin duplicados |
| D11 | A·h | Filtro de realismo pediátrico; ventana terapéutica explícita para cualquier citotóxico |
| D12 | A·m | Top 10 → instanciar los slots Z |

## Slots dinámicos (40)

Z01–Z10 A·m evidence card por candidato (a favor, en contra, estado regulatorio, seguridad
pediátrica, experimento falsable) · Z11–Z20 A·l Owl de precedentes por candidato ·
Z21–Z25 A·l Phoenix: química/ADMET del top 5 · Z26–Z30 A·m experimento in vitro detallado del
top 5 (línea, readout, controles, criterio de falsación) · Z31–Z40 A·m deep-dive de las 10
variantes top del Track 1 (lecturas, literatura del gen, ACMG final).

## Revisión adversarial (40, A·m, Opus)

W01–W40: una de cada cinco salidas A, muestreo estratificado por rol; **100 %** de P10–P14, Z01–Z10
y reportes finales. Busca citas inventadas (re-resuelve cada PMID), `external_calls` que violen
§3, saltos lógicos y lenguaje de consejo clínico.

## Convergencia (12, A·h)

C01–C12, una cada 2 h: consolidar, deduplicar, detectar contradicciones, repriorizar la cola,
reescribir `STATUS.md` (candidatos, confianza, bloqueos, submissions y costo restantes).

## Reportes y entrega (12)

O01 Rigor · O02 Impacto · O03 Innovación · O04 Escalabilidad (desde E06) · O05 resumen ejecutivo +
disclaimer no clínico · O06 guion del video de 3 min · O07 README + CC BY 4.0 · O08 checklist de
privacidad (`git log --all --stat`, grep de coordenadas) · O09 STATUS final · O12 pasada
anti-consejo-médico (todas A) · **O10 H aprobar borrador T2** · **O11 H aprobar push público**.

---

## §3. Presupuesto estimado

### LLM (tareas A)

Perfil por tarea, estimado para un `claude -p` con herramientas en Claude Opus 5 ($5 / $25 por MTok;
lecturas de caché ≈ 0.1× de input, escritura ≈ 1.25×):

| Peso | Tokens procesados | Output | Costo/tarea |
|---|---|---|---|
| l | ~150k (80 % caché) | ~6k | ~$0.40 |
| m | ~500k (80 % caché) | ~15k | ~$1.20 |
| h | ~1.5M (80 % caché) | ~40k | ~$3.50 |

| Fase | Tareas A | Tokens (aprox.) | Costo Opus 5 |
|---|---|---|---|
| 0 Gates/infra | 11 | 7.6M | $18 |
| 1–5 QC/align/call/mosaico/CNV (interpretación) | 12 | 6.2M | $14 |
| 6–7b Priorización + ACMG + barridos | 45 | 16M | $40 |
| 8 Sintético/calibración + escalabilidad | 13 | 11M | $23 |
| T2 literatura (+ créditos Edison aparte) | 30 | 7.4M | $18 |
| T2 mecanismo + reposicionamiento | 13 | 11M | $23 |
| Slots dinámicos | 40 | 15M | $36 |
| Adversarial | 40 | 21M | $48 |
| Convergencia | 12 | 18M | $42 |
| Reportes | 10 | 10M | $24 |
| **Total** | **226** | **~123M (~100M de caché)** | **~$285** |
| + 30 % de reintentos y sobrecosto | | | **~$370** |

Si usas Sonnet 5 ($2 / $10) en las tareas `l` y `m` y dejas en Opus solo ACMG, adversarial,
convergencia, mecanismo y reportes, baja a unos **$200**. La decisión es tuya.

Ojo: si Claude Code corre con suscripción y no con API key, el límite es de rate y no de dólares.
15 workers concurrentes con suscripción se van a trabar.

### Cómputo AWS (tareas C)

| Bloque | Horas-instancia (64 vCPU) | Costo on-demand aprox. |
|---|---|---|
| Alineamiento + DV + HC | ~13 h | ~$40 |
| Mosaico (4 callers) + SV/CNV + MoChA | ~15 h | ~$45 |
| Anotación + priorización | ~4 h | ~$12 |
| Set sintético (HG002, pipeline completo) | ~25 h | ~$70 |
| Escalabilidad (30 casos en modo VCF) | ~3 h | ~$9 |
| Almacenamiento (EBS 1.5 TB + S3 de recursos, unas 2 semanas) | — | ~$70 |
| **Total** | **~60 h** | **~$250** (spot: ~$120) |

Las tarifas horarias son órdenes de magnitud y se miden de verdad en L01.

### FutureHouse / Edison

Se compran créditos y el precio por consulta no está en la documentación que revisé. G16/G17
miden el costo real de 1 Crow y 1 Robin, y con eso se fija un tope en `state.json` antes de lanzar
las 30 consultas F y los 10 Owl Z. Robin usa `o4-mini` (OpenAI) por defecto vía LiteLLM: su input
es solo el nombre de la enfermedad (T2), pero hay que decidir si lo apuntamos a Claude.

### Tiempo

El camino crítico (FASTQ → CRAM → callers → anotación → ranking) son unas 15–20 h de pared en una
sola instancia grande. Con el fast-path sobre el VCF entregado (G21) hay candidatos en la hora 2.
La ventana cierra el 2026-10-24: conviene correr **sprints** de 24 h, no una sola corrida de 24 h.

---

## §4. Los tres supuestos más frágiles y cómo los pruebo en la hora 1

### 1. "Podemos procesar los datos del paciente con Claude y en AWS, desde esta máquina"

- **Por qué es frágil:** la regla 2 del prompt prohíbe mandar coordenadas a APIs de terceros, y la
  API de Anthropic es un tercero. Si un worker lee un VCF, las coordenadas salen a Anthropic. Además
  el repo está en `OneDrive\Desktop`: `data/` se sincronizaría a Microsoft apenas se descargue.
  La ficha del dataset en HF solo pide aceptar las reglas del hackathon (licencia CC-BY-4.0) y no
  dice nada sobre LLMs ni nube.
- **Prueba:** G01 (leer las reglas; si no son claras, preguntar a los organizadores en Discussions
  sin dar detalles del caso) + G04 (sacar el repo de OneDrive antes de descargar un solo byte) +
  G05/G06 (egress-guard con 30/30 canarios bloqueados).
- **Plan B:** si el DUA no cubre LLMs, los agentes solo ven derivados agregados T2/T3 y el ranking
  lo producen scripts deterministas escritos por agentes a partir de datos sintéticos.

### 2. "Tenemos cómputo para reprocesar 85 GB de FASTQ con 4+ callers en la ventana"

- **Por qué es frágil:** esta máquina tiene 4 hilos, 7.7 GB de RAM, 105 GB libres y ninguna distro
  WSL. El índice de bwa-mem2 para GRCh38 no cabe en esa RAM; 85 GB de input más CRAM e intermedios
  no caben en ese disco; DeepVariant WGS con 4 hilos tarda días. `run.sh` además asume bash GNU y
  tmux.
- **Prueba:** G03 + L01 (1M pares en la instancia → extrapolar horas y dólares). En paralelo, G09
  confirma si entre los 11 archivos viene un CRAM/BAM (evita realinear) y G21 arranca el fast-path.
- **Criterio de corte:** si la extrapolación da más de 12 h de alineamiento, se usa el VCF entregado
  como fuente principal para germinal, y el cómputo se reserva para mosaico y aneuploidía.

### 3. "La señal causal está en los VAF intermedios (5–35 %)"

- **Por qué es frágil:** en MVA1–3 las variantes causales publicadas son constitucionales y
  bialélicas (VAF ~50/100 %). Lo mosaico es la aneuploidía, no la variante. Si tratamos los VAF
  intermedios como señal causal, la lista se llena de artefactos y F-max cae, que es justo lo que
  la métrica castiga. Riesgo adicional: hay reportes de MVA1 con un alelo truncante más un alelo
  hipomórfico no codificante o de baja expresión en BUB1B `[UNVERIFIED — resolver en G18]`; un
  filtro solo codificante no lo vería.
- **Prueba:** G22 (balance alélico por cromosoma desde el VCF entregado: ¿hay aneuploidía visible en
  este tejido?) + G21 (candidatos recesivos germinales en SAC) + S01–S03 con evidencia negativa
  explícita. Si aparece un compuesto germinal plausible en BUB1B/CEP57/TRIP13 en la hora 2, el
  canal de mosaico queda como evidencia ortogonal y no como fuente de candidatos.

---

## Defectos del scaffold actual que se corrigen con G19 (después del visto bueno)

1. Todas las tareas de `queue.jsonl` tienen `depends_on: []`: el dispatcher lanzaría Exomiser antes
   de la anotación y DeepVariant antes del alineamiento.
2. `dispatch.py` espera a que termine todo el lote: una tarea de 30 min bloquea 14 slots libres.
3. `--permission-mode acceptEdits` en headless deniega Bash: los workers no pueden correr nada.
   Hace falta `--allowedTools` explícito por rol.
4. Jobs de horas (alinear, DeepVariant) están modelados como tareas de LLM con timeout de 30 min.
5. La fuga de datos es autodeclarada (`leaked_patient_data`); no hay un control real de egreso.
6. 24 tareas de agente, una por cromosoma, para algo que es un solo script vectorizado.
7. "RZZ" en la lista de genes no es un gen.
8. El backoff se dispara ante cualquier fallo, no solo ante 429.
9. No hay tracker de costo ni tope de presupuesto.

---

# Actualización tras leer las reglas (2026-09-12)

Detalle completo en `RULES-FINDINGS.md`. Deltas sobre el plan de arriba:

## Gates

- **G01 resuelto.** Los LLM de terceros están permitidos bajo la prueba de procesador. Se reemplaza
  por tres tareas: **G01a** fijar proveedor/plan (API comercial recomendada), **G01b** verificar que
  el entrenamiento está desactivado y prohibir calificar outputs en toda la corrida, **G01c**
  redactar la línea obligatoria de divulgación para el methods writeup.
- **Nuevo grupo DEL (5).** DEL01 layout con todo lo purgable bajo un root único, incluidos los
  transcripts locales de Claude Code · DEL02 `make purge` + verificación · DEL03 auditoría de que
  ningún artefacto conservado permite reconstruir el genoma · DEL04 ensayo de purga en seco ·
  DEL05 (H) atestación por correo antes del 2026-11-23.

## Track 1

- **G11/G12 se fusionan en una sola tarea de validación**: el fenotipo ya viene como términos HPO
  estandarizados, no hay que extraerlos de texto libre.
- **L02–L06 quedan condicionadas.** Las reglas dicen que los datos crudos vienen opcionalmente en
  BAM/CRAM. Si hay CRAM, no se realinea: se ahorran unas 13 horas de instancia y unos 40 dólares.
  G08 lo decide.
- **G13 cambia de fasta.** El VCF viene de un analysis set con decoys hs38d1 y contigs sin prefijo
  `chr`. Si realineamos, usamos ese fasta exacto.
- **Nueva tarea de nomenclatura de contigs** en el formatter, con test: el VCF usa `15`, el CSV de
  submission usa `chr15`. Es el error tonto que cuesta una submission entera.
- **Nueva tarea de cis/trans** con PGT/PID del VCF más phasing por lecturas. Es el muro donde está
  hoy la comunidad (discusión #22) y es material directo para el writeup.
- **Mosaico se recorta**: Mutect2 tumor-only y MosaicHunter, más MoChA para la aneuploidía.
  DeepMosaic y MosaicForecast salen del camino crítico y quedan como verificación opcional.
- **X10 cambia de objetivo.** El techo son 10 filas y F-max toma el máximo sobre umbrales: llenar
  las 10 filas es seguro si el orden es correcto. La tarea pasa a calibrar el *orden* y la
  separación de los EPCR, no el tamaño de la lista.
- **Presupuesto de submissions: 5, no 6** (discrepancia entre `config.py` y la pestaña de envío).

## Reequilibrio

El leaderboard de Track 1 ya tiene puntajes perfectos y los organizadores lo llaman un track
fundacional. El peso se mueve del cómputo exótico hacia el methods writeup y Track 2, que es donde
se reparte el puntaje del panel.

**Total actualizado: ~336 tareas.** Cómputo AWS estimado baja a unos 140-250 dólares según haya o no
CRAM. El presupuesto de LLM no cambia: todo en Opus 5, unos 370 dólares con reintentos.

## Supuestos frágiles, revisados

1. ~~Permiso de datos~~ **resuelto**. El nuevo riesgo es de *configuración*: si corremos con una
   suscripción de consumidor en vez de API comercial, caemos en la pregunta abierta de la discusión
   #21, que nadie ha respondido. Prueba: G01a/G01b antes de tocar un solo dato.
2. **Cómputo**: sin cambios, pero el riesgo baja mucho si hay CRAM o si trabajamos sobre el VCF.
3. **VAF intermedios**: sin cambios, y se refuerza. Que el track sea "alcanzable" y ya tenga
   puntajes perfectos apunta a una variante germinal recuperable del VCF entregado.
