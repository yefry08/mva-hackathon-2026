# MVA Hackathon 2026 — predicción de variante y reposicionamiento

Participación en *Rare Disease, Real Kid: The MVA Hackathon 2026* (Sage
Bionetworks, MVA Society, Hugging Face, BEACON, con AWS y Anthropic).

## Qué hay aquí

Un pipeline de priorización de variantes que corre entero con Python de la
biblioteca estándar, en un portátil, sin nube ni contenedores, y sin que los
datos del paciente salgan de la máquina.

**Track 1.** El candidato es un heterocigoto compuesto en **BUB1B** (MVA1):
p.Leu737Ter, nonsense con NMD predicho y clasificada Pathogenic/Likely_pathogenic
en ClinVar para *Mosaic variegated aneuploidy syndrome 1*, más p.Asn1002Lys,
missense presente en 1 de 1.461.878 alelos de gnomAD.

**Track 2.** Caracterización del mecanismo con 16 citas verificadas y tres
candidatos de reposicionamiento planteados **como hipótesis falsables**, con su
experimento de refutación, su contraevidencia y su filtro de realismo pediátrico.

**Nada de esto es consejo médico.**

![Embudo de siete etapas del genoma completo a un gen, sin panel y sin ClinVar: de 5.012.204 variantes a BUB1B](reports/figures/fig1_embudo_ciego.svg)

## Documentos

| Archivo | Qué contiene |
|---|---|
| `reports/TRACK1_METHODS.md` | Methods writeup: pipeline, decisiones, límites |
| `reports/ACMG_BUB1B.md` | Clasificación ACMG criterio por criterio, con los no aplicados |
| `reports/MECHANISM_BUB1B.md` | Cadena mecanística, ocho eslabones, y dónde falla |
| `reports/TRACK2_REPORT.md` | Reporte contra la rúbrica del panel |
| `reports/REPURPOSING.md` | Fichas de evidencia de los candidatos |
| `reports/SCALABILITY.md` | El mismo pipeline sobre tres genomas sanos, como controles negativos |
| `reports/drug_landscape.tsv` | Esencialidad (DepMap) y fármacos conocidos para los 32 genes del panel |
| `reports/cmap_signature.tsv` | Firma transcripcional de la pérdida de BUB1B (LINCS), contra 5.211 knockouts |
| `reports/cmap_reversal.tsv`, `cmap_summary.json` | Búsqueda por conectividad, sus controles y por qué ningún candidato sobrevive |
| `reports/pitch/index.html` | Presentación para grabar el video de 3 minutos |
| `reports/STATUS.md` | Estado y bitácora |

## Reproducir

Requisitos: Python 3.12 o superior, biblioteca estándar. Nada más.

```bash
python scripts/gene_sweep.py data/<tu>.vcf.gz   # barrido por panel de genes
python scripts/annotate_gnomad.py               # gnomAD por rangos de bytes remotos
python scripts/candidates.py                    # consecuencia y fase
python scripts/clinvar_annot.py                 # cruce dirigido con ClinVar
python scripts/clinvar_genomewide.py            # ranking ciego, sin prior de genes
python scripts/protein_consequence.py           # traducción local
python scripts/make_submission.py in.tsv out.csv --proband PROBAND01
```

Recursos de referencia que se descargan solos o con `curl` (MANE, ClinVar,
chr15). gnomAD **no** se descarga: `scripts/tabix_remote.py` implementa un cliente
tabix por HTTP que pide solo los bloques que cubren cada gen, unos pocos MB en vez
de 7.4 GB por cromosoma.

Correrlo sobre otro genoma:

```bash
export MVA_WORK=work_otro
python scripts/gene_sweep.py ruta/otro.vcf.gz
```

Pruebas:

```bash
python tests/test_submission.py      # 12 casos, incluida la réplica del scoring
python tests/test_egress_guard.py    # 30 canarios, 10 consultas legítimas
python scripts/privacy_check.py      # checklist bloqueante antes de publicar
```

## Qué NO hay aquí, y por qué

Ningún dato del paciente. Ni VCF, ni BAM, ni FASTQ, ni tablas de genotipos, ni
identificadores de muestra, ni siquiera en el historial de git. Las reglas del
hackathon prohíben redistribuir los datos por cualquier canal y obligan a
borrarlos al cierre.

Eso se comprueba con `scripts/privacy_check.py`, que es bloqueante y revisa el
índice, **todo el historial de commits**, el `.gitignore` y el contenido de cada
archivo rastreado. La única excepción declarada es `submissions/`, que contiene la
lista rankeada de candidatos: ese es el entregable, y las reglas permiten
conservarlo.

Durante el análisis, un hook `PreToolUse` (`.claude/hooks/egress_guard.py`)
bloqueó automáticamente cualquier llamada a un servicio externo que llevara
coordenadas, rsIDs, HGVS con posición o identificadores de muestra.

## Uso de IA

El análisis lo ejecutaron scripts deterministas. Un modelo de lenguaje (Claude,
Anthropic) escribió esos scripts, interpretó salidas agregadas y redactó los
informes. El modelo recibió símbolos de genes, consecuencias a nivel de proteína,
métricas de calidad, frecuencias y clasificaciones; no recibió coordenadas
genómicas, tablas de genotipos ni identificadores de muestra. Todas las citas se
resolvieron contra PubMed durante la sesión.

## Licencia

CC BY 4.0. Ver `LICENSE`.

## Agradecimiento

> *This work was made possible through the Hackathon, organized by Sage
> Bionetworks in partnership with the MVA Society, Hugging Face, and BEACON (The
> Benchmarking, Evaluation, and Assessment Consortium for Science), with prize
> sponsorship from AWS and Anthropic. We are deeply grateful to the child and
> their family who generously contributed their data and their story to advance
> research into this rare disease. We acknowledge their trust in making this
> Hackathon possible.*
