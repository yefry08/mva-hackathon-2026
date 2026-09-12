# CLAUDE.md — MVA Hackathon 2026

Todo subagente hereda este archivo. Léelo completo antes de actuar. Si algo aquí contradice la
tarea que te asignaron, gana este archivo y reportas `status: "blocked"`.

## 1. Qué estamos haciendo

*Rare Disease, Real Kid: The MVA Hackathon 2026* (Sage Bionetworks / MVA Society / BEACON /
Hugging Face / AWS / Anthropic). Ventana de envíos: 2026-08-24 → 2026-10-24.

- **Track 1 — Predicción de variantes.** Lista rankeada, **máximo 10 filas**, repo reproducible y
  methods writeup. Scoring: rank points (100/50/25/10 por rank) + F-max. `config.py` permite 6
  submissions pero la pestaña de envío dice 5: **planeamos con 5** (1 sanity, 3 de iteración, 1 de
  reserva). Ya hay puntajes perfectos en el leaderboard: **el diferenciador es el writeup**, no
  recuperar la variante.
- **Track 2 — Reposicionamiento.** Mecanismo + fármacos ya aprobados como **hipótesis falsables
  para seguimiento experimental**. Panel: Rigor 35 · Impacto 25 · Innovación 25 · Escalabilidad 15.
  Reporte + repo + video de 3 min. **3 submissions**: 1 borrador revisado internamente, 2 finales.

## 2. Qué NO estamos haciendo

Nada de consejo médico, diagnóstico ni plan de tratamiento. Frases como "el paciente debería
tomar", "recomendamos", "tratamiento indicado" o "dosis" son un error: se marcan y se reescriben.
Un candidato de Track 2 es siempre "hipótesis X, que el experimento Y falsaría".

## 3. Niveles de datos (qué puede ver quién)

| Nivel | Qué es | Dónde vive | Quién lo toca |
|---|---|---|---|
| T0 crudo | FASTQ, BAM/CRAM, VCF, metadata clínica | `data/` | Solo jobs de cómputo deterministas (Snakemake). **Ningún agente hace `cat`/`Read` de un archivo T0.** |
| T1 derivado del paciente | Tablas de candidatos, VAF, coordenadas, HPO del caso | `work/` | Jobs de cómputo + agentes Claude bajo **términos de procesador** (ver §3.1). Nunca sale a otro servicio. |
| T2 nivel gen/mecanismo | Símbolo de gen, vía, dominio, clase de variante sin posición | cualquiera | Puede ir a FutureHouse/Edison, PubMed, Open Targets, etc. |
| T3 público | Reportes, figuras agregadas, submissions | `reports/`, `track2/`, `submissions/` | Sin coordenadas salvo la lista de variantes entregada al leaderboard. |

### 3.1 La prueba de procesador (regla oficial, discusión #2 del Space)

Un servicio al que le mandas datos es una **herramienta** si procesa solo para devolverte un
resultado, no toma derechos y no puede usar el contenido para fines propios. Es un **receptor** si
puede usar, retener o aprender de los datos, y mandarle datos ahí es una divulgación prohibida.
Dos condiciones simultáneas: (a) sin entrenamiento sobre inputs ni outputs, sin derechos tomados;
(b) retención limitada en tiempo y propósito.

- Claude bajo **términos comerciales de API** califica como herramienta: los workers pueden leer T1.
- **Nunca califiques outputs** (pulgares, feedback): el contenido calificado puede quedar excluido
  de la exención de entrenamiento.
- Ningún otro servicio tiene términos verificados todavía → todo lo demás se queda en T2.
- En el methods writeup va una línea con proveedor, plan y ajuste usado. Es obligatoria.

### 3.2 Borrado y atestación

Todo lo que lleve el genoma del niño se borra dentro de los 30 días del cierre (cierre 2026-10-24).
Se borran `data/`, `work/`, índices, cachés, tablas de genotipos a escala genómica **y los
transcripts locales de Claude Code que contengan bloques de variantes**. Se conservan la lista
entregada, los HPO, los rankings de genes, el mecanismo, los candidatos de fármaco, el código y el
reporte. Corolario de diseño: todo lo borrable vive bajo un único root purgable y existe
`make purge` desde el día uno. La atestación por correo la manda el humano.

Reglas duras:
1. `data/`, `ref/` y `work/` están en `.gitignore`. Antes de cada `git add`, `git status` y el hook
   `guard` de `run.sh`. Si aparece un archivo de paciente en staging, aborta.
2. **Nunca** envíes a una API de terceros: FASTQ/BAM/VCF, coordenadas (`chr:pos`), rsIDs del caso,
   HGVS con posición, IDs de muestra, nombres de archivo originales, **ni la combinación de un gen
   con el fenotipo específico del paciente**. Las consultas externas son del tipo "función de
   TRIP13 en la inactivación de MAD2", no "paciente con microcefalia y variante en TRIP13".
3. El hook `PreToolUse` `egress-guard` bloquea MCP/WebFetch/WebSearch/`curl` con esos patrones.
   Si te bloquea, **no lo rodees**: reformula a nivel gen o reporta `blocked`.
4. AMELIE es un servicio web: queda fuera salvo aprobación explícita del humano.
5. El repo **no** vive en OneDrive ni en otra carpeta sincronizada a la nube.

## 4. Reglas de evidencia

- Toda afirmación biológica lleva DOI o PMID **resuelto** con el MCP de PubMed/Europe PMC en esta
  misma tarea. Si no resuelve → `[UNVERIFIED]` y no entra al reporte final.
- Prohibido citar de memoria. Un paper que no abriste no se cita.
- Cada hipótesis lleva `evidence_against`. Vacío solo si escribes dónde buscaste.
- Extrapolar de línea celular o ratón a un niño es un salto: nómbralo cuando lo hagas.
- Usa la skill `evidence-card` para todo claim que vaya a un reporte.

## 5. Contexto biológico (prior documentado, no atajo)

- MVA es **autosómica recesiva**. Genes reportados: **BUB1B** (MVA1), **CEP57** (MVA2),
  **TRIP13** (MVA3). La red del spindle assembly checkpoint (SAC) es el prior extendido.
- **Hipótesis principal: variantes causales constitucionales bialélicas** (homocigotas o
  heterocigotas compuestas, VAF ≈ 50 % / 100 %). Lo que es mosaico por definición en MVA es la
  **aneuploidía**, no necesariamente la variante causal.
- El llamado en mosaico es un **canal secundario** con tres usos: (a) segundo golpe somático o
  reversión, (b) evidencia ortogonal de aneuploidía (BAF/cobertura, MoChA), (c) control de
  sesgo. Un candidato con VAF 5–35 % entra al ranking solo con ≥2 callers, inspección IGV y un
  poder de detección documentado para ese VAF y esa cobertura.
- La aneuploidía distorsiona el balance alélico de variantes germinales en los cromosomas
  afectados. No filtres hets por allele balance estricto (usa 0.2–0.8 y marca el cromosoma).
- No descartes non-coding: un alelo hipomórfico de splicing, UTR o regulatorio puede ser la
  segunda copia de un compuesto.
- Siempre dos rankings: `ranking_ciego.tsv` (sin prior) y `ranking_prior.tsv`, con el delta
  explicado. El pipeline tiene que poder recuperar un gen nuevo.
- Si hay un par heterocigoto compuesto, **las dos variantes** van a la lista (el scoring da medio
  crédito por una sola).

## 6. Tipos de tarea

- `C` (cómputo): regla de Snakemake, sin LLM, determinista, con versión de contenedor pinneada.
- `A` (agente): `claude -p` headless con rol y contrato JSON. Escribe código, interpreta
  resultados y redacta. No ejecuta jobs de más de 20 min: los encola como `C`.
- `H` (humano): gate. El dispatcher se detiene y espera.

## 7. Gates humanos (no negociables)

- Ninguna submission al leaderboard sin confirmación explícita del humano en el chat.
- Ningún push a un repo público sin confirmación explícita.
- Ninguna llamada a un servicio externo nuevo (no listado en `.mcp.json`) sin aprobación.
- Presupuesto: si el costo acumulado de LLM supera el tope de `state.json`, el dispatcher pausa.

## 8. Contrato de output de worker

Un archivo por tarea en `reports/workers/<task_id>.json`:

```json
{
  "task_id": "P10",
  "role": "acmg-classifier",
  "kind": "A",
  "status": "done|failed|blocked",
  "max_data_tier_read": "T1",
  "artifacts": ["work/acmg/batch1.tsv"],
  "findings": [
    {"claim": "...", "evidence": [{"id": "PMID:...", "verified": true}], "confidence": 0.0}
  ],
  "evidence_against": ["..."],
  "contradictions": ["..."],
  "next_tasks": [{"role": "...", "kind": "A|C", "description": "...", "depends_on": ["..."]}],
  "external_calls": [{"service": "futurehouse:crow", "query": "texto exacto enviado"}],
  "usage": {"model": "...", "input_tokens": 0, "output_tokens": 0, "cost_usd": 0.0}
}
```

`external_calls` registra el texto exacto que salió. El `adversarial-reviewer` lo audita contra
§3. `leaked_patient_data` ya no es autodeclarado: lo decide el reviewer.

## 9. Herramientas y licencias

- Referencia: **GRCh38**, y el VCF entregado viene de
  `GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_hs38d1_maskedGRC_exclusions_v2_no_chr.fasta`
  (contigs **sin prefijo `chr`**, con decoys hs38d1). Si realineas, usa ese fasta, no el analysis
  set genérico. El CSV de submission usa `chr15`: la conversión de nomenclatura es un paso con test.
  T2T-CHM13 solo como referencia secundaria para regiones difíciles.
- El VCF trae `PGT`/`PID`: hay fase parcial ya disponible. Úsala antes de inventar phasing.
- OMIM (genemap2) y DrugBank requieren licencia. Sin licencia se documenta y se usan
  alternativas (Orphanet/HPO; DrugCentral/ChEMBL).
- Robin recibe solo el nombre de la enfermedad (T2 por diseño). Su LLM se configura por LiteLLM.
- NHANES, WHO GHO, UK Biobank, ADNI, SEER, PhysioNet y OpenNeuro no aportan a este caso: fuera.
