#!/usr/bin/env python3
"""Genera la cola de tareas atomicas para la corrida multiagente.

Uso:  python orchestrator/generate_queue.py > orchestrator/queue.jsonl
"""
import json
import sys

SAC_GENES = [
    "BUB1B", "CEP57", "TRIP13", "BUB1", "BUB3", "MAD1L1", "MAD2L1", "TTK",
    "CENPE", "AURKB", "AURKA", "PLK1", "PLK4", "CDC20", "ESPL1", "PTTG1",
    "NDC80", "CENPF", "KNL1", "ZW10", "ZWILCH", "RZZ", "SKA1", "KIF11",
    "CENPA", "INCENP", "SGO1", "CDCA8", "MASTL", "CHMP4C",
]

CHROMS = [str(c) for c in range(1, 23)] + ["X", "Y"]

tasks = []
n = 0


def add(role, phase, desc, blocking=None, priority=5):
    global n
    n += 1
    tasks.append({
        "task_id": f"T{n:03d}",
        "role": role,
        "phase": phase,
        "priority": priority,
        "description": desc,
        "depends_on": blocking or [],
        "status": "pending",
    })


# ---------------------------------------------------------------- Fase 0: setup
add("qc", 0, "Inventariar data/: listar todos los archivos, tamanos, checksums, y determinar tipo (FASTQ R1/R2, VCF, metadata clinica).", priority=0)
add("qc", 0, "Determinar el build de referencia del VCF entregado (GRCh37 vs GRCh38) desde el header. NO alinear nada hasta confirmarlo.", priority=0)
add("qc", 0, "Parsear la descripcion clinica y convertirla a una lista de terminos HPO con IDs. Guardar en pipelines/07_prioritize/hpo_terms.txt", priority=0)
add("qc", 0, "Verificar .gitignore: confirmar que data/ y ref/ no aparecen en git status. Abortar la corrida si aparecen.", priority=0)
add("qc", 0, "Descargar y indexar el analysis set GRCh38 no-alt (BWA, samtools faidx, GATK dict).", priority=1)

# ---------------------------------------------------------------- Fase 1: QC
add("qc", 1, "FastQC sobre todos los FASTQ + MultiQC consolidado.", priority=1)
add("qc", 1, "Estimar cobertura media y uniformidad; reportar regiones con cobertura <10x.", priority=2)
add("qc", 1, "VerifyBamID: estimar contaminacion cruzada de muestra.", priority=2)
add("qc", 1, "Determinar sexo genetico desde cobertura X/Y y contrastar con metadata clinica.", priority=2)
add("qc", 1, "Estimar tasa de duplicados y complejidad de libreria.", priority=3)

# ---------------------------------------------------------------- Fase 2: alineamiento
add("aligner", 2, "Alinear FASTQ con BWA-MEM2 contra el no-alt analysis set. Read groups correctos.", priority=1)
add("aligner", 2, "Alinear en paralelo con DRAGMAP como control ortogonal.", priority=3)
add("aligner", 2, "MarkDuplicates + BQSR sobre el BAM de BWA-MEM2.", priority=1)
add("aligner", 2, "Comparar metricas de alineamiento BWA-MEM2 vs DRAGMAP; documentar discrepancias.", priority=4)
add("aligner", 2, "Generar BAM de cobertura por ventana de 10kb para uso en aneuploidia.", priority=2)

# ---------------------------------------------------------------- Fase 3: llamado germinal
add("germline-caller", 3, "DeepVariant sobre el BAM final. Generar VCF + gVCF.", priority=1)
add("germline-caller", 3, "GATK HaplotypeCaller sobre el BAM final.", priority=1)
add("germline-caller", 3, "Consenso DeepVariant/HaplotypeCaller: variantes en ambos, y set de discordantes para revision manual.", priority=2)
add("germline-caller", 3, "Comparar nuestro VCF contra el VCF entregado por el hackathon; cuantificar concordancia y explicar diferencias.", priority=2)
add("germline-caller", 3, "Hard-filter y VQSR; reportar cuantas variantes sobreviven cada paso.", priority=3)

# ---------------------------------------------------------------- Fase 4: mosaico (el track lo decide aqui)
add("mosaic-caller", 4, "MuTect2 en modo tumor-only sobre el BAM, con panel of normals publico.", priority=1)
add("mosaic-caller", 4, "MosaicHunter con el modelo single-sample.", priority=1)
add("mosaic-caller", 4, "DeepMosaic sobre el BAM final.", priority=2)
add("mosaic-caller", 4, "Calcular VAF por variante y construir el histograma de VAF. Identificar el pico de mosaicismo.", priority=1)
add("mosaic-caller", 4, "Extraer el set de variantes con VAF 5-35% en genes del SAC. Este es el set de mayor interes.", priority=1)
add("mosaic-caller", 4, "Validar visualmente (IGV screenshots) las 20 variantes en mosaico mejor rankeadas.", priority=3)
add("mosaic-caller", 4, "Estimar la tasa de falsos positivos del llamado en mosaico usando regiones control.", priority=3)
add("mosaic-caller", 4, "Consenso entre los tres callers de mosaico; reportar union e interseccion.", priority=2)

# ---------------------------------------------------------------- Fase 5: CNV / SV / aneuploidia
add("sv-cnv", 5, "Manta para SVs estructurales.", priority=2)
add("sv-cnv", 5, "GRIDSS como caller ortogonal de SV.", priority=3)
add("sv-cnv", 5, "CNVnator y cn.mops para CNVs.", priority=2)
add("sv-cnv", 5, "Anotar SVs contra gnomAD-SV y DGV para filtrar polimorfismos comunes.", priority=3)
for c in CHROMS:
    add("aneuploidy-profiler", 5, f"Perfilar cobertura normalizada y B-allele frequency en chr{c}; cuantificar evidencia de aneuploidia en mosaico y estimar fraccion celular.", priority=2)
add("aneuploidy-profiler", 5, "Consolidar el perfil de aneuploidia genome-wide en una figura y una tabla. Esta es evidencia ortogonal fuerte del fenotipo variegado.", priority=1)
add("aneuploidy-profiler", 5, "Contrastar el patron de aneuploidia observado contra los patrones reportados en literatura de MVA.", priority=2)

# ---------------------------------------------------------------- Fase 6: anotacion
add("annotator", 6, "VEP con plugin LOFTEE sobre el VCF de consenso.", priority=1)
add("annotator", 6, "Anotar con dbNSFP (CADD, REVEL, MetaSVM, MutationTaster).", priority=1)
add("annotator", 6, "SpliceAI sobre todas las variantes intronicas y sinonimas cercanas a sitios de splicing.", priority=1)
add("annotator", 6, "AlphaMissense sobre todas las missense.", priority=1)
add("annotator", 6, "Anotar frecuencias gnomAD v4 (global y por poblacion).", priority=1)
add("annotator", 6, "Cruzar contra ClinVar; extraer clasificaciones y conflictos de interpretacion.", priority=1)
add("annotator", 6, "Anotar constraint scores (pLI, LOEUF, o/e) por gen.", priority=2)
add("annotator", 6, "Anotar expresion GTEx en tejidos relevantes (cerebro, medula osea, fibroblasto).", priority=3)
add("annotator", 6, "Anotar variantes en regiones regulatorias y UTR; no descartar non-coding por defecto.", priority=3)

# ---------------------------------------------------------------- Fase 7: priorizacion
add("phenotype-matcher", 7, "Exomiser con los terminos HPO, modelo AR. Guardar ranking completo.", priority=1)
add("phenotype-matcher", 7, "LIRICAL con los mismos terminos HPO.", priority=1)
add("phenotype-matcher", 7, "AMELIE sobre el set filtrado.", priority=2)
add("phenotype-matcher", 7, "Phen2Gene / GADO como cuarto ranking independiente.", priority=3)
add("phenotype-matcher", 7, "Comparar los cuatro rankings; calcular concordancia y construir un ranking de consenso ponderado.", priority=1)
add("acmg-classifier", 7, "Filtrar por modelo recesivo: homocigotos y heterocigotos compuestos con AF gnomAD <0.001.", priority=1)
add("acmg-classifier", 7, "Para cada par heterocigoto compuesto candidato, verificar fase si los datos lo permiten (read-backed phasing).", priority=2)
add("acmg-classifier", 7, "Clasificar ACMG/AMP las 50 variantes top; documentar cada criterio aplicado.", priority=1)
add("acmg-classifier", 7, "Generar el ranking CIEGO (sin prior de genes MVA) y guardarlo aparte.", priority=1)
add("acmg-classifier", 7, "Generar el ranking CON PRIOR sobre genes MVA y red SAC. Documentar el delta contra el ciego.", priority=1)
add("acmg-classifier", 7, "Revisar explicitamente variantes fuera del prior que rankeen alto: candidatas a gen nuevo.", priority=2)
add("acmg-classifier", 7, "Identificar hallazgos secundarios / incidentales reportables (no penalizan el score automatico).", priority=3)

# ---------------------------------------------------------------- Fase 7b: barrido por gen
for g in SAC_GENES:
    add("annotator", 7, f"Barrido dirigido de {g}: todas las variantes (codificantes, splicing, UTR, CNV) con VAF y cobertura. Reportar aunque no pasen filtros.", priority=3)

# ---------------------------------------------------------------- Fase 8: ranking y calibracion
add("rank", 8, "Construir el score compuesto: consenso fenotipico + ACMG + VAF mosaico + constraint + evidencia de aneuploidia.", priority=1)
add("rank", 8, "Calibrar el umbral de confianza sobre un set sintetico de validacion (variantes causales conocidas insertadas en datos publicos).", priority=1)
add("rank", 8, "Optimizar F-max: barrer el tamano de la lista entregada (5, 10, 20, 50) y estimar el trade-off precision/recall.", priority=1)
add("rank", 8, "Verificar que si hay un par heterocigoto compuesto AMBAS variantes esten en la lista. Medio credito por una sola.", priority=1)
add("rank", 8, "Generar submission Track 1 #1 (sanity check). NO enviar sin confirmacion humana.", priority=2)
add("rank", 8, "Redactar el methods writeup del Track 1 (evaluado por panel aparte del score automatico).", priority=2)

# ---------------------------------------------------------------- Track 2: literatura
LIT_QUERIES = [
    "Funcion normal de BUB1B en el spindle assembly checkpoint",
    "Mecanismo molecular de CEP57 en el ensamblaje del centrosoma",
    "Rol de TRIP13 en la inactivacion de MAD2 y silenciamiento del checkpoint",
    "Correlacion genotipo-fenotipo en MVA reportada a la fecha",
    "Espectro tumoral en pacientes con MVA (Wilms, rabdomiosarcoma, leucemia)",
    "Modelos celulares y animales de deficiencia de BUB1B",
    "Consecuencias metabolicas y de estres proteotoxico de la aneuploidia",
    "Respuesta a estres oxidativo en celulas aneuploides",
    "Senescencia celular inducida por inestabilidad cromosomica",
    "Activacion de cGAS-STING por micronucleos en celulas con CIN",
    "Letalidad sintetica con inestabilidad cromosomica: dianas conocidas",
    "Autofagia y aclaramiento de proteinas en desbalance estequiometrico por aneuploidia",
    "Terapias en investigacion para sindromes de inestabilidad cromosomica",
    "Farmacos aprobados que modulan el checkpoint mitotico",
    "Seguridad pediatrica de agentes antimitoticos: ventana terapeutica",
    "Vigilancia oncologica en sindromes de predisposicion tumoral pediatrica",
    "Reposicionamiento de farmacos en enfermedades ultra-raras: casos exitosos",
    "Firmas transcripcionales de aneuploidia y conectividad LINCS L1000",
    "Dependencias DepMap en lineas celulares con alta inestabilidad cromosomica",
    "Chaperonas y HSF1 como buffer del estres por aneuploidia",
]
for q in LIT_QUERIES:
    add("literature-crow", 9, f"Crow (busqueda concisa): {q}. Solo nivel gen/mecanismo, NUNCA datos del paciente.", priority=3)

for q in LIT_QUERIES[:8]:
    add("literature-falcon", 9, f"Falcon (revision profunda): {q}. Reporte estructurado con citas verificables.", priority=3)

PRECEDENT_QUERIES = [
    "Ha intentado alguien reposicionar farmacos aprobados para MVA",
    "Ha rescatado alguien fenotipos de BUB1B con moleculas pequenas",
    "Existe algun ensayo clinico en sindromes de aneuploidia en mosaico",
    "Ha usado alguien conectividad LINCS para enfermedades de inestabilidad cromosomica",
    "Se ha aplicado letalidad sintetica a enfermedad germinal (no oncologica)",
    "Ha modelado alguien MVA en organoides o iPSC de paciente",
]
for q in PRECEDENT_QUERIES:
    add("precedent-owl", 9, f"Owl (precedentes): {q}", priority=3)

# ---------------------------------------------------------------- Track 2: mecanismo
add("mechanism-modeler", 10, "Cadena causal completa: funcion normal -> efecto de la variante -> checkpoint mitotico -> missegregacion -> aneuploidia en mosaico -> CIN -> fenotipo del desarrollo y predisposicion tumoral. Cada flecha con cita.", priority=2)
add("mechanism-modeler", 10, "Modelar estructuralmente el efecto de la variante (AlphaFold/PDB): dominio afectado, interfaz de interaccion perdida.", priority=3)
add("mechanism-modeler", 10, "Mapear la red de interaccion afectada con STRING y Reactome; identificar nodos druggables.", priority=2)
add("mechanism-modeler", 10, "Identificar biomarcadores medibles del mecanismo (micronucleos, indice mitotico, marcadores de senescencia).", priority=3)
add("mechanism-modeler", 10, "Contrastar el mecanismo propuesto contra la evidencia que lo contradice. Seccion obligatoria.", priority=2)

# ---------------------------------------------------------------- Track 2: reposicionamiento
add("repurposing-robin", 11, "Correr Robin (Future-House/robin) sobre el mecanismo caracterizado. Guardar el ranking crudo de candidatos.", priority=2)
add("repurposing-robin", 11, "Open Targets: asociaciones gen-farmaco para el gen afectado y su red.", priority=2)
add("repurposing-robin", 11, "DGIdb + DrugBank + ChEMBL: interacciones conocidas, estado regulatorio, indicacion actual.", priority=2)
add("repurposing-robin", 11, "LINCS L1000 / CMap: buscar farmacos aprobados que reviertan la firma transcripcional asociada.", priority=2)
add("repurposing-robin", 11, "DepMap: dependencias selectivas en lineas con alta CIN; traducir a oportunidades de letalidad sintetica.", priority=3)
add("repurposing-robin", 11, "TCGA/GDC: firmas de aneuploidia y CIN como contexto (NO como cohorte de pacientes comparables).", priority=4)
add("repurposing-robin", 11, "Filtro de realismo pediatrico: descartar o justificar explicitamente cualquier citotoxico mitotico.", priority=2)
add("repurposing-robin", 11, "Para cada candidato top-5: ficha con mecanismo, evidencia a favor, evidencia en contra, estado regulatorio, seguridad pediatrica y EXPERIMENTO IN VITRO FALSABLE. Sin experimento falsable no entra.", priority=2)

# ---------------------------------------------------------------- Revision adversarial
for i in range(1, 16):
    add("adversarial-reviewer", 12, f"Auditoria #{i}: muestrear 1 de cada 5 outputs de worker. Buscar (a) citas fabricadas, (b) fuga de datos del paciente, (c) saltos logicos sin evidencia, (d) lenguaje que suene a consejo clinico.", priority=2)

# ---------------------------------------------------------------- Reporte y entrega
add("report-writer", 13, "Redactar el reporte Track 2 contra los cuatro encabezados de la rubrica: Rigor 35, Impacto 25, Innovacion 25, Escalabilidad 15.", priority=2)
add("report-writer", 13, "Seccion de Escalabilidad: demostrar el pipeline corriendo sobre OTRO caso no diagnosticado publico. Es 15% del score y casi nadie lo hace.", priority=2)
add("report-writer", 13, "Guion del video pitch de 3 minutos: problema, metodo, hallazgo, que se validaria despues.", priority=3)
add("report-writer", 13, "README del repo publico: reproducibilidad completa, licencia CC BY 4.0, SIN datos de paciente.", priority=2)
add("report-writer", 13, "Checklist final de privacidad antes de cualquier push publico.", priority=1)
add("report-writer", 13, "Consolidar reports/STATUS.md con candidatos, confianza, bloqueos y submissions restantes.", priority=1)

for t in tasks:
    print(json.dumps(t, ensure_ascii=False))

print(f"# total: {len(tasks)} tareas", file=sys.stderr)
