# Reglas del hackathon — hallazgos (leídas 2026-09-12)

Fuentes (todas públicas):
- Reglas oficiales: `tabs/rules.py` del Space `SageBio/rare-disease-real-kid-mva-hackathon-2026`
- FAQ: `tabs/faq.py` · Scoring: `evaluation.py` · Formato: `tabs/submit_track1.py`, `tabs/submit_track2.py`, `config.py`
- Discusión #2 (uso de LLMs de terceros), #17 (build de referencia), #21 (suscripción de consumidor, **sin responder**), #22 (cis/trans)

---

## 1. LLMs de terceros: permitidos bajo condiciones

La respuesta la dio la Chief Privacy and Compliance Officer de Sage en la discusión #2, y el
presidente de Sage la reafirmó. La prueba que aplican es **Processor vs Recipient**:

- Un servicio que procesa los datos solo para devolverte un resultado, no toma derechos sobre
  ellos y no puede usarlos para fines propios, es una **herramienta** — igual que un almacenamiento
  en la nube, un alineador hospedado o una API de anotación.
- Un servicio que gana el derecho de usar, retener o aprender de los datos es un **receptor**, y
  mandarle datos es una divulgación.

Dos condiciones que deben cumplirse a la vez:
1. Sin entrenamiento sobre tus inputs ni outputs, y sin que el proveedor tome derechos sobre ellos.
2. Retención limitada en tiempo y propósito. Cero retención no es obligatorio: logs cortos para
   abuso, depuración o calidad de servicio son aceptables si el contenido no se usa para otra cosa.

Advertencias explícitas de los organizadores:
- Desactivar el entrenamiento no siempre es completo: algunos proveedores sí usan el contenido
  sobre el que das feedback. **Hay que optar por salir y no calificar outputs** (nada de pulgares).
- Los programas de créditos pueden traer términos propios que anulan los de tu cuenta. Hay que leer
  los términos de los créditos que uses.
- Hay que registrar proveedor, plan y ajuste en la descripción de métodos. El ejemplo que dan es una
  línea del tipo: "Anthropic API, Claude <modelo>, términos comerciales, sin entrenamiento sobre
  contenido del cliente".

La obligación de borrado **no** alcanza a los logs del proveedor que no controlas. No se exige
inferencia local.

### Lo que esto implica para nosotros
- Los workers de Claude **sí pueden** ver datos a nivel de variante, siempre que corran bajo términos
  de procesador. La vía limpia es una **API key comercial** (los términos comerciales de Anthropic no
  entrenan sobre contenido del cliente).
- La discusión #21 pregunta justamente si una suscripción de consumidor (Claude Max) con
  "Help improve our AI models" apagada cumple. **Nadie la ha respondido.** Si corremos con
  suscripción, estamos en la zona sin resolver; con API key comercial, no.
- FutureHouse/Edison y cualquier otro servicio pasan la misma prueba. Como no verifiqué sus
  términos, se mantienen en nivel gen/mecanismo, que además es lo que Robin necesita.

## 2. Borrado y atestación (esto hay que diseñarlo desde ahora)

Todo lo que lleve el genoma del niño se borra dentro de los 30 días del cierre (cierre 2026-10-24
→ borrado antes del ~2026-11-23), de **todos** los entornos, y se confirma por correo a
`RarediseaserealkidMVAhackathon2026@synapse.org`.

**Se borra:** VCF, BAM, CRAM y cualquier copia, subconjunto, recorte o reformateo, índices incluidos ·
archivos intermedios con genotipos por variante a escala de tabla · cachés y estado de notebooks ·
**prompts o logs guardados en nuestros sistemas que contengan bloques de datos de variantes** ·
pesos, embeddings o fine-tunes entrenados sobre los datos crudos.

**Se conserva:** la lista rankeada de candidatos (es la submission) · los términos HPO · rankings de
genes y vías, mecanismo, candidatos de fármacos · código, reporte, pitch y entrada del leaderboard.

Regla práctica que dan: un puñado de variantes nombradas en un reporte es un hallazgo; una tabla de
genotipos a escala genómica es el dataset en otro formato.

> **Consecuencia operativa que casi nadie va a ver:** Claude Code guarda los transcripts de sesión en
> disco (`~/.claude/projects/...`). Si un worker lee tablas de variantes, esos transcripts son
> "logs con bloques de datos de variantes" y entran en la lista de borrado. Hay que mantenerlos en el
> mismo volumen purgable que `data/` y `work/`, y meterlos en el `make purge`.

Otras obligaciones: no re-contactar a la familia ni a la MVA Society · no redistribuir los datos por
ningún canal · el repo puede ser privado durante el hackathon pero debe ser público al cerrar ·
las submissions salen con licencia CC BY · embargo de publicación revisada por pares hasta que los
organizadores publiquen su reporte; código y salidas derivadas se pueden compartir cuando queramos ·
hay un texto de agradecimiento obligatorio en cualquier publicación.

## 3. Datos y referencia

- Formato: VCF; datos crudos **opcionalmente** en BAM/CRAM; fenotipo ya viene como términos HPO
  estandarizados. ~85 GB, un solo sujeto. En la discusión #17 se menciona que hay FASTQ.
- Build: **GRCh38**. El header del VCF apunta a
  `GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_hs38d1_maskedGRC_exclusions_v2_no_chr.fasta`:
  contigs **sin prefijo `chr`** (1, 2, X), con decoys hs38d1 y exclusiones enmascaradas.
- El VCF trae campos **PGT/PID**, es decir, fase parcial ya calculada.
- El formulario de submission usa `chrom` con prefijo (`chr15`). **Hay que convertir la
  nomenclatura de contigs al generar el CSV**, y no mezclar el analysis set genérico con el fasta
  exacto del header si llegamos a realinear.

## 4. Scoring y formato exacto (Track 1)

CSV con: `proband_id, chrom_1, pos_1, ref_1, alt_1, chrom_2, pos_2, ref_2, alt_2, epcr` y los
opcionales `finding_type` (`primary`/`secondary`) y `notes`. Coordenadas GRCh38.

- **Máximo 10 filas por probando.** Variante simple: los campos `_2` van vacíos. Par heterocigoto
  compuesto: las dos variantes en la misma fila.
- `epcr` en el rango (0, 1]. Las filas se ordenan por EPCR descendente antes de puntuar; empates por
  orden de envío.
- **Rank points:** rank 1 = 100 · ranks 2–3 = 50 · ranks 4–5 = 25 · ranks 6–10 = 10 · >10 = 0.
  Coincidencia parcial en un compuesto (recuperas una de las dos) = mitad de los puntos de ese rank.
- **F-max:** a nivel de variante individual, barriendo cada umbral de EPCR presente en el envío y
  quedándose con el mejor F.
- Hallazgos secundarios no penalizan el score automático; van a revisión cualitativa del panel.

**Consecuencia de estrategia:** como F-max toma el *máximo* sobre los umbrales, llenar las 10 filas
no hace daño **siempre que las variantes verdaderas tengan los EPCR más altos**. Lo que hunde el
score es una variante equivocada por encima de la correcta. El trabajo está en el *orden* y en la
separación de los EPCR, no en el volumen.

- Límite de submissions: `config.py` dice 6 (Track 1) y 3 (Track 2). El texto de la pestaña de envío
  dice 5 y muestra el mensaje de "ya usaste las 5". **Planeamos con 5** y confirmamos antes del
  último envío.

## 5. Track 2: entregables

Reporte (PDF o Markdown, con el nombre de usuario en el archivo) + repo de GitHub público al cierre +
video de 3 minutos en YouTube/Vimeo. Solo revisan la **última** entrega. El reporte debe caracterizar
el mecanismo de la variante (pérdida o ganancia de función, vía afectada, consecuencia biológica) y
**declarar el uso de LLM/IA**. Rúbrica: Rigor 35 · Impacto 25 · Innovación 25 · Escalabilidad 15.
Premios: 50.000 dólares en total; 1er lugar 12.000 en efectivo más 12.000 en créditos de Claude.

## 6. El dato incómodo

El FAQ admite que ya hay **puntajes perfectos** en el leaderboard de Track 1 y lo llama un track
"fundacional", diseñado para ser alcanzable. Recuperar la variante no es la meta: el panel juzga el
methods write-up. Traducción: no gastemos el presupuesto en cómputo exótico para Track 1; gastémoslo
en rigor documentado y en Track 2, que es donde se decide.
