# Track 2 — Mecanismo y reposicionamiento

Rare Disease, Real Kid: The MVA Hackathon 2026.

> **Esto no es consejo médico ni un plan de tratamiento.** Ninguna molécula de
> este informe tiene evidencia de eficacia en MVA. Son hipótesis derivadas del
> mecanismo, cada una acompañada del experimento que la refutaría. Nada de esto
> debe llegar a un paciente sin pasar antes por ese experimento y por un ensayo
> clínico.

Estructurado contra los cuatro encabezados de la rúbrica.

---

## 1. Rigor científico (35%)

### El mecanismo, eslabón por eslabón

Cadena completa en `MECHANISM_BUB1B.md`, con 16 citas verificadas contra PubMed.
En resumen:

**Variante → función.** BUB1B codifica BUBR1, componente del complejo de
checkpoint mitótico que se cataliza en el cinetocoro (PMID 36518060). El alelo
p.Leu737Ter introduce un codón de parada en el exón 17 de 23, con degradación por
NMD predicha: es **pérdida de función**, y elimina por completo el dominio
766-1050.

**El segundo alelo.** p.Asn1002Lys cae dentro de ese dominio. UniProt lo anota
como "Protein kinase", pero está demostrado que es un **pseudoquinasa sin
actividad catalítica** cuyo papel es reclutar la fosfatasa PP2A-B56 al cinetocoro
(PMID 33207204, PMID 35525552). La hipótesis mecanística, por tanto, no es pérdida
de catálisis: es interferencia con ese andamiaje.

**Función → fenotipo.** La dosis de BUBR1 se comporta como un reóstato, no como un
interruptor: niveles reducidos se asocian a aneuploidía (PMID 18699967) y
aumentarlos protege frente a aneuploidía y cáncer en ratón (PMID 23242215). Un
checkpoint débil produce errores de segregación independientes célula a célula,
lo que genera un mosaico de cromosomas distintos —aneuploidía variegada— en lugar
de una trisomía uniforme. El cuadro clínico resultante incluye microcefalia y
retraso del crecimiento (PMID 18548531), y riesgo tumoral pediátrico con
rabdomiosarcoma embrionario y tumor de Wilms descritos en MVA (PMID 9916837,
PMID 42595739). La causa genética bialelica está establecida desde
PMID 15475955, con el patrón de alelo truncante más alelo hipomórfico documentado
en siete familias (PMID 16411201).

### Qué sostiene la hipótesis y qué no

**A favor:** el gen correcto para el fenotipo; un alelo patogénico conocido con
NMD verificado contra la estructura exónica; un segundo alelo prácticamente
privado en el mismo gen; ausencia de homocigotos en gnomAD para ambos; calidad
alta en las dos llamadas (GQ 99).

**En contra, o sin resolver:**

1. **La fase no está demostrada.** Los bloques de fase del VCF llegan a 206 pb y
   las variantes están a 10.911 pb. Si estuvieran en cis, el paciente tendría un
   alelo intacto y toda la hipótesis se cae.
2. **El segundo alelo es VUS** por ACMG, y lo sigue siendo aunque se resuelva la
   fase: hace falta evidencia funcional.
3. **La aneuploidía no se observa en estos datos.** El barrido de BAF no detecta
   eventos por encima de ~10-15% de fracción celular.
4. **Buena parte de la biología de consecuencias viene de modelos tumorales**, y
   extrapolar de una célula seleccionada por tolerar aneuploidía a un tejido en
   desarrollo es un salto que nombramos cada vez que lo damos.

### Control de calidad del propio método

Se documentan dos errores encontrados y corregidos durante el análisis —un test
de BAF que producía fracciones celulares idénticas en los 24 cromosomas, y un
cribado de splicing que daba positivo en todo— junto con la tasa de fondo medida
del método que los sustituyó. Detalle en `TRACK1_METHODS.md`.

---

## 2. Impacto potencial (25%)

### Lo que ya existe

La única intervención con respaldo hoy en MVA no es un fármaco: son los
protocolos de **vigilancia oncológica** en niños con predisposición genética
(PMID 39264246). Cualquier hipótesis farmacológica compite por atención y
recursos con ese estándar, y no lo sustituye. Decirlo primero es parte del rigor.

### Qué cambiaría para este niño

Un diagnóstico molecular confirmado no es un tratamiento, pero cambia cosas
concretas: fija la intensidad y el calendario de la vigilancia tumoral, permite
consejo genético a la familia, y habilita el reclutamiento en registros y
ensayos. La confirmación de fase por genotipado parental es una prueba barata con
rendimiento diagnóstico alto, y es la recomendación operativa de este trabajo.

### Los candidatos, con su confianza

Fichas completas en `REPURPOSING.md`, con evidencia a favor, evidencia en contra,
seguridad pediátrica y experimento falsable.

| Candidato | Eje | Confianza | Por qué |
|---|---|---|---|
| Eje senolítico | carga de células senescentes | 0.35 | la intervención se probó en el ratón progeroide **BubR1** (PMID 22048312): el modelo del gen correcto |
| Sirolimus | recambio proteico y tolerancia a aneuploidía | 0.25 | tolerancia a aneuploidía ligada a recambio de proteína (PMID 38778096); amplio historial pediátrico (PMID 26783326) |
| N-acetilcisteína | carga oxidativa | 0.15 | el estrés oxidativo induce CIN en fibroblastos (PMID 37288672); vínculo indirecto con BUB1B |

Ninguna confianza llega a la mitad, y es deliberado: **ninguno tiene un solo dato
en MVA**.

El contraargumento más fuerte es específico y merece destacarse: el sirolimus
inhibe mTOR y frena el crecimiento, y el retraso del crecimiento ya forma parte
del fenotipo. El fármaco empujaría en la misma dirección que la enfermedad.

### Filtro de realismo pediátrico

| Clase descartada | Motivo |
|---|---|
| Inhibidores de Aurora A/B | Se usan **contra** tumores con CIN. Dárselos a un niño cuyas células ya missegregan empuja hacia el daño |
| Inhibidores de MPS1/TTK | TTK pertenece al mismo checkpoint que ya está debilitado: agrava el defecto primario |
| Antimitóticos citotóxicos | Toxicidad de quimioterapia en indicación no oncológica, en un paciente con predisposición tumoral |
| Letalidad sintética con CIN | Estrategia para matar células tumorales; aquí las células diana son las del propio niño |

Estos candidatos aparecen espontáneamente en cualquier búsqueda de literatura
sobre aneuploidía —alisertib salió en las nuestras— y son exactamente los que un
panel experto rechazaría. El filtro forma parte del método, no del prólogo.

---

## 3. Innovación (25%)

**No hay precedente de reposicionamiento en MVA.** La búsqueda cruzando MVA con
tratamiento, terapia o fármaco no devuelve ningún intento previo. Es terreno
virgen, y conviene leer eso en las dos direcciones: es una oportunidad, y también
una señal de que justificarlo es difícil.

Tres aportes metodológicos que sobreviven al caso concreto:

1. **Anotación genómica sin descargar los recursos.** Un cliente tabix por HTTP en
   Python puro consulta gnomAD por rangos de bytes: unos pocos MB en vez de 7.4 GB
   por cromosoma. Elimina la necesidad de clúster o nube, y de paso hace que lo
   que sale de la máquina sea la coordenada de un gen del panel, nunca una del
   paciente.
2. **Control de egreso verificable.** Un hook que bloquea coordenadas, rsIDs, HGVS
   e identificadores de muestra hacia servicios externos, con prueba de red-team
   de 30 canarios y 10 consultas legítimas. La privacidad deja de ser una promesa
   en el README y pasa a ser un test que corre.
3. **Reencuadre del segundo alelo.** Tratar el dominio 766-1050 como pseudoquinasa
   de andamiaje en vez de como quinasa cambia la hipótesis funcional y cambia el
   experimento que la falsa: no se mide actividad catalítica, se mide
   reclutamiento de PP2A-B56.

---

## 4. Escalabilidad (15%)

Demostrada corriendo, no prometida. Detalle en `SCALABILITY.md`.

El pipeline completo se ejecutó sobre **HG002 (GIAB)**, un genoma público de un
adulto sano, cambiando solo el VCF de entrada y el directorio de salida. Ningún
umbral, ningún parámetro, ningún criterio se ajustó.

| | Paciente | HG002 |
|---|---|---|
| Raras o ausentes en el panel | 44 | **45** |
| Genes con 2 o más raras | 14 | 10 |
| Genes con al menos una rara codificante | 2 | **0** |
| Coincidencias P/LP en ClinVar | 1 | **0** |
| Candidato final | BUB1B | **ninguno** |

El pipeline no fabrica un diagnóstico cuando no lo hay.

Y el control enseña algo sobre el método que no se ve mirando solo al paciente:
un adulto sano carga **45 variantes raras** en genes del checkpoint mitótico, una
más que el niño. A ese nivel los dos genomas son indistinguibles. Toda la
discriminación vive en la anotación de consecuencia y en el cruce con ClinVar. Un
pipeline que se detuviera en el filtro de frecuencia habría entregado 14 genes
candidatos para el paciente y 10 para una persona sana. **El esfuerzo va en la
anotación funcional, no en afinar el umbral de frecuencia.**

Requisitos para reutilizarlo: Python de la biblioteca estándar, un portátil de 4
núcleos y 8 GB de RAM. Sin nube, sin contenedores, sin GPU. Coste de cómputo: 0
dólares. Un laboratorio sin infraestructura puede correrlo sobre el genoma de su
paciente sin que los datos salgan de su máquina.

**Límites de la demostración:** un solo control no estima una tasa de falsos
positivos; HG002 es un adulto sano, no un caso sin diagnosticar; y se probó la
generalización a otro individuo, no a otra enfermedad.

---

## Qué haría falta a continuación

En orden de coste y de rendimiento esperado:

1. **Genotipar a los padres.** Barato, resuelve la fase, y es lo que convierte la
   hipótesis en diagnóstico.
2. **SpliceAI o Pangolin** sobre la variante intrónica profunda.
3. **Ensayo funcional del missense**: reclutamiento de PP2A-B56 al cinetocoro en
   células que expresen p.Asn1002Lys.
4. **Caracterización celular del paciente**: carga senescente, frecuencia de
   micronúcleos, tasa de separación prematura de cromátidas. Es el experimento que
   discrimina entre los tres candidatos de reposicionamiento y el que, si sale
   negativo, mata el más fuerte de los tres antes de gastar un euro en fármacos.

## Divulgación de uso de IA

El análisis lo ejecutaron scripts deterministas. Un modelo de lenguaje (Claude,
Anthropic) escribió los scripts, interpretó salidas agregadas y redactó los
informes. El modelo recibió símbolos de genes, consecuencias a nivel de proteína,
métricas de calidad, frecuencias y clasificaciones; no recibió coordenadas
genómicas, tablas de genotipos ni identificadores de muestra. Todas las citas de
este informe se resolvieron contra PubMed durante la propia sesión; ninguna se
escribió de memoria.
