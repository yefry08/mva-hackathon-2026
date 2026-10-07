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

**La esencialidad predice qué segundo alelo cabe esperar.** Cruzando los tres
genes MVA con DepMap (1.258 líneas celulares): CEP57 y TRIP13 no son esenciales
(0% y 1% de las líneas dependen de ellos) y en sus pacientes se describen dobles
nulos (PMID 30035751, PMID 28553959). BUB1B es esencial en el 77% de las líneas, y
lo descrito es un nulo más un alelo que conserva función. El paciente tiene
exactamente eso: un nonsense y un missense. Son tres genes y líneas tumorales, así
que es coherencia mecanística y no una ley, ni un criterio ACMG; pero es una
predicción que el caso cumple y que no se buscó para que la cumpliera.

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

![Histograma logarítmico de los bloques de fase del cromosoma 15: el más largo mide 206 pares de bases y las dos variantes de BUB1B están a 10.911, fuera del alcance de cualquier lectura](figures/fig2_fase.svg)

### Control de calidad del propio método

Se documentan cuatro errores encontrados y corregidos durante el análisis: un test
de BAF que producía fracciones celulares idénticas en los 24 cromosomas; un
cribado de splicing que daba positivo en todo; un fallo de red que se convertía en
silencio en "variante rara", destapado por los genomas de control; y un control de
privacidad ciego a los CSV. Ninguno cambió el hallazgo del paciente. Detalle en
`TRACK1_METHODS.md`.

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

**Y el filtro se comprobó con datos, no solo con argumentos.** Consultando DGIdb y
Open Targets para los 32 genes del panel, solo con símbolos de gen: **123
interacciones fármaco-gen, todas inhibidoras salvo una vacuna, y ningún
activador**. Solo 5 genes tienen candidatos en fase clínica, todos inhibidores
oncológicos, y BUB1B no tiene ninguno. No existe nada, aprobado ni experimental,
que empuje esta vía en la dirección que necesita un paciente con poco checkpoint.
Por eso los tres candidatos actúan aguas abajo, sobre las consecuencias de la
aneuploidía: aguas arriba no hay dónde apoyarse. Detalle en `REPURPOSING.md`.

**El filtro también se probó contra una búsqueda que no parte de la literatura.**
La conectividad transcripcional (LINCS L1000) busca compuestos que inviertan la
firma de una célula sin BUB1B. La corrimos con tres controles: 40 knockouts al azar
como nulo, una réplica por shRNA limitada a las líneas donde BUB1B baja de verdad,
y la comprobación de en qué genes se apoya cada coincidencia. **Ningún candidato
sobrevive:**

- Los más repetidos, inhibidores de HSP90 y de CDK, revierten igual la pérdida de
  genes ajenos al checkpoint (15-42% del nulo).
- Los dos aprobados que parecían específicos, pentobarbital y naltrexona, salen de
  una sola placa y coinciden por dos genes: es un efecto de lote.

Lo revelador es lo que sale **sin** esos controles: mebendazol, albendazol,
vincristina, docetaxel. Todos aprobados, y el mebendazol es barato, pediátrico y
un clásico del reposicionamiento. Son venenos del huso, y una célula con BubR1
bajo no sostiene la parada mitótica ante ellos (PMID 23812934). Un pipeline sin
filtro le habría propuesto a este niño un fármaco que empuja hacia más
missegregación. El resultado negativo, con sus razones, es el aporte.

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
4. **Una firma de referencia de la pérdida de BUB1B, con su nulo.** El knockout
   CRISPR activa interferón, p53 y apoptosis con una fuerza que solo alcanza el
   0,1-3% de las 5.211 firmas de knockout restantes. El interferón lo comparten
   BUB1, ZWINT y TTK, que señalizan el checkpoint, y no AURKB ni PLK1. No se
   replica con shRNA, lo que deja abierta una pregunta de dosis que el RNA-seq de
   fibroblastos del paciente puede responder: ¿se parecen al knockout o al
   knockdown? Detalle en `MECHANISM_BUB1B.md`.

---

## 4. Escalabilidad (15%)

Demostrada corriendo, no prometida. Detalle en `SCALABILITY.md`.

El pipeline completo se ejecutó sobre **tres genomas públicos de adultos sanos**
de GIAB, no emparentados y de tres ascendencias (askenazí, europea y china),
cambiando solo el VCF de entrada y el directorio de salida. Ningún umbral, ningún
parámetro, ningún criterio se ajustó.

![Cuatro paneles que comparan al paciente con tres genomas sanos. Variantes raras y genes con una rara codificante no separan al paciente; solo él tiene dos raras codificantes en un mismo gen y una patogénica exacta en ClinVar](figures/fig4_controles.svg)

| | Paciente | HG002 | HG001 | HG005 |
|---|---|---|---|---|
| Raras o ausentes en el panel | 44 | 45 | 35 | 51 |
| Genes con al menos una rara codificante | 2 | 0 | 1 | **2** |
| Genes con 2 o más raras codificantes | **1** | 0 | 0 | 0 |
| Coincidencias P/LP en ClinVar | **1** | 0 | 0 | 0 |
| Candidato final | BUB1B | **ninguno** | **ninguno** | **ninguno** |

El pipeline no fabrica un diagnóstico cuando no lo hay, en ninguno de los tres.

Los controles enseñan algo que no se ve mirando solo al paciente, y además
corrigen lo que se afirmó con uno solo. A nivel de variante rara los cuatro genomas
son indistinguibles: 44 en el paciente, entre 35 y 51 en los controles. Con un
único control parecía que tener una variante rara codificante ya discriminaba;
**HG005 tiene dos, igual que el paciente.** Lo que discrimina es la conjunción de
dos variantes raras codificantes en el mismo gen con patogenicidad conocida. **El
esfuerzo va en la anotación funcional y el modelo de herencia, no en afinar el
umbral de frecuencia.**

Los controles también destaparon un defecto real del pipeline —un fallo de red
que se convertía en silencio en "variante rara"—, ya corregido y documentado en
`SCALABILITY.md`. El hallazgo del paciente no estaba afectado.

Requisitos para reutilizarlo: Python de la biblioteca estándar, un portátil de 4
núcleos y 8 GB de RAM. Sin nube, sin contenedores, sin GPU. Coste de cómputo: 0
dólares. Un laboratorio sin infraestructura puede correrlo sobre el genoma de su
paciente sin que los datos salgan de su máquina.

**Límites de la demostración:** tres controles no estiman con precisión una tasa
de falsos positivos; son adultos sanos, no casos sin diagnosticar; y se probó la
generalización a otros individuos y ascendencias, no a otra enfermedad.

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
   Añadir RNA-seq de esas mismas células cuesta poco y se puntúa directamente
   contra la firma de `cmap_signature.tsv`.

## Divulgación de uso de IA

El análisis lo ejecutaron scripts deterministas. Un modelo de lenguaje (Claude,
Anthropic) escribió los scripts, interpretó salidas agregadas y redactó los
informes. El modelo recibió símbolos de genes, consecuencias a nivel de proteína,
métricas de calidad, frecuencias y clasificaciones; no recibió coordenadas
genómicas, tablas de genotipos ni identificadores de muestra. Todas las citas de
este informe se resolvieron contra PubMed durante la propia sesión; ninguna se
escribió de memoria.
