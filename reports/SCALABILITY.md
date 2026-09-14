# Escalabilidad: el mismo pipeline sobre otros genomas

La sección de escalabilidad suele resolverse con una promesa: "el método se
generaliza". Aquí se resuelve corriéndolo.

Se ejecutó el pipeline completo sobre **tres genomas públicos de referencia de
GIAB**, todos de adultos sanos con datos abiertos y consentidos para uso público:

| Muestra | Ascendencia | Por qué |
|---|---|---|
| HG002 (NA24385) | askenazí | el genoma de referencia más usado en benchmarking |
| HG001 (NA12878) | europea | independiente de HG002, sin parentesco |
| HG005 (NA24631) | china (han) | otra ascendencia: las frecuencias de gnomAD cambian por población |

Se eligieron a propósito muestras **no emparentadas**. Los padres de HG002 habrían
sido la opción más cómoda, pero serían controles correlacionados entre sí.

Dos objetivos: demostrar que el método corre sobre otros individuos sin tocar
código, y usarlos como **controles negativos**, que es la pregunta que casi nadie
se hace sobre su propio priorizador: ¿qué devuelve cuando no hay nada que
encontrar?

![Cuatro paneles que comparan al paciente con tres genomas sanos. Variantes raras en el panel y genes con una rara codificante no separan al paciente; genes con dos raras codificantes y patogénicas exactas en ClinVar sí, con 1 en el paciente y 0 en los tres controles](figures/fig4_controles.svg)

## Qué hubo que cambiar

Dos cosas, ninguna de lógica:

```bash
export MVA_WORK=work_hg001                        # separar la salida
python scripts/gene_sweep.py scale/HG001.vcf.gz   # otro VCF de entrada
python scripts/annotate_gnomad.py
python scripts/candidates.py
python scripts/clinvar_annot.py
```

Más un parche para aceptar contigs con prefijo `chr`: el VCF del paciente usa `15`
y los de GIAB usan `chr15`. Es exactamente el tipo de fragilidad que un caso nuevo
saca a la luz y que un pipeline probado sobre un solo genoma nunca revela.

Umbrales, panel de genes y criterios: idénticos en los cuatro genomas.

## Resultado

| | Paciente | HG002 | HG001 | HG005 |
|---|---|---|---|---|
| Variantes en el panel de 32 genes | 2369 | 1968 | 2276 | 2054 |
| Sin anotar por fallo de red | 0 | 0 | 0 | 0 |
| Raras o ausentes de gnomAD | 44 | 45 | 35 | 51 |
| Genes con 2 o más raras | 14 | 10 | 7 | 10 |
| Genes con al menos una rara codificante | 2 | 0 | 1 | **2** |
| Genes con **2 o más** raras codificantes | **1** | 0 | 0 | 0 |
| Coincidencias exactas P/LP en ClinVar | **1** | 0 | 0 | 0 |
| Candidato final | BUB1B | **ninguno** | **ninguno** | **ninguno** |

**Sobre tres genomas sanos, el pipeline devuelve cero candidatos**, sin ajustar
ningún parámetro.

## Lo que los controles enseñan sobre el propio método

No se ve mirando solo al paciente, y además corrige lo que se afirmó con un único
control.

**A nivel de variante rara, los cuatro genomas son indistinguibles.** El paciente
tiene 44; los controles, 35, 45 y 51. Todos tienen varios genes con dos o más
variantes raras.

**Tampoco discrimina tener una variante rara codificante.** Con un solo control
(HG002, que tenía 0) parecía que sí. HG005 tiene 2, igual que el paciente. Una
variante rara en región codificante es algo normal en un genoma sano.

**Lo que discrimina es la conjunción de dos cosas:** dos variantes raras
codificantes en el mismo gen, y patogenicidad conocida de al menos una. Solo el
paciente reúne las dos.

Consecuencia práctica para quien reutilice el método: el esfuerzo va en la
anotación funcional y en el modelo de herencia, no en afinar el umbral de
frecuencia. Y una lección de proceso: una conclusión sacada con n = 1 se cayó en
cuanto hubo n = 3.

## Un fallo que solo apareció gracias a los controles

Al correr HG001 y HG005 por primera vez, salieron con **1.735 y 1.760 variantes
"raras"** en el panel: el 76%, algo biológicamente imposible. La causa era un
defecto real del pipeline. Las consultas a gnomAD fallaban por red, el script lo
avisaba por la salida de errores y seguía, y cada variante del gen quedaba
clasificada como "ausente en gnomAD", es decir, rara. **Un corte de red se
convertía en silencio en evidencia de variante rara**, que es el peor modo de fallo
posible para un priorizador.

Corregido en los tres scripts que consultan gnomAD: reintentos con espera
exponencial, una consulta solo cuenta si termina completa, las variantes de un gen
sin respuesta quedan como `sin_anotacion` y nunca como raras, y el script termina
con código de error para que ningún paso posterior las use. La fila "sin anotar por
fallo de red" de la tabla existe por eso.

El hallazgo del paciente no estaba afectado: sus 32 genes tenían anotación
completa desde la primera corrida, y la repetición con detección de fallos da el
mismo resultado.

## Límites de esta demostración

1. **n = 3 controles** no estima una tasa de falsos positivos con precisión. Para
   eso harían falta decenas de controles, o casos sintéticos con variantes
   causales conocidas insertadas.
2. **Los VCF no son comparables en detalle.** El del paciente viene de GATK sobre
   el genoma completo; los de GIAB son conjuntos de alta confianza restringidos a
   regiones benchmark. Los conteos absolutos no se comparan; lo comparable son los
   puntos de decisión.
3. **Son adultos sanos, no casos sin diagnosticar.** Demostrar utilidad real
   exigiría correr el pipeline sobre casos sin resolver y medir cuántos resuelve.
4. **Mismo panel.** Se probó la generalización a otros individuos y ascendencias,
   no a otra enfermedad. Cambiar de enfermedad es cambiar `ref/sac_genes.txt`, pero
   eso no está demostrado aquí.

## Coste y requisitos

El pipeline completo sobre un genoma nuevo corre en un portátil con 4 núcleos y
8 GB de RAM, con Python de la biblioteca estándar. Sin clúster, sin nube, sin
contenedores, sin GPU.

- Coste de cómputo: **0 dólares**, para los cuatro genomas.
- Coste de LLM sobre datos del paciente: **0 dólares** en la vía determinista.
- Descargas: recursos públicos de referencia (MANE, ClinVar, chr15), más consultas
  por rango de bytes a gnomAD que traen unos pocos MB en vez de los 7.4 GB por
  cromosoma del archivo completo.

Un laboratorio sin infraestructura puede correrlo sobre el genoma de su paciente
sin pedirle presupuesto a nadie y sin que los datos salgan de su máquina.
