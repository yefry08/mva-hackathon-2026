# Escalabilidad: el mismo pipeline sobre otro genoma

La sección de escalabilidad suele resolverse con una promesa: "el método se
generaliza". Aquí se resuelve corriéndolo.

Se ejecutó el pipeline completo sobre **HG002 (NA24385)**, el genoma de referencia
público de GIAB, un adulto sano cuyos datos son abiertos y consentidos para uso
público. Dos motivos para elegirlo: demuestra que el método corre sobre otro
individuo sin tocar código, y funciona como **control negativo**, que es la
pregunta que nadie se hace sobre su propio priorizador: ¿qué devuelve cuando no
hay nada que encontrar?

## Qué hubo que cambiar

Dos cosas, ninguna de lógica:

```bash
export MVA_WORK=work_hg002                       # separar la salida
python scripts/gene_sweep.py scale/HG002.vcf.gz  # otro VCF de entrada
python scripts/annotate_gnomad.py
python scripts/candidates.py
python scripts/clinvar_annot.py
```

Más un parche de una línea para aceptar contigs con prefijo `chr`, porque el VCF
del paciente usa `15` y el de GIAB usa `chr15`. Ese parche es exactamente el tipo
de fragilidad que un caso nuevo saca a la luz y que un pipeline probado sobre un
solo genoma nunca revela.

Umbrales, panel de genes y criterios: idénticos. Ningún parámetro se ajustó.

## Resultado

| | Paciente | HG002 (control) |
|---|---|---|
| Variantes en el panel de 32 genes | 2369 | 1968 |
| Raras o ausentes de gnomAD | 44 | **45** |
| Genes con 2 o más raras | 14 | 10 |
| Genes con al menos una rara **codificante** | 2 | **0** |
| Genes con 2 o más raras codificantes | 1 | **0** |
| Coincidencias exactas P/LP en ClinVar | 1 | **0** |
| Candidato final | BUB1B | **ninguno** |

**El pipeline no inventa un diagnóstico cuando no lo hay.** Sobre un genoma sano
devuelve cero candidatos, sin necesidad de ajustar nada.

## Lo que el control enseña sobre el propio método

Es el hallazgo más útil de este ejercicio, y no se ve mirando solo al paciente.

Un adulto sano de referencia carga **45 variantes raras** en genes del checkpoint
mitótico, una más que el paciente. En 10 de esos genes tiene dos o más. A nivel de
"variante rara en gen candidato", los dos genomas son **indistinguibles**.

Toda la capacidad de discriminación vive en dos pasos posteriores: la anotación
de consecuencia, que separa codificante de intrónico, y el cruce con ClinVar. Un
pipeline que se detuviera en el filtro de frecuencia habría entregado 14 genes
candidatos para el paciente, y habría entregado 10 para una persona sana.

Esto tiene una consecuencia práctica para quien reutilice el método: **el
presupuesto de esfuerzo va en la anotación funcional, no en afinar el umbral de
frecuencia.**

## Límites de esta demostración

1. **n = 1 control.** Un solo genoma sano no estima una tasa de falsos positivos.
   Para eso harían falta decenas de controles, o casos sintéticos con variantes
   causales conocidas insertadas, que es el diseño que quedó planificado y sin
   correr.
2. **Los VCF no son comparables en detalle.** El del paciente viene de GATK sobre
   el genoma completo; el de HG002 es un conjunto de alta confianza restringido a
   regiones benchmark. Los conteos absolutos no se pueden comparar entre sí; lo
   comparable son los puntos de decisión, que es lo que la tabla usa.
3. **HG002 es un adulto sano.** No es un caso no diagnosticado. Demostrar utilidad
   real exigiría correrlo sobre casos sin resolver y medir cuántos se resuelven.
4. **Mismo panel.** Se probó la generalización a otro individuo, no a otra
   enfermedad. Cambiar de enfermedad es cambiar `ref/sac_genes.txt`, pero eso no
   está demostrado aquí.

## Coste y requisitos

El pipeline completo sobre un genoma nuevo corre en un portátil con 4 núcleos y
8 GB de RAM, con Python de la biblioteca estándar. Sin clúster, sin nube, sin
contenedores, sin GPU.

- Coste de cómputo: **0 dólares**.
- Coste de LLM: **0 dólares** en la vía determinista.
- Descargas: recursos públicos de referencia (MANE, ClinVar, chr15), más consultas
  por rango de bytes a gnomAD que traen unos pocos MB en vez de los 7.4 GB por
  cromosoma del archivo completo.

Esa es la parte escalable de verdad: un laboratorio sin infraestructura puede
correr esto sobre el genoma de su paciente sin pedirle presupuesto a nadie y sin
que los datos salgan de su máquina.
