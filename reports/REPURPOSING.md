# Reposicionamiento: hipótesis falsables, no recomendaciones

**Esto no es consejo médico.** Ninguna de estas moléculas tiene evidencia de
eficacia en MVA. Son hipótesis derivadas del mecanismo, cada una con el
experimento concreto que la refutaría. Nada de esto debe llegar a un paciente sin
pasar antes por ese experimento y por un ensayo clínico.

Citas verificadas contra PubMed el 2026-09-12.

---

## Lo que ya existe, y que hay que decir primero

La única intervención con respaldo hoy en MVA no es un fármaco: son los
**protocolos de vigilancia oncológica** en niños con predisposición tumoral
genética (PMID 39264246). Cualquier hipótesis farmacológica compite con ese
estándar en atención y recursos, y no lo sustituye.

**Precedente de reposicionamiento en MVA: no encontramos ninguno.** La búsqueda
en PubMed cruzando MVA con tratamiento, terapia o fármaco no devuelve ningún
intento previo. Es terreno virgen, lo que es a la vez la oportunidad y la señal
de alarma: puede que nadie lo haya intentado porque es difícil justificarlo.

---

## Candidato 1 — Eliminación de células senescentes (eje senolítico)

```yaml
claim: "La carga de células senescentes es un efector del fenotipo por
        hipomorfismo de BubR1, y reducirla podría mitigar manifestaciones
        somáticas sin tocar la tasa de missegregación."
evidence_for:
  - source: PMID:22048312
    finding: "En el fondo del ratón progeroide BubR1, eliminar células positivas
              para p16Ink4a retrasa trastornos asociados al envejecimiento en
              tejido adiposo, músculo esquelético y ojo."
    strength: strong
    model: mouse
    nota: "Es el modelo del gen correcto. Esta es la coincidencia más cercana que
           existe entre BUB1B y una intervención."
  - source: PMID:34995493
    finding: "La aneuploidía constitucional (trisomía 21) induce senescencia en
              progenitores neurales y altera su arquitectura nuclear."
    strength: moderate
    model: human_cells
    nota: "Aneuploidía constitucional, no tumoral: el modelo correcto para un
           fenotipo del neurodesarrollo."
evidence_against:
  - finding: "INK-ATTAC es una herramienta genética inducible, no un fármaco. El
              salto de ahí a un senolítico (dasatinib más quercetina,
              navitoclax) es un salto real y no está hecho en este modelo."
  - finding: "Eliminar células senescentes en un organismo en desarrollo no es lo
              mismo que hacerlo en uno envejecido. La senescencia tiene papeles
              fisiológicos en el desarrollo, y aquí el paciente es un niño."
  - finding: "Dasatinib en pediatría se usa en leucemia Ph+, con toxicidad
              relevante. El perfil riesgo-beneficio en una indicación no
              oncológica es otro problema."
confidence: 0.35
falsifiable_by: >
  Fibroblastos del paciente frente a control pareado por edad: cuantificar carga
  senescente (SA-beta-gal, p16, p21) y perfil SASP. Si la carga no está elevada
  respecto del control, la hipótesis muere ahí, antes de hablar de fármacos. Si
  lo está, probar si un senolítico elimina selectivamente esa población sin
  afectar la fracción proliferante, y confirmar que la frecuencia de micronúcleos
  NO cambia: la senolisis no debería corregir la missegregación, y si lo hiciera,
  el mecanismo que propusimos está mal.
```

---

## Candidato 2 — Sirolimus (inhibición de mTOR, recambio proteico)

```yaml
claim: "Potenciar el recambio proteico podría aumentar la tolerancia celular al
        desbalance estequiométrico que impone la aneuploidía."
evidence_for:
  - source: PMID:38778096
    finding: "La diversidad natural del proteoma vincula la tolerancia a
              aneuploidía con el recambio de proteínas."
    strength: moderate
    model: cell_line
  - source: PMID:39247952
    finding: "Contrarrestar la carga transcripcional de la aneuploidía exige
              aumentar la degradación de ARN y proteína."
    strength: moderate
    model: cell_line
  - source: PMID:27308438
    finding: "La aneuploidía impone estrés proteotóxico."
    strength: moderate
    model: cell_line
  - source: PMID:26783326
    finding: "Sirolimus en anomalías vasculares complicadas pediátricas: eficacia
              y seguridad documentadas."
    strength: strong
    model: human
  - source: PMID:34524406
    finding: "Sirolimus en malformaciones de flujo lento en niños, fase
              observacional."
    strength: moderate
    model: human
evidence_against:
  - finding: "Es inmunosupresor, y este niño tiene predisposición tumoral. Reducir
              la vigilancia inmunitaria en ese contexto es un riesgo direccional,
              no teórico."
  - finding: "La inhibición de mTOR frena el crecimiento, y el retraso del
              crecimiento ya es parte del fenotipo de MVA. El fármaco empujaría en
              la misma dirección que la enfermedad."
  - finding: "Toda la evidencia de tolerancia a aneuploidía viene de líneas
              celulares seleccionadas por sobrevivir con aneuploidía. Un tejido en
              desarrollo no pasó por esa selección."
confidence: 0.25
falsifiable_by: >
  Fibroblastos del paciente: medir flujo autofágico (LC3-II con y sin bafilomicina)
  e inducción de chaperonas basal frente a control. Con rapamicina, comprobar si la
  viabilidad y el estrés proteotóxico mejoran a concentraciones que no frenen la
  proliferación. Si el beneficio solo aparece a dosis que detienen el ciclo, no hay
  ventana terapéutica y el candidato se cae.
```

---

## Candidato 3 — N-acetilcisteína (carga oxidativa)

```yaml
claim: "Reducir el estrés oxidativo podría bajar el componente de inestabilidad
        cromosómica que depende de estrés replicativo, sin tocar el checkpoint."
evidence_for:
  - source: PMID:37288672
    finding: "El estrés oxidativo induce inestabilidad cromosómica a través de
              estrés replicativo en fibroblastos."
    strength: moderate
    model: human_cells
    nota: "Fibroblastos, no tumor. Modelo mejor emparejado que la mayoría."
  - source: PMID:36384314
    finding: "N-acetilcisteína en niños y adolescentes: uso crónico con
              tolerabilidad documentada en indicaciones psiquiátricas."
    strength: moderate
    model: human
evidence_against:
  - finding: "El vínculo con BUB1B es indirecto. La inestabilidad de MVA es de
              origen mitótico, por checkpoint débil, no principalmente
              replicativo."
  - finding: "Los antioxidantes acumulan un historial largo de fracasos al pasar
              de célula a clínica."
confidence: 0.15
falsifiable_by: >
  Linfocitos y fibroblastos del paciente con y sin NAC: frecuencia de micronúcleos
  y tasa de separación prematura de cromátidas. Si la PCS no baja, que es lo
  esperable si el defecto es puramente de checkpoint, el candidato queda
  descartado para este mecanismo.
```

---

## Descartados de forma explícita, y por qué

El filtro de realismo pediátrico no es decorativo. Estos son los candidatos que
un pipeline ingenuo propondría y que un panel experto rechazaría en la primera
lectura:

| Clase | Por qué se descarta |
|---|---|
| Inhibidores de Aurora A/B (alisertib y similares) | Aparecen en las búsquedas de aneuploidía porque se usan **contra** tumores con CIN. Administrárselos a un niño cuyas células ya missegregan empuja en la dirección del daño, no en contra. Sin ventana terapéutica defendible. |
| Inhibidores de MPS1/TTK | Peor aún: TTK es parte del mismo checkpoint que ya está debilitado. Inhibirlo agrava el defecto primario. |
| Antimitóticos citotóxicos en general | Toxicidad de quimioterapia en una indicación no oncológica, en un paciente con predisposición tumoral por daño genómico. |
| Letalidad sintética con CIN | Es una estrategia para matar células tumorales. En una enfermedad constitucional, las células diana son las del propio niño. |

### El filtro, con datos: qué existe contra el checkpoint

La tabla anterior es un argumento. Esto es la comprobación. Se consultaron DGIdb y
Open Targets para los 32 genes del panel, solo con símbolos e identificadores de
gen (`scripts/drug_landscape.py`, resultado en `drug_landscape.tsv`).

| Qué se midió | Resultado |
|---|---|
| Interacciones fármaco-gen en DGIdb, todo el panel | **123 inhibidores, 1 vacuna, 0 activadores** |
| Genes con candidatos en fase clínica (Open Targets) | 5 de 32: AURKB (18), PLK1 (8), PLK4 (3), TTK (2), BIRC5 (2) |
| Candidatos clínicos contra **BUB1B** | **0**, y 0 interacciones en DGIdb |

Todo lo que la farmacología conoce sobre esta vía **inhibe** el checkpoint. Un
paciente cuyo problema es tener demasiado poco checkpoint no tiene nada que ganar
ahí: no hay ningún activador, aprobado o experimental, de ninguno de los 32 genes.

Eso convierte en obligada una decisión que parecía de estilo. Los tres candidatos de
arriba actúan **aguas abajo**, sobre las consecuencias celulares de la aneuploidía,
porque aguas arriba, sobre la vía misma, no hay nada que empuje en la dirección
correcta.

Una advertencia de lectura: DGIdb también lista 65 fármacos aprobados asociados a
algún gen del panel, pero son en su mayoría citotóxicos oncológicos (cisplatino,
doxorrubicina, docetaxel, citarabina) y asociaciones de evidencia débil extraídas
de la literatura. El dato robusto es la dirección de las interacciones, no ese
recuento.

---

## Lectura honesta del conjunto

Las tres confianzas son bajas a propósito: 0.35, 0.25 y 0.15. Ninguna llega a la
mitad porque **ninguna tiene un solo dato en MVA**, y decir lo contrario sería
inventar.

El candidato 1 es el más fuerte por una razón concreta y no por entusiasmo: es el
único donde la intervención se probó en un modelo del gen correcto. Su debilidad
también es concreta: la prueba se hizo con una herramienta genética, en un
organismo envejecido, no con un fármaco en uno en desarrollo.

**Pendiente y declarado, con la razón correcta.** Una versión anterior de este
texto decía que Robin, LINCS y DepMap quedaban fuera "por falta de credenciales".
Solo es cierto para Robin, que usa créditos de pago de Edison. DepMap y LINCS son
gratuitos:

- **DepMap** se usó al final, a través de Open Targets, pero para una pregunta
  distinta de la letalidad sintética: la esencialidad de los genes MVA, que está
  en `MECHANISM_BUB1B.md`. La letalidad sintética sigue descartada por la razón de
  la tabla de arriba, no por acceso.
- **LINCS/CMap no se corrió.** Un cruce de firmas transcripcionales contra el
  knockdown de BUB1B podría proponer candidatos que esta ruta, basada en
  literatura, no ve. Es el pendiente con más potencial de Track 2.
