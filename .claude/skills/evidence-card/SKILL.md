---
name: evidence-card
description: Formato estandar de ficha de evidencia con cita obligatoria. Usar en toda tarea de literatura, mecanismo o reposicionamiento de farmacos.
---

# Ficha de evidencia

Toda afirmacion biologica que entre a un reporte usa este formato:

```yaml
claim: "<una afirmacion, una sola>"
evidence_for:
  - source: "PMID:12345678"   # o DOI
    finding: "<que mostro ese trabajo>"
    strength: strong|moderate|weak
    model: human|mouse|cell_line|in_silico
evidence_against:
  - source: "PMID:..."
    finding: "..."
confidence: 0.0-1.0
falsifiable_by: "<experimento concreto que refutaria esto>"
```

## Reglas
- Sin DOI/PMID verificable → marca `[UNVERIFIED]` y no entra al reporte final.
- `evidence_against` vacio es sospechoso. Si de verdad no hay contraevidencia, escribe
  "no encontrada tras busqueda en X" y di donde buscaste.
- Nunca cites un paper que no leiste. El adversarial-reviewer verifica citas al azar
  y una cita fabricada invalida el output completo del worker.
- Extrapolar de linea celular a paciente pediatrico es un salto; nombralo cuando lo hagas.

## Para candidatos de farmaco, ademas
- Estado regulatorio y indicacion aprobada actual
- Perfil de seguridad conocido **en poblacion pediatrica**
- Ventana terapeutica argumentada
- Sin `falsifiable_by` concreto, el candidato no entra
