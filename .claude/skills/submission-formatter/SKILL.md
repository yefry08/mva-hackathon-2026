---
name: submission-formatter
description: Genera y valida los archivos de submission para Track 1 y Track 2. Usar antes de cualquier entrega al leaderboard.
---

# Formateo de submissions

## Antes de generar cualquier archivo
Lee las specs vigentes en la pagina del hackathon (pestanas "Submit - Track 1" y
"Submit - Track 2"). El formato exacto manda sobre lo que diga esta skill.

## Track 1
- Lista rankeada de variantes (rank 1 = mas confiable)
- Coordenadas en el build correcto, con REF/ALT explicitos
- Score de confianza por variante — F-max depende de la calibracion, no del volumen
- Pares heterocigotos compuestos: ambas variantes presentes
- Hallazgos secundarios permitidos (no penalizan el score automatico)
- Link a GitHub + methods writeup

## Track 2
- Reporte estructurado contra la rubrica: Rigor 35%, Impacto 25%, Innovacion 25%,
  Escalabilidad 15%
- La seccion de Escalabilidad debe mostrar el pipeline corriendo sobre OTRO caso
  no diagnosticado publico. Es 15% del score y casi nadie lo hace.
- Link a GitHub + video de 3 minutos

## Checklist de privacidad (bloqueante)
- [ ] Ningun FASTQ/BAM/VCF/CRAM en el repo
- [ ] Ninguna coordenada cruda fuera de la lista de variantes entregada
- [ ] Ningun identificador de muestra ni nombre de archivo original del paciente
- [ ] `git log --all --stat | grep -Ei '\.(bam|vcf|fastq|cram)'` sin resultados
- [ ] Licencia CC BY 4.0 declarada

## Presupuesto
Track 1: 6 submissions (1 sanity + 3 iteracion + 2 reserva)
Track 2: 3 submissions (1 borrador interno + 2 finales)

**Nunca envies sin confirmacion explicita del humano.**
