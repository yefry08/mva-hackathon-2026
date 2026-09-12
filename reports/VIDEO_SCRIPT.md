# Guion del pitch — 3 minutos

Ritmo objetivo: unas 145 palabras por minuto. Total ~430 palabras, que deja aire
para respirar y para que las imágenes asienten.

---

## [0:00 – 0:25] El problema, y lo que no vamos a prometer

> Hay un niño con aneuploidía variegada en mosaico. Menos de cincuenta personas
> en el mundo tienen este diagnóstico. Su familia hizo públicos sus datos para
> que alguien mirara.
>
> No vamos a decir que lo curamos. Vamos a enseñar qué encontramos, cómo lo
> comprobamos, y en qué puntos exactos el método se queda corto.

*Imagen: el título del hackathon. Nada de stock de laboratorios.*

## [0:25 – 1:05] El hallazgo

> Buscamos bajo modelo recesivo en los treinta y dos genes del checkpoint
> mitótico. Cero homocigotos raros. Pero **un solo gen del panel tiene dos
> variantes raras en región codificante: BUB1B**, el gen de MVA tipo 1.
>
> Una de las dos ya está en ClinVar como patogénica para esta enfermedad
> exacta: introduce un codón de parada, y verificamos contra la estructura
> exónica que el transcrito va a degradación. La otra es un missense presente en
> **un alelo entre 1,46 millones**.
>
> Un alelo nulo y un alelo hipomórfico: es la arquitectura que la literatura
> describe para esta enfermedad desde 2004.

*Imagen: la tabla de los 32 genes, con BUB1B como única fila con dos codificantes.*

## [1:05 – 1:35] Y sale sin que se lo digamos

> Un pipeline que solo encuentra la respuesta cuando le dices dónde mirar no
> demuestra nada. Así que lo corrimos ciego: **todo el genoma contra ClinVar, sin
> prior de genes**. Siete variantes patogénicas en total, ninguna homocigota. Y
> BUB1B es la única que además tiene un segundo alelo raro en el mismo gen.

*Imagen: las siete, con BUB1B destacada.*

## [1:35 – 2:10] Lo que no podemos afirmar

> Aquí es donde muchos pipelines callan. Nosotros lo medimos.
>
> Para que las dos variantes causen la enfermedad tienen que estar en cromosomas
> distintos. Los bloques de fase de estos datos llegan como máximo a **206 pares
> de bases**. Las dos variantes están a **diez mil novecientas once**. Cincuenta y
> tres veces más lejos. Ninguna lectura de este experimento puede cubrir ambas:
> no es un fallo del análisis, la información no está en los datos.
>
> Y el segundo alelo, por criterios ACMG, es una variante de significado
> incierto. Lo seguirá siendo aunque se resuelva la fase.
>
> Lo que lo resolvería es barato: genotipar a los padres.

*Imagen: la distribución de tamaños de bloque de fase, con la distancia marcada.*

## [2:10 – 2:40] Que no inventa diagnósticos

> Corrimos el mismo pipeline, sin tocar un parámetro, sobre HG002: un genoma
> público de un adulto sano. Resultado: **cero candidatos**.
>
> Y algo que solo se ve con el control: esa persona sana carga **cuarenta y cinco
> variantes raras** en genes del checkpoint mitótico. Una más que el paciente. A
> nivel de variante rara, los dos genomas son indistinguibles. Toda la
> discriminación vive en la anotación funcional, no en el umbral de frecuencia.

*Imagen: la tabla comparativa de dos columnas.*

## [2:40 – 3:00] Cierre

> Todo esto corre en un portátil, con Python de la biblioteca estándar, sin nube,
> y sin que los datos del niño salgan de la máquina: un hook bloquea cualquier
> coordenada que intente salir, y lo probamos con treinta canarios.
>
> Coste de cómputo: cero.
>
> Para Track 2 proponemos tres hipótesis de reposicionamiento, cada una con el
> experimento que la refutaría. La más fuerte llega a una confianza de 0.35,
> porque ninguna tiene todavía un solo dato en esta enfermedad. Decirlo también
> es parte del trabajo.

*Imagen: el checklist de privacidad en verde y la línea de coste.*
