###### Emillium v1.3 - Medidor de centroides estelares ######
###### Hecho por Plutarco Martínez y Emilio Ruiz con asistencia de Claude Code ##########

Programa de escritorio en Python para abrir imágenes astronómicas **FITS**, girarlas
o invertirlas como un espejo, y medir la posición de las estrellas **en píxeles con
precisión subpíxel** (centroide). Sirve para comparar posiciones entre distintas fotos.

No usa WCS ni bibliotecas de astronomía (astropy, photutils…). El lector FITS, el giro
de la imagen y el cálculo del centroide están programados en `Emillium.py` y comentados
paso a paso. Solo necesita Python (con tkinter, que ya viene incluido) y **numpy**.

## Contenido de esta carpeta

- `Emillium.exe`: el programa listo para usar (doble clic; no necesita Python).
- `codigo\`: el código fuente (`Emillium.py`), el icono y lo necesario para crear el .exe.

Las versiones anteriores están en `..\Emillium1.2` y `..\Emillium1.1`; la 1.0 y la 0.9 (con WCS), en `..\Emillium`.

## Novedades de la 1.3

- **Medidas más precisas en estrellas débiles**: el «Umbral sobre el ruido» (antes
  «Sensibilidad») decide qué píxeles cuentan para el centro, así el ruido del cielo ya
  no lo desvía (en estrellas débiles, de ~0,8 px de error a ~0,03 px con el valor 3).
  Al cambiar una opción avanzada, la medida se recalcula sola.
- **Estrellas juntas**: se mide la estrella más cercana al clic (no la más brillante) y
  se quita la luz de las vecinas brillantes. Si el centro coincide con una estrella ya
  medida con otra letra, sale un aviso y no se guarda.
- **Giro con ángulo personalizado** con signo bien visible (+ antihorario, − horario) y
  botones −10° / +10°, en el paso «Modifica la foto».
- **Las estrellas medidas siguen marcadas al girar** la foto o ponerla en espejo
  (círculo magenta discontinuo).
- **Botones nuevos**: «Cerrar la foto» (vuelve al inicio sin borrar medidas), «Borrar la
  medida» (aparece al seleccionar filas; también tecla Supr) y «Limpiar medidas»
  (pequeño y rojo, pide confirmación).
- **Panel derecho desplazable** con la rueda del ratón o una barra, si no cabe.
- **Easter eggs nuevos** (ver la tabla de combinaciones secretas), entre ellos `2001`
  con música sintetizada, y Pilares de la Creación rehechos, con vaciado y relleno del cielo.

## Novedades de la 1.2

- **Idioma**: catalán por defecto; se puede cambiar a castellano con el desplegable
  **Idioma** de arriba a la derecha (no se pierde la foto ni las medidas).
- **Estrellas con letra** (A, B, C…): antes de hacer clic se elige qué estrella se mide.
  Por defecto, después de medir pasa sola a la letra siguiente y, al abrir otra foto,
  vuelve a la A. Si se mide otra vez la misma letra en la misma foto, sustituye a la anterior.
- **Tabla** reducida a: N.º de foto, foto, estrella, orientación, X, Y y flujo.
- **Excel estrella por estrella**: una fila por foto y, para cada estrella, sus columnas
  `A_X (px)`, `A_Y (px)`, `A_Flux`, `B_X (px)`…

## Usarlo desde el código

Desde la carpeta `codigo`:

```bash
pip install -r requirements.txt
python Emillium.py
```

## Cómo se usa

La ventana tiene un panel a la derecha con los pasos numerados. Si dejas el ratón
encima de un botón, aparece una explicación.

| Acción | Cómo |
|---|---|
| Abrir foto | **Abrir foto FITS…** (o Ctrl+O). Se lee la primera imagen 2D del archivo. |
| Elegir la estrella | Desplegable **Estrella** (A, B, C…). Usa la misma letra para la misma estrella en todas las fotos. |
| Medir una estrella | **Clic izquierdo** cerca de la estrella. El resultado sale en grande a la derecha. |
| Cambiar de idioma | Desplegable **Idioma** (Català / Castellano). |
| Zoom | Rueda del ratón, o los botones **−**, **+** y **Ver foto entera** de arriba. |
| Mover la foto | Arrastrar con el **botón derecho** (o la rueda pulsada). |
| Girar | En **2. Modifica la foto**, escribir el ángulo con su signo (`+12,5`, `-7`…) y pulsar Intro o **Aplicar el giro**. **+** = sentido antihorario, **−** = sentido horario. Los botones **−10°** y **+10°** lo cambian de 10 en 10. El ángulo es siempre respecto a la foto original. |
| Espejo | **Horizontal** (izquierda ↔ derecha) o **Vertical** (arriba ↔ abajo). |
| Volver al original | **Volver a la foto original**. |
| Ver una medición | Doble clic en su fila de la tabla. |
| Guardar | **Guardar para Excel…** (Ctrl+S): guarda toda la tabla en `.xlsx` (o `.csv`). |
| Borrar una medida | Seleccionarla en la tabla: aparece el botón **Borrar la medida** (también tecla Supr). Se pueden seleccionar varias con Ctrl o Mayús. Ctrl+Z borra la última. |
| Empezar de cero | Botón pequeño y rojo **Limpiar medidas**, debajo de la tabla a la derecha: borra todas las medidas (pide confirmación) y vuelve a numerar las fotos desde el 1. |
| Quitar la foto | **Cerrar la foto** (debajo del nombre de la foto): vuelve a la pantalla de inicio sin borrar las medidas. |

Opciones avanzadas (botón **Mostrar opciones avanzadas**):

- **Distancia máx. al clic**: hasta qué distancia del clic se busca la estrella.
- **Radio para el centro**: radio del círculo que se usa para calcular el centroide.
  Conviene que sea unas 2–3 veces la FWHM de las estrellas.
- **Umbral sobre el ruido (σ)** (antes "Sensibilidad"): cuántas veces el ruido del cielo
  tiene que superar un píxel para contar al calcular el centro (y para detectar la estrella).
  Recomendado: 3. Con valores bajos, el ruido del cielo desvía el centro de las estrellas
  débiles; con estrellas simuladas débiles el error pasa de ~0,4 px (umbral 0) a ~0,05–0,08 px
  (umbral 3).
- Al cambiar cualquiera de estas opciones, la medida seleccionada en la tabla (o, si no hay
  ninguna, la última) se vuelve a calcular sola con los valores nuevos.

El desplegable **Brillo** solo cambia cómo se ve la foto, no afecta a las medidas.

## Combinaciones secretas (easter eggs)

Se teclean seguidas (menos de 1,2 s entre tecla y tecla), sin estar escribiendo en una casilla.

| Teclas | Qué pasa | Dónde |
|---|---|---|
| `M42` | La Nebulosa de Orión cruza el cielo | Pantalla de inicio |
| `618` | TON 618: agujero negro que hace orbitar las estrellas | Pantalla de inicio |
| `M1` o `CRAB` | Supernova: fogonazo, onda de luz de lado a lado, Nebulosa del Cangrejo y un púlsar que se queda parpadeando en la esquina de la cabecera | Pantalla de inicio |
| `M16` | Cortina de gas con los Pilares de la Creación; nacen estrellas nuevas que se quedan en el cielo | Pantalla de inicio |
| `1977` | La sonda Voyager cruza despacio dejando ondas de radio, con sonido de espacio profundo | Pantalla de inicio |
| `ISS` | La Estación Espacial cruza el cielo | Pantalla de inicio |
| `2001` | Una nave en homenaje a «2001: una odisea del espacio» cruza despacio con el comienzo de «Así habló Zaratustra» (R. Strauss, 1896, dominio público), sintetizado por el propio programa y a volumen bajo | Pantalla de inicio |
| `EARTH` | Zoom hacia atrás desde la Tierra hasta el "punto azul pálido" de la Voyager 1 | Siempre |
| `6EQUJ5` o `WOW` | Monitor de radiotelescopio con la señal Wow! | Siempre |
| `SPOCK` o `VULCAN` | Saludo vulcano | Siempre |
| `GARGANTUA` o `INTERSTELLAR` | Lente gravitatoria durante 3,6 s (deforma las estrellas o la foto abierta, no las medidas) | Siempre |
| `42` | Pensamiento Profundo calcula la respuesta definitiva | Siempre |

Las de "Siempre" tapan la foto un momento; un clic las cierra antes de tiempo.

## Coordenadas

- X crece hacia la derecha e Y hacia arriba. El centro del píxel de abajo a la
  izquierda es (0, 0).
- Las coordenadas se dan en la imagen **tal como se ve** (después del giro y el espejo).
  En la tabla y en el Excel se guarda también la orientación de cada medición, descrita
  (por ejemplo "Espejo vertical, girada +25° (antihorario)"). Primero se aplican los
  espejos y después el giro.
- Los giros múltiplos de 90° son exactos. Los demás ángulos interpolan los píxeles
  (bilineal); las esquinas vacías quedan sin datos.
- Las diferencias y distancias entre mediciones se calculan en Excel a partir de X e Y.

## Método del centroide

1. Se buscan los picos (píxeles más brillantes que sus 8 vecinos) que destacan sobre el
   ruido cerca del clic y se elige **el más cercano al clic**, no el más brillante. Así una
   estrella brillante al lado no "roba" el clic a una débil.
   Si hay una estrella brillante muy cerca, se le quita su luz aprovechando que es
   simétrica (su luz en un punto es igual a la del punto opuesto respecto a su centro),
   y el círculo para el centro se hace más pequeño.
   Si el centro calculado coincide (a menos de 1,5 px) con una estrella ya medida con
   otra letra en la misma foto, sale un aviso y la medida no se guarda.
2. Se mide el fondo del cielo en un anillo alrededor de la estrella (mediana con
   recorte a 3σ).
3. Se comprueba que el pico supera el umbral en σ.
4. Se calcula el centro de masas de la luz (brillo − fondo − umbral·σ, solo los píxeles
   que superan el umbral) dentro de la apertura, y se repite recentrando la apertura hasta
   que se mueve menos de 0,0001 px.
5. Se dan también el flujo, la FWHM aproximada y el pico/σ.

Con estrellas simuladas, el error del centroide fue de entre 0,003 y 0,025 px.

## Exportar a Excel

- **`.xlsx` (recomendado)**: libro de Excel con dos hojas:
  - **Coordenades píxels**: una fila por foto con `N.º de foto`, `Foto`, `Orientació` y,
    para cada estrella, `A_X (px)`, `A_Y (px)`, `A_Flux`, `B_X (px)`… Si alguna foto
    tiene la fecha de la toma en la cabecera FITS (`DATE-OBS`), sale también la columna
    `Data i hora`. Si una estrella no se midió en una foto, sus casillas quedan vacías.
  - **Totes les mesures**: la misma información en forma de lista (una fila por medida).
  Los números se guardan con todos sus decimales (en pantalla se ven con 3), así que
  se puede calcular con ellos directamente. Funciona con Excel en cualquier idioma,
  LibreOffice y Google Sheets. La fila de títulos es azul oscuro, queda fija al bajar
  y tiene los botones de filtro.
- **`.csv`**: la tabla por estrellas, con `;` como separador y coma decimal. Excel lo
  abre bien solo si está configurado en catalán o español.
- Los títulos de las columnas salen en el idioma elegido en el programa.
- El número de foto es el orden en que se abrió cada foto por primera vez.

## Crear el .exe

```bash
pip install pyinstaller
pyinstaller Emillium.spec
```

(desde la carpeta `codigo`; el .exe nuevo aparece en `codigo\dist\`).


