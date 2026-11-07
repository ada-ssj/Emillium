<a name="catala"></a>

<p align="center">
  <img src="docs/logo.png" width="112" alt="Logo d'Emillium">
</p>

<h1 align="center">Emillium</h1>

<p align="center">
  <b>Mesurador de centroides estel·lars per a fotos astronòmiques FITS</b><br>
  <sub>Complete Edition · v1.4</sub>
</p>

<p align="center">
  <img src="https://img.shields.io/badge/Complete-1.4-ffd27a" alt="Complete Edition 1.4">
  <img src="https://img.shields.io/badge/Lite-1.3-86c8ff" alt="Edició Lite 1.3">
  <img src="https://img.shields.io/badge/Python-3.9%2B-3776AB?logo=python&logoColor=white" alt="Python 3.9+">
  <img src="https://img.shields.io/badge/depend%C3%A8ncies-numpy-013243?logo=numpy&logoColor=white" alt="Dependències: numpy">
  <img src="https://img.shields.io/badge/Windows-.exe-0078D6" alt="Windows .exe">
  <img src="https://img.shields.io/badge/llic%C3%A8ncia-MIT-2ea44f" alt="Llicència MIT">
  <img src="https://img.shields.io/badge/idiomes-catal%C3%A0%20%7C%20castellano-8a96bb" alt="Idiomes: català i castellano">
</p>

<p align="center">
  <b>Català</b> · <a href="#castellano">Castellano</a>
</p>

<p align="center">
  <img src="docs/mesures.png" width="880" alt="Emillium mesurant un asteroide en una sèrie de fotos">
</p>

**Emillium** és un programa d'escriptori, lleuger en català i castellà, per mesurar la posició de
les estrelles en fotos **FITS** amb **precisió de centèsimes de píxel**: obres una foto (o una sèrie de fotos),
fas clic a les estrelles i pots exportar tot en un Excel, amb la seva incertesa, magnitud i coordenades.

Està pensat per a **astrònoms aficionats**, per a **treballs escolars** o **activitats educatives**:
moviments propis, estrelles dobles, satèl·lits de Júpiter, seguiment d'asteroides i cometes,
estrelles variables...

> **Lleuger** — un sol `.exe` de ~25 MB, sense instal·lar res.<br>
> **Privat** — funciona sense internet i no envia cap dada.<br>
> **Obert i explicat** — tot (lector FITS, centroide, coordenades, Excel...) està programat a
> [`codigo/Emillium.py`](codigo/Emillium.py) i comentat pas a pas.

---

## Índex

- [Dues edicions: Complete i Lite](#dues-edicions-complete-i-lite)
- [Funcions](#funcions)
- [Captures](#captures)
- [Instal·lació](#installació)
- [Guia ràpida](#guia-ràpida)
- [Precisió i validació](#precisió-i-validació)
- [Limitacions](#limitacions)
- [Com es calcula](#com-es-calcula)
- [Estructura del projecte](#estructura-del-projecte)
- [Autors i llicència](#autors-i-llicència)

---

## Dues edicions: Complete i Lite

Emillium es distribueix en dues edicions. Les dues mesuren les estrelles amb el mateix centroide
subpíxel; la diferència és tot el que hi ha al voltant.

- **Complete Edition (v1.4)** — per a astrònoms aficionats: incerteses, AR/Dec, temps, fotometria,
  sèries de fotos, fitxers `.fz` i informes per al MPC. És la que es descriu en aquest README.
- **Edició Lite (v1.3)** — la versió senzilla, pensada per a activitats educatives: obres una foto, fas clic
  a les estrelles i obtens les seves coordenades **X i Y en píxels** i el flux. Menys columnes i menys
  opcions, per no distreure.

| | Lite 1.3 | Complete 1.4 |
|---|:---:|:---:|
| Centroide subpíxel (estrella més propera al clic, veïnes brillants, avís d'estrella repetida) | ✓ | ✓ |
| Gir amb angle i signe, mirall, marques en girar | ✓ | ✓ |
| Estrelles amb lletra i Excel «estrella per estrella» (X, Y, flux) | ✓ | ✓ |
| Català i castellà, *easter eggs* | ✓ | ✓ |
| Incertesa ± de cada mesura | | ✓ |
| Magnitud instrumental i fotometria diferencial | | ✓ |
| AR/Dec (solució astromètrica TAN/SIP) | | ✓ |
| Hora UTC del mig de l'exposició i dia julià | | ✓ |
| Sèries de fotos: ◀ ▶ i «Mesurar com a la foto anterior» | | ✓ |
| FITS comprimits (`.fz`) | | ✓ |
| Informe per al Minor Planet Center | | ✓ |
| Recorda l'idioma, les opcions i l'última carpeta | | ✓ |
| Protecció contra fitxers FITS manipulats | | ✓ |

**Quina triar?** Si només necessites posicions en píxels per a un treball o una pràctica, la **Lite**.
Si fas astrometria o fotometria (asteroides, estrelles variables, sèries de fotos), la **Complete**.

---

## Funcions

Funcions de la **Complete Edition**:

| | |
|---|---|
| **Centroide subpíxel** | Centre de massa iteratiu amb llindar sobre el soroll. Agafa l'estrella **més propera al clic** i, si en té una de brillant enganxada, li treu la llum. |
| **Incertesa** | ± en X, Y i en el flux de cada mesura (soroll del cel i, si el FITS porta `GAIN`, de la pròpia estrella). |
| **Coordenades celestes** | **AR i Dec** de cada estrella si la foto té solució astromètrica (per exemple d'[Astrometry.net](https://nova.astrometry.net)): projecció TAN, amb o sense distorsió **SIP**. |
| **Temps** | Hora UTC del **mig de l'exposició** i **dia julià (JD)**, a partir de `DATE-OBS` (+ `TIME-OBS`) i `EXPTIME`. |
| **Fotometria** | Flux d'obertura i **magnitud instrumental**; amb una **estrella de comparació** de magnitud coneguda, **magnitud real** de totes les altres (fotometria diferencial). |
| **Sèries de fotos** | Obre-ne moltes de cop, passa d'una a l'altra amb ◀ ▶ mantenint gir i zoom, i **mesura automàticament les mateixes estrelles que a la foto anterior**. |
| **FITS comprimits** | Obre directament fitxers **`.fz`** (RICE_1, GZIP_1, GZIP_2, també en coma flotant quantitzada), sense programes externs. |
| **Informe per al MPC** | Assistent que prepara les observacions d'un asteroide o cometa en el **format de 80 columnes del Minor Planet Center**, amb vista prèvia. |
| **Gir i mirall** | Angle amb signe, botons −10° / +10°, mirall horitzontal i vertical. Les estrelles ja mesurades segueixen marcades. |
| **Excel** | `.xlsx` real (no depèn de l'idioma d'Excel): una fila per foto i columnes per estrella, més un full amb totes les mesures. També CSV. |
| **Comoditat** | Català per defecte i castellà; recorda l'idioma, les opcions i l'última carpeta. |

## Captures

<table>
  <tr>
    <td width="58%"><img src="docs/inici.png" alt="Pantalla d'inici"></td>
    <td><img src="docs/informe_mpc.png" alt="Assistent de l'informe per al MPC"></td>
  </tr>
  <tr>
    <td align="center"><sub>Pantalla d'inici</sub></td>
    <td align="center"><sub>Assistent de l'informe per al MPC</sub></td>
  </tr>
</table>

---

## Instal·lació

### Opció 1 · Executable per a Windows (recomanat)

1. A l'apartat [**Releases**](../../releases), descarrega l'executable de l'edició que vulguis
   (**Complete** o **Lite**).
2. Fes-hi doble clic. No cal instal·lar res. Les dues edicions es poden tenir alhora.

> [!NOTE]
> Com que el programa no està signat digitalment, la primera vegada Windows pot mostrar
> «Windows ha protegit l'equip». Prem **Més informació → Executa igualment**.

També pots arrossegar una o diverses fotos FITS sobre el `.exe` per obrir-les directament.

### Opció 2 · Des del codi

Necessita **Python 3.9 o superior** (provat amb Python 3.14), amb `tkinter` (ja ve inclòs a Windows)
i `numpy`:

```bash
pip install -r codigo/requirements.txt
python codigo/Emillium.py
```

Està provat a Windows. A Linux i macOS hauria de funcionar des del codi, però encara no s'ha provat.

---

## Guia ràpida

1. **Obrir foto FITS…** (`Ctrl+O`). Pots triar-ne diverses alhora per treballar amb una sèrie.
2. *(Opcional)* Gira-la o posa-la en mirall.
3. Tria la **lletra** de l'estrella (A, B, C...) i fes-hi **clic**. Fes servir la mateixa lletra per a la
   mateixa estrella a totes les fotos.
4. Passa a la foto següent (▶) i prem **Mesurar com a la foto anterior**.
5. *(Opcional)* A **Fotometria diferencial**, tria l'estrella de comparació i escriu la seva magnitud.
6. **Desar per a Excel…** (`Ctrl+S`) o, per a un asteroide, **Informe MPC…**

<details>
<summary><b>Ratolí i dreceres de teclat</b></summary>

| Acció | Ratolí / tecla |
|---|---|
| Mesurar una estrella | Clic esquerre |
| Zoom | Roda del ratolí (a sobre del panell dret, el desplaça) |
| Moure la foto | Arrossegar amb el botó dret |
| Foto anterior / següent | `RePàg` / `AvPàg` |
| Esborrar les mesures seleccionades | `Supr` |
| Esborrar l'última mesura | `Ctrl+Z` |
| Obrir / desar | `Ctrl+O` / `Ctrl+S` |

</details>

<details>
<summary><b>Opcions avançades</b></summary>

- **Distància màx. al clic (px)**: fins on es busca l'estrella.
- **Radi per al centre (px)**: unes 2–3 vegades la FWHM de les estrelles.
- **Llindar sobre el soroll (σ)**: quins píxels compten per al centre. **Recomanat: 3.**
  Més baix, el soroll desvia el centre de les estrelles febles; molt alt (5 o més), en estrelles
  febles gairebé només compta el píxel central.

En canviar una opció, la mesura seleccionada (o l'última) es torna a calcular sola.

</details>

**Coordenades en píxels:** X cap a la dreta i Y cap amunt; el centre del píxel de baix a l'esquerra és
(0, 0). Es donen sobre la foto **tal com es veu** (després del gir i el mirall). Les AR/Dec sempre es
calculen sobre la foto original, així que no depenen del gir.

---

## Precisió i validació

Tot s'ha comprovat amb estrelles simulades en posicions conegudes i comparant amb
[astropy](https://www.astropy.org), la biblioteca de referència en astronomia:

| Què | Resultat |
|---|---|
| Centre d'estrelles brillants / febles | error de ~0,003–0,01 px / ~0,05 px |
| Incertesa ± predita | coincideix amb la dispersió real (simulació de Monte Carlo) amb un ~10–15 % |
| AR/Dec (TAN, CDELT+CROTA i SIP) | diferència màxima de **0,004″** amb astropy en 3.000 punts |
| Dia julià | diferència màxima de **0,5 ms** amb astropy en 500 dates |
| Fitxers `.fz` | resultat **idèntic** a astropy en 20 combinacions de compressió, tipus de dades i quantització |
| Fotometria diferencial | magnitud real d'una estrella simulada recuperada amb un error de 0,02 mag |
| Informe MPC | línies de 80 columnes; posicions a 0,02–0,11″ de les reals |

Per repetir les proves bàsiques (només necessita `numpy`):

```bash
python proves/prova_precisio.py
```

---

## Limitacions

**Fotos i coordenades**

- **No resol l'astrometria**: per tenir AR/Dec, la foto ha de venir ja resolta (Astrometry.net, ASTAP, Siril...).
- Les AR/Dec només surten amb solucions de projecció **TAN** (amb o sense SIP). Amb altres projeccions
  (TPV, ZPN...) la foto s'obre i es mesura en píxels, però sense coordenades celestes.
- Considera que `DATE-OBS` és l'**inici** de l'exposició (i hi suma `EXPTIME`/2). Si la teva càmera hi
  desa el mig de l'exposició, l'hora quedarà desplaçada mitja exposició.
- Llegeix només la **primera imatge** del fitxer FITS; si és una imatge en color (diversos plans), els
  fa la mitjana.
- No llegeix els formats comprimits poc habituals **HCOMPRESS_1** i **PLIO_1** (avisa; cal fer servir `funpack`).
- No calibra (biaix, foscos, flats) ni apila fotos.

**Mesures i resultats**

- Les incerteses ± són **estadístiques** (soroll): no inclouen errors sistemàtics ni l'error de la
  solució astromètrica de la foto.
- La fotometria diferencial fa servir **una** estrella de comparació, sense correcció de color ni d'extinció.
- L'informe del MPC és en el format clàssic de 80 columnes, que el MPC encara accepta, però des de 2018
  el format recomanat és **ADES**, i Emillium no el genera.

**Sistema**

- Només s'ha provat a **Windows**.

> [!IMPORTANT]
> Abans d'enviar cap informe al Minor Planet Center, revisa'l i consulta les instruccions vigents
> a la [web del MPC](https://minorplanetcenter.net). Comprova també que la solució astromètrica de
> les fotos sigui bona.

---

## Com es calcula

<details>
<summary><b>Resum del mètode</b> (tot està explicat també al codi)</summary>

1. **Pic més proper al clic** entre els que destaquen sobre el soroll.
2. **Fons del cel** en un anell al voltant de l'estrella (mediana amb retall a 3σ).
3. Si hi ha una **veïna brillant**, se li treu la llum aprofitant que és simètrica.
4. **Centre de massa** dels píxels que superen el fons + llindar·σ, recentrant fins que es mou menys
   de 0,0001 px.
5. **Incerteses**: σ<sub>X</sub>² = Σ (x<sub>i</sub> − X)² · σ<sub>i</sub>² / (Σ w<sub>i</sub>)², amb
   σ<sub>i</sub>² = σ<sub>cel</sub>² (+ senyal / guany).
6. **AR/Dec**: píxel → correcció SIP → matriu CD → projecció gnomònica inversa.
7. **Dia julià**: algorisme de Meeus (*Astronomical Algorithms*).
8. **Fotometria diferencial**: m = m<sub>inst</sub> − m<sub>inst, comp</sub> + m<sub>catàleg, comp</sub>,
   amb error √(σ² + σ<sub>comp</sub>²).
9. **`.fz`**: lectura de la taula binària de «teseles» i descompressió RICE (decodificada en paral·lel
   per a totes les files amb numpy) o gzip; a les imatges en coma flotant es desfà la quantització amb
   la taula de nombres aleatoris estàndard del FITS.
10. **MPC**: data UTC del mig de l'exposició amb 5 decimals de dia (~0,9 s), AR a 0,01 s i Dec a 0,1″.

</details>

<details>
<summary><b>Seguretat i privacitat</b></summary>

- No fa servir internet ni executa programes externs.
- Només llegeix les fotos que obres i escriu on tu dius (Excel, CSV, informe), més els ajustos a
  `%APPDATA%\Emillium\ajustes.json` i dos sons a la carpeta temporal.
- Els FITS manipulats (mides negatives o enormes) es rebutgen, i els textos que van al CSV es
  protegeixen contra la «injecció de fórmules».

</details>

<details>
<summary><b>Combinacions secretes</b></summary>

Hi ha uns quants *easter eggs* astronòmics. Prova d'escriure a la pantalla d'inici `M42`, `618`, `M1`,
`M16`, `1977`, `2001`, `ISS`, `EARTH`, `WOW`, `SPOCK`, `GARGANTUA` o `42`...

</details>

---

## Estructura del projecte

```
Emillium/
├── codigo/
│   ├── Emillium.py          # tot el programa, comentat
│   ├── Emillium.spec        # configuració per crear el .exe
│   ├── requirements.txt     # numpy
│   └── icono.ico
├── proves/
│   └── prova_precisio.py    # prova de precisió amb estrelles simulades
├── docs/                    # imatges d'aquest README
├── LICENSE
└── README.md
```

### Crear l'executable

```bash
pip install pyinstaller
cd codigo
pyinstaller Emillium.spec
```

L'`Emillium.exe` apareix a `codigo/dist/`.

---

## Autors i llicència

Fet per **Plutarco Martínez** i **Emilio Ruiz**.

Distribuït sota la [llicència MIT](LICENSE): el pots fer servir, modificar i compartir lliurement,
mantenint l'avís de llicència.

---

<br>

<a name="castellano"></a>

<p align="center">
  <img src="docs/logo.png" width="112" alt="Logo de Emillium">
</p>

<h1 align="center">Emillium</h1>

<p align="center">
  <b>Medidor de centroides estelares para fotos astronómicas FITS</b><br>
  <sub>Complete Edition · v1.4</sub>
</p>

<p align="center">
  <img src="https://img.shields.io/badge/Complete-1.4-ffd27a" alt="Complete Edition 1.4">
  <img src="https://img.shields.io/badge/Lite-1.3-86c8ff" alt="Edición Lite 1.3">
  <img src="https://img.shields.io/badge/Python-3.9%2B-3776AB?logo=python&logoColor=white" alt="Python 3.9+">
  <img src="https://img.shields.io/badge/dependencias-numpy-013243?logo=numpy&logoColor=white" alt="Dependencias: numpy">
  <img src="https://img.shields.io/badge/Windows-.exe-0078D6" alt="Windows .exe">
  <img src="https://img.shields.io/badge/licencia-MIT-2ea44f" alt="Licencia MIT">
  <img src="https://img.shields.io/badge/idiomas-catal%C3%A0%20%7C%20castellano-8a96bb" alt="Idiomas: catalán y castellano">
</p>

<p align="center">
  <a href="#catala">Català</a> · <b>Castellano</b>
</p>

<p align="center">
  <img src="docs/mesures.png" width="880" alt="Emillium midiendo un asteroide en una serie de fotos">
</p>

**Emillium** es un programa de escritorio, ligero en catalán y castellano, para medir la posición de
las estrellas en fotos **FITS** con **precisión de centésimas de píxel**: abres una foto (o una serie de fotos),
haces clic en las estrellas y lo exportas todo en un Excel, con su incertidumbre, magnitud y coordenadas.

Está pensado para **astrónomos aficionados**, para **trabajos escolares** o **actividades educativas**:
movimientos propios, estrellas dobles, satélites de Júpiter, seguimiento de asteroides y cometas,
estrellas variables...

> **Ligero** — un solo `.exe` de ~25 MB, sin instalar nada.<br>
> **Privado** — funciona sin internet y no envía ningún dato.<br>
> **Abierto y explicado** — todo (lector FITS, centroide, coordenadas, Excel...) está programado a mano en
> [`codigo/Emillium.py`](codigo/Emillium.py) y comentado paso a paso.

---

## Índice

- [Dos ediciones: Complete y Lite](#dos-ediciones-complete-y-lite)
- [Funciones](#funciones)
- [Capturas](#capturas)
- [Instalación](#instalación)
- [Guía rápida](#guía-rápida)
- [Precisión y validación](#precisión-y-validación)
- [Limitaciones](#limitaciones)
- [Cómo se calcula](#cómo-se-calcula)
- [Estructura del proyecto](#estructura-del-proyecto)
- [Autores y licencia](#autores-y-licencia)

---

## Dos ediciones: Complete y Lite

Emillium se divide en dos ediciones. Las dos miden las estrellas con el mismo centroide
subpíxel; la diferencia es todo lo que hay alrededor.

- **Complete Edition (v1.4)** — para astrónomos aficionados: incertidumbres, AR/Dec, tiempo, fotometría,
  series de fotos, archivos `.fz` e informes para el MPC. Es la que se describe en este README.
- **Edición Lite (v1.3)** — la versión sencilla, pensada para actividades educativas: abres una foto, haces clic
  en las estrellas y obtienes sus coordenadas **X e Y en píxeles** y el flujo. Menos columnas y menos
  opciones, para no distraer.

| | Lite 1.3 | Complete 1.4 |
|---|:---:|:---:|
| Centroide subpíxel (estrella más cercana al clic, vecinas brillantes, aviso de estrella repetida) | ✓ | ✓ |
| Giro con ángulo y signo, espejo, marcas al girar | ✓ | ✓ |
| Estrellas con letra y Excel «estrella por estrella» (X, Y, flujo) | ✓ | ✓ |
| Catalán y castellano, *easter eggs* | ✓ | ✓ |
| Incertidumbre ± de cada medida | | ✓ |
| Magnitud instrumental y fotometría diferencial | | ✓ |
| AR/Dec (solución astrométrica TAN/SIP) | | ✓ |
| Hora UTC del centro de la exposición y día juliano | | ✓ |
| Series de fotos: ◀ ▶ y «Medir como en la foto anterior» | | ✓ |
| FITS comprimidos (`.fz`) | | ✓ |
| Informe para el Minor Planet Center | | ✓ |
| Recuerda el idioma, las opciones y la última carpeta | | ✓ |
| Protección contra archivos FITS manipulados | | ✓ |

**¿Cuál elegir?** Si solo necesitas posiciones en píxeles para un trabajo o una práctica, la **Lite**.
Si haces astrometría o fotometría (asteroides, estrellas variables, series de fotos), la **Complete**.

---

## Funciones

Funciones de la **Complete Edition**:

| | |
|---|---|
| **Centroide subpíxel** | Centro de masas iterativo con umbral sobre el ruido. Coge la estrella **más cercana al clic** y, si tiene una brillante pegada, le quita su luz. |
| **Incertidumbre** | ± en X, Y y en el flujo de cada medida (ruido del cielo y, si el FITS lleva `GAIN`, de la propia estrella). |
| **Coordenadas celestes** | **AR y Dec** de cada estrella si la foto tiene solución astrométrica (por ejemplo de [Astrometry.net](https://nova.astrometry.net)): proyección TAN, con o sin distorsión **SIP**. |
| **Tiempo** | Hora UTC del **centro de la exposición** y **día juliano (JD)**, a partir de `DATE-OBS` (+ `TIME-OBS`) y `EXPTIME`. |
| **Fotometría** | Flujo de apertura y **magnitud instrumental**; con una **estrella de comparación** de magnitud conocida, **magnitud real** de todas las demás (fotometría diferencial). |
| **Series de fotos** | Abre muchas a la vez, pasa de una a otra con ◀ ▶ manteniendo giro y zoom, y **mide automáticamente las mismas estrellas que en la foto anterior**. |
| **FITS comprimidos** | Abre directamente archivos **`.fz`** (RICE_1, GZIP_1, GZIP_2, también en coma flotante cuantizada), sin programas externos. |
| **Informe para el MPC** | Asistente que prepara las observaciones de un asteroide o cometa en el **formato de 80 columnas del Minor Planet Center**, con vista previa. |
| **Giro y espejo** | Ángulo con signo, botones −10° / +10°, espejo horizontal y vertical. Las estrellas ya medidas siguen marcadas. |
| **Excel** | `.xlsx` real (no depende del idioma de Excel): una fila por foto y columnas por estrella, más una hoja con todas las medidas. También CSV. |
| **Comodidad** | Catalán por defecto y castellano; recuerda el idioma, las opciones y la última carpeta. |

## Capturas

<table>
  <tr>
    <td width="58%"><img src="docs/inici.png" alt="Pantalla de inicio"></td>
    <td><img src="docs/informe_mpc.png" alt="Asistente del informe para el MPC"></td>
  </tr>
  <tr>
    <td align="center"><sub>Pantalla de inicio</sub></td>
    <td align="center"><sub>Asistente del informe para el MPC</sub></td>
  </tr>
</table>

<sub>La interfaz se muestra en catalán; se cambia a castellano con el desplegable «Idioma».</sub>

---

## Instalación

### Opción 1 · Ejecutable para Windows (recomendado)

1. En el apartado [**Releases**](../../releases), descarga el ejecutable de la edición que quieras
   (**Complete** o **Lite**).
2. Haz doble clic. No hace falta instalar nada. Puedes tener las dos ediciones a la vez.

> [!NOTE]
> Como el programa no está firmado digitalmente, la primera vez Windows puede mostrar
> «Windows protegió su PC». Pulsa **Más información → Ejecutar de todas formas**.

También puedes arrastrar una o varias fotos FITS sobre el `.exe` para abrirlas directamente.

### Opción 2 · Desde el código

Necesita **Python 3.9 o superior** (probado con Python 3.14), con `tkinter` (ya viene incluido en Windows)
y `numpy`:

```bash
pip install -r codigo/requirements.txt
python codigo/Emillium.py
```

Está probado en Windows. En Linux y macOS debería funcionar desde el código, pero todavía no se ha probado.

---

## Guía rápida

1. **Abrir foto FITS…** (`Ctrl+O`). Puedes elegir varias a la vez para trabajar con una serie.
2. *(Opcional)* Gírala o ponla en espejo.
3. Elige la **letra** de la estrella (A, B, C...) y haz **clic** en ella. Usa la misma letra para la
   misma estrella en todas las fotos.
4. Pasa a la foto siguiente (▶) y pulsa **Medir como en la foto anterior**.
5. *(Opcional)* En **Fotometría diferencial**, elige la estrella de comparación y escribe su magnitud.
6. **Guardar para Excel…** (`Ctrl+S`) o, para un asteroide, **Informe MPC…**

<details>
<summary><b>Ratón y atajos de teclado</b></summary>

| Acción | Ratón / tecla |
|---|---|
| Medir una estrella | Clic izquierdo |
| Zoom | Rueda del ratón (encima del panel derecho, lo desplaza) |
| Mover la foto | Arrastrar con el botón derecho |
| Foto anterior / siguiente | `RePág` / `AvPág` |
| Borrar las medidas seleccionadas | `Supr` |
| Borrar la última medida | `Ctrl+Z` |
| Abrir / guardar | `Ctrl+O` / `Ctrl+S` |

</details>

<details>
<summary><b>Opciones avanzadas</b></summary>

- **Distancia máx. al clic (px)**: hasta dónde se busca la estrella.
- **Radio para el centro (px)**: unas 2–3 veces la FWHM de las estrellas.
- **Umbral sobre el ruido (σ)**: qué píxeles cuentan para el centro. **Recomendado: 3.**
  Más bajo, el ruido desvía el centro de las estrellas débiles; muy alto (5 o más), en estrellas
  débiles casi solo cuenta el píxel central.

Al cambiar una opción, la medida seleccionada (o la última) se vuelve a calcular sola.

</details>

**Coordenadas en píxeles:** X hacia la derecha e Y hacia arriba; el centro del píxel de abajo a la izquierda es
(0, 0). Se dan sobre la foto **tal como se ve** (después del giro y el espejo). Las AR/Dec siempre se
calculan sobre la foto original, así que no dependen del giro.

---

## Precisión y validación

Todo se ha comprobado con estrellas simuladas en posiciones conocidas y comparando con
[astropy](https://www.astropy.org), la biblioteca de referencia en astronomía:

| Qué | Resultado |
|---|---|
| Centro de estrellas brillantes / débiles | error de ~0,003–0,01 px / ~0,05 px |
| Incertidumbre ± predicha | coincide con la dispersión real (simulación de Monte Carlo) con un ~10–15 % |
| AR/Dec (TAN, CDELT+CROTA y SIP) | diferencia máxima de **0,004″** con astropy en 3.000 puntos |
| Día juliano | diferencia máxima de **0,5 ms** con astropy en 500 fechas |
| Archivos `.fz` | resultado **idéntico** a astropy en 20 combinaciones de compresión, tipo de datos y cuantización |
| Fotometría diferencial | magnitud real de una estrella simulada recuperada con un error de 0,02 mag |
| Informe MPC | líneas de 80 columnas; posiciones a 0,02–0,11″ de las reales |

Para repetir las pruebas básicas (solo necesita `numpy`):

```bash
python proves/prova_precisio.py
```

---

## Limitaciones

**Fotos y coordenadas**

- **No resuelve la astrometría**: para tener AR/Dec, la foto tiene que venir ya resuelta (Astrometry.net, ASTAP, Siril...).
- Las AR/Dec solo salen con soluciones de proyección **TAN** (con o sin SIP). Con otras proyecciones
  (TPV, ZPN...) la foto se abre y se mide en píxeles, pero sin coordenadas celestes.
- Considera que `DATE-OBS` es el **inicio** de la exposición (y le suma `EXPTIME`/2). Si tu cámara guarda
  ahí el centro de la exposición, la hora quedará desplazada media exposición.
- Lee solo la **primera imagen** del archivo FITS; si es una imagen en color (varios planos), hace la
  media.
- No lee los formatos comprimidos poco habituales **HCOMPRESS_1** y **PLIO_1** (avisa; hay que usar `funpack`).
- No calibra (bias, darks, flats) ni apila fotos.

**Medidas y resultados**

- Las incertidumbres ± son **estadísticas** (ruido): no incluyen errores sistemáticos ni el error de la
  solución astrométrica de la foto.
- La fotometría diferencial usa **una** estrella de comparación, sin corrección de color ni de extinción.
- El informe del MPC está en el formato clásico de 80 columnas, que el MPC todavía acepta, pero desde 2018
  el formato recomendado es **ADES**, y Emillium no lo genera.

**Sistema**

- Solo se ha probado en **Windows**.

> [!IMPORTANT]
> Antes de enviar ningún informe al Minor Planet Center, revísalo y consulta las instrucciones vigentes
> en la [web del MPC](https://minorplanetcenter.net). Comprueba también que la solución astrométrica de
> las fotos sea buena.

---

## Cómo se calcula

<details>
<summary><b>Resumen del método</b> (todo está explicado también en el código)</summary>

1. **Pico más cercano al clic** entre los que destacan sobre el ruido.
2. **Fondo del cielo** en un anillo alrededor de la estrella (mediana con recorte a 3σ).
3. Si hay una **vecina brillante**, se le quita la luz aprovechando que es simétrica.
4. **Centro de masas** de los píxeles que superan el fondo + umbral·σ, recentrando hasta que se mueve menos
   de 0,0001 px.
5. **Incertidumbres**: σ<sub>X</sub>² = Σ (x<sub>i</sub> − X)² · σ<sub>i</sub>² / (Σ w<sub>i</sub>)², con
   σ<sub>i</sub>² = σ<sub>cielo</sub>² (+ señal / ganancia).
6. **AR/Dec**: píxel → corrección SIP → matriz CD → proyección gnomónica inversa.
7. **Día juliano**: algoritmo de Meeus (*Astronomical Algorithms*).
8. **Fotometría diferencial**: m = m<sub>inst</sub> − m<sub>inst, comp</sub> + m<sub>catálogo, comp</sub>,
   con error √(σ² + σ<sub>comp</sub>²).
9. **`.fz`**: lectura de la tabla binaria de «teselas» y descompresión RICE (decodificada en paralelo
   para todas las filas con numpy) o gzip; en las imágenes en coma flotante se deshace la cuantización con
   la tabla de números aleatorios estándar del FITS.
10. **MPC**: fecha UTC del centro de la exposición con 5 decimales de día (~0,9 s), AR a 0,01 s y Dec a 0,1″.

</details>

<details>
<summary><b>Seguridad y privacidad</b></summary>

- No usa internet ni ejecuta programas externos.
- Solo lee las fotos que abres y escribe donde tú dices (Excel, CSV, informe), además de los ajustes en
  `%APPDATA%\Emillium\ajustes.json` y dos sonidos en la carpeta temporal.
- Los FITS manipulados (tamaños negativos o enormes) se rechazan, y los textos que van al CSV se
  protegen contra la «inyección de fórmulas».

</details>

<details>
<summary><b>Combinaciones secretas</b></summary>

Hay unos cuantos *easter eggs* astronómicos. Prueba a escribir en la pantalla de inicio `M42`, `618`, `M1`,
`M16`, `1977`, `2001`, `ISS`, `EARTH`, `WOW`, `SPOCK`, `GARGANTUA` o `42`...

</details>

---

## Estructura del proyecto

```
Emillium/
├── codigo/
│   ├── Emillium.py          # todo el programa, comentado
│   ├── Emillium.spec        # configuración para crear el .exe
│   ├── requirements.txt     # numpy
│   └── icono.ico
├── proves/
│   └── prova_precisio.py    # prueba de precisión con estrellas simuladas
├── docs/                    # imágenes de este README
├── LICENSE
└── README.md
```

### Crear el ejecutable

```bash
pip install pyinstaller
cd codigo
pyinstaller Emillium.spec
```

El `Emillium.exe` aparece en `codigo/dist/`.

---

## Autores y licencia

Hecho por **Plutarco Martínez** y **Emilio Ruiz**.

Distribuido bajo la [licencia MIT](LICENSE): puedes usarlo, modificarlo y compartirlo libremente,
manteniendo el aviso de licencia.
