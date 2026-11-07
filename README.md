###### Emillium v1.3 - Mesurador de centroides estel·lars ######
###### Fet per Plutarco Martínez i Emilio Ruiz amb l'assistència de Claude Code ##########

**Emillium** és un programa d'escriptori, lleuger i en català i castellà, per mesurar la
posició de les estrelles en fotos astronòmiques **FITS** amb **precisió de centèsimes de
píxel**: obres una foto (o una sèrie), fas clic a les estrelles i ho desa tot en un Excel.

Pensat per a astrònoms aficionats i per a treballs escolars (com un *treball de recerca*):
moviments propis, estrelles dobles, satèl·lits de Júpiter, desplaçament d'asteroides i
cometes, corbes de llum aproximades...

- **Lleuger i ràpid**: un sol `.exe` de ~25 MB, sense instal·lar res. Només necessita Python
  i `numpy` si es fa servir des del codi.
- **Sense internet**: no es connecta enlloc ni envia cap dada.
- **Codi obert i explicat**: el lector FITS, el gir, el centroide, les coordenades celestes i
  l'exportació a Excel estan programats a `codigo/Emillium.py`, comentats pas a pas.

---

## Què fa

| | |
|---|---|
| **Centroide subpíxel** | Centre de massa iteratiu amb llindar sobre el soroll (error de ~0,01 px en estrelles brillants i ~0,05 px en estrelles febles, mesurat amb estrelles simulades). Agafa l'estrella **més propera al clic**, i si en té una de brillant enganxada li treu la llum. |
| **Incertesa de cada mesura** | ± en X i Y (propagació d'errors amb el soroll del cel i, si el FITS porta `GAIN`, el soroll de la pròpia estrella). Comprovada amb simulacions de Monte Carlo. |
| **Coordenades celestes** | Si la foto té **solució astromètrica** (per exemple d'[Astrometry.net](https://nova.astrometry.net)), mostra **AR i Dec** de cada estrella. Projecció TAN, amb o sense distorsió **SIP**; coincideix amb astropy a menys de 0,01″. |
| **Temps** | Hora UTC del **mig de l'exposició** i **dia julià (JD)** a partir de `DATE-OBS` (+ `TIME-OBS`) i `EXPTIME`. |
| **Fotometria** | Flux d'obertura i **magnitud instrumental** amb el seu error (−2,5·log₁₀ del flux per segon + 25). |
| **Sèries de fotos** | Obre-ne moltes de cop, passa d'una a l'altra amb ◀ ▶ (o RePàg/AvPàg) mantenint gir i zoom, i **mesura automàticament les mateixes estrelles que a la foto anterior** (calcula sol quant s'ha desplaçat el camp). |
| **Girar i mirall** | Angle amb signe (+ antihorari, − horari), botons −10°/+10°, mirall horitzontal i vertical. Les estrelles ja mesurades segueixen marcades. |
| **Excel** | `.xlsx` de veritat (sense dependre de l'idioma d'Excel): una fila per foto i columnes per estrella (A_X, A_Y, ±, flux, magnitud, AR, Dec), més un full amb totes les mesures. També CSV. |
| **Català i castellà** | Català per defecte. Recorda l'idioma, les opcions i l'última carpeta. |

### Què NO fa (encara)

- No resol l'astrometria: per tenir AR/Dec, la foto ha de venir ja resolta (Astrometry.net, ASTAP, Siril...).
- No llegeix FITS comprimits (`.fz`): cal descomprimir-los abans (per exemple amb `funpack`).
- No calibra (biaix, foscos, flats) ni apila fotos.
- La magnitud és **instrumental**: per passar-la a magnitud real cal restar-la d'una estrella de magnitud coneguda (fotometria diferencial, fàcil de fer a l'Excel).
- No genera informes per al Minor Planet Center.

---

## Descàrrega i ús

**Amb el `.exe` (Windows):** descarrega `Emillium.exe` de l'apartat *Releases* i fes-hi doble clic.
La primera vegada Windows pot mostrar «Windows ha protegit l'equip» perquè el programa no està
signat: *Més informació → Executa igualment*.

**Des del codi (Windows, Linux o macOS):**

```bash
pip install -r codigo/requirements.txt
python codigo/Emillium.py
```

També pots arrossegar una o diverses fotos sobre el `.exe`.

## Guia ràpida

1. **Obrir foto FITS…** (Ctrl+O). Pots triar-ne diverses alhora per treballar amb una sèrie.
2. *(Opcional)* Gira-la o posa-la en mirall.
3. Tria la **lletra** de l'estrella (A, B, C...) i fes-hi **clic**. Fes servir la mateixa lletra
   per a la mateixa estrella a totes les fotos.
4. Passa a la foto següent (▶) i prem **«Mesurar com a la foto anterior»**.
5. **Desar per a Excel…** (Ctrl+S).

| Ratolí / tecla | |
|---|---|
| Clic esquerre | Mesurar l'estrella |
| Roda | Zoom (a sobre del panell dret, el desplaça) |
| Arrossegar amb el botó dret | Moure la foto |
| RePàg / AvPàg | Foto anterior / següent de la sèrie |
| Supr | Esborrar les mesures seleccionades a la taula |
| Ctrl+Z | Esborrar l'última mesura |

### Opcions avançades

- **Distància màx. al clic (px)**: fins on es busca l'estrella.
- **Radi per al centre (px)**: unes 2–3 vegades la FWHM de les estrelles.
- **Llindar sobre el soroll (σ)**: quins píxels compten per al centre. **Recomanat: 3.**
  Més baix, el soroll desvia el centre de les estrelles febles; molt alt (5 o més), en
  estrelles febles gairebé només compta el píxel central.

En canviar una opció, la mesura seleccionada (o l'última) es torna a calcular sola.

---

## Com es calcula (resum)

1. **Pic més proper al clic** entre els que destaquen sobre el soroll.
2. **Fons del cel** en un anell al voltant (mediana amb retall a 3σ).
3. Si hi ha una **veïna brillant**, se li treu la llum aprofitant que és simètrica.
4. **Centre de massa** dels píxels que superen el fons + llindar·σ, recentrant fins que es mou
   menys de 0,0001 px.
5. **Incerteses**: σ_X² = Σ (x_i − X)²·σ_i² / (Σ w_i)², amb σ_i² = σ_cel² (+ senyal/guany).
6. **AR/Dec**: píxel → (correcció SIP) → matriu CD → projecció gnomònica inversa.
7. **JD**: algorisme de Meeus (*Astronomical Algorithms*).

Coordenades en píxels: X cap a la dreta i Y cap amunt; el centre del píxel de baix a l'esquerra
és (0, 0) i es donen sobre la foto **tal com es veu** (després del gir i el mirall). Les AR/Dec
sempre es calculen sobre la foto original, així que no depenen del gir.

## Proves

```bash
python proves/prova_precisio.py
```

Crea estrelles simulades en posicions conegudes i comprova el lector FITS, el centroide, el gir
i les coordenades. Durant el desenvolupament també es va comparar amb astropy (AR/Dec: diferència
màxima de 0,004″ en 3.000 punts; JD: 0,5 ms en 500 dates).

## Seguretat i privacitat

- No fa servir internet ni executa programes externs.
- Només llegeix les fotos que obres i escriu on tu dius (Excel/CSV), més els ajustos a
  `%APPDATA%\Emillium\ajustes.json` i dos sons a la carpeta temporal.
- Els FITS manipulats (mides negatives o enormes) es rebutgen; els textos que van al CSV es
  protegeixen contra la «injecció de fórmules».

## Combinacions secretes

Hi ha uns quants *easter eggs* astronòmics: prova d'escriure `M42`, `618`, `M1`, `M16`, `1977`,
`2001`, `ISS`, `EARTH`, `WOW`, `SPOCK`, `GARGANTUA` o `42`...

## Crear el `.exe`

```bash
pip install pyinstaller
cd codigo
pyinstaller Emillium.spec
```

## Llicència

[MIT](LICENSE): pots fer-lo servir, modificar-lo i compartir-lo lliurement, mantenint l'avís de
llicència.

---

### Resum en castellano

Emillium mide la posición de las estrellas en fotos FITS con precisión subpíxel, con su
incertidumbre, magnitud instrumental, AR/Dec (si la foto tiene solución astrométrica TAN/SIP),
hora UTC y día juliano del centro de la exposición. Trabaja con series de fotos (incluso mide
automáticamente las mismas estrellas que en la foto anterior) y exporta a Excel. Funciona sin
internet, está en catalán y castellano (desplegable «Idioma») y el código está comentado en
castellano paso a paso.
