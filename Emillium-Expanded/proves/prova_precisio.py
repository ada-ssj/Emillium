"""
Prova de precisió d'Emillium (només necessita numpy).

Crea fotos FITS amb estrelles simulades en posicions conegudes i comprova:
  - el lector FITS (escrivint els FITS a mà),
  - el centroide i la seva incertesa,
  - el gir i el mirall (i la seva inversa),
  - les coordenades celestes (AR/Dec) i el dia julià.

Ús:  python proves/prova_precisio.py
"""
import math
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'codigo'))
import Emillium as E  # noqa: E402

generador = np.random.default_rng(2026)
errors = []


def comprova(condicio, text):
    print(('  OK   ' if condicio else '  ERROR ') + text)
    if not condicio:
        errors.append(text)


def escriu_fits(ruta, imatge, capcalera=()):
    """FITS mínim de 32 bits en coma flotant, escrit a mà."""
    targetes = ['SIMPLE  =                    T', 'BITPIX  =                  -32', 'NAXIS   =                    2',
                f'NAXIS1  = {imatge.shape[1]:20d}', f'NAXIS2  = {imatge.shape[0]:20d}', *capcalera, 'END']
    text = ''.join(t.ljust(80) for t in targetes)
    text = text.ljust((len(text) + 2879) // 2880 * 2880)
    dades = imatge.astype('>f4').tobytes()
    with open(ruta, 'wb') as fitxer:
        fitxer.write(text.encode('ascii') + dades + b'\0' * ((-len(dades)) % 2880))


# --- 1. Estrelles simulades i lector FITS
alt, ample = 300, 400
yy, xx = np.mgrid[0:alt, 0:ample]
estrelles = [(100.37, 80.81, 5000, 1.6), (250.12, 200.55, 1500, 2.0), (330.9, 40.2, 400, 1.3)]
imatge = np.full((alt, ample), 1000.0)
for x, y, amp, s in estrelles:
    imatge += amp * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2 * s * s))
imatge += generador.normal(0, 10, imatge.shape)
ruta = os.path.join(tempfile.gettempdir(), 'emillium_prova.fits')
escriu_fits(ruta, imatge, ["DATE-OBS= '2000-01-01T11:59:30'", 'EXPTIME =                 60.0',
                           "CTYPE1  = 'RA---TAN'", "CTYPE2  = 'DEC--TAN'",
                           'CRVAL1  =                 83.8', 'CRVAL2  =                 -5.4',
                           'CRPIX1  =                200.5', 'CRPIX2  =                150.5',
                           'CD1_1   =             -0.00030', 'CD1_2   =                  0.0',
                           'CD2_1   =                  0.0', 'CD2_2   =              0.00030'])
llegida, capcalera, _nota = E.leer_fits(ruta)
print('Lector FITS')
comprova(np.allclose(llegida, imatge.astype(np.float32)), 'les dades llegides coincideixen amb les escrites')

# --- 2. Centroide i incertesa
print('Centroide')
for x, y, amp, s in estrelles:
    r = E.medir_estrella(llegida, x + 2.5, y - 1.5, 15, 3 * s, 3)
    error = math.hypot(r['x'] - x, r['y'] - y)
    comprova(error < 0.1, f'estrella a ({x}, {y}): error {error:.4f} px (± predita {max(r["sigma_x"], r["sigma_y"]):.4f})')

# --- 3. Gir i mirall
print('Gir i mirall')
for angle, mirall in ((90, False), (-37, True), (12.5, False)):
    girada = E.transformar_imagen(llegida, angle, mirall)
    x, y = estrelles[0][:2]
    gx, gy = E.transformar_punto(x, y, ample, alt, angle, mirall)
    r = E.medir_estrella(girada, gx, gy, 6, 5, 3)
    tornada = E.destransformar_punto(r['x'], r['y'], ample, alt, angle, mirall)
    comprova(math.hypot(tornada[0] - x, tornada[1] - y) < 0.1, f'angle {angle}°, mirall {mirall}: torna a la posició original')

# --- 4. Coordenades celestes i temps
print('Coordenades i temps')
wcs = E.SolucionWCS.desde_cabecera(capcalera)
ar, dec = wcs.a_cielo(200.5, 150.5)
comprova(abs(ar - 83.8) < 1e-9 and abs(dec + 5.4) < 1e-9, 'al píxel de referència (CRPIX) dona exactament CRVAL')
ar2, _ = wcs.a_cielo(210.5, 150.5)
separacio = (ar - ar2) * math.cos(math.radians(dec)) * 3600   # segons d'arc
comprova(abs(separacio - 10 * 0.00030 * 3600) < 0.01,
         f'10 píxels cap a la dreta = 10 × 1,08″ = {separacio:.3f}″ (cap a l\'oest: AR més petita)')
_text, jd = E.instante_observacion(capcalera)
comprova(abs(jd - 2451545.0) < 1e-8, f'J2000.0 (2000-01-01 12:00 UTC, mig de l\'exposició) = JD 2451545,0 (dona {jd})')

# --- 5. Fotometria i informe del MPC
print('Fotometria i MPC')
m1, _ = E.magnitud_instrumental(5000.0, 10.0, 60)
m2, _ = E.magnitud_instrumental(2000.0, 10.0, 60)
comprova(abs((m1 - m2) + 2.5 * math.log10(5000 / 2000)) < 1e-12, 'diferència de magnituds = −2,5·log10(F1/F2)')
linia = E.linea_mpc('EML0001', True, 'C', 2451545.0, 83.633212, 22.014501, 12.34, 'R', 'J95')
comprova(len(linia) == 80, 'la línia del MPC té exactament 80 columnes')
comprova(linia[5:12] == 'EML0001' and linia[12] == '*' and linia[14] == 'C', 'designació, descobriment i tipus al seu lloc')
comprova(linia[15:32] == '2000 01 01.50000 ', 'data a les columnes 16-32')
comprova(linia[32:44] == '05 34 31.97 ' and linia[44:56] == '+22 00 52.2 ', 'AR i Dec a les columnes 33-56')
comprova(linia[65:71] == '12.3 R' and linia[77:80] == 'J95', 'magnitud, banda i codi d\'observatori')
comprova(E.ar_mpc(359.99999999) == '00 00 00.00', 'els segons de AR s\'arrodoneixen bé (59,999… -> 00)')

print()
print('TOT CORRECTE' if not errors else f'{len(errors)} PROVES HAN FALLAT')
sys.exit(1 if errors else 0)
