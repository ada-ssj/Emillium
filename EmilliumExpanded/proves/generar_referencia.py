"""
Genera el conjunt de fotos FITS de referència d'Emillium Expanded (només numpy).

Cada foto és simulada amb una "veritat" coneguda (posicions, magnituds,
moviments) que es desa a referencia/referencia.json. Les proves de regressió
(prova_regressio.py) mesuren aquestes fotos i comparen amb la veritat.

Condicions comunes (també al JSON):
  - Estrelles gaussianes circulars; fons del cel constant; soroll de Poisson
    (amb el GAIN de la capçalera) + soroll de lectura.
  - Imatges de 16 bits (BITPIX=16, BZERO=32768) com les d'una càmera real.
  - Solució astromètrica TAN a la capçalera (excepte on s'indica).
  - Magnituds: flux = 10^(-0,4·(m − 20,5)) · temps d'exposició (en ADU): una estrella\n    de magnitud 11 té un pic d'uns 38.000 ADU (sense saturar).

Casos (taula del document de requisits):
  brillants        precisió del centre en estrelles brillants
  febles           robustesa davant del soroll (S/R 10–30)
  dobles           separació i contaminació (parelles a 3–8 px)
  variable_1..6    precisió fotomètrica (estrella variable + 3 comparacions)
  asteroide_1..4   posició i seguiment d'un objecte mòbil
  desenfocada      detecció d'imatges dolentes (FWHM 3 vegades més gran)
  camp_dens        identificació de fonts properes (400 estrelles)
  escala_fina/gruixuda  robustesa davant de diferents escales (0,6 i 3″/px)

Ús:  python proves/generar_referencia.py   (sobreescriu proves/referencia/)
"""
import json
import math
import os

import numpy as np

CARPETA = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'referencia')
LLAVOR = 20261108
JD0 = 2461110.40


def targeta(clau, valor, comentari=''):
    if isinstance(valor, bool):
        text = 'T' if valor else 'F'
        t = f'{clau:<8}= {text:>20}'
    elif isinstance(valor, (int, np.integer)):
        t = f'{clau:<8}= {int(valor):>20}'
    elif isinstance(valor, float):
        t = f'{clau:<8}= {valor:>20.12G}'
    else:
        t = f"{clau:<8}= '{str(valor):<8}'"
    if comentari:
        t += f' / {comentari}'
    return t[:80].ljust(80)


def escriu_fits(ruta, imatge, capcalera):
    """FITS de 16 bits sense signe (BZERO = 32768), escrit a mà."""
    dades = np.clip(np.round(imatge), 0, 65535).astype(np.int64) - 32768
    targetes = [targeta('SIMPLE', True), targeta('BITPIX', 16), targeta('NAXIS', 2),
                targeta('NAXIS1', imatge.shape[1]), targeta('NAXIS2', imatge.shape[0]),
                targeta('BZERO', 32768.0), targeta('BSCALE', 1.0)]
    targetes += [targeta(k, v) for k, v in capcalera.items()]
    targetes.append('END'.ljust(80))
    text = ''.join(targetes)
    text = text.ljust((len(text) + 2879) // 2880 * 2880)
    binari = dades.astype('>i2').tobytes()
    with open(ruta, 'wb') as fitxer:
        fitxer.write(text.encode('ascii') + binari + b'\0' * ((-len(binari)) % 2880))


def wcs(escala_arcsec, ample, alt, rotacio=12.0, centre=(150.12, 22.45)):
    """Capçalera TAN (CD) i funció píxel(0) -> cel per a la veritat."""
    e = escala_arcsec / 3600
    r = math.radians(rotacio)
    cd = [[-e * math.cos(r), e * math.sin(r)], [e * math.sin(r), e * math.cos(r)]]
    crpix = ((ample + 1) / 2, (alt + 1) / 2)
    cap = {'CTYPE1': 'RA---TAN', 'CTYPE2': 'DEC--TAN', 'CRPIX1': float(crpix[0]), 'CRPIX2': float(crpix[1]),
           'CRVAL1': float(centre[0]), 'CRVAL2': float(centre[1]), 'CD1_1': cd[0][0], 'CD1_2': cd[0][1],
           'CD2_1': cd[1][0], 'CD2_2': cd[1][1]}

    def a_cel(x, y):   # píxel amb origen 0 (Emillium) -> AR, Dec en graus
        u, v = x + 1 - crpix[0], y + 1 - crpix[1]
        xi = math.radians(cd[0][0] * u + cd[0][1] * v)
        eta = math.radians(cd[1][0] * u + cd[1][1] * v)
        a0, d0 = math.radians(centre[0]), math.radians(centre[1])
        den = math.cos(d0) - eta * math.sin(d0)
        ar = a0 + math.atan2(xi, den)
        dec = math.atan2(math.sin(d0) + eta * math.cos(d0), math.hypot(xi, den))
        return math.degrees(ar) % 360, math.degrees(dec)
    return cap, a_cel


def foto(ample, alt, estrelles, fwhm, generador, fons=800.0, guany=1.5, lectura=6.0, exposicio=60.0):
    """estrelles: [(x, y, magnitud)] amb el centre del píxel (0, 0) a (0, 0)."""
    s = fwhm / 2.3548
    imatge = np.full((alt, ample), fons, dtype=np.float64)
    yy, xx = np.mgrid[0:alt, 0:ample]
    for x, y, m in estrelles:
        flux = 10 ** (-0.4 * (m - 20.5)) * exposicio
        r = int(6 * s) + 2
        x0, x1, y0, y1 = max(int(x) - r, 0), min(int(x) + r + 2, ample), max(int(y) - r, 0), min(int(y) + r + 2, alt)
        # integral exacta de la gaussiana sobre cada píxel (amb la funció d'error)
        ex = _integral_1d(np.arange(x0, x1), x, s)
        ey = _integral_1d(np.arange(y0, y1), y, s)
        imatge[y0:y1, x0:x1] += flux * np.outer(ey, ex)
    imatge = generador.poisson(imatge * guany) / guany + generador.normal(0, lectura / guany, imatge.shape)
    return imatge


def _integral_1d(pixels, centre, s):
    from math import erf
    a = (pixels - 0.5 - centre) / (s * math.sqrt(2))
    b = (pixels + 0.5 - centre) / (s * math.sqrt(2))
    return np.array([(erf(bb) - erf(aa)) / 2 for aa, bb in zip(a, b)])


def capcalera_base(cap_wcs, k=0, exposicio=60.0, minuts=10.0, filtre='V'):
    from datetime import datetime, timedelta, timezone
    inici = datetime(2026, 3, 10, 21, 35, 30, tzinfo=timezone.utc) + timedelta(minutes=minuts * k)
    c = dict(cap_wcs)
    c.update({'DATE-OBS': inici.strftime('%Y-%m-%dT%H:%M:%S.000'), 'EXPTIME': exposicio, 'GAIN': 1.5,
              'SATURATE': 65000.0, 'FILTER': filtre, 'OBJECT': 'EMILLIUM-REF'})
    return c


def main():
    os.makedirs(CARPETA, exist_ok=True)
    g = np.random.default_rng(LLAVOR)
    veritat = {'llavor': LLAVOR, 'condicions': {'fons_adu': 800.0, 'guany_e_adu': 1.5, 'lectura_e': 6.0,
                                                'exposicio_s': 60.0, 'punt_zero': 20.5, 'bitpix': 16}, 'casos': {}}
    A, H = 240, 200

    def desa(nom, imatge, cap, dades):
        escriu_fits(os.path.join(CARPETA, nom + '.fits'), imatge, cap)
        veritat['casos'][nom] = dades

    # --- brillants i febles
    for nom, mags, fwhm in (('brillants', (11.0, 11.5, 12.0, 12.5), 3.0), ('febles', (15.5, 15.9, 16.3, 16.7), 3.0)):
        cap, a_cel = wcs(1.2, A, H)
        est = [(x + g.uniform(-0.5, 0.5), y + g.uniform(-0.5, 0.5), m)
               for (x, y), m in zip(((60, 50), (180, 60), (70, 150), (170, 140)), mags)]
        desa(nom, foto(A, H, est, fwhm, g), capcalera_base(cap),
             {'fwhm': fwhm, 'estrelles': [{'x': x, 'y': y, 'mag': m, 'ar': a_cel(x, y)[0], 'dec': a_cel(x, y)[1]}
                                          for x, y, m in est]})
    # --- dobles
    cap, a_cel = wcs(1.2, A, H)
    parelles = []
    for (x, y), sep, angle, dm in (((60, 60), 3.0, 30, 0.5), ((170, 60), 5.0, 120, 1.0), ((60, 150), 8.0, 250, 1.5),
                                   ((170, 150), 6.0, 340, 0.0)):
        x2, y2 = x + sep * math.cos(math.radians(angle)), y + sep * math.sin(math.radians(angle))
        parelles.append(((x + 0.3, y - 0.2, 12.0), (x2, y2, 12.0 + dm)))
    est = [e for p in parelles for e in p]
    desa('dobles', foto(A, H, est, 2.8, g), capcalera_base(cap),
         {'fwhm': 2.8, 'parelles': [{'a': {'x': a[0], 'y': a[1], 'mag': a[2], 'ar': a_cel(a[0], a[1])[0], 'dec': a_cel(a[0], a[1])[1]},
                                     'b': {'x': b[0], 'y': b[1], 'mag': b[2], 'ar': a_cel(b[0], b[1])[0], 'dec': a_cel(b[0], b[1])[1]}}
                                    for a, b in parelles]})
    # --- variable (6 fotos)
    cap, a_cel = wcs(1.2, A, H)
    base = {'V': (120.4, 100.3, 12.8), 'C1': (50.2, 40.7, 12.0), 'C2': (190.6, 50.1, 12.6), 'C3': (60.3, 160.2, 13.1)}
    for k in range(6):
        delta = 0.40 * math.sin(2 * math.pi * k / 6)
        est = [(x, y, m + (delta if nom == 'V' else 0)) for nom, (x, y, m) in base.items()]
        desa(f'variable_{k + 1}', foto(A, H, est, 3.0, g), capcalera_base(cap, k),
             {'fwhm': 3.0, 'estrelles': {nom: {'x': x, 'y': y, 'mag': m} for nom, (x, y, m) in
                                         zip(base, est)}, 'jd_mig': JD0 + 10 / 1440 * k})
    # --- asteroide (4 fotos): 3 px per foto
    cap, a_cel = wcs(1.2, A, H)
    fixes = [(40.2, 40.7, 12.5), (200.4, 60.3, 12.8), (110.6, 170.2, 13.0)]
    for k in range(4):
        x, y = 70.3 + 3.0 * k, 80.6 + 1.2 * k
        desa(f'asteroide_{k + 1}', foto(A, H, fixes + [(x, y, 14.0)], 3.0, g), capcalera_base(cap, k),
             {'fwhm': 3.0, 'fixes': [{'x': a, 'y': b} for a, b, _m in fixes],
              'asteroide': {'x': x, 'y': y, 'ar': a_cel(x, y)[0], 'dec': a_cel(x, y)[1]}, 'jd_mig': JD0 + 10 / 1440 * k})
    # --- desenfocada (comparar amb brillants)
    cap, a_cel = wcs(1.2, A, H)
    est = [(60.2, 50.4, 11.0), (180.7, 60.1, 11.5), (70.4, 150.6, 12.0), (170.1, 140.3, 12.5)]
    desa('desenfocada', foto(A, H, est, 9.0, g), capcalera_base(cap), {'fwhm': 9.0, 'estrelles': [{'x': x, 'y': y} for x, y, _ in est]})
    # --- camp dens
    cap, a_cel = wcs(1.2, A, H)
    est = [(g.uniform(8, A - 8), g.uniform(8, H - 8), g.uniform(12.5, 16.0)) for _ in range(400)]
    desa('camp_dens', foto(A, H, est, 2.6, g), capcalera_base(cap),
         {'fwhm': 2.6, 'estrelles': [{'x': x, 'y': y, 'mag': m} for x, y, m in est]})
    # --- escales diferents
    for nom, escala, fwhm in (('escala_fina', 0.6, 5.5), ('escala_gruixuda', 3.0, 1.6)):
        cap, a_cel = wcs(escala, A, H, rotacio=-25.0)
        est = [(x + 0.37, y + 0.81, m) for (x, y), m in zip(((50, 50), (180, 50), (50, 150), (180, 150), (120, 100)),
                                                           (12.5, 13.0, 13.5, 14.0, 13.2))]
        desa(nom, foto(A, H, est, fwhm, g), capcalera_base(cap),
             {'fwhm': fwhm, 'escala': escala, 'estrelles': [{'x': x, 'y': y, 'mag': m, 'ar': a_cel(x, y)[0],
                                                              'dec': a_cel(x, y)[1]} for x, y, m in est]})
    with open(os.path.join(CARPETA, 'referencia.json'), 'w', encoding='utf-8') as fitxer:
        json.dump(veritat, fitxer, indent=1, ensure_ascii=False)
    print(f"{len(veritat['casos'])} fotos de referència a {CARPETA}")


if __name__ == '__main__':
    main()
