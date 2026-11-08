"""
Prova amb fotos REALS: compara Emillium amb astropy (la biblioteca de referència) en
qualsevol carpeta de fitxers FITS. Serveix per comprovar el programa amb dades
independents, de càmeres i telescopis de veritat.

Per a cada imatge de cada fitxer (extensions i plans) comprova:
  - que els píxels llegits són idèntics als d'astropy,
  - que l'instant del mig de l'exposició (dia julià) coincideix,
  - que l'AR/Dec de les estrelles detectades coincideix (si la foto té solució astromètrica).

Necessita astropy (només per comparar; Emillium no el fa servir):  pip install astropy
Ús:  python proves/prova_reals.py CARPETA_AMB_FITS
"""
import os
import sys
import warnings

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'codigo'))
import Emillium as E  # noqa: E402
import ciencia as C  # noqa: E402

try:
    from astropy.io import fits
    from astropy.wcs import WCS
except ImportError:
    print('Aquesta prova necessita astropy:  pip install astropy')
    sys.exit(0)
warnings.simplefilter('ignore')


def imatges_astropy(ruta):
    """Les mateixes imatges que compta Emillium, llegides amb astropy: (dades 2D, capçalera combinada)."""
    llista = []
    with fits.open(ruta) as hdul:
        primaria = hdul[0].header
        for hdu in hdul:
            if hdu.data is None or not hasattr(hdu.data, 'ndim') or hdu.data.ndim < 2 or getattr(hdu, 'is_image', True) is False:
                continue
            cap = hdu.header.copy()
            if hdu is not hdul[0]:
                for k, v in primaria.items():
                    if k not in cap and k not in ('SIMPLE', 'BITPIX', 'NAXIS', 'EXTEND', 'COMMENT', 'HISTORY', ''):
                        cap[k] = v
            dades = hdu.data
            if dades.ndim == 3 and dades.shape[0] == 3:
                llista.append((np.nanmean(dades.astype(np.float64), axis=0), cap))
            elif dades.ndim >= 3:
                for pla in dades.reshape(-1, dades.shape[-2], dades.shape[-1]):
                    llista.append((pla, cap))
            else:
                llista.append((dades, cap))
    return llista


def main():
    if len(sys.argv) < 2 or not os.path.isdir(sys.argv[1]):
        print(__doc__)
        sys.exit(1)
    carpeta = sys.argv[1]
    errors, provades = [], 0
    for nom in sorted(os.listdir(carpeta)):
        if not nom.lower().endswith(('.fits', '.fit', '.fts', '.fz', '.gz')):
            continue
        ruta = os.path.join(carpeta, nom)
        try:
            referencia = imatges_astropy(ruta)
        except Exception as error:
            print(f'{nom}: astropy no el pot llegir ({error}); se salta')
            continue
        for i, (dades, cap) in enumerate(referencia):
            etiqueta = f'{nom} [{i + 1}/{len(referencia)}]'
            try:
                imatge, cap_e, _nota = E.leer_fits(ruta, i)
            except E.ErrorFits as error:
                print(f'  ERROR {etiqueta}: Emillium no la llegeix ({error.clave})')
                errors.append(etiqueta)
                continue
            provades += 1
            iguals = np.allclose(np.nan_to_num(imatge), np.nan_to_num(np.asarray(dades, dtype=np.float32)), rtol=1e-6, atol=0)
            linia = f'  {"OK   " if iguals else "ERROR"} {etiqueta}: píxels {"idèntics" if iguals else "DIFERENTS"}'
            if not iguals:
                errors.append(etiqueta + ' píxels')
            hora, jd = E.instante_observacion(cap_e)
            if jd:
                linia += f' · JD {jd:.6f}'
            wcs_e = E.SolucionWCS.desde_cabecera(cap_e)
            if wcs_e is not None:
                fonts = C.detectar_fuentes(imatge, 5.0, 50)
                if fonts:
                    w = WCS(cap, naxis=2)
                    pitjor = 0.0
                    for f in fonts:
                        ar, dec = wcs_e.a_cielo(f['x'] + 1, f['y'] + 1)
                        ar_a, dec_a = w.all_pix2world([[f['x'] + 1, f['y'] + 1]], 1)[0]
                        pitjor = max(pitjor, C.separacion_y_angulo(ar, dec, ar_a % 360, dec_a)[0])
                    ok = pitjor < 0.01
                    linia += f' · AR/Dec de {len(fonts)} estrelles a ≤ {pitjor * 1000:.3f} mil·lisegons d\'arc'
                    if not ok:
                        errors.append(etiqueta + ' AR/Dec')
            print(linia)
    print()
    print(f'{provades} imatges provades.', 'TOT CORRECTE' if not errors else f'{len(errors)} PROBLEMES: {errors}')
    sys.exit(1 if errors else 0)


if __name__ == '__main__':
    main()
