"""
Proves de regressió d'Emillium Expanded amb el conjunt de fotos de referència.

Mesura les fotos de proves/referencia (creades per generar_referencia.py) i les
compara amb la veritat. Separa tres coses que no són el mateix:

  EXACTITUD           error respecte a la veritat (on és de debò l'estrella).
  REPETIBILITAT       dispersió quan es repeteix la mesura amb soroll diferent
                      (i si l'error ± que dona Emillium la prediu bé).
  PRECISIÓ NUMÈRICA   el càlcul en si: el mateix resultat cada vegada i el
                      mateix punt quan la foto es gira 90°.

Si alguna prova falla, el programa acaba amb codi 1 (útil per a proves
automàtiques). Només necessita numpy.

Ús:  python proves/prova_regressio.py
"""
import json
import math
import os
import sys

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(AQUI, '..', 'codigo'))
import Emillium as E  # noqa: E402
import ades  # noqa: E402
import ciencia as C  # noqa: E402

REF = os.path.join(AQUI, 'referencia')
errors = []


def comprova(condicio, categoria, text):
    print(f"  {'OK   ' if condicio else 'ERROR'} [{categoria}] {text}")
    if not condicio:
        errors.append(f'[{categoria}] {text}')


def llegeix(nom):
    imatge, cap, _nota = E.leer_fits(os.path.join(REF, nom + '.fits'))
    return imatge, cap


def mesura(imatge, cap, x, y, radi=5.0):
    m = E.medir_estrella(imatge, x, y, 6, radi, 3.0, E.ganancia_de(cap))
    if m is None:
        return None
    fot = C.fotometria_apertura(imatge, m['x'], m['y'], 5.0, 8.0, 14.0, E.ganancia_de(cap), C.nivel_saturacion(cap))
    m['flux_ap'], m['sigma_ap'] = fot['flujo'], fot['sigma_flujo']
    m['mag_ap'], m['sigma_mag_ap'] = C.magnitud(fot['flujo'], fot['sigma_flujo'], E.exposicion_de(cap))
    return m


def distancia_cel(a1, d1, a2, d2):
    return C.separacion_y_angulo(a1, d1, a2, d2)[0]


def main():
    if not os.path.exists(os.path.join(REF, 'referencia.json')):
        print('Falta el conjunt de referència: executa primer  python proves/generar_referencia.py')
        sys.exit(1)
    with open(os.path.join(REF, 'referencia.json'), encoding='utf-8') as fitxer:
        veritat = json.load(fitxer)['casos']

    print('== Estrelles brillants i febles (exactitud del centre, AR/Dec i magnituds)')
    for nom in ('brillants', 'febles'):
        imatge, cap = llegeix(nom)
        wcs = E.SolucionWCS(cap)
        errs, errs_cel, sigmes, dmags = [], [], [], []
        mesures = []
        for e in veritat[nom]['estrelles']:
            m = mesura(imatge, cap, e['x'] + 0.6, e['y'] - 0.4)
            mesures.append((m, e))
            errs.append(math.hypot(m['x'] - e['x'], m['y'] - e['y']))
            sigmes.append(math.hypot(m['sigma_x'], m['sigma_y']))
            ar, dec = wcs.a_cielo(m['x'] + 1, m['y'] + 1)
            errs_cel.append(distancia_cel(ar, dec, e['ar'], e['dec']))
        ref_m, ref_e = mesures[0]
        for m, e in mesures[1:]:
            dmags.append((m['mag_ap'] - ref_m['mag_ap']) - (e['mag'] - ref_e['mag']))
        # cada error ha de ser menys de 3 vegades el seu ± predit (i com a mínim 0,005 px)
        comprova(all(e < 3 * max(s, 0.005) for e, s in zip(errs, sigmes)), 'exactitud',
                 f'{nom}: centre a {max(errs):.4f} px de la veritat (màx.; ± predit {np.mean(sigmes):.4f} px)')
        comprova(all(e < 3 * max(s, 0.005) * 1.2 + 0.005 for e, s in zip(errs_cel, sigmes)), 'exactitud',
                 f'{nom}: AR/Dec a {max(errs_cel):.4f}″ de la veritat')
        # l'obertura de 5 px no recull tota la llum, però la fracció és la mateixa per a totes: les
        # diferències de magnitud han de ser correctes
        sig_m = [math.hypot(m['sigma_mag_ap'], ref_m['sigma_mag_ap']) for m, _e in mesures[1:]]
        comprova(all(abs(d) < 3 * max(s, 0.003) for d, s in zip(dmags, sig_m)), 'exactitud',
                 f'{nom}: diferències de magnitud a {max(abs(d) for d in dmags):.4f} mag (± predit {np.mean(sig_m):.4f})')

    print('== Repetibilitat (la mateixa estrella amb 60 sorolls diferents)')
    g = np.random.default_rng(7)
    xs, ys, sig, flux, sig_f = [], [], [], [], []
    yy, xx = np.mgrid[0:41, 0:41]
    s = 3.0 / 2.3548
    net = 800 + 4000 / (2 * math.pi * s * s) * np.exp(-((xx - 20.3) ** 2 + (yy - 19.8) ** 2) / (2 * s * s))
    for _ in range(60):
        im = g.poisson(net * 1.5) / 1.5
        m = E.medir_estrella(im, 20, 20, 6, 5.0, 3.0, 1.5)
        fot = C.fotometria_apertura(im, m['x'], m['y'], 5.0, 8.0, 14.0, 1.5)
        xs.append(m['x']); ys.append(m['y']); sig.append(m['sigma_x']); flux.append(fot['flujo']); sig_f.append(fot['sigma_flujo'])
    r_x = np.std(xs) / np.mean(sig)
    r_f = np.std(flux) / np.mean(sig_f)
    comprova(0.7 < r_x < 1.4, 'repetibilitat', f'dispersió del centre {np.std(xs):.4f} px vs ± predit {np.mean(sig):.4f} px (quocient {r_x:.2f})')
    comprova(0.7 < r_f < 1.4, 'repetibilitat', f'dispersió del flux {np.std(flux):.1f} vs ± predit {np.mean(sig_f):.1f} (quocient {r_f:.2f})')

    print('== Precisió numèrica')
    imatge, cap = llegeix('brillants')
    e = veritat['brillants']['estrelles'][0]
    m1 = E.medir_estrella(imatge, e['x'], e['y'], 6, 5.0, 3.0)
    m2 = E.medir_estrella(imatge, e['x'], e['y'], 6, 5.0, 3.0)
    comprova(m1['x'] == m2['x'] and m1['y'] == m2['y'], 'numèrica', 'el mateix resultat exacte si es repeteix')
    alt, ample = imatge.shape
    girada = E.girar_imagen(imatge, 90.0)
    xg, yg = E.transformar_punto(m1['x'], m1['y'], ample, alt, 90.0, False)
    mg = E.medir_estrella(girada, xg, yg, 6, 5.0, 3.0)
    xb, yb = E.destransformar_punto(mg['x'], mg['y'], ample, alt, 90.0, False)
    comprova(math.hypot(xb - m1['x'], yb - m1['y']) < 1e-6, 'numèrica',
             f'girant la foto 90° s\'obté el mateix punt (diferència {math.hypot(xb - m1["x"], yb - m1["y"]):.1e} px)')
    wcs = E.SolucionWCS(cap)
    xr, yr = wcs.a_pixel(*wcs.a_cielo(123.456, 78.9))
    comprova(abs(xr - 123.456) < 1e-7 and abs(yr - 78.9) < 1e-7, 'numèrica', 'píxel -> cel -> píxel torna al mateix punt')

    print('== Estrelles dobles (separació, angle i contaminació)')
    imatge, cap = llegeix('dobles')
    wcs = E.SolucionWCS(cap)
    for p in veritat['dobles']['parelles']:
        a, b = p['a'], p['b']
        sep_px = math.hypot(b['x'] - a['x'], b['y'] - a['y'])
        # radi de cerca 10 px: així es veu la veïna i se li treu la llum (com a l'aplicació)
        ma = E.medir_estrella(imatge, a['x'], a['y'], 10, 5.0, 3.0)
        mb = E.medir_estrella(imatge, b['x'], b['y'], 10, 5.0, 3.0)
        if sep_px < 4:
            comprova(True, 'límit', f'parella a {sep_px:.1f} px: massa junta per a la FWHM de 2,8 px (no es comprova)')
            continue
        sep_v, pa_v = C.separacion_y_angulo(a['ar'], a['dec'], b['ar'], b['dec'])
        sep, pa = C.separacion_y_angulo(*wcs.a_cielo(ma['x'] + 1, ma['y'] + 1), *wcs.a_cielo(mb['x'] + 1, mb['y'] + 1))
        comprova(abs(sep - sep_v) < 0.12 and abs((pa - pa_v + 180) % 360 - 180) < 1.5, 'exactitud',
                 f'parella a {sep_px:.0f} px: ρ {sep:.3f}″ (veritat {sep_v:.3f}), θ {pa:.2f}° (veritat {pa_v:.2f})')
        fot = C.fotometria_apertura(imatge, ma['x'], ma['y'], 5.0, 8.0, 14.0)
        fraccio, _v = C.contaminacion(imatge, ma['x'], ma['y'], 5.0, 2.8, fot['flujo'], fondo=fot['fondo'], sigma=fot['sigma_fondo'])
        if sep_px <= 6:
            comprova(fraccio > 0.01, 'qualitat', f'parella a {sep_px:.0f} px marcada com a contaminada ({fraccio * 100:.1f} %)')

    print('== Estrella variable (fotometria diferencial amb 3 comparacions)')
    mags, errs_m, errs = [], [], []
    for k in range(1, 7):
        nom = f'variable_{k}'
        imatge, cap = llegeix(nom)
        est = veritat[nom]['estrelles']
        mesures = {n: mesura(imatge, cap, e['x'], e['y']) for n, e in est.items()}
        comps = [(mesures[n]['mag_ap'], mesures[n]['sigma_mag_ap'], est[n]['mag'], None) for n in ('C1', 'C2', 'C3')]
        mag, err = C.magnitud_diferencial(mesures['V']['mag_ap'], mesures['V']['sigma_mag_ap'], C.punto_cero(comps))
        mags.append(mag); errs_m.append(err); errs.append(mag - est['V']['mag'])
    comprova(max(abs(e) for e in errs) < 0.02, 'exactitud', f'magnitud de la variable a {max(abs(e) for e in errs):.4f} mag de la veritat')
    comprova(C.analisis_variabilidad(mags, errs_m)['variable'], 'exactitud', 'la variable es detecta com a variable')

    print('== Asteroide (posició i moviment)')
    jds, ars, decs = [], [], []
    for k in range(1, 5):
        nom = f'asteroide_{k}'
        imatge, cap = llegeix(nom)
        v = veritat[nom]
        hora, jd = E.instante_observacion(cap)
        comprova(abs(jd - v['jd_mig']) * 86400 < 0.01, 'exactitud', f'{nom}: instant del mig de l\'exposició ({hora})')
        m = E.medir_estrella(imatge, v['asteroide']['x'], v['asteroide']['y'], 3, 5.0, 3.0)
        ar, dec = E.SolucionWCS(cap).a_cielo(m['x'] + 1, m['y'] + 1)
        jds.append(jd); ars.append(ar); decs.append(dec)
    vel, ang, rms, _r = C.movimiento_en_cielo(jds, ars, decs)
    v_ver, ang_ver, _rms, _r = C.movimiento_en_cielo([veritat[f'asteroide_{k}']['jd_mig'] for k in range(1, 5)],
                                                     [veritat[f'asteroide_{k}']['asteroide']['ar'] for k in range(1, 5)],
                                                     [veritat[f'asteroide_{k}']['asteroide']['dec'] for k in range(1, 5)])
    comprova(abs(vel - v_ver) / v_ver < 0.02 and abs(ang - ang_ver) < 1.0, 'exactitud',
             f'velocitat {vel:.2f}″/h (veritat {v_ver:.2f}) cap a {ang:.1f}° (veritat {ang_ver:.1f}°), rms {rms:.3f}″')
    context = {'designacion': 'EMLREF1', 'codigo': 'XXX', 'modo': 'CCD', 'remitente': 'Proves', 'observadores': 'Proves',
               'medidores': 'Proves', 'apertura': '0.2', 'diseno': 'Reflector', 'catalogo_astrometrico': 'Gaia3', 'banda': 'V'}
    obs = [{'jd': j, 'ar': a, 'dec': d, 'rms_ar': 0.1, 'rms_dec': 0.1} for j, a, d in zip(jds, ars, decs)]
    problemes = [p for p in ades.validar(context, obs) if p[0] == 'error']
    comprova(not problemes, 'format', f'informe ADES de l\'asteroide vàlid {problemes or ""}')

    print('== Foto desenfocada, camp dens i escales')
    stats = [C.estadisticas_imagen(*[llegeix(n)[0], 65000.0]) for n in ('brillants', 'variable_1', 'asteroide_1', 'desenfocada')]
    motius = C.fotos_defectuosas(stats)
    comprova('desenfocada' in motius[3] and not motius[0], 'qualitat', f'foto desenfocada detectada (FWHM {stats[3]["fwhm"]:.1f} px vs {stats[0]["fwhm"]:.1f})')
    imatge, cap = llegeix('camp_dens')
    fonts = C.detectar_fuentes(imatge)
    est = veritat['camp_dens']['estrelles']
    v = np.array([(e['x'], e['y']) for e in est])
    brillants = [e for e in est if e['mag'] < 15.0]
    trobada = lambda e: min(math.hypot(f['x'] - e['x'], f['y'] - e['y']) for f in fonts) < 1.0   # noqa: E731
    # Estrelles "separables": sense cap altra estrella a menys d'1,5 vegades la FWHM. Les que en
    # tenen una de més a prop es fonen en una sola taca (límit de la detecció automàtica).
    limit = 1.5 * veritat['camp_dens']['fwhm']
    separables = [e for e in brillants if np.sort(np.hypot(v[:, 0] - e['x'], v[:, 1] - e['y']))[1] > limit]
    trobades = sum(1 for e in separables if trobada(e))
    totals = sum(1 for e in brillants if trobada(e))
    falses = sum(1 for f in fonts if np.min(np.hypot(v[:, 0] - f['x'], v[:, 1] - f['y'])) > 2.0)
    comprova(trobades / len(separables) > 0.90 and falses < 0.05 * len(fonts), 'exactitud',
             f'camp dens: {trobades}/{len(separables)} estrelles separables trobades ({totals}/{len(brillants)} de '
             f'totes les de mag < 15), {falses} deteccions falses de {len(fonts)}')
    for nom, tol in (('escala_fina', 0.05), ('escala_gruixuda', 0.25)):
        imatge, cap = llegeix(nom)
        wcs = E.SolucionWCS(cap)
        err = [distancia_cel(*wcs.a_cielo(m['x'] + 1, m['y'] + 1), e['ar'], e['dec'])
               for e in veritat[nom]['estrelles'] for m in [E.medir_estrella(imatge, e['x'], e['y'], 4, 5.0, 3.0)]]
        comprova(max(err) < tol, 'exactitud', f"{nom} ({veritat[nom]['escala']}″/px, FWHM {veritat[nom]['fwhm']} px): AR/Dec a {max(err):.3f}″ (màx.)")

    print('== Resoldre una foto sense WCS, cercar mòbils i calibrar amb el camp')
    imatge, cap = llegeix('camp_dens')
    wcs = E.SolucionWCS(cap)
    est = sorted(veritat['camp_dens']['estrelles'], key=lambda e: e['mag'])
    cel = [wcs.a_cielo(e['x'] + 1, e['y'] + 1) for e in est]
    alt, ample = imatge.shape
    fonts = C.detectar_fuentes(imatge)
    # apuntat desviat 30″ i sense escala
    r = C.resolver_campo(fonts, ample, alt, [c[0] for c in cel], [c[1] for c in cel], 150.12 + 0.008, 22.45 - 0.004)
    if r is None:
        comprova(False, 'astrometria', 'resoldre el camp dens sense WCS')
    else:
        nova = E.SolucionWCS(r['cabecera'])
        err = max(C.separacion_y_angulo(*nova.a_cielo(x, y), *wcs.a_cielo(x, y))[0]
                  for x, y in ((1, 1), (ample, 1), (1, alt), (ample, alt)))
        comprova(err < 0.15, 'astrometria', f"camp resolt sense escala: {r['n_usadas']} estrelles, rms {r['rms']:.3f}″, "
                                            f"escala {r['escala']:.4f}″/px, error a les cantonades {err:.3f}″")
    llistes, desplaçaments, temps = [], [], []
    for k in range(1, 5):
        im, c = llegeix(f'asteroide_{k}')
        llistes.append(C.detectar_fuentes(im))
        temps.append(veritat[f'asteroide_{k}']['jd_mig'])
    for k in range(4):
        d = C.desplazamiento_entre(llistes[0], llistes[k]) if k else (0.0, 0.0, 0)
        desplaçaments.append(d[:2] if d else (0.0, 0.0))
    pistes = C.buscar_moviles(llistes, desplaçaments, temps)
    if pistes:
        p = max(pistes, key=lambda p: p['n'])
        err = max(math.hypot(x - veritat[f'asteroide_{k + 1}']['asteroide']['x'], y - veritat[f'asteroide_{k + 1}']['asteroide']['y'])
                  for k, (x, y, _f) in p['posiciones'].items())
        comprova(len(pistes) == 1 and p['n'] == 4 and err < 0.5, 'exactitud',
                 f"asteroide trobat sol: {len(pistes)} pista, {p['n']} fotos, a {err:.2f} px de la veritat")
    else:
        comprova(False, 'exactitud', 'asteroide trobat sol')
    # Calibració amb estrelles aïllades de diverses fotos de referència (mateixa càmera i exposició)
    m_inst, m_cat, errs = [], [], []
    for nom, estrelles in (('brillants', veritat['brillants']['estrelles']), ('febles', veritat['febles']['estrelles']),
                           ('variable_1', [veritat['variable_1']['estrelles'][n] for n in ('V', 'C1', 'C2', 'C3')]),
                           ('escala_fina', veritat['escala_fina']['estrelles'])):
        im, c = llegeix(nom)
        for e in estrelles:
            fot = C.fotometria_apertura(im, e['x'], e['y'], 8.0, 11.0, 17.0, E.ganancia_de(c), C.nivel_saturacion(c))
            m, s = C.magnitud(fot['flujo'], fot['sigma_flujo'], E.exposicion_de(c))
            m_inst.append(m); m_cat.append(e['mag']); errs.append(s)
    cal = C.calibrar_campo(m_inst, m_cat, None, errs)
    # punt zero esperat: 20,5 − 25 = −4,5 (l'obertura de 8 px recull tota la llum)
    comprova(abs(cal['zp'] + 4.5) < 0.01 and cal['rms'] < 0.05, 'exactitud',
             f"calibració: ZP {cal['zp']:.4f} ± {cal['sigma_zp']:.4f} (esperat −4,50), rms {cal['rms']:.3f} mag, {cal['n']} estrelles")

    print('== Temps heliocèntric, període, mínim d\'eclipsi i distorsió SIP')
    # valors de referència calculats amb astropy (correcció heliocèntrica, en segons)
    for jd, ar, dec, segons in ((2461110.4, 150.12, 22.45, 440.738), (2451545.0, 10.0, -60.0, -171.346),
                                (2455000.25, 279.23, 38.78, 227.366)):
        diferencia = (C.jd_heliocentrico(jd, ar, dec) - jd) * 86400 - segons
        comprova(abs(diferencia) < 0.5, 'exactitud', f'HJD per a JD {jd}: {diferencia:+.2f} s respecte a astropy')
    g = np.random.default_rng(11)
    t = np.concatenate([n + np.sort(g.uniform(0, 0.25, 35)) for n in range(0, 10, 2)])
    m = 12 + 0.2 * np.sin(2 * np.pi * t / 0.3712) + g.normal(0, 0.01, len(t))
    _p, _pot, periode = C.periodograma(t, m, np.full(len(t), 0.01), 0.05, 5)
    comprova(abs(periode - 0.3712) / 0.3712 < 0.003, 'exactitud', f'període trobat {periode:.5f} d (veritat 0,3712 d)')
    t0 = 2461110.4567
    t = np.sort(2461110.40 + g.uniform(0, 0.12, 120))
    m = 12 + 0.6 * np.exp(-((t - t0) / 0.012) ** 2) + g.normal(0, 0.006, len(t))
    minim = C.minimo_kvw(t, m)
    comprova(minim is not None and abs(minim[0] - t0) * 86400 < 3 * max(minim[1] * 86400, 5),
             'exactitud', f"mínim d'eclipsi a {(minim[0] - t0) * 86400:+.1f} s de la veritat (± predit {minim[1] * 86400:.1f} s)")
    # distorsió de barril coneguda: la solució lineal no basta i la SIP sí
    u = g.uniform(-1000, 1000, 400)
    v = g.uniform(-700, 700, 400)
    du, dv = 1e-6 * u * u + 4e-9 * (u ** 3 + u * v * v), -1e-6 * v * v + 4e-9 * (v ** 3 + u * u * v)
    escala = 3.5 / 3600
    ar, dec = C.desproyectar(-(u + du) * escala, (v + dv) * escala, 83.8, -5.4)
    lineal = C.ajustar_wcs(u + 1024.5, v + 768.5, ar, dec, crpix=(1024.5, 768.5))
    sip = C.ajustar_sip(u + 1024.5, v + 768.5, ar, dec, orden=3, crpix=(1024.5, 768.5))
    comprova(sip is not None and sip['rms'] < 0.02 and lineal['rms'] > 1, 'astrometria',
             f"distorsió: lineal rms {lineal['rms']:.2f}″ → SIP rms {sip['rms'] if sip else float('nan'):.4f}″")

    print('== Detecció automàtica en fotos difícils')
    g = np.random.default_rng(11)
    yy, xx = np.mgrid[0:600, 0:700]
    def estrella(cx, cy, a, s, sy=None):
        return a * np.exp(-((xx - cx) ** 2 / (2 * s * s) + (yy - cy) ** 2 / (2 * (sy or s) ** 2)))
    foto = g.normal(1000, 10, yy.shape) + 0.6 * xx                          # degradat fort
    veritat_det = [(300, 300), (450, 450), (200, 450), (560, 180)]
    foto += estrella(300, 300, 3000, 1.5) + estrella(450, 450, 80, 1.3)     # normal i feble (~8σ)
    foto += estrella(200, 450, 3000, 6) + estrella(560, 180, 1500, 4, 1.3)  # desenfocada i moguda
    foto += estrella(150, 150, 2e5, 4)                                      # saturada (centre pla)
    foto = np.minimum(foto, 30000)
    veritat_det.append((150, 150))
    foto[440:480, 20:80] = 9000                                             # rètol / bloc brillant
    foto[:, 640:643] += 400                                                 # línia (vora del sensor)
    foto[:60, :] = 0                                                        # zona sense dades
    foto[60:, :40] = 0
    fonts = C.detectar_fuentes(foto.astype(np.float32), saturacion=30000)
    trobades = [min((math.hypot(f['x'] - x, f['y'] - y) for f in fonts), default=99) for x, y in veritat_det]
    falses = [f for f in fonts if min(math.hypot(f['x'] - x, f['y'] - y) for x, y in veritat_det) > 3]
    comprova(max(trobades) < 0.6 and not falses, 'exactitud',
             f'{len(veritat_det)} estrelles trobades (la saturada a {trobades[-1]:.2f} px del centre), '
             f'{len(falses)} falses (degradat, bloc brillant, línia i zona sense dades)')
    # Estrella molt brillant amb anell i aspes (com les del Hubble) i una companya de veritat
    foto = g.normal(1000, 10, yy.shape) + estrella(350, 300, 40000, 1.5) + estrella(377, 291, 600, 1.5)
    rr, aa = np.hypot(xx - 350, yy - 300), np.arctan2(yy - 300, xx - 350)
    foto += 400 * np.exp(-(rr - 12) ** 2 / 2) * (1 + np.cos(8 * aa)) / 2          # anell amb 8 grumolls
    for angle in (45, 135, 225, 315):                                           # 4 aspes amb grumolls
        dist = (xx - 350) * math.cos(math.radians(angle)) + (yy - 300) * math.sin(math.radians(angle))
        perp = -(xx - 350) * math.sin(math.radians(angle)) + (yy - 300) * math.cos(math.radians(angle))
        foto += np.where(dist > 8, 3000 / np.maximum(dist, 1) * (1 + np.cos(dist)) / 2, 0) * np.exp(-perp ** 2 / 2)
    fonts = C.detectar_fuentes(foto.astype(np.float32))
    reals = [(350, 300), (377, 291)]
    trobades = [min((math.hypot(f['x'] - x, f['y'] - y) for f in fonts), default=99) for x, y in reals]
    falses = [f for f in fonts if min(math.hypot(f['x'] - x, f['y'] - y) for x, y in reals) > 3]
    comprova(max(trobades) < 1 and not falses, 'exactitud',
             f'estrella brillant amb anell i aspes: {len(fonts)} deteccions ({len(falses)} falses en l\'anell o les aspes)')
    print('== Temps baricèntric (BJD_TDB)')
    # Valors d'astropy (Time.tdb + light_travel_time 'barycentric', efemèride 'builtin')
    for (jd, ar, dec), bjd_astropy in (((2461110.4, 83.8, -5.4), 2461110.40103487),
                                       ((2451545.0, 10.0, 41.0), 2451545.00207582),
                                       ((2455000.25, 280.0, -20.0), 2455000.25645312)):
        diferencia = (C.bjd_tdb(jd, ar, dec) - bjd_astropy) * 86400
        comprova(abs(diferencia) < 0.5, 'exactitud', f'BJD_TDB per a JD {jd}: {diferencia:+.2f} s respecte a astropy')

    print('== Alinear i apilar fotos')
    yy, xx = np.mgrid[0:300, 0:400]
    cel = list(zip(g.uniform(-100, 500, 70), g.uniform(-100, 400, 70), 10 ** g.uniform(2.8, 4.6, 70)))
    def foto_vista(A, t, rai=False):
        inv = np.linalg.inv(A)
        im = g.normal(500, 10, yy.shape)
        for x, y, f in cel:
            px, py = inv @ (np.array([x, y]) - t)
            im += f / (2 * math.pi * 1.7 ** 2) * np.exp(-((xx - px) ** 2 + (yy - py) ** 2) / (2 * 1.7 ** 2))
        if rai:
            im[150, 60:64] += 4000
        return im.astype(np.float32)
    gir = lambda a: np.array([[math.cos(a), -math.sin(a)], [math.sin(a), math.cos(a)]])   # noqa: E731
    casos = [(np.eye(2), np.zeros(2)), (np.eye(2), np.array([6.4, -3.7])), (gir(0.06), np.array([-9.0, 8.0])),
             (-np.eye(2), np.array([401.3, 297.6])), (np.array([[-1.0, 0], [0, 1.0]]), np.array([396.2, 1.9]))]
    fotos = [foto_vista(A, t, rai=(k == 2)) for k, (A, t) in enumerate(casos)]
    referencia = C.fuentes_para_alinear(fotos[0])
    transformacions, pitjor = [], 0.0
    punts = np.array([[30, 30], [370, 270], [200, 150]], dtype=np.float64)
    for (A, t), im in zip(casos, fotos):
        r = C.alinear_fuentes(referencia, C.fuentes_para_alinear(im))
        if r is None:
            pitjor = 99.0
            break
        pitjor = max(pitjor, float(np.abs(punts @ r['A'].T + r['t'] - (punts @ A.T + t)).max()))
        transformacions.append((r['A'], r['t']))
    comprova(pitjor < 0.1, 'exactitud', f'5 fotos (desplaçades, girades, a 180° i en mirall) alineades a ≤ {pitjor:.3f} px')
    if len(transformacions) == len(fotos):
        apilada, _n = C.apilar(fotos, transformacions, fotos[0].shape, 'media')
        s1, s5 = C.fondo_y_ruido(fotos[0])[1], C.fondo_y_ruido(apilada)[1]
        comprova(s1 / s5 > 2.0, 'exactitud', f'apilant 5 fotos el soroll baixa {s1 / s5:.2f} vegades (√5 = 2,24)')
        comprova(apilada[150, 61] - C.fondo_y_ruido(apilada)[0] < 6 * s5, 'qualitat', 'el raig còsmic desapareix de la foto apilada')
        fluxos = [(E.medir_estrella(fotos[0], x, y), E.medir_estrella(apilada, x, y)) for x, y, f in cel
                  if 30 < x < 370 and 30 < y < 270 and f > 8000]
        dif = [2.5 * math.log10(a['flujo'] / b['flujo']) for a, b in fluxos if a and b]
        comprova(dif and abs(float(np.median(dif))) < 0.01, 'exactitud',
                 f'la foto apilada conserva el flux de les estrelles ({float(np.median(dif)):+.4f} mag)')
    soroll = C.detectar_fuentes(g.normal(100, 10, (1000, 1000)).astype(np.float32))
    comprova(len(soroll) <= 5, 'exactitud', f'foto de només soroll (1 megapíxel): {len(soroll)} deteccions falses')

    print()
    if errors:
        print(f'{len(errors)} PROVES FALLADES:')
        for e in errors:
            print('  -', e)
        sys.exit(1)
    print('TOT CORRECTE')


if __name__ == '__main__':
    main()
