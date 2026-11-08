"""
Emillium Expanded - informes ADES para el Minor Planet Center
=============================================================

ADES (Astrometry Data Exchange Standard) es el formato que el MPC recomienda desde
2018 para enviar observaciones de asteroides y cometas. Tiene dos formas con la
misma información:
  - PSV ("pipe-separated values"): texto con columnas separadas por "|".
  - XML.

Un informe tiene una cabecera (quién observa, con qué telescopio, desde qué
observatorio) y una fila por observación (objeto, instante, AR, Dec, errores,
catálogo de referencia, magnitud...).

Este módulo solo escribe y valida el texto; no envía nada a ningún sitio.
Validar un informe NO garantiza que el MPC acepte las observaciones: también cuentan
la identificación del objeto, la calidad de las medidas y sus normas de envío.

Las comprobaciones devuelven una lista de (nivel, código, datos), donde nivel es
'error' (el informe no es válido) o 'aviso' (es válido pero conviene revisarlo).
Los textos de cada código, traducidos, están en Emillium.py (claves 'ades_<código>').
"""

import math
import re
from datetime import datetime, timezone
from xml.sax.saxutils import escape

VERSION_ADES = '2017'

# Catálogos de referencia astrométrica y fotométrica más habituales (códigos ADES)
CATALOGOS_ASTROMETRICOS = ['Gaia3', 'Gaia3E', 'Gaia2', 'Gaia1', 'Gaia_2016', 'UCAC4', 'URAT1', 'PPMXL',
                           'ATLAS2', 'PS1_DR1', 'PS1_DR2', 'Tycho2', 'Hipparcos2']
CATALOGOS_FOTOMETRICOS = ['Gaia3', 'Gaia2', 'Gaia1', 'ATLAS2', 'PS1_DR2', 'PS1_DR1', 'UCAC4', 'URAT1', 'APASS']
BANDAS = ['G', 'Gb', 'Gr', 'V', 'B', 'R', 'I', 'U', 'C', 'o', 'c', 'w', 'g', 'r', 'i', 'z', 'y', 'J', 'H', 'K']
MODOS = ['CCD', 'CMO']
DISENOS = ['Reflector', 'Refractor', 'Schmidt', 'Astrograph', 'Catadioptric']

_RE_ESTACION = re.compile(r'^[0-9A-Z]{3}$')
_RE_PROVISIONAL = re.compile(r'^\d{4} [A-Z]{2}\d{0,4}$')       # 2024 AB12
_RE_COMETA = re.compile(r'^[PCDXA]/\d{4} [A-Z]{1,2}\d{0,3}$')   # C/2023 A3
_RE_PERMANENTE = re.compile(r'^\d{1,7}[PCDXI]?$')              # 433, 1P
_RE_TRKSUB = re.compile(r'^[A-Za-z0-9]{1,8}$')


def tipo_designacion(designacion):
    """Clasifica la designación: 'permID' (número), 'provID' (provisional del MPC) o
    'trkSub' (nombre propio del observador para un objeto aún sin designación)."""
    d = designacion.strip()
    if _RE_PERMANENTE.match(d):
        return 'permID'
    if _RE_PROVISIONAL.match(d) or _RE_COMETA.match(d):
        return 'provID'
    return 'trkSub'


def instante_iso(jd, decimales=2):
    """Día juliano -> 'AAAA-MM-DDThh:mm:ss.ssZ' (UTC)."""
    segundos_totales = round((jd - 2440587.5) * 86400, decimales)        # desde 1970-01-01
    entero = math.floor(segundos_totales)
    fecha = datetime.fromtimestamp(entero, tz=timezone.utc)
    fraccion = segundos_totales - entero
    texto_seg = f'{fecha.second + fraccion:0{3 + decimales}.{decimales}f}'
    return fecha.strftime('%Y-%m-%dT%H:%M:') + texto_seg + 'Z'


def _numero(valor, decimales, signo=False):
    if valor is None:
        return ''
    return f'{valor:+.{decimales}f}' if signo else f'{valor:.{decimales}f}'


# Campos de una observación en el orden que pide ADES
CAMPOS = ['permID', 'provID', 'trkSub', 'mode', 'stn', 'obsTime', 'ra', 'dec', 'rmsRA', 'rmsDec', 'rmsCorr',
          'astCat', 'mag', 'rmsMag', 'band', 'photCat', 'disc', 'notes', 'remarks']


def filas_observaciones(contexto, observaciones):
    """Convierte las observaciones en filas de texto (diccionario campo -> texto).

    observaciones: lista de diccionarios con jd, ar, dec (grados), rms_ar y rms_dec
    (segundos de arco; rms_ar ya multiplicado por cos(Dec)), y opcionalmente mag y rms_mag."""
    tipo = tipo_designacion(contexto['designacion'])
    filas = []
    for i, o in enumerate(observaciones):
        f = dict.fromkeys(CAMPOS, '')
        f[tipo] = contexto['designacion'].strip()
        f['mode'] = contexto['modo']
        f['stn'] = contexto['codigo'].strip().upper()
        f['obsTime'] = instante_iso(o['jd'])
        f['ra'] = _numero(o['ar'] % 360.0, 6)
        f['dec'] = _numero(o['dec'], 6, signo=True)
        f['rmsRA'] = _numero(o.get('rms_ar'), 3)
        f['rmsDec'] = _numero(o.get('rms_dec'), 3)
        f['astCat'] = contexto['catalogo_astrometrico']
        if o.get('mag') is not None:
            f['mag'] = _numero(o['mag'], 2)
            f['rmsMag'] = _numero(o.get('rms_mag'), 2)
            f['band'] = contexto['banda']
            f['photCat'] = contexto.get('catalogo_fotometrico', '')
        if contexto.get('descubrimiento') and i == 0:
            f['disc'] = '*'     # marca de descubrimiento: solo en una observación (la primera)
        filas.append(f)
    return filas


def _cabecera_psv(contexto):
    lineas = [f'# version={VERSION_ADES}', '# observatory', f"! mpcCode {contexto['codigo'].strip().upper()}"]
    if contexto.get('nombre_observatorio', '').strip():
        lineas.append(f"! name {contexto['nombre_observatorio'].strip()}")
    lineas += ['# submitter', f"! name {contexto['remitente'].strip()}"]
    if contexto.get('institucion', '').strip():
        lineas.append(f"! institution {contexto['institucion'].strip()}")
    lineas.append('# observers')
    lineas += [f'! name {n}' for n in _nombres(contexto['observadores'])]
    lineas.append('# measurers')
    lineas += [f'! name {n}' for n in _nombres(contexto['medidores'])]
    lineas += ['# telescope', f"! aperture {contexto['apertura']}", f"! design {contexto['diseno']}",
               f"! detector {contexto['modo']}"]
    if contexto.get('comentario', '').strip():
        lineas += ['# comment', f"! line {contexto['comentario'].strip()}"]
    return lineas


def _nombres(texto):
    return [n.strip() for n in re.split(r'[;,]', texto or '') if n.strip()]


def crear_psv(contexto, observaciones):
    """Informe ADES en formato PSV (columnas alineadas, separadas por "|")."""
    filas = filas_observaciones(contexto, observaciones)
    usados = [c for c in CAMPOS if c in ('mode', 'stn', 'obsTime', 'ra', 'dec', 'astCat')
              or any(f[c] for f in filas)]
    anchos = {c: max([len(c)] + [len(f[c]) for f in filas]) for c in usados}
    lineas = _cabecera_psv(contexto)
    lineas.append('|'.join(c.ljust(anchos[c]) for c in usados).rstrip())
    for f in filas:
        lineas.append('|'.join(f[c].ljust(anchos[c]) for c in usados).rstrip())
    return '\n'.join(lineas) + '\n'


def crear_xml(contexto, observaciones):
    """El mismo informe en formato XML."""
    filas = filas_observaciones(contexto, observaciones)
    e = lambda t: escape(str(t).strip())                                        # noqa: E731
    x = ['<?xml version="1.0" encoding="UTF-8"?>', f'<ades version="{VERSION_ADES}">', '  <obsBlock>',
         '    <obsContext>', '      <observatory>', f"        <mpcCode>{e(contexto['codigo'].upper())}</mpcCode>"]
    if contexto.get('nombre_observatorio', '').strip():
        x.append(f"        <name>{e(contexto['nombre_observatorio'])}</name>")
    x += ['      </observatory>', '      <submitter>', f"        <name>{e(contexto['remitente'])}</name>"]
    if contexto.get('institucion', '').strip():
        x.append(f"        <institution>{e(contexto['institucion'])}</institution>")
    x.append('      </submitter>')
    for bloque, clave in (('observers', 'observadores'), ('measurers', 'medidores')):
        x.append(f'      <{bloque}>')
        x += [f'        <name>{e(n)}</name>' for n in _nombres(contexto[clave])]
        x.append(f'      </{bloque}>')
    x += ['      <telescope>', f"        <aperture>{e(contexto['apertura'])}</aperture>",
          f"        <design>{e(contexto['diseno'])}</design>", f"        <detector>{e(contexto['modo'])}</detector>",
          '      </telescope>']
    if contexto.get('comentario', '').strip():
        x += ['      <comment>', f"        <line>{e(contexto['comentario'])}</line>", '      </comment>']
    x += ['    </obsContext>', '    <obsData>']
    for f in filas:
        x.append('      <optical>')
        x += [f'        <{c}>{e(f[c])}</{c}>' for c in CAMPOS if f[c]]
        x.append('      </optical>')
    x += ['    </obsData>', '  </obsBlock>', '</ades>']
    return '\n'.join(x) + '\n'


# =============================================================================
# Validación
# =============================================================================


def validar(contexto, observaciones, ahora_jd=None):
    """Comprueba el informe antes de exportarlo. Devuelve [(nivel, código, datos)]."""
    p = []
    codigo = contexto.get('codigo', '').strip().upper()
    if not _RE_ESTACION.match(codigo):
        p.append(('error', 'codigo_formato', {'c': codigo}))
    elif codigo == 'XXX':
        p.append(('aviso', 'codigo_xxx', {}))
    if not contexto.get('remitente', '').strip():
        p.append(('error', 'falta_remitente', {}))
    if not _nombres(contexto.get('observadores')):
        p.append(('error', 'falta_observadores', {}))
    if not _nombres(contexto.get('medidores')):
        p.append(('error', 'falta_medidores', {}))
    try:
        apertura = float(str(contexto.get('apertura', '')).replace(',', '.'))
        if not 0.01 <= apertura <= 20:
            p.append(('error', 'apertura_rango', {'a': apertura}))
        elif apertura > 3:
            p.append(('aviso', 'apertura_metros', {'a': apertura}))
    except ValueError:
        p.append(('error', 'apertura_numero', {}))
    if contexto.get('diseno') not in DISENOS:
        p.append(('aviso', 'diseno_raro', {'d': contexto.get('diseno', '')}))
    if contexto.get('modo') not in MODOS:
        p.append(('error', 'modo', {'m': contexto.get('modo', '')}))
    if contexto.get('catalogo_astrometrico') not in CATALOGOS_ASTROMETRICOS:
        p.append(('aviso', 'catalogo_raro', {'c': contexto.get('catalogo_astrometrico', '')}))
    designacion = contexto.get('designacion', '').strip()
    tipo = tipo_designacion(designacion) if designacion else None
    if not designacion:
        p.append(('error', 'falta_designacion', {}))
    elif tipo == 'trkSub' and not _RE_TRKSUB.match(designacion):
        p.append(('error', 'trksub_formato', {'d': designacion}))

    if not observaciones:
        p.append(('error', 'sin_observaciones', {}))
        return p
    if len(observaciones) < 2:
        p.append(('error', 'pocas_observaciones', {'n': len(observaciones)}))
    elif len(observaciones) < 3:
        p.append(('aviso', 'recomendable_3', {'n': len(observaciones)}))
    vistos = set()
    if ahora_jd is None:
        ahora_jd = datetime.now(timezone.utc).timestamp() / 86400 + 2440587.5
    for i, o in enumerate(observaciones, 1):
        if not (0 <= o['ar'] < 360 or math.isclose(o['ar'], 360)):
            p.append(('error', 'ar_rango', {'i': i}))
        if not -90 <= o['dec'] <= 90:
            p.append(('error', 'dec_rango', {'i': i}))
        if o['jd'] > ahora_jd + 0.01:
            p.append(('error', 'fecha_futura', {'i': i}))
        if o['jd'] < 2415020.5:        # antes de 1900
            p.append(('error', 'fecha_antigua', {'i': i}))
        instante = instante_iso(o['jd'])
        if instante in vistos:
            p.append(('error', 'fecha_repetida', {'i': i, 't': instante}))
        vistos.add(instante)
        for clave in ('rms_ar', 'rms_dec'):
            valor = o.get(clave)
            if valor is None:
                p.append(('aviso', 'sin_rms', {'i': i}))
                break
            if valor <= 0:
                p.append(('error', 'rms_cero', {'i': i}))
                break
            if valor > 5:
                p.append(('aviso', 'rms_grande', {'i': i, 'r': valor}))
                break
        if o.get('mag') is not None:
            if not -2 < o['mag'] < 30:
                p.append(('error', 'mag_rango', {'i': i}))
            if contexto.get('banda') not in BANDAS:
                p.append(('error', 'banda', {'b': contexto.get('banda', '')}))
            if not contexto.get('catalogo_fotometrico'):
                p.append(('aviso', 'sin_photcat', {}))
    # ¿Las posiciones siguen un movimiento razonable (casi en línea recta y a velocidad
    # constante durante una noche)? Si no, alguna medida puede ser de otro objeto.
    if len(observaciones) >= 3:
        orden = sorted(observaciones, key=lambda o: o['jd'])
        if orden[-1]['jd'] - orden[0]['jd'] < 1.0:
            residuo = _residuo_lineal(orden)
            if residuo is not None and residuo > max(1.5, 4 * max(o.get('rms_ar') or 0.3 for o in orden)):
                p.append(('aviso', 'movimiento_irregular', {'r': residuo}))
    return p


def _residuo_lineal(observaciones):
    """Máximo residuo (″) respecto a un movimiento lineal en el cielo."""
    jd0 = observaciones[0]['jd']
    ar0 = observaciones[0]['ar']
    dec0 = observaciones[0]['dec']
    ts = [o['jd'] - jd0 for o in observaciones]
    xs = [((o['ar'] - ar0 + 180) % 360 - 180) * math.cos(math.radians(dec0)) * 3600 for o in observaciones]
    ys = [(o['dec'] - dec0) * 3600 for o in observaciones]
    n = len(ts)
    tm = sum(ts) / n
    var_t = sum((t - tm) ** 2 for t in ts)
    if var_t == 0:
        return None
    peor = 0.0
    for valores in (xs, ys):
        vm = sum(valores) / n
        pendiente = sum((t - tm) * (v - vm) for t, v in zip(ts, valores)) / var_t
        peor = max(peor, max(abs(v - (vm + pendiente * (t - tm))) for t, v in zip(ts, valores)))
    return peor


def hay_errores(problemas):
    return any(nivel == 'error' for nivel, _c, _d in problemas)
