"""
Emillium Expanded - catálogos astronómicos con caché local
==========================================================

REGLAS (requisito del proyecto):
  * Emillium funciona sin internet. Nada de este módulo se ejecuta solo: las
    descargas empiezan únicamente cuando el usuario pulsa un botón de descarga.
  * Lo descargado se guarda en una carpeta local y se puede usar para siempre sin
    volver a conectarse (%APPDATA%\\Emillium\\catalogos en Windows).
  * Sin telemetría: solo se piden las estrellas de la zona del cielo de la foto.
    No se envía nada más.

Catálogos
---------
  gaia    Gaia DR3 (ESA): posiciones de altísima precisión (época 2016,0),
          movimientos propios, paralajes y magnitud G y color BP−RP.
  apass   APASS DR9 (AAVSO): magnitudes B, V, g', r', i' para fotometría.
  vsx     VSX (AAVSO): catálogo de estrellas variables conocidas.
  skybot  SkyBoT (IMCCE): asteroides y cometas conocidos en el campo en el instante
          de la foto.
  local   Catálogos propios importados de un archivo CSV (sin conexión).
"""

import csv
import io
import json
import math
import os
import re
import time
import urllib.parse
import urllib.request

AGENTE = 'Emillium/1.5 (desktop astronomy software)'
ESPEJOS_VIZIER = ['https://vizier.cds.unistra.fr', 'https://vizier.cfa.harvard.edu',
                  'https://vizier.idia.ac.za', 'https://vizier.nao.ac.jp']
TAP_GAIA = 'https://gea.esac.esa.int/tap-server/tap/sync'
SKYBOT = 'https://ssp.imcce.fr/webservices/skybot/api/conesearch.php'

NOMBRES = {'gaia': 'Gaia DR3', 'apass': 'APASS DR9', 'vsx': 'AAVSO VSX', 'skybot': 'SkyBoT', 'ngc': 'NGC/IC 2000'}


class ErrorCatalogo(Exception):
    """Error al descargar o leer un catálogo (el mensaje se traduce en Emillium.py)."""

    def __init__(self, clave, **datos):
        self.clave, self.datos = clave, datos
        super().__init__(clave)


# =============================================================================
# Carpeta de la caché
# =============================================================================


def carpeta_cache():
    base = os.environ.get('APPDATA') or os.path.join(os.path.expanduser('~'), '.config')
    carpeta = os.path.join(base, 'Emillium', 'catalogos')
    os.makedirs(carpeta, exist_ok=True)
    return carpeta


def _guardar(datos):
    """Guarda un trozo de catálogo descargado (o importado) y devuelve su ruta."""
    nombre = '{}_{:.4f}_{:+.4f}_{:.3f}_{}.json'.format(datos['catalogo'], datos['ar'], datos['dec'],
                                                       datos['radio'], int(time.time()))
    nombre = re.sub(r'[^A-Za-z0-9_.+-]', '_', nombre)
    ruta = os.path.join(carpeta_cache(), nombre)
    temporal = ruta + '.tmp'
    with open(temporal, 'w', encoding='utf-8') as archivo:
        json.dump(datos, archivo, ensure_ascii=False)
    os.replace(temporal, ruta)
    _escribir_ficha(ruta, datos)
    return ruta


def _ruta_ficha(ruta):
    return ruta[:-5] + '.meta'


def _escribir_ficha(ruta, datos):
    """Ficha pequeña con lo que hay en un trozo de catálogo (sin las estrellas): así, para
    saber qué zonas hay guardadas, no hace falta leer los archivos grandes."""
    ficha = {k: v for k, v in datos.items() if k != 'filas'} | {'n': len(datos['filas'])}
    with open(_ruta_ficha(ruta), 'w', encoding='utf-8') as archivo:
        json.dump(ficha, archivo, ensure_ascii=False)
    return ficha


_MEMORIA = {}          # últimos trozos leídos: ruta -> (fecha de modificación, datos)


def _cargar(ruta):
    """Lee un trozo de catálogo (guardando en memoria los últimos para no releerlos)."""
    momento = os.path.getmtime(ruta)
    guardado = _MEMORIA.get(ruta)
    if guardado and guardado[0] == momento:
        return guardado[1]
    with open(ruta, encoding='utf-8') as archivo:
        datos = json.load(archivo)
    if len(_MEMORIA) >= 6:
        _MEMORIA.pop(next(iter(_MEMORIA)))
    _MEMORIA[ruta] = (momento, datos)
    return datos


def listar_cache():
    """Trozos de catálogo guardados: lista de (ruta, metadatos sin las filas)."""
    resultado = []
    for nombre in sorted(os.listdir(carpeta_cache())):
        if not nombre.endswith('.json'):
            continue
        ruta = os.path.join(carpeta_cache(), nombre)
        try:
            ficha_ruta = _ruta_ficha(ruta)
            if os.path.exists(ficha_ruta) and os.path.getmtime(ficha_ruta) >= os.path.getmtime(ruta):
                with open(ficha_ruta, encoding='utf-8') as archivo:
                    ficha = json.load(archivo)
            else:      # archivos guardados antes de que existieran las fichas
                ficha = _escribir_ficha(ruta, _cargar(ruta))
            resultado.append((ruta, ficha))
        except (OSError, ValueError, KeyError):
            continue
    return resultado


def borrar_de_cache(ruta):
    """Borra un trozo de catálogo guardado. Por seguridad, solo archivos .json de la
    carpeta de la caché (nunca otra cosa, aunque se pida)."""
    carpeta = os.path.realpath(carpeta_cache())
    ruta = os.path.realpath(ruta)
    if os.path.dirname(ruta) != carpeta or not ruta.endswith('.json'):
        raise ValueError('ruta fuera de la caché')
    os.remove(ruta)
    _MEMORIA.pop(ruta, None)
    if os.path.exists(_ruta_ficha(ruta)):
        os.remove(_ruta_ficha(ruta))


def tamano_cache():
    """Espacio que ocupan los catálogos guardados (bytes)."""
    total = 0
    for nombre in os.listdir(carpeta_cache()):
        if nombre.endswith(('.json', '.meta')):
            total += os.path.getsize(os.path.join(carpeta_cache(), nombre))
    return total


def _distancia_grados(ar1, dec1, ar2, dec2):
    a1, d1, a2, d2 = (math.radians(v) for v in (ar1, dec1, ar2, dec2))
    c = math.sin(d1) * math.sin(d2) + math.cos(d1) * math.cos(d2) * math.cos(a1 - a2)
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


def desde_cache(catalogo, ar, dec, radio, jd=None):
    """Estrellas del catálogo en el círculo (ar, dec, radio) usando SOLO lo guardado
    en local. Sirve un trozo guardado que cubra todo el círculo (para SkyBoT, además,
    de un instante a menos de 30 minutos). Devuelve (filas, metadatos) o None."""
    mejor = None
    for ruta, meta in listar_cache():
        if meta.get('catalogo') != catalogo and not (catalogo == 'local' and str(meta.get('catalogo', '')).startswith('local')):
            continue
        if catalogo == 'skybot' and (jd is None or abs(meta.get('jd', 0) - jd) > 30 / 1440):
            continue
        d = _distancia_grados(ar, dec, meta['ar'], meta['dec'])
        if d + radio <= meta['radio'] * 1.0001 or (str(meta.get('catalogo', '')).startswith('local') and d <= meta['radio']):
            if mejor is None or meta.get('fecha', 0) > mejor[1].get('fecha', 0):
                mejor = (ruta, meta)
    if mejor is None:
        return None
    datos = _cargar(mejor[0])
    filas = [f for f in datos['filas'] if _distancia_grados(ar, dec, f['ar'], f['dec']) <= radio]
    return filas, {k: v for k, v in datos.items() if k != 'filas'}


# =============================================================================
# Descargas (SOLO se llaman desde un botón)
# =============================================================================


def _pedir(url, datos=None, tiempo=60):
    peticion = urllib.request.Request(url, data=datos, headers={'User-Agent': AGENTE})
    with urllib.request.urlopen(peticion, timeout=tiempo) as respuesta:
        return respuesta.read().decode('utf-8', 'replace')


def _leer_tsv_vizier(texto):
    """Tabla TSV de VizieR -> lista de diccionarios (columna -> texto)."""
    lineas = [l for l in texto.splitlines() if l.strip() and not l.startswith('#')]
    if len(lineas) < 3:
        return []
    columnas = [c.strip() for c in lineas[0].split('\t')]
    filas = []
    for linea in lineas[3:]:          # detrás: unidades y una línea de guiones
        valores = [v.strip() for v in linea.split('\t')]
        if len(valores) == len(columnas):
            filas.append(dict(zip(columnas, valores)))
    return filas


def _vizier(fuente, ar, dec, radio, columnas, orden, max_filas, extra=None):
    parametros = {'-source': fuente, '-c': f'{ar:.6f} {dec:+.6f}', '-c.rd': f'{radio:.5f}',
                  '-out': ','.join(columnas), '-out.max': str(max_filas), '-sort': orden}
    parametros.update(extra or {})
    ultimo_error = None
    for espejo in ESPEJOS_VIZIER:
        try:
            filas = _leer_tsv_vizier(_pedir(espejo + '/viz-bin/asu-tsv?' + urllib.parse.urlencode(parametros), tiempo=45))
        except Exception as error:   # este espejo no responde: se prueba el siguiente
            ultimo_error = error
            continue
        if filas:
            return filas, espejo
    if ultimo_error is not None:
        raise ErrorCatalogo('sin_conexion', detalle=str(ultimo_error)[:200])
    return [], ESPEJOS_VIZIER[-1]


def _num(texto):
    try:
        return float(texto)
    except (TypeError, ValueError):
        return None


def descargar(catalogo, ar, dec, radio, max_filas=5000, jd=None, codigo_observatorio='500'):
    """Descarga las estrellas (u objetos) del catálogo en un círculo y las guarda en la
    caché. ESTA FUNCIÓN SE CONECTA A INTERNET: llamarla solo por una acción del usuario.
    Devuelve (filas, ruta del archivo guardado)."""
    radio = min(max(radio, 0.01), 3.0)
    if catalogo == 'gaia':
        filas, fuente = _descargar_gaia(ar, dec, radio, max_filas)
        epoca = 2016.0
    elif catalogo == 'apass':
        crudas, fuente = _vizier('II/336/apass9', ar, dec, radio,
                                 ['RAJ2000', 'DEJ2000', 'Vmag', 'e_Vmag', 'Bmag', 'e_Bmag', "g'mag", "r'mag", "i'mag"],
                                 'Vmag', max_filas)
        filas = [{'id': f'APASS {i + 1}', 'ar': _num(c['RAJ2000']), 'dec': _num(c['DEJ2000']), 'mag': _num(c['Vmag']),
                  'V': _num(c['Vmag']), 'e_V': _num(c['e_Vmag']), 'B': _num(c['Bmag']), 'e_B': _num(c['e_Bmag']),
                  'g': _num(c["g'mag"]), 'r': _num(c["r'mag"]), 'i': _num(c["i'mag"])} for i, c in enumerate(crudas)]
        epoca = 2000.0
    elif catalogo == 'vsx':
        crudas, fuente = _vizier('B/vsx/vsx', ar, dec, radio, ['Name', 'RAJ2000', 'DEJ2000', 'Type', 'max', 'min', 'Period'],
                                 'max', max_filas)
        filas = [{'id': c['Name'], 'ar': _num(c['RAJ2000']), 'dec': _num(c['DEJ2000']), 'mag': _num(c['max']),
                  'tipo': c['Type'], 'max': _num(c['max']), 'min': _num(c['min']), 'periodo': _num(c['Period'])}
                 for c in crudas]
        epoca = 2000.0
    elif catalogo == 'ngc':
        # NGC 2000.0: galaxias, nebulosas y cúmulos (NGC e IC). Se piden también las
        # coordenadas calculadas en grados (_RAJ2000, _DEJ2000).
        crudas, fuente = _vizier('VII/118/ngc2000', ar, dec, radio, ['Name', 'Type', 'size', 'mag'], 'mag', max_filas,
                                 extra={'-out.add': '_RAJ,_DEJ', '-oc.form': 'd'})
        filas = []
        for c in crudas:
            nombre = c['Name'].strip()
            nombre = ('IC ' + nombre[1:].strip()) if nombre.startswith('I') else ('NGC ' + nombre)
            filas.append({'id': nombre, 'ar': _num(c.get('_RAJ2000')), 'dec': _num(c.get('_DEJ2000')),
                          'tipo': c.get('Type', '').strip(), 'tamano': _num(c.get('size')), 'mag': _num(c.get('mag'))})
        epoca = 2000.0
    elif catalogo == 'skybot':
        if jd is None:
            raise ErrorCatalogo('skybot_sin_fecha')
        filas, fuente = _descargar_skybot(jd, ar, dec, radio, codigo_observatorio)
        epoca = None
    else:
        raise ValueError(catalogo)
    filas = [f for f in filas if f['ar'] is not None and f['dec'] is not None]
    datos = {'catalogo': catalogo, 'nombre': NOMBRES[catalogo], 'ar': ar, 'dec': dec, 'radio': radio,
             'fecha': time.time(), 'fuente': fuente, 'epoca': epoca, 'filas': filas}
    if catalogo == 'skybot':
        datos['jd'] = jd
    return filas, _guardar(datos)


def _descargar_gaia(ar, dec, radio, max_filas):
    """Gaia DR3: primero el archivo oficial de la ESA; si no responde, VizieR."""
    adql = (f'SELECT TOP {int(max_filas)} source_id, ra, dec, ra_error, dec_error, pmra, pmdec, parallax, '
            f'phot_g_mean_mag, bp_rp FROM gaiadr3.gaia_source '
            f'WHERE 1=CONTAINS(POINT(ra, dec), CIRCLE({ar:.6f}, {dec:.6f}, {radio:.5f})) '
            f'AND phot_g_mean_mag IS NOT NULL ORDER BY phot_g_mean_mag')
    try:
        texto = _pedir(TAP_GAIA, urllib.parse.urlencode({'REQUEST': 'doQuery', 'LANG': 'ADQL', 'FORMAT': 'csv',
                                                         'QUERY': adql}).encode(), tiempo=90)
        filas = []
        for c in csv.DictReader(io.StringIO(texto)):
            filas.append({'id': f"Gaia DR3 {c['source_id']}", 'ar': _num(c['ra']), 'dec': _num(c['dec']),
                          'mag': _num(c['phot_g_mean_mag']), 'e_ar_mas': _num(c['ra_error']),
                          'e_dec_mas': _num(c['dec_error']), 'pmra': _num(c['pmra']), 'pmdec': _num(c['pmdec']),
                          'plx': _num(c['parallax']), 'color': _num(c['bp_rp'])})
        if filas or 'source_id' in texto:
            return filas, TAP_GAIA
    except Exception:
        pass   # el archivo de la ESA no responde: se prueba VizieR
    crudas, fuente = _vizier('I/355/gaiadr3', ar, dec, radio,
                             ['Source', 'RA_ICRS', 'DE_ICRS', 'e_RA_ICRS', 'e_DE_ICRS', 'pmRA', 'pmDE', 'Plx', 'Gmag', 'BP-RP'],
                             'Gmag', max_filas)
    return [{'id': f"Gaia DR3 {c['Source']}", 'ar': _num(c['RA_ICRS']), 'dec': _num(c['DE_ICRS']),
             'mag': _num(c['Gmag']), 'e_ar_mas': _num(c['e_RA_ICRS']), 'e_dec_mas': _num(c['e_DE_ICRS']),
             'pmra': _num(c['pmRA']), 'pmdec': _num(c['pmDE']), 'plx': _num(c['Plx']), 'color': _num(c['BP-RP'])}
            for c in crudas], fuente


def _sexagesimal(texto, horas):
    partes = texto.replace(':', ' ').split()
    signo = -1 if partes[0].startswith('-') else 1
    g, m, s = abs(float(partes[0])), float(partes[1]), float(partes[2])
    valor = signo * (g + m / 60 + s / 3600)
    return valor * 15 if horas else valor


def _descargar_skybot(jd, ar, dec, radio, codigo):
    parametros = {'-ep': f'{jd:.6f}', '-ra': f'{ar:.6f}', '-dec': f'{dec:.6f}', '-rd': f'{radio:.5f}',
                  '-mime': 'text', '-output': 'object', '-loc': codigo or '500'}
    try:
        texto = _pedir(SKYBOT + '?' + urllib.parse.urlencode(parametros), tiempo=90)
    except Exception as error:
        raise ErrorCatalogo('sin_conexion', detalle=str(error)[:200])
    filas = []
    for linea in texto.splitlines():
        if linea.startswith('#') or '|' not in linea:
            continue
        p = [c.strip() for c in linea.split('|')]
        if len(p) < 8:
            continue
        try:
            filas.append({'id': p[1] if p[0] in ('-', '') else f'({p[0]}) {p[1]}', 'numero': p[0], 'nombre': p[1],
                          'ar': _sexagesimal(p[2], True), 'dec': _sexagesimal(p[3], False), 'clase': p[4],
                          'mag': _num(p[5]), 'error_pos': _num(p[6]), 'distancia': _num(p[7])})
        except (ValueError, IndexError):
            continue
    return filas, SKYBOT


# =============================================================================
# Catálogos locales (CSV) - sin conexión
# =============================================================================

_NOMBRES_AR = ('ra', 'ra_deg', 'raj2000', 'ra_icrs', 'ar', 'alpha', 'radeg')
_NOMBRES_DEC = ('dec', 'dec_deg', 'dej2000', 'de_icrs', 'de', 'delta', 'decdeg')
_NOMBRES_MAG = ('mag', 'gmag', 'vmag', 'rmag', 'magnitude', 'phot_g_mean_mag', 'v')
_NOMBRES_ID = ('id', 'name', 'nom', 'nombre', 'source', 'source_id', 'designation')


def importar_csv(ruta):
    """Importa un catálogo propio (CSV con columnas de AR y Dec en grados y, si se
    quiere, magnitud y nombre). Se detecta el separador. Se guarda en la caché como
    catálogo 'local' para usarlo siempre. Devuelve (número de estrellas, ruta guardada)."""
    with open(ruta, encoding='utf-8-sig', errors='replace') as archivo:
        texto = archivo.read()
    try:
        dialecto = csv.Sniffer().sniff(texto[:4096], delimiters=',;\t|')
    except csv.Error:
        dialecto = csv.excel
    lector = csv.reader(io.StringIO(texto), dialecto)
    filas_texto = [f for f in lector if f and any(c.strip() for c in f)]
    if len(filas_texto) < 2:
        raise ErrorCatalogo('csv_vacio')
    cabecera = [c.strip().lower() for c in filas_texto[0]]

    def columna(nombres):
        return next((i for i, c in enumerate(cabecera) if c in nombres), None)
    i_ar, i_dec, i_mag, i_id = columna(_NOMBRES_AR), columna(_NOMBRES_DEC), columna(_NOMBRES_MAG), columna(_NOMBRES_ID)
    if i_ar is None or i_dec is None:
        raise ErrorCatalogo('csv_columnas', columnas=', '.join(filas_texto[0][:8]))
    filas = []
    for n, f in enumerate(filas_texto[1:], 1):
        try:
            ar_v = float(f[i_ar].replace(',', '.')) if dialecto.delimiter != ',' else float(f[i_ar])
            dec_v = float(f[i_dec].replace(',', '.')) if dialecto.delimiter != ',' else float(f[i_dec])
        except (ValueError, IndexError):
            continue
        if not (0 <= ar_v <= 360 and -90 <= dec_v <= 90):
            continue
        mag = None
        if i_mag is not None and i_mag < len(f):
            mag = _num(f[i_mag].replace(',', '.'))
        filas.append({'id': f[i_id].strip() if i_id is not None and i_id < len(f) else f'#{n}', 'ar': ar_v,
                      'dec': dec_v, 'mag': mag})
    if not filas:
        raise ErrorCatalogo('csv_vacio')
    # círculo que contiene todas las estrellas
    x = sum(math.cos(math.radians(f['dec'])) * math.cos(math.radians(f['ar'])) for f in filas)
    y = sum(math.cos(math.radians(f['dec'])) * math.sin(math.radians(f['ar'])) for f in filas)
    z = sum(math.sin(math.radians(f['dec'])) for f in filas)
    ar_c = math.degrees(math.atan2(y, x)) % 360
    dec_c = math.degrees(math.atan2(z, math.hypot(x, y)))
    radio = max(_distancia_grados(ar_c, dec_c, f['ar'], f['dec']) for f in filas) + 1e-4
    nombre = os.path.basename(ruta)
    datos = {'catalogo': 'local', 'nombre': nombre, 'ar': ar_c, 'dec': dec_c, 'radio': radio, 'fecha': time.time(),
             'fuente': ruta, 'epoca': None, 'filas': filas}
    return len(filas), _guardar(datos)
