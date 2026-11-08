"""
Emillium v1.5 (Expanded) - Astrometría y fotometría en imágenes FITS
====================================================================

Programa de escritorio para abrir imágenes astronómicas en formato FITS,
girarlas o invertirlas como un espejo, y medir la posición de las estrellas
en PÍXELES con precisión subpíxel (centroide), su brillo (fotometría de
apertura) y, si la foto tiene solución astrométrica, sus coordenadas celestes.

Archivos del programa:
  - Emillium.py   la ventana y todo lo de la versión anterior (este archivo).
  - ciencia.py    cálculos nuevos: fotometría de apertura, detección de
                  estrellas, estrellas dobles, ajuste WCS, objetos móviles...
  - ades.py       informes ADES para el Minor Planet Center y su validación.
  - catalogos.py  Gaia DR3, APASS, VSX y SkyBoT con caché local permanente.

Solo utiliza:
  - La biblioteca estándar de Python (tkinter para la ventana, csv, math...).
  - numpy, para hacer los cálculos con matrices de forma rápida.
  - matplotlib, solo para dibujar los gráficos (curvas de luz, residuos...).

No usa astropy ni photutils: el lector de archivos FITS, el giro de la imagen,
el cálculo del centroide, la fotometría y las coordenadas están programados
aquí mismo, paso a paso.

Funciona SIN INTERNET. Solo se conecta cuando el usuario pulsa un botón de
descarga de catálogos (y lo descargado se guarda en local para siempre).

Convenio de coordenadas
-----------------------
  - X crece hacia la derecha e Y hacia arriba (como en astronomía).
  - El píxel (0, 0) es el de la esquina inferior izquierda y su centro tiene
    coordenadas exactamente (0, 0). Por eso un valor como X = 120,37 quiere
    decir "un poco a la derecha del centro del píxel 120".
  - Las coordenadas se dan siempre en la imagen TAL COMO SE VE, es decir,
    después de aplicar el giro y el espejo.
"""

import csv
import json
import math
import os
import base64
import random
import re
import struct
import tempfile
import threading
import time
import wave
import zipfile
import zlib
import sys
import tkinter as tk
from datetime import datetime
from tkinter import filedialog, messagebox, ttk
from xml.sax.saxutils import escape

import numpy as np

import ades
import catalogos
import ciencia as C

VERSION = '1.5'
TITULO_VENTANA = 'Emillium'  # texto de la barra de la ventana (siempre el mismo)


class ErrorFits(ValueError):
    """Error al leer o transformar una foto. Guarda la clave del mensaje para
    poder mostrarlo traducido (ver TEXTOS más abajo)."""

    def __init__(self, clave, **datos):
        self.clave, self.datos = clave, datos
        super().__init__(clave)


# =============================================================================
# 1. LECTOR DE ARCHIVOS FITS
# =============================================================================
# Un archivo FITS está formado por uno o varios "HDU" (Header Data Unit).
# Cada HDU tiene:
#   - Una CABECERA de texto en bloques de 2880 bytes. Cada bloque contiene
#     36 "tarjetas" de 80 caracteres del tipo:  NAXIS1  =   1024 / comentario
#     La cabecera termina con la tarjeta END.
#   - Los DATOS de la imagen en binario (big-endian), también rellenados
#     hasta un múltiplo de 2880 bytes.
# Las palabras clave importantes son:
#   BITPIX : tipo de número de cada píxel (8, 16, 32, 64 enteros; -32, -64 reales)
#   NAXIS  : número de dimensiones; NAXIS1 = ancho, NAXIS2 = alto
#   BZERO, BSCALE : valor_real = BZERO + BSCALE * valor_guardado
#   BLANK  : valor entero que indica "píxel sin dato"

TAMANO_BLOQUE = 2880
TAMANO_TARJETA = 80

TIPOS_BITPIX = {
    8: '>u1',    # entero sin signo de 8 bits
    16: '>i2',   # entero de 16 bits
    32: '>i4',   # entero de 32 bits
    64: '>i8',   # entero de 64 bits
    -32: '>f4',  # real de 32 bits
    -64: '>f8',  # real de 64 bits
}


def _interpretar_valor(texto):
    """Convierte el valor de una tarjeta FITS (texto) en un dato de Python."""
    texto = texto.strip()
    if texto.startswith("'"):
        # Cadena de texto entre comillas simples ('' dentro = una comilla).
        resultado = []
        i = 1
        while i < len(texto):
            if texto[i] == "'":
                if i + 1 < len(texto) and texto[i + 1] == "'":
                    resultado.append("'")
                    i += 2
                    continue
                break
            resultado.append(texto[i])
            i += 1
        return ''.join(resultado).rstrip()

    valor = texto.split('/')[0].strip()  # quitar el comentario
    if valor == 'T':
        return True
    if valor == 'F':
        return False
    try:
        return int(valor)
    except ValueError:
        pass
    try:
        return float(valor.replace('D', 'E'))
    except ValueError:
        return valor


def _leer_cabecera(archivo):
    """Lee una cabecera FITS completa (hasta END). Devuelve None si el archivo
    ya no tiene más datos."""
    cabecera = {}
    primer_bloque = True
    while True:
        bloque = archivo.read(TAMANO_BLOQUE)
        if primer_bloque and len(bloque) == 0:
            return None
        if len(bloque) < TAMANO_BLOQUE:
            raise ErrorFits('err_cabecera')
        primer_bloque = False
        texto = bloque.decode('ascii', errors='replace')
        for i in range(0, TAMANO_BLOQUE, TAMANO_TARJETA):
            tarjeta = texto[i:i + TAMANO_TARJETA]
            clave = tarjeta[:8].strip()
            if clave == 'END':
                return cabecera
            if clave and tarjeta[8] == '=' and clave not in cabecera:
                cabecera[clave] = _interpretar_valor(tarjeta[10:])


# Palabras que describen la estructura de cada HDU: nunca se heredan de la cabecera principal.
_ESTRUCTURALES = re.compile(r'^(SIMPLE|XTENSION|BITPIX|NAXIS\d*|EXTEND|PCOUNT|GCOUNT|GROUPS|BZERO|BSCALE|BLANK|'
                            r'TFIELDS|TFORM\d+|TTYPE\d+|TUNIT\d+|TDIM\d+|THEAP|Z[A-Z0-9_]*|EXTNAME|EXTVER|'
                            r'CHECKSUM|DATASUM|END)$')


def _heredar_cabecera(primaria, extension):
    """Muchos FITS (HST, por ejemplo) guardan la imagen en una extensión y dejan la
    fecha, la exposición o el filtro en la cabecera principal. Como dice la norma
    (INHERIT), la extensión se completa con lo que le falte de la principal."""
    if extension.get('INHERIT') is False:
        return extension
    combinada = {k: v for k, v in primaria.items() if not _ESTRUCTURALES.match(k)}
    combinada.update(extension)
    return combinada


def separar_ruta(ruta):
    """Una foto se identifica por su archivo y, si el archivo tiene varias imágenes
    (extensiones o planos), por el número de imagen: 'foto.fits::2' es la tercera
    (se cuenta desde 0; la primera es simplemente 'foto.fits'). Devuelve (archivo, número)."""
    m = re.match(r'^(.*)::(\d+)$', ruta)
    return (m.group(1), int(m.group(2))) if m else (ruta, 0)


def unir_ruta(archivo, indice):
    return archivo if indice == 0 else f'{archivo}::{indice}'


_CLAVES_WCS = re.compile(r'^(CTYPE\d|CRPIX\d|CRVAL\d|CDELT\d|CROTA\d|CD\d_\d|PC\d_\d|CUNIT\d|[AB]P?_\d+_\d+|[AB]P?_ORDER|'
                         r'WCSAXES|LONPOLE|LATPOLE)$')


LIMITE_GZ = 1_600_000_000   # un .fits.gz no puede descomprimirse en más de 1,6 GB (protección)


def _abrir_fits_crudo(ruta):
    """Abre un FITS para leerlo. Si es un .fits.gz (todo el archivo comprimido con gzip,
    como los de muchos archivos astronómicos), se descomprime en memoria, a trozos y con
    un límite de tamaño (un archivo manipulado podría ocupar muchísimo al descomprimirse)."""
    archivo = open(ruta, 'rb')
    if archivo.read(2) != b'\x1f\x8b':
        archivo.seek(0)
        return archivo
    archivo.close()
    import gzip
    import io
    datos = bytearray()
    try:
        with gzip.open(ruta, 'rb') as comprimido:
            while True:
                trozo = comprimido.read(16 * 1024 * 1024)
                if not trozo:
                    break
                datos += trozo
                if len(datos) > LIMITE_GZ:
                    raise ErrorFits('err_enorme')
    except (OSError, EOFError, zlib.error):
        raise ErrorFits('err_gz')        # .gz roto o que no es gzip de verdad
    return io.BytesIO(bytes(datos))


def leer_fits(ruta, indice=0):
    """Lee una imagen 2D de un archivo FITS: la número 'indice' (desde 0) entre todas
    las que tiene (cada extensión de imagen y cada plano de un cubo cuentan como una;
    una imagen en color de 3 planos cuenta como una sola).

    Devuelve (imagen, cabecera, nota):
      imagen   -> matriz numpy float32 con forma (alto, ancho); la fila 0 es
                  la de abajo.
      cabecera -> diccionario con las palabras clave. Además, '_N_IMAGENES' dice
                  cuántas imágenes tiene el archivo.
      nota     -> texto con avisos (por ejemplo, si era una imagen en color).
    """
    resultado = None
    contadas = 0
    with _abrir_fits_crudo(ruta) as archivo:
        archivo.seek(0, 2)
        tamano_archivo = archivo.tell()
        archivo.seek(0)
        if archivo.read(6) != b'SIMPLE':
            raise ErrorFits('err_no_fits')
        archivo.seek(0)

        comprimido = False
        primaria = None   # cabecera del primer HDU: las extensiones "heredan" sus datos (fecha, exposición...)
        while True:
            cabecera = _leer_cabecera(archivo)
            if cabecera is None:
                break
            if primaria is None:
                primaria = cabecera
            else:
                cabecera = _heredar_cabecera(primaria, cabecera)

            # Se comprueba que la cabecera tenga sentido antes de creérsela: un archivo
            # manipulado podría pedir tamaños negativos o enormes y colgar el programa.
            try:
                bitpix = int(cabecera.get('BITPIX', 0))
                naxis = int(cabecera.get('NAXIS', 0))
                ejes = [int(cabecera.get(f'NAXIS{i}', 0)) for i in range(1, min(naxis, 999) + 1)]
                pcount = int(cabecera.get('PCOUNT', 0))
                gcount = int(cabecera.get('GCOUNT', 1))
            except (TypeError, ValueError):
                raise ErrorFits('err_cabecera_rara')
            if not 0 <= naxis <= 999 or any(e < 0 for e in ejes) or pcount < 0 or gcount < 0:
                raise ErrorFits('err_cabecera_rara')
            num_valores = math.prod(ejes) if naxis > 0 else 0   # entero de Python: no se desborda
            bytes_por_valor = abs(bitpix) // 8
            bytes_datos = bytes_por_valor * gcount * (pcount + num_valores) if naxis > 0 else 0
            if bytes_datos > tamano_archivo:  # dice tener más datos de los que caben en el archivo
                raise ErrorFits('err_datos')
            relleno = (-bytes_datos) % TAMANO_BLOQUE
            inicio_datos = archivo.tell()

            if cabecera.get('ZIMAGE') is True and str(cabecera.get('XTENSION', '')).strip() == 'BINTABLE':
                # Imagen comprimida por "teselas" (.fz): se descomprime aquí mismo.
                comprimido = True
                if resultado is None and indice == contadas:
                    tabla = archivo.read(bytes_datos)
                    if len(tabla) < bytes_datos:
                        raise ErrorFits('err_datos')
                    datos, bitpix_z = descomprimir_fz(cabecera, tabla)
                    imagen, nota = _convertir_datos(datos, cabecera, bitpix_z, blank=cabecera.get('ZBLANK'))
                    resultado = (imagen, dict(cabecera), nota)
                contadas += 1
                archivo.seek(inicio_datos + bytes_datos + relleno)
                continue

            extension = cabecera.get('XTENSION')
            es_imagen = extension is None or str(extension).strip() == 'IMAGE'
            if es_imagen and naxis >= 2 and num_valores > 0:
                if bitpix not in TIPOS_BITPIX:
                    raise ErrorFits('err_bitpix', bitpix=bitpix)
                planos = math.prod(e for e in ejes[2:] if e > 1) if naxis > 2 else 1
                color = naxis >= 3 and planos == 3
                aqui = 1 if color else planos
                if resultado is None and contadas <= indice < contadas + aqui:
                    if num_valores > 400_000_000 and not planos > 1:  # más de 400 megapíxeles: no cabría
                        raise ErrorFits('err_enorme')
                    plano_pixeles = ejes[0] * ejes[1]
                    if plano_pixeles > 400_000_000:
                        raise ErrorFits('err_enorme')
                    if color:
                        cuantos, desde = num_valores, 0
                    else:
                        cuantos, desde = plano_pixeles, (indice - contadas) * plano_pixeles
                    archivo.seek(inicio_datos + desde * bytes_por_valor)
                    crudo = archivo.read(bytes_por_valor * cuantos)
                    if len(crudo) < bytes_por_valor * cuantos:
                        raise ErrorFits('err_datos')
                    datos = np.frombuffer(crudo, dtype=TIPOS_BITPIX[bitpix])
                    datos = datos.reshape(ejes[::-1] if color else (ejes[1], ejes[0]))   # numpy: (alto, ancho)
                    imagen, nota = _convertir_datos(datos, cabecera, bitpix)
                    cab = dict(cabecera)
                    if not color and indice - contadas > 0:
                        # los otros planos de un cubo suelen ser otras cámaras o filtros: la
                        # solución astrométrica de la cabecera (la del primer plano) no les sirve
                        cab = {k: v for k, v in cab.items() if not _CLAVES_WCS.match(k)}
                        nota = 'nota_plano_sin_wcs'
                    resultado = (imagen, cab, nota)
                contadas += aqui
            # se salta lo que quede de los datos de este HDU y se sigue
            archivo.seek(inicio_datos + bytes_datos + relleno)

    if resultado is not None:
        resultado[1]['_N_IMAGENES'] = contadas
        return resultado
    if contadas:
        raise ErrorFits('err_indice', n=contadas)
    if comprimido:
        raise ErrorFits('err_comprimido')
    raise ErrorFits('err_sin_imagen')

# -----------------------------------------------------------------------------
# FITS comprimidos (.fz): la imagen se guarda partida en "teselas" (normalmente
# una fila cada una), cada tesela comprimida por separado dentro de una tabla
# binaria. Aquí se descomprimen RICE_1 (el más habitual), GZIP_1 y GZIP_2,
# incluidas las imágenes de coma flotante "cuantizadas" (con o sin dithering).
# -----------------------------------------------------------------------------
_BYTES_TFORM = {'L': 1, 'B': 1, 'I': 2, 'J': 4, 'K': 8, 'A': 1, 'E': 4, 'D': 8, 'C': 8, 'M': 16, 'P': 8, 'Q': 16}
_TIPOS_TFORM = {'B': '>u1', 'I': '>i2', 'J': '>i4', 'K': '>i8', 'E': '>f4', 'D': '>f8'}
_TABLA_ALEATORIA = None


def _tabla_aleatoria():
    """Los 10.000 números pseudoaleatorios del estándar FITS para el "dithering"
    (generador de Park y Miller), iguales en todos los programas."""
    global _TABLA_ALEATORIA
    if _TABLA_ALEATORIA is None:
        a, m, semilla = 16807.0, 2147483647.0, 1.0
        valores = np.empty(10000)
        for i in range(10000):
            temp = a * semilla
            semilla = temp - m * int(temp / m)
            valores[i] = semilla / m
        _TABLA_ALEATORIA = valores
    return _TABLA_ALEATORIA


def _columnas_tabla(cab):
    """Posición, repetición y tipo de cada columna de una tabla binaria FITS."""
    columnas, posicion = {}, 0
    for i in range(1, int(cab.get('TFIELDS', 0)) + 1):
        forma = str(cab.get(f'TFORM{i}', '')).strip().upper()
        nombre = str(cab.get(f'TTYPE{i}', f'COL{i}')).strip().upper()
        trozo = re.match(r'(\d*)([LXBIJKAEDCMPQ])([A-Z]?)', forma)
        if not trozo:
            raise ErrorFits('err_cabecera_rara')
        repeticion = int(trozo.group(1)) if trozo.group(1) else 1
        letra = trozo.group(2)
        tamano = (repeticion + 7) // 8 if letra == 'X' else repeticion * _BYTES_TFORM[letra]
        columnas[nombre] = (posicion, repeticion, letra, trozo.group(3))
        posicion += tamano
    return columnas


def _rice_decodificar(teselas, n, bytes_pixel, bloque):
    """Descomprime varias teselas RICE de n píxeles a la vez.

    El código RICE guarda cada píxel como la diferencia con el anterior: por cada
    bloque de 32 píxeles se elige un número k de bits; cada diferencia se escribe
    como (diferencia >> k) en "unario" (tantos ceros como indique, y un 1) seguida
    de los k bits bajos. En lugar de ir bit a bit, se avanza a la vez en todas las
    teselas con numpy (muy rápido)."""
    bits_fs, fs_max, bits_valor = {1: (3, 6, 8), 2: (4, 14, 16), 4: (5, 25, 32)}[bytes_pixel]
    mascara = (1 << bits_valor) - 1
    total = len(teselas)
    primero = np.array([int.from_bytes(t[:bytes_pixel], 'big') for t in teselas], dtype=np.int64)
    flujos = [t[bytes_pixel:] for t in teselas]
    inicios = np.cumsum([0] + [len(f) * 8 for f in flujos[:-1]]).astype(np.int64)
    bits = np.unpackbits(np.frombuffer(b''.join(flujos) + b'\0' * 16, dtype=np.uint8))
    indices = np.where(bits == 1, np.arange(bits.size, dtype=np.int64), bits.size - 1)
    siguiente_uno = np.minimum.accumulate(indices[::-1])[::-1]       # posición del próximo bit 1
    columnas = np.arange(32)

    def leer(posicion, k, k_max):
        """Lee k bits (k distinto en cada tesela) a partir de posicion."""
        if k_max <= 0:
            return np.zeros(total, dtype=np.int64)
        j = columnas[:k_max]
        trozo = bits[np.minimum(posicion[:, None] + j, bits.size - 1)].astype(np.int64)
        pesos = np.where(j < k[:, None], np.left_shift(1, np.maximum(k[:, None] - 1 - j, 0)), 0)
        return (trozo * pesos).sum(axis=1)

    posicion = inicios.copy()
    ultimo = primero.copy()
    salida = np.empty((total, n), dtype=np.int64)
    for i0 in range(0, n, bloque):
        fs = leer(posicion, np.full(total, bits_fs), bits_fs) - 1
        posicion = posicion + bits_fs
        cero, alto = fs < 0, fs == fs_max
        fs_normal = np.clip(fs, 0, fs_max - 1)
        hay_alto = bool(alto.any())
        k_max = int(fs_normal.max())
        for i in range(i0, min(i0 + bloque, n)):
            uno = siguiente_uno[posicion]
            pos_bajos = uno + 1
            diferencia = ((uno - posicion) << fs_normal) | leer(pos_bajos, fs_normal, k_max)
            nueva_posicion = pos_bajos + fs_normal
            if hay_alto:  # bloques de "alta entropía": cada diferencia va entera, sin comprimir
                cruda = leer(posicion, np.full(total, bits_valor), bits_valor)
                diferencia = np.where(alto, cruda, diferencia)
                nueva_posicion = np.where(alto, posicion + bits_valor, nueva_posicion)
            diferencia = np.where(cero, 0, diferencia)
            posicion = np.where(cero, posicion, nueva_posicion)
            # las diferencias negativas se guardan "en zigzag": 0, -1, 1, -2, 2...
            diferencia = np.where(diferencia & 1, ~(diferencia >> 1), diferencia >> 1)
            ultimo = (ultimo + diferencia) & mascara
            salida[:, i] = ultimo
    if bytes_pixel > 1:  # enteros con signo
        salida = np.where(salida >= 1 << (bits_valor - 1), salida - (1 << bits_valor), salida)
    return salida


def _desordenar_bytes(crudo, bytes_valor):
    """GZIP_2 guarda primero todos los primeros bytes, luego todos los segundos..."""
    n = len(crudo) // bytes_valor
    return np.frombuffer(crudo, dtype=np.uint8)[:n * bytes_valor].reshape(bytes_valor, n).T.tobytes()


def descomprimir_fz(cab, datos):
    """Reconstruye la imagen de un HDU comprimido (ZIMAGE = T).
    Devuelve (matriz numpy, BITPIX original)."""
    ancho_fila, filas = int(cab['NAXIS1']), int(cab['NAXIS2'])
    inicio_monton = int(cab.get('THEAP', ancho_fila * filas))
    columnas = _columnas_tabla(cab)
    tipo = str(cab.get('ZCMPTYPE', '')).strip().upper()
    if tipo not in ('RICE_1', 'GZIP_1', 'GZIP_2', 'NOCOMPRESS'):
        raise ErrorFits('err_fz_tipo', tipo=tipo or '?')
    zbitpix = int(cab['ZBITPIX'])
    ejes = [int(cab[f'ZNAXIS{i}']) for i in range(1, int(cab['ZNAXIS']) + 1)]
    if zbitpix not in TIPOS_BITPIX or any(e <= 0 for e in ejes) or math.prod(ejes) > 400_000_000:
        raise ErrorFits('err_cabecera_rara')
    tesela = [int(cab.get(f'ZTILE{i + 1}', ejes[0] if i == 0 else 1)) for i in range(len(ejes))]
    if any(t <= 0 for t in tesela):
        raise ErrorFits('err_cabecera_rara')
    parametros = {}
    for i in range(1, 20):
        if f'ZNAME{i}' in cab:
            parametros[str(cab[f'ZNAME{i}']).strip().upper()] = cab.get(f'ZVAL{i}')
    bloque = int(parametros.get('BLOCKSIZE', 32))
    cuantizada = zbitpix < 0 and ('ZSCALE' in columnas or 'ZSCALE' in cab)
    metodo = str(cab.get('ZQUANTIZ', 'NO_DITHER' if cuantizada else '')).strip().upper()
    if cuantizada and metodo not in ('NO_DITHER', 'SUBTRACTIVE_DITHER_1', 'SUBTRACTIVE_DITHER_2'):
        raise ErrorFits('err_fz_tipo', tipo=metodo)
    semilla = int(cab.get('ZDITHER0', 1))
    # tipo de los números guardados en cada tesela
    if cuantizada:
        tipo_guardado, bytes_pixel = '>i4', 4
    else:
        tipo_guardado = TIPOS_BITPIX[zbitpix]
        bytes_pixel = abs(zbitpix) // 8
    if tipo == 'RICE_1':
        bytes_pixel = int(parametros.get('BYTEPIX', 4 if cuantizada else bytes_pixel))
        if bytes_pixel not in (1, 2, 4):
            raise ErrorFits('err_fz_tipo', tipo=f'RICE_1 BYTEPIX={bytes_pixel}')

    def celda_variable(fila, nombre):
        if nombre not in columnas:
            return b'', None
        posicion, _rep, letra, subtipo = columnas[nombre]
        base = fila * ancho_fila + posicion
        if letra == 'P':
            n, desplazamiento = struct.unpack('>ii', datos[base:base + 8])
        elif letra == 'Q':
            n, desplazamiento = struct.unpack('>qq', datos[base:base + 16])
        else:
            raise ErrorFits('err_cabecera_rara')
        tam = _BYTES_TFORM.get(subtipo, 1)
        inicio = inicio_monton + desplazamiento
        if n < 0 or desplazamiento < 0 or inicio + n * tam > len(datos):
            raise ErrorFits('err_datos')
        return datos[inicio:inicio + n * tam], subtipo

    def celda_numero(fila, nombre, defecto):
        if nombre in columnas:
            posicion, _rep, letra, _sub = columnas[nombre]
            base = fila * ancho_fila + posicion
            return float(np.frombuffer(datos[base:base + _BYTES_TFORM[letra]], dtype=_TIPOS_TFORM[letra])[0])
        return float(cab.get(nombre, defecto))

    # Teselas: posición y tamaño de cada una (el primer eje varía más deprisa)
    cuantas = [math.ceil(e / t) for e, t in zip(ejes, tesela)]
    if math.prod(cuantas) != filas:
        raise ErrorFits('err_cabecera_rara')
    imagen = np.zeros(ejes[::-1], dtype=np.float64 if zbitpix < 0 else np.int64)
    pendientes_rice = {}   # n píxeles -> [(fila, trozos de la imagen, bytes)]
    for fila in range(filas):
        indice, cortes = fila, []
        for eje, (e, t, c) in enumerate(zip(ejes, tesela, cuantas)):
            k = indice % c
            indice //= c
            cortes.append(slice(k * t, min((k + 1) * t, e)))
        cortes = tuple(reversed(cortes))
        forma = tuple(s.stop - s.start for s in cortes)
        n = math.prod(forma)
        comprimido, _sub = celda_variable(fila, 'COMPRESSED_DATA')
        if not comprimido:  # tesela que no se pudo comprimir/cuantizar: va tal cual (o con gzip)
            crudo, subtipo = celda_variable(fila, 'GZIP_COMPRESSED_DATA')
            if crudo:
                crudo = zlib.decompress(crudo, 47)
                valores = np.frombuffer(crudo, dtype='>f4' if zbitpix == -32 else '>f8')
            else:
                crudo, subtipo = celda_variable(fila, 'UNCOMPRESSED_DATA')
                valores = np.frombuffer(crudo, dtype=_TIPOS_TFORM.get(subtipo, tipo_guardado))
            imagen[cortes] = valores[:n].reshape(forma)
            continue
        if tipo == 'RICE_1':
            pendientes_rice.setdefault(n, []).append((fila, cortes, forma, comprimido))
            continue
        if tipo == 'NOCOMPRESS':
            crudo = comprimido
        else:
            crudo = zlib.decompress(comprimido, 47)
            if tipo == 'GZIP_2':
                crudo = _desordenar_bytes(crudo, np.dtype(tipo_guardado).itemsize)
        enteros = np.frombuffer(crudo, dtype=tipo_guardado)[:n].astype(np.float64 if zbitpix < 0 and not cuantizada
                                                                       else np.int64)
        imagen[cortes] = _valores_tesela(enteros, fila, cuantizada, metodo, semilla, celda_numero, cab).reshape(forma)
    for n, grupo in pendientes_rice.items():
        # se descomprimen por lotes, para no gastar demasiada memoria con imágenes grandes
        lote = max(1, int(30_000_000 // max(1, n * bytes_pixel * 8)))
        for k in range(0, len(grupo), lote):
            parte = grupo[k:k + lote]
            enteros = _rice_decodificar([g[3] for g in parte], n, bytes_pixel, bloque)
            for (fila, cortes, forma, _b), valores in zip(parte, enteros):
                imagen[cortes] = _valores_tesela(valores, fila, cuantizada, metodo, semilla, celda_numero,
                                                 cab).reshape(forma)
    return imagen, zbitpix


def _valores_tesela(enteros, fila, cuantizada, metodo, semilla, celda_numero, cab):
    """Pasa los enteros guardados a valores reales (si la imagen estaba cuantizada)."""
    if not cuantizada:
        return enteros
    escala = celda_numero(fila, 'ZSCALE', 1.0)
    cero = celda_numero(fila, 'ZZERO', 0.0)
    nulo = int(cab.get('ZBLANK', -2147483647))
    valores = enteros.astype(np.float64)
    if metodo == 'NO_DITHER':
        resultado = valores * escala + cero
    else:
        # dithering: se restó un número aleatorio (de la tabla estándar) antes de redondear
        tabla = _tabla_aleatoria()
        iseed = (fila + semilla - 1) % 10000
        aleatorios = np.empty(len(valores))
        hechos, siguiente = 0, int(tabla[iseed] * 500)
        while hechos < len(valores):
            cuantos = min(len(valores) - hechos, 10000 - siguiente)
            aleatorios[hechos:hechos + cuantos] = tabla[siguiente:siguiente + cuantos]
            hechos += cuantos
            iseed = (iseed + 1) % 10000
            siguiente = int(tabla[iseed] * 500)
        resultado = (valores - aleatorios + 0.5) * escala + cero
        if metodo == 'SUBTRACTIVE_DITHER_2':
            resultado[enteros == -2147483646] = 0.0
    resultado[enteros == nulo] = np.nan
    return resultado


def _convertir_datos(datos, cabecera, bitpix, blank=None):
    """Aplica BZERO/BSCALE/BLANK y deja una imagen 2D en float32."""
    bzero = float(cabecera.get('BZERO', 0.0))
    bscale = float(cabecera.get('BSCALE', 1.0))
    blank = cabecera.get('BLANK') if blank is None else blank

    sin_dato = None
    if bitpix > 0 and isinstance(blank, int):
        sin_dato = datos == blank

    # Con datos de 8 o 16 bits (o ya en float32) la precisión simple es exacta: se evita
    # una copia en doble precisión (con una foto de 24 MP son cientos de MB menos).
    tipo = np.float32 if bitpix in (8, 16, -32) else np.float64
    imagen = datos.astype(tipo)
    if bscale != 1.0 or bzero != 0.0:
        imagen *= tipo(bscale)     # en el mismo sitio: sin copias de más
        imagen += tipo(bzero)
    if sin_dato is not None:
        imagen[sin_dato] = np.nan

    nota = ''
    imagen = np.squeeze(imagen)  # quita dimensiones de tamaño 1
    if imagen.ndim == 3 and imagen.shape[0] == 3:
        imagen = np.nanmean(imagen, axis=0)
        nota = 'nota_color'  # clave del texto traducido (ver TEXTOS)
    while imagen.ndim > 2:
        imagen = imagen[0]
        nota = 'nota_planos'
    if imagen.ndim < 2:
        raise ErrorFits('err_no_2d')

    return np.ascontiguousarray(imagen, dtype=np.float32), nota


# =============================================================================
# 2. GIRO Y ESPEJO
# =============================================================================
# El estado de la imagen se guarda como (espejo, ángulo):
#   1) si espejo = True, primero se invierte izquierda <-> derecha;
#   2) después se gira el ángulo indicado en sentido ANTIHORARIO.
# Siempre se calcula a partir de la imagen original, así los giros
# repetidos no van perdiendo calidad.

def normalizar_angulo(angulo):
    """Deja el ángulo entre -180 y 180 grados."""
    angulo = angulo % 360.0
    if angulo > 180.0:
        angulo -= 360.0
    angulo = round(angulo, 9)
    return 0.0 if angulo == 0 else angulo


def es_multiplo_de_90(angulo):
    return abs(angulo / 90.0 - round(angulo / 90.0)) < 1e-9


def _seno_coseno(angulo):
    """Seno y coseno del ángulo; exactos (0, 1, -1) en múltiplos de 90°."""
    if es_multiplo_de_90(angulo):
        k = int(round(angulo / 90.0)) % 4
        return (0, 1, 0, -1)[k], (1, 0, -1, 0)[k]
    a = math.radians(angulo)
    return math.sin(a), math.cos(a)


def tamano_girado(ancho, alto, angulo):
    """Tamaño (ancho, alto) de la imagen girada, para que quepa entera."""
    s, c = _seno_coseno(angulo)
    nuevo_ancho = int(math.ceil(abs(ancho * c) + abs(alto * s) - 1e-6))
    nuevo_alto = int(math.ceil(abs(ancho * s) + abs(alto * c) - 1e-6))
    return nuevo_ancho, nuevo_alto


def girar_imagen(imagen, angulo):
    """Gira la imagen `angulo` grados en sentido antihorario.

    - Si el ángulo es múltiplo de 90°, el giro es exacto (solo se reordenan
      los píxeles, sin modificar ningún valor).
    - Si no, se usa interpolación bilineal: para cada píxel de la imagen
      nueva se calcula de qué punto de la original viene y se mezclan los
      4 píxeles vecinos según la distancia. Las esquinas que quedan vacías
      se marcan como NaN (sin dato).
    """
    if es_multiplo_de_90(angulo):
        k = int(round(angulo / 90.0)) % 4
        # Como la fila 0 se dibuja abajo, rot90 con -k gira en sentido antihorario.
        return np.ascontiguousarray(np.rot90(imagen, -k))

    alto, ancho = imagen.shape
    if alto < 2 or ancho < 2:
        raise ErrorFits('err_pequena')
    s, c = _seno_coseno(angulo)
    nuevo_ancho, nuevo_alto = tamano_girado(ancho, alto, angulo)
    cx, cy = (ancho - 1) / 2.0, (alto - 1) / 2.0                # centro original
    ncx, ncy = (nuevo_ancho - 1) / 2.0, (nuevo_alto - 1) / 2.0  # centro nuevo

    salida = np.empty((nuevo_alto, nuevo_ancho), dtype=np.float32)
    plana = np.ascontiguousarray(imagen, dtype=np.float32).ravel()   # índices planos: más rápido
    dx = (np.arange(nuevo_ancho) - ncx)[None, :]
    filas_por_trozo = max(1, 1_000_000 // nuevo_ancho)  # por trozos, para no gastar mucha memoria

    for fila_ini in range(0, nuevo_alto, filas_por_trozo):
        fila_fin = min(fila_ini + filas_por_trozo, nuevo_alto)
        dy = (np.arange(fila_ini, fila_fin) - ncy)[:, None]
        # Giro inverso: punto de la imagen original del que viene cada píxel nuevo.
        x_orig = cx + c * dx + s * dy
        y_orig = cy - s * dx + c * dy
        dentro = (x_orig >= 0) & (x_orig <= ancho - 1) & (y_orig >= 0) & (y_orig <= alto - 1)

        x0 = np.clip(np.floor(x_orig).astype(np.int64), 0, ancho - 2)
        y0 = np.clip(np.floor(y_orig).astype(np.int64), 0, alto - 2)
        fx = (x_orig - x0).astype(np.float32)
        fy = (y_orig - y0).astype(np.float32)
        i00 = y0 * ancho + x0
        a, b = plana[i00], plana[i00 + 1]
        abajo = a + (b - a) * fx
        a, b = plana[i00 + ancho], plana[i00 + ancho + 1]
        arriba = a + (b - a) * fx
        valor = abajo + (arriba - abajo) * fy
        valor[~dentro] = np.nan
        salida[fila_ini:fila_fin] = valor
    return salida


def transformar_imagen(original, angulo, espejo):
    """Aplica el espejo (si toca) y después el giro."""
    imagen = original[:, ::-1] if espejo else original
    return girar_imagen(imagen, angulo)


def transformar_punto(x, y, ancho, alto, angulo, espejo):
    """Dónde acaba el punto (x, y) de la imagen original después de aplicar
    espejo y giro. Es la misma cuenta que hace girar_imagen pero al revés."""
    if espejo:
        x = (ancho - 1) - x
    s, c = _seno_coseno(angulo)
    nuevo_ancho, nuevo_alto = tamano_girado(ancho, alto, angulo)
    dx = x - (ancho - 1) / 2.0
    dy = y - (alto - 1) / 2.0
    return ((nuevo_ancho - 1) / 2.0 + c * dx - s * dy,
            (nuevo_alto - 1) / 2.0 + s * dx + c * dy)


def destransformar_punto(x, y, ancho, alto, angulo, espejo):
    """Lo contrario de transformar_punto: de un punto de la foto girada/invertida
    vuelve a sus coordenadas en la foto original."""
    s, c = _seno_coseno(angulo)
    nuevo_ancho, nuevo_alto = tamano_girado(ancho, alto, angulo)
    dx = x - (nuevo_ancho - 1) / 2.0
    dy = y - (nuevo_alto - 1) / 2.0
    x_orig = (ancho - 1) / 2.0 + c * dx + s * dy   # giro inverso
    y_orig = (alto - 1) / 2.0 - s * dx + c * dy
    if espejo:
        x_orig = (ancho - 1) - x_orig
    return x_orig, y_orig


# =============================================================================
# 3. CÁLCULO DEL CENTROIDE (PRECISIÓN SUBPÍXEL)
# =============================================================================
# Pasos:
#   1) Buscar el píxel más brillante cerca del clic (suavizando 3x3 para que
#      un píxel "caliente" suelto no engañe).
#   2) Medir el FONDO del cielo en un anillo alrededor de la estrella, con
#      recorte sigma (se descartan valores raros, como otras estrellas).
#   3) Comprobar que la estrella destaca sobre el ruido (umbral en sigmas).
#   4) Centroide iterativo ("centro de masas" de la luz):
#          Xc = suma(I * x) / suma(I)      Yc = suma(I * y) / suma(I)
#      donde I = brillo del píxel - fondo, usando solo los píxeles dentro de
#      un círculo (apertura). El círculo se vuelve a centrar en el nuevo
#      centroide y se repite hasta que se mueve menos de 0,0001 píxeles.
#   5) Extras: flujo total, FWHM aproximada (anchura de la estrella).

def _centro_rapido(imagen, x, y, fondo, radio=2.5):
    """Centro aproximado de una estrella vecina: centro de masas de su luz en un
    círculo pequeño, recentrado unas cuantas veces."""
    alto, ancho = imagen.shape
    m = int(math.ceil(radio)) + 2
    x0, x1, y0, y1 = max(x - m, 0), min(x + m + 1, ancho), max(y - m, 0), min(y + m + 1, alto)
    trozo = np.clip(np.nan_to_num(imagen[y0:y1, x0:x1].astype(np.float64) - fondo), 0, None)
    yy, xx = np.mgrid[y0:y1, x0:x1]
    cx, cy = float(x), float(y)
    for _ in range(10):
        pesos = np.where((xx - cx) ** 2 + (yy - cy) ** 2 <= radio ** 2, trozo, 0.0)
        total = pesos.sum()
        if total <= 0:
            break
        nx, ny = float((pesos * xx).sum() / total), float((pesos * yy).sum() / total)
        if math.hypot(nx - cx, ny - cy) < 1e-3:
            cx, cy = nx, ny
            break
        cx, cy = nx, ny
    return cx, cy


def _luz_reflejada(imagen, xs, ys, fondo):
    """Brillo (menos el fondo) de la imagen en los puntos (xs, ys), con interpolación
    bilineal. Fuera de la imagen vale 0."""
    alto, ancho = imagen.shape
    dentro = (xs >= 0) & (xs <= ancho - 1) & (ys >= 0) & (ys <= alto - 1)
    x0 = np.clip(np.floor(xs).astype(int), 0, ancho - 2)
    y0 = np.clip(np.floor(ys).astype(int), 0, alto - 2)
    fx, fy = xs - x0, ys - y0
    valor = (imagen[y0, x0] * (1 - fx) * (1 - fy) + imagen[y0, x0 + 1] * fx * (1 - fy)
             + imagen[y0 + 1, x0] * (1 - fx) * fy + imagen[y0 + 1, x0 + 1] * fx * fy)
    return np.where(dentro, np.clip(np.nan_to_num(valor - fondo), 0, None), 0.0)


def recorte_sigma(valores, n_sigma=3.0, iteraciones=5):
    """Mediana y desviación típica del fondo, descartando valores atípicos."""
    valores = valores[np.isfinite(valores)]
    for _ in range(iteraciones):
        mediana = np.median(valores)
        sigma = np.std(valores)
        if sigma == 0:
            break
        quedan = valores[np.abs(valores - mediana) < n_sigma * sigma]
        if quedan.size == valores.size or quedan.size < 5:
            break
        valores = quedan
    return float(np.median(valores)), float(np.std(valores))


def medir_estrella(imagen, x_clic, y_clic, radio_busqueda=15, radio_apertura=5.0,
                   umbral_sigma=3.0, ganancia=None):
    """Mide la estrella más cercana al clic (x_clic, y_clic).

    `ganancia` (electrones por unidad del FITS, palabra GAIN de la cabecera) sirve
    para incluir el ruido de la propia estrella en las incertidumbres; si no se
    conoce, solo se tiene en cuenta el ruido del fondo del cielo.
    Devuelve un diccionario con el resultado, o None si no hay estrella.
    """
    alto, ancho = imagen.shape
    radio_busqueda = int(radio_busqueda)

    # --- Paso 1: la estrella más cercana al clic ------------------------------
    # Se buscan los "picos" (píxeles más brillantes que sus 8 vecinos) que
    # destacan sobre el ruido y se elige el más cercano al clic. Así, si hay
    # una estrella brillante al lado, no "roba" el clic a la estrella débil.
    cx, cy = int(round(x_clic)), int(round(y_clic))
    x0, x1 = max(cx - radio_busqueda, 0), min(cx + radio_busqueda + 1, ancho)
    y0, y1 = max(cy - radio_busqueda, 0), min(cy + radio_busqueda + 1, alto)
    if x1 - x0 < 3 or y1 - y0 < 3:
        return None
    ventana = imagen[y0:y1, x0:x1].astype(np.float64)
    finitos = np.isfinite(ventana)
    if not finitos.any():
        return None
    ventana = np.where(finitos, ventana, np.nanmedian(ventana))

    h, w = ventana.shape
    ampliada = np.pad(ventana, 1, mode='edge')
    suavizada = sum(ampliada[j:j + h, i:i + w] for j in range(3) for i in range(3)) / 9.0
    yy, xx = np.mgrid[y0:y1, x0:x1]
    distancia2 = (xx - x_clic) ** 2 + (yy - y_clic) ** 2
    fuera = distancia2 > radio_busqueda ** 2
    # Los picos se buscan en la imagen sin suavizar: si se suaviza, la luz de una
    # estrella brillante "tapa" a una débil que esté a pocos píxeles.
    mediana = np.median(ventana[~fuera])
    sigma = 1.4826 * np.median(np.abs(ventana[~fuera] - mediana))  # ruido (robusto a las estrellas)
    if sigma <= 0:
        sigma = np.std(ventana[~fuera]) or 1e-12
    suavizada[fuera] = -np.inf
    cruda = np.where(fuera, -np.inf, ventana)
    rodeada = np.pad(cruda, 1, mode='constant', constant_values=-np.inf)
    es_pico = np.ones_like(cruda, dtype=bool)
    for j in range(3):
        for i in range(3):
            if (j, i) != (1, 1):
                es_pico &= cruda >= rodeada[j:j + h, i:i + w]
    picos = es_pico & ~fuera & (cruda - mediana > max(umbral_sigma, 3.0) * sigma)
    otros_picos = []
    if picos.any():
        filas, columnas = np.nonzero(picos)
        k = int(np.argmin(distancia2[filas, columnas]))
        fila, columna = filas[k], columnas[k]
        # Los "picos" pegados al elegido, o del mismo valor a pocos píxeles, son la MISMA
        # estrella: pasa con las estrellas saturadas, que tienen el centro plano (una
        # "meseta" de píxeles iguales). No son vecinas y no hay que restarles la luz.
        valor = cruda[fila, columna]
        otros_picos = [(x0 + c, y0 + f) for f, c in zip(filas, columnas)
                       if (f, c) != (fila, columna) and math.hypot(f - fila, c - columna) > 2.0
                       and not (abs(cruda[f, c] - valor) <= 1e-6 * abs(valor) + 1e-12
                                and math.hypot(f - fila, c - columna) <= 8.0)]
    else:  # no hay ningún pico claro: el píxel más brillante (luego se comprueba si destaca)
        fila, columna = np.unravel_index(np.argmax(suavizada), suavizada.shape)
    px, py = x0 + columna, y0 + fila

    # Si hay otra estrella muy cerca, el círculo para el centro se hace más pequeño
    # (como mucho, algo más de la mitad de la distancia entre las dos) para no tocarla.
    if otros_picos:
        mas_cercana = min(math.hypot(ox - px, oy - py) for ox, oy in otros_picos)
        radio_apertura = max(1.5, min(radio_apertura, 0.55 * mas_cercana))

    # --- Paso 2: fondo del cielo en un anillo --------------------------------
    r_int = radio_apertura + 3
    r_ext = r_int + 6
    m = int(math.ceil(r_ext)) + 1
    sx0, sx1 = max(px - m, 0), min(px + m + 1, ancho)
    sy0, sy1 = max(py - m, 0), min(py + m + 1, alto)
    recorte = imagen[sy0:sy1, sx0:sx1].astype(np.float64)
    gy, gx = np.mgrid[sy0:sy1, sx0:sx1]
    dist2 = (gx - px) ** 2 + (gy - py) ** 2
    anillo = recorte[(dist2 >= r_int ** 2) & (dist2 <= r_ext ** 2)]
    anillo = anillo[np.isfinite(anillo)]
    if anillo.size < 8:
        return None
    fondo, sigma_fondo = recorte_sigma(anillo)

    # --- Paso 3: ¿la estrella destaca sobre el ruido? ------------------------
    pico = float(imagen[py, px]) - fondo
    if not np.isfinite(pico) or pico <= 0 or pico < umbral_sigma * sigma_fondo:
        return None

    # --- Paso 4: centroide iterativo -----------------------------------------
    luz = np.nan_to_num(recorte - fondo, nan=0.0)
    # Si hay otra estrella muy cerca, se le quita su luz para que no arrastre el centro.
    # Truco: una estrella es simétrica, así que su luz en un punto es la misma que en el
    # punto opuesto respecto a su centro (que queda lejos de nuestra estrella).
    for ox, oy in otros_picos:
        if math.hypot(ox - px, oy - py) > radio_apertura + 8:
            continue
        if imagen[oy, ox] - fondo < 0.5 * (imagen[py, px] - fondo):
            continue  # solo molestan las vecinas brillantes; las débiles casi no influyen
        cxv, cyv = _centro_rapido(imagen, ox, oy, fondo)
        luz -= _luz_reflejada(imagen, 2 * cxv - gx, 2 * cyv - gy, fondo)
    xc, yc = float(px), float(py)
    # Solo cuentan para el centro los píxeles que superan el fondo en más de
    # `umbral_sigma` veces el ruido (la «Sensibilidad»). Si contaran todos, el ruido
    # del cielo que hay en el círculo arrastraría el centro de las estrellas débiles.
    corte = max(umbral_sigma, 0.0) * sigma_fondo
    luz_util = np.clip(luz - corte, 0.0, None)
    iteraciones = 0
    for iteraciones in range(1, 101):
        dentro = (gx - xc) ** 2 + (gy - yc) ** 2 <= radio_apertura ** 2
        pesos = np.where(dentro, luz_util, 0.0)
        total = pesos.sum()
        if total <= 0:
            return None
        nuevo_x = float((pesos * gx).sum() / total)
        nuevo_y = float((pesos * gy).sum() / total)
        movimiento = math.hypot(nuevo_x - xc, nuevo_y - yc)
        xc, yc = nuevo_x, nuevo_y
        if movimiento < 1e-4:
            break

    # --- Paso 5: datos extra --------------------------------------------------
    dentro = (gx - xc) ** 2 + (gy - yc) ** 2 <= radio_apertura ** 2
    pesos = np.where(dentro, luz_util, 0.0)
    total = pesos.sum()
    momento2 = (pesos * ((gx - xc) ** 2 + (gy - yc) ** 2)).sum() / total
    fwhm = 2.3548 * math.sqrt(momento2 / 2.0)  # FWHM = 2,3548 * sigma (gaussiana)
    flujo = float(luz[dentro].sum())

    # --- Paso 6: incertidumbres (propagación de errores) ----------------------
    # Cada píxel tiene un ruido σ_i² = ruido del cielo² (+ señal/ganancia, si se conoce).
    # Como Xc = Σ w_i·x_i / Σ w_i, su error es σ_X² = Σ (x_i − Xc)²·σ_i² / (Σ w_i)².
    var_pixel = np.full_like(luz, sigma_fondo ** 2)
    if ganancia and ganancia > 0:
        var_pixel = var_pixel + np.clip(luz, 0, None) / ganancia
    cuentan = pesos > 0
    sigma_x = math.sqrt(float((((gx - xc) ** 2) * var_pixel)[cuentan].sum())) / total
    sigma_y = math.sqrt(float((((gy - yc) ** 2) * var_pixel)[cuentan].sum())) / total
    n_ap = int(dentro.sum())
    var_flujo = n_ap * sigma_fondo ** 2 * (1 + n_ap / max(anillo.size, 1))  # fondo (y su error)
    if ganancia and ganancia > 0:
        var_flujo += max(flujo, 0.0) / ganancia                              # ruido de la estrella
    sigma_flujo = math.sqrt(var_flujo)

    toca_borde = (xc - radio_apertura < 0 or yc - radio_apertura < 0
                  or xc + radio_apertura > ancho - 1 or yc + radio_apertura > alto - 1
                  or not np.isfinite(recorte[dentro]).all())

    return {
        'x': xc,
        'y': yc,
        'fwhm': fwhm,
        'flujo': flujo,
        'pico': pico,
        'fondo': fondo,
        'sigma_fondo': sigma_fondo,
        'pico_sigma': pico / sigma_fondo if sigma_fondo > 0 else float('inf'),
        'iteraciones': iteraciones,
        'borde': toca_borde,
        'sigma_x': sigma_x,
        'sigma_y': sigma_y,
        'sigma_flujo': sigma_flujo,
    }


def magnitud_instrumental(flujo, sigma_flujo, exposicion=None):
    """m = 25 − 2,5·log10(flujo / tiempo de exposición), y su error 1,0857·σF/F.
    (El 25 es solo un "punto cero" para que salgan números cómodos; para
    magnitudes reales se compara con una estrella de magnitud conocida.)"""
    if flujo is None or flujo <= 0:
        return None, None
    por_segundo = flujo / exposicion if exposicion and exposicion > 0 else flujo
    return 25.0 - 2.5 * math.log10(por_segundo), 1.0857 * sigma_flujo / flujo


# =============================================================================
# 3b. COORDENADAS CELESTES (WCS) Y TIEMPO
# =============================================================================

class SolucionWCS:
    """Convierte píxeles en coordenadas celestes (AR, Dec) con la "solución
    astrométrica" de la cabecera FITS (la que añaden programas como Astrometry.net).
    Solo la proyección gnomónica (TAN), con o sin corrección de distorsión SIP,
    que es la habitual en fotos de telescopio."""

    def __init__(self, cabecera):
        c = cabecera
        tipo1, tipo2 = str(c.get('CTYPE1', '')), str(c.get('CTYPE2', ''))
        if not (tipo1.startswith('RA---TAN') and tipo2.startswith('DEC--TAN')):
            raise ValueError('sin TAN')
        self.crpix = (float(c['CRPIX1']), float(c['CRPIX2']))
        self.crval = (math.radians(float(c['CRVAL1'])), math.radians(float(c['CRVAL2'])))
        if 'CD1_1' in c:
            self.cd = ((float(c.get('CD1_1', 0)), float(c.get('CD1_2', 0))),
                       (float(c.get('CD2_1', 0)), float(c.get('CD2_2', 0))))
        else:
            d1, d2 = float(c.get('CDELT1', 1)), float(c.get('CDELT2', 1))
            if 'PC1_1' in c:
                pc = ((float(c.get('PC1_1', 1)), float(c.get('PC1_2', 0))),
                      (float(c.get('PC2_1', 0)), float(c.get('PC2_2', 1))))
            else:
                rot = math.radians(float(c.get('CROTA2', 0)))
                pc = ((math.cos(rot), -math.sin(rot) * d2 / d1), (math.sin(rot) * d1 / d2, math.cos(rot)))
            self.cd = ((d1 * pc[0][0], d1 * pc[0][1]), (d2 * pc[1][0], d2 * pc[1][1]))
        if abs(self.cd[0][0] * self.cd[1][1] - self.cd[0][1] * self.cd[1][0]) == 0:
            raise ValueError('matriz CD singular')
        # Distorsión SIP (polinomios A y B), si la hay
        self.sip = None
        if tipo1.endswith('-SIP'):
            a = {(p, q): float(c[f'A_{p}_{q}']) for p in range(10) for q in range(10) if f'A_{p}_{q}' in c}
            b = {(p, q): float(c[f'B_{p}_{q}']) for p in range(10) for q in range(10) if f'B_{p}_{q}' in c}
            self.sip = (a, b)
        self.escala = 3600 * math.sqrt(abs(self.cd[0][0] * self.cd[1][1] - self.cd[0][1] * self.cd[1][0]))

    @classmethod
    def desde_cabecera(cls, cabecera):
        try:
            return cls(cabecera)
        except (KeyError, ValueError, TypeError):
            return None

    def a_cielo(self, x_fits, y_fits):
        """Píxel FITS (empieza en 1) -> (AR, Dec) en grados."""
        u, v = x_fits - self.crpix[0], y_fits - self.crpix[1]
        if self.sip:
            a, b = self.sip
            du = sum(coef * u ** p * v ** q for (p, q), coef in a.items())
            dv = sum(coef * u ** p * v ** q for (p, q), coef in b.items())
            u, v = u + du, v + dv
        xi = math.radians(self.cd[0][0] * u + self.cd[0][1] * v)    # coordenadas "estándar" en el plano
        eta = math.radians(self.cd[1][0] * u + self.cd[1][1] * v)
        a0, d0 = self.crval
        denominador = math.cos(d0) - eta * math.sin(d0)               # proyección gnomónica inversa
        ar = a0 + math.atan2(xi, denominador)
        dec = math.atan2(math.sin(d0) + eta * math.cos(d0), math.hypot(xi, denominador))
        return math.degrees(ar) % 360.0, math.degrees(dec)

    def a_pixel(self, ar, dec):
        """(AR, Dec) en grados -> píxel FITS (empieza en 1). Es la inversa de a_cielo:
        primero sin distorsión (invirtiendo la matriz CD) y luego unas iteraciones de
        Newton para tener en cuenta la distorsión SIP, si la hay."""
        xi, eta = C.proyectar(ar, dec, math.degrees(self.crval[0]), math.degrees(self.crval[1]))
        (a, b), (c, d) = self.cd
        det = a * d - b * c
        u = (d * float(xi) - b * float(eta)) / det
        v = (-c * float(xi) + a * float(eta)) / det
        x, y = u + self.crpix[0], v + self.crpix[1]
        if self.sip:
            for _ in range(8):
                ar1, dec1 = self.a_cielo(x, y)
                xi1, eta1 = C.proyectar(ar1, dec1, math.degrees(self.crval[0]), math.degrees(self.crval[1]))
                dxi, deta = float(xi) - float(xi1), float(eta) - float(eta1)
                x += (d * dxi - b * deta) / det
                y += (-c * dxi + a * deta) / det
                if abs(dxi) + abs(deta) < 1e-10:
                    break
        return x, y

    def como_cabecera(self):
        """La solución como palabras clave FITS (para guardarla en una sesión)."""
        cab = {'CTYPE1': 'RA---TAN' + ('-SIP' if self.sip else ''), 'CTYPE2': 'DEC--TAN' + ('-SIP' if self.sip else ''),
               'CRPIX1': self.crpix[0], 'CRPIX2': self.crpix[1],
               'CRVAL1': math.degrees(self.crval[0]), 'CRVAL2': math.degrees(self.crval[1]),
               'CD1_1': self.cd[0][0], 'CD1_2': self.cd[0][1], 'CD2_1': self.cd[1][0], 'CD2_2': self.cd[1][1]}
        if self.sip:
            for letra, coef in zip('AB', self.sip):
                for (p, q), valor in coef.items():
                    cab[f'{letra}_{p}_{q}'] = valor
        return cab


def ar_hms(grados):
    """AR en horas, minutos y segundos (texto)."""
    segundos = round((grados % 360) / 15 * 3600, 2)
    h, resto = divmod(segundos, 3600)
    m, s = divmod(resto, 60)
    return f'{int(h) % 24:02d}h {int(m):02d}m {s:05.2f}s'.replace('.', ',')


def dec_dms(grados):
    """Declinación en grados, minutos y segundos (texto)."""
    signo = '-' if grados < 0 else '+'
    segundos = round(abs(grados) * 3600, 1)
    g, resto = divmod(segundos, 3600)
    m, s = divmod(resto, 60)
    return f'{signo}{int(g):02d}° {int(m):02d}′ {s:04.1f}″'.replace('.', ',')


def dia_juliano(anyo, mes, dia, horas=0.0):
    """Día juliano (JD) de una fecha del calendario gregoriano en tiempo universal
    (algoritmo de Meeus, «Astronomical Algorithms»)."""
    if mes <= 2:
        anyo, mes = anyo - 1, mes + 12
    a = anyo // 100
    b = 2 - a + a // 4
    return int(365.25 * (anyo + 4716)) + int(30.6001 * (mes + 1)) + dia + b - 1524.5 + horas / 24.0


def instante_observacion(cabecera):
    """Momento central de la exposición a partir de DATE-OBS (+ TIME-OBS/UT) y
    EXPTIME. Devuelve (texto UTC, JD) o (texto original, None) si no se entiende."""
    fecha = re.sub(r'\s*(GMT|UTC|UT)\s*$', '', str(cabecera.get('DATE-OBS', '') or '').strip())
    if not fecha:
        return '', None
    # Formato antiguo de la norma FITS (antes de 1999): 'DD/MM/AA', siempre del siglo XX
    antiguo = re.match(r'^(\d{1,2})/(\d{1,2})/(\d{2})$', fecha)
    if antiguo:
        fecha = f'19{antiguo.group(3)}-{int(antiguo.group(2)):02d}-{int(antiguo.group(1)):02d}'
    if 'T' not in fecha:
        hora = cabecera.get('TIME-OBS', '') or cabecera.get('UT', '') or cabecera.get('UTSTART', '') or ''
        if isinstance(hora, (int, float)) and not isinstance(hora, bool) and 0 <= hora < 24:
            hora = f'{int(hora):02d}:{int(hora * 60 % 60):02d}:{hora * 3600 % 60:06.3f}'   # horas con decimales
        hora = re.sub(r'\s*(GMT|UTC|UT)\s*$', '', str(hora).strip())
        if hora:
            fecha = f'{fecha}T{hora}'
    try:
        parte_fecha, _, parte_hora = fecha.partition('T')
        anyo, mes, dia = (int(p) for p in parte_fecha.split('-')[:3])
        h, m, s = (parte_hora.split(':') + ['0', '0', '0'])[:3] if parte_hora else ('0', '0', '0')
        segundos = int(h) * 3600 + int(m) * 60 + float(s.rstrip('Z'))
        exposicion = float(cabecera.get('EXPTIME', cabecera.get('EXPOSURE', 0)) or 0)
        segundos += max(exposicion, 0) / 2                               # centro de la exposición
        jd = dia_juliano(anyo, mes, dia, segundos / 3600.0)
        dias_extra, resto = divmod(segundos, 86400)
        hh, resto = divmod(resto, 3600)
        mm, ss = divmod(resto, 60)
        texto = f'{anyo:04d}-{mes:02d}-{dia:02d} {int(hh):02d}:{int(mm):02d}:{ss:06.3f}'
        if dias_extra:  # la exposición acaba pasada la medianoche: se toma la fecha del JD
            texto = fecha_desde_jd(jd)
        return texto, jd
    except (ValueError, TypeError):
        return fecha, None


def fecha_desde_jd(jd):
    """Texto UTC 'AAAA-MM-DD hh:mm:ss.sss' a partir de un día juliano (Meeus)."""
    z, f = divmod(jd + 0.5, 1)
    z = int(z)
    alfa = int((z - 1867216.25) / 36524.25)
    a = z + 1 + alfa - alfa // 4
    b = a + 1524
    c = int((b - 122.1) / 365.25)
    d = int(365.25 * c)
    e = int((b - d) / 30.6001)
    dia = b - d - int(30.6001 * e)
    mes = e - 1 if e < 14 else e - 13
    anyo = c - 4716 if mes > 2 else c - 4715
    segundos = f * 86400
    hh, resto = divmod(segundos, 3600)
    mm, ss = divmod(resto, 60)
    return f'{anyo:04d}-{mes:02d}-{dia:02d} {int(hh):02d}:{int(mm):02d}:{ss:06.3f}'


# -----------------------------------------------------------------------------
# Informe para el Minor Planet Center (formato de 80 columnas, "obs80")
# -----------------------------------------------------------------------------
# Cada observación es una línea de exactamente 80 caracteres:
#   1-5 número (vacío si no tiene)   6-12 designación provisional   13 descubrimiento (*)
#   15 tipo de observación (C = CCD, B = CMOS)   16-32 fecha UTC "AAAA MM DD.dddddd"
#   33-44 AR "HH MM SS.ss"   45-56 Dec "sDD MM SS.s"   66-70 magnitud   71 banda
#   78-80 código del observatorio

def _calendario_desde_jd(jd):
    """(año, mes, día con decimales) a partir del día juliano (Meeus)."""
    z, f = divmod(jd + 0.5, 1)
    z = int(z)
    alfa = int((z - 1867216.25) / 36524.25)
    a = z + 1 + alfa - alfa // 4
    b = a + 1524
    c = int((b - 122.1) / 365.25)
    d = int(365.25 * c)
    e = int((b - d) / 30.6001)
    dia = b - d - int(30.6001 * e) + f
    mes = e - 1 if e < 14 else e - 13
    anyo = c - 4716 if mes > 2 else c - 4715
    return anyo, mes, dia


def fecha_mpc(jd):
    """Fecha con 5 decimales de día (0,9 s): la precisión habitual para un reloj de ordenador
    sincronizado (el MPC pide no dar más precisión de la que se tiene)."""
    anyo, mes, dia = _calendario_desde_jd(round(jd * 1e5) / 1e5)
    return f'{anyo:04d} {mes:02d} {dia:08.5f}'


def ar_mpc(grados):
    centesimas = round((grados % 360) / 15 * 3600 * 100) % (24 * 3600 * 100)
    h, resto = divmod(centesimas, 360000)
    m, resto = divmod(resto, 6000)
    return f'{h:02d} {m:02d} {resto / 100:05.2f}'


def dec_mpc(grados):
    signo = '-' if grados < 0 else '+'
    decimas = round(abs(grados) * 36000)
    g, resto = divmod(decimas, 36000)
    m, resto = divmod(resto, 600)
    return f'{signo}{g:02d} {m:02d} {resto / 10:04.1f}'


def linea_mpc(designacion, descubrimiento, tipo, jd, ar, dec, magnitud, banda, codigo):
    """Una línea de observación en el formato de 80 columnas del MPC."""
    mag = f'{magnitud:4.1f} ' if magnitud is not None and 0 < magnitud < 100 else '     '
    linea = (' ' * 5 + f'{designacion.strip()[:7]:<7}' + ('*' if descubrimiento else ' ') + ' ' + tipo[:1]
             + f'{fecha_mpc(jd):<17}' + f'{ar_mpc(ar):<12}' + f'{dec_mpc(dec):<12}' + ' ' * 9
             + mag + (banda[:1] if mag.strip() else ' ') + ' ' * 6 + f'{codigo.strip()[:3]:>3}')
    assert len(linea) == 80
    return linea


def ganancia_de(cabecera):
    """Ganancia (electrones por unidad) de la cabecera FITS, si viene."""
    for clave in ('GAIN', 'EGAIN', 'CCDGAIN'):
        try:
            valor = float(cabecera.get(clave))
            if 0 < valor < 1000:
                return valor
        except (TypeError, ValueError):
            pass
    return None


def exposicion_de(cabecera):
    try:
        valor = float(cabecera.get('EXPTIME', cabecera.get('EXPOSURE')))
        return valor if valor > 0 else None
    except (TypeError, ValueError):
        return None


def _tarjeta_fits(clave, valor, comentario=''):
    """Una línea (tarjeta) de 80 caracteres de una cabecera FITS, o None si no se puede escribir."""
    if isinstance(valor, bool):
        texto = ('T' if valor else 'F').rjust(20)
    elif isinstance(valor, (int, np.integer)):
        texto = str(int(valor)).rjust(20)
    elif isinstance(valor, (float, np.floating)):
        if not math.isfinite(valor):
            return None
        texto = f'{float(valor):.16G}'
        if '.' not in texto and 'E' not in texto:
            texto += '.0'
        texto = texto.rjust(20)
    elif isinstance(valor, str):
        limpio = ''.join(c if 32 <= ord(c) < 127 else '?' for c in valor)[:66].replace("'", "''")[:66]
        texto = f"'{limpio.ljust(8)}'"
    else:
        return None
    linea = f'{clave.ljust(8)}= {texto}'
    if comentario:
        linea += ' / ' + comentario
    return linea[:80].ljust(80)


def escribir_fits(ruta, imagen, cabecera, comentarios=()):
    """Guarda una imagen 2D (fila 0 = la de abajo, como al leerla) como FITS de 32 bits en
    coma flotante, con las palabras clave de 'cabecera' (las de estructura se ponen nuevas)."""
    alto, ancho = imagen.shape
    tarjetas = [_tarjeta_fits('SIMPLE', True), _tarjeta_fits('BITPIX', -32), _tarjeta_fits('NAXIS', 2),
                _tarjeta_fits('NAXIS1', ancho), _tarjeta_fits('NAXIS2', alto)]
    for clave, valor in cabecera.items():
        if (not isinstance(clave, str) or clave.startswith('_') or _ESTRUCTURALES.match(clave)
                or not re.match(r'^[A-Z0-9_-]{1,8}$', clave) or clave in ('COMMENT', 'HISTORY')):
            continue
        tarjeta = _tarjeta_fits(clave, valor)
        if tarjeta:
            tarjetas.append(tarjeta)
    for texto in comentarios:
        limpio = ''.join(c if 32 <= ord(c) < 127 else '?' for c in texto)
        for i in range(0, max(len(limpio), 1), 72):
            tarjetas.append(('HISTORY ' + limpio[i:i + 72]).ljust(80))
    tarjetas.append('END'.ljust(80))
    texto = ''.join(tarjetas)
    texto += ' ' * ((-len(texto)) % TAMANO_BLOQUE)
    datos = np.ascontiguousarray(imagen, dtype='>f4').tobytes()
    temporal = ruta + '.tmp'
    with open(temporal, 'wb') as archivo:
        archivo.write(texto.encode('ascii'))
        archivo.write(datos)
        archivo.write(b'\0' * ((-len(datos)) % TAMANO_BLOQUE))
    os.replace(temporal, ruta)


# =============================================================================
# 4. VENTANA DEL PROGRAMA
# =============================================================================

COLOR_FUERA = 10  # casi negro para la zona fuera de la imagen

# Paleta "estelar": azul noche, luz de estrella y dorado
FONDO = '#0b1022'          # fondo de la ventana (cielo nocturno)
ESPACIO = '#04060e'        # fondo del visor de la foto
CAMPO = '#141c36'          # cajas de texto, tabla
BORDE = '#26335a'
BOTON = '#1a2444'
BOTON_ENCIMA = '#25335e'
BOTON_PULSADO = '#31437a'
PRINCIPAL = '#2c4f9e'      # botón principal
PRINCIPAL_ENCIMA = '#3a66c4'
TEXTO = '#e4eaff'
TEXTO_SUAVE = '#8a96bb'
COLOR_ACENTO = '#86c8ff'   # azul estrella
DORADO = '#ffd27a'
ROJO_SUAVE = '#ff8a9a'
FUENTE_TITULO = 'Segoe UI Light'


# -----------------------------------------------------------------------------
# Decoración animada: cielo que se desplaza, título flotante, estrellas fugaces
# -----------------------------------------------------------------------------
# Todo se dibuja con figuras simples del Canvas de tkinter (óvalos, líneas y
# polígonos) y se mueve unas 30 veces por segundo con root.after().

def mezclar(color1, color2, t):
    """Color intermedio entre color1 (t = 0) y color2 (t = 1)."""
    t = min(max(t, 0.0), 1.0)
    a = [int(color1[i:i + 2], 16) for i in (1, 3, 5)]
    b = [int(color2[i:i + 2], 16) for i in (1, 3, 5)]
    return '#%02x%02x%02x' % tuple(int(round(p + (q - p) * t)) for p, q in zip(a, b))


def puntos_destello(x, y, r, angulo=0.0):
    """Vértices de una estrella de 4 puntas (destello), girada `angulo` radianes."""
    puntos = []
    for k in range(8):
        a = angulo - math.pi / 2 + k * math.pi / 4
        radio = r if k % 2 == 0 else r * 0.3
        puntos += [x + radio * math.cos(a), y + radio * math.sin(a)]
    return puntos


def dibujar_destello(lienzo, x, y, r, color, angulo=0.0, etiquetas=()):
    return lienzo.create_polygon(*puntos_destello(x, y, r, angulo), fill=color, outline='',
                                 tags=etiquetas)


class CampoEstelar:
    """Estrellas que se desplazan despacio hacia la izquierda, como si el
    telescopio hiciera un paneo del cielo. Las estrellas "cercanas" (más
    grandes y brillantes) se mueven más rápido que las lejanas (paralaje),
    y de vez en cuando alguna titila."""

    def __init__(self, lienzo, ancho, alto, semilla, densidad, velocidad, x_inicio=0, destellos=False):
        self.lienzo = lienzo
        self.ancho, self.alto = ancho, alto
        self.x_inicio = x_inicio
        self.velocidad = velocidad  # píxeles por segundo de las estrellas más cercanas
        self.azar = random.Random(semilla)
        # Cada estrella: [objeto, x, y, profundidad, color_base, desvío_x, desvío_y]
        # (x, y) es su posición "real"; el desvío es lo que la mueve una lente gravitatoria.
        self.estrellas = []
        self.lente = None  # (centro_x, centro_y, radio_de_Einstein) mientras dura GARGANTUA
        self.titileo = True  # se apaga mientras M16 vacía el cielo
        azar = self.azar
        for _ in range(int(max(ancho - x_inicio, 0) * alto / densidad)):
            x, y = azar.uniform(x_inicio, ancho), azar.uniform(0, alto)
            profundidad = azar.choice((0.3, 0.3, 0.3, 0.55, 0.55, 1.0))
            r = {0.3: 0.6, 0.55: 0.85, 1.0: 1.25}[profundidad] * azar.uniform(0.85, 1.2)
            nivel = min(255, int(azar.randint(60, 140) + 110 * profundidad))
            tono = azar.choice(((1, 1, 1), (0.8, 0.88, 1), (1, 0.93, 0.8)))
            color = '#%02x%02x%02x' % tuple(int(nivel * t) for t in tono)
            objeto = lienzo.create_oval(x - r, y - r, x + r, y + r, fill=color, outline='', tags=('cielo',))
            self.estrellas.append([objeto, x, y, profundidad, color, 0.0, 0.0])
        if destellos:
            for _ in range(max(2, int(ancho * alto / 120000))):
                x, y = azar.uniform(x_inicio, ancho), azar.uniform(0, alto)
                color = azar.choice(('#cfe3ff', '#fff1cf'))
                objeto = dibujar_destello(lienzo, x, y, azar.uniform(4, 7), color, etiquetas=('cielo',))
                self.estrellas.append([objeto, x, y, 1.3, color, 0.0, 0.0])

    def _desvio_lente(self, x, y):
        """Cuánto desplaza la lente gravitatoria una estrella que está en (x, y).
        Ecuación de una lente puntual: la imagen se ve a una distancia
        θ = (β + √(β² + 4·θE²)) / 2 del centro, siendo β la distancia real."""
        if self.lente is None:
            return 0.0, 0.0
        cx, cy, theta_e = self.lente
        dx, dy = x - cx, y - cy
        beta = math.hypot(dx, dy)
        if beta < 1e-6 or theta_e <= 0:
            return 0.0, 0.0
        theta = (beta + math.sqrt(beta * beta + 4 * theta_e * theta_e)) / 2
        factor = theta / beta - 1
        return dx * factor, dy * factor

    def mover(self, dt):
        lienzo = self.lienzo
        salto = self.ancho - self.x_inicio + 16
        for estrella in self.estrellas:
            objeto, x, y, profundidad, _color, ox, oy = estrella
            nueva_x = x - self.velocidad * profundidad * dt
            nueva_y = y
            if nueva_x < self.x_inicio - 8:  # sale por la izquierda: vuelve a entrar por la derecha
                nueva_x += salto
                nueva_y = self.azar.uniform(0, self.alto)
            nox, noy = self._desvio_lente(nueva_x, nueva_y)
            lienzo.move(objeto, (nueva_x + nox) - (x + ox), (nueva_y + noy) - (y + oy))
            estrella[1], estrella[2], estrella[5], estrella[6] = nueva_x, nueva_y, nox, noy
        # titileo: unas pocas estrellas cambian de brillo en cada fotograma
        for _ in range(max(1, len(self.estrellas) // 80) if self.titileo else 0):
            estrella = self.azar.choice(self.estrellas)
            objeto, profundidad, color = estrella[0], estrella[3], estrella[4]
            if profundidad < 1.3:
                if self.azar.random() < 0.5:
                    nuevo = mezclar(color, '#000000', self.azar.uniform(0.2, 0.7))
                else:
                    nuevo = mezclar(color, '#ffffff', self.azar.uniform(0.0, 0.5))
                lienzo.itemconfig(objeto, fill=nuevo)

    def orbitar(self, agujero, dt):
        """Las estrellas cercanas a un agujero negro giran a su alrededor.

        Cuanto más cerca, más rápido giran (como en las órbitas de Kepler,
        velocidad angular proporcional a 1/r^1,5), van cayendo poco a poco
        hacia dentro y las que llegan al horizonte desaparecen ("se las traga")
        y vuelven a aparecer por la derecha.
        """
        cx, cy = agujero.x, agujero.y
        r_horizonte, r_influencia = agujero.r_horizonte, agujero.r_influencia
        for estrella in self.estrellas:
            objeto, x, y = estrella[0], estrella[1], estrella[2]
            ex, ey = x - cx, y - cy
            r = math.hypot(ex, ey)
            if r >= r_influencia or r == 0:
                continue
            if r < r_horizonte * 1.1:  # tragada: reaparece por la derecha
                nueva_x, nueva_y = self.ancho + self.azar.uniform(2, 10), self.azar.uniform(0, self.alto)
            else:
                f = min(1.0, (r_influencia - r) / (r_influencia * 0.35))  # suaviza el borde de la zona
                omega = min(3.0, agujero.omega * (3 * r_horizonte / r) ** 1.5 * f)
                angulo = -omega * dt  # sentido antihorario en pantalla
                nuevo_r = r * (1 - 0.05 * f * dt * (3 * r_horizonte / r))  # espiral hacia dentro
                c, s = math.cos(angulo), math.sin(angulo)
                nueva_x = cx + nuevo_r * (ex * c - ey * s) / r + agujero.vx * dt * f * 0.7
                nueva_y = cy + nuevo_r * (ex * s + ey * c) / r
            self.lienzo.move(objeto, nueva_x - x, nueva_y - y)
            estrella[1], estrella[2] = nueva_x, nueva_y


class Galaxia:
    """Una galaxia lejana muy tenue (elipse inclinada) que cruza el cielo despacio."""

    def __init__(self, lienzo, ancho, alto, azar, velocidad):
        self.lienzo, self.ancho, self.velocidad = lienzo, ancho, velocidad
        self.tamano = azar.uniform(60, 90)
        self.x = azar.uniform(ancho * 0.6, ancho * 0.9)
        y = alto * azar.choice((0.18, 0.82))
        inclinacion = math.radians(azar.uniform(20, 40))
        capas = ((1.0, 0.10), (0.75, 0.20), (0.5, 0.32), (0.3, 0.48), (0.14, 0.7))
        for escala, intensidad in capas:
            a, b = self.tamano * escala, self.tamano * escala * 0.32
            puntos = []
            for k in range(36):
                ang = 2 * math.pi * k / 36
                px, py = a * math.cos(ang), b * math.sin(ang)
                puntos += [self.x + px * math.cos(inclinacion) - py * math.sin(inclinacion),
                           y + px * math.sin(inclinacion) + py * math.cos(inclinacion)]
            lienzo.create_polygon(*puntos, fill=mezclar(ESPACIO, '#9aa3ff', intensidad), outline='',
                                  smooth=True, tags=('galaxia',))
        lienzo.create_oval(self.x - 2, y - 2, self.x + 2, y + 2, fill='#f2ecff', outline='', tags=('galaxia',))

    def mover(self, dt):
        dx = -self.velocidad * dt
        self.x += dx
        if self.x < -self.tamano * 1.2:
            dx += self.ancho + self.tamano * 2.4
            self.x += self.ancho + self.tamano * 2.4
        self.lienzo.move('galaxia', dx, 0)


class EstrellaFugaz:
    """Un meteoro: una cabeza brillante con una estela que se apaga."""
    SEGMENTOS = 10

    def __init__(self, lienzo, ancho, alto, azar):
        self.lienzo = lienzo
        self.x = azar.uniform(ancho * 0.3, ancho * 1.05)
        self.y = azar.uniform(-10, alto * 0.45)
        angulo = math.radians(azar.uniform(18, 38))
        velocidad = azar.uniform(520, 820)
        self.ux, self.uy = -math.cos(angulo), math.sin(angulo)  # dirección del movimiento
        self.vx, self.vy = self.ux * velocidad, self.uy * velocidad
        self.largo = azar.uniform(90, 170)
        self.vida = azar.uniform(0.8, 1.3)
        self.edad = 0.0
        self.color = azar.choice(('#ffffff', '#d9e6ff', '#fff3d6'))
        self.trozos = []
        for i in range(self.SEGMENTOS):
            ancho_linea = 2.4 - 1.6 * i / self.SEGMENTOS
            self.trozos.append(lienzo.create_line(0, 0, 0, 0, fill=ESPACIO, width=ancho_linea,
                                                  capstyle=tk.ROUND, tags=('fugaz',)))
        self.cabeza = lienzo.create_oval(0, 0, 0, 0, fill=ESPACIO, outline='', tags=('fugaz',))
        for objeto in self.trozos + [self.cabeza]:
            lienzo.tag_lower(objeto, 'frente')  # por detrás de los textos

    def paso(self, dt):
        """Avanza la animación. Devuelve False cuando ya ha desaparecido."""
        self.edad += dt
        if self.edad >= self.vida:
            for objeto in self.trozos + [self.cabeza]:
                self.lienzo.delete(objeto)
            return False
        self.x += self.vx * dt
        self.y += self.vy * dt
        visible = min(1.0, self.edad / 0.12) * min(1.0, (self.vida - self.edad) / 0.4)
        for i, objeto in enumerate(self.trozos):
            t0, t1 = i / self.SEGMENTOS, (i + 1) / self.SEGMENTOS
            self.lienzo.coords(objeto,
                               self.x - self.ux * self.largo * t0, self.y - self.uy * self.largo * t0,
                               self.x - self.ux * self.largo * t1, self.y - self.uy * self.largo * t1)
            self.lienzo.itemconfig(objeto, fill=mezclar(ESPACIO, self.color, visible * (1 - t0) ** 1.6))
        r = 2.0
        self.lienzo.coords(self.cabeza, self.x - r, self.y - r, self.x + r, self.y + r)
        self.lienzo.itemconfig(self.cabeza, fill=mezclar(ESPACIO, '#ffffff', visible))
        return True


class Satelite:
    """Un puntito que cruza el cielo en línea recta y despacio, como un satélite."""

    def __init__(self, lienzo, ancho, alto, azar):
        self.lienzo, self.ancho, self.alto = lienzo, ancho, alto
        if azar.random() < 0.5:
            self.x, self.vx = -5.0, azar.uniform(35, 60)
        else:
            self.x, self.vx = ancho + 5.0, -azar.uniform(35, 60)
        self.y = azar.uniform(alto * 0.1, alto * 0.9)
        self.vy = azar.uniform(-12, 12)
        self.objeto = lienzo.create_oval(self.x - 1.2, self.y - 1.2, self.x + 1.2, self.y + 1.2,
                                         fill='#d8dcef', outline='', tags=('satelite',))
        lienzo.tag_lower(self.objeto, 'frente')

    def paso(self, dt):
        dx, dy = self.vx * dt, self.vy * dt
        self.x += dx
        self.y += dy
        if self.x < -10 or self.x > self.ancho + 10:
            self.lienzo.delete(self.objeto)
            return False
        self.lienzo.move(self.objeto, dx, dy)
        return True


def png_rgba(rgba, compresion=6):
    """Escribe a mano un archivo PNG (en memoria) a partir de una matriz numpy
    (alto, ancho, 4) de bytes RGBA: cabecera + datos comprimidos con zlib."""
    alto, ancho, _ = rgba.shape
    filas = np.concatenate([np.zeros((alto, 1), np.uint8), rgba.reshape(alto, ancho * 4)], axis=1)

    def trozo(tipo, datos):
        return (struct.pack('>I', len(datos)) + tipo + datos
                + struct.pack('>I', zlib.crc32(tipo + datos) & 0xFFFFFFFF))

    return (b'\x89PNG\r\n\x1a\n'
            + trozo(b'IHDR', struct.pack('>IIBBBBB', ancho, alto, 8, 6, 0, 0, 0))
            + trozo(b'IDAT', zlib.compress(filas.tobytes(), compresion))
            + trozo(b'IEND', b''))


def foto_rgba(rgba):
    """Imagen de tkinter CON transparencia a partir de una matriz RGBA
    (tkinter sabe leer PNG, así que se le pasa un PNG hecho a mano)."""
    return tk.PhotoImage(data=base64.b64encode(png_rgba(rgba)).decode('ascii'), format='png')


def ruido_fractal(alto, ancho, generador, celdas=(3, 6, 12, 24, 48)):
    """Ruido suave tipo "nubes": suma de rejillas aleatorias cada vez más finas,
    interpoladas suavemente. Devuelve valores entre 0 y 1."""
    total = np.zeros((alto, ancho))
    peso_total = 0.0
    for i, n in enumerate(celdas):
        rejilla = generador.random((n + 1, n + 1))
        ys, xs = np.linspace(0, n, alto), np.linspace(0, n, ancho)
        y0, x0 = np.minimum(ys.astype(int), n - 1), np.minimum(xs.astype(int), n - 1)
        fy, fx = (ys - y0)[:, None], (xs - x0)[None, :]
        fy, fx = fy * fy * (3 - 2 * fy), fx * fx * (3 - 2 * fx)  # transición suave
        valor = (rejilla[y0][:, x0] * (1 - fx) * (1 - fy) + rejilla[y0][:, x0 + 1] * fx * (1 - fy)
                 + rejilla[y0 + 1][:, x0] * (1 - fx) * fy + rejilla[y0 + 1][:, x0 + 1] * fx * fy)
        peso = 0.55 ** i
        total += peso * valor
        peso_total += peso
    return total / peso_total


class NebulosaOrion:
    """Easter egg (teclas M-4-2): la Nebulosa de Orión (M42) cruzando el cielo.

    La imagen se genera con numpy: gas de hidrógeno rojo-rosado por fuera,
    zona verde-azulada y blanca en el centro, una lengua oscura de polvo,
    la pequeña M43 encima y las 4 estrellas del Trapecio en el corazón.
    Todo con transparencia, así se ven las estrellas y el cielo a través.
    """

    def __init__(self, lienzo, ancho, alto, texto='M42', velocidad=26):
        self.lienzo, self.velocidad = lienzo, velocidad
        self.r = r = min(alto * 0.34, 180)
        self.x = x = ancho + r * 0.1  # asoma ya por la derecha
        y = alto * 0.42
        lado = int(r * 2.4)
        generador = np.random.default_rng(42)
        v, u = np.mgrid[-1.2:1.2:lado * 1j, -1.2:1.2:lado * 1j]  # coordenadas en "radios"
        nubes = ruido_fractal(lado, lado, generador)
        nubes2 = ruido_fractal(lado, lado, generador)
        # hidrógeno: forma general irregular
        hidrogeno = np.exp(-(u ** 2 / 0.62 ** 2 + (v - 0.05) ** 2 / 0.55 ** 2)) * (0.25 + 1.3 * nubes)
        hidrogeno = np.clip(hidrogeno * 1.25 - 0.3, 0, None) ** 1.3
        hidrogeno += 0.55 * np.exp(-((u - 0.12) ** 2 + (v + 0.8) ** 2) / 0.16 ** 2) * (0.4 + nubes2)  # M43
        # oxígeno (verde-azulado) y núcleo blanco alrededor del Trapecio
        cx0, cy0 = 0.08, -0.04
        d2 = (u - cx0) ** 2 + (v - cy0) ** 2
        oxigeno = np.exp(-d2 / 0.2 ** 2) * (0.5 + 0.8 * nubes2)
        nucleo = np.exp(-d2 / 0.07 ** 2)
        # polvo oscuro que se mete desde la izquierda ("la boca del pez")
        polvo = 1 - 0.9 * np.exp(-(((u + 0.32) / 0.25) ** 2 + ((v + 0.1) / 0.11) ** 2))
        polvo *= 1 - 0.5 * np.exp(-(((u - 0.35) / 0.12) ** 2 + ((v - 0.3) / 0.2) ** 2)) * nubes
        emision = (hidrogeno[..., None] * np.array([235, 75, 110]) * 0.9
                   + oxigeno[..., None] * np.array([120, 225, 205]) * 0.8
                   + nucleo[..., None] * np.array([255, 250, 238])) * polvo[..., None]
        brillo = emision.max(axis=2) / 255.0
        alfa = np.clip(brillo * 1.15, 0, 0.92)
        color = np.clip(emision / np.maximum(brillo, 1e-6)[..., None], 0, 255)
        rgba = np.dstack([color, alfa * 255]).astype(np.uint8)
        self.foto = foto_rgba(rgba)  # hay que guardar la imagen o tkinter la borra
        lienzo.create_image(x, y, image=self.foto, tags=('m42',))
        azar = random.Random(42)
        # el Trapecio: 4 estrellas jóvenes muy brillantes
        for dx, dy in ((-0.035, -0.02), (0.03, -0.03), (-0.02, 0.03), (0.04, 0.025)):
            sx, sy = x + 0.08 * r + dx * r, y - 0.04 * r + dy * r
            lienzo.create_oval(sx - 2, sy - 2, sx + 2, sy + 2, fill='#ffffff', outline='', tags=('m42',))
        # algunas estrellas dentro de la nebulosa
        for _ in range(25):
            ang, dist = azar.uniform(0, 2 * math.pi), r * azar.uniform(0.1, 0.9) ** 0.5
            sx, sy = x + dist * math.cos(ang), y + dist * 0.85 * math.sin(ang)
            t = azar.uniform(0.6, 1.4)
            lienzo.create_oval(sx - t, sy - t, sx + t, sy + t, fill='#f4ecff', outline='', tags=('m42',))
        lienzo.create_text(x, y + r * 0.75, text=texto,fill=mezclar(ESPACIO, '#f0a0b4', 0.7),
                           font=('Segoe UI', 10, 'italic'), tags=('m42',))
        lienzo.tag_lower('m42', 'cielo')  # por detrás de las estrellas del fondo

    def paso(self, dt, _t=0.0):
        dx = -self.velocidad * dt
        self.x += dx
        self.lienzo.move('m42', dx, 0)
        if self.x < -self.r * 1.4:
            self.lienzo.delete('m42')
            return False
        return True


class AgujeroNegro:
    """Easter egg (teclas 6-1-8): TON 618, uno de los agujeros negros más
    masivos conocidos (unos 66.000 millones de soles), cruzando el
    cielo. Tiene un horizonte de sucesos negro, un anillo de fotones, un disco
    de acreción incandescente (con la parte de detrás curvada por encima por
    la gravedad) y hace orbitar a las estrellas del fondo (ver CampoEstelar.orbitar)."""

    def __init__(self, lienzo, ancho, alto, velocidad=20):
        self.lienzo = lienzo
        self.r_horizonte = rh = max(18.0, min(alto * 0.06, 32.0))
        self.r_influencia = rh * 8
        self.omega = 0.9               # velocidad de giro (rad/s) a 3 radios del horizonte
        self.vx = -velocidad
        self.x = ancho + rh * 0.5  # asoma ya por la derecha
        self.y = max(rh * 3.0, alto * 0.2)
        self.a = rh * 3.3              # semieje del disco de acreción
        self.b = self.a * 0.26
        x, y = self.x, self.y
        etiquetas = ('agujero',)
        # resplandor anaranjado difuso alrededor (imagen con transparencia, detrás de las estrellas)
        lado = int(rh * 9)
        v, u = np.mgrid[-1:1:lado * 1j, -1:1:lado * 1j]
        d = np.hypot(u, v * 1.6) * lado / 2 / rh  # distancia en radios del horizonte, achatado
        alfa = 0.55 * np.exp(-((d - 1.2) / 1.6) ** 2)
        rgba = np.zeros((lado, lado, 4), np.uint8)
        rgba[..., 0], rgba[..., 1], rgba[..., 2] = 255, 140, 60
        rgba[..., 3] = (np.clip(alfa, 0, 1) * 255).astype(np.uint8)
        self.foto = foto_rgba(rgba)
        lienzo.create_image(x, y, image=self.foto, tags=('agujero', 'agujero_fondo'))
        lienzo.tag_lower('agujero_fondo', 'cielo')
        colores_disco = ('#5a1d06', '#8f3410', '#c5561a', '#e8832c', '#ffb85c', '#ffe3a3')
        semiejes = [self.a * (1 - 0.1 * i) for i in range(len(colores_disco))]
        # disco completo (la mitad de delante se vuelve a dibujar después, encima del agujero)
        for semieje, color in zip(semiejes, colores_disco):
            sb = semieje * 0.26
            lienzo.create_oval(x - semieje, y - sb, x + semieje, y + sb, outline=color, width=3, tags=etiquetas)
        # luz del disco de detrás, curvada por encima del agujero (lente gravitatoria)
        for factor, color, grosor in ((1.75, '#b8501a', 2), (1.5, '#f0a050', 3), (1.25, '#ffe0a0', 2)):
            r = rh * factor
            lienzo.create_arc(x - r, y - r, x + r, y + r, start=0, extent=180, style=tk.ARC,
                              outline=color, width=grosor, tags=etiquetas)
        # horizonte de sucesos y anillo de fotones
        lienzo.create_oval(x - rh * 1.06, y - rh * 1.06, x + rh * 1.06, y + rh * 1.06,
                           outline='#fff2cf', width=2, tags=etiquetas)
        lienzo.create_oval(x - rh, y - rh, x + rh, y + rh, fill='#000000', outline='', tags=etiquetas)
        # mitad delantera del disco, por delante del agujero
        for semieje, color in zip(semiejes, colores_disco):
            sb = semieje * 0.26
            lienzo.create_arc(x - semieje, y - sb, x + semieje, y + sb, start=180, extent=180,
                              style=tk.ARC, outline=color, width=3, tags=etiquetas)
        # grumos calientes que giran con el disco
        self.grumos = []
        for i in range(7):
            objeto = lienzo.create_oval(0, 0, 0, 0, fill='#fff4d6', outline='', tags=etiquetas)
            self.grumos.append((objeto, 2 * math.pi * i / 7, random.uniform(0.62, 0.95)))
        self.texto = lienzo.create_text(x, y - rh * 2.3, text='TON 618',
                                        fill='#d9a070', font=('Segoe UI', 10, 'italic'), tags=etiquetas)
        for objeto in lienzo.find_withtag('agujero'):
            if 'agujero_fondo' not in lienzo.gettags(objeto):
                lienzo.tag_lower(objeto, 'frente')  # por delante de las estrellas, detrás de los textos
        self.t = 0.0

    def paso(self, dt, _t=0.0):
        self.t += dt
        dx = self.vx * dt
        self.x += dx
        self.lienzo.move('agujero', dx, 0)
        # los grumos del disco dan vueltas (más rápido cuanto más cerca del centro)
        for objeto, fase, radio in self.grumos:
            ang = fase + self.t * 1.6 / radio ** 1.5
            gx = self.x + self.a * radio * math.cos(ang)
            gy = self.y + self.b * radio * math.sin(ang)
            detras = math.sin(ang) < 0 and abs(gx - self.x) < self.r_horizonte
            r = 0 if detras else 1.8
            self.lienzo.coords(objeto, gx - r, gy - r, gx + r, gy + r)
        if self.x < -self.a * 1.3:
            self.lienzo.delete('agujero')
            return False
        return True


# -----------------------------------------------------------------------------
# Más easter eggs
# -----------------------------------------------------------------------------

def contorno_irregular(cx, cy, radio, azar, irregularidad=0.3, puntos=36, limite=None):
    """Contorno de una "mancha" irregular (continentes, nubes...). Si se da
    `limite` = (centro_x, centro_y, radio_max), no se sale de ese círculo."""
    ondas = [(azar.randint(2, 6), azar.uniform(0, 2 * math.pi), azar.uniform(0.3, 1.0)) for _ in range(4)]
    total = sum(a for _f, _p, a in ondas)
    coordenadas = []
    for k in range(puntos):
        ang = 2 * math.pi * k / puntos
        r = radio * (1 + irregularidad * sum(a * math.sin(f * ang + p) for f, p, a in ondas) / total)
        x, y = cx + r * math.cos(ang), cy + r * math.sin(ang)
        if limite:
            lx, ly, lr = limite
            d = math.hypot(x - lx, y - ly)
            if d > lr:
                x, y = lx + (x - lx) * lr / d, ly + (y - ly) * lr / d
        coordenadas += [x, y]
    return coordenadas


class Pulsar:
    """Púlsar diminuto en la esquina de la cabecera: parpadea con un ritmo
    fijo (pulso principal y otro más débil, como el púlsar del Cangrejo)."""
    PERIODO = 0.7

    def __init__(self, lienzo, x, y):
        self.lienzo, self.t = lienzo, 0.0
        self.haces = [lienzo.create_line(x - 14, y + 9, x + 14, y - 9, fill='#060a18', width=1),
                      lienzo.create_line(x - 5, y - 13, x + 5, y + 13, fill='#060a18', width=1)]
        self.brillo = lienzo.create_oval(x - 6, y - 6, x + 6, y + 6, fill='#060a18', outline='')
        self.nucleo = lienzo.create_oval(x - 2, y - 2, x + 2, y + 2, fill='#8fa8e0', outline='')

    def paso(self, dt):
        self.t += dt
        fase = (self.t % self.PERIODO) / self.PERIODO
        pulso = math.exp(-(min(fase, 1 - fase) / 0.05) ** 2) + 0.45 * math.exp(-((fase - 0.42) / 0.04) ** 2)
        self.lienzo.itemconfig(self.nucleo, fill=mezclar('#6f84c0', '#ffffff', pulso))
        self.lienzo.itemconfig(self.brillo, fill=mezclar('#060a18', '#5f8cff', pulso * 0.8))
        for haz in self.haces:
            self.lienzo.itemconfig(haz, fill=mezclar('#060a18', '#cfe4ff', pulso))
        return True


class Supernova:
    """Easter egg M1 / CRAB: fogonazo de supernova, una onda de luz que barre
    la pantalla de lado a lado y, al final, la Nebulosa del Cangrejo y un
    púlsar parpadeando en la esquina de la cabecera."""

    def __init__(self, lienzo, ancho, alto, texto, al_terminar):
        self.lienzo, self.ancho, self.alto, self.t = lienzo, ancho, alto, 0.0
        self.al_terminar, self.pulsar_lanzado = al_terminar, False
        self.x0, self.y0 = ancho * 0.74, alto * 0.3
        self.R = R = min(alto * 0.19, 110)
        # Restos de la supernova (se generan ya, pero quedan tapados por el fogonazo)
        lado = int(R * 2.4)
        generador = np.random.default_rng(1054)  # año en que se vio la supernova
        v, u = np.mgrid[-1.2:1.2:lado * 1j, -1.2:1.2:lado * 1j]
        forma = np.exp(-1.4 * ((u / 1.0) ** 2 + (v / 0.72) ** 2) ** 2)
        filamentos = (1 - np.abs(2 * ruido_fractal(lado, lado, generador, (4, 8, 16, 32, 64)) - 1)) ** 6
        nubes = ruido_fractal(lado, lado, generador)
        rojo = forma * (0.15 + 2.2 * filamentos) * (0.5 + 0.7 * nubes)
        azul = np.exp(-(u ** 2 + (v / 0.8) ** 2) / 0.33 ** 2) * 0.75 * (0.5 + 0.6 * nubes)
        emision = rojo[..., None] * np.array([255, 100, 45]) + azul[..., None] * np.array([140, 185, 255])
        brillo = emision.max(axis=2) / 255.0
        color = np.clip(emision / np.maximum(brillo, 1e-6)[..., None], 0, 255)
        alfa = np.clip(brillo * 1.2, 0, 0.93)
        self.foto = foto_rgba(np.dstack([color, alfa * 255]).astype(np.uint8))
        lienzo.create_image(self.x0, self.y0, image=self.foto, tags=('crab',))
        lienzo.create_text(self.x0, self.y0 + R * 0.95, text=texto, fill=mezclar(ESPACIO, '#ffb08a', 0.75),
                           font=('Segoe UI', 10, 'italic'), tags=('crab',))
        lienzo.tag_lower('crab', 'cielo')
        # Ondas de luz que barren la pantalla hacia los dos lados (franjas finas en degradado)
        self.bandas = []
        for lado_onda in (-1, 1):
            for i in range(14):
                color = mezclar('#ffffff', '#122040', (i / 13) ** 0.8)
                objeto = lienzo.create_rectangle(0, 0, 0, 0, fill=color, outline='', tags=('sn',))
                self.bandas.append((objeto, lado_onda, i))
        # Fogonazo: halo, estrella y puntas de luz
        self.halos = [lienzo.create_oval(0, 0, 0, 0, fill=ESPACIO, outline='', tags=('sn',)) for _ in range(3)]
        self.estrella = lienzo.create_oval(0, 0, 0, 0, fill='#ffffff', outline='', tags=('sn',))
        self.puntas = [lienzo.create_line(0, 0, 0, 0, fill='#ffffff', width=2, tags=('sn',)) for _ in range(2)]

    def paso(self, dt):
        c, t = self.lienzo, self.t + dt
        self.t = t
        # 1) La estrella explota: crece muy deprisa y luego se apaga
        if t < 2.2:
            r = 2 + 18 * min(t / 0.4, 1) ** 0.5
            apagado = max(0.0, (t - 0.6) / 1.6)
            color = mezclar('#ffffff', ESPACIO, apagado ** 0.7)
            c.coords(self.estrella, self.x0 - r, self.y0 - r, self.x0 + r, self.y0 + r)
            c.itemconfig(self.estrella, fill=color)
            for k, halo in enumerate(self.halos):  # resplandor azulado alrededor
                rh = r * (1.8 + 1.1 * k)
                c.coords(halo, self.x0 - rh, self.y0 - rh, self.x0 + rh, self.y0 + rh)
                vida_halo = max(0.0, 1 - t / 1.1) ** 2  # el resplandor dura poco
                c.itemconfig(halo, fill=mezclar(ESPACIO, '#9cc4ff', (0.4 - 0.12 * k) * vida_halo),
                             state=tk.NORMAL if vida_halo > 0.02 else tk.HIDDEN)
            largo = r * 5
            c.coords(self.puntas[0], self.x0 - largo, self.y0, self.x0 + largo, self.y0)
            c.coords(self.puntas[1], self.x0, self.y0 - largo, self.x0, self.y0 + largo)
            for punta in self.puntas:
                c.itemconfig(punta, fill=color)
            c.tag_raise('sn')
        elif self.estrella:
            c.delete('sn')
            self.estrella, self.bandas = None, []
        # 2) La onda de luz barre la pantalla de lado a lado
        if 0.15 < t < 2.2 and self.bandas:
            frente = (t - 0.15) * self.ancho * 0.75
            for objeto, lado_onda, i in self.bandas:
                x = self.x0 + lado_onda * (frente - i * 6)
                c.coords(objeto, x - 3, 0, x + 3, self.alto)
        # 3) El púlsar se queda en la esquina
        if t > 1.8 and not self.pulsar_lanzado:
            self.pulsar_lanzado = True
            self.al_terminar()
        # 4) Los restos se alejan despacio con el paneo
        dx = -6 * dt
        self.x0 += dx
        c.move('crab', dx, 0)
        if self.x0 < -self.R * 1.4:
            c.delete('crab')
            return False
        return True


def suave(a, b, x):
    """Transición suave de 0 a 1 entre a y b (si a > b, de 1 a 0)."""
    t = np.clip((x - a) / (b - a), 0, 1)
    return t * t * (3 - 2 * t)


def ruido_vertical(alto, ancho, generador, factor, celdas):
    """Ruido con vetas verticales: se calcula con menos filas y se estira hacia abajo."""
    corto = ruido_fractal(max(alto // factor, 4), ancho, generador, celdas)
    ys = np.linspace(0, corto.shape[0] - 1, alto)
    y0 = np.floor(ys).astype(int)
    y1 = np.minimum(y0 + 1, corto.shape[0] - 1)
    f = (ys - y0)[:, None]
    return corto[y0] * (1 - f) + corto[y1] * f


def difuminar(a, r):
    """Media en una ventana cuadrada de lado 2r+1 (con sumas acumuladas, muy rápido)."""
    if r < 1:
        return a
    p = np.pad(a, r, mode='edge')
    c = np.pad(np.cumsum(np.cumsum(p, axis=0), axis=1), ((1, 0), (1, 0)))
    n = 2 * r + 1
    return (c[n:, n:] - c[:-n, n:] - c[n:, :-n] + c[:-n, :-n]) / (n * n)


PILARES = [  # (x abajo, x en la cima, altura de la cima, anchura abajo, anchura arriba), en fracciones de la ventana
    (0.40, 0.335, 0.07, 0.075, 0.034),
    (0.535, 0.565, 0.24, 0.054, 0.026),
    (0.63, 0.655, 0.38, 0.036, 0.019),
]


def generar_pilares(ancho, alto, semilla=16):
    """Dibuja con numpy los Pilares de la Creación: gas de fondo, tres columnas
    de polvo con textura y vetas, iluminadas desde arriba, y gas brillante que
    se evapora de las cimas. Devuelve (color, transparencia)."""
    g = np.random.default_rng(semilla)
    escala = alto / 600
    v, u = np.mgrid[0:1:alto * 1j, 0:1:ancho * 1j]
    n_gas = ruido_fractal(alto, ancho, g, (3, 6, 12, 24, 48, 96))
    n_gas2 = ruido_fractal(alto, ancho, g, (2, 4, 8, 16))
    n_fino = ruido_fractal(alto, ancho, g, (12, 24, 48, 96, 192))
    n_det = ruido_fractal(alto, ancho, g, (40, 80, 160, 320))
    vetas = ruido_vertical(alto, ancho, g, 6, (4, 8, 16, 32, 64, 128))
    # Gas de fondo: azul verdoso arriba, cálido abajo a la izquierda, más oscuro a la derecha
    arriba, abajo, derecha = np.array([46, 110, 140]), np.array([170, 145, 80]), np.array([22, 40, 80])
    m_abajo = suave(0.45, 1.15, v * 0.85 + (1 - u) * 0.45)[..., None]
    m_der = (suave(0.55, 1.0, u) * 0.65)[..., None]
    fondo = arriba * (1 - m_abajo) + abajo * m_abajo
    fondo = fondo * (1 - m_der) + derecha * m_der
    fondo = fondo * (0.3 + 0.95 * n_gas + 0.25 * (n_det - 0.5))[..., None]
    # Silueta de los pilares (bordes deshilachados con ruido)
    dens = np.zeros_like(u)
    halo = np.zeros_like(u)
    for (xb, xc, vc, wb, wc) in PILARES:
        s = np.clip((v - vc) / (1 - vc), 0, 1)
        centro = xc + (xb - xc) * s + 0.010 * np.sin(v * 17 + xb * 40) + 0.015 * (n_gas2 - 0.5)
        anchura = (wc + (wb - wc) * s ** 0.7) * (1 + 0.5 * (vetas - 0.5))
        cabeza = 1 + 0.6 * np.exp(-((v - vc - 0.05) / 0.04) ** 2)
        d = np.abs(u - centro) / (anchura * cabeza) - 1 + 0.55 * (n_fino - 0.5) + 0.35 * (n_det - 0.5)
        cima = suave(vc - 0.03, vc + 0.05, v + 0.10 * (n_fino - 0.5) + 0.05 * (n_det - 0.5))
        dens = np.maximum(dens, suave(0.12, -0.12, d) * cima)
        cerca = np.exp(-(np.maximum(d, 0) / 1.0) ** 2)
        halo = np.maximum(halo, cerca * np.exp(-np.maximum(v - vc, 0) / 0.18) * suave(vc - 0.16, vc + 0.02, v))
    base = (suave(0.70, 0.88, v + 0.12 * (n_gas2 - 0.5) + 0.08 * (n_fino - 0.5))
            * suave(0.24, 0.33, u) * suave(0.74, 0.64, u))
    dens = np.maximum(dens, base)
    # Relieve: luz que viene de arriba (zonas junto al borde y superficies que miran hacia arriba)
    k = max(1, int(5 * escala))
    db = difuminar(dens, max(1, int(5 * escala)))
    contorno = np.clip((db - difuminar(dens, max(2, int(16 * escala)))) * 3.5, 0, 1) * dens
    hacia_arriba = np.clip((db - np.roll(db, k, axis=0)) * 6, 0, 1)
    de_arriba = 0.35 + 0.65 * np.exp(-np.clip(v - 0.07, 0, 1) / 0.35)
    luz = np.clip((contorno * 0.9 + hacia_arriba) * de_arriba + 0.25 * np.clip(vetas - 0.55, 0, 1), 0, 1)
    # Color del polvo: marrón muy oscuro con zonas rojizas, vetas y bordes claros
    textura = np.clip((n_fino - 0.5) * 2.4 + (vetas - 0.5) * 2.0 + 0.45, 0, 1)
    oscuro, medio, claro = np.array([22, 11, 11]), np.array([128, 56, 36]), np.array([245, 185, 130])
    color_p = oscuro + (medio - oscuro) * (textura ** 1.3)[..., None]
    color_p = color_p * (0.75 + 0.5 * n_det)[..., None]
    color_p = color_p + (claro - color_p) * (luz ** 1.2)[..., None] * 0.85
    color = fondo * (1 - dens[..., None]) + color_p * dens[..., None]
    color = color + (halo * (1 - dens * 0.8))[..., None] * np.array([120, 205, 215]) * 0.55
    color = np.clip(color, 0, 255)
    # Transparencia: el gas se funde con el espacio hacia los bordes de la ventana
    r = np.hypot((u - 0.5) / 0.6, (v - 0.5) / 0.7)
    alfa = np.clip(0.4 + 0.5 * n_gas, 0, 1) * np.clip(1.3 - r, 0, 1) ** 1.2
    alfa = np.maximum(alfa, np.clip(dens, 0, 1) * 0.98)
    return color, alfa


class PilaresCreacion:
    """Easter egg M16: primero el cielo se vacía de estrellas; después aparecen
    los Pilares de la Creación mientras nacen estrellas nuevas; al final la
    nube se desvanece y el cielo se vuelve a llenar poco a poco.
    La imagen (pesada de calcular) se prepara en segundo plano mientras el
    cielo se vacía, para que la animación no se congele."""
    VACIAR, ENTRADA, DURACION, SALIDA, LLENAR = 1.8, 2.5, 9.0, 3.0, 6.0
    NIVELES = 6

    def __init__(self, lienzo, ancho, alto, texto, cielo):
        self.lienzo, self.ancho, self.alto, self.cielo = lienzo, ancho, alto, cielo
        self.t, self.t_cortina = 0.0, None
        self.azar = random.Random()
        cielo.titileo = False
        self.orden = [(estrella, self.azar.random()) for estrella in cielo.estrellas]  # cuándo se apaga cada una
        self.niveles_png, self.niveles = None, []
        threading.Thread(target=self._preparar, daemon=True).start()
        self.cortina = lienzo.create_image(0, 0, anchor=tk.NW, state=tk.HIDDEN)
        self.texto = lienzo.create_text(ancho * 0.5, alto * 0.05, text=texto, state=tk.HIDDEN,
                                        fill='#f0d2a0', font=('Segoe UI', 11, 'italic'))
        lienzo.tag_lower(self.cortina, 'frente')
        lienzo.tag_lower(self.texto, 'frente')
        self.nacientes = []
        self.proximo_nacimiento = 0.5

    def _preparar(self):
        """En segundo plano: calcula la imagen y la convierte en PNG con varias transparencias."""
        try:
            color, alfa = generar_pilares(self.ancho, self.alto)
            niveles = []
            for k in range(1, self.NIVELES + 1):
                rgba = np.dstack([color, alfa * 255 * k / self.NIVELES]).astype(np.uint8)
                niveles.append(base64.b64encode(png_rgba(rgba, 3)).decode('ascii'))
            self.niveles_png = niveles
        except Exception:
            self.niveles_png = []

    def _brillo_cielo(self, f):
        """f = 1: todas las estrellas normales; f = 0: cielo vacío (se apagan/encienden escalonadas)."""
        for estrella, retraso in self.orden:
            g = min(1.0, max(0.0, (f - retraso * 0.6) / 0.4))
            self.lienzo.itemconfig(estrella[0], fill=mezclar(ESPACIO, estrella[4], g))

    def _nivel(self, i):
        while len(self.niveles) <= i:  # las imágenes se crean cuando se necesitan
            self.niveles.append(tk.PhotoImage(data=self.niveles_png[len(self.niveles)], format='png'))
        return self.niveles[i]

    def paso(self, dt):
        c, self.t = self.lienzo, self.t + dt
        t = self.t
        if self.t_cortina is None:
            self._brillo_cielo(max(0.0, 1 - t / self.VACIAR))
            if t >= self.VACIAR and self.niveles_png is not None:
                self.t_cortina = t
                if not self.niveles_png:  # no se ha podido generar: se vuelve a llenar el cielo
                    self.t_cortina = t - self.ENTRADA - self.DURACION - self.SALIDA
            return True
        tc = t - self.t_cortina
        fin_cortina = self.ENTRADA + self.DURACION
        if tc < self.ENTRADA:
            nivel = tc / self.ENTRADA
        elif tc < fin_cortina:
            nivel = 1.0
        else:
            nivel = max(0.0, 1 - (tc - fin_cortina) / self.SALIDA)
        indice = int(round(nivel * self.NIVELES)) - 1
        if indice < 0 or not self.niveles_png:
            c.itemconfig(self.cortina, state=tk.HIDDEN)
            c.itemconfig(self.texto, state=tk.HIDDEN)
        else:
            c.itemconfig(self.cortina, image=self._nivel(indice), state=tk.NORMAL)
            c.itemconfig(self.texto, state=tk.NORMAL)
        # Al final el cielo se vuelve a llenar poco a poco
        llenado = max(0.0, (tc - fin_cortina) / self.LLENAR)
        if tc >= fin_cortina:
            self._brillo_cielo(min(llenado, 1.0))
        # Nacen estrellas (sobre todo en las puntas de los pilares)
        if 0.6 < tc < fin_cortina and tc >= self.proximo_nacimiento and self.niveles_png:
            self.proximo_nacimiento = tc + self.azar.uniform(0.2, 0.5)
            self._nacer()
        for naciente in list(self.nacientes):
            naciente[4] += dt
            objeto, chispa, x, y, edad, color, radio = naciente
            # nace creciendo desde un punto (con su color), así se ve bien también sobre el gas brillante
            r = radio * min(1.0, edad / 1.0) ** 0.7
            if chispa:   # estrella pequeña con un destello que crece y se apaga
                c.coords(objeto, x - r, y - r, x + r, y + r)
                rc = 7 * math.sin(min(edad / 1.4, 1) * math.pi)
                c.coords(chispa, *puntos_destello(x, y, max(rc, 0.1)))
            else:        # estrella brillante con puntas
                c.coords(objeto, *puntos_destello(x, y, max(r, 0.1)))
            if edad >= 1.4:  # ya ha nacido: pasa a ser una estrella más del cielo
                if chispa:
                    c.delete(chispa)
                self.cielo.estrellas.append([objeto, x, y, 0.55, color, 0.0, 0.0])
                self.nacientes.remove(naciente)
        if llenado >= 1 and not self.nacientes:
            c.delete(self.cortina)
            c.delete(self.texto)
            self.cielo.titileo = True
            return False
        return True

    def _nacer(self):
        c, azar = self.lienzo, self.azar
        xb, xc, vc, wb, wc = azar.choice(PILARES)
        if azar.random() < 0.65:
            x = (xc + azar.uniform(-2.5, 2.5) * wc) * self.ancho
            y = (vc + azar.uniform(-0.08, 0.10)) * self.alto
        else:
            x, y = azar.uniform(0.1, 0.9) * self.ancho, azar.uniform(0.05, 0.9) * self.alto
        if azar.random() < 0.3:  # estrella brillante con puntas de luz, como en las fotos del Hubble
            color = azar.choice(('#ffd2a0', '#ffffff', '#ffbf86', '#ffe6c8'))
            radio = azar.uniform(4, 7)
            objeto = dibujar_destello(c, x, y, 0.1, color, etiquetas=('cielo',))
            chispa = None
        else:
            color = azar.choice(('#ffffff', '#dbe8ff', '#fff0d8', '#ffd9b0'))
            radio = azar.uniform(0.9, 1.6)
            objeto = c.create_oval(x, y, x, y, fill=color, outline='', tags=('cielo',))
            chispa = dibujar_destello(c, x, y, 0.1, color)
            c.tag_lower(chispa, 'frente')
        c.tag_lower(objeto, 'frente')
        self.nacientes.append([objeto, chispa, x, y, 0.0, color, radio])


def dibujar_voyager(c, x, y, s, etiquetas):
    """Sonda Voyager vista de lado, hecha con figuras simples (escala s)."""
    gris, oscuro = '#c8ccd6', '#6d7383'
    c.create_line(x - 4 * s, y - 6 * s, x - 70 * s, y - 60 * s, fill=oscuro, width=1, tags=etiquetas)  # magnetómetro
    c.create_line(x + 2 * s, y + 8 * s, x + 20 * s, y + 70 * s, fill=oscuro, width=1, tags=etiquetas)  # antenas
    c.create_line(x - 2 * s, y + 8 * s, x - 28 * s, y + 66 * s, fill=oscuro, width=1, tags=etiquetas)
    c.create_line(x - 6 * s, y + 4 * s, x - 42 * s, y + 26 * s, fill=gris, width=2, tags=etiquetas)    # brazo RTG
    for i in range(3):
        cx, cy = x - (30 + 5 * i) * s, y + (19 + 3 * i) * s
        c.create_rectangle(cx - 4 * s, cy - 3 * s, cx + 4 * s, cy + 3 * s, fill='#4b4f5a', outline='', tags=etiquetas)
    c.create_line(x + 6 * s, y + 6 * s, x + 40 * s, y + 22 * s, fill=gris, width=2, tags=etiquetas)    # brazo de cámaras
    c.create_rectangle(x + 36 * s, y + 17 * s, x + 48 * s, y + 27 * s, fill='#9aa0ac', outline='', tags=etiquetas)
    puntos = []
    for k in range(10):  # cuerpo decagonal dorado
        a = 2 * math.pi * k / 10
        puntos += [x + 13 * s * math.cos(a), y + 10 + 8 * s * math.sin(a)]
    c.create_polygon(*puntos, fill='#b48a32', outline='#e2bd62', tags=etiquetas)
    c.create_oval(x + 3 * s, y + 6 * s, x + 11 * s, y + 14 * s, fill='#f2cf6a', outline='#fff1b8', tags=etiquetas)  # disco de oro
    c.create_oval(x - 26 * s, y - 34 * s, x + 26 * s, y - 2 * s, fill='#e9ecf2', outline='#9aa2b4', tags=etiquetas)  # antena parabólica
    c.create_oval(x - 14 * s, y - 26 * s, x + 14 * s, y - 10 * s, fill='#d3d8e2', outline='', tags=etiquetas)
    c.create_line(x, y - 18 * s, x, y - 36 * s, fill=oscuro, width=2, tags=etiquetas)


class Voyager:
    """Easter egg 1977: la sonda Voyager cruza despacio la pantalla dejando
    ondas de radio detrás, con un sonido de espacio profundo."""

    def __init__(self, lienzo, ancho, alto, texto):
        self.lienzo, self.ancho, self.t = lienzo, ancho, 0.0
        self.x, self.y, self.s = -60.0, alto * 0.8, max(0.8, min(alto / 600, 1.3))
        dibujar_voyager(lienzo, self.x, self.y, self.s, ('voyager',))
        lienzo.create_text(self.x, self.y + 82 * self.s, text=texto, fill='#e2c67a',
                           font=('Segoe UI', 10, 'italic'), tags=('voyager',))
        lienzo.tag_lower('voyager', 'frente')
        self.ondas = []  # [objeto, radio]
        self.proxima_onda = 0.0
        reproducir_sonido_espacio()

    def paso(self, dt):
        c, self.t = self.lienzo, self.t + dt
        dx = 55 * dt  # cruza la pantalla en unos 20-25 segundos
        self.x += dx
        c.move('voyager', dx, 0)
        if self.t >= self.proxima_onda and self.x < self.ancho + 40:
            self.proxima_onda = self.t + 0.6
            objeto = c.create_arc(0, 0, 0, 0, start=120, extent=120, style=tk.ARC, outline='#7fd0ff', width=2)
            c.tag_lower(objeto, 'frente')
            self.ondas.append([objeto, 10.0, self.x, self.y - 18 * self.s])
        for onda in list(self.ondas):
            onda[1] += 45 * dt
            objeto, r, ox, oy = onda
            if r > 170:
                c.delete(objeto)
                self.ondas.remove(onda)
                continue
            c.coords(objeto, ox - r, oy - r, ox + r, oy + r)
            c.itemconfig(objeto, outline=mezclar('#7fd0ff', ESPACIO, r / 170))
        if self.x > self.ancho + 90 and not self.ondas:
            c.delete('voyager')
            return False
        return True


def reproducir_sonido_espacio(segundos=10.0):
    """Genera un sonido de "espacio profundo" (silbidos que bajan de tono, como
    las ondas de plasma que grabaron las Voyager, sobre un rumor grave) y lo
    reproduce sin parar el programa. Solo en Windows; si falla, no pasa nada."""
    try:
        import winsound
        muestras = 22050
        t = np.arange(int(segundos * muestras)) / muestras
        generador = np.random.default_rng()
        ruido = np.convolve(generador.normal(0, 1, t.size), np.ones(80) / 80, mode='same')
        senal = 0.5 * ruido / (ruido.std() + 1e-9) * 0.12
        senal += 0.10 * np.sin(2 * np.pi * 55 * t) * (0.6 + 0.4 * np.sin(2 * np.pi * 0.15 * t))
        for inicio in (0.6, 2.4, 4.7, 6.3, 8.0):  # silbidos ("whistlers")
            duracion = 1.7
            tramo = (t >= inicio) & (t < inicio + duracion)
            tt = t[tramo] - inicio
            frecuencia = 2600 * np.exp(-2.4 * tt) + 280
            fase = 2 * np.pi * np.cumsum(frecuencia) / muestras
            senal[tramo] += 0.32 * np.sin(np.pi * tt / duracion) ** 2 * np.sin(fase)
        fundido = np.minimum(1, np.minimum(t / 0.8, (segundos - t) / 1.5))
        senal = senal * fundido
        volumen = 0.16  # suave: de 0 (silencio) a 1 (máximo)
        datos = (np.clip(senal / (np.abs(senal).max() + 1e-9) * volumen, -1, 1) * 32767).astype('<i2')
        ruta = os.path.join(tempfile.gettempdir(), 'emillium_voyager.wav')
        with wave.open(ruta, 'wb') as archivo:
            archivo.setnchannels(1)
            archivo.setsampwidth(2)
            archivo.setframerate(muestras)
            archivo.writeframes(datos.tobytes())
        winsound.PlaySound(ruta, winsound.SND_FILENAME | winsound.SND_ASYNC | winsound.SND_NODEFAULT)
    except Exception:
        pass


def _guardar_wav(senal, muestras, nombre, volumen):
    """Normaliza la señal, la baja al volumen pedido y la guarda como WAV temporal."""
    datos = (np.clip(senal / (np.abs(senal).max() + 1e-9) * volumen, -1, 1) * 32767).astype('<i2')
    ruta = os.path.join(tempfile.gettempdir(), nombre)
    with wave.open(ruta, 'wb') as archivo:
        archivo.setnchannels(1)
        archivo.setsampwidth(2)
        archivo.setframerate(muestras)
        archivo.writeframes(datos.tobytes())
    return ruta


def sintetizar_zaratustra(muestras=22050):
    """Sintetiza el famoso comienzo de «Así habló Zaratustra» (Richard Strauss, 1896,
    obra de dominio público): pedal grave de órgano, trompetas Do-Sol-Do, acorde
    mayor-menor (y al revés), timbales... tres veces, hasta el gran acorde final.
    No es una grabación: cada nota se construye sumando sus armónicos."""
    duracion = 27.5
    t_total = np.arange(int(duracion * muestras)) / muestras
    senal = np.zeros_like(t_total)
    nota = {'C1': 32.70, 'C2': 65.41, 'G1': 49.00, 'C3': 130.81, 'F3': 174.61, 'G3': 196.00, 'A3': 220.00,
            'B3': 246.94, 'C4': 261.63, 'D4': 293.66, 'Eb4': 311.13, 'E4': 329.63, 'F4': 349.23, 'G4': 392.00,
            'A4': 440.00, 'B4': 493.88, 'C5': 523.25, 'D5': 587.33, 'Eb5': 622.25, 'E5': 659.25}
    trompeta = (1, .85, .7, .55, .45, .35, .25, .18, .12, .08)
    cuerdas = (1, .5, .33, .25, .15, .1, .06)

    def tocar(inicio, dur, frecuencias, timbre, amp, ataque=0.06, caida=0.35):
        i0 = int(inicio * muestras)
        n = int((dur + caida) * muestras)
        t = np.arange(n) / muestras
        envolvente = np.minimum(1, t / ataque) * np.where(t < dur, 1, np.exp(-(t - dur) / (caida / 3)))
        onda = np.zeros(n)
        for f in frecuencias:
            vibrato = 1 + 0.003 * np.sin(2 * np.pi * 5.2 * t) * np.minimum(1, t / 0.4)
            fase = 2 * np.pi * f * np.cumsum(vibrato) / muestras
            for k, a in enumerate(timbre, start=1):
                if f * k < muestras / 2.2:
                    onda += a * np.sin(k * fase)
        trozo = senal[i0:i0 + n]
        trozo += amp * (onda * envolvente)[:len(trozo)] / max(1, len(frecuencias)) ** 0.5

    def timbales(inicio, golpes=10, separacion=0.19):
        for g in range(golpes):
            i0 = int((inicio + g * separacion) * muestras)
            n = int(0.9 * muestras)
            t = np.arange(n) / muestras
            f = (65.41 if g % 2 == 0 else 49.00) * (1 + 0.04 * np.exp(-t / 0.05))  # el parche baja un poco de tono
            golpe = (np.sin(2 * np.pi * np.cumsum(f) / muestras) + 0.3 * np.sin(2 * np.pi * 2.76 * np.cumsum(f) / muestras))
            golpe *= np.exp(-t / 0.35) * (1.0 - 0.04 * g)
            trozo = senal[i0:i0 + n]
            trozo += 0.9 * golpe[:len(trozo)]

    mayor = ['C3', 'G3', 'C4', 'E4', 'G4', 'C5', 'E5']
    menor = ['C3', 'G3', 'C4', 'Eb4', 'G4', 'C5', 'Eb5']
    # pedal de órgano muy grave durante toda la pieza (crece al final)
    organo = np.minimum(1, t_total / 3.0) * (0.35 + 0.35 * np.clip((t_total - 22.5) / 1.0, 0, 1))
    organo *= np.where(t_total < 26.0, 1, np.exp(-(t_total - 26.0) / 0.5))
    senal += organo * (np.sin(2 * np.pi * 32.70 * t_total) + 0.6 * np.sin(2 * np.pi * 65.41 * t_total)
                       + 0.3 * np.sin(2 * np.pi * 130.81 * t_total))
    for inicio, primero, segundo in ((2.0, mayor, menor), (10.2, menor, mayor)):
        tocar(inicio, 1.25, [nota['C4']], trompeta, 1.0, ataque=0.12)
        tocar(inicio + 1.4, 1.25, [nota['G4']], trompeta, 1.0, ataque=0.10)
        tocar(inicio + 2.8, 0.85, [nota['C5']], trompeta, 1.1, ataque=0.06)
        tocar(inicio + 3.8, 0.35, [nota[n] for n in primero], cuerdas, 1.6, ataque=0.03, caida=0.15)
        tocar(inicio + 4.2, 1.4, [nota[n] for n in segundo], cuerdas, 1.4, ataque=0.03, caida=0.6)
        timbales(inicio + 5.8)
    inicio = 18.4  # tercera vez: sube hasta el gran acorde final
    tocar(inicio, 1.0, [nota['C4']], trompeta, 1.0, ataque=0.1)
    tocar(inicio + 1.2, 1.0, [nota['G4']], trompeta, 1.0, ataque=0.1)
    tocar(inicio + 2.4, 0.8, [nota['C5']], trompeta, 1.15, ataque=0.06)
    tocar(inicio + 3.4, 0.55, [nota[n] for n in ('F3', 'A3', 'C4', 'F4', 'A4', 'C5')], cuerdas, 1.5, ataque=0.04)
    tocar(inicio + 4.0, 0.55, [nota[n] for n in ('G3', 'B3', 'D4', 'G4', 'B4', 'D5')], cuerdas, 1.6, ataque=0.04)
    tocar(inicio + 4.6, 3.6, [nota[n] for n in mayor] + [nota['C2']], cuerdas + (0.04,), 2.0, ataque=0.05, caida=1.4)
    tocar(inicio + 4.6, 3.6, [nota['C5'], nota['E5']], trompeta, 0.9, ataque=0.05, caida=1.4)
    fundido = np.minimum(1, (duracion - t_total) / 1.0)
    return senal * np.clip(fundido, 0, 1), muestras


def reproducir_zaratustra(volumen=0.14):
    """Genera la música (en segundo plano, para no parar la animación) y la reproduce."""
    def tarea():
        try:
            import winsound
            senal, muestras = sintetizar_zaratustra()
            ruta = _guardar_wav(senal, muestras, 'emillium_2001.wav', volumen)
            winsound.PlaySound(ruta, winsound.SND_FILENAME | winsound.SND_ASYNC | winsound.SND_NODEFAULT)
        except Exception:
            pass
    threading.Thread(target=tarea, daemon=True).start()


def dibujar_discovery(c, x, y, s, etiquetas):
    """Nave de exploración (homenaje a «2001»): módulo esférico delante, espina
    larga y la zona de propulsión detrás. (x, y) es el centro de la esfera; la nave
    mira hacia la derecha."""
    blanco, gris, oscuro = '#e9ecf2', '#c4c9d4', '#7d8494'
    largo = 205 * s
    # espina con "vértebras"
    c.create_rectangle(x - largo, y - 1.6 * s, x - 14 * s, y + 1.6 * s, fill=gris, outline='', tags=etiquetas)
    for k in range(14):
        xv = x - 30 * s - k * 12.5 * s
        c.create_rectangle(xv - 2.2 * s, y - 3.4 * s, xv + 2.2 * s, y + 3.4 * s, fill=blanco, outline=oscuro,
                           tags=etiquetas)
    # cuello y módulo esférico con ventanas
    c.create_rectangle(x - 26 * s, y - 5 * s, x - 10 * s, y + 5 * s, fill=gris, outline=oscuro, tags=etiquetas)
    c.create_oval(x - 14 * s, y - 14 * s, x + 14 * s, y + 14 * s, fill=blanco, outline='#9aa2b4', tags=etiquetas)
    for k in range(5):
        a = math.radians(-35 + k * 17)
        wx, wy = x + 10.5 * s * math.cos(a), y + 10.5 * s * math.sin(a)
        c.create_oval(wx - 1.3 * s, wy - 1.3 * s, wx + 1.3 * s, wy + 1.3 * s, fill='#2b3346', outline='',
                      tags=etiquetas)
    # antena (mástil y plato) en la parte delantera de la espina
    c.create_line(x - 42 * s, y - 3 * s, x - 48 * s, y - 16 * s, fill=oscuro, width=max(1, int(1.5 * s)), tags=etiquetas)
    c.create_oval(x - 56 * s, y - 22 * s, x - 42 * s, y - 14 * s, fill=gris, outline=oscuro, tags=etiquetas)
    # zona de propulsión: bloque, depósitos y toberas
    c.create_rectangle(x - largo - 32 * s, y - 9 * s, x - largo, y + 9 * s, fill=gris, outline=oscuro, tags=etiquetas)
    for k in range(3):
        xd = x - largo + 3 * s - k * 10 * s
        c.create_rectangle(xd - 4 * s, y - 12 * s, xd + 4 * s, y + 12 * s, fill=blanco, outline=oscuro,
                           tags=etiquetas)
    for dy in (-6, 0, 6):
        c.create_polygon(x - largo - 32 * s, y + (dy - 2.2) * s, x - largo - 42 * s, y + (dy - 3.6) * s,
                         x - largo - 42 * s, y + (dy + 3.6) * s, x - largo - 32 * s, y + (dy + 2.2) * s,
                         fill=oscuro, outline='', tags=etiquetas)
    return largo + 42 * s  # longitud total de la nave


class Odisea2001:
    """Easter egg 2001: la nave cruza despacio (a la velocidad de la Voyager)
    mientras suena el comienzo de «Así habló Zaratustra»."""

    def __init__(self, lienzo, ancho, alto, texto):
        self.lienzo, self.ancho = lienzo, ancho
        self.s = max(0.8, min(alto / 600, 1.3))
        self.x, self.y = -10.0, alto * 0.17
        self.largo = dibujar_discovery(lienzo, self.x, self.y, self.s, ('odisea',))
        lienzo.create_text(self.x - self.largo / 2, self.y + 26 * self.s, text=texto, fill='#c9d3f0',
                           font=('Segoe UI', 10, 'italic'), tags=('odisea',))
        lienzo.tag_lower('odisea', 'frente')
        reproducir_zaratustra()

    def paso(self, dt):
        dx = 55 * dt  # igual que la Voyager
        self.x += dx
        self.lienzo.move('odisea', dx, 0)
        if self.x - self.largo > self.ancho + 20:
            self.lienzo.delete('odisea')
            return False
        return True


def dibujar_iss(c, x, y, s, etiquetas):
    """Estación Espacial Internacional vista desde abajo (escala s)."""
    panel, borde = '#b7832f', '#5e4316'
    c.create_line(x - 80 * s, y, x + 80 * s, y, fill='#d9dde6', width=max(2, int(3 * s)), tags=etiquetas)  # viga
    for px in (-70, -52, 52, 70):  # 8 paneles solares, 4 a cada lado
        for lado in (-1, 1):
            y0, y1 = (y + lado * 4 * s, y + lado * 34 * s)
            c.create_rectangle(x + (px - 7) * s, min(y0, y1), x + (px + 7) * s, max(y0, y1),
                               fill=panel, outline=borde, tags=etiquetas)
    for dx, dy, w, h in ((0, 0, 9, 30), (0, 0, 26, 7), (-12, -10, 5, 9), (14, 10, 5, 9)):  # módulos
        c.create_rectangle(x + (dx - w / 2) * s, y + (dy - h / 2) * s, x + (dx + w / 2) * s, y + (dy + h / 2) * s,
                           fill='#eef1f6', outline='#9aa2b4', tags=etiquetas)
    for dx in (-26, 26):  # radiadores
        c.create_rectangle(x + (dx - 3) * s, y + 4 * s, x + (dx + 3) * s, y + 20 * s, fill='#f5f7fb',
                           outline='', tags=etiquetas)


class EstacionEspacial:
    """Easter egg ISS: la estación cruza el cielo de izquierda a derecha muy rápido."""

    def __init__(self, lienzo, ancho, alto, texto):
        self.lienzo, self.ancho = lienzo, ancho
        self.x, self.y = -120.0, alto * random.uniform(0.12, 0.25)
        escala = max(0.7, min(alto / 650, 1.2))
        dibujar_iss(lienzo, self.x, self.y, escala, ('iss',))
        lienzo.create_text(self.x, self.y + 48 * escala, text=texto, fill='#d6dcea',
                           font=('Segoe UI', 10, 'italic'), tags=('iss',))
        lienzo.tag_lower('iss', 'frente')
        self.velocidad = (ancho + 240) / 7.5  # cruza en unos 7-8 segundos (da tiempo a leer el nombre)

    def paso(self, dt):
        dx = self.velocidad * dt
        self.x += dx
        self.lienzo.move('iss', dx, 0)
        if self.x > self.ancho + 130:
            self.lienzo.delete('iss')
            return False
        return True


class PuntoAzulPalido:
    """Easter egg EARTH (en una capa encima de la foto): zoom hacia atrás desde
    la Tierra hasta que solo queda un píxel azul dentro de un rayo de luz,
    como la foto "Un punto azul pálido" de la Voyager 1 (1990)."""
    ZOOM = 6.0

    def __init__(self, c, ancho, alto, texto):
        self.c, self.t = c, 0.0
        self.cx, self.cy = ancho * 0.6, alto * 0.55
        azar = random.Random(1990)
        self.estrellas = []
        for _ in range(int(ancho * alto / 3500)):
            x, y = azar.uniform(0, ancho), azar.uniform(0, alto)
            objeto = c.create_oval(x - 0.8, y - 0.8, x + 0.8, y + 0.8, fill='#0a1a38', outline='')
            self.estrellas.append((objeto, azar.choice(('#ffffff', '#cfdcff', '#ffe9d0')), azar.uniform(0.3, 0.9)))
        self.R = R = max(ancho, alto) * 0.9
        self.k = math.log(R / 1.0) / self.ZOOM
        c.create_oval(self.cx - R, self.cy - R, self.cx + R, self.cy + R, fill='#1f5ea8', outline='', tags=('tierra',))
        for _ in range(8):   # continentes
            a, d = azar.uniform(0, 2 * math.pi), R * azar.uniform(0, 0.65)
            puntos = contorno_irregular(self.cx + d * math.cos(a), self.cy + d * math.sin(a), R * azar.uniform(0.12, 0.3),
                                        azar, 0.45, limite=(self.cx, self.cy, R * 0.98))
            c.create_polygon(*puntos, fill=azar.choice(('#3d7a3a', '#6f7b3c', '#8a7550')), outline='', smooth=True,
                             tags=('tierra',))
        for _ in range(11):  # nubes
            a, d = azar.uniform(0, 2 * math.pi), R * azar.uniform(0, 0.8)
            puntos = contorno_irregular(self.cx + d * math.cos(a), self.cy + d * math.sin(a), R * azar.uniform(0.05, 0.16),
                                        azar, 0.5, limite=(self.cx, self.cy, R * 0.98))
            c.create_polygon(*puntos, fill='#e9f0fb', outline='', smooth=True, tags=('tierra',))
        self.punto = None
        self.rayos = []
        self.ancho, self.alto = ancho, alto
        self.texto = c.create_text(ancho / 2, alto - 30, text=texto, fill='#000000',
                                   font=('Segoe UI Light', 15))

    def paso(self, dt):
        c, self.t = self.c, self.t + dt
        t = self.t
        if self.R > 1.2:
            f = math.exp(-self.k * dt)
            c.scale('tierra', self.cx, self.cy, f, f)
            self.R *= f
        elif self.punto is None:
            c.delete('tierra')
            # rayos de luz del Sol dispersada en la cámara (franjas suaves, en degradado)
            for desfase, ancho_r, color in ((-110, 70, '#6e4f33'), (-25, 34, '#8a6a48'),
                                            (0, 52, '#a7835a'), (80, 40, '#5e4430')):
                x = self.cx + desfase
                for capa in range(5):
                    w = ancho_r * (1 - capa / 5)
                    objeto = c.create_polygon(x - w + 200, -10, x + w + 200, -10,
                                              x + w - 200, self.alto + 10, x - w - 200, self.alto + 10,
                                              fill='#000000', outline='')
                    self.rayos.append((objeto, color, 0.06 + 0.06 * capa))
            self.punto = c.create_oval(self.cx - 1.6, self.cy - 1.6, self.cx + 1.6, self.cy + 1.6,
                                       fill='#b4d8ff', outline='')
            c.tag_raise(self.texto)
        c.configure(background=mezclar('#0a1a38', '#000000', min(t / 4.0, 1.0)))
        for objeto, color, intensidad in self.estrellas:
            c.itemconfig(objeto, fill=mezclar('#000000', color, intensidad * min(max(t - 1.0, 0) / 3.0, 1.0)))
        if self.rayos:
            luz = min((t - self.ZOOM) / 1.8, 1.0)
            for objeto, color, intensidad in self.rayos:
                c.itemconfig(objeto, fill=mezclar('#000000', color, luz * intensidad))
        if t > self.ZOOM + 1.5:
            c.itemconfig(self.texto, fill=mezclar('#000000', '#c9d6f2', min((t - self.ZOOM - 1.5) / 1.5, 1)))
        return t < self.ZOOM + 8.0


class SenalWow:
    """Easter egg 6EQUJ5 (en una capa encima de la foto): monitor verde de
    radiotelescopio antiguo. Pasa un pico de señal, la impresora saca la
    columna 6EQUJ5, se rodea en rojo y aparece un ¡WOW! enorme."""

    def __init__(self, c, ancho, alto, texto_wow, texto_sub):
        self.c, self.ancho, self.alto, self.t = c, ancho, alto, 0.0
        c.configure(background='#010a03')
        verde, verde_oscuro = '#3cff6e', '#0b2a12'
        for y in range(0, alto, 3):
            c.create_line(0, y, ancho, y, fill='#021006')
        self.x0, self.x1, self.y0, self.y1 = 40, ancho * 0.72, alto * 0.18, alto * 0.82
        for i in range(11):
            x = self.x0 + (self.x1 - self.x0) * i / 10
            c.create_line(x, self.y0, x, self.y1, fill=verde_oscuro)
        for i in range(7):
            y = self.y0 + (self.y1 - self.y0) * i / 6
            c.create_line(self.x0, y, self.x1, y, fill=verde_oscuro)
        c.create_text(self.x0, self.y0 - 30, anchor=tk.W, text='BIG EAR  ·  1420,4556 MHz  ·  15-08-1977',
                      fill=verde, font=('Consolas', 13, 'bold'))
        self.valores = [0.0] * 220
        self.traza = c.create_line(*([self.x0, self.y1] * 2), fill=verde, width=2)
        self.azar = random.Random(6)
        # Columna de la impresora (un carácter por línea)
        self.columna = list('1 1 1 2 1 1 1 1 1 3 1 1 6 E Q U J 5 1 1 2 1 1 1 1'.split())
        self.indice_wow = 12
        self.impresos = 0
        self.xi = ancho * 0.84
        c.create_text(self.xi, self.y0 - 30, text='S/N', fill=verde, font=('Consolas', 13, 'bold'))
        self.circulo = None
        self.fondo_wow = c.create_rectangle(0, 0, 0, 0, fill='#010a03', outline='#ff4a4a', width=2)
        self.wow = c.create_text(ancho * 0.42, alto * 0.45, text='', fill='#ff4a4a', font=('Segoe UI Black', 10))
        self.sub = c.create_text(ancho * 0.42, alto * 0.45 + 68, text='', fill='#ff9d9d', font=('Consolas', 13))
        self.texto_wow, self.texto_sub = texto_wow, texto_sub

    def paso(self, dt):
        c, self.t = self.c, self.t + dt
        t = self.t
        # señal: ruido de fondo y, a los 3 s, el paso de la fuente por el haz de la antena
        for _ in range(3):
            pico = math.exp(-((t - 3.0) / 0.55) ** 2)
            self.valores.append(0.08 + 0.05 * self.azar.random() + 0.85 * pico * (0.9 + 0.1 * self.azar.random()))
        self.valores = self.valores[-220:]
        puntos = []
        for i, valor in enumerate(self.valores):
            puntos += [self.x0 + (self.x1 - self.x0) * i / 219, self.y1 - (self.y1 - self.y0) * valor]
        c.coords(self.traza, *puntos)
        # la impresora saca un carácter cada 0,2 s
        while self.impresos < len(self.columna) and t > 0.2 * self.impresos:
            y = self.y0 + 4 + 18 * self.impresos
            c.create_text(self.xi, y, text=self.columna[self.impresos], fill='#3cff6e', font=('Consolas', 14))
            self.impresos += 1
        if self.circulo is None and self.impresos > self.indice_wow + 6:
            ya, yb = self.y0 + 4 + 18 * self.indice_wow - 12, self.y0 + 4 + 18 * (self.indice_wow + 5) + 12
            self.circulo = c.create_oval(self.xi - 22, ya, self.xi + 22, yb, outline='#ff4a4a', width=3)
        if t > 4.6:
            tamano = int(10 + 80 * min((t - 4.6) / 0.35, 1))
            c.itemconfig(self.wow, text=self.texto_wow, font=('Segoe UI Black', tamano))
            c.itemconfig(self.sub, text=self.texto_sub)
            x0, y0, x1, y1 = c.bbox(self.wow)
            sx0, _sy0, sx1, sy1 = c.bbox(self.sub)
            c.coords(self.fondo_wow, min(x0, sx0) - 18, y0 - 6, max(x1, sx1) + 18, sy1 + 12)
        return t < 11.0


class SaludoVulcano:
    """Easter egg SPOCK / VULCAN (en una capa encima de la foto): aparece la
    silueta del saludo vulcano y una frase, y luego se desvanece."""

    def __init__(self, c, ancho, alto, frase):
        self.c, self.t = c, 0.0
        azar = random.Random(7)
        for _ in range(int(ancho * alto / 4000)):
            x, y = azar.uniform(0, ancho), azar.uniform(0, alto)
            c.create_oval(x - 0.7, y - 0.7, x + 0.7, y + 0.7, fill=mezclar(ESPACIO, '#ffffff', azar.uniform(0.2, 0.7)),
                          outline='')
        s = min(ancho, alto) / 420
        cx, cy = ancho / 2, alto / 2 - 30 * s
        grosor = 25 * s
        dedos = (((-40, -36), (-62, -116)), ((-16, -42), (-40, -162)),   # meñique y anular, juntos
                 ((16, -42), (40, -166)), ((40, -36), (62, -148)),       # corazón e índice, juntos
                 ((44, 22), (100, -20)))                                 # pulgar
        # (entre el anular y el corazón queda la "V" del saludo)
        for (x0, y0), (x1, y1) in dedos:
            c.create_line(cx + x0 * s, cy + y0 * s, cx + x1 * s, cy + y1 * s, width=grosor,
                          capstyle=tk.ROUND, fill=ESPACIO, tags=('mano',))
        palma = [(-54, -44), (54, -44), (58, 28), (32, 74), (-32, 74), (-58, 28)]
        c.create_polygon(*[v for x, y in palma for v in (cx + x * s, cy + y * s)], smooth=True,
                         fill=ESPACIO, outline='', tags=('mano',))
        brazo = [(-34, 60), (34, 60), (30, 150), (-30, 150)]
        c.create_polygon(*[v for x, y in brazo for v in (cx + x * s, cy + y * s)], fill=ESPACIO, outline='',
                         tags=('mano',))
        self.frase = c.create_text(cx, cy + 195 * s, text=frase, fill=ESPACIO, font=('Segoe UI Light', 26))

    def paso(self, dt):
        self.t += dt
        t = self.t
        nivel = min(t / 1.0, 1.0) if t < 3.8 else max(0.0, 1 - (t - 3.8) / 1.0)
        self.c.itemconfig('mano', fill=mezclar(ESPACIO, '#e7d6a6', nivel))
        self.c.itemconfig(self.frase, fill=mezclar(ESPACIO, '#cfdcff', max(0.0, min((t - 0.6) / 1.0, 1.0)) * nivel))
        return t < 4.9


class Logo:
    """'Emillium' con un '1.0' pequeño y dorado al lado (como un exponente)
    y un destello que gira. Puede flotar suavemente, como en ingravidez."""

    def __init__(self, lienzo, x, y, tamano, centrado=False, etiqueta='logo'):
        self.lienzo, self.tamano, self.etiqueta = lienzo, tamano, etiqueta
        nombre = lienzo.create_text(x, y, text='Emillium', fill=TEXTO, anchor=tk.CENTER if centrado else tk.W,
                                    font=(FUENTE_TITULO, tamano), tags=(etiqueta, 'frente'))
        if centrado:  # recentrar contando también el "1.0"
            lienzo.move(nombre, -int(tamano * 0.45), 0)
        x0, y0, x1, y1 = lienzo.bbox(nombre)
        self.version = lienzo.create_text(x1 + 3, y0 + (y1 - y0) * 0.22, text=f'v{VERSION}', fill=DORADO,
                                          anchor=tk.W, font=('Segoe UI', max(8, int(tamano * 0.36)), 'bold'),
                                          tags=(etiqueta, 'frente'))
        self.cx, self.cy = x0 - tamano * 0.45, y0 + (y1 - y0) * 0.5
        self.destello = dibujar_destello(lienzo, self.cx, self.cy, tamano * 0.32, DORADO,
                                         etiquetas=('frente',))
        self.caja = (x0, y0, x1, y1)
        self.dx = self.dy = self.dy_version = 0.0

    def animar(self, t, amplitud):
        # Flotar: dos senos lentos con periodos distintos para que no parezca un vaivén mecánico.
        dx = amplitud * 0.45 * math.sin(t * 0.53)
        dy = amplitud * math.sin(t * 1.05)
        self.lienzo.move(self.etiqueta, dx - self.dx, dy - self.dy)
        self.dx, self.dy = dx, dy
        # El "1.0" flota un poquito por su cuenta.
        dy_version = amplitud * 0.25 * math.sin(t * 1.7 + 1.0)
        self.lienzo.move(self.version, 0, dy_version - self.dy_version)
        self.dy_version = dy_version
        # El destello late y gira despacio.
        pulso = 0.5 + 0.5 * math.sin(t * 2.4)
        r = self.tamano * 0.32 * (0.88 + 0.24 * pulso)
        self.lienzo.coords(self.destello, *puntos_destello(self.cx + dx, self.cy + dy, r, t * 0.5))
        self.lienzo.itemconfig(self.destello, fill=mezclar(DORADO, '#fff8e6', pulso * 0.7))


# -----------------------------------------------------------------------------
# Textos en catalán y castellano
# -----------------------------------------------------------------------------
IDIOMAS = (('ca', 'Català'), ('es', 'Castellano'))
IDIOMA_POR_DEFECTO = 'ca'

TEXTOS = {
    'ca': {
        'subtitulo': 'Mesurador de centroides estel·lars',
        'barra_ayuda': 'Clic esquerre: mesurar estrella   ·   Roda: zoom   ·   '
                       'Arrossegar amb el botó dret: moure la foto',
        'paso1': ' 1. Obre una foto ',
        'abrir': 'Obrir foto FITS…',
        'ayuda_abrir': 'Obre una o diverses imatges FITS (.fits, .fit, .fts, també comprimides .fz). Si en tries diverses, '
                       'podràs passar d\'una a l\'altra amb ◀ ▶. També: Ctrl+O.',
        'ninguna_foto': 'Cap foto oberta',
        'info_foto': '{nombre}\nFoto n.º {n} · {ancho} × {alto} píxels',
        'info_hora': 'UTC (mig): {hora}',
        'info_exposicion': ' · {s} s',
        'info_wcs': '✓ Solució astromètrica: AR/Dec disponibles ({escala}″/px)',
        'foto_anterior': '◀ Anterior',
        'foto_siguiente': 'Següent ▶',
        'ayuda_serie': "Passa a la foto anterior o següent de la sèrie (també: tecles RePàg i AvPàg). "
                       'Es manté el gir, el mirall i el zoom.',
        'como_anterior': 'Mesurar com a la foto anterior',
        'paso4': ' 4. Fotometria diferencial (opcional) ',
        'fotometria_texto': "Tria una estrella de comparació de magnitud coneguda (d'un catàleg) i escriu-la: "
                            'la resta tindran la seva magnitud real.',
        'comparacion': 'Comparació:',
        'magnitud_conocida': 'Magnitud:',
        'ayuda_comparacion': "Lletra de l'estrella de comparació. Ha d'estar mesurada a cada foto.",
        'ayuda_magnitud_conocida': "Magnitud de catàleg de l'estrella de comparació, en la mateixa banda que el "
                                   'filtre de la foto (per exemple, 12,35).',
        'fotometria_falta_mag': "Escriu la magnitud coneguda de l'estrella {letra}.",
        'fotometria_ok': 'Magnituds calibrades amb {letra} = {mag}.',
        'fotometria_sin_comp': "Falta mesurar l'estrella {letra} a les fotos: {fotos}.",
        'detalle_mag_cal': 'Magnitud calibrada: {mag} ± {sm} (amb {comp})',
        'col_mag_inst': 'Mag (inst.)',
        'col_mag_cal': 'Mag (cal.)',
        'mag_cal_col': '{l}_Mag cal. ({c}={m})',
        'mpc_boton': 'Informe MPC…',
        'ayuda_mpc': 'Prepara un informe d\'astrometria per al Minor Planet Center (format de 80 columnes) amb '
                     "l'AR/Dec i l'hora de l'objecte a cada foto.",
        'mpc_titulo': 'Informe per al Minor Planet Center',
        'mpc_sin_datos': "Per fer l'informe cal que l'objecte estigui mesurat en fotos amb solució astromètrica "
                         '(AR/Dec) i amb data (DATE-OBS).',
        'mpc_explicacion': "Tria l'objecte (asteroide o cometa) i omple les dades. Cada línia és una mesura en el "
                           'format de 80 columnes del MPC. REVISA-HO abans d\'enviar-ho: comprova a la web del MPC el '
                           'format i les instruccions vigents (ara recomanen el format ADES), i que la solució '
                           'astromètrica de les fotos sigui bona.',
        'mpc_objeto': 'Objecte (lletra):',
        'mpc_designacion': 'Designació provisional (7):',
        'mpc_codigo': "Codi d'observatori:",
        'mpc_tipo': 'Detector (C=CCD, B=CMOS):',
        'mpc_observadores': 'Observadors (OBS):',
        'mpc_medidores': 'Mesuradors (MEA):',
        'mpc_telescopio': 'Telescopi (TEL):',
        'mpc_catalogo': 'Catàleg de referència (NET):',
        'mpc_contacto': 'Contacte (CON, opcional):',
        'mpc_banda': 'Banda de la magnitud:',
        'mpc_descubrimiento': 'Primera observació d\'un objecte nou (*)',
        'mpc_con_magnitud': 'Incloure la magnitud calibrada (fotometria diferencial)',
        'mpc_n_obs': '{n} observacions.',
        'mpc_sin_codigo': 'Sense codi d\'observatori del MPC, posa «XXX» i explica on és al contacte.',
        'mpc_pocas': 'El MPC demana com a mínim 2–3 observacions de la mateixa nit.',
        'mpc_guardado': 'Informe desat a:\n{ruta}',
        'mpc_copiado': 'Copiat al porta-retalls.',
        'mpc_desar': 'Desar .txt…',
        'mpc_copiar': 'Copiar',
        'mpc_cerrar': 'Tancar',
        'mpc_comentario': 'Mesurat amb Emillium {v} (centroide; AR/Dec de la solucio WCS del FITS)',
        'ayuda_como_anterior': "Busca en aquesta foto totes les estrelles mesurades a l'anterior, amb la "
                               "mateixa lletra. Primer calcula quant s'ha desplaçat el camp amb l'estrella més brillant.",
        'como_anterior_ok': "S'han mesurat {n} de {total} estrelles de la foto anterior.",
        'como_anterior_faltan': "No s'han trobat: {letras}. Mesura-les a mà.",
        'detalle_mag': 'Magnitud instrumental: {mag} ± {sm}',
        'col_error': '± (px)',
        'col_ar': 'AR',
        'col_dec': 'Dec',
        'col_jd': 'JD (mig de l\'exposició)',
        'col_exposicion': 'Exposició (s)',
        'col_sigma_flujo': 'σ Flux',
        'col_ar_grados': 'AR (°)',
        'col_dec_grados': 'Dec (°)',
        'col_ar_hms': 'AR (h m s)',
        'col_dec_dms': 'Dec (° ′ ″)',
        'ar_col': '{l}_AR (°)',
        'dec_col': '{l}_Dec (°)',
        'paso2': ' 2. Modifica la foto (opcional) ',
        'giro_label': 'Gir de la foto (graus, amb signe):',
        'aplicar_giro': 'Aplicar el gir',
        'ayuda_giro': "Escriu l'angle amb el seu signe (per exemple +12,5 o -7) i prem Retorn o "
                      "«Aplicar el gir». L'angle és sempre respecte a la foto original. "
                      'Els botons −10° i +10° el canvien de 10 en 10.',
        'sentido_pos': "+ : sentit antihorari (cap a l'esquerra)",
        'sentido_neg': '− : sentit horari (cap a la dreta)',
        'sentido_cero': 'Sense girar',
        'angulo_invalido': 'Angle no vàlid',
        'espejo': 'Mirall:',
        'horizontal': 'Horitzontal',
        'vertical': 'Vertical',
        'ayuda_espejo_h': "Inverteix la foto d'esquerra a dreta, com un mirall.",
        'ayuda_espejo_v': 'Inverteix la foto de dalt a baix.',
        'volver': 'Tornar a la foto original',
        'ahora_esta': 'Ara la foto està:',
        'original': 'Original',
        'o_espejo_h': 'mirall horitzontal',
        'o_espejo_v': 'mirall vertical',
        'o_giro': 'girada {g}° ({sentido})',
        'antihorario': 'antihorari',
        'horario': 'horari',
        'paso3': ' 3. Mesura estrelles ',
        'medir_texto': 'Tria quina estrella mesures i fes-hi clic. A cada foto, '
                       'fes servir la mateixa lletra per a la mateixa estrella.',
        'estrella': 'Estrella:',
        'avanzar': 'Passar a la lletra següent',
        'ayuda_estrella': "La lletra identifica l'estrella a totes les fotos. Al full de càlcul "
                          'cada estrella té les seves columnes (A_X, A_Y, B_X…).',
        'ultima': 'Última estrella: -',
        'titulo_medida': 'Estrella {letra} · foto {n}:',
        'detalle': "Flux: {flujo}\nMida (FWHM): {fwhm} px\nPíxel més brillant: {ps} vegades el soroll d'un píxel",
        'sin_estrella': "Aquí no hi ha cap estrella clara. Fes clic més a prop d'una estrella "
                        '(o baixa el «Llindar sobre el soroll» a les opcions avançades).',
        'borde': "Atenció: l'estrella és a la vora o en una zona sense dades; "
                 'el centre pot ser menys precís.',
        'sustituida': "S'ha substituït la mesura anterior de l'estrella {letra} en aquesta foto.",
        'duplicada_t': 'Estrella ja mesurada',
        'duplicada': "Aquest centre és el de l'estrella {otra}, que ja està mesurada en aquesta foto "
                     '(a {d} px).\n\nNo s\'ha desat com a estrella {letra}. Si volies mesurar una altra '
                     "estrella, fes clic més a prop d'ella o fes més petita la «Distància màx. al clic» "
                     'a les opcions avançades.',
        'otra_foto': "Aquesta mesura és d'una altra foto o orientació.",
        'mostrar_avanzadas': 'Mostrar les opcions avançades',
        'ocultar_avanzadas': 'Amagar les opcions avançades',
        'adv_busqueda': 'Distància màx. al clic (px)',
        'ayuda_busqueda': "Fins a quina distància del clic es busca l'estrella. "
                          'Si hi ha estrelles molt juntes, abaixa-la.',
        'adv_apertura': 'Radi per al centre (px)',
        'ayuda_apertura': 'Radi del cercle de píxels que es fa servir per calcular el centre. '
                          "L'ideal és unes 2 o 3 vegades la mida (FWHM) de les estrelles.",
        'recalculada': "S'ha tornat a calcular l'estrella {letra} amb les opcions noves.",
        'adv_umbral': 'Llindar sobre el soroll (σ)',
        'ayuda_umbral': 'Quantes vegades el soroll del cel ha de superar un píxel perquè compti per '
                        "calcular el centre (i perquè es detecti una estrella). Recomanat: 3.\n"
                        "Més baix: el soroll del cel desvia el centre de les estrelles febles.\n"
                        "Més alt: només compten els píxels més brillants (amb 5 o més, en estrelles "
                        'febles gairebé només el píxel central).',
        'zoom': 'Zoom:',
        'ver_entera': 'Veure la foto sencera',
        'zoom_pct': 'Zoom {p} %',
        'brillo': 'Brillantor:',
        'modos_brillo': ('Suau', 'Normal', 'Ressaltar estrelles febles'),
        'ayuda_brillo': 'Només canvia com es veu la foto a la pantalla. No afecta les mesures.',
        'marcas': 'Marcar les estrelles mesurades',
        'idioma': 'Idioma:',
        'mediciones': 'Mesures',
        'ninguna_todavia': '(encara cap)',
        'contador_1': '(1 mesura · doble clic en una fila per veure-la)',
        'contador_n': '({n} mesures · doble clic en una fila per veure-la)',
        'borrar_medida': 'Esborrar la mesura',
        'borrar_medidas': 'Esborrar les {n} mesures',
        'ayuda_borrar_medida': 'Esborra les mesures seleccionades a la taula. També: tecla Supr.',
        'borrar_todas': 'Netejar mesures',
        'ayuda_borrar_todas': 'Esborra TOTES les mesures i torna a numerar les fotos des de l\'1. '
                              'Abans demana confirmació.',
        'borrar_msg_n': "Segur que vols esborrar les {n} mesures? Això no es pot desfer.\n\n"
                        "(Si encara no les has desat, prem «No» i desa-les abans per a Excel.)",
        'cerrar_foto': 'Tancar la foto',
        'ayuda_cerrar_foto': "Treu la foto i torna a la pantalla d'inici. Les mesures no s'esborren.",
        'borrar_ultima': "Esborrar l'última",
        'guardar': 'Desar per a Excel…',
        'ayuda_guardar': 'Desa les mesures en un full de càlcul (.xlsx): una fila per foto i '
                         'columnes X, Y i flux per a cada estrella. També: Ctrl+S.',
        'col_foto_n': 'N.º foto',
        'col_foto': 'Foto',
        'col_estrella': 'Estrella',
        'col_orientacion': 'Orientació',
        'col_flujo': 'Flux',
        'cursor': 'Ratolí: X = {x}   Y = {y}   brillantor = {b}',
        'sin_dato': 'sense dades',
        'fuera': 'Ratolí: fora de la foto',
        'lema': 'Mesura la posició de les estrelles a les teves fotos FITS',
        'pasos_bienvenida': ('Obre una foto FITS  (botó de la dreta o Ctrl+O)',
                             'Si cal, gira-la o posa-la en mirall',
                             'Tria la lletra i fes clic a les estrelles',
                             'Desa les mesures per a Excel'),
        'abrir_titulo': 'Obrir foto FITS',
        'tipo_fits': 'Fotos FITS',
        'tipo_todos': 'Tots els fitxers',
        'error_abrir': "No s'ha pogut obrir la foto",
        'foto_vacia_t': 'Foto buida',
        'foto_vacia': 'La foto no conté cap valor vàlid.',
        'angulo_t': 'Angle',
        'angulo_msg': 'Escriu un nombre de graus, per exemple 12,5.',
        'borrar_t': 'Esborrar',
        'borrar_msg': 'Segur que vols esborrar totes les mesures?',
        'guardar_t': 'Desar',
        'sin_medidas': 'Encara no hi ha mesures. Fes clic en alguna estrella primer.',
        'guardar_titulo': 'Desar les mesures',
        'archivo_defecto': 'mesures',
        'tipo_xlsx': "Llibre d'Excel (.xlsx)",
        'tipo_csv': 'CSV per a Excel en català o castellà (.csv)',
        'error_guardar': "No s'ha pogut desar",
        'abierto_otro': "El fitxer està obert en un altre programa (Excel?). Tanca'l i torna-ho a provar.",
        'guardado_t': 'Desat',
        'guardado': "S'han desat {n} mesures de {f} fotos a:\n{ruta}",
        'hoja_estrellas': 'Coordenades píxels',
        'hoja_todas': 'Totes les mesures',
        'x_col': '{l}_X (px)',
        'y_col': '{l}_Y (px)',
        'flujo_col': '{l}_Flux',
        'col_foto_n_largo': 'N.º de foto',
        'col_hora': "Data i hora UTC (mig de l'exposició)",
        'm42': "M42 · Nebulosa d'Orió",
        'crab': 'M1 · Nebulosa del Cranc',
        'm16': 'M16 · Pilars de la Creació',
        'voyager': 'Voyager · 1977',
        'iss': 'ISS · Estació Espacial Internacional',
        'odisea': "2001: una odissea de l'espai",
        'pale_blue_dot': 'Un punt blau pàl·lid · Voyager 1, 14 de febrer de 1990',
        'wow': 'WOW!',
        'wow_sub': 'Senyal Wow! · 6EQUJ5 · 72 segons que mai no es van repetir',
        'vulcano': 'Llarga vida i prosperitat',
        't42_titulo': 'Pensament Profund',
        't42_calc': "Calculant la resposta definitiva a la vida, l'univers i tot plegat",
        't42_res': 'Resultat: 42',
        'aceptar': "D'acord",
        'err_no_fits': 'No és un fitxer FITS vàlid (no comença per SIMPLE).',
        'err_cabecera': 'El fitxer està incomplet: la capçalera FITS es talla.',
        'err_datos': 'El fitxer està incomplet: falten dades de la imatge.',
        'err_bitpix': 'BITPIX = {bitpix} no és un tipus FITS vàlid.',
        'err_comprimido': 'La imatge FITS està comprimida (.fz). Descomprimeix-la primer '
                          '(per exemple amb el programa «funpack»).',
        'err_sin_imagen': 'El fitxer FITS no conté cap imatge de 2 dimensions.',
        'err_no_2d': 'Les dades del FITS no formen una imatge 2D.',
        'err_pequena': 'La imatge és massa petita per girar-la.',
        'err_cabecera_rara': 'La capçalera FITS té valors que no tenen sentit (mides negatives o no vàlides).',
        'err_fz_tipo': 'Aquest FITS està comprimit amb «{tipo}», que Emillium encara no sap descomprimir '
                       '(sí: RICE_1, GZIP_1 i GZIP_2). Descomprimeix-lo abans amb el programa «funpack».',
        'err_enorme': 'La imatge és massa gran (més de 400 megapíxels) per obrir-la.',
        'nota_color': 'Imatge en color (3 plans): es fa servir la mitjana dels tres.',
        'nota_planos': 'La imatge té diversos plans: es fa servir el primer.',
    },
    'es': {
        'subtitulo': 'Medidor de centroides estelares',
        'barra_ayuda': 'Clic izquierdo: medir estrella   ·   Rueda: zoom   ·   '
                       'Arrastrar con botón derecho: mover la foto',
        'paso1': ' 1. Abre una foto ',
        'abrir': 'Abrir foto FITS…',
        'ayuda_abrir': 'Abre una o varias imágenes FITS (.fits, .fit, .fts, también comprimidas .fz). Si eliges varias, '
                       'podrás pasar de una a otra con ◀ ▶. También: Ctrl+O.',
        'ninguna_foto': 'Ninguna foto abierta',
        'info_foto': '{nombre}\nFoto n.º {n} · {ancho} × {alto} píxeles',
        'info_hora': 'UTC (centro): {hora}',
        'info_exposicion': ' · {s} s',
        'info_wcs': '✓ Solución astrométrica: AR/Dec disponibles ({escala}″/px)',
        'foto_anterior': '◀ Anterior',
        'foto_siguiente': 'Siguiente ▶',
        'ayuda_serie': 'Pasa a la foto anterior o siguiente de la serie (también: teclas RePág y AvPág). '
                       'Se mantienen el giro, el espejo y el zoom.',
        'como_anterior': 'Medir como en la foto anterior',
        'paso4': ' 4. Fotometría diferencial (opcional) ',
        'fotometria_texto': 'Elige una estrella de comparación de magnitud conocida (de un catálogo) y escríbela: '
                            'las demás tendrán su magnitud real.',
        'comparacion': 'Comparación:',
        'magnitud_conocida': 'Magnitud:',
        'ayuda_comparacion': 'Letra de la estrella de comparación. Tiene que estar medida en cada foto.',
        'ayuda_magnitud_conocida': 'Magnitud de catálogo de la estrella de comparación, en la misma banda que el '
                                   'filtro de la foto (por ejemplo, 12,35).',
        'fotometria_falta_mag': 'Escribe la magnitud conocida de la estrella {letra}.',
        'fotometria_ok': 'Magnitudes calibradas con {letra} = {mag}.',
        'fotometria_sin_comp': 'Falta medir la estrella {letra} en las fotos: {fotos}.',
        'detalle_mag_cal': 'Magnitud calibrada: {mag} ± {sm} (con {comp})',
        'col_mag_inst': 'Mag (inst.)',
        'col_mag_cal': 'Mag (cal.)',
        'mag_cal_col': '{l}_Mag cal. ({c}={m})',
        'mpc_boton': 'Informe MPC…',
        'ayuda_mpc': 'Prepara un informe de astrometría para el Minor Planet Center (formato de 80 columnas) con '
                     'la AR/Dec y la hora del objeto en cada foto.',
        'mpc_titulo': 'Informe para el Minor Planet Center',
        'mpc_sin_datos': 'Para hacer el informe el objeto tiene que estar medido en fotos con solución astrométrica '
                         '(AR/Dec) y con fecha (DATE-OBS).',
        'mpc_explicacion': 'Elige el objeto (asteroide o cometa) y rellena los datos. Cada línea es una medida en el '
                           'formato de 80 columnas del MPC. REVÍSALO antes de enviarlo: comprueba en la web del MPC el '
                           'formato y las instrucciones vigentes (ahora recomiendan el formato ADES), y que la solución '
                           'astrométrica de las fotos sea buena.',
        'mpc_objeto': 'Objeto (letra):',
        'mpc_designacion': 'Designación provisional (7):',
        'mpc_codigo': 'Código de observatorio:',
        'mpc_tipo': 'Detector (C=CCD, B=CMOS):',
        'mpc_observadores': 'Observadores (OBS):',
        'mpc_medidores': 'Medidores (MEA):',
        'mpc_telescopio': 'Telescopio (TEL):',
        'mpc_catalogo': 'Catálogo de referencia (NET):',
        'mpc_contacto': 'Contacto (CON, opcional):',
        'mpc_banda': 'Banda de la magnitud:',
        'mpc_descubrimiento': 'Primera observación de un objeto nuevo (*)',
        'mpc_con_magnitud': 'Incluir la magnitud calibrada (fotometría diferencial)',
        'mpc_n_obs': '{n} observaciones.',
        'mpc_sin_codigo': 'Sin código de observatorio del MPC, pon «XXX» y explica dónde está en el contacto.',
        'mpc_pocas': 'El MPC pide como mínimo 2–3 observaciones de la misma noche.',
        'mpc_guardado': 'Informe guardado en:\n{ruta}',
        'mpc_copiado': 'Copiado al portapapeles.',
        'mpc_desar': 'Guardar .txt…',
        'mpc_copiar': 'Copiar',
        'mpc_cerrar': 'Cerrar',
        'mpc_comentario': 'Medido con Emillium {v} (centroide; AR/Dec de la solucion WCS del FITS)',
        'ayuda_como_anterior': 'Busca en esta foto todas las estrellas medidas en la anterior, con la misma '
                               'letra. Primero calcula cuánto se ha desplazado el campo con la estrella más brillante.',
        'como_anterior_ok': 'Se han medido {n} de {total} estrellas de la foto anterior.',
        'como_anterior_faltan': 'No se han encontrado: {letras}. Mídelas a mano.',
        'detalle_mag': 'Magnitud instrumental: {mag} ± {sm}',
        'col_error': '± (px)',
        'col_ar': 'AR',
        'col_dec': 'Dec',
        'col_jd': 'JD (centro de la exposición)',
        'col_exposicion': 'Exposición (s)',
        'col_sigma_flujo': 'σ Flujo',
        'col_ar_grados': 'AR (°)',
        'col_dec_grados': 'Dec (°)',
        'col_ar_hms': 'AR (h m s)',
        'col_dec_dms': 'Dec (° ′ ″)',
        'ar_col': '{l}_AR (°)',
        'dec_col': '{l}_Dec (°)',
        'paso2': ' 2. Modifica la foto (opcional) ',
        'giro_label': 'Giro de la foto (grados, con signo):',
        'aplicar_giro': 'Aplicar el giro',
        'ayuda_giro': 'Escribe el ángulo con su signo (por ejemplo +12,5 o -7) y pulsa Intro o '
                      '«Aplicar el giro». El ángulo es siempre respecto a la foto original. '
                      'Los botones −10° y +10° lo cambian de 10 en 10.',
        'sentido_pos': '+ : sentido antihorario (hacia la izquierda)',
        'sentido_neg': '− : sentido horario (hacia la derecha)',
        'sentido_cero': 'Sin girar',
        'angulo_invalido': 'Ángulo no válido',
        'espejo': 'Espejo:',
        'horizontal': 'Horizontal',
        'vertical': 'Vertical',
        'ayuda_espejo_h': 'Invierte la foto de izquierda a derecha, como un espejo.',
        'ayuda_espejo_v': 'Invierte la foto de arriba a abajo.',
        'volver': 'Volver a la foto original',
        'ahora_esta': 'Ahora la foto está:',
        'original': 'Original',
        'o_espejo_h': 'espejo horizontal',
        'o_espejo_v': 'espejo vertical',
        'o_giro': 'girada {g}° ({sentido})',
        'antihorario': 'antihorario',
        'horario': 'horario',
        'paso3': ' 3. Mide estrellas ',
        'medir_texto': 'Elige qué estrella mides y haz clic en ella. En cada foto, '
                       'usa la misma letra para la misma estrella.',
        'estrella': 'Estrella:',
        'avanzar': 'Pasar a la letra siguiente',
        'ayuda_estrella': 'La letra identifica la estrella en todas las fotos. En la hoja de cálculo '
                          'cada estrella tiene sus columnas (A_X, A_Y, B_X…).',
        'ultima': 'Última estrella: -',
        'titulo_medida': 'Estrella {letra} · foto {n}:',
        'detalle': 'Flujo: {flujo}\nTamaño (FWHM): {fwhm} px\nPíxel más brillante: {ps} veces el ruido de un píxel',
        'sin_estrella': 'No hay ninguna estrella clara ahí. Haz clic más cerca de una estrella '
                        '(o baja el «Umbral sobre el ruido» en las opciones avanzadas).',
        'borde': 'Atención: la estrella está en el borde o en una zona sin datos; '
                 'el centro puede ser menos preciso.',
        'sustituida': 'Se ha sustituido la medida anterior de la estrella {letra} en esta foto.',
        'duplicada_t': 'Estrella ya medida',
        'duplicada': 'Este centro es el de la estrella {otra}, que ya está medida en esta foto '
                     '(a {d} px).\n\nNo se ha guardado como estrella {letra}. Si querías medir otra '
                     'estrella, haz clic más cerca de ella o reduce la «Distancia máx. al clic» '
                     'en las opciones avanzadas.',
        'otra_foto': 'Esta medida es de otra foto u orientación.',
        'mostrar_avanzadas': 'Mostrar opciones avanzadas',
        'ocultar_avanzadas': 'Ocultar opciones avanzadas',
        'adv_busqueda': 'Distancia máx. al clic (px)',
        'ayuda_busqueda': 'Hasta qué distancia del clic se busca la estrella. '
                          'Si hay estrellas muy juntas, bájalo.',
        'adv_apertura': 'Radio para el centro (px)',
        'ayuda_apertura': 'Radio del círculo de píxeles que se usa para calcular el centro. '
                          'Lo ideal es unas 2 o 3 veces el tamaño (FWHM) de las estrellas.',
        'recalculada': 'Se ha vuelto a calcular la estrella {letra} con las opciones nuevas.',
        'adv_umbral': 'Umbral sobre el ruido (σ)',
        'ayuda_umbral': 'Cuántas veces el ruido del cielo tiene que superar un píxel para contar al '
                        'calcular el centro (y para detectar una estrella). Recomendado: 3.\n'
                        'Más bajo: el ruido del cielo desvía el centro de las estrellas débiles.\n'
                        'Más alto: solo cuentan los píxeles más brillantes (con 5 o más, en estrellas '
                        'débiles casi solo el píxel central).',
        'zoom': 'Zoom:',
        'ver_entera': 'Ver foto entera',
        'zoom_pct': 'Zoom {p} %',
        'brillo': 'Brillo:',
        'modos_brillo': ('Suave', 'Normal', 'Resaltar estrellas débiles'),
        'ayuda_brillo': 'Solo cambia cómo se ve la foto en pantalla. No afecta a las medidas.',
        'marcas': 'Marcar las estrellas medidas',
        'idioma': 'Idioma:',
        'mediciones': 'Mediciones',
        'ninguna_todavia': '(ninguna todavía)',
        'contador_1': '(1 medida · doble clic en una fila para verla)',
        'contador_n': '({n} medidas · doble clic en una fila para verla)',
        'borrar_medida': 'Borrar la medida',
        'borrar_medidas': 'Borrar las {n} medidas',
        'ayuda_borrar_medida': 'Borra las medidas seleccionadas en la tabla. También: tecla Supr.',
        'borrar_todas': 'Limpiar medidas',
        'ayuda_borrar_todas': 'Borra TODAS las medidas y vuelve a numerar las fotos desde el 1. '
                              'Antes pide confirmación.',
        'borrar_msg_n': '¿Seguro que quieres borrar las {n} medidas? No se puede deshacer.\n\n'
                        '(Si todavía no las has guardado, pulsa «No» y guárdalas antes para Excel.)',
        'cerrar_foto': 'Cerrar la foto',
        'ayuda_cerrar_foto': 'Quita la foto y vuelve a la pantalla de inicio. Las medidas no se borran.',
        'borrar_ultima': 'Borrar la última',
        'guardar': 'Guardar para Excel…',
        'ayuda_guardar': 'Guarda las medidas en una hoja de cálculo (.xlsx): una fila por foto y '
                         'columnas X, Y y flujo para cada estrella. También: Ctrl+S.',
        'col_foto_n': 'N.º foto',
        'col_foto': 'Foto',
        'col_estrella': 'Estrella',
        'col_orientacion': 'Orientación',
        'col_flujo': 'Flujo',
        'cursor': 'Ratón: X = {x}   Y = {y}   brillo = {b}',
        'sin_dato': 'sin dato',
        'fuera': 'Ratón: fuera de la foto',
        'lema': 'Mide la posición de las estrellas en tus fotos FITS',
        'pasos_bienvenida': ('Abre una foto FITS  (botón de la derecha o Ctrl+O)',
                             'Si hace falta, gírala o ponla en espejo',
                             'Elige la letra y haz clic en las estrellas',
                             'Guarda las mediciones para Excel'),
        'abrir_titulo': 'Abrir foto FITS',
        'tipo_fits': 'Fotos FITS',
        'tipo_todos': 'Todos los archivos',
        'error_abrir': 'No se ha podido abrir la foto',
        'foto_vacia_t': 'Foto vacía',
        'foto_vacia': 'La foto no contiene ningún valor válido.',
        'angulo_t': 'Ángulo',
        'angulo_msg': 'Escribe un número de grados, por ejemplo 12,5.',
        'borrar_t': 'Borrar',
        'borrar_msg': '¿Seguro que quieres borrar todas las mediciones?',
        'guardar_t': 'Guardar',
        'sin_medidas': 'Todavía no hay mediciones. Haz clic en alguna estrella primero.',
        'guardar_titulo': 'Guardar mediciones',
        'archivo_defecto': 'mediciones',
        'tipo_xlsx': 'Libro de Excel (.xlsx)',
        'tipo_csv': 'CSV para Excel en español o catalán (.csv)',
        'error_guardar': 'No se ha podido guardar',
        'abierto_otro': 'El archivo está abierto en otro programa (¿Excel?). Ciérralo y vuelve a probar.',
        'guardado_t': 'Guardado',
        'guardado': '{n} medidas de {f} fotos guardadas en:\n{ruta}',
        'hoja_estrellas': 'Coordenadas píxeles',
        'hoja_todas': 'Todas las medidas',
        'x_col': '{l}_X (px)',
        'y_col': '{l}_Y (px)',
        'flujo_col': '{l}_Flujo',
        'col_foto_n_largo': 'N.º de foto',
        'col_hora': 'Fecha y hora UTC (centro de la exposición)',
        'm42': 'M42 · Nebulosa de Orión',
        'crab': 'M1 · Nebulosa del Cangrejo',
        'm16': 'M16 · Pilares de la Creación',
        'voyager': 'Voyager · 1977',
        'iss': 'ISS · Estación Espacial Internacional',
        'odisea': '2001: una odisea del espacio',
        'pale_blue_dot': 'Un punto azul pálido · Voyager 1, 14 de febrero de 1990',
        'wow': '¡WOW!',
        'wow_sub': 'Señal Wow! · 6EQUJ5 · 72 segundos que nunca se repitieron',
        'vulcano': 'Larga vida y prosperidad',
        't42_titulo': 'Pensamiento Profundo',
        't42_calc': 'Calculando la respuesta definitiva a la vida, el universo y todo lo demás',
        't42_res': 'Resultado: 42',
        'aceptar': 'Aceptar',
        'err_no_fits': 'No es un archivo FITS válido (no empieza por SIMPLE).',
        'err_cabecera': 'El archivo está incompleto: la cabecera FITS se corta.',
        'err_datos': 'El archivo está incompleto: faltan datos de la imagen.',
        'err_bitpix': 'BITPIX = {bitpix} no es un tipo FITS válido.',
        'err_comprimido': 'La imagen FITS está comprimida (.fz). Descomprímela primero '
                          '(por ejemplo con el programa «funpack»).',
        'err_sin_imagen': 'El archivo FITS no contiene ninguna imagen de 2 dimensiones.',
        'err_no_2d': 'Los datos del FITS no forman una imagen 2D.',
        'err_pequena': 'La imagen es demasiado pequeña para girarla.',
        'err_cabecera_rara': 'La cabecera FITS tiene valores sin sentido (tamaños negativos o no válidos).',
        'err_fz_tipo': 'Este FITS está comprimido con «{tipo}», que Emillium todavía no sabe descomprimir '
                       '(sí: RICE_1, GZIP_1 y GZIP_2). Descomprímelo antes con el programa «funpack».',
        'err_enorme': 'La imagen es demasiado grande (más de 400 megapíxeles) para abrirla.',
        'nota_color': 'Imagen en color (3 planos): se usa la media de los tres.',
        'nota_planos': 'La imagen tiene varios planos: se usa el primero.',
    },
}

# Textos nuevos de la versión 1.5 (y algunos de antes que cambian). Se añaden a TEXTOS.
TEXTOS_15 = {
    'ca': {
        # --- fotometria
        'paso4': ' 4. Fotometria diferencial (opcional) ',
        'fotometria_texto': "Escriu les estrelles de comparació. Amb la magnitud de catàleg (B=12,35) s'obtenen "
                            'magnituds reals; només amb la lletra (B C D), magnituds relatives per a corbes de llum.',
        'comparaciones': 'Comparació:',
        'ayuda_comparaciones': "Lletres de les estrelles de comparació separades per espais, amb la seva magnitud de "
                               'catàleg i, si vols, el seu índex de color: «B=12,35 C=11,90/0,65 D=13,1». Si alguna no '
                               "té magnitud, tot es fa en mode relatiu. Han d'estar mesurades a cada foto.",
        'boton_fotometria': 'Fotometria i corbes de llum…',
        'ayuda_boton_fotometria': 'Corbes de llum, estrelles variables, comparacions automàtiques i magnituds dels '
                                  'catàlegs guardats.',
        'fot_formato': 'No entenc les comparacions. Exemple: B=12,35 C=11,90/0,65 D',
        'fot_relativa': 'Fotometria RELATIVA amb {letras} (falta la magnitud de catàleg d\'alguna).',
        'fot_calibrada': 'Magnituds calibrades amb {n} comparació/ons: {letras}.',
        'fot_sin_comp': 'Encara no hi ha estrelles de comparació: es mostren magnituds instrumentals.',
        'col_mag_rel': 'Mag (rel.)',
        'mag_cal_col': '{l}_Mag cal.',
        'mag_rel_col': '{l}_Mag rel.',
        'sigma_cal_col': '{l}_σMag dif.',
        'detalle_mag_cal': 'Magnitud calibrada: {mag} ± {sm} (amb {comp})',
        'detalle_mag_rel': 'Magnitud relativa: {mag} ± {sm} (amb {comp})',
        'detalle_snr': "Senyal/soroll: {snr} (de tota l'estrella, obertura de {r} px)",
        'fot_titulo': 'Fotometria i corbes de llum',
        'fot_explicacion': "Magnitud de cada estrella al llarg de la sèrie. «χ²» compara la dispersió amb la que "
                           'esperaríem només pels errors: molt més gran que 1 vol dir que l\'estrella canvia de brillantor '
                           "(o que alguna mesura és dolenta). Selecciona files per veure'n la corba de llum.",
        'fot_coef_color': 'Coef. de color T:',
        'ayuda_coef_color': 'Coeficient de color de la teva càmera i filtre (normalment entre −0,2 i 0,2). Amb 0 no es '
                            'fa correcció de color. Cal posar també els índexs de color de les estrelles.',
        'fot_colores': 'Índexs de color:',
        'ayuda_colores': "Índex de color (per exemple B−V) de les estrelles objectiu: «A=0,65 D=1,10». S'utilitza amb "
                         'el coeficient T.',
        'fot_n': 'Fotos',
        'fot_media': 'Mag mitjana',
        'fot_desv': 'Desviació',
        'fot_amplitud': 'Amplitud',
        'fot_chi2': 'χ² red.',
        'fot_diagnostico': 'Diagnòstic',
        'fot_es_comp': 'Comparació (es calcula amb les altres)',
        'fot_saturada': 'Saturada: la magnitud no és fiable',
        'fot_pocas': 'Massa poques fotos (mínim 3)',
        'fot_variable': '⚠ Sembla VARIABLE',
        'fot_constante': 'Constant (dins dels errors)',
        'fot_actualizar': 'Actualitzar',
        'fot_auto': 'Triar comparacions automàticament',
        'fot_auto_pocas': 'Calen almenys 3 estrelles mesurades a totes les fotos (sense saturar) per triar comparacions.',
        'fot_banda': 'Banda:',
        'fot_catalogo': 'Magnituds del catàleg',
        'fot_cat_sin_comps': 'Escriu primer les lletres de les estrelles de comparació.',
        'fot_cat_sin_cache': "No hi ha {cat} guardat per a aquest camp. Descarrega'l a «Catàlegs i astrometria».",
        'fot_cat_ok': "S'han posat les magnituds de {n} estrelles amb {cat}.",
        'fot_cat_faltan': 'Sense parella al catàleg (o sense AR/Dec): {letras}.',
        'fot_exportar': 'Exportar corba (CSV)…',
        'fot_archivo': 'corba_de_llum.csv',
        'eje_horas': 'Hores des de {fecha} UTC',
        'eje_foto': 'Número de foto',
        # --- opcions avançades
        'adv_r_fot': 'Radi fotomètric (px)',
        'adv_r_int': 'Anell del cel: interior (px)',
        'adv_r_ext': 'Anell del cel: exterior (px)',
        'ayuda_r_fot': "Radi del cercle on se suma la llum de l'estrella (fotometria d'obertura). Fes servir el mateix a "
                       'totes les estrelles. Un bon valor: 1,5–2 vegades la FWHM.',
        'ayuda_r_anillo': "L'anell on es mesura el fons del cel, al voltant de l'estrella. Ha de quedar fora de "
                          "l'estrella i, si pot ser, lluny d'altres estrelles.",
        # --- qualitat
        'col_calidad': 'Qualitat',
        'calidad_col': '{l}_Qualitat',
        'calidad_bien': 'Qualitat: ✓ cap problema detectat',
        'calidad_titulo': 'Qualitat:',
        'calidad_SAT': "Saturada: l'estrella arriba al límit de la càmera; el flux i el centre no són fiables.",
        'calidad_CON': "Contaminada: una estrella veïna posa més d'un 1 % de llum dins de l'obertura.",
        'calidad_BOR': 'Toca la vora de la foto o una zona sense dades.',
        'calidad_SNR': 'Senyal/soroll baix (< 10): magnitud amb més de 0,1 mag d\'error.',
        'calidad_FWHM': "Mida estranya comparada amb les altres estrelles (galàxia, raig còsmic o dues estrelles juntes?).",
        'calidad_ELO': 'Estrella allargada (moguda o doble sense separar).',
        'calidad_corto_SAT': 'saturada',
        'calidad_corto_CON': 'contaminada',
        'calidad_corto_BOR': 'a la vora',
        'calidad_corto_SNR': 'S/R baix',
        'calidad_corto_FWHM': 'mida estranya',
        'calidad_corto_ELO': 'allargada',
        'col_snr': 'S/R',
        'col_r_fot': 'Radi fot. (px)',
        'col_fondo': 'Fons del cel',
        # --- panell 5
        'paso5': ' 5. Anàlisi ',
        'boton_detectar': 'Detectar estrelles',
        'ayuda_detectar': 'Busca totes les estrelles de la foto i les marca (cian; vermell = saturada). Torna a prémer per '
                          "treure les marques. També serveix per comprovar l'astrometria.",
        'detectadas': "S'han detectat {n} estrelles (FWHM típica {fwhm} px).",
        'boton_dobles': 'Estrelles dobles…',
        'ayuda_dobles': "Separació i angle de posició entre dues estrelles mesurades, a cada foto i de mitjana.",
        'boton_moviles': 'Objectes mòbils…',
        'ayuda_moviles': "Quines estrelles mesurades es mouen respecte a les altres (asteroides, cometes), a quina "
                         'velocitat i cap on.',
        'boton_serie': 'Revisar la sèrie…',
        'ayuda_serie_rev': 'Revisa totes les fotos de la sèrie: fons, estrelles, nitidesa i desplaçament, i marca les '
                           'dolentes (desenfocades, mogudes, amb núvols...).',
        'boton_catalogos': 'Catàlegs i astrometria…',
        'ayuda_catalogos': "Gaia DR3, APASS, VSX i asteroides coneguts (SkyBoT). Comprova i millora la solució "
                           'astromètrica i identifica les estrelles mesurades.',
        'boton_wcs_manual': 'WCS manual…',
        'ayuda_wcs_manual': "Calcula la solució astromètrica (AR/Dec) a partir de 3 o més estrelles mesurades de les "
                            'quals coneixes les coordenades.',
        'boton_guardar_sesion': 'Desar sessió…',
        'boton_abrir_sesion': 'Obrir sessió…',
        'ayuda_sesion': 'Desa (o recupera) tota la feina: fotos, mesures, comparacions, opcions i solucions '
                        'astromètriques, per continuar un altre dia.',
        'necesita_foto': 'Obre primer una foto.',
        # --- dobles
        'dobles_titulo': 'Estrelles dobles',
        'dobles_explicacion': "Separació i angle de posició (des del nord cap a l'est) de la secundària vista des de la "
                              'principal. En segons d\'arc cal que les fotos tinguin solució astromètrica. L\'error és '
                              "l'error de la mitjana (desviació / √n).",
        'dobles_principal': 'Principal:',
        'dobles_secundaria': 'Secundària:',
        'dobles_calcular': 'Calcular',
        'dobles_sep_px': 'Separació (px)',
        'dobles_sep_arc': 'Separació (″)',
        'dobles_pa': 'Angle (°)',
        'dobles_n': '{n} fotos amb les dues estrelles.',
        'dobles_media_px': 'Separació mitjana: {m} ± {e} px',
        'dobles_media_arc': 'ρ = {m} ± {e}″    θ = {pa} ± {epa}°',
        'dobles_sin_wcs': 'Sense solució astromètrica: només separació en píxels.',
        'dobles_pocas': 'Cal haver mesurat almenys dues estrelles.',
        # --- mòbils
        'moviles_titulo': 'Objectes mòbils',
        'moviles_explicacion': "Per a cada lletra mesurada en 2 fotos o més: quant s'ha mogut respecte a les estrelles. "
                               'Amb solució astromètrica, en segons d\'arc; si no, en píxels respecte a l\'estrella més '
                               "brillant. «rms» és la distància típica a un moviment en línia recta.",
        'moviles_total': 'Desplaçament',
        'moviles_velocidad': 'Velocitat',
        'moviles_angulo': 'Direcció',
        'moviles_rms': 'rms',
        'moviles_estado': 'Estat',
        'moviles_se_mueve': '➜ Es mou',
        'moviles_fija': 'Fixa',
        'moviles_eje_este': '← Est   (″)',
        'moviles_eje_norte': 'Nord ↑ (″)',
        'moviles_trayectoria': 'Trajectòria de {letra}',
        'moviles_pocas': 'Cal haver mesurat alguna estrella en 2 fotos o més.',
        'moviles_skybot': 'Identificar amb SkyBoT (asteroides coneguts)…',
        'moviles_elige': 'Selecciona primer un objecte a la taula.',
        'moviles_skybot_sin_wcs': "Per identificar-lo calen fotos amb solució astromètrica i data.",
        'moviles_skybot_nada': 'SkyBoT no coneix cap asteroide ni cometa a prop.',
        'moviles_skybot_cabecera': 'Objectes coneguts més propers a {letra} (SkyBoT):',
        'moviles_skybot_fila': '  {nombre}  ({clase})  V={mag}  a {d}″',
        'como_anterior_moviles': 'Objectes en moviment buscats on havien de ser: {letras}.',
        # --- revisió de la sèrie
        'serie_titulo': 'Revisió de la sèrie',
        'serie_pocas': 'Obre primer una sèrie de 2 fotos o més.',
        'serie_leyendo': 'Llegint {i} de {n}: {nombre}',
        'serie_resumen': "{n} fotos revisades. Fotos amb possibles problemes: {malas}. El desplaçament és respecte a "
                         'la primera foto.',
        'serie_fondo': 'Fons del cel',
        'serie_ruido': 'Soroll',
        'serie_estrellas': 'Estrelles',
        'serie_elong': 'Allargament',
        'serie_diagnostico': 'Diagnòstic',
        'serie_sin_registro': 'no alineable',
        'motivo_desenfocada': 'desenfocada',
        'motivo_movida': 'moguda (estrelles allargades)',
        'motivo_pocas_estrellas': 'poques estrelles (núvols?)',
        'motivo_fondo_alto': 'fons alt (Lluna, núvols, alba?)',
        'hoja_revision': 'Revisió de la sèrie',
        # --- catàlegs i astrometria
        'cat_titulo': 'Catàlegs i astrometria',
        'cat_explicacion': "Emillium funciona sense internet. Només es connecta quan prems «Descarregar», i només demana "
                           "les estrelles de la zona del cel de la foto. El que es descarrega es guarda a l'ordinador per "
                           'sempre i després es fa servir sense connexió.',
        'cat_campo': 'Camp de la foto: AR {ar}, Dec {dec}, radi {radio}′.',
        'cat_sin_wcs': "Aquesta foto no té solució astromètrica. Resol-la aquí sota amb Gaia o amb «WCS manual» "
                       '(3 estrelles conegudes); després ja es podran fer servir els catàlegs.',
        'cat_locales': ' Catàlegs guardats per a aquest camp ',
        'cat_propio': 'Catàleg propi (CSV)',
        'cat_descargar': 'Descarregar (internet)…',
        'cat_importar': 'Importar CSV…',
        'cat_no_guardado': 'no guardat',
        'cat_guardado': '✓ {n} objectes (descarregat el {fecha})',
        'cat_confirmar_t': 'Connexió a internet',
        'cat_confirmar': "Emillium es connectarà a {servidor} per descarregar {cat}.\n\nNomés s'envia la zona del cel: "
                         'AR {ar}°, Dec {dec}°, radi {radio}°. No s\'envia cap altra dada.\n\nEl resultat es guardarà en '
                         'aquest ordinador i es podrà fer servir sense connexió.\n\nVols continuar?',
        'cat_ya_descargando': 'Ja hi ha una descàrrega en marxa.',
        'cat_error_t': 'Catàleg',
        'cat_err_sin_conexion': "No s'ha pogut connectar amb el servidor (sense internet o servidor caigut).\n\n{detalle}",
        'cat_err_skybot_sin_fecha': 'SkyBoT necessita la data de la foto (DATE-OBS).',
        'cat_err_csv_vacio': 'El fitxer no té cap estrella vàlida.',
        'cat_err_csv_columnas': "No trobo les columnes d'AR i Dec (en graus). Columnes del fitxer: {columnas}",
        'cat_importado': "S'han importat {n} estrelles. Ja es poden fer servir sense connexió.",
        'ast_titulo': ' Solució astromètrica ',
        'ast_comprobar': 'Comprovar amb Gaia',
        'ast_ajustar': 'Ajustar una solució nova',
        'ast_aplicar': 'Fer servir aquesta solució',
        'ast_sin_wcs': 'La foto no té solució astromètrica.',
        'ast_sin_gaia': "No hi ha Gaia DR3 guardat per a aquest camp. Descarrega'l primer (a dalt).",
        'ast_pocas_estrellas': "S'han trobat massa poques estrelles a la foto.",
        'ast_sin_parejas': "Només {n} estrelles coincideixen amb Gaia: la solució pot estar molt malament o el camp és "
                           'massa buit.',
        'ast_resultado': '{n} estrelles emparellades amb Gaia (de {nf} detectades i les {nc} més brillants del catàleg; ~{azar} per atzar).\n'
                         'Desplaçament sistemàtic: ΔAR·cosδ = {dar}″  ΔDec = {ddec}″\n'
                         'Error típic (rms): {rms}″  ·  sense el desplaçament: {rms2}″\n'
                         "(Aquest rms s'afegeix com a error sistemàtic als informes del MPC.)",
        'ast_ajuste': 'Solució nova (TAN lineal): {n} de {total} estrelles, rms {rms}″\n'
                      'Escala {escala}″/px (abans {escala0}″/px) · eix Y cap a {rot}° (des del nord cap a l\'est){sip}',
        'ast_aviso_sip': '\nAtenció: la solució original té distorsió SIP i la nova és lineal; compara els rms abans d\'aplicar-la.',
        'ast_aplicada': '✓ Aplicada: les AR/Dec de les mesures d\'aquesta foto s\'han recalculat.',
        'ast_grafico': 'Residus Gaia − foto (fletxes exagerades)',
        'ast_grafico_ajuste': 'Residus de la solució nova',
        'cmp_titulo': ' Estrelles mesurades i catàlegs ',
        'cmp_boton': 'Identificar les estrelles mesurades',
        'cmp_otros': 'Variables (VSX) i catàleg propi',
        'info_rms': 'Astrometria comprovada amb Gaia: rms {rms}″{ajustada}',
        'info_ajustada': ' (solució ajustada)',
        'info_filtro': 'Filtre: {f}',
        # --- WCS manual
        'wcsm_titulo': 'Solució astromètrica amb estrelles conegudes',
        'wcsm_explicacion': "Escriu l'AR i la Dec d'almenys 3 estrelles mesurades en aquesta foto (d'un atles, "
                            'Stellarium, Aladin...). AR en hores (10 00 29,6) o graus (150,123); Dec en graus (+22 27 '
                            "24 o 22,456). Amb 4 o més es pot calcular l'error.",
        'wcsm_ar': 'AR',
        'wcsm_dec': 'Dec',
        'wcsm_calcular': 'Calcular',
        'wcsm_pocas': 'Calen almenys 3 estrelles mesurades en aquesta foto amb AR i Dec.',
        'wcsm_fuera': "Les coordenades de l'estrella {letra} no són vàlides.",
        'wcsm_resultado': '{n} estrelles · rms {rms}″ · escala {escala}″/px · eix Y cap a {rot}° · mirall: {espejo}',
        'wcsm_solo3': "Amb només 3 estrelles la solució passa exactament per elles: no se'n pot saber l'error.",
        # --- sessions
        'sesion_archivo': 'sessio',
        'sesion_tipo': "Sessió d'Emillium",
        'sesion_guardada': 'Sessió desada a:\n{ruta}',
        'sesion_reemplazar': 'Obrir una sessió substitueix les mesures actuals. Vols continuar?',
        'sesion_no_valida': "No és una sessió d'Emillium Expanded vàlida.",
        'sesion_abierta': "Sessió oberta: {n} mesures de {f} fotos.",
        'sesion_faltan': 'No es troben aquestes fotos (les mesures sí que hi són): {fotos}',
        # --- Excel
        'hoja_configuracion': 'Configuració',
        'conf_parametro': 'Paràmetre',
        'conf_valor': 'Valor',
        'conf_version': 'Programa',
        'conf_fecha': 'Data del full',
        'conf_astrometria': 'Astrometria foto {n}',
        'conf_rms': 'rms Gaia {rms}″ · solució ajustada: {ajustada}',
        'si': 'sí',
        'no': 'no',
        'ayuda_guardar': 'Desa les mesures en un full de càlcul (.xlsx): una fila per foto i columnes per estrella, '
                         'totes les mesures, la revisió de la sèrie i la configuració. També: Ctrl+S.',
        # --- MPC / ADES
        'ayuda_mpc': "Prepara un informe d'astrometria per al Minor Planet Center en format ADES (el recomanat) o "
                     'de 80 columnes, amb validació abans de desar-lo.',
        'mpc_explicacion': "Tria l'objecte i omple les dades. ADES és el format que recomana el MPC; el de 80 columnes "
                           "és l'antic. L'informe es comprova abans de desar-lo, però REVISA-HO: un informe correcte no "
                           'garanteix que el MPC accepti les observacions (també compten la identificació de l\'objecte, la '
                           'qualitat de les mesures i les seves normes vigents).',
        'mpc_formato': 'Format:',
        'mpc_formato_80': '80 columnes (antic)',
        'mpc_designacion': 'Designació / trkSub:',
        'mpc_remitente': 'Remitent (ADES):',
        'mpc_apertura': 'Obertura del telescopi (m):',
        'mpc_diseno': 'Disseny (ADES):',
        'mpc_cat_ast': 'Catàleg astromètric (ADES):',
        'mpc_cat_fot': 'Catàleg fotomètric (ADES):',
        'mpc_desar': 'Desar…',
        'mpc_valido': '✓ Cap problema detectat.',
        'mpc_con_errores': "L'informe té errors i no es pot desar:",
        'mpc_no_garantia': 'Validar no vol dir que el MPC l\'accepti.',
        'mpc_comentario': 'Mesurat amb Emillium {v} (centroide; AR/Dec de la solucio WCS)',
        'ades_codigo_formato': "Codi d'observatori «{c}» no vàlid (han de ser 3 caràcters, com 568 o C51).",
        'ades_codigo_xxx': 'Sense codi d\'observatori (XXX): el MPC demanarà la ubicació. Demana un codi si envies sovint.',
        'ades_falta_remitente': 'Falta el remitent (qui envia l\'informe).',
        'ades_falta_observadores': 'Falten els observadors.',
        'ades_falta_medidores': 'Falten els mesuradors.',
        'ades_apertura_numero': "L'obertura del telescopi ha de ser un número en metres (per exemple 0,20).",
        'ades_apertura_rango': "L'obertura ({a} m) no té sentit.",
        'ades_apertura_metros': "L'obertura és {a} m: segur que és en metres i no en mil·límetres?",
        'ades_diseno_raro': 'Disseny de telescopi poc habitual: «{d}».',
        'ades_modo': 'Tipus de detector «{m}» no vàlid.',
        'ades_catalogo_raro': 'Catàleg astromètric «{c}» poc habitual.',
        'ades_falta_designacion': "Falta la designació de l'objecte.",
        'ades_trksub_formato': 'La designació pròpia «{d}» ha de tenir 1–8 lletres o xifres, sense espais.',
        'ades_designacion_80': 'En el format de 80 columnes la designació té com a màxim 7 caràcters.',
        'ades_sin_observaciones': 'No hi ha cap observació amb AR/Dec i data.',
        'ades_pocas_observaciones': 'Només hi ha {n} observació: el MPC en demana com a mínim 2.',
        'ades_recomendable_3': 'Hi ha {n} observacions: el MPC recomana 3 o més per nit.',
        'ades_ar_rango': 'Observació {i}: AR fora de rang.',
        'ades_dec_rango': 'Observació {i}: Dec fora de rang.',
        'ades_fecha_futura': 'Observació {i}: la data és del futur.',
        'ades_fecha_antigua': 'Observació {i}: la data és massa antiga.',
        'ades_fecha_repetida': 'Observació {i}: data repetida ({t}).',
        'ades_sin_rms': 'Observació {i}: sense error de posició.',
        'ades_rms_cero': 'Observació {i}: error de posició zero o negatiu.',
        'ades_rms_grande': 'Observació {i}: error de posició gran ({r}″).',
        'ades_mag_rango': 'Observació {i}: magnitud fora de rang.',
        'ades_banda': 'Banda «{b}» no vàlida.',
        'ades_sin_photcat': 'Hi ha magnituds però no el catàleg fotomètric.',
        'ades_movimiento_irregular': "Les posicions no segueixen un moviment uniforme (desviació de {r}″): potser alguna "
                                     "mesura és d'un altre objecte.",
    },
    'es': {
        # --- fotometría
        'paso4': ' 4. Fotometría diferencial (opcional) ',
        'fotometria_texto': 'Escribe las estrellas de comparación. Con su magnitud de catálogo (B=12,35) se obtienen '
                            'magnitudes reales; solo con la letra (B C D), magnitudes relativas para curvas de luz.',
        'comparaciones': 'Comparación:',
        'ayuda_comparaciones': 'Letras de las estrellas de comparación separadas por espacios, con su magnitud de '
                               'catálogo y, si quieres, su índice de color: «B=12,35 C=11,90/0,65 D=13,1». Si alguna no '
                               'tiene magnitud, todo se hace en modo relativo. Tienen que estar medidas en cada foto.',
        'boton_fotometria': 'Fotometría y curvas de luz…',
        'ayuda_boton_fotometria': 'Curvas de luz, estrellas variables, comparaciones automáticas y magnitudes de los '
                                  'catálogos guardados.',
        'fot_formato': 'No entiendo las comparaciones. Ejemplo: B=12,35 C=11,90/0,65 D',
        'fot_relativa': 'Fotometría RELATIVA con {letras} (falta la magnitud de catálogo de alguna).',
        'fot_calibrada': 'Magnitudes calibradas con {n} comparación/es: {letras}.',
        'fot_sin_comp': 'Todavía no hay estrellas de comparación: se muestran magnitudes instrumentales.',
        'col_mag_rel': 'Mag (rel.)',
        'mag_cal_col': '{l}_Mag cal.',
        'mag_rel_col': '{l}_Mag rel.',
        'sigma_cal_col': '{l}_σMag dif.',
        'detalle_mag_cal': 'Magnitud calibrada: {mag} ± {sm} (con {comp})',
        'detalle_mag_rel': 'Magnitud relativa: {mag} ± {sm} (con {comp})',
        'detalle_snr': 'Señal/ruido: {snr} (de toda la estrella, apertura de {r} px)',
        'fot_titulo': 'Fotometría y curvas de luz',
        'fot_explicacion': 'Magnitud de cada estrella a lo largo de la serie. «χ²» compara la dispersión con la que '
                           'esperaríamos solo por los errores: mucho mayor que 1 significa que la estrella cambia de brillo '
                           '(o que alguna medida es mala). Selecciona filas para ver su curva de luz.',
        'fot_coef_color': 'Coef. de color T:',
        'ayuda_coef_color': 'Coeficiente de color de tu cámara y filtro (normalmente entre −0,2 y 0,2). Con 0 no se '
                            'corrige el color. Hay que poner también los índices de color de las estrellas.',
        'fot_colores': 'Índices de color:',
        'ayuda_colores': 'Índice de color (por ejemplo B−V) de las estrellas objetivo: «A=0,65 D=1,10». Se usa con el '
                         'coeficiente T.',
        'fot_n': 'Fotos',
        'fot_media': 'Mag media',
        'fot_desv': 'Desviación',
        'fot_amplitud': 'Amplitud',
        'fot_chi2': 'χ² red.',
        'fot_diagnostico': 'Diagnóstico',
        'fot_es_comp': 'Comparación (se calcula con las demás)',
        'fot_saturada': 'Saturada: la magnitud no es fiable',
        'fot_pocas': 'Demasiado pocas fotos (mínimo 3)',
        'fot_variable': '⚠ Parece VARIABLE',
        'fot_constante': 'Constante (dentro de los errores)',
        'fot_actualizar': 'Actualizar',
        'fot_auto': 'Elegir comparaciones automáticamente',
        'fot_auto_pocas': 'Hacen falta al menos 3 estrellas medidas en todas las fotos (sin saturar) para elegir comparaciones.',
        'fot_banda': 'Banda:',
        'fot_catalogo': 'Magnitudes del catálogo',
        'fot_cat_sin_comps': 'Escribe primero las letras de las estrellas de comparación.',
        'fot_cat_sin_cache': 'No hay {cat} guardado para este campo. Descárgalo en «Catálogos y astrometría».',
        'fot_cat_ok': 'Se han puesto las magnitudes de {n} estrellas con {cat}.',
        'fot_cat_faltan': 'Sin pareja en el catálogo (o sin AR/Dec): {letras}.',
        'fot_exportar': 'Exportar curva (CSV)…',
        'fot_archivo': 'curva_de_luz.csv',
        'eje_horas': 'Horas desde {fecha} UTC',
        'eje_foto': 'Número de foto',
        # --- opciones avanzadas
        'adv_r_fot': 'Radio fotométrico (px)',
        'adv_r_int': 'Anillo del cielo: interior (px)',
        'adv_r_ext': 'Anillo del cielo: exterior (px)',
        'ayuda_r_fot': 'Radio del círculo donde se suma la luz de la estrella (fotometría de apertura). Usa el mismo en '
                       'todas las estrellas. Un buen valor: 1,5–2 veces la FWHM.',
        'ayuda_r_anillo': 'El anillo donde se mide el fondo del cielo, alrededor de la estrella. Tiene que quedar fuera '
                          'de la estrella y, si puede ser, lejos de otras estrellas.',
        # --- calidad
        'col_calidad': 'Calidad',
        'calidad_col': '{l}_Calidad',
        'calidad_bien': 'Calidad: ✓ ningún problema detectado',
        'calidad_titulo': 'Calidad:',
        'calidad_SAT': 'Saturada: la estrella llega al límite de la cámara; el flujo y el centro no son fiables.',
        'calidad_CON': 'Contaminada: una estrella vecina mete más de un 1 % de luz en la apertura.',
        'calidad_BOR': 'Toca el borde de la foto o una zona sin datos.',
        'calidad_SNR': 'Señal/ruido baja (< 10): magnitud con más de 0,1 mag de error.',
        'calidad_FWHM': 'Tamaño raro comparado con las demás estrellas (¿galaxia, rayo cósmico o dos estrellas juntas?).',
        'calidad_ELO': 'Estrella alargada (movida o doble sin separar).',
        'calidad_corto_SAT': 'saturada',
        'calidad_corto_CON': 'contaminada',
        'calidad_corto_BOR': 'en el borde',
        'calidad_corto_SNR': 'S/R baja',
        'calidad_corto_FWHM': 'tamaño raro',
        'calidad_corto_ELO': 'alargada',
        'col_snr': 'S/R',
        'col_r_fot': 'Radio fot. (px)',
        'col_fondo': 'Fondo del cielo',
        # --- panel 5
        'paso5': ' 5. Análisis ',
        'boton_detectar': 'Detectar estrellas',
        'ayuda_detectar': 'Busca todas las estrellas de la foto y las marca (cian; rojo = saturada). Vuelve a pulsar '
                          'para quitar las marcas. También sirve para comprobar la astrometría.',
        'detectadas': 'Se han detectado {n} estrellas (FWHM típica {fwhm} px).',
        'boton_dobles': 'Estrellas dobles…',
        'ayuda_dobles': 'Separación y ángulo de posición entre dos estrellas medidas, en cada foto y de media.',
        'boton_moviles': 'Objetos móviles…',
        'ayuda_moviles': 'Qué estrellas medidas se mueven respecto a las demás (asteroides, cometas), a qué velocidad y '
                         'hacia dónde.',
        'boton_serie': 'Revisar la serie…',
        'ayuda_serie_rev': 'Revisa todas las fotos de la serie: fondo, estrellas, nitidez y desplazamiento, y marca las '
                           'malas (desenfocadas, movidas, con nubes...).',
        'boton_catalogos': 'Catálogos y astrometría…',
        'ayuda_catalogos': 'Gaia DR3, APASS, VSX y asteroides conocidos (SkyBoT). Comprueba y mejora la solución '
                           'astrométrica e identifica las estrellas medidas.',
        'boton_wcs_manual': 'WCS manual…',
        'ayuda_wcs_manual': 'Calcula la solución astrométrica (AR/Dec) a partir de 3 o más estrellas medidas cuyas '
                            'coordenadas conoces.',
        'boton_guardar_sesion': 'Guardar sesión…',
        'boton_abrir_sesion': 'Abrir sesión…',
        'ayuda_sesion': 'Guarda (o recupera) todo el trabajo: fotos, medidas, comparaciones, opciones y soluciones '
                        'astrométricas, para seguir otro día.',
        'necesita_foto': 'Abre primero una foto.',
        # --- dobles
        'dobles_titulo': 'Estrellas dobles',
        'dobles_explicacion': 'Separación y ángulo de posición (desde el norte hacia el este) de la secundaria vista desde '
                              'la principal. En segundos de arco hace falta que las fotos tengan solución astrométrica. El '
                              'error es el error de la media (desviación / √n).',
        'dobles_principal': 'Principal:',
        'dobles_secundaria': 'Secundaria:',
        'dobles_calcular': 'Calcular',
        'dobles_sep_px': 'Separación (px)',
        'dobles_sep_arc': 'Separación (″)',
        'dobles_pa': 'Ángulo (°)',
        'dobles_n': '{n} fotos con las dos estrellas.',
        'dobles_media_px': 'Separación media: {m} ± {e} px',
        'dobles_media_arc': 'ρ = {m} ± {e}″    θ = {pa} ± {epa}°',
        'dobles_sin_wcs': 'Sin solución astrométrica: solo separación en píxeles.',
        'dobles_pocas': 'Hay que haber medido al menos dos estrellas.',
        # --- móviles
        'moviles_titulo': 'Objetos móviles',
        'moviles_explicacion': 'Para cada letra medida en 2 fotos o más: cuánto se ha movido respecto a las estrellas. '
                               'Con solución astrométrica, en segundos de arco; si no, en píxeles respecto a la estrella '
                               'más brillante. «rms» es la distancia típica a un movimiento en línea recta.',
        'moviles_total': 'Desplazamiento',
        'moviles_velocidad': 'Velocidad',
        'moviles_angulo': 'Dirección',
        'moviles_rms': 'rms',
        'moviles_estado': 'Estado',
        'moviles_se_mueve': '➜ Se mueve',
        'moviles_fija': 'Fija',
        'moviles_eje_este': '← Este   (″)',
        'moviles_eje_norte': 'Norte ↑ (″)',
        'moviles_trayectoria': 'Trayectoria de {letra}',
        'moviles_pocas': 'Hay que haber medido alguna estrella en 2 fotos o más.',
        'moviles_skybot': 'Identificar con SkyBoT (asteroides conocidos)…',
        'moviles_elige': 'Selecciona primero un objeto en la tabla.',
        'moviles_skybot_sin_wcs': 'Para identificarlo hacen falta fotos con solución astrométrica y fecha.',
        'moviles_skybot_nada': 'SkyBoT no conoce ningún asteroide ni cometa cerca.',
        'moviles_skybot_cabecera': 'Objetos conocidos más cercanos a {letra} (SkyBoT):',
        'moviles_skybot_fila': '  {nombre}  ({clase})  V={mag}  a {d}″',
        'como_anterior_moviles': 'Objetos en movimiento buscados donde tenían que estar: {letras}.',
        # --- revisión de la serie
        'serie_titulo': 'Revisión de la serie',
        'serie_pocas': 'Abre primero una serie de 2 fotos o más.',
        'serie_leyendo': 'Leyendo {i} de {n}: {nombre}',
        'serie_resumen': '{n} fotos revisadas. Fotos con posibles problemas: {malas}. El desplazamiento es respecto a la '
                         'primera foto.',
        'serie_fondo': 'Fondo del cielo',
        'serie_ruido': 'Ruido',
        'serie_estrellas': 'Estrellas',
        'serie_elong': 'Alargamiento',
        'serie_diagnostico': 'Diagnóstico',
        'serie_sin_registro': 'no alineable',
        'motivo_desenfocada': 'desenfocada',
        'motivo_movida': 'movida (estrellas alargadas)',
        'motivo_pocas_estrellas': 'pocas estrellas (¿nubes?)',
        'motivo_fondo_alto': 'fondo alto (¿Luna, nubes, amanecer?)',
        'hoja_revision': 'Revisión de la serie',
        # --- catálogos y astrometría
        'cat_titulo': 'Catálogos y astrometría',
        'cat_explicacion': 'Emillium funciona sin internet. Solo se conecta cuando pulsas «Descargar», y solo pide las '
                           'estrellas de la zona del cielo de la foto. Lo que se descarga se guarda en el ordenador para '
                           'siempre y después se usa sin conexión.',
        'cat_campo': 'Campo de la foto: AR {ar}, Dec {dec}, radio {radio}′.',
        'cat_sin_wcs': 'Esta foto no tiene solución astrométrica. Resuélvela aquí debajo con Gaia o con «WCS manual» '
                       '(3 estrellas conocidas); después ya se podrán usar los catálogos.',
        'cat_locales': ' Catálogos guardados para este campo ',
        'cat_propio': 'Catálogo propio (CSV)',
        'cat_descargar': 'Descargar (internet)…',
        'cat_importar': 'Importar CSV…',
        'cat_no_guardado': 'no guardado',
        'cat_guardado': '✓ {n} objetos (descargado el {fecha})',
        'cat_confirmar_t': 'Conexión a internet',
        'cat_confirmar': 'Emillium se conectará a {servidor} para descargar {cat}.\n\nSolo se envía la zona del cielo: '
                         'AR {ar}°, Dec {dec}°, radio {radio}°. No se envía ningún otro dato.\n\nEl resultado se guardará en '
                         'este ordenador y se podrá usar sin conexión.\n\n¿Quieres continuar?',
        'cat_ya_descargando': 'Ya hay una descarga en marcha.',
        'cat_error_t': 'Catálogo',
        'cat_err_sin_conexion': 'No se ha podido conectar con el servidor (sin internet o servidor caído).\n\n{detalle}',
        'cat_err_skybot_sin_fecha': 'SkyBoT necesita la fecha de la foto (DATE-OBS).',
        'cat_err_csv_vacio': 'El archivo no tiene ninguna estrella válida.',
        'cat_err_csv_columnas': 'No encuentro las columnas de AR y Dec (en grados). Columnas del archivo: {columnas}',
        'cat_importado': 'Se han importado {n} estrellas. Ya se pueden usar sin conexión.',
        'ast_titulo': ' Solución astrométrica ',
        'ast_comprobar': 'Comprobar con Gaia',
        'ast_ajustar': 'Ajustar una solución nueva',
        'ast_aplicar': 'Usar esta solución',
        'ast_sin_wcs': 'La foto no tiene solución astrométrica.',
        'ast_sin_gaia': 'No hay Gaia DR3 guardado para este campo. Descárgalo primero (arriba).',
        'ast_pocas_estrellas': 'Se han encontrado demasiado pocas estrellas en la foto.',
        'ast_sin_parejas': 'Solo {n} estrellas coinciden con Gaia: la solución puede estar muy mal o el campo está '
                           'demasiado vacío.',
        'ast_resultado': '{n} estrellas emparejadas con Gaia (de {nf} detectadas y las {nc} más brillantes del catálogo; ~{azar} por azar).\n'
                         'Desplazamiento sistemático: ΔAR·cosδ = {dar}″  ΔDec = {ddec}″\n'
                         'Error típico (rms): {rms}″  ·  sin el desplazamiento: {rms2}″\n'
                         '(Este rms se añade como error sistemático a los informes del MPC.)',
        'ast_ajuste': 'Solución nueva (TAN lineal): {n} de {total} estrellas, rms {rms}″\n'
                      'Escala {escala}″/px (antes {escala0}″/px) · eje Y hacia {rot}° (desde el norte hacia el este){sip}',
        'ast_aviso_sip': '\nAtención: la solución original tiene distorsión SIP y la nueva es lineal; compara los rms antes de aplicarla.',
        'ast_aplicada': '✓ Aplicada: las AR/Dec de las medidas de esta foto se han recalculado.',
        'ast_grafico': 'Residuos Gaia − foto (flechas exageradas)',
        'ast_grafico_ajuste': 'Residuos de la solución nueva',
        'cmp_titulo': ' Estrellas medidas y catálogos ',
        'cmp_boton': 'Identificar las estrellas medidas',
        'cmp_otros': 'Variables (VSX) y catálogo propio',
        'info_rms': 'Astrometría comprobada con Gaia: rms {rms}″{ajustada}',
        'info_ajustada': ' (solución ajustada)',
        'info_filtro': 'Filtro: {f}',
        # --- WCS manual
        'wcsm_titulo': 'Solución astrométrica con estrellas conocidas',
        'wcsm_explicacion': 'Escribe la AR y la Dec de al menos 3 estrellas medidas en esta foto (de un atlas, '
                            'Stellarium, Aladin...). AR en horas (10 00 29,6) o grados (150,123); Dec en grados (+22 27 24 '
                            'o 22,456). Con 4 o más se puede calcular el error.',
        'wcsm_ar': 'AR',
        'wcsm_dec': 'Dec',
        'wcsm_calcular': 'Calcular',
        'wcsm_pocas': 'Hacen falta al menos 3 estrellas medidas en esta foto con AR y Dec.',
        'wcsm_fuera': 'Las coordenadas de la estrella {letra} no son válidas.',
        'wcsm_resultado': '{n} estrellas · rms {rms}″ · escala {escala}″/px · eje Y hacia {rot}° · espejo: {espejo}',
        'wcsm_solo3': 'Con solo 3 estrellas la solución pasa exactamente por ellas: no se puede saber su error.',
        # --- sesiones
        'sesion_archivo': 'sesion',
        'sesion_tipo': 'Sesión de Emillium',
        'sesion_guardada': 'Sesión guardada en:\n{ruta}',
        'sesion_reemplazar': 'Abrir una sesión sustituye las medidas actuales. ¿Quieres continuar?',
        'sesion_no_valida': 'No es una sesión de Emillium Expanded válida.',
        'sesion_abierta': 'Sesión abierta: {n} medidas de {f} fotos.',
        'sesion_faltan': 'No se encuentran estas fotos (las medidas sí están): {fotos}',
        # --- Excel
        'hoja_configuracion': 'Configuración',
        'conf_parametro': 'Parámetro',
        'conf_valor': 'Valor',
        'conf_version': 'Programa',
        'conf_fecha': 'Fecha de la hoja',
        'conf_astrometria': 'Astrometría foto {n}',
        'conf_rms': 'rms Gaia {rms}″ · solución ajustada: {ajustada}',
        'si': 'sí',
        'no': 'no',
        'ayuda_guardar': 'Guarda las medidas en una hoja de cálculo (.xlsx): una fila por foto y columnas por estrella, '
                         'todas las medidas, la revisión de la serie y la configuración. También: Ctrl+S.',
        # --- MPC / ADES
        'ayuda_mpc': 'Prepara un informe de astrometría para el Minor Planet Center en formato ADES (el recomendado) o '
                     'de 80 columnas, con validación antes de guardarlo.',
        'mpc_explicacion': 'Elige el objeto y rellena los datos. ADES es el formato que recomienda el MPC; el de 80 '
                           'columnas es el antiguo. El informe se comprueba antes de guardarlo, pero REVÍSALO: un informe '
                           'correcto no garantiza que el MPC acepte las observaciones (también cuentan la identificación '
                           'del objeto, la calidad de las medidas y sus normas vigentes).',
        'mpc_formato': 'Formato:',
        'mpc_formato_80': '80 columnas (antiguo)',
        'mpc_designacion': 'Designación / trkSub:',
        'mpc_remitente': 'Remitente (ADES):',
        'mpc_apertura': 'Apertura del telescopio (m):',
        'mpc_diseno': 'Diseño (ADES):',
        'mpc_cat_ast': 'Catálogo astrométrico (ADES):',
        'mpc_cat_fot': 'Catálogo fotométrico (ADES):',
        'mpc_desar': 'Guardar…',
        'mpc_valido': '✓ Ningún problema detectado.',
        'mpc_con_errores': 'El informe tiene errores y no se puede guardar:',
        'mpc_no_garantia': 'Validar no significa que el MPC lo acepte.',
        'mpc_comentario': 'Medido con Emillium {v} (centroide; AR/Dec de la solucion WCS)',
        'ades_codigo_formato': 'Código de observatorio «{c}» no válido (tienen que ser 3 caracteres, como 568 o C51).',
        'ades_codigo_xxx': 'Sin código de observatorio (XXX): el MPC pedirá la ubicación. Pide un código si envías a menudo.',
        'ades_falta_remitente': 'Falta el remitente (quién envía el informe).',
        'ades_falta_observadores': 'Faltan los observadores.',
        'ades_falta_medidores': 'Faltan los medidores.',
        'ades_apertura_numero': 'La apertura del telescopio tiene que ser un número en metros (por ejemplo 0,20).',
        'ades_apertura_rango': 'La apertura ({a} m) no tiene sentido.',
        'ades_apertura_metros': 'La apertura es {a} m: ¿seguro que está en metros y no en milímetros?',
        'ades_diseno_raro': 'Diseño de telescopio poco habitual: «{d}».',
        'ades_modo': 'Tipo de detector «{m}» no válido.',
        'ades_catalogo_raro': 'Catálogo astrométrico «{c}» poco habitual.',
        'ades_falta_designacion': 'Falta la designación del objeto.',
        'ades_trksub_formato': 'La designación propia «{d}» tiene que tener 1–8 letras o cifras, sin espacios.',
        'ades_designacion_80': 'En el formato de 80 columnas la designación tiene como máximo 7 caracteres.',
        'ades_sin_observaciones': 'No hay ninguna observación con AR/Dec y fecha.',
        'ades_pocas_observaciones': 'Solo hay {n} observación: el MPC pide como mínimo 2.',
        'ades_recomendable_3': 'Hay {n} observaciones: el MPC recomienda 3 o más por noche.',
        'ades_ar_rango': 'Observación {i}: AR fuera de rango.',
        'ades_dec_rango': 'Observación {i}: Dec fuera de rango.',
        'ades_fecha_futura': 'Observación {i}: la fecha es del futuro.',
        'ades_fecha_antigua': 'Observación {i}: la fecha es demasiado antigua.',
        'ades_fecha_repetida': 'Observación {i}: fecha repetida ({t}).',
        'ades_sin_rms': 'Observación {i}: sin error de posición.',
        'ades_rms_cero': 'Observación {i}: error de posición cero o negativo.',
        'ades_rms_grande': 'Observación {i}: error de posición grande ({r}″).',
        'ades_mag_rango': 'Observación {i}: magnitud fuera de rango.',
        'ades_banda': 'Banda «{b}» no válida.',
        'ades_sin_photcat': 'Hay magnitudes pero no el catálogo fotométrico.',
        'ades_movimiento_irregular': 'Las posiciones no siguen un movimiento uniforme (desviación de {r}″): quizá alguna '
                                     'medida es de otro objeto.',
    },
}
for _idioma, _textos in TEXTOS_15.items():
    TEXTOS[_idioma].update(_textos)


# Textos de la segunda parte de la 1.5 (resolver, series, calibración, interfaz, PDF...)
TEXTOS_15B = {
    'ca': {
        'nota_plano_sin_wcs': "Aquest pla del fitxer no fa servir la solució astromètrica de la capçalera (és la del primer pla).",
        'err_indice': 'El fitxer només té {n} imatges.',
        'imagen_n': 'Imatge {i} de {n}',
        'ayuda_imagen': 'Aquest fitxer FITS té diverses imatges (extensions o plans, per exemple diverses càmeres o '
                        'filtres). Canvia d\'imatge amb ◀ ▶.',
        'boton_medir_serie': 'Mesurar tota la sèrie',
        'ayuda_medir_serie': 'Mesura d\'un clic totes les fotos següents de la sèrie, amb les mateixes lletres que aquesta '
                             '(els objectes mòbils es busquen on han de ser).',
        'serie_medir_falta': 'Mesura primer les estrelles en aquesta foto de la sèrie.',
        'serie_midiendo': 'Mesurant la foto {i} de {n}…',
        'serie_medida': "S'han mesurat {n} fotos més ({medidas} mesures). Revisa la columna Qualitat.",
        'boton_parpadeo': 'Parpelleig',
        'ayuda_parpadeo': 'Alterna aquesta foto i l\'anterior de la sèrie, alineades: el que es mou «salta» i es veu de '
                          'seguida. Torna a prémer (o fes clic a la foto) per parar.',
        'parpadeo_falta': 'Cal una sèrie de 2 fotos o més.',
        'parpadeo_tamano': 'Les dues fotos no tenen la mateixa mida.',
        'parpadeo_mostrando': 'Parpelleig: {nombre}',
        'buscar_mov_boton': 'Buscar objectes mòbils',
        'ayuda_buscar_mov': 'Busca automàticament objectes que es moguin en línia recta entre les fotos (asteroides, '
                            'cometes), també els que no has mesurat.',
        'buscar_mov_falta': 'Calen almenys 3 fotos revisades.',
        'buscar_mov_n': '{n} objecte/s mòbil/s trobat/s.',
        'buscar_mov_vel': 'Velocitat',
        'buscar_mov_pos': 'Posició',
        'buscar_mov_letra': 'Mesurar amb la lletra:',
        'buscar_mov_medir': 'Mesurar a totes les fotos',
        'buscar_mov_medido': "S'ha mesurat l'objecte com a {letra} en {n} fotos.",
        'cal_boton': 'Calibrar amb tot el camp…',
        'ayuda_cal': 'Ajusta el punt zero i el coeficient de color amb totes les estrelles de la foto que són al catàleg '
                     '(APASS o Gaia, guardat), no només amb les comparacions.',
        'cal_titulo': 'Calibració amb el camp',
        'cal_pocas': "Només {n} estrelles útils (no saturades, S/R > 20 i al catàleg): en calen més.",
        'cal_resultado': '{n} de {total} estrelles ({cat}) · punt zero ZP = {zp} ± {szp} · rms {rms} mag',
        'cal_color': 'Coeficient de color T = {t} ± {st}',
        'cal_usadas': 'fetes servir',
        'cal_descartadas': 'descartades (> 3σ)',
        'cal_eje_color': 'Índex de color',
        'cal_eje_mag': 'Magnitud del catàleg',
        'cal_usar_t': 'Fer servir aquest coeficient de color',
        'cal_t_puesto': "Coeficient de color T = {t}. Escriu els índexs de color de les estrelles objectiu a «Índexs de color».",
        'ext_boton': 'Extinció de la sèrie…',
        'ayuda_ext': "Calcula quantes magnituds perd la llum per cada massa d'aire (coeficient k), amb el punt zero de "
                     'cada foto i la seva massa d\'aire (AIRMASS).',
        'ext_titulo': 'Extinció atmosfèrica',
        'ext_falta': "Calen 3 fotos o més amb massa d'aire (AIRMASS) i catàleg guardat (n'hi ha {n}).",
        'ext_resultado': 'k = {k} ± {sk} mag per massa d\'aire · {n} fotos · rms {rms} mag',
        'ext_eje': "Massa d'aire",
        'ropt_boton': 'Radi òptim…',
        'ayuda_ropt': 'Corba de creixement de les estrelles mesurades: quin radi d\'obertura dona més senyal/soroll.',
        'ropt_titulo': "Radi d'obertura òptim",
        'ropt_falta': 'Mesura primer alguna estrella (no saturada) en aquesta foto.',
        'ropt_resultado': 'Radi amb més senyal/soroll: {r} px ({f} vegades la FWHM, {fwhm} px). Per a fotometria precisa '
                          "d'estrelles brillants es pot fer servir una mica més gran.",
        'ropt_flujo': 'Flux (fracció)',
        'ropt_snr': 'S/R (fracció)',
        'ropt_radio': 'Radi (px)',
        'ropt_usar': 'Fer servir {r} px',
        'capas_mostrar': 'Mostrar a la foto:',
        'capa_vsx': 'variables',
        'capa_ngc': 'NGC/IC',
        'capa_skybot': 'asteroides',
        'capa_gaia': 'Gaia',
        'capas_n': '{n} marques',
        'resolver_titulo': ' Resoldre la foto amb Gaia ',
        'resolver_explicacion': "La foto no té solució astromètrica. Si saps més o menys cap on apuntava el telescopi (i, "
                                "si pot ser, l'escala), Emillium la pot trobar comparant les figures que formen les estrelles "
                                'amb les de Gaia. Primer descarrega Gaia per a la zona (cal internet, un cop) i després prem '
                                '«Resoldre» (sense internet).',
        'resolver_escala': 'Escala (″/px, opcional):',
        'resolver_descargar': 'Descarregar Gaia per a la zona…',
        'resolver_boton': 'Resoldre',
        'resolver_falta_pos': "Escriu l'AR i la Dec aproximades.",
        'resolver_descargado': '✓ {n} estrelles de Gaia guardades. Ara prem «Resoldre».',
        'resolver_buscando': 'Buscant… (pot tardar uns segons)',
        'resolver_sin_gaia': 'No hi ha Gaia guardat per a aquesta zona: descarrega\'l primer.',
        'resolver_no': "No s'ha trobat la solució. Comprova la posició aproximada, prova sense escala o descarrega una "
                       'zona més gran.',
        'resolver_ok': 'Foto resolta amb {n} estrelles · rms {rms}″ · escala {escala}″/px · eix Y cap a {rot}° · mirall: {espejo}',
        'aperturas': 'Obertures',
        'ayuda_aperturas': "Dibuixa l'obertura fotomètrica i l'anell del cel de cada estrella mesurada. Les estrelles amb "
                           'algun problema de qualitat es marquen en taronja.',
        'modo': 'Mode:',
        'modo_avanzado': 'Mode avançat',
        'modo_sencillo': 'Mode senzill',
        'ayuda_modo': 'Senzill: només el bàsic (obrir, girar, mesurar, fotometria i desar). Avançat: totes les eines.',
        'boton_info_foto': 'Informació de la foto…',
        'ayuda_info_foto': 'Estadístiques de la foto (fons, soroll, estrelles, FWHM), histograma i capçalera FITS.',
        'info_titulo': 'Informació · {nombre}',
        'info_tamano': 'Mida',
        'info_fondo': 'Fons del cel',
        'info_ruido': 'Soroll',
        'info_estrellas': 'Estrelles detectades',
        'info_saturacion': 'Nivell de saturació',
        'info_saturadas': 'Estrelles saturades',
        'info_escala': 'Escala',
        'info_masa': "Massa d'aire",
        'info_histograma': 'Valor dels píxels (escala logarítmica)',
        'info_buscar': 'Buscar a la capçalera:',
        'info_clave': 'Paraula clau',
        'boton_historial': 'Historial…',
        'ayuda_historial': "Totes les operacions fetes (fotos obertes, mesures, solucions, descàrregues, desats...).",
        'hist_titulo': 'Historial',
        'hist_hora': 'Hora',
        'hist_operacion': 'Operació',
        'hist_abrir': 'Foto oberta: {nombre}',
        'hist_medida': 'Mesura: estrella {letra} a la foto {n} (X {x}, Y {y})',
        'hist_serie': 'Sèrie mesurada: {n} fotos',
        'hist_parpadeo': 'Parpelleig entre {a} i {b}',
        'hist_buscar_mov': 'Cerca d\'objectes mòbils: {n} trobats',
        'hist_medir_pista': 'Objecte mòbil mesurat com a {letra} en {n} fotos',
        'hist_calibracion': 'Calibració amb el camp: ZP {zp} ({n} estrelles)',
        'hist_coef_color': 'Coeficient de color T = {t}',
        'hist_extincion': 'Extinció: k = {k} ({n} fotos)',
        'hist_radio': 'Radi fotomètric: {r} px',
        'hist_resuelta': 'Foto resolta amb Gaia: {n} estrelles, rms {rms}″',
        'hist_wcs': 'Solució astromètrica nova a la foto {n} (rms {rms}″)',
        'hist_descarga': 'Descàrrega de {cat}: {n} objectes (correcta: {ok})',
        'hist_revision': 'Sèrie revisada: {n} fotos, {malas} amb problemes',
        'hist_excel': 'Desat per a Excel: {ruta} ({n} mesures)',
        'hist_mpc': 'Informe MPC ({formato}): {ruta}',
        'hist_pdf': 'Informe PDF: {ruta}',
        'hist_sesion_guardada': 'Sessió desada: {ruta}',
        'hist_sesion_abierta': 'Sessió oberta: {ruta}',
        'boton_pdf': 'Informe PDF…',
        'ayuda_pdf': 'Desa un informe en PDF: resum, taula de mesures, corbes de llum, objectes mòbils i revisió de la sèrie.',
        'pdf_titulo': 'Informe PDF',
        'pdf_archivo': 'informe_emillium',
        'pdf_resumen': "Informe d'Emillium Expanded",
        'pdf_tabla': 'Mesures',
        'pdf_curvas': 'Corbes de llum',
        'pdf_medidas': 'mesures',
        'pdf_fotos': 'fotos',
        'boton_guia': 'Guia pas a pas…',
        'guia_titulo': 'Guia pas a pas',
        'guia_explicacion': "Els passos per treure profit d'Emillium Expanded. Es marquen sols a mesura que els fas "
                            '(prem «Actualitzar» per tornar-los a mirar).',
        'guia_1': 'Obre una foto o una sèrie amb «Obrir foto FITS…».',
        'guia_2': 'Tria la lletra i fes clic a les estrelles (o «Detectar estrelles» per veure-les totes).',
        'guia_3': 'Amb una sèrie: passa a la foto següent i prem «Mesurar com a la foto anterior» (o «Mesurar tota la sèrie»).',
        'guia_4': 'Mira la columna Qualitat: evita les estrelles saturades (SAT) i contaminades (CON).',
        'guia_5': 'Fotometria: escriu les estrelles de comparació (per exemple B=12,35 C D) i obre «Fotometria i corbes de llum».',
        'guia_6': 'Coordenades: la foto ha de tenir solució astromètrica (o resol-la a «Catàlegs i astrometria»).',
        'guia_7': 'Comprova la solució astromètrica amb Gaia a «Catàlegs i astrometria».',
        'guia_8': 'Desa la feina: Excel, informe PDF, informe MPC o sessió.',
    },
    'es': {
        'nota_plano_sin_wcs': 'Este plano del archivo no usa la solución astrométrica de la cabecera (es la del primer plano).',
        'err_indice': 'El archivo solo tiene {n} imágenes.',
        'imagen_n': 'Imagen {i} de {n}',
        'ayuda_imagen': 'Este archivo FITS tiene varias imágenes (extensiones o planos, por ejemplo varias cámaras o '
                        'filtros). Cambia de imagen con ◀ ▶.',
        'boton_medir_serie': 'Medir toda la serie',
        'ayuda_medir_serie': 'Mide de un clic todas las fotos siguientes de la serie, con las mismas letras que esta '
                             '(los objetos móviles se buscan donde tienen que estar).',
        'serie_medir_falta': 'Mide primero las estrellas en esta foto de la serie.',
        'serie_midiendo': 'Midiendo la foto {i} de {n}…',
        'serie_medida': 'Se han medido {n} fotos más ({medidas} medidas). Revisa la columna Calidad.',
        'boton_parpadeo': 'Parpadeo',
        'ayuda_parpadeo': 'Alterna esta foto y la anterior de la serie, alineadas: lo que se mueve «salta» y se ve '
                          'enseguida. Vuelve a pulsar (o haz clic en la foto) para parar.',
        'parpadeo_falta': 'Hace falta una serie de 2 fotos o más.',
        'parpadeo_tamano': 'Las dos fotos no tienen el mismo tamaño.',
        'parpadeo_mostrando': 'Parpadeo: {nombre}',
        'buscar_mov_boton': 'Buscar objetos móviles',
        'ayuda_buscar_mov': 'Busca automáticamente objetos que se muevan en línea recta entre las fotos (asteroides, '
                            'cometas), también los que no has medido.',
        'buscar_mov_falta': 'Hacen falta al menos 3 fotos revisadas.',
        'buscar_mov_n': '{n} objeto/s móvil/es encontrado/s.',
        'buscar_mov_vel': 'Velocidad',
        'buscar_mov_pos': 'Posición',
        'buscar_mov_letra': 'Medir con la letra:',
        'buscar_mov_medir': 'Medir en todas las fotos',
        'buscar_mov_medido': 'Se ha medido el objeto como {letra} en {n} fotos.',
        'cal_boton': 'Calibrar con todo el campo…',
        'ayuda_cal': 'Ajusta el punto cero y el coeficiente de color con todas las estrellas de la foto que están en el '
                     'catálogo (APASS o Gaia, guardado), no solo con las comparaciones.',
        'cal_titulo': 'Calibración con el campo',
        'cal_pocas': 'Solo {n} estrellas útiles (no saturadas, S/R > 20 y en el catálogo): hacen falta más.',
        'cal_resultado': '{n} de {total} estrellas ({cat}) · punto cero ZP = {zp} ± {szp} · rms {rms} mag',
        'cal_color': 'Coeficiente de color T = {t} ± {st}',
        'cal_usadas': 'usadas',
        'cal_descartadas': 'descartadas (> 3σ)',
        'cal_eje_color': 'Índice de color',
        'cal_eje_mag': 'Magnitud del catálogo',
        'cal_usar_t': 'Usar este coeficiente de color',
        'cal_t_puesto': 'Coeficiente de color T = {t}. Escribe los índices de color de las estrellas objetivo en «Índices de color».',
        'ext_boton': 'Extinción de la serie…',
        'ayuda_ext': 'Calcula cuántas magnitudes pierde la luz por cada masa de aire (coeficiente k), con el punto cero '
                     'de cada foto y su masa de aire (AIRMASS).',
        'ext_titulo': 'Extinción atmosférica',
        'ext_falta': 'Hacen falta 3 fotos o más con masa de aire (AIRMASS) y catálogo guardado (hay {n}).',
        'ext_resultado': 'k = {k} ± {sk} mag por masa de aire · {n} fotos · rms {rms} mag',
        'ext_eje': 'Masa de aire',
        'ropt_boton': 'Radio óptimo…',
        'ayuda_ropt': 'Curva de crecimiento de las estrellas medidas: qué radio de apertura da más señal/ruido.',
        'ropt_titulo': 'Radio de apertura óptimo',
        'ropt_falta': 'Mide primero alguna estrella (no saturada) en esta foto.',
        'ropt_resultado': 'Radio con más señal/ruido: {r} px ({f} veces la FWHM, {fwhm} px). Para fotometría precisa de '
                          'estrellas brillantes se puede usar uno algo mayor.',
        'ropt_flujo': 'Flujo (fracción)',
        'ropt_snr': 'S/R (fracción)',
        'ropt_radio': 'Radio (px)',
        'ropt_usar': 'Usar {r} px',
        'capas_mostrar': 'Mostrar en la foto:',
        'capa_vsx': 'variables',
        'capa_ngc': 'NGC/IC',
        'capa_skybot': 'asteroides',
        'capa_gaia': 'Gaia',
        'capas_n': '{n} marcas',
        'resolver_titulo': ' Resolver la foto con Gaia ',
        'resolver_explicacion': 'La foto no tiene solución astrométrica. Si sabes más o menos hacia dónde apuntaba el '
                                'telescopio (y, si puede ser, la escala), Emillium la puede encontrar comparando las figuras '
                                'que forman las estrellas con las de Gaia. Primero descarga Gaia para la zona (hace falta '
                                'internet, una vez) y después pulsa «Resolver» (sin internet).',
        'resolver_escala': 'Escala (″/px, opcional):',
        'resolver_descargar': 'Descargar Gaia para la zona…',
        'resolver_boton': 'Resolver',
        'resolver_falta_pos': 'Escribe la AR y la Dec aproximadas.',
        'resolver_descargado': '✓ {n} estrellas de Gaia guardadas. Ahora pulsa «Resolver».',
        'resolver_buscando': 'Buscando… (puede tardar unos segundos)',
        'resolver_sin_gaia': 'No hay Gaia guardado para esta zona: descárgalo primero.',
        'resolver_no': 'No se ha encontrado la solución. Comprueba la posición aproximada, prueba sin escala o descarga '
                       'una zona más grande.',
        'resolver_ok': 'Foto resuelta con {n} estrellas · rms {rms}″ · escala {escala}″/px · eje Y hacia {rot}° · espejo: {espejo}',
        'aperturas': 'Aperturas',
        'ayuda_aperturas': 'Dibuja la apertura fotométrica y el anillo del cielo de cada estrella medida. Las estrellas con '
                           'algún problema de calidad se marcan en naranja.',
        'modo': 'Modo:',
        'modo_avanzado': 'Modo avanzado',
        'modo_sencillo': 'Modo sencillo',
        'ayuda_modo': 'Sencillo: solo lo básico (abrir, girar, medir, fotometría y guardar). Avanzado: todas las herramientas.',
        'boton_info_foto': 'Información de la foto…',
        'ayuda_info_foto': 'Estadísticas de la foto (fondo, ruido, estrellas, FWHM), histograma y cabecera FITS.',
        'info_titulo': 'Información · {nombre}',
        'info_tamano': 'Tamaño',
        'info_fondo': 'Fondo del cielo',
        'info_ruido': 'Ruido',
        'info_estrellas': 'Estrellas detectadas',
        'info_saturacion': 'Nivel de saturación',
        'info_saturadas': 'Estrellas saturadas',
        'info_escala': 'Escala',
        'info_masa': 'Masa de aire',
        'info_histograma': 'Valor de los píxeles (escala logarítmica)',
        'info_buscar': 'Buscar en la cabecera:',
        'info_clave': 'Palabra clave',
        'boton_historial': 'Historial…',
        'ayuda_historial': 'Todas las operaciones hechas (fotos abiertas, medidas, soluciones, descargas, guardados...).',
        'hist_titulo': 'Historial',
        'hist_hora': 'Hora',
        'hist_operacion': 'Operación',
        'hist_abrir': 'Foto abierta: {nombre}',
        'hist_medida': 'Medida: estrella {letra} en la foto {n} (X {x}, Y {y})',
        'hist_serie': 'Serie medida: {n} fotos',
        'hist_parpadeo': 'Parpadeo entre {a} y {b}',
        'hist_buscar_mov': 'Búsqueda de objetos móviles: {n} encontrados',
        'hist_medir_pista': 'Objeto móvil medido como {letra} en {n} fotos',
        'hist_calibracion': 'Calibración con el campo: ZP {zp} ({n} estrellas)',
        'hist_coef_color': 'Coeficiente de color T = {t}',
        'hist_extincion': 'Extinción: k = {k} ({n} fotos)',
        'hist_radio': 'Radio fotométrico: {r} px',
        'hist_resuelta': 'Foto resuelta con Gaia: {n} estrellas, rms {rms}″',
        'hist_wcs': 'Solución astrométrica nueva en la foto {n} (rms {rms}″)',
        'hist_descarga': 'Descarga de {cat}: {n} objetos (correcta: {ok})',
        'hist_revision': 'Serie revisada: {n} fotos, {malas} con problemas',
        'hist_excel': 'Guardado para Excel: {ruta} ({n} medidas)',
        'hist_mpc': 'Informe MPC ({formato}): {ruta}',
        'hist_pdf': 'Informe PDF: {ruta}',
        'hist_sesion_guardada': 'Sesión guardada: {ruta}',
        'hist_sesion_abierta': 'Sesión abierta: {ruta}',
        'boton_pdf': 'Informe PDF…',
        'ayuda_pdf': 'Guarda un informe en PDF: resumen, tabla de medidas, curvas de luz, objetos móviles y revisión de la serie.',
        'pdf_titulo': 'Informe PDF',
        'pdf_archivo': 'informe_emillium',
        'pdf_resumen': 'Informe de Emillium Expanded',
        'pdf_tabla': 'Medidas',
        'pdf_curvas': 'Curvas de luz',
        'pdf_medidas': 'medidas',
        'pdf_fotos': 'fotos',
        'boton_guia': 'Guía paso a paso…',
        'guia_titulo': 'Guía paso a paso',
        'guia_explicacion': 'Los pasos para sacar partido a Emillium Expanded. Se marcan solos a medida que los haces '
                            '(pulsa «Actualizar» para volver a mirarlos).',
        'guia_1': 'Abre una foto o una serie con «Abrir foto FITS…».',
        'guia_2': 'Elige la letra y haz clic en las estrellas (o «Detectar estrellas» para verlas todas).',
        'guia_3': 'Con una serie: pasa a la foto siguiente y pulsa «Medir como en la foto anterior» (o «Medir toda la serie»).',
        'guia_4': 'Mira la columna Calidad: evita las estrellas saturadas (SAT) y contaminadas (CON).',
        'guia_5': 'Fotometría: escribe las estrellas de comparación (por ejemplo B=12,35 C D) y abre «Fotometría y curvas de luz».',
        'guia_6': 'Coordenadas: la foto tiene que tener solución astrométrica (o resuélvela en «Catálogos y astrometría»).',
        'guia_7': 'Comprueba la solución astrométrica con Gaia en «Catálogos y astrometría».',
        'guia_8': 'Guarda el trabajo: Excel, informe PDF, informe MPC o sesión.',
    },
}
for _idioma, _textos in TEXTOS_15B.items():
    TEXTOS[_idioma].update(_textos)


# Textos de la tercera parte de la 1.5 (periodos, mínimos, SIP, recortes, gestor...)
TEXTOS_15C = {
    'ca': {
        'err_gz': "El fitxer .gz està fet malbé o no és un FITS comprimit amb gzip.",
        'error_inesperado': "Hi ha hagut un error inesperat, però el programa continua obert.\n\n{error}\n\n"
                            "El detall s'ha desat a %APPDATA%\\Emillium\\errores.log.",
        'adv_corregir': 'Corregir la contaminació de les veïnes',
        'ayuda_corregir': "Resta del flux la llum estimada de les estrelles veïnes que cau dins de l'obertura (suposant "
                          'estrelles gaussianes). L\'error de la mesura augmenta en proporció al que es resta.',
        'boton_recortes': 'Retalls d\'una estrella…',
        'ayuda_recortes': "La mateixa estrella a totes les fotos on s'ha mesurat, una al costat de l'altra: per veure d'un "
                          'cop d\'ull si alguna mesura és dolenta.',
        'rec_titulo': 'Retalls',
        'per_boton': 'Període…',
        'ayuda_per': "Busca el període d'una estrella variable (periodograma de Lomb-Scargle) i en mostra la corba de llum en fase.",
        'per_titulo': 'Període',
        'per_falta': 'Cal una estrella mesurada en 6 fotos o més.',
        'per_explicacion': "Per a cada període de prova s'ajusta una sinusoide a les magnituds: la potència (0 a 1) diu quina "
                           'part de la variació explica. Compte: amb fotos d\'una sola nit, o cada nit a la mateixa hora, '
                           'poden sortir períodes falsos (àlies).',
        'per_min': 'Període mín. (dies):',
        'per_max': 'màx.:',
        'per_eje_periodo': 'Període (dies)',
        'per_eje_potencia': 'Potència',
        'per_eje_fase': 'Fase',
        'per_resultado': '{letra}: període {p} dies ({h} h) · potència {pot} · temps {tiempo}',
        'per_sin_tiempo': 'Calen fotos amb data (DATE-OBS) per buscar el període.',
        'min_boton': "Mínim d'eclipsi…",
        'ayuda_min': "Instant del mínim d'un eclipsi (mètode de Kwee i van Woerden), amb el seu error.",
        'min_titulo': "Mínim d'eclipsi",
        'min_falta': 'Cal una estrella mesurada en 7 fotos o més.',
        'min_explicacion': "Fes servir fotos de la mateixa nit al voltant del mínim (amb punts abans i després). Es busca "
                           "l'instant on la baixada i la pujada de la corba són més simètriques.",
        'min_eje': 'Hores des de {t} {t0}',
        'min_no': "No s'ha trobat el mínim: cal que el punt més feble tingui mesures abans i després.",
        'min_resultado': 'Mínim: {tiempo} {t} ± {e} s  ({fecha} UTC)',
        'col_hjd': "HJD (mig de l'exposició)",
        'ast_ajustar_sip': 'Ajustar amb distorsió (SIP)',
        'ast_sip_pocas': 'Per ajustar la distorsió calen més estrelles emparellades (ara {n}; mínim 18).',
        'ast_sip': 'Solució amb distorsió SIP d\'ordre {orden}: {n} de {total} estrelles, rms {rms}″',
        'ast_grafico_sip': 'Residus de la solució amb distorsió',
        'gestor_boton': 'Gestionar els catàlegs guardats…',
        'gestor_titulo': 'Catàlegs guardats',
        'gestor_catalogo': 'Catàleg',
        'gestor_radio': 'Radi (°)',
        'gestor_objetos': 'Objectes',
        'gestor_fecha': 'Descarregat',
        'gestor_resumen': '{n} trossos de catàleg · {mb} MB · {carpeta}',
        'gestor_borrar': 'Esborrar els seleccionats',
        'gestor_confirmar': "Esborrar {n} trossos de catàleg? Per tornar-los a tenir caldrà descarregar-los (amb internet).",
        'hist_periodo': 'Període de {letra}: {p} dies',
        'hist_minimo': 'Mínim de {letra}: {tipo} {t}',
        'hist_borrar_cat': '{n} catàlegs guardats esborrats',
    },
    'es': {
        'err_gz': 'El archivo .gz está estropeado o no es un FITS comprimido con gzip.',
        'error_inesperado': 'Ha habido un error inesperado, pero el programa sigue abierto.\n\n{error}\n\n'
                            'El detalle se ha guardado en %APPDATA%\\Emillium\\errores.log.',
        'adv_corregir': 'Corregir la contaminación de las vecinas',
        'ayuda_corregir': 'Resta del flujo la luz estimada de las estrellas vecinas que cae dentro de la apertura '
                          '(suponiendo estrellas gaussianas). El error de la medida aumenta en proporción a lo que se resta.',
        'boton_recortes': 'Recortes de una estrella…',
        'ayuda_recortes': 'La misma estrella en todas las fotos donde se ha medido, una al lado de otra: para ver de un '
                          'vistazo si alguna medida es mala.',
        'rec_titulo': 'Recortes',
        'per_boton': 'Periodo…',
        'ayuda_per': 'Busca el periodo de una estrella variable (periodograma de Lomb-Scargle) y muestra su curva de luz en fase.',
        'per_titulo': 'Periodo',
        'per_falta': 'Hace falta una estrella medida en 6 fotos o más.',
        'per_explicacion': 'Para cada periodo de prueba se ajusta una sinusoide a las magnitudes: la potencia (0 a 1) dice '
                           'qué parte de la variación explica. Cuidado: con fotos de una sola noche, o cada noche a la misma '
                           'hora, pueden salir periodos falsos (alias).',
        'per_min': 'Periodo mín. (días):',
        'per_max': 'máx.:',
        'per_eje_periodo': 'Periodo (días)',
        'per_eje_potencia': 'Potencia',
        'per_eje_fase': 'Fase',
        'per_resultado': '{letra}: periodo {p} días ({h} h) · potencia {pot} · tiempo {tiempo}',
        'per_sin_tiempo': 'Hacen falta fotos con fecha (DATE-OBS) para buscar el periodo.',
        'min_boton': 'Mínimo de eclipse…',
        'ayuda_min': 'Instante del mínimo de un eclipse (método de Kwee y van Woerden), con su error.',
        'min_titulo': 'Mínimo de eclipse',
        'min_falta': 'Hace falta una estrella medida en 7 fotos o más.',
        'min_explicacion': 'Usa fotos de la misma noche alrededor del mínimo (con puntos antes y después). Se busca el '
                           'instante en el que la bajada y la subida de la curva son más simétricas.',
        'min_eje': 'Horas desde {t} {t0}',
        'min_no': 'No se ha encontrado el mínimo: hace falta que el punto más débil tenga medidas antes y después.',
        'min_resultado': 'Mínimo: {tiempo} {t} ± {e} s  ({fecha} UTC)',
        'col_hjd': 'HJD (centro de la exposición)',
        'ast_ajustar_sip': 'Ajustar con distorsión (SIP)',
        'ast_sip_pocas': 'Para ajustar la distorsión hacen falta más estrellas emparejadas (ahora {n}; mínimo 18).',
        'ast_sip': 'Solución con distorsión SIP de orden {orden}: {n} de {total} estrellas, rms {rms}″',
        'ast_grafico_sip': 'Residuos de la solución con distorsión',
        'gestor_boton': 'Gestionar los catálogos guardados…',
        'gestor_titulo': 'Catálogos guardados',
        'gestor_catalogo': 'Catálogo',
        'gestor_radio': 'Radio (°)',
        'gestor_objetos': 'Objetos',
        'gestor_fecha': 'Descargado',
        'gestor_resumen': '{n} trozos de catálogo · {mb} MB · {carpeta}',
        'gestor_borrar': 'Borrar los seleccionados',
        'gestor_confirmar': '¿Borrar {n} trozos de catálogo? Para volver a tenerlos habrá que descargarlos (con internet).',
        'hist_periodo': 'Periodo de {letra}: {p} días',
        'hist_minimo': 'Mínimo de {letra}: {tipo} {t}',
        'hist_borrar_cat': '{n} catálogos guardados borrados',
    },
}
for _idioma, _textos in TEXTOS_15C.items():
    TEXTOS[_idioma].update(_textos)

# Textos de la última parte de la 1.5: BJD_TDB, alinear y apilar
TEXTOS_15D = {
    'ca': {
        'escala_tiempo': 'Temps:',
        'col_bjd': "BJD_TDB (mig de l'exposició)",
        'parpadeo_tamano': "No s'han pogut alinear les dues fotos (no s'hi troben prou estrelles en comú).",
        'boton_alinear': 'Alinear i apilar…',
        'ayuda_alinear': "Alinea totes les fotos de la sèrie amb la que tens oberta (desplaçament, gir i mirall, "
                         "amb precisió de centèsimes de píxel) i les suma en una foto nova amb menys soroll.",
        'ali_titulo': 'Alinear i apilar fotos',
        'ali_falta': 'Obre una sèrie de 2 fotos o més (pots triar-ne diverses alhora a «Obrir») i posa\'t a la '
                     'que vulguis fer servir de referència.',
        'ali_explicacion': "Es busquen les mateixes estrelles a cada foto i es calcula com s'ha de moure, girar o "
                           "capgirar (com en un mirall) perquè quedin a sobre de les de la referència ({ref}). Després es "
                           "poden apilar: la mitjana de N fotos té √N vegades menys soroll, així que hi surten "
                           "estrelles i objectes més febles. La foto apilada es guarda en un FITS nou (les originals "
                           "no es toquen) amb la solució astromètrica de la referència.",
        'ali_alineando': 'Alineant {i} de {n}: {nombre}',
        'ali_estrellas': 'Estrelles',
        'ali_giro': 'Gir (°)',
        'ali_no': 'no alineada',
        'ali_referencia': 'referència',
        'ali_espejo': '(mirall)',
        'ali_resumen': "{n} de {total} fotos alineades. Tria a la taula les que vols apilar (Ctrl+clic) i el mètode.",
        'ali_metodo': 'Mètode:',
        'ali_media': 'Mitjana amb retall (recomanat)',
        'ali_mediana': 'Mediana',
        'ali_seleccion': "La mitjana amb retall descarta en cada píxel els valors estranys (raigs còsmics, satèl·lits, "
                         "avions) i fa la mitjana de la resta: és la millor per mesurar. La mediana encara és més "
                         "robusta, però té una mica més de soroll.",
        'ali_apilar': 'Apilar i guardar…',
        'ayuda_apilar': 'Crea una foto nova amb la mitjana de les fotos triades, alineades amb la referència.',
        'ali_pocas': 'Tria almenys una foto alineada a més de la referència.',
        'ali_sufijo': '_apilada_{n}',
        'ali_leyendo': 'Llegint {nombre}…',
        'ali_apilando': 'Apilant… {p} %',
        'ali_hecho': "Foto apilada amb {n} fotos: {ruta}",
        'ali_abrir': 'Vols obrir-la ara? (S\'afegirà al final de la sèrie.)',
        'hist_alinear': 'Alineades {n} de {total} fotos de la sèrie',
        'hist_apilar': 'Apilades {n} fotos a {nombre}',
    },
    'es': {
        'escala_tiempo': 'Tiempo:',
        'col_bjd': 'BJD_TDB (centro de la exposición)',
        'parpadeo_tamano': 'No se han podido alinear las dos fotos (no tienen suficientes estrellas en común).',
        'boton_alinear': 'Alinear y apilar…',
        'ayuda_alinear': 'Alinea todas las fotos de la serie con la que tienes abierta (desplazamiento, giro y espejo, '
                         'con precisión de centésimas de píxel) y las suma en una foto nueva con menos ruido.',
        'ali_titulo': 'Alinear y apilar fotos',
        'ali_falta': 'Abre una serie de 2 fotos o más (puedes elegir varias a la vez en «Abrir») y ponte en la '
                     'que quieras usar de referencia.',
        'ali_explicacion': 'Se buscan las mismas estrellas en cada foto y se calcula cómo hay que moverla, girarla o '
                           'darle la vuelta como un espejo para que queden encima de las de la referencia ({ref}). '
                           'Después se pueden apilar: la media de N fotos tiene √N veces menos ruido, así que salen '
                           'estrellas y objetos más débiles. La foto apilada se guarda en un FITS nuevo (las '
                           'originales no se tocan) con la solución astrométrica de la referencia.',
        'ali_alineando': 'Alineando {i} de {n}: {nombre}',
        'ali_estrellas': 'Estrellas',
        'ali_giro': 'Giro (°)',
        'ali_no': 'no alineada',
        'ali_referencia': 'referencia',
        'ali_espejo': '(espejo)',
        'ali_resumen': '{n} de {total} fotos alineadas. Elige en la tabla las que quieres apilar (Ctrl+clic) y el método.',
        'ali_metodo': 'Método:',
        'ali_media': 'Media con recorte (recomendado)',
        'ali_mediana': 'Mediana',
        'ali_seleccion': 'La media con recorte descarta en cada píxel los valores raros (rayos cósmicos, satélites, '
                         'aviones) y hace la media del resto: es la mejor para medir. La mediana es aún más robusta, '
                         'pero tiene algo más de ruido.',
        'ali_apilar': 'Apilar y guardar…',
        'ayuda_apilar': 'Crea una foto nueva con la media de las fotos elegidas, alineadas con la referencia.',
        'ali_pocas': 'Elige al menos una foto alineada además de la referencia.',
        'ali_sufijo': '_apilada_{n}',
        'ali_leyendo': 'Leyendo {nombre}…',
        'ali_apilando': 'Apilando… {p} %',
        'ali_hecho': 'Foto apilada con {n} fotos: {ruta}',
        'ali_abrir': '¿Quieres abrirla ahora? (Se añadirá al final de la serie.)',
        'hist_alinear': 'Alineadas {n} de {total} fotos de la serie',
        'hist_apilar': 'Apiladas {n} fotos en {nombre}',
    },
}
for _idioma, _textos in TEXTOS_15D.items():
    TEXTOS[_idioma].update(_textos)

ESCALAS_TIEMPO = ['BJD_TDB', 'HJD', 'JD']


LETRAS_ESTRELLA = [chr(c) for c in range(ord('A'), ord('Z') + 1)]

# Combinaciones secretas: lo que hay que teclear seguido -> easter egg
CODIGOS_SECRETOS = {
    'm42': 'nebulosa', '618': 'agujero', 'm1': 'supernova', 'crab': 'supernova', 'm16': 'pilares',
    'earth': 'tierra', '1977': 'voyager', 'iss': 'iss', '6equj5': 'wow', 'wow': 'wow', 'spock': 'vulcano',
    'vulcan': 'vulcano', 'gargantua': 'lente', 'interstellar': 'lente', '42': 'respuesta',
    '2001': 'odisea',
}
EFECTOS_DE_CIELO = ('nebulosa', 'agujero', 'supernova', 'pilares', 'voyager', 'iss',
                    'odisea')  # solo en la pantalla de inicio

COLUMNAS_TABLA = [
    # (clave, clave del título en TEXTOS, ancho, alineación)
    ('foto_n', 'col_foto_n', 58, tk.CENTER),
    ('archivo', 'col_foto', 130, tk.W),
    ('estrella', 'col_estrella', 58, tk.CENTER),
    ('orientacion', 'col_orientacion', 120, tk.W),
    ('x', None, 80, tk.E),
    ('y', None, 80, tk.E),
    ('error', 'col_error', 58, tk.E),
    ('flujo', 'col_flujo', 82, tk.E),
    ('mag', 'col_mag_inst', 78, tk.E),
    ('ar', 'col_ar', 108, tk.E),
    ('dec', 'col_dec', 108, tk.E),
    ('calidad', 'col_calidad', 82, tk.CENTER),
]


def resource_path(relative_path):
    """Ruta a un archivo del programa (compatible con el .exe de PyInstaller)."""
    if getattr(sys, 'frozen', False):
        base_path = getattr(sys, '_MEIPASS', os.path.dirname(sys.executable))
    else:
        base_path = os.path.abspath(os.path.dirname(__file__))
    return os.path.join(base_path, relative_path)


def numero(valor, decimales):
    """Número con coma decimal ('' si no hay valor)."""
    if valor is None or not np.isfinite(valor):
        return ''
    if valor != 0 and (abs(valor) < min(1e-6, 0.5 * 10.0 ** -decimales) or abs(valor) >= 1e12):
        # Fotos en unidades físicas (p. ej. erg/s/cm²/Å ≈ 1e-13): con decimales fijos
        # saldría «0,0». Se muestran 3 cifras significativas: 1,23e-13.
        return f'{valor:.3g}'.replace('.', ',')
    return f'{valor:.{decimales}f}'.replace('.', ',')


def texto_angulo(grados):
    """Ángulo con su signo siempre visible: +12,5 / -7 / 0."""
    if grados == 0:
        return '0'
    return f'{grados:+g}'.replace('.', ',')


def leer_angulo(texto):
    """Lee un ángulo escrito a mano ('+12,5', '-7', '−7°'...). None si no es un número."""
    texto = str(texto).replace('−', '-').replace('°', '').replace(',', '.').replace(' ', '')
    try:
        valor = float(texto)
    except ValueError:
        return None
    return valor if math.isfinite(valor) else None


def leer_numero(texto, defecto, minimo, maximo):
    try:
        valor = float(str(texto).replace(',', '.'))
    except ValueError:
        return defecto
    return min(max(valor, minimo), maximo)


# -----------------------------------------------------------------------------
# Exportar a Excel (.xlsx) sin bibliotecas externas
# -----------------------------------------------------------------------------
# Cada hoja es (nombre, columnas, filas). Columnas: (título, decimales, ancho);
# decimales = None -> texto; 0 -> entero; 1-4 -> número con esos decimales.

def sin_caracteres_control(texto):
    """Quita caracteres invisibles de control (no válidos en un .xlsx)."""
    return ''.join(ch for ch in str(texto) if ch in '\t\n' or ord(ch) >= 32)


def texto_seguro_csv(texto):
    """Evita la "inyección de fórmulas": si un texto (por ejemplo, la fecha de la cabecera
    FITS o un nombre de archivo) empieza por = + - @, Excel lo podría ejecutar como fórmula
    al abrir el CSV. Se le pone un apóstrofo delante para que sea solo texto."""
    texto = sin_caracteres_control(texto)
    return "'" + texto if texto[:1] in ('=', '+', '-', '@') else texto


def _letra_columna(indice):
    """0 -> A, 1 -> B, ..., 26 -> AA (nombres de columna de Excel)."""
    letras = ''
    indice += 1
    while indice:
        indice, resto = divmod(indice - 1, 26)
        letras = chr(65 + resto) + letras
    return letras


def _hoja_xml(columnas, filas):
    estilo_decimales = {0: 2, 1: 3, 2: 4, 3: 5, 4: 6, 5: 7, 6: 8}  # índice del estilo según decimales

    def celda_texto(ref, texto, estilo=0):
        texto = escape(sin_caracteres_control(texto))
        return f'<c r="{ref}" t="inlineStr" s="{estilo}"><is><t xml:space="preserve">{texto}</t></is></c>'

    filas_xml = ['<row r="1" ht="22" customHeight="1">' + ''.join(
        celda_texto(f'{_letra_columna(j)}1', titulo, 1) for j, (titulo, _d, _a) in enumerate(columnas))
        + '</row>']
    for i, fila in enumerate(filas, start=2):
        celdas = []
        for j, (valor, (_t, decimales, _a)) in enumerate(zip(fila, columnas)):
            ref = f'{_letra_columna(j)}{i}'
            if valor is None or valor == '':
                continue
            if decimales is None:
                celdas.append(celda_texto(ref, valor))
            elif isinstance(valor, (int, float, np.integer, np.floating)) and np.isfinite(valor):
                numero_txt = str(int(valor)) if decimales == 0 else repr(float(valor))
                celdas.append(f'<c r="{ref}" s="{estilo_decimales[decimales]}"><v>{numero_txt}</v></c>')
        filas_xml.append(f'<row r="{i}">' + ''.join(celdas) + '</row>')
    anchos = ''.join(f'<col min="{j + 1}" max="{j + 1}" width="{ancho}" customWidth="1"/>'
                     for j, (_t, _d, ancho) in enumerate(columnas))
    rango = f'A1:{_letra_columna(len(columnas) - 1)}{len(filas) + 1}'
    ns = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'
    return ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
            f'<worksheet xmlns="{ns}"><sheetViews><sheetView workbookViewId="0">'
            '<pane ySplit="1" topLeftCell="A2" activePane="bottomLeft" state="frozen"/>'
            '</sheetView></sheetViews>'
            f'<cols>{anchos}</cols><sheetData>{"".join(filas_xml)}</sheetData>'
            f'<autoFilter ref="{rango}"/></worksheet>'), rango


def escribir_xlsx(ruta, hojas):
    """Crea un libro de Excel (.xlsx) con una o varias hojas.

    Un .xlsx es en realidad un archivo ZIP con varios XML dentro. Aquí se
    escriben los mínimos: las hojas con los datos, los estilos (cabecera azul
    oscuro con letra blanca y número de decimales) y los archivos que los
    enlazan. Los números se guardan como números de verdad, así que Excel los
    entiende igual en cualquier idioma (con coma o con punto decimal).
    Cada hoja lleva los botones de filtro en la fila de títulos.
    """
    ns = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'
    ns_r = 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'
    ns_paquete = 'http://schemas.openxmlformats.org/package/2006/relationships'
    cabecera_xml = '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'

    hojas_xml, nombres, rangos = [], [], []
    for nombre, columnas, filas in hojas:
        xml_hoja, rango = _hoja_xml(columnas, filas)
        hojas_xml.append(xml_hoja)
        for prohibido in '[]:*?/\\':
            nombre = nombre.replace(prohibido, ' ')
        nombres.append(nombre[:31])
        rangos.append(rango)

    estilos = (cabecera_xml + f'<styleSheet xmlns="{ns}">'
               '<numFmts count="6"><numFmt numFmtId="164" formatCode="0.0"/>'
               '<numFmt numFmtId="165" formatCode="0.00"/><numFmt numFmtId="166" formatCode="0.000"/>'
               '<numFmt numFmtId="167" formatCode="0.0000"/><numFmt numFmtId="168" formatCode="0.00000"/>'
               '<numFmt numFmtId="169" formatCode="0.000000"/></numFmts>'
               '<fonts count="2"><font><sz val="11"/><name val="Calibri"/></font>'
               '<font><b/><sz val="11"/><color rgb="FFFFFFFF"/><name val="Calibri"/></font></fonts>'
               '<fills count="3"><fill><patternFill patternType="none"/></fill>'
               '<fill><patternFill patternType="gray125"/></fill>'
               '<fill><patternFill patternType="solid"><fgColor rgb="FF203864"/><bgColor indexed="64"/>'
               '</patternFill></fill></fills>'
               '<borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders>'
               '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>'
               '<cellXfs count="9">'
               '<xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>'
               '<xf numFmtId="0" fontId="1" fillId="2" borderId="0" xfId="0" applyFont="1" applyFill="1" '
               'applyAlignment="1"><alignment vertical="center"/></xf>'
               '<xf numFmtId="1" fontId="0" fillId="0" borderId="0" xfId="0" applyNumberFormat="1"/>'
               '<xf numFmtId="164" fontId="0" fillId="0" borderId="0" xfId="0" applyNumberFormat="1"/>'
               '<xf numFmtId="165" fontId="0" fillId="0" borderId="0" xfId="0" applyNumberFormat="1"/>'
               '<xf numFmtId="166" fontId="0" fillId="0" borderId="0" xfId="0" applyNumberFormat="1"/>'
               '<xf numFmtId="167" fontId="0" fillId="0" borderId="0" xfId="0" applyNumberFormat="1"/>'
               '<xf numFmtId="168" fontId="0" fillId="0" borderId="0" xfId="0" applyNumberFormat="1"/>'
               '<xf numFmtId="169" fontId="0" fillId="0" borderId="0" xfId="0" applyNumberFormat="1"/>'
               '</cellXfs>'
               '<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>'
               '</styleSheet>')

    def nombre_formula(nombre):
        return "'" + nombre.replace("'", "''") + "'"

    def absoluto(rango):  # A1:K19 -> $A$1:$K$19
        partes = []
        for ref in rango.split(':'):
            letras = ''.join(ch for ch in ref if ch.isalpha())
            partes.append(f'${letras}${ref[len(letras):]}')
        return ':'.join(partes)

    hojas_libro = ''.join(f'<sheet name="{escape(n, {chr(34): "&quot;"})}" sheetId="{i + 1}" r:id="rId{i + 1}"/>'
                          for i, n in enumerate(nombres))
    nombres_definidos = ''.join(
        f'<definedName name="_xlnm._FilterDatabase" localSheetId="{i}" hidden="1">'
        f'{escape(nombre_formula(n))}!{absoluto(r)}</definedName>'
        for i, (n, r) in enumerate(zip(nombres, rangos)))
    libro = (cabecera_xml + f'<workbook xmlns="{ns}" xmlns:r="{ns_r}"><sheets>{hojas_libro}</sheets>'
             f'<definedNames>{nombres_definidos}</definedNames></workbook>')
    enlaces_libro = (cabecera_xml + f'<Relationships xmlns="{ns_paquete}">'
                     + ''.join(f'<Relationship Id="rId{i + 1}" Type="{ns_r}/worksheet" '
                               f'Target="worksheets/sheet{i + 1}.xml"/>' for i in range(len(hojas_xml)))
                     + f'<Relationship Id="rId{len(hojas_xml) + 1}" Type="{ns_r}/styles" Target="styles.xml"/>'
                     '</Relationships>')
    enlaces = (cabecera_xml + f'<Relationships xmlns="{ns_paquete}">'
               f'<Relationship Id="rId1" Type="{ns_r}/officeDocument" Target="xl/workbook.xml"/>'
               '</Relationships>')
    tipo = 'application/vnd.openxmlformats-officedocument.spreadsheetml'
    tipos = (cabecera_xml + '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
             '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
             '<Default Extension="xml" ContentType="application/xml"/>'
             f'<Override PartName="/xl/workbook.xml" ContentType="{tipo}.sheet.main+xml"/>'
             + ''.join(f'<Override PartName="/xl/worksheets/sheet{i + 1}.xml" ContentType="{tipo}.worksheet+xml"/>'
                       for i in range(len(hojas_xml)))
             + f'<Override PartName="/xl/styles.xml" ContentType="{tipo}.styles+xml"/>'
             '</Types>')

    with zipfile.ZipFile(ruta, 'w', zipfile.ZIP_DEFLATED) as paquete:
        paquete.writestr('[Content_Types].xml', tipos)
        paquete.writestr('_rels/.rels', enlaces)
        paquete.writestr('xl/workbook.xml', libro)
        paquete.writestr('xl/_rels/workbook.xml.rels', enlaces_libro)
        paquete.writestr('xl/styles.xml', estilos)
        for i, xml_hoja in enumerate(hojas_xml):
            paquete.writestr(f'xl/worksheets/sheet{i + 1}.xml', xml_hoja)


def _ruta_ajustes():
    # Windows: %APPDATA%\Emillium · Linux y macOS: ~/.config/Emillium (igual que los catálogos)
    base = os.environ.get('APPDATA') or os.path.join(os.path.expanduser('~'), '.config')
    return os.path.join(base, 'Emillium', 'ajustes.json')


def cargar_ajustes():
    """Ajustes guardados la vez anterior (idioma, opciones, última carpeta). Si no hay
    o el archivo está dañado, se empieza con los valores por defecto."""
    try:
        with open(_ruta_ajustes(), encoding='utf-8') as archivo:
            datos = json.load(archivo)
        return datos if isinstance(datos, dict) else {}
    except (OSError, ValueError):
        return {}


def guardar_ajustes(ajustes):
    try:
        os.makedirs(os.path.dirname(_ruta_ajustes()), exist_ok=True)
        with open(_ruta_ajustes(), 'w', encoding='utf-8') as archivo:
            json.dump(ajustes, archivo, ensure_ascii=False, indent=1)
    except OSError:
        pass  # si no se puede guardar, no pasa nada


class Ayuda:
    """Recuadro de ayuda que aparece al dejar el ratón encima de un control."""

    def __init__(self, widget, texto):
        self.widget = widget
        self.texto = texto
        self.ventana = None
        self.espera = None
        widget.bind('<Enter>', self._programar, add='+')
        widget.bind('<Leave>', self._ocultar, add='+')
        widget.bind('<ButtonPress>', self._ocultar, add='+')

    def _programar(self, _evento):
        self.espera = self.widget.after(500, self._mostrar)

    def _mostrar(self):
        x = self.widget.winfo_rootx() + 10
        y = self.widget.winfo_rooty() + self.widget.winfo_height() + 4
        self.ventana = tk.Toplevel(self.widget)
        self.ventana.wm_overrideredirect(True)
        self.ventana.wm_geometry(f'+{x}+{y}')
        tk.Label(self.ventana, text=self.texto, background=CAMPO, foreground=TEXTO,
                 highlightbackground=COLOR_ACENTO, highlightthickness=1, borderwidth=0,
                 wraplength=320, justify=tk.LEFT, padx=8, pady=6).pack()

    def _ocultar(self, _evento=None):
        if self.espera:
            self.widget.after_cancel(self.espera)
            self.espera = None
        if self.ventana:
            self.ventana.destroy()
            self.ventana = None


class EmilliumApp:
    def __init__(self, root):
        self.root = root
        self.root.title(TITULO_VENTANA)
        self.root.configure(background=FONDO)
        ancho_pantalla, alto_pantalla = root.winfo_screenwidth(), root.winfo_screenheight()
        self.root.geometry(f'{min(1500, int(ancho_pantalla * 0.9))}x{min(950, int(alto_pantalla * 0.85))}')
        self.root.minsize(900, 600)
        self.ajustes = cargar_ajustes()   # idioma, opciones y última carpeta de la vez anterior
        self.idioma = self.ajustes.get('idioma') if self.ajustes.get('idioma') in TEXTOS else IDIOMA_POR_DEFECTO
        self.serie = []                   # fotos abiertas juntas (para ◀ ▶)
        self.alineaciones = {}            # (referencia, otra) -> cómo se alinea una foto con otra
        self.ruta_anterior = None         # foto de la que se viene (para «medir como en la anterior»)

        # Imagen
        self.ruta = None
        self.original = None     # imagen tal como está en el archivo
        self.imagen = None       # imagen después de espejo y giro (la que se ve)
        self.nota = ''           # aviso al leer el FITS (clave de TEXTOS)
        # Lo que elige el usuario: espejos y un ángulo de giro (con signo, respecto a la original).
        self.mirall_h = False
        self.mirall_v = False
        self.rotacion = 0.0      # + = sentido antihorario, - = sentido horario
        # Estado interno equivalente que usa transformar_imagen: espejo horizontal y luego giro
        self.angulo = 0.0
        self.espejo = False
        self.nivel_bajo = 0.0    # niveles para convertir a gris
        self.nivel_alto = 1.0
        self.fotos = {}          # ruta -> {'n': número de foto, 'archivo', 'hora', 'ancho', 'alto'}

        # Vista (zoom y desplazamiento)
        self.zoom = 1.0          # píxeles de pantalla por píxel de imagen
        self.centro_x = 0.0      # punto de la imagen que está en el centro de la ventana
        self.centro_y = 0.0
        self.foto = None
        self._dibujo_pendiente = False
        self._arrastre = None
        self.modo_brillo = int(self.ajustes.get('brillo', 1)) % 3  # 0 suave, 1 normal, 2 resaltar débiles

        # Mediciones
        self.mediciones = []
        self._contador = 0

        # Opciones que se conservan al cambiar de idioma
        self.var_angulo = tk.StringVar(value='0')
        self.var_angulo.trace_add('write', lambda *_: self._actualizar_sentido())
        self.var_busqueda = tk.StringVar(value=str(self.ajustes.get('busqueda', '15')))
        self.var_apertura = tk.StringVar(value=str(self.ajustes.get('apertura', '5')))
        self.var_umbral = tk.StringVar(value=str(self.ajustes.get('umbral', '3')))
        self.var_marcas = tk.BooleanVar(value=True)
        self.var_estrella = tk.StringVar(value='A')
        self.var_avanzar = tk.BooleanVar(value=True)
        # Fotometría diferencial: estrellas de comparación ("B=12,35 C=11,90/0,65 D"),
        # coeficiente de color T e índices de color de las estrellas objetivo
        self.var_comps = tk.StringVar(value='')
        self.var_color_t = tk.StringVar(value='0')
        self.var_colores = tk.StringVar(value='')
        for var in (self.var_comps, self.var_color_t, self.var_colores):
            var.trace_add('write', lambda *_: self._actualizar_fotometria())
        # Fotometría de apertura: radio del círculo y del anillo del cielo (píxeles)
        self.var_r_fot = tk.StringVar(value=str(self.ajustes.get('r_fot', '5')))
        self.var_r_int = tk.StringVar(value=str(self.ajustes.get('r_int', '8')))
        self.var_r_ext = tk.StringVar(value=str(self.ajustes.get('r_ext', '14')))
        self.avanzadas_visibles = False
        # Al cambiar una opción avanzada se vuelve a calcular la medida (la seleccionada o la última)
        self._recalculo_pendiente = None
        for var in (self.var_busqueda, self.var_apertura, self.var_umbral, self.var_r_fot, self.var_r_int, self.var_r_ext):
            var.trace_add('write', lambda *_: self._programar_recalculo())
        # Herramientas de la versión 1.5
        self.fuentes = None         # estrellas detectadas en la foto abierta (coordenadas de la foto original)
        self.revision = None        # resultado de «Revisar la sèrie»
        self._descargas = {}        # descargas de catálogos en marcha (las pide el usuario)
        self.historial = []         # operaciones hechas (se guarda en la sesión y en el Excel)
        self.cabecera = {}          # cabecera FITS de la foto abierta
        self._parpadeo = None
        self._guardado = False
        self.var_aperturas = tk.BooleanVar(value=bool(self.ajustes.get('aperturas', False)))
        self.var_corregir = tk.BooleanVar(value=bool(self.ajustes.get('corregir_contaminacion', False)))
        self.var_corregir.trace_add('write', lambda *_: self._cambio_correccion())
        self.root.report_callback_exception = self._error_inesperado

        # Animaciones
        self.escena = None          # cielo animado de la pantalla de bienvenida
        self.cielo_cabecera = None
        self.logo_cabecera = None
        self._t_animacion = time.perf_counter()
        self._teclas = []           # últimas teclas pulsadas (para los easter eggs)
        self._codigo_pendiente = None
        self.capa = None            # (lienzo, efecto) que tapa la foto un momento
        self.pulsar = None
        self.pulsar_activo = False
        self.lente = None           # GARGANTUA / INTERSTELLAR en marcha
        self.lente_theta = 0.0
        self.ventana_42 = None

        self._estilos()
        self._construir_interfaz()
        self.root.bind('<Control-o>', lambda e: self.abrir_fits())
        self.root.bind('<Control-s>', lambda e: self.guardar_csv())
        self.root.bind('<Control-z>', lambda e: self.borrar_ultima())
        self.root.bind('<Prior>', lambda e: self.ir_a_foto(-1))   # RePág: foto anterior de la serie
        self.root.bind('<Next>', lambda e: self.ir_a_foto(+1))    # AvPág: foto siguiente
        self.root.bind('<KeyPress>', self._tecla_secreta, add='+')
        self.root.protocol('WM_DELETE_WINDOW', self._al_cerrar)
        for rueda in ('<MouseWheel>', '<Button-4>', '<Button-5>'):
            self.root.bind_all(rueda, self._rueda_panel, add='+')
        self.root.after(100, self._animar)

    def _guardar_ajustes(self):
        self.ajustes.update(idioma=self.idioma, brillo=self.modo_brillo, busqueda=self.var_busqueda.get(),
                            apertura=self.var_apertura.get(), umbral=self.var_umbral.get(),
                            r_fot=self.var_r_fot.get(), r_int=self.var_r_int.get(), r_ext=self.var_r_ext.get())
        guardar_ajustes(self.ajustes)

    def _al_cerrar(self):
        self._guardar_ajustes()
        self.root.destroy()

    def t(self, clave, **datos):
        """Texto traducido al idioma actual."""
        texto = TEXTOS[self.idioma][clave]
        return texto.format(**datos) if datos else texto

    # ------------------------------------------------------------------ interfaz
    def _estilos(self):
        from tkinter import font
        for nombre in ('TkDefaultFont', 'TkTextFont', 'TkHeadingFont'):
            try:
                font.nametofont(nombre).configure(family='Segoe UI', size=10)
            except tk.TclError:
                pass
        estilo = ttk.Style()
        estilo.theme_use('clam')  # el tema 'clam' deja cambiar todos los colores
        estilo.configure('.', background=FONDO, foreground=TEXTO, fieldbackground=CAMPO,
                         bordercolor=BORDE, lightcolor=FONDO, darkcolor=FONDO, troughcolor=FONDO,
                         focuscolor=COLOR_ACENTO, selectbackground=PRINCIPAL, selectforeground=TEXTO,
                         insertcolor=TEXTO, arrowcolor=TEXTO, font=('Segoe UI', 10))
        estilo.configure('TButton', background=BOTON, foreground=TEXTO, bordercolor=BORDE,
                         lightcolor=BOTON, darkcolor=BOTON, padding=(8, 5), relief=tk.FLAT)
        estilo.map('TButton',
                   background=[('pressed', BOTON_PULSADO), ('active', BOTON_ENCIMA)],
                   lightcolor=[('pressed', BOTON_PULSADO), ('active', BOTON_ENCIMA)],
                   darkcolor=[('pressed', BOTON_PULSADO), ('active', BOTON_ENCIMA)],
                   bordercolor=[('active', COLOR_ACENTO)])
        estilo.configure('Grande.TButton', font=('Segoe UI', 11, 'bold'), padding=(10, 9),
                         background=PRINCIPAL, lightcolor=PRINCIPAL, darkcolor=PRINCIPAL,
                         bordercolor=COLOR_ACENTO, foreground='#ffffff')
        estilo.map('Grande.TButton',
                   background=[('pressed', BOTON_PULSADO), ('active', PRINCIPAL_ENCIMA)],
                   lightcolor=[('pressed', BOTON_PULSADO), ('active', PRINCIPAL_ENCIMA)],
                   darkcolor=[('pressed', BOTON_PULSADO), ('active', PRINCIPAL_ENCIMA)])
        estilo.configure('TLabelframe', background=FONDO, bordercolor=BORDE, lightcolor=FONDO, darkcolor=FONDO)
        estilo.configure('TLabelframe.Label', background=FONDO, foreground=COLOR_ACENTO,
                         font=('Segoe UI Semibold', 11))
        estilo.configure('Gris.TLabel', foreground=TEXTO_SUAVE)
        estilo.configure('Resultado.TLabel', font=('Consolas', 17, 'bold'), foreground=DORADO)
        estilo.configure('Orientacion.TLabel', foreground=COLOR_ACENTO, font=('Segoe UI', 10, 'bold'))
        estilo.configure('Letra.TCombobox', font=('Segoe UI', 12, 'bold'))
        for campo in ('TEntry', 'TSpinbox', 'TCombobox'):
            estilo.configure(campo, fieldbackground=CAMPO, foreground=TEXTO, background=BOTON,
                             bordercolor=BORDE, lightcolor=CAMPO, darkcolor=CAMPO, arrowcolor=TEXTO)
            estilo.map(campo, bordercolor=[('focus', COLOR_ACENTO)],
                       fieldbackground=[('readonly', CAMPO)], foreground=[('readonly', TEXTO)],
                       selectbackground=[('readonly', CAMPO)], selectforeground=[('readonly', TEXTO)])
        self.root.option_add('*TCombobox*Listbox.background', CAMPO)
        self.root.option_add('*TCombobox*Listbox.foreground', TEXTO)
        self.root.option_add('*TCombobox*Listbox.selectBackground', PRINCIPAL)
        self.root.option_add('*TCombobox*Listbox.selectForeground', TEXTO)
        estilo.configure('TCheckbutton', background=FONDO, foreground=TEXTO,
                         indicatorbackground=CAMPO, indicatorforeground=DORADO)
        estilo.map('TCheckbutton', background=[('active', FONDO)],
                   indicatorbackground=[('selected', CAMPO)])
        estilo.configure('Treeview', background=CAMPO, fieldbackground=CAMPO, foreground=TEXTO,
                         bordercolor=BORDE, lightcolor=CAMPO, darkcolor=CAMPO, rowheight=24)
        estilo.map('Treeview', background=[('selected', PRINCIPAL)], foreground=[('selected', '#ffffff')])
        estilo.configure('Treeview.Heading', background=BOTON, foreground=COLOR_ACENTO,
                         bordercolor=BORDE, lightcolor=BOTON, darkcolor=BOTON,
                         font=('Segoe UI', 10, 'bold'), relief=tk.FLAT)
        estilo.map('Treeview.Heading', background=[('active', BOTON_ENCIMA)])
        estilo.configure('Vertical.TScrollbar', background=BOTON, troughcolor=FONDO, bordercolor=FONDO,
                         lightcolor=BOTON, darkcolor=BOTON, arrowcolor=TEXTO_SUAVE)
        estilo.map('Vertical.TScrollbar', background=[('active', BOTON_ENCIMA)])
        estilo.configure('Peligro.TButton', background='#7a1f2b', foreground='#ffe9ec', bordercolor='#b8394b',
                         lightcolor='#7a1f2b', darkcolor='#7a1f2b', font=('Segoe UI', 8, 'bold'), padding=(6, 1))
        estilo.map('Peligro.TButton',
                   background=[('pressed', '#5c1520'), ('active', '#a32838')],
                   lightcolor=[('pressed', '#5c1520'), ('active', '#a32838')],
                   darkcolor=[('pressed', '#5c1520'), ('active', '#a32838')],
                   bordercolor=[('active', '#ff8a9a')])
        estilo.configure('Pequeno.TButton', font=('Segoe UI', 9), padding=(6, 2))
        estilo.configure('Emillium.Horizontal.TProgressbar', background=PRINCIPAL, troughcolor=CAMPO,
                         bordercolor=BORDE, lightcolor=PRINCIPAL_ENCIMA, darkcolor=PRINCIPAL)

    def _cabecera(self):
        """Franja de arriba con cielo estrellado y el nombre del programa."""
        self.lienzo_cabecera = tk.Canvas(self.root, height=58, background='#060a18', highlightthickness=0)
        self.lienzo_cabecera.pack(side=tk.TOP, fill=tk.X)
        self.lienzo_cabecera.bind('<Configure>', self._dibujar_cabecera)

    def _dibujar_cabecera(self, evento):
        c = self.lienzo_cabecera
        c.delete('all')
        ancho, alto = evento.width, evento.height
        # degradado suave de azul noche
        for i in range(0, ancho, 6):
            t = i / max(ancho, 1)
            color = '#%02x%02x%02x' % (int(6 + 14 * t), int(10 + 14 * t), int(24 + 30 * t))
            c.create_rectangle(i, 0, i + 6, alto, fill=color, outline='')
        self.logo_cabecera = Logo(c, 40, alto / 2, 20)
        x1 = self.logo_cabecera.caja[2]
        subtitulo = c.create_text(x1 + 42, alto / 2 + 2, anchor=tk.W, fill=TEXTO_SUAVE,
                                  font=('Segoe UI', 10), text=self.t('subtitulo'))
        # estrellas solo a la derecha del texto, para que no lo tapen
        self.cielo_cabecera = CampoEstelar(c, ancho, alto, semilla=3, densidad=700, velocidad=7,
                                           x_inicio=c.bbox(subtitulo)[2] + 30)
        c.create_line(0, alto - 1, ancho, alto - 1, fill=BORDE)
        if self.pulsar_activo:  # el púlsar de la supernova se queda para siempre
            self.pulsar = Pulsar(c, ancho - 26, alto / 2)

    def _animar(self):
        """Bucle de animación (unas 30 veces por segundo)."""
        ahora = time.perf_counter()
        dt = min(ahora - self._t_animacion, 0.1)
        self._t_animacion = ahora
        try:
            if self.cielo_cabecera is not None:
                self.cielo_cabecera.mover(dt)
                self.logo_cabecera.animar(ahora, 1.5)
            escena = self.escena
            if escena is not None and self.imagen is None:
                escena['cielo'].mover(dt)
                escena['galaxia'].mover(dt)
                # easter eggs del cielo
                for cual in ('nebulosa', 'supernova', 'pilares', 'voyager', 'iss', 'odisea'):
                    if escena.get(cual) and not escena[cual].paso(dt):
                        escena[cual] = None
                if escena.get('agujero'):
                    escena['cielo'].orbitar(escena['agujero'], dt)
                    if not escena['agujero'].paso(dt):
                        escena['agujero'] = None
                escena['logo'].animar(ahora, 7)
                azar = escena['azar']
                if ahora >= escena['proxima_fugaz']:
                    escena['fugaces'].append(EstrellaFugaz(self.lienzo, escena['ancho'], escena['alto'], azar))
                    escena['proxima_fugaz'] = ahora + azar.uniform(1.8, 5.0)
                if ahora >= escena['proximo_satelite']:
                    escena['fugaces'].append(Satelite(self.lienzo, escena['ancho'], escena['alto'], azar))
                    escena['proximo_satelite'] = ahora + azar.uniform(18, 35)
                escena['fugaces'] = [f for f in escena['fugaces'] if f.paso(dt)]
        except tk.TclError:
            pass  # algún dibujo se ha borrado (cambio de idioma, ventana cerrándose...)
        for nombre, paso in (('pulsar', lambda: self.pulsar.paso(dt) if self.pulsar else None),
                             ('capa', lambda: self.capa and (self.capa[1].paso(dt) or self._cerrar_capa())),
                             ('lente', lambda: self.lente is not None and self._paso_lente(dt))):
            try:
                paso()
            except tk.TclError:
                if nombre == 'capa':
                    self.capa = None
                elif nombre == 'lente':
                    self.lente, self.lente_theta = None, 0.0
        try:
            self.root.after(25, self._animar)  # en Windows esto da unos 30 fotogramas por segundo
        except tk.TclError:
            pass  # la ventana ya se ha cerrado

    def _construir_interfaz(self):
        self._cabecera()

        # Barra de estado (abajo del todo)
        barra = ttk.Frame(self.root, padding=(8, 3))
        barra.pack(side=tk.BOTTOM, fill=tk.X)
        self.lbl_cursor = ttk.Label(barra, text='', width=44)
        self.lbl_cursor.pack(side=tk.LEFT)
        ttk.Label(barra, style='Gris.TLabel', text=self.t('barra_ayuda')).pack(side=tk.RIGHT)

        # Panel lateral derecho con los pasos. Va dentro de un lienzo con barra de
        # desplazamiento (que solo aparece si no cabe todo) para poder subir y bajar.
        self.contenedor_panel = ttk.Frame(self.root, width=345)
        self.contenedor_panel.pack(side=tk.RIGHT, fill=tk.Y)
        self.contenedor_panel.pack_propagate(False)
        self.lienzo_panel = tk.Canvas(self.contenedor_panel, background=FONDO, highlightthickness=0,
                                      yscrollincrement=24)
        self.barra_panel = ttk.Scrollbar(self.contenedor_panel, orient=tk.VERTICAL, command=self.lienzo_panel.yview)
        self.lienzo_panel.configure(yscrollcommand=self.barra_panel.set)
        self.lienzo_panel.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        panel = ttk.Frame(self.lienzo_panel, padding=(4, 8, 8, 4))
        self._ventana_panel = self.lienzo_panel.create_window(0, 0, window=panel, anchor=tk.NW)
        panel.bind('<Configure>', lambda e: self._ajustar_panel())
        self.lienzo_panel.bind('<Configure>', lambda e: (
            self.lienzo_panel.itemconfig(self._ventana_panel, width=e.width), self._ajustar_panel()))
        self._panel_foto(panel)
        self._panel_orientacion(panel)
        self._panel_medir(panel)
        self._panel_fotometria(panel)
        self._panel_analisis(panel)

        # Tabla de mediciones (debajo de la imagen) y barra de vista (encima)
        self._panel_tabla()
        self._panel_vista()
        self._actualizar_fotometria()   # cabecera de la columna de magnitud y avisos de la fotometría
        # Lo que se esconde en el modo sencillo
        self._widgets_avanzados = [(self.caja_analisis, {'fill': tk.X, 'pady': (0, 8)}),
                                   (self.boton_fotometria, {'fill': tk.X, 'pady': (6, 0)}),
                                   (self.boton_avanzadas, {'fill': tk.X, 'pady': (8, 0), 'after': self.lbl_mensaje})]
        self.root.after_idle(self._aplicar_modo)

        # Lienzo de la imagen (ocupa el resto)
        self.lienzo = tk.Canvas(self.root, background=ESPACIO, highlightthickness=1,
                                highlightbackground=BORDE, cursor='crosshair')
        self.lienzo.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=(8, 4), pady=(8, 4))
        self.lienzo.bind('<Configure>', lambda e: self._programar_dibujo())
        self.lienzo.bind('<Button-1>', self._clic_izquierdo)
        for boton in (2, 3):
            self.lienzo.bind(f'<ButtonPress-{boton}>', self._empezar_arrastre)
            self.lienzo.bind(f'<B{boton}-Motion>', self._arrastrar)
            self.lienzo.bind(f'<ButtonRelease-{boton}>', self._terminar_arrastre)
        self.lienzo.bind('<MouseWheel>', self._rueda)           # Windows y macOS
        self.lienzo.bind('<Button-4>', lambda e: self._zoom_en(e.x, e.y, 1.25))  # Linux
        self.lienzo.bind('<Button-5>', lambda e: self._zoom_en(e.x, e.y, 1 / 1.25))
        self.lienzo.bind('<Motion>', self._mover_cursor)
        self.lienzo.bind('<Leave>', lambda e: self.lbl_cursor.config(text=''))

    # ------------------------------------------------------------------ panel desplazable
    def _panel_desborda(self):
        contenido = self.lienzo_panel.bbox(self._ventana_panel)
        return contenido is not None and contenido[3] > self.lienzo_panel.winfo_height() + 1

    def _ajustar_panel(self):
        """Actualiza la zona desplazable y muestra la barra solo si hace falta."""
        try:
            self.lienzo_panel.configure(scrollregion=self.lienzo_panel.bbox(self._ventana_panel))
            if self._panel_desborda():
                if not self.barra_panel.winfo_ismapped():
                    self.barra_panel.pack(side=tk.RIGHT, fill=tk.Y, before=self.lienzo_panel)
            else:
                if self.barra_panel.winfo_ismapped():
                    self.barra_panel.pack_forget()
                self.lienzo_panel.yview_moveto(0)
        except tk.TclError:
            pass

    def _rueda_panel(self, evento):
        """Rueda del ratón encima del panel derecho: lo sube o lo baja."""
        try:
            widget = str(evento.widget)
            if not widget.startswith(str(self.contenedor_panel)) or not self._panel_desborda():
                return
        except (tk.TclError, AttributeError):
            return
        pasos = -1 if (getattr(evento, 'delta', 0) > 0 or getattr(evento, 'num', 0) == 4) else 1
        self.lienzo_panel.yview_scroll(pasos * 2, 'units')

    # ------------------------------------------------------------------ idioma
    def cambiar_idioma(self, idioma):
        """Cambia el idioma rehaciendo la ventana, sin perder la foto ni las medidas."""
        if idioma == self.idioma:
            return
        self.idioma = idioma
        self.capa = None
        self.lente, self.lente_theta = None, 0.0
        self.escena = None
        self.cielo_cabecera = None
        self.logo_cabecera = None
        for widget in self.root.winfo_children():
            widget.destroy()
        self._construir_interfaz()
        self._mostrar_avanzadas()
        # volver a poner la información de la foto, la orientación, la tabla y el último resultado
        self._mostrar_info_foto()
        self._actualizar_orientacion()
        for medida in self.mediciones:
            self._anadir_a_tabla(medida)
        self._actualizar_contador()
        if self.mediciones:
            self._mostrar_resultado(self.mediciones[-1])
        self._actualizar_botones_serie()
        self._guardar_ajustes()
        self._programar_dibujo()

    # ------------------------------------------------------------------ easter eggs
    def _tecla_secreta(self, evento):
        """Detecta las combinaciones secretas (CODIGOS_SECRETOS) tecleadas seguidas,
        con menos de 1,2 s entre tecla y tecla."""
        try:
            foco = self.root.focus_get()
        except (KeyError, tk.TclError):
            foco = None
        if isinstance(foco, (tk.Entry, ttk.Entry, tk.Spinbox)):
            return  # se está escribiendo un número en una casilla
        letra = (evento.char or '').lower()
        if not letra or not letra.isalnum():
            return
        ahora = time.perf_counter()
        if self._teclas and ahora - self._teclas[-1][1] > 1.2:
            self._teclas = []
        self._teclas = (self._teclas + [(letra, ahora)])[-12:]
        if self._codigo_pendiente:
            self.root.after_cancel(self._codigo_pendiente)
            self._codigo_pendiente = None
        escrito = ''.join(l for l, _t in self._teclas)
        coincidencias = [codigo for codigo in CODIGOS_SECRETOS if escrito.endswith(codigo)]
        if not coincidencias:
            return
        codigo = max(coincidencias, key=len)  # "m42" gana a "42"

        def disparar():
            self._codigo_pendiente = None
            self._teclas = []
            self.lanzar_easter_egg(CODIGOS_SECRETOS[codigo])
        if any(len(otro) > len(codigo) and otro.startswith(codigo) for otro in CODIGOS_SECRETOS):
            # "m1" podría ser el principio de "m16": se espera un momento por si sigue
            self._codigo_pendiente = self.root.after(800, disparar)
        else:
            disparar()

    def lanzar_easter_egg(self, cual):
        if cual in EFECTOS_DE_CIELO:
            escena = self.escena
            if escena is None or self.imagen is not None or escena.get(cual):
                return  # solo en la pantalla de inicio, y uno de cada a la vez
            l, ancho, alto = self.lienzo, escena['ancho'], escena['alto']
            if cual == 'nebulosa':
                escena[cual] = NebulosaOrion(l, ancho, alto, self.t('m42'))
            elif cual == 'agujero':
                escena[cual] = AgujeroNegro(l, ancho, alto)
            elif cual == 'supernova':
                escena[cual] = Supernova(l, ancho, alto, self.t('crab'), self.activar_pulsar)
            elif cual == 'pilares':
                escena[cual] = PilaresCreacion(l, ancho, alto, self.t('m16'), escena['cielo'])
            elif cual == 'voyager':
                escena[cual] = Voyager(l, ancho, alto, self.t('voyager'))
            elif cual == 'iss':
                escena[cual] = EstacionEspacial(l, ancho, alto, self.t('iss'))
            elif cual == 'odisea':
                escena[cual] = Odisea2001(l, ancho, alto, self.t('odisea'))
        elif cual == 'lente':
            if self.lente is None:
                self.lente = {'t': 0.0}
        elif cual == 'respuesta':
            self.mostrar_42()
        else:
            self._mostrar_capa(cual)

    # --- capa que tapa la foto un momento (EARTH, 6EQUJ5, SPOCK)
    def _mostrar_capa(self, cual):
        self._cerrar_capa()
        capa = tk.Canvas(self.root, background=ESPACIO, highlightthickness=0, cursor='hand2')
        capa.place(in_=self.lienzo, x=0, y=0, relwidth=1, relheight=1)
        capa.bind('<Button-1>', lambda e: self._cerrar_capa())  # un clic la cierra
        self.root.update_idletasks()
        ancho, alto = max(capa.winfo_width(), 50), max(capa.winfo_height(), 50)
        if cual == 'tierra':
            efecto = PuntoAzulPalido(capa, ancho, alto, self.t('pale_blue_dot'))
        elif cual == 'wow':
            efecto = SenalWow(capa, ancho, alto, self.t('wow'), self.t('wow_sub'))
        else:
            efecto = SaludoVulcano(capa, ancho, alto, self.t('vulcano'))
        self.capa = (capa, efecto)

    def _cerrar_capa(self):
        if self.capa:
            try:
                self.capa[0].destroy()
            except tk.TclError:
                pass
            self.capa = None

    # --- púlsar en la cabecera (queda después de la supernova)
    def activar_pulsar(self):
        self.pulsar_activo = True
        c = self.lienzo_cabecera
        self.pulsar = Pulsar(c, c.winfo_width() - 26, c.winfo_height() / 2)

    # --- lente gravitatoria (GARGANTUA / INTERSTELLAR)
    def _paso_lente(self, dt):
        lente = self.lente
        lente['t'] += dt
        duracion = 3.6
        fraccion = lente['t'] / duracion
        ancho_v, alto_v = self._tamano_lienzo()
        theta = 0.0 if fraccion >= 1 else 0.17 * min(ancho_v, alto_v) * math.sin(math.pi * fraccion) ** 2
        self.lente_theta = theta
        escena = self.escena
        if self.imagen is not None:
            self._programar_dibujo()  # la foto se vuelve a dibujar deformada
        elif escena is not None:
            cx, cy = escena['ancho'] * 0.5, escena['alto'] * 0.24  # encima del título, para que se vea
            escena['cielo'].lente = (cx, cy, theta) if theta > 0 else None
            if 'sombra' not in lente:
                lente['sombra'] = self.lienzo.create_oval(0, 0, 0, 0, fill='#000000', outline='')
                lente['anillo'] = self.lienzo.create_oval(0, 0, 0, 0, outline='#ffcf8a', width=2)
                lente['halo'] = self.lienzo.create_oval(0, 0, 0, 0, outline='#ff8a3d', width=5)
                for clave in ('halo', 'sombra', 'anillo'):
                    self.lienzo.tag_lower(lente[clave], 'frente')
            for clave, factor in (('sombra', 0.62), ('anillo', 1.0), ('halo', 1.06)):
                r = theta * factor
                self.lienzo.coords(lente[clave], cx - r, cy - r, cx + r, cy + r)
        if fraccion >= 1:
            if escena is not None and 'sombra' in lente:
                for clave in ('sombra', 'anillo', 'halo'):
                    self.lienzo.delete(lente[clave])
            self.lente = None
            self.lente_theta = 0.0

    # --- 42
    def mostrar_42(self):
        if self.ventana_42 is not None and self.ventana_42.winfo_exists():
            return
        v = tk.Toplevel(self.root)
        self.ventana_42 = v
        v.title(self.t('t42_titulo'))
        v.configure(background=FONDO)
        v.transient(self.root)
        v.resizable(False, False)
        marco = ttk.Frame(v, padding=22)
        marco.pack(fill=tk.BOTH, expand=True)
        mensaje = ttk.Label(marco, text=self.t('t42_calc'), wraplength=440, justify=tk.CENTER,
                            font=('Segoe UI', 11))
        mensaje.pack(pady=(0, 14))
        barra = ttk.Progressbar(marco, mode='determinate', maximum=100, length=440, style='Emillium.Horizontal.TProgressbar')
        barra.pack()
        resultado = ttk.Label(marco, text=' ', font=('Consolas', 30, 'bold'), foreground=DORADO)
        resultado.pack(pady=(16, 10))
        boton = ttk.Button(marco, text=self.t('aceptar'), command=v.destroy, state=tk.DISABLED)
        boton.pack()
        v.update_idletasks()
        x = self.root.winfo_rootx() + (self.root.winfo_width() - v.winfo_width()) // 2
        y = self.root.winfo_rooty() + (self.root.winfo_height() - v.winfo_height()) // 2
        v.geometry(f'+{x}+{y}')
        texto_base = self.t('t42_calc')

        def avanzar(paso=0):
            if not v.winfo_exists():
                return
            barra['value'] = min(100, paso * 100 / 75)
            mensaje.config(text=texto_base + '.' * (paso // 6 % 4))
            if paso < 75:
                v.after(40, avanzar, paso + 1)
            else:
                mensaje.config(text=texto_base + '...')
                resultado.config(text=self.t('t42_res'))
                boton.config(state=tk.NORMAL)
                boton.focus_set()
        avanzar()

    # ------------------------------------------------------------------ paneles
    def _panel_foto(self, panel):
        caja = ttk.LabelFrame(panel, text=self.t('paso1'), padding=8)
        caja.pack(fill=tk.X, pady=(0, 8))
        boton = ttk.Button(caja, text=self.t('abrir'), style='Grande.TButton', command=self.abrir_fits)
        boton.pack(fill=tk.X)
        Ayuda(boton, self.t('ayuda_abrir'))
        self.lbl_archivo = ttk.Label(caja, text=self.t('ninguna_foto'), style='Gris.TLabel', wraplength=300)
        self.lbl_archivo.pack(fill=tk.X, pady=(6, 0))
        # Navegación por una serie de fotos: ◀ 3 / 12 ▶ (solo si se abrieron varias)
        self.fila_serie = ttk.Frame(caja)
        self.boton_anterior = ttk.Button(self.fila_serie, text=self.t('foto_anterior'), command=lambda: self.ir_a_foto(-1))
        self.boton_anterior.pack(side=tk.LEFT)
        self.boton_siguiente = ttk.Button(self.fila_serie, text=self.t('foto_siguiente'), command=lambda: self.ir_a_foto(+1))
        self.boton_siguiente.pack(side=tk.RIGHT)
        self.lbl_serie = ttk.Label(self.fila_serie, text='', anchor=tk.CENTER, font=('Segoe UI', 10, 'bold'))
        self.lbl_serie.pack(side=tk.LEFT, expand=True)
        Ayuda(self.boton_anterior, self.t('ayuda_serie'))
        Ayuda(self.boton_siguiente, self.t('ayuda_serie'))
        # Herramientas de serie (solo con varias fotos): medir todas y parpadeo
        self.fila_serie2 = ttk.Frame(caja)
        self.fila_serie2.columnconfigure((0, 1), weight=1)
        b = ttk.Button(self.fila_serie2, text=self.t('boton_medir_serie'), style='Pequeno.TButton', command=self.medir_serie)
        b.grid(row=0, column=0, sticky=tk.EW, padx=(0, 2))
        Ayuda(b, self.t('ayuda_medir_serie'))
        b = ttk.Button(self.fila_serie2, text=self.t('boton_parpadeo'), style='Pequeno.TButton', command=self.parpadeo)
        b.grid(row=0, column=1, sticky=tk.EW, padx=(2, 0))
        Ayuda(b, self.t('ayuda_parpadeo'))
        b = ttk.Button(self.fila_serie2, text=self.t('boton_alinear'), style='Pequeno.TButton', command=self.ventana_alinear)
        b.grid(row=1, column=0, columnspan=2, sticky=tk.EW, pady=(4, 0))
        Ayuda(b, self.t('ayuda_alinear'))
        # Varias imágenes en un mismo archivo (extensiones o planos)
        self.fila_imagen = ttk.Frame(caja)
        self.boton_img_ant = ttk.Button(self.fila_imagen, text='◀', width=3, command=lambda: self.cambiar_imagen(-1))
        self.boton_img_ant.pack(side=tk.LEFT)
        self.boton_img_sig = ttk.Button(self.fila_imagen, text='▶', width=3, command=lambda: self.cambiar_imagen(+1))
        self.boton_img_sig.pack(side=tk.RIGHT)
        self.lbl_imagen = ttk.Label(self.fila_imagen, text='', anchor=tk.CENTER)
        self.lbl_imagen.pack(side=tk.LEFT, expand=True)
        Ayuda(self.lbl_imagen, self.t('ayuda_imagen'))
        self.boton_cerrar_foto = ttk.Button(caja, text=self.t('cerrar_foto'), style='Pequeno.TButton',
                                            command=self.cerrar_foto)
        Ayuda(self.boton_cerrar_foto, self.t('ayuda_cerrar_foto'))
        if self.ruta is not None:  # solo se ve cuando hay una foto abierta
            self.boton_cerrar_foto.pack(anchor=tk.E, pady=(6, 0))

    def _panel_orientacion(self, panel):
        caja = ttk.LabelFrame(panel, text=self.t('paso2'), padding=8)
        caja.pack(fill=tk.X, pady=(0, 8))
        caja.columnconfigure((0, 1), weight=1)

        # Giro: un ángulo con signo, grande y bien visible, con botones de -10° y +10°
        ttk.Label(caja, text=self.t('giro_label')).grid(row=0, column=0, columnspan=2, sticky=tk.W)
        fila = ttk.Frame(caja)
        fila.grid(row=1, column=0, columnspan=2, sticky=tk.EW, pady=(4, 0))
        fila.columnconfigure(1, weight=1)
        ttk.Button(fila, text='−10°', width=6, command=lambda: self.girar(-10)).grid(row=0, column=0, sticky=tk.NS)
        entrada = tk.Entry(fila, textvariable=self.var_angulo, width=8, justify=tk.CENTER,
                           font=('Consolas', 18, 'bold'), background=CAMPO, foreground=DORADO,
                           insertbackground=TEXTO, relief=tk.FLAT, highlightthickness=1,
                           highlightbackground=BORDE, highlightcolor=COLOR_ACENTO)
        entrada.grid(row=0, column=1, sticky=tk.EW, padx=6, ipady=2)
        entrada.bind('<Return>', lambda e: self.aplicar_angulo_escrito())
        entrada.bind('<KP_Enter>', lambda e: self.aplicar_angulo_escrito())
        Ayuda(entrada, self.t('ayuda_giro'))
        ttk.Button(fila, text='+10°', width=6, command=lambda: self.girar(+10)).grid(row=0, column=2, sticky=tk.NS)
        self.lbl_sentido = ttk.Label(caja, text='', font=('Segoe UI', 10, 'bold'))
        self.lbl_sentido.grid(row=2, column=0, columnspan=2, pady=(4, 0))
        ttk.Button(caja, text=self.t('aplicar_giro'), command=self.aplicar_angulo_escrito).grid(
            row=3, column=0, columnspan=2, sticky=tk.EW, pady=(4, 2))
        self._actualizar_sentido()

        ttk.Label(caja, text=self.t('espejo')).grid(row=4, column=0, columnspan=2, sticky=tk.W, pady=(6, 0))
        b = ttk.Button(caja, text=self.t('horizontal'), command=self.espejo_horizontal)
        b.grid(row=5, column=0, sticky=tk.EW, padx=(0, 3), pady=2)
        Ayuda(b, self.t('ayuda_espejo_h'))
        b = ttk.Button(caja, text=self.t('vertical'), command=self.espejo_vertical)
        b.grid(row=5, column=1, sticky=tk.EW, padx=(3, 0), pady=2)
        Ayuda(b, self.t('ayuda_espejo_v'))

        ttk.Button(caja, text=self.t('volver'), command=self.restablecer).grid(
            row=6, column=0, columnspan=2, sticky=tk.EW, pady=(6, 2))
        ttk.Label(caja, text=self.t('ahora_esta'), style='Gris.TLabel').grid(
            row=7, column=0, columnspan=2, sticky=tk.W, pady=(6, 0))
        self.lbl_orientacion = ttk.Label(caja, text=self.t('original'), style='Orientacion.TLabel',
                                         wraplength=300)
        self.lbl_orientacion.grid(row=8, column=0, columnspan=2, sticky=tk.W)

    def _actualizar_sentido(self):
        """Debajo del ángulo escrito, explica qué significa su signo (antes de aplicarlo)."""
        if not hasattr(self, 'lbl_sentido'):
            return
        valor = leer_angulo(self.var_angulo.get())
        try:
            if valor is None:
                self.lbl_sentido.config(text=self.t('angulo_invalido'), foreground=ROJO_SUAVE)
            elif valor > 0:
                self.lbl_sentido.config(text=self.t('sentido_pos'), foreground=COLOR_ACENTO)
            elif valor < 0:
                self.lbl_sentido.config(text=self.t('sentido_neg'), foreground='#ffb36b')
            else:
                self.lbl_sentido.config(text=self.t('sentido_cero'), foreground=TEXTO_SUAVE)
        except tk.TclError:
            pass  # la etiqueta se está rehaciendo (cambio de idioma)

    def _panel_medir(self, panel):
        caja = ttk.LabelFrame(panel, text=self.t('paso3'), padding=8)
        caja.pack(fill=tk.X, pady=(0, 8))
        ttk.Label(caja, text=self.t('medir_texto'), wraplength=300, style='Gris.TLabel').pack(fill=tk.X)

        # Qué estrella se está midiendo (A, B, C...)
        fila = ttk.Frame(caja)
        fila.pack(fill=tk.X, pady=(6, 0))
        etiqueta = ttk.Label(fila, text=self.t('estrella'), font=('Segoe UI', 11, 'bold'))
        etiqueta.pack(side=tk.LEFT)
        selector = ttk.Combobox(fila, textvariable=self.var_estrella, values=LETRAS_ESTRELLA,
                                state='readonly', width=3, font=('Segoe UI', 12, 'bold'), justify=tk.CENTER)
        selector.pack(side=tk.LEFT, padx=6)
        selector.bind('<<ComboboxSelected>>', lambda e: self.lienzo.focus_set())
        ttk.Checkbutton(fila, text=self.t('avanzar'), variable=self.var_avanzar).pack(side=tk.LEFT, padx=(4, 0))
        Ayuda(etiqueta, self.t('ayuda_estrella'))
        Ayuda(selector, self.t('ayuda_estrella'))
        self.fila_estrella = fila
        # Aparece al pasar a otra foto: mide de golpe las estrellas medidas en la foto anterior
        self.boton_como_anterior = ttk.Button(caja, text=self.t('como_anterior'), command=self.medir_como_anterior)
        Ayuda(self.boton_como_anterior, self.t('ayuda_como_anterior'))

        self.lbl_titulo_medida = ttk.Label(caja, text=self.t('ultima'), font=('Segoe UI', 10, 'bold'))
        self.lbl_titulo_medida.pack(anchor=tk.W, pady=(8, 0))
        self.lbl_x = ttk.Label(caja, text='X = -', style='Resultado.TLabel')
        self.lbl_x.pack(anchor=tk.W)
        self.lbl_y = ttk.Label(caja, text='Y = -', style='Resultado.TLabel')
        self.lbl_y.pack(anchor=tk.W)
        self.lbl_cielo = ttk.Label(caja, text='', foreground=COLOR_ACENTO, font=('Consolas', 11, 'bold'), justify=tk.LEFT)
        self.lbl_cielo.pack(anchor=tk.W)
        self.lbl_detalle = ttk.Label(caja, text='', style='Gris.TLabel', wraplength=300, justify=tk.LEFT)
        self.lbl_detalle.pack(fill=tk.X, pady=(2, 0))
        self.lbl_mensaje = ttk.Label(caja, text='', foreground=ROJO_SUAVE, wraplength=300, justify=tk.LEFT)
        self.lbl_mensaje.pack(fill=tk.X)

        # Opciones avanzadas (plegadas)
        self.boton_avanzadas = ttk.Button(caja, text=self.t('mostrar_avanzadas'),
                                          command=self._alternar_avanzadas)
        self.boton_avanzadas.pack(fill=tk.X, pady=(8, 0))
        self.caja_avanzadas = ttk.Frame(caja, padding=(0, 6, 0, 0))
        for fila, (clave, var, desde, hasta, paso, ayuda) in enumerate((
                ('adv_busqueda', self.var_busqueda, 3, 200, 1, 'ayuda_busqueda'),
                ('adv_apertura', self.var_apertura, 1, 50, 0.5, 'ayuda_apertura'),
                ('adv_umbral', self.var_umbral, 0, 50, 0.5, 'ayuda_umbral'),
                ('adv_r_fot', self.var_r_fot, 1, 60, 0.5, 'ayuda_r_fot'),
                ('adv_r_int', self.var_r_int, 2, 80, 0.5, 'ayuda_r_anillo'),
                ('adv_r_ext', self.var_r_ext, 3, 120, 0.5, 'ayuda_r_anillo'))):
            etiqueta = ttk.Label(self.caja_avanzadas, text=self.t(clave))
            etiqueta.grid(row=fila, column=0, sticky=tk.W, pady=2)
            caja_num = ttk.Spinbox(self.caja_avanzadas, textvariable=var, from_=desde, to=hasta,
                                   increment=paso, width=6, command=self._programar_dibujo)
            caja_num.grid(row=fila, column=1, sticky=tk.E, padx=(8, 0), pady=2)
            Ayuda(etiqueta, self.t(ayuda))
            Ayuda(caja_num, self.t(ayuda))
        c = ttk.Checkbutton(self.caja_avanzadas, text=self.t('adv_corregir'), variable=self.var_corregir)
        c.grid(row=10, column=0, columnspan=2, sticky=tk.W, pady=(4, 0))
        Ayuda(c, self.t('ayuda_corregir'))
        self.caja_avanzadas.columnconfigure(0, weight=1)

    def _panel_vista(self):
        """Barra fina encima de la foto: zoom, brillo, marcas e idioma."""
        fila = ttk.Frame(self.root, padding=(8, 8, 4, 0))
        fila.pack(side=tk.TOP, fill=tk.X)
        ttk.Label(fila, text=self.t('zoom')).pack(side=tk.LEFT)
        ttk.Button(fila, text='−', width=3, command=lambda: self._zoom_centro(1 / 1.5)).pack(side=tk.LEFT, padx=(4, 0))
        ttk.Button(fila, text='+', width=3, command=lambda: self._zoom_centro(1.5)).pack(side=tk.LEFT, padx=3)
        ttk.Button(fila, text=self.t('ver_entera'), command=self.ajustar_vista).pack(side=tk.LEFT, padx=3)
        self.lbl_zoom = ttk.Label(fila, text='', style='Gris.TLabel', width=11)
        self.lbl_zoom.pack(side=tk.LEFT, padx=4)

        ttk.Label(fila, text=self.t('brillo')).pack(side=tk.LEFT, padx=(10, 0))
        combo = ttk.Combobox(fila, state='readonly', width=15, values=self.t('modos_brillo'))
        combo.current(self.modo_brillo)
        combo.pack(side=tk.LEFT, padx=6)

        def elegir_brillo(_evento):
            self.modo_brillo = combo.current()
            self._programar_dibujo()
        combo.bind('<<ComboboxSelected>>', elegir_brillo)
        Ayuda(combo, self.t('ayuda_brillo'))
        ttk.Checkbutton(fila, text=self.t('marcas'), variable=self.var_marcas,
                        command=self._programar_dibujo).pack(side=tk.LEFT, padx=12)
        c = ttk.Checkbutton(fila, text=self.t('aperturas'), variable=self.var_aperturas,
                            command=lambda: (self.ajustes.update(aperturas=self.var_aperturas.get()), self._programar_dibujo()))
        c.pack(side=tk.LEFT)
        Ayuda(c, self.t('ayuda_aperturas'))
        # Modo sencillo / avanzado (a la derecha, antes del idioma)
        modos = [self.t('modo_avanzado'), self.t('modo_sencillo')]
        combo_modo = ttk.Combobox(fila, state='readonly', width=13, values=modos)
        combo_modo.current(1 if self.ajustes.get('modo') == 'sencillo' else 0)

        def elegir_modo(_evento):
            self.ajustes['modo'] = 'sencillo' if combo_modo.current() == 1 else 'avanzado'
            guardar_ajustes(self.ajustes)
            self._aplicar_modo()
            self._actualizar_botones_serie()
        combo_modo.bind('<<ComboboxSelected>>', elegir_modo)
        Ayuda(combo_modo, self.t('ayuda_modo'))

        # Selector de idioma (a la derecha)
        idiomas = ttk.Combobox(fila, state='readonly', width=11, values=[nombre for _c, nombre in IDIOMAS])
        idiomas.current([c for c, _n in IDIOMAS].index(self.idioma))
        idiomas.pack(side=tk.RIGHT)
        idiomas.bind('<<ComboboxSelected>>',
                     lambda e: self.root.after_idle(self.cambiar_idioma, IDIOMAS[idiomas.current()][0]))
        ttk.Label(fila, text=self.t('idioma')).pack(side=tk.RIGHT, padx=(0, 4))
        combo_modo.pack(side=tk.RIGHT, padx=(0, 12))

    def _panel_tabla(self):
        marco = ttk.Frame(self.root, padding=(8, 0, 4, 4))
        marco.pack(side=tk.BOTTOM, fill=tk.X)
        cabecera = ttk.Frame(marco)
        cabecera.pack(fill=tk.X, pady=(0, 4))
        ttk.Label(cabecera, text=self.t('mediciones'), font=('Segoe UI', 11, 'bold'),
                  foreground=COLOR_ACENTO).pack(side=tk.LEFT)
        self.lbl_num_medidas = ttk.Label(cabecera, text=self.t('ninguna_todavia'), style='Gris.TLabel')
        self.lbl_num_medidas.pack(side=tk.LEFT, padx=8)
        b = ttk.Button(cabecera, text=self.t('guardar'), command=self.guardar_csv)
        b.pack(side=tk.RIGHT)
        Ayuda(b, self.t('ayuda_guardar'))
        b = ttk.Button(cabecera, text=self.t('boton_pdf'), command=self.guardar_pdf)
        b.pack(side=tk.RIGHT, padx=(0, 6))
        Ayuda(b, self.t('ayuda_pdf'))
        b = ttk.Button(cabecera, text=self.t('mpc_boton'), command=self.abrir_informe_mpc)
        b.pack(side=tk.RIGHT, padx=(0, 6))
        Ayuda(b, self.t('ayuda_mpc'))
        # Solo aparece cuando hay alguna fila seleccionada en la tabla
        self.boton_borrar_sel = ttk.Button(cabecera, text=self.t('borrar_medida'), command=self.borrar_seleccionadas)
        Ayuda(self.boton_borrar_sel, self.t('ayuda_borrar_medida'))
        self._cabecera_tabla = cabecera

        cuerpo = ttk.Frame(marco)
        cuerpo.pack(fill=tk.X)
        self.tabla = ttk.Treeview(cuerpo, columns=[c[0] for c in COLUMNAS_TABLA], show='headings', height=6)
        for clave, clave_titulo, ancho, alineacion in COLUMNAS_TABLA:
            titulo = self.t(clave_titulo) if clave_titulo else f'{clave.upper()} (px)'
            self.tabla.heading(clave, text=titulo)
            self.tabla.column(clave, width=ancho, minwidth=40, anchor=alineacion,
                              stretch=clave in ('archivo', 'orientacion'))
        self.tabla.tag_configure('par', background='#19234a')
        barra = ttk.Scrollbar(cuerpo, orient=tk.VERTICAL, command=self.tabla.yview)
        self.tabla.configure(yscrollcommand=barra.set)
        barra.pack(side=tk.RIGHT, fill=tk.Y)
        self.tabla.pack(side=tk.LEFT, fill=tk.X, expand=True)
        self.tabla.bind('<Double-1>', self._ir_a_medicion)
        self.tabla.bind('<<TreeviewSelect>>', lambda e: self._actualizar_boton_borrar())
        self.tabla.bind('<Delete>', lambda e: self.borrar_seleccionadas())  # tecla Supr
        # Botón para borrar todo: pequeño, rojo y apartado debajo de la tabla (lejos de los demás)
        pie = ttk.Frame(marco)
        pie.pack(fill=tk.X, pady=(3, 0))
        boton = ttk.Button(pie, text=self.t('borrar_todas'), style='Peligro.TButton', command=self.borrar_todas)
        boton.pack(side=tk.RIGHT)
        Ayuda(boton, self.t('ayuda_borrar_todas'))

    def _alternar_avanzadas(self):
        self.avanzadas_visibles = not self.avanzadas_visibles
        self._mostrar_avanzadas()
        if self.avanzadas_visibles:  # bajar el panel para que se vean las opciones recién abiertas
            def bajar():
                self.root.update_idletasks()
                self._ajustar_panel()
                if self._panel_desborda():
                    self.lienzo_panel.yview_moveto(1.0)
            self.root.after_idle(bajar)

    def _mostrar_avanzadas(self):
        if self.avanzadas_visibles:
            self.caja_avanzadas.pack(fill=tk.X)
            self.boton_avanzadas.config(text=self.t('ocultar_avanzadas'))
        else:
            self.caja_avanzadas.pack_forget()
            self.boton_avanzadas.config(text=self.t('mostrar_avanzadas'))

    def _mensaje(self, texto):
        """Aviso en rojo debajo del resultado (vacío para borrarlo)."""
        self.lbl_mensaje.config(text=texto)

    # ------------------------------------------------------------------ archivo
    def abrir_fits(self):
        """Abre una foto o una serie de fotos (se pueden elegir varias a la vez)."""
        rutas = filedialog.askopenfilenames(
            title=self.t('abrir_titulo'), initialdir=self.ajustes.get('carpeta') or None,
            filetypes=[(self.t('tipo_fits'), '*.fits *.fit *.fts *.fz *.gz *.FITS *.FIT *.FTS *.FZ *.GZ'),
                       (self.t('tipo_todos'), '*.*')])
        if rutas:
            self.abrir_serie(list(rutas))

    def abrir_serie(self, rutas):
        """Guarda la serie ordenada por nombre (así suelen quedar en orden de tiempo) y abre la primera."""
        self.serie = sorted(rutas, key=lambda r: os.path.basename(r).lower())
        self.ajustes['carpeta'] = os.path.dirname(self.serie[0])
        guardar_ajustes(self.ajustes)
        self.cargar_fits(self.serie[0])

    def ir_a_foto(self, paso):
        """Botones ◀ / ▶ (o RePág / AvPág): foto anterior o siguiente de la serie,
        manteniendo el giro y el espejo (en una serie todas las fotos suelen ir igual)."""
        self._parar_parpadeo()
        if not self.serie or self.ruta not in self.serie:
            return
        i = self.serie.index(self.ruta) + paso
        if 0 <= i < len(self.serie):
            self.cargar_fits(self.serie[i], mantener_orientacion=True)

    def cargar_fits(self, ruta, mantener_orientacion=False):
        self._parar_parpadeo()
        self.root.config(cursor='watch')
        self.root.update_idletasks()
        try:
            original, cabecera, nota = leer_fits(*separar_ruta(ruta))
        except ErrorFits as error:
            messagebox.showerror(self.t('error_abrir'), self.t(error.clave, **error.datos))
            return
        except Exception as error:
            messagebox.showerror(self.t('error_abrir'), str(error))
            return
        finally:
            self.root.config(cursor='')

        paso = max(1, int(math.sqrt(original.size / 1_000_000)))
        muestra = original[::paso, ::paso]
        muestra = muestra[np.isfinite(muestra)]
        if muestra.size == 0:
            messagebox.showerror(self.t('foto_vacia_t'), self.t('foto_vacia'))
            return

        if self.ruta and self.ruta != ruta:
            self.ruta_anterior = self.ruta  # para «Medir como en la foto anterior»
        if ruta not in self.serie:
            self.serie = [ruta]
        self.ruta = ruta
        self.original = original
        self.nota = nota
        if not mantener_orientacion:
            self.mirall_h = self.mirall_v = False
            self.rotacion = 0.0
        # Número de foto: el orden en que se abrió por primera vez (si se vuelve a abrir, conserva el suyo)
        archivo, indice = separar_ruta(ruta)
        n_imagenes = int(cabecera.get('_N_IMAGENES', 1))
        if ruta not in self.fotos:
            hora, jd = instante_observacion(cabecera)
            nombre = os.path.basename(archivo) + (f' [{indice + 1}/{n_imagenes}]' if n_imagenes > 1 else '')
            self.fotos[ruta] = {'n': len(self.fotos) + 1, 'archivo': nombre, 'hora': hora, 'jd': jd,
                                'exposicion': exposicion_de(cabecera), 'ganancia': ganancia_de(cabecera),
                                'wcs': SolucionWCS.desde_cabecera(cabecera)}
        alto, ancho = original.shape
        self.fotos[ruta].update(ancho=ancho, alto=alto, saturacion=C.nivel_saturacion(cabecera),
                                filtro=str(cabecera.get('FILTER', cabecera.get('FILTNAM1', '')) or '').strip(),
                                objeto=str(cabecera.get('OBJECT', '') or '').strip(),
                                n_imagenes=n_imagenes, indice=indice, airmass=C.masa_de_aire(cabecera),
                                apuntado=C.apuntado_aproximado(cabecera))
        self.cabecera = cabecera
        self.fuentes = None   # las estrellas detectadas eran de la foto anterior
        self._registrar('hist_abrir', nombre=self.fotos[ruta]['archivo'])
        self.var_estrella.set('A')  # en cada foto nueva se empieza por la estrella A
        # Niveles de gris: percentiles 1 % y 99,8 %
        self.nivel_bajo, self.nivel_alto = (float(v) for v in np.percentile(muestra, [1.0, 99.8]))
        if self.nivel_alto <= self.nivel_bajo:
            self.nivel_alto = self.nivel_bajo + 1.0

        self.imagen = original
        if mantener_orientacion and (self.mirall_h or self.mirall_v or self.rotacion):
            self._aplicar_estado()      # se vuelve a aplicar el mismo giro/espejo a la foto nueva
        else:
            self.angulo, self.espejo = 0.0, False
            self.var_angulo.set('0')
        self._actualizar_orientacion()
        self.root.title(TITULO_VENTANA)
        self._mostrar_info_foto()
        if not self.boton_cerrar_foto.winfo_ismapped():
            self.boton_cerrar_foto.pack(anchor=tk.E, pady=(6, 0))
        self._actualizar_botones_serie()
        self._mensaje('')
        if mantener_orientacion and self.zoom > 0:
            self._programar_dibujo()   # misma vista (zoom y posición) que la foto anterior
        else:
            self.ajustar_vista()

    def _mostrar_info_foto(self):
        if self.ruta is None:
            return
        foto = self.fotos[self.ruta]
        texto = self.t('info_foto', nombre=foto['archivo'], n=foto['n'], ancho=foto['ancho'], alto=foto['alto'])
        if foto.get('hora'):
            texto += '\n' + self.t('info_hora', hora=foto['hora'])
        if foto.get('exposicion'):
            texto += self.t('info_exposicion', s=numero(foto['exposicion'], 1).rstrip('0').rstrip(','))
        if foto.get('wcs'):
            texto += '\n' + self.t('info_wcs', escala=numero(foto['wcs'].escala, 2))
            if foto.get('rms_astrometrico') is not None:
                texto += '\n' + self.t('info_rms', rms=numero(foto['rms_astrometrico'], 2),
                                       ajustada=self.t('info_ajustada') if foto.get('wcs_ajustada') else '')
        if foto.get('filtro'):
            texto += '\n' + self.t('info_filtro', f=foto['filtro'])
        if self.nota:
            texto += '\n' + self.t(self.nota)
        self.lbl_archivo.config(text=texto, style='TLabel')

    def _actualizar_botones_serie(self):
        """Muestra ◀ n/N ▶ solo si hay una serie de varias fotos abierta."""
        if len(self.serie) > 1 and self.ruta in self.serie:
            i = self.serie.index(self.ruta)
            self.lbl_serie.config(text=f'{i + 1} / {len(self.serie)}')
            self.boton_anterior.state(['!disabled'] if i > 0 else ['disabled'])
            self.boton_siguiente.state(['!disabled'] if i < len(self.serie) - 1 else ['disabled'])
            if not self.fila_serie.winfo_ismapped():
                self.fila_serie.pack(fill=tk.X, pady=(6, 0), before=self.boton_cerrar_foto)
            if not self.fila_serie2.winfo_ismapped() and self.ajustes.get('modo') != 'sencillo':
                self.fila_serie2.pack(fill=tk.X, pady=(4, 0), before=self.boton_cerrar_foto)
            elif self.fila_serie2.winfo_ismapped() and self.ajustes.get('modo') == 'sencillo':
                self.fila_serie2.pack_forget()
        else:
            if self.fila_serie.winfo_ismapped():
                self.fila_serie.pack_forget()
            if self.fila_serie2.winfo_ismapped():
                self.fila_serie2.pack_forget()
        self._actualizar_selector_imagen()
        hay_anterior = any(m['ruta'] == self.ruta_anterior for m in self.mediciones) and self.ruta is not None
        if hay_anterior and self.ruta_anterior != self.ruta:
            if not self.boton_como_anterior.winfo_ismapped():
                self.boton_como_anterior.pack(fill=tk.X, pady=(6, 0), after=self.fila_estrella)
        elif self.boton_como_anterior.winfo_ismapped():
            self.boton_como_anterior.pack_forget()

    # ------------------------------------------------------------------ giro / espejo
    # El usuario elige: espejo horizontal (sí/no), espejo vertical (sí/no) y un
    # ángulo de giro con signo respecto a la foto original. Primero se aplican
    # los espejos y después el giro (+ = antihorario, - = horario).
    # transformar_imagen solo sabe hacer "espejo horizontal y luego giro", así que
    # se traduce: el espejo vertical es igual que un espejo horizontal + 180°.

    def _aplicar_estado(self):
        self.espejo = self.mirall_h != self.mirall_v
        self.angulo = normalizar_angulo(self.rotacion + (180.0 if self.mirall_v else 0.0))
        self.var_angulo.set(texto_angulo(self.rotacion))
        self._aplicar_transformacion()

    def girar_a(self, grados):
        """Deja la foto girada exactamente `grados` respecto a la original."""
        if self.original is None:
            return
        self.rotacion = normalizar_angulo(grados)
        self._aplicar_estado()

    def girar(self, grados):
        """Botones -10° / +10°: suman al giro actual."""
        self.girar_a(self.rotacion + grados)

    def aplicar_angulo_escrito(self):
        grados = leer_angulo(self.var_angulo.get())
        if grados is None:
            messagebox.showwarning(self.t('angulo_t'), self.t('angulo_msg'))
            return
        self.girar_a(grados)

    def espejo_horizontal(self):
        if self.original is None:
            return
        self.mirall_h = not self.mirall_h
        self._aplicar_estado()

    def espejo_vertical(self):
        if self.original is None:
            return
        self.mirall_v = not self.mirall_v
        self._aplicar_estado()

    def restablecer(self):
        if self.original is None:
            return
        self.mirall_h = self.mirall_v = False
        self.rotacion = 0.0
        self._aplicar_estado()

    def orientacion_actual(self):
        return (self.mirall_h, self.mirall_v, self.rotacion)

    def describir_orientacion(self, orientacion):
        """Explica con palabras (en el idioma actual) cómo se ha modificado la foto."""
        mirall_h, mirall_v, rotacion = orientacion
        partes = []
        if mirall_h:
            partes.append(self.t('o_espejo_h'))
        if mirall_v:
            partes.append(self.t('o_espejo_v'))
        if rotacion:
            sentido = self.t('antihorario') if rotacion > 0 else self.t('horario')
            partes.append(self.t('o_giro', g=texto_angulo(rotacion), sentido=sentido))
        if not partes:
            return self.t('original')
        texto = ', '.join(partes)
        return texto[0].upper() + texto[1:]

    def descripcion_orientacion(self):
        return self.describir_orientacion(self.orientacion_actual())

    def _aplicar_transformacion(self):
        self.root.config(cursor='watch')
        self.root.update_idletasks()
        try:
            self.imagen = transformar_imagen(self.original, self.angulo, self.espejo)
        finally:
            self.root.config(cursor='')
        alto, ancho = self.imagen.shape
        self.centro_x = (ancho - 1) / 2.0
        self.centro_y = (alto - 1) / 2.0
        self._actualizar_orientacion()
        self._programar_dibujo()

    def _actualizar_orientacion(self):
        self.lbl_orientacion.config(text=self.descripcion_orientacion())

    # ------------------------------------------------------------------ vista
    def ajustar_vista(self):
        if self.imagen is None:
            return
        self.root.update_idletasks()
        ancho_v, alto_v = self._tamano_lienzo()
        if min(ancho_v, alto_v) < 50 and getattr(self, '_esperas_vista', 0) < 40:
            # la ventana aún no tiene su tamaño (se acaba de abrir): se encuadra en cuanto lo tenga
            self._esperas_vista = getattr(self, '_esperas_vista', 0) + 1
            self.root.after(50, self.ajustar_vista)
            return
        self._esperas_vista = 0
        alto, ancho = self.imagen.shape
        self.zoom = min(ancho_v / ancho, alto_v / alto) * 0.98
        self.centro_x = (ancho - 1) / 2.0
        self.centro_y = (alto - 1) / 2.0
        self._programar_dibujo()

    def _tamano_lienzo(self):
        return max(self.lienzo.winfo_width(), 10), max(self.lienzo.winfo_height(), 10)

    def lienzo_a_imagen(self, u, v):
        """Pasa de un punto de la ventana a coordenadas (x, y) de la imagen."""
        ancho_v, alto_v = self._tamano_lienzo()
        x = self.centro_x + (u + 0.5 - ancho_v / 2.0) / self.zoom
        y = self.centro_y - (v + 0.5 - alto_v / 2.0) / self.zoom
        return x, y

    def imagen_a_lienzo(self, x, y):
        ancho_v, alto_v = self._tamano_lienzo()
        u = (x - self.centro_x) * self.zoom + ancho_v / 2.0 - 0.5
        v = alto_v / 2.0 - (y - self.centro_y) * self.zoom - 0.5
        return u, v

    def _rueda(self, evento):
        self._zoom_en(evento.x, evento.y, 1.25 if evento.delta > 0 else 1 / 1.25)

    def _zoom_centro(self, factor):
        ancho_v, alto_v = self._tamano_lienzo()
        self._zoom_en(ancho_v / 2, alto_v / 2, factor)

    def _zoom_en(self, u, v, factor):
        if self.imagen is None:
            return
        x, y = self.lienzo_a_imagen(u, v)  # el punto bajo el ratón se queda quieto
        self.zoom = min(max(self.zoom * factor, 0.01), 64.0)
        ancho_v, alto_v = self._tamano_lienzo()
        self.centro_x = x - (u + 0.5 - ancho_v / 2.0) / self.zoom
        self.centro_y = y + (v + 0.5 - alto_v / 2.0) / self.zoom
        self._programar_dibujo()

    def _empezar_arrastre(self, evento):
        self._arrastre = (evento.x, evento.y, self.centro_x, self.centro_y)
        self.lienzo.config(cursor='fleur')

    def _arrastrar(self, evento):
        if self._arrastre is None or self.imagen is None:
            return
        u0, v0, cx0, cy0 = self._arrastre
        self.centro_x = cx0 - (evento.x - u0) / self.zoom
        self.centro_y = cy0 + (evento.y - v0) / self.zoom
        self._programar_dibujo()

    def _terminar_arrastre(self, _evento):
        self._arrastre = None
        self.lienzo.config(cursor='crosshair')

    def _mover_cursor(self, evento):
        if self.imagen is None:
            return
        x, y = self.lienzo_a_imagen(evento.x, evento.y)
        alto, ancho = self.imagen.shape
        ix, iy = int(math.floor(x + 0.5)), int(math.floor(y + 0.5))
        if 0 <= ix < ancho and 0 <= iy < alto:
            valor = self.imagen[iy, ix]
            brillo = f'{valor:.0f}' if np.isfinite(valor) else self.t('sin_dato')
            self.lbl_cursor.config(text=self.t('cursor', x=numero(x, 1), y=numero(y, 1), b=brillo))
        else:
            self.lbl_cursor.config(text=self.t('fuera'))

    def _ir_a_medicion(self, _evento):
        """Doble clic en la tabla: centra la vista en esa estrella (si se ve en esta foto)."""
        seleccion = self.tabla.selection()
        if not seleccion:
            return
        n = int(seleccion[0])
        for m in self.mediciones:
            if m['n'] == n:
                posicion = self.posicion_en_vista(m)
                if posicion:
                    self.centro_x, self.centro_y = posicion[0], posicion[1]
                    self.zoom = max(self.zoom, 4.0)
                    self._programar_dibujo()
                else:
                    self._mensaje(self.t('otra_foto'))
                return

    def _se_ve_ahora(self, m):
        return m['ruta'] == self.ruta and m['giro'] == self.angulo and m['espejo'] == self.espejo

    def posicion_en_vista(self, m):
        """Dónde se ve ahora la estrella medida `m`, aunque se midiera con la foto
        girada o en espejo de otra manera. Devuelve (x, y, misma_orientacion)
        o None si es de otra foto."""
        if m['ruta'] != self.ruta:
            return None
        if self._se_ve_ahora(m):
            return m['x'], m['y'], True
        foto = self.fotos[self.ruta]
        x0, y0 = destransformar_punto(m['x'], m['y'], foto['ancho'], foto['alto'], m['giro'], m['espejo'])
        x, y = transformar_punto(x0, y0, foto['ancho'], foto['alto'], self.angulo, self.espejo)
        return x, y, False

    # ------------------------------------------------------------------ dibujo
    def _programar_dibujo(self):
        if not self._dibujo_pendiente:
            self._dibujo_pendiente = True
            self.root.after_idle(self._dibujar)

    def _a_gris(self, valores):
        """Convierte brillos en niveles de gris 0-255 según el brillo elegido."""
        t = (valores - self.nivel_bajo) / (self.nivel_alto - self.nivel_bajo)
        sin_dato = ~np.isfinite(t)
        t = np.clip(np.nan_to_num(t, nan=0.0), 0.0, 1.0)
        if self.modo_brillo == 1:
            t = np.sqrt(t)
        elif self.modo_brillo == 2:
            t = np.log10(1.0 + 999.0 * t) / 3.0
        gris = (t * 255.0 + 0.5).astype(np.uint8)
        gris[sin_dato] = COLOR_FUERA
        return gris

    def _dibujar(self):
        self._dibujo_pendiente = False
        try:
            self.lienzo.delete('all')
        except tk.TclError:
            return  # el lienzo se está rehaciendo (cambio de idioma)
        ancho_v, alto_v = self._tamano_lienzo()
        if self.imagen is None:
            self._dibujar_bienvenida(ancho_v, alto_v)
            return
        self.escena = None  # con una foto abierta ya no hay cielo animado

        # Para cada píxel de la ventana, qué píxel de la imagen le toca.
        alto, ancho = self.imagen.shape
        if self.lente is not None and self.lente_theta > 0.5:
            gris = self._dibujo_con_lente(ancho_v, alto_v)
        else:
            xs = self.centro_x + (np.arange(ancho_v) + 0.5 - ancho_v / 2.0) / self.zoom
            ys = self.centro_y - (np.arange(alto_v) + 0.5 - alto_v / 2.0) / self.zoom
            ix = np.floor(xs + 0.5).astype(np.int64)
            iy = np.floor(ys + 0.5).astype(np.int64)
            valido_x = (ix >= 0) & (ix < ancho)
            valido_y = (iy >= 0) & (iy < alto)
            muestra = self.imagen[np.clip(iy, 0, alto - 1)[:, None], np.clip(ix, 0, ancho - 1)[None, :]]
            gris = self._a_gris(muestra)
            gris[~valido_y, :] = COLOR_FUERA
            gris[:, ~valido_x] = COLOR_FUERA

        # Imagen PGM (formato de grises que tkinter sabe mostrar)
        cabecera = f'P5 {ancho_v} {alto_v} 255\n'.encode('ascii')
        self.foto = tk.PhotoImage(data=cabecera + gris.tobytes(), format='PPM')
        self.lienzo.create_image(0, 0, anchor=tk.NW, image=self.foto)
        self.lbl_zoom.config(text=self.t('zoom_pct', p=f'{self.zoom * 100:.0f}'))

        if self.var_marcas.get():
            self._dibujar_marcas()

    def _dibujo_con_lente(self, ancho_v, alto_v):
        """La foto vista a través de un agujero negro en el centro de la ventana
        (lente gravitatoria). Solo cambia lo que se ve, no las medidas."""
        alto, ancho = self.imagen.shape
        u, v = np.meshgrid(np.arange(ancho_v) + 0.5, np.arange(alto_v) + 0.5)
        cx, cy, theta = ancho_v / 2.0, alto_v / 2.0, self.lente_theta
        dx, dy = u - cx, v - cy
        r = np.maximum(np.hypot(dx, dy), 1e-6)
        factor = (r - theta * theta / r) / r      # ecuación de la lente: β = θ - θE²/θ
        su, sv = cx + dx * factor, cy + dy * factor
        ix = np.floor(self.centro_x + (su - ancho_v / 2.0) / self.zoom + 0.5).astype(np.int64)
        iy = np.floor(self.centro_y - (sv - alto_v / 2.0) / self.zoom + 0.5).astype(np.int64)
        valido = (ix >= 0) & (ix < ancho) & (iy >= 0) & (iy < alto)
        gris = self._a_gris(self.imagen[np.clip(iy, 0, alto - 1), np.clip(ix, 0, ancho - 1)])
        gris[~valido] = COLOR_FUERA
        gris[r < theta * 0.62] = 0                                # sombra del agujero
        gris[np.abs(r - theta * 0.66) < 1.6] = 255                # anillo de fotones
        return gris

    def _dibujar_bienvenida(self, ancho_v, alto_v):
        cx, cy = ancho_v / 2, alto_v / 2
        # halo suave detrás del título (las estrellas pasan por encima)
        for i, color in enumerate(('#070b1a', '#0a1024', '#0d152e', '#101a38')):
            r = 300 - i * 55
            self.lienzo.create_oval(cx - r * 1.4, cy - 60 - r * 0.5, cx + r * 1.4, cy - 60 + r * 0.5,
                                    fill=color, outline='')
        azar = random.Random()
        galaxia = Galaxia(self.lienzo, ancho_v, alto_v, azar, velocidad=4)
        cielo = CampoEstelar(self.lienzo, ancho_v, alto_v, semilla=7, densidad=2200, velocidad=16,
                             destellos=True)
        logo = Logo(self.lienzo, cx, cy - 75, 44, centrado=True)
        ahora = time.perf_counter()
        self.escena = {'cielo': cielo, 'galaxia': galaxia, 'logo': logo, 'azar': azar,
                       'fugaces': [], 'ancho': ancho_v, 'alto': alto_v,
                       'proxima_fugaz': ahora + 1.2, 'proximo_satelite': ahora + azar.uniform(6, 12)}
        # tarjeta oscura detrás de las instrucciones, para que se lean bien
        self.lienzo.create_rectangle(cx - 270, cy - 36, cx + 270, cy + 150, fill='#0a0f22',
                                     outline=BORDE, tags=('frente',))
        self.lienzo.create_text(cx, cy - 18, fill=TEXTO_SUAVE, font=('Segoe UI', 12), tags=('frente',),
                                text=self.t('lema'))
        for i, texto in enumerate(self.t('pasos_bienvenida')):
            y = cy + 35 + i * 30
            self.lienzo.create_text(cx - 225, y, text=str(i + 1), fill=DORADO, anchor=tk.CENTER,
                                    font=('Segoe UI', 12, 'bold'))
            self.lienzo.create_oval(cx - 237, y - 12, cx - 213, y + 12, outline=DORADO)
            self.lienzo.create_text(cx - 200, y, text=texto, anchor=tk.W, fill='#c9d3f0',
                                    font=('Segoe UI', 12))

    def _dibujar_marcas(self):
        """Dibuja las estrellas medidas en esta foto, con su letra. Si se midieron
        con otra orientación, se recalcula dónde quedan ahora y se marcan con
        un círculo discontinuo (como referencia)."""
        self._dibujar_detecciones()   # estrellas detectadas y objetos de catálogo (si se han pedido)
        radio_apertura = leer_numero(self.var_apertura.get(), 5.0, 1.0, 50.0)
        for medida in self.mediciones:
            posicion = self.posicion_en_vista(medida)
            if posicion is None:
                continue
            x, y, misma_orientacion = posicion
            u, v = self.imagen_a_lienzo(x, y)
            r = max(radio_apertura * self.zoom, 6.0)
            if misma_orientacion:
                color = DORADO if medida is self.mediciones[-1] else COLOR_ACENTO
                if set(medida.get('calidad') or []) - {'BOR'}:
                    color = '#ffb36b'   # naranja: la medida tiene algún problema de calidad
                if self.var_aperturas.get() and medida.get('r_fot'):
                    for radio, trazo in ((medida['r_fot'], ()), (medida.get('r_int'), (3, 3)), (medida.get('r_ext'), (3, 3))):
                        if radio:
                            rr = radio * self.zoom
                            self.lienzo.create_oval(u - rr, v - rr, u + rr, v + rr, outline='#7be0a0', dash=trazo)
                self.lienzo.create_oval(u - r, v - r, u + r, v + r, outline=color, width=2)
                self.lienzo.create_line(u - r - 6, v, u - r * 0.4, v, fill=color)
                self.lienzo.create_line(u + r * 0.4, v, u + r + 6, v, fill=color)
                self.lienzo.create_line(u, v - r - 6, u, v - r * 0.4, fill=color)
                self.lienzo.create_line(u, v + r * 0.4, u, v + r + 6, fill=color)
                self.lienzo.create_text(u + r + 4, v - r - 4, text=medida['estrella'], fill=DORADO,
                                        anchor=tk.SW, font=('Segoe UI', 13, 'bold'))
            else:
                color = '#ff5cf4'  # magenta: medida con otra orientación
                r = r * 1.4 + 6    # algo más grande, para que no lo tape la estrella
                self.lienzo.create_oval(u - r, v - r, u + r, v + r, outline='#1a0018', width=5)  # contorno oscuro
                self.lienzo.create_oval(u - r, v - r, u + r, v + r, outline=color, width=3, dash=(7, 4))
                tx, ty = u + r * 0.75 + 2, v - r * 0.75 - 2
                self.lienzo.create_text(tx + 1, ty + 1, text=medida['estrella'], fill='#1a0018',
                                        anchor=tk.SW, font=('Segoe UI', 14, 'bold'))  # sombra
                self.lienzo.create_text(tx, ty, text=medida['estrella'], fill=color,
                                        anchor=tk.SW, font=('Segoe UI', 14, 'bold'))

    # ------------------------------------------------------------------ medir
    def _clic_izquierdo(self, evento):
        self.lienzo.focus_set()
        if self._parpadeo:
            self._parar_parpadeo()
            return
        if self.imagen is None:
            self.abrir_fits()
            return
        x, y = self.lienzo_a_imagen(evento.x, evento.y)
        self.medir_en(x, y)

    def _opciones_medida(self):
        return (int(leer_numero(self.var_busqueda.get(), 15, 3, 200)),
                leer_numero(self.var_apertura.get(), 5.0, 1.0, 50.0),
                leer_numero(self.var_umbral.get(), 3.0, 0.0, 50.0))

    def _cambio_correccion(self):
        """Al activar o desactivar la corrección de la contaminación se recalculan las medidas
        de la foto abierta (las de otras fotos, cuando se vuelvan a medir o recalcular)."""
        self.ajustes['corregir_contaminacion'] = self.var_corregir.get()
        if self.original is None:
            return
        for m in self.mediciones:
            if m['ruta'] == self.ruta and m.get('giro') is not None:
                self._completar_medida(m, self.original)
                if self.tabla.exists(str(m['n'])):
                    self.tabla.item(str(m['n']), values=self._valores_fila(m))
        self._actualizar_fotometria()

    def _radios_fotometria(self):
        """Radio de la apertura y del anillo del cielo (corregidos si no tienen sentido)."""
        r = leer_numero(self.var_r_fot.get(), 5.0, 1.0, 60.0)
        r_int = max(leer_numero(self.var_r_int.get(), 8.0, 2.0, 80.0), r + 1.0)
        r_ext = max(leer_numero(self.var_r_ext.get(), 14.0, 3.0, 120.0), r_int + 3.0)
        return r, r_int, r_ext

    def _original_xy(self, m):
        """Posición de una medida en la foto original (sin giro ni espejo)."""
        foto = self.fotos[m['ruta']]
        return destransformar_punto(m['x'], m['y'], foto['ancho'], foto['alto'], m['giro'], m['espejo'])

    def _fwhm_referencia(self, m):
        """Tamaño típico de las estrellas de la foto: mediana de las otras medidas (si hay 3 o más)."""
        otras = [o['fwhm'] for o in self.mediciones if o['ruta'] == m['ruta'] and o is not m and o.get('fwhm')]
        return float(np.median(otras)) if len(otras) >= 3 else None

    def _completar_medida(self, m, original=None):
        """Añade a una medida lo que depende de la foto: fotometría de apertura e
        indicadores de calidad (si se da la imagen original), magnitud instrumental y,
        si la foto tiene solución astrométrica, sus coordenadas celestes (AR, Dec)."""
        foto = self.fotos[m['ruta']]
        m['x0'], m['y0'] = self._original_xy(m)
        if original is not None:
            r, r_int, r_ext = self._radios_fotometria()
            fotometria = C.fotometria_apertura(original, m['x0'], m['y0'], r, r_int, r_ext,
                                               foto.get('ganancia'), foto.get('saturacion'))
            if fotometria is not None:
                m.update(flujo=fotometria['flujo'], sigma_flujo=fotometria['sigma_flujo'], snr=fotometria['snr'],
                         saturada=fotometria['saturada'], fondo_cielo=fotometria['fondo'], r_fot=r, r_int=r_int,
                         r_ext=r_ext, borde=bool(m.get('borde') or fotometria['borde']))
                m['contaminacion'], _vecinas = C.contaminacion(original, m['x0'], m['y0'], r, m.get('fwhm'),
                                                               fotometria['flujo'], fondo=fotometria['fondo'],
                                                               sigma=fotometria['sigma_fondo'])
                m['elongacion'] = C.elongacion_local(original, m['x0'], m['y0'], fotometria['fondo'],
                                                     radio=max(2.5, 1.5 * (m.get('fwhm') or 2.5)))
                m['corregida'] = False
                if self.var_corregir.get() and m['contaminacion'] > 0.01:
                    m['flujo'], m['sigma_flujo'] = C.flujo_corregido(m['flujo'], m['sigma_flujo'], m['contaminacion'])
                    m['corregida'] = True
            m['calidad'] = C.indicadores_calidad(m.get('saturada'), m.get('contaminacion'), m.get('borde'),
                                                 m.get('snr'), m.get('fwhm'), self._fwhm_referencia(m),
                                                 m.get('elongacion'))
        m['mag'], m['sigma_mag'] = magnitud_instrumental(m['flujo'], m['sigma_flujo'], foto.get('exposicion'))
        self._actualizar_cielo(m)

    def _actualizar_cielo(self, m):
        """AR y Dec de una medida con la solución astrométrica de su foto (si la tiene)."""
        foto = self.fotos[m['ruta']]
        m['ar'] = m['dec'] = None
        if foto.get('wcs'):
            x0, y0 = self._original_xy(m)
            m['ar'], m['dec'] = foto['wcs'].a_cielo(x0 + 1, y0 + 1)   # el FITS cuenta los píxeles desde 1

    def medir_en(self, x, y, letra=None, silencioso=False, radio_busqueda=None):
        busqueda, radio_apertura, umbral = self._opciones_medida()
        if radio_busqueda is not None:
            busqueda = radio_busqueda
        ganancia = self.fotos[self.ruta].get('ganancia')
        resultado = medir_estrella(self.imagen, x, y, busqueda, radio_apertura, umbral, ganancia)
        if resultado is None:
            if not silencioso:
                self._mensaje(self.t('sin_estrella'))
                self.root.bell()
            return None

        letra = letra or self.var_estrella.get() or 'A'
        # ¿Este centro ya es el de otra estrella medida (con otra letra) en esta foto? -> aviso y no se guarda
        for m in self.mediciones:
            if m['ruta'] != self.ruta or m['estrella'] == letra:
                continue
            ox, oy, _misma = self.posicion_en_vista(m)
            distancia = math.hypot(resultado['x'] - ox, resultado['y'] - oy)
            if distancia < 1.5:
                if not silencioso:
                    self.root.bell()
                    messagebox.showwarning(self.t('duplicada_t'),
                                           self.t('duplicada', otra=m['estrella'], letra=letra, d=numero(distancia, 2)))
                return None
        # Si esta estrella ya estaba medida en esta foto, la medida nueva sustituye a la anterior.
        anterior = [m for m in self.mediciones if m['ruta'] == self.ruta and m['estrella'] == letra]
        for m in anterior:
            self.mediciones.remove(m)
            self.tabla.delete(str(m['n']))

        self._contador += 1
        medida = dict(resultado)
        medida.update({
            'n': self._contador,
            'foto_n': self.fotos[self.ruta]['n'],
            'estrella': letra,
            'archivo': self.fotos[self.ruta]['archivo'],
            'ruta': self.ruta,
            'giro': self.angulo,
            'espejo': self.espejo,
            'orient': self.orientacion_actual(),
            'radio_apertura': radio_apertura,
            'fecha': datetime.now().isoformat(sep=' ', timespec='seconds'),
        })
        self._completar_medida(medida, self.original)
        self.mediciones.append(medida)
        self._registrar('hist_medida', letra=letra, n=medida['foto_n'], x=numero(medida['x'], 3), y=numero(medida['y'], 3))
        self._anadir_a_tabla(medida)
        self._actualizar_contador()
        if silencioso:
            return medida
        self._mostrar_resultado(medida)
        avisos = []
        if anterior:
            avisos.append(self.t('sustituida', letra=letra))
        if medida['borde']:
            avisos.append(self.t('borde'))
        for codigo in medida.get('calidad', []):
            if codigo != 'BOR':
                avisos.append('⚠ ' + self.t('calidad_' + codigo))
        self._mensaje('\n'.join(avisos))
        # pasar a la siguiente letra
        if self.var_avanzar.get() and letra != LETRAS_ESTRELLA[-1]:
            self.var_estrella.set(LETRAS_ESTRELLA[LETRAS_ESTRELLA.index(letra) + 1])
        self._actualizar_botones_serie()
        self._programar_dibujo()
        return medida

    def medir_como_anterior(self):
        """Mide en esta foto, de golpe, todas las estrellas que se midieron en la foto
        anterior, buscando cada una cerca de donde estaba allí (con la misma letra).
        Muy cómodo en series de fotos del mismo campo."""
        if self.imagen is None or not self.ruta_anterior:
            return
        de_antes = sorted((m for m in self.mediciones if m['ruta'] == self.ruta_anterior),
                          key=lambda m: LETRAS_ESTRELLA.index(m['estrella']))
        # Las estrellas que no se encontraron en la foto anterior (por ejemplo, una foto
        # movida) se buscan desde la última foto en la que sí se midieron.
        ya = {m['estrella'] for m in de_antes}
        for m in sorted(self.mediciones, key=lambda m: -m['foto_n']):
            if m['estrella'] not in ya and m['ruta'] != self.ruta:
                de_antes.append(m)
                ya.add(m['estrella'])
        foto_ahora = self.fotos[self.ruta]
        busqueda = self._opciones_medida()[0]
        # Primero se calcula cuánto se ha movido el campo entre las dos fotos (con la estrella
        # más brillante), y luego se busca cada estrella en su sitio corregido.
        # Los objetos que se mueven respecto a las estrellas (asteroides, cometas) se buscan
        # donde deberían estar ahora según su velocidad en las fotos anteriores.
        desplazamiento = (0.0, 0.0)
        referencia = None          # (letra, x, y en la foto original) de la estrella más brillante
        medidas_ok, fallidas, moviles = [], [], []
        for intento, m in enumerate(sorted(de_antes, key=lambda m: -m['flujo'])):
            prediccion = self._prediccion_movil(m['estrella'], referencia) if referencia else None
            if prediccion is not None:
                x, y = transformar_punto(prediccion[0], prediccion[1], foto_ahora['ancho'], foto_ahora['alto'],
                                         self.angulo, self.espejo)
                nueva = self.medir_en(x, y, letra=m['estrella'], silencioso=True, radio_busqueda=min(busqueda, 8))
                if nueva:
                    moviles.append(m['estrella'])
            else:
                x0, y0 = self._original_xy(m)
                x, y = transformar_punto(x0, y0, foto_ahora['ancho'], foto_ahora['alto'], self.angulo, self.espejo)
                nueva = self.medir_en(x + desplazamiento[0], y + desplazamiento[1], letra=m['estrella'],
                                      silencioso=True, radio_busqueda=busqueda if intento == 0 else min(busqueda, 8))
            if nueva:
                if intento == 0:
                    desplazamiento = (nueva['x'] - x, nueva['y'] - y)
                    referencia = (m['estrella'], *self._original_xy(nueva))
                medidas_ok.append(m['estrella'])
            else:
                fallidas.append(m['estrella'])
        texto = self.t('como_anterior_ok', n=len(medidas_ok), total=len(de_antes))
        if moviles:
            texto += '\n' + self.t('como_anterior_moviles', letras=', '.join(sorted(moviles)))
        if fallidas:
            texto += '\n' + self.t('como_anterior_faltan', letras=', '.join(sorted(fallidas)))
        if self.mediciones and medidas_ok:
            self._mostrar_resultado(self.mediciones[-1])
        self._mensaje(texto)
        self._actualizar_botones_serie()
        self._programar_dibujo()

    def _prediccion_movil(self, letra, referencia):
        """Si la estrella 'letra' se ha movido respecto a la estrella de referencia en las
        fotos anteriores (con hora), dónde debería estar ahora en la foto original.
        Devuelve (x, y) o None si no se mueve o no hay datos suficientes."""
        jd_ahora = self.fotos[self.ruta].get('jd')
        letra_ref, x_ref, y_ref = referencia
        if jd_ahora is None or letra == letra_ref:
            return None
        puntos = []
        for m in self.mediciones:
            if m['estrella'] != letra or m['ruta'] == self.ruta:
                continue
            ref = next((o for o in self.mediciones if o['ruta'] == m['ruta'] and o['estrella'] == letra_ref), None)
            jd = self.fotos[m['ruta']].get('jd')
            if ref is None or jd is None:
                continue
            (xm, ym), (xr, yr) = self._original_xy(m), self._original_xy(ref)
            puntos.append((jd, xm - xr, ym - yr, math.hypot(m['sigma_x'], m['sigma_y'])))
        if len(puntos) < 2:
            return None
        puntos.sort()
        tiempos, xs, ys, errores = zip(*puntos)
        if not C.parece_movil(tiempos, xs, ys, list(errores)):
            return None
        ajuste = C.ajuste_lineal_movimiento(tiempos, xs, ys)
        if ajuste is None:
            return None
        x0, y0, vx, vy = ajuste[:4]
        dt = jd_ahora - tiempos[0]
        return x_ref + x0 + vx * dt, y_ref + y0 + vy * dt

    def _mostrar_resultado(self, medida):
        self.lbl_titulo_medida.config(text=self.t('titulo_medida', letra=medida['estrella'], n=medida['foto_n']))
        self.lbl_x.config(text=f'X = {numero(medida["x"], 3)} ± {numero(medida["sigma_x"], 3)} px')
        self.lbl_y.config(text=f'Y = {numero(medida["y"], 3)} ± {numero(medida["sigma_y"], 3)} px')
        if medida.get('ar') is not None:
            self.lbl_cielo.config(text=f'AR  {ar_hms(medida["ar"])}\nDec {dec_dms(medida["dec"])}')
        else:
            self.lbl_cielo.config(text='')
        detalle = self.t('detalle', flujo=numero(medida['flujo'], 1),
                         fwhm=numero(medida['fwhm'], 2), ps=numero(medida['pico_sigma'], 0))
        if medida.get('snr') is not None and np.isfinite(medida['snr']):
            detalle += '\n' + self.t('detalle_snr', snr=numero(medida['snr'], 0), r=numero(medida.get('r_fot'), 1))
        if medida.get('mag') is not None:
            detalle += '\n' + self.t('detalle_mag', mag=numero(medida['mag'], 3), sm=numero(medida['sigma_mag'], 3))
        calibrada, error = self.magnitud_calibrada(medida)
        if calibrada is not None:
            clave = 'detalle_mag_rel' if self._fotometria_relativa() else 'detalle_mag_cal'
            detalle += '\n' + self.t(clave, mag=numero(calibrada, 3), sm=numero(error, 3),
                                     comp=' '.join(c[0] for c in self._comparaciones() or [] if c[0] != medida['estrella']))
        calidad = medida.get('calidad') or []
        detalle += '\n' + (self.t('calidad_bien') if not calidad else
                           self.t('calidad_titulo') + ' ' + ', '.join(self.t('calidad_corto_' + c) for c in calidad))
        self.lbl_detalle.config(text=detalle)

    def _texto_calidad(self, m):
        """Columna «Qualitat» de la tabla: ✓ o los códigos de los problemas."""
        if 'calidad' not in m:
            return ''
        return '✓' if not m['calidad'] else ' '.join(m['calidad'])

    def _valores_fila(self, m):
        error = max(m['sigma_x'], m['sigma_y'])
        calibrada, _err = self.magnitud_calibrada(m)
        if self._fotometria_activa():
            mag = numero(calibrada, 3) if calibrada is not None else ''
        else:
            mag = numero(m.get('mag'), 3)
        return (m['foto_n'], m['archivo'], m['estrella'], self.describir_orientacion(m['orient']),
                numero(m['x'], 3), numero(m['y'], 3), numero(error, 3), numero(m['flujo'], 1), mag,
                ar_hms(m['ar']) if m.get('ar') is not None else '',
                dec_dms(m['dec']) if m.get('dec') is not None else '', self._texto_calidad(m))

    # ------------------------------------------------------------------ fotometría diferencial
    def _panel_fotometria(self, panel):
        """4. Fotometría diferencial con una o varias estrellas de comparación.
        Con magnitudes de catálogo: magnitud real (calibrada). Sin ellas: magnitud
        relativa a la media de las comparaciones (sirve para curvas de luz)."""
        caja = ttk.LabelFrame(panel, text=self.t('paso4'), padding=8)
        caja.pack(fill=tk.X, pady=(0, 8))
        ttk.Label(caja, text=self.t('fotometria_texto'), wraplength=300, style='Gris.TLabel').pack(fill=tk.X)
        fila = ttk.Frame(caja)
        fila.pack(fill=tk.X, pady=(6, 0))
        ttk.Label(fila, text=self.t('comparaciones')).pack(side=tk.LEFT)
        entrada = tk.Entry(fila, textvariable=self.var_comps, justify=tk.LEFT,
                           font=('Consolas', 11, 'bold'), background=CAMPO, foreground=DORADO,
                           insertbackground=TEXTO, relief=tk.FLAT, highlightthickness=1,
                           highlightbackground=BORDE, highlightcolor=COLOR_ACENTO)
        entrada.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(6, 0), ipady=2)
        Ayuda(entrada, self.t('ayuda_comparaciones'))
        self.lbl_fotometria = ttk.Label(caja, text='', style='Gris.TLabel', wraplength=300, justify=tk.LEFT)
        self.lbl_fotometria.pack(fill=tk.X, pady=(4, 0))
        self.boton_fotometria = ttk.Button(caja, text=self.t('boton_fotometria'), command=self.ventana_fotometria)
        self.boton_fotometria.pack(fill=tk.X, pady=(6, 0))
        Ayuda(self.boton_fotometria, self.t('ayuda_boton_fotometria'))

    @staticmethod
    def _leer_lista_letras(texto, con_color=True):
        """«B=12,35 C=11,90/0,65 D» -> [(letra, valor o None, color o None)]; None si está mal escrito."""
        resultado = []
        for trozo in re.split(r'[\s;]+', texto.strip().upper()):
            if not trozo:
                continue
            m = re.match(r'^([A-Z])(?:=([-+]?\d+(?:[.,]\d*)?)(?:/([-+]?\d+(?:[.,]\d*)?))?)?$', trozo)
            if not m or (m.group(3) and not con_color):
                return None
            valor = float(m.group(2).replace(',', '.')) if m.group(2) else None
            color = float(m.group(3).replace(',', '.')) if m.group(3) else None
            if m.group(1) in [r[0] for r in resultado]:
                return None   # letra repetida
            resultado.append((m.group(1), valor, color))
        return resultado

    def _comparaciones(self):
        """Estrellas de comparación escritas: [(letra, magnitud de catálogo o None, índice de color o None)]."""
        return self._leer_lista_letras(self.var_comps.get())

    def _coef_color(self):
        valor = leer_angulo(self.var_color_t.get())   # (acepta coma o punto decimal)
        return valor or 0.0

    def _indices_color(self):
        lista = self._leer_lista_letras(self.var_colores.get(), con_color=False) or []
        return {letra: valor for letra, valor, _c in lista if valor is not None}

    def _fotometria_activa(self):
        return bool(self._comparaciones())

    def _fotometria_relativa(self):
        """Sin la magnitud de catálogo de alguna comparación: magnitudes relativas."""
        comps = self._comparaciones()
        return bool(comps) and any(c[1] is None for c in comps)

    def magnitud_calibrada(self, m):
        """Magnitud de la medida m por fotometría diferencial con las estrellas de
        comparación medidas en la misma foto (ver ciencia.punto_cero). Si m es una de las
        comparaciones, se calcula con las demás (así sirve de comprobación)."""
        comps = self._comparaciones()
        if not comps or m.get('mag') is None:
            return None, None
        relativa = any(c[1] is None for c in comps)
        lista = []
        for letra, catalogo, color in comps:
            if letra == m['estrella']:
                continue
            c = next((c for c in self.mediciones if c['ruta'] == m['ruta'] and c['estrella'] == letra), None)
            if c is None or c.get('mag') is None:
                continue
            lista.append((c['mag'], c['sigma_mag'], None if relativa else catalogo, color))
        if not lista:
            return None, None
        t = self._coef_color()
        return C.magnitud_diferencial(m['mag'], m['sigma_mag'], C.punto_cero(lista, t), t,
                                      self._indices_color().get(m['estrella']))

    def _actualizar_fotometria(self):
        """Al cambiar las comparaciones: rehace la columna de magnitud y los avisos."""
        if not hasattr(self, 'lbl_fotometria') or not hasattr(self, 'tabla'):
            return
        try:
            comps = self._comparaciones()
            texto = ''
            if comps is None:
                texto = self.t('fot_formato')
            elif comps:
                letras = ' '.join(c[0] for c in comps)
                if any(c[1] is None for c in comps):
                    texto = self.t('fot_relativa', letras=letras)
                else:
                    texto = self.t('fot_calibrada', letras=letras, n=len(comps))
                fotos = {m['ruta'] for m in self.mediciones}
                for letra, _cat, _col in comps:
                    sin = [self.fotos[r]['n'] for r in fotos
                           if not any(c['ruta'] == r and c['estrella'] == letra for c in self.mediciones)]
                    if sin:
                        texto += '\n' + self.t('fotometria_sin_comp', letra=letra,
                                               fotos=', '.join(str(n) for n in sorted(sin)))
            self.lbl_fotometria.config(text=texto, foreground=ROJO_SUAVE if comps is None else TEXTO_SUAVE)
            if self._fotometria_relativa():
                cabecera = self.t('col_mag_rel')
            else:
                cabecera = self.t('col_mag_cal') if self._fotometria_activa() else self.t('col_mag_inst')
            self.tabla.heading('mag', text=cabecera)
            for m in self.mediciones:
                if self.tabla.exists(str(m['n'])):
                    self.tabla.item(str(m['n']), values=self._valores_fila(m))
            if self.mediciones:
                self._mostrar_resultado(self.mediciones[-1])
        except tk.TclError:
            pass

    def _anadir_a_tabla(self, m):
        etiquetas = ('par',) if m['foto_n'] % 2 == 0 else ()  # color alterno por foto
        self.tabla.insert('', tk.END, iid=str(m['n']), values=self._valores_fila(m), tags=etiquetas)
        self.tabla.see(str(m['n']))

    def _programar_recalculo(self):
        """Espera un momento (por si se siguen tocando las flechas) y recalcula."""
        if self._recalculo_pendiente:
            self.root.after_cancel(self._recalculo_pendiente)
        self._recalculo_pendiente = self.root.after(350, self.recalcular_medida)

    def recalcular_medida(self):
        """Vuelve a calcular con las opciones avanzadas actuales la medida seleccionada en
        la tabla (o, si no hay ninguna, la última), si es de la foto que se está viendo."""
        self._recalculo_pendiente = None
        if self.imagen is None or not self.mediciones:
            return
        seleccion = self.tabla.selection()
        m = next((m for m in self.mediciones if seleccion and str(m['n']) == seleccion[0]), self.mediciones[-1])
        posicion = self.posicion_en_vista(m)
        if posicion is None or not posicion[2]:
            return  # es de otra foto u orientación: no se puede recalcular aquí
        radio_busqueda, radio_apertura, umbral = self._opciones_medida()
        # se busca muy cerca de donde estaba, para no saltar a otra estrella
        resultado = medir_estrella(self.imagen, m['x'], m['y'], min(radio_busqueda, 4), radio_apertura, umbral,
                                   self.fotos[m['ruta']].get('ganancia'))
        if resultado is None:
            self._mensaje(self.t('sin_estrella'))
            return
        m.update(resultado)
        m['radio_apertura'] = radio_apertura
        self._completar_medida(m, self.original)
        self.tabla.item(str(m['n']), values=self._valores_fila(m))
        self._mostrar_resultado(m)
        self._mensaje(self.t('recalculada', letra=m['estrella']))
        self._programar_dibujo()

    def _actualizar_contador(self):
        n = len(self.mediciones)
        if n == 0:
            texto = self.t('ninguna_todavia')
        elif n == 1:
            texto = self.t('contador_1')
        else:
            texto = self.t('contador_n', n=n)
        self.lbl_num_medidas.config(text=texto)
        if self._fotometria_activa():  # una medida nueva (o borrada) puede cambiar las magnitudes calibradas
            self._actualizar_fotometria()

    def borrar_ultima(self):
        if not self.mediciones:
            return
        medida = self.mediciones.pop()
        self.tabla.delete(str(medida['n']))
        if medida['ruta'] == self.ruta:
            self.var_estrella.set(medida['estrella'])  # así el próximo clic vuelve a medir esa estrella
        self._actualizar_contador()
        self._programar_dibujo()

    def _actualizar_boton_borrar(self):
        """Muestra el botón «Borrar la medida» solo si hay filas seleccionadas."""
        seleccion = self.tabla.selection()
        if seleccion:
            texto = self.t('borrar_medida') if len(seleccion) == 1 else self.t('borrar_medidas', n=len(seleccion))
            self.boton_borrar_sel.config(text=texto)
            if not self.boton_borrar_sel.winfo_ismapped():
                self.boton_borrar_sel.pack(side=tk.RIGHT, padx=(0, 6))
        elif self.boton_borrar_sel.winfo_ismapped():
            self.boton_borrar_sel.pack_forget()

    def borrar_seleccionadas(self):
        """Borra las medidas seleccionadas en la tabla."""
        seleccion = set(self.tabla.selection())
        if not seleccion:
            return
        borradas = [m for m in self.mediciones if str(m['n']) in seleccion]
        for m in borradas:
            self.mediciones.remove(m)
            self.tabla.delete(str(m['n']))
        if len(borradas) == 1 and borradas[0]['ruta'] == self.ruta:
            self.var_estrella.set(borradas[0]['estrella'])  # así el próximo clic vuelve a medir esa estrella
        self._actualizar_boton_borrar()
        self._actualizar_contador()
        self._programar_dibujo()

    def borrar_todas(self):
        """Botón rojo «Limpiar medidas»: borra todo y empieza de cero (pide confirmación)."""
        if not self.mediciones:
            return
        if not messagebox.askyesno(self.t('borrar_t'), self.t('borrar_msg_n', n=len(self.mediciones)),
                                   icon=messagebox.WARNING, default=messagebox.NO):
            return
        self.mediciones.clear()
        self.tabla.delete(*self.tabla.get_children())
        self.var_estrella.set('A')
        # las fotos se vuelven a numerar desde el 1 (la abierta ahora pasa a ser la n.º 1)
        if self.ruta is not None:
            foto = dict(self.fotos[self.ruta], n=1)
            self.fotos = {self.ruta: foto}
            self._mostrar_info_foto()
        else:
            self.fotos = {}
        self.lbl_titulo_medida.config(text=self.t('ultima'))
        self.lbl_x.config(text='X = -')
        self.lbl_y.config(text='Y = -')
        self.lbl_detalle.config(text='')
        self._mensaje('')
        self._actualizar_contador()
        self._programar_dibujo()

    def cerrar_foto(self):
        """Quita la foto y vuelve a la pantalla de inicio (las medidas se conservan)."""
        self._parar_parpadeo()
        if self.ruta is None:
            return
        self.ruta = self.original = self.imagen = None
        self.nota = ''
        self.mirall_h = self.mirall_v = False
        self.rotacion = self.angulo = 0.0
        self.espejo = False
        self.var_angulo.set('0')
        self.lente, self.lente_theta = None, 0.0
        self.root.title(TITULO_VENTANA)
        self.lbl_archivo.config(text=self.t('ninguna_foto'), style='Gris.TLabel')
        self.boton_cerrar_foto.pack_forget()
        self.lbl_cielo.config(text='')
        self._actualizar_botones_serie()
        self._actualizar_orientacion()
        self.lbl_zoom.config(text='')
        self.lbl_cursor.config(text='')
        self._mensaje('')
        self._programar_dibujo()  # vuelve a dibujar el cielo animado de bienvenida

    # ------------------------------------------------------------------ exportar
    def guardar_csv(self):
        """Botón 'Desar per a Excel…': guarda en .xlsx (recomendado) o en .csv."""
        if not self.mediciones:
            messagebox.showinfo(self.t('guardar_t'), self.t('sin_medidas'))
            return
        ruta = filedialog.asksaveasfilename(
            title=self.t('guardar_titulo'), defaultextension='.xlsx',
            initialfile=self.t('archivo_defecto') + '.xlsx',
            filetypes=[(self.t('tipo_xlsx'), '*.xlsx'), (self.t('tipo_csv'), '*.csv')])
        if not ruta:
            return
        if not ruta.lower().endswith(('.xlsx', '.csv')):
            ruta += '.xlsx'
        try:
            if ruta.lower().endswith('.csv'):
                self.escribir_csv(ruta)
            else:
                self.escribir_xlsx(ruta)
        except PermissionError:
            messagebox.showerror(self.t('error_guardar'), self.t('abierto_otro'))
            return
        except Exception as error:
            messagebox.showerror(self.t('error_guardar'), str(error))
            return
        fotos = len({m['ruta'] for m in self.mediciones})
        self._guardado = True
        self._registrar('hist_excel', ruta=ruta, n=len(self.mediciones))
        messagebox.showinfo(self.t('guardado_t'), self.t('guardado', n=len(self.mediciones), f=fotos, ruta=ruta))

    # ------------------------------------------------------------------ informe MPC (80 columnas y ADES)
    def _letras_para_mpc(self):
        """Estrellas que tienen al menos una medida con AR/Dec y hora (se pueden enviar)."""
        letras = {m['estrella'] for m in self.mediciones
                  if m.get('ar') is not None and self.fotos[m['ruta']].get('jd')}
        return sorted(letras, key=LETRAS_ESTRELLA.index)

    def _observaciones_mpc(self, letra, con_magnitud):
        """Observaciones del objeto 'letra' ordenadas por tiempo, con sus errores en ″.

        El error de posición junta el estadístico del centroide (pasado a ″ con la escala
        de la foto) y, si se ha comprobado la astrometría con Gaia, el sistemático de la
        solución astrométrica de esa foto."""
        usar_mag = con_magnitud and self._fotometria_activa() and not self._fotometria_relativa()
        observaciones = []
        for m in sorted((m for m in self.mediciones if m['estrella'] == letra),
                        key=lambda m: self.fotos[m['ruta']].get('jd') or 0):
            foto = self.fotos[m['ruta']]
            if m.get('ar') is None or not foto.get('jd') or not foto.get('wcs'):
                continue
            estadistico = math.hypot(m['sigma_x'], m['sigma_y']) / math.sqrt(2) * foto['wcs'].escala
            rms = max(math.hypot(estadistico, foto.get('rms_astrometrico') or 0.0), 0.01)
            magnitud, error = self.magnitud_calibrada(m) if usar_mag else (None, None)
            observaciones.append({'jd': foto['jd'], 'ar': m['ar'], 'dec': m['dec'], 'rms_ar': rms, 'rms_dec': rms,
                                  'mag': magnitud, 'rms_mag': error})
        return observaciones

    @staticmethod
    def _contexto_ades(datos):
        return {'designacion': datos['designacion'], 'codigo': datos['codigo'] or 'XXX',
                'modo': 'CMO' if datos['tipo'] == 'B' else 'CCD', 'remitente': datos['remitente'],
                'observadores': datos['observadores'], 'medidores': datos['medidores'],
                'apertura': datos['apertura'].replace(',', '.'), 'diseno': datos['diseno'],
                'catalogo_astrometrico': datos['cat_ast'], 'catalogo_fotometrico': datos['cat_fot'],
                'banda': datos['banda'], 'descubrimiento': datos['descubrimiento'],
                'comentario': f'Measured with Emillium Expanded {VERSION}'}

    def texto_informe_mpc(self, datos):
        """Informe completo en el formato elegido. Devuelve (texto, número de observaciones)."""
        observaciones = self._observaciones_mpc(datos['letra'], datos['con_magnitud'])
        formato = datos.get('formato', 'obs80')
        if formato == 'psv':
            return ades.crear_psv(self._contexto_ades(datos), observaciones), len(observaciones)
        if formato == 'xml':
            return ades.crear_xml(self._contexto_ades(datos), observaciones), len(observaciones)
        cabecera = [f"COD {datos['codigo'].strip()[:3] or 'XXX'}"]
        for clave, campo in (('CON', 'contacto'), ('OBS', 'observadores'), ('MEA', 'medidores'),
                             ('TEL', 'telescopio'), ('NET', 'catalogo')):
            if datos[campo].strip():
                cabecera.append(f'{clave} {datos[campo].strip()}'[:80])
        cabecera.append(f'COM {self.t("mpc_comentario", v=VERSION)}'[:80])
        lineas = [linea_mpc(datos['designacion'], datos['descubrimiento'] and i == 0,   # el * solo en una línea
                            datos['tipo'], o['jd'], o['ar'], o['dec'], o['mag'], datos['banda'], datos['codigo'] or 'XXX')
                  for i, o in enumerate(observaciones)]
        return '\n'.join(cabecera + lineas) + '\n', len(observaciones)

    def validar_informe_mpc(self, datos):
        """Comprobaciones antes de exportar (ver ades.validar). En el formato de 80
        columnas no cuentan las que solo pide ADES."""
        problemas = ades.validar(self._contexto_ades(datos), self._observaciones_mpc(datos['letra'], datos['con_magnitud']))
        if datos.get('formato', 'obs80') == 'obs80':
            solo_ades = {'falta_remitente', 'falta_observadores', 'falta_medidores', 'apertura_rango', 'apertura_metros',
                         'apertura_numero', 'diseno_raro', 'catalogo_raro', 'sin_photcat', 'trksub_formato', 'sin_rms',
                         'rms_grande', 'rms_cero'}
            problemas = [p for p in problemas if p[1] not in solo_ades]
            if len(datos['designacion'].strip()) > 7:
                problemas.append(('error', 'designacion_80', {}))
        return problemas

    def texto_problemas(self, problemas):
        lineas = []
        for nivel, codigo, datos in problemas:
            datos = {k: (numero(v, 2) if isinstance(v, float) else v) for k, v in datos.items()}
            lineas.append(('✖ ' if nivel == 'error' else '⚠ ') + self.t('ades_' + codigo, **datos))
        return '\n'.join(lineas)

    def abrir_informe_mpc(self):
        """Asistente: genera un informe de astrometría para el Minor Planet Center."""
        letras = self._letras_para_mpc()
        if not letras:
            messagebox.showinfo(self.t('mpc_titulo'), self.t('mpc_sin_datos'))
            return
        v, marco = self._ventana(self.t('mpc_titulo'), 'ventana_mpc')
        ttk.Label(marco, text=self.t('mpc_explicacion'), wraplength=900, style='Gris.TLabel',
                  justify=tk.LEFT).grid(row=0, column=0, columnspan=4, sticky=tk.W, pady=(0, 8))
        guardados = self.ajustes.get('mpc', {})
        formatos = [('obs80', self.t('mpc_formato_80')), ('psv', 'ADES PSV'), ('xml', 'ADES XML')]
        var_formato = tk.StringVar(value=dict(formatos).get(guardados.get('formato', 'psv'), 'ADES PSV'))
        variables = {
            'letra': tk.StringVar(value=self.var_estrella.get() if self.var_estrella.get() in letras else letras[0]),
            'designacion': tk.StringVar(value=guardados.get('designacion', 'EML0001')),
            'codigo': tk.StringVar(value=guardados.get('codigo', 'XXX')),
            'tipo': tk.StringVar(value=guardados.get('tipo', 'C')),
            'observadores': tk.StringVar(value=guardados.get('observadores', '')),
            'medidores': tk.StringVar(value=guardados.get('medidores', '')),
            'remitente': tk.StringVar(value=guardados.get('remitente', '')),
            'contacto': tk.StringVar(value=guardados.get('contacto', '')),
            'telescopio': tk.StringVar(value=guardados.get('telescopio', '')),
            'apertura': tk.StringVar(value=guardados.get('apertura', '0.20')),
            'diseno': tk.StringVar(value=guardados.get('diseno', 'Reflector')),
            'catalogo': tk.StringVar(value=guardados.get('catalogo', 'Gaia DR3')),
            'cat_ast': tk.StringVar(value=guardados.get('cat_ast', 'Gaia3')),
            'cat_fot': tk.StringVar(value=guardados.get('cat_fot', 'Gaia3')),
            'banda': tk.StringVar(value=guardados.get('banda', 'G')),
        }
        descubrimiento = tk.BooleanVar(value=False)
        con_magnitud = tk.BooleanVar(value=self._fotometria_activa() and not self._fotometria_relativa())
        campos = (('mpc_formato', None, 'formato', [n for _c, n in formatos]),
                  ('mpc_objeto', 'letra', 'combo', letras), ('mpc_designacion', 'designacion', 'entrada', None),
                  ('mpc_codigo', 'codigo', 'entrada', None), ('mpc_tipo', 'tipo', 'combo', ['C', 'B']),
                  ('mpc_observadores', 'observadores', 'entrada', None), ('mpc_medidores', 'medidores', 'entrada', None),
                  ('mpc_remitente', 'remitente', 'entrada', None), ('mpc_contacto', 'contacto', 'entrada', None),
                  ('mpc_telescopio', 'telescopio', 'entrada', None), ('mpc_apertura', 'apertura', 'entrada', None),
                  ('mpc_diseno', 'diseno', 'combo', ades.DISENOS), ('mpc_catalogo', 'catalogo', 'entrada', None),
                  ('mpc_cat_ast', 'cat_ast', 'combo', ades.CATALOGOS_ASTROMETRICOS),
                  ('mpc_cat_fot', 'cat_fot', 'combo', ades.CATALOGOS_FOTOMETRICOS),
                  ('mpc_banda', 'banda', 'combo', ades.BANDAS))
        for i, (clave_texto, clave, clase, valores) in enumerate(campos):
            fila, columna = 1 + i // 2, (i % 2) * 2
            ttk.Label(marco, text=self.t(clave_texto)).grid(row=fila, column=columna, sticky=tk.W, padx=(0, 6), pady=3)
            if clase == 'formato':
                widget = ttk.Combobox(marco, textvariable=var_formato, values=valores, state='readonly', width=22)
            elif clase == 'combo':
                widget = ttk.Combobox(marco, textvariable=variables[clave], values=valores, state='readonly', width=12)
            else:
                widget = ttk.Entry(marco, textvariable=variables[clave], width=34)
            widget.grid(row=fila, column=columna + 1, sticky=tk.W, padx=(0, 16), pady=3)
        fila = 2 + len(campos) // 2
        ttk.Checkbutton(marco, text=self.t('mpc_descubrimiento'), variable=descubrimiento).grid(
            row=fila, column=0, columnspan=2, sticky=tk.W, pady=(4, 0))
        ttk.Checkbutton(marco, text=self.t('mpc_con_magnitud'), variable=con_magnitud).grid(
            row=fila, column=2, columnspan=2, sticky=tk.W, pady=(4, 0))
        vista = tk.Text(marco, width=118, height=12, font=('Consolas', 9), background=CAMPO, foreground=TEXTO,
                        insertbackground=TEXTO, relief=tk.FLAT, highlightthickness=1, highlightbackground=BORDE, wrap=tk.NONE)
        vista.grid(row=fila + 1, column=0, columnspan=4, sticky=tk.EW, pady=(10, 4))
        aviso = ttk.Label(marco, text='', foreground=TEXTO_SUAVE, wraplength=900, justify=tk.LEFT)
        aviso.grid(row=fila + 2, column=0, columnspan=4, sticky=tk.W)
        problemas_lbl = ttk.Label(marco, text='', foreground=ROJO_SUAVE, wraplength=900, justify=tk.LEFT)
        problemas_lbl.grid(row=fila + 3, column=0, columnspan=4, sticky=tk.W)

        def datos_actuales():
            datos = {k: var.get() for k, var in variables.items()}
            datos.update(descubrimiento=descubrimiento.get(), con_magnitud=con_magnitud.get(),
                         formato=next(c for c, n in formatos if n == var_formato.get()))
            return datos

        def actualizar(*_):
            datos = datos_actuales()
            texto, n = self.texto_informe_mpc(datos)
            vista.delete('1.0', tk.END)
            vista.insert('1.0', texto)
            aviso.config(text=self.t('mpc_n_obs', n=n) + '  ' + self.t('mpc_no_garantia'))
            problemas = self.validar_informe_mpc(datos)
            problemas_lbl.config(text=self.texto_problemas(problemas) or self.t('mpc_valido'),
                                 foreground=ROJO_SUAVE if ades.hay_errores(problemas) else
                                 ('#ffb36b' if problemas else '#7be0a0'))
        for var in list(variables.values()) + [descubrimiento, con_magnitud, var_formato]:
            var.trace_add('write', actualizar)

        def recordar():
            guardado = {k: var.get() for k, var in variables.items() if k != 'letra'}
            guardado['formato'] = datos_actuales()['formato']
            self.ajustes['mpc'] = guardado
            guardar_ajustes(self.ajustes)

        def guardar():
            datos = datos_actuales()
            problemas = self.validar_informe_mpc(datos)
            if ades.hay_errores(problemas):
                messagebox.showerror(self.t('mpc_titulo'), self.t('mpc_con_errores') + '\n\n' +
                                     self.texto_problemas([p for p in problemas if p[0] == 'error']), parent=v)
                return
            extension = {'obs80': '.txt', 'psv': '.psv', 'xml': '.xml'}[datos['formato']]
            ruta = filedialog.asksaveasfilename(parent=v, title=self.t('mpc_titulo'), defaultextension=extension,
                                                initialfile=f"MPC_{datos['designacion'].strip() or 'obs'}{extension}",
                                                filetypes=[(extension[1:].upper(), '*' + extension)])
            if not ruta:
                return
            try:
                with open(ruta, 'w', encoding='utf-8' if datos['formato'] == 'xml' else 'ascii', errors='replace',
                          newline='\n') as archivo:
                    archivo.write(self.texto_informe_mpc(datos)[0])
            except OSError as error:
                messagebox.showerror(self.t('error_guardar'), str(error), parent=v)
                return
            recordar()
            self._registrar('hist_mpc', ruta=ruta, formato=datos['formato'])
            messagebox.showinfo(self.t('guardado_t'), self.t('mpc_guardado', ruta=ruta), parent=v)

        def copiar():
            v.clipboard_clear()
            v.clipboard_append(self.texto_informe_mpc(datos_actuales())[0])
            recordar()
            aviso.config(text=self.t('mpc_copiado'))

        botones = ttk.Frame(marco)
        botones.grid(row=fila + 4, column=0, columnspan=4, sticky=tk.E, pady=(8, 0))
        ttk.Button(botones, text=self.t('mpc_desar'), style='Grande.TButton', command=guardar).pack(side=tk.RIGHT)
        ttk.Button(botones, text=self.t('mpc_copiar'), command=copiar).pack(side=tk.RIGHT, padx=6)
        ttk.Button(botones, text=self.t('mpc_cerrar'), command=lambda: (recordar(), v.destroy())).pack(side=tk.RIGHT)
        actualizar()
        self._centrar(v)

    # ------------------------------------------------------------------ tablas para Excel
    def tabla_por_estrellas(self):
        """Tabla "estrella por estrella": una fila por foto y, para cada estrella
        (A, B, C...), sus columnas (X, Y, errores, flujo, magnitud y, si la foto tiene
        solución astrométrica, AR y Dec). Devuelve (columnas, filas)."""
        letras = sorted({m['estrella'] for m in self.mediciones}, key=LETRAS_ESTRELLA.index)
        rutas = sorted({m['ruta'] for m in self.mediciones}, key=lambda r: self.fotos[r]['n'])
        con_hora = any(self.fotos[r].get('hora') for r in rutas)
        con_jd = any(self.fotos[r].get('jd') for r in rutas)
        con_expo = any(self.fotos[r].get('exposicion') for r in rutas)
        con_cielo = any(m.get('ar') is not None for m in self.mediciones)
        con_cal = self._fotometria_activa()
        nombre_cal = 'mag_rel_col' if self._fotometria_relativa() else 'mag_cal_col'
        columnas = [(self.t('col_foto_n_largo'), 0, 11), (self.t('col_foto'), None, 28)]
        if con_hora:
            columnas.append((self.t('col_hora'), None, 24))
        if con_jd:
            columnas.append((self.t('col_jd'), 6, 17))
        con_hjd = con_jd and any(self._hjd_foto(r) for r in rutas)
        if con_hjd:
            columnas.append((self.t('col_hjd'), 6, 17))
            columnas.append((self.t('col_bjd'), 6, 17))
        if con_expo:
            columnas.append((self.t('col_exposicion'), 1, 10))
        columnas.append((self.t('col_orientacion'), None, 30))
        for letra in letras:
            columnas += [(self.t('x_col', l=letra), 3, 12), (self.t('y_col', l=letra), 3, 12),
                         (f'{letra}_σX (px)', 3, 9), (f'{letra}_σY (px)', 3, 9),
                         (self.t('flujo_col', l=letra), 1, 12), (f'{letra}_Mag inst.', 3, 9), (f'{letra}_σMag', 3, 9)]
            if con_cal:
                columnas += [(self.t(nombre_cal, l=letra), 3, 12), (self.t('sigma_cal_col', l=letra), 3, 10)]
            if con_cielo:
                columnas += [(self.t('ar_col', l=letra), 6, 13), (self.t('dec_col', l=letra), 6, 13)]
            columnas.append((self.t('calidad_col', l=letra), None, 12))
        filas = []
        for ruta in rutas:
            foto = self.fotos[ruta]
            de_esta_foto = {m['estrella']: m for m in self.mediciones if m['ruta'] == ruta}
            orientaciones = []
            for m in sorted(de_esta_foto.values(), key=lambda m: LETRAS_ESTRELLA.index(m['estrella'])):
                texto = self.describir_orientacion(m['orient'])
                if texto not in orientaciones:
                    orientaciones.append(texto)
            fila = [foto['n'], foto['archivo']]
            if con_hora:
                fila.append(foto.get('hora'))
            if con_jd:
                fila.append(foto.get('jd'))
            if con_hjd:
                fila += [self._hjd_foto(ruta), self._bjd_foto(ruta)]
            if con_expo:
                fila.append(foto.get('exposicion'))
            fila.append(' / '.join(orientaciones))
            for letra in letras:
                m = de_esta_foto.get(letra)
                if m:
                    fila += [m['x'], m['y'], m['sigma_x'], m['sigma_y'], m['flujo'], m.get('mag'), m.get('sigma_mag')]
                    if con_cal:
                        fila += list(self.magnitud_calibrada(m))
                    if con_cielo:
                        fila += [m.get('ar'), m.get('dec')]
                    fila.append(self._texto_calidad(m) or None)
                else:
                    fila += [None] * (8 + (2 if con_cal else 0) + (2 if con_cielo else 0))
            filas.append(fila)
        return columnas, filas

    def tabla_todas(self):
        """Lista completa: una fila por medida, ordenada por foto y estrella."""
        columnas = [(self.t('col_foto_n_largo'), 0, 11), (self.t('col_foto'), None, 28),
                    (self.t('col_estrella'), None, 9), (self.t('col_orientacion'), None, 30),
                    ('X (px)', 3, 12), ('Y (px)', 3, 12), ('σX (px)', 3, 9), ('σY (px)', 3, 9),
                    (self.t('col_flujo'), 1, 12), (self.t('col_sigma_flujo'), 1, 10), (self.t('col_snr'), 1, 8),
                    (self.t('col_r_fot'), 1, 9), (self.t('col_fondo'), 2, 11),
                    ('Mag inst.', 3, 9), ('σMag', 3, 9), (self.t('col_mag_cal'), 3, 12), ('σMag cal.', 3, 10),
                    ('FWHM (px)', 2, 10), (self.t('col_calidad'), None, 14),
                    (self.t('col_ar_grados'), 6, 13), (self.t('col_dec_grados'), 6, 13),
                    (self.t('col_ar_hms'), None, 17), (self.t('col_dec_dms'), None, 17),
                    (self.t('col_hora'), None, 24), (self.t('col_jd'), 6, 17)]
        orden = sorted(self.mediciones, key=lambda m: (m['foto_n'], LETRAS_ESTRELLA.index(m['estrella'])))
        filas = []
        for m in orden:
            foto = self.fotos[m['ruta']]
            hay_cielo = m.get('ar') is not None
            calibrada, error_cal = self.magnitud_calibrada(m)
            filas.append([m['foto_n'], m['archivo'], m['estrella'], self.describir_orientacion(m['orient']),
                          m['x'], m['y'], m['sigma_x'], m['sigma_y'], m['flujo'], m['sigma_flujo'],
                          m.get('snr') if m.get('snr') is not None and np.isfinite(m['snr']) else None,
                          m.get('r_fot'), m.get('fondo_cielo'),
                          m.get('mag'), m.get('sigma_mag'), calibrada, error_cal, m['fwhm'], self._texto_calidad(m) or None,
                          m.get('ar'), m.get('dec'),
                          ar_hms(m['ar']) if hay_cielo else None, dec_dms(m['dec']) if hay_cielo else None,
                          foto.get('hora') or None, foto.get('jd')])
        return columnas, filas

    def tabla_configuracion(self):
        """Hoja con cómo se ha medido todo (para poder repetirlo): versión, opciones,
        comparaciones, soluciones astrométricas comprobadas..."""
        r, r_int, r_ext = self._radios_fotometria()
        busqueda, apertura, umbral = self._opciones_medida()
        filas = [[self.t('conf_version'), f'Emillium Expanded {VERSION}'],
                 [self.t('conf_fecha'), datetime.now().isoformat(sep=' ', timespec='seconds')],
                 [self.t('adv_busqueda'), busqueda], [self.t('adv_apertura'), apertura], [self.t('adv_umbral'), umbral],
                 [self.t('adv_r_fot'), r], [self.t('adv_r_int'), r_int], [self.t('adv_r_ext'), r_ext],
                 [self.t('comparaciones'), self.var_comps.get() or '—'],
                 [self.t('fot_coef_color'), self.var_color_t.get()], [self.t('fot_colores'), self.var_colores.get() or '—']]
        for ruta, foto in sorted(self.fotos.items(), key=lambda par: par[1]['n']):
            if foto.get('rms_astrometrico') is not None:
                filas.append([self.t('conf_astrometria', n=foto['n']),
                              self.t('conf_rms', rms=numero(foto['rms_astrometrico'], 3),
                                     ajustada=self.t('si') if foto.get('wcs_ajustada') else self.t('no'))])
        for codigo in ('SAT', 'CON', 'BOR', 'SNR', 'FWHM', 'ELO'):
            filas.append([codigo, self.t('calidad_' + codigo)])
        return [(self.t('conf_parametro'), None, 42), (self.t('conf_valor'), None, 70)], filas

    def escribir_xlsx(self, ruta):
        hojas = [(self.t('hoja_estrellas'), *self.tabla_por_estrellas()),
                 (self.t('hoja_todas'), *self.tabla_todas())]
        if self.revision:
            hojas.append((self.t('hoja_revision'), *self.tabla_revision()))
        hojas.append((self.t('hoja_configuracion'), *self.tabla_configuracion()))
        if self.historial:
            hojas.append((self.t('hist_titulo'), *self.tabla_historial()))
        escribir_xlsx(ruta, hojas)

    def escribir_csv(self, ruta):
        # Separador ";" y coma decimal: así Excel configurado en catalán o español lo abre bien.
        columnas, filas = self.tabla_por_estrellas()
        with open(ruta, 'w', newline='', encoding='utf-8-sig') as archivo:
            escritor = csv.writer(archivo, delimiter=';')
            escritor.writerow([c[0] for c in columnas])
            for fila in filas:
                escritor.writerow(['' if valor is None else
                                   (texto_seguro_csv(valor) if decimales is None else numero(valor, decimales))
                                   for valor, (_t, decimales, _a) in zip(fila, columnas)])

    # =================================================================== HERRAMIENTAS DE LA 1.5
    # ------------------------------------------------------------------ utilidades de ventana
    def _ventana(self, titulo, nombre):
        """Ventana secundaria con el estilo del programa (si ya estaba abierta, se rehace)."""
        anterior = getattr(self, nombre, None)
        try:
            if isinstance(anterior, tk.Misc) and anterior.winfo_exists():
                anterior.destroy()
        except tk.TclError:
            pass
        v = tk.Toplevel(self.root)
        v.title(titulo)
        v.configure(background=FONDO)
        v.transient(self.root)
        setattr(self, nombre, v)
        marco = ttk.Frame(v, padding=12)
        marco.pack(fill=tk.BOTH, expand=True)
        return v, marco

    def _centrar(self, v):
        v.update_idletasks()
        x = self.root.winfo_rootx() + max(0, (self.root.winfo_width() - v.winfo_width()) // 2)
        y = self.root.winfo_rooty() + max(0, (self.root.winfo_height() - v.winfo_height()) // 4)
        v.geometry(f'+{x}+{y}')

    def _tabla_simple(self, padre, columnas, alto=8, seleccion=tk.EXTENDED):
        """Tabla (Treeview) con barra de desplazamiento. columnas: [(título, ancho, alineación)]."""
        marco = ttk.Frame(padre)
        tabla = ttk.Treeview(marco, columns=[str(i) for i in range(len(columnas))], show='headings',
                             height=alto, selectmode=seleccion)
        for i, (titulo, ancho, alineacion) in enumerate(columnas):
            tabla.heading(str(i), text=titulo)
            tabla.column(str(i), width=ancho, minwidth=30, anchor=alineacion, stretch=True)
        barra = ttk.Scrollbar(marco, orient=tk.VERTICAL, command=tabla.yview)
        tabla.configure(yscrollcommand=barra.set)
        barra.pack(side=tk.RIGHT, fill=tk.Y)
        tabla.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        return marco, tabla

    def _grafico(self, padre, ancho=7.6, alto=3.6, filas=1):
        """Gráfico de matplotlib dentro de una ventana, con los colores del programa y la
        barra de herramientas (zoom, desplazar, guardar como imagen). matplotlib solo se
        carga aquí, la primera vez que hace falta un gráfico."""
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
        figura = Figure(figsize=(ancho, alto), dpi=96, facecolor=FONDO)
        ejes = [figura.add_subplot(filas, 1, i + 1) for i in range(filas)]
        for ax in ejes:
            self._estilo_ejes(ax)
        contenedor = ttk.Frame(padre)
        lienzo = FigureCanvasTkAgg(figura, master=contenedor)
        barra = NavigationToolbar2Tk(lienzo, contenedor, pack_toolbar=False)
        barra.update()
        barra.pack(side=tk.BOTTOM, fill=tk.X)
        lienzo.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        return contenedor, figura, (ejes if filas > 1 else ejes[0]), lienzo

    @staticmethod
    def _estilo_ejes(ax, enteros=False):
        """Colores del programa en unos ejes de matplotlib (hay que volver a ponerlos
        después de ax.clear(), que deja los de por defecto)."""
        ax.set_facecolor(CAMPO)
        ax.tick_params(colors=TEXTO_SUAVE, labelsize=8)
        for borde in ax.spines.values():
            borde.set_color(BORDE)
        ax.grid(True, color=BORDE, alpha=0.6, linewidth=0.6)
        ax.xaxis.label.set_color(TEXTO_SUAVE)
        ax.yaxis.label.set_color(TEXTO_SUAVE)
        ax.title.set_color(TEXTO)
        ax.ticklabel_format(axis='both', style='plain', useOffset=False)   # 389,37 y no "+3,893e2 + 0,07"
        if enteros:
            from matplotlib.ticker import MaxNLocator
            ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    def _necesita_foto(self):
        if self.original is None:
            messagebox.showinfo(TITULO_VENTANA, self.t('necesita_foto'))
            return True
        return False

    def _letras_medidas(self):
        return sorted({m['estrella'] for m in self.mediciones}, key=LETRAS_ESTRELLA.index)

    def _eje_tiempo(self, rutas):
        """Eje X de los gráficos: horas desde la primera foto (si todas tienen hora) o número de foto."""
        jds = [self.fotos[r].get('jd') for r in rutas]
        if jds and all(jds):
            jd0 = min(jds)
            return {r: (self.fotos[r]['jd'] - jd0) * 24 for r in rutas}, self.t('eje_horas', fecha=fecha_desde_jd(jd0)[:16])
        return {r: self.fotos[r]['n'] for r in rutas}, self.t('eje_foto')

    # ------------------------------------------------------------------ panel 5: análisis
    def _panel_analisis(self, panel):
        caja = ttk.LabelFrame(panel, text=self.t('paso5'), padding=8)
        caja.pack(fill=tk.X, pady=(0, 8))
        self.caja_analisis = caja
        caja.columnconfigure((0, 1), weight=1, uniform='botones')
        botones = (('boton_guia', self.ventana_guia, 'guia_explicacion'),
                   ('boton_detectar', self.detectar_estrellas, 'ayuda_detectar'),
                   ('boton_dobles', self.ventana_dobles, 'ayuda_dobles'),
                   ('boton_moviles', self.ventana_moviles, 'ayuda_moviles'),
                   ('boton_serie', self.revisar_serie, 'ayuda_serie_rev'),
                   ('boton_catalogos', self.ventana_catalogos, 'ayuda_catalogos'),
                   ('boton_wcs_manual', self.ventana_wcs_manual, 'ayuda_wcs_manual'),
                   ('boton_info_foto', self.ventana_info_foto, 'ayuda_info_foto'),
                   ('boton_recortes', self.ventana_recortes, 'ayuda_recortes'),
                   ('boton_historial', self.ventana_historial, 'ayuda_historial'),
                   ('boton_guardar_sesion', self.guardar_sesion, 'ayuda_sesion'),
                   ('boton_abrir_sesion', self.abrir_sesion, 'ayuda_sesion'))
        for i, (clave, orden, ayuda) in enumerate(botones):
            b = ttk.Button(caja, text=self.t(clave), command=orden, style='Pequeno.TButton')
            b.grid(row=i // 2, column=i % 2, sticky=tk.EW, padx=2, pady=2)
            Ayuda(b, self.t(ayuda))

    # ------------------------------------------------------------------ detectar estrellas
    def detectar_estrellas(self):
        """Busca todas las estrellas de la foto y las marca (otra pulsación las quita)."""
        if self._necesita_foto():
            return
        if self.fuentes is not None:
            self.fuentes = None
            self._mensaje('')
            self._programar_dibujo()
            return
        self.root.config(cursor='watch')
        self.root.update_idletasks()
        try:
            self.fuentes = C.detectar_fuentes(self.original, 5.0, 3000, saturacion=self.fotos[self.ruta].get('saturacion'))
        finally:
            self.root.config(cursor='')
        estadisticas = [f['fwhm'] for f in self.fuentes[:50] if not f['saturada']]
        self._mensaje(self.t('detectadas', n=len(self.fuentes),
                             fwhm=numero(float(np.median(estadisticas)), 2) if estadisticas else '—'))
        self._programar_dibujo()

    def _dibujar_detecciones(self):
        """Marcas de las estrellas detectadas y de los objetos de catálogo en la foto."""
        if self.ruta is None or not self.var_marcas.get():
            return
        foto = self.fotos[self.ruta]
        ancho_v, alto_v = self.lienzo.winfo_width(), self.lienzo.winfo_height()
        for f in (self.fuentes or [])[:3000]:
            x, y = transformar_punto(f['x'], f['y'], foto['ancho'], foto['alto'], self.angulo, self.espejo)
            u, v = self.imagen_a_lienzo(x, y)
            if -10 <= u <= ancho_v + 10 and -10 <= v <= alto_v + 10:
                r = max(f['fwhm'] * self.zoom, 4.0)
                self.lienzo.create_oval(u - r, v - r, u + r, v + r, outline='#ff5c5c' if f['saturada'] else '#4fd1c5')
        for marca in foto.get('marcas_catalogo', []):
            x, y = transformar_punto(marca['x'], marca['y'], foto['ancho'], foto['alto'], self.angulo, self.espejo)
            u, v = self.imagen_a_lienzo(x, y)
            if -40 <= u <= ancho_v + 40 and -40 <= v <= alto_v + 40:
                color = marca.get('color', '#ff9cf0')
                self.lienzo.create_line(u - 9, v, u - 3, v, fill=color, width=2)
                self.lienzo.create_line(u + 3, v, u + 9, v, fill=color, width=2)
                self.lienzo.create_line(u, v - 9, u, v - 3, fill=color, width=2)
                self.lienzo.create_line(u, v + 3, u, v + 9, fill=color, width=2)
                self.lienzo.create_text(u + 10, v + 10, text=marca['texto'], fill=color, anchor=tk.NW,
                                        font=('Segoe UI', 9, 'bold'))

    # ------------------------------------------------------------------ fotometría y curvas de luz
    def _serie_fotometrica(self):
        """{letra: [(ruta, magnitud, error)]} con la magnitud calibrada/relativa (si hay
        comparaciones) o la instrumental."""
        activa = self._fotometria_activa()
        series = {}
        for m in sorted(self.mediciones, key=lambda m: m['foto_n']):
            mag, err = self.magnitud_calibrada(m) if activa else (m.get('mag'), m.get('sigma_mag'))
            if mag is not None:
                series.setdefault(m['estrella'], []).append((m['ruta'], mag, err))
        return series

    def ventana_fotometria(self):
        if not self.mediciones:
            messagebox.showinfo(TITULO_VENTANA, self.t('sin_medidas'))
            return
        v, marco = self._ventana(self.t('fot_titulo'), 'ventana_fot')
        ttk.Label(marco, text=self.t('fot_explicacion'), wraplength=860, style='Gris.TLabel', justify=tk.LEFT).pack(fill=tk.X)
        ajustes = ttk.Frame(marco)
        ajustes.pack(fill=tk.X, pady=(8, 4))
        for i, (clave, var, ancho, ayuda) in enumerate((('comparaciones', self.var_comps, 40, 'ayuda_comparaciones'),
                                                         ('fot_coef_color', self.var_color_t, 8, 'ayuda_coef_color'),
                                                         ('fot_colores', self.var_colores, 26, 'ayuda_colores'))):
            ttk.Label(ajustes, text=self.t(clave)).grid(row=0, column=2 * i, sticky=tk.W, padx=(0 if i == 0 else 12, 4))
            e = ttk.Entry(ajustes, textvariable=var, width=ancho)
            e.grid(row=0, column=2 * i + 1, sticky=tk.W)
            Ayuda(e, self.t(ayuda))
        estado = ttk.Label(marco, text='', style='Gris.TLabel', wraplength=860, justify=tk.LEFT)
        estado.pack(fill=tk.X)
        acciones = ttk.Frame(marco)
        acciones.pack(fill=tk.X, pady=(4, 6))
        bandas = ['V (APASS)', 'B (APASS)', "g' (APASS)", "r' (APASS)", "i' (APASS)", 'G (Gaia DR3)']
        var_banda = tk.StringVar(value=self.ajustes.get('banda_catalogo', bandas[0]))
        columnas = [(self.t('col_estrella'), 60, tk.CENTER), (self.t('fot_n'), 55, tk.CENTER),
                    (self.t('fot_media'), 90, tk.E), (self.t('fot_desv'), 80, tk.E), (self.t('fot_amplitud'), 80, tk.E),
                    (self.t('fot_chi2'), 70, tk.E), (self.t('fot_diagnostico'), 300, tk.W)]
        marco_tabla, tabla = self._tabla_simple(marco, columnas, alto=7)
        marco_tabla.pack(fill=tk.X)
        contenedor, figura, ax, lienzo = self._grafico(marco, 8.6, 3.4)
        contenedor.pack(fill=tk.BOTH, expand=True, pady=(8, 0))

        def recalcular(*_):
            series = self._serie_fotometrica()
            comps = {c[0] for c in self._comparaciones() or []}
            estado.config(text=self.lbl_fotometria.cget('text') or self.t('fot_sin_comp'))
            tabla.delete(*tabla.get_children())
            for letra in sorted(series, key=LETRAS_ESTRELLA.index):
                datos = series[letra]
                info = C.analisis_variabilidad([d[1] for d in datos], [d[2] for d in datos])
                saturada = any('SAT' in (m.get('calidad') or []) for m in self.mediciones if m['estrella'] == letra)
                if letra in comps:
                    diagnostico = self.t('fot_es_comp')
                elif saturada:
                    diagnostico = self.t('fot_saturada')
                elif info['n'] < 3:
                    diagnostico = self.t('fot_pocas')
                elif info['variable']:
                    diagnostico = self.t('fot_variable')
                else:
                    diagnostico = self.t('fot_constante')
                tabla.insert('', tk.END, iid=letra, values=(letra, info['n'], numero(info['media'], 3),
                                                             numero(info['desviacion'], 3), numero(info['amplitud'], 3),
                                                             numero(info['chi2_red'], 1), diagnostico))
            dibujar()

        def dibujar(*_):
            series = self._serie_fotometrica()
            comps = {c[0] for c in self._comparaciones() or []}
            elegidas = list(tabla.selection()) or [l for l in sorted(series, key=LETRAS_ESTRELLA.index)
                                                   if l not in comps][:6]
            rutas = sorted({r for l in elegidas for r, _m, _e in series.get(l, [])}, key=lambda r: self.fotos[r]['n'])
            eje, titulo_x = self._eje_tiempo(rutas)
            ax.clear()
            self._estilo_ejes(ax)
            colores = ['#ffd27a', '#86c8ff', '#7be0a0', '#ff9cf0', '#ffb36b', '#c9b6ff']
            for i, letra in enumerate(elegidas):
                datos = series.get(letra, [])
                if not datos:
                    continue
                ax.errorbar([eje[r] for r, _m, _e in datos], [m for _r, m, _e in datos],
                            yerr=[e or 0 for _r, _m, e in datos], fmt='o-', ms=4, lw=1, capsize=2,
                            color=colores[i % len(colores)], label=letra)
            ax.invert_yaxis()
            ax.set_xlabel(titulo_x)
            ax.set_ylabel(self.tabla.heading('mag')['text'])
            if elegidas:
                leyenda = ax.legend(fontsize=8, facecolor=FONDO, edgecolor=BORDE)
                for texto in leyenda.get_texts():
                    texto.set_color(TEXTO)
            figura.tight_layout()
            lienzo.draw_idle()

        def automaticas():
            elegidas = self.comparaciones_automaticas(excluir=tuple(tabla.selection()))
            if len(elegidas) < 2:
                messagebox.showinfo(self.t('fot_titulo'), self.t('fot_auto_pocas'), parent=v)
                return
            conocidas = {c[0]: c for c in self._comparaciones() or []}
            self.var_comps.set(' '.join(self._formato_comp(conocidas.get(l, (l, None, None))) for l in elegidas))
            recalcular()

        def del_catalogo():
            self.ajustes['banda_catalogo'] = var_banda.get()
            mensaje = self.magnitudes_de_catalogo(bandas.index(var_banda.get()))
            messagebox.showinfo(self.t('fot_titulo'), mensaje, parent=v)
            recalcular()

        def exportar():
            series = self._serie_fotometrica()
            ruta = filedialog.asksaveasfilename(parent=v, defaultextension='.csv', initialfile=self.t('fot_archivo'),
                                                filetypes=[('CSV', '*.csv')])
            if not ruta:
                return
            letras = sorted(series, key=LETRAS_ESTRELLA.index)
            rutas = sorted({r for l in letras for r, _m, _e in series[l]}, key=lambda r: self.fotos[r]['n'])
            try:
                with open(ruta, 'w', newline='', encoding='utf-8-sig') as archivo:
                    escritor = csv.writer(archivo, delimiter=';')
                    escritor.writerow([self.t('col_foto_n_largo'), self.t('col_foto'), self.t('col_jd'), 'HJD', 'BJD_TDB'] +
                                      [x for l in letras for x in (f'{l}_Mag', f'{l}_σMag')])
                    for r in rutas:
                        hjd, bjd = self._hjd_foto(r), self._bjd_foto(r)
                        fila = [self.fotos[r]['n'], texto_seguro_csv(self.fotos[r]['archivo']),
                                numero(self.fotos[r].get('jd'), 6) if self.fotos[r].get('jd') else '',
                                numero(hjd, 6) if hjd else '', numero(bjd, 6) if bjd else '']
                        for l in letras:
                            dato = next(((m, e) for rr, m, e in series[l] if rr == r), (None, None))
                            fila += [numero(dato[0], 4) if dato[0] is not None else '',
                                     numero(dato[1], 4) if dato[1] is not None else '']
                        escritor.writerow(fila)
            except OSError as error:
                messagebox.showerror(self.t('error_guardar'), str(error), parent=v)
                return
            messagebox.showinfo(self.t('guardado_t'), self.t('mpc_guardado', ruta=ruta), parent=v)

        ttk.Button(acciones, text=self.t('fot_actualizar'), command=recalcular).pack(side=tk.LEFT)
        ttk.Button(acciones, text=self.t('fot_auto'), command=automaticas).pack(side=tk.LEFT, padx=6)
        ttk.Label(acciones, text=self.t('fot_banda')).pack(side=tk.LEFT, padx=(12, 4))
        ttk.Combobox(acciones, textvariable=var_banda, values=bandas, state='readonly', width=14).pack(side=tk.LEFT)
        ttk.Button(acciones, text=self.t('fot_catalogo'), command=del_catalogo).pack(side=tk.LEFT, padx=6)
        ttk.Button(acciones, text=self.t('fot_exportar'), command=exportar).pack(side=tk.RIGHT)
        mas = ttk.Frame(marco)
        mas.pack(fill=tk.X, pady=(0, 6), before=marco_tabla)
        for clave, orden, ayuda in (('cal_boton', lambda: self.ventana_calibracion(var_banda.get()), 'ayuda_cal'),
                                    ('ext_boton', lambda: self.extincion_serie(var_banda.get()), 'ayuda_ext'),
                                    ('ropt_boton', self.radio_optimo, 'ayuda_ropt'),
                                    ('per_boton', lambda: self.ventana_periodo((tabla.selection() or [None])[0]), 'ayuda_per'),
                                    ('min_boton', self.ventana_minimo, 'ayuda_min')):
            b = ttk.Button(mas, text=self.t(clave), command=orden)
            b.pack(side=tk.LEFT, padx=(0, 6))
            Ayuda(b, self.t(ayuda))
        tabla.bind('<<TreeviewSelect>>', dibujar)
        vigilancia = self.var_comps.trace_add('write', lambda *_: v.after_idle(recalcular) if v.winfo_exists() else None)
        v.bind('<Destroy>', lambda e: self.var_comps.trace_remove('write', vigilancia) if e.widget is v else None)
        recalcular()
        self._centrar(v)

    def comparaciones_automaticas(self, excluir=(), maximo=4):
        """Las mejores estrellas de comparación (ver ciencia.elegir_comparaciones). No se
        usan las saturadas ni las contaminadas por una vecina en ninguna foto."""
        malas = {m['estrella'] for m in self.mediciones if {'SAT', 'CON'} & set(m.get('calidad') or [])}
        series = {}
        for m in self.mediciones:
            if m.get('mag') is not None and m['estrella'] not in malas:
                series.setdefault(m['estrella'], {})[m['foto_n']] = (m['mag'], m['sigma_mag'])
        return C.elegir_comparaciones(series, excluir=tuple(excluir), maximo=maximo)

    @staticmethod
    def _formato_comp(comp):
        letra, mag, color = comp
        texto = letra
        if mag is not None:
            texto += '=' + numero(mag, 3)
            if color is not None:
                texto += '/' + numero(color, 3)
        return texto

    def magnitudes_de_catalogo(self, banda):
        """Pone a las estrellas de comparación su magnitud (y su índice de color) del
        catálogo guardado en local (APASS o Gaia). No se conecta a internet."""
        catalogo, clave, color = [('apass', 'V', ('B', 'V')), ('apass', 'B', ('B', 'V')), ('apass', 'g', ('g', 'r')),
                                  ('apass', 'r', ('g', 'r')), ('apass', 'i', ('r', 'i')), ('gaia', 'mag', None)][banda]
        comps = self._comparaciones()
        if not comps:
            return self.t('fot_cat_sin_comps')
        nuevas, faltan = [], []
        for letra, mag, ic in comps:
            m = next((m for m in self.mediciones if m['estrella'] == letra and m.get('ar') is not None), None)
            if m is None:
                faltan.append(letra)
                nuevas.append((letra, mag, ic))
                continue
            datos = catalogos.desde_cache(catalogo, m['ar'], m['dec'], 10 / 3600)
            if datos is None:
                return self.t('fot_cat_sin_cache', cat=catalogos.NOMBRES[catalogo])
            filas = [f for f in datos[0] if f.get(clave) is not None]
            idx, _d = C.emparejar([m['ar']], [m['dec']], [f['ar'] for f in filas], [f['dec'] for f in filas], 3.0)
            if idx[0] < 0:
                faltan.append(letra)
                nuevas.append((letra, mag, ic))
                continue
            fila = filas[idx[0]]
            if color is None:
                indice = fila.get('color')
            else:
                indice = (fila[color[0]] - fila[color[1]]) if fila.get(color[0]) is not None and fila.get(color[1]) is not None else None
            nuevas.append((letra, fila[clave], indice))
        self.var_comps.set(' '.join(self._formato_comp(c) for c in nuevas))
        texto = self.t('fot_cat_ok', n=len(comps) - len(faltan), cat=catalogos.NOMBRES[catalogo])
        if faltan:
            texto += '\n' + self.t('fot_cat_faltan', letras=', '.join(faltan))
        return texto

    # ------------------------------------------------------------------ estrellas dobles
    def ventana_dobles(self):
        letras = self._letras_medidas()
        if len(letras) < 2:
            messagebox.showinfo(TITULO_VENTANA, self.t('dobles_pocas'))
            return
        v, marco = self._ventana(self.t('dobles_titulo'), 'v_dobles')
        ttk.Label(marco, text=self.t('dobles_explicacion'), wraplength=760, style='Gris.TLabel', justify=tk.LEFT).pack(fill=tk.X)
        fila = ttk.Frame(marco)
        fila.pack(fill=tk.X, pady=8)
        var_a, var_b = tk.StringVar(value=letras[0]), tk.StringVar(value=letras[1])
        for texto, var in ((self.t('dobles_principal'), var_a), (self.t('dobles_secundaria'), var_b)):
            ttk.Label(fila, text=texto).pack(side=tk.LEFT, padx=(0, 4))
            ttk.Combobox(fila, textvariable=var, values=letras, state='readonly', width=4).pack(side=tk.LEFT, padx=(0, 14))
        resumen = ttk.Label(marco, text='', font=('Consolas', 11, 'bold'), foreground=DORADO, justify=tk.LEFT)
        resumen.pack(fill=tk.X)
        columnas = [(self.t('col_foto_n'), 60, tk.CENTER), (self.t('col_foto'), 200, tk.W), (self.t('dobles_sep_px'), 110, tk.E),
                    (self.t('dobles_sep_arc'), 110, tk.E), (self.t('dobles_pa'), 100, tk.E)]
        marco_tabla, tabla = self._tabla_simple(marco, columnas, alto=7)
        marco_tabla.pack(fill=tk.X, pady=(6, 0))
        contenedor, figura, (ax1, ax2), lienzo = self._grafico(marco, 7.6, 4.0, filas=2)
        contenedor.pack(fill=tk.BOTH, expand=True, pady=(8, 0))
        self._dobles_resultado = []

        def calcular(*_):
            a, b = var_a.get(), var_b.get()
            tabla.delete(*tabla.get_children())
            filas = []
            for ruta in sorted({m['ruta'] for m in self.mediciones}, key=lambda r: self.fotos[r]['n']):
                ma = next((m for m in self.mediciones if m['ruta'] == ruta and m['estrella'] == a), None)
                mb = next((m for m in self.mediciones if m['ruta'] == ruta and m['estrella'] == b), None)
                if ma is None or mb is None or a == b:
                    continue
                (xa, ya), (xb, yb) = self._original_xy(ma), self._original_xy(mb)
                sep_px = math.hypot(xb - xa, yb - ya)
                sep_arc = pa = None
                if ma.get('ar') is not None and mb.get('ar') is not None:
                    sep_arc, pa = C.separacion_y_angulo(ma['ar'], ma['dec'], mb['ar'], mb['dec'])
                filas.append((ruta, sep_px, sep_arc, pa))
                tabla.insert('', tk.END, values=(self.fotos[ruta]['n'], self.fotos[ruta]['archivo'], numero(sep_px, 3),
                                                 numero(sep_arc, 3) if sep_arc is not None else '—',
                                                 numero(pa, 2) if pa is not None else '—'))
            self._dobles_resultado = [(a, b, filas)]
            texto = self.t('dobles_n', n=len(filas))
            media_px = C.media_y_error([f[1] for f in filas])
            if media_px[0] is not None:
                texto += '\n' + self.t('dobles_media_px', m=numero(media_px[0], 3), e=numero(media_px[2], 3))
            media_arc = C.media_y_error([f[2] for f in filas])
            media_pa = C.media_y_error([f[3] for f in filas], circular=True)
            if media_arc[0] is not None:
                texto += '\n' + self.t('dobles_media_arc', m=numero(media_arc[0], 3), e=numero(media_arc[2], 3),
                                       pa=numero(media_pa[0], 2), epa=numero(media_pa[2], 2))
            elif filas:
                texto += '\n' + self.t('dobles_sin_wcs')
            resumen.config(text=texto)
            eje, titulo_x = self._eje_tiempo([f[0] for f in filas])
            for ax in (ax1, ax2):
                ax.clear()
                self._estilo_ejes(ax)
            con_arc = [f for f in filas if f[2] is not None]
            if con_arc:
                ax1.plot([eje[f[0]] for f in con_arc], [f[2] for f in con_arc], 'o-', color='#ffd27a', ms=4)
                ax1.set_ylabel(self.t('dobles_sep_arc'))
                ax2.plot([eje[f[0]] for f in con_arc], [f[3] for f in con_arc], 'o-', color='#86c8ff', ms=4)
                ax2.set_ylabel(self.t('dobles_pa'))
            else:
                ax1.plot([eje[f[0]] for f in filas], [f[1] for f in filas], 'o-', color='#ffd27a', ms=4)
                ax1.set_ylabel(self.t('dobles_sep_px'))
            ax2.set_xlabel(titulo_x)
            figura.tight_layout()
            lienzo.draw_idle()

        ttk.Button(fila, text=self.t('dobles_calcular'), style='Grande.TButton', command=calcular).pack(side=tk.LEFT)
        for var in (var_a, var_b):
            var.trace_add('write', calcular)
        calcular()
        self._centrar(v)

    # ------------------------------------------------------------------ objetos móviles
    def analisis_moviles(self):
        """Para cada letra medida en 2 fotos o más: ¿se mueve respecto a las estrellas?
        Con AR/Dec: velocidad en ″/hora. Sin ellas: en píxeles respecto a la estrella de
        referencia (la más brillante medida en todas las fotos)."""
        resultado = []
        letras = self._letras_medidas()
        rutas_todas = {m['ruta'] for m in self.mediciones}
        candidatas = [l for l in letras if {m['ruta'] for m in self.mediciones if m['estrella'] == l} >= rutas_todas]
        for letra in letras:
            medidas = sorted((m for m in self.mediciones if m['estrella'] == letra),
                             key=lambda m: (self.fotos[m['ruta']].get('jd') or 0, m['foto_n']))
            if len(medidas) < 2:
                continue
            jds = [self.fotos[m['ruta']].get('jd') for m in medidas]
            info = {'letra': letra, 'n': len(medidas), 'medidas': medidas}
            if all(jds) and all(m.get('ar') is not None for m in medidas):
                ars, decs = self._cielo_respecto_al_campo(letra, medidas)
                mov = C.movimiento_en_cielo(jds, ars, decs)
                if mov:
                    escala = np.median([self.fotos[m['ruta']]['wcs'].escala for m in medidas])
                    error = float(np.median([math.hypot(m['sigma_x'], m['sigma_y']) for m in medidas])) * escala
                    total = mov[0] * (max(jds) - min(jds)) * 24
                    info.update(unidad='″', velocidad=mov[0], angulo=mov[1], rms=mov[2], total=total,
                                se_mueve=total > max(3 * error * math.sqrt(2), 0.5))
            else:
                otras = [l for l in candidatas if l != letra]
                if otras:
                    ref = max(otras, key=lambda l: np.median([m['flujo'] for m in self.mediciones if m['estrella'] == l]))
                    puntos = []
                    for m in medidas:
                        r = next((o for o in self.mediciones if o['ruta'] == m['ruta'] and o['estrella'] == ref), None)
                        if r is not None:
                            (xm, ym), (xr, yr) = self._original_xy(m), self._original_xy(r)
                            puntos.append((self.fotos[m['ruta']].get('jd') or m['foto_n'], xm - xr, ym - yr,
                                           math.hypot(m['sigma_x'], m['sigma_y'])))
                    if len(puntos) >= 2:
                        t, xs, ys, es = zip(*puntos)
                        total = math.hypot(xs[-1] - xs[0], ys[-1] - ys[0])
                        ajuste = C.ajuste_lineal_movimiento(t, xs, ys)
                        horas = all(jds)
                        info.update(unidad='px', referencia=ref, total=total, rms=ajuste[5] if ajuste else None,
                                    velocidad=(math.hypot(ajuste[2], ajuste[3]) / 24 if horas else None) if ajuste else None,
                                    angulo=None, se_mueve=C.parece_movil(t, xs, ys, list(es)))
            resultado.append(info)
        return resultado

    def _cielo_respecto_al_campo(self, letra, medidas):
        """AR/Dec de las medidas de 'letra' corregidas foto a foto con las demás estrellas
        medidas: si la solución astrométrica de una foto está un poco desplazada, todas sus
        estrellas "saltan" lo mismo, y ese salto se resta. Así un error del WCS de una foto
        no parece un movimiento. (Sin otras estrellas, se devuelven tal cual.)"""
        ar0, dec0 = medidas[0]['ar'], medidas[0]['dec']
        otras = {}
        for m in self.mediciones:
            if m['estrella'] != letra and m.get('ar') is not None:
                otras.setdefault(m['estrella'], []).append(m)
        desvios = {}   # ruta -> lista de (dξ, dη) de las otras estrellas respecto a su posición media
        for lista in otras.values():
            if len(lista) < 2:
                continue
            xi, eta = C.proyectar(np.array([m['ar'] for m in lista]), np.array([m['dec'] for m in lista]), ar0, dec0)
            for m, x, y in zip(lista, xi - xi.mean(), eta - eta.mean()):
                desvios.setdefault(m['ruta'], []).append((x, y))
        xi, eta = C.proyectar(np.array([m['ar'] for m in medidas]), np.array([m['dec'] for m in medidas]), ar0, dec0)
        for i, m in enumerate(medidas):
            if desvios.get(m['ruta']):
                d = np.median(np.array(desvios[m['ruta']]), axis=0)
                xi[i] -= d[0]
                eta[i] -= d[1]
        ars, decs = C.desproyectar(xi, eta, ar0, dec0)
        return list(np.atleast_1d(ars)), list(np.atleast_1d(decs))

    def ventana_moviles(self):
        datos = self.analisis_moviles()
        if not datos:
            messagebox.showinfo(TITULO_VENTANA, self.t('moviles_pocas'))
            return
        v, marco = self._ventana(self.t('moviles_titulo'), 'v_moviles')
        ttk.Label(marco, text=self.t('moviles_explicacion'), wraplength=820, style='Gris.TLabel', justify=tk.LEFT).pack(fill=tk.X)
        columnas = [(self.t('col_estrella'), 60, tk.CENTER), (self.t('fot_n'), 50, tk.CENTER),
                    (self.t('moviles_total'), 120, tk.E), (self.t('moviles_velocidad'), 130, tk.E),
                    (self.t('moviles_angulo'), 110, tk.E), (self.t('moviles_rms'), 100, tk.E),
                    (self.t('moviles_estado'), 200, tk.W)]
        marco_tabla, tabla = self._tabla_simple(marco, columnas, alto=6, seleccion=tk.BROWSE)
        marco_tabla.pack(fill=tk.X, pady=(8, 0))
        for d in datos:
            u = d.get('unidad', '')
            tabla.insert('', tk.END, iid=d['letra'], values=(
                d['letra'], d['n'], (numero(d['total'], 2) + ' ' + u) if d.get('total') is not None else '—',
                (numero(d['velocidad'], 2) + f' {u}/h') if d.get('velocidad') is not None else '—',
                numero(d['angulo'], 1) + '°' if d.get('angulo') is not None else '—',
                (numero(d['rms'], 3) + ' ' + u) if d.get('rms') is not None else '—',
                self.t('moviles_se_mueve') if d.get('se_mueve') else self.t('moviles_fija')))
        contenedor, figura, ax, lienzo = self._grafico(marco, 7.6, 3.8)
        contenedor.pack(fill=tk.BOTH, expand=True, pady=(8, 0))
        info_sky = ttk.Label(marco, text='', style='Gris.TLabel', wraplength=820, justify=tk.LEFT)
        info_sky.pack(fill=tk.X, pady=(4, 0))

        def dibujar(*_):
            seleccion = tabla.selection()
            d = next((d for d in datos if seleccion and d['letra'] == seleccion[0]), None)
            if d is None:
                d = next((d for d in datos if d.get('se_mueve')), datos[0])
            ax.clear()
            self._estilo_ejes(ax)
            medidas = d['medidas']
            if d.get('unidad') == '″':
                ar0, dec0 = medidas[0]['ar'], medidas[0]['dec']
                xi, eta = C.proyectar(np.array([m['ar'] for m in medidas]), np.array([m['dec'] for m in medidas]), ar0, dec0)
                xs, ys = -xi * 3600, eta * 3600       # este a la izquierda, como en el cielo
                ax.set_xlabel(self.t('moviles_eje_este'))
                ax.set_ylabel(self.t('moviles_eje_norte'))
            else:
                xs = [self._original_xy(m)[0] for m in medidas]
                ys = [self._original_xy(m)[1] for m in medidas]
                ax.set_xlabel('X (px)')
                ax.set_ylabel('Y (px)')
            ax.plot(xs, ys, 'o-', color='#ffd27a', ms=5)
            for x, y, m in zip(xs, ys, medidas):
                ax.annotate(str(m['foto_n']), (x, y), color=TEXTO, fontsize=8, xytext=(4, 4), textcoords='offset points')
            ax.set_title(self.t('moviles_trayectoria', letra=d['letra']), fontsize=10)
            ax.set_aspect('equal', adjustable='datalim')
            figura.tight_layout()
            lienzo.draw_idle()

        def skybot():
            seleccion = tabla.selection()
            d = next((d for d in datos if seleccion and d['letra'] == seleccion[0]), None)
            if d is None:
                messagebox.showinfo(self.t('moviles_titulo'), self.t('moviles_elige'), parent=v)
                return
            m = next((m for m in d['medidas'] if m.get('ar') is not None and self.fotos[m['ruta']].get('jd')), None)
            if m is None:
                messagebox.showinfo(self.t('moviles_titulo'), self.t('moviles_skybot_sin_wcs'), parent=v)
                return
            jd = self.fotos[m['ruta']]['jd']

            def mostrar(filas):
                cerca = sorted(((C.separacion_y_angulo(m['ar'], m['dec'], f['ar'], f['dec'])[0], f) for f in filas),
                               key=lambda par: par[0])[:5]
                if not cerca:
                    info_sky.config(text=self.t('moviles_skybot_nada'))
                    return
                lineas = [self.t('moviles_skybot_cabecera', letra=d['letra'])]
                for distancia, f in cerca:
                    lineas.append(self.t('moviles_skybot_fila', nombre=f['id'], clase=f.get('clase', ''),
                                         mag=numero(f.get('mag'), 1), d=numero(distancia, 1)))
                info_sky.config(text='\n'.join(lineas))
            guardado = catalogos.desde_cache('skybot', m['ar'], m['dec'], 0.05, jd=jd)
            if guardado is not None:
                mostrar(guardado[0])
                return
            self._preguntar_y_descargar('skybot', m['ar'], m['dec'], 0.1, jd=jd, padre=v,
                                        al_terminar=lambda filas: mostrar(filas))

        botones = ttk.Frame(marco)
        botones.pack(fill=tk.X, pady=(6, 0))
        ttk.Button(botones, text=self.t('moviles_skybot'), command=skybot).pack(side=tk.LEFT)
        tabla.bind('<<TreeviewSelect>>', dibujar)
        dibujar()
        self._centrar(v)

    # ------------------------------------------------------------------ revisar la serie
    def revisar_serie(self):
        """Mira todas las fotos de la serie abierta: fondo, ruido, número de estrellas,
        FWHM, elongación y desplazamiento respecto a la primera. Marca las fotos raras
        (desenfocadas, movidas, con nubes o luz de la Luna...)."""
        rutas = list(self.serie) if len(self.serie) > 1 else sorted(self.fotos, key=lambda r: self.fotos[r]['n'])
        rutas = [r for r in rutas if os.path.exists(separar_ruta(r)[0])]
        if len(rutas) < 2:
            messagebox.showinfo(TITULO_VENTANA, self.t('serie_pocas'))
            return
        v, marco = self._ventana(self.t('serie_titulo'), 'ventana_serie')
        progreso = ttk.Label(marco, text='', style='Gris.TLabel')
        progreso.pack(fill=tk.X)
        barra = ttk.Progressbar(marco, maximum=len(rutas), style='Emillium.Horizontal.TProgressbar')
        barra.pack(fill=tk.X, pady=4)
        resultados, primeras = [], None
        for i, ruta in enumerate(rutas):
            try:
                progreso.config(text=self.t('serie_leyendo', i=i + 1, n=len(rutas), nombre=self._nombre_ruta(ruta)))
                barra['value'] = i
                v.update()
            except tk.TclError:
                return   # se ha cerrado la ventana: se deja de revisar
            try:
                imagen, cabecera, _nota = leer_fits(*separar_ruta(ruta))
            except Exception as error:
                resultados.append({'ruta': ruta, 'error': str(error)})
                continue
            est = C.estadisticas_imagen(imagen, C.nivel_saturacion(cabecera))
            if primeras is None:
                primeras = est['fuentes']
                desplazamiento = (0.0, 0.0, len(primeras))
            else:
                desplazamiento = C.desplazamiento_entre(primeras, est['fuentes'])
            est['desplazamiento'] = desplazamiento
            est['ruta'] = ruta
            est['jd'] = instante_observacion(cabecera)[1]
            # se guardan las 400 más brillantes (x, y, flujo) para buscar objetos móviles
            est['fuentes_xyf'] = [(f['x'], f['y'], f['flujo']) for f in est['fuentes'][:400] if not f['saturada']]
            del est['fuentes']
            resultados.append(est)
            del imagen
        validos = [r for r in resultados if 'error' not in r]
        motivos = C.fotos_defectuosas(validos)
        for r, lista in zip(validos, motivos):
            r['motivos'] = lista
        self.revision = resultados
        barra.destroy()
        progreso.config(text=self.t('serie_resumen', n=len(validos), malas=sum(1 for r in validos if r['motivos'])))
        columnas = [('#', 40, tk.CENTER), (self.t('col_foto'), 200, tk.W), (self.t('serie_fondo'), 115, tk.E),
                    (self.t('serie_ruido'), 80, tk.E), (self.t('serie_estrellas'), 80, tk.E), ('FWHM (px)', 80, tk.E),
                    (self.t('serie_elong'), 100, tk.E), ('dx, dy (px)', 120, tk.E), (self.t('serie_diagnostico'), 260, tk.W)]
        marco_tabla, tabla = self._tabla_simple(marco, columnas, alto=min(12, len(resultados)))
        marco_tabla.pack(fill=tk.X)
        tabla.tag_configure('mala', foreground='#ffb36b')
        for i, r in enumerate(resultados, 1):
            if 'error' in r:
                tabla.insert('', tk.END, values=(i, self._nombre_ruta(r['ruta']), '', '', '', '', '', '', r['error']),
                             tags=('mala',))
                continue
            d = r['desplazamiento']
            tabla.insert('', tk.END, values=(
                i, self._nombre_ruta(r['ruta']), numero(r['fondo'], 1), numero(r['ruido'], 2), r['n_fuentes'],
                numero(r['fwhm'], 2), numero(r['elongacion'], 2),
                f'{numero(d[0], 1)}, {numero(d[1], 1)}' if d else self.t('serie_sin_registro'),
                ', '.join(self.t('motivo_' + x) for x in r['motivos']) or '✓'), tags=('mala',) if r['motivos'] else ())
        contenedor, figura, (ax1, ax2), lienzo = self._grafico(marco, 8.4, 3.8, filas=2)
        contenedor.pack(fill=tk.BOTH, expand=True, pady=(8, 0))
        x = list(range(1, len(validos) + 1))
        for ax in (ax1, ax2):
            self._estilo_ejes(ax, enteros=True)
        ax1.plot(x, [r['fwhm'] for r in validos], 'o-', color='#ffd27a', ms=4)
        ax1.set_ylabel('FWHM (px)')
        ax2.plot(x, [r['fondo'] for r in validos], 'o-', color='#86c8ff', ms=4)
        ax2.set_ylabel(self.t('serie_fondo'))
        ax2.set_xlabel(self.t('eje_foto'))
        for i, r in enumerate(validos, 1):
            if r['motivos']:
                for ax, clave in ((ax1, 'fwhm'), (ax2, 'fondo')):
                    ax.plot([i], [r[clave]], 'o', color='#ff5c5c', ms=8, mfc='none')
        figura.tight_layout()
        lienzo.draw_idle()
        self._registrar('hist_revision', n=len(validos), malas=sum(1 for r in validos if r['motivos']))
        self._caja_buscar_moviles(marco, v)
        self._centrar(v)

    def _caja_buscar_moviles(self, marco, v):
        """Botón y resultados de «Buscar objectes mòbils» en la ventana de revisión de la serie."""
        fila = ttk.Frame(marco)
        fila.pack(fill=tk.X, pady=(8, 0))
        resultado = ttk.Label(fila, text='', style='Gris.TLabel')
        columnas = [('#', 40, tk.CENTER), (self.t('fot_n'), 60, tk.CENTER), (self.t('buscar_mov_vel'), 140, tk.E),
                    (self.t('buscar_mov_pos'), 260, tk.W)]
        marco_tabla, tabla = self._tabla_simple(marco, columnas, alto=4, seleccion=tk.BROWSE)
        var_letra = tk.StringVar(value=next((l for l in LETRAS_ESTRELLA if l not in self._letras_medidas()), 'Z'))
        estado = {'pistas': []}
        acciones = ttk.Frame(marco)

        def buscar():
            pistas = self.buscar_moviles_serie(v)
            if pistas is None:
                return
            estado['pistas'] = pistas
            tabla.delete(*tabla.get_children())
            for i, p in enumerate(pistas):
                k0 = min(p['posiciones'])
                x, y, _f = p['posiciones'][k0]
                tabla.insert('', tk.END, iid=str(i), values=(i + 1, p['n'], numero(p['rapidez_px_dia'] / 24, 2) + ' px/h',
                                                              f'X {numero(x, 1)}, Y {numero(y, 1)} (foto {k0 + 1})'))
            resultado.config(text=self.t('buscar_mov_n', n=len(pistas)))
            if not marco_tabla.winfo_ismapped():
                marco_tabla.pack(fill=tk.X, pady=(4, 0))
                acciones.pack(fill=tk.X, pady=(4, 0))

        def medir():
            seleccion = tabla.selection()
            if not seleccion:
                messagebox.showinfo(self.t('serie_titulo'), self.t('moviles_elige'), parent=v)
                return
            n = self.medir_pista(estado['pistas'][int(seleccion[0])], var_letra.get())
            messagebox.showinfo(self.t('serie_titulo'), self.t('buscar_mov_medido', n=n, letra=var_letra.get()), parent=v)
        b = ttk.Button(fila, text=self.t('buscar_mov_boton'), command=buscar)
        b.pack(side=tk.LEFT)
        Ayuda(b, self.t('ayuda_buscar_mov'))
        resultado.pack(side=tk.LEFT, padx=8)
        ttk.Label(acciones, text=self.t('buscar_mov_letra')).pack(side=tk.LEFT)
        ttk.Combobox(acciones, textvariable=var_letra, values=LETRAS_ESTRELLA, state='readonly', width=3).pack(side=tk.LEFT, padx=4)
        ttk.Button(acciones, text=self.t('buscar_mov_medir'), command=medir).pack(side=tk.LEFT, padx=6)

    def tabla_revision(self):
        columnas = [(self.t('col_foto'), None, 30), (self.t('col_jd'), 6, 17), (self.t('serie_fondo'), 2, 12),
                    (self.t('serie_ruido'), 3, 10), (self.t('serie_estrellas'), 0, 10), ('FWHM (px)', 2, 10),
                    (self.t('serie_elong'), 2, 10), ('dx (px)', 2, 9), ('dy (px)', 2, 9), (self.t('serie_diagnostico'), None, 40)]
        filas = []
        for r in self.revision or []:
            if 'error' in r:
                filas.append([self._nombre_ruta(r['ruta'])] + [None] * 8 + [r['error']])
                continue
            d = r['desplazamiento'] or (None, None)
            filas.append([self._nombre_ruta(r['ruta']), r.get('jd'), r['fondo'], r['ruido'], r['n_fuentes'], r['fwhm'],
                          r['elongacion'], d[0], d[1], ', '.join(self.t('motivo_' + x) for x in r.get('motivos', [])) or 'OK'])
        return columnas, filas

    # ------------------------------------------------------------------ catálogos y astrometría
    def _campo_foto(self, ruta=None):
        """Centro (AR, Dec) y radio (grados) del campo de la foto, o None si no tiene WCS."""
        foto = self.fotos[ruta or self.ruta]
        wcs = foto.get('wcs')
        if wcs is None:
            return None
        ar, dec = wcs.a_cielo((foto['ancho'] + 1) / 2, (foto['alto'] + 1) / 2)
        radio = math.hypot(foto['ancho'], foto['alto']) / 2 * wcs.escala / 3600 * 1.05
        return ar, dec, min(max(radio, 0.01), 3.0)

    def _codigo_observatorio(self):
        codigo = (self.ajustes.get('mpc', {}).get('codigo') or '500').strip().upper()
        return codigo if re.match(r'^[0-9A-Z]{3}$', codigo) and codigo != 'XXX' else '500'

    def _preguntar_y_descargar(self, catalogo, ar, dec, radio, jd=None, padre=None, al_terminar=None):
        """Pide permiso y descarga (en segundo plano) un trozo de catálogo. Es la ÚNICA
        forma en que Emillium se conecta a internet: siempre por una acción del usuario."""
        if catalogo in self._descargas:
            messagebox.showinfo(TITULO_VENTANA, self.t('cat_ya_descargando'), parent=padre)
            return
        servidor = {'gaia': 'ESA Gaia Archive / VizieR (CDS)', 'apass': 'VizieR (CDS)', 'vsx': 'VizieR (CDS)',
                    'skybot': 'SkyBoT (IMCCE)'}[catalogo]
        if not messagebox.askyesno(self.t('cat_confirmar_t'),
                                   self.t('cat_confirmar', cat=catalogos.NOMBRES[catalogo], servidor=servidor,
                                          ar=numero(ar, 4), dec=numero(dec, 4), radio=numero(radio, 3)),
                                   parent=padre):
            return
        resultado = {}

        def tarea():
            try:
                resultado['filas'], resultado['ruta'] = catalogos.descargar(
                    catalogo, ar, dec, radio, jd=jd, codigo_observatorio=self._codigo_observatorio())
            except catalogos.ErrorCatalogo as error:
                resultado['error'] = (error.clave, error.datos)
            except Exception as error:   # error inesperado: se muestra tal cual
                resultado['error'] = ('sin_conexion', {'detalle': str(error)[:200]})
        hilo = threading.Thread(target=tarea, daemon=True)
        self._descargas[catalogo] = hilo
        hilo.start()

        def comprobar():
            if hilo.is_alive():
                self.root.after(300, comprobar)
                return
            del self._descargas[catalogo]
            self._registrar('hist_descarga', cat=catalogos.NOMBRES[catalogo], ok='error' not in resultado,
                            n=len(resultado.get('filas') or []))
            if 'error' in resultado:
                clave, datos = resultado['error']
                messagebox.showerror(self.t('cat_error_t'), self.t('cat_err_' + clave, **datos), parent=padre)
                return
            if al_terminar is not None:
                try:
                    al_terminar(resultado['filas'])
                except tk.TclError:
                    pass   # la ventana se ha cerrado mientras se descargaba
        self.root.after(300, comprobar)

    def _gaia_en_epoca(self, ruta):
        """Estrellas de Gaia del campo (desde la caché local) llevadas a la fecha de la foto."""
        campo = self._campo_foto(ruta)
        if campo is None:
            return None
        datos = catalogos.desde_cache('gaia', *campo)
        if datos is None:
            return None
        filas, meta = datos
        foto = self.fotos[ruta]
        epoca = C.epoca_juliana(foto['jd']) if foto.get('jd') else (meta.get('epoca') or 2016.0)
        ar, dec = C.aplicar_movimiento_propio([f['ar'] for f in filas], [f['dec'] for f in filas],
                                              [f.get('pmra') or 0 for f in filas], [f.get('pmdec') or 0 for f in filas],
                                              meta.get('epoca') or 2016.0, epoca)
        return filas, ar, dec

    def comprobar_astrometria(self):
        """Compara la solución astrométrica de la foto abierta con Gaia DR3 (de la caché):
        detecta las estrellas de la foto, las pasa a AR/Dec con la solución y busca su
        pareja en Gaia. Devuelve un diccionario con los residuos o un texto de error."""
        foto = self.fotos[self.ruta]
        wcs = foto.get('wcs')
        if wcs is None:
            return self.t('ast_sin_wcs')
        gaia = self._gaia_en_epoca(self.ruta)
        if gaia is None:
            return self.t('ast_sin_gaia')
        filas, ar_cat, dec_cat = gaia
        fuentes = [f for f in (self.fuentes or C.detectar_fuentes(self.original, 5.0, 1500,
                                                                   saturacion=foto.get('saturacion')))
                   if not f['saturada']][:600]
        if len(fuentes) < 3:
            return self.t('ast_pocas_estrellas')
        # Solo las estrellas de Gaia más brillantes (unas 3 por cada estrella detectada): con
        # todas, en un campo con muchas estrellas cualquier punto encontraría "pareja" por azar.
        orden = sorted(range(len(filas)), key=lambda i: filas[i].get('mag') if filas[i].get('mag') is not None else 99)
        orden = np.array(orden[:max(3 * len(fuentes), 50)], dtype=int)
        filas = [filas[i] for i in orden]
        ar_cat, dec_cat = np.asarray(ar_cat)[orden], np.asarray(dec_cat)[orden]
        campo = self._campo_foto()
        densidad = len(filas) / (math.pi * campo[2] ** 2)          # estrellas del catálogo por grado²
        cielo = [wcs.a_cielo(f['x'] + 1, f['y'] + 1) for f in fuentes]
        ar_med, dec_med = np.array([c[0] for c in cielo]), np.array([c[1] for c in cielo])
        # 1) emparejado amplio para encontrar un posible desplazamiento sistemático
        radio = max(3.0, 4 * wcs.escala)
        idx, _d = C.emparejar(ar_med, dec_med, ar_cat, dec_cat, radio)
        ok = idx >= 0
        if ok.sum() < 3:
            return self.t('ast_sin_parejas', n=int(ok.sum()))
        cosd = np.cos(np.radians(dec_med))
        d_ar = np.median((ar_cat[idx[ok]] - ar_med[ok] + 180) % 360 - 180) * 1.0
        d_dec = np.median(dec_cat[idx[ok]] - dec_med[ok])
        # 2) emparejado fino tras quitar el desplazamiento
        radio_fino = max(1.5, 2 * wcs.escala)
        idx, dist = C.emparejar(ar_med + d_ar, dec_med + d_dec, ar_cat, dec_cat, radio_fino)
        # cada estrella de Gaia solo puede emparejarse una vez (la más cercana)
        usadas = {}
        for i, (j, d) in enumerate(zip(idx, dist)):
            if j >= 0 and (j not in usadas or d < dist[usadas[j]]):
                usadas[j] = i
        pares = sorted((i, j) for j, i in usadas.items())
        if len(pares) < 3:
            return self.t('ast_sin_parejas', n=len(pares))
        i_med = np.array([p[0] for p in pares])
        i_cat = np.array([p[1] for p in pares])
        res_ar = ((ar_cat[i_cat] - ar_med[i_med] + 180) % 360 - 180) * cosd[i_med] * 3600
        res_dec = (dec_cat[i_cat] - dec_med[i_med]) * 3600
        return {'fuentes': [fuentes[i] for i in i_med], 'ar_cat': ar_cat[i_cat], 'dec_cat': dec_cat[i_cat],
                'filas_cat': [filas[j] for j in i_cat], 'res_ar': res_ar, 'res_dec': res_dec,
                'media_ar': float(res_ar.mean()), 'media_dec': float(res_dec.mean()),
                'rms': float(np.sqrt(np.mean(res_ar ** 2 + res_dec ** 2))),
                'rms_sin_media': float(np.sqrt(np.mean((res_ar - res_ar.mean()) ** 2 + (res_dec - res_dec.mean()) ** 2))),
                'n': len(pares), 'n_fuentes': len(fuentes), 'n_catalogo': len(filas),
                # parejas que saldrían por pura casualidad (estrellas detectadas × probabilidad
                # de que caiga una del catálogo dentro del círculo de búsqueda)
                'azar': len(fuentes) * min(1.0, densidad * math.pi * (radio_fino / 3600) ** 2)}

    def aplicar_wcs(self, ruta, cabecera_wcs, rms, ajustada=True):
        """Usa una solución astrométrica nueva en una foto y recalcula AR/Dec de sus medidas."""
        foto = self.fotos[ruta]
        foto['wcs'] = SolucionWCS(cabecera_wcs)
        foto['wcs_ajustada'] = ajustada
        foto['rms_astrometrico'] = rms
        self._registrar('hist_wcs', n=foto['n'], rms=numero(rms, 3) if rms is not None else '—')
        for m in self.mediciones:
            if m['ruta'] == ruta:
                self._actualizar_cielo(m)
                if self.tabla.exists(str(m['n'])):
                    self.tabla.item(str(m['n']), values=self._valores_fila(m))
        self._mostrar_info_foto()
        if self.mediciones:
            self._mostrar_resultado(self.mediciones[-1])

    def comparar_con_catalogos(self):
        """Para cada estrella medida en la foto abierta, su pareja en los catálogos guardados
        (Gaia, APASS, VSX y catálogos propios). Devuelve una lista de filas de texto."""
        foto = self.fotos[self.ruta]
        medidas = [m for m in self.mediciones if m['ruta'] == self.ruta and m.get('ar') is not None]
        gaia = self._gaia_en_epoca(self.ruta)
        campo = self._campo_foto()
        apass = catalogos.desde_cache('apass', *campo) if campo else None
        vsx = catalogos.desde_cache('vsx', *campo) if campo else None
        local = catalogos.desde_cache('local', *campo) if campo else None
        filas = []
        for m in sorted(medidas, key=lambda m: LETRAS_ESTRELLA.index(m['estrella'])):
            fila = [m['estrella'], '', '', '', '', '', '']
            if gaia is not None:
                idx, dist = C.emparejar([m['ar']], [m['dec']], gaia[1], gaia[2], max(2.0, 2 * foto['wcs'].escala))
                if idx[0] >= 0:
                    g = gaia[0][idx[0]]
                    d_ar = ((gaia[1][idx[0]] - m['ar'] + 180) % 360 - 180) * math.cos(math.radians(m['dec'])) * 3600
                    d_dec = (gaia[2][idx[0]] - m['dec']) * 3600
                    fila[1:5] = [g['id'].replace('Gaia DR3 ', ''), numero(g.get('mag'), 2), numero(d_ar, 2), numero(d_dec, 2)]
            if apass is not None:
                filas_a = apass[0]
                idx, _d = C.emparejar([m['ar']], [m['dec']], [f['ar'] for f in filas_a], [f['dec'] for f in filas_a], 3.0)
                if idx[0] >= 0:
                    fila[5] = numero(filas_a[idx[0]].get('V'), 2)
            extra = []
            for datos, radio in ((vsx, 10.0), (local, 3.0)):
                if datos is not None and datos[0]:
                    lista = datos[0]
                    idx, _d = C.emparejar([m['ar']], [m['dec']], [f['ar'] for f in lista], [f['dec'] for f in lista], radio)
                    if idx[0] >= 0:
                        f = lista[idx[0]]
                        extra.append(f"{f['id']}" + (f" ({f['tipo']})" if f.get('tipo') else ''))
            fila[6] = '; '.join(extra)
            filas.append(fila)
        return filas

    def ventana_catalogos(self):
        if self._necesita_foto():
            return
        v, marco = self._ventana(self.t('cat_titulo'), 'v_catalogos')
        foto = self.fotos[self.ruta]
        campo = self._campo_foto()
        texto = self.t('cat_explicacion') + '\n\n'
        texto += (self.t('cat_campo', ar=ar_hms(campo[0]), dec=dec_dms(campo[1]), radio=numero(campo[2] * 60, 1))
                  if campo else self.t('cat_sin_wcs'))
        ttk.Label(marco, text=texto, wraplength=900, style='Gris.TLabel', justify=tk.LEFT).pack(fill=tk.X)

        # --- catálogos guardados y botones de descarga
        caja = ttk.LabelFrame(marco, text=self.t('cat_locales'), padding=8)
        caja.pack(fill=tk.X, pady=(8, 4))
        etiquetas = {}
        for fila, catalogo in enumerate(('gaia', 'apass', 'vsx', 'ngc', 'skybot', 'local')):
            nombre = catalogos.NOMBRES.get(catalogo, self.t('cat_propio'))
            ttk.Label(caja, text=nombre, font=('Segoe UI', 10, 'bold')).grid(row=fila, column=0, sticky=tk.W, padx=(0, 10))
            etiquetas[catalogo] = ttk.Label(caja, text='', style='Gris.TLabel')
            etiquetas[catalogo].grid(row=fila, column=1, sticky=tk.W)
            if catalogo == 'local':
                ttk.Button(caja, text=self.t('cat_importar'), style='Pequeno.TButton',
                           command=lambda: importar()).grid(row=fila, column=2, sticky=tk.E, pady=1)
            elif campo:
                ttk.Button(caja, text=self.t('cat_descargar'), style='Pequeno.TButton',
                           command=lambda c=catalogo: descargar(c)).grid(row=fila, column=2, sticky=tk.E, pady=1)
        ttk.Button(caja, text=self.t('gestor_boton'), style='Pequeno.TButton',
                   command=lambda: self.ventana_gestor_catalogos(v)).grid(row=11, column=2, sticky=tk.E, pady=(6, 0))
        caja.columnconfigure(1, weight=1)

        def estado():
            for catalogo, etiqueta in etiquetas.items():
                if campo is None:
                    etiqueta.config(text='—')
                    continue
                datos = catalogos.desde_cache(catalogo, *campo, jd=foto.get('jd'))
                if datos is None:
                    etiqueta.config(text=self.t('cat_no_guardado'))
                else:
                    etiqueta.config(text=self.t('cat_guardado', n=len(datos[0]),
                                                fecha=datetime.fromtimestamp(datos[1]['fecha']).strftime('%Y-%m-%d %H:%M')))

        def descargar(catalogo):
            radio = campo[2] * (1.0 if catalogo not in ('vsx', 'ngc') else 1.2)

            def terminado(_filas):
                estado()
                if catalogo in ('skybot', 'vsx', 'ngc'):
                    capas[catalogo].set(True)
                    marcar()
            self._preguntar_y_descargar(catalogo, campo[0], campo[1], radio, jd=foto.get('jd'), padre=v,
                                        al_terminar=terminado)

        def importar():
            ruta = filedialog.askopenfilename(parent=v, title=self.t('cat_importar'),
                                              filetypes=[('CSV', '*.csv *.txt *.tsv'), (self.t('tipo_todos'), '*.*')])
            if not ruta:
                return
            try:
                n, _ruta = catalogos.importar_csv(ruta)
            except catalogos.ErrorCatalogo as error:
                messagebox.showerror(self.t('cat_error_t'), self.t('cat_err_' + error.clave, **error.datos), parent=v)
                return
            except OSError as error:
                messagebox.showerror(self.t('cat_error_t'), str(error), parent=v)
                return
            messagebox.showinfo(self.t('cat_titulo'), self.t('cat_importado', n=n), parent=v)
            estado()

        # --- marcas de los catálogos guardados sobre la foto
        capas = {c: tk.BooleanVar(value=c in foto.get('capas', set())) for c in ('vsx', 'ngc', 'skybot', 'gaia')}

        def marcar():
            elegidas = {c for c, var in capas.items() if var.get()}
            foto['capas'] = elegidas
            n = self.actualizar_marcas_catalogo(self.ruta, elegidas)
            lbl_capas.config(text=self.t('capas_n', n=n) if elegidas else '')
        if campo:
            fila_capas = ttk.Frame(caja)
            fila_capas.grid(row=10, column=0, columnspan=3, sticky=tk.W, pady=(6, 0))
            ttk.Label(fila_capas, text=self.t('capas_mostrar')).pack(side=tk.LEFT, padx=(0, 6))
            for c, clave in (('vsx', 'capa_vsx'), ('ngc', 'capa_ngc'), ('skybot', 'capa_skybot'), ('gaia', 'capa_gaia')):
                ttk.Checkbutton(fila_capas, text=self.t(clave), variable=capas[c], command=marcar).pack(side=tk.LEFT, padx=4)
            lbl_capas = ttk.Label(fila_capas, text='', style='Gris.TLabel')
            lbl_capas.pack(side=tk.LEFT, padx=8)

        # --- resolver la foto con Gaia (si no tiene solución astrométrica)
        if campo is None:
            self._caja_resolver(marco, v, foto)

        # --- astrometría
        caja_ast = ttk.LabelFrame(marco, text=self.t('ast_titulo'), padding=8)
        caja_ast.pack(fill=tk.X, pady=4)
        resultado_lbl = ttk.Label(caja_ast, text='', font=('Consolas', 10), justify=tk.LEFT, foreground=TEXTO)
        resultado_lbl.pack(fill=tk.X)
        botones_ast = ttk.Frame(caja_ast)
        botones_ast.pack(fill=tk.X, pady=(6, 0))
        contenedor, figura, ax, lienzo = self._grafico(marco, 8.0, 2.7)
        estado_ast = {'comprobacion': None, 'ajuste': None}

        def comprobar():
            self.root.config(cursor='watch')
            v.config(cursor='watch')
            v.update_idletasks()
            try:
                r = self.comprobar_astrometria()
            finally:
                self.root.config(cursor='')
                v.config(cursor='')
            if isinstance(r, str):
                resultado_lbl.config(text=r, foreground=ROJO_SUAVE)
                return
            estado_ast['comprobacion'] = r
            foto['rms_astrometrico'] = r['rms']
            resultado_lbl.config(foreground=TEXTO, text=self.t(
                'ast_resultado', n=r['n'], nf=r['n_fuentes'], nc=r['n_catalogo'], dar=numero(r['media_ar'], 3),
                ddec=numero(r['media_dec'], 3), rms=numero(r['rms'], 3), rms2=numero(r['rms_sin_media'], 3),
                azar=numero(r['azar'], 1)))
            dibujar_residuos(r)
            boton_ajustar.state(['!disabled'])
            boton_sip.state(['!disabled'])

        def dibujar_residuos(r, titulo=None):
            if not contenedor.winfo_ismapped():
                contenedor.pack(fill=tk.BOTH, expand=True, pady=(8, 0))
            ax.clear()
            self._estilo_ejes(ax)
            xs = [f['x'] for f in r['fuentes']]
            ys = [f['y'] for f in r['fuentes']]
            q = ax.quiver(xs, ys, r['res_ar'], r['res_dec'], color='#ffd27a', angles='xy', scale_units='xy',
                          scale=max(float(np.hypot(r['res_ar'], r['res_dec']).max()) / (0.08 * max(foto['ancho'], foto['alto'])), 1e-6))
            ax.quiverkey(q, 0.92, -0.16, 1.0, '1″', labelcolor=TEXTO, color='#ffd27a')
            ax.set_xlim(0, foto['ancho'])
            ax.set_ylim(0, foto['alto'])
            ax.set_aspect('equal')
            ax.set_xlabel('X (px)')
            ax.set_ylabel('Y (px)')
            ax.set_title(titulo or self.t('ast_grafico'), fontsize=9)
            figura.tight_layout()
            lienzo.draw_idle()

        def ajustar():
            r = estado_ast['comprobacion']
            if r is None:
                return
            wcs = foto['wcs']
            aj = C.ajustar_wcs([f['x'] + 1 for f in r['fuentes']], [f['y'] + 1 for f in r['fuentes']],
                               r['ar_cat'], r['dec_cat'], crpix=wcs.crpix)
            if aj is None:
                resultado_lbl.config(text=self.t('ast_sin_parejas', n=r['n']), foreground=ROJO_SUAVE)
                return
            estado_ast['ajuste'] = aj
            resultado_lbl.config(foreground=TEXTO, text=resultado_lbl.cget('text') + '\n\n' + self.t(
                'ast_ajuste', n=aj['n_usadas'], total=aj['n'], rms=numero(aj['rms'], 3), escala=numero(aj['escala'], 4),
                escala0=numero(wcs.escala, 4), rot=numero(aj['rotacion'], 2),
                sip=self.t('ast_aviso_sip') if wcs.sip else ''))
            usadas = aj['usadas']
            dibujar_residuos({'fuentes': [f for f, u in zip(r['fuentes'], usadas) if u],
                              'res_ar': -aj['residuo_x'][usadas], 'res_dec': -aj['residuo_y'][usadas]},
                             self.t('ast_grafico_ajuste'))
            boton_aplicar.state(['!disabled'])

        def aplicar():
            aj = estado_ast['ajuste']
            if aj is None:
                return
            self.aplicar_wcs(self.ruta, aj['cabecera'], aj['rms'])
            boton_aplicar.state(['disabled'])
            resultado_lbl.config(text=resultado_lbl.cget('text') + '\n' + self.t('ast_aplicada'))

        ttk.Button(botones_ast, text=self.t('ast_comprobar'), command=comprobar).pack(side=tk.LEFT)
        boton_ajustar = ttk.Button(botones_ast, text=self.t('ast_ajustar'), command=ajustar)
        boton_ajustar.pack(side=tk.LEFT, padx=6)
        def ajustar_sip():
            r = estado_ast['comprobacion']
            if r is None:
                return
            orden = 3 if r['n'] >= 30 else 2
            aj = C.ajustar_sip([f['x'] + 1 for f in r['fuentes']], [f['y'] + 1 for f in r['fuentes']],
                               r['ar_cat'], r['dec_cat'], orden=orden, crpix=foto['wcs'].crpix)
            if aj is None:
                resultado_lbl.config(text=resultado_lbl.cget('text') + '\n\n' + self.t('ast_sip_pocas', n=r['n']),
                                     foreground=TEXTO)
                return
            estado_ast['ajuste'] = aj
            resultado_lbl.config(foreground=TEXTO, text=resultado_lbl.cget('text') + '\n\n' + self.t(
                'ast_sip', orden=orden, n=aj['n_usadas'], total=aj['n'], rms=numero(aj['rms'], 3)))
            usadas = aj['usadas']
            dibujar_residuos({'fuentes': [f for f, u in zip(r['fuentes'], usadas) if u],
                              'res_ar': -aj['residuo_x'][usadas], 'res_dec': -aj['residuo_y'][usadas]},
                             self.t('ast_grafico_sip'))
            boton_aplicar.state(['!disabled'])
        boton_sip = ttk.Button(botones_ast, text=self.t('ast_ajustar_sip'), command=ajustar_sip)
        boton_sip.pack(side=tk.LEFT, padx=(0, 6))
        boton_sip.state(['disabled'])
        boton_aplicar = ttk.Button(botones_ast, text=self.t('ast_aplicar'), command=aplicar)
        boton_aplicar.pack(side=tk.LEFT)
        boton_ajustar.state(['disabled'])
        boton_aplicar.state(['disabled'])
        if campo is None:
            ttk.Button(botones_ast, text=self.t('boton_wcs_manual'), command=self.ventana_wcs_manual).pack(side=tk.LEFT, padx=6)

        # --- identificación cruzada de las estrellas medidas
        caja_cmp = ttk.LabelFrame(marco, text=self.t('cmp_titulo'), padding=8)
        caja_cmp.pack(fill=tk.X, pady=4)
        columnas = [(self.t('col_estrella'), 60, tk.CENTER), ('Gaia DR3', 170, tk.W), ('G', 60, tk.E),
                    ('ΔAR·cosδ (″)', 100, tk.E), ('ΔDec (″)', 90, tk.E), ('V (APASS)', 80, tk.E),
                    (self.t('cmp_otros'), 260, tk.W)]
        marco_tabla, tabla = self._tabla_simple(caja_cmp, columnas, alto=5)

        def comparar():
            if campo is None:
                return
            tabla.delete(*tabla.get_children())
            for fila in self.comparar_con_catalogos():
                tabla.insert('', tk.END, values=fila)
            if not marco_tabla.winfo_ismapped():
                marco_tabla.pack(fill=tk.X, pady=(6, 0))
        ttk.Button(caja_cmp, text=self.t('cmp_boton'), command=comparar).pack(anchor=tk.W)
        estado()
        self._centrar(v)

    def _caja_resolver(self, marco, v, foto):
        """En la ventana de catálogos, para una foto sin solución astrométrica: resolverla con
        Gaia a partir de la posición y la escala aproximadas (de la cabecera o escritas)."""
        caja = ttk.LabelFrame(marco, text=self.t('resolver_titulo'), padding=8)
        caja.pack(fill=tk.X, pady=4)
        ttk.Label(caja, text=self.t('resolver_explicacion'), wraplength=860, style='Gris.TLabel',
                  justify=tk.LEFT).grid(row=0, column=0, columnspan=6, sticky=tk.W)
        apuntado = foto.get('apuntado') or {}
        guardado = self.ajustes.get('resolver', {})
        var_ar = tk.StringVar(value=ar_hms(apuntado['ar']).replace('h', '').replace('m', '').replace('s', '')
                              if apuntado.get('ar') is not None else guardado.get('ar', ''))
        var_dec = tk.StringVar(value=dec_dms(apuntado['dec']).replace('°', '').replace('′', '').replace('″', '')
                               if apuntado.get('dec') is not None else guardado.get('dec', ''))
        var_escala = tk.StringVar(value=numero(apuntado['escala'], 3) if apuntado.get('escala') else guardado.get('escala', ''))
        for i, (clave, var) in enumerate((('wcsm_ar', var_ar), ('wcsm_dec', var_dec), ('resolver_escala', var_escala))):
            ttk.Label(caja, text=self.t(clave)).grid(row=1, column=2 * i, sticky=tk.W, padx=(0 if i == 0 else 12, 4), pady=(6, 0))
            ttk.Entry(caja, textvariable=var, width=16).grid(row=1, column=2 * i + 1, sticky=tk.W, pady=(6, 0))
        resultado = ttk.Label(caja, text='', font=('Consolas', 10), justify=tk.LEFT)
        resultado.grid(row=3, column=0, columnspan=6, sticky=tk.W, pady=(6, 0))
        botones = ttk.Frame(caja)
        botones.grid(row=2, column=0, columnspan=6, sticky=tk.W, pady=(6, 0))

        def leer():
            ar = self._leer_coordenada(var_ar.get(), True)
            dec = self._leer_coordenada(var_dec.get(), False)
            escala = leer_angulo(var_escala.get()) if var_escala.get().strip() else None
            if ar is None or dec is None or not (0 <= ar < 360 and -90 <= dec <= 90):
                resultado.config(text=self.t('resolver_falta_pos'), foreground=ROJO_SUAVE)
                return None
            self.ajustes['resolver'] = {'ar': var_ar.get(), 'dec': var_dec.get(), 'escala': var_escala.get()}
            # radio de búsqueda: 1,5 veces el campo (con la escala) o 1° si no se sabe la escala
            radio = (math.hypot(foto['ancho'], foto['alto']) / 2 * escala / 3600 * 1.5) if escala else 1.0
            return ar, dec, escala, min(max(radio, 0.05), 3.0)

        def descargar():
            datos = leer()
            if datos:
                self._preguntar_y_descargar('gaia', datos[0], datos[1], datos[3], padre=v,
                                            al_terminar=lambda filas: resultado.config(
                                                text=self.t('resolver_descargado', n=len(filas)), foreground=TEXTO))

        def resolver():
            datos = leer()
            if not datos:
                return
            resultado.config(text=self.t('resolver_buscando'), foreground=TEXTO_SUAVE)
            v.config(cursor='watch')
            v.update()
            try:
                r = self.resolver_con_gaia(*datos)
            finally:
                v.config(cursor='')
            if isinstance(r, str):
                resultado.config(text=r, foreground=ROJO_SUAVE)
                return
            self.aplicar_wcs(self.ruta, r['cabecera'], r['rms'])
            self._registrar('hist_resuelta', n=r['n_usadas'], rms=numero(r['rms'], 2))
            messagebox.showinfo(self.t('resolver_titulo'), self.t(
                'resolver_ok', n=r['n_usadas'], rms=numero(r['rms'], 2), escala=numero(r['escala'], 3),
                rot=numero(r['rotacion'], 1), espejo=self.t('si') if r['espejada'] else self.t('no')), parent=v)
            self.ventana_catalogos()          # se rehace la ventana, ahora con solución astrométrica
        ttk.Button(botones, text=self.t('resolver_descargar'), command=descargar).pack(side=tk.LEFT)
        b = ttk.Button(botones, text=self.t('resolver_boton'), style='Grande.TButton', command=resolver)
        b.pack(side=tk.LEFT, padx=6)

    # ------------------------------------------------------------------ WCS a partir de estrellas identificadas
    @staticmethod
    def _leer_coordenada(texto, horas):
        """'10 00 29.6' / '10:00:29.6' (horas o grados) o '150.1234' (grados) -> grados."""
        texto = texto.strip().replace(',', '.')
        if not texto:
            return None
        partes = re.split(r'[\s:hmsd°′″\'"]+', texto)
        partes = [p for p in partes if p]
        try:
            if len(partes) == 1:
                return float(partes[0])
            signo = -1 if partes[0].startswith('-') else 1
            valor = abs(float(partes[0])) + float(partes[1]) / 60 + (float(partes[2]) if len(partes) > 2 else 0) / 3600
            return signo * valor * (15 if horas else 1)
        except ValueError:
            return None

    def ventana_wcs_manual(self):
        """Solución astrométrica a partir de 3 o más estrellas medidas cuyas coordenadas
        se conocen (de un atlas, Stellarium, Aladin...). Útil cuando la foto no tiene WCS."""
        if self._necesita_foto():
            return
        medidas = sorted((m for m in self.mediciones if m['ruta'] == self.ruta),
                         key=lambda m: LETRAS_ESTRELLA.index(m['estrella']))
        if len(medidas) < 3:
            messagebox.showinfo(TITULO_VENTANA, self.t('wcsm_pocas'))
            return
        v, marco = self._ventana(self.t('wcsm_titulo'), 'v_wcs_manual')
        ttk.Label(marco, text=self.t('wcsm_explicacion'), wraplength=640, style='Gris.TLabel', justify=tk.LEFT).grid(
            row=0, column=0, columnspan=5, sticky=tk.W, pady=(0, 8))
        for c, titulo in enumerate((self.t('col_estrella'), 'X (px)', 'Y (px)', self.t('wcsm_ar'), self.t('wcsm_dec'))):
            ttk.Label(marco, text=titulo, font=('Segoe UI', 10, 'bold')).grid(row=1, column=c, padx=4, sticky=tk.W)
        guardadas = self.fotos[self.ruta].setdefault('coordenadas_manuales', {})
        entradas = []
        for i, m in enumerate(medidas):
            fila = 2 + i
            ttk.Label(marco, text=m['estrella']).grid(row=fila, column=0, padx=4)
            x0, y0 = self._original_xy(m)
            ttk.Label(marco, text=numero(x0, 2)).grid(row=fila, column=1, padx=4, sticky=tk.E)
            ttk.Label(marco, text=numero(y0, 2)).grid(row=fila, column=2, padx=4, sticky=tk.E)
            var_ar = tk.StringVar(value=guardadas.get(m['estrella'], ('', ''))[0])
            var_dec = tk.StringVar(value=guardadas.get(m['estrella'], ('', ''))[1])
            ttk.Entry(marco, textvariable=var_ar, width=16).grid(row=fila, column=3, padx=4, pady=1)
            ttk.Entry(marco, textvariable=var_dec, width=16).grid(row=fila, column=4, padx=4, pady=1)
            entradas.append((m, var_ar, var_dec))
        resultado = ttk.Label(marco, text='', font=('Consolas', 10), justify=tk.LEFT)
        resultado.grid(row=3 + len(medidas), column=0, columnspan=5, sticky=tk.W, pady=(8, 0))
        estado = {'ajuste': None}

        def calcular():
            xs, ys, ars, decs = [], [], [], []
            for m, var_ar, var_dec in entradas:
                guardadas[m['estrella']] = (var_ar.get(), var_dec.get())
                ar, dec = self._leer_coordenada(var_ar.get(), True), self._leer_coordenada(var_dec.get(), False)
                if ar is None or dec is None:
                    continue
                if not (0 <= ar < 360 and -90 <= dec <= 90):
                    resultado.config(text=self.t('wcsm_fuera', letra=m['estrella']), foreground=ROJO_SUAVE)
                    return
                x0, y0 = self._original_xy(m)
                xs.append(x0 + 1)
                ys.append(y0 + 1)
                ars.append(ar)
                decs.append(dec)
            if len(xs) < 3:
                resultado.config(text=self.t('wcsm_pocas'), foreground=ROJO_SUAVE)
                return
            aj = C.ajustar_wcs(xs, ys, ars, decs, recorte=1e9)   # pocas estrellas: no se descarta ninguna
            if aj is None:
                resultado.config(text=self.t('wcsm_pocas'), foreground=ROJO_SUAVE)
                return
            estado['ajuste'] = aj
            texto = self.t('wcsm_resultado', n=aj['n'], rms=numero(aj['rms'], 2), escala=numero(aj['escala'], 3),
                           rot=numero(aj['rotacion'], 1), espejo=self.t('si') if aj['espejada'] else self.t('no'))
            if aj['n'] == 3:
                texto += '\n' + self.t('wcsm_solo3')
            resultado.config(text=texto, foreground=TEXTO)
            boton_aplicar.state(['!disabled'])

        def aplicar():
            aj = estado['ajuste']
            if aj is None:
                return
            self.aplicar_wcs(self.ruta, aj['cabecera'], aj['rms'] if aj['n'] > 3 else None)
            resultado.config(text=resultado.cget('text') + '\n' + self.t('ast_aplicada'))
            boton_aplicar.state(['disabled'])

        botones = ttk.Frame(marco)
        botones.grid(row=4 + len(medidas), column=0, columnspan=5, sticky=tk.E, pady=(8, 0))
        ttk.Button(botones, text=self.t('wcsm_calcular'), command=calcular).pack(side=tk.LEFT)
        boton_aplicar = ttk.Button(botones, text=self.t('ast_aplicar'), command=aplicar, style='Grande.TButton')
        boton_aplicar.pack(side=tk.LEFT, padx=6)
        boton_aplicar.state(['disabled'])
        self._centrar(v)

    # ------------------------------------------------------------------ sesiones
    def datos_sesion(self):
        """Todo lo necesario para seguir otro día: fotos, medidas y opciones."""
        fotos = {}
        for ruta, foto in self.fotos.items():
            copia = {k: v for k, v in foto.items() if k not in ('wcs', 'marcas_catalogo')}
            copia['wcs'] = foto['wcs'].como_cabecera() if foto.get('wcs') else None
            fotos[ruta] = copia
        return {'programa': 'Emillium Expanded', 'version': VERSION, 'fecha': datetime.now().isoformat(timespec='seconds'),
                'fotos': fotos, 'mediciones': self.mediciones, 'serie': self.serie, 'ruta': self.ruta,
                'contador': self._contador, 'historial': self.historial,
                'revision': [{k: v for k, v in r.items() if k != 'fuentes_xyf'} for r in self.revision or []] or None,
                'opciones': {'comps': self.var_comps.get(), 'color_t': self.var_color_t.get(),
                             'colores': self.var_colores.get(), 'busqueda': self.var_busqueda.get(),
                             'apertura': self.var_apertura.get(), 'umbral': self.var_umbral.get(),
                             'r_fot': self.var_r_fot.get(), 'r_int': self.var_r_int.get(), 'r_ext': self.var_r_ext.get(),
                             'rotacion': self.rotacion, 'mirall_h': self.mirall_h, 'mirall_v': self.mirall_v}}

    def guardar_sesion(self):
        if not self.mediciones and not self.fotos:
            messagebox.showinfo(TITULO_VENTANA, self.t('sin_medidas'))
            return
        ruta = filedialog.asksaveasfilename(title=self.t('boton_guardar_sesion'), defaultextension='.emillium',
                                            initialfile=self.t('sesion_archivo') + '.emillium',
                                            filetypes=[(self.t('sesion_tipo'), '*.emillium'), ('JSON', '*.json')])
        if not ruta:
            return
        try:
            temporal = ruta + '.tmp'
            with open(temporal, 'w', encoding='utf-8') as archivo:
                json.dump(self.datos_sesion(), archivo, ensure_ascii=False, indent=1,
                          default=lambda o: o.item() if hasattr(o, 'item') else list(o) if isinstance(o, (tuple, set)) else str(o))
            os.replace(temporal, ruta)
        except OSError as error:
            messagebox.showerror(self.t('error_guardar'), str(error))
            return
        self._guardado = True
        self._registrar('hist_sesion_guardada', ruta=ruta)
        messagebox.showinfo(self.t('guardado_t'), self.t('sesion_guardada', ruta=ruta))

    def abrir_sesion(self, ruta=None):
        if ruta is None:
            if self.mediciones and not messagebox.askyesno(self.t('boton_abrir_sesion'), self.t('sesion_reemplazar')):
                return
            ruta = filedialog.askopenfilename(title=self.t('boton_abrir_sesion'),
                                              filetypes=[(self.t('sesion_tipo'), '*.emillium *.json'), (self.t('tipo_todos'), '*.*')])
            if not ruta:
                return
        try:
            with open(ruta, encoding='utf-8') as archivo:
                datos = json.load(archivo)
            if datos.get('programa') != 'Emillium Expanded' or not isinstance(datos.get('fotos'), dict):
                raise ValueError(self.t('sesion_no_valida'))
            fotos = {}
            for r, foto in datos['fotos'].items():
                foto = dict(foto)
                foto['wcs'] = SolucionWCS.desde_cabecera(foto['wcs']) if foto.get('wcs') else None
                fotos[r] = foto
            mediciones = []
            for m in datos['mediciones']:
                m = dict(m)
                m['orient'] = tuple(m['orient']) if isinstance(m.get('orient'), list) else m.get('orient')
                if m.get('ruta') not in fotos:
                    raise ValueError(self.t('sesion_no_valida'))
                mediciones.append(m)
        except (OSError, ValueError, KeyError, TypeError) as error:
            detalle = str(error)
            if detalle != self.t('sesion_no_valida'):
                detalle = self.t('sesion_no_valida') + f'\n\n({detalle})'
            messagebox.showerror(self.t('error_abrir'), detalle)
            return
        # se cierra lo que hubiera y se pone lo de la sesión
        self.cerrar_foto()
        self.fotos, self.mediciones = fotos, mediciones
        self._contador = int(datos.get('contador') or max([m['n'] for m in mediciones] + [0]))
        self.revision = datos.get('revision')
        self.historial = list(datos.get('historial') or [])
        opciones = datos.get('opciones', {})
        for clave, var in (('comps', self.var_comps), ('color_t', self.var_color_t), ('colores', self.var_colores),
                           ('busqueda', self.var_busqueda), ('apertura', self.var_apertura), ('umbral', self.var_umbral),
                           ('r_fot', self.var_r_fot), ('r_int', self.var_r_int), ('r_ext', self.var_r_ext)):
            if clave in opciones:
                var.set(opciones[clave])
        if self._recalculo_pendiente:   # cambiar opciones no debe recalcular nada al abrir una sesión
            self.root.after_cancel(self._recalculo_pendiente)
            self._recalculo_pendiente = None
        self.tabla.delete(*self.tabla.get_children())
        for m in self.mediciones:
            self._anadir_a_tabla(m)
        self._actualizar_contador()
        self._actualizar_fotometria()
        serie = [r for r in datos.get('serie') or [] if os.path.exists(separar_ruta(r)[0])]
        actual = datos.get('ruta')
        if actual and os.path.exists(separar_ruta(actual)[0]):
            self.serie = serie if actual in serie else [actual]
            self.rotacion = float(opciones.get('rotacion', 0.0))
            self.mirall_h = bool(opciones.get('mirall_h'))
            self.mirall_v = bool(opciones.get('mirall_v'))
            self.cargar_fits(actual, mantener_orientacion=True)
            self.ajustar_vista()
        faltan = [self.fotos[r]['archivo'] for r in self.fotos if not os.path.exists(separar_ruta(r)[0])]
        texto = self.t('sesion_abierta', n=len(self.mediciones), f=len(self.fotos))
        if faltan:
            texto += '\n' + self.t('sesion_faltan', fotos=', '.join(faltan[:6]) + ('…' if len(faltan) > 6 else ''))
        self._mensaje(texto)
        self._registrar('hist_sesion_abierta', ruta=ruta)
        if self.mediciones:
            self._mostrar_resultado(self.mediciones[-1])
        self._actualizar_botones_serie()
        self._programar_dibujo()

    # =================================================================== HERRAMIENTAS DE LA 1.5 (II)
    # ------------------------------------------------------------------ historial de operaciones
    def _registrar(self, clave, **datos):
        """Apunta una operación en el historial (se guarda en la sesión y en el Excel)."""
        if not hasattr(self, 'historial'):
            self.historial = []
        self.historial.append({'hora': datetime.now().isoformat(sep=' ', timespec='seconds'), 'clave': clave,
                               'datos': {k: (v if isinstance(v, (str, int, float)) or v is None else str(v))
                                         for k, v in datos.items()}})
        if len(self.historial) > 5000:
            del self.historial[:1000]

    def _texto_historial(self, entrada):
        try:
            return self.t(entrada['clave'], **entrada['datos'])
        except (KeyError, IndexError, ValueError):
            return entrada['clave']

    def ventana_historial(self):
        v, marco = self._ventana(self.t('hist_titulo'), 'v_historial')
        columnas = [(self.t('hist_hora'), 150, tk.W), (self.t('hist_operacion'), 620, tk.W)]
        marco_tabla, tabla = self._tabla_simple(marco, columnas, alto=18)
        marco_tabla.pack(fill=tk.BOTH, expand=True)
        for entrada in reversed(getattr(self, 'historial', [])):
            tabla.insert('', tk.END, values=(entrada['hora'], self._texto_historial(entrada)))
        self._centrar(v)

    def tabla_historial(self):
        return ([(self.t('hist_hora'), None, 20), (self.t('hist_operacion'), None, 90)],
                [[e['hora'], self._texto_historial(e)] for e in getattr(self, 'historial', [])])

    def _nombre_ruta(self, ruta):
        if ruta in self.fotos:
            return self.fotos[ruta]['archivo']
        archivo, indice = separar_ruta(ruta)
        return os.path.basename(archivo) + (f' [{indice + 1}]' if indice else '')

    # ------------------------------------------------------------------ varias imágenes en un archivo
    def cambiar_imagen(self, paso):
        """Otra imagen del mismo archivo FITS (otra extensión u otro plano)."""
        if self.ruta is None:
            return
        foto = self.fotos[self.ruta]
        archivo, indice = separar_ruta(self.ruta)
        nuevo = indice + paso
        if 0 <= nuevo < foto.get('n_imagenes', 1):
            self.cargar_fits(unir_ruta(archivo, nuevo))

    def _actualizar_selector_imagen(self):
        if not hasattr(self, 'fila_imagen'):
            return
        foto = self.fotos.get(self.ruta) if self.ruta else None
        n = foto.get('n_imagenes', 1) if foto else 1
        if n > 1:
            i = foto.get('indice', 0)
            self.lbl_imagen.config(text=self.t('imagen_n', i=i + 1, n=n))
            self.boton_img_ant.state(['!disabled'] if i > 0 else ['disabled'])
            self.boton_img_sig.state(['!disabled'] if i < n - 1 else ['disabled'])
            if not self.fila_imagen.winfo_ismapped():
                self.fila_imagen.pack(fill=tk.X, pady=(6, 0), before=self.boton_cerrar_foto)
        elif self.fila_imagen.winfo_ismapped():
            self.fila_imagen.pack_forget()

    # ------------------------------------------------------------------ medir toda la serie
    def medir_serie(self):
        """Mide de un clic todas las fotos siguientes de la serie con «Medir como en la
        foto anterior» (los objetos móviles se buscan donde tienen que estar)."""
        if not self.serie or self.ruta not in self.serie or not any(m['ruta'] == self.ruta for m in self.mediciones):
            messagebox.showinfo(TITULO_VENTANA, self.t('serie_medir_falta'))
            return
        inicio = self.serie.index(self.ruta)
        resumen = []
        for i in range(inicio + 1, len(self.serie)):
            self._mensaje(self.t('serie_midiendo', i=i + 1, n=len(self.serie)))
            self.root.update()
            self.cargar_fits(self.serie[i], mantener_orientacion=True)
            if self.ruta != self.serie[i]:
                break                       # no se ha podido abrir
            self.medir_como_anterior()
            resumen.append(sum(1 for m in self.mediciones if m['ruta'] == self.ruta))
        self._registrar('hist_serie', n=len(resumen))
        self._mensaje(self.t('serie_medida', n=len(resumen), medidas=sum(resumen)))

    # ------------------------------------------------------------------ parpadeo
    def parpadeo(self):
        """Alterna la foto abierta con la anterior de la serie, alineadas: lo que se mueve
        «salta» y se ve enseguida (así se descubrieron muchos asteroides)."""
        if getattr(self, '_parpadeo', None):
            self._parar_parpadeo()
            return
        if self.original is None or len(self.serie) < 2 or self.ruta not in self.serie:
            messagebox.showinfo(TITULO_VENTANA, self.t('parpadeo_falta'))
            return
        i = self.serie.index(self.ruta)
        otra = self.serie[i - 1] if i > 0 else self.serie[i + 1]
        try:
            imagen_otra, _cab, _n = leer_fits(*separar_ruta(otra))
        except Exception as error:
            messagebox.showerror(self.t('error_abrir'), str(error))
            return
        # Se lleva la otra foto encima de esta (desplazamiento, giro y espejo, con precisión
        # de centésimas de píxel). Si no se consigue, se avisa y se usa un desplazamiento simple.
        self.root.config(cursor='watch')
        self.root.update_idletasks()
        try:
            alineacion = self._alineacion(self.ruta, otra, self.original, imagen_otra)
        finally:
            self.root.config(cursor='')
        if alineacion is not None:
            alineada = C.remuestrear(imagen_otra, alineacion['A'], alineacion['t'], self.original.shape)
        elif imagen_otra.shape == self.original.shape:
            d = C.desplazamiento_entre(C.detectar_fuentes(self.original, 5.0, 300),
                                       C.detectar_fuentes(imagen_otra, 5.0, 300)) or (0.0, 0.0, 0)
            alineada = np.roll(imagen_otra, (-int(round(d[1])), -int(round(d[0]))), axis=(0, 1))
        else:
            messagebox.showinfo(TITULO_VENTANA, self.t('parpadeo_tamano'))
            return
        alineada = transformar_imagen(alineada, self.angulo, self.espejo)
        self._parpadeo = {'imagenes': [self.imagen, alineada], 'i': 0, 'nombres': [self.fotos[self.ruta]['archivo'],
                                                                                  self._nombre_ruta(otra)]}
        self._registrar('hist_parpadeo', a=self.fotos[self.ruta]['archivo'], b=self._nombre_ruta(otra))
        self._paso_parpadeo()

    def _paso_parpadeo(self):
        estado = getattr(self, '_parpadeo', None)
        if not estado:
            return
        estado['i'] = 1 - estado['i']
        self.imagen = estado['imagenes'][estado['i']]
        self._mensaje(self.t('parpadeo_mostrando', nombre=estado['nombres'][estado['i']]))
        self._programar_dibujo()
        estado['after'] = self.root.after(650, self._paso_parpadeo)

    def _parar_parpadeo(self):
        estado = getattr(self, '_parpadeo', None)
        if not estado:
            return
        if estado.get('after'):
            self.root.after_cancel(estado['after'])
        self.imagen = estado['imagenes'][0]
        self._parpadeo = None
        self._mensaje('')
        self._programar_dibujo()

    # ------------------------------------------------------------------ alinear y apilar
    def _alineacion(self, referencia, otra, imagen_ref=None, imagen_otra=None):
        """Cómo se lleva la foto 'otra' encima de la 'referencia' (C.alinear_fuentes), o None.
        Se recuerda para no repetir el cálculo."""
        clave = (referencia, otra)
        if clave in self.alineaciones:
            return self.alineaciones[clave]
        if referencia == otra:
            resultado = {'A': np.eye(2), 't': np.zeros(2), 'n': 0, 'rms': 0.0, 'giro': 0.0, 'escala': 1.0, 'espejo': False}
        else:
            fuentes = []
            for ruta, imagen in ((referencia, imagen_ref), (otra, imagen_otra)):
                cabecera = {}
                if imagen is None:
                    imagen, cabecera, _n = leer_fits(*separar_ruta(ruta))
                saturacion = self.fotos.get(ruta, {}).get('saturacion') or C.nivel_saturacion(cabecera)
                fuentes.append(C.fuentes_para_alinear(imagen, saturacion=saturacion))
            resultado = C.alinear_fuentes(fuentes[0], fuentes[1])
        self.alineaciones[clave] = resultado
        return resultado

    def ventana_alinear(self):
        """Alinea todas las fotos de la serie con la abierta (la referencia) y las apila en una
        foto nueva más profunda (menos ruido: se ven objetos más débiles)."""
        rutas = [r for r in self.serie if os.path.exists(separar_ruta(r)[0])]
        if self.original is None or len(rutas) < 2 or self.ruta not in rutas:
            messagebox.showinfo(TITULO_VENTANA, self.t('ali_falta'))
            return
        referencia = self.ruta
        v, marco = self._ventana(self.t('ali_titulo'), 'v_alinear')
        ttk.Label(marco, text=self.t('ali_explicacion', ref=self._nombre_ruta(referencia)), wraplength=820,
                  style='Gris.TLabel', justify=tk.LEFT).pack(fill=tk.X)
        progreso = ttk.Label(marco, text='', style='Gris.TLabel')
        progreso.pack(fill=tk.X, pady=(6, 0))
        barra = ttk.Progressbar(marco, maximum=len(rutas), style='Emillium.Horizontal.TProgressbar')
        barra.pack(fill=tk.X, pady=4)
        columnas = [('#', 40, tk.CENTER), (self.t('col_foto'), 230, tk.W), (self.t('ali_estrellas'), 90, tk.E),
                    ('dx, dy (px)', 130, tk.E), (self.t('ali_giro'), 80, tk.E), ('rms (px)', 80, tk.E),
                    (self.t('serie_diagnostico'), 170, tk.W)]
        marco_tabla, tabla = self._tabla_simple(marco, columnas, alto=min(12, len(rutas)), seleccion=tk.EXTENDED)
        marco_tabla.pack(fill=tk.X)
        tabla.tag_configure('mala', foreground='#ffb36b')
        alineadas = {}
        for i, ruta in enumerate(rutas):
            try:
                progreso.config(text=self.t('ali_alineando', i=i + 1, n=len(rutas), nombre=self._nombre_ruta(ruta)))
                barra['value'] = i
                v.update()
            except tk.TclError:
                return      # se ha cerrado la ventana
            try:
                r = self._alineacion(referencia, ruta, self.original if ruta == referencia else None)
            except Exception as error:
                tabla.insert('', tk.END, iid=ruta, values=(i + 1, self._nombre_ruta(ruta), '', '', '', '', str(error)[:60]),
                             tags=('mala',))
                continue
            if r is None:
                tabla.insert('', tk.END, iid=ruta, values=(i + 1, self._nombre_ruta(ruta), '', '', '', '', self.t('ali_no')),
                             tags=('mala',))
                continue
            alineadas[ruta] = r
            alto, ancho = self.original.shape
            centro = np.array([(ancho - 1) / 2, (alto - 1) / 2])
            dx, dy = r['A'] @ centro + r['t'] - centro
            giro = (r['giro'] + 180) % 360 - 180
            giro = 0.0 if abs(giro) < 0.005 else giro          # sin «-0,00»
            estado = self.t('ali_referencia') if ruta == referencia else ('✓' + (' ' + self.t('ali_espejo') if r['espejo'] else ''))
            tabla.insert('', tk.END, iid=ruta, values=(
                i + 1, self._nombre_ruta(ruta), r['n'] if ruta != referencia else '',
                f'{numero(dx, 2)}, {numero(dy, 2)}', numero(giro, 2), numero(r['rms'], 3) if ruta != referencia else '', estado))
        barra.destroy()
        progreso.config(text=self.t('ali_resumen', n=len(alineadas), total=len(rutas)))
        tabla.selection_set([r for r in rutas if r in alineadas])
        self._registrar('hist_alinear', n=len(alineadas), total=len(rutas))

        fila = ttk.Frame(marco)
        fila.pack(fill=tk.X, pady=(8, 0))
        ttk.Label(fila, text=self.t('ali_metodo')).pack(side=tk.LEFT)
        metodos = [self.t('ali_media'), self.t('ali_mediana')]
        var_metodo = tk.StringVar(value=metodos[0])
        ttk.Combobox(fila, textvariable=var_metodo, values=metodos, state='readonly', width=33).pack(side=tk.LEFT, padx=6)
        resultado = ttk.Label(marco, text=self.t('ali_seleccion'), style='Gris.TLabel', wraplength=820, justify=tk.LEFT)

        def apilar():
            elegidas = [r for r in rutas if r in tabla.selection() and r in alineadas]
            if referencia not in elegidas:
                elegidas.insert(0, referencia)
            if len(elegidas) < 2:
                messagebox.showinfo(self.t('ali_titulo'), self.t('ali_pocas'), parent=v)
                return
            base = os.path.splitext(separar_ruta(referencia)[0])[0]
            ruta_salida = filedialog.asksaveasfilename(
                parent=v, defaultextension='.fits', initialdir=os.path.dirname(base),
                initialfile=os.path.basename(base) + self.t('ali_sufijo', n=len(elegidas)) + '.fits',
                filetypes=[(self.t('tipo_fits'), '*.fits *.fit')])
            if not ruta_salida:
                return
            imagenes, transformaciones, mitades, exposiciones, ganancias = [], [], [], [], []
            cab_ref = None
            try:
                for ruta in elegidas:
                    progreso.config(text=self.t('ali_leyendo', nombre=self._nombre_ruta(ruta)))
                    v.update()
                    imagen, cabecera, _n = leer_fits(*separar_ruta(ruta))
                    if ruta == referencia:
                        cab_ref = cabecera
                    imagenes.append(imagen)
                    transformaciones.append((alineadas[ruta]['A'], alineadas[ruta]['t']))
                    mitades.append(instante_observacion(cabecera)[1])
                    exposiciones.append(exposicion_de(cabecera))
                    ganancias.append(ganancia_de(cabecera))

                def avance(fraccion):
                    progreso.config(text=self.t('ali_apilando', p=int(fraccion * 100)))
                    v.update()
                apilada, _cuantas = C.apilar(imagenes, transformaciones, self.original.shape,
                                             'mediana' if var_metodo.get() == metodos[1] else 'media', progreso=avance)
                del imagenes
                cabecera = self._cabecera_apilada(cab_ref, referencia, len(elegidas), mitades, exposiciones, ganancias,
                                                  var_metodo.get())
                escribir_fits(ruta_salida, apilada, cabecera,
                              [f'Emillium {VERSION}: {var_metodo.get()} de {len(elegidas)} fotos alineadas con '
                               f'{self._nombre_ruta(referencia)}'] + [self._nombre_ruta(r) for r in elegidas])
            except tk.TclError:
                return
            except (OSError, ErrorFits, MemoryError) as error:
                messagebox.showerror(self.t('error_guardar'), str(error), parent=v)
                return
            texto = self.t('ali_hecho', n=len(elegidas), ruta=os.path.basename(ruta_salida))
            progreso.config(text=texto)
            self._registrar('hist_apilar', n=len(elegidas), nombre=os.path.basename(ruta_salida))
            if messagebox.askyesno(self.t('ali_titulo'), texto + '\n\n' + self.t('ali_abrir'), parent=v):
                if ruta_salida not in self.serie:
                    self.serie.append(ruta_salida)
                self.fotos.pop(ruta_salida, None)          # por si se ha sobrescrito una apilada anterior
                self.cargar_fits(ruta_salida, mantener_orientacion=True)
        b = ttk.Button(fila, text=self.t('ali_apilar'), style='Grande.TButton', command=apilar)
        b.pack(side=tk.LEFT, padx=6)
        Ayuda(b, self.t('ayuda_apilar'))
        resultado.pack(fill=tk.X, pady=(6, 0))
        self._centrar(v)

    def _cabecera_apilada(self, cab_ref, referencia, n, mitades, exposiciones, ganancias, metodo):
        """Cabecera de la foto apilada: la de la referencia (su WCS sirve, porque la apilada está
        en su misma rejilla) con el instante, la exposición y la ganancia de la pila."""
        cabecera = {k: v for k, v in (cab_ref or {}).items()
                    if k not in ('DATE-OBS', 'TIME-OBS', 'UT', 'UTSTART', 'UTC-OBS', 'MJD-OBS', 'MJD-AVG', 'DATE-AVG',
                                 'DATE-END', 'JD', 'JD-OBS', 'EXPSTART', 'EXPEND', 'EXPTIME', 'EXPOSURE', 'GAIN',
                                 'EGAIN', 'CCDGAIN', 'NCOMBINE', 'BZERO', 'BSCALE', 'BLANK', 'DATAMIN', 'DATAMAX')}
        foto = self.fotos.get(referencia, {})
        if foto.get('wcs') is not None:          # la solución que usa Emillium (puede estar ajustada)
            cabecera = {k: v for k, v in cabecera.items() if not _CLAVES_WCS.match(k)}
            cabecera.update(foto['wcs'].como_cabecera())
        # Es una MEDIA: el brillo de cada estrella es el de una foto normal (no la suma), así
        # que la exposición equivalente es la media; la ganancia se multiplica por n (hay n
        # veces más electrones detrás de cada valor) y el ruido baja √n veces.
        if all(e for e in exposiciones):
            cabecera['EXPTIME'] = float(np.mean(exposiciones))
            cabecera['TOTEXP'] = float(np.sum(exposiciones))
        if all(g for g in ganancias):
            cabecera['GAIN'] = float(np.mean(ganancias)) * n
        if all(m for m in mitades):
            mitad = float(np.mean(mitades))
            inicio = mitad - (cabecera.get('EXPTIME') or 0.0) / 2 / 86400
            cabecera['DATE-OBS'] = fecha_desde_jd(inicio).replace(' ', 'T')
        cabecera['NCOMBINE'] = n
        cabecera['STACKMET'] = 'median' if metodo == self.t('ali_mediana') else 'sigma-clipped mean'
        return cabecera

    # ------------------------------------------------------------------ buscar objetos móviles en la serie
    def buscar_moviles_serie(self, padre=None):
        """Busca objetos desconocidos que se muevan en línea recta en la serie revisada."""
        revision = [r for r in (self.revision or []) if 'error' not in r and r.get('fuentes_xyf') is not None]
        if len(revision) < 3:
            messagebox.showinfo(TITULO_VENTANA, self.t('buscar_mov_falta'), parent=padre)
            return None
        listas = [[{'x': x, 'y': y, 'flujo': f} for x, y, f in r['fuentes_xyf']] for r in revision]
        desplazamientos = [r['desplazamiento'][:2] if r['desplazamiento'] else None for r in revision]
        tiempos = [r.get('jd') if r.get('jd') else i * 0.01 for i, r in enumerate(revision)]
        pistas = C.buscar_moviles(listas, desplazamientos, tiempos)
        for p in pistas:
            p['rutas'] = {revision[k]['ruta']: v for k, v in p['posiciones'].items()}
        self._registrar('hist_buscar_mov', n=len(pistas))
        return pistas

    def medir_pista(self, pista, letra):
        """Mide un objeto móvil encontrado en todas las fotos donde aparece, con la letra dada."""
        actual = self.ruta
        medidas = 0
        for ruta, (x, y, _f) in sorted(pista['rutas'].items(), key=lambda par: self.serie.index(par[0])
                                       if par[0] in self.serie else 0):
            if ruta != self.ruta:
                self.cargar_fits(ruta, mantener_orientacion=True)
            if self.ruta != ruta:
                continue
            foto = self.fotos[ruta]
            xv, yv = transformar_punto(x, y, foto['ancho'], foto['alto'], self.angulo, self.espejo)
            if self.medir_en(xv, yv, letra=letra, silencioso=True, radio_busqueda=4):
                medidas += 1
        if actual and actual != self.ruta:
            self.cargar_fits(actual, mantener_orientacion=True)
        self._registrar('hist_medir_pista', letra=letra, n=medidas)
        self._programar_dibujo()
        return medidas

    # ------------------------------------------------------------------ fotometría con todo el campo
    def _fotometria_campo(self, imagen, ruta, catalogo, banda):
        """Fotometría de apertura de las estrellas detectadas en una foto y su pareja en el
        catálogo guardado (APASS o Gaia). Devuelve (m_inst, σ, m_cat, color) o un texto de error."""
        foto = self.fotos.get(ruta, {})
        wcs = foto.get('wcs')
        if wcs is None:
            return self.t('ast_sin_wcs')
        campo = self._campo_foto(ruta)
        datos = catalogos.desde_cache(catalogo, *campo)
        if datos is None:
            return self.t('fot_cat_sin_cache', cat=catalogos.NOMBRES[catalogo])
        clave, color = {'V': ('V', ('B', 'V')), 'B': ('B', ('B', 'V')), 'g': ('g', ('g', 'r')), 'r': ('r', ('g', 'r')),
                        'i': ('i', ('r', 'i')), 'G': ('mag', None)}[banda]
        filas = [f for f in datos[0] if f.get(clave) is not None]
        if len(filas) < 3:
            return self.t('cal_pocas', n=len(filas))
        r, r_int, r_ext = self._radios_fotometria()
        fuentes = [f for f in C.detectar_fuentes(imagen, 8.0, 600, saturacion=foto.get('saturacion')) if not f['saturada']]
        cielo = [wcs.a_cielo(f['x'] + 1, f['y'] + 1) for f in fuentes]
        idx, _d = C.emparejar([c[0] for c in cielo], [c[1] for c in cielo], [f['ar'] for f in filas],
                              [f['dec'] for f in filas], max(2.0, 1.5 * wcs.escala))
        m_inst, sig, m_cat, colores = [], [], [], []
        for f, j in zip(fuentes, idx):
            if j < 0:
                continue
            fot = C.fotometria_apertura(imagen, f['x'], f['y'], r, r_int, r_ext, foto.get('ganancia'), foto.get('saturacion'))
            if fot is None or fot['saturada'] or fot['borde'] or not fot['snr'] > 20:
                continue
            # las estrellas con una vecina dentro de la apertura parecerían más brillantes
            ajena, _vecinas = C.contaminacion(imagen, f['x'], f['y'], r, f['fwhm'], fot['flujo'], umbral_sigma=3.0,
                                              fondo=fot['fondo'], sigma=fot['sigma_fondo'])
            if ajena > 0.01:
                continue
            m, s = C.magnitud(fot['flujo'], fot['sigma_flujo'], foto.get('exposicion'))
            fila = filas[j]
            if color is None:
                indice = fila.get('color')
            else:
                indice = (fila[color[0]] - fila[color[1]]) if fila.get(color[0]) is not None and fila.get(color[1]) is not None else None
            m_inst.append(m)
            sig.append(s)
            m_cat.append(fila[clave])
            colores.append(indice if indice is not None else np.nan)
        if len(m_inst) < 3:
            return self.t('cal_pocas', n=len(m_inst))
        return np.array(m_inst), np.array(sig), np.array(m_cat), np.array(colores, dtype=np.float64)

    def ventana_calibracion(self, banda_texto):
        """Calibración fotométrica con todo el campo: punto cero y término de color."""
        if self._necesita_foto():
            return
        banda = banda_texto[0].replace("'", '')
        catalogo = 'gaia' if banda == 'G' else 'apass'
        datos = self._fotometria_campo(self.original, self.ruta, catalogo, banda)
        if isinstance(datos, str):
            messagebox.showinfo(self.t('cal_titulo'), datos)
            return
        m_inst, sig, m_cat, color = datos
        con_color = np.isfinite(color).sum() >= 5
        cal = C.calibrar_campo(m_inst, m_cat, color if con_color else None, sig)
        if cal is None:
            messagebox.showinfo(self.t('cal_titulo'), self.t('cal_pocas', n=len(m_inst)))
            return
        v, marco = self._ventana(self.t('cal_titulo'), 'v_calibracion')
        texto = self.t('cal_resultado', n=cal['n'], total=len(m_inst), zp=numero(cal['zp'], 3), szp=numero(cal['sigma_zp'], 3),
                       rms=numero(cal['rms'], 3), cat=f"{catalogos.NOMBRES[catalogo]} {banda}")
        if con_color:
            texto += '\n' + self.t('cal_color', t=numero(cal['t'], 3), st=numero(cal['sigma_t'], 3))
        ttk.Label(marco, text=texto, font=('Consolas', 10), justify=tk.LEFT).pack(fill=tk.X)
        contenedor, figura, ax, lienzo = self._grafico(marco, 7.2, 3.6)
        contenedor.pack(fill=tk.BOTH, expand=True, pady=(8, 0))
        x = color if con_color else m_cat
        usadas = cal['usadas']
        ax.plot(x[usadas], (m_cat - m_inst)[usadas], 'o', color='#86c8ff', ms=4, label=self.t('cal_usadas'))
        ax.plot(x[~usadas], (m_cat - m_inst)[~usadas], 'x', color='#ff5c5c', ms=6, label=self.t('cal_descartadas'))
        if con_color:
            xs = np.linspace(np.nanmin(x), np.nanmax(x), 20)
            ax.plot(xs, cal['zp'] + cal['t'] * xs, '-', color='#ffd27a')
            ax.set_xlabel(self.t('cal_eje_color'))
        else:
            ax.axhline(cal['zp'], color='#ffd27a')
            ax.set_xlabel(self.t('cal_eje_mag'))
        ax.set_ylabel('m_cat − m_inst')
        leyenda = ax.legend(fontsize=8, facecolor=FONDO, edgecolor=BORDE)
        for t in leyenda.get_texts():
            t.set_color(TEXTO)
        figura.tight_layout()
        lienzo.draw_idle()
        if con_color:
            def usar_t():
                self.var_color_t.set(numero(cal['t'], 3))
                self._registrar('hist_coef_color', t=numero(cal['t'], 3))
                messagebox.showinfo(self.t('cal_titulo'), self.t('cal_t_puesto', t=numero(cal['t'], 3)), parent=v)
            ttk.Button(marco, text=self.t('cal_usar_t'), command=usar_t).pack(anchor=tk.E, pady=(6, 0))
        self._registrar('hist_calibracion', zp=numero(cal['zp'], 3), n=cal['n'])
        self._centrar(v)

    def extincion_serie(self, banda_texto):
        """Coeficiente de extinción con el punto cero de cada foto de la serie y su masa de aire."""
        banda = banda_texto[0].replace("'", '')
        catalogo = 'gaia' if banda == 'G' else 'apass'
        rutas = [r for r in (self.serie if len(self.serie) > 1 else list(self.fotos)) if os.path.exists(separar_ruta(r)[0])]
        masas, zps, errs, nombres = [], [], [], []
        for i, ruta in enumerate(rutas):
            self._mensaje(self.t('serie_leyendo', i=i + 1, n=len(rutas), nombre=self._nombre_ruta(ruta)))
            self.root.update()
            try:
                imagen, cabecera, _n = leer_fits(*separar_ruta(ruta))
            except Exception:
                continue
            masa = C.masa_de_aire(cabecera)
            if masa is None:
                continue
            if ruta not in self.fotos:   # foto de la serie aún no abierta: lo mínimo para calibrarla
                self.fotos[ruta] = {'n': len(self.fotos) + 1, 'archivo': self._nombre_ruta(ruta),
                                    'wcs': SolucionWCS.desde_cabecera(cabecera), 'ancho': imagen.shape[1],
                                    'alto': imagen.shape[0], 'exposicion': exposicion_de(cabecera),
                                    'ganancia': ganancia_de(cabecera), 'saturacion': C.nivel_saturacion(cabecera),
                                    'hora': instante_observacion(cabecera)[0], 'jd': instante_observacion(cabecera)[1]}
            datos = self._fotometria_campo(imagen, ruta, catalogo, banda)
            if isinstance(datos, str):
                continue
            cal = C.calibrar_campo(datos[0], datos[2], None, datos[1])
            if cal:
                masas.append(masa)
                zps.append(cal['zp'])
                errs.append(cal['sigma_zp'])
                nombres.append(self._nombre_ruta(ruta))
        self._mensaje('')
        ext = C.ajustar_extincion(masas, zps, errs)
        if ext is None:
            messagebox.showinfo(self.t('ext_titulo'), self.t('ext_falta', n=len(masas)))
            return
        v, marco = self._ventana(self.t('ext_titulo'), 'v_extincion')
        ttk.Label(marco, text=self.t('ext_resultado', k=numero(ext['k'], 3), sk=numero(ext['sigma_k'], 3), n=ext['n'],
                                     rms=numero(ext['rms'], 3)), font=('Consolas', 10), justify=tk.LEFT).pack(fill=tk.X)
        contenedor, figura, ax, lienzo = self._grafico(marco, 6.8, 3.4)
        contenedor.pack(fill=tk.BOTH, expand=True, pady=(8, 0))
        ax.errorbar(masas, zps, yerr=errs, fmt='o', color='#86c8ff', ms=4, capsize=2)
        xs = np.linspace(min(masas), max(masas), 10)
        ax.plot(xs, ext['zp0'] - ext['k'] * xs, '-', color='#ffd27a')
        ax.set_xlabel(self.t('ext_eje'))
        ax.set_ylabel('ZP')
        figura.tight_layout()
        lienzo.draw_idle()
        self._registrar('hist_extincion', k=numero(ext['k'], 3), n=ext['n'])
        self._centrar(v)

    def radio_optimo(self):
        """Radio de apertura con más señal/ruido para las estrellas medidas en esta foto."""
        if self._necesita_foto():
            return
        medidas = [m for m in self.mediciones if m['ruta'] == self.ruta and 'SAT' not in (m.get('calidad') or [])]
        if not medidas:
            messagebox.showinfo(self.t('ropt_titulo'), self.t('ropt_falta'))
            return
        _r, r_int, r_ext = self._radios_fotometria()
        foto = self.fotos[self.ruta]
        radios = np.arange(1.0, r_int - 0.4, 0.5)
        curvas = [C.curva_crecimiento(self.original, m['x0'], m['y0'], radios, r_int, r_ext, foto.get('ganancia'))
                  for m in medidas]
        mejor = C.radio_optimo(curvas)
        if mejor is None:
            return
        v, marco = self._ventana(self.t('ropt_titulo'), 'v_ropt')
        fwhm = float(np.median([m['fwhm'] for m in medidas]))
        ttk.Label(marco, text=self.t('ropt_resultado', r=numero(mejor, 1), f=numero(mejor / fwhm, 2), fwhm=numero(fwhm, 2)),
                  wraplength=640, justify=tk.LEFT).pack(fill=tk.X)
        contenedor, figura, (ax1, ax2), lienzo = self._grafico(marco, 6.8, 4.4, filas=2)
        contenedor.pack(fill=tk.BOTH, expand=True, pady=(8, 0))
        for m, curva in zip(medidas, curvas):
            if not curva:
                continue
            rs, fl, sn = zip(*curva)
            ax1.plot(rs, np.array(fl) / max(fl), '-', lw=1, label=m['estrella'])
            ax2.plot(rs, np.array(sn) / max(sn), '-', lw=1)
        for ax in (ax1, ax2):
            ax.axvline(mejor, color='#ffd27a', ls='--')
        ax1.set_ylabel(self.t('ropt_flujo'))
        ax2.set_ylabel(self.t('ropt_snr'))
        ax2.set_xlabel(self.t('ropt_radio'))
        figura.tight_layout()
        lienzo.draw_idle()

        def usar():
            self.var_r_fot.set(numero(mejor, 1))
            self.var_r_int.set(numero(max(r_int, mejor + 3), 1))
            self.var_r_ext.set(numero(max(r_ext, mejor + 9), 1))
            self._registrar('hist_radio', r=numero(mejor, 1))
            v.destroy()
        ttk.Button(marco, text=self.t('ropt_usar', r=numero(mejor, 1)), command=usar).pack(anchor=tk.E, pady=(6, 0))
        self._centrar(v)

    # ------------------------------------------------------------------ marcas de catálogo en la foto
    def actualizar_marcas_catalogo(self, ruta, cuales):
        """Pone en la foto las marcas de los catálogos guardados que se piden
        (cuales: conjunto de 'vsx', 'ngc', 'skybot', 'gaia'). No se conecta a nada."""
        foto = self.fotos[ruta]
        wcs, campo = foto.get('wcs'), self._campo_foto(ruta)
        marcas = []
        if wcs is None or campo is None:
            foto['marcas_catalogo'] = marcas
            return 0
        estilos = {'vsx': '#ff9cf0', 'ngc': '#7be0a0', 'skybot': '#ffb36b', 'gaia': '#8fa8ff'}
        for catalogo in ('ngc', 'vsx', 'skybot', 'gaia'):
            if catalogo not in cuales:
                continue
            datos = catalogos.desde_cache(catalogo, *campo, jd=foto.get('jd'))
            if datos is None:
                continue
            filas = datos[0][:150] if catalogo == 'gaia' else datos[0]
            for f in filas:
                x, y = wcs.a_pixel(f['ar'], f['dec'])
                if -20 <= x <= foto['ancho'] + 20 and -20 <= y <= foto['alto'] + 20:
                    if catalogo == 'gaia':
                        texto = f"G={numero(f.get('mag'), 1)}"
                    elif catalogo == 'vsx':
                        texto = f"{f['id']} ({f.get('tipo', '').strip()})"
                    elif catalogo == 'ngc':
                        texto = f"{f['id']} {f.get('tipo', '')}"
                    else:
                        texto = f"{f['id']} V={numero(f.get('mag'), 1)}"
                    marcas.append({'x': x - 1, 'y': y - 1, 'texto': texto, 'color': estilos[catalogo]})
        foto['marcas_catalogo'] = marcas
        self._programar_dibujo()
        return len(marcas)

    # ------------------------------------------------------------------ resolver la foto con Gaia
    def resolver_con_gaia(self, ar, dec, escala, radio):
        """Busca la solución astrométrica de la foto abierta con las estrellas de Gaia
        guardadas (sin conexión). Devuelve el resultado o un texto de error."""
        datos = catalogos.desde_cache('gaia', ar, dec, radio)
        if datos is None:
            return self.t('resolver_sin_gaia')
        filas = sorted(datos[0], key=lambda f: f.get('mag') if f.get('mag') is not None else 99)
        foto = self.fotos[self.ruta]
        epoca = C.epoca_juliana(foto['jd']) if foto.get('jd') else 2016.0
        ars, decs = C.aplicar_movimiento_propio([f['ar'] for f in filas], [f['dec'] for f in filas],
                                                [f.get('pmra') or 0 for f in filas], [f.get('pmdec') or 0 for f in filas],
                                                datos[1].get('epoca') or 2016.0, epoca)
        fuentes = [f for f in (self.fuentes or C.detectar_fuentes(self.original, 5.0, 1500, saturacion=foto.get('saturacion')))
                   if not f['saturada']]
        if len(fuentes) < 5:
            return self.t('ast_pocas_estrellas')
        resultado = C.resolver_campo(fuentes, foto['ancho'], foto['alto'], ars, decs, ar, dec, escala)
        if resultado is None:
            return self.t('resolver_no')
        return resultado

    # ------------------------------------------------------------------ información de la foto
    def ventana_info_foto(self):
        if self._necesita_foto():
            return
        foto = self.fotos[self.ruta]
        v, marco = self._ventana(self.t('info_titulo', nombre=foto['archivo']), 'v_info')
        est = C.estadisticas_imagen(self.original, foto.get('saturacion'))
        datos = [(self.t('info_tamano'), f"{foto['ancho']} × {foto['alto']} px"),
                 (self.t('info_fondo'), numero(est['fondo'], 2)), (self.t('info_ruido'), numero(est['ruido'], 3)),
                 (self.t('info_estrellas'), est['n_fuentes']), ('FWHM', numero(est['fwhm'], 2) + ' px'),
                 (self.t('serie_elong'), numero(est['elongacion'], 2)),
                 (self.t('info_saturacion'), numero(foto.get('saturacion'), 0) if foto.get('saturacion') else '—'),
                 (self.t('info_saturadas'), sum(1 for f in est['fuentes'] if f['saturada'])),
                 (self.t('info_escala'), (numero(foto['wcs'].escala, 3) + '″/px') if foto.get('wcs') else '—'),
                 (self.t('info_masa'), numero(foto.get('airmass'), 3) if foto.get('airmass') else '—'),
                 (self.t('col_hora'), foto.get('hora') or '—')]
        resumen = ttk.Frame(marco)
        resumen.pack(fill=tk.X)
        for i, (clave, valor) in enumerate(datos):
            ttk.Label(resumen, text=clave + ':', style='Gris.TLabel').grid(row=i // 2, column=(i % 2) * 2, sticky=tk.W, padx=(0, 6))
            ttk.Label(resumen, text=str(valor), font=('Segoe UI', 10, 'bold')).grid(row=i // 2, column=(i % 2) * 2 + 1,
                                                                                   sticky=tk.W, padx=(0, 24))
        contenedor, figura, ax, lienzo = self._grafico(marco, 7.4, 2.6)
        contenedor.pack(fill=tk.X, pady=(8, 4))
        valores = self.original[np.isfinite(self.original)].ravel()
        if valores.size > 2_000_000:
            valores = valores[::max(1, valores.size // 2_000_000)]
        lo, hi = np.percentile(valores, [0.5, 99.9])
        ax.hist(valores, bins=200, range=(lo, hi), color='#86c8ff', log=True)
        ax.axvline(est['fondo'], color='#ffd27a', ls='--')
        ax.set_xlabel(self.t('info_histograma'))
        figura.tight_layout()
        lienzo.draw_idle()
        # cabecera FITS, con buscador
        fila = ttk.Frame(marco)
        fila.pack(fill=tk.X, pady=(4, 2))
        ttk.Label(fila, text=self.t('info_buscar')).pack(side=tk.LEFT)
        var = tk.StringVar()
        ttk.Entry(fila, textvariable=var, width=24).pack(side=tk.LEFT, padx=6)
        marco_tabla, tabla = self._tabla_simple(marco, [(self.t('info_clave'), 120, tk.W), (self.t('conf_valor'), 520, tk.W)],
                                                alto=10)
        marco_tabla.pack(fill=tk.BOTH, expand=True)
        cabecera = getattr(self, 'cabecera', {}) or {}

        def llenar(*_):
            filtro = var.get().strip().upper()
            tabla.delete(*tabla.get_children())
            for clave, valor in cabecera.items():
                if clave.startswith('_'):
                    continue
                if filtro and filtro not in clave and filtro not in str(valor).upper():
                    continue
                tabla.insert('', tk.END, values=(clave, valor))
        var.trace_add('write', llenar)
        llenar()
        self._centrar(v)

    # ------------------------------------------------------------------ guía paso a paso
    def ventana_guia(self):
        """Lista de pasos para empezar, que se marca sola según lo que ya está hecho."""
        v, marco = self._ventana(self.t('guia_titulo'), 'v_guia')
        ttk.Label(marco, text=self.t('guia_explicacion'), wraplength=560, style='Gris.TLabel', justify=tk.LEFT).pack(fill=tk.X)
        lista = ttk.Frame(marco)
        lista.pack(fill=tk.X, pady=8)

        def estado():
            for w in lista.winfo_children():
                w.destroy()
            hay_foto = self.ruta is not None
            medidas = len(self.mediciones)
            fotos_medidas = len({m['ruta'] for m in self.mediciones})
            con_wcs = any(f.get('wcs') for f in self.fotos.values())
            pasos = [(hay_foto, 'guia_1'), (medidas > 0, 'guia_2'), (fotos_medidas > 1 or len(self.serie) < 2, 'guia_3'),
                     (all('calidad' in m and not (set(m['calidad']) & {'SAT', 'CON'}) for m in self.mediciones) and medidas > 0,
                      'guia_4'),
                     (self._fotometria_activa(), 'guia_5'), (con_wcs, 'guia_6'),
                     (any(f.get('rms_astrometrico') is not None for f in self.fotos.values()), 'guia_7'),
                     (bool(getattr(self, '_guardado', False)), 'guia_8')]
            for i, (hecho, clave) in enumerate(pasos):
                ttk.Label(lista, text='✓' if hecho else '○', foreground='#7be0a0' if hecho else TEXTO_SUAVE,
                          font=('Segoe UI', 13, 'bold')).grid(row=i, column=0, sticky=tk.N, padx=(0, 8), pady=3)
                ttk.Label(lista, text=self.t(clave), wraplength=520, justify=tk.LEFT,
                          foreground=TEXTO_SUAVE if hecho else TEXTO).grid(row=i, column=1, sticky=tk.W, pady=3)
        ttk.Button(marco, text=self.t('fot_actualizar'), command=estado).pack(anchor=tk.E)
        estado()
        self._centrar(v)

    # ------------------------------------------------------------------ modo sencillo / avanzado
    def _aplicar_modo(self):
        """Modo sencillo: solo lo esencial (se esconden el análisis y las opciones avanzadas)."""
        sencillo = self.ajustes.get('modo') == 'sencillo'
        try:
            for widget, opciones in getattr(self, '_widgets_avanzados', []):
                if sencillo and widget.winfo_ismapped():
                    widget.pack_forget()
                elif not sencillo and not widget.winfo_ismapped():
                    widget.pack(**opciones)
            if sencillo and self.avanzadas_visibles:
                self._alternar_avanzadas()
        except tk.TclError:
            pass

    # ------------------------------------------------------------------ informe PDF
    def guardar_pdf(self):
        if not self.mediciones:
            messagebox.showinfo(self.t('guardar_t'), self.t('sin_medidas'))
            return
        ruta = filedialog.asksaveasfilename(title=self.t('pdf_titulo'), defaultextension='.pdf',
                                            initialfile=self.t('pdf_archivo') + '.pdf', filetypes=[('PDF', '*.pdf')])
        if not ruta:
            return
        try:
            self.escribir_pdf(ruta)
        except PermissionError:
            messagebox.showerror(self.t('error_guardar'), self.t('abierto_otro'))
            return
        except Exception as error:
            messagebox.showerror(self.t('error_guardar'), str(error))
            return
        self._registrar('hist_pdf', ruta=ruta)
        messagebox.showinfo(self.t('guardado_t'), self.t('mpc_guardado', ruta=ruta))

    def escribir_pdf(self, ruta):
        """Informe en PDF (con matplotlib): resumen, tabla de medidas, curvas de luz,
        objetos móviles y revisión de la serie."""
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_pdf import PdfPages
        a4 = (8.27, 11.69)

        def pagina(titulo):
            fig = Figure(figsize=a4)
            fig.text(0.07, 0.955, titulo, fontsize=15, weight='bold', color='#1a2444')
            fig.text(0.07, 0.94, f'Emillium Expanded {VERSION} · ' + datetime.now().strftime('%Y-%m-%d %H:%M'),
                     fontsize=8, color='#666666')
            return fig
        with PdfPages(ruta) as pdf:
            # 1) resumen
            fig = pagina(self.t('pdf_resumen'))
            lineas = []
            for r, foto in sorted(self.fotos.items(), key=lambda par: par[1]['n']):
                n = sum(1 for m in self.mediciones if m['ruta'] == r)
                if not n:
                    continue
                linea = f"{foto['n']:>3}. {foto['archivo'][:42]:<42} {foto.get('hora') or '':<24} {n} {self.t('pdf_medidas')}"
                if foto.get('rms_astrometrico') is not None:
                    linea += f"  rms {foto['rms_astrometrico']:.2f}″"
                lineas.append(linea)
            parametros = self.tabla_configuracion()[1]
            texto = '\n'.join(lineas[:45]) + '\n\n' + '\n'.join(f'{p}: {v}' for p, v in parametros[:12])
            fig.text(0.07, 0.9, texto, fontsize=7.5, family='monospace', va='top')
            pdf.savefig(fig)
            # 2) medidas (tabla, varias páginas)
            columnas, filas = self.tabla_todas()
            elegidas = [0, 2, 4, 5, 6, 7, 8, 15, 18, 21, 22]
            cabecera = [columnas[i][0] for i in elegidas]
            for inicio in range(0, len(filas), 42):
                fig = pagina(self.t('pdf_tabla'))
                ax = fig.add_axes([0.04, 0.04, 0.92, 0.88])
                ax.axis('off')
                celdas = [[('' if fila[i] is None else (numero(fila[i], d) if isinstance(fila[i], float) and d is not None else str(fila[i])))
                           for i, d in ((i, columnas[i][1]) for i in elegidas)] for fila in filas[inicio:inicio + 42]]
                tabla = ax.table(cellText=celdas, colLabels=cabecera, loc='upper center', cellLoc='right')
                tabla.auto_set_font_size(False)
                tabla.set_fontsize(5.6)
                tabla.scale(1, 1.15)
                pdf.savefig(fig)
            # 3) curvas de luz
            series = self._serie_fotometrica()
            if series and len({m['ruta'] for m in self.mediciones}) > 1:
                fig = pagina(self.t('pdf_curvas'))
                ax = fig.add_axes([0.1, 0.5, 0.85, 0.38])
                rutas = sorted({r for l in series for r, _m, _e in series[l]}, key=lambda r: self.fotos[r]['n'])
                eje, titulo_x = self._eje_tiempo(rutas)
                for letra in sorted(series, key=LETRAS_ESTRELLA.index)[:10]:
                    datos = series[letra]
                    ax.errorbar([eje[r] for r, _m, _e in datos], [m for _r, m, _e in datos],
                                yerr=[e or 0 for _r, _m, e in datos], fmt='o-', ms=3, lw=0.8, capsize=1.5, label=letra)
                ax.invert_yaxis()
                ax.set_xlabel(titulo_x)
                ax.set_ylabel(self.tabla.heading('mag')['text'])
                ax.legend(fontsize=7, ncol=5)
                ax.grid(alpha=0.3)
                texto = []
                for letra in sorted(series, key=LETRAS_ESTRELLA.index):
                    info = C.analisis_variabilidad([d[1] for d in series[letra]], [d[2] for d in series[letra]])
                    if info['n'] >= 3:
                        texto.append(f"{letra}: n={info['n']}  {self.t('fot_media')} {info['media']:.3f}  σ {info['desviacion']:.3f}  "
                                     f"χ² {info['chi2_red']:.1f}  {'VARIABLE' if info['variable'] else ''}")
                fig.text(0.1, 0.44, '\n'.join(texto[:30]), fontsize=7.5, family='monospace', va='top')
                pdf.savefig(fig)
            # 4) objetos móviles
            moviles = [d for d in self.analisis_moviles() if d.get('se_mueve')]
            if moviles:
                fig = pagina(self.t('moviles_titulo'))
                ax = fig.add_axes([0.12, 0.45, 0.8, 0.45])
                texto = []
                for d in moviles:
                    medidas = d['medidas']
                    if d.get('unidad') == '″':
                        xi, eta = C.proyectar(np.array([m['ar'] for m in medidas]), np.array([m['dec'] for m in medidas]),
                                              medidas[0]['ar'], medidas[0]['dec'])
                        ax.plot(-xi * 3600, eta * 3600, 'o-', ms=3, label=d['letra'])
                    texto.append(f"{d['letra']}: {d['n']} {self.t('pdf_fotos')} · {numero(d.get('velocidad'), 2)} {d.get('unidad', '')}/h · "
                                 f"{numero(d.get('angulo'), 1)}° · rms {numero(d.get('rms'), 3)}")
                ax.set_xlabel(self.t('moviles_eje_este'))
                ax.set_ylabel(self.t('moviles_eje_norte'))
                ax.legend(fontsize=7)
                ax.grid(alpha=0.3)
                ax.set_aspect('equal', adjustable='datalim')
                fig.text(0.12, 0.38, '\n'.join(texto), fontsize=8, family='monospace', va='top')
                pdf.savefig(fig)
            # 5) revisión de la serie
            if self.revision:
                validos = [r for r in self.revision if 'error' not in r]
                if validos:
                    fig = pagina(self.t('hoja_revision'))
                    ax1 = fig.add_axes([0.12, 0.62, 0.8, 0.26])
                    ax2 = fig.add_axes([0.12, 0.3, 0.8, 0.26])
                    x = list(range(1, len(validos) + 1))
                    ax1.plot(x, [r['fwhm'] for r in validos], 'o-', ms=3)
                    ax1.set_ylabel('FWHM (px)')
                    ax2.plot(x, [r['fondo'] for r in validos], 'o-', ms=3, color='tab:orange')
                    ax2.set_ylabel(self.t('serie_fondo'))
                    ax2.set_xlabel(self.t('eje_foto'))
                    for ax in (ax1, ax2):
                        ax.grid(alpha=0.3)
                    malas = [f"{i}. {self._nombre_ruta(r['ruta'])}: {', '.join(self.t('motivo_' + m) for m in r['motivos'])}"
                             for i, r in enumerate(validos, 1) if r.get('motivos')]
                    fig.text(0.12, 0.24, '\n'.join(malas[:20]) or '✓', fontsize=8, va='top')
                    pdf.savefig(fig)
            info = pdf.infodict()
            info['Title'] = self.t('pdf_resumen')
            info['Creator'] = f'Emillium Expanded {VERSION}'


    # =================================================================== HERRAMIENTAS DE LA 1.5 (III)
    # ------------------------------------------------------------------ errores inesperados
    def _error_inesperado(self, tipo, valor, traza):
        """Si algo falla dentro de un botón, el programa no se cierra: se avisa y se guarda el
        detalle en %APPDATA%\\Emillium\\errores.log (para poder arreglarlo)."""
        import traceback
        texto = ''.join(traceback.format_exception(tipo, valor, traza))
        try:
            carpeta = os.path.dirname(_ruta_ajustes())
            os.makedirs(carpeta, exist_ok=True)
            with open(os.path.join(carpeta, 'errores.log'), 'a', encoding='utf-8') as archivo:
                archivo.write(f'--- {datetime.now().isoformat(timespec="seconds")} · Emillium {VERSION}\n{texto}\n')
        except OSError:
            pass
        try:
            self.root.config(cursor='')
            messagebox.showerror(TITULO_VENTANA, self.t('error_inesperado', error=f'{tipo.__name__}: {valor}'))
        except Exception:
            pass

    # ------------------------------------------------------------------ tiempo heliocéntrico
    def _hjd_foto(self, ruta):
        """HJD del centro de la exposición de una foto (para el centro del campo o, sin solución
        astrométrica, para donde apuntaba el telescopio). None si falta algo."""
        return self._tiempo_foto(ruta, 'HJD')

    def _bjd_foto(self, ruta):
        """BJD_TDB del centro de la exposición de una foto (o None)."""
        return self._tiempo_foto(ruta, 'BJD_TDB')

    def _tiempo_foto(self, ruta, escala):
        foto = self.fotos.get(ruta, {})
        if not foto.get('jd'):
            return None
        campo = self._campo_foto(ruta) if foto.get('wcs') else None
        if campo:
            ar, dec = campo[0], campo[1]
        else:
            apuntado = foto.get('apuntado') or {}
            ar, dec = apuntado.get('ar'), apuntado.get('dec')
        if ar is None or dec is None:
            return None
        if escala == 'BJD_TDB':
            return C.bjd_tdb(foto['jd'], ar, dec)
        return C.jd_heliocentrico(foto['jd'], ar, dec)

    def _selector_tiempo(self, fila):
        """Desplegable con la escala de tiempo (BJD_TDB, HJD o JD). Devuelve su variable."""
        variable = tk.StringVar(value=self.ajustes.get('escala_tiempo', 'BJD_TDB'))
        ttk.Label(fila, text=self.t('escala_tiempo')).pack(side=tk.LEFT, padx=(8, 2))
        ttk.Combobox(fila, textvariable=variable, values=ESCALAS_TIEMPO, state='readonly', width=9).pack(side=tk.LEFT)

        def recordar(*_):
            self.ajustes['escala_tiempo'] = variable.get()
            self._guardar_ajustes()
        variable.trace_add('write', recordar)
        return variable

    def _tiempos_curva(self, rutas, escala):
        """Instante de cada foto en la escala pedida ('BJD_TDB', 'HJD' o 'JD'; si falta algo
        para calcularla, la siguiente que se pueda). Devuelve ({ruta: tiempo}, nombre del eje)."""
        if escala is True:
            escala = 'HJD'
        if escala in ('BJD_TDB', 'HJD'):
            t = {r: self._tiempo_foto(r, escala) for r in rutas}
            if all(v is not None for v in t.values()) and t:
                return t, escala
        t = {r: self.fotos[r].get('jd') for r in rutas}
        if all(v is not None for v in t.values()) and t:
            return t, 'JD'
        return {r: float(self.fotos[r]['n']) for r in rutas}, self.t('eje_foto')

    # ------------------------------------------------------------------ periodo y mínimos
    def ventana_periodo(self, letra=None):
        """Periodograma de Lomb-Scargle de una estrella y su curva de luz en fase."""
        series = self._serie_fotometrica()
        letras = [l for l in sorted(series, key=LETRAS_ESTRELLA.index) if len(series[l]) >= 6]
        if not letras:
            messagebox.showinfo(self.t('per_titulo'), self.t('per_falta'))
            return
        v, marco = self._ventana(self.t('per_titulo'), 'v_periodo')
        ttk.Label(marco, text=self.t('per_explicacion'), wraplength=760, style='Gris.TLabel', justify=tk.LEFT).pack(fill=tk.X)
        fila = ttk.Frame(marco)
        fila.pack(fill=tk.X, pady=6)
        var_letra = tk.StringVar(value=letra if letra in letras else letras[0])
        var_min, var_max = tk.StringVar(value=''), tk.StringVar(value='')
        ttk.Label(fila, text=self.t('col_estrella') + ':').pack(side=tk.LEFT)
        ttk.Combobox(fila, textvariable=var_letra, values=letras, state='readonly', width=4).pack(side=tk.LEFT, padx=(4, 12))
        ttk.Label(fila, text=self.t('per_min')).pack(side=tk.LEFT)
        ttk.Entry(fila, textvariable=var_min, width=8).pack(side=tk.LEFT, padx=(4, 8))
        ttk.Label(fila, text=self.t('per_max')).pack(side=tk.LEFT)
        ttk.Entry(fila, textvariable=var_max, width=8).pack(side=tk.LEFT, padx=(4, 8))
        var_hjd = self._selector_tiempo(fila)
        resultado = ttk.Label(marco, text='', font=('Consolas', 10, 'bold'), foreground=DORADO)
        resultado.pack(fill=tk.X)
        contenedor, figura, (ax1, ax2), lienzo = self._grafico(marco, 7.6, 5.0, filas=2)
        contenedor.pack(fill=tk.BOTH, expand=True, pady=(6, 0))

        def calcular(*_):
            datos = series[var_letra.get()]
            tiempos, nombre = self._tiempos_curva([d[0] for d in datos], var_hjd.get())
            t = np.array([tiempos[d[0]] for d in datos], dtype=np.float64)
            m = np.array([d[1] for d in datos])
            e = np.array([d[2] or 0.01 for d in datos])
            pmin = leer_angulo(var_min.get()) if var_min.get().strip() else None
            pmax = leer_angulo(var_max.get()) if var_max.get().strip() else None
            r = C.periodograma(t, m, e, pmin, pmax)
            for ax in (ax1, ax2):
                ax.clear()
                self._estilo_ejes(ax)
            if r is None or nombre not in ESCALAS_TIEMPO:
                resultado.config(text=self.t('per_sin_tiempo'))
                lienzo.draw_idle()
                return
            periodos, potencia, mejor = r
            ax1.plot(periodos, potencia, '-', color='#86c8ff', lw=0.8)
            ax1.axvline(mejor, color='#ffd27a', ls='--')
            ax1.set_xscale('log')
            ax1.set_xlabel(self.t('per_eje_periodo'))
            ax1.set_ylabel(self.t('per_eje_potencia'))
            fases = C.fase(t, mejor)
            for desplazamiento in (0, 1):
                ax2.errorbar(fases + desplazamiento, m, yerr=e, fmt='o', ms=3, color='#ffd27a', alpha=0.85 if desplazamiento == 0 else 0.35)
            ax2.invert_yaxis()
            ax2.set_xlabel(self.t('per_eje_fase'))
            ax2.set_ylabel(self.tabla.heading('mag')['text'])
            resultado.config(text=self.t('per_resultado', letra=var_letra.get(), p=f'{mejor:.6f}', h=numero(mejor * 24, 3),
                                         pot=numero(float(potencia.max()), 2), tiempo=nombre))
            figura.tight_layout()
            lienzo.draw_idle()
            self._registrar('hist_periodo', letra=var_letra.get(), p=f'{mejor:.6f}')
        ttk.Button(fila, text=self.t('dobles_calcular'), style='Grande.TButton', command=calcular).pack(side=tk.LEFT, padx=6)
        var_letra.trace_add('write', calcular)
        calcular()
        self._centrar(v)

    def ventana_minimo(self):
        """Instante del mínimo de un eclipse (Kwee-van Woerden) de la estrella elegida."""
        series = self._serie_fotometrica()
        letras = [l for l in sorted(series, key=LETRAS_ESTRELLA.index) if len(series[l]) >= 7]
        if not letras:
            messagebox.showinfo(self.t('min_titulo'), self.t('min_falta'))
            return
        v, marco = self._ventana(self.t('min_titulo'), 'v_minimo')
        ttk.Label(marco, text=self.t('min_explicacion'), wraplength=720, style='Gris.TLabel', justify=tk.LEFT).pack(fill=tk.X)
        fila = ttk.Frame(marco)
        fila.pack(fill=tk.X, pady=6)
        var_letra = tk.StringVar(value=letras[0])
        ttk.Label(fila, text=self.t('col_estrella') + ':').pack(side=tk.LEFT)
        ttk.Combobox(fila, textvariable=var_letra, values=letras, state='readonly', width=4).pack(side=tk.LEFT, padx=(4, 12))
        var_hjd = self._selector_tiempo(fila)
        resultado = ttk.Label(marco, text='', font=('Consolas', 10, 'bold'), foreground=DORADO, justify=tk.LEFT)
        resultado.pack(fill=tk.X)
        contenedor, figura, ax, lienzo = self._grafico(marco, 7.2, 3.4)
        contenedor.pack(fill=tk.BOTH, expand=True, pady=(6, 0))

        def calcular(*_):
            datos = series[var_letra.get()]
            tiempos, nombre = self._tiempos_curva([d[0] for d in datos], var_hjd.get())
            t = np.array([tiempos[d[0]] for d in datos], dtype=np.float64)
            m = np.array([d[1] for d in datos])
            ax.clear()
            self._estilo_ejes(ax)
            r = C.minimo_kvw(t, m) if nombre in ESCALAS_TIEMPO else None
            t0 = float(t.min())
            ax.plot((t - t0) * 24, m, 'o', ms=4, color='#86c8ff')
            ax.invert_yaxis()
            ax.set_xlabel(self.t('min_eje', t=nombre, t0=f'{t0:.5f}'))
            ax.set_ylabel(self.tabla.heading('mag')['text'])
            if r is None:
                resultado.config(text=self.t('min_no'), foreground=ROJO_SUAVE)
            else:
                ax.axvline((r[0] - t0) * 24, color='#ffd27a', ls='--')
                resultado.config(text=self.t('min_resultado', tiempo=nombre, t=f'{r[0]:.6f}', e=numero(r[1] * 86400, 1),
                                             fecha=fecha_desde_jd(r[0])), foreground=DORADO)
                self._registrar('hist_minimo', letra=var_letra.get(), t=f'{r[0]:.6f}', tipo=nombre)
            figura.tight_layout()
            lienzo.draw_idle()
        var_letra.trace_add('write', calcular)
        var_hjd.trace_add('write', calcular)
        calcular()
        self._centrar(v)

    # ------------------------------------------------------------------ galería de recortes
    def ventana_recortes(self):
        """La misma estrella en todas las fotos donde se ha medido, una al lado de otra:
        para ver de un vistazo si alguna medida es mala (nube, rayo cósmico, estrella movida...)."""
        letras = self._letras_medidas()
        if not letras:
            messagebox.showinfo(TITULO_VENTANA, self.t('sin_medidas'))
            return
        v, marco = self._ventana(self.t('rec_titulo'), 'v_recortes')
        fila = ttk.Frame(marco)
        fila.pack(fill=tk.X)
        var_letra = tk.StringVar(value=self.var_estrella.get() if self.var_estrella.get() in letras else letras[0])
        ttk.Label(fila, text=self.t('col_estrella') + ':').pack(side=tk.LEFT)
        ttk.Combobox(fila, textvariable=var_letra, values=letras, state='readonly', width=4).pack(side=tk.LEFT, padx=4)
        contenedor, figura, _ax, lienzo = self._grafico(marco, 8.4, 5.6)
        contenedor.pack(fill=tk.BOTH, expand=True, pady=(6, 0))
        cache_imagenes = {}

        def dibujar(*_):
            medidas = sorted((m for m in self.mediciones if m['estrella'] == var_letra.get()), key=lambda m: m['foto_n'])[:30]
            figura.clear()
            n = len(medidas)
            columnas = min(6, max(n, 1))
            filas = max(1, math.ceil(n / columnas))
            r_px = int(max(8, 3 * max((m.get('r_fot') or 5) for m in medidas))) if medidas else 10
            for i, m in enumerate(medidas):
                ruta = m['ruta']
                if ruta == self.ruta and self.original is not None:
                    imagen = self.original
                elif ruta in cache_imagenes:
                    imagen = cache_imagenes[ruta]
                else:
                    try:
                        imagen = leer_fits(*separar_ruta(ruta))[0]
                    except Exception:
                        imagen = None
                    if len(cache_imagenes) < 40:   # como mucho 40 fotos en memoria
                        cache_imagenes[ruta] = imagen
                ax = figura.add_subplot(filas, columnas, i + 1)
                ax.set_xticks([])
                ax.set_yticks([])
                for borde in ax.spines.values():
                    borde.set_color('#ffb36b' if set(m.get('calidad') or []) - {'BOR'} else BORDE)
                if imagen is None:
                    ax.set_facecolor(CAMPO)
                    continue
                trozo = C.recorte(imagen, m['x0'], m['y0'], r_px)
                finitos = trozo[np.isfinite(trozo)]
                lo, hi = (np.percentile(finitos, [5, 99.5]) if finitos.size else (0, 1))
                ax.imshow(trozo, origin='lower', cmap='gray', vmin=lo, vmax=hi if hi > lo else lo + 1, interpolation='nearest')
                ax.plot([r_px], [r_px], '+', color='#ffd27a', ms=8)
                ax.set_title(f"{m['foto_n']} · {' '.join(m.get('calidad') or []) or '✓'}", fontsize=7, color=TEXTO)
            figura.tight_layout()
            lienzo.draw_idle()
        var_letra.trace_add('write', dibujar)
        dibujar()
        self._centrar(v)

    # ------------------------------------------------------------------ gestor de catálogos guardados
    def ventana_gestor_catalogos(self, padre=None):
        v, marco = self._ventana(self.t('gestor_titulo'), 'v_gestor')
        resumen = ttk.Label(marco, text='', style='Gris.TLabel')
        resumen.pack(fill=tk.X)
        columnas = [(self.t('gestor_catalogo'), 140, tk.W), ('AR (°)', 90, tk.E), ('Dec (°)', 90, tk.E),
                    (self.t('gestor_radio'), 80, tk.E), (self.t('gestor_objetos'), 80, tk.E),
                    (self.t('gestor_fecha'), 130, tk.W), ('KB', 70, tk.E)]
        marco_tabla, tabla = self._tabla_simple(marco, columnas, alto=12)
        marco_tabla.pack(fill=tk.BOTH, expand=True, pady=(6, 0))
        rutas = {}

        def llenar():
            tabla.delete(*tabla.get_children())
            rutas.clear()
            for i, (ruta, meta) in enumerate(catalogos.listar_cache()):
                rutas[str(i)] = ruta
                tabla.insert('', tk.END, iid=str(i), values=(
                    meta.get('nombre', meta.get('catalogo')), numero(meta.get('ar'), 3), numero(meta.get('dec'), 3),
                    numero(meta.get('radio'), 3), meta.get('n'),
                    datetime.fromtimestamp(meta.get('fecha', 0)).strftime('%Y-%m-%d %H:%M'),
                    os.path.getsize(ruta) // 1024))
            resumen.config(text=self.t('gestor_resumen', n=len(rutas), mb=numero(catalogos.tamano_cache() / 1e6, 1),
                                       carpeta=catalogos.carpeta_cache()))

        def borrar():
            seleccion = tabla.selection()
            if not seleccion:
                return
            if not messagebox.askyesno(self.t('gestor_titulo'), self.t('gestor_confirmar', n=len(seleccion)), parent=v,
                                       icon=messagebox.WARNING, default=messagebox.NO):
                return
            for iid in seleccion:
                try:
                    catalogos.borrar_de_cache(rutas[iid])
                except (OSError, ValueError) as error:
                    messagebox.showerror(self.t('gestor_titulo'), str(error), parent=v)
            self._registrar('hist_borrar_cat', n=len(seleccion))
            llenar()
        ttk.Button(marco, text=self.t('gestor_borrar'), style='Peligro.TButton', command=borrar).pack(anchor=tk.E, pady=(6, 0))
        llenar()
        self._centrar(v)


def main():
    if sys.platform == 'win32':
        try:  # texto nítido en pantallas con escalado (125 %, 150 %...)
            import ctypes
            ctypes.windll.shcore.SetProcessDpiAwareness(1)
        except Exception:
            pass
    root = tk.Tk()
    app = EmilliumApp(root)
    ruta_icono = resource_path('icono.ico')
    if os.path.exists(ruta_icono):
        try:
            root.iconbitmap(ruta_icono)
        except Exception:
            pass
    # Si se abre con un archivo (por ejemplo arrastrándolo sobre el .exe), cargarlo.
    # (si se arrastran varias fotos a la vez, se abren como una serie)
    archivos = [a for a in sys.argv[1:] if os.path.isfile(a)]
    if archivos:
        root.after(100, lambda: app.abrir_serie(archivos))
    root.mainloop()


if __name__ == '__main__':
    main()

