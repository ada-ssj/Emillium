"""
Emillium v1.3 - Medidor de centroides estelares en imágenes FITS
=================================================================

Programa de escritorio para abrir imágenes astronómicas en formato FITS,
girarlas o invertirlas como un espejo, y medir la posición de las estrellas
en PÍXELES con precisión subpíxel (centroide).

Solo utiliza:
  - La biblioteca estándar de Python (tkinter para la ventana, csv, math...).
  - numpy, para hacer los cálculos con matrices de forma rápida.

No usa astropy, photutils ni ninguna solución WCS: el lector de archivos
FITS, el giro de la imagen y el cálculo del centroide están programados
aquí mismo, paso a paso.

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
import math
import os
import base64
import random
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

VERSION = '1.3'


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


def leer_fits(ruta):
    """Lee la primera imagen 2D de un archivo FITS.

    Devuelve (imagen, cabecera, nota):
      imagen   -> matriz numpy float32 con forma (alto, ancho); la fila 0 es
                  la de abajo.
      cabecera -> diccionario con las palabras clave.
      nota     -> texto con avisos (por ejemplo, si era una imagen en color).
    """
    with open(ruta, 'rb') as archivo:
        if archivo.read(6) != b'SIMPLE':
            raise ErrorFits('err_no_fits')
        archivo.seek(0)

        comprimido = False
        while True:
            cabecera = _leer_cabecera(archivo)
            if cabecera is None:
                break

            bitpix = int(cabecera.get('BITPIX', 0))
            naxis = int(cabecera.get('NAXIS', 0))
            ejes = [int(cabecera.get(f'NAXIS{i}', 0)) for i in range(1, naxis + 1)]
            num_valores = int(np.prod(ejes)) if naxis > 0 else 0
            pcount = int(cabecera.get('PCOUNT', 0))
            gcount = int(cabecera.get('GCOUNT', 1))
            bytes_por_valor = abs(bitpix) // 8
            bytes_datos = bytes_por_valor * gcount * (pcount + num_valores) if naxis > 0 else 0
            relleno = (-bytes_datos) % TAMANO_BLOQUE

            if cabecera.get('ZIMAGE') is True:
                comprimido = True

            extension = cabecera.get('XTENSION')
            es_imagen = extension is None or str(extension).strip() == 'IMAGE'
            if es_imagen and naxis >= 2 and num_valores > 0:
                if bitpix not in TIPOS_BITPIX:
                    raise ErrorFits('err_bitpix', bitpix=bitpix)
                crudo = archivo.read(bytes_por_valor * num_valores)
                if len(crudo) < bytes_por_valor * num_valores:
                    raise ErrorFits('err_datos')
                datos = np.frombuffer(crudo, dtype=TIPOS_BITPIX[bitpix])
                datos = datos.reshape(ejes[::-1])  # numpy: (planos, alto, ancho)
                imagen, nota = _convertir_datos(datos, cabecera, bitpix)
                return imagen, cabecera, nota

            # Este HDU no es una imagen: saltamos sus datos y seguimos.
            archivo.seek(bytes_datos + relleno, 1)

    if comprimido:
        raise ErrorFits('err_comprimido')
    raise ErrorFits('err_sin_imagen')


def _convertir_datos(datos, cabecera, bitpix):
    """Aplica BZERO/BSCALE/BLANK y deja una imagen 2D en float32."""
    bzero = float(cabecera.get('BZERO', 0.0))
    bscale = float(cabecera.get('BSCALE', 1.0))
    blank = cabecera.get('BLANK')

    sin_dato = None
    if bitpix > 0 and isinstance(blank, int):
        sin_dato = datos == blank

    imagen = datos.astype(np.float64)
    if bscale != 1.0 or bzero != 0.0:
        imagen = imagen * bscale + bzero
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
    dx = (np.arange(nuevo_ancho) - ncx)[None, :]
    filas_por_trozo = max(1, 500_000 // nuevo_ancho)  # por trozos, para no gastar mucha memoria

    for fila_ini in range(0, nuevo_alto, filas_por_trozo):
        fila_fin = min(fila_ini + filas_por_trozo, nuevo_alto)
        dy = (np.arange(fila_ini, fila_fin) - ncy)[:, None]
        # Giro inverso: punto de la imagen original del que viene cada píxel nuevo.
        x_orig = cx + c * dx + s * dy
        y_orig = cy - s * dx + c * dy
        dentro = (x_orig >= 0) & (x_orig <= ancho - 1) & (y_orig >= 0) & (y_orig <= alto - 1)

        x0 = np.clip(np.floor(x_orig).astype(np.int64), 0, ancho - 2)
        y0 = np.clip(np.floor(y_orig).astype(np.int64), 0, alto - 2)
        fx = x_orig - x0
        fy = y_orig - y0
        valor = (imagen[y0, x0] * (1 - fx) * (1 - fy)
                 + imagen[y0, x0 + 1] * fx * (1 - fy)
                 + imagen[y0 + 1, x0] * (1 - fx) * fy
                 + imagen[y0 + 1, x0 + 1] * fx * fy)
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
                   umbral_sigma=3.0):
    """Mide la estrella más brillante cerca de (x_clic, y_clic).

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
        otros_picos = [(x0 + c, y0 + f) for f, c in zip(filas, columnas) if (f, c) != (fila, columna)]
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
    }


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
        'ayuda_abrir': 'Obre una imatge en format FITS (.fits, .fit, .fts). També: Ctrl+O.',
        'ninguna_foto': 'Cap foto oberta',
        'info_foto': '{nombre}\nFoto n.º {n} · {ancho} × {alto} píxels',
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
        'detalle': 'Flux: {flujo}\nMida (FWHM): {fwhm} px\nPic: {ps} vegades el soroll del cel',
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
        'col_hora': 'Data i hora (DATE-OBS)',
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
        'nota_color': 'Imatge en color (3 plans): es fa servir la mitjana dels tres.',
        'nota_planos': 'La imatge té diversos plans: es fa servir el primer.',
    },
    'es': {
        'subtitulo': 'Medidor de centroides estelares',
        'barra_ayuda': 'Clic izquierdo: medir estrella   ·   Rueda: zoom   ·   '
                       'Arrastrar con botón derecho: mover la foto',
        'paso1': ' 1. Abre una foto ',
        'abrir': 'Abrir foto FITS…',
        'ayuda_abrir': 'Abre una imagen en formato FITS (.fits, .fit, .fts). También: Ctrl+O.',
        'ninguna_foto': 'Ninguna foto abierta',
        'info_foto': '{nombre}\nFoto n.º {n} · {ancho} × {alto} píxeles',
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
        'detalle': 'Flujo: {flujo}\nTamaño (FWHM): {fwhm} px\nPico: {ps} veces el ruido del cielo',
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
        'col_hora': 'Fecha y hora (DATE-OBS)',
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
        'nota_color': 'Imagen en color (3 planos): se usa la media de los tres.',
        'nota_planos': 'La imagen tiene varios planos: se usa el primero.',
    },
}

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
    ('foto_n', 'col_foto_n', 70, tk.CENTER),
    ('archivo', 'col_foto', 200, tk.W),
    ('estrella', 'col_estrella', 70, tk.CENTER),
    ('orientacion', 'col_orientacion', 280, tk.W),
    ('x', None, 100, tk.E),
    ('y', None, 100, tk.E),
    ('flujo', 'col_flujo', 120, tk.E),
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

def _letra_columna(indice):
    """0 -> A, 1 -> B, ..., 26 -> AA (nombres de columna de Excel)."""
    letras = ''
    indice += 1
    while indice:
        indice, resto = divmod(indice - 1, 26)
        letras = chr(65 + resto) + letras
    return letras


def _hoja_xml(columnas, filas):
    estilo_decimales = {0: 2, 1: 3, 2: 4, 3: 5, 4: 6}  # índice del estilo según decimales

    def celda_texto(ref, texto, estilo=0):
        texto = escape(str(texto))
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
               '<numFmts count="4"><numFmt numFmtId="164" formatCode="0.0"/>'
               '<numFmt numFmtId="165" formatCode="0.00"/><numFmt numFmtId="166" formatCode="0.000"/>'
               '<numFmt numFmtId="167" formatCode="0.0000"/></numFmts>'
               '<fonts count="2"><font><sz val="11"/><name val="Calibri"/></font>'
               '<font><b/><sz val="11"/><color rgb="FFFFFFFF"/><name val="Calibri"/></font></fonts>'
               '<fills count="3"><fill><patternFill patternType="none"/></fill>'
               '<fill><patternFill patternType="gray125"/></fill>'
               '<fill><patternFill patternType="solid"><fgColor rgb="FF203864"/><bgColor indexed="64"/>'
               '</patternFill></fill></fills>'
               '<borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders>'
               '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>'
               '<cellXfs count="7">'
               '<xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>'
               '<xf numFmtId="0" fontId="1" fillId="2" borderId="0" xfId="0" applyFont="1" applyFill="1" '
               'applyAlignment="1"><alignment vertical="center"/></xf>'
               '<xf numFmtId="1" fontId="0" fillId="0" borderId="0" xfId="0" applyNumberFormat="1"/>'
               '<xf numFmtId="164" fontId="0" fillId="0" borderId="0" xfId="0" applyNumberFormat="1"/>'
               '<xf numFmtId="165" fontId="0" fillId="0" borderId="0" xfId="0" applyNumberFormat="1"/>'
               '<xf numFmtId="166" fontId="0" fillId="0" borderId="0" xfId="0" applyNumberFormat="1"/>'
               '<xf numFmtId="167" fontId="0" fillId="0" borderId="0" xfId="0" applyNumberFormat="1"/>'
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
        self.root.title('Emillium')
        self.root.configure(background=FONDO)
        ancho_pantalla, alto_pantalla = root.winfo_screenwidth(), root.winfo_screenheight()
        self.root.geometry(f'{min(1500, int(ancho_pantalla * 0.9))}x{min(950, int(alto_pantalla * 0.85))}')
        self.root.minsize(900, 600)
        self.idioma = IDIOMA_POR_DEFECTO

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
        self.modo_brillo = 1     # 0 suave, 1 normal, 2 resaltar estrellas débiles

        # Mediciones
        self.mediciones = []
        self._contador = 0

        # Opciones que se conservan al cambiar de idioma
        self.var_angulo = tk.StringVar(value='0')
        self.var_angulo.trace_add('write', lambda *_: self._actualizar_sentido())
        self.var_busqueda = tk.StringVar(value='15')
        self.var_apertura = tk.StringVar(value='5')
        self.var_umbral = tk.StringVar(value='3')
        self.var_marcas = tk.BooleanVar(value=True)
        self.var_estrella = tk.StringVar(value='A')
        self.var_avanzar = tk.BooleanVar(value=True)
        self.avanzadas_visibles = False
        # Al cambiar una opción avanzada se vuelve a calcular la medida (la seleccionada o la última)
        self._recalculo_pendiente = None
        for var in (self.var_busqueda, self.var_apertura, self.var_umbral):
            var.trace_add('write', lambda *_: self._programar_recalculo())

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
        self.root.bind('<KeyPress>', self._tecla_secreta, add='+')
        for rueda in ('<MouseWheel>', '<Button-4>', '<Button-5>'):
            self.root.bind_all(rueda, self._rueda_panel, add='+')
        self.root.after(100, self._animar)

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

        # Tabla de mediciones (debajo de la imagen) y barra de vista (encima)
        self._panel_tabla()
        self._panel_vista()

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

        self.lbl_titulo_medida = ttk.Label(caja, text=self.t('ultima'), font=('Segoe UI', 10, 'bold'))
        self.lbl_titulo_medida.pack(anchor=tk.W, pady=(8, 0))
        self.lbl_x = ttk.Label(caja, text='X = -', style='Resultado.TLabel')
        self.lbl_x.pack(anchor=tk.W)
        self.lbl_y = ttk.Label(caja, text='Y = -', style='Resultado.TLabel')
        self.lbl_y.pack(anchor=tk.W)
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
                ('adv_umbral', self.var_umbral, 0, 50, 0.5, 'ayuda_umbral'))):
            etiqueta = ttk.Label(self.caja_avanzadas, text=self.t(clave))
            etiqueta.grid(row=fila, column=0, sticky=tk.W, pady=2)
            caja_num = ttk.Spinbox(self.caja_avanzadas, textvariable=var, from_=desde, to=hasta,
                                   increment=paso, width=6, command=self._programar_dibujo)
            caja_num.grid(row=fila, column=1, sticky=tk.E, padx=(8, 0), pady=2)
            Ayuda(etiqueta, self.t(ayuda))
            Ayuda(caja_num, self.t(ayuda))
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
        combo = ttk.Combobox(fila, state='readonly', width=24, values=self.t('modos_brillo'))
        combo.current(self.modo_brillo)
        combo.pack(side=tk.LEFT, padx=6)

        def elegir_brillo(_evento):
            self.modo_brillo = combo.current()
            self._programar_dibujo()
        combo.bind('<<ComboboxSelected>>', elegir_brillo)
        Ayuda(combo, self.t('ayuda_brillo'))
        ttk.Checkbutton(fila, text=self.t('marcas'), variable=self.var_marcas,
                        command=self._programar_dibujo).pack(side=tk.LEFT, padx=12)

        # Selector de idioma (a la derecha)
        idiomas = ttk.Combobox(fila, state='readonly', width=11, values=[nombre for _c, nombre in IDIOMAS])
        idiomas.current([c for c, _n in IDIOMAS].index(self.idioma))
        idiomas.pack(side=tk.RIGHT)
        idiomas.bind('<<ComboboxSelected>>',
                     lambda e: self.root.after_idle(self.cambiar_idioma, IDIOMAS[idiomas.current()][0]))
        ttk.Label(fila, text=self.t('idioma')).pack(side=tk.RIGHT, padx=(0, 4))

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
        ruta = filedialog.askopenfilename(
            title=self.t('abrir_titulo'),
            filetypes=[(self.t('tipo_fits'), '*.fits *.fit *.fts *.FITS *.FIT *.FTS'),
                       (self.t('tipo_todos'), '*.*')])
        if ruta:
            self.cargar_fits(ruta)

    def cargar_fits(self, ruta):
        self.root.config(cursor='watch')
        self.root.update_idletasks()
        try:
            original, cabecera, nota = leer_fits(ruta)
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

        self.ruta = ruta
        self.original = original
        self.nota = nota
        self.angulo = 0.0
        self.espejo = False
        self.mirall_h = self.mirall_v = False
        self.rotacion = 0.0
        self.var_angulo.set('0')
        # Número de foto: el orden en que se abrió por primera vez (si se vuelve a abrir, conserva el suyo)
        if ruta not in self.fotos:
            hora = cabecera.get('DATE-OBS', '')
            self.fotos[ruta] = {'n': len(self.fotos) + 1, 'archivo': os.path.basename(ruta),
                                'hora': str(hora).strip() if hora else ''}
        alto, ancho = original.shape
        self.fotos[ruta].update(ancho=ancho, alto=alto)
        self.var_estrella.set('A')  # en cada foto nueva se empieza por la estrella A
        # Niveles de gris: percentiles 1 % y 99,8 %
        self.nivel_bajo, self.nivel_alto = (float(v) for v in np.percentile(muestra, [1.0, 99.8]))
        if self.nivel_alto <= self.nivel_bajo:
            self.nivel_alto = self.nivel_bajo + 1.0

        self.imagen = original
        self._actualizar_orientacion()
        self.root.title('Emillium')
        self._mostrar_info_foto()
        if not self.boton_cerrar_foto.winfo_ismapped():
            self.boton_cerrar_foto.pack(anchor=tk.E, pady=(6, 0))
        self._mensaje('')
        self.ajustar_vista()

    def _mostrar_info_foto(self):
        if self.ruta is None:
            return
        foto = self.fotos[self.ruta]
        texto = self.t('info_foto', nombre=foto['archivo'], n=foto['n'], ancho=foto['ancho'], alto=foto['alto'])
        if self.nota:
            texto += '\n' + self.t(self.nota)
        self.lbl_archivo.config(text=texto, style='TLabel')

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
        if self.imagen is None:
            self.abrir_fits()
            return
        x, y = self.lienzo_a_imagen(evento.x, evento.y)
        self.medir_en(x, y)

    def medir_en(self, x, y):
        radio_busqueda = int(leer_numero(self.var_busqueda.get(), 15, 3, 200))
        radio_apertura = leer_numero(self.var_apertura.get(), 5.0, 1.0, 50.0)
        umbral = leer_numero(self.var_umbral.get(), 3.0, 0.0, 50.0)
        resultado = medir_estrella(self.imagen, x, y, radio_busqueda, radio_apertura, umbral)
        if resultado is None:
            self._mensaje(self.t('sin_estrella'))
            self.root.bell()
            return None

        letra = self.var_estrella.get() or 'A'
        # ¿Este centro ya es el de otra estrella medida (con otra letra) en esta foto? -> aviso y no se guarda
        for m in self.mediciones:
            if m['ruta'] != self.ruta or m['estrella'] == letra:
                continue
            ox, oy, _misma = self.posicion_en_vista(m)
            distancia = math.hypot(resultado['x'] - ox, resultado['y'] - oy)
            if distancia < 1.5:
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
            'archivo': os.path.basename(self.ruta),
            'ruta': self.ruta,
            'giro': self.angulo,
            'espejo': self.espejo,
            'orient': self.orientacion_actual(),
            'radio_apertura': radio_apertura,
            'fecha': datetime.now().isoformat(sep=' ', timespec='seconds'),
        })
        self.mediciones.append(medida)
        self._anadir_a_tabla(medida)
        self._actualizar_contador()
        self._mostrar_resultado(medida)
        avisos = []
        if anterior:
            avisos.append(self.t('sustituida', letra=letra))
        if medida['borde']:
            avisos.append(self.t('borde'))
        self._mensaje('\n'.join(avisos))
        # pasar a la siguiente letra
        if self.var_avanzar.get() and letra != LETRAS_ESTRELLA[-1]:
            self.var_estrella.set(LETRAS_ESTRELLA[LETRAS_ESTRELLA.index(letra) + 1])
        self._programar_dibujo()
        return medida

    def _mostrar_resultado(self, medida):
        self.lbl_titulo_medida.config(text=self.t('titulo_medida', letra=medida['estrella'], n=medida['foto_n']))
        self.lbl_x.config(text=f'X = {numero(medida["x"], 3)} px')
        self.lbl_y.config(text=f'Y = {numero(medida["y"], 3)} px')
        self.lbl_detalle.config(text=self.t('detalle', flujo=numero(medida['flujo'], 1),
                                            fwhm=numero(medida['fwhm'], 2), ps=numero(medida['pico_sigma'], 0)))

    def _valores_fila(self, m):
        return (m['foto_n'], m['archivo'], m['estrella'], self.describir_orientacion(m['orient']),
                numero(m['x'], 3), numero(m['y'], 3), numero(m['flujo'], 1))

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
        radio_busqueda = int(leer_numero(self.var_busqueda.get(), 15, 3, 200))
        radio_apertura = leer_numero(self.var_apertura.get(), 5.0, 1.0, 50.0)
        umbral = leer_numero(self.var_umbral.get(), 3.0, 0.0, 50.0)
        # se busca muy cerca de donde estaba, para no saltar a otra estrella
        resultado = medir_estrella(self.imagen, m['x'], m['y'], min(radio_busqueda, 4), radio_apertura, umbral)
        if resultado is None:
            self._mensaje(self.t('sin_estrella'))
            return
        m.update(resultado)
        m['radio_apertura'] = radio_apertura
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
        if self.ruta is None:
            return
        self.ruta = self.original = self.imagen = None
        self.nota = ''
        self.mirall_h = self.mirall_v = False
        self.rotacion = self.angulo = 0.0
        self.espejo = False
        self.var_angulo.set('0')
        self.lente, self.lente_theta = None, 0.0
        self.root.title('Emillium')
        self.lbl_archivo.config(text=self.t('ninguna_foto'), style='Gris.TLabel')
        self.boton_cerrar_foto.pack_forget()
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
        messagebox.showinfo(self.t('guardado_t'), self.t('guardado', n=len(self.mediciones), f=fotos, ruta=ruta))

    def tabla_por_estrellas(self):
        """Tabla "estrella por estrella": una fila por foto y, para cada estrella
        (A, B, C...), sus columnas X, Y y flujo. Devuelve (columnas, filas)."""
        letras = sorted({m['estrella'] for m in self.mediciones}, key=LETRAS_ESTRELLA.index)
        rutas = sorted({m['ruta'] for m in self.mediciones}, key=lambda r: self.fotos[r]['n'])
        con_hora = any(self.fotos[r]['hora'] for r in rutas)
        columnas = [(self.t('col_foto_n_largo'), 0, 11), (self.t('col_foto'), None, 28)]
        if con_hora:
            columnas.append((self.t('col_hora'), None, 24))
        columnas.append((self.t('col_orientacion'), None, 30))
        for letra in letras:
            columnas += [(self.t('x_col', l=letra), 3, 12), (self.t('y_col', l=letra), 3, 12),
                         (self.t('flujo_col', l=letra), 1, 12)]
        filas = []
        for ruta in rutas:
            de_esta_foto = {m['estrella']: m for m in self.mediciones if m['ruta'] == ruta}
            orientaciones = []
            for m in sorted(de_esta_foto.values(), key=lambda m: LETRAS_ESTRELLA.index(m['estrella'])):
                texto = self.describir_orientacion(m['orient'])
                if texto not in orientaciones:
                    orientaciones.append(texto)
            fila = [self.fotos[ruta]['n'], self.fotos[ruta]['archivo']]
            if con_hora:
                fila.append(self.fotos[ruta]['hora'])
            fila.append(' / '.join(orientaciones))
            for letra in letras:
                m = de_esta_foto.get(letra)
                fila += [m['x'], m['y'], m['flujo']] if m else [None, None, None]
            filas.append(fila)
        return columnas, filas

    def tabla_todas(self):
        """Lista simple: una fila por medida, ordenada por foto y estrella."""
        columnas = [(self.t('col_foto_n_largo'), 0, 11), (self.t('col_foto'), None, 28),
                    (self.t('col_estrella'), None, 9), (self.t('col_orientacion'), None, 30),
                    ('X (px)', 3, 12), ('Y (px)', 3, 12), (self.t('col_flujo'), 1, 12)]
        orden = sorted(self.mediciones, key=lambda m: (m['foto_n'], LETRAS_ESTRELLA.index(m['estrella'])))
        filas = [[m['foto_n'], m['archivo'], m['estrella'], self.describir_orientacion(m['orient']),
                  m['x'], m['y'], m['flujo']] for m in orden]
        return columnas, filas

    def escribir_xlsx(self, ruta):
        escribir_xlsx(ruta, [(self.t('hoja_estrellas'), *self.tabla_por_estrellas()),
                             (self.t('hoja_todas'), *self.tabla_todas())])

    def escribir_csv(self, ruta):
        # Separador ";" y coma decimal: así Excel configurado en catalán o español lo abre bien.
        columnas, filas = self.tabla_por_estrellas()
        with open(ruta, 'w', newline='', encoding='utf-8-sig') as archivo:
            escritor = csv.writer(archivo, delimiter=';')
            escritor.writerow([c[0] for c in columnas])
            for fila in filas:
                escritor.writerow(['' if valor is None else (valor if decimales is None else numero(valor, decimales))
                                   for valor, (_t, decimales, _a) in zip(fila, columnas)])


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
    if len(sys.argv) > 1 and os.path.isfile(sys.argv[1]):
        root.after(100, lambda: app.cargar_fits(sys.argv[1]))
    root.mainloop()


if __name__ == '__main__':
    main()

