"""
Emillium Expanded - cálculos científicos (solo numpy)
=====================================================

Este módulo no tiene ventanas: solo cálculos, para poder probarlos por separado
(ver la carpeta «proves»). Lo usa Emillium.py.

Convenio de píxeles (el mismo que en Emillium.py): imagen[fila, columna], con
x = columna e y = fila, y el centro del píxel (0, 0) en las coordenadas (0, 0).

Contenido
---------
  1. Fotometría de apertura (con área exacta de píxel, fondo en un anillo,
     saturación y contaminación por estrellas vecinas).
  2. Detección automática de fuentes y estadísticas de una imagen.
  3. Fotometría diferencial con varias estrellas de comparación ("ensemble"),
     término de color y detección de variabilidad.
  4. Geometría en el cielo: separación y ángulo de posición (estrellas dobles),
     proyección gnomónica, ajuste de una solución WCS y emparejado con catálogos.
  5. Objetos móviles: movimiento lineal, predicción y detección.
  6. Series de fotos: desplazamiento entre fotos y fotos defectuosas.
  7. Indicadores de calidad de cada medida.
"""

import math

import numpy as np

# =============================================================================
# Utilidades
# =============================================================================


def recorte_sigma(valores, n_sigma=3.0, iteraciones=5):
    """Mediana y desviación típica descartando valores atípicos (a más de n_sigma)."""
    valores = np.asarray(valores, dtype=np.float64)
    valores = valores[np.isfinite(valores)]
    if valores.size == 0:
        return float('nan'), float('nan')
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


def sigma_robusta(valores):
    """Ruido estimado con la mediana de las desviaciones (MAD), que no se deja
    engañar por las estrellas: 1,4826·MAD = σ para un ruido gaussiano."""
    valores = np.asarray(valores, dtype=np.float64)
    valores = valores[np.isfinite(valores)]
    if valores.size == 0:
        return float('nan')
    return float(1.4826 * np.median(np.abs(valores - np.median(valores))))


def nivel_saturacion(cabecera):
    """Valor a partir del cual la cámara está saturada: SATURATE (o SATLEVEL) de la
    cabecera; si no viene, el máximo que cabe en los datos enteros (por ejemplo
    65535 en una cámara de 16 bits). En imágenes de coma flotante no se sabe (None).
    (DATAMAX no sirve: muchas veces es simplemente el valor más alto de la imagen.)"""
    for clave in ('SATURATE', 'SATLEVEL'):
        try:
            valor = float(cabecera.get(clave))
            if valor > 0:
                return valor
        except (TypeError, ValueError):
            pass
    try:
        bitpix = int(cabecera.get('ZBITPIX', cabecera.get('BITPIX', 0)))
        bzero = float(cabecera.get('BZERO', 0) or 0)
        bscale = float(cabecera.get('BSCALE', 1) or 1)
    except (TypeError, ValueError):
        return None
    maximos = {8: 255, 16: 32767, 32: 2147483647}
    if bitpix not in maximos:
        return None
    return maximos[bitpix] * bscale + bzero


# =============================================================================
# 1. FOTOMETRÍA DE APERTURA
# =============================================================================
# La luz de la estrella se suma dentro de un círculo de radio r (la "apertura").
# Para que el resultado no dependa de dónde cae el centro respecto a la cuadrícula,
# cada píxel del borde cuenta solo la parte que queda dentro del círculo (se calcula
# dividiendo el píxel en sub×sub trocitos). El fondo del cielo se mide en un anillo
# (de r_int a r_ext) y se resta.
#
# Error del flujo (en las unidades de la imagen):
#   σF² = N·σ²·(1 + N/N_anillo)  +  F/ganancia
#   (ruido del cielo en la apertura, error de la estimación del fondo y, si se
#    conoce la ganancia, el ruido de la propia luz de la estrella)


def fraccion_en_circulo(xs, ys, cx, cy, r, sub=5):
    """Fracción de cada píxel (centros xs, ys) que cae dentro del círculo."""
    paso = (np.arange(sub) + 0.5) / sub - 0.5
    dx = xs[..., None, None] + paso[None, None, :] - cx
    dy = ys[..., None, None] + paso[None, :, None] - cy
    return ((dx * dx + dy * dy) <= r * r).mean(axis=(-1, -2))


def fotometria_apertura(imagen, x, y, r, r_int, r_ext, ganancia=None, saturacion=None, sub=5):
    """Fotometría de apertura en (x, y). Devuelve un diccionario o None si no se puede."""
    alto, ancho = imagen.shape
    if not (r > 0 and r_ext > r_int >= r):
        raise ValueError('radios incorrectos')
    m = int(math.ceil(r_ext)) + 1
    cx, cy = int(round(x)), int(round(y))
    x0, x1 = max(cx - m, 0), min(cx + m + 1, ancho)
    y0, y1 = max(cy - m, 0), min(cy + m + 1, alto)
    if x1 - x0 < 3 or y1 - y0 < 3:
        return None
    trozo = imagen[y0:y1, x0:x1].astype(np.float64)
    yy, xx = np.mgrid[y0:y1, x0:x1].astype(np.float64)
    d2 = (xx - x) ** 2 + (yy - y) ** 2
    anillo = trozo[(d2 >= r_int ** 2) & (d2 <= r_ext ** 2) & np.isfinite(trozo)]
    if anillo.size < 10:
        return None
    fondo, sigma = recorte_sigma(anillo)
    # Píxeles que tocan la apertura: completos o en parte
    cerca = d2 <= (r + 0.75) ** 2
    fraccion = np.zeros_like(trozo)
    fraccion[cerca] = fraccion_en_circulo(xx[cerca], yy[cerca], x, y, r, sub)
    validos = np.isfinite(trozo)
    area = float(fraccion.sum())
    area_valida = float(fraccion[validos].sum())
    if area <= 0:
        return None
    toca_borde = (x - r < -0.5 or y - r < -0.5 or x + r > ancho - 0.5 or y + r > alto - 0.5
                  or area_valida < 0.999 * area)
    senal = np.where(validos, trozo - fondo, 0.0)
    flujo = float((senal * fraccion).sum())
    varianza = area * sigma ** 2 * (1 + area / anillo.size)
    if ganancia and ganancia > 0:
        varianza += max(flujo, 0.0) / ganancia
    dentro = (fraccion > 0) & validos
    maximo = float(np.max(trozo[dentro])) if dentro.any() else float('nan')
    saturada = bool(saturacion is not None and maximo >= 0.97 * saturacion)
    return {
        'flujo': flujo,
        'sigma_flujo': math.sqrt(max(varianza, 0.0)),
        'fondo': fondo,
        'sigma_fondo': sigma,
        'area': area,
        'n_anillo': int(anillo.size),
        'maximo': maximo,
        'saturada': saturada,
        'borde': bool(toca_borde),
        'snr': flujo / math.sqrt(varianza) if varianza > 0 else float('inf'),
    }


def magnitud(flujo, sigma_flujo, exposicion=None, punto_cero=25.0):
    """m = punto_cero − 2,5·log10(flujo/t) y su error 1,0857·σF/F."""
    if flujo is None or not flujo > 0:
        return None, None
    por_segundo = flujo / exposicion if exposicion and exposicion > 0 else flujo
    return punto_cero - 2.5 * math.log10(por_segundo), 1.0857 * sigma_flujo / flujo


def contaminacion(imagen, x, y, r, fwhm, flujo, umbral_sigma=5.0, fondo=None, sigma=None):
    """Estrellas vecinas que pueden meter luz en la apertura de radio r.

    Busca picos claros (más brillantes que sus 8 vecinos y a más de umbral_sigma
    sobre el fondo) a menos de r + 2·FWHM del centro, sin contar la propia estrella.
    Para cada una estima qué fracción de su luz cae dentro de la apertura (suponiendo
    una estrella gaussiana de esa FWHM) y la compara con el flujo de la estrella.
    Devuelve (fracción de luz ajena estimada, lista de vecinas [(x, y, distancia)])."""
    fwhm = max(float(fwhm or 2.5), 1.0)
    alcance = r + 2.0 * fwhm
    m = int(math.ceil(alcance)) + 2
    alto, ancho = imagen.shape
    cx, cy = int(round(x)), int(round(y))
    x0, x1 = max(cx - m, 0), min(cx + m + 1, ancho)
    y0, y1 = max(cy - m, 0), min(cy + m + 1, alto)
    trozo = imagen[y0:y1, x0:x1].astype(np.float64)
    if trozo.shape[0] < 3 or trozo.shape[1] < 3:
        return 0.0, []
    if fondo is None or sigma is None:
        fondo, sigma = float(np.nanmedian(trozo)), sigma_robusta(trozo)
    if not sigma or not np.isfinite(sigma):
        return 0.0, []
    lleno = np.where(np.isfinite(trozo), trozo, fondo)
    rodeado = np.pad(lleno, 1, mode='constant', constant_values=-np.inf)
    h, w = lleno.shape
    es_pico = np.ones_like(lleno, dtype=bool)
    for j in range(3):
        for i in range(3):
            if (j, i) != (1, 1):
                es_pico &= lleno > rodeado[j:j + h, i:i + w]
    es_pico &= lleno - fondo > umbral_sigma * sigma
    filas, columnas = np.nonzero(es_pico)
    s = fwhm / 2.3548
    vecinas, ajena = [], 0.0
    for f, c in zip(filas, columnas):
        vx, vy = x0 + c, y0 + f
        d = math.hypot(vx - x, vy - y)
        if d < max(1.5, 0.6 * fwhm) or d > alcance:
            continue
        altura = lleno[f, c] - fondo
        flujo_vecina = 2 * math.pi * s * s * altura
        # fracción aproximada de la luz de la vecina dentro del círculo: gaussiana 2D
        # desplazada d (aproximación de Rice/Marcum con integración numérica sencilla)
        fr = _fraccion_gaussiana_en_circulo(d, r, s)
        ajena += flujo_vecina * fr
        vecinas.append((float(vx), float(vy), float(d)))
    return (ajena / flujo if flujo and flujo > 0 else 0.0), vecinas


def _fraccion_gaussiana_en_circulo(d, r, s, n=24):
    """Fracción de una gaussiana circular (σ = s) centrada a distancia d que cae
    dentro de un círculo de radio r (integración numérica en anillos)."""
    radios = (np.arange(n) + 0.5) / n * r
    angulos = (np.arange(48) + 0.5) / 48 * 2 * math.pi
    rr, aa = np.meshgrid(radios, angulos)
    px, py = rr * np.cos(aa) - d, rr * np.sin(aa)
    densidad = np.exp(-(px * px + py * py) / (2 * s * s)) / (2 * math.pi * s * s)
    return float((densidad * rr).sum() * (r / n) * (2 * math.pi / 48))


# =============================================================================
# 2. DETECCIÓN DE FUENTES Y ESTADÍSTICAS DE UNA IMAGEN
# =============================================================================


def fondo_y_ruido(imagen, muestra_max=400_000):
    """Fondo (mediana con recorte) y ruido robusto de toda la imagen, con una muestra."""
    paso = max(1, int(math.sqrt(imagen.size / muestra_max)))
    muestra = imagen[::paso, ::paso].astype(np.float64).ravel()
    muestra = muestra[np.isfinite(muestra)]
    if muestra.size == 0:
        return float('nan'), float('nan')
    fondo, _ = recorte_sigma(muestra)
    sigma = sigma_robusta(muestra)
    if not sigma > 0:
        # Imágenes de recuento de fotones con casi todos los píxeles a 0 (la MAD da 0):
        # desviación típica de los píxeles sin las estrellas (por debajo del percentil 99)
        sin_estrellas = muestra[muestra <= np.percentile(muestra, 99)]
        sigma = float(np.std(sin_estrellas)) if sin_estrellas.size else 0.0
    return fondo, sigma


def _cuadrado_bool(mascara, radio, erosion):
    """Erosión (erosion=True) o dilatación de una máscara con un cuadrado de lado 2·radio+1."""
    salida = mascara
    for eje in (0, 1):
        nueva = salida.copy()
        for d in range(1, radio + 1):
            if eje == 0:
                destino_a, origen_a, destino_b, origen_b = nueva[:-d], salida[d:], nueva[d:], salida[:-d]
            else:
                destino_a, origen_a, destino_b, origen_b = nueva[:, :-d], salida[:, d:], nueva[:, d:], salida[:, :-d]
            if erosion:
                destino_a &= origen_a
                destino_b &= origen_b
            else:
                destino_a |= origen_a
                destino_b |= origen_b
        salida = nueva
    return salida


def zona_sin_datos(imagen, muestra_max=400_000):
    """Máscara de los píxeles SIN DATOS: NaN, infinitos y zonas grandes con exactamente el
    mismo valor (lo que queda fuera del campo expuesto, que muchos telescopios rellenan con 0).
    Una zona cuenta si forma cuadrados de 5×5 píxeles idénticos: en el cielo de verdad el
    ruido hace que eso no pase nunca."""
    malo = ~np.isfinite(imagen)
    paso = max(1, int(math.sqrt(imagen.size / muestra_max)))
    muestra = imagen[::paso, ::paso].ravel()
    muestra = muestra[np.isfinite(muestra)]
    if muestra.size == 0 or min(imagen.shape) < 5:
        return malo
    valores, cuentas = np.unique(muestra, return_counts=True)
    for v in valores[np.argsort(-cuentas)][:3]:
        if (muestra == v).sum() < 0.02 * muestra.size:
            break
        zona = _cuadrado_bool(_cuadrado_bool(imagen == v, 2, True), 2, False)
        if (malo | zona).mean() < 0.97:     # si es casi toda la foto, el valor sí es el cielo
            malo |= zona
    return malo


def _malla_fondo(datos, malo, caja):
    """Fondo (mediana con recorte a 3σ) y ruido (rango intercuartílico) en cajas de
    caja×caja píxeles. Las cajas casi sin datos quedan en NaN."""
    alto, ancho = datos.shape
    s = max(1, caja // 32)                           # basta con ~32×32 píxeles por caja
    d = datos[::s, ::s].astype(np.float32, copy=True)
    d[malo[::s, ::s]] = np.nan
    cs = caja // s
    ny, nx = -(-alto // caja), -(-ancho // caja)
    relleno = np.full((ny * cs, nx * cs), np.nan, dtype=np.float32)
    relleno[:d.shape[0], :d.shape[1]] = d[:ny * cs, :nx * cs]
    cajas = np.sort(relleno.reshape(ny, cs, nx, cs).transpose(0, 2, 1, 3).reshape(ny, nx, cs * cs), axis=2)
    n = np.isfinite(cajas).sum(axis=2)
    indice = lambda k: np.take_along_axis(cajas, np.clip(k, 0, cs * cs - 1)[..., None], axis=2)[..., 0]
    mediana = indice((n - 1) // 2)
    # ruido: rango intercuartílico / 1,349 (como la desviación típica de una gaussiana)
    sigma = (indice(np.floor(0.75 * (n - 1)).astype(int)) - indice(np.floor(0.25 * (n - 1)).astype(int))) / 1.349
    # Una pasada de recorte a ±3σ (quita las estrellas) y mediana otra vez
    abajo = (cajas < (mediana - 3 * sigma)[..., None]).sum(axis=2)
    dentro = (cajas <= (mediana + 3 * sigma)[..., None]).sum(axis=2)
    mediana = np.where(dentro > abajo, indice(abajo + (dentro - abajo - 1) // 2), mediana)
    vacias = n < 0.3 * cs * cs
    mediana[vacias] = np.nan
    sigma[vacias] = np.nan
    return mediana.astype(np.float64), sigma.astype(np.float64)


def _completar_malla(malla):
    """Rellena las cajas sin datos con sus vecinas y suaviza con una mediana 3×3."""
    malla = malla.copy()
    if not np.isfinite(malla).any():
        return None
    for _ in range(max(malla.shape)):
        huecos = ~np.isfinite(malla)
        if not huecos.any():
            break
        amp = np.pad(malla, 1, mode='constant', constant_values=np.nan)
        vecinas = np.stack([amp[j:j + malla.shape[0], i:i + malla.shape[1]] for j in range(3) for i in range(3)])
        cuantas = np.isfinite(vecinas).sum(axis=0)
        media = np.where(cuantas > 0, np.nansum(vecinas, axis=0) / np.maximum(cuantas, 1), np.nan)
        malla[huecos] = media[huecos]
    amp = np.pad(malla, 1, mode='edge')
    return np.median(np.stack([amp[j:j + malla.shape[0], i:i + malla.shape[1]]
                               for j in range(3) for i in range(3)]), axis=0)


def _interpolar_malla(malla, caja, filas, columnas):
    """Valor de la malla (interpolación bilineal entre centros de caja) en unos píxeles.
    filas y columnas pueden ser listas de puntos o, con forma (alto,1) y (1,ancho), una rejilla."""
    ny, nx = malla.shape
    def eje(p, n):
        u = np.clip((np.asarray(p, dtype=np.float64) + 0.5) / caja - 0.5, 0, n - 1)
        i0 = np.minimum(np.floor(u).astype(int), max(n - 2, 0))
        return i0, np.minimum(i0 + 1, n - 1), (u - i0).astype(np.float32)
    j0, j1, wy = eje(filas, ny)
    i0, i1, wx = eje(columnas, nx)
    m = malla.astype(np.float32)
    return ((m[j0, i0] * (1 - wx) + m[j0, i1] * wx) * (1 - wy)
            + (m[j1, i0] * (1 - wx) + m[j1, i1] * wx) * wy)


def _suavizar3(datos):
    """Media 3×3 en dos pasadas (horizontal y vertical): menos memoria que sumar 9 copias."""
    amp = np.pad(datos, 1, mode='edge')
    horizontal = amp[:, :-2] + amp[:, 1:-1] + amp[:, 2:]
    del amp
    suave = horizontal[:-2] + horizontal[1:-1]
    suave += horizontal[2:]
    del horizontal
    suave /= np.float32(9.0)
    return suave


def _sin_halos(suave, filas, columnas, ruido):
    """Quita los grumos de los HALOS DE ESTRELLAS BRILLANTES (devuelve una máscara de las
    detecciones que se quedan). Los anillos y los rayos (aspas) que el telescopio pone
    alrededor de una estrella brillante tienen grumos que parecen estrellitas. Se reconocen
    porque se REPITEN: en el lado opuesto de la estrella brillante, a la misma distancia,
    hay otro grumo parecido. Una estrella de verdad que está al lado no tiene esa copia."""
    alto, ancho = suave.shape
    n = filas.size
    valor = suave[filas, columnas].astype(np.float64)
    ruido = np.broadcast_to(np.asarray(ruido, dtype=np.float64), (n,))   # ruido del suavizado
    queda = np.ones(n, dtype=bool)

    # Alrededor de las 50 más brillantes, grumos con "copias"
    angulos = np.linspace(-math.pi, math.pi, 256, endpoint=False)
    for j in np.argsort(-valor)[:50]:
        if not queda[j]:
            continue
        dy, dx = (filas - filas[j]).astype(np.float64), (columnas - columnas[j]).astype(np.float64)
        d = np.hypot(dx, dy)
        cand = np.nonzero(queda & (d >= 4) & (d <= 150) & (valor * 5 <= valor[j]))[0]
        if cand.size == 0:
            continue
        r = d[cand][:, None]
        ff = np.clip(np.rint(filas[j] + r * np.sin(angulos)[None]).astype(int), 0, alto - 1)
        cc = np.clip(np.rint(columnas[j] + r * np.cos(angulos)[None]).astype(int), 0, ancho - 1)
        circulo = suave[ff, cc].astype(np.float64)
        mediana = np.median(circulo, axis=1)
        exceso = valor[cand] - mediana                        # lo que el grumo sobresale
        # ángulo de cada punto del círculo respecto a la dirección del grumo
        delta = (angulos[None] - np.arctan2(dy[cand], dx[cand])[:, None] + math.pi) % (2 * math.pi) - math.pi
        # Las figuras del telescopio son simétricas: hay una copia JUSTO en el lado opuesto
        # (la otra punta del aspa, el otro lado del anillo), a menos de 2,5 px de ese punto.
        # Solo cuenta DENTRO del halo (donde la luz de la estrella brillante aún se nota), y
        # la copia tiene que parecerse (ni la mitad ni el doble): en un campo con muchas
        # estrellas, enfrente puede haber otra estrella de verdad por casualidad.
        cerca_opuesto = np.abs(delta) >= math.pi - np.maximum(2.5 / r, math.radians(2))
        opuesto = np.where(cerca_opuesto, circulo, -np.inf).max(axis=1) - mediana
        en_halo = mediana > ruido[cand]
        # ... o unido a la estrella brillante por una línea de luz (un aspa), mirando solo el
        # tramo a más de 4 px de las dos (para no confundirlo con la luz de ellas mismas)
        t = np.linspace(0.0, 1.0, 16)[None]
        largo = np.maximum(r - 8.0, 0.0)
        fr = np.rint(filas[j] + dy[cand][:, None] * ((4.0 + largo * t) / r)).astype(int)
        cr = np.rint(columnas[j] + dx[cand][:, None] * ((4.0 + largo * t) / r)).astype(int)
        en_halo |= (d[cand] >= 12) & (suave[fr, cr].min(axis=1) > 2 * ruido[cand])
        # ... o la estrella es tan brillante (50 veces o más) que sus figuras llegan lejos
        en_halo |= valor[cand] * 50 <= valor[j]
        copia = (opuesto >= 0.5 * exceso) & (opuesto <= 2.0 * exceso)
        queda[cand[en_halo & copia]] = False
    return queda


def detectar_fuentes(imagen, umbral_sigma=5.0, max_fuentes=2000, borde=4, saturacion=None):
    """Busca estrellas en toda la imagen.

    1) Zonas sin datos (NaN, o fuera del campo expuesto rellenas con un valor fijo).
    2) Fondo y ruido LOCALES: la foto se divide en cajas y en cada una se mide el cielo,
       así un degradado, una nebulosa o la luz de la Luna no crean estrellas falsas.
    3) Suavizado 3×3 (para que el ruido no cree picos falsos).  4) Picos: máximos
    locales 5×5 que superan el fondo en umbral_sigma veces el ruido (del suavizado).
    5) Se descartan los objetos extensos o planos (rótulos, bordes, bloques brillantes):
    una estrella cae alrededor de su pico; un objeto así sigue brillante a 6, 10 y 15 px.
    6) Para cada pico, centro de masas en 5×5, flujo, FWHM y elongación por momentos.

    Devuelve una lista de diccionarios ordenada de más a menos brillante."""
    # Se trabaja en precisión simple y sin copias de más: con una foto de 24 MP cada copia
    # en doble precisión son 192 MB.
    lleno = np.asarray(imagen, dtype=np.float32)
    if lleno.ndim != 2:
        return []
    alto, ancho = lleno.shape
    borde = max(borde, 2)
    if alto <= 2 * borde or ancho <= 2 * borde:
        return []
    malo = zona_sin_datos(lleno)
    hay_malo = bool(malo.any())
    if hay_malo:
        validos = lleno[~malo]
        if validos.size < 25:
            return []
        fondo_global, sigma = fondo_y_ruido(validos[None, :])
        lleno = np.where(malo, np.float32(fondo_global), lleno)
    else:
        fondo_global, sigma = fondo_y_ruido(lleno)
    if not np.isfinite(sigma) or sigma <= 0:
        return []

    # Fondo local (malla de cajas). En fotos pequeñas, un solo valor para toda la foto.
    caja = 64 if min(alto, ancho) >= 512 else 32 if min(alto, ancho) >= 128 else 0
    malla = None
    if caja:
        malla, ruido_bruto = _malla_fondo(lleno, malo, caja)
        malla = _completar_malla(malla)
        # Con un degradado, el ruido de toda la foto junta sale enorme: vale más el de las cajas
        ruido_bruto = ruido_bruto[np.isfinite(ruido_bruto) & (ruido_bruto > 0)]
        if ruido_bruto.size:
            sigma = min(sigma, float(np.median(ruido_bruto)))
    suave = _suavizar3(lleno)
    if malla is not None:
        # Interpolación separable: primero a lo ancho (una fila por fila de cajas) y luego,
        # por bandas de filas, entre las dos filas de cajas más cercanas.
        por_filas = np.stack([_interpolar_malla(malla[j:j + 1], caja, np.zeros(1, int), np.arange(ancho))
                              for j in range(malla.shape[0])])
        u = np.clip((np.arange(alto) + 0.5) / caja - 0.5, 0, malla.shape[0] - 1)
        j0 = np.minimum(np.floor(u).astype(int), max(malla.shape[0] - 2, 0))
        j1 = np.minimum(j0 + 1, malla.shape[0] - 1)
        wy = (u - j0).astype(np.float32)[:, None]
        for a in range(0, alto, 1024):                # por bandas: menos memoria
            b = slice(a, a + 1024)
            suave[b] -= por_filas[j0[b]] * (1 - wy[b]) + por_filas[j1[b]] * wy[b]
    else:
        suave -= np.float32(fondo_global)
    # Ruido de la imagen suavizada, medido directamente (si los píxeles fueran independientes
    # bajaría 3 veces, pero en imágenes remuestreadas los vecinos se parecen y baja menos),
    # también por cajas: en las zonas brillantes (nebulosas) el ruido es mayor.
    _f, sigma_suave = fondo_y_ruido(suave[~malo][None, :] if hay_malo else suave)
    sigma_suave = max(sigma_suave, sigma / 3.0) if sigma_suave > 0 else sigma / 3.0
    malla_ruido = None
    if caja:
        _m, ruido_cajas = _malla_fondo(suave, malo, caja)
        ruido_cajas[~(ruido_cajas > 0)] = np.nan
        malla_ruido = _completar_malla(ruido_cajas)
        if malla_ruido is not None:
            malla_ruido = np.maximum(malla_ruido, 0.5 * sigma_suave)
    # Candidatos: píxeles que destacan (lejos del borde y de las zonas sin datos). Solo con
    # ellos se mira si son máximos locales en 5×5: así no hay que comparar la imagen 24 veces.
    minimo = umbral_sigma * (float(malla_ruido.min()) if malla_ruido is not None else sigma_suave)
    interior = suave[borde:-borde, borde:-borde]
    filas, columnas = np.nonzero(interior > np.float32(minimo))
    if filas.size == 0:
        return []
    filas += borde
    columnas += borde
    if malla_ruido is not None:
        ok = suave[filas, columnas] > umbral_sigma * _interpolar_malla(malla_ruido, caja, filas, columnas)
        filas, columnas = filas[ok], columnas[ok]
    if hay_malo:
        cerca = _cuadrado_bool(malo, 6, False)     # el borde del campo suele ser raro
        ok = ~cerca[filas, columnas]
        filas, columnas = filas[ok], columnas[ok]
        del cerca
    valor = suave[filas, columnas]
    es_pico = np.ones(filas.size, dtype=bool)
    for dj in range(-2, 3):
        for di in range(-2, 3):
            if (dj, di) == (0, 0):
                continue
            vecino = suave[filas + dj, columnas + di]
            # empates: solo cuenta el primero (para no duplicar un pico plano)
            es_pico &= (valor > vecino) | ((valor == vecino) & ((dj, di) > (0, 0)))
    filas, columnas = filas[es_pico], columnas[es_pico]
    if filas.size > 20 * max_fuentes:   # imagen llena de ruido: nos quedamos con los más altos
        orden = np.argsort(-suave[filas, columnas])[:20 * max_fuentes]
        filas, columnas = filas[orden], columnas[orden]
    # Objetos extensos o planos, en 3 círculos (radio 6, 10 y 15 px) de 16 puntos alrededor:
    #  - ¿cuántos puntos siguen por encima de la mitad del pico? En una estrella, casi ninguno
    #    (aunque esté saturada, a 15 px ya ha bajado); en un rótulo, un bloque o un borde,
    #    muchos en los tres círculos.
    #  - una estrella es COMPACTA: su pico destaca sobre lo que tiene alrededor (la mediana del
    #    círculo) al menos en uno de los círculos. Una nebulosa, un arco o un halo brillante
    #    supera el fondo de la caja, pero no destaca sobre sus propios alrededores.
    if filas.size:
        angulos = np.arange(16) * (2 * math.pi / 16)
        pico_s = suave[filas, columnas]
        ruido_punto = (umbral_sigma * _interpolar_malla(malla_ruido, caja, filas, columnas)
                       if malla_ruido is not None else np.float32(umbral_sigma * sigma_suave))
        extenso = np.ones(filas.size, dtype=bool)
        compacto = np.zeros(filas.size, dtype=bool)
        linea = np.ones(filas.size, dtype=bool)
        for r in (6, 10, 15):
            ff = np.clip(filas[:, None] + np.round(r * np.sin(angulos)).astype(int)[None], 0, alto - 1)
            cc = np.clip(columnas[:, None] + np.round(r * np.cos(angulos)).astype(int)[None], 0, ancho - 1)
            circulo = suave[ff, cc]
            alto_circ = circulo > 0.5 * pico_s[:, None]
            extenso &= alto_circ.sum(axis=1) >= 4
            compacto |= pico_s - np.median(circulo, axis=1) > ruido_punto
            # - una LÍNEA recta (borde del sensor, rayo, satélite) pasa por el pico y sale por
            #   dos puntos opuestos del círculo, en los tres radios
            if r < 15:
                linea &= (alto_circ[:, :8] & alto_circ[:, 8:]).any(axis=1)
        ok = compacto & ~extenso & ~linea
        filas, columnas = filas[ok], columnas[ok]
        if filas.size > 1:
            ruido_punto = np.broadcast_to(ruido_punto, ok.shape)[ok] / umbral_sigma
            ok = _sin_halos(suave, filas, columnas, ruido_punto)
            filas, columnas = filas[ok], columnas[ok]
    del suave
    if filas.size == 0:
        return []
    if malla is not None:
        fondo_punto = _interpolar_malla(malla, caja, filas, columnas).astype(np.float64)
    else:
        fondo_punto = np.full(filas.size, fondo_global)

    def ventana(k):
        off = np.arange(-k, k + 1)
        oy, ox = np.meshgrid(off, off, indexing='ij')
        ff = np.clip(filas[:, None, None] + oy[None], 0, alto - 1)
        cc = np.clip(columnas[:, None, None] + ox[None], 0, ancho - 1)
        return ox, oy, np.clip(lleno[ff, cc].astype(np.float64) - fondo_punto[:, None, None], 0, None)
    # Posición: centro de masas en una ventana 5×5 (pequeña, para que no entren estrellas
    # vecinas en los campos densos).
    ox, oy, luz = ventana(2)
    suma = luz.sum(axis=(1, 2))
    ok = suma > 0
    filas, columnas, fondo_punto = filas[ok], columnas[ok], fondo_punto[ok]
    mx = (luz[ok] * ox[None]).sum(axis=(1, 2)) / suma[ok]
    my = (luz[ok] * oy[None]).sum(axis=(1, 2)) / suma[ok]
    # Centro plano (estrella saturada): 3 píxeles o más exactamente al máximo, o el pico por
    # encima del nivel de saturación
    maximo = luz[ok].max(axis=(1, 2))
    plana = (luz[ok] >= maximo[:, None, None] - 1e-6 * np.abs(maximo)[:, None, None]).sum(axis=(1, 2)) >= 3
    if saturacion is not None:
        plana |= lleno[filas, columnas] >= 0.97 * saturacion
    # Forma (FWHM y elongación) y flujo: ventana 11×11. Una más pequeña recortaría las
    # estrellas grandes (desenfocadas) y daría una FWHM demasiado pequeña.
    ox, oy, luz = ventana(5)
    total = np.maximum(luz.sum(axis=(1, 2)), 1e-30)
    # Estrellas saturadas (centro plano): el pico puede quedar en cualquier punto de la
    # meseta y el 5×5 no sabe hacia dónde está el centro. El centro es el de la meseta:
    # el centro de los píxeles casi al máximo en una ventana 21×21 (recentrada si hace falta).
    if plana.any():
        indices = np.nonzero(plana)[0]
        for _ in range(4):
            sub_f, sub_c = filas[indices], columnas[indices]
            off = np.arange(-10, 11)
            oyy, oxx = np.meshgrid(off, off, indexing='ij')
            grande = lleno[np.clip(sub_f[:, None, None] + oyy[None], 0, alto - 1),
                           np.clip(sub_c[:, None, None] + oxx[None], 0, ancho - 1)]
            meseta = grande >= grande[:, 10, 10][:, None, None] - 0.02 * np.abs(grande[:, 10, 10] - fondo_punto[indices])[:, None, None]
            n = np.maximum(meseta.sum(axis=(1, 2)), 1)
            cx, cy = (meseta * oxx[None]).sum(axis=(1, 2)) / n, (meseta * oyy[None]).sum(axis=(1, 2)) / n
            paso_x, paso_y = np.round(cx).astype(int), np.round(cy).astype(int)
            filas[indices] = np.clip(sub_f + paso_y, 0, alto - 1)
            columnas[indices] = np.clip(sub_c + paso_x, 0, ancho - 1)
            mx[indices], my[indices] = cx - paso_x, cy - paso_y
            if not ((paso_x != 0) | (paso_y != 0)).any():
                break
        ox, oy, luz = ventana(5)
        total = np.maximum(luz.sum(axis=(1, 2)), 1e-30)
    dx, dy = ox[None] - mx[:, None, None], oy[None] - my[:, None, None]
    mxx = (luz * dx * dx).sum(axis=(1, 2)) / total
    myy = (luz * dy * dy).sum(axis=(1, 2)) / total
    mxy = (luz * dx * dy).sum(axis=(1, 2)) / total
    traza, det = mxx + myy, np.sqrt(np.maximum((mxx - myy) ** 2 + 4 * mxy ** 2, 0))
    l1, l2 = (traza + det) / 2, np.maximum((traza - det) / 2, 1e-6)
    fwhm = 2.3548 * np.sqrt(np.maximum(traza / 2, 0))
    elongacion = np.sqrt(l1 / l2)
    pico = lleno[filas, columnas] - fondo_punto
    fuentes = []
    for i in np.argsort(-total)[:max_fuentes]:
        fuentes.append({
            'x': float(columnas[i] + mx[i]), 'y': float(filas[i] + my[i]),
            'flujo': float(total[i]), 'pico': float(pico[i]), 'fwhm': float(fwhm[i]),
            'elongacion': float(elongacion[i]),
            'saturada': bool(saturacion is not None and lleno[filas[i], columnas[i]] >= 0.97 * saturacion),
        })
    return fuentes


def estadisticas_imagen(imagen, saturacion=None):
    """Resumen de la calidad de una foto: fondo, ruido, número de estrellas, FWHM y
    elongación medianas (de las estrellas no saturadas más brillantes)."""
    fondo, ruido = fondo_y_ruido(imagen)
    fuentes = detectar_fuentes(imagen, 5.0, 1500, saturacion=saturacion)
    buenas = [f for f in fuentes if not f['saturada'] and 0.8 < f['fwhm'] < 30][:60]
    return {
        'fondo': fondo, 'ruido': ruido, 'n_fuentes': len(fuentes),
        'fwhm': float(np.median([f['fwhm'] for f in buenas])) if buenas else float('nan'),
        'elongacion': float(np.median([f['elongacion'] for f in buenas])) if buenas else float('nan'),
        'fuentes': fuentes,
    }


def fotos_defectuosas(estadisticas):
    """Marca las fotos raras de una serie comparándolas con la mediana de todas.

    estadisticas: lista de diccionarios de estadisticas_imagen (uno por foto).
    Devuelve una lista paralela de listas de motivos ('desenfocada', 'movida',
    'pocas_estrellas', 'fondo_alto')."""
    if not estadisticas:
        return []
    def mediana(clave):
        valores = [e[clave] for e in estadisticas if np.isfinite(e[clave])]
        return float(np.median(valores)) if valores else float('nan')
    fwhm_m, n_m = mediana('fwhm'), mediana('n_fuentes')
    fondos = np.array([e['fondo'] for e in estadisticas], dtype=np.float64)
    fondo_m = float(np.nanmedian(fondos))
    ruido_m = mediana('ruido')
    mad_fondo = sigma_robusta(fondos) if len(estadisticas) >= 3 else 0.0
    motivos = []
    for e in estadisticas:
        lista = []
        if np.isfinite(e['fwhm']) and np.isfinite(fwhm_m) and e['fwhm'] > 1.5 * fwhm_m:
            lista.append('desenfocada')
        if np.isfinite(e['elongacion']) and e['elongacion'] > 1.6:
            lista.append('movida')
        if n_m > 0 and e['n_fuentes'] < 0.5 * n_m:
            lista.append('pocas_estrellas')
        limite = max(5 * mad_fondo, 3 * (ruido_m if np.isfinite(ruido_m) else 0))
        if np.isfinite(e['fondo']) and limite > 0 and e['fondo'] - fondo_m > limite:
            lista.append('fondo_alto')
        motivos.append(lista)
    return motivos


def desplazamiento_entre(fuentes_a, fuentes_b, n=40, tolerancia=2.0, max_desplazamiento=None):
    """Desplazamiento (dx, dy) que lleva las estrellas de la foto A a la foto B
    (solo traslación). Se prueban todas las parejas entre las n más brillantes de
    cada foto y se queda el desplazamiento que más parejas comparten ("votación");
    luego se afina con la media de las parejas que coinciden.
    Devuelve (dx, dy, número de estrellas emparejadas) o None."""
    a = np.array([(f['x'], f['y']) for f in fuentes_a[:n]], dtype=np.float64)
    b = np.array([(f['x'], f['y']) for f in fuentes_b[:n]], dtype=np.float64)
    if len(a) < 3 or len(b) < 3:
        return None
    dif = (b[None, :, :] - a[:, None, :]).reshape(-1, 2)
    if max_desplazamiento is not None:
        dif = dif[np.hypot(dif[:, 0], dif[:, 1]) <= max_desplazamiento]
        if len(dif) == 0:
            return None
    # votación: para cada candidato, cuántos otros caen a menos de "tolerancia"
    mejor, votos_max = None, 0
    for d in dif:
        votos = int((np.hypot(dif[:, 0] - d[0], dif[:, 1] - d[1]) < tolerancia).sum())
        if votos > votos_max:
            mejor, votos_max = d, votos
    if mejor is None or votos_max < 3:
        return None
    cerca = dif[np.hypot(dif[:, 0] - mejor[0], dif[:, 1] - mejor[1]) < tolerancia]
    return float(cerca[:, 0].mean()), float(cerca[:, 1].mean()), int(len(cerca))


# =============================================================================
# 3. FOTOMETRÍA DIFERENCIAL CON VARIAS ESTRELLAS DE COMPARACIÓN
# =============================================================================
# En cada foto, cada estrella de comparación c con magnitud de catálogo M_c da una
# estimación del "punto cero" de esa foto:   ZP_c = M_c − m_inst,c − T·IC_c
# (T = coeficiente de color de la cámara y el filtro, IC = índice de color; si no se
# conocen, T = 0). Se hace la media ponderada (por 1/σ²) de todas y:
#   M_objetivo = m_inst,objetivo + ZP + T·IC_objetivo
# El error junta el del objetivo y el de la media, y la dispersión entre las
# comparaciones avisa si alguna es variable o está mal.
#
# Si no se conoce la magnitud de ninguna comparación, se hace fotometría RELATIVA:
# ZP = −(media de las magnitudes instrumentales de las comparaciones), así que
# el resultado es la diferencia respecto a la media de las comparaciones (útil para
# curvas de luz: la forma de la curva es correcta, aunque no su nivel absoluto).
# La extinción atmosférica de primer orden es la misma para todas las estrellas de una
# misma foto, así que se cancela en la fotometría diferencial.


def punto_cero(comparaciones, color_t=0.0):
    """comparaciones: lista de (m_inst, σ_inst, M_catalogo o None, IC o None) de UNA foto.
    Devuelve (ZP, σ_ZP, dispersión, número usado, relativa) o None."""
    usadas = [c for c in comparaciones if c[0] is not None and np.isfinite(c[0])]
    if not usadas:
        return None
    relativa = any(c[2] is None for c in usadas)
    if relativa:
        # Sin magnitudes de catálogo: el punto cero es menos la media simple de las
        # comparaciones. Que unas sean más brillantes que otras no es un error, así que
        # su diferencia no cuenta; solo cuentan sus errores de medida.
        m = np.array([c[0] for c in usadas], dtype=np.float64)
        s = np.array([max(c[1] or 0.0, 0.001) for c in usadas], dtype=np.float64)
        return -float(m.mean()), float(math.sqrt((s ** 2).sum()) / len(m)), 0.0, len(m), True
    estimaciones, pesos = [], []
    for m_inst, sigma, catalogo, ic in usadas:
        termino_color = color_t * ic if (color_t and ic is not None) else 0.0
        valor = (catalogo - m_inst - termino_color) if not relativa else -m_inst
        estimaciones.append(valor)
        pesos.append(1.0 / max(sigma or 0.01, 0.001) ** 2)
    estimaciones, pesos = np.array(estimaciones), np.array(pesos)
    zp = float((estimaciones * pesos).sum() / pesos.sum())
    sigma_zp = float(1.0 / math.sqrt(pesos.sum()))
    dispersion = float(np.std(estimaciones, ddof=1)) if len(estimaciones) > 1 else 0.0
    if len(estimaciones) > 1:  # error de la media realista: el mayor entre el teórico y el observado
        sigma_zp = max(sigma_zp, dispersion / math.sqrt(len(estimaciones)))
    return zp, sigma_zp, dispersion, len(estimaciones), relativa


def magnitud_diferencial(m_inst, sigma_inst, zp_info, color_t=0.0, ic=None):
    """Magnitud del objetivo a partir del punto cero de su foto."""
    if m_inst is None or zp_info is None:
        return None, None
    zp, sigma_zp = zp_info[0], zp_info[1]
    termino_color = color_t * ic if (color_t and ic is not None) else 0.0
    return m_inst + zp + termino_color, math.hypot(sigma_inst or 0.0, sigma_zp)


def analisis_variabilidad(magnitudes, errores):
    """¿Es variable una estrella? Compara la dispersión de sus magnitudes con lo que
    esperaríamos solo por los errores (χ² reducido de una recta horizontal).
    Devuelve n, media, desviación, amplitud, χ²_red y si parece variable."""
    m = np.asarray([v for v in magnitudes if v is not None], dtype=np.float64)
    e = np.asarray([max(err or 0.0, 1e-4) for v, err in zip(magnitudes, errores) if v is not None])
    n = len(m)
    if n < 3:
        return {'n': n, 'media': float(m.mean()) if n else None, 'desviacion': None, 'amplitud': None,
                'chi2_red': None, 'variable': False}
    pesos = 1 / e ** 2
    media = float((m * pesos).sum() / pesos.sum())
    chi2_red = float((((m - media) / e) ** 2).sum() / (n - 1))
    amplitud = float(m.max() - m.min())
    variable = bool(n >= 5 and chi2_red > 5.0 and amplitud > 5 * float(np.median(e)))
    return {'n': n, 'media': media, 'desviacion': float(np.std(m, ddof=1)), 'amplitud': amplitud,
            'chi2_red': chi2_red, 'variable': variable}


def elegir_comparaciones(series, excluir=(), maximo=4):
    """Elige automáticamente las mejores estrellas de comparación.

    series: {letra: {foto: (m_inst, σ)}}. Se quedan las estrellas medidas en todas
    las fotos y se ordenan por lo constantes que son respecto a la media de las demás
    (las variables y las mal medidas quedan las últimas). Devuelve una lista de letras."""
    fotos = set.intersection(*[set(s) for s in series.values()]) if series else set()
    candidatas = [l for l, s in series.items() if l not in excluir and set(s) >= fotos and fotos]
    if len(candidatas) < 2:
        return candidatas[:maximo]
    fotos = sorted(fotos)
    matriz = np.array([[series[l][f][0] for f in fotos] for l in candidatas], dtype=np.float64)
    errores = np.array([[series[l][f][1] or 0.01 for f in fotos] for l in candidatas], dtype=np.float64)
    puntuacion = []
    for i, letra in enumerate(candidatas):
        otras = np.delete(matriz, i, axis=0).mean(axis=0)
        dif = matriz[i] - otras
        dispersion = float(np.std(dif, ddof=1)) if len(fotos) > 1 else 0.0
        puntuacion.append((dispersion + float(np.median(errores[i])), letra))
    return [l for _p, l in sorted(puntuacion)[:maximo]]


# =============================================================================
# 4. GEOMETRÍA EN EL CIELO, WCS Y CATÁLOGOS
# =============================================================================


def separacion_y_angulo(ar1, dec1, ar2, dec2):
    """Separación (segundos de arco) y ángulo de posición (grados, desde el norte
    hacia el este) de la estrella 2 vista desde la 1. Fórmula de Vincenty (exacta
    también para separaciones muy pequeñas)."""
    a1, d1, a2, d2 = (math.radians(v) for v in (ar1, dec1, ar2, dec2))
    da = a2 - a1
    num = math.hypot(math.cos(d2) * math.sin(da),
                     math.cos(d1) * math.sin(d2) - math.sin(d1) * math.cos(d2) * math.cos(da))
    den = math.sin(d1) * math.sin(d2) + math.cos(d1) * math.cos(d2) * math.cos(da)
    separacion = math.degrees(math.atan2(num, den)) * 3600
    angulo = math.degrees(math.atan2(math.sin(da) * math.cos(d2),
                                     math.cos(d1) * math.sin(d2) - math.sin(d1) * math.cos(d2) * math.cos(da)))
    return separacion, angulo % 360.0


def media_y_error(valores, circular=False):
    """Media, desviación típica y error de la media. Con circular=True (ángulos en
    grados) se promedia como vectores, para que 359° y 1° den 0° y no 180°."""
    v = np.asarray([x for x in valores if x is not None and np.isfinite(x)], dtype=np.float64)
    if v.size == 0:
        return None, None, None
    if circular:
        rad = np.radians(v)
        media = math.degrees(math.atan2(np.sin(rad).mean(), np.cos(rad).mean())) % 360.0
        desv = np.degrees(np.angle(np.exp(1j * (rad - math.radians(media)))))
        s = float(np.sqrt((desv ** 2).sum() / (v.size - 1))) if v.size > 1 else 0.0
    else:
        media = float(v.mean())
        s = float(v.std(ddof=1)) if v.size > 1 else 0.0
    return media, s, s / math.sqrt(v.size)


def proyectar(ar, dec, ar0, dec0):
    """Proyección gnomónica (TAN): (AR, Dec) -> coordenadas estándar (ξ, η) en grados
    sobre el plano tangente en (ar0, dec0). Acepta arrays."""
    a, d = np.radians(ar), np.radians(dec)
    a0, d0 = math.radians(ar0), math.radians(dec0)
    cosc = np.sin(d0) * np.sin(d) + np.cos(d0) * np.cos(d) * np.cos(a - a0)
    xi = np.cos(d) * np.sin(a - a0) / cosc
    eta = (np.cos(d0) * np.sin(d) - np.sin(d0) * np.cos(d) * np.cos(a - a0)) / cosc
    return np.degrees(xi), np.degrees(eta)


def desproyectar(xi, eta, ar0, dec0):
    """Inversa de proyectar: (ξ, η) en grados -> (AR, Dec) en grados."""
    x, y = np.radians(xi), np.radians(eta)
    a0, d0 = math.radians(ar0), math.radians(dec0)
    den = np.cos(d0) - y * np.sin(d0)
    ar = a0 + np.arctan2(x, den)
    dec = np.arctan2(np.sin(d0) + y * np.cos(d0), np.hypot(x, den))
    return np.degrees(ar) % 360.0, np.degrees(dec)


def ajustar_wcs(x_fits, y_fits, ar, dec, crpix=None, iteraciones=6, recorte=3.0):
    """Ajusta una solución astrométrica TAN lineal (CRVAL, CRPIX y matriz CD) a partir
    de estrellas con posición en la foto (píxeles FITS, empiezan en 1) y en el cielo.

    Método: se fija CRPIX (el centro de las estrellas o el que se dé), se proyectan
    las estrellas del catálogo en el plano tangente y se ajusta por mínimos cuadrados
        ξ = CD11·u + CD12·v + ξ0      η = CD21·u + CD22·v + η0     (u, v = píxel − CRPIX)
    El punto (ξ0, η0) dice cuánto hay que mover el punto de tangencia (CRVAL); se
    repite hasta que no se mueve. Se descartan las parejas con residuo > recorte·σ
    (estrellas mal emparejadas).

    Devuelve un diccionario con la cabecera WCS, los residuos en segundos de arco y
    qué estrellas se han usado, o None si hay menos de 3 estrellas."""
    x_fits, y_fits, ar, dec = (np.asarray(v, dtype=np.float64) for v in (x_fits, y_fits, ar, dec))
    n = len(x_fits)
    if n < 3:
        return None
    if crpix is None:
        crpix = (float(np.mean(x_fits)), float(np.mean(y_fits)))
    usar = np.ones(n, dtype=bool)
    # punto de tangencia inicial: el centro de las estrellas en el cielo
    vx = np.cos(np.radians(dec)) * np.cos(np.radians(ar))
    vy = np.cos(np.radians(dec)) * np.sin(np.radians(ar))
    vz = np.sin(np.radians(dec))
    ar0 = math.degrees(math.atan2(vy.mean(), vx.mean())) % 360.0
    dec0 = math.degrees(math.atan2(vz.mean(), math.hypot(vx.mean(), vy.mean())))
    u, v = x_fits - crpix[0], y_fits - crpix[1]
    cd = None
    for _ in range(iteraciones):
        for _ in range(4):   # recorte de parejas malas
            if usar.sum() < 3:
                return None
            xi, eta = proyectar(ar, dec, ar0, dec0)
            matriz = np.column_stack([u, v, np.ones(n)])
            cx, *_ = np.linalg.lstsq(matriz[usar], xi[usar], rcond=None)
            cy, *_ = np.linalg.lstsq(matriz[usar], eta[usar], rcond=None)
            res = np.hypot(matriz @ cx - xi, matriz @ cy - eta) * 3600
            sigma = math.sqrt(np.mean(res[usar] ** 2)) if usar.sum() > 3 else float('inf')
            nuevo = res < max(recorte * sigma, 0.05)
            if usar.sum() <= 3 or np.array_equal(nuevo & usar, usar) or (nuevo & usar).sum() < 3:
                break
            usar &= nuevo
        cd = ((cx[0], cx[1]), (cy[0], cy[1]))
        ar0, dec0 = (float(t) for t in desproyectar(cx[2], cy[2], ar0, dec0))
        if abs(cx[2]) < 1e-10 and abs(cy[2]) < 1e-10:
            break
    # residuos finales (con el CRVAL ya movido)
    xi, eta = proyectar(ar, dec, ar0, dec0)
    xi_m = cd[0][0] * u + cd[0][1] * v
    eta_m = cd[1][0] * u + cd[1][1] * v
    res_x = (xi_m - xi) * 3600
    res_y = (eta_m - eta) * 3600
    res = np.hypot(res_x, res_y)
    det = cd[0][0] * cd[1][1] - cd[0][1] * cd[1][0]
    escala = 3600 * math.sqrt(abs(det))
    # Hacia dónde apunta el eje +Y de la foto en el cielo (ángulo desde el norte hacia el este)
    rotacion = math.degrees(math.atan2(cd[0][1], cd[1][1])) % 360.0
    cabecera = {'CTYPE1': 'RA---TAN', 'CTYPE2': 'DEC--TAN', 'CRPIX1': crpix[0], 'CRPIX2': crpix[1],
                'CRVAL1': ar0, 'CRVAL2': dec0, 'CD1_1': cd[0][0], 'CD1_2': cd[0][1],
                'CD2_1': cd[1][0], 'CD2_2': cd[1][1]}
    return {
        'cabecera': cabecera, 'usadas': usar, 'n_usadas': int(usar.sum()), 'n': n,
        'residuo_x': res_x, 'residuo_y': res_y, 'residuo': res,
        'rms': float(math.sqrt(np.mean(res[usar] ** 2))),
        'rms_x': float(math.sqrt(np.mean(res_x[usar] ** 2))), 'rms_y': float(math.sqrt(np.mean(res_y[usar] ** 2))),
        'escala': escala, 'rotacion': rotacion, 'espejada': det > 0,
    }


def emparejar(ar1, dec1, ar2, dec2, radio_arcsec):
    """Para cada estrella de la lista 1 busca la más cercana de la lista 2 a menos
    de radio_arcsec. Devuelve (índices en 2 o −1, distancias en segundos de arco)."""
    ar1, dec1, ar2, dec2 = (np.asarray(v, dtype=np.float64) for v in (ar1, dec1, ar2, dec2))
    indices = np.full(len(ar1), -1, dtype=int)
    distancias = np.full(len(ar1), np.inf)
    if len(ar1) == 0 or len(ar2) == 0:
        return indices, distancias
    v2 = np.column_stack([np.cos(np.radians(dec2)) * np.cos(np.radians(ar2)),
                          np.cos(np.radians(dec2)) * np.sin(np.radians(ar2)), np.sin(np.radians(dec2))])
    for inicio in range(0, len(ar1), 500):          # a trozos, para no gastar demasiada memoria
        a, d = np.radians(ar1[inicio:inicio + 500]), np.radians(dec1[inicio:inicio + 500])
        v1 = np.column_stack([np.cos(d) * np.cos(a), np.cos(d) * np.sin(a), np.sin(d)])
        coseno = np.clip(v1 @ v2.T, -1, 1)
        j = np.argmax(coseno, axis=1)
        dist = np.degrees(np.arccos(coseno[np.arange(len(j)), j])) * 3600
        ok = dist <= radio_arcsec
        indices[inicio:inicio + 500] = np.where(ok, j, -1)
        distancias[inicio:inicio + 500] = dist
    return indices, distancias


def epoca_juliana(jd):
    """Año juliano (por ejemplo 2026,19) a partir del día juliano."""
    return 2000.0 + (jd - 2451545.0) / 365.25


def aplicar_movimiento_propio(ar, dec, pmra, pmdec, epoca_catalogo, epoca):
    """Lleva posiciones de catálogo a la fecha de la foto (aproximación lineal, de sobra
    para unos años). pmra ya incluye cos(Dec), en mas/año, como en Gaia."""
    ar, dec = np.asarray(ar, dtype=np.float64), np.asarray(dec, dtype=np.float64)
    pmra = np.nan_to_num(np.asarray(pmra, dtype=np.float64))
    pmdec = np.nan_to_num(np.asarray(pmdec, dtype=np.float64))
    dt = epoca - epoca_catalogo
    nueva_dec = dec + pmdec * dt / 3.6e6
    nueva_ar = ar + pmra * dt / 3.6e6 / np.cos(np.radians(dec))
    return nueva_ar % 360.0, nueva_dec


# =============================================================================
# 5. OBJETOS MÓVILES
# =============================================================================


def ajuste_lineal_movimiento(tiempos, xs, ys, errores=None):
    """Ajusta x(t) = x0 + vx·t, y(t) = y0 + vy·t (t en días respecto al primero).
    Devuelve (x0, y0, vx, vy, residuos, rms) o None si hay menos de 2 puntos."""
    t = np.asarray(tiempos, dtype=np.float64)
    if len(t) < 2:
        return None
    t = t - t[0]
    xs, ys = np.asarray(xs, dtype=np.float64), np.asarray(ys, dtype=np.float64)
    if np.ptp(t) == 0:
        return None
    matriz = np.column_stack([np.ones_like(t), t])
    pesos = None if errores is None else 1 / np.maximum(np.asarray(errores, dtype=np.float64), 1e-6)
    a = matriz if pesos is None else matriz * pesos[:, None]
    bx = xs if pesos is None else xs * pesos
    by = ys if pesos is None else ys * pesos
    cx, *_ = np.linalg.lstsq(a, bx, rcond=None)
    cy, *_ = np.linalg.lstsq(a, by, rcond=None)
    residuos = np.hypot(matriz @ cx - xs, matriz @ cy - ys)
    rms = float(math.sqrt(np.mean(residuos ** 2))) if len(t) > 2 else 0.0
    return float(cx[0]), float(cy[0]), float(cx[1]), float(cy[1]), residuos, rms


def movimiento_en_cielo(jds, ars, decs):
    """Velocidad aparente de un objeto: (″/hora, ángulo de posición del movimiento,
    rms de los residuos en ″, residuos) a partir de sus posiciones en varias fotos."""
    jds = np.asarray(jds, dtype=np.float64)
    if len(jds) < 2:
        return None
    ar0, dec0 = float(np.mean(ars)), float(np.mean(decs))
    xi, eta = proyectar(np.asarray(ars), np.asarray(decs), ar0, dec0)
    ajuste = ajuste_lineal_movimiento(jds, xi * 3600, eta * 3600)
    if ajuste is None:
        return None
    _x0, _y0, vx, vy, residuos, rms = ajuste
    velocidad = math.hypot(vx, vy) / 24.0                    # ″/día -> ″/hora
    angulo = math.degrees(math.atan2(vx, vy)) % 360.0        # ξ hacia el este, η hacia el norte
    return velocidad, angulo, rms, residuos


def parece_movil(tiempos, xs, ys, errores, minimo_px=0.8):
    """¿Se mueve un objeto respecto a las estrellas? Sí si su desplazamiento total es
    claramente mayor que sus errores (3σ) y que minimo_px."""
    if len(tiempos) < 2:
        return False
    total = math.hypot(xs[-1] - xs[0], ys[-1] - ys[0])
    error = math.sqrt(errores[0] ** 2 + errores[-1] ** 2) if errores else 0.0
    return total > max(3 * error, minimo_px)


# =============================================================================
# 6. INDICADORES DE CALIDAD DE UNA MEDIDA
# =============================================================================
# Códigos (los textos traducidos están en Emillium.py):
#   SAT  saturada: el píxel más brillante llega al límite de la cámara (flujo y
#        centro poco fiables)
#   CON  contaminada: una estrella vecina mete más de un 1 % de luz en la apertura
#   BOR  toca el borde de la foto
#   SNR  señal/ruido baja (< 10): magnitud con más de 0,1 mag de error
#   FWHM tamaño raro (más de 1,6 veces o menos de 0,6 veces la mediana de la foto):
#        puede ser una galaxia, un rayo cósmico o dos estrellas juntas
#   ELO  estrella alargada (movida o doble sin separar)


def indicadores_calidad(saturada=False, fraccion_ajena=0.0, borde=False, snr=None, fwhm=None,
                        fwhm_referencia=None, elongacion=None):
    codigos = []
    if saturada:
        codigos.append('SAT')
    if fraccion_ajena and fraccion_ajena > 0.01:
        codigos.append('CON')
    if borde:
        codigos.append('BOR')
    if snr is not None and np.isfinite(snr) and snr < 10:
        codigos.append('SNR')
    if fwhm and fwhm_referencia and np.isfinite(fwhm_referencia) and fwhm_referencia > 0:
        if fwhm > 1.6 * fwhm_referencia or fwhm < 0.6 * fwhm_referencia:
            codigos.append('FWHM')
    if elongacion is not None and np.isfinite(elongacion) and elongacion > 1.6:
        codigos.append('ELO')
    return codigos


def elongacion_local(imagen, x, y, fondo, radio=4):
    """Elongación (eje mayor / eje menor) de la estrella en (x, y) por momentos."""
    alto, ancho = imagen.shape
    m = int(math.ceil(radio)) + 1
    cx, cy = int(round(x)), int(round(y))
    x0, x1 = max(cx - m, 0), min(cx + m + 1, ancho)
    y0, y1 = max(cy - m, 0), min(cy + m + 1, alto)
    trozo = np.nan_to_num(imagen[y0:y1, x0:x1].astype(np.float64) - fondo)
    yy, xx = np.mgrid[y0:y1, x0:x1]
    dentro = (xx - x) ** 2 + (yy - y) ** 2 <= radio ** 2
    luz = np.where(dentro, np.clip(trozo, 0, None), 0.0)
    total = luz.sum()
    if total <= 0:
        return float('nan')
    dx, dy = xx - x, yy - y
    mxx, myy, mxy = (luz * dx * dx).sum() / total, (luz * dy * dy).sum() / total, (luz * dx * dy).sum() / total
    traza, det = mxx + myy, math.sqrt(max((mxx - myy) ** 2 + 4 * mxy ** 2, 0))
    l1, l2 = (traza + det) / 2, max((traza - det) / 2, 1e-6)
    return math.sqrt(l1 / l2)


# =============================================================================
# 8. DATOS DE LA CABECERA: MASA DE AIRE Y APUNTADO APROXIMADO
# =============================================================================


def masa_de_aire(cabecera):
    """Masa de aire (cuánta atmósfera atraviesa la luz; 1 en el cénit) de la cabecera:
    AIRMASS o SECZ, o calculada con la altura del objeto (OBJCTALT/ALTITUDE)."""
    for clave in ('AIRMASS', 'SECZ', 'AIRMAS'):
        try:
            valor = float(cabecera.get(clave))
            if 1.0 <= valor < 40:
                return valor
        except (TypeError, ValueError):
            pass
    for clave in ('OBJCTALT', 'ALTITUDE', 'ALT-OBJ', 'ELEVATIO'):
        try:
            altura = float(cabecera.get(clave))
            if 3 < altura <= 90:
                # fórmula de Pickering (2002): buena también cerca del horizonte
                return 1 / math.sin(math.radians(altura + 244 / (165 + 47 * altura ** 1.1)))
        except (TypeError, ValueError):
            pass
    return None


def _angulo_texto(valor, horas):
    """'10 00 29.6' / '10:00:29.6' / '+22 27 24' / 150.12 -> grados (o None)."""
    if valor is None or isinstance(valor, bool):
        return None
    if isinstance(valor, (int, float)):
        return float(valor)
    texto = str(valor).strip().replace(',', '.')
    partes = [p for p in texto.replace(':', ' ').replace('h', ' ').replace('d', ' ').replace('m', ' ').replace('s', ' ').split() if p]
    try:
        if len(partes) == 1:
            return float(partes[0])
        signo = -1 if partes[0].startswith('-') else 1
        v = abs(float(partes[0])) + float(partes[1]) / 60 + (float(partes[2]) if len(partes) > 2 else 0) / 3600
        return signo * v * (15 if horas else 1)
    except (ValueError, IndexError):
        return None


def apuntado_aproximado(cabecera):
    """Dónde apuntaba el telescopio y la escala aproximada de la foto, si la cabecera lo
    dice (RA/DEC u OBJCTRA/OBJCTDEC; PIXSCALE, o XPIXSZ y FOCALLEN). Sirve para
    resolver la foto con Gaia. Devuelve {'ar', 'dec', 'escala'} (lo que no se sepa, None)."""
    ar = dec = None
    for clave_ar, clave_dec in (('RA', 'DEC'), ('OBJCTRA', 'OBJCTDEC'), ('RA_OBJ', 'DEC_OBJ'), ('TELRA', 'TELDEC'),
                                ('RA_TARG', 'DEC_TARG'), ('CAT-RA', 'CAT-DEC')):
        if clave_ar in cabecera and clave_dec in cabecera:
            valor_ar = cabecera.get(clave_ar)
            # en texto con espacios o ':' la AR va en horas; un número suelto, en grados
            horas = isinstance(valor_ar, str) and (':' in valor_ar or ' ' in valor_ar.strip())
            ar, dec = _angulo_texto(valor_ar, horas), _angulo_texto(cabecera.get(clave_dec), False)
            if ar is not None and dec is not None and 0 <= ar < 360 and -90 <= dec <= 90:
                break
            ar = dec = None
    if ar is None and str(cabecera.get('CTYPE1', '')).startswith('RA'):
        try:
            ar, dec = float(cabecera['CRVAL1']) % 360, float(cabecera['CRVAL2'])
        except (KeyError, TypeError, ValueError):
            pass
    escala = None
    for clave in ('PIXSCALE', 'SECPIX', 'SCALE', 'PIXSCAL1'):
        try:
            valor = float(cabecera.get(clave))
            if 0.01 < valor < 3600:
                escala = valor
                break
        except (TypeError, ValueError):
            pass
    if escala is None:
        try:
            pixel_um = float(cabecera.get('XPIXSZ', cabecera.get('PIXSIZE1')))
            focal_mm = float(cabecera.get('FOCALLEN'))
            if pixel_um > 0 and focal_mm > 0:
                escala = 206.265 * pixel_um / focal_mm
        except (TypeError, ValueError):
            pass
    return {'ar': ar, 'dec': dec, 'escala': escala}


# =============================================================================
# 9. RESOLVER LA FOTO CON UN CATÁLOGO (TRIÁNGULOS)
# =============================================================================
# Sin solución astrométrica, pero sabiendo más o menos hacia dónde apuntaba el
# telescopio, se buscan en el catálogo las mismas figuras que forman las estrellas de
# la foto. Un triángulo de estrellas tiene la misma FORMA aunque la foto esté girada,
# en espejo o a otra escala: sus lados ordenados (a ≤ b ≤ c) dan dos números, a/c y
# b/c, que no cambian. Se comparan los triángulos de las estrellas más brillantes de la
# foto con los del catálogo; cada pareja con la misma forma propone una transformación
# (giro, escala, espejo y desplazamiento), y se queda la que hace coincidir más
# estrellas. Después se afina con todas las parejas (ajustar_wcs).


def _triangulos(puntos, lado_minimo):
    """Todos los triángulos de una lista de puntos: invariantes (a/c, b/c) y vértices
    ordenados (opuesto al lado corto, al mediano y al largo)."""
    from itertools import combinations
    n = len(puntos)
    if n < 3:
        return np.zeros((0, 2)), np.zeros((0, 3), dtype=int)
    idx = np.array(list(combinations(range(n), 3)), dtype=int)
    p = puntos[idx]                                           # (T, 3, 2)
    lados = np.stack([np.hypot(*(p[:, 1] - p[:, 2]).T),       # opuesto al vértice 0
                      np.hypot(*(p[:, 0] - p[:, 2]).T),       # opuesto al 1
                      np.hypot(*(p[:, 0] - p[:, 1]).T)], axis=1)
    orden = np.argsort(lados, axis=1)
    ordenados = np.take_along_axis(lados, orden, axis=1)
    c = ordenados[:, 2]
    ok = (c >= lado_minimo) & (ordenados[:, 0] / np.maximum(c, 1e-9) > 0.12)
    # se evitan los casi isósceles: los vértices no se podrían ordenar sin ambigüedad
    ok &= (ordenados[:, 1] - ordenados[:, 0]) / np.maximum(c, 1e-9) > 0.03
    ok &= (ordenados[:, 2] - ordenados[:, 1]) / np.maximum(c, 1e-9) > 0.03
    invariantes = np.stack([ordenados[:, 0] / c, ordenados[:, 1] / c], axis=1)
    vertices = np.take_along_axis(idx, orden, axis=1)          # vértice opuesto al lado corto, mediano, largo
    return invariantes[ok], vertices[ok]


def _afin_3(origen, destino):
    """Transformación afín exacta que lleva 3 puntos a otros 3: destino = A·origen + t."""
    m = np.column_stack([origen, np.ones(3)])
    try:
        sol = np.linalg.solve(m, destino)
    except np.linalg.LinAlgError:
        return None, None
    return sol[:2].T, sol[2]


def resolver_campo(fuentes, ancho, alto, cat_ar, cat_dec, ar0, dec0, escala=None, tolerancia_escala=0.3,
                   n_foto=25, n_catalogo=60, tiempo_max=40.0, margen=1.5):
    """Busca la solución astrométrica de una foto sin WCS.

    fuentes: estrellas detectadas (x, y con el píxel 0 en 0), de más a menos brillante.
    cat_ar, cat_dec: estrellas del catálogo de la zona, de más a menos brillante.
    ar0, dec0: hacia dónde apuntaba el telescopio, más o menos.
    escala: ″/px aproximada. Si no se conoce, se prueban escalas de 0,3 a 12″/px.
    margen: el catálogo se busca hasta 'margen' veces el radio del campo (para tolerar
    que el apuntado no sea exacto).
    Devuelve el resultado de ajustar_wcs (con 'cabecera') o None si no lo encuentra."""
    import time
    cat_ar, cat_dec = np.asarray(cat_ar, dtype=np.float64), np.asarray(cat_dec, dtype=np.float64)
    escalas = [escala] if escala else [0.35 * 1.45 ** i for i in range(10)]     # 0,35 ... 9,8 ″/px
    inicio = time.perf_counter()
    encontradas = []
    for e in escalas:
        # Solo las estrellas del catálogo que pueden caer en la foto: si no, las más brillantes
        # de una zona grande quedarían casi todas fuera del campo.
        radio = math.hypot(ancho, alto) / 2 * e / 3600 * margen
        dist = np.degrees(np.arccos(np.clip(np.sin(np.radians(dec0)) * np.sin(np.radians(cat_dec))
                                            + np.cos(np.radians(dec0)) * np.cos(np.radians(cat_dec))
                                            * np.cos(np.radians(cat_ar - ar0)), -1, 1)))
        cerca = dist <= radio
        if cerca.sum() < 5:
            continue
        resto = max(tiempo_max - (time.perf_counter() - inicio), 2.0)
        r = _resolver_con_escala(fuentes, ancho, alto, cat_ar[cerca], cat_dec[cerca], ar0, dec0, e,
                                 tolerancia_escala if escala else 0.2, n_foto, n_catalogo, resto)
        if r is not None:
            if escala:
                return r
            encontradas.append(r)
        if time.perf_counter() - inicio > tiempo_max:
            break
    # Sin escala conocida se prueban todas y se queda la que empareja más estrellas: en un
    # campo con muchas estrellas, una escala equivocada puede juntar unas pocas por azar.
    return max(encontradas, key=lambda r: (r['n_usadas'], -r['rms'])) if encontradas else None


def _resolver_con_escala(fuentes, ancho, alto, cat_ar, cat_dec, ar0, dec0, escala, tolerancia_escala,
                         n_foto, n_catalogo, tiempo_max):
    import time
    inicio = time.perf_counter()
    det = np.array([(f['x'], f['y']) for f in fuentes[:n_foto]], dtype=np.float64)
    if len(det) < 5 or len(cat_ar) < 5:
        return None
    centro = np.array([(ancho - 1) / 2, (alto - 1) / 2])
    uv = det - centro
    xi, eta = proyectar(np.asarray(cat_ar[:n_catalogo]), np.asarray(cat_dec[:n_catalogo]), ar0, dec0)
    cat = np.column_stack([xi, eta]) * 3600.0                      # segundos de arco en el plano tangente
    cat_todo_xi, cat_todo_eta = proyectar(np.asarray(cat_ar[:4 * n_catalogo]), np.asarray(cat_dec[:4 * n_catalogo]), ar0, dec0)
    cat_todo = np.column_stack([cat_todo_xi, cat_todo_eta]) * 3600.0
    inv_d, ver_d = _triangulos(uv, lado_minimo=8.0)
    inv_c, ver_c = _triangulos(cat, lado_minimo=1.0)
    if len(inv_d) == 0 or len(inv_c) == 0:
        return None
    orden = np.argsort(inv_c[:, 0])
    inv_c, ver_c = inv_c[orden], ver_c[orden]
    tol = 0.008
    mejor, mejor_votos = None, 0
    radio_px = 2.5
    candidatos = 0
    for t in range(len(inv_d)):
        a, b = inv_d[t]
        i0, i1 = np.searchsorted(inv_c[:, 0], [a - tol, a + tol])
        if i1 <= i0:
            continue
        cercanos = np.nonzero(np.abs(inv_c[i0:i1, 1] - b) < tol)[0] + i0
        for k in cercanos:
            A, tr = _afin_3(uv[ver_d[t]], cat[ver_c[k]])
            if A is None:
                continue
            sv = np.linalg.svd(A, compute_uv=False)
            if sv[1] <= 0 or sv[0] / sv[1] > 1.04:                # tiene que ser casi un giro + escala
                continue
            s = math.sqrt(abs(np.linalg.det(A)))
            if escala is not None and not (1 - tolerancia_escala) < s / escala < (1 + tolerancia_escala):
                continue
            candidatos += 1
            proyectadas = uv @ A.T + tr
            d = np.hypot(proyectadas[:, None, 0] - cat_todo[None, :, 0], proyectadas[:, None, 1] - cat_todo[None, :, 1])
            votos = int((d.min(axis=1) < radio_px * s).sum())
            if votos > mejor_votos:
                mejor, mejor_votos = (A, tr, s), votos
            if candidatos > 30000 or time.perf_counter() - inicio > tiempo_max:
                break
        if candidatos > 30000 or time.perf_counter() - inicio > tiempo_max:
            break
    if mejor is None or mejor_votos < max(5, int(0.3 * len(uv))):
        return None
    A, tr, s = mejor
    # Afinado: todas las estrellas detectadas contra todo el catálogo, varias veces
    todas = np.array([(f['x'], f['y']) for f in fuentes[:400]], dtype=np.float64)
    xi_t, eta_t = proyectar(np.asarray(cat_ar), np.asarray(cat_dec), ar0, dec0)
    cat_t = np.column_stack([xi_t, eta_t]) * 3600.0
    resultado = None
    for radio in (3.0, 2.0, 1.5):
        proyectadas = (todas - centro) @ A.T + tr
        d = np.hypot(proyectadas[:, None, 0] - cat_t[None, :, 0], proyectadas[:, None, 1] - cat_t[None, :, 1])
        j = d.argmin(axis=1)
        ok = d[np.arange(len(j)), j] < radio * s
        # cada estrella del catálogo, una sola vez
        usados, pares = set(), []
        for i in np.argsort(d[np.arange(len(j)), j]):
            if ok[i] and j[i] not in usados:
                usados.add(j[i])
                pares.append((i, j[i]))
        if len(pares) < 5:
            return resultado
        ii = np.array([p[0] for p in pares])
        jj = np.array([p[1] for p in pares])
        resultado = ajustar_wcs(todas[ii, 0] + 1, todas[ii, 1] + 1, np.asarray(cat_ar)[jj], np.asarray(cat_dec)[jj],
                                crpix=((ancho + 1) / 2, (alto + 1) / 2))
        if resultado is None:
            return None
        # transformación afín actualizada (para el siguiente emparejado), en el plano de (ar0, dec0)
        cab = resultado['cabecera']
        cd = np.array([[cab['CD1_1'], cab['CD1_2']], [cab['CD2_1'], cab['CD2_2']]]) * 3600
        x0, y0 = proyectar(cab['CRVAL1'], cab['CRVAL2'], ar0, dec0)
        A, tr = cd, np.array([float(x0), float(y0)]) * 3600
        s = math.sqrt(abs(np.linalg.det(cd)))
    # ¿Es creíble? Una buena solución empareja una parte importante de las estrellas que
    # deberían verse: las detectadas o las del catálogo que caen dentro de la foto.
    pix = (cat_t - tr) @ np.linalg.inv(A).T + centro
    en_foto = int(((pix[:, 0] >= 0) & (pix[:, 0] < ancho) & (pix[:, 1] >= 0) & (pix[:, 1] < alto)).sum())
    if resultado['n_usadas'] < max(6, 0.2 * min(len(todas), en_foto)):
        return None
    resultado['votos'] = mejor_votos
    resultado['candidatos'] = candidatos
    return resultado


# =============================================================================
# 10. CALIBRACIÓN FOTOMÉTRICA CON TODO EL CAMPO Y EXTINCIÓN
# =============================================================================


def calibrar_campo(m_inst, m_cat, color=None, errores=None, recorte=3.0):
    """Ajusta  m_catálogo − m_instrumental = ZP + T·color  con muchas estrellas del campo
    (mínimos cuadrados ponderados, quitando las que se apartan más de 3σ: variables,
    estrellas mal emparejadas...). Sin color, solo el punto cero.
    Devuelve ZP, T, sus errores, el rms y qué estrellas se han usado."""
    m_inst, m_cat = np.asarray(m_inst, dtype=np.float64), np.asarray(m_cat, dtype=np.float64)
    con_color = color is not None
    c = np.asarray(color, dtype=np.float64) if con_color else np.zeros_like(m_inst)
    e = np.maximum(np.asarray(errores, dtype=np.float64), 0.005) if errores is not None else np.full_like(m_inst, 0.02)
    usar = np.isfinite(m_inst) & np.isfinite(m_cat) & np.isfinite(c)
    y = m_cat - m_inst
    for _ in range(6):
        if usar.sum() < (3 if con_color else 2):
            return None
        a = np.column_stack([np.ones(usar.sum()), c[usar]]) if con_color else np.ones((usar.sum(), 1))
        w = 1 / e[usar]
        coef, *_ = np.linalg.lstsq(a * w[:, None], y[usar] * w, rcond=None)
        modelo = coef[0] + (coef[1] * c if con_color else 0)
        res = y - modelo
        rms = float(np.sqrt(np.mean(res[usar] ** 2)))
        # el recorte usa una σ robusta (MAD): unas pocas estrellas muy malas no la inflan
        s_robusta = sigma_robusta(res[usar]) or rms
        nuevo = usar & (np.abs(res) < max(recorte * s_robusta, 0.02))
        if nuevo.sum() == usar.sum():
            break
        usar = nuevo
    a = np.column_stack([np.ones(usar.sum()), c[usar]]) if con_color else np.ones((usar.sum(), 1))
    cov = np.linalg.pinv(a.T @ a) * rms ** 2
    return {'zp': float(coef[0]), 'sigma_zp': float(math.sqrt(cov[0, 0])),
            't': float(coef[1]) if con_color else 0.0, 'sigma_t': float(math.sqrt(cov[1, 1])) if con_color else 0.0,
            'rms': rms, 'n': int(usar.sum()), 'usadas': usar, 'residuos': y - modelo}


def ajustar_extincion(masas, puntos_cero, errores=None):
    """Coeficiente de extinción k (mag por masa de aire): ZP(X) = ZP0 − k·X, con el punto
    cero de cada foto y su masa de aire. (Una estrella pierde k magnitudes por cada
    masa de aire de atmósfera que atraviesa su luz.)"""
    x, y = np.asarray(masas, dtype=np.float64), np.asarray(puntos_cero, dtype=np.float64)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3 or np.ptp(x[ok]) < 0.05:
        return None
    w = 1 / np.maximum(np.asarray(errores, dtype=np.float64)[ok], 0.003) if errores is not None else np.ones(ok.sum())
    a = np.column_stack([np.ones(ok.sum()), x[ok]])
    coef, *_ = np.linalg.lstsq(a * w[:, None], y[ok] * w, rcond=None)
    res = y[ok] - a @ coef
    gl = max(ok.sum() - 2, 1)
    cov = np.linalg.pinv((a * w[:, None]).T @ (a * w[:, None])) * float((res ** 2 * w ** 2).sum() / gl)
    return {'k': float(-coef[1]), 'sigma_k': float(math.sqrt(max(cov[1, 1], 0))), 'zp0': float(coef[0]),
            'rms': float(np.sqrt(np.mean(res ** 2))), 'n': int(ok.sum())}


def curva_crecimiento(imagen, x, y, radios, r_int, r_ext, ganancia=None):
    """Flujo y señal/ruido en aperturas de varios radios (curva de crecimiento).
    El radio con más señal/ruido es el mejor para fotometría de estrellas débiles."""
    filas = []
    for r in radios:
        if r >= r_int:
            break
        f = fotometria_apertura(imagen, x, y, r, r_int, r_ext, ganancia)
        if f is not None:
            filas.append((float(r), f['flujo'], f['snr']))
    return filas


def radio_optimo(curvas):
    """Radio que da más señal/ruido, como mediana de las curvas de varias estrellas."""
    mejores = [max(c, key=lambda fila: fila[2])[0] for c in curvas if c]
    return float(np.median(mejores)) if mejores else None


# =============================================================================
# 11. BUSCAR OBJETOS MÓVILES AUTOMÁTICAMENTE EN UNA SERIE
# =============================================================================


def buscar_moviles(listas, desplazamientos, tiempos, tolerancia=1.5, minimo_fotos=3, velocidad_minima=0.5):
    """Busca objetos que se mueven en línea recta entre las fotos de una serie.

    listas: para cada foto, sus estrellas detectadas (x, y, flujo).
    desplazamientos: (dx, dy) de cada foto respecto a la primera (desplazamiento_entre).
    tiempos: instante de cada foto (días).

    1) Todo se pasa a las coordenadas de la primera foto.
    2) Las fuentes que están quietas (en el mismo sitio en la mitad de las fotos o más)
       son estrellas: se quitan.
    3) Con lo que queda, cada pareja (una de la primera foto donde aparece y otra de la
       última) define una velocidad; se busca el objeto en las fotos intermedias donde
       debería estar. Si aparece en suficientes fotos, es un objeto móvil.
    Devuelve una lista de pistas: posiciones por foto, velocidad (px/día) y número de fotos."""
    n = len(listas)
    if n < minimo_fotos:
        return []
    puntos = []
    for k, (lista, d) in enumerate(zip(listas, desplazamientos)):
        if d is None:
            puntos.append(np.zeros((0, 3)))
            continue
        p = np.array([(f['x'] - d[0], f['y'] - d[1], f['flujo']) for f in lista], dtype=np.float64).reshape(-1, 3)
        puntos.append(p)
    # 2) quietas
    todas = np.vstack([p[:, :2] for p in puntos if len(p)]) if any(len(p) for p in puntos) else np.zeros((0, 2))
    transitorias = []
    for k, p in enumerate(puntos):
        if not len(p):
            transitorias.append(np.zeros((0, 3)))
            continue
        cuenta = np.zeros(len(p), dtype=int)
        for j, q in enumerate(puntos):
            if j == k or not len(q):
                continue
            d = np.hypot(p[:, None, 0] - q[None, :, 0], p[:, None, 1] - q[None, :, 1])
            cuenta += (d.min(axis=1) < tolerancia).astype(int)
        quieta = cuenta >= max(1, (n - 1) // 2)
        transitorias.append(p[~quieta])
    del todas
    # 3) pistas: de la primera foto (o la segunda) a la última (o la penúltima)
    t = np.asarray(tiempos, dtype=np.float64)
    pistas = []
    for a in range(min(2, n)):
        for b in range(n - 1, max(n - 3, a), -1):
            if t[b] == t[a]:
                continue
            for pa in transitorias[a]:
                for pb in transitorias[b]:
                    v = (pb[:2] - pa[:2]) / (t[b] - t[a])
                    rapidez = math.hypot(*v)
                    if rapidez * abs(t[b] - t[a]) < velocidad_minima:
                        continue
                    posiciones = {}
                    for k in range(n):
                        esperada = pa[:2] + v * (t[k] - t[a])
                        q = transitorias[k]
                        if not len(q):
                            continue
                        d = np.hypot(q[:, 0] - esperada[0], q[:, 1] - esperada[1])
                        i = int(d.argmin())
                        if d[i] < tolerancia * 1.5:
                            posiciones[k] = (float(q[i, 0]), float(q[i, 1]), float(q[i, 2]))
                    if len(posiciones) >= max(minimo_fotos, int(0.6 * n)):
                        pistas.append({'posiciones': posiciones, 'velocidad': v.tolist(), 'n': len(posiciones),
                                       'rapidez_px_dia': rapidez})
    # quitar repetidas (la misma pista encontrada desde parejas distintas)
    unicas = []
    for p in sorted(pistas, key=lambda p: -p['n']):
        k0 = min(p['posiciones'])
        x0, y0 = p['posiciones'][k0][:2]
        if any(k0 in q['posiciones'] and math.hypot(q['posiciones'][k0][0] - x0, q['posiciones'][k0][1] - y0) < 3
               for q in unicas):
            continue
        unicas.append(p)
    # posiciones de vuelta a las coordenadas de cada foto
    for p in unicas:
        p['posiciones'] = {k: (x + desplazamientos[k][0], y + desplazamientos[k][1], f)
                           for k, (x, y, f) in p['posiciones'].items()}
    return unicas


# =============================================================================
# 12. TIEMPO HELIOCÉNTRICO (HJD)
# =============================================================================
# La luz de una estrella llega a la Tierra antes o después según dónde esté la Tierra en
# su órbita (hasta ±8,3 minutos). Para comparar instantes de fotos de distintas fechas
# (el mínimo de un eclipse, el periodo de una variable) se usa el día juliano
# heliocéntrico: el instante en que la luz habría llegado al centro del Sol.
#   HJD = JD − (r/c)·[cos δ cos α cos L + (sen ε sen δ + cos ε cos δ sen α) sen L]
# con L la longitud geométrica del Sol, r la distancia Tierra-Sol y ε la oblicuidad.
# Sol de baja precisión (Meeus, cap. 25): ~0,01°, de sobra (error < 1 s en el HJD).


def posicion_sol(jd):
    """Longitud geométrica del Sol (grados), distancia Tierra-Sol (UA) y oblicuidad (grados)."""
    t = (jd - 2451545.0) / 36525.0
    l0 = (280.46646 + 36000.76983 * t + 0.0003032 * t * t) % 360
    m = math.radians((357.52911 + 35999.05029 * t - 0.0001537 * t * t) % 360)
    e = 0.016708634 - 0.000042037 * t - 0.0000001267 * t * t
    c = ((1.914602 - 0.004817 * t - 0.000014 * t * t) * math.sin(m) + (0.019993 - 0.000101 * t) * math.sin(2 * m)
         + 0.000289 * math.sin(3 * m))
    # La longitud sale respecto al equinoccio de la fecha; las coordenadas de las estrellas son
    # J2000, así que se pasa a J2000 quitando la precesión (Meeus: 0,01397° por año).
    longitud = (l0 + c - 1.397 * t) % 360
    anomalia = m + math.radians(c)
    r = 1.000001018 * (1 - e * e) / (1 + e * math.cos(anomalia))
    oblicuidad = 23.4392911          # oblicuidad de J2000
    return longitud, r, oblicuidad


def jd_heliocentrico(jd, ar, dec):
    """Día juliano heliocéntrico para un objeto en (ar, dec) en grados."""
    longitud, r, eps = posicion_sol(jd)
    a, d, l, e = (math.radians(v) for v in (ar, dec, longitud, eps))
    correccion = 0.0057755183 * r * (math.cos(d) * math.cos(a) * math.cos(l)
                                     + (math.sin(e) * math.sin(d) + math.cos(e) * math.cos(d) * math.sin(a)) * math.sin(l))
    return jd - correccion


# -----------------------------------------------------------------------------
# BJD_TDB: el estándar actual para cronometrar eclipses y tránsitos.
# Dos diferencias con el HJD:
#  1) El centro de referencia es el BARICENTRO del sistema solar (el centro de masas, que
#     por culpa de Júpiter y Saturno está hasta a ~1 radio solar del centro del Sol: hasta
#     ±4 s de diferencia).
#  2) La escala de tiempo es TDB (uniforme), no UTC: TDB ≈ UTC + segundos intercalares +
#     32,184 s (unos 69 s en 2026). Sin esto, dos medidas separadas por un segundo
#     intercalar no se podrían comparar.
#   BJD_TDB = JD_UTC + (TT − UTC) + (TDB − TT) + (posición de la Tierra respecto al
#             baricentro · dirección de la estrella) / c
# Los planetas gigantes salen de sus elementos orbitales medios (JPL, Standish 1992),
# de sobra para mover el Sol con un error de milisegundos.
# -----------------------------------------------------------------------------

# Segundos intercalares: (día juliano UTC desde el que valen, TAI − UTC en segundos)
_TAI_UTC = [(2441317.5, 10), (2441499.5, 11), (2441683.5, 12), (2442048.5, 13), (2442413.5, 14), (2442778.5, 15),
            (2443144.5, 16), (2443509.5, 17), (2443874.5, 18), (2444239.5, 19), (2444786.5, 20), (2445151.5, 21),
            (2445516.5, 22), (2446247.5, 23), (2447161.5, 24), (2447892.5, 25), (2448257.5, 26), (2448804.5, 27),
            (2449169.5, 28), (2449534.5, 29), (2450083.5, 30), (2450630.5, 31), (2451179.5, 32), (2453736.5, 33),
            (2454832.5, 34), (2456109.5, 35), (2457204.5, 36), (2457754.5, 37)]

# Elementos medios (época J2000 y variación por siglo): a (UA), e, i, L, ϖ, Ω (grados);
# y masa del Sol / masa del planeta.
_PLANETAS = {
    'Júpiter': ((5.20288700, -0.00011607), (0.04838624, -0.00013253), (1.30439695, -0.00183714),
                (34.39644051, 3034.74612775), (14.72847983, 0.21252668), (100.47390909, 0.20469106), 1047.3486),
    'Saturno': ((9.53667594, -0.00125060), (0.05386179, -0.00050991), (2.48599187, 0.00193609),
                (49.95424423, 1222.49362201), (92.59887831, -0.41897216), (113.66242448, -0.28867794), 3497.898),
    'Urano': ((19.18916464, -0.00196176), (0.04725744, -0.00004397), (0.77263783, -0.00242939),
              (313.23810451, 428.48202785), (170.95427630, 0.40805281), (74.01692503, 0.04240589), 22902.98),
    'Neptuno': ((30.06992276, 0.00026291), (0.00859048, 0.00005105), (1.77004347, 0.00035372),
                (-55.12002969, 218.45945325), (44.96476227, -0.32241464), (131.78422574, -0.00508664), 19412.24),
}
_OBLICUIDAD_J2000 = math.radians(23.4392911)


def tt_menos_utc(jd_utc):
    """TT − UTC en segundos (segundos intercalares + 32,184 s)."""
    tai = _TAI_UTC[0][1]
    for inicio, valor in _TAI_UTC:
        if jd_utc >= inicio:
            tai = valor
    return tai + 32.184


def _posicion_planeta(elementos, t):
    """Posición heliocéntrica eclíptica J2000 (UA) de un planeta, t en siglos desde J2000."""
    a, e, i, l, w, o = (v0 + v1 * t for v0, v1 in elementos[:6])
    i, o = math.radians(i), math.radians(o)
    perihelio = math.radians(w) - o
    m = math.radians((l - w) % 360)
    excentrica = m
    for _ in range(12):                                   # ecuación de Kepler (Newton)
        excentrica -= (excentrica - e * math.sin(excentrica) - m) / (1 - e * math.cos(excentrica))
    xp, yp = a * (math.cos(excentrica) - e), a * math.sqrt(1 - e * e) * math.sin(excentrica)
    cw, sw, co, so, ci, si = (math.cos(perihelio), math.sin(perihelio), math.cos(o), math.sin(o),
                              math.cos(i), math.sin(i))
    return ((cw * co - sw * so * ci) * xp + (-sw * co - cw * so * ci) * yp,
            (cw * so + sw * co * ci) * xp + (-sw * so + cw * co * ci) * yp,
            (sw * si) * xp + (cw * si) * yp)


def posicion_tierra_baricentrica(jd):
    """Posición de la Tierra respecto al baricentro del sistema solar, en UA, en
    coordenadas ecuatoriales J2000 (x, y, z)."""
    t = (jd - 2451545.0) / 36525.0
    longitud, r, _eps = posicion_sol(jd)
    lon = math.radians(longitud)
    tierra = [-r * math.cos(lon), -r * math.sin(lon), 0.0]          # Tierra respecto al Sol (eclíptica)
    suma_masas = 1.0 + sum(1.0 / p[6] for p in _PLANETAS.values())
    for elementos in _PLANETAS.values():                             # el Sol respecto al baricentro
        px, py, pz = _posicion_planeta(elementos, t)
        fraccion = (1.0 / elementos[6]) / suma_masas
        tierra[0] -= fraccion * px
        tierra[1] -= fraccion * py
        tierra[2] -= fraccion * pz
    ce, se = math.cos(_OBLICUIDAD_J2000), math.sin(_OBLICUIDAD_J2000)
    return tierra[0], ce * tierra[1] - se * tierra[2], se * tierra[1] + ce * tierra[2]


def bjd_tdb(jd_utc, ar, dec):
    """BJD_TDB (día juliano baricéntrico en tiempo dinámico) de un instante JD en UTC
    para un objeto en (ar, dec) en grados."""
    jd_tt = jd_utc + tt_menos_utc(jd_utc) / 86400.0
    g = math.radians(357.53 + 0.98560028 * (jd_tt - 2451545.0))
    jd_tdb = jd_tt + (0.001657 * math.sin(g) + 0.000014 * math.sin(2 * g)) / 86400.0
    x, y, z = posicion_tierra_baricentrica(jd_tdb)
    a, d = math.radians(ar), math.radians(dec)
    proyeccion = x * math.cos(d) * math.cos(a) + y * math.cos(d) * math.sin(a) + z * math.sin(d)
    return jd_tdb + 0.0057755183 * proyeccion


# =============================================================================
# 13. PERIODOS Y MÍNIMOS DE ESTRELLAS VARIABLES
# =============================================================================


def periodograma(tiempos, magnitudes, errores=None, periodo_min=None, periodo_max=None, n=20000):
    """Periodograma de Lomb-Scargle generalizado (con media flotante y pesos).

    Para cada periodo de prueba ajusta una sinusoide a los datos y mide qué parte de la
    variación explica (potencia de 0 a 1). Los periodos con potencia alta son candidatos.
    Devuelve (periodos, potencia, mejor periodo)."""
    t = np.asarray(tiempos, dtype=np.float64)
    y = np.asarray(magnitudes, dtype=np.float64)
    w = 1 / np.maximum(np.asarray(errores, dtype=np.float64), 1e-4) ** 2 if errores is not None else np.ones_like(y)
    ok = np.isfinite(t) & np.isfinite(y) & np.isfinite(w)
    t, y, w = t[ok], y[ok], w[ok]
    if len(t) < 6:
        return None
    w = w / w.sum()
    duracion = t.max() - t.min()
    if duracion <= 0:
        return None
    diferencias = np.diff(np.sort(t))
    paso = np.median(diferencias[diferencias > 0]) if (diferencias > 0).any() else duracion / len(t)
    periodo_min = periodo_min or max(2 * paso, duracion / 2000)
    periodo_max = periodo_max or duracion
    frecuencias = np.linspace(1 / periodo_max, 1 / periodo_min, n)
    y_media = (w * y).sum()
    yy = (w * (y - y_media) ** 2).sum()
    potencia = np.empty(n)
    for inicio in range(0, n, 2000):                     # a trozos, para no gastar mucha memoria
        f = frecuencias[inicio:inicio + 2000, None]
        arg = 2 * np.pi * f * t[None, :]
        c, s = np.cos(arg), np.sin(arg)
        C_ = (w * c).sum(axis=1); S_ = (w * s).sum(axis=1)
        YC = (w * (y - y_media) * c).sum(axis=1); YS = (w * (y - y_media) * s).sum(axis=1)
        CC = (w * c * c).sum(axis=1) - C_ ** 2; SS = (w * s * s).sum(axis=1) - S_ ** 2
        CS = (w * c * s).sum(axis=1) - C_ * S_
        D = CC * SS - CS ** 2
        potencia[inicio:inicio + 2000] = (SS * YC ** 2 + CC * YS ** 2 - 2 * CS * YC * YS) / (yy * np.where(D > 0, D, np.inf))
    periodos = 1 / frecuencias
    return periodos, potencia, float(periodos[int(np.argmax(potencia))])


def fase(tiempos, periodo, epoca=None):
    """Fase (0 a 1) de cada instante para un periodo dado."""
    t = np.asarray(tiempos, dtype=np.float64)
    epoca = t.min() if epoca is None else epoca
    return ((t - epoca) / periodo) % 1.0


def minimo_kvw(tiempos, magnitudes, puntos=None):
    """Instante del mínimo de un eclipse por el método de Kwee y van Woerden (1956).

    Se prueban instantes alrededor del punto más débil: para cada uno se refleja la curva
    respecto a ese instante y se mide cuánto se parecen la rama de bajada y la de subida
    (suma de cuadrados S). El mínimo de S está en el centro del eclipse; una parábola en
    S(T) da el instante y su error.
    Las magnitudes más grandes son las más débiles (el mínimo de luz es el máximo de magnitud).
    Devuelve (instante, error) o None."""
    t = np.asarray(tiempos, dtype=np.float64)
    m = np.asarray(magnitudes, dtype=np.float64)
    ok = np.isfinite(t) & np.isfinite(m)
    t, m = t[ok], m[ok]
    orden = np.argsort(t)
    t, m = t[orden], m[orden]
    if len(t) < 7:
        return None
    k = int(np.argmax(m))                                    # punto más débil
    if k < 2 or k > len(t) - 3:
        return None                                          # el mínimo no está dentro de los datos
    n = min(k, len(t) - 1 - k)
    paso = np.median(np.diff(t))
    candidatos = t[k] + np.linspace(-2 * paso, 2 * paso, 41)
    sumas = []
    for T in candidatos:
        reflejados = np.interp(2 * T - t, t, m, left=np.nan, right=np.nan)
        cerca = np.abs(t - T) <= n * paso
        validos = cerca & np.isfinite(reflejados)
        if validos.sum() < 4:
            sumas.append(np.nan)
            continue
        sumas.append(float(((m[validos] - reflejados[validos]) ** 2).sum() / validos.sum()))
    sumas = np.array(sumas)
    ok = np.isfinite(sumas)
    if ok.sum() < 5:
        return None
    i = int(np.nanargmin(sumas))
    lo, hi = max(i - 5, 0), min(i + 6, len(sumas))
    x, s = candidatos[lo:hi][ok[lo:hi]], sumas[lo:hi][ok[lo:hi]]
    if len(x) < 3:
        return None
    a, b, c = np.polyfit(x - candidatos[i], s, 2)
    if a <= 0:
        return None
    t0 = candidatos[i] - b / (2 * a)
    s_min = c - b * b / (4 * a)
    z = max(n // 2, 2)                                      # puntos independientes (aprox. de KvW)
    error = math.sqrt(max(s_min, 1e-12) / (a * (z - 1))) if z > 1 else float('nan')
    return float(t0), float(error)


# =============================================================================
# 14. SOLUCIÓN ASTROMÉTRICA CON DISTORSIÓN (SIP)
# =============================================================================
# En fotos de campo ancho (teleobjetivos, reductores) la óptica deforma un poco el cielo:
# una solución lineal deja residuos que crecen hacia los bordes. SIP (Simple Imaging
# Polynomial) corrige los píxeles con polinomios antes de la matriz CD:
#   u' = u + Σ A_pq·u^p·v^q     v' = v + Σ B_pq·u^p·v^q     (2 ≤ p + q ≤ orden)


def ajustar_sip(x_fits, y_fits, ar, dec, orden=3, crpix=None, recorte=3.0, iteraciones=6):
    """Ajusta una solución TAN con distorsión SIP.

    Una solución TAN-SIP equivale exactamente a dos polinomios en el plano tangente:
        ξ = CD11·(u + f) + CD12·(v + g)     η = CD21·(u + f) + CD22·(v + g)
    con f y g los polinomios SIP. Así que se ajustan directamente por mínimos cuadrados
    ξ(u, v) y η(u, v) como polinomios completos (de grado 0 a 'orden'): la parte lineal da
    la matriz CD, la parte constante dice cuánto mover el punto de tangencia, y los
    coeficientes de grado ≥ 2, pasados por la inversa de CD, son los de SIP.
    Necesita bastantes estrellas (unas 3 por coeficiente). Devuelve un diccionario como
    ajustar_wcs (con la cabecera -SIP) o None."""
    x_fits, y_fits, ar, dec = (np.asarray(v, dtype=np.float64) for v in (x_fits, y_fits, ar, dec))
    terminos = [(p, q) for p in range(orden + 1) for q in range(orden + 1) if p + q <= orden]
    altos = [(p, q) for p, q in terminos if p + q >= 2]
    n = len(x_fits)
    if n < 3 * len(terminos):
        return None
    lineal = ajustar_wcs(x_fits, y_fits, ar, dec, crpix=crpix)
    if lineal is None:
        return None
    crpix = (lineal['cabecera']['CRPIX1'], lineal['cabecera']['CRPIX2'])
    ar0, dec0 = lineal['cabecera']['CRVAL1'], lineal['cabecera']['CRVAL2']
    u, v = x_fits - crpix[0], y_fits - crpix[1]
    escala = max(np.abs(u).max(), np.abs(v).max(), 1.0)          # para que los números sean manejables
    un, vn = u / escala, v / escala
    matriz = np.column_stack([un ** p * vn ** q for p, q in terminos])
    usar = np.ones(n, dtype=bool)
    for _ in range(iteraciones):
        xi, eta = proyectar(ar, dec, ar0, dec0)
        for _ in range(5):          # recorte de parejas malas
            cx, *_ = np.linalg.lstsq(matriz[usar], xi[usar], rcond=None)
            cy, *_ = np.linalg.lstsq(matriz[usar], eta[usar], rcond=None)
            res = np.hypot(matriz @ cx - xi, matriz @ cy - eta) * 3600
            sigma = math.sqrt(np.mean(res[usar] ** 2))
            nuevo = res < max(recorte * sigma, 0.05)
            if (nuevo & usar).sum() < 3 * len(terminos) or np.array_equal(nuevo & usar, usar):
                break
            usar &= nuevo
        # el término constante mueve el punto de tangencia (CRVAL); se repite hasta que no se mueve
        i0 = terminos.index((0, 0))
        if abs(cx[i0]) < 1e-11 and abs(cy[i0]) < 1e-11:
            break
        ar0, dec0 = (float(t) for t in desproyectar(cx[i0], cy[i0], ar0, dec0))
    i10, i01 = terminos.index((1, 0)), terminos.index((0, 1))
    cd = np.array([[cx[i10], cx[i01]], [cy[i10], cy[i01]]]) / escala
    inv = np.linalg.inv(cd)
    cabecera = {'CTYPE1': 'RA---TAN-SIP', 'CTYPE2': 'DEC--TAN-SIP', 'CRPIX1': crpix[0], 'CRPIX2': crpix[1],
                'CRVAL1': ar0, 'CRVAL2': dec0, 'CD1_1': cd[0, 0], 'CD1_2': cd[0, 1], 'CD2_1': cd[1, 0], 'CD2_2': cd[1, 1],
                'A_ORDER': orden, 'B_ORDER': orden}
    for p, q in altos:
        k = terminos.index((p, q))
        alto = np.array([cx[k], cy[k]]) / escala ** (p + q)         # coeficientes en píxeles sin normalizar
        a, b = inv @ alto
        cabecera[f'A_{p}_{q}'] = float(a)
        cabecera[f'B_{p}_{q}'] = float(b)
    xi, eta = proyectar(ar, dec, ar0, dec0)
    res_x = (matriz @ cx - xi) * 3600
    res_y = (matriz @ cy - eta) * 3600
    res = np.hypot(res_x, res_y)
    det = np.linalg.det(cd)
    return {'cabecera': cabecera, 'usadas': usar, 'n_usadas': int(usar.sum()), 'n': n,
            'residuo_x': res_x, 'residuo_y': res_y, 'residuo': res,
            'rms': float(math.sqrt(np.mean(res[usar] ** 2))),
            'rms_x': float(math.sqrt(np.mean(res_x[usar] ** 2))), 'rms_y': float(math.sqrt(np.mean(res_y[usar] ** 2))),
            'escala': 3600 * math.sqrt(abs(det)), 'rotacion': math.degrees(math.atan2(cd[0, 1], cd[1, 1])) % 360.0,
            'espejada': det > 0, 'orden': orden}


# =============================================================================
# 15. CORRECCIÓN DE LA CONTAMINACIÓN Y RECORTES
# =============================================================================


def flujo_corregido(flujo, sigma_flujo, fraccion_ajena):
    """Flujo sin la luz estimada de las vecinas. La estimación (estrellas gaussianas)
    no es perfecta: se añade al error la mitad de lo que se ha restado."""
    if not fraccion_ajena or fraccion_ajena <= 0 or flujo is None or flujo <= 0:
        return flujo, sigma_flujo
    # la fracción es respecto al flujo medido (que ya incluye la luz de la vecina)
    restado = flujo * min(fraccion_ajena, 0.9)
    return flujo - restado, math.hypot(sigma_flujo, 0.5 * restado)


def recorte(imagen, x, y, radio):
    """Trozo cuadrado de la imagen centrado en (x, y), rellenado con NaN fuera de la foto."""
    alto, ancho = imagen.shape
    cx, cy = int(round(x)), int(round(y))
    lado = 2 * radio + 1
    salida = np.full((lado, lado), np.nan, dtype=np.float32)
    x0, x1 = max(cx - radio, 0), min(cx + radio + 1, ancho)
    y0, y1 = max(cy - radio, 0), min(cy + radio + 1, alto)
    if x0 < x1 and y0 < y1:
        salida[y0 - (cy - radio):y1 - (cy - radio), x0 - (cx - radio):x1 - (cx - radio)] = imagen[y0:y1, x0:x1]
    return salida


# =============================================================================
# 16. ALINEAR Y APILAR FOTOS
# =============================================================================
# Alinear: se buscan las mismas estrellas en las dos fotos con triángulos (como al
# resolver una foto con un catálogo, pero de píxeles a píxeles). Cada pareja de triángulos
# con la misma forma propone una transformación (desplazamiento, giro, espejo y un poco de
# escala); se queda la que hace coincidir más estrellas y se afina por mínimos cuadrados
# con todas las parejas (transformación afín: 6 números). Precisión: centésimas de píxel.
#
# Remuestrear: para cada píxel de la foto de referencia se calcula dónde cae en la otra foto
# y se interpola (bilineal: con los 4 píxeles de alrededor, que conserva el flujo total).
#
# Apilar: se suman muchas fotos alineadas para que el ruido baje (√N veces) y salgan
# objetos débiles. Antes se igualan los fondos del cielo. Con la "media con recorte" se
# descartan, píxel a píxel, los valores raros (rayos cósmicos, satélites, píxeles calientes)
# antes de hacer la media; la "mediana" es aún más robusta pero algo más ruidosa.


def fuentes_para_alinear(imagen, n=100, saturacion=None):
    """Las estrellas más brillantes (no saturadas) de una foto con un centro FINO: centro de
    masas iterativo con peso gaussiano (como hace el centrado "windowed" de los programas
    profesionales), mucho más preciso que el centro rápido de la detección."""
    fuentes = [f for f in detectar_fuentes(imagen, 8.0, 3 * n, saturacion=saturacion) if not f['saturada']][:n]
    alto, ancho = imagen.shape
    r = 6
    off = np.arange(-r, r + 1)
    oy, ox = np.meshgrid(off, off, indexing='ij')
    anillo = (ox ** 2 + oy ** 2 >= 25) & (ox ** 2 + oy ** 2 <= 36)
    finas = []
    for f in fuentes:
        x, y = f['x'], f['y']
        s = max(0.6, min(f['fwhm'], 8.0) / 2.3548)
        for _ in range(10):
            cx, cy = int(round(x)), int(round(y))
            if cx - r < 0 or cy - r < 0 or cx + r >= ancho or cy + r >= alto:
                break
            v = imagen[cy - r:cy + r + 1, cx - r:cx + r + 1].astype(np.float64)
            if not np.isfinite(v).all():
                break
            luz = v - np.median(v[anillo])
            peso = luz * np.exp(-((cx + ox - x) ** 2 + (cy + oy - y) ** 2) / (2 * s * s))
            total = peso.sum()
            if total <= 0:
                break
            nx, ny = cx + (peso * ox).sum() / total, cy + (peso * oy).sum() / total
            mov = math.hypot(nx - x, ny - y)
            x, y = float(nx), float(ny)
            if mov < 1e-3 or mov > 3:
                break
        if abs(x - f['x']) < 2 and abs(y - f['y']) < 2:
            finas.append(dict(f, x=x, y=y))
    return finas


def _semejanza(origen, destino, espejo=False):
    """Mejor giro + escala + desplazamiento (y espejo, si se pide) que lleva unos puntos a
    otros, por mínimos cuadrados: 4 números en vez de los 6 de una afín, porque las dos
    fotos son de la misma cámara (los píxeles no se deforman). Devuelve (A, t)."""
    x, y = origen[:, 0], origen[:, 1]
    s = -1.0 if espejo else 1.0
    # sin espejo: X = a·x − b·y + tx,  Y = b·x + a·y + ty
    # con espejo: X = a·x + b·y + tx,  Y = b·x − a·y + ty
    uno, cero = np.ones_like(x), np.zeros_like(x)
    m = np.vstack([np.column_stack([x, -s * y, uno, cero]), np.column_stack([s * y, x, cero, uno])])
    sol, *_ = np.linalg.lstsq(m, np.concatenate([destino[:, 0], destino[:, 1]]), rcond=None)
    a, b, tx, ty = sol
    A = np.array([[a, b], [b, -a]]) if espejo else np.array([[a, -b], [b, a]])
    return A, np.array([tx, ty])


def alinear_fuentes(referencia, otra, n=25, tolerancia=2.0, tiempo_max=20.0):
    """Transformación que lleva las estrellas de 'otra' foto encima de las de la
    'referencia': (x_ref, y_ref) = A·(x, y) + t.

    referencia, otra: estrellas detectadas (detectar_fuentes), de más a menos brillante.
    Devuelve {'A', 't', 'n' (estrellas emparejadas), 'rms' (px), 'giro' (°), 'escala',
    'espejo'} o None si no se encuentra."""
    import time
    inicio = time.perf_counter()
    ref_todas = np.array([(f['x'], f['y']) for f in referencia[:4 * n]], dtype=np.float64)
    otra_todas = np.array([(f['x'], f['y']) for f in otra[:4 * n]], dtype=np.float64)
    if len(ref_todas) < 3 or len(otra_todas) < 3:
        return None
    inv_r, ver_r = _triangulos(ref_todas[:n], lado_minimo=10.0)
    inv_o, ver_o = _triangulos(otra_todas[:n], lado_minimo=10.0)
    if len(inv_r) == 0 or len(inv_o) == 0:
        return None
    orden = np.argsort(inv_r[:, 0])
    inv_r, ver_r = inv_r[orden], ver_r[orden]
    objetivo = min(len(ref_todas), len(otra_todas))

    def votos_de(A, t):
        movidas = otra_todas @ A.T + t
        d = np.hypot(movidas[:, None, 0] - ref_todas[None, :, 0], movidas[:, None, 1] - ref_todas[None, :, 1])
        return int((d.min(axis=1) < tolerancia).sum())

    mejor, mejor_votos, probadas = None, 0, 0
    tol = 0.006
    for k in np.argsort(-np.abs(np.diff(inv_o, axis=1)[:, 0])):        # primero los triángulos más "irregulares"
        a, b = inv_o[k]
        i0, i1 = np.searchsorted(inv_r[:, 0], [a - tol, a + tol])
        for j in np.nonzero(np.abs(inv_r[i0:i1, 1] - b) < tol)[0] + i0:
            A, t = _afin_3(otra_todas[ver_o[k]], ref_todas[ver_r[j]])
            if A is None:
                continue
            sv = np.linalg.svd(A, compute_uv=False)
            if sv[1] <= 0 or sv[0] / sv[1] > 1.05 or not 0.5 < math.sqrt(abs(np.linalg.det(A))) < 2.0:
                continue
            probadas += 1
            votos = votos_de(A, t)
            if votos > mejor_votos:
                mejor, mejor_votos = (A, t), votos
        if mejor_votos >= max(6, 0.6 * objetivo) or probadas > 20000 or time.perf_counter() - inicio > tiempo_max:
            break
    if mejor is None or mejor_votos < max(4, min(8, int(0.25 * objetivo))):
        return None
    A, t = mejor
    # Afinar: mínimos cuadrados con todas las parejas, quitando las que no encajan
    for _ in range(4):
        movidas = otra_todas @ A.T + t
        d = np.hypot(movidas[:, None, 0] - ref_todas[None, :, 0], movidas[:, None, 1] - ref_todas[None, :, 1])
        cercana = d.argmin(axis=1)
        dist = d[np.arange(len(otra_todas)), cercana]
        usar = dist < tolerancia
        if usar.sum() >= 6:
            rms_actual = math.sqrt(float(np.mean(dist[usar] ** 2)))
            usar &= dist < max(3.0 * rms_actual, 0.15)
        if usar.sum() < 3:
            break
        A, t = _semejanza(otra_todas[usar], ref_todas[cercana[usar]], espejo=np.linalg.det(A) < 0)
    movidas = otra_todas[usar] @ A.T + t
    residuos = np.hypot(*(movidas - ref_todas[cercana[usar]]).T)
    rms = math.sqrt(float(np.mean(residuos ** 2))) if residuos.size else float('nan')
    det = float(np.linalg.det(A))
    if det >= 0:      # giro: A ≈ s·[[cos, −sin], [sin, cos]]
        giro = math.atan2(A[1, 0] - A[0, 1], A[0, 0] + A[1, 1])
    else:             # con espejo: A ≈ s·[[cos, sin], [sin, −cos]]
        giro = math.atan2(A[1, 0] + A[0, 1], A[0, 0] - A[1, 1])
    return {'A': A, 't': np.asarray(t, dtype=np.float64), 'n': int(usar.sum()), 'rms': rms,
            'giro': math.degrees(giro) % 360, 'escala': math.sqrt(abs(det)), 'espejo': det < 0}


def remuestrear(imagen, A, t, forma, filas=None):
    """La foto 'imagen' llevada a la rejilla de la referencia (forma = (alto, ancho) de la
    referencia) con la transformación (A, t) de alinear_fuentes. Interpolación bilineal;
    lo que cae fuera de la foto queda en NaN. 'filas' = (desde, hasta) para hacerlo a trozos."""
    alto, ancho = forma
    f0, f1 = filas if filas else (0, alto)
    inversa = np.linalg.inv(np.asarray(A, dtype=np.float64))
    # Las coordenadas son lineales: x = a·columna + b·fila + c (se suman una fila y una columna)
    columnas = np.arange(ancho, dtype=np.float64) - t[0]
    filas_v = np.arange(f0, f1, dtype=np.float64) - t[1]
    x = (inversa[0, 0] * columnas)[None, :].astype(np.float32) + (inversa[0, 1] * filas_v)[:, None].astype(np.float32)
    y = (inversa[1, 0] * columnas)[None, :].astype(np.float32) + (inversa[1, 1] * filas_v)[:, None].astype(np.float32)
    alto_o, ancho_o = imagen.shape
    fuera = (x < 0) | (y < 0) | (x > ancho_o - 1) | (y > alto_o - 1)
    x0 = np.clip(x.astype(np.int32), 0, max(ancho_o - 2, 0))    # (dentro, x ≥ 0: truncar = suelo)
    y0 = np.clip(y.astype(np.int32), 0, max(alto_o - 2, 0))
    wx = np.clip(x - x0, 0, 1)
    wy = np.clip(y - y0, 0, 1)
    del x, y
    plana = np.ascontiguousarray(imagen, dtype=np.float32).ravel()
    i00 = y0 * ancho_o + x0
    del x0, y0
    paso_x = 1 if ancho_o > 1 else 0
    paso_y = ancho_o if alto_o > 1 else 0
    arriba = plana.take(i00) * (1 - wx) + plana.take(i00 + paso_x) * wx
    abajo = plana.take(i00 + paso_y) * (1 - wx) + plana.take(i00 + paso_y + paso_x) * wx
    salida = arriba * (1 - wy) + abajo * wy
    salida[fuera] = np.nan
    return salida


def apilar(imagenes, transformaciones, forma, metodo='media', recorte_sigma_px=3.0, progreso=None, banda=200):
    """Apila fotos ya emparejadas con la referencia.

    imagenes: lista de fotos (la referencia incluida); transformaciones: lista paralela de
    (A, t) que las lleva a la referencia (la de la referencia, la identidad).
    metodo: 'media' (media con recorte de los valores raros, a 3σ) o 'mediana'.
    Devuelve (apilada, cuantas): la foto apilada (float32, con el fondo de la referencia)
    y cuántas fotos han contribuido a cada píxel."""
    import warnings
    alto, ancho = forma
    n = len(imagenes)
    fondos, ruidos = [], []
    for im in imagenes:
        validos = im[np.isfinite(im)]
        f, r = fondo_y_ruido(validos[None, :]) if validos.size else (0.0, 0.0)
        fondos.append(f)
        ruidos.append(r)
    ruido_cielo = float(np.median([r for r in ruidos if r > 0])) if any(r > 0 for r in ruidos) else 0.0
    apilada = np.empty((alto, ancho), dtype=np.float32)
    cuantas = np.zeros((alto, ancho), dtype=np.int16)
    for f0 in range(0, alto, banda):
        f1 = min(f0 + banda, alto)
        pila = np.stack([remuestrear(im, A, t, forma, (f0, f1)) - np.float32(fondo - fondos[0])
                         for im, (A, t), fondo in zip(imagenes, transformaciones, fondos)])
        validos = np.isfinite(pila)
        cuantas[f0:f1] = validos.sum(axis=0)
        # (np.median es mucho más rápida que np.nanmedian: se usa si en este trozo no hay huecos)
        mediana_de = np.median if validos.all() else np.nanmedian
        with np.errstate(all='ignore'), warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)       # píxeles sin ninguna foto
            mediana = mediana_de(pila, axis=0)
            if metodo == 'mediana':
                resultado = mediana
            elif n < 3:
                resultado = np.nanmean(pila, axis=0)
            else:
                # σ de cada píxel: la dispersión entre fotos (MAD), y como mínimo el ruido
                # del cielo; con pocas fotos la MAD es poco fiable y se deja más margen.
                mad = 1.4826 * mediana_de(np.abs(pila - mediana[None]), axis=0)
                sigma = np.maximum(mad, ruido_cielo) * (1.0 if n >= 5 else 1.5)
                pila[np.abs(pila - mediana[None]) > recorte_sigma_px * sigma[None]] = np.nan
                resultado = np.nanmean(pila, axis=0)
        apilada[f0:f1] = resultado
        if progreso:
            progreso(f1 / alto)
    apilada[cuantas == 0] = np.nan
    return apilada, cuantas
