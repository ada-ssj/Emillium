# -*- mode: python ; coding: utf-8 -*-
# Emillium Expanded 1.5: numpy, tkinter y matplotlib (solo para los gráficos).
# Los módulos ciencia.py, ades.py y catalogos.py se incluyen solos (Emillium.py los importa).

a = Analysis(
    ['Emillium.py'],
    pathex=[],
    binaries=[],
    datas=[('icono.ico', '.')],
    hiddenimports=['matplotlib.backends.backend_tkagg'],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=['astropy', 'photutils', 'scipy', 'pandas', 'IPython', 'jupyter', 'notebook',
              'matplotlib.backends.backend_qt5agg', 'matplotlib.backends.backend_qtagg', 'PyQt5', 'PyQt6',
              'PySide2', 'PySide6', 'wx', 'gi'],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='EmilliumExpanded',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=['icono.ico'],
)
