# -*- mode: python ; coding: utf-8 -*-
# Emillium v1.0: solo necesita numpy y tkinter.

a = Analysis(
    ['Emillium.py'],
    pathex=[],
    binaries=[],
    datas=[('icono.ico', '.')],
    hiddenimports=[],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=['astropy', 'photutils', 'matplotlib', 'scipy', 'PIL', 'pandas'],
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
    name='Emillium',
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
