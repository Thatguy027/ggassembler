# PyInstaller spec for the double-clickable app.
#
# Build it from the repository root:
#
#     pyinstaller packaging/GGAssembler.spec --noconfirm
#
# One spec for every platform. macOS gets a .app from the BUNDLE block at the
# bottom; on Windows and Linux that block is skipped and COLLECT is the answer.
# Kept as a spec rather than a pile of command-line flags because the data
# files and the hidden imports below are not guesses - each one is something
# the analysis cannot see, and a flag list that long belongs in a file.

from pathlib import Path

from PyInstaller.utils.hooks import collect_data_files, collect_submodules

SPEC_DIR = Path(SPECPATH)
ROOT = SPEC_DIR.parent

datas = [
    # The whole frontend. ``api.main`` finds it with ``__file__``, which under
    # PyInstaller points into the unpacked bundle, so the layout has to match
    # the package: ggassembler/web/..., not web/... at the top level.
    (str(ROOT / "ggassembler" / "web"), "ggassembler/web"),
]
datas += collect_data_files("Bio")

hiddenimports = [
    # Imported by name inside uvicorn rather than at the top of a module, so
    # the analysis never sees them.
    "uvicorn.logging",
    "uvicorn.loops.auto",
    "uvicorn.loops.asyncio",
    "uvicorn.protocols.http.auto",
    "uvicorn.protocols.http.h11_impl",
    "uvicorn.protocols.websockets.auto",
    "uvicorn.lifespan.on",
]
hiddenimports += collect_submodules("ggassembler")

a = Analysis(
    [str(SPEC_DIR / "launch.py")],
    pathex=[str(ROOT)],
    binaries=[],
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    runtime_hooks=[],
    # Nothing in the app draws a plot or opens a notebook, and these arrive as
    # transitive weight if they happen to be in the build environment.
    excludes=["tkinter", "matplotlib", "numpy.testing", "pytest", "IPython"],
    noarchive=False,
)

pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="GG Assembler",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=False,
)

coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=False,
    name="GG Assembler",
)

icon = ROOT / "packaging" / "icon.icns"

app = BUNDLE(
    coll,
    name="GG Assembler.app",
    icon=str(icon) if icon.exists() else None,
    bundle_identifier="edu.ucla.ggassembler",
    info_plist={
        "CFBundleName": "GG Assembler",
        "CFBundleDisplayName": "GG Assembler",
        "CFBundleShortVersionString": "0.1.1",
        "CFBundleVersion": "0.1.1",
        "NSHighResolutionCapable": True,
        # The app is one window; closing it should quit rather than leave a
        # menu bar behind with no way back to a window.
        "LSApplicationCategoryType": "public.app-category.education",
        "NSHumanReadableCopyright": "MIT",
    },
)
