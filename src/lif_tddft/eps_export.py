"""Checked EPS export for Matplotlib manuscript figures.

Type 3 fonts contain vector glyph outlines; they do not rasterize a figure.
Some Type 42 export stacks emit glyphshow names absent from the embedded font.
Checking the process exit code alone does not detect the resulting blank text.
"""

from __future__ import annotations

from pathlib import Path
import re
import shutil
import subprocess
import tempfile


class EPSFontError(ValueError):
    """An embedded-font glyph reference cannot be resolved."""


_FONT = re.compile(
    r"/FontName /(?P<name>[^\s]+) def.*?"
    r"/CharStrings \d+ dict dup begin(?P<glyphs>.*?)end readonly def",
    re.DOTALL,
)
_GLYPH_DEFINITION = re.compile(r"/([^\s/{}]+)(?=\s+\d+\s+def|\s*\{)")
_DRAW = re.compile(
    r"/(?P<font>[^\s]+)\s+[-+.0-9]+\s+selectfont"
    r"|/(?P<glyph>[^\s]+)\s+glyphshow"
)


def validate_eps_fonts(path: str | Path) -> dict:
    """Check named glyph references in Matplotlib Type 3/42 EPS output.

    This is a check of Matplotlib's embedded-font/glyphshow convention, not a
    general PostScript interpreter. Raster-only EPS files have zero references.
    Unknown fonts and undefined glyphs fail explicitly. Syntax and visual
    correctness still require rendering. The optional DSC %%EOF is not required.
    """
    path = Path(path)
    text = path.read_bytes().decode("latin1")
    if not text.startswith("%!PS-Adobe-") or "EPSF-" not in text.splitlines()[0]:
        raise EPSFontError(f"Not an EPS file: {path}")
    fonts = {
        match["name"]: set(_GLYPH_DEFINITION.findall(match["glyphs"]))
        for match in _FONT.finditer(text)
    }
    active_font = None
    count = 0
    failures = set()
    for match in _DRAW.finditer(text):
        if match["font"]:
            active_font = match["font"]
        else:
            count += 1
            glyph = match["glyph"]
            if glyph not in fonts.get(active_font, set()):
                failures.add((str(active_font), glyph))
    if failures:
        examples = ", ".join(f"{font}: /{glyph}" for font, glyph in sorted(failures)[:8])
        raise EPSFontError(f"{path.name}: unresolved glyph references ({examples})")
    return {"file": str(path), "embedded_fonts": sorted(fonts), "glyph_references": count}


def save_eps(figure, path: str | Path, *, dpi: int = 600) -> dict:
    """Export with vector glyph outlines and validate before replacing a file."""
    import matplotlib

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, suffix=".eps", delete=False) as stream:
        temporary = Path(stream.name)
    try:
        with matplotlib.rc_context({"ps.fonttype": 3, "ps.useafm": False}):
            figure.savefig(temporary, format="eps", dpi=dpi)
        audit = validate_eps_fonts(temporary)
        temporary.replace(path)
        audit["file"] = str(path)
        return audit
    finally:
        temporary.unlink(missing_ok=True)


def find_ghostscript(executable: str | Path | None = None) -> str:
    """Resolve an explicit executable or an installed Ghostscript command."""
    if executable:
        found = shutil.which(str(executable))
        if not found:
            raise FileNotFoundError(f"Ghostscript executable not found: {executable}")
        return found
    for name in ("gs", "gswin64c", "gswin32c"):
        found = shutil.which(name)
        if found:
            return found
    raise FileNotFoundError("Install Ghostscript or provide --ghostscript /path/to/executable")


def render_eps_for_review(
    path: str | Path,
    output_dir: str | Path,
    *,
    ghostscript: str | Path | None = None,
    dpi: int = 200,
) -> dict:
    """Convert the checked EPS to PDF and render that PDF, not the source plot."""
    audit = validate_eps_fonts(path)
    executable = find_ghostscript(ghostscript)
    path = Path(path).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    pdf = output_dir / f"{path.stem}_from_eps.pdf"
    png = output_dir / f"{path.stem}_from_eps.png"
    common = [executable, "-dSAFER", "-dBATCH", "-dNOPAUSE"]
    commands = [
        common + ["-dEPSCrop", "-sDEVICE=pdfwrite", "-dAutoRotatePages=/None",
                  "-dEmbedAllFonts=true", "-dDownsampleColorImages=false",
                  "-dDownsampleGrayImages=false", f"-sOutputFile={pdf}", str(path)],
        common + ["-sDEVICE=png16m", "-dTextAlphaBits=4", "-dGraphicsAlphaBits=4",
                  f"-r{dpi}", f"-sOutputFile={png}", str(pdf)],
    ]
    logs = []
    for command in commands:
        result = subprocess.run(command, capture_output=True, text=True, timeout=120)
        logs.append(result.stdout + result.stderr)
        if result.returncode:
            raise RuntimeError(f"Ghostscript failed for {path.name}:\n{logs[-1]}")
    if not pdf.is_file() or not png.is_file() or not pdf.stat().st_size or not png.stat().st_size:
        raise RuntimeError(f"Missing PDF/PNG review output for {path.name}")
    (output_dir / f"{path.stem}_render.log").write_text("\n".join(logs), encoding="utf-8")
    return {**audit, "review_pdf": str(pdf), "review_png": str(png), "render_dpi": dpi}
