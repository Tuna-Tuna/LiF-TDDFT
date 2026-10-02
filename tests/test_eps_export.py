import os
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from lif_tddft.eps_export import (
    EPSFontError, find_ghostscript, render_eps_for_review, save_eps,
    validate_eps_fonts,
)

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from PIL import Image
    HAS_FIGURES = True
except ImportError:
    HAS_FIGURES = False


class EPSGlyphTests(unittest.TestCase):
    def check_text(self, text):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "test.eps"
            path.write_text(text, encoding="ascii")
            return validate_eps_fonts(path)

    def font(self, name, definitions):
        return (f"/FontName /{name} def\n/CharStrings 2 dict dup begin\n"
                f"{definitions}\nend readonly def\n")

    def test_type42_defined_glyph(self):
        text = "%!PS-Adobe-3.0 EPSF-3.0\n" + self.font("Times", "/two 12 def")
        result = self.check_text(text + "/Times 12 selectfont\n/two glyphshow\nshowpage\n")
        self.assertEqual(result["glyph_references"], 1)

    def test_type42_mismatched_original_glyph_id_is_rejected(self):
        text = "%!PS-Adobe-3.0 EPSF-3.0\n" + self.font("Times", "/two 12 def")
        with self.assertRaisesRegex(EPSFontError, "uni00000015"):
            self.check_text(text + "/Times 12 selectfont\n/uni00000015 glyphshow\n%%EOF")

    def test_glyph_must_exist_in_active_font(self):
        text = "%!PS-Adobe-3.0 EPSF-3.0\n"
        text += self.font("Roman", "/two 12 def") + self.font("Italic", "/v 9 def")
        with self.assertRaisesRegex(EPSFontError, "Italic: /two"):
            self.check_text(text + "/Italic 12 selectfont\n/two glyphshow\n%%EOF")

    def test_type3_outline_names_without_optional_eof(self):
        text = "%!PS-Adobe-3.0 EPSF-3.0\n"
        text += self.font("Outline", "/uni00000015{100 0 0 0 100 100 setcachedevice} def")
        result = self.check_text(text + "/Outline 12 selectfont\n/uni00000015 glyphshow\nshowpage")
        self.assertEqual(result["glyph_references"], 1)

    def test_raster_only_and_non_eps(self):
        result = self.check_text("%!PS-Adobe-3.0 EPSF-3.0\nshowpage\n%%EOF")
        self.assertEqual(result["glyph_references"], 0)
        with self.assertRaises(EPSFontError):
            self.check_text("not an EPS")


@unittest.skipUnless(HAS_FIGURES, "Install the figures extra")
class EPSRenderTests(unittest.TestCase):
    def make_figure(self):
        fig = plt.figure(figsize=(4, 2))
        fig.text(.1, .75, "0123456789 regular", family="DejaVu Serif", fontsize=20)
        fig.text(.1, .5, "Bold glyphs", family="DejaVu Serif", weight="bold", fontsize=20)
        fig.text(.1, .25, "Italic glyphs", family="DejaVu Serif", style="italic", fontsize=20)
        return fig

    def test_export_uses_valid_outlines_and_restores_settings(self):
        with tempfile.TemporaryDirectory() as folder:
            fig = self.make_figure()
            try:
                with matplotlib.rc_context({"ps.fonttype": 42}):
                    path = Path(folder) / "figure.eps"
                    audit = save_eps(fig, path)
                    self.assertEqual(matplotlib.rcParams["ps.fonttype"], 42)
                self.assertGreater(audit["glyph_references"], 20)
                self.assertGreaterEqual(len(audit["embedded_fonts"]), 3)
                self.assertIn("/FontType 3 def", path.read_text(encoding="latin1"))
            finally:
                plt.close(fig)

    def test_text_survives_eps_to_pdf_to_png(self):
        try:
            gs = find_ghostscript(os.environ.get("GHOSTSCRIPT"))
        except FileNotFoundError:
            self.skipTest("Ghostscript is not installed")
        with tempfile.TemporaryDirectory() as folder:
            fig = self.make_figure()
            try:
                path = Path(folder) / "figure.eps"
                save_eps(fig, path)
                review = render_eps_for_review(path, Path(folder) / "review", ghostscript=gs)
                with Image.open(review["review_png"]) as im:
                    gray = np.asarray(im.convert("L"))
                # This figure contains only text: each independent line must
                # leave visible ink after conversion, including bold/italic.
                height = gray.shape[0]
                for lo, hi in ((.1, .3), (.35, .55), (.6, .8)):
                    band = gray[int(height * lo):int(height * hi)]
                    self.assertGreater(int((band < 128).sum()), 1000)
            finally:
                plt.close(fig)


if __name__ == "__main__":
    unittest.main()
