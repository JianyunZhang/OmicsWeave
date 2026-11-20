from pathlib import Path
from xml.etree import ElementTree

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest

from omicsweave.visualization import GLYPH_STYLES, draw_glyph, glyph_ids, read_glyph_svg
from omicsweave.visualization.graphical_abstract import load_abstract_spec, render_graphical_abstract

_NS = {"svg": "http://www.w3.org/2000/svg"}


@pytest.mark.parametrize("style", GLYPH_STYLES)
def test_packaged_glyphs_draw_as_vectors(style):
    fig, ax = plt.subplots()
    try:
        for gid in glyph_ids():
            root = ElementTree.fromstring(read_glyph_svg(gid, style=style))
            assert root.attrib["viewBox"] == "0 0 32 32"
            assert root.findtext("svg:title", namespaces=_NS) == gid
            assert root.findtext("svg:desc", namespaces=_NS)
            assert {child.tag.split("}")[-1] for child in root} == {"title", "desc", "path"}
            patches = draw_glyph(ax, gid, (0.5, 0.5), size=0.2, style=style)
            assert patches
            if style == "pixel_abstract":
                for patch in patches:
                    assert not patch.get_antialiased()
                    assert (patch.get_path().vertices % 2 == 0).all()
    finally:
        plt.close(fig)


def test_pixel_graphical_abstract_svg_has_no_raster_images(tmp_path):
    example = Path(__file__).resolve().parents[2] / "examples" / "graphical_abstract.example.yaml"
    (path,) = render_graphical_abstract(
        load_abstract_spec(example), tmp_path / "pixel", style="pixel_abstract", formats=("svg",),
    )
    root = ElementTree.parse(path)
    assert root.find(".//svg:path", _NS) is not None
    assert root.find(".//svg:image", _NS) is None
    assert root.find(".//svg:text", _NS) is not None
