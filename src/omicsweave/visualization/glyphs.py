"""Packaged schematic SVGs for standalone reuse and vector Matplotlib rendering."""
from __future__ import annotations

import json
from functools import lru_cache
from importlib.resources import files
from math import isfinite
from xml.etree import ElementTree

from matplotlib.axes import Axes
from matplotlib.patches import PathPatch
from matplotlib.path import Path
from matplotlib.transforms import Affine2D

_ASSETS = files("omicsweave.visualization").joinpath("assets", "glyphs")
_SVG_NS = "{http://www.w3.org/2000/svg}"
_VIEW_SIZE = 32
_COMMANDS = {
    "M": (Path.MOVETO, 1),
    "L": (Path.LINETO, 1),
    "Q": (Path.CURVE3, 2),
    "C": (Path.CURVE4, 3),
    "Z": (Path.CLOSEPOLY, 0),
}


@lru_cache(maxsize=1)
def _manifest() -> dict:
    return json.loads(_ASSETS.joinpath("manifest.json").read_text(encoding="utf-8"))


GLYPH_STYLES = tuple(_manifest()["styles"])


def glyph_ids(*, category: str | None = None) -> tuple[str, ...]:
    """List glyph IDs in manifest order, optionally filtering by semantic category."""
    return tuple(
        glyph_id for glyph_id, item in _manifest()["glyphs"].items()
        if category is None or item["category"] == category
    )


def read_glyph_svg(glyph_id: str, *, style: str = "nature_editorial") -> str:
    """Read a self-contained SVG by semantic ID from the installed package."""
    if style not in GLYPH_STYLES:
        raise ValueError(f"Unknown glyph style {style!r}. Available styles: {', '.join(GLYPH_STYLES)}")
    try:
        relative = _manifest()["glyphs"][glyph_id]["files"][style]
    except KeyError:
        raise ValueError(f"Unknown schematic glyph {glyph_id!r}") from None
    return _ASSETS.joinpath(*relative.split("/")).read_text(encoding="utf-8")


def _path(data: str) -> Path:
    """Decode the explicit absolute commands used by the bundled SVG paths."""
    tokens = iter(data.split())
    vertices = []
    codes = []
    for command in tokens:
        code, count = _COMMANDS[command]
        if command == "Z":
            vertices.append((0.0, 0.0))
            codes.append(code)
        else:
            for _ in range(count):
                vertices.append((float(next(tokens)), float(next(tokens))))
                codes.append(code)
    return Path(vertices, codes, readonly=True)


@lru_cache(maxsize=None)
def _layers(glyph_id: str, style: str) -> tuple[tuple[Path, str, str, float], ...]:
    root = ElementTree.fromstring(read_glyph_svg(glyph_id, style=style))
    return tuple(
        (_path(layer.attrib["d"]), layer.attrib["fill"], layer.get("stroke", "none"),
         float(layer.get("stroke-width", "0")))
        for layer in root.findall(f"{_SVG_NS}path")
    )


def draw_glyph(
    ax: Axes,
    glyph_id: str,
    xy: tuple[float, float],
    *,
    style: str = "nature_editorial",
    size: float = 4.0,
    zorder: float = 5.0,
) -> tuple[PathPatch, ...]:
    """Draw an SVG glyph as vector patches on linear schematic axes.

    `xy` is the centre and `size` is the side of its square view box in data units.
    The orientation follows the axes' current y direction. Set the axes limits
    and aspect before drawing so that stroke widths match the rendered size.
    """
    if not isfinite(size) or size <= 0:
        raise ValueError("Glyph size must be finite and positive")
    layers = _layers(glyph_id, style)
    scale = size / _VIEW_SIZE
    y_sign = 1 if ax.yaxis_inverted() else -1
    transform = (
        Affine2D().translate(-_VIEW_SIZE / 2, -_VIEW_SIZE / 2)
        .scale(scale, y_sign * scale).translate(*xy) + ax.transData
    )
    origin, unit = transform.transform(((0.0, 0.0), (1.0, 0.0)))
    points_per_unit = abs(unit[0] - origin[0]) * 72 / ax.figure.dpi
    patches = []
    for path, fill, stroke, width in layers:
        patch = PathPatch(
            path, transform=transform, facecolor=fill, edgecolor=stroke,
            linewidth=width * points_per_unit, capstyle="round", joinstyle="round",
            antialiased=style != "pixel_abstract", zorder=zorder,
        )
        ax.add_patch(patch)
        patches.append(patch)
    return tuple(patches)
