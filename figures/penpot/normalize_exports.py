#!/usr/bin/env python3
"""Make Penpot Figure 1 exports portable and publication-sized.

Penpot exports at CSS-pixel dimensions and references its self-hosted Inter
font endpoint.  This script embeds those font files, sets a physical width of
180 mm, and renders deterministic PDF/PNG derivatives.  It changes no layout
or scientific content; the native Penpot file remains the design source.
"""

from __future__ import annotations

import base64
import re
from pathlib import Path
from urllib.request import urlopen

import cairosvg

ROOT = Path(__file__).resolve().parents[2]
EXPORTS = ROOT / "figures" / "penpot" / "exports"
WIDTH_MM = 180.0
VIEWBOX_WIDTH = 1800.0
VIEWBOX_HEIGHT = 1120.0
HEIGHT_MM = WIDTH_MM * VIEWBOX_HEIGHT / VIEWBOX_WIDTH
FONT_URL = re.compile(r"url\((http://localhost:9001/internal/gfonts/[^)]+\.woff2)\)")


def normalize_svg(path: Path) -> None:
    svg = path.read_text(encoding="utf-8")
    svg = re.sub(r'<svg width="[^"]+"', f'<svg width="{WIDTH_MM:g}mm"', svg, count=1)
    svg = re.sub(r'height="[^"]+"', f'height="{HEIGHT_MM:g}mm"', svg, count=1)

    cache: dict[str, str] = {}

    def embed_font(match: re.Match[str]) -> str:
        url = match.group(1)
        if url not in cache:
            with urlopen(url, timeout=30) as response:
                payload = response.read()
            cache[url] = base64.b64encode(payload).decode("ascii")
        return f"url(data:font/woff2;base64,{cache[url]})"

    svg = FONT_URL.sub(embed_font, svg)
    if "http://localhost:9001" in svg:
        raise RuntimeError(f"Unresolved local Penpot URL in {path}")
    path.write_text(svg, encoding="utf-8")

    cairosvg.svg2pdf(
        bytestring=svg.encode("utf-8"), write_to=str(path.with_suffix(".pdf"))
    )
    cairosvg.svg2png(
        bytestring=svg.encode("utf-8"),
        write_to=str(path.with_suffix(".png")),
        output_width=3600,
        output_height=2240,
    )


def main() -> None:
    paths = sorted(EXPORTS.glob("fig1_variant_*.svg"))
    if not paths:
        raise SystemExit(f"No Penpot SVG exports found in {EXPORTS}")
    for path in paths:
        normalize_svg(path)
        print(path.relative_to(ROOT))


if __name__ == "__main__":
    main()
