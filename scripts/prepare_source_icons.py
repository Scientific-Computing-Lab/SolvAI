#!/usr/bin/env python3
"""Verify and crop the white-background conceptual source-icon atlas.

The atlas is a non-evidential illustration already processed to have a pure
#FFFFFF background. Source labels, connections and counts remain native SVG.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "paper/figures/journal/teaser_components"
SOURCE = ASSETS / "source_icon_atlas.png"


def main():
    atlas = Image.open(SOURCE).convert("RGB")
    if atlas.size != (1536, 1024):
        raise ValueError(f"Expected a 1536x1024 atlas, got {atlas.size}")
    rgb = np.asarray(atlas, dtype=np.uint8)
    edge = np.concatenate([rgb[0], rgb[-1], rgb[:, 0], rgb[:, -1]], axis=0)
    if not np.all(edge == 255):
        raise AssertionError("Source atlas has a non-white unillustrated edge")
    outputs = []
    for index in range(6):
        col, row = index % 3, index // 3
        crop = atlas.crop((col * 512, row * 512, (col + 1) * 512, (row + 1) * 512))
        path = ASSETS / f"source_icon_{index+1}.png"
        crop.save(path, optimize=True)
        corners = [crop.getpixel(p) for p in [(0, 0), (511, 0), (0, 511), (511, 511)]]
        if any(pixel != (255, 255, 255) for pixel in corners):
            raise AssertionError(f"Cell {index+1} has a nonwhite corner: {corners}")
        outputs.append({"path": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                        "corner_rgb": [255, 255, 255]})
    record = {
        "role": "conceptual source-category icons; not molecular structures, model outputs or data",
        "source_atlas": SOURCE.name,
        "source_atlas_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "background": "unillustrated background and cell edges #FFFFFF",
        "outputs": outputs,
    }
    (ASSETS / "source_icon_processing.json").write_text(json.dumps(record, indent=2) + "\n")


if __name__ == "__main__":
    main()
