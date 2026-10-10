"""Put generated conceptual cutouts on an exact white publication background.

The Azure raster may contain a near-white paper tint despite a white prompt.
Only high-value, nearly neutral pixels are reset; colored illustration pixels
remain unchanged. Scientific labels and geometry are added later as native SVG.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
ASSET = ROOT / "paper/figures/journal/teaser_components"
NAMES = ("solvent_halo", "endpoint_lens")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(name: str) -> None:
    source = ASSET / f"{name}_v1.png"
    output = ASSET / f"{name}_white.png"
    pixels = np.asarray(Image.open(source).convert("RGB")).copy()
    hi = pixels.max(axis=2).astype(np.int16)
    lo = pixels.min(axis=2).astype(np.int16)
    background = (lo >= 238) & (hi - lo <= 12)
    pixels[background] = (255, 255, 255)
    Image.fromarray(pixels, "RGB").save(output, optimize=True)

    height, width = pixels.shape[:2]
    border = np.concatenate((pixels[:24].reshape(-1, 3),
                             pixels[-24:].reshape(-1, 3),
                             pixels[:, :24].reshape(-1, 3),
                             pixels[:, -24:].reshape(-1, 3)))
    assert np.all(border == 255), f"{name}: background edge is not pure white"
    center = pixels[height // 2 - 50:height // 2 + 50,
                    width // 2 - 50:width // 2 + 50]
    assert np.all(center == 255), f"{name}: central opening is not pure white"

    record = {
        "role": "conceptual_illustration_not_scientific_data",
        "source": str(source.relative_to(ROOT)),
        "source_sha256": sha256(source),
        "output": str(output.relative_to(ROOT)),
        "output_sha256": sha256(output),
        "background_rgb": [255, 255, 255],
        "background_rule": "RGB minimum >= 238 and channel range <= 12",
        "reset_pixels": int(background.sum()),
        "total_pixels": int(width * height),
        "verified_border_px": 24,
        "verified_center_px": 100,
    }
    output.with_suffix(".provenance.json").write_text(
        json.dumps(record, indent=2) + "\n")
    print(f"{output.relative_to(ROOT)}: white border and opening verified")


if __name__ == "__main__":
    for component in NAMES:
        prepare(component)
