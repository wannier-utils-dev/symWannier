#!/usr/bin/env python3
"""Build an external Ni projector file from the PP_CHI arrays in the UPF."""

from __future__ import annotations

import math
import os
from pathlib import Path
import re


CASE_DIR = Path(__file__).resolve().parent
UPF_PATH = Path(os.environ.get(
    "NI_UPF", str(CASE_DIR / "Ni.pbe-n-rrkjus_psl.0.1.UPF")
))
OUTPUT_PATH = CASE_DIR / "atom_proj" / "Ni.dat"


upf_text = UPF_PATH.read_text(encoding="utf-8")


def extract_tag(name: str) -> tuple[str, list[float]]:
    escaped_name = re.escape(name)
    match = re.search(
        rf"<{escaped_name}\b([^>]*)>(.*?)</{escaped_name}>",
        upf_text,
        flags=re.DOTALL,
    )
    if match is None:
        raise RuntimeError(f"{name} not found in {UPF_PATH}")
    attributes, body = match.groups()
    return attributes, [float(value) for value in body.split()]


_, rgrid = extract_tag("PP_R")
channels: list[tuple[int, list[float]]] = []
for index in range(1, 4):
    attributes, radial = extract_tag(f"PP_CHI.{index}")
    l_match = re.search(r'\bl="([0-9]+)"', attributes)
    if l_match is None:
        raise RuntimeError(f"l attribute missing from PP_CHI.{index}")
    channels.append((int(l_match.group(1)), radial))

if [channel[0] for channel in channels] != [0, 1, 2]:
    raise RuntimeError("expected Ni PP_CHI channels in s, p, d order")
if any(len(radial) != len(rgrid) for _, radial in channels):
    raise RuntimeError("PP_R and PP_CHI mesh sizes do not match")

# Duplicate each radial channel. The second s/p/d copies are excluded in
# pw2wan.in, so this sample exercises atom_proj_ext + atom_proj_exclude while
# retaining the same nine s/p/d orbitals as the ordinary Ni UPF test.
expanded_l = [0, 0, 1, 1, 2, 2]
expanded_radial = [
    channels[0][1],
    channels[0][1],
    channels[1][1],
    channels[1][1],
    channels[2][1],
    channels[2][1],
]

OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
with OUTPUT_PATH.open("w", encoding="utf-8") as output:
    output.write("# Generated from Ni PP_R and PP_CHI.1-3 by make_atom_proj.py\n")
    output.write(f"{len(rgrid)} {len(expanded_l)}\n")
    output.write(" ".join(str(l_value) for l_value in expanded_l) + "\n")
    for row, radius in enumerate(rgrid):
        values = [math.log(radius), radius]
        values.extend(radial[row] for radial in expanded_radial)
        output.write(" ".join(f"{value:.16e}" for value in values) + "\n")

print(OUTPUT_PATH)
