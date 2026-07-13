"""Regenerate the README telemetry plot from the bundled historical sample.

This script intentionally does not import the unfinished historical decoder. It
implements only the scale factors documented in ``tlm_transcription2.py`` and
accepts the ``%DATE@TIME;`` sample rows that file was written to decode.
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "tlm.3"
OUTPUT = ROOT / "docs" / "assets" / "tlm-sensor-timeseries.png"
ROW = re.compile(r"%(?P<date>\d{8})@(?P<time>\d{6});(?P<payload>[0-9a-fA-F]{24})")
COLORS = ["#56A3FF", "#35C759", "#FF4D4D", "#F2A900", "#A970FF", "#34C6D9"]


def decode_rows(text: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    payloads = [bytes.fromhex(match.group("payload")) for match in ROW.finditer(text)]
    if not payloads:
        raise ValueError(f"No %DATE@TIME telemetry rows found in {SOURCE}")
    values = np.asarray([list(payload) for payload in payloads], dtype=float)
    gyro = (values[:, 0:3] - 2**7) * 0.14
    magnetometer = (values[:, 3:6] - 2**7) * 0.29
    solar = values[:, 6:12] * 12.89e-3
    return gyro, magnetometer, solar


def style_axis(axis: Axes) -> None:
    axis.set_facecolor("#171C24")
    axis.grid(True, color="#303641", linewidth=0.8, alpha=0.8)
    axis.tick_params(colors="#CBD3DF")
    for spine in axis.spines.values():
        spine.set_color("#3A414C")


def main() -> None:
    gyro, magnetometer, solar = decode_rows(SOURCE.read_text(encoding="utf-8", errors="replace"))
    sample_index = np.arange(len(gyro))

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})
    figure, axes = plt.subplots(3, 1, figsize=(12.8, 8.53), sharex=True)
    figure.patch.set_facecolor("#0B0F14")
    figure.suptitle(
        "PSC CubeSat telemetry decoder sample — didactic only",
        color="#EEF2F8",
        fontsize=20,
        fontweight="bold",
        y=0.985,
    )
    figure.text(
        0.5,
        0.946,
        f"Decoded {len(gyro)} sample rows from tlm.3 using tlm_transcription2.py scale factors; not flight validation evidence.",
        ha="center",
        color="#98A1AF",
        fontsize=11,
    )

    for axis in axes:
        style_axis(axis)

    for column, label, color in zip(range(3), ("gyro x", "gyro y", "gyro z"), COLORS[:3], strict=True):
        axes[0].plot(sample_index, gyro[:, column], color=color, linewidth=1.8, label=label)
    axes[0].set_ylabel("deg/s", color="#CBD3DF")
    axes[0].legend(loc="upper right", ncol=3, frameon=True, facecolor="#0B0F14", edgecolor="#3A414C", labelcolor="#EEF2F8")

    for column, label, color in zip(range(3), ("mag x", "mag y", "mag z"), COLORS[:3], strict=True):
        axes[1].plot(sample_index, magnetometer[:, column], color=color, linewidth=1.8, label=label)
    axes[1].set_ylabel("µT", color="#CBD3DF")
    axes[1].legend(loc="upper right", ncol=3, frameon=True, facecolor="#0B0F14", edgecolor="#3A414C", labelcolor="#EEF2F8")

    for column, color in enumerate(COLORS):
        axes[2].plot(sample_index, solar[:, column], color=color, linewidth=1.6, label=f"solar {column + 1}")
    axes[2].set_ylabel("scaled units", color="#CBD3DF")
    axes[2].set_xlabel("decoded sample row index", color="#CBD3DF")
    axes[2].legend(loc="upper right", ncol=3, frameon=True, facecolor="#0B0F14", edgecolor="#3A414C", labelcolor="#EEF2F8")

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    figure.tight_layout(rect=(0.02, 0.04, 0.995, 0.92), h_pad=1.1)
    figure.savefig(OUTPUT, dpi=300, facecolor=figure.get_facecolor(), bbox_inches="tight")
    plt.close(figure)
    print(f"wrote {OUTPUT} ({len(gyro)} decoded rows)")


if __name__ == "__main__":
    main()
