#!/usr/bin/env python3
"""Build an H.264 MP4 from a consistently named PNG frame sequence."""

from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--frames", required=True, help="ffmpeg glob, for example frames/*.png")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--fps", type=float, default=6.0)
    args = parser.parse_args()
    if shutil.which("ffmpeg") is None:
        raise SystemExit("ffmpeg is not on PATH")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    command = [
        "ffmpeg", "-y", "-framerate", str(args.fps), "-pattern_type", "glob",
        "-i", args.frames, "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2",
        "-c:v", "libx264", "-pix_fmt", "yuv420p", str(args.output)
    ]
    subprocess.run(command, check=True)
    print(args.output.resolve(), flush=True)


if __name__ == "__main__":
    main()
