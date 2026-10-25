"""Movie helpers for turning saved figure frames into animations."""

from __future__ import annotations

from pathlib import Path
import shutil
import subprocess

import numpy as np


MOVIE_QUALITY_PRESETS = {
    "simple": {
        "description": "Small preview file, fast encoding.",
        "crf": "28",
        "preset": "veryfast",
    },
    "high": {
        "description": "Higher visual quality with larger output files.",
        "crf": "12",
        "preset": "slow",
    },
}


def make_movie_from_images(
    image_paths,
    output: str | Path,
    *,
    fps: float = 6.0,
    backend: str = "auto",
    quality: str = "simple",
    overwrite: bool = False,
):
    """Create a movie from an ordered list of image paths.

    Parameters
    ----------
    image_paths
        Ordered image paths. The order is preserved in the output movie.
    output
        Output movie path, usually `.mp4` or `.gif`.
    fps
        Frames per second.
    backend
        `auto`, `ffmpeg`, or `imageio`. `auto` uses `ffmpeg` when available and
        falls back to `imageio`.
    quality
        Encoding quality preset. `simple` creates a compact preview movie;
        `high` uses lower compression for clearer frames.
    overwrite
        Replace an existing output file.
    """

    if quality not in MOVIE_QUALITY_PRESETS:
        choices = ", ".join(sorted(MOVIE_QUALITY_PRESETS))
        raise ValueError(f"quality must be one of: {choices}")

    image_paths = [Path(path) for path in image_paths]
    if not image_paths:
        raise ValueError("image_paths must contain at least one frame")
    missing = [path for path in image_paths if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing movie frame: {missing[0]}")

    output = Path(output)
    if output.exists() and not overwrite:
        raise FileExistsError(f"Output already exists: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)

    if backend == "auto":
        backend = "ffmpeg" if shutil.which("ffmpeg") else "imageio"
    if backend == "ffmpeg":
        return _make_movie_with_ffmpeg(image_paths, output, fps=fps, quality=quality, overwrite=overwrite)
    if backend == "imageio":
        return _make_movie_with_imageio(image_paths, output, fps=fps, quality=quality, overwrite=overwrite)
    raise ValueError("backend must be 'auto', 'ffmpeg', or 'imageio'")


def _make_movie_with_ffmpeg(
    image_paths: list[Path],
    output: Path,
    *,
    fps: float,
    quality: str,
    overwrite: bool,
) -> Path:
    ffmpeg = shutil.which("ffmpeg")
    if ffmpeg is None:
        raise RuntimeError("ffmpeg executable was not found")

    list_path = output.with_suffix(output.suffix + ".frames.txt")
    duration = 1.0 / fps
    with list_path.open("w", encoding="utf-8") as handle:
        for path in image_paths:
            handle.write(f"file '{path.resolve()}'\n")
            handle.write(f"duration {duration:.8f}\n")
        handle.write(f"file '{image_paths[-1].resolve()}'\n")

    command = [
        ffmpeg,
        "-y" if overwrite else "-n",
        "-f",
        "concat",
        "-safe",
        "0",
        "-i",
        str(list_path),
        "-c:v",
        "libx264",
        "-preset",
        MOVIE_QUALITY_PRESETS[quality]["preset"],
        "-crf",
        MOVIE_QUALITY_PRESETS[quality]["crf"],
        "-vf",
        "format=yuv420p",
        "-movflags",
        "+faststart",
        str(output),
    ]
    try:
        subprocess.run(command, check=True, capture_output=True, text=True)
    finally:
        list_path.unlink(missing_ok=True)
    return output


def _make_movie_with_imageio(
    image_paths: list[Path],
    output: Path,
    *,
    fps: float,
    quality: str,
    overwrite: bool,
) -> Path:
    if output.exists() and overwrite:
        output.unlink()
    try:
        import imageio.v2 as iio
    except ImportError as exc:
        raise RuntimeError("imageio is required for the imageio movie backend") from exc

    suffix = output.suffix.lower()
    if suffix == ".gif":
        with iio.get_writer(output, mode="I", duration=1000.0 / fps, loop=0) as writer:
            for path in image_paths:
                writer.append_data(iio.imread(path))
        return output
    if suffix != ".mp4":
        raise ValueError("the imageio backend supports only .mp4 and .gif outputs")

    try:
        import imageio_ffmpeg  # noqa: F401
    except ImportError as exc:
        raise RuntimeError(
            "MP4 output with the imageio backend requires imageio-ffmpeg; "
            "install solarutils[movie]"
        ) from exc

    ffmpeg_params = [
        "-crf",
        MOVIE_QUALITY_PRESETS[quality]["crf"],
        "-preset",
        MOVIE_QUALITY_PRESETS[quality]["preset"],
        "-pix_fmt",
        "yuv420p",
        "-movflags",
        "+faststart",
    ]
    with iio.get_writer(
        output,
        format="FFMPEG",
        fps=fps,
        codec="libx264",
        macro_block_size=1,
        ffmpeg_params=ffmpeg_params,
    ) as writer:
        for path in image_paths:
            writer.append_data(_pad_frame_to_even_shape(iio.imread(path)))
    return output


def _pad_frame_to_even_shape(frame):
    """Pad image frames to even pixel dimensions for common MP4 encoders."""

    height, width = frame.shape[:2]
    pad_height = height % 2
    pad_width = width % 2
    if not pad_height and not pad_width:
        return frame
    pad_spec = [(0, pad_height), (0, pad_width)]
    pad_spec.extend((0, 0) for _ in frame.shape[2:])
    return np.pad(frame, pad_spec, mode="edge")
