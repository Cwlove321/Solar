from pathlib import Path

import numpy as np
import pytest
import imageio.v2 as imageio

from solarutils.plotting.movie import _pad_frame_to_even_shape, make_movie_from_images


def test_make_movie_from_images_rejects_empty_frame_list(tmp_path):
    with pytest.raises(ValueError, match="at least one"):
        make_movie_from_images([], tmp_path / "movie.mp4")


def test_make_movie_from_images_rejects_missing_frame(tmp_path):
    with pytest.raises(FileNotFoundError, match="Missing movie frame"):
        make_movie_from_images([tmp_path / "missing.png"], tmp_path / "movie.mp4")


def test_make_movie_from_images_rejects_existing_output(tmp_path):
    frame = tmp_path / "frame.png"
    frame.write_bytes(b"not a real png")
    output = tmp_path / "movie.mp4"
    output.write_bytes(b"existing")

    with pytest.raises(FileExistsError, match="already exists"):
        make_movie_from_images([frame], output, backend="imageio")


def test_make_movie_from_images_rejects_unknown_quality(tmp_path):
    with pytest.raises(ValueError, match="quality must be one of"):
        make_movie_from_images([], tmp_path / "movie.mp4", quality="giant")


def test_make_movie_from_images_supports_high_quality_imageio(tmp_path):
    frame_paths = []
    for index in range(2):
        frame = np.full((5, 7, 3), index * 80, dtype=np.uint8)
        path = tmp_path / f"frame_{index:02d}.png"
        imageio.imwrite(path, frame)
        frame_paths.append(path)

    output = tmp_path / "movie.mp4"
    make_movie_from_images(frame_paths, output, fps=2, backend="imageio", quality="high")

    assert output.exists()
    assert output.stat().st_size > 0


def test_pad_frame_to_even_shape_pads_odd_dimensions():
    frame = np.ones((3, 5, 4), dtype=np.uint8)

    padded = _pad_frame_to_even_shape(frame)

    assert padded.shape == (4, 6, 4)
    assert np.array_equal(padded[:3, :5], frame)
