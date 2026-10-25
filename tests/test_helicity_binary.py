from pathlib import Path

import numpy as np
import pytest

from workflows.Helicity.helicity_binary import (
    FLOAT_DTYPE,
    INT_DTYPE,
    boundary_faces,
    read_boundary_data,
    read_magnetic_data,
    write_boundary_data,
    write_magnetic_data,
)


def sample_fields() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    base = np.arange(2 * 3 * 4, dtype=np.float64).reshape(2, 3, 4)
    return base, base + 100.0, base + 200.0


def test_magnetic_stream_matches_fortran_layout(tmp_path: Path) -> None:
    b1, b2, b3 = sample_fields()
    output = tmp_path / "nlfff_0000.binary"

    write_magnetic_data(output, 0.25, b1, b2, b3)

    expected = b"".join(
        (
            np.asarray(0.25, dtype=FLOAT_DTYPE).tobytes(),
            np.asarray((4, 3, 2), dtype=INT_DTYPE).tobytes(),
            np.asarray(b1, dtype=FLOAT_DTYPE, order="C").tobytes(order="C"),
            np.asarray(b2, dtype=FLOAT_DTYPE, order="C").tobytes(order="C"),
            np.asarray(b3, dtype=FLOAT_DTYPE, order="C").tobytes(order="C"),
        )
    )
    assert output.read_bytes() == expected
    spacing, actual_b1, actual_b2, actual_b3 = read_magnetic_data(output)
    assert spacing == 0.25
    np.testing.assert_array_equal(actual_b1, b1)
    np.testing.assert_array_equal(actual_b2, b2)
    np.testing.assert_array_equal(actual_b3, b3)


def test_boundary_stream_has_expected_signs_order_and_shapes(tmp_path: Path) -> None:
    b1, b2, b3 = sample_fields()
    output = tmp_path / "bd_0000.binary"

    flux_sum = write_boundary_data(output, b1, b2, b3)
    expected_faces = boundary_faces(b1, b2, b3)
    actual_faces = read_boundary_data(output, b1.shape)

    assert flux_sum == pytest.approx(sum(np.sum(face) for face in expected_faces))
    for actual, expected in zip(actual_faces, expected_faces, strict=True):
        np.testing.assert_array_equal(actual, expected)
    assert output.stat().st_size == 8 * 2 * (3 * 2 + 4 * 2 + 4 * 3)


@pytest.mark.parametrize("kind", ["truncated", "trailing"])
def test_magnetic_reader_rejects_invalid_file_length(tmp_path: Path, kind: str) -> None:
    b1, b2, b3 = sample_fields()
    output = tmp_path / "nlfff_0000.binary"
    write_magnetic_data(output, 0.25, b1, b2, b3)
    content = output.read_bytes()
    output.write_bytes(content[:-1] if kind == "truncated" else content + b"x")

    with pytest.raises(ValueError):
        read_magnetic_data(output)


def test_writer_rejects_mismatched_or_non_3d_fields(tmp_path: Path) -> None:
    b1, b2, b3 = sample_fields()
    with pytest.raises(ValueError, match="形状必须一致"):
        write_magnetic_data(tmp_path / "bad.binary", 1.0, b1, b2[:, :, :-1], b3)
    with pytest.raises(ValueError, match="必须是三维"):
        write_boundary_data(tmp_path / "bad-bd.binary", b1[0], b2[0], b3[0])


def test_failed_write_does_not_replace_existing_file(tmp_path: Path) -> None:
    b1, b2, b3 = sample_fields()
    output = tmp_path / "nlfff_0000.binary"
    output.write_bytes(b"existing")

    with pytest.raises(ValueError):
        write_magnetic_data(output, 1.0, b1, b2[:, :, :-1], b3)

    assert output.read_bytes() == b"existing"
