"""Prevent the standalone NLFFF and TMF bundles from drifting apart."""

from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
NLFFF_CORE = REPO_ROOT / "workflows/Cartesian_NLFFF/Preprocess/nlfff_core"
TMF_CORE = REPO_ROOT / "workflows/Cartesian_TMF/Preprocess/nlfff_core"
SPHERICAL_CORE = REPO_ROOT / "workflows/Spherical_RBSL/Preprocess/rbsl_core"


def source_files(root: Path) -> dict[Path, bytes]:
    return {
        path.relative_to(root): path.read_bytes()
        for path in root.rglob("*")
        if path.is_file() and path.suffix == ".py"
    }


def test_standalone_nlfff_cores_are_identical() -> None:
    """Both handoff bundles must ship the same tested NLFFF/downsampling core."""

    nlfff_files = source_files(NLFFF_CORE)
    tmf_files = source_files(TMF_CORE)

    assert nlfff_files, f"No Python sources found in {NLFFF_CORE}"
    assert nlfff_files.keys() == tmf_files.keys()

    different = [path for path in nlfff_files if nlfff_files[path] != tmf_files[path]]
    assert not different, (
        "Cartesian_NLFFF and Cartesian_TMF contain different nlfff_core sources: "
        + ", ".join(str(path) for path in different)
    )


def test_standalone_spherical_core_matches_solarutils() -> None:
    """The portable spherical bundle must stay synchronized with solarutils."""

    pairs = {
        REPO_ROOT / "src/solarutils/sim/spherical.py": SPHERICAL_CORE / "spherical.py",
        REPO_ROOT / "src/solarutils/obs/maps.py": SPHERICAL_CORE / "maps.py",
    }
    different = [
        f"{source.relative_to(REPO_ROOT)} != {vendored.relative_to(REPO_ROOT)}"
        for source, vendored in pairs.items()
        if source.read_bytes() != vendored.read_bytes()
    ]
    assert not different, "Standalone Spherical_RBSL sources have drifted: " + "; ".join(different)
