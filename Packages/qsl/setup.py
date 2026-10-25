import re
import subprocess
import sys
from pathlib import Path
from typing import Optional

from setuptools import setup
from wheel.bdist_wheel import bdist_wheel


def bundled_platform_tag() -> Optional[str]:
    executable = Path(__file__).parent / "src" / "qsl" / "bin" / "qfactor.x"
    if sys.platform != "darwin" or not executable.is_file():
        return None
    architecture = subprocess.check_output(
        ["lipo", "-archs", executable], text=True
    ).strip()
    load_commands = subprocess.check_output(["otool", "-l", executable], text=True)
    match = re.search(r"\bminos\s+(\d+)\.(\d+)", load_commands)
    if not match:
        return None
    major, minor = match.groups()
    return f"macosx_{major}_{minor}_{architecture}"


class PlatformWheel(bdist_wheel):
    """Mark the wheel as platform-specific because it bundles qfactor.x."""

    def finalize_options(self) -> None:
        super().finalize_options()
        self.root_is_pure = False

    def get_tag(self) -> tuple[str, str, str]:
        _, _, platform = super().get_tag()
        return "py3", "none", bundled_platform_tag() or platform


setup(cmdclass={"bdist_wheel": PlatformWheel})
