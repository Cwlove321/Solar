#!/usr/bin/env python3
"""将 AMRVAC ``.dat`` 快照转换为磁螺度程序所需的数据文件。"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np


# 优先使用从原作者仓库拉取到 solar_code/Packages 的新版 simesh。
# 安装后运行时，下面的路径注入也不会改变 API。
REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
LOCAL_SIMESH_SRC = REPOSITORY_ROOT / "Packages" / "simesh" / "src"
if LOCAL_SIMESH_SRC.is_dir():
    sys.path.insert(0, str(LOCAL_SIMESH_SRC))

try:
    from simesh.amrvac import open_dataset, read_uniform
except (ImportError, ModuleNotFoundError) as exc:  # pragma: no cover - 环境错误提示
    raise ModuleNotFoundError(
        "无法导入新版 simesh。请在 solar_code 根目录执行：\n"
        "  bash bootstrap_dependencies.sh\n"
        "  python3 -m pip install -e Packages/simesh\n"
        "原作者仓库：https://github.com/Astery0502/simesh"
    ) from exc

try:
    from .helicity_binary import (
        read_magnetic_data,
        write_boundary_data,
        write_helicity_inputs,
        write_magnetic_data,
    )
except ImportError:  # 直接执行 process_dat.py 时使用同目录模块
    from helicity_binary import (
        read_magnetic_data,
        write_boundary_data,
        write_helicity_inputs,
        write_magnetic_data,
    )


def init_boundary(filename: str | Path, b1: np.ndarray, b2: np.ndarray, b3: np.ndarray) -> float:
    """按 Fortran 程序要求写出六个外法向磁场边界面。

    输入磁场数组的轴顺序必须为 ``(z, y, x)``。返回六个面的离散通量和，
    可用作边界通量平衡的快速诊断值。
    """
    return write_boundary_data(filename, b1, b2, b3)


def write_binary_data(
    filename: str | Path,
    dx: float,
    b1: np.ndarray,
    b2: np.ndarray,
    b3: np.ndarray,
) -> None:
    """写出 ``dx, nx, ny, nz, Bx, By, Bz`` 的 Fortran stream 文件。"""
    write_magnetic_data(filename, dx, b1, b2, b3)


def read_binary_data(filename: str | Path) -> tuple[float, np.ndarray, np.ndarray, np.ndarray]:
    """读取 :func:`write_binary_data` 产生的文件，主要用于检查与测试。"""
    return read_magnetic_data(filename)


def curl(
    b1: np.ndarray,
    b2: np.ndarray,
    b3: np.ndarray,
    dx: float = 1.0,
    dy: float = 1.0,
    dz: float = 1.0,
    dim: str = "zyx",
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """计算磁场旋度；本工作流默认数组顺序为 ``(z, y, x)``。"""
    if dim == "zyx":
        curl_x = np.gradient(b3, axis=1) / dy - np.gradient(b2, axis=0) / dz
        curl_y = np.gradient(b1, axis=0) / dz - np.gradient(b3, axis=2) / dx
        curl_z = np.gradient(b2, axis=2) / dx - np.gradient(b1, axis=1) / dy
    else:
        curl_x = np.gradient(b3, axis=1) / dy - np.gradient(b2, axis=2) / dz
        curl_y = np.gradient(b1, axis=2) / dz - np.gradient(b3, axis=0) / dx
        curl_z = np.gradient(b2, axis=0) / dx - np.gradient(b1, axis=1) / dy
    return curl_x, curl_y, curl_z


def load_uniform_snapshot(
    filename: str | Path,
    resolution_multiplier: int,
    interpolation: str = "zero",
    field_names: set[str] | None = None,
    allow_anisotropic: bool = False,
) -> tuple[SimpleNamespace, tuple[int, int, int], np.ndarray, np.ndarray, np.ndarray]:
    """使用新版 simesh 读取快照，并转换为下游使用的 ``(z,y,x)`` 布局。"""
    if resolution_multiplier < 1:
        raise ValueError("level 必须是大于等于 1 的整数")

    dataset = open_dataset(str(filename))
    if int(dataset.ndim) != 3:
        raise ValueError(f"Helicity 工作流只支持 Cartesian 3D 数据，实际 ndim={dataset.ndim}")

    resolution_array = resolution_multiplier * np.asarray(dataset.domain_nx, dtype=np.int64)
    resolution = tuple(int(value) for value in resolution_array)
    xmin = np.asarray(dataset.physical_domain[0], dtype=np.float64)
    xmax = np.asarray(dataset.physical_domain[1], dtype=np.float64)

    available_names = list(dataset.wnames)
    if field_names is None:
        selected_indices = list(range(len(available_names)))
    else:
        missing_requested = field_names.difference(available_names)
        if missing_requested:
            raise KeyError(f"快照缺少所需变量：{sorted(missing_requested)}")
        selected_indices = [index for index, name in enumerate(available_names) if name in field_names]
    selected_names = [available_names[index] for index in selected_indices]

    ghost_width = 1 if interpolation == "linear" else 0
    uniform = read_uniform(
        str(filename),
        resolution=resolution,
        bounds=(xmin, xmax),
        field_indices=selected_indices,
        ghost_width=ghost_width,
        interpolation=interpolation,
    )
    if uniform.shape != (*resolution, len(selected_names)):
        raise ValueError(
            "simesh 返回了非预期形状："
            f"{uniform.shape}，预期 {(*resolution, len(selected_names))}"
        )

    fields = {
        name: np.ascontiguousarray(np.transpose(uniform[..., index], (2, 1, 0)))
        for index, name in enumerate(selected_names)
    }
    missing = {"b1", "b2", "b3"}.difference(fields)
    if missing:
        raise KeyError(f"快照缺少 Helicity 必需的磁场变量：{sorted(missing)}")

    spacing = (xmax - xmin) / resolution_array
    if not np.allclose(spacing, spacing[0]):
        message = (
            "helicity.f90 只保存一个 dx；"
            f"实际 (dx,dy,dz)={spacing.tolist()}"
        )
        if not allow_anisotropic:
            raise ValueError(message + "。若要沿用旧流程统一使用 dx，请添加 --allow-anisotropic")
        print(f"警告：{message}；Fortran 输入将统一使用 dx={spacing[0]:.12g}")
    return SimpleNamespace(**fields), resolution, xmin, xmax, spacing


def process_data(
    path: str | Path,
    out_path: str | Path,
    start: int,
    end: int | None,
    level: int = 2,
    save_nlfff: bool = False,
    save_bd: bool = False,
    save_vtk: bool = False,
    save_npz: bool = False,
    interpolation: str = "zero",
    velocity_unit: float = 116.45,
    pattern: str = "*.dat",
    allow_anisotropic: bool = False,
) -> None:
    """转换编号区间 ``[start, end)`` 内的 AMRVAC 快照。"""
    input_path = Path(path).expanduser().resolve()
    output_path = Path(out_path).expanduser().resolve()
    files = sorted(input_path.glob(pattern))
    if not files:
        raise FileNotFoundError(f"输入目录中没有匹配 {pattern!r} 的文件：{input_path}")
    if end is None:
        end = len(files)
    if start < 0 or end < start or end > len(files):
        raise ValueError(f"索引范围应满足 0 <= start <= end <= {len(files)}，实际为 {start}:{end}")
    if not any((save_nlfff, save_bd, save_vtk, save_npz)):
        raise ValueError("至少选择一个输出选项：--save-nlfff/--save-bd/--save-vtk/--save-npz")

    selected_directories = {
        "nlfff": save_nlfff,
        "bd": save_bd,
        "vtu": save_vtk,
        "npz": save_npz,
    }
    for name, selected in selected_directories.items():
        if selected:
            (output_path / name).mkdir(parents=True, exist_ok=True)

    if save_vtk:
        try:
            from pyevtk.hl import gridToVTK
        except ModuleNotFoundError as exc:
            raise ModuleNotFoundError("保存 VTK 需要安装 pyevtk：python -m pip install pyevtk") from exc
    else:
        gridToVTK = None

    for index in range(start, end):
        snapshot = files[index]
        requested_fields = None if (save_vtk or save_npz) else {"b1", "b2", "b3"}
        w, (nx, ny, nz), xmin, xmax, spacing = load_uniform_snapshot(
            snapshot,
            resolution_multiplier=level,
            interpolation=interpolation,
            field_names=requested_fields,
            allow_anisotropic=allow_anisotropic,
        )
        dx, dy, dz = (float(value) for value in spacing)
        number = f"{index:04d}"

        if save_nlfff and save_bd:
            flux_sum = write_helicity_inputs(
                output_path / "nlfff" / f"nlfff_{number}.binary",
                output_path / "bd" / f"bd_{number}.binary",
                dx,
                w.b1,
                w.b2,
                w.b3,
            )
            print(f"[{number}] 边界离散通量和：{flux_sum:.8e}")
        elif save_nlfff:
            write_binary_data(output_path / "nlfff" / f"nlfff_{number}.binary", dx, w.b1, w.b2, w.b3)
        elif save_bd:
            flux_sum = init_boundary(output_path / "bd" / f"bd_{number}.binary", w.b1, w.b2, w.b3)
            print(f"[{number}] 边界离散通量和：{flux_sum:.8e}")

        if save_vtk:
            required_fluid_fields = {"m1", "m2", "m3", "rho", "e"}
            missing = sorted(name for name in required_fluid_fields if not hasattr(w, name))
            if missing:
                raise KeyError(f"保存 VTK 所需变量缺失：{missing}")

            j1, j2, j3 = curl(w.b1, w.b2, w.b3, dx=dx, dy=dy, dz=dz, dim="zyx")
            current = np.sqrt(j1**2 + j2**2 + j3**2)
            v1 = w.m1 / w.rho * velocity_unit
            v2 = w.m2 / w.rho * velocity_unit
            v3 = w.m3 / w.rho * velocity_unit

        if save_vtk and gridToVTK is not None:
            # pyevtk 需要 (x,y,z) 数组，因此将内部 (z,y,x) 布局转置回去。
            x = np.linspace(xmin[0], xmax[0], nx)
            y = np.linspace(xmin[1], xmax[1], ny)
            z = np.linspace(xmin[2], xmax[2], nz)
            gridToVTK(
                str(output_path / "vtu" / f"eru{number}"),
                x,
                y,
                z,
                pointData={
                    "b1": np.ascontiguousarray(w.b1.T),
                    "b2": np.ascontiguousarray(w.b2.T),
                    "b3": np.ascontiguousarray(w.b3.T),
                    "j": np.ascontiguousarray(current.T),
                    "e": np.ascontiguousarray(w.e.T),
                    "rho": np.ascontiguousarray(w.rho.T),
                    "v1": np.ascontiguousarray(v1.T),
                    "v2": np.ascontiguousarray(v2.T),
                    "v3": np.ascontiguousarray(v3.T),
                },
            )

        if save_npz:
            excluded = {"b1", "b2", "b3"}
            filtered = {name: value for name, value in vars(w).items() if name not in excluded}
            np.savez_compressed(output_path / "npz" / f"eru{number}.npz", **filtered)

        print(f"[{number}] 已处理 {snapshot.name}，分辨率={nx}x{ny}x{nz}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="用新版 simesh 将 AMRVAC .dat 快照转换为 Helicity 输入文件。"
    )
    parser.add_argument("path", help="包含 .dat 快照的输入目录")
    parser.add_argument("out_path", help="输出目录")
    parser.add_argument("start", type=int, nargs="?", default=0, help="起始索引（包含），默认 0")
    parser.add_argument("end", type=int, nargs="?", default=None, help="结束索引（不包含），默认处理到最后")
    parser.add_argument("--level", type=int, default=2, help="相对 domain_nx 的分辨率倍数，默认 2")
    parser.add_argument(
        "--interpolation",
        choices=("zero", "linear"),
        default="zero",
        help="simesh 重采样方法；默认 zero（分片常数，稳健且省内存）",
    )
    parser.add_argument("--save-nlfff", "--save_nlfff", action="store_true", help="写出 nlfff/*.binary")
    parser.add_argument("--save-bd", "--save_bd", action="store_true", help="写出 bd/*.binary")
    parser.add_argument("--save-vtk", "--save_vtk", action="store_true", help="写出 vtu/*.vtr")
    parser.add_argument("--save-npz", "--save_npz", action="store_true", help="写出 npz/*.npz")
    parser.add_argument(
        "--velocity-unit",
        type=float,
        default=116.45,
        help="VTK 速度换算系数，默认沿用原流程的 116.45",
    )
    parser.add_argument(
        "--pattern",
        default="*.dat",
        help="输入文件匹配模式，默认 *.dat；例如 driven[0-9][0-9][0-9][0-9].dat",
    )
    parser.add_argument(
        "--allow-anisotropic",
        action="store_true",
        help="允许 dx、dy、dz 不完全相等，并沿用旧 Fortran 流程统一写入 dx",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    selected_any_output = any((args.save_nlfff, args.save_bd, args.save_vtk, args.save_npz))
    # 最常用的无开关调用直接生成 helicity.f90 所需的两类输入。
    save_nlfff = args.save_nlfff or not selected_any_output
    save_bd = args.save_bd or not selected_any_output
    process_data(
        args.path,
        args.out_path,
        args.start,
        args.end,
        level=args.level,
        save_nlfff=save_nlfff,
        save_bd=save_bd,
        save_vtk=args.save_vtk,
        save_npz=args.save_npz,
        interpolation=args.interpolation,
        velocity_unit=args.velocity_unit,
        pattern=args.pattern,
        allow_anisotropic=args.allow_anisotropic,
    )


if __name__ == "__main__":
    main()
