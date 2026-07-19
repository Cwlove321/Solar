#!/usr/bin/env python3
"""读取 ``helicity`` 生成的 out.binary，并可导出 CSV。"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


BASE_COLUMNS = ("energy", "free_energy", "h_v", "h_j", "b_flux")
ERROR_COLUMNS = ("divergence_error", "curl_a_error")


def read_results(path: str | Path, error_estimation: bool = False) -> tuple[tuple[str, ...], np.ndarray]:
    columns = BASE_COLUMNS + (ERROR_COLUMNS if error_estimation else ())
    values = np.fromfile(Path(path), dtype=np.float64)
    if values.size % len(columns) != 0:
        raise ValueError(
            f"文件包含 {values.size} 个 float64，不能按每帧 {len(columns)} 个量解析；"
            "请检查 --error-estimation 是否与计算时一致。"
        )
    return columns, values.reshape(-1, len(columns))


def main() -> None:
    parser = argparse.ArgumentParser(description="读取 Helicity out.binary 并显示或导出 CSV。")
    parser.add_argument("input", help="out.binary 路径")
    parser.add_argument("--start", type=int, default=0, help="第一行对应的帧编号，默认 0")
    parser.add_argument("--error-estimation", action="store_true", help="每帧按 7 个量解析")
    parser.add_argument("--csv", help="可选 CSV 输出路径")
    args = parser.parse_args()

    columns, values = read_results(args.input, args.error_estimation)
    frames = np.arange(args.start, args.start + len(values), dtype=np.int64)[:, None]
    table = np.column_stack((frames, values))
    header = "frame," + ",".join(columns)

    if args.csv:
        destination = Path(args.csv)
        destination.parent.mkdir(parents=True, exist_ok=True)
        np.savetxt(destination, table, delimiter=",", header=header, comments="", fmt=["%d"] + ["%.16e"] * len(columns))
        print(destination)
    else:
        print(header)
        for row in table:
            print(f"{int(row[0])}," + ",".join(f"{value:.16e}" for value in row[1:]))


if __name__ == "__main__":
    main()
