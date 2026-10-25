#!/usr/bin/env python3
"""Write a compact JSON/CSV health report for a running AMRVAC case."""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from simulation_analysis.config import load_config
from simulation_analysis.logs import build_report, write_report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--job-id")
    parser.add_argument("--watch", action="store_true")
    parser.add_argument("--interval", type=int, default=600)
    args = parser.parse_args()
    config = load_config(args.config)
    while True:
        report = build_report(config, args.job_id)
        write_report(config, report)
        print(json.dumps(report, indent=2, allow_nan=True), flush=True)
        if not args.watch or report["job_status"] not in {"PEND", "RUN"}:
            break
        time.sleep(args.interval)


if __name__ == "__main__":
    main()
