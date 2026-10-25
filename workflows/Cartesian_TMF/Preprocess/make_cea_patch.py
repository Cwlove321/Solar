#!/usr/bin/env python3
"""Generate CEA Br/Bt/Bp with the shared, validated official JSOC core.

Keep Packages/sharp alongside workflows when copying this workflow.
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from Packages.sharp.workflow import infer_p_angle_degree, main, parse_args

if __name__ == "__main__":
    raise SystemExit(main())
