#!/usr/bin/env bash
# Reference implementation installer.
#
# HTB_SKIP=<difficulty id>[,<difficulty id>...] reverts selected fixes while
# installing all others. Each replacement must match exactly once so a changed
# reference implementation cannot silently invalidate the leave-one-out runs.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
install -m 0644 "$HERE/app/main.py" /app/main.py

SKIP="${HTB_SKIP:-}" /usr/local/bin/python3 - /app/main.py <<'PY'
import os
import sys
from pathlib import Path

REVERT = {
    "difficulty_1": (
        '''        conversion = [
            (initial - float(point["mass_mg"])) / removable for point in run["points"]
        ]''',
        '''        conversion = [
            (initial - float(point["mass_mg"])) / initial for point in run["points"]
        ]''',
    ),
    "difficulty_2": (
        '''        beta_k_s = float(run["beta_k_per_min"]) / 60.0
        extents = [
            value / beta_k_s
            for value in _cumulative_integral(temperatures, activation_energy_kj_mol)
        ]''',
        '''        beta_k_s = float(run["beta_k_per_min"])
        extents = [
            value / beta_k_s
            for value in _cumulative_integral(temperatures, activation_energy_kj_mol)
        ]''',
    ),
    "difficulty_3": (
        '''MODEL_NAMES = (
    "nth_order_1",
    "nth_order_2",
    "avrami_2",
    "contracting_sphere",
)''',
        '''MODEL_NAMES = ("nth_order_1",)''',
    ),
}
skip = [item for item in os.environ.get("SKIP", "").split(",") if item]
unknown = sorted(set(skip) - set(REVERT))
if unknown:
    sys.exit(f"HTB_SKIP contains unknown difficulty ids: {unknown}")
path = Path(sys.argv[1])
text = path.read_text(encoding="utf-8")
for key in skip:
    fixed, original = REVERT[key]
    if text.count(fixed) != 1:
        sys.exit(f"{key} replacement fragment occurs {text.count(fixed)} times, expected 1")
    text = text.replace(fixed, original)
    print(f"leave-one-out: reverted {key}", file=sys.stderr)
path.write_text(text, encoding="utf-8")
PY
chown runner:runner /app/main.py 2>/dev/null || true
