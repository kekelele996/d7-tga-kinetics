"""TGA decomposition analysis prototype.

Reads one JSON request and writes one JSON response. The laboratory tool was
originally written for a single first-order process and is being extended for
multi-heating-rate solid-state decomposition work.
"""

from __future__ import annotations

import json
import math
import sys

R = 8.31446261815324


def _clamp_unit(value: float) -> float:
    return min(max(float(value), 1.0e-14), 1.0 - 1.0e-14)


def _first_order_extent(alpha: float) -> float:
    return -math.log1p(-_clamp_unit(alpha))


def _first_order_inverse(extent: float) -> float:
    return 1.0 - math.exp(-max(float(extent), 0.0))


def _arrhenius(temperature_k: float, activation_energy_kj_mol: float) -> float:
    return math.exp(-(activation_energy_kj_mol * 1000.0) / (R * temperature_k))


def _temperature_integral(
    start_k: float, end_k: float, activation_energy_kj_mol: float, intervals: int = 8
) -> float:
    if end_k == start_k:
        return 0.0
    n = max(2, int(intervals))
    if n % 2:
        n += 1
    step = (end_k - start_k) / n
    total = _arrhenius(start_k, activation_energy_kj_mol)
    total += _arrhenius(end_k, activation_energy_kj_mol)
    for index in range(1, n):
        total += (4.0 if index % 2 else 2.0) * _arrhenius(
            start_k + index * step, activation_energy_kj_mol
        )
    return total * step / 3.0


def _cumulative_integral(
    temperatures_k: list[float], activation_energy_kj_mol: float
) -> list[float]:
    values = [0.0]
    for index in range(1, len(temperatures_k)):
        values.append(
            values[-1]
            + _temperature_integral(
                temperatures_k[index - 1],
                temperatures_k[index],
                activation_energy_kj_mol,
            )
        )
    return values


def _best_prefactor(runs: list[dict], activation_energy: float) -> tuple[float, float]:
    prepared = []
    for run in runs:
        temperatures = [float(value) + 273.15 for value in run["temperature_c"]]
        beta = float(run["beta_k_per_min"])
        extents = [value / beta for value in _cumulative_integral(temperatures, activation_energy)]
        prepared.append((extents, [float(value) for value in run["conversion"]]))

    def objective(log_prefactor: float) -> float:
        prefactor = 10.0 ** log_prefactor
        error = 0.0
        for extents, observed_values in prepared:
            for extent, observed in zip(extents, observed_values):
                residual = _first_order_inverse(prefactor * extent) - observed
                error += residual * residual
        return error

    coarse = min((objective(-10.0 + 0.25 * index), -10.0 + 0.25 * index) for index in range(145))
    low = max(-12.0, coarse[1] - 0.5)
    high = min(28.0, coarse[1] + 0.5)
    ratio = (math.sqrt(5.0) - 1.0) / 2.0
    left = high - ratio * (high - low)
    right = low + ratio * (high - low)
    f_left = objective(left)
    f_right = objective(right)
    for _ in range(30):
        if f_left < f_right:
            high = right
            right = left
            f_right = f_left
            left = high - ratio * (high - low)
            f_left = objective(left)
        else:
            low = left
            left = right
            f_left = f_right
            right = low + ratio * (high - low)
            f_right = objective(right)
    log_prefactor = (low + high) / 2.0
    return objective(log_prefactor), log_prefactor


def _fit_first_order(runs: list[dict]) -> tuple[float, float, float]:
    candidates = []
    for activation_energy in range(40, 401, 10):
        error, log_prefactor = _best_prefactor(runs, float(activation_energy))
        candidates.append((error, float(activation_energy), log_prefactor))
    seeds = sorted(candidates)[:3]
    low = max(30.0, min(item[1] for item in seeds) - 10.0)
    high = min(430.0, max(item[1] for item in seeds) + 10.0)
    ratio = (math.sqrt(5.0) - 1.0) / 2.0
    left = high - ratio * (high - low)
    right = low + ratio * (high - low)
    left_error, left_log = _best_prefactor(runs, left)
    right_error, right_log = _best_prefactor(runs, right)
    for _ in range(25):
        if left_error < right_error:
            high = right
            right = left
            right_error = left_error
            right_log = left_log
            left = high - ratio * (high - low)
            left_error, left_log = _best_prefactor(runs, left)
        else:
            low = left
            left = right
            left_error = right_error
            left_log = right_log
            right = low + ratio * (high - low)
            right_error, right_log = _best_prefactor(runs, right)
    candidates.append((left_error, left, left_log))
    candidates.append((right_error, right, right_log))
    return min(candidates)


def normalize_runs(runs: list[dict]) -> dict:
    output = []
    for run in runs:
        initial = float(run["initial_mass_mg"])
        conversion = [
            (initial - float(point["mass_mg"])) / initial for point in run["points"]
        ]
        output.append({"run_id": run.get("run_id", ""), "conversion": conversion})
    return {"runs": output}


def identify_runs(runs: list[dict]) -> dict:
    error, activation_energy, log_prefactor = _fit_first_order(runs)
    return {
        "model": "nth_order_1",
        "activation_energy_kj_mol": activation_energy,
        "log10_prefactor_s": log_prefactor,
        "model_scores": {"nth_order_1": error},
    }


def predict_conversion(request: dict) -> dict:
    activation_energy = float(request["activation_energy_kj_mol"])
    prefactor = 10.0 ** float(request["log10_prefactor_s"])
    schedule = request["schedule"]
    kind = schedule["type"]
    if kind == "linear":
        start_k = float(schedule["start_temperature_c"]) + 273.15
        end_k = float(schedule["end_temperature_c"]) + 273.15
        beta_k_s = float(schedule["beta_k_per_min"]) / 60.0
        extent = prefactor * _temperature_integral(start_k, end_k, activation_energy) / beta_k_s
    elif kind == "isothermal":
        temperature_k = float(schedule["temperature_c"]) + 273.15
        extent = (
            prefactor
            * _arrhenius(temperature_k, activation_energy)
            * float(schedule["time_min"])
            * 60.0
        )
    else:
        raise ValueError(f"unknown schedule type: {kind}")
    return {"conversion": _first_order_inverse(extent)}


def dispatch(request: dict) -> dict:
    operation = request["operation"]
    if operation == "normalize":
        return normalize_runs(request["runs"])
    if operation == "identify":
        return identify_runs(request["runs"])
    if operation == "predict":
        return predict_conversion(request)
    raise ValueError(f"unknown operation: {operation}")


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print("usage: python3 /app/main.py <request.json>", file=sys.stderr)
        return 2
    try:
        with open(argv[1], "r", encoding="utf-8") as handle:
            payload = json.load(handle)
        result = dispatch(payload)
    except (OSError, json.JSONDecodeError, KeyError, TypeError, ValueError, OverflowError) as exc:
        print(f"input error: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(result, ensure_ascii=False, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
