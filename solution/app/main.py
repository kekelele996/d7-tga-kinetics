"""TGA decomposition analysis utilities.

The command reads one JSON request and writes one JSON response. Temperatures are
converted to kelvin and heating rates from K/min to K/s at the boundary.
"""

from __future__ import annotations

import json
import math
import sys

R = 8.31446261815324
MODEL_NAMES = (
    "nth_order_1",
    "nth_order_2",
    "avrami_2",
    "contracting_sphere",
)


def _clamp_unit(value: float) -> float:
    return min(max(float(value), 1.0e-14), 1.0 - 1.0e-14)


def inverse_extent(model: str, extent: float) -> float:
    z = max(float(extent), 0.0)
    if model == "nth_order_1":
        return 1.0 - math.exp(-z)
    if model == "nth_order_2":
        return z / (1.0 + z)
    if model == "avrami_2":
        return 1.0 - math.exp(-z * z)
    if model == "contracting_sphere":
        return 1.0 - (1.0 - min(z, 1.0)) ** 3
    raise ValueError(f"unknown reaction model: {model}")


def _arrhenius(temperature_k: float, activation_energy_kj_mol: float) -> float:
    return math.exp(-(activation_energy_kj_mol * 1000.0) / (R * temperature_k))


def _temperature_integral(
    start_k: float, end_k: float, activation_energy_kj_mol: float, intervals: int = 32
) -> float:
    if end_k < start_k:
        return -_temperature_integral(end_k, start_k, activation_energy_kj_mol, intervals)
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


def _best_log_prefactor(
    model: str, runs: list[dict], activation_energy_kj_mol: float
) -> tuple[float, float]:
    prepared = []
    for run in runs:
        temperatures = [float(value) + 273.15 for value in run["temperature_c"]]
        beta_k_s = float(run["beta_k_per_min"]) / 60.0
        extents = [
            value / beta_k_s
            for value in _cumulative_integral(temperatures, activation_energy_kj_mol)
        ]
        prepared.append((extents, [float(value) for value in run["conversion"]]))

    def objective(log_prefactor: float) -> float:
        prefactor = 10.0 ** log_prefactor
        squared_error = 0.0
        for extents, observed_values in prepared:
            for extent, observed in zip(extents, observed_values):
                residual = inverse_extent(model, prefactor * extent) - observed
                squared_error += residual * residual
        return squared_error

    coarse = min((objective(-8.0 + 0.5 * index), -8.0 + 0.5 * index) for index in range(65))
    low = max(-10.0, coarse[1] - 1.0)
    high = min(26.0, coarse[1] + 1.0)
    ratio = (math.sqrt(5.0) - 1.0) / 2.0
    left = high - ratio * (high - low)
    right = low + ratio * (high - low)
    f_left = objective(left)
    f_right = objective(right)
    for _ in range(38):
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


def fit_reaction_model(model: str, runs: list[dict]) -> tuple[float, float, float]:
    coarse = []
    for activation_energy in range(40, 401, 10):
        error, log_prefactor = _best_log_prefactor(model, runs, float(activation_energy))
        coarse.append((error, float(activation_energy), log_prefactor))
    seeds = sorted(coarse)[:4]
    lower = max(30.0, min(item[1] for item in seeds) - 12.0)
    upper = min(430.0, max(item[1] for item in seeds) + 12.0)
    ratio = (math.sqrt(5.0) - 1.0) / 2.0
    left = upper - ratio * (upper - lower)
    right = lower + ratio * (upper - lower)
    left_error, left_log = _best_log_prefactor(model, runs, left)
    right_error, right_log = _best_log_prefactor(model, runs, right)
    for _ in range(30):
        if left_error < right_error:
            upper = right
            right = left
            right_error = left_error
            right_log = left_log
            left = upper - ratio * (upper - lower)
            left_error, left_log = _best_log_prefactor(model, runs, left)
        else:
            lower = left
            left = right
            left_error = right_error
            left_log = right_log
            right = lower + ratio * (upper - lower)
            right_error, right_log = _best_log_prefactor(model, runs, right)
    candidates = list(coarse)
    candidates.extend([(left_error, left, left_log), (right_error, right, right_log)])
    error, activation_energy, log_prefactor = min(candidates)
    return float(activation_energy), float(log_prefactor), float(error)


def normalize_runs(runs: list[dict]) -> dict:
    result = []
    for run in runs:
        initial = float(run["initial_mass_mg"])
        residual = float(run["residual_mass_mg"])
        removable = initial - residual
        if not removable > 0.0:
            raise ValueError("initial_mass_mg must exceed residual_mass_mg")
        conversion = [
            (initial - float(point["mass_mg"])) / removable for point in run["points"]
        ]
        result.append({"run_id": run.get("run_id", ""), "conversion": conversion})
    return {"runs": result}


def identify_runs(runs: list[dict]) -> dict:
    if len(runs) < 2:
        raise ValueError("identify requires at least two heating-rate runs")
    for run in runs:
        beta = float(run["beta_k_per_min"])
        temperatures = [float(value) for value in run["temperature_c"]]
        conversions = run["conversion"]
        if beta <= 0.0:
            raise ValueError("beta_k_per_min must be positive")
        if len(temperatures) < 2 or len(temperatures) != len(conversions):
            raise ValueError("temperature_c and conversion must have equal nonzero length")
        if any(right <= left for left, right in zip(temperatures, temperatures[1:])):
            raise ValueError("temperature_c must be strictly increasing")
    fits = {}
    for model in MODEL_NAMES:
        activation_energy, log_prefactor, error = fit_reaction_model(model, runs)
        fits[model] = {
            "activation_energy_kj_mol": activation_energy,
            "log10_prefactor_s": log_prefactor,
            "sse": error,
        }
    selected = min(MODEL_NAMES, key=lambda name: fits[name]["sse"])
    best = fits[selected]
    return {
        "model": selected,
        "activation_energy_kj_mol": best["activation_energy_kj_mol"],
        "log10_prefactor_s": best["log10_prefactor_s"],
        "model_scores": {name: fits[name]["sse"] for name in MODEL_NAMES},
    }


def predict_conversion(request: dict) -> dict:
    model = request["model"]
    if model not in MODEL_NAMES:
        raise ValueError(f"unknown reaction model: {model}")
    activation_energy = float(request["activation_energy_kj_mol"])
    prefactor = 10.0 ** float(request["log10_prefactor_s"])
    schedule = request["schedule"]
    kind = schedule["type"]
    if kind == "linear":
        start_k = float(schedule["start_temperature_c"]) + 273.15
        end_k = float(schedule["end_temperature_c"]) + 273.15
        beta_k_s = float(schedule["beta_k_per_min"]) / 60.0
        if beta_k_s <= 0.0 or end_k <= start_k:
            raise ValueError("linear schedule requires beta_k_per_min > 0 and end > start")
        extent = prefactor * _temperature_integral(start_k, end_k, activation_energy) / beta_k_s
    elif kind == "isothermal":
        temperature_k = float(schedule["temperature_c"]) + 273.15
        time_s = float(schedule["time_min"]) * 60.0
        if time_s < 0.0:
            raise ValueError("time_min must be nonnegative")
        extent = prefactor * _arrhenius(temperature_k, activation_energy) * time_s
    else:
        raise ValueError(f"unknown schedule type: {kind}")
    return {"conversion": float(inverse_extent(model, extent))}


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
            request = json.load(handle)
        result = dispatch(request)
    except (OSError, json.JSONDecodeError, KeyError, TypeError, ValueError, OverflowError) as exc:
        print(f"input error: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(result, ensure_ascii=False, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
