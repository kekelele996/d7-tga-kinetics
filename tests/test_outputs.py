"""Hidden acceptance tests for the TGA decomposition analysis task.

Expected parameters are ground truth from independently generated multirate
curves. The candidate is always executed as runner in a clean environment.
"""

import json
import math
import os
import signal
import subprocess
import tempfile
from pathlib import Path

ENTRY = "/app/main.py"
PYTHON = "/usr/local/bin/python3"
CANDIDATE_USER = "runner"
CANDIDATE_ENV = {
    "PATH": "/usr/local/bin:/usr/bin:/bin",
    "HOME": "/home/runner",
    "LANG": "C.UTF-8",
    "TZ": "UTC",
}
CALL_TIMEOUT_SEC = 45
R = 8.31446261815324


def run_candidate(*args, timeout=CALL_TIMEOUT_SEC):
    """Run the candidate entry point as runner and kill its process group on timeout."""
    as_root = os.getuid() == 0
    proc = subprocess.Popen(
        [PYTHON, ENTRY, *[str(value) for value in args]],
        cwd="/app",
        env=CANDIDATE_ENV,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        start_new_session=True,
        user=CANDIDATE_USER if as_root else None,
        group=CANDIDATE_USER if as_root else None,
        extra_groups=[] if as_root else None,
    )
    try:
        out, err = proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        os.killpg(proc.pid, signal.SIGKILL)
        proc.communicate()
        raise TimeoutError(f"{ENTRY} exceeded {timeout}s")
    return subprocess.CompletedProcess(proc.args, proc.returncode, out, err)


def request_candidate(payload):
    """Write a private request file and return the candidate's parsed JSON output."""
    folder = Path(tempfile.mkdtemp(prefix="tga-in-"))
    folder.chmod(0o755)
    path = folder / "request.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    path.chmod(0o644)
    proc = run_candidate(path)
    assert proc.returncode == 0, f"exit={proc.returncode} stderr={proc.stderr[-700:]}"
    return json.loads(proc.stdout)


def request_failure(payload):
    """Run an invalid request and return the completed process."""
    folder = Path(tempfile.mkdtemp(prefix="tga-bad-"))
    folder.chmod(0o755)
    path = folder / "request.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    path.chmod(0o644)
    return run_candidate(path)


def clamp_unit(value):
    return min(max(float(value), 1.0e-14), 1.0 - 1.0e-14)


def extent(model, alpha):
    a = clamp_unit(alpha)
    if model == "nth_order_1":
        return -math.log1p(-a)
    if model == "nth_order_2":
        return a / (1.0 - a)
    if model == "avrami_2":
        return math.sqrt(-math.log1p(-a))
    if model == "contracting_sphere":
        return 1.0 - (1.0 - a) ** (1.0 / 3.0)
    raise ValueError(model)


def inverse_extent(model, value):
    z = max(float(value), 0.0)
    if model == "nth_order_1":
        return 1.0 - math.exp(-z)
    if model == "nth_order_2":
        return z / (1.0 + z)
    if model == "avrami_2":
        return 1.0 - math.exp(-z * z)
    if model == "contracting_sphere":
        return 1.0 - (1.0 - min(z, 1.0)) ** 3
    raise ValueError(model)


def thermal_integral(start_c, end_c, activation_energy, intervals=32):
    start_k = start_c + 273.15
    end_k = end_c + 273.15
    step = (end_k - start_k) / intervals
    total = math.exp(-activation_energy * 1000.0 / (R * start_k))
    total += math.exp(-activation_energy * 1000.0 / (R * end_k))
    for index in range(1, intervals):
        temperature = start_k + index * step
        weight = 4.0 if index % 2 else 2.0
        total += weight * math.exp(-activation_energy * 1000.0 / (R * temperature))
    return total * step / 3.0


def multirate_runs(
    model,
    activation_energy,
    log_prefactor,
    betas,
    start_c,
    end_c,
    initial_conversion=0.0,
):
    runs = []
    points = 160
    for run_index, beta in enumerate(betas):
        temperatures = [
            start_c + (end_c - start_c) * index / points for index in range(points + 1)
        ]
        cumulative = [0.0]
        for index in range(1, len(temperatures)):
            cumulative.append(
                cumulative[-1]
                + thermal_integral(temperatures[index - 1], temperatures[index], activation_energy)
            )
        selected_temperature = []
        selected_conversion = []
        for index, (temperature, integral) in enumerate(zip(temperatures, cumulative)):
            alpha = inverse_extent(
                model,
                extent(model, initial_conversion)
                + 10.0 ** log_prefactor * integral / (beta / 60.0),
            )
            alpha += 4.0e-5 * math.sin((run_index + 1) * 0.71 + index * 1.3)
            if 0.01 < alpha < 0.99 and index % 2 == 0:
                selected_temperature.append(temperature)
                selected_conversion.append(alpha)
        runs.append(
            {
                "run_id": f"r{run_index + 1}",
                "beta_k_per_min": beta,
                "temperature_c": selected_temperature,
                "conversion": selected_conversion,
                "initial_conversion": initial_conversion,
            }
        )
    return runs


def identify_case(case):
    model, activation_energy, log_prefactor, betas, start_c, end_c = case[:6]
    initial_conversion = case[6] if len(case) > 6 else 0.0
    return request_candidate(
        {
            "operation": "identify",
            "runs": multirate_runs(
                model,
                activation_energy,
                log_prefactor,
                betas,
                start_c,
                end_c,
                initial_conversion,
            ),
        }
    )


def normalize_case(runs):
    return request_candidate({"operation": "normalize", "runs": runs})


def close_list(got, want, tolerance=1.0e-9):
    return len(got) == len(want) and all(
        abs(float(a) - b) <= tolerance for a, b in zip(got, want)
    )


# ---------- gate: entry and response contract ----------

def test_normalize_response_shape():
    """「返回 `{"runs": [{"run_id": "...", "conversion": [...]}, ...]}`」"""
    got = normalize_case(
        [
            {
                "run_id": "g1",
                "initial_mass_mg": 10.0,
                "residual_mass_mg": 2.0,
                "points": [{"time_min": 0, "temperature_c": 25, "mass_mg": 10.0}],
            }
        ]
    )
    assert set(got) == {"runs"}
    assert got["runs"][0]["run_id"] == "g1"
    assert isinstance(got["runs"][0]["conversion"], list)


def test_bad_operation_exits_2_without_stdout():
    """「JSON 读不了、`operation` 不认识、模型名不认识或字段缺失时，退出码为 2，标准输出留空。」"""
    folder = Path(tempfile.mkdtemp(prefix="tga-bad-"))
    folder.chmod(0o755)
    path = folder / "request.json"
    path.write_text('{"operation":"unknown"}', encoding="utf-8")
    path.chmod(0o644)
    proc = run_candidate(path)
    assert proc.returncode == 2
    assert proc.stdout == ""


# ---------- difficulty 1: residual-mass conversion baseline ----------

def test_complete_run_reaches_unit_conversion():
    """「转化率在初始质量时为 0，到残余质量时为 1。」"""
    got = normalize_case(
        [
            {
                "run_id": "complete",
                "initial_mass_mg": 12.0,
                "residual_mass_mg": 3.0,
                "points": [
                    {"time_min": i, "temperature_c": 30 + 10 * i, "mass_mg": mass}
                    for i, mass in enumerate((12.0, 10.5, 9.0, 7.5, 6.0, 4.5, 3.0))
                ],
            }
        ]
    )
    assert close_list(
        got["runs"][0]["conversion"], [0.0, 1 / 6, 1 / 3, 0.5, 2 / 3, 5 / 6, 1.0]
    )


def test_partial_run_uses_declared_residual_mass():
    """「转化率在初始质量时为 0，到残余质量时为 1。」"""
    got = normalize_case(
        [
            {
                "run_id": "partial",
                "initial_mass_mg": 20.0,
                "residual_mass_mg": 8.0,
                "points": [
                    {"time_min": i, "temperature_c": 100 + 20 * i, "mass_mg": mass}
                    for i, mass in enumerate((20.0, 18.0, 15.5, 13.0, 10.0))
                ],
            }
        ]
    )
    assert close_list(got["runs"][0]["conversion"], [0.0, 1 / 6, 0.375, 7 / 12, 5 / 6])


def test_multiple_runs_preserve_order_and_ids():
    """「输出顺序按输入来，`conversion` 跟同一次实验的 `points` 对齐。」"""
    got = normalize_case(
        [
            {
                "run_id": "z-run",
                "initial_mass_mg": 10.0,
                "residual_mass_mg": 2.0,
                "points": [
                    {"time_min": 0, "temperature_c": 20, "mass_mg": 10.0},
                    {"time_min": 1, "temperature_c": 30, "mass_mg": 6.0},
                ],
            },
            {
                "run_id": "a-run",
                "initial_mass_mg": 20.0,
                "residual_mass_mg": 5.0,
                "points": [
                    {"time_min": 0, "temperature_c": 20, "mass_mg": 20.0},
                    {"time_min": 1, "temperature_c": 30, "mass_mg": 12.5},
                ],
            },
        ]
    )
    assert [item["run_id"] for item in got["runs"]] == ["z-run", "a-run"]
    assert close_list(got["runs"][0]["conversion"], [0.0, 0.5])
    assert close_list(got["runs"][1]["conversion"], [0.0, 0.5])


def test_first_point_uses_declared_initial_mass_not_sample_maximum():
    """「输出顺序按输入来，`conversion` 跟同一次实验的 `points` 对齐。」"""
    got = normalize_case(
        [
            {
                "run_id": "offset",
                "initial_mass_mg": 15.0,
                "residual_mass_mg": 4.0,
                "points": [
                    {"time_min": 0, "temperature_c": 40, "mass_mg": 14.9},
                    {"time_min": 1, "temperature_c": 60, "mass_mg": 12.0},
                    {"time_min": 2, "temperature_c": 80, "mass_mg": 9.0},
                    {"time_min": 3, "temperature_c": 100, "mass_mg": 6.0},
                    {"time_min": 4, "temperature_c": 120, "mass_mg": 4.0},
                ],
            }
        ]
    )
    assert close_list(
        got["runs"][0]["conversion"],
        [0.1 / 11.0, 3.0 / 11.0, 6.0 / 11.0, 9.0 / 11.0, 1.0],
    )


def test_normalize_preserves_values_outside_endpoints():
    """「称量抖动让质量略越过两个端点时，保留带符号的数值，不要截断到 0 或 1。」"""
    got = normalize_case(
        [
            {
                "run_id": "drift",
                "initial_mass_mg": 10.0,
                "residual_mass_mg": 2.0,
                "points": [
                    {"time_min": i, "temperature_c": 20 + 20 * i, "mass_mg": mass}
                    for i, mass in enumerate((10.2, 10.0, 6.0, 2.0, 1.8))
                ],
            }
        ]
    )
    assert close_list(got["runs"][0]["conversion"], [-0.025, 0.0, 0.5, 1.0, 1.025])


def test_equal_initial_and_residual_mass_is_rejected():
    """「`initial_mass_mg` 必须大于 `residual_mass_mg`，否则按坏请求处理。」"""
    proc = request_failure(
        {
            "operation": "normalize",
            "runs": [
                {
                    "run_id": "bad-baseline",
                    "initial_mass_mg": 5.0,
                    "residual_mass_mg": 5.0,
                    "points": [{"time_min": 0, "temperature_c": 25, "mass_mg": 5.0}],
                }
            ],
        }
    )
    assert proc.returncode == 2
    assert proc.stdout == ""


# ---------- difficulty 2: second-based prefactor from K/min heating rates ----------

def test_prefactor_is_second_based_for_first_order_case_a():
    """「秒制要算进指前因子，不能把分钟制数值直接交上去。」"""
    got = identify_case(("nth_order_1", 115.0, 9.2, (1.5, 3.0, 8.0, 18.0), 120.0, 340.0))
    assert abs(got["activation_energy_kj_mol"] - 115.0) <= 3.0
    assert abs(got["log10_prefactor_s"] - 9.2) <= 0.32


def test_prefactor_is_second_based_for_first_order_case_b():
    """「秒制要算进指前因子，不能把分钟制数值直接交上去。」"""
    got = identify_case(("nth_order_1", 168.0, 12.8, (2.0, 5.0, 10.0, 20.0), 200.0, 430.0))
    assert abs(got["activation_energy_kj_mol"] - 168.0) <= 3.0
    assert abs(got["log10_prefactor_s"] - 12.8) <= 0.32


def test_prefactor_is_second_based_for_first_order_case_c():
    """「秒制要算进指前因子，不能把分钟制数值直接交上去。」"""
    got = identify_case(("nth_order_1", 88.0, 6.8, (0.8, 2.0, 6.0, 15.0), 100.0, 320.0))
    assert abs(got["activation_energy_kj_mol"] - 88.0) <= 3.0
    assert abs(got["log10_prefactor_s"] - 6.8) <= 0.32


def test_prefactor_is_second_based_for_first_order_case_d():
    """「秒制要算进指前因子，不能把分钟制数值直接交上去。」"""
    got = identify_case(("nth_order_1", 205.0, 15.3, (2.5, 6.0, 12.0, 24.0), 240.0, 450.0))
    assert abs(got["activation_energy_kj_mol"] - 205.0) <= 3.0
    assert abs(got["log10_prefactor_s"] - 15.3) <= 0.32


def test_identify_rejects_non_monotonic_temperature():
    """「温度倒退的记录不收。」"""
    proc = request_failure(
        {
            "operation": "identify",
            "runs": [
                {
                    "run_id": "bad-order",
                    "beta_k_per_min": 5.0,
                    "temperature_c": [100.0, 120.0, 110.0, 140.0],
                    "conversion": [0.1, 0.2, 0.3, 0.4],
                },
                {
                    "run_id": "ok-order",
                    "beta_k_per_min": 10.0,
                    "temperature_c": [100.0, 120.0, 140.0, 160.0],
                    "conversion": [0.1, 0.2, 0.3, 0.4],
                },
            ],
        }
    )
    assert proc.returncode == 2
    assert proc.stdout == ""


def test_prefactor_with_nonzero_initial_conversion():
    """「有些样品送到仪器前已经反应了一部分，`initial_conversion` 就是那个起点。」"""
    got = identify_case(
        ("nth_order_1", 132.0, 10.1, (1.2, 3.5, 9.0, 20.0), 130.0, 360.0, 0.18)
    )
    assert abs(got["activation_energy_kj_mol"] - 132.0) <= 3.5
    assert abs(got["log10_prefactor_s"] - 10.1) <= 0.32


# ---------- difficulty 3: model discrimination and schedule prediction ----------

def test_model_selects_second_order_reaction():
    """「同一批样品的所有升温速率曲线要一起参与模型判断。」"""
    got = identify_case(("nth_order_2", 142.0, 10.5, (2.0, 5.0, 10.0, 20.0), 190.0, 420.0))
    assert got["model"] == "nth_order_2"


def test_model_selects_avrami_erofeev_reaction():
    """「同一批样品的所有升温速率曲线要一起参与模型判断。」"""
    got = identify_case(("avrami_2", 125.0, 8.7, (2.0, 5.0, 10.0, 20.0), 190.0, 430.0))
    assert got["model"] == "avrami_2"


def test_model_selects_contracting_sphere_reaction():
    """「同一批样品的所有升温速率曲线要一起参与模型判断。」"""
    got = identify_case(("contracting_sphere", 155.0, 9.7, (2.0, 5.0, 10.0, 20.0), 270.0, 500.0))
    assert got["model"] == "contracting_sphere"


def test_avrami_linear_schedule_prediction():
    """「升温段从起始温度一路走到终止温度」"""
    got = request_candidate(
        {
            "operation": "predict",
            "model": "avrami_2",
            "activation_energy_kj_mol": 125.0,
            "log10_prefactor_s": 8.7,
            "schedule": {
                "type": "linear",
                "start_temperature_c": 100.0,
                "end_temperature_c": 300.0,
                "beta_k_per_min": 4.2,
            },
        }
    )
    assert abs(got["conversion"] - 0.294921) <= 0.001


def test_contracting_sphere_isothermal_prediction():
    """「保温段从零转化率开始计时，`time_min` 不能为负。」"""
    got = request_candidate(
        {
            "operation": "predict",
            "model": "contracting_sphere",
            "activation_energy_kj_mol": 155.0,
            "log10_prefactor_s": 9.7,
            "schedule": {"type": "isothermal", "temperature_c": 370.0, "time_min": 2.0},
        }
    )
    assert abs(got["conversion"] - 0.397014) <= 0.001

def test_predict_rejects_nonpositive_linear_rate():
    """「`beta_k_per_min` 必须为正且终止温度高于起始温度。」"""
    proc = request_failure(
        {
            "operation": "predict",
            "model": "nth_order_1",
            "activation_energy_kj_mol": 100.0,
            "log10_prefactor_s": 8.0,
            "schedule": {
                "type": "linear",
                "start_temperature_c": 100.0,
                "end_temperature_c": 300.0,
                "beta_k_per_min": 0.0,
            },
        }
    )
    assert proc.returncode == 2
    assert proc.stdout == ""


def test_predict_rejects_negative_isothermal_time():
    """「保温段从零转化率开始计时，`time_min` 不能为负。」"""
    proc = request_failure(
        {
            "operation": "predict",
            "model": "nth_order_1",
            "activation_energy_kj_mol": 100.0,
            "log10_prefactor_s": 8.0,
            "schedule": {"type": "isothermal", "temperature_c": 150.0, "time_min": -1.0},
        }
    )
    assert proc.returncode == 2
    assert proc.stdout == ""

