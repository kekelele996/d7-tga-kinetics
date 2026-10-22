# 参考解说明

## 三个难点

| 编号 | 错误认知 | 正确知识 | 对应测试（tests/scoring.json） |
|---|---|---|---|
| difficulty_1 | 把当前失重除以初始质量，忽略残余质量 | α=(m0-m)/(m0-mf) | `test_complete_run_reaches_unit_conversion`、`test_partial_run_uses_declared_residual_mass`、`test_multiple_runs_preserve_order_and_ids`、`test_first_point_uses_declared_initial_mass_not_sample_maximum` |
| difficulty_2 | 把 K/min 直接当成 K/s，导致秒制 A 偏大 60 倍 | 积分前使用 β_s=β_Kmin/60 | `test_prefactor_is_second_based_for_first_order_case_a` 至 `_d` |
| difficulty_3 | 只保留一级模型，或者只用单一曲线判别，并且预测时总用一级模型 | 对四条曲线联合拟合 F1/F2/A2/R3，比较 SSE 后代入请求模型预测 | `test_model_selects_*`、`test_avrami_linear_schedule_prediction`、`test_contracting_sphere_isothermal_prediction` |

## 参考解文件

- `solution/app/main.py`：完整实现，含质量归一、温度积分、四模型拟合、模型选择、线性和恒温预测。
- `solve.sh`：安装参考实现，并通过 `HTB_SKIP` 精确回退三个难点中的任意一个。

## 关键步骤

1. 用声明的初始质量和残余质量归一化每次实验的质量序列。
2. 把摄氏度转成开尔文，把 K/min 的升温速率转成 K/s。
3. 对任意 Ea，先计算每个温度区间的 Arrhenius 温度积分，再在固定 Ea 下优化 log10(A/s⁻¹)。
4. 对 F1、F2、A2、R3 分别拟合，比较全部升温速率上的残差平方和，选 SSE 最小的模型。
5. 用请求中的模型、Ea、A 与 schedule 计算线性升温或恒温后的转化率；R3 在 g≥1 时按完全转化处理。

## 核心公式与独立复算

- 解析推导：`research/derivations.md` 第 1–5 节。
- 独立数值复核：`python3 -I research/verify_kinetics.py`。
- 两条路径的四类指前因子误差为 0，最大线性化残差约 1.14e-13。
- A2 线性工况：32 段积分的参考实现与 256 段独立复核差值约 6.4e-7，采用容差 0.001。
- R3 恒温工况：独立解析值为 0.3970144875，采用容差 0.001。

## 容差依据

- 质量归一：纯代数，容差 1e-9。
- 秒制指前因子：4e-5 确定性扰动和有限网格造成约 0.23 个对数单位偏移，取 0.32；漏除 60 的典型错误偏差为 1.778，明显可分离。
- 活化能：多曲线网格拟合误差小于 3 kJ/mol；错用反应模型通常偏差 20 kJ/mol 以上。
- 转化率：差分小于 1e-6，取 0.001；忽略模型造成的偏差超过 0.15。

## 权重依据

- gate 测试（不计分）：入口能跑、错误 operation 的退出码与空标准输出。
- difficulty_1 合计 0.34：质量归一决定后续所有动力学参数，权重最高。
- difficulty_2 合计 0.33：秒制换算决定指前因子的物理量纲与外推可信度。
- difficulty_3 合计 0.33：模型选择和外推是最终业务结果，覆盖三个模型辨识与两个预测工况。

## 数据来源

- 真实数据：无。
- 生成数据：四类固体反应模型的已知参数曲线，使用 32 段复合 Simpson 积分和确定性 4e-5 正弦扰动，无随机种子；复现脚本为 `research/verify_kinetics.py`。
- 生成规则与来源说明只位于 `research/`，不进入 `environment/`。

## 答案隔离

- [x] 镜像只复制 `environment/app/`，不含 `tests/`、`solution/`、`.git/`、缓存和备份。
- [x] 环境代码只说明原型限制，不含正确参数、隐藏曲线或测试期望。
- [x] Dockerfile 从原始原型直接构建，没有复制正确版本后再改坏。
