# 热重分解动力学分析

树脂样品做热重时，会按几种不同速度升温。实验员想把曲线整理成一组能复用的动力学参数，并估算新批次的升温或保温会烧掉多少。现在 `/app/main.py` 报出的模型和指前因子经常对不上复核单。

## 现状

每次实验记录初始质量、残余质量，以及一串时间、温度、质量读数。样品先恒温吹扫，称量总有小抖动。分析只认从初始质量走到残余质量的那一段，仪器漂移不该算成反应。

入口能读 JSON。最近不同机理的样品都被报成同一种模型；换一套升温速率，指前因子又整体差了数量级。程序只用容器里的标准库，不能联网。

夜班交接时，分析员只拿到一份请求文件，不会回头翻原始仪器日志。程序给出的数值要能直接写进复核单。同一份请求重复运行时，各字段也不能变化。

## 要交付的内容

修好 `/app/main.py`，调用方式是 `python3 /app/main.py <请求文件>`。请求是 JSON 对象，靠 `operation` 选处理方式。

### `normalize`

`runs` 里每次实验带 `run_id`、`initial_mass_mg`、`residual_mass_mg` 和 `points`。`points` 的每个元素有 `time_min`、`temperature_c`、`mass_mg`。

转化率在初始质量时为 0，到残余质量时为 1。中间点按本次实验的总失重量归一。拿 12.0 mg 降到 3.0 mg 的实验来说，读到 7.5 mg 时正好是一半。返回 `{"runs": [{"run_id": "...", "conversion": [...]}, ...]}`。输出顺序照输入走，`conversion` 与同一次实验的 `points` 一一对应。

### `identify`

`runs` 至少有两次实验。每条记录自己的升温速度写在 `beta_k_per_min`，单位 K/min。另外两组是等长的 `temperature_c` 和 `conversion`，从前到后按升温顺序排。

返回对象带四个字段。`model` 只能是 `nth_order_1`、`nth_order_2`、`avrami_2` 或 `contracting_sphere`。`activation_energy_kj_mol` 是表观活化能，单位 kJ/mol。`log10_prefactor_s` 是指前因子按 s⁻¹ 计时后的常用对数。`model_scores` 放四个模型的拟合残差平方和，键用模型名。

秒制要算进指前因子，不能把分钟制数值直接交上去。同一批样品的所有升温速率曲线要一起参与模型判断。

### `predict`

请求带 `model`、`activation_energy_kj_mol`、`log10_prefactor_s` 和 `schedule`。升温段用 `type=linear`，再给 `start_temperature_c`、`end_temperature_c`、`beta_k_per_min`；保温段用 `type=isothermal`，再给 `temperature_c`、`time_min`。返回 `{"conversion": 数值}`。

升温段从起始温度一路走到终止温度，保温段从零转化率开始计时。计算时使用请求给出的模型、活化能和指前因子。

## 约束

- 温度输摄氏度，升温速率输 K/min，时间用分钟，指前因子按 s⁻¹ 输出。
- 只用 Python 标准库。不要联网，也不要留下后台进程。
- 同一请求重复执行，结果要一致。
- JSON 读不了、`operation` 不认识、模型名不认识或字段缺失时，退出码为 2，标准输出留空。
