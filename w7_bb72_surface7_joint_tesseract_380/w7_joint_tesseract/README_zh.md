# Weight-7 / BB72 / surface d=7：同一 shot 的联合 X/Z 电路级模拟

本程序比较 weight-7 `[[56,12,7]]`、IBM BB `[[72,12,6]]` 和 rotated surface `[[49,1,7]]`，采用完整含噪电路采样，以及 Tesseract 的联合 detector error model（DEM）解码。每个 shot 同时产生 X、Z 错误及相关测量记录；主结果是该 shot 中任意逻辑错误或解码失败的概率。不是把两次独立 memory 实验的统计量相乘。

## 安装与运行

建议使用 Linux 服务器，在解压后的本目录执行。安装脚本会建立 `.venv` 并安装固定版本的 Stim 和已验证的 Tesseract：

```bash
bash install.sh
source .venv/bin/activate
python main.py --smoke
```

完成快速验证后正式运行：

```bash
bash run_linux.sh --workers 376
```

若需要使用全部 380 核：

```bash
bash run_linux.sh --workers 380
```

后台运行并保存日志：

```bash
nohup bash run_linux.sh --workers 376 > run.log 2>&1 &
```

查看日志：

```bash
tail -f run.log
```

程序按物理错误率由高到低处理各点，在当前点内分配并行 shot 批次；一个 worker 使用一个 CPU 线程。376 个 worker 为系统保留约 4 核。实际内存占用与每个进程的 Tesseract 搜索队列有关，启动时会实测解码器内存，并按 CPU affinity 和可用内存限制并发，日志显示实际 worker 数。初始批次自动缩小，避免高错误率点一次提交过多 shots。

## 默认设置

| 项目 | 默认值 |
|---|---|
| 编码 | `w7,bb72,surface7` |
| noisy rounds | 按静态码距 `R=d`：weight-7 为 7，BB72 为 6，surface 为 7 |
| 收尾 | 2 个无噪 syndrome cycles；不计入 R |
| 物理错误率 | `0.006, 0.005, 0.004, 0.0035, 0.003, 0.0025, 0.002, 0.0015, 0.001` |
| CNOT | 理想门后以 p 概率施加 15 种非平凡二比特 Pauli 之一 |
| 初始化 | 以 p 概率初始化到相反本征态 |
| 测量 | 以 p 概率翻转读出 |
| idle | `IDLE_SCALE=0`，即无 idle noise |
| 解码器 | Tesseract，联合 X/Z DEM，`det_beam=2` |
| worker | 376，可改成 380 |
| 最少 shots/点 | 1,000 |
| 停止精度 | 联合失败概率 P 的相对二项标准差 `sqrt(P(1-P)/N)/P ≤10%`，且已有失败样本 |
| shots 上限/点 | 100,000,000 |
| 运行时间上限 | 无 |
| 随机数 | 固定任务种子及独立批次种子，保存批次记录以支持续跑 |

`det_beam=2` 固定启发式剪枝阈值，`beam_climbing=False`，`num_det_orders=20`，`pqlimit=200000`。它不保证得到最优解码结果；队列上限也不是时间上限，少数困难 shots 可能明显慢于平均。父进程每 30 秒打印运行状态。程序另外保存低置信结果与无有效解码结果的计数；不丢弃这些 shots。

初始化和读出使用原生 X/Z 基操作，即 `RX/MX` 和 `R/M`；surface circuit 也使用这一约定。没有把 X 基操作额外拆成带独立噪声的 Hadamard 门。idle 位置按 IBM 的 data-qubit 时序定义；默认系数为 0，因此不会增加任何 idle fault。

可只运行一个码或选定几个点，例如：

```bash
python main.py --codes w7 --p 0.006,0.004,0.002 --workers 376
```

在相同实验设置下提高精度并继续已有结果：

```bash
python main.py --workers 376 --relative-error 0.05 --max-shots 300000000
```

也可编辑 `campaign_config.py` 中的 `RELATIVE_ERROR`、`MIN_SHOTS`、`MAX_SHOTS`、`WORKERS` 等运行控制项。修改 decoder、rounds、噪声、码或 schedule 会改变实验身份，不能把改变前后的计数合并。调整 worker 数和 `--batch-size` 不改变实验身份，也可继续已保存的批次。增加 shots 上限或收紧停止精度属于同一实验的续跑。

本版本将默认轮数由统一 R=8 改为 R=d，并将相对 1σ 目标从 20% 改为 10%。旧 R=8 的结果仍保留，但因电路轮数不同，不会与新结果合并。建议新版本使用新的结果目录（或解压到新目录）；若在原目录运行，可添加 `--output results_Rd`。仅改变停止精度时仍支持原实验续跑。

## IBM 电路与收尾边界

BB72 使用 `ell=m=6`、`A=x^3+y+y^2`、`B=y^3+x+x^2`。其 CNOT 顺序参考 IBM 提供的 `BivariateBicycleCodes` 源码：

```text
sX = ['idle', 1, 4, 3, 5, 0, 2]
sZ = [3, 5, 0, 1, 2, 4, 'idle']
```

每周期是 7 层 CNOT，加上与初始化、测量交错的第 8 个时间步。X ancilla 在周期开始初始化；Z ancilla 在周期结束初始化，供下一个周期使用。第一个 noisy cycle 之前的 Z ancilla 是理想状态。最后一个 noisy cycle 末尾发生的 Z 初始化错误仍然保留，并进入后面的无噪收尾周期。实现使用两个无噪收尾周期，与用户提供的 IBM `decoder_setup.py` 和 `decoder_run.py` 一致。

detector 使用相邻轮 syndrome 的差分，并包括初始参考和最终无噪边界。不能遗漏 terminal syndrome，也不能把每轮读出当成彼此独立的 detector。

IBM 原论文采用 idle error=p，本任务按用户设置采用 idle error=0；本程序的数值因此不应被称为原论文噪声设置的直接复现。另外，本程序使用联合 Tesseract DEM，保留 X/Z detector 的相关信息；IBM 原代码使用同一条物理错误历史，但分别执行 X 和 Z 的 BP+OSD。两者的解码器不同。

已在本程序实际使用的 `R=d`（分别 7、6、7）、`idle=0`、两个无噪收尾周期下严格核验：

| 码 | 静态 d | CNOT 层/周期 | 时间步/周期 | CNOT 数/周期 | 电路距离 |
|---|---:|---:|---:|---:|---:|
| weight-7 `[[56,12,7]]` | 7 | 8 | 10 | 392 | **d_circ=6** |
| IBM BB `[[72,12,6]]` | 6 | 7 | 8 | 432 | **d_circ=6** |
| rotated surface `[[49,1,7]]` | 7 | 4 | 6 | 168 | **d_circ=7** |

weight-7 使用 `x^2=y^14=1`、`A=1+y+x*y^10`、`B=1+y^3+x*y^11+x*y^2`，`HX=[A|B]`、`HZ=[B^T|A^T]`。从 624 个合法候选中选取了当前 8 层排程；最初 7 层排程的精确距离为 5，保存在 `distance7layer.json`。这是有限范围的排程改进，不宣称全局最优。

weight-7 与 BB72 的严格下界来自完整的低重量故障排除，且各有实际 6 故障逻辑实例给出上界。surface 的下界通过完整 CSS 投影的精确图搜索获得，并用实际 7 故障实例核对上界。距离验证与 Tesseract 的启发式解码性能分开；证明所用投影只是放松约束来计算下界，不是将 Monte Carlo 拆成两次独立实验。

这些结论适用于报告中的排程、轮数和噪声支持。更改相关参数后，程序会隐藏不再适用的距离标签。实际操作、证据和验证方法保存在 `distance_bounds.json`、`schedule_search_summary.json` 及生成的电路文件中。三个静态码的 X/Z distance 已分别精确验证为 7、6、7。`distance_bounds.json` 是绘图的电路距离标注来源。只有下界、上界一致且获得证明时，才标 `d_circ=...`；仅找到某个 fault witness 时只能给出上界。静态 distance 或 Monte Carlo 曲线的斜率都不能替代这一验证。

## 同一 shot 如何得到 X/Z 逻辑错误

对于残余 Pauli `E=X(e_X)Z(e_Z)`，逻辑 X 分量由 `L_Z e_X` 检测，逻辑 Z 分量由 `L_X e_Z` 检测。任意一位为 1 即表示该 sector 的 block failure。物理 depolarizing fault 可以在同一 shot 中同时贡献 X 与 Z 分量。

实现用理想的 encoded Bell reference 记录所有逻辑 Pauli 分量。reference 是模拟和记账手段，无噪声，不参加被比较的硬件开销，也不在真实 data qubit 上同时测量反对易的逻辑 X、Z。Bell 扩展后的相关观测量彼此对易；这与 IBM 源码跟踪经典 Pauli 错误历史后计算逻辑奇偶性具有相同作用。

硬件数只计 data 与 syndrome ancilla：

| 编码 | data | syndrome ancilla | 硬件总数 | k |
|---|---:|---:|---:|---:|
| weight-7 | 56 | 56 | 112 | 12 |
| BB72 | 72 | 72 | 144 | 12 |
| surface d=7，1 patch | 49 | 48 | 97 | 1 |
| surface d=7，12 patches | 588 | 576 | 1,164 | 12 |

计数含义：

| 字段 | 含义 |
|---|---|
| `N` | 总 shots；不是总逻辑比特数 |
| `FX` | 得到有效解码结果后，残余逻辑 X 分量非零的 shots |
| `FZ` | 得到有效解码结果后，残余逻辑 Z 分量非零的 shots |
| `Fboth` | 同一个 shot 同时发生上述 X 和 Z 逻辑错误 |
| `decodefail` | 没有有效解码结果的 shots，与有效解码后的逻辑错误计数分开 |
| `Fjoint` | `FX+FZ-Fboth+decodefail`，主结果使用的操作失败计数 |
| `lowconfidence` | 解码器报告低置信的 shots；不等同于失败，也不丢弃 |
| `logical_X`, `logical_Z` | 各逻辑比特的已知残余错误计数 |

未能返回有效解码时，不能声称存在某个具体 X/Z 残余，因此不把这些 shots 伪造为 `Fboth`。绘图的 sector 曲线使用 `FX+decodefail` 或 `FZ+decodefail` 作为保守的操作失败数；只要 `decodefail>0`，图上会明确标为上界。当 `decodefail=0` 时它们就是通常的 sector failure 曲线。

## 概率、误差棒和 surface 对比

每个点的整段 block failure probability 为 `P=Fjoint/N`。依 IBM 的归一化约定，绘制每 noisy round 的失败率

```text
p_round = 1 - (1-P)^(1/R).
```

小 P 时接近 `P/R`。这是一种归一化约定，不是从最终错误概率严格反演任意时刻的独立逻辑跳变率。

主图提供两种对比：一个编码块的失败率，以及相同 `k=12` 的失败率。后一张图从单个 surface patch 的联合失败率 `P_surface` 推出 12 个独立 patch 的整段失败率：

```text
P_surface_12 = 1 - (1-P_surface)^12
p_surface_12_round = 1 - (1-P_surface)^(12/R).
```

这是独立 patch 假设下的解析换算，不额外模拟 12 倍 shots，也不把单 patch 的样本计数乘以 12。同一 patch 内的 X/Z 相关性已经包含在 `P_surface` 中。

误差棒由原始二项计数的 Wilson 区间、`z=1` 得到，名义覆盖率约 68.27%，再对上下端点做同样的单调概率换算。零失败样本只画该区间的上界和向下箭头；不会当成已知的零逻辑错误率，也不参与伪造低错误率幂律拟合。达到 shots 上限却未达到精度的点仍保留，其状态会记录为未满足目标精度。

## 续跑、原始数据与只重画图

重新运行相同命令即可继续已有结果。配置和电路身份匹配的完成批次会保留；已达到当前停止条件的点会跳过。程序通过输出目录锁阻止两个主进程同时写入。第一次 Ctrl+C 停止提交新任务并等待正在运行的批次保存；再次 Ctrl+C 终止 worker，尚未记账的批次下次会重新计算，已完成批次不会重复累计。

每点主要文件在：

```text
results/points/<code>_p..._R..._<experiment hash>/result.json
```

`result.json` 保存原始计数、实验设置、置信区间、解码状态和批次进度。可随时只读取已有结果重画图；不修改运行状态：

```bash
python plot_results.py --results results
```

若使用了自定义输出目录：

```bash
python plot_results.py --results my_results --out my_plots
```

若结果根目录同时包含同一码、同一 p/R 的多个不同实验，绘图会拒绝静默合并；请为不同实验指定单独的输出目录，或选择对应结果文件绘图。

输出包括 5 组 PDF/PNG：

- `joint_one_block_1sigma`、`joint_one_block_points`：每个编码块的联合失败率。
- `joint_equal_k12_1sigma`、`joint_equal_k12_points`：统一存储 12 个逻辑比特的对比。
- `sectors_X_then_Z_1sigma`：X 在左、Z 在右的同一 shot sector 统计。

`plot_summary.csv` 保存作图数值和原始计数。对于 `patches=12` 的派生行，计数仍然指实际采样的单 patch；只有概率及区间进行了 12-patch 换算。

论文参考：S. Bravyi et al., *High-threshold and low-overhead fault-tolerant quantum memory*, Nature 627, 778–782 (2024), https://doi.org/10.1038/s41586-024-07107-7 。

## 单独复核（不运行生产 Monte Carlo）

```bash
python validate_circuits.py
python static_distance.py
python distance_audit.py --codes w7 bb72 surface7 --seconds 0 --projected-five-seconds 120 --output distance_bounds_recheck.json
```

低重量故障枚举可能使用约 1 GB 内存，并按码顺序执行；不需要启动 380 个距离验证进程。若验证预算耗尽，输出会明确保留当前上下界，不把超时当作证明。随包 `validation_circuits.json` 记录独立经典 Pauli 传播与 Stim 的 517 个故障历史对照，以及对用户上传 IBM 原电路的逐操作核对；`validation_decoder.json` 仅记录小规模功能测试，不能用作生产逻辑错误率数据。
