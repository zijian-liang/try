# Weight-7 电路调度对照：先筛电路距离，再做联合 X/Z 仿真

本程序用于检查：同一个 [[56,12,7]] 码，只改变 stabilizer extraction 的 CNOT 调度，是否会使当前 Tesseract beam=2 的逻辑错误率曲线上下移动。它包含原先使用的8层电路 `w7_baseline`，并生成新的8层、9层及少量14层调度。

只跑物理错误率 **p=0.003（0.3%）和 p=0.0025（0.25%）**。X/Z 来自同一个物理 shot，使用完整 correlated DEM 一次解码。不同调度独立采样，保持同样的码、噪声模型、轮数和解码器参数。

## 默认设置

| 参数 | 默认值 / 作用 |
|---|---|
| 码 | weight-7 [[56,12,7]]；56 data + 56 syndrome ancilla = 112 物理比特 |
| 候选数 | 64；8层30个（含原基线）、9层30个、14层4个 |
| 仿真电路数上限 | 12，包含基线；从距离筛选通过者中按调度类别选择，不根据仿真结果预挑 |
| 噪声轮数 | R=7，随后2轮理想收尾；沿用原IBM边界约定 |
| CNOT 数 | 所有调度每轮392个；整段2744个有噪声 CNOT |
| 准备 / 测量 | 各392个有噪声位置；native RX/MX、R/M |
| 噪声 | CNOT 后总概率 p 的双比特 depolarizing；准备/测量翻转概率 p；idle=0 |
| 解码 | Tesseract 0.1.1.dev20260910235247；beam=2；beam_climbing=False |
| 其余解码参数 | num_det_orders=20；no_revisit_dets=True；merge_errors=True；pqlimit=200000；sparsify_errors=False；decoder seed=2384753 |
| 仿真进程 | 请求380；根据 CPU affinity、资源限制和内存估计向下调整；每个进程单线程 |
| 距离进程 | 请求32；每个证明需要较多内存，额外自动限制并发 |
| 距离搜索预算 | 每候选上界搜索最多30秒；下界证明总预算300秒；超时保留未证实状态 |
| 停止条件 | 每点至少1000 shots，联合失败数>0，相对二项近似1σ≤10%；最多1亿 shots；仿真无时间上限 |
| 判断停止的时机 | 一整批预先安排的并行任务全部完成后；运行中显示 provisional |
| 顺序 | 全部候选先完成距离筛选，然后 p=0.003 → 0.0025；每个 p 内先基线，再其他调度 |

## 距离筛选的含义

`upper_bound >= 6` 本身不能保证真实电路距离至少为6。例如，`3 <= d_circ <= 7` 的真实距离仍可能为3、4或5。程序因此始终分别报告 lower_bound、upper_bound、exact。

默认执行以下流程：

1. 检查每条 check-data 边恰好出现一次、无硬件冲突、交错顺序正确，并验证无噪声情况下 detector 和24个逻辑观测量确定为0。
2. 在 **R=7、idle=0、p=0.003** 的实际电路上计算完整 DEM。相同的带 detector/logical 标签的 DEM 去重；不声称解决了所有图同构。
3. 搜索不可检测逻辑故障的物理 witness。找到5个及以下故障位置的 witness，直接排除。
4. 对其余方案做下界证明：在保持逻辑标签的CSS投影中，完整排除权重≤4，再在适用的已证前提下完整排除权重5。只有证明完成，才报告 `lower_bound >= 6`。超时不是证明。
5. `lower_bound >= 6` 且有 `upper_bound >= 6` 的物理 witness 才默认通过。`[6,6]` 标为 exact 6；`[6,7]` 标为区间，绝不自动写成7。

每个上界 witness 都映射回实际允许的物理故障，并核对 detector=0、logical≠0。下界来自完整物理故障模型的放松问题，完整排除才能用于证明。由于两个仿真 p 均为正且故障支持相同，距离筛查在 p=0.003 做一次即可用于这两个点。

如果你希望严格按“只要求搜索所得上界≥6”做更宽松的探索，可以在**新的输出目录**增加 `--allow-uncertified`。这时未证明 lower>=6 的候选会显式标为未认证；不能据它们的曲线宣称所有被比较电路都有 d_circ>=6。建议首先使用默认严格模式。

8/9层候选来自有限的 translation-invariant 调度搜索。历史搜索只提供调度种子，旧的R2/R8上界不会继承为本次距离结论。14层分开执行 X/Z 的方案也必须通过距离筛选，才会进入仿真。此包不保证全局最优调度，也不保证一定找到 d_circ=7。

## 在你已有的 Python 3.12 环境中运行

把压缩包解压，进入包含 `main.py` 的目录。你之前已成功安装依赖的解释器可以直接复用，不需要重新安装 Tesseract。例如：

```bash
cd /你的解压路径/w7_schedule_scan_380

QEC_PYTHON=/home/zijianliang/Mycode/w7_bb72_surface7_joint_tesseract_380/w7_joint_tesseract/.venv-qec312/bin/python

"$QEC_PYTHON" -c 'import stim; from tesseract_decoder import tesseract; print("Dependencies OK")'

nohup "$QEC_PYTHON" -u main.py --workers 380 > scan.log 2>&1 &
echo "PID: $!"
```

请按你服务器上的实际位置修改目录。如果原仿真仍占用380核，等其退出后再启动此任务，或者为两个任务分别分配核数。输出默认位于这个新目录的 `scan_results/`，不读取或合并旧仿真的结果。

查看进度：

```bash
tail -f scan.log
```

此时 Ctrl+C 仅退出 `tail`。要请求后台主程序停止，可对上面打印的主进程 PID 发送 SIGINT：

```bash
kill -INT 你的主进程PID
```

仿真阶段第一次 SIGINT 完成当前已安排的整批任务并保存后退出；第二次 SIGINT 终止工作进程，未计入的批次在下次运行时重做。重新执行相同启动命令即可续跑。不要同时向同一个输出目录启动两个主程序，程序有写锁检查。

## 常用调整

扩大候选集合和实际对照数量：

```bash
"$QEC_PYTHON" -u main.py --workers 380 --candidates 128 --max-circuits 24
```

只生成调度并查看清单：

```bash
"$QEC_PYTHON" -u main.py --stage generate --candidates 64
```

只做距离筛查：

```bash
"$QEC_PYTHON" -u main.py --stage screen --distance-workers 32 --proof-seconds 300
```

下界证明超时，增加预算重试：

```bash
"$QEC_PYTHON" -u main.py --stage screen --proof-seconds 1200 --retry-inconclusive
```

然后进入仿真：

```bash
"$QEC_PYTHON" -u main.py --stage simulate --workers 380 --max-circuits 12
```

仿真入口仍会核查距离身份并复用完整报告，避免误用另一条电路的证书。如果合格方案少于指定上限，会跑现有合格方案；如果只有基线，则停止并提示增加候选数或证明预算，避免把单条曲线当成电路对照。

只跑其中一个指定 p：`--p 0.003` 或 `--p 0.0025`。程序拒绝这两个值以外的 p。

提高精度继续累积：

```bash
"$QEC_PYTHON" -u main.py --stage simulate --workers 380 --relative-error 0.05
```

统计精度、worker 数和 batch size 不进入物理实验身份，允许续跑；电路、p、解码器设置/版本、抽样种子不同则不合并。

## 输出与绘图

- `scan_results/candidates.json`：全部调度、稳定ID、8/9/14层类别和去重信息。
- `scan_results/distance/<ID>.json`：每条调度的距离上/下界、证明状态、物理 witness 及身份信息。
- `scan_results/screening_summary.json`：距离筛查汇总，包括通过、淘汰、超时未证实等状态。
- `scan_results/selected_schedules.json`：本次主对照的固定候选名单，含原基线。
- `scan_results/points/<ID>_p..._<hash>/`：实际 circuit.stim、model.dem、config.json、result.json。
- `result.json` 保留 N、Fjoint、FX、FZ、Fboth、decodefail、lowconfidence 和每逻辑比特的X/Z分量计数；不丢弃困难样本。
- `scan_results/plots/comparison.pdf/png`：两个 p 分别成图，显示每条调度的每块每轮联合错误率和距离标签。
- `scan_results/plots/comparison.csv`：读取到的原始计数和统计量。
- `scan_results/plots/ranking.csv`：描述性排序及相对基线的比值；不是“最优调度”的证明。

程序结束会自动绘图。运行中也可在另一终端执行：

```bash
"$QEC_PYTHON" plot_scan.py --results scan_results
```

默认只画完成的点。`max_shots` 达上限却未达到精度的结果会特别标记；`--include-incomplete` 可以额外显示明确标注的运行中临时结果。绘图不修改仿真数据。

所有被比较的码均 k=12、R=7，主图使用 `1-(1-Fjoint/N)^(1/7)`，即 per block per round；无需通过 per-logical 归一化来比较这些调度。误差棒为 Wilson z=1 区间的单调变换，是常规近似统计表示，不宣称在自适应停止和多重比较下具有严格覆盖保证。

## 如何判断是否来自电路

先看原基线的新结果是否与此前曲线在统计误差内一致。若同码、同噪声、同解码参数下，另一条已认证调度在两个 p 都明显更低，说明调度确实会影响当前解码配置的表现。其原因可能涉及故障传播方式、致错故障组合数量及解码器对新 DEM 的处理，单凭这个对照不能进一步把这些机制全部分离。

比较了很多调度后选出最低点，会有挑选产生的偏差。可把初步优胜调度与基线放在新目录，用独立抽样种子做确认，例如：

```bash
"$QEC_PYTHON" -u main.py --workers 380 --output confirm_results \
  --sample-seed 20260922 --only-ids 在这里填写候选ID \
  --min-shots 100000 --max-shots 100000
```

这里候选生成种子 `--seed` 保持默认不变；`--sample-seed` 只改变Monte Carlo抽样。确认实验固定每点100000 shots；若原搜索用了超过64个候选，也应保留对应的 `--candidates` 设置以重新生成相同候选。该命令是可选的后续确认，默认扫描不会自动启动它。

## 新建环境（已有环境可跳过）

已安装 uv 时：

```bash
"$HOME/.local/bin/uv" venv --python 3.12 --seed .venv-qec312
.venv-qec312/bin/python -m pip install -r requirements.txt
bash run_linux.sh --workers 380
```

或 `PYTHON_BIN=/实际路径/python3.12 bash install.sh`。固定版本的 Linux wheel 已核实有 CPython3.12 + x86_64 + glibc>=2.35 的组合；你的服务器此前满足系统条件。

## 本包验证

`validate_scan.py` 独立检查不同深度调度的理想电路、全部15类CNOT Pauli错误、随机多故障历史及最后有噪声轮的Z准备错误，对比独立经典Pauli传播与Stim detector/observable结果。交付的验证报告说明实际测试范围。短shot运行只验证解码、保存和续跑流程，不作为调度性能数据；生产仿真结果由你在服务器运行后生成。

交付前已实际完成前8条候选的距离认证：原基线精确为6；另外7条8/9层调度均证明 `6 <= d_circ <= 7`，未宣称精确为7。完整证明记录、源码哈希与物理 witness 位于 `certificate_examples.json`。默认64条候选中其余方案由服务器运行时逐条筛选。
