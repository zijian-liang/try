# Weight-7 高斜率电路：精确距离判定与独立仿真

直接运行 `python main.py`。默认选取之前两点斜率最高的三个未决电路，加原电路对照，重新验证距离，并运行两个物理错误率点。所有可直接修改的参数集中在 **settings.py**。

## 本次选择

只读取你之前日志中最终的 `point_done`，没有使用 provisional 进度。按

\[
q=1-(1-F_{\mathrm{joint}}/N)^{1/(12\times7)},\qquad
\alpha=\frac{\log(q_{0.003}/q_{0.0025})}{\log(0.003/0.0025)}
\]

排序。这里的 q 是联合 block 失败概率的**等效 per logical per round 归一化**，不是逐逻辑比特边际失败率。

| 电路 | 历史两点斜率及名义 1σ | 历史距离范围 |
|---|---:|---:|
| parallel9_a2b5cacfbb84 | 6.6687 ± 1.3864 | [6,7] |
| parallel9_4ec1ca467df4 | 6.0208 ± 1.5513 | [6,7] |
| parallel8_4c4fb6c67ebc | 5.8322 ± 1.2733 | [6,7] |
| w7_baseline（对照） | 5.5930 ± 1.3025 | [6,6] |

这些误差不包含选择多个候选及自适应停止带来的统计影响。排序仅用于分配计算资源，不能推出距离为7或统计显著优于对照。data/ 中附有原始日志、12个实际测量电路的完整 schedule、排名及原始计数；schedule 与原来64候选的确定性生成结果逐一核对过。

**本次实际求解新增结果：** `parallel8_4c4fb6c67ebc` 已找到6个不同物理位置的故障见证，结合重新完成的≤5故障排除证明，精确距离为 **6**。前两名在本地短预算下仍为 `[6,7]`。包内 `data/known_physical_witnesses.json` 保存该见证与baseline见证；运行时会重新核对电路hash、完整DEM和物理位置，验证成功后不再对已知≤6的电路浪费SAT计算。第三名仍保留在独立仿真中作为额外对照。

## 运行

建议直接使用服务器上**已经能成功运行 Tesseract 的同一个 Python 环境**。本包继续使用原版本，不需要重装已有的 Tesseract。

```bash
unzip w7_high_slope_dc7_380.zip
cd w7_high_slope_dc7_380
python3 -m pip install 'python-sat==1.9.dev15'
nohup python3 -u main.py > run.log 2>&1 &
tail -f run.log
```

若你之前使用的是虚拟环境，请将上述所有 `python3` 换成该环境的绝对 Python 路径，或先激活环境。`requirements.txt` 列出了完整依赖，只有在新环境才需要整体安装。

默认结果写入本包内部 `results/`；不会读取或改写旧 simulation 的断点。重复相同命令会续跑。改变随机种子或候选集合时请指定新的 `--output`。

```bash
# 只求解距离，不跑仿真
python3 -u main.py --stage distance

# 距离任务超时后，增加每个 SAT 分片的预算；已有完成任务继续复用
nohup python3 -u main.py --stage distance --sat-seconds 3600 --retry-unknown > distance_retry.log 2>&1 &

# 使用当前距离报告，单独跑独立仿真
nohup python3 -u main.py --stage simulate > simulation.log 2>&1 &

# 将仿真的目标误差改成10%，同一个结果目录会继续积累shots
python3 -u main.py --stage simulate --relative-error 0.10

# 查看/重新生成结果图；不进行仿真
python3 main.py --stage plot

# 扩大到所有11个已测量的未决候选，加baseline；使用独立目录
nohup python3 -u main.py --top 11 --output results_top11 > top11.log 2>&1 &
```

不要在同一结果目录同时启动多个协调进程。进程锁会拒绝这种写入。只求解距离和 simulation 也是先后运行。

## 默认参数，可在 settings.py 直接修改

| 参数 | 默认值 | 含义 |
|---|---:|---|
| TOP_CANDIDATES | 3 | 历史斜率最大的未决候选数，额外加入baseline |
| SIMULATION_WORKERS | 380 | 仿真进程上限，每进程一个CPU线程 |
| DISTANCE_WORKERS | 380 | SAT分片进程上限 |
| SCREENING_WORKERS | 32 | 重新验证小权重下界的进程上限 |
| DISTANCE_SHARDS | 16 | 每个逻辑bit的分片数；每电路最多384个任务 |
| SAT_SECONDS_PER_TASK | 300 | 每个SAT任务的搜索秒数，超时记UNKNOWN |
| SAT_SOLVER | minisat22 | 也可选 glucose3 / glucose4，改变求解器会重新分配任务 |
| SCREENING_PROOF_SECONDS | 300 | 每电路重新验证六故障下界的预算 |
| P_VALUES | [0.003, 0.0025] | 先0.30%，后0.25% |
| RELATIVE_1SIGMA | 0.20 | 联合失败概率的目标相对标准误差 |
| MIN_SHOTS | 1000 | 每点至少的shots数 |
| MAX_SHOTS | 100000000 | 每点shots上限；达到上限不冒充达到精度 |
| SAMPLE_SEED | 20260923 | 与历史筛选数据独立的新种子 |

实际并行数会根据CPU配额和可用内存自动降低。距离筛查、SAT、simulation依次运行，进程数不叠加。380核不意味着每个阶段都恰好能用满380核：初步筛查仅有4个电路，而SAT有足够多的分片。内存估计不是严格峰值上限；其他作业占用大量内存时可手动减小worker数。

## 怎么严格判定是否为7

使用当前完整电路、完整未拆分的DEM和真实单物理故障响应，固定 R=7、idle=0、联合X/Z。历史日志里的 `[6,7]` 不直接充当证明：程序先重新验证小权重下界，并重新构造物理上可实现的7故障上界。

对每个CSS投影，SAT求解如下**不超过六故障**的问题：

\[
H_D e=0\pmod2,\qquad (H_L e)_j=1\pmod2,\qquad |e|\le6.
\]

枚举两种CSS分量、每种分量的12个逻辑bit，按第一个被选中的故障响应索引分片。分片并集覆盖所有非零故障选择，允许分片之间逻辑条件重叠。X/Z投影是完整故障问题的松弛，因此每种投影都排除六故障后，下界对完整相关噪声仍成立；并未把Monte Carlo解码拆成独立CSS解码。

- **全部分片 UNSAT + 实际7故障见证：** 得到严格 `d_circ=7`。
- **发现实际6故障见证 + 已证明没有≤5故障见证：** 得到 `d_circ=6`。
- **某分片超时/出错/未完成：** 保留已证上下界，通常为 `[6,7]`，不把“没找到”写成“不存在”。
- 投影里的SAT解必须映射到完整DEM并核对物理位置，才接受为上界；无法映射的投影解只说明该松弛尚不能证明下界。

由于最低权重求解的复杂度，300秒预算不保证所有候选都能定论。缓存基于电路/DEM/算法/任务身份；增加预算可以重试UNKNOWN任务。求解器返回UNSAT作为精确判定，程序不声称附有独立DRAT证明检查器。

即使距离仍是 `[6,7]` 或已确定为6，默认也继续仿真所选电路，并在图中如实标注。这使得距离验证暂时未完成时，仍能独立检查之前的高斜率现象。

## 仿真与输出

物理模型保持一致：W7 [[56,12,7]]；R=7 noisy rounds及原来的两个ideal closing rounds；CNOT后总概率p的双比特去极化；制备/测量错误率p；idle=0；同一shot的联合X/Z；完整相关DEM；Tesseract固定beam=2。失败数是发生任一逻辑错误的shots数，不是逻辑bit错误数；解码异常单独统计并计入joint失败。

每个采样wave在运行前固定，全部收齐后才判断停止，以避免先完成的易解码shots被优先计入。名义相对误差是 `sqrt((1-Pjoint)/Fjoint)`；要求正失败数、至少1000shots。无仿真wall-time截止，仅有shots上限。首次Ctrl+C会在当前wave完整保存后停止；重新执行继续断点。

主要文件：

- `selected_candidates.json`、`historical_ranking.csv`：本次选择及历史排名。
- `distance_records.json`：每个电路的新距离上下界、精确性及任务状态。
- 距离工作目录内的分片JSON：已完成/未决任务及可复用断点。
- `points/*/result.json`：原始计数、seed、完整波次记录和断点；同目录保存电路与DEM。
- `simulation_summary.csv`、`simulation_summary.json`：新的独立统计结果。
- `two_point_slopes.csv`：新的两点斜率及名义误差；未完成或零失败点不伪造斜率。
- `independent_confirmation.png/.pdf`：仅等效per logical per round，附caption、Wilson 1σ误差条与距离标注；旧数据不混入新曲线。

开放符号表示达到shots上限但尚未达到目标精度；零失败点使用上限箭头。自适应停止下的Wilson区间是名义区间，不是anytime-valid置信序列。

依赖与算法参考：[Stim](https://github.com/quantumlib/Stim)、[PySAT](https://pysathq.github.io/docs/html/api/solvers.html)、[IBM BB memory paper](https://arxiv.org/abs/2308.07915)。这里只复用已经采用的噪声模型和电路构建规则，W7候选是原扫描中的具体schedule。

## 随包验证记录

`validation/` 保存本次本地运行记录。实际四个电路重新完成小权重下界验证；baseline精确6，第三名候选精确6，前两名未定。局部SAT测试使用每逻辑bit一个分片、每任务2秒，生产默认则为16分片、300秒/任务。不能用短预算UNKNOWN推断不存在六故障逻辑错误。

另使用真实Stim采样和Tesseract beam=2，对第一名候选试跑了640shots：p=.003有384shots，p=.0025有256shots，均零joint失败及零解码异常；验证了中断后完整波次续跑及无重复计数。样本量不足以比较性能，未达到20%精度，也没有计算两点斜率。这些验证样本不会混入你服务器的新结果。
