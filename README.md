# Weight-7 QEC simulation projects

本仓库保存以下两个项目的主要源码、运行脚本、依赖清单、说明文档和运行所需的小型输入数据，保留原有目录结构。

| 项目 | 说明与运行入口 |
| --- | --- |
| Weight-7 / BB72 / surface d=7 联合模拟 | [README](w7_bb72_surface7_joint_tesseract_380/w7_joint_tesseract/README_zh.md) |
| Weight-7 电路调度扫描 | [README](w7_schedule_scan_380/w7_schedule_scan_380/README_zh.md) |
| 调度扫描内的 high-slope 实验程序 | [README](w7_schedule_scan_380/w7_schedule_scan_380/w7_high_slope_dc7_380/w7_high_slope_dc7_380/README_zh.md) |

## 上传范围

- 包含 Python 源码、Shell 脚本、`requirements.txt`、原有文档和第三方许可证。
- 保留调度种子、基线距离及候选数据等程序输入。`data/old_scan.log` 是 high-slope 程序默认读取的候选选择输入，因此也保留。
- 不包含虚拟环境、Python 缓存、批量实验结果、生成的图表、验证输出以及原归档的 `SHA256SUMS.txt`。原始文档中提到的历史输出可在本地原始项目中查看。
- `.gitignore` 用于防止后续提交虚拟环境、缓存和实验输出。

运行前请进入对应程序目录，按照其 README 安装依赖并运行。实验输出需自行生成。本次整理保留原始程序文件内容，未执行完整实验或重新验证历史数值结论。
