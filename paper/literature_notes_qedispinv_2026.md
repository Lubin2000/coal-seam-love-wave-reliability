# QEDispInv 与三篇关联论文核验记录

## 软件基本信息

- **名称：** QEDispInv
- **地址：** https://github.com/pan3rock/QEDispInv
- **默认分支：** `main`
- **检查版本：** 提交 `ed1b5dd7`，提交时间 2026-05-31
- **主要语言：** C++17、Fortran，Python 脚本用于绘图和结果转换
- **依赖：** Eigen3、OpenMP、CMake、fmt、toml11、HighFive、LBFGSpp
- **仓库许可证：** 根目录 `LICENSE` 为 LGPL-3.0
- **许可证风险：** C++ 源文件头部仍写为 GPL-3.0，与根目录 LGPL-3.0 不一致；复制或修改代码前应要求作者澄清
- **成熟度：** 37 stars、8 forks；无 release、无 tag、未见独立测试套件；截至检查日有一个关于 GCC 15 未限定 `abs()` 导致 NaN 和漏根的开放问题

## 三篇关联论文

### 1. Pan & Chen (2025)

- **题目：** Efficient Computation of Dispersion Curves in Low-Velocity-Layered Half-Spaces
- **期刊：** *Bulletin of the Seismological Society of America*
- **卷期页：** 115(4), 1381–1391
- **DOI：** https://doi.org/10.1785/0120250003
- **核心：** 针对低速层模型中相邻模态靠近的 mode-kissing 区域，改进频散函数采样与有效根计算。
- **对本文价值：** 如果高阶根在正演阶段被漏掉，就不能把“模态不可用”解释为物理缺失。可用于设计根值交叉验证和漏根率统计。

### 2. Pan, Wang & Chen (2026)

- **题目：** Surface-Wave Dispersion Curve Computation and Inversion: A Framework Integrating Quadratic Extrema Interpolation and Randomized Layering with Multiple Initial Models
- **期刊：** *Bulletin of the Seismological Society of America*
- **卷期页：** 116(3), 902–922
- **DOI：** https://doi.org/10.1785/0120250207
- **核心：** 二次极值插值、随机分层、按波长设置层厚和多初始模型 L-BFGS-B 反演。
- **对本文价值：** 可作为独立正演和“多初始模型”不确定性的对照，但不等同于煤厚—速度—\(Q\) 联合反演。

### 3. Pan, Chen, Wang, Yang & Zhang (2019)

- **题目：** Sensitivity Analysis of Dispersion Curves of Rayleigh Waves with Fundamental and Higher Modes
- **期刊：** *Geophysical Journal International*
- **卷期页：** 216(2), 1276–1303
- **DOI：** https://doi.org/10.1093/gji/ggy479
- **核心：** 比较基阶与高阶 Rayleigh 频散对层参数的敏感性。
- **对本文价值：** 支持高阶信息价值，同时提醒多参数权衡；不能把高阶数量增加直接等同于厚度、速度和 \(Q\) 都可独立识别。

## 仓库实际功能核验

### 支持

- 任意一维分层模型的 Rayleigh/Love 相速度计算；Love 波使用 `--sh`。
- 计算 0 到指定最高阶模态，适合检查煤层次级模态和 mode-kissing。
- 输出格式为 `频率(Hz) 相速度(km/s) 模态号`。
- 可计算 Vp、Vs、密度灵敏度核并保存为 HDF5。
- 反演：随机层厚/参考模型、多初始模型、L-BFGS-B、一阶或自适应平滑正则化。
- 反演数据可以指定不同模态和模态权重。

### 不支持或不适合直接替代本文

- 不含衰减、复频散或 \(Q\) 参数；Fortran 内核明确注明 `no attenuation`。
- 反演变量主要是分层 \(V_S\)；Vp 和密度由固定参考模型或经验关系生成。
- 煤厚不是显式参数，层界面是随机层厚和正则化平滑对象。
- 频散类型主要是相速度。底层程序含群速度/灵敏度计算，但当前 CLI 的主要输出不是群速度曲线。
- `P10/P50/P90` 是多初始反演集合按拟合度加权的分布统计，不等同于经过覆盖率校准的贝叶斯可信区间。
- Python 文件是绘图和后处理脚本，不是完整 Python API 绑定。
- 未见自动化测试；GCC 15 已报告数值求根问题，使用前应固定编译器/提交版本并做基准检查。

## 对第一篇论文的正确用法

### 1. 相速度独立基准

构造与本文相同的三层非对称煤层模型：

- 顶板 \(\beta_0=2400\) m/s；
- 煤层 \(h=2.4\) m、\(\beta_1=1120\) m/s；
- 底板 \(\beta_2=2650\) m/s；
- 频率 60–920 Hz；
- 计算模态 0、1、2 的 Love 相速度。

比较：

- 本文解析方程的根；
- QEDispInv 的根；
- 合成能量图峰位。

统计相速度相对误差、模态存在/缺失一致性以及 mode-kissing 频段。

### 2. 漏模态诊断

对随机煤层模型，区分三类“模态缺失”：

1. 物理截止：解析方程在可追踪速度范围内无实根；
2. 数值漏根：解析方程存在根，但 QEDispInv 的常规扫描漏掉；
3. 可观测性缺失：根存在，但合成能量图或现场记录中无可拾取能量。

这一分类可以直接提升本文模态可用率研究的可信度。

### 3. 不作为联合反演的直接替代

QEDispInv 不包含衰减，不能完成“频散—衰减联合反演”；其 \(V_S\) 剖面反演也不能自动恢复三层模型的具体煤厚和物理界面。若把其 P10–P90 当作本文参数后验，会混淆“多初始模型集成范围”和“校准可信区间”。

## 建议的交叉验证表

| 项目 | 本文解析正演 | QEDispInv | 合成能量图 |
|---|---|---|---|
| 模态 0–2 是否存在 | 解析根 | 独立求根 | 峰值与跟踪结果 |
| 相速度根值 | 方程解 | 交叉基准 | 从能量图拾取 |
| 群速度 | 频散差分 | 当前 CLI 不宜直接使用 | 可作为独立近似 |
| \(Q\)/衰减 | 本文一阶模型 | 不支持 | 仅本文流程 |
| 界面/煤厚 | 显式参数 | 难与煤厚一一对应 | 反演结果 |
| 不确定度 | Fisher/局部近似 | 多初始集统计 | 本文稳健反演/保形校准 |

## 复现注意事项

1. 固定提交 `ed1b5dd7`，不要直接跟随 `main`。
2. 在编译前检查 `std::abs` 修正和 GCC 版本；避免使用已报告的 GCC 15 组合而不做回归测试。
3. 模型文件单位为 km、密度 g/cm³、速度 km/s；与本文的 m、m/s、kg/m³ 不同，需统一单位。
4. 先复现仓库自带 `demo/lvl-l4`，再运行三层煤层模型。
5. 保存程序版本、编译器和编译选项，避免结果不可追溯。
6. 如果只调用独立可执行文件，按 LGPL 分发边界处理；如果复制源码进本项目，先澄清 LGPL/GPL 头文件冲突。

## Windows 编译与实际验证结果

2026-11-06 已在便携式 MSYS2 UCRT64 环境中成功编译提交 `ed1b5dd7` 的 `forward` 模块：

- GCC 16.2.0；
- gfortran 16.2.0；
- CMake 4.4.5；
- Eigen 5.0.1；
- HDF5 2.2.0；
- OpenMP 关闭。

为适配当前依赖，只在构建副本中做了以下平台修补：

1. 移除上游过严的 Eigen 3.3 版本限制；
2. 为 Eigen 5 显式传递 include 目录；
3. 在 HighFive 中取消 Windows 宏 `FILE_CREATE` 的冲突。

没有修改频散方程或 Fortran 数值内核。

煤层验证模型采用 2.4 m 煤层、1.12 km/s 煤体横波速度，100 m 顶板近似上半空间。QEDispInv 同时产生顶板自由表面模态，因此对每个解析煤层模态匹配全部 QEDispInv 输出中最近的根。

100–560 Hz 结果：

- 34 个煤层模态点全部匹配；
- 27/34 误差小于 1 m/s；
- 34/34 误差小于 10 m/s；
- 中位误差 0.477 m/s；
- 平均误差 0.972 m/s；
- P95 误差 2.582 m/s；
- 最大误差 9.104 m/s。

最大误差位于 100 Hz，原因是解析基阶接近顶板速度，对“顶板半空间替换为有限厚层”的边界近似最敏感。
