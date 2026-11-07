# 文献核验记录：Schott and Waclawik (2015)

## 基本信息

- **作者：** W. Schott, P. Waclawik
- **题目：** On the quantitative determination of coal seam thickness by means of in-seam seismic surveys
- **期刊：** Canadian Geotechnical Journal
- **卷页：** 52, 1496–1504
- **年份：** 2015
- **DOI：** https://doi.org/10.1139/cgj-2014-0466
- **文件：** 用户提供的 9 页 PDF

## 研究背景与目标

研究地点为捷克 Karvina 煤矿的两个相邻盘区。煤层受侵蚀影响，厚度可从约 0.3 m 变化到约 4.4 m。研究目标是利用槽波透射数据圈定侵蚀区，并定量给出煤层厚度分布，为采掘设计提供依据。

## 数据与方法

- 使用井下炸药震源和双分量检波器记录 SH/Love 型槽波。
- 每个探头包含相互垂直的 X、Y 分量，能够记录 Love 槽波的完整振动。
- 震源和检波器间距约 20 m。
- 使用多重滤波分析确定瞬时群速度或常频率群走时。
- 选择 200 Hz 作为目标频率，对群走时执行 SIRT 层析，获得群速度分布。
- 利用巷道已知厚度与层析网格群速度建立三阶多项式经验关系。
- 将群速度层析结果转换为厚度图。
- 结果用于规划三个盘区的采掘，并经过后续采掘揭露验证。

## 主要结果

- 200 Hz 的群速度—厚度关系主要能够分辨约 1–3 m 厚度范围。
- 群速度层析成功圈定多个侵蚀区及其向盘区内部的延伸范围。
- 采掘时在预测边界附近遇到约 2.3–2.6 m 的煤层厚度，与预测目标基本一致。
- 另一区域实际采掘遇到 0.9–2.4 m 的较大厚度变化。
- 作者指出：如果要覆盖更大厚度范围或完整厚度分布，应在多个常数频率重复处理并拼接结果。

## 对本文最有价值的发现

### 1. 群速度具有现场可操作性

该研究使用群走时而没有直接使用相速度，因为现场相速度更难稳定提取。本文因此将群速度作为独立信息源，并证明它与高阶模态联合后能够明显缩短厚度区间。

### 2. 单一频率只覆盖有限的厚度范围

200 Hz 只适用于有限厚度区间；完整厚度范围需要多频率。本文的信息预算从理论上解释了这种频率—厚度敏感性，并给出不同观测组合的理论区间宽度。

### 3. 现场经验关系与理论关系存在显著差异

作者将巷道已知厚度与层析群速度拟合为三阶多项式，但该经验曲线与理论群速度—厚度曲线差异较大。作者认为可能原因包括：

- 巷道局部厚度与层析网格平均速度不匹配；
- 层析迭代中的平滑处理；
- 更复杂的、尚未识别的模型误差和传播效应。

这与本文的设计原则高度一致：不能把确定性正演误差忽略，也不能把窄后验区间直接解释为可信度。

### 4. 现场验证为合成研究提供了外部约束

本文的基线速度、密度和频率范围参考了该研究报道的量级：

- 顶板横波速度约 2450 m/s；
- 煤体横波速度约 1000 m/s；
- 顶底板密度约 2.5 g/cm³；
- 煤密度约 1.3 g/cm³；
- 常用槽波工作频率约 100–600 Hz。

## 与该研究的明确区别

本文不是对该现场层析方法的替代，而是补充：

- 从单一基阶和单一常数频率扩展到多模态、群速度和衰减；
- 从经验厚度标定扩展到基于精确频散的可辨识性信息预算；
- 从人工识别频散扩展到自动模态跟踪和伪拾取分析；
- 从确定性结果扩展到稳健参数反演和不确定度评价；
- 从现场空间层析扩展到受控合成机理实验。

## 适合写入论文的位置

1. **引言：** 说明 ISS 已能通过群速度层析进行现场厚度制图，但单频、基阶和经验标定存在限制。
2. **方法参数：** 用该研究的煤、岩速度和密度量级约束合成模型。
3. **讨论：** 对比现场经验关系与理论关系的差异，支持模型误差和稳健后验研究。
4. **局限性：** 明确本文尚无 Karvina 式采掘验证，不能把合成精度外推为现场精度。

## 可引用观点

可以在论文中写为：

> Schott and Waclawik (2015) demonstrated field-scale quantitative seam-thickness mapping using group-velocity tomography at a fixed frequency and an empirical velocity–thickness calibration. Their study also reported a substantial discrepancy between the empirically fitted and theoretically predicted group-velocity–thickness relations, highlighting unresolved model and upscaling errors.

不要将该研究描述为“使用高阶模态”或“进行了概率不确定度分析”，因为原文并未这样处理。

## 仍需补读的相关文献

该文参考文献中可优先查找：

- Räder et al. (1985): Love channel wave dispersion and amplitude–depth distributions.
- Dresen and Rüter (1994): In-seam seismic handbook.
- Schott and Brandt (2003): seam thickness and seam splitting from ISS.
- Buchanan and Jackson (1986): Coal geophysics.
- Dobróka (2005): tomographic imaging of 3D waveguide structures.
- Dziewonski et al. (1969): multiple-filter analysis.

