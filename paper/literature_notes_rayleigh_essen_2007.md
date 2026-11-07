# 文献核验记录：Essen et al. (2007)

## 正式题录

- **作者：** Katja Essen, Thomas Bohlen, Wolfgang Friederich, Thomas Meier
- **标题：** Modelling of Rayleigh-type seam waves in disturbed coal seams and around a coal mine roadway
- **期刊：** Geophysical Journal International
- **卷期：** 170, 511–526
- **年份：** 2007
- **DOI：** https://doi.org/10.1111/j.1365-246X.2007.03436.x
- 用户另外提供了 16 页正式出版版；此前三页作者稿已由正式版替代。

## 方法与模型

- 简单横向均匀模型使用 Green 函数方法；
- 复杂二维/三维几何使用有限差分波场模拟；
- 研究煤层冲刷、尖灭、分叉、变薄、断层和巷道；
- 波场包含 P-SV 体波和 Rayleigh 型槽波。

## 主要结论

1. 简单三层煤层中，煤层内波场以基阶 Rayleigh 槽波为主；基阶模的垂直分量关于煤层中心对称，水平分量反对称。
2. 煤层中含高速、高密度夹矸时，高阶模式可能主导传播，具体模态取决于夹矸厚度。
3. 冲刷、尖灭或煤层完全中断会产生模式转换：一部分成为反射槽波，大量能量转化为 S 波并向顶板围岩辐射。
4. 煤层分叉时反射槽波很弱，可能被噪声掩盖；仅依靠反射槽波不能可靠探测分叉。
5. 扰动反射振幅依赖断层落差、变薄或冲刷程度；高频相位更容易反射，低频相位更容易透射。
6. 巷道端面会把槽波转换为沿巷道表面传播的 Rayleigh 波，其频率略低于煤层槽波。
7. 煤层槽波在垂直—径向平面内呈椭圆偏振，巷道 Rayleigh 波在径向—横向平面内呈椭圆偏振，因此可通过三分量偏振分析区分。

## 补充方法参数

- Green 函数方法：横向均匀模型，频率域计算完整地震图，包含泄漏模、体波和衰减；
- 有限差分方法：并行粘弹性速度—应力方程，空间四阶或二阶、时间二阶交错网格；
- 二维/三维扰动模型源为 500 Hz Ricker 子波，垂直单力；
- 巷道采用真空形式近似自由表面；
- 模型涵盖 1–2 m 夹矸、冲刷、尖灭、变薄、断层、分叉和巷道端面。

## 对第一篇论文的支持

- 解释高阶模态并不只来自厚度连续变化，也可由夹矸和结构扰动产生；
- 说明能量图中的额外波包可能代表模式转换，而不一定是基阶伪峰；
- 支持多候选模态和物理连续性约束；
- 说明仅使用反射信息会漏掉某些工况；
- 与姬广忠等（2012）的巷道波结论相互印证。

## 边界

- 该文重点研究 Rayleigh 型槽波，而第一篇论文采用 Love 型 SH 槽波；
- 其目标是波场与散射建模，不是自动多模态参数反演；
- 当前第一篇不模拟巷道、围岩 S 波辐射、煤层分叉或 Rayleigh 巷道波。

## 论文中的谨慎表述

可以写：

> Essen et al. (2007) showed that seam disturbances can cause strong mode conversion and that reflected seam waves may be weak or absent for seam splitting, implying that higher-mode and transmitted-wave information should be interpreted jointly rather than through reflection alone.

不要写：

> Essen et al. already solved automatic higher-mode inversion for thickness.

该文没有完成这一工作。
