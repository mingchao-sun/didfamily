# didfamily:在 jamovi 中比较各种双重差分(DiD)估计量

用同一份面板数据一次设定,并排比较多种现代 DiD 估计量的结果,并提供诊断、敏感性分析和可复现的 R 代码。版本 0.5.0(开发中)。

## 十一个分析(jamovi 菜单:DiD,按工作流程排列并分组)

1. **面板结构检查**:估计前检查面板(平衡性、重复行、缺失值、处理是否吸收、从未处理组占比、各队列规模),并画出处理状态图。问题只列出,不会中断。
2. **经典DiD(2×2)**:两组两期,可带或不带个体 ID;给出分组均值、回归 DiD,以及双重稳健 / IPW / 结果回归估计量(`DRDID`),图中画出隐含的反事实路径。
3. **估计量**(子菜单):每个估计量一个独立分析,各有自己的选项,输出总体 ATT、动态效应表、事件研究图,并可生成 R 代码:**双向固定效应(TWFE)**、**Callaway & Sant'Anna**(对照组、双重稳健/IPW/回归、乘数自助法、组-时热图、分队列与分日历时间效应)、**Sun & Abraham**、**Borusyak, Jaravel & Spiess(插补法)**、**Gardner 两阶段(did2s)**、**堆叠 DiD**。它们与"比较"使用同一份代码,数字完全一致。
4. **DiD估计量比较**:TWFE、Callaway & Sant'Anna、Sun & Abraham、Borusyak-Jaravel-Spiess(插补法)、Gardner 两阶段、堆叠 DiD 并排比较。
   输出总体 ATT 表、动态效应表、叠加事件研究图、森林图、Callaway & Sant'Anna 的组-时效应热图与分队列/分日历时间效应,
   支持协变量、聚类变量、Bootstrap,并可生成复现全部数字的 R 脚本。
5. **Bacon、事前趋势与 HonestDiD**:Goodman-Bacon 分解(需平衡面板)、TWFE 权重、事前趋势联合检验、HonestDiD(较慢,默认关闭,结果逐行显示)。
6. **连续处理DiD**:剂量-反应曲线 ATT(d)、边际效应 ACRT(d) 及总体平均(`contdid`)。
7. **三重差分(DDD)**:三重差分估计,并与"仅符合条件者"和"仅不符合条件者"的普通 DiD 对照。

处理变量可以是 0/1 指示变量(须为吸收型)或首次处理时期。标准误默认按个体聚类,也可指定聚类变量。可选的**抽样权重**变量(须在个体内恒定,如基期人口)用于 TWFE、Callaway & Sant'Anna、Sun & Abraham、Gardner 和堆叠 DiD;使用权重时 BJS 会被跳过(`didimputation` 不支持权重)。

**推断**:TWFE、Sun & Abraham 和堆叠 DiD 与 `fixest`、Stata `reghdfe` 一样使用 t(G-1) 分布(12 个聚类时,正态近似给出 p = 0.0009,t 分布给出 0.0069);其余估计量沿用各自包的正态推断。总表给出各估计量的观测数和聚类数,聚类很少或队列很小时面板概况会提示。
面板不平衡时,Callaway & Sant'Anna、DDD 和连续处理估计量只使用每期都有观测的个体(`did` 的不平衡面板代码路径会使 R 崩溃,故不使用);固定效应类估计量使用全部行。

## 在 R 中使用
每个分析也是一个 R 函数(Rj 编辑器、脚本均可),例如 `didfamily::didcompare(data = d, dep = "y", id = "id", time = "year", treat = "D")`;"显示可复现的 R 代码"可生成只依赖 `fixest`、`did`、`didimputation`、`did2s` 的独立脚本。各分析背后的论文和 R 包列在 jamovi 输出末尾的参考文献中。

## 数据格式
长面板:每个个体每期一行。必需:结果变量、个体 ID、整数时期、处理指示变量(0/1)或首次处理时期(0 或缺失 = 从未处理)。
第一期就已处理的个体会被剔除(没有处理前时期);含缺失值的行会被剔除并提示。

示例数据以预设了变量类型和说明的 `.omv` 文件随模块安装(ID 设为 ID 类型,0/1 指示变量设为名义,时期与队列设为连续),在 jamovi 的 **打开 > 数据库 > Difference-in-Differences Family** 中可以找到,每条说明里写明用哪个分析、变量怎么放。对应的 CSV 在 `dev/csv/`,`.omv` 可用 `dev/make_omv.R` 重新生成:`sim_staggered.csv`(交错处理,TWFE 明显有偏)、`sim_2x2.csv`(两期,趋势依赖协变量)、`sim_ddd.csv`(三重差分)、
`sim_dose.csv`(连续处理)、`mpdta.csv`(`did` 包的最低工资示例,处理变量类型选"首次处理时期",使用 `first.treat`)、`castle.csv`(Castle Doctrine 数据,交错处理的经典例子)。

**图的样式**:各分析沿用 jamovi"绘图"标签页的主题;模块另有"Plot style"选项:彩色(色盲友好)或灰度(灰度加不同点形状和线型,便于黑白印刷),事件研究图还可选误差线或阴影带及估计量间距。

## 验证
`dev/test_*.R` 中约 200 项检查:已知真实效应的模拟设计、精确恒等式(Bacon 权重、TWFE 权重、Sun & Abraham 与 `fixest::sunab`、生成的 R 代码与模块输出)、
`mpdta` 示例,以及 Castle Doctrine 数据(TWFE 0.0818,91% 的权重在"处理组 vs 从未处理组"的比较上)。

## 已知限制
- 抽样权重仅在"DiD估计量比较"中提供;只支持面板数据(个体重复截面须先汇总成"单位×时期"单元);没有野 Bootstrap 聚类、没有一致置信带。
- 没有 dCDH(`DIDmultiplegtDYN` 依赖 jamovi 软件包快照中没有的 `polars`)、没有合成 DiD(`synthdid` 仅在 GitHub)、没有 `fect`(其 FEct 估计量与已有的插补法相同)。
- Wooldridge 扩展 TWFE 未单列:对从未处理组,它与插补法的处理后效应相同;加入前置项后等同于 Sun & Abraham。
- 没有预期效应期数,也不能自选参考期(事件时间相对于 e = -1)。
- 界面标签、表格与选项已提供日文和简体中文翻译,说明性文字和提示信息仅有英文。

更多细节见英文版 [README.md](README.md) 与 `dev/NOTES.md`。
