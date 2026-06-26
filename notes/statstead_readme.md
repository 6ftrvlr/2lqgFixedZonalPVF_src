# 基于相邻时间窗的一阶、二阶统计量偏移的统计稳态判定方法

## 1. 目的

本文档记录一种用于数值模拟后处理的统计稳态判定方法。该方法面向若干诊断量时间序列，例如

$$
X_j(t),\quad j=1,2,\dots,N_d,
$$

通过比较相邻时间窗内的一阶统计量和二阶统计量的变化，判断每个诊断量进入统计稳态的时间，并以所有诊断量判稳时间的最大值作为系统进入统计稳态的时间。

该方法的基本思想是：

> 若某个诊断量在连续多个相邻时间窗中，其均值和标准差的相对偏移均小于指定阈值，则认为该诊断量已经进入统计稳态。

---

## 2. 输入数据结构

程序读入 `cases_stat.mat`，其中包含 `CaseMap`。每个 case 由一个名称 `cn` 索引：

```matlab
dataDir = './fumo/fumo';
matFile = fullfile(dataDir, "cases_stat.mat");
```

`CaseMap(cn)` 为一个结构体，典型字段包括：

```matlab
sourceFile
sourceName
nRows
number
t
EHF
U
KE1
KE2
PE
type
constantKind
constantValue
name
baseName
tStart
tEnd
```

其中：

* `t` 为时间序列自变量；
* `EHF, U, KE1, KE2, PE` 为待检测的诊断量；
* `constantValue` 用作后续作图中的自变量；
* `type, constantKind` 等字段可用于定义不同 case 的类别标签。

---

## 3. 待判稳诊断量

在当前实现中，标注为 `o` 的五个诊断量为：

```matlab
diagFields = ["EHF", "U", "KE1", "KE2", "PE"];
```

对每个 case，程序分别判断这五个诊断量进入统计稳态的时间：

$$
T_{\mathrm{st},j}.
$$

系统整体进入统计稳态的时间定义为：

$$
T_{\mathrm{st}}
===============

\max_j T_{\mathrm{st},j}.
$$

也就是说，只有所有指定诊断量都通过判稳后，才认为系统进入统计稳态。

---

## 4. 时间窗设置

用户指定两个统一参数：

```matlab
windowLength = W;
windowStep   = \Delta W;
```

例如：

```matlab
windowLength = 0.2;
windowStep   = 0.05;
```

程序采用相邻两时间窗比较。对某个起点 (s)，定义：

$$
W_1=[s,s+W),
$$

$$
W_2=[s+W,s+2W).
$$

其中 (W_1) 为前窗，(W_2) 为后窗。

每次窗口整体向前移动 (\Delta W)。

---

## 5. 一阶与二阶统计量

对给定诊断量 (X(t))，在两个相邻窗口中分别计算均值：

$$
\mu_1 = \langle X\rangle_{W_1},
$$

$$
\mu_2 = \langle X\rangle_{W_2}.
$$

二阶统计量取标准差：

$$
\sigma_1 = \mathrm{std}(X)_{W_1},
$$

$$
\sigma_2 = \mathrm{std}(X)_{W_2}.
$$

---

## 6. 偏移判据

均值相对偏移定义为：

$$
\delta_\mu
==========

\frac{|\mu_2-\mu_1|}{S_\mu}.
$$

标准差相对偏移定义为：

$$
\delta_\sigma
=============

\frac{|\sigma_2-\sigma_1|}{S_\sigma}.
$$

一个窗口对通过判稳的条件为：

$$
\delta_\mu \le \varepsilon_\mu
\quad \text{and} \quad
\delta_\sigma \le \varepsilon_\sigma.
$$

例如：

```matlab
meanTol = 0.01;   % 1 %
sdTol   = 0.03;   % 3 %
```

---

## 7. 归一化尺度设置

当前设计保留两种归一化尺度选项。

### 7.1 使用参考时间后的整体统计量

设参考时间为：

```matlab
refTime = 3.0;
```

对 (t\ge t_{\mathrm{ref}}) 后的全部样本计算参考尺度。

均值尺度取：

$$
S_\mu = \left| \langle X\rangle_{t\ge t_{\mathrm{ref}}} \right|.
$$

标准差尺度取：

$$
S_\sigma = \mathrm{std}(X)*{t\ge t*{\mathrm{ref}}}.
$$

### 7.2 不设置参考尺度

也可以不设置参考尺度，此时：

$$
S_\mu = |\mu_2|,
$$

$$
S_\sigma = \sigma_2.
$$

这种方式要求用户保证相关诊断量的均值不会接近零。若程序检测到归一化尺度接近零，则应报错中止，避免得到没有意义的相对偏移。

---

## 8. 连续通过次数

为了避免偶然通过，用户指定连续通过次数：

```matlab
nPass = 100;
```

如果某个诊断量连续 (N_{\mathrm{pass}}) 个窗口对均满足：

$$
\delta_\mu \le \varepsilon_\mu,
$$

$$
\delta_\sigma \le \varepsilon_\sigma,
$$

则认为该诊断量已经进入统计稳态。

---

## 9. 最早允许判稳时间

用户指定最早允许判断稳定的时间：

```matlab
tEarliest = 4.0;
```

在此时间之前，即使窗口统计量满足偏移判据，也不允许程序判定系统已经进入统计稳态。

---

## 10. 稳态时间定义

对某个通过判据的窗口对：

$$
W_1=[s,s+W),
$$

$$
W_2=[s+W,s+2W),
$$

稳态时间定义为后窗 (W_2) 的前端：

$$
T_{\mathrm{st}} = s+W.
$$

若需要连续 (N_{\mathrm{pass}}) 次通过，则程序在首次达到连续通过次数时，取对应窗口对中后窗前端作为该诊断量的判稳时间。

---

## 11. 系统整体判稳流程

对一个 case，分别得到五个诊断量的判稳时间：

$$
T_{\mathrm{st,EHF}},
\quad
T_{\mathrm{st,U}},
\quad
T_{\mathrm{st,KE1}},
\quad
T_{\mathrm{st,KE2}},
\quad
T_{\mathrm{st,PE}}.
$$

系统整体稳态时间定义为：

$$
T_{\mathrm{st}}
===============

\max
\left(
T_{\mathrm{st,EHF}},
T_{\mathrm{st,U}},
T_{\mathrm{st,KE1}},
T_{\mathrm{st,KE2}},
T_{\mathrm{st,PE}}
\right).
$$

采用最大值的原因是：若某个诊断量仍未达到统计稳定，则整个系统仍不应被判定为统计稳态。

---

## 12. 稳态段统计量计算

得到系统整体稳态时间 (T_{\mathrm{st}}) 后，程序对每个诊断量在区间：

$$
t\in [T_{\mathrm{st}}, T_{\mathrm{end}}]
$$

内计算统计量。

需要记录的统计量包括：

1. 样本数；
2. 均值；
3. 标准差；
4. 方差；
5. 三阶无量纲中心矩，即偏度；
6. 四阶无量纲中心矩，即峰度；
7. 最小值；
8. 第一四分位值；
9. 中位数；
10. 第三四分位值；
11. 最大值。

---

## 13. 前四阶统计量定义

设稳态段样本为：

$$
X_i,\quad i=1,2,\dots,N.
$$

均值为：

$$
\mu = \frac{1}{N}\sum_{i=1}^{N}X_i.
$$

标准差为：

$$
\sigma =
\left[
\frac{1}{N-1}\sum_{i=1}^{N}(X_i-\mu)^2
\right]^{1/2}.
$$

方差为：

$$
\sigma^2.
$$

三阶统计量采用偏度：

$$
\mathrm{skewness}
=================

\frac{1}{N}
\sum_{i=1}^{N}
\left(
\frac{X_i-\mu}{\sigma}
\right)^3.
$$

四阶统计量采用峰度：

$$
\mathrm{kurtosis}
=================

\frac{1}{N}
\sum_{i=1}^{N}
\left(
\frac{X_i-\mu}{\sigma}
\right)^4.
$$

这里三阶和四阶统计量使用前两阶统计量进行无量纲化，因此可用于比较不同量纲、不同幅值的诊断量的分布形态。

---

## 14. 输出数据

统计后处理结果追加保存到原文件：

```matlab
cases_stat.mat
```

建议保存以下对象：

```matlab
CaseMap
CaseStatMap
SteadyStatTable
p
```

其中：

* `CaseMap`：原 case 信息，并可追加每个 case 的系统稳态时间；
* `CaseStatMap`：按 case 保存完整判稳细节和稳态段统计量；
* `SteadyStatTable`：整理后的表格，便于作图、筛选和导出；
* `p`：本次后处理所使用的参数结构体。

---

## 15. pointType 字段

为了便于后续作图时对数据点分门别类，统计结果中预留至少两个具有 `pointType` 含义的字段：

```matlab
pointType1
pointType2
```

默认可设置为：

```matlab
pointType1 = type;
pointType2 = constantKind;
```

用户可根据具体课题自行修改，例如按 forcing 类型、控制参数类型、边界条件类型、数值实验组别等重新定义。

---

## 16. 作图方式

遍历 case 名称集合 `cns`，以：

$$
| \mathrm{constantValue} |
$$

作为横坐标。

当前实现可支持至少一种图形。

### 16.1 均值-标准差图

横坐标为：

$$
|\mathrm{constantValue}|.
$$

纵坐标为稳态段均值，并用误差棒表示标准差：

$$
\mu \pm \sigma.
$$

### 16.2 箱线图

横坐标为：

$$
|\mathrm{constantValue}|.
$$

每个横坐标位置对应相应 case 在稳态段的样本分布。

箱线图用于观察中位数、四分位区间、极值范围和不同 case 之间的波动差异。

---

## 17. 注意事项

### 17.1 窗口长度不能过短

窗口长度 (W) 应足以覆盖主要振荡周期或主要相关时间。如果窗口太短，均值和标准差会受到瞬时波动影响，导致误判。

### 17.2 窗口长度不能过长

如果窗口过长，判稳时间会被推迟，并且短暂的慢漂移可能被窗口平均掩盖。

### 17.3 连续通过次数不宜过小

若 `nPass` 太小，可能把偶然平缓段误判为统计稳态。使用连续通过次数可以提高判稳的鲁棒性。

### 17.4 归一化尺度接近零时必须谨慎

当均值或标准差接近零时，相对偏移可能失去意义。因此程序应在归一化尺度小于指定容忍值时停止，并提示用户改用参考尺度或更换判据。

### 17.5 系统稳态应由多个诊断量共同决定

单一诊断量进入平台区并不意味着整个系统进入统计稳态。例如，某个平均流量可能已经稳定，但能量分量、通量或关键模态仍可能缓慢演化。

### 17.6 累计统计量不作为主判据

从 (0) 到 (t) 的累计均值和累计标准差可以用于检查统计量是否收敛，但不建议作为判定进入稳态时间的唯一依据。

原因是累计统计量会长期保留初始暂态的影响，并且随着 (t) 增大天然变得越来越平滑，可能掩盖后期慢漂移。

---

## 18. 方法总结

本方法最终采用：

> 多诊断量 + 相邻时间窗 + 均值偏移 + 标准差偏移 + 连续通过次数

作为统计稳态判据。

单个诊断量的判稳条件为：

$$
\frac{|\mu_2-\mu_1|}{S_\mu}
\le
\varepsilon_\mu,
$$

$$
\frac{|\sigma_2-\sigma_1|}{S_\sigma}
\le
\varepsilon_\sigma,
$$

并且该条件连续满足 (N_{\mathrm{pass}}) 次。

系统整体稳态时间为所有诊断量判稳时间的最大值：

$$
T_{\mathrm{st}}
===============

\max_j T_{\mathrm{st},j}.
$$

之后，所有统计量均在：

$$
[T_{\mathrm{st}},T_{\mathrm{end}}]
$$

上计算。

该方法适合用于批量数值实验的自动化后处理，能够将“看图判断是否达到统计稳态”的经验过程转化为可重复、可记录、可编码实现的定量流程。
::: 
