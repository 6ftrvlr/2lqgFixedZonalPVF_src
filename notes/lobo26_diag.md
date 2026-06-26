# 2LQG 系统能量传递诊断方案：基于 EAPE、BC EKE 与 BT EKE 分解

## 1. 目标

本文档给出两层准地转系统（2LQG）的能量诊断方案。目标是从模型输出的两层扰动流函数 (\psi_1,\psi_2) 出发，计算并统计以下三类能量及其谱空间预算项：

1. EAPE：eddy available potential energy，涡动有效位能；
2. BC EKE：baroclinic eddy kinetic energy，斜压涡动动能；
3. BT EKE：barotropic eddy kinetic energy，正压涡动动能。

诊断重点包括：

* EAPE、BC EKE、BT EKE 的二维谱和各向同性谱；
* EAPE 与 BC EKE 之间由界面垂直速度控制的转换项 (T^W)；
* BC/BT 模态之间的线性与非线性能量转换；
* 各能量类型内部的非线性尺度再分配项；
* 底摩擦对 BT 与 BC 模态的作用；
* 统计稳态窗口内的时间平均；
* 可选的 (T^W) 机制分解与相干性分析。

---

## 2. 基本变量与模态分解

两层扰动流函数为

[
\psi_1,\qquad \psi_2.
]

定义正压模态和斜压模态：

[
\psi=\frac{\psi_1+\psi_2}{2},
\qquad
\tau=\frac{\psi_1-\psi_2}{2}.
]

因此

[
\psi_1=\psi+\tau,
\qquad
\psi_2=\psi-\tau.
]

相对涡度为

[
\zeta=\nabla^2\psi,
\qquad
\zeta_\tau=\nabla^2\tau.
]

界面扰动高度与斜压流函数的关系为

[
\eta=\frac{f_0(\psi_2-\psi_1)}{g_r}
=-\frac{2f_0\tau}{g_r}.
]

两层等厚时，变形半径满足

[
\lambda^{-2}=\frac{4f_0^2}{Hg_r}.
]

其中：

* (f_0)：参考 Coriolis 参数；
* (g_r)：约化重力；
* (H)：总水深；
* (\lambda)：两层 QG 变形半径；
* (\kappa)：线性底摩擦系数；
* (S)：背景界面坡度；
* (U)：背景垂向剪切对应的背景流参数；
* (J(a,b)=a_x b_y-a_y b_x)：Jacobian。

---

## 3. 首先诊断界面垂直速度 (w_{3/2})

### 3.1 为什么需要 (w_{3/2})

EAPE 与 BC EKE 之间的转换项是

[
T^W=-\frac{2f_0}{H}\tau^\dagger_{\boldsymbol{k}}w_{3/2,\boldsymbol{k}}.
]

在 EAPE 预算中，(T^W>0) 表示 BC EKE 转换为 EAPE；在 BC EKE 预算中，对应项为 (-T^W)。

因此，(w_{3/2}) 是整个能量循环诊断的核心变量之一。

### 3.2 (w_{3/2}) 不是预报变量，而是诊断变量

在 QG 中，年龄散度项为

[
A_k=\partial_xu_{ag,k}+\partial_yv_{ag,k}.
]

层积分连续方程给出

[
w_{3/2}=A_1H/2,
]

[
w_b-w_{3/2}=A_2H/2.
]

底 Ekman 垂直速度为

[
w_b=\frac{\kappa H\nabla^2\psi_2}{2f_0}
=\frac{\kappa H\nabla^2(\psi-\tau)}{2f_0}.
]

在诊断程序中不需要显式预报 (u_{ag},v_{ag})。应通过 omega 方程从瞬时 (\psi,\tau) 诊断 (w_{3/2})。

### 3.3 本文使用的 (w_{3/2}) omega 方程

若采用上层背景流 (U_1=U)、下层 (U_2=0) 的记号，界面坡度为

[
S=\frac{f_0U}{g_r}.
]

此时 (w_{3/2}) 满足

[
(\nabla^2-\lambda^{-2})w_{3/2}
==============================

-\frac{f_0}{g_r}
\left[
4J(\tau,\nabla^2\psi)
+2J(\tau,f)
+J(Uy,\nabla^2\psi)
+\frac{2f_0w_b}{H}
\right].
]

如果 (f=f_0+\beta y)，则

[
J(\tau,f)=\beta \partial_x\tau.
]

并且

[
J(Uy,\nabla^2\psi)=-U\partial_x\nabla^2\psi
]

在本文的 Jacobian 定义下成立。

### 3.4 谱空间求解

在双周期区域中，对上式做二维 Fourier 变换。记

[
K^2=k_x^2+k_y^2.
]

因为

[
\widehat{\nabla^2 w}=-K^2\widehat{w},
]

所以

[
\widehat{w}_{3/2}(\boldsymbol{k})
=================================

\frac{\widehat{\mathrm{RHS}}(\boldsymbol{k})}
{-(K^2+\lambda^{-2})}.
]

数值实现建议：

1. 用伪谱法计算各 Jacobian；
2. 对 RHS 做 FFT；
3. 用 (-(K^2+\lambda^{-2})) 除 RHS；
4. 对 (\boldsymbol{k}=0) 模式，理论上 RHS 平均值应为 0；实际编码中可将 (\widehat{w}_{3/2}(0,0)=0) 作为质量守恒约束；
5. 对非线性项使用 (2/3) dealiasing 或与主模式一致的滤波策略。

---

## 4. 三类能量定义

### 4.1 BT EKE

[
\widehat{\mathrm{EKE}}_{BT}
===========================

\frac{|\boldsymbol{k}|^2|\widehat{\psi}_{\boldsymbol{k}}|^2}{2}.
]

### 4.2 BC EKE

[
\widehat{\mathrm{EKE}}_{BC}
===========================

\frac{|\boldsymbol{k}|^2|\widehat{\tau}_{\boldsymbol{k}}|^2}{2}.
]

### 4.3 EAPE

[
\widehat{\mathrm{EAPE}}
=======================

\frac{\lambda^{-2}|\widehat{\tau}_{\boldsymbol{k}}|^2}{2}.
]

---

## 5. EAPE 谱预算

EAPE 预算为

[
\partial_t\widehat{\mathrm{EAPE}}
=================================

\widehat{P}
+\widehat{R}^{PE}_{BC}
+\widehat{T}^W
+\widehat{\mathrm{ssd}}
+\mathrm{c.c.}
]

其中

### 5.1 背景位能释放 / EAPE production

[
\widehat{P}
===========

\frac{2f_0S}{H}
\widehat{\tau}*{\boldsymbol{k}}^\dagger
\widehat{\partial_x\psi}*{\boldsymbol{k}}.
]

该项表示正压流对斜压模态的搅拌，从背景斜压位能中产生 EAPE。

### 5.2 EAPE 与 BC EKE 转换

[
\widehat{T}^W
=============

-\frac{2f_0}{H}
\widehat{\tau}*{\boldsymbol{k}}^\dagger
\widehat{w}*{3/2,\boldsymbol{k}}.
]

符号解释：

* (T^W<0)：EAPE 转换为 BC EKE，常对应斜压不稳定；
* (T^W>0)：BC EKE 转换回 EAPE，常出现在大尺度，特别是强底摩擦或有 (\beta) 和喷流时。

### 5.3 EAPE 非线性再分配

[
\widehat{R}^{PE}_{BC}
=====================

-\lambda^{-2}
\widehat{\tau}*{\boldsymbol{k}}^\dagger
\widehat{J(\psi,\tau)}*{\boldsymbol{k}}.
]

该项只在尺度之间重分配 EAPE；对所有波数积分后应近似为 0。

---

## 6. BC EKE 谱预算

BC EKE 预算为

[
\partial_t\widehat{\mathrm{EKE}}_{BC}
=====================================

\widehat{T}^{L}*{BT\rightarrow BC}
-\widehat{T}^W
+\widehat{T}^{N}*{BT\rightarrow BC}
+\widehat{R}^{KE}*{BC}
+\widehat{D}*{BC}
+\widehat{\mathrm{ssd}}
+\mathrm{c.c.}
]

### 6.1 线性 BT 到 BC 转换

[
\widehat{T}^{L}_{BT\rightarrow BC}
==================================

\frac{U}{2}
\widehat{\tau}*{\boldsymbol{k}}^\dagger
\widehat{\partial_x\nabla^2\psi}*{\boldsymbol{k}}.
]

### 6.2 EAPE 到 BC EKE 转换

BC EKE 预算中对应项为

[
-\widehat{T}^W.
]

若 (T^W<0)，则 (-T^W>0)，表示 EAPE 正在转化为 BC EKE。

### 6.3 非线性 BT 到 BC 转换

[
\widehat{T}^{N}_{BT\rightarrow BC}
==================================

\widehat{\tau}*{\boldsymbol{k}}^\dagger
\widehat{J(\psi,\nabla^2\tau)}*{\boldsymbol{k}}.
]

### 6.4 BC EKE 非线性再分配

[
\widehat{R}^{KE}_{BC}
=====================

\widehat{\tau}*{\boldsymbol{k}}^\dagger
\widehat{J(\tau,\nabla^2\psi)}*{\boldsymbol{k}}.
]

该项描述 BC EKE 在尺度之间的非线性再分配。对所有波数积分后应近似为 0。

### 6.5 BC 摩擦项

[
\widehat{D}_{BC}
================

\frac{\kappa}{2}
\widehat{\tau}*{\boldsymbol{k}}^\dagger
\widehat{\nabla^2(\tau-\psi)}*{\boldsymbol{k}}.
]

该项来自底层线性摩擦在 BC 模态上的投影。

---

## 7. BT EKE 谱预算

BT EKE 预算为

[
\partial_t\widehat{\mathrm{EKE}}_{BT}
=====================================

\widehat{T}^{L}*{BC\rightarrow BT}
+\widehat{T}^{N}*{BC\rightarrow BT}
+\widehat{R}*{BT}
+\widehat{D}*{BT}
+\widehat{\mathrm{ssd}}
+\mathrm{c.c.}
]

### 7.1 线性 BC 到 BT 转换

[
\widehat{T}^{L}_{BC\rightarrow BT}
==================================

\frac{U}{2}
\widehat{\psi}*{\boldsymbol{k}}^\dagger
\widehat{\partial_x\nabla^2\tau}*{\boldsymbol{k}}.
]

### 7.2 非线性 BC 到 BT 转换

[
\widehat{T}^{N}_{BC\rightarrow BT}
==================================

\widehat{\psi}*{\boldsymbol{k}}^\dagger
\widehat{J(\tau,\nabla^2\tau)}*{\boldsymbol{k}}.
]

### 7.3 BT EKE 非线性再分配

[
\widehat{R}_{BT}
================

\widehat{\psi}*{\boldsymbol{k}}^\dagger
\widehat{J(\psi,\nabla^2\psi)}*{\boldsymbol{k}}.
]

该项只重分配 BT EKE，对所有波数积分后应近似为 0。

### 7.4 BT 摩擦项

[
\widehat{D}_{BT}
================

\frac{\kappa}{2}
\widehat{\psi}*{\boldsymbol{k}}^\dagger
\widehat{\nabla^2(\psi-\tau)}*{\boldsymbol{k}}.
]

---

## 8. 复共轭与实值化

表中每个谱预算项都需要加上 complex conjugate，即

[
\widehat{X}+\mathrm{c.c.}=2\operatorname{Re}(\widehat{X}).
]

编码时建议统一使用

```python
term_real = 2 * np.real(np.conj(field_hat) * tendency_hat * coefficient)
```

具体归一化取决于 FFT 约定。最重要的是所有能量、转换项和 tendency 使用一致的 Parseval 归一化。

---

## 9. 推荐的诊断流程

### Step 1：读取模型输出

输入变量：

* (\psi_1(x,y,t))
* (\psi_2(x,y,t))

以及参数：

* (f_0)
* (\beta)
* (g_r)
* (H)
* (\lambda)
* (U)
* (S)
* (\kappa)
* 网格长度 (L_x,L_y)
* 网格数 (N_x,N_y)

### Step 2：模态分解

```python
psi = 0.5 * (psi1 + psi2)
tau = 0.5 * (psi1 - psi2)
psi2 = psi - tau
```

### Step 3：谱导数

构造

```python
kx, ky = spectral_wavenumbers(...)
K2 = kx**2 + ky**2
lap = -K2
```

并使用

```python
d_dx_hat = 1j * kx * field_hat
d_dy_hat = 1j * ky * field_hat
lap_hat = -K2 * field_hat
```

### Step 4：计算 Jacobian

推荐伪谱法：

```python
def jacobian(a, b):
    ax = dx(a)
    ay = dy(a)
    bx = dx(b)
    by = dy(b)
    return ax * by - ay * bx
```

其中导数用谱法，乘法在物理空间完成，非线性项需要 dealiasing。

### Step 5：诊断底部垂直速度

```python
lap_psi2 = laplacian(psi - tau)
w_b = kappa * H * lap_psi2 / (2 * f0)
```

### Step 6：诊断 (w_{3/2})

计算 RHS：

```python
zeta_bt = laplacian(psi)

rhs = -(f0 / gr) * (
    4 * jacobian(tau, zeta_bt)
    + 2 * jacobian(tau, f0 + beta * y)
    + jacobian(U * y, zeta_bt)
    + 2 * f0 * w_b / H
)
```

谱空间求解：

```python
rhs_hat = fft2(rhs)
w_hat = rhs_hat / (-(K2 + lambda_inv2))
w_hat[0, 0] = 0.0
w = ifft2(w_hat).real
```

其中

```python
lambda_inv2 = 4 * f0**2 / (H * gr)
```

### Step 7：计算所有谱预算项

需要计算：

EAPE：

* (P)
* (T^W)
* (R^{PE}_{BC})

BC EKE：

* (-T^W)
* (T^L_{BT\rightarrow BC})
* (T^N_{BT\rightarrow BC})
* (R^{KE}_{BC})
* (D_{BC})

BT EKE：

* (T^L_{BC\rightarrow BT})
* (T^N_{BC\rightarrow BT})
* (R_{BT})
* (D_{BT})

每一项都先得到 Fourier 空间表达，再取

[
2\operatorname{Re}(\cdot).
]

### Step 8：做二维谱图与各向同性壳平均

保留两类结果：

1. 二维谱：((k_x,k_y)) 平面；
2. 各向同性谱：按 (|\boldsymbol{k}|) 做 shell integral。

注意：如果系统含 (\beta) 并出现 zonal jets，各向同性谱可能掩盖重要的各向异性结构。因此建议同时输出二维谱，尤其关注 (k_x=0) 处的峰值。

### Step 9：统计稳态时间平均

推荐流程：

1. 先运行模型到统计稳态；
2. 只在稳态窗口内计算能量诊断；
3. 每个 snapshot 计算一次完整预算；
4. 对预算项、能谱、二维谱进行时间平均；
5. 同时保存 budget residual，检查闭合误差。

如果有 (N_t) 个诊断时刻，则平均为

[
\overline{X}(\boldsymbol{k})
============================

\frac{1}{N_t}
\sum_{n=1}^{N_t}X_n(\boldsymbol{k}).
]

---

## 10. 预算闭合检查

需要检查以下关系：

### 10.1 非线性再分配项积分为零

[
\sum_{\boldsymbol{k}}R_{BT}\approx 0,
]

[
\sum_{\boldsymbol{k}}R^{KE}_{BC}\approx 0,
]

[
\sum_{\boldsymbol{k}}R^{PE}_{BC}\approx 0.
]

### 10.2 (T^W) 在 EAPE 与 BC EKE 之间抵消

EAPE 预算中是 (+T^W)，BC EKE 预算中是 (-T^W)。两者相加后应抵消。

### 10.3 线性 BT/BC 转换逐波数抵消

[
T^L_{BT\rightarrow BC}
+
T^L_{BC\rightarrow BT}
\approx 0
]

应在每个波数上近似成立，前提是符号和 FFT 归一化一致。

### 10.4 非线性 BT/BC 转换全局积分抵消

[
\sum_{\boldsymbol{k}}
T^N_{BT\rightarrow BC}
+
\sum_{\boldsymbol{k}}
T^N_{BC\rightarrow BT}
\approx 0.
]

---

## 11. 可选：分解 (T^W) 的物理来源

因为 (w_{3/2}) 的诊断方程是线性的，可以把 (w_{3/2}) 分解为不同物理过程的贡献：

[
w_{3/2}
=======

w[J(\tau,\nabla^2\psi)]
+
w[J(\tau,f)]
+
w[J(Uy,\nabla^2\psi)]
+
w[w_b].
]

对应地，

[
T^W
===

T^W[J(\tau,\nabla^2\psi)]
+
T^W[J(\tau,f)]
+
T^W[J(Uy,\nabla^2\psi)]
+
T^W[w_b].
]

每一项的计算方式：

1. 只保留 RHS 中对应的一个 forcing；
2. 单独求解 Helmholtz 方程得到对应的 (w)；
3. 用

[
T^W[\cdot]
==========

-\frac{2f_0}{H}
\widehat{\tau}*{\boldsymbol{k}}^\dagger
\widehat{w[\cdot]}*{\boldsymbol{k}}
+\mathrm{c.c.}
]

得到该过程对 EAPE/BC EKE 转换的贡献。

该分解有助于判断：

* 大尺度 BC EKE 转回 EAPE 是否主要来自 BC 逆级串；
* 是否由底摩擦诱导；
* 是否与背景流对 BT 涡度的平流有关；
* 是否与 (\beta) 项直接有关。

---

## 12. 可选：相干性诊断

若要判断两个能量传递项是否在空间结构上相关，可以计算谱相干性：

[
C(A,B;\boldsymbol{k})
=====================

\frac{
\left|\sum_{n=1}^{N}A_n(\boldsymbol{k})B_n^\dagger(\boldsymbol{k})\right|^2
}{
\left(\sum_{n=1}^{N}|A_n(\boldsymbol{k})|^2\right)
\left(\sum_{n=1}^{N}|B_n(\boldsymbol{k})|^2\right)
}.
]

可重点比较：

* (C(D_{BC},T^W))
* (C(R^{KE}_{BC},T^W))
* (C(R^{KE}*{BC},T^N*{BT\rightarrow BC}))
* (C(D_{BC},T^N_{BT\rightarrow BC}))

这些相干性可帮助判断大尺度 BC EKE 是被 barotropization 转为 BT EKE，还是被 (T^W) 转回 EAPE。

---

## 13. 物理解读建议

### 13.1 弱底摩擦、无 (\beta)

典型能量循环是 BC-BT dual cascade：

1. 大尺度 (P) 产生 EAPE；
2. EAPE 前向级串到变形尺度；
3. (T^W<0)，EAPE 转为 BC EKE；
4. BC EKE 通过 barotropization 转为 BT EKE；
5. BT EKE 逆级串到大尺度；
6. BT EKE 由底摩擦耗散。

### 13.2 强底摩擦、无 (\beta)

大尺度可能出现 (T^W>0)，即 BC EKE 转回 EAPE。此时大尺度 EAPE 来源可由摩擦诱导的 BC energization 支持，但不一定形成闭合的纯 BC 能量环。

### 13.3 有 (\beta) 且形成 zonal jets

需要特别关注二维谱结构：

* 很多能量转换项会在 (k_x=0)、(k_y) 对应 jet-spacing scale 的位置出现尖峰；
* (P) 往往不属于这个 (k_x=0) 能量循环；
* (T^W>0) 可成为大尺度 EAPE 的主源；
* BC EKE 逆级串可把能量送到大尺度，再经 (T^W) 转回 EAPE，形成闭合的纯斜压能量环。

因此，在有喷流的 2LQG 系统中，不应只看各向同性谱。必须同时分析 ((k_x,k_y)) 二维谱。

---

## 14. 最小输出清单

一个完整的诊断程序至少应输出：

1. (\widehat{\mathrm{EAPE}}(k_x,k_y))
2. (\widehat{\mathrm{EKE}}_{BC}(k_x,k_y))
3. (\widehat{\mathrm{EKE}}_{BT}(k_x,k_y))
4. (P)
5. (T^W)
6. (R^{PE}_{BC})
7. (T^L_{BT\rightarrow BC})
8. (T^N_{BT\rightarrow BC})
9. (R^{KE}_{BC})
10. (D_{BC})
11. (T^L_{BC\rightarrow BT})
12. (T^N_{BC\rightarrow BT})
13. (R_{BT})
14. (D_{BT})
15. budget residuals
16. isotropic shell-integrated spectra
17. two-dimensional spectra
18. optional (T^W) decomposition
19. optional coherence diagnostics

---

## 15. 编码注意事项

1. 保持 FFT 归一化一致，否则能量和 tendency 的量纲会错。
2. 所有 nonlinear Jacobian 建议用伪谱法计算，并进行 dealiasing。
3. 所有谱预算项都要取 (2\operatorname{Re})。
4. 不要只保存 isotropic spectra；有 (\beta) 或 jets 时，二维谱非常关键。
5. (w_{3/2}) 应从 omega 方程诊断，不建议用界面高度方程的时间差分反推。
6. 对每个 snapshot 可独立诊断完整预算，但最终图应使用统计稳态窗口时间平均。
7. 必须检查预算闭合误差，尤其是：

   * (R) 项全波数积分是否为 0；
   * (T^W) 是否在 EAPE 与 BC EKE 预算之间抵消；
   * 线性 BT/BC 转换是否逐波数抵消；
   * 非线性 BT/BC 转换是否全局抵消。
8. 如果使用对称背景流 (U_1=+U, U_2=-U)，请重新核对 (S)、背景流项和热风关系。最安全的做法是从自己模型的原始两层方程重新推导 (w_{3/2}) 的 RHS，而不是直接套用单一记号。
