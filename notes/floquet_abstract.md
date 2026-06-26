
# 可直接实现的矩阵形式

## 线性化方程

原波动方程为
\[
(q_1)_t+U\big((q_1)_x+k_D^2(\psi_1)_x\big)+[\psi_1,q_1]-\overline{[\psi_1,q_1]}=\Delta q_1,
\tag{1}
\]
\[
(q_2)_t-U\big((q_2)_x+k_D^2(\psi_2)_x\big)+[\psi_2,q_2]-\overline{[\psi_2,q_2]}=\Delta q_2.
\tag{2}
\]
注意原方程的 $U$ 满足一个代数关系
\[
  \overline{(\psi_1)_xq_1}-k_D^2 U = F_1
\]
其中 $\psi_1$ 和 $q_1$ 由背景态和扰动态联合确定，这使得 $U$ 极容易发生围绕前述饱和值的偏移，变为 $U+\delta U$。因此$\delta U$的具体表达式可由各假设推出来。

令
\[
\psi_i=\psi_i^{(0)}+\tilde\psi_i,\qquad q_i=q_i^{(0)}+\tilde q_i.
\]
保留一阶项，得

\[
(\tilde q_1)_t+U\big((\tilde q_1)_x+k_D^2(\tilde\psi_1)_x\big)
+\delta U\big((q_1^{(0)})_x+k_D^2(\psi_1^{(0)})_x\big)
+[\psi_1^{(0)},\tilde q_1]+[\tilde\psi_1,q_1^{(0)}]
-\overline{[\psi_1^{(0)},\tilde q_1]+[\tilde\psi_1,q_1^{(0)}]}
=\Delta\tilde q_1,
\tag{L1}
\]

\[
(\tilde q_2)_t-U\big((\tilde q_2)_x+k_D^2(\tilde\psi_2)_x\big)
-\delta U\big((q_2^{(0)})_x+k_D^2(\psi_2^{(0)})_x\big)
+[\psi_2^{(0)},\tilde q_2]+[\tilde\psi_2,q_2^{(0)}]
-\overline{[\psi_2^{(0)},\tilde q_2]+[\tilde\psi_2,q_2^{(0)}]}
=\Delta\tilde q_2.
\tag{L2}
\]

并且
\[
\tilde q_1=\Delta \tilde\psi_1+\frac12 k_D^2(\tilde\psi_2-\tilde\psi_1),
\qquad
\tilde q_2=\Delta \tilde\psi_2+\frac12 k_D^2(\tilde\psi_1-\tilde\psi_2).
\tag{L3}
\]

## 基准态

基准态
\[
\psi_1^{(0)}=\Re(A_1 e^{i\theta}),\qquad
\psi_2^{(0)}=\Re(A_2 e^{i\theta}),\qquad
\theta=kx+ly.
\]

记
\[
\psi_i^{(0)}=\frac12\Big(A_i e^{i\theta}+A_i^* e^{-i\theta}\Big).
\]

由于
\[
q_i=\Delta\psi_i+\frac12k_D^2(\psi_{3-i}-\psi_i),
\]
对单 Fourier 模 \(e^{i\theta}\) 有
\[
\Delta e^{i\theta}=-(k^2+l^2)e^{i\theta}.
\]
记
\[
\kappa^2:=k^2+l^2.
\]
则基准态 PV 的复振幅为
\[
Q_1^+ = -\Big(\kappa^2+\frac12k_D^2\Big)A_1+\frac12k_D^2 A_2,\qquad
Q_2^+ = -\Big(\kappa^2+\frac12k_D^2\Big)A_2+\frac12k_D^2 A_1.
\]
于是
\[
q_1^{(0)}=\frac12\Big(Q_1^+ e^{i\theta}+(Q_1^+)^*e^{-i\theta}\Big),\qquad
q_2^{(0)}=\frac12\Big(Q_2^+ e^{i\theta}+(Q_2^+)^*e^{-i\theta}\Big).
\]

## 饱和态振幅 $A_i$ 的定法

给定
\[
S=\frac{k_D^2+k^2+l^2}{k_D^2-k^2-l^2},\qquad
U=\sqrt{\frac{k_D^2+k^2+l^2}{k_D^2-k^2-l^2}}\frac{k^2+l^2}{k}
=\sqrt{S}\,\frac{k^2+l^2}{k},
\]
以及稳态模振幅关系
\[
A_1(1+iS)+A_2(1-iS)=0.
\]
故
\[
A_2=-A_1\frac{1+iS}{1-iS}.
\]
再用固定 \(F\) 的 PV flux 条件
\[
\frac{i k k_D^2}{8}\left(A_1A_2^*-A_1^*A_2\right)-k_D^2U=F
\]
定出 \(A_1\) 的模长。不妨让 \(A_1\in\mathbb R\)，可得
\[
A_1A_2^*-A_1^*A_2
= A_1^2\!\left(\frac{-1+iS}{1+iS}-\frac{-1-iS}{1-iS}\right)
= \frac{4iS}{1+S^2}A_1^2.
\]
于是
\[
F = -\frac{k k_D^2 S}{2(1+S^2)}A_1^2-k_D^2U.
\]
故
\[
A_1^2 = -\,\frac{2(1+S^2)}{k\,k_D^2\,S}\,\bigl(F+k_D^2U\bigr),
\]
该算式右端项中的 $F$ 会是一个比较大的负数，而 $U$ 是不随之而动的定值，这保证了它一定是正数，因此可以在实数范围内开根号。程序中取实正根
\[
A_1=\sqrt{A_1^2},\qquad
A_2=-A_1\frac{1+iS}{1-iS}.
\]

## Floquet/Bloch 展开

取
\[
\tilde\psi_1=e^{\lambda t}\sum_{n=-N}^N a_n e^{i\phi_n},
\qquad
\tilde\psi_2=e^{\lambda t}\sum_{n=-N}^N b_n e^{i\phi_n},
\]
其中
\[
\phi_n=(p+nk)x+(q+nl)y.
\]
记第 \(n\) 个 Floquet 波矢为
\[
\alpha_n:=p+nk,\qquad \beta_n:=q+nl,
\qquad m_n^2:=\alpha_n^2+\beta_n^2.
\]

于是
\[
\Delta e^{i\phi_n}=-m_n^2 e^{i\phi_n},\qquad
\partial_x e^{i\phi_n}=i\alpha_n e^{i\phi_n}.
\]

由反演关系 (L3)，定义每个 \(n\) 上的 PV 振幅
\[
\hat q_{1,n}=c_n a_n+d\, b_n,
\qquad
\hat q_{2,n}=d\, a_n+c_n b_n,
\]
其中
\[
c_n:=-m_n^2-\frac12 k_D^2,\qquad d:=\frac12 k_D^2.
\]
即
\[
\tilde q_1=e^{\lambda t}\sum_n \hat q_{1,n}e^{i\phi_n},\qquad
\tilde q_2=e^{\lambda t}\sum_n \hat q_{2,n}e^{i\phi_n}.
\]

这一步非常关键，因为之后所有矩阵元都可以直接写成 \(a_n,b_n\) 的线性组合。




## 广义特征值问题


令
\[
X=(a_{-N},b_{-N},a_{-N+1},b_{-N+1},\dots,a_N,b_N)^T.
\]
广义特征值问题
\[
(\lambda M+K)X=0
\]
 \(K\) 是块三对角矩阵：

\[
K=
\begin{pmatrix}
K_{-N}^0 & K_{-N}^+ & & & \\
K_{-N+1}^- & K_{-N+1}^0 & K_{-N+1}^+ & & \\
& \ddots & \ddots & \ddots & \\
& & K_{N-1}^- & K_{N-1}^0 & K_{N-1}^+ \\
& & & K_N^- & K_N^0
\end{pmatrix},
\]
而
\[
M=\operatorname{diag}(M_{-N},M_{-N+1},\dots,M_N).
\]

块矩阵 \(M\)仅主对角非零，主对角第 \(m\) 块
\[
M_m=
\begin{pmatrix}
-\sigma_m^2-\frac12k_D^2 & \frac12k_D^2\\
\frac12k_D^2 & -\sigma_m^2-\frac12k_D^2
\end{pmatrix}.
\]

块矩阵 \(K\)

第 \(m\) 行块满足
\[
K_{m,m-1}=K_m^-,
\qquad
K_{m,m}=K_m^0,
\qquad
K_{m,m+1}=K_m^+.
\]

其中
\[
\sigma_m^2=(p+mk)^2+(q+ml)^2,\qquad
\alpha_m=p+mk,\qquad
\gamma=kq-lp,\qquad d=\frac12k_D^2.
\]

主对角块
\[
K_m^0=
\begin{pmatrix}
-\sigma_m^4-\frac12k_D^2\sigma_m^2+i\alpha_m U(-\sigma_m^2+\frac12k_D^2)
&
d(\sigma_m^2+i\alpha_m U)
\\[1mm]
d(\sigma_m^2-i\alpha_m U)
&
-\sigma_m^4-\frac12k_D^2\sigma_m^2-i\alpha_m U(-\sigma_m^2+\frac12k_D^2)
\end{pmatrix}.
\]

下对角块
\[
K_m^-=
-\frac{\gamma}{2}
\begin{pmatrix}
A_1 c_{m-1}-Q_1^+ & A_1 d\\
A_2 d & A_2 c_{m-1}-Q_2^+
\end{pmatrix},
\]
上对角块
\[
K_m^+=
\frac{\gamma}{2}
\begin{pmatrix}
A_1^* c_{m+1}-(Q_1^+)^* & A_1^* d\\
A_2^* d & A_2^* c_{m+1}-(Q_2^+)^*
\end{pmatrix},
\]
其中
\[
c_n=-\sigma_n^2-\frac12k_D^2.
\]


### 4. Resonant corrections to the Floquet matrix

When the perturbation contains a zero-wavenumber component, the Jacobian terms and the induced mean-flow correction \(\delta U\) generate additional couplings that are not included in the generic nearest-neighbour block structure above. These terms must be added separately under the exact resonance conditions.

We denote
\[
\alpha_n = p + nk,\qquad \beta_n = q + nl,
\]
and recall
\[
\hat q_{1,n} = c_n a_n + d\, b_n,\qquad
\hat q_{2,n} = c_n b_n + d\, a_n,
\]
with
\[
c_n = -(\alpha_n^2+\beta_n^2)-d,\qquad d=\frac12 k_D^2.
\]

We also define
\[
R_1^+ = Q_1^+ + k_D^2 A_1,\qquad
R_2^+ = Q_2^+ + k_D^2 A_2.
\]

---

#### 4.1. Jacobian resonance correction

The generic Jacobian contribution
\[
J(\psi_i',q_i^{(0)})+J(\psi_i^{(0)},q_i')
\]
produces only the nearest-neighbour couplings \(n\leftrightarrow n\pm1\) as long as no zero-wavenumber mode is present. However, when one of the sideband wavenumbers vanishes, the coefficients multiplying the resonant harmonics must be reinterpreted as contributions to the corresponding Floquet modes.

There are two possible resonance conditions:

- backward resonance:
\[
\alpha_{n-1}=0,\qquad \beta_{n-1}=0,
\]
- forward resonance:
\[
\alpha_{n+1}=0,\qquad \beta_{n+1}=0.
\]

Under these conditions, the Jacobian terms contribute additional on-site corrections to the \(n\)-th Floquet equation.

For the first layer, the resonant correction is
\[
\Delta J_{1,n}^{(-)}
= -\frac12 (kq-lp)\,Q_1^+\, a_n
+\frac12 (kq-lp)\,A_1\, \hat q_{1,n},
\qquad \text{if } \alpha_{n-1}=\beta_{n-1}=0,
\]
and
\[
\Delta J_{1,n}^{(+)}
= -\frac12 (kq-lp)\,(Q_1^+)^*\, a_n
+\frac12 (kq-lp)\,A_1^*\, \hat q_{1,n},
\qquad \text{if } \alpha_{n+1}=\beta_{n+1}=0.
\]

Using \(\hat q_{1,n}=c_n a_n+d\,b_n\), these become
\[
\Delta J_{1,n}^{(-)}
=
\frac12 (kq-lp)\Big[(A_1c_n-Q_1^+)a_n + A_1 d\, b_n\Big],
\]
\[
\Delta J_{1,n}^{(+)}
=
\frac12 (kq-lp)\Big[(A_1^*c_n-(Q_1^+)^*)a_n + A_1^* d\, b_n\Big].
\]

Similarly, for the second layer,
\[
\Delta J_{2,n}^{(-)}
= -\frac12 (kq-lp)\,Q_2^+\, b_n
+\frac12 (kq-lp)\,A_2\, \hat q_{2,n},
\qquad \text{if } \alpha_{n-1}=\beta_{n-1}=0,
\]
\[
\Delta J_{2,n}^{(+)}
= -\frac12 (kq-lp)\,(Q_2^+)^*\, b_n
+\frac12 (kq-lp)\,A_2^*\, \hat q_{2,n},
\qquad \text{if } \alpha_{n+1}=\beta_{n+1}=0,
\]
that is,
\[
\Delta J_{2,n}^{(-)}
=
\frac12 (kq-lp)\Big[A_2 d\, a_n + (A_2c_n-Q_2^+) b_n\Big],
\]
\[
\Delta J_{2,n}^{(+)}
=
\frac12 (kq-lp)\Big[A_2^* d\, a_n + (A_2^*c_n-(Q_2^+)^*) b_n\Big].
\]

Therefore the Jacobian resonance correction modifies only the diagonal block \(K_{n,n}\), with
\[
\Delta K_{n,n}^{J,-}
=
\frac12 (kq-lp)
\begin{pmatrix}
A_1c_n-Q_1^+ & A_1 d\\
A_2 d & A_2c_n-Q_2^+
\end{pmatrix},
\]
when \(\alpha_{n-1}=\beta_{n-1}=0\), and
\[
\Delta K_{n,n}^{J,+}
=
\frac12 (kq-lp)
\begin{pmatrix}
A_1^*c_n-(Q_1^+)^* & A_1^* d\\
A_2^* d & A_2^*c_n-(Q_2^+)^*
\end{pmatrix},
\]
when \(\alpha_{n+1}=\beta_{n+1}=0\).

If neither resonance condition is satisfied, no extra Jacobian correction is present.

---

#### 4.2. Mean-flow correction \(\delta U\)

The perturbation also induces a correction to the zonal mean flow,
\[
\delta U
=
\frac{ik}{2k_D^2}
\left\{
[A_1\hat q_{1,n}-a_nQ_1^+]\,\delta(\alpha_{n+1},\beta_{n+1})
+
[-A_1^*\hat q_{1,n}+a_n(Q_1^+)^*]\,\delta(\alpha_{n-1},\beta_{n-1})
\right\}.
\]

This mean-flow correction enters the two perturbation equations through
\[
+\delta U\big[(q_1^{(0)})_x+k_D^2(\psi_1^{(0)})_x\big],
\qquad
-\delta U\big[(q_2^{(0)})_x+k_D^2(\psi_2^{(0)})_x\big].
\]
Using
\[
q_i^{(0)}=\frac12(Q_i^+e^{i\theta}+\text{c.c.}),\qquad
\psi_i^{(0)}=\frac12(A_i e^{i\theta}+\text{c.c.}),\qquad
\theta=kx+ly,
\]
we obtain
\[
(q_i^{(0)})_x+k_D^2(\psi_i^{(0)})_x
=
\frac{ik}{2}R_i^+ e^{i\theta}
-\frac{ik}{2}(R_i^+)^* e^{-i\theta}.
\]

Because \(\delta U\) is itself generated by the \(n\)-th perturbation mode under a resonance condition, it produces couplings from the \(n\)-th column to two different Floquet rows.

---

##### 4.2.1. Case \(\alpha_{n-1}=\beta_{n-1}=0\)

In this case,
\[
\delta U=\frac{ik}{2k_D^2} Z_n^-,
\qquad
Z_n^-=-A_1^*\hat q_{1,n}+a_n(Q_1^+)^*,
\]
or
\[
Z_n^-=\big[(Q_1^+)^*-A_1^*c_n\big]a_n-A_1^* d\, b_n.
\]

Since \(\phi_{n-1}=0\), the factors \(e^{\pm i\theta}\) correspond to the Floquet modes \(n\) and \(n-2\), respectively. Hence \(\delta U\) modifies the blocks \(K_{n,n}\) and \(K_{n-2,n}\).

The contribution to the \(n\)-th row block is
\[
\Delta K_{n,n}^{U,-}
=
\frac{k^2}{4k_D^2}
\begin{pmatrix}
- R_1^+\big[(Q_1^+)^*-A_1^*c_n\big]
&
+ R_1^+A_1^*d
\\[1ex]
+ R_2^+\big[(Q_1^+)^*-A_1^*c_n\big]
&
- R_2^+A_1^*d
\end{pmatrix}.
\]

The contribution to the \((n-2)\)-th row block is
\[
\Delta K_{n-2,n}^{U,-}
=
\frac{k^2}{4k_D^2}
\begin{pmatrix}
(R_1^+)^*\big[(Q_1^+)^*-A_1^*c_n\big]
&
-(R_1^+)^*A_1^*d
\\[1ex]
-(R_2^+)^*\big[(Q_1^+)^*-A_1^*c_n\big]
&
+(R_2^+)^*A_1^*d
\end{pmatrix}.
\]

---

##### 4.2.2. Case \(\alpha_{n+1}=\beta_{n+1}=0\)

In this case,
\[
\delta U=\frac{ik}{2k_D^2} Z_n^+,
\qquad
Z_n^+=A_1\hat q_{1,n}-a_nQ_1^+,
\]
or
\[
Z_n^+=(A_1c_n-Q_1^+)a_n + A_1 d\, b_n.
\]

Since \(\phi_{n+1}=0\), the factors \(e^{\pm i\theta}\) now correspond to the Floquet modes \(n+2\) and \(n\), respectively. Hence \(\delta U\) modifies the blocks \(K_{n,n}\) and \(K_{n+2,n}\).

The contribution to the \(n\)-th row block is
\[
\Delta K_{n,n}^{U,+}
=
\frac{k^2}{4k_D^2}
\begin{pmatrix}
(R_1^+)^*(A_1c_n-Q_1^+)
&
(R_1^+)^*A_1 d
\\[1ex]
-(R_2^+)^*(A_1c_n-Q_1^+)
&
-(R_2^+)^*A_1 d
\end{pmatrix}.
\]

The contribution to the \((n+2)\)-th row block is
\[
\Delta K_{n+2,n}^{U,+}
=
\frac{k^2}{4k_D^2}
\begin{pmatrix}
- R_1^+(A_1c_n-Q_1^+)
&
- R_1^+A_1 d
\\[1ex]
+ R_2^+(A_1c_n-Q_1^+)
&
+ R_2^+A_1 d
\end{pmatrix}.
\]

---

#### 4.3. Summary of the resonant block corrections

When the resonance option is enabled, the full Floquet matrix consists of the generic tridiagonal block structure plus the additional resonant corrections described above.

For a given mode \(n\),

- if \(\alpha_{n-1}=\beta_{n-1}=0\), then add
\[
\Delta K_{n,n}^{J,-}+\Delta K_{n,n}^{U,-}
\]
to the diagonal block, and add
\[
\Delta K_{n-2,n}^{U,-}
\]
to the off-diagonal block two rows below;

- if \(\alpha_{n+1}=\beta_{n+1}=0\), then add
\[
\Delta K_{n,n}^{J,+}+\Delta K_{n,n}^{U,+}
\]
to the diagonal block, and add
\[
\Delta K_{n+2,n}^{U,+}
\]
to the off-diagonal block two rows above.

If both conditions fail, no resonant correction is applied.

These terms are sparse and affect only isolated blocks associated with exact zero-wavenumber sidebands.