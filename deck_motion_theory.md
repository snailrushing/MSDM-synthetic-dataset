# 甲板运动仿真原理与数学推导

> 本文档描述用 MSS（Marine Systems Simulator）对 **Tanker** 船在 **6 级海况、135° 斜浪、固定航向 30 节**工况下，仿真甲板上某一点（`r_p = [-123, 0, -5]`）的四自由度运动（升沉、升沉速度、横摇、纵摇）的完整原理与数学推导。
> 所有公式均与 MSS 源码实现一致，具体源文件在正文中引用。

---

## 目录

1. [概述](#1-概述)
2. [坐标系与符号约定](#2-坐标系与符号约定)
3. [运动学方程](#3-运动学方程)
4. [6-DOF 动力学方程](#4-6-dof-动力学方程)
5. [波浪建模与一阶波浪力](#5-波浪建模与一阶波浪力)
6. [功率等效水动力模型](#6-功率等效水动力模型)
7. [航向自动舵](#7-航向自动舵)
8. [数值积分（RK4）](#8-数值积分rk4)
9. [甲板点运动转换（核心）](#9-甲板点运动转换核心)
10. [本工况参数](#10-本工况参数)
11. [参考文献](#11-参考文献)

---

## 1. 概述

甲板运动仿真的完整流程为：

1. **建模**：建立船体 6 自由度刚体运动方程（`hydroVessel.m`），水动力系数由势流计算（WAMIT）给出，并通过"功率等效"方法转化为常值矩阵（`computeManeuveringModel.m`）。
2. **激励**：用波浪谱合成不规则波，通过一阶波浪力 RAO 计算作用在船体上的 6 自由度广义波浪力（`waveForceRAO.m`）。
3. **控制**：航向自动舵锁定艏向，恒定纵向推力维持航速（`headingAutopilot.m`）。
4. **积分**：用四阶龙格-库塔（RK4）积分 12 维状态方程，得到船体参考点（原点 CO）处的 6 自由度运动 `η, ν`。
5. **甲板点换算**：用力臂变换，把 CO 处的运动换算到甲板点 `r_p`，输出升沉、升沉速度、横摇、纵摇。

其中第 5 步是本项目"甲板运动"区别于普通船舶仿真的关键。

---

## 2. 坐标系与符号约定

### 2.1 NED 与体坐标 FSD

MSS 使用两套右手坐标系：

- **NED（North-East-Down）惯性坐标系**：

  $$
  \eta = [x,\ y,\ z,\ \phi,\ \theta,\ \psi]^T
  $$

  其中 $x$ 北向、$y$ 东向、$z$ **向下为正**，$\phi,\theta,\psi$ 为欧拉角。
- **体坐标（body-fixed）FSD（Front-Starboard-Down）**：

  $$
  \nu = [u,\ v,\ w,\ p,\ q,\ r]^T
  $$

  $u,v,w$ 为线速度（m/s），$p,q,r$ 为角速度（rad/s）。轴向约定：

  | 轴      | 正方向   |
  | ------- | -------- |
  | $x_b$ | 指向船首 |
  | $y_b$ | 指向右舷 |
  | $z_b$ | 指向下方 |

  角速度 $p$（横摇率）、$q$（纵摇率）、$r$（艏摇率）分别绕 $x_b,y_b,z_b$，按右手定则。

### 2.2 欧拉角与旋转矩阵

采用 ZYX 欧拉角约定（先绕 $z$ 转 $\psi$，再绕 $y$ 转 $\theta$，最后绕 $x$ 转 $\phi$），体坐标到 NED 的旋转矩阵为（`Rzyx.m`）：

$$
R(\phi,\theta,\psi) = R_z(\psi)\,R_y(\theta)\,R_x(\phi)
$$

$$
R =
\begin{bmatrix}
c\psi\,c\theta & c\psi\,s\theta\,s\phi - s\psi\,c\phi & c\psi\,s\theta\,c\phi + s\psi\,s\phi \\
s\psi\,c\theta & s\psi\,s\theta\,s\phi + c\psi\,c\phi & s\psi\,s\theta\,c\phi - c\psi\,s\phi \\
-s\theta & c\theta\,s\phi & c\theta\,c\phi
\end{bmatrix}
$$

其中 $c(\cdot)=\cos(\cdot)$，$s(\cdot)=\sin(\cdot)$。

> **欧拉角符号约定（FSD）**：$\phi>0$ 为右舷下沉（right-hand rule 绕 $x_b$）；$\theta>0$ 为船首上抬；$\psi$ 为艏向角（相对北）。

### 2.3 Tanker 的体坐标定义

Tanker（WAMIT 数据）的体坐标原点 **CO** 位于：

$$
\text{CO} = \big(\underbrace{L_{pp}/2}_{\text{舯}},\ \underbrace{B/2}_{\text{中线}},\ \underbrace{\text{WL}}_{\text{静水面}}\big)
$$

即 $x=0$ 在舯、$y=0$ 在中线、$z=0$ 在静水面。主尺度和质量参数：

| 参数     | 符号                     | 值                                                            |
| -------- | ------------------------ | ------------------------------------------------------------- |
| 垂线间长 | $L_{pp}$               | 246 m                                                         |
| 型宽     | $B$                    | 46 m                                                          |
| 吃水     | $T$                    | 10 m                                                          |
| 排水体积 | $\nabla$               | ≈ 92,312 m³                                                 |
| 质量     | $m$                    | ≈ 94,620 t                                                   |
| 重心     | $r_g=[x_g,y_g,z_g]^T$  | $[3.93,\ 0,\ -2.5]^T$（FSD，即舯前 3.93 m、水线上方 2.5 m） |
| 惯性半径 | $k_{44},k_{55},k_{66}$ | $0.37B,\ 0.26L_{pp},\ 0.27L_{pp}$                           |

---

## 3. 运动学方程

位置/姿态对时间的导数为（`eulerang.m`）：

$$
\dot\eta = J(\eta)\,\nu
$$

$$
J(\eta) = \begin{bmatrix} R(\phi,\theta,\psi) & 0_{3\times3} \\ 0_{3\times3} & T_\Theta(\phi,\theta) \end{bmatrix}
$$

其中姿态变换矩阵（`Tzyx.m`）：

$$
T_\Theta =
\begin{bmatrix}
1 & s\phi\,\tan\theta & c\phi\,\tan\theta \\
0 & c\phi & -s\phi \\
0 & s\phi/\cos\theta & c\phi/\cos\theta
\end{bmatrix}
$$

线速度关系：$\dot\eta_{1:3} = R\,\nu_{1:3}$；角速度关系：$\dot\Theta = T_\Theta\,\nu_{4:6}$，其中 $\Theta=[\phi,\theta,\psi]^T$。

---

## 4. 6-DOF 动力学方程

MSS 的 12 状态非线性模型（`hydroVessel.m`）求解：

$$
\dot\eta = J(\eta)\,\nu
$$

$$
\dot\nu = \dot\nu_c + M^{-1}\Big[\tau_{\text{ext}} + \tau_{\text{crossflow}} + \tau_{\text{liftdrag}}
- \big(C_{RB}(\nu) + C_A(\nu_r) + D(\nu_r)\big)\nu_r - G\,\eta\Big]
$$

其中各符号定义如下。

### 4.1 系统惯性矩阵

$$
M = M_{RB} + M_A
$$

- $M_{RB}$：刚体质量矩阵（关于原点 CO，重心位于 $r_g=[x_g,y_g,z_g]^T$）：

$$
M_{RB} =
\begin{bmatrix}
m & 0 & 0 & 0 & m z_g & -m y_g \\
0 & m & 0 & -m z_g & 0 & m x_g \\
0 & 0 & m & m y_g & -m x_g & 0 \\
0 & -m z_g & m y_g & I_x & -I_{xy} & -I_{xz} \\
m z_g & 0 & -m x_g & -I_{yx} & I_y & -I_{yz} \\
-m y_g & m x_g & 0 & -I_{zx} & -I_{zy} & I_z
\end{bmatrix}
$$

  转动惯量 $I_x = m k_{44}^2,\ I_y = m k_{55}^2,\ I_z = m k_{66}^2$（Tanker 取 $I_{xy}=I_{xz}=I_{yz}=0$，主惯量轴与体轴对齐）。

- $M_A$：附加质量矩阵（由第 6 节的功率等效 $A_{eq}$ 给出）。

### 4.2 刚体科里奥利-向心力矩阵

MSS 采用"线性速度无关"表示（`rbody.m`），关于 CG 写作：

$$
C_{RB}(\nu_2) =
\begin{bmatrix}
m\,S(\nu_2) & 0_{3\times3} \\
0_{3\times3} & -S(I_G\,\nu_2)
\end{bmatrix},
\qquad
\nu_2 = [p,q,r]^T
$$

其中 $S(\cdot)$ 为斜对称矩阵（`Smtrx.m`）：

$$
S(a) = \begin{bmatrix} 0 & -a_3 & a_2 \\ a_3 & 0 & -a_1 \\ -a_2 & a_1 & 0 \end{bmatrix}
\quad\Longrightarrow\quad S(a)\,b = a\times b
$$

> 注：代码中 $C_{RB}$ 用 `rbody(...,r_bG=[0,0,0])` 计算（即假设 CG 与 CO 重合），这是对 CG 偏移的科里奥利项的近似处理；而 $M_{RB}$ 保留了完整的 CG 偏移。

### 4.3 附加质量科里奥利矩阵

由 `m2c.m` 计算。记 $M_A$ 的分块：

$$
M_A = \begin{bmatrix} M_{11} & M_{12} \\ M_{12}^T & M_{22} \end{bmatrix},
\qquad
\nu = \begin{bmatrix} \nu_1 \\ \nu_2 \end{bmatrix}
$$

令动量导数 $\dot\nu_1 = M_{11}\nu_1 + M_{12}\nu_2$，$\dot\nu_2 = M_{12}^T\nu_1 + M_{22}\nu_2$，则

$$
C_A(\nu) =
\begin{bmatrix}
0_{3\times3} & -S(\dot\nu_1) \\
-S(\dot\nu_1) & -S(\dot\nu_2)
\end{bmatrix}
$$

### 4.4 阻尼矩阵

总阻尼由势流等效阻尼、粘性阻尼和非线性项组成（`computeManeuveringModel.m` + `hydroVessel.m`）：

$$
D(\nu_r) = B_{eq} + B_v + \Delta D_{\text{nl}}(\nu_r)
$$

非线性对角线项：

$$
\begin{aligned}
D(1,1) &= D(1,1)\,e^{-k_u |u_r|} - X_{uu}\,|u_r| & \text{(纵荡)} \\
D(4,4) &= D(4,4)\,(1+\kappa_4 |p|) & \text{(横摇)} \\
D(5,5) &= D(5,5)\,(1+\kappa_5 |q|) & \text{(纵摇)} \\
D(6,6) &= D(6,6)\,(1+\kappa_6 |r|) & \text{(艏摇)}
\end{aligned}
$$

其中 $X_{uu}$ 由 ITTC-1957 摩擦阻力给出（`XuuITTC.m`）：

$$
X_{uu} = -\tfrac12\,\rho\,S\,(1+k)\,C_f,
\qquad
C_f = \frac{0.075}{(\log_{10}Re - 2)^2},
\qquad
Re = \frac{L\,|u_r|}{\nu_{\text{kin}}}
$$

湿表面积用 Mumford 公式：$S = 1.025\,L\,(C_B B + 1.7 T)$。

### 4.5 恢复力矩阵

$$
G = \begin{bmatrix} 0 & 0 \\ 0 & G_{3:5,3:5} \end{bmatrix},
\qquad
G_{3:5,3:5} = C(3{:}5,\,3{:}5,1)
$$

即恢复力只作用于垂荡、横摇、纵摇三个自由度：

$$
\begin{aligned}
G_{33} &= \rho g A_{wp} & \text{(水线面面积 → 垂荡恢复)} \\
G_{44} &= \rho g \nabla\,GM_T & \text{(横稳心 → 横摇恢复)} \\
G_{55} &= \rho g \nabla\,GM_L & \text{(纵稳心 → 纵摇恢复)}
\end{aligned}
$$

恢复力项为 $-G\eta$：$F_z = -G_{33} z$，$M_x = -G_{44}\phi$，$M_y = -G_{55}\theta$。

### 4.6 海流

无旋水平海流（`hydroVessel.m`）：

$$
\nu_c = \big[V_c\cos(\beta_c-\psi),\ V_c\sin(\beta_c-\psi),\ 0,\ 0,\ 0,\ 0\big]^T
$$

相对速度 $\nu_r = \nu - \nu_c$，其导数为

$$
\dot\nu_c = \begin{bmatrix} -S(\nu_{4:6})\,v_c \\ 0_{3\times1} \end{bmatrix},
\qquad
v_c = [V_c\cos(\beta_c-\psi),\ V_c\sin(\beta_c-\psi),\ 0]^T
$$

### 4.7 外力

$$
\tau_{\text{ext}} = \tau_{\text{control}} + \tau_{\text{wave1}}
$$

- $\tau_{\text{control}}$：控制力（航向自动舵 + 恒定推力，见第 7 节）。
- $\tau_{\text{wave1}}$：一阶波浪力（见第 5.5 节）。
- $\tau_{\text{crossflow}},\ \tau_{\text{liftdrag}}$：横流阻力、升阻力项（Tanker 用 Hoerner 横流阻力，见 `hydroVessel.m`）。

---

## 5. 波浪建模与一阶波浪力

### 5.1 波浪谱（JONSWAP）

本工况用 JONSWAP 谱（`waveSpectrum.m` case 7）：

$$
S(\omega) = \big[1 - 0.287\ln\gamma\big]\,
\alpha\, g^2\, \omega^{-5}\,
\exp\!\Big(-\tfrac54 (\omega_0/\omega)^4\Big)\,
\gamma^{\exp\!\big(-\frac{(\omega-\omega_0)^2}{2\sigma^2\omega_0^2}\big)}
$$

$$
\alpha = \frac{5}{16}\frac{H_s^2\,\omega_0^4}{g^2},
\qquad
\sigma = \begin{cases} 0.07, & \omega < \omega_0 \\ 0.09, & \omega \ge \omega_0 \end{cases},
\qquad
\gamma = 3.3
$$

谱峰频率与周期的关系：

$$
\omega_0 = \frac{2\pi}{T_0},
\qquad
T_z = 0.710\,T_0
$$

### 5.2 方向扩散

无扩散（长峰波）时方向函数 $M(\mu)=1$，$\mu=0$。有扩散时（本工况关闭）用 $\cos^2$ 扩散：

$$
M(\mu) = \frac{2}{\pi}\cos^2\mu,\quad \mu\in[-\pi/2,\pi/2]
$$

方向谱：$S_M(\omega,\mu) = S(\omega)\,M(\mu)$。

### 5.3 波浪高程合成

把频谱离散成 $N$ 个频率区间、$J$ 个方向，每分量幅值：

$$
A_{kj} = \sqrt{2\,S_M(\omega_k,\mu_j)\,\Delta\omega\,\Delta\mu}
$$

不规则波高程为各分量余弦叠加（`waveForceRAO.m`）：

$$
\zeta(t) = \sum_k \sum_j A_{kj}\,\cos\big(\omega_{e,kj}\,t + \varepsilon_{kj}\big)
$$

其中 $\varepsilon_{kj}\in[0,2\pi)$ 为随机相位。

### 5.4 遭遇频率

考虑航速与浪向的多普勒频移，遭遇频率为：

$$
\omega_e = \Big|\omega - \frac{\omega^2}{g}\,U\cos\beta_r\Big|
$$

其中 $U=\sqrt{u^2+v^2}$ 为船速，$\beta_r = \beta_{\text{wave}} - \psi$ 为相对浪向角（$\beta_r=0$ 顺浪、$90°$ 横浪、$180°$ 迎浪）。

### 5.5 一阶波浪力（力 RAO）

一阶波浪激励力由**力传递函数（力 RAO）**描述：对频率为 $\omega$、浪向为 $\beta$ 的单位波幅规则波，作用在第 $i$ 个自由度的广义波浪力是复数

$$
X_i(\omega,\beta) = |X_i(\omega,\beta)|\,e^{j\delta_i(\omega,\beta)} = \Re_i(\omega,\beta) + j\,\Im_i(\omega,\beta)
$$

- $|X_i|$：单位波幅对应的力/力矩幅值（纵荡/横荡/垂荡为力，单位 N/m；横摇/纵摇/艏摇为力矩，单位 N·m/m）；
- $\delta_i$：波浪力相对入射波高程的相位差，由波浪绕射效应决定；
- $\Re_i,\Im_i$：复数力 RAO 的实部与虚部，$\Re_i=|X_i|\cos\delta_i$，$\Im_i=|X_i|\sin\delta_i$。

该传递函数由线性势流边界元软件 WAMIT 在频域求解绕射边值问题得到，输出离散频率与浪向点上的幅值 $|X_i|$ 与相位 $\delta_i$（原始输出 `tanker.3`）；经 `wamit2vessel.m` 解析、量纲化并变换到体坐标 FSD 后存入 `tanker.mat` 的 `vessel.forceRAO`，再由 `waveInitialization.m` 按上式分解为实部/虚部并插值到波浪谱频率网格（存入 `waveModel.RAO_re` / `RAO_im`）。

时域 6 自由度广义波浪力为各谱分量的线性叠加（`waveForceRAO.m`，Fossen 2027 Eq. 10.96）：

$$
\tau_{\text{wave1},i}(t)
= \sum_k\sum_j A_{kj}\Big[
\Re_i(\omega_k,\beta_{r,kj})\cos(\omega_{e,kj}t+\varepsilon_{kj})
- \Im_i(\omega_k,\beta_{r,kj})\sin(\omega_{e,kj}t+\varepsilon_{kj})
\Big],
\qquad i=1,\dots,6
$$

即 $\tau_{\text{wave1},i}=\sum_k\sum_j A_{kj}|X_i|\cos(\omega_{e,kj}t+\varepsilon_{kj}+\delta_i)$，其中 $\delta_i=\mathrm{atan2}(\Im_i,\Re_i)$。

> 区分：本文用**力 RAO**（波浪激励力传递函数）。**运动 RAO**（`waveMotionRAO.m`）是船体运动响应传递函数，可跳过动力学方程直接合成波频运动 $\eta_w,\nu_w,\dot\nu_w$：
> $\eta_{w,i}=\sum_k\sum_j A_{kj}|H_i|\cos(\omega_e t+\angle H_i+\varepsilon)$。本项目采用力 RAO + 时域积分的完整路线。

---

## 6. 功率等效水动力模型

WAMIT 给出的附加质量 $A(\omega)$ 与势流阻尼 $B(\omega)$ 是频率相关的。为得到常值模型，用归一化波谱做加权平均（`computeManeuveringModel.m`）：

$$
A_{eq} = \int_0^\infty A(\omega)\,S_N(\omega)\,d\omega,
\qquad
B_{eq} = \int_0^\infty B(\omega)\,S_N(\omega)\,d\omega
$$

其中归一化谱 $S_N(\omega) = S(\omega)/m_0$，$m_0 = \int_0^\infty S(\omega)d\omega$ 为零阶谱矩。

粘性阻尼矩阵 $B_v$（对角）：

- 无约束模态（纵荡、横荡、艏摇，DOF 1,2,6）：

$$
B_v(i,i) = \kappa_i\,B_{eq}(i,i),\qquad i=1,2,6
$$

- 有恢复模态（垂荡、横摇、纵摇，DOF 3,4,5）：

$$
B_v(i,i) = 2\,\Delta\zeta_i\,\sqrt{M(i,i)\,G(i,i)},\qquad i=3,4,5
$$

总阻尼 $D = B_{eq} + B_v$，惯性矩阵 $M = M_{RB} + A_{eq}$。

---

## 7. 航向自动舵

SISO PID 极点配置算法（`headingAutopilot.m`，Fossen 2027 Algorithm 15.1）：

$$
\tau_N = -K_p\,\operatorname{ssa}(\psi-\psi_d) - K_i\,z_{int} - K_d\,(r - r_d)
$$

$$
K_p = M_{66}\,\omega_n^2,\qquad
K_d = M_{66}\,2\zeta\omega_n,\qquad
K_i = \tfrac{1}{10}K_p\,\omega_n
$$

积分状态：$\dot z_{int} = \operatorname{ssa}(\psi-\psi_d)$；$\operatorname{ssa}(\cdot)$ 为最小符号角（映射到 $[-\pi,\pi]$）。参考信号 $\psi_d$ 由三阶参考模型生成（比 $\omega_n$ 慢 20 倍），并受最大转艏率 $r_{max}$ 约束。

控制力向量：

$$
\tau_{\text{control}} = [\tau_X,\ 0,\ 0,\ 0,\ 0,\ \tau_N]^T
$$

其中 $\tau_X$ 为恒定纵向推力，用于维持航速。

---

## 8. 数值积分（RK4）

对 12 维状态 $x=[\nu^T,\ \eta^T]^T$，方程 $\dot x = f(x,\tau_{\text{ext}})$ 用四阶龙格-库塔积分（`rk4.m`）：

$$
x_{k+1} = x_k + \frac{h}{6}\big(k_1 + 2k_2 + 2k_3 + k_4\big)
$$

$$
\begin{aligned}
k_1 &= f(x_k) \\
k_2 &= f(x_k + \tfrac{h}{2}k_1) \\
k_3 &= f(x_k + \tfrac{h}{2}k_2) \\
k_4 &= f(x_k + h\,k_3)
\end{aligned}
$$

---

## 9. 甲板点运动转换（核心）

船体方程给出的 $\eta,\nu$ 是原点 CO（舯/中线/静水面）处的运动。甲板上一点 $P$ 在体坐标中固定于力臂 $r_p=[x_p,y_p,z_p]^T$（从 CO 量起），其运动由 CO 的运动加上杠杆效应得到。

### 9.1 位置

$P$ 在 NED 中的位置：

$$
p_{NED}(t) = \eta_{1:3}(t) + R\big(\phi(t),\theta(t),\psi(t)\big)\,r_p
$$

### 9.2 速度

$P$ 在体坐标中的速度（CO 速度 + 转动引起的牵连速度）：

$$
v_p^b = \nu_{1:3} + S(\nu_{4:6})\,r_p = \nu_{1:3} + \omega \times r_p
$$

变换到 NED：$\dot p_{NED} = R\,v_p^b$。显式地，

$$
\dot p_{NED} = R\,\nu_{1:3} + R\,S(\omega)\,r_p
$$

### 9.3 加速度

利用 $\dot R = R\,S(\omega)$，对位置求二阶导：

$$
\ddot p_{NED} = R\Big[\dot\nu_{1:3} + S(\omega)\,\nu_{1:3} + S(\dot\omega)\,r_p + S(\omega)^2\,r_p\Big]
$$

其中 $S(\omega)^2 r_p = \omega\times(\omega\times r_p)$ 为向心加速度项。

### 9.4 小角度线性化

小角度下 $R \approx I + S(\Theta)$（$\Theta=[\phi,\theta,\psi]^T$），于是

$$
(R - I)\,r_p \approx \Theta \times r_p
$$

其第 3 分量（垂向）：

$$
(\Theta\times r_p)_3 = \phi\,y_p - \theta\,x_p
$$

故甲板点垂向位移（向上为正，相对初始位置）：

$$
\boxed{\,z_p = -\eta(3) - \phi\,y_p + \theta\,x_p\,}
$$

即

$$
z_p = \text{heave}_{CO} - x_p\,\theta + y_p\,\phi
$$

其中 $\text{heave}_{CO}=-\eta(3)$ 为 CO 处升沉（向上为正）。

对本工况 $r_p=[-123,\ 0,\ -5]$（$x_p=-123,\ y_p=0$）：

$$
\boxed{z_p = \text{heave}_{CO} - 123\,\theta}
$$

可见**舯后 123 m 处，纵摇被 123 m 的力臂放大**：即使纵摇只有 $0.05$ rad（约 2.9°），也会在甲板点产生 $123\times0.05 = 6.15$ m 的垂向运动——这是甲板点升沉的主导项。

### 9.5 四自由度输出

| 输出量   | 符号         | 定义                               | 单位            |
| -------- | ------------ | ---------------------------------- | --------------- |
| 升沉     | $z_p$      | $-\big[p_{NED,3}(t)-r_p(3)\big]$ | m（向上为正）   |
| 升沉速度 | $\dot z_p$ | $-\dot p_{NED,3}(t)$             | m/s（向上为正） |
| 横摇     | $\phi$     | $\eta(4)$                        | °              |
| 纵摇     | $\theta$   | $\eta(5)$                        | °              |

> 注意横摇/纵摇是全船各处相同的姿态角，直接取 $\eta(4),\eta(5)$；只有升沉和升沉速度需要做上述力臂变换。

---

## 10. 本工况参数

| 参数          | 符号                    | 值                                                     |
| ------------- | ----------------------- | ------------------------------------------------------ |
| 船舶          | —                      | Tanker（WAMIT）                                        |
| 海况等级      | —                      | 6 级（$H_s = 4\sim6$ m）                             |
| 有义波高      | $H_s$                 | 5.0 m                                                  |
| 谱峰频率      | $\omega_0$            | 0.5 rad/s（$T_0\approx12.6$ s，$T_z\approx8.9$ s） |
| 波浪谱        | —                      | JONSWAP（$\gamma=3.3$，无方向扩散）                  |
| 浪向角        | $\beta_{\text{wave}}$ | 135°（斜浪）                                          |
| 海流          | $V_c$                 | 0                                                      |
| 航速          | $U$                   | 30 节 = 15.43 m/s                                      |
| 艏向          | $\psi$                | 0°（固定）                                            |
| 甲板点力臂    | $r_p$                 | $[-123,\ 0,\ -5]^T$ m                                |
| 仿真时长/步长 | $T_f,\ h$             | 600 s / 0.05 s                                         |
| 恒定推力      | $\tau_X$              | $2.4\times10^6$ N                                    |

---

## 11. 参考文献

1. Fossen, Thor. I. (2027). *Handbook of Marine Craft Hydrodynamics and Motion Control*, 3rd ed., John Wiley & Sons Ltd.
2. Fossen, T. I. (2025). *Maneuvering Coefficient Estimation from Frequency-Dependent Added Mass and Damping: A Power-Based Approach*. Ocean Engineering, 341, 122494.
3. DNV-RP-C205 (2025). *Environmental Conditions and Environmental Loads*. Recommended Practice, DNV, Høvik, Norway.
4. Perez, T. & Fossen, T. I. (2009). *A Matlab Tool for Parametric Identification of Radiation-Force Models of Ships and Offshore Structures*. MIC-30(1):1-15.
5. WAMIT User Manual (Chapter 4: body coordinate system, mass matrix, RAO conventions).

---

## 附录：与本源码的对应关系

| 公式/模块                                          | 源文件                                                 |
| -------------------------------------------------- | ------------------------------------------------------ |
| 旋转矩阵$R$、姿态矩阵 $T_\Theta$、斜对称 $S$ | `LIBRARY/kinematics/Rzyx.m`, `Tzyx.m`, `Smtrx.m` |
| 运动学$J(\eta)$                                  | `LIBRARY/kinematics/eulerang.m`                      |
| 6-DOF 动力学                                       | `CRAFT/hydroVessel.m`                                |
| 刚体质量/科里奥利$M_{RB},C_{RB}$                 | `LIBRARY/modeling/rbody.m`, `Hmtrx.m`              |
| 附加质量科里奥利$C_A$                            | `LIBRARY/modeling/m2c.m`                             |
| ITTC 二次阻力$X_{uu}$                            | `LIBRARY/modeling/XuuITTC.m`                         |
| 功率等效$A_{eq},B_{eq},B_v$                      | `HYDRO/computeManeuveringModel.m`                    |
| 波浪谱                                             | `LIBRARY/environment/waveSpectrum.m`                 |
| 方向谱                                             | `LIBRARY/environment/waveDirectionalSpectrum.m`      |
| 一阶波浪力                                         | `LIBRARY/environment/waveForceRAO.m`                 |
| 航向自动舵                                         | `GNC/headingAutopilot.m`                             |
| RK4 积分                                           | `LIBRARY/numericalMethods/rk4.m`                     |
| 甲板点换算 + CSV 输出                              | `SIMtankerDeckMotion.m`                              |
