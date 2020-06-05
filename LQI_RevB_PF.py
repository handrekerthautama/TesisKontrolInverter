import numpy as np
import matplotlib.pyplot as plt
from Functions import ss, nl, l2h1, l2h2, m11, m12, m21, m22, pf, qf, il, transform
from matplotlib.ticker import (MultipleLocator)
from math import exp, sqrt
from numpy.linalg import inv


w = 2 * 3.142 * 50
L = 0.0008
C = 0.000075
t_init = 0.0
t_final = 4.01
t_delta = 0.001
r1 = 230
r2 = 0

t = np.arange(t_init, t_final, t_delta)
x = np.zeros((4, len(t)))
z = np.zeros((4, len(t)))
dx = np.zeros((4, len(t)))
y = np.zeros((2, len(t)))
v1 = np.zeros(len(t))
v2 = np.zeros(len(t))
u1 = np.zeros(len(t))
u2 = np.zeros(len(t))
iL = np.zeros((2, len(t)))
pqL = np.zeros((2, len(t)))
sOut = np.zeros(len(t))
sL = np.zeros(len(t))
PFOut = np.zeros(len(t))
PFL = np.zeros(len(t))
p = np.zeros(len(t))
q = np.zeros(len(t))
p_out = np.zeros(len(t))
q_out = np.zeros(len(t))
e = np.zeros((2, len(t)))
e_sum = np.zeros((2, len(t)))
pint = 0
ki = 2236.1
k1 = 365.3
k2 = 27

x[0, 0] = 0.24734592
x[1, 0] = 0.0
x[2, 0] = 0.00364869
x[3, 0] = 0.01165741
dx[0, 0] = 48.64925313
dx[1, 0] = 0.0
e[0, 0] = 269.44
e_sum[0, 0] = 5.3888
y[0, 0] = 0.24734592
y[1, 0] = 0.0

for i in range(len(t) - 1):
    if 1000 < i < 2000:
        x[2, i] = 0.28 * (1 - exp(-i * t_delta / 0.001)) + C * dx[0, i] - C * w * x[1, i]
        x[3, i] = -0.22 * (1 - exp(-i * t_delta / 0.001)) + C * dx[1, i] + C * w * x[0, 1]
    if 2000 < i < 3000:
        x[2, i] = 0.32 * (1 - exp(-i * t_delta / 0.001)) + C * dx[0, i] - C * w * x[1, i]
        x[3, i] = -0.29 * (1 - exp(-i * t_delta / 0.001)) + C * dx[1, i] + C * w * x[0, 1]
    if i > 3000:
        x[2, i] = 0.35 * (1 - exp(-i * t_delta / 0.001)) + C * dx[0, i] - C * w * x[1, i]
        x[3, i] = -0.17 * (1 - exp(-i * t_delta / 0.001)) + C * dx[1, i] + C * w * x[0, 1]
    # controller
    e[0, i] = r1 - x[0, i]
    e[1, i] = r2 - x[1, i]
    for k in range(i):
        e_sum[0, i] = e_sum[0, i] + e[0, k] * t_delta
        e_sum[1, i] = e_sum[1, i] + e[1, k] * t_delta
    v1[i] = ki * e_sum[0, i] - k1 * x[0, i] - k2 * dx[0, i]
    v2[i] = ki * e_sum[1, i] - k1 * x[1, i] - k2 * dx[1, i]
    # calculate pf & qf
    p[i] = pf(x[0, i], x[1, i], x[2, i], x[3, i], dx[0, i], dx[1, i], C)
    q[i] = qf(x[0, i], x[1, i], x[2, i], x[3, i], dx[0, i], dx[1, i], w, L, C)
    # calculate u1 & u2
    L2H1 = l2h1(x[0, i], x[1, i], x[2, i], x[3, i], dx[0, i], dx[1, i], p[i], q[i], w, L, C)
    L2H2 = l2h2(x[0, i], x[1, i], x[2, i], x[3, i], dx[0, i], dx[1, i], p[i], q[i], w, L, C)
    M11 = m11(x[0, i], x[1, i], x[2, i], w, L, C)
    M12 = m12(x[0, i], x[1, i], x[3, i], w, L, C)
    M21 = m21(x[0, i], x[1, i], x[2, i], w, L, C)
    M22 = m22(x[0, i], x[1, i], x[3, i], w, L, C)
    detA = M11 * M22 - M12 * M21
    u1[i] = (M22 * (v1[i] - L2H1) - M12 * (v2[i] - L2H2)) / detA
    u2[i] = (-M21 * (v1[i] - L2H1) + M11 * (v2[i] - L2H2)) / detA
    # updating state variable
    dx_temp = nl(x[0, i], x[1, i], x[2, i], x[3, i], p[i], q[i], w, L, C, u1[i], u2[i])
    dx[2, i] = dx_temp[2]
    dx[3, i] = dx_temp[3]
    for j in range(4):
        x[j, i + 1] = x[j, i] + dx_temp[j] * t_delta
    y[0, i + 1] = x[0, i + 1]
    y[1, i + 1] = x[1, i + 1]
    dz_temp = ss(dx[0, i], dx[1, i], v1[i], v2[i])
    for j in range(4):
        z[j, i + 1] = z[j, i] + dz_temp[j] * t_delta
    dz_temp[0] = dz_temp[0] + 0.0559 * (x[0, i + 1] - z[0, i + 1])
    dz_temp[1] = dz_temp[1] - 0.5000 * (x[0, i + 1] - z[0, i + 1])
    dz_temp[2] = dz_temp[2] + 0.0559 * (x[1, i + 1] - z[2, i + 1])
    dz_temp[3] = dz_temp[3] - 0.5000 * (x[1, i + 1] - z[2, i + 1])
    for j in range(4):
        z[j, i + 1] = z[j, i] + dz_temp[j] * t_delta
    # updating value
    if i == 1000:
        x[0, i + 1] = 225.2
        z[0, i + 1] = 225.2
    if i == 2000:
        x[0, i + 1] = 224
        z[0, i + 1] = 224
    if i == 3000:
        x[0, i + 1] = 225
        z[0, i + 1] = 225
    '''else:
        x[0, i + 1] = z[0, i + 1]
    x[1, i + 1] = z[2, i + 1]'''
    dx[0, i + 1] = z[1, i + 1]
    dx[1, i + 1] = z[3, i + 1]
    [iL[0, i + 1], iL[1, i + 1]] = il(x[0, i], x[1, i], x[2, i], x[3, i], p[i], q[i], w, C, L)
    p_out[i] = 2.1 * (u1[i] * x[2, i] + u2[i] * x[3, i])
    pint = pint + p_out[i] * 0.01
    q_out[i] = 2.1 * (u2[i] * x[2, i] - u1[i] * x[3, i])
    sOut[i] = sqrt(p_out[i] ** 2 + q_out[i] ** 2)
    PFOut[i] = p_out[i] / sOut[i]
    pqL[0, i] = 2.1 * (x[0, i] * iL[0, i] + x[1, i] * iL[1, i])
    pqL[1, i] = 2.1 * (x[1, i] * iL[0, i] - x[0, i] * iL[1, i])
    sL[i] = sqrt(pqL[0, i]**2 + pqL[1, i]**2)
    PFL[i] = pqL[0, i] / sL[i]

t1 = 1800
t2 = 4000
fig, (ax1, ax2) = plt.subplots(2)
ax1.plot(t[t1:t2], x[0, t1:t2], color='red', linestyle='-', linewidth='2', label='$V_{Ld}$')
ax1.plot(t[t1:t2], x[1, t1:t2], color='blue', linestyle='-', linewidth='2', label='$V_{Lq}$')
ax1.yaxis.set_minor_locator(MultipleLocator(50))
# ax1.xaxis.set_minor_locator(MultipleLocator(0.001))
ax1.set(title='Keluaran Tegangan dan PF Inverter', ylabel='Tegangan (V)')
ax1.legend()
ax1.grid(True, linestyle='-')


ax2.plot(t[t1:t2], PFL[t1:t2], color='green', linestyle='-', linewidth='2', label='$PF$')
ax2.yaxis.set_minor_locator(MultipleLocator(0.5))
ax2.axhline(y=1, color="grey", linestyle="-")
# ax2.axhline(y=0, color="grey", linestyle="-")
ax2.set(xlabel='Waktu (detik)', ylabel='PF')
ax2.legend()
ax2.grid(True, linestyle='-')
plt.show()
