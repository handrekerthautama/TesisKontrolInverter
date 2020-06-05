# import module
import numpy as np
import matplotlib.pyplot as plt
from Functions import ss, nl, l2h1, l2h2, m11, m12, m21, m22, pf, qf, il  # transform
from matplotlib.ticker import (MultipleLocator)
# from numpy.linalg import inv

# simulation parameters
w = 2 * 3.142 * 50
L = 0.0008
C = 0.000075
t_init = 0.0
t_final = 2.01
t_delta = 0.01
r = np.array([230, 0])

# control gain
ki = 2236.1
k1 = 365.3
k2 = 27

# initialize variable
t = np.arange(t_init, t_final, t_delta)
x = np.zeros((4, len(t)))
z = np.zeros((4, len(t)))
dx = np.zeros((4, len(t)))
y = np.zeros((2, len(t)))
v = np.zeros((2, len(t)))
u = np.zeros((2, len(t)))
iL = np.zeros((2, len(t)))
pqf = np.zeros((2, len(t)))
e = np.zeros((2, len(t)))
e_sum = np.zeros((2, len(t)))

# variable assignment
x[0, 0] = 0.24734592
x[1, 0] = 0.0
z[0, 0] = 0.24734592
z[2, 0] = 0.0
x[2, 0] = 0.00364869
x[3, 0] = 0.01165741
dx[0, 0] = 48.64925313
dx[1, 0] = 0.0
z[1, 0] = 48.64925313
z[3, 0] = 0.0
e[0, 0] = 269.44
e_sum[0, 0] = 5.3888
y[0, 0] = 0.24734592
y[1, 0] = 0.0
ise = 0
iae = 0
itae = 0

for i in range(len(t) - 1):
    # load
    if i > 100:
        x[2, i] = 13.0 + C * z[1, i] - C * w * z[1, i]
        x[3, i] = 0.00 + C * z[3, i] + C * w * z[2, 1]
    # calculate error
    e[0, i] = r[0] - z[0, i]
    e[1, i] = r[1] - z[2, i]
    ise = ise + ((e[0, i]) ** 2) * t_delta
    iae = iae + abs(e[0, i]) * t_delta
    itae = itae + (i * abs(e[0, i])) ** t_delta
    for k in range(i):
        e_sum[0, i] = e_sum[0, i] + e[0, k] * t_delta
        e_sum[1, i] = e_sum[1, i] + e[1, k] * t_delta
    # calculate v1 dan v2
    v[0, i] = ki * e_sum[0, i] - k1 * z[0, i] - k2 * z[1, i]
    v[1, i] = ki * e_sum[1, i] - k1 * z[2, i] - k2 * z[3, i]
    # calculate u1 & u2
    pqf[0, i] = pf(z[0, i], z[2, i], x[2, i], x[3, i], z[1, i], z[3, i], C)
    pqf[1, i] = qf(z[0, i], z[2, i], x[2, i], x[3, i], z[1, i], z[3, i], w, L, C)
    L2H1 = l2h1(z[0, i], z[2, i], x[2, i], x[3, i], z[1, i], z[3, i], pqf[0, i], pqf[1, i], w, L, C)
    L2H2 = l2h2(z[0, i], z[2, i], x[2, i], x[3, i], z[1, i], z[3, i], pqf[0, i], pqf[1, i], w, L, C)
    M11 = m11(z[0, i], z[2, i], x[2, i], w, L, C)
    M12 = m12(z[0, i], z[2, i], x[3, i], w, L, C)
    M21 = m21(z[0, i], z[2, i], x[2, i], w, L, C)
    M22 = m22(z[0, i], z[2, i], x[3, i], w, L, C)
    detA = M11 * M22 - M12 * M21
    u[0, i] = (M22 * (v[0, i] - L2H1) - M12 * (v[1, i] - L2H2)) / detA
    u[1, i] = (-M21 * (v[0, i] - L2H1) + M11 * (v[1, i] - L2H2)) / detA
    # inverter model
    dx_temp = nl(x[0, i], x[1, i], x[2, i], x[3, i], pqf[0, i], pqf[1, i], w, L, C, u[0, i], u[1, i])
    dx[0, i] = dx_temp[0]
    dx[1, i] = dx_temp[1]
    dx[2, i] = dx_temp[2]
    dx[3, i] = dx_temp[3]
    for j in range(4):
        x[j, i + 1] = x[j, i] + dx_temp[j] * t_delta
    # measurement
    y[0, i + 1] = x[0, i + 1]
    y[1, i + 1] = x[1, i + 1]
    [iL[0, i + 1], iL[1, i + 1]] = il(x[0, i], x[1, i], x[2, i], x[3, i], pqf[0, i], pqf[1, i], w, C, L)
    # estimate states
    dz_temp = ss(z[1, i], z[3, i], v[0, i], v[1, i])
    for j in range(4):
        z[j, i + 1] = z[j, i] + dz_temp[j] * t_delta
    dz_temp[0] = dz_temp[0] + 0.0559 * (y[0, i + 1] - z[0, i + 1])
    dz_temp[1] = dz_temp[1] - 0.5000 * (y[0, i + 1] - z[0, i + 1])
    dz_temp[2] = dz_temp[2] + 0.0559 * (y[1, i + 1] - z[2, i + 1])
    dz_temp[3] = dz_temp[3] - 0.5000 * (y[1, i + 1] - z[2, i + 1])
    for j in range(4):
        z[j, i + 1] = z[j, i] + dz_temp[j] * t_delta

print([ise, iae, itae])

# plot result
fig, ax = plt.subplots()
t1 = 0
t2 = 200
ax.plot(t[t1:t2], x[0, t1:t2], color='red', linestyle='-', linewidth='2', label='$V_{Ld}$')
ax.plot(t[t1:t2], x[1, t1:t2], color='blue', linestyle='-', linewidth='2', label='$V_{Lq}$')
ax.axhline(y=230, color="grey", linestyle="--")
ax.tick_params(labelcolor='black', labelsize='medium', width=3)
# ax.yaxis.set_minor_locator(MultipleLocator(10))
ax.set_xlabel('Waktu (detik)')
ax.set_ylabel('Tegangan (V)')
ax.set_title('Tegangan Keluaran Inverter [Koordinat dq]')
ax.grid(True, linestyle='-')
ax.legend()
plt.show()
