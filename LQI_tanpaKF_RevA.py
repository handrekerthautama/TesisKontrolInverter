import numpy as np
from matplotlib import pyplot as plt
from Functions import ss, nl, l2h1, l2h2, m11, m12, m21, m22, pf, qf  # transform
from matplotlib.ticker import (MultipleLocator)

w = 2 * 3.142 * 50
L = 0.0008
C = 0.000075
t_init = 0.0
t_final = 20.0
t_delta = 0.01
r1 = 230
r2 = 0
t = np.arange(t_init, t_final, t_delta)
ki = [3.1623, 3.1623, 3.1623]
k1 = [4.6054, 5.1187, 5.0587]
k2 = [3.1954, 3.3522, 3.8881]
y = np.zeros((len(ki), len(t)))
param = np.zeros((3, len(ki)))

for k in range(len(ki)):
    x = np.zeros((4, len(t)))
    z = np.zeros((4, len(t)))
    dx = np.zeros((4, len(t)))
    y1 = np.zeros(len(t))
    v1 = np.zeros(len(t))
    v2 = np.zeros(len(t))
    u1 = np.zeros(len(t))
    u2 = np.zeros(len(t))
    p = np.zeros(len(t))
    q = np.zeros(len(t))
    e1 = np.zeros(len(t))
    e2 = np.zeros(len(t))
    e1sum = np.zeros(len(t))
    e2sum = np.zeros(len(t))
    x[0, 0] = 0.24734592
    x[1, 0] = 0.0
    x[2, 0] = 0.00364869
    x[3, 0] = 0.01165741
    dx[0, 0] = 48.64925313
    dx[1, 0] = 0.0
    e1[0] = 269.44
    e1sum[0] = 5.3888
    ise_tmp = 0
    iae_tmp = 0
    itae_tmp = 0

    for i in range(len(t) - 1):
        # controller
        e1[i] = r1 - x[0, i]
        e2[i] = r2 - x[1, i]
        ise_tmp = ise_tmp + ((e1[i]) ** 2) * t_delta
        iae_tmp = iae_tmp + abs(e1[i]) * t_delta
        itae_tmp = itae_tmp + (i * abs(e1[i])) ** t_delta
        e1sum[i + 1] = e1sum[i] + e1[i] * t_delta
        e2sum[i + 1] = e2sum[i] + e2[i] * t_delta
        v1[i] = ki[k] * e1sum[i + 1] - k1[k] * x[0, i] - k2[k] * dx[0, i]
        v2[i] = ki[k] * e2sum[i + 1] - k1[k] * x[1, i] - k2[k] * dx[1, i]
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
        dz_temp = ss(dx[0, i], dx[1, i], v1[i], v2[i])
        for j in range(4):
            z[j, i + 1] = z[j, i] + dz_temp[j] * t_delta
        # updating value
        x[0, i + 1] = z[0, i + 1]
        x[1, i + 1] = z[2, i + 1]
        dx[0, i + 1] = z[1, i + 1]
        dx[1, i + 1] = z[3, i + 1]

    param[0, k] = ise_tmp
    param[1, k] = iae_tmp
    param[2, k] = itae_tmp

    y[k, :] = x[0, :]

fig, ax = plt.subplots()
t1 = 0
t2 = 1400
ax.plot(t[t1:t2], y[0, t1:t2], color='red', linestyle='-', linewidth='2', label='$y_1$')
ax.plot(t[t1:t2], y[1, t1:t2], color='green', linestyle='-', linewidth='2', label='$y_2$')
ax.plot(t[t1:t2], y[2, t1:t2], color='blue', linestyle='-', linewidth='2', label='$y_3$')
# ax.plot(t[t1:t2], y[3, t1:t2], color='blue', linestyle='-', linewidth='2', label='$dV_{Ld}, dV_{Lq}$')
# ax.plot(t[t1:t2], y[4, t1:t2], color='gold', linestyle='-', linewidth='2', label='$V_{Id}, V_{Iq}$')
ax.axhline(y=230, color="grey", linestyle="--")
ax.grid(True, linestyle='-')
ax.tick_params(labelcolor='black', labelsize='medium', width=3)
ax.yaxis.set_minor_locator(MultipleLocator(10))
ax.xaxis.set_minor_locator(MultipleLocator(0.5))
ax.set_xlabel('Waktu (detik)')
ax.set_ylabel('Tegangan (V)')
ax.set_title('Tegangan Keluaran Inverter [Koordinat dq]')
ax.legend()
plt.show()

print(param)
