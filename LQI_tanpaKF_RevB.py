import numpy as np
from matplotlib import pyplot as plt
from Functions import ss, nl, l2h1, l2h2, m11, m12, m21, m22, pf, qf  # transform

w = 2 * 3.142 * 50
L = 0.0008
C = 0.000075
t_init = 0.0
t_final = 15.0
t_delta = 0.01
r = np.array([230, 0])
t = np.arange(t_init, t_final, t_delta)
ki = [3162.3, 3162.3, 1414.2, 2236.1, 2236.1]
k1 = [431.40, 446.10, 257.20, 351.60, 365.30]
k2 = [29.400, 29.900, 22.700, 26.500, 27.000]
y = np.zeros((len(ki), len(t)))

for k in range(len(ki)):
    x = np.zeros((4, len(t)))
    z = np.zeros((4, len(t)))
    dx = np.zeros((4, len(t)))
    v = np.zeros((2, len(t)))
    u = np.zeros((2, len(t)))
    iL = np.zeros((2, len(t)))
    pqf = np.zeros((2, len(t)))
    e = np.zeros((2, len(t)))
    e_sum = np.zeros((2, len(t)))
    e1 = np.zeros(len(t))
    e2 = np.zeros(len(t))
    e1sum = np.zeros(len(t))
    e2sum = np.zeros(len(t))

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

    for i in range(len(t) - 1):
        # calculate v1 dan v2
        e1[i] = r[0] - z[0, i]
        e2[i] = r[1] - z[2, i]
        e1sum[i + 1] = e1sum[i] + e1[i] * t_delta
        e2sum[i + 1] = e2sum[i] + e2[i] * t_delta
        v[0, i] = ki[k] * e_sum[0, i] - k1[k] * z[0, i] - k2[k] * z[1, i]
        v[1, i] = ki[k] * e_sum[1, i] - k1[k] * z[2, i] - k2[k] * z[3, i]
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

    y[k, :] = x[0, :]

# plot
t1 = 0
t2 = 200

plt.plot(t[t1:t2], y[0, t1:t2], color='black', linestyle='-', linewidth='2', label='$y_1$')
plt.plot(t[t1:t2], y[1, t1:t2], color='red', linestyle='-', linewidth='2', label='$y_2$')
plt.plot(t[t1:t2], y[2, t1:t2], color='green', linestyle='-', linewidth='2', label='$y_3$')
plt.plot(t[t1:t2], y[3, t1:t2], color='blue', linestyle='-', linewidth='2', label='$y_4$')
plt.plot(t[t1:t2], y[4, t1:t2], color='gold', linestyle='-', linewidth='2', label='$y_5$')
# plt.axhline(y=274.8288, color="grey", linestyle="--")
# plt.axhline(y=269.4400, color="grey", linestyle="--")
# plt.axhline(y=264.0512, color="grey", linestyle="--")
plt.xlabel('Waktu (detik)')
plt.ylabel('Tegangan (V)')
plt.title('Tegangan Keluaran Inverter [Koordinat dq]')
plt.legend()
plt.grid(True)
plt.show()
