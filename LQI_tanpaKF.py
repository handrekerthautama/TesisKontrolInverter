import numpy as np
# from numpy.linalg import inv
import matplotlib.pyplot as plt
from Functions import ss, nl, l2h1, l2h2, m11, m12, m21, m22, pf, qf  # transform


w = 2 * 3.142 * 50
L = 0.0008
C = 0.000075
t_init = 0.0
t_final = 15.0
t_delta = 0.01
r1 = 269.44
r2 = 0

ki = 1
k1 = 2.4142
k2 = 2.4142
t = np.arange(t_init, t_final, t_delta)
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

for i in range(len(t) - 1):
    # controller
    e1[i] = r1 - x[0, i]
    e2[i] = r2 - x[1, i]
    e1sum[i+1] = e1sum[i] + e1[i] * t_delta
    e2sum[i+1] = e2sum[i] + e2[i] * t_delta
    v1[i] = ki * e1sum[i+1] - k1 * x[0, i] - k2 * dx[0, i]
    v2[i] = ki * e2sum[i+1] - k1 * x[1, i] - k2 * dx[1, i]
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

y1 = x[0, :]

ki = 3.1623
k1 = 4.6054
k2 = 3.1954
t = np.arange(t_init, t_final, t_delta)
x = np.zeros((4, len(t)))
z = np.zeros((4, len(t)))
dx = np.zeros((4, len(t)))
y2 = np.zeros(len(t))
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

for i in range(len(t) - 1):
    # controller
    e1[i] = r1 - x[0, i]
    e2[i] = r2 - x[1, i]
    e1sum[i+1] = e1sum[i] + e1[i] * t_delta
    e2sum[i+1] = e2sum[i] + e2[i] * t_delta
    v1[i] = ki * e1sum[i+1] - k1 * x[0, i] - k2 * dx[0, i]
    v2[i] = ki * e2sum[i+1] - k1 * x[1, i] - k2 * dx[1, i]
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

y2 = x[0, :]

ki = 1
k1 = 3.2352
k2 = 2.7332
t = np.arange(t_init, t_final, t_delta)
x = np.zeros((4, len(t)))
z = np.zeros((4, len(t)))
dx = np.zeros((4, len(t)))
y3 = np.zeros(len(t))
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

for i in range(len(t) - 1):
    # controller
    e1[i] = r1 - x[0, i]
    e2[i] = r2 - x[1, i]
    e1sum[i+1] = e1sum[i] + e1[i] * t_delta
    e2sum[i+1] = e2sum[i] + e2[i] * t_delta
    v1[i] = ki * e1sum[i+1] - k1 * x[0, i] - k2 * dx[0, i]
    v2[i] = ki * e2sum[i+1] - k1 * x[1, i] - k2 * dx[1, i]
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

y3 = x[0, :]

ki = 1
k1 = 2.7332
k2 = 3.2352
t = np.arange(t_init, t_final, t_delta)
x = np.zeros((4, len(t)))
z = np.zeros((4, len(t)))
dx = np.zeros((4, len(t)))
y4 = np.zeros(len(t))
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

for i in range(len(t) - 1):
    # controller
    e1[i] = r1 - x[0, i]
    e2[i] = r2 - x[1, i]
    e1sum[i+1] = e1sum[i] + e1[i] * t_delta
    e2sum[i+1] = e2sum[i] + e2[i] * t_delta
    v1[i] = ki * e1sum[i+1] - k1 * x[0, i] - k2 * dx[0, i]
    v2[i] = ki * e2sum[i+1] - k1 * x[1, i] - k2 * dx[1, i]
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

y4 = x[0, :]

plt.plot(t[:], y1, color='black', linestyle='-', linewidth='2', label='y1')
plt.plot(t[:], y2, color='red', linestyle='-', linewidth='2', label='y2')
plt.plot(t[:], y3, color='green', linestyle='-', linewidth='2', label='y3')
plt.plot(t[:], y4, color='blue', linestyle='-', linewidth='2', label='y4')
plt.xlabel('Waktu (detik)')
plt.ylabel('Tegangan (V)')
plt.title('Tegangan Keluaran Inverter [Koordinat dq]')
plt.legend()
plt.grid(True)
plt.show()
