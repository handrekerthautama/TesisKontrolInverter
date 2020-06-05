import numpy as np
# from numpy.linalg import inv
import matplotlib.pyplot as plt
from Functions import ss, nl, l2h1, l2h2, m11, m12, m21, m22, pf, qf  # transform
# import csv

w = 2 * 3.142 * 50
L = 0.0008
C = 0.000075
t_init = 0.0
t_final = 15.0
t_delta = 0.01
r1 = 230
r2 = 0

t = np.arange(t_init, t_final, t_delta)
x = np.zeros((4, len(t)))
z = np.zeros((4, len(t)))
dx = np.zeros((4, len(t)))
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
kp = 1
ki = 10
kd = 0
x[0, 0] = 0.24734592
x[1, 0] = 0.0
x[2, 0] = 0.00364869
x[3, 0] = 0.01165741
dx[0, 0] = 48.64925313
dx[1, 0] = 0.0
e1[0] = 269.44
e1sum[0] = 5.3888
dz_temp = [0, 48.64925313, 0, 0]

for i in range(len(t) - 1):
    # controller
    e1[i] = r1 - x[0, i]
    e2[i] = r2 - x[1, i]
    for k in range(i):
        e1sum[i] = e1sum[i] + e1[k] * t_delta
        e2sum[i] = e2sum[i] + e2[k] * t_delta
    v1[i] = kp * e1[i] + ki * e1sum[i] + kd * (e1[i] + e1[i-1]) / 0.01
    v2[i] = kp * e2[i] + ki * e2sum[i] + kd * (e2[i] + e2[i-1]) / 0.01
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
    dz_temp[0] = dz_temp[0] + 0.0838 * (x[0, i + 1] - z[0, i + 1])
    dz_temp[1] = dz_temp[1] - 0.5000 * (x[0, i + 1] - z[0, i + 1])
    dz_temp[2] = dz_temp[2] + 0.0838 * (x[1, i + 1] - z[2, i + 1])
    dz_temp[3] = dz_temp[3] - 0.5000 * (x[1, i + 1] - z[2, i + 1])
    for j in range(4):
        z[j, i + 1] = z[j, i] + dz_temp[j] * t_delta
    # updating value
    x[0, i + 1] = z[0, i + 1]
    x[1, i + 1] = z[2, i + 1]
    dx[0, i + 1] = z[1, i + 1]
    dx[1, i + 1] = z[3, i + 1]

fig, ax = plt.subplots()
t1 = 0
t2 = 1400
ax.plot(t[t1:t2], x[0, t1:t2], color='red', linestyle='-', linewidth='2', label='$V_{Ld}$')
ax.plot(t[t1:t2], x[1, t1:t2], color='blue', linestyle='-', linewidth='2', label='$V_{Lq}$')
# ax.plot(t[t1:t2], z[0, t1:t2], color='red', linestyle='-', marker='o', linewidth='2', label='$I_{Ld}$')
# ax.plot(t[t1:t2], z[2, t1:t2], color='blue', linestyle='-', marker='o', linewidth='2', label='$I_{Lq}$')
# ax.plot(t[t1:t2], pqL[0, t1:t2], color='red', linestyle='-.', linewidth='2', label='$p_{L}$')
# ax.plot(t[t1:t2], pqL[1, t1:t2], color='blue', linestyle='-.', linewidth='2', label='$q_{L}$')
ax.axhline(y=230, color="grey", linestyle="--")
ax.grid(True, linestyle='-')
ax.tick_params(labelcolor='black', labelsize='medium', width=3)
# ax.xaxis.set_minor_locator(MultipleLocator(0.05))
# ax.yaxis.set_minor_locator(MultipleLocator(10))
ax.set_xlabel('Waktu (detik)')
ax.set_ylabel('Tegangan (V)')
ax.set_title('Tegangan Keluaran Inverter [Koordinat dq]')
ax.legend()
plt.show()

'''ts_init = 0.0
ts_final = 5.0
ts_delta = 0.001

ts = np.arange(ts_init, ts_final, ts_delta)
xs = np.zeros((4, len(ts)))
v_abc = np.zeros((3, len(ts)))
vdq0 = np.zeros((3, len(ts)))
v_abc_tmp = np.zeros((3, 1))

# Transformation from dq to abc
for i in range(len(t)):
    var = i * 10
    for j in range(10):
        xs[0, var + j] = x[0, i]
        xs[1, var + j] = x[1, i]
        xs[2, var + j] = x[2, i]
        xs[3, var + j] = x[3, i]

for k in range(len(ts)):
    vdq0[0, k] = xs[0, k]
    vdq0[1, k] = xs[1, k]
    vdq0[2, k] = 0.000
    T = transform(w, ts[k])
    T_inv = inv(T)
    vdq0_tmp = np.array(([[vdq0[0, k]], [vdq0[1, k]], [vdq0[2, k]]]))
    v_abc_tmp = np.dot(T_inv, vdq0_tmp)
    v_abc[0, k] = v_abc_tmp[0, 0]
    v_abc[1, k] = v_abc_tmp[1, 0]
    v_abc[2, k] = v_abc_tmp[2, 0]

plt.plot(ts[0:4999], v_abc[0, 0:4999], color='red', linestyle='-', linewidth='2', label='Va')
plt.plot(ts[0:4999], v_abc[1, 0:4999], color='blue', linestyle='-', linewidth='2', label='Vb')
plt.plot(ts[0:4999], v_abc[2, 0:4999], color='green', linestyle='-', linewidth='2', label='Vc')
plt.xlabel('Waktu (detik)')
plt.ylabel('Tegangan (V)')
plt.title('Tegangan Keluaran Inverter [Koordinat abc]')
plt.legend()
plt.grid(True)
plt.show()'''

'''with open('data.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    for k in range(len(ts)):
        writer.writerow([v_abc[0, k]])'''