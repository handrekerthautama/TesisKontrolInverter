# import csv
import matplotlib.pyplot as plt
from Functions import nl, l2h1, l2h2, m11, m12, m21, m22, pf, qf, il, transform
from matplotlib.ticker import (MultipleLocator)
from numpy.linalg import inv
import numpy as np
import serial
import time

w = 2 * 3.142 * 50
L = 0.0008
C = 0.000075
t_init = 0.0
t_final = 2.00
t_delta = 0.01
t = np.arange(t_init, t_final, t_delta)
ser = serial.Serial('COM3', baudrate=9600, timeout=1)
x = np.zeros((4, len(t)))
dx = np.zeros((4, len(t)))
z = np.zeros((4, len(t)))
p = np.zeros(len(t))
q = np.zeros(len(t))
u1 = np.zeros(len(t))
u2 = np.zeros(len(t))
iL = np.zeros((2, len(t)))
tmp = [' ', ' ', ' ', ' ']
z[0, 1] = 0.247346
z[1, 1] = 48.64925
z[2, 1] = 0
z[3, 1] = 0
x[0, 1] = 0.247346
x[1, 1] = 0
x[2, 1] = 0.00364869
x[3, 1] = 0.01165741
dx[0, 1] = 48.64925313
dx[1, 1] = 0
dx[2, 1] = 0.234600327
dx[3, 1] = -2.34669971

for i in range(len(t)-1):
    if i > 100:
        x[2, i] = 0.35 + C * dx[0, i] - C * w * x[1, i]
        x[3, i] = 0.00 + C * dx[1, i] + C * w * x[0, 1]
    if i == 100:
        x[0, i] = 225.5
    tmp[0] = str(round(1000 * x[0, i], 2)) + ","
    tmp[1] = str(round(1000 * x[1, i], 2)) + ","
    # tmp[2] = str(round(1000 * z[1, idx], 2)) + ","
    # tmp[3] = str(round(1000 * z[3, idx], 2)) + ","
    send = tmp[0] + tmp[1]
    ser.write(send.encode())
    time.sleep(1.0)
    data = ser.readline().decode('ascii')
    if len(data) > 0:
        split = [0, 0, 0]
        idx = 0
        for j in range(len(data)):
            if data[j] == ',':
                split[idx] = j
                idx = idx + 1

        v1 = float(data[0:split[0]])
        v2 = float(data[split[0] + 1:split[1]])
        z2 = float(data[split[1] + 1:split[2]])
        z4 = float(data[split[2] + 1:len(data)+1])
        dx[0, i] = z2
        dx[1, i] = z4

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
        u1[i] = (M22 * (v1 - L2H1) - M12 * (v2 - L2H2)) / detA
        u2[i] = (-M21 * (v1 - L2H1) + M11 * (v2 - L2H2)) / detA
        # updating state variable
        dx_temp = nl(x[0, i], x[1, i], x[2, i], x[3, i], p[i], q[i], w, L, C, u1[i], u2[i])
        dx[2, i] = dx_temp[2]
        dx[3, i] = dx_temp[3]
        for j in range(4):
            x[j, i + 1] = x[j, i] + dx_temp[j] * t_delta
        [iL[0, i + 1], iL[1, i + 1]] = il(x[0, i], x[1, i], x[2, i], x[3, i], p[i], q[i], w, C, L)
        # dz_temp = ss(z[1, idx], z[3, idx], v1, v2)
        # for k in range(4):
        #    z[k, idx + 1] = z[k, idx] + dz_temp[k] * t_delta

        # print([z2, z4, z[1, idx+1], z[3, idx+1]])

        print(v1, u1[i], x[0, i+1])

fig, ax = plt.subplots()
t1 = 0
t2 = 200
ax.plot(t[t1:t2], x[0, t1:t2], color='red', linestyle='-', linewidth='2', label='$V_{Ld}$')
ax.plot(t[t1:t2], x[1, t1:t2], color='blue', linestyle='-', linewidth='2', label='$V_{Lq}$')
# ax.plot(t[t1:t2], iL[0, t1:t2], color='red', linestyle='-.', linewidth='2', label='$I_{Ld}$')
# ax.plot(t[t1:t2], iL[1, t1:t2], color='blue', linestyle='-.', linewidth='2', label='$I_{Lq}$')
# ax.plot(t[t1:t2], pqL[0, t1:t2], color='red', linestyle='-.', linewidth='2', label='$p_{L}$')
# ax.plot(t[t1:t2], pqL[1, t1:t2], color='blue', linestyle='-.', linewidth='2', label='$q_{L}$')
ax.axhline(y=230, color="grey", linestyle="--")
ax.grid(True, linestyle='-')
ax.tick_params(labelcolor='black', labelsize='medium', width=3)
# ax.xaxis.set_minor_locator(MultipleLocator(0.05))
ax.yaxis.set_minor_locator(MultipleLocator(50))
ax.set_xlabel('Waktu (detik)')
ax.set_ylabel('Tegangan (V)')
ax.set_title('Tegangan Keluaran Inverter [Koordinat dq]')
ax.legend()
plt.show()

ts_init = t_init
ts_final = t_final
ts_delta = 0.001

ts = np.arange(ts_init, ts_final, ts_delta)
xs = np.zeros((4, len(ts)))
v_abc = np.zeros((3, len(ts)))
i_abc = np.zeros((3, len(ts)))
vdq0 = np.zeros((3, len(ts)))
idq0 = np.zeros((3, len(ts)))
v_abc_tmp = np.zeros((3, 1))
i_abc_tmp = np.zeros((3, 1))

# Transformation from dq to abc
for i in range(len(t)):
    var = i * 10
    for j in range(10):
        xs[0, var + j] = x[0, i]
        xs[1, var + j] = x[1, i]
        xs[2, var + j] = iL[0, i]
        xs[3, var + j] = iL[1, i]

for k in range(len(ts)):
    vdq0[0, k] = xs[0, k]
    vdq0[1, k] = xs[1, k]
    vdq0[2, k] = 0.000
    idq0[0, k] = xs[2, k]
    idq0[1, k] = xs[3, k]
    idq0[2, k] = 0.000
    T = transform(w, ts[k])
    T_inv = inv(T)
    vdq0_tmp = np.array(([[vdq0[0, k]], [vdq0[1, k]], [vdq0[2, k]]]))
    idq0_tmp = np.array(([[idq0[0, k]], [idq0[1, k]], [idq0[2, k]]]))
    v_abc_tmp = np.dot(T_inv, vdq0_tmp)
    i_abc_tmp = np.dot(T_inv, idq0_tmp)
    v_abc[0, k] = v_abc_tmp[0, 0]
    v_abc[1, k] = v_abc_tmp[1, 0]
    v_abc[2, k] = v_abc_tmp[2, 0]
    i_abc[0, k] = i_abc_tmp[0, 0]
    i_abc[1, k] = i_abc_tmp[1, 0]
    i_abc[2, k] = i_abc_tmp[2, 0]

t1 = 900
t2 = 1200
fig, (ax1, ax2) = plt.subplots(2)
ax1.plot(ts[t1:t2], v_abc[0, t1:t2], color='red', linestyle='-', linewidth='2', label='$V_{La}$')
ax1.plot(ts[t1:t2], v_abc[1, t1:t2], color='blue', linestyle='-', linewidth='2', label='$V_{Lb}$')
ax1.plot(ts[t1:t2], v_abc[2, t1:t2], color='green', linestyle='-', linewidth='2', label='$V_{Lc}$')
ax1.yaxis.set_minor_locator(MultipleLocator(50))
# ax1.xaxis.set_minor_locator(MultipleLocator(0.001))
ax1.set(title='Keluaran Inverter', ylabel='Tegangan (V)')
ax1.legend()
ax1.grid(True, linestyle='-')


ax2.plot(ts[t1:t2], i_abc[0, t1:t2], color='red', linestyle='-', linewidth='2', label='$I_{La}$')
ax2.plot(ts[t1:t2], i_abc[1, t1:t2], color='blue', linestyle='-', linewidth='2', label='$I_{Lb}$')
ax2.plot(ts[t1:t2], i_abc[2, t1:t2], color='green', linestyle='-', linewidth='2', label='$I_{Lc}$')
ax2.yaxis.set_minor_locator(MultipleLocator(0.1))
# ax2.xaxis.set_minor_locator(MultipleLocator(0.001))
# ax2.axhline(y=2, color="grey", linestyle="-")
ax2.set(xlabel='Waktu (detik)', ylabel='Arus (A)')
ax2.legend()
ax2.grid(True, linestyle='-')

plt.show()
