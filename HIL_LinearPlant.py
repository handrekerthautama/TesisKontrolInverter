# import csv
import matplotlib.pyplot as plt
from Functions import ss
import numpy as np
import serial
import time

t_init = 0.0
t_final = 4.00
t_delta = 0.01
t = np.arange(t_init, t_final, t_delta)
ser = serial.Serial('COM3', baudrate=9600, timeout=1)
z = np.zeros((4, len(t)))
tmp = [' ', ' ', ' ', ' ']
z[0, 1] = 0.247346
z[1, 1] = 48.64925
z[2, 1] = 0
z[3, 1] = 0


for idx in range(len(t)-1):
    for i in range(4):
        tmp[i] = str(round(1000 * z[i, idx], 2)) + ","
    send = tmp[0] + tmp[2] + tmp[1] + tmp[3]
    ser.write(send.encode())
    time.sleep(1.0)
    data = ser.readline().decode('ascii')
    if len(data) > 0:
        split = 0
        for j in range(len(data)):
            if data[j] == ',':
                split = j

        v1 = float(data[0:split])
        v2 = float(data[split + 1:len(data) + 1])

        dz_temp = ss(z[1, idx], z[3, idx], v1, v2)
        for k in range(4):
            z[k, idx + 1] = z[k, idx] + dz_temp[k] * t_delta

plt.plot(t, z[0, :])
plt.show()
