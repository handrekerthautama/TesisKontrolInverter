import serial
import time

ser = serial.Serial('COM3', baudrate=9600, timeout=1)

while 1:
    i = "22,"
    j = "33,"
    ser.write(i.encode())
    ser.write(j.encode())
    time.sleep(0.5)
    data = ser.readline().decode('ascii')
    if len(data) > 0:
        print(data[0:2])
        print(data[2:4])
