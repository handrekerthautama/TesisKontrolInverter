import matplotlib.pyplot as plt

x1 = [-1.0000, -0.7071, -0.7071]
y1 = [0.00000, 0.70710, -0.7071]
x2 = [-1.3716, -0.7713, -0.7713]
y2 = [0.00000, 1.10176, -1.10176]
x3 = [-0.4547, -1.1392, -1.1392]
y3 = [0.00000, 0.9493, -0.9493]
x4 = [-2.1991, -0.5180, -0.5180]
y4 = [0.00000, 0.43170, -0.4317]
x5 = [-0.7291, -0.4731, -0.4731]
y5 = [0.00000, 0.62420, -0.6242]

plt.scatter(x1, y1, label='$I$', marker="x", color="black")
plt.scatter(x2, y2, label='$e_I$', marker="x", color="red")
plt.scatter(x3, y3, label='$V_{Ld},V_{Lq}$', marker="x", color="green")
plt.scatter(x4, y4, label='$dV_{Ld},dV_{Lq}$', marker="x", color="blue")
plt.scatter(x5, y5, label='$V_{Id},V_{Iq}$', marker="x", color="gold")
plt.axhline(y=0, color="grey", linestyle="-")
plt.axvline(x=0, color="grey", linestyle="-")
plt.xlabel('sumbu riil')
plt.ylabel('sumbu imajiner')
plt.title('Peta Poles')
plt.legend()
plt.grid()
plt.show()
