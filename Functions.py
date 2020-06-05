from math import pi, sin, cos, sqrt
import numpy as np


def f1(x1, x2, p, cf):
    val = p * x1 / (cf * (x1 ** 2 + x2 ** 2))
    return val


def f2(x1, x2, x3, x4, q, w, lf, cf):
    val = (q - w * lf * (x3 ** 2 + x4 ** 2)) * x2 / (cf * (x1 ** 2 + x2 ** 2))
    return val


def dx1dt(x1, x2, x3, x4, p, q, w, lf, cf):
    val = x3 / cf - f1(x1, x2, p, cf) - f2(x1, x2, x3, x4, q, w, lf, cf)
    return val


def dx2dt(x1, x2, x3, x4, p, q, w, lf, cf):
    val = x4 / cf - f1(x2, x1, p, cf) + f2(x2, x1, x3, x4, q, w, lf, cf)
    return val


def dx3dt(x1, x4, w, lf):
    val = w * x4 - x1 / lf
    return val


def dx4dt(x2, x3, w, lf):
    val = -w * x3 - x2 / lf
    return val


def nl(x1, x2, x3, x4, p, q, w, lf, cf, u1, u2):
    dx1 = dx1dt(x1, x2, x3, x4, p, q, w, lf, cf)
    dx2 = dx2dt(x1, x2, x3, x4, p, q, w, lf, cf)
    dx3 = dx3dt(x1, x4, w, lf) + u1 / lf
    dx4 = dx4dt(x2, x3, w, lf) + u2 / lf
    return [dx1, dx2, dx3, dx4]


def ff1(x1, x2, x3, x4, p, q, w, lf, cf):
    val = (p * (x1 ** 2 - x2 ** 2) + 2 * x1 * x2 * (q - w * lf * (x3 ** 2 + x4 ** 2))) / (cf * (x1 ** 2 + x2 ** 2) ** 2)
    return val


def ff2(x1, x2, x3, x4, p, q, w, lf, cf):
    val = (2 * p * x1 * x2 - (q - w * lf * (x3 ** 2 + x4 ** 2)) * (x1 ** 2 - x2 ** 2)) / (cf * (x1 ** 2 + x2 ** 2) ** 2)
    return val


def ff3(x1, x2, x3, w, lf, cf):
    val = (2 * w * lf * x2 * x3) / (cf * (x1 ** 2 + x2 ** 2))
    return val


def ff4(x1, x2, x4, w, lf, cf):
    val = (2 * w * lf * x2 * x4) / (cf * (x1 ** 2 + x2 ** 2))
    return val


def l2h1(x1, x2, x3, x4, dx1, dx2, p, q, w, lf, cf):
    val = dx1 * ff1(x1, x2, x3, x4, p, q, w, lf, cf) + dx2 * ff2(x1, x2, x3, x4, p, q, w, lf, cf) + \
          (1 / cf + ff3(x1, x2, x3, w, lf, cf)) * dx3dt(x1, x4, w, lf) + \
          ff4(x1, x2, x4, w, lf, cf) * dx4dt(x2, x3, w, lf)
    return val


def l2h2(x1, x2, x3, x4, dx1, dx2, p, q, w, lf, cf):
    val = dx1 * ff2(x1, x2, x3, x4, p, q, w, lf, cf) - dx2 * ff1(x1, x2, x3, x4, p, q, w, lf, cf) - \
          ff3(x2, x1, x3, w, lf, cf) * dx3dt(x1, x4, w, lf) + \
          (1 / cf - ff4(x2, x1, x4, w, lf, cf)) * dx4dt(x2, x3, w, lf)
    return val


def m11(x1, x2, x3, w, lf, cf):
    val = (1 / cf + ff3(x1, x2, x3, w, lf, cf)) / lf
    return val


def m12(x1, x2, x4, w, lf, cf):
    val = ff4(x1, x2, x4, w, lf, cf) / lf
    return val


def m21(x1, x2, x3, w, lf, cf):
    val = -1 * ff3(x2, x1, x3, w, lf, cf) / lf
    return val


def m22(x1, x2, x4, w, lf, cf):
    val = (1 / cf - ff4(x2, x1, x4, w, lf, cf)) / lf
    return val


def ss(z2, z4, v1, v2):
    dz1dt = 1 * z2
    dz2dt = 1 * v1
    dz3dt = 1 * z4
    dz4dt = 1 * v2
    return [dz1dt, dz2dt, dz3dt, dz4dt]


def kf(z1, z2, z3, z4, v1, v2, y1, y2):
    dz1dt = -0.05092 * z1 + 1 * z2 + 0.05092 * y1
    dz2dt = -2.205 * z1 - 2.115 * z2 + 1 * v1 - 0.0307 * y1
    dz3dt = -0.05092 * z3 + 1 * z4 + 0.05092 * y2
    dz4dt = -2.205 * z3 - 2.115 * z4 + 1 * v2 - 0.0307 * y2
    return [dz1dt, dz2dt, dz3dt, dz4dt]


def pf(x1, x2, x3, x4, dx1, dx2, cf):
    val = x1 * x3 + x2 * x4 - cf * (x1 * dx1 + x2 * dx2)
    return val


def qf(x1, x2, x3, x4, dx1, dx2, w, lf, cf):
    val = x2 * x3 - x1 * x4 + cf * (x1 * dx2 - dx1 * x2) + w * lf * (x3 ** 2 + x4 ** 2)
    return val


def transform(w, t):
    a1 = np.array([sin(w * t), sin(w * t - 2 * pi / 3), sin(w * t + 2 * pi / 3)])
    a2 = np.array([cos(w * t), cos(w * t - 2 * pi / 3), cos(w * t + 2 * pi / 3)])
    a3 = np.array([1/2, 1/2, 1/2])
    a = np.array(([2 / 3 * a1, 2 / 3 * a2, 2 / 3 * a3]))
    return a


def il(x1, x2, x3, x4, p, q, w, cf, lf):
    sq_v = x1 ** 2 + x2 ** 2
    sq_i = x3 ** 2 + x4 ** 2
    ild = (p * x1 + q * x2) / sq_v + w * cf * x2 - (w * lf * x2 * sq_i) / sq_v
    ilq = (p * x2 - q * x1) / sq_v - w * cf * x1 + (w * lf * x1 * sq_i) / sq_v + 5.419939231630957
    return [ild, ilq]
