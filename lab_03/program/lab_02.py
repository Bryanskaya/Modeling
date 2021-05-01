from math import *
import matplotlib.pyplot as plt

table_I = [0.5, 1, 5, 10, 50, 200, 400, 800, 1200]
table_T0 = [6730, 6790, 7150, 7270, 8010, 9185, 10010, 11140, 12010]
table_m = [0.5, 0.55, 1.7, 3, 11, 32, 40, 41, 39]
table_T = [4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000]
table_sigma = [0.031, 0.27, 2.05, 6.06, 12, 19.9, 29.6, 41.1, 54.1, 67.7, 81.5]

R = 0.35
le = 12
Lk = 187 * 10 ** (-6)
Ck = 268 * 10 ** (-6)
Rk = 0.25
Tw = 2000

T0 = 0
Rp = 0

res_I = []
res_U = []
res_t = []
res_Rp = []
res_IRp = []
res_T0 = []

Ck_arr = [i for i in range(200, 650, 50)]

res_t_task = []
res_I_task = []
res_Ck_lenp = []

temp1 = [i for i in range(200, 700, 50)]
temp2 = [0.000485, 0.000546, 0.000602, 0.000654, 0.000703, 0.000750, 0.0007946999999998214,\
         0.000837999999999798, 0.0008796999999997755, 0.0009203999999997616]
temp3 = [100, 150, 200, 250, 300, 350, 400, 450, 500, 550]
temp4 = [425.50, 511.80, 584.70, 649.30, 708.20, 762.50, 813.20, 861.00, \
         906.20, 949.20]
temp5 = [0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
temp6 = [556.20, 562.90, 580.2, 604.5, 636.3, 675.3, 723.9, 781.01, 842.6,\
         903.2, 959.7]


def find_sigma(z, I, Tw):
    t0 = interpolation(I, 2, table_I, table_T0)

    global T0
    T0 = t0

    m = interpolation(I, 2, table_I, table_m)
    T = T0 + (Tw - T0) * z ** m
    sigma = interpolation(T, 2, table_T, table_sigma)

    return sigma


def find_Rp(I, Tw):
    a, b = 0, 1
    n = 100
    dz = (b - a) / n
    intgr = 0
    z = 0

    for j in range(n):
        intgr += z * dz * find_sigma(z, I, Tw)
        z += dz

    return le / (2 * pi * R * R * intgr)


def f(I, U):
    global Rp
    Rp = find_Rp(I, Tw)
    return (U - (Rk + Rp) * I) / Lk


def phi(I):
    return -I / Ck


def find_I_U(I, U, h):
    k1 = h * f(I, U)
    p1 = h * phi(I)

    k2 = h * f(I + k1 / 2, U + p1 / 2)
    p2 = h * phi(I + k1 / 2)

    k3 = h * f(I + k2 / 2, U + p2 / 2)
    p3 = h * phi(I + k2 / 2)

    k4 = h * f(I + k3, U + p3)
    p4 = h * phi(I + k3)

    return I + 1 / 6 * (k1 + 2 * k2 + 2 * k3 + k4), \
           U + 1 / 6 * (p1 + 2 * p2 + 2 * p3 + p4)


def find_position(x, data_x, number):
    left = 0
    right = number - 1

    if data_x[0] < data_x[-1]:
        if x < data_x[left] or x > data_x[right]: return -1
        while left + 1 != right:
            pos_x = int((left + right)/2)
            if data_x[pos_x] <= x:
                left = pos_x
            else:
                right = pos_x
        return left
    else:
        if x < data_x[right] or x > data_x[left]: return -1
        while left + 1 != right:
            pos_x = int((left + right)/2)
            if data_x[pos_x] >= x:
                left = pos_x
            else:
                right = pos_x
        return left


def find_range(pos_x, degree, data_x, data_y, number):
    half = int(degree / 2)
    left = pos_x - half
    right = pos_x + (degree - half)
    if left < 0:
        right += -left
        left = 0
    elif right > number - 1:
        left -= right - (number - 1)
        right = number - 1

    return data_x[left: right + 1], data_y[left: right + 1]


def find_polynomial(degree, nodes_x, nodes_y):
    div_diff, old_arr = [0] * (degree + 1), nodes_y
    new_arr = [0] * degree

    div_diff[0] = old_arr[0]
    for i in range(degree):
        for j in range(degree - i):
            new_arr[j] = (old_arr[j] - old_arr[j+1]) / \
                         (nodes_x[j] - nodes_x[j+i+1])
        div_diff[i+1] = new_arr[0]
        old_arr = new_arr

    return lambda x: nwtn_polynom(x, degree, div_diff, nodes_x)


def nwtn_polynom(x, degree, div_diff, nodes_x):
    y = div_diff[0]
    x_pl = 1
    for i in range(1, degree + 1):
        x_pl *= (x - nodes_x[i-1])
        y += x_pl * div_diff[i]

    return y


def interpolation(x, degree, data_x, data_y):
    pos = find_position(x, data_x, len(data_x))
    nodes_x, nodes_y = find_range(pos, degree, data_x, data_y, len(data_x))
    f = find_polynomial(degree, nodes_x, nodes_y)

    return f(x)


def show_graph():
    '''plt.subplot(2, 2, 1)
    plt.xlabel('Ck')
    plt.ylabel('tимп')
    plt.plot(temp1, temp2)

    plt.subplot(2, 2, 2)
    plt.xlabel('Lk')
    plt.ylabel('tимп')
    plt.plot(temp3, temp4)

    plt.subplot(2, 2, 3)
    plt.xlabel('Rk')
    plt.ylabel('tимп')
    plt.plot(temp5, temp6)'''

    plt.subplot(3, 2, 1)
    plt.xlabel('t')
    plt.ylabel('I')
    plt.plot(res_t, res_I)

    plt.subplot(3, 2, 2)
    plt.xlabel('t')
    plt.ylabel('U')
    plt.plot(res_t, res_U)

    plt.subplot(3, 2, 3)
    plt.xlabel('t')
    plt.ylabel('Rp')
    plt.plot(res_t, res_Rp)

    plt.subplot(3, 2, 4)
    plt.xlabel('t')
    plt.ylabel('I * Rp')
    plt.plot(res_t, res_IRp)

    plt.subplot(3, 2, 5)
    plt.xlabel('t')
    plt.ylabel('T0')
    plt.plot(res_t, res_T0)

    plt.show()


def find_length_p():
    max_I = max(res_I)
    i, t = 0, 0

    while res_I[i] < 0.35 * max_I:
        i += 1

    t_start = res_t[i]

    while i < len(res_I) and res_I[i] >= 0.35 * max_I:
        i += 1

    t_end = res_t[i]

    return t_end - t_start


def main():
    Uc = 1400
    I = 0
    h = 10 ** (-7)
    t = 0

    while t < 20 * 10 ** (-6): #-3
        I, Uc = find_I_U(I, Uc, h)
        t += h

        res_I.append(I)
        res_U.append(Uc)
        res_t.append(t)
        res_T0.append(T0)
        res_Rp.append(Rp)

    for i in range(len(res_Rp)):
        res_IRp.append(res_I[i] * res_Rp[i])

    show_graph()

    #print(Rk, find_length_p())


if __name__ == '__main__':
    main()