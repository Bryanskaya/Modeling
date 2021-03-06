from math import *
import matplotlib.pyplot as plt

table_I = [0.5, 1, 5, 10, 50, 200, 400, 800, 1200]
table_T0 = [6730, 6790, 7150, 7270, 8010, 9185, 10010, 11140, 12010]
table_m = [0.5, 0.55, 1.7, 3, 11, 32, 40, 41, 39]
table_T = [4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000]
table_sigma = [0.031, 0.27, 2.05, 6.06, 12, 19.9, 29.6, 41.1, 54.1, 67.7, 81.5]

res_I = []
res_U = []


def find_sigma(z, I, Tw):
    T0 = interpolation(I, 4, table_I, table_T0)
    m = interpolation(I, 4, table_I, table_m)
    T = T0 + (Tw - T0) * z ** int(m)
    sigma = interpolation(T, 4, table_T, table_sigma)

    return sigma


def find_Rp(le, R, I, Tw):
    a, b = 0, 1
    n = 100
    dz = (b - a) / n
    intgr = 0
    z = 0

    for j in range(n):
        intgr += z * dz * find_sigma(z, I, Tw)
        z += dz

    return le / (2 * pi * R * R * intgr)


def f(I, U, Rk, Lk, le, R, Tw):
    return (U - (Rk + find_Rp(le, R, I, Tw)) * I * 0) / Lk


def phi(I, Ck):
    return -I / Ck


def find_I(I, U, Rk, Lk, h, Ck, le, R, Tw):
    k1 = h * f(I, U, Rk, Lk, le, R, Tw)
    p1 = h * phi(I, Ck)

    k2 = h * f(I + k1 / 2, U + p1 / 2, Rk, Lk, le, R, Tw)
    p2 = h * phi(I + k1 / 2, Ck)

    k3 = h * f(I + k2 / 2, U + p2 / 2, Rk, Lk, le, R, Tw)
    p3 = h * phi(I + k2 / 2, Ck)

    k4 = h * f(I + k3, U + p3, Rk, Lk, le, R, Tw)
    p4 = h * phi(I + k3, Ck)

    return I + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def find_U(I, U, Rk, Lk, h, Ck, le, R, Tw):
    k1 = h * f(I, U, Rk, Lk, le, R, Tw)
    p1 = h * phi(I, Ck)

    k2 = h * f(I + k1 / 2, U + p1 / 2, Rk, Lk, le, R, Tw)
    p2 = h * phi(I + k1 / 2, Ck)

    k3 = h * f(I + k2 / 2, U + p2 / 2, Rk, Lk, le, R, Tw)
    p3 = h * phi(I + k2 / 2, Ck)

    k4 = h * f(I + k3, U + p3, Rk, Lk, le, R, Tw)
    p4 = h * phi(I + k3, Ck)

    return U + h / 6 * (p1 + 2 * p2 + 2 * p3 + p4)


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


def find_polynomial(x, y0, degree, nodes_x, result):
    Pn = y0
    prev = 1
    for i in range(degree):
        prev *= (x - nodes_x[i])
        Pn += prev * result[i]

    return Pn


def find_diff(degree, nodes_x, nodes_y):
    result = []

    for i in range(degree):
        for j in range(degree - i):
            nodes_y[j] = (nodes_y[j] - nodes_y[j + 1]) / \
                         (nodes_x[j] - nodes_x[j + i + 1])
        result.append(nodes_y[0])

    return result


def interpolation(x, degree, data_x, data_y):
    pos = find_position(x, data_x, len(data_x))
    nodes_x, nodes_y = find_range(pos, degree, data_x, data_y, len(data_x))
    res = find_diff(degree, nodes_x, nodes_y)
    return round(find_polynomial(x, nodes_y[0], degree, nodes_x, res), 3)


def show_graph():
    plt.plot(res_I)
    #plt.plot(res_U)

    plt.show()

def main():
    R = 0.35
    le = 12
    Lk = 187 * 10 ** (-6)
    Ck = 268 * 10 ** (-6)
    Rk = 0.25
    Uc = 1400
    I = 0  # 0..3
    Tw = 2000
    h = 10 ** (-6)

    for i in range(1200):
        I = find_I(I, Uc, Rk, Lk, h, Ck, le, R, Tw)
        Uc = find_U(I, Uc, Rk, Lk, h, Ck, le, R, Tw)

        res_I.append(I)
        res_U.append(Uc)

    show_graph()


if __name__ == '__main__':
    main()