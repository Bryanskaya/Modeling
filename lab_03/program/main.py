import lab_02 as interpolate
import matplotlib.pyplot as plt

np = 1.4
l = 0.2
T0 = 300
sigma = 5.668e-12
F0 = 100
alpha = 0.05 * 10

step = 2e-4
N = round(l / step) + 1
max_iter = 20
e1, e2 = 1e-3, 1e-3

t = [0] * N
t_balanced = [0] * N

k_relax = 0.3

table_k = [2.0e-2, 5.0e-2, 7.8e-2, 1.0e-1, 1.3e-1, 2.0e-1]
table_kt = [293, 1278, 1528, 1677, 2000, 2400]

table_lambda = [1.36e-2, 1.63e-2, 1.81e-2, 1.98e-2, 2.50e-2, 2.74e-2]
table_lambdat = [300, 500, 800, 1100, 2000, 2400]


def draw_graphs():
    x = [i * step for i in range(0, N)]

    plt.plot(x, t)
    plt.xlabel("x, см")
    plt.ylabel("T, K")
    plt.grid()
    plt.show()



def k(t):
    return interpolate.interpolation(t, 2, table_kt, table_k)


def lmbd(t):
    return interpolate.interpolation(t, 2, table_lambdat, table_lambda)


def hee_minus(n):
    return (lmbd(t_balanced[n]) + lmbd(t_balanced[n - 1])) / 2


def hee_plus(n):
    return (lmbd(t_balanced[n]) + lmbd(t_balanced[n + 1])) / 2


def f(n):
    return -4 * k(t_balanced[n]) * np**2 * sigma * (t_balanced[n]**4 - T0**4)


def A(n):
    return hee_minus(n) / step
def C(n):
    return hee_plus(n) / step
def B(n):
    return A(n) + C(n)
def D(n):
    return f(n) * step


def K0():
    return hee_plus(0)
def M0():
    return -hee_plus(0)
def P0():
    return step * F0 + step**2 / 4 * ((f(0) + f(1)) / 2 + f(0))
def KN():
    return -lmbd(N - 1)
def MN():
    return alpha * step + lmbd(N - 1)
def PN():
    return alpha * T0 * step


def forward_move():
    eps = [0] * N
    eta = [0] * N
    #TODO index = 1 => index = 0
    eps[1] = - M0() / K0()
    eta[1] = P0() / K0()

    for i in range(1, N - 1):
        eps[i + 1] = C(i) / (B(i) - A(i) * eps[i])
        eta[i + 1] = (A(i) * eta[i] + D(i)) / (B(i) - A(i) * eps[i])

    return eps, eta


def back_move(eps, eta):
    t2 = [0] * N
    t2[N - 1] = (PN() - MN() * eta[N - 1]) / (KN() + MN() * eps[N - 1])
    for i in range(N - 1, 0, -1):
        t2[i - 1] = eps[i] * t2[i] + eta[i]

    return t2


def init_t_values(begin, end):
    global t, t_balanced
    temp_step = (end - begin) / (N - 1)

    for i in range(N):
        t[i] = begin + temp_step * i
        t_balanced[i] = t[i]


def temperature_check(t_old, t_new):
    for i in range(N):
        if abs((t_new[i] - t_old[i]) / t_new[i]) > e1:
            return False
    return True


def energy_check():
    f1 = F0 - alpha * (t_balanced[N - 1] - T0)
    s = 0
    for i in range(N):
        #s += 4 * k(t_balanced[i]) * np ** 2 * sigma * (t_new[i] ** 4 - T0 ** 4)
        s += -f(i)
    f2 = s * l / N
    print("&${:.4f}$ &${:.4f}$ \\\\".format(f1, f2))
    return abs((f1 - f2) / f1) <= e2


def do_balance():
    global t_balanced
    for i in range(N):
        t_balanced[i] += k_relax * (t[i] - t_balanced[i])


def main():
    global t
    init_t_values(1200, 300)

    for j in range(max_iter):
        eps, eta = forward_move()
        t2 = back_move(eps, eta)

        t_old, t = t, t2
        if temperature_check(t_old, t) and energy_check():
            break

        do_balance()

    draw_graphs()


if __name__ == '__main__':
    main()