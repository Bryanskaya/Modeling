import matplotlib.pyplot as plt

a1 = 0.0134
b1 = 1
c1 = 4.35e-4
m1 = 1
a2 = 2.049 #* 20
b2 = 0.563e-3
c2 = 0.528e5
m2 = 1
alpha_0 = 0.05
alpha_N = 0.01
l = 10
T0 = 300
R = 0.5
Fс = -15

tau = 0.01
step = 0.1 #/ 2
N = round(l / step) + 1
eps = 1e-3

cur_temp = [T0] * N
x = [0] * N
for i in range(N):
    x[i] = i * step

def k(t):
    return a1 * (b1 + c1 * t**m1)

def c(t):
    return a2 + b2 * t**m2 - c2 / (t**2)

def alpha(x):
    d = alpha_N * l / (alpha_N - alpha_0)
    c = -(alpha_N * alpha_0 * l) / (alpha_N - alpha_0)

    return c / (x - d)

def p(x):
    return 2 / R * alpha(x)

def f(x):
    return 2 * T0 / R * alpha(x)

def hee_minus(n):
    return (k(cur_temp[n]) + k(cur_temp[n - 1])) / 2

def hee_plus(n):
    return (k(cur_temp[n]) + k(cur_temp[n + 1])) / 2

def A(n):
    return hee_minus(n) * tau / step

def D(n):
    return hee_plus(n) * tau / step

def B(n):
    return A(n) + D(n) + c(cur_temp[n]) * step + p(x[n]) * step * tau

def F(n, t):
    return f(x[n]) * step * tau + c(cur_temp[n]) * t * step

def M0():
    c_12 = c((cur_temp[0] + cur_temp[1]) / 2)
    p_12 = p(x[0] + step/2)
    hee_12 = hee_plus(0)

    return step / 8 * c_12 - hee_12 * tau / step + tau * step / 8 * p_12

def K0():
    p0 = p(x[0])
    p_12 = p(x[0] + step / 2)
    c0 = c(cur_temp[0])
    c_12 = c((cur_temp[0] + cur_temp[1]) / 2)
    hee_12 = hee_plus(0)

    return step / 8 * (c_12 + 2 * c0) + hee_12 * tau / step + tau * step / 8 * (p_12 + 2 * p0)

def P0():
    c0 = c(cur_temp[0])
    c_12 = c((cur_temp[0] + cur_temp[1]) / 2)
    f0 = f(x[0])
    f_12 = f(x[0] + step/2)

    return step/8 * c_12 * (cur_temp[0] + cur_temp[1]) + \
           step/4 * c0 * cur_temp[0] + Fс*tau + tau*step/4 * (f_12 + f0)

def MN():
    c_N = c(cur_temp[N - 1])
    c_12 = c((cur_temp[N - 1] + cur_temp[N - 2]) / 2)
    p_N = p(l)
    p_12 = p(l - step/2)
    hee_12 = hee_minus(N - 1)

    return step / 8 * (2*c_N + c_12) + tau*step/8 * (2*p_N + p_12) + tau * (alpha_N + hee_12/step)

def KN():
    hee_12 = hee_minus(N - 1)
    c_12 = c((cur_temp[N-1] + cur_temp[N-2]) / 2)
    p_12 = p(l - step/2)

    return c_12*step/8 + p_12*tau*step/8 - tau*hee_12/step

def PN():
    c_N = c(cur_temp[N - 1])
    c_12 = c((cur_temp[N - 1] + cur_temp[N - 2]) / 2)
    f_N = f(l)
    f_12 = f(l - step/2)

    return step/4 * (c_N*cur_temp[N-1] + c_12/2 * (cur_temp[N-1] + cur_temp[N-2])) + step*tau/4 * (f_N + f_12) + tau*alpha_N*T0

def calc_iter(prev_temp):
    m0 = M0()
    k0 = K0()
    p0 = P0()
    eps = [0] * N
    eta = [0] * N
    eps[1] = -m0 / k0
    eta[1] = p0 / k0

    for i in range(1, N - 1):
        eps[i + 1] = D(i) / (B(i) - A(i) * eps[i])
        eta[i + 1] = (F(i, prev_temp[i]) + A(i) * eta[i]) / (B(i) - A(i) * eps[i])

    mN = MN()
    kN = KN()
    pN = PN()
    new_temp = [0] * N
    new_temp[-1] = (pN - kN * eta[-1]) / (mN + kN * eps[-1])
    for i in range(N - 1, 0, -1):
        new_temp[i - 1] = eps[i] * new_temp[i] + eta[i]

    return new_temp

def temperature_check(old_t, new_t):
    for i in range(N):
        if abs((new_t[i] - old_t[i]) / new_t[i]) > eps:
            return False
    return True

def draw(x_arr, t_arr, temp_arr, step):
    x = x_arr
    t = t_arr[::step]
    temp = temp_arr[::step]

    t.append(t_arr[-1])
    temp.append(temp_arr[-1])

    for i in range(len(t)):
        #if i % 2 == 0:
            plt.plot(x,
                     temp[i],
                     label="{:.2f} с".format(t[i]))
            plt.grid()

    plt.xlabel("x, см")
    plt.ylabel("T, K")
    plt.legend()
    plt.show()

def draw_2(x_arr, t_arr, temp_arr, x_t):
    x_t.sort()
    x_ind = 0
    for i in range(len(x_arr)):
        if x_arr[i] + step / 2 > x_t[x_ind]:
            plt.plot(t_arr,
                     [tmp[i] for tmp in temp_arr],
                     label="{:.2f} см".format(x_arr[i]))
            plt.grid()

            x_ind += 1
            if x_ind == len(x_t):
                break

    plt.xlabel("t, с")
    plt.ylabel("T, K")
    plt.legend()
    plt.show()

def main():
    global cur_temp

    t_arr = [0]
    temp_arr = [cur_temp.copy()]

    while True:
        prev_temp = cur_temp
        while True:
            new_temp = calc_iter(prev_temp)

            if temperature_check(cur_temp, new_temp):
                break
            cur_temp = new_temp

        t_arr.append(t_arr[-1] + tau)
        temp_arr.append(new_temp)

        if temperature_check(prev_temp, new_temp):
            break
        cur_temp = new_temp

    cur_temp = new_temp

    draw(x, t_arr, temp_arr, round(1/tau))
    #draw_2(x, t_arr, temp_arr, [0.0, 0.2, 0.4, 0.6, 0.8, 1.0])


if __name__ == '__main__':
    main()