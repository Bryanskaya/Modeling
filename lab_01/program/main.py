# Метод Пикара, метод Эйлера, метод Рунге-Кутта
from math import sqrt


def f(x, y):
    return x ** 2 + y ** 2


# Метод Пикара
def picard_1(x_args):
    res = []

    for x in x_args:
        res.append(x ** 3 / 3)

    return res


def picard_2(x_args):
    res = []

    for x in x_args:
        res.append(x ** 3 / 3 + x ** 7 / 63)

    return res


def picard_3(x_args):
    res = []

    for x in x_args:
        res.append(x ** 3 / 3 + x ** 7 / 63 + x ** 15 / 59535 + 2 * x ** 11 / 2079)

    return res


def picard_4(x_args):
    res = []

    for x in x_args:
        res.append(x ** 3 / 3 + x ** 7 / 63 + x ** 15 / 59535 + 2 * x ** 11 / 2079 + \
                   x ** 31 / 109876902975 + 4 * x ** 23 / 99411543 + \
                   4 * x ** 27 / 3341878155 + 2 * x ** 23 / 86266215 + \
                   2 * x ** 19 / 3393495 + 4 * x ** 19 / 2488563 + 4 * x ** 15 / 93555)

    return res


def runge_kutta(x, y, h, num):
    alpha = 0.5
    res = []

    temp = h / (2 * alpha)

    for i in range(num):
        res.append(y)

        k1 = f(x, y)
        k2 = f(x + temp, y + temp * k1)

        y += h * ((1 - alpha) * k1 + alpha * k2)
        x += h  # or temp

    return res


def euler_explicit(x, y, h, num):
    res = []

    for i in range(num):
        res.append(y)

        try:
            y += h * f(x, y)
            x += h
        except OverFlowError:
            for k in range(i, num):
                res.append('---')
            break

    return res


def euler_implicit(x, y, h, num):
    res = []

    for i in range(num):
        res.append(y)

        x += h
        D = 1 - 4 * h * y - 4 * h ** 2 * x ** 2
        if D < 0:
            for k in range(i, n):
                res.append('---')
            break
        y = (1 - sqrt(D)) / (2 * h)

    return res


def count_x_args(x, x_max, h):
    x_args = []

    while x <= x_max:
        x_args.append(x)
        x += h

    return x_args
    

def print_head():
    print(' ' * 4 + 'x' + ' ' * 4 + '|' + \
          ' ' * 17 + 'Метод Пикара' + ' ' * 18 + '|' + \
          ' ' * 6 + 'Метод Эйлера' + ' ' * 5 + '|' + \
          ' ' * 3 + 'Метод Рунге-Кутта' + ' ' * 3 +
          '\n' + \
          ' ' * 9 + '|' + \
          ' ' * 5 + '1' + ' ' * 5 + '|' + \
          ' ' * 5 + '2' + ' ' * 5 + '|' + \
          ' ' * 5 + '3' + ' ' * 5 + '|' + \
          ' ' * 5 + '4' + ' ' * 5 + '|' + \
          ' ' * 9 + 'Явный' + ' ' * 9 + '|' + \
          '\n' + '̅ ' * 55)


def main():
    print_head()

    x, x_max, y = 0, 2 + 10e-5, 0
    h = 10 ** -6

    num = int((x_max - x) / h)
    x_args = count_x_args(x, x_max, h)

    res_runge_kutta = runge_kutta(x, y, h, num)
    res_euler_explicit = euler_explicit(x, y, h, num)
    #res_euler_implicit = euler_implicit(x, y, h, num)
    res_picard_1 = picard_1(x_args)
    res_picard_2 = picard_2(x_args)
    res_picard_3 = picard_3(x_args)
    res_picard_4 = picard_4(x_args)

    f = open('results.csv', 'w')
    for i in range(len(x_args)):
        if not i % 50000: #100000
            print('{:9.2f}|{:11.2e}|{:11.2e}|{:11.2e}|{:11.2e}|{:11.2e}|{:11.4f}'.format(x_args[i],
                                                                                      res_picard_1[i],
                                                                                      res_picard_2[i],
                                                                                      res_picard_3[i],
                                                                                      res_picard_4[i],
                                                                                      res_euler_explicit[i],
                                                                                      res_runge_kutta[i]))
            f.write('{:9.2f},{:11.2e},{:11.2e},{:11.2e},{:11.2e},{:11.2e},{:11.2e}\n'.format(x_args[i],
                                                                                      res_picard_1[i],
                                                                                      res_picard_2[i],
                                                                                      res_picard_3[i],
                                                                                      res_picard_4[i],
                                                                                      res_euler_explicit[i],
                                                                                      res_runge_kutta[i]))


if __name__ == '__main__':
    main()
