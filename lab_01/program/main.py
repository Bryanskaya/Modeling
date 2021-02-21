# Метод Пикара, метод Эйлера, метод Рунге-Кутта

def f(x, y):
    return x ** 2 + y ** 2


# Метод Пикара
def picard_1(x):
    return x ** 3 / 3


def picard_2(x):
    return picard_1(x) + x ** 7 / 63


def picard_3(x):
    return picard_2(x) + x ** 15 / 59535 + 2 * x ** 11 / 2079


def picard_4(x):
    return picard_3(x) + x ** 31 / 109876902975 + 4 * x ** 23 / 99411543 + \
           4 * x ** 27 / 3341878155 + 2 * x ** 23 / 86266215 + \
           2 * x ** 19 / 3393495 + 4 * x ** 19 / 2488563 + 4 * x ** 15 / 93555


def runge_kutta(x, y, h, num):
    alpha = 0.5
    res = []

    temp = h / (2 * alpha)

    for i in range(num):
        res.append([x, y])

        k1 = f(x, y)
        k2 = f(x + temp, y + temp * k1)

        y += h * ((1 - alpha) * k1 + alpha * k2)
        x += h  # or temp

    return res


def euler_explicit(x, y, h, num):
    res = []

    for i in range(num):
        res.append([x, y])

        y += h * f(x, y)
        x += h

    return res


def euler_implicit(x, y, h, num):
    res = []

    



def print_head():
    print(' ' * 4 + 'x' + ' ' * 4 + '|' + \
          ' ' * 9 + 'Метод Пикара' + ' ' * 10 + '|' + \
          ' ' * 3 + 'Метод Эйлера' + ' ' * 3 + '|' + \
          ' ' * 3 + 'Метод Рунге-Кутта' + ' ' * 3 +
          '\n' + \
          ' ' * 9 + '|' + \
          ' ' * 3 + '1' + ' ' * 3 + '|' + \
          ' ' * 3 + '2' + ' ' * 3 + '|' + \
          ' ' * 3 + '3' + ' ' * 3 + '|' + \
          ' ' * 3 + '4' + ' ' * 3 + '|' + \
          ' ' * 18 + '|' + \
          '\n' + '̅ ' * 50)


def main():
    print_head()

    x, y = 0, 0
    h = 0.1

    print('Кутт', runge_kutta(x, y, h, 10))
    print('Явный Эйлер', euler_explicit(x, y, h, 10))

    for x in range(0, 10):
        print('{:9.3f}|{:7.3f}|{:7.3f}|{:7.3f}|{:7.3f}'.format(x,
                                                               picard_1(x),
                                                               picard_2(x),
                                                               picard_3(x),
                                                               picard_4(x)))


if __name__ == '__main__':
    main()
