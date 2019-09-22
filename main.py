import math
from scipy import sin, cos, sqrt
from scipy.integrate import quad
a = -1
b = 0
c = 0.5
d = 1.5
m = 15
eps = 10**-3
t_range = [c + j*(d - c)/m for j in range(m+1)]


def f(x: float, t: float) -> float:
    return cos(t + x**2 + 0.5) / (2 + sin(x**2 + 1))


def s_n(f: 'func', t: float, a: float, b: float, n: int) -> float:
    """
    Вычесляет квадратуру функции f методом трапеций
    на отрезке [a, b] в точке t
    с разбиением [a, b] на n точек
    """
    h = (b - a) / n
    return sum([h*(f(a + (j - 1)*h, t) + f(a + j*h, t))/2
               for j in range(1, n+1)])


def double_n(t: float) -> tuple:
    """
    Метод удвоения числа шагов для стандартной квадратуры
    """
    n = 1
    while abs(s_n(f, t, a, b, n) - s_n(f, t, a, b, 2*n)) > eps:
        n += 1
    ans = s_n(f, t, a, b, n)
    return ans, n


def Gauss(f: 'func', t: float, a: float, b: float) -> float:
    """
    Квадратруа Гаусса
    """
    _c = [i*(b - a)/2 for i in [5/9, 8/9, 5/9]]
    _x = [i*(b - a)/2 + (b + a)/2 for i in [-(3/5)**0.5, 0, (3/5)**0.5]]
    return sum([c*f(x, t) for c, x in zip(_c, _x)])


for t in t_range:
    I_1 = double_n(t)
    I_2 = quad(f, a, b, args=(t,))
    I_3 = Gauss(f, t, a, b)
    print(('|' + '-'*10)*3 + '|')
    print('|' + '{0:10.7f}'.format(I_1[0]), end='')
    print('|' + '{0:10.7f}'.format(I_2[0]), end='')
    print('|' + '{0:10.7f}'.format(I_3) + '|')
