import math
from scipy import sin, cos
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


def Gauss(t: float) -> None:
    pass


for t in t_range:
    I_1 = double_n(t)
    I_2 = quad(f, a, b, args=(t,))
    print(I_1[0], I_2[0])


exit(0)

print('-'*8)
for _ in range(8):
    print('|' + '_'*6 + '|')
print('-'*8)
