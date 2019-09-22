import math
from scipy import sin, cos, sqrt
from scipy.integrate import quad
import matplotlib.pyplot as plt
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


def Gauss(f: 'func', t: float, a: float, b: float, n4: bool = False) -> float:
    """
    Квадратруа Гаусса
    """
    if n4:
        _c = [0.347855, 0.652145, 0.652145, 0.347855]
        _x = [-0.861136, -0.339981, 0.339981, 0.861136]
    else:
        _c = [5/9, 8/9, 5/9]
        _x = [-(3/5)**0.5, 0, (3/5)**0.5]
    _c = [i*(b - a)/2 for i in _c]
    _x = [i*(b - a)/2 + (b + a)/2 for i in _x]
    return sum([c*f(x, t) for c, x in zip(_c, _x)])


def main():
    plot1 = []
    plot2 = []
    plot3 = []
    print('|' + 'Trapezoid'.rjust(10, ' '), end='')
    print('|' + 'Gauss'.rjust(10, ' '), end='')
    print('|' + 'SciPy'.rjust(10, ' ') + '|')
    for t in t_range:
        I_1 = double_n(t)
        I_2 = Gauss(f, t, a, b)
        I_3 = quad(f, a, b, args=(t,))
        plot1.append(I_1[0])
        plot2.append(I_2)
        plot3.append(I_3[0])
        print(('|' + '-'*10)*3 + '|')
        print('|' + '{0:10.7f}'.format(I_1[0]), end='')
        print('|' + '{0:10.7f}'.format(I_2), end='')
        print('|' + '{0:10.7f}'.format(I_3[0]) + '|')
    d = plt.scatter(t_range, plot1, s=1, color='r')
    g = plt.scatter(t_range, plot2, s=1, color='b')
    q = plt.scatter(t_range, plot3, s=1, color='g')
    plt.legend((d, g, q),
               ('Trapezoid', 'Gauss', 'SciPy'),
               scatterpoints=1,
               loc='lower left',
               ncol=3,
               fontsize=8)
    plt.xlabel('t')
    plt.ylabel('F(t)')
    plt.show()

if __name__ == '__main__':
    main()
