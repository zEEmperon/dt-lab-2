import math
from tabulate import tabulate
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Variant 14
Mx = 15
My = 25
sigma_x = 1
sigma_y = 0.5
R1 = 0.5
R2 = -0.9
X_1 = 16
X_2 = 15
X_3 = 14
X_4 = 13


def m_y_div_x(x, r):
    return Mx + r * (sigma_y / sigma_x) * (x - Mx)


def d_delta_y(r):
    return sigma_y ** 2 * (1 - r ** 2)


def d_y_div_x(r):
    return sigma_y ** 2 * (1 - r ** 2)


def get_W_x_y(x, y, r):
    return (1 / (2 * math.pi * sigma_x * sigma_y * math.sqrt(1 - r**2))) \
            * math.exp((-1 / 2 * (1 - r ** 2)
                       * abs((x - Mx) ** 2 / sigma_x ** 2
                             - ((2 * r * (x - Mx) * (y - My)) / (sigma_x * sigma_y))
                             + ((y - My) ** 2 / sigma_y ** 2)
                             )))


def get_W_x(x):
    return (1 / (sigma_x * math.sqrt(2 * math.pi))) \
           * math.e ** (-((x - Mx)**2) / 2 * sigma_x ** 2)


def get_W_y(y):
    return (1 / (sigma_y * math.sqrt(2 * math.pi))) \
           * math.e ** (-((y - My)**2) / 2 * sigma_y ** 2)


def get_W_y_div_x(x, y, r):
    return (1 / (sigma_x * math.sqrt(1 - r**2) * math.sqrt(2 * math.pi))) \
           * math.exp(-(y * ((sigma_x * My + sigma_y * r * (x - Mx)
                              / sigma_x) / sigma_x)) / (2 * sigma_y ** 2 * (1 - r ** 2)))


def main():
    # bisection method
    # https://docs.python.org/3/library/bisect.html
    # use scipy optimize algorithms

    # D[y]
    label = "Залежність дисперсії похибки прогнозування від зміни коефіцієнту кореляції"
    r_arr = list(map(lambda r: r / 10, range(0, 11)))
    d_delta_y_arr = list(map(lambda r: round(r * d_delta_y(r / 10), 1), range(0, 11)))

    col_names = ["r", "D(delta y)"]
    table_data = np.vstack((r_arr, d_delta_y_arr)).T

    print("\n", label)
    print(tabulate(table_data, headers=col_names, tablefmt="fancy_grid"))

    plt.plot(r_arr, d_delta_y_arr)
    plt.title(label)
    plt.show()

    # W(y)


if __name__ == '__main__':
    main()
