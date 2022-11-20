import math
from tabulate import tabulate
import matplotlib.pyplot as plt
import numpy as np

# Variant 14
Mx = 15
My = 25
sigma_x = 1
sigma_y = 0.5
R1 = 0.5
R2 = -0.9
X_arr = range(13, 17)


def get_m_y_div_x(x, r):
    return sigma_x * My + (r * sigma_y * (x - Mx)) / sigma_x


def get_d_delta_y(r):
    return sigma_y ** 2 * (1 - r ** 2)


def get_d_y_div_x(r):
    return sigma_y ** 2 * (1 - r ** 2)


def get_W_x_y(x, y, r):
    return (1 / (2 * math.pi * sigma_x * sigma_y * math.sqrt(1 - r ** 2))) \
           * math.exp((-1) / (2 * (1 - (r ** 2)))) * (
                   ((x - Mx) * (x - Mx)) / (sigma_x ** 2) - (2 * r * (x - Mx) * (y - My)) / (sigma_x * sigma_y) + (
                   (y - My) * (y - My)) / (sigma_y ** 2))


def get_W_x(x):
    return (1 / (sigma_x * math.sqrt(2 * math.pi))) * math.exp(-1 * ((x - Mx) * (x - Mx)) / (2 * sigma_x ** 2))


def get_W_y(y):
    return (1 / (sigma_y * math.sqrt(2 * math.pi))) * math.exp(-1 * ((y - My) * (y - My)) / (2 * sigma_y ** 2))


def get_W_y_div_x(x, y, r):
    return (1 / (sigma_y * math.sqrt(2 * math.pi + (1 - r ** 2)))) \
           * math.exp((-1 / (2 * sigma_y ** 2 * (1 - r ** 2))) * (y - My - (sigma_y * (x - Mx) / sigma_x)))


def main():
    # D[y]
    label = "Залежність дисперсії похибки прогнозування від зміни коефіцієнту кореляції"
    r_arr = list(map(lambda r: r / 10, range(0, 11)))
    d_delta_y_arr = list(map(lambda r: round(r * get_d_delta_y(r / 10), 1), range(0, 11)))

    col_names = ["r", "D(delta y)"]
    table_data = np.vstack((r_arr, d_delta_y_arr)).T

    print("\n", label)
    print(tabulate(table_data, headers=col_names, tablefmt="fancy_grid"))

    plt.plot(r_arr, d_delta_y_arr)
    plt.title(label)
    plt.xlabel('r')
    plt.ylabel('D(delta y)')
    plt.show()

    # utils
    sigma_y_coefs = list(set(np.array([*map(lambda x: [-x, x], [0, 1, 2, 3, 4, 5, 10, 15])]).flatten()))
    sigma_y_coefs.sort()
    sigma_y_label_parts = list(map(
        lambda coef: (' + ' if coef > 0 else ' - ') + str(abs(round(coef * 0.2, 1))) + ' * sigma_y',
        sigma_y_coefs
    ))
    get_y_labels = lambda expression: list(map(lambda x:
                                               expression + (x if x != ' - 0.0 * sigma_y' else ''),
                                               sigma_y_label_parts))

    # W(y)
    label = "Значення кривої безумовної густини розподілу прогнозованого параметра W(y)"

    y_arr = list(map(lambda x: My + x * 0.2, sigma_y_coefs))
    W_y_arr = list(map(lambda y: get_W_y(y), y_arr))
    y_labels = get_y_labels("My")

    col_names = ["Формула для y", "y", "W(y)"]
    table_data = np.vstack((y_labels, y_arr, W_y_arr)).T

    print()
    print(label)
    print(tabulate(table_data, headers=col_names, tablefmt="fancy_grid"))

    plt.plot(y_arr, W_y_arr)
    plt.title(label)
    plt.xlabel('y')
    plt.ylabel('W(y)')
    plt.show()

    y_arr_for_non_conditional_distribution = y_arr
    w_y_arr_to_compare = W_y_arr

    # W(y/x(j))
    label = "Значення кривих густин умовного розподілу W(y/x(j))"

    calculate_m_y_div_x = lambda x, r: list(map(lambda sigma_y_coef:
                                                round(get_m_y_div_x(x, r) + sigma_y_coef * 0.2, 2),
                                                sigma_y_coefs))

    x_and_r_arr = list(map(lambda x:
                           [x, R1, x, R2],
                           X_arr))

    x_and_r_arr = np.array(x_and_r_arr).flatten().reshape((-1, 2)).tolist()
    x_and_r_arr = [*x_and_r_arr[::2], *x_and_r_arr[::-2][::-1]]

    y_arr = list(map(lambda params:
                     calculate_m_y_div_x(x=params[0], r=params[1])
                     , x_and_r_arr))

    w_y_div_x_results = list()
    for i in range(len(x_and_r_arr)):
        x = x_and_r_arr[i][0]
        r = x_and_r_arr[i][1]
        w_y_div_x_results.append(list(map(lambda y: get_W_y_div_x(x, y, r), y_arr[i])))

    y_labels = ["", *get_y_labels("M[y/x(j)]")]
    x_and_r_labels = list(map(lambda params: "x = {}, r = {}".format(*params), x_and_r_arr))

    col_names = y_labels
    table_data = list()

    for i in range(len(w_y_div_x_results)):
        table_data.append([x_and_r_labels[i], *w_y_div_x_results[i]])

    print()
    print(label)
    print(tabulate(table_data, headers=col_names, tablefmt="fancy_grid"))

    display_x_arr = X_arr

    result_column = np.array(w_y_div_x_results).T[0]
    w_y_div_x_r1 = result_column[0:len(result_column) // 2]

    plt.plot(display_x_arr, w_y_div_x_r1)
    plt.show()

    y_arr_for_conditional_distribution = calculate_m_y_div_x(x=X_arr[0], r=R1)
    w_y_div_x_to_compare = list(map(lambda y: get_W_y_div_x(X_arr[0], y, R1), y_arr_for_conditional_distribution))

    plt.plot(y_arr_for_non_conditional_distribution, w_y_arr_to_compare, label="W(y)")
    plt.plot(y_arr_for_conditional_distribution, w_y_div_x_to_compare, label="W(y/x(0))")
    plt.xlabel("y")
    plt.ylabel("Густина розподілу")
    plt.title("Порівняння умовної і безумовної густин розподілу")
    plt.legend()
    plt.show()

    # W(x,y)
    x = X_arr[2]
    r = R1
    y_arr = range(15, 36)
    w_x_y_results = list(map(lambda y: get_W_x_y(x, y, r), y_arr))

    plt.plot(y_arr, w_x_y_results, label="W(x, y), де x = {} і r = {}".format(x, r))
    plt.legend()
    plt.show()


if __name__ == '__main__':
    main()
