from tabulate import tabulate
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def load_csv(filename: str) -> pd.DataFrame:
    return pd.read_csv(filename)


def main():
    # bisection method
    # https://docs.python.org/3/library/bisect.html
    # use scipy optimize algorithms

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

    # loading .csv file
    filename = "data.csv"
    df = load_csv(filename)

    # dividing dataset into X and Y sets
    X = df['U']
    Y = df['a']


if __name__ == '__main__':
    main()
