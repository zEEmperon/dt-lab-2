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

    # loading .csv file
    filename = "data.csv"
    df = load_csv(filename)

    # dividing dataset into X and Y sets
    X = df['U']
    Y = df['a']


if __name__ == '__main__':
    main()
