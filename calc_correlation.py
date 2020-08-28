#!/usr/bin/env python
#SBATCH -t 7-00:00:00
#SBATCH -N 1

import sys
import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns


def calc_mean(data_1, data_2):
    df_1 = pd.read_csv(data_1)
    df_2 = pd.read_csv(data_2)
    print(type(df_1))

    print(df_1)
    print(df_2)


    #heatmap = sns.heatmap(df, linewidth = 0.5, cmap = "Blues").set_title(title)


    array_1 = df_1.to_numpy()
    array_2 = df_2.to_numpy()


    array_1 = np.delete(array_1, 0, axis = 1)
    array_1 = array_1.astype(float)
    array_2 = np.delete(array_2, 0, axis = 1)
    array_2 = array_2.astype(float)


    m1 = np.asmatrix(array_1)
    m2 = np.asmatrix(array_2)

    v1 = m1.flatten()
    v2 = m2.flatten()

    print(v1)
    print(v2)
    print(np.shape(v1))

    similarity_matrix = np.corrcoef(v1, v2)
    print(similarity_matrix)




    return








def main():

    try:
        args = sys.argv[1:]
        data_1 = args[0] #my results
        data_2 = args[1]

    except IndexError:
        print("INPUT ERROR: <file1.csv> <file2.csv>")
        return

    calc_mean(data_1, data_2)











if __name__ == '__main__':
    main()
