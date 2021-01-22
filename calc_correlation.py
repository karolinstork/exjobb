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
    print(data_1, data_2)

    if ("glaser" not in data_1) and ("glaser" not in data_2): #works for pro vs pep and pro vs whole database
        print("whole database vs true protein protein")
        df_1 = pd.read_csv(data_1)
        df_2 = pd.read_csv(data_2)

        if "pro" not in data_1:
            del df_1['Unnamed: 0']
        if "pro" not in data_2:
            del df_2['Unnamed: 0']

        upper_tri1 = df_1.where(np.tril(np.ones(df_1.shape), 0).astype(bool)).stack()
        upper_tri2 = df_2.where(np.tril(np.ones(df_2.shape), 0).astype(bool)).stack()

        print(upper_tri1)
        print(upper_tri2)

        arr1 = upper_tri1.to_numpy()
        arr2 = upper_tri2.to_numpy()

        print(arr1)
        print(arr2)
        print(arr1.shape)
        print(arr2.shape)


        similarity_matrix = np.corrcoef(arr1, arr2)
        print(similarity_matrix)




    else: #comparing glaser and my results with the whole database
        print("comparing glaser with something...")

        df_1 = pd.read_csv(data_1) #my first results
        df_2 = pd.read_csv(data_2)

        if "glaser" in data_1:
            del df_1['Unnamed: 0'] #removing letters
        if "glaser" in data_2:
            del df_2['Unnamed: 0']

        print(df_1)
        print(df_2)

        upper_tri1 = df_1.where(np.tril(np.ones(df_1.shape), 0).astype(bool)).stack()
        upper_tri2 = df_2.where(np.tril(np.ones(df_2.shape), 0).astype(bool)).stack()

        print(upper_tri1)
        print(upper_tri2)

        arr1 = upper_tri1.to_numpy()
        arr2 = upper_tri2.to_numpy()

        print(arr1)
        print(arr2)
        print(arr1.shape)
        print(arr2.shape)


        similarity_matrix = np.corrcoef(arr1, arr2)
        print(similarity_matrix)




    return








def main():

    try:
        args = sys.argv[1:]
        data_1 = args[0]
        data_2 = args[1]

    except IndexError:
        print("INPUT ERROR: <file1.csv> <file2.csv>")
        return

    calc_mean(data_1, data_2)











if __name__ == '__main__':
    main()
