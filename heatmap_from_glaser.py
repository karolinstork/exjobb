#!/usr/bin/env python

import sys
import os
import pandas as pd
import seaborn as sns
import csv
import matplotlib.pyplot as plt
from collections import OrderedDict

def main():
    residue_order = ["R", "K", "N", "D", "Q", "E", "H", "P", "Y", "W", "S", "T", "G", "A", "M", "C", "F", "L", "V", "I"] # to sort the same way the article did
    ordered_dict = OrderedDict({key: OrderedDict({key: None for key in residue_order}) for key in residue_order})
    position_dict = {0:"I", 1: "V", 2:"L", 3:"F", 4:"C", 5: "M", 6: "A", 7: "G", 8: "T", 9: "S", 10: "W", 11: "Y", 12:"P", 13: "H", 14: "E", 15:"Q", 16: "D" , 17: "N", 18:"K", 19: "R"}

    file =  open('/proj/wallner/users/x_karst/exjobb/tables/gij_v_glaser_table.csv', 'r')
    lines = file.readlines()
    lines = lines[1:]

    for line in lines:
        line = line.strip("\n")
        line = line.lstrip(",")
        print(line)

        master_res = line.split(",")[0]
        all_values_for_that_row = line.split(",")[1:]
        column = 0
        for value in all_values_for_that_row:
            mini_res = position_dict[column]
            value = float(value)

            print(master_res, mini_res, value, type(value))


            ordered_dict[master_res][mini_res] = value
            column = column + 1


    df= pd.DataFrame(ordered_dict)
    print(df)

    title = "Gij(v) results from Glaser et al study"
    heatmap= sns.heatmap(df, linewidth = 0.5, cmap = "seismic", center = 0, vmin = -15, vmax = 25).set_title(title)
    plt.savefig("/proj/wallner/users/x_karst/exjobb/pictures/gijv_glaser.png", bbox_inches='tight')

    plt.show()















if __name__ == '__main__':
    main()
