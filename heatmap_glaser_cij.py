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

    tot_numb_interactions_dict = {}

    file =  open('/proj/wallner/users/x_karst/exjobb/tables/cij_glaser_table.csv', 'r')
    lines = file.readlines()
    lines = lines[1:] #remove labels in top row

    for line in lines:
        line = line.strip("\n")
        line = line.lstrip(",")
        master_res = line.split(",")[0]
        all_values_for_that_row = line.split(",")[1:]
        column = 0

        for value in all_values_for_that_row:
            mini_res = position_dict[column]
            value = int(value)

            ordered_dict[master_res][mini_res] = value
            column = column + 1

    print(ordered_dict)

    df= pd.DataFrame(ordered_dict)
    print("Number of contacts")
    print(df)

    title = "Cij results from Glaser et al study"
    heatmap= sns.heatmap(df, linewidth = 0.5, cmap = "Blues").set_title(title)
    plt.savefig("/proj/wallner/users/x_karst/exjobb/pictures/cij_glaser.png", bbox_inches='tight')
    plt.show()


    ########################################################################################################################
    normalized_dict = OrderedDict({key: OrderedDict({key: None for key in residue_order}) for key in residue_order})

    for master_residue, small_dict in ordered_dict.items():
        all_interactions_with_master_residue  = 0
        for mini_res, contacts in small_dict.items():
            if mini_res == master_residue:
                all_interactions_with_master_residue = all_interactions_with_master_residue + contacts
            else:
                all_interactions_with_master_residue = all_interactions_with_master_residue + (contacts/2)
        tot_numb_interactions_dict[master_residue] = all_interactions_with_master_residue




    tot_numb_interactions = sum(tot_numb_interactions_dict.values())
    print("Total number of interactions:", tot_numb_interactions)

    for master_residue, small_dict in ordered_dict.items():
        for mini_residue, contacts in small_dict.items():
            norm_contacts = float(contacts) / tot_numb_interactions
            normalized_dict[master_residue][mini_residue] = norm_contacts




    df_2= pd.DataFrame(normalized_dict)
    print("Normalized number of contacts")
    print(df_2)
    title_2 = "Qij results from Glaser et al study"
    heatmap= sns.heatmap(df_2, linewidth = 0.5, cmap = "Blues").set_title(title_2)
    plt.savefig("/proj/wallner/users/x_karst/exjobb/pictures/qij_glaser.png", bbox_inches='tight')
    plt.show()


















if __name__ == '__main__':
    main()
