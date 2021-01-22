#!/usr/bin/env python
#SBATCH -t 7-00:00:00
#SBATCH -N 1

import sys
import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

import pickle
from collections import defaultdict
import math
import operator




def remove_none(protein_protein_file, resolution_dict):
    protein_protein_file = open(protein_protein_file, "r")
    lines = protein_protein_file.readlines()
    cluster_complexes_dict = defaultdict(list)

    for line in lines:
        if "None" not in line:
            chain_chain_area_cluster_cluster_list = line.split()
            chain1 = chain_chain_area_cluster_cluster_list[0]
            chain2 = chain_chain_area_cluster_cluster_list[1]
            area = chain_chain_area_cluster_cluster_list[2]
            cluster1 = chain_chain_area_cluster_cluster_list[3]
            cluster2 = chain_chain_area_cluster_cluster_list[4]

            pdb_chain = chain1.split(".")[0]
            resolution_key = pdb_chain[:4] + ".pdb"
            #print(resolution_key)

            try:
                resolution = resolution_dict[resolution_key]
                if resolution == None:
                    resolution = float(math.inf) #if file exists but has None as resolution

            except:
                resolution = float(math.inf) #if file doesnt exist in pickle dict


            complex_data = (chain1, chain2, float(area), resolution)

            if (cluster2, cluster1) in cluster_complexes_dict.keys():
                clusters = (cluster2, cluster1)
                cluster_complexes_dict[clusters].append(complex_data)
            else:
                clusters = (cluster1, cluster2)
                cluster_complexes_dict[clusters].append(complex_data)




    return cluster_complexes_dict




def find_representative(cluster_complexes_dict, output_file):

    for cluster_combo in cluster_complexes_dict.keys():
        print("-----------------------------------------")
        print("Family: ", cluster_combo)

        area_list = []
        dict_90 = defaultdict(list)
        list_of_complexes = cluster_complexes_dict[cluster_combo]

        for complex_data in list_of_complexes:
            area = complex_data[2]
            area_list.append(area)

        largest_area = max(area_list)
        numb_of_interfaces = len(area_list)
        threshold = 0.9 * largest_area

        print("Number of interfaces: ", numb_of_interfaces)
        print("Largest area: ", largest_area)
        print("Threshold: >=", threshold)

        for complex_data in list_of_complexes:
            if complex_data[2] >= threshold:
                dict_90[cluster_combo].append(complex_data)


        try:
            representative_complex = min(dict_90[cluster_combo], key=operator.itemgetter(3)) #lowest resolution aka best reoslution
            print("Representative complex: ", representative_complex)
            first_chain = representative_complex[0]
            second_chain = representative_complex[1]
            representative_complex = str(first_chain)+'\t'+str(second_chain)+'\n'
            output_file.write(representative_complex)
        except TypeError:
            print("ERROR")
            print(dict_90[cluster_combo])






    return






def main():

    try:
        args = sys.argv[1:]
        protein_protein_file = args[0]


    except IndexError:
        print("INPUT ERROR: <protein_protein_complexes_list>")
        return



    resolution_dict = pickle.load(open("/proj/wallner/users/x_karst/exjobb/pickles/resolutions_sbatch.p", "rb"))


    cluster_complexes_dict = remove_none(protein_protein_file, resolution_dict)

    #print(cluster_complexes_dict)

    output_file = open("/proj/wallner/users/x_karst/exjobb/protein_peptide_data/protein_protein_dataset.out", "w")
    find_representative(cluster_complexes_dict, output_file)
    list1 = list(cluster_complexes_dict.keys())
    myset = set(list1)
    print(len(list1))
    print(len(myset))

    output_file.close()








if __name__ == '__main__':
    main()
