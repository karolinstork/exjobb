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




def remove_none(protein_peptide_file, cluster_dict, resolution_dict):
    protein_peptide_file = open(protein_peptide_file, "r")
    lines = protein_peptide_file.readlines()
    cluster_complexes_dict = defaultdict(list)

    for line in lines:
        if "None" not in line:
            chain_chain_area_cluster_cluster_ligand_list = line.split()
            ligand = chain_chain_area_cluster_cluster_ligand_list[-1]
            chain1 = chain_chain_area_cluster_cluster_ligand_list[0]
            chain2 = chain_chain_area_cluster_cluster_ligand_list[1]
            area = chain_chain_area_cluster_cluster_ligand_list[2]
            receptor = None

            if ligand == chain1:
                receptor = chain2
            if ligand == chain2:
                receptor = chain1

            pdb_chain_model = receptor.split(".")[0]
            file_extension = receptor.split(".")[1]
            pdb_chain = pdb_chain_model[:5]

            print("receptor:", pdb_chain)
            print("ligand:", ligand)

            cluster = cluster_dict[pdb_chain]


            resolution_key = pdb_chain[:4] + ".pdb"
            print(resolution_key)

            try:
                resolution = resolution_dict[resolution_key]
                if resolution == None:
                    resolution = float(math.inf) #if file exists but has None as resolution

            except:
                resolution = float(math.inf) #if file doesnt exist in pickle dict

            print(type(resolution))

            complex_data = (chain1, chain2, float(area), resolution)
            print(complex_data)
            cluster_complexes_dict[cluster].append(complex_data)


    return cluster_complexes_dict




def find_representative(cluster_complexes_dict, output_file):

    for cluster in cluster_complexes_dict.keys():
        print("-----------------------------------------")
        print("Cluster: ", cluster)

        area_list = []
        dict_90 = defaultdict(list)
        list_of_complexes = cluster_complexes_dict[cluster]

        for complex_data in list_of_complexes:
            area = complex_data[2] #apparently i made the areas in string ..
            area_list.append(area)

        largest_area = max(area_list)
        numb_of_interfaces = len(area_list)
        threshold = 0.9 * largest_area

        print("Number of interfaces: ", numb_of_interfaces)
        print("Largest area: ", largest_area)
        print("Threshold: >=", threshold)

        for complex_data in list_of_complexes:
            if complex_data[2] >= threshold:
                dict_90[cluster].append(complex_data)


        try:
            representative_complex = min(dict_90[cluster], key=operator.itemgetter(3)) #lowest resolution aka best reoslution
            print("Representative complex: ", representative_complex)
            first_chain = representative_complex[0]
            second_chain = representative_complex[1]
            representative_complex = str(first_chain)+'\t'+str(second_chain)+'\n'
            output_file.write(representative_complex)
        except TypeError:
            print("ERROR")
            print(dict_90[cluster])






    return






def main():

    try:
        args = sys.argv[1:]
        protein_peptide_file = args[0]


    except IndexError:
        print("INPUT ERROR: <protein_peptide_complexes_list>")
        return


    cluster_dict = pickle.load(open("/proj/wallner/users/x_karst/exjobb/pickles/file_cluster_dict.p", "rb"))
    resolution_dict = pickle.load(open("/proj/wallner/users/x_karst/exjobb/pickles/resolutions_sbatch.p", "rb"))


    cluster_complexes_dict = remove_none(protein_peptide_file, cluster_dict, resolution_dict)
    #print(cluster_complexes_dict)

    output_file = open("/proj/wallner/users/x_karst/exjobb/protein_peptide_data/protein_peptide_dataset.out", "w")
    find_representative(cluster_complexes_dict, output_file)

    output_file.close()








if __name__ == '__main__':
    main()
