#!/usr/bin/env python3
#SBATCH -t 7-00:00:00
#SBATCH -N 1

import sys
import os
import fnmatch
import re
from operator import itemgetter
from collections import defaultdict
import itertools
import string
from Bio.PDB import *
import Bio
import pickle


import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning) #ignores if chains are non continous


def find_representative(cluster_file, naccess_file):
    files_1 = []
    files_2 = []
    cluster_file = open(cluster_file, 'r')
    lines = cluster_file.readlines()
    tuples_of_families = itertools.combinations_with_replacement(lines, 2)
    tuples_of_families = list(tuples_of_families) #allocates a lot of memory
    numb_families = len(tuples_of_families)
    print(numb_families, "families created")
    counter=0
    numb_repr = 0

    for family_tuple in tuples_of_families:
        i = tuples_of_families.index(family_tuple)
        i= i+1
        print("------------family", i, "--------------")
        progress = (float(i/numb_families))*100
        dict={}
        trunc_dict={}
        reso_dict = {}

        counter = counter + 1
        if counter == 10000:
          print("PROGRESS REPORT:",progress, "%")
          print("Number of representatives found:", numb_repr)
          counter = 0

        line1 =family_tuple[0]
        line2 =family_tuple[1]

        files_group1_list =line1.split()
        files_group2_list =line2.split()
        string_combinations = list(itertools.product(files_group1_list, files_group2_list)) #making all binary combinations 1quuA.pdb1-1quuA.pdb1


        for string_combination in string_combinations:
            file1 = string_combination[0]
            file2 = string_combination[1]
            is_complex = compare_strings(file1, file2)

            if is_complex:
                biggest_area_tuple = get_interface_area(file1, file2)
                if biggest_area_tuple[2] != None:
                    print("Complex with confirmed interface found. ", biggest_area_tuple)
                    complex = (biggest_area_tuple[0], biggest_area_tuple[1])
                    area = biggest_area_tuple[2]
                    dict[complex] = area

        if len(dict)>0:
            biggest_area = max(dict.values())
            print("Biggest interface in family: ", biggest_area)
            threshold = biggest_area*0.7
            print("threshold: ", threshold)
            trunc_dict = {k: v for k, v in dict.items() if v > threshold} #only keeping the biggest interfaces

            for complex in trunc_dict.keys():
                file1 = complex[0]
            #    resolution = get_resolution(file1, resolutions)
                resolution = 0.0
                print("resolution: ", resolution)
                reso_dict[complex] = resolution

            representative_tuple = max(reso_dict.items(), key=itemgetter(1))[0]
            print("family representative found: ", representative_tuple)
            print("---------------------------------------------")
            numb_repr = numb_repr + 1
            result = str(representative_tuple[0])+ " "+str(representative_tuple[1])+'\n'
            #dataset_file.write(result)

    cluster_file.close()

    return

def compare_strings(file1, file2):
    is_complex = False
    file1 = (file1.split("."))[0] #files from blastclust
    file2 = (file2.split("."))[0]

    if file1[:-1] == file2[:-1]:
        is_complex = True

    return is_complex

def get_interface_area(file1, file2):
    naccess_file = open("/proj/wallner/users/x_karst/exjobb/naccess_results.txt", "r")
    text = naccess_file.read()
    list=[]

    file1_parts = file1.split(".")
    file2_parts = file2.split(".")

    file1_file2_area_list = re.findall(rf"({file1_parts[0]}[\d].{file1_parts[1]})[\t]({file2_parts[0]}[\d].{file1_parts[1]})[\t]([\d.]+)", text)

    for file1_file2_area in file1_file2_area_list:
        file1 = file1_file2_area[0]
        file2 = file1_file2_area[1]
        file1_split = (file1.split("."))[0]
        file2_split = (file2.split("."))[0]

        if file1_split[-2:] == file2_split[-2:]:
            file1_file2_area_list.remove(file1_file2_area)
            print(file1_file2_area," removed due to same chain, same model")


    if len(file1_file2_area_list)>1:
        #print("interfaces found: ", file1_file2_area_list)
        for tuple in file1_file2_area_list:
            file1_file2_area =(tuple[0], tuple[1], float(tuple[2]))
            list.append(file1_file2_area)

        biggest_area_tuple = max(list,key=lambda item:item[2])

    else:
        biggest_area_tuple = (file1, file2, None)

    naccess_file.close()

    return biggest_area_tuple


def get_resolution(file1, resolutions):
    dir = file1[1:3]
    pdb_id = (file1.split("."))[0]
    pdb_id = pdb_id[:-2]
    pdb_file = pdb_id+"."+(file1.split("."))[1]

    resolution = resolutions[pdb_file]


    return resolution

def main():
    args = sys.argv[1:]
    try:
        cluster_file = args[0]
        naccess_file = args[1]
    except IndexError:
        print("ERROR: <file> <cluster_file> <naccess_file>")


    #resolutions = pickle.load(open("resolutions.p", "rb"))

    #dataset_file = open("dataset_file.txt", "w")
#    find_representative(cluster_file, naccess_file, dataset_file, resolutions)

    find_representative(cluster_file, naccess_file)

    #dataset_file.close()







if __name__ == '__main__':
    main()
