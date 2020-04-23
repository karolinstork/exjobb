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

import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning) #ignores if chains are non continous



def find_representative(cluster_file, naccess_file):
    files_1 = []
    files_2 = []
    cluster_file = open(cluster_file, 'r')
    lines = cluster_file.readlines()
    tuples_of_families = itertools.combinations_with_replacement(lines, 2)
    print("combinations created")

    for family_tuple in tuples_of_families:
        dict={}
        trunc_dict={}
        reso_dict = {}

        print("checking for matches . . . ")
        string1 =family_tuple[0]
        string2 =family_tuple[1]

        files_group1_list =string1.split()
        files_group2_list =string2.split()
        string_combinations = list(itertools.product(files_group1_list, files_group2_list)) #making all binary combinations 1quuA.pdb1

        for string_combination in string_combinations:
            file1 = string_combination[0]
            file2 = string_combination[1]
            complex = compare_strings(file1, file2)

            if complex:
                area = get_interface_area(file1, file2)
                if area != None:
                    print("appending to dict ", file1, file2, area)
                    dict[(file1, file2)] = area

        biggest_area = max(dict.values())
        threshold = biggest_area*0.7
        trunc_dict = {k: v for k, v in dict.items() if v < threshold}

        for (file1, file2) in dict.keys():
            resolution = get_resolution(file1)
            reso_dict[(file1, file2)] = resolution

        representative = max(reso_dict.items(), key=operator.itemgetter(1))[0]

    cluster_file.close()

    return

def compare_strings(file1, file2):
    complex = False
    file1 = (file1.split("."))[0]
    file2 = (file2.split("."))[0]
    if file1[:-1] == file2[:-1] and file1[-1] != file2[-1]:
        complex = True

    return complex

def get_interface_area(file1, file2):
    naccess_file = open("/proj/wallner/users/x_karst/exjobb/result_naccess_hc.txt", "r")
    text = naccess_file.read()

    file1_parts = file1.split(".")
    file2_parts = file2.split(".")


    file1_file2_area_list = re.findall(rf"({file1_parts[0]}[\d].{file1_parts[1]})[\t]({file2_parts[0]}[\d].{file1_parts[1]})[\t]([\d.]+)", text)

    if len(file1_file2_area_list)<1:
        print("no interface area was found between ", file1, file2)

    #else:
        #find biggest area?


    naccess_file.close()

    return


def get_resolution(file1):
    dir = file1[1:3]
    path = "/proj/wallner/share/PDB/191015_biounit/"+dir+file1
    parser=Bio.PDB.PDBParser(PERMISSIVE=1)
    structure = parser.get_structure(file1, path)
    resolution = structure.header["resolution"]

    return resolution

def main():
    args = sys.argv[1:]
    try:
        cluster_file = args[0]
        naccess_file = args[1]
    except IndexError:
        print("ERROR: <file> <cluster_file> <naccess_file>")



    dataset_file = open("dataset_file.txt", "w")

    find_representative(cluster_file, naccess_file)

    dataset_file.close()







if __name__ == '__main__':
    main()
