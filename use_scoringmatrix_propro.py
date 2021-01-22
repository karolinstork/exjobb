#!/usr/bin/env python
#SBATCH -t 7-00:00:00
#SBATCH -N 1

import sys
import os
import pandas as pd
from os import walk
import subprocess
import Bio
from Bio.PDB import *
import string
from collections import defaultdict
from collections import Counter
import re
import numpy as np


import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning) #ignores if chains are non continous

def read_csv(matrix_file):
    likelihood_score_dict = {}
    df = pd.read_csv(matrix_file)

    for column in df.columns:
        i = 0
        for value in df[column]:
            if ((column, df.columns[i]) and (df.columns[i], column)) not in likelihood_score_dict.keys():
                likelihood_score_dict[(column, df.columns[i])] = value
            i = i + 1

    return likelihood_score_dict





def make_list_of_files(target_path):
    list_of_files = []

    for (dir_path, dirnames, filenames) in walk(target_path):
        for filename in filenames:
            if filename[-3:]=="pdb":
                target_path = dir_path+"/"+filename
                list_of_files.append(target_path)

    print("Number of files in", target_path, ":", len(list_of_files), '\n')
    return list_of_files


def check_file(file):
    pdb_file = open(file, "r")
    s = pd.Series(pdb_file)
    orig_chains = list(s.str.get(21).unique()) #find chain names

    if len(orig_chains)>2:
        print("WARNING: TOO MANY CHAINS IN THIS FILE ")
        should_continue = False
        return


    pdb_text = pdb_file.read()

    numb_models = 0
    all_matches= re.findall("MODEL", pdb_text)
    if len(all_matches)>0:
        print("WARNING: TOO MANY MODELS")
        should_continue = False
        return

    should_continue = True


    return should_continue


def find_contacts(file):

    pdb_file = open(file, "r")
    pdb_text = pdb_file.read()

    filename = os.path.basename(os.path.normpath(file))
    if "MOAL" in file:
        output_file = open(f'/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/glaser/MOAL/{filename}.output', "w")

    if "CAPRI" in file:
        output_file = open(f'/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/glaser/CAPRI/{filename}.output', "w")

    subprocess.run(["/proj/wallner/users/x_bjowa/git/pdb_codes/contact-maps/./contact_map_chain", file, "-chains", "A", "B", "-cb"], stdout = output_file)

    return


def calc_score(file, likelihood_score_dict, cut_off):
    pro_interface_dict = defaultdict(Counter)
    total_score = 0
    translation_dict = {}
    string_position_to_chain_dict = {}

    filename = os.path.basename(os.path.normpath(file))

    if "MOAL" in file:
        output_read = open(f"/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/glaser/MOAL/{filename}.output", "r")
    if "CAPRI" in file:
        output_read = open(f"/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/glaser/CAPRI/{filename}.output", "r")


    lines = output_read.readlines()

    string_pos = 0

    for line in lines: #translation_dict
        if line[:3] == "RES":

            line_split_by_space = line.split()
            pos = (line_split_by_space[2])
            res = (line_split_by_space[3])
            chain = (line_split_by_space[4]).strip(":")
            pos_chain = pos+chain

            string_position_to_chain_dict[string_pos] = pos_chain #tex row 0 = 1C
            translation_dict[pos_chain] = res #tex 1C = glycine

            string_pos = string_pos + 1



    for line in lines:
        contact_with_anything_across = None
        if line[:3] == "RES":
            line_split_by_space = line.split()
            pos = (line_split_by_space[2])
            res = (line_split_by_space[3])
            chain = (line_split_by_space[4]).strip(":")

            pattern = f":(.*)"
            all_contacts_in_string = re.findall(pattern, line) #list with length 1, all distances as a string
            list_of_distances = all_contacts_in_string[0].split() #distances in string format in list

            column = 0

            for distance in list_of_distances:
                distance = float(distance)
                if distance < cut_off:
                    position_chain_id = string_position_to_chain_dict[column]
                    contact_res = translation_dict[position_chain_id]
                    if position_chain_id[-1] != chain: #only intrested in interactions between two different chains
                        contact_residue = translation_dict[position_chain_id]
                        pro_interface_dict[res].update(contact_residue)

                        if contact_residue != res: #to avoid adding extra bonds between two identical residues ie glycin - glycin
                            pro_interface_dict[contact_residue].update(res)


                column = column + 1

    output_read.close()

    if "MOAL" in file:
        os.remove(f"/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/glaser/MOAL/{filename}.output")

    if "CAPRI" in file:
        os.remove(f"/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/glaser/CAPRI/{filename}.output")


    asymmetric_numb_contacts_dict = {} # to NOT have the symmetric heatmap ie same information twice
    for master_residue, small_dict in pro_interface_dict.items():
        for mini_residue, contacts in small_dict.items():
            if ((master_residue, mini_residue) not in asymmetric_numb_contacts_dict.keys()) and ((mini_residue, master_residue) not in asymmetric_numb_contacts_dict.keys()):
                    asymmetric_numb_contacts_dict[(master_residue, mini_residue)] = contacts


    tot_numb_contacts =  sum(asymmetric_numb_contacts_dict.values())
    if tot_numb_contacts == 0:
        total_score = None
        return total_score



    for (master_residue, mini_residue) in asymmetric_numb_contacts_dict.keys(): #calc score
        occurence = asymmetric_numb_contacts_dict[(master_residue, mini_residue)]
        if (master_residue, mini_residue) in likelihood_score_dict.keys(): #likelihood dict has either (A, B) or (B, A)
            value = likelihood_score_dict[(master_residue,mini_residue)]
            score = occurence*value
            total_score = total_score + score
        elif (mini_residue, master_residue) in likelihood_score_dict.keys():
            value = likelihood_score_dict[(mini_residue, master_residue)]
            score = occurence * value
            total_score = total_score + score

    print("Total number of contacts:", tot_numb_contacts)
    total_score = total_score / tot_numb_contacts
    print("Total score:", total_score, '\n')


    return total_score



def write_to_resultfile(file, total_score, result_file):
    filename = os.path.basename(os.path.normpath(file))
    result = filename + '\t' + str(total_score) + '\n'
    result_file.write(result)

    return










def main():

    try:
        args = sys.argv[1:]
        cut_off = int(args[0]) #for the c script
        project_path = args[1] #target directory either MOAL or CAPRI: "/proj/wallner/CnM-dataset/MOAL_Benchmark/" or "/proj/wallner/CnM-dataset/CAPRI_ScoreSet/"


    except IndexError:
        print("USAGE: <c script cutoff> <target directory>")
        return


    matrix_file = "/proj/wallner/users/x_karst/exjobb/tables/gijv_whole.csv"
    likelihood_score_dict = read_csv(matrix_file)


    for (dir_path, dirnames, filenames) in walk(project_path):
        for dir in dirnames:  #for every target in CAPRI or MOAL
            target_path = (dir_path+dir)

            list_of_files = make_list_of_files(target_path)

            break_point = 0
            correct_files = 0
            no_contacts = 0


            if "MOAL" in target_path:
                result_file = open(f"/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/glaser/MOAL/{dir}/myscore_{dir}.txt", "w")
            if "CAPRI" in target_path:
                result_file = open(f"/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/glaser/CAPRI/{dir}/myscore_{dir}.txt", "w")


            for file in list_of_files: #for every model in that target
                should_continue = None
                print(file, break_point)

                should_continue = check_file(file) #checks if only binary file ie 1 model and 2 chains
                if should_continue == True:
                    find_contacts(file)
                    total_score = calc_score(file, likelihood_score_dict, cut_off)
                    if total_score != None: #checks for tot numb contacts not = 0
                        write_to_resultfile(file, total_score, result_file)
                        correct_files = correct_files + 1
                    else:
                        print("Total number of contacts less than 1")
                        no_contacts = no_contacts + 1


                break_point = break_point + 1

            result_file.close()
            print("Number of files:", len(list_of_files))
            print("Correct files for this target:", correct_files)
            print("Number of files without close contacts:", no_contacts)
            print("-------------------------------")












if __name__ == '__main__':
    main()
