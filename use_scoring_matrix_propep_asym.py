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
import traceback

import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning) #ignores if chains are non continous




def read_csv(matrix_file):
    likelihood_score_dict = {}
    df = pd.read_csv(matrix_file) #proteins at x axis and peptides at y axis
    print(df)

    for column in df.columns:
        i = 0
        for value in df[column]:
            if (column, df.columns[i]) not in likelihood_score_dict.keys():
                likelihood_score_dict[(column, df.columns[i])] = value
            i = i + 1

    #tuple in dictionary (protein res, peptide res)
    return likelihood_score_dict





def make_list_of_files(target_directory):
    list_of_files = []
    for (dir_path, dirnames, filenames) in walk(f"/proj/wallner/users/x_isaak/ModelsKarolin/models/{target_directory}"):
        for filename in filenames:
            print(dir_path)
            if os.path.basename(dir_path) != "3gz1QA03": # dont have a native file so was not possible to run dockq
                if filename[-3:]=="pdb":
                    path = dir_path+"/"+filename
                    list_of_files.append(path)
                #    print(path)
    print("Number of files:", len(list_of_files), '\n')
    return list_of_files




def check_file(file):
    pdb_file = open(file, "r")
    s = pd.Series(file)
    pdb_text = pdb_file.read()

    orig_chains = list(s.str.get(21).unique()) #find chain names
    if len(orig_chains)>2:
        print("WARNING: TOO MANY CHAINS IN THIS FILE ")
        should_continue = False
        return

    numb_models = 0
    all_matches= re.findall("MODEL", pdb_text) #find number of models
    if len(all_matches)>0:
        print("WARNING: TOO MANY MODELS")
        should_continue = False
        return

    should_continue = True
    return should_continue




def find_contacts(file, peptide, receptor):
    pdb_file = open(file, "r")
    pdb_text = pdb_file.read()

    first_name = file.split("/")[-2]
    last_name = os.path.basename(os.path.normpath(file))

    output_file = open(f'/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/pro_pep/c_script_output/{first_name}/{last_name}.output', "w")
    subprocess.run(["/proj/wallner/users/x_bjowa/git/pdb_codes/contact-maps/./contact_map_chain", file, "-chains", "A", "B", "-cb"], stdout = output_file)

    return




def calc_score(file, likelihood_score_dict, cut_off, peptide, receptor):
    print(file)
    interface_counter = defaultdict(Counter)

    total_score = 0
    translation_dict = {}
    string_position_to_chain_dict = {}

    first_name = file.split("/")[-2]
    last_name = os.path.basename(os.path.normpath(file))

    output_read = open(f"/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/pro_pep/c_script_output/{first_name}/{last_name}.output", "r")
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
                    if chain == receptor: #look only from protein side to avoid adding duplicates in counter
                        if position_chain_id[-1] != chain: #only intrested in interactions between two different chains
                            contact_residue = translation_dict[position_chain_id]
                            interface_counter[res].update(contact_residue)


                column = column + 1

    output_read.close()

    os.remove(f"/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/pro_pep/c_script_output/{first_name}/{last_name}.output")

    interface_dict = {}

    for protein_residue, small_dict in interface_counter.items():
        for peptide_residue, numb_contacts in small_dict.items():
            interface_dict[(protein_residue, peptide_residue)] = numb_contacts




    tot_numb_contacts =  sum(interface_dict.values())
    print("Tot numb contacts with cutoff", cut_off, "Å:", tot_numb_contacts)

    if tot_numb_contacts == 0:
        total_score = None
        return total_score

    for (protein_residue, peptide_residue) in interface_dict.keys(): #calc score
        try:
            value = likelihood_score_dict[(protein_residue, peptide_residue)]
        except KeyError:
            print("ERROR: Unknown residues in this file. File removed.") #some files had an X residue
            total_score = None
            return
        occurence = interface_dict[(protein_residue, peptide_residue)]
        score = occurence*value
        total_score = total_score + score


    print("Total number of contacts:", tot_numb_contacts)
    total_score = total_score / tot_numb_contacts
    print("Total score:", total_score, '\n')

    return total_score


def write_to_resultfile(file, total_score, result_file, target_directory):
    id = target_directory[:4] #4xev
    filename = os.path.basename(os.path.normpath(file)) #output.833.pdb
    numb = filename.split(".")[1] #833

    model = id+"_"+numb


    result = model + '\t' + str(total_score) + '\n'
    result_file.write(result)

    return













def main():
    try:
        args = sys.argv[1:]
        cut_off = int(args[0]) #the threshold distance for a contact between two amino acids
        target_directory = args[1]

    except IndexError:
        print("USAGE: <c script cutoff> <target directory>")
        return


    matrix_file = "/proj/wallner/users/x_karst/exjobb/protein_peptide_data/gij_v_results_pro_pep_asymmetric.csv"
    likelihood_score_dict = read_csv(matrix_file)

    list_of_files = make_list_of_files(target_directory)

    result_file = open(f"/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/pro_pep/myscores_asym/myscore_propep_{target_directory}.txt", "w")

    for file in list_of_files:
        should_continue = None
        should_continue = check_file(file) #checks if only binary file ie 1 model and 2 chains


        if should_continue == True:
            first_name = file.split("/")[-2]
            peptide = first_name[4]
            receptor = first_name[5]

            find_contacts(file, peptide, receptor)
            total_score = calc_score(file, likelihood_score_dict, cut_off, peptide, receptor)

        if total_score != None: #checks for tot numb contacts not = 0 and not unknown residues
            write_to_resultfile(file, total_score, result_file, target_directory)




    result_file.close()













if __name__ == '__main__':
    main()
