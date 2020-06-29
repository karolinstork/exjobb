#!/usr/bin/env python
#SBATCH -t 7-00:00:00
#SBATCH -N 1

import sys
import os
import subprocess
import string
import re
from collections import defaultdict
from collections import Counter
import pandas as pd
import seaborn as sns



def calculate_contacts(complex, res_dict):
    line = complex.strip('\n')
    files = line.split('\t')
    file_1 = files[0]
    file_2 = files[1]

    file1_parts = file_1.split(".")
    file2_parts = file_2.split(".")
    base_file = (file1_parts[0])[:4]+"."+file1_parts[1]
    pdb_id = (file1_parts[0])[:4]

    chain_1 = (file1_parts[0])[4:]
    chain_2 = (file2_parts[0])[4:]

    only_1_model = False

    if (chain_1[1]=="1") and (chain_2[1]=="1") and (len(file1_parts[0])==6) and (len(file2_parts [0]) == 6):
        only_1_model = True

    if only_1_model:

        #print(base_file, file_1, file_2, chain_1, chain_2)
        output_file = open(f"/proj/wallner/users/x_karst/exjobb/c_script_output/{pdb_id}.output", "w")
        dir = base_file[1:3]
        target = f'/proj/wallner/share/PDB/191015_biounit/{dir}/{base_file}'
        subprocess.run(["/proj/wallner/users/x_bjowa/git/pdb_codes/contact-maps/./contact_map_chain", target, "-chains", chain_1[0], chain_2[0], "-c", "6", "-cb"], stdout = output_file)
        output_file.close()

        output_read = open(f"/proj/wallner/users/x_karst/exjobb/c_script_output/{pdb_id}.output", "r")
        #lines = output_read.readlines()
        text = output_read.read()
        lines = text.split('\n')


        for line in lines:
            if line[:3] == "RES":
                line_split_by_space = line.split()
                res = (line_split_by_space[3])
                chain = (line_split_by_space[4]).strip(":")
                #print(res, chain, chain_1[:-1], chain_2[:-1])
                #print(line)

                if chain == chain_1[:-1]:       #only interesed in contacts between 2 chains not in the same chain
                    pattern = f"[\d]+{chain_2[:-1]}"
                    contacts = re.findall(pattern, line)
                    if len(contacts)>0:
                        #print(contacts)
                        for contact in contacts:    #find out what residue the contacts are
                            contact_chain = contact[-1:]
                            contact_pos = contact[:-1]
                            #print(contact_chain, contact_pos)
                            pattern1 = f'RES[\s]+[\d]+[\s]+{contact_pos}[\s]+([\w])[\s]+{contact_chain}:'
                            #print(pattern1)
                            match = re.search(pattern1, text)
                            if match:
                                #print("REAL AMINO ACID IS: ", match.group(1))
                                res_dict[res].update(match.group(1))



                        #res_dict[res].append(contacts)

                if chain == chain_2[:-1]:
                    pattern = f"[\d]+{chain_1[:-1]}"
                    contacts = re.findall(pattern, line)
                    if len(contacts)>0:
                        #print(contacts)
                        for contact in contacts:    #find out what residue the contacts are
                            contact_chain = contact[-1:]
                            contact_pos = contact[:-1]
                            #print(contact_chain, contact_pos)
                            pattern2 = f'RES[\s]+[\d]+[\s]+{contact_pos}[\s]+([\w])[\s]+{contact_chain}:'
                            #print(pattern1)
                            match = re.search(pattern2, text)
                            if match:
                                #print("REAL AMINO ACID IS: ", match.group(1))
                                res_dict[res].update(match.group(1))


    return



def heatmap(res_dict):
    #x_axis = list(res_dict.keys())
    df = pd.DataFrame(res_dict)
    df = df.sort_index(0, ascending=False)
    df = df.sort_index(1)
    print(df)
    map = sns.heatmap(df, linewidth = 0.5, cmap="YlGnBu")

    import matplotlib.pyplot as plt
    plt.show()

    return


def calc_composition(res_dict):
    tot_interactions = 0
    for key, value in res_dict.items():
        list_of_contacts = list(value.elements())
        print("Residue", key, "has contact with", len(list_of_contacts), "residues")

        tot_interactions = tot_interactions + len(list_of_contacts)

    print("Total number of interactions in this complex: ", tot_interactions)



    return




def main():
    args = sys.argv[1:]
    dataset_file = args[0]
    dataset_file= open(dataset_file, "r")

    complexes = dataset_file.readlines()

    res_dict = defaultdict(Counter)

    check = 0
    for complex in complexes:
        calculate_contacts(complex, res_dict)
        check = check + 1
        calc_composition(res_dict)
        print("File ", check, "complete")
        if check == 10:
            print(res_dict)
            heatmap(res_dict)

            return












    dataset_file.close()







if __name__ == '__main__':
    main()
