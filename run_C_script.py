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





def find_all_contacts(complex, res_dict, composition_dict):
    translation_dict = {}

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
        complex = complex.strip('\n')
        print("----------------", complex, "------------------")

        output_file = open(f"/proj/wallner/users/x_karst/exjobb/c_script_output/{pdb_id}.output", "w")
        dir = base_file[1:3]
        target = f'/proj/wallner/share/PDB/191015_biounit/{dir}/{base_file}'
        subprocess.run(["/proj/wallner/users/x_bjowa/git/pdb_codes/contact-maps/./contact_map_chain", target, "-chains", chain_1[0], chain_2[0], "-c", "6", "-cb"], stdout = output_file)
        output_file.close()

        output_read = open(f"/proj/wallner/users/x_karst/exjobb/c_script_output/{pdb_id}.output", "r")

        text = output_read.read()
        lines = text.split('\n')

        for line in lines: #translation_dict
            if line[:3] == "RES":
                line_split_by_space = line.split()
                pos = (line_split_by_space[2])
                res = (line_split_by_space[3])
                chain = (line_split_by_space[4]).strip(":")
                pos_chain = pos+chain
                translation_dict[pos_chain] = res



        for line in lines:
            if line[:3] == "RES":
                line_split_by_space = line.split()
                res = (line_split_by_space[3])
                chain = (line_split_by_space[4]).strip(":")
                #print(res, chain, chain_1[:-1], chain_2[:-1])

                pattern = f":(.*)"
                all_contacts_in_string = re.findall(pattern, line)

                if len(all_contacts_in_string)>1:
                    print("ERROR")
                    return

                list_of_contacts = all_contacts_in_string[0].split()
                if len(list_of_contacts)>0:
                    for contact in list_of_contacts:
                        if contact[-1:] != chain:  #only interested in contacts between 2 chains not in the same chain
                            contact_residue = translation_dict[contact]
                            res_dict[res].update(contact_residue)
                            res_dict[contact_residue].update(res)

                            composition_dict[res].update(contact_residue)
                            composition_dict[contact_residue].update(res)
                            #print("Inter-chain contact with ",res ,"found:", contact, "=", contact_residue)


    return







def calc_composition(composition_dict, check):
    tot_interactions = 0
    number_of_contacts_dict = {}
    percentages_dict = {}
    plot_list = []
    x_list = []
    y_list = []

    for key, value in composition_dict.items():
        list_of_contacts = list(value.elements())
        number_of_contacts_dict[key] = len(list_of_contacts)
        tot_interactions = tot_interactions + len(list_of_contacts)

    print("Total number of interactions: ", tot_interactions)
    print("Prevalence of residue in interactions:")
    for residue, interactions in number_of_contacts_dict.items():
        ratio = interactions/tot_interactions
        ratio = (ratio*100)
        ratio = round(ratio, 2)
        print(residue+":", ratio, "%")
        percentages_dict[residue] = ratio
        plot_list.append((residue, ratio))


    plot_list = sorted(plot_list, key = lambda x:x[1], reverse = True)
    for item in plot_list:
        x_list.append(item[0])
        y_list.append(item[1])



    # df = pd.DataFrame(plot_list)
    # df = df.transpose()
    # print(df)
    # percent_plot = sns.scatterplot(data = df)



    import matplotlib.pyplot as plt
    plt.scatter(x_list, y_list)
    plt.plot(x_list, y_list)
    plt.xlabel("Residues")
    plt.ylabel("Wi [%]")
    title = "Residue composition for " + str(check)+ " complexes"
    plt.title(title)
    plt.show()


    return



def plot_heatmap(res_dict, check):
    df = pd.DataFrame(res_dict)
    df = df.sort_index(0, ascending=False)
    df = df.sort_index(1)
    print(df)
    title = "Number of interactions in "+ str(check)+ " complexes"
    heatmap = sns.heatmap(df, linewidth = 0.5, cmap="YlGnBu").set_title(title)

    import matplotlib.pyplot as plt
    plt.show()

    return






def main():
    args = sys.argv[1:]
    dataset_file = args[0]
    dataset_file= open(dataset_file, "r")

    complexes = dataset_file.readlines()

    res_dict = defaultdict(Counter)
    composition_dict = defaultdict(Counter)

    check = 0

    for complex in complexes:
        try:
            find_all_contacts(complex, res_dict, composition_dict)
        except KeyError:
            print("KeyError in ", complex)



        check = check + 1
        if check == 621:
            print(res_dict)
            plot_heatmap(res_dict, check)
            calc_composition(composition_dict, check)

            return





    dataset_file.close()







if __name__ == '__main__':
    main()
