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





def find_all_contacts(complex, res_dict, composition_dict, check):
    translation_dict = {}

    line = complex.strip('\n')
    files = line.split('\t')
    file_1 = files[0]
    file_2 = files[1]

    file1_parts = file_1.split(".")
    file2_parts = file_2.split(".")
    base_file = (file1_parts[0])[:4]+"."+file1_parts[1]
    dir = base_file[1:3]
    pdb_id = (file1_parts[0])[:4]

    chain_1 = (file1_parts[0])[4:]
    chain_2 = (file2_parts[0])[4:]

    only_1_model = False
    numb_of_models = 0

    pdbfile = open(f"/proj/wallner/share/PDB/191015_biounit/{dir}/{base_file}", "r")
    textlines_pdbfile = pdbfile.readlines()
    for line in textlines_pdbfile:
        if line[:5]== "MODEL":
            numb_of_models = numb_of_models + 1

    if numb_of_models == 1:
        check = check +1
        complex = complex.strip('\n')

        print("----------------", complex, "------------------")

        output_file = open(f"/proj/wallner/users/x_karst/exjobb/c_script_output/{pdb_id}_{chain_1}_{chain_2}.output", "w")
        target = f'/proj/wallner/share/PDB/191015_biounit/{dir}/{base_file}'
        subprocess.run(["/proj/wallner/users/x_bjowa/git/pdb_codes/contact-maps/./contact_map_chain", target, "-chains", chain_1[0], chain_2[0], "-c", "6", "-cb"], stdout = output_file)
        output_file.close()


    if numb_of_models >1:
        complex = complex.strip('\n')
        print("----------------", complex, "------------------")
        different_models = many_models_find_contacts(complex)
        if different_models:
            print("Interactions between different models in", complex)
            return check
        check = check + 1


    output_read = open(f"/proj/wallner/users/x_karst/exjobb/c_script_output/{pdb_id}_{chain_1}_{chain_2}.output", "r")

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
                    #    print("Inter-chain contact with ",res ,"found:", contact, "=", contact_residue)





    return check


def many_models_find_contacts(complex):
    different_models = False
    line = complex.strip('\n')
    files = line.split('\t')
    file_1 = files[0]
    file_2 = files[1]

    file1_parts = file_1.split(".")
    file2_parts = file_2.split(".")
    base_file = (file1_parts[0])[:4]+"."+file1_parts[1] #abcd.pdb1
    pdb_id = base_file.split(".")[0] #abcd
    filetype = base_file.split(".")[1] #.pdb1
    dir = base_file[1:3] #bc

    chain_1 = (file1_parts[0])[4:] #A2
    chain_2 = (file2_parts[0])[4:]

    file_path = f"/proj/wallner/share/PDB/191015_biounit/{dir}/{base_file}"
    ##################################################################################################
    file = open(file_path, 'r')
    s = pd.Series(file)

    number_of_models = sum(s.str.startswith("MODEL"))
    row_models= s.index[s.str.startswith('MODEL')].tolist()
    row_endmdls=s.index[s.str.startswith('ENDMDL')].tolist()
    models_start_end= list((zip(row_models, row_endmdls)))

    model_chain1 = chain_1[1:] #2
    model_chain2 = chain_2[1:]

    if model_chain1 == model_chain2:
        model_numb = int(model_chain1) - 1
        model_start = (models_start_end[model_numb])[0] #atom records
        model_end = (models_start_end[model_numb])[1] #last row of atom (or hetatm) records
        model_series = s[model_start:model_end]
        temp_file = open(f'/proj/wallner/users/x_karst/exjobb/c_script_temp/{pdb_id}_{chain_1}_{chain_2}.{filetype}', "w")
        for line in model_series:
            temp_file.write(line)
        temp_file.close()

        output_file = open(f"/proj/wallner/users/x_karst/exjobb/c_script_output/{pdb_id}_{chain_1}_{chain_2}.output", "w")
        target = f'/proj/wallner/users/x_karst/exjobb/c_script_temp/{pdb_id}_{chain_1}_{chain_2}.{filetype}'
        subprocess.run(["/proj/wallner/users/x_bjowa/git/pdb_codes/contact-maps/./contact_map_chain", target, "-chains", chain_1[0], chain_2[0], "-c", "6", "-cb"], stdout = output_file)
        output_file.close()

    else:
        different_models = True



    return different_models




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



    import matplotlib.pyplot as plt
    plt.scatter(x_list, y_list)
    plt.plot(x_list, y_list)
    plt.xlabel("Residues")
    plt.ylabel("Wi [%]")
    plt.ylim(0,12)
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
            check = find_all_contacts(complex, res_dict, composition_dict, check)
        except KeyError:
            print("KeyError in ", complex)


        if check == 100:
            print(res_dict)
            #plot_heatmap(res_dict, check)
            #calc_composition(composition_dict, check)

            return




    dataset_file.close()







if __name__ == '__main__':
    main()
