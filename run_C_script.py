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
import math
from matplotlib import colors
import traceback



def find_all_contacts(complex, res_dict, composition_dict, check, one_model, within_model, between_models, number_of_insertion_codes):
    translation_dict = {}
    insertion_residues_list = []


    line = complex.strip('\n')
    files = line.split('\t')
    file_1 = files[0]
    file_2 = files[1]

    file1_parts = file_1.split(".")
    file2_parts = file_2.split(".")
    base_file = (file1_parts[0])[:4]+"."+file1_parts[1]
    filetype = base_file.split(".")[1]
    dir = base_file[1:3]
    pdb_id = (file1_parts[0])[:4]

    chain_1 = (file1_parts[0])[4:]
    chain_2 = (file2_parts[0])[4:]

    only_1_model = False
    numb_of_models = 0
    new_chain_name = None

    pdbfile = open(f"/proj/wallner/share/PDB/191015_biounit/{dir}/{base_file}", "r")
    textlines_pdbfile = pdbfile.readlines()
    for line in textlines_pdbfile:
        if line[:5]== "MODEL":
            numb_of_models = numb_of_models + 1

    #############################################################################################

    if numb_of_models == 1:
        check = check +1
        complex = complex.strip('\n')
        print("----------------", complex, "------------------")

        output_file = open(f"/proj/wallner/users/x_karst/exjobb/c_script_output/{pdb_id}_{chain_1}_{chain_2}_{filetype}.output", "w")
        target = f'/proj/wallner/share/PDB/191015_biounit/{dir}/{base_file}'
        subprocess.run(["/proj/wallner/users/x_bjowa/git/pdb_codes/contact-maps/./contact_map_chain", target, "-chains", chain_1[0], chain_2[0], "-c", "6", "-cb"], stdout = output_file)
        output_file.close()
        one_model = one_model + 1

    if numb_of_models >1:
        complex = complex.strip('\n')
        print("----------------", complex, "------------------")
        new_chain_name, within_model, between_models = many_models_find_contacts(complex, within_model, between_models)
        check = check + 1


    if new_chain_name== None:
        output_read = open(f"/proj/wallner/users/x_karst/exjobb/c_script_output/{pdb_id}_{chain_1}_{chain_2}_{filetype}.output", "r")
    else:
        output_read = open(f"/proj/wallner/users/x_karst/exjobb/c_script_output/{pdb_id}_{chain_1}_{chain_2}_{new_chain_name}_{filetype}.output", "r")


    text = output_read.read()
    lines = text.split('\n')


    for line in lines: #translation_dict
        if line[:3] == "RES":
            line_split_by_space = line.split()
            pos = (line_split_by_space[2]) #should be able to handle negative positions?
            res = (line_split_by_space[3])
            chain = (line_split_by_space[4]).strip(":")
            pos_chain = pos+chain
            if pos[-1].isalpha(): #dont add if insertion code
                print(pos, "Insertion code")
                only_pos_chain = pos[:-1]+chain
                insertion_residues_list.append(only_pos_chain)
            else:
                translation_dict[pos_chain] = res


    for line in lines:
        if line[:3] == "RES":
            line_split_by_space = line.split()
            pos = (line_split_by_space[2])
            res = (line_split_by_space[3])
            chain = (line_split_by_space[4]).strip(":")

            if not pos[-1].isalpha():
                pattern = f":(.*)"
                all_contacts_in_string = re.findall(pattern, line) #list with length 1, all matches as a string
                list_of_contacts = all_contacts_in_string[0].split()

                if len(list_of_contacts)>0:
                    for contact in list_of_contacts:
                        if contact[-1:] != chain:  #only interested in contacts between 2 chains not in the same chain
                            if contact not in insertion_residues_list: #not interested in insertion codes
                                contact_residue = translation_dict[contact]

                                res_dict[res].update(contact_residue)
                                res_dict[contact_residue].update(res)

                                composition_dict[res].update(contact_residue)
                                composition_dict[contact_residue].update(res)
                            #    print("Inter-chain contact with ",res ,"found:", contact, "=", contact_residue)


    output_read.close()
    number_of_insertion_codes = len(insertion_residues_list) + number_of_insertion_codes



    return check, one_model, within_model, between_models, number_of_insertion_codes


def many_models_find_contacts(complex, within_model, between_models):
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
    temp_dict = {}

    number_of_models = sum(s.str.startswith("MODEL"))
    row_models= s.index[s.str.startswith('MODEL')].tolist()
    row_endmdls=s.index[s.str.startswith('ENDMDL')].tolist()
    models_start_end= list((zip(row_models, row_endmdls)))


    model_chain1 = chain_1[1:] #2
    model_chain2 = chain_2[1:]

    new_chain_name = None

    if model_chain1 == model_chain2:
        model_numb = int(model_chain1) - 1
        model_start = (models_start_end[model_numb])[0] #atom records
        model_end = (models_start_end[model_numb])[1] #last row of atom (or hetatm) records
        model_series = s[model_start:model_end]
        temp_file = open(f'/proj/wallner/users/x_karst/exjobb/c_script_temp/{pdb_id}_{chain_1}_{chain_2}.{filetype}', "w")
        for line in model_series:
            temp_file.write(line)
        temp_file.close()

        output_file = open(f"/proj/wallner/users/x_karst/exjobb/c_script_output/{pdb_id}_{chain_1}_{chain_2}_{filetype}.output", "w")
        target = f'/proj/wallner/users/x_karst/exjobb/c_script_temp/{pdb_id}_{chain_1}_{chain_2}.{filetype}'
        subprocess.run(["/proj/wallner/users/x_bjowa/git/pdb_codes/contact-maps/./contact_map_chain", target, "-chains", chain_1[0], chain_2[0], "-c", "6", "-cb"], stdout = output_file)
        output_file.close()
        within_model = within_model + 1

    else:
        #write first model to temp file
        between_models = between_models + 1
        model_numb = int(model_chain1) - 1
        model_start = (models_start_end[model_numb])[0] +1 #atom records
        model_end = (models_start_end[model_numb])[1] #last row of atom (or hetatm) records
        model_series = s[model_start:model_end]
        orig_chains = list(model_series.str.get(21).unique())
        #print(orig_chains)
        temp_file = open(f'/proj/wallner/users/x_karst/exjobb/c_script_temp/{pdb_id}_{chain_1}_{chain_2}.{filetype}', "w")
        for line in model_series:
            #print(line,)
            temp_file.write(line)


        #change chain id:s in second model and write to temp file
        model_numb = int(model_chain2) - 1 #find list index for model
        model_start = (models_start_end[model_numb])[0] +1 #atom records
        model_end = (models_start_end[model_numb])[1] #last row of atom (or hetatm) records
        model_series = s[model_start:model_end]

        for line in model_series:
            curr_chain_id = line[21]
            if curr_chain_id in temp_dict:
                new_chain_id = temp_dict[curr_chain_id]
            else:
                new_chain_id = next_available(curr_chain_id, orig_chains)
                temp_dict[curr_chain_id] = new_chain_id
                orig_chains.append(new_chain_id)

            line= line[:21]+ new_chain_id+line[22:]
            temp_file.write(line)

        new_chain_name = temp_dict[chain_2[0]]
        #print(new_chain_name)

        temp_file.close()
        os.rename(f'/proj/wallner/users/x_karst/exjobb/c_script_temp/{pdb_id}_{chain_1}_{chain_2}.{filetype}', f'/proj/wallner/users/x_karst/exjobb/c_script_temp/{pdb_id}_{chain_1}_{chain_2}_{new_chain_name}.{filetype}' )


        output_file = open(f"/proj/wallner/users/x_karst/exjobb/c_script_output/{pdb_id}_{chain_1}_{chain_2}_{new_chain_name}_{filetype}.output", "w")
        target = f'/proj/wallner/users/x_karst/exjobb/c_script_temp/{pdb_id}_{chain_1}_{chain_2}_{new_chain_name}.{filetype}'
        subprocess.run(["/proj/wallner/users/x_bjowa/git/pdb_codes/contact-maps/./contact_map_chain", target, "-chains", chain_1[0], new_chain_name, "-c", "6", "-cb"], stdout = output_file)
        output_file.close()



    return new_chain_name, within_model, between_models


def next_available(curr_chain_id, orig_chains):
    stepsize= 1
    while curr_chain_id in orig_chains or curr_chain_id.isalpha() == False:
        curr_chain_id= chr(ord(curr_chain_id)+stepsize)

    new_chain_id = curr_chain_id

    return new_chain_id




def calc_composition(composition_dict, check):
    tot_number_involved_residues = 0
    number_of_contacts_dict = {}
    fraction_dict = {}
    plot_list = []
    x_list = []
    y_list = []

    for key, value in composition_dict.items(): #tex G:2, A:4... but in Counter form
        number_of_residues_list = list(value.elements()) # tex G G och A A A A
        number_of_contacts_dict[key] = len(number_of_residues_list) # saves G: 2
        tot_number_involved_residues = tot_number_involved_residues + len(number_of_residues_list)

    print("Total number of residues in interface: ", tot_number_involved_residues)
    print("Prevalence of residue in interactions:")
    for residue, interactions in number_of_contacts_dict.items():
        ratio = interactions/tot_number_involved_residues
        ratio = (ratio*100) #percent form for graph
    #    ratio = round(ratio, 2)
        print(residue+":", ratio, "%")
        plot_list.append((residue, ratio))
        ratio = ratio/100 #decimal form for further calculations
        fraction_dict[residue] = ratio


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
    title = "Interface composition for " + str(check)+ " complexes"
    plt.title(title)
    plt.show()


    return fraction_dict



def plot_number_of_interactions(res_dict, check):
    df = pd.DataFrame(res_dict)
    df = df.sort_index(0, ascending=False)
    df = df.sort_index(1)
    print('\n', "Number of interactions")
    print(df)
    title = "Number of interactions in "+ str(check)+ " complexes"
    heatmap = sns.heatmap(df, linewidth = 0.5, cmap="Blues").set_title(title)

    import matplotlib.pyplot as plt
    plt.show()

    return



def plot_normalized_contacts_and_likelihood(res_dict, check, fraction_dict):
    normalized_dict = defaultdict(dict)
    numb_all_interactions = {}
    res_dict = dict(res_dict)
    tot_numb_interactions = 0

    for master_residue, small_dict in res_dict.items():
        all_interactions_with_master_residue = sum(small_dict.values())
        numb_all_interactions[master_residue] = all_interactions_with_master_residue
        tot_numb_interactions = tot_numb_interactions + all_interactions_with_master_residue
        #print(all_interactions)
    for master_residue, small_dict in res_dict.items():
        for mini_residue, number in small_dict.items():
            normalized_number_contacts = number/tot_numb_interactions
            normalized_dict[master_residue][mini_residue] =normalized_number_contacts

    #print(normalized_dict)
    df = pd.DataFrame(normalized_dict)
    df = df.sort_index(0, ascending=False)
    df = df.sort_index(1)
    print("Normalized number of contacts")
    print(df)

    title = "Normalized interactions in "+ str(check)+ " complexes"
    heatmap = sns.heatmap(df, linewidth = 0.5, cmap = "Blues").set_title(title)

    import matplotlib.pyplot as plt
    plt.show()


    likelihood_dict = defaultdict(dict)

    for master_residue, small_dict in normalized_dict.items():
        for mini_residue, number in small_dict.items():
            w_i = fraction_dict[master_residue]
            w_j = fraction_dict[mini_residue]
            q_ij = normalized_dict[master_residue][mini_residue]

            likelihood = 10 * (math.log(q_ij/(w_i * w_j)))

            likelihood_dict[master_residue][mini_residue] = likelihood

    df_2 = pd.DataFrame(likelihood_dict)
    df_2 = df_2.sort_index(0, ascending=False)
    df_2 = df_2.sort_index(1)
    print("Likelihood")
    print(df_2)

    title_2 = "Likelihood of interaction Gij based on "+ str(check)+ " complexes"

    heatmap = sns.heatmap(df_2, linewidth = 0.5, cmap = "seismic", center = 0).set_title(title_2)
    import matplotlib.pyplot as plt
    plt.show()

    ###########################################
    volume_dict = {"I": 166.1, "F":189.7, "V":138.8, "L":168.0 , "W":227.9, "M":165.2, "A":87.8, "G":59.9, "C":105.4, "Y":191.2, "P":123.3, "T":118.3, "S":91.7, "H":156.3, "E":140.9, "N":120.1, "Q":145.1, "D":115.4, "K":172.7, "R":188.2}
    normalized_volume_dict = defaultdict(dict)
    likelihood_volume_dict = defaultdict(dict)
    sum_of_contacts_times_volume = 0

    for master_residue, small_dict in res_dict.items():
        for mini_residue, number in small_dict.items():
            product = number * volume_dict[master_residue]*volume_dict[mini_residue]
            sum_of_contacts_times_volume = sum_of_contacts_times_volume + product

    for master_residue, small_dict in res_dict.items():
        for mini_residue, number in small_dict.items():
            c_ij = number
            v_i = volume_dict[master_residue]
            v_j = volume_dict[mini_residue]
            norm_vol = c_ij * v_i * v_j /sum_of_contacts_times_volume
            normalized_volume_dict[master_residue][mini_residue] = norm_vol #Q(v)ij

    for master_residue, small_dict in res_dict.items():
        for mini_residue, number in small_dict.items():
            w_i = fraction_dict[master_residue]
            w_j = fraction_dict[mini_residue]
            q_v_ij = normalized_volume_dict[master_residue][mini_residue]
            likelihood_vol = 10 * math.log(q_v_ij/(w_i* w_j))
            likelihood_volume_dict[master_residue][mini_residue] = likelihood_vol


    df_3 = pd.DataFrame(likelihood_volume_dict)
    df_3 = df_3.sort_index(0, ascending = False)
    df_3 = df_3.sort_index(1)

    title_3 = "Likelihood of interaction Gij(v) based on "+ str(check)+" complexes"
    heatmap_3 = sns.heatmap(df_3, linewidth = 0.5, cmap = "seismic", center = 0).set_title(title_3)


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

    number_of_insertion_codes = 0
    one_model = 0
    within_model = 0
    between_models = 0
    error_files = 0
    wrongs_file = open("/proj/wallner/users/x_karst/exjobb/c_script_errors.txt", "w")

    check = 0


    for complex in complexes:
        try:
            check, one_model, within_model, between_models, number_of_insertion_codes = find_all_contacts(complex, res_dict, composition_dict, check, one_model, within_model, between_models, number_of_insertion_codes)
        except KeyError:
            print("KeyError in ", complex)
            error = traceback.format_exc()
            print(error)
            wrongs_file.write(complex)
            wrongs_file.write(error)
            error_files = error_files + 1


        #
        # if check == 1000:
        #     #print(res_dict)
        #     print("Only one model: ", one_model)
        #     print("Interactions within same model but more models exist: ", within_model)
        #     print("Interactions between models: ", between_models )
        #     plot_number_of_interactions(res_dict, check)
        #     fraction_dict = calc_composition(composition_dict, check)
        #     plot_normalized_contacts_and_likelihood(res_dict, check, fraction_dict)
        #     dataset_file.close()
        #
        #     return


    print("Number of complexes: ", check)
    print("Only one model: ", one_model)
    print("Interactions within same model but more models exist: ", within_model)
    print("Interactions between models: ", between_models )
    print("Error files: ", error_files) # 11 st
    print("Number of insertion codes in total: ", number_of_insertion_codes)

    average_insert_codes = number_of_insertion_codes/check
    print("Average number of insertion codes in a file: ", average_insert_codes)

    plot_number_of_interactions(res_dict, check)
    fraction_dict = calc_composition(composition_dict, check)
    plot_normalized_contacts_and_likelihood(res_dict, check, fraction_dict)
    dataset_file.close()


    wrongs_file.close()








if __name__ == '__main__':
    main()
