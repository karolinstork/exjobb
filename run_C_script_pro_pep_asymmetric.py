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
from collections import OrderedDict


def find_all_contacts(complex, pro_interface_dict, pep_interface_dict, composition_counter_pro, composition_counter_pep, check, one_model, within_model, between_models, cut_off, no_interchain_contacts_file, which_is_protein_dict):
    translation_dict = {}
    string_position_to_chain_dict = {}
    insertion_residues_list = []
    number_of_close_contacts = 0
    number_of_contacts_between_chains = 0

    line = complex.strip('\n')
    files = line.split('\t')
    file_1 = files[0]
    file_2 = files[1]

    protein = which_is_protein_dict[(file_1, file_2)]
    protein_chain = protein[4]


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
        print("PROTEIN CHAIN = ", protein_chain)

        output_file = open(f"/proj/wallner/users/x_karst/exjobb/c_script_output/{pdb_id}_{chain_1}_{chain_2}_{filetype}.output", "w")
        target = f'/proj/wallner/share/PDB/191015_biounit/{dir}/{base_file}'
        #subprocess.run(["/proj/wallner/users/x_bjowa/git/pdb_codes/contact-maps/./contact_map_chain", target, "-chains", chain_1[0], chain_2[0], "-c", "6", "-cb"], stdout = output_file)
        subprocess.run(["/proj/wallner/users/x_bjowa/git/pdb_codes/contact-maps/./contact_map_chain", target, "-chains", chain_1[0], chain_2[0], "-cb"], stdout = output_file)
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
                    number_of_close_contacts = number_of_close_contacts + 1
                    position_chain_id = string_position_to_chain_dict[column]
                    contact_res = translation_dict[position_chain_id]
                    if position_chain_id[-1] != chain: #only intrested in interactions between two different chains
                        number_of_contacts_between_chains = number_of_contacts_between_chains + 1
                        print(res, pos+chain, "is in contact with", contact_res, "at", position_chain_id, "Distance: ", distance)
                        contact_residue = translation_dict[position_chain_id]

                        if chain==protein_chain: #if contact is from protein side save in protein dict
                            pro_interface_dict[res].update(contact_residue)
                        else:
                            pep_interface_dict[res].update(contact_residue)

                        contact_with_anything_across = True

                column = column + 1

        if contact_with_anything_across == True: #if residue is closer than 6Å to anything across interface, add to counter
            if chain==protein_chain:
                composition_counter_pro.update(res) #composition on protein side
            else:
                composition_counter_pep.update(res) #composition on peptide side




    output_read.close()
    print("Total number of contacts below given cut off ("+ str(cut_off)+" Å):", number_of_close_contacts)
    print("Number of inter-chain contacts:", number_of_contacts_between_chains)

    if number_of_contacts_between_chains == 0:
        mysterious_complex = str(complex)+'\n'
        no_interchain_contacts_file.write(mysterious_complex)
        check = check - 1




    return check, one_model, within_model, between_models


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
        #subprocess.run(["/proj/wallner/users/x_bjowa/git/pdb_codes/contact-maps/./contact_map_chain", target, "-chains", chain_1[0], chain_2[0], "-c", "6", "-cb"], stdout = output_file)
        subprocess.run(["/proj/wallner/users/x_bjowa/git/pdb_codes/contact-maps/./contact_map_chain", target, "-chains", chain_1[0], chain_2[0], "-cb"], stdout = output_file)
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

        temp_file = open(f'/proj/wallner/users/x_karst/exjobb/c_script_temp/{pdb_id}_{chain_1}_{chain_2}.{filetype}', "w")
        for line in model_series:
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


        temp_file.close()
        os.rename(f'/proj/wallner/users/x_karst/exjobb/c_script_temp/{pdb_id}_{chain_1}_{chain_2}.{filetype}', f'/proj/wallner/users/x_karst/exjobb/c_script_temp/{pdb_id}_{chain_1}_{chain_2}_{new_chain_name}.{filetype}' )


        output_file = open(f"/proj/wallner/users/x_karst/exjobb/c_script_output/{pdb_id}_{chain_1}_{chain_2}_{new_chain_name}_{filetype}.output", "w")
        target = f'/proj/wallner/users/x_karst/exjobb/c_script_temp/{pdb_id}_{chain_1}_{chain_2}_{new_chain_name}.{filetype}'
        #subprocess.run(["/proj/wallner/users/x_bjowa/git/pdb_codes/contact-maps/./contact_map_chain", target, "-chains", chain_1[0], new_chain_name, "-c", "6", "-cb"], stdout = output_file)
        subprocess.run(["/proj/wallner/users/x_bjowa/git/pdb_codes/contact-maps/./contact_map_chain", target, "-chains", chain_1[0], new_chain_name, "-cb"], stdout = output_file)
        output_file.close()



    return new_chain_name, within_model, between_models


def next_available(curr_chain_id, orig_chains):
    stepsize= 1
    while curr_chain_id in orig_chains or curr_chain_id.isalpha() == False:
        curr_chain_id= chr(ord(curr_chain_id)+stepsize)

    new_chain_id = curr_chain_id

    return new_chain_id




def calc_composition(composition_counter_pro, composition_counter_pep, check):
    ######################################################################################## protein side
    number_of_contacts_dict_pro = {}
    fraction_dict_pro = {}
    plot_list_pro = []
    x_list_pro = []
    y_list_pro = []


    number_of_contacts_dict_pro = dict(composition_counter_pro) #tex G: 72
    number_of_residues_list_pro = list(composition_counter_pro.elements()) # tex G G och A A A A
    tot_number_involved_residues_pro = len(number_of_residues_list_pro)

    print("Number of residues on protein side of interface:", tot_number_involved_residues_pro)

    for residue, interactions in number_of_contacts_dict_pro.items():
        ratio = interactions/tot_number_involved_residues_pro
        ratio = (ratio*100) #percent form for graph
        plot_list_pro.append((residue, ratio))
        ratio = ratio/100 #decimal form for further calculations
        fraction_dict_pro[residue] = ratio

    plot_list_pro = sorted(plot_list_pro, key = lambda x:x[1], reverse = True)

    print("Composition protein side")
    for (residue, ratio) in plot_list_pro:
        print(str(residue) + ":", ratio)

    for item in plot_list_pro:
        x_list_pro.append(item[0])
        y_list_pro.append(item[1])
    ############################################################################################## peptide side
    number_of_contacts_dict_pep = {}
    fraction_dict_pep = {}
    plot_list_pep = []
    x_list_pep = []
    y_list_pep = []


    number_of_contacts_dict_pep = dict(composition_counter_pep) #tex G: 72
    number_of_residues_list_pep = list(composition_counter_pep.elements()) # tex G G och A A A A
    tot_number_involved_residues_pep = len(number_of_residues_list_pep) #DUBBELKOLLA

    print("Number of residues on peptide side of interface:", tot_number_involved_residues_pep)

    for residue, interactions in number_of_contacts_dict_pep.items():
        ratio = interactions/tot_number_involved_residues_pep
        ratio = (ratio*100) #percent form for graph
        plot_list_pep.append((residue, ratio))
        ratio = ratio/100 #decimal form for further calculations
        fraction_dict_pep[residue] = ratio

    plot_list_pep = sorted(plot_list_pep, key = lambda x:x[1], reverse = True)

    print("Composition peptide side")
    for (residue, ratio) in plot_list_pep:
        print(str(residue) + ":", ratio)

    for item in plot_list_pep:
        x_list_pep.append(item[0])
        y_list_pep.append(item[1])
    ##############################################################################################


    import matplotlib.pyplot as plt

    plt.scatter(x_list_pep, y_list_pep, label = "Peptide side")
    plt.scatter(x_list_pro, y_list_pro, label = "Protein side")


    plt.xlabel("Residues")
    plt.ylabel("Wi [%]")
    plt.ylim(0,13)
    title = "Interface composition for " + str(check)+ " protein peptide complexes"
    plt.title(title)
    legend = plt.legend(loc='best')
    plt.savefig("/proj/wallner/users/x_karst/exjobb/pictures/pro_pep_interface_composition_asymmetric.png")
    plt.show()


    return fraction_dict_pro, fraction_dict_pep



def plot_number_of_interactions(pro_interface_dict, pep_interface_dict, check):
    print('\n')
    print("from protein-side to pep side")
    print(pro_interface_dict)
    print("pep side to pro side")
    print(pep_interface_dict)

    df = pd.DataFrame(pro_interface_dict, pep_interface_dict)
    df = df.sort_index(0, ascending=False)
    df = df.sort_index(1)


    print('\n', "Number of interactions between pairs (x=pro, y=pep)")
    print(df)


    title = "Number of interactions in "+ str(check)+ " complexes"

    heatmap = sns.heatmap(df, linewidth = 0.5, cmap="Blues").set_title(title)


    import matplotlib.pyplot as plt
    plt.xlabel("Protein residues")
    plt.ylabel("Peptide residues")
    plt.savefig("/proj/wallner/users/x_karst/exjobb/pictures/numb_pro_pep_asymmetric.png", bbox_inches='tight')
    plt.show()

    return



def plot_normalized_contacts_and_likelihood(pro_interface_dict, pep_interface_dict, check, fraction_dict_pro, fraction_dict_pep):
    ######################################################################################### Normalized number of contacts Qij
    df = pd.DataFrame(pro_interface_dict, pep_interface_dict)
    df = df.sort_index(0, ascending=False)
    df = df.sort_index(1)


    tot_numb_interactions = df.values.sum()
    tot_pro_interface_interactions = 0
    tot_pep_interface_interactions = 0

    for res, counter in pro_interface_dict.items(): #calc number of contacts from protein side
        all_pro_interactions_list = list(counter.elements()) # tex G G och A A A A
        numb_interactions_pro = len(all_pro_interactions_list)
        tot_pro_interface_interactions = tot_pro_interface_interactions + numb_interactions_pro

    for res, counter in pep_interface_dict.items(): #calc number of contacts peptide side
        all_pep_interactions_list = list(counter.elements()) # tex G G och A A A A
        numb_interactions_pep = len(all_pep_interactions_list)
        tot_pep_interface_interactions = tot_pep_interface_interactions + numb_interactions_pep

    print("Interactions from protein side", tot_pro_interface_interactions)
    print("Interactions from peptide side", tot_pep_interface_interactions)
    print('\n'+ "Total number of interactions:", tot_numb_interactions)

    df = df.div(tot_numb_interactions)
    print('\n'+ "Normalized number of contacts")
    print(df)

    title = "Normalized interactions in "+ str(check)+ " complexes"
    heatmap = sns.heatmap(df, linewidth = 0.5, cmap = "Blues").set_title(title)

    import matplotlib.pyplot as plt
    plt.xlabel("Protein residues")
    plt.ylabel("Peptide residues")
    plt.savefig("/proj/wallner/users/x_karst/exjobb/pictures/norm_numb_pro_pep_asymmetric.png", bbox_inches='tight')
    plt.show()

######################################################################################################### Likelihood Gij
    residue_order = ["R", "K", "N", "D", "Q", "E", "H", "P", "Y", "W", "S", "T", "G", "A", "M", "C", "F", "L", "V", "I"] # to sort the same way the article did
    likelihood_dict_pro = OrderedDict({key: OrderedDict({key: None for key in residue_order}) for key in residue_order})
    likelihood_dict_pep = OrderedDict({key: OrderedDict({key: None for key in residue_order}) for key in residue_order})
    pro_interface_dict = dict(pro_interface_dict)
    pep_interface_dict = dict(pep_interface_dict)


    for master_residue, small_dict in pro_interface_dict.items(): #protein side
        for mini_residue, number in small_dict.items():
            w_i = fraction_dict_pro[master_residue]  #dubbelkolla detta
            w_j = fraction_dict_pep[mini_residue]
            c_ij = pro_interface_dict[master_residue][mini_residue]
            q_ij = c_ij / tot_numb_interactions
            likelihood = 10 * (math.log(q_ij/(w_i * w_j)))
            likelihood_dict_pro[master_residue][mini_residue] = likelihood

    for master_residue, small_dict in pep_interface_dict.items(): #peptide side
        for mini_residue, number in small_dict.items():
            w_i = fraction_dict_pep[master_residue]
            w_j = fraction_dict_pro[mini_residue]
            c_ij = pep_interface_dict[master_residue][mini_residue]
            q_ij = c_ij / tot_numb_interactions
            likelihood = 10 * (math.log(q_ij/(w_i * w_j)))
            likelihood_dict_pep[master_residue][mini_residue] = likelihood


    df_2 = pd.DataFrame(likelihood_dict_pro, likelihood_dict_pep)
    print('\n'+ "Likelihood Gij")
    print(df_2)

    title_2 = "Likelihood of interaction Gij based on "+ str(check)+ " complexes"

    heatmap = sns.heatmap(df_2, linewidth = 0.5, cmap = "seismic", center = 0).set_title(title_2)

    import matplotlib.pyplot as plt
    plt.xlabel("Protein residues")
    plt.ylabel("Peptide residues")
    plt.savefig("/proj/wallner/users/x_karst/exjobb/pictures/gij_pro_pep_asymmetric.png", bbox_inches='tight')
    plt.show()

    ########################################### Likelihood Gij(v) 
    volume_dict = {"I": 166.1, "F":189.7, "V":138.8, "L":168.0 , "W":227.9, "M":165.2, "A":87.8, "G":59.9, "C":105.4, "Y":191.2, "P":123.3, "T":118.3, "S":91.7, "H":156.3, "E":140.9, "N":120.1, "Q":145.1, "D":115.4, "K":172.7, "R":188.2}


    normalized_volume_dict_pro = defaultdict(dict)
    normalized_volume_dict_pep = defaultdict(dict)
    likelihood_volume_dict_pro = OrderedDict({key: OrderedDict({key: None for key in residue_order}) for key in residue_order})
    likelihood_volume_dict_pep = OrderedDict({key: OrderedDict({key: None for key in residue_order}) for key in residue_order})
    sum_of_contacts_times_volume = 0 #sum(Ckl * Vk * Vl)
    onedim_contacts_dict ={}

    for master_residue, small_dict in pro_interface_dict.items():
        for mini_residue, number_of_contacts in small_dict.items():
            onedim_contacts_dict[(master_residue, mini_residue)] = number_of_contacts

    for master_residue, small_dict in pep_interface_dict.items():
        for mini_residue, number_of_contacts in small_dict.items():
            onedim_contacts_dict[(master_residue, mini_residue)] = number_of_contacts

    print(onedim_contacts_dict)
    print('\n'+ "Number of interactions TEST:", sum(onedim_contacts_dict.values()))



    for (res1, res2) in onedim_contacts_dict.keys(): #calc sum(Ckl * Vk * Vl)
        numb_of_contacts_res1_res2 = onedim_contacts_dict[(res1, res2)]
        v_i = volume_dict[res1]
        v_j = volume_dict[res2]
        product = numb_of_contacts_res1_res2 * v_i * v_j
        sum_of_contacts_times_volume = sum_of_contacts_times_volume + product

    print("Sum(Ckl * Vk * Vl):", sum_of_contacts_times_volume)

    for master_residue, small_dict in pro_interface_dict.items(): #calculate new Q(v)ij with regard to volume protein side
        for mini_residue, contacts in small_dict.items():
            c_ij = contacts
            v_i = volume_dict[master_residue]
            v_j = volume_dict[mini_residue]
            norm_vol = c_ij * v_i * v_j /sum_of_contacts_times_volume
            normalized_volume_dict_pro[master_residue][mini_residue] = norm_vol

    for master_residue, small_dict in pep_interface_dict.items(): #calculate new Q(v)ij with regard to volume peptide side
        for mini_residue, contacts in small_dict.items():
            c_ij = contacts
            v_i = volume_dict[master_residue]
            v_j = volume_dict[mini_residue]
            norm_vol = c_ij * v_i * v_j /sum_of_contacts_times_volume
            normalized_volume_dict_pep[master_residue][mini_residue] = norm_vol



    for master_residue, small_dict in pro_interface_dict.items(): #calculate G(v)ij protein side
        for mini_residue, number in small_dict.items():
            w_i = fraction_dict_pro[master_residue]
            w_j = fraction_dict_pep[mini_residue]
            q_v_ij = normalized_volume_dict_pro[master_residue][mini_residue]
            likelihood_vol = 10 * math.log(q_v_ij/(w_i* w_j))
            likelihood_volume_dict_pro[master_residue][mini_residue] = likelihood_vol



    for master_residue, small_dict in pep_interface_dict.items(): #calculate G(v)ij peptide side
        for mini_residue, number in small_dict.items():
            w_i = fraction_dict_pep[master_residue]
            w_j = fraction_dict_pro[mini_residue]
            q_v_ij = normalized_volume_dict_pep[master_residue][mini_residue]
            likelihood_vol = 10 * math.log(q_v_ij/(w_i* w_j))
            likelihood_volume_dict_pep[master_residue][mini_residue] = likelihood_vol



    #print(likelihood_volume_dict)
    df_3 = pd.DataFrame(likelihood_volume_dict_pro, likelihood_volume_dict_pep)

    df_3_reorder = df_3.iloc[::-1]
    df_3_reorder = df_3_reorder.iloc[:, ::-1]
    print('\n'+ "Likelihood Gij(v)")
    print(df_3_reorder)


    df_3_reorder.to_csv('/proj/wallner/users/x_karst/exjobb/protein_peptide_data/gij_v_results_pro_pep_asymmetric.csv', index = False)

    title_3 = "Likelihood of interaction Gij(v) based on "+ str(check)+" complexes"
    heatmap_3 = sns.heatmap(df_3, linewidth = 0.5, cmap = "seismic", center = 0).set_title(title_3)

    import matplotlib.pyplot as plt
    plt.xlabel("Protein residues")
    plt.ylabel("Peptide residues")
    plt.savefig("/proj/wallner/users/x_karst/exjobb/pictures/gijv_pro_pep.png", bbox_inches='tight')
    plt.show()

    print("likelihood protein dict")
    print(likelihood_volume_dict_pro)
    print("likelihood peptide dict")
    print(likelihood_volume_dict_pep)
    gijv_list1 =[]
    gijv_list2 =[]
    diff_dict = defaultdict(dict)

    for res, small_dict in likelihood_volume_dict_pro.items():
        for partner_res in small_dict:
            value1 = small_dict[partner_res]
            print(res, partner_res, value1)
            gijv_list1.append(tuple([res, partner_res, value1]))


    for res, small_dict in likelihood_volume_dict_pep.items():
        for partner_res in small_dict:
            value2 = small_dict[partner_res]
            print(res, partner_res, value2)
            gijv_list2.append(tuple([res, partner_res, value2]))

    for (res11, res12, value1), (res21, res22, value2) in zip(gijv_list1, gijv_list2):
        diff = abs(value1 - value2)
        print(res11, res12, res21, res22)
        diff_dict[res11][res12] = diff

    print(diff_dict)
    df = pd.DataFrame(diff_dict)
    print(df)

    title = "Likelihood differences in two way interactions"
    heatmap = sns.heatmap(df, linewidth = 0.5, cmap = "seismic", center = 0).set_title(title)

    import matplotlib.pyplot as plt
    plt.savefig("/proj/wallner/users/x_karst/exjobb/pictures/same_res_different_sides.png", bbox_inches='tight')
    plt.show()

    return



def main():
    args = sys.argv[1:]
    dataset_file = args[0]



    cut_off=args[1]
    cut_off = float(cut_off)
    dataset_file= open(dataset_file, "r")


    complexes = dataset_file.readlines()

    pro_interface_dict = defaultdict(Counter)
    pep_interface_dict = defaultdict(Counter)

    composition_counter_pro = Counter()
    composition_counter_pep = Counter()


    which_is_protein_dict = {}
    pro_pep_list = open("/proj/wallner/users/x_karst/exjobb/protein_peptide_data/protein_peptide_list.txt", "r")
    list_of_pro_peps = pro_pep_list.readlines()
    for pro_pep in list_of_pro_peps:
        chain1 = pro_pep.split()[0]
        chain2 = pro_pep.split()[1]
        protein = pro_pep.split()[5]
        which_is_protein_dict[(chain1, chain2)] = protein


    one_model = 0
    within_model = 0
    between_models = 0
    error_files = 0
    wrongs_file = open("/proj/wallner/users/x_karst/exjobb/protein_peptide_data/c_script_errors.txt", "w")
    no_interchain_contacts_file = open("/proj/wallner/users/x_karst/exjobb/protein_peptide_data/c_script_no_interchain_contacts_peptide.txt", "w")

    check = 0


    for complex in complexes:
        try:
            check, one_model, within_model, between_models = find_all_contacts(complex, pro_interface_dict, pep_interface_dict, composition_counter_pro, composition_counter_pep, check, one_model, within_model, between_models, cut_off, no_interchain_contacts_file, which_is_protein_dict)
        except KeyError:
            print("KeyError in ", complex)
            error = traceback.format_exc()
            print(error)
            wrongs_file.write(complex)
            wrongs_file.write(error)
            error_files = error_files + 1







    plot_number_of_interactions(pro_interface_dict, pep_interface_dict, check)
    fraction_dict_pro, fraction_dict_pep = calc_composition(composition_counter_pro, composition_counter_pep, check)
    plot_normalized_contacts_and_likelihood(pro_interface_dict, pep_interface_dict, check, fraction_dict_pro, fraction_dict_pep)

    print('\n')
    print("Number of complexes: ", check)
    print("Only one model: ", one_model)
    print("Interactions within same model but more models exist: ", within_model)
    print("Interactions between models: ", between_models )
    print("Error files: ", error_files)



    dataset_file.close()
    wrongs_file.close()
    no_interchain_contacts_file.close()








if __name__ == '__main__':
    main()
