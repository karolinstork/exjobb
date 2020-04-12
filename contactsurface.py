#!/usr/bin/env python
import sys
import os
import fnmatch
import re
from Bio.PDB import *
import Bio
from operator import itemgetter
import pandas as pd
import subprocess
import string
import itertools


import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning) #ignores if chains are non continous



def find_pdbfiles(dir_name):
    list_of_files=[]
    includes = ['*.pdb?']
    includes = r'|'.join([fnmatch.translate(x) for x in includes]) #from fnmatch to regex

    for root, dirs, files in os.walk(dir_name):
        for file in files:
            match = re.search(includes, file)
            if match:
                list_of_files.append(os.path.join(root,file))

    return list_of_files


def count_models_chains(file):
    pdb_id=file[-8:]
    parser=Bio.PDB.PDBParser(PERMISSIVE=1)
    structure = parser.get_structure(pdb_id, file)

    numb_chains = 0
    numb_models = 0

    for model in structure:
        numb_models = numb_models +1
        for chain in model:
            numb_chains = numb_chains+1

    numb_chains=numb_chains/numb_models
    numb_models_and_chains = (numb_models, numb_chains)

    return numb_models_and_chains


def rename_chains(file_path):

    print(file_path[-9:] + "-----------------------------------")
    file = open(file_path, 'r')
    s = pd.Series(file)

    number_of_models = sum(s.str.startswith("MODEL"))
    row_models= s.index[s.str.startswith('MODEL')].tolist()
    row_endmdls=s.index[s.str.startswith('ENDMDL')].tolist()
    models_start_end= list((zip(row_models, row_endmdls)))

    model1_start = (models_start_end[0])[0]+1
    model1_end = (models_start_end[0])[1]
    model1_series=s[model1_start:model1_end]

    orig_chains = list(model1_series.str.get(21).unique())
    new_chains=orig_chains
    numb_orig_chains = len(orig_chains)
    stepsize=numb_orig_chains

    print("number of models: ", number_of_models)
    print("chains in model1 : ", orig_chains)

    if number_of_models*numb_orig_chains>40:
        print("WARNING: Too many models and/or chains")
        return



    if number_of_models>1:
        i=1
        for i in range(1,len(models_start_end)):
            model_numb=i+1
            #print("MODEL", model_numb)

            model_start = (models_start_end[i])[0]+1 #start of atom records
            model_end = (models_start_end[i])[1] #last row of atom (or hetatm) records
            model_series = s[model_start:model_end]

            for line in model_series: #renaming chains
                curr_chain_id = line[21]

                if curr_chain_id in orig_chains:
                    print("I RECOGNIZE THIS ID! ", orig_chains)
                    new_chain_id = next_available(curr_chain_id, orig_chains)
                    line= line[:21]+ new_chain_id+line[22:]

                    pdbfile = os.path.basename(file_path)
                    filename_split = pdbfile.split('.')
                    path_chainfile = "files_for_naccess/"+filename_split[0]+new_chain_id+'.'+filename_split[1]

                    if not os.path.exists(path_chainfile):
                        chain_file = open(path_chainfile, "w")
                        print("------------START OF FILE >1 models-----------")

                    if line.startswith("ATOM"):# or line.startswith("HETATM"):
                        chain_file.write(line)
                        #print(line.strip('\n'))

                    if line.startswith("TER"):

                        chain_file.close()
                        print("-------------END OF FILE----------------")





    if number_of_models<=1:  #no need to change chain ids

        atom_rows= s.index[s.str.startswith('ATOM')].tolist()
        first_atom_row = atom_rows[0]
        ter_rows=s.index[s.str.startswith('TER')].tolist()
        last_atom_row = int(ter_rows[-1])-1
        model_series = s[first_atom_row:last_atom_row]

        for line in model_series:
            chain_id = line[21]
            pdbfile = os.path.basename(file_path)
            filename_split = pdbfile.split('.')
            path_chainfile = "files_for_naccess/"+filename_split[0]+chain_id+'.'+filename_split[1]

            if not os.path.exists(path_chainfile):
                chain_file = open(path_chainfile, "w")
                print("------------START OF FILE-----------")

            if line.startswith("ATOM"):# or line.startswith("HETATM"):
                chain_file.write(line)
                #print(line.strip('\n'))

            if line.startswith("TER"):
                chain_file.close()
                print("------------END OF FILE-----------")



        #    for chain in chains_in_one_model:
        #        print("CHAIN "+ chain)
        #        for line in model_i_series:
        #            if line[21] == chain:
        #                print(line.strip('\n'))

        #chain_ids = slized_s.str.get(21)
        #print(chain_ids)

        #duplicates = slized_s.duplicated(keep='first')
        #print(s.loc[duplicates])

    #    for line in s:
    #        curr_res = line[17:19]
    #        if line.startswith("ATOM") or line.startswith("HETATM"):
    #            if line[21] not in chain_ids:
    #                orig_chains.append(line[21])
    #            else:
    #                temp_chains.append(i)
    ##
    #        prev_res = curr_res
    ##            i=i+1


            #if not (line.startswith("HETATM") and line[17:19]=="HOH"):
        #    print(line.strip('\n'))

            #else:
                #continue


            #renamed_chains_file.write(line)

        #row_models= s.index[s.str.startswith('MODEL')].tolist()
        #row_endmdls=s.index[s.str.startswith('ENDMDL')].tolist()
        #models_start_end= list((zip(row_models, row_endmdls)))
        #print(number_of_models," models in total")

        #print(models_start_end)
        #iterator = len(models_start_end) -1

        #for i in range(iterator):
        #    last_atomrecord_row = (models_start_end[i])[1] -1
        #    line = s.iloc[last_atomrecord_row]
        #    print(line.strip('\n'))
        #    last_chainid = line[21]
        #    print("Last chain id found was ",last_chainid," on row " ,last_atomrecord_row)


            #print("changing chain ID in model", i+2)

        #    nextmodel_start=(((models_start_end[i])[0])+1)
        #    nextmodel_end=(((models_start_end[i])[1])-1)
        #    next_model =s[nextmodel_start:nextmodel_end]

        #    for row in next_model:
        #        if row[21] == last_chainid:
        #            row[21] == last_chainid


            #NUMBER OF MODELS = 2
        #    if number_of_models==100:
            #    for line in model_series: #renaming chains
            #        if line[21] in orig_chains:
            #            curr_chain_id = line[21]
            #            new_chain_id = lower_upper(curr_chain_id)
#
#                        pdbfile = os.path.basename(file_path)
#                        filename_split = pdbfile.split('.')
#                        pdbfile_chain_name = "files_for_naccess/"+filename_split[0]+new_chain_id+'.'+filename_split[1]
#
#                        if not os.path.exists(pdbfile_chain_name):
#                            pdb_chain_file = open(pdbfile_chain_name, "w")
#                            line= line[:21]+ new_chain_id+line[22:]
#
#                        if line.startswith("ATOM") or line.startswith("HETATM"):
#                            pdb_chain_file.write(line.strip('\n'))

#                        if line.startswith("TER"):
#                            pdb_chain_file.close()


    return



def lower_upper(curr_chain_id):
    if curr_chain_id.lower():
        new_chain_id = curr_chain_id.upper()

    if curr_chain_id.upper():
        new_chain_id = curr_chain_id.lower()

    if not curr_chain_id.isalpha():
        print("ERROR: chain id is not a letter")

    return new_chain_id


def next_available(curr_chain_id, orig_chains):
    stepsize= 1

    while curr_chain_id in orig_chains:
        curr_chain_id= chr(ord(curr_chain_id)+stepsize)
        stepsize = stepsize + 1

    new_chain_id = curr_chain_id
    orig_chains.append(new_chain_id)
    #stepsize = stepsize + numb_orig_chains

    return new_chain_id



def calc_interaction_area():
    list_of_files=os.listdir("/proj/wallner/users/x_karst/files_for_naccess")
    combinations_list= itertools.combinations(list_of_files, 2)
    print(list(combination_of_files))

    for combination_tuple in combinations_list:
        largest_area=0
        individual_chain_area = []
        for file in combination_tuple:
            #subprocess.run(["/proj/wallner/users/x_bjowa/local/naccess/./naccess", "/proj/wallner/user/x_karst/files_for_naccess/file"])
            #output_file = file[:-4]+".rsa"
            #output_file = open(output_file, "r")
            #text = output_file.read()
            #chain_area = re.findall("TOTAL[\s]+([\d]+)", text)
            #individual_chain_area.append(chain_area)
            #remove outputfile
            largest_area=0

        file1 = open(combination[0], "r")
        file2 = open(combination[1], "r")
        file3 = open("binary_file.pdb", "w")

        file1_text = file1.read()
        file2_text = file2.read()
        file3_text= file1_text+file2_text
        file3.write(file3_text)

        file1.close()
        file2.close()
        file3.close()


        #subprocess.run(["/proj/wallner/users/x_bjowa/local/naccess/./naccess", "/proj/wallner/user/x_karst/files_for_naccess/binary_file.pdb"])
        #output_file = binary_file.pdb.rsa"
        #output_file = open(output_file, "r")
        #text = output_file.read()
        #chain_area = re.findall("", text)
        #individual_chain_area.append(chain_area)

        #calculate difference

        #save difference with the right chain name



    return





def main():
    args = sys.argv[1:]
    dir_name = args[0]

    list_of_files = find_pdbfiles(dir_name)
    filtered_pdb_list=[]


    for file in list_of_files:
        try:
            numb_models_and_chains = count_models_chains(file)
            if numb_models_and_chains != (1,1):
                filtered_pdb_list.append(file)
        except Exception as error:
            print(error)
            print(" in ")
            print(file)



    for file_path in filtered_pdb_list:
    #file_path = "/proj/wallner/share/PDB/191015_biounit/mq/1mqn.pdb4"
        rename_chains(file_path)
        #calc_interaction_area()





if __name__ == '__main__':
    main()
