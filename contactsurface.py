#!/usr/bin/env python
#SBATCH -t 7-00:00:00
#SBATCH -N 1

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
import glob


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

def create_path_renamed_chain(file_path, curr_chain_id, curr_model, new_chain_id):
    pdbfile = os.path.basename(file_path)
    filename_split = pdbfile.split('.')
    path_chainfile = f'files_for_naccess/{filename_split[0]}{curr_chain_id}{curr_model}_{new_chain_id}.{filename_split[1]}'

    return path_chainfile

def create_path_chain(file_path, chain_id, curr_model):
    pdbfile = os.path.basename(file_path)
    filename_split = pdbfile.split('.')
    path_chainfile = f'files_for_naccess/{filename_split[0]}{chain_id}{curr_model}.{filename_split[1]}'

    return path_chainfile


def rename_xchains(file_path):
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


    if number_of_models<=1:
        for line in model1_series:
            chain_id = line[21]
            path_chainfile = create_path_chain(file_path, chain_id, curr_model)

            if not os.path.exists(path_chainfile) and not line.startswith("ANISOU"):
                chain_file = open(path_chainfile, "w")
                open_files_dict[chain_id] = chain_file
                opened_files.append(chain_file)
                print("Creating file "+ path_chainfile)

            if line.startswith("ATOM") or line.startswith("HETATM"):
                open_files_dict[chain_id].write(line)
            #    print(line.strip('\n'))
        for file in opened_files:
            file.close()

    if number_of_models>1:
        text = s.to_string()
        print(text)
        chains = re.findall("TER[\s]+[\d.]+[\s]+[\w]+[\s]+([\w])",text)




    return

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
    numb_orig_chains = len(orig_chains)

    previous_chain = None
    new_model = None
    prev_model = None
    new_chain_id= None
    #stepsize=numb_orig_chains
    open_files_dict={}
    temp_dict={}

    opened_files=[]

    print("Number of models: ", number_of_models)
    print("Chains in model 1: ", orig_chains)

    if number_of_models*numb_orig_chains>40:
        print("WARNING: Too many models and/or chains")
        return

    if number_of_models>1:

        curr_model = 1
        for line in model1_series: #model1 no need to change chain ids
            chain_id = line[21]
            path_chainfile = create_path_chain(file_path, chain_id, curr_model)

            if not os.path.exists(path_chainfile) and not line.startswith("ANISOU"):
                chain_file = open(path_chainfile, "w")
                open_files_dict[chain_id] = chain_file
                opened_files.append(chain_file)
                print("Creating file "+ path_chainfile)

            if line.startswith("ATOM") or line.startswith("HETATM"):
                open_files_dict[chain_id].write(line)
                #print(line.strip('\n'))

        for file in opened_files:
            file.close()


        for i in range(1,len(models_start_end)): #rest of models, change chain ids
            new_model = True
            curr_model=i+1

            print("Changing the chain id:s in MODEL", curr_model, ". . .")

            model_start = (models_start_end[i])[0]+1 #atom records
            model_end = (models_start_end[i])[1] #last row of atom (or hetatm) records
            model_series = s[model_start:model_end]

            for line in model_series: #renaming chains
                curr_chain_id = line[21]
                if curr_chain_id in temp_dict:
                    new_chain_id = temp_dict[curr_chain_id]
                else:
                    new_chain_id = next_available(curr_chain_id, orig_chains)
                    temp_dict[curr_chain_id] = new_chain_id
                    orig_chains.append(new_chain_id)
                    path_chainfile = create_path_renamed_chain(file_path, curr_chain_id, curr_model, new_chain_id)

                line= line[:21]+ new_chain_id+line[22:]
                chain_id = new_chain_id

                if not os.path.exists(path_chainfile) and not line.startswith("ANISOU"):
                    chain_file = open(path_chainfile, "w")
                    open_files_dict[chain_id] = chain_file
                    opened_files.append(chain_file)
                    print("Creating file "+path_chainfile)

                if line.startswith("ATOM") or line.startswith("HETATM"):
                    chain_id = line[21]
                    open_files_dict[chain_id].write(line)
                    #print(line.strip('\n'))

                previous_chain = curr_chain_id
                #new_model = False

            #    else: #some files with multiple models already have unique chain ids
            #        chain_id=line[21]

                    #line= line[:21]+ new_chain_id+line[22:]

                #    path_chainfile = create_path_chain(file_path, chain_id, curr_model)

                    # if not os.path.exists(path_chainfile) and not line.startswith("ANISOU"):
                    #     chain_file = open(path_chainfile, "w")
                    #     open_files_dict[chain_id] = chain_file
                    #     opened_files.append(chain_file)
                    #     print("Creating ~special~ file "+path_chainfile)
                    #
                    # if line.startswith("ATOM") or line.startswith("HETATM"):
                    #     chain_id = line[21]
                    #     open_files_dict[chain_id].write(line)
                    #     #print(line.strip('\n'))
                    # if line.startswith("TER"):
                    #     chain =line[21]
                    #     orig_chains.append(chain)


            print(orig_chains)
            temp_dict = {}
            prev_model = curr_model
            for file in opened_files:
                file.close()



    if number_of_models<=1:  #no need to change chain ids
        curr_model = 1

        atom_rows= s.index[s.str.startswith('ATOM')].tolist()
        first_atom_row = atom_rows[0]
        ter_rows=s.index[s.str.startswith('TER')].tolist()
        last_atom_row = int(ter_rows[-1])-1
        model_series = s[first_atom_row:last_atom_row]

        for line in model_series:
            chain_id = line[21]
            path_chainfile = create_path_chain(file_path, chain_id, curr_model)

            if not os.path.exists(path_chainfile) and not line.startswith("ANISOU"):
                chain_file = open(path_chainfile, "w")
                open_files_dict[chain_id] = chain_file
                opened_files.append(chain_file)
                print("Creating file "+ path_chainfile)

            if line.startswith("ATOM") or line.startswith("HETATM"):
                open_files_dict[chain_id].write(line)
            #    print(line.strip('\n'))

        for file in opened_files:
            file.close()

    return




def next_available(curr_chain_id, orig_chains):
    stepsize= 1
    while curr_chain_id in orig_chains:
        curr_chain_id= chr(ord(curr_chain_id)+stepsize)
    new_chain_id = curr_chain_id

    return new_chain_id


def get_orig_chain_id(filename):
    if "_" in filename:
        seperation_by_underscore = filename.split("_")
        separation_by_dot = seperation_by_underscore[1].split(".")
        orig_chain_id = seperation_by_underscore[0]+"."+separation_by_dot[1]
    else:
        orig_chain_id = filename

    return orig_chain_id


def run_naccess(result_naccess_file):
    list_of_files=os.listdir("/proj/wallner/users/x_karst/exjobb/files_for_naccess")

    if len(list_of_files)==0:
        return

    combination_of_files= itertools.combinations(list_of_files, 2)
    combinations_list = list(combination_of_files)
    print("Doing binary comparisons . . . ")

    for combination_tuple in combinations_list:
        print(combination_tuple)
        individual_chain_area = []
        for file in combination_tuple: #each chain area for themselves
            target_file=f'/proj/wallner/users/x_karst/exjobb/files_for_naccess/{file}'
            try:
                subprocess.run(["/proj/wallner/users/x_bjowa/local/naccess/./naccess", target_file], timeout = 10)
            except subprocess.TimeoutExpired:
                print("Time limit exceeded for ", file)
                files_to_be_removed = glob.glob("/proj/wallner/users/x_karst/exjobb/files_for_naccess/*")
                for f in files_to_be_removed:
                    os.remove(f)
                    

                return


            filename_split = file.split('.')
            file=filename_split[0]
            rsa_file = f'{file}.rsa'

            rsa_file = open(rsa_file, "r")
            text = rsa_file.read()
            chain_area = re.findall("TOTAL[\s]+([\d.]+)", text)

            chain_area = float(chain_area[0])
            individual_chain_area.append(chain_area)

            os.remove(file+".rsa")
            os.remove(file+".log")
            os.remove(file+".asa")

        sum = individual_chain_area[0] + individual_chain_area[1]

        path = "/proj/wallner/users/x_karst/exjobb/files_for_naccess/"
        file1 = open(path+combination_tuple[0], "r")
        file2 = open(path+combination_tuple[1], "r")
        binary_file = open("binary_file.pdb", "w")

        file1_text = file1.read()
        file2_text = file2.read()
        text= file1_text+file2_text
        binary_file.write(text)

        file1.close()
        file2.close()
        binary_file.close()

        try:
            subprocess.run(["/proj/wallner/users/x_bjowa/local/naccess/./naccess", "/proj/wallner/users/x_karst/exjobb/binary_file.pdb"], timeout = 10)
        except subprocess.TimeoutExpired:
            print("Time limit exceeded for binary_file.pdb")
            files_to_be_removed = glob.glob("/proj/wallner/users/x_karst/exjobb/files_for_naccess/*")
            for f in files_to_be_removed:
                os.remove(f)
            os.remove(file+".rsa")
            os.remove(file+".log")
            os.remove(file+".asa")

            return


        output_file = "binary_file.rsa"
        output_file = open(output_file, "r")
        text = output_file.read()

        complex_areas = re.findall("CHAIN[\s]+[\d]+[\s]+[\w]+[\s]+([\d.]+)", text)
        complex_area1 = float(complex_areas[0])
        complex_area2 = float(complex_areas[1])

        size_contactarea=(sum-complex_area1-complex_area2)/2
        size_contactarea = round(size_contactarea, 4)

        if size_contactarea>0:
            input_file1 = combination_tuple[0]
            input_file2 = combination_tuple[1]

            input_file1 = get_orig_chain_id(input_file1)
            input_file2 = get_orig_chain_id(input_file2)

            result = f'{input_file1}\t{input_file2}\t{size_contactarea}\n'
            result_naccess_file.write(result)

        os.remove("binary_file.rsa")
        os.remove("binary_file.log")
        os.remove("binary_file.asa")
        os.remove("binary_file.pdb")

    files_to_be_removed = glob.glob("/proj/wallner/users/x_karst/exjobb/files_for_naccess/*")

    for f in files_to_be_removed:
        os.remove(f)

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


    result_naccess_file = open("/proj/wallner/users/x_karst/exjobb/result_naccess.txt", "w")


    for file_path in filtered_pdb_list:
        rename_chains(file_path)
        run_naccess(result_naccess_file)



    result_naccess_file.close()


if __name__ == '__main__':
    main()
