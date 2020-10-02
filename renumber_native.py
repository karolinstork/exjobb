import os
import sys
import subprocess
from os import walk
import string
import Bio.PDB

from collections import Counter

import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)




def renumber(file):

    native_path = "/proj/wallner/users/x_karst/exjobb/natives_propep/"+file
    new_native_path = "/proj/wallner/users/x_karst/exjobb/renumbered_natives/"+file
    new_native = open(new_native_path, "w")


    for (dir_path, dirnames, filenames) in walk("/proj/wallner/users/x_isaak/ModelsKarolin/models"): #find a corresponding model
        for filename in filenames:
            dir = dir_path.split("/")[-1]
            if dir[:4] == file[:4]:
                print("")
                print("Model found for native", file)
                print(dir_path, filename)
                model_path = dir_path+"/"+filename
                break


    should_continue = check_similarity(native_path, model_path)

    if should_continue == True:
        try:
            native = open(native_path, "r")
            native_lines = native.readlines()
            model = open(model_path, "r")
            model_lines = model.readlines()
            first_line = model_lines[0]
            start_pos_chain1 = int(first_line[22:27])

            new_chain = False

            for line in model_lines:
                if new_chain == True:
                    start_pos_chain2 = int(line[22:27])
                    new_chain = False
                if line[:3]== "TER":
                    new_chain = True


            print( start_pos_chain1, start_pos_chain2)

            curr_pos = None
            prev_pos = None
            new_chain = True
            chain_numb = 0


            should_break = 0

            for line in native_lines:
                should_break = should_break + 1
                if line[:4] == "ATOM":

                    curr_pos = int(line[24:27])


                    if new_chain == True: #find out start position
                        chain_numb = chain_numb + 1
                        if chain_numb == 1:
                            start_number = str(start_pos_chain1)
                        if chain_numb == 2:
                            start_number = str(start_pos_chain2)

                        new_line = line[:22] + "{:>4}".format(start_number) + line[26:]
                    #    print(new_line)
                        new_chain = False
                    else:
                        if curr_pos != prev_pos:
                            start_number = int(start_number) + 1
                            start_number = str(start_number)
                            new_line = new_line = line[:22] + "{:>4}".format(start_number) + line[26:]
                        #    print(new_line)
                        new_line = new_line = line[:22] + "{:>4}".format(start_number) + line[26:]

                    prev_pos = curr_pos

                if line[:3] == "TER":
                    new_chain = True
                    new_line = line



                new_native.write(new_line)

            native.close()
            model.close()
            new_native.close()

        except:
            print("ERROR in file:", file)
            old_path = "/proj/wallner/users/x_karst/exjobb/renumbered_natives/"+file
            ERROR_path = "/proj/wallner/users/x_karst/exjobb/renumbered_natives/"+file+"_ERROR"
            os.rename(old_path, ERROR_path)

    else:
        print(file, "not corresponding with model")

    return
#########################################################################################

def check_similarity(native_path, model_path):
    #### Native
    pdb_id = "native.pdb"
    parser=Bio.PDB.PDBParser(PERMISSIVE=1)
    structure = parser.get_structure(pdb_id, native_path)

    native_list = []


    for model in structure:
        for chain in model:
            residues = chain.get_residues()
            native_list.append((len(list(residues))))

    #### Model
    pdb_id = "model.pdb"
    parser=Bio.PDB.PDBParser(PERMISSIVE=1)
    structure = parser.get_structure(pdb_id, model_path)

    model_list = []


    for model in structure:
        for chain in model:
            residues = chain.get_residues()
            model_list.append((len(list(residues))))

    print(native_path, model_path)
    print(native_list, model_list)
    if native_list == model_list:
        should_continue = True
    else:
        should_continue = False

    return should_continue



def main():

    list_of_files = []

    for (dir_path, dirnames, filenames) in walk("/proj/wallner/users/x_karst/exjobb/natives_propep"):
        for filename in filenames:
            list_of_files.append(filename)



    for file in list_of_files:
        renumber(file)
    














if __name__ == '__main__':
    main()
