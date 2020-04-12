#!/usr/bin/env python3
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
from collections import defaultdict


import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)


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


def convert_to_fasta(filtered_pdb_file):
    #print(filtered_pdb_file+"--------------------------------------")
    fasta_seq=""
    three_to_one={"ALA":'A',"ARG":'R',"ASN":'N',"ASP":'D',"CYS":'C',"GLU":'E',"GLN":'Q',"GLY":'G',"HIS":'H',"ILE":'I',"LEU":'L',"LYS":'K',"MET":'M',"PHE":'F',"PRO":'P',"SER":'S', "THR":'T',"TRP":'W',"TYR":'Y',"VAL":'V'}

    chain_dict=defaultdict(list)

    file = open(filtered_pdb_file, 'r')
    s = pd.Series(file)

    number_of_models = sum(s.str.contains("MODEL", regex= False))


    if number_of_models>1:
        #Change the chain identifier for NACCESS

        #for model in models_start_end:
        #    endmdl=model[1]
        #    last_row = endmdl-1
        #    print(last_row)
        #    print(atoms.iloc[last_row])


        #print("Multiple models:", number_of_models)

        row_models= s.index[s.str.startswith('MODEL')].tolist()
        row_endmdls=s.index[s.str.startswith('ENDMDL')].tolist()
        models_start_end= list((zip(row_models, row_endmdls)))

        #Keeping only model 1
        row_model1=row_models[0]
        row_endmdl1=row_endmdls[0]
        s=s.truncate(before=row_model1+1, after=row_endmdl1-1) #keeping only model 1


    #Convert chains to FASTA format
    atoms = s.loc[s.str.startswith(('ATOM'), na=False)] #keep only atom records
    res_chain_pos = atoms.str.slice(start=17, stop=26) #keep only residue, chain, position
    c_alphas = res_chain_pos.drop_duplicates() #keep only c alpha

    for c_alpha in c_alphas:
        res = c_alpha[0:3]
        chain = c_alpha[4]
        pos= c_alpha[5:9]
        chain_dict[chain].append((res, pos))

    #print(chain_dict.keys())



    for chain in chain_dict.keys():
        list_of_residues=[]
    #    print(chain)
        for residue_pos in chain_dict[chain]:
            if residue_pos[0] in three_to_one:
                fasta_res = three_to_one[residue_pos[0]]
            else:
                fasta_res = 'X'
            list_of_residues.append(fasta_res)

        if 'X' not in list_of_residues:
            if len(''.join(list_of_residues)) > 6:
                fasta_seq = ">"+filtered_pdb_file[-9:-5]+chain+filtered_pdb_file[-5:]+'\n'+ ''.join(list_of_residues)+'\n'
                print("SUCCESS: File "+filtered_pdb_file[-9:-5]+chain+filtered_pdb_file[-5:]+" was converted to fasta format")
            else:
                print("WARNING: Too short sequence! File "+filtered_pdb_file[-9:-5]+chain+filtered_pdb_file[-5:]+" was removed.")
                fasta_seq=""
        else:
            print("WARNING: Unknown content in sequence. File "+filtered_pdb_file[-9:-5]+chain+filtered_pdb_file[-5:]+" was removed.")
            fasta_seq=""

    return fasta_seq




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



    fasta_file=open("fasta_file.txt", "w")

    for filtered_pdb_file in filtered_pdb_list:
        fasta_seq= convert_to_fasta(filtered_pdb_file)
        fasta_file.write(fasta_seq)



    fasta_file.close()




if __name__ == '__main__':
    main()
