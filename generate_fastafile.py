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
    pdb_id=file[-8:] #temp name
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


def convert_to_fasta(pdb_file, fasta_file, entry_list):
    #print(filtered_pdb_file+"--------------------------------------")
    fasta_seq=""
    three_to_one={"ALA":'A',"ARG":'R',"ASN":'N',"ASP":'D',"CYS":'C',"GLU":'E',"GLN":'Q',"GLY":'G',"HIS":'H',"ILE":'I',"LEU":'L',"LYS":'K',"MET":'M',"PHE":'F',"PRO":'P',"SER":'S', "THR":'T',"TRP":'W',"TYR":'Y',"VAL":'V'}
    file = open(pdb_file, 'r')
    s = pd.Series(file)

    number_of_models = sum(s.str.contains("MODEL", regex= False))

    row_models= s.index[s.str.startswith('MODEL')].tolist()
    row_endmdls=s.index[s.str.startswith('ENDMDL')].tolist()
    models_start_end= list((zip(row_models, row_endmdls)))

    i = 0
    for model in models_start_end:
        chain_dict=defaultdict(list)
        model_numb = i + 1
    #    print("MODEL", model_numb)
        row_model=row_models[i]
        row_endmdl=row_endmdls[i]

        model=s.truncate(before=row_model, after=row_endmdl)

        #Slicing out ATOM records from PDB file
        atoms = model.loc[model.str.startswith(('ATOM'), na=False)] #keep only atom records
        res_chain_pos = atoms.str.slice(start=17, stop=26) #keep only residue, chain, position
        c_alphas = res_chain_pos.drop_duplicates() #keep only c alpha

        #Building dictionary where chain is key and (res,pos) is value
        for c_alpha in c_alphas:
            res = c_alpha[0:3]
            chain = c_alpha[4]+str(model_numb)
            pos= c_alpha[5:9]
            chain_dict[chain].append((res, pos))


        for chain in chain_dict.keys(): #translate each chain to FASTA
            entry = pdb_file[-9:-5]+chain[:-1] #only write if not already in file

            if entry not in entry_list:
                entry_list.append(entry)
                fasta_seq=""
                list_of_residues=[]

                for residue_pos in chain_dict[chain]:
                    if residue_pos[0] in three_to_one:
                        fasta_res = three_to_one[residue_pos[0]] #check if known residue ie amino acids
                    else:
                        fasta_res = 'X'

                    list_of_residues.append(fasta_res)


                if 'X' not in list_of_residues:
                    if len(''.join(list_of_residues)) > 6:
                        fasta_seq = ">"+pdb_file[-9:-5]+chain+pdb_file[-5:]+'\n'+ ''.join(list_of_residues)+'\n'
                        print("Converting "+ chain+ " ...")

                    else:
                        print("WARNING: Too short sequence! File "+pdb_file[-9:-5]+chain+pdb_file[-5:]+" was removed.")
                        fasta_seq=""
                else:
                        print("WARNING: Unknown content in sequence. File "+pdb_file[-9:-5]+chain+str(model_numb)+pdb_file[-5:]+" was removed.")
                        fasta_seq=""
                print(fasta_seq)
                #fasta_file.write(fasta_seq)


            else:
                print("WARNING:", entry, "already exists in FASTA file!", pdb_file[-9:-5]+chain, "will not be added." )
        i = i + 1

    return




def main():
    args = sys.argv[1:]
    dir_name = args[0]

    list_of_files = find_pdbfiles(dir_name)


    fasta_file=open("fasta_file_test.txt", "w")
    entry_list = []
    x=0

    for file in list_of_files:
        x = x+1
        try:
            numb_models_and_chains = count_models_chains(file)

            if numb_models_and_chains != (1,1.0):
                print("----------------------------------------------------")
                print(file, "Number of models and chains:" ,numb_models_and_chains)
                convert_to_fasta(file, fasta_file, entry_list) #måste entry list skickas ut igen eller spars allt i dennu tills nästa varv i loopen?

        except Exception as error:
            print(error)
            print(file)

        if x == 20 :
            break





    fasta_file.close()




if __name__ == '__main__':
    main()
