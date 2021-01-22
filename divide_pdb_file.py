#!/usr/bin/env python
#SBATCH -t 7-00:00:00
#SBATCH -N 1

import sys
import os
import fnmatch
import re
from operator import itemgetter
import pandas as pd
import subprocess
import string





def rename_chains(complex):
    line = complex.strip('\n')
    files = line.split('\t')
    file_1 = files[0]
    file_2 = files[1]

    file1_parts = file_1.split(".")
    file2_parts = file_2.split(".")
    base_file = (file1_parts[0])[:4]+"."+file1_parts[1] #abcd.pdb1
    dir = base_file[1:3]

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

    model_chain1 = chain_1[-1]
    model_chain2 = chain_2[-1] 




def main():

    args = sys.argv[1:]
    file = args[0]

    rename_chains(complex)






if __name__ == '__main__':
    main()
