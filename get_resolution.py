#!/usr/bin/env python3
#SBATCH -t 7-00:00:00
#SBATCH -N 1



import sys
import os
import fnmatch
import re
from operator import itemgetter
from collections import defaultdict
import itertools
import string
from Bio.PDB import *
import Bio
import glob
import pickle

import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning) #ignores if chains are non continous

def main():
    args = sys.argv[1:]
    path = args[0]

    pickle_dict={}
    all_files =[]

    # all_files = glob.glob(path+"/*")
    # dir = path[-2:]

    list_of_dirs = glob.glob(path+"/*")

    for dir in list_of_dirs:
        files_in_dir = glob.glob(dir+"/*")
        for file in files_in_dir:
            all_files.append(file)

    print("looking for resolutions . . . ")


    for file in all_files:
        try:
            pdb_id = os.path.basename(os.path.normpath(file))
            if pdb_id[-3:] == "pdb":
                parser=Bio.PDB.PDBParser(PERMISSIVE=1)
                structure = parser.get_structure(pdb_id, file)
                resolution = structure.header["resolution"]
                pickle_dict[pdb_id]=resolution
                print(pdb_id, resolution)

        except Exception as error:
            print("ERROR!", file, "caused: ")
            print(error)





    pickle.dump(pickle_dict, open("resolutions_sbatch.p", "wb"))


if __name__ == '__main__':
    main()
