
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
import glob

import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning) #ignores if chains are non continous

def remove_duplicates(line, output_file):
    unique_names = set()


    words = line.split(" ")
    for word in words:
        if len(word)>1:
            #print(word)
            split_file=word.split(".")
            name =split_file[0]
            #print(split_file[0], split_file[1])
            #print(name)
            if name not in unique_names:
                unique_names.add(name)
                name = name + "." + split_file[1]
                output_file.write(name+ " ")
    output_file.write('\n')



    return


def main():
    args = sys.argv[1:]
#    filename = args[0]

    filename=open("blastclust.out", "r")
    output_file=open("blastclust_version2.out", "w")
    text = filename.readlines()


    for line in text:
        remove_duplicates(line, output_file)

    filename.close()
    output_file.close()



if __name__ == '__main__':
        main()
