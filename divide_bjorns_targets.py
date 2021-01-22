import os
import sys
from scipy.stats.stats import pearsonr
import pandas as pd

from os import walk





def main():

    bjornsfile = open("/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/pro_pro/bjorns/bjorns_scores.txt", "r")
    lines = bjornsfile.readlines()
    dirs = []


    for (dir_path, dirnames, filenames) in walk("/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/pro_pro/MOAL/"):
        for dir in dirnames:
            dirs.append(dir)


    open_files_paths = []

    for line in lines:
        file = line.split()[0]
        target = file.split("-")[0]
        if target in dirs:
            target_file = open(f"/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/pro_pro/bjorns/bjorns_{target}", "a")
        #    path = f"/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/pro_pro/bjorns/bjorns_{target}"


            target_file.write(line)





















if __name__ == '__main__':
    main()
