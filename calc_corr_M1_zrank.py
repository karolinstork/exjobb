import os
import sys
from scipy.stats.stats import pearsonr
import pandas as pd







def main():
    f = open("/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/pro_pro/bjorns/bjorns_scores.txt", "r")
    lines = f.readlines()
    lines = lines[1:]


    zrank_list = []
    dockq_list = []



    for line in lines:
        print(line)

        zrank = line.split()[7]
        dockq = line.split()[16]

        print(zrank, dockq)
        zrank = float(zrank)
        dockq = float(dockq)


        zrank_list.append(zrank)
        dockq_list.append(dockq)



    correlation, p_value = pearsonr(dockq_list, zrank_list)
    print("Correlation:", correlation)

















if __name__ == '__main__':
    main()
