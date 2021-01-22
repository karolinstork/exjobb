import os
import sys
import pickle
from collections import defaultdict




def create_dict(cluster_file):
    open_file = open(cluster_file, "r")
    lines = open_file.readlines()

    cluster_dict = defaultdict(list)
    file_dict={}
    i=0

    for line in lines:


        if line[:3] =="pdb": #one row with header info
            continue


        line = line.strip('\n')


        file_and_cluster=line.split('\t')
        file = file_and_cluster[0]
        file_chain = (file.split(".")[0])[:5] #1ix2A
        cluster = int(file_and_cluster[1])
        file_dict[file_chain] = cluster

        print(file_chain)
        print(cluster)

        cluster_dict[cluster].append(file_chain)



    pickle.dump(file_dict, open("pickles/file_cluster_dict.p", "wb"))

    pickle.dump(cluster_dict, open("pickles/cluster_files_dict.p", "wb"))


    return















def main():
    args = sys.argv[1:]
    cluster_file = args[0]

    create_dict(cluster_file)





if __name__ == '__main__':
    main()
