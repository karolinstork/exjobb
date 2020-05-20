import os
import sys
import pickle
from collections import defaultdict




def create_dict(file):
    open_file = open(file, "r")
    lines = open_file.readlines()

    cluster_dict = {}
    file_dict=defaultdict(list)
    i=0

    for line in lines:
        line = line.strip('\n')
        i = i + 1

        print("on row number: ", i)
        files_list = line.split(" ")

        filtered_list = [file for file in files_list if file != ""]
        cluster_dict[i]=filtered_list #clusternumber gives all file in that cluster. whole file name: xxxxA.pdb1

        for file in files_list:
            pdb_id = file.split(".")[0]
            file_dict[pdb_id] = i #file as key gives which row it is on 



    pickle.dump(file_dict, open("pdbfile_clusters_dict.p", "wb"))

    pickle.dump(cluster_dict, open("cluster_files_dict.p", "wb"))


    return















def main():
    args = sys.argv[1:]
    file = args[0]

    create_dict(file)





if __name__ == '__main__':
    main()
