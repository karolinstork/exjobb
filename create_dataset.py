#!/usr/bin/env python3
#SBATCH -t 7-00:00:00
#SBATCH -N 1

import sys
import os
import pickle
import string
from collections import defaultdict
import operator
import math


def create_dict(complex_area_clusters, resolutions):

    text = complex_area_clusters.readlines()
    dict = defaultdict(list)

    for line in text:
        line = line.strip('\n')
        line = line.split('\t')

        file1 = line[0]
        file2 = line[1]
        area = line[2]
        cluster1 = line[3]
        cluster2 = line[4]

        file1_split = file1.split('.')[0] # in dictionary only 2zzn.pdb, no biounits
        file1_base = file1_split[:4]
        base_file = file1_base+".pdb"

        try:
            resolution = resolutions[base_file]
            if resolution == None: #resolution not applicable in original file, tex NMR
                print("resolution is None in: ", base_file)
                resolution = float(math.inf)
        except KeyError:
            resolution = float(math.inf) #resolutions have been collected from mirabellos directory not 191015 biounit



    #    print(file1, file2, area, cluster1, cluster2, resolution)
        entry = (file1, file2, area, resolution)
        if (cluster2, cluster1) in dict.keys():
            clusters = (cluster2, cluster1)
            dict[clusters].append(entry)
        else:
            clusters = (cluster1, cluster2)
            dict[clusters].append(entry)

    return dict




def find_representative(dict, dataset_file):

    for cluster_pair in dict.keys():
        print("-----------------------------------------")
        print("Family: ", cluster_pair)
        if not "None" in cluster_pair: #ignoring None:s, None:s might be because chains is shorter than 6 aa -> removed from FASTA file

            area_list = []
            dict_90 = defaultdict(list)
            values = dict[cluster_pair]

            for tuple in values:
            #    print(tuple)
                area = float(tuple[2]) #apparently i made the areas in string ..
                area_list.append(area)

            largest_area = max(area_list)
            numb_of_interfaces = len(area_list)
            threshold = 0.9 * largest_area

            print("Number of interfaces: ", numb_of_interfaces)
            print("Largest area: ", largest_area)
            print("Threshold: >=", threshold)

            for tuple in values:
                if float(tuple[2]) >= threshold:
                    dict_90[cluster_pair].append(tuple)


            try:
                representative_complex = min(dict_90[cluster_pair], key=operator.itemgetter(3)) #lowest resolution aka best reoslution
                print("Representative complex: ", representative_complex)
                first_file = representative_complex[0]
                second_file = representative_complex[1]
                representative_complex = str(first_file)+'\t'+str(second_file)+'\n'
                dataset_file.write(representative_complex)
            except TypeError:
                print("ERROR")
                print(dict_90[cluster_pair])



        else:
            print("ERROR: No cluster found for file1 chain and/or file2 chain.")




    return








def main():
    args = sys.argv[1:]
    complex_area_clusters=args[0]

    complex_area_clusters = open(complex_area_clusters, "r")


    resolutions = pickle.load(open("pickles/resolutions_sbatch.p", "rb"))
    dataset_file = open("dataset_file.txt", "w")

    dictionary = create_dict(complex_area_clusters, resolutions)

    find_representative(dictionary, dataset_file)




    dataset_file.close()







if __name__ == '__main__':
    main()
