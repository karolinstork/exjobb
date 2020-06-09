import sys
import os
import pickle
import re





def edit(dict, naccess_file, complex_area_clusters):
    lines = naccess_file.readlines()
    numb_of_missing = 0
    numb_of_weirdos = 0
    nones = 0

    for line in lines:
        print("---------------------------")
        file_file_area = line.split()

        file1_full = file_file_area[0] #3go6A1.pdb1
        file1 = file1_full.split(".")[0]
        file1 = file1_full[:5] #3go6A

        file2_full = file_file_area[1]
        file2 = file2_full.split(".")[0]
        file2 = file2[:5] #3go6A

        area = file_file_area[2]



        try:
            cluster1 = dict[file1] #3go6A1.pdb1
        except KeyError:
            cluster1= None
            nones = nones + 1





        try:
            cluster2 = dict[file2]
        except KeyError:
            cluster2= None
            nones = nones +1




        string = str(file1_full)+'\t'+str(file2_full)+'\t'+str(area)+'\t'+ str(cluster1)+'\t'+str(cluster2)+'\n'
        print(string)


        complex_area_clusters.write(string)


    print("number of nones: ", nones)





    return




def main():
    args = sys.argv[1:]
    naccess_file = args[0]

    dict = pickle.load(open("pickles/file_cluster_dict.p", "rb"))
    naccess_file = open(naccess_file, "r")
    complex_area_clusters = open("complex_area_clusters.txt", "w")

    edit(dict, naccess_file, complex_area_clusters)

    complex_area_clusters.close()









if __name__ == '__main__':
    main()
