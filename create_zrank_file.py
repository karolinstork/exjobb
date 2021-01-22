
import os
import sys

from os import walk




def find_files():
    list_of_files = []
    for (dir_path, dirnames, filenames) in walk("/proj/wallner/users/x_isaak/ModelsKarolin/zrank"):
        for filename in filenames:
            path = dir_path+"/"+filename

            list_of_files.append(path)


    return list_of_files


def main():

    pickle_dict = {}
    zrank_file = open("/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/pro_pep/zranks.txt", "w")

    list_of_files = find_files()

    for file in list_of_files:
        name = os.path.basename(file)
        first_name = name.split(".")[0]

        file = open(file, "r")
        lines = file.readlines()
        for line in lines:
            z_value = line.split()[1]
            id = line.split()[0]
            id = os.path.basename(id)
            last_name = id.split(".")[1]
            full_name = first_name+"_"+last_name
            result = str(full_name)+ '\t'+str(z_value)+ '\n'
            print(result)
            zrank_file.write(result)




    zrank_file.close()











if __name__ == '__main__':
    main()
