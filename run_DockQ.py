#!/usr/bin/env python
#SBATCH -t 7-00:00:00
#SBATCH -N 1



import os
import sys
import subprocess
from os import walk
import string
import re




def find_files():
    list_of_files = []
    for (dir_path, dirnames, filenames) in walk("/proj/wallner/users/x_isaak/ModelsKarolin/models/"):
        for filename in filenames:
            path = dir_path+"/"+filename
            list_of_files.append(path)



    print(list_of_files[0])
    print(len(list_of_files))


    return list_of_files



def run_dockq(file):
    first_name = file.split("/")[-2] #tex 4xevCA01
    first_name = first_name[:4] # 4xev
    last_name = file.split("/")[-1] #tex output.2488.pdb
    last_name = last_name.split(".")[1] # 2488

    for (dir_path, dirnames, filenames) in walk("/proj/wallner/users/x_karst/exjobb/renumbered_natives/"): #find the corresponding native file
        for filename in filenames: #list of native files
            if filename[:4] == first_name[:4]:
                native_file = "/proj/wallner/users/x_karst/exjobb/renumbered_natives/"+filename
                print("Model:", file)
                print("Native:", native_file)


    output_file = open(f'/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/pro_pep/dockq_output/{first_name}/{last_name}.dockq', "w")
    output_path = f"/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/pro_pep/dockq_output/{first_name}/{last_name}.dockq"
    subprocess.call(["/proj/wallner/users/x_isaak/DockQ/DockQ.py", file, native_file], stdout = output_file)
    output_file.close()


    return output_path


def save_result(output_path, result_file, file):
    if output_path != None:
        output_file = open(output_path, "r")
        output_end = output_file.readlines()[-1]

        matches = re.findall(r'DockQ[\s]+([\d\.]+)', output_end)

        if len(matches) == 1:
            dockq_value = float(matches[0])


        model_file = file
        model_name = model_file.split("/")[-1]
        native_name = model_file.split("/")[-2]
        native_name = native_name[:4]

        model_id = (native_name + "_"+model_name.split(".")[1])
        full_res = (model_id + '\t'+str(dockq_value)+'\n')
        print(full_res)

        result_file.write(full_res)
        output_file.close()


    return



def main():
    result_file = open("/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/pro_pep/dock_q.txt", "w")
    result_file.write("<MODEL>          <DockQ RESULT>"+'\n')

    list_of_files = find_files() #list of model files

    for file in list_of_files:
        output_path = run_dockq(file) #find corresponding, re-numbered, native file and run DockQ
        save_result(output_path, result_file, file)

















if __name__ == '__main__':
    main()
