import os
import sys
from os import walk
import pickle
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import csv







def find_targets():

    targets_list = []
    for (dir_path, dirnames, filenames) in walk(f"/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/propro_matrix/proteinpeptide_scored/"):
        for file in filenames:
            target_path = dir_path+file
            targets_list.append(target_path)


    for entry in targets_list:
        print(entry)

    return targets_list #find all scoring files for targets.



def sort_scores_023(targets_list):
    dict_023_performance = {}

    correct_models_dict = pickle.load(open("/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/propep_matrices/dockq_natives_0.23.p", "rb"))


    for target in targets_list:
        print("------------------------------------------")
        print(target, "Benchmark: 0.23")

        target_file = target
        target_name = os.path.basename(target)

        score_file = open(target_file, "r")
        score_file_lines = score_file.readlines()
        sorted_lines = sorted(score_file_lines, key=lambda line: float(line.split()[1]), reverse = True)

        try:
            top1 = sorted_lines[0]
            top10 = sorted_lines[:10]
            top100 = sorted_lines[:100]
        except IndexError:
            print("ERROR in", target, "File is empty due to no close contacts <6Å")
            continue


        performance_1 = check_top1(top1, correct_models_dict, )
        performance_10 = check_top10(top10, correct_models_dict, )
        performance_100 = check_top100(top100, correct_models_dict, )


        dict_023_performance[target_name] = (performance_1, performance_10, performance_100)



    return dict_023_performance



def sort_scores_049(targets_list):
    dict_049_performance = {}
    correct_models_dict = pickle.load(open("/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/propep_matrices/dockq_natives_0.49.p", "rb"))


    for target in targets_list:
        print("------------------------------------------")
        print(target, "Benchmark: 0.49")

        target_file = target
        target_name = os.path.basename(target)


        score_file = open(target_file, "r")
        score_file_lines = score_file.readlines()
        sorted_lines = sorted(score_file_lines, key=lambda line: float(line.split()[1]), reverse = True)

        try:
            top1 = sorted_lines[0]
            top10 = sorted_lines[:10]
            top100 = sorted_lines[:100]
        except IndexError:
            print("ERROR in", target)
            continue


        performance_1 = check_top1(top1, correct_models_dict, )
        performance_10 = check_top10(top10, correct_models_dict, )
        performance_100 = check_top100(top100, correct_models_dict, )


        dict_049_performance[target_name] = (performance_1, performance_10, performance_100)

    return dict_049_performance




def sort_scores_080(targets_list):
    dict_080_performance = {}


    correct_models_dict = pickle.load(open("/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/propep_matrices/dockq_natives_0.8.p", "rb"))


    for target in targets_list:
        print("------------------------------------------")
        print(target, "Benchmark: 0.80")

        target_file = target
        target_name = os.path.basename(target)

        score_file = open(target_file, "r")
        score_file_lines = score_file.readlines()
        sorted_lines = sorted(score_file_lines, key=lambda line: float(line.split()[1]), reverse = True)

        try:
            top1 = sorted_lines[0]
            top10 = sorted_lines[:10]
            top100 = sorted_lines[:100]
        except IndexError:
            print("ERROR in", target)
            continue


        performance_1 = check_top1(top1, correct_models_dict, )
        performance_10 = check_top10(top10, correct_models_dict, )
        performance_100 = check_top100(top100, correct_models_dict, )


        dict_080_performance[target_name] = (performance_1, performance_10, performance_100)

    return dict_080_performance



def check_top100(top100, correct_models_dict):
    correct_models = 0
    total_models = 0

    for entry in top100:
        model = entry.split()[0]

        if model in correct_models_dict.keys():
            correct_models = correct_models + 1
        total_models = total_models + 1


    performance_100 = (float(correct_models)/total_models) * 100
    print("Number of correct models:", correct_models, "of", total_models, "=> Performance:", performance_100, "%" )

    return performance_100




def check_top10(top10, correct_models_dict):
    correct_models = 0
    total_models = 0

    for entry in top10:
        model = entry.split()[0]
        if model in correct_models_dict.keys():
            correct_models = correct_models + 1
        total_models = total_models + 1

    performance_10 = (float(correct_models)/total_models) * 100
    print("Number of correct models:", correct_models, "of", total_models, "=> Performance:", performance_10, "%" )

    return performance_10





def check_top1(top1, correct_models_dict):
    correct_models = 0
    total_models = 0

    model = top1.split()[0]


    if model in correct_models_dict.keys():
        correct_models = correct_models + 1

    total_models = total_models + 1

    performance_1 = (float(correct_models)/total_models) * 100
    print("Number of correct models:", correct_models, "of", total_models, "=> Performance:", performance_1, "%" )


    return performance_1




def join_dicts(dict_023_performance, dict_049_performance, dict_080_performance):
    dict_023_049_080 = defaultdict(list)

    for target, tuple in dict_023_performance.items():
        dict_023_049_080[target].append(tuple)
        dict_023_049_080[target].append(dict_049_performance[target])
        dict_023_049_080[target].append(dict_080_performance[target])

    return dict_023_049_080








def create_csv(dict_023_049_080):
    print("Creating csv file...")
    csv_dict = {}

    f = open('/proj/wallner/users/x_karst/exjobb/tables/myperformance_propro_on_propep.csv','w')

    f.write("Target,Top1 0.23 [%],Top10 0.23 [%],Top100 0.23 [%],Top1 0.49 [%],Top10 0.49 [%],Top100 0.49 [%],Top1 0.80 [%],Top10 0.80 [%],Top100 0.80 [%]\n")

    for target, tuples in dict_023_049_080.items():

        target= target.split("_")[2] # myscore_propep_4uqyBA11.txt -> 4uqyBA11.txt
        target = target.split(".")[0] #4uqyBA11.txt -> 4uqyBA11



        first = tuples[0]
        second = tuples[1]
        third = tuples[2]

        top1_first = int(first[0])
        top10_first = int(first[1])
        top100_first = int(first[2])

        top1_second = int(second[0])
        top10_second = int(second[1])
        top100_second = int(second[2])

        top1_third = int(third[0])
        top10_third = int(third[1])
        top100_third = int(third[2])

        f.write(str(target)+","+str(top1_first)+","+str(top10_first)+","+str(top100_first)+","+str(top1_second)+","+str(top10_second)+","+str(top100_second)+","+str(top1_third)+","+str(top10_third)+","+str(top100_third)+'\n')



    f.close()
    print("Done")


    return


def main():

    args = sys.argv[1:]



    targets_list = find_targets()   #list of files, one file contains scores for all models of that target


    print("----------------------BENCHMARK 0.23 ----------------------------")
    dict_023_performance = sort_scores_023(targets_list)
    print("----------------------BENCHMARK 0.49 ----------------------------")
    dict_049_performance = sort_scores_049(targets_list)
    print("----------------------BENCHMARK 0.80 ----------------------------")
    dict_080_performance = sort_scores_080(targets_list)



    dict_023_049_080 = join_dicts(dict_023_performance, dict_049_performance, dict_080_performance)


    create_csv(dict_023_049_080)




















if __name__ == '__main__':
    main()
