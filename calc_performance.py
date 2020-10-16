import os
import sys
from os import walk
import pickle
import numpy as np
import matplotlib.pyplot as plt









def find_targets(which_set):
    if which_set == "propro":
        targets_list = []
        for (root, dirs, files) in walk(f"/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/pro_pro/CAPRI/"):
            for dir in dirs:
                target_path = (root+dir)
                targets_list.append(target_path)


        for (dir_path, dirnames, filenames) in walk(f"/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/pro_pro/MOAL/"):
            for dir in dirnames:
                target_path = (dir_path+dir)
                targets_list.append(target_path)



    if which_set == "propep":
        targets_list = []
        for (dir_path, dirnames, filenames) in walk(f"/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/pro_pep/my_scores_output_symmetric/"):
            for file in filenames:
                target_path = dir_path+file
                targets_list.append(target_path)



    return targets_list #find all targets. for propro: CAPRI and MOAL. For propep it is in either asymmetric or symmetric scores.



def sort_scores_023(targets_list, which_set):
    dict_023_performance = {}

    if which_set == "propro":
        correct_models_dict = pickle.load(open("/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/pro_pro/dockq_natives_0.23.p", "rb"))
    if which_set == "propep":
        correct_models_dict = pickle.load(open("/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/pro_pep/dockq_natives_0.23.p", "rb"))

    for target in targets_list:
        print("------------------------------------------")
        print(target)
        target_name = os.path.basename(target)
        target_file = target+"/myscore_"+target_name+".txt"

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


        performance_1 = check_top1(top1, correct_models_dict, which_set)
        performance_10 = check_top10(top10, correct_models_dict, which_set)
        performance_100 = check_top100(top100, correct_models_dict, which_set)


        dict_023_performance[target_name] = (performance_1, performance_10, performance_100)



    return dict_023_performance



def sort_scores_049(targets_list, which_set):
    dict_049_performance = {}

    if which_set == "propro":
        correct_models_dict = pickle.load(open("/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/pro_pro/dockq_natives_0.49.p", "rb"))
    if which_set == "propep":
        correct_models_dict = pickle.load(open("/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/pro_pep/dockq_natives_0.49.p", "rb"))

    for target in targets_list:
        print("------------------------------------------")
        print(target)
        target_name = os.path.basename(target)
        target_file = target+"/myscore_"+target_name+".txt"

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


        performance_1 = check_top1(top1, correct_models_dict, which_set)
        performance_10 = check_top10(top10, correct_models_dict, which_set)
        performance_100 = check_top100(top100, correct_models_dict, which_set)


        dict_049_performance[target_name] = (performance_1, performance_10, performance_100)

    return dict_049_performance




def sort_scores_080(targets_list, which_set):
    dict_080_performance = {}

    if which_set == "propro":
        correct_models_dict = pickle.load(open("/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/pro_pro/dockq_natives_0.8.p", "rb"))
    if which_set == "propep":
        correct_models_dict = pickle.load(open("/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/pro_pep/dockq_natives_0.8.p", "rb"))

    for target in targets_list:
        print("------------------------------------------")
        print(target)
        target_name = os.path.basename(target)
        target_file = target+"/myscore_"+target_name+".txt"

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


        performance_1 = check_top1(top1, correct_models_dict, which_set)
        performance_10 = check_top10(top10, correct_models_dict, which_set)
        performance_100 = check_top100(top100, correct_models_dict, which_set)


        dict_080_performance[target_name] = (performance_1, performance_10, performance_100)

    return dict_080_performance



def check_top100(top100, correct_models_dict, which_set):
    correct_models = 0
    total_models = 0

    for entry in top100:
        model = entry.split()[0]
        if which_set == "propro":
            if model[0] == "T": #from CAPRI data set
                first_name = model.split("-")[0]
                first_name = first_name[6:]
                first_name = "T"+first_name
                last_name = model.split("-")[1]
                last_name = last_name.split(".")[0]
                last_name = last_name[5:]
                name = first_name+"-"+last_name

            if model[0] == "D": #from MOAL data set
                name_parts = model.split("-")
                name = name_parts[0] +"-"+ name_parts[1]


            if name in correct_models_dict.keys():
                correct_models = correct_models + 1

            total_models = total_models + 1

    performance_100 = (float(correct_models)/total_models) * 100
    print("Number of correct models:", correct_models, "of", total_models, "=> Performance:", performance_100, "%" )

    return performance_100




def check_top10(top10, correct_models_dict, which_set):
    correct_models = 0
    total_models = 0
    for entry in top10:
        model = entry.split()[0]
        if which_set == "propro":
            if model[0] == "T": #from CAPRI data set
                first_name = model.split("-")[0]
                first_name = first_name[6:]
                first_name = "T"+first_name
                last_name = model.split("-")[1]
                last_name = last_name.split(".")[0]
                last_name = last_name[5:]
                name = first_name+"-"+last_name

            if model[0] == "D": #from MOAL data set
                name_parts = model.split("-")
                name = name_parts[0] + "-"+ name_parts[1]

            if name in correct_models_dict.keys():
                correct_models = correct_models + 1


            total_models = total_models + 1

    performance_10 = (float(correct_models)/total_models) * 100
    print("Number of correct models:", correct_models, "of", total_models, "=> Performance:", performance_10, "%" )

    return performance_10





def check_top1(top1, correct_models_dict, which_set):
    correct_models = 0
    total_models = 0

    model = top1.split()[0]

    if which_set == "propro":
        if model[0] == "T": #from CAPRI data set
            first_name = model.split("-")[0]
            first_name = first_name[6:]
            first_name = "T"+first_name
            last_name = model.split("-")[1]
            last_name = last_name.split(".")[0]
            last_name = last_name[5:]
            name = first_name+"-"+last_name

        elif model[0] == "D": #from MOAL data set ie D2QFW-a117b-merged.pdb. In facit D2QFW-a117b-merged.pdb
            name_parts = model.split("-")
            name = name_parts[0] + "-" + name_parts[1]


        if name in correct_models_dict.keys():
            correct_models = correct_models + 1


        total_models = total_models + 1

    performance_1 = (float(correct_models)/total_models) * 100
    print("Number of correct models:", correct_models, "of", total_models, "=> Performance:", performance_1, "%" )


    return performance_1


def plot_bars(performance_dict):
    top1_bars_MOAL = []
    top10_bars_MOAL = []
    top100_bars_MOAL = []
    targets_list_MOAL = []

    top1_bars_CAPRI = []
    top10_bars_CAPRI = []
    top100_bars_CAPRI = []
    targets_list_CAPRI = []


    for key, value in performance_dict.items():
        top1_value = value[0]
        top10_value = value[1]
        top100_value = value[2]

        if key[0] == "D":
            key = key[0]
            top1_bars_MOAL.append(top1_value)
            top10_bars_MOAL.append(top10_value)
            top100_bars_MOAL.append(top100_value)
            targets_list_MOAL.append(key)


        if key[0] == "T":
            top1_bars_CAPRI.append(top1_value)
            top10_bars_CAPRI.append(top10_value)
            top100_bars_CAPRI.append(top100_value)
            targets_list_CAPRI.append(key)



    barWidth = 0.25
    r1 = np.arange(len(top1_bars_MOAL))
    r2 = [x + barWidth for x in r1]
    r3 = [x + barWidth for x in r2]
    plt.bar(r1, top1_bars_MOAL, color='#FCD484', width=barWidth, edgecolor='white', label='Top1')
    plt.bar(r2, top10_bars_MOAL, color='#FC9F84', width=barWidth, edgecolor='white', label='Top10')
    plt.bar(r3, top100_bars_MOAL, color='#FC84B5', width=barWidth, edgecolor='white', label='Top100')

    plt.xlabel('MOAL set with DockQ benchmark XX', fontweight='bold')
    plt.ylabel('Performance [%]', fontweight = 'bold')

    plt.xticks([r + barWidth for r in range(len(top1_bars_MOAL))], targets_list_MOAL)

    plt.legend()
    plt.show()




    barWidth = 0.25
    r1 = np.arange(len(top1_bars_CAPRI))
    r2 = [x + barWidth for x in r1]
    r3 = [x + barWidth for x in r2]
    plt.bar(r1, top1_bars_CAPRI, color='#FCD484', width=barWidth, edgecolor='white', label='Top1')
    plt.bar(r2, top10_bars_CAPRI, color='#FC9F84', width=barWidth, edgecolor='white', label='Top10')
    plt.bar(r3, top100_bars_CAPRI, color='#FC84B5', width=barWidth, edgecolor='white', label='Top100')

    plt.xlabel('CAPRI set with DockQ benchmark XX', fontweight='bold')
    plt.ylabel('Performance [%]', fontweight = 'bold')
    plt.xticks([r + barWidth for r in range(len(top1_bars_CAPRI))], targets_list_CAPRI)

    plt.legend()
    plt.show()

    return



def main():
    try:
        args = sys.argv[1:]
        which_set = args[0] #propro or propep


    except IndexError:
        print("USAGE: <file.csv> [propro] [propep]")
        return

    targets_list = find_targets(which_set)   #list of files, one file contains scores for all models of that target


    print("----------------------BENCHMARK 0.23 ----------------------------")
    dict_023_performance = sort_scores_023(targets_list, which_set)
    #plot_bars(dict_023_performance)
    print('\n')
    print('\n')
    print('\n')
    print('\n')





    print("----------------------BENCHMARK 0.49 ----------------------------")
    dict_049_performance = sort_scores_049(targets_list, which_set)
    plot_bars(dict_049_performance)
    print('\n')
    print('\n')
    print('\n')
    print('\n')

    print("----------------------BENCHMARK 0.80 ----------------------------")


    dict_080_performance = sort_scores_080(targets_list, which_set)



#
#




















if __name__ == '__main__':
    main()
