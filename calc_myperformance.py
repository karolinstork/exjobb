import os
import sys
from os import walk
import pickle
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import csv







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



    if which_set == "propep_sym":
        targets_list = []
        for (dir_path, dirnames, filenames) in walk(f"/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/pro_pep/myscores_sym/"):
            for file in filenames:
                target_path = dir_path+file
                #print(target_path)
                targets_list.append(target_path)


    if which_set == "propep_asym":
        targets_list = []
        for (dir_path, dirnames, filenames) in walk(f"/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/pro_pep/myscores_asym/"):
            for file in filenames:
                target_path = dir_path+file
                targets_list.append(target_path)


    return targets_list #find all targets. for propro: CAPRI and MOAL. For propep it is in either asymmetric or symmetric scores.



def sort_scores_023(targets_list, which_set):
    dict_023_performance = {}

    if which_set == "propro":
        correct_models_dict = pickle.load(open("/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/pro_pro/dockq_natives_0.23.p", "rb"))
    if which_set == "propep_asym" or which_set == "propep_sym":
        correct_models_dict = pickle.load(open("/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/pro_pep/dockq_natives_0.23.p", "rb"))

    for target in targets_list:
        print("------------------------------------------")
        print(target, "Benchmark: 0.23")

        if which_set == "propro":
            target_name = os.path.basename(target)
            target_file = target+"/myscore_"+target_name+".txt"

        if which_set == "propep_asym" or which_set == "propep_sym":
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


        performance_1 = check_top1(top1, correct_models_dict, which_set)
        performance_10 = check_top10(top10, correct_models_dict, which_set)
        performance_100 = check_top100(top100, correct_models_dict, which_set)


        dict_023_performance[target_name] = (performance_1, performance_10, performance_100)



    return dict_023_performance



def sort_scores_049(targets_list, which_set):
    dict_049_performance = {}

    if which_set == "propro":
        correct_models_dict = pickle.load(open("/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/pro_pro/dockq_natives_0.49.p", "rb"))
    if which_set == "propep_sym" or which_set == "propep_asym":
        correct_models_dict = pickle.load(open("/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/pro_pep/dockq_natives_0.49.p", "rb"))

    for target in targets_list:
        print("------------------------------------------")
        print(target, "Benchmark: 0.49")


        if which_set == "propro":
            target_name = os.path.basename(target)
            target_file = target+"/myscore_"+target_name+".txt"

        if which_set == "propep_asym" or which_set == "propep_sym":
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


        performance_1 = check_top1(top1, correct_models_dict, which_set)
        performance_10 = check_top10(top10, correct_models_dict, which_set)
        performance_100 = check_top100(top100, correct_models_dict, which_set)


        dict_049_performance[target_name] = (performance_1, performance_10, performance_100)

    return dict_049_performance




def sort_scores_080(targets_list, which_set):
    dict_080_performance = {}

    if which_set == "propro":
        correct_models_dict = pickle.load(open("/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/pro_pro/dockq_natives_0.8.p", "rb"))
    if which_set == "propep_sym" or which_set == "propep_asym":
        correct_models_dict = pickle.load(open("/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/pro_pep/dockq_natives_0.8.p", "rb"))

    for target in targets_list:
        print("------------------------------------------")
        print(target, "Benchmark: 0.80")


        if which_set == "propro":
            target_name = os.path.basename(target)
            target_file = target+"/myscore_"+target_name+".txt"

        if which_set == "propep_asym" or which_set == "propep_sym":
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

        if which_set == "propep_sym" or which_set == "propep_asym":
            if model in correct_models_dict.keys():
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


        if which_set == "propep_sym" or which_set == "propep_asym":
            if model in correct_models_dict.keys():
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

    if which_set == "propep_sym" or which_set == "propep_asym":
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


def plot_bars_propro(dict_023_049_080):
    top1_023_MOAL = []
    top1_049_MOAL = []
    top1_080_MOAL = []

    top10_023_MOAL = []
    top10_049_MOAL = []
    top10_080_MOAL = []

    top100_023_MOAL = []
    top100_049_MOAL = []
    top100_080_MOAL = []

    targets_list_MOAL = []

    top1_023_CAPRI = []
    top1_049_CAPRI = []
    top1_080_CAPRI = []

    top10_023_CAPRI = []
    top10_049_CAPRI = []
    top10_080_CAPRI = []

    top100_023_CAPRI = []
    top100_049_CAPRI = []
    top100_080_CAPRI = []

    targets_list_CAPRI = []

    check = 0
    edition = 0

    for target, tuples in dict_023_049_080.items():
        tuple023 = tuples[0]
        tuple049 = tuples[1]
        tuple080 = tuples[2]

        top1_023 = tuple023[0]
        top10_023 = tuple023[1]
        top100_023 = tuple023[2]

        top1_049 = tuple049[0]
        top10_049 = tuple049[1]
        top100_049 = tuple049[2]

        top1_080 = tuple080[0]
        top10_080 = tuple080[1]
        top100_080 = tuple080[2]

        if target[0] == "D":
            top1_023_MOAL.append(top1_023)
            top10_023_MOAL.append(top10_023)
            top100_023_MOAL.append(top100_023)

            top1_049_MOAL.append(top1_049)
            top10_049_MOAL.append(top10_049)
            top100_049_MOAL.append(top100_049)

            top1_080_MOAL.append(top1_080)
            top10_080_MOAL.append(top10_080)
            top100_080_MOAL.append(top100_080)

            targets_list_MOAL.append(target)

            check = check + 1

            if (check == 20) or (check == 40) or (check == 60) or (check == 80) or (check == 100) or check == 118 :
                edition = edition + 1

                barWidth = 0.25

                r = np.arange(len(top1_023_MOAL))
                plt.figure(figsize = (15.4, 7.0), dpi = 300)

                bars_023 = np.add(top1_023_MOAL, top10_023_MOAL).tolist()
                bars_049 = np.add(top1_049_MOAL, top10_049_MOAL).tolist()
                bars_080 = np.add(top1_080_MOAL, top10_080_MOAL).tolist()

                r1 = np.arange(len(top1_023_MOAL))
                r2 = [x + barWidth for x in r1]
                r3 = [x + barWidth for x in r2]

                plt.bar(r1, top1_023_MOAL, color='#FFD435', edgecolor='white', width=barWidth, label = "Top1")
                plt.bar(r1, top10_023_MOAL, bottom=top1_023_MOAL, color='#FF9735', edgecolor='white', width=barWidth, label = "Top10")
                plt.bar(r1, top100_023_MOAL, bottom=bars_023, color='#FC84B5', edgecolor='white', width=barWidth, label = "Top100")

                plt.bar(r2, top1_049_MOAL, color='#FFD435', edgecolor='white', hatch = "////", width=barWidth)
                plt.bar(r2, top10_049_MOAL, bottom=top1_049_MOAL, color='#FF9735', edgecolor='white', hatch = "////", width=barWidth)
                plt.bar(r2, top100_049_MOAL, bottom=bars_049, color='#FC84B5', edgecolor='white', hatch = "////", width=barWidth)


                plt.bar(r3, top1_080_MOAL, color='#FFD435', edgecolor='white', hatch = "//", width=barWidth)
                plt.bar(r3, top10_080_MOAL, bottom=top1_080_MOAL, color='#FF9735', edgecolor='white', hatch = "//", width=barWidth)
                plt.bar(r3, top100_080_MOAL, bottom=bars_080, color='#FC84B5', edgecolor='white', hatch = "//", width=barWidth)



                plt.ylim((0,305))
                plt.xticks([r + barWidth for r in range(len(top1_023_MOAL))], targets_list_MOAL)
                plt.title(f'Performance on MOAL set ({edition}/6) with DockQ benchmark 0.23, 0.49 and 0.80', fontweight='bold')
                plt.xlabel("Target complexes")
                plt.ylabel('Accuracy [%]')
                plt.savefig(f"/proj/wallner/users/x_karst/exjobb/pictures/my_performance_MOAL_{edition}.png")

                plt.clf()


                top1_023_MOAL = []
                top1_049_MOAL = []
                top1_080_MOAL = []

                top10_023_MOAL = []
                top10_049_MOAL = []
                top10_080_MOAL = []

                top100_023_MOAL = []
                top100_049_MOAL = []
                top100_080_MOAL = []

                targets_list_MOAL = []


        if target[0] == "T":
            top1_023_CAPRI.append(top1_023)
            top10_023_CAPRI.append(top10_023)
            top100_023_CAPRI.append(top100_023)

            top1_049_CAPRI.append(top1_049)
            top10_049_CAPRI.append(top10_049)
            top100_049_CAPRI.append(top100_049)

            top1_080_CAPRI.append(top1_080)
            top10_080_CAPRI.append(top10_080)
            top100_080_CAPRI.append(top100_080)

            targets_list_CAPRI.append(target)


    barWidth = 0.25
    r = np.arange(len(top1_023_CAPRI))
    plt.figure(figsize = (9.4, 6.0))
    bars_023 = np.add(top1_023_CAPRI, top10_023_CAPRI).tolist()
    bars_049 = np.add(top1_049_CAPRI, top10_049_CAPRI).tolist()
    bars_080 = np.add(top1_080_CAPRI, top10_080_CAPRI).tolist()

    r1 = np.arange(len(top1_023_CAPRI))
    r2 = [x + barWidth for x in r1]
    r3 = [x + barWidth for x in r2]

    plt.bar(r1, top1_023_CAPRI, color='#FFD435', edgecolor='white', width=barWidth, label = "Top1")
    plt.bar(r1, top10_023_CAPRI, bottom=top1_023_CAPRI, color='#FF9735', edgecolor='white', width=barWidth, label = "Top10")
    plt.bar(r1, top100_023_CAPRI, bottom=bars_023, color='#FC84B5', edgecolor='white', width=barWidth, label = "Top100")

    plt.bar(r2, top1_049_CAPRI, color='#FFD435', edgecolor='white', hatch = "////", width=barWidth)
    plt.bar(r2, top10_049_CAPRI, bottom=top1_049_CAPRI, color='#FF9735', edgecolor='white', hatch = "////", width=barWidth)
    plt.bar(r2, top100_049_CAPRI, bottom=bars_049, color='#FC84B5', edgecolor='white', hatch = "////", width=barWidth)


    plt.bar(r3, top1_080_CAPRI, color='#FFD435', edgecolor='white', hatch = "//", width=barWidth)
    plt.bar(r3, top10_080_CAPRI, bottom=top1_080_CAPRI, color='#FF9735', edgecolor='white', hatch = "//", width=barWidth)
    plt.bar(r3, top100_080_CAPRI, bottom=bars_080, color='#FC84B5', edgecolor='white', hatch = "//", width=barWidth)



    plt.ylim((0,305))
    plt.xticks([r + barWidth for r in range(len(top1_023_CAPRI))], targets_list_CAPRI)
    plt.xlabel("Target complexes")
    plt.title(f'Performance on CAPRI set with DockQ benchmark 0.23, 0.49 and 0.80', fontweight='bold')
    plt.ylabel('Accuracy [%]')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"/proj/wallner/users/x_karst/exjobb/pictures/my_performance_CAPRI.png")

    plt.clf()


    return



def plot_bars_propep(dict_023_049_080, which_set):
    top1_023_list = []
    top1_049_list = []
    top1_080_list = []

    top10_023_list = []
    top10_049_list = []
    top10_080_list = []

    top100_023_list = []
    top100_049_list = []
    top100_080_list = []

    targets_list = []
    check = 0
    edition = 0

    for target, tuples in dict_023_049_080.items():
        tuple023 = tuples[0]
        tuple049 = tuples[1]
        tuple080 = tuples[2]

        top1_023 = tuple023[0]
        top10_023 = tuple023[1]
        top100_023 = tuple023[2]

        top1_049 = tuple049[0]
        top10_049 = tuple049[1]
        top100_049 = tuple049[2]

        top1_080 = tuple080[0]
        top10_080 = tuple080[1]
        top100_080 = tuple080[2]

        top1_023_list.append(top1_023)
        top10_023_list.append(top10_023)
        top100_023_list.append(top100_023)

        top1_049_list.append(top1_049)
        top10_049_list.append(top10_049)
        top100_049_list.append(top100_049)

        top1_080_list.append(top1_080)
        top10_080_list.append(top10_080)
        top100_080_list.append(top100_080)

        target= target.split("_")[2] # myscore_propep_4uqyBA11.txt
        target = target.split(".")[0] #4uqyBA11.txt
        targets_list.append(target) #4uqyBA11

        check = check + 1

        if check == 18 or check == 36:
            edition = edition + 1


            barWidth = 0.25
            r = np.arange(len(top1_023_list))
            plt.figure(figsize = (12.4, 6.0))
            bars_023 = np.add(top1_023_list, top10_023_list).tolist()
            bars_049 = np.add(top1_049_list, top10_049_list).tolist()
            bars_080 = np.add(top1_080_list, top10_080_list).tolist()

            r1 = np.arange(len(top1_023_list))
            r2 = [x + barWidth for x in r1]
            r3 = [x + barWidth for x in r2]

            plt.bar(r1, top1_023_list, color='#FFD435', edgecolor='white', width=barWidth, label = "Top1")
            plt.bar(r1, top10_023_list, bottom=top1_023_list, color='#FF9735', edgecolor='white', width=barWidth, label = "Top10")
            plt.bar(r1, top100_023_list, bottom=bars_023, color='#FC84B5', edgecolor='white', width=barWidth, label = "Top100")

            plt.bar(r2, top1_049_list, color='#FFD435', edgecolor='white', hatch = "////", width=barWidth)
            plt.bar(r2, top10_049_list, bottom=top1_049_list, color='#FF9735', edgecolor='white', hatch = "////", width=barWidth)
            plt.bar(r2, top100_049_list, bottom=bars_049, color='#FC84B5', edgecolor='white', hatch = "////", width=barWidth)


            plt.bar(r3, top1_080_list, color='#FFD435', edgecolor='white', hatch = "//", width=barWidth)
            plt.bar(r3, top10_080_list, bottom=top1_080_list, color='#FF9735', edgecolor='white', hatch = "//", width=barWidth)
            plt.bar(r3, top100_080_list, bottom=bars_080, color='#FC84B5', edgecolor='white', hatch = "//", width=barWidth)



            plt.ylim((0,305))
            plt.xticks([r + barWidth for r in range(len(top1_023_list))], targets_list, fontsize = 7.5)
            plt.xlabel("Target complexes")

            if which_set == "propep_sym":
                plt.title(f'Performance of symmetric propep matrix with DockQ benchmark 0.23, 0.49 and 0.80 ({edition}/2)', fontweight='bold')
            if which_set == "propep_asym":
                plt.title(f'Performance of asymmetric propep matrix with DockQ benchmark 0.23, 0.49 and 0.80 ({edition}/2)', fontweight='bold')

            plt.ylabel('Accuracy [%]')
            plt.legend()
            plt.tight_layout()
            plt.savefig(f"/proj/wallner/users/x_karst/exjobb/pictures/my_performance_{which_set}_{edition}.png")
            plt.show()
            plt.clf()


            top1_023_list = []
            top1_049_list = []
            top1_080_list = []

            top10_023_list = []
            top10_049_list = []
            top10_080_list = []

            top100_023_list = []
            top100_049_list = []
            top100_080_list = []

            targets_list = []

    return



def create_csv(dict_023_049_080, which_set):
    csv_dict = {}
    if which_set == "propro":
        f = open('/proj/wallner/users/x_karst/exjobb/tables/myperformance_CnM_023_049_080.csv','w')

    if which_set == "propep_sym":
        f = open('/proj/wallner/users/x_karst/exjobb/tables/myperformance_sympeptide_023_049_080.csv','w')

    if which_set == "propep_asym":
        f = open('/proj/wallner/users/x_karst/exjobb/tables/myperformance_asympeptide_023_049_080.csv','w')


    f.write("Target,Top1 0.23 [%],Top10 0.23 [%],Top100 0.23 [%],Top1 0.49 [%],Top10 0.49 [%],Top100 0.80 [%],Top1 0.80 [%],Top10 0.80 [%],Top100 0.80 [%]\n")

    for target, tuples in dict_023_049_080.items():
        if which_set == "propep_sym" or which_set == "propep_asym":
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



    return


def main():
    try:
        args = sys.argv[1:]
        which_set = args[0] #propro or propep_asym or propep_sym


    except IndexError:
        print("USAGE: [propro] [propep_sym] [propep_asym]")
        return

    targets_list = find_targets(which_set)   #list of files, one file contains scores for all models of that target


    print("----------------------BENCHMARK 0.23 ----------------------------")
    dict_023_performance = sort_scores_023(targets_list, which_set)
    print("----------------------BENCHMARK 0.49 ----------------------------")
    dict_049_performance = sort_scores_049(targets_list, which_set)
    print("----------------------BENCHMARK 0.80 ----------------------------")
    dict_080_performance = sort_scores_080(targets_list, which_set)



    dict_023_049_080 = join_dicts(dict_023_performance, dict_049_performance, dict_080_performance)

    if which_set == "propro":
        plot_bars_propro(dict_023_049_080)
    if which_set == "propep_sym" or which_set == "propep_asym":
        plot_bars_propep(dict_023_049_080, which_set)

    create_csv(dict_023_049_080, which_set)




















if __name__ == '__main__':
    main()
