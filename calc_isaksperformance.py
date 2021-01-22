import os
import sys
from os import walk
import pickle
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import csv




def find_targets():
    unavailable_targets = ["1gy3", "4imi", "4xvj", "1kcr", "4h3h", "3gz1", "4l1u", "4ajy", "5fgc", "1z9o"] #targets that were not possible to use for the scoring matrices
    targets_list = []

    for (dir_path, dirnames, filenames) in walk(f"/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/propep_matrices/isaks_zrank"):
        for file in filenames:
            if file[:4] not in unavailable_targets:
                target_path = dir_path+"/"+file
                targets_list.append(target_path)

    return targets_list


def sort_file(target):
    file = open(target, "r")
    lines = file.readlines()

    sorted_zrank = sorted(lines,key=lambda line: line.split()[1]) #low value is better quality model than high value in zrank, default for sorted function is ascending

    return sorted_zrank



def calc_performance(sorted_zrank, target):
    top100_zrank = sorted_zrank[:100] #Generates top lists from the full list that is sorted according to zrank values
    top10_zrank = sorted_zrank[:10]
    top1_zrank = sorted_zrank[0]

    correct_100_zrank_023 = 0
    correct_100_zrank_049 = 0
    correct_100_zrank_080 = 0

    correct_10_zrank_023 = 0
    correct_10_zrank_049 = 0
    correct_10_zrank_080 = 0

    correct_1_zrank_023 = 0
    correct_1_zrank_049 = 0
    correct_1_zrank_080 = 0

    for model in top100_zrank:
        correct_100_zrank_023, correct_100_zrank_049, correct_100_zrank_080 = check_pickledict(model, target, correct_100_zrank_023, correct_100_zrank_049, correct_100_zrank_080) #walks through top100 and checks how many models are in pickle dicts

    for model in top10_zrank:
        correct_10_zrank_023, correct_10_zrank_049, correct_10_zrank_080 = check_pickledict(model, target, correct_10_zrank_023, correct_10_zrank_049, correct_10_zrank_080) #walks through top10 and checks how many models are in pickle dicts


    correct_1_zrank_023, correct_1_zrank_049, correct_1_zrank_080 = check_pickledict(top1_zrank, target, correct_1_zrank_023, correct_1_zrank_049, correct_1_zrank_080)  #takes top1 and checks if it is in pickle dicts






    ####################################################################################


    zrank_023 = (correct_1_zrank_023, correct_10_zrank_023, correct_100_zrank_023)
    zrank_049 = (correct_1_zrank_049, correct_10_zrank_049, correct_100_zrank_049)
    zrank_080 = (correct_1_zrank_080, correct_10_zrank_080, correct_100_zrank_080)


    return zrank_023, zrank_049, zrank_080





def check_pickledict(model, target, correct_023, correct_049, correct_080): #checks if model is a correct model according to all three benchmarks
    correct_models_023 = pickle.load(open("/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/propep_matrices/dockq_natives_0.23.p", "rb"))
    correct_models_049 = pickle.load(open("/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/propep_matrices/dockq_natives_0.49.p", "rb"))
    correct_models_080 = pickle.load(open("/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/propep_matrices/dockq_natives_0.8.p", "rb"))

    target = os.path.basename(target)
    target = target.split(".")[0] #1elwCA01
    target= target[:-4] #1elw

    split1 = model.split(".")[1]
    model_numb = split1.split("_")[0]

    model_name = target+"_"+model_numb

    if model_name in correct_models_023.keys():
        #print("MATCH in 0.23", model_name)
        correct_023 = correct_023 + 1

    if model_name in correct_models_049.keys():
        #print("MATCH in 0.49", model_name)
        correct_049 = correct_049 + 1

    if model_name in correct_models_080.keys():
        #print("MATCH in 0.80", model_name)
        correct_080 = correct_080 + 1


    return correct_023, correct_049, correct_080






def plot_performance(zrank_targets023, zrank_targets049, zrank_targets080, tot_targets):
    top1_023 = [zrank_targets023[0], 7, 4, 8, 7]
    top10_023 = [zrank_targets023[1], 10, 13, 13, 16]
    top100_023 = [zrank_targets023[2], 20, 19, 23, 22]

    top1_049 = [zrank_targets049[0], 1, 0, 1, 2]
    top10_049 = [zrank_targets049[1], 2, 3, 4, 4]
    top100_049 = [zrank_targets049[2], 8, 7, 10, 9]

    top1_080 = [zrank_targets080[0], 0, 0, 0, 0]
    top10_080 = [zrank_targets080[1], 0, 0, 0, 0]
    top100_080 = [zrank_targets080[2], 0, 0, 1, 1]

    bars10 = top1_023 + top10_023 + top100_023 + top1_049 + top10_049 + top100_049 + top1_080 + top10_080 + top100_080


    x_list = ["ZRANK", "Glaser matrix", "True propro matrix", "Symmetric propep matrix", "Asymmetric propep matrix"]


    barWidth = 0.08

    r1 = list(np.arange(len(top1_023)))

    r2 = [x + barWidth for x in r1]
    r3 = [x + barWidth for x in r2]
    r4 = [x + barWidth for x in r3]
    r5 = [x + barWidth for x in r4]
    r6 = [x + barWidth for x in r5]
    r7 = [x + barWidth for x in r6]
    r8 = [x + barWidth for x in r7]
    r9 = [x + barWidth for x in r8]

    r10 = r1+r2+r3+r4+r5+r6+r7+r8+r9

    plt.figure(figsize = (9.4, 6.0), dpi = 190)
    plt.bar(r1, top1_023, color='#FECD1C', edgecolor='white', width=barWidth, label = "Acceptable top 1")
    plt.bar(r2, top10_023, color='#FFDC5D', edgecolor='white', width=barWidth, label = "Acceptable top 10")
    plt.bar(r3, top100_023, color='#FFE99A', edgecolor='white', width=barWidth, label = "Acceptable top 100")

    plt.bar(r4, top1_049, color='#FF7D00', edgecolor='white', width=barWidth, label = "Medium top 1")
    plt.bar(r5, top10_049, color='#FF9735', edgecolor='white', width=barWidth, label = "Medium top 10")
    plt.bar(r6, top100_049, color='#FCB775', edgecolor='white', width=barWidth, label = "Medium top 100")

    plt.bar(r7, top1_080, color='#FF62A2', edgecolor='white', width=barWidth, label = "High top 1")
    plt.bar(r8, top10_080, color='#FC84B5', edgecolor='white', width=barWidth, label = "High top 10")
    plt.bar(r9, top100_080, color='#FFC7DE', edgecolor='white', width=barWidth, label = "High top 100")

    plt.ylim((0,40))


    labels = list(top1_023 + top10_023 + top100_023 + top1_049 + top10_049 + top100_049 + top1_080 + top10_080 + top100_080)


    for i in range(len(r10)):
        plt.text(x = r10[i]-0.03, y = bars10[i]+0.5, s = labels[i], size = 6)



    plt.xticks(r5, x_list, size = 8)


    plt.xlabel("Methods", fontweight = "bold")
    plt.title(f'Performance on {tot_targets} targets with 3 DockQ quality benchmarks', fontweight='bold')
    plt.ylabel('Number of targets with >=1 model hits', fontweight = "bold")
    plt.tight_layout()
    plt.legend(fontsize = "small")
    plt.savefig(f"/proj/wallner/users/x_karst/exjobb/pictures/isaks_performance.png")
    plt.show()
    plt.clf()



    return






def main():

    args = sys.argv[1:]
    targets_list = find_targets()




    tot_targets = 0
    zrank_targets023 = [0, 0, 0]
    zrank_targets049 = [0, 0, 0]
    zrank_targets080 = [0, 0, 0]


    for target in targets_list:
        print(target)
        sorted_zrank = sort_file(target)
        zrank_023, zrank_049, zrank_080 = calc_performance(sorted_zrank, target)

        if int(zrank_023[0])>0:
            zrank_targets023[0] = zrank_targets023[0] + 1
        if int(zrank_023[1])>0:
            zrank_targets023[1] = zrank_targets023[1] + 1
        if int(zrank_023[2])>0:
            zrank_targets023[2] = zrank_targets023[2] + 1
        if int(zrank_049[0])>0:
            zrank_targets049[0] = zrank_targets049[0] + 1
        if int(zrank_049[1])>0:
            zrank_targets049[1] = zrank_targets049[1] + 1
        if int(zrank_049[2])>0:
            zrank_targets049[2] = zrank_targets049[2] + 1
        if int(zrank_080[0])>0:
            zrank_targets080[0] = zrank_targets080[0] + 1
        if int(zrank_080[1])>0:
            zrank_targets080[1] = zrank_targets080[1] + 1
        if int(zrank_080[2])>0:
            zrank_targets080[2] = zrank_targets080[2] + 1

        tot_targets = tot_targets + 1





    print("-------ZRANK TOTAL PERFORMANCE OUT OF", tot_targets, "-------")
    print("0.23 <top1, top10, top100>:", zrank_targets023)
    print("0.49 <top1, top10, top100>:", zrank_targets049)
    print("0.80 <top1, top10, top100>:", zrank_targets080)



    plot_performance(zrank_targets023, zrank_targets049, zrank_targets080, tot_targets)




if __name__ == '__main__':
    main()
