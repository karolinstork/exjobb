import os
import sys
from os import walk
import pickle
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import csv




def find_targets():
    unavailable_targets = ["1gy3", "4imi", "4xvj", "1kcr", "4h3h", "3gz1", "4l1u", "4ajy", "5fgc", "1z9o"]
    targets_list = []

    for (dir_path, dirnames, filenames) in walk(f"/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/pro_pep/isaks"):
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




    print("------------ZRANK PERFORMANCE ON TARGET-------------")
    print("\tTop1\tTop 10\tTop 100")
    print("0.23:"+'\t'+ str((correct_1_zrank_023*100))+"%"+'\t'+ str((correct_10_zrank_023/10)*100)+"%"+'\t'+ str(correct_100_zrank_023)+"%")
    print("0.49:"+'\t'+ str((correct_1_zrank_049*100))+"%"+'\t'+ str((correct_10_zrank_049/10)*100)+"%"+'\t'+ str(correct_100_zrank_049)+"%")
    print("0.80:"+'\t'+ str((correct_1_zrank_080*100))+"%"+'\t'+ str((correct_10_zrank_080/10)*100)+"%"+'\t'+ str(correct_100_zrank_080)+"%")
    print('\n')



    ####################################################################################


    if correct_1_zrank_049>0:
        zrank_target_accurate_top1 = True
    else:
        zrank_target_accurate_top1 = None


    if correct_10_zrank_049>0:
        zrank_target_accurate_top10 = True
    else:
        zrank_target_accurate_top10 = None


    if correct_100_zrank_049>0:
        zrank_target_accurate_top100 = True
    else:
        zrank_target_accurate_top100 = None



    return zrank_target_accurate_top1, zrank_target_accurate_top10, zrank_target_accurate_top100






def check_pickledict(model, target, correct_023, correct_049, correct_080): #checks if model is a correct model according to all three benchmarks
    correct_models_023 = pickle.load(open("/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/pro_pep/dockq_natives_0.23.p", "rb"))
    correct_models_049 = pickle.load(open("/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/pro_pep/dockq_natives_0.49.p", "rb"))
    correct_models_080 = pickle.load(open("/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/pro_pep/dockq_natives_0.8.p", "rb"))

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






def plot_performance(succeded_predictions_zrank_top1, succeded_predictions_zrank_top10, succeded_predictions_zrank_top100, tot_targets):

    top1 = [succeded_predictions_zrank_top1, 2]
    top10 = [succeded_predictions_zrank_top10, 4]
    top100 = [succeded_predictions_zrank_top100, 9]
    bars4 = top1+top10+top100

    x_list = ["ZRANK", "Asymmetric scoring matrix"]

    r = np.arange(len(top1))
    barWidth = 0.25
    r1 = list(np.arange(len(top1)))
    r2 = [x + barWidth for x in r1]
    r3 = [x + barWidth for x in r2]
    r4 = r1+r2+r3

    print(r1)
    print(r2)
    print(r3)
    print(r4)

    plt.bar(r1, top1, color='#FFD435', edgecolor='white', width=barWidth, label = "Top1 hits")
    plt.bar(r2, top10, color='#FF9735', edgecolor='white', width=barWidth, label = "Top10 hits")
    plt.bar(r3, top100, color='#FC84B5', edgecolor='white', width=barWidth, label = "Top100 hits")


    plt.ylim((0,60))


    label = [succeded_predictions_zrank_top1, 2, succeded_predictions_zrank_top10, 4, succeded_predictions_zrank_top100, 9]
    print(label)

    for i in range(len(r4)):
        plt.text(x = r4[i]-0.03, y = bars4[i]+0.5, s = label[i], size = 8)



    plt.xticks([r + barWidth for r in range(len(x_list))], x_list)
    plt.xlabel("Methods")
    plt.title(f'Performance on {tot_targets} targets', fontweight='bold')
    plt.ylabel('Number of targets with >=1 medium quality model hits')
    plt.tight_layout()
    plt.legend()
    plt.savefig(f"/proj/wallner/users/x_karst/exjobb/pictures/isaks_performance.png")
    plt.show()
    plt.clf()



    return






def main():

    args = sys.argv[1:]
    targets_list = find_targets()


    succeded_predictions_zrank_top1 = 0
    succeded_predictions_zrank_top10 = 0
    succeded_predictions_zrank_top100 = 0

    tot_targets = 0


    for target in targets_list:
        sorted_zrank = sort_file(target)
        zrank_target_accurate_top1, zrank_target_accurate_top10, zrank_target_accurate_top100 = calc_performance(sorted_zrank, target)

        if zrank_target_accurate_top1 == True:
            succeded_predictions_zrank_top1 = succeded_predictions_zrank_top1 + 1

        if zrank_target_accurate_top10 == True:
            succeded_predictions_zrank_top10 = succeded_predictions_zrank_top10 + 1

        if zrank_target_accurate_top100 == True:
            succeded_predictions_zrank_top100 = succeded_predictions_zrank_top100 + 1

        tot_targets = tot_targets + 1




    print("------ZRANK TOTAL PERFORMANCE 0.49------")
    print("Top1:", succeded_predictions_zrank_top1, "out of", tot_targets, "successful")
    print("Top10:", succeded_predictions_zrank_top10, "out of", tot_targets, "successful")
    print("Top100:", succeded_predictions_zrank_top100, "out of", tot_targets, "successful")



    plot_performance(succeded_predictions_zrank_top1, succeded_predictions_zrank_top10, succeded_predictions_zrank_top100, tot_targets)




if __name__ == '__main__':
    main()
