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
    for (dir_path, dirnames, filenames) in walk(f"/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/propro_matrix/bjorns"):
        for file in filenames:
            if file[-4:] != ".txt":
                target_path = dir_path+"/"+file
                targets_list.append(target_path)



    return targets_list



def sort_file(target):
    file = open(target, "r")
    lines = file.readlines()

    sorted_pqd = sorted(lines,key=lambda line: line.split()[18], reverse = True)
    sorted_pqdz = sorted(lines,key=lambda line: line.split()[19], reverse = True) #reverse=true gives descending list
    sorted_zrank = sorted(lines,key=lambda line: line.split()[7]) #low value is better quality model than high value in zrank, default for sorted function is ascending


    return sorted_pqd, sorted_pqdz, sorted_zrank



def calc_performance(sorted_pqd, sorted_pqdz, sorted_zrank):

############# PQD #####################
    top100_pqd = sorted_pqd[:100] #Generates top lists from the full list that is sorted according to pqd values
    top10_pqd = sorted_pqd[:10]
    top1_pqd = sorted_pqd[0]

    correct_100_pqd_023 = 0
    correct_100_pqd_049 = 0
    correct_100_pqd_080 = 0

    correct_10_pqd_023 = 0
    correct_10_pqd_049 = 0
    correct_10_pqd_080 = 0

    correct_1_pqd_023 = 0
    correct_1_pqd_049 = 0
    correct_1_pqd_080 = 0

    for model in top100_pqd:
        correct_100_pqd_023, correct_100_pqd_049, correct_100_pqd_080 = check_pickledict(model, correct_100_pqd_023, correct_100_pqd_049, correct_100_pqd_080)

    for model in top10_pqd:
        correct_10_pqd_023, correct_10_pqd_049, correct_10_pqd_080 = check_pickledict(model, correct_10_pqd_023, correct_10_pqd_049, correct_10_pqd_080)


    correct_1_pqd_023, correct_1_pqd_049, correct_1_pqd_080 = check_pickledict(top1_pqd, correct_1_pqd_023, correct_1_pqd_049, correct_1_pqd_080)

    print("------------PQD PERFORMANCE ON TARGET-------------")
    print("\tTop1\tTop 10\tTop 100")
    print("0.23:"+'\t'+ str((correct_1_pqd_023*100))+"%"+'\t'+ str((correct_10_pqd_023/10)*100)+"%"+'\t'+ str(correct_100_pqd_023)+"%")
    print("0.49:"+'\t'+ str((correct_1_pqd_049*100))+"%"+'\t'+ str((correct_10_pqd_049/10)*100)+"%"+'\t'+ str(correct_100_pqd_049)+"%")
    print("0.80:"+'\t'+ str((correct_1_pqd_080*100))+"%"+'\t'+ str((correct_10_pqd_080/10)*100)+"%"+'\t'+ str(correct_100_pqd_080)+"%")
    print('\n')



##########  PQDZ #############
    top100_pqdz = sorted_pqdz[:100] #Generates top lists from the full list that is sorted according to pqdz values
    top10_pqdz = sorted_pqdz[:10]
    top1_pqdz = sorted_pqdz[0]

    correct_100_pqdz_023 = 0
    correct_100_pqdz_049 = 0
    correct_100_pqdz_080 = 0

    correct_10_pqdz_023 = 0
    correct_10_pqdz_049 = 0
    correct_10_pqdz_080 = 0

    correct_1_pqdz_023 = 0
    correct_1_pqdz_049 = 0
    correct_1_pqdz_080 = 0

    for model in top100_pqdz:
        correct_100_pqdz_023, correct_100_pqdz_049, correct_100_pqdz_080 = check_pickledict(model, correct_100_pqdz_023, correct_100_pqdz_049, correct_100_pqdz_080)

    for model in top10_pqdz:
        correct_10_pqdz_023, correct_10_pqdz_049, correct_10_pqdz_080 = check_pickledict(model, correct_10_pqdz_023, correct_10_pqdz_049, correct_10_pqdz_080)


    correct_1_pqdz_023, correct_1_pqdz_049, correct_1_pqdz_080 = check_pickledict(top1_pqdz, correct_1_pqdz_023, correct_1_pqdz_049, correct_1_pqdz_080)

    print("------------PQDZ PERFORMANCE ON TARGET-------------")
    print("\tTop1\tTop 10\tTop 100")
    print("0.23:"+'\t'+ str((correct_1_pqdz_023*100))+"%"+'\t'+ str((correct_10_pqdz_023/10)*100)+"%"+'\t'+ str(correct_100_pqdz_023)+"%")
    print("0.49:"+'\t'+ str((correct_1_pqdz_049*100))+"%"+'\t'+ str((correct_10_pqdz_049/10)*100)+"%"+'\t'+ str(correct_100_pqdz_049)+"%")
    print("0.80:"+'\t'+ str((correct_1_pqdz_080*100))+"%"+'\t'+ str((correct_10_pqdz_080/10)*100)+"%"+'\t'+ str(correct_100_pqdz_080)+"%")
    print('\n')

################## ZRANK ##################

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
        correct_100_zrank_023, correct_100_zrank_049, correct_100_zrank_080 = check_pickledict(model, correct_100_zrank_023, correct_100_zrank_049, correct_100_zrank_080) #walks through top100 and checks how many models are in pickle dicts

    for model in top10_zrank:
        correct_10_zrank_023, correct_10_zrank_049, correct_10_zrank_080 = check_pickledict(model, correct_10_zrank_023, correct_10_zrank_049, correct_10_zrank_080) #walks through top10 and checks how many models are in pickle dicts


    correct_1_zrank_023, correct_1_zrank_049, correct_1_zrank_080 = check_pickledict(top1_zrank, correct_1_zrank_023, correct_1_zrank_049, correct_1_zrank_080)  #takes top1 and checks if it is in pickle dicts




    print("------------ZRANK PERFORMANCE ON TARGET-------------")
    print("\tTop1\tTop 10\tTop 100")
    print("0.23:"+'\t'+ str((correct_1_zrank_023*100))+"%"+'\t'+ str((correct_10_zrank_023/10)*100)+"%"+'\t'+ str(correct_100_zrank_023)+"%")
    print("0.49:"+'\t'+ str((correct_1_zrank_049*100))+"%"+'\t'+ str((correct_10_zrank_049/10)*100)+"%"+'\t'+ str(correct_100_zrank_049)+"%")
    print("0.80:"+'\t'+ str((correct_1_zrank_080*100))+"%"+'\t'+ str((correct_10_zrank_080/10)*100)+"%"+'\t'+ str(correct_100_zrank_080)+"%")
    print('\n')



    ####################################################################################

    pqd_023 = (correct_1_pqd_023, correct_10_pqd_023, correct_100_pqd_023)
    pqd_049 = (correct_1_pqd_049, correct_10_pqd_049, correct_100_pqd_049)
    pqd_080 = (correct_1_pqd_080, correct_10_pqd_080, correct_100_pqd_080)


    pqdz_023 = (correct_1_pqdz_023, correct_10_pqdz_023, correct_100_pqdz_023)
    pqdz_049 = (correct_1_pqdz_049, correct_10_pqdz_049, correct_100_pqdz_049)
    pqdz_080 = (correct_1_pqdz_080, correct_10_pqdz_080, correct_100_pqdz_080)

    zrank_023 = (correct_1_zrank_023, correct_10_zrank_023, correct_100_zrank_023)
    zrank_049 = (correct_1_zrank_049, correct_10_zrank_049, correct_100_zrank_049)
    zrank_080 = (correct_1_zrank_080, correct_10_zrank_080, correct_100_zrank_080)








    return pqd_023, pqd_049, pqd_080, pqdz_023, pqdz_049, pqdz_080, zrank_023, zrank_049, zrank_080






def check_pickledict(model, correct_023, correct_049, correct_080): #checks if model is a correct model according to all three benchmarks
    correct_models_023 = pickle.load(open("/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/propro_matrix/dockq_natives_0.23.p", "rb"))
    correct_models_049 = pickle.load(open("/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/propro_matrix/dockq_natives_0.49.p", "rb"))
    correct_models_080 = pickle.load(open("/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/propro_matrix/dockq_natives_0.8.p", "rb"))

    model_name = (model.split()[0])

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



def plot_performance(pqd_targets023, pqd_targets049, pqd_targets080, pqdz_targets023, pqdz_targets049, pqdz_targets080, zrank_targets023, zrank_targets049, zrank_targets080):

    top1_023 = [pqd_targets023[0], pqdz_targets023[0], zrank_targets023[0], 6, 6]
    top10_023 = [pqd_targets023[1], pqdz_targets023[1], zrank_targets023[1], 26, 24]
    top100_023 = [pqd_targets023[2], pqdz_targets023[2], zrank_targets023[2], 62, 70]

    top1_049 = [pqd_targets049[0], pqdz_targets049[0], zrank_targets049[0], 3, 2]
    top10_049 = [pqd_targets049[1], pqdz_targets049[1], zrank_targets049[1], 12, 12]
    top100_049 = [pqd_targets049[2], pqdz_targets049[2], zrank_targets049[2], 33, 37]

    top1_080 = [pqd_targets080[0], pqdz_targets080[0], zrank_targets080[0], 1, 1]
    top10_080 = [pqd_targets080[1], pqdz_targets080[1], zrank_targets080[1], 3, 2]
    top100_080 = [pqd_targets080[2], pqdz_targets080[2], zrank_targets080[2], 5, 5]




    bars10 = top1_023 + top10_023 + top100_023 + top1_049 + top10_049 + top100_049 + top1_080 + top10_080 + top100_080

    x_list = ["ProQDock", "ProQDockZ", "ZRANK", "Glaser matrix", "True propro matrix"]
    #x_list = ["Scoring matrix", "Glaser scoring matrix"]


    barWidth = 0.08 #0.1

    r1 = list(np.arange(len(top1_023)))
    print(r1)
    r2 = [x + barWidth for x in r1]
    r3 = [x + barWidth for x in r2]
    r4 = [x + barWidth for x in r3]
    r5 = [x + barWidth for x in r4]
    r6 = [x + barWidth for x in r5]
    r7 = [x + barWidth for x in r6]
    r8 = [x + barWidth for x in r7]
    r9 = [x + barWidth for x in r8]

    r10 = r1+r2+r3+r4+r5+r6+r7+r8+r9


    plt.figure(figsize = (9.4, 6.0), dpi = 190) #9.4, 6.0



    plt.bar(r1, top1_023, color='#FECD1C', edgecolor='white', width=barWidth, label = "Acceptable top 1")
    plt.bar(r2, top10_023, color='#FFDC5D', edgecolor='white', width=barWidth, label = "Acceptable top 10")
    plt.bar(r3, top100_023, color='#FFE99A', edgecolor='white', width=barWidth, label = "Acceptable top 100")

    plt.bar(r4, top1_049, color='#FF7D00', edgecolor='white', width=barWidth, label = "Medium top 1")
    plt.bar(r5, top10_049, color='#FF9735', edgecolor='white', width=barWidth, label = "Medium top 10")
    plt.bar(r6, top100_049, color='#FCB775', edgecolor='white', width=barWidth, label = "Medium top 100")

    plt.bar(r7, top1_080, color='#FF62A2', edgecolor='white', width=barWidth, label = "High top 1")
    plt.bar(r8, top10_080, color='#FC84B5', edgecolor='white', width=barWidth, label = "High top 10")
    plt.bar(r9, top100_080, color='#FFC7DE', edgecolor='white', width=barWidth, label = "High top 100")


    plt.ylim((0,134))


    labels = list(top1_023 + top10_023 + top100_023 + top1_049 + top10_049 + top100_049 + top1_080 + top10_080 + top100_080)


    for i in range(len(r10)):
        plt.text(x = r10[i]-0.03, y = bars10[i]+0.5, s = labels[i], size = 6)



    plt.xticks(r5, x_list)

    plt.xlabel("Methods", fontweight = "bold")
    plt.title(f'Performance on 124 targets with 3 DockQ quality benchmarks', fontweight='bold')
    plt.ylabel('Number of targets with >=1 model hits', fontweight = "bold")
    plt.tight_layout()
    plt.legend(fontsize = "small")
    plt.savefig(f"/proj/wallner/users/x_karst/exjobb/pictures/bjorns_performance_CnM.png")
    plt.show()
    plt.clf()



    return



def main():

    args = sys.argv[1:]
    targets_list = find_targets()


    pqd_targets023 = [0, 0, 0]
    pqd_targets049 = [0, 0, 0]
    pqd_targets080 = [0, 0, 0]

    pqdz_targets023 = [0, 0, 0]
    pqdz_targets049 = [0, 0, 0]
    pqdz_targets080 = [0, 0, 0]

    zrank_targets023 = [0, 0, 0]
    zrank_targets049 = [0, 0, 0]
    zrank_targets080 = [0, 0, 0]

    tot_targets = 0

    for target in targets_list:
        tot_targets = tot_targets + 1
        print(target)

        sorted_pqd, sorted_pqdz, sorted_zrank = sort_file(target)
        pqd_023, pqd_049, pqd_080, pqdz_023, pqdz_049, pqdz_080, zrank_023, zrank_049, zrank_080 = calc_performance(sorted_pqd, sorted_pqdz, sorted_zrank)



        if int(pqd_023[0])>0:
            pqd_targets023[0] = pqd_targets023[0] + 1
        if int(pqd_023[1])>0:
            pqd_targets023[1] = pqd_targets023[1] + 1
        if int(pqd_023[2])>0:
            pqd_targets023[2] = pqd_targets023[2] + 1
        if int(pqd_049[0])>0:
            pqd_targets049[0] = pqd_targets049[0] + 1
        if int(pqd_049[1])>0:
            pqd_targets049[1] = pqd_targets049[1] + 1
        if int(pqd_049[2])>0:
            pqd_targets049[2] = pqd_targets049[2] + 1
        if int(pqd_080[0])>0:
            pqd_targets080[0] = pqd_targets080[0] + 1
        if int(pqd_080[1])>0:
            pqd_targets080[1] = pqd_targets080[1] + 1
        if int(pqd_080[2])>0:
            pqd_targets080[2] = pqd_targets080[2] + 1

        if int(pqdz_023[0])>0:
            pqdz_targets023[0] = pqdz_targets023[0] + 1
        if int(pqdz_023[1])>0:
            pqdz_targets023[1] = pqdz_targets023[1] + 1
        if int(pqdz_023[2])>0:
            pqdz_targets023[2] = pqdz_targets023[2] + 1
        if int(pqdz_049[0])>0:
            pqdz_targets049[0] = pqdz_targets049[0] + 1
        if int(pqdz_049[1])>0:
            pqdz_targets049[1] = pqdz_targets049[1] + 1
        if int(pqdz_049[2])>0:
            pqdz_targets049[2] = pqdz_targets049[2] + 1
        if int(pqdz_080[0])>0:
            pqdz_targets080[0] = pqdz_targets080[0] + 1
        if int(pqdz_080[1])>0:
            pqdz_targets080[1] = pqdz_targets080[1] + 1
        if int(pqdz_080[2])>0:
            pqdz_targets080[2] = pqdz_targets080[2] + 1

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







    print('\n\n\n\n\n')
    print("-------PQD TOTAL PERFORMANCE OUT OF", tot_targets, "-------")
    print("0.23 <top1, top10, top100>:", pqd_targets023)
    print("0.49 <top1, top10, top100>:", pqd_targets049)
    print("0.80 <top1, top10, top100>:", pqd_targets080)

    print("-------PQDZ TOTAL PERFORMANCE OUT OF", tot_targets, "-------")
    print("0.23 <top1, top10, top100>:", pqdz_targets023)
    print("0.49 <top1, top10, top100>:", pqdz_targets049)
    print("0.80 <top1, top10, top100>:", pqdz_targets080)

    print("-------ZRANK TOTAL PERFORMANCE OUT OF", tot_targets, "-------")
    print("0.23 <top1, top10, top100>:", zrank_targets023)
    print("0.49 <top1, top10, top100>:", zrank_targets049)
    print("0.80 <top1, top10, top100>:", zrank_targets080)






    plot_performance(pqd_targets023, pqd_targets049, pqd_targets080, pqdz_targets023, pqdz_targets049, pqdz_targets080, zrank_targets023, zrank_targets049, zrank_targets080)


















if __name__ == '__main__':
    main()
