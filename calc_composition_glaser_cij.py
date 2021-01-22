#!/usr/bin/env python
#SBATCH -t 7-00:00:00
#SBATCH -N 1




import os
import sys



def main():


    #dict_res = {'I': 1697, 'V':2370, 'L': 2258, 'F': 1337, 'C':623, 'M': 793, 'A': 3056, 'G': 2860, 'T':2289 , 'S':2333 , 'W':544 , 'Y': 1234, 'P': 2243, 'H': 854 , 'E':1568, 'Q':1137 , 'D': 1686, 'N': 1630, 'K':1242 , 'R': 1485}
    dict = {'I': 1775, 'V':2549, 'L': 2438, 'F': 1399, 'C': 666, 'M': 821, 'A': 3284, 'G': 3052, 'T':2412 , 'S':2480 , 'W':551 , 'Y':1289, 'P': 2340, 'H': 882, 'E':1624, 'Q': 1173, 'D':1741, 'N': 1723, 'K': 1278 , 'R': 1523}


    tot_numb_of_res = sum(dict.values())
    tot_numb_of_res = float(tot_numb_of_res)
    print(tot_numb_of_res)

    plot_list = []
    x_list = []
    y_list = []

    for res in dict.keys():
        print(dict[res])
        ratio = dict[res]/tot_numb_of_res
        ratio = ratio * 100 #percent
        plot_list.append((res, ratio))

    plot_list = sorted(plot_list, key = lambda x:x[1], reverse = True)
    for item in plot_list:
        x_list.append(item[0])
        y_list.append(item[1])

    print(x_list)
    print(y_list)

    import matplotlib.pyplot as plt
    plt.scatter(x_list, y_list)
    plt.plot(x_list, y_list)
    plt.xlabel("Residues")
    plt.ylabel("Wi [%]")
    plt.ylim(0,12)
    title = "Interface composition Glaser et al"
    plt.title(title)
    plt.show()






























if __name__ == '__main__':
    main()
