import pickle
import os
import sys
import matplotlib.pyplot as plt





def main():

    dataset = open("/proj/wallner/users/x_karst/exjobb/dataset_whole.out", "r")
    lines = dataset.readlines()

    resolution_dict = pickle.load(open("/proj/wallner/users/x_karst/exjobb/pickles/resolutions_sbatch.p", "rb"))
    not_in_dict = []
    x = []
    y = []

    key = 0
    for line in lines:
        id = line.split(".")[0] #3vskA1.pdb1
        id = id[:4]
        id = str(id) + ".pdb"
        try:
            resolution = resolution_dict[id]
            if resolution != None:
                key = key + 1
                print(id, resolution)
                x.append(resolution)
                y.append(key)


        except KeyError:
            not_in_dict.append(id)


    print("NOT IN RESOLUTION DICT:")
    print(not_in_dict)


    plt.title("Resolution distribution in data set")

    print(len(x))
    plt.scatter(x, y, s=1)
    plt.show()














if __name__ == '__main__':
    main()
