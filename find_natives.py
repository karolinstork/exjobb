
import os
import sys
import pickle





def main():
    try:
        args = sys.argv[1:]
        cut_off = float(args[0]) #dockq cut offs
        which_set = args[1] #if pickle should be saved in propro or propep
    except IndexError:
        print("INPUT ERROR")
        print("USAGE: find_natives.py <dockq cutoff> [propep] [propro]")
        return

    if which_set == "propep":
        scoring_file = open("/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/pro_pep/dockq_results.txt", "r")
        scoring_file = scoring_file.readlines()
        pickle_dict = {}

        for line in scoring_file:
            if line[0] != "<":
                line_split = line.split()
                filename = line_split[0]
                dockq = line_split[1]

                if float(dockq) >= cut_off:
                    if float(dockq)< 1:
                        pickle_dict[filename] = dockq
                        print(filename, dockq)


        pickle.dump(pickle_dict, open(f"/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/pro_pep/dockq_natives_{cut_off}.p", "wb"))
        print("Pickle saved in /proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/pro_pep/")
        print("Threshold:", cut_off)




    if which_set == "propro":
        scoring_file = open("/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/pro_pro/bjorns_scores.txt", "r")
        scoring_file = scoring_file.readlines()
        pickle_dict = {}

        for line in scoring_file:
            if line[0] != "#":
                line_split = line.split()
                filename = line_split[0]
                dockq = line_split[16]

                if float(dockq) >= cut_off:
                    if float(dockq)< 1: #some files had higher than 1 which should not be possible
                        pickle_dict[filename] = dockq
                        print(filename, dockq)


        pickle.dump(pickle_dict, open(f"/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/pro_pro/dockq_natives_{cut_off}.p", "wb"))
        print("Pickle saved in /proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/pro_pro/")
        print("Threshold:", cut_off)








if __name__ == '__main__':
    main()
