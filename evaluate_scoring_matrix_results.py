import os
import sys
import pickle





def main():

    args = sys.argv[1:]
    score_file = args[0]
    scoring_file = open(score_file, "r")
    scoring_file = scoring_file.readlines()
    correct_models_dict = pickle.load(open("/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/pro_pro/m1_natives.p", "rb"))


    tot_numb_models = 0
    tot_correct_models = 0.0


    if "bjorns" not in score_file:
        for line in scoring_file:
            tot_numb_models = tot_numb_models + 1
            if line[:6] == "Target": #for the files from target directories
                score = line.split()[1]
                filename = line.split()[0]

                target_numb = filename.split("-")[0]
                first_numb = target_numb[6:]

                model_numb = filename.split("-")[1]
                model_numb = model_numb.split(".")[0] #remvoing file extension
                second_numb = model_numb[5:]

                filename = "T"+str(first_numb)+"-"+str(second_numb)
                print(filename)

                if filename in correct_models_dict.keys():
                    print("MATCH!")
                    tot_correct_models = tot_correct_models + 1

            else: #for files from the other directory
                filename = line.split()[0]
                first_part = filename.split("-")[0]
                second_part = filename.split("-")[1]
                filename = first_part+"-"+second_part
                print(filename)
                if filename in correct_models_dict.keys():
                    print("MATCH!")
                    tot_correct_models = tot_correct_models + 1

    else: #
        for line in scoring_file:
            tot_numb_models = tot_numb_models + 1
            filename = line.split()[0]
            print(filename)
            if filename in correct_models_dict.keys():
                print("MATCH!")
                tot_correct_models = tot_correct_models + 1

    percent = (tot_correct_models / tot_numb_models) *100
    print("----------------------------")
    print("Correctness:", percent, "%")



















if __name__ == '__main__':
    main()
