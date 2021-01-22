import os
import sys
import pickle





def main():

    args = sys.argv[1:]
    try:
        score_file = args[0]
        whos_result = args[1] #my bjorns or isaks
        which_set = args[2] #propro or propep

    except IndexError:
        print("INPUT ERROR")
        print("USAGE: compare_scoring_to_dockq.py <score_file> [my] [isaks] [bjorns] [propro] [propep]")
        return

    scoring_file = open(score_file, "r")
    scoring_file = scoring_file.readlines()

    if which_set == "propro":
        correct_models_dict = pickle.load(open("/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/pro_pro/dockq_natives.p", "rb"))
    if which_set == "propep":
        correct_models_dict = pickle.load(open("/proj/wallner/users/x_karst/exjobb/evaluate_scoring_matrix/pro_pep/dockq_natives.p", "rb"))



    tot_numb_models = 0
    tot_correct_models = 0.0


    if whos_result == "my" and which_set == "propro": #my filenames are not in same format as propro facit
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


                if filename in correct_models_dict.keys():
                    print(filename, "MATCH!")
                    tot_correct_models = tot_correct_models + 1
                else:
                    print(filename)

            else: #for files from the other directory
                filename = line.split()[0]
                first_part = filename.split("-")[0]
                second_part = filename.split("-")[1]
                filename = first_part+"-"+second_part

                if filename in correct_models_dict.keys():
                    print(filename, "MATCH!")
                    tot_correct_models = tot_correct_models + 1
                else:
                    print(filename)


    if whos_result == "my" and which_set == "propep": #in the propep set, my filenames has the same name in both facit and score file
        for line in scoring_file:
            if line[0] != "<":
                tot_numb_models = tot_numb_models + 1
                filename = line.split()[0]

                if filename in correct_models_dict.keys():
                    print(filename, "MATCH!")
                    tot_correct_models = tot_correct_models + 1
                else:
                    print(filename)




    if whos_result == "bjorns": #for bjorns propro results, the files has the same name in facit and score file
        for line in scoring_file:
            tot_numb_models = tot_numb_models + 1
            filename = line.split()[0]

            if filename in correct_models_dict.keys():
                print(filename, "MATCH!")
                tot_correct_models = tot_correct_models + 1
            else:
                print(filename)



    if whos_result == "isaks": #results from isaks zrank for propep
        for line in scoring_file:
            tot_numb_models = tot_numb_models + 1
            full_name = line.split()[0]
            first_name = full_name.split("_")[0]
            first_name = first_name[:4]
            numb = full_name.split("_")[1]
            name = first_name+"_"+numb

            if name in correct_models_dict.keys():
                print(name, "MATCH!")
                tot_correct_models = tot_correct_models + 1
            else:
                print(name)




    percent = (tot_correct_models / tot_numb_models) *100
    print("----------------------------")
    print(tot_correct_models, "/",tot_numb_models)
    print("Correctness:", percent, "%")



















if __name__ == '__main__':
    main()
