import os
import sys



def main():

    args = sys.argv[1:]
    csv_file = args[0]
    dockq_benchmark = float(args[1])
    top = int(args[2])


    numb_success_predictions = 0
    f = open(csv_file, "r")
    lines = f.readlines()
    lines = lines[1:]



    for line in lines: #number of hits for that target
    #    print(line)
        if dockq_benchmark == 0.23:
            if top == 1:
                hit = line.split(",")[1]
            if top == 10:
                hit = line.split(",")[2]
            if top == 100:
                hit = line.split(",")[3]

        if dockq_benchmark == 0.49:
            if top == 1:
                hit = line.split(",")[4]
            if top == 10:
                hit = line.split(",")[5]
            if top == 100:
                hit = line.split(",")[6]

        if dockq_benchmark == 0.80:
            if top == 1:
                hit = line.split(",")[7]
            if top == 10:
                hit = line.split(",")[8]
            if top == 100:
                hit = line.split(",")[9]



        if int(hit)>0:
            print("SUCCESS FOR", line)
            numb_success_predictions = numb_success_predictions + 1



    print("Number of successful predictions for", dockq_benchmark, "top"+str(top)+":", numb_success_predictions)





    f.close()







if __name__ == '__main__':
    main()
