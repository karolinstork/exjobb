import os
import sys



def main():

    args = sys.argv[1:]
    csv_file = args[0]
    dockq_benchmark = float(args[1])



    hits_1 = 0
    hits_10 = 0
    hits_100 = 0



    f = open(csv_file, "r")
    lines = f.readlines()
    lines = lines[1:]



    for line in lines: #number of hits for that target
        if dockq_benchmark == 0.23:
            hit_1 = line.split(",")[1]
            hit_10 = line.split(",")[2]
            hit_100 = line.split(",")[3]

        if dockq_benchmark == 0.49:
            hit_1 = line.split(",")[4]
            hit_10 = line.split(",")[5]
            hit_100 = line.split(",")[6]

        if dockq_benchmark == 0.80:
            hit_1 = line.split(",")[7]
            hit_10 = line.split(",")[8]
            hit_100 = line.split(",")[9]



        if int(hit_1)>0:
            #print("SUCCESS FOR TOP1", line)
            hits_1 = hits_1 + 1

        if int(hit_10)>0:
            #print("SUCCESS FOR TOP10", line)
            hits_10 = hits_10 + 1


        if int(hit_100)>0:
            #print("SUCCESS FOR TOP100", line)
            hits_100 = hits_100 + 1


    print("Number of successful predictions for", dockq_benchmark)
    print("Top 1:", hits_1, "   Top 10:", hits_10, "    Top 100:", hits_100)





    f.close()







if __name__ == '__main__':
    main()
