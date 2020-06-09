import os
import glob

def main():

    files =  glob.glob("/proj/wallner/users/x_karst/exjobb/results_naccess/*")

    big_file = open("naccess_results_0906.txt", "w")

    for file in files:
        print(file)
        text = open(file, "r")
        text = text.read()
        big_file.write(text)


    print("combined: ", len(files), "files")



if __name__ == '__main__':
    main()
