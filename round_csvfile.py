import sys
import os
import pandas as pd






def round_data(csv_file):

    df = pd.read_csv(csv_file)
    print(df)

    rounded_df = df.round(decimals=2)
    print(rounded_df)

    title = (os.path.basename(csv_file))
    full_path = "/proj/wallner/users/x_karst/exjobb/rounded_csv/"+ str(title)
    rounded_df.to_csv(full_path, index = False)
    return








def main():

    args = sys.argv[1:]
    csv_file = args[0]


    round_data(csv_file)



















if __name__ == '__main__':
    main()
