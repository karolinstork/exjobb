import sys
import os







def main():
  args = sys.argv[1:]
  cluster_info_file = args[0]
  new_file = open("cdhit/clusters.txt", "w")

  clusterinfo = open(cluster_info_file, "r")
  lines = clusterinfo.readlines()



  for line in lines:
      line = line.strip('\n')
      words = line.split(',')
      print(words)
      file = words[0]
      cluster_thirtypercent = words[-1]
      string = str(file)+'\t'+str(cluster_thirtypercent)+'\n'
      print(string)
      new_file.write(string)




  new_file.close()
  clusterinfo.close()



if __name__ == '__main__':
    main()
