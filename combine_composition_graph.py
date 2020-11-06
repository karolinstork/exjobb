import sys
import os
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr






def main():
  print("Plotting interface composition . . .")
  x = ["L", "A", "G", "V", "S", "I", "T", "E", "R", "P", "D", "N", "K", "F", "Y", "Q", "H", "M", "W", "C"]
  y_propro = [9.96, 8.36, 7.31, 7.29, 6.93, 5.99, 5.92, 5.49, 5.05, 5.04, 4.74, 4.48, 4.33, 4.25, 3.91, 3.79, 2.45, 2.15, 1.36, 1.21]
  y_propep = [9.81, 7.55, 6.48, 7.01, 6.50, 5.82, 5.62, 5.11, 4.98, 5.15, 4.61, 4.88, 5.04, 4.76, 4.35, 3.42, 2.42, 2.56, 2.32, 1.61]


  corr, p = pearsonr(y_propro, y_propep)
  print("Correlation:", corr)

  diff =[]

  # for element1,element2 in zip(y_propro,y_propep):
  #      diff.append(abs(element1 - element2))
  # print(diff)


  plt.scatter(x, y_propro, s = 20, color = "C1")
  plt.plot(x, y_propro, label = "Protein protein (16286 complexes)", color = "C1")

  plt.scatter(x, y_propep, s = 20, color= "C0")
  plt.plot(x, y_propep, label = "Protein peptide (1397 complexes)", color = "C0")
  plt.xlabel("Residues")
  plt.ylabel("Wi [%]")
  plt.ylim(0,12)
  plt.title("Interface composition")
  plt.legend()

  plt.savefig("/proj/wallner/users/x_karst/exjobb/pictures/interface_composition_comparison.png", bbox_inches='tight')
  plt.show()


#############################################################################################################

  print("Plotting leucine pairing patterns . . .")
  x = ["R", "K", "N", "D", "Q", "E", "H", "P", "Y", "W", "S", "T", "G", "A", "M", "C", "F", "L", "V", "I"]
  y_glaser = [4.99, 3.15, 2.31, 1.40, 3.46, 3.12, 4.88, 2.50, 4.19, 5.77, 1.41, 2.07, -0.37, 2.77, 5.32, 2.93, 4.86, 4.03, 4.20, 4.59]
  y_my = [8.93, 8.28, 3.79, 2.68, 7.46, 6.48, 7.63, 6.40, 11.04, 13.02, 3.51, 4.64, -1.70, 8.04, 10.36, 5.53, 12.50, 9.02, 9.84, 11.25]

  corr, p = pearsonr(y_glaser, y_my)
  print("Correlation:", corr)

  plt.scatter(x, y_my, s = 20)
  plt.plot(x, y_my, label = "Thesis results")
  plt.scatter(x, y_glaser, s = 20)
  plt.plot(x, y_glaser, label = "Glaser et al results")
  plt.xlabel("Partner residue")
  plt.ylabel("Likelihood")
  plt.ylim(-5,15)
  plt.title("Leucine pairing patterns")
  legend = plt.legend(loc='best')
  plt.axhline(y=0, color = "grey")

  plt.savefig("/proj/wallner/users/x_karst/exjobb/pictures/leucine_pairings.png", bbox_inches='tight')
  plt.show()




if __name__ == '__main__':
    main()
