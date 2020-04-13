#!/bin/bash
#SBATCH -N 1
#SBATCH -t 1440
#SBATCH snic-2019-9-224

/proj/wallner/users/x_bjowa/local/blast-2.2.26/bin/./blastclust -a 32 -i fasta_file.txt -S 30 -L 0.5 -o blastclust.out
