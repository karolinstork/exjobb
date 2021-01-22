#!/usr/bin/env python
#SBATCH -t 7-00:00:00
#SBATCH -N 1

import sys
import os
import Bio

import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)



def check_length(chain_1, chain_2, cut_off):
    protein_peptide = None
    protein_protein = None
    peptide_peptide = None
    receptor = None

    from Bio.PDB.PDBParser import PDBParser
    parser = PDBParser(PERMISSIVE=1)

    file_extension = chain_1.split(".")[-1]
    pdb_file = chain_1[0:4]+"."+file_extension
    dir= pdb_file[1:3]
    path = f"/proj/wallner/share/PDB/191015_biounit/{dir}/{pdb_file}"
    print(path)
    structure = parser.get_structure(pdb_file, path)
    first_number_residues = 0
    second_number_residues = 0


    ####################################################
    chain1_name = chain_1[4]
    piece1 = chain_1.split(".")[0]
    model_for_chain1 = piece1[5:]

    chain2_name = chain_2[4]
    piece2 = chain_2.split(".")[0]
    model_for_chain2 = piece2[5:]

    model_for_chain1 = int(model_for_chain1) - 1 #convert for python format
    model_for_chain2 = int(model_for_chain2) -1

    first_model = structure[model_for_chain1]
    second_model = structure[model_for_chain2]
    first_chain = first_model[chain1_name]
    second_chain = second_model[chain2_name]

    first_residues = first_chain.get_residues()
    for residue in first_residues:
        if Bio.PDB.is_aa(residue):
            first_number_residues = first_number_residues + 1


    second_residues = second_chain.get_residues()
    for residue in second_residues:
        if Bio.PDB.is_aa(residue): #atom records must only have amino acids, cut off applies to atom records
            second_number_residues = second_number_residues + 1


    print(chain_1, "Length:", first_number_residues)
    print(chain_2, "Length:", second_number_residues)

    if (first_number_residues>cut_off) and (second_number_residues>cut_off):
        protein_protein = True

    elif (first_number_residues <= cut_off) and (second_number_residues <=cut_off):
        peptide_peptide = True

    elif ((first_number_residues>cut_off) and (second_number_residues<=cut_off)) or ((first_number_residues<=cut_off) and (second_number_residues>cut_off)):
        protein_peptide = True
        if first_number_residues<=cut_off:
            receptor = chain_2
        elif second_number_residues<=cut_off:
            receptor = chain_1

    else:
        print("ERROR")
        return




    return protein_peptide, protein_protein, peptide_peptide, receptor









def main():
    try:
        args = sys.argv[1:]
        complexes_file = args[0]
        cut_off = args[1]
        cut_off = int(cut_off)

    except IndexError:
        print("INPUT ERROR: <categorize_protein_peptide_data.py> <file with cluster data>")
        return

    protein_peptide_list = open("/proj/wallner/users/x_karst/exjobb/protein_peptide_data/protein_peptide_list_part2.txt", "w")
    protein_protein_list = open("/proj/wallner/users/x_karst/exjobb/protein_peptide_data/protein_protein_list_part2.txt", "w")
    peptide_peptide_list = open("/proj/wallner/users/x_karst/exjobb/protein_peptide_data/peptide_peptide_list_part2.txt", "w")

    complexes_file = open(complexes_file, "r")
    lines_of_complexes = complexes_file.readlines()
    print(len(lines_of_complexes))

    for complex_data in lines_of_complexes:
        print("--------------------------------------------------------------------------")
        chain_1 = complex_data.split()[0]
        chain_2 = complex_data.split()[1]

        protein_peptide, protein_protein, peptide_peptide, receptor = check_length(chain_1, chain_2, cut_off)
        print("PROTEIN PROTEIN:", protein_protein)
        print("PEPTIDE PEPTIDE:", peptide_peptide)
        print("PROTEIN PEPTIDE:", protein_peptide)
        print("Receptor:", receptor)

        if protein_peptide == True:
            data = complex_data.strip("\n")+"\t"+str(receptor)+"\n"
            print(data)
            protein_peptide_list.write(data)

        elif protein_protein == True:
            data = complex_data
            print("PROTEIN PROTEIN CONFIRMED")
            print(data)
            protein_protein_list.write(data)

        elif peptide_peptide == True:
            data = complex_data
            print("PEPTIDE PEPTIDE CONFIRMED")
            print(data)
            peptide_peptide_list.write(data)







    protein_peptide_list.close()
    protein_protein_list.close()
    peptide_peptide_list.close()






if __name__ == '__main__':
    main()
