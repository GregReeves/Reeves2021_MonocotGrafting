import Supermatrix
import sys, os

if len(sys.argv) != 3:
	print sys.argv[0] + " super.matrix super.model"
	sys.exit(0)

file1 = open(sys.argv[1],"r")
file2 = open(sys.argv[2],"r")
Super = Supermatrix.Alignment()
	
Super.read_fasta_as_hash(file1)
Super.get_part_info(file2)
#Supermatrix.missing_check()
	
Genes = Supermatrix.GeneSet()
Genes = Super.part_to_genes()
Genes.print_all()
