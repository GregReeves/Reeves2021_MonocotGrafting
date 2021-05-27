import re,sys
import random

#An object to store sets of alignment objects
class GeneSet:
	def __init__(self):
		
		self.gene_set = []

	def print_all(self):
		
		for seq in self.gene_set:
			seq.print_to_file(seq.name)
		



# An object designed to get info about sequences
class Alignment:

	def __init__(self):
		
		self.Supermatrix = {}
		self.length = 0
		self.amino_acid = False
		self.partition_info = [] #can do an array of tuples
		self.name = ""  
		self.taxa = 0 #total taxa that exist
	
		
	#takes in a file and turns it into a hash
	def read_fasta_as_hash(self,file1):
		
		for line in file1:
			line = line.strip("\n\r")
			if line[0] == ">":
				seq_name = line[1:]
				self.Supermatrix[seq_name] = ""
			else:
				self.Supermatrix[seq_name] += line
		
		self.length = len(self.Supermatrix[seq_name])
		
		ambigs = "EFHILPQ" #All the IUPAC letters that don't double as ambiguous for nucleotide
		
		if any(c in self.Supermatrix[seq_name] for c in ambigs):
			self.amino_acid = True
		
		self.taxa = len(self.Supermatrix)
	
	#Gets the partition information
	def get_part_info(self,partition_file):
		
		for line in partition_file:
			
			m = re.search(".*?,(.*?)=(.*?)-(.*)",line.strip("\n\r").replace(" ",""))
			tup = (m.group(1),m.group(2),m.group(3))
			self.partition_info.append(tup)
			
	
	#Jackknife some number of sites and return the results in a HASH
	def jackknife(self,sites):
		
		knife = set()
		rep = Alignment()
		rep.amino_acid = self.amino_acid
		rep.length = sites
		
		if self.length < sites:
			 print "Specified more sites than length"
			 sys.exit()
		
		while len(knife) < sites:
			#correct for sites starting at 0
			a = (random.randint(1,self.length) - 1)
			knife.add(a)
		
		for i in knife:
			for seq in self.Supermatrix:
				
				try: 
					rep.Supermatrix[seq] += self.Supermatrix[seq][i]
				except:
					rep.Supermatrix[seq] = ""
					rep.Supermatrix[seq] += self.Supermatrix[seq][i]
		
		rep.taxa = len(rep.Supermatrix)
		return rep
	
	#Check if taxa has just missing data
	def MissingCheck(self):
	
		miss_array = []
		for i in self.Supermatrix:
			
			#this is the number of sites that are not missing data
			InfoData = 0 

			if self.amino_acid == False:
				InfoData = (self.length - (self.Supermatrix[i].count("N") + self.Supermatrix[i].count("-") + self.Supermatrix[i].count("?")))
			else:
				InfoData = (self.length - (self.Supermatrix[i].count("X") + self.Supermatrix[i].count("-") + self.Supermatrix[i].count("?")))
			if InfoData == 0:
				#can't change a map during an iterative procedure so story in an array the
				#values to be removed
				miss_array.append(i)
		
		for i in miss_array:
			self.Supermatrix.pop(i)

		self.taxa = len(self.Supermatrix)
	
	#Prints the results to a file on the command line				
	def print_to_file(self,out_name):
		
		outw = open(out_name + ".fa", "w")
		for x in self.Supermatrix:
			outw.write(">"+x+"\n"+self.Supermatrix[x]+"\n")
		outw.close()

	#remove taxa from partitions and move to a GeneSet file
	def part_to_genes(self):
		
		if len(self.partition_info) == 0:
			print "No partition file available"
			sys.exit()
		
		Genes = GeneSet()
		for part in self.partition_info:

			Aln = Alignment()
			Aln.name = part[0]
			for taxa in self.Supermatrix:
				Aln.Supermatrix[taxa] = self.Supermatrix[taxa][(int(part[1])-1):(int(part[2]))]
			Aln.length = ((int(part[2]) - int(part[1])+1))

			Aln.MissingCheck()
			Genes.gene_set.append(Aln)
		
		return Genes



if __name__ == "__main__":
	
	#test_string = ">Fasta1\nATCG\n>Fasta2\nATTT\nFasta3\nCATA\nFasta4\nA--T"
	#print test_string
	file1 = open(sys.argv[1], "r")
	Matrix = Alignment()
	Matrix.name = sys.argv[1]
	Matrix.read_fasta_as_hash(file1)
	file2 = open(sys.argv[2], "r")
	Matrix.get_part_info(file2)
	Replicate = Alignment()
	Replicate = Matrix.jackknife(50)
	Replicate.MissingCheck()
	Replicate.print_to_file("rep")
	
	Genes = GeneSet()
	Genes = Matrix.part_to_genes()
	print Genes.gene_set
	Genes.print_all()
	

	

	
