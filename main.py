#import the requierd packages
from Bio.Seq import Seq
from Bio import Entrez
from Bio import SeqIO
import argparse
import sys
import os

#implement a function to be used by the user to choose either running the whole data
#or only the test files that are provided
def main(args = None):
	parser = argparse.ArgumentParser(description='perform HCMV analysis')
	parser.add_argument('-test',action='store_true',help='use -test if you want to run the program using testdata')
	return parser.parse_args(args)

#check command line
args = main(sys.argv[1:])
test=args.test

#use Entrez here:
Entrez.email = "faltunusi1@luc.edu"

#get to directory 
dirr = os.getcwd()

#store the samll dataset
dataset = dirr+"/dataset"

#data IDs
SRR = ['SRR5660030','SRR5660033','SRR5660044','SRR5660045']

#if -test flag was not used by user, the whole data will be downloaded 
if (test == False): 
	os.system("mkdir "+dirr+"/miniProject_Feras_Altunusi")
	outdir = dirr+"/miniProject_Feras_Altunusi"
	outfile = open(outdir+"/miniProject.log","w")

	os.chdir("dataset")
	path = {} 
	for s in SRR:
		os.system("wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/"+s+"/"+s+".1")
		fastqDump = "fastq-dump --split-files "+s+".1 -O "+dataset
		os.system(fastqDump)
		path[s] = [dataset+"/"+s+".1_1.fastq",dataset+"/"+s+".1_2.fastq"]

        #go back to dirr 
	os.chdir(dirr)

#otherwise, stored dataset will be used in this case
else:
	os.system("mkdir "+dirr+"/test_outputs")
	outdir = dirr+"/test_outputs"
	
	#Opening the output file
	outfile = open(outdir+"/miniProject.log","w")

	path = {}
	for s in SRR:
        	path[s] = [dataset+"/"+s+".1_1.fastq",dataset+"/"+s+".1_2.fastq"]

#----------------------------------------------------------------
        	
#use kallisto in this first case
os.chdir("dataset")

#initialize the CDS
CDS = 0 
handle = Entrez.efetch(db="nucleotide",id="EF999921",rettype="gb",retmode="text")

#read in genbank entry
record = SeqIO.read(handle,"genbank")
handle.close()

#get the path to transcriptome file
transPath = dataset+"/HCMV_transcriptome.fa"
transcriptome = open(transPath,"w")

#use for loop to go over this list:
for rec in record.features: #iterate through list of features
        if rec.type == "CDS":
                CDS +=1

                #get the secuence IDs
                transcriptome.write(">"+str(record.id)+"\n")

                #get the sequence!
                transcriptome.write(str(rec.extract(record.seq))+"\n")

#close the transcriptome file
transcriptome.close()

#Print number of CDS to outfile file
outfile.write("The HCMV genome (EF999921) has "+str(CDS)+" CDS.\n")

#------------------------------------------------------

#second case, reference genome for bowtie2
handle = Entrez.efetch(db="nucleotide", id="EF999921", rettype="fasta", retmode="text")
record=SeqIO.read(handle,"fasta")
handle.close()

#get the genome file:
genPath = dataset+"/HCMV_genome.fa" 
genome = open(genPath,"w")
genome.write(">"+str(record.id)+"\n")
genome.write(str(record.seq)+"\n")
genome.close()

#get the betaherpes file
os.system("gunzip sequence.fasta.gz")
sub = dataset+"/sequence.fasta"

#done for now, go back to main directory
os.chdir(dirr) 

#third case, build kallisto index
os.system("mkdir "+outdir+"/kallisto") 
kallistoCommand = "kallisto index -i "+outdir+"/kallisto/HCMV_transcriptome.idx --make-unique "+transPath
os.system(kallistoCommand)

#use deafult k-mer size which is 31
sleuthIn = open(outdir+"/kallisto/sample_table.tsv","w") 
sleuthIn.write("sample\tdays post-infection\tpath\n")

#use for loop to go over kallisto quantification
for s in SRR:
        kallistoQuant = "kallisto quant -i "+outdir+"/kallisto/HCMV_transcriptome.idx -o "+outdir+"/kallisto/"+s+" -b 30 -t 2 "+path[s][0]+" "+path[s][1]
        os.system(kallistoQuant)
        if s == SRR[0] or s == SRR[2]: #1 and 3 are 2 d
                sleuthIn.write(s+"\t2\t"+outdir+"/kallisto/"+s+"\n")
        else: #others are 6 d
                sleuthIn.write(s+"\t6\t"+outdir+"/kallisto/"+s+"\n")
sleuthIn.close()
outfile.close()

#--------------------------------------------------------------------------------------------

#analyze the kallisto quant output with sleuth in R
os.system("Rscript sleuth.R "+str(test)) 
outfile = open(outdir+"/miniProject.log","a")

#fourth case, create bowtie2 index 
os.system("mkdir "+outdir+"/bowtie")
bowtie2 = "bowtie2-build "+genPath+" "+outdir+"/bowtie/HCMV"
os.system(bowtie2)

#aligne transcriptomes to HCMV
for s in SRR:
        
        #get number of read pairs before filtering
	before = len(open(path[s][0],"r").readlines())/4 

	#Perform bowtie2 alignment and store the data 
	bowtie1 = "bowtie2 --quiet -x "+outdir+"/bowtie/HCMV -1 "+path[s][0]+" -2 "+path[s][1]+" -S "+outdir+"/bowtie/HCMVmap.sam --al-conc "+outdir+"/bowtie/"+s+".1_mapped_%.fq"
	os.system(bowtie1)

	#add the reads to the paths dictionary
	path[s].append(outdir+"/bowtie/"+s+".1_mapped_1.fq") #at 2
	path[s].append(outdir+"/bowtie/"+s+".1_mapped_2.fq") #at 3

	#get number of read pairs after filtering
	after = len(open(path[s][2],"r").readlines())/4

        #writing to output file
	if s == SRR[0] or s == SRR[1]: #1 and 2 are from donor 1
		if s == SRR[0]: #first one is 2dpi
			outfile.write("Donor 1 (2dpi) had "+str(before)+" read pairs before Bowtie2 filtering and "+str(after)+" read pairs after.\n")
		else: #2 is 6dpi
			outfile.write("Donor 2 (6dpi) had "+str(before)+" read pairs before Bowtie2 filtering and "+str(after)+" read pairs after.\n")
	else: #3 and 4 are from donor 3
		if s == SRR[2]: #3 is 2dpi
			outfile.write("Donor 3 (2dpi) had "+str(before)+" read pairs before Bowtie2 filtering and "+str(after)+" read pairs after.\n")
		else: #4 is 6dpi
			outfile.write("Donor 4 (6dpi) had "+str(before)+" read pairs before Bowtie2 filtering and "+str(after)+" read pairs after.\n")

#assmbly the four transcriptomes with SPAdes
os.system("mkdir "+outdir+"/SPAdes") 
SPAdesC = "spades -k 55,77,99,127 -t 2 --only-assembler --pe1-1 "+path[SRR[0]][2]+" --pe1-2 "+path[SRR[0]][3]+" --pe2-1 "+path[SRR[1]][2]+" --pe2-2 "+path[SRR[1]][3]+" --pe3-1 "+path[SRR[2]][2]+" --pe3-2 "+path[SRR[2]][3]+" --pe4-1 "+path[SRR[3]][2]+" --pe4-2 "+path[SRR[3]][3]+" -o "+outdir+"/SPAdes/HCMV_assembly"
os.system(SPAdesC)
outfile.write(SPAdesC+"\n")

#evaluate the assembly with reads 
contigs = open(outdir+"/SPAdes/HCMV_assembly/contigs.fasta")

#try to get the longest one
longestContig = Seq("") 
assembly_length = 0

#counting contigs > 1000bp
long_contigs = 0

#Parsing the contigs file
records = SeqIO.parse(contigs,"fasta")

#go over the records
for record in records:
	assembly_length = assembly_length+len(record.seq) #add sequence length to total
	if len(record.seq) > 1000:
		long_contigs+=1
	if len(record.seq) > len(longestContig):
		longestContig=record.seq
outfile.write("There are "+str(long_contigs)+" contigs > 1000 bp in the assembly.\n")
outfile.write("There are "+str(assembly_length)+" bp in the assembly.\n")
contigs.close()

#get query file for blast with the longest contig
os.system("mkdir "+outdir+"/BLAST")
queryFile = open(outdir+"/BLAST/query.fasta","w")
queryFile.write(">HCMV\n")
queryFile.write(str(longestContig)+"\n")
queryFile.close()

#---------------------------------------------------------

#Creating the BLAST db
make_blast_command = "makeblastdb -in "+sub+" -out "+outdir+"/BLAST/Betaherpesvirinae  -title Betaherpesvirinae  -dbtype nucl"
os.system(make_blast_command)

#Running BLAST
blastC = 'blastn -query '+outdir+'/BLAST/query.fasta -db '+outdir+'/BLAST/Betaherpesvirinae  -out '+outdir+'/BLAST/myresults.csv -outfmt "10 sacc pident length qstart qend sstart send bitscore evalue stitle"'
os.system(blastC)

#Parsing the BLAST output
outfile.write("sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle\n")
blast_output = open(outdir+"/BLAST/myresults.csv","r")
lines = blast_output.readlines()

#get the top 10 hits from the top of the file
for i in range(0,10): 
	line = lines[i][:-1].split(",")
	for i in range(0,10): 
		outfile.write(str(line[i])+"\t")
	outfile.write("\n")

outfile.close() #done!
