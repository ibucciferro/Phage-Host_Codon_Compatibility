from Bio import SeqIO
import pandas
import os
import argparse
from Bio import Entrez, SeqIO
Entrez.email = "ecrum@luc.edu"


# supply the arguments providing the bacterial (-s flag) and phage (-q flag) genomes for the pipeline
phcc = argparse.ArgumentParser(description='Use GenBank Accession (-a) or Local File (-l)')
phcc.add_argument('-s', '--subject', nargs='?', default='CP014272.1',
                 help='"-s" flag specifies the GenBank accession or local CDS file of the bacteria to be used.')

phcc.add_argument('-q', '--query', nargs='?', default='J02459.1',
                 help='"-q" flag specifies the GenBank accession or local CDS file of the phage to be used.')

args = phcc.parse_args()
bacteria = vars(args)['subject']
print("Using %s for the bacterial input" % bacteria)
phage = vars(args)['query']
print("Using %s for the phage input" % phage)

# Keep track of current directory
cwd = os.getcwd()

#-----------------------------------------------------------

'''
0) (optional) Fetches CDS fasta sequence from a provided accession
1) Parses inputted CDS fasta file and returns {geneID: sequence}
2) Determines codon frequency for each phage gene and writes output to phageGeneCodons.txt file
'''

# Creates a fasta CDS file from a provided genbank accession via the Entrez database
def getFasta(accession):
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()

        inSeq = cwd + '/' + accession + '.fasta'
        getSeq = open(inSeq, 'w')
        c = 1
        for rec in record.features:  # iterate through list of features
            if rec.type == "CDS":
                try:
                    getSeq.write(">" + str(rec.qualifiers['protein_id'])[2:-2] + "\n")  # protein IDs
                    getSeq.write(str(rec.extract(record.seq)) + "\n")  # gene sequence
                except:
                    getSeq.write(">" + str(record.id) + '.' + str(c) + "\n")  # protein IDs
                    getSeq.write(str(rec.extract(record.seq)) + "\n")  # gene sequence
                    c += 1
        getSeq.close()
        return inSeq
    except:
        return "Not a valid NCBI Accession"


# Opens fasta CDS files, parses them, returns a dictionary {geneID: sequence}
def parseAccessionFasta(file_path):
    f = open(file_path, 'r')
    fs = f.read().strip().split('>')
    final_genes = {}
    for i in fs[1:]:
        si = i.rstrip().split('\n')
        final_genes[si[0]] = si[1]
    return final_genes


def parseLocalFasta(file_path):
    f = open(file_path, 'r')
    fs = f.read().strip().split('>')
    final_genes = {}
    for i in fs[1:]:
        si = i.rstrip().split(' ')
        for comp in si:
            if 'gbkey' in comp:
                final_genes[si[0]] = comp[12:].replace('\n', '')
    return final_genes


def parseHEGFasta(file_path):
    f = open(file_path, 'r')
    fs = f.read().strip().split('>')
    final_genes = {}
    for i in fs[1:]:
        si = i.rstrip().split('.')
        for comp in si:
            final_genes[si[0] + '.1'] = comp[2:].replace('\n', '')
    return final_genes


# Creates a tab delineated file documenting each phage gene's codon frequencies.
def phageCodons(phageGeneDict):
    codonsDict = {
        'TTT': 0, 'TTC': 0, 'TTA': 0, 'TTG': 0, 'CTT': 0,
        'CTC': 0, 'CTA': 0, 'CTG': 0, 'ATT': 0, 'ATC': 0,
        'ATA': 0, 'ATG': 0, 'GTT': 0, 'GTC': 0, 'GTA': 0,
        'GTG': 0, 'TAT': 0, 'TAC': 0, 'TAA': 0, 'TAG': 0,
        'CAT': 0, 'CAC': 0, 'CAA': 0, 'CAG': 0, 'AAT': 0,
        'AAC': 0, 'AAA': 0, 'AAG': 0, 'GAT': 0, 'GAC': 0,
        'GAA': 0, 'GAG': 0, 'TCT': 0, 'TCC': 0, 'TCA': 0,
        'TCG': 0, 'CCT': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0,
        'ACT': 0, 'ACC': 0, 'ACA': 0, 'ACG': 0, 'GCT': 0,
        'GCC': 0, 'GCA': 0, 'GCG': 0, 'TGT': 0, 'TGC': 0,
        'TGA': 0, 'TGG': 0, 'CGT': 0, 'CGC': 0, 'CGA': 0,
        'CGG': 0, 'AGT': 0, 'AGC': 0, 'AGA': 0, 'AGG': 0,
        'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0}

    of = open("phageGeneCodons.txt", 'w')

    for gene in phageGeneDict:
        seq = phageGeneDict[gene]
        sseq = [seq[i:i + 3] for i in range(0, len(seq), 3)]

        of.write(gene + '\n')

        freshdict = codonsDict.copy()

        for c in sseq:
            curr = freshdict[c]
            curr += 1
            freshdict[c] = curr

        [of.write('%s:%d\t' % (a, freshdict[a])) for a in freshdict]
        of.write('\n')


# Driver Code (Phage)
print("Fetching the %s's CDS..." % phage)
yap = getFasta(phage)
if 'Not a' in yap and '/' in phage:
    phagedict = parseLocalFasta(phage)
elif 'Not a' in yap:
    try:
        phagedict = parseLocalFasta(phage)
    except:
        print("Not a valid phage accession or input file")
else:
    phagedict = parseAccessionFasta(yap)

try:
    phageCodons(phagedict)
    print("Success, the phage codon frequencies are within phageGeneCodons.txt")
except:
    print("problem :(")


#--------------------------------------------------------------------------------------------------


'''
1)  Uses DIAMOND to run a local BLASTx on the database created from the protein_db.fasta file
        (for more info about the protein_db.fasta file see the git repo wiki)
    
    Outputs a file named matches : sequence title, bitscore, query sequence id
        K is specified to avoid redundancy and get the top hit for each query
        Diamond runs via a binary in this application

2)  Trims the number of HEG matches from the protein_db.fasta to 40 (if possible)
    Outputs a file named HEGs.fasta : : sequenceID, sequence
'''

# Run the Diamond binary using wget functions below
def _get_hegs(file):
    os.chdir(cwd)
    currdir = os.listdir()
    if 'diamond' not in currdir:
        os.system('wget http://github.com/bbuchfink/diamond/releases/download/v0.9.24/diamond-linux64.tar.gz')
        os.system('tar xzf diamond-linux64.tar.gz')
    os.system("nohup ./diamond makedb --in protein_database.fasta -d dmndDB")
    os.system('nohup ./diamond blastx -d dmndDB.dmnd -q %s -o matches -f 6 stitle bitscore qtitle -p 1 -k 1' % file)

    def _get_hegs_to_forty():
        df = pandas.read_table("matches", names=["Subject", "Bit", "Query"], skipinitialspace=True)
        df = df.replace('\[.*\]', '', regex=True)
        df["Subject"] = df["Subject"].str.strip()
        df["Subject"] = df["Subject"].apply(lambda x: ' '.join(x.split(' ')[1:]))
        df["Subject"] = df["Subject"].str.lower()
        df = df.replace("elongation factor ef-2", "elongation factor g")
        df2 = df.sort_values(["Subject", "Bit"], ascending=[True, False])
        df2 = df2.loc[df2.groupby('Subject')["Bit"].idxmax()].reset_index(drop=True)
        items = df2.Query.unique()
        newSeqs = []
        for seq_record in SeqIO.parse(file, "fasta"):
            if seq_record.id in items:
                newSeqs.append(seq_record)
        if len(newSeqs) < 38:
            print("WARNING there are fewer than 38 HEGs.")
        with open("HEGS.fasta", "w") as handle:
            SeqIO.write(newSeqs, handle, "fasta")
        return "HEGS.fasta"
    _get_hegs_to_forty()

# Creates a tab delineated file documenting all HEGs codon frequencies.
def bacteriaCodons(bacteriaGeneDict):
    hegDict = {
        'TTT': 0, 'TTC': 0, 'TTA': 0, 'TTG': 0, 'CTT': 0,
        'CTC': 0, 'CTA': 0, 'CTG': 0, 'ATT': 0, 'ATC': 0,
        'ATA': 0, 'ATG': 0, 'GTT': 0, 'GTC': 0, 'GTA': 0,
        'GTG': 0, 'TAT': 0, 'TAC': 0, 'TAA': 0, 'TAG': 0,
        'CAT': 0, 'CAC': 0, 'CAA': 0, 'CAG': 0, 'AAT': 0,
        'AAC': 0, 'AAA': 0, 'AAG': 0, 'GAT': 0, 'GAC': 0,
        'GAA': 0, 'GAG': 0, 'TCT': 0, 'TCC': 0, 'TCA': 0,
        'TCG': 0, 'CCT': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0,
        'ACT': 0, 'ACC': 0, 'ACA': 0, 'ACG': 0, 'GCT': 0,
        'GCC': 0, 'GCA': 0, 'GCG': 0, 'TGT': 0, 'TGC': 0,
        'TGA': 0, 'TGG': 0, 'CGT': 0, 'CGC': 0, 'CGA': 0,
        'CGG': 0, 'AGT': 0, 'AGC': 0, 'AGA': 0, 'AGG': 0,
        'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0}

    of = open("bacteriaHEGCodons.txt", 'w')

    f = ''
    for heg in bacteriaGeneDict:
        f = heg
        bseq = bacteriaGeneDict[heg]
        bsseq = [bseq[b:b + 3] for b in range(0, len(bseq), 3)]
        for q in bsseq:
            bcurr = hegDict[q]
            bcurr += 1
            hegDict[q] = bcurr

    of.write(f[:10] + '\n')
    [of.write('%s:%d\t' % (a, hegDict[a])) for a in hegDict]


#Driver Code (Bacteria)
print("Fetching the bacteria's CDS...")
yab = getFasta(bacteria)
if 'Not a' in yab and '/' in bacteria:
    _get_hegs(bacteria)
elif 'Not a' in yab:
    try:
        _get_hegs(bacteria)
    except:
        print("Not a valid accession or input file format")
else:
    _get_hegs(yab)

try:
    bacteriaCodons(parseHEGFasta('HEGS.fasta'))
    print("Success, the bacterial codon frequencies are within bacteriaHEGCodons.txt")
except:
    print("problem :(")

