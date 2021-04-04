from Bio import SeqIO
import pandas
import os
import subprocess

cwd = os.getcwd()
## Uses DIAMOND on the database called testDB that contains a database assembled from the identical protein groups NCBI database of the
## 40 highly expressed genes in bacteria. Outputs a file called matches with the sequence title, bitscore, and query sequence id. K is specified to avoid redundancy and get
## the top hit for each query. Diamond has a binary in the parent working directory in this application. Diamond could be installed to avoid this change in directory.

def _get_hegs(file):
    os.chdir(cwd)
    os.system("./diamond makedb --in protein_database.fasta -d dmndDB")
    os.system('./diamond -d dmndDB.dmnd -q %s -out /matches -f 6 stitle bitscore qseqid' % file)


_get_hegs('CP014272.fasta')


def _get_hegs_to_forty(bact_genes):
    df = pandas.read_table("matches", names=["Subject", "Bit", "SeqID"], skipinitialspace=True)
    df = df.replace('\[.*\]', '', regex=True)
    df["Subject"] = df["Subject"].str.strip()
    df["Subject"] = df["Subject"].apply(lambda x: ' '.join(x.split(' ')[1:]))
    df["Subject"] = df["Subject"].str.lower()
    df = df.replace("elongation factor ef-2", "elongation factor g")
    df2 = df.sort_values(["Subject", "Bit"], ascending=[True, False])
    df2 = df2.loc[df2.groupby('Subject')["Bit"].idxmax()].reset_index(drop=True)
    items = df2.SeqID.unique()
    newSeqs = []
    for seq_record in SeqIO.parse(self.file, "fasta"):
        if (seq_record.id in items):
            newSeqs.append(seq_record)
    if len(newSeqs) < 38:
        print("WARNING there are less than 38 sequences.")
    with open("HEGS.fasta", "w") as handle:
        SeqIO.write(newSeqs, handle, "fasta")
    return len(newSeqs)



