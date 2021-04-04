from Bio import SeqIO
import pandas
import os
cwd = os.getcwd()

'''
1)  Uses DIAMOND to run a local BLASTx on the database created from the protein_db.fasta file
        (for more info about the protein_db.fasta file see the git repo wiki)
    
    Outputs a file named matches : sequence title, bitscore, query sequence id
        K is specified to avoid redundancy and get the top hit for each query
        Diamond runs via a binary in this application
2)  Trims the number of HEG matches from the protein_db.fasta to 40 (if possible)

    Outputs a file named HEGs.fasta : : sequenceID, sequence
'''


def _get_hegs(file):
    os.chdir(cwd)
    os.system("./diamond makedb --in protein_database.fasta -d dmndDB")
    os.system('./diamond blastx -d dmndDB.dmnd -q %s -o matches -f 6 stitle bitscore qseqid -p 1 -k 1' % file)

    def _get_hegs_to_forty():
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
        for seq_record in SeqIO.parse(file, "fasta"):
            if seq_record.id in items:
                newSeqs.append(seq_record)
        if len(newSeqs) < 38:
            print("WARNING there are fewer than 38 sequences.")
        with open("HEGS.fasta", "w") as handle:
            SeqIO.write(newSeqs, handle, "fasta")
        return len(newSeqs)

    _get_hegs_to_forty()


_get_hegs('Bacteria.txt')



