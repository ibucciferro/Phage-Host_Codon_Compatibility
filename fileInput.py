import os


def parse(file_path):
    # Opens fasta CDS files, parses them and returns a dictionary containing the gene name and sequence
    f = open(file_path, 'r')
    fs = f.read().strip().split('>')
    final_genes = {}
    for i in fs[1:]:
        si = i.split(' ')
        for comp in si:
            if 'gbkey' in comp:
                final_genes[si[0][4:]] = comp[12:].replace('\n', '')
    return final_genes

# Test
i = 'Phage.txt'
j = 'Bacteria.txt'
parse(i)
parse(j)

#create a definition to take the parsed file from the Phage, split the codons into 



