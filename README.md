# Phage-Host_Codon_Compatibility
Repo for Phage-Host Codon Compatibility Pipeline


**Instructions**

1. Download all files into a directory

2. Navigate to the directory containing the the downloaded repo files

3. Execute the command: python3 wrapper.py -s XXX -q XXX
        
        Arguments:
        
        -s AA012345.1    (Accession of bacterial genome or file path to bacterial CDS fasta file)
        -q A01234.1      (Accession of phage genome or file path to phage CDS fasta file)


**Files used for testing:**
    
    Bacteria -- E. coli K12 strain C3026 (GCF_001559675.1); Accession: CP014272.1
    Positive Control Phage -- Escherichia phage lambda; Accession: J02459.1
    Negative Control Phage -- Lactobacillus phage ATCC 8014-B1; Accession: NC_019916.1


**Testing**

1) Try to run wrapper.py with the files provided 

        python3 wrapper.py -s Bacteria.txt -q PCPhage.txt
        python3 wrapper.py -s Bacteria.txt -q NCPhage.txt

2) Try to run wrapper.py using record accessions

        python3 wrapper.py -s CP014272.1 -q J02459.1

3) Go on NCBI, find a random bacterial + phage sequence and run them (using accessions would be easiest)

        python3 wrapper.py -s XXXX -q XXXX
