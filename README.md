# Phage-Host_Codon_Compatibility
Repo for Phage-Host Codon Compatibility Pipeline


**Instructions**

1. Download all files into a directory

2. Navigate to the directory containing the the downloaded repo files

3. Execute the command: python3 pipeline.py -s XXX -q XXX
        
        Arguments:
        
        -s AA012345.1    (Accession of bacterial genome or file path to bacterial CDS fasta file)
        -q J02459.1      (Accession of phage genome or file path to phage CDS fasta file)

**Files used for testing:**
    
    Bacteria -- E. coli K12 strain C3026 (GCF_001559675.1): CP014272.1
    Positive Control Phage -- Escherichia phage lambda: J02459.1
    Negative Control Phage -- Lactobacillus phage ATCC 8014-B1: NC_019916.1
