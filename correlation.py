from scipy import stats

# Reads the bacterial HEG codon frequencies into an array
def readBactCodonFeqs(bacteriaf):
    b = open(bacteriaf, 'r')

    b_s = b.read().strip().split('\n')

    name = b_s.pop(0)
    bcodon_nums = []
    for d in sorted(b_s):
        c_s = d.split('\t')

        for h in c_s:
            n_s = h.split(':')
            bcodon_nums.append(float(n_s[1]))
    return name, bcodon_nums


# Reads the phage gene codon frequencies into an array
def readPhageCodonFeqs(phagef):
    p = open(phagef, 'r')
    p_s = p.read().strip().split('\n')

    genenames = []
    pcodon_nums = []
    codons = []

    for r in p_s:
        if 'NP' in r:
            genenames.append(r)
        else:
            codons.append(r)

    count = 0
    for k in codons:
        cs = k.split('\t')
        try:
            cs.remove('')
        except:
            pass
        curr_nums = []

        if len(cs) != 64:               # Ask Dr. Wheeler
            print('The phage gene %s contains an abnormal codon! The codon %s will not be included in the correlation calculation.'
                  % (genenames[count], cs[-1][:3]))
            cs.pop()

        for s in cs:
            s_s = s.split(':')
            curr_nums.append(float(s_s[1]))
        count += 1
        pcodon_nums.append(curr_nums)

    return genenames, pcodon_nums


def calculation(bcodons, bname, pcodons, pgenes):
    out = open('phage_host_codon_correlation.txt', 'w')
    out.write('%s\tphage accession\n' % bname)


    def genes():
        for gs in range(len(pcodons)-1):
            gcorr = stats.pearsonr(bcodons, pcodons[gs])[0]
            gn = pgenes[gs]
            out.write('\n%s\t%f' % (gn, gcorr))

    def allgenes():
        curr_pos = 0
        pglobcodons = []
        while curr_pos < len(pcodons[1]):
            value = 0
            for ag in range(len(pcodons)-1):
                value += pcodons[ag][curr_pos]
            pglobcodons.append(value)
            curr_pos += 1

        globcorr = stats.pearsonr(bcodons, pglobcodons)[0]
        out.write('\n%f\n' % globcorr)

    allgenes()
    genes()

bact, bnums = readBactCodonFeqs('bacteriaHEGCodons.txt')
phage_genes, pnums = readPhageCodonFeqs('phageCodons.txt')
calculation(bnums, bact, pnums, phage_genes)
