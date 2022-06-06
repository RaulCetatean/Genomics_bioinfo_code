from urllib.request import urlopen
import random
import math
import sys

file = "chr22.fa"

def analyze(file):

    '''analyze(file) analyzes the FASTA file, counts the nucleotides in the file and creates two strings; one representing the
    cpg islands and the other strings simulates a sequence containing non cpg regions. '''

    fasta = open(file)
    content = fasta.read().upper()
    lines = content.splitlines()
    f = open("chr22q.txt", "a+")
    
    for e in lines: 
        if ">" in e or e[0] == "N" and e[-1] == "N":
            lines.remove(e)
        else:
            f.write(e)

    nucleotides = f.read()       
    all_chr = open("all_chr.txt", "a+")
    for l in lines:
        if l[0] != '>':
            all_chr.write(l)

    r_all = all_chr.read()
        
    request_url = "http://www.haowulab.org/software/makeCGI/model-based-cpg-islands-hg19.txt"
    with urlopen(request_url) as response:
        content = response.read().decode()
        content = content.splitlines()
    
    CpG_islands = []
    for lines in content:
        lines = str(lines).split()
        if lines[0] == 'chr22':
            CpG_islands += [[int(lines[1])-1,int(lines[2]), int(lines[3])]]
                            #start          #end            #lenght
    
    tot_cpg = 0
    for e in CpG_islands:
        tot_cpg += e[2]

    avg = tot_cpg//len(CpG_islands)
    '''Since I already have here the lengths of every CpG island, I compute the the average length necessary
    for the scanning of the genome with the sliding window'''
    
    cpg = ''    #the string will contain all the CpG islands
    for p in CpG_islands:
        cpg += r_all[p[0]-1:p[1]]
   
    out_region = []
    for e in CpG_islands:
        point = random.randint(0, len(nucleotides)-e[2])
        out_region += [[point, e[2] + point]]

    non_cpg = ''    #the string will contain all the non CpG sequences

    for p in out_region:
        non_cpg += r_all[p[0]-1:p[1]]

    fasta.close()
    f.close()
    return cpg, non_cpg, avg

cpg,non_cpg, avg = analyze(file)

def probability(cpg, non_cpg):
        
    '''The probability(cpg, non_cpg) function computes the probabilities for the inside and the outside models based on
    the cpg and non_cpg sequences respectively.'''

    dn_i = counter(cpg)
    
    n_i = {"A":cpg.count("A") ,"C":cpg.count("C")  ,"G":cpg.count("G") ,"T":cpg.count("T")}
    
    prob_in = {'AA':dn_i["AA"]/n_i['A'], 'AC':dn_i["AC"]/n_i['A'], 'AG':dn_i["AG"]/n_i['A'], 'AT':dn_i["AT"]/n_i['A'],
              'CA':dn_i["CA"]/n_i['C'], 'CC':dn_i["CC"]/n_i['C'], 'CG':dn_i["CG"]/n_i['C'], 'CT':dn_i["CT"]/n_i['C'],
              'GA':dn_i["GA"]/n_i['G'], 'GC':dn_i["GC"]/n_i['G'], 'GG':dn_i["GG"]/n_i['G'], 'GT':dn_i["GT"]/n_i['G'],
              'TA':dn_i["TA"]/n_i['T'], 'TC':dn_i["TC"]/n_i['T'], 'TG':dn_i["TG"]/n_i['T'], 'TT':dn_i["TT"]/n_i['T']}

    '''AA should not be seen as A|A as normal but the inverse: e.g. prob_in["CA"] == p(A|C) because it easier for computation
    now i should compute the probability of being inside a cpg island and this should be the inside model
    remember that firt character has probability 0.25'''

    dn_o = counter(non_cpg)

    n_o = {"A":non_cpg.count("A") ,"C":non_cpg.count("C")  ,"G":non_cpg.count("G") ,"T":non_cpg.count("T")}
    
    prob_out = {'AA':dn_o["AA"]/n_o['A'], 'AC':dn_o["AC"]/n_o['A'], 'AG':dn_o["AG"]/n_o['A'], 'AT':dn_o["AT"]/n_o['A'],
               'CA':dn_o["CA"]/n_o['C'], 'CC':dn_o["CC"]/n_o['C'], 'CG':dn_o["CG"]/n_o['C'], 'CT':dn_o["CT"]/n_o['C'],
               'GA':dn_o["GA"]/n_o['G'], 'GC':dn_o["GC"]/n_o['G'], 'GG':dn_o["GG"]/n_o['G'], 'GT':dn_o["GT"]/n_o['G'],
               'TA':dn_o["TA"]/n_o['G'], 'TC':dn_o["TC"]/n_o['G'], 'TG':dn_o["TG"]/n_o['G'], 'TT':dn_o["TT"]/n_o['G']}
    
    #the probabilities should be seen as already stated for the prob_in

    return prob_in, prob_out
    

def prediction(sequence,cpg, non_cpg):

    '''prediction(sequence) takes a sequence and, using the probabilities of the inside and outside models from the
    prediction function, considers every dinucleotide of the sequence and computes the probabilities of the sequence
    to be inside and outside of a CpG island. Then it takes the logarithm of their ratio and predicts if the sequence
    belongs inside or outside a CpG island. '''

    prob_in, prob_out = probability(cpg, non_cpg)
    p_in = 0.25 #unconditioned probability of the first nuclueotide
    p_out = 0.25
    i = 0
    while i+2 < len(sequence): 
        p_in = p_in * prob_in[sequence[i:i+2]]
        p_out = p_out * prob_out[sequence[i:i+2]]
        i += 1
        
    s = math.log10(p_in/p_out)
    if s > 0 :
        return s, ' > 0, inside model is more probable than outside model'
    elif s < 0:
        return s, ' < 0, outside model is more probable than inside model'
    elif s == 0:
        return 's is ', s,'the probability of being inside is equal of the probability of being outside' 

def create(n):
    
    bases = ['A','C','G','T']
    sequence = ''
    for b in range(n):
        sequence += random.choice(bases)
    return sequence

def counter(s):

    '''counter(s) takes a string s as input and counts every dinucleotide, and considers also the overlapping cases.'''
    
    c = 0
    dinucleotides = {}
    while c+2 < len(s):
        if s[c:c+2] not in dinucleotides:
            dinucleotides[s[c:c+2]] = 1
            c += 1
        else:
            dinucleotides[s[c:c+2]] += 1
            c += 1
    return dinucleotides 

def scan(G):

    '''scan(G) takes a long sequence G and analyzes, using a window that is the average length of a CpG island,
    the score of being a sequence from the inside or the outside model. If the score is positive, the program 
    prints the score. '''

    prob_in, prob_out = probability(cpg,non_cpg)
    slide = avg
    i = 0
    while i + slide < len(G):
        score = prediction(G[i:slide +1])
        
        if score > 0:
            print(score)
        i += 1

def another_pred(seq):

    '''another_pred(seq) has the same function of prediction(sequence), but works better with longer
    iterations.'''

    prob_in, prob_out = probability(cpg,non_cpg)
    p_in = math.log10(0.25)
    p_out = math.log10(0.25)
    i = 0
    while i+2 < len(seq):
        p_in = p_in + math.log10(prob_in[seq[i:i+2]])
        p_out = p_out + math.log10(prob_out[seq[i:i+2]])
        i +=1
    print(p_in,p_out)
    s = p_in - p_out
    return s
    
if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: %s sequence, cpg sequence, non_cpg sequence" % sys.argv[0] )
	sys.exit()
    sequence = sys.argv[1]
    cpg = sys.argv[2]
    non_cpg = sys.argv[3]
    prediction(sequence, cpg, non_cpg)
    
#for the scanning of the genome
'''if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: genome G" % sys.argv[0])
        sys.exit()
    G = sys.argv[1]
    scan(G)'''    
    





    
