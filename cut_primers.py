from fastq import *

import sys


PRIMERS_FOR  = {"c1a1_for": "CCCTGACCSAGACCTG", 
                "dqba1_for": "AGKCTTTGCGGATCCC", 
                "drba1_for": "CTGAGCTCCCSACTGG", 
                "c1a2_for": "CGACGGCAARGATTAC", 
                "c2a2_for": "GGAACAGCCAGAAGGA"
               }
PRIMERS_REV  = {"c1a1_rev": "CCMTCCAGGTAGRCTCT", 
                "c2a1_rev": "YCAGHAGGTTGTGGTG", 
                "c1a2_rev": "YGGTGGYCTGGGAAGA", 
                "c2a2_rev": "CCACKTGGCAGGTGTA"
               }
PRIMERS_PAIR = {("c1a1_for", "c1a1_rev"): "class1_amp1", 
                ("c1a2_for", "c1a2_rev"): "class1_amp2",
                ("dqba1_for", "c2a1_rev"): "DQB_amp1",
                ("drba1_for", "c2a1_rev"): "DRB_amp1",
                ("c2a2_for", "c2a2_rev"): "class2_amp2"
               }

MAX_ERRORS = 4


def cut_seq(fqrec, n):
    return FASTQData(name = fqrec.name, attr = fqrec.attr, seq = fqrec.seq[n:], qual = fqrec.qual[n:])


def cut_tuple(fqtuple, n1, n2):
    return cut_seq(fqtuple[0], n1), cut_seq(fqtuple[1], n2)


def hamm(alpha, beta, max_err = MAX_ERRORS):
    if len(alpha) != len(beta): 
        print("!!!!")
    err = 0
    for a, b in zip(alpha, beta):
        err += (a != b) and (a != 'N') and (b != 'N')
        if err > max_err:
            return False
    return True


def check_primers(fq_tuple):
    left_pr = "undefined"
    for key, val in PRIMERS_FOR.items():
        if hamm(fq_tuple[0].seq[:len(val)], val):
            left_pr = key
            break
    
    right_pr = "undefined"
    for key, val in PRIMERS_REV.items():
        if hamm(fq_tuple[1].seq[:len(val)], val):
            right_pr = key
            break
    
    if left_pr != "undefined" and right_pr !=  "undefined":
        return (left_pr, right_pr), cut_tuple(fq_tuple, len(PRIMERS_FOR[left_pr]), len(PRIMERS_REV[right_pr]))
    else:
        return (left_pr, right_pr), 0
    

if __name__ == "__main__":
    f1 = sys.argv[1]
    f2 = f1.replace("R1.fastq", "R2.fastq")
    parser = PairendFASTQParser(f1, f2)
    
    stats = {}
    writers = {pair: PairendFASTQWriter(f1.replace("R1.fastq", "__" + PRIMERS_PAIR[pair] + "__.R1.fastq"),
                                        f2.replace("R2.fastq", "__" + PRIMERS_PAIR[pair] + "__.R2.fastq")
                                       ) for pair in PRIMERS_PAIR}
    
    for fq_tuple in parser:
        pr_pair, fq_tuple_cut = check_primers(fq_tuple)
        stats[pr_pair] = stats.get(pr_pair, 0) + 1
        if pr_pair in PRIMERS_PAIR:
            writers[pr_pair].write(fq_tuple_cut)
    
    for pair, cnt in sorted(stats.items(), reverse = True, key = lambda x: x[1]):
        print(pair[0], " ", pair[1], "\t", cnt, sep = "", end = "\t")
        if pair in PRIMERS_PAIR:
            print(PRIMERS_PAIR[pair])
        else:
            print("~")
            
    for key, wr in writers.items():
        wr.close()