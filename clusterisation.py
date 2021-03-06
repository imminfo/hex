import sys
import os
import random
import math
import statistics
from collections import namedtuple

import seaborn as sns
import numpy as np

from fasta import *
from fastq import *
from log_progress import *


Read = namedtuple('Read', 'index seq qual')


def replace_low_quality(data, q_threshold):
    res = ''
    res_q = ''
    for nuc, q in zip(data.seq, data.qual):
        if ord(q) - ord('!') <= q_threshold:
            res += 'N'
            res_q += '!'
        else:
            res += nuc
            res_q += q

    return FASTQData(seq = res, qual = res_q, name = data.name, attr = data.attr)


def median(data):
    return statistics.median([ord(q) - ord('!') for q in data.qual])


def hamm(alpha, beta, max_err = 0):
    # if len(alpha) != len(beta): 
    #     print("!!!!")
    err = 0
    for a, b in zip(alpha, beta):
        err += (a != b) and (a != 'N') and (b != 'N')
        if err > max_err:
            return False
    return True


def make_best_pos_consensus(orig_seq, orig_num, seq_list):
    res = ''
    for i in range(len(orig_seq)):
        d = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}
        d[orig_seq[i]] += orig_num
        for seq, cnt in seq_list:
            d[seq[i]] += cnt
        d.pop('N')
        res += max(d.items(), key = lambda x: x[1])[0]
    return res


def classify(fqdata, low_threshold, high_threshold):
    qmed = median(fqdata)
    if qmed > low_threshold and qmed < high_threshold:
        return 1
    elif qmed >= high_threshold:
        return 2
    else:
        return 0

    
def cls_stats(cls_seq, which_class):
    print("\t", which_class, " seq:\t", len(cls_seq[which_class]), sep = "", end = "")
    print("(", round(len(cls_seq[which_class]) / sum(map(len, cls_seq)), 4), ")", sep = "")
    print("\t", which_class, " rds:\t", sum(cls_seq[which_class].values()), sep = "", end = "")
    print("(", round(len(cls_seq[which_class]) / sum(map(lambda x: sum(x.values()), cls_seq)), 4), ")", sep = "")
        

def classify_fastq_file(filepath, replace_threshold, low_threshold, high_threshold):    
    print("Fastq records:\t", int(num_lines(filepath) / 4), sep = "")
    
    q_distr = []
    
    # median - sequences
    seqmed = {}
    cls_seq = [{}, {}, {}]
    for fqdata in FASTQParser(filepath):
        fqdata = replace_low_quality(fqdata, replace_threshold)
        med = median(fqdata)
        q_distr.append(med)
        if med not in seqmed: seqmed[med] = []
        seqmed[med].append(fqdata.seq)
        
    lo, hi = np.trunc(np.percentile(q_distr, [low_threshold, high_threshold]))
    print("Lo / hi percentile median values:", (lo, hi), sep = "\t")
        
    kmers = {}
    for med in seqmed:
        cls = 2
        if med < lo:
            cls = 0
        elif med >= lo and med < hi:
            cls = 1
        for seq in seqmed[med]:
            cls_seq[cls][seq] = cls_seq[cls].get(seq, 0) + 1
            if cls == 0 or cls == 1:
                kmers[seq[10:20]] = kmers.get(seq[10:20], 0) + 1
    print(sorted(kmers.items(), reverse = True, key = lambda x: x[1]))
        
    print("Classes:")
    cls_stats(cls_seq, 0)
    cls_stats(cls_seq, 1)
    cls_stats(cls_seq, 2)
    print()
    
#     sns.distplot(list(kmers.values()), bins = max(list(kmers.values())) - min(list(kmers.values())) + 1, axlabel = "Kmer distribution");
    sns.distplot(q_distr, bins = max(q_distr) - min(q_distr) + 1, axlabel = "Median quality distribution");
    
    return cls_seq


def X_clust(minors, majors, x_clust_hamm):
    """
    N-clust: merge bad sequences with all other sequences.
    H-clust: merge medium quality sequences with the high quality ones.
    """
    
    def _get_keys(seq_dict):
        keys = list(seq_dict.keys())
        key_ind = dict([(keys[i], i) for i in range(len(keys))])
        return keys, key_ind
        
    cands = {}
    merged_seq = {x: set([]) for x in majors}
    minor_keys, minor_key_inds = _get_keys(minors)
    major_keys, major_keys_inds = _get_keys(majors)
    
    
    print("Computing distances...")
    for i in range(len(minors)):
        for j in range(len(majors)):
            if hamm(minor_keys[i], major_keys[j], x_clust_hamm):
                if i not in cands: cands[i] = []
                cands[i].append(j)

                
    print("Merging sequences...")
    for minor_key, targets in cands.items():
        for target in targets:
            merged_seq[major_keys[target]].add((minor_key, minors[minor_keys[minor_key]] / len(targets)))

    cand_keys = set([minor_keys[k] for k in cands])
    new_minors = {x: minors[x] for x in minor_keys if x not in cand_keys}
    print("# candidates:", len(cands))
    print("# distants", len(new_minors))
        
    print("Making consensuses...")
#     sum_pre = sum(majors.values())
    n_merged = 0
    new_seq_dict = {}
    for seq, seq_ls in merged_seq.items():
        new_seq = make_best_pos_consensus(seq, majors[seq], [(minor_keys[x[0]], x[1]) for x in seq_ls])
        if new_seq not in new_seq_dict:
            new_seq_dict[new_seq] = 0
        else:
            n_merged += 1
        new_seq_dict[new_seq] += majors[seq] + sum([x[1] for x in merged_seq[seq]])

    print("# merged:", n_merged)
#     sum_post = sum(new_seq_dict.values())
#     if sum_pre != sum_post:
#         print("Sums are not equal!", sum_pre, "vs.", sum_post)

    for seq in new_seq_dict:
        new_seq_dict[seq] = round(new_seq_dict[seq], 3)

    return new_minors, new_seq_dict


def write_and_blast(seq_dict, f1, out_seq, out_blast, max_sequences):
    with open(out_seq, 'w') as file:
        i = 0
        for key, val in reversed(sorted(seq_dict.items(), key = lambda x: x[1])):
            # print(val, " (", round(100 * val / sum(final_seq.values()), 4), "%)", sep = '')
            print(val, key, sep = '\t', file = file)
            i += 1
            if i == max_sequences: break
                
    ls = []
    i = 0
    for key, val in reversed(sorted(seq_dict.items(), key = lambda x: x[1])):
        ls.append(faseq(name = "sequence" + str(i) + "_" + str(val) + "_(" + str(round(100 * val / sum(seq_dict.values()), 4)) + ")", seq = key, comm = ''))
        i += 1
        if i == max_sequences: break
    write_fasta(f1 + ".seq.fasta.txt", ls)
    os.system("blastn -query " + f1 + ".seq.fasta.txt" + " -db hlabase/hlabase.fasta -outfmt 6 -num_alignments 4 > " + out_blast)


def clusterise_sequences(f1, replace_threshold, low_threshold, high_threshold, n_clust_hamm, h_clust_hamm, max_sequences, out_seq_prefix = "tmp.topseq1", out_blast_prefix = "tmp.blast1"):
    prefix = f1[:f1.find(".fastq")]
    
    
    print("*** Searching for unique sequences ***")
    cls_seq = classify_fastq_file(f1, replace_threshold, low_threshold, high_threshold)

    
    print("*** N-clusterisation ***")    
    cls_seq[0], cls_seq[1] = X_clust(cls_seq[0], cls_seq[1], n_clust_hamm)
    print("Intermediate statistics by class:")
    cls_stats(cls_seq, 0)
    cls_stats(cls_seq, 1)
    print()
    
    
    print("*** H-clusterisation ***")
    cls_seq[1], cls_seq[2] = X_clust(cls_seq[1], cls_seq[2], h_clust_hamm)    
    print("Final statistics by class:")
    cls_stats(cls_seq, 1)
    cls_stats(cls_seq, 2)
    print()
    
    cls_seq[2].update(cls_seq[1])

    print("Move minors to the major class.")
    print("Final number of sequences:", len(cls_seq[2]), sep = "\t")
          
    print("*** Writing results ***")
#     write_and_blast(cls_seq[1], f1 + ".minor", out_seq_prefix + ".minor.txt", out_blast_prefix + ".minor.txt", max_sequences)
    write_and_blast(cls_seq[2], f1 + ".major", out_seq_prefix + ".major.txt", out_blast_prefix + ".major.txt", max_sequences)
    print("\n*** DONE ***")
    

if __name__ == '__main__':
    aggregate_sequences(sys.argv[1], 50, 7, 5, "tmp.topseq1.txt", "tmp.blast1.txt")
    aggregate_sequences(sys.argv[2], 50, 7, 5, "tmp.topseq2.txt", "tmp.blast2.txt")