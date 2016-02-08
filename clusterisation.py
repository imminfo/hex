import sys
import os
import random
import math
import statistics
from collections import namedtuple

from fasta import *
from fastq import *
from log_progress import *


Read = namedtuple('Read', 'index seq qual')


def replace_low_quality(data, q_threshold):
    res = ''
    for nuc, q in zip(data.seq, data.qual):
        if ord(q) - ord('!') <= q_threshold:
            res += 'N'
        else:
            res += nuc

    return FASTQData(seq = res, qual = data.qual, name = data.name, attr = data.attr)



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
    if qmed >= low_threshold and qmed < high_threshold:
        return 1
    elif qmed >= high_threshold:
        return 2
    return 0


def classify_fastq_file(filepath, replace_threshold, low_threshold, high_threshold):
    def _stats(cls_seq, which_class):
        print("\t", which_class, " seq:\t", len(cls_seq[which_class]))
        print("(", round(len(cls_seq[which_class]) / sum(map(len, cls_seq), 4), ")", sep = "", end = ""))
        print("\t", which_class, " rds:\t", len(sum(cls_seq[which_class].values())))
        print("(", round(len(cls_seq[which_class]) / sum(map(lambda x: sum(x.values()), cls_seq), 4), ")", sep = "", end = ""))
    
    cls_seq = [{}, {}, {}]
    for fqdata in FASTQParser(filepath):
        fqdata = replace_low_quality(fqdata, replace_threshold)
        cls = classify(fqdata, low_threshold, high_threshold)
        cls_seq[cls][fqdata.seq] = cls_seq[cls].get(fqdata.seq, 0) + 1
    print("Classes:")
    _stats(cls_seq, 0)
    _stats(cls_seq, 1)
    _stats(cls_seq, 2)
    print()
    return cls_seq


def N_clust(cls_seq, n_clust_hamm):
    pass


def H_clust(cls_seq, h_clust_hamm):
    pass
    

def merge_with_clusters(seq_dict, seq_threshold, qmed_threshold):
    print("sum pre filter", sum(seq_dict.values()))
    seq_dict = dict(filter(lambda x: x[1] >= seq_threshold, seq_dict.items()))
    print("sum post filter", sum(seq_dict.values()))

    keys = list(seq_dict.keys())
    key_ind = dict([(keys[i], i) for i in range(len(keys))])
    cands = {}
    merged_seq = {x: set([]) for x in keys}
    for i in range(len(seq_dict) - 1):
        for j in range(i + 1, len(seq_dict)):
            if hamm(keys[i], keys[j]):
                key = i
                append_to = j
                if seq_dict[keys[i]] > seq_dict[keys[j]]:
                    key = j
                    append_to = i

                if key not in cands: cands[key] = []
                cands[key].append(append_to)

    while cands:
        for seq, cnt in sorted(seq_dict.items(), key = lambda x: x[1]):
            if key_ind[seq] in cands:
                for append_to in cands[key_ind[seq]]:
                    # seq_dict[keys[append_to]] += math.ceil(seq_dict[seq] / len(cands[key_ind[seq]]))
                    seq_dict[keys[append_to]] += seq_dict[seq] / len(cands[key_ind[seq]])
                    if merged_seq[seq]:
                        merged_seq[keys[append_to]] = merged_seq[keys[append_to]].union(merged_seq[seq])
                    # merged_seq[keys[append_to]].add((key_ind[seq], math.ceil(seq_dict[seq] / len(cands[key_ind[seq]]))))
                    merged_seq[keys[append_to]].add((key_ind[seq], seq_dict[seq] / len(cands[key_ind[seq]])))
                seq_dict[seq] = 0

        to_pop = set([x for x in cands])
        for source in cands:
            for to_append in cands[source]:
                if seq_dict[keys[source]] and source in to_pop:
                    to_pop.remove(source)

        for p in to_pop:
            cands.pop(p)
            merged_seq.pop(keys[p])
            seq_dict.pop(keys[p])


    print("Making consensuses...")

    print("sum pre merging", sum(seq_dict.values()))
    n_merged = 0
    new_seq_dict = {}
    for seq, seq_ls in merged_seq.items():
        new_seq = make_best_pos_consensus(seq, seq_dict[seq], [(keys[x[0]], x[1]) for x in seq_ls])
        if new_seq not in new_seq_dict:
            new_seq_dict[new_seq] = 0
        else:
            n_merged += 1
        new_seq_dict[new_seq] += seq_dict[seq]

    print("# merged:", n_merged)
    print("sum post merging", sum(new_seq_dict.values()))
    return new_seq_dict


def clusterise_sequences(f1, replace_threshold, low_threshold, high_threshold, out_seq = "tmp.topseq1.txt", out_blast = "tmp.blast1.txt"):
    prefix = f1[:f1.find(".fastq")]
    print("Searching for unique sequences..."); d1 = {}

    classify_fastq_file(f1, replace_threshold, low_threshold, high_threshold)


    print("Clustering error sequences...")
#     d1 = merge_with_clusters(d1, seq_threshold, qmed_threshold)


#     print("Writing results...")

#     with open(out_seq, 'w') as file:
#         i = 0
#         for key, val in reversed(sorted(d1.items(), key = lambda x: x[1])):
#             # print(val, " (", round(100 * val / sum(d1.values()), 4), "%)", sep = '')
#             print(val, key, sep = '\t', file = file)
#             i += 1
#             if i == max_sequences: break
#             # if val < 2: break

#     ls = []
#     i = 0
#     for key, val in reversed(sorted(d1.items(), key = lambda x: x[1])):
#         ls.append(faseq(name = "sequence" + str(i) + "_" + str(val) + "_(" + str(round(100 * val / sum(d1.values()), 4)) + ")", seq = key, comm = ''))
#         i += 1
#         if i == max_sequences: break
#     write_fasta(f1 + ".seq.fasta.txt", ls)
#     os.system("blastn -query " + f1 + ".seq.fasta.txt" + " -db hlabase/hlabase.fasta -outfmt 6 -num_alignments 4 > " + out_blast)


if __name__ == '__main__':
    aggregate_sequences(sys.argv[1], 50, 7, 5, "tmp.topseq1.txt", "tmp.blast1.txt")
    aggregate_sequences(sys.argv[2], 50, 7, 5, "tmp.topseq2.txt", "tmp.blast2.txt")