import sys
import argparse
from math import log
import itertools
import os


from fasta import *
from fastq import *


def read_clustalw(filename):
	sequences = []
	with open(filename) as inp:
		for line in inp:
			words = line.split()
			if words:
				sequences.append(words[1])
	return sequences


def read_fasta_alignment(filename):
	sequences = []
	for data in read_fasta(filename):
		sequences.append(data.seq)
	return sequences


def nucl_count_matrix(sequences):
	nucl_poses = [{} for i in range(max(map(len, sequences)))]

	for seq in sequences:
		for i in range(len(seq)):
			nucl_poses[i][seq[i]] = nucl_poses[i].get(seq[i], 0) + 1

	return nucl_poses


def make_consensus(sequences):
	cons = ''

	nucl_poses = nucl_count_matrix(sequences)

	for nucls in nucl_poses:
		cons += max(nucls.items(), key = lambda x: x[1])[0]

	return faseq(seq = cons, comm = 'algo:cons', name = 'Cons:1')


def count_table(kmers):
	kmer_prop = {}
	for kmer in kmers:
		kmer_prop[kmer] = kmer_prop.get(kmer, 0) + 1
	return kmer_prop


def prop_table(kmers):
	kmer_prop = {}
	all_kmers = len(kmers)
	for kmer in kmers:
		kmer_prop[kmer] = kmer_prop.get(kmer, 0) + 1
	for kmer in kmer_prop:
		kmer_prop[kmer] = kmer_prop[kmer] / all_kmers
	return kmer_prop


def split_to_kmers(sequences, k = 10, steps_back = 3):
	nkmers = len(sequences[0])

	if k != 1:
		nkmers = len(sequences[0]) // k + (1 if len(sequences[0]) % k else 0)

	# inds = [x * k for x in range(nkmers)] + [len(sequences[0]) + 1]
	inds = [x * k for x in range(nkmers)]
	inds[-1] = len(sequences[0]) + 1

	return [['N' for i in range(len(sequences))] for j in range(steps_back)] + [[seq[inds[i-1]:inds[i]] for seq in sequences] for i in range(1, len(inds))]


def entropy(sequences, ind):
	distr = {'A': 1, 'C': 1, 'G': 1, 'T': 1}
	for i in range(len(sequences)):
		distr[sequences[i][ind]] += 1
	distr_sum = sum(distr.values())
	return sum(map(lambda x: x / distr_sum * log(x / distr_sum, 2), distr.values()))


def maximise_information(alpha, beta, noise_level):
	"""
	Given two distributions of equal lengths,
	get those pairs {a from alpha: bs from beta} which
	maximises transmitted information, i.e. drop all pairs / bs
	which has noise >= noise_level in percentage
	"""
	if len(alpha) != len(beta):
		return -1

	# Compute counts of next kmer for each previous kmer.
	prev = {}
	X_only = set()
	for i in range(len(alpha)):
		if alpha[i] not in prev:
			prev[alpha[i]] = {}
		prev[alpha[i]][beta[i]] = prev[alpha[i]].get(beta[i], 0) + 1

		if 'X' in beta[i]:
			X_only.add(alpha[i])

	for xo in X_only:
		has_not_x = False
		xs_items = set()
		for item in prev[xo]:
			if item.count('X') == 0:
				has_not_x = True
			else:
				xs_items.add(item)
		if has_not_x:
			for it in xs_items:
				prev[xo].pop(it)

	# Compute frequency for each kmer pair
	# and drop which have frequency <= noise_level
	# for p in prev:
	# 	nsum = sum(prev[p].values())
	# 	prev[p] = dict(filter(lambda x: x[1] / nsum >= noise_level, prev[p].items()))
	# 	for n in prev[p]:
	# 		prev[p][n] = log(prev[p][n] / nsum, 2)

	# Choose the top 2 nucleotides by frequency.
	# print(prev)
	if len(prev) == 1:
		srt = sorted(prev[("N",)].items(), reverse = True, key = lambda x: x[1])
		if len(srt) > 1:
			if srt[0][1] * .1 <= srt[1][1]:
				prev[("N",)] = {}
				prev[("N",)][srt[0][0]] = srt[0][1]
				prev[("N",)][srt[1][0]] = srt[1][1]
			else:
				prev[("N",)] = {}
				prev[("N",)][srt[0][0]] = srt[0][1]

	else:
		for p in prev:
			top = max(prev[p].items(), key = lambda x: x[1])
			prev[p] = {}
			prev[p][top[0]] = top[1]
			# prev[p] = dict(filter(lambda x: x[0] != top[0], prev[p].items()))

	return prev


def make_kmer_consensus(sequences, k = 10, steps_back = 1, noise_level = .1, verbose = 1):
	"""
	Make a consensuses using kmers.

	Invariant:
		1. We know current good nucleotides.
		2. Using mutual information for each current nucleotide
		we find next nucleotides.

	Arguments:
	sequences -- list of nucleotide sequences
	k -- length of kmers for use
	noise_level -- ???
	"""
	# Column - sequence
	# Row - i-th kmer of sequences
	# Replace all Ns with Xs
	# DICT = set(["A", "C", "G", "T"])
	kmers_mat = split_to_kmers(sequences, k, steps_back)

	# Dictionary of consensus sequences:
	# {previous non-conservative kmer: [(consensus sequence, loglikelihood), ...]}
	prev_kmer_cons = {tuple('N' for i in range(steps_back)): [('', 0)]}
	prev_i = [i for i in range(steps_back)]
	# Each step:
	for cur_i in range(steps_back, len(kmers_mat)):
		if verbose == 2: print('\n', cur_i - steps_back + 1, '/', len(kmers_mat) - steps_back)

		# Make two lists (of equal length) with previous kmers and related current kmers.
		_A = {}

		prev_kmer_list = []
		next_kmer_list = []
		for seq_i in range(len(sequences)):
			# Get previous kmers for sequence with index seq_i
			prev_kmer_tuple = tuple(kmers_mat[prev_i[step]][seq_i] for step in range(steps_back))
			if prev_kmer_tuple in prev_kmer_cons:
				prev_kmer_list.append(prev_kmer_tuple)
				next_kmer_list.append(kmers_mat[cur_i][seq_i])
			_A[kmers_mat[prev_i[0]][seq_i] + kmers_mat[cur_i][seq_i]] = _A.get(kmers_mat[prev_i[0]][seq_i] + kmers_mat[cur_i][seq_i], 0) + 1
		if verbose == 2:
			print(list(sorted(_A.items(), key = lambda x: x[1], reverse = True)))


		# Compute frequency of pairs (previous kmer, current kmer) and drop the noisy ones.
		# {prev: [(next1, loglikelihood1), ...]}
		double_kmers_prop = maximise_information(prev_kmer_list, next_kmer_list, noise_level)
		if verbose == 2: print('-- Kmer chains')
		if verbose == 2: print(double_kmers_prop)

		# Generate new consensus sequences.
		# Check for conservative kmer first.
		new_prev_kmer_cons = {}
		is_conservative = len(set(itertools.chain(*map(lambda x: tuple(x[1].keys()), double_kmers_prop.items())))) == 1
		for prev_kmer_tuple in double_kmers_prop:
			for cur_kmer in double_kmers_prop[prev_kmer_tuple]:
				for i in range(len(prev_kmer_cons[prev_kmer_tuple])):
					work_kmer = prev_kmer_tuple
					if not is_conservative:
						work_kmer = prev_kmer_tuple[1:] + tuple([cur_kmer])
					if work_kmer not in new_prev_kmer_cons:
						new_prev_kmer_cons[work_kmer] = []
					# new_prev_kmer_cons[work_kmer].append((prev_kmer_cons[prev_kmer_tuple][i][0] + cur_kmer,
					# 	prev_kmer_cons[prev_kmer_tuple][i][1] + double_kmers_prop[prev_kmer_tuple][cur_kmer]))
					new_prev_kmer_cons[work_kmer].append((prev_kmer_cons[prev_kmer_tuple][i][0] + cur_kmer,
						double_kmers_prop[prev_kmer_tuple][cur_kmer]))

		prev_kmer_cons = new_prev_kmer_cons
		if not is_conservative:
			prev_i = prev_i[1:] + [cur_i]

		if verbose == 2:
			# print(prev_kmer_cons)
			c_i = 0
			for c in prev_kmer_cons:
				for elem in prev_kmer_cons[c]:
					# print(c_i, ":\t", elem[0])
					c_i += 1
			print('-- Consensuses', end = '\t')
			print(c_i)


	cons_list = []
	comments = []
	for kmer in prev_kmer_cons:
		for i in range(len(prev_kmer_cons[kmer])):
			cons_list.append(prev_kmer_cons[kmer][i])
	cons_list.sort(reverse = False, key = lambda x: x[1])
	for i in range(len(cons_list)):
		comments.append(cons_list[i][0].replace("X", ""))
		# cons_list[i] = cons_list[i][0].replace("X", "").replace("D", ""), cons_list[i][1]


	if verbose:
		print('\nConsensus sequences:')
		for i in range(len(cons_list)):
			print('>Cons:' + str(i + 1), '\tSequences:', cons_list[i][1], '\n', cons_list[i][0], sep = '')

	# return [faseq(seq = cons_list[i][0], comm = 'algo:kmer' + ';k:' + str(k) + ';noise:' + str(noise_level) + ';steps:' + str(steps_back) + ';score:' + str(cons_list[i][1]), name = 'Cons:' + str(i+1)) for i in range(len(cons_list))]
	return [faseq(seq = cons_list[i][0], comm = comments[i], name = 'Cons:' + str(i+1)) for i in range(len(cons_list))]


def find_centre_pos(sequences):
	centre_pos = -1
	best_sum_top2 = 0
	best_second = 0

	tmp = {}

	for pos_i in range(2, len(sequences[0]) - 2):
		d = {'A': 0, 'C': 0, 'T': 0, 'G': 0, 'X': 0, 'D': 0}
		for seq in sequences:
			d[seq[pos_i]] += 1
		d.pop('X')
		d_sorted = list(sorted(d.items(), reverse = True, key = lambda x: x[1]))
		# print(d_sorted)
		# if (d_sorted[0][1] + d_sorted[1][1] > best_sum_top2) and (d_sorted[1][1] >= best_second):
		if (d_sorted[0][1] + d_sorted[1][1] > .5 * len(sequences)) and (d_sorted[1][1] >= best_second):
			best_sum_top2 = d_sorted[0][1] + d_sorted[1][1]
			best_second = d_sorted[1][1]
			centre_pos = pos_i
			tmp = d_sorted

	# input()
	print(tmp)
	# input()
	return centre_pos


def merge_consensuses(forw_seq, back_seq):
	res = []

	if len(forw_seq) == 1:
		forw_seq.append(faseq(name = forw_seq[0].name + "$", comm = forw_seq[0].comm, seq = forw_seq[0].seq))

	if len(back_seq) == 1:
		back_seq.append(faseq(name = back_seq[0].name + "$", comm = back_seq[0].comm, seq = back_seq[0].seq))

	if forw_seq[0].seq[0] == back_seq[0].seq[0]:
		res.append(faseq(name = forw_seq[0].name + "=" + back_seq[0].name, comm = back_seq[0].comm[::-1] + forw_seq[0].comm[1:], seq = back_seq[0].seq[::-1] + forw_seq[0].seq[1:]))
		res.append(faseq(name = forw_seq[1].name + "=" + back_seq[1].name, comm = back_seq[1].comm[::-1] + forw_seq[1].comm[1:], seq = back_seq[1].seq[::-1] + forw_seq[1].seq[1:]))
	else:
		res.append(faseq(name = forw_seq[1].name + "=" + back_seq[0].name, comm = back_seq[0].comm[::-1] + forw_seq[1].comm[1:], seq = back_seq[0].seq[::-1] + forw_seq[1].seq[1:]))
		res.append(faseq(name = forw_seq[0].name + "=" + back_seq[1].name, comm = back_seq[1].comm[::-1] + forw_seq[0].comm[1:], seq = back_seq[1].seq[::-1] + forw_seq[0].seq[1:]))

	return res


# def find_centre_pos(sequences):
# 	# find three neighbourng positions with minimal number of Xs
# 	min_x = len(sequences)
# 	centre_pos = 2
# 	for pos_i in range(2, len(sequences[0]) - 2):
# 		Xs = 0
# 		for seq in sequences:
# 			Xs += seq[pos_i-1 : pos_i+2].count('X')
# 		if Xs < min_x:
# 			min_x = Xs
# 			centre_pos = pos_i
#
# 	return centre_pos
#
#
# def merge_consensuses(forw_seq, back_seq):
# 	res = []
# 	res.append(faseq(name = forw_seq[0].name + "=" + back_seq[0].name, comm = "", seq = forw_seq[0].seq + back_seq[0].seq[::-1][3:]))
# 	if len(forw_seq) == 2:
# 		res.append(faseq(name = forw_seq[1].name + "=" + back_seq[0].name, comm = "", seq = forw_seq[1].seq + back_seq[0].seq[::-1][3:]))
# 	if len(back_seq) == 2:
# 		res.append(faseq(name = forw_seq[0].name + "=" + back_seq[1].name, comm = "", seq = forw_seq[0].seq + back_seq[1].seq[::-1][3:]))
# 	if len(forw_seq) == 2 and len(back_seq) == 2:
# 		res.append(faseq(name = forw_seq[1].name + "=" + back_seq[1].name, comm = "", seq = forw_seq[1].seq + back_seq[1].seq[::-1][3:]))
# 	return res


if __name__ == '__main__':
	"""
	1: Make consensus from the ClustalW file:
		$python3 cons.py --mode clw -i <input file> -o <optional output file>

	2: Make consensus from the FASTA file:
		$python3 cons.py --mode fasta -i <input file> -o <optional output file>

	3: Make consensus using kmers:
		$python3 cons.py --mode kmer -i <input file> -k <optional k (default 10)> -o <optional output file>
	"""
	p = argparse.ArgumentParser()
	p.add_argument('-i', help = 'input file', type = str)
	p.add_argument('-o', help = 'optional output file', type = str)
	p.add_argument('--mode', help = 'consensus sequence input file mode: clw, fasta or kmer', type = str)
	p.add_argument('-k', help = 'kmer size', default = 1, type = int)
	p.add_argument('--noise', help = 'noise level in proportion', default = .1, type = float)
	p.add_argument('--steps', help = 'number of steps to look back', default = 1, type = int)
	p.add_argument('--verbose', help = '2 for full output of the process, 1 for consensus only output, 0 for no output', default = 1, type = int)
	args = p.parse_args()

	cons = []
	if args.mode == 'clw':
		cons = make_consensus(read_clustalw(args.i))
	elif args.mode == 'fasta':
		cons = make_consensus(read_fasta_alignment(args.i))
	elif args.mode == 'kmer':
		lens = {}
		sequences = []
		for fs in read_fasta(args.i):
			sequences.append(fs.seq)
			lens[len(fs.seq)] = lens.get(len(fs.seq), 0) + 1
		LEN = max(lens.items(), key = lambda x: x[1])[0]
		cons = make_kmer_consensus([x for x in sequences if len(x) == LEN], k = args.k, steps_back = args.steps, noise_level = args.noise, verbose = args.verbose)

	elif args.mode == 'kmer-aligned':
		allele_list = {}
		cons = []
		# for d in read_fasta(args.i):
		for d in FASTQParser(args.i):
			# nameall = d.name.split('///')[0][1:]
			nameall = d.name.split('///')[0]
			if nameall not in allele_list:
				allele_list[nameall] = []
			allele_list[nameall].append(d.seq.replace("N", "X"))
		print("Found alleles:")
		print(", ".join(allele_list.keys()))
		with open(args.i + '.' + 'with_dels' + '.txt', 'w') as f:
			for allele in allele_list:
				# if allele == "HLA_A_nuc":
				# if allele in ["HLA_A_nuc", "HLA_B_nuc", "HLA_C_nuc", "DPB1_nuc"]:
				# if allele in ["HLA_A_nuc", "HLA_B_nuc", "HLA_C_nuc", "DQB1_nuc", "DPB1_nuc", "DRB1_nuc"]:
				if True:
					print("============================")
					print("New allele:", allele, "(", len(allele_list[allele]), " reads)")

					centre_pos = find_centre_pos(allele_list[allele])
					print("Centre position =", centre_pos, "( length = ", len(allele_list[allele][0]), ")")

					forw_seq = [x[centre_pos:] for x in allele_list[allele]]
					forw_cons_res = make_kmer_consensus(forw_seq, k = args.k, steps_back = args.steps, noise_level = args.noise, verbose = args.verbose)

					back_seq = [x[:centre_pos+1][::-1] for x in allele_list[allele]]
					back_cons_res = make_kmer_consensus(back_seq, k = args.k, steps_back = args.steps, noise_level = args.noise, verbose = args.verbose)

					cons_res = merge_consensuses(forw_cons_res, back_cons_res)

					for cn in cons_res:
						f.write(allele + '\t' + cn.name + '\n' + cn.seq + '\n')
						new_cn = faseq(name = allele + '(len=' + str(len(cn.seq)) + ')' + '|' + cn.name, comm = "", seq = cn.seq.replace('X', '').replace('D', ''))
						cons.append(new_cn)
					print(allele, ':', len(cons_res), 'consensuses (', len(allele_list[allele]), "sequences )")
					print()


	# if args.mode != 'kmer' or args.mode != 'kmer-aligned':
	# 	args.o = args.i + '.cons'
	# else:
	# 	args.o = args.i + '.k' + str(args.k) + '.noise' + str(args.noise) + '.steps' + str(args.steps) + '.cons'
	directory = '/'.join(args.o.split('/')[:-1])
	if not os.path.exists(directory):
		os.makedirs(directory)

	write_fasta(args.o, cons)
