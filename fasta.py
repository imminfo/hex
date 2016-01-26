import sys
from collections import namedtuple


faseq = namedtuple('faseq', 'name comm seq')


def fastq2fasta(filename, outname = ''):
	if not outname:
		outname = filename.replace('.fastq', '.fasta')
		outname = outname.replace('_R1', '')

	with open(outname, 'w') as out:
		with open(filename) as inp:
			write_line = False
			for line in inp:
				if write_line:
					out.write(line)
					write_line = False
				elif line.find(' ') != -1:
					out.write('>' + line[1:])
					write_line = True

	return outname


def read_fasta(filename):
	with open(filename) as inp:
		name = ''
		comm = ''
		seq = ''
		for line in inp:
			if line[0] == '>':
				if name:
					yield faseq(name, comm, seq)
				words = line.split()
				name = words[0][1:]
				comm = ''
				if len(words) > 1:
					comm = words[1]
				seq = ''
			else:
				seq += line.strip()
		yield faseq(name, comm, seq)


def search_fasta(filename, targets):
	return [x for x in read_fasta(filename) if x.name in targets]


def write_fasta(filename, faseq_list, sep = ' ', length = 80):
	with open(filename, 'w') as out:
		for fs in faseq_list:
			out.write('>' + fs.name + sep + fs.comm + '\n')
			i = 0
			while i < len(fs.seq):
				out.write(fs.seq[i:i + length] + '\n')
				i += length


def filter_fasta(filename, target_file, out_postfix = ''):
	'''
	Read file with targets, search for them in filename and write them to the output file.
	'''
	if not out_postfix:
		out_postfix = target_file + '.search'
	outfile = filename + out_postfix
	targets = []
	with open(target_file) as tf:
		targets = tf.read().split()
	print('Searching for', len(targets), 'sequences...')
	fl = search_fasta(filename, targets)
	print('Found:', len(fl), 'sequences.')
	write_fasta(outfile, fl)


if __name__ == '__main__':
	"""
	1: Transform .fastq file to .fasta file:
		$python3 fasta.py fastq <input .fasta> <optional output name of .fastq>

	2: Search for targets:
		$python3 fasta.py search <input .fasta> <file with names> <optional postfix>

	3: Switch names and attributes in .fasta file:
		$python3 fasta.py switch <input .fasta>
	"""
	if sys.argv[1] == 'fastq':
		output = ''
		if len(sys.argv) > 3:
			output = sys.argv[3]
		fastq2fasta(sys.argv[2], output)
	elif sys.argv[1] == 'search':
		postfix = ''
		if len(sys.argv) > 4:
			postfix = sys.argv[4]
		filter_fasta(sys.argv[2], sys.argv[3], postfix)
	elif sys.argv[1] == 'switch':
		write_fasta(sys.argv[2], [faseq(name = ','.join(map(lambda y: 'HLA:HLA' + y, x.comm.split(','))), comm = x.name, seq = x.seq) for x in read_fasta(sys.argv[2])])
