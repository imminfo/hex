from collections import namedtuple
from writer import *
from parser import *

FASTQData = namedtuple('FASTQData', 'name attr seq qual')

DEFAULT_FASTAQP_BUFFER_SIZE = 100000
DEFAULT_FASTAQW_BUFFER_SIZE = 100000


def subseq(fastq_data, n, last = 0):
	'''
	Return the subsequence of the sequence of the given FASTQData starting on the pos n and ending in the pos last and return new FASTQData.
	'''
	if last:
		return FASTQData(seq = fastq_data.seq[n:last], attr = fastq_data.attr, qual = fastq_data.qual[n:last], name = fastq_data.name)
	else:
		return FASTQData(seq = fastq_data.seq[n:], attr = fastq_data.attr, qual = fastq_data.qual[n:], name = fastq_data.name)


class FASTQParser(Parser):
	'''
	Class for parsing files written in the FASTQ format.
	'''

	def __init__(self, filepath, buffer_size = DEFAULT_FASTAQP_BUFFER_SIZE):
		'''
		FASTQParser(filename, buffer_size) -> new parser
		filename - name of the file for parsing.
		buffer_size - size of the inner buffer.
		'''
		Parser.__init__(self, filepath, buffer_size = DEFAULT_FASTAQP_BUFFER_SIZE)

	def _parse_data(self):
		'''
		P._parse_data() -> FASTQData or None

		Load one data record from the file or return None if EOF.
		'''
		line = self._file.readline().strip()
		if line:
			words = tuple(line.split())
			if len(words) == 1:
				rec_name = words[0]
				rec_attr = ''
			else:
				rec_name, rec_attr = words
			if rec_name[0] == '>':
				rec_name = rec_name[1:]
			rec_seq = self._file.readline().strip()
			self._file.readline()
			rec_qual = self._file.readline().strip()
			return FASTQData(rec_name, rec_attr, rec_seq, rec_qual)
		else:
			return None


class PairendFASTQParser(PairendParser):
	'''
	Class for parsing FASTQ files with mates - pairs of FASTQData with similar names from the given pair of files.
	After processing all mates returns data records with no mate in the other file.
	Parses files with the FASTQParsers.
	'''

	def __init__(self, filepath1, filepath2, buffer_size = DEFAULT_FASTAQP_BUFFER_SIZE, max_mates = -1):
		'''
		PairwiseFASTQParser(filepath1, filepath2, buffer_size) -> new PairwiseFASTQParser		
		filepath1, filepath2 - 2 strs naming two filenames for parsing.
		buffer_size - buffer size of the FASTQParsers.
		'''
		PairendParser.__init__(self, filepath1, filepath2, buffer_size)
		# Dictionary for mates {str: (number of file, FASTQData)}. If the string S is in the dictionary,
		# then it has a mate in the self._mates[S] file (0 for the left or alpha file, 1 for the right or beta file)
		self._mates = {}
		self._nodata_left = False
		self._nodata_right = False

	def _open_parser(self, filepath, buffer_size):
		'''
		P._open_parser(filepath, buffer) -> FASTQParser

		Open the FASTQParser to the given filepath and with the given buffer_size.
		'''
		return FASTQParser(filepath, buffer_size)

	def open(self, filepath1, filepath2):
		'''
		P.open(filpath1, filepath2)

		Open the given files for parsing.
		'''
		PairendParser.open(self, filepath1, filepath2)
		self._nodata_left = False
		self._nodata_right = False

	def read(self):
		'''
		P.read() -> (FASTQData from the left or alpha file, FASTQData from the right or beta file)

		Return pair of FASTQData corresponding to the given files,
		(None, FASTQData) or (FASTQData, None) if no corresponding mate has been found,
		or (None, None) if no data left.
		'''
		while not (self._nodata_left and self._nodata_right):
			if not self._nodata_left:
				tmp = left_data = self._parser_alpha.read()
				if tmp:
					left_name = left_data.name

					if left_name in self._mates:
						return (left_data, self._mates.pop(left_name)[1])
					else:
						self._mates[left_name] = (0, left_data)
				else:
					self._nodata_left = True

			if not self._nodata_right:
				if tmp:
					right_data = self._parser_beta.read()
					right_name = right_data.name

					if right_name in self._mates:
						return (self._mates.pop(right_name)[1], right_data)
					else:
						self._mates[right_name] = (1, right_data)			
				else:
					self._nodata_right = True

		rec = [None, None]
		if self._nodata_left and self._nodata_right:
			try:
				item = self._mates.popitem()
				rec[item[1][0]] = item[1][1]
				rec[1 - item[1][0]] = None
			except KeyError:
				return None
		return tuple(rec)
		

class FASTQWriter(Writer):
	'''
	Class for writing files in FASTQ format. Receives FASTQData named tuples and writes them to the file.
	'''

	def __init__(self, filepath, buffer_size = DEFAULT_FASTAQW_BUFFER_SIZE):
		'''
		FASTQWriter(filename, buffer_size) -> new FASTQWriter
		filepath - name of the file to writing.
		buffer_size - size of the inner buffer.
		'''
		Writer.__init__(self, filepath, buffer_size)

	def _write_data(self, some_data):
		'''
		W._write_data(some_data)

		Write the given FASTQData to the output file.
		'''
		if some_data:
			self._file.write('@' + some_data.name + ' ' + some_data.attr + '\n' + some_data.seq + '\n' + '+\n' + some_data.qual + '\n')


class PairendFASTQWriter(PairendWriter):
	'''
	Class for writing the given FASTQData named tuple pairs (i.e., mates) to the given output files.
	'''

	def __init__(self, filepath1, filepath2, buffer_size = DEFAULT_FASTAQW_BUFFER_SIZE):
		'''
		PairwiseFASTQWriter(filenames, writers_buffer_size) -> new PairwiseFASTQWriter
		filepath1, filepath2 - pair of the output files.
		buffer_size - buffer of FASTQWriters.
		'''
		PairendWriter.__init__(self, filepath1, filepath2, buffer_size)

	def _open_writer(self, filepath, buffer_size):
		'''
		W.open(filepath, buffer_size)

		Open the output file for writing with the given filename.
		'''
		return FASTQWriter(filepath, buffer_size)


def filter_mates(f, input_l, input_r):
	'''
	filter_mates(f, input_l, input_r)

	Read mates from the given files and write filtered to the output files.
	Given filter is a lambda or callable object that receives order of the read FASTQData.
	'''
	parser = PairendFASTQParser(input_l, input_r)
	writer = PairendFASTQWriter(input_l + '.filtered', input_r + '.filtered')

	order = 1
	count = 0
	fcount = 0
	for d in parser:
		if f(order):
			writer.write(d)
			fcount += 1
		order += 1
		count += 1 if d[0] else 0
		count += 1 if d[1] else 0
	print('Done.')
	print(count, 'mates has been processed.')
	print(fcount, 'pairs of mates has been written to the output files:')
	print(input_l + '.filtered')
	print(input_r + '.filtered')

	parser.close()
	writer.close()


if __name__ == '__main__':
	p = PairendFASTQParser('data/P83_1_ATCACG_L003_R1_001_sample.fastq', 'data/P83_1_ATCACG_L003_R2_001_sample.fastq')
	for d in p:
		print(d[0], d[1], sep = ' | ')