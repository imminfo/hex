from collections import deque

class Parser:
	''' 
	Abstract parser for some filetypes with buffering.
	'''

	def __init__(self, filepath, buffer_size):
		'''
		Parser(filepath, buffer_size)

		Open the parser of file with name filepath and with the given buffer size.
		'''
		self._buffer = deque(maxlen = buffer_size)
		self._buffer_size = buffer_size
		self.open(filepath)

	def __iter__(self):
		'''
		iter(P) -> return itself
		'''
		return self

	def __next__(self):
		'''
		next(P) -> data or None
		'''
		# Try to get data.
		data = self.read()
		# If no data left in both buffer and input file, stop iteration.
		if data:
			return data
		else:
			raise StopIteration

	def _parse_data(self):
		'''
		P._parse_data() -> one data unit or None if no data left in the file.

		Parse the input file and return one data unit or None if no data left.
		'''
		raise NotImplementedError

	def _parse_to_buffer(self):
		'''
		P._parse_to_buffer() -> None

		Parse data from the input file to the buffer or raise BufferError if no data left in the file.
		'''
		for i in range(self._buffer_size):
			data = self._parse_data()
			if data:
				self._buffer.append(data)
			else:
				break;

	def read(self):
		'''
		P.read() -> one data unit or None if no data left in file.

		Get one data record or None if no data left.
		'''
		# Try to get data from the buffer.
		# If no data left in buffer, try to parse new data to the buffer.
		if len(self._buffer) == 0:
			self._parse_to_buffer()
			# If no data left in file, return void.
			if len(self._buffer) == 0:
				return None
		# Return the oldest data from file.
		return self._buffer.popleft()

	def open(self, filepath):
		'''
		P.open(filepath)

		Open the input file for parsing.
		'''
		self._file = open(filepath)

	def is_eof(self):
		'''
		Check if the parser has reached the end of the inpiut file.
		'''
		if not self._buffer:
			self._parse_to_buffer()
			if not self._buffer:
				return False
		return True

	def data(self):
		'''
		'''
		tmp = self._buffer[:]
		self._buffer.clear()
		self._parse_to_buffer()
		return tmp

	def close(self):
		'''
		P.close()

		Close the input file.
		'''
		self._file.close()
		self._buffer.clear()


class PairendParser:
	'''
	Abstract parser for pair of files based with underlying Parser objects.
	'''

	def __init__(self, filepath1, filepath2, buffer_size):
		'''
		P(filepath1, filepath2, buffer_size)

		Open the given two files with the given buffer size for underlying Parser objects.
		'''
		self._buffer_size = buffer_size
		self.open(filepath1, filepath2)

	def __iter__(self):
		'''
		iter(P) -> return itself
		'''
		return self

	def __next__(self):
		'''
		next(P) -> data or None
		'''
		# Try to get data.
		data = self.read()
		# If no data left, stop iteration.
		if data:
			return data
		else:
			raise StopIteration

	def _open_parser(self, filepath, buffer_size):
		'''
		P._open_parser(filepath, buffer_size) -> Parser

		Open and return Parser to the given file with the given buffer_size
		'''
		raise NotImplementedError

	def read(self):
		'''
		P.read() -> (left_data, right_data)

		Get the pair of data where is left_data from the first file and right_data from the second file.
		'''
		raise NotImplementedError

	def open(self, filepath1, filepath2):
		'''
		P.open(self, filepath1, filepath2)

		Open the given two files.
		'''
		self._parser_alpha = self._open_parser(filepath1, self._buffer_size)
		self._parser_beta = self._open_parser(filepath2, self._buffer_size)

	def close(self):
		'''
		P.close()

		Close writers.
		'''
		self._parser_alpha.close()
		self._parser_beta.close()

	def is_eof(self):	
		'''
		Check if the parser has reached the end of the inpiut file.
		'''
		return self._parser_alpha.is_eof() and self._parser_beta.is_eof()