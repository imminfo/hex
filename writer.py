from collections import deque

class Writer:
	'''
	Abstract class for writing to some file with buffering.
	'''

	def __init__(self, filepath, buffer_size):
		'''
		Writer(filepath, buffer_size)

		Open the writer to the filepath with the given buffer_size.
		'''
		self._buffer = deque(maxlen = buffer_size)
		self._buffer_size = buffer_size
		self.open(filepath)

	def _write_data(self, some_data):
		'''
		W._write_data(some_data)

		Write the given data to the output file.
		'''
		raise NotImplementedError

	def _write_buffer(self):
		'''
		W._write_buffer()

		Write all data from the buffer to the output file.
		'''
		for some_data in self._buffer:
			self._write_data(some_data)
		self._buffer.clear()

	def write(self, some_data):
		'''
		W.write(some_data)

		Write the given data to the buffer and (later) to the output file.
		'''
		if len(self._buffer) >= self._buffer_size:	# Just in case.
			self._write_buffer()
			self._buffer.append(some_data)
		elif len(self._buffer) == (self._buffer_size - 1):
			self._buffer.append(some_data)
			self._write_buffer()
		else:
			self._buffer.append(some_data)

	def open(self, filepath):
		'''
		W.open(filepath)

		Open the given file.
		'''
		self._file = open(filepath, 'w')

	def close(self):
		'''
		W.close()

		Write buffer and close the output file.
		'''
		self._write_buffer();
		self._file.close()		


class PairendWriter:
	'''
	Abstract class for writing pairs of data to two output files with underlying Writer objects.
	'''

	def __init__(self, filepath1, filepath2, buffer_size):
		'''
		PairendWriter(filepath1, filepath2, buffer_size)

		Open the pairend writer to the given two files with the buffer size of underlying Writer objects.
		'''
		self._buffer_size = buffer_size
		self.open(filepath1, filepath2)

	def _open_writer(self, filepath, buffer_size):
		'''
		W._open_writer(filepath, buffer_size) -> Writer

		Open Writer to the output filepath with buffer size and return it.
		'''
		raise NotImplementedError

	def write(self, some_data_pair):
		'''
		W.write(some_data_pair)

		Write the given pair of data to corresponding writers.
		'''
		self._writer_alpha.write(some_data_pair[0])
		self._writer_beta.write(some_data_pair[1])

	def open(self, filepath1, filepath2):
		'''
		W.open(filepath1, filepath2)

		Open Writer objects to the output files.
		'''
		self._writer_alpha = self._open_writer(filepath1, self._buffer_size)
		self._writer_beta = self._open_writer(filepath2, self._buffer_size)

	def close(self):
		'''
		W.close()

		Close the output files.
		'''
		self._writer_alpha.close()
		self._writer_beta.close()