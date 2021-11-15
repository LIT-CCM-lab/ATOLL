import sys
import os

import initialize as init

logger = init.get_logger(__file__)

class FileChecker:
	_SUPPORTED_FILEEXTS = ()

	def __init__(self, filepath):
		self.filepath = filepath
		self.check_path()    

	def check_path(self):
		self.dirpath = os.path.dirname(self.filepath)
		self.filename = os.path.basename(self.filepath)
		self.fileprefix, self.fileext = os.path.splitext(self.filename)

		if not os.path.exists(self.filepath):
			logger.critical(f'{os.path.basename(self.filepath)} path is not found')

		if not os.path.isfile(self.filepath):
			logger.critical(f'{os.path.basename(self.filepath)} path is not a file')

		if not os.access(self.filepath, os.R_OK):
			logger.critical(f'{os.path.basename(self.filepath)} is not readable')
		
		if self.fileext not in self._SUPPORTED_FILEEXTS:
			logger.critical(f'{self.fileext} file extention from file {os.path.basename(self.filepath)} is not a supported format')


class DirChecker:
	def __init__(self, dirpath):
		self.dirpath = dirpath
		self.check_path()

	def check_path(self):
		if not os.path.exists(self.dirpath):
			logger.critical(f'{os.path.basename(self.dirpath)} path is not found')

		if not os.path.isdir(self.dirpath):
			logger.critical(f'{os.path.basename(self.dirpath)} path is not a directory')

		if not os.access(self.dirpath, os.R_OK):
			logger.critical(f'{os.path.basename(self.dirpath)} directory is not readable')

		if not os.access(self.dirpath, os.X_OK):
			logger.critical(f'{os.path.basename(self.dirpath)} directory is not accessible')