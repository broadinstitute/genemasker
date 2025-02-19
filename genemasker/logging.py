import sys
import logging

class MultiLogger(logging.StreamHandler):
	def __init__(self, stream, file_path):
		super().__init__(stream)
		self.file = open(file_path, "w", buffering=1)

	def emit(self, record):
		log_entry = self.format(record)
		self.stream.write(log_entry + "\n")
		self.stream.flush()
		self.file.write(log_entry + "\n")
		self.file.flush()

	def close(self):
		self.file.close()
		super().close()

def setup_logger(logfile):
	logger = logging.getLogger()
	logger.setLevel(logging.INFO)
	
	logger_handler = MultiLogger(sys.stdout, logfile)
	
	formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s", datefmt = "%Y-%m-%d %H:%M:%S %Z")
	logger_handler.setFormatter(formatter)
	
	logger.addHandler(logger_handler)

	return logger, logger_handler
