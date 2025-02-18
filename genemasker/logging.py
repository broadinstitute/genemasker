import sys
import logging

# Define a custom class to log to stdout and file simultaneously
class MultiLogger(logging.StreamHandler):
	def __init__(self, stream, file_path):
		super().__init__(stream)
		self.file = open(file_path, "w", buffering=1)  # Line-buffered mode (minimizes buffering)

	def emit(self, record):
		log_entry = self.format(record)
		self.stream.write(log_entry + "\n")  # Log to stdout
		self.stream.flush()  # Ensure immediate output to stdout
		self.file.write(log_entry + "\n")  # Log to file
		self.file.flush()  # Ensure immediate output to file

	def close(self):
		self.file.close()
		super().close()

def setup_logger(logfile):
	# Create a logger
	logger = logging.getLogger()
	logger.setLevel(logging.INFO)
	
	# Create a handler that writes to both stdout and a file with no buffering
	logger_handler = MultiLogger(sys.stdout, logfile)
	
	# Define a simple log format
	formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s", datefmt = "%Y-%m-%d %H:%M:%S")
	logger_handler.setFormatter(formatter)
	
	# Add handler to the logger
	logger.addHandler(logger_handler)

	return logger, logger_handler
