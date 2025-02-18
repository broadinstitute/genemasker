import os
import time
import psutil
import functools

def memory_usage():
	"""Returns the current memory usage in MB."""
	process = psutil.Process(os.getpid())
	return process.memory_info().rss / 1024 ** 2


def resource_tracker(logger):
	"""Decorator to track memory usage before and after a function call."""
	def decorator(func):
		@functools.wraps(func)
		def wrapper(*args, **kwargs):
			start_time = time.time()
			mem_before = memory_usage()
			logger.info(f"[MEMORY] {func.__name__} - Memory before: {mem_before:.2f} MB")
			result = func(*args, **kwargs)
			elapsed_time = time.time() - start_time
			mem_after = memory_usage()
			logger.info(f"[TIME] {func.__name__} completed in {elapsed_time:.2f} seconds")
			logger.info(f"[MEMORY] {func.__name__} - Memory after: {mem_after:.2f} MB")
			logger.info(f"[MEMORY] {func.__name__} - Memory change: {mem_after - mem_before:.2f} MB")
			return result
		return wrapper
	return decorator
