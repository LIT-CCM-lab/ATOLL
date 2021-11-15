import sys
import os
import inspect
import logging

def get_logger(caller_module_name):
    log_formatter = logging.Formatter("[%(levelname)s]\t%(module)s\t%(message)s")

    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setLevel(logging.INFO)
    stream_handler.setFormatter(log_formatter)
    logger = logging.getLogger(caller_module_name)

    logger.addHandler(stream_handler)
    logger.setLevel(logging.INFO)

    return logger