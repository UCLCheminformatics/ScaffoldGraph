"""
scaffoldgraph.scripts.misc
"""

import logging
import os

import tqdm


class TqdmHandler(logging.Handler):
    """Logging handler for use with tqdm (used in CLI)."""

    def __init__(self, level=logging.NOTSET):
        super().__init__(level)

    def emit(self, record):
        try:
            msg = self.format(record)
            tqdm.tqdm.write(msg)
            self.flush()
        except (KeyboardInterrupt, SystemExit):
            raise
        except Exception:
            self.handleError(record)


def file_format(path):
    """Determine an input file format from a path."""
    split_path, extension = os.path.splitext(path)
    if extension == '.sdf':
        return 'SDF', False
    elif extension == '.smi':
        return 'SMI', False
    elif extension == '.gz' or extension == '.gzip':
        new_extension = file_format(split_path)
        if new_extension[0] is not None:
            return new_extension[0], True
        else:
            return None, False
    else:
        return None, False
