"""
scaffoldgraph.utils.logging

Utilities for dealing with rdkit logging.
"""

import functools
import warnings

from rdkit import __version__ as rdversion
from rdkit import RDLogger, rdBase


DEFAULT_RDLOGGER_STATUS = {
    'rdApp.debug': True,
    'rdApp.info': True,
    'rdApp.warning': True,
    'rdApp.error': True
}

QUIET_RDLOGGER_STATUS = {
    'rdApp.debug': False,
    'rdApp.info': False,
    'rdApp.warning': False,
    'rdApp.error': True
}

UNKNOWN_RDLOGGER_STATUS = DEFAULT_RDLOGGER_STATUS


def get_rdlogger_status():
    """dict : Return the status of the rdlogger."""
    status_dict = {}
    if rdversion < '2020.09.01':
        warnings.warn('Failed to get status of rdlogger')
        return UNKNOWN_RDLOGGER_STATUS
    for status in rdBase.LogStatus().split('\n'):
        level, state = status.split(':')
        status_dict[level] = True if state == 'enabled' else False
    return status_dict


def set_rdlogger_status(status_dict):
    """Set the state of the rdlogger."""
    for level, state in status_dict.items():
        if state is True:
            rdBase.EnableLog(level)
        else:
            rdBase.DisableLog(level)


def set_rdlogger_quiet():
    """Set the rdlogger to quiet status."""
    set_rdlogger_status(QUIET_RDLOGGER_STATUS)


def reset_rdlogger():
    """Reset the rdlogger status to default."""
    set_rdlogger_status(DEFAULT_RDLOGGER_STATUS)


def set_rdlogger_debug_status(status):
    """Set status of rdlogger: debug."""
    set_rdlogger_status({'rdApp.debug': status})


def set_rdlogger_info_status(status):
    """Set status of rdlogger: info."""
    set_rdlogger_status({'rdApp.info': status})


def set_rdlogger_warning_status(status):
    """Set status of rdlogger: warning."""
    set_rdlogger_status({'rdApp.warning': status})


def set_rdlogger_error_status(status):
    """Set status of rdlogger: error."""
    set_rdlogger_status({'rdApp.error': status})


def suppress_rdlogger(
        suppress_info=True,
        suppress_warning=True,
        suppress_error=True,
        suppress_debug=True
):
    """Decorator for controlling the output level of the rdkit logger.

    Useful for supressing the output of noisy functions related to
    the rdkit logger. The previous status of the logger is returned
    after the function has been executed.

    Parameters
    ----------
    suppress_info : bool, optional
        Suppress logs from rdApp.info. The default is True.
    suppress_warning : bool, optional
        Suppress logs from rdApp.warning. The default is True.
    suppress_error : bool, optional
        Suppress logs from rdApp.error. The default is True.
    suppress_debug : bool, optional
        Suppress logs from rdApp.debug. The default is True.

    Returns
    -------
    decorator : function

    Notes
    -----
    The prior state of the logger can only be returned in the newer
    versions of rdkit (>= '2020.09.01'). In previous versions the
    logger status is returned to its default state.

    """
    rdlogger, altered_status = RDLogger.logger(), {}
    altered_status['rdApp.info'] = not suppress_info
    altered_status['rdApp.warning'] = not suppress_warning
    altered_status['rdApp.error'] = not suppress_error
    altered_status['rdApp.debug'] = not suppress_debug

    def decorator(func):
        @functools.wraps(func)
        def wrap_suppress(*args, **kwargs):
            # rdkit version compatability.
            prior_status = DEFAULT_RDLOGGER_STATUS
            if rdversion >= '2020.09.01':
                prior_status = get_rdlogger_status()
            set_rdlogger_status(altered_status)
            try:  # restore status of rdlogger on failure.
                result = func(*args, **kwargs)
            except Exception as e:
                set_rdlogger_status(prior_status)
                raise e
            set_rdlogger_status(prior_status)
            return result
        return wrap_suppress
    return decorator
