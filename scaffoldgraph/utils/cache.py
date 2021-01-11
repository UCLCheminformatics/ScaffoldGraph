"""
scaffoldgraph.utils.cache
"""

from collections import OrderedDict
from operator import eq as _eq


class Cache(OrderedDict):
    """A basic implementation of an LRU cache using OrderedDict.

    Adapted (slightly) from the collections ``OrderedDict``
    documentation.

    .. _collections OrderedDict Documentation:
   https://docs.python.org/3/library/collections.html#collections.OrderedDict

    """
    def __init__(self, maxsize=None, *args, **kwargs):
        """
        Parameters
        ----------
        maxsize : int, None, optional
            Set the maximum size of the cache, if None the cache
            has no size limitation. The default is None.
        *args
            Variable length argument list.
            Passed to OrderedDict.
        **kwargs
            Arbitrary keyword arguments.
            Passed to OrderedDict.

        """
        self._maxsize = maxsize
        super(Cache, self).__init__(*args, **kwargs)

    @property
    def maxsize(self):
        """int: The maximum size of the cache."""
        return self._maxsize

    def __getitem__(self, key):
        value = super().__getitem__(key)
        self.move_to_end(key)
        return value

    def __setitem__(self, key, value):
        super().__setitem__(key, value)
        if self.maxsize and len(self) > self.maxsize:
            oldest = next(iter(self))
            del self[oldest]

    def __eq__(self, other):
        if isinstance(other, Cache):
            return dict.__eq__(self, other) and all(map(_eq, self, other))
        return dict.__eq__(self, other)

    def __repr__(self):
        return '{}(maxsize={})'.format(
            self.__class__.__name__,
            self.maxsize
        )
