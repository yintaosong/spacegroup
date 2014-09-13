from __future__ import print_function, division, absolute_import
import numpy as np
import numpy.linalg as la

try:
    from numba import jit
except ImportError:
    def jit(sig=''):
        def jit_decorator(func):
            def func_wrapper(*args, **kwargs):
                return func(*args, **kwargs)
            return func_wrapper
        return jit_decorator