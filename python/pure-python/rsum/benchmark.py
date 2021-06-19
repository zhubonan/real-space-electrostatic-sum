"""
Benchmarking
"""

from .test_rsum import _al_base as al
from .test_rsum import _nacl_base as nacl
from .rsum import energy
from time import time

def timeit(timeout=10):
    """
    Run tests for a function
    """
    def _wrapped(f):
        def __wrapped(*args, **kwargs):
            start = time()
            count = 0
            while time() - start < timeout:
                res = f(*args, **kwargs)
                count += 1
            finished = time()
            elapsed = finished - start
            avg = elapsed / count
            print(f'Average execuation time: {avg:.2g} s, over {count} runs')
            return res
        return __wrapped
    return _wrapped


energy = timeit(timeout=5)(energy)

if __name__ == '__main__':
    energy(*al()[:-1])
    energy(*nacl()[:-1])