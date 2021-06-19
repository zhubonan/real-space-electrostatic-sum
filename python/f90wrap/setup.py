import setuptools
from numpy.distutils.core import Extension
from numpy.distutils.core import setup
import os

source_folder = os.path.join(os.path.dirname(__file__), '../../source')

ext1 = Extension(name='rsum.extension._rsumf',
                 sources=[
                     'rsum/extension/f90wrap_real_space_electrostatic_sum.f90',
                     os.path.join(source_folder,
                                  'real_space_electrostatic_sum.f90')
                 ])

setup(name='rsum',
      install_requires=['f90wrap>=0.2.3'],
      description="Real-space electrostatic summation",
      packages=['rsum', 'rsum.extension'], 
      ext_modules=[ext1])
