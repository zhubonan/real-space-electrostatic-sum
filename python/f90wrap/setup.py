from numpy.distutils.core import Extension
import os

source_folder = os.path.join(os.path.dirname(__file__), '../source')

ext1 = Extension(name='_rsumf',
                 sources=[
                     'f90wrap_real_space_electrostatic_sum.f90',
                     os.path.join(source_folder,
                                  'real_space_electrostatic_sum.f90')
                 ])

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(name='rsumf',
          install_requires=['f90wrap>=0.2.3'],
          description="Real-space electrostatic summation",
          ext_modules=[ext1])