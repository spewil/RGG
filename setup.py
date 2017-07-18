from distutils.core import setup
from Cython.Build import cythonize
from setuptools import setup, find_packages, Extension
from Cython.Distutils import build_ext

# for profiling code 

#this is only for pyx / cythonized files 

# basic way:

# setup(
#     ext_modules = cythonize("network_generation/HopcroftKarp.pyx")
# )

# advanced way that ensures the .so files end up in the correct location:

ext_modules=[Extension("network_generation.generation", # location of the resulting .so
             ["network_generation/generation.pyx"]) 
			]

setup(name='package',
      packages=find_packages(),
      cmdclass = {'build_ext': build_ext},
      ext_modules = ext_modules,
     )

# ext_modules=[Extension("network_generation.HopcroftKarp", # location of the resulting .so
#              ["network_generation/HopcroftKarp.pyx"]) 
# 			]

# setup(name='package',
#       packages=find_packages(),
#       cmdclass = {'build_ext': build_ext},
#       ext_modules = ext_modules,
#       compiler_directives={'linetrace': True},
#      )

#Some info about compiling cython:

#http://stackoverflow.com/questions/8106258/cc1plus-warning-command-line-option-wstrict-prototypes-is-valid-for-ada-c-o
#http://stackoverflow.com/questions/18423512/calling-c-code-from-python-using-cython-whith-the-distutilis-approach
#https://github.com/cython/cython/wiki/WrappingCPlusPlus
#http://stackoverflow.com/questions/9846182/errors-when-compiling-first-cython-program


