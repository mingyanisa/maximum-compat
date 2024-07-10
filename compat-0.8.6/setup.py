#!/usr/bin/env python


from setuptools import setup, Extension
from distutils.sysconfig import get_python_inc,get_python_lib
print(get_python_inc()) 
print(get_python_lib())
setup(name='compat',
      version='0.8.6',
      description='Maximum compatibility phylogenetic algorithm',
      author='Joshua L. Cherry',
      author_email='jcherry@ncbi.nlm.nih.gov',
      url='https://ftp.ncbi.nlm.nih.gov/pub/jcherry/compat/',
      py_modules=['compat'],
      ext_modules=[Extension('_compat',
                             ['compat.cpp',
                              'compat_wrap.cpp'],
                             depends=['compat.hpp'],
                             extra_compile_args=['--std=gnu++11'])],
      scripts=['compat'],
     )
# conda install conda-forge/label/h6f9ffa1_0::gcc
# pip3 install -e .