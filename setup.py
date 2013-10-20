# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 23:45:07 2013

@author: Wiliam
"""

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy

ext = Extension("momentum", ["momentum.pyx"],
    include_dirs = [numpy.get_include()])
                
setup(ext_modules=[ext],
      cmdclass = {'build_ext': build_ext})