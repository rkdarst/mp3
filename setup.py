#!/usr/bin/env python
from distutils.core import setup
import sys

py_modules = [ ]
#py_modules.append(["src/bin/"])
#py_modules.append(["mp3.system"])
#py_modules.append([""])
#py_modules.append([""])
#py_modules.append([""])
#py_modules.append([""])
#py_modules.append([""])

data_files = []
#data_files.append(('bin', ['src/bin/catdcd.py', 'src/bin/verify.py']))
#ext_modules=[Extension("lj.lj_c", ["src/lj/lj_c.c"]
                                              #include_dirs=C_include_dirs,
                                              #library_dirs=C_library_dirs,
                                              #libraries=["m", "argtable2", "gsl", "gslcblas"],
                                              #extra_compile_args=[]
#                                              )]
scripts = []
scripts.append('src/bin/mp3_verify')
scripts.append('src/bin/mp3_catdcd')
scripts.append('src/bin/mp3_geom')
scripts.append('src/bin/mp3_namdawk')
scripts.append('src/bin/mp3_tinkerawk')
scripts.append('src/bin/mp3_xstboxsize')
scripts.append('src/bin/mp3_xyztodcd')
scripts.append('src/bin/mp3_interactivesystem')


setup(name="mp3",
      version="0.4",
      description="Molecular Processing / data extraction tools",
      author="Richard Darst",
      author_email="rkd@zgib.net",
      url="http://www.cm.utexas.edu/rossky/",
      packages=['mp3'],
      py_modules=py_modules,
#      ext_modules=ext_modules,
      data_files=data_files,
      scripts=scripts,
      package_dir={'': 'src', 'bin':'src/bin'},  #where the root of the packages are
      )
