import os
import numpy as np

from setuptools import setup
from setuptools.extension import Extension

include_dirs = [np.get_include()]
library_dirs = []
if os.name == 'nt':  # Windows, assumming MSVC compiler
  libraries = []
  compiler_args = ['/Ox', '/fp:fast']
elif os.name == 'posix':  # UNIX, assumming GCC compiler
  libraries = ['m']
  compiler_args = ['-O3', '-ffast-math']


"""
Allow users to install the module even if they do not have cython.
If cython is not found the c sources are compiled instead. More details at:
http://docs.cython.org/en/latest/src/reference/compilation.html
"""
try:
    from Cython.Build import cythonize
    USE_CYTHON = True
except ImportError:
    import warnings
    USE_CYTHON = False
    warnings.warn('Cython not found, compiling from c sources')


def no_cythonize(extensions, **_ignore):
    for extension in extensions:
        sources = []
        for sfile in extension.sources:
            path, ext = os.path.splitext(sfile)
            if ext in ('.pyx', '.py'):
                if extension.language == 'c++':
                    ext = '.cpp'
                else:
                    ext = '.c'
                sfile = path + ext
            sources.append(sfile)
        extension.sources[:] = sources
    return extensions

extensions = [
Extension("cypdm", 
	  sources=["cypdm.pyx", "cpdm.c"],
	  extra_compile_args=compiler_args,
    include_dirs=include_dirs,
    libraries=libraries,
    library_dirs=library_dirs
    ),
  ]

if USE_CYTHON:
    extensions = cythonize(extensions, annotate=False)
else:
    extensions = no_cythonize(extensions)

setup(
	name="cypdm",
	packages="cypdm",
	ext_modules = extensions,
  description = 'A fast package to apply the phase disperion minimization (PDM) algorithm',
  author = 'Yan-Rong Li',
  author_email = 'liyanrong@mail.ihep.ac.cn',
	)
