from distutils.core import setup, Extension
from Cython.Build import cythonize

ext = Extension("cypdm", sources=["cypdm.pyx", "cpdm.c"])

setup(
	name="cypdm",
	ext_modules = cythonize(ext)
	)