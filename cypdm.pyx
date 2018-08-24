#!/usr/bin/python
#cython: initializedcheck=False, boundscheck=False, wraparound=False, cdivision=False, profile=False

cimport cython
from cpython.mem cimport PyMem_Malloc, PyMem_Free

import numpy as np
cimport numpy as np

ctypedef double DTYPE_t

cdef extern from "cpdm.h":
  ctypedef struct TypePDM
  TypePDM * cmkPDM(double *jd, double *fs, unsigned int nd, unsigned int nbins, unsigned int covers)
  void cfreePDM(TypePDM * pdm)
  void cpdmEquiBin(TypePDM *pdm, double *periods, double *thetas, unsigned int nper)
  void cpdmEquiBinCover(TypePDM *pdm, double *periods, double *thetas, unsigned int nper)
  void cpdm(TypePDM *pdm, double *periods, double *thetas, unsigned int nper)


cdef class PyPDM:
  cdef TypePDM * _thisptr
  cdef double *jd
  cdef double *fs
  cdef unsigned int n
  cdef unsigned int nbins
  cdef unsigned int covers

  def __cinit__(self, DTYPE_t [:] jd, DTYPE_t [:] fs, int nbins=10, int covers=3):

    if nbins <= 0:
      msg = "nbins=%d incorrect."%nbins
      raise ValueError(msg)
    if covers <= 0:
      msg = "covers=%d incorrect."%covers
      raise ValueError(msg)

    self.nbins = nbins
    self.covers = covers
    self.n = len(jd)

    self.jd = <double *>PyMem_Malloc(self.n * sizeof(DTYPE_t))
    self.fs = <double *>PyMem_Malloc(self.n * sizeof(DTYPE_t))

    for i in range(self.n):
      self.jd[i] = jd[i]
      self.fs[i] = fs[i]   

    self._thisptr = cmkPDM(self.jd, self.fs, self.n, self.nbins, self.covers)

    if self._thisptr == NULL or self.jd == NULL or self.fs == NULL:
      msg = "Fail to creat PDM instance."
      raise MemoryError(msg)


  def __dealloc__(self):
    if self._thisptr != NULL:
      cfreePDM(self._thisptr)

    if self.jd != NULL:
      PyMem_Free(self.jd)

    if self.fs != NULL:
      PyMem_Free(self.fs)

  cpdef setData(self, DTYPE_t [:] jd, DTYPE_t [:] fs):
    for i in range(self.n):
      self.jd[i] = jd[i]
      self.fs[i] = fs[i]

  cpdef getPDM(self, np.ndarray[DTYPE_t, ndim=1] periods):
    thetas = np.empty_like(periods)
    cpdm(self._thisptr, <double *>np.PyArray_DATA(periods), <double *>np.PyArray_DATA(thetas), len(periods))
    return thetas

  cpdef getPDM_EquiBin(self, np.ndarray[DTYPE_t, ndim=1] periods):
    thetas = np.empty_like(periods)
    cpdmEquiBin(self._thisptr, <double *>np.PyArray_DATA(periods), <double *>np.PyArray_DATA(thetas), len(periods))
    return thetas

  cpdef getPDM_EquiBinCover(self, np.ndarray[DTYPE_t, ndim=1] periods):
    thetas = np.empty_like(periods)
    cpdmEquiBinCover(self._thisptr, <double *>np.PyArray_DATA(periods), <double *>np.PyArray_DATA(thetas), len(periods))
    return thetas