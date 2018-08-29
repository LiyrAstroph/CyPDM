#!/usr/bin/python
#cython: initializedcheck=False, boundscheck=False, wraparound=False, cdivision=True, profile=False

cimport cython
from cpython.mem cimport PyMem_Malloc, PyMem_Free, PyMem_Realloc

import numpy as np
cimport numpy as np

ctypedef double DTYPE_t


"""
 * CyPDM
 * 
 * A fast package to apply the phase disperion minimization (PDM) algorithm, 
 * based on PyAstronomy module PyPDM (https://github.com/sczesla/PyAstronomy).
 *
 * This is a Cython verison with improved computation speed.
 * 
 * The PDM alogrithm refers to the reference:
 *    http://adsabs.harvard.edu/abs/1978ApJ...224..953S
 *
 *
 * Author:
 *   Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * 
"""

cdef extern from "cpdm.h":
  ctypedef struct TypePDM

  TypePDM * cmkPDM(unsigned int nbins, unsigned int covers, unsigned int nd)
  
  void cgetScan(TypePDM * pdm, unsigned int *nbins, unsigned int *covers)

  void csetScan(TypePDM * pdm, unsigned int nbins, unsigned int covers)

  void cfreePDM(TypePDM * pdm)

  void cpdmEquiBin(TypePDM *pdm, DTYPE_t *datax, DTYPE_t *datay, unsigned int nd, \
    DTYPE_t *periods, DTYPE_t *thetas, unsigned int nper)

  void cpdmEquiBinCover(TypePDM *pdm, DTYPE_t *datax, DTYPE_t *datay, unsigned int nd, \
    DTYPE_t *periods, DTYPE_t *thetas, unsigned int nper)

  void cpdm(TypePDM *pdm, DTYPE_t *datax, DTYPE_t *datay, unsigned int nd, \
    DTYPE_t *periods, DTYPE_t *thetas, unsigned int nper)


cdef class CyPDM:
  cdef TypePDM * _thisptr
  cdef DTYPE_t *jd
  cdef DTYPE_t *fs
  cdef unsigned int n

  def __cinit__(self, DTYPE_t [:] jd, DTYPE_t [:] fs, int nbins=10, int covers=3):

    if nbins <= 0:
      msg = "nbins=%d incorrect."%nbins
      raise ValueError(msg)
    if covers <= 0:
      msg = "covers=%d incorrect."%covers
      raise ValueError(msg)

    self.n = len(jd)
    self.jd = <DTYPE_t *>PyMem_Malloc(self.n * sizeof(DTYPE_t))
    self.fs = <DTYPE_t *>PyMem_Malloc(self.n * sizeof(DTYPE_t))

    for i in range(self.n):
      self.jd[i] = jd[i]
      self.fs[i] = fs[i]   

    self._thisptr = cmkPDM(nbins, covers, self.n)

    if self._thisptr == NULL or self.jd == NULL or self.fs == NULL:
      msg = "Fail to create PDM instance."
      raise MemoryError(msg)


  def __dealloc__(self):
    if self._thisptr != NULL:
      cfreePDM(self._thisptr)

    if self.jd != NULL:
      PyMem_Free(self.jd)

    if self.fs != NULL:
      PyMem_Free(self.fs)

  cpdef setData(self, DTYPE_t [:] jd, DTYPE_t [:] fs):
    if len(jd) != self.n:
      print "realloc self.jd and self.fs."
      self.n = len(jd)
      PyMem_Realloc(self.jd, len(jd) * sizeof(double))
      PyMem_Realloc(self.fs, len(jd) * sizeof(double))
      if self.jd == NULL or self.fs == NULL:
        msg = "Fail to realloc self.jd or self.fs."
        raise MemoryError(msg)

    for i in range(self.n):
      self.jd[i] = jd[i]
      self.fs[i] = fs[i]

  cpdef getScan(self):
    cdef unsigned int nbins
    cdef unsigned int covers
    cgetScan(self._thisptr, &nbins, &covers)

    return nbins, covers

  cpdef setScan(self, unsigned int nbins, unsigned int covers):
    if nbins <= 0:
      msg = "nbins=%d incorrect."%nbins
      raise ValueError(msg)
    if covers <= 0:
      msg = "covers=%d incorrect."%covers
      raise ValueError(msg)

    csetScan(self._thisptr, nbins, covers)

    return

  cpdef getPDM(self, np.ndarray[DTYPE_t, ndim=1] periods):
    thetas = np.empty_like(periods)
    cpdm(self._thisptr, self.jd, self.fs, self.n,  \
      <DTYPE_t *>np.PyArray_DATA(periods), <DTYPE_t *>np.PyArray_DATA(thetas), len(periods))
    return thetas

  cpdef getPDM_EquiBin(self, np.ndarray[DTYPE_t, ndim=1] periods):
    thetas = np.empty_like(periods)
    cpdmEquiBin(self._thisptr, self.jd, self.fs, self.n, \
      <DTYPE_t *>np.PyArray_DATA(periods), <DTYPE_t *>np.PyArray_DATA(thetas), len(periods))
    return thetas

  cpdef getPDM_EquiBinCover(self, np.ndarray[DTYPE_t, ndim=1] periods):
    thetas = np.empty_like(periods)
    cpdmEquiBinCover(self._thisptr, self.jd, self.fs, self.n, \
      <DTYPE_t *>np.PyArray_DATA(periods), <DTYPE_t *>np.PyArray_DATA(thetas), len(periods))
    return thetas