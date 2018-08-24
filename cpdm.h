/*
 * CyPDM
 * 
 * A fast package to apply the phase disperion minimization (PDM) algorithm. 
 *
 * Based on PyAstronomy module PyPDM (https://github.com/sczesla/PyAstronomy).
 * 
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * 
 */

#ifndef _PDM_H
#include <stdlib.h>
#include <stdio.h>

#define minBinPoints (3)

typedef struct
{
  unsigned int i;
  double value;
}TypeSorter;

typedef struct 
{
  int n;
  double *x, *y;
  unsigned int nbins, covers;
}TypePDM;

void cdophase(double *x, unsigned int n, double period, double *phase, unsigned int *order);
int ccmp_sorter(const void *a, const void *b);
void cargsort(const double *x, unsigned int *order, unsigned int n);
void csetUpEquiBlocks(unsigned int nbins, double *phase, unsigned int n, double *bbeg, double *bend, unsigned int *nb);
void csetUpEquiBlocksCover(unsigned int nbins, unsigned int covers, double *bbeg, double *bend);
void cpdmEquiBin(TypePDM *pdm, double *periods, double *thetas, unsigned int np);
void cpdmEquiBinCover(TypePDM *pdm, double *periods, double *thetas, unsigned int np);
double cgetTheta(double *phase, double *y, unsigned int n, double *bbeg, double *bend, unsigned int nb);
void cpdm(TypePDM *pdm, double *periods, double *thetas, unsigned int np);

TypePDM * cmkPDM(double *jd, double *fs, unsigned int nd, unsigned int mode, unsigned int type);
void cfreePDM(TypePDM *pdm);
#endif