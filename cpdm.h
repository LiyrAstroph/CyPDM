/*
 * CyPDM
 * 
 * A fast package to apply the phase disperion minimization (PDM) algorithm, 
 * based on PyAstronomy module PyPDM (https://github.com/sczesla/PyAstronomy).
 *
 * This is a C verison with improved computation speed.
 * 
 * The PDM alogrithm refers to the reference:
 *    http://adsabs.harvard.edu/abs/1978ApJ...224..953S
 *
 *
 * Author:
 *   Yan-Rong Li, liyanrong@mail.ihep.ac.cn
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
  unsigned int nbins, covers;
}TypePDM;

void cdophase(double *x, unsigned int n, double period, double *phase, unsigned int *order);
int ccmp_sorter(const void *a, const void *b);
void cargsort(const double *x, unsigned int *order, unsigned int n);
void csetUpEquiBlocks(unsigned int nbins, double *phase, unsigned int n, double *bbeg, double *bend, unsigned int *nb);
void csetUpEquiBlocksCover(unsigned int nbins, unsigned int covers, double *bbeg, double *bend);
void cpdmEquiBin(TypePDM *pdm, double *datax, double *datay, unsigned int nd, double *periods, double *thetas, unsigned int np);
void cpdmEquiBinCover(TypePDM *pdm, double *datax, double *datay, unsigned int nd, double *periods, double *thetas, unsigned int np);
double cgetTheta(double *phase, double *y, unsigned int n, double *bbeg, double *bend, unsigned int nb);
void cpdm(TypePDM *pdm, double *datax, double *datay, unsigned int nd, double *periods, double *thetas, unsigned int np);

TypePDM * cmkPDM(unsigned int mode, unsigned int type);
void cfreePDM(TypePDM *pdm);
#endif