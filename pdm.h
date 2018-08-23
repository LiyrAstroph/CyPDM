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
  int i;
  double value;
}TypeSorter;

typedef struct
{
  double *x, *y;
  int n;
}TypeData;
TypeData data;

typedef struct
{
  double minVal, maxVal, dVal;
  int nVal, mode;
}TypeScanner;
TypeScanner scanner;

void pdmInit(double *jd, double *fs, int nd, double minVal, double maxVal, double dVal, int mode);
void setData(double *jd, double *fs, int nd);
void setScanner(double minVal, double maxVal, double dVal, int mode);
void pdmEnd();
void dophase(double period, double *phase, int *order);
int cmp(const void *a, const void *b);
int cmp_sorter(const void *a, const void *b);
void argsort(const double *x, int *order, int n);
void setUpEquiBlocks(int nbins, double *phase, int n, double *bbeg, double *bend, int *nb);
void setUpEquiBlocksCover(int nbins, int covers, double *bbeg, double *bend);
void pdmEquiBin(const int nbins, double *periods, double *thetas);
void pdmEquiBinCover(const int nbins, const int covers, double *periods, double *thetas);
double getTheta(double *phase, double *y, int n, double *bbeg, double *bend, int nb);
#endif