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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>

#include "cpdm.h"

/*!
 * create a TypePDM instance.  
 *
 */
TypePDM * cmkPDM(unsigned int nbins, unsigned int covers)
{
  TypePDM *pdm = (TypePDM *)malloc(sizeof(TypePDM));

  pdm->nbins = nbins;
  pdm->covers = covers; 

  return pdm;
}

/*!
 * free the TypePDM instance.
 */
void cfreePDM(TypePDM *pdm)
{
  free(pdm);
  return;
}


void cpdm(TypePDM *pdm, double *datax, double *datay, unsigned int nd, double *periods, double *thetas, unsigned int np)
{
  if(pdm->covers == 0)
  {
    cpdmEquiBin(pdm, datax, datay, nd, periods, thetas, np);
  }
  else
  {
    cpdmEquiBinCover(pdm, datax, datay, nd, periods, thetas, np);
  }
}

double cgetTheta(double *phase, double *y, unsigned int n, double *bbeg, double *bend, unsigned int nb)
{
  unsigned int i, j, ic;
  double mean, sigmaSqr, sSqrUp, sSqrDown, sigmaj;

  mean = 0.0;
  for(i=0; i<n; i++)
  {
    mean += y[i];
  }
  mean /= n;

  sigmaSqr = 0.0;
  for(i=0; i<n; i++)
  {
    sigmaSqr += (y[i] - mean) * (y[i] - mean);
  }
  sigmaSqr /= (n-1.0);
  
  sSqrUp = 0.0;
  sSqrDown = 0.0;
  for(i=0; i<nb; i++)
  {
    ic = 0;
    mean = 0.0;
    for(j=0; j<n; j++)
    {
      if(phase[j] >= bbeg[i] && phase[j] < bend[i])
      {
        mean += y[j];
        ic++;
      }
    }

    if(ic > 1)
    {
      mean /= ic;
      sigmaj = 0.0;
      for(j=0; j<n; j++)
      {
        if(phase[j] >= bbeg[i] && phase[j] < bend[i])
        {
          sigmaj += (y[j] - mean)*(y[j] - mean);
        }
      }
      sigmaj /= (ic-1.0);
      sSqrUp += (ic-1.0)*sigmaj;
      sSqrDown += ic;
    }
  }

  sSqrDown -= nb;

  return sSqrUp/sSqrDown / sigmaSqr;
}

void cpdmEquiBinCover(TypePDM *pdm, double *datax, double *datay, unsigned int nd, double *periods, double *thetas, unsigned int np)
{
  unsigned int i, j;
  double *phase, *phaseSort, *tmpy;
  double *bbeg, *bend;
  unsigned int *order;

  phase = malloc(nd*sizeof(double));
  phaseSort = malloc(2*nd*sizeof(double));
  tmpy = malloc(2*nd*sizeof(double));
  order = malloc(nd*sizeof(unsigned int));


  bbeg = malloc(pdm->nbins*pdm->covers*sizeof(double));
  bend = malloc(pdm->nbins*pdm->covers*sizeof(double));

  for(i=0; i<np; i++)
  {
    cdophase(datax, nd, periods[i], phase, order);

    for(j=0; j<nd; j++)
    {
      phaseSort[j] = phase[order[j]];
      phaseSort[j+nd] = phaseSort[j] + 1.0;
      tmpy[j] = datay[order[j]];
      tmpy[j+nd] = tmpy[j];
    }
        
    csetUpEquiBlocksCover(pdm->nbins, pdm->covers, bbeg,  bend);
    thetas[i]=cgetTheta(phaseSort, tmpy, 2*nd, bbeg, bend, pdm->nbins*pdm->covers);
  }

  free(phase);
  free(tmpy);
  free(phaseSort);
  free(order);
  free(bbeg);
  free(bend);
}

void cpdmEquiBin(TypePDM *pdm, double *datax, double *datay, unsigned int nd, double *periods, double *thetas, unsigned int np)
{
  unsigned int i, j;
  double *phase, *phaseSort;
  double *bbeg, *bend;
  unsigned int *order, nb;

  phase = malloc(nd*sizeof(double));
  phaseSort = malloc(nd*sizeof(double));
  order = malloc(nd*sizeof(unsigned int));

  bbeg = malloc(pdm->nbins*sizeof(double));
  bend = malloc(pdm->nbins*sizeof(double));

  for(i=0; i<np; i++)
  {
    cdophase(datax, nd, periods[i], phase, order);
    
    for(j=0; j<nd;j++)
    {
      phaseSort[j] = phase[order[j]];
    }
    
    csetUpEquiBlocks(pdm->nbins, phaseSort, nd, bbeg,  bend, &nb);
    thetas[i]=cgetTheta(phase, datay, nd, bbeg, bend, nb);
  }
  
  free(phase);
  free(phaseSort);
  free(order);
  free(bbeg);
  free(bend);
}


void csetUpEquiBlocks(unsigned int nbins, double *phaseSort, unsigned int n, double *bbeg, double *bend, unsigned int *nb)
{
  unsigned int i, j, ic;
  unsigned int *Ns;
  int nBlock, iPlus, iMinu, NPlus, NMinu;
  double *blockBegin, *blockEnd;
  int badBlock;
  blockBegin = bbeg;
  blockEnd = bend;

  nBlock = nbins;
  for(i=0; i<nbins; i++)
  {
    blockBegin[i] = i*1.0/nbins;
    blockEnd[i] = (i+1)*1.0/nbins;
  }
  
  Ns = malloc(nbins * sizeof(unsigned int));

  for(i=0; i<nBlock; i++)
  {
    ic=0;
    for(j=0; j<n; j++)
    {
      if(phaseSort[j] >= blockBegin[i] && phaseSort[j] < blockEnd[i])
      {
        ic++;
      }
    }
    Ns[i] = ic;
    //printf("%d %d\n", i, Ns[i]);
  }

  i=0;
  while(i<nBlock)
  {
    badBlock = -1;
    if(Ns[i] < minBinPoints)
    {
      iPlus = i + 1;
      iMinu = i - 1;
      NPlus=0;
      NMinu=0;

      if(iPlus == nBlock)
        NPlus =  INT_MAX;
      else
        NPlus = Ns[iPlus];

      if(iMinu == -1)
        NMinu = INT_MAX;
      else
        NMinu = Ns[iMinu];

      if(NMinu <= NPlus)
      {
        badBlock = i;
      }
      else
      {
        badBlock = i+1;
      }
    }

    if(badBlock != -1)
    {
      nBlock--;

      Ns[badBlock] += Ns[badBlock-1];

      for(i=badBlock; i<nBlock; i++)
      {
        blockBegin[i] = blockBegin[i+1];
      }
      for(i=badBlock-1; i<nBlock; i++)
      {
        blockEnd[i] = blockEnd[i+1];
        Ns[i] = Ns[i+1];
      }

      for(i=nBlock; i<nbins; i++)
      {
        blockBegin[i]=-1;
        blockEnd[i]=-1;
      }

      i=0; 
      continue;
    }
    i++;
  }

  *nb = nBlock;
  free(Ns);
  return;
}

void csetUpEquiBlocksCover(unsigned int nbins, unsigned int covers, double *bbeg, double *bend)
{
  unsigned int i, j;
  double offset;

  for(i=0; i<nbins; i++)
  {
    bbeg[i]  = i*1.0/nbins;
    bend[i] =  (i+1)*1.0/nbins;
  }

  for(i=1; i<covers; i++)
  {
    offset = i*1.0/(nbins*covers);
    for(j=0; j<nbins; j++)
    {
      bbeg[i*nbins + j] = bbeg[j] + offset;
      bend[i*nbins + j] = bend[j] + offset;
    }
  }

  return;
}

/*!
 *  folding phase of a light curve with a given period. 
 */
void cdophase(double *x, unsigned int n, double period, double *phase, unsigned int *order)
{
  unsigned int i;

  for(i=0; i<n; i++)
  {
    phase[i] = x[i]/period - floor(x[i]/period);
  }

  cargsort(phase, order, n);

  return;
}

/*!
 * comparison function
 */
int ccmp_sorter(const void *a, const void *b)
{
  return ((TypeSorter *)a)->value>=((TypeSorter *)b)->value?1:0;
}

/*!
 *  generate indices that would sort an array.
 */
void cargsort(const double *x, unsigned int *order, unsigned int n)
{
  unsigned int i;
  TypeSorter *sort = malloc(n*sizeof(TypeSorter));
  for(i=0; i<n; i++)
  {
    sort[i].i=i;
    sort[i].value=x[i];
  }
  qsort(sort, n, sizeof(TypeSorter), ccmp_sorter);

  for(i=0; i<n; i++)
    order[i] = sort[i].i;

  free(sort);
  return;
}