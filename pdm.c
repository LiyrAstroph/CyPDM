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

#include "pdm.h"

int main(int argc, char **argv)
{
  FILE *fp;
  int nd_max = 5000, nd;
  double jd[nd_max],mag[nd_max];
  int i, nbins, covers;
  char buf[256];
  double *periods, *thetas;

  nbins=10;
  covers=3;

  fp = fopen("test/con_all.txt", "r");
  if(fp==NULL)
  {
    printf("Cannot open con_all.txt.\n");
    exit(-1);
  }

  for(i=0; i<nd_max; i++)
  {
    fgets(buf, 256, fp);
    sscanf(buf, "%lf %lf", &jd[i], &mag[i]);
    if(feof(fp)!=0)
      break;
  }
  fclose(fp);
  nd = i;
  
  
  pdmInit(jd, mag, nd, 1.0/(50.0*365.0), 1.0e-2+1.0/(50.0*365.0), 1.0e-5, 0);

  periods = malloc(scanner.nVal * sizeof(double));
  thetas = malloc(scanner.nVal * sizeof(double));

  //pdmEquiBin(nbins, periods, thetas);

  pdmEquiBinCover(nbins, covers, periods, thetas);

  for(i=0; i<scanner.nVal; i++)
  {
    printf("%f %f\n", periods[i], thetas[i]);
  }

  pdmEnd();
}

void pdmInit(double *jd, double *fs, int nd, double minVal, double maxVal, double dVal, int mode)
{
  setData(jd, fs, nd);
  setScanner(minVal, maxVal, dVal, mode);
  return;
}

void setData(double *jd, double *fs, int nd)
{
  int i;
  data.n = nd;
  data.x = (double *)malloc(nd*sizeof(double));
  data.y = (double *)malloc(nd*sizeof(double));

  memcpy(data.x, jd, nd*sizeof(double));
  memcpy(data.y, fs, nd*sizeof(double));

  for(i=nd-1; i>=0; i--)
  {
    data.x[i] -= data.x[0];
  }
  
  return;
}

void setScanner(double minVal, double maxVal, double dVal, int mode)
{
  scanner.nVal = ceil((maxVal - minVal)/dVal + 1);
  scanner.minVal = minVal;
  scanner.maxVal = maxVal;
  scanner.dVal = dVal;
  scanner.mode = mode;
}

void pdmEnd()
{
  free(data.x);
  free(data.y);
}

double getTheta(double *phase, double *y, int n, double *bbeg, double *bend, int nb)
{
  int i, j, N, M, ic;
  double mean, sigmaSqr, sSqrUp, sSqrDown, sigmaj;

  mean = 0.0;
  for(i=0; i<data.n; i++)
  {
    mean += y[i];
  }
  mean /= data.n;

  N = n;
  M = nb;

  sigmaSqr = 0.0;
  for(i=0; i<N; i++)
  {
    sigmaSqr += (y[i] - mean) * (y[i] - mean);
  }
  sigmaSqr /= (N-1.0);
  
  sSqrUp = 0.0;
  sSqrDown = 0.0;
  for(i=0; i<M; i++)
  {
    ic = 0;
    mean = 0.0;
    for(j=0; j<N; j++)
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
      for(j=0; j<N; j++)
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

  sSqrDown -= M;

  return sSqrUp/sSqrDown / sigmaSqr;
}

void pdmEquiBinCover(const int nbins, const int covers, double *periods, double *thetas)
{
  int i, j;
  double *phase, *phaseSort, *tmpy;
  double *bbeg, *bend;
  int *order;

  phase = malloc(data.n*sizeof(double));
  phaseSort = malloc(2*data.n*sizeof(double));
  tmpy = malloc(2*data.n*sizeof(double));
  order = malloc(data.n*sizeof(int));

  if(scanner.mode == 0)
  {
    for(i=0; i<scanner.nVal; i++)
    {
      periods[i] = 1.0/(scanner.minVal + i*scanner.dVal);
    }
  }
  else
  {
    for(i=0; i<scanner.nVal; i++)
    {
      periods[i] = scanner.minVal + i*scanner.dVal;
    }
  }

  bbeg = malloc(nbins*covers*sizeof(double));
  bend = malloc(nbins*covers*sizeof(double));

  for(i=0; i<scanner.nVal; i++)
  {
    dophase(periods[i], phase, order);

    for(j=0; j<data.n; j++)
    {
      phaseSort[j] = phase[order[j]];
      phaseSort[j+data.n] = phaseSort[j] + 1.0;
      tmpy[j] = data.y[order[j]];
      tmpy[j+data.n] = tmpy[j];
    }
        
    setUpEquiBlocksCover(nbins, covers, bbeg,  bend);
    thetas[i]=getTheta(phaseSort, tmpy, 2*data.n, bbeg, bend, nbins*covers);
  }

  free(phase);
  free(tmpy);
  free(phaseSort);
  free(order);
  free(bbeg);
  free(bend);
}

void pdmEquiBin(const int nbins, double *periods, double *thetas)
{
  int i, j;
  double *phase, *phaseSort;
  double *bbeg, *bend;
  int *order, nb;

  phase = malloc(data.n*sizeof(double));
  phaseSort = malloc(data.n*sizeof(double));
  order = malloc(data.n*sizeof(int));

  if(scanner.mode == 0)
  {
    for(i=0; i<scanner.nVal; i++)
    {
      periods[i] = 1.0/(scanner.minVal + i*scanner.dVal);
    }
  }
  else
  {
    for(i=0; i<scanner.nVal; i++)
    {
      periods[i] = scanner.minVal + i*scanner.dVal;
    }
  }

  bbeg = malloc(nbins*sizeof(double));
  bend = malloc(nbins*sizeof(double));

  for(i=0; i<scanner.nVal; i++)
  {
    dophase(periods[i], phase, order);
    
    for(j=0; j<data.n;j++)
    {
      phaseSort[j] = phase[order[j]];
    }
    
    setUpEquiBlocks(nbins, phaseSort, data.n, bbeg,  bend, &nb);
    thetas[i]=getTheta(phase, data.y, data.n, bbeg, bend, nb);
  }
  
  free(phase);
  free(phaseSort);
  free(order);
  free(bbeg);
  free(bend);
}


void setUpEquiBlocks(int nbins, double *phaseSort, int n, double *bbeg, double *bend, int *nb)
{
  int i, j, ic;
  int *Ns;
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
  
  Ns = malloc(nbins * sizeof(int));

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

  /*printf("nb: %d\n", *nb);
  for(i=0; i<*nb; i++)
  {
    printf("%d %f %f %d\n", i, bbeg[i], bend[i], Ns[i]);
  }*/

  free(Ns);
}

void setUpEquiBlocksCover(int nbins, int covers, double *bbeg, double *bend)
{
  int i, j;
  double offset;

  for(i=0; i<nbins; i++)
  {
    bbeg[i]  = i*1.0/nbins;
    bend[i] =  (i+1)*1.0/nbins;
  }
  
  offset = 0.0;
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

void dophase(double period, double *phase, int *order)
{
  int i;

  for(i=0; i<data.n; i++)
  {
    phase[i] = data.x[i]/period - floor(data.x[i]/period);
  }

  argsort(phase, order, data.n);

  /*for(i=0; i<10; i++)
  {
    printf("%f %f %d\n", phase[i], phase[order[i]], order[i]);
  }*/
  return;
}

int cmp(const void *a, const void *b)
{
  return (*(double *)a)>=(*(double *)b)?1:0;
}

int cmp_sorter(const void *a, const void *b)
{
  return ((TypeSorter *)a)->value>=((TypeSorter *)b)->value?1:0;
}

void argsort(const double *x, int *order, int n)
{
  int i;
  TypeSorter *sort = malloc(n*sizeof(TypeSorter));
  for(i=0; i<n; i++)
  {
    sort[i].i=i;
    sort[i].value=x[i];
  }
  qsort(sort, n, sizeof(TypeSorter), cmp_sorter);

  for(i=0; i<n; i++)
    order[i] = sort[i].i;

  free(sort);
  return;
}