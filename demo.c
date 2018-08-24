#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>

#include "cpdm.h"

int nd_max = 5000;

int readData(double *jd, double *mag, int *nd)
{
  FILE *fp;
  int i;
  char buf[256];

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
  *nd = i;

  for(i=*nd-1; i>=0; i--)
  {
    jd[i] -= jd[0];
  }
  return 0;
}

int main(int argc, char **argv)
{
  int nd, np = 995;
  int i, nbins, covers;
  double jd[nd_max], mag[nd_max];
  double *periods, *thetas;
  TypePDM * pdm;

  periods = malloc(np * sizeof(double));
  thetas = malloc(np * sizeof(double));

  nbins=10;
  covers=0;

  for(i=0; i<np; i++)
  {
    periods[i] = (1.0/365.0/50.0) + i*1.0e-5;
    periods[i] = 1.0/periods[i];
  }

  readData(jd, mag, &nd);

  pdm = cmkPDM(nbins, covers);

  cpdm(pdm, jd, mag, nd, periods, thetas, np);

  for(i=0; i<np; i++)
  {
    printf("%f %f\n", periods[i], thetas[i]);
  }

  free(periods);
  free(thetas);
}