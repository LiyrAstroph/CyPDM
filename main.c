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

  TypePDM * pdm;

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

  pdm = mkPDM(jd, mag, nd, 1.0/(50.0*365.0), 1.0e-2+1.0/(50.0*365.0), 1.0e-5, 0);

  periods = malloc(pdm->scanner->nVal * sizeof(double));
  thetas = malloc(pdm->scanner->nVal * sizeof(double));

  //pdmEquiBin(pdm, nbins, periods, thetas);

  pdmEquiBinCover(pdm, nbins, covers, periods, thetas);

  for(i=0; i<pdm->scanner->nVal; i++)
  {
    printf("%f %f\n", periods[i], thetas[i]);
  }
}