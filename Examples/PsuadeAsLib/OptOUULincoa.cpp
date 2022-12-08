#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "OUUOptimizer.h"
#include "Globals.h"

void evalF(int nInps, double *X, int nOuts, double *Y)
{
   Y[0] = 2 * X[0] * X[0] + X[1] * X[1];
}

int main(int argc, char **argv)
{
   int    maxF=10000, nInps=2, ii, nOuts=1;
   double tol=1e-6, *X, *L, *U, *optX;
   char   cfile[1000];
   OUUOptimizer *opt;

   X = new double[nInps];
   L = new double[nInps];
   U = new double[nInps];
   for (ii = 0; ii < nInps; ii++)
   {
      X[ii] = -0.5;
      L[ii] = -2;
      U[ii] = 2;
   }
   for (ii = 0; ii < nInps; ii++)
     printf("OUU Lincoa initial X %d = %12.4e\n",ii+1,X[ii]);
   opt = new OUUOptimizer();
   for (ii = 0; ii < nInps; ii++)
      opt->setVariableType(ii+1, OUUDesignCont);
   opt->setObjectiveFunction(evalF);
   strcpy(cfile, "psuade_lincoa_constraints");
   opt->setConstraintFile(cfile);
   psConfig_.InteractiveOn();
   opt->optimize(nInps, X, L, U, nOuts, maxF, tol);
   optX = opt->getOptimalX();
   for (ii = 0; ii < nInps; ii++)
     printf("OUU Lincoa final   X %d = %12.4e\n",ii+1,optX[ii]);
   printf("Optimal Y = %e\n", opt->getOptimalY());
   delete opt;
   delete [] X;
   delete [] L;
   delete [] U;
}


