#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "LBFGSOptimizer.h"
#include "Globals.h"

void evalF(int nInps, double *X, int nOuts, double *Y)
{
   int ii;
   Y[0] = 0.0;
   for (ii = 0; ii < nInps-1; ii++)
   {
      Y[0] += (pow(1.0 - X[ii],2.0) + 100.0 * pow(X[ii+1] - X[ii]*X[ii],2.0));
      Y[ii+1] = -2.0*(1.0-X[ii]) + 200.0*(2.0*pow(X[ii],3.0)-2.0*X[ii]*X[ii+1]);
      if (ii > 0) Y[ii+1] += 200.0 * (X[ii] - X[ii-1] * X[ii-1]);
   }
   Y[nInps] = 200.0 * (X[nInps-1] - X[nInps-2] * X[nInps-2]);
}

int main(int argc, char **argv)
{
   int    maxF=10000, nInps=4, ii;
   double tol=1e-6, *X, *L, *U, *optX;
   LBFGSOptimizer *opt;

   X = new double[nInps];
   L = new double[nInps];
   U = new double[nInps];
   for (ii = 0; ii < nInps; ii++)
   {
      X[ii] = 0;
      L[ii] = -2;
      U[ii] = +2;
   }
   for (ii = 0; ii < nInps; ii++)
     printf("LBFGS initial X %d = %12.4e\n",ii+1,X[ii]);
   opt = new LBFGSOptimizer();
   opt->setObjectiveFunction(evalF);
   psConfig_.InteractiveOn();
   opt->optimize(nInps, X, L, U, maxF, tol);
   optX = opt->getOptimalX();
   for (ii = 0; ii < nInps; ii++)
      printf("LBFGS final   X %d = %12.4e\n",ii+1,optX[ii]);
   printf("Optimal Y = %e\n", opt->getOptimalY());
   delete opt;
   delete [] X;
   delete [] L;
   delete [] U;
}


