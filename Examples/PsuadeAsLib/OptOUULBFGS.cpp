#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "OUUOptimizer.h"
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
   int    maxF=10000, nInps=4, ii, nOuts=5;
   double tol=1e-6, *X, *L, *U, *optX;
   OUUOptimizer *opt;

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
     printf("OUU LBFGS initial X %d = %12.4e\n",ii+1,X[ii]);
   opt = new OUUOptimizer();
   for (ii = 0; ii < nInps; ii++)
      opt->setVariableType(ii+1, OUUDesignCont);
   opt->setObjectiveFunction(evalF);
   opt->setLocalOptimizer(3);
   psConfig_.InteractiveOn();
   opt->optimize(nInps, X, L, U, nOuts, maxF, tol);
   optX = opt->getOptimalX();
   for (ii = 0; ii < nInps; ii++)
      printf("OUU LBFGS final   X %d = %12.4e\n",ii+1,optX[ii]);
   printf("Optimal Y = %e\n", opt->getOptimalY());
   delete opt;
   delete [] X;
   delete [] L;
   delete [] U;
}


