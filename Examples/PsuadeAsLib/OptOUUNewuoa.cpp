#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "OUUOptimizer.h"
#include "Globals.h"

void evalF(int nInps, double *X, int nOuts, double *Yout)
{
   double X1, X2, X3, X4, Y;
   X1 = X[0];
   X2 = X[1];
   X3 = X[2];
   X4 = X[3];
   Y = 100.0 * (X2 - X1 * X1) * (X2 - X1 * X1)  + (1.0 - X1) * (1.0 - X1);
   Y += 100.0 * (X3 - X2 * X2) * (X3 - X2 * X2)  + (1.0 - X2) * (1.0 - X2);
   Y += 100.0 * (X4 - X3 * X3) * (X4 - X3 * X3)  + (1.0 - X3) * (1.0 - X3);
   Yout[0] = Y;
}

int main(int argc, char **argv)
{
   int    maxF=10000, nInps=4, ii, nOuts=1;
   double tol=1e-6, *lbnds, *ubnds, *X, *optX;
   OUUOptimizer *opt;

   lbnds = new double[nInps];
   ubnds = new double[nInps];
   X     = new double[nInps];
   for (ii = 0; ii < nInps; ii++) 
   {
      X[ii] = -1.0;
      lbnds[ii] = -2;
      ubnds[ii] =  2;
   }
   for (ii = 0; ii < nInps; ii++)
     printf("OUU Newuoa initial X %d = %12.4e\n",ii+1,X[ii]);
   opt = new OUUOptimizer();
   opt->setLocalOptimizer(1);
   opt->setObjectiveFunction(evalF);
   for (ii = 0; ii < nInps; ii++)
      opt->setVariableType(ii+1, OUUDesignCont);
   psConfig_.InteractiveOn();
   opt->optimize(nInps, X, lbnds, ubnds, nOuts, maxF, tol);
   optX = opt->getOptimalX();
   for (ii = 0; ii < nInps; ii++)
     printf("OUU Newuoa final   X %d = %12.4e\n",ii+1,optX[ii]);
   printf("Optimal Y = %e\n", opt->getOptimalY());
   delete opt;
   delete [] X;
   delete [] lbnds;
   delete [] ubnds;
}


