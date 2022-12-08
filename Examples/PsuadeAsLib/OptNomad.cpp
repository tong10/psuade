#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "NomadOptimizer.h"
#include "Globals.h"

int counter=0;
void evalF(int nInps, double *X, int nOuts, double *Yout)
{
   Yout[0] = X[0]*X[0] + X[1]*X[1] + 2.0*X[2]*X[2] + X[3]*X[3] -
             5.0*X[0] - 5.0*X[1] - 21.0*X[2] + 7.0*X[3];
   Yout[1] = -(-X[0]*X[0] - X[1]*X[1] - X[2]*X[2] - X[3]*X[3] -
             X[0] + X[1] - X[2] + X[3] + 8.0);
   Yout[2] = -(-X[0]*X[0] - 2.0*X[1]*X[1] - X[2]*X[2] - 2.0*X[3]*X[3] +
             X[0] + X[3] + 10.0);
   Yout[3] = -(-2.0*X[0]*X[0] - X[1]*X[1] - X[2]*X[2] -
             2.0*X[0] + X[1] + X[3] + 5.0);
   counter++;
   printf("Iteration %3d: Y = %e\n", counter, Yout[0]);
}

int main(int argc, char **argv)
{
   int    maxF=10000, nInps=4, ii, nOuts=4;
   double tol=1e-6, *lbnds, *ubnds, *X, *optX;
   NomadOptimizer *opt;

   lbnds = new double[nInps];
   ubnds = new double[nInps];
   X     = new double[nInps];
   for (ii = 0; ii < nInps; ii++) 
   {
      X[ii] = 0;
      lbnds[ii] = 0;
      ubnds[ii] = 4;
   }
   for (ii = 0; ii < nInps; ii++)
     printf("Nomad initial X %d = %12.4e\n",ii+1,X[ii]);
   opt = new NomadOptimizer();
   opt->setObjectiveFunction(evalF);
   opt->setDiscreteVariable(1);
   opt->setDiscreteVariable(2);
   opt->setDiscreteVariable(4);
   psConfig_.InteractiveOn();
   opt->optimize(nInps, X, lbnds, ubnds, nOuts, maxF, tol);
   optX = opt->getOptimalX();
   for (ii = 0; ii < nInps; ii++)
     printf("Nomad final   X %d = %12.4e\n",ii+1,optX[ii]);
   printf("Optimal Y = %e\n", opt->getOptimalY());
   delete opt;
   delete [] X;
   delete [] lbnds;
   delete [] ubnds;
}

