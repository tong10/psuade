#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "CobylaOptimizer.h"
#include "Globals.h"

void evalF(int nInps, double *X, int nOuts, double *Y)
{
   Y[0] = X[0]*X[0] + X[1]*X[1] + X[0]*X[1] - 14*X[0] - 16*X[1] +
          pow(X[2]-10,2) + 4 * pow(X[3]-5,2) + pow(X[4]-3,2) +
          2*pow(X[5]-1,2) + 5*X[6]*X[6] + 7 * pow(X[7]-11,2) +
          2*pow(X[8]-10,2) + pow(X[9]-7,2) + 45;
   Y[1] = 105 - 4*X[0] - 5*X[1] + 3*X[6] - 9*X[7];
   Y[2] = - 10*X[0] + 8*X[1] + 17*X[6] - 2*X[7];
   Y[3] = 8*X[0] - 2*X[1] - 5*X[8] + 2*X[9] + 12;
   Y[4] = -3*pow(X[0]-2,2) - 4*pow(X[1]-3,2) - 2*pow(X[2],2) +
          7*X[3] + 120;
   Y[5] = -5*X[0]*X[0] - 8*X[1] - pow(X[2]-6,2) + 2*X[3] + 40;
   Y[6] = -0.5*pow(X[0]-8,2) - 2*pow(X[1]-4,2) - 3*X[4]*X[4] + X[5] + 30;
   Y[7] = -X[0]*X[0] - 2*pow(X[1]-2,2) + 2*X[0]*X[1] - 14*X[4] + 6*X[5];
   Y[8] = 3*X[0] - 6*X[1] - 12*pow(X[8]-8,2) + 7*X[9];
}

int main(int argc, char **argv)
{
   int    maxF=10000, nInps=10, ii, nOuts=9;
   double tol=1e-6, *X, *L, *U, *optX;
   CobylaOptimizer *opt;

   X = new double[nInps];
   L = new double[nInps];
   U = new double[nInps];
   for (ii = 0; ii < nInps; ii++) X[ii] = L[ii] = 0;
   for (ii = 0; ii < nInps; ii++) U[ii] = 10;
   for (ii = 0; ii < nInps; ii++)
     printf("Cobyla initial X %d = %12.4e\n",ii+1,X[ii]);
   opt = new CobylaOptimizer();
   opt->setObjectiveFunction(evalF);
   psConfig_.InteractiveOn();
   opt->optimize(nInps, X, L, U, nOuts, maxF, tol);
   optX = opt->getOptimalX();
   for (ii = 0; ii < nInps; ii++)
      printf("Cobyla final   X %d = %12.4e\n",ii+1,optX[ii]);
   printf("Optimal Y = %e\n", opt->getOptimalY());
   delete opt;
   delete [] X;
   delete [] L;
   delete [] U;
}


