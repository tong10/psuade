#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "OUUOptimizer.h"
#include "Globals.h"

int counter=0;
void evalF(int nInps, double *X, int nOuts, double *Yout)
{
   int flag=0;
   Yout[0] = X[0]*X[0] + X[1]*X[1] + 2.0*X[2]*X[2] + X[3]*X[3] -
             5.0*X[0] - 5.0*X[1] - 21.0*X[2] + 7.0*X[3];
   Yout[1] = -(-X[0]*X[0] - X[1]*X[1] - X[2]*X[2] - X[3]*X[3] -
             X[0] + X[1] - X[2] + X[3] + 8.0);
   if (Yout[1] > 0) flag++;
   Yout[2] = -(-X[0]*X[0] - 2.0*X[1]*X[1] - X[2]*X[2] - 2.0*X[3]*X[3] +
             X[0] + X[3] + 10.0);
   if (Yout[2] > 0) flag++;
   Yout[3] = -(-2.0*X[0]*X[0] - X[1]*X[1] - X[2]*X[2] -
             2.0*X[0] + X[1] + X[3] + 5.0);
   if (Yout[3] > 0) flag++;
   counter++;
   printf("Iteration %3d: \n", counter);
   printf("   X = %14.6e %14.6e %14.6e %14.6e\n",X[0],X[1],X[2],X[3]);
   printf("   Y = %e ", Yout[0]);
   if (flag == 0) printf(" - feasible\n");
   else           printf(" - infeasible\n");
}

int main(int argc, char **argv)
{
   int    maxF=10000, nInps=4, ii, nOuts=4;
   double tol=1e-6, *lbnds, *ubnds, *X;
   OUUOptimizer *opt;

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
     printf("OUU Nomad initial X %d = %12.4e\n",ii+1,X[ii]);
   opt = new OUUOptimizer();
   for (ii = 0; ii < nInps; ii++)
      opt->setVariableType(ii+1, OUUDesignCont);
   opt->setObjectiveFunction(evalF);
   opt->setDiscreteVariable(1);
   opt->setDiscreteVariable(2);
   opt->setDiscreteVariable(4);
   psConfig_.InteractiveOn();
   opt->optimize(nInps, X, lbnds, ubnds, nOuts, maxF, tol);
   for (ii = 0; ii < nInps; ii++)
     printf("OUU Nomad final   X %d = %12.4e\n",ii+1,X[ii]);
   printf("Optimal Y = %e\n", opt->getOptimalY());
   delete opt;
   delete [] X;
   delete [] lbnds;
   delete [] ubnds;
}

