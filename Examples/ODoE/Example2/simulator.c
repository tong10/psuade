#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
   int    ii, nInputs;
   double X[3], Y;
   FILE   *fIn  = fopen(argv[1], "r");
   FILE   *fOut;

   if (fIn == NULL)
   {
      printf("Simulator ERROR - cannot open in/out files.\n");
      exit(1);
   }
   fscanf(fIn, "%d", &nInputs);
   for (ii = 0; ii < nInputs; ii++) fscanf(fIn, "%lg", &X[ii]);
   Y = 1.0 + X[0] * X[0] + X[1] * X[1] + 0.5 * 
       exp(-10*(X[0]-0.8)*(X[0]-0.8)) * 
       exp(-10*(X[1]-0.8)*(X[1]-0.8)) * (X[2] + 0.2); 
   fOut = fopen(argv[2], "w");
   fprintf(fOut, " %24.16e\n", Y);
   //Y = 1.0 + (1 - X[0]) * (1 - X[0]) + (1 - X[1]) * (1 - X[1]) + 0.5 * 
   //    exp(-2*(X[0]-0.2)*(X[0]-0.2)) * 
   //    exp(-2*(X[1]-0.2)*(X[1]-0.2)) * (X[2] + 0.2); 
   //fprintf(fOut, " %24.16e\n", Y);
   fclose(fOut);   
}

