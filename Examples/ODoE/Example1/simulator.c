#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
   int    ii, nInputs;
   double X[2], Y, pi=3.1416;
   FILE   *fIn  = fopen(argv[1], "r");
   FILE   *fOut;

   if (fIn == NULL)
   {
      printf("Simulator ERROR - cannot open in/out files.\n");
      exit(1);
   }
   fscanf(fIn, "%d", &nInputs);
   for (ii = 0; ii < nInputs; ii++) fscanf(fIn, "%lg", &X[ii]);
   Y = pow(X[0] - 0.5, 2.0) + pow(X[1] - 0.5, 2.0) + 0.1;
   fOut = fopen(argv[2], "w");
   fprintf(fOut, " %24.16e\n", Y);
   fclose(fOut);   
}

