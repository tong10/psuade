#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
   int    i, nInputs;
   double *X, Y;
   FILE   *fIn  = fopen(argv[1], "r");
   FILE   *fOut;

   if (fIn == NULL)
   {
      printf("Simulator ERROR - cannot open in/out files.\n");
      exit(1);
   }
   fscanf(fIn, "%d", &nInputs);
   X = (double *) malloc(nInputs * sizeof(double));
   for (i = 0; i < nInputs; i++) fscanf(fIn, "%lg", &X[i]);
   Y = 0.0;
   for (i = 0; i < nInputs; i++) Y += X[i] * X[i];
   fOut = fopen(argv[2], "w");
   fprintf(fOut, " %24.16e\n", Y);
   for (i = 0; i < nInputs; i++) fprintf(fOut, " %24.16e\n",2*X[i]);
   fclose(fOut);   
   free(X);
}

