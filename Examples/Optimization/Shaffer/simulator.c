#include <math.h>
#include <stdio.h>
#include <stdlib.h>

main(int argc, char **argv)
{
   int    n, ii;
   double *X, Y, mult;
   char   lineIn[100], stringPtr[100], equal[2];

   FILE  *fIn  = fopen(argv[1], "r");
   FILE  *fOut;
   if (fIn == NULL)
   {
      printf("Griewank ERROR - cannot open in/out files.\n");
      exit(1);
   }
   fscanf(fIn, "%d",  &n);
   X = (double *) malloc(n * sizeof(double));
   for (ii = 0; ii < n; ii++) fscanf(fIn, "%lg", &X[ii]);
   fclose(fIn);

   Y = 0;
   for (ii = 0; ii < 2; ii++)
   {
      mult = X[ii] * X[ii] + X[ii+1] * X[ii+1];
      Y = Y + pow(mult,0.25) * (pow(sin(50*pow(mult,0.1)),2.0) + 1);
   }
   free(X);
   fOut = fopen(argv[2], "w");
   fprintf(fOut, "%24.16e\n", Y);
   fclose(fOut);
   return 0;
}
