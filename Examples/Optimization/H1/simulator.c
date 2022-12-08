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

   Y = pow(sin(X[0] - X[1]/8), 2.0);
   Y += pow(sin(X[1] + X[0]/8), 2.0);
   mult = pow(X[0] - 8.6998, 2.0) + pow(X[1] - 6.7665, 2.0);
   Y = Y / sqrt(mult+1);
   Y = - Y;
   free(X);
   fOut = fopen(argv[2], "w");
   fprintf(fOut, "%24.16e\n", Y);
   fclose(fOut);
   return 0;
}
