#include <math.h>
#include <stdio.h>
#include <stdlib.h>

main(int argc, char **argv)
{
   int    n, ii;
   double *X, Y;
   char   lineIn[100], stringPtr[100], equal[2];

   FILE  *fIn  = fopen(argv[1], "r");
   FILE  *fOut;
   if (fIn == NULL)
   {
      printf("Schwefel ERROR - cannot open in/out files.\n");
      exit(1);
   }
   fscanf(fIn, "%d",  &n);
   X = (double *) malloc(n * sizeof(double));
   for (ii = 0; ii < n; ii++) fscanf(fIn, "%lg", &X[ii]);
   fclose(fIn);

   Y = 418.9828872724339 * 3;
   for (ii = 0; ii < 3; ii++)
      Y = Y - X[ii] * sin(sqrt(fabs(X[ii])));
   free(X);
   fOut = fopen(argv[2], "w");
   fprintf(fOut, "%24.16e\n", Y);
   fclose(fOut);
   return 0;
}
