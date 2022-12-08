#include <math.h>
#include <stdio.h>
#include <stdlib.h>

main(int argc, char **argv)
{
   int    i, nInputs, count;
   double X[2], Y;
   FILE   *fIn  = fopen(argv[1], "r");
   FILE   *fOut;

   if (fIn == NULL)
   {
      printf("Simulator ERROR - cannot open in/out files.\n");
      exit(1);
   }
   fscanf(fIn, "%d", &nInputs);
   for (i = 0; i < nInputs; i++) fscanf(fIn, "%lg", &X[i]);
/*
   if (X[0] < 0.5 && X[1] < 0.5) Y = 0;
   if (X[0] < 0.5 && X[1] >= 0.5) Y = X[1] * X[1];
   if (X[0] >= 0.5 && X[1] < 0.5) Y = X[0] * X[0];
   if (X[0] >= 0.5 && X[1] >= 0.5) Y = X[0] * X[0] + X[1] * X[1];
*/
   if ((X[0]+X[1]) < 1.3) Y = 0;
   else                   Y = 1.0 * X[0] * X[1];
   fOut = fopen(argv[2], "w");
   fprintf(fOut, " %24.16e\n", Y);
   fclose(fOut);   
}

