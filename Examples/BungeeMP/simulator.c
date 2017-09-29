#include <math.h>
#include <stdio.h>
#include <stdlib.h>

main(int argc, char **argv)
{
   int    i, nInputs, count;
   double X[3], Y, kel=1.5, g=9.8;
   FILE   *fIn  = fopen("psuadeApps.in", "r");
   FILE   *fOut;

   if (fIn == NULL)
   {
      printf("Simulator ERROR - cannot open in/out files.\n");
      exit(1);
   }
   fscanf(fIn, "%d", &nInputs);
   for (i = 0; i < nInputs; i++) fscanf(fIn, "%lg", &X[i]);
   Y = X[0] - 2.0 * X[1] * g / (kel * X[2]);
   system("sleep 60");
   fOut = fopen("psuadeApps.out", "w");
   fprintf(fOut, " %24.16e\n", Y);
   fclose(fOut);   
}

