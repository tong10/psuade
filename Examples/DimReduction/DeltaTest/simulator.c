#include <math.h>
#include <stdio.h>
#include <stdlib.h>

main(int argc, char **argv)
{
   int    i, nInputs, count;
   double X[6], Y, pi=3.1416;
   FILE   *fIn  = fopen(argv[1], "r");
   FILE   *fOut;

   if (fIn == NULL)
   {
      printf("Simulator ERROR - cannot open in/out files.\n");
      exit(1);
   }
   fscanf(fIn, "%d", &nInputs);
   for (i = 0; i < nInputs; i++) fscanf(fIn, "%lg", &X[i]);
   Y = cos(2*pi*X[0]) * cos(4*pi*X[1]) * exp(X[1]) * exp(2*X[2]) + drand48()*10.0;
   fOut = fopen(argv[2], "w");
   fprintf(fOut, " %24.16e\n", Y);
   fclose(fOut);   
}

