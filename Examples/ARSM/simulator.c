#include <math.h>
#include <stdio.h>
#include <stdlib.h>

main(int argc, char **argv)
{
   int    i, nInputs, count;
   double X[2], Y;
   FILE   *fIn;
   FILE   *fOut;

   fIn  = fopen(argv[1], "r");
   if (fIn == NULL)
   {
      printf("Simulator ERROR - cannot open in/out files.\n");
      exit(1);
   }
   fscanf(fIn, "%d", &nInputs);
   for (i = 0; i < nInputs; i++) fscanf(fIn, "%lg", &X[i]);
   Y = 0;
   if (X[0] > 0.5) 
      Y = Y + pow(X[0]-0.5,2.0);
   if (X[1] > 0.5)
      Y = Y + pow(X[1]-0.5,2.0);
/*
   Y = X[0] + X[1];
*/
   fOut = fopen(argv[2], "w");
   fprintf(fOut, " %24.16e\n", Y);
   fclose(fOut);   
}

