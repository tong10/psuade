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
   Y = X[0] + X[1];
   fOut = fopen(argv[2], "w");
   fprintf(fOut, " %24.16e\n", Y);
   fclose(fOut);   
}

