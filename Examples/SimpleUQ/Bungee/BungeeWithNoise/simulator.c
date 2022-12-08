#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
   int    i, nInputs, count;
   double X[40], Y, kel=1.5, g=9.8, ddata;
   FILE   *fIn  = fopen(argv[1], "r");
   FILE   *fOut;

   if (fIn == NULL)
   {
      printf("Simulator ERROR - cannot open in/out files.\n");
      exit(1);
   }
   srand(5415417);
   fscanf(fIn, "%d", &nInputs);
   for (i = 0; i < nInputs; i++) fscanf(fIn, "%lg", &X[i]);
   Y = X[0] - 2.0 * X[1] * g / (kel * X[2]);
   for (i = 3; i < nInputs; i++)
   {
      ddata = drand48();
      /*printf("%2d : %e\n", i+1, ddata);*/
      Y += 0.01 * (ddata - 0.5) * X[i];
   }
   fOut = fopen(argv[2], "w");
   fprintf(fOut, " %24.16e\n", Y);
   fclose(fOut);   
}

