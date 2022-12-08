#include <math.h>
#include <stdio.h>
#include <stdlib.h>

main(int argc, char **argv)
{
   int    count, i;
   double X[100], Y, Xm;
   FILE   *fIn  = fopen(argv[1], "r");
   FILE   *fOut;

   if (fIn == NULL)
   {
      printf("Simulator ERROR - cannot open in/out files.\n");
      exit(1);
   }
   fscanf(fIn, "%d", &count);
   if (count != 100)
   {
      printf("Simulator ERROR - invalid nInputs.\n");
      exit(1);
   }
   for (i = 0; i < 100; i++) fscanf(fIn, "%lg", &X[i]);
   Xm = 0.0;
   for (i = 0; i < 30; i++) Xm += X[i];
   Xm /= 30.0;
   Y = 0.0;
   for (i = 0; i < 30; i++)
      Y += exp(5.5 * X[i] - 1.5 * Xm);

   fOut = fopen(argv[2], "w");
   fprintf(fOut, "%24.16e\n", Y);
   fclose(fIn);   
   fclose(fOut);   
}

