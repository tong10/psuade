#include <math.h>
#include <stdio.h>
#include <stdlib.h>

main(int argc, char **argv)
{
   int    count, i;
   double X[3], Y, pi=3.14159;
   FILE   *fIn  = fopen(argv[1], "r");
   FILE   *fOut;

   if (fIn == NULL)
   {
      printf("Simulator ERROR - cannot open in/out files.\n");
      exit(1);
   }
   fscanf(fIn, "%d", &count);
   if (count != 3)
   {
      printf("Simulator ERROR - invalid nInputs.\n");
      exit(1);
   }
   for (i = 0; i < 3; i++) fscanf(fIn, "%lg", &X[i]);
   Y = sin(3*pi*X[0]) + sin(3*pi*X[1]) + sin(3*pi*X[2]);
   fOut = fopen(argv[2], "w");
   fprintf(fOut, "%24.16e\n", Y);
   fclose(fIn);   
   fclose(fOut);   
}

