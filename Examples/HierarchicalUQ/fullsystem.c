#include <math.h>
#include <stdio.h>
#include <stdlib.h>

main(int argc, char **argv)
{
   int    count, i;
   double theta[2], Y, R, gamma, pi=3.1459;
   FILE   *fIn  = fopen(argv[1], "r");
   FILE   *fOut = fopen(argv[2], "w");
   char   lineIn[100], stringPtr[100], equal[2];

   if (fIn == NULL || fOut == NULL)
   {
      printf("Simulator ERROR - cannot open in/out files.\n");
      exit(1);
   }
   fscanf(fIn, "%d", &count);
   for (i = 0; i < 2; i++) fscanf(fIn, "%lg", &theta[i]);
   gamma = atan2(theta[1],theta[0]);
   R = sqrt(theta[0] * theta[0] + theta[1] * theta[1]);
   Y = (0.8 * R + 0.35 * sin(2.4*pi*R/sqrt(2.0))) * (1.5 * sin(1.3*gamma));
   fprintf(fOut, "%e\n", Y);
   fclose(fIn);
   fclose(fOut);
}

