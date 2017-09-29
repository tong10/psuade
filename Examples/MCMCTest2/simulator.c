#include <math.h>
#include <stdio.h>
#include <stdlib.h>

main(int argc, char **argv)
{
   int    nInputs;
   double X1, X2, Y;
   FILE   *fIn  = fopen(argv[1], "r");
   FILE   *fOut;
                                                                                
   if (fIn == NULL)
   {
      printf("Simulator ERROR - cannot open in/out files.\n");
      exit(1);
   }
   fscanf(fIn, "%d", &nInputs);
   fscanf(fIn, "%lg", &X1);
   fscanf(fIn, "%lg", &X2);
   /* X1 = [0 2] */
   /* X2 = [0 10] */
   Y = pow(X1 - 0.3,2.0) * pow(X1 - 1.7,2) + pow(X2 - 0.3, 2) * pow(X2 - 1.7,2);
   fOut = fopen(argv[2], "w");
   fprintf(fOut, " %24.16e\n", Y);
   fclose(fOut);
}

