#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
   int    nInputs;
   double X, Y, error;
   FILE   *fIn  = fopen(argv[1], "r");
   FILE   *fOut;
                                                                                
   if (fIn == NULL)
   {
      printf("Simulator ERROR - cannot open in/out files.\n");
      exit(1);
   }
   fscanf(fIn, "%d", &nInputs);
   fscanf(fIn, "%lg", &X);
   Y = X;
   fOut = fopen(argv[2], "w");
   error = Y - 0.5;
   fprintf(fOut, " %24.16e\n", error);
   fclose(fOut);
}

