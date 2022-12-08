#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

main(int argc, char **argv)
{
   int    nInputs;
   double Y, theta2;
   FILE   *fIn  = fopen(argv[1], "r");
   FILE   *fOut;
                                                                                
   if (fIn == NULL)
   {
      printf("Simulator2 ERROR - cannot open in/out files.\n");
      exit(1);
   }
   fscanf(fIn, "%d", &nInputs);
   fscanf(fIn, "%lg", &theta2);
   fclose(fIn);
   Y = theta2 + theta2 * theta2;
   fOut  = fopen(argv[2], "w");
   fprintf(fOut, " %24.16e\n", Y);
   fclose(fOut);
}

