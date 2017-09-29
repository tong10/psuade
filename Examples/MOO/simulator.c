#include <math.h>
#include <stdio.h>
#include <stdlib.h>

main(int argc, char **argv)
{
   int    count, i;
   double X[5], F1, G, H, F2;
   FILE   *fIn  = fopen(argv[1], "r");
   FILE   *fOut;

   if (fIn == NULL)
   {
      printf("Simulator ERROR - cannot open in/out files.\n");
      exit(1);
   }
   fscanf(fIn, "%d", &count);
   if (count != 5)
   {
      printf("Simulator ERROR - wrong nInputs.\n");
      exit(1);
   }
   for (i = 0; i < count; i++) fscanf(fIn, "%lg", &X[i]);
   F1 = X[0];
   fOut = fopen(argv[2], "w");
   fprintf(fOut, "%24.16e\n", F1);
   G = 1.0;
   for (i = 1; i < count; i++) G += 9.0 * X[i] / 4.0;
   H = 1 - (F1 / G) * (F1 / G);
   F2 = G * H;
   fprintf(fOut, "%24.16e\n", F2);
   fclose(fIn);   
   fclose(fOut);   
}

