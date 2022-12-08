#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
   int    count, ii;
   double X[100], Y, d=2, ddata, t;
   FILE   *fIn  = fopen(argv[1], "r");
   FILE   *fOut = fopen(argv[2], "w");

   if (fIn == NULL || fOut == NULL)
   {
      printf("Simulator ERROR - cannot open in/out files.\n");
      exit(1);
   }
   fscanf(fIn, "%d", &count);
   if (count != 2)
   {
      printf("Simulator ERROR - wrong nInputs.\n");
      exit(1);
   }
   for (ii = 0; ii < count; ii++) fscanf(fIn, "%lg", &X[ii]);
   Y = 0.0;
   for (ii = 0; ii < count; ii++)
   {
      t = X[ii];
      if ( d == 1 )
         ddata = t;
      else if ( d == 2 )
         ddata = 0.5 * (3 * t * t - 1);
      else if ( d == 3 )
         ddata = 0.5 * (5 * t * t * t - 3 * t);
      else if ( d == 4 )
         ddata = 0.125 * (35 * t * t * t * t - 30 * t * t + 3);
      else if ( d == 5 )
         ddata = 0.125 * (63 * t * t * t * t *t - 70 * t * t * t + 15 * t);
      Y += ddata;
   }
   fprintf(fOut, "%24.16e\n", Y);
   fclose(fIn);   
   fclose(fOut);   
}

