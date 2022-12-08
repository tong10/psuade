#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
   int    count, i;
   double X[2], Y, d1, scale=4, pi=3.1415928;
   FILE   *fIn  = fopen(argv[1], "r");
   FILE   *fOut = fopen(argv[2], "w");

   if (fIn == NULL || fOut == NULL)
   {
      printf("Simulator ERROR - cannot open in/out files.\n");
      exit(1);
   }
   fscanf(fIn, "%d", &count);
   for (i = 0; i < 2; i++) fscanf(fIn, "%lg", &X[i]);
   if (X[1] < 3*pi) Y = 0;
   else
   {
     d1 = 2*pi-0.25*sin(scale*(X[1]-3*pi));
     if (X[0] - d1 > 0) Y = 1;
     else               Y = 0;
   }
   fprintf(fOut, "%e\n", Y);
   fclose(fIn);   
   fclose(fOut);   
}

