#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
   int    ii, jj, nSam, nInp;
   double X1, X2, Y;
   FILE   *fIn, *fOut;
                                                                                
   if (argc < 3)
   {
     printf("Experiment ERROR: not enough arguments.\n");
     exit(1);
   }
   fIn = fopen(argv[1], "r");
   fOut = fopen(argv[2], "w");
   fscanf(fIn, "%d %d", &nSam, &nInp);
   fprintf(fOut, "%d %d\n", nSam, nInp+2);
   for (ii = 0; ii < nSam; ii++)
   {
     fscanf(fIn, "%d %lg %lg", &jj, &X1, &X2);
     Y = 1.0 + X1 * X1 + X2 * X2 +
         0.5 * exp(-10*(X1-0.75)*(X1-0.75)) *
               exp(-10*(X2-0.75)*(X2-0.75)) * 0.5 +
         0.25 * exp(-10*(X1-0.25)*(X1-0.25)) *
                exp(-10*(X2-0.25)*(X2-0.25)) * 0.5;
     Y = Y * (1 + 0.01 * (drand48() - 0.5));
     //fprintf(fOut, "%d %18.10e %18.10e %18.10e %18.10e\n",
     //        ii+1, X1, X2, Y, Y*0.01);   
     fprintf(fOut, "%d %18.10e %18.10e %18.10e 0.05\n",
             ii+1, X1, X2, Y);   
   }  
   fclose(fIn);
   fclose(fOut);
}

