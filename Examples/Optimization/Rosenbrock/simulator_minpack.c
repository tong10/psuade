#include <math.h>
#include <stdio.h>
#include <stdlib.h>
int main(int argc, char **argv)
{
   int    n, i;
   double *X, Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8, Y9, YA;
   char   lineIn[100], stringPtr[100], equal[2];

   FILE  *fIn  = fopen(argv[1], "r");
   FILE  *fOut;
   if (fIn == NULL)
   {
      printf("RosenbrockM ERROR - cannot open in/out files.\n");
      exit(1);
   }
   fscanf(fIn, "%d",  &n);
   X = (double *) malloc(n * sizeof(double));
   for (i = 0; i < n; i++) fscanf(fIn, "%lg", &X[i]);
   Y1 = pow(1.0 - X[0],2.0);
   Y2 = 100.0 * pow(X[1] - X[0]*X[0],2.0);
   Y3 = pow(1.0 - X[1],2.0);
   Y4 = 100.0 * pow(X[2] - X[1]*X[1],2.0);
   Y5 = pow(1.0 - X[2],2.0);
   Y6 = 100.0 * pow(X[3] - X[2]*X[2],2.0);
   Y7 = pow(1.0 - X[3],2.0);
   Y8 = 100.0 * pow(X[4] - X[3]*X[3],2.0);
   Y9 = pow(1.0 - X[4],2.0);
   YA = 100.0 * pow(X[5] - X[3]*X[3],2.0);
   fOut = fopen(argv[2], "w");
   fprintf(fOut, "%24.16e\n", Y1);
   fprintf(fOut, "%24.16e\n", Y2);
   fprintf(fOut, "%24.16e\n", Y3);
   fprintf(fOut, "%24.16e\n", Y4);
   fprintf(fOut, "%24.16e\n", Y5);
   fprintf(fOut, "%24.16e\n", Y6);
   fprintf(fOut, "%24.16e\n", Y7);
   fprintf(fOut, "%24.16e\n", Y8);
   fprintf(fOut, "%24.16e\n", Y9);
   fprintf(fOut, "%24.16e\n", YA);
   free(X);
   fclose(fIn);   
   fclose(fOut);   
}

