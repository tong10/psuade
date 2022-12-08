#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
   int    ii, jj, kk, nPts=5, count;
   double X1, X2, X3, width;
   FILE   *fOut;
                                                                                
   fOut = fopen("CandidateSet0", "w");
   fprintf(fOut, "%2d 3\n", nPts*nPts*nPts);   
   width = 1.0 / (nPts - 1);
   count = 1;
   for (ii = 0; ii < nPts; ii++)
   {
     X1 = width * ii;
     for (jj = 0; jj < nPts; jj++)
     {
       X2 = width * jj;
       for (kk = 0; kk < nPts; kk++)
       {
         X3 = width * kk;
         fprintf(fOut, "%2d %18.10e %18.10e %18.10e\n",
                 count, X1, X2, X3);   
         count++;
       }
     }
   }  
   fclose(fOut);
}

