#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
   int    ii, jj, nPts=5, count;
   double X1, X2, width;
   FILE   *fOut;
                                                                                
   fOut = fopen("CandidateSet0", "w");
   fprintf(fOut, "%2d 2\n", nPts*nPts);   
   width = 1.0 / (nPts - 1);
   count = 1;
   for (ii = 0; ii < nPts; ii++)
   {
     X1 = width * ii;
     for (jj = 0; jj < nPts; jj++)
     {
       X2 = width * jj;
       fprintf(fOut, "%2d %18.10e %18.10e\n",
               count, X1, X2);   
       count++;
     }
   }  
   fclose(fOut);
}

