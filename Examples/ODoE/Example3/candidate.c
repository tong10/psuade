#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
   int    ii, jj, nPts=5, count;
   double X1, X2, Y, width;
   FILE   *fOut;
                                                                                
   fOut = fopen("CandidateSet0", "w");
   //fprintf(fOut, "%2d 2\n", nPts*nPts);   
   fprintf(fOut, "%2d 4\n", nPts*nPts);   
   width = 1.0 / (nPts - 1);
   count = 1;
   for (ii = 0; ii < nPts; ii++)
   {
     X1 = width * ii;
     for (jj = 0; jj < nPts; jj++)
     {
       X2 = width * jj;
       Y = 1.0 + X1 * X1 + X2 * X2 +
           0.5 * exp(-10*(X1-0.75)*(X1-0.75)) *
                 exp(-10*(X2-0.75)*(X2-0.75)) * 0.5 +
           0.25 * exp(-10*(X1-0.25)*(X1-0.25)) *
                  exp(-10*(X2-0.25)*(X2-0.25)) * 0.5;
       //fprintf(fOut, "%2d %18.10e %18.10e\n",
       //        count, X1, X2);
       //fprintf(fOut, "%2d %18.10e %18.10e %18.10e %18.10e\n",
       //        count, X1, X2, Y, Y*0.05);
       fprintf(fOut, "%2d %18.10e %18.10e %18.10e 0.05\n",
               count, X1, X2, Y);
       count++;
     }
   }  
   fclose(fOut);
}

