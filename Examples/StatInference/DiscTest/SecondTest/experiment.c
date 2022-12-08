#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

double getClock();

int main(int argc, char **argv)
{
   int    ii, nInputs, ntime;
   double X1, X2, X3, Y, dtime;
   FILE   *fOut;
                                                                                
   dtime = getClock();
   ntime = (int) (dtime * 1.0E6);
   ntime = ntime % 10000;
   dtime = getClock();
   ntime = ntime * 10000 + (int) (dtime * 1.0E6);
   srand48((long) ntime);
   X2 = 0.9;
   fOut = fopen("expdata", "w");
   fprintf(fOut, "PSUADE_BEGIN\n");
   fprintf(fOut, "21 1 1 1\n");
   for (ii = 0; ii < 21; ii++)
   {
     fprintf(fOut, " %d ", ii+1);
     X1 = 0.6 + 0.015 * ii;
     X3 = drand48();
     Y = X1 + X2 + X1 * X2 + X2 * X3 * 0.5;
     fprintf(fOut, " %24.16e %24.16e 0.5\n", X1, Y);
   }
   fprintf(fOut, "PSUADE_END\n");
   fclose(fOut);
}

double getClock()
{
   double time_i;
   double time_d;
   struct timeval tp;
   struct timezone tzp;
   gettimeofday(&tp,&tzp);
   time_i = tp.tv_sec % 10;
   time_d = (double) time_i;
   time_i = tp.tv_usec;
   time_d = time_d + (double) time_i / 1000000.0;
   return(time_d);
}


