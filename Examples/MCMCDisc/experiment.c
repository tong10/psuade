#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

double getClock();

main(int argc, char **argv)
{
   int    ii, nInputs, ntime;
   double X[2], Y, X1, X2, X3, dtime;
   FILE   *fIn  = fopen(argv[1], "r");
   FILE   *fOut;
                                                                                
   if (fIn == NULL)
   {
      printf("Experiment ERROR - cannot open in/out files.\n");
      exit(1);
   }
   fscanf(fIn, "%d", &nInputs);
   for (ii = 0; ii < nInputs; ii++) fscanf(fIn, "%lg", &X[ii]);
   fclose(fIn);
   dtime = getClock();
   ntime = (int) (dtime * 1.0E6);
   ntime = ntime % 10000;
   dtime = getClock();
   ntime = ntime * 10000 + (int) (dtime * 1.0E6);
   srand48((long) ntime);
   X1 = X[0];
   X2 = 0.9;
   X3 = drand48();
   fOut = fopen(argv[2], "w");
   /*Y = X1 + X2 + X1 * X2 + X2 * X3 * 0.5;
   */
   Y = X1 + X2 + X1 * X2;
   fprintf(fOut, " %24.16e\n", Y);
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


