#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double gaussrand();
double genNormal();      
#define PABS(x) ((x > 0) ? (x) : -(x))

int main(int argc, char **argv)
{
   int    ii, nInputs;
   double X[2], Y, beta=0.65, A=20.0;
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
   fOut = fopen(argv[2], "w");
   /* Y = beta * X[0] / (1.0 + X[0] / A) + 0.01 * gaussrand(); */
   Y = beta * X[0] / (1.0 + X[0] / A) + genNormal();
   fprintf(fOut, " %24.16e\n", Y);
   fclose(fOut);
}


#define PI 3.141592654

double gaussrand()
{
   static double U, V;
   static int phase = 0;
   double Z;

   if(phase == 0) {
   	U = (rand() + 1.) / (RAND_MAX + 2.);
   	V = rand() / (RAND_MAX + 1.);
   	Z = sqrt(-2 * log(U)) * sin(2 * PI * V);
   } else
   	Z = sqrt(-2 * log(U)) * cos(2 * PI * V);
   
   phase = 1 - phase;

   return Z;
}

double genNormal()      
{      
   double ddata, xlo=-5, xhi=5,ylo,yhi, iroot2, scale=10;
   double outdata, mean=0.0, stdev=0.1, xmi, ymi;

   iroot2 = sqrt(0.5)/stdev;
   ddata = drand48();
   ylo = 0.5 * (1.0 + erf((xlo-mean)*iroot2));
   yhi = 0.5 * (1.0 + erf((xhi-mean)*iroot2));
   ddata = (ddata + 5) / scale  * (yhi - ylo) + ylo;
   if      (ddata <= ylo) outdata = xlo;
   else if (ddata >= yhi) outdata = xhi;
   else
   {
      while (PABS(ddata-ylo) > 1.0e-12 || PABS(ddata-yhi) > 1.0e-12)
      {
         xmi = 0.5 * (xhi + xlo);
         ymi = 0.5 * (1.0 + erf((xmi-mean)*iroot2));
         if (ddata > ymi)
         {
            xlo = xmi;
            ylo = ymi;
         }
         else
         {
            xhi = xmi;
            yhi = ymi;
         }
      }
      if (PABS(ddata-ylo) < PABS(ddata-yhi)) outdata = xlo;
      else                                   outdata = xhi;
   }
   return outdata;
}

