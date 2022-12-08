#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

int Normal(int nInputs, double *x, double *x2);
double getClock();


main(int argc, char **argv)
{
   int    count, i, j, k, l, ntime;
   double X[20], W[20], B[20], B2[20], B3[20], Y, B0, dtime;
   FILE   *fIn  = fopen(argv[1], "r");
   FILE   *fOut;

   if (fIn == NULL)
   {
      printf("Simulator ERROR - cannot open in/out files.\n");
      exit(1);
   }
   fscanf(fIn, "%d", &count);
   if (count != 20)
   {
      printf("Simulator ERROR - invalid nInputs.\n");
      exit(1);
   }
   for (i = 0; i < count; i++) fscanf(fIn, "%lg", &X[i]);

   dtime = getClock();
   ntime = (int) (dtime * 1.0E6);
   ntime = ntime % 10000;
   dtime = getClock();
   ntime = ntime * 10000 + (int) (dtime * 1.0E6);
   srand48((long) ntime);

   for (i = 0; i < count; i++) W[i] = 2 * (X[i] - 0.5);
   W[2] = 2.0 * ((1.1 * X[2])/(X[2] + 0.1) - 0.5);
   W[4] = 2.0 * ((1.1 * X[4])/(X[4] + 0.1) - 0.5);
   W[6] = 2.0 * ((1.1 * X[6])/(X[6] + 0.1) - 0.5);

   for (i = 0; i < 10; i++) B[i] = 20.0;
   for (i = 10; i < 20; i++) B[i] = drand48();
   Normal(10, &B[10], B2);
   for (i = 10; i < 20; i++) B[i] = B2[i-10];

   B0 = drand48();
   Normal(1, &B0, &Y);
   for (i = 0; i < 20; i++) Y += B[i] * W[i];
   for (i = 0; i < 6; i++)
      for (j = i+1; j < 6; j++) Y += -15.0 * W[i] * W[j];
   for (i = 0; i < 6; i++)
   {
      for (j = 0; j < 14; j++) B2[j] = drand48(); 
      Normal(14, B2, B3);
      for (j = 6; j < 20; j++) Y += B3[j-6] * W[i] * W[j];
   }
   for (i = 6; i < 20; i++)
   {
      for (j = 0; j < 15; j++) B2[j] = drand48(); 
      Normal(15, B2, B3);
      for (j = i+1; j < 20; j++) Y += B3[j-6] * W[i] * W[j];
   }
   for (i = 0; i < 5; i++)
      for (j = i+1; j < 5; j++)
         for (k = j+1; k < 5; k++)
            Y += -10.0 * W[i] * W[j] * W[k];
   for (i = 0; i < 4; i++)
      for (j = i+1; j < 4; j++)
         for (k = j+1; k < 4; k++)
            for (l = k+1; l < 4; l++)
               Y += 5.0 * W[i] * W[j] * W[k] * W[l];

   fOut = fopen(argv[2], "w");
   fprintf(fOut, "%24.16e\n", Y);
   fclose(fIn);   
   fclose(fOut);   
}

int Normal(int nInputs, double *x, double *x2)
{
   int    i, nParts, finish;
   double step, mu=0.0, sigma=1.0, pi=3.1415928, exp1;
   double left, right, area, yAcc, total, coef, exp2, xplus, xminus;

   coef   = 1.0 / ( sigma * sqrt(2.0*pi) );
   left   = mu - 5.0 * sigma; 
   right  = mu + 5.0 * sigma; 
   step   = (right - left) / 10000.0;
   nParts = 10001;
   finish = 0;
   total   = 0.0;
   for ( i = 0; i < nParts+1; i++ )
   {
      xplus  = left + step * ( i + 0.5 );
      xminus = left + step * ( i - 0.5 );
      exp1 = - (xminus - mu) * (xminus - mu) / (2.0 * sigma * sigma);
      exp2 = - (xplus  - mu) * (xplus  - mu) / (2.0 * sigma * sigma);
      area = (exp(exp1) + exp(exp2)) * coef * step * 0.5;
      total += area;
   }
   yAcc = 0.0;
   for ( i = 0; i < nParts; i++ )
   {
      xplus  = left + step * ( i + 0.5 );
      xminus = left + step * ( i - 0.5 );
      exp1 = - (xminus - mu) * (xminus - mu) / (2.0 * sigma * sigma);
      exp2 = - (xplus  - mu) * (xplus  - mu) / (2.0 * sigma * sigma);
      area = (exp(exp1) + exp(exp2)) * coef * step * 0.5;
      yAcc += area;
      if ( yAcc/total >= x[0] && ((finish & 1) == 0)) 
      {
         x2[0] = left + step * i;
         finish |= 1;
      }
      if ( nInputs > 1 && yAcc/total >= x[1] && ((finish & 2) == 0)) 
      {
         x2[1] = left + step * i;
         finish |= 2;
      } else if (nInputs <= 1) finish |= 2;
      if ( nInputs > 2 && yAcc/total >= x[2] && ((finish & 4) == 0)) 
      {
         x2[2] = left + step * i;
         finish |= 4;
      } else if (nInputs <= 2) finish |= 4;
      if ( nInputs > 3 && yAcc/total >= x[3] && ((finish & 8) == 0)) 
      {
         x2[3] = left + step * i;
         finish |= 8;
      } else if (nInputs <= 3) finish |= 8;
      if ( nInputs > 4 && yAcc/total >= x[4] && ((finish & 16) == 0)) 
      {
         x2[4] = left + step * i;
         finish |= 16;
      } else if (nInputs <= 4) finish |= 16;
      if ( nInputs > 5 && yAcc/total >= x[5] && ((finish & 32) == 0)) 
      {
         x2[5] = left + step * i;
         finish |= 32;
      } else if (nInputs <= 5) finish |= 32;
      if ( nInputs > 6 && yAcc/total >= x[6] && ((finish & 64) == 0)) 
      {
         x2[6] = left + step * i;
         finish |= 64;
      } else if (nInputs <= 6) finish |= 64;
      if ( nInputs > 7 && yAcc/total >= x[7] && ((finish & 128) == 0)) 
      {
         x2[7] = left + step * i;
         finish |= 128;
      } else if (nInputs <= 7) finish |= 128;
      if ( nInputs > 8 && yAcc/total >= x[8] && ((finish & 256) == 0)) 
      {
         x2[8] = left + step * i;
         finish |= 256;
      }
      else if (nInputs <= 8) finish |= 256;
      if ( nInputs > 9 && yAcc/total >= x[9] && ((finish & 512) == 0)) 
      {
         x2[9] = left + step * i;
         finish |= 512;
      }
      else if (nInputs <= 9) finish |= 512;
      if ( finish == 1023 ) break;
   }
   if ( nInputs > 0 && (finish & 1) == 0 ) x2[0]  = right;
   if ( nInputs > 1 && (finish & 2) == 0 ) x2[1]  = right;
   if ( nInputs > 2 && (finish & 4) == 0 ) x2[2]  = right;
   if ( nInputs > 3 && (finish & 8) == 0 ) x2[3]  = right;
   if ( nInputs > 4 && (finish & 16) == 0 ) x2[4] = right;
   if ( nInputs > 5 && (finish & 32) == 0 ) x2[5] = right;
   if ( nInputs > 6 && (finish & 64) == 0 ) x2[6] = right;
   if ( nInputs > 7 && (finish & 128) == 0 ) x2[7] = right;
   if ( nInputs > 8 && (finish & 256) == 0 ) x2[8] = right;
   if ( nInputs > 9 && (finish & 512) == 0 ) x2[9] = right;
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

