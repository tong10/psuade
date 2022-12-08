#include <math.h>
#include <stdio.h>

main()
{
   int    ii, jj, kk, n=9;
   double dx1, dx2, dx3, Y;
   double dx1a, dx2a, dx3a;
   FILE   *fp;
   double pi2 = 3.14159, pi1=0.5*3.14159;

   fp = fopen("train.dat", "w");
   for (ii = 0; ii < n; ii++)
   {
      dx1  = 1.0 / (n-1.0) * ii;
      dx1a = dx1 * pi2 - pi1;
      for (jj = 0; jj < n; jj++)
      {
         dx2 = 1.0 / (n-1.0) * jj;
         dx2a = dx2 * pi2 - pi1;
         for (kk = 0; kk < n; kk++)
         {
            dx3 = 1.0 / (n-1.0) * kk;
            dx3a = dx3 * pi2 - pi1;
            Y=sin(dx1a)+7.0*sin(dx2a)*sin(dx2a)+0.1*dx3a*dx3a*dx3a*dx3a*sin(dx1a);
            fprintf(fp, "%e 1:%e 2:%e 3:%e\n", Y, dx1, dx2, dx3);
         }
      }
   }
   fclose(fp);

   n = 5;
   Y = 0.0;
   fp = fopen("test.dat", "w");
   FILE *fp2 = fopen("test.m", "w");
   for (ii = 0; ii < n; ii++)
   {
      dx1  = 1.0 / (n-1.0) * ii;
      dx1a = dx1 * pi2 - pi1;
      for (jj = 0; jj < n; jj++)
      {
         dx2  = 1.0 / (n-1.0) * jj;
         dx2a = dx2 * pi2 - pi1;
         for (kk = 0; kk < n; kk++)
         {
            dx3  = 1.0 / (n-1.0) * kk;
            dx3a = dx3 * pi2 - pi1;
            Y=sin(dx1a)+7.0*sin(dx2a)*sin(dx2a)+0.1*dx3a*dx3a*dx3a*dx3a*sin(dx1a);
            fprintf(fp, "%e 1:%e 2:%e 3:%e\n", Y, dx1, dx2, dx3);
            fprintf(fp2, "%e %e %e %e\n", dx1, dx2, dx3, Y);
         }
      }
   }
   fclose(fp);
   fclose(fp2);
}

