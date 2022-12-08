#include <math.h>
#include <stdio.h>

main()
{
   int    ii, jj, kk, n=5;
   double dx1, dx2, dx3, Y, ddata, err, error;
   double dx1a, dx2a, dx3a, maxerr;
   FILE   *fp = fopen("pred.dat", "r");
   double pi2 = 3.14159, pi1=0.5*3.14159;

   error = 0.0;
   maxerr = 0.0;
   for (ii = 0; ii < n; ii++)
   {
      dx1  = 1.0 / (n-1) * ii;
      dx1a = pi2 * dx1 - pi1;
      for (jj = 0; jj < n; jj++)
      {
         dx2  = 1.0 / (n-1) * jj;
         dx2a = pi2 * dx2 - pi1;
         for (kk = 0; kk < n; kk++)
         {
            dx3  = 1.0 / (n-1) * kk;
            dx3a = pi2 * dx3 - pi1;
            Y=sin(dx1a)+7.0*sin(dx2a)*sin(dx2a)+0.1*dx3a*dx3a*dx3a*dx3a*sin(dx1a);
            fscanf(fp, "%lg \n", &ddata);
            err = (Y - ddata) * (Y - ddata);
            if (err > maxerr) maxerr = err;
            error += err;
            printf("%5d : %12.4e    %12.4e   %12.4e \n", ii*n*n+jj*n+kk+1, Y, ddata, err);
         }
      }
   }
   printf("L-1 norm squared error = %e\n", maxerr);
   printf("L-2 norm squared error = %e\n", error);
   fclose(fp);
}

