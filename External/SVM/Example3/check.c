#include <math.h>
#include <stdio.h>

main()
{
   int    ii, n=101;
   double dx1, Y, ddata, err, error, maxerr;
   FILE   *fp = fopen("pred.dat", "r");

   error = 0.0;
   maxerr = 0.0;
   for (ii = 0; ii < n; ii++)
   {
      dx1 = 1.0 / (n-1.0) * ii;
      Y   = dx1 * dx1 * dx1 * dx1 * dx1;
      fscanf(fp, "%lg \n", &ddata);
      err = (Y - ddata) * (Y - ddata);
      if (err > maxerr) maxerr = err;
      error += err;
      printf("%5d : %12.4e    %12.4e   %12.4e \n", ii, Y, ddata, err);
   }
   printf("L-1 norm squared error = %e\n", maxerr);
   printf("L-2 norm squared error = %e\n", error);
   fclose(fp);
}

