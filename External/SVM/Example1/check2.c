#include <stdio.h>

main()
{
   int    ii, jj;
   double dx1, dx2, Y, ddata, error;
   FILE   *fp = fopen("pred2.dat", "r");

   error = 0.0;
   for (ii = 0; ii < 33; ii++)
   {
      dx1 = 1.0 / 32.0 * ii;
      for (jj = 0; jj < 33; jj++)
      {
         fscanf(fp, "%lg \n", &ddata);
         dx2 = 1.0 / 32.0 * jj;
         Y = dx1 * dx1 * dx1 + dx2 * dx2 * dx2;
         printf("%4d : %12.4e     %12.4e      %12.4e\n", ii*33+jj+1, Y, ddata, error);
         error += (Y - ddata) * (Y - ddata);
      }
   }
   printf("L-2 norm squared error = %e\n", error);
   fclose(fp);
}

