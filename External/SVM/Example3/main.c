#include <math.h>
#include <stdio.h>

main()
{
   int    ii, n=5;
   double dx1, Y, pi2=2.0*3.14159;
   FILE   *fp;

   fp = fopen("train.dat", "w");
   for (ii = 0; ii < n; ii++)
   {
      dx1  = 1.0 / (n-1.0) * ii;
      Y = sin(pi2 *dx1) * sin(pi2 *dx1);
      fprintf(fp, "%e 1:%e\n", Y, dx1);
   }
   fclose(fp);

   n = 101;
   Y = 0.0;
   fp = fopen("test.dat", "w");
   for (ii = 0; ii < n; ii++)
   {
      dx1  = 1.0 / (n-1.0) * ii;
      Y = sin(pi2 *dx1) * sin(pi2 *dx1);
      fprintf(fp, "%e 1:%e\n", dx1, dx1);
   }
   fclose(fp);
}

