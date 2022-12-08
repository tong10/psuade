#include <stdio.h>

main()
{
   int    ii, jj;
   double dx1, dx2, Y;
   FILE   *fp = fopen("train2.dat", "w");

   for (ii = 0; ii < 10; ii++)
   {
      dx1 = 1.0 / 9.0 * ii;
      for (jj = 0; jj < 10; jj++)
      {
         dx2 = 1.0 / 9.0 * jj;
         Y = dx1 * dx1 * dx1 + dx2 * dx2 * dx2;
         fprintf(fp, "%e 1:%e 2:%e\n", Y, dx1, dx2);
      }
   }
   fclose(fp);

   fp = fopen("test2.dat", "w");
   for (ii = 0; ii < 33; ii++)
   {
      dx1 = 1.0 / 32.0 * ii;
      for (jj = 0; jj < 33; jj++)
      {
         dx2 = 1.0 / 32.0 * jj;
         Y = dx1 * dx1 * dx1 + dx2 * dx2 * dx2;
         fprintf(fp, "%e 1:%e 2:%e\n", Y, dx1, dx2);
      }
   }
   fclose(fp);
}

