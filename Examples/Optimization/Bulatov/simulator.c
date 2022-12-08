#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define PABS(x)  ((x) > 0 ? x : -(x))

main(int argc, char **argv)
{
   int    nInputs, i, j, npts, *CC;
   double X[2], Y, *XX, *YY;
   FILE  *fIn = fopen(  argv[1], "r");
   FILE  *fOut, *fp;

   double a1=-1.4763;
   double a11=-6.2554;
   double a111=28.6679;
   double a1111=-27.5092;
   double a11111=2.4055;
   double a111111=4.0179;
   double a2=-0.1594;
   double a22= 5.9116;
   double a222=0;
   double a2222=-34.0007;
   double a22222=43.5076;
   double a222222=-15.2657;
   if (fIn == NULL)
   {
      printf("low_fidelity ERROR - cannot open in/out files.\n");
      exit(1);
   }
   fscanf(fIn, "%d",  &nInputs);
   if (nInputs != 2)
   {
      printf("ERROR : nInputs != 2.\n");
      fclose(fIn);
      fclose(fOut);
      exit(1);
   }
   for (i = 0; i < 2; i++) fscanf(fIn, "%lg", &X[i]);

   Y = a1*X[0]+a11*pow(X[0],2)+a111*pow(X[0],3)+a1111*pow(X[0],4)+
       a11111*pow(X[0],5)+a111111*pow(X[0],6)+
       a2*X[1]+a22*pow(X[1],2)+a222*pow(X[1],3)+a2222*pow(X[1],4)+
       a22222*pow(X[1],5)+a222222*pow(X[1],6);
   fOut = fopen(argv[2], "w");
   fprintf(fOut, "%24.16f\n", Y);
   fclose(fOut);
}


