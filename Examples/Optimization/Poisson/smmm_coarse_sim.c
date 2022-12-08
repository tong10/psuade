#include <math.h>
#include <stdio.h>
#include <stdlib.h>


int matvec(int N, int *IA, int *JA, double *AA, double *X, double *Y);
double inner(int N, double *X, double *Y);

int main(int argc, char **argv)
{
   int    nInputs, n, N, i, j, row, nnz, *IA, *JA, iter;
   double S[4], *AA, *B, *R, *P, *AP, *U, h;
   double rho, rhom1, alpha, beta, sigma, rnorm, rnorm0;
   double y[4], cx[4];

   FILE  *fIn     = fopen(    argv[1], "r");
   FILE  *fOut    = fopen(    argv[2], "w");

   if (fIn == NULL || fOut == NULL)
   {
      printf("Poisson ERROR - cannot open in/out files.\n");
      exit(1);
   }
   fscanf(fIn, "%d",  &nInputs);
   for (i = 0; i < nInputs; i++) fscanf(fIn, "%lg", &S[i]);

   {
      /* set up problem */
      n = 33;
      N = n * n;
      IA = (int *) malloc((N + 1) * sizeof(int));
      JA = (int *) malloc(5 * N * sizeof(int));
      AA = (double *) malloc(5 * N * sizeof(double));
      row = 0;
      nnz = 0;
      for (j = 0; j < n; j++)
      {
         for (i = 0; i < n; i++)
         {
            IA[row] = nnz;
            if (j > 0)
            {
               JA[nnz] = row - n;
               AA[nnz++] = -1.0;
            }
            if (i > 0)
            {
               JA[nnz] = row - 1;
               AA[nnz++] = -1.0;
            }
            JA[nnz] = row;
            AA[nnz++] = 4.0;
            if (i < n-1)
            {
               JA[nnz] = row + 1;
               AA[nnz++] = -1.0;
            }
            if (j < n-1)
            {
               JA[nnz] = row + n;
               AA[nnz++] = -1.0;
            }
            row++;
         }
      }
      IA[N] = nnz;
      B = (double *) malloc(N * sizeof(double));
      h = 1.0 / (double) (n - 1);

      for (i = 0; i < N; i++) B[i] = 0.0;
      //B[ 2] = S[0];
      //B[11] = S[1];
      //B[(n-1)/2] = S[0];
      //B[n*(n-1)/2+(n-1)/4] = S[1];
      //B[n*n-1] = S[2];
      B[0] = S[0];
      B[n-1] = S[1];
      B[n*(n-1)] = S[2];
      B[n*n-1] = S[3];

      /* CG */

      R = (double *) malloc(N * sizeof(double));
      P = (double *) malloc(N * sizeof(double));
      AP = (double *) malloc(N * sizeof(double));
      U = (double *) malloc(N * sizeof(double));
      for (i = 0; i < N; i++) R[i] = B[i];
      for (i = 0; i < N; i++) P[i] = 0.0;
      iter = 0;
      rho = inner(N, R, R);
      rnorm0 = sqrt(rho);
      rnorm  = rnorm0;
      /* printf("CG initial rnorm = %e\n", rnorm0); */
      while (rnorm/rnorm0 > 1.0e-8)
      {
         iter++;
         if (iter == 1) beta = 0.0;
         else           beta = rho / rhom1;
         for (i = 0; i < N; i++) P[i] = R[i] + beta * P[i];
         matvec(N, IA, JA, AA, P, AP);
         sigma = inner(N, P, AP);
         alpha = rho / sigma;
         for (i = 0; i < N; i++) U[i] += alpha * P[i];
         for (i = 0; i < N; i++) R[i] -= alpha * AP[i];
         rhom1 = rho;
         rho = inner(N, R, R);
         rnorm = sqrt(rho);
         /* printf("CG iter = %5d, rnorm = %e\n", iter, rnorm); */
      }

      /* clean up */
      free(B);
      free(IA);
      free(JA);
      free(AA);
      free(R);
      free(P);
      free(AP);
      free(U);
   }

   /* OUTPUT: MODEL RESPONSE */
   double scale = h * h;
   cx[0] = U[n*(n-1)/4+(n-1)/4] / scale;
   cx[1] = U[n*(n-1)/4+(n-1)/4*3] / scale;
   cx[2] = U[n*(n-1)/4*3+(n-1)/4] / scale;
   cx[3] = U[n*(n-1)/4*3+(n-1)/4*3] / scale;
//printf("%e %e %e %e\n",cx[0],cx[1],cx[2],cx[3]);
   for (i = 0; i < 4 ; i++) fprintf(fOut, "%24.16e\n", cx[i]);

   fclose(fIn);
   fclose(fOut);
}

int matvec(int N, int *IA, int *JA, double *AA, double *X, double *Y)
{
   int    i, j;
   double dtemp;

   for (i = 0; i < N; i++)
   {
      dtemp = 0.0;
      for (j = IA[i]; j < IA[i+1]; j++)
         dtemp += AA[j] * X[JA[j]];
      Y[i] = dtemp;
   }
   return 0;
}

double inner(int N, double *X, double *Y)
{
   int i;
   double dtemp;

   dtemp = 0.0;
   for (i = 0; i < N; i++) dtemp += X[i] * Y[i];
   return dtemp;
}
 
