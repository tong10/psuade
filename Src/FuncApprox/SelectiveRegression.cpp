// ************************************************************************
// Copyright (c) 2007   Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory.
// Written by the PSUADE team.
// All rights reserved.
//
// Please see the COPYRIGHT_and_LICENSE file for the copyright notice,
// disclaimer, contact information and the GNU Lesser General Public License.
//
// PSUADE is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License (as published by the Free Software
// Foundation) version 2.1 dated February 1999.
//
// PSUADE is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the terms and conditions of the GNU General
// Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
// ************************************************************************
// Functions for the class SelectiveRegression
// AUTHOR : CHARLES TONG
// DATE   : 2005
// ************************************************************************
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "SelectiveRegression.h"
#include "Psuade.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "PDFManager.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

extern "C" {
   void dgels_(char *, int *, int *, int *, double *, int *,
               double *, int *, double *, int *, int *);
   void dgesvd_(char *, char *, int *, int *, double *, int *, double *,
               double *, int *, double *, int *, double *, int *, int *);
   void dgetrf_(int *, int *, double *, int *, int *, int *);
   void dgetri_(int *, double *, int *, int *, double*, int *, int *);
}

// ************************************************************************
// Constructor 
// ------------------------------------------------------------------------
SelectiveRegression::SelectiveRegression(int nInputs,int nSamples):
                                         FuncApprox(nInputs,nSamples)
{
   int  i, j, termCheck, length;
   char inStr[300];
   FILE *fp = fopen("selective_regression_file", "r");

   faID_ = PSUADE_RS_REGRU;
   if (fp == NULL)
   {
      printf("SelectiveRegression ERROR: selective_regression_file not found.\n");
      printf("This file is needed to identify selected terms.\n");
      printf("This file should have the following format: \n");
      printf("line   1: PSUADE_BEGIN\n");
      printf("line   2: the number of terms.\n");
      printf("line   3: 1  number of inputs, input list (1-based).\n");
      printf("line   4: 2  number of inputs, input list (1-based).\n");
      printf(".........\n");
      printf("line   m: 1 0 (meaning constant term).\n");
      printf("lastline: END\n");
      exit(1);
   }
   else
   {
      printf("SelectiveRegression INFO: using selective_regression_file.\n");
   }
   fscanf(fp, "%s", inStr);
   if (strcmp(inStr, "PSUADE_BEGIN"))
   {
      printf("SelectiveRegression: wrong format in selective_regression_file.\n");
      printf("line   1: PSUADE_BEGIN\n");
      printf("line   2: the number of terms.\n");
      printf("line   3: 1  number of inputs, input list (1-based).\n");
      printf("line   4: 2  number of inputs, input list (1-based).\n");
      printf(".........\n");
      printf("line   m: 1 0 (meaning constant term).\n");
      printf("lastline: END\n");
      fclose(fp);
      exit(1);
   }
   fscanf(fp, "%d", &numTerms_);
   if (numTerms_ >= nSamples)
   {
      printf("SelectiveRegression ERROR: no. of terms should be < nSamples.\n");
      printf("                    nTerms   = %d\n", numTerms_);
      printf("                    nSamples = %d\n", nSamples);
      printf("Check your selective_regression_file file format as follow: \n");
      printf("line   1: PSUADE_BEGIN\n");
      printf("line   2: the number of terms.\n");
      printf("line   3: 1  number of inputs, input list (1-based).\n");
      printf("line   4: 2  number of inputs, input list (1-based).\n");
      printf(".........\n");
      printf("line   m: 1 0 (meaning constant term).\n");
      printf("lastline: END\n");
      exit(1);
   }
   regCoeffs_ = NULL;
   regStdevs_ = NULL;
   fuzzyC_ = NULL;
   coefTerms_ = new int*[numTerms_];
   for (i = 0; i < numTerms_; i++)
   {
      fscanf(fp, "%d %d", &termCheck, &length);
      if (termCheck != (i+1))
      {
         printf("SelectiveRegression ERROR: %d-th term has index %d.\n",i+1,
                termCheck);
         exit(1);
      }
      if (length < 1)
      {
         printf("SelectiveRegression ERROR: each term should have >1 element.\n");
         exit(1);
      }
      coefTerms_[i] = new int[length+1];
      coefTerms_[i][0] = length;
      for (j = 1; j <= length; j++) coefTerms_[i][j] = 0;
      for (j = 0; j < length; j++) 
      {
         fscanf(fp, "%d", &coefTerms_[i][j+1]);
         if (coefTerms_[i][j+1] == 0 && length != 1)
         {
            printf("SelectiveRegression ERROR : constant term should have length=1\n");
            exit(1);
         }
         if (coefTerms_[i][j+1] < 0 || coefTerms_[i][j+1] > nInputs)
         {
            printf("SelectiveRegression ERROR : nTerms should be <= %d\n",
                   nInputs);
            exit(1);
         }
         coefTerms_[i][j+1]--;
      }
   }
   fscanf(fp, "%s", inStr);
   fclose(fp);
   if (strcmp(inStr, "PSUADE_END"))
   {
      printf("SelectiveRegression ERROR: wrong format in user regression file.\n");
      printf("The file should end with PSUADE_END\n");
      printf("Check your selective_regression_file file format as follow: \n");
      printf("line   1: PSUADE_BEGIN\n");
      printf("line   2: the number of terms.\n");
      printf("line   3: 1  number of inputs, input list (1-based).\n");
      printf("line   4: 2  number of inputs, input list (1-based).\n");
      printf(".........\n");
      printf("line   m: 1 0 (meaning constant term).\n");
      printf("lastline: END\n");
      exit(1);
   }
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
SelectiveRegression::~SelectiveRegression()
{
   int ii;
   if (regCoeffs_ != NULL) delete [] regCoeffs_;
   if (regStdevs_ != NULL) delete [] regStdevs_;
   if (coefTerms_ != NULL)
   {
      for (int i = 0; i < numTerms_; i++)
         if (coefTerms_[i] != NULL) delete [] coefTerms_[i];
      delete [] coefTerms_;
   }
   if (fuzzyC_ != NULL)
   {
      for (ii = 0; ii < numTerms_; ii++)
         if (fuzzyC_[ii] != NULL) delete [] fuzzyC_[ii];
      delete [] fuzzyC_;
   }
}

// ************************************************************************
// initialize
// ------------------------------------------------------------------------
int SelectiveRegression::initialize(double *X, double *Y)
{
   int    ii, status;
   double *XX;
   if (nSamples_ <= nInputs_)
   {
      printf("SelectiveRegression::initialize INFO- not enough points.\n");
      printf("                nSamples should be larger than nInputs.\n");
      return -1;
   }
   
   if (fuzzyC_ != NULL)
   {
      for (ii = 0; ii < numTerms_; ii++)
         if (fuzzyC_[ii] != NULL) delete [] fuzzyC_[ii];
      delete [] fuzzyC_;
      fuzzyC_ = NULL;
   }
   if (regCoeffs_ != NULL) delete [] regCoeffs_;
   if (regStdevs_ != NULL) delete [] regStdevs_;
   regCoeffs_ = NULL;
   regStdevs_ = NULL;

   printEquals(PL_INFO, 0);
   printf("* Selective Regression Analysis\n");
   printf("* You have the option to scale the sample matrix for\n");
   printf("* conditioning, but it may require more terms than\n");
   printf("* what you have provided.\n");
   printf("* To turn on scaling, first set rs_expert mode.\n");
   printEquals(PL_INFO, 0);
   XX = new double[nSamples_*nInputs_];
   initInputScaling(X, XX, 0);
   status = analyze(XX,Y);
   delete [] XX;
   if (status != 0)
   {
      printf("SelectiveRegression::initialize - ERROR detected.\n");
      return -1;
   }
   return 0;
}

// ************************************************************************
// Generate lattice data based on the input set
// ------------------------------------------------------------------------
int SelectiveRegression::genNDGridData(double *X, double *Y, int *NN,
                                       double **XX, double **YY)
{
   int totPts, ss;

   if (initialize(X,Y) != 0)
   {
      printf("SelectiveRegression::genNDGridData - ERROR detected.\n");
      return -1;
   }

   if ((*NN) == -999) return 0;

   genNDGrid(NN, XX);
   if ((*NN) == 0) return 0;
   totPts = (*NN);

   (*YY) = new double[totPts];
   for (ss = 0; ss < totPts; ss++)
      (*YY)[ss] = evaluatePoint(&((*XX)[ss*nInputs_]));

   return 0;
}

// ************************************************************************
// Generate 1D mesh results (set all other inputs to nominal values)
// ------------------------------------------------------------------------
int SelectiveRegression::gen1DGridData(double *X, double *Y, int ind1,
                                       double *settings, int *NN, 
                                       double **XX, double **YY)
{
   int    totPts, mm, nn;
   double HX, *Xloc;

   if (initialize(X,Y) != 0)
   {
      printf("SelectiveRegression::gen1DGridData - ERROR detected.\n");
      return -1;
   }

   totPts = nPtsPerDim_;
   HX = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 

   (*NN) = totPts;
   (*XX) = new double[totPts];
   (*YY) = new double[totPts];
   Xloc  = new double[nInputs_];
   for (nn = 0; nn < nInputs_; nn++) Xloc[nn] = settings[nn]; 
    
   for (mm = 0; mm < nPtsPerDim_; mm++) 
   {
      Xloc[ind1] = HX * mm + lowerBounds_[ind1];
      (*XX)[mm] = Xloc[ind1];
      (*YY)[mm] = evaluatePoint(Xloc);
   }

   delete [] Xloc;
   return 0;
}

// ************************************************************************
// Generate 2D mesh results (set all other inputs to nominal values)
// ------------------------------------------------------------------------
int SelectiveRegression::gen2DGridData(double *X, double *Y, int ind1,
                                       int ind2, double *settings, int *NN, 
                                       double **XX, double **YY)
{
   int    totPts, mm, nn, index;
   double *HX, *Xloc;

   if (initialize(X,Y) != 0)
   {
      printf("SelectiveRegression::gen2DGridData - ERROR detected.\n");
      return -1;
   }

   totPts = nPtsPerDim_ * nPtsPerDim_;
   HX    = new double[2];
   HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
   HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 

   (*NN) = totPts;
   (*XX) = new double[totPts * 2];
   (*YY) = new double[totPts];
   Xloc  = new double[nInputs_];
   for (nn = 0; nn < nInputs_; nn++) Xloc[nn] = settings[nn]; 
    
   for (mm = 0; mm < nPtsPerDim_; mm++) 
   {
      for (nn = 0; nn < nPtsPerDim_; nn++)
      {
         index = mm * nPtsPerDim_ + nn;
         Xloc[ind1] = HX[0] * mm + lowerBounds_[ind1];
         Xloc[ind2] = HX[1] * nn + lowerBounds_[ind2];
         (*XX)[index*2]   = Xloc[ind1];
         (*XX)[index*2+1] = Xloc[ind2];
         (*YY)[index] = evaluatePoint(Xloc);
      }
   }

   delete [] Xloc;
   delete [] HX;
   return 0;
}

// ************************************************************************
// Generate 3D mesh results (setting others to some nominal values) 
// ------------------------------------------------------------------------
int SelectiveRegression::gen3DGridData(double *X, double *Y, int ind1,
                                  int ind2, int ind3, double *settings, 
                                  int *NN, double **XX, double **YY)
{
   int    totPts, mm, nn, pp, index;
   double *HX, *Xloc;

   if (initialize(X,Y) != 0)
   {
      printf("SelectiveRegression::gen3DGridData - ERROR detected.\n");
      return -1;
   }

   totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
   HX    = new double[3];
   HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
   HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 
   HX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_ - 1); 

   (*NN) = totPts;
   (*XX) = new double[totPts * 3];
   (*YY) = new double[totPts];
   Xloc  = new double[nInputs_];
   for (nn = 0; nn < nInputs_; nn++) Xloc[nn] = settings[nn]; 
    
   for (mm = 0; mm < nPtsPerDim_; mm++) 
   {
      for (nn = 0; nn < nPtsPerDim_; nn++)
      {
         for (pp = 0; pp < nPtsPerDim_; pp++)
         {
            index = mm * nPtsPerDim_ * nPtsPerDim_ + nn * nPtsPerDim_ + pp;
            Xloc[ind1] = HX[0] * mm + lowerBounds_[ind1];
            Xloc[ind2] = HX[1] * nn + lowerBounds_[ind2];
            Xloc[ind3] = HX[2] * pp + lowerBounds_[ind3];
            (*XX)[index*3]   = Xloc[ind1];
            (*XX)[index*3+1] = Xloc[ind2];
            (*XX)[index*3+2] = Xloc[ind3];
            (*YY)[index] = evaluatePoint(Xloc);
         }
      }
   }

   delete [] Xloc;
   delete [] HX;
   return 0;
}

// ************************************************************************
// Generate 4D mesh results (setting others to some nominal values) 
// ------------------------------------------------------------------------
int SelectiveRegression::gen4DGridData(double *X, double *Y, int ind1, 
                            int ind2, int ind3, int ind4, double *settings, 
                            int *NN, double **XX, double **YY)
{
   int    totPts, mm, nn, pp, qq, index;
   double *HX, *Xloc;

   if (initialize(X,Y) != 0)
   {
      printf("SelectiveRegression::gen4DGridData - ERROR detected.\n");
      return -1;
   }

   totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
   HX    = new double[4];
   HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
   HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 
   HX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_ - 1); 
   HX[3] = (upperBounds_[ind4] - lowerBounds_[ind4]) / (nPtsPerDim_ - 1); 

   (*NN) = totPts;
   (*XX) = new double[totPts * 4];
   (*YY) = new double[totPts];
   Xloc  = new double[nInputs_];
   for (nn = 0; nn < nInputs_; nn++) Xloc[nn] = settings[nn]; 
    
   for (mm = 0; mm < nPtsPerDim_; mm++) 
   {
      for (nn = 0; nn < nPtsPerDim_; nn++)
      {
         for (pp = 0; pp < nPtsPerDim_; pp++)
         {
            for (qq = 0; qq < nPtsPerDim_; qq++)
            {
               index = mm*nPtsPerDim_*nPtsPerDim_*nPtsPerDim_ +
                       nn*nPtsPerDim_*nPtsPerDim_ + pp * nPtsPerDim_ + qq;
               Xloc[ind1] = HX[0] * mm + lowerBounds_[ind1];
               Xloc[ind2] = HX[1] * nn + lowerBounds_[ind2];
               Xloc[ind3] = HX[2] * pp + lowerBounds_[ind3];
               Xloc[ind4] = HX[3] * qq + lowerBounds_[ind4];
               (*XX)[index*4]   = Xloc[ind1];
               (*XX)[index*4+1] = Xloc[ind2];
               (*XX)[index*4+2] = Xloc[ind3];
               (*XX)[index*4+3] = Xloc[ind4];
               (*YY)[index] = evaluatePoint(Xloc);
            }
         }
      }
   }

   delete [] Xloc;
   delete [] HX;
   return 0;
}

// ************************************************************************
// Evaluate a given point
// ------------------------------------------------------------------------
double SelectiveRegression::evaluatePoint(double *X)
{
   int    ii, jj, kk;
   double Y, multiplier, x;

   if (regCoeffs_ == NULL) return 0.0;
   Y = 0.0;
   for (ii = 0; ii < numTerms_; ii++)
   {
      multiplier = 1.0;
      for (jj = 0; jj < coefTerms_[ii][0]; jj++)
      {
         kk = coefTerms_[ii][jj+1]; 
         if (kk < 0) x = 1.0;
         else        x = (X[kk] - XMeans_[kk]) / XStds_[kk];
         multiplier *= x;
      }
      Y += regCoeffs_[ii] * multiplier;
   }
   return Y;
}

// ************************************************************************
// Evaluate a number of points
// ------------------------------------------------------------------------
double SelectiveRegression::evaluatePoint(int npts, double *X, double *Y)
{
   for (int kk = 0; kk < npts; kk++)
      Y[kk] = evaluatePoint(&X[kk*nInputs_]);
   return 0.0;
}

// ************************************************************************
// Evaluate a given point and also its standard deviation
// ------------------------------------------------------------------------
double SelectiveRegression::evaluatePointFuzzy(double *X, double &std)
{
   int    iOne=1;
   double Y, ddata;
   evaluatePointFuzzy(iOne, X, &Y, &ddata);
   std = ddata;
   return Y;
}

// ************************************************************************
// Evaluate a number of points and also their standard deviations
// ------------------------------------------------------------------------
double SelectiveRegression::evaluatePointFuzzy(int npts, double *X, double *Y,
                                               double *Ystd)
{
   int    tt, kk, ii, ntimes=100;
   double *regStore, *Ys, *Xs, mean, std;

   if (regCoeffs_ == NULL)
   {
      printf("SelectiveRegression ERROR: initialize has not been called.\n");
      exit(1);
   }

   regStore = new double[numTerms_];
   for (kk = 0; kk < numTerms_; kk++) regStore[kk] = regCoeffs_[kk];
   Ys = new double[ntimes];
   for (kk = 0; kk < npts; kk++)
   {
      Xs = &(X[kk*nInputs_]);
      for (tt = 0; tt < ntimes; tt++)
      {
         for (ii = 0; ii < numTerms_; ii++) regCoeffs_[ii] = fuzzyC_[ii][tt];
         Ys[tt] = evaluatePoint(Xs);
      }
      mean = 0.0;
      for (tt = 0; tt < ntimes; tt++) mean += Ys[tt];
      mean /= (double) ntimes;
      std = 0.0;
      for (tt = 0; tt < ntimes; tt++) std += (Ys[tt]-mean)*(Ys[tt]-mean);
      std = sqrt(std/(ntimes-1));
      Y[kk] = mean;
      Ystd[kk] = std;
   }

   for (kk = 0; kk < numTerms_; kk++) regCoeffs_[kk] = regStore[kk];
   delete [] regStore;
   delete [] Ys;
   return 0.0;
}

// ************************************************************************
// perform mean/variance analysis
// ------------------------------------------------------------------------
int SelectiveRegression::analyze(double *X, double *Y)
{
   int    N, M, mm, nn, wlen, info, NRevised;
   double *B, *XX, SSresid, SStotal, R2, *XTX, var, *Bvar;
   double esum, ymax, *WW, *SS, *AA, *UU, *VV;
   char   jobu  = 'S', jobvt = 'S';
   char   pString[100];
   FILE   *fp;

   if (nSamples_ <= nInputs_) 
   {
      printf("SelectiveRegression::analyze ERROR - sample size too small.\n");
      return -1;
   }
   
   if (outputLevel_ >= 1)
   {
      printAsterisks(PL_INFO, 0);
      printf("*          Selective Regression Analysis\n");
      printf("* R-square gives a measure of the goodness of the model.\n");
      printf("* R-square should be close to 1 if it is a good model.\n");
      printEquals(PL_INFO, 0);
   }
   N = loadXMatrix(X, &XX);
   M = nSamples_;

   wlen = 5 * M;
   AA = new double[M*N];
   UU = new double[M*N];
   SS = new double[N];
   VV = new double[M*N];
   WW = new double[wlen];
   B  = new double[N];
   for (mm = 0; mm < M; mm++)
      for (nn = 0; nn < N; nn++)
         AA[mm+nn*M] = sqrt(weights_[mm]) * XX[mm+nn*M];

   if (outputLevel_ > 3) printf("Running SVD ...\n");
   dgesvd_(&jobu, &jobvt, &M, &N, AA, &M, SS, UU, &M, VV, &N, WW,
           &wlen, &info);
   if (outputLevel_ > 3) 
      printf("SVD completed: status = %d (should be 0).\n",info);

   if (info != 0)
   {
      printf("* SelectiveRegression Info: dgesvd returns a nonzero (%d).\n",
             info);
      printf("* SelectiveRegression terminates further processing.\n");
      delete [] XX;
      delete [] AA;
      delete [] UU;
      delete [] SS;
      delete [] VV;
      delete [] WW;
      delete [] B;
      return -1;
   }

   mm = 0;
   for (nn = 0; nn < N; nn++) if (SS[nn] < 0) mm++;
   if (mm > 0)
   {
      printf("* SelectiveRegression WARNING: some of the singular values\n");
      printf("*            are < 0. May spell trouble but will\n");
      printf("*            proceed anyway (%d).\n", mm);
   }
   if (SS[0] == 0.0) NRevised = 0;
   else
   {
      NRevised = N;
      for (nn = 1; nn < N; nn++)
         if (SS[nn-1] > 0 && SS[nn]/SS[nn-1] < 1.0e-8) NRevised--;
   }
   if (NRevised < N)
   {
      printf("* SelectiveRegression ERROR: sample rank = %d (need %d)\n",
             NRevised, N);
      delete [] XX;
      delete [] AA;
      delete [] UU;
      delete [] SS;
      delete [] VV;
      delete [] WW;
      delete [] B;
      return -1;
   }
   if (psMasterMode_ == 1)
   {
      printf("SelectiveRegression: matrix singular values \n");
      printf("The VERY small ones may cause poor numerical accuracy,\n");
      printf("but not keeping them may ruin the approximation power.\n");
      printf("So, select them judiciously.\n");
      for (nn = 0; nn < N; nn++)
         printf("Singular value %5d = %e\n", nn+1, SS[nn]);
      sprintf(pString, "How many to keep (1 - %d, 0 - all) ? ", N);
      NRevised = getInt(0,N,pString);
      if (NRevised == 0) NRevised = N;
      for (nn = NRevised; nn < N; nn++) SS[nn] = 0.0;
   }
   else
   {
      NRevised = N;
      for (nn = 1; nn < N; nn++)
      {
         if (SS[nn-1] > 0.0 && SS[nn]/SS[nn-1] < 1.0e-8)
         {
            SS[nn] = 0.0;
            NRevised--;
         }
      }
   }
   if (NRevised < N)
      printf("Number of singular values deleted = %d\n",N-NRevised);
   for (mm = 0; mm < NRevised; mm++)
   {
      WW[mm] = 0.0;
      for (nn = 0; nn < M; nn++)
         WW[mm] += UU[mm*M+nn] * sqrt(weights_[nn]) * Y[nn];
   }
   for (nn = 0; nn < NRevised; nn++) WW[nn] /= SS[nn];
   for (nn = NRevised; nn < N; nn++) WW[nn] = 0.0;
   for (mm = 0; mm < N; mm++)
   {
      B[mm] = 0.0;
      for (nn = 0; nn < NRevised; nn++) B[mm] += VV[mm*N+nn] * WW[nn];
   }
   delete [] AA;
   delete [] SS;
   delete [] UU;
   delete [] VV;

   if (psMasterMode_ == 1)
   {
      fp = fopen("selective_regression_error.m", "w");
      if(fp == NULL)
      {
         printf("fopen returned NULL in file %s line %d, exiting\n",
                 __FILE__, __LINE__);
         exit(1);
      }
      fprintf(fp, "%% This file contains errors of each data point.\n");
   }
   else fp = NULL;

   esum = ymax = 0.0;
   for (mm = 0; mm < nSamples_; mm++)
   {
      WW[mm] = 0.0;
      for (nn = 0; nn < N; nn++)
         WW[mm] = WW[mm] + XX[mm+nn*nSamples_] * B[nn];
      WW[mm] = WW[mm] - Y[mm];
      esum = esum + WW[mm] * WW[mm] * weights_[mm];
      if (fp != NULL) 
         fprintf(fp, "%6d %24.16e\n",mm+1,WW[mm]*sqrt(weights_[mm]));
      if (PABS(Y[mm]) > ymax) ymax = PABS(Y[mm]);
   }
   esum /= (double) nSamples_;
   esum = sqrt(esum);
   printf("* SelectiveRegression:: LS mean error = %11.4e (max=%11.4e)\n",
          esum, ymax); 

   if (fp != NULL)
   {
      fclose(fp);
      printf("FILE selective_regression_error.m contains data errors.\n");
   }

   computeSS(N, XX, Y, B, SSresid, SStotal);
   if (SStotal == 0) R2 = 1.0;
   else              R2  = 1.0 - SSresid / SStotal;
   if (nSamples_ > N) var = SSresid / (double) (nSamples_ - N);
   else               var = 0.0;
   if (var < 0)
   { 
      if (PABS(var) > 1.0e-12)
           printf("SelectiveRegression WARNING: var < 0.\n");
      else var = 0;
   }

   Bvar = new double[N];
   computeXTX(N, XX, &XTX);
   computeCoeffVariance(N, XTX, var, Bvar);
   regCoeffs_ = B;
   regStdevs_ = Bvar;

   PDFManager *pdfman = new PDFManager();
   int    cc, nTimes=100;
   int    *inPDFs = new int[numTerms_];
   double *inMeans = new double[numTerms_];
   double *inStds = new double[numTerms_];
   double *inUppers = new double[numTerms_];
   double *inLowers = new double[numTerms_];
   for (nn = 0; nn < numTerms_; nn++)
   {
      inPDFs[nn] = PSUADE_PDF_NORMAL;
      inMeans[nn] = regCoeffs_[nn];
      inStds[nn] = regStdevs_[nn];
      inUppers[nn] = inMeans[nn] + 4.0 * inStds[nn];
      inLowers[nn] = inMeans[nn] - 4.0 * inStds[nn];
      if (inUppers[nn] == inLowers[nn])
      {
         if (inUppers[nn] > 0) inUppers[nn] *= (1.0 + 1.0e-14);
         else                  inUppers[nn] *= (1.0 - 1.0e-14);
         if (inLowers[nn] > 0) inLowers[nn] *= (1.0 - 1.0e-14);
         else                  inLowers[nn] *= (1.0 + 1.0e-14);
         if (inUppers[nn] == 0.0) inUppers[nn] = 1e-14;
         if (inLowers[nn] == 0.0) inLowers[nn] = -1e-14;
      }
   }
   pdfman->initialize(numTerms_,inPDFs,inMeans,inStds,covMatrix_,NULL,NULL);
   psVector vLower, vUpper, vOut;
   vLower.load(numTerms_, inLowers);
   vUpper.load(numTerms_, inUppers);
   vOut.setLength(numTerms_*nTimes);
   pdfman->genSample(nTimes, vOut, vLower, vUpper);
   fuzzyC_ = new double*[numTerms_];
   for (nn = 0; nn < numTerms_; nn++)
   {
      fuzzyC_[nn] = new double[nTimes];
      for (cc = 0; cc < nTimes; cc++)
         fuzzyC_[nn][cc] = vOut[cc*numTerms_+nn];
   }
   delete pdfman;
   delete [] inPDFs;
   delete [] inStds;
   delete [] inMeans;
   delete [] inLowers;
   delete [] inUppers;

   if (outputLevel_ >= 0)
   {
      if (outputLevel_ >= 0) printRC(N, B, Bvar, XX, Y);
      printf("* SelectiveRegression model R-square = %12.4e\n",R2);
      printf("* adjusted   R-square = %12.4e\n",
             1.0 - (1.0 - R2) * ((M - 1) / (M - N - 1)));
      if (outputLevel_ > 1) printSRC(X, B, SStotal);
   }
   fp = NULL;
   if (psRSCodeGen_ == 1) fp = fopen("psuade_rs.info", "w");
   if (fp != NULL)
   {
      fprintf(fp,"/* ************************************************/\n");
      fprintf(fp,"/* Selective regression interpolator from PSUADE. */\n");
      fprintf(fp,"/* ===============================================*/\n");
      fprintf(fp,"/* This file contains information for interpolation\n");
      fprintf(fp,"   using response surface. Follow the steps below:\n");
      fprintf(fp,"   1. move this file to *.c file (e.g. main.c)\n");
      fprintf(fp,"   2. Compile main.c (cc -o main main.c -lm) \n");
      fprintf(fp,"   3. run: main input output\n");
      fprintf(fp,"          where input has the number of inputs and\n");
      fprintf(fp,"          the input values\n");
      fprintf(fp,"*/\n");
      fprintf(fp,"/* ===============================================*/\n");
      fprintf(fp,"#include <math.h>\n");
      fprintf(fp,"#include <stdlib.h>\n");
      fprintf(fp,"#include <stdio.h>\n");
      fprintf(fp,"int interpolate(int,double*,double*,double*);\n");
      fprintf(fp,"main(int argc, char **argv) {\n");
      fprintf(fp,"  int    i, iOne=1, nInps;\n");
      fprintf(fp,"  double X[%d], Y, S;\n",nInputs_);
      fprintf(fp,"  FILE   *fIn=NULL, *fOut=NULL;\n");
      fprintf(fp,"  if (argc != 3) {\n");
      fprintf(fp,"     printf(\"ERROR: not enough argument.\\n\");\n");
      fprintf(fp,"     exit(1);\n");
      fprintf(fp,"  }\n");
      fprintf(fp,"  fIn = fopen(argv[1], \"r\");\n");
      fprintf(fp,"  if (fIn == NULL) {\n");
      fprintf(fp,"     printf(\"ERROR: cannot open input file.\\n\");\n");
      fprintf(fp,"     exit(1);\n");
      fprintf(fp,"  }\n");
      fprintf(fp,"  fscanf(fIn, \"%%d\", &nInps);\n");
      fprintf(fp,"  if (nInps != %d) {\n", nInputs_);
      fprintf(fp,"    printf(\"ERROR - wrong nInputs.\\n\");\n");
      fprintf(fp,"    exit(1);\n");
      fprintf(fp,"  }\n");
      fprintf(fp,"  for (i=0; i<%d; i++) fscanf(fIn, \"%%lg\", &X[i]);\n",nInputs_);
      fprintf(fp,"  fclose(fIn);\n");
      fprintf(fp,"  interpolate(iOne, X, &Y, &S);\n");
      fprintf(fp,"  printf(\"Y = %%e\\n\", Y);\n");
      fprintf(fp,"  printf(\"S = %%e\\n\", S);\n");
      fprintf(fp,"  fOut = fopen(argv[2], \"w\");\n");
      fprintf(fp,"  if (fOut == NULL) {\n");
      fprintf(fp,"     printf(\"ERROR: cannot open output file.\\n\");\n");
      fprintf(fp,"     exit(1);\n");
      fprintf(fp,"  }\n");
      fprintf(fp,"  fprintf(fOut,\" %%e\\n\", Y);\n");
      fprintf(fp,"  fclose(fOut);\n");
      fprintf(fp,"}\n");
      fprintf(fp,"/* ===============================================*/\n");
      fprintf(fp,"/* Selective regression interpolation function    */\n");
      fprintf(fp,"/* X[0], X[1],   .. X[m-1]   - first point\n");
      fprintf(fp," * X[m], X[m+1], .. X[2*m-1] - second point\n");
      fprintf(fp," * ... */\n");
      fprintf(fp,"/* ===============================================*/\n");
      fprintf(fp,"static int\n");
      fprintf(fp,"CoefTerms[%d][%d] = \n", numTerms_, nInputs_+1);
      fprintf(fp,"{\n");
      for (mm = 0; mm < numTerms_; mm++)
      {
         fprintf(fp,"  { %d", coefTerms_[mm][0]);
         for (nn = 0; nn < coefTerms_[mm][0]; nn++)
            fprintf(fp," , %d", coefTerms_[mm][nn+1]);
         for (nn = coefTerms_[mm][0]; nn < nInputs_; nn++)
            fprintf(fp," , 0");
         fprintf(fp," },\n");
      }
      fprintf(fp,"};\n");
      fprintf(fp,"static double\n");
      fprintf(fp,"Coefs[%d][2] = \n", N);
      fprintf(fp,"{\n");
      for (mm = 0; mm < numTerms_; mm++)
         fprintf(fp,"  { %24.16e , %24.16e},\n",B[mm],Bvar[mm]);
      fprintf(fp,"};\n");
      fprintf(fp,"static double\n");
      fprintf(fp,"XStat[%d][2] = \n", nInputs_);
      fprintf(fp,"{\n");
      for (mm = 0; mm < nInputs_; mm++)
         fprintf(fp,"  { %24.16e , %24.16e},\n",XMeans_[mm],XStds_[mm]);
      fprintf(fp,"};\n");
      fprintf(fp,"void getCoefs(int kk, double *coefs);\n");
      fprintf(fp,"/* ===============================================*/\n");
      fprintf(fp,"int interpolate(int npts,double *X,double *Y,double *S){\n");
      fprintf(fp,"  int    ii,jj,kk,ll,mm,nterms=%d,nInps=%d;\n",
              numTerms_,nInputs_);
      fprintf(fp,"  int    ntimes=%d;\n",nTimes);
      fprintf(fp,"  double y, d, u, *yt, *x, *coefs, mean, std;\n");
      fprintf(fp,"  coefs = (double *) malloc(nterms *sizeof(double));\n");
      fprintf(fp,"  yt = (double *) malloc(ntimes * sizeof(double));\n");
      fprintf(fp,"  for (ii = 0; ii < npts; ii++) {\n");
      fprintf(fp,"    x = &(X[ii* nInps]);\n");
      fprintf(fp,"    for (kk = 0; kk < ntimes; kk++) {\n");
      fprintf(fp,"      getCoefs(kk, coefs);\n");
      fprintf(fp,"      y = 0.0;\n");
      fprintf(fp,"      for (jj = 0; jj < nterms; jj++) {\n");
      fprintf(fp,"        u = 1;\n");
      fprintf(fp,"        for (ll = 0; ll < CoefTerms[jj][0]; ll++){\n");
      fprintf(fp,"          mm = CoefTerms[jj][ll+1];\n");
      fprintf(fp,"          if (mm >= 0)\n");
      fprintf(fp,"            u = u * (x[mm]-XStat[mm][0])/XStat[mm][1];\n");
      fprintf(fp,"        }\n");
      fprintf(fp,"        y += coefs[jj] * u;\n");
      fprintf(fp,"      }\n");
      fprintf(fp,"      yt[kk] = y * %e + %e;\n", YStd_, YMean_);
      fprintf(fp,"    }\n");
      fprintf(fp,"    mean = 0.0;\n");
      fprintf(fp,"    for (kk = 0; kk < ntimes; kk++) mean += yt[kk];\n");
      fprintf(fp,"    mean /= ntimes;\n");
      fprintf(fp,"    Y[ii] = mean;\n");
      fprintf(fp,"    std = 0.0;\n");
      fprintf(fp,"    for (kk = 0; kk < ntimes; kk++) ");
      fprintf(fp," std += (yt[kk] - mean) * (yt[kk] - mean);\n");
      fprintf(fp,"    std = sqrt(std / (ntimes - 1.0));\n");
      fprintf(fp,"    S[ii] = std;\n");
      fprintf(fp,"  }\n");
      fprintf(fp,"  free(yt);\n");
      fprintf(fp,"  free(coefs);\n");
      fprintf(fp,"  return 0;\n");
      fprintf(fp,"}\n");
      fprintf(fp,"/* ===============================================*/\n");
      fprintf(fp,"static double\n");
      fprintf(fp,"CoefEnsemble[%d][%d] = \n", numTerms_, nTimes);
      fprintf(fp,"{\n");
      for (mm = 0; mm < numTerms_; mm++)
      {
         fprintf(fp," { %24.16e", fuzzyC_[mm][0]);
         for (nn = 1; nn < nTimes; nn++)
            fprintf(fp,", %24.16e", fuzzyC_[mm][nn]);
         fprintf(fp," },\n");
      }
      fprintf(fp," };\n");
      fprintf(fp,"void getCoefs(int kk, double *coefs) {\n");
      fprintf(fp,"  int mm;\n");
      fprintf(fp,"  for (mm = 0; mm < %d; mm++)\n",numTerms_);
      fprintf(fp,"    coefs[mm] = CoefEnsemble[mm][kk];\n");
      fprintf(fp,"}\n");
      fprintf(fp,"/* ===============================================*/\n");
      fclose(fp);
      printf("FILE psuade_rs.info contains the polynomial functional\n");
      printf("     form.\n");
   }
   fp = NULL;
   if (psRSCodeGen_ == 1) fp = fopen("psuade_rs.py", "w");
   if (fp != NULL)
   {
      fwriteRSPythonHeader(fp);
      fprintf(fp,"#==================================================\n");
      fprintf(fp,"# Selective regression interpolation\n");
      fprintf(fp,"#==================================================\n");
      fwriteRSPythonCommon(fp);
      fprintf(fp,"CoefTerms = [\n");
      for (mm = 0; mm < numTerms_; mm++)
      {
         fprintf(fp," [ %d", coefTerms_[mm][0]);
         for (nn = 0; nn < coefTerms_[mm][0]; nn++)
            fprintf(fp,", %d", coefTerms_[mm][nn+1]);
         for (nn = coefTerms_[mm][0]; nn < nInputs_; nn++)
            fprintf(fp,", 0");
         fprintf(fp," ],\n");
      }
      fprintf(fp,"]\n");
      fprintf(fp,"XStat = [\n");
      for (mm = 0; mm < nInputs_; mm++)
         fprintf(fp," [ %24.16e, %24.16e ],\n", XMeans_[mm], XStds_[mm]);
      fprintf(fp,"]\n");
      fprintf(fp,"CoefEnsemble = [\n");
      for (mm = 0; mm < numTerms_; mm++)
      {
         fprintf(fp," [ %24.16e", fuzzyC_[mm][0]);
         for (nn = 1; nn < nTimes; nn++)
            fprintf(fp,", %24.16e", fuzzyC_[mm][nn]);
         fprintf(fp," ],\n");
      }
      fprintf(fp,"]\n");
      fprintf(fp,"###################################################\n");
      fprintf(fp,"def getCoefs(kk) :\n");
      fprintf(fp,"  coefs = %d * [0.0]\n",numTerms_);
      fprintf(fp,"  for mm in range(%d) : \n", numTerms_);
      fprintf(fp,"    coefs[mm] = CoefEnsemble[mm][kk];\n");
      fprintf(fp,"  return coefs\n");
      fprintf(fp,"###################################################\n");
      fprintf(fp,"# Regression interpolation function  \n");
      fprintf(fp,"# X[0], X[1],   .. X[m-1]   - first point\n");
      fprintf(fp,"# X[m], X[m+1], .. X[2*m-1] - second point\n");
      fprintf(fp,"# ... \n");
      fprintf(fp,"#==================================================\n");
      fprintf(fp,"def interpolate(X): \n");
      fprintf(fp,"  nSamp = int(len(X) / %d)\n",nInputs_);
      fprintf(fp,"  Yt = %d * [0.0]\n", nTimes);
      fprintf(fp,"  Xt = %d * [0.0]\n", nInputs_);
      fprintf(fp,"  Ys = 2 * nSamp * [0.0]\n");
      fprintf(fp,"  for ss in range(nSamp) : \n");
      fprintf(fp,"    for ii in range(%d) : \n", nInputs_);
      fprintf(fp,"      Xt[ii] = X[ss*%d+ii]\n",nInputs_);
      fprintf(fp,"    for kk in range(%d) : \n", nTimes);
      fprintf(fp,"      coefs = getCoefs(kk)\n");
      fprintf(fp,"      Y = 0.0\n");
      fprintf(fp,"      for jj in range(%d) : \n", numTerms_);
      fprintf(fp,"        u = 1;\n");
      fprintf(fp,"        for ll in range(CoefTerms[jj][0]) : \n");
      fprintf(fp,"          mm = CoefTerms[jj][ll+1]\n");
      fprintf(fp,"          if (mm >= 0) :\n");
      fprintf(fp,"            u = u * (Xt[mm]-XStat[mm][0])/XStat[mm][1]\n");
      fprintf(fp,"        Y += coefs[jj] * u\n");
      fprintf(fp,"      Yt[kk] = Y * %e + %e;\n", YStd_, YMean_);
      fprintf(fp,"    mean = 0.0;\n");
      fprintf(fp,"    for kk in range(%d) : \n", nTimes);
      fprintf(fp,"      mean += Yt[kk]\n");
      fprintf(fp,"    mean /= %d\n", nTimes);
      fprintf(fp,"    std = 0.0;\n");
      fprintf(fp,"    for kk in range(%d) : \n", nTimes);
      fprintf(fp,"      std += (Yt[kk] - mean) * (Yt[kk] - mean)\n");
      fprintf(fp,"    std = math.sqrt(std / (%d - 1.0))\n",nTimes);
      fprintf(fp,"    Ys[2*ss] = mean\n");
      fprintf(fp,"    Ys[2*ss+1] = std\n");
      fprintf(fp,"  return Ys\n");
      fprintf(fp,"###################################################\n");
      fprintf(fp,"# main program\n");
      fprintf(fp,"#==================================================\n");
      fprintf(fp,"infileName  = sys.argv[1]\n");
      fprintf(fp,"outfileName = sys.argv[2]\n");
      fprintf(fp,"inputs = getInputData(infileName)\n");
      fprintf(fp,"outputs = interpolate(inputs)\n");
      fprintf(fp,"genOutputFile(outfileName, outputs)\n");
      fprintf(fp,"###################################################\n");
      printf("FILE psuade_rs.py contains the final selective polynomial\n");
      printf("     functional form.\n");
      fclose(fp);
   }
 
   delete [] XX;
   delete [] XTX;
   delete [] WW;
   return 0;
}

// *************************************************************************
// load the X matrix
// -------------------------------------------------------------------------
int SelectiveRegression::loadXMatrix(double *X, double **XXOut)
{
   int    M, N=0, m, n, k;
   double *XX=NULL, multiplier;

   M = nSamples_;
   N = numTerms_;
   XX = new double[M*N];
   for (m = 0; m < M; m++)
   {
      for (n = 0; n < N; n++)
      {
         multiplier = 1.0;
         for (k = 0; k < coefTerms_[n][0]; k++)
            if (coefTerms_[n][k+1] >= 0)
               multiplier *= X[m*nInputs_+coefTerms_[n][k+1]];
         XX[M*n+m] = multiplier;
      }
   }
   (*XXOut) = XX;
   return N;
}

// *************************************************************************
// form X^T X 
// -------------------------------------------------------------------------
int SelectiveRegression::computeXTX(int N, double *X, double **XXOut)
{
   int    nn, nn2, mm;
   double *XX, coef;

   XX = new double[nSamples_*N];
   for (nn = 0; nn < N; nn++)
   {
      for (nn2 = 0; nn2 < N; nn2++)
      {
         coef = 0.0;
         for (mm = 0; mm < nSamples_; mm++)
            coef += X[nn*nSamples_+mm] * weights_[mm] * X[nn2*nSamples_+mm];
         XX[nn*N+nn2] = coef;
      }
   }
   (*XXOut) = XX;
   return 0;
}

// *************************************************************************
// compute SS (sum of squares) statistics
// -------------------------------------------------------------------------
int SelectiveRegression::computeSS(int N, double *XX, double *Y,
                              double *B, double &SSresid, double &SStotal)
{
   int    nn, mm;
   double ymean, SSreg, SSresidCheck, ddata, rdata;
                                                                                
   SSresid = SStotal = ymean = SSreg = 0.0;
   for (mm = 0; mm < nSamples_; mm++) ymean += (sqrt(weights_[mm]) * Y[mm]);
   ymean /= (double) nSamples_;
   for (mm = 0; mm < nSamples_; mm++)
   {
      ddata = 0.0;
      for (nn = 0; nn < N; nn++) ddata += (XX[mm+nn*nSamples_] * B[nn]);
      rdata = Y[mm] - ddata;
      SSresid += rdata * Y[mm] * weights_[mm];
      SSresidCheck += rdata * rdata * weights_[mm];
      SSreg += (ddata - ymean) * (ddata - ymean);
   }
   for (mm = 0; mm < nSamples_; mm++)
      SStotal += weights_[mm] * (Y[mm] - ymean) * (Y[mm] - ymean);
   if (outputLevel_ > 0)
   {
      printf("* SelectiveRegression: SStot  = %24.16e\n", SStotal);
      printf("* SelectiveRegression: SSreg  = %24.16e\n", SSreg);
      printf("* SelectiveRegression: SSres  = %24.16e\n", SSresid);
      printf("* SelectiveRegression: SSres  = %24.16e (true)\n", SSresidCheck);
   }
   SSresid = SSresidCheck;
   if (outputLevel_ > 0 && nSamples_ != N)
   {
      printf("* SelectiveRegression: eps(Y) = %24.16e\n",
             SSresidCheck/(nSamples_-N));
   }
   return 0;
}

// *************************************************************************
// compute coefficient variance ((diagonal of sigma^2 (X' X)^(-1))
// -------------------------------------------------------------------------
int SelectiveRegression::computeCoeffVariance(int N,double *XX,double var,
                                              double *B)
{
   int    ii, jj, lwork, iOne=1, info, errCnt=0;
   double *B2, *work, *XT;
   char   trans[1];
   FILE   *fp;

   (*trans) = 'N';
   B2 = new double[N];
   XT = new double[N*N];
   lwork = 2 * N * N;
   work  = new double[lwork];
   for (ii = 0; ii < N; ii++)
   {
      for (jj = 0; jj < N*N; jj++) XT[jj] = XX[jj];
      for (jj = 0; jj < N; jj++) B2[jj] = 0.0;
      B2[ii] = var;
      dgels_(trans, &N, &N, &iOne, XT, &N, B2, &N, work, &lwork, &info);
      if (info != 0)
         printf("SelectiveRegression WARNING: dgels returns error %d.\n",info);
      if (B2[ii] < 0) errCnt++;
      if (B2[ii] < 0) B[ii] = sqrt(-B2[ii]);
      else            B[ii] = sqrt(B2[ii]);
   }
   if (errCnt > 0)
   {
      printf("* SelectiveRegression WARNING: some of the coefficient\n");
      printf("*            variances are < 0. May spell trouble but\n");
      printf("*            will proceed anyway (%d).\n",errCnt);
   }
   delete [] B2;
   delete [] XT;

   int    *ipiv = new int[N+1];
   double *invA = new double[lwork];
   double ddata, ddata2;
   for (ii = 0; ii < N*N; ii++) invA[ii] = XX[ii];
   dgetrf_(&N, &N, invA, &N, ipiv, &info);
   if (info != 0)
      printf("LegendreRegression WARNING: dgels returns error %d.\n",info);
   dgetri_(&N, invA, &N, ipiv, work, &lwork, &info);
   covMatrix_.setDim(N,N);
   for (ii = 0; ii < N; ii++)
   {
      for (jj = 0; jj < N; jj++)
      {
         ddata = invA[ii*N+jj] * var;
         covMatrix_.setEntry(ii,jj,ddata);
      }
   }
   for (ii = 0; ii < N; ii++)
   {
      ddata = covMatrix_.getEntry(ii,ii);
      ddata = sqrt(ddata);
      for (jj = 0; jj < N; jj++)
      {
         if (ii != jj)
         {
            ddata2 = covMatrix_.getEntry(ii,jj);
            if (ddata != 0) ddata2 /= ddata;
            covMatrix_.setEntry(ii,jj,ddata2);
         }
      }
   }
   for (jj = 0; jj < N; jj++)
   {
      ddata = covMatrix_.getEntry(jj,jj);
      ddata = sqrt(ddata);
      for (ii = 0; ii < N; ii++)
      {
         if (ii != jj)
         {
            ddata2 = covMatrix_.getEntry(ii,jj);
            if (ddata != 0) ddata2 /= ddata;
            covMatrix_.setEntry(ii,jj,ddata2);
         }
      }
   }
   ddata = 1.0;
   for (ii = 0; ii < N; ii++) covMatrix_.setEntry(ii,ii,ddata);
   for (ii = 0; ii < N; ii++)
   {
      for (jj = 0; jj < ii; jj++)
      {
         ddata  = covMatrix_.getEntry(ii,jj);
         ddata2 = covMatrix_.getEntry(jj,ii);
         ddata  = 0.5 * (ddata + ddata2);
         covMatrix_.setEntry(ii,jj,ddata);
         covMatrix_.setEntry(jj,ii,ddata);
      }
   }
   if (psMasterMode_ == 1)
   {
      fp = fopen("user_regression_correlation_matrix","w");
      for (ii = 0; ii < N; ii++)
      {
          for (jj = 0; jj < N; jj++)
             fprintf(fp, "%e ", covMatrix_.getEntry(ii,jj));
          fprintf(fp, "\n");
      }
      fclose(fp);
   }
   errCnt = 0;
   for (ii = 0; ii < N; ii++)
   {
      for (jj = 0; jj < N; jj++)
      {
         ddata = covMatrix_.getEntry(ii,jj);
         if (ii != jj && (ddata >=1 || ddata <= -1))
         {
            errCnt++;
            covMatrix_.setEntry(ii,jj,0.0);
         }
      }
   }
   char inStr[1001];
   if (errCnt > 0)
   {
      printf("UserRegression WARNING:\n");
      printf("  Correlation matrix has invalid entries (%d out of %d).\n",
             errCnt, N*(N-1));
      printf("  EVALUATION MAY BE INCORRECT.\n");
      printf("  CONTINUE ANYWAY (will set them to zeros)? (y or n)");
      scanf("%s", inStr);
      if (inStr[0] != 'y') exit(1);
      fgets(inStr, 100, stdin);
   }
   delete [] work;
   delete [] ipiv;
   delete [] invA;
   return info;
}

// *************************************************************************
// print statistics
// -------------------------------------------------------------------------
int SelectiveRegression::printRC(int N,double *B,double *Bvar,double *XX,
                                 double *Y)
{
   int    ii, jj, maxTerms;
   double coef;
   char   fname[200];
   FILE   *fp;

   maxTerms = 0;
   for (jj = 0; jj < numTerms_; jj++) 
      if (coefTerms_[jj][0] > maxTerms) maxTerms = coefTerms_[jj][0];

   printf("*----------------------------------------------------------*\n");
   printf("*      ");
   for (jj = 0; jj < maxTerms; jj++) printf("    ");
   printf("   coefficient   std. error   t-value\n");

   for (ii = 0; ii < numTerms_; ii++)
   {
      if (PABS(Bvar[ii]) < 1.0e-15) coef = 0.0;
      else                          coef = B[ii] / Bvar[ii]; 
      if (PABS(coef) > 0.0)
      {
         printf("* Input");
         for (jj = 0; jj < coefTerms_[ii][0]; jj++)
            printf(" %2d ", coefTerms_[ii][jj+1]+1);
         for (jj = coefTerms_[ii][0]; jj < maxTerms; jj++) printf("    ");
         printf("= %12.4e %12.4e %12.4e\n", B[ii], Bvar[ii], coef);
      }
   }
   strcpy(fname, "dataVariance");
   printDashes(PL_INFO, 0);

   if (psMasterMode_ == 1)
   {
      fp = fopen(fname, "w");
      if(fp == NULL)
      {
         printf("fopen returned NULL in file %s line %d, exiting\n", 
                __FILE__, __LINE__);
         exit(1);
      }
      fprintf(fp, "data number     data         standard error\n");
      for (ii = 0; ii < nSamples_; ii++)
      {
         coef = 0.0;
         for (jj = 0; jj < N; jj++) coef += XX[jj*nSamples_+ii] * Bvar[jj];
         if (coef < 0)
            fprintf(fp,"%7d         %12.4e %12.4e\n",ii+1,Y[ii],-sqrt(-coef));
         else
            fprintf(fp,"%7d         %12.4e %12.4e\n",ii+1,Y[ii],sqrt(coef));
      }
      fclose(fp);
   }
   return 0;
}

// *************************************************************************
// print standardized regression coefficients
// -------------------------------------------------------------------------
int SelectiveRegression::printSRC(double *X, double *B, double SStotal)
{
   int    nn, mm, ii, index, length, maxTerms;
   double denom, xmean, coef, Bmax, coef1, *B2;

   printEquals(PL_INFO, 0);
   printf("* Standardized Regression Coefficients (SRC)\n");
   printf("* When R-square is acceptable (order assumption holds), the*\n");
   printf("* absolute values of SRCs provide variable importance.     *\n"); 
   printDashes(PL_INFO, 0);
   printf("* based on nSamples = %d\n", nSamples_);

   maxTerms = 0;
   for (ii = 0; ii < numTerms_; ii++) 
      if (coefTerms_[ii][0] > maxTerms) maxTerms = coefTerms_[ii][0];
   B2 = new double[nSamples_];
   denom = sqrt(SStotal / (double) (nSamples_ - 1));
   Bmax  = 0.0;
   for (nn = 0; nn < numTerms_; nn++)
   {
      coef = 1.0;
      length = coefTerms_[nn][0];
      for (ii = 0; ii < length; ii++)
      {
         xmean = 1.0;
         coef1 = 0.0;
         index = coefTerms_[nn][ii+1];
         if (index >= 0)
         {
            for (mm = 0; mm < nSamples_; mm++) xmean += X[mm*nInputs_+index];
            xmean /= (double) nSamples_;
            for (mm = 0; mm < nSamples_; mm++)
               coef1 += (X[mm*nInputs_+index]-xmean) * (X[mm*nInputs_+index]-xmean);
            coef1 = sqrt(coef1 / (double) (nSamples_ - 1));
            coef *= coef1;
         }
      }
      B2[nn] = B[nn] * coef / denom;
      if (PABS(B2[nn]) > Bmax) Bmax = PABS(B2[nn]);
   }
   for (nn = 0; nn < numTerms_; nn++)
   {
      if (PABS(B2[nn]) > 1.0e-12 * Bmax)
      {
         printf("* Input");
         for (ii = 0; ii < coefTerms_[nn][0]; ii++)
            printf(" %2d ",coefTerms_[nn][ii+1]+1);
         for (ii = coefTerms_[nn][0]; ii < maxTerms; ii++) printf("    ");
         printf("= %12.4e\n",B2[nn]);
      }
   }
   delete [] B2;
   printAsterisks(PL_INFO, 0);
   return 0;
}

