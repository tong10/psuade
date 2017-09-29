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
#include "Main/Psuade.h"
#include "Util/sysdef.h"
#include "Util/PsuadeUtil.h"
#include "PDFLib/PDFBase.h"
#include "PDFLib/PDFNormal.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

extern "C" {
   void dgels_(char *, int *, int *, int *, double *A, int *LDA,
               double *B, int *LDB, double *WORK, int *LWORK, int *INFO);
   void dgesvd_(char *, char *, int *, int *, double *, int *, double *,
               double *, int *, double *, int *, double *, int *, int *);
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
      printf("  file format : first line - numTerms.\n");
      printf("                2nd line on - term, length, indices(1-based).\n");
      exit(1);
   }
   else
   {
      printf("SelectiveRegression INFO: using selective_regression_file.\n");
   }
   fscanf(fp, "%s", inStr);
   if (strcmp(inStr, "BEGIN"))
   {
      printf("SelectiveRegression: wrong format in selective_regression_file.\n");
      printf("First  line: BEGIN\n");
      printf("Second line: the number of terms.\n");
      printf("Third  line: 1  number of inputs, input list (1-based).\n");
      printf("Fourth line: 2  number of inputs, input list (1-based).\n");
      printf(".........\n");
      printf("Last   line : END\n");
      fclose(fp);
      exit(1);
   }
   fscanf(fp, "%d", &numTerms_);
   if (numTerms_ >= nSamples)
   {
      printf("SelectiveRegression ERROR: no. of terms should be < nSamples.\n");
      exit(1);
   }
   regCoeffs_ = new double[numTerms_+1];
   for (i = 0; i <= numTerms_; i++) regCoeffs_[i] = 0.0;
   regStdevs_ = new double[numTerms_+1];
   for (i = 0; i <= numTerms_; i++) regStdevs_[i] = 0.0;
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
         if (coefTerms_[i][j+1] < 1 || coefTerms_[i][j+1] > nInputs)
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
   if (strcmp(inStr, "END"))
   {
      printf("SelectiveRegression: wrong format in user regression file.\n");
      printf("The file should end with END\n");
      exit(1);
   }
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
SelectiveRegression::~SelectiveRegression()
{
   if (regCoeffs_ != NULL) delete [] regCoeffs_;
   if (regStdevs_ != NULL) delete [] regStdevs_;
   if (coefTerms_ != NULL)
   {
      for (int i = 0; i < numTerms_; i++)
         if (coefTerms_[i] != NULL) delete [] coefTerms_[i];
      delete [] coefTerms_;
   }
}

// ************************************************************************
// Generate lattice data based on the input set
// ------------------------------------------------------------------------
int SelectiveRegression::genNDGridData(double *X, double *Y, int *N2,
                                       double **X2, double **Y2)
{
   int     totPts, sampleID, inputID, status;
   double  *HX, *Xloc;

   if (nInputs_ <= 0 || nSamples_ <= 0)
   {
      printf("SelectiveRegression::genNDGridData ERROR- invalid argument.\n");
      exit(1);
   } 
   if (nSamples_ <= nInputs_)
   {
      printf("SelectiveRegression::genNDGridData INFO- not enough points.\n");
      printf("                nSamples should be larger than nInputs.\n");
      return -1;
   }
   
   status = analyze(X, Y);
   if (status != 0)
   {
      printf("SelectiveRegression::genNDGridData - ERROR detected.\n");
      return -1;
   }

   if ((*N2) == -999) return 0;
   if (nInputs_ > 21)
   {
      printf("SelectiveRegression::genNDGridData INFO - nInputs > 21.\n");
      printf("                No lattice points generated.\n");
      (*N2) = 0;
      return -1;
   }

   if (nInputs_ == 21 && nPtsPerDim_ >    2) nPtsPerDim_ =  2;
   if (nInputs_ == 20 && nPtsPerDim_ >    2) nPtsPerDim_ =  2;
   if (nInputs_ == 19 && nPtsPerDim_ >    2) nPtsPerDim_ =  2;
   if (nInputs_ == 18 && nPtsPerDim_ >    2) nPtsPerDim_ =  2;
   if (nInputs_ == 17 && nPtsPerDim_ >    2) nPtsPerDim_ =  2;
   if (nInputs_ == 16 && nPtsPerDim_ >    2) nPtsPerDim_ =  2;
   if (nInputs_ == 15 && nPtsPerDim_ >    2) nPtsPerDim_ =  2;
   if (nInputs_ == 14 && nPtsPerDim_ >    2) nPtsPerDim_ =  2;
   if (nInputs_ == 13 && nPtsPerDim_ >    3) nPtsPerDim_ =  3;
   if (nInputs_ == 12 && nPtsPerDim_ >    3) nPtsPerDim_ =  3;
   if (nInputs_ == 11 && nPtsPerDim_ >    3) nPtsPerDim_ =  3;
   if (nInputs_ == 10 && nPtsPerDim_ >    4) nPtsPerDim_ =  4;
   if (nInputs_ ==  9 && nPtsPerDim_ >    5) nPtsPerDim_ =  5;
   if (nInputs_ ==  8 && nPtsPerDim_ >    6) nPtsPerDim_ =  6;
   if (nInputs_ ==  7 && nPtsPerDim_ >    8) nPtsPerDim_ =  8;
   if (nInputs_ ==  6 && nPtsPerDim_ >   10) nPtsPerDim_ = 10;
   if (nInputs_ ==  5 && nPtsPerDim_ >   16) nPtsPerDim_ = 16;
   if (nInputs_ ==  4 && nPtsPerDim_ >   32) nPtsPerDim_ = 32;
   if (nInputs_ ==  3 && nPtsPerDim_ >   64) nPtsPerDim_ = 64;
   if (nInputs_ ==  2 && nPtsPerDim_ > 1024) nPtsPerDim_ = 1024;
   if (nInputs_ ==  1 && nPtsPerDim_ > 8192) nPtsPerDim_ = 8192;
   totPts = nPtsPerDim_;
   for (inputID = 1; inputID < nInputs_; inputID++)
      totPts = totPts * nPtsPerDim_;
   HX = new double[nInputs_];
   for (inputID = 0; inputID < nInputs_; inputID++)
      HX[inputID] = (upperBounds_[inputID] - lowerBounds_[inputID]) /
                    (double) (nPtsPerDim_ - 1);

   (*X2) = new double[nInputs_ * totPts];
   (*Y2) = new double[totPts];
   (*N2) = totPts;
   Xloc  = new double[nInputs_];

   for (inputID = 0; inputID < nInputs_; inputID++)
      Xloc[inputID] = lowerBounds_[inputID];

   for (sampleID = 0; sampleID < totPts; sampleID++)
   {
      for (inputID = 0; inputID < nInputs_; inputID++ )
         (*X2)[sampleID*nInputs_+inputID] = Xloc[inputID];
      for (inputID = 0; inputID < nInputs_; inputID++ )
      {
         Xloc[inputID] += HX[inputID];
         if (Xloc[inputID] < upperBounds_[inputID] ||
              PABS(Xloc[inputID] - upperBounds_[inputID]) < 1.0E-7) break;
         else Xloc[inputID] = lowerBounds_[inputID];
      }
      (*Y2)[sampleID] = evaluatePoint(&((*X2)[sampleID*nInputs_]));
   }

   delete [] Xloc;
   delete [] HX;
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

   (*NN) = -999;
   genNDGridData(X, Y, NN, NULL, NULL);

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

   (*NN) = -999;
   genNDGridData(X, Y, NN, NULL, NULL);

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

   (*NN) = -999;
   genNDGridData(X, Y, NN, NULL, NULL);

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

   (*NN) = -999;
   genNDGridData(X, Y, NN, NULL, NULL);

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
   int    i, j;
   double Y, multiplier;

   if (regCoeffs_ == NULL) return 0.0;
   Y = regCoeffs_[numTerms_];
   for (i = 0; i < numTerms_; i++)
   {
      multiplier = 1.0;
      for (j = 0; j < coefTerms_[i][0]; j++)
         multiplier *= X[coefTerms_[i][j+1]]; 
      Y += regCoeffs_[i] * multiplier;
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
   int    mm, nn, cc, nTimes=100;
   double dtmp, *Ys, C1, C2, mean, stds, *uppers, *lowers, multiplier;
   PDFBase **PDFPtrs;

   if (regCoeffs_ == NULL) return 0.0;

   PDFPtrs = new PDFBase*[numTerms_+1];
   uppers = new double[numTerms_+1];
   lowers = new double[numTerms_+1];
   for (mm = 0; mm <= numTerms_; mm++) 
   {
      mean = regCoeffs_[mm];
      stds = regStdevs_[mm];
      uppers[mm] = mean + 2.0 * std;
      lowers[mm] = mean - 2.0 * std;
      PDFPtrs[mm] = (PDFBase *) new PDFNormal(mean, std);
   }
   Ys = new double[nTimes];

   mean = 0.0;
   for (cc = 0; cc < nTimes; cc++)
   {
      C1 = PSUADE_drand();
      PDFPtrs[numTerms_]->invCDF(1, &C1, &C2, lowers[numTerms_], 
                                  uppers[numTerms_]);
      dtmp = C2;
      for (mm = 0; mm < numTerms_; mm++)
      {
         multiplier = 1.0;
         for (nn = 0; nn < coefTerms_[mm][0]; nn++)
            multiplier *= X[coefTerms_[mm][nn+1]]; 
         C1 = PSUADE_drand();
         PDFPtrs[mm]->invCDF(1, &C1, &C2, lowers[mm], uppers[mm]);
         dtmp += C2 * multiplier;
      }
      Ys[cc] = dtmp;
      mean += dtmp;
   }
   mean /= (double) nTimes;
   stds = 0.0;
   for (cc = 0; cc < nTimes; cc++)
      stds += (Ys[cc] - mean) * (Ys[cc] - mean);
   stds = sqrt(stds / (nTimes - 1));

   for (mm = 0; mm <= numTerms_; mm++) delete [] PDFPtrs[mm];
   delete [] PDFPtrs;
   delete [] uppers;
   delete [] lowers;
   delete [] Ys;
   std = stds; 
   return mean;
}

// ************************************************************************
// Evaluate a number of points and also their standard deviations
// ------------------------------------------------------------------------
double SelectiveRegression::evaluatePointFuzzy(int npts, double *X, double *Y,
                                          double *Ystd)
{
   for (int kk = 0; kk < npts; kk++)
      Y[kk] = evaluatePointFuzzy(&X[kk*nInputs_], Ystd[kk]);
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
   char   jobu  = 'A', jobvt = 'A';
   char   pString[100];
   FILE   *fp;

   if (nInputs_ <= 0 || nSamples_ <= 0)
   {
      printf("SelectiveRegression::analyze ERROR - invalid arguments.\n");
      exit(1);
   } 
   if (nSamples_ <= nInputs_) 
   {
      printf("SelectiveRegression::analyze ERROR - sample size too small.\n");
      return -1;
   }
   
   if (outputLevel_ >= 1)
   {
      printAsterisks(0);
      printf("*          Selective Regression Analysis\n");
      printf("* R-square gives a measure of the goodness of the model.\n");
      printf("* R-square should be close to 1 if it is a good model.\n");
      printEquals(0);
   }
   N = loadXMatrix(X, &XX);
   M = nSamples_;

   wlen = 5 * M;
   AA = new double[M*N];
   UU = new double[M*M];
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
   if (outputLevel_ > 3) printf("SVD completed: status = %d (should be 0).\n",info);

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

   if (SS[0] == 0.0) NRevised = 0;
   else
   {
      NRevised = N;
      for (nn = 1; nn < N; nn++)
         if (SS[nn-1] > 0 && SS[nn]/SS[nn-1] < 1.0e-8) NRevised--;
   }
   if (NRevised < N)
   {
      printf("* SelectiveRegression ERROR: sample true rank = %d (need %d)\n",
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
   if (psRSExpertMode_ == 1)
   {
      printf("SelectiveRegression: singular values for the Vandermonde matrix\n");
      printf("The VERY small ones may cause poor numerical accuracy,\n");
      printf("but not keeping them may ruin the approximation power.\n");
      printf("So, select them judiciously.\n");
      for (nn = 0; nn < N; nn++)
         printf("Singular value %5d = %e\n", nn+1, SS[nn]);
      sprintf(pString, "How many to keep (1 - %d) ? ", N);
      NRevised = getInt(1,N,pString);
      for (nn = NRevised; nn < N; nn++) SS[nn] = 0.0;
   }
   else
   {
      NRevised = N;
      for (nn = 1; nn < N; nn++)
      {
         if (SS[nn-1] > 0.0 && SS[nn]/SS[nn-1] < 1.0e-8) SS[nn] = 0.0;
         else NRevised++;
      }
   }
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

   if (outputLevel_ > 5)
   {
      fp = fopen("selective_regression_error.m", "w");
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
      if (fp != NULL) fprintf(fp, "%6d %24.16e\n",mm+1,WW[mm]*sqrt(weights_[mm]));
      if (PABS(Y[mm]) > ymax) ymax = PABS(Y[mm]);
   }
   esum /= (double) nSamples_;
   esum = sqrt(esum);
   printf("* SelectiveRegression:: LS mean error = %11.4e (max=%11.4e)\n",
          esum, ymax); 

   if (fp != NULL)
   {
      fclose(fp);
      printf("FILE selective_regression_error contains data errors.\n");
   }

   computeSS(N, XX, Y, B, SSresid, SStotal);
   R2  = (SStotal - SSresid) / SStotal;
   if (nSamples_ > N) var = SSresid / (double) (nSamples_ - N);
   else               var = 0.0;

   Bvar = new double[N];
   computeXTX(N, XX, &XTX);
   computeCoeffVariance(N, XTX, var, Bvar);

   if (outputLevel_ >= 0)
   {
      if (outputLevel_ >= 0) printRC(N, B, Bvar, XX, Y);
      printf("* SelectiveRegression model R-square = %12.4e\n",R2);
      printf("* adjusted   R-square = %12.4e\n",
             1.0 - (1.0 - R2) * ((M - 1) / (M - N - 1)));
      if (outputLevel_ > 1) printSRC(X, B, SStotal);
   }
 
   if (regCoeffs_ != NULL) delete [] regCoeffs_;
   regCoeffs_ = B;
   delete [] XX;
   delete [] XTX;
   delete [] Bvar;
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
   N = numTerms_ + 1;
   XX = new double[M*N];
   for (m = 0; m < M; m++)
   {
      for (n = 0; n < N-1; n++)
      {
         multiplier = 1.0;
         for (k = 0; k < coefTerms_[n][0]; k++)
            multiplier *= X[m*nInputs_+coefTerms_[n][k+1]];
         XX[M*n+m] = multiplier;
      }
      XX[M*(N-1)+m] = 1.0;
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
   double *R, ymean;
                                                                                
   SSresid = SStotal = ymean = 0.0;
   R = new double[nSamples_];
   for (mm = 0; mm < nSamples_; mm++)
   {
      R[mm] = Y[mm];
      for (nn = 0; nn < N; nn++) R[mm] -= (XX[mm+nn*nSamples_] * B[nn]);
      SSresid += R[mm] * Y[mm] * weights_[mm];
      ymean += (sqrt(weights_[mm]) * Y[mm]);
   }
   ymean /= (double) nSamples_;
   SStotal = - ymean * ymean * (double) nSamples_;
   for (mm = 0; mm < nSamples_; mm++)
      SStotal += weights_[mm] * Y[mm] * Y[mm];
   delete [] R;
   return 0;
}

// *************************************************************************
// compute coefficient variance ((diagonal of sigma^2 (X' X)^(-1))
// -------------------------------------------------------------------------
int SelectiveRegression::computeCoeffVariance(int N,double *XX,double var,
                                              double *B)
{
   int    ii, jj, lwork, iOne=1, info;
   double *B2, *work, *XT;
   char   trans[1];

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
      if (B2[ii] < 0) B[ii] = -sqrt(-B2[ii]);
      else            B[ii] = sqrt(B2[ii]);
   }
   delete [] B2;
   delete [] work;
   delete [] XT;
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

   if (PABS(Bvar[numTerms_]) < 1.0e-15) coef = 0.0;
   else                                 coef = B[numTerms_] / Bvar[numTerms_]; 
   if (PABS(coef) > 0.0)
   {
      printf("* Constant ");
      for (jj = 1; jj < maxTerms; jj++) printf("    ");
      printf("= %12.4e %12.4e %12.4e\n", B[numTerms_], Bvar[numTerms_], coef);
   }
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
   printDashes(0);

   if (outputLevel_ >= 3)
   {
      fp = fopen(fname, "w");
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

   printEquals(0);
   printf("* Standardized Regression Coefficients (SRC)\n");
   printf("* When R-square is acceptable (order assumption holds), the*\n");
   printf("* absolute values of SRCs provide variable importance.     *\n"); 
   printDashes(0);
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
         index = coefTerms_[nn][ii+1];
         xmean = 0.0;
         for (mm = 0; mm < nSamples_; mm++) xmean += X[mm*nInputs_+index];
         xmean /= (double) nSamples_;
         coef1 = 0.0;
         for (mm = 0; mm < nSamples_; mm++)
            coef1 += (X[mm*nInputs_+index]-xmean) * (X[mm*nInputs_+index]-xmean);
         coef1 = sqrt(coef1 / (double) (nSamples_ - 1));
         coef *= coef1;
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
   printAsterisks(0);
   return 0;
}

