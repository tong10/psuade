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
// Functions for the class Regression
// AUTHOR : CHARLES TONG
// DATE   : 2005
// ************************************************************************
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "Util/sysdef.h"
#include "Main/Psuade.h"
#include "Regression.h"
#include "PDFLib/PDFBase.h"
#include "PDFLib/PDFNormal.h"
#include "Util/PsuadeUtil.h"

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
Regression::Regression(int nInputs,int nSamples):FuncApprox(nInputs,nSamples)
{
   faID_ = PSUADE_RS_REGR4;
   order_     = 4;
   numTerms_  = 0;
   regCoeffs_ = NULL;
   regStdevs_ = NULL;
   fuzzyC_    = NULL;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
Regression::~Regression()
{
   if (regCoeffs_ != NULL) delete [] regCoeffs_;
   if (regStdevs_ != NULL) delete [] regStdevs_;
   if (fuzzyC_ != NULL)
   {
      for (int ii = 0; ii < numTerms_; ii++) delete [] fuzzyC_[ii];
      delete [] fuzzyC_;
   }
}

// ************************************************************************
// Generate lattice data based on the input set
// ------------------------------------------------------------------------
int Regression::genNDGridData(double *X, double *Y, int *N2, double **X2, 
                              double **Y2)
{
   int     totPts, mm, ii, status;
   double  *HX, *Xloc;

   if (nInputs_ <= 0 || nSamples_ <= 0)
   {
      printf("Regression ERROR: parameters not selected.\n");
      printf("    nInputs  = %6d (should be > 0).\n", nInputs_);
      printf("    nSamples = %6d (should be > 0).\n", nSamples_);
      return 0;
   } 
   if (nSamples_ <= nInputs_)
   {
      printf("Regression ERROR: no regression, not enough points.\n");
      printf("   First order regression needs >= nInputs+1 sample points.\n");
      return 0;
   }
   
   status = analyze(X, Y);
   if (status != 0)
   {
      printf("Regression: ERROR detected in regression analysis.\n");
      (*N2) = 0;
      return -1;
   }

   if ((*N2) == -999) return 0;
   if (nInputs_ > 21)
   {
      printf("Regression INFO: nInputs > 21 => no lattice point generated.\n");
      (*N2) = 0;
      return 0;
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
   for (ii = 1; ii < nInputs_; ii++) totPts = totPts * nPtsPerDim_;
   HX = new double[nInputs_];
   for (ii = 0; ii < nInputs_; ii++)
      HX[ii] = (upperBounds_[ii] - lowerBounds_[ii]) /
               (double) (nPtsPerDim_ - 1);

   (*X2) = new double[nInputs_ * totPts];
   (*Y2) = new double[totPts];
   (*N2) = totPts;
   Xloc  = new double[nInputs_];

   for (ii = 0; ii < nInputs_; ii++) Xloc[ii] = lowerBounds_[ii];

   for (mm = 0; mm < totPts; mm++)
   {
      for (ii = 0; ii < nInputs_; ii++ ) (*X2)[mm*nInputs_+ii] = Xloc[ii];
      for (ii = 0; ii < nInputs_; ii++ )
      {
         Xloc[ii] += HX[ii];
         if (Xloc[ii] < upperBounds_[ii] ||
             PABS(Xloc[ii] - upperBounds_[ii]) < 1.0E-7) break;
         else Xloc[ii] = lowerBounds_[ii];
      }
      (*Y2)[mm] = evaluatePoint(&((*X2)[mm*nInputs_]));
   }

   delete [] Xloc;
   delete [] HX;
   return 0;
}

// ************************************************************************
// Generate 1D mesh results (setting others to some nominal values) 
// ------------------------------------------------------------------------
int Regression::gen1DGridData(double *X, double *Y, int ind1,
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
// Generate 2D mesh results (setting others to some nominal values) 
// ------------------------------------------------------------------------
int Regression::gen2DGridData(double *X, double *Y, int ind1,
                              int ind2, double *settings, int *NN, 
                              double **XX, double **YY)
{
   int    totPts, mm, nn, ind;
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
         ind = mm * nPtsPerDim_ + nn;
         Xloc[ind1] = HX[0] * mm + lowerBounds_[ind1];
         Xloc[ind2] = HX[1] * nn + lowerBounds_[ind2];
         (*XX)[ind*2]   = Xloc[ind1];
         (*XX)[ind*2+1] = Xloc[ind2];
         (*YY)[ind] = evaluatePoint(Xloc);
      }
   }

   delete [] Xloc;
   delete [] HX;
   return 0;
}

// ************************************************************************
// Generate 3D mesh results (setting others to some nominal values) 
// ------------------------------------------------------------------------
int Regression::gen3DGridData(double *X, double *Y, int ind1,
                              int ind2, int ind3, double *settings, 
                              int *NN, double **XX, double **YY)
{
   int    totPts, mm, nn, pp, ind;
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
            ind = mm * nPtsPerDim_ * nPtsPerDim_ + nn * nPtsPerDim_ + pp;
            Xloc[ind1] = HX[0] * mm + lowerBounds_[ind1];
            Xloc[ind2] = HX[1] * nn + lowerBounds_[ind2];
            Xloc[ind3] = HX[2] * pp + lowerBounds_[ind3];
            (*XX)[ind*3]   = Xloc[ind1];
            (*XX)[ind*3+1] = Xloc[ind2];
            (*XX)[ind*3+2] = Xloc[ind3];
            (*YY)[ind] = evaluatePoint(Xloc);
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
int Regression::gen4DGridData(double *X, double *Y, int ind1, int ind2,
                              int ind3, int ind4, double *settings, 
                              int *NN, double **XX, double **YY)
{
   int    totPts, mm, nn, pp, qq, ind;
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
               ind = mm*nPtsPerDim_*nPtsPerDim_*nPtsPerDim_ +
                       nn*nPtsPerDim_*nPtsPerDim_ + pp*nPtsPerDim_ + qq;
               Xloc[ind1] = HX[0] * mm + lowerBounds_[ind1];
               Xloc[ind2] = HX[1] * nn + lowerBounds_[ind2];
               Xloc[ind3] = HX[2] * pp + lowerBounds_[ind3];
               Xloc[ind4] = HX[3] * qq + lowerBounds_[ind4];
               (*XX)[ind*4]   = Xloc[ind1];
               (*XX)[ind*4+1] = Xloc[ind2];
               (*XX)[ind*4+2] = Xloc[ind3];
               (*XX)[ind*4+3] = Xloc[ind4];
               (*YY)[ind] = evaluatePoint(Xloc);
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
double Regression::evaluatePoint(double *X)
{
   int    mm, nn, pp, qq, offset;
   double Y;

   if (regCoeffs_ == NULL) return 0.0;
   Y = regCoeffs_[0];
   for (mm = 0; mm < nInputs_; mm++) Y += regCoeffs_[mm+1] * X[mm];
   if (order_ >= 2)
   {
      offset = nInputs_ + 1;
      for (mm = 0; mm < nInputs_; mm++)
         for (nn = mm; nn < nInputs_; nn++)
            Y += (regCoeffs_[offset++] * X[mm] * X[nn]);
   }
   if (order_ >= 3)
   {
      for (mm = 0; mm < nInputs_; mm++)
         for (nn = mm; nn < nInputs_; nn++)
            for (pp = nn; pp < nInputs_; pp++)
               Y += regCoeffs_[offset++] * X[mm] * X[nn] * X[pp];
   }
   if (order_ >= 4)
   {
      for (mm = 0; mm < nInputs_; mm++)
         for (nn = mm; nn < nInputs_; nn++)
            for (pp = nn; pp < nInputs_; pp++)
               for (qq = pp; qq < nInputs_; qq++)
                  Y += regCoeffs_[offset++] * X[mm] * X[nn] * X[pp] * X[qq];
   }
   return Y;
}

// ************************************************************************
// Evaluate a number of points
// ------------------------------------------------------------------------
double Regression::evaluatePoint(int npts, double *X, double *Y)
{
   for (int kk = 0; kk < npts; kk++)
      Y[kk] = evaluatePoint(&X[kk*nInputs_]);
   return 0.0;
}

// ************************************************************************
// Evaluate a given point and also its standard deviation
// ------------------------------------------------------------------------
double Regression::evaluatePointFuzzy(double *X, double &std)
{
   int    mm, nn, pp, qq, cc, offset, nTimes=100;
   double dtmp, *Ys, C1, C2, mean, stds, *uppers, *lowers;
   PDFBase **PDFPtrs;

   if (regCoeffs_ == NULL) return 0.0;

   if (fuzzyC_ == NULL)
   {
      fuzzyC_ = new double*[numTerms_];
      for (mm = 0; mm < numTerms_; mm++) 
         fuzzyC_[mm] = new double[nTimes]; 
      PDFPtrs = new PDFBase*[numTerms_];
      uppers = new double[numTerms_];
      lowers = new double[numTerms_];
      for (mm = 0; mm < numTerms_; mm++) 
      {
         mean = regCoeffs_[mm];
         stds = regStdevs_[mm];
         uppers[mm] = mean + 2.0 * stds;
         lowers[mm] = mean - 2.0 * stds;
         PDFPtrs[mm] = (PDFBase *) new PDFNormal(mean, std);
      }
      Ys = new double[nTimes];

      mean = 0.0;
      for (cc = 0; cc < nTimes; cc++)
      {
         C1 = PSUADE_drand();
         PDFPtrs[0]->invCDF(1, &C1, &C2, lowers[0], uppers[0]);
         dtmp = C2;
         fuzzyC_[0][cc] = C2;

         for (mm = 0; mm < nInputs_; mm++)
         {
            C1 = PSUADE_drand();
            PDFPtrs[mm+1]->invCDF(1, &C1, &C2, lowers[mm+1], uppers[mm+1]);
            dtmp += C2 * X[mm];
            fuzzyC_[mm+1][cc] = C2;
         }
         if (order_ >= 2)
         {
            offset = nInputs_ + 1;
            for (mm = 0; mm < nInputs_; mm++)
            {
               for (nn = mm; nn < nInputs_; nn++)
               {
                  C1 = PSUADE_drand();
                  PDFPtrs[offset]->invCDF(1,&C1,&C2,lowers[offset],uppers[offset]);
                  dtmp += (C2 * X[mm] * X[nn]);
                  fuzzyC_[offset][cc] = C2;
                  offset++;
               }
            }
         }
         if (order_ >= 3)
         {
            for (mm = 0; mm < nInputs_; mm++)
            {
               for (nn = mm; nn < nInputs_; nn++)
               {
                  for (pp = nn; pp < nInputs_; pp++)
                  {
                     C1 = PSUADE_drand();
                     PDFPtrs[offset]->invCDF(1,&C1,&C2,lowers[offset],uppers[offset]);
                     fuzzyC_[offset][cc] = C2;
                     dtmp += (C2 * X[mm] * X[nn] * X[pp]);
                     offset++;
                  }
               }
            }
         }
         if (order_ >= 4)
         {
            for (mm = 0; mm < nInputs_; mm++)
            {
               for (nn = mm; nn < nInputs_; nn++)
               {
                  for (pp = nn; pp < nInputs_; pp++)
                  {
                     for (qq = pp; qq < nInputs_; qq++)
                     {
                        C1 = PSUADE_drand();
                        PDFPtrs[offset]->invCDF(1, &C1, &C2, lowers[offset],
                                                 uppers[offset]);
                        fuzzyC_[offset][cc] = C2;
                        dtmp += (C2 * X[mm] * X[nn] * X[pp] * X[qq]);
                        offset++;
                     }
                  }
               }
            }
         }
         Ys[cc] = dtmp;
         mean += dtmp;
      }
      mean /= (double) nTimes;
      stds = 0.0;
      for (cc = 0; cc < nTimes; cc++)
         stds += (Ys[cc] - mean) * (Ys[cc] - mean);
      stds = sqrt(stds / (nTimes - 1));

      for (mm = 0; mm < numTerms_; mm++) delete PDFPtrs[mm];
      delete [] PDFPtrs;
      delete [] uppers;
      delete [] lowers;
      delete [] Ys;
   }
   else
   {
      Ys = new double[nTimes];

      mean = 0.0;
      for (cc = 0; cc < nTimes; cc++)
      {
         C2 = fuzzyC_[0][cc];
         dtmp = C2;

         for (mm = 0; mm < nInputs_; mm++)
         {
            C2 = fuzzyC_[mm+1][cc];
            dtmp += C2 * X[mm];
         }
         if (order_ >= 2)
         {
            offset = nInputs_ + 1;
            for (mm = 0; mm < nInputs_; mm++)
            {
               for (nn = mm; nn < nInputs_; nn++)
               {
                  C2 = fuzzyC_[offset][cc];
                  dtmp += (C2 * X[mm] * X[nn]);
                  offset++;
               }
            }
         }
         if (order_ >= 3)
         {
            for (mm = 0; mm < nInputs_; mm++)
            {
               for (nn = mm; nn < nInputs_; nn++)
               {
                  for (pp = nn; pp < nInputs_; pp++)
                  {
                     C2 = fuzzyC_[offset][cc];
                     dtmp += (C2 * X[mm] * X[nn] * X[pp]);
                     offset++;
                  }
               }
            }
         }
         if (order_ >= 4)
         {
            for (mm = 0; mm < nInputs_; mm++)
            {
               for (nn = mm; nn < nInputs_; nn++)
               {
                  for (pp = nn; pp < nInputs_; pp++)
                  {
                     for (qq = pp; qq < nInputs_; qq++)
                     {
                        C2 = fuzzyC_[offset][cc];
                        dtmp += (C2 * X[mm] * X[nn] * X[pp] * X[qq]);
                        offset++;
                     }
                  }
               }
            }
         }
         Ys[cc] = dtmp;
         mean += dtmp;
      }
      mean /= (double) nTimes;
      stds = 0.0;
      for (cc = 0; cc < nTimes; cc++)
         stds += (Ys[cc] - mean) * (Ys[cc] - mean);
      stds = sqrt(stds / (nTimes - 1));
      delete [] Ys;
   }
   std = stds; 
   return mean;
}

// ************************************************************************
// Evaluate a number of points and also their standard deviations
// ------------------------------------------------------------------------
double Regression::evaluatePointFuzzy(int npts, double *X, double *Y,
                                      double *Ystd)
{
   for (int kk = 0; kk < npts; kk++)
      Y[kk] = evaluatePointFuzzy(&(X[kk*nInputs_]), Ystd[kk]);
   return 0.0;
}

// ************************************************************************
// set parameters
// ------------------------------------------------------------------------
double Regression::setParams(int targc, char **targv)
{
   order_ = *(int *) targv[0];
   if (order_ <= 0 || order_ > 4) order_ = 1;
   switch (order_)
   {
      case 1: faID_ = PSUADE_RS_REGR1; break;
      case 2: faID_ = PSUADE_RS_REGR2; break;
      case 3: faID_ = PSUADE_RS_REGR3; break;
      case 4: faID_ = PSUADE_RS_REGR4; break;
   }
   return 0.0;
}

// ************************************************************************
// perform mean/variance analysis
// ------------------------------------------------------------------------
int Regression::analyze(double *X, double *Y)
{
   int    M, N, ii, mm, nn, nn2, nn3, nn4, info, wlen, ind, last, NRevised;
   double *B, *XX, SSresid, SStotal, R2, *XTX, var, *Bstd;
   double esum, ymax, *txArray, *AA, *SS, *UU, *VV, *WW;
   char   jobu  = 'A', jobvt = 'A';
   char   pString[100];
   FILE   *fp;

   if (nInputs_ <= 0 || nSamples_ <= 0)
   {
      printf("Regression ERROR: consult PSUADE developers.\n");
      exit( 1 );
   } 
   if (nSamples_ <= nInputs_) return 0;
   
   txArray = new double[nSamples_];
   for (ii = 0; ii < nInputs_; ii++)
   {
      for (mm = 0; mm < nSamples_; mm++)
         txArray[mm] = X[mm*nInputs_+ii];
      sortDbleList(nSamples_, txArray);
      last = 1;
      for (mm = 1; mm < nSamples_; mm++)
      {
         if (txArray[mm] != txArray[last-1])
         {
            txArray[last] = txArray[mm];
            last++;
         }
      }
      if (order_ >= last)
      {
         order_ = last - 1;
         printf("* Regression: order reduced to %d - not enough levels.\n",
               order_);
         delete [] txArray;
         return -1;
      }
   }
   delete [] txArray;

   N = loadXMatrix(X, &XX); 
   if (N == 0)
   {
      printf("* Regression ERROR: loadXMatrix returns NULL.\n");
      return -1;
   }
   while (N > nSamples_ && order_ > 1)
   {
      if (XX != NULL) delete [] XX;
      printf("* Regression ERROR: sample too small for order %d.\n",order_);
      printf("                    Try lower order.\n");
      return -1;
   }
   if (N > nSamples_)
   {
      printf("* Regression: sample size too small (%d > %d).\n", N, nSamples_);
      if (XX != NULL) delete [] XX;
      return -1;
   }
   M = nSamples_;

   if (outputLevel_ >= 0)
   {
      printAsterisks(0);
      if (order_ == 1)
        printf("*               Linear Regression Analysis\n");
      else if (order_ == 2)
        printf("*             Quadratic Regression Analysis\n");
      else if (order_ == 3)
        printf("*               cubic Regression Analysis\n");
      else if (order_ == 4)
        printf("*               quartic Regression Analysis\n");
      printf("* R-squared gives a measure of the goodness of the model.\n");
      printf("* R-squared should be close to 1 if it is a good model.\n");
      printf("* TURN ON rs_expert mode to output regression matrix.\n");
      printf("* TURN ON rs_expert mode to condition covariance matrix.\n");
      printf("* SET print level to 5 to output regression error splot.\n");
      printf("* SET print level to 4 to output data standard deviations.\n");
      printDashes(0);
      printf("* Suggestion: if your parameter ranges are too high, SCALE\n");
      printf("*             them first using 'irerange' command in PSUADE\n");
      printf("*             command line mode.\n");
      printEquals(0);
   }

   wlen = 5 * M;
   AA = new double[M*N];
   UU = new double[M*M];
   SS = new double[N];
   VV = new double[M*N];
   WW = new double[wlen];
   for (mm = 0; mm < M; mm++) 
      for (nn = 0; nn < N; nn++) 
         AA[mm+nn*M] = sqrt(weights_[mm]) * XX[mm+nn*M];
   if (psRSExpertMode_ == 1)
   {
      fp = fopen("regression_matrix.m", "w");
      fprintf(fp, "%% the sample matrix where svd is computed\n");
      fprintf(fp, "%% the last column is the right hand side\n");
      fprintf(fp, "AA = [\n");
      for (mm = 0; mm < M; mm++) 
      {
         for (nn = 0; nn < N; nn++) 
            fprintf(fp, "%16.6e ", AA[mm+nn*M]);
         fprintf(fp, "%16.6e \n",Y[mm]);
      }
      fprintf(fp, "];\n");
      fprintf(fp, "A = AA(:,1:%d);\n", N);
      fprintf(fp, "B = AA(:,%d);\n", N+1);
      fclose(fp);
      printf("Regression: regression matrix (X^TX) now in regression_matrix.m\n");
   }
   if (outputLevel_ > 3) printf("Running SVD ...\n"); 
   dgesvd_(&jobu, &jobvt, &M, &N, AA, &M, SS, UU, &M, VV, &N, WW,
           &wlen, &info);
   if (outputLevel_ > 3) printf("SVD completed: status = %d (should be 0).\n",info); 

   if (info != 0)
   {
      printf("* Regression ERROR: dgesvd returns a nonzero (%d).\n",info);
      printf("* Regression terminates further processing.\n");
      delete [] XX;
      delete [] AA;
      delete [] UU;
      delete [] SS;
      delete [] VV;
      delete [] WW;
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
      printf("* Regression ERROR: true rank of sample = %d (need %d)\n",
             NRevised, N);
      printf("* This can be due to the quality of the sample.\n");
      delete [] XX;
      delete [] AA;
      delete [] UU;
      delete [] SS;
      delete [] VV;
      delete [] WW;
      return -1;
   }
   if (psRSExpertMode_ == 1)
   {
      printf("Regression: singular values for the Vandermonde matrix\n");
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
         if (SS[nn-1] > 0 && SS[nn]/SS[nn-1] < 1.0e-8)
         {
            SS[nn] = 0.0;
            NRevised--;
         }
      }
   }
   for (mm = 0; mm < N; mm++) 
   {
      WW[mm] = 0.0;
      for (nn = 0; nn < M; nn++) 
         WW[mm] += UU[mm*M+nn] * sqrt(weights_[nn]) * Y[nn]; 
   }
   for (nn = 0; nn < NRevised; nn++) WW[nn] /= SS[nn];
   for (nn = NRevised; nn < N; nn++) WW[nn] = 0.0;
   B = new double[N];
   for (mm = 0; mm < N; mm++) 
   {
      B[mm] = 0.0;
      for (nn = 0; nn < N; nn++) B[mm] += VV[mm*N+nn] * WW[nn]; 
   }
   delete [] SS;
   delete [] UU;
   delete [] VV;

   fp = NULL;
   if (outputLevel_ >= 5)
   {
      fp = fopen("regression_err_splots.m", "w");
      fprintf(fp, "%% This file contains scatter plots of\n");
      fprintf(fp, "%% errors with respect to each input.\n");
      fprintf(fp, "%% This is useful to see which input fits the worst.\n");
      fprintf(fp,"X = [\n");
      for (mm = 0; mm < nSamples_; mm++)
      {
         for (nn = 0; nn < nInputs_; nn++) 
            fprintf(fp,"%16.8e ", X[mm*nInputs_+nn]);
         fprintf(fp,"\n");
      }
      fprintf(fp,"];\n");
      fprintf(fp,"R = [\n");
   }

   esum = ymax = 0.0;
   for (mm = 0; mm < M; mm++)
   {
      WW[mm] = 0.0;
      for (nn = 0; nn < N; nn++) WW[mm] += XX[mm+nn*nSamples_] * B[nn];
      WW[mm] -= Y[mm];
      esum = esum + WW[mm] * WW[mm] * weights_[mm];
      if (fp != NULL) 
         fprintf(fp, "%6d %24.16e\n", mm+1, WW[mm]*sqrt(weights_[mm]));
      if (PABS(Y[mm]) > ymax) ymax = PABS(Y[mm]);
   }
   esum /= (double) nSamples_;
   esum = sqrt(esum);
   if (outputLevel_ > 1)
      printf("* Regression:: LS average error = %11.4e (Ymax=%9.2e)\n",
             esum, ymax); 
   if (fp != NULL) 
   {
      fprintf(fp,"];\n");
      fprintf(fp,"for ii = 1 : size(X,2)\n");
      fprintf(fp,"   plot(X(:,ii), R(:,2), 'x')\n");
      fprintf(fp,"   title('Residual Plot')\n");
      fprintf(fp,"   disp(['Residual versus input ' int2str(ii)])\n");
      fprintf(fp,"   disp('Press Enter to continue')\n");
      fprintf(fp,"   pause\n");
      fprintf(fp,"end\n");
      fclose(fp);
      printf("FILE regression_err_splots.m contains fitting errors.\n");
      printf("     (resubstitution test) with respect to each input.\n");
   }

   computeSS(N, XX, Y, B, SSresid, SStotal);
   R2  = (SStotal - SSresid) / SStotal;
   if (nSamples_ > N) var = SSresid / (double) (nSamples_ - N);
   else               var = 0.0;

   Bstd = new double[N];
   computeXTX(N, XX, &XTX);
   computeCoeffVariance(N, XTX, var, Bstd);

   if (outputLevel_ >= 0)
   {
      printRC(N, B, Bstd, XX, Y);
      printf("* Regression R-square = %12.4e\n", R2);
      if ((M - N - 1) > 0)
         printf("* adjusted   R-square = %12.4e\n",
                1.0 - (1.0 - R2) * ((M - 1) / (M - N - 1)));
      if (outputLevel_ > 1) printSRC(X, B, SStotal);
   }
   if (psRSExpertMode_ == 1)
   {
      fp = fopen("regression_function.m", "w");
      fprintf(fp, "Y = %24.16e ", B[0]);
      for (nn = 1; nn <= nInputs_; nn++)
         fprintf(fp, "+ %24.16e * X(%d) ", B[nn], nn);
      if (order_ >= 2)
      {
         ind = nInputs_ + 1;
         for (nn = 0; nn < nInputs_; nn++)
         {
            for (nn2 = nn; nn2 < nInputs_; nn2++)
            {
               fprintf(fp, "+ %24.16e * X(%d) * X(%d) ", B[ind], nn+1, nn2+1);
               ind++;
            }
         }
      }
      if (order_ >= 3)
      {
         for (nn = 0; nn < nInputs_; nn++)
         {
            for (nn2 = nn; nn2 < nInputs_; nn2++)
            {
               for (nn3 = nn2; nn3 < nInputs_; nn3++)
               {
                  fprintf(fp, "+ %24.16e * X(%d) * X(%d) * X(%d) ", B[ind], 
                          nn+1, nn2+1, nn3+1);
                  ind++;
               }
            }
         }
      }
      if (order_ >= 4)
      {
         for (nn = 0; nn < nInputs_; nn++)
         {
            for (nn2 = nn; nn2 < nInputs_; nn2++)
            {
               for (nn3 = nn2; nn3 < nInputs_; nn3++)
               {
                  for (nn4 = nn3; nn4 < nInputs_; nn4++)
                  {
                     fprintf(fp, "+ %24.16e * X(%d) * X(%d) * X(%d) *X(%d) ", 
                             B[ind], nn+1, nn2+1, nn3+1, nn4+1);
                     ind++;
                  }
               }
            }
         }
      }
      fprintf(fp, ";\n");
      printf("FILE regression_function.m contains the final polynomial\n");
      printf("     functional form.\n");
      fclose(fp);
   }
 
   if (regCoeffs_ != NULL) delete [] regCoeffs_;
   if (regStdevs_ != NULL) delete [] regStdevs_;
   numTerms_  = N;
   regCoeffs_ = B;
   regStdevs_ = Bstd;
   delete [] AA;
   delete [] XX;
   delete [] XTX;
   delete [] WW;
   return 0;
}

// *************************************************************************
// load the X matrix
// -------------------------------------------------------------------------
int Regression::loadXMatrix(double *X, double **XXOut)
{
   int    M, N=0, mm, nn, nn2, nn3, nn4, ind;
   double *XX=NULL;

   (*XXOut) = NULL;
   if (order_ <= 0 || order_ > 4) return 0;

   M = nSamples_;
   N = nInputs_ + 1;
   if (order_ >= 2) N += nInputs_ * (nInputs_ + 1) / 2;
   if (order_ >= 3)
   {
      for (nn = 0; nn < nInputs_; nn++)
         for (nn2 = nn; nn2 < nInputs_; nn2++)
            for (nn3 = nn2; nn3 < nInputs_; nn3++) N++;
   }
   if (order_ >= 4)
   {
      for (nn = 0; nn < nInputs_; nn++)
         for (nn2 = nn; nn2 < nInputs_; nn2++)
            for (nn3 = nn2; nn3 < nInputs_; nn3++)
               for (nn4 = nn3; nn4 < nInputs_; nn4++) N++;
   }
   if (N > M) return N;
   XX = new double[M*N];
   for (mm = 0; mm < M; mm++)
   {
      XX[mm] = 1.0;
      for (nn = 0; nn < nInputs_; nn++)
         XX[M*(nn+1)+mm] = X[mm*nInputs_+nn];
   }
   if (order_ >= 2)
   {
      ind = nInputs_ + 1;
      for (nn = 0; nn < nInputs_; nn++)
      {
         for (nn2 = nn; nn2 < nInputs_; nn2++)
         {
            for (mm = 0; mm < M; mm++)
               XX[M*ind+mm] = X[mm*nInputs_+nn] * X[mm*nInputs_+nn2];
            ind++;
         }
      }
   }
   if (order_ >= 3)
   {
      for (nn = 0; nn < nInputs_; nn++)
      {
         for (nn2 = nn; nn2 < nInputs_; nn2++)
         {
            for (nn3 = nn2; nn3 < nInputs_; nn3++)
            {
               for (mm = 0; mm < M; mm++)
                  XX[M*ind+mm] = X[mm*nInputs_+nn] * X[mm*nInputs_+nn2] * 
                                   X[mm*nInputs_+nn3];
               ind++;
            }
         }
      }
   }
   if (order_ >= 4)
   {
      for (nn = 0; nn < nInputs_; nn++)
      {
         for (nn2 = nn; nn2 < nInputs_; nn2++)
         {
            for (nn3 = nn2; nn3 < nInputs_; nn3++)
            {
               for (nn4 = nn3; nn4 < nInputs_; nn4++)
               {
                  for (mm = 0; mm < M; mm++)
                     XX[M*ind+mm] = X[mm*nInputs_+nn] * X[mm*nInputs_+nn2] * 
                                      X[mm*nInputs_+nn3] * X[mm*nInputs_+nn4];
                  ind++;
               }
            }
         }
      }
   }
   (*XXOut) = XX;
   return N;
}

// *************************************************************************
// form X^T X 
// -------------------------------------------------------------------------
int Regression::computeXTX(int N, double *X, double **XXOut)
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
int Regression::computeSS(int N, double *XX, double *Y,
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
// compute coefficient variances (diagonal of sigma^2 (X' X)^(-1))
// -------------------------------------------------------------------------
int Regression::computeCoeffVariance(int N,double *XX,double var,double *B)
{
   int    nn, nn2, lwork, iOne=1, info;
   double *B2, *work, *XT;
   char   trans[1];

   (*trans) = 'N';
   B2 = new double[N];
   XT = new double[N*N];
   lwork = 2 * N * N;
   work  = new double[lwork];
   for (nn = 0; nn < N; nn++)
   {
      for (nn2 = 0; nn2 < N*N; nn2++) XT[nn2] = XX[nn2];
      for (nn2 = 0; nn2 < N; nn2++) B2[nn2] = 0.0;
      B2[nn] = var;
      dgels_(trans, &N, &N, &iOne, XT, &N, B2, &N, work, &lwork, &info);
      if (info != 0)
         printf("Regression WARNING: dgels returns error %d.\n",info);
      if (B2[nn] < 0) B[nn] = -sqrt(-B2[nn]);
      else            B[nn] = sqrt(B2[nn]);
   }

   delete [] B2;
   delete [] work;
   delete [] XT;
   return info;
}

// *************************************************************************
// print statistics
// -------------------------------------------------------------------------
int Regression::printRC(int N,double *B,double *Bvar,double *XX,double *Y)
{
   int    nn, ii, ind, nn2, nn3, nn4, mm;
   double coef, Bmax;
   char   fname[200];
   FILE   *fp;

   if (order_ <= 0 || order_ > 4) return 0;
   printf("* ======== Note: entries with t-value < 1 suppressed ========\n");
   printf("*            ");
   for (ii = 1; ii < order_; ii++) printf("    ");
   printf("  coefficient   std. error   t-value\n");

   if (PABS(Bvar[0]) < 1.0e-15) coef = 0.0;
   else                         coef = B[0] / Bvar[0]; 
   if (PABS(coef) > 1.0)
   {
      printf("* Constant  ");
      for (ii = 1; ii < order_; ii++) printf("    ");
      printf("= %12.4e %12.4e %12.4e\n", B[0], Bvar[0], coef);
   }
   for (nn = 1; nn <= nInputs_; nn++)
   {
      if (PABS(Bvar[nn]) < 1.0e-15) coef = 0.0;
      else                          coef = B[nn] / Bvar[nn]; 
      {
         printf("* Input %3d ", nn);
         for (ii = 1; ii < order_; ii++) printf("    ");
         printf("= %12.4e %12.4e %12.4e\n", B[nn], Bvar[nn], coef);
      }
      strcpy(fname, "dataVariance1");
   }
   if (order_ >= 2)
   {
      Bmax = 0.0;
      for (nn = nInputs_+1; nn < N; nn++)
         if (PABS(B[nn]) > Bmax) Bmax = PABS(B[nn]); 
      ind = nInputs_ + 1;
      for (nn = 0; nn < nInputs_; nn++)
      {
         for (nn2 = nn; nn2 < nInputs_; nn2++)
         {
            if (PABS(B[ind]) > 1.0e-6 * Bmax) 
            {
               if (PABS(Bvar[ind]) < 1.0e-15) coef = 0.0;
               else coef = B[ind] / Bvar[ind]; 
               {
                  printf("* Input %3d %3d ", nn+1, nn2+1);
                  for (ii = 2; ii < order_; ii++) printf("    ");
                  printf("= %12.4e %12.4e %12.4e\n",B[ind],Bvar[ind],coef); 
               }
            }
            ind++;
         }
      }
      strcpy(fname, "dataVariance2");
   }
   if (order_ >= 3)
   {
      Bmax = 0.0;
      for (nn = ind; nn < N; nn++)
         if (PABS(B[nn]) > Bmax) Bmax = PABS(B[nn]); 
      for (nn = 0; nn < nInputs_; nn++)
      {
         for (nn2 = nn; nn2 < nInputs_; nn2++)
         {
            for (nn3 = nn2; nn3 < nInputs_; nn3++)
            {
               if (PABS(B[ind]) > 1.0e-6 * Bmax) 
               {
                  if (PABS(Bvar[ind]) < 1.0e-15) coef = 0.0;
                  else coef = B[ind] / Bvar[ind]; 
                  {
                     printf("* Input %3d %3d %3d ", nn+1, nn2+1, nn3+1);
                     for (ii = 3; ii < order_; ii++) printf("    ");
                     printf("= %12.4e %12.4e %12.4e\n",B[ind],Bvar[ind],
                            coef);
                  }
               }
               ind++;
            }
         }
      }
      strcpy(fname, "dataVariance3");
   }
   if (order_ >= 4)
   {
      for (nn = 0; nn < nInputs_; nn++)
      {
         for (nn2 = nn; nn2 < nInputs_; nn2++)
         {
            for (nn3 = nn2; nn3 < nInputs_; nn3++)
            {
               for (nn4 = nn3; nn4 < nInputs_; nn4++)
               {
                  if (PABS(B[ind]) > 1.0e-6 * Bmax) 
                  {
                     if (PABS(Bvar[ind]) < 1.0e-15) coef = 0.0;
                     else coef = B[ind] / Bvar[ind]; 
                     {
                        printf("* Input %3d %3d %3d %3d ",nn+1,nn2+1,nn3+1,
                               nn4+1);
                        for (ii = 4; ii < order_; ii++) printf("    ");
                        printf("= %12.4e %12.4e %12.4e\n",B[ind],Bvar[ind],
                               coef); 
                     }
                  }
                  ind++;
               }
            }
         }
      }
      strcpy(fname, "dataVariance4");
   }
   printDashes(0);

   if (order_ >= 1 && order_ <= 4 && outputLevel_ > 3)
   {
      fp = fopen(fname, "w");
      fprintf(fp, "data number     data         standard error\n");
      for (mm = 0; mm < nSamples_; mm++)
      {
         coef = 0.0;
         for (nn = 0; nn < N; nn++) coef += PABS(XX[nn*nSamples_+mm]*Bvar[nn]);
         fprintf(fp,"%7d         %12.4e %12.4e\n",mm+1,Y[mm],sqrt(coef));
      }
      fclose(fp);
      printf("FILE %s contains output data standard errors.\n",fname);
   }
   return 0;
}

// *************************************************************************
// print standardized regression coefficients
// -------------------------------------------------------------------------
int Regression::printSRC(double *X, double *B, double SStotal)
{
   int    nn, mm, ind, ii, nn2, itemp, nn3, nn4, *iArray;
   double denom, xmean, coef, Bmax, coef1, coef2, xmean1, xmean2;
   double xmean3, coef3, xmean4, coef4, *B2;

   if (order_ <= 0 || order_ > 4) return 0;

   printAsterisks(0);
   printf("*       Standardized Regression Coefficients (SRC)\n");
   printf("* When R-square is acceptable (order assumption holds), the\n");
   printf("* absolute values of SRCs provide variable importance.\n"); 
   printEquals(0);
   printf("* based on nSamples = %d\n", nSamples_);

   B2 = new double[nSamples_];
   denom = sqrt(SStotal / (double) (nSamples_ - 1));
   Bmax = 0.0;
   for (nn = 0; nn < nInputs_; nn++)
   {
      xmean = 0.0;
      for (mm = 0; mm < nSamples_; mm++) 
         xmean += X[mm*nInputs_+nn] * sqrt(weights_[mm]);
      xmean /= (double) nSamples_;
      coef = 0.0;
      for (mm = 0; mm < nSamples_; mm++)
      {
         coef += (sqrt(weights_[mm])*X[mm*nInputs_+nn] - xmean) * 
                 (sqrt(weights_[mm])*X[mm*nInputs_+nn] - xmean);
      }
      coef = sqrt(coef / (double) (nSamples_ - 1)) / denom;
      printf("* Input %3d ", nn+1);
      for (ii = 1; ii < order_; ii++) printf("    ");
      if ((B[nn+1]*coef) > Bmax) Bmax = B[nn+1]*coef;
      printf("= %12.4e\n", B[nn+1]*coef);
      B2[nn] = PABS(B[nn+1] * coef);
   }
   printDashes(0);
   printf("*    ordered list of SRCs\n");
   printDashes(0);
   iArray = new int[nInputs_];
   for (nn = 0; nn < nInputs_; nn++) iArray[nn] = nn;
   sortDbleList2a(nInputs_, B2, iArray); 
   for (nn = nInputs_-1; nn >= 0; nn--)
   {
      printf("* Input %3d ", iArray[nn]+1);
      for (ii = 1; ii < order_; ii++) printf("    ");
      printf("= %12.4e\n", B2[nn]);
   }
   delete [] iArray; 

   if (order_ >= 2)
   {
      ind = nInputs_ + 1;
      for (nn = 0; nn < nInputs_; nn++)
      {
         xmean1 = 0.0;
         for (mm = 0; mm < nSamples_; mm++) 
            xmean1 += X[mm*nInputs_+nn] * sqrt(weights_[mm]);
         xmean1 /= (double) nSamples_;
         coef1 = 0.0;
         for (mm = 0; mm < nSamples_; mm++)
            coef1 += (sqrt(weights_[mm])*X[mm*nInputs_+nn] - xmean1) * 
                     (sqrt(weights_[mm])*X[mm*nInputs_+nn] - xmean1);
         coef1 = sqrt(coef1 / (double) (nSamples_ - 1));
         for (nn2 = nn; nn2 < nInputs_; nn2++)
         {
            xmean2 = 0.0;
            for (mm = 0; mm < nSamples_; mm++)
               xmean2 += X[mm*nInputs_+nn2] * sqrt(weights_[mm]);
            xmean2 /= (double) nSamples_;
            coef2 = 0.0;
            for (mm = 0; mm < nSamples_; mm++)
               coef2 += (sqrt(weights_[mm])*X[mm*nInputs_+nn2] - xmean2) * 
                        (sqrt(weights_[mm])*X[mm*nInputs_+nn2] - xmean2);
            coef2 = sqrt(coef2 / (double) (nSamples_ - 1));
            B2[ind] = B[ind] * coef1 * coef2 / denom;
            if (PABS(B2[ind]) > Bmax) Bmax = PABS(B2[ind]);
            ind++;
         }
      }
      ind = nInputs_ + 1;
      for (nn = 0; nn < nInputs_; nn++)
      {
         for (nn2 = nn; nn2 < nInputs_; nn2++)
         {
            if (PABS(B2[ind]) > 1.0e-3 * Bmax)
            {
               printf("* Input %3d,%3d ", nn+1, nn2+1);
               for (ii = 2; ii < order_; ii++) printf("    ");
               printf("= %12.4e\n", B2[ind]);
            }
            ind++;
         }
      }
   }
   if (order_ >= 3)
   {
      itemp = ind;
      for (nn = 0; nn < nInputs_; nn++)
      {
         xmean1 = 0.0;
         for (mm = 0; mm < nSamples_; mm++)
            xmean1 += X[mm*nInputs_+nn] * sqrt(weights_[mm]);
         xmean1 /= (double) nSamples_;
         coef1 = 0.0;
         for (mm = 0; mm < nSamples_; mm++)
            coef1 += (sqrt(weights_[mm])*X[mm*nInputs_+nn] - xmean1) * 
                     (sqrt(weights_[mm])*X[mm*nInputs_+nn] - xmean1);
         coef1 = sqrt(coef1 / (double) (nSamples_ - 1));
         for (nn2 = nn; nn2 < nInputs_; nn2++)
         {
            xmean2 = 0.0;
            for (mm = 0; mm < nSamples_; mm++)
               xmean2 += X[mm*nInputs_+nn2] * sqrt(weights_[mm]);
            xmean2 /= (double) nSamples_;
            coef2 = 0.0;
            for (mm = 0; mm < nSamples_; mm++)
               coef2 += (sqrt(weights_[mm])*X[mm*nInputs_+nn2] - xmean2) * 
                        (sqrt(weights_[mm])*X[mm*nInputs_+nn2] - xmean2);
            coef2 = sqrt(coef2 / (double) (nSamples_ - 1));
            for (nn3 = nn2; nn3 < nInputs_; nn3++)
            {
               xmean3 = 0.0;
               for (mm = 0; mm < nSamples_; mm++)
                  xmean3 += X[mm*nInputs_+nn3] * sqrt(weights_[mm]);
               xmean3 /= (double) nSamples_;
               coef3 = 0.0;
               for (mm = 0; mm < nSamples_; mm++)
                  coef3 += (sqrt(weights_[mm])*X[mm*nInputs_+nn3] - xmean3) * 
                           (sqrt(weights_[mm])*X[mm*nInputs_+nn3] - xmean3);
               coef3 = sqrt(coef3 / (double) (nSamples_ - 1));
               B2[ind] = B[ind] * coef1 * coef2 *coef3 / denom;
               if (PABS(B2[ind]) > Bmax) Bmax = PABS(B2[ind]);
               ind++;
            }
         }
      }
      ind = itemp;
      for (nn = 0; nn < nInputs_; nn++)
      {
         for (nn2 = nn; nn2 < nInputs_; nn2++)
         {
            for (nn3 = nn2; nn3 < nInputs_; nn3++)
            {
               if (PABS(B2[ind]) > 1.0e-3 * Bmax)
               {
                  printf("* Input %3d,%3d,%3d ", nn+1, nn2+1, nn3+1);
                  for (ii = 3; ii < order_; ii++) printf("    ");
                  printf("= %12.4e\n", B2[ind]);
               }
               ind++;
            }
         }
      }
   }
   if (order_ >= 4)
   {
      itemp = ind;
      for (nn = 0; nn < nInputs_; nn++)
      {
         xmean1 = 0.0;
         for (mm = 0; mm < nSamples_; mm++)
            xmean1 += X[mm*nInputs_+nn] * sqrt(weights_[mm]);
         xmean1 /= (double) nSamples_;
         coef1 = 0.0;
         for (mm = 0; mm < nSamples_; mm++)
            coef1 += (sqrt(weights_[mm])*X[mm*nInputs_+nn] - xmean1) * 
                     (sqrt(weights_[mm])*X[mm*nInputs_+nn] - xmean1);
         coef1 = sqrt(coef1 / (double) (nSamples_ - 1));
         for (nn2 = nn; nn2 < nInputs_; nn2++)
         {
            xmean2 = 0.0;
            for (mm = 0; mm < nSamples_; mm++)
               xmean2 += X[mm*nInputs_+nn2] * sqrt(weights_[mm]);
            xmean2 /= (double) nSamples_;
            coef2 = 0.0;
            for (mm = 0; mm < nSamples_; mm++)
               coef2 += (sqrt(weights_[mm])*X[mm*nInputs_+nn2] - xmean2) * 
                        (sqrt(weights_[mm])*X[mm*nInputs_+nn2] - xmean2);
            coef2 = sqrt(coef2 / (double) (nSamples_ - 1));
            for (nn3 = nn2; nn3 < nInputs_; nn3++)
            {
               xmean3 = 0.0;
               for (mm = 0; mm < nSamples_; mm++)
                  xmean3 += X[mm*nInputs_+nn3] * sqrt(weights_[mm]);
               xmean3 /= (double) nSamples_;
               coef3 = 0.0;
               for (mm = 0; mm < nSamples_; mm++)
                  coef3 += (sqrt(weights_[mm])*X[mm*nInputs_+nn3] - xmean3) * 
                           (sqrt(weights_[mm])*X[mm*nInputs_+nn3] - xmean3);
               coef3 = sqrt(coef3 / (double) (nSamples_ - 1));
               for (nn4 = nn3; nn4 < nInputs_; nn4++)
               {
                  xmean4 = 0.0;
                  for (mm = 0; mm < nSamples_; mm++)
                     xmean3 += X[mm*nInputs_+nn4] * sqrt(weights_[mm]);
                  xmean4 /= (double) nSamples_;
                  coef4 = 0.0;
                  for (mm = 0; mm < nSamples_; mm++)
                     coef4 += (sqrt(weights_[mm])*X[mm*nInputs_+nn4] - xmean4) * 
                              (sqrt(weights_[mm])*X[mm*nInputs_+nn4] - xmean4);
                  coef4 = sqrt(coef4 / (double) (nSamples_ - 1));
                  B2[ind] = B[ind] * coef1 * coef2 *coef3 *coef4 / denom;
                  if (PABS(B2[ind]) > Bmax) Bmax = PABS(B2[ind]);
                  ind++;
               }
            }
         }
      }
      ind = itemp;
      for (nn = 0; nn < nInputs_; nn++)
      {
         for (nn2 = nn; nn2 < nInputs_; nn2++)
         {
            for (nn3 = nn2; nn3 < nInputs_; nn3++)
            {
               for (nn4 = nn3; nn4 < nInputs_; nn4++)
               {
                  if (PABS(B2[ind]) > 1.0e-3 * Bmax)
                  {
                     printf("* Input %3d,%3d,%3d,%3d ",nn+1,nn2+1,nn3+1,nn4+1);
                     for (ii = 4; ii < order_; ii++) printf("    ");
                     printf("= %12.4e\n",B2[ind]);
                  }
                  ind++;
               }
            }
         }
      }
   }
   delete [] B2;
   printAsterisks(0);
   return 0;
}

