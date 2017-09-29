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
// Functions for the class GP1
// AUTHOR : CHARLES TONG
// DATE   : 2005
// ************************************************************************

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "GP1.h"
#include "Util/sysdef.h"
#include "Util/PsuadeUtil.h"

#define PABS(x) ((x) > 0 ? (x) : (-(x)))

#ifdef HAVE_TPROS
extern "C" 
{
  void TprosTrain(int nInputs, int nTrains, double *trainInputs,
                  double *trainOutput, int, double *, double *);

  void TprosInterp(int nTests, double *inputs, double *output, double *stds);

  void TprosGetLengthScales(int nInputs, double *lengthScales);
}
#endif

// ************************************************************************
// Constructor for object class GP1
// ------------------------------------------------------------------------
GP1::GP1(int nInputs,int nSamples) : FuncApprox(nInputs,nSamples)
{
   faID_ = PSUADE_RS_GP1;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
GP1::~GP1()
{
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int GP1::genNDGridData(double *X, double *Y, int *N, double **X2, 
                      double **Y2)
{
#ifdef HAVE_TPROS
   int    totPts, ss, ii;
   double *HX, *Xloc, *XX, *YY, *stds;

   stds = new double[nSamples_];
   if (outputLevel_ >= 1) printf("GP1 training begins....\n");
   TprosTrain(nInputs_, nSamples_, X, Y, 0, NULL, stds);
   for (ss = 0; ss < nSamples_; ss++) stds[ss] = 0.0;
   if (outputLevel_ >= 1) printf("GP1 training completed.\n");
   delete [] stds;
   if ((*N) == -999 || X2 == NULL || Y2 == NULL) return 0;
  
   totPts = nPtsPerDim_;
   for (ii = 1; ii < nInputs_; ii++)
      totPts = totPts * nPtsPerDim_;
   HX = new double[nInputs_];
   for (ii = 0; ii < nInputs_; ii++) 
      HX[ii] = (upperBounds_[ii] - lowerBounds_[ii]) /
               (double) (nPtsPerDim_ - 1); 
 
   XX = new double[totPts*nInputs_];
   YY = new double[totPts];
   Xloc = new double[nInputs_];
 
   for (ii = 0; ii < nInputs_; ii++) Xloc[ii] = lowerBounds_[ii];
 
   for (ss = 0; ss < totPts; ss++)
   {
      for (ii = 0; ii < nInputs_; ii++ ) XX[ss*nInputs_+ii] = Xloc[ii];
      for (ii = 0; ii < nInputs_; ii++ ) 
      {
         Xloc[ii] += HX[ii];
         if (Xloc[ii] < upperBounds_[ii] || 
             PABS(Xloc[ii] - upperBounds_[ii]) < 1.0E-7) break;
         else Xloc[ii] = lowerBounds_[ii];
      }
   }
 
   if (outputLevel_ >= 1) printf("GP1 interpolation begins....\n");
   TprosInterp(totPts, XX, YY, NULL);
   if (outputLevel_ >= 1) printf("GP1 interpolation completed.\n");
   (*N) = totPts;
   (*X2) = XX;
   (*Y2) = YY;
   delete [] Xloc;
   delete [] HX;
#else
   printf("PSUADE ERROR : GP1 not installed.\n");
   return -1;
#endif
   return 0;
}

// ************************************************************************
// Generate 1D results for display
// ------------------------------------------------------------------------
int GP1::gen1DGridData(double *X, double *Y, int ind1,
                       double *settings, int *n, double **X2, double **Y2)
{
#ifdef HAVE_TPROS
   int    ii, kk, totPts;
   double HX, *XX, *YY;

   if (outputLevel_ >= 1) printf("GP1 training begins....\n");
   TprosTrain(nInputs_, nSamples_, X, Y, 0, NULL, NULL);
   if (outputLevel_ >= 1) printf("GP1 training completed.\n");
  
   totPts = nPtsPerDim_;
   HX = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 

   (*X2) = new double[totPts];
   XX = new double[totPts*nInputs_];
   for (ii = 0; ii < nPtsPerDim_; ii++) 
   {
      for (kk = 0; kk < nInputs_; kk++) 
         XX[ii*nInputs_+kk] = settings[kk]; 
      XX[ii*nInputs_+ind1] = HX * ii + lowerBounds_[ind1];
      (*X2)[ii] = HX * ii + lowerBounds_[ind1];
   }
    
   YY = new double[totPts];
   if (outputLevel_ >= 1) printf("GP1 interpolation begins....\n");
   TprosInterp(totPts, XX, YY, NULL);
   if (outputLevel_ >= 1) printf("GP1 interpolation completed.\n");
   (*n) = totPts;
   (*Y2) = YY;
   delete [] XX;
#else
   printf("PSUADE ERROR : GP1 not installed.\n");
#endif
   return 0;
}

// ************************************************************************
// Generate 2D results for display
// ------------------------------------------------------------------------
int GP1::gen2DGridData(double *X, double *Y, int ind1, int ind2, 
                       double *settings, int *n, double **X2, double **Y2)
{
#ifdef HAVE_TPROS
   int    ii, jj, kk, totPts, index;
   double *HX, *XX, *YY;

   if (outputLevel_ >= 1) printf("GP1 training begins....\n");
   TprosTrain(nInputs_, nSamples_, X, Y, 0, NULL, NULL);
   if (outputLevel_ >= 1) printf("GP1 training completed.\n");
  
   totPts = nPtsPerDim_ * nPtsPerDim_;
   HX    = new double[2];
   HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
   HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 

   XX = new double[totPts*nInputs_];
   (*X2) = new double[2*totPts];
   for (ii = 0; ii < nPtsPerDim_; ii++) 
   {
      for (jj = 0; jj < nPtsPerDim_; jj++) 
      {
         index = ii * nPtsPerDim_ + jj;
         for (kk = 0; kk < nInputs_; kk++) 
            XX[index*nInputs_+kk] = settings[kk]; 
         XX[index*nInputs_+ind1]  = HX[0] * ii + lowerBounds_[ind1];
         XX[index*nInputs_+ind2]  = HX[1] * jj + lowerBounds_[ind2];
         (*X2)[index*2]   = HX[0] * ii + lowerBounds_[ind1];
         (*X2)[index*2+1] = HX[1] * jj + lowerBounds_[ind2];
      }
   }
    
   YY = new double[totPts];
   if (outputLevel_ >= 1) printf("GP1 interpolation begins....\n");
   TprosInterp(totPts, XX, YY, NULL);
   if (outputLevel_ >= 1) printf("GP1 interpolation completed.\n");
   (*n) = totPts;
   (*Y2) = YY;
   delete [] XX;
   delete [] HX;
#else
   printf("PSUADE ERROR : GP1 not installed.\n");
#endif
   return 0;
}

// ************************************************************************
// Generate 3D results for display
// ------------------------------------------------------------------------
int GP1::gen3DGridData(double *X, double *Y, int ind1, int ind2, int ind3,
                       double *settings, int *n, double **X2, double **Y2)
{
#ifdef HAVE_TPROS
   int    ii, jj, ll, kk, totPts, index;
   double *HX, *XX, *YY;

   if (outputLevel_ >= 1) printf("GP1 training begins....\n");
   TprosTrain(nInputs_, nSamples_, X, Y, 0, NULL, NULL);
   if (outputLevel_ >= 1) printf("GP1 training completed.\n");
  
   totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
   HX    = new double[3];
   HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
   HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 
   HX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_ - 1); 

   XX = new double[totPts*nInputs_];
   (*X2) = new double[3*totPts];
   for (ii = 0; ii < nPtsPerDim_; ii++) 
   {
      for (jj = 0; jj < nPtsPerDim_; jj++) 
      {
         for (ll = 0; ll < nPtsPerDim_; ll++) 
         {
            index = ii * nPtsPerDim_ * nPtsPerDim_ + jj * nPtsPerDim_ + ll;
            for (kk = 0; kk < nInputs_; kk++) 
               XX[index*nInputs_+kk] = settings[kk]; 
            XX[index*nInputs_+ind1]  = HX[0] * ii + lowerBounds_[ind1];
            XX[index*nInputs_+ind2]  = HX[1] * jj + lowerBounds_[ind2];
            XX[index*nInputs_+ind3]  = HX[2] * ll + lowerBounds_[ind3];
            (*X2)[index*3]   = HX[0] * ii + lowerBounds_[ind1];
            (*X2)[index*3+1] = HX[1] * jj + lowerBounds_[ind2];
            (*X2)[index*3+2] = HX[2] * ll + lowerBounds_[ind3];
         }
      }
   }
    
   YY = new double[totPts];
   if (outputLevel_ >= 1) printf("GP1 interpolation begins....\n");
   TprosInterp(totPts, XX, YY, NULL);
   if (outputLevel_ >= 1) printf("GP1 interpolation completed.\n");
   (*n) = totPts;
   (*Y2) = YY;
   delete [] XX;
   delete [] HX;
#else
   printf("PSUADE ERROR : GP1 not installed.\n");
#endif
   return 0;
}

// ************************************************************************
// Generate 4D results for display
// ------------------------------------------------------------------------
int GP1::gen4DGridData(double *X, double *Y, int ind1, int ind2, int ind3,
                       int ind4, double *settings, int *n, double **X2, 
                       double **Y2)
{
#ifdef HAVE_TPROS
   int    ii, jj, ll, mm, kk, totPts, index;
   double *HX, *XX, *YY;

   if (outputLevel_ >= 1) printf("GP1 training begins....\n");
   TprosTrain(nInputs_, nSamples_, X, Y, 0, NULL, NULL);
   if (outputLevel_ >= 1) printf("GP1 training completed.\n");
  
   totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
   HX    = new double[4];
   HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
   HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 
   HX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_ - 1); 
   HX[3] = (upperBounds_[ind4] - lowerBounds_[ind4]) / (nPtsPerDim_ - 1); 

   XX = new double[totPts*nInputs_];
   (*X2) = new double[4*totPts];
   for (ii = 0; ii < nPtsPerDim_; ii++) 
   {
      for (jj = 0; jj < nPtsPerDim_; jj++) 
      {
         for (ll = 0; ll < nPtsPerDim_; ll++) 
         {
            for (mm = 0; mm < nPtsPerDim_; mm++) 
            {
               index = ii*nPtsPerDim_*nPtsPerDim_*nPtsPerDim_ +
                       jj*nPtsPerDim_*nPtsPerDim_ + ll*nPtsPerDim_ + mm;
               for (kk = 0; kk < nInputs_; kk++) 
                  XX[index*nInputs_+kk] = settings[kk]; 
               XX[index*nInputs_+ind1]  = HX[0] * ii + lowerBounds_[ind1];
               XX[index*nInputs_+ind2]  = HX[1] * jj + lowerBounds_[ind2];
               XX[index*nInputs_+ind3]  = HX[2] * ll + lowerBounds_[ind3];
               XX[index*nInputs_+ind4]  = HX[3] * mm + lowerBounds_[ind4];
               (*X2)[index*4]   = HX[0] * ii + lowerBounds_[ind1];
               (*X2)[index*4+1] = HX[1] * jj + lowerBounds_[ind2];
               (*X2)[index*4+2] = HX[2] * ll + lowerBounds_[ind3];
               (*X2)[index*4+3] = HX[3] * mm + lowerBounds_[ind4];
            }
         }
      }
   }
    
   YY = new double[totPts];
   if (outputLevel_ >= 1) printf("GP1 interpolation begins....\n");
   TprosInterp(totPts, XX, YY, NULL);
   if (outputLevel_ >= 1) printf("GP1 interpolation completed.\n");
   (*n) = totPts;
   (*Y2) = YY;
   delete [] XX;
   delete [] HX;
#else
   printf("PSUADE ERROR : GP1 not installed.\n");
#endif
   return 0;
}

// ************************************************************************
// Evaluate a given point 
// ------------------------------------------------------------------------
double GP1::evaluatePoint(double *X)
{
   double Y=0.0;
#ifdef HAVE_TPROS
   int    iOne=1;
   TprosInterp(iOne, X, &Y, NULL);
#else
   printf("PSUADE ERROR : GP1 not installed.\n");
#endif
   return Y;
}

// ************************************************************************
// Evaluate a number of points 
// ------------------------------------------------------------------------
double GP1::evaluatePoint(int npts, double *X, double *Y)
{
#ifdef HAVE_TPROS
   TprosInterp(npts, X, Y, NULL);
#else
   printf("PSUADE ERROR : GP1 not installed.\n");
#endif
   return 0.0;
}

// ************************************************************************
// Evaluate a given point, return also the standard deviation 
// ------------------------------------------------------------------------
double GP1::evaluatePointFuzzy(double *X, double &std)
{
   double Y=0.0;
#ifdef HAVE_TPROS
   TprosInterp(1, X, &Y, &std);
   if (std < 0) printf("GP1 ERROR: variance < 0\n");
   else         std = sqrt(std);
#else
   printf("PSUADE ERROR : GP1 not installed.\n");
#endif
   return Y;
}

// ************************************************************************
// Evaluate a number of points, return also the standard deviation 
// ------------------------------------------------------------------------
double GP1::evaluatePointFuzzy(int npts,double *X, double *Y, double *Ystd)
{
#ifdef HAVE_TPROS
   TprosInterp(npts, X, Y, Ystd);
   for (int ii = 0; ii < npts; ii++)
   {
      if (Ystd[ii] < 0) printf("GP1 ERROR: variance < 0\n");
      else              Ystd[ii] = sqrt(Ystd[ii]);
   }
#else
   printf("PSUADE ERROR : GP1 not installed.\n");
#endif
   return 0.0;
}

// ************************************************************************
// set parameters
// ------------------------------------------------------------------------
double GP1::setParams(int targc, char **targv)
{
   int    ii, *iArray, ind;
   double *lengthScales=NULL, mmax, ddata=0.0, range;
                                                                                
   if (targc > 0 && !strcmp(targv[0], "rank"))
   {
      lengthScales = new double[nInputs_];
#ifdef HAVE_TPROS
      TprosGetLengthScales(nInputs_, lengthScales);
#else
      printf("PSUADE ERROR : GP1 not installed.\n");
#endif
      mmax = 0.0;
      for (ii = 0; ii < nInputs_; ii++)
      {
         lengthScales[ii] = 1.0/lengthScales[ii];
         range = upperBounds_[ii] - lowerBounds_[ii];
         lengthScales[ii] *= range;
         if (lengthScales[ii] > mmax) mmax = lengthScales[ii];
      }
      for (ii = 0; ii < nInputs_; ii++)
         lengthScales[ii] = lengthScales[ii] / mmax * 100.0;
      iArray = new int[nInputs_];
      for (ii = 0; ii < nInputs_; ii++) iArray[ii] = ii;
      sortDbleList2a(nInputs_, lengthScales, iArray);
      if (targc == 1)
      {
         printAsterisks(0);
         printf("* GP1 screening rankings \n");
         printAsterisks(0);
         for (ii = nInputs_-1; ii >= 0; ii--)
            printf("*  Rank %3d : Input = %3d (score = %4.1f) (ref = %e)\n", 
                   nInputs_-ii, iArray[ii]+1, lengthScales[ii], lengthScales[ii]*mmax*0.01);
         printAsterisks(0);
      }
      if (targc > 1)
      {
         ind = *(int *) targv[1];
         if (ind >= 0 && ind < nInputs_) ddata = lengthScales[ind];
      }
      delete [] lengthScales;
   }
   return ddata;
}

