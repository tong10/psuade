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
// Functions for the class Mars
// AUTHOR : CHARLES TONG
// DATE   : 2003
// ************************************************************************
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "Mars.h"
#include "Util/sysdef.h"
#include "Util/PsuadeUtil.h"
#include "Main/Psuade.h"

#define PABS(x) ((x) > 0 ? (x) : (-(x)))

#ifdef HAVE_MARS
extern "C" 
{
  void mars_process(int,int,double**,double*,float*,int,int,int*,
                 float*,int*);
  void mars_fmod(int,int,double**,double*,float*,int*);
}
#endif

// ************************************************************************
// Constructor for object class Mars
// ------------------------------------------------------------------------
Mars::Mars(int nInputs,int nSamples) : FuncApprox(nInputs,nSamples)
{
#ifdef HAVE_MARS
   int  ss, ii;
   char pString[501];

   faID_ = PSUADE_RS_MARS;

   nBasisFcns_ = 50;
   if (nBasisFcns_ > nSamples) nBasisFcns_ = nSamples - 1;

   if (nInputs >= 8) maxVarPerBasis_ = 8;
   else              maxVarPerBasis_ = nInputs;

   varFlags_ = new int[nInputs];
   for (ii = 0; ii < nInputs; ii++) varFlags_[ii] = 1;
   for (ii = 0; ii < nSamples_; ii++) weights_[ii] = 1000.0;
   chooseWght_ = 0;

   if (psRSExpertMode_ == 1)
   {
      printf("MARS: Current number of basis functions = %d\n", nBasisFcns_);
      sprintf(pString,"Enter the number of basis functions (>=10, <= %d): ",
              nSamples);
      nBasisFcns_ = getInt(10, nSamples, pString);
      printf("MARS: Current degree of interactions    = %d\n",maxVarPerBasis_);
      sprintf(pString, "Enter the degree of interactions (<=%d) : ", nInputs);
      maxVarPerBasis_ = getInt(1, nInputs, pString);
      printf("MARS: current weight of each sample point is set to 1.\n");
      printf("      Other option (1) is to set weight = abs(output)\n");
      sprintf(pString, "Change to option (1) ? (1 means yes, otherwise means no) ");
      chooseWght_ = getInt(0, 1, pString);
   }

   wgts_ = new float[nSamples];
   for (ss = 0; ss < nSamples; ss++) wgts_[ss] = 1.0;
   fm_ = NULL;
   im_ = NULL;
#else
   printf("PSUADE ERROR : Mars not installed.\n");
   exit(1);
#endif
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
Mars::~Mars()
{
   if (wgts_     != NULL) delete [] wgts_;
   if (varFlags_ != NULL) delete [] varFlags_;
   if (fm_       != NULL) delete [] fm_;
   if (im_       != NULL) delete [] im_;
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int Mars::genNDGridData(double *XX, double *Y, int *N, double **XX2, 
                        double **Y2)
{
#ifdef HAVE_MARS
   int    totPts, length, ss, ii, iOne=1;
   double *HX, *Xloc, **X, **X2;
   char   *targv[1], *cStr;

   if (fm_ != NULL) delete [] fm_;
   if (im_ != NULL) delete [] im_;
   length = 3 + 5 * nBasisFcns_ * (5 * maxVarPerBasis_ + nSamples_ + 6) +
            2 * nInputs_ + nSamples_;
   fm_    = new float[length];
   length = 21 + 10 * nBasisFcns_ * (3 * maxVarPerBasis_ + 8);
   im_ = new int[length];
   if (chooseWght_ == 1)
   {
      for (ss = 0; ss < nSamples_; ss++)
      {
         if (PABS(Y[ss]) < 1.0e-12) wgts_[ss] = 1.0;
         else                       wgts_[ss] = (float) PABS(Y[ss]);
      }
   }
   else
   {
      for (ss = 0; ss < nSamples_; ss++) wgts_[ss] = (float) weights_[ss];
   }

   X = new double*[nSamples_];
   for (ss = 0; ss < nSamples_; ss++) X[ss] = new double[nInputs_];
   for (ss = 0; ss < nSamples_; ss++) 
      for (ii = 0; ii < nInputs_; ii++) X[ss][ii] = XX[ss*nInputs_+ii];

   if (outputLevel_ >= 2) printf("Entering Mars (process)\n");
   if (outputLevel_ >= 2) 
      printf("If it crashes here, it is mars_process problem (try different nSamples).\n");
   if (outputLevel_ >= 2) 
      printf("Mars: nBasis, maxVarPerBasis = %d %d\n", 
             nBasisFcns_, maxVarPerBasis_);
   mars_process(nSamples_, nInputs_, X, Y, wgts_, nBasisFcns_,
                maxVarPerBasis_, varFlags_, fm_, im_);
   if (outputLevel_ >= 2) 
      printf("Exit from mars_process.\n");
   if (outputLevel_ >= 2) printf("Returning from Mars\n");
   for (ss = 0; ss < nSamples_; ss++) delete [] X[ss];
   delete [] X;
   if (outputLevel_ >= 4) 
   {
      cStr = new char[10];
      strcpy(cStr, "rank");
      targv[0] = (char *) cStr;
      setParams(iOne, targv);
      delete [] cStr;
   }
   if ((*N) == -999) return 0;
  
   totPts = nPtsPerDim_;
   for (ii = 1; ii < nInputs_; ii++) totPts = totPts * nPtsPerDim_;
   HX = new double[nInputs_];
   for (ii = 0; ii < nInputs_; ii++) 
      HX[ii] = (upperBounds_[ii]-lowerBounds_[ii]) / (double) (nPtsPerDim_-1); 
 
   X2 = new double*[totPts];
   for (ss = 0; ss < totPts; ss++) X2[ss] = new double[nInputs_];
   (*Y2) = new double[totPts];
   Xloc  = new double[nInputs_];
 
   for (ii = 0; ii < nInputs_; ii++) Xloc[ii] = lowerBounds_[ii];
 
   for (ss = 0; ss < totPts; ss++)
   {
      for (ii = 0; ii < nInputs_; ii++ ) X2[ss][ii] = Xloc[ii];
      for (ii = 0; ii < nInputs_; ii++ ) 
      {
         Xloc[ii] += HX[ii];
         if (Xloc[ii] < upperBounds_[ii] || 
              PABS(Xloc[ii] - upperBounds_[ii]) < 1.0E-7) break;
         else Xloc[ii] = lowerBounds_[ii];
      }
   }
 
   if (outputLevel_ >= 2) printf("Entering Mars (fmod)\n");
   if (outputLevel_ >= 2) 
      printf("If it crashes here, it is mars_fmod problem.\n");
   mars_fmod(totPts, nInputs_, X2, *Y2, fm_, im_);
   if (outputLevel_ >= 2) 
      printf("Exit from mars_fmod.\n");
   if (outputLevel_ >= 2) printf("Returning from Mars\n");
   (*N) = totPts;
   delete [] Xloc;

   (*XX2) = new double[nInputs_ * totPts];
   for (ss = 0; ss < totPts; ss++) 
      for (ii = 0; ii < nInputs_; ii++) (*XX2)[ss*nInputs_+ii] = X2[ss][ii];
   for (ss = 0; ss < totPts; ss++ ) delete [] X2[ss];
   delete [] X2;
   delete [] HX;
#else
   printf("PSUADE ERROR : Mars not used.\n");
   return -1;
#endif
   return 0;
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int Mars::gen1DGridData(double *X, double *Y, int ind1, 
                        double *settings, int *N, double **XX, double **YY)
{
#ifdef HAVE_MARS
   int    ss, ii, totPts, length;
   double HX, **X2, **X1;

   if (fm_ != NULL) delete [] fm_;
   if (im_ != NULL) delete [] im_;
   length = 3 + 5 * nBasisFcns_ * (5 * maxVarPerBasis_ + nSamples_ + 6) +
            2 * nInputs_ + nSamples_;
   fm_    = new float[length];
   length = 21 + 10 * nBasisFcns_ * (3 * maxVarPerBasis_ + 8);
   im_ = new int[length];
   if (chooseWght_ == 1)
   {
      for (ss = 0; ss < nSamples_; ss++)
      {
         if (PABS(Y[ss]) < 1.0e-12) wgts_[ss] = 1.0;
         else                       wgts_[ss] = (float) PABS(Y[ss]);
      }
   }
   else
   {
      for (ss = 0; ss < nSamples_; ss++) wgts_[ss] = (float) weights_[ss];
   }

   X1 = new double*[nSamples_];
   for (ss = 0; ss < nSamples_; ss++) X1[ss] = new double[nInputs_];
   for (ss = 0; ss < nSamples_; ss++) 
      for (ii = 0; ii < nInputs_; ii++) X1[ss][ii] = X[ss*nInputs_+ii];
   if (outputLevel_ >= 2) printf("Entering Mars (process)\n");
   if (outputLevel_ >= 2) 
      printf("If it crashes here, it is mars_process problem (try different nSamples).\n");
   mars_process(nSamples_, nInputs_, X1, Y, wgts_, nBasisFcns_,
                maxVarPerBasis_, varFlags_, fm_, im_);
   if (outputLevel_ >= 2) 
      printf("Exit from mars_process.\n");
   if (outputLevel_ >= 2) printf("Returning from Mars\n");
   for (ss = 0; ss < nSamples_; ss++) delete [] X1[ss];
   delete [] X1;

   totPts = nPtsPerDim_;
   HX = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 

   X2 = new double*[totPts];
   for (ss = 0; ss < totPts; ss++) 
   {
      X2[ss] = new double[nInputs_];
      for (ii = 0; ii < nInputs_; ii++) X2[ss][ii] = settings[ii]; 
   }
    
   (*XX) = new double[totPts];
   (*YY) = new double[totPts];

   for (ii = 0; ii < totPts; ii++) 
   {
      X2[ii][ind1]  = HX * ii + lowerBounds_[ind1];
      (*XX)[ii]     = HX * ii + lowerBounds_[ind1];
      (*YY)[ii]     = 0.0;
   }

   if (outputLevel_ >= 2) printf("Entering Mars (fmod)\n");
   if (outputLevel_ >= 2) 
      printf("If it crashes here, it is mars_fmod problem.\n");
   mars_fmod(totPts, nInputs_, X2, *YY, fm_, im_);
   if (outputLevel_ >= 2) 
      printf("Exit from mars_fmod.\n");
   if (outputLevel_ >= 2) printf("Returning from Mars\n");
   (*N) = totPts;
   for (ss = 0; ss < totPts; ss++) delete [] X2[ss];
   delete [] X2;
#else
   printf("PSUADE ERROR : Mars not used.\n");
   return -1;
#endif
   return 0;
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int Mars::gen2DGridData(double *X, double *Y, int ind1, int ind2, 
                        double *settings, int *N, double **XX, double **YY)
{
#ifdef HAVE_MARS
   int    ss, ii, jj, totPts, index, length;
   double *HX, **X2, **X1;

   if (fm_ != NULL) delete [] fm_;
   if (im_ != NULL) delete [] im_;
   length = 3 + 5 * nBasisFcns_ * (5 * maxVarPerBasis_ + nSamples_ + 6) +
            2 * nInputs_ + nSamples_;
   fm_    = new float[length];
   length = 21 + 10 * nBasisFcns_ * (3 * maxVarPerBasis_ + 8);
   im_ = new int[length];
   if (chooseWght_ == 1)
   {
      for (ss = 0; ss < nSamples_; ss++)
      {
         if (PABS(Y[ss]) < 1.0e-12) wgts_[ss] = 1.0;
         else                       wgts_[ss] = (float) PABS(Y[ss]);
      }
   }
   else
   {
      for (ss = 0; ss < nSamples_; ss++) wgts_[ss] = (float) weights_[ss];
   }

   X1 = new double*[nSamples_];
   for (ss = 0; ss < nSamples_; ss++) X1[ss] = new double[nInputs_];
   for (ss = 0; ss < nSamples_; ss++) 
      for (ii = 0; ii < nInputs_; ii++) X1[ss][ii] = X[ss*nInputs_+ii];
   if (outputLevel_ >= 2) printf("Entering Mars (process)\n");
   if (outputLevel_ >= 2) 
      printf("If it crashes here, it is mars_process problem (try different nSamples).\n");
   mars_process(nSamples_, nInputs_, X1, Y, wgts_, nBasisFcns_,
                maxVarPerBasis_, varFlags_, fm_, im_);
   if (outputLevel_ >= 2) 
      printf("Exit from mars_process.\n");
   if (outputLevel_ >= 2) printf("Returning from Mars\n");
   for (ss = 0; ss < nSamples_; ss++) delete [] X1[ss];
   delete [] X1;

   totPts = nPtsPerDim_ * nPtsPerDim_;
   HX = new double[2];
   HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
   HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 

   X2 = new double*[totPts];
   for (ss = 0; ss < totPts; ss++) 
   {
      X2[ss] = new double[nInputs_];
      for (ii = 0; ii < nInputs_; ii++) X2[ss][ii] = settings[ii]; 
   }
    
   (*XX) = new double[2*totPts];
   (*YY) = new double[totPts];

   for (ii = 0; ii < nPtsPerDim_; ii++) 
   {
      for (jj = 0; jj < nPtsPerDim_; jj++)
      {
         index = ii * nPtsPerDim_ + jj;
         X2[index][ind1]  = HX[0] * ii + lowerBounds_[ind1];
         X2[index][ind2]  = HX[1] * jj + lowerBounds_[ind2];
         (*XX)[index*2]   = HX[0] * ii + lowerBounds_[ind1];
         (*XX)[index*2+1] = HX[1] * jj + lowerBounds_[ind2];
      }
   }

   if (outputLevel_ >= 2) printf("Entering Mars (fmod)\n");
   if (outputLevel_ >= 2) 
      printf("If it crashes here, it is mars_fmod problem.\n");
   mars_fmod(totPts, nInputs_, X2, *YY, fm_, im_);
   if (outputLevel_ >= 2) 
      printf("Exit from mars_fmod.\n");
   if (outputLevel_ >= 2) printf("Returning from Mars\n");
   (*N) = totPts;
   for (ss = 0; ss < totPts; ss++) delete [] X2[ss];
   delete [] X2;
   delete [] HX;
#else
   printf("PSUADE ERROR : Mars not used.\n");
   return -1;
#endif
   return 0;
}

// ************************************************************************
// Generate 3D results for display
// ------------------------------------------------------------------------
int Mars::gen3DGridData(double *X, double *Y, int ind1, int ind2, int ind3,
                        double *settings, int *N, double **XX, double **YY)
{
#ifdef HAVE_MARS
   int    ss, ii, jj, ll, totPts, index, length;
   double *HX, **X2, **X1;

   if (fm_ != NULL) delete [] fm_;
   if (im_ != NULL) delete [] im_;
   length = 3 + 5 * nBasisFcns_ * (5 * maxVarPerBasis_ + nSamples_ + 6) +
            2 * nInputs_ + nSamples_;
   fm_    = new float[length];
   length = 21 + 10 * nBasisFcns_ * (3 * maxVarPerBasis_ + 8);
   im_ = new int[length];
   if (chooseWght_ == 1)
   {
      for (ss = 0; ss < nSamples_; ss++)
      {
         if (PABS(Y[ss]) < 1.0e-12) wgts_[ss] = 1.0;
         else                       wgts_[ss] = (float) PABS(Y[ss]);
      }
   }
   else
   {
      for (ss = 0; ss < nSamples_; ss++) wgts_[ss] = (float) weights_[ss];
   }

   X1 = new double*[nSamples_];
   for (ss = 0; ss < nSamples_; ss++) X1[ss] = new double[nInputs_];
   for (ss = 0; ss < nSamples_; ss++) 
      for (ii = 0; ii < nInputs_; ii++) X1[ss][ii] = X[ss*nInputs_+ii];
   if (outputLevel_ >= 2) printf("Entering Mars (process)\n");
   if (outputLevel_ >= 2) 
      printf("If it crashes here, it is mars_process problem (try different nSamples).\n");
   mars_process(nSamples_, nInputs_, X1, Y, wgts_, nBasisFcns_,
                maxVarPerBasis_, varFlags_, fm_, im_);
   if (outputLevel_ >= 2) 
      printf("Exit from mars_process.\n");
   if (outputLevel_ >= 2) printf("Returning from Mars\n");
   for (ss = 0; ss < nSamples_; ss++) delete [] X1[ss];
   delete [] X1;

   totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
   HX = new double[3];
   HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
   HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 
   HX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_ - 1); 

   X2 = new double*[totPts];
   for (ss = 0; ss < totPts; ss++) 
   {
      X2[ss] = new double[nInputs_];
      for (ii = 0; ii < nInputs_; ii++) X2[ss][ii] = settings[ii]; 
   }
    
   (*XX) = new double[3*totPts];
   (*YY) = new double[totPts];

   for (ii = 0; ii < nPtsPerDim_; ii++) 
   {
      for (jj = 0; jj < nPtsPerDim_; jj++)
      {
         for (ll = 0; ll < nPtsPerDim_; ll++)
         {
            index = ii * nPtsPerDim_ * nPtsPerDim_ + jj * nPtsPerDim_ + ll;
            X2[index][ind1]  = HX[0] * ii + lowerBounds_[ind1];
            X2[index][ind2]  = HX[1] * jj + lowerBounds_[ind2];
            X2[index][ind3]  = HX[2] * ll + lowerBounds_[ind3];
            (*XX)[index*3]   = HX[0] * ii + lowerBounds_[ind1];
            (*XX)[index*3+1] = HX[1] * jj + lowerBounds_[ind2];
            (*XX)[index*3+2] = HX[2] * ll + lowerBounds_[ind3];
         }
      }
   }

   if (outputLevel_ >= 2) printf("Entering Mars (fmod)\n");
   if (outputLevel_ >= 2) 
      printf("If it crashes here, it is mars_fmod problem.\n");
   mars_fmod(totPts, nInputs_, X2, *YY, fm_, im_);
   if (outputLevel_ >= 2) 
      printf("Exit from mars_fmod.\n");
   if (outputLevel_ >= 2) printf("Returning from Mars\n");
   (*N) = totPts;
   for (ss = 0; ss < totPts; ss++) delete [] X2[ss];
   delete [] X2;
   delete [] HX;
#else
   printf("PSUADE ERROR : Mars not used.\n");
   return -1;
#endif
   return 0;
}

// ************************************************************************
// Generate 4D results for display
// ------------------------------------------------------------------------
int Mars::gen4DGridData(double *X, double *Y, int ind1, int ind2, int ind3,
                        int ind4, double *settings, int *N, double **XX, 
                        double **YY)
{
#ifdef HAVE_MARS
   int    ss, ii, jj, ll, mm, totPts, index, length;
   double *HX, **X2, **X1;

   if (fm_ != NULL) delete [] fm_;
   if (im_ != NULL) delete [] im_;
   length = 3 + 5 * nBasisFcns_ * (5 * maxVarPerBasis_ + nSamples_ + 6) +
            2 * nInputs_ + nSamples_;
   fm_    = new float[length];
   length = 21 + 10 * nBasisFcns_ * (3 * maxVarPerBasis_ + 8);
   im_ = new int[length];
   if (chooseWght_ == 1)
   {
      for (ss = 0; ss < nSamples_; ss++)
      {
         if (PABS(Y[ss]) < 1.0e-12) wgts_[ss] = 1.0;
         else                       wgts_[ss] = (float) PABS(Y[ss]);
      }
   }
   else
   {
      for (ss = 0; ss < nSamples_; ss++) wgts_[ss] = (float) weights_[ss];
   }

   X1 = new double*[nSamples_];
   for (ss = 0; ss < nSamples_; ss++) X1[ss] = new double[nInputs_];
   for (ss = 0; ss < nSamples_; ss++) 
      for (ii = 0; ii < nInputs_; ii++) X1[ss][ii] = X[ss*nInputs_+ii];
   if (outputLevel_ >= 2) printf("Entering Mars (process)\n");
   if (outputLevel_ >= 2) 
      printf("If it crashes here, it is mars_process problem (try different nSamples).\n");
   mars_process(nSamples_, nInputs_, X1, Y, wgts_, nBasisFcns_,
                maxVarPerBasis_, varFlags_, fm_, im_);
   if (outputLevel_ >= 2) 
      printf("Exit from mars_process.\n");
   if (outputLevel_ >= 2) printf("Returning from Mars\n");
   for (ss = 0; ss < nSamples_; ss++) delete [] X1[ss];
   delete [] X1;

   totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
   HX = new double[4];
   HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
   HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 
   HX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_ - 1); 
   HX[3] = (upperBounds_[ind4] - lowerBounds_[ind4]) / (nPtsPerDim_ - 1); 

   X2 = new double*[totPts];
   for (ss = 0; ss < totPts; ss++) 
   {
      X2[ss] = new double[nInputs_];
      for (ii = 0; ii < nInputs_; ii++) X2[ss][ii] = settings[ii]; 
   }
    
   (*XX) = new double[4*totPts];
   (*YY) = new double[totPts];

   for (ii = 0; ii < nPtsPerDim_; ii++) 
   {
      for (jj = 0; jj < nPtsPerDim_; jj++)
      {
         for (ll = 0; ll < nPtsPerDim_; ll++)
         {
            for (mm = 0; mm < nPtsPerDim_; mm++)
            {
               index = ii*nPtsPerDim_*nPtsPerDim_ * nPtsPerDim_ +
                       jj*nPtsPerDim_*nPtsPerDim_ + ll*nPtsPerDim_ + mm;
               X2[index][ind1]  = HX[0] * ii + lowerBounds_[ind1];
               X2[index][ind2]  = HX[1] * jj + lowerBounds_[ind2];
               X2[index][ind3]  = HX[2] * ll + lowerBounds_[ind3];
               X2[index][ind4]  = HX[3] * mm + lowerBounds_[ind4];
               (*XX)[index*4]   = HX[0] * ii + lowerBounds_[ind1];
               (*XX)[index*4+1] = HX[1] * jj + lowerBounds_[ind2];
               (*XX)[index*4+2] = HX[2] * ll + lowerBounds_[ind3];
               (*XX)[index*4+3] = HX[3] * mm + lowerBounds_[ind4];
            }
         }
      }
   }

   if (outputLevel_ >= 2) printf("Entering Mars (fmod)\n");
   if (outputLevel_ >= 2) 
      printf("If it crashes here, it is mars_fmod problem.\n");
   mars_fmod(totPts, nInputs_, X2, *YY, fm_, im_);
   if (outputLevel_ >= 2) 
      printf("Exit from mars_fmod.\n");
   if (outputLevel_ >= 2) printf("Returning from Mars\n");
   (*N) = totPts;
   for (ss = 0; ss < totPts; ss++) delete [] X2[ss];
   delete [] X2;
   delete [] HX;
#else
   printf("PSUADE ERROR : Mars not used.\n");
   return -1;
#endif
   return 0;
}

// ************************************************************************
// Generate results to a file called psuade_grid_data
// ------------------------------------------------------------------------
int Mars::writeToFileGrid2DData(double *XX, double *Y, int ind1,
                                int ind2, double *settings)
{
#ifdef HAVE_MARS
   int    ss, ii, jj, totPts, index;
   double *HX, **X2, dtmp, *YY, **X;
   FILE   *fp;

   X = new double*[nSamples_];
   for (ss = 0; ss < nSamples_; ss++) X[ss] = new double[nInputs_];
   for (ss = 0; ss < nSamples_; ss++) 
      for (ii = 0; ii < nInputs_; ii++) X[ss][ii] = XX[nInputs_*ss+ii];

   if (chooseWght_ == 1)
   {
      for (ss = 0; ss < nSamples_; ss++)
      {
         if (PABS(Y[ss]) < 1.0e-12) wgts_[ss] = 1.0;
         else                       wgts_[ss] = (float) PABS(Y[ss]);
      }
   }
   else
   {
      for (ss = 0; ss < nSamples_; ss++) wgts_[ss] = (float) weights_[ss];
   }

   if (outputLevel_ >= 2) printf("Entering Mars (process)\n");
   if (outputLevel_ >= 2)
      printf("If it crashes here, it is mars_process problem (try different nSamples).\n");
   mars_process(nSamples_, nInputs_, X, Y, wgts_, nBasisFcns_,
                maxVarPerBasis_, varFlags_, fm_, im_);
   if (outputLevel_ >= 2) printf("Returning from Mars\n");
   for (ss = 0; ss < nSamples_; ss++) delete [] X[ss];
   delete [] X;

   totPts = nPtsPerDim_ * nPtsPerDim_;
   HX = new double[2];
   HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
   HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 

   X2 = new double*[totPts];
   for (ss = 0; ss < totPts; ss++) 
   {
      X2[ss] = new double[nInputs_];
      for (ii = 0; ii < nInputs_; ii++) X2[ss][ii] = settings[ii]; 
   }
   YY = new double[2 * totPts];

   for (ii = 0; ii < nPtsPerDim_; ii++) 
   {
      for (jj = 0; jj < nPtsPerDim_; jj++) 
      {
        index = ii * nPtsPerDim_ + jj;
        X2[index][ind1] = HX[0] * ii + lowerBounds_[ind1];
        X2[index][ind2] = HX[1] * jj + lowerBounds_[ind2];
      }
   }

   if (outputLevel_ >= 2) printf("Entering Mars (fmod)\n");
   if (outputLevel_ >= 2)
      printf("If it crashes here, it is mars_fmod problem.\n");
   mars_fmod(totPts, nInputs_, X2, YY, fm_, im_);
   if (outputLevel_ >= 2) printf("Returning from Mars\n");

   fp = fopen("psuade_grid_data", "w");
   fprintf(fp, "%d\n",nPtsPerDim_);
   dtmp = lowerBounds_[ind1];
   for (ii = 0; ii < nPtsPerDim_; ii++) 
   {
      fprintf(fp, "%e\n", dtmp);
      dtmp += HX[0];
   }
   fprintf(fp, "%d\n",nPtsPerDim_);
   dtmp = lowerBounds_[ind2];
   for (ii = 0; ii < nPtsPerDim_; ii++) 
   {
      fprintf(fp, "%e\n", dtmp);
      dtmp += HX[1];
   }
   for (ii = 0; ii < nPtsPerDim_; ii++)
      for (jj = 0; jj < nPtsPerDim_; jj++)
         fprintf(fp, "%e\n", YY[ii*nPtsPerDim_+jj]);
   fclose(fp);

   for (ii = 0; ii < totPts; ii++) delete [] X2[ii];
   delete [] X2;
   delete [] YY;
   delete [] HX;
#else
   printf("PSUADE ERROR : Mars not used.\n");
#endif
   return 0;
}

// ************************************************************************
// Evaluate a given point 
// ------------------------------------------------------------------------
double Mars::evaluatePoint(double *X)
{
   double Y=0.0;
#ifdef HAVE_MARS
   int    ii;
   double **X2;

   X2    = new double*[1];
   X2[0] = new double[nInputs_];
   for (ii = 0; ii < nInputs_; ii++) X2[0][ii] = X[ii]; 
   mars_fmod(1, nInputs_, X2, &Y, fm_, im_);
   delete [] X2[0];
   delete [] X2;
#else
   printf("PSUADE ERROR : Mars not used.\n");
#endif
   return Y;
}

// ************************************************************************
// Evaluate a number of points
// ------------------------------------------------------------------------
double Mars::evaluatePoint(int npts, double *X, double *Y)
{
#ifdef HAVE_MARS
   int    ii, kk;
   double **X2;

   X2 = new double*[npts];
   for (kk = 0; kk < npts; kk++)
   {
      X2[kk] = new double[nInputs_];
      for (ii = 0; ii < nInputs_; ii++) X2[kk][ii] = X[kk*nInputs_+ii]; 
   }
   mars_fmod(npts, nInputs_, X2, Y, fm_, im_);
   for (kk = 0; kk < npts; kk++) delete [] X2[kk];
   delete [] X2;
#else
   printf("PSUADE ERROR : Mars not used.\n");
#endif
   return 0.0;
}

// ************************************************************************
// Evaluate a given point with standard deviation
// ------------------------------------------------------------------------
double Mars::evaluatePointFuzzy(double *X, double &std)
{
   double Y=0.0;
#ifdef HAVE_MARS
   int    ii;
   double **X2;

   X2    = new double*[1];
   X2[0] = new double[nInputs_];
   for (ii = 0; ii < nInputs_; ii++) X2[0][ii] = X[ii]; 
   mars_fmod(1, nInputs_, X2, &Y, fm_, im_);
   delete [] X2[0];
   delete [] X2;
#else
   printf("PSUADE ERROR : Mars not used.\n");
#endif
   std = 0.0;
   return Y;
}

// ************************************************************************
// Evaluate a number of points with standard deviations
// ------------------------------------------------------------------------
double Mars::evaluatePointFuzzy(int npts,double *X, double *Y,double *Ystd)
{
#ifdef HAVE_MARS
   int    ii, kk;
   double **X2;

   X2 = new double*[npts];
   for (kk = 0; kk < npts; kk++)
   {
      X2[kk] = new double[nInputs_];
      for (ii = 0; ii < nInputs_; ii++) X2[kk][ii] = X[kk*nInputs_+ii]; 
      Ystd[kk] = 0.0;
   }
   mars_fmod(npts, nInputs_, X2, Y, fm_, im_);
   for (kk = 0; kk < npts; kk++) delete [] X2[kk];
   delete [] X2;
#else
   printf("PSUADE ERROR : Mars not used.\n");
#endif
   return 0.0;
}

// ************************************************************************
// print MARS output information
// ------------------------------------------------------------------------
void Mars::printStat()
{
   int    lineLeng=500, i1, i2, nk=0, ind1, ind2;
   double *coefs;
   char   line[500], word1[100], word2[100], word3[100];
   FILE   *fp;
#if 0
   int    nGroups, *groupMembers, ivars, parent;
   double *vars, *knots, dtmp1, dtmp2, dtmp3, dtmp4, dtmp5, dtmp6;
#endif

   // first get the nk parameter
   fp = fopen(".psuade_mars", "r");
   if (fp == NULL) 
   {
      printf("Mars importance information not available.\n");
      return;
   }
   while ((fgets(line, lineLeng, fp) != NULL) && (feof(fp) == 0))
   {
      sscanf(line,"%s %s", word1, word2);
      if (!strcmp(word1,"input") && !strcmp(word2,"parameters"))
      {
         fgets(line, lineLeng, fp);
         fgets(line, lineLeng, fp);
         sscanf(line,"%d %d %d", &i1, &i2, &nk);
         break;
      }
   }
   fclose(fp);
   if (nk <= 0) 
   {
      printf("Mars::printStat ERROR - number of basis functions = %d <= 0.\n",nk);
      return;
   }

   coefs = new double[nk+1];
   for (i1 = 0; i1 <= nk; i1++) coefs[i1] = -9999.0;
   while ((fgets(line, lineLeng, fp) != NULL) && (feof(fp) == 0))
   {
      sscanf(line,"%s %s", word1, word2);
      if (!strcmp(word1,"final") && !strcmp(word2,"model"))
      {
         for (i1 = 0; i1 <= nk; i1+=6)
         {
            fgets(line, lineLeng, fp);
            fscanf(fp, "%s", word3);
            if (strcmp(word3, "bsfn:")) 
            {
               nk = i1 - 1;
               break;
            }
            else
            {
               for (i2 = 0; i2 < 6; i2++)
               {
                  fscanf(fp, "%d", &ind1);
                  if (ind1 != (i1 + i2))
                  {
                     printf("Mars : nk reset from %d to %d\n", nk, i1+i2-1);
                     nk = i1 + i2 - 1; 
                     break;
                  }
               }
               if (i2 == 6) fgets(line, lineLeng, fp);
            } 
            fscanf(fp, "%s", word3);
            if (strcmp(word3, "coef:")) 
            {
               printf("Mars::printStat WARNING reading coefficients (%s).\n",
                      word3);
               printf("%s\n", line);
            }
            for (i2 = 0; i2 < 6; i2++) 
               if (i1+i2 <= nk) fscanf(fp, "%lg", &coefs[i1+i2]);
            fgets(line, lineLeng, fp);
         }
         break;
      }
   }
   fclose(fp);
   if (coefs[nk] == -9999.0)
   {
      printf("Mars::printStat ERROR - invalid Mars coefficients.\n");
      delete [] coefs;
      return;
   }
 
#if 0
   fp = fopen(".psuade_mars", "r");
   while ((fgets(line, lineLeng, fp) != NULL) && (feof(fp) == 0))
   {
      sscanf(line,"%s %s %s", word1, word2, word3);
      if (!strcmp(word1,"forward") && !strcmp(word2,"stepwise"))
      {
         nGroups = nk;
         groupMembers = new int*[nk];
         for (i1 = 0; i1 < nk; i1++) 
         {
            groupMembers[i1] = new int[2];
            for (i2 = 0; i2 < 2; i2++) groupMembers[i1][i2] = -1; 
         }
         knots = new double[nk+1];
         vars = new double[nk+1];
         ivars = new int[nk+1];
         parent = new int[nk+1];
         fgets(line, lineLeng, fp);
         fgets(line, lineLeng, fp);
         ind1 = 0;
         nGroups = 0;
         while (ind1 < nk)
         {
            fgets(line, lineLeng, fp);
            sscanf(line, "%d %lg", &i1, &dtmp1);
            i2 = (int) dtmp1;
            if ((dtmp1 - (double) i2) != 0) i2 = -1;
            if ((i1 - i2) != 1)
            {
               if (i1 == 0)
               {
                  sscanf(line, "%d %lg %lg %lg",&groupMembers[nGroups][0], 
                         &dtmp1, &dtmp2, &dtmp3); 
                  ivars[i1] = 0;
                  knots[i1] = 0.0;
                  parent[i1] = 0;
               }
               else
               {
                  sscanf(line, "%d %lg %lg %lg %lg %lg %lg", 
                      &groupMembers[nGroups][0], &dtmp1, &dtmp2, &dtmp3, 
                      &dtmp4, &dtmp5, &dtmp6);
                  ivars[i1] = (int) dtmp4;
                  knots[i1] = dtmp5;
                  parent[i1] = dtmp6;
               }
               nGroups++;
            }
            else
            {   
               sscanf(line, "%d %d %lg %lg %lg %lg %lg %lg", 
                      &groupMembers[nGroups][0], &groupMembers[nGroups][1],
                      &dtmp1, &dtmp2, &dtmp3, &dtmp4, &dtmp5, &dtmp6);
               ivars[i1] = ivars[ind1+1] = (int) dtmp4;
               knots[i1] = knots[ind1+1] = dtmp5;
               parent[i1] = parent[ind1+1] = (int) dtmp6;
               nGroups++;
            }
            ind1 = groupMembers[nGroups-1][0];
         }

         // analyze and print out

         for (i1 = 0; i1 < nGroups; i1++)
         {
            ind1 = groupMembers[i1][0];
            ind2 = groupMembers[i1][1];
            if (parent[ind1] == 0)
            {
               if (ind2 > 0 && coefs[ind2] != 0.0)
               {
                  printf("BF%4d : variable %3d > %9.3e has coefficient = %9.3e\n",
                         ind2, ivars[ind2], knots[ind2], coefs[ind2]);
               }
               if (ind1 > 0 && ind2 > 0 && coefs[ind1] != 0.0)
               {
                  printf("BF%4d : variable %3d < %9.3e has coefficient = %9.3e\n",
                         ind1, ivars[ind1], knots[ind1], coefs[ind1]);
               }
               else if (ind1 > 0 && ind2 < 0 && coefs[ind1] != 0.0)
               {
                  printf("BF%4d : variable %3d > %9.3e has coefficient = %9.3e\n",
                         ind1, ivars[ind1], knots[ind1], coefs[ind1]);
               }
            }
         }
         for (i1 = 0; i1 < nGroups; i1++)
         {
            ind1 = groupMembers[i1][0];
            if (parent[ind1] != 0)
            {
               if (coefs[ind1] != 0.0)
                  printf("Interaction: variables (%3d, %3d) with strength = %9.3e\n",
                         ivars[ind1], ivars[parent[ind1]], PABS(coefs[ind1]));
               ind1 = groupMembers[i1][1];
               if (ind1 > 0 && coefs[ind1] != 0.0)
                  printf("Interaction: variables (%3d, %3d) with strength = %9.3e\n",
                         ivars[ind1], ivars[parent[ind1]], PABS(coefs[ind1]));
            }
         }
         for (i1 = 0; i1 < nk; i1++) delete [] groupMembers[i1];
         delete [] groupMembers;
         delete [] knots;
         delete [] vars;
         delete [] ivars;
         delete [] parent;
         break;
      }
   }
   delete [] coefs;
   fclose(fp);
#endif

   fp = fopen(".psuade_mars", "r");
   while ((fgets(line, lineLeng, fp) != NULL) && (feof(fp) == 0))
   {
      sscanf(line,"%s %s", word1, word2);
      if (!strcmp(word1,"anova") && !strcmp(word2,"decomposition"))
      {
         printf("\n Mars ANOVA Decomposition\n");
         fgets(line, lineLeng, fp);
         printf("%s", line);
         fgets(line, lineLeng, fp);
         ind2 = 1;
         sscanf(line, "%d", &ind1);
         while (ind2 == ind1)
         {
            printf("%s", line);
            fgets(line, lineLeng, fp);
            ind1 = -1;
            sscanf(line, "%d", &ind1);
            ind2++;
         }
         break;
      }
   }
   fclose(fp);
}

// ************************************************************************
// set parameters
// ------------------------------------------------------------------------
double Mars::setParams(int targc, char **targv)
{
   int    lineLeng=500, *iArray, ii, jj, nCount;
   double *dArray, dmax;
   char   line[500], word1[500], word2[500], word3[500];
   FILE   *fp, *fp2;

   if (targc > 0 && !strcmp(targv[0], "rank"))
   {
      fp = fopen(".psuade_mars", "r");
      if (fp != NULL)
      {
         nCount = 0;
         strcpy(word1, "none");
         dArray = new double[nInputs_];
         iArray = new int[nInputs_];
         while ((fgets(line, lineLeng, fp) != NULL) && (feof(fp) == 0))
         {
            sscanf(line,"%s %s %s", word1, word2, word3);
            if (!strcmp(word1,"relative") && !strcmp(word2,"variable") && 
                !strcmp(word3,"importance:"))
            {
               fgets(line, lineLeng, fp);
               fgets(line, lineLeng, fp);
               if (feof(fp) == 0)
               {
                  for (ii = 0; ii < nInputs_; ii+=6)
                  {
                     for (jj = 0; jj < 6; jj++)
                     {
                        if (ii+jj < nInputs_)
                        {
                           fscanf(fp,"%lg", &(dArray[ii+jj]));
                           nCount++;
                        }
                     }
                     fgets(line, lineLeng, fp);
                     fgets(line, lineLeng, fp);
                     fgets(line, lineLeng, fp);
                  }
               }
            }
         }
         if (nCount != nInputs_)
         {
             printf("* Mars gives insufficient importance measures.\n");
             printf("* To inspect Mars output, examine .psuade_mars.\n");
         }
         else
         {
            dmax = dArray[0];
            for (ii = 1; ii < nInputs_; ii++) 
               if (dArray[ii] > dmax) dmax = dArray[ii];
            for (ii = 0; ii < nInputs_; ii++) dArray[ii] /= (dmax / 100.0);

            fp2 = fopen("matlabmars.m", "w");
            fwritePlotCLF(fp2);
            fprintf(fp2, "A = [\n");
            for (ii = 0; ii < nInputs_; ii++)
               fprintf(fp2, "%e\n", 0.01 * dArray[ii]);
            fprintf(fp2, "];\n");
            fprintf(fp2, "bar(A, 0.8);\n");
            fwritePlotAxes(fp2);
            fwritePlotTitle(fp2, "MARS Rankings");
            fwritePlotXLabel(fp2, "Input parameters");
            fwritePlotYLabel(fp2, "MARS ranks");
            fclose(fp2);
            printf("MARS ranking is now in matlabmars.m.\n");

            for (ii = 0; ii < nInputs_; ii++) iArray[ii] = ii;
            sortDbleList2a(nInputs_, dArray, iArray);
            printAsterisks(0);
            printf("* Mars screening rankings \n");
            printAsterisks(0);
            for (ii = nInputs_-1; ii >= 0; ii--)
               printf("*  Rank %3d : Input = %3d (score = %4.1f)\n", 
                      nInputs_-ii, iArray[ii]+1, dArray[ii]);
            printAsterisks(0);

         }
         delete [] dArray;
         delete [] iArray;
         fclose(fp);
      }
      else printf("Mars: cannot rank - .psuade_mars file not found.\n");
   }
   else if (targc == 3 && !strcmp(targv[0], "mars_params"))
   {
      nBasisFcns_ = *(int *) targv[1];
      if (nBasisFcns_ < 20) nBasisFcns_ = 20;
      maxVarPerBasis_ = *(int *) targv[2];
      if (maxVarPerBasis_ < 1) maxVarPerBasis_ = 1;
      if (psRSExpertMode_ == 1)
      {
         printf("MARS: numBasis    set to = %d.\n", nBasisFcns_);
         printf("MARS: varPerBasis set to = %d.\n", maxVarPerBasis_);
      }
   }
   else 
   {
      printf("Mars setParams ERROR: invalid command.\n");
   }
   return 0.0;
}

