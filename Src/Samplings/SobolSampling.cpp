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
// Functions for the Sobol's one-at-a-time class 
// AUTHOR : CHARLES TONG
// DATE   : 2004
// ************************************************************************
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "SobolSampling.h"
#include "Sampling.h"
#include "MCSampling.h"
#include "LHSampling.h"

// fix delta h for Sobol instead of random
//#define PSUADE_SAL_GRID

// ************************************************************************
// Constructor
// ------------------------------------------------------------------------
SobolSampling::SobolSampling() : Sampling()
{
   samplingID_ = PSUADE_SAMP_SOBOL;
}

// *******************************************************************
// destructor 
// -------------------------------------------------------------------
SobolSampling::~SobolSampling()
{
}

// ************************************************************************
// initialize the sampling data
// ------------------------------------------------------------------------
int SobolSampling::initialize(int initLevel)
{
   int    iD, iD2, inputID, nReps, sampleCount, iR;
   double **M1Mat, **M2Mat, *ranges, ddata;

   if (nInputs_ == 0 || lowerBounds_ == NULL || upperBounds_ == NULL)
   {
      printf("SobolSampling::initialize ERROR - input not set up.\n");
      exit(1);
   }

   deleteSampleData();
   if (nInputs_ > 25)
   {
      printf("SobolSampling : can only handle nInputs <= 25\n");
      exit(1);
   }
   if (nSamples_ / (nInputs_ + 2) * (nInputs_ + 2) != nSamples_) 
   {
      printf("SobolSampling : nSamples must be multiples of nInputs+2.\n");
      exit(1);
   }
   if (initLevel != 0) return 0;
   allocSampleData();
   nReps = nSamples_ / (nInputs_ + 2);

   if (printLevel_ > 4)
   {
      printf("SobolSampling::initialize: nSamples = %d\n", nSamples_);
      printf("SobolSampling::initialize: nInputs  = %d\n", nInputs_);
      printf("SobolSampling::initialize: nOutputs = %d\n", nOutputs_);
      for (inputID = 0; inputID < nInputs_; inputID++)
         printf("    SobolSampling input %3d = [%e %e]\n", inputID+1,
                lowerBounds_[inputID], upperBounds_[inputID]);
   }

   M1Mat = new double*[nReps];
   for (iD = 0; iD < nReps; iD++) M1Mat[iD] = new double[nInputs_];
   M2Mat = new double*[nReps];
   for (iD = 0; iD < nReps; iD++) M2Mat[iD] = new double[nInputs_];

   ranges = new double[nInputs_];
   for (iD = 0;  iD < nInputs_;  iD++) 
      ranges[iD] = upperBounds_[iD] - lowerBounds_[iD];

   sampleCount = 0;
   generate(M2Mat, nReps);
#ifdef PSUADE_SAL_GRID
   ddata  = 1.0 / (double) (nReps/2 - 1);
   for (iD = 0; iD < nReps; iD++)
   {
      for (iD2 = 0; iD2 < nInputs_; iD2++)
      {
         ddata2 = PSUADE_drand();
         if (ddata2 > 0.5) ddata2 = 1.0;
         else              ddata2 = -1.0;
         if ((M2Mat[iD][iD2]+ddata2*ddata)>1.0 || (M2Mat[iD][iD2]+ddata2*ddata)< 0.0)
              M1Mat[iD][iD2] = M2Mat[iD][iD2] - ddata2*ddata;
         else M1Mat[iD][iD2] = M2Mat[iD][iD2] + ddata2*ddata;
      }
   }
#endif
#ifdef PSUADE_SAL_GRID
   for (iD = 0; iD < nReps; iD++)
      for (iD2 = 0; iD2 < nInputs_; iD2++)
         M1Mat[iD][iD2] = M2Mat[iD][iD2];
   for (iD = 0; iD < nReps; iD++)
   {
      ddata = 0.25 * PSUADE_drand();
      for (iD2 = 0; iD2 < nInputs_; iD2++)
      {
         ddata2 = PSUADE_drand();
         if (ddata2 > 0.5) ddata2 = 1.0;
         else              ddata2 = -1.0;
         if ((M2Mat[iD][iD2]+ddata2*ddata)>1.0 || 
             (M2Mat[iD][iD2]+ddata2*ddata)< 0.0)
              M1Mat[iD][iD2] = M2Mat[iD][iD2] - ddata2*ddata;
         else M1Mat[iD][iD2] = M2Mat[iD][iD2] + ddata2*ddata;
      }
   }
#endif
#if 1
   int    iOne=1;
   double *lbounds = new double[2*nInputs_];
   double *ubounds = new double[2*nInputs_];
   for (iD = 0; iD < 2*nInputs_; iD++) lbounds[iD] = 0.0;
   for (iD = 0; iD < 2*nInputs_; iD++) ubounds[iD] = 1.0;
   Sampling *sampler = (Sampling *) new MCSampling();
   sampler->setInputBounds(2*nInputs_, lbounds, ubounds);
   sampler->setOutputParams(iOne);
   sampler->setSamplingParams(nReps, iOne, iOne);
   sampler->initialize(0);
   double *sInputs  = new double[nReps*2*nInputs_];
   double *sOutputs = new double[nReps];
   int    *sStates  = new int[nReps];
   sampler->getSamples(nReps,2*nInputs_,iOne,sInputs, sOutputs, 
                       sStates);
   for (iD = 0; iD < nReps; iD++)
      for (iD2 = 0; iD2 < nInputs_; iD2++) 
         M1Mat[iD][iD2] = sInputs[iD*2*nInputs_+iD2];
   for (iD = 0; iD < nReps; iD++)
      for (iD2 = 0; iD2 < nInputs_; iD2++) 
         M2Mat[iD][iD2] = sInputs[iD*2*nInputs_+nInputs_+iD2];
   delete [] sInputs;
   delete [] sOutputs;
   delete [] sStates;
   delete [] lbounds;
   delete [] ubounds;
   delete sampler;
#endif
   for (iR = 0; iR < nReps; iR++)
   {
      for (iD2 = 0; iD2 < nInputs_; iD2++)
      {
         ddata = M2Mat[iR][iD2];
         ddata = ddata * ranges[iD2] + lowerBounds_[iD2];
         sampleMatrix_[sampleCount][iD2] = ddata;
      }
      sampleCount++;
      for (iD = 0; iD < nInputs_; iD++)
      {
         for (iD2 = 0; iD2 < nInputs_; iD2++)
         {
            ddata = M2Mat[iR][iD2];
            ddata = ddata * ranges[iD2] + lowerBounds_[iD2];
            sampleMatrix_[sampleCount][iD2] = ddata;
         }
         ddata = M1Mat[iR][iD];
         ddata = ddata * ranges[iD] + lowerBounds_[iD];
         sampleMatrix_[sampleCount][iD] = ddata;
         sampleCount++;
      }
      for (iD2 = 0; iD2 < nInputs_; iD2++)
      {
         ddata = M1Mat[iR][iD2];
         ddata = ddata * ranges[iD2] + lowerBounds_[iD2];
         sampleMatrix_[sampleCount][iD2] = ddata;
      }
      sampleCount++;
   }

   delete [] ranges;
   for (iD = 0;  iD < nReps; iD++)
   {
      delete [] M1Mat[iD];
      delete [] M2Mat[iD];
   }
   delete [] M1Mat;
   delete [] M2Mat;
   return 0;
}

// ************************************************************************
// generate the BS matrix
// ------------------------------------------------------------------------
int SobolSampling::generate(double **inMat, int size)
{
   int iD, iD2, nmax;
#ifdef PSUADE_SAL_GRID
   int idata;
#endif

   nmax = size / 2;
   for (iD = 0; iD < size; iD++)
   {
      for (iD2 = 0; iD2 < nInputs_; iD2++)
      {
#ifdef PSUADE_SAL_GRID
         idata = (int) (PSUADE_drand() * nmax);
         if (idata == nmax) idata = nmax - 1;
         inMat[iD][iD2] = (double) idata / (double) (nmax - 1);
#else
         inMat[iD][iD2] = PSUADE_drand();
#endif
      }
   }
   return 0;
}

// ************************************************************************
// refine the sample space
// ------------------------------------------------------------------------
int SobolSampling::refine(int refineRatio, int randomize, double thresh,
                          int nSamples, double *sampleErrors)
{
   int    iD, iD2, sampleCount, *newStates, nReps, iR, nLevels;
   double **M1Mat, **M2Mat, **newSamples, *newOutputs, *ranges, ddata;

   (void) randomize;
   (void) thresh;
   (void) nSamples;
   (void) sampleErrors;

   nLevels = refineRatio;
   nReps = nSamples_ * (nLevels - 1) / (nInputs_ + 2);

   // First do some defensive programming and range checking by Bill Oliver
   if(nSamples_*nLevels <= 0)
   {
      printf("nSamples_*nLevels <= 0 in file %s line %d\n",__FILE__,__LINE__);
      exit(1);
   }
   newSamples = new double*[nSamples_*nLevels];
   for (iD = 0;  iD < nSamples_*nLevels; iD++)
      newSamples[iD] = new double[nInputs_];
   newOutputs = new double[nSamples_*nLevels*nOutputs_];
   newStates = new int[nSamples_*nLevels];
   for (iD = 0;  iD < nSamples_*nLevels; iD++)
   {
      newStates[iD] = 0;
      for (iD2 = 0; iD2 < nOutputs_; iD2++)
         newOutputs[iD*nOutputs_+iD2] = PSUADE_UNDEFINED;
   }

   for (iD = 0;  iD < nSamples_; iD++) 
   {
      for (iD2 = 0; iD2 < nInputs_; iD2++)
         newSamples[iD][iD2] = sampleMatrix_[iD][iD2];
      for (iD2 = 0; iD2 < nOutputs_; iD2++)
         newOutputs[iD*nOutputs_+iD2] = sampleOutput_[iD*nOutputs_+iD2];
      newStates[iD] = 1;
   }

   M1Mat = new double*[nReps];
   for (iD = 0; iD < nReps; iD++) M1Mat[iD] = new double[nInputs_];
   M2Mat = new double*[nReps];
   for (iD = 0; iD < nReps; iD++) M2Mat[iD] = new double[nInputs_];

   ranges = new double[nInputs_];
   for (iD = 0;  iD < nInputs_;  iD++) 
      ranges[iD] = upperBounds_[iD] - lowerBounds_[iD];

   sampleCount = nSamples_;
   generate(M2Mat, nReps);
#ifdef PSUADE_SAL_GRID
   ddata  = 1.0 / (double) (nReps/2 - 1);
   for (iD = 0; iD < nReps; iD++)
   {
      for (iD2 = 0; iD2 < nInputs_; iD2++)
      {
         ddata2 = PSUADE_drand();
         if (ddata2 > 0.5) ddata2 = 1.0;
         else              ddata2 = -1.0;
         if ((M2Mat[iD][iD2]+ddata2*ddata)>1.0 || 
             (M2Mat[iD][iD2]+ddata2*ddata)< 0.0)
              M1Mat[iD][iD2] = M2Mat[iD][iD2] - ddata2*ddata;
         else M1Mat[iD][iD2] = M2Mat[iD][iD2] + ddata2*ddata;
      }
   }
#endif
#ifdef PSUADE_SAL_GRID
   for (iD = 0; iD < nReps; iD++)
      for (iD2 = 0; iD2 < nInputs_; iD2++)
         M1Mat[iD][iD2] = M2Mat[iD][iD2];
   for (iD = 0; iD < nReps; iD++)
   {
      ddata = 0.25 * PSUADE_drand();
      for (iD2 = 0; iD2 < nInputs_; iD2++)
      {
         ddata2 = PSUADE_drand();
         if (ddata2 > 0.5) ddata2 = 1.0;
         else              ddata2 = -1.0;
         if ((M2Mat[iD][iD2]+ddata2*ddata)>1.0 || 
             (M2Mat[iD][iD2]+ddata2*ddata)< 0.0)
              M1Mat[iD][iD2] = M2Mat[iD][iD2] - ddata2*ddata;
         else M1Mat[iD][iD2] = M2Mat[iD][iD2] + ddata2*ddata;
      }
   }
#endif
#if 1
   int iOne=1;
   double *lbounds  = new double[2*nInputs_];
   double *ubounds  = new double[2*nInputs_];
   for (iD = 0; iD < 2*nInputs_; iD++) lbounds[iD] = 0.0;
   for (iD = 0; iD < 2*nInputs_; iD++) ubounds[iD] = 1.0;
   Sampling *sampler = (Sampling *) new MCSampling();
   sampler->setInputBounds(2*nInputs_, lbounds, ubounds);
   sampler->setOutputParams(iOne);
   sampler->setSamplingParams(nReps*nLevels, iOne, iOne);
   sampler->initialize(0);
   double *sInputs  = new double[nReps*nLevels*nInputs_];
   double *sOutputs = new double[nReps*nLevels];
   int    *sStates  = new int[nReps*nLevels];
   sampler->getSamples(nReps*nLevels, nInputs_, iOne, sInputs, 
                       sOutputs, sStates);
   for (iD = 0; iD < nReps*(nLevels-1); iD++)
      for (iD2 = 0; iD2 < nInputs_; iD2++) 
         M1Mat[iD][iD2] = sInputs[(nReps+iD)*2*nInputs_+iD2];
   for (iD = 0; iD < nReps*(nLevels-1); iD++)
      for (iD2 = 0; iD2 < nInputs_; iD2++) 
         M2Mat[iD][iD2] = sInputs[(nReps+iD)*2*nInputs_+nInputs_+iD2];
   delete [] sInputs;
   delete [] sOutputs;
   delete [] sStates;
   delete [] lbounds;
   delete [] ubounds;
   delete sampler;
#endif
   for (iR = 0; iR < nReps; iR++)
   {
      for (iD2 = 0; iD2 < nInputs_; iD2++)
      {
         ddata = M2Mat[iR][iD2];
         ddata = ddata * ranges[iD2] + lowerBounds_[iD2];
         newSamples[sampleCount][iD2] = ddata;
      }
      sampleCount++;
      for (iD = 0; iD < nInputs_; iD++)
      {
         for (iD2 = 0; iD2 < nInputs_; iD2++)
         {
            ddata = M2Mat[iR][iD2];
            ddata = ddata * ranges[iD2] + lowerBounds_[iD2];
            newSamples[sampleCount][iD2] = ddata;
         }
         ddata = M1Mat[iR][iD];
         ddata = ddata * ranges[iD] + lowerBounds_[iD];
         newSamples[sampleCount][iD] = ddata;
         sampleCount++;
      }
      for (iD2 = 0; iD2 < nInputs_; iD2++)
      {
         ddata = M1Mat[iR][iD2];
         ddata = ddata * ranges[iD2] + lowerBounds_[iD2];
         newSamples[sampleCount][iD2] = ddata;
      }
      sampleCount++;
   }

   deleteSampleData();
   delete [] ranges;
   for (iD = 0;  iD < nReps; iD++)
   {
      delete [] M1Mat[iD];
      delete [] M2Mat[iD];
   }
   delete [] M1Mat;
   delete [] M2Mat;

   nSamples_ = nSamples_ * nLevels;
   sampleMatrix_ = newSamples;
   sampleOutput_ = newOutputs;
   sampleStates_ = newStates;

   if (printLevel_ > 4)
   {
      printf("SobolSampling::refine: nSamples = %d\n", nSamples_);
      printf("SobolSampling::refine: nInputs  = %d\n", nInputs_);
      printf("SobolSampling::refine: nOutputs = %d\n", nOutputs_);
      for (iD2 = 0; iD2 < nInputs_; iD2++)
         printf("    SobolSampling input %3d = [%e %e]\n", iD2+1,
                lowerBounds_[iD2], upperBounds_[iD2]);
   }
   return 0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
SobolSampling& SobolSampling::operator=(const SobolSampling &)
{
   printf("SobolSampling operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

