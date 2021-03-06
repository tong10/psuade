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
// Functions for the Fourier Amplitude Sampling Test (FAST) class 
// AUTHOR : CHARLES TONG
// DATE   : 2005 (verified against SobolG problem - Jan 2006)
// ************************************************************************
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "sysdef.h"
#include "PsuadeUtil.h"
#include "FASTSampling.h"

#define PABS(x) ((x) > 0 ? x : -(x))

// ------------------------------------------------------------------------
// local definitions
// ------------------------------------------------------------------------
#define  PSUADE_FAST_MaxDimension  50

static unsigned long
PSUADE_FAST_OMEGA[PSUADE_FAST_MaxDimension] =
{
      1,   3,     1,    5,   11,    1,   17,   23,   19,   25,   
     41,   31,   23,   87,   67,   73,   85,  143,  149,   99, 
    119,  237,  267,  283,  151,  385,  157,  215,  449,  163,
    337,  253,  375,  441,  673,  773,  875,  873,  587,  849, 
    623,  637,  891,  943, 1171, 1225, 1335, 1725, 1663, 2019
};

static unsigned long
PSUADE_FAST_DELTA[PSUADE_FAST_MaxDimension] =
{
      4,    8,    6,   10,   20,   22,   32,   40,   38,   26,
     56,   62,   46,   76,   96,   60,   86,  126,  134,  112,
     92,  128,  154,  196,   34,  416,  106,  208,  328,  198,
    382,   88,  348,  186,  140,  170,  284,  568,  302,  438,
    410,  248,  448,  388,  596,  216,  100,  488,  166,    0
};

// ************************************************************************
// Constructor 
// ------------------------------------------------------------------------
FASTSampling::FASTSampling() : Sampling()
{
}

// *******************************************************************
// destructor 
// -------------------------------------------------------------------
FASTSampling::~FASTSampling()
{
}

// ************************************************************************
// initialize the sampling data (default M = 4)
// ------------------------------------------------------------------------
int FASTSampling::initialize(int initLevel)
{
   int    M=4, *omegas, inputID, ii;
   double *ranges, ds, ss, ps_pi=3.14159, ddata;

   if (nSamples_ == 0)
   {
      printf("FASTSampling::initialize ERROR - nSamples = 0.\n");
      exit(1);
   }
   if (nInputs_ == 0 || lowerBounds_ == NULL || upperBounds_ == NULL)
   {
      printf("FASTSampling::initialize ERROR - input not set up.\n");
      exit(1);
   }

   deleteSampleData();
   if (nInputs_ > 50) 
   {
      printf("FASTSampling : cannot handle nInputs > 50.\n");
      exit(1);
   }
   if (initLevel != 0) return 0;

   omegas = new int[nInputs_];
   calculateOmegas(nInputs_, omegas);

   ii = (nSamples_ - 1) / omegas[nInputs_-1] / 2;
   if ((2*ii*omegas[nInputs_-1]+1) == nSamples_ && ii >= 2) M = ii;
   else nSamples_ = 2 * M * omegas[nInputs_-1] + 1;

   if (printLevel_ > 4)
   {
      printf("FASTSampling::initialize: nSamples = %d\n", nSamples_);
      printf("FASTSampling::initialize: nInputs  = %d\n", nInputs_);
      printf("FASTSampling::initialize: nOutputs = %d\n", nOutputs_);
      for (inputID = 0; inputID < nInputs_; inputID++)
         printf("    FASTSampling input %3d = [%e %e]\n", inputID+1,
                lowerBounds_[inputID], upperBounds_[inputID]);
      for (inputID = 0; inputID < nInputs_; inputID++)
         printf("    FASTSampling::input %4d fundamental frequency = %d\n",
                inputID+1, omegas[inputID]);
   }

   allocSampleData();
   ds = ps_pi / (double) (2 * nSamples_);
   for (ii = 0; ii < nSamples_; ii++)
   {
      ss = - 0.5 * ps_pi + ds * (2 * ii + 1);
      for (inputID = 0; inputID < nInputs_; inputID++)
      {
         ddata = ss * omegas[inputID];
         sampleMatrix_[ii][inputID] = 0.5 + 1.0 / ps_pi * asin(sin(ddata));
      }
   }

   ranges = new double[nInputs_];
   for (inputID = 0; inputID < nInputs_; inputID++) 
      ranges[inputID] = upperBounds_[inputID] - lowerBounds_[inputID];

   for (ii = 0; ii < nSamples_; ii++) 
   {
      for (inputID = 0; inputID < nInputs_; inputID++)
      {
         ddata = sampleMatrix_[ii][inputID];
         sampleMatrix_[ii][inputID] = ddata * ranges[inputID] + 
                                      lowerBounds_[inputID];
      }
   } 

   delete [] ranges;
   delete [] omegas;
   return 0;
}

// ************************************************************************
// refine the sample space
// ------------------------------------------------------------------------
int FASTSampling::refine(int refineRatio, int randomize, double thresh,
                         int nSamples, double *sampleErrors)
{
   int    *omegas, M, oldNSamples, *oldStates, ii, inputID, fail;
   double *ranges, **oldSamples, *oldOutputs, ds, ss, ps_pi=3.14159, ddata;

   (void) randomize;
   (void) thresh;
   (void) nSamples;
   (void) sampleErrors;

   omegas = new int[nInputs_];
   calculateOmegas(nInputs_, omegas);
   M = (nSamples_ - 1) / omegas[nInputs_-1];
   if ((M*omegas[nInputs_-1]+1) != nSamples_)
   {
      printf("FASTSampling refine ERROR : invalid sample size.\n");
      delete [] omegas;
      exit(1);
   } 
   
   oldNSamples = nSamples_;
   oldSamples  = sampleMatrix_;
   oldOutputs  = sampleOutput_;
   oldStates   = sampleStates_;
   sampleMatrix_ = NULL;
   sampleOutput_ = NULL;
   sampleStates_ = NULL;
   deleteSampleData();
   nSamples_ = (nSamples_ - 1) * 2 + 1;
   allocSampleData();

   if (printLevel_ > 4)
   {
      printf("FASTSampling::refine: nSamples = %d\n", nSamples_);
      printf("FASTSampling::refine: nInputs  = %d\n", nInputs_);
      printf("FASTSampling::refine: nOutputs = %d\n", nOutputs_);
      for (inputID = 0; inputID < nInputs_; inputID++)
         printf("    FASTSampling input %3d = [%e %e]\n", inputID+1,
                lowerBounds_[inputID], upperBounds_[inputID]);
   }

   for (ii = 0; ii < nSamples_; ii++)
   {
      ds = ps_pi / (double) (2 * nSamples_);
      ss = - 0.5 * ps_pi + ds * (2 * ii + 1);
      for (inputID = 0; inputID < nInputs_; inputID++)
      {
         ddata = ss * omegas[inputID];
         sampleMatrix_[ii][inputID] = 0.5 + 1.0 / ps_pi * asin(sin(ddata));
      }
   }

   ranges = new double[nInputs_];
   for (inputID = 0; inputID < nInputs_; inputID++) 
      ranges[inputID] = upperBounds_[inputID] - lowerBounds_[inputID];

   for (ii = 0; ii < nSamples_; ii++) 
   {
      for (inputID = 0; inputID < nInputs_; inputID++)
      {
         ddata = sampleMatrix_[ii][inputID];
         sampleMatrix_[ii][inputID] = ddata * ranges[inputID] + 
                                      lowerBounds_[inputID];
      }
   } 

   fail = 0;
   for (ii = 0; ii < oldNSamples; ii++)
   {
      for (inputID = 0; inputID < nInputs_; inputID++)
      {
         ddata = oldSamples[ii][inputID] - sampleMatrix_[ii*2][inputID];
         if (PABS(ddata) > 1.0e-13) fail = 1;
      }
      for (inputID = 0; inputID < nOutputs_; inputID++)
         sampleOutput_[ii*2*nOutputs_+inputID] = 
                                     oldOutputs[ii*nOutputs_+inputID];
      sampleStates_[ii*2] = oldStates[ii];
   }
   if (fail == 1) printf("FASTSampling fails checking.\n");

   delete [] omegas;
   delete [] ranges;
   return 0;
}

// ************************************************************************
// calculate frequencies
// ------------------------------------------------------------------------
int FASTSampling::calculateOmegas(int nInputs, int *omegas)
{
   int ii;

   omegas[0] = PSUADE_FAST_OMEGA[nInputs-1];
   for (ii = 1; ii < nInputs; ii++)
      omegas[ii] = omegas[ii-1] + PSUADE_FAST_DELTA[nInputs-1-ii];
   return 0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
FASTSampling& FASTSampling::operator=(const FASTSampling &)
{
   printf("FASTSampling operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

