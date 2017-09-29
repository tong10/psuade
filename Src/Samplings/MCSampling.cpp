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
// Functions for the Monte Carlo sampling class
// AUTHOR : CHARLES TONG
// DATE   : 2003
// ************************************************************************
using namespace std;
#include "Util/sysdef.h"
#include "Util/PsuadeUtil.h"
#include "MCSampling.h"

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
MCSampling::MCSampling() : Sampling()
{
   samplingID_ = PSUADE_SAMP_MC;
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
MCSampling::~MCSampling()
{
}

// ************************************************************************
// initialization 
// ------------------------------------------------------------------------
int MCSampling::initialize(int initLevel)
{
   int    ii, inputID;
   double range;

   if (nInputs_ == 0 || lowerBounds_ == NULL || upperBounds_ == NULL)
   {
      printf("MCSampling::initialize ERROR - input not set up.\n");
      exit(1);
   }

   if (printLevel_ > 4)
   {
      printf("MCSampling::initialize: nSamples = %d\n", nSamples_);
      printf("MCSampling::initialize: nInputs  = %d\n", nInputs_);
      printf("MCSampling::initialize: nOutputs = %d\n", nOutputs_);
      for (inputID = 0; inputID < nInputs_; inputID++)
         printf("    MCSampling input %3d = [%e %e]\n", inputID+1,
                lowerBounds_[inputID], upperBounds_[inputID]);
   }

   deleteSampleData();
   if (initLevel != 0) return 0;

   allocSampleData();
   for (inputID = 0; inputID < nInputs_; inputID++)
   { 
      range = upperBounds_[inputID] - lowerBounds_[inputID];
      for (ii = 0; ii < nSamples_; ii++) 
         sampleMatrix_[ii][inputID] = PSUADE_drand() * range + 
                                      lowerBounds_[inputID];
   }
   return 0;
}

// ************************************************************************
// refine 
// ------------------------------------------------------------------------
int MCSampling::refine(int refineRatio, int randomize, double thresh,
                       int nSamples, double *sampleErrors)
{
   int    ii, inputID, oldNumSamples, *oldSampleStates;
   int    nLevels;
   double **oldSampleMatrix, *oldSampleOutput, range;

   (void) randomize;
   (void) thresh;
   (void) nSamples;
   (void) sampleErrors;

   nLevels         = refineRatio;
   oldNumSamples   = nSamples_;
   oldSampleMatrix = sampleMatrix_;
   oldSampleOutput = sampleOutput_;
   oldSampleStates = sampleStates_;

   nSamples_ *= nLevels;
   allocSampleData();

   if (printLevel_ > 4)
   {
      printf("MCSampling::refine: nSamples = %d\n", nSamples_);
      printf("MCSampling::refine: nInputs  = %d\n", nInputs_);
      printf("MCSampling::refine: nOutputs = %d\n", nOutputs_);
      for (inputID = 0; inputID < nInputs_; inputID++)
         printf("    MCSampling input %3d = [%e %e]\n", inputID+1,
                lowerBounds_[inputID], upperBounds_[inputID]);
   }

   for (ii = 0; ii < oldNumSamples; ii++) 
   {
      for (inputID = 0; inputID < nInputs_; inputID++) 
         sampleMatrix_[ii][inputID] = oldSampleMatrix[ii][inputID];
      for (inputID = 0; inputID < nOutputs_; inputID++) 
         sampleOutput_[ii*nOutputs_+inputID] = 
                    oldSampleOutput[ii*nOutputs_+inputID];
      sampleStates_[ii] = oldSampleStates[ii];
   }
   for (inputID = 0; inputID < nInputs_; inputID++) 
   {
      range = upperBounds_[inputID] - lowerBounds_[inputID];
      for (ii = oldNumSamples; ii < nSamples_; ii++)
         sampleMatrix_[ii][inputID] = PSUADE_drand() * range + 
                                      lowerBounds_[inputID];
   }
   for (ii = 0; ii < oldNumSamples; ii++) delete [] oldSampleMatrix[ii];
   delete [] oldSampleMatrix;
   delete [] oldSampleOutput;
   delete [] oldSampleStates;
   return 0;
}

