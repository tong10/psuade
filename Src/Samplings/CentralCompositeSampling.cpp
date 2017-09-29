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
// Functions for the Central Composite sampling class 
// AUTHOR : CHARLES TONG
// DATE   : 2005
// ************************************************************************
#include <stdio.h>
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "CentralCompositeSampling.h"

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
CentralCompositeSampling::CentralCompositeSampling() : Sampling()
{
   samplingID_ = PSUADE_SAMP_CCI4;
   scheme_ = 0;
   resolution_ = 4;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
CentralCompositeSampling::~CentralCompositeSampling()
{
}

// ************************************************************************
// initialize the sampler data
// ------------------------------------------------------------------------
int CentralCompositeSampling::initialize(int initLevel)
{
   int      nSamples, inputID, inputID2, sampleID, *sampleStates;
   double   *sampleInputs, *sampleOutputs, alpha;
   Sampling *samplePtr;

   if (nSamples_ == 0)
   {
      printf("CentralCompositeSampling::initialize ERROR: nSamples = 0.\n");
      exit(1);
   }
   if (nInputs_ == 0 || lowerBounds_ == NULL || upperBounds_ == NULL)
   {
      printf("CentralCompositeSampling::initialize ERROR: no inputs.\n");
      exit(1);
   }

   deleteSampleData();
   if (initLevel == 1) return 0;

   if (resolution_ == 4)
        samplePtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_FF4);
   else if (resolution_ == 5)
        samplePtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_FF5);
   else
   {
      samplePtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_FACT);
      nSamples = 1;
      for (inputID = 0; inputID < nInputs_; inputID++) nSamples *= 2;
   }
   samplePtr->setInputBounds(nInputs_, lowerBounds_, upperBounds_);
   samplePtr->setOutputParams(nOutputs_);
   nSamples = samplePtr->getNumSamples();
   samplePtr->setSamplingParams(nSamples, 1, 0);
   samplePtr->initialize(0);
   // moved this statement above the statement that uses nSamples
   // nSamples = samplePtr->getNumSamples();
   sampleInputs  = new double[nSamples * nInputs_];
   sampleOutputs = new double[nSamples * nOutputs_];
   sampleStates  = new int[nSamples];
   samplePtr->getSamples(nSamples, nInputs_, nOutputs_, sampleInputs,
                         sampleOutputs, sampleStates);
   delete [] sampleOutputs;
   delete [] sampleStates;

   nSamples_ = nSamples + 2 * nInputs_ + 1;
   allocSampleData();
   if (scheme_ == 0) alpha = pow((double) nSamples, 0.25);
   else              alpha = 1.0;

   for (sampleID = 0; sampleID < nSamples; sampleID++)
   {
      for (inputID = 0; inputID < nInputs_; inputID++)
      {
         if (alpha == 1.0)
            sampleMatrix_[sampleID][inputID] = 
               sampleInputs[sampleID*nInputs_+inputID];
         else
         {
            if (sampleInputs[sampleID*nInputs_+inputID]==lowerBounds_[inputID])
               sampleMatrix_[sampleID][inputID] = 0.5 * 
                  (lowerBounds_[inputID] + upperBounds_[inputID]) + 
                  0.5/alpha*(lowerBounds_[inputID]-upperBounds_[inputID]); 
            else
               sampleMatrix_[sampleID][inputID] = 0.5 *  
                  (lowerBounds_[inputID] + upperBounds_[inputID]) + 
                  0.5/alpha*(upperBounds_[inputID]-lowerBounds_[inputID]); 
         }
      }
   }
   delete [] sampleInputs;

   if (printLevel_ > 4)
   {
      printf("CentralCompositeSampling::initialize: nSamples = %d\n",nSamples_);
      printf("CentralCompositeSampling::initialize: nInputs  = %d\n",nInputs_);
      printf("CentralCompositeSampling::initialize: nOutputs = %d\n",nOutputs_);
      for (inputID = 0; inputID < nInputs_; inputID++)
         printf("   CentralCompositeSampling input %3d = [%e %e]\n",
                inputID+1, lowerBounds_[inputID], upperBounds_[inputID]);
   }

   if (scheme_ == 2) alpha = pow((double) nSamples, 0.25);
   else              alpha = 1.0;

   for (inputID = 0; inputID < nInputs_; inputID++)
   {
      for (inputID2 = 0; inputID2 < nInputs_; inputID2++)
      {
         if (inputID == inputID2)
         {
            sampleMatrix_[nSamples+2*inputID][inputID2] = 0.5 *  
                  (lowerBounds_[inputID] + upperBounds_[inputID]) + 
                  0.5*alpha*(lowerBounds_[inputID]-upperBounds_[inputID]); 
            sampleMatrix_[nSamples+2*inputID+1][inputID2] = 0.5 * 
                  (lowerBounds_[inputID] + upperBounds_[inputID]) + 
                  0.5*alpha*(upperBounds_[inputID]-lowerBounds_[inputID]); 
         }
         else
         {
            sampleMatrix_[nSamples+2*inputID][inputID2] = 0.5 * 
                             (lowerBounds_[inputID2] + upperBounds_[inputID2]); 
            sampleMatrix_[nSamples+2*inputID+1][inputID2] = 0.5 * 
                             (lowerBounds_[inputID2] + upperBounds_[inputID2]); 
         }
      }
   }
   for (inputID = 0; inputID < nInputs_; inputID++)
       sampleMatrix_[nSamples_-1][inputID] = 0.5 * 
                    (lowerBounds_[inputID] + upperBounds_[inputID]);
   // Cleanup by Bill Oliver
   delete samplePtr;
   return 0;
}

// ************************************************************************
// refine the sample space
// ------------------------------------------------------------------------
int CentralCompositeSampling::refine(int,int,double, int, double *)
{
   printf("CentralCompositeSampling::refine ERROR - not available.\n");
   exit(1);
   return 0;
}

// ************************************************************************
// set internal scheme
// ------------------------------------------------------------------------
int CentralCompositeSampling::setParam(char *sparam)
{
   char winput[501];
   sscanf(sparam, "%s", winput);
   if (!strcmp(winput, "setResolution"))
   {
      sscanf(sparam, "%s %d", winput, &resolution_);
      if (resolution_ != 4 && resolution_ != 5 && resolution_ != 0)
         resolution_ = 4;
   }
   else if (!strcmp(winput, "setScheme"))
   {
      sscanf(sparam, "%s %d", winput, &scheme_);
      if (scheme_ < 0 && scheme_ > 2) scheme_ = 0;
   }
   if (resolution_ == 0 && scheme_ == 0) samplingID_ = PSUADE_SAMP_CCIF;
   if (resolution_ == 4 && scheme_ == 0) samplingID_ = PSUADE_SAMP_CCI4;
   if (resolution_ == 5 && scheme_ == 0) samplingID_ = PSUADE_SAMP_CCI5;
   if (resolution_ == 0 && scheme_ == 1) samplingID_ = PSUADE_SAMP_CCFF;
   if (resolution_ == 4 && scheme_ == 1) samplingID_ = PSUADE_SAMP_CCF4;
   if (resolution_ == 5 && scheme_ == 1) samplingID_ = PSUADE_SAMP_CCF5;
   if (resolution_ == 0 && scheme_ == 2) samplingID_ = PSUADE_SAMP_CCCF;
   if (resolution_ == 4 && scheme_ == 2) samplingID_ = PSUADE_SAMP_CCC4;
   if (resolution_ == 5 && scheme_ == 2) samplingID_ = PSUADE_SAMP_CCC5;
   return 0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
CentralCompositeSampling& CentralCompositeSampling::operator=
                                  (const CentralCompositeSampling &)
{
   printf("CentralCompositeSampling operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

