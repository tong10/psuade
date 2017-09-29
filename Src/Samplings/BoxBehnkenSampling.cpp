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
// Functions for the Box-Behnken sampling class 
// This sample design is compatible with quadratic regression
// AUTHOR : CHARLES TONG
// DATE   : 2005
// ************************************************************************
#include <stdio.h>
using namespace std;

#include "Util/sysdef.h"
#include "Util/PsuadeUtil.h"
#include "BoxBehnkenSampling.h"

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
BoxBehnkenSampling::BoxBehnkenSampling() : Sampling()
{
   samplingID_ = PSUADE_SAMP_BBD;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
BoxBehnkenSampling::~BoxBehnkenSampling()
{
}

// ************************************************************************
// initialize the sampler data
// ------------------------------------------------------------------------
int BoxBehnkenSampling::initialize(int initLevel)
{
   int inputID, sampleID, inputID2;

   if (nInputs_ == 0 || lowerBounds_ == NULL || upperBounds_ == NULL)
   {
      printf("BoxBehnkenSampling::initialize ERROR - input not set up.\n");
      exit(1);
   }
   if (nInputs_ < 2)
   {
      printf("BoxBehnkenSampling::initialize ERROR - nInputs < 2.\n");
      printf("                    nInputs should be [2,7].\n");
      exit(1);
   }
   if (nInputs_ > 7)
   {
      printf("BoxBehnkenSampling::initialize ERROR - nInputs > 7.\n");
      printf("                    nInputs should be [2,7].\n");
      exit(1);
   }
   
   deleteSampleData();
   if (initLevel == 1) return 0;
   if      (nInputs_ == 2) nSamples_ = 5;
   else if (nInputs_ == 3) nSamples_ = 13;
   else if (nInputs_ == 4) nSamples_ = 25;
   else if (nInputs_ == 5) nSamples_ = 41;
   else if (nInputs_ == 6) nSamples_ = 48;
   else if (nInputs_ == 7) nSamples_ = 57;
   allocSampleData();

   if (printLevel_ > 4)
   {
      printf("BoxBehnkenSampling::initialize: nSamples = %d\n", nSamples_);
      printf("BoxBehnkenSampling::initialize: nInputs  = %d\n", nInputs_);
      printf("BoxBehnkenSampling::initialize: nOutputs = %d\n", nOutputs_);
      for (inputID = 0; inputID < nInputs_; inputID++)
         printf("    BoxBehnkenSampling input %3d = [%e %e]\n", inputID+1,
                lowerBounds_[inputID], upperBounds_[inputID]);
   }

   for (sampleID = 0; sampleID < nSamples_; sampleID++)
      for (inputID = 0; inputID < nInputs_; inputID++)
         sampleMatrix_[sampleID][inputID] = 
                0.5 * (lowerBounds_[inputID] + upperBounds_[inputID]);
   if (nInputs_ >= 2 && nInputs_ <= 5)
   {
      sampleID = 0;
      for (inputID = 0; inputID < nInputs_; inputID++)
      {
         for (inputID2 = inputID+1; inputID2 < nInputs_; inputID2++)
         {
            sampleMatrix_[sampleID][inputID]  = lowerBounds_[inputID];
            sampleMatrix_[sampleID][inputID2] = lowerBounds_[inputID2];
            sampleMatrix_[sampleID+1][inputID]  = lowerBounds_[inputID];
            sampleMatrix_[sampleID+1][inputID2] = upperBounds_[inputID2];
            sampleMatrix_[sampleID+2][inputID]  = upperBounds_[inputID];
            sampleMatrix_[sampleID+2][inputID2] = lowerBounds_[inputID2];
            sampleMatrix_[sampleID+3][inputID]  = upperBounds_[inputID];
            sampleMatrix_[sampleID+3][inputID2] = upperBounds_[inputID2];
            sampleID += 4;
         }
      }
   }
   if (nInputs_ == 6)
   {
      sampleID = 0;
      for (inputID = 0; inputID < nInputs_; inputID++)
      {
         sampleMatrix_[sampleID][inputID] = lowerBounds_[inputID];
         inputID2 = (inputID + 1) % nInputs_;
         sampleMatrix_[sampleID][inputID2]  = lowerBounds_[inputID2];
         inputID2 = (inputID + 3) % nInputs_;
         sampleMatrix_[sampleID][inputID2] = lowerBounds_[inputID2];
         sampleMatrix_[sampleID+1][inputID] = lowerBounds_[inputID];
         inputID2 = (inputID + 1) % nInputs_;
         sampleMatrix_[sampleID+1][inputID2]  = lowerBounds_[inputID2];
         inputID2 = (inputID + 3) % nInputs_;
         sampleMatrix_[sampleID+1][inputID2] = upperBounds_[inputID2];
         sampleMatrix_[sampleID+2][inputID] = lowerBounds_[inputID];
         inputID2 = (inputID + 1) % nInputs_;
         sampleMatrix_[sampleID+2][inputID2]  = upperBounds_[inputID2];
         inputID2 = (inputID + 3) % nInputs_;
         sampleMatrix_[sampleID+2][inputID2] = lowerBounds_[inputID2];
         sampleMatrix_[sampleID+3][inputID] = lowerBounds_[inputID];
         inputID2 = (inputID + 1) % nInputs_;
         sampleMatrix_[sampleID+3][inputID2]  = upperBounds_[inputID2];
         inputID2 = (inputID + 3) % nInputs_;
         sampleMatrix_[sampleID+3][inputID2] = upperBounds_[inputID2];

         sampleMatrix_[sampleID+4][inputID] = upperBounds_[inputID];
         inputID2 = (inputID + 1) % nInputs_;
         sampleMatrix_[sampleID+4][inputID2]  = lowerBounds_[inputID2];
         inputID2 = (inputID + 3) % nInputs_;
         sampleMatrix_[sampleID+4][inputID2] = lowerBounds_[inputID2];
         sampleMatrix_[sampleID+5][inputID] = upperBounds_[inputID];
         inputID2 = (inputID + 1) % nInputs_;
         sampleMatrix_[sampleID+5][inputID2]  = lowerBounds_[inputID2];
         inputID2 = (inputID + 3) % nInputs_;
         sampleMatrix_[sampleID+5][inputID2] = upperBounds_[inputID2];
         sampleMatrix_[sampleID+6][inputID] = upperBounds_[inputID];
         inputID2 = (inputID + 1) % nInputs_;
         sampleMatrix_[sampleID+6][inputID2]  = upperBounds_[inputID2];
         inputID2 = (inputID + 3) % nInputs_;
         sampleMatrix_[sampleID+6][inputID2] = lowerBounds_[inputID2];
         sampleMatrix_[sampleID+7][inputID] = upperBounds_[inputID];
         inputID2 = (inputID + 1) % nInputs_;
         sampleMatrix_[sampleID+7][inputID2]  = upperBounds_[inputID2];
         inputID2 = (inputID + 3) % nInputs_;
         sampleMatrix_[sampleID+7][inputID2] = upperBounds_[inputID2];
         sampleID += 8;
      }
   }
   if (nInputs_ == 7)
   {
printf("check 1\n");
      sampleID = 0;
      sampleMatrix_[sampleID][3] = lowerBounds_[3];
      sampleMatrix_[sampleID][4] = lowerBounds_[4];
      sampleMatrix_[sampleID][5] = lowerBounds_[5];
      sampleMatrix_[sampleID+1][3] = lowerBounds_[3];
      sampleMatrix_[sampleID+1][4] = lowerBounds_[4];
      sampleMatrix_[sampleID+1][5] = upperBounds_[5];
      sampleMatrix_[sampleID+2][3] = lowerBounds_[3];
      sampleMatrix_[sampleID+2][4] = upperBounds_[4];
      sampleMatrix_[sampleID+2][5] = lowerBounds_[5];
      sampleMatrix_[sampleID+3][3] = lowerBounds_[3];
      sampleMatrix_[sampleID+3][4] = upperBounds_[4];
      sampleMatrix_[sampleID+3][5] = upperBounds_[5];
      sampleMatrix_[sampleID+4][3] = upperBounds_[3];
      sampleMatrix_[sampleID+4][4] = lowerBounds_[4];
      sampleMatrix_[sampleID+4][5] = lowerBounds_[5];
      sampleMatrix_[sampleID+5][3] = upperBounds_[3];
      sampleMatrix_[sampleID+5][4] = lowerBounds_[4];
      sampleMatrix_[sampleID+5][5] = upperBounds_[5];
      sampleMatrix_[sampleID+6][3] = upperBounds_[3];
      sampleMatrix_[sampleID+6][4] = upperBounds_[4];
      sampleMatrix_[sampleID+6][5] = lowerBounds_[5];
      sampleMatrix_[sampleID+7][3] = upperBounds_[3];
      sampleMatrix_[sampleID+7][4] = upperBounds_[4];
      sampleMatrix_[sampleID+7][5] = upperBounds_[5];
      sampleID += 8;
      sampleMatrix_[sampleID][0] = lowerBounds_[0];
      sampleMatrix_[sampleID][5] = lowerBounds_[5];
      sampleMatrix_[sampleID][6] = lowerBounds_[6];
      sampleMatrix_[sampleID+1][0] = lowerBounds_[0];
      sampleMatrix_[sampleID+1][5] = lowerBounds_[5];
      sampleMatrix_[sampleID+1][6] = upperBounds_[6];
      sampleMatrix_[sampleID+2][0] = lowerBounds_[0];
      sampleMatrix_[sampleID+2][5] = upperBounds_[5];
      sampleMatrix_[sampleID+2][6] = lowerBounds_[6];
      sampleMatrix_[sampleID+3][0] = lowerBounds_[0];
      sampleMatrix_[sampleID+3][5] = upperBounds_[5];
      sampleMatrix_[sampleID+3][6] = upperBounds_[6];
      sampleMatrix_[sampleID+4][0] = upperBounds_[0];
      sampleMatrix_[sampleID+4][5] = lowerBounds_[5];
      sampleMatrix_[sampleID+4][6] = lowerBounds_[6];
      sampleMatrix_[sampleID+5][0] = upperBounds_[0];
      sampleMatrix_[sampleID+5][5] = lowerBounds_[5];
      sampleMatrix_[sampleID+5][6] = upperBounds_[6];
      sampleMatrix_[sampleID+6][0] = upperBounds_[0];
      sampleMatrix_[sampleID+6][5] = upperBounds_[5];
      sampleMatrix_[sampleID+6][6] = lowerBounds_[6];
      sampleMatrix_[sampleID+7][0] = upperBounds_[0];
      sampleMatrix_[sampleID+7][5] = upperBounds_[5];
      sampleMatrix_[sampleID+7][6] = upperBounds_[6];
      sampleID += 8;
      sampleMatrix_[sampleID][1] = lowerBounds_[1];
      sampleMatrix_[sampleID][4] = lowerBounds_[4];
      sampleMatrix_[sampleID][6] = lowerBounds_[6];
      sampleMatrix_[sampleID+1][1] = lowerBounds_[1];
      sampleMatrix_[sampleID+1][4] = lowerBounds_[4];
      sampleMatrix_[sampleID+1][6] = upperBounds_[6];
      sampleMatrix_[sampleID+2][1] = lowerBounds_[1];
      sampleMatrix_[sampleID+2][4] = upperBounds_[4];
      sampleMatrix_[sampleID+2][6] = lowerBounds_[6];
      sampleMatrix_[sampleID+3][1] = lowerBounds_[1];
      sampleMatrix_[sampleID+3][4] = upperBounds_[4];
      sampleMatrix_[sampleID+3][6] = upperBounds_[6];
      sampleMatrix_[sampleID+4][1] = upperBounds_[1];
      sampleMatrix_[sampleID+4][4] = lowerBounds_[4];
      sampleMatrix_[sampleID+4][6] = lowerBounds_[6];
      sampleMatrix_[sampleID+5][1] = upperBounds_[1];
      sampleMatrix_[sampleID+5][4] = lowerBounds_[4];
      sampleMatrix_[sampleID+5][6] = upperBounds_[6];
      sampleMatrix_[sampleID+6][1] = upperBounds_[1];
      sampleMatrix_[sampleID+6][4] = upperBounds_[4];
      sampleMatrix_[sampleID+6][6] = lowerBounds_[6];
      sampleMatrix_[sampleID+7][1] = upperBounds_[1];
      sampleMatrix_[sampleID+7][4] = upperBounds_[4];
      sampleMatrix_[sampleID+7][6] = upperBounds_[6];
      sampleID += 8;
      sampleMatrix_[sampleID][0] = lowerBounds_[0];
      sampleMatrix_[sampleID][1] = lowerBounds_[1];
      sampleMatrix_[sampleID][3] = lowerBounds_[3];
      sampleMatrix_[sampleID+1][0] = lowerBounds_[0];
      sampleMatrix_[sampleID+1][1] = lowerBounds_[1];
      sampleMatrix_[sampleID+1][3] = upperBounds_[3];
      sampleMatrix_[sampleID+2][0] = lowerBounds_[0];
      sampleMatrix_[sampleID+2][1] = upperBounds_[1];
      sampleMatrix_[sampleID+2][3] = lowerBounds_[3];
      sampleMatrix_[sampleID+3][0] = lowerBounds_[0];
      sampleMatrix_[sampleID+3][1] = upperBounds_[1];
      sampleMatrix_[sampleID+3][3] = upperBounds_[3];
      sampleMatrix_[sampleID+4][0] = upperBounds_[0];
      sampleMatrix_[sampleID+4][1] = lowerBounds_[1];
      sampleMatrix_[sampleID+4][3] = lowerBounds_[3];
      sampleMatrix_[sampleID+5][0] = upperBounds_[0];
      sampleMatrix_[sampleID+5][1] = lowerBounds_[1];
      sampleMatrix_[sampleID+5][3] = upperBounds_[3];
      sampleMatrix_[sampleID+6][0] = upperBounds_[0];
      sampleMatrix_[sampleID+6][1] = upperBounds_[1];
      sampleMatrix_[sampleID+6][3] = lowerBounds_[3];
      sampleMatrix_[sampleID+7][0] = upperBounds_[0];
      sampleMatrix_[sampleID+7][1] = upperBounds_[1];
      sampleMatrix_[sampleID+7][3] = upperBounds_[3];
      sampleID += 8;
      sampleMatrix_[sampleID][2] = lowerBounds_[2];
      sampleMatrix_[sampleID][3] = lowerBounds_[3];
      sampleMatrix_[sampleID][6] = lowerBounds_[6];
      sampleMatrix_[sampleID+1][2] = lowerBounds_[2];
      sampleMatrix_[sampleID+1][3] = lowerBounds_[3];
      sampleMatrix_[sampleID+1][6] = upperBounds_[6];
      sampleMatrix_[sampleID+2][2] = lowerBounds_[2];
      sampleMatrix_[sampleID+2][3] = upperBounds_[3];
      sampleMatrix_[sampleID+2][6] = lowerBounds_[6];
      sampleMatrix_[sampleID+3][2] = lowerBounds_[2];
      sampleMatrix_[sampleID+3][3] = upperBounds_[3];
      sampleMatrix_[sampleID+3][6] = upperBounds_[6];
      sampleMatrix_[sampleID+4][2] = upperBounds_[2];
      sampleMatrix_[sampleID+4][3] = lowerBounds_[3];
      sampleMatrix_[sampleID+4][6] = lowerBounds_[6];
      sampleMatrix_[sampleID+5][2] = upperBounds_[2];
      sampleMatrix_[sampleID+5][3] = lowerBounds_[3];
      sampleMatrix_[sampleID+5][6] = upperBounds_[6];
      sampleMatrix_[sampleID+6][2] = upperBounds_[2];
      sampleMatrix_[sampleID+6][3] = upperBounds_[3];
      sampleMatrix_[sampleID+6][6] = lowerBounds_[6];
      sampleMatrix_[sampleID+7][2] = upperBounds_[2];
      sampleMatrix_[sampleID+7][3] = upperBounds_[3];
      sampleMatrix_[sampleID+7][6] = upperBounds_[6];
      sampleID += 8;
      sampleMatrix_[sampleID][0] = lowerBounds_[0];
      sampleMatrix_[sampleID][2] = lowerBounds_[2];
      sampleMatrix_[sampleID][4] = lowerBounds_[4];
      sampleMatrix_[sampleID+1][0] = lowerBounds_[0];
      sampleMatrix_[sampleID+1][2] = lowerBounds_[2];
      sampleMatrix_[sampleID+1][4] = upperBounds_[4];
      sampleMatrix_[sampleID+2][0] = lowerBounds_[0];
      sampleMatrix_[sampleID+2][2] = upperBounds_[2];
      sampleMatrix_[sampleID+2][4] = lowerBounds_[4];
      sampleMatrix_[sampleID+3][0] = lowerBounds_[0];
      sampleMatrix_[sampleID+3][2] = upperBounds_[2];
      sampleMatrix_[sampleID+3][4] = upperBounds_[4];
      sampleMatrix_[sampleID+4][0] = upperBounds_[0];
      sampleMatrix_[sampleID+4][2] = lowerBounds_[2];
      sampleMatrix_[sampleID+4][4] = lowerBounds_[4];
      sampleMatrix_[sampleID+5][0] = upperBounds_[0];
      sampleMatrix_[sampleID+5][2] = lowerBounds_[2];
      sampleMatrix_[sampleID+5][4] = upperBounds_[4];
      sampleMatrix_[sampleID+6][0] = upperBounds_[0];
      sampleMatrix_[sampleID+6][2] = upperBounds_[2];
      sampleMatrix_[sampleID+6][4] = lowerBounds_[4];
      sampleMatrix_[sampleID+7][0] = upperBounds_[0];
      sampleMatrix_[sampleID+7][2] = upperBounds_[2];
      sampleMatrix_[sampleID+7][4] = upperBounds_[4];
      sampleID += 8;
      sampleMatrix_[sampleID][1] = lowerBounds_[1];
      sampleMatrix_[sampleID][2] = lowerBounds_[2];
      sampleMatrix_[sampleID][5] = lowerBounds_[5];
      sampleMatrix_[sampleID+1][1] = lowerBounds_[1];
      sampleMatrix_[sampleID+1][2] = lowerBounds_[2];
      sampleMatrix_[sampleID+1][5] = upperBounds_[5];
      sampleMatrix_[sampleID+2][1] = lowerBounds_[1];
      sampleMatrix_[sampleID+2][2] = upperBounds_[2];
      sampleMatrix_[sampleID+2][5] = lowerBounds_[5];
      sampleMatrix_[sampleID+3][1] = lowerBounds_[1];
      sampleMatrix_[sampleID+3][2] = upperBounds_[2];
      sampleMatrix_[sampleID+3][5] = upperBounds_[5];
      sampleMatrix_[sampleID+4][1] = upperBounds_[1];
      sampleMatrix_[sampleID+4][2] = lowerBounds_[2];
      sampleMatrix_[sampleID+4][5] = lowerBounds_[5];
      sampleMatrix_[sampleID+5][1] = upperBounds_[1];
      sampleMatrix_[sampleID+5][2] = lowerBounds_[2];
      sampleMatrix_[sampleID+5][5] = upperBounds_[5];
      sampleMatrix_[sampleID+6][1] = upperBounds_[1];
      sampleMatrix_[sampleID+6][2] = upperBounds_[2];
      sampleMatrix_[sampleID+6][5] = lowerBounds_[5];
      sampleMatrix_[sampleID+7][1] = upperBounds_[1];
      sampleMatrix_[sampleID+7][2] = upperBounds_[2];
      sampleMatrix_[sampleID+7][5] = upperBounds_[5];
printf("check 2\n");
   }
   return 0;
}

// ************************************************************************
// refine the sample space
// ------------------------------------------------------------------------
int BoxBehnkenSampling::refine(int ratio, int randomize, double thresh,
                               int nSamples, double *sampleErrors)
{
   (void) ratio;
   (void) randomize;
   (void) thresh;
   (void) nSamples;
   (void) sampleErrors;
   printf("BoxBehnkenSampling::refine ERROR - not available.\n");
   exit(1);
   return 0;
}

