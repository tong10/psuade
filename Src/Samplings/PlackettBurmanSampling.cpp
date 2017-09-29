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
// Functions for the Plackett-Burman sampling class 
// AUTHOR : CHARLES TONG
// DATE   : 2005
// ************************************************************************
#include <stdio.h>
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "PlackettBurmanSampling.h"

// ************************************************************************
// generating vectors for the Plackett Burman designs
// ------------------------------------------------------------------------
static int PSUADE_PB_GV12[12] = {1,1,0,1,1,1,0,0,0,1,0};
static int PSUADE_PB_GV20[20] = {1,1,0,0,1,1,1,1,0,1,0,1,0,0,0,0,1,1,0};
static int PSUADE_PB_GV24[24] = 
{1,1,1,1,1,0,1,0,1,1,0,0,1,1,0,0,1,0,1,0,0,0,0};
static int PSUADE_PB_GV36[36] = 
{0,1,0,1,1,1,0,0,0,1,1,1,1,1,0,1,1,1,0,0,1,0,0,0,0,1,0,1,0,1,1,0,0,1,0};
static int PSUADE_PB_GV44[44] = 
{1,1,0,0,1,0,1,0,0,1,1,1,0,1,1,1,1,1,0,0,0,1,0,1,1,1,0,0,0,0,0,1,0,0,0,
 1,1,0,1,0,1,1,0};

static int PSUADE_PB_GV28X[9][9] = 
{
   {1,0,1,1,1,1,0,0,0}, {1,1,0,1,1,1,0,0,0}, {0,1,1,1,1,1,0,0,0},
   {0,0,0,1,0,1,1,1,1}, {0,0,0,1,1,0,1,1,1}, {0,0,0,0,1,1,1,1,1},
   {1,1,1,0,0,0,1,1,0}, {1,1,1,0,0,0,1,1,0}, {1,1,1,0,0,0,0,1,1},
};

static int PSUADE_PB_GV28Y[9][9] = 
{
   {0,1,0,0,0,1,0,0,1}, {0,0,1,1,0,0,1,0,0}, {1,0,0,0,1,0,0,1,0},
   {0,0,1,0,1,0,0,1,0}, {1,0,0,0,0,1,1,0,0}, {0,1,0,1,0,0,0,1,0},
   {0,0,1,0,0,1,0,1,0}, {1,0,0,1,0,0,0,0,1}, {0,1,0,0,1,0,1,0,0},
};

static int PSUADE_PB_GV28Z[9][9] = 
{
   {1,1,0,1,0,1,1,0,1}, {0,1,1,1,1,0,1,1,0}, {1,0,1,0,1,1,0,1,1},
   {1,0,1,1,1,0,1,0,1}, {1,1,0,0,1,1,1,1,0}, {0,1,1,1,0,1,0,1,1},
   {1,0,1,1,0,1,1,1,0}, {1,1,0,1,1,0,0,1,1}, {0,1,1,0,1,1,1,0,1},
};

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
PlackettBurmanSampling::PlackettBurmanSampling() : Sampling()
{
   samplingID_ = PSUADE_SAMP_PBD;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
PlackettBurmanSampling::~PlackettBurmanSampling()
{
}

// ************************************************************************
// initialize the sampler data
// ------------------------------------------------------------------------
int PlackettBurmanSampling::initialize(int initLevel)
{
   int ii, jj, kk, ipower, whichPlan, checkN;
   int **patternMatrix, **tempMatrix, rowCnt, *genVec;

   if (nSamples_ = 0)
   {
      printf("PlackettBurmanSampling::initialize ERROR - nSamples = 0.\n");
      exit(1);
   }
   if (nInputs_ == 0 || lowerBounds_ == NULL || upperBounds_ == NULL)
   {
      printf("PlackettBurmanSampling::initialize ERROR - input not set up.\n");
      exit(1);
   }
   if (nInputs_ < 2)
   {
      printf("PlackettBurmanSampling::initialize ERROR - nInputs < 2.\n");
      exit(1);
   }
   if (nInputs_ > 47)
   {
      printf("PlackettBurmanSampling::initialize ERROR - nInputs > 47 not");
      printf(" supported.\n");
      exit(1);
   }
   

   deleteSampleData();
   if (initLevel != 0) return 0;
   ii = nInputs_ / 4;
   nSamples_ = (ii + 1) * 4;
   allocSampleData();

   patternMatrix = new int*[nSamples_];
   for (ii = 0; ii < nSamples_; ii++)
      patternMatrix[ii] = new int[nSamples_];
   whichPlan = -1;
   ipower = int (log((double) nSamples_) / log(2.0e0) + 1.0e-8);
   checkN = (int) pow(2.0, (double) ipower);
   if (checkN == nSamples_) whichPlan = 0;
   else if (nSamples_ == 12) whichPlan = 1;
   else if (nSamples_ == 20) whichPlan = 2;
   else if (nSamples_ == 24) whichPlan = 3;
   else if (nSamples_ == 28) whichPlan = 4;
   else if (nSamples_ == 36) whichPlan = 5;
   else if (nSamples_ == 40) whichPlan = 6;
   else if (nSamples_ == 44) whichPlan = 7;
   else if (nSamples_ == 48) whichPlan = 8;

   if (whichPlan == 0)
   {
      patternMatrix[0][0] =  1; patternMatrix[0][1] = 1;
      patternMatrix[0][2] =  1; patternMatrix[0][3] = 1;
      patternMatrix[1][0] =  1; patternMatrix[1][1] = -1;
      patternMatrix[1][2] =  1; patternMatrix[1][3] = -1;
      patternMatrix[2][0] =  1; patternMatrix[2][1] = 1;
      patternMatrix[2][2] = -1; patternMatrix[2][3] = -1;
      patternMatrix[3][0] =  1; patternMatrix[3][1] = -1;
      patternMatrix[3][2] = -1; patternMatrix[3][3] = 1;
      if (nSamples_ > 4)
      {
         ipower = int (log((double) nSamples_) / log(2.0e0) + 1.0e-8);
         tempMatrix = new int*[nSamples_];
         for (ii = 0; ii < nSamples_; ii++)
            tempMatrix[ii] = new int[nSamples_];
         rowCnt = 4;
         for (ii = 0; ii < ipower-2; ii++)
         {
            for (jj = 0; jj < rowCnt; jj++)
               for (kk = 0; kk < rowCnt; kk++)
                  tempMatrix[jj][kk] = patternMatrix[jj][kk];
            for (jj = 0; jj < rowCnt; jj++)
            {
               for (kk = 0; kk < rowCnt; kk++)
               {
                  patternMatrix[rowCnt+jj][kk] = tempMatrix[jj][kk];
                  patternMatrix[jj][rowCnt+kk] = tempMatrix[jj][kk];
                  patternMatrix[rowCnt+jj][rowCnt+kk] = -tempMatrix[jj][kk];
               }
            }
            rowCnt *= 2;
         }
         for (ii = 0; ii < nSamples_; ii++) delete [] tempMatrix[ii];
         delete [] tempMatrix;
      }

      for (ii = 0; ii < nSamples_; ii++)
         for (jj = 0; jj < nSamples_-1; jj++)
            patternMatrix[ii][jj] = patternMatrix[ii][jj+1];
   }
   else if (whichPlan != 4 && whichPlan != 6 && whichPlan != 8)
   {
      switch(whichPlan)
      {
         case 1: genVec = PSUADE_PB_GV12; break;
         case 2: genVec = PSUADE_PB_GV20; break;
         case 3: genVec = PSUADE_PB_GV24; break;
         case 5: genVec = PSUADE_PB_GV36; break;
         case 7: genVec = PSUADE_PB_GV44; break;
      }

      for (jj = 0; jj < nSamples_-1; jj++)
         patternMatrix[0][jj] = genVec[jj];

      for (ii = 1; ii < nSamples_-1; ii++)
      {
         for (jj = 1; jj < nSamples_-1; jj++)
            patternMatrix[ii][jj] = patternMatrix[ii-1][jj-1];
         patternMatrix[ii][0] = patternMatrix[ii-1][nSamples_-2];
      }

      for (jj = 0; jj < nSamples_-1; jj++)
         patternMatrix[nSamples_-1][jj] = 0;
   }
   else if (whichPlan == 4)
   {
      for (ii = 0; ii < 9; ii++)
      {
         for (jj = 0; jj < 9; jj++)
         {
            patternMatrix[ii][jj]       = PSUADE_PB_GV28X[ii][jj]; 
            patternMatrix[ii+9][jj+9]   = PSUADE_PB_GV28X[ii][jj]; 
            patternMatrix[ii+18][jj+18] = PSUADE_PB_GV28X[ii][jj]; 
            patternMatrix[ii][jj+9]     = PSUADE_PB_GV28Y[ii][jj]; 
            patternMatrix[ii+9][jj+18]  = PSUADE_PB_GV28Y[ii][jj]; 
            patternMatrix[ii+18][jj]    = PSUADE_PB_GV28Y[ii][jj]; 
            patternMatrix[ii][jj+18]    = PSUADE_PB_GV28Z[ii][jj]; 
            patternMatrix[ii+9][jj]     = PSUADE_PB_GV28Z[ii][jj]; 
            patternMatrix[ii+18][jj+9]  = PSUADE_PB_GV28Z[ii][jj]; 
         }
      }
      for (jj = 0; jj < 27; jj++) patternMatrix[27][jj] = 0;
   }
   else if (whichPlan == 6)
   {
      genVec = PSUADE_PB_GV20;
      for (jj = 1; jj < 20; jj++) patternMatrix[0][jj] = genVec[jj];
      for (ii = 0; ii < 19; ii++)
      {
         patternMatrix[ii][0] = 1;
         for (jj = 2; jj < 20; jj++)
            patternMatrix[ii][jj] = patternMatrix[ii-1][jj-1];
         patternMatrix[ii][1] = patternMatrix[ii-1][19];
      }
      for (jj = 1; jj < 20; jj++) patternMatrix[19][jj] = 0;
      patternMatrix[19][0] = 1;

      for (ii = 0; ii < 20; ii++) 
      {
         for (jj = 0; jj < 20; jj++) 
         {
            patternMatrix[ii+20][jj] = patternMatrix[ii][jj]; 
            patternMatrix[ii][jj+20] = patternMatrix[ii][jj]; 
            patternMatrix[ii+20][jj+20] = -patternMatrix[ii][jj]; 
         }
      }

      for (ii = 0; ii < nSamples_; ii++)
         for (jj = 0; jj < nSamples_-1; jj++)
            patternMatrix[ii][jj] = patternMatrix[ii][jj+1];
   } 
   else if (whichPlan == 8)
   {
      genVec = PSUADE_PB_GV24;
      for (jj = 1; jj < 24; jj++) patternMatrix[0][jj] = genVec[jj];
      for (ii = 0; ii < 23; ii++)
      {
         patternMatrix[ii][0] = 1;
         for (jj = 2; jj < 24; jj++)
            patternMatrix[ii][jj] = patternMatrix[ii-1][jj-1];
         patternMatrix[ii][1] = patternMatrix[ii-1][23];
      }
      for (jj = 1; jj < 24; jj++) patternMatrix[23][jj] = 0;
      patternMatrix[23][0] = 1;

      for (ii = 0; ii < 24; ii++) 
      {
         for (jj = 0; jj < 24; jj++) 
         {
            patternMatrix[ii+24][jj] = patternMatrix[ii][jj]; 
            patternMatrix[ii][jj+24] = patternMatrix[ii][jj]; 
            patternMatrix[ii+24][jj+24] = -patternMatrix[ii][jj]; 
         }
      }

      for (ii = 0; ii < nSamples_; ii++)
         for (jj = 0; jj < nSamples_-1; jj++)
            patternMatrix[ii][jj] = patternMatrix[ii][jj+1];
   }
      
   for (ii = 0; ii < nSamples_; ii++)
      for (jj = 0; jj < nInputs_; jj++)
         sampleMatrix_[ii][jj] = (upperBounds_[jj] - lowerBounds_[jj]) *
               0.5 * (patternMatrix[ii][jj] + 1.0) + lowerBounds_[jj];

   for (ii = 0; ii < nSamples_; ii++) delete [] patternMatrix[ii];
   delete [] patternMatrix;

   if (printLevel_ > 4)
   {
      printf("PlackettBurmanSampling::initialize: nSamples = %d\n", nSamples_);
      printf("PlackettBurmanSampling::initialize: nInputs  = %d\n", nInputs_);
      printf("PlackettBurmanSampling::initialize: nOutputs = %d\n", nOutputs_);
      for (jj = 0; jj < nInputs_; jj++)
         printf("    PlackettBurmanSampling input %3d = [%e %e]\n", jj+1,
                lowerBounds_[jj], upperBounds_[jj]);
   }
   return 0;
}

// ************************************************************************
// refine the sample space
// ------------------------------------------------------------------------
int PlackettBurmanSampling::refine(int, int, double, int, double *)
{
   printf("PlackettBurmanSampling::refine ERROR - not supported.\n");
   exit(1);
   return 0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
PlackettBurmanSampling& PlackettBurmanSampling::operator=
                                     (const PlackettBurmanSampling &)
{
   printf("PlackettBurmanSampling operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

