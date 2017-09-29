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
// Functions for the Latin hypercube class 
// AUTHOR : CHARLES TONG
// DATE   : 2003
// ************************************************************************
#include <stdlib.h>
#include <string.h>
#include <math.h>
using namespace std;

#include "Main/Psuade.h"
#include "Util/sysdef.h"
#include "Util/PsuadeUtil.h"
#include "LHSampling.h"
#define PABS(x) ((x) > 0 ? (x) : -(x))

// ************************************************************************
// external functions
// ************************************************************************
extern "C" 
{
   void OA_strength(int q,int nrow,int ncol,int** A,int *str,int verbose);
}

// ************************************************************************
// Constructor 
// ------------------------------------------------------------------------
LHSampling::LHSampling() : Sampling()
{
   samplingID_ = PSUADE_SAMP_LHS;
   trueRandom_ = 0;
}

// *******************************************************************
// destructor 
// -------------------------------------------------------------------
LHSampling::~LHSampling()
{
}

// ************************************************************************
// initialize the sampling data
// ------------------------------------------------------------------------
int LHSampling::initialize(int initLevel)
{
   int    ss, jj, ii, repID, *iArray1, *iArray2;
   int    nReps, strength, **permMatrix=NULL, **storeMatrix, index;
   int    ir, maxMinDist, minDist, dist, ss2, ntimes;
   double *ranges, scale, ddata, **perturbMatrix=NULL;
   char   pString[1000];

   if (nInputs_ == 0 || lowerBounds_ == NULL || upperBounds_ == NULL)
   {
      printf("LHSampling::initialize ERROR - input not set up.\n");
      exit(1);
   }

   deleteSampleData();
   nReps = nReplications_;
   trueRandom_ = 0;
   if (randomize_ & 2) trueRandom_ = 1;
   if (nSamples_ / nReps * nReps != nSamples_) 
   {
      printf("LHSampling : nSamples must be multiples of replications.\n");
      exit(1);
   }
   nSymbols_ = nSamples_ / nReps;
   if (initLevel != 0) return 0;

   if (printLevel_ > 4)
   {
      printf("LHSampling::initialize: nSamples = %d\n", nSamples_);
      printf("LHSampling::initialize: nInputs  = %d\n", nInputs_);
      printf("LHSampling::initialize: nOutputs = %d\n", nOutputs_);
      if (randomize_ != 0)
           printf("LHSampling::initialize: randomize on\n");
      else printf("LHSampling::initialize: randomize off\n");
      if (trueRandom_ != 0)
           printf("LHSampling::initialize: more randomize on\n");
      else printf("LHSampling::initialize: more randomize off\n");
      for (ii = 0; ii < nInputs_; ii++)
         printf("    LHSampling input %3d = [%e %e]\n", ii+1,
                lowerBounds_[ii], upperBounds_[ii]);
      if (inputNumSettings_ != NULL)
      {
         for (ii = 0; ii < nInputs_; ii++)
         {
            if (inputNumSettings_[ii] != nSymbols_)
            {
               printf("LHSampling ERROR: inputSetting not compabible.\n");
               exit(1);
            }
         }
      }
   }

   if (nSymbols_ == 1)
   {
      allocSampleData();
      for (ss = 0; ss < nSamples_; ss++)
      {
         if (trueRandom_ != 0)
         {
            for (ii = 0; ii < nInputs_; ii++) 
            {
               ddata = upperBounds_[ii] - lowerBounds_[ii];
               sampleMatrix_[ss][ii] = PSUADE_drand() * ddata +
                                       lowerBounds_[ii];
            }
         }
         else
         {
            for (ii = 0; ii < nInputs_; ii++)
               sampleMatrix_[ss][ii] = 0.5 * 
                           (upperBounds_[ii] + lowerBounds_[ii]);
         }
      }
      return 0;
   }

   permMatrix = new int*[nSamples_];
   storeMatrix = new int*[nSamples_];
   for (ss = 0; ss < nSamples_; ss++) 
   {
      permMatrix[ss] = new int[nInputs_];
      storeMatrix[ss] = new int[nInputs_];
   }
   for (ss = 0; ss < nSamples_; ss+=nSymbols_) 
      for (jj = 0; jj < nSymbols_; jj++) 
         for (ii = 0; ii < nInputs_; ii++) permMatrix[ss+jj][ii] = jj;

   iArray1 = new int[nSymbols_];
   iArray2 = new int[nSymbols_];
   maxMinDist = 0;
   if      (nSamples_ >= 10000) ntimes = 1;
   else if (nSamples_ >= 9000)  ntimes = 2;
   else if (nSamples_ >= 4000)  ntimes = 5;
   else if (nSamples_ >= 1000)  ntimes = 10;
   else                         ntimes = 50;
   
   if (psSamExpertMode_ == 1 && printLevel_ > 0)
   {
      sprintf(pString,
          "LHSampling: number of times to search for maxi-min (1-100): \n");
      ntimes = getInt(1, 1000, pString);
   }

   for (ir = 0; ir < ntimes; ir++)
   {
      if (printLevel_ >= 5 && (ir % (ntimes/10+1) == 0)) 
         printf("LHSampling::initialize - max-min cycle %d (out of %d).\n",
                ir+1, ntimes); 
      for (ss = 0; ss < nSamples_; ss += nSymbols_) 
      {
         for (ii = 0; ii < nInputs_; ii++) 
         { 
            generateRandomIvector(nSymbols_, iArray1);
            for (jj = 0; jj < nSymbols_; jj++) 
               iArray2[jj] = permMatrix[ss+iArray1[jj]][ii];
            for (jj = 0; jj < nSymbols_; jj++) 
               permMatrix[ss+jj][ii] = iArray2[jj];
         }
      }
      if (ntimes == 1)
      {
         for (ss = 0; ss < nSamples_; ss++) 
            for (ii = 0; ii < nInputs_; ii++) 
               storeMatrix[ss][ii] = permMatrix[ss][ii];
         break;
      }
      minDist = nInputs_ * nSymbols_;
      for (ss = 0; ss < nSamples_; ss++) 
      {
         for (ss2 = ss+1; ss2 < nSamples_; ss2++) 
         {
            dist = 0;
            for (ii = 0; ii < nInputs_; ii++) 
               dist += PABS(permMatrix[ss][ii] - permMatrix[ss2][ii]);
            if (dist > 0 && dist < minDist) minDist = dist;
         }
      }
      if (minDist > maxMinDist) 
      {
         for (ss = 0; ss < nSamples_; ss++) 
            for (ii = 0; ii < nInputs_; ii++) 
               storeMatrix[ss][ii] = permMatrix[ss][ii];
         maxMinDist = minDist; 
         if (printLevel_ >= 1)
         {
            printf("LHSampling::current max-min distance    = %d (sample %d)\n",
                   maxMinDist, ir+1);
            printf("            each input has max distance = %d\n",
                   nSamples_-1);
            printf("            number of inputs            = %d\n",nInputs_);
         }
      }
   }
   for (ss = 0; ss < nSamples_; ss++) delete [] permMatrix[ss];
   delete [] permMatrix;
   permMatrix = storeMatrix;
   delete [] iArray1;
   delete [] iArray2;

   if (nSamples_ < 1000)
   {
      if (printLevel_ >= 1 && ntimes > 1) 
         printf("LHSampling::initialize - checking.\n");
      for (repID = 0; repID < nReps; repID++)
      {
         OA_strength(nSymbols_, nSymbols_, nInputs_, 
                     &(permMatrix[repID*nSymbols_]), &strength, 0);
         if (strength < 1)
         {
            printf("LHSampling : failed strength test (%d,%d).\n",
                   strength,repID);
            printf("   ==> Please consult PSUADE developers.\n");
            exit(1);
         }
      }
   }

   allocSampleData();
   ranges = new double[nInputs_];
   for (ii = 0; ii < nInputs_; ii++) 
      ranges[ii] = upperBounds_[ii] - lowerBounds_[ii];


   if (trueRandom_ != 0)
   {
      scale = 1.0 / ((double) nSymbols_);
      for (ss = 0; ss < nSamples_; ss++) 
      {
         for (ii = 0; ii < nInputs_; ii++)
         {
            index = permMatrix[ss][ii];
            ddata = (PSUADE_drand() + index) * scale;
            sampleMatrix_[ss][ii] = ddata * ranges[ii] + lowerBounds_[ii];
         }
      }
   } 

   else if (randomize_ != 0)
   {
      perturbMatrix = new double*[nInputs_];
      for (ii = 0; ii < nInputs_; ii++)
      {
         perturbMatrix[ii] = new double[nSymbols_];
         for (jj = 0; jj < nSymbols_; jj++)
            perturbMatrix[ii][jj] = PSUADE_drand() - 0.5;
      }
      scale = 1.0 / ((double) nSymbols_);
      for (ss = 0; ss < nSamples_; ss++)
      {
         for (ii = 0; ii < nInputs_; ii++)
         {
            index = permMatrix[ss][ii];
            ddata = (perturbMatrix[ii][index] + index + 0.5) * scale;
            sampleMatrix_[ss][ii] = ddata * ranges[ii] + lowerBounds_[ii];
         }
      }
   }

   else 
   {
      scale = 1.0 / ((double) (nSymbols_ - 1));
      for (ss = 0; ss < nSamples_; ss++) 
      {
         for (ii = 0; ii < nInputs_; ii++)
         {
            if (inputSettings_ != NULL && inputSettings_[ii] != NULL)
            {
               index = permMatrix[ss][ii];
               sampleMatrix_[ss][ii] = inputSettings_[ii][index];
            }
            else
            {
               ddata = (double) permMatrix[ss][ii];
               ddata = ddata * scale;
               sampleMatrix_[ss][ii] = ddata * ranges[ii] + 
                                       lowerBounds_[ii];
            }
         }
      }
   }

   delete [] ranges;
   if (permMatrix != NULL)
   {
      for (ss = 0; ss < nSamples_; ss++) delete [] permMatrix[ss];
      delete [] permMatrix;
   }
   if (perturbMatrix != NULL)
   {
      for (ii = 0; ii < nInputs_; ii++) delete [] perturbMatrix[ii];
      delete [] perturbMatrix;
   }
   return 0;
}

// ************************************************************************
// refine the sample space
// ------------------------------------------------------------------------
int LHSampling::refine(int refineRatio, int randomize, double thresh,
                       int nSamples, double *sampleErrors)
{
   int    ss, ii, jj, newNSymbols, **permMatrix=NULL;   
   int    *iArray1, *iArray2, **binFlags, sampleOffset, index, nReps;
   int    newSampleOffset, repID, strength, *newStates, outputID;
   int    binCount, addNSymbols, nLevels;
   double scale, **newSamples, *ranges, **bounds, **perturbMatrix=NULL;
   double *newOutputs, ddata;

   (void) randomize;
   (void) thresh;
   (void) nSamples;
   (void) sampleErrors;

   nLevels = refineRatio;
   ranges = new double[nInputs_];
   for (ii = 0; ii < nInputs_; ii++) 
      ranges[ii] = upperBounds_[ii] - lowerBounds_[ii];
   nReps = nSamples_ / nSymbols_;

   if (randomize_ != 0 || trueRandom_ != 0)
   {
      newSamples = new double*[nSamples_*nLevels];
      for (ss = 0; ss < nSamples_*nLevels; ss++) 
         newSamples[ss] = new double[nInputs_];
      newOutputs = new double[nSamples_*nLevels*nOutputs_];
      newStates = new int[nSamples_*nLevels];
      for (ss = 0; ss < nSamples_*nLevels; ss++)
      {
         newStates[ss] = 0;
         for (outputID = 0; outputID < nOutputs_; outputID++)
            newOutputs[ss*nOutputs_+outputID] = PSUADE_UNDEFINED;
      }

      newNSymbols = nSymbols_ * nLevels;
      bounds = new double*[nInputs_];
      for (ii = 0; ii < nInputs_; ii++) 
      {
         bounds[ii] = new double[newNSymbols+1];
         for (jj = 0; jj <= newNSymbols; jj++) 
            bounds[ii][jj] = ranges[ii] / newNSymbols * 
                             jj + lowerBounds_[ii];
      }

      perturbMatrix = new double*[nInputs_];
      for (ii = 0; ii < nInputs_; ii++)
      {
         perturbMatrix[ii] = new double[nSymbols_*nLevels];
         for (jj = 0; jj < nSymbols_*nLevels; jj++)
            perturbMatrix[ii][jj] = PSUADE_drand();
      }

      for (ss = 0; ss < nSamples_; ss++)
      {
         for (ii = 0; ii < nInputs_; ii++) 
         {
            ddata = sampleMatrix_[ss][ii]; 
            if (ddata == bounds[ii][newNSymbols]) jj = newNSymbols;
            else
            {
               for (jj = 1; jj <= newNSymbols; jj++) 
                  if (ddata < bounds[ii][jj]) break;
            }
            jj--;
            if (jj >= newNSymbols)
            {
               printf("LHSampling::refine ERROR (3) - %d %d %e\n",jj,
                      newNSymbols,ddata);
               exit(1);
            }
            ddata -= lowerBounds_[ii];
            ddata  = ddata / ranges[ii] * (double) newNSymbols;
            ddata -= (double) jj;
            perturbMatrix[ii][jj] = ddata;
         }
      }

      permMatrix  = new int*[nSymbols_];
      for (jj = 0; jj < nSymbols_; jj++) 
         permMatrix[jj] = new int[nInputs_];
      iArray1 = new int[nSymbols_];
      iArray2 = new int[nSymbols_];
      binFlags = new int*[nInputs_];
      for (ii = 0; ii < nInputs_; ii++) 
         binFlags[ii] = new int[newNSymbols];

      sampleOffset = 0;
      newSampleOffset = 0;
      for (repID = 0; repID < nReps; repID++)
      {
         for (jj = 0; jj < nSymbols_; jj++) 
         {
            for (ii = 0; ii < nInputs_; ii++)
               newSamples[newSampleOffset+jj][ii] = 
                          sampleMatrix_[sampleOffset+jj][ii];
            for (outputID = 0; outputID < nOutputs_; outputID++)
               newOutputs[(newSampleOffset+jj)*nOutputs_+outputID] = 
                  sampleOutput_[(sampleOffset+jj)*nOutputs_+outputID];
            newStates[newSampleOffset+jj] = 
                        sampleStates_[sampleOffset+jj]; 
         }
         newSampleOffset += nSymbols_;

         for (ii = 0; ii < nInputs_; ii++) 
            for (jj = 0; jj < newNSymbols; jj++)
               binFlags[ii][jj] = 0;

         for (ss = sampleOffset; ss < sampleOffset+nSymbols_; ss++) 
         {
            for (ii = 0; ii < nInputs_; ii++) 
            {
               ddata = sampleMatrix_[ss][ii]; 
               if (ddata == bounds[ii][newNSymbols])
                  jj = newNSymbols;
               else
               {
                  for (jj = 1; jj < (newNSymbols+1); jj++) 
                     if (ddata < bounds[ii][jj]) break;
               }
               binFlags[ii][jj-1] = -1;
            }
         }

         for (ii = 0; ii < nInputs_; ii++)
         {
            for (jj = 0; jj < newNSymbols; jj++) 
               if (binFlags[ii][jj] == 0)
                  binFlags[ii][jj] = jj;
            binCount = 0;
            for (jj = 0; jj < newNSymbols; jj++) 
               if (binFlags[ii][jj] >= 0)
                  binFlags[ii][binCount++] = binFlags[ii][jj];
         }

         for (jj = 0; jj < binCount; jj++) 
            for (ii = 0; ii < nInputs_; ii++) 
               permMatrix[jj][ii] = binFlags[ii][jj];
         for (ii = 0; ii < nInputs_; ii++) 
         {
            generateRandomIvector(binCount, iArray1);
            for (jj = 0; jj < binCount; jj++)
               iArray2[jj] = permMatrix[iArray1[jj]][ii];
            for (jj = 0; jj < binCount; jj++)
               permMatrix[jj][ii] = iArray2[jj];
         }

         if (trueRandom_ != 0)
         {
            scale = 1.0 / ((double) (nSymbols_*nLevels));
            for (ss = 0; ss < binCount; ss++) 
            {
               for (ii = 0; ii < nInputs_; ii++)
               {
                  index =  permMatrix[ss][ii];
                  ddata = (PSUADE_drand() + index) * scale;
                  newSamples[newSampleOffset+ss][ii] = 
                            ddata * ranges[ii] + lowerBounds_[ii];
               }
            }
         }
         else
         {
            scale = 1.0 / ((double) (nSymbols_*nLevels));
            for (ss = 0; ss < binCount; ss++) 
            {
               for (ii = 0; ii < nInputs_; ii++)
               {
                  index =  permMatrix[ss][ii];
                  ddata = (perturbMatrix[ii][index] + index) * scale;
                  newSamples[newSampleOffset+ss][ii] = 
                            ddata * ranges[ii] + lowerBounds_[ii];
               }
            }
         }

         newSampleOffset += binCount;
         sampleOffset += nSymbols_;
      }

      deleteSampleData();
      delete [] iArray1;
      delete [] iArray2;
      for (jj = 0; jj < nSymbols_; jj++) delete [] permMatrix[jj];
      delete [] permMatrix;
      for (ii = 0; ii < nInputs_; ii++) delete [] binFlags[ii];
      delete [] binFlags;
      for (ii = 0; ii < nInputs_; ii++) delete [] bounds[ii];
      delete [] bounds;

      nSamples_ = nSamples_ * nLevels;
      nSymbols_ = nSymbols_ * nLevels;
      sampleMatrix_ = newSamples;
      sampleOutput_ = newOutputs;
      sampleStates_ = newStates;
   }
   else
   {

      newNSymbols = (nSymbols_ - 1) * nLevels + 1;
      addNSymbols = newNSymbols - nSymbols_;
      newSamples = new double*[newNSymbols*nReps];
      for (ss = 0; ss < newNSymbols*nReps; ss++) 
         newSamples[ss] = new double[nInputs_];
      newOutputs = new double[newNSymbols*nReps*nOutputs_];
      newStates = new int[newNSymbols*nReps];
      for (ss = 0; ss < newNSymbols*nReps; ss++)
      {
         newStates[ss] = 0;
         for (outputID = 0; outputID < nOutputs_; outputID++)
            newOutputs[ss*nOutputs_+outputID] = PSUADE_UNDEFINED;
      }

      permMatrix  = new int*[nSymbols_*nLevels];
      for (jj = 0; jj < nSymbols_*nLevels; jj++) 
         permMatrix[jj] = new int[nInputs_];
      iArray1 = new int[nSymbols_];
      iArray2 = new int[nSymbols_];

      sampleOffset = 0;
      newSampleOffset = 0;
      for (repID = 0; repID < nReps; repID++)
      {
         for (jj = 0; jj < nSymbols_; jj++) 
         {
            for (ii = 0; ii < nInputs_; ii++)
               newSamples[newSampleOffset+jj][ii] = 
                          sampleMatrix_[sampleOffset+jj][ii];
            for (outputID = 0; outputID < nOutputs_; outputID++)
               newOutputs[(newSampleOffset+jj)*nOutputs_+outputID] = 
                  sampleOutput_[(sampleOffset+jj)*nOutputs_+outputID];
            newStates[newSampleOffset+jj] = 
                        sampleStates_[sampleOffset+jj]; 
         }
         newSampleOffset += nSymbols_;

         for (ii = 0; ii < nInputs_; ii++)
         {
            addNSymbols = 0;
            for (jj = 0; jj < newNSymbols; jj++) 
               if ((jj % nLevels) != 0)
                  permMatrix[addNSymbols++][ii] = jj;
            generateRandomIvector(addNSymbols, iArray1);
            for (jj = 0; jj < addNSymbols; jj++)
               iArray2[jj] = permMatrix[iArray1[jj]][ii];
            for (jj = 0; jj < addNSymbols; jj++)
               permMatrix[jj][ii] = iArray2[jj];
         }

         scale = 1.0 / ((double) (newNSymbols - 1));
         for (jj = 0; jj < addNSymbols; jj++) 
         {
            for (ii = 0; ii < nInputs_; ii++)
            {
               ddata = (double) permMatrix[jj][ii];
               ddata = ddata * scale;
               newSamples[newSampleOffset+jj][ii] = 
                       ddata * ranges[ii] + lowerBounds_[ii];
            }
         }
         newSampleOffset += addNSymbols;
         sampleOffset += nSymbols_;
      }

      deleteSampleData();
      delete [] iArray1;
      delete [] iArray2;
      for (jj = 0; jj < nSymbols_; jj++) delete [] permMatrix[jj];
      delete [] permMatrix;
      nSymbols_ = newNSymbols;
      nSamples_ = nSymbols_ * nReps;
      sampleMatrix_ = newSamples;
      sampleOutput_ = newOutputs;
      sampleStates_ = newStates;
   }

   if (printLevel_ > 4)
   {
      printf("LHSampling refine: nSamples = %d\n", nSamples_);
      printf("LHSampling refine: nInputs  = %d\n", nInputs_);
      printf("LHSampling refine: nOutputs = %d\n", nOutputs_);
      if (randomize_ != 0)
           printf("LHSampling refine: randomize on\n");
      else printf("LHSampling refine: randomize off\n");
      if (trueRandom_ != 0)
           printf("LHSampling refine: more randomize on\n");
      else printf("LHSampling refine: more randomize off\n");
      for (ii = 0; ii < nInputs_; ii++)
         printf("    LHSampling input %3d = [%e %e]\n", ii+1,
                lowerBounds_[ii], upperBounds_[ii]);
      if (inputNumSettings_ != NULL || symTable_ != NULL)
         printf("LHSampling refine: diable input settings, symbol table.\n");
   }

   if (nSamples_ < 1000)
   {
      bounds = new double*[nInputs_];
      for (ii = 0; ii < nInputs_; ii++) 
      {
         bounds[ii] = new double[nSymbols_+1];
         for (jj = 0; jj <= nSymbols_; jj++) 
            bounds[ii][jj] = ranges[ii] / nSymbols_ * 
                             jj + lowerBounds_[ii];
      }
      permMatrix = new int*[nSamples_/nReps];
      for (ss = 0; ss < nSamples_/nReps; ss++)
         permMatrix[ss] = new int[nInputs_];

      sampleOffset = 0;
      for (repID = 0; repID < nReps; repID++)
      {
         for (ss = 0; ss < nSamples_/nReps; ss++) 
         {
#if 0
            printf("sample %4d : ", ss);
#endif
            for (ii = 0; ii < nInputs_; ii++) 
            {
               ddata = sampleMatrix_[sampleOffset+ss][ii]; 
               if (ddata == bounds[ii][nSymbols_]) jj = nSymbols_;
               else
               {
                  for (jj = 1; jj < (nSymbols_+1); jj++) 
                     if (ddata < bounds[ii][jj]) break;
               }
               permMatrix[ss][ii] = jj-1;
#if 0
               printf(" %4d ", jj-1);
#endif
            }
#if 0
            printf("\n");
#endif
         }
         OA_strength(nSamples_/nReps, nSymbols_, nInputs_,
                     permMatrix, &strength, 0);
         if (strength != 1)
            printf("LHS refine ERROR : replication %d : OA_strength = %d\n",
                   repID, strength);
         sampleOffset += (nSamples_/nReps);
      }
      for (ss = 0; ss < nSamples_/nReps; ss++)
         delete [] permMatrix[ss];
      delete [] permMatrix;
      for (ii = 0; ii < nInputs_; ii++)
          delete [] bounds[ii];
      delete [] bounds;
   }

   delete [] ranges;
   if (perturbMatrix != NULL)
   {
      for (ii = 0; ii < nInputs_; ii++) delete [] perturbMatrix[ii];
      delete [] perturbMatrix;
   }
   return 0;
}

// ************************************************************************
// set input settings
// ------------------------------------------------------------------------
int LHSampling::setInputParams(int nInputs, int *counts, double **settings,
                               int *symtable)
{
   int ii, inputCnt, ss;

   if (nInputs_ != 0 && nInputs != nInputs_) 
   { 
      printf("LHSampling::setInputParams - nInputs mismatch.\n");
      exit(1);
   }
   nInputs_ = nInputs;
   if (symtable != NULL)
   {
      symTable_ = new int[nInputs];
      for (ii = 0; ii < nInputs_; ii++) symTable_[ii] = symtable[ii];
   }
   if (counts != NULL)
   {
      inputCnt = 0;
      for (ii = 0; ii < nInputs_; ii++)
      {
         if (counts[ii] != 0 && counts[ii] != nSymbols_)
         {
            printf("LHSampling::setInputParams - counts mismatch.\n");
            printf("            count data = %d %d\n", nSymbols_,
                   counts[ii]);
            exit(1);
         }
         else if (counts[ii] == nSymbols_) inputCnt++;
      }
      if (inputCnt > 0)
      {
         inputSettings_ = new double*[nInputs_];
         inputNumSettings_ = new int[nInputs_];
         for (ii = 0; ii < nInputs_; ii++)
         {
            if (counts[ii] == nSymbols_)
            {
               inputSettings_[ii] = new double[counts[ii]];
               inputNumSettings_[ii] = counts[ii];
               for (ss = 0; ss < counts[ii]; ss++)
                  inputSettings_[ii][ss] = settings[ii][ss];
            }
            else
            {
               inputSettings_[ii] = NULL;
               inputNumSettings_[ii] = 0;
            }
         } 
      }
   }
   return 0;
}

