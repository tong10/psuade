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
// Functions for the factorial sampling class 
// AUTHOR : CHARLES TONG
// DATE   : 2003
// ************************************************************************
#include <stdio.h>
using namespace std;

#include "Util/sysdef.h"
#include "Util/PsuadeUtil.h"
#include "FactorialSampling.h"

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
FactorialSampling::FactorialSampling() : Sampling()
{
   samplingID_ = PSUADE_SAMP_FACT;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
FactorialSampling::~FactorialSampling()
{
}

// ************************************************************************
// initialize the sampler data
// ------------------------------------------------------------------------
int FactorialSampling::initialize(int initLevel)
{
   int    inputID, sampleID, *iCounts, *iCntMax, curSymbol, *iSymbol, nsym;
   double scale, *ranges=NULL, ddata, dpower;

   if (nInputs_ == 0 || lowerBounds_ == NULL || upperBounds_ == NULL)
   {
      printf("FactorialSampling::initialize ERROR - input not set up.\n");
      exit(1);
   }

   deleteSampleData();

   if (inputSettings_ != NULL)
   {
      if (printLevel_ > 4)
         printf("FactorialSampling::initialize: inputSettings used.\n");
      if (printLevel_ > 4 && symTable_ != NULL)
         printf("FactorialSampling::initialize: symbol table overwritten.\n");
      if (symTable_ != NULL) delete [] symTable_;
      symTable_ = new int[nInputs_];
      for (inputID = 0; inputID < nInputs_; inputID++) 
         symTable_[inputID] = inputNumSettings_[inputID];
      nSamples_ = 1;
      for (inputID = 0; inputID < nInputs_; inputID++) 
         nSamples_ *= symTable_[inputID];
      if (nSamples_ == 0)
      {
         printf("FactorialSampling::initialize ERROR - ");
         printf("incomplete inputSettings.\n");
         exit(1);
      }
   } 
   else if (symTable_ == NULL)
   {
      ddata  = (double) nSamples_;
      dpower = 1.0 / (double) nInputs_;
      ddata  = pow(ddata, dpower+1.0E-12);
      nsym   = (int) ddata;
      if (symTable_ != NULL) delete [] symTable_;
      symTable_ = new int[nInputs_];
      for (inputID = 0; inputID < nInputs_; inputID++) 
         symTable_[inputID] = nsym;
      nSamples_ = 1;
      for (inputID = 0; inputID < nInputs_; inputID++) 
         nSamples_ *= symTable_[inputID];
   }
   allocSampleData();
   if (initLevel == 1) return 0;

   if (printLevel_ > 4)
   {
      printf("FactorialSampling::initialize: nSamples = %d\n", nSamples_);
      printf("FactorialSampling::initialize: nInputs  = %d\n", nInputs_);
      printf("FactorialSampling::initialize: nOutputs = %d\n", nOutputs_);
      if (randomize_ != 0)
           printf("FactorialSampling::initialize: randomize on\n");
      else printf("FactorialSampling::initialize: randomize off\n");
      for (inputID = 0; inputID < nInputs_; inputID++)
         printf("    FactorialSampling input %3d = [%e %e]\n", inputID+1,
                lowerBounds_[inputID], upperBounds_[inputID]);
      if (symTable_ != NULL)
         for (inputID = 0; inputID < nInputs_; inputID++) 
            printf("    FactorialSampling symbol table %2d = %d\n",
                   inputID+1, symTable_[inputID]);
   }

   if (nSamples_ == 1)
   {
      for (inputID = 0; inputID < nInputs_; inputID++) 
      {
         sampleMatrix_[0][inputID] = 0.5 * (lowerBounds_[inputID] +
                                            upperBounds_[inputID]);
      }
      return 0;
   }
   ranges = new double[nInputs_];
   for (inputID = 0; inputID < nInputs_; inputID++)
      ranges[inputID] = upperBounds_[inputID] - lowerBounds_[inputID];
   iCntMax = new int[nInputs_];
   iCounts = new int[nInputs_];
   iSymbol = new int[nInputs_];
   iCntMax[0] = 1;
   for (inputID = 1; inputID < nInputs_; inputID++)
      iCntMax[inputID] = iCntMax[inputID-1] * symTable_[inputID-1];
   for (inputID = 0; inputID < nInputs_; inputID++)
      iCounts[inputID] = iSymbol[inputID] = 0;

   if (randomize_ != 0)
   {
      for (sampleID = 0; sampleID < nSamples_; sampleID++)
      {
         for (inputID = 0; inputID < nInputs_; inputID++) 
         {
            curSymbol = iSymbol[inputID];
            iCounts[inputID]++;
            if (iCounts[inputID] == iCntMax[inputID])
            {
               iCounts[inputID] = 0;
               iSymbol[inputID]++;
               if (iSymbol[inputID] >= symTable_[inputID]) 
                  iSymbol[inputID] = 0;
            }
            scale = ranges[inputID] / ((double) (symTable_[inputID]));
            sampleMatrix_[sampleID][inputID] = 
                            curSymbol * scale + lowerBounds_[inputID];
            sampleMatrix_[sampleID][inputID] += (PSUADE_drand() * scale);
         }
      }
   }
   else
   {
      for (sampleID = 0; sampleID < nSamples_; sampleID++)
      {
         for (inputID = 0; inputID < nInputs_; inputID++) 
         {
            curSymbol = iSymbol[inputID];
            iCounts[inputID]++;
            if (iCounts[inputID] == iCntMax[inputID])
            {
               iCounts[inputID] = 0;
               iSymbol[inputID]++;
               if (iSymbol[inputID] >= symTable_[inputID])
                  iSymbol[inputID] = 0;
            }
            scale = ranges[inputID] / ((double) (symTable_[inputID] - 1));
            if (inputSettings_ != NULL && inputSettings_[inputID] != NULL)
               sampleMatrix_[sampleID][inputID] =
                                 inputSettings_[inputID][curSymbol]; 
            else
               sampleMatrix_[sampleID][inputID] = 
                             curSymbol * scale + lowerBounds_[inputID];
         }
      }
   }

   delete [] ranges;
   delete [] iCounts;
   delete [] iCntMax;
   delete [] iSymbol;
   return 0;
}

// ************************************************************************
// refine the sample space
// ------------------------------------------------------------------------
int FactorialSampling::refine(int refineRatio, int randomize, double thresh,
                              int nSamples, double *sampleErrors)
{
   int    newNSamples, inputID, *useArray, *iCounts, *iCntMax, *iSymbol;
   int    sampleID, index, curSymbol, *newSampleStates, outputID, nLevels;
   double *ranges=NULL, **newSampleMatrix, stepSize, ddata;
   double *newSampleOutput;

   (void) randomize;
   (void) thresh;
   (void) nSamples;
   (void) sampleErrors;

   nLevels = refineRatio;
   if (nLevels != 2)
   {
      printf("FactorialSampling::refine WARNING - nLevels set to 2.\n");
      nLevels = 2;
   }

   newNSamples = 1;
   if (randomize_ != 0)
   {
      for (inputID = 0; inputID < nInputs_; inputID++)
      {
         symTable_[inputID] *= 2;
         newNSamples *= symTable_[inputID];
      }
   }
   else
   {
      for (inputID = 0; inputID < nInputs_; inputID++)
      {
         symTable_[inputID] = symTable_[inputID] * 2 - 1;
         newNSamples *= symTable_[inputID];
      }
   }

   ranges = new double[nInputs_];
   for (inputID = 0; inputID < nInputs_; inputID++)
      ranges[inputID] = upperBounds_[inputID] - lowerBounds_[inputID];
   useArray = new int[newNSamples];
   for (sampleID = 0; sampleID < newNSamples; sampleID++)
      useArray[sampleID] = 0;
   iCntMax = new int[nInputs_];
   iCounts = new int[nInputs_];
   iSymbol = new int[nInputs_];
   iCntMax[0] = 1;
   for (inputID = 1; inputID < nInputs_; inputID++)
      iCntMax[inputID] = iCntMax[inputID-1] * symTable_[inputID-1];
   for (inputID = 0; inputID < nInputs_; inputID++)
      iCounts[inputID] = iSymbol[inputID] = 0;

   for (sampleID = 0; sampleID < nSamples_; sampleID++)
   {
      index = 0;
      for (inputID = 0; inputID < nInputs_; inputID++) 
      {
         if (randomize_ != 0)
            stepSize = ((double) (symTable_[inputID])) / ranges[inputID];
         else
            stepSize = ((double) (symTable_[inputID]-1)) / ranges[inputID];
         ddata = sampleMatrix_[sampleID][inputID] - lowerBounds_[inputID];
         curSymbol = (int) (ddata * stepSize * 1.00000000001);
         index += curSymbol * iCntMax[inputID];
      } 
      useArray[index] = sampleID + 1;
      //printf("Fact: use %d = %d\n", index, useArray[index]);
   }

   newSampleMatrix = new double*[newNSamples];
   for (sampleID = 0; sampleID < newNSamples; sampleID++)
      newSampleMatrix[sampleID] = new double[nInputs_];
   newSampleOutput = new double[newNSamples*nOutputs_];
   for (sampleID = 0; sampleID < newNSamples*nOutputs_; sampleID++)
      newSampleOutput[sampleID] = PSUADE_UNDEFINED;
   newSampleStates = new int[newNSamples];
   for (sampleID = 0; sampleID < newNSamples; sampleID++)
      newSampleStates[sampleID] = 0;
   
   for (sampleID = 0; sampleID < newNSamples; sampleID++)
   {
      index = -1;
      if (useArray[sampleID] > 0) index = useArray[sampleID] - 1;

      if (index >= 0)
      {
         for (inputID = 0; inputID < nInputs_; inputID++) 
            newSampleMatrix[sampleID][inputID] = 
                           sampleMatrix_[index][inputID]; 
         for (outputID = 0; outputID < nOutputs_; outputID++) 
            newSampleOutput[sampleID*nOutputs_+outputID] = 
                           sampleOutput_[index*nOutputs_+outputID]; 
         newSampleStates[sampleID] = sampleStates_[index];
      }
      for (inputID = 0; inputID < nInputs_; inputID++) 
      {
         curSymbol = iSymbol[inputID];
         iCounts[inputID]++;
         if (iCounts[inputID] == iCntMax[inputID])
         {
            iCounts[inputID] = 0;
            iSymbol[inputID]++;
            if (iSymbol[inputID] >= symTable_[inputID]) 
               iSymbol[inputID] = 0;
         }
         if (index < 0)
         {
            if (randomize_ != 0)
            {
               stepSize = ranges[inputID] / ((double) (symTable_[inputID]));
               newSampleMatrix[sampleID][inputID] = lowerBounds_[inputID] +
                                    (curSymbol + PSUADE_drand()) * stepSize;
            }
            else
            {
               stepSize = ranges[inputID] / ((double) (symTable_[inputID]-1));
               newSampleMatrix[sampleID][inputID] = 
                                curSymbol * stepSize + lowerBounds_[inputID];
            }
         }
      }
   }
   deleteSampleData();
   nSamples_ = newNSamples;
   sampleMatrix_ = newSampleMatrix;
   sampleOutput_ = newSampleOutput;
   sampleStates_ = newSampleStates;

   if (printLevel_ > 4)
   {
      printf("FactorialSampling::refine: nSamples = %d\n", nSamples_);
      printf("FactorialSampling::refine: nInputs  = %d\n", nInputs_);
      printf("FactorialSampling::refine: nOutputs = %d\n", nOutputs_);
      if (randomize_ != 0)
           printf("FactorialSampling::refine: randomize on\n");
      else printf("FactorialSampling::refine: randomize off\n");
      for (inputID = 0; inputID < nInputs_; inputID++)
         printf("    FactorialSampling input %3d = [%e %e]\n", inputID+1,
                lowerBounds_[inputID], upperBounds_[inputID]);
      if (symTable_ != NULL)
         for (inputID = 0; inputID < nInputs_; inputID++) 
            printf("    FactorialSampling symbol table %2d = %d\n",
                   inputID+1, symTable_[inputID]);
   }

   delete [] useArray;
   delete [] ranges;
   delete [] iCounts;
   delete [] iCntMax;
   delete [] iSymbol;
   return 0;
}

// ************************************************************************
// set input settings
// ------------------------------------------------------------------------
int FactorialSampling::setInputParams(int nInputs, int *counts, 
                                      double **settings, int *symtable)
{
   int inputID, sID;

   if (nInputs_ != 0 && nInputs != nInputs_)
   {
      printf("FactorialSampling::setInputParams - nInputs mismatch.\n");
      exit(1);
   }
   nInputs_ = nInputs;
   if (symtable != NULL)
   {
      symTable_ = new int[nInputs];
      for (inputID = 0; inputID < nInputs_; inputID++)
         symTable_[inputID] = symtable[inputID];
   }
   if (counts != NULL)
   {
      inputSettings_ = new double*[nInputs_];
      inputNumSettings_ = new int[nInputs_];
      for (inputID = 0; inputID < nInputs_; inputID++)
      {
         inputSettings_[inputID] = new double[counts[inputID]];
         inputNumSettings_[inputID] = counts[inputID];
         for (sID = 0; sID < counts[inputID]; sID++)
            inputSettings_[inputID][sID] = settings[inputID][sID];
      }
   }
   return 0;
}

