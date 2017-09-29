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
// Functions for the OA-based Latin Hypercube sampling class
// AUTHOR : CHARLES TONG
// DATE   : 2004
// ************************************************************************
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "OALHSampling.h"

// ************************************************************************
// link to external OA generator
// ************************************************************************
extern "C" 
{
   int  bose_link(int n, int ninputs, int str, int ***AA);
   void OA_strength(int q,int nrow,int ncol,int** A,int *str,int verbose);
}

// ************************************************************************
// Constructor
// ------------------------------------------------------------------------
OALHSampling::OALHSampling() : Sampling()
{
   samplingID_ = PSUADE_SAMP_OALH;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
OALHSampling::~OALHSampling()
{
}

// ************************************************************************
// initialize sampling data
// ------------------------------------------------------------------------
int OALHSampling::initialize(int initLevel)
{
   int    inputID, sampleID, *iArray1, strength, offset, **OAMatrix, nReps;
   int    repID, status, **tempMatrix, **iMatrix, *iArray2, nSamples1;
   int    *permute, symbolID, symbolID2, index, nSubSamples, nSamples2;
   double *ranges, scale, ddata, **perturbMatrix=NULL;

   if (nInputs_ == 0 || lowerBounds_ == NULL || upperBounds_ == NULL)
   {
      printf("OALHSampling ERROR: input not set up.\n");
      exit(1);
   }

   deleteSampleData();
   nReps = nReplications_;
   if (nSamples_ / nReps * nReps != nSamples_)
   {
      printf("OALHSampling ERROR: nSamples not multiples of replications.\n");
      exit(1);
   }
   ddata     = (double) (nSamples_ / nReps);
   ddata     = pow(ddata, 0.5000001);
   nSymbols_ = (int) ddata;
   nSamples1  = nSymbols_ * nSymbols_ * nReps;
   if (nSamples1 < nSamples_)
   {
      nSamples2 = (nSymbols_ + 1) * (nSymbols_ + 1) * nReps;
      if ((nSamples_ - nSamples1) < (nSamples2 - nSamples_))
         nSamples_ = nSamples1;
      else
      {
         nSamples_ = nSamples2;
         nSymbols_++;
      }
   }
   if (initLevel != 0) return 0;
   allocSampleData();

   if (printLevel_ > 4)
   {
      printf("OALHSampling initialize: nSamples = %d\n", nSamples_);
      printf("OALHSampling initialize: nInputs  = %d\n", nInputs_);
      printf("OALHSampling initialize: nOutputs = %d\n", nOutputs_);
      if (randomize_ != 0)
           printf("OALHSampling initialize: randomize on\n");
      else printf("OALHSampling initialize: randomize off\n");
      for (inputID = 0; inputID < nInputs_; inputID++)
         printf("    OALHSampling input %3d = [%e %e]\n", inputID+1,
                lowerBounds_[inputID], upperBounds_[inputID]);
      if (inputNumSettings_ != NULL)
      {
         for (inputID = 0; inputID < nInputs_; inputID++)
         {
            if (inputNumSettings_[inputID] != nSymbols_)
            {
               printf("OALHSampling ERROR: inputSetting not compabible.\n");
               exit(1);
            }
         }
      }
   }

   OAMatrix = new int*[nSamples_];
   for (sampleID = 0; sampleID < nSamples_; sampleID++)
      OAMatrix[sampleID] = new int[nInputs_];

   nReps       = nSamples_ / (nSymbols_ * nSymbols_);
   nSubSamples = nSamples_ / nReps;
   iArray1     = new int[nSymbols_];
   strength    = 2;
   offset      = 0;

   for (repID = 0; repID < nReps; repID++)
   {
      status = bose_link(nSubSamples, nInputs_, strength, &tempMatrix);
      if ((status >= 0 && status != nSubSamples) || status < 0)
      {
         printf("OALHSampling ERROR: Bose failure.\n");
         printf("    ==> Consult PSUADE developers for help.\n");
         exit(1);
      }
      for (sampleID = 0; sampleID < nSubSamples; sampleID++) 
      {
         for (inputID = 0; inputID < nInputs_; inputID++) 
            OAMatrix[offset+sampleID][inputID] = 
                               tempMatrix[sampleID][inputID];
         free(tempMatrix[sampleID]);
      }
      free(tempMatrix);

      for (inputID = 0; inputID < nInputs_; inputID++) 
      {
         generateRandomIvector(nSymbols_, iArray1);
         for (sampleID = 0; sampleID < nSubSamples; sampleID++) 
            OAMatrix[offset+sampleID][inputID] = 
                     iArray1[OAMatrix[offset+sampleID][inputID]];
      }

      OA_strength(nSymbols_,nSubSamples, nInputs_, &(OAMatrix[offset]),
                  &strength, 0);
      if (strength != 2)
      {
         printf("OALHSampling ERROR: failed strength 2 test.\n");
         printf("   ==> Please consult PSUADE developers.\n");
         exit(1);
      }
      offset += nSubSamples;
   }
   delete [] iArray1;

   iMatrix = new int*[nSymbols_];
   for (symbolID = 0; symbolID < nSymbols_; symbolID++) 
      iMatrix[symbolID] = new int[nSymbols_];
   iArray1 = new int[nSymbols_];
   iArray2 = new int[nSymbols_];
   permute = new int[nSubSamples];
   offset  = 0;

   for (repID = 0; repID < nReps; repID++)
   {
      for (inputID = 0; inputID < nInputs_; inputID++)
      {
         for (symbolID = 0; symbolID < nSymbols_; symbolID++)
         {
            for (symbolID2 = 0; symbolID2 < nSymbols_; symbolID2++)
                iArray1[symbolID2] = symbolID * nSymbols_ + symbolID2;
            generateRandomIvector(nSymbols_, iArray2);
            for (symbolID2 = 0; symbolID2 < nSymbols_; symbolID2++ )
               iMatrix[symbolID][symbolID2] = iArray1[iArray2[symbolID2]];
         }
         for (symbolID = 0; symbolID < nSymbols_; symbolID++) 
            iArray1[symbolID] = 0;
         for (sampleID = 0; sampleID < nSubSamples; sampleID++)
         {
            index = OAMatrix[offset+sampleID][inputID];
            permute[sampleID] = iMatrix[index][iArray1[index]++];
         }
         for (sampleID = 0; sampleID < nSubSamples; sampleID++) 
            OAMatrix[offset+sampleID][inputID] = permute[sampleID];
      }
      offset += nSubSamples;
   }
   for (symbolID = 0; symbolID < nSymbols_; symbolID++) 
      delete [] iMatrix[symbolID];
   delete [] iMatrix;
   delete [] iArray1;
   delete [] iArray2;
   delete [] permute;


   offset = 0;
   for (repID = 0; repID < nReps; repID++)
   {
      OA_strength(nSubSamples,nSubSamples,nInputs_,
                  &(OAMatrix[nSubSamples*repID]),&strength,0);
      if (strength != 1)
      {
         printf("OALHSampling ERROR: failed strength 1 test.\n");
         printf("   ==> Please consult PSUADE developers.\n");
         exit(1);
      }
      offset += nSubSamples;
   }

   iMatrix = new int*[nSamples_];
   for (sampleID = 0; sampleID < nSamples_; sampleID++) 
      iMatrix[sampleID] = new int[nInputs_];
   offset = 0;

   for (repID = 0; repID < nReps; repID++)
   {
      for (sampleID = 0; sampleID < nSubSamples; sampleID++)
         for (inputID = 0; inputID < nInputs_; inputID++)
            iMatrix[sampleID][inputID] = 
                    OAMatrix[offset+sampleID][inputID] / nSymbols_;
      OA_strength(nSymbols_, nSubSamples, nInputs_, iMatrix, &strength, 0);
      if (strength != 2)
      {
         printf("OALHSampling ERROR: failed strength 2 test.\n");
         printf("   ==> Please consult PSUADE developers.\n");
         exit(1);
      }
      offset += nSubSamples;
   }
   for (sampleID = 0; sampleID < nSamples_; sampleID++) 
      delete [] iMatrix[sampleID];
   delete [] iMatrix;

   ranges = new double[nInputs_];
   for (inputID = 0; inputID < nInputs_; inputID++) 
      ranges[inputID] = upperBounds_[inputID] - lowerBounds_[inputID];
   if (randomize_ != 0)
   {
      perturbMatrix = new double*[nInputs_];
      for (inputID = 0; inputID < nInputs_; inputID++)
      {
         perturbMatrix[inputID] = new double[nSymbols_*nSymbols_];
         for (symbolID = 0; symbolID < nSymbols_*nSymbols_; symbolID++)
            perturbMatrix[inputID][symbolID] = PSUADE_drand();
      }
      scale = 1.0 / (double) nSubSamples;
      for (sampleID = 0; sampleID < nSamples_; sampleID++)
      {
         for (inputID = 0; inputID < nInputs_; inputID++)
         {
            index = OAMatrix[sampleID][inputID];
            ddata = (perturbMatrix[inputID][index] + index) * scale;
            sampleMatrix_[sampleID][inputID] = 
                       ddata * ranges[inputID] + lowerBounds_[inputID];
         }
      }
   }
   else
   {
      scale = 1.0 / (double) (nSubSamples - 1);
      for (sampleID = 0; sampleID < nSamples_; sampleID++)
      {
         for (inputID = 0; inputID < nInputs_; inputID++)
         {
            if (inputSettings_ != NULL && inputSettings_[inputID] != NULL)
            {
               index = OAMatrix[sampleID][inputID];
               sampleMatrix_[sampleID][inputID] = 
                        inputSettings_[inputID][index];
            }
            else
            {
               ddata = (double) OAMatrix[sampleID][inputID];
               ddata = ddata * scale;
               sampleMatrix_[sampleID][inputID] = 
                       ddata * ranges[inputID] + lowerBounds_[inputID];
            }
         }
      }

   }
   for (sampleID = 0; sampleID < nSamples_; sampleID++)
      delete [] OAMatrix[sampleID];
   delete [] OAMatrix;
   if (perturbMatrix != NULL)
   {
      for (inputID = 0; inputID < nInputs_; inputID++)
         delete [] perturbMatrix[inputID];
      delete [] perturbMatrix;
   }
   delete [] ranges;
   return 0;
}

// ************************************************************************
// perform refinements
// ------------------------------------------------------------------------
int OALHSampling::refine(int refineRatio, int randomize, double thresh,
                         int nSamples, double *sampleErrors)
{
   int    nReps, symMult, *factors, ncount, ind1, samMult, strength;
   int    status, **tempMatrix, newNSymbols, newNSamples, *newStates;
   int    sampleID, repID, inputID, symbolID, nSubSamples, **newOAMatrix;
   int    nFactors, sampleOffset, newSampleOffset, ind2, sampleID2;
   int    newNSubSamples, **OAMatrix, *symbolArray, *iArray1, *iArray2;
   int    *iArray3, **iMatrix, symbolID2, outputID, nLevels;
   double *ranges, **newSamples, *newOutputs, **bounds, **perturbMatrix;
   double scale, dtemp;

   (void) randomize;
   (void) thresh;
   (void) nSamples;
   (void) sampleErrors;

   printf("OALHSampling refine - not implemented yet.\n");
   return -1;
   nLevels = refineRatio;
   nReps = nSamples_ / (nSymbols_ * nSymbols_);
   nSubSamples = nSamples_ / nReps;
   symMult = nLevels;
   factors  = factorize(nSamples_/nReps);
   nFactors = factors[0];
   ncount = 0;
   for (ind1 = 2; ind1 < nFactors; ind1++)
      if (factors[ind1] != factors[ind1-1]) ncount++;
   if (ncount > 1)
   {
      printf("OALHSampling refine ERROR: %d not prime power.\n",
             nSamples_/nReps);
      exit(1);
   }
   nLevels = nInputs_ - 1;
   printf("OALHSampling refine INFO: nLevels set to %d.\n",nLevels);
   symMult = factors[1] * nLevels;
   delete [] factors;

   samMult = symMult * symMult;
   strength = 2;
   status = bose_link(samMult, nInputs_, strength, &tempMatrix);
   if (status < 0)
   {
      printf("OASampling refine ERROR: cannot refine.\n");
      printf("         ==> Consult PSUADE developers for help.\n");
      exit(1);
   }
   for (sampleID = 0; sampleID < samMult; sampleID++)
      if (tempMatrix[sampleID] != NULL) free(tempMatrix[sampleID]);
   free(tempMatrix);

   randomize_ = 1;
   newNSymbols    = nSymbols_ * symMult;
   newNSamples    = newNSymbols * newNSymbols * nReps;
   newNSubSamples = newNSamples / nReps;

   ranges = new double[nInputs_];
   for (inputID = 0; inputID < nInputs_; inputID++)
      ranges[inputID] = upperBounds_[inputID] - lowerBounds_[inputID];
   bounds = new double*[nInputs_];
   for (inputID = 0; inputID < nInputs_; inputID++)
   {
      bounds[inputID] = new double[newNSubSamples+1];
      for (symbolID = 0; symbolID <= newNSubSamples; symbolID++)
         bounds[inputID][symbolID] = ranges[inputID] / newNSubSamples *
                                     symbolID + lowerBounds_[inputID];
   }
   perturbMatrix = new double*[nInputs_];
   for (inputID = 0; inputID < nInputs_; inputID++)
   {
      perturbMatrix[inputID] = new double[newNSubSamples];
      for (symbolID = 0; symbolID < newNSubSamples; symbolID++)
         perturbMatrix[inputID][symbolID] = PSUADE_drand();
   }
   for (sampleID = 0; sampleID < nSubSamples; sampleID++)
   {
      for (inputID = 0; inputID < nInputs_; inputID++)
      {
         dtemp = sampleMatrix_[sampleID][inputID];
         for (symbolID = 1; symbolID <= newNSubSamples; symbolID++)
            if (dtemp < bounds[inputID][symbolID]) break;
         symbolID--;
         if (symbolID >= newNSubSamples)
         {
            printf("OALHSampling refine ERROR (3)\n");
            exit(1);
         }
         dtemp -= lowerBounds_[inputID];
         dtemp  = dtemp / ranges[inputID] * newNSubSamples;
         dtemp -= (double) symbolID;
         perturbMatrix[inputID][symbolID] = dtemp;
      }
   }
   for (inputID = 0; inputID < nInputs_; inputID++) delete [] bounds[inputID];
   delete [] bounds;

   bounds = new double*[nInputs_];
   for (inputID = 0; inputID < nInputs_; inputID++)
   {
      bounds[inputID] = new double[newNSymbols+1];
      for (symbolID = 0; symbolID <= newNSymbols; symbolID++)
         bounds[inputID][symbolID] = ranges[inputID] / newNSymbols *
                                     symbolID + lowerBounds_[inputID];
   }
   OAMatrix = new int*[nSamples_];
   for (sampleID = 0; sampleID < nSamples_; sampleID++)
   {
      OAMatrix[sampleID] = new int[nInputs_];
      for (inputID = 0; inputID < nInputs_; inputID++)
      {
         dtemp = sampleMatrix_[sampleID][inputID];
         for (symbolID = 1; symbolID <= newNSymbols; symbolID++)
            if (dtemp < bounds[inputID][symbolID]) break;
         OAMatrix[sampleID][inputID] = symbolID - 1;
      }
   }
   for (inputID = 0; inputID < nInputs_; inputID++) delete [] bounds[inputID];
   delete [] bounds;
   delete [] ranges;

   newOAMatrix = new int*[newNSamples];
   for (sampleID = 0; sampleID < newNSamples; sampleID++)
      newOAMatrix[sampleID] = new int[nInputs_];

   sampleOffset = 0;
   newSampleOffset = 0;

   for (repID = 0; repID < nReps; repID++)
   {
      for (sampleID = 0; sampleID < nSubSamples; sampleID++)
         for (inputID = 0; inputID < nInputs_; inputID++)
            newOAMatrix[newSampleOffset+sampleID][inputID] =
               OAMatrix[sampleOffset+sampleID][inputID];
      newSampleOffset += nSubSamples;

      for (sampleID = 0; sampleID < nSubSamples; sampleID++)
      {
         status = bose_link(samMult, nInputs_, strength, &tempMatrix);

         for (inputID = 0; inputID < nInputs_; inputID++)
         {
            ind1 = tempMatrix[0][inputID];
            ind2 = OAMatrix[sampleOffset+sampleID][inputID] % symMult;
            if (ind1 != ind2)
            {
               for (sampleID2 = 1; sampleID2 < samMult; sampleID2++)
               {
                  ind1 = tempMatrix[sampleID2][inputID];
                  if (ind1 == ind2) break;
               }
            }
            if (ind1 == ind2)
            {
               ind1 = tempMatrix[0][inputID];
               for (sampleID2 = 0; sampleID2 < samMult; sampleID2++)
               {
                  if (tempMatrix[sampleID2][inputID] == ind1)
                     tempMatrix[sampleID2][inputID] = ind2;
                  else if (tempMatrix[sampleID2][inputID] == ind2)
                     tempMatrix[sampleID2][inputID] = ind1;
               }
            }
         }

         for (sampleID2 = 1; sampleID2 < samMult; sampleID2++)
         {
            for (inputID = 0; inputID < nInputs_; inputID++)
            {
               ind1 = OAMatrix[sampleOffset+sampleID][inputID] / symMult; 
               newOAMatrix[newSampleOffset+sampleID2-1][inputID] = 
                  tempMatrix[sampleID2][inputID] + ind1 * symMult;
            }
         }
         newSampleOffset += (samMult - 1);

         for (sampleID2 = 0; sampleID2 < samMult; sampleID2++)
            if (tempMatrix[sampleID2] != NULL) free(tempMatrix[sampleID2]);
         free(tempMatrix);
      }
      sampleOffset += nSubSamples;
   }

   newSampleOffset = 0;

   for (repID = 0; repID < nReps; repID++)
   {
      for (inputID = 0; inputID < nInputs_; inputID++)
      {
         for (symbolID = 0; symbolID < newNSymbols; symbolID++)
            symbolArray[symbolID] = symbolID;
         for (sampleID = 0; sampleID < nSubSamples; sampleID++)
         {
            ind1 = newOAMatrix[newSampleOffset+sampleID][inputID];
            symbolArray[ind1] = -1;
         }
         ncount = 0;
         for (symbolID = 0; symbolID < newNSymbols; symbolID++)
            if (symbolArray[symbolID] >= 0) ncount++;
         if (ncount > 1)
         {
            ncount = 0;
            for (symbolID = 0; symbolID < newNSymbols; symbolID++)
               if (symbolArray[symbolID] >= 0)
                  iArray3[ncount++] = symbolArray[symbolID];
            generateRandomIvector(ncount, iArray1);
            for (symbolID = 0; symbolID < ncount; symbolID++)
               iArray2[symbolID] = iArray3[iArray1[symbolID]];
            ncount = 0;
            for (symbolID = 0; symbolID < newNSymbols; symbolID++)
            {
               if (symbolArray[symbolID] >= 0)
                  symbolArray[symbolID] = iArray2[ncount++];
               else
                  symbolArray[symbolID] = symbolID;
            }
            for (sampleID = 0; sampleID < nSubSamples; sampleID++)
            {
               ind1 = newOAMatrix[newSampleOffset+sampleID][inputID];
               newOAMatrix[newSampleOffset+sampleID][inputID] = 
                                                   symbolArray[ind1];
            }
         }
      }

      OA_strength(newNSymbols, newNSubSamples, nInputs_,
                  &(newOAMatrix[newSampleOffset]), &strength, 0);
#if 0
      for (sampleID = 0; sampleID < newNSubSamples; sampleID++)
      {
         printf("sample %3d (%3d) : ",sampleID, currOffset);
         for (inputID = 0; inputID < nInputs_; inputID++)
            printf(" %3d ", OASamples[newSampleOffset+sampleID][inputID]);
         printf("\n");
      }
#endif
      if (strength != 2)
      {
         printf("OALHSampling refine ERROR: strength %d != 2\n",strength);
         exit(1);
      }
      newSampleOffset += newNSubSamples;
   }

   iMatrix = new int*[newNSymbols];
   for (symbolID = 0; symbolID < newNSymbols; symbolID++) 
      iMatrix[symbolID] = new int[newNSymbols];
   iArray1 = new int[newNSymbols];
   iArray2 = new int[newNSymbols];
   ranges  = new double[nInputs_];
   for (inputID = 0; inputID < nInputs_; inputID++)
      ranges[inputID] = upperBounds_[inputID] - lowerBounds_[inputID];
   bounds = new double*[nInputs_];
   for (inputID = 0; inputID < nInputs_; inputID++)
   {
      bounds[inputID] = new double[newNSubSamples+1];
      for (symbolID = 0; symbolID <= newNSubSamples; symbolID++)
         bounds[inputID][symbolID] = ranges[inputID] / newNSubSamples *
                                     symbolID + lowerBounds_[inputID];
   }

   sampleOffset = 0;
   newSampleOffset = 0;
   for (repID = 0; repID < nReps; repID++)
   {
      for (inputID = 0; inputID < nInputs_; inputID++)
      {
         for (symbolID = 0; symbolID < newNSymbols; symbolID++)
         {
            for (symbolID2 = 0; symbolID2 < newNSymbols; symbolID2++ )
               iMatrix[symbolID][symbolID2] = symbolID * newNSymbols +
                                              symbolID2;
         }
         for (sampleID = 0; sampleID < nSubSamples; sampleID++)
         {
            dtemp = sampleMatrix_[sampleOffset+sampleID][inputID];
            for (symbolID = 1; symbolID <= newNSubSamples; symbolID++)
               if (dtemp < bounds[inputID][symbolID]) break;
            symbolID--;
            iMatrix[symbolID/newNSubSamples][symbolID%newNSubSamples] =
               - iMatrix[symbolID/newNSubSamples][symbolID%newNSubSamples] - 1;
         }
         for (symbolID = 0; symbolID < newNSymbols; symbolID++)
         {
            ncount = 0;
            for (symbolID2 = 0; symbolID2 < newNSymbols; symbolID2++)
                if (iMatrix[symbolID][symbolID2] >= 0) 
                   iMatrix[symbolID][ncount++] = iMatrix[symbolID][symbolID2];
            if (ncount > 1)
            {
               for (symbolID2 = 0; symbolID2 < ncount; symbolID2++)
                     iArray1[symbolID2] = iMatrix[symbolID][symbolID2];
               generateRandomIvector(ncount, iArray2);
               for (symbolID2 = 0; symbolID2 < ncount; symbolID2++ )
                  iMatrix[symbolID][symbolID2] = iArray1[iArray2[symbolID2]];
            }
         }
         for (symbolID = 0; symbolID < newNSymbols; symbolID++) 
            iArray1[symbolID] = 0;
         for (sampleID = nSubSamples; sampleID < newNSubSamples; sampleID++) 
         {
            ind1 = newOAMatrix[newSampleOffset+sampleID][inputID];
            ind2  = iMatrix[ind1][iArray1[ind1]++];
            if (ind2 < 0 || ind2 >= newNSubSamples)
            {
               printf("OALHSampling::refine ERROR : ind2 problem.\n");
               exit(1);
            }
            newOAMatrix[newSampleOffset+sampleID][inputID] = ind2;
         }
      }
      newSampleOffset += newNSubSamples;
      sampleOffset += nSubSamples;
   }
   for (symbolID = 0; symbolID < nSymbols_; symbolID++) 
      delete [] iMatrix[symbolID];
   delete [] iMatrix;
   delete [] iArray1;
   delete [] iArray2;
   for (inputID = 0; inputID < nInputs_; inputID++) delete [] bounds[inputID];
   delete [] bounds;
   delete [] ranges;
   for (sampleID = 0; sampleID < nSamples_; sampleID++)
      delete [] OAMatrix[sampleID];
   delete [] OAMatrix;

   newSamples = new double*[newNSamples];
   newOutputs = new double[newNSamples*nOutputs_];
   newStates  = new int[newNSamples];
   for (sampleID = 0; sampleID < newNSamples; sampleID++)
   {
      newSamples[sampleID] = new double[nInputs_];
      newStates[sampleID] = 0;
      for (outputID = 0; outputID < nOutputs_; outputID++)
         newOutputs[sampleID*nOutputs_+outputID] = PSUADE_UNDEFINED;
   }

   ranges = new double[nInputs_];
   for (inputID = 0; inputID < nInputs_; inputID++) 
      ranges[inputID] = upperBounds_[inputID] - lowerBounds_[inputID];
   scale = 1.0 / (double) newNSubSamples;
   newSampleOffset = 0;
   sampleOffset = 0;
   for (repID = 0; repID < nReps; repID++)
   {
      for (sampleID = 0; sampleID < nSubSamples; sampleID++)
      {
         for (inputID = 0; inputID < nInputs_; inputID++)
            newSamples[newSampleOffset+sampleID][inputID] =
                sampleMatrix_[sampleOffset+sampleID][inputID];
         for (outputID = 0; outputID < nOutputs_; outputID++)
            newOutputs[(newSampleOffset+sampleID)*nOutputs_+outputID] =
               sampleOutput_[(sampleOffset+sampleID)*nOutputs_+outputID];
         newStates[newSampleOffset+sampleID] =
               sampleStates_[sampleOffset+sampleID];
      }
      for (sampleID = nSubSamples; sampleID < newNSamples; sampleID++)
      {
         for (inputID = 0; inputID < nInputs_; inputID++)
         {
            ind1  = newOAMatrix[newSampleOffset+sampleID][inputID];
            dtemp = (perturbMatrix[inputID][ind1 ] + ind1 ) * scale;
            newSamples[newSampleOffset+sampleID][inputID] = 
                       dtemp * ranges[inputID] + lowerBounds_[inputID];
         }
      }
      newSampleOffset += newNSubSamples; 
      sampleOffset += nSubSamples; 
   }
   delete [] ranges;
   for (sampleID = 0; sampleID < newNSamples; sampleID++)
      delete [] newOAMatrix;
   delete [] newOAMatrix;
   if (perturbMatrix != NULL)
   {
      for (inputID = 0; inputID < nInputs_; inputID++)
         delete [] perturbMatrix[inputID];
      delete [] perturbMatrix;
   }

   deleteSampleData();
   nSamples_ = newNSamples;
   nSymbols_ = newNSymbols;
   sampleMatrix_ = newSamples;
   sampleOutput_ = newOutputs;
   sampleStates_ = newStates;

   if (printLevel_ > 4)
   {
      printf("OALHSampling refine: nSamples = %d\n", nSamples_);
      printf("OALHSampling refine: nInputs  = %d\n", nInputs_);
      printf("OALHSampling refine: nOutputs = %d\n", nOutputs_);
      if (randomize_ != 0)
           printf("OALHSampling refine: randomize on\n");
      else printf("OALHSampling refine: randomize off\n");
      for (inputID = 0; inputID < nInputs_; inputID++)
         printf("    OALHSampling input %3d = [%e %e]\n", inputID+1,
                lowerBounds_[inputID], upperBounds_[inputID]);
      if (inputNumSettings_ != NULL || symTable_ != NULL)
         printf("OALHSampling refine: diable input settings, symbol table.\n");
   }
   return 0;
}

// ************************************************************************
// set input settings
// ------------------------------------------------------------------------
int OALHSampling::setInputParams(int nInputs, int *counts,
                                 double **settings, int *symtable)
{
   int inputID, inputCnt, sID;

   if (nInputs_ != 0 && nInputs != nInputs_)
   {
      printf("OALHSampling setInputParams ERROR: nInputs mismatch.\n");
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
      inputCnt = 0;
      for (inputID = 0; inputID < nInputs_; inputID++)
      {
         if (counts[inputID] != 0 && counts[inputID] != nSymbols_*nSymbols_)
         {
            printf("OALHSampling setInputParams ERROR: counts mismatch %d %d.\n",
                   counts[inputID], nSymbols_*nSymbols_);
            exit(1);
         }
         else if (counts[inputID] == nSymbols_*nSymbols_) inputCnt++;
      }
      if (inputCnt > 0)
      {
         inputSettings_ = new double*[nInputs_];
         inputNumSettings_ = new int[nInputs_];
         for (inputID = 0; inputID < nInputs_; inputID++)
         {
            if (counts[inputID] == nSymbols_*nSymbols_)
            {
               inputSettings_[inputID] = new double[counts[inputID]];
               inputNumSettings_[inputID] = counts[inputID];
               for (sID = 0; sID < counts[inputID]; sID++)
                  inputSettings_[inputID][sID] = settings[inputID][sID];
            }
            else
            {
               inputSettings_[inputID] = NULL;
               inputNumSettings_[inputID] = 0;
            }
         }
      }
   }
   return 0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
OALHSampling& OALHSampling::operator=(const OALHSampling &)
{
   printf("OALHSampling operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

