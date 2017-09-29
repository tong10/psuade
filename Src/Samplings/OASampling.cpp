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
// Functions for the orthogonal array sampling class
// AUTHOR : CHARLES TONG
// DATE   : 2003
// ************************************************************************

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "OASampling.h"

#define PABS(x) ((x) > 0 ? (x) : -(x))

// ************************************************************************
// external functions
// ************************************************************************

extern "C" 
{
   int  bose_link(int n, int ninputs, int str, int ***AA);
   void OA_strength(int q,int nrow,int ncol,int** A,int *str,int verbose);
}

// ************************************************************************
// Constructor
// ------------------------------------------------------------------------
OASampling::OASampling() : Sampling()
{
   samplingID_ = PSUADE_SAMP_OA;
   trueRandom_ = 0;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
OASampling::~OASampling()
{
}

// ************************************************************************
// initialize sampling data
// ------------------------------------------------------------------------
int OASampling::initialize(int initLevel)
{
   int    ii, jj, kk, status, *intVec, strength, nReps, maxMinDist, dist2;
   int    **permMatrix, **tempMatrix, index, repID, offset, nSamples1;
   int    nSamples2, ntimes=1, ll, minDist, dist, ss, ss2, **storeMatrix;
   double *ranges, scale, ddata, **perturbMatrix=NULL;

   if (nInputs_ == 0 || lowerBounds_ == NULL || upperBounds_ == NULL)
   {
      printf("OASampling::initialize ERROR - input not set up.\n");
      exit(1);
   }
   if (nSamples_ == 0)
   {
      printf("OASampling::initialize ERROR - nSamples = 0.\n");
      exit(1);
   }

   deleteSampleData();
   trueRandom_ = 0;
   if (randomize_ & 2) trueRandom_ = 1;
   if (randomize_ & 1) randomize_ = 1;
   nReps = nReplications_;
   if (nSamples_ / nReps * nReps != nSamples_)
   {
      printf("OASampling : nSamples must be multiples of replications.\n");
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
      printf("OASampling: initialize: nSamples = %d\n", nSamples_);
      printf("OASampling: initialize: nInputs  = %d\n", nInputs_);
      printf("OASampling: initialize: nOutputs = %d\n", nOutputs_);
      if (randomize_ != 0)
           printf("OASampling: initialize: randomize on\n");
      else printf("OASampling: initialize: randomize off\n");
      if (trueRandom_ != 0)
           printf("OASampling: initialize: more randomize on\n");
      else printf("OASampling: initialize: more randomize off\n");
      for (jj = 0; jj < nInputs_; jj++)
         printf("    OASampling input %3d = [%e %e]\n", jj+1,
                lowerBounds_[jj], upperBounds_[jj]);
      if (inputNumSettings_ != NULL)
      {
         for (jj = 0; jj < nInputs_; jj++)
         {
            if (inputNumSettings_[jj] != nSymbols_)
            {
               printf("OASampling ERROR: inputSetting not compabible.\n");
               exit(1);
            }
         }
      }
   }

   nReps = nSamples_ / (nSymbols_ * nSymbols_);
   permMatrix = new int*[nSamples_+nSamples_/nReps];
   for (ii = 0; ii < nSamples_+nSamples_/nReps; ii++) 
      permMatrix[ii] = new int[nInputs_];
   storeMatrix = new int*[nSamples_/nReps];
   for (ii = 0; ii < nSamples_/nReps; ii++) storeMatrix[ii] = new int[nInputs_];
   intVec = new int[nSymbols_];

   strength = 2;
   offset = 0;
   maxMinDist = 0;
   for (repID = 0; repID < nReps; repID++)
   {
      if (printLevel_ > 4)
         printf("OASampling: creating the %d-th (of %d) replication.\n",
                repID+1,nReps);

      status = bose_link(nSamples_/nReps, nInputs_, strength, &tempMatrix);
      if ((status >= 0 && status != (nSamples_/nReps)) || status < 0)
      {
         printf("OASampling ERROR: Bose failure.\n");
         printf("    ==> Consult PSUADE developers for help.\n");
         exit(1);
      }
      for (ii = 0; ii < nSamples_/nReps; ii++) 
      {
         for (jj = 0; jj < nInputs_; jj++) 
            permMatrix[nSamples_+ii][jj] = tempMatrix[ii][jj];
         free(tempMatrix[ii]);
      }
      free(tempMatrix);
     
      for (ll = 0; ll < ntimes; ll++)
      {

         for (jj = 0; jj < nInputs_; jj++) 
         {
            generateRandomIvector(nSymbols_, intVec);
            for (ii = 0; ii < nSamples_/nReps; ii++) 
               permMatrix[offset+ii][jj] = 
                        intVec[permMatrix[nSamples_+ii][jj]];
         }
         if (ntimes > 1)
         {
            minDist = 2 * nSymbols_ * nInputs_;
            for (ss = 0; ss < nSamples_/nReps; ss++)
            {
               for (ss2 = ss+1; ss2 < nSamples_/nReps; ss2++)
               {
                  dist = 0;
                  for (ii = 0; ii < nInputs_; ii++)
                  {
                     dist2 = permMatrix[offset+ss][ii]-permMatrix[offset+ss2][ii];
                     dist += PABS(dist2);
                  }
                  if (dist > 0 && dist < minDist) minDist = dist;
               }
            }
         }
         else minDist = maxMinDist + 1;
         if (minDist > maxMinDist)
         {
            for (ss = 0; ss < nSamples_/nReps; ss++)
               for (ii = 0; ii < nInputs_; ii++)
                  storeMatrix[ss][ii] = permMatrix[offset+ss][ii];
            maxMinDist = minDist;
         }
      }
      for (ss = 0; ss < nSamples_/nReps; ss++)
         for (ii = 0; ii < nInputs_; ii++)
            permMatrix[offset+ss][ii] = storeMatrix[ss][ii];

      if (nSamples_/nReps < 1000 && printLevel_ > 5)
      {
         OA_strength(nSymbols_,nSamples_/nReps, nInputs_, 
                     &(permMatrix[offset]), &strength, 0);
         if (strength != 2)
         {
            printf("OASampling ERROR: fail strength 2 test.\n");
            printf("   ==> Please consult PSUADE developers.\n");
            exit(1);
         }
      }
      offset += nSamples_/nReps;
   }
   delete [] intVec;
   for (ii = 0; ii < nSamples_/nReps; ii++)
      if (storeMatrix[ii] != NULL) delete [] storeMatrix[ii];
   delete [] storeMatrix;

#if 0
   for (ii = 0; ii < nSamples_; ii++) 
   {
      printf("OA sample %3d = ", ii);
      for (jj = 0; jj < nInputs_; jj++) 
         printf("%d ", permMatrix[ii][jj]);
      printf("\n");
   }
#endif

   // ----------------------------------------------------------------
   // generate sample data
   // ----------------------------------------------------------------
   ranges = new double[nInputs_];
   for (jj = 0; jj < nInputs_; jj++) 
      ranges[jj] = upperBounds_[jj] - lowerBounds_[jj];

   if (trueRandom_ != 0)
   {
      scale = 1.0 / (double) nSymbols_;
      for (ii = 0; ii < nSamples_; ii++)
      {
         for (jj = 0; jj < nInputs_; jj++)
         {
            index = permMatrix[ii][jj];
            ddata = (PSUADE_drand() + index) * scale;
            sampleMatrix_[ii][jj] = ddata * ranges[jj] + 
                                               lowerBounds_[jj];
         }
      }
   }

   else if (randomize_ != 0)
   {
      perturbMatrix = new double*[nInputs_];
      for (jj = 0; jj < nInputs_; jj++)
      {
         perturbMatrix[jj] = new double[nSymbols_];
         for (kk = 0; kk < nSymbols_; kk++)
            perturbMatrix[jj][kk] = PSUADE_drand();
      }
      scale = 1.0 / (double) nSymbols_;
      for (ii = 0; ii < nSamples_; ii++)
      {
         for (jj = 0; jj < nInputs_; jj++)
         {
            index = permMatrix[ii][jj];
            ddata = (perturbMatrix[jj][index] + index) * scale;
            sampleMatrix_[ii][jj] = ddata * ranges[jj] + lowerBounds_[jj];
         }
      }
   }

   else
   {
      scale = 1.0 / (double) (nSymbols_ - 1);
      for (ii = 0; ii < nSamples_; ii++)
      {
         for (jj = 0; jj < nInputs_; jj++)
         {
            if (inputSettings_ != NULL && inputSettings_[jj] != NULL)
            {
               index = permMatrix[ii][jj];
               sampleMatrix_[ii][jj] = inputSettings_[jj][index];
            }
            else
            {
               ddata = (double) permMatrix[ii][jj];
               ddata = ddata * scale;
               sampleMatrix_[ii][jj] = ddata * ranges[jj] + lowerBounds_[jj];
            }
         }
      }
   }

   for (ii = 0; ii < nSamples_+nSamples_/nReps; ii++)
      if (permMatrix[ii] != NULL) delete [] permMatrix[ii];
   delete [] permMatrix;
   delete [] ranges;
   if (perturbMatrix != NULL)
   {
      for (jj = 0; jj < nInputs_; jj++) delete [] perturbMatrix[jj];
      delete [] perturbMatrix;
   }
   return 0;
}

// ************************************************************************
// perform refinements
// ------------------------------------------------------------------------
int OASampling::refine(int refineRatio, int randomize, double threshold,
                       int nSamples, double *sampleErrors)
{
   int    nReps, symMult, *factors, nFactors, ncount, ind1, samMult;
   int    strength, status, **tempMatrix, newNSamples, newNSymbols;
   int    *newStates, ii, kk, jj, repID, outputID;
   int    **OASamples, *symbolArray, sampleOffset, newSampleOffset;
   int    ii2, *iArray1, *iArray2, *iArray3, currOffset, ind2;
   int    nLevels;
   double *ranges, **newSamples, *newOutputs, **perturbMatrix=NULL;
   double **bounds, scale, ddata;

   (void) refineRatio;
   (void) randomize;
   (void) threshold;
   (void) nSamples;
   (void) sampleErrors;

   nReps = nSamples_ / (nSymbols_ * nSymbols_);
   factors  = factorize(nSamples_/nReps);
   nFactors = factors[0];
   ncount = 0;
   for (ind1 = 2; ind1 < nFactors; ind1++)
      if (factors[ind1] != factors[ind1-1]) ncount++;
   if (ncount > 1)
   {
      printf("OASampling refine ERROR: %d not prime power.\n",
             nSamples_/nReps);
      exit(1);
   }
   if (nInputs_ > 2) nLevels = nInputs_ - 1;
   else              nLevels = 2;
   printf("OASampling refine INFO: nLevels set to %d.\n",nInputs_-1);
   symMult = nLevels;
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
   for (ii = 0; ii < samMult; ii++)
      if (tempMatrix[ii] != NULL) free(tempMatrix[ii]);
   free(tempMatrix);

   ranges = new double[nInputs_];
   for (jj = 0; jj < nInputs_; jj++)
      ranges[jj] = upperBounds_[jj] - lowerBounds_[jj];
   randomize_ = 1;


   newNSymbols = nSymbols_ * symMult;
   newNSamples = newNSymbols * newNSymbols * nReps;

   newSamples = new double*[newNSamples];
   newOutputs = new double[newNSamples*nOutputs_];
   newStates  = new int[newNSamples];
   for (ii = 0; ii < newNSamples; ii++)
   {
      newSamples[ii] = new double[nInputs_];
      newStates[ii] = 0;
      for (outputID = 0; outputID < nOutputs_; outputID++)
         newOutputs[ii*nOutputs_+outputID] = PSUADE_UNDEFINED;
   }

   bounds = new double*[nInputs_];
   for (jj = 0; jj < nInputs_; jj++)
   {
      bounds[jj] = new double[newNSymbols+1];
      for (kk = 0; kk <= newNSymbols; kk++)
         bounds[jj][kk] = ranges[jj] / newNSymbols * kk + lowerBounds_[jj];
   }

   perturbMatrix = new double*[nInputs_];
   for (jj = 0; jj < nInputs_; jj++)
   {
      perturbMatrix[jj] = new double[newNSymbols];
      for (kk = 0; kk < newNSymbols; kk++)
         perturbMatrix[jj][kk] = PSUADE_drand();
   }

   for (ii = 0; ii < nSamples_/nReps; ii++)
   {
      for (jj = 0; jj < nInputs_; jj++)
      {
         ddata = sampleMatrix_[ii][jj];
         if (ddata == bounds[jj][newNSymbols]) kk = newNSymbols;
         else
         {
            for (kk = 1; kk < (newNSymbols+1); kk++)
               if (ddata < bounds[jj][kk]) break;
         }
         kk--;
         if (kk >= newNSymbols)
         {
            printf("OASampling refine ERROR (3): %d %d %e\n",kk,
                   newNSymbols,ddata);
            exit(1);
         }
         ddata -= lowerBounds_[jj];
         ddata  = ddata / ranges[jj] * (double) newNSymbols;
         ddata -= (double) kk;
         perturbMatrix[jj][kk] = ddata;
      }
   }

   symbolArray = new int[newNSymbols];
   iArray1 = new int[newNSymbols];
   iArray2 = new int[newNSymbols];
   iArray3 = new int[newNSymbols];
   OASamples = new int*[newNSamples];
   for (ii = 0; ii < newNSamples; ii++)
      OASamples[ii] = new int[nInputs_];

   sampleOffset = 0;
   newSampleOffset = 0;
   for (repID = 0; repID < nReps; repID++)
   {
      for (ii = 0; ii < nSamples_/nReps; ii++)
      {
         for (jj = 0; jj < nInputs_; jj++)
            newSamples[newSampleOffset+ii][jj] =
                          sampleMatrix_[sampleOffset+ii][jj];
         for (outputID = 0; outputID < nOutputs_; outputID++)
            newOutputs[(newSampleOffset+ii)*nOutputs_+outputID] =
                  sampleOutput_[(sampleOffset+ii)*nOutputs_+outputID];
         newStates[newSampleOffset+ii] = 
            sampleStates_[sampleOffset+ii];
      }

      for (ii = 0; ii < nSamples_/nReps; ii++)
      {
         for (jj = 0; jj < nInputs_; jj++)
         {
            ddata = sampleMatrix_[sampleOffset+ii][jj];
            if (ddata == bounds[jj][newNSymbols]) kk = newNSymbols;
            else
            {
               for (kk = 1; kk < (newNSymbols+1); kk++)
                  if (ddata < bounds[jj][kk]) break;
            }
            OASamples[newSampleOffset+ii][jj] = kk - 1;
         }
      }
      currOffset = newSampleOffset;
      newSampleOffset += (nSamples_ / nReps);

      for (ii = 0; ii < nSamples_/nReps; ii++)
      {
         status = bose_link(samMult, nInputs_, strength, &tempMatrix);

         for (jj = 0; jj < nInputs_; jj++)
         {
            ind1 = tempMatrix[0][jj];
            ind2 = OASamples[currOffset+ii][jj] % symMult;
            if (ind1 != ind2)
            {
               for (ii2 = 1; ii2 < samMult; ii2++)
               {
                  ind1 = tempMatrix[ii2][jj];
                  if (ind1 == ind2) break;
               }
            }
            if (ind1 == ind2)
            {
               ind1 = tempMatrix[0][jj];
               for (ii2 = 0; ii2 < samMult; ii2++)
               {
                  if (tempMatrix[ii2][jj] == ind1)
                     tempMatrix[ii2][jj] = ind2;
                  else if (tempMatrix[ii2][jj] == ind2)
                     tempMatrix[ii2][jj] = ind1;
               }
            }
         }

         for (ii2 = 1; ii2 < samMult; ii2++)
         {
            for (jj = 0; jj < nInputs_; jj++)
            {
               ind1 = OASamples[currOffset+ii][jj] / symMult; 
               OASamples[newSampleOffset+ii2-1][jj] = 
                  tempMatrix[ii2][jj] + ind1 * symMult;
            }
         }
         newSampleOffset += (samMult - 1);

         for (ii2 = 0; ii2 < samMult; ii2++)
            if (tempMatrix[ii2] != NULL) free(tempMatrix[ii2]);
         free(tempMatrix);
      }

      for (jj = 0; jj < nInputs_; jj++)
      {
         for (kk = 0; kk < newNSymbols; kk++) symbolArray[kk] = kk;
         for (ii = 0; ii < nSamples_/nReps; ii++)
         {
            ind1 = OASamples[currOffset+ii][jj];
            symbolArray[ind1] = -1;
         }
         ncount = 0;
         for (kk = 0; kk < newNSymbols; kk++)
            if (symbolArray[kk] >= 0) ncount++;
         if (ncount > 1)
         {
            ncount = 0;
            for (kk = 0; kk < newNSymbols; kk++)
               if (symbolArray[kk] >= 0)
                  iArray3[ncount++] = symbolArray[kk];
            generateRandomIvector(ncount, iArray1);
            for (kk = 0; kk < ncount; kk++)
               iArray2[kk] = iArray3[iArray1[kk]];
            ncount = 0;
            for (kk = 0; kk < newNSymbols; kk++)
            {
               if (symbolArray[kk] >= 0)
                  symbolArray[kk] = iArray2[ncount++];
               else
                  symbolArray[kk] = kk;
            }
            for (ii = 0; ii < newNSamples/nReps; ii++)
            {
               ind1 = OASamples[currOffset+ii][jj];
               OASamples[currOffset+ii][jj] = symbolArray[ind1];
            }
         }
      }

      if (newNSamples/nReps < 1000)
      {
         OA_strength(newNSymbols, newNSamples/nReps, nInputs_,
                     &(OASamples[currOffset]), &strength, 0);
         if (strength != 2)
         {
            printf("OASampling::refine ERROR : strength %d != 2\n",
                   strength);
            exit(1);
         }
      }
#if 0
      for (ii = 0; ii < newNSamples/nReps; ii++)
      {
         printf("sample %3d (%3d) : ",ii, currOffset);
         for (jj = 0; jj < nInputs_; jj++)
            printf(" %3d ", OASamples[currOffset+ii][jj]);
         printf("\n");
      }
#endif

      currOffset += nSamples_ / nReps; 
      for (ii = currOffset; ii < newSampleOffset; ii++)
      {
         scale = 1.0 / ((double) newNSymbols);
         for (jj = 0; jj < nInputs_; jj++)
         {
            ind1 = OASamples[ii][jj];
            ddata = (perturbMatrix[jj][ind1] + ind1) * scale;
            newSamples[ii][jj] = ddata * ranges[jj] + lowerBounds_[jj];
         }
      }
      sampleOffset += (nSymbols_ * nSymbols_);
   }

   deleteSampleData();
   nSamples_ = newNSamples;
   nSymbols_ = newNSymbols;
   sampleMatrix_ = newSamples;
   sampleOutput_ = newOutputs;
   sampleStates_ = newStates;

   if (printLevel_ > 4)
   {
      printf("OASampling refine: nSamples = %d\n", nSamples_);
      printf("OASampling refine: nInputs  = %d\n", nInputs_);
      printf("OASampling refine: nOutputs = %d\n", nOutputs_);
      if (randomize_ != 0)
           printf("OASampling refine: randomize on\n");
      else printf("OASampling refine: randomize off\n");
      if (trueRandom_ != 0)
           printf("OASampling refine: more randomize on\n");
      else printf("OASampling refine: more randomize off\n");
      for (jj = 0; jj < nInputs_; jj++)
         printf("    OASampling input %3d = [%e %e]\n", jj+1,
                lowerBounds_[jj], upperBounds_[jj]);
   }

   delete [] ranges;
   for (jj = 0; jj < nInputs_; jj++)
       delete [] bounds[jj];
   delete [] bounds;
   for (ii = 0; ii < newNSamples; ii++) delete [] OASamples[ii];
   for (jj = 0; jj < nInputs_; jj++) delete [] perturbMatrix[jj];
   delete [] perturbMatrix;
   delete [] OASamples;
   delete [] symbolArray;
   delete [] iArray1;
   delete [] iArray2;
   delete [] iArray3;
   return 0;
}

// ************************************************************************
// set input settings
// ------------------------------------------------------------------------
int OASampling::setInputParams(int nInputs, int *counts, 
                               double **settings, int *symtable)
{
   int jj, inputCnt, sID;

   if (nInputs_ != 0 && nInputs != nInputs_)
   {
      printf("OASampling setInputParams ERROR: nInputs mismatch.\n");
      exit(1);
   }
   nInputs_ = nInputs;
   if (symtable != NULL)
   {
      symTable_ = new int[nInputs];
      for (jj = 0; jj < nInputs_; jj++) symTable_[jj] = symtable[jj];
   }
   if (counts != NULL)
   {
      inputCnt = 0;
      for (jj = 0; jj < nInputs_; jj++)
      {
         if (counts[jj] != 0 && counts[jj] != nSymbols_)
         {
            printf("OASampling setInputParams ERROR: counts mismatch.\n");
            exit(1);
         }
         else if (counts[jj] == nSymbols_) inputCnt++;
      }
      if (inputCnt > 0)
      {
         inputSettings_ = new double*[nInputs_];
         inputNumSettings_ = new int[nInputs_];
         for (jj = 0; jj < nInputs_; jj++)
         {
            if (counts[jj] == nSymbols_)
            {
               inputSettings_[jj] = new double[counts[jj]];
               inputNumSettings_[jj] = counts[jj];
               for (sID = 0; sID < counts[jj]; sID++)
                  inputSettings_[jj][sID] = settings[jj][sID];
            }
            else
            {
               inputSettings_[jj] = NULL;
               inputNumSettings_[jj] = 0;
            }
         }
      }
   }
   return 0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
OASampling& OASampling::operator=(const OASampling &)
{
   printf("OASampling operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

