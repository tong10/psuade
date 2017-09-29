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
// Functions for the Morris one-at-a-time class (improved) 
// AUTHOR : CHARLES TONG
// DATE   : 2004 (updated in 2006 and 2007)
//*------------------------------------------------------------------------
// ************************************************************************
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "pData.h"
#include "pData.h"
#include "PsuadeData.h"
#include "FuncApprox.h"
#include "Psuade.h"
#include "MOATSampling.h"
#include "PrintingTS.h"
#define PABS(x) ((x) > 0 ? (x) : -(x))

//*************************************************************************
//* Constructor
//*------------------------------------------------------------------------
MOATSampling::MOATSampling() : Sampling()
{
   samplingID_ = PSUADE_SAMP_MOAT;
   P_ = 4;
   inputSubset_ = NULL;
}

//*************************************************************************
//* Copy Constructor added by Bill Oliver
//*------------------------------------------------------------------------
MOATSampling::MOATSampling(const MOATSampling & ms) : Sampling()
{
   samplingID_ = ms.samplingID_;
   P_ = ms.P_;
   nInputs_ = ms.nInputs_;
   inputSubset_ = new int[nInputs_];
   for(int i = 0; i < nInputs_; i++)
     inputSubset_[i] = ms.inputSubset_[i];

   printLevel_ = ms.printLevel_;
   samplingID_ = ms.samplingID_;
   nSamples_ = ms.nSamples_;
   nOutputs_ = ms.nOutputs_;
   randomize_ = ms.randomize_;
   nReplications_ = ms.nReplications_;
   lowerBounds_ = new double[nInputs_];
   upperBounds_ = new double[nInputs_];
   for (int i = 0; i < nInputs_; i++)
   {
      lowerBounds_[i] = ms.lowerBounds_[i];
      upperBounds_[i] = ms.upperBounds_[i];
   }
   sampleMatrix_ = new double*[nSamples_];
   for (int i = 0; i < nSamples_; i++)
   {
      sampleMatrix_[i] = new double[nInputs_];
      for(int j = 0; j < nInputs_; j++)
         sampleMatrix_[i][j] = ms.sampleMatrix_[i][j];
   }
   sampleOutput_ = new double[nSamples_*nOutputs_];
   for (int i = 0; i < nSamples_*nOutputs_; i++)
      sampleOutput_[i] = ms.sampleOutput_[i];
   sampleStates_ = new int[nSamples_];
   for (int i = 0; i < nSamples_; i++)
      sampleStates_[i] = ms.sampleStates_[i];
}

//********************************************************************
//* destructor 
//*-------------------------------------------------------------------
MOATSampling::~MOATSampling()
{
   if (inputSubset_ != NULL) delete [] inputSubset_;
}

//*************************************************************************
//* initialize the sampling data
//*------------------------------------------------------------------------
int MOATSampling::initialize(int initLevel)
{
   int    ii, ii2, rr, ss, nReps, nn, randomize, *bins, currBin, nSub=0;
   int    kk1, kk2, maxReps=500, maxSamples, index, base1, base2, setFlag=0;
   double **BS, *ranges, ddata, *tempX, maxDist, dDist;
   char   *cString, partitionFile[200], winput1[200], winput2[200];
   FILE*  fp;

   if (nSamples_ == 0)
   {
      printOutTS(PL_ERROR,
           "MOATSampling::initialize ERROR - nSamples = 0.\n");
      exit(1);
   }
   if (nInputs_ == 0 || lowerBounds_ == NULL || upperBounds_ == NULL)
   {
      printOutTS(PL_ERROR,
           "MOATSampling::initialize ERROR - input not set up.\n");
      exit(1);
   }

   deleteSampleData();
   randomize = (randomize_ & 1);
   if (nSamples_/(nInputs_+1) * (nInputs_+1) != nSamples_) 
   {
      printOutTS(PL_INFO,
           "MOATSampling: nSamples should be multiples of nInputs+1.\n");
      printOutTS(PL_INFO,
           "              nSamples reset to be 10*(nInputs+1).\n");
      nSamples_ = 10 * (nInputs_ + 1);
   }

   if (initLevel != 0) return 0;

   if (psConfig_ != NULL)
   {
      cString = psConfig_->getParameter("MOAT_P");
      if (cString != NULL)
      {
         sscanf(cString, "%s %s %d",winput1,winput2,&P_);
         P_ = P_ / 2 * 2;
         if (P_ <= 0 || P_ > 100) P_ = 4;
         printOutTS(PL_INFO,"MOATSampling: P set to %d (config)\n", P_);
         setFlag = 1;
      }
      cString = psConfig_->getParameter("MOAT_partition_file");
      if (cString != NULL)
      {
         sscanf(cString, "%s %s %s",winput1,winput2,partitionFile);
         printOutTS(PL_INFO,
              "MOATSampling: use MOAT input partition file %s.\n",
              partitionFile);
         fp = fopen(partitionFile, "r");
         if (fp != NULL)
         {
            if (inputSubset_ != NULL) delete [] inputSubset_;
            inputSubset_ = NULL;
            fscanf(fp, "%d", &ss);
            if (ss <= 0 || ss >= nInputs_)
            {
               printOutTS(PL_INFO,
                    "MOATSampling: invalid MOAT input partition file.\n");
               printOutTS(PL_INFO,
                    "              The first line should be nInputs.\n");
               fclose(fp);
               fp = NULL;
            }
            else
            {
               inputSubset_ = new int[nInputs_];
               for (ii = 0; ii < nInputs_; ii++) inputSubset_[ii] = 0;
               for (ii = 0; ii < ss; ii++)
	       {
		  fscanf(fp, "%d", &ii2);
                  if (ii2 < 1 || ii2 > nInputs_)
                  {
                     printOutTS(PL_INFO,
                          "MOATSampling: invalid input partition file.\n");
                     printOutTS(PL_INFO,
                          "               invalid input index %d.\n", ii);
                     delete [] inputSubset_;
                     inputSubset_ = NULL;
                     if(fp != NULL) fclose(fp);
                     fp = NULL;
                     break;
                  }
                  else inputSubset_[ii2-1] = 1;
               }
            }
            if (fp != NULL) fclose(fp);
         }
      }
   }
   if (psSamExpertMode_ == 1 && setFlag == 0)
   {
      printOutTS(PL_INFO,"MOATSampling: the current P is %d.\n", P_);
      sprintf(winput1, "Please choose a new P: (4 - 10, even) ");
      P_ = getInt(4, 10, winput1);
      P_ = P_ / 2 * 2;
   }

   if (printLevel_ > 4)
   {
      printOutTS(PL_INFO,"MOATSampling: initialize: nSamples  = %d\n", 
                 nSamples_);
      printOutTS(PL_INFO,"MOATSampling: initialize: nInputs   = %d\n", 
                 nInputs_);
      printOutTS(PL_INFO,"MOATSampling: initialize: nOutputs  = %d\n", 
                 nOutputs_);
      printOutTS(PL_INFO,"MOATSampling: initialize: numLevels = %d\n",P_);
      if (randomize != 0)
           printOutTS(PL_INFO,"MOATSampling: initialize: randomize on\n");
      else printOutTS(PL_INFO,"MOATSampling: initialize: randomize off\n");
   }

   if (nInputs_ > 100)
   {
      printOutTS(PL_INFO,"MOATSampling: nInputs > 100, use fast version.\n");
      initializeHighDimension();
      return 0;
   }

   nReps = nSamples_ / (nInputs_ + 1);
   if (nReps > maxReps) maxReps = nReps;
   maxSamples = (nInputs_ + 1) * maxReps;
   BS = new double*[maxSamples];
   for (ii = 0; ii < maxSamples; ii++) BS[ii] = new double[nInputs_];

   for (rr = 0; rr < maxReps; rr++) generate(&BS[rr*(nInputs_+1)]);

   for (rr = 1; rr < nReps; rr++) 
   {
      if (printLevel_ > 0)
         printOutTS(PL_INFO,
              "MOATSampling::generate: finding path %d (out of %d)\n", 
              rr+1, nReps);
      maxDist = 0;
      index = rr;
      base1 = (rr - 1) * (nInputs_ + 1);
      for (ss = rr; ss < maxReps; ss++)
      {
         dDist = 0;
         base2 = ss * (nInputs_ + 1);
         for (kk1 = 0; kk1 <= nInputs_; kk1++) 
         {
            for (kk2 = 0; kk2 <= nInputs_; kk2++) 
            {
               for (ii = 0; ii < nInputs_; ii++) 
               {
                  ddata = BS[base1+kk1][ii] - BS[base2+kk2][ii];
                  dDist += ddata * ddata;
               }
            }
         }
         if (dDist > maxDist)
         {
            maxDist = dDist;
            index = ss;
         }
      }
      if (index != rr)
      { 
         base1 = rr * (nInputs_ + 1);
         base2 = index * (nInputs_ + 1);
         for (kk1 = 0; kk1 <= nInputs_; kk1++) 
         {
            for (ii = 0; ii < nInputs_; ii++)
            {
               ddata = BS[base1+kk1][ii];
               BS[base1+kk1][ii] = BS[base2+kk1][ii];
               BS[base2+kk1][ii] = ddata;
            }
         }
      }
   }

   allocSampleData();
   ranges = new double[nInputs_];
   for (ii = 0;  ii < nInputs_;  ii++) 
      ranges[ii] = upperBounds_[ii] - lowerBounds_[ii];

   if (inputSubset_ != NULL)
   {
      nSub = 0;
      for (ii = 0; ii < nInputs_; ii++)
         if (inputSubset_[ii] == 1) nSub++;
   }
   else nSub = nInputs_;

   for (ss = 0; ss < nSamples_; ss+=(nInputs_+1))
   {
      for (ii = 0; ii <= nInputs_; ii++)
      {
         for (ii2 = 0; ii2 < nInputs_; ii2++)
         {
            ddata = BS[ss+ii][ii2];
            ddata = ddata * ranges[ii2] + lowerBounds_[ii2];
            sampleMatrix_[ss+ii][ii2] = ddata;
         }
      }
      for (ii = nSub+1; ii <= nInputs_; ii++)
         sampleStates_[ss+ii] = 1;
   }

   if (repair(NULL,0) != 0)
   {
      if (checkSample(nInputs_, nSamples_, sampleMatrix_) != 0)
      {
         printOutTS(PL_ERROR,
              "MOATSampling: generated sample is not MOAT.\n");
         exit(1);
      }
   }

   if (printLevel_ > 2)
   {
      tempX = new double[2*nReps];
      bins  = new int[P_];
      for (ii = 0; ii < P_; ii++) bins[ii] = 0;
      for (ii = 0; ii < nInputs_; ii++)
      {
         for (rr = 0; rr < nReps; rr++)
         {
            tempX[rr*2] = sampleMatrix_[rr*(nInputs_+1)][ii];
            for (ii2 = 1; ii2 < nInputs_+1; ii2++)
            {
               if (sampleMatrix_[rr*(nInputs_+1)+ii2][ii] != tempX[2*rr])
               {
                  tempX[rr*2+1] = sampleMatrix_[rr*(nInputs_+1)+ii2][ii];
                  break;
               }
            }
         }
         sortDbleList(2*nReps, tempX);
         nn = 1;
         currBin = 0;
         for (ii2 = 1; ii2 < 2*nReps; ii2++)
         {
            if (tempX[ii2] != tempX[ii2-1])
            {
               printf("MOAT: input %3d - level = %12.4e, # times = %4d\n",
                      ii+1, tempX[ii2-1], nn);
               bins[currBin++] += nn;
               nn = 1;
            } else nn++;
         }
         printf("MOAT: input %3d - level = %12.4e, # times = %4d\n",
                ii+1, tempX[ii2-1], nn);
         bins[currBin++] += nn;
      }
      for (ii = 0; ii < P_; ii++) 
         printf("MOAT: frequency of visit to bin %5d = %d\n",ii+1,bins[ii]);
      delete [] bins;
      delete [] tempX;
   }
   
   delete [] ranges;
   for (ii = 0;  ii < maxSamples; ii++) delete [] BS[ii];
   delete [] BS;
   return 0;
}

//*************************************************************************
//* generate the BS matrix
//*------------------------------------------------------------------------
int MOATSampling::generate(double **BS)
{
   int    *permute, ss, ii, ii2, idata, *subset, nSub, imax;
   double **B, *D, *X, **B2, delta;

   delta = P_ / ((double) (2*P_) - 2.0);

   B = new double*[nInputs_+1];
   for (ii = 0; ii <= nInputs_; ii++)
   {
      B[ii] = new double[nInputs_];
      for (ii2 = 0; ii2 < ii; ii2++) B[ii][ii2] = 1.0;
      for (ii2 = ii; ii2 < nInputs_; ii2++) B[ii][ii2] = 0.0;
   }
   D = new double[nInputs_];
   X = new double[nInputs_];
   permute = new int[nInputs_];
   B2   = new double*[nInputs_+1];
   for (ii = 0; ii <= nInputs_; ii++) B2[ii] = new double[nInputs_];

   for (ii = 0; ii < nInputs_; ii++)
   {
      D[ii] = PSUADE_drand();
      if (D[ii] > 0.5) D[ii] = 1.0;
      else             D[ii] = -1.0;
   }

   imax = (P_ - 1) / 2;
   for (ii = 0; ii < nInputs_; ii++)
   {
       X[ii] = PSUADE_drand();
       idata = (int) (X[ii] * (imax + 1));
       if (idata > imax) idata--;
       X[ii] = (double) idata / (double) (P_ - 1);
   }
   
   if (inputSubset_ == NULL)
   {
      generateRandomIvector(nInputs_, permute);
   }
   else
   {
      subset = new int[nInputs_];
      nSub = 0;
      for (ii = 0; ii < nInputs_; ii++)
         if (inputSubset_[ii] == 1) nSub++;
      generateRandomIvector(nSub, subset);
      generateRandomIvector(nInputs_-nSub, permute);
      for (ii = 0; ii < nInputs_; ii++)
         subset[nSub+ii] = permute[ii] + nSub;

      ss = 0;
      for (ii = 0; ii < nInputs_; ii++)
      {
         if (inputSubset_[ii] == 1)
         {
            permute[ii] = subset[ss];
            ss++;
         }
      }
      for (ii = 0; ii < nInputs_; ii++)
      {
         if (inputSubset_[ii] != 1)
         {
            permute[ii] = subset[ss];
            ss++;
         }
      }
      delete [] subset;
   } 
   
   for (ii = 0; ii <= nInputs_; ii++)
      for (ii2 = 0; ii2 < nInputs_; ii2++)
         B2[ii][ii2] = X[ii2]+delta/2*((B[ii][ii2]*2-1.0)*D[ii2]+1.0);
   for (ii = 0; ii <= nInputs_; ii++)
      for (ii2 = 0; ii2 < nInputs_; ii2++)
         BS[ii][ii2] = B2[ii][permute[ii2]];

   for (ii = 0;  ii <= nInputs_; ii++)
   {
      delete [] B[ii];
      delete [] B2[ii];
   }
   delete [] B;
   delete [] B2;
   delete [] D;
   delete [] X;
   delete [] permute;
   return 0;
}

//*************************************************************************
//* refine the sample space
//*------------------------------------------------------------------------
int MOATSampling::refine(int refineRatio, int randomize, double thresh,
                         int nSamples, double *sampleErrors)
{
   int    ss, ii, ii2, jj, rr, *newStates, nReps, nTimes, newNSamples;
   double **BS, **newSamples, *newOutputs, *ranges, ddata;

   (void) randomize;
   (void) thresh;
   (void) nSamples;
   (void) sampleErrors;

   if (inputSubset_ != NULL)
   {
      printOutTS(PL_INFO,"MOATSampling: refine is not available due to the\n");
      printOutTS(PL_INFO,"              use of selective replications.\n");
      return 0;
   }
   if (refineRatio != 2)
      printOutTS(PL_INFO,"MOATSampling WARNING: refinement ratio set to 2.\n");

   nTimes = 2;

   newNSamples = nSamples_ * nTimes;
   newSamples = new double*[newNSamples];
   for (ss = 0;  ss < newNSamples; ss++)
      newSamples[ss] = new double[nInputs_];
   newOutputs = new double[newNSamples*nOutputs_];
   newStates = new int[newNSamples];
   for (ss = 0;  ss < newNSamples; ss++)
   {
      newStates[ss] = 0;
      for (jj = 0; jj < nOutputs_; jj++)
         newOutputs[ss*nOutputs_+jj] = PSUADE_UNDEFINED;
   }

   for (ss = 0; ss < nSamples_; ss++) 
   {
      for (ii = 0; ii < nInputs_; ii++)
         newSamples[ss][ii] = sampleMatrix_[ss][ii];
      for (jj = 0; jj < nOutputs_; jj++)
         newOutputs[ss*nOutputs_+jj] = sampleOutput_[ss*nOutputs_+jj];
      newStates[ss] = 1;
   }

   BS = new double*[nSamples_];
   for (ii = 0; ii < nSamples_; ii++) BS[ii] = new double[nInputs_];

   nReps  = nSamples_ / (nInputs_ + 1);
   ranges = new double[nInputs_];
   for (ii = 0;  ii < nInputs_;  ii++) 
      ranges[ii] = upperBounds_[ii] - lowerBounds_[ii];

   for (rr = 0; rr < nReps; rr++) generate(&BS[rr*(nInputs_+1)]);
   for (ss = nSamples_; ss < newNSamples; ss+=(nInputs_+1))
   {
      for (ii = 0; ii <= nInputs_; ii++)
      {
         for (ii2 = 0; ii2 < nInputs_; ii2++)
         {
            ddata = BS[ss-nSamples_+ii][ii2];
            ddata = ddata * ranges[ii2] + lowerBounds_[ii2];
            newSamples[ss+ii][ii2] = ddata;
         }
      }
   }


   deleteSampleData();
   delete [] ranges;
   for (ii = 0;  ii < nSamples_; ii++) delete [] BS[ii];
   delete [] BS;

   nSamples_ = nSamples_ * nTimes;
   sampleMatrix_ = newSamples;
   sampleOutput_ = newOutputs;
   sampleStates_ = newStates;

   if (checkSample(nInputs_, nSamples_, sampleMatrix_) != 0)
   {
      printOutTS(PL_ERROR,"MOATSampling: refined sample is not MOAT.\n");
      exit(1);
   }

   if (printLevel_ > 4)
   {
      printOutTS(PL_INFO,"MOATSampling::refine: nSamples = %d\n",nSamples_);
      printOutTS(PL_INFO,"MOATSampling::refine: nInputs  = %d\n",nInputs_);
      printOutTS(PL_INFO,"MOATSampling::refine: nOutputs = %d\n",nOutputs_);
      if (randomize != 0)
           printOutTS(PL_INFO,"MOATSampling::refine: randomize on\n");
      else printOutTS(PL_INFO,"MOATSampling::refine: randomize off\n");
   }

   if (repair(NULL, nSamples_/nTimes) != 0)
   {
      if (checkSample(nInputs_, nSamples_, sampleMatrix_) != 0)
      {
         printOutTS(PL_ERROR,"MOATSampling: refined sample is not MOAT.\n");
         exit(1);
      }
   }
   return 0;
}

//*************************************************************************
//* repair a MOAT design due to constraints
//*------------------------------------------------------------------------
int MOATSampling::repair(char *fname, int start)
{
   int    nPatterns, nInps, nSets, ii, jj, kk, *inpList, pindex, iindex;
   double **patterns, *tempSample;
   char   inStr[200], *cString, winput1[200], winput2[200], repairFile[200];
   FILE   *fp=NULL;

   if (start / (nInputs_+1) * (nInputs_ + 1) != start)
   {
      printOutTS(PL_ERROR,
           "MOATSampling: start should be multiples of nInputs+1.\n");
      exit(1);
   }
  
   fp = NULL;
   if (fname != NULL)
   {
      fp = fopen(fname, "r");
      if (fp != NULL)
            printf("MOAT repair file found = %s.\n", fname);
   }
   if (fp == NULL && psConfig_ != NULL)
   {
      cString = psConfig_->getParameter("MOAT_repair_file");
      if (cString != NULL)
      {
         sscanf(cString, "%s %s %s",winput1,winput2,repairFile);
         fp = fopen(repairFile, "r");
         if (fp != NULL)
            printf("MOAT repair file found = %s.\n", repairFile);
      }
   }
   if (fp == NULL)
   {
      if (fname != NULL)
         printOutTS(PL_INFO,"MOAT repair: repair file not found.\n");
      return 1;
   }
  
   fscanf(fp, "%s", inStr);
   if (strcmp(inStr, "BEGIN"))
   {
      printOutTS(PL_ERROR,"MOATSampling: wrong format in repair file.\n");
      printOutTS(PL_ERROR,"First  line : BEGIN\n");
      printOutTS(PL_ERROR,"Second line : nPatterns nInputs.\n");
      printOutTS(PL_ERROR,"Third  line : a list of input IDs (1-based).\n");
      printOutTS(PL_ERROR,"Fourth line : (and on) set of patterns.\n");
      printOutTS(PL_ERROR,"Last   line : END\n");
      fclose(fp);
      exit(1);
   }

   fscanf(fp, "%d %d", &nPatterns, &nInps);
   if (nPatterns <= 0 || nInps <= 0)
   {
      printOutTS(PL_ERROR,"MOATSampling: nPatterns or nInps <= 0.\n");
      fclose(fp);
      exit(1);
   }
   nSets = nPatterns / (nInps + 1);
   if (nSets*(nInps+1) != nPatterns)
   {
      printOutTS(PL_ERROR,
           "MOATSampling: nPatterns should be multiples of nInputs+1.\n");
      fclose(fp);
      exit(1);
   }
   printf("nPatterns = %d involving %d inputs\n", nPatterns, nInps);

   inpList = new int[nInps];
   for (ii = 0;  ii < nInps; ii++)
   {
      fscanf(fp, "%d", &inpList[ii]);
      if (inpList[ii] <= 0 || inpList[ii] > nInputs_)
      {
         printOutTS(PL_ERROR,
              "MOATSampling ERROR: input index out of range (%d,%d)\n",
              inpList[ii], nInputs_);
         fclose(fp);
         exit(1);
      }
      for (jj = 0; jj < ii; jj++)
      {
         if (inpList[ii] == inpList[jj])
         {
            printOutTS(PL_ERROR,
                 "MOATSampling ERROR: repeated index (%d)\n",inpList[ii]);
            fclose(fp);
            return 1;
         }
      }
   }

   patterns = new double*[nPatterns];
   for (ii = 0;  ii < nPatterns; ii++)
   {
      patterns[ii] = new double[nInps];
      for (jj = 0; jj < nInps; jj++) fscanf(fp, "%lg", &patterns[ii][jj]);
   }
   fscanf(fp, "%s", inStr);
   fclose(fp);

   if (strcmp(inStr, "END"))
   {
      printOutTS(PL_ERROR,"MOATSampling ERROR: wrong format in repair file.\n");
      printOutTS(PL_ERROR,"The file should end with END\n");
      exit(1);
   }

   if (checkSample(nInps, nPatterns, patterns) != 0)
   {
      printOutTS(PL_ERROR,
           "MOATSampling ERROR: pattern in repair file is not MOAT.\n");
      exit(1);
   }

   tempSample = new double[nInps];
   for (ii = start; ii < nSamples_; ii+=nInputs_+1)
   {
      pindex = ii / (nInputs_ + 1) * (nInps + 1);

      for (kk = 0; kk < nInps; kk++)
      {
         iindex = inpList[kk] - 1;
         tempSample[kk] = sampleMatrix_[ii][iindex];
         sampleMatrix_[ii][iindex] = patterns[pindex][kk];
      }

      for (jj = 1; jj <= nInputs_; jj++)
      {
         for (kk = 0; kk < nInps; kk++)
         {
            iindex = inpList[kk] - 1;
            if (sampleMatrix_[ii+jj][iindex] != tempSample[kk]) break; 
         }

         if (kk < nInps) pindex++;

         for (kk = 0; kk < nInps; kk++)
         {
            iindex = inpList[kk] - 1;
            tempSample[kk] = sampleMatrix_[ii+jj][iindex];
            sampleMatrix_[ii+jj][iindex] = patterns[pindex][kk]; 
         }
      }
   }

   if (checkSample(nInputs_, nSamples_, sampleMatrix_) != 0)
   {
      printOutTS(PL_ERROR,"MOATSampling ERROR: repaired file is not MOAT.\n");
      exit(1);
   }

   delete [] inpList;
   for (ii = 0; ii < nPatterns; ii++) delete [] patterns[ii];
   delete [] patterns;
   delete [] tempSample;
   return 0;
}

//*************************************************************************
//* merge two MOAT designs (FOR 2 SETS OF INPUTS)
//*------------------------------------------------------------------------
int MOATSampling::merge()
{
   int        nInps1, nOuts1, nSamp1, nInps2, nSamp2, nOutputs, count; 
   int        samplingMethod, ii, ii2, rr, nSamples, nInputs, cnt1, cnt2;
   int        *sampleStates, nReps;
   double     *sampleInputs1;
   double     *sampleInputs2, *sampleInputs, *sampleOutputs;
   char       **inpNames1, **inpNames2, **inpNames, file1[500], file2[500];
   char       pString[500];
   FILE       *fp1, *fp2;
   PsuadeData *psuadeIO1, *psuadeIO2;
   pData      pPtr1, pINames1;
   pData      pPtr2, pINames2;

   sprintf(pString,"Please enter the name of the first MOAT datafile: ");
   getString(pString, file1);
   file1[strlen(file1)-1] = '\0';
   if ((fp1=fopen(file1,"r")) == NULL)
   {
      printOutTS(PL_ERROR,"ERROR : File %s not found.\n", file1);
      return 1;
   }
   else fclose(fp1);
   psuadeIO1 = new PsuadeData();
   psuadeIO1->setOutputLevel(0);
   if (psuadeIO1->readPsuadeFile(file1) != 0)
   {
      printOutTS(PL_ERROR,
           "MOAT ERROR : problem with reading file %s.\n", file1);
      delete psuadeIO1;
      return 1;
   }
   assert(psuadeIO1->getParameter("input_ninputs", pPtr1) == 0);
   nInps1 = pPtr1.intData_;
   assert(psuadeIO1->getParameter("input_names", pINames1) == 0);
   inpNames1 = pINames1.strArray_;
   assert(psuadeIO1->getParameter("output_noutputs", pPtr1) == 0);
   nOuts1 = pPtr1.intData_;
   assert(psuadeIO1->getParameter("method_sampling", pPtr1) == 0);
   samplingMethod = pPtr1.intData_;
   if (samplingMethod != PSUADE_SAMP_MOAT &&
       samplingMethod != PSUADE_SAMP_GMOAT)
   {
      printOutTS(PL_ERROR,"MOAT Merge ERROR : data1 is not MOAT.\n");
      delete psuadeIO1;
      pINames1.clean();
      return 1;
   }
   assert(psuadeIO1->getParameter("method_nsamples", pPtr1) == 0);
   nSamp1 = pPtr1.intData_;
   assert(psuadeIO1->getParameter("input_sample", pPtr1) == 0);
   sampleInputs1 = pPtr1.dbleArray_;
   if (checkSample2(nInps1, nSamp1, sampleInputs1) != 0)
   {
      printOutTS(PL_ERROR,
           "MOAT Merge ERROR : first sample is not MOAT.\n");
      delete psuadeIO1;
      pINames1.clean();
      return 1;
   }

   sprintf(pString,"Please enter the name of the second MOAT datafile: ");
   getString(pString, file2);
   file2[strlen(file2)-1] = '\0';
   if ((fp2=fopen(file2,"r")) == NULL)
   {
      printOutTS(PL_ERROR,"MOAT ERROR : File %s not found.\n", file2);
      delete psuadeIO1;
      pINames1.clean();
      return 1;
   }
   else fclose(fp2);
                                                                  
   psuadeIO2 = new PsuadeData();
   psuadeIO2->setOutputLevel(0);
   if (psuadeIO2->readPsuadeFile(file2) != 0)
   {
      printOutTS(PL_ERROR,
           "MOAT ERROR : problem with reading file %s.\n", file2);
      delete psuadeIO1;
      delete psuadeIO2;
      pINames1.clean();
      return 1;
   }
   assert(psuadeIO2->getParameter("input_ninputs", pPtr2) == 0);
   nInps2 = pPtr2.intData_;
   assert(psuadeIO2->getParameter("input_names", pINames2) == 0);
   inpNames2 = pINames2.strArray_;
   assert(psuadeIO2->getParameter("method_sampling", pPtr2) == 0);
   samplingMethod = pPtr2.intData_;
   if (samplingMethod != PSUADE_SAMP_MOAT &&
       samplingMethod != PSUADE_SAMP_GMOAT)
   {
      printOutTS(PL_ERROR,"MOAT Merge ERROR : data2 is not MOAT.\n");
      delete psuadeIO1;
      delete psuadeIO2;
      pINames1.clean();
      pINames2.clean();
      return 1;
   }
   assert(psuadeIO2->getParameter("method_nsamples", pPtr2) == 0);
   nSamp2 = pPtr2.intData_;
   assert(psuadeIO2->getParameter("input_sample", pPtr2) == 0);
   sampleInputs2 = pPtr2.dbleArray_;
   if (checkSample2(nInps2, nSamp2, sampleInputs2) != 0)
   {
      printOutTS(PL_ERROR,
           "MOAT Merge ERROR : second sample is not MOAT.\n");
      delete psuadeIO1;
      delete psuadeIO2;
      pINames1.clean();
      pINames2.clean();
      return 1;
   }

   nReps = nSamp1 / (nInps1 + 1);
   if (nReps != (nSamp2 / (nInps2 + 1)))
   {
      printOutTS(PL_ERROR,
           "MOAT Merge ERROR : different number of replications.\n");
      delete psuadeIO1;
      delete psuadeIO2;
      pINames1.clean();
      pINames2.clean();
      return 1;
   }

   nSamples = nSamp1 + nSamp2 - nReps;
   nInputs  = nInps1 + nInps2;
   nOutputs = nOuts1;
   sampleInputs  = new double[nSamples*nInputs];
   sampleOutputs = new double[nSamples*nOutputs];
   sampleStates  = new int[nSamples];
   for (rr = 0; rr < nReps; rr++)
   {
      count = rr * (nInputs + 1);
      cnt1  = rr * (nInps1 + 1);
      cnt2  = rr * (nInps2 + 1);
      for (ii = 0; ii < nInps1; ii++)
         sampleInputs[count*nInputs+ii] = sampleInputs1[cnt1*nInps1+ii];
      for (ii = 0; ii < nInps2; ii++)
         sampleInputs[count*nInputs+nInps1+ii] = sampleInputs2[cnt2*nInps2+ii];
      for (ii = 0; ii < nInps1; ii++)
      {
         count++;
         cnt1++;
         for (ii2 = 0; ii2 < nInps1; ii2++)
            sampleInputs[count*nInputs+ii2] = sampleInputs1[cnt1*nInps1+ii2];
         for (ii2 = 0; ii2 < nInps2; ii2++)
            sampleInputs[count*nInputs+nInps1+ii2] = 
               sampleInputs2[cnt2*nInps2+ii2];
      }
      for (ii = 0; ii < nInps2; ii++)
      {
         count++;
         cnt2++;
         for (ii2 = 0; ii2 < nInps1; ii2++)
            sampleInputs[count*nInputs+ii2] = sampleInputs1[cnt1*nInps1+ii2];
         for (ii2 = 0; ii2 < nInps2; ii2++)
            sampleInputs[count*nInputs+nInps1+ii2] = 
               sampleInputs2[cnt2*nInps2+ii2];
      }
   }
   for (ii = 0; ii < nSamples; ii++) sampleStates[ii] = 0;
   for (ii = 0; ii < nSamples*nOutputs; ii++)
      sampleOutputs[ii] = PSUADE_UNDEFINED;
   inpNames = new char*[nInputs];
   for (ii = 0; ii < nInps1; ii++)
   {
      inpNames[ii] = new char[100]; 
      strcpy(inpNames[ii], inpNames1[ii]);
   }
   for (ii = 0; ii < nInps1; ii++)
   {
      inpNames[nInps1+ii] = new char[100]; 
      strcpy(inpNames[nInps1+ii], inpNames2[ii]);
   }
   psuadeIO1->updateInputSection(nSamples,nInputs,NULL,NULL,NULL,
                                 sampleInputs,inpNames,NULL,NULL,NULL,NULL);
   psuadeIO1->updateOutputSection(nSamples,nOutputs,sampleOutputs,
                                  sampleStates,NULL);
   psuadeIO1->updateMethodSection(-1,nSamples,-1,-1,-1);
   psuadeIO1->writePsuadeFile(NULL,0);
   
   delete psuadeIO1;
   delete psuadeIO2;
   pINames1.clean();
   pINames2.clean();
   delete [] sampleInputs;
   delete [] sampleOutputs;
   delete [] sampleStates;
   for (ii = 0; ii < nInputs; ii++) delete [] inpNames[ii];
   delete [] inpNames;
   return 0;
}

//*************************************************************************
//* check whether the sample is MOAT (return 0 if yes)
//*------------------------------------------------------------------------
int MOATSampling::checkSample(int nInputs, int nSamples, double **X)
{
   int    ss, ii, ii2, *errArray, nDiff;
   double xtemp1, xtemp2;

   errArray = new int[nInputs];
   for (ss = 0; ss < nSamples; ss+=(nInputs+1))
   {
      for (ii = 0; ii < nInputs; ii++) errArray[ii] = 0;
      for (ii = 1; ii <= nInputs; ii++)
      {
         for (ii2 = 0; ii2 < nInputs; ii2++)
         {
            xtemp1 = X[(ss+ii-1)][ii2];
            xtemp2 = X[(ss+ii)][ii2];
            if (xtemp1 != xtemp2) errArray[ii2]++;
         }
      }
      nDiff = 0;
      for (ii = 0; ii < nInputs; ii++) nDiff += errArray[ii];
      if (nDiff != nInputs)
      {
         // Cleanup by Bill Oliver
	 delete [] errArray;
	 return 1;
      }
   }
   delete [] errArray;
   return 0;
}

//*************************************************************************
//* check whether the sample is MOAT (return 0 if yes)
//*------------------------------------------------------------------------
int MOATSampling::checkSample2(int nInputs, int nSamples, double *X)
{
   int    ss, ii, ii2, *errArray, nDiff;
   double xtemp1, xtemp2;

   errArray = new int[nInputs];
   for (ss = 0; ss < nSamples; ss+=(nInputs+1))
   {
      for (ii = 0; ii < nInputs; ii++) errArray[ii] = 0;
      for (ii = 1; ii <= nInputs; ii++)
      {
         for (ii2 = 0; ii2 < nInputs; ii2++)
         {
            xtemp1 = X[(ss+ii-1)*nInputs+ii2];
            xtemp2 = X[(ss+ii)*nInputs+ii2];
            if (xtemp1 != xtemp2) errArray[ii2]++;
         }
      }
      nDiff = 0;
      for (ii = 0; ii < nInputs; ii++) nDiff += errArray[ii];
      if (nDiff != nInputs)
      {
         // Cleanup by Bill Oliver
	 delete [] errArray;
	 return 1;
      }
   }
   delete [] errArray;
   return 0;
}

//*************************************************************************
//* generate sample with constraints
//*------------------------------------------------------------------------
int MOATSampling::genRepair(int nInputs, double *lbounds, double *ubounds)
{
   int        status, nFiles, faFlag, kk, jj, ii, nPaths, currP, nTrials;
   int        *indSet, count, trial, iInd, ind, ind2, *states;
   double     *threshLs, *threshUs, **moatSample, *tempW, dtemp, ddata;
   double     currX1, currX2, currX1F, currX2F, filterRange, currY1, currY2;
   double     sLo, sHi;
   char       pString[1000], winput[1000];
   FILE       *fp;
   FuncApprox **faPtrs;
   PsuadeData *ioPtr=NULL;
   pData      pPtr, pLower, pUpper;

   status = 0;
   sprintf(pString,"How many constraint data files are there (1-10)? ");
   nFiles = getInt(1, 10, pString);
   faFlag = 2;
   faPtrs = new FuncApprox*[nFiles];
   threshLs = new double[nFiles];
   threshUs = new double[nFiles];
   for (kk = 0; kk < nFiles; kk++)
   {
      sprintf(pString,"Enter name of file #%d : ", kk+1);
      getString(pString, winput);
      ioPtr = new PsuadeData;
      status = ioPtr->readPsuadeFile(winput);
      if (status != 0)
      {
         printOutTS(PL_ERROR,"moatgen READ ERROR: file = %s\n", winput);
         exit(1);
      }
      ioPtr->getParameter("input_ninputs", pPtr);
      jj = pPtr.intData_;
      if (jj != nInputs)
      {
         printOutTS(PL_ERROR,"moatgen ERROR: nInputs mismatch.\n");
         exit(1);
      }
      pLower.clean();
      ioPtr->getParameter("input_lbounds", pLower);
      for (ii = 0; ii < nInputs; ii++)
      {
         if (lbounds[ii] != pLower.dbleArray_[ii])
         {
            printOutTS(PL_ERROR,
                 "MOAT genRepair ERROR: lower bound mismatch.\n");
            exit(1);
         }
      }
      pUpper.clean();
      ioPtr->getParameter("input_ubounds", pUpper);
      for (ii = 0; ii < nInputs; ii++)
      {
         if (ubounds[ii] != pUpper.dbleArray_[ii])
         {
            printOutTS(PL_ERROR,"moatgen ERROR: upper bound mismatch.\n");
            exit(1);
         }
      }
      faPtrs[kk] = genFAInteractive(ioPtr, faFlag);
      if (faPtrs[kk] == NULL) {printf("ERROR detected.\n"); exit(1);}
      faPtrs[kk]->setOutputLevel(printLevel_);
      sprintf(pString,"Constraint %d lower bound : ",kk+1);
      threshLs[kk] = getDouble(pString);
      sprintf(pString,"Constraint %d upper bound : ",kk+1);
      threshUs[kk] = getDouble(pString);
      if (threshLs[kk] >= threshUs[kk])
      {
         printf("ERROR : lower bound >= upper bound.\n");
         exit(1);;
      }
      delete ioPtr;
   }
   sprintf(pString,"Please enter the number of paths to search: ");
   nPaths = getInt(1, 1000, pString);
   sprintf(pString, "Please enter P (resolution: try 4-10) : ");
   currP = getInt(4, 10, pString);
   sprintf(pString, "Please enter the number of trials (> 100) : ");
   nTrials = getInt(101, 10000000, pString);
   moatSample = new double*[nPaths*(nInputs+1)];
   for (ii = 0; ii < nPaths*(nInputs+1); ii++)
      moatSample[ii] = new double[nInputs];
   tempW = new double[nInputs];
   indSet = new int[nInputs];
   count = 0;
   for (ii = 0; ii < nPaths; ii++)
   {
      trial = 0; 
      while (trial < nTrials)
      {
         iInd = count;
         trial++;
         for (jj = 0; jj < nInputs; jj++)
         {
            ind = PSUADE_rand() % currP;
            dtemp = ind * (ubounds[jj] - lbounds[jj]) / (currP - 1.0);
            moatSample[iInd][jj] = lbounds[jj] + dtemp;
         }
         for (kk = 0; kk < nFiles; kk++)
         {
            dtemp = faPtrs[kk]->evaluatePoint(moatSample[iInd]);
            if (dtemp < threshLs[ii] || dtemp > threshUs[ii]) break;
         }
         if (kk != nFiles) continue;

         generateRandomIvector(nInputs, indSet);

         for (jj = 0; jj < nInputs; jj++)
         {
            iInd++;
            for (kk = 0; kk < nInputs; kk++)
               moatSample[iInd][kk] = moatSample[iInd-1][kk];

            ind2 = indSet[jj];

            ddata = moatSample[iInd][ind2]; 
            currX1F = - PSUADE_UNDEFINED;
            currX2F =   PSUADE_UNDEFINED;
            for (kk = 0; kk < nFiles; kk++)
            {
               filterRange = threshUs[kk] - threshLs[kk];
               moatSample[iInd][ind2] = lbounds[ind2];
               currY1 = faPtrs[kk]->evaluatePoint(moatSample[iInd]);
               moatSample[iInd][ind2] = ubounds[ind2];
               currY2 = faPtrs[kk]->evaluatePoint(moatSample[iInd]);
               currX1 = lbounds[ind2];
               currX2 = ubounds[ind2];

               if (currY2 >= threshUs[kk] && currY1 >= threshUs[kk])
                  currX1 = currX2 = 0.0;
               else if (currY2 <= threshLs[kk] && currY1 <= threshLs[kk])
                  currX1 = currX2 = 0.0;
               else if (currY2 > currY1)
               {
                  if (currY2 <= threshUs[kk]) currX2 = ubounds[ind2];
                  else
                  {
                     sLo = lbounds[ind2];
                     sHi = ubounds[ind2];
                     while (PABS((currY2-threshUs[kk])/filterRange)>1e-4)
                     {
                        moatSample[iInd][ind2] = 0.5 * (sLo + sHi);
                        currY2=faPtrs[kk]->evaluatePoint(moatSample[iInd]);
                        if (currY2 > threshUs[kk]) sHi = 0.5 * (sLo+sHi);
                        else                       sLo = 0.5 * (sLo+sHi);
                     }
                     currX2 = moatSample[iInd][ind2];
                  }
                  if (currY1 >= threshLs[kk]) currX1 = lbounds[ind2];
                  else
                  {
                     sLo = lbounds[ind2];
                     sHi = ubounds[ind2];
                     while (PABS((currY1-threshLs[kk])/filterRange)>1e-4)
                     {
                        moatSample[iInd][ind2] = 0.5 * (sLo + sHi);
                        currY1=faPtrs[kk]->evaluatePoint(moatSample[iInd]);
                        if (currY1 < threshLs[kk]) sLo = 0.5 * (sLo+sHi);
                        else                       sHi = 0.5 * (sLo+sHi);
                     }
                     currX1 = moatSample[iInd][ind2];
                  }
               }
               else
               {
                  if (currY1 <= threshUs[kk]) currX1 = lbounds[ind2];
                  else
                  {
                     sLo = lbounds[ind2];
                     sHi = ubounds[ind2];
                     while (PABS((currY1-threshUs[kk])/filterRange)>1e-4)
                     {
                        moatSample[iInd][ind2] = 0.5 * (sLo + sHi);
                        currY1=faPtrs[kk]->evaluatePoint(moatSample[iInd]);
                        if (currY1 > threshUs[kk]) sLo = 0.5 * (sLo+sHi);
                        else                       sHi = 0.5 * (sLo+sHi);
                     }
                     currX1 = moatSample[iInd][ind2];
                  }
                  if (currY2 >= threshLs[kk]) currX2 = ubounds[ind2];
                  else
                  {
                     sLo = lbounds[ind2];
                     sHi = ubounds[ind2];
                     while (PABS((currY2-threshLs[kk])/filterRange)>1e-4)
                     {
                        moatSample[iInd][ind2] = 0.5 * (sLo + sHi);
                        currY2=faPtrs[kk]->evaluatePoint(moatSample[iInd]);
                        if (currY2 < threshLs[kk]) sHi = 0.5 * (sLo+sHi);
                        else                       sLo = 0.5 * (sLo+sHi);
                     }
                     currX2 = moatSample[iInd][ind2];
                  }
               }
               if (PABS(currX2-currX1)<0.1*(ubounds[ind2]-lbounds[ind2])) 
                  break;
               if (currX1 > currX1F) currX1F = currX1;
               if (currX2 < currX2F) currX2F = currX2;
            }
            moatSample[iInd][ind2] = ddata;
            if (kk != nFiles) break;
            tempW[ind2] = PABS(currX2F - currX1F) / (currP-1.0);
            moatSample[iInd][ind2] += tempW[ind2];
            if (moatSample[iInd][ind2] > ubounds[ind2])
                 dtemp = threshLs[kk] - 1.0;
            else
            {
               for (kk = 0; kk < nFiles; kk++)
               {
                  moatSample[iInd][ind2] = ddata;
                  moatSample[iInd][ind2] += tempW[ind2];
                  dtemp = faPtrs[kk]->evaluatePoint(moatSample[iInd]);
                  if (dtemp < threshLs[kk] || dtemp > threshUs[kk]) 
                  {
                     moatSample[iInd][ind2] -= 2.0 * tempW[ind2];
                     if (moatSample[iInd][ind2] < lbounds[ind2])
                        break;
                     dtemp = faPtrs[kk]->evaluatePoint(moatSample[iInd]);
                     if (dtemp < threshLs[kk] || dtemp > threshUs[kk]) 
                        break;
                  }
               }
               moatSample[iInd][ind2] = ddata;
               if (kk != nFiles) break;
            }
         }
         if (jj == nInputs) 
         {
            count += (nInputs + 1);
            if (printLevel_ > 2)
               printOutTS(PL_INFO,
                    "MOAT genRepair: path %d (out of %d) found.\n", ii+1,
                    nPaths);
            break; 
         }
         else
         {
            if (printLevel_ > 2)
               printOutTS(PL_INFO,
                    "Current path fails (%d out of max %d).\n",
                    trial, nTrials); 
         }
      }
      if (trial >= nTrials)
      {
         printOutTS(PL_INFO,"moatgen fails to find all possible paths.\n");
         printOutTS(PL_INFO,"Suggestion: try a larger P than %d.\n", currP);
         break;
      }
   }
   for (ii = 0; ii < nPaths; ii++)
   {
      for (kk = 0; kk < nFiles; kk++)
      {
         dtemp = faPtrs[kk]->evaluatePoint(moatSample[ii]);
         if (dtemp < threshLs[kk] || dtemp > threshUs[kk])
         printOutTS(PL_ERROR,
            "MOAT genRepair:sample %d fails final test (%e <? %e <? %e).\n",
            ii, threshLs[kk], dtemp, threshUs[kk]);
      }
   }
   delete [] tempW;
   delete [] indSet;
   delete [] threshLs;
   delete [] threshUs;
   for (kk = 0; kk < nFiles; kk++) delete faPtrs[kk];
   delete [] faPtrs;
   if (trial >= nTrials)
   {
      for (ii = 0; ii < nPaths*(nInputs+1); ii++)
         delete [] moatSample[ii];
      delete [] moatSample;
      printOutTS(PL_ERROR,"MOAT genRepair FAILS.\n");
      return 0; 
   }
   fp = fopen("MOAT_repair_file", "w");
   // Add a check for NULL by Bill Oliver
   if(fp != NULL)
   {
      fprintf(fp, "BEGIN\n");
      fprintf(fp, "%d %d\n", nPaths*(nInputs+1), nInputs);
      for (ii = 0; ii < nInputs; ii++) fprintf(fp, "%d ", ii+1);
      fprintf(fp, "\n");
      for (ii = 0; ii < nPaths*(nInputs+1); ii++)
      {
         for (jj = 0; jj < nInputs; jj++)
            fprintf(fp, "%e ", moatSample[ii][jj]);
         fprintf(fp, "\n");
      }
      fprintf(fp, "END\n");
      fclose(fp);
   }
   count = nPaths * (nInputs + 1); 
   tempW = new double[count*nInputs];
   states = new int[count];
   for (ii = 0; ii < count; ii++) states[ii] = 1;
   for (ii = 0; ii < count; ii++)
      for (jj = 0; jj < nInputs; jj++)
         tempW[ii*nInputs+jj] = moatSample[ii][jj];
   printOutTS(PL_INFO,"MOAT genRepair: check for repeated sample points.\n");
   for (ii = 0; ii < count; ii++)
   {
      status = compareSamples(ii,count,nInputs, tempW, states);
      if (status >= 0)
         printOutTS(PL_INFO,
              "MOAT genRepair check: sample %d and %d are identical.\n",
              ii+1,status+1);
   }
   printOutTS(PL_INFO,
        "MOAT genRepair: repair file created in MOAT_repair_file.\n");
   printOutTS(PL_INFO,"         Make sure to change the input indices.\n");
   for (ii = 0; ii < nPaths*(nInputs+1); ii++) delete [] moatSample[ii];
   delete [] moatSample;
   delete [] tempW;
   delete [] states;
   return 0;
}

//*************************************************************************
//* generate sample with response/surface constraints
//*------------------------------------------------------------------------
int MOATSampling::genRepair(PsuadeData *psIO)
{
   int        status, sInd, faFlag, jj, ii, nPaths, currP, nTrials;
   int        outputID, *indSet, count, trial, iInd, ind, ind2, *states;
   int        nInputs, nOutputs, kk, nSamples;
   double     threshL, threshU, **moatSample, *tempW, dtemp, ddata;
   double     currX1, currX2, filterRange, currY1, currY2;
   double     sLo, sHi, Ymax, Ymin, *iUpperB, *iLowerB, *sampleOutputs;
   char       pString[1000];
   FILE       *fp;
   FuncApprox *faPtr;
   pData      pPtr, pLower, pUpper;

   faFlag = 3;
   faPtr = genFAInteractive(psIO, faFlag);
   if (faPtr == NULL) 
   {
      printOutTS(PL_ERROR,"ERROR detected.\n"); 
      return 0;
   }
   faPtr->setOutputLevel(printLevel_);
   psIO->getParameter("ana_outputid", pPtr);
   outputID = pPtr.intData_;
   Ymax = - 1.0e35;
   Ymin =   1.0e35;
   psIO->getParameter("input_ninputs", pPtr);
   nInputs = pPtr.intData_;
   psIO->getParameter("output_noutputs", pPtr);
   nOutputs = pPtr.intData_;
   psIO->getParameter("method_nsamples", pPtr);
   nSamples = pPtr.intData_;
   psIO->getParameter("output_sample", pPtr);
   sampleOutputs  = pPtr.dbleArray_;

   pLower.clean();
   psIO->getParameter("input_lbounds", pLower);
   iLowerB = pLower.dbleArray_;
   pUpper.clean();
   psIO->getParameter("input_ubounds", pUpper);
   iUpperB = pUpper.dbleArray_;
   for (sInd = 0; sInd < nSamples; sInd++)
   {
      if (sampleOutputs[sInd*nOutputs+outputID] > Ymax)
         Ymax = sampleOutputs[sInd*nOutputs+outputID];
      if (sampleOutputs[sInd*nOutputs+outputID] < Ymin)
         Ymin = sampleOutputs[sInd*nOutputs+outputID];
   }
   sprintf(pString,
           "Please enter the lower bound constraint (Ymin=%e) : ",Ymin);
   threshL = getDouble(pString);
   sprintf(pString,
           "Please enter the upper bound constraint (Ymax=%e) : ",Ymax);
   threshU = getDouble(pString);
   if (threshL >= threshU)
   {
      printf("ERROR : lower bound >= upper bound.\n");
      // Cleanup by Bill Oliver
      delete faPtr;
      return 0;
   }
   sprintf(pString,"Please enter the number of paths to search: ");
   nPaths = getInt(1, 1000, pString);
   sprintf(pString, "Please enter P (resolution: try 4-10) : ");
   currP = getInt(4, 10, pString);
   sprintf(pString, "Please enter the number of trials (> 100) : ");
   nTrials = getInt(101, 10000000, pString);
   moatSample = new double*[nPaths*(nInputs+1)];
   for (ii = 0; ii < nPaths*(nInputs+1); ii++)
      moatSample[ii] = new double[nInputs];
   tempW = new double[nInputs];
   indSet = new int[nInputs];
   count = 0;
   filterRange = threshU - threshL;
   for (ii = 0; ii < nPaths; ii++)
   {
      trial = 0; 
      while (trial < nTrials)
      {
         iInd = count;
         trial++;
         for (jj = 0; jj < nInputs; jj++)
         {
            ind = PSUADE_rand() % currP;
            dtemp = ind * (iUpperB[jj] - iLowerB[jj]) / (currP - 1.0);
            moatSample[iInd][jj] = iLowerB[jj] + dtemp;
         }
         dtemp = faPtr->evaluatePoint(moatSample[iInd]);
         if (dtemp < threshL || dtemp > threshU) continue;

         generateRandomIvector(nInputs, indSet);

         for (jj = 0; jj < nInputs; jj++)
         {
            iInd++;
            for (kk = 0; kk < nInputs; kk++)
               moatSample[iInd][kk] = moatSample[iInd-1][kk];

            ind2 = indSet[jj];

            ddata = moatSample[iInd][ind2]; 
            moatSample[iInd][ind2] = iLowerB[ind2];
            currY1 = faPtr->evaluatePoint(moatSample[iInd]);
            moatSample[iInd][ind2] = iUpperB[ind2];
            currY2 = faPtr->evaluatePoint(moatSample[iInd]);
            currX1 = iLowerB[ind2];
            currX2 = iUpperB[ind2];
            moatSample[iInd][ind2] = ddata;

            if (currY2 >= threshU && currY1 >= threshU)
               currX1 = currX2 = 0.0;
            else if (currY2 <= threshL && currY1 <= threshL)
               currX1 = currX2 = 0.0;
            else if (currY2 > currY1)
            {
               if (currY2 <= threshU) currX2 = iUpperB[ind2];
               else
               {
                  sLo = iLowerB[ind2];
                  sHi = iUpperB[ind2];
                  while (PABS((currY2-threshU)/filterRange)>0.0001)
                  {
                     moatSample[iInd][ind2] = 0.5 * (sLo + sHi);
                     currY2 = faPtr->evaluatePoint(moatSample[iInd]);
                     if (currY2 > threshU) sHi = 0.5 * (sLo + sHi);
                     else                  sLo = 0.5 * (sLo + sHi);
                  }
                  currX2 = moatSample[iInd][ind2];
               }
               if (currY1 >= threshL) currX1 = iLowerB[ind2];
               else
               {
                  sLo = iLowerB[ind2];
                  sHi = iUpperB[ind2];
                  while (PABS((currY1-threshL)/filterRange)>0.0001)
                  {
                     moatSample[iInd][ind2] = 0.5 * (sLo + sHi);
                     currY1 = faPtr->evaluatePoint(moatSample[iInd]);
                     if (currY1 < threshL) sLo = 0.5 * (sLo + sHi);
                     else                  sHi = 0.5 * (sLo + sHi);
                  }
                  currX1 = moatSample[iInd][ind2];
               }
            }
            else
            {
               if (currY1 <= threshU) currX1 = iLowerB[ind2];
               else
               {
                  sLo = iLowerB[ind2];
                  sHi = iUpperB[ind2];
                  while (PABS((currY1-threshU)/filterRange)>0.0001)
                  {
                     moatSample[iInd][ind2] = 0.5 * (sLo + sHi);
                     currY1 = faPtr->evaluatePoint(moatSample[iInd]);
                     if (currY1 > threshU) sLo = 0.5 * (sLo + sHi);
                     else                  sHi = 0.5 * (sLo + sHi);
                  }
                  currX1 = moatSample[iInd][ind2];
               }
               if (currY2 >= threshL) currX2 = iUpperB[ind2];
               else
               {
                  sLo = iLowerB[ind2];
                  sHi = iUpperB[ind2];
                  while (PABS((currY2-threshL)/filterRange)>0.0001)
                  {
                     moatSample[iInd][ind2] = 0.5 * (sLo + sHi);
                     currY2 = faPtr->evaluatePoint(moatSample[iInd]);
                     if (currY2 < threshL) sHi = 0.5 * (sLo + sHi);
                     else                  sLo = 0.5 * (sLo + sHi);
                  }
                  currX2 = moatSample[iInd][ind2];
               }
            }
            if (PABS(currX2-currX1)<0.1*(iUpperB[ind2]-iLowerB[ind2])) 
               break;
            moatSample[iInd][ind2] = ddata;
            tempW[ind2] = PABS(currX2 - currX1) / (currP-1.0);
            moatSample[iInd][ind2] += tempW[ind2];
            if (moatSample[iInd][ind2] > iUpperB[ind2])
                 dtemp = threshL - 1.0;
            else dtemp = faPtr->evaluatePoint(moatSample[iInd]);
            if (dtemp < threshL || dtemp > threshU) 
            {
               moatSample[iInd][ind2] -= 2.0 * tempW[ind2];
               if (moatSample[iInd][ind2] < iLowerB[ind2])
                  break;
               dtemp = faPtr->evaluatePoint(moatSample[iInd]);
               if (dtemp < threshL || dtemp > threshU) 
                  break;
            }
         }
         if (jj == nInputs) 
         {
            count += (nInputs + 1);
            if (printLevel_ > 2)
               printOutTS(PL_INFO,"moatgen: path %d (out of %d) found.\n", 
                    ii+1, nPaths);
            break; 
         }
         else
         {
            if (printLevel_ > 2)
               printOutTS(PL_INFO,"Current path fails (%d out of max %d).\n",
                          trial, nTrials); 
         }
      }
      if (trial >= nTrials)
      {
         printOutTS(PL_INFO,
              "MOAT genRepair FAILS to find all possible paths.\n");
         printOutTS(PL_INFO,
              "Suggestion: try a larger P than %d.\n", currP);
      }
   }
   for (ii = 0; ii < nPaths; ii++)
   {
      dtemp = faPtr->evaluatePoint(moatSample[ii]);
      if (dtemp < threshL || dtemp > threshU)
         printOutTS(PL_INFO,
            "MOAT genRepair:sample %d fails final test (%e <? %e <? %e).\n",
            ii, threshL, dtemp, threshU);
   }
   delete [] tempW;
   delete [] indSet;
   if (trial >= nTrials)
   {
      for (ii = 0; ii < nPaths*(nInputs+1); ii++)
         delete [] moatSample[ii];
      delete [] moatSample;
      delete faPtr;
      return 0;
   }
   fp = fopen("MOAT_repair_file", "w");
   // add a check for NULL by Bill Oliver
   if(fp != NULL)
   {
      fprintf(fp, "BEGIN\n");
      fprintf(fp, "%d %d\n", nPaths*(nInputs+1), nInputs);
      for (ii = 0; ii < nInputs; ii++) fprintf(fp, "%d ", ii+1);
      fprintf(fp, "\n");
      for (ii = 0; ii < nPaths*(nInputs+1); ii++)
      {
         for (jj = 0; jj < nInputs; jj++)
            fprintf(fp, "%e ", moatSample[ii][jj]);
         fprintf(fp, "\n");
      }
      fprintf(fp, "END\n");
      fclose(fp);
   }
   count = nPaths * (nInputs + 1); 
   tempW = new double[count*nInputs];
   states = new int[count];
   for (ii = 0; ii < count; ii++) states[ii] = 1;
   for (ii = 0; ii < count; ii++)
      for (jj = 0; jj < nInputs; jj++)
         tempW[ii*nInputs+jj] = moatSample[ii][jj];
   printOutTS(PL_INFO,
        "MOAT genRepair: check for repeated sample points.\n");
   for (ii = 0; ii < count; ii++)
   {
      status = compareSamples(ii,count,nInputs, tempW, states);
      if (status >= 0)
         printf("moatgen check: sample %d and %d are identical.\n",
                ii+1,status+1);
   }
   printOutTS(PL_INFO,
        "MOAT genRepair: repair file created in MOAT_repair_file.\n");
   printOutTS(PL_INFO,
        "                Make sure to change the input indices.\n");
   for (ii = 0; ii < nPaths*(nInputs+1); ii++) delete [] moatSample[ii];
   delete [] moatSample;
   delete [] tempW;
   delete [] states;
   delete faPtr;
   return 0;
}

//*************************************************************************
//* initialize the sampling data
//*------------------------------------------------------------------------
int MOATSampling::initializeHighDimension()
{
   int    nReps, ii, rr, ind1, ind2, kk, *SI, count, *bins;
   int    ii2, nn, currBin;
   double *ranges, *S1, *S2, *tempX;

   nReps = nSamples_ / (nInputs_ + 1);
   allocSampleData();
   ranges = new double[nInputs_];
   for (ii = 0;  ii < nInputs_;  ii++) 
      ranges[ii] = upperBounds_[ii] - lowerBounds_[ii];

   S1 = new double[nInputs_];
   S2 = new double[nInputs_];
   SI = new int[nInputs_];
   for (rr = 0; rr < nReps; rr++)
   {
      for (ii = 0; ii < nInputs_; ii++)
      {
         ind1 = PSUADE_rand() % P_; 
         S1[ii] = ranges[ii] / (P_ - 1) * ind1 + lowerBounds_[ii]; 
         kk = PSUADE_rand() % 2;
         if (kk == 0) ind2 = ind1 - P_ / 2;
         else         ind2 = ind1 + P_ / 2;
         if      (ind2 < 0)    ind2 = ind2 + P_;
         else if (ind2 > P_-1) ind2 = ind2 - P_;
         S2[ii] = ranges[ii] / (P_ - 1) * ind2 + lowerBounds_[ii]; 
      }
      for (ii = 0; ii < nInputs_; ii++)
         sampleMatrix_[rr*(nInputs_+1)][ii] = S1[ii];
      for (ii = 0; ii < nInputs_; ii++) SI[ii] = ii;
      count = nInputs_;
      for (kk = 0; kk < nInputs_; kk++)
      {
         ind1 = PSUADE_rand() % count; 
         ind2 = SI[ind1];
         for (ii = ind1; ii < count-1; ii++) SI[ii] = SI[ii+1];
         count--;
         for (ii = 0; ii < nInputs_; ii++)
            sampleMatrix_[rr*(nInputs_+1)+kk+1][ii] = 
               sampleMatrix_[rr*(nInputs_+1)+kk][ii]; 
         sampleMatrix_[rr*(nInputs_+1)+kk+1][ind2] = S2[ind2]; 
      }
   }

   delete [] ranges;
   delete [] S1;
   delete [] S2;
   delete [] SI;

   if (printLevel_ > 2)
   {
      tempX = new double[2*nReps];
      bins  = new int[P_];
      for (ii = 0; ii < P_; ii++) bins[ii] = 0;
      for (ii = 0; ii < nInputs_; ii++)
      {
         for (rr = 0; rr < nReps; rr++)
         {
            tempX[rr*2] = sampleMatrix_[rr*(nInputs_+1)][ii];
            for (ii2 = 1; ii2 < nInputs_+1; ii2++)
            {
               if (sampleMatrix_[rr*(nInputs_+1)+ii2][ii] != tempX[2*rr])
               {
                  tempX[rr*2+1] = sampleMatrix_[rr*(nInputs_+1)+ii2][ii];
                  break;
               }
            }
         }
         sortDbleList(2*nReps, tempX);
         nn = 1;
         currBin = 0;
         for (ii2 = 1; ii2 < 2*nReps; ii2++)
         {
            if (tempX[ii2] != tempX[ii2-1])
            {
               printOutTS(PL_INFO,
                    "MOAT: input %3d - level = %12.4e, # times = %4d\n",
                    ii+1, tempX[ii2-1], nn);
               bins[currBin++] += nn;
               nn = 1;
            } else nn++;
         }
         printOutTS(PL_INFO,
              "MOAT: input %3d - level = %12.4e, # times = %4d\n",
              ii+1, tempX[ii2-1], nn);
         bins[currBin++] += nn;
      }
      for (ii = 0; ii < P_; ii++) 
         printOutTS(PL_INFO,
              "MOAT: frequency of visit to bin %5d = %d\n",ii+1,bins[ii]);
      delete [] bins;
      delete [] tempX;
   }
   return 0;
}

// ************************************************************************
// equal operator modified by Bill Oliver
// ------------------------------------------------------------------------
MOATSampling& MOATSampling::operator=(const MOATSampling & ms)
{
  if(this == &ms) return *this;
   samplingID_ = ms.samplingID_;
   P_ = ms.P_;
   nInputs_ = ms.nInputs_;
   inputSubset_ = new int[nInputs_];
   for(int i = 0; i < nInputs_; i++)
     inputSubset_[i] = ms.inputSubset_[i];

   printLevel_ = ms.printLevel_;
   samplingID_ = ms.samplingID_;
   nSamples_ = ms.nSamples_;
   nOutputs_ = ms.nOutputs_;
   randomize_ = ms.randomize_;
   nReplications_ = ms.nReplications_;
   lowerBounds_ = new double[nInputs_];
   upperBounds_ = new double[nInputs_];
   for (int i = 0; i < nInputs_; i++)
   {
      lowerBounds_[i] = ms.lowerBounds_[i];
      upperBounds_[i] = ms.upperBounds_[i];
   }
   sampleMatrix_ = new double*[nSamples_];
   for (int i = 0; i < nSamples_; i++)
   {
      sampleMatrix_[i] = new double[nInputs_];
      for(int j = 0; j < nInputs_; j++)
         sampleMatrix_[i][j] = ms.sampleMatrix_[i][j];
   }
   sampleOutput_ = new double[nSamples_*nOutputs_];
   for (int i = 0; i < nSamples_*nOutputs_; i++)
      sampleOutput_[i] = ms.sampleOutput_[i];
   sampleStates_ = new int[nSamples_];
   for (int i = 0; i < nSamples_; i++)
      sampleStates_[i] = ms.sampleStates_[i];
   return (*this);
}

