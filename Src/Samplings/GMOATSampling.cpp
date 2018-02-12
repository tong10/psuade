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
// Functions for the Generalized Morris one-at-a-time class 
// AUTHOR : CHARLES TONG
// DATE   : 2006
// ************************************************************************
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "sysdef.h"
#include "PsuadeUtil.h"
#include "pData.h"
#include "PsuadeData.h"
#include "Psuade.h"
#include "GMOATSampling.h"
#define PABS(x) ((x) > 0 ? (x) : -(x))

//*************************************************************************
//* Constructor
//*------------------------------------------------------------------------
GMOATSampling::GMOATSampling() : Sampling()
{
   samplingID_  = PSUADE_SAMP_GMOAT;
   inpLevels_   = NULL;
   inputSubset_ = NULL;
   P_  = 4;
   initX_ = NULL;
}

//*************************************************************************
//* Copy Constructor added by Bill Oliver
//*------------------------------------------------------------------------
GMOATSampling::GMOATSampling(const GMOATSampling & gms) : Sampling()
{
   int maxReps = 5000;
   int nlevels;
   samplingID_  = gms.samplingID_;
   P_ = gms.P_;
   nInputs_ = gms.nInputs_;
   inpLevels_   = new int[nInputs_];
   for(int i = 0; i < nInputs_; i++)
     inpLevels_[i] = gms.inpLevels_[i];
   
   inputSubset_ = new int[nInputs_];
   for(int i = 0; i < nInputs_; i++)
     inputSubset_[i] = gms.inputSubset_[i];
   
   initX_ = new int*[nInputs_];
   for(int i = 0; i < nInputs_; i++){
     nlevels = inpLevels_[i];
     initX_[i] = new int[maxReps+nlevels];
     for(int j = 0; j < maxReps; j++)
       initX_[i][j] = gms.initX_[i][j];
   }

   // MetisSampling inherits from Sampling so include the parent 
   // class data members
   printLevel_ = gms.printLevel_;
   samplingID_ = gms.samplingID_;
   nSamples_ = gms.nSamples_;
   nOutputs_ = gms.nOutputs_;
   randomize_ = gms.randomize_;
   nReplications_ = gms.nReplications_;
   lowerBounds_ = new double[nInputs_];
   upperBounds_ = new double[nInputs_];
   for (int i = 0; i < nInputs_; i++)
   {
      lowerBounds_[i] = gms.lowerBounds_[i];
      upperBounds_[i] = gms.upperBounds_[i];
   }
   sampleMatrix_ = new double*[nSamples_];
   for (int i = 0; i < nSamples_; i++)
   {
      sampleMatrix_[i] = new double[nInputs_];
      for(int j = 0; j < nInputs_; j++)
         sampleMatrix_[i][j] = gms.sampleMatrix_[i][j];
   }
   sampleOutput_ = new double[nSamples_*nOutputs_];
   for (int i = 0; i < nSamples_*nOutputs_; i++)
      sampleOutput_[i] = gms.sampleOutput_[i];
   sampleStates_ = new int[nSamples_];
   for (int i = 0; i < nSamples_; i++)
      sampleStates_[i] = gms.sampleStates_[i];
}

//*************************************************************************
//* destructor 
//*------------------------------------------------------------------------
GMOATSampling::~GMOATSampling()
{
   if (inpLevels_   != NULL) delete [] inpLevels_;
   if (inputSubset_ != NULL) delete [] inputSubset_;
   if (initX_ != NULL)
   {
      for (int ii = 0; ii < nInputs_; ii++)
         if (initX_[ii] != NULL) delete [] initX_[ii];
      delete [] initX_;
   }
}

//*************************************************************************
//* initialize the sampling data
//*------------------------------------------------------------------------
int GMOATSampling::initialize(int initLevel)
{
   int    ii, ii2, rr, ss, nReps, nn, index, idata, nlevels;
   int    *counts, kk1, kk2, maxSamples;
   int    maxReps=5000, base1, base2, nSub;
   double **BS, *ranges, ddata, delta, rdata, **ranArray, *tempX, maxDist;
   double dDist;
   char   *cString, winput1[500], winput2[500];
   char   partitionFile[500], pString[500], pString2[500];
   FILE   *fp;

   if (nSamples_ == 0)
   {
      printf("GMOATSampling::initialize ERROR - nSamples = 0.\n");
      exit(1);
   }
   if (nInputs_ == 0 || lowerBounds_ == NULL || upperBounds_ == NULL)
   {
      printf("GMOATSampling::initialize ERROR - input not set up.\n");
      exit(1);
   }

   deleteSampleData();
   nReps = nSamples_ / (nInputs_ + 1);
   if ((nReps * (nInputs_+1)) != nSamples_) 
   {
      printf("GMOATSampling : nSamples should be multiples of nInputs+1.\n");
      printf("                nSamples reset to be 10*(nInputs+1).\n");
      nSamples_ = 10 * (nInputs_ + 1);
   }
   if (initLevel != 0) return 0;

   if (printLevel_ > 4)
   {
      printf("GMOATSampling::initialize: nSamples  = %d\n", nSamples_);
      printf("GMOATSampling::initialize: nInputs   = %d\n", nInputs_);
      printf("GMOATSampling::initialize: nOutputs  = %d\n", nOutputs_);
      printf("GMOATSampling: initialize: numLevels = %d\n", P_);
      if (randomize_ != 0)
           printf("GMOATSampling::initialize: randomize on\n");
      else printf("GMOATSampling::initialize: randomize off\n");
      for (ii2 = 0; ii2 < nInputs_; ii2++)
         printf("    GMOATSampling input %3d = [%e %e]\n", ii2+1,
                lowerBounds_[ii2], upperBounds_[ii2]);
   }

   if (nInputs_ > 100)
   {
      printf("GMOATSampling: nInputs > 100, use fast version.\n");
      initializeHighDimension();
      return 0;
   }

   if (nReps > maxReps) maxReps = nReps;
   if (inpLevels_ != NULL) delete [] inpLevels_;
   inpLevels_ = new int[nInputs_];
   for (ii = 0; ii < nInputs_; ii++) inpLevels_[ii] = P_;

   if (psConfig_ != NULL)
   {
      cString = psConfig_->getParameter("GMOAT_P");
      if (cString != NULL)
      {
         sscanf(cString, "%s %s %d",winput1,winput2,&P_);
         P_ = P_ / 2 * 2;
         if (P_ <= 0 || P_ > 100) P_ = 4;
         printf("GMOATSampling: P set to %d (config)\n", P_);
         for (ii = 0; ii < nInputs_; ii++) inpLevels_[ii] = P_;
      }
      cString = psConfig_->getParameter("GMOAT_partition_file");
      if (cString != NULL)
      {
         sscanf(cString, "%s %s %s",winput1,winput2,partitionFile);
         printf("GMOATSampling: use GMOAT input partition file %s.\n",
                partitionFile);
         fp = fopen(partitionFile, "r");
         if (fp != NULL)
         {
            if (inputSubset_ != NULL) delete [] inputSubset_;
            inputSubset_ = NULL;
            fscanf(fp, "%d", &ss);
            if (ss <= 0 || ss >= nInputs_)
            {
               printf("GMOATSampling: invalid GMOAT input partition file.\n");
               printf("              The first line should be nInputs.\n");
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
                     printf("GMOATSampling: invalid input partition file.\n");
                     printf("               invalid input index %d.\n", ii);
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
   if (psSamExpertMode_ == 1)
   {
      printf("GMOATSampling: default number of levels = %d.\n", P_);
      sprintf(pString,
              "Change number of levels for individual inputs ? (y or n) ");
      getString(pString, winput1);
      if (winput1[0] == 'y')
      {
         ii = nInputs_+1;
         sprintf(pString,
                 "Which input to change ? (1 - %d, 0 to end) ",nInputs_);
         while (ii > 0)
         {
            ii = getInt(0, nInputs_, pString);
            if (ii == 0) break;
            sprintf(pString2,
                    "Enter number of levels for input %d (2 - 1000) : ",ii);
            inpLevels_[ii-1] = getInt(2, 1000, pString2);
         }
      }
   }
   if (printLevel_ > 2) 
   {
      for (ii = 0; ii < nInputs_; ii++)
         printf("GMOATSampling: input %3d has %4d levels\n",ii+1,
                inpLevels_[ii]);
   }
 
   if (psSamExpertMode_ == 1)
   {
      printf("GMOAT has an option for an input partition file.\n");
      printf("It is used in case you need to separate the inputs\n");
      printf("into two sets so that for each MOAT path members of\n");
      printf("the selected subset will be varied first.\n");
      printf("This option is generally not needed, but maybe useful\n");
      printf("for parallel computing.\n");
      sprintf(pString,"Do you have a MOAT input partition file ? (y or n) ");
      getString(pString, winput1);
      if (winput1[0] == 'y')
         printf("GMOAT INFO: this option is now available only in configFile.\n");
   }

   if (initX_ != NULL)
   {
      for (ii = 0; ii < nInputs_; ii++)
         if (initX_[ii] != NULL) delete [] initX_[ii];
      delete [] initX_;
   }
   initX_ = new int*[nInputs_];
   for (ii = 0; ii < nInputs_; ii++)
   {
      nlevels = inpLevels_[ii] - 1;
      initX_[ii] = new int[maxReps+nlevels];
      for (ii2 = 0; ii2 < maxReps; ii2+=nlevels)
         generateRandomIvector(nlevels, &(initX_[ii][ii2]));
   }

   maxSamples = maxReps * (nInputs_ + 1);
   if (nSamples_ > maxSamples) maxSamples = nSamples_;
   BS = new double*[maxSamples];
   for (ii = 0; ii < maxSamples; ii++) BS[ii] = new double[nInputs_];

   allocSampleData();
   ranges = new double[nInputs_];
   for (ii = 0;  ii < nInputs_;  ii++) 
      ranges[ii] = upperBounds_[ii] - lowerBounds_[ii];

   for (rr = 0; rr < maxReps; rr++) generate(&BS[rr*(nInputs_+1)],rr);
   for (rr = 1; rr < nReps; rr++)
   {
      printf("GMOATSampling::generate: finding path %d (out of %d)\n",
             rr+1, nReps);
      maxDist = 0.0;
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

   // ----------------------------------------------------------------
   if (printLevel_ > 1)
   {
      counts = new int[nSamples_];
      tempX  = new double[nSamples_];
      for (ii = 0; ii < nInputs_; ii++)
      {
         nlevels = inpLevels_[ii] - 1;
         for (ss = 0; ss <= nlevels; ss++) counts[ss] = 0;
         for (ss = 0; ss < nSamples_; ss+=(nInputs_+1))
         {
            for (ii2 = 0; ii2 < nInputs_+1; ii2++)
               tempX[ii2] = BS[ii2+ss][ii];
            sortDbleList(nInputs_+1, tempX);
            nn = (int) ((tempX[0] + 1.0e-8) * nlevels);
            counts[nn]++;
            nn = (int) ((tempX[nInputs_] + 1.0e-8) * nlevels);
            counts[nn]++;
         }
         for (ii2 = 0; ii2 <= nlevels; ii2++)
            printf("GMOAT: input %3d - level = %2d, # replications = %3d\n",
                   ii+1, ii2, counts[ii2]);
      }
      delete [] tempX;
      delete [] counts;
   }
   
   if (inputSubset_ != NULL)
   {
      nSub = 0;
      for (ii = 0; ii < nInputs_; ii++)
         if (inputSubset_[ii] == 1) nSub++;
   }
   else nSub = nInputs_;

   if (randomize_ == 1)
   {
      printf("You are turning on randomFlag in GMOAT.\n");
      printf("Please confirm (y or n): \n");
      scanf("%s", winput1);
      if (winput1[0] != 'y')
      {
         printf("Turning off randomFlag.\n");
         randomize_ = 0;
      }
   }
   if (randomize_ == 0)
   {
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
   }    
   else
   {
      ranArray = new double*[nInputs_];
      for (ii = 0; ii < nInputs_; ii++)
      {
         ranArray[ii] = new double[inpLevels_[ii]];
         delta = 0.5 / (double) (inpLevels_[ii] - 1);
         for (ii2 = 0; ii2 < inpLevels_[ii]; ii2++)
            ranArray[ii][ii2] = (PSUADE_drand() - 0.5) * delta;
      } 
      for (ss = 0; ss < nSamples_; ss+=(nInputs_+1))
      {
         for (ii = 0; ii < nInputs_; ii++)
         {
            ddata = BS[ss][ii];
            ddata = ddata * ranges[ii] + lowerBounds_[ii];
            sampleMatrix_[ss][ii] = ddata;
         }
         for (ii = 0; ii < nInputs_; ii++)
         {
            for (ii2 = 0; ii2 < nInputs_; ii2++)
            {
               if (BS[ss+ii+1][ii2] != BS[ss+ii][ii2])
               {
                  index = ii2; 
                  break;
               }
            }
            for (ii2 = 0; ii2 < nInputs_; ii2++)
               sampleMatrix_[ss+ii+1][ii2] = sampleMatrix_[ss+ii][ii2];
            ddata = BS[ss+ii+1][index];
            idata = (int) (ddata * (inpLevels_[ii] - 1) + 1.0e-5);
            rdata = ranArray[ii][idata];
            if (ddata > 0.0 && ddata < 1.0) ddata += rdata;
            sampleMatrix_[ss+ii+1][index] = ddata * ranges[index] + 
                                            lowerBounds_[index];
         }
         for (ii = nSub+1; ii <= nInputs_; ii++)
            sampleStates_[ss+ii] = 1;
      }
      for (ii = 0; ii < nInputs_; ii++) delete [] ranArray[ii]; 
      delete [] ranArray;
   }    

   if (repair(NULL, 0) != 0)
   {
      if (checkSample(nInputs_, nSamples_, sampleMatrix_) != 0)
      {
         printf("GMOATSampling : generated sample is not MOAT.\n");
         exit(1);
      }
   }

   if (printLevel_ > 1)
   {
      tempX = new double[nSamples_];
      for (ii = 0; ii < nInputs_; ii++)
      {
         for (ii2 = 0; ii2 < nSamples_; ii2++)
            tempX[ii2] = sampleMatrix_[ii2][ii];
         sortDbleList(nSamples_, tempX);
         nn = 1;
         for (ii2 = 1; ii2 < nSamples_; ii2++)
         {
            if (tempX[ii2] != tempX[ii2-1])
            {
               printf("GMOAT: input %3d - level = %12.4e, # visits = %4d\n",
                      ii+1, tempX[ii2-1], nn);
               nn = 1;
            } else nn++;
         }
         printf("GMOAT: input %3d - level = %12.4e, # visits = %4d\n",
             ii+1, tempX[ii2-1], nn);
      }
      delete [] tempX;
   }
   
   for (ii = 0;  ii < maxSamples; ii++) delete [] BS[ii];
   delete [] BS;
   delete [] ranges;
   return 0;
}

//*************************************************************************
//* generate the BS matrix
//*------------------------------------------------------------------------
int GMOATSampling::generate(double **BS, int index)
{
   int    *permute, ss, ii, ii2, idata, nSub, *subset;
   double **B, *D, *X, **B2, delta;

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

#if 1 
   for (ii = 0; ii < nInputs_; ii++)
   {
      idata = initX_[ii][index];
      X[ii] = ((double) idata) / (double) (inpLevels_[ii] - 1);
   }
#else
   int imax = P_ - 1;
   for (ii = 0; ii < nInputs_; ii++)
   {
      X[ii] = PSUADE_drand();
      idata = (int) (X[ii] * imax);
      if (idata >= imax) idata--;
      X[ii] = (double) idata / (double) (P_ - 1);
   }
#endif

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
         B2[ii][ii2] = B[ii][permute[ii2]];
   for (ii = 0; ii <= nInputs_; ii++)
   {
      for (ii2 = 0; ii2 < nInputs_; ii2++)
      {
         delta = 1.0 / ((double) inpLevels_[ii2] - 1.0);
         BS[ii][ii2] = X[ii2]+delta/2*((2.0*B2[ii][ii2]-1.0)*D[ii2]+1.0);
      }
   }

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
// refine the sample space (not supported)
//-------------------------------------------------------------------------
int GMOATSampling::refine(int, int, double, int, double *)
{
   printf("GMOATSampling: refine not implemented and not needed.\n");
   printf("               You can create a new GMOAT sample and\n");
   printf("               concatenate with the old one.\n");
   return 1;
}

//*************************************************************************
//* repair a MOAT design due to constraints
//  (substitute a subset of inputs with another set of patterns from a
//   repair file) 
//*------------------------------------------------------------------------
int GMOATSampling::repair(char *fname, int start)
{
   int    nPatterns, nInps, nSets, ii, jj, kk, *inpList, pindex, iindex;
   double **patterns, *tempSample;
   char   inStr[500], winput1[500], winput2[500], *cString, repairFile[500];
   FILE   *fp;

   if (start / (nInputs_+1) * (nInputs_ + 1) != start)
   {
      printf("GMOATSampling : start should be multiples of nInputs+1.\n");
      exit(1);
   }

   fp = NULL;
   if (fname != NULL)
   {
      fp = fopen(fname, "r");
      if (fp != NULL)
         printf("GMOAT repair file found = %s.\n", fname);
   }
   if (fp == NULL && psConfig_ != NULL)
   {
      cString = psConfig_->getParameter("GMOAT_repair_file");
      if (cString != NULL)
      {
         sscanf(cString, "%s %s %s",winput1,winput2,repairFile);
         fp = fopen(repairFile, "r");
         if (fp != NULL)
            printf("GMOAT repair file found = %s.\n", repairFile);
      }
   }
   if (fp == NULL)
   {
      printf("GMOAT repair: no repair file found.\n");
      return 1;
   }

   fscanf(fp, "%s", inStr);
   if (strcmp(inStr, "BEGIN"))
   {
      printf("GMOATSampling : wrong format in repair file.\n");
      printf("First  line : BEGIN\n");
      printf("Second line : nPatterns nInputs.\n");
      printf("Third  line : a list of input IDs (1-based).\n");
      printf("Fourth line : (and on) set of patterns.\n");
      printf("Last   line : END\n");
      fclose(fp);
      exit(1);
   }

   fscanf(fp, "%d %d", &nPatterns, &nInps);
   if (nPatterns <= 0 || nInps <= 0)
   {
      printf("GMOATSampling : nPatterns or nInps <= 0.\n");
      fclose(fp);
      exit(1);
   }
   nSets = nPatterns / (nInps + 1);
   if (nSets*(nInps+1) != nPatterns)
   {
      printf("GMOATSampling ERROR: nPatterns must be multiples of nInputs+1\n");
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
         printf("GMOATSampling ERROR: input index out of range (%d,%d)\n",
                inpList[ii], nInputs_);
         fclose(fp);
         exit(1);
      }
      for (jj = 0; jj < ii; jj++)
      {
         if (inpList[ii] == inpList[jj])
         {
            printf("GOATSampling ERROR: repeated index (%d)\n",inpList[ii]);
            fclose(fp);
            exit(1);
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
      printf("GMOATSampling : wrong format in repair file.\n");
      printf("The file should end with END\n");
      exit(1);
   }

   if (checkSample(nInps, nPatterns, patterns) != 0)
   {
      printf("GMOATSampling : pattern in repair file is not MOAT.\n");
      exit(1);
   }

   tempSample = new double[nInps];
   for (ii = start; ii < nSamples_; ii+=nInputs_+1)
   {
      pindex = (ii / (nInputs_+1)) * (nInps + 1);

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
      printf("GMOATSampling : repaired file is not MOAT.\n");
      exit(1);
   }

   delete [] inpList;
   for (ii = 0; ii < nPatterns; ii++) delete [] patterns[ii];
   delete [] patterns;
   delete [] tempSample;
   return 0;
}

//*************************************************************************
//* merge two MOAT designs
//*------------------------------------------------------------------------
int GMOATSampling::merge()
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
      printf("GMOATSampling ERROR: File %s not found.\n", file1);
      return 1;
   }
   else fclose(fp1);
   psuadeIO1 = new PsuadeData();
   psuadeIO1->setOutputLevel(0);
   if (psuadeIO1->readPsuadeFile(file1) != 0)
   {
      printf("GMOAT ERROR: problem with reading file %s.\n", file1);
      delete psuadeIO1;
      return 1;
   }
   psuadeIO1->getParameter("input_ninputs", pPtr1);
   nInps1 = pPtr1.intData_;
   psuadeIO1->getParameter("input_names", pINames1);
   inpNames1 = pINames1.strArray_;
   psuadeIO1->getParameter("output_noutputs", pPtr1);
   nOuts1 = pPtr1.intData_;
   psuadeIO1->getParameter("method_sampling", pPtr1);
   samplingMethod = pPtr1.intData_;
   if (samplingMethod != PSUADE_SAMP_GMOAT)
   {
      printf("GMOAT Merge ERROR: data1 is not GMOAT.\n");
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
      printf("GMOAT Merge ERROR: first sample is not GMOAT.\n");
      delete psuadeIO1;
      pINames1.clean();
      return 1;
   }

   sprintf(pString,"Please enter the name of the second MOAT datafile: ");
   getString(pString, file2);
   file2[strlen(file2)-1] = '\0';
   if ((fp2=fopen(file2,"r")) == NULL)
   {
      printf("GMOAT ERROR : File %s not found.\n", file2);
      delete psuadeIO1;
      pINames1.clean();
      return 1;
   }
   else fclose(fp2);
                                                                  
   psuadeIO2 = new PsuadeData();
   psuadeIO2->setOutputLevel(0);
   if (psuadeIO2->readPsuadeFile(file2) != 0)
   {
      printf("GMOAT ERROR : problem with reading file %s.\n", file2);
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
   if (samplingMethod != PSUADE_SAMP_GMOAT)
   {
      printf("GMOAT Merge ERROR : data2 is not GMOAT.\n");
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
      printf("GMOAT Merge ERROR : second sample is not GMOAT.\n");
      delete psuadeIO1;
      delete psuadeIO2;
      pINames1.clean();
      pINames2.clean();
      return 1;
   }

   nReps = nSamp1 / (nInps1 + 1);
   if (nReps != (nSamp2 / (nInps2 + 1)))
   {
      printf("GMOAT Merge ERROR : different number of replications.\n");
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
int GMOATSampling::checkSample(int nInputs, int nSamples, double **X)
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
      if (nDiff != nInputs){
	// Clean up by Bill Oliver
	delete [] errArray;
	return 1;
      }
   }
   delete [] errArray;
   return 0;
}

//*************************************************************************
//* check whether the sample is MOAT (return 0 if yes)
//  (input X is in a different format)
//*------------------------------------------------------------------------
int GMOATSampling::checkSample2(int nInputs, int nSamples, double *X)
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
      if (nDiff != nInputs){
	// Cleanup by Bill Oliver
	delete [] errArray;
	return 1;
      }
   }
   delete [] errArray;
   return 0;
}

//*************************************************************************
//* initialize the sampling data
//*------------------------------------------------------------------------
int GMOATSampling::initializeHighDimension()
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
         if (kk == 0) ind2 = ind1 - 1;
         else         ind2 = ind1 + 1;
         if      (ind2 < 0)    ind2 = ind2 + 2;
         else if (ind2 > P_-1) ind2 = ind2 - 2;
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
               printf("GMOAT: input %3d - level = %12.4e, # times = %4d\n",
                      ii+1, tempX[ii2-1], nn);
               bins[currBin++] += nn;
               nn = 1;
            } else nn++;
         }
         printf("GMOAT: input %3d - level = %12.4e, # times = %4d\n",
                ii+1, tempX[ii2-1], nn);
         bins[currBin++] += nn;
      }
      for (ii = 0; ii < P_; ii++) 
         printf("GMOAT: frequency of visit to bin %5d = %d\n",ii+1,bins[ii]);
      delete [] bins;
      delete [] tempX;
   }
   return 0;
}

// ************************************************************************
// equal operator modified by Bill Oliver
// ------------------------------------------------------------------------
GMOATSampling& GMOATSampling::operator=(const GMOATSampling &gms)
{
   if(this == &gms) return *this;

   int maxReps = 5000;
   int nlevels;
   samplingID_  = gms.samplingID_;
   P_ = gms.P_;
   nInputs_ = gms.nInputs_;
   inpLevels_   = new int[nInputs_];
   for(int i = 0; i < nInputs_; i++)
      inpLevels_[i] = gms.inpLevels_[i];
   
   inputSubset_ = new int[nInputs_];
   for(int i = 0; i < nInputs_; i++)
      inputSubset_[i] = gms.inputSubset_[i];
   
   initX_ = new int*[nInputs_];
   for (int i = 0; i < nInputs_; i++)
   {
      nlevels = inpLevels_[i];
      initX_[i] = new int[maxReps+nlevels];
      for(int j = 0; j < maxReps; j++)
         initX_[i][j] = gms.initX_[i][j];
   }

   // MetisSampling inherits from Sampling so include the parent 
   // class data members
   printLevel_ = gms.printLevel_;
   samplingID_ = gms.samplingID_;
   nSamples_ = gms.nSamples_;
   nOutputs_ = gms.nOutputs_;
   randomize_ = gms.randomize_;
   nReplications_ = gms.nReplications_;
   lowerBounds_ = new double[nInputs_];
   upperBounds_ = new double[nInputs_];
   for (int i = 0; i < nInputs_; i++)
   {
      lowerBounds_[i] = gms.lowerBounds_[i];
      upperBounds_[i] = gms.upperBounds_[i];
   }
   sampleMatrix_ = new double*[nSamples_];
   for (int i = 0; i < nSamples_; i++)
   {
      sampleMatrix_[i] = new double[nInputs_];
      for(int j = 0; j < nInputs_; j++)
         sampleMatrix_[i][j] = gms.sampleMatrix_[i][j];
   }
   sampleOutput_ = new double[nSamples_*nOutputs_];
   for (int i = 0; i < nSamples_*nOutputs_; i++)
      sampleOutput_[i] = gms.sampleOutput_[i];
   sampleStates_ = new int[nSamples_];
   for (int i = 0; i < nSamples_; i++)
      sampleStates_[i] = gms.sampleStates_[i];
   return (*this);
}

