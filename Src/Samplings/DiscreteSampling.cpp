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
// Functions for the discrete sampling class 
// AUTHOR : CHARLES TONG
// DATE   : 2010
// ************************************************************************
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "sysdef.h"
#include "PsuadeUtil.h"
#include "FunctionInterface.h"
#include "DiscreteSampling.h"

#define PABS(x) ((x) > 0 ? x : -(x))
// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
DiscreteSampling::DiscreteSampling() : Sampling()
{
   samplingID_   = PSUADE_SAMP_DISCRETE;
   inputValCnts_ = NULL;
   inputValues_  = NULL;
   inputProbs_   = NULL;
}

// ************************************************************************
// copy constructor by Bill Oliver
// ------------------------------------------------------------------------
DiscreteSampling::DiscreteSampling(const DiscreteSampling & ds) : Sampling()
{
   printLevel_ = ds.printLevel_;
   samplingID_ = ds.samplingID_;
   nSamples_ = ds.nSamples_;
   nInputs_ = ds.nInputs_;
   nOutputs_ = ds.nOutputs_;
   randomize_ = ds.randomize_;
   nReplications_ = ds.nReplications_;
   lowerBounds_ = new double[nInputs_];
   upperBounds_ = new double[nInputs_];
   for (int ii = 0; ii < nInputs_; ii++)
   {
      lowerBounds_[ii] = ds.lowerBounds_[ii];
      upperBounds_[ii] = ds.upperBounds_[ii];
   }
   sampleMatrix_ = new double*[nSamples_];
   for (int ii = 0; ii < nSamples_; ii++)
   {
      sampleMatrix_[ii] = new double[nInputs_];
      for(int jj = 0; jj < nInputs_; jj++)
         sampleMatrix_[ii][jj] = ds.sampleMatrix_[ii][jj];
   }
   sampleOutput_ = new double[nSamples_*nOutputs_];
   for (int ii = 0; ii < nSamples_*nOutputs_; ii++)
      sampleOutput_[ii] = ds.sampleOutput_[ii];
   sampleStates_ = new int[nSamples_];
   for (int ii = 0; ii < nSamples_; ii++)
      sampleStates_[ii] = ds.sampleStates_[ii];
   inputValCnts_ = NULL;
   inputValues_  = NULL;
   inputProbs_   = NULL;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
DiscreteSampling::~DiscreteSampling()
{
   if (inputValues_ != NULL)
   {
      for (int ii = 0; ii < nInputs_; ii++)
         if (inputValues_[ii] != NULL) delete [] inputValues_[ii];
      delete [] inputValues_;
   }
   if (inputValCnts_ != NULL) delete [] inputValCnts_;
   if (inputProbs_ != NULL)
   {
      for (int ii = 0; ii < nInputs_; ii++)
         if (inputProbs_[ii] != NULL) delete [] inputProbs_[ii];
      delete [] inputProbs_;
   }
   if (upperBounds_  != NULL) delete [] upperBounds_;
   if (lowerBounds_  != NULL) delete [] lowerBounds_;
   if (sampleMatrix_ != NULL)
   {
      for (int ii = 0; ii < nSamples_; ii++) delete [] sampleMatrix_[ii];
      delete [] sampleMatrix_;
   }
   if (sampleOutput_ != NULL) delete [] sampleOutput_;
   if (sampleStates_ != NULL) delete [] sampleStates_;
}

// ************************************************************************
// initialize the sampler data
// ------------------------------------------------------------------------
int DiscreteSampling::initialize(int initLevel)
{
   int    inLeng, ii, jj, kk, ll, sampleID, **probArrays, cnt, num;
   int    **inputValues_=NULL, *inputValCnts_=NULL;
   double dsum, **inputProbs_ = NULL;
   char   lineIn[1001], filename[1001], cword[1001];
   FILE   *fp=NULL;

   deleteSampleData();
   if (inputValues_ != NULL)
   {
      for (ii = 0; ii < nInputs_; ii++)
         if (inputValues_[ii] != NULL) delete [] inputValues_[ii];
      delete [] inputValues_;
   }
   if (inputProbs_ != NULL)
   {
      for (ii = 0; ii < nInputs_; ii++)
         if (inputProbs_[ii] != NULL) delete [] inputProbs_[ii];
      delete [] inputProbs_;
   }
   if (inputValCnts_ != NULL) delete [] inputValCnts_;
   inputValues_  = NULL;
   inputValCnts_ = NULL;
   inputProbs_   = NULL;

   printf("To use discrete sampling, you will have to provide a\n");
   printf("parameter file. The parameter file should be of the\n");
   printf("following format:\n");
   printf("PSUADE_BEGIN\n");
   printf("<number of inputs> <number of sample points>\n");
   printf("1 numLevels(n) L1 P1 ... Ln Pn <for input 1>\n");
   printf("2 numLevels(n) L1 P1 ... Ln Pn <for input 2>\n");
   printf("...\n");
   printf("PSUADE_END\n");
   printf("Note: numLevels - number of discrete values for the input.\n");
   printf("      L1        - value of the first level\n");
   printf("      P1        - probability of the first level\n");
   printf("      Ln        - value of the last level\n");
   printf("      Pn        - probability of the last level\n");
   printf("Note: the sum of P's for an input should be 1.\n");
   printf("Enter the name of the parameter file: ");
   scanf("%s", filename);
   fgets(lineIn, 500, stdin);
   inLeng = strlen(filename);
   if (inLeng < 500)
   {
      filename[inLeng] = '\0';
      fp = fopen(filename, "r");
      if (fp == NULL)
      {
         printf("DiscreteSampling ERROR: cannot open file %s.\n",filename);
         return -1;
      }
   }
   else
   {
      printf("DiscreteSampling ERROR: file name too long.\n");
      return -1;
   }
   fgets(lineIn, 1000, stdin);
   sscanf(lineIn, "%s", cword);
   if (!strcmp(cword, "PSUADE_BEGIN"))
   {
      scanf("%d %d", &nInputs_, &nSamples_);
      if (nInputs_ <= 0)
      {
         printf("DiscreteSampling ERROR : nInputs <= 0.\n");
         fclose(fp);
         return -1;
      }
      if (nSamples_ <= 0)
      {
         printf("DiscreteSampling ERROR : nSamples <= 0.\n");
         fclose(fp);
         return -1;
      }
      inputValCnts_ = new int[nInputs_];
      inputValues_ = new int*[nInputs_];
      inputProbs_  = new double*[nInputs_];
      for (ii = 0; ii < nInputs_; ii++)
      {
         inputValues_[ii] = NULL;
         inputProbs_[ii]  = NULL;
         inputValCnts_[ii] = 0;
      }
      for (ii = 0; ii < nInputs_; ii++)
      {
         scanf("%d", &kk);
         if (kk != ii+1)
         {
            printf("DiscreteSampling ERROR: invalid input index %d (!= %d).\n",
                   kk,ii+1);
            fclose(fp);
            if (inputValues_ != NULL)
            {
               for (int ii = 0; ii < nInputs_; ii++)
                  if (inputValues_[ii] != NULL) delete [] inputValues_[ii];
                     delete [] inputValues_;
            }
            if (inputProbs_ != NULL)
            {
               for (int ii = 0; ii < nInputs_; ii++)
                  if (inputProbs_[ii] != NULL) delete [] inputProbs_[ii];
                     delete [] inputProbs_;
            }
            if (inputValCnts_ != NULL) delete [] inputValCnts_;
            inputValues_  = NULL;
            inputValCnts_ = NULL;
            inputProbs_   = NULL;
            return -1;
         }
         scanf("%d", &num);
         if (num <= 0)
         {
            printf("DiscreteSampling ERROR: invalid numLevels %d (input %d).\n",
                   num,ii+1);
            fclose(fp);
            for(int ii = 0; ii < nInputs_; ii++) delete [] inputValues_[ii];
	    delete [] inputValues_;
            for(int ii = 0; ii < nInputs_; ii++) delete [] inputProbs_[ii];
	    delete [] inputProbs_;
	    delete [] inputValCnts_;
	    inputValues_ = NULL;
	    inputValCnts_ = NULL;
	    inputProbs_ = NULL;
            return -1;
         }
         inputValCnts_[ii] = num;
         inputValues_[ii] = new int[num];
         inputProbs_[ii] = new double[num];
         dsum = 0.0;
         if (printLevel_ > 1) printf("Input %4d: \n", ii+1);
         for (jj = 0; jj < num; jj++)
         {
            scanf("%d %lg",&inputValues_[ii][jj],&inputProbs_[ii][jj]);
            if (inputProbs_[ii][jj] <= 0.0)
            {
               fclose(fp);
               printf("DiscreteSampling ERROR: probability should be > 0.\n");
	       for(int ii = 0; ii < nInputs_; ii++) delete [] inputValues_[ii];
	       delete [] inputValues_;
               for(int ii = 0; ii < nInputs_; ii++) delete [] inputProbs_[ii];
	       delete [] inputProbs_;
	       delete [] inputValCnts_;
	       inputValues_ = NULL;
	       inputValCnts_ = NULL;
	       inputProbs_ = NULL;
               return -1;
            }
            if (printLevel_ > 1)
               printf("     Level = %8d, Probability = %e\n", 
                      inputValues_[ii][jj], inputProbs_[ii][jj]);
            dsum += inputProbs_[ii][jj];
         } 
         if (dsum != 1.0)
         {
            fclose(fp);
            printf("DiscreteSampling ERROR: sum of probability should be 1.\n");
	    for(int ii = 0; ii < nInputs_; ii++) delete [] inputValues_[ii];
	    delete [] inputValues_;
            for(int ii = 0; ii < nInputs_; ii++) delete [] inputProbs_[ii];
	    delete [] inputProbs_;
	    delete [] inputValCnts_;
	    inputValues_ = NULL;
	    inputValCnts_ = NULL;
	    inputProbs_ = NULL;
            return -1;
         }
      }
   }
   else
   {
      printf("DiscreteSampling ERROR: PSUADE_BEGIN not found.\n");
      fclose(fp);
      return -1;
   }
   fgets(lineIn, 1000, stdin);
   fgets(lineIn, 1000, stdin);
   sscanf(lineIn, "%s", cword);
   if (strcmp(cword, "PSUADE_END"))
   {
      printf("DiscreteSampling ERROR: PSUADE_END not found.\n");
      fclose(fp);
      for(int ii = 0; ii < nInputs_; ii++) delete [] inputValues_[ii];
      delete [] inputValues_;
      for(int ii = 0; ii < nInputs_; ii++) delete [] inputProbs_[ii];
      delete [] inputProbs_;
      delete [] inputValCnts_;
      inputValues_ = NULL;
      inputValCnts_ = NULL;
      inputProbs_ = NULL;
      return -1;
   }
   fclose(fp);

   allocSampleData();
   probArrays = new int*[nInputs_];
   for (ii = 0; ii < nInputs_; ii++)
   {
      probArrays[ii] = new int[101];
      cnt = 0;
      kk = inputValCnts_[ii];
      for (jj = 0; jj < kk; jj++)
      {
         num = (int) ((inputProbs_[ii][jj] + 1.0e-12) * 100);
         for (ll = 0; ll < num; ll++)
            probArrays[ii][cnt+ll] = inputValues_[ii][jj];
         cnt += num;
      }
   }
   for (sampleID = 0; sampleID < nSamples_; sampleID++)
   {
      for (ii = 0; ii < nInputs_; ii++)
      {
         kk = PSUADE_rand() % 100;
         sampleMatrix_[sampleID][ii] = (double) probArrays[ii][kk];
      }
   }

   if (printLevel_ > 4)
   {
      printf("DiscreteSampling::initialize: nSamples = %d\n", nSamples_);
      printf("DiscreteSampling::initialize: nInputs  = %d\n", nInputs_);
   }

   for (ii = 0; ii < nInputs_; ii++) delete [] probArrays[ii];
   delete [] probArrays;
   for(int ii = 0; ii < nInputs_; ii++) delete [] inputValues_[ii];
   delete [] inputValues_;
   for(int ii = 0; ii < nInputs_; ii++) delete [] inputProbs_[ii];
   delete [] inputProbs_;
   delete [] inputValCnts_;
   inputValues_ = NULL;
   inputValCnts_ = NULL;
   inputProbs_ = NULL;
   return 0;
}

// ************************************************************************
// refine the sample space
// ------------------------------------------------------------------------
int DiscreteSampling::refine(int refineRatio,int randomize,double thresh,
                              int nSamples, double *sampleErrors)
{
   (void) refineRatio;
   (void) randomize;
   (void) thresh;
   (void) nSamples;
   (void) sampleErrors;
   printf("DiscreteSampling refine: not implemented yet.\n");
   return -1;
} 

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
DiscreteSampling& DiscreteSampling::operator=(const DiscreteSampling & ds)
{
     if(this == & ds) return *this;

   printLevel_ = ds.printLevel_;
   samplingID_ = ds.samplingID_;
   nSamples_ = ds.nSamples_;
   nInputs_ = ds.nInputs_;
   nOutputs_ = ds.nOutputs_;
   randomize_ = ds.randomize_;
   nReplications_ = ds.nReplications_;
   lowerBounds_ = new double[nInputs_];
   upperBounds_ = new double[nInputs_];
   for (int i = 0; i < nInputs_; i++)
   {
      lowerBounds_[i] = ds.lowerBounds_[i];
      upperBounds_[i] = ds.upperBounds_[i];
   }
   sampleMatrix_ = new double*[nSamples_];
   for (int i = 0; i < nSamples_; i++)
   {
      sampleMatrix_[i] = new double[nInputs_];
      for(int j = 0; j < nInputs_; j++)
         sampleMatrix_[i][j] = ds.sampleMatrix_[i][j];
   }
   sampleOutput_ = new double[nSamples_*nOutputs_];
   for (int i = 0; i < nSamples_*nOutputs_; i++)
      sampleOutput_[i] = ds.sampleOutput_[i];
   sampleStates_ = new int[nSamples_];
   for (int i = 0; i < nSamples_; i++)
      sampleStates_[i] = ds.sampleStates_[i];

   return *this;
}

