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
// Functions for the User-specified Metis sampling class 
// (arbitrary domain - determined by UserMetisDriver) 
// AUTHOR : CHARLES TONG
// DATE   : 2006
// ************************************************************************
#include <assert.h>
#include <stdio.h>
#include <unistd.h>

#include "sysdef.h"
#include "PsuadeUtil.h"
#include "FunctionInterface.h"
#include "UserMetisSampling.h"

#define PABS(x) ((x) > 0 ? x : -(x))

#ifdef HAVE_METIS

extern "C" 
{
void METIS_PartGraphRecursive(int *, int *, int *, int *, int *,
                              int *, int *, int *, int *, int *, int *);
}

#endif

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
UserMetisSampling::UserMetisSampling() : Sampling()
{
   samplingID_ = PSUADE_SAMP_UMETIS;
   n1d_        = -1;
   nAggrs_     = 0;
   aggrCnts_   = NULL;
   aggrLabels_ = NULL;
   cellsOccupied_ = NULL;
}

// ************************************************************************
// copy constructor added by Bill Oliver 
// ------------------------------------------------------------------------
UserMetisSampling::UserMetisSampling(const UserMetisSampling & ums) : Sampling()
{
   int graphN = 1;
   n1d_ = ums.n1d_;
   nAggrs_ = ums.nAggrs_;
  
   if(nAggrs_ > 0) 
   {
      aggrCnts_ = new int[nAggrs_];
      for( int i = 0; i < nAggrs_; i++) aggrCnts_[i] = ums.aggrCnts_[i];
  
      aggrLabels_ = new int*[nAggrs_];
      for(int i = 0; i < nAggrs_; i++)
      {
         aggrLabels_[i] = new int[aggrCnts_[i]];
         for(int j = 0; j < aggrCnts_[i]; j++)
	    aggrLabels_[i][j] = ums.aggrLabels_[i][j];
      }
   } 
   else 
   {
      initialize(0);
   }
   nInputs_ = ums.nInputs_;

   for(int i = 0; i < nInputs_; i++) graphN *= n1d_;

   cellsOccupied_ = new int[graphN];
   for(int i = 0; i < graphN; i++)
      cellsOccupied_[i] = ums.cellsOccupied_[i];
  
   // MetisSampling inherits from Sampling so include the parent 
   // class data members
   printLevel_ = ums.printLevel_;
   samplingID_ = ums.samplingID_;
   nSamples_ = ums.nSamples_;
   nOutputs_ = ums.nOutputs_;
   randomize_ = ums.randomize_;
   nReplications_ = ums.nReplications_;
   lowerBounds_ = new double[nInputs_];
   upperBounds_ = new double[nInputs_];
   for (int i = 0; i < nInputs_; i++)
   {
      lowerBounds_[i] = ums.lowerBounds_[i];
      upperBounds_[i] = ums.upperBounds_[i];
   }
   sampleMatrix_ = new double*[nSamples_];
   for (int i = 0; i < nSamples_; i++)
   {
      sampleMatrix_[i] = new double[nInputs_];
      for(int j = 0; j < nInputs_; j++)
         sampleMatrix_[i][j] = ums.sampleMatrix_[i][j];
   }
   sampleOutput_ = new double[nSamples_*nOutputs_];
   for (int i = 0; i < nSamples_*nOutputs_; i++)
      sampleOutput_[i] = ums.sampleOutput_[i];
   sampleStates_ = new int[nSamples_];
   for (int i = 0; i < nSamples_; i++)
      sampleStates_[i] = ums.sampleStates_[i];
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
UserMetisSampling::~UserMetisSampling()
{
   if (aggrCnts_ != NULL && nAggrs_ > 0) delete [] aggrCnts_;
   if (aggrLabels_ != NULL && nAggrs_ > 0)
   {
      for (int ii = 0; ii < nAggrs_; ii++)
         if (aggrLabels_[ii] != NULL) delete [] aggrLabels_[ii];
      delete [] aggrLabels_;
   }
   if (cellsOccupied_ != NULL) delete [] cellsOccupied_;
}

// ************************************************************************
// initialize the sampler data
// ------------------------------------------------------------------------
int UserMetisSampling::initialize(int initLevel)
{
   int    *incrs, inputID, ii, jj, kk, itmp, jtmp, nnz, sampleID, randFlag;
   int    options[10], index, haveExec;
   int    nAppFiles=5, newGraphN, graphN, *graphI, *graphJ, *aggrMap;
   double *ranges=NULL, dtmp, *tempSample, *tempOutput, *coeffs = NULL; 
   double constraints[2], *lbounds, *ubounds;
   char   **inputNames, **outputNames, **driverNames, command[500];
   char   pString[500];
   FILE   *fp;
   FunctionInterface *funcIO;
                                                                                
   if (nSamples_ == 0)
   {
      printf("UserMetisSampling::initialize ERROR - nSamples = 0.\n");
      exit(1);
   }
   if (nInputs_ == 0 || lowerBounds_ == NULL || upperBounds_ == NULL)
   {
      printf("UserMetisSampling::initialize ERROR - input not set up.\n");
      exit(1);
   }

   haveExec = 1;
   fp = fopen("UserMetisDriver", "r");
   if (fp == NULL)
   {
      haveExec = 0;
      printf("UserMetisSampling missing : \n");
      sprintf(pString,
              "Is it a linear regression with simple bounds? (y or n) ");
      getString(pString, command);
      if (command[0] == 'y')
      {
         coeffs = new double[nInputs_+1];
         sprintf(pString,"Enter the constant coefficient : ");
         coeffs[0] = getDouble(pString);
         for (ii = 1; ii <= nInputs_; ii++)
         {
            sprintf(pString,"Enter the coefficient for input %d : ",ii);
            coeffs[ii] = getDouble(pString);
         }
         sprintf(pString,"Enter the lower bound : ");
         constraints[0] = getDouble(pString);
         sprintf(pString,"Enter the upper bound : ");
         constraints[1] = getDouble(pString);
      }
      else
      {
         printf("UserMetisSampling ERROR : missing program UserMetisDriver.\n");
         printf("      UserMetisDriver should be an executable taking\n");
         printf("      the first argument as input file and second argument\n");
         printf("      as output file. The input file has the first line\n");
         printf("      the number of inputs, followed by all the inputs.\n");
         printf("        The output file just contains all outputs.\n");
         exit(1);
      }
   }
   else fclose(fp);

   deleteSampleData();
   if (aggrCnts_ != NULL && nAggrs_ > 0) delete [] aggrCnts_;
   if (aggrLabels_ != NULL && nAggrs_ > 0)
   {
      for (ii = 0; ii < nAggrs_; ii++)
         if (aggrLabels_[ii] != NULL) delete [] aggrLabels_[ii];
      delete [] aggrLabels_;
   }
   if (cellsOccupied_ != NULL) delete [] cellsOccupied_;

   randFlag = 0;
   if (randomize_ & 1) randFlag = 1;
   if (nInputs_ > 12)
   {
      printf("UserMetisSampling ERROR : does not support nInputs > 12.\n");
      exit(1);
   }
   if (nInputs_ == 1 ) n1d_ = nSamples_;
   if (nInputs_ == 2 ) n1d_ = 512;
   if (nInputs_ == 3 ) n1d_ = 64;
   if (nInputs_ == 4 ) n1d_ = 23;
   if (nInputs_ == 5 ) n1d_ = 12;
   if (nInputs_ == 6 ) n1d_ = 8;
   if (nInputs_ == 7 ) n1d_ = 6;
   if (nInputs_ == 8 ) n1d_ = 5;
   if (nInputs_ == 9 ) n1d_ = 4;
   if (nInputs_ == 10) n1d_ = 3;
   if (nInputs_ >= 11) n1d_ = 3;

   ranges  = new double[nInputs_];
   lbounds = new double[nInputs_];
   ubounds = new double[nInputs_];
   if (printLevel_ > 1)
      printf("UserMetisSampling:: ranges expanded by 10 %% on each side.\n");
   for (inputID = 0; inputID < nInputs_; inputID++)
   {
      ranges[inputID]  = (upperBounds_[inputID] - lowerBounds_[inputID])*0.1;
      lbounds[inputID] = lowerBounds_[inputID] - ranges[inputID];
      ubounds[inputID] = upperBounds_[inputID] + ranges[inputID];
   }
   for (inputID = 0; inputID < nInputs_; inputID++)
      ranges[inputID]  = ubounds[inputID] - lbounds[inputID];

   incrs    = new int[nInputs_+1];
   graphN   = 1;
   incrs[0] = graphN;
   for (inputID = 1; inputID <= nInputs_; inputID++)
   {
      graphN *= n1d_;
      incrs[inputID] = graphN;
   }
   if (nSamples_ > 2*graphN)
   {
      printf("UserMetisSampling ERROR : nSamples %d too large.\n",
             nSamples_);
      exit(1);
   }

   if (haveExec == 1)
   {
      funcIO = new FunctionInterface();
      inputNames = new char*[nInputs_];
      for (ii = 0; ii < nInputs_; ii++)
      {
         inputNames[ii] = new char[200];
         sprintf(inputNames[ii], "X%d", ii+1);
      }
      outputNames = new char*[nInputs_];
      for (ii = 0; ii < nOutputs_; ii++)
      {
         outputNames[ii] = new char[200];
         sprintf(outputNames[ii], "Y%d", ii+1);
      }
      driverNames = new char*[5];
      for (ii = 0; ii < 5; ii++)
      {
         driverNames[ii] = new char[200];
         strcpy(driverNames[ii], "NONE");
      }
      strcpy(driverNames[0], "UserMetisDriver");
      funcIO->loadInputData(nInputs_, inputNames);
      funcIO->loadOutputData(nOutputs_,outputNames);
      funcIO->loadFunctionData(nAppFiles, driverNames);
      for (ii = 0; ii < nInputs_; ii++) delete [] inputNames[ii];
      delete [] inputNames;
      for (ii = 0; ii < nOutputs_; ii++) delete [] outputNames[ii];
      delete [] outputNames;
      for (ii = 0; ii < 5; ii++) delete [] driverNames[ii];
      delete [] driverNames;
   }

   cellsOccupied_ = new int[graphN];
   tempSample = new double[nInputs_];
   tempOutput = new double[nOutputs_];
   newGraphN  = 0;
   for (ii = 0; ii < graphN; ii++)
   {
      itmp = ii;
      for (inputID = 0; inputID < nInputs_; inputID++)
      {
         jtmp = itmp % n1d_;
         itmp = itmp / n1d_;
         dtmp = ((double) (jtmp + 0.5)) / (double) n1d_;
         tempSample[inputID] = dtmp * ranges[inputID] +
                               lowerBounds_[inputID];
      }
      if (ii % 100 == 0) 
         printf("UserMetisSampling::initialize - running %d out of %d.\n",
                ii+1, graphN);
      if (haveExec == 1)
         funcIO->evaluate(1,nInputs_,tempSample,nOutputs_,tempOutput,0);
      else
      {
         for (kk = 0; kk < nOutputs_; kk++)
         {
            dtmp = coeffs[0];
            for (jj = 0; jj < nOutputs_; jj++)
               dtmp += coeffs[jj+1]*tempSample[jj];
            if (dtmp < constraints[0] || dtmp > constraints[1])
                 tempOutput[kk] = 0.0;
            else tempOutput[kk] = 1.0;
         }
      }
      dtmp = 0.0;
      for (jj = 0; jj < nOutputs_; jj++) dtmp += tempOutput[jj];
      if (dtmp == 0.0) cellsOccupied_[ii] = newGraphN++;
      else             cellsOccupied_[ii] = 9999999;
   }
   delete [] tempOutput;
   delete [] tempSample;
   if (haveExec == 1) delete funcIO;
   if (coeffs != NULL) delete [] coeffs;
 
   graphI = new int[newGraphN+1];
   graphJ = new int[newGraphN*nInputs_*2];
   nnz = 0;
   graphI[0] = nnz;
   newGraphN = 0;
   for (ii = 0; ii < graphN; ii++)
   {
      if (cellsOccupied_[ii] != 9999999)
      {
         itmp = ii;
         for (inputID = 0; inputID < nInputs_; inputID++)
         {
            jtmp = itmp % n1d_;
            itmp = itmp / n1d_;
            if (jtmp > 0)
            {
               jj = cellsOccupied_[ii-incrs[inputID]];
               if (jj != 9999999) graphJ[nnz++] = jj;
            }
            if (jtmp < n1d_-1)
            {
               jj = cellsOccupied_[ii+incrs[inputID]];
               if (jj != 9999999) graphJ[nnz++] = jj;
            }
         }
         newGraphN++;
         graphI[newGraphN] = nnz;
      }
   }
   delete [] incrs;

   options[0] = 0;
   aggrMap = new int[newGraphN];
   if (printLevel_ > 1)
      printf("UserMetisSampling:: calling domain partitioner.\n");
#ifdef HAVE_METIS
   int wgtflag=0, numflag=0, edgeCut=0;
   METIS_PartGraphRecursive(&newGraphN, graphI, graphJ, NULL, NULL,
             &wgtflag,&numflag,&nSamples_,options,&edgeCut,aggrMap);
#else
   printf("UserMetisSampling ERROR : METIS not installed.\n");
#endif
   if (printLevel_ > 1)
      printf("UserMetisSampling:: subdomains created.\n");

   nAggrs_   = nSamples_;
   // Bill Oliver added range check for nAggrs_ for defensive programming
   if(nAggrs_ <= 0)
   {
      printf("nAggrs_ is <= 0 in file %s line %d\n", __FILE__, __LINE__);
      exit(1);
   }
   aggrCnts_ = new int[nAggrs_];
   for (ii = 0; ii < nAggrs_; ii++) aggrCnts_[ii] = 0;
   for (ii = 0; ii < graphN; ii++)
   {
      if (cellsOccupied_[ii] != 9999999)
      {
         index = aggrMap[cellsOccupied_[ii]];
         aggrCnts_[index]++;  
      }
   }
   
   aggrLabels_ = new int*[nAggrs_];
   for (ii = 0; ii < nAggrs_; ii++)
   {
      aggrLabels_[ii] = new int[aggrCnts_[ii]];
      aggrCnts_[ii] = 0;
   }
   for (ii = 0; ii < graphN; ii++)
   {
      if (cellsOccupied_[ii] != 9999999)
      {
         index = aggrMap[cellsOccupied_[ii]];
         aggrLabels_[index][aggrCnts_[index]++] = ii;  
      }
   }
   if (initLevel != 0)
   {
      delete [] ubounds;
      delete [] ranges;
      delete [] lbounds;
      delete [] aggrMap;
      return 0;
   }

   allocSampleData();

   for (sampleID = 0; sampleID < nSamples_; sampleID++)
   {
      if (randFlag == 1)
           index = (int) (PSUADE_drand() * aggrCnts_[sampleID]);
      else index = aggrCnts_[sampleID] / 2;
      if (index == aggrCnts_[sampleID]) index--;
      index = aggrLabels_[sampleID][index];
      cellsOccupied_[index] = -(cellsOccupied_[index] + 1);
      itmp = index;
      for (inputID = 0; inputID < nInputs_; inputID++)
      {
         jtmp = itmp % n1d_;
         itmp = itmp / n1d_;
         if (randFlag == 1)
              dtmp = ((double) (jtmp + 0.999*PSUADE_drand())) / (double) n1d_;
         else dtmp = ((double) (jtmp + 0.5)) / (double) n1d_;
         sampleMatrix_[sampleID][inputID] = dtmp * ranges[inputID] +
                                            lbounds[inputID];
         if (sampleMatrix_[sampleID][inputID] < lowerBounds_[inputID])
            sampleMatrix_[sampleID][inputID] = lowerBounds_[inputID];
         if (sampleMatrix_[sampleID][inputID] > upperBounds_[inputID])
            sampleMatrix_[sampleID][inputID] = upperBounds_[inputID];
      }
   }

   if (printLevel_ > 4)
   {
      printf("UserMetisSampling::initialize: nSamples = %d\n", nSamples_);
      printf("UserMetisSampling::initialize: nInputs  = %d\n", nInputs_);
      printf("UserMetisSampling::initialize: nOutputs = %d\n", nOutputs_);
      if (randFlag != 0)
           printf("UserMetisSampling::initialize: randomize on\n");
      else printf("UserMetisSampling::initialize: randomize off\n");
      for (inputID = 0; inputID < nInputs_; inputID++)
         printf("    UserMetisSampling input %3d = [%e %e]\n", inputID+1,
                lowerBounds_[inputID], upperBounds_[inputID]);
   }

   delete [] graphI;
   delete [] graphJ;
   delete [] ranges;
   delete [] ubounds;
   delete [] lbounds;
   delete [] aggrMap;
   return 0;
}

// ************************************************************************
// refine the sample space
// ------------------------------------------------------------------------
int UserMetisSampling::refine(int,int,double, int, double *)
{
   printf("UserMetisSampling refine: not implemented yet.\n");
   return -1;
} 

// ************************************************************************
// set internal scheme
// ------------------------------------------------------------------------
int UserMetisSampling::setParam(char *sparam)
{
   char winput[1001];
   FILE *fp;

   sscanf(sparam, "%s", winput);
   if (!strcmp(winput, "reset"))
   {
      fp = fopen("psuadeMetisInfo", "r");
      if (fp != NULL)
      {
         fclose(fp);
         unlink("psuadeMetisInfo");
      }
   }
   return 0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
UserMetisSampling& UserMetisSampling::operator=(const UserMetisSampling & ums)
{
   if(this == &ums) return *this;

   int graphN = 1;
   n1d_ = ums.n1d_;
   nAggrs_ = ums.nAggrs_;
   if(nAggrs_ > 0)
   {
      aggrCnts_ = new int[nAggrs_];
      for( int i = 0; i < nAggrs_; i++)
         aggrCnts_[i] = ums.aggrCnts_[i];
  
      aggrLabels_ = new int*[nAggrs_];
      for(int i= 0; i < nAggrs_; i++)
      {
         aggrLabels_[i] = new int[aggrCnts_[i]];
         for(int j = 0; j < aggrCnts_[i]; j++)
            aggrLabels_[i][j] = ums.aggrLabels_[i][j];
      }
   }
   nInputs_ = ums.nInputs_;
   for(int i = 0; i < nInputs_; i++) graphN *= n1d_;

   cellsOccupied_ = new int[graphN];
   for(int i = 0; i < graphN; i++)
      cellsOccupied_[i] = ums.cellsOccupied_[i];
  
   // MetisSampling inherits from Sampling so include the parent '
   // class data members
   printLevel_ = ums.printLevel_;
   samplingID_ = ums.samplingID_;
   nSamples_ = ums.nSamples_;
   nOutputs_ = ums.nOutputs_;
   randomize_ = ums.randomize_;
   nReplications_ = ums.nReplications_;
   lowerBounds_ = new double[nInputs_];
   upperBounds_ = new double[nInputs_];
   for (int i = 0; i < nInputs_; i++)
   {
      lowerBounds_[i] = ums.lowerBounds_[i];
      upperBounds_[i] = ums.upperBounds_[i];
   }
   sampleMatrix_ = new double*[nSamples_];
   for (int i = 0; i < nSamples_; i++)
   {
      sampleMatrix_[i] = new double[nInputs_];
      for(int j = 0; j < nInputs_; j++)
         sampleMatrix_[i][j] = ums.sampleMatrix_[i][j];
   }
   sampleOutput_ = new double[nSamples_*nOutputs_];
   for (int i = 0; i < nSamples_*nOutputs_; i++)
      sampleOutput_[i] = ums.sampleOutput_[i];
   sampleStates_ = new int[nSamples_];
   for (int i = 0; i < nSamples_; i++)
      sampleStates_[i] = ums.sampleStates_[i];

   return (*this);
}

