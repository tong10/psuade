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
// Functions for the Metis sampling class 
// AUTHOR : CHARLES TONG
// DATE   : 2004
// ************************************************************************
#include <stdio.h>
#include <assert.h>
#include <sstream>
#include <unistd.h>

#include "Psuade.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "FuncApprox.h"
#include "Mars.h"
#include "MetisSampling.h"

#define PABS(x) ((x) > 0 ? x : -(x))

#ifdef HAVE_METIS
extern "C" 
{
void METIS_PartGraphRecursive(int *, int *, int *, int *, int *,
                              int *, int *, int *, int *, int *, int *);
}
#endif

int modeSet_ = 0;
int pruneSet_ = 0;

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
MetisSampling::MetisSampling() : Sampling()
{
   samplingID_ = PSUADE_SAMP_METIS;
   refineType_ = 0;
   refineSize_ = 1000000;
   n1d_        = -1;
   nAggrs_     = 0;
   aggrCnts_   = NULL;
   aggrLabels_ = NULL;
   graphN_     = 0;
   graphI_     = NULL;
   graphJ_     = NULL;
   cellsOccupied_ = NULL;
   changeInfoName_ = 0;
}

// ************************************************************************
// copy constructor added by Bill Oliver
// ------------------------------------------------------------------------
MetisSampling::MetisSampling(const MetisSampling & ms) : Sampling()
{
   refineType_ = ms.refineType_;
   refineSize_ = ms.refineSize_;
   n1d_ = ms.n1d_;
   nAggrs_ = ms.nAggrs_;
   graphN_ = ms.graphN_;
 
   aggrLabels_ = new int*[nAggrs_]; 
   aggrCnts_ = new int[nAggrs_];
   for(int i = 0; i < nAggrs_; i++)
   {
      aggrCnts_[i] = ms.aggrCnts_[i];
      for(int j = 0; j < aggrCnts_[j]; j++)
      {
         aggrLabels_[j] = new int[aggrCnts_[i]];
         aggrLabels_[i][j] = ms.aggrLabels_[i][j];
      }
   }
   
   // MetisSampling inherits from Sampling so include the parent 
   // class data members
   printLevel_ = ms.printLevel_;
   samplingID_ = ms.samplingID_;
   nSamples_ = ms.nSamples_;
   nInputs_ = ms.nInputs_;
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

   graphI_ = new int[graphN_+1];
   for(int i = 0; i <= graphN_; i++)
      graphI_[i] = ms.graphI_[i];

   graphJ_ = new int[graphN_*nInputs_*2 + 1];
   for(int i = 0; i <= graphN_*nInputs_*2; i++)
      graphJ_[i] = ms.graphJ_[i];

   cellsOccupied_ = new int[graphN_];
   for(int i = 0; i < graphN_; i++)
      cellsOccupied_[i] = ms.cellsOccupied_[i];
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
MetisSampling::~MetisSampling()
{
   if (aggrCnts_ != NULL && nAggrs_ > 0) delete [] aggrCnts_;
   if (aggrLabels_ != NULL && nAggrs_ > 0)
   {
      for (int ii = 0; ii < nAggrs_; ii++)
         if (aggrLabels_[ii] != NULL) delete [] aggrLabels_[ii];
      delete [] aggrLabels_;
   }
   if (graphI_ != NULL) delete [] graphI_;
   if (graphJ_ != NULL) delete [] graphJ_;
   if (cellsOccupied_ != NULL) delete [] cellsOccupied_;
}

// ************************************************************************
// initialize the sampler data
// ------------------------------------------------------------------------
int MetisSampling::initialize(int initLevel)
{
   int    *incrs, inputID, ii, jj, itmp, jtmp, nnz, sampleID, randFlag;
   int    options[10], index, count, kk, saveFlag=1;
#ifdef HAVE_METIS
   int    wgtflag=0, numflag=0, edgeCut=0;
#endif
   double *ranges=NULL, dtmp, *lbounds, *ubounds, expand;
   char   response[1001], pString[1001], filename[1001];
   FILE   *fp;
 
   if( nSamples_ <= 0)
   {
      printf("nSamples_ in file %s line %d is <= 0\n", __FILE__, __LINE__);
      exit(1);
   }
   if (nInputs_ == 0 || lowerBounds_ == NULL || upperBounds_ == NULL)
   {
      printf("MetisSampling::initialize ERROR - input not set up.\n");
      exit(1);
   }

   deleteSampleData();
   if (aggrCnts_ != NULL && nAggrs_ > 0) delete [] aggrCnts_;
   if (aggrLabels_ != NULL && nAggrs_ > 0)
   {
      for (ii = 0; ii < nAggrs_; ii++)
         if (aggrLabels_[ii] != NULL) delete [] aggrLabels_[ii];
      delete [] aggrLabels_;
   }
   if (graphI_ != NULL) delete [] graphI_;
   if (graphJ_ != NULL) delete [] graphJ_;
   if (cellsOccupied_ != NULL) delete [] cellsOccupied_;
   aggrCnts_      = NULL;
   aggrLabels_    = NULL;
   cellsOccupied_ = NULL;
   randFlag  = randomize_;

   if (nInputs_ > 21)
   {
      printf("MetisSampling ERROR : nInputs > 21 currently not supported.\n");
      exit(1);
   }
   if (nInputs_ == 1 ) n1d_ = nSamples_*10;
   if (nInputs_ == 2 ) n1d_ = 1024;
   if (nInputs_ == 3 ) n1d_ = 100;
   if (nInputs_ == 4 ) n1d_ = 36;
   if (nInputs_ == 5 ) n1d_ = 16;
   if (nInputs_ == 6 ) n1d_ = 11;
   if (nInputs_ == 7 ) n1d_ = 8;
   if (nInputs_ == 8 ) n1d_ = 6;
   if (nInputs_ == 9 ) n1d_ = 5;
   if (nInputs_ == 10) n1d_ = 4;
   if (nInputs_ == 11) n1d_ = 3;
   if (nInputs_ == 12) n1d_ = 3;
   if (nInputs_ == 13) n1d_ = 3;
   if (nInputs_ >= 14) n1d_ = 2;

   incrs   = new int[nInputs_+1];
   graphN_ = 1;
   incrs[0] = graphN_;
   for (inputID = 1; inputID <= nInputs_; inputID++)
   {
      graphN_ *= n1d_;
      incrs[inputID] = graphN_;
   }
   if (nSamples_ > 2*graphN_)
   {
      printf("MetisSampling ERROR : nSamples %d too large.\n",nSamples_);
      exit(1);
   }

   graphI_ = new int[graphN_+1];
   graphJ_ = new int[graphN_*nInputs_*2+1];
   nnz = 0;
   graphI_[0] = nnz;
   for (ii = 0; ii < graphN_; ii++)
   {
      itmp = ii;
      for (inputID = 0; inputID < nInputs_; inputID++)
      {
         jtmp = itmp % n1d_;
         itmp = itmp / n1d_;
         if (jtmp > 0     ) graphJ_[nnz++] = ii - incrs[inputID];
         if (jtmp < n1d_-1) graphJ_[nnz++] = ii + incrs[inputID];
      }
      graphI_[ii+1] = nnz;
   }
   delete [] incrs;
   cellsOccupied_ = new int[graphN_];

   if (changeInfoName_ == 0) fp = fopen("psuadeMetisInfo", "r");
   else                      fp = fopen("psuadeMetisInfo.tmp", "r");
   if (fp != NULL)
   {
      printf("INFO: psuadeMetisInfo file found. Reading it in ...\n");
      fscanf(fp, "%d %d %d", &jj, &itmp, &jtmp);
      if (itmp != nSamples_ || jtmp != nInputs_ || jj != nSamples_)
      {
         fclose(fp);
         printf("MetisSampling INFO: a partition file is found but\n");
         printf("      the data is not consistent with this setup\n");
         printf("      (The file name is psuadeMetisInfo).\n");
         if (itmp != nSamples_ || jj != nSamples_)
            printf("      nSamples : %d != %d.\n", nSamples_, itmp);
         if (jtmp != nInputs_)
            printf("      nInputs  : %d != %d.\n", nInputs_, jtmp);
         sprintf(pString,"Would you like to provide another file? (y or n) ");
         getString(pString, response);
         if (response[0] == 'y')
         {
            sprintf(pString,"Enter partition file name : ");
            getString(pString, filename);
            fp = fopen(filename, "r");
            if (fp == NULL)
            {
               printf("MetisSampling ERROR: partition file not found.\n");
               exit(1);
            }
         }
         else
         {
            printf("MetisSampling INFO: delete psuadeMetisInfo file and.\n");
            printf("                    re-launch.\n");
            exit(1);
         }
         fscanf(fp, "%d %d %d", &jj, &itmp, &jtmp);
         if (itmp != nSamples_ || jtmp != nInputs_)
         {
            printf("MetisSampling ERROR: partition file not valid.\n");
            exit(1);
         }
         printf("Metis INFO: partition file has %d subdomains\n",jj);
         printf("            and %d sample points.\n",itmp);
         printf("Metis INFO: reconstructing the partitioning.\n");
      }
      printf("      Incoming nSamples : %d.\n", itmp);
      printf("      Incoming nInputs  : %d.\n", jtmp);
      nAggrs_ = jj;
      aggrCnts_ = new int[nAggrs_];
      aggrLabels_ = new int*[nAggrs_];
      for (ii = 0; ii < nAggrs_; ii++)
      {
         fscanf(fp, "%d", &count);
         if (printLevel_ > 4) 
            printf("Metis read: aggr %8d, size = %d\n", ii+1, count);
         if (count > 0)
         {
            aggrCnts_[ii] = count;
            aggrLabels_[ii] = new int[count];
         }
         else 
         {
            aggrCnts_[ii] = count = 0;
            aggrLabels_[ii] = NULL;
         }
         for (jj = 0; jj < count; jj++)
         {
            fscanf(fp, "%d", &(aggrLabels_[ii][jj]));
            if (aggrLabels_[ii][jj] < 0 || aggrLabels_[ii][jj] >= graphN_)
            {
               printf("Metis ERROR: psuadeMetisInfo file has invalid info.\n");
               exit(1);
            }
            cellsOccupied_[aggrLabels_[ii][jj]] = ii;
         }
      }
      fscanf(fp, "%d", &count);
      for (ii = 0; ii < count; ii++)
      {
         fscanf(fp, "%d %d", &jj, &kk);
         if (jj < 0 || jj >= graphN_)
         {
            printf("Metis ERROR: psuadeMetisInfo file has invalid info (2).\n");
            exit(1);
         }
         cellsOccupied_[jj] = - kk - 1;
      }
      fclose(fp);
      saveFlag = 0;
   }

   if (aggrCnts_ == NULL)
   {
      options[0] = 0;
      if (printLevel_ > 1)
         printf("MetisSampling:: calling domain partitioner.\n");
#ifdef HAVE_METIS
      METIS_PartGraphRecursive(&graphN_, graphI_, graphJ_, NULL, NULL,
                &wgtflag,&numflag,&nSamples_,options,&edgeCut,cellsOccupied_);
#else
      printf("MetisSampling ERROR : METIS not installed.\n");
      exit(1);
#endif
      if (printLevel_ > 1)
         printf("MetisSampling:: subdomains created.\n");

      nAggrs_   = nSamples_;
      aggrCnts_ = new int[nAggrs_];
      for (ii = 0; ii < nAggrs_; ii++) aggrCnts_[ii] = 0;
      for (ii = 0; ii < graphN_; ii++)
      {
         if (cellsOccupied_[ii] < 0 || cellsOccupied_[ii] >= nAggrs_)
         {
            printf("MetisSampling INTERNAL ERROR.\n");
            exit(1);
         }
	 aggrCnts_[cellsOccupied_[ii]]++;
      }  
      aggrLabels_ = new int*[nAggrs_];
      for (ii = 0; ii < nAggrs_; ii++)
      {
         if (aggrCnts_[ii] <= 0)
         {
            printf("MetisSampling INTERNAL ERROR (2).\n");
            exit(1);
         }
         aggrLabels_[ii] = new int[aggrCnts_[ii]];
         aggrCnts_[ii] = 0;
      }
      for (ii = 0; ii < graphN_; ii++)
      {
         index = cellsOccupied_[ii];
         if (index < 0 || index >= nAggrs_)
         {
            printf("MetisSampling INTERNAL ERROR (3).\n");
            exit(1);
         }
         aggrLabels_[index][aggrCnts_[index]++] = ii;  
      }
   }

   if (saveFlag == 1)
   {
      if (changeInfoName_ == 0) fp = fopen("psuadeMetisInfo", "w");
      else                      fp = fopen("psuadeMetisInfo.tmp", "w");
      if (fp != NULL)
      {
         fprintf(fp, "%d %d %d\n", nAggrs_, nSamples_, nInputs_);
         for (ii = 0; ii < nAggrs_; ii++)
         {
            fprintf(fp, "%d\n", aggrCnts_[ii]);
            for (jj = 0; jj < aggrCnts_[ii]; jj++)
            {
               fprintf(fp, "%d ", aggrLabels_[ii][jj]);
               if (jj != 0 && jj % 10 == 0) fprintf(fp, "\n");
            }
            fprintf(fp, "\n");
         }
         jj = 0;
         for (ii = 0; ii < graphN_; ii++) if (cellsOccupied_[ii] < 0) jj++;
         fprintf(fp, "%d\n", jj);
         for (ii = 0; ii < graphN_; ii++)
            if (cellsOccupied_[ii] < 0) 
               fprintf(fp, "%d %d\n",ii,-(cellsOccupied_[ii]+1));
         fclose(fp);
      }
   }
   if (initLevel != 0) return 0;

   char *inStr, winput1[500], winput2[500];

   expand = 0.0;
   if (psSamExpertMode_ == 1)
   {
      printf("For the Metis sampling, you can divert more data points\n");
      printf("to the surface of the parameter space by sampling from\n");
      printf("the expanded parameter space folowed by projecting the\n");
      printf("sample points back to the original parameter space.\n");
      printf("In order to do so, an expansion ratio needs to be given.\n");
      printf("An expansion ratio of 0.0 means no expansion.\n");
      printf("Usually the expansion ratio should be no more than 0.1-0.2.\n");
      sprintf(pString, "Enter an expansion ratio (default=0): ");
      expand = -0.1;
      while (expand < 0.0) expand = getDouble(pString);
   }
   else if (psConfig_ != NULL)
   {
      inStr = psConfig_->getParameter("METIS_expand_ratio");
      if (inStr != NULL)
      {
         sscanf(inStr, "%s %s %lg\n", winput1, winput2, &expand);
         if (winput2[0] != '=')
         {
            printf("METIS read config file syntax error : %s\n", inStr);
            expand = 0.0;
         }
         else if (expand < 0.0 || expand >= 1.0)
         {
            printf("METIS read config file read expand error : %e\n", expand);
            expand = 0.0;
         }
         else printf("METIS : expansion ratio set to %e\n", expand);
      }
   }

   allocSampleData();
   ranges  = new double[nInputs_];
   lbounds = new double[nInputs_];
   ubounds = new double[nInputs_];
   for (inputID = 0; inputID < nInputs_; inputID++)
   {
      ranges[inputID]  = (upperBounds_[inputID] - lowerBounds_[inputID]) *
                         0.5 * expand;
      lbounds[inputID] = lowerBounds_[inputID] - ranges[inputID];
      ubounds[inputID] = upperBounds_[inputID] + ranges[inputID];
   }
   for (inputID = 0; inputID < nInputs_; inputID++)
      ranges[inputID]  = ubounds[inputID] - lbounds[inputID];

   for (sampleID = 0; sampleID < nSamples_; sampleID++)
   {
      if ((randFlag & 1) == 1)
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
         if ((randFlag & 1) == 1)
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
      printf("MetisSampling::initialize: nSamples = %d\n", nSamples_);
      printf("MetisSampling::initialize: nInputs  = %d\n", nInputs_);
      printf("MetisSampling::initialize: nOutputs = %d\n", nOutputs_);
      if (randomize_ != 0)
           printf("MetisSampling::initialize: randomize on\n");
      else printf("MetisSampling::initialize: randomize off\n");
      for (inputID = 0; inputID < nInputs_; inputID++)
         printf("    MetisSampling input %3d = [%e %e]\n", inputID+1,
                lowerBounds_[inputID], upperBounds_[inputID]);
   }

   delete [] ranges;
   delete [] ubounds;
   delete [] lbounds;
   return 0;
}

// ************************************************************************
// refine the sample space
// ------------------------------------------------------------------------
int MetisSampling::refine(int nLevels, int randFlag, double threshold,
                          int nSamples, double *sampleErrors)
{
   int    inputID, ii, jj, count, itmp, jtmp, localN, *localIA, *localJA;
   int    maxN, rowInd, colInd, localNNZ, *subLabels, index, count0, kk;
   int    options[10], count1, iOne=1;
#ifdef HAVE_METIS
   int    wgtflag=0, numflag=0, edgeCut=0, itwo=2;
#endif
   int    *tmpCnts, **tmpLabels, currNAggr, newCell, oldNumSamples, mode;
   int    *oldSampleStates, *labels, status, outputID, ss, maxNNZ, marsCnt;
   int    *refineArray=NULL, *node2Aggr=NULL, aggr1, splitCount, *tags;
   int    splitSuccess, maxCellSize, minCellSize, avgCellSize, cellCnt=0;
   int    *tagArray, aggr2, rowInd2, colInd2, localN2, ii2, jj2;
   int    *subLabels2, limit;
   double *ranges, dtmp, **oldSampleMatrix, *oldSampleOutput;
   double *lbounds, *ubounds, *diffArray=NULL, ddmax, *sortList2, diffCnt;
   double *sortList1, ddmin, ddmean, thresh, *marsIn, *marsOut, ddsd;
   char   pString[500];
   FILE   *fp;
   FuncApprox *faPtr;

   if (nSamples != 0 && nSamples != nAggrs_)
   {
      printf("MetisSampling::refine ERROR - sampleErrors length mismatch.\n");
      printf("               This may be due to the use of wrong\n");
      printf("               psuadeMetisInfo file.\n");
      exit(1);
   }
   if (cellsOccupied_ == NULL || aggrLabels_ == NULL)
   {
      printf("MetisSampling::refine ERROR - need to call initialize first.\n");
      exit(1);
   }
   if (printLevel_ > 0)
   {
      printf("MetisSampling::refine(1): nSamples = %d\n", nSamples_);
      printf("MetisSampling::refine(1): nInputs  = %d\n", nInputs_);
      printf("MetisSampling::refine(1): nOutputs = %d\n", nOutputs_);
   }

   ranges  = new double[nInputs_];
   lbounds = new double[nInputs_];
   ubounds = new double[nInputs_];
   for (inputID = 0; inputID < nInputs_; inputID++)
   {
      ranges[inputID]  = (upperBounds_[inputID] - lowerBounds_[inputID])*0.0;
      lbounds[inputID] = lowerBounds_[inputID] - ranges[inputID];
      ubounds[inputID] = upperBounds_[inputID] + ranges[inputID];
   }
   for (inputID = 0; inputID < nInputs_; inputID++)
   {
      ranges[inputID] = ubounds[inputID] - lbounds[inputID];
      if (ranges[inputID] <= 0.0)
      {
         printf("MetisSampling::refine ERROR - lbound/ubound mismatch.\n");
         exit(1);
      }
   }

   count = 0;
   for (ii = 0; ii < graphN_; ii++) if (cellsOccupied_[ii] < 0) count++;

   if (count == 0)
   {
      for (ss = 0; ss < nSamples_; ss++)
      {
         itmp = 0;
         for (inputID = nInputs_-1; inputID >= 0; inputID--)
         {
            itmp = itmp * n1d_;
            dtmp = sampleMatrix_[ss][inputID];
            dtmp = (dtmp - lowerBounds_[inputID]) / ranges[inputID];
            if (dtmp == 1.0) jtmp = n1d_ - 1;
            else             jtmp = (int) (dtmp * n1d_);
            itmp += jtmp;
         }
         if (itmp < 0 || itmp >= graphN_)
         {
            printf("MetisSampling::refine INTERNAL ERROR.\n");
            printf("               Consult PSUADE developer.\n");
         }
         cellsOccupied_[itmp] = -(cellsOccupied_[itmp] + 1); 
      }
   }
   else if (count != nSamples_)
   {
      printf("MetisSampling::refine ERROR - in re-initialize.\n");
      printf("               Consult PSUADE developer.\n");
      printf("               Sample size mismatch: %d vs %d\n",
             count, nSamples_);
      printf("                between psuadeMetisInfo and current data set.\n");
      exit(1);
   }
   
   count = 0;
   for (ii = 0; ii < graphN_; ii++) 
      if (cellsOccupied_[ii] < 0) count -= (cellsOccupied_[ii] + 1);
   itmp = (nSamples_ - 1) * nSamples_ / 2;
   if (nSamples_ < 40000 && count != itmp)
   {
      printf("MetisSampling::refine ERROR - METIS not used in initialize.\n");
      printf("               so Metis cannot be used in refine (%d,%d).\n",
             count, itmp);
      printf("Note: It can be due to psuadeMetisInfo file being modified.\n");
      exit(1);
   }
   else if (nSamples_ >= 40000)
   {
      printf("MetisSampling::refine INFO - less error checking (N>40000).\n");
   }

   tmpCnts = aggrCnts_;
   tmpLabels = aggrLabels_;
   aggrCnts_ = new int[2*nAggrs_];
   aggrLabels_ = new int*[2*nAggrs_];
   for (ii = 0; ii < nAggrs_; ii++)
   {
      aggrCnts_[ii] = tmpCnts[ii];
      aggrLabels_[ii] = new int[aggrCnts_[ii]];
      for (jj = 0; jj < aggrCnts_[ii]; jj++) 
         aggrLabels_[ii][jj] = tmpLabels[ii][jj];
   }
   delete [] tmpCnts;
   for (ii = 0; ii < nAggrs_; ii++)
      if (tmpLabels[ii] != NULL) delete [] tmpLabels[ii];
   delete [] tmpLabels;
   for (ii = nAggrs_; ii < 2*nAggrs_; ii++) 
   {
      aggrCnts_[ii] = 0;
      aggrLabels_[ii] = NULL;
   }

   maxN = 0;
   for (ss = 0; ss < nAggrs_; ss++)
      if (aggrCnts_[ss] > maxN) maxN = aggrCnts_[ss];
   labels  = new int[maxN];
   assert(labels != NULL);
   localIA = new int[2*maxN+1];
   assert(localIA != NULL);
   maxNNZ = maxN * (2 * nInputs_ + 1);
   localJA = new int[maxNNZ];
   assert(localJA != NULL);

   refineArray = NULL;

   if (refineType_ == 1 && sampleErrors == NULL)
   {
      printf("MetisSampling::refine ERROR - error based but no error given.\n");
      exit(1);
   }
   if (refineType_ == 1 && sampleErrors != NULL)
   {
      printf("MetisSampling::refine - maximum number of new points = %d.\n",
             refineSize_);
      printf("Note: error scaled by the cell size.\n");
      node2Aggr = new int[graphN_];
      for (ii = 0; ii < graphN_; ii++) node2Aggr[ii] = -1;
      for (ii = 0; ii < nAggrs_; ii++)
      {
         localN = aggrCnts_[ii];
         subLabels = aggrLabels_[ii];
         for (jj = 0; jj < localN; jj++)
         {
            rowInd = subLabels[jj];
            if (rowInd < 0 || rowInd >= graphN_)
            {
               printf("MetisSampling::refine ERROR (index out of bound.)\n");
               printf("               rowInd = %d (should be in [0,%d])\n",
                      rowInd, graphN_-1);
               printf("               Aggregate = %d (loc=%d)\n", ii, jj);
               printf("               Please consult PSUADE developers.\n");
               exit(1);
            }
            node2Aggr[rowInd] = ii;
         }
      }
      for (ii = 0; ii < graphN_; ii++)
      {
         if (node2Aggr[ii] == -1)
         {
            printf("MetisSampling::refine ERROR - node2Aggr not correct.\n");
            exit(1);
         }
      }

      sortList1 = new double[nAggrs_];
      sortList2 = new double[nAggrs_];
      for (ii = 0; ii < nAggrs_; ii++)
      {
         localN = aggrCnts_[ii];
         sortList1[ii] = PABS(sampleErrors[ii]);
         sortList1[ii] *= pow(1.0*localN,1.0/nInputs_);
      }
      for (ii = 0; ii < nAggrs_; ii++) sortList2[ii] = (double) ii;
      sortDbleList2(nAggrs_, sortList1, sortList2);
      refineArray = new int[nAggrs_];
      for (ii = nAggrs_-1; ii >= 0; ii--) refineArray[ii] = 0;
      cellCnt = 0;

      for (ii = nAggrs_-1; ii >= 0; ii--)
      {
         index = (int) sortList2[ii];
         localN = aggrCnts_[index];
         if (sampleErrors[index] != 0.0) 
         {
            //printf("MetisSampling::refine - chosen aggregate, error = %e\n",
            //        sampleErrors[index]);
            printf("MetisSampling: selected for refinement: error = %13.5e, size = %d\n",
                   sampleErrors[index], localN);
            if (refineArray[index] == 0)
            {
               cellCnt++;
               refineArray[index] = 1;
               if (cellCnt >= refineSize_) break;
            }
#if 0
            subLabels = aggrLabels_[index];
            for (jj = 0; jj < localN; jj++)
            {
               rowInd = subLabels[jj];
               for (kk = graphI_[rowInd]; kk < graphI_[rowInd+1]; kk++)
               {
                  colInd = graphJ_[kk];
                  aggr1 = node2Aggr[colInd];
                  if (refineArray[aggr1] == 0)
                  {
                     cellCnt++;
                     refineArray[aggr1] = 1;
                     if (cellCnt >= refineSize_) break;
                  }
               }
               if (cellCnt >= refineSize_) break;
            }
#endif
            if (cellCnt >= refineSize_) break;
         }
         if (cellCnt >= refineSize_) break;
      }
      delete [] sortList1;
      delete [] sortList2;
      delete [] node2Aggr;
   }

   if (refineType_ == 2)
   {
      printf("MetisSampling::refine - maximum number of new points = %d.\n",
             refineSize_);
      tags = new int[nAggrs_];
      for (ii = 0; ii < nAggrs_; ii++) tags[ii] = 0;
      node2Aggr = new int[graphN_];
      for (ii = 0; ii < nAggrs_; ii++)
      {
         localN = aggrCnts_[ii];
         subLabels = aggrLabels_[ii];
         for (jj = 0; jj < localN; jj++)
         {
            rowInd = subLabels[jj];
            if (rowInd < 0 || rowInd >= graphN_)
            {
               printf("MetisSampling::refine ERROR (index out of bound.)\n");
               printf("               rowInd = %d (should be in [0,%d])\n",
                      rowInd, graphN_-1);
               printf("               Aggregate = %d (loc=%d)\n", ii, jj);
               printf("               Please consult PSUADE developers.\n");
               printf("               next loc index = %d\n", subLabels[jj+1]);
               exit(1);
            }
            node2Aggr[rowInd] = ii;
         }
      }

      diffArray = new double[nAggrs_];
      for (ss = 0; ss < nAggrs_; ss++)
      {
         diffArray[ss] = 0;
         localN = aggrCnts_[ss];
         subLabels = aggrLabels_[ss];

         ddmean = sampleOutput_[ss];
         count0 = 1;
         for (ii = 0; ii < localN; ii++)
         {
            rowInd = subLabels[ii];
            for (jj = graphI_[rowInd]; jj < graphI_[rowInd+1]; jj++)
            {
               colInd = graphJ_[jj];
               aggr1 = node2Aggr[colInd];
               if (aggr1 != ss && tags[aggr1] == 0)
               {
                  tags[aggr1] = 1;
                  ddmean += sampleOutput_[aggr1];
                  count0++;
                  localN2 = aggrCnts_[aggr1];
                  subLabels2 = aggrLabels_[aggr1];
                  for (ii2 = 0; ii2 < localN2; ii2++)
                  {
                     rowInd2 = subLabels2[ii2];
                     for (jj2 = graphI_[rowInd2];jj2 < graphI_[rowInd2+1];jj2++)
                     {
                        colInd2 = graphJ_[jj2];
                        aggr2 = node2Aggr[colInd2];
                        if (aggr2 != ss && tags[aggr2] == 0)
                        {
                           tags[aggr2] = 1;
                           ddmean += sampleOutput_[aggr2];
                           count0++;
                        }
                     }
                  }
               }
            }
         }
         ddmean /= count0;
         ddsd = 0.0;
         if (count0 > 1)
         {
            ddsd += (sampleOutput_[ss] - ddmean) * (sampleOutput_[ss] - ddmean);
            for (ii = 0; ii < localN; ii++)
            {
               rowInd = subLabels[ii];
               for (jj = graphI_[rowInd]; jj < graphI_[rowInd+1]; jj++)
               {
                  colInd = graphJ_[jj];
                  aggr1 = node2Aggr[colInd];
                  if (aggr1 != ss && tags[aggr1] == 1)
                  {
                     dtmp = sampleOutput_[aggr1];
                     ddsd += (dtmp - ddmean) * (dtmp - ddmean);
                     tags[aggr1] = 0;
                     localN2 = aggrCnts_[aggr1];
                     subLabels2 = aggrLabels_[aggr1];
                     for (ii2 = 0; ii2 < localN2; ii2++)
                     {
                        rowInd2 = subLabels2[ii2];
                        for (jj2=graphI_[rowInd2];jj2<graphI_[rowInd2+1];jj2++)
                        {
                           colInd2 = graphJ_[jj2];
                           aggr2 = node2Aggr[colInd2];
                           if (aggr2 != ss && tags[aggr2] == 1)
                           {
                              tags[aggr2] = 0;
                              dtmp = sampleOutput_[aggr2];
                              ddsd += (dtmp - ddmean) * (dtmp - ddmean);
                           }
                        }
                     }
                  }
               }
            }
         }
         if (count0 > 0) diffArray[ss] = sqrt(ddsd) / count0;

         if (printLevel_ > 4 && nInputs_ == 2)
         {
            for (ii = 0; ii < nAggrs_; ii++) tags[ii] = 0;
            printf(" ..... Sample %6d inputs : %e %e = %e\n", ss+1, 
                   sampleMatrix_[ss][0],sampleMatrix_[ss][1],
                   sampleOutput_[ss]);
            for (ii = 0; ii < localN; ii++)
            {
               rowInd = subLabels[ii];
               for (jj = graphI_[rowInd]; jj < graphI_[rowInd+1]; jj++)
               {
                  colInd = graphJ_[jj];
                  aggr1 = node2Aggr[colInd];
                  if (aggr1 != ss && tags[aggr1] == 0)
                  {
                     printf(" .....          neighbor : %e %e = %e\n", 
                         sampleMatrix_[aggr1][0],sampleMatrix_[aggr1][1],
                         sampleOutput_[aggr1]);
                     tags[aggr1] = 1;
                  }
               }
            }
         }
      }
      delete [] tags;

      diffCnt = 0.0;
      for (ii = 0; ii < nAggrs_; ii++) diffCnt += diffArray[ii];

      refineArray = new int[nAggrs_];
      for (ii = 0; ii < nAggrs_; ii++) refineArray[ii] = 0;
      maxCellSize = cellCnt = 0;
      avgCellSize = 0;
      minCellSize = 1000000000;
      if (diffCnt == 0 || diffCnt == nAggrs_)
      { 
         for (ii = 0; ii < nAggrs_; ii++) 
         {
            refineArray[ii] = 0;
            if (localN > 1)
            {
               refineArray[ii] = 1;
               localN = aggrCnts_[ii];
               maxCellSize = (localN > maxCellSize) ? localN : maxCellSize;
               minCellSize = (localN < maxCellSize) ? localN : minCellSize;
               avgCellSize += localN;
               cellCnt++;
            }
         }
      }
      else
      {
         mode = 2;
         if (psSamExpertMode_ == 1 && modeSet_ == 0)
         {
            modeSet_ = 1;
         }
         limit = 0;
         if (mode == 2)
         {
            for (ss = 0; ss < nAggrs_; ss++) diffArray[ss] *= aggrCnts_[ss];
         }
         else if (mode == 3 || mode == 4)
         {
            if (mode == 4)
            {
               limit = nInputs_ + 1 + nInputs_ * (nInputs_ + 1) / 2;
            }
            marsIn  = new double[nAggrs_*nInputs_];
            marsOut = new double[nAggrs_];
            tagArray = new int[nAggrs_];
            for (ss = 0; ss < nAggrs_; ss++)
            {
               for (ii = 0; ii < nAggrs_; ii++) tagArray[ii] = 0;
               tagArray[ss] = 1;
               marsCnt = 0;
               localN = aggrCnts_[ss];
               subLabels = aggrLabels_[ss];
               for (ii = 0; ii < localN; ii++)
               {
                  rowInd = subLabels[ii];
                  for (jj = graphI_[rowInd]; jj < graphI_[rowInd+1]; jj++)
                  {
                     colInd = graphJ_[jj];
                     aggr1 = node2Aggr[colInd];
                     if (tagArray[aggr1] == 0)
                     {
                        for (kk = 0; kk < nInputs_; kk++)
                           marsIn[marsCnt*nInputs_+kk]=sampleMatrix_[aggr1][kk];
                        marsOut[marsCnt++] = sampleOutput_[aggr1];
                        tagArray[aggr1] = 1;
                     }
                  }
               }
               if (marsCnt < limit)
               {
                  for (ii = 0; ii < localN; ii++)
                  {
                     rowInd = subLabels[ii];
                     for (jj = graphI_[rowInd]; jj < graphI_[rowInd+1]; jj++)
                     {
                        colInd = graphJ_[jj];
                        aggr1 = node2Aggr[colInd];
                        localN2 = aggrCnts_[aggr1];
                        subLabels2 = aggrLabels_[aggr1];
                        for (ii2 = 0; ii2 < localN2; ii2++)
                        {
                           rowInd2 = subLabels2[ii2];
                           for (jj2=graphI_[rowInd2];jj2<graphI_[rowInd2+1];jj2++)
                           {
                              colInd2 = graphJ_[jj2];
                              aggr2 = node2Aggr[colInd2];
                           }
                           if (tagArray[aggr2] == 0)
                           {
                              for (kk = 0; kk < nInputs_; kk++)
                                 marsIn[marsCnt*nInputs_+kk] = 
                                                   sampleMatrix_[aggr2][kk];
                               marsOut[marsCnt++] = sampleOutput_[aggr2];
                               tagArray[aggr2] = 1;
                           }
                        }
                     }
                  }
               }
               if (printLevel_ > 2)
                  printf("Metis refine aggr %d (%d): sample size = %d\n",
                         ss+1,nAggrs_,marsCnt);
               if (mode == 3)
                  faPtr = genFA(0, nInputs_, iOne, marsCnt);
               else
                  faPtr = genFA(2, nInputs_, iOne, marsCnt);
               faPtr->setBounds(lbounds, ubounds);
               faPtr->setOutputLevel(1);
               faPtr->initialize(marsIn, marsOut);
               diffArray[ss] = faPtr->evaluatePoint(sampleMatrix_[ss]);
               diffArray[ss] = PABS(diffArray[ss] - sampleOutput_[ss]);
               delete faPtr;
            }
            delete [] marsIn;
            delete [] marsOut;
            delete [] tagArray;
            dtmp = 0.0;
            for (ss = 0; ss < nAggrs_; ss++) dtmp += diffArray[ss];
            printf("Average interpolation error = %e\n", dtmp/nAggrs_);
         }

         thresh = 0.0;
         ddmax = -PSUADE_UNDEFINED;
         ddmin =  PSUADE_UNDEFINED;
         for (ii = 0; ii < nAggrs_; ii++)
         {
            if (diffArray[ii] > ddmax) ddmax = diffArray[ii];
            if (diffArray[ii] < ddmin) ddmin = diffArray[ii];
         }
         thresh = 0.0;
         if (psSamExpertMode_ == 1 && pruneSet_ == 0)
         {
            printf("Maximum discrepancy between aggregates = %e\n",ddmax);
            printf("Minimum discrepancy between aggregates = %e\n",ddmin);
            fp = fopen("arsmPDF.m","w");
            // add check for NULL by Bill Oliver
            if(fp != NULL)
            {
               fprintf(fp,"A = [\n");
               for (ii = 0; ii < nAggrs_; ii++)
                  fprintf(fp,"%e\n",diffArray[ii]);
               fprintf(fp,"];\n");
               fprintf(fp,"hist(A,10)\n");
               fclose(fp);
            }
            printf("A PDF for error has been given to you in arsmPDF.m.\n");
            printf("Default threshold for pruning = 0.\n");
            sprintf(pString,"Enter thresh to prune refinement candidates : ");
            thresh = getDouble(pString);
            if (thresh < 0.0)
            {
               printf("WARNING: threshold < 0 not allowed. Reset to 0.\n");
               thresh = 0.0;
            }
            pruneSet_ = 1;
         }
         for (ii = 0; ii < nAggrs_; ii++)
            if (diffArray[ii] < thresh) diffArray[ii] = 0.0;

         sortList2 = new double[nAggrs_];
         for (ii = 0; ii < nAggrs_; ii++) sortList2[ii] = (double) ii;
         sortDbleList2(nAggrs_, diffArray, sortList2);
         cellCnt = 0;
         ddmax = -PSUADE_UNDEFINED;
         for (ii = nAggrs_-1; ii >= 0; ii--) refineArray[ii] = 0;
         for (ii = nAggrs_-1; ii >= 0; ii--)
         {
            jj = (int) sortList2[ii];
            localN = aggrCnts_[jj];
            if (localN > 1 && diffArray[ii] > 0.0)
            {
               refineArray[jj] = 1;
               maxCellSize = (localN > maxCellSize) ? localN : maxCellSize;
               minCellSize = (localN < maxCellSize) ? localN : minCellSize;
               avgCellSize += localN;
               if (diffArray[ii] > ddmax) ddmax = diffArray[ii];
               cellCnt++;
            }
            if (cellCnt >= refineSize_) break;
         }
         //printf("MetisSampling:: maximum neighbor discrepancy = %e\n",ddmax);
         delete [] sortList2;
      }
      delete [] diffArray;
      delete [] node2Aggr;
      printf("MetisSampling::refine - number of new sample points = %d.\n",
             cellCnt);
   }

   options[0] = 0;
   currNAggr = nAggrs_;
   splitCount = splitSuccess = 0;
   for (ss = 0; ss < nAggrs_; ss++)
   {
      localN = aggrCnts_[ss];
      if (localN > maxN) 
      {
         printf("MetisSampling INTERNAL ERROR (1): consult PSUADE developers\n");
         exit(1);
      }
      if (refineArray == NULL || refineArray[ss] == 1) splitCount++;
      if (localN == 1 && (refineArray == NULL || refineArray[ss] == 1))
         printf("MetisSampling::refine INFO- cannot split cell (too small).\n");
      if (localN > 1 && (refineArray == NULL || refineArray[ss] == 1))
      {
         if (printLevel_ > 4)
            printf("MetisSampling::refine - split sample %d\n",ss+1);
         splitSuccess++;
         localNNZ = 0;
         subLabels = aggrLabels_[ss];
         localIA[0] = localNNZ;
         for (ii = 0; ii < localN; ii++)
         {
            rowInd = subLabels[ii];
            for (jj = graphI_[rowInd]; jj < graphI_[rowInd+1]; jj++)
            {
               colInd = graphJ_[jj];
               status = binarySearchInt(colInd,subLabels,localN);
               if (status >= 0) localJA[localNNZ++] = status;
            }
            if (ii >= maxN)
            {
               printf("MetisSampling INTERNAL ERROR (2): consult developers\n");
               exit(1);
            }
            localIA[ii+1] = localNNZ;
         }
         if (localNNZ > maxNNZ) 
         {
            printf("MetisSampling INTERNAL ERROR (3): consult developers\n");
            exit(1);
         }

#ifdef HAVE_METIS
         METIS_PartGraphRecursive(&localN, localIA, localJA, NULL, NULL,
                &wgtflag,&numflag,&itwo,options,&edgeCut,labels);
#else
         printf("MetisSampling ERROR : METIS not installed.\n");
         exit(1);
#endif

         count0 = 0;
         for (ii = 0; ii < localN; ii++) if (labels[ii] == 0) count0++;
         count1 = localN - count0;

         for (ii = 0; ii < localN; ii++) 
         {
            index = aggrLabels_[ss][ii];
            if (cellsOccupied_[index] < 0) break;
         }
         newCell = 0;
         if (labels[ii] == 0) newCell = 1;

         if (newCell == 0)
         {
            aggrCnts_[currNAggr] = count0;
            aggrLabels_[currNAggr] = new int[count0]; 
            count = 0;
            for (ii = 0; ii < localN; ii++) 
            {
               if (labels[ii] == 0) 
               {
                  aggrLabels_[currNAggr][count++] = aggrLabels_[ss][ii];
                  cellsOccupied_[aggrLabels_[ss][ii]] = currNAggr;
               }
            }
            if (count != count0)
            {
               printf("MetisSampling INTERNAL ERROR (4): consult developers\n");
               exit(1);
            }
            aggrCnts_[ss] = count1; 
            count = 0;
            for (ii = 0; ii < localN; ii++) 
               if (labels[ii] == 1) 
                  aggrLabels_[ss][count++] = aggrLabels_[ss][ii];
            if (count != count1)
            {
               printf("MetisSampling INTERNAL ERROR (5): consult developers\n");
               exit(1);
            }
         }
         else
         {
            aggrCnts_[currNAggr] = count1;
            aggrLabels_[currNAggr] = new int[count1]; 
            count = 0;
            for (ii = 0; ii < localN; ii++) 
            {
               if (labels[ii] == 1) 
               {
                  aggrLabels_[currNAggr][count++] = aggrLabels_[ss][ii];
                  cellsOccupied_[aggrLabels_[ss][ii]] = currNAggr;
               }
            }
            if (count != count1)
            {
               printf("MetisSampling INTERNAL ERROR (6): consult developers\n");
               exit(1);
            }
            aggrCnts_[ss] = count0; 
            count = 0;
            for (ii = 0; ii < localN; ii++) 
               if (labels[ii] == 0) 
                  aggrLabels_[ss][count++] = aggrLabels_[ss][ii];
            if (count != count0)
            {
               printf("MetisSampling INTERNAL ERROR (7): consult developers\n");
               exit(1);
            }
         }
         currNAggr++;
      }
   }
   if (printLevel_ > 4 && refineType_ == 2 && maxCellSize != 0)
      printf("MetisSampling:: maximum resolution = %d \n", maxCellSize);
   if (printLevel_ > 4 && refineType_ == 2 && minCellSize != 1000000000)
      printf("MetisSampling:: minimum resolution = %d \n", minCellSize);
   if (printLevel_ > 4 && refineType_ == 2 && cellCnt != 0)
      printf("MetisSampling:: average resolution = %e (%d)\n", 
             1.0*avgCellSize/cellCnt, cellCnt);
           
   if (printLevel_ > 4 && splitSuccess != splitCount)
      printf("MetisSampling:: number of successful splits = %d (out of %d)\n",
             splitSuccess, splitCount);

   oldNumSamples   = nSamples_;
   oldSampleMatrix = sampleMatrix_;
   oldSampleOutput = sampleOutput_;
   oldSampleStates = sampleStates_;
   nSamples_ = currNAggr;
   allocSampleData();
   for (ss = 0; ss < oldNumSamples; ss++)
   {
      for (inputID = 0; inputID < nInputs_; inputID++)
         sampleMatrix_[ss][inputID] = oldSampleMatrix[ss][inputID];
      for (outputID = 0; outputID < nOutputs_; outputID++)
         sampleOutput_[ss*nOutputs_+outputID] =
                    oldSampleOutput[ss*nOutputs_+outputID];
      sampleStates_[ss] = oldSampleStates[ss];
   }
   for (ss = 0; ss < oldNumSamples; ss++) delete [] oldSampleMatrix[ss];
   delete [] oldSampleMatrix;
   delete [] oldSampleOutput;
   delete [] oldSampleStates;

   for (ss = oldNumSamples; ss < currNAggr; ss++)
   {
      if ((randFlag & 1) == 1)
           index = (int) (PSUADE_drand() * aggrCnts_[ss]);
      else index = aggrCnts_[ss] / 2;
      if (index == aggrCnts_[ss]) index--;
      index = aggrLabels_[ss][index];
      cellsOccupied_[index] = -(cellsOccupied_[index] + 1);
      itmp = index;
      for (inputID = 0; inputID < nInputs_; inputID++)
      {
         jtmp = itmp % n1d_;
         itmp = itmp / n1d_;
         if ((randFlag & 1) == 1)
              dtmp = ((double) (jtmp + 0.999*PSUADE_drand())) / (double) n1d_;
         else dtmp = ((double) (jtmp + 0.5)) / (double) n1d_;
         sampleMatrix_[ss][inputID] = dtmp * ranges[inputID] +
                                            lbounds[inputID];
         if (sampleMatrix_[ss][inputID] < lowerBounds_[inputID])
            sampleMatrix_[ss][inputID] = lowerBounds_[inputID];
         if (sampleMatrix_[ss][inputID] > upperBounds_[inputID])
            sampleMatrix_[ss][inputID] = upperBounds_[inputID];
      }
   }
   nAggrs_ = currNAggr;

   if (changeInfoName_ == 0) fp = fopen("psuadeMetisInfo", "w");
   else                      fp = fopen("psuadeMetisInfo.tmp", "w");
   if (fp != NULL)
   {
      fprintf(fp, "%d %d %d\n", nAggrs_, nSamples_, nInputs_);
      for (ii = 0; ii < nAggrs_; ii++)
      {
         fprintf(fp, "%d\n", aggrCnts_[ii]);
         for (jj = 0; jj < aggrCnts_[ii]; jj++)
         {
            fprintf(fp, "%d ", aggrLabels_[ii][jj]);
            if (jj % 10 == 0) fprintf(fp, "\n");
         }
         fprintf(fp, "\n");
      }
      jj = 0;
      for (ii = 0; ii < graphN_; ii++) if (cellsOccupied_[ii] < 0) jj++;
      fprintf(fp, "%d\n", jj);
      for (ii = 0; ii < graphN_; ii++)
         if (cellsOccupied_[ii] < 0) 
            fprintf(fp, "%7d %7d\n",ii,-(cellsOccupied_[ii]+1));
      fclose(fp);
   }

   delete [] ranges;
   delete [] ubounds;
   delete [] lbounds;
   delete [] labels;
   delete [] localIA;
   delete [] localJA;
   if (refineArray != NULL) delete [] refineArray;

   if (printLevel_ > 0)
   {
      printEquals(PL_INFO, 0);
      printf("MetisSampling::refine(2): nSamples = %d\n", nSamples_);
      printf("MetisSampling::refine(2): nInputs  = %d\n", nInputs_);
      printf("MetisSampling::refine(2): nOutputs = %d\n", nOutputs_);
      if (randomize_ != 0)
           printf("MetisSampling::refine: randomize on\n");
      else printf("MetisSampling::refine: randomize off\n");
      printEquals(PL_INFO, 0);
   }
   return 0;
}

// ************************************************************************
// set internal scheme
// ------------------------------------------------------------------------
int MetisSampling::setParam(char * sparam)
{
   int  ii,curVol,count, *fGridFlags, sID, localN, *subLabels, sID2, rowInd;
   char  winput[501];
   FILE *fp;

   sscanf(sparam, "%s", winput);
   if (!strcmp(winput, "reset"))
   {
      if (changeInfoName_ == 0) fp = fopen("psuadeMetisInfo", "r");
      else                      fp = fopen("psuadeMetisInfo.tmp", "r");
      if (fp != NULL)
      {
         fclose(fp);
         if (changeInfoName_ == 0) unlink("psuadeMetisInfo");
         else                      unlink("psuadeMetisInfo.tmp");
      }
   }
   else if (!strcmp(winput, "changeInfoName"))
   {
      changeInfoName_ = 1;
   }
   else if (!strcmp(winput, "setUniformRefinement"))
   {
      refineType_ = 0;
   }
   else if (!strcmp(winput, "setAdaptiveRefinementBasedOnErrors"))
   {
      refineType_ = 1;
   }
   else if (!strcmp(winput, "setAdaptiveRefinementBasedOnOutputs"))
   {
      refineType_ = 2;
   }
   else if (!strcmp(winput, "calVolumes"))
   {
      curVol = count = 0;
      for (ii = 0; ii < nAggrs_; ii++)
      {
         if (sampleOutput_[ii] == 1)
         {
            curVol += aggrCnts_[ii];
            count++;
         }
      }
      return curVol;
   }
   else if (!strcmp(winput, "totalVolumes"))
   {
      curVol = 0;
      for (ii = 0; ii < nAggrs_; ii++) curVol += aggrCnts_[ii];
      return curVol;
   }
   else if (!strcmp(winput, "setRefineSize"))
   {
      sscanf(sparam, "%s %d", winput, &refineSize_);
   }
   else if (!strcmp(winput, "genMeshPlot") && nInputs_ == 2)
   {
      fGridFlags = new int[graphN_];
      for (sID = 0; sID < nAggrs_; sID++)
      {
         localN = aggrCnts_[sID];
         subLabels = aggrLabels_[sID];
         for (ii = 0; ii < localN; ii++)
         {
            rowInd = subLabels[ii];
            if (sampleOutput_[sID] == 1) fGridFlags[rowInd] = 1;
            else                         fGridFlags[rowInd] = 0;
         }
      }
      fp = fopen("metisMeshPlot.m", "w");
      // add check for NULL by Bill Oliver
      if(fp != NULL)
      {
         fprintf(fp, "meshA = [ \n");
         for (sID = 0; sID < n1d_; sID++)
         {
            for (sID2 = 0; sID2 < n1d_; sID2++)
               fprintf(fp, "%d ", fGridFlags[sID*n1d_+sID2]);
            fprintf(fp, "\n");
         }
         fprintf(fp, "];\n");
         fprintf(fp, "contour(meshA)\n");
         fclose(fp);
      }
      // Cleanup by Bill Oliver
      delete [] fGridFlags;
      return 0;
   }
   return 0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
MetisSampling& MetisSampling::operator=(const MetisSampling & ms)
{
   if (this == &ms) return *this;
   refineType_ = ms.refineType_;
   refineSize_ = ms.refineSize_;
   n1d_ = ms.n1d_;
   nAggrs_ = ms.nAggrs_;
   graphN_ = ms.graphN_;
 
   aggrLabels_ = new int*[nAggrs_]; 
   aggrCnts_ = new int[nAggrs_];
   for(int i = 0; i < nAggrs_; i++)
   {
      aggrCnts_[i] = ms.aggrCnts_[i];
      for(int j = 0; j < aggrCnts_[j]; j++)
      {
         aggrLabels_[j] = new int[aggrCnts_[i]];
         aggrLabels_[i][j] = ms.aggrLabels_[i][j];
      }
   }
   
   // MetisSampling inherits from Sampling so include the parent 
   // class data members
   printLevel_ = ms.printLevel_;
   samplingID_ = ms.samplingID_;
   nSamples_ = ms.nSamples_;
   nInputs_ = ms.nInputs_;
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

   graphI_ = new int[graphN_+1];
   for(int i = 0; i <= graphN_; i++)
      graphI_[i] = ms.graphI_[i];

   graphJ_ = new int[graphN_*nInputs_*2 + 1];
   for(int i = 0; i <= graphN_*nInputs_*2; i++)
      graphJ_[i] = ms.graphJ_[i];

   cellsOccupied_ = new int[graphN_];
   for(int i = 0; i < graphN_; i++)
      cellsOccupied_[i] = ms.cellsOccupied_[i];
   return (*this);
}

