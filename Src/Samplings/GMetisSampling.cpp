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
// Functions for the GMetis sampling class (generalized Metis) 
// AUTHOR : CHARLES TONG
// DATE   : 2010
// ************************************************************************
#include <stdio.h>
#include <assert.h>
#include <sstream>
#include <unistd.h>
using namespace std;

#include "Psuade.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "FuncApprox.h"
#include "Mars.h"
#include "GMetisSampling.h"

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
GMetisSampling::GMetisSampling() : Sampling()
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
}

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
GMetisSampling::GMetisSampling(const GMetisSampling & gms) : Sampling()
{
   refineType_ = gms.refineType_;
   refineSize_ = gms.refineSize_;
   n1d_ = gms.n1d_;
   nAggrs_ = gms.nAggrs_;
   graphN_ = gms.graphN_;
 
   aggrLabels_ = new int*[nAggrs_]; 
   aggrCnts_ = new int[nAggrs_];
   for(int i = 0; i < nAggrs_; i++)
   {
      aggrCnts_[i] = gms.aggrCnts_[i];
      for(int j = 0; j < aggrCnts_[j]; j++)
      {
         aggrLabels_[j] = new int[aggrCnts_[i]];
         aggrLabels_[i][j] = gms.aggrLabels_[i][j];
      }
   }
   
   // MetisSampling inherits from Sampling so include the parent 
   // class data members
   printLevel_ = gms.printLevel_;
   samplingID_ = gms.samplingID_;
   nSamples_ = gms.nSamples_;
   nInputs_ = gms.nInputs_;
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

   graphI_ = new int[graphN_+1];
   for(int i = 0; i <= graphN_; i++)
      graphI_[i] = gms.graphI_[i];

   graphJ_ = new int[graphN_*nInputs_*2 + 1];
   for(int i = 0; i <= graphN_*nInputs_*2; i++)
      graphJ_[i] = gms.graphJ_[i];

   cellsOccupied_ = new int[graphN_];
   for(int i = 0; i < graphN_; i++)
      cellsOccupied_[i] = gms.cellsOccupied_[i];

}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
GMetisSampling::~GMetisSampling()
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
int GMetisSampling::initialize(int initLevel)
{
   int    *incrs, inputID, ii, jj, itmp, jtmp, nnz, sampleID;
   int    options[10], index, count, kk;
#ifdef HAVE_METIS
   int    wgtflag=0, numflag=0, edgeCut=0;
#endif
   double *ranges=NULL, dtmp, *lbounds, *ubounds, expand=0.0;
   char   response[500], pString[500];
   FILE   *fp;

   if (nInputs_ == 0 || lowerBounds_ == NULL || upperBounds_ == NULL)
   {
      printf("GMetisSampling::initialize ERROR - input not set up.\n");
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
   if (nInputs_ > 21)
   {
      printf("GMetisSampling ERROR : nInputs > 21 currently not supported.\n");
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

   fp = fopen("psuadeGMetisInfo", "r");
   if (fp != NULL)
   {
      printf("INFO: psuadGeMetisInfo file found. Reading it in ...\n");
      fscanf(fp, "%d", &jj);
      printf("GMetisSampling Info: a partition file is found.\n");
      printf("          The file name is psuadeGMetisInfo.\n");
      printf("          It may have been left behind by previous\n");
      printf("          call to GMetis (It has %d subdomains).\n",jj);
      sprintf(pString,"Reuse the file ? (y or n) ");
      getString(pString, response);

      if (response[0] == 'y')
      {
         printf("GMetisSampling Info: reconstructing the partitioning.\n");
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
               cellsOccupied_[aggrLabels_[ii][jj]] = ii;
            }
         }
         fscanf(fp, "%d", &count);
         for (ii = 0; ii < count; ii++)
         {
            fscanf(fp, "%d %d", &jj, &kk);
            cellsOccupied_[jj] = - kk - 1;
         }
      }
   }

   if (aggrCnts_ == NULL)
   {
      options[0] = 0;
#ifdef HAVE_METIS
      if (nSamples_ <= 0)
      {
         strcpy(pString, "Enter the number of sample points: ");
         nSamples_ = getInt(1, 1000000, pString);
         strcpy(pString, "Enter the number of partitions: ");
         nAggrs_ = getInt(1, nSamples_, pString);
      }
      else nAggrs_ = nSamples_;
      if (printLevel_ > 1)
         printf("GMetisSampling: creating %d partitions...\n", nAggrs_);
      METIS_PartGraphRecursive(&graphN_, graphI_, graphJ_, NULL, NULL,
                &wgtflag,&numflag,&nAggrs_,options,&edgeCut,cellsOccupied_);
#else
      printf("GMetisSampling ERROR : METIS not installed.\n");
      exit(1);
#endif
      if (printLevel_ > 1)
         printf("GMetisSampling:: %d subdomains created.\n", nAggrs_);

      aggrCnts_ = new int[nAggrs_];
      for (ii = 0; ii < nAggrs_; ii++) aggrCnts_[ii] = 0;
      for (ii = 0; ii < graphN_; ii++) aggrCnts_[cellsOccupied_[ii]]++;  
      aggrLabels_ = new int*[nAggrs_];
      for (ii = 0; ii < nAggrs_; ii++)
      {
         aggrLabels_[ii] = new int[aggrCnts_[ii]];
         aggrCnts_[ii] = 0;
      }
      for (ii = 0; ii < graphN_; ii++)
      {
         index = cellsOccupied_[ii];
         aggrLabels_[index][aggrCnts_[index]++] = ii;  
      }
   }
  
   if(fp != NULL)  fclose(fp);
   fp = fopen("psuadeGMetisInfo", "w");
   if (fp != NULL)
   {
      fprintf(fp, "%d\n", nAggrs_);
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
   if (initLevel != 0) return 0;

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

   sampleID = 0;
   while (sampleID < nSamples_)
   {
      index = (int) (PSUADE_drand() * aggrCnts_[sampleID%nAggrs_]);
      if (index == aggrCnts_[sampleID%nAggrs_]) index--;
      index = aggrLabels_[sampleID%nAggrs_][index];
      cellsOccupied_[index] = -(cellsOccupied_[index] + 1);
      itmp = index;
      for (inputID = 0; inputID < nInputs_; inputID++)
      {
         jtmp = itmp % n1d_;
         itmp = itmp / n1d_;
         dtmp = ((double) (jtmp + 0.999*PSUADE_drand())) / (double) n1d_;
         sampleMatrix_[sampleID][inputID] = dtmp * ranges[inputID] +
                                            lbounds[inputID];
      }
      sampleID++;
   }

   if (printLevel_ > 4)
   {
      printf("GMetisSampling::initialize: nSamples = %d\n", nSamples_);
      printf("GMetisSampling::initialize: nInputs  = %d\n", nInputs_);
      printf("GMetisSampling::initialize: nOutputs = %d\n", nOutputs_);
      if (randomize_ != 0)
           printf("GMetisSampling::initialize: randomize on\n");
      else printf("GMetisSampling::initialize: randomize off\n");
      for (inputID = 0; inputID < nInputs_; inputID++)
         printf("    GMetisSampling input %3d = [%e %e]\n", inputID+1,
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
int GMetisSampling::refine(int nLevels, int randFlag, double threshold,
                           int nSamples, double *sampleErrors)
{
   int    inputID, ii, jj, count, itmp, jtmp, localN, *localIA, *localJA;
   int    maxN, rowInd, colInd, localNNZ, *subLabels, index, count0;
   int    options[10], count1;
#ifdef HAVE_METIS
   int    wgtflag=0, numflag=0, edgeCut=0, itwo=2;
#endif
   int    *tmpCnts, **tmpLabels, currNAggr, newCell, oldNumSamples;
   int    *oldSampleStates, *labels, status, outputID, ss, maxNNZ;
   int    *refineArray=NULL, *node2Aggr=NULL, splitCount;
   int    splitSuccess, cellCnt=0;
   double *ranges, dtmp, **oldSampleMatrix, *oldSampleOutput;
   double *lbounds, *ubounds, *sortList2;
   double *sortList1, *aggrErrs=NULL;
   FILE   *fp;

   if (cellsOccupied_ == NULL || aggrLabels_ == NULL)
   {
      printf("GMetisSampling::refine ERROR - need to call initialize first.\n");
      exit(1);
   }
   ranges  = new double[nInputs_];
   lbounds = lowerBounds_;
   ubounds = upperBounds_;
   for (inputID = 0; inputID < nInputs_; inputID++)
   {
      ranges[inputID] = ubounds[inputID] - lbounds[inputID];
      if (ranges[inputID] <= 0.0)
      {
         printf("GMetisSampling::refine ERROR - lbound/ubound mismatch.\n");
         exit(1);
      }
   }

   for (ss = 0; ss < nSamples_; ss++)
   {
      itmp = 0;
      for (inputID = nInputs_-1; inputID >= 0; inputID--)
      {
         itmp = itmp * n1d_;
         dtmp = sampleMatrix_[ss][inputID];
         dtmp = (dtmp - lowerBounds_[inputID]) / ranges[inputID];
         jtmp = (int) (dtmp * n1d_);
         itmp += jtmp;
      }
      cellsOccupied_[itmp] = -(cellsOccupied_[itmp] + 1); 
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
   localIA = new int[maxN+1];
   assert(localIA != NULL);
   maxNNZ = maxN * (2 * nInputs_ + 1);
   localJA = new int[maxNNZ];
   assert(localJA != NULL);

   refineArray = NULL;
   aggrErrs = NULL;

   if (refineType_ == 1 && sampleErrors == NULL)
   {
      printf("GMetisSampling::refine ERROR- error based but no error given.\n");
      exit(1);
   }

   if (refineType_ == 1 && sampleErrors != NULL)
   {
      printf("GMetisSampling::refine - maximum number of new points = %d.\n",
             refineSize_);
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
               printf("GMetisSampling::refine ERROR (index out of bound.)\n");
               printf("               rowInd = %d (should be in [0,%d]\n",
                      rowInd, graphN_-1);
               printf("               Aggregate = %d \n", ii);
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
            printf("GMetisSampling::refine ERROR - node2Aggr not correct.\n");
            printf("                Please consult PSUADE developers.\n");
            exit(1);
         }
      }

      aggrErrs = new double[nAggrs_];
      for (ii = 0; ii < nAggrs_; ii++)
          aggrErrs[ii] = PABS(sampleErrors[ii]);

      sortList1 = new double[nAggrs_];
      sortList2 = new double[nAggrs_];
      for (ii = 0; ii < nAggrs_; ii++) sortList1[ii] = PABS(aggrErrs[ii]);
      for (ii = 0; ii < nAggrs_; ii++) sortList2[ii] = (double) ii;
      sortDbleList2(nAggrs_, sortList1, sortList2);
      refineArray = new int[nAggrs_];
      for (ii = nAggrs_-1; ii >= 0; ii--) refineArray[ii] = 0;
      cellCnt = 0;

      for (ii = nAggrs_-1; ii >= 0; ii--)
      {
         index = (int) sortList2[ii];
         localN = aggrCnts_[index];
         aggrErrs[index] *= pow(1.0*localN,1.0/nInputs_);
         if (aggrErrs[index] != 0.0) 
         {
            //printf("GMetisSampling::refine - chosen aggregate, error = %e\n",
            //        aggrErrs[index]);
            if (refineArray[index] == 0)
            {
               cellCnt++;
               refineArray[index] = 1;
               if (cellCnt >= refineSize_) break;
            }
            if (cellCnt >= refineSize_) break;
         }
         if (cellCnt >= refineSize_) break;
      }
      delete [] sortList1;
      delete [] sortList2;
      delete [] node2Aggr;
   }

   options[0] = 0;
   currNAggr = nAggrs_;
   splitCount = splitSuccess = 0;
   for (ss = 0; ss < nAggrs_; ss++)
   {

      localN = aggrCnts_[ss];
      if (localN > maxN) 
      {
         printf("GMetisSampling INTERNAL ERROR : consult PSUADE developers\n");
         exit(1);
      }
      if (refineArray == NULL || refineArray[ss] == 1) splitCount++;
      if (localN == 1 && (refineArray == NULL || refineArray[ss] == 1))
         printf("GMetisSampling::refine INFO- cannot split cell (too small).\n");
      if (localN > 1 && (refineArray == NULL || refineArray[ss] == 1))
      {
         if (printLevel_ > 4)
            printf("GMetisSampling::refine - split sample %d\n",ss+1);
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
            localIA[ii+1] = localNNZ;
         }
         if (localNNZ > maxNNZ) 
         {
           printf("GMetisSampling INTERNAL ERROR: consult PSUADE developers\n");
           exit(1);
         }

#ifdef HAVE_METIS
         METIS_PartGraphRecursive(&localN, localIA, localJA, NULL, NULL,
                &wgtflag,&numflag,&itwo,options,&edgeCut,labels);
#else
         printf("GMetisSampling ERROR : METIS not installed.\n");
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
            aggrCnts_[ss] = count1; 
            count = 0;
            for (ii = 0; ii < localN; ii++) 
               if (labels[ii] == 1) 
                  aggrLabels_[ss][count++] = aggrLabels_[ss][ii];
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
            aggrCnts_[ss] = count0; 
            count = 0;
            for (ii = 0; ii < localN; ii++) 
               if (labels[ii] == 0) 
                  aggrLabels_[ss][count++] = aggrLabels_[ss][ii];
         }
         currNAggr++;
      }
   }
           
   if (printLevel_ > 4 && splitSuccess != splitCount)
      printf("GMetisSampling:: number of successful splits = %d (out of %d)\n",
             splitSuccess, splitCount);

   oldNumSamples   = nSamples_;
   oldSampleMatrix = sampleMatrix_;
   oldSampleOutput = sampleOutput_;
   oldSampleStates = sampleStates_;
   nSamples_ = nSamples_ + currNAggr - nAggrs_;
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

   for (ss = nAggrs_; ss < currNAggr; ss++)
   {
      index = (int) (PSUADE_drand() * aggrCnts_[ss]);
      if (index == aggrCnts_[ss]) index--;
      index = aggrLabels_[ss][index];
      cellsOccupied_[index] = -(cellsOccupied_[index] + 1);
      itmp = index;
      for (inputID = 0; inputID < nInputs_; inputID++)
      {
         jtmp = itmp % n1d_;
         itmp = itmp / n1d_;
         dtmp = ((double) (jtmp + 0.999*PSUADE_drand())) / (double) n1d_;
         sampleMatrix_[oldNumSamples][inputID] = dtmp * ranges[inputID] +
                                                 lbounds[inputID];
      }
      oldNumSamples++;
   }
   nAggrs_ = currNAggr;

   fp = fopen("psuadeGMetisInfo", "w");
   if (fp != NULL)
   {
      fprintf(fp, "%d\n", nAggrs_);
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
   delete [] labels;
   delete [] localIA;
   delete [] localJA;
   if (refineArray != NULL) delete [] refineArray;
   if (aggrErrs    != NULL) delete [] aggrErrs;

   if (printLevel_ > 4)
   {
      printf("GMetisSampling::refine: nSamples = %d\n", nSamples_);
      printf("GMetisSampling::refine: nInputs  = %d\n", nInputs_);
      printf("GMetisSampling::refine: nOutputs = %d\n", nOutputs_);
      for (inputID = 0; inputID < nInputs_; inputID++)
         printf("    GMetisSampling input %3d = [%e %e]\n", inputID+1,
                lowerBounds_[inputID], upperBounds_[inputID]);
   }
   return 0;
}

// ************************************************************************
// set internal scheme
// ------------------------------------------------------------------------
int GMetisSampling::setParam(string sparam)
{
   int           pos, ii, curVol, count;
   istringstream buffer;
   string        substr;
   FILE *fp;

   pos = sparam.find("reset");
   if (pos >= 0)
   {
      fp = fopen("psuadeGMetisInfo", "r");
      if (fp != NULL)
      {
         fclose(fp);
         unlink("psuadeGMetisInfo");
      }
      return 0;
   }

   pos = sparam.find("setUniformRefinement");
   if (pos >= 0)
   {
      refineType_ = 0;
      return 0;
   }

   pos = sparam.find("setAdaptiveRefinementBasedOnErrors");
   if (pos >= 0)
   {
      refineType_ = 1;
      return 0;
   }

   pos = sparam.find("calVolumes");
   if (pos >= 0)
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

   pos = sparam.find("totalVolumes");
   if (pos >= 0)
   {
      curVol = 0;
      for (ii = 0; ii < nAggrs_; ii++) curVol += aggrCnts_[ii];
      return curVol;
   }

   pos = sparam.find("setRefineSize");
   if (pos >= 0)
   {
      substr = sparam.substr(14);
      buffer.str(substr);
      buffer >> refineSize_;
      return 0;
   }
   printf("GMetisSampling ERROR:: setParam - invalid param.\n");
   return -1;
}

// ************************************************************************
// equal operator  Modified by Bill Oliver
// ------------------------------------------------------------------------
GMetisSampling& GMetisSampling::operator=(const GMetisSampling & gms)
{
   if(this == &gms) return *this;
   refineType_ = gms.refineType_;
   refineSize_ = gms.refineSize_;
   n1d_ = gms.n1d_;
   nAggrs_ = gms.nAggrs_;
   graphN_ = gms.graphN_;
 
   aggrLabels_ = new int*[nAggrs_]; 
   aggrCnts_ = new int[nAggrs_];
   for(int i = 0; i < nAggrs_; i++)
   {
      aggrCnts_[i] = gms.aggrCnts_[i];
      for(int j = 0; j < aggrCnts_[j]; j++)
      {
         aggrLabels_[j] = new int[aggrCnts_[i]];
         aggrLabels_[i][j] = gms.aggrLabels_[i][j];
      }
   }
   
   // MetisSampling inherits from Sampling so include the parent 
   // class data members
   printLevel_ = gms.printLevel_;
   samplingID_ = gms.samplingID_;
   nSamples_ = gms.nSamples_;
   nInputs_ = gms.nInputs_;
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

   graphI_ = new int[graphN_+1];
   for(int i = 0; i <= graphN_; i++)
      graphI_[i] = gms.graphI_[i];

   graphJ_ = new int[graphN_*nInputs_*2 + 1];
   for(int i = 0; i <= graphN_*nInputs_*2; i++)
      graphJ_[i] = gms.graphJ_[i];

   cellsOccupied_ = new int[graphN_];
   for(int i = 0; i < graphN_; i++)
      cellsOccupied_[i] = gms.cellsOccupied_[i];
   return (*this);
}

