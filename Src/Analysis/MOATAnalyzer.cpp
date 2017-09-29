// ************************************************************************
// Copyright (c) 2007   Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory.
// Written by the PSUADE team. 
// All rights reserved.
//
// Please see the COPYRIGHT and LICENSE file for the copyright notice,
// disclaimer, contact information and the GNU Lesser General Public License.
//
// PSUADE is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free 
// Software Foundation) version 2.1 dated February 1999.
//
// PSUADE is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the terms and conditions of the GNU Lesser
// General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
// ************************************************************************
// Functions for the class MOATAnalyzer  
// AUTHOR : CHARLES TONG
// DATE   : 2004
// ************************************************************************
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>

#include "MOATConstraints.h"
#include "PsuadeUtil.h"
#include "sysdef.h"
#include "MOATAnalyzer.h"
#include "BinomialAnalyzer.h"
#include "BootstrapAnalyzer.h"
#include "Psuade.h"
#include "PrintingTS.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
MOATAnalyzer::MOATAnalyzer() : Analyzer(), nInputs_(0), nOutputs_(0), 
                  nSamples_(0), outputID_(0), analysisError_(0), means_(0), 
                  stds_(0), modifiedMeans_(0), modifiedStds_(0),
                  indexesSortedByModifiedMeans_(0), 
                  sigma1LowerBounds_(0), sigma1UpperBounds_(0), 
                  sigma2LowerBounds_(0), sigma2UpperBounds_(0)
{
   setName("MOAT");
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
MOATAnalyzer::~MOATAnalyzer()
{
   if (means_) delete [] means_;
   if (modifiedMeans_) delete [] modifiedMeans_;
   if (stds_) delete [] stds_;
   if (modifiedStds_) delete [] modifiedStds_;
   if (indexesSortedByModifiedMeans_) delete [] indexesSortedByModifiedMeans_;
   if (sigma1LowerBounds_) delete [] sigma1LowerBounds_;
   if (sigma1UpperBounds_) delete [] sigma1UpperBounds_;
   if (sigma2LowerBounds_) delete [] sigma2LowerBounds_;
   if (sigma2UpperBounds_) delete [] sigma2UpperBounds_;
}

// ************************************************************************
// perform analysis
// ------------------------------------------------------------------------
double MOATAnalyzer::analyze(aData &adata)
{
   int    *S, ii, jj, diffIndex;
   int    iD, *indexTrack, index, n1, n2, flag, actualPaths, sigFlag;
   int    diffCnt, *counts, itemp, iaCnt, iD2, nPaths, ip, printLevel;
   double *X, *Y, *YY, *YG, *XG, *Xbase;
   double xtemp1, xtemp2;
   double ytemp1, ytemp2, scale, *xLower, *xUpper, dtemp, dtemp2, thresh;
   double *XSort, *YSort, dsum, dstdev, *YB, binMax=-1.0e35, binMin=1.0e35;
   char   winput[500];
   char   pString[500];
   PsuadeData *ioPtr;
   BinomialAnalyzer  binAnalyzer;
   BootstrapAnalyzer bsAnalyzer;
   aData             bData;
   pData             qData;
   MOATConstraints   *constrPtr;

   // ---------------------------------------------------------------
   // clean up first 
   // ---------------------------------------------------------------
   if (means_) delete [] means_;
   if (modifiedMeans_) delete [] modifiedMeans_;
   if (stds_) delete [] stds_;
   if (modifiedStds_) delete [] modifiedStds_;
   if (indexesSortedByModifiedMeans_) delete [] indexesSortedByModifiedMeans_;
   if (sigma1LowerBounds_) delete [] sigma1LowerBounds_;
   if (sigma1UpperBounds_) delete [] sigma1UpperBounds_;
   if (sigma2LowerBounds_) delete [] sigma2LowerBounds_;
   if (sigma2UpperBounds_) delete [] sigma2UpperBounds_;

   // ---------------------------------------------------------------
   // extract test information and data
   // ---------------------------------------------------------------
   printLevel = adata.printLevel_;
   nInputs_   = adata.nInputs_;
   nOutputs_  = adata.nOutputs_;
   nSamples_  = adata.nSamples_;
   xLower     = adata.iLowerB_;
   xUpper     = adata.iUpperB_;
   X          = adata.sampleInputs_;
   Y          = adata.sampleOutputs_;
   S          = adata.sampleStates_;
   outputID_  = adata.outputID_;
   ioPtr      = adata.ioPtr_;
   if (adata.inputPDFs_ != NULL)
   {
      n1 = 0;
      for (ii = 0; ii < nInputs_; ii++) n1 += adata.inputPDFs_[ii];
      if (n1 > 0)
      {
         printOutTS(1, "MOATAnalysis INFO: some inputs have non-uniform PDFs,\n");
         printOutTS(1, "             but they are not relevant in this analysis\n");
         printOutTS(1, "             (the sample should have been generated\n");
         printOutTS(1, "              with the desired distributions.)\n");
      }
   }
   if (ioPtr != NULL) ioPtr->getParameter("input_names", qData);

   if (nInputs_ <= 0 || nSamples_ <= 0 || nOutputs_ <= 0 || 
       outputID_ < 0 || outputID_ >= nOutputs_)
   {
      printOutTS(0, "MOATAnalysis ERROR: invalid arguments.\n");
      printOutTS(0, "    nInputs  = %d\n", nInputs_);
      printOutTS(0, "    nOutputs = %d\n", nOutputs_);
      printOutTS(0, "    nSamples = %d\n", nSamples_);
      printOutTS(0, "    outputID = %d\n", outputID_+1);
      exit(1);
   } 
   n1 = 0;
   for (iD = 0; iD < nSamples_; iD++)
      if (S[iD] != 1 || Y[iD*nOutputs_+outputID_] == PSUADE_UNDEFINED) n1++; 
   if (n1 > 0) printOutTS(1, "MOATAnalysis INFO: %d invalid data points.\n",n1);

   XSort = new double[nSamples_];
   for (ii = 0; ii < nInputs_; ii++)
   {
      for (iD = 0; iD < nSamples_; iD++)
         XSort[iD] = X[iD*nInputs_+ii];
      sortDbleList(nSamples_, XSort);
      dtemp = XSort[nSamples_-1] - XSort[0];
      dtemp2 = xUpper[ii] - xLower[ii];
      if (PABS(dtemp-dtemp2) > 1.0e-6)
      {
         printOutTS(1, "MOATAnalyzer WARNING: input and data range mismatch.\n");
         printOutTS(1, "    (This can be the result of applying constraints.)\n");
         printOutTS(1, "    diagnostics: \n");
         printOutTS(1, "    Input %3d: original vs new ranges = %e %e\n", ii+1,
                dtemp2, dtemp);
      }
   }
   delete [] XSort;

   printAsterisks(0);
   for (ii = 0; ii < 23; ii++) printOutTS(2, "*");
   printOutTS(2, " MOAT Analysis ");
   for (ii = 0; ii < 23; ii++) printOutTS(2, "*");
   printOutTS(2, "\n");
   printDashes(0);

   constrPtr = new MOATConstraints();
   if(ioPtr != NULL) constrPtr->initialize(ioPtr);

   YY = new double[nSamples_];
   YG = new double[nSamples_];
   XG = new double[nSamples_];
   for (iD = 0; iD < nSamples_; iD++) YY[iD] = Y[nOutputs_*iD+outputID_];
   counts = new int[nInputs_];
   for (ii = 0; ii < nInputs_; ii++) counts[ii] = 0;
   means_ = new double[nInputs_];
   modifiedMeans_ = new double[nInputs_];
   indexesSortedByModifiedMeans_ = new int[nInputs_];
   stds_ = new double[nInputs_];
   modifiedStds_ = new double[nInputs_];
   for (ii = 0; ii < nInputs_; ii++)
      means_[ii] = modifiedMeans_[ii] = stds_[ii] = modifiedStds_[ii] = 0.0;
   indexTrack = new int[nSamples_];
   for (iD = 0; iD < nSamples_; iD++) indexTrack[iD] = -1;
   Xbase = new double[nSamples_];
   for (iD = 0; iD < nSamples_; iD++) Xbase[iD] = 0.0;

   indexTrack[0] = -1;
   for (iD = 1; iD < nSamples_; iD++)
   {
      if (printLevel > 3 && (iD % 10 == 0))
         printOutTS(2, "MOATAnalysis: processing sample %d\n", iD+1);

      Xbase[iD] = 0.0;
      ytemp1 = YY[iD-1]; 
      ytemp2 = YY[iD]; 
      diffCnt = 0;
      for (ii = 0; ii < nInputs_; ii++)
      {
         xtemp1 = X[(iD-1)*nInputs_+ii]; 
         xtemp2 = X[iD*nInputs_+ii]; 
         if (xtemp1 != xtemp2 && ytemp1 !=  PSUADE_UNDEFINED &&
             ytemp2 != PSUADE_UNDEFINED && S[iD-1] != 0 &&
             S[iD] != 0)
         {
            diffCnt++;
            diffIndex = ii;
         }
      }
      if (diffCnt == 1)
      {
         indexTrack[iD] = diffIndex;
         xtemp1 = X[(iD-1)*nInputs_+diffIndex]; 
         xtemp2 = X[iD*nInputs_+diffIndex]; 
         scale = constrPtr->getScale(&X[iD*nInputs_],diffIndex,flag);
         if (flag == 1) scale = xUpper[diffIndex] - xLower[diffIndex];
         else if (printLevel > 3)
         {
            printOutTS(2, "MOATAnalysis: sample group %d, ", iD+1);
            printOutTS(2, "input %d, scale = %e\n", diffIndex, scale);
         }
         YG[iD] = (ytemp2-ytemp1)/(xtemp2-xtemp1)*scale;
         if (xtemp2 > xtemp1) XG[iD] = xtemp2;
         else                 XG[iD] = xtemp1;
         counts[diffIndex]++;
         if (xtemp2 > xtemp1) Xbase[iD] = xtemp1;
         else                 Xbase[iD] = xtemp2;
      }
      else
      {
         YG[iD] = PSUADE_UNDEFINED;
         indexTrack[iD] = -1;
      }
   }

   if (nSamples_ / (nInputs_+1) * (nInputs_+1) == nSamples_)
   {
      for (iD = 0; iD < nSamples_; iD+=(nInputs_+1))
         indexTrack[iD] = -1;
   }

   for (iD = 0; iD < nSamples_; iD++)
   {
      if (YG[iD] != PSUADE_UNDEFINED)
      {
         index = indexTrack[iD];
         if (index >= 0)
         {
            means_[index] += YG[iD];
            modifiedMeans_[index] += PABS(YG[iD]);
         }
      }
   }
   for (ii = 0; ii < nInputs_; ii++)
   {
      if (counts[ii] > 0)
      {
         means_[ii] /= (double) counts[ii];
         modifiedMeans_[ii] /= (double) counts[ii];
      }
      else 
      {
         printOutTS(2, "MOATAnalysis analyze: zero data points for input %d\n",
                ii+1);
         means_[ii] = 0.0;
         modifiedMeans_[ii] = 0.0;
      }
   }
   for (iD = 0; iD < nSamples_; iD++)
   {
      if (YG[iD] != PSUADE_UNDEFINED)
      {
         index = indexTrack[iD];
         if (index >= 0)
         {
            stds_[index] += (YG[iD] - means_[index])*(YG[iD] - means_[index]);
            modifiedStds_[index] += (PABS(YG[iD]) - modifiedMeans_[index])*
                                    (PABS(YG[iD]) - modifiedMeans_[index]);
         }
      }
   }
   for (ii = 0; ii < nInputs_; ii++)
   {
      if (counts[ii] > 1)
      {
         stds_[ii] /= (double) (counts[ii] - 1);
         modifiedStds_[ii] /= (double) (counts[ii] - 1);
      }
      else
      {
         printOutTS(2, "MOATAnalysis analyze: %d data points for input %d\n",
                counts[ii], ii+1);
         stds_[ii] = 0.0;
         modifiedStds_[ii] = 0.0;
      }
      if (stds_[ii] < 0.0) stds_[ii] = -sqrt(-stds_[ii]);
      else                 stds_[ii] = sqrt(stds_[ii]);
      if (modifiedStds_[ii] < 0.0) modifiedStds_[ii] = -sqrt(-modifiedStds_[ii]);
      else                         modifiedStds_[ii] = sqrt(modifiedStds_[ii]);
   }
   n1 = 0;
   for (ii = 0; ii < nInputs_; ii++) n1 += counts[ii];
   if (n1 <= 0)
   {
      delete [] counts;
      delete constrPtr;
      delete [] YY;
      delete [] YG;
      delete [] XG;
      delete [] indexTrack;
      delete [] Xbase;
      printOutTS(1, "MOATAnalysis INFO: probably not an MOAT sample.\n");
      return 1;
   }

   printDashes(0);
   for (ii = 0; ii < nInputs_; ii++)
      printOutTS(2, "Input %3d (unmod. mean & std) = %12.4e %12.4e \n",
             ii+1, means_[ii], stds_[ii]);

   printAsterisks(0);
   for (ii = 0; ii < 23; ii++) printOutTS(2, "*");
   printOutTS(2, " MOAT Analysis ");
   for (ii = 0; ii < 23; ii++) printOutTS(2, "*");
   printOutTS(2, "\n");
   printDashes(0);
   for (ii = 0; ii < nInputs_; ii++)
      printOutTS(2, "Input %3d (mod. mean & std) = %12.4e %12.4e \n",
             ii+1, modifiedMeans_[ii], stds_[ii]);

   // ---------------------------------------------------------------
   // create matlab files
   // ---------------------------------------------------------------
   if (psAnalysisInteractive_ == 1 || psAnaExpertMode_ == 1)
      createScreenDiagramFile(nSamples_, nInputs_, YG, indexTrack, 
                        modifiedMeans_,stds_,outputID_,qData.strArray_);
   if (psAnalysisInteractive_ == 1 || psAnaExpertMode_ == 1)
      createScatterFile(nSamples_,nInputs_,YG,Xbase,indexTrack,
                        qData.strArray_);

   if (psAnalysisInteractive_ == 1 || psAnaExpertMode_ == 1)
      createBootstrapFile(nSamples_,nInputs_,YG,Xbase,indexTrack,
                        qData.strArray_);

   // ---------------------------------------------------------------
   // sorted list 
   // ---------------------------------------------------------------
   double indexes[nInputs_];
   double sortedModifiedMeans[nInputs_];
   std::copy(modifiedMeans_, modifiedMeans_+nInputs_, sortedModifiedMeans);
   for (ii = 0; ii < nInputs_; ii++) indexes[ii] = ii;
   sortDbleList2(nInputs_, sortedModifiedMeans, indexes);
   for (ii = 0; ii < nInputs_; ii++) 
      indexesSortedByModifiedMeans_[ii] = (int) indexes[ii];
   printAsterisks(0);
   for (ii = 0; ii < 18; ii++) printOutTS(2, "*");
   printOutTS(2, " MOAT Analysis (ordered) ");
   for (ii = 0; ii < 18; ii++) printOutTS(2, "*");
   printOutTS(2, "\n");
   printDashes(0);
   for (ii = nInputs_-1; ii >= 0; ii--)
   {
     iD =  indexesSortedByModifiedMeans_[ii];
      itemp = counts[iD] - 1;
      printOutTS(2, "%6d: Input %3d (mu*, sigma, dof) = %12.4e %12.4e %d\n",
             nInputs_-ii,iD+1, sortedModifiedMeans[ii], stds_[iD], itemp);
   }
   printAsterisks(0);
   if (printLevel > 1)
   {
      sigma1LowerBounds_ = new double[nInputs_];
      sigma1UpperBounds_ = new double[nInputs_];
      for (ii = 0; ii < 11; ii++) printOutTS(2, "-");
      printOutTS(2, " MOAT Analysis (ordered) : +- 1 sigma  ");
      for (ii = 0; ii < 11; ii++) printOutTS(2, "-");
      printOutTS(2, "\n");
      printDashes(0);
      for (ii = nInputs_-1; ii >= 0; ii--)
      {
         iD = indexesSortedByModifiedMeans_[ii];
         if (counts[iD] > 1) dtemp = sqrt((double) counts[iD] - 1.0);
         else                dtemp = 1.0;

         sigma1LowerBounds_[ii] = sortedModifiedMeans[ii]-modifiedStds_[iD]/dtemp;
         sigma1UpperBounds_[ii] = sortedModifiedMeans[ii]+modifiedStds_[iD]/dtemp;
         printOutTS(2, "Input %3d bounds = %12.4e %12.4e\n",
                iD+1, sigma1LowerBounds_[ii], sigma1UpperBounds_[ii]);
         if (printLevel > 2 && ii > 0)
         {
            iD2 = indexesSortedByModifiedMeans_[ii-1];
            if (counts[iD2-1] > 1) dtemp2 = sqrt((double) counts[iD2-1] - 1.0);
            else                   dtemp2 = 1.0;
            if ((sortedModifiedMeans[ii]-modifiedStds_[iD]/dtemp) >
                (sortedModifiedMeans[ii-1]+modifiedStds_[iD2]/dtemp2))
               printOutTS(2, "=============> Input %3d is different from input %3d.\n",
                  iD+1, iD2+1);
         }
      }
      printAsterisks(0);
      sigma2LowerBounds_ = new double[nInputs_];
      sigma2UpperBounds_ = new double[nInputs_];

      for (ii = 0; ii < 11; ii++) printOutTS(2, "-");
      printOutTS(2, " MOAT Analysis (ordered) : +- 2 sigma  ");
      for (ii = 0; ii < 11; ii++) printOutTS(2, "-");
      printOutTS(2, "\n");
      printDashes(0);
      for (ii = nInputs_-1; ii >= 0; ii--)
      {
         iD = indexesSortedByModifiedMeans_[ii];
         if (counts[iD] > 1) dtemp = sqrt((double) counts[iD] - 1.0);
         else                dtemp = 1.0;
         sigma2LowerBounds_[ii] = sortedModifiedMeans[ii]-2.0*modifiedStds_[iD]/dtemp;
         sigma2UpperBounds_[ii] = sortedModifiedMeans[ii]+2.0*modifiedStds_[iD]/dtemp;
         printOutTS(2, "Input %3d bounds = %12.4e %12.4e\n",
                iD+1,  sigma2LowerBounds_[ii], sigma2UpperBounds_[ii]);

         if (printLevel > 2 && ii > 0)
         {
            iD2 = indexesSortedByModifiedMeans_[ii-1];
            if (counts[iD2-1] > 1) dtemp2 = sqrt((double) counts[iD2-1] - 1.0);
            else                   dtemp2 = 1.0;
            if ((sortedModifiedMeans[ii]-2.0*modifiedStds_[iD]/dtemp) >
                (sortedModifiedMeans[ii-1]+2.0*modifiedStds_[iD2]/dtemp2))
               printOutTS(2, "=============> Input %3d is different from input %3d.\n",
                  iD+1, iD2+1);
         }
      }
   }

   if (psAnaExpertMode_ == 1)
   {
      printAsterisks(0);
      sprintf(pString, "Perform MOAT interaction study ? (y or n) ");
      getString(pString, winput);
      if (winput[0] == 'y')
      {
         printAsterisks(0);
         for (ii = 0; ii < 17; ii++) printOutTS(2, "*");
         printOutTS(2, " MOAT Interaction Analysis ");
         for (ii = 0; ii < 17; ii++) printOutTS(2, "*");
         printOutTS(2, "\n");
         printDashes(0);
         printOutTS(1, "** No data for an input means no interaction info.\n");
         printOutTS(1, "** Small stdev means (crudely) little interaction.\n");
         nPaths = nSamples_ / (nInputs_ + 1);
         XSort = new double[nPaths];
         YSort = new double[nPaths];
         for (ii = 0; ii < nInputs_; ii++)
         {
            actualPaths = 0;
            for (iD = 0; iD < nSamples_; iD++)
            {
               if (indexTrack[iD] == ii && YG[iD] != PSUADE_UNDEFINED)
               { 
                  XSort[actualPaths] = XG[iD]; 
                  YSort[actualPaths] = PABS(YG[iD]); 
                  actualPaths++;
               }
            }
            sortDbleList2(actualPaths, XSort, YSort);
            ip = 0; 
            iaCnt = 0;
            dstdev = 0.0;
            while (ip < actualPaths)
            {
               n1 = ip;
               for (ip = n1+1; ip < actualPaths; ip++)
                  if (XSort[ip] != XSort[ip-1]) break;
               n2 = ip;
               if ((n2 - n1) >= 2)
               {
                  dsum = 0.0;
                  for (iD = n1; iD < n2; iD++) dsum += YSort[iD];
                  dsum /= (double) (n2-n1);
                  dtemp = 0.0;
                  for (iD = n1; iD < n2; iD++) 
                     dtemp += (YSort[iD] - dsum) * (YSort[iD] - dsum);
                  dtemp /= (double) (n2-n1-1);
                  dtemp = sqrt(dtemp);
                  dstdev += dtemp;
                  iaCnt++;
                  if (printLevel > 3)
                     printOutTS(2, "Input %3d stdev at %12.4e = %12.4e (%d) \n",
                            ii+1, XSort[n1], dtemp, n2-n1);
               }
               else
                  if (printLevel > 3)
                     printOutTS(2, "Input %3d stdev at %12.4e = not available (%d).\n",
                            ii+1, XSort[n1], n2-n1);
            }
            if (iaCnt > 0)
               printOutTS(2, "Input %3d average interaction measure (std dev) = %11.3e\n",
                      ii+1, dstdev/iaCnt);
         }
         delete [] XSort;
         delete [] YSort;
      }
   }

   if (psAnaExpertMode_ == 1)
   {
      printAsterisks(0);
      sprintf(pString, "Perform hypothesis tests ? (y or n) ");
      getString(pString, winput);
      if (winput[0] == 'y')
      {
         printAsterisks(0);
         for (ii = 0; ii < 19; ii++) printOutTS(2, "*");
         printOutTS(2, "MOAT Hypothesis Tests");
         for (ii = 0; ii < 19; ii++) printOutTS(2, "*");
         printOutTS(2, "\n");
         printDashes(0);
         printOutTS(2, "* This consists of 2 tests: \n");
         printOutTS(2, "* (1) bootstrap confidence interval test (95 %%) \n");
         printOutTS(2, "* (2) binomial test (90 %%) \n");
         printOutTS(2, "* If either test indicates significance, the corresponding\n");
         printOutTS(2, "*   input will be declared significant.\n");
         printDashes(0);
         for (ii = 0; ii < nInputs_; ii++)
         {
            if (sortedModifiedMeans[ii] > binMax) binMax = sortedModifiedMeans[ii];
            if (sortedModifiedMeans[ii] < binMin) binMin = sortedModifiedMeans[ii];
         }
         sprintf(pString,"Enter the significance threshold (min,max = %e %e): ",
                 binMin, binMax);
         thresh = binMin - 1.0;
         while (thresh < binMin || thresh > binMax)
            thresh = getDouble(pString);
         YB = new double[nSamples_/(nInputs_+1)*100];
         bData.nOutputs_ = 1;
         bData.outputID_ = 0;
         bData.sampleOutputs_ = YB;
         bData.analysisThreshold_ = thresh;
         bData.printLevel_ = 0;
         for (ii = 0; ii < nInputs_; ii++)
         {
            printOutTS(2, ">>>>>> MOAT hypothesis tests for INPUT %3d : \n", ii+1);
            sigFlag = 0;
            actualPaths = 0;
            for (iD = 0; iD < nSamples_; iD++)
            {
               if (YG[iD] != PSUADE_UNDEFINED)
               {
                  index = indexTrack[iD];
                  if (index == ii) YB[actualPaths++] = PABS(YG[iD]);
               }
            }
            bData.nSamples_ = actualPaths;
            dtemp = bsAnalyzer.analyze(bData);
            if (dtemp > thresh) sigFlag = 1; 
            if (printLevel > 1)
            {
               printOutTS(2, "   ---- Bootstrap confidence interval = [0, %e] (>? %e)",
                      dtemp, thresh);
               if (sigFlag == 1) printOutTS(2, " **");
               printOutTS(2, "\n");
            } 
            jj = actualPaths;
            for (iD = 0; iD < actualPaths; iD++)
            {
               if (YB[iD] > 1.5*thresh)   
               {
                  dtemp = YB[iD];
                  while (dtemp > 1.5*thresh)
                  {
                     YB[jj++] = 1.5 * thresh;
                     dtemp -= (1.5 * thresh);
                     if (jj >= (nSamples_/nInputs_+1)*100)
                     {
                        printOutTS(2, "MOATAnalysis ERROR: binomial test.\n");
                        exit(1);
                     }
                  }
               }
            }
            actualPaths = jj;
            bData.nSamples_ = actualPaths;
            dtemp = binAnalyzer.analyze(bData);
            if (dtemp > 0.1) sigFlag++;
            if (printLevel > 1)
               printOutTS(2,"   ---- Binomial test ERROR = %4.1f %% if declared insignificant.\n",
                      100.0*dtemp);
            if (sigFlag > 0) 
                 printOutTS(2, "<<<<<< INPUT %3d IS SIGNIFICANT (%d) *****\n", 
                            ii+1,sigFlag); 
            else printOutTS(2, "<<<<<< INPUT %3d IS NOT SIGNIFICANT\n", ii+1); 
            
         }
         bData.sampleOutputs_ = NULL;
      }
   }

   delete [] counts;
   delete constrPtr;
   delete [] YY;
   delete [] YG;
   delete [] XG;
   delete [] indexTrack;
   delete [] Xbase;
   return 0.0;
}

// ************************************************************************
// functions for getting results 
// ------------------------------------------------------------------------
double MOATAnalyzer::get_analysisError()
{
   return analysisError_;
}

double* MOATAnalyzer::get_means() 
{
   double* retVal = NULL;
   if (means_) 
   {
      retVal = new double[nInputs_];
      std::copy(means_, means_+nInputs_, retVal);
   } 
   return retVal;
}

double* MOATAnalyzer::get_stds()
{
   double* retVal = NULL;
   if (stds_) 
   {
     retVal = new double[nInputs_];
     std::copy(stds_, stds_+nInputs_, retVal);
   } 
   return retVal;
}

double* MOATAnalyzer::get_modifiedMeans() 
{
   double* retVal = NULL;
   if(modifiedMeans_)
   {
      retVal = new double[nInputs_];
      std::copy(modifiedMeans_, modifiedMeans_+nInputs_, retVal);
   } 
   return retVal;
}

int* MOATAnalyzer::get_indexesSortedByModifiedMeans()
{
   int* retVal = NULL;
   if(indexesSortedByModifiedMeans_)
   {
      retVal = new int[nInputs_];
      std::copy(indexesSortedByModifiedMeans_, 
                indexesSortedByModifiedMeans_+nInputs_, retVal);
   } 
   return retVal;
}

double* MOATAnalyzer::get_sigma1LowerBounds() 
{
   double* retVal = NULL;
   if (sigma1LowerBounds_)
   {
      retVal = new double[nInputs_];
      std::copy(sigma1LowerBounds_, sigma1LowerBounds_+nInputs_, retVal);
   } 
   return retVal;
}

double* MOATAnalyzer::get_sigma1UpperBounds() 
{
   double* retVal = NULL;
   if (sigma1UpperBounds_)
   {
      retVal = new double[nInputs_];
      std::copy(sigma1UpperBounds_, sigma1UpperBounds_+nInputs_, retVal);
   } 
   return retVal;
}

double* MOATAnalyzer::get_sigma2LowerBounds() 
{
   double* retVal = NULL;
   if(sigma2LowerBounds_)
   {
      retVal = new double[nInputs_];
      std::copy(sigma2LowerBounds_, sigma2LowerBounds_+nInputs_, retVal);
   } 
   return retVal;
}

double* MOATAnalyzer::get_sigma2UpperBounds()
{
   double* retVal = NULL;
   if(sigma2UpperBounds_)
   {
      retVal = new double[nInputs_];
      std::copy(sigma2UpperBounds_, sigma2UpperBounds_+nInputs_, retVal);
   } 
   return retVal;
}

// ************************************************************************
// create Morris diagram matlab/scilab file
// ------------------------------------------------------------------------
int MOATAnalyzer::createScreenDiagramFile(int nSamples, int nInputs, 
                        double *Y, int *indices, double *modifiedMeans,
                        double *stds, int outputID, char **iNames)
{
   int  iD, ii, index, cnt;
   char winput[500], moatFile[500], pString[500];
   FILE *fp;

   printAsterisks(0);
   sprintf(pString,"Create screening diagram? (y or n) ");
   getString(pString, winput);
   if (winput[0] != 'y') return 0;
   sprintf(pString, "matlab/scilab screening diagram file name (no extension): ");
   getString(pString, moatFile);
   cnt = strlen(moatFile);
   if (cnt > 500)
   {
      printOutTS(0, "ERROR: file name too long.\n");
      exit(1);
   }
   moatFile[cnt-1] = '.';
   if (psPlotTool_ == 1)
   {
      moatFile[cnt] = 's';
      moatFile[cnt+1] = 'c';
      moatFile[cnt+2] = 'i';
      moatFile[cnt+3] = '\0';
   }
   else
   {
      moatFile[cnt] = 'm';
      moatFile[cnt+1] = '\0';
   }
   fp = fopen(moatFile, "w");

   if (fp == NULL)
   {
      printOutTS(0, "MOATAnalyzer: cannot open MOAT plot file %s\n",moatFile);
      return 0;
   }

   if (psPlotTool_ == 1)
   {
      fprintf(fp, "// Morris one-at-a-time screening plots\n");
      fprintf(fp, "// Z contains the gradient data\n");
      fprintf(fp, "// C contains the number of data for each input\n");
      fprintf(fp, "// turn plotAll on or off\n");
   }
   else
   {
      fprintf(fp, "%% Morris one-at-a-time screening plots\n");
      fprintf(fp, "%% Z contains the gradient data\n");
      fprintf(fp, "%% C contains the number of data for each input\n");
      fprintf(fp, "%% turn plotAll on or off\n");
   }
   fprintf(fp, "maxFlag = 0;\n");
   fprintf(fp, "plotAll = 1;\n");
   fprintf(fp, "Z = zeros(%d,%d);\n",nInputs,nSamples/(nInputs+1));
   fprintf(fp, "C = zeros(%d,1);\n",nInputs);
   fwritePlotCLF(fp);
   for (iD = 0; iD < nSamples; iD++)
   {
      if (Y[iD] != PSUADE_UNDEFINED)
      {
         index = indices[iD];
         if (index >= 0)
         {
            fprintf(fp, "C(%d) = C(%d) + 1;\n", index+1, index+1);
            fprintf(fp, "Z(%d,C(%d)) = %24.16e;\n",index+1,index+1,Y[iD]);
         }
      }
   }
   if (iNames == NULL)
   {
      fprintf(fp, "Str = {");
      for (ii = 0; ii < nInputs-1; ii++) fprintf(fp,"'X%d',",ii+1);
      fprintf(fp,"'X%d'};\n",nInputs);
   }
   else
   {
      fprintf(fp, "Str = {");
      for (ii = 0; ii < nInputs-1; ii++) 
      {
         if (iNames[ii] != NULL) fprintf(fp,"'%s',",iNames[ii]);
         else                    fprintf(fp,"'X%d',",ii+1);
      }
      if (iNames[nInputs-1] != NULL) fprintf(fp,"'%s'};\n",iNames[nInputs-1]);
      else                           fprintf(fp,"'X%d'};\n",nInputs);
   }
   if (psPlotTool_ == 1)
   {
      fprintf(fp, "// compute max and min for axis \n");
      fprintf(fp, "// XM : for computing mean \n");
      fprintf(fp, "// Xm : for computing max \n");
      fprintf(fp, "// XX : for computing modified mean \n");
      fprintf(fp, "// VV : counts for each input \n");
      fprintf(fp, "// YY : standard deviation for each input \n");
   }
   else
   {
      fprintf(fp, "%% compute max and min for axis \n");
      fprintf(fp, "%% XM : for computing mean \n");
      fprintf(fp, "%% Xm : for computing max \n");
      fprintf(fp, "%% XX : for computing modified mean \n");
      fprintf(fp, "%% VV : counts for each input \n");
      fprintf(fp, "%% YY : standard deviation for each input \n");
   }
   fprintf(fp, "nn = %d;\n",nInputs);
   fprintf(fp, "XX = zeros(nn,1);\n");
   fprintf(fp, "XM = zeros(nn,1);\n");
   fprintf(fp, "Xm = zeros(nn,1);\n");
   fprintf(fp, "YY = zeros(nn,1);\n");
   fprintf(fp, "for jj = 1 : nn\n");
   fprintf(fp, "   VV = zeros(nn,1);\n");
   fprintf(fp, "   for kk = 1 : C(jj)\n");
   fprintf(fp, "      if (VV(jj) <= C(jj))\n");
   fprintf(fp, "         if (abs(Z(jj,kk)) > Xm(jj))\n");
   fprintf(fp, "            Xm(jj) = abs(Z(jj,kk));\n");
   fprintf(fp, "         end;\n");
   fprintf(fp, "         XM(jj) = XM(jj) + Z(jj,kk);\n");
   fprintf(fp, "         XX(jj) = XX(jj) + abs(Z(jj,kk));\n");
   fprintf(fp, "         VV(jj) = VV(jj) + 1;\n");
   fprintf(fp, "      end;\n");
   fprintf(fp, "   end;\n");
   fprintf(fp, "   if (VV(jj) > 0)\n");
   fprintf(fp, "      XM(jj) = XM(jj) / VV(jj);\n");
   fprintf(fp, "      XX(jj) = XX(jj) / VV(jj);\n");
   fprintf(fp, "   end;\n");
   fprintf(fp, "   VV = zeros(nn,1);\n");
   fprintf(fp, "   for kk = 1 : C(jj)\n");
   fprintf(fp, "      if (VV(jj) <= C(jj))\n");
   fprintf(fp, "         YY(jj) = YY(jj) + (Z(jj,kk)-XM(jj))^2;\n");
   fprintf(fp, "         VV(jj) = VV(jj) + 1;\n");
   fprintf(fp, "      end;\n");
   fprintf(fp, "   end;\n");
   fprintf(fp, "   if (VV(jj) > 1)\n");
   fprintf(fp, "      YY(jj) = sqrt(YY(jj) / (VV(jj)-1));\n");
   fprintf(fp, "   end;\n");
   fprintf(fp, "end;\n");
   if (psPlotTool_ == 1)
        fprintf(fp, "// plot sequence of Morris plot\n");
   else fprintf(fp, "%% plot sequence of Morris plot\n");
   fprintf(fp, "last  = max(C);\n");
   fprintf(fp, "inc = floor((last - 2)/3);\n");
   fprintf(fp, "if (inc < 1)\n");
   fprintf(fp, "   inc = 1;\n");
   fprintf(fp, "end;\n");
   fprintf(fp, "list = [last : -inc : 2];\n");
   fprintf(fp, "if (length(list) > 4)\n");
   fprintf(fp, "   list = list(1:4);\n");
   fprintf(fp, "end;\n");
   if (psPlotTool_ == 1)
        fprintf(fp, "list = gsort(list,'g','i');\n");
   else fprintf(fp, "list = sort(list);\n");
   fprintf(fp, "list = unique(list);\n");
   fprintf(fp, "count = length(list);\n");
   fprintf(fp, "if (plotAll == 1)\n");
   if (psPlotTool_ == 1)
        fprintf(fp, "scf(1)\n");
   else fprintf(fp, "figure(1)\n");
   fprintf(fp, "for mm = 1 : count\n");
   fprintf(fp, "   ii = list(mm);\n");
   fprintf(fp, "   XX = zeros(nn,1);\n");
   fprintf(fp, "   XM = zeros(nn,1);\n");
   fprintf(fp, "   YY = zeros(nn,1);\n");
   fprintf(fp, "   VV = zeros(nn,1);\n");
   fprintf(fp, "   for jj = 1 : nn\n");
   fprintf(fp, "      cnt = 0;\n");
   fprintf(fp, "      for kk = 1 : ii\n");
   fprintf(fp, "         if (cnt <= C(jj))\n");
   fprintf(fp, "            XM(jj) = XM(jj) + Z(jj,kk);\n");
   fprintf(fp, "            XX(jj) = XX(jj) + abs(Z(jj,kk));\n");
   fprintf(fp, "            cnt = cnt + 1;\n");
   fprintf(fp, "         end;\n");
   fprintf(fp, "      end;\n");
   fprintf(fp, "      if (cnt > 0)\n");
   fprintf(fp, "         XM(jj) = XM(jj) / cnt;\n");
   fprintf(fp, "         XX(jj) = XX(jj) / cnt;\n");
   fprintf(fp, "      end;\n");
   fprintf(fp, "      cnt = 0;\n");
   fprintf(fp, "      for kk = 1 : ii\n");
   fprintf(fp, "         if (cnt <= C(jj))\n");
   fprintf(fp, "            YY(jj) = YY(jj) + (Z(jj,kk)-XM(jj))^2;\n");
   fprintf(fp, "            cnt = cnt + 1;\n");
   fprintf(fp, "         end;\n");
   fprintf(fp, "      end;\n");
   fprintf(fp, "      if (cnt > 1)\n");
   fprintf(fp, "         YY(jj) = sqrt(YY(jj) / (cnt-1));\n");
   fprintf(fp, "      end;\n");
   fprintf(fp, "   end;\n");
   fprintf(fp, "   subplot(2,2,mm) \n");
   fprintf(fp, "   bar(XX,0.8)\n");
   fwritePlotAxes(fp);
   fprintf(fp,"   pStr=sprintf('%%d replications',list(mm));\n");
   if (psPlotTool_ == 1)
   {
      fprintf(fp, "   a = gca();\n");
      fprintf(fp, "   a.title.text = pStr;\n");
      fprintf(fp, "   a.title.font_size = 3;\n");
      fprintf(fp, "   a.title.font_style = 4;\n");
   }
   else fprintf(fp,"   title(pStr);\n");
   fprintf(fp,"Xmin = min(XX) - (max(XX) - min(XX)) * 0.1;\n");
   fprintf(fp,"Xmax = max(XX) + (max(XX) - min(XX)) * 0.1;\n");
   if (psPlotTool_ == 1)
   {
      fprintf(fp, "   a=gca();\n");
      fprintf(fp, "   a.data_bounds=[0, Xmin; nn+1, Xmax];\n");
      fprintf(fp, "   a.x_ticks(2) = [1:nn]';\n");
      fprintf(fp, "   a.x_ticks(3) = Str';\n");
      fprintf(fp, "   a.x_label.font_size = 3;\n");
      fprintf(fp, "   a.x_label.font_style = 4;\n");
   }
   else
   {
      fprintf(fp,"    axis([0 nn+1 Xmin Xmax])\n");
      fprintf(fp, "   set(gca,'XTickLabel',[]);\n");
      fprintf(fp, "   th=text(1:nn, repmat(Xmin-0.1*(Xmax-Xmin),nn,1),Str,");
      fprintf(fp, "'HorizontalAlignment','left','rotation',90);\n");
      fprintf(fp, "   set(th, 'fontsize', 12)\n");
      fprintf(fp, "   set(th, 'fontweight', 'bold')\n");
   }
   fwritePlotYLabel(fp, "Modified Means (of gradients)");
   fprintf(fp, "end;\n");
   if (psPlotTool_ == 1)
        fprintf(fp, "scf(2)\n");
   else fprintf(fp, "figure(2)\n");
   fprintf(fp, "for mm = 1 : count\n");
   fprintf(fp, "   ii = list(mm);\n");
   fprintf(fp, "   XX = zeros(nn,1);\n");
   fprintf(fp, "   XM = zeros(nn,1);\n");
   fprintf(fp, "   YY = zeros(nn,1);\n");
   fprintf(fp, "   VV = zeros(nn,1);\n");
   fprintf(fp, "   for jj = 1 : nn\n");
   fprintf(fp, "      cnt = 0;\n");
   fprintf(fp, "      for kk = 1 : ii\n");
   fprintf(fp, "         if (cnt <= C(jj))\n");
   fprintf(fp, "            XM(jj) = XM(jj) + Z(jj,kk);\n");
   fprintf(fp, "            XX(jj) = XX(jj) + abs(Z(jj,kk));\n");
   fprintf(fp, "            cnt = cnt + 1;\n");
   fprintf(fp, "         end;\n");
   fprintf(fp, "      end;\n");
   fprintf(fp, "      if (cnt > 0)\n");
   fprintf(fp, "         XM(jj) = XM(jj) / cnt;\n");
   fprintf(fp, "         XX(jj) = XX(jj) / cnt;\n");
   fprintf(fp, "      end;\n");
   fprintf(fp, "      cnt = 0;\n");
   fprintf(fp, "      for kk = 1 : ii\n");
   fprintf(fp, "         if (cnt <= C(jj))\n");
   fprintf(fp, "            YY(jj) = YY(jj) + (Z(jj,kk)-XM(jj))^2;\n");
   fprintf(fp, "            cnt = cnt + 1;\n");
   fprintf(fp, "         end;\n");
   fprintf(fp, "      end;\n");
   fprintf(fp, "      if (cnt > 1)\n");
   fprintf(fp, "         YY(jj) = sqrt(YY(jj) / (cnt-1));\n");
   fprintf(fp, "      end;\n");
   fprintf(fp, "   end;\n");
   fprintf(fp, "   subplot(2,2,mm) \n");
   fprintf(fp, "   bar(YY,0.8)\n");
   fprintf(fp,"    Ymin = 0.0;\n");
   fprintf(fp,"    Ymax = max(YY) + (max(YY) - min(YY)) * 0.1;\n");
   fwritePlotAxes(fp);
   if (psPlotTool_ == 1)
   {
      fprintf(fp, "a=gca();\n");
      fprintf(fp, "a.data_bounds=[0, Ymin; nn+1, Ymax];\n");
      fprintf(fp, "   a.x_ticks(2) = [1:nn]';\n");
      fprintf(fp, "   a.x_ticks(3) = Str';\n");
      fprintf(fp, "   a.x_label.font_size = 3;\n");
      fprintf(fp, "   a.x_label.font_style = 4;\n");
   }
   else
   {
      fprintf(fp,"    axis([0 nn+1 Ymin Ymax])\n");
      fprintf(fp, "   set(gca,'XTickLabel',[]);\n");
      fprintf(fp, "   th=text(1:nn, repmat(Ymin-0.1*(Ymax-Ymin),nn,1),Str,");
      fprintf(fp, "'HorizontalAlignment','left','rotation',90);\n");
      fprintf(fp, "   set(th, 'fontsize', 12)\n");
      fprintf(fp, "   set(th, 'fontweight', 'bold')\n");
   }
   fwritePlotYLabel(fp, "Std. Devs.");
   fprintf(fp,"    pStr=sprintf('%%d replications',list(mm));\n");
   if (psPlotTool_ == 1)
   {
      fprintf(fp, "a = gca();\n");
      fprintf(fp, "a.title.text = pStr;\n");
      fprintf(fp, "a.title.font_size = 3;\n");
      fprintf(fp, "a.title.font_style = 4;\n");
   }
   else fprintf(fp,"    title(pStr);\n");
   fprintf(fp, "end;\n");
   fprintf(fp, "end;\n");
   fprintf(fp, "if (plotAll == 1)\n");
   if (psPlotTool_ == 1)
        fprintf(fp, "scf(3)\n");
   else fprintf(fp, "figure(3)\n");
   fprintf(fp, "end;\n");
   fprintf(fp, "Y = [\n");
   for (ii = 0; ii < nInputs; ii++) fprintf(fp, "%24.16e\n", stds[ii]);
   fprintf(fp, "];\n");
   fprintf(fp, "X = [\n");
   for (ii = 0; ii < nInputs; ii++)
      fprintf(fp, "%24.16e\n", modifiedMeans[ii]);
   fprintf(fp, "];\n");
   fprintf(fp, "plot(X,Y,'*','MarkerSize',12)\n");
   if (psPlotTool_ == 1)
   {
      fprintf(fp, "a=gca();\n");
      printf("Morris scilab plot: labels to be put in later.\n");
   }
   else
   {
      fprintf(fp, "text(X*1.01,Y,{");
      if (iNames != NULL)
      {
         for (ii = 0; ii < nInputs-1; ii++) fprintf(fp, "'%s',",iNames[ii]);
         fprintf(fp, "'%s'},'FontWeight','bold','FontSize',12)\n",
                 iNames[nInputs-1]);
      }
      else
      {
         for (ii = 0; ii < nInputs-1; ii++) fprintf(fp, "'X%d',",ii+1);
         fprintf(fp, "'X%d'},'FontWeight','bold','FontSize',12)\n",nInputs);
      }
   }
   fprintf(fp, "Xmin = min(X) - (max(X) - min(X)) * 0.1;\n");
   fprintf(fp, "Ymin = min(Y) - (max(Y) - min(Y)) * 0.1;\n");
   fprintf(fp, "Xmax = max(X) + (max(X) - min(X)) * 0.1;\n");
   fprintf(fp, "Ymax = max(Y) + (max(Y) - min(Y)) * 0.1;\n");
   if (psPlotTool_ == 1)
   {
      fprintf(fp, "a=gca();\n");
      fprintf(fp, "a.data_bounds=[Xmin, Ymin; Xmax, Ymax];\n");
   }
   else fprintf(fp, "axis([Xmin Xmax Ymin Ymax])\n");
   fwritePlotAxes(fp);
   fwritePlotXLabel(fp, "Modified Means (of gradients)");
   fwritePlotYLabel(fp, "Std Deviations (of gradients)");
   sprintf(pString, "Modified Morris Diagram for Output %d", outputID+1);
   fwritePlotTitle(fp, pString);
   fprintf(fp, "if (maxFlag == 1)\n");
   if (psPlotTool_ == 1)
        fprintf(fp, "   scf(4)\n");
   else fprintf(fp, "   figure(4)\n");
   fprintf(fp, "   bar(Xm, 0.8)\n");
   fprintf(fp,"    Ymin = 0.0;\n");
   fprintf(fp,"    Ymax = max(YY) + (max(YY) - min(YY)) * 0.1;\n");
   fwritePlotAxes(fp);
   if (psPlotTool_ == 1)
   {
      fprintf(fp, "   a=gca();\n");
      fprintf(fp, "   a.data_bounds=[0, Ymin; nn+1, Ymax];\n");
      fprintf(fp, "   a.x_ticks(2) = [1:nn]';\n");
      fprintf(fp, "   a.x_ticks(3) = Str';\n");
      fprintf(fp, "   a.x_label.font_size = 3;\n");
      fprintf(fp, "   a.x_label.font_style = 4;\n");
   }
   else
   {
      fprintf(fp, "   axis([0 nn+1 Ymin Ymax])\n");
      fprintf(fp, "   set(gca,'XTickLabel',[]);\n");
      fprintf(fp, "   th=text(1:nn, repmat(Ymin-0.1*(Ymax-Ymin),nn,1),Str,");
      fprintf(fp, "'HorizontalAlignment','left','rotation',90);\n");
      fprintf(fp, "   set(th, 'fontsize', 12)\n");
      fprintf(fp, "   set(th, 'fontweight', 'bold')\n");
   }
   fwritePlotYLabel(fp, "   Max Absolute Gradients");
   fprintf(fp, "end;\n");
   printOutTS(2, "MOAT screening diagram file = %s\n", moatFile);
   fclose(fp); 
   printAsterisks(0);
   return 0;
}

// ************************************************************************
// create scatter plot matlab file
// ------------------------------------------------------------------------
int MOATAnalyzer::createScatterFile(int nSamples, int nInputs, double *Y, 
                                    double *X, int *indices, char **iNames)
{
   int  iD, cnt, ii;
   char winput[500], scatterFile[500], pString[500];
   FILE *fp;

   printAsterisks(0);
   sprintf(pString, "Create scatter plot ? (y or n) ");
   getString(pString, winput);
   if (winput[0] != 'y')
   {
      return 0;
   }
   sprintf(pString,"Enter matlab/scilab scatter plot file name (no extension): ");
   getString(pString, scatterFile);
   cnt = strlen(scatterFile);
   if (cnt > 500)
   {
     printOutTS(0, "ERROR: file name too long.\n");
      exit(1);
   }
   scatterFile[cnt-1] = '.';
   if (psPlotTool_ == 1)
   {
      scatterFile[cnt] = 's';
      scatterFile[cnt+1] = 'c';
      scatterFile[cnt+2] = 'i';
      scatterFile[cnt+3] = '\0';
   }
   else
   {
      scatterFile[cnt] = 'm';
      scatterFile[cnt+1] = '\0';
   }

   fp = fopen(scatterFile, "w");
   if (fp == NULL)
   {
      printOutTS(0, "MOATAnalysis: cannot open scatterplot file %s\n",scatterFile);
      return 0;
   }

   if (psPlotTool_ == 1)
   {
      fprintf(fp, "// This file contains individual gradient info.\n");
      fprintf(fp, "// The gradients are normalized for the range.\n");
      fprintf(fp, "// To display only the important inputs, first\n");
      fprintf(fp, "// create and run the bootstrapped analysis\n");
      fprintf(fp, "// which creates an I2 array. \n");
   }
   else
   {
      fprintf(fp, "%% This file contains individual gradient info.\n");
      fprintf(fp, "%% The gradients are normalized for the range.\n");
      fprintf(fp, "%% To display only the important inputs, first\n");
      fprintf(fp, "%% create and run the bootstrapped analysis\n");
      fprintf(fp, "%% which creates an I2 array. \n");
   }
   fprintf(fp, "sortFlag = 0;\n");
   fprintf(fp, "if (sortFlag == 1);\n");
   fprintf(fp, "  nn = length(I2);\n");
   fprintf(fp, "else\n");
   fprintf(fp, "  nn = %d;\n", nInputs);
   fprintf(fp, "end\n");
   fprintf(fp, "A = [\n");
   for (iD = 1; iD < nSamples; iD++)
   {
      if (Y[iD] != PSUADE_UNDEFINED)
         fprintf(fp,"%4d %24.16e %24.16e %d\n",indices[iD]+1,Y[iD],X[iD],iD+1);
   }
   fprintf(fp, "];\n");
   if (iNames == NULL)
   {
      fprintf(fp, "Str = {");
      for (ii = 0; ii < nInputs-1; ii++) fprintf(fp,"'X%d',",ii+1);
      fprintf(fp,"'X%d'};\n",nInputs);
   }
   else
   {
      fprintf(fp, "Str = {");
      for (ii = 0; ii < nInputs-1; ii++) 
      {
         if (iNames[ii] != NULL) fprintf(fp,"'%s',",iNames[ii]);
         else                    fprintf(fp,"'X%d',",ii+1);
      }
      if (iNames[nInputs-1] != NULL) fprintf(fp,"'%s'};\n",iNames[nInputs-1]);
      else                           fprintf(fp,"'X%d'};\n",nInputs);
   }
   fprintf(fp, "if (sortFlag == 1);\n");
   fprintf(fp, "  Str = Str(I2);\n");
   fprintf(fp, "end\n");

   fwritePlotCLF(fp);
   if (psPlotTool_ == 1)
        fprintf(fp, "set(gca(),\"auto_clear\",\"off\")\n");
   else fprintf(fp, "hold on\n");
   fprintf(fp, "ymin = min(A(:,2));\n");
   fprintf(fp, "ymax = max(A(:,2));\n");
   if (psPlotTool_ == 1)
   {
      fprintf(fp, "a=gca();\n");
      fprintf(fp, "a.data_bounds=[0, ymin; nn+1, ymax];\n");
      fprintf(fp, "newtick = a.x_ticks;\n");
      fprintf(fp, "newtick(2) = [1:nn]';\n");
      fprintf(fp, "newtick(3) = Str';\n");
      fprintf(fp, "a.x_ticks = newtick;\n");
      fprintf(fp, "a.x_label.font_size = 3;\n");
      fprintf(fp, "a.x_label.font_style = 4;\n");
   }
   else
   {
      fprintf(fp, "axis([0  nn+1 ymin ymax])\n");
      fprintf(fp, "set(gca,'XTickLabel',[]);\n");
      fprintf(fp, "th=text(1:nn, repmat(ymin-0.05*(ymax-ymin),nn,1),Str,");
      fprintf(fp, "'HorizontalAlignment','left','rotation',90);\n");
      fprintf(fp, "set(th, 'fontsize', 12)\n");
      fprintf(fp, "set(th, 'fontweight', 'bold')\n");
   }
   fprintf(fp, "for ii2 = 1 : nn\n");
   fprintf(fp, "  if (sortFlag == 1)\n");
   fprintf(fp, "     ii = I2(ii2);\n");
   fprintf(fp, "  else\n");
   fprintf(fp, "     ii = ii2;\n");
   fprintf(fp, "  end;\n");
   fprintf(fp, "  inds = find(A(:,1)==ii);\n");
   fprintf(fp, "  leng = length(inds);\n");
   fprintf(fp, "  if (leng > 0)\n");
   if (psPlotTool_ == 1)
        fprintf(fp, "    [AA,II] = gsort(A(inds,3),'g','i');\n");
   else fprintf(fp, "    [AA,II] = sort(A(inds,3));\n");
   fprintf(fp, "    if (sortFlag == 1)\n");
   fprintf(fp, "       xx = ii2;\n");
   fprintf(fp, "    else\n");
   fprintf(fp, "       xx = A(inds(II(1)),1);\n");
   fprintf(fp, "    end;\n");
   fprintf(fp, "    yy = A(inds(II(1)),2);\n");
   fprintf(fp, "    plot(xx,yy,'ro','MarkerSize',11,'MarkerFaceColor','r',");
   fprintf(fp, "'MarkerEdgeColor','k')\n");
   fprintf(fp, "    for jj = 2 : leng\n");
   fprintf(fp, "      x1 = A(inds(II(jj)),3);\n");
   fprintf(fp, "      x2 = A(inds(II(jj-1)),3);\n");
   fprintf(fp, "      if x1 ~= x2\n");
   fprintf(fp, "        break;\n");
   fprintf(fp, "      else\n");
   fprintf(fp, "        if (sortFlag == 1)\n");
   fprintf(fp, "           xx = ii2;\n");
   fprintf(fp, "        else\n");
   fprintf(fp, "           xx = A(inds(II(jj)),1);\n");
   fprintf(fp, "        end;\n");
   fprintf(fp, "        yy = A(inds(II(jj)),2);\n");
   fprintf(fp, "        plot(xx,yy,'ro','MarkerSize',11,'MarkerFaceColor','r',");
   fprintf(fp, "'MarkerEdgeColor','k')\n");
   fprintf(fp, "        if (jj == leng)\n");
   fprintf(fp, "          jj = jj + 1;\n");
   fprintf(fp, "        end;\n");
   fprintf(fp, "      end;\n");
   fprintf(fp, "    end;\n");
   fprintf(fp, "    if (jj <= leng)\n");
   fprintf(fp, "      if (sortFlag == 1)\n");
   fprintf(fp, "         xx = ii2;\n");
   fprintf(fp, "      else\n");
   fprintf(fp, "         xx = A(inds(II(jj)),1);\n");
   fprintf(fp, "      end;\n");
   fprintf(fp, "      yy = A(inds(II(jj)),2);\n");
   fprintf(fp, "      plot(xx+0.1,yy,'go','MarkerSize',11,'MarkerFaceColor','g',");
   fprintf(fp, "'MarkerEdgeColor','k')\n");
   fprintf(fp, "    end;\n");
   fprintf(fp, "    for kk = jj+1 : leng\n");
   fprintf(fp, "      x1 = A(inds(II(kk)),3);\n");
   fprintf(fp, "      x2 = A(inds(II(kk-1)),3);\n");
   fprintf(fp, "      if x1 ~= x2\n");
   fprintf(fp, "        break;\n");
   fprintf(fp, "      else\n");
   fprintf(fp, "        if (sortFlag == 1)\n");
   fprintf(fp, "           xx = ii2;\n");
   fprintf(fp, "        else\n");
   fprintf(fp, "           xx = A(inds(II(kk)),1);\n");
   fprintf(fp, "        end;\n");
   fprintf(fp, "        yy = A(inds(II(kk)),2);\n");
   fprintf(fp, "        plot(xx+0.1,yy,'go','MarkerSize',11,'MarkerFaceColor','g',");
   fprintf(fp, "'MarkerEdgeColor','k')\n");
   fprintf(fp, "        if (kk == leng)\n");
   fprintf(fp, "          kk = kk + 1;\n");
   fprintf(fp, "        end;\n");
   fprintf(fp, "      end;\n");
   fprintf(fp, "    end;\n");
   fprintf(fp, "    if (kk <= leng)\n");
   fprintf(fp, "      if (sortFlag == 1)\n");
   fprintf(fp, "         xx = ii2;\n");
   fprintf(fp, "      else\n");
   fprintf(fp, "         xx = A(inds(II(kk)),1);\n");
   fprintf(fp, "      end;\n");
   fprintf(fp, "      yy = A(inds(II(kk)),2);\n");
   fprintf(fp, "      plot(xx+0.2,yy,'bo','MarkerSize',11,'MarkerFaceColor','b',");
   fprintf(fp, "'MarkerEdgeColor','k')\n");
   fprintf(fp, "    end;\n");
   fprintf(fp, "    for jj = kk+1 : leng\n");
   fprintf(fp, "      x1 = A(inds(II(jj)),3);\n");
   fprintf(fp, "      x2 = A(inds(II(jj-1)),3);\n");
   fprintf(fp, "      if x1 ~= x2\n");
   fprintf(fp, "        break;\n");
   fprintf(fp, "      else\n");
   fprintf(fp, "        if (sortFlag == 1)\n");
   fprintf(fp, "           xx = ii2;\n");
   fprintf(fp, "        else\n");
   fprintf(fp, "           xx = A(inds(II(jj)),1);\n");
   fprintf(fp, "        end;\n");
   fprintf(fp, "        yy = A(inds(II(jj)),2);\n");
   fprintf(fp, "        plot(xx+0.2,yy,'bo','MarkerSize',11,'MarkerFaceColor','b',");
   fprintf(fp, "'MarkerEdgeColor','k')\n");
   fprintf(fp, "        if (jj == leng)\n");
   fprintf(fp, "          jj = jj + 1;\n");
   fprintf(fp, "        end;\n");
   fprintf(fp, "      end;\n");
   fprintf(fp, "    end;\n");
   fprintf(fp, "    if (jj <= leng)\n");
   fprintf(fp, "      xx = A(inds(II(jj)),1);\n");
   fprintf(fp, "      if (sortFlag == 1)\n");
   fprintf(fp, "         xx = ii2;\n");
   fprintf(fp, "      else\n");
   fprintf(fp, "         xx = A(inds(II(jj)),1);\n");
   fprintf(fp, "      end;\n");
   fprintf(fp, "      yy = A(inds(II(jj)),2);\n");
   fprintf(fp, "      plot(xx+0.3,yy,'mo','MarkerSize',11,'MarkerFaceColor','m',");
   fprintf(fp, "'MarkerEdgeColor','k')\n");
   fprintf(fp, "    end;\n");
   fprintf(fp, "    for kk = jj+1 : leng\n");
   fprintf(fp, "      x1 = A(inds(II(kk)),3);\n");
   fprintf(fp, "      x2 = A(inds(II(kk-1)),3);\n");
   fprintf(fp, "      if x1 ~= x2\n");
   fprintf(fp, "        break;\n");
   fprintf(fp, "      else\n");
   fprintf(fp, "        if (sortFlag == 1)\n");
   fprintf(fp, "           xx = ii2;\n");
   fprintf(fp, "        else\n");
   fprintf(fp, "           xx = A(inds(II(kk)),1);\n");
   fprintf(fp, "        end;\n");
   fprintf(fp, "        yy = A(inds(II(kk)),2);\n");
   fprintf(fp, "        plot(xx+0.3,yy,'mo','MarkerSize',11,'MarkerFaceColor','m',");
   fprintf(fp, "'MarkerEdgeColor','k')\n");
   fprintf(fp, "        if (kk == leng)\n");
   fprintf(fp, "          kk = kk + 1;\n");
   fprintf(fp, "        end;\n");
   fprintf(fp, "      end;\n");
   fprintf(fp, "    end;\n");
   fprintf(fp, "    if (kk <= leng)\n");
   fprintf(fp, "      for jj = kk : leng\n");
   fprintf(fp, "        if (sortFlag == 1)\n");
   fprintf(fp, "           xx = ii2;\n");
   fprintf(fp, "        else\n");
   fprintf(fp, "           xx = A(inds(II(jj)),1);\n");
   fprintf(fp, "        end;\n");
   fprintf(fp, "        yy = A(inds(II(jj)),2);\n");
   fprintf(fp, "        plot(xx+0.4,yy,'co','MarkerSize',11,'MarkerFaceColor','c',");
   fprintf(fp, "'MarkerEdgeColor','k')\n");
   fprintf(fp, "      end;\n");
   fprintf(fp, "    end;\n");
   fprintf(fp, "  end;\n");
   fprintf(fp, "end;\n");
   fprintf(fp, "for ii2 = 1 : nn\n");
   fprintf(fp, "   if (sortFlag == 1)\n");
   fprintf(fp, "      ii = I2(ii2);\n");
   fprintf(fp, "   else\n");
   fprintf(fp, "      ii = ii2;\n");
   fprintf(fp, "   end;\n");
   fprintf(fp, "   inds = find(A(:,1)==ii);\n");
   fprintf(fp, "   if length(inds) > 0\n");
   fprintf(fp, "      asum = sum(abs(A(inds,2)))/length(inds);\n");
   if (psPlotTool_ == 1)
        fprintf(fp, "      plot(ii2,asum,'kp','MarkerSize',14);\n");
   else fprintf(fp, "      plot(ii2,asum,'kh','MarkerSize',14);\n");
   fprintf(fp, "   end;\n");
   fprintf(fp, "end;\n");
   fwritePlotAxes(fp);
   fwritePlotYLabel(fp, "Individual Gradients");
   fwritePlotTitle(fp, "Scatter Plots for the Gradients (colored)");
   fprintf(fp, "disp('from lo to hi : red,green,blue,magenta,cyan')\n");
   fprintf(fp, "disp('hexagrams: means of the absolute gradients')\n");
   if (psPlotTool_ == 1)
        fprintf(fp, "set(gca(),\"auto_clear\",\"on\")\n");
   else fprintf(fp, "hold off\n");
   fclose(fp);
   printOutTS(2, "MOAT scatter plot file = %s.\n", scatterFile);
   printAsterisks(0);
   return 0;
}

// ************************************************************************
// create bootstrap mean plot matlab file
// ------------------------------------------------------------------------
int MOATAnalyzer::createBootstrapFile(int nSamples, int nInputs, double *Y, 
                                      double *X, int *indices, char **iNames)
{
   int    nReps, *validCnts, ii, is, jj, ib, input, count, index;
   double **YGs, *means, *stds, *bArray;
   char   winput[500], bootstrapFile[500], pString[500];
   FILE   *fp;

   // range check nInputs by Bill Oliver
   if(nInputs <= 0)
   {
      printOutTS(0, "nInputs is <= 0 in file %s, line %d. returning -1\n",
             __FILE__, __LINE__);
      return -1;
   }

   printAsterisks(0);
   sprintf(pString, "Create bootstrap modified mean plot ? (y or n) ");
   getString(pString, winput);
   if (winput[0] != 'y') return 0;
   sprintf(pString,"Enter matlab/scilab bootstrap file name (no extension): ");
   getString(pString, bootstrapFile);
   count = strlen(bootstrapFile);
   if (count > 500)
   {
      printOutTS(0, "ERROR: file name too long.\n");
      exit(1);
   }
   bootstrapFile[count-1] = '.';
   if (psPlotTool_ == 1)
   {
      bootstrapFile[count] = 's';
      bootstrapFile[count+1] = 'c';
      bootstrapFile[count+2] = 'i';
      bootstrapFile[count+3] = '\0';
   }
   else
   {
      bootstrapFile[count] = 'm';
      bootstrapFile[count+1] = '\0';
   }
   fp = fopen(bootstrapFile, "w");
   if (fp == NULL)
   {
      printOutTS(0, "MOATAnalysis: cannot write to plot file %s\n",
                 bootstrapFile);
      return 0;
   }

   nReps = nSamples / (nInputs + 1);
   validCnts = new int[nInputs];
   for (ii = 0; ii < nInputs; ii++) validCnts[ii] = 0;
   for (is = 0; is < nSamples; is++)
   {
      if (Y[is] != PSUADE_UNDEFINED && indices[is] >= 0)
         validCnts[indices[is]]++;
   }
   for (ii = 0; ii < nInputs; ii++)
      if (validCnts[ii] > nReps) nReps = validCnts[ii];
   YGs = new double*[nInputs]; 
   for (ii = 0; ii < nInputs; ii++)
   {
      if (validCnts[ii] > 0) YGs[ii] = new double[validCnts[ii]];
      else                   YGs[ii] = NULL;
   }
   for (ii = 0; ii < nInputs; ii++) validCnts[ii] = 0;
   for (is = 0; is < nSamples; is++)
   {
      if (Y[is] != PSUADE_UNDEFINED && indices[is] >= 0)
      {
         input = indices[is];
         count = validCnts[input]++;
         YGs[input][count] = PABS(Y[is]);
      }
   }

   bArray = new double[1000];
   means = new double[nInputs];
   stds = new double[nInputs];
   for (ii = 0; ii < nInputs; ii++)
   {
      if (validCnts[ii] > 5)
      {
         for (ib = 0; ib < 1000; ib++)
         {
            bArray[ib] = 0.0;
            for (jj = 0; jj < validCnts[ii]; jj++)
            {
               index = PSUADE_rand() % validCnts[ii];
               bArray[ib] += YGs[ii][index];
            }
            bArray[ib] /= validCnts[ii];
         }
         means[ii] = 0.0;
         for (ib = 0; ib < 1000; ib++) means[ii] += bArray[ib];
         means[ii] /= 1000.0;
         stds[ii] = 0.0;
         for (ib = 0; ib < 1000; ib++)
            stds[ii] += pow(bArray[ib] - means[ii], 2.0);
         stds[ii] /= (1000.0 - 1.0);
         stds[ii] = sqrt(stds[ii]);
      }
      else
      {
         printOutTS(1,"MOATAnalysis WARNING: input %d needs > 5 replications\n",
                    ii+1);
         printOutTS(1,"                      to perform bootstrapping.\n");
         stds[ii] = 0.0;
         means[ii] = 0.0;
         for (jj = 0; jj < validCnts[ii]; jj++) means[ii] += YGs[ii][jj];
         if (validCnts[ii] > 0) means[ii] /= (double) validCnts[ii];
      }
   }

   if (psPlotTool_ == 1)
   {
      fprintf(fp,"// This file contains modified means of gradients\n");
      fprintf(fp,"// and also their spreads based on bootstraping.\n");
      fprintf(fp,"// to select the most important ones to display,\n");
      fprintf(fp,"// set sortFlag = 1 and set nn to be the number\n");
      fprintf(fp,"// of inputs to display.\n");
   }
   else
   {
      fprintf(fp,"%% This file contains modified means of gradients\n");
      fprintf(fp,"%% and also their spreads based on bootstraping.\n");
      fprintf(fp,"%% to select the most important ones to display,\n");
      fprintf(fp,"%% set sortFlag = 1 and set nn to be the number\n");
      fprintf(fp,"%% of inputs to display.\n");
   }
   fprintf(fp, "sortFlag = 0;\n");
   fprintf(fp, "nn = %d;\n", nInputs);
   fprintf(fp, "Means = [\n");
   for (ii = 0; ii < nInputs; ii++) fprintf(fp,"%24.16e\n", means[ii]);
   fprintf(fp, "];\n");
   fprintf(fp, "Stds = [\n");
   for (ii = 0; ii < nInputs; ii++) fprintf(fp,"%24.16e\n", stds[ii]);
   fprintf(fp, "];\n");
   if (iNames == NULL)
   {
      fprintf(fp, "  Str = {");
      for (ii = 0; ii < nInputs-1; ii++) fprintf(fp,"'X%d',",ii+1);
      fprintf(fp,"'X%d'};\n",nInputs);
   }
   else
   {
      fprintf(fp, "  Str = {");
      for (ii = 0; ii < nInputs-1; ii++) 
      {
         if (iNames[ii] != NULL) fprintf(fp,"'%s',",iNames[ii]);
         else                    fprintf(fp,"'X%d',",ii+1);
      }
      if (iNames[nInputs-1] != NULL) fprintf(fp,"'%s'};\n",iNames[nInputs-1]);
      else                           fprintf(fp,"'X%d'};\n",nInputs);
   }
   fwritePlotCLF(fp);
   fprintf(fp, "if (sortFlag == 1)\n");
   if (psPlotTool_ == 1)
        fprintf(fp, "  [Means, I2] = gsort(Means,'g','d');\n");
   else fprintf(fp, "  [Means, I2] = sort(Means,'descend');\n");
   fprintf(fp, "  Stds = Stds(I2);\n");
   fprintf(fp, "  I2 = I2(1:nn);\n");
   fprintf(fp, "  Means = Means(1:nn);\n");
   fprintf(fp, "  Stds = Stds(1:nn);\n");
   fprintf(fp, "  Str  = Str(I2);\n");
   fprintf(fp, "end\n");
   fprintf(fp, "ymin = min(Means-Stds);\n");
   fprintf(fp, "ymax = max(Means+Stds);\n");
   fprintf(fp, "h2 = 0.05 * (ymax - ymin);\n");
   fprintf(fp, "bar(Means,0.8);\n");
   fprintf(fp, "for ii = 1:nn\n");
   fprintf(fp, "   if (ii == 1)\n");
   if (psPlotTool_ == 1)
        fprintf(fp, "   set(gca(),\"auto_clear\",\"off\")\n");
   else fprintf(fp, "   hold on\n");
   fprintf(fp, "   end;\n");
   fprintf(fp, "   XX = [ii ii];\n");
   fprintf(fp, "   YY = [Means(ii)-Stds(ii) Means(ii)+Stds(ii)];\n");
   fprintf(fp, "   plot(XX,YY,'-ko','LineWidth',3.0,'MarkerEdgeColor',");
   fprintf(fp, "'k','MarkerFaceColor','g','MarkerSize',12)\n");
   fprintf(fp, "end;\n");
   if (psPlotTool_ == 1)
   {
      fprintf(fp, "a=gca();\n");
      fprintf(fp, "a.data_bounds=[0, ymin; nn+1, ymax];\n");
      fprintf(fp, "a.x_ticks(2) = [1:nn]';\n");
      fprintf(fp, "a.x_ticks(3) = Str';\n");
      fprintf(fp, "a.x_label.font_size = 3;\n");
      fprintf(fp, "a.x_label.font_style = 4;\n");
   }
   else
   {
      fprintf(fp, "axis([0 nn+1 ymin ymax])\n");
      fprintf(fp, "set(gca,'XTickLabel',[]);\n");
      fprintf(fp, "th=text(1:nn, repmat(ymin-0.05*(ymax-ymin),nn,1),Str,");
      fprintf(fp, "'HorizontalAlignment','left','rotation',90);\n");
      fprintf(fp, "set(th, 'fontsize', 12)\n");
      fprintf(fp, "set(th, 'fontweight', 'bold')\n");
   }
   fwritePlotAxes(fp);
   fwritePlotTitle(fp,"Modified Means Plot (bootstrap)");
   fwritePlotYLabel(fp, "Modified Means (of gradients)");
   if (psPlotTool_ == 1)
        fprintf(fp, "   set(gca(),\"auto_clear\",\"on\")\n");
   else fprintf(fp, "   hold off\n");
   fclose(fp);
   printOutTS(2, "MOAT bootstrap plot file = %s.\n", bootstrapFile);
   printAsterisks(0);

   for (ii = 0; ii < nInputs; ii++) if (YGs[ii] != NULL) delete [] YGs[ii];
   delete [] YGs;
   delete [] validCnts;
   delete [] means;
   delete [] stds;
   delete [] bArray;
   return 0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
MOATAnalyzer& MOATAnalyzer::operator=(const MOATAnalyzer &)
{
   printOutTS(0, "MOATAnalysis operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

