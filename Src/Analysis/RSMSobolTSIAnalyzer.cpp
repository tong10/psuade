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
// Functions for the class RSMSobolTSIAnalyzer  
// AUTHOR : CHARLES TONG
// DATE   : 2006
// ************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Util/PsuadeUtil.h"
#include "Util/sysdef.h"
#include "Util/Matrix.h"
#include "Util/Vector.h"
#include "Main/Psuade.h"
#include "FuncApprox/FuncApprox.h"
#include "Samplings/Sampling.h"
#include "Samplings/RSConstraints.h"
#include "PDFLib/PDFManager.h"
#include "DataIO/PsuadeData.h"
#include "RSMSobolTSIAnalyzer.h"

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
RSMSobolTSIAnalyzer::RSMSobolTSIAnalyzer() : Analyzer()
{
   setName("RSMSOBOLTSI");
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
RSMSobolTSIAnalyzer::~RSMSobolTSIAnalyzer()
{
}

// ************************************************************************
// perform analysis
// ------------------------------------------------------------------------
double RSMSobolTSIAnalyzer::analyze(aData &adata)
{
   int    nInputs, nOutputs, nSamples, ii, jj, kk, *S, status, outputID;
   int    length, nSubSamples=20000, *sampleStates, iL, sCnt, nLevels=200;
   int    *pdfFlags, printLevel, *bins, totalCnt, nSamp, *pdfFlags1;
   int    corFlag;
   double *xLower, *xUpper, *X, *Y, *cLower, *cUpper, *sampleInputs;
   double *sampleOutputs, *oneSamplePt, *Y2, *means,*tsi, variance;
   double *vars, dmean, dvar, *tsiRes, ddata, *mSamplePts, *inputMeans;
   double *inputStdevs, *inputMeans1, *inputStdevs1, *samplePts1D, *Y3;
   char   pString[500];
   Sampling      *sampler;
   PsuadeData    *ioPtr;
   FuncApprox    *faPtr;
   RSConstraints *constrPtr;
   Vector        vecIn, vecOut, vecUB, vecLB;
   pData         pCorMat;
   Matrix        *corMatp, corMat, corMat1, corMat2;
   PDFManager    *pdfman, *pdfman1, *pdfman2;

   printLevel  = adata.printLevel_;
   nInputs     = adata.nInputs_;
   nOutputs    = adata.nOutputs_;
   nSamples    = adata.nSamples_;
   xLower      = adata.iLowerB_;
   xUpper      = adata.iUpperB_;
   X           = adata.sampleInputs_;
   Y3          = adata.sampleOutputs_;
   S           = adata.sampleStates_;
   outputID    = adata.outputID_;
   ioPtr       = adata.ioPtr_;
   pdfFlags    = adata.inputPDFs_;
   inputMeans  = adata.inputMeans_;
   inputStdevs = adata.inputStdevs_;
   if (pdfFlags != NULL)
   {
      printf("RSMSobolTSI INFO: non-uniform distributions detected which\n");
      printf("                  will be used in this analysis.\n");
   }

   if (nInputs <= 0 || nSamples <= 0 || nOutputs <= 0) 
   {
      printf("RSMSobolTSI ERROR: invalid arguments.\n");
      printf("   nInputs  = %d\n", nInputs);
      printf("   nOutputs = %d\n", nOutputs);
      printf("   nSamples = %d\n", nSamples);
      return PSUADE_UNDEFINED;
   } 
   if (nInputs <= 1)
   {
      printf("RSMSobolTSI: nInputs<=1 does not need this analysis.\n");
      return PSUADE_UNDEFINED;
   }
   if (outputID >= nOutputs || outputID < 0)
   {
      printf("RSMSobolTSI ERROR: invalid output ID (%d).\n", outputID);
      return PSUADE_UNDEFINED;
   }
   if (ioPtr == NULL)
   {
      printf("RSMSobolTSI ERROR: no data.\n");
      return PSUADE_UNDEFINED;
   }
   ioPtr->getParameter("input_cor_matrix", pCorMat);
   corMatp = (Matrix *) pCorMat.psObject_;
   for (ii = 0; ii < nInputs; ii++)
   {
      for (jj = 0; jj < ii; jj++)
      {
         if (corMatp->getEntry(ii,jj) != 0.0)
         {
            printf("RSMSobolTSI INFO: this method cannot handle correlated\n");
            printf("          inputs using joint PDFs yet. Use group\n");
            printf("          variance-based method.\n");
            return PSUADE_UNDEFINED;
         }
      }
   }
   Y = new double[nSamples];
   for (ii = 0; ii < nSamples; ii++) Y[ii] = Y3[ii*nOutputs+outputID];

   constrPtr = new RSConstraints();
   constrPtr->genConstraints(ioPtr);

   faPtr = genFAInteractive(ioPtr, 0);
   length = -999;
   status = faPtr->genNDGridData(X, Y, &length, NULL, NULL);

   if (psAnaExpertMode_ == 1)
   {
      printAsterisks(0);
      printf("\n");
      printf("* For each input, this RSMSobolTSIAnalyzer creates M sample\n");
      printf("* points by varying all other inputs. For each of the M \n");
      printf("* sample points, K levels of each input is used.\n");
      printf("* The total sample size = nInputs * M * K.\n");
      printf("* As a user, you are asked to decide on M and K.\n");
      printEquals(0);

      sprintf(pString,"Enter M (suggestion: > 2000) : ");
      nSubSamples = getInt(1000,1000000,pString);
      if (nSubSamples > 100000)
         printf("A M of %d may take very long time.\n", nSubSamples);

      sprintf(pString,"Enter K (suggestion: > 500) : ");
      nLevels = getInt(100,100000,pString);
      if (nLevels > 5000)
         printf("A K of %d may take very long time.\n", nLevels);
      printAsterisks(0);
   }
   else
   {  
      nSubSamples = 2000;
      nLevels = 500;
      if (printLevel > 0)
      {
         printf("RSMSobolTSI: default M = 2000.\n");
         printf("RSMSobolTSI: default K = 500.\n");
         printf("To change these settings, turn on ana_expert mode and rerun.\n");
      }
   }

   nSamp = 100000;
   if (printLevel > 1)
   {
      printf("RSMSobolTSI INFO: creating a sample for basic statistics.\n");
      printf("                  sample size = %d\n", nSamp);
   }

   sampleInputs  = new double[nSamp*nInputs];
   sampleOutputs = new double[nSamp];

   pdfman = new PDFManager();
   pdfman->initialize(nInputs,pdfFlags,inputMeans,
                      inputStdevs,*corMatp);
   vecLB.load(nInputs, xLower);
   vecUB.load(nInputs, xUpper);
   vecOut.setLength(nSamp*nInputs);
   pdfman->genSample(nSamp, vecOut, vecLB, vecUB);
   for (ii = 0; ii < nSamp*nInputs; ii++) sampleInputs[ii] = vecOut[ii];

   if (printLevel > 1)
   {
      printf("RSMSobolTSI: running the sample with response surface...\n");
   }
   faPtr->evaluatePoint(nSamp, sampleInputs, sampleOutputs);
   if (printLevel > 1)
   {
      printf("RSMSobolTSI: done running the sample with response surface.\n");
   }

   for (ii = 0; ii < nSamp; ii++)
   {
      oneSamplePt = &(sampleInputs[nInputs*ii]);
      ddata = constrPtr->evaluate(oneSamplePt,sampleOutputs[ii],status);
      if (status == 0) sampleOutputs[ii] = PSUADE_UNDEFINED;
   }

   dmean = 0.0;
   sCnt = 0;
   for (ii = 0; ii < nSamp; ii++)
   {
      if (sampleOutputs[ii] != PSUADE_UNDEFINED)
      {
         dmean += sampleOutputs[ii];
         sCnt++;
      }
   }
   if (sCnt > 1) dmean /= (double) sCnt;
   else
   {
      printf("RSMSobolTSI ERROR: too few samples that satisify the ");
      printf("constraints (%d out of %d).\n", sCnt, nSubSamples);
      delete [] sampleInputs;
      delete [] sampleOutputs;
      delete faPtr;
      delete pdfman;
      return PSUADE_UNDEFINED;
   }
   variance = 0.0;
   for (ii = 0; ii < nSamp; ii++)
   {
      if (sampleOutputs[ii] != PSUADE_UNDEFINED)
         variance += (sampleOutputs[ii] - dmean) *
                     (sampleOutputs[ii] - dmean) ;
   }
   variance /= (double) sCnt;
   printf("RSMSobolTSI: mean     (based on N = %d) = %10.3e\n",
          sCnt, dmean);
   printf("RSMSobolTSI: std dev. (based on N = %d) = %10.3e\n",
          sCnt, sqrt(variance));
   if (variance == 0.0) variance = 1.0;
   delete [] sampleInputs;
   delete [] sampleOutputs;

   cLower = new double[nInputs];
   cUpper = new double[nInputs];
   sampleInputs  = new double[nSubSamples*nInputs];
   Y2            = new double[nLevels];
   means         = new double[nSubSamples];
   bins          = new int[nSubSamples];
   vars          = new double[nSubSamples];
   tsi           = new double[nInputs];
   tsiRes        = new double[nInputs];
   mSamplePts    = new double[nLevels*nInputs];
   samplePts1D   = new double[nLevels];
   inputMeans1   = new double[nInputs];
   inputStdevs1  = new double[nInputs];
   pdfFlags1     = new int[nInputs];

   for (ii = 0; ii < nInputs; ii++)
   {
      if (printLevel > 1)
         printf("RSMSobolTSI: processing input %d\n", ii+1);

      corFlag = 0;
      for (jj = 0; jj < nInputs; jj++)
         if (ii != jj && corMatp->getEntry(ii,jj) != 0.0) corFlag = 1;
      if (corFlag == 0)
      {
         for (jj = 0; jj < nInputs; jj++) corFlag += pdfFlags[jj];
         if (corFlag == 0) corFlag = 1;
      }
      
      if (corFlag == 0)
      {
         if (printLevel > 1)
            printf("RSMSobolTSI: create samples (1)\n");
         corMat1.setDim(nInputs-1, nInputs-1);
         for (jj = 0; jj < ii; jj++)
         {
            cLower[jj] = xLower[jj];
            cUpper[jj] = xUpper[jj];
            pdfFlags1[jj] = pdfFlags[jj];
            inputMeans1[jj] = inputMeans[jj];
            inputStdevs1[jj] = inputStdevs[jj];
            for (kk = 0; kk < ii; kk++)
               corMat1.setEntry(jj, kk, corMatp->getEntry(jj,kk));
            for (kk = ii+1; kk < nInputs; kk++)
               corMat1.setEntry(jj, kk-1, corMatp->getEntry(jj,kk));
         }
         for (jj = ii+1; jj < nInputs; jj++)
         {
            cLower[jj-1] = xLower[jj];
            cUpper[jj-1] = xUpper[jj];
            pdfFlags1[jj-1] = pdfFlags[jj];
            inputMeans1[jj-1] = inputMeans[jj];
            inputStdevs1[jj-1] = inputStdevs[jj];
            for (kk = 0; kk < ii; kk++)
               corMat1.setEntry(jj-1, kk, corMatp->getEntry(jj,kk));
            for (kk = ii+1; kk < nInputs; kk++)
               corMat1.setEntry(jj-1, kk-1, corMatp->getEntry(jj,kk));
         }
         pdfman1 = new PDFManager();
         pdfman1->initialize(nInputs-1,pdfFlags1,inputMeans1,
                             inputStdevs1,corMat1);
         vecLB.load(nInputs-1, cLower);
         vecUB.load(nInputs-1, cUpper);
         vecOut.setLength(nSubSamples*(nInputs-1));
         pdfman1->genSample(nSubSamples, vecOut, vecLB, vecUB);
         for (jj = 0; jj < nSubSamples*(nInputs-1); jj++)
            sampleInputs[jj] = vecOut[jj];
         delete pdfman1;
         corMat2.setDim(1, 1);
         corMat2.setEntry(0, 0, corMatp->getEntry(0,0));
         pdfman2 = new PDFManager();
         pdfman2->initialize(1, &pdfFlags[ii], &inputMeans[ii],
                             &inputStdevs[ii],corMat2);
         vecLB.load(1, &xLower[ii]);
         vecUB.load(1, &xUpper[ii]);
         vecOut.setLength(nLevels);
         pdfman2->genSample(nLevels, vecOut, vecLB, vecUB);
         for (iL = 0; iL < nLevels; iL++) samplePts1D[iL] = vecOut[iL];
         delete pdfman2;
      }
      else
      {
         if (printLevel > 1)
            printf("RSMSobolTSI: create samples (2)\n");
         if (nInputs > 51)
            sampler = SamplingCreateFromID(PSUADE_SAMP_LHS);
         else
            sampler = SamplingCreateFromID(PSUADE_SAMP_LPTAU);
         for (jj = 0; jj < ii; jj++)
         {
            cLower[jj] = xLower[jj];
            cUpper[jj] = xUpper[jj];
         }
         for (jj = ii+1; jj < nInputs; jj++)
         {
            cLower[jj-1] = xLower[jj];
            cUpper[jj-1] = xUpper[jj];
         }
         sampler->setInputBounds(nInputs-1, cLower, cUpper);
         sampler->setOutputParams(1);
         sampler->setSamplingParams(nSubSamples, 1, 0);
         sampler->initialize(0);
         sampleInputs  = new double[nSubSamples*(nInputs-1)];
         sampleOutputs = new double[nSubSamples];
         sampleStates  = new int[nSubSamples];
         sampler->getSamples(nSubSamples, nInputs-1, 1, sampleInputs,
                             sampleOutputs, sampleStates);
         delete [] sampleStates;
         delete [] sampleOutputs;
         delete sampler;
         for (iL = 0; iL < nLevels; iL++)
            samplePts1D[iL] = (xUpper[ii] - xLower[ii]) / (nLevels-1) * 
                              iL + xLower[ii];
      }

      for (jj = 0; jj < nSubSamples; jj++)
      {
         for (iL = 0; iL < nLevels; iL++)
         {
            for (kk = 0; kk < ii; kk++)
               mSamplePts[iL*nInputs+kk] = sampleInputs[jj*(nInputs-1)+kk];
            for (kk = ii+1; kk < nInputs; kk++)
               mSamplePts[iL*nInputs+kk] = sampleInputs[jj*(nInputs-1)+kk-1];
            mSamplePts[iL*nInputs+ii] = samplePts1D[iL];
         } 

         if (corFlag == 1)
         {
            vecIn.load(nLevels*nInputs, mSamplePts);
            vecOut.setLength(nLevels*nInputs);
            pdfman->invCDF(nLevels, vecIn, vecOut, vecLB, vecUB);
            for (kk = 0; kk < nLevels*nInputs; kk++) 
               mSamplePts[kk] = vecOut[kk];
         }

         faPtr->evaluatePoint(nLevels,mSamplePts,Y2);

         for (iL = 0; iL < nLevels; iL++)
         {
            oneSamplePt = &(mSamplePts[iL*nInputs]);
            ddata = constrPtr->evaluate(oneSamplePt, Y2[iL],status);
            if (status == 0) Y2[iL] = PSUADE_UNDEFINED;

            if (ddata == PSUADE_UNDEFINED) Y2[iL] = PSUADE_UNDEFINED;
         } 

         means[jj] = 0.0;
         vars[jj] = 0.0;
         sCnt = 0;
         for (iL = 0; iL < nLevels; iL++)
         {
            if (Y2[iL] != PSUADE_UNDEFINED)
            {
               means[jj] += Y2[iL];
               sCnt++;
            }
         }
         bins[jj] = sCnt;
         if (sCnt >= 1) means[jj] /= (double) sCnt;
         else           means[jj] = PSUADE_UNDEFINED;
         if (means[jj] == PSUADE_UNDEFINED) vars[jj] = PSUADE_UNDEFINED;
         else
         {
            for (iL = 0; iL < nLevels; iL++)
            {
               if (Y2[iL] != PSUADE_UNDEFINED)
                  vars[jj] += (Y2[iL] - means[jj]) * (Y2[iL] - means[jj]);
            }
            if (sCnt == 1) vars[jj] = 0.0;
            else           vars[jj] = vars[jj] / (double) sCnt;
         }
         if (sCnt < nSubSamples/10 && printLevel >= 5)
            printf("RSMSobolTSI WARNING: subsample size = %d\n", sCnt);
      }

      totalCnt = 0;
      for (jj = 0; jj < nSubSamples; jj++) totalCnt += bins[jj];
      dvar = 0.0;
      for (jj = 0; jj < nSubSamples; jj++)
      {
         if (vars[jj] != PSUADE_UNDEFINED)
            dvar += vars[jj] * bins[jj] / totalCnt;
      }
      tsi[ii] = dvar;

      dmean = 0.0;
      for (jj = 0; jj < nSubSamples; jj++)
      {
         if (means[jj] != PSUADE_UNDEFINED)
            dmean += means[jj] * bins[jj] / totalCnt;
      }
      dvar = 0.0;
      for (jj = 0; jj < nSubSamples; jj++)
         if (means[jj] != PSUADE_UNDEFINED)
            dvar += (means[jj] - dmean) * (means[jj] - dmean) * 
                    bins[jj] / totalCnt;
      tsiRes[ii] = dvar;

      if (printLevel > 1)
      {
         printf("RSMSobolTSI (unnormalized) for input %3d = %12.4e\n",
                ii+1, tsi[ii]);
         printf("RSMSobolTSI (  normalized) for input %3d = %12.4e\n",
                ii+1, tsi[ii]/variance);
      }
   }
   printAsterisks(0);
   for (ii = 0; ii < nInputs; ii++) tsi[ii] /= variance;
   for (ii = 0; ii < nInputs; ii++)
      printf("RSMSobolTSI (normalized) for input %3d = %12.4e\n",ii+1,tsi[ii]);
   printEquals(0);
   if (printLevel > 1)
   {
      for (ii = 0; ii < nInputs; ii++) means[ii] = (double) ii;
      sortDbleList2(nInputs, tsi, means);
      for (ii = nInputs-1; ii >= 0; ii--)
         printf("RSMSobolTSI (normalized,ordered) for input %3d = %12.4e\n",
                (int) means[ii]+1,tsi[ii]);
      printAsterisks(0);
   }
   if (printLevel > 2)
   {
      printAsterisks(0);
      for (ii = 0; ii < nInputs; ii++)
         printf("RSMSobolTSI residual (normalized) for input %3d = %12.4e\n",
                ii+1, tsiRes[ii]/variance);
   }
   printAsterisks(0);
    
   pData *pPtr = ioPtr->getAuxData();
   pPtr->nDbles_ = nInputs;
   pPtr->dbleArray_ = new double[nInputs];
   for (ii = 0; ii < nInputs; ii++) pPtr->dbleArray_[ii] = tsi[ii] * variance;
   pPtr->dbleData_ = variance;

   delete constrPtr;
   delete faPtr;
   delete [] cLower;
   delete [] cUpper;
   delete [] sampleInputs;
   delete [] Y2;
   delete [] Y;
   delete [] vars;
   delete [] means;
   delete [] tsi;
   delete [] tsiRes;
   delete [] bins;
   delete [] mSamplePts;
   delete pdfman;
   delete [] inputMeans1;
   delete [] inputStdevs1;
   delete [] pdfFlags1;
   return 0.0;
}

