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

#include "PsuadeUtil.h"
#include "sysdef.h"
#include "Matrix.h"
#include "Vector.h"
#include "Psuade.h"
#include "FuncApprox.h"
#include "Sampling.h"
#include "RSConstraints.h"
#include "PDFManager.h"
#include "PDFNormal.h"
#include "PsuadeData.h"
#include "PsuadeConfig.h"
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
   int    corFlag, method=1, noPDF;
   double *xLower, *xUpper, *X, *Y, *cLower, *cUpper, *sampleInputs;
   double *sampleOutputs, *oneSamplePt, *Y2, *means,*tsi, variance;
   double *vars, dmean, dvar, *tsiRes, ddata, *mSamplePts, *inputMeans;
   double *inputStdevs, *inputMeans1, *inputStdevs1, *samplePts1D, *Y3;
   char   pString[500], *cString, winput1[500], winput2[500];
   Sampling      *sampler;
   PsuadeData    *ioPtr;
   FuncApprox    *faPtr;
   RSConstraints *constrPtr;
   Vector        vecIn, vecOut, vecUB, vecLB;
   pData         pCorMat;
   Matrix        *corMatp, corMat, corMat1, corMat2;
   PDFManager    *pdfman, *pdfman1, *pdfman2;

   if (method == 1)
   {
      analyze3(adata);
      return 0;
   }

   printAsterisks(0);
   printf("*          RS-based Total Order Sobol' Indices \n");
   printEquals(0); 
   printf("* TO GAIN ACCESS TO DIFFERENT OPTIONS: SET \n");
   printf("*\n");
   printf("* - ana_expert mode to finetune RSMSobolTSI parameters, \n");
   printf("*   (e.g. sample size for integration can be adjusted).\n");
   printf("* - rs_expert to mode finetune response surface for RSMSobolTSI,\n");
   printf("* - printlevel to 1 to display more information.\n");
   printf("* - ntimes to 100 to compute also error bars for the results\n");
   printEquals(0);

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
   noPDF = 1;
   if (pdfFlags != NULL)
   {
      for (ii = 0; ii < nInputs; ii++)
         if (pdfFlags[ii] != 0) noPDF = 0;
   }
   if (noPDF == 1) printf("* RSMSobolTSI INFO: all uniform distributions.\n");
   else
   {
      printf("* RSMSobolTSI INFO: non-uniform distributions detected\n");
      printf("                    which will be used in this analysis.\n");
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
            printf("* RSMSobolTSI INFO: this method cannot handle correlated\n");
            printf("*           inputs using joint PDFs yet. Use group\n");
            printf("*           variance-based method.\n");
            return PSUADE_UNDEFINED;
         }
      }
   }
   status = 0;
   for (ii = 0; ii < nSamples; ii++)
      if (Y3[nOutputs*ii+outputID] > 0.9*PSUADE_UNDEFINED) status = 1;
   if (status == 1)
   {
      printf("RSMSobolTSI ERROR: Some outputs are undefined. Prune the\n");
      printf("                   undefined sample points first.\n");
      return PSUADE_UNDEFINED;
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
      nSubSamples = 1000;
      nLevels = 100;
      printf("* nInputs m = %d\n", nInputs);
      printf("* default M = %d\n", nSubSamples);
      printf("* default K = %d\n", nLevels);
      printf("* As a user, you are asked to decide on M and K.\n");
      printEquals(0);
      sprintf(pString,"Enter M (suggestion: >= 1000) : ");
      nSubSamples = getInt(1000,100000,pString);
      sprintf(pString,"Enter K (suggestion: >= 100) : ");
      nLevels = getInt(100,1000,pString);
      printAsterisks(0);
   }
   else
   {  
      nSubSamples = 1000;
      nLevels = 100;
      if (psConfig_ != NULL)
      {
         cString = psConfig_->getParameter("RSMSoboltsi_nsubsamples");
         if (cString != NULL)
         {
            sscanf(cString, "%s %s %d", winput1, winput2, &nSubSamples);
            if (nSubSamples < 1000)
            {
               printf("RSMSobolTSI INFO: nSubSamples should be >= 1000.\n");
               nSubSamples = 1000;
            }
            else
            {
               printf("RSMSobolTSI INFO: nSubSamples = %d (config).\n",nSubSamples);
            }
         }
         cString = psConfig_->getParameter("RSMSoboltsi_nlevels");
         if (cString != NULL)
         {
            sscanf(cString, "%s %s %d", winput1, winput2, &nLevels);
            if (nLevels < 100)
            {
               printf("RSMSobolTSI INFO: nLevels should be >= 100.\n");
               nLevels = 100;
            }
            else
            {
               printf("RSMSobolTSI INFO: nLevels = %d (config).\n",nLevels);
            }
         }
      }
      if (printLevel > 0)
      {
         printAsterisks(0);
         printf("\n");
         printf("* For each input, this RSMSobolTSIAnalyzer creates M sample\n");
         printf("* points by varying all other inputs. For each of the M \n");
         printf("* sample points, K levels of each input is used.\n");
         printf("* The total sample size = nInputs * M * K.\n");
         printf("* nInputs m = %d\n", nInputs);
         printf("* default M = %d\n", nSubSamples);
         printf("* default K = %d\n", nLevels);
         printf("To change these settings, turn on ana_expert mode and rerun.\n");
         printEquals(0);
      }
   }

   nSamp = 100000;
   if (printLevel > 1)
   {
      printf("* RSMSobolTSI INFO: creating a sample for basic statistics.\n");
      printf("*                   sample size = %d\n", nSamp);
   }

   sampleInputs  = new double[nSamp*nInputs];
   sampleOutputs = new double[nSamp];

   pdfman = new PDFManager();
   if (pdfFlags == NULL)
   {
      printf("pdfFlags is NULL in file %s line %d aborting\n",__FILE__,
             __LINE__);
      abort();
   }
   pdfman->initialize(nInputs,pdfFlags,inputMeans,
                      inputStdevs,*corMatp);
   vecLB.load(nInputs, xLower);
   vecUB.load(nInputs, xUpper);
   vecOut.setLength(nSamp*nInputs);
   pdfman->genSample(nSamp, vecOut, vecLB, vecUB);
   for (ii = 0; ii < nSamp*nInputs; ii++) sampleInputs[ii] = vecOut[ii];

   if (printLevel > 1)
   {
      printf("* RSMSobolTSI: running the sample with response surface...\n");
   }
   faPtr->evaluatePoint(nSamp, sampleInputs, sampleOutputs);
   if (printLevel > 1)
   {
      printf("* RSMSobolTSI: done running the sample with response surface.\n");
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
   printf("* RSMSobolTSI: mean     (based on N = %d) = %10.3e\n",
          sCnt, dmean);
   printf("* RSMSobolTSI: std dev. (based on N = %d) = %10.3e\n",
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

// ************************************************************************
// perform analysis (version with error bars 12/2012)
// ------------------------------------------------------------------------
double RSMSobolTSIAnalyzer::analyze3(aData &adata)
{
   int    nInputs, nOutputs, nSamples, ii, jj, kk, *S, status, outputID;
   int    length, nSubSamples=20000, *sampleStates, iL, sCnt, nLevels=200;
   int    *pdfFlags, printLevel, *bins, totalCnt, nSamp, *pdfFlags1;
   int    corFlag, ntimes=1, noPDF;
   double *xLower, *xUpper, *X, *Y, *cLower, *cUpper, *sampleInputs;
   double *sampleOutputs, *oneSamplePt, *means,*tsi, variance;
   double *vars, dmean, dvar, *tsiRes, ddata, *mSamplePts, *inputMeans;
   double *inputStdevs, *inputMeans1, *inputStdevs1, *samplePts1D, *Y3;
   char   pString[500], *cString, winput1[500], winput2[500];
   Sampling      *sampler;
   PsuadeData    *ioPtr;
   FuncApprox    *faPtr;
   RSConstraints *constrPtr;
   Vector        vecIn, vecOut, vecUB, vecLB;
   pData         pCorMat;
   Matrix        *corMatp, corMat, corMat1, corMat2;
   PDFManager    *pdfman, *pdfman1, *pdfman2;

   printAsterisks(0);
   printf("*          RS-based Total Order Sobol' Indices \n");
   printEquals(0); 
   printf("* TO GAIN ACCESS TO DIFFERENT OPTIONS: SET \n");
   printf("*\n");
   printf("* - ana_expert mode to finetune RSMSobolTSI parameters, \n");
   printf("*   (e.g. sample size for integration can be adjusted).\n");
   printf("* - rs_expert to mode finetune response surface for RSMSobolTSI,\n");
   printf("* - printlevel to 1 to display more information.\n");
   printf("* - ntimes to 100 to compute also error bars for the results\n");
   printEquals(0);

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
   noPDF = 1;
   if (pdfFlags != NULL)
   {
      for (ii = 0; ii < nInputs; ii++)
         if (pdfFlags[ii] != 0) noPDF = 0;
   }
   if (noPDF == 1) printf("* RSMSobolTSI INFO: all uniform distributions.\n");
   else
   {
      printf("* RSMSobolTSI INFO: non-uniform distributions detected\n");
      printf("                    which will be used in this analysis.\n");
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
            printf("            inputs using joint PDFs yet. Use group\n");
            printf("            variance-based method.\n");
            return PSUADE_UNDEFINED;
         }
      }
   }
   status = 0;
   for (ii = 0; ii < nSamples; ii++)
      if (Y3[nOutputs*ii+outputID] > 0.9*PSUADE_UNDEFINED) status = 1;
   if (status == 1)
   {
      printf("RSMSobolTSI ERROR: Some outputs are undefined. Prune the\n");
      printf("                   undefined sample points first.\n");
      return PSUADE_UNDEFINED;
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
      nSubSamples = 1000;
      nLevels = 100;
      printf("* nInputs m = %d\n", nInputs);
      printf("* default M = %d\n", nSubSamples);
      printf("* default K = %d\n", nLevels);
      printf("* As a user, you are asked to decide on M and K.\n");
      printf("* Note: large M and K may take a long time to compute.\n");
      printEquals(0);
      sprintf(pString,"Enter M (1000 - 100000) : ");
      nSubSamples = getInt(1000,100000,pString);
      sprintf(pString,"Enter K (100 - 1000) : ");
      nLevels = getInt(100,1000,pString);
      printf("* This analysis is to be done ntimes=%d times\n",ntimes);
      printf("* to also error bars. If you do not need error bars,\n");
      printf("* Keep ntimes=1. Otherwise use 100-500.\n");
      sprintf(pString,"Enter ntimes (suggestion: 1 or 100 - 500) : ");
      ntimes = getInt(1, 500, pString);
      printAsterisks(0);
   }
   else
   {  
      nSubSamples = 1000;
      nLevels = 100;
      if (psConfig_ != NULL)
      {
         cString = psConfig_->getParameter("RSMSoboltsi_nsubsamples");
         if (cString != NULL)
         {
            sscanf(cString, "%s %s %d", winput1, winput2, &nSubSamples);
            if (nSubSamples < 1000)
            {
               printf("RSMSobolTSI INFO: nSubSamples should be >= 1000.\n");
               nSubSamples = 1000;
            }
            else
            {
               printf("RSMSobolTSI INFO: nSubSamples = %d (config).\n",nSubSamples);
            }
         }
         cString = psConfig_->getParameter("RSMSoboltsi_nlevels");
         if (cString != NULL)
         {
            sscanf(cString, "%s %s %d", winput1, winput2, &nLevels);
            if (nLevels < 100)
            {
               printf("RSMSobolTSI INFO: nLevels should be >= 100.\n");
               nLevels = 100;
            }
            else
            {
               printf("RSMSobolTSI INFO: nLevels = %d (config).\n",nLevels);
            }
         }
      }
      if (printLevel > 0)
      {
         printAsterisks(0);
         printf("\n");
         printf("* For each input, this RSMSobolTSIAnalyzer creates M sample\n");
         printf("* points by varying all other inputs. For each of the M \n");
         printf("* sample points, K levels of each input is used.\n");
         printf("* The total sample size = nInputs * M * K.\n");
         printf("* nInputs m = %d\n", nInputs);
         printf("* default M = %d\n", nSubSamples);
         printf("* default K = %d\n", nLevels);
         printf("To change these settings, turn on ana_expert mode and rerun.\n");
         printAsterisks(0);
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
   if (pdfFlags == NULL)
   {
      printf("pdfFlags is NULL in file %s line %d aborting\n", __FILE__, 
             __LINE__);
      abort();
   }
   pdfman->initialize(nInputs,pdfFlags,inputMeans,
                      inputStdevs,*corMatp);
   vecLB.load(nInputs, xLower);
   vecUB.load(nInputs, xUpper);
   vecOut.setLength(nSamp*nInputs);
   pdfman->genSample(nSamp, vecOut, vecLB, vecUB);
   for (ii = 0; ii < nSamp*nInputs; ii++) sampleInputs[ii] = vecOut[ii];

   double *sampleStdevs = new double[nSamp];
   if (printLevel > 1)
   {
      printf("RSMSobolTSI: running the sample with response surface...\n");
   }
   if (ntimes <= 1)
   {
      faPtr->evaluatePoint(nSamp,sampleInputs,sampleOutputs);
      for (ii = 0; ii < nSamp; ii++) sampleStdevs[ii] = 0.0;
   }
   else
   {
      faPtr->evaluatePointFuzzy(nSamp,sampleInputs,sampleOutputs,sampleStdevs);
   }
   if (printLevel > 1)
   {
      printf("RSMSobolTSI: done running the sample with response surface.\n");
   }

   int    nn, count;
   int    *Ycnts  = new int[ntimes];
   double *Ymeans = new double[ntimes];
   double *Ystds  = new double[ntimes];
   double *newY   = new double[ntimes];
   PDFNormal **rsPDFs = new PDFNormal*[nSamp];
   for (ii = 0; ii < nSamp; ii++)
   {
      if (sampleStdevs[ii] != 0) 
           rsPDFs[ii] = new PDFNormal(sampleOutputs[ii],sampleStdevs[ii]); 
      else rsPDFs[ii] = NULL;
   }
   for (nn = 0; nn < ntimes; nn++)
   {
      Ymeans[nn] = 0.0;
      Ystds[nn] = 0.0;
      Ycnts[nn] = 0;
   }
   for (ii = 0; ii < nSamp; ii++)
   {
      if (rsPDFs[ii] != NULL)
      {
         rsPDFs[ii]->genSample(ntimes,newY,sampleOutputs[ii]-4*sampleStdevs[ii],
                               sampleOutputs[ii]+4*sampleStdevs[ii]);
      }
      else 
      {
         for (nn = 0; nn < ntimes; nn++) newY[nn] = sampleOutputs[ii];
      }
      oneSamplePt = &(sampleInputs[ii*nInputs]);
      for (nn = 0; nn < ntimes; nn++)
      {
         ddata = constrPtr->evaluate(oneSamplePt,newY[nn],status);
         if (status != 0)
         {
            Ymeans[nn] += newY[nn];
            Ycnts[nn]++;
         }
      }
   }

   for (nn = 0; nn < ntimes; nn++)
      if (Ycnts[nn] != 0) Ymeans[nn] /= (double) Ycnts[nn];
   for (ii = 0; ii < nSamp; ii++)
   {
      if (rsPDFs[ii] != NULL)
      {
         rsPDFs[ii]->genSample(ntimes,newY,sampleOutputs[ii]-4*sampleStdevs[ii],
                               sampleOutputs[ii]+4*sampleStdevs[ii]);
      } 
      else 
      {
         for (nn = 0; nn < ntimes; nn++) newY[nn] = sampleOutputs[ii];
      }
      oneSamplePt = &(sampleInputs[ii*nInputs]);
      for (nn = 0; nn < ntimes; nn++)
      {
         ddata = constrPtr->evaluate(oneSamplePt,newY[nn],status);
         if (status != 0)
            Ystds[nn] += pow(newY[nn] - Ymeans[nn], 2.0);
      }
   }
   for (ii = 0; ii < nSamp; ii++)
      if (rsPDFs[ii] != NULL) delete rsPDFs[ii];
   delete rsPDFs;

   count  = 0;
   for (nn = 0; nn < ntimes; nn++) count += Ycnts[nn];
   printf("* RSMSobolTSI INFO: %6.2f percent passes the contrained filters.\n",
          (double) count * 100.0 /((double) ntimes*nSamp));
   if (100.0 * count / ((double) ntimes*nSamp) < 10.0)
   {
      printf("RSMSobolTSI ERROR: too few samples that satisify the ");
      printf("constraints (%d out of %d).\n", count, nSubSamples);
      delete [] sampleInputs;
      delete [] sampleOutputs;
      delete [] sampleStdevs;
      delete [] Ycnts;
      delete [] Ystds;
      delete [] Ymeans;
      delete [] newY;
      delete faPtr;
      delete pdfman;
      return PSUADE_UNDEFINED;
   }
   for (nn = 0; nn < ntimes; nn++)
   {
      if (Ycnts[nn] > 0) Ystds[nn] = sqrt(Ystds[nn]/nSamp);
      else               Ystds[nn] = 0.0;
   }
   double smean = 0.0;
   for (nn = 0; nn < ntimes; nn++) smean += Ystds[nn];
   smean /= (double) ntimes;
   double sstd = 0.0;
   for (nn = 0; nn < ntimes; nn++)
      sstd += (Ystds[nn] - smean) * (Ystds[nn] - smean) ;
   sstd = sqrt(sstd/ntimes);

   dmean = 0.0;
   for (nn = 0; nn < ntimes; nn++) dmean += Ymeans[nn];
   dmean /= (double) ntimes;
   variance = 0.0;
   for (nn = 0; nn < ntimes; nn++)
      variance += (Ymeans[nn] - dmean) * (Ymeans[nn] - dmean) ;
   variance /= (double) ntimes;

   printf("* RSMSobolTSI: sample mean (std dev of mean) = %10.3e (%10.3e)\n",
          dmean, sqrt(variance));
   printf("* RSMSobolTSI: std dev (std dev of std dev)  = %10.3e (%10.3e)\n",
          smean, sstd);
   if (smean == 0.0) smean = 1.0;
   delete [] sampleInputs;
   delete [] sampleOutputs;
   delete [] sampleStdevs;
   delete [] Ycnts;
   delete [] Ystds;
   delete [] Ymeans;
   delete [] newY;

   double *samplePtsND  = new double[nSubSamples*nInputs*nInputs];
   cLower = new double[nInputs];
   cUpper = new double[nInputs];
   sampleOutputs = new double[nSubSamples];
   sampleStates  = new int[nSubSamples];
   samplePts1D   = new double[nLevels*nInputs];
   inputMeans1   = new double[nInputs];
   inputStdevs1  = new double[nInputs];
   pdfFlags1     = new int[nInputs];

   for (ii = 0; ii < nInputs; ii++)
   {
      if (printLevel > 0)
         printf("RSMSobolTSI: creating sample for input %d\n", ii+1);

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
            samplePtsND[nSubSamples*ii*nInputs+jj] = vecOut[jj];
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
         for (iL = 0; iL < nLevels; iL++)
            samplePts1D[ii*nLevels+iL] = vecOut[iL];
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
         sampler->getSamples(nSubSamples, nInputs-1, 1, 
                             &(samplePtsND[ii*nSubSamples*nInputs]),
                             sampleOutputs, sampleStates);
         delete sampler;
         for (iL = 0; iL < nLevels; iL++)
            samplePts1D[iL+ii*nLevels] = (xUpper[ii]-xLower[ii])/(nLevels-1) * 
                              iL + xLower[ii];
      }
   }
   delete [] inputMeans1;
   delete [] inputStdevs1;
   delete [] pdfFlags1;
   delete [] cLower;
   delete [] cUpper;
   delete [] sampleStates;
   delete [] sampleOutputs;

   double *YZ      = new double[nLevels*ntimes];
   double *tsiMeds = new double[nInputs];
   double *tsiMins = new double[nInputs];
   double *tsiMaxs = new double[nInputs];
   double *YFuzzy = new double[nLevels];
   double *SFuzzy = new double[nLevels];
   PDFNormal *rsPDF;
   mSamplePts = new double[nLevels*nInputs];
   means      = new double[nSubSamples*ntimes];
   bins       = new int[nSubSamples*ntimes];
   vars       = new double[nSubSamples*ntimes];
   tsi        = new double[nInputs*ntimes];
   tsiRes     = new double[nInputs*ntimes];

   for (ii = 0; ii < nInputs; ii++)
   {
      if (printLevel > 0)
         printf("RSMSobolTSI: processing input %d\n", ii+1);
      for (jj = 0; jj < nSubSamples; jj++)
      {
         for (iL = 0; iL < nLevels; iL++)
         {
            for (kk = 0; kk < ii; kk++)
               mSamplePts[iL*nInputs+kk] = 
                  samplePtsND[ii*nSubSamples*nInputs+jj*(nInputs-1)+kk];
            for (kk = ii+1; kk < nInputs; kk++)
               mSamplePts[iL*nInputs+kk] = 
                  samplePtsND[ii*nSubSamples*nInputs+jj*(nInputs-1)+kk-1];
            mSamplePts[iL*nInputs+ii] = samplePts1D[iL+ii*nLevels];
         } 

         if (corFlag == 1)
         {
            vecIn.load(nLevels*nInputs, mSamplePts);
            vecOut.setLength(nLevels*nInputs);
            pdfman->invCDF(nLevels, vecIn, vecOut, vecLB, vecUB);
            for (kk = 0; kk < nLevels*nInputs; kk++) 
               mSamplePts[kk] = vecOut[kk];
         }

         if (ntimes == 1)
         {
            faPtr->evaluatePoint(nLevels,mSamplePts,YFuzzy);
            for (iL = 0; iL < nLevels; iL++) SFuzzy[iL] = 0.0;
         }
         else
         {
            faPtr->evaluatePointFuzzy(nLevels,mSamplePts,YFuzzy,SFuzzy);
         }
         for (iL = 0; iL < nLevels; iL++)
         {
            if (SFuzzy[iL] != 0)
            {
               rsPDF = new PDFNormal(YFuzzy[iL],SFuzzy[iL]);
               rsPDF->genSample(ntimes,&(YZ[iL*ntimes]),YFuzzy[iL]-4*SFuzzy[iL],
                                YFuzzy[iL]+4*SFuzzy[iL]);
               delete rsPDF;
            }
            else
            {
               for (nn = 0; nn < ntimes; nn++) YZ[iL*ntimes+nn] = YFuzzy[iL];
            }
         }

         for (iL = 0; iL < nLevels; iL++)
         {
            oneSamplePt = &(mSamplePts[iL*nInputs]);
            for (nn = 0; nn < ntimes; nn++)
            {
               ddata = constrPtr->evaluate(oneSamplePt,YZ[iL*ntimes+nn],status);
               if (status == 0) YZ[iL*ntimes+nn] = PSUADE_UNDEFINED;
            }
         }

         for (nn = 0; nn < ntimes; nn++)
         {
            means[jj*ntimes+nn] = 0.0;
            vars[jj*ntimes+nn] = 0.0;
            sCnt = 0;
            for (iL = 0; iL < nLevels; iL++)
            {
               if (YZ[iL*ntimes+nn] != PSUADE_UNDEFINED)
               {
                  means[jj*ntimes+nn] += YZ[iL*ntimes+nn];
                  sCnt++;
               }
            }
            bins[jj*ntimes+nn] = sCnt;
            if (sCnt >= 1) means[jj*ntimes+nn] /= (double) sCnt;
            else           means[jj*ntimes+nn] = PSUADE_UNDEFINED;
            if (means[jj*ntimes+nn] == PSUADE_UNDEFINED)
               vars[jj*ntimes+nn] = PSUADE_UNDEFINED;
            else
            {
               for (iL = 0; iL < nLevels; iL++)
               {
                  if (YZ[iL*ntimes+nn] != PSUADE_UNDEFINED)
                     vars[jj*ntimes+nn] += 
                       pow(YZ[iL*ntimes+nn]-means[jj*ntimes+nn],2.0);
               }
               if (sCnt == 1) vars[jj*ntimes+nn] = 0.0;
               else           vars[jj*ntimes+nn] /= (double) sCnt;
            }
         }
      }

      for (nn = 0; nn < ntimes; nn++)
      {
         totalCnt = 0;
         for (jj = 0; jj < nSubSamples; jj++) totalCnt += bins[jj*ntimes+nn];
         if (totalCnt == 0)
         {
            printf("RSMSobolTSI ERROR: no feasible region.\n");
            exit(1);
         }
         dvar = 0.0;
         for (jj = 0; jj < nSubSamples; jj++)
         {
            if (vars[jj*ntimes+nn] != PSUADE_UNDEFINED)
               dvar += vars[jj*ntimes+nn] * bins[jj*ntimes+nn] / totalCnt;
         }
         tsi[ii*ntimes+nn] = dvar;

         dmean = 0.0;
         for (jj = 0; jj < nSubSamples; jj++)
         {
            if (means[jj*ntimes+nn] != PSUADE_UNDEFINED)
               dmean += means[jj*ntimes+nn]*bins[jj*ntimes+nn]/totalCnt;
         }
         dvar = 0.0;
         for (jj = 0; jj < nSubSamples; jj++)
            if (means[jj*ntimes+nn] != PSUADE_UNDEFINED)
               dvar += pow(means[jj*ntimes+nn]-dmean, 2.0)*bins[jj]/totalCnt;
         tsiRes[ii*ntimes+nn] = dvar;
      }
   }

   for (ii = 0; ii < nInputs; ii++)
   {
      tsiMaxs[ii] = -PSUADE_UNDEFINED;
      tsiMins[ii] =  PSUADE_UNDEFINED;
      tsiMeds[ii] =  0.0;
      for (nn = 0; nn < ntimes; nn++)
      {
         if (tsi[ii*ntimes+nn] > tsiMaxs[ii]) tsiMaxs[ii] = tsi[ii*ntimes+nn];
         if (tsi[ii*ntimes+nn] < tsiMins[ii]) tsiMins[ii] = tsi[ii*ntimes+nn];
         tsiMeds[ii] += tsi[ii*ntimes+nn];
      }
      tsiMeds[ii] /= (double) ntimes;
      if (smean != 0)
      {
         tsiMeds[ii] /= (smean * smean);
         tsiMaxs[ii] /= (smean * smean);
         tsiMins[ii] /= (smean * smean);
      }
      else tsiMeds[ii] = tsiMaxs[ii] = tsiMins[ii] = 0.0;

      if (printLevel >= 0)
      {
         printf("RSMSobolTSI (normalized) for input %3d = %12.4e",
                ii+1,tsiMeds[ii]);
         if (ntimes > 1)
              printf(", bounds = [%12.4e, %12.4e]\n",tsiMins[ii],tsiMaxs[ii]);
         else printf("\n");
      }
   }
   if (ntimes == 1)
   {
      pData *pPtr = ioPtr->getAuxData();
      pPtr->nDbles_ = nInputs;
      pPtr->dbleArray_ = new double[nInputs];
      for (ii = 0; ii < nInputs; ii++) 
         pPtr->dbleArray_[ii] = tsiMeds[ii] * smean * smean;
      pPtr->dbleData_ = smean*smean;
   }
   if (printLevel > 1)
   {
      printEquals(0);
      for (ii = 0; ii < nInputs; ii++) means[ii] = (double) ii;
      sortDbleList2(nInputs, tsiMeds, means);
      for (ii = nInputs-1; ii >= 0; ii--)
         printf("RSMSobolTSI (normalized,ordered) for input %3d = %12.4e\n",
                (int) means[ii]+1,tsiMeds[ii]);
      printAsterisks(0);
   }
   printAsterisks(0);
    
   if (ntimes == 1)
   {
      for (ii = 0; ii < nInputs; ii++) tsiMeds[ii] *= smean * smean;
      printResults(nInputs, smean*smean, tsiMeds, NULL, NULL, ioPtr, 0);
   }
   else
   {
      for (ii = 0; ii < nInputs; ii++)
      {
         tsiMeds[ii] *= smean * smean;
         tsiMins[ii] *= smean * smean;
         tsiMaxs[ii] *= smean * smean;
      }
      printResults(nInputs,smean*smean,tsiMeds,tsiMins,tsiMaxs,ioPtr,1);
   }

   delete constrPtr;
   delete faPtr;
   delete [] Y;
   delete [] YZ;
   delete [] YFuzzy;
   delete [] SFuzzy;
   delete [] vars;
   delete [] means;
   delete [] tsi;
   delete [] tsiMeds;
   delete [] tsiMaxs;
   delete [] tsiMins;
   delete [] tsiRes;
   delete [] bins;
   delete [] mSamplePts;
   delete pdfman;
   return 0.0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
RSMSobolTSIAnalyzer& RSMSobolTSIAnalyzer::operator=(const RSMSobolTSIAnalyzer &)
{
   printf("RSMSobolTSI operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

// ************************************************************************
// print result
// ------------------------------------------------------------------------
int RSMSobolTSIAnalyzer::printResults(int nInputs, double variance,
                                      double *tsi, double *tsiMins, 
                                      double *tsiMaxs, PsuadeData *ioPtr,
                                      int flag)
{
   int   ii;
   FILE  *fp;
   char  **iNames;
   pData qData;

   if (ioPtr != NULL) ioPtr->getParameter("input_names", qData);
   if (qData.strArray_ != NULL) iNames = qData.strArray_;
   else                         iNames = NULL;
   printEquals(0);
   if (variance == 0.0)
   {
      printf("Total variance = 0. Hence, no total effect plot.\n");
      return 0;
   }

   printEquals(0);
   printf("Total Effect Statistics: \n");
   for (ii = 0; ii < nInputs; ii++)
   {
      printf("Input %4d: Sobol' total sensitivity = %12.4e",ii+1,tsi[ii]/variance);
      if (flag == 1)
           printf(", bounds = [%12.4e, %12.4e]\n",tsiMins[ii]/variance,
                  tsiMaxs[ii]/variance);
      else printf("\n");
   }
   if (psPlotTool_ == 1) fp = fopen("scilabrssoboltsi.sci", "w");
   else                  fp = fopen("matlabrssoboltsi.m", "w");
   if (fp != NULL)
   {
      if (psPlotTool_ == 1)
      {
         fprintf(fp, "// This file contains Sobol' total indices\n");
         fprintf(fp, "// set sortFlag = 1 and set nn to be the number\n");
         fprintf(fp, "// of inputs to display.\n");
      }
      else
      {
         fprintf(fp, "%% This file contains Sobol' total indices\n");
         fprintf(fp, "%% set sortFlag = 1 and set nn to be the number\n");
         fprintf(fp, "%% of inputs to display.\n");
      }
      fprintf(fp, "sortFlag = 0;\n");
      fprintf(fp, "nn = %d;\n", nInputs);
      fprintf(fp, "var = %e;\n", variance);
      fprintf(fp, "Mids = [\n");
      for (ii = 0; ii < nInputs; ii++) fprintf(fp,"%24.16e\n",tsi[ii]/variance);
      fprintf(fp, "];\n");
      if (flag == 1)
      {
         fprintf(fp, "Mins = [\n");
         for (ii = 0; ii < nInputs; ii++) fprintf(fp,"%24.16e\n",tsiMins[ii]/variance);
         fprintf(fp, "];\n");
         fprintf(fp, "Maxs = [\n");
         for (ii = 0; ii < nInputs; ii++) fprintf(fp,"%24.16e\n",tsiMaxs[ii]/variance);
         fprintf(fp, "];\n");
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
      fwritePlotCLF(fp);
      fprintf(fp, "if (sortFlag == 1)\n");
      if (psPlotTool_ == 1)
           fprintf(fp, "  [Mids, I2] = gsort(Mids);\n");
      else fprintf(fp, "  [Mids, I2] = sort(Mids,'descend');\n");
      if (flag == 1)
      {
         fprintf(fp, "  Maxs = Maxs(I2);\n");
         fprintf(fp, "  Mins = Mins(I2);\n");
      }
      fprintf(fp, "  Str  = Str(I2);\n");
      fprintf(fp, "  I2 = I2(1:nn);\n");
      fprintf(fp, "  Mids = Mids(1:nn);\n");
      if (flag == 1)
      {
         fprintf(fp, "  Maxs = Maxs(1:nn);\n");
         fprintf(fp, "  Mins = Mins(1:nn);\n");
      }
      fprintf(fp, "  Str  = Str(1:nn);\n");
      fprintf(fp, "end\n");
      if (flag == 1)
      {
         fprintf(fp, "ymin = min(Mins);\n");
         fprintf(fp, "ymax = max(Maxs);\n");
      }
      else
      {
         fprintf(fp, "ymin = min(Mids);\n");
         fprintf(fp, "ymax = max(Mids);\n");
      }
      fprintf(fp, "h2 = 0.05 * (ymax - ymin);\n");
      fprintf(fp, "bar(Mids,0.8);\n");
      if (flag == 1)
      {
         fprintf(fp, "for ii = 1:nn\n");
         fprintf(fp, "   if (ii == 1)\n");
         if (psPlotTool_ == 1)
              fprintf(fp, "      set(gca(),\"auto_clear\",\"off\")\n");
         else fprintf(fp, "      hold on\n");
         fprintf(fp, "   end;\n");
         fprintf(fp, "   XX = [ii ii];\n");
         fprintf(fp, "   YY = [Mins(ii) Maxs(ii)];\n");
         fprintf(fp, "   plot(XX,YY,'-ko','LineWidth',3.0,'MarkerEdgeColor',");
         fprintf(fp, "'k','MarkerFaceColor','g','MarkerSize',12)\n");
         fprintf(fp, "end;\n");
      }
      fwritePlotAxes(fp);
      fprintf(fp, "ymin=0;\n");
      if (psPlotTool_ == 1)
      {
         fprintf(fp, "a=gca();\n");
         fprintf(fp, "a.data_bounds=[0, 0; nn+1, ymax];\n");
         fprintf(fp, "newtick = a.x_ticks;\n");
         fprintf(fp, "newtick(2) = [1:nn]';\n");
         fprintf(fp, "newtick(3) = Str';\n");
         fprintf(fp, "a.x_ticks = newtick;\n");
         fprintf(fp, "a.x_label.font_size = 3;\n");
         fprintf(fp, "a.x_label.font_style = 4;\n");
      }
      else
      {
         fprintf(fp, "axis([0 nn+1 0 ymax])\n");
         fprintf(fp, "set(gca,'XTickLabel',[]);\n");
         fprintf(fp, "th=text(1:nn, repmat(ymin-0.07*(ymax-ymin),nn,1),Str,");
         fprintf(fp, "'HorizontalAlignment','left','rotation',90);\n");
         fprintf(fp, "set(th, 'fontsize', 12)\n");
         fprintf(fp, "set(th, 'fontweight', 'bold')\n");
      }
      fwritePlotTitle(fp,"Sobol Total Order Indices");
      fwritePlotYLabel(fp, "Sobol Indices");
      if (psPlotTool_ == 1)
           fprintf(fp, "      set(gca(),\"auto_clear\",\"on\")\n");
      else fprintf(fp, "      hold off\n");
      fclose(fp);
      if (psPlotTool_ == 1)
           printf("RSMSobolTSI plot file = scilabrssoboltsi.sci\n");
      else printf("RSMSobolTSI plot file = matlabrssoboltsi.m\n");
      return 0;
   }
   else
   {
      printf("RSMSobolTSI ERROR: cannot create rssoboltsi graphics file.\n");
      return 0;
   }
}

