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
#include "sysdef.h"
#include "PrintingTS.h"

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
RSMSobolTSIAnalyzer::RSMSobolTSIAnalyzer() : Analyzer(), nInputs_(0), 
                  outputMean_(0), outputStd_(0), totalSensitivity_(0)
{
   setName("RSMSOBOLTSI");
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
RSMSobolTSIAnalyzer::~RSMSobolTSIAnalyzer()
{
   if (totalSensitivity_) delete[] totalSensitivity_;
}

// ************************************************************************
// perform analysis
// ------------------------------------------------------------------------
double RSMSobolTSIAnalyzer::analyze(aData &adata)
{
   int    nInputs, nOutputs, nSamples, ii, jj, kk, *S, status, outputID;
   int    nSubSamples=10000, *sampleStates, iL, sCnt, nLevels=50;
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

   if (totalSensitivity_) delete[] totalSensitivity_;
   totalSensitivity_  = NULL;
   if (method == 1)
   {
      analyze3(adata);
      return 0;
   }

   printAsterisks(PL_INFO, 0);
   printOutTS(PL_INFO, "*          RS-based Total Order Sobol' Indices \n");
   printEquals(PL_INFO, 0); 
   printOutTS(PL_INFO,"* TO GAIN ACCESS TO DIFFERENT OPTIONS: SET \n");
   printOutTS(PL_INFO,"*\n");
   printOutTS(PL_INFO,
        "* - ana_expert mode to finetune RSMSobolTSI parameters, \n");
   printOutTS(PL_INFO,
        "*   (e.g. sample size for integration can be adjusted).\n");
   printOutTS(PL_INFO,
        "* - rs_expert mode to finetune response surface for RSMSobolTSI,\n");
   printOutTS(PL_INFO,"* - printlevel to 1 to display more information.\n");
   printOutTS(PL_INFO,
        "* - ntimes to 100 to compute also error bars for the results\n");
   printEquals(PL_INFO,0);

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
      for (ii = 0; ii < nInputs; ii++)
      {
         if (pdfFlags[ii] == PSUADE_PDF_USER)
         {
            printOutTS(PL_ERROR,
              "* RSMSobolTSI ERROR: S PDF type currently not supported.\n");
            return PSUADE_UNDEFINED;
         }
      }
   }
   if (noPDF == 1) 
      printOutTS(PL_INFO, "* RSMSobolTSI INFO: all uniform distributions.\n");
   else
   {
      printOutTS(PL_INFO,"* RSMSobolTSI INFO: non-uniform distributions\n");
      printOutTS(PL_INFO,"     detected which will be used in this analysis.\n");
   }

   if (nInputs <= 0 || nSamples <= 0 || nOutputs <= 0) 
   {
      printOutTS(PL_ERROR, "RSMSobolTSI ERROR: invalid arguments.\n");
      printOutTS(PL_ERROR, "   nInputs  = %d\n", nInputs);
      printOutTS(PL_ERROR, "   nOutputs = %d\n", nOutputs);
      printOutTS(PL_ERROR, "   nSamples = %d\n", nSamples);
      return PSUADE_UNDEFINED;
   } 
   if (nInputs <= 1)
   {
      printOutTS(PL_ERROR,
                 "RSMSobolTSI: nInputs<=1 does not need this analysis.\n");
      return PSUADE_UNDEFINED;
   }
   if (outputID >= nOutputs || outputID < 0)
   {
      printOutTS(PL_ERROR,"RSMSobolTSI ERROR: invalid output ID (%d).\n", 
                 outputID);
      return PSUADE_UNDEFINED;
   }
   if (ioPtr == NULL)
   {
      printOutTS(PL_ERROR, "RSMSobolTSI ERROR: no data.\n");
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
            printOutTS(PL_ERROR,
               "* RSMSobolTSI INFO: this method cannot handle correlated\n");
            printOutTS(PL_ERROR,
               "*           inputs using joint PDFs yet. Use group\n");
            printOutTS(PL_ERROR, "*           variance-based method.\n");
            return PSUADE_UNDEFINED;
         }
      }
   }
   status = 0;
   for (ii = 0; ii < nSamples; ii++)
      if (Y3[nOutputs*ii+outputID] > 0.9*PSUADE_UNDEFINED) status = 1;
   if (status == 1)
   {
      printOutTS(PL_ERROR, 
           "RSMSobolTSI ERROR: Some outputs are undefined. Prune the\n");
      printOutTS(PL_ERROR,
           "                   undefined sample points first.\n");
      return PSUADE_UNDEFINED;
   }
   Y = new double[nSamples];
   for (ii = 0; ii < nSamples; ii++) Y[ii] = Y3[ii*nOutputs+outputID];

   constrPtr = new RSConstraints();
   constrPtr->genConstraints(ioPtr);

   faPtr = genFAInteractive(ioPtr, 0);
   status = faPtr->initialize(X, Y);

   if (psAnaExpertMode_ == 1)
   {
      printAsterisks(PL_INFO, 0);
      printOutTS(PL_INFO, "\n");
      printOutTS(PL_INFO,
           "* For each input, this RSMSobolTSIAnalyzer creates M sample\n");
      printOutTS(PL_INFO,
           "* points by varying all other inputs. For each of the M \n");
      printOutTS(PL_INFO,"* sample points, K levels of each input is used.\n");
      printOutTS(PL_INFO,"* The total sample size = nInputs * M * K.\n");
      nSubSamples = 10000;
      nLevels = 50;
      printOutTS(PL_INFO,"* nInputs m = %d\n", nInputs);
      printOutTS(PL_INFO,"* default M = %d\n", nSubSamples);
      printOutTS(PL_INFO,"* default K = %d\n", nLevels);
      printOutTS(PL_INFO,"* As a user, you are asked to decide on M and K.\n");
      printEquals(PL_INFO, 0);
      sprintf(pString,"Enter M (suggestion: >= 1000) : ");
      nSubSamples = getInt(1000,100000,pString);
      sprintf(pString,"Enter K (suggestion: >= 50) : ");
      nLevels = getInt(50,1000,pString);
      printAsterisks(PL_INFO, 0);
   }
   else
   {  
      nSubSamples = 10000;
      nLevels = 50;
      if (psConfig_ != NULL)
      {
         cString = psConfig_->getParameter("RSMSoboltsi_nsubsamples");
         if (cString != NULL)
         {
            sscanf(cString, "%s %s %d", winput1, winput2, &nSubSamples);
            if (nSubSamples < 1000)
            {
               printOutTS(PL_INFO, 
                    "RSMSobolTSI INFO: nSubSamples should be >= 1000.\n");
               nSubSamples = 10000;
            }
            else
            {
               printOutTS(PL_INFO,
                     "RSMSobolTSI INFO: nSubSamples = %d (config).\n",
                     nSubSamples);
            }
         }
         cString = psConfig_->getParameter("RSMSoboltsi_nlevels");
         if (cString != NULL)
         {
            sscanf(cString, "%s %s %d", winput1, winput2, &nLevels);
            if (nLevels < 50)
            {
               printOutTS(PL_INFO,
                    "RSMSobolTSI INFO: nLevels should be >= 50.\n");
               nLevels = 50;
            }
            else
            {
               printOutTS(PL_INFO, 
                    "RSMSobolTSI INFO: nLevels = %d (config).\n",nLevels);
            }
         }
      }
      if (printLevel > 0)
      {
         printAsterisks(PL_INFO, 0);
         printOutTS(PL_INFO,"\n");
         printOutTS(PL_INFO,
              "* For each input, this RSMSobolTSIAnalyzer creates M sample\n");
         printOutTS(PL_INFO,
              "* points by varying all other inputs. For each of the M \n");
         printOutTS(PL_INFO,"* sample points, K levels of each input is used.\n");
         printOutTS(PL_INFO,"* The total sample size = nInputs * M * K.\n");
         printOutTS(PL_INFO,"* nInputs m = %d\n", nInputs);
         printOutTS(PL_INFO,"* default M = %d\n", nSubSamples);
         printOutTS(PL_INFO,"* default K = %d\n", nLevels);
         printOutTS(PL_INFO,
              "To change these settings, turn on ana_expert mode and rerun.\n");
         printEquals(PL_INFO, 0);
      }
   }

   nSamp = 100000;
   printOutTS(PL_INFO,
        "* RSMSobolTSI INFO: creating a sample for basic statistics.\n");
   printOutTS(PL_INFO,"*                   sample size = %d\n", nSamp);

   sampleInputs  = new double[nSamp*nInputs];
   sampleOutputs = new double[nSamp];

   pdfman = new PDFManager();
   if (pdfFlags == NULL)
   {
      printOutTS(PL_ERROR,"pdfFlags is NULL in file %s line %d aborting\n",
                 __FILE__, __LINE__);
      abort();
   }
   pdfman->initialize(nInputs,pdfFlags,inputMeans,inputStdevs,*corMatp,
                      NULL,NULL);
   vecLB.load(nInputs, xLower);
   vecUB.load(nInputs, xUpper);
   vecOut.setLength(nSamp*nInputs);
   pdfman->genSample(nSamp, vecOut, vecLB, vecUB);
   for (ii = 0; ii < nSamp*nInputs; ii++) sampleInputs[ii] = vecOut[ii];

   printOutTS(PL_INFO, 
        "* RSMSobolTSI: running the sample with response surface...\n");
   faPtr->evaluatePoint(nSamp, sampleInputs, sampleOutputs);
   printOutTS(PL_INFO,
        "* RSMSobolTSI: done running the sample with response surface.\n");

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
      printOutTS(PL_ERROR,
           "RSMSobolTSI ERROR: too few samples that satisify the ");
      printOutTS(PL_ERROR, "constraints (%d out of %d).\n", 
                 sCnt, nSubSamples);
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
   printOutTS(PL_INFO,
        "* RSMSobolTSI: mean     (based on N = %d) = %10.3e\n",sCnt,dmean);
   printOutTS(PL_INFO, "* RSMSobolTSI: std dev. (based on N = %d) = %10.3e\n",
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
      printOutTS(PL_INFO, "RSMSobolTSI: processing input %d\n", ii+1);

      corFlag = 0;
      for (jj = 0; jj < nInputs; jj++)
         if (ii != jj && corMatp->getEntry(ii,jj) != 0.0) corFlag = 1;
      if (corFlag == 0)
      {
         for (jj = 0; jj < nInputs; jj++) corFlag += pdfFlags[jj];
         if (corFlag == 0) corFlag = 1;
         else              corFlag = 0;
      }
      
      if (corFlag == 0)
      {
         printOutTS(PL_INFO, "RSMSobolTSI: create samples (1)\n");
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
                             inputStdevs1,corMat1,NULL,NULL);
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
                             &inputStdevs[ii],corMat2,NULL,NULL);
         vecLB.load(1, &xLower[ii]);
         vecUB.load(1, &xUpper[ii]);
         vecOut.setLength(nLevels);
         pdfman2->genSample(nLevels, vecOut, vecLB, vecUB);
         for (iL = 0; iL < nLevels; iL++) samplePts1D[iL] = vecOut[iL];
         delete pdfman2;
      }
      else
      {
         printOutTS(PL_INFO, "RSMSobolTSI: create samples (2)\n");
         if (nInputs > 51)
              sampler = SamplingCreateFromID(PSUADE_SAMP_LHS);
         else sampler = SamplingCreateFromID(PSUADE_SAMP_LPTAU);
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
            printOutTS(PL_DUMP,
                 "RSMSobolTSI WARNING: subsample size = %d\n", sCnt);
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

      printOutTS(PL_INFO,
           "RSMSobolTSI (unnormalized) for input %3d = %12.4e\n",
           ii+1,tsi[ii]);
      printOutTS(PL_INFO,
           "RSMSobolTSI (  normalized) for input %3d = %12.4e\n",
           ii+1, tsi[ii]/variance);

   }
   printAsterisks(PL_INFO, 0);
   for (ii = 0; ii < nInputs; ii++) tsi[ii] /= variance;
   for (ii = 0; ii < nInputs; ii++)
      printOutTS(PL_INFO,
           "RSMSobolTSI (normalized) for input %3d = %12.4e\n",
           ii+1,tsi[ii]);
   printEquals(PL_INFO, 0);
   if (printLevel > 1)
   {
      for (ii = 0; ii < nInputs; ii++)
         printOutTS(PL_INFO, 
              "RSMSobolTSI (unnormalized) for input %3d = %12.4e\n",
              ii+1, tsi[ii] * variance);
      for (ii = 0; ii < nInputs; ii++) means[ii] = (double) ii;
      sortDbleList2(nInputs, tsi, means);
      for (ii = nInputs-1; ii >= 0; ii--)
         printOutTS(PL_INFO,
              "RSMSobolTSI (normalized,ordered) for input %3d = %12.4e\n",
              (int) means[ii]+1,tsi[ii]);
      printAsterisks(PL_INFO, 0);
   }
   if (printLevel > 2)
   {
      printAsterisks(PL_INFO, 0);
      for (ii = 0; ii < nInputs; ii++)
         printOutTS(PL_INFO, 
              "RSMSobolTSI residual (normalized) for input %3d = %12.4e\n",
              ii+1, tsiRes[ii]/variance);
   }
   printAsterisks(PL_INFO, 0);
    
   pData *pPtr = ioPtr->getAuxData();
   pPtr->nDbles_ = nInputs;
   pPtr->dbleArray_ = new double[nInputs];
   for (ii = 0; ii < nInputs; ii++) pPtr->dbleArray_[ii] = tsi[ii]*variance;
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
   int    nSubSamples=10000, *sampleStates, iL, sCnt, nLevels=50;
   int    *pdfFlags, printLevel, *bins, totalCnt, nSamp, *pdfFlags1;
   int    corFlag, ntimes=1, noPDF, offset;
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

   printAsterisks(PL_INFO, 0);
   printOutTS(PL_INFO, "*          RS-based Total Order Sobol' Indices \n");
   printEquals(PL_INFO, 0); 
   printOutTS(PL_INFO,"* TO GAIN ACCESS TO DIFFERENT OPTIONS: SET \n");
   printOutTS(PL_INFO,"*\n");
   printOutTS(PL_INFO,"* - ana_expert mode to finetune internal parameters\n");
   printOutTS(PL_INFO,"*   (e.g. adjust sample size for integration).\n");
   printOutTS(PL_INFO,"* - rs_expert to mode finetune response surface\n");
   printOutTS(PL_INFO,"* - printlevel to 1 to display more information.\n");
   printOutTS(PL_INFO,"* - ntimes=100 to compute error bars for the results\n");
   printEquals(PL_INFO, 0);

   printLevel  = adata.printLevel_;
   nInputs     = adata.nInputs_;
   nInputs_    = nInputs;
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
      for (ii = 0; ii < nInputs; ii++)
      {
         if (pdfFlags[ii] == PSUADE_PDF_USER)
         {
            printOutTS(PL_ERROR,
                 "* RSMSobolTSI ERROR: S PDF not supported.\n");
            return PSUADE_UNDEFINED;
         }
      }
   }
   if (noPDF == 1) 
      printOutTS(PL_INFO,"* RSMSobolTSI INFO: all uniform distributions.\n");
   else
   {
      printOutTS(PL_INFO,"* RSMSobolTSI INFO: non-uniform distributions\n");
      printOutTS(PL_INFO,"* detected.\n");
   }

   if (nInputs <= 0 || nSamples <= 0 || nOutputs <= 0) 
   {
      printOutTS(PL_ERROR, "RSMSobolTSI ERROR: invalid arguments.\n");
      printOutTS(PL_ERROR, "   nInputs  = %d\n", nInputs);
      printOutTS(PL_ERROR, "   nOutputs = %d\n", nOutputs);
      printOutTS(PL_ERROR, "   nSamples = %d\n", nSamples);
      return PSUADE_UNDEFINED;
   } 
   if (nInputs <= 1)
   {
      printOutTS(PL_ERROR,"RSMSobolTSI: analysis not needed for nInputs<=1.\n");
      return PSUADE_UNDEFINED;
   }
   if (outputID >= nOutputs || outputID < 0)
   {
      printOutTS(PL_ERROR, "RSMSobolTSI ERROR: invalid output ID (%d).\n", 
                 outputID);
      return PSUADE_UNDEFINED;
   }
   if (ioPtr == NULL)
   {
      printOutTS(PL_ERROR, "RSMSobolTSI ERROR: no data.\n");
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
            printOutTS(PL_ERROR,
                 "RSMSobolTSI INFO: this method cannot handle\n");
            printOutTS(PL_ERROR,
                 "                  correlated inputs yet. \n");
            return PSUADE_UNDEFINED;
         }
      }
   }
   status = 0;
   for (ii = 0; ii < nSamples; ii++)
      if (Y3[nOutputs*ii+outputID] > 0.9*PSUADE_UNDEFINED) status = 1;
   if (status == 1)
   {
      printOutTS(PL_ERROR,
           "RSMSobolTSI ERROR: Some outputs are undefined. Prune\n");
      printOutTS(PL_ERROR,
           "                   the undefined sample points first.\n");
      return PSUADE_UNDEFINED;
   }
   Y = new double[nSamples];
   for (ii = 0; ii < nSamples; ii++) Y[ii] = Y3[ii*nOutputs+outputID];

   constrPtr = new RSConstraints();
   constrPtr->genConstraints(ioPtr);

   faPtr = genFAInteractive(ioPtr, 0);
   status = faPtr->initialize(X, Y);

   if (psAnaExpertMode_ == 1)
   {
      printAsterisks(PL_INFO, 0);
      printOutTS(PL_INFO,"\n");
      printOutTS(PL_INFO,"* For each input, this analyzer creates M sample\n");
      printOutTS(PL_INFO,"* points by varying all other inputs. For each of\n");
      printOutTS(PL_INFO,"* the M sample points, K levels of each input is\n");
      printOutTS(PL_INFO,"* used. The total sample size = nInputs * M * K.\n");
      nSubSamples = 10000;
      nLevels = 50;
      printOutTS(PL_INFO,"* nInputs m = %d\n", nInputs);
      printOutTS(PL_INFO,"* default M = %d\n", nSubSamples);
      printOutTS(PL_INFO,"* default K = %d\n", nLevels);
      printOutTS(PL_INFO,"* Please enter your desired M and K below.\n");
      printOutTS(PL_INFO,"* NOTE: large M and K may take a long time\n");
      printEquals(PL_INFO, 0);
      sprintf(pString,"Enter M (1000 - 50000) : ");
      nSubSamples = getInt(1000,510000,pString);
      sprintf(pString,"Enter K (50 - 500) : ");
      nLevels = getInt(50,500,pString);
      printOutTS(PL_INFO,"* To include response surface uncertainties in\n");
      printOutTS(PL_INFO,"* this analysis, the sensitivity calculation\n");
      printOutTS(PL_INFO,"* is to be repeated a number of times using\n");
      printOutTS(PL_INFO,"* different bootstrapped samples. Please specify\n");
      printOutTS(PL_INFO,"* the number of bootstrapped samples below.\n");
      printOutTS(PL_INFO,"* If you do not need error bars, set it to 1.\n");
      sprintf(pString,"Enter the number of bootstrapped samples (1 - 500) : ");
      ntimes = getInt(1, 500, pString);
      printAsterisks(PL_INFO, 0);
   }
   else
   {  
      nSubSamples = 10000;
      nLevels = 50;
      if (psConfig_ != NULL)
      {
         cString = psConfig_->getParameter("RSMSoboltsi_nsubsamples");
         if (cString != NULL)
         {
            sscanf(cString, "%s %s %d", winput1, winput2, &nSubSamples);
            if (nSubSamples < 1000)
            {
               printOutTS(PL_INFO,
                    "RSMSobolTSI INFO: nSubSamples should be >= 1000.\n");
               nSubSamples = 1000;
            }
            else
            {
               printOutTS(PL_INFO, 
                    "RSMSobolTSI INFO: nSubSamples = %d (config).\n",
                    nSubSamples);
            }
         }
         cString = psConfig_->getParameter("RSMSoboltsi_nlevels");
         if (cString != NULL)
         {
            sscanf(cString, "%s %s %d", winput1, winput2, &nLevels);
            if (nLevels < 50)
            {
               printOutTS(PL_INFO,"RSMSobolTSI INFO: nLevels should be >= 50.\n");
               nLevels = 50;
            }
            else
            {
               printOutTS(PL_INFO, "RSMSobolTSI INFO: nLevels = %d (config).\n",
                          nLevels);
            }
         }
      }
      if (printLevel > 0)
      {
         printAsterisks(PL_INFO, 0);
         printOutTS(PL_INFO,"\n");
         printOutTS(PL_INFO,"* For each input, this analyzer creates M sample\n");
         printOutTS(PL_INFO,"* points by varying all other inputs. For each of\n");
         printOutTS(PL_INFO,"* the M sample points, K levels of each input is\n");
         printOutTS(PL_INFO,"* used. The total sample size = nInputs * M * K.\n");
         printOutTS(PL_INFO,"* nInputs m = %d\n", nInputs);
         printOutTS(PL_INFO,"* default M = %d\n", nSubSamples);
         printOutTS(PL_INFO,"* default K = %d\n", nLevels);
         printOutTS(PL_INFO,
              "To change settings, re-run with ana_expert mode on.\n");
         printAsterisks(PL_INFO, 0);
      }
   }

   nSamp = 200000;

   printOutTS(PL_INFO,
        "RSMSobolTSI INFO: creating a sample for basic statistics.\n");
   printOutTS(PL_INFO,"                  sample size = %d\n", nSamp);

   sampleInputs  = new double[nSamp*nInputs];
   sampleOutputs = new double[nSamp];

   pdfman = new PDFManager();
   if (pdfFlags == NULL)
   {
      printOutTS(PL_INFO, "pdfFlags is NULL in file %s line %d aborting\n", 
                 __FILE__, __LINE__);
      abort();
   }
   pdfman->initialize(nInputs,pdfFlags,inputMeans,
                      inputStdevs,*corMatp,NULL,NULL);
   vecLB.load(nInputs, xLower);
   vecUB.load(nInputs, xUpper);
   vecOut.setLength(nSamp*nInputs);
   pdfman->genSample(nSamp, vecOut, vecLB, vecUB);
   for (ii = 0; ii < nSamp*nInputs; ii++) sampleInputs[ii] = vecOut[ii];

   double *sampleStdevs = new double[nSamp];

   if (printLevel > 1)
      printOutTS(PL_INFO,
           "RSMSobolTSI: response surface evaluation begins...\n");

   if (ntimes <= 1)
   {
      faPtr->evaluatePoint(nSamp,sampleInputs,sampleOutputs);
      for (ii = 0; ii < nSamp; ii++) sampleStdevs[ii] = 0.0;
   }
   else
   {
      faPtr->evaluatePointFuzzy(nSamp,sampleInputs,sampleOutputs,
                                sampleStdevs);
   }
   if (printLevel > 1)
      printOutTS(PL_INFO,
           "RSMSobolTSI: response surface evaluation ends...\n");

   dmean = 0.0;
   for (ii = 0; ii < nSamp; ii++) dmean += sampleOutputs[ii];
   dmean /= (double) nSamp;
   dvar = 0.0;
   for (ii = 0; ii < nSamp; ii++)
      dvar += (sampleOutputs[ii] - dmean) * (sampleOutputs[ii] - dmean) ;
   dvar /= (double) nSamp;
   if (printLevel > 1)
   {
      printOutTS(PL_INFO,"RSMSobolTSI: sample mean (N=%d) = %e.\n",
                 nSamp, dmean);
      printOutTS(PL_INFO,"RSMSobolTSI: sample std  (N=%d) = %e.\n",
                 nSamp, sqrt(dvar));
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
         rsPDFs[ii]->genSample(ntimes,newY,
                               sampleOutputs[ii]-4*sampleStdevs[ii],
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
         rsPDFs[ii]->genSample(ntimes,newY,
                               sampleOutputs[ii]-4*sampleStdevs[ii],
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
   printOutTS(PL_INFO,
        "* RSMSobolTSI INFO: %6.2f percent passes the contraints\n",
        (double) count * 100.0 /((double) ntimes*nSamp));
   if (100.0 * count / ((double) ntimes*nSamp) < 10.0)
   {
      printOutTS(PL_ERROR,
           "RSMSobolTSI ERROR: too few samples that satisify the ");
      printOutTS(PL_ERROR, 
           "constraints (%d out of %d).\n", count, nSubSamples);
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

   printOutTS(PL_INFO,
        "* RSMSobolTSI: sample mean (std dev of mean) = %10.3e (%10.3e)\n",
        dmean, sqrt(variance));
   printOutTS(PL_INFO,
        "* RSMSobolTSI: std dev (std dev of std dev)  = %10.3e (%10.3e)\n",
        smean, sstd);

   //save mean & std
   outputMean_ = dmean;
   outputStd_ = smean;
   //cout << outputMean_ << ", " << outputStd_ << endl;

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

      printOutTS(PL_INFO, "RSMSobolTSI: processing input %d\n", ii+1);

      corFlag = 0;
      for (jj = 0; jj < nInputs; jj++)
         if (ii != jj && corMatp->getEntry(ii,jj) != 0.0) corFlag = 1;
      if (corFlag == 0)
      {
         for (jj = 0; jj < nInputs; jj++) corFlag += pdfFlags[jj];
         if (corFlag == 0) corFlag = 1;
         else              corFlag = 0;
      }
      
      if (corFlag == 0)
      {
         printOutTS(PL_DETAIL,
              "RSMSobolTSI: create samples (no input correlation)\n");
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
                             inputStdevs1,corMat1,NULL,NULL);
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
                             &inputStdevs[ii],corMat2,NULL,NULL);
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
         printOutTS(PL_DETAIL,
              "RSMSobolTSI: create samples (has input correlation)\n");
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
   double *YFuzzy = new double[nLevels*nSubSamples];
   double *SFuzzy = new double[nLevels*nSubSamples];
   PDFNormal *rsPDF;
   mSamplePts = new double[nSubSamples*nLevels*nInputs];
   means      = new double[nSubSamples*ntimes];
   bins       = new int[nSubSamples*ntimes];
   vars       = new double[nSubSamples*ntimes];
   tsi        = new double[nInputs*ntimes];
   tsiRes     = new double[nInputs*ntimes];
   if (tsiRes == NULL)
   {
      printOutTS(PL_ERROR,"ERROR: memory allocation in file %s line %d\n",
                 __FILE__, __LINE__);
      abort();
   }
   for (ii = 0; ii < nInputs; ii++)
   {

      printOutTS(PL_DETAIL, "RSMSobolTSI: processing input %d (phase 2)\n", 
                 ii+1);
      for (jj = 0; jj < nSubSamples; jj++)
      {
         offset = jj * nLevels * nInputs;
         for (iL = 0; iL < nLevels; iL++)
         {
            for (kk = 0; kk < ii; kk++)
               mSamplePts[offset+iL*nInputs+kk] = 
                  samplePtsND[ii*nSubSamples*nInputs+jj*(nInputs-1)+kk];
            for (kk = ii+1; kk < nInputs; kk++)
               mSamplePts[offset+iL*nInputs+kk] = 
                  samplePtsND[ii*nSubSamples*nInputs+jj*(nInputs-1)+kk-1];
            mSamplePts[offset+iL*nInputs+ii] = samplePts1D[iL+ii*nLevels];
         } 

         if (corFlag == 1)
         {
            vecIn.load(nLevels*nInputs, &mSamplePts[offset]);
            vecOut.setLength(nLevels*nInputs);
            pdfman->invCDF(nLevels, vecIn, vecOut, vecLB, vecUB);
            for (kk = 0; kk < nLevels*nInputs; kk++) 
               mSamplePts[offset+kk] = vecOut[kk];
         }
      }

      if (ntimes == 1)
      {
         faPtr->evaluatePoint(nLevels*nSubSamples,mSamplePts,YFuzzy);
         for (iL = 0; iL < nLevels*nSubSamples; iL++) SFuzzy[iL] = 0.0;
      }
      else
      {
         faPtr->evaluatePointFuzzy(nLevels*nSubSamples,mSamplePts,
                                   YFuzzy,SFuzzy);
      }

      for (jj = 0; jj < nSubSamples; jj++)
      {
         offset = jj * nLevels;
         for (iL = 0; iL < nLevels; iL++)
         {
            if (SFuzzy[offset+iL] != 0)
            {
               rsPDF = new PDFNormal(YFuzzy[offset+iL],SFuzzy[offset+iL]);
               rsPDF->genSample(ntimes,&(YZ[iL*ntimes]),
                                YFuzzy[offset+iL]-4*SFuzzy[offset+iL],
                                YFuzzy[offset+iL]+4*SFuzzy[offset+iL]);
               delete rsPDF;
            }
            else
            {
               for (nn = 0; nn < ntimes; nn++) 
                  YZ[iL*ntimes+nn] = YFuzzy[offset+iL];
            }
         }

         for (iL = 0; iL < nLevels; iL++)
         {
            oneSamplePt = &(mSamplePts[jj*nLevels*nInputs+iL*nInputs]);
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
            printOutTS(PL_ERROR, "RSMSobolTSI ERROR: no feasible region.\n");
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

      printOutTS(PL_INFO, "RSMSobolTSI (normalized) for input %3d = %12.4e",
                 ii+1,tsiMeds[ii]);
      if (ntimes > 1)
           printOutTS(PL_INFO, ", bounds = [%12.4e, %12.4e]\n",
                      tsiMins[ii],tsiMaxs[ii]);
      else printOutTS(PL_INFO, "\n");
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
   printAsterisks(PL_INFO, 0);
    
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
   if (printLevel > 1)
   {
      printEquals(PL_INFO, 0);
      for (ii = 0; ii < nInputs; ii++)
         printOutTS(PL_INFO, 
              "RSMSobolTSI (unnormalized) for input %3d = %12.4e\n",
              ii+1, tsiMeds[ii]);
      for (ii = 0; ii < nInputs; ii++) means[ii] = (double) ii;
      sortDbleList2(nInputs, tsiMeds, means);
      for (ii = nInputs-1; ii >= 0; ii--)
         printOutTS(PL_INFO, 
              "RSMSobolTSI (unnormalized,ordered) for input %3d = %12.4e\n",
              (int) means[ii]+1,tsiMeds[ii]);
      printAsterisks(PL_INFO, 0);
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
   printOutTS(PL_ERROR,
              "RSMSobolTSI operator= ERROR: operation not allowed.\n");
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
   printEquals(PL_INFO, 0);
   if (variance == 0.0)
   {
      printOutTS(PL_INFO, "Total variance = 0. Hence, no total effect plot.\n");
      return 0;
   }

   printEquals(PL_INFO, 0);
   printOutTS(PL_INFO, "Total Effect Statistics: \n");
   totalSensitivity_ = new double[nInputs_];
   for (ii = 0; ii < nInputs; ii++)
   {
      printOutTS(PL_INFO, "Input %4d: Sobol' total sensitivity = %12.4e",
                 ii+1,tsi[ii]/variance);
      //save total sensitivity
      totalSensitivity_[ii] = tsi[ii]/variance;

      if (flag == 1)
           printOutTS(PL_INFO, ", bounds = [%12.4e, %12.4e]\n",
                      tsiMins[ii]/variance, tsiMaxs[ii]/variance);
      else printOutTS(PL_INFO, "\n");
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
      for (ii = 0; ii < nInputs; ii++) 
         fprintf(fp,"%24.16e\n",tsi[ii]/variance);
      fprintf(fp, "];\n");
      if (flag == 1)
      {
         fprintf(fp, "Mins = [\n");
         for (ii = 0; ii < nInputs; ii++)
            fprintf(fp,"%24.16e\n",tsiMins[ii]/variance);
         fprintf(fp, "];\n");
         fprintf(fp, "Maxs = [\n");
         for (ii = 0; ii < nInputs; ii++) 
            fprintf(fp,"%24.16e\n",tsiMaxs[ii]/variance);
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
         if (iNames[nInputs-1] != NULL) 
              fprintf(fp,"'%s'};\n",iNames[nInputs-1]);
         else fprintf(fp,"'X%d'};\n",nInputs);
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
           printOutTS(PL_INFO, 
                "RSMSobolTSI plot file = scilabrssoboltsi.sci\n");
      else printOutTS(PL_INFO, 
                "RSMSobolTSI plot file = matlabrssoboltsi.m\n");
      return 0;
   }
   else
   {
      printOutTS(PL_ERROR,
           "RSMSobolTSI ERROR: cannot create rssoboltsi graphics.\n");
      return 0;
   }
}

// ************************************************************************
// functions for getting results
// ------------------------------------------------------------------------
int RSMSobolTSIAnalyzer::get_nInputs()
{
   return nInputs_;
}
double RSMSobolTSIAnalyzer::get_outputMean()
{
   return outputMean_;
}
double RSMSobolTSIAnalyzer::get_outputStd()
{
   return outputStd_;
}
double *RSMSobolTSIAnalyzer::get_totalSensitivity()
{
  double* retVal = NULL;
  if (totalSensitivity_)
  {
     retVal = new double[nInputs_];
     for (int ii = 0; ii < nInputs_; ii++) 
        retVal[ii] = totalSensitivity_[ii];
  }
  return retVal;
}

