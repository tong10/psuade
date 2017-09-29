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
// Functions for the class RSMSobol2Analyzer  
// (Sobol' second order sensitivity analysis - with response surface)
// AUTHOR : CHARLES TONG
// DATE   : 2006
// ************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "PsuadeUtil.h"
#include "sysdef.h"
#include "Vector.h"
#include "Matrix.h"
#include "Psuade.h"
#include "FuncApprox.h"
#include "Sampling.h"
#include "RSConstraints.h"
#include "PDFManager.h"
#include "pData.h"
#include "PsuadeData.h"
#include "PsuadeConfig.h"
#include "RSMSobol2Analyzer.h"

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
RSMSobol2Analyzer::RSMSobol2Analyzer() : Analyzer()
{
   setName("RSMSOBOL2");
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
RSMSobol2Analyzer::~RSMSobol2Analyzer()
{
}

// ************************************************************************
// perform analysis
// ------------------------------------------------------------------------
double RSMSobol2Analyzer::analyze(aData &adata)
{
   int    nInputs, nOutputs, nSamples, ii, jj, kk, *S, status, outputID;
   int    length, nSubSamples=500, iL, sCnt, *SS, nLevels=200, noPDF=1;
   int    ii2, iL2, iR, currNLevels, nSamp, printLevel, *pdfFlags;
   int    *bins, totalCnt, *pdfFlags1, *pdfFlags2, selectedInput=-1;
   double *xLower, *xUpper, *X, *Y, *Y2, *YY, *cLower, *cUpper, *XX;
   double *oneSamplePt, *means, *vars, vce, ddata, variance, *ZZ;
   double dmean, ecv, *inputMeans, *inputStdevs, *mSamplePts;
   double *samplePts2D, *inputMeans1, *inputMeans2, *inputStdevs1;
   double *inputStdevs2;
   char   pString[500], *cString, winput[500], winput2[500];
   PsuadeData    *ioPtr;
   FuncApprox    *faPtr;
   RSConstraints *constrPtr;
   Vector        vecIn, vecOut, vecUB, vecLB;
   pData         pCorMat;
   Matrix        *corMatp, corMat;
   PDFManager    *pdfman, *pdfman1, *pdfman2;
   Sampling      *sampler;

   printAsterisks(0);
   printf("*          RS-based Second Order Sobol' Indices \n");
   printEquals(0);
   printf("* TO GAIN ACCESS TO DIFFERENT OPTIONS: SET\n");
   printf("*\n");
   printf("* - ana_expert mode to finetune RSMSobol2 parameters, \n");
   printf("*   (e.g. sample size for integration can be adjusted).\n");
   printf("* - rs_expert to mode finetune response surface for RSMSobol2,\n");
   printf("* - printlevel to 1 or higher to display more information.\n");
   printEquals(0);
   
   nInputs     = adata.nInputs_;
   nOutputs    = adata.nOutputs_;
   nSamples    = adata.nSamples_;
   xLower      = adata.iLowerB_;
   xUpper      = adata.iUpperB_;
   X           = adata.sampleInputs_;
   Y2          = adata.sampleOutputs_;
   S           = adata.sampleStates_;
   outputID    = adata.outputID_;
   ioPtr       = adata.ioPtr_;
   printLevel  = adata.printLevel_;
   pdfFlags    = adata.inputPDFs_;
   inputMeans  = adata.inputMeans_;
   inputStdevs = adata.inputStdevs_;
   if (pdfFlags != NULL)
   {
      for (ii = 0; ii < nInputs; ii++)
         if (pdfFlags[ii] != 0) noPDF = 0;
   }
   if (noPDF == 1) printf("* RSMSobol2 INFO: all uniform distributions.\n");
   else
   {
      printf("RSMSobol2 INFO: non-uniform distributions detected which\n");
      printf("                will be used in this analysis.\n");
   }

   if (nInputs <= 1 || nSamples <= 0 || nOutputs <= 0)
   {
      printf("RSMSobol2 ERROR: invalid arguments.\n");
      printf("   nInputs  = %d\n", nInputs);
      printf("   nOutputs = %d\n", nOutputs);
      printf("   nSamples = %d\n", nSamples);
      return PSUADE_UNDEFINED;
   }
   if (nInputs <= 2)
   {
      printf("RSMSobol2 ERROR: no need for this analysis (nInputs=2).\n"); 
      return PSUADE_UNDEFINED;
   }
   if (outputID < 0 || outputID >= nOutputs)
   {
      printf("RSMSobol2 ERROR: invalid outputID (%d).\n", outputID);
      return PSUADE_UNDEFINED;
   }
   if (ioPtr == NULL)
   {
      printf("RSMSobol2 ERROR: no data.\n");
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
            printf("RSMSobol2 INFO: this method cannot handle correlated\n");
            printf("          inputs using joint PDFs. PSUADE will try\n");
            printf("          a variant of this method, or you can re-run\n");
            printf("          using the group variance-based method.\n");
            return analyze2(adata);
         }
      }
   }
   status = 0;
   for (ii = 0; ii < nSamples; ii++)
      if (Y2[nOutputs*ii+outputID] > 0.9*PSUADE_UNDEFINED) status = 1;
   if (status == 1)
   {
      printf("RSMSobol2 ERROR: Some outputs are undefined. Prune the\n");
      printf("                 undefined sample points first.\n");
      return PSUADE_UNDEFINED;
   }
   Y = new double[nSamples];
   for (ii = 0; ii < nSamples; ii++) Y[ii] = Y2[ii*nOutputs+outputID];

   constrPtr = new RSConstraints();
   constrPtr->genConstraints(ioPtr);

   faPtr = genFAInteractive(ioPtr, 0);
   length = -999;
   status = faPtr->genNDGridData(X, Y, &length, NULL, NULL);

   printAsterisks(0);
   if (psAnaExpertMode_ == 1)
   {
      printf("* RSMSobol2 creates a sample of size M for each\n");
      printf("* input pair of K levels. The total sample size is \n");
      printf("* N = M * K * K * nInputs * (nInputs - 1) / 2.\n");
      printf("* NOW, nInputs = %d\n", nInputs);
      printf("* As a user, please decide on M and K.\n");
      printf("* Default M = %d\n", nSubSamples);
      printf("* Default K = %d\n", nLevels);
      printEquals(0);

      sprintf(pString,"Enter M (suggestion: 200 - 1000) : ");
      nSubSamples = getInt(200, 100000, pString);
      if (nSubSamples > 10000)
         printf("A M of %d may take very long time.\n",nSubSamples);

      sprintf(pString, "Enter K (suggestion: 200 - 500) : ");
      nLevels = getInt(200, 10000, pString);
      if (nLevels > 1000)
      {
         printf("A K of %d may take very long time.\n",nLevels);
      }
      sprintf(pString,
              "If analyze one input only, enter input number (0 if all): ");
      selectedInput = getInt(0, nInputs, pString);
      selectedInput--;
      printAsterisks(0);
   }
   else
   {
      nSubSamples = 1000;
      nLevels = 200;
      if (psConfig_ != NULL)
      {
         cString = psConfig_->getParameter("RSMSobol2_nsubsamples");
         if (cString != NULL)
         {
            sscanf(cString, "%s %s %d", winput, winput2, &nSubSamples);
            if (nSubSamples < 1000)
            {
               printf("RSMSobol2 INFO: nSubSamples should be >= 1000.\n");
               nSubSamples = 1000;
            }
            else
            {
               printf("RSMSobol2 INFO: nSubSamples = %d (config).\n",nSubSamples);
            }
         }
         cString = psConfig_->getParameter("RSMSobol2_nlevels");
         if (cString != NULL)
         {
            sscanf(cString, "%s %s %d", winput, winput2, &nLevels);
            if (nLevels < 200)
            {
               printf("RSMSobol2 INFO: nLevels should be >= 200.\n");
               nLevels = 200;
            }
            else
            {
               printf("RSMSobol2 INFO: nLevels = %d (config).\n",nLevels);
            }
         }
      }
      if (printLevel > 1)
      {
         printf("RSMSobol2: default M = %d.\n", nSubSamples);
         printf("RSMSobol2: default K = %d.\n", nLevels);
         printf("To change these settings, turn on ana_expert mode and rerun.\n");
      }
   }
   printEquals(0);

   nSamp = 100000;
   if (printLevel > 1)
   {
      printf("RSMSobol2 INFO: creating a sample for basic statistics.\n");
      printf("                sample size = %d\n", nSamp);
   }

   XX = new double[nSamp*nInputs];
   YY = new double[nSamp];

   if (noPDF == 0)
   {
      pdfman = new PDFManager();
      pdfman->initialize(nInputs,pdfFlags,inputMeans,
                         inputStdevs,*corMatp);
      vecLB.load(nInputs, xLower);
      vecUB.load(nInputs, xUpper);
      vecOut.setLength(nSamp*nInputs);
      pdfman->genSample(nSamp, vecOut, vecLB, vecUB);
      for (ii = 0; ii < nSamp*nInputs; ii++) XX[ii] = vecOut[ii];
      delete pdfman;
   }
   else
   {
      if (nInputs > 51)
           sampler = SamplingCreateFromID(PSUADE_SAMP_LHS);
      else sampler = SamplingCreateFromID(PSUADE_SAMP_LPTAU);
      sampler->setInputBounds(nInputs, xLower, xUpper);
      sampler->setOutputParams(1);
      sampler->setSamplingParams(nSamp, 1, 1);
      sampler->initialize(0);
      SS = new int[nSamp];
      sampler->getSamples(nSamp, nInputs, 1, XX, YY, SS);
      delete [] SS;
      delete sampler;
   }

   if (printLevel > 1)
   {
      printf("RSMSobol2: running the sample with response surface...\n");
   }
   faPtr->evaluatePoint(nSamp, XX, YY);
   if (printLevel > 1)
   {
      printf("RSMSobol2: done running the sample with response surface.\n");
   }
   
   for (ii = 0; ii < nSamp; ii++)
   {
      oneSamplePt = &(XX[ii*nInputs]);
      ddata = constrPtr->evaluate(oneSamplePt,YY[ii],status);
      if (status == 0) YY[ii] = PSUADE_UNDEFINED;
   }
   
   dmean = 0.0;
   sCnt = 0;
   for (ii = 0; ii < nSamp; ii++)
   {
      if (YY[ii] != PSUADE_UNDEFINED)
      {
         dmean += YY[ii];
         sCnt++;
      }
   }
   if (sCnt > 1) dmean /= (double) sCnt;
   else
   {
      printf("RSMSobol2 ERROR: too few samples that satisify\n");
      printf("constraints (%d out of %d).\n",sCnt,nSamp);
      delete [] XX;
      delete [] YY;
      delete faPtr;
      return PSUADE_UNDEFINED;
   }
   variance = 0.0;
   for (ii = 0; ii < nSamp; ii++)
   {
      if (YY[ii] != PSUADE_UNDEFINED)
         variance += (YY[ii] - dmean) * (YY[ii] - dmean) ;
   }
   variance /= (double) sCnt;
   printf("RSMSobol2: sample mean    (based on N = %d) = %10.3e\n",
          sCnt, dmean);
   printf("RSMSobol2: sample std dev (based on N = %d) = %10.3e\n",
          sCnt, sqrt(variance));
   if (variance == 0.0) variance = 1.0;
   delete [] XX;
   delete [] YY;

   cLower = new double[nInputs];
   cUpper = new double[nInputs];
   nSamp  = nSubSamples;
   XX     = new double[nSubSamples*nInputs];
   YY     = new double[nSubSamples];
   means  = new double[nLevels*nLevels];
   vars   = new double[nLevels*nLevels];
   bins   = new int[nLevels*nLevels];
   pdfFlags1    = new int[2];
   inputMeans1  = new double[2];
   inputStdevs1 = new double[2];
   pdfFlags2    = new int[nInputs-2];
   inputMeans2  = new double[nInputs-2];
   inputStdevs2 = new double[nInputs-2];
   samplePts2D  = new double[nLevels*2];
   mSamplePts   = new double[nInputs*nSubSamples];

   pData *pPtr = ioPtr->getAuxData();
   pPtr->nDbles_ = nInputs * nInputs;
   pPtr->dbleArray_ = new double[nInputs * nInputs];
   for (ii = 0; ii < nInputs*nInputs; ii++) pPtr->dbleArray_[ii] = 0.0;
   pPtr->dbleData_ = variance;

   printAsterisks(0);
   for (ii = 0; ii < nInputs; ii++)
   {
      for (ii2 = 0; ii2 < nInputs; ii2++)
      {
         vce = 0.0;
         if (ii2 <= ii) continue;
         if (selectedInput != -1 && ii != selectedInput && ii2 != selectedInput)
            continue;
         if (printLevel > 1)
           printf("RSMSobol2: processing input pair %d, %d\n",
                  ii+1, ii2+1);

         currNLevels = nLevels / 8;
         for (iR = 0; iR < 4; iR++)
         {
            if (printLevel > 3)
               printf("RSMSobol2: processing refinement %d\n",iR);

            cLower[0] = xLower[ii];
            cUpper[0] = xUpper[ii];
            cLower[1] = xLower[ii2];
            cUpper[1] = xUpper[ii2];
            if (noPDF == 0)
            {
               corMat.setDim(2,2);
               corMat.setEntry(0, 0, corMatp->getEntry(ii,ii));
               corMat.setEntry(1, 1, corMatp->getEntry(ii2,ii2));
               pdfFlags1[0] = pdfFlags[ii];
               pdfFlags1[1] = pdfFlags[ii2];
               inputMeans1[0] = inputMeans[ii];
               inputMeans1[1] = inputMeans[ii2];
               inputStdevs1[0] = inputStdevs[ii];
               inputStdevs1[1] = inputStdevs[ii2];
               pdfman1 = new PDFManager();
               pdfman1->initialize(2,pdfFlags1,inputMeans1,
                                   inputStdevs1,corMat);
               vecLB.load(2, cLower);
               vecUB.load(2, cUpper);
               vecOut.setLength(currNLevels*2);
               pdfman1->genSample(currNLevels, vecOut, vecLB, vecUB);
               for (jj = 0; jj < currNLevels*2; jj++) 
                  samplePts2D[jj] = vecOut[jj];
               delete pdfman1;
            }
            else
            {
               sampler = SamplingCreateFromID(PSUADE_SAMP_LHS);
               sampler->setInputBounds(2, cLower, cUpper);
               sampler->setOutputParams(1);
               sampler->setSamplingParams(currNLevels, 1, 0);
               sampler->initialize(0);
               SS = new int[currNLevels];
               ZZ = new double[currNLevels];
               sampler->getSamples(currNLevels,2,1,samplePts2D,ZZ,SS);
               delete [] SS;
               delete [] ZZ;
               delete sampler;
            }

            if (noPDF == 0)
            {
               corMat.setDim(nInputs-2, nInputs-2);
               for (jj = 0; jj < ii; jj++)
               {
                  cLower[jj] = xLower[jj];
                  cUpper[jj] = xUpper[jj];
                  corMat.setEntry(jj, jj, corMatp->getEntry(jj,jj));
                  pdfFlags2[jj] = pdfFlags[jj];
                  inputMeans2[jj] = inputMeans[jj];
                  inputStdevs2[jj] = inputStdevs[jj];
               }
               for (jj = ii+1; jj < ii2; jj++)
               {
                  cLower[jj-1] = xLower[jj];
                  cUpper[jj-1] = xUpper[jj];
                  corMat.setEntry(jj-1, jj-1, corMatp->getEntry(jj,jj));
                  pdfFlags2[jj-1] = pdfFlags[jj];
                  inputMeans2[jj-1] = inputMeans[jj];
                  inputStdevs2[jj-1] = inputStdevs[jj];
               }
               for (jj = ii2+1; jj < nInputs; jj++)
               {
                  cLower[jj-2] = xLower[jj];
                  cUpper[jj-2] = xUpper[jj];
                  corMat.setEntry(jj-2, jj-2, corMatp->getEntry(jj,jj));
                  pdfFlags2[jj-2] = pdfFlags[jj];
                  inputMeans2[jj-2] = inputMeans[jj];
                  inputStdevs2[jj-2] = inputStdevs[jj];
               }
               pdfman2 = new PDFManager();
               pdfman2->initialize(nInputs-2,pdfFlags2,inputMeans2,
                                   inputStdevs2,corMat);
               vecLB.load(nInputs-2, cLower);
               vecUB.load(nInputs-2, cUpper);
               vecOut.setLength(nSubSamples*(nInputs-2));
               pdfman2->genSample(nSubSamples, vecOut, vecLB, vecUB);
               for (jj = 0; jj < nSubSamples*(nInputs-2); jj++)
                  XX[jj] = vecOut[jj];
               delete pdfman2;
            }
            else
            {
               for (jj = 0; jj < ii; jj++)
               {
                  cLower[jj] = xLower[jj];
                  cUpper[jj] = xUpper[jj];
               }
               for (jj = ii+1; jj < ii2; jj++)
               {
                  cLower[jj-1] = xLower[jj];
                  cUpper[jj-1] = xUpper[jj];
               }
               for (jj = ii2+1; jj < nInputs; jj++)
               {
                  cLower[jj-2] = xLower[jj];
                  cUpper[jj-2] = xUpper[jj];
               }
               if (nInputs-1 > 51)
                    sampler = SamplingCreateFromID(PSUADE_SAMP_LHS);
               else sampler = SamplingCreateFromID(PSUADE_SAMP_LPTAU);
               sampler->setInputBounds(nInputs-2, cLower, cUpper);
               sampler->setOutputParams(1);
               sampler->setSamplingParams(nSubSamples, 1, 1);
               sampler->initialize(0);
               SS = new int[nSubSamples];
               ZZ = new double[nSubSamples];
               sampler->getSamples(nSubSamples,nInputs-2,1,XX,ZZ,SS);
               delete [] SS;
               delete [] ZZ;
               delete sampler;
            }

            for (iL = 0; iL < currNLevels; iL++)
            {
               for (iL2 = 0; iL2 < currNLevels; iL2++)
               {
                  for (jj = 0; jj < nSubSamples; jj++)
                  {
                     oneSamplePt = &(XX[jj*(nInputs-2)]);
                     for (kk = 0; kk < ii; kk++)
                        mSamplePts[jj*nInputs+kk] = oneSamplePt[kk];
                     for (kk = ii+1; kk < ii2; kk++)
                        mSamplePts[jj*nInputs+kk] = oneSamplePt[kk-1];
                     for (kk = ii2+1; kk < nInputs; kk++)
                        mSamplePts[jj*nInputs+kk] = oneSamplePt[kk-2];
                     mSamplePts[jj*nInputs+ii] = samplePts2D[iL*2];
                     mSamplePts[jj*nInputs+ii2] = samplePts2D[iL2*2+1];
                  }

                  faPtr->evaluatePoint(nSubSamples,mSamplePts,YY);

                  for (jj = 0; jj < nSubSamples; jj++)
                  {
                     oneSamplePt = &(mSamplePts[jj*nInputs]);
                     ddata = constrPtr->evaluate(oneSamplePt,YY[jj],status);
                     if (status == 0) YY[jj] = PSUADE_UNDEFINED;
                  }

                  means[iL*currNLevels+iL2] = 0.0;
                  sCnt = 0;
                  for (jj = 0; jj < nSubSamples; jj++)
                  {
                     if (YY[jj] != PSUADE_UNDEFINED)
                     {
                        means[iL*currNLevels+iL2] += YY[jj];
                        sCnt++;
                     }
                  }
                  bins[iL*currNLevels+iL2] = sCnt;
                  if (sCnt < 1 && printLevel >= 5)
                     printf("RSMSobol2 WARNING: subsample size = 0.\n");
                  if (sCnt < 1) means[iL*currNLevels+iL2] = PSUADE_UNDEFINED;
                  else          means[iL*currNLevels+iL2] /= (double) sCnt;

                  vars[iL*currNLevels+iL2] = 0.0;
                  ddata = means[iL*currNLevels+iL2];
                  for (jj = 0; jj < nSubSamples; jj++)
                  {
                     if (YY[jj] != PSUADE_UNDEFINED)
                        vars[iL*currNLevels+iL2] += 
                                      (YY[jj]-ddata)*(YY[jj]-ddata);
                  }
                  if (sCnt < 1) vars[iL*currNLevels+iL2] = PSUADE_UNDEFINED;
                  else          vars[iL*currNLevels+iL2] /= (double) sCnt;

                  if (printLevel > 4)
                  {
                     printf("RSMSobol2: inputs (%d,%d)\n",ii+1,ii2+1);
                     printf("  refinement %d, levels (%d,%d), size %d (%d)\n",
                            iR, iL+1, iL2+1, sCnt, nSubSamples);
                     printf("     mean = %e\n", means[iL*currNLevels+iL2]);
                     printf("     var  = %e\n", vars[iL*currNLevels+iL2]);
                  }
               }
            }

            dmean = 0.0;
            totalCnt = 0;
            for (iL = 0; iL < currNLevels; iL++)
               for (iL2 = 0; iL2 < currNLevels; iL2++)
                  totalCnt += bins[iL*currNLevels+iL2];;
            if (totalCnt == 0)
            {
               printf("RSMSobol2 ERROR: empty constrained space.\n");
               printf("          Either try larger sample size or\n");
               printf("          use looser constraints.\n");
               exit(1);
            }

            for (iL = 0; iL < currNLevels; iL++)
            {
               for (iL2 = 0; iL2 < currNLevels; iL2++)
               {
                  if (means[iL*currNLevels+iL2] != PSUADE_UNDEFINED)
                     dmean += means[iL*currNLevels+iL2] * 
                              bins[iL*currNLevels+iL2] / totalCnt;
               }
            }

            vce = 0.0;
            for (iL = 0; iL < currNLevels; iL++)
               for (iL2 = 0; iL2 < currNLevels; iL2++)
                  if (means[iL*currNLevels+iL2] != PSUADE_UNDEFINED)
                     vce += (means[iL*currNLevels+iL2] - dmean) * 
                            (means[iL*currNLevels+iL2] - dmean) *
                            bins[iL*currNLevels+iL2] / totalCnt;

            ecv = 0.0;
            for (iL = 0; iL < currNLevels; iL++)
            {
               for (iL2 = 0; iL2 < currNLevels; iL2++)
               {
                  if (vars[iL*currNLevels+iL2] != PSUADE_UNDEFINED)
                     ecv += vars[iL*currNLevels+iL2] *
                            bins[iL*currNLevels+iL2] / totalCnt;
               }
            }

            if (printLevel > 2 || iR == 3)
            {
               printf("VCE(%3d,%3d) = %12.4e, (normalized) = %10.3e\n",
                      ii+1, ii2+1, vce, vce/variance);
            }
            if (printLevel > 3)
            {
               printf("ECV(%3d,%3d) = %12.4e, (normalized) = %10.3e\n",
                      ii+1, ii2+1, ecv, ecv/variance);
            }
            currNLevels *= 2;
         }
         pPtr->dbleArray_[ii*nInputs+ii2] = vce;
      }
   }
   printAsterisks(0);
    
   delete constrPtr;
   delete faPtr;
   delete [] cLower;
   delete [] cUpper;
   delete [] XX;
   delete [] YY;
   delete [] Y;
   delete [] means;
   delete [] vars;
   delete [] bins;
   delete [] mSamplePts;
   delete [] samplePts2D;
   delete [] pdfFlags1;
   delete [] pdfFlags2;
   delete [] inputMeans1;
   delete [] inputMeans2;
   delete [] inputStdevs1;
   delete [] inputStdevs2;
   return 0.0;
}

// ************************************************************************
// perform analysis (for problems with joint PDFs)
// ------------------------------------------------------------------------
double RSMSobol2Analyzer::analyze2(aData &adata)
{
   int    nInputs, nOutputs, nSamples, ii, ii2, jj, iR, *S, status;
   int    length, nSubSamples=500, iL, iL2, sCnt, nLevels=200, outputID;
   int    currNLevels, nSamp, printLevel, *bins, totalCnt, bin1, bin2;
   double *xLower, *xUpper, *X, *Y, *XX, *YY, *Y2, dmean, ecv, ddata;
   double *oneSamplePt, *means, *vars, vce, variance, width1, width2;
   char   pString[500];
   PsuadeData    *ioPtr;
   FuncApprox    *faPtr;
   RSConstraints *constrPtr;
   Vector        vecIn, vecOut, vecUB, vecLB;
   PDFManager    *pdfman;

   printAsterisks(0);
   printf("RSMSobol2: since joint PDFs have been specified, a different\n");
   printf("           interaction analysis will be performed.\n");
   nInputs     = adata.nInputs_;
   nOutputs    = adata.nOutputs_;
   nSamples    = adata.nSamples_;
   xLower      = adata.iLowerB_;
   xUpper      = adata.iUpperB_;
   X           = adata.sampleInputs_;
   Y2          = adata.sampleOutputs_;
   S           = adata.sampleStates_;
   outputID    = adata.outputID_;
   ioPtr       = adata.ioPtr_;
   printLevel  = adata.printLevel_;
   Y = new double[nSamples];
   for (ii = 0; ii < nSamples; ii++) Y[ii] = Y2[ii*nOutputs+outputID];

   constrPtr = new RSConstraints();
   constrPtr->genConstraints(ioPtr);

   faPtr = genFAInteractive(ioPtr, 0);
   if (faPtr == NULL)
   {
      printf("RSMSobol2 ERROR: cannot create function approximator.\n");
      delete [] Y;
      delete constrPtr;
      return 1.0e12;
   }
   length = -999;
   status = faPtr->genNDGridData(X, Y, &length, NULL, NULL);

   printAsterisks(0);
   if (psAnaExpertMode_ == 1)
   {
      printf("* RSMSobol2 creates a sample of size M for each\n");
      printf("* input pair of K levels each. The total sample size is \n");
      printf("* N = M * K * K.\n");
      printf("* NOW, nInputs = %d\n", nInputs);
      printf("* As a user, please decide on M and K.\n");
      printEquals(0);

      sprintf(pString,"Enter M (suggestion: 200 - 1000) : ");
      nSubSamples = getInt(200, 100000, pString);
      if (nSubSamples > 10000)
         printf("A M of %d may take very long time.\n",nSubSamples);

      sprintf(pString, "Enter nLevels (suggestion: 200 - 500) : ");
      nLevels = getInt(200, 10000, pString);
      if (nLevels > 1000)
      {
         printf("* A K of %d may take very long time.\n",nLevels);
      }
      printAsterisks(0);
   }
   else
   {
      nSubSamples = 500;
      nLevels = 200;
      printf("* RSMSobol2: default M = %d.\n", nSubSamples);
      printf("* RSMSobol2: default K = %d.\n", nLevels);
      printf("* To change these settings, turn on ana_expert mode and rerun.\n");
   }
   printEquals(0);

   nSamp = nSubSamples * nLevels * nLevels;
   if (printLevel > 1)
   {
      printf("* RSMSobol2 INFO: creating a sample for basic statistics.\n");
      printf("*                 sample size = %d\n", nSamp);
   }

   XX  = new double[nSamp*nInputs];
   YY = new double[nSamp];

   pdfman = new PDFManager();
   pdfman->initialize(ioPtr);
   vecLB.load(nInputs, xLower);
   vecUB.load(nInputs, xUpper);
   vecOut.setLength(nSamp*nInputs);
   pdfman->genSample(nSamp, vecOut, vecLB, vecUB);
   for (ii = 0; ii < nSamp*nInputs; ii++) XX[ii] = vecOut[ii];

   if (printLevel > 1)
   {
      printf("RSMSobol2: running the sample with response surface...\n");
   }
   faPtr->evaluatePoint(nSamp, XX, YY);
   if (printLevel > 1)
   {
      printf("RSMSobol2: done running the sample with response surface.\n");
   }
   
   for (ii = 0; ii < nSamp; ii++)
   {
      oneSamplePt = &(XX[ii*nInputs]);
      ddata = constrPtr->evaluate(oneSamplePt,YY[ii],status);
      if (status == 0) YY[ii] = PSUADE_UNDEFINED;
   }
   
   dmean = 0.0;
   sCnt = 0;
   for (ii = 0; ii < nSamp; ii++)
   {
      if (YY[ii] != PSUADE_UNDEFINED)
      {
         dmean += YY[ii];
         sCnt++;
      }
   }
   if (sCnt > 1) dmean /= (double) sCnt;
   else
   {
      printf("RSMSobol2 ERROR: too few samples that satisify\n");
      printf("constraints (%d out of %d).\n",sCnt,nSamp);
      delete [] XX;
      delete [] YY;
      delete faPtr;
      delete pdfman;
      delete constrPtr;
      return PSUADE_UNDEFINED;
   }
   variance = 0.0;
   for (ii = 0; ii < nSamp; ii++)
   {
      if (YY[ii] != PSUADE_UNDEFINED)
         variance += (YY[ii] - dmean) * (YY[ii] - dmean) ;
   }
   variance /= (double) sCnt;
   printf("* RSMSobol2: sample mean    (based on %d points) = %e\n",
          sCnt, dmean);
   printf("* RSMSobol2: total std dev  (based on %d points) = %e\n",
          sCnt, sqrt(variance));
   if (variance == 0.0) variance = 1.0;
   delete pdfman;
   printAsterisks(0);

   means = new double[nLevels*nLevels];
   vars  = new double[nLevels*nLevels];
   bins  = new int[nLevels*nLevels];

   for (ii = 0; ii < nInputs; ii++)
   {
      for (ii2 = ii+1; ii2 < nInputs; ii2++)
      {
         if (printLevel > 1)
           printf("RSMSobol2: processing input pair %d, %d\n",
                  ii+1, ii2+1);

         currNLevels = nLevels / 8;
         for (iR = 0; iR < 4; iR++)
         {
            if (printLevel > 3)
               printf("RSMSobol2: processing refinement %d\n",iR);

            width1 = (xUpper[ii] - xLower[ii]) / currNLevels;
            width2 = (xUpper[ii2] - xLower[ii2]) / currNLevels;
            for (iL = 0; iL < currNLevels*currNLevels; iL++)
            {
               means[iL] = 0.0;
               vars[iL] = 0.0;
               bins[iL] = 0;
            }
            for (jj = 0; jj < nSamp; jj++)
            {
               if (YY[jj] != PSUADE_UNDEFINED)
               {
                  ddata = XX[jj*nInputs+ii];
                  bin1 = (int) ((ddata - xLower[ii]) / width1);
                  ddata = XX[jj*nInputs+ii2];
                  bin2 = (int) ((ddata - xLower[ii2]) / width2);
                  if (bin1 == currNLevels) bin1 = currNLevels - 1;
                  if (bin2 == currNLevels) bin2 = currNLevels - 1;
                  means[bin1*currNLevels+bin2] += YY[jj];
                  bins[bin1*currNLevels+bin2]++;
               }
            }
            for (iL = 0; iL < currNLevels; iL++)
            {
               for (iL2 = 0; iL2 < currNLevels; iL2++)
               {
                  sCnt = bins[iL*currNLevels+iL2];
                  if (sCnt < 1 && printLevel >= 5)
                     printf("RSMSobol2 WARNING: subsample size = 0.\n");
                  if (sCnt < 1) means[iL*currNLevels+iL2] = PSUADE_UNDEFINED;
                  else          means[iL*currNLevels+iL2] /= (double) sCnt;
               }
            }
            for (jj = 0; jj < nSamp; jj++)
            {
               if (YY[jj] != PSUADE_UNDEFINED)
               {
                  ddata = XX[jj*nInputs+ii];
                  bin1 = (int) ((ddata - xLower[ii]) / width1);
                  ddata = XX[jj*nInputs+ii2];
                  bin2 = (int) ((ddata - xLower[ii2]) / width2);
                  if (bin1 == currNLevels) bin1 = currNLevels - 1;
                  if (bin2 == currNLevels) bin2 = currNLevels - 1;
                  ddata = means[bin1*currNLevels+bin2];
                  vars[bin1*currNLevels+bin2] += 
                         (YY[jj]-ddata)*(YY[jj]-ddata);
               }
            }
            for (iL = 0; iL < currNLevels; iL++)
            {
               for (iL2 = 0; iL2 < currNLevels; iL2++)
               {
                  sCnt = bins[iL*currNLevels+iL2];
                  if (sCnt < 1) vars[iL*currNLevels+iL2] = PSUADE_UNDEFINED;
                  else          vars[iL*currNLevels+iL2] /= (double) sCnt;
                  if (printLevel > 3)
                  {
                     printf("RSMSobol2: inputs (%d,%d)\n",ii+1,ii2+1);
                     printf("  refinement %d, levels (%d,%d), size %d (%d)\n",
                            iR, iL+1, iL2+1, sCnt, nSubSamples);
                     printf("     mean = %e\n", means[iL*currNLevels+iL2]);
                     printf("     var  = %e\n", vars[iL*currNLevels+iL2]);
                  }
               }
            }

            totalCnt = 0;
            for (iL = 0; iL < currNLevels; iL++)
               for (iL2 = 0; iL2 < currNLevels; iL2++)
                  totalCnt += bins[iL*currNLevels+iL2];;
            if (totalCnt == 0)
            {
               printf("RSMSobol2 ERROR: empty constrained space.\n");
               printf("          Either try larger sample size or\n");
               printf("          use looser constraints.\n");
               exit(1);
            }

            dmean = 0.0;
            for (iL = 0; iL < currNLevels; iL++)
            {
               for (iL2 = 0; iL2 < currNLevels; iL2++)
               {
                  if (means[iL*currNLevels+iL2] != PSUADE_UNDEFINED)
                     dmean += means[iL*currNLevels+iL2] * 
                              bins[iL*currNLevels+iL2] / totalCnt;
               }
            }
            vce = 0.0;
            for (iL = 0; iL < currNLevels; iL++)
               for (iL2 = 0; iL2 < currNLevels; iL2++)
                  if (means[iL*currNLevels+iL2] != PSUADE_UNDEFINED)
                     vce += (means[iL*currNLevels+iL2] - dmean) * 
                            (means[iL*currNLevels+iL2] - dmean) *
                            bins[iL*currNLevels+iL2] / totalCnt;

            ecv = 0.0;
            for (iL = 0; iL < currNLevels; iL++)
            {
               for (iL2 = 0; iL2 < currNLevels; iL2++)
               {
                  if (vars[iL*currNLevels+iL2] != PSUADE_UNDEFINED)
                     ecv += vars[iL*currNLevels+iL2] *
                            bins[iL*currNLevels+iL2] / totalCnt;
               }
            }
            if (printLevel > 2 || iR == 3)
            {
               printf("VCE(%3d,%3d) = %12.4e, (normalized) = %10.3e\n",
                      ii+1, ii2+1, vce, vce/variance);
            }
            if (printLevel > 3 || iR == 3)
            {
               printf("ECV(%3d,%3d) = %12.4e, (normalized) = %10.3e\n",
                      ii+1, ii2+1, ecv, ecv/variance);
            }
            currNLevels *= 2;
         }
      }
   }
   printAsterisks(0);
    
   delete constrPtr;
   delete faPtr;
   delete [] XX;
   delete [] YY;
   delete [] Y;
   delete [] means;
   delete [] vars;
   delete [] bins;
   return 0.0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
RSMSobol2Analyzer& RSMSobol2Analyzer::operator=(const RSMSobol2Analyzer &)
{
   printf("RSMSobol2 operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

