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
// Functions for the class RSMSobol1Analyzer  
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
#include "pData.h"
#include "RSMSobol1Analyzer.h"
#include "Sampling.h"
#include "PDFManager.h"
#include "PDFNormal.h"
#include "RSConstraints.h"
#include "Psuade.h"
#include "PsuadeData.h"
#include "PsuadeConfig.h"

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
RSMSobol1Analyzer::RSMSobol1Analyzer() : Analyzer()
{
   setName("RSMSOBOL1");
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
RSMSobol1Analyzer::~RSMSobol1Analyzer()
{
}

// ************************************************************************
// perform analysis
// ------------------------------------------------------------------------
double RSMSobol1Analyzer::analyze(aData &adata)
{
   int    nInputs, nOutputs, nSamples, *S, *SS, nSteps, outputID, noPDF=1;
   int    ii, jj, kk, iL, iR, sCnt, length, status, currNLevels, method=1;
   int    nSubSamples=2000, nLevels=500, pdfFileFlag=0, *bins, *pdfFlags2;
   int    printLevel, nSamp, scatterFileFlag=0, totalCnt, *pdfFlags;
   double *X, *Y, *Y2, *XX, *YY, *ZZ, *oneSamplePt, dmean;
   double *xLower, *xUpper, *cLower, *cUpper, *means, variance, *vces;
   double *vars, *ecvs, ddata, *mSamplePts, *inputMeans;
   double *inputStdevs, *inputMeans2, *inputStdevs2, *samplePts1D;
   char   winput1[500], winput2[500], pdfFile[500], scatterFile[500];
   char   *cString, pString[500];
   FILE   *fp;
   PsuadeData    *ioPtr;
   FuncApprox    *faPtr;
   RSConstraints *constrPtr;
   Vector        vecIn, vecOut, vecUB, vecLB;
   pData         pCorMat;
   Matrix        *corMatp, corMat;
   PDFManager    *pdfman, *pdfman1, *pdfman2;
   Sampling      *sampler;

   if (method == 1)
   {
      analyze3(adata);
      return 0;
   }
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
   printAsterisks(0);
   printf("* RSMSobol1 constructor\n");
   printDashes(0);
   if (noPDF == 1) printf("* RSMSobol1 INFO: all uniform distributions.\n");
   else
   {
      printf("* RSMSobol1 INFO: non-uniform distributions detected\n");
      printf("                  which will be used in this analysis.\n");
   }

   if (nInputs <= 0 || nSamples <= 0 || nOutputs <= 0) 
   {
      printf("RSMSobol1 ERROR: invalid arguments.\n");
      printf("   nInputs  = %d\n", nInputs);
      printf("   nOutputs = %d\n", nOutputs);
      printf("   nSamples = %d\n", nSamples);
      return PSUADE_UNDEFINED;
   } 
   if (outputID >= nOutputs || outputID < 0)
   {
      printf("RSMSobol1 ERROR: invalid output ID (%d).\n", outputID);
      return PSUADE_UNDEFINED;
   }
   if (nInputs <= 1)
   {
      printf("RSMSobol1 ERROR: nInputs=1 does not need this analysis.\n");
      return PSUADE_UNDEFINED;
   }
   if (ioPtr == NULL)
   {
      printf("RSMSobol1 ERROR: no data (PsuadeData).\n");
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
            printf("RSMSobol1 INFO: this method cannot handle correlated\n");
            printf("          inputs using joint PDFs. PSUADE will use\n");
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
      printf("RSMSobol1 ERROR: Some outputs are undefined. Prune the\n");
      printf("                 undefined sample point first.\n");
      return PSUADE_UNDEFINED;
   }
   Y = new double[nSamples];
   for (ii = 0; ii < nSamples; ii++) Y[ii] = Y2[ii*nOutputs+outputID];

   constrPtr = new RSConstraints();
   constrPtr->genConstraints(ioPtr);

   faPtr = genFAInteractive(ioPtr, 0);
   length = -999;
   status = faPtr->genNDGridData(X, Y, &length, NULL, NULL);

   if (psAnaExpertMode_ == 1)
   {
      printAsterisks(0);
      printf("* RSMSobol1 creates one sample of size M \n");
      printf("* for each of the K input levels. Therefore,\n");
      printf("* the total sample size is\n");
      printf("* N = M * K * nInputs\n");
      nSubSamples = 1000;
      nLevels = 200;
      printf("* nInputs m = %d, and\n", nInputs);
      printf("* default M = %d\n", nSubSamples);
      printf("* default K = %d\n", nLevels);
      printf("* As a user, please decide on M and K.\n");
      printf("* Note: large M and K may take a long time\n");
      printEquals(0);
      sprintf(pString,"Enter M (suggestion: 1000-10000) : ");
      nSubSamples = getInt(1000, 50000, pString);
      sprintf(pString,"Enter K (suggestion: 100 - 1000) : ");
      nLevels = getInt(100, 5000, pString);
      printAsterisks(0);
   }
   else
   {
      nSubSamples = 1000;
      nLevels = 200;
      if (psConfig_ != NULL)
      {
         cString = psConfig_->getParameter("RSMSobol1_nsubsamples");
         if (cString != NULL)
         {
            sscanf(cString, "%s %s %d", winput1, winput2, &nSubSamples);
            if (nSubSamples < 1000)
            {
               printf("RSMSobol1 INFO: nSubSamples should be >= 1000.\n");
               nSubSamples = 1000;
            }
            else
            {
               printf("RSMSobol1 INFO: nSubSamples = %d (config).\n",nSubSamples);
            }
         }
         cString = psConfig_->getParameter("RSMSobol1_nlevels");
         if (cString != NULL)
         {
            sscanf(cString, "%s %s %d", winput1, winput2, &nLevels);
            if (nLevels < 200)
            {
               printf("RSMSobol1 INFO: nLevels should be >= 200.\n");
               nLevels = 200;
            }
            else
            {
               printf("RSMSobol1 INFO: nLevels = %d (config).\n",nLevels);
            }
         }
      }
      if (printLevel > 1)
      {
         printAsterisks(0);
         printf("* RSMSobol1 creates one sample of size M \n");
         printf("* for each of the K input levels. Therefore,\n");
         printf("* the total sample size is\n");
         printf("* N = M * K * nInputs.\n");
         printf("* nInputs m = %d \n", nInputs);
         printf("* default M = %d\n", nSubSamples);
         printf("* default K = %d\n", nLevels);
         printf("* To change these settings, turn on ana_expert mode and rerun.\n");
         printAsterisks(0);
      }
   }

   if (psAnaExpertMode_ == 0 && psConfig_ != NULL)
   {
      cString = psConfig_->getParameter("RSMSobol1_pdffile");
      if (cString != NULL)
      {
         strcpy(pdfFile, cString);
         pdfFileFlag = 1; 
      }
   }
   if (psAnaExpertMode_ == 1)
   {
      printf("RSMSobol1 will create a sample for basic statistics. You\n");
      printf("have the option to plot the probability density function.\n");
      sprintf(pString, "Create a pdf (probability) bar graph? (y or n) ");
      getString(pString, winput1);
      if (winput1[0] == 'y')
      {
         sprintf(pString,"Enter the file name for the bar graph : ");
         getString(pString, pdfFile);
         pdfFile[strlen(pdfFile)-1] = '\0';
         fp = fopen(pdfFile, "w");
         if (fp != NULL)
         {
            fclose(fp);
            pdfFileFlag = 1; 
         }
         else 
         {
            printf("RSMSobol1 WARNING: cannot open file %s\n", pdfFile);
            pdfFileFlag = 0; 
         }
      }
      printEquals(0);
   }
   else pdfFileFlag = 0;

   if (psAnaExpertMode_ == 0 && psConfig_ != NULL)
   {
      cString = psConfig_->getParameter("RSMSobol1_scatterfile");
      if (cString != NULL)
      {
         strcpy(pdfFile, cString);
         scatterFileFlag = 1; 
      }
   }
   if (psAnaExpertMode_ == 1)
   {
      printf("RSMSobol1 will create many samples for Sobol1 analysis. You\n");
      printf("have the option to plot these sample data.\n");
      sprintf(pString, "Create a scatter plot for RSMSobol1? (y or n) ");
      getString(pString, winput1);
      if (winput1[0] == 'y')
      {
         sprintf(pString,"Enter the file name for the scatter plot : ");
         getString(pString, scatterFile);
         scatterFile[strlen(scatterFile)-1] = '\0';
         fp = fopen(scatterFile, "w");
         if (fp != NULL)
         {
            fclose(fp);
            scatterFileFlag = 1;
         }
         else 
         {
            printf("RSMSobol1 ERROR: cannot open file %s\n", 
                   scatterFile); 
            scatterFileFlag = 0;
         }
      }
      printEquals(0);
   }
   else scatterFileFlag = 0;

   nSamp = 100000;
   if (printLevel > 1)
   {
      printf("RSMSobol1 INFO: creating a sample for basic statistics.\n");
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
      else
         sampler = SamplingCreateFromID(PSUADE_SAMP_LPTAU);
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
      printf("RSMSobol1: running the sample with response surface...\n");
   }
   faPtr->evaluatePoint(nSamp, XX, YY);
   if (printLevel > 1)
   {
      printf("RSMSobol1: done running the sample with response surface.\n");
   }

   for (ii = 0; ii < nSamp; ii++)
   {
      oneSamplePt = &(XX[ii*nInputs]);
      ddata = constrPtr->evaluate(oneSamplePt,YY[ii],status);
      if (status == 0) YY[ii] = PSUADE_UNDEFINED;
   }

   sCnt  = 0;
   dmean = 0.0;
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
      printf("RSMSobol1 ERROR: too few samples that satisify the ");
      printf("constraints (%d out of %d).\n", sCnt, nSamp);
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
   printf("RSMSobol1: sample mean    (based on N = %d) = %10.3e\n",
          sCnt, dmean);
   printf("RSMSobol1: sample std dev (based on N = %d) = %10.3e\n",
          sCnt, sqrt(variance));
   if (variance == 0.0) variance = 1.0;

   nSteps = 1;
   if (nSamp > 2000) nSteps = nSamp / 1000;
   if (pdfFileFlag == 1)
   {
      fp = fopen(pdfFile, "w");
      if (fp != NULL)
      {
         if (psPlotTool_ == 1)
              fprintf(fp, "// Sample PDF based on response surface.\n");
         else fprintf(fp, "%% Sample PDF based on response surface.\n");
         fprintf(fp, "A = [\n");
         for (ii = 0; ii < nSamp; ii++)
         {
            if (YY[ii] != PSUADE_UNDEFINED)
            {
               if ((pdfFileFlag == 1) && (ii % nSteps) == 0) 
                  fprintf(fp,"   %e\n", YY[ii]);
            }
         }
         fprintf(fp, "];\n");
         fprintf(fp, "amax = max(A);\n");
         fprintf(fp, "amin = min(A);\n");
         fprintf(fp, "step = (amax - amin) / 10;\n");
         fprintf(fp, "numa = size(A,1);\n");
         fprintf(fp, "B = zeros(numa,1);\n");
         fprintf(fp, "for ii = 1 : numa\n");
         fprintf(fp, "   B(ii) = ceil(((A(ii) - amin))/step);\n");
         fprintf(fp, "end;\n");
         fprintf(fp, "X = zeros(10,1);\n");
         fprintf(fp, "Y = zeros(10,1);\n");
         fprintf(fp, "for ii = 1 : 10\n");
         fprintf(fp, "   X(ii) = amin + (ii - 1 + 0.5) * step;\n");
         fprintf(fp, "   [ia,ja,aa] = find(B==ii);\n");
         fprintf(fp, "   Y(ii) = length(ia);\n");
         fprintf(fp, "end;\n");
         fprintf(fp, "Y = Y / numa;\n");
         fprintf(fp, "bar(X,Y,1.0)\n");
         fwritePlotXLabel(fp, "Output Values");
         fwritePlotYLabel(fp, "Probabilities");
         fwritePlotAxes(fp);
         fclose(fp);
      }
   }
   delete [] XX;
   delete [] YY;

   cLower = new double[nInputs];
   cUpper = new double[nInputs];
   nSamp  = nSubSamples;
   XX            = new double[nSubSamples*nInputs];
   samplePts1D   = new double[nLevels];
   mSamplePts    = new double[nInputs*nSubSamples];
   YY            = new double[nSubSamples];
   means         = new double[nLevels];
   vars          = new double[nLevels];
   vces          = new double[nInputs];
   ecvs          = new double[nInputs];
   bins          = new int[nLevels];
   pdfFlags2     = new int[nInputs-1];
   inputMeans2   = new double[nInputs-1];
   inputStdevs2  = new double[nInputs-1];

   fp = NULL;
   if (scatterFileFlag == 1)
   {
      fp = fopen(scatterFile, "w");
      if (psPlotTool_ == 1)
           fprintf(fp, "// scatter plot for RSMSobol1 sample\n");
      else fprintf(fp, "%% scatter plot for RSMSobol1 sample\n");
      fprintf(fp, "A = [ \n");
   }

   for (ii = 0; ii < nInputs; ii++)
   {
      if (printLevel >= 1)
         printf("RSMSobol1: processing input %d\n", ii+1);

      currNLevels = nLevels / 8;
      for (iR = 0; iR < 4; iR++)
      {
         if (noPDF == 0)
         {
            corMat.setDim(1,1);
            corMat.setEntry(0, 0, corMatp->getEntry(ii,ii));
            pdfman1 = new PDFManager();
            pdfman1->initialize(1,&pdfFlags[ii],&inputMeans[ii],
                                &inputStdevs[ii],corMat);
            vecLB.load(1, &xLower[ii]);
            vecUB.load(1, &xUpper[ii]);
            vecOut.setLength(currNLevels);
            pdfman1->genSample(currNLevels, vecOut, vecLB, vecUB);
            for (jj = 0; jj < currNLevels; jj++)
               samplePts1D[jj] = vecOut[jj];
            delete pdfman1;
         }
         else
         {
            sampler = SamplingCreateFromID(PSUADE_SAMP_LHS);
            sampler->setInputBounds(1, &xLower[ii], &xUpper[ii]);
            sampler->setOutputParams(1);
            sampler->setSamplingParams(currNLevels, 1, 1);
            sampler->initialize(0);
            SS = new int[currNLevels];
            ZZ = new double[currNLevels];
            sampler->getSamples(currNLevels, 1, 1, samplePts1D, ZZ, SS);
            delete [] SS;
            delete [] ZZ;
            delete sampler;
         }

         if (noPDF == 0)
         {
            corMat.setDim(nInputs-1, nInputs-1);
            for (jj = 0; jj < ii; jj++)
            {
               cLower[jj] = xLower[jj];
               cUpper[jj] = xUpper[jj];
               pdfFlags2[jj] = pdfFlags[jj];
               inputMeans2[jj] = inputMeans[jj];
               inputStdevs2[jj] = inputStdevs[jj];
               for (kk = 0; kk < ii; kk++)
                  corMat.setEntry(jj, kk, corMatp->getEntry(jj,kk));
               for (kk = ii+1; kk < nInputs; kk++)
                  corMat.setEntry(jj, kk-1, corMatp->getEntry(jj,kk));
            }
            for (jj = ii+1; jj < nInputs; jj++)
            {
               cLower[jj-1] = xLower[jj];
               cUpper[jj-1] = xUpper[jj];
               pdfFlags2[jj-1] = pdfFlags[jj];
               inputMeans2[jj-1] = inputMeans[jj];
               inputStdevs2[jj-1] = inputStdevs[jj];
               for (kk = 0; kk < ii; kk++)
                  corMat.setEntry(jj-1, kk, corMatp->getEntry(jj,kk));
               for (kk = ii+1; kk < nInputs; kk++)
                  corMat.setEntry(jj-1, kk-1, corMatp->getEntry(jj,kk));
            }
            pdfman2 = new PDFManager();
            pdfman2->initialize(nInputs-1,pdfFlags2,inputMeans2,
                                inputStdevs2,corMat);
            vecLB.load(nInputs-1, cLower);
            vecUB.load(nInputs-1, cUpper);
            vecOut.setLength(nSubSamples*(nInputs-1));
            pdfman2->genSample(nSubSamples, vecOut, vecLB, vecUB);
            for (jj = 0; jj < nSubSamples*(nInputs-1); jj++) XX[jj] = vecOut[jj];
            delete pdfman2;
         }
         else
         {
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
            if (nInputs-1 > 51)
                 sampler = SamplingCreateFromID(PSUADE_SAMP_LHS);
            else sampler = SamplingCreateFromID(PSUADE_SAMP_LPTAU);
            sampler->setInputBounds(nInputs-1, cLower, cUpper);
            sampler->setOutputParams(1);
            sampler->setSamplingParams(nSubSamples, 1, 1);
            sampler->initialize(0);
            SS = new int[nSubSamples];
            ZZ = new double[nSubSamples];
            sampler->getSamples(nSubSamples, nInputs-1, 1, XX, ZZ, SS);
            delete [] SS;
            delete [] ZZ;
            delete sampler;
         }

         for (iL = 0; iL < currNLevels; iL++)
         {
            for (jj = 0; jj < nSubSamples; jj++)
            {
               oneSamplePt = &(XX[jj*(nInputs-1)]);
               for (kk = 0; kk < ii; kk++)
                  mSamplePts[jj*nInputs+kk] = oneSamplePt[kk];
               for (kk = ii+1; kk < nInputs; kk++)
                  mSamplePts[jj*nInputs+kk] = oneSamplePt[kk-1];
               mSamplePts[jj*nInputs+ii] = samplePts1D[iL];
            }

            faPtr->evaluatePoint(nSubSamples, mSamplePts, YY);

            for (jj = 0; jj < nSubSamples; jj++)
            {
               oneSamplePt = &(mSamplePts[jj*nInputs]);
               ddata = constrPtr->evaluate(oneSamplePt, YY[jj], status);

               if (status == 0) YY[jj] = PSUADE_UNDEFINED;

               if (scatterFileFlag == 1 && iR == 0 && ddata != PSUADE_UNDEFINED)
               {
                  for (kk = 0; kk < nInputs; kk++)
                     fprintf(fp, "%e ", mSamplePts[jj*nInputs+kk]);
                  fprintf(fp, " %e\n", YY[jj]);
               }
            }

            means[iL] = 0.0;
            sCnt = 0;
            for (jj = 0; jj < nSubSamples; jj++)
            {
               if (YY[jj] != PSUADE_UNDEFINED)
               {
                  means[iL] += YY[jj];
                  sCnt++;
               }
            }
            bins[iL] = sCnt;
            if (sCnt < nSubSamples/10 && printLevel >= 5)
              printf("RSMSobol1 WARNING: subsample size = %d\n",sCnt); 
            if (sCnt >= 1) means[iL] /= (double) sCnt;
            else           means[iL] = PSUADE_UNDEFINED;
            if (printLevel > 3)
            {
               printf("RSMSobol1: input %d :\n", ii+1);
               printf("  refinement %2d, level %3d, size %d), mean = %e\n",iR,
                      iL+1, nSubSamples, means[iL]);
            }
            vars[iL] = 0.0;
            if (sCnt > 1)
            {
               for (jj = 0; jj < nSubSamples; jj++)
                  if (YY[jj] != PSUADE_UNDEFINED)
                     vars[iL] += (YY[jj] - means[iL]) * (YY[jj] - means[iL]);
               vars[iL] /= (double) sCnt;
            }
            else vars[iL] = PSUADE_UNDEFINED;
         }

         totalCnt = 0;
         for (iL = 0; iL < currNLevels; iL++) totalCnt += bins[iL];
         if (totalCnt == 0) 
         {
            printf("RSMSobol1 ERROR: no feasible region.\n"); 
            exit(1);
         }

         ddata = 0.0;
         for (iL = 0; iL < currNLevels; iL++) 
         {
            if (means[iL] != PSUADE_UNDEFINED)
               ddata += means[iL] * bins[iL] / totalCnt;
         }
         vces[ii] = 0.0;
         for (iL = 0; iL < currNLevels; iL++)
            if (means[iL] != PSUADE_UNDEFINED)
               vces[ii] += (means[iL] - ddata) * (means[iL] - ddata) * 
                           bins[iL] / totalCnt;
         if (printLevel > 2 || iR == 3)
         {
            printf("VCE(%3d) (refinement=%3d) = %e, (normalized) = %e\n",
                   ii+1, iR, vces[ii], vces[ii]/variance);
         }
         if (iR == 3) vces[ii] = vces[ii] / variance;

         ecvs[ii] = 0.0;
         for (iL = 0; iL < currNLevels; iL++)
         {
            if (vars[iL] != PSUADE_UNDEFINED)
               ecvs[ii] += vars[iL] * bins[iL] / totalCnt;
         }
         if (printLevel > 2)
         {
            printf("ECV(%3d) (refinement=%3d) = %e, (normalized) = %e\n",
                   ii+1, iR, ecvs[ii], ecvs[ii]/variance);
         }
         currNLevels *= 2;
      }
   }
   if (printLevel > 0)
   {
      printAsterisks(0);
      for (ii = 0; ii < nInputs; ii++)
         printf("RSMSobol1: Normalized VCE for input %3d = %e\n",ii+1,vces[ii]);
      for (ii = 0; ii < nInputs; ii++) means[ii] = (double) ii;
      printEquals(0);
      sortDbleList2(nInputs, vces, means);
      for (ii = nInputs-1; ii >= 0; ii--)
         printf("RSMSobol1: Normalized VCE (ordered) for input %3d = %e\n",
                (int) means[ii]+1,vces[ii]);
      printAsterisks(0);
   }

   if (scatterFileFlag == 1)
   {
      fprintf(fp, "];\n");
      fprintf(fp, "for ii = 1 : %d\n", nInputs);
      fprintf(fp, "   plot(A(:,ii),A(:,%d),'x')\n", nInputs+1);
      fwritePlotAxes(fp);
      fwritePlotXLabel(fp, "['input ' int2str(ii)]");
      fwritePlotYLabel(fp, "Y");
      fprintf(fp, "   disp('Press enter to continue')\n");
      fprintf(fp, "   pause\n");
      fprintf(fp, "end\n");
      fclose(fp);
   }
    
   pData *pPtr = ioPtr->getAuxData();
   pPtr->nDbles_ = nInputs;
   pPtr->dbleArray_ = new double[nInputs];
   for (ii = 0; ii < nInputs; ii++)
      pPtr->dbleArray_[ii] = vces[ii] * variance;
   pPtr->dbleData_ = variance;

   delete [] pdfFlags2;
   delete faPtr;
   delete constrPtr;
   delete [] cLower;
   delete [] cUpper;
   delete [] XX;
   delete [] YY;
   delete [] Y;
   delete [] samplePts1D;
   delete [] mSamplePts;
   delete [] inputMeans2;
   delete [] inputStdevs2;
   delete [] means;
   delete [] vars;
   delete [] vces;
   delete [] ecvs;
   delete [] bins;
   return 0.0;
}

// ************************************************************************
// perform analysis (for problems with joint PDFs)
// ------------------------------------------------------------------------
double RSMSobol1Analyzer::analyze2(aData &adata)
{
   int    nInputs, nOutputs, nSamples, *S, nSteps, outputID;
   int    ii, jj, iL, iR, sCnt, length, status, currNLevels;
   int    nSubSamples=2000, nLevels=500, pdfFileFlag=0, *bins, bin;
   int    printLevel, nSamp, scatterFileFlag=0, totalCnt;
   double *X, *Y, *Y2, *XX, *YY, dmean, variance, *xLower, *xUpper, *means;
   double *vces, *vars, *ecvs, ddata, width, *oneSamplePt;
   char   pdfFile[500], scatterFile[500], pString[500], winput1[500];
   FILE   *fp;
   PsuadeData    *ioPtr;
   FuncApprox    *faPtr;
   RSConstraints *constrPtr;
   Vector        vecIn, vecOut, vecUB, vecLB;
   PDFManager    *pdfman;

   printAsterisks(0);
   printf("RSMSobol1: since joint PDFs have been specified, a different\n");
   printf("           interaction analysis will be performed.\n");
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
   Y           = new double[nSamples];
   for (ii = 0; ii < nSamples; ii++) Y[ii] = Y2[ii*nOutputs+outputID];

   constrPtr = new RSConstraints();
   constrPtr->genConstraints(ioPtr);

   faPtr = genFAInteractive(ioPtr, 0);
   length = -999;
   if(faPtr == NULL){
     printf("faPtr is NULL in file %s line %d aborting. \n", __FILE__, __LINE__);
     abort();
   }else{
     status = faPtr->genNDGridData(X, Y, &length, NULL, NULL);
   }

   printAsterisks(0);
   if (psAnaExpertMode_ == 1)
   {
      printAsterisks(0);
      printf("* RSMSobol1 creates one sample of size M \n");
      printf("* for each of the K input levels. Therefore,\n");
      printf("* the total sample size is\n");
      printf("* N = M * K.\n");
      printf("* NOW, nInputs = %d\n", nInputs);
      printf("* As a user, please decide on M and K.\n\n");
      printEquals(0);

      sprintf(pString,"Enter M (suggestion: 1000-10000) : ");
      nSubSamples = getInt(100, 1000000, pString);
      if (nSubSamples > 100000)
         printf("An M of %d may take very long time.\n", nSubSamples);

      sprintf(pString,"Enter K (suggestion: 100 - 1000) : ");
      nLevels = getInt(10, 10000, pString);
      if (nLevels > 1000)
         printf("* A K of %d may take very long time.\n", nLevels);
      printAsterisks(0);
   }
   else
   {
      nSubSamples = 10000;
      nLevels     = 1000;
      printf("* RSMSobol1: default M = %d.\n", nSubSamples);
      printf("* RSMSobol1: default K = %d.\n", nLevels);
      printf("* To change these settings, turn on ana_expert mode and rerun.\n");
   }
   printEquals(0);

   if (psAnaExpertMode_ == 1)
   {
      printf("* RSMSobol1 will create a sample for basic statistics. You\n");
      printf("* have the option to plot the probability density function.\n");
      sprintf(pString, "Create a pdf (probability) bar graph? (y or n) ");
      getString(pString, winput1);
      if (winput1[0] == 'y')
      {
         sprintf(pString,"Enter the file name for the bar graph : ");
         getString(pString, pdfFile);
         pdfFile[strlen(pdfFile)-1] = '\0';
         fp = fopen(pdfFile, "w");
         if (fp != NULL)
         {
            fclose(fp);
            if (psPlotTool_ == 0) pdfFileFlag = 1; 
         }
         else 
         {
            printf("RSMSobol1 WARNING: cannot open file %s\n", pdfFile);
            pdfFileFlag = 0; 
         }
      }
      printEquals(0);
   }
   else pdfFileFlag = 0;

   if (psAnaExpertMode_ == 1)
   {
      printf("RSMSobol1 will create many samples for Sobol1 analysis. You\n");
      printf("have the option to plot these sample data.\n");
      sprintf(pString, "Create a scatter plot for RSMSobol1? (y or n) ");
      getString(pString, winput1);
      if (winput1[0] == 'y')
      {
         sprintf(pString,"Enter the file name for the scatter plot : ");
         getString(pString, scatterFile);
         scatterFile[strlen(scatterFile)-1] = '\0';
         fp = fopen(scatterFile, "w");
         if (fp != NULL)
         {
            fclose(fp);
            scatterFileFlag = 1;
         }
         else 
         {
            printf("RSMSobol1 ERROR: cannot open file %s\n", 
                   scatterFile); 
            scatterFileFlag = 0;
         }
      }
      printEquals(0);
   }
   else scatterFileFlag = 0;

   nSamp = nSubSamples * nLevels;
   if (nSamp < 100000) nSamp = 100000;
   if (printLevel > 1)
   {
      printf("* RSMSobol1 INFO: creating a sample for basic statistics.\n");
      printf("*                 sample size = %d\n", nSamp);
   }

   XX = new double[nSamp*nInputs];
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
      printf("RSMSobol1: running the sample with response surface...\n");
   }
   faPtr->evaluatePoint(nSamp, XX, YY);
   if (printLevel > 1)
   {
      printf("RSMSobol1: done running the sample with response surface.\n");
   }

   for (ii = 0; ii < nSamp; ii++)
   {
      oneSamplePt = &(XX[ii*nInputs]);
      ddata = constrPtr->evaluate(oneSamplePt,YY[ii],status);
      if (status == 0) YY[ii] = PSUADE_UNDEFINED;
   }

   sCnt  = 0;
   dmean = 0.0;
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
      printf("RSMSobol1 ERROR: too few samples that satisify the ");
      printf("constraints (%d out of %d).\n", sCnt, nSamp);
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
   printf("* RSMSobol1: sample mean     (based on %d points) = %e\n",
          sCnt, dmean);
   printf("* RSMSobol1: sample std dev  (based on %d points) = %e\n",
          sCnt, sqrt(variance));
   if (variance == 0.0) variance = 1.0;
   delete pdfman;
   printAsterisks(0);

   if (pdfFileFlag == 1 && (fp = fopen(pdfFile, "w")) != NULL)
   {
      if (psPlotTool_ == 1)
           fprintf(fp,"// sample PDF based on the response surface.\n");
      else fprintf(fp,"%% sample PDF based on the response surface.\n");
      fprintf(fp, "A = [\n");
      nSteps = 1;
      if (nSamp > 2000) nSteps = nSamp / 1000;
      for (ii = 0; ii < nSamp; ii++)
      {
         if (YY[ii] != PSUADE_UNDEFINED)
         {
            if ((pdfFileFlag == 1) && (ii % nSteps) == 0) 
               fprintf(fp,"   %e\n", YY[ii]);
         }
      }
      fprintf(fp, "];\n");
      fprintf(fp, "amax = max(A);\n");
      fprintf(fp, "amin = min(A);\n");
      fprintf(fp, "step = (amax - amin) / 10;\n");
      fprintf(fp, "numa = size(A,1);\n");
      fprintf(fp, "B = zeros(numa,1);\n");
      fprintf(fp, "for ii = 1 : numa\n");
      fprintf(fp, "   B(ii) = ceil(((A(ii) - amin))/step);\n");
      fprintf(fp, "end;\n");
      fprintf(fp, "X = zeros(10,1);\n");
      fprintf(fp, "Y = zeros(10,1);\n");
      fprintf(fp, "for ii = 1 : 10\n");
      fprintf(fp, "   X(ii) = amin + (ii - 1 + 0.5) * step;\n");
      fprintf(fp, "   [ia,ja,aa] = find(B==ii);\n");
      fprintf(fp, "   Y(ii) = length(ia);\n");
      fprintf(fp, "end;\n");
      fprintf(fp, "Y = Y / numa;\n");
      fprintf(fp, "bar(X,Y,1.0)\n");
      fwritePlotXLabel(fp, "Output Values");
      fwritePlotYLabel(fp, "Probabilities");
      fwritePlotAxes(fp);
      fclose(fp);
   }

   fp = NULL;
   if (scatterFileFlag == 1 && (fp = fopen(scatterFile, "w")) != NULL)
   {
     //fp = fopen(scatterFile, "w");
      fprintf(fp, "%% scatter plot for RSMSobol1 sample\n");
      fprintf(fp, "A = [ \n");
      nSteps = 1;
      if (nSamp > 2000) nSteps = nSamp / 1000;
      for (ii = 0; ii < nSamp; ii++)
      {
         if (YY[ii] != PSUADE_UNDEFINED)
         {
            if ((pdfFileFlag == 1) && (ii % nSteps) == 0) 
               fprintf(fp,"   %e\n", YY[ii]);
         }
      }
      fprintf(fp, "];\n");
      fprintf(fp, "for ii = 1 : %d\n", nInputs);
      fprintf(fp, "   plot(A(:,ii),A(:,%d),'x')\n", nInputs+1);
      fwritePlotXLabel(fp, "['input ' int2str(ii)]");
      fwritePlotYLabel(fp, "Y");
      fwritePlotAxes(fp);
      fprintf(fp, "   disp('Press enter to continue')\n");
      fprintf(fp, "   pause\n");
      fprintf(fp, "end\n");
      fclose(fp);
   }

   means = new double[nLevels];
   vars  = new double[nLevels];
   vces  = new double[nInputs];
   ecvs  = new double[nInputs];
   bins  = new int[nLevels];

   for (ii = 0; ii < nInputs; ii++)
   {
      if (printLevel >= 1)
         printf("RSMSobol1: processing input %d\n", ii+1);

      currNLevels = nLevels / 8;
      for (iR = 0; iR < 4; iR++)
      {
         width = (xUpper[ii] - xLower[ii]) / currNLevels;
         for (iL = 0; iL < currNLevels; iL++)
         {
            bins[iL] = 0;
            means[iL] = 0.0;
            vars[iL] = 0.0;
         }
         for (jj = 0; jj < nSamp; jj++)
         {
            if (YY[jj] != PSUADE_UNDEFINED)
            {
               ddata = XX[jj*nInputs+ii];
               bin = (int) ((ddata - xLower[ii]) / width); 
               if (bin == currNLevels) bin = currNLevels - 1;
               means[bin] += YY[jj];
               bins[bin]++;
            }
         }
         for (iL = 0; iL < currNLevels; iL++)
         {
            if (bins[iL] < nSubSamples/10 && printLevel >= 5)
               printf("RSMSobol1 WARNING: subsample size = %d\n",bins[iL]); 
            if (bins[iL] >= 1) means[iL] /= (double) bins[iL];
            else               means[iL] = PSUADE_UNDEFINED;
            if (printLevel > 3)
            {
               printf("RSMSobol1: input %d :\n", ii+1);
               printf("  refinement %2d, level %3d, size %d), mean = %e\n",iR,
                      iL, nSubSamples, means[iL]);
            }
         }
         for (jj = 0; jj < nSamp; jj++)
         {
            ddata = XX[jj*nInputs+ii];
            bin = (int) ((ddata - xLower[ii]) / width); 
            if (bin == currNLevels) bin = currNLevels - 1;
            if (YY[jj] != PSUADE_UNDEFINED)
            {
               vars[bin] += (YY[jj] - means[bin]) * (YY[jj] - means[bin]);
            }
         }
         for (iL = 0; iL < currNLevels; iL++)
         {
            if (bins[iL] >= 1) vars[iL] /= (double) bins[iL];
            else               vars[iL] = PSUADE_UNDEFINED;
         }

         totalCnt = 0;
         for (iL = 0; iL < currNLevels; iL++) totalCnt += bins[iL];
         if (totalCnt == 0) 
         {
            printf("RSMSobol1 ERROR: no feasible region.\n"); 
            exit(1);
         }

         ddata = 0.0;
         for (iL = 0; iL < currNLevels; iL++) 
         {
            if (means[iL] != PSUADE_UNDEFINED)
               ddata += means[iL] * bins[iL] / totalCnt;
         }
         vces[ii] = 0.0;
         for (iL = 0; iL < currNLevels; iL++)
            if (means[iL] != PSUADE_UNDEFINED)
               vces[ii] += (means[iL] - ddata) * (means[iL] - ddata) * 
                           bins[iL] / totalCnt;
         if (printLevel > 2 || iR == 3)
         {
            printf("VCE(%3d) (refinement=%3d) = %e, (normalized) = %e\n",
                   ii+1, iR, vces[ii], vces[ii]/variance);
         }
         if (iR == 3) vces[ii] = vces[ii] / variance;

         ecvs[ii] = 0.0;
         for (iL = 0; iL < currNLevels; iL++)
         {
            if (vars[iL] != PSUADE_UNDEFINED)
               ecvs[ii] += vars[iL] * bins[iL] / totalCnt;
         }
         if (printLevel > 2)
         {
            printf("ECV(%3d) (refinement=%3d) = %e, (normalized) = %e\n",
                   ii+1, iR, ecvs[ii], ecvs[ii]/variance);
         }
         currNLevels *= 2;
      }
   }
   for (ii = 0; ii < nInputs; ii++)
      printf("RSMSobol1: Normalized VCE for input %3d = %e\n",ii+1,vces[ii]);

   delete faPtr;
   delete constrPtr;
   delete [] XX;
   delete [] YY;
   delete [] Y;
   delete [] means;
   delete [] vars;
   delete [] vces;
   delete [] ecvs;
   delete [] bins;
   return 0.0;
}

// ************************************************************************
// perform analysis using fuzzy evaluation (12/2012)
// ------------------------------------------------------------------------
double RSMSobol1Analyzer::analyze3(aData &adata)
{
   int    ii, jj, nn, count, length, status, *SS, nSamp;
   int    nSubSamples=2000, nLevels=500, ntimes=1;
   double *Y, dmean, variance, ddata;
   char   pString[500], *cString, winput1[500], winput2[500];;
   FuncApprox    *faPtr;
   RSConstraints *constrPtr;
   Vector        vecIn, vecOut, vecUB, vecLB;
   pData         pCorMat;
   Matrix        *corMatp, corMat;
   PDFManager    *pdfman;
   Sampling      *sampler;

   printAsterisks(0);
   printf("*          RS-based First Order Sobol' Indices \n");
   printEquals(0);
   printf("* TO GAIN ACCESS TO DIFFERENT OPTIONS: SET\n");
   printf("*\n");
   printf("* - ana_expert mode to finetune RSMSobol1 parameters, \n");
   printf("*   (e.g. sample size for integration can be adjusted).\n");
   printf("* - rs_expert to mode finetune response surface for RSMSobol1,\n");
   printf("* - printlevel to 1 to display more information.\n");
   printf("* - ntimes to 100 to compute also error bars for the results\n");
   printEquals(0);
   
   int    nInputs, nOutputs, nSamples, *S, outputID, printLevel;
   int    noPDF=1, *pdfFlags;
   double *xLower, *xUpper, *X, *Y2, *inputMeans, *inputStdevs;
   PsuadeData *ioPtr;

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
   if (noPDF == 1) printf("* RSMSobol1 INFO: all uniform distributions.\n");
   else
   {
      printf("RSMSobol1 INFO: non-uniform distributions detected\n");
      printf("                which will be used in this analysis.\n");
   }

   if (nInputs <= 0 || nSamples <= 0 || nOutputs <= 0) 
   {
      printf("RSMSobol1 ERROR: invalid arguments.\n");
      printf("   nInputs  = %d\n", nInputs);
      printf("   nOutputs = %d\n", nOutputs);
      printf("   nSamples = %d\n", nSamples);
      return PSUADE_UNDEFINED;
   } 
   if (outputID >= nOutputs || outputID < 0)
   {
      printf("RSMSobol1 ERROR: invalid output ID (%d).\n", outputID);
      return PSUADE_UNDEFINED;
   }
   if (nInputs <= 1)
   {
      printf("RSMSobol1 ERROR: nInputs=1 does not need this analysis.\n");
      return PSUADE_UNDEFINED;
   }
   if (ioPtr == NULL)
   {
      printf("RSMSobol1 ERROR: no data (PsuadeData).\n");
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
            printf("RSMSobol1 INFO: this method cannot handle correlated\n");
            printf("          inputs using joint PDFs. You can re-run\n");
            printf("          using the group variance-based method.\n");
            return -1;
         }
      }
   }
   status = 0;
   for (ii = 0; ii < nSamples; ii++)
      if (Y2[nOutputs*ii+outputID] > 0.9*PSUADE_UNDEFINED) status = 1;
   if (status == 1)
   {
      printf("RSMSobol1 ERROR: Some outputs are undefined. Prune the\n");
      printf("                 undefined sample point first.\n");
      return PSUADE_UNDEFINED;
   }
   Y = new double[nSamples];
   for (ii = 0; ii < nSamples; ii++) Y[ii] = Y2[ii*nOutputs+outputID];

   constrPtr = new RSConstraints();
   constrPtr->genConstraints(ioPtr);

   faPtr = genFAInteractive(ioPtr, 0);
   length = -999;
   status = faPtr->genNDGridData(X, Y, &length, NULL, NULL);
   delete [] Y;

   if (psAnaExpertMode_ == 1)
   {
      printAsterisks(0);
      printf("* RSMSobol1 creates one sample of size M \n");
      printf("* for each of the K input levels. Therefore,\n");
      printf("* the total sample size is\n");
      printf("* N = M * K * nInputs\n");
      nSubSamples = 1000;
      nLevels = 200;
      printf("* nInputs m = %d\n", nInputs);
      printf("* default M = %d\n", nSubSamples);
      printf("* default K = %d\n", nLevels);
      printf("* As a user, please decide on M and K.\n");
      printf("* Note: large M and K may take very long time.\n");
      printEquals(0);
      sprintf(pString,"Enter M (suggestion: 1000-10000) : ");
      nSubSamples = getInt(1000, 10000, pString);
      sprintf(pString,"Enter K (suggestion: 100 - 1000) : ");
      nLevels = getInt(100, 1000, pString);
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
      nLevels = 200;
      if (psConfig_ != NULL)
      {
         cString = psConfig_->getParameter("RSMSobol1_nsubsamples");
         if (cString != NULL)
         {
            sscanf(cString, "%s %s %d", winput1, winput2, &nSubSamples);
            if (nSubSamples < 1000)
            {
               printf("RSMSobol1 INFO: nSubSamples should be >= 1000.\n");
               nSubSamples = 1000;
            }
            else
            {
               printf("RSMSobol1 INFO: nSubSamples = %d (config).\n",nSubSamples);
            }
         }
         cString = psConfig_->getParameter("RSMSobol1_nlevels");
         if (cString != NULL)
         {
            sscanf(cString, "%s %s %d", winput1, winput2, &nLevels);
            if (nLevels < 200)
            {
               printf("RSMSobol1 INFO: nLevels should be >= 200.\n");
               nLevels = 200;
            }
            else
            {
               printf("RSMSobol1 INFO: nLevels = %d (config).\n",nLevels);
            }
         }
      }
      if (printLevel > 1)
      {
         printf("* RSMSobol1 creates one sample of size M \n");
         printf("* for each of the K input levels. Therefore,\n");
         printf("* the total sample size is\n");
         printf("* N = M * K * nInputs\n");
         printf("* nInputs m = %d\n", nInputs);
         printf("* default M = %d\n", nSubSamples);
         printf("* default K = %d\n", nLevels);
         printf("* To change these settings, turn on ana_expert mode and rerun.\n");
      }
      printEquals(0);
   }

   nSamp = 100000;
   if (printLevel > 1)
   {
      printf("RSMSobol1 INFO: creating a sample for basic statistics.\n");
      printf("                sample size = %d\n", nSamp);
   }

   double *XX, *YY;
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
      else
         sampler = SamplingCreateFromID(PSUADE_SAMP_LPTAU);
      sampler->setInputBounds(nInputs, xLower, xUpper);
      sampler->setOutputParams(1);
      sampler->setSamplingParams(nSamp, 1, 1);
      sampler->initialize(0);
      SS = new int[nSamp];
      sampler->getSamples(nSamp, nInputs, 1, XX, YY, SS);
      delete [] SS;
      delete sampler;
   }

   double *YS = new double[nSamp];
   if (printLevel > 1)
   {
      printf("RSMSobol1: running the sample with response surface...\n");
   }
   if (ntimes == 1)
   {
      faPtr->evaluatePoint(nSamp, XX, YY);
      for (ii = 0; ii < nSamp; ii++) YS[ii] = 0.0;
   }
   else faPtr->evaluatePointFuzzy(nSamp, XX, YY, YS);
   if (printLevel > 1)
   {
      printf("RSMSobol1: done running the sample with response surface.\n");
   }

   double    *oneSamplePt;
   double    *YYS = new double[ntimes];
   double    *YYm = new double[ntimes];
   double    *YYs = new double[ntimes];
   int       *YYc = new int[ntimes];
   PDFNormal **rsPDFs = new PDFNormal*[nSamp];
   for (ii = 0; ii < nSamp; ii++)
   {
      if (YS[ii] != 0) rsPDFs[ii] = new PDFNormal(YY[ii],YS[ii]); 
      else             rsPDFs[ii] = NULL;
   }
   
   for (nn = 0; nn < ntimes; nn++)
   {
      YYm[nn] = 0.0;
      YYc[nn] = 0;
      YYs[nn] = 0;
   }
   for (ii = 0; ii < nSamp; ii++)
   {
      if (rsPDFs[ii] != NULL)
      {
         rsPDFs[ii]->genSample(ntimes,YYS,YY[ii]-4*YS[ii],YY[ii]+4*YS[ii]);
      }
      else 
      {
         for (nn = 0; nn < ntimes; nn++) YYS[nn] = YY[ii];
      }
      oneSamplePt = &(XX[ii*nInputs]);
      for (nn = 0; nn < ntimes; nn++)
      {
         ddata = constrPtr->evaluate(oneSamplePt,YYS[nn],status);
         if (status != 0)
         {
            YYm[nn] += YYS[nn];
            YYc[nn]++;
         }
      }
   }
   for (nn = 0; nn < ntimes; nn++)
      if (YYc[nn] != 0) YYm[nn] /= (double) YYc[nn];
   for (ii = 0; ii < nSamp; ii++)
   {
      if (rsPDFs[ii] != NULL)
      {
         rsPDFs[ii]->genSample(ntimes,YYS,YY[ii]-4*YS[ii],YY[ii]+4*YS[ii]);
      }
      else
      {
         for (nn = 0; nn < ntimes; nn++) YYS[nn] = YY[ii];
      }
      oneSamplePt = &(XX[ii*nInputs]);
      for (nn = 0; nn < ntimes; nn++)
      {
         ddata = constrPtr->evaluate(oneSamplePt,YYS[nn],status);
         if (status != 0)
            YYs[nn] += pow(YYS[nn] - YYm[nn], 2.0);
      }
   }
   for (ii = 0; ii < nSamp; ii++)
      if (rsPDFs[ii] != NULL) delete rsPDFs[ii];
   delete rsPDFs;
   
   count  = 0;
   for (nn = 0; nn < ntimes; nn++) count += YYc[nn];
   printf("RSMSobol1 INFO: %6.2f percent passes the contrained filters.\n",
          (double) count * 100.0 /((double) ntimes*nSamp));
   for (nn = 0; nn < ntimes; nn++)
   {
      if (YYc[nn] > 0) YYs[nn] = sqrt(YYs[nn]/YYc[nn]);
      else             YYs[nn] = 0.0;
   }
   double smean = 0.0;
   for (nn = 0; nn < ntimes; nn++) smean += YYs[nn];
   smean /= (double) ntimes;
   double sstd = 0.0;
   for (nn = 0; nn < ntimes; nn++)
      sstd += (YYs[nn] - smean) * (YYs[nn] - smean) ;
   sstd = sqrt(sstd/ntimes);

   dmean = 0.0;
   for (nn = 0; nn < ntimes; nn++) dmean += YYm[nn];
   dmean /= (double) ntimes;
   variance = 0.0;
   for (nn = 0; nn < ntimes; nn++)
      variance += (YYm[nn] - dmean) * (YYm[nn] - dmean) ;
   variance /= (double) ntimes;

   printf("RSMSobol1: sample mean (std dev of mean) = %10.3e (%10.3e)\n",
          dmean, sqrt(variance));
   printf("RSMSobol1: std dev (std dev of std dev)  = %10.3e (%10.3e)\n",
          smean, sstd);
   if (smean == 0.0) smean = 1.0;
   delete [] XX;
   delete [] YY;
   delete [] YS;
   delete [] YYm;
   delete [] YYc;
   delete [] YYs;

   int    kk, *pdfFlags2, nSteps=1;
   double *inputMeans2, *inputStdevs2, *samplePts1D, *samplePtsND;
   double *cLower, *cUpper, *ZZ;
   PDFManager *pdfman1, *pdfman2;

   if (nSamp > 2000) nSteps = nSamp / 1000;
   cLower = new double[nInputs];
   cUpper = new double[nInputs];
   nSamp  = nSubSamples;
   samplePts1D   = new double[nLevels*nInputs];
   samplePtsND   = new double[nSubSamples*nInputs*nInputs];
   pdfFlags2     = new int[nInputs-1];
   inputMeans2   = new double[nInputs-1];
   inputStdevs2  = new double[nInputs-1];

   if (nLevels > nSubSamples) 
   {
      SS = new int[nLevels];
      ZZ = new double[nLevels];
   }
   else
   {
      SS = new int[nSubSamples];
      ZZ = new double[nSubSamples];
   }
   for (ii = 0; ii < nInputs; ii++)
   {
      if (noPDF == 0)
      {
         corMat.setDim(1,1);
         corMat.setEntry(0, 0, corMatp->getEntry(ii,ii));
         pdfman1 = new PDFManager();
         pdfman1->initialize(1,&pdfFlags[ii],&inputMeans[ii],
                             &inputStdevs[ii],corMat);
         vecLB.load(1, &xLower[ii]);
         vecUB.load(1, &xUpper[ii]);
         vecOut.setLength(nLevels);
         pdfman1->genSample(nLevels, vecOut, vecLB, vecUB);
         for (jj = 0; jj < nLevels; jj++)
            samplePts1D[ii*nLevels+jj] = vecOut[jj];
         delete pdfman1;
      }
      else
      {
         sampler = SamplingCreateFromID(PSUADE_SAMP_LHS);
         sampler->setInputBounds(1, &xLower[ii], &xUpper[ii]);
         sampler->setOutputParams(1);
         sampler->setSamplingParams(nLevels, 1, 0);
         sampler->initialize(0);
         sampler->getSamples(nLevels, 1, 1, 
                      &(samplePts1D[ii*nLevels]), ZZ, SS);
         delete sampler;
      }

      if (noPDF == 0)
      {
         corMat.setDim(nInputs-1, nInputs-1);
         for (jj = 0; jj < ii; jj++)
         {
            cLower[jj] = xLower[jj];
            cUpper[jj] = xUpper[jj];
            pdfFlags2[jj] = pdfFlags[jj];
            inputMeans2[jj] = inputMeans[jj];
            inputStdevs2[jj] = inputStdevs[jj];
            for (kk = 0; kk < ii; kk++)
               corMat.setEntry(jj, kk, corMatp->getEntry(jj,kk));
            for (kk = ii+1; kk < nInputs; kk++)
               corMat.setEntry(jj, kk-1, corMatp->getEntry(jj,kk));
         }
         for (jj = ii+1; jj < nInputs; jj++)
         {
            cLower[jj-1] = xLower[jj];
            cUpper[jj-1] = xUpper[jj];
            pdfFlags2[jj-1] = pdfFlags[jj];
            inputMeans2[jj-1] = inputMeans[jj];
            inputStdevs2[jj-1] = inputStdevs[jj];
            for (kk = 0; kk < ii; kk++)
               corMat.setEntry(jj-1, kk, corMatp->getEntry(jj,kk));
            for (kk = ii+1; kk < nInputs; kk++)
               corMat.setEntry(jj-1, kk-1, corMatp->getEntry(jj,kk));
         }
         pdfman2 = new PDFManager();
         pdfman2->initialize(nInputs-1,pdfFlags2,inputMeans2,
                             inputStdevs2,corMat);
         vecLB.load(nInputs-1, cLower);
         vecUB.load(nInputs-1, cUpper);
         vecOut.setLength(nSubSamples*(nInputs-1));
         pdfman2->genSample(nSubSamples, vecOut, vecLB, vecUB);
         for (jj = 0; jj < nSubSamples*(nInputs-1); jj++) 
            samplePtsND[ii*nSubSamples*nInputs+jj] = vecOut[jj];
         delete pdfman2;
      }
      else
      {
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
         if (nInputs-1 > 51)
              sampler = SamplingCreateFromID(PSUADE_SAMP_LHS);
         else sampler = SamplingCreateFromID(PSUADE_SAMP_LPTAU);
         sampler->setInputBounds(nInputs-1, cLower, cUpper);
         sampler->setOutputParams(1);
         sampler->setSamplingParams(nSubSamples, 1, 0);
         sampler->initialize(0);
         sampler->getSamples(nSubSamples, nInputs-1, 1, 
                      &(samplePtsND[ii*nSubSamples*nInputs]), ZZ, SS);
         delete sampler;
      }
   }
   delete [] SS;
   delete [] ZZ;
   delete [] cLower;
   delete [] cUpper;
   delete [] pdfFlags2;
   delete [] inputMeans2;
   delete [] inputStdevs2;

   int       iL, sCnt, totalCnt, *bins;
   double    *YZ, *vceMaxs, *vceMins, *vceMeds;
   double    *mSamplePts, *vces, *ecvs, *vars, *means;
   PDFNormal *rsPDF;
   YY            = new double[nSubSamples];
   YS            = new double[nSubSamples];
   YZ            = new double[nSubSamples*ntimes];
   vceMaxs       = new double[nInputs];
   vceMins       = new double[nInputs];
   vceMeds       = new double[nInputs];
   vars          = new double[nLevels*ntimes];
   vces          = new double[nInputs*ntimes];
   ecvs          = new double[nInputs*ntimes];
   bins          = new int[nLevels*ntimes];
   means         = new double[nLevels*ntimes];
   mSamplePts    = new double[nInputs*nSubSamples];
   for (ii = 0; ii < nInputs; ii++)
   {
      if (printLevel >= 0)
         printf("RSMSobol1: processing input %d\n", ii+1);

      for (iL = 0; iL < nLevels; iL++)
      {
         if (printLevel > 3 && (iL % (nLevels/10) == 0))
            printf("RSMSobol1: level %d (%d) \n", iL+1, nLevels);
         for (jj = 0; jj < nSubSamples; jj++)
         {
            oneSamplePt = &(samplePtsND[ii*nSubSamples*nInputs+jj*(nInputs-1)]);
            for (kk = 0; kk < ii; kk++)
               mSamplePts[jj*nInputs+kk] = oneSamplePt[kk];
            for (kk = ii+1; kk < nInputs; kk++)
               mSamplePts[jj*nInputs+kk] = oneSamplePt[kk-1];
            mSamplePts[jj*nInputs+ii] = samplePts1D[ii*nLevels+iL];
         }

         if (printLevel > 3 && (iL % (nLevels/10) == 0))
            printf("RSMSobol1: fuzzy function evaluation %d (%d)\n",iL+1,nLevels);
         if (ntimes == 1)
         {
            faPtr->evaluatePoint(nSubSamples,mSamplePts,YY);
            for (jj = 0; jj < nSubSamples; jj++) YS[jj] = 0.0;
         }
         else faPtr->evaluatePointFuzzy(nSubSamples,mSamplePts,YY,YS);

         if (printLevel > 3 && (iL % (nLevels/10) == 0))
            printf("RSMSobol1: output perturbation\n");
         for (jj = 0; jj < nSubSamples; jj++)
         {
            if (YS[jj] != 0)
            {
               rsPDF = new PDFNormal(YY[jj],YS[jj]); 
               rsPDF->genSample(ntimes,&(YZ[jj*ntimes]),YY[jj]-4*YS[jj],
                                YY[jj]+4*YS[jj]);
               delete rsPDF;
            }
            else 
            {
               for (nn = 0; nn < ntimes; nn++) YZ[jj*ntimes+nn] = YY[jj];
            }
         }

         if (printLevel > 3 && (iL % (nLevels/10) == 0))
            printf("RSMSobol1: further processing\n");
         for (jj = 0; jj < nSubSamples; jj++)
         {
            oneSamplePt = &(mSamplePts[jj*nInputs]);
            for (nn = 0; nn < ntimes; nn++)
            {
               ddata = constrPtr->evaluate(oneSamplePt, YZ[jj*ntimes+nn], status);
               if (status == 0) YZ[jj*ntimes+nn] = PSUADE_UNDEFINED;
            }
         }

         if (printLevel > 3 && (iL % (nLevels/10) == 0))
            printf("RSMSobol1: compute mean and std dev\n");
         for (nn = 0; nn < ntimes; nn++)
         {
            means[iL*ntimes+nn] = 0.0;
            sCnt = 0;
            for (jj = 0; jj < nSubSamples; jj++)
            {
               if (YZ[jj*ntimes+nn] != PSUADE_UNDEFINED)
               {
                  means[iL*ntimes+nn] += YZ[jj*ntimes+nn];
                  sCnt++;
               }
            }
            bins[iL*ntimes+nn] = sCnt;
            if (sCnt < nSubSamples/10 && printLevel >= 5)
               printf("RSMSobol1 WARNING: subsample size = %d\n",sCnt); 
            if (sCnt >= 1) means[iL*ntimes+nn] /= (double) sCnt;
            else           means[iL*ntimes+nn] = PSUADE_UNDEFINED;
            vars[iL*ntimes+nn] = 0.0;
            if (sCnt > 1)
            {
               for (jj = 0; jj < nSubSamples; jj++)
                  if (YZ[jj*ntimes+nn] != PSUADE_UNDEFINED)
                     vars[iL*ntimes+nn] += pow(YZ[jj*ntimes+nn]-means[iL*ntimes+nn],2.0);
               vars[iL*ntimes+nn] /= (double) sCnt;
            }
            else vars[iL*ntimes+nn] = PSUADE_UNDEFINED;
         }
      }

      for (nn = 0; nn < ntimes; nn++)
      {
         totalCnt = 0;
         for (iL = 0; iL < nLevels; iL++) totalCnt += bins[iL*ntimes+nn];
         if (totalCnt == 0) 
         {
            printf("RSMSobol1 ERROR: no feasible region.\n"); 
            exit(1);
         }

         ddata = 0.0;
         for (iL = 0; iL < nLevels; iL++) 
         {
            if (means[iL*ntimes+nn] != PSUADE_UNDEFINED)
               ddata += means[iL*ntimes+nn] * bins[iL*ntimes+nn] / totalCnt;
         }
         vces[ii*ntimes+nn] = 0.0;
         for (iL = 0; iL < nLevels; iL++)
            if (means[iL*ntimes+nn] != PSUADE_UNDEFINED)
               vces[ii*ntimes+nn] += pow(means[iL*ntimes+nn]-ddata,2.0) * 
                           bins[iL*ntimes+nn] / totalCnt;
         ecvs[ii*ntimes+nn] = 0.0;
         for (iL = 0; iL < nLevels; iL++)
         {
            if (vars[iL*ntimes+nn] != PSUADE_UNDEFINED)
               ecvs[ii*ntimes+nn] += vars[iL*ntimes+nn]*bins[iL*ntimes+nn]/totalCnt;
         }
      }
   }
   for (ii = 0; ii < nInputs; ii++)
   {
      vceMaxs[ii] = -PSUADE_UNDEFINED;
      vceMins[ii] =  PSUADE_UNDEFINED;
      vceMeds[ii] =  0.0;
      for (nn = 0; nn < ntimes; nn++)
      {
         if (vces[ii*ntimes+nn] > vceMaxs[ii]) vceMaxs[ii] = vces[ii*ntimes+nn];
         if (vces[ii*ntimes+nn] < vceMins[ii]) vceMins[ii] = vces[ii*ntimes+nn];
         vceMeds[ii] += vces[ii*ntimes+nn];
      }
      vceMeds[ii] /= (double) ntimes;
      if (smean != 0)
      {
         vceMeds[ii] /= (smean * smean);
         vceMaxs[ii] /= (smean * smean);
         vceMins[ii] /= (smean * smean);
      }
      else vceMeds[ii] = vceMaxs[ii] = vceMins[ii] = 0.0;
   }

   if (printLevel >= 0)
   {
      printAsterisks(0);
      for (ii = 0; ii < nInputs; ii++)
      {
         printf("RSMSobol1: Normalized mean VCE for input %3d = %12.4e",
                ii+1, vceMeds[ii]);
         if (ntimes > 1)
              printf(",bounds = [%12.4e, %12.4e]\n",vceMins[ii],vceMaxs[ii]);
         else printf("\n");
      }
   }

   pData *pPtr = ioPtr->getAuxData();
   if (ntimes == 1)
   {
      pPtr->nDbles_ = nInputs;
      pPtr->dbleArray_ = new double[nInputs];
      for (ii = 0; ii < nInputs; ii++)
         pPtr->dbleArray_[ii] = vceMeds[ii] * smean * smean;
      pPtr->dbleData_ = smean * smean;
   }
   else
   {
      pPtr->nDbles_ = 3*nInputs;
      pPtr->dbleArray_ = new double[nInputs*3];
      for (ii = 0; ii < nInputs; ii++)
         pPtr->dbleArray_[ii] = vceMeds[ii] * smean * smean;
      for (ii = 0; ii < nInputs; ii++)
         pPtr->dbleArray_[nInputs+ii] = vceMins[ii] * smean * smean;
      for (ii = 0; ii < nInputs; ii++)
         pPtr->dbleArray_[2*nInputs+ii] = vceMaxs[ii] * smean * smean;
      pPtr->dbleData_ = smean * smean;
   }

   if (printLevel >= 0)
   {
      for (ii = 0; ii < nInputs; ii++) means[ii] = (double) ii;
      printEquals(0);
      sortDbleList2(nInputs, vceMeds, means);
      for (ii = nInputs-1; ii >= 0; ii--)
         printf("RSMSobol1: Normalized mean VCE (ordered) for input %3d = %e\n",
                (int) means[ii]+1,vceMeds[ii]);
      printAsterisks(0);
   }

   delete [] mSamplePts;
   delete [] samplePts1D;
   delete [] samplePtsND;
   delete [] YY;
   delete [] YZ;
   delete [] YS;
   delete [] vceMaxs;
   delete [] vceMins;
   delete [] vceMeds;
   delete [] means;
   delete [] vars;
   delete [] vces;
   delete [] ecvs;
   delete [] bins;
   delete faPtr;
   delete constrPtr;
   return 0.0;
}
// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
RSMSobol1Analyzer& RSMSobol1Analyzer::operator=(const RSMSobol1Analyzer &)
{
   printf("RSMSobol1 operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

