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

#include "Util/PsuadeUtil.h"
#include "Util/sysdef.h"
#include "Util/Matrix.h"
#include "Util/Vector.h"
#include "DataIO/pData.h"
#include "RSMSobol1Analyzer.h"
#include "Samplings/Sampling.h"
#include "PDFLib/PDFManager.h"
#include "Samplings/RSConstraints.h"
#include "Main/Psuade.h"
#include "DataIO/PsuadeData.h"

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
   int    ii, jj, kk, iL, iR, sCnt, length, status, currNLevels;
   int    nSubSamples=2000, nLevels=500, pdfFileFlag=0, *bins, *pdfFlags2;
   int    printLevel, nSamp, scatterFileFlag=0, totalCnt, *pdfFlags;
   double *X, *Y, *Y2, *XX, *YY, *ZZ, *oneSamplePt, dmean;
   double *xLower, *xUpper, *cLower, *cUpper, *means, variance, *vces;
   double *vars, *ecvs, ddata, *mSamplePts, *inputMeans;
   double *inputStdevs, *inputMeans2, *inputStdevs2, *samplePts1D;
   char   winput1[500], pdfFile[500], scatterFile[500];
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
   if (noPDF == 1) printf("RSMSobol1 INFO: all uniform distributions.\n");
   else
   {
      printf("RSMSobol1 INFO: non-uniform distributions detected which\n");
      printf("                will be used in this analysis.\n");
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
   Y = new double[nSamples];
   for (ii = 0; ii < nSamples; ii++) Y[ii] = Y2[ii*nOutputs+outputID];
   
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

   constrPtr = new RSConstraints();
   constrPtr->genConstraints(ioPtr);

   faPtr = genFAInteractive(ioPtr, 0);
   length = -999;
   status = faPtr->genNDGridData(X, Y, &length, NULL, NULL);

   if (psAnaExpertMode_ == 1)
   {
      printAsterisks(0);
      printf("* RSMSobol1Analyzer creates one sample of size M \n");
      printf("* for each of the K input levels. Therefore,\n");
      printf("* the total sample size is\n");
      printf("* N = M * K * nInputs.\n");
      printf("* NOW, nInputs = %d\n", nInputs);
      printf("* As a user, please decide on M and K.\n");
      printEquals(0);

      sprintf(pString,"Enter M (suggestion: 5000-10000) : ");
      nSubSamples = getInt(1000, 1000000, pString);
      if (nSubSamples > 100000)
         printf("An M of %d may take very long time.\n", nSubSamples);

      sprintf(pString,"Enter K (suggestion: 100 - 1000) : ");
      nLevels = getInt(10, 10000, pString);
      if (nLevels > 500)
         printf("A K of %d may take very long time.\n", nLevels);
      printAsterisks(0);
   }
   else
   {
      nSubSamples = 10000;
      nLevels = 1000;
      if (printLevel > 1)
      {
         printf("RSMSobol1: default M = %d.\n", nSubSamples);
         printf("RSMSobol1: default K = %d.\n", nLevels);
         printf("To change these settings, turn on ana_expert mode and rerun.\n");
      }
      printEquals(0);
   }

   if (psConfig_ != NULL)
   {
      cString = psConfig_->getParameter("RSMSobol1_pdffile");
      if (cString != NULL)
      {
         strcpy(pdfFile, cString);
         pdfFileFlag = 1; 
      }
   }
   else if (psAnaExpertMode_ == 1)
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

   if (psConfig_ != NULL)
   {
      cString = psConfig_->getParameter("RSMSobol1_scatterfile");
      if (cString != NULL)
      {
         strcpy(pdfFile, cString);
         scatterFileFlag = 1; 
      }
   }
   else if (psAnaExpertMode_ == 1)
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
         fprintf(fp, "%% PDF for the sample based on the response surface.\n");
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
         fprintf(fp, "grid on\n");
         fprintf(fp, "box on\n");
         fprintf(fp,"xlabel('Output values','FontSize',12,'FontWeight','bold')\n");
         fprintf(fp,"ylabel('Probabilities','Fontsize',12,'FontWeight','bold')\n");
         fprintf(fp, "set(gca,'linewidth',2)\n");
         fprintf(fp, "set(gca,'fontweight','bold')\n");
         fprintf(fp, "set(gca,'fontsize',12)\n");
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
      fprintf(fp, "%% scatter plot for RSMSobol1 sample\n");
      fprintf(fp, "A = [ \n");
   }

   for (ii = 0; ii < nInputs; ii++)
   {
      if (printLevel > 1)
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
      fprintf(fp, "   grid on\n");
      fprintf(fp, "   box on\n");
      fprintf(fp,"    xlabel(['input ' int2str(ii)],");
      fprintf(fp, "'FontSize',12,'FontWeight','bold')\n");
      fprintf(fp,"    ylabel('Y','Fontsize',12,'FontWeight','bold')\n");
      fprintf(fp, "   set(gca,'linewidth',2)\n");
      fprintf(fp, "   set(gca,'fontweight','bold')\n");
      fprintf(fp, "   set(gca,'fontsize',12)\n");
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
   status = faPtr->genNDGridData(X, Y, &length, NULL, NULL);

   printAsterisks(0);
   if (psAnaExpertMode_ == 1)
   {
      printAsterisks(0);
      printf("* RSMSobol1Analyzer creates one sample of size M \n");
      printf("* for each of the K input intervals. Therefore,\n");
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
         printf("A K of %d may take very long time.\n", nLevels);
      printAsterisks(0);
   }
   else
   {
      nSubSamples = 10000;
      nLevels     = 1000;
      printf("RSMSobol1: default M = %d.\n", nSubSamples);
      printf("RSMSobol1: default K = %d.\n", nLevels);
      printf("To change these settings, turn on ana_expert mode and rerun.\n");
   }
   printEquals(0);

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
      printf("RSMSobol1 INFO: creating a sample for basic statistics.\n");
      printf("                sample size = %d\n", nSamp);
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
      return PSUADE_UNDEFINED;
   }
   variance = 0.0;
   for (ii = 0; ii < nSamp; ii++)
   {
      if (YY[ii] != PSUADE_UNDEFINED)
         variance += (YY[ii] - dmean) * (YY[ii] - dmean) ;
   }
   variance /= (double) sCnt;
   printf("RSMSobol1: sample mean     (based on %d points) = %e\n",
          sCnt, dmean);
   printf("RSMSobol1: sample std dev  (based on %d points) = %e\n",
          sCnt, sqrt(variance));
   if (variance == 0.0) variance = 1.0;
   delete pdfman;
   printAsterisks(0);

   if (pdfFileFlag == 1)
   {
      fp = fopen(pdfFile, "w");
      fprintf(fp, "%% PDF for the sample based on the response surface.\n");
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
      fprintf(fp, "grid on\n");
      fprintf(fp, "box on\n");
      fprintf(fp,"xlabel('Output values','FontSize',12,'FontWeight','bold')\n");
      fprintf(fp,"ylabel('Probabilities','Fontsize',12,'FontWeight','bold')\n");
      fprintf(fp, "set(gca,'linewidth',2)\n");
      fprintf(fp, "set(gca,'fontweight','bold')\n");
      fprintf(fp, "set(gca,'fontsize',12)\n");
      fclose(fp);
   }

   fp = NULL;
   if (scatterFileFlag == 1)
   {
      fp = fopen(scatterFile, "w");
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
      fprintf(fp, "   grid on\n");
      fprintf(fp, "   box on\n");
      fprintf(fp,"    xlabel(['input ' int2str(ii)],");
      fprintf(fp, "'FontSize',12,'FontWeight','bold')\n");
      fprintf(fp,"    ylabel('Y','Fontsize',12,'FontWeight','bold')\n");
      fprintf(fp, "   set(gca,'linewidth',2)\n");
      fprintf(fp, "   set(gca,'fontweight','bold')\n");
      fprintf(fp, "   set(gca,'fontsize',12)\n");
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
      if (printLevel > 1)
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

