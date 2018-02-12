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
#include "PrintingTS.h"

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
RSMSobol1Analyzer::RSMSobol1Analyzer() : Analyzer(),nInputs_(0),
                     outputMean_(0), outputStd_(0), vces_(0)
{
  setName("RSMSOBOL1");
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
RSMSobol1Analyzer::~RSMSobol1Analyzer()
{
  if (vces_ != NULL) delete [] vces_;
}

// ************************************************************************
// perform analysis
// ------------------------------------------------------------------------
void RSMSobol1Analyzer::analyze(int nInps, int nSamp, double *lbs, 
                                double *ubs, double *X, double *Y)
{
  aData adata;
  adata.nInputs_ = nInps;
  adata.nOutputs_ = 1;
  adata.nSamples_ = nSamp;
  adata.iLowerB_ = lbs;
  adata.iUpperB_ = ubs;
  adata.sampleInputs_ = X;
  adata.sampleOutputs_ = Y;
  adata.outputID_ = 0;
  adata.printLevel_ = 0;
  analyze3(adata);
}

// ************************************************************************
// perform analysis
// ------------------------------------------------------------------------
double RSMSobol1Analyzer::analyze(aData &adata)
{
  int    nInputs, nOutputs, nSamples, *SS, nSteps, outputID, noPDF=1;
  int    ii, jj, kk, iL, iR, sCnt, status, currNLevels, method=1, corFlag;
  int    nSubSamples=1000, nLevels=200, pdfFileFlag=0, *bins, *pdfFlags2;
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
  psVector      vecIn, vecOut, vecUB, vecLB;
  pData         pCorMat;
  psMatrix      *corMatp, corMat;
  PDFManager    *pdfman, *pdfman1, *pdfman2;
  Sampling      *sampler;

  if (vces_ != NULL) delete [] vces_;
  vces_ = NULL;
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
  outputID    = adata.outputID_;
  ioPtr       = adata.ioPtr_;
  printLevel  = adata.printLevel_;
  pdfFlags    = adata.inputPDFs_;
  inputMeans  = adata.inputMeans_;
  inputStdevs = adata.inputStdevs_;
  corFlag = 0;
  if (pdfFlags != NULL)
  {
    for (ii = 0; ii < nInputs; ii++)
      if (pdfFlags[ii] != 0) noPDF = 0;
    for (ii = 0; ii < nInputs; ii++)
      if (pdfFlags[ii] == PSUADE_PDF_SAMPLE) corFlag = 1;
  } 
  printAsterisks(PL_INFO, 0);
  printOutTS(PL_INFO, "* RSMSobol1 constructor\n");
  printDashes(PL_INFO, 0);
  if (noPDF == 1) 
    printOutTS(PL_INFO,"* RSMSobol1 INFO: all uniform distributions.\n");
  else
  {
    printOutTS(PL_INFO,"* RSMSobol1 INFO: non-uniform distributions");
    printOutTS(PL_INFO," detected.\n");
  }

  if (nInputs <= 0 || nSamples <= 0 || nOutputs <= 0) 
  {
    printOutTS(PL_ERROR, "RSMSobol1 ERROR: invalid arguments.\n");
    printOutTS(PL_ERROR, "   nInputs  = %d\n", nInputs);
    printOutTS(PL_ERROR, "   nOutputs = %d\n", nOutputs);
    printOutTS(PL_ERROR, "   nSamples = %d\n", nSamples);
    return PSUADE_UNDEFINED;
  } 
  if (outputID >= nOutputs || outputID < 0)
  {
    printOutTS(PL_ERROR,"RSMSobol1 ERROR: invalid output ID (%d).\n", 
               outputID);
    return PSUADE_UNDEFINED;
  }
  if (nInputs <= 1)
  {
    printOutTS(PL_ERROR,"RSMSobol1 ERROR: nInputs=1 does not need");
    printOutTS(PL_ERROR," this analysis.\n");
    return PSUADE_UNDEFINED;
  }
  if (ioPtr == NULL)
  {
    printOutTS(PL_ERROR,
         "RSMSobol1 ERROR: no data object (PsuadeData) found.\n");
    return PSUADE_UNDEFINED;
  } 
  ioPtr->getParameter("input_cor_matrix", pCorMat);
  corMatp = (psMatrix *) pCorMat.psObject_;
  for (ii = 0; ii < nInputs; ii++)
  {
    for (jj = 0; jj < ii; jj++)
    {
      if (corMatp->getEntry(ii,jj) != 0.0)
      {
        printOutTS(PL_INFO,"RSMSobol1 INFO: Correlated inputs have\n");
        printOutTS(PL_INFO,"   been detected. PSUADE will use a\n");
        printOutTS(PL_INFO,"   variant of this method.\n");
        corFlag = 1;
      }
    }
  }
  status = 0;
  for (ii = 0; ii < nSamples; ii++)
    if (Y2[nOutputs*ii+outputID] > 0.9*PSUADE_UNDEFINED) status = 1;
  if (status == 1)
  {
    printOutTS(PL_ERROR, 
       "RSMSobol1 ERROR: Some outputs are undefined. Prune the\n");
    printOutTS(PL_ERROR,
       "                 undefined sample point first.\n");
    return PSUADE_UNDEFINED;
  }
  if (corFlag != 0) return analyze2(adata);

  Y = new double[nSamples];
  checkAllocate(Y, "Y in RSMSobol1::analyze");
  for (ii = 0; ii < nSamples; ii++) Y[ii] = Y2[ii*nOutputs+outputID];

  constrPtr = new RSConstraints();
  constrPtr->genConstraints(ioPtr);

  faPtr = genFAInteractive(ioPtr, 0);
  status = faPtr->initialize(X, Y);

  if (psAnaExpertMode_ == 1)
  {
    printAsterisks(PL_INFO, 0);
    printOutTS(PL_INFO, "* RSMSobol1 creates one sample of size M \n");
    printOutTS(PL_INFO, "* for each of the K input levels. Therefore,\n");
    printOutTS(PL_INFO, "* the total sample size is\n");
    printOutTS(PL_INFO, "* N = M * K * nInputs\n");
    nSubSamples = 1000;
    nLevels = 200;
    printOutTS(PL_INFO, "* nInputs m = %d, and\n", nInputs);
    printOutTS(PL_INFO, "* default M = %d\n", nSubSamples);
    printOutTS(PL_INFO, "* default K = %d\n", nLevels);
    printOutTS(PL_INFO, "* As a user, please decide on M and K.\n");
    printOutTS(PL_INFO, "* Note: large M and K may take a long time\n");
    printEquals(PL_INFO, 0);
    sprintf(pString,"Enter M (suggestion: 1000-10000) : ");
    nSubSamples = getInt(1000, 50000, pString);
    sprintf(pString,"Enter K (suggestion: 100 - 1000) : ");
    nLevels = getInt(100, 5000, pString);
    printAsterisks(PL_INFO, 0);
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
          printOutTS(PL_INFO,
             "RSMSobol1 INFO: nSubSamples should be >= 1000.\n");
          nSubSamples = 1000;
        }
        else
        {
          printOutTS(PL_INFO, 
             "RSMSobol1 INFO: nSubSamples = %d (config).\n",nSubSamples);
        }
      }
      cString = psConfig_->getParameter("RSMSobol1_nlevels");
      if (cString != NULL)
      {
        sscanf(cString, "%s %s %d", winput1, winput2, &nLevels);
        if (nLevels < 200)
        {
          printOutTS(PL_INFO,
             "RSMSobol1 INFO: nLevels should be >= 200.\n");
          nLevels = 200;
        }
        else
        {
          printOutTS(PL_INFO,
             "RSMSobol1 INFO: nLevels = %d (config).\n",nLevels);
        }
      }
    }

    printAsterisks(PL_INFO, 0);
    printOutTS(PL_INFO,"* RSMSobol1 creates one sample of size M \n");
    printOutTS(PL_INFO,"* for each of the K input levels. Therefore,\n");
    printOutTS(PL_INFO,"* the total sample size is\n");
    printOutTS(PL_INFO,"* N = M * K * nInputs.\n");
    printOutTS(PL_INFO,"* nInputs m = %d \n", nInputs);
    printOutTS(PL_INFO,"* default M = %d\n", nSubSamples);
    printOutTS(PL_INFO,"* default K = %d\n", nLevels);
    printOutTS(PL_INFO,
         "* To change these settings, turn on ana_expert mode and rerun.\n");
    printAsterisks(PL_INFO, 0);
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
    printOutTS(PL_INFO,
         "RSMSobol1 will create a sample for basic statistics. You\n");
    printOutTS(PL_INFO,
         "have the option to plot the probability density function.\n");
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
        printOutTS(PL_INFO, 
             "RSMSobol1 WARNING: cannot open file %s\n", pdfFile);
        pdfFileFlag = 0; 
      }
    }
    printEquals(PL_INFO, 0);
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
    printOutTS(PL_INFO,
         "RSMSobol1 will create many samples for Sobol1 analysis.\n");
    printOutTS(PL_INFO,"You have the option to plot these sample data.\n");
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
        printOutTS(PL_INFO, "RSMSobol1 ERROR: cannot open file %s\n",
               scatterFile); 
        scatterFileFlag = 0;
      }
    }
    printEquals(PL_INFO, 0);
  }
  else scatterFileFlag = 0;

  nSamp = 100000;
  if (printLevel > 1)
  {
    printOutTS(PL_INFO,
         "RSMSobol1 INFO: creating a sample for basic statistics.\n");
    printOutTS(PL_INFO,"                sample size = %d\n", nSamp);
  }

  XX = new double[nSamp*nInputs];
  YY = new double[nSamp];
  checkAllocate(YY, "YY in RSMSobol1::analyze");

  if (noPDF == 0)
  {
    pdfman = new PDFManager();
    pdfman->initialize(nInputs,pdfFlags,inputMeans,
                       inputStdevs,*corMatp,NULL,NULL);
    vecLB.load(nInputs, xLower);
    vecUB.load(nInputs, xUpper);
    vecOut.setLength(nSamp*nInputs);
    pdfman->genSample(nSamp, vecOut, vecLB, vecUB);
    for (ii = 0; ii < nSamp*nInputs; ii++) XX[ii] = vecOut[ii];
    delete pdfman;
  }
  else
  {
    if (nInputs > 51) sampler = SamplingCreateFromID(PSUADE_SAMP_LHS);
    else              sampler = SamplingCreateFromID(PSUADE_SAMP_LPTAU);
    sampler->setInputBounds(nInputs, xLower, xUpper);
    sampler->setOutputParams(1);
    sampler->setSamplingParams(nSamp, 1, 1);
    sampler->initialize(0);
    SS = new int[nSamp];
    sampler->getSamples(nSamp, nInputs, 1, XX, YY, SS);
    delete [] SS;
    delete sampler;
  }

  printOutTS(PL_INFO,
          "RSMSobol1: running the sample with response surface...\n");
  faPtr->evaluatePoint(nSamp, XX, YY);
  printOutTS(PL_INFO,
          "RSMSobol1: done running the sample with response surface.\n");

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
    printOutTS(PL_ERROR, 
         "RSMSobol1 ERROR: too few samples that satisify the ");
    printOutTS(PL_ERROR, "constraints (%d out of %d).\n", sCnt, nSamp);
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
  printOutTS(PL_INFO, 
       "RSMSobol1: sample mean    (based on N = %d) = %10.3e\n",
       sCnt, dmean);
  printOutTS(PL_INFO, 
       "RSMSobol1: sample std dev (based on N = %d) = %10.3e\n",
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
  checkAllocate(inputStdevs2, "inputStdevs2 in RSMSobol1::analyze");

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
    printOutTS(PL_INFO, "RSMSobol1: processing input %d\n", ii+1);

    currNLevels = nLevels / 8;
    for (iR = 0; iR < 4; iR++)
    {
      if (noPDF == 0)
      {
        corMat.setDim(1,1);
        corMat.setEntry(0, 0, corMatp->getEntry(ii,ii));
        pdfman1 = new PDFManager();
        pdfman1->initialize(1,&pdfFlags[ii],&inputMeans[ii],
                            &inputStdevs[ii],corMat,NULL,NULL);
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
        checkAllocate(ZZ, "ZZ in RSMSobol1::analyze");
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
                            inputStdevs2,corMat,NULL,NULL);
        vecLB.load(nInputs-1, cLower);
        vecUB.load(nInputs-1, cUpper);
        vecOut.setLength(nSubSamples*(nInputs-1));
        pdfman2->genSample(nSubSamples, vecOut, vecLB, vecUB);
        for (jj = 0; jj < nSubSamples*(nInputs-1); jj++) 
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
        checkAllocate(ZZ, "ZZ(2) in RSMSobol1::analyze");
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
          printOutTS(PL_WARN, 
               "RSMSobol1 WARNING: subsample size = %d\n",sCnt);
        if (sCnt >= 1) means[iL] /= (double) sCnt;
        else           means[iL] = PSUADE_UNDEFINED;
        if (printLevel > 3)
        {
          printOutTS(PL_INFO,"RSMSobol1: input %d :\n", ii+1);
          printOutTS(PL_INFO,
               "  refinement %2d, level %3d, size %d), mean = %e\n",iR,
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
        printOutTS(PL_ERROR, "RSMSobol1 ERROR: no feasible region.\n");
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
        printOutTS(PL_INFO,
             "VCE(%3d) (refinement=%3d) = %e, (normalized) = %e\n",
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
        printOutTS(PL_INFO,
             "ECV(%3d) (refinement=%3d) = %e, (normalized) = %e\n",
             ii+1, iR, ecvs[ii], ecvs[ii]/variance);
      }
      currNLevels *= 2;
    }
  }
  if (printLevel > 0)
  {
    printAsterisks(PL_INFO, 0);
    for (ii = 0; ii < nInputs; ii++)
      printOutTS(PL_INFO, "RSMSobol1: Normalized VCE for input %3d = %e\n",
                 ii+1,vces[ii]);
    for (ii = 0; ii < nInputs; ii++) means[ii] = (double) ii;
    printEquals(PL_INFO, 0);
    sortDbleList2(nInputs, vces, means);
    for (ii = nInputs-1; ii >= 0; ii--)
      printOutTS(PL_INFO,
           "RSMSobol1: Normalized VCE (ordered) for input %3d = %e\n",
           (int) means[ii]+1,vces[ii]);
    printAsterisks(PL_INFO, 0);
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
    
  pData *pPtr;
  if (ioPtr != NULL)
  {
    pPtr = ioPtr->getAuxData();
    pPtr->nDbles_ = nInputs;
    pPtr->dbleArray_ = new double[nInputs];
    for (ii = 0; ii < nInputs; ii++)
      pPtr->dbleArray_[ii] = vces[ii] * variance;
    pPtr->dbleData_ = variance;
  }

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
   int    nInputs, nOutputs, nSamples, nSteps, outputID;
   int    ii, jj, iL, iR, sCnt, status, currNLevels;
   int    nSubSamples=2000, nLevels=500, pdfFileFlag=0, *bins, bin;
   int    printLevel, nSamp, scatterFileFlag=0, totalCnt;
   double *X, *Y, *Y2, *XX, *YY, dmean, variance, *xLower, *xUpper, *means;
   double *vces, *vars, *ecvs, ddata, width, *oneSamplePt;
   char   pdfFile[500], scatterFile[500], pString[500], winput1[500];
   FILE   *fp;
   PsuadeData    *ioPtr;
   FuncApprox    *faPtr;
   RSConstraints *constrPtr;
   psVector      vecIn, vecOut, vecUB, vecLB;
   PDFManager    *pdfman;

   printAsterisks(PL_INFO, 0);
   printOutTS(PL_INFO,"RSMSobol1: since joint PDFs have been specified,\n");
   printOutTS(PL_INFO,"           a different interaction analysis will\n");
   printOutTS(PL_INFO,"           be performed.\n");
   printEquals(PL_INFO, 0);
   nInputs     = adata.nInputs_;
   nOutputs    = adata.nOutputs_;
   nSamples    = adata.nSamples_;
   xLower      = adata.iLowerB_;
   xUpper      = adata.iUpperB_;
   X           = adata.sampleInputs_;
   Y2          = adata.sampleOutputs_;
   outputID    = adata.outputID_;
   ioPtr       = adata.ioPtr_;
   printLevel  = adata.printLevel_;
   Y           = new double[nSamples];
   checkAllocate(Y, "Y in RSMSobol1::analyze2");
   for (ii = 0; ii < nSamples; ii++) Y[ii] = Y2[ii*nOutputs+outputID];

   constrPtr = new RSConstraints();
   constrPtr->genConstraints(ioPtr);

   faPtr = genFAInteractive(ioPtr, 0);
   if(faPtr == NULL)
   {
      printOutTS(PL_INFO, "faPtr is NULL in file %s line %d aborting. \n", 
                 __FILE__, __LINE__);
      abort();
   }
   else status = faPtr->initialize(X, Y);

   printAsterisks(PL_INFO, 0);
   if (psAnaExpertMode_ == 1)
   {
      printAsterisks(PL_INFO, 0);
      printf("* RSMSobol1 creates one sample of size M \n");
      printf("* for each of the K input levels. Therefore,\n");
      printf("* the total sample size is\n");
      printf("* N = M * K.\n");
      printf("* NOW, nInputs = %d\n", nInputs);
      printf("* As a user, please decide on M and K.\n\n");
      printEquals(PL_INFO, 0);

      sprintf(pString,"Enter M (suggestion: 100-10000) : ");
      nSubSamples = getInt(100, 100000, pString);
      if (nSubSamples > 100000)
         printOutTS(PL_INFO, "An M of %d may take very long time.\n", 
                    nSubSamples);

      sprintf(pString,"Enter K (suggestion: 1000 - 10000) : ");
      nLevels = getInt(1000, 10000, pString);
      if (nLevels > 5000)
         printOutTS(PL_INFO, "* A K of %d may take very long time.\n", 
                    nLevels);
      printAsterisks(PL_INFO, 0);
   }
   else
   {
      nSubSamples = 100;
      nLevels     = 1000;
      printOutTS(PL_INFO,"* RSMSobol1: default M = %d.\n", nSubSamples);
      printOutTS(PL_INFO,"* RSMSobol1: default K = %d.\n", nLevels);
      printOutTS(PL_INFO,"* To change settings, rerun with ana_expert on.\n");
   }
   printEquals(PL_INFO, 0);

   if (psMasterMode_ == 1)
   {
      printOutTS(PL_INFO,"* RSMSobol1 will create a sample for basic\n");
      printOutTS(PL_INFO,"* statistics. You have the option to plot\n");
      printOutTS(PL_INFO,"* the PDF.\n");
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
            printOutTS(PL_INFO,"RSMSobol1 WARNING: cannot open file %s\n", 
                       pdfFile);
            pdfFileFlag = 0; 
         }
      }
      printEquals(PL_INFO, 0);
   }
   else pdfFileFlag = 0;

   if (psMasterMode_ == 1)
   {
      printOutTS(PL_INFO,"RSMSobol1 will create many samples for Sobol1\n");
      printOutTS(PL_INFO,"analysis. You have the option to plot these\n");
      printOutTS(PL_INFO,"sample data.\n");
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
            printOutTS(PL_ERROR, "RSMSobol1 ERROR: cannot open file %s\n",
                   scatterFile); 
            scatterFileFlag = 0;
         }
      }
      printEquals(PL_INFO, 0);
   }
   else scatterFileFlag = 0;

   nSamp = nSubSamples * nLevels;
   if (nSamp < 100000) nSamp = 100000;
   if (printLevel > 1)
   {
      printOutTS(PL_INFO,"* RSMSobol1 INFO: creating a sample for basic\n");
      printOutTS(PL_INFO,"*           statistics. Sample size = %d\n",nSamp);
   }

   XX = new double[nSamp*nInputs];
   YY = new double[nSamp];
   checkAllocate(YY, "YY in RSMSobol1::analyze2");

   pdfman = new PDFManager();
   pdfman->initialize(ioPtr);
   vecLB.load(nInputs, xLower);
   vecUB.load(nInputs, xUpper);
   vecOut.setLength(nSamp*nInputs);
   pdfman->genSample(nSamp, vecOut, vecLB, vecUB);
   for (ii = 0; ii < nSamp*nInputs; ii++) XX[ii] = vecOut[ii];

   printOutTS(PL_INFO,"RSMSobol1: response surface evaluation begins...\n");
   faPtr->evaluatePoint(nSamp, XX, YY);
   printOutTS(PL_INFO,"RSMSobol1: response surface evaluation ends...\n");

   for (ii = 0; ii < nSamp; ii++)
   {
      oneSamplePt = &(XX[ii*nInputs]);
      status = 1;
      if (constrPtr != NULL)
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
      printOutTS(PL_ERROR,"RSMSobol1 ERROR: too few samples that satisify");
      printOutTS(PL_ERROR,"the constraints (%d out of %d).\n", sCnt, nSamp);
      delete [] XX;
      delete [] YY;
      delete faPtr;
      delete pdfman;
      if (constrPtr != NULL) delete constrPtr;
      return PSUADE_UNDEFINED;
   }
   variance = 0.0;
   for (ii = 0; ii < nSamp; ii++)
   {
      if (YY[ii] != PSUADE_UNDEFINED)
         variance += (YY[ii] - dmean) * (YY[ii] - dmean) ;
   }
   variance /= (double) sCnt;
   printOutTS(PL_INFO,
        "* RSMSobol1: sample mean    (based on %d points) = %e\n",
        sCnt,dmean);
   printOutTS(PL_INFO,
        "* RSMSobol1: sample std dev (based on %d points) = %e\n",
        sCnt, sqrt(variance));
   if (variance == 0.0) variance = 1.0;
   delete pdfman;
   printAsterisks(PL_INFO, 0);

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
   checkAllocate(bins, "bins in RSMSobol1::analyze2");

   for (ii = 0; ii < nInputs; ii++)
   {

      printOutTS(PL_INFO, "RSMSobol1: processing input %d\n", ii+1);

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
               printOutTS(PL_INFO,
                    "RSMSobol1 WARNING: subsample size = %d\n",bins[iL]);
            if (bins[iL] >= 1) means[iL] /= (double) bins[iL];
            else               means[iL] = PSUADE_UNDEFINED;
            if (printLevel > 3)
            {
               printOutTS(PL_INFO,"RSMSobol1: input %d :\n", ii+1);
               printOutTS(PL_INFO,
                    "  refinement %2d, level %3d, size %d), mean = %e\n",
                    iR, iL, nSubSamples, means[iL]);
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
            printOutTS(PL_ERROR, "RSMSobol1 ERROR: no feasible region.\n");
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
            printOutTS(PL_INFO,
                 "VCE(%3d) (refinement=%3d) = %e, (normalized) = %e\n",
                 ii+1, iR, vces[ii], vces[ii]/variance);
         }
         if (iR == 3) vces[ii] = vces[ii] / variance;

         ecvs[ii] = 0.0;
         for (iL = 0; iL < currNLevels; iL++)
         {
            if (vars[iL] != PSUADE_UNDEFINED)
               ecvs[ii] += vars[iL] * bins[iL] / totalCnt;
         }

         printOutTS(PL_INFO, 
              "ECV(%3d) (refinement=%3d) = %e, (normalized) = %e\n",
              ii+1, iR, ecvs[ii], ecvs[ii]/variance);

         currNLevels *= 2;
      }
   }
   for (ii = 0; ii < nInputs; ii++)
      printOutTS(PL_INFO,"RSMSobol1: Normalized VCE for input %3d = %e\n",
                 ii+1,vces[ii]);

   delete faPtr;
   if (constrPtr != NULL) delete constrPtr;
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
   int    ii, jj, kk, nn, ss, count, status, nSamp, rstype;
   int    nSubSamples=1000, nLevels=200, ntimes=1, nConstr, pdfNull=0;
   double *Y, dmean, dstds, ddata, ddata2, frac=0.8;
   char   pString[500], *cString, winput1[500], winput2[500];
   FuncApprox    *faPtr;
   RSConstraints *constrPtr;
   psVector      vecIn, vecOut, vecUB, vecLB;
   pData         pCorMat, pPtr;
   psMatrix      *corMatp=NULL, corMat;
   PDFManager    *pdfman;
   Sampling      *sampler;

   if (isScreenDumpModeOn())
   {
     printAsterisks(PL_INFO, 0);
     printOutTS(PL_INFO,"*          RS-based First Order Sobol' Analysis \n");
     printEquals(PL_INFO, 0);
     printOutTS(PL_INFO,"* TO GAIN ACCESS TO DIFFERENT OPTIONS: SET\n");
     printOutTS(PL_INFO,"* - ana_expert mode to finetune RSMSobol1 parameters\n");
     printOutTS(PL_INFO,"*   (e.g. to adjust sample size for integration).\n");
     printOutTS(PL_INFO,"* - rs_expert mode to finetune response surface\n");
     printOutTS(PL_INFO,"* - printlevel to display more information\n");
     printOutTS(PL_INFO,"* Or, use configure file to finetune parameters\n");
     printEquals(PL_INFO, 0);
   }
 
   int    nInputs, nOutputs, nSamples, outputID, printLevel;
   int    noPDF=1, *pdfFlags, corFlag=0;
   double *xLower, *xUpper, *X, *Y2, *inputMeans, *inputStdevs;
   PsuadeData *ioPtr;

   nInputs     = adata.nInputs_;
   nInputs_    = nInputs;
   nOutputs    = adata.nOutputs_;
   nSamples    = adata.nSamples_;
   xLower      = adata.iLowerB_;
   xUpper      = adata.iUpperB_;
   X           = adata.sampleInputs_;
   Y2          = adata.sampleOutputs_;
   outputID    = adata.outputID_;
   ioPtr       = adata.ioPtr_;
   printLevel  = adata.printLevel_;
   pdfFlags    = adata.inputPDFs_;
   inputMeans  = adata.inputMeans_;
   inputStdevs = adata.inputStdevs_;
   if (inputMeans == NULL || pdfFlags == NULL || inputStdevs == NULL)
   {
      pdfNull = 1;
      pdfFlags    = new int[nInputs];
      inputMeans  = new double[nInputs];
      inputStdevs = new double[nInputs];
      checkAllocate(inputStdevs, "inputStdevs in RSMSobol1::analyze3");
      for (ii = 0; ii < nInputs; ii++)
      {
         pdfFlags[ii] = 0;
         inputMeans[ii]  = 0;
         inputStdevs[ii] = 0;
      }
   }
   if (pdfFlags != NULL)
   {
      for (ii = 0; ii < nInputs; ii++)
         if (pdfFlags[ii] != 0) noPDF = 0;
      for (ii = 0; ii < nInputs; ii++)
         if (pdfFlags[ii] == PSUADE_PDF_SAMPLE) corFlag = 1;
   } 
   if (isScreenDumpModeOn())
   {
      if (noPDF == 1) 
         printOutTS(PL_INFO,"* RSMSobol1 INFO: all uniform distributions.\n");
      else
      {
         printOutTS(PL_INFO,"RSMSobol1 INFO: non-uniform distributions");
         printOutTS(PL_INFO," detected.\n");
      }
   }

   if (nInputs <= 0 || nSamples <= 0 || nOutputs <= 0) 
   {
      printOutTS(PL_ERROR, "RSMSobol1 ERROR: invalid arguments.\n");
      printOutTS(PL_ERROR, "   nInputs  = %d\n", nInputs);
      printOutTS(PL_ERROR, "   nOutputs = %d\n", nOutputs);
      printOutTS(PL_ERROR, "   nSamples = %d\n", nSamples);
      return PSUADE_UNDEFINED;
   } 
   if (outputID >= nOutputs || outputID < 0)
   {
      printOutTS(PL_ERROR,"RSMSobol1 ERROR: invalid output ID (%d).\n", 
                 outputID);
      return PSUADE_UNDEFINED;
   }
   if (nInputs <= 1)
   {
      printOutTS(PL_ERROR,
           "RSMSobol1 INFO: analysis not needed for nInputs=1\n");
      return PSUADE_UNDEFINED;
   }
   if (ioPtr == NULL)
   {
      if (isScreenDumpModeOn())
      {
         printOutTS(PL_INFO,
           "RSMSobol1 INFO: no data object (PsuadeData) found.\n");
         printOutTS(PL_INFO,"   Several features will be turned off.\n");
      }
      corMatp = new psMatrix();
      corMatp->setDim(nInputs, nInputs);
      for (ii = 0; ii < nInputs; ii++) corMatp->setEntry(ii,ii,1.0e0);
   } 
   else
   {
      ioPtr->getParameter("input_cor_matrix", pCorMat);
      corMatp = (psMatrix *) pCorMat.psObject_;
      for (ii = 0; ii < nInputs; ii++)
      {
         for (jj = 0; jj < ii; jj++)
         {
            if (corMatp->getEntry(ii,jj) != 0.0)
            {
               if (isScreenDumpModeOn())
               {
                  printOutTS(PL_INFO, 
                    "RSMSobol1 INFO: Correlated inputs detected.\n");
                  printOutTS(PL_INFO, 
                    "   Alternative analyis is to be performed.\n");
               }
               corFlag = 1;
            }
         }
      }
   }
   status = 0;
   for (ii = 0; ii < nSamples; ii++)
      if (Y2[nOutputs*ii+outputID] > 0.9*PSUADE_UNDEFINED) status = 1;
   if (status == 1)
   {
      printOutTS(PL_ERROR,"RSMSobol1 ERROR: Some outputs are undefined.\n");
      printOutTS(PL_ERROR,"     Prune the undefined sample point first.\n");
      return PSUADE_UNDEFINED;
   }
   if (corFlag != 0) 
   {
      if (ioPtr == NULL) delete corMatp;
      if (pdfNull == 1)
      {
         delete [] pdfFlags;
         delete [] inputMeans;
         delete [] inputStdevs;
      }
      return analyze2(adata);
   }

   if (ioPtr != NULL)
   {
      constrPtr = new RSConstraints();
      constrPtr->genConstraints(ioPtr);
      nConstr = constrPtr->getNumConstraints();
   }
   else
   {
      nConstr = 0;
      constrPtr = NULL;
      if (isScreenDumpModeOn())
        printf("RSMSobolTSI INFO: no PsuadeData ==> no constraints.\n");
   }

   if (ioPtr == NULL)
   {
      if (rstype_ < 0)
      {
         printf("Select response surface. Options are: \n");
         writeFAInfo(0);
         strcpy(pString, "Choose response surface: ");
         rstype = getInt(0, PSUADE_NUM_RS, pString);
      }
      else rstype = rstype_;
   }
   else
   {
      ioPtr->getParameter("ana_rstype", pPtr);
      rstype = pPtr.intData_;
   }

   if (isScreenDumpModeOn() && psAnaExpertMode_ == 1)
   {
      printAsterisks(PL_INFO, 0);
      printOutTS(PL_INFO,"* RSMSobol1 generates K levels for each input\n");
      printOutTS(PL_INFO,"* and creates a sample of size M for each level.\n");
      printOutTS(PL_INFO,"* Therefore, the total sample size for this\n");
      printOutTS(PL_INFO,"* analysis is:\n");
      printOutTS(PL_INFO,"*      N = M * K * nInputs\n");
      nSubSamples = 100;
      nLevels = 1000;
      printOutTS(PL_INFO,"* nInputs   = %d\n", nInputs);
      printOutTS(PL_INFO,"* default M = %d\n", nSubSamples);
      printOutTS(PL_INFO,"* default K = %d\n", nLevels);
      printOutTS(PL_INFO,"* Recommendation: K >> M\n");
      printOutTS(PL_INFO,"* Please select different M and K.\n");
      printOutTS(PL_INFO,"* NOTE: large M and K may take very long time.\n");
      printEquals(PL_INFO, 0);
      sprintf(pString,"Enter M (suggestion: 100 - 10000) : ");
      nSubSamples = getInt(100, 50000, pString);
      sprintf(pString,"Enter K (suggestion: 1000 - 10000) : ");
      nLevels = getInt(1000, 50000, pString);
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
      nSubSamples = 100;
      nLevels = 1000;
      if (psConfig_ != NULL)
      {
         cString = psConfig_->getParameter("RSMSobol1_nsubsamples");
         if (cString != NULL)
         {
            sscanf(cString, "%s %s %d", winput1, winput2, &nSubSamples);
            if (nSubSamples < 100)
            {
               printOutTS(PL_INFO,
                    "RSMSobol1 INFO: nSubSamples should be >= 100.\n");
               nSubSamples = 100;
            }
            else
            {
               printOutTS(PL_INFO,
                    "RSMSobol1 INFO: nSubSamples = %d (config).\n",
                    nSubSamples);
            }
         }
         cString = psConfig_->getParameter("RSMSobol1_nlevels");
         if (cString != NULL)
         {
            sscanf(cString, "%s %s %d", winput1, winput2, &nLevels);
            if (nLevels < 1000)
            {
               printOutTS(PL_INFO,
                    "RSMSobol1 INFO: nLevels should be >= 1000.\n");
               nLevels = 1000;
            }
            else
            {
               printOutTS(PL_INFO,
                    "RSMSobol1 INFO: nLevels = %d (config).\n",nLevels);
            }
         }
      }
      if (isScreenDumpModeOn() && printLevel > 1)
      {
         printOutTS(PL_INFO,"* RSMSobol1 creates one sample of size M \n");
         printOutTS(PL_INFO,"* for each of the K levels of each input.\n");
         printOutTS(PL_INFO,"* Therefore, the total sample size is\n");
         printOutTS(PL_INFO,"* N = M * K * nInputs\n");
         printOutTS(PL_INFO,"* nInputs m = %d\n", nInputs);
         printOutTS(PL_INFO,"* default M = %d\n", nSubSamples);
         printOutTS(PL_INFO,"* default K = %d\n", nLevels);
         printOutTS(PL_INFO,
              "* To make changes, re-run with ana_expert mode on.\n");
         printEquals(PL_INFO, 0);
      }
   }
   if (ntimes > 1)
   {
     if (isScreenDumpModeOn())
     {
       printOutTS(PL_INFO,"RSMSobol1 INFO: number of bootstrapped samples\n");
       printOutTS(PL_INFO,"          greater than 1 is not recommended.\n");
       printOutTS(PL_INFO,"          Use the rssobol1b command instead.\n");
     }
     ntimes = 1;
   }
   if (ntimes == 1) frac = 1.0;

   nSamp = nLevels * nSubSamples;
   if (isScreenDumpModeOn())
   {
      printOutTS(PL_INFO,
        "RSMSobol1 INFO: creating a sample for basic statistics.\n");
      printOutTS(PL_INFO,"                sample size = %d\n", nSamp);
   }

   int    *SS;
   double *XX, *YY;
   XX = new double[nSamp*nInputs];
   YY = new double[nSamp];
   checkAllocate(YY, "YY in RSMSobol1::analyze3");
   if (noPDF == 0)
   {
      pdfman = new PDFManager();
      pdfman->initialize(nInputs,pdfFlags,inputMeans,
                         inputStdevs,*corMatp,NULL,NULL);
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
      checkAllocate(SS, "SS in RSMSobol1::analyze3");
      sampler->getSamples(nSamp, nInputs, 1, XX, YY, SS);
      delete [] SS;
      delete sampler;
   }

   int    *bsFlags = new int[nSamples];
   double *bsX = new double[nSamples*nInputs];
   double *bsY = new double[nSamples];
   double *bsMeans, *bsStds;
   bsMeans = new double[ntimes];
   bsStds  = new double[ntimes];
   count = (int) (frac * nSamples);
   faPtr = genFA(rstype, nInputs, 0, count);
   faPtr->setBounds(xLower, xUpper);
   faPtr->setOutputLevel(0);
   for (nn = 0; nn < ntimes; nn++)
   {
      if (ntimes == 1)
      {
         for (ss = 0; ss < nSamples*nInputs; ss++) bsX[ss] = X[ss];
         for (ss = 0; ss < nSamples; ss++) 
            bsY[ss] = Y2[ss*nOutputs+outputID];
      }
      else
      {
         for (ss = 0; ss < nSamples; ss++) bsFlags[ss] = 0;
         count = 0;
         while (count < frac * nSamples)
         {
            jj = PSUADE_rand() % nSamples;
            if (bsFlags[jj] == 0)
            {
               for (ii = 0; ii < nInputs; ii++)
                  bsX[count*nInputs+ii] = X[jj*nInputs+ii];
               bsY[count] = Y2[jj*nOutputs+outputID];
               bsFlags[jj] = 1;
               count++;
            }
         }
      }
      status = faPtr->initialize(bsX, bsY);
      if (status != 0)
      {
         printf("RSMSobol1 ERROR: in initializing response surface.\n");
         delete [] bsX;
         delete [] bsY;
         delete [] bsFlags;
         delete [] XX;
         delete [] YY;
         if (pdfNull == 1)
         {
            delete [] pdfFlags;
            delete [] inputMeans;
            delete [] inputStdevs;
         }
         return -1;
      }
      faPtr->evaluatePoint(nSamp, XX, YY);
      count = 0;
      if (nConstr > 0)
      {
         for (kk = 0; kk < nSamp; kk++)
         {
            ddata = constrPtr->evaluate(&XX[kk*nInputs],YY[kk],status);
            if (status == 0) 
            {
               YY[kk] = PSUADE_UNDEFINED;
               count++;
            }
         }
      }
      count = nSamp - count;
      if (nConstr > 0)
      {
         printOutTS(PL_INFO,
              "RSMSobol1 INFO: %6.2f percent passes the contraints.\n",
              (double) count * 100.0 /((double) nSamp));
      }
      if (count <= 1)
      {
         printf("RSMSobol1 ERROR: too few samples left after filtering\n");
         delete [] bsX;
         delete [] bsY;
         delete [] bsFlags;
         delete [] XX;
         delete [] YY;
         if (pdfNull == 1)
         {
            delete [] pdfFlags;
            delete [] inputMeans;
            delete [] inputStdevs;
         }
         return -1;
      }
      dmean = 0.0;
      for (kk = 0; kk < nSamp; kk++)
         if (YY[kk] != PSUADE_UNDEFINED) dmean += YY[kk];
      dmean /= (double) count;
      dstds = 0.0;
      for (kk = 0; kk < nSamp; kk++)
         if (YY[kk] != PSUADE_UNDEFINED) dstds += pow(YY[kk]-dmean,2.0);
      dstds /= (double) (count - 1);
      bsMeans[nn] = dmean;
      bsStds[nn]  = sqrt(dstds);
   }   
   dmean = 0.0; 
   for (nn = 0; nn < ntimes; nn++) dmean += bsMeans[nn];
   dmean /= (double) ntimes;
   dstds = 0.0;
   if (ntimes > 1)  
   {
      for (nn = 0; nn < ntimes; nn++) dstds += pow(bsMeans[nn]-dmean,2.0);
      dstds = sqrt(dstds/(double) (ntimes - 1));
   } 
   double smean=0, sstd=0;
   smean = 0.0; 
   for (nn = 0; nn < ntimes; nn++) smean += bsStds[nn];
   smean /= (double) ntimes;
   sstd = 0.0;
   if (ntimes > 1)  
   {
      for (nn = 0; nn < ntimes; nn++) sstd += pow(bsStds[nn]-smean,2.0);
      sstd = sqrt(sstd/(double) (ntimes - 1));
   } 
   if (isScreenDumpModeOn())
   {
      printOutTS(PL_INFO,
        "RSMSobol1: sample mean (std dev of mean) = %10.3e (%10.3e)\n",
        dmean, dstds);
      printOutTS(PL_INFO,
        "RSMSobol1: std dev (std dev of std dev)  = %10.3e (%10.3e)\n",
        smean, sstd);
   }
   if (smean == 0.0) smean = 1.0;
   delete [] XX;
   delete [] YY;
   delete [] bsX;
   delete [] bsY;
   delete [] bsFlags;
   delete [] bsMeans;
   delete [] bsStds;
   //save mean & std
   outputMean_ = dmean;
   outputStd_  = smean;
   //cout << outputMean_ << ", " << outputStd_ << endl;

   int    *pdfFlags2, nSteps=1;
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
   checkAllocate(ZZ, "ZZ in RSMSobol1::analyze3");
   for (ii = 0; ii < nInputs; ii++)
   {
      if (noPDF == 0)
      {
         corMat.setDim(1,1);
         corMat.setEntry(0, 0, corMatp->getEntry(ii,ii));
         pdfman1 = new PDFManager();
         pdfman1->initialize(1,&pdfFlags[ii],&inputMeans[ii],
                             &inputStdevs[ii],corMat,NULL,NULL);
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
                             inputStdevs2,corMat,NULL,NULL);
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
   if (pdfNull == 1)
   {
      delete [] pdfFlags;
      delete [] inputMeans;
      delete [] inputStdevs;
   }
   if (ioPtr == NULL) delete corMatp;

   int       iL, sCnt, totalCnt, *bins, offset;
   double    *vceMaxs, *vceMins, *vceMeds;
   double    *mSamplePts, *vces, *ecvs, *vars, *means, *oneSamplePt;
   YY            = new double[nSubSamples*nLevels];
   bsX           = new double[nSamples*nInputs];
   bsY           = new double[nSamples];
   bsFlags       = new int[nSamples];
   vceMaxs       = new double[nInputs];
   vceMins       = new double[nInputs];
   vceMeds       = new double[nInputs];
   vars          = new double[nLevels*ntimes];
   vces          = new double[nInputs*ntimes];
   ecvs          = new double[nInputs*ntimes];
   bins          = new int[nLevels*ntimes];
   means         = new double[nLevels*ntimes];
   mSamplePts    = new double[nInputs*nSubSamples*nLevels];
   checkAllocate(mSamplePts, "mSamplePts in RSMSobol1::analyze3");

   for (ii = 0; ii < nInputs; ii++)
   {
      if (isScreenDumpModeOn())
         printOutTS(PL_INFO, "RSMSobol1: processing input %d\n",ii+1);

      for (iL = 0; iL < nLevels; iL++)
      {
         offset = iL * nSubSamples * nInputs;
         for (jj = 0; jj < nSubSamples; jj++)
         {
            oneSamplePt = &(samplePtsND[ii*nSubSamples*nInputs+jj*(nInputs-1)]);
            for (kk = 0; kk < ii; kk++)
               mSamplePts[offset+jj*nInputs+kk] = oneSamplePt[kk];
            for (kk = ii+1; kk < nInputs; kk++)
               mSamplePts[offset+jj*nInputs+kk] = oneSamplePt[kk-1];
            mSamplePts[offset+jj*nInputs+ii] = samplePts1D[ii*nLevels+iL];
         }
      }

      for (nn = 0; nn < ntimes; nn++)
      {
         if (ntimes == 1)
         {
            for (ss = 0; ss < nSamples*nInputs; ss++) bsX[ss] = X[ss];
            for (ss = 0; ss < nSamples; ss++) 
               bsY[ss] = Y2[ss*nOutputs+outputID];
         }
         else
         {
            for (ss = 0; ss < nSamples; ss++) bsFlags[ss] = 0;
            count = 0;
            while (count < frac * nSamples)
            {
               jj = PSUADE_rand() % nSamples;
               if (bsFlags[jj] == 0)
               {
                  for (kk = 0; kk < nInputs; kk++)
                     bsX[count*nInputs+kk] = X[jj*nInputs+kk];
                  bsY[count] = Y2[jj*nOutputs+outputID];
                  bsFlags[jj] = 1;
                  count++;
               }
            }
         }
         if (ii == 0 || ntimes > 1) 
            status = faPtr->initialize(bsX, bsY);
         if (isScreenDumpModeOn() && printLevel > 3)
            printOutTS(PL_INFO,"RSMSobol1: function evaluations\n");
         faPtr->evaluatePoint(nSubSamples*nLevels,mSamplePts,YY);
         if (nConstr > 0)
         {
            for (kk = 0; kk < nLevels*nSubSamples; kk++)
            {
               ddata = constrPtr->evaluate(&mSamplePts[kk*nInputs],
                                           YY[kk],status);
               if (status == 0) YY[kk] = PSUADE_UNDEFINED;
            }
         }
   
         if (isScreenDumpModeOn() && printLevel > 3 &&
             (iL % (nLevels/10) == 0))
            printOutTS(PL_INFO, "RSMSobol1: compute mean and std dev\n");
         for (iL = 0; iL < nLevels; iL++)
         {
            means[iL*ntimes+nn] = 0.0;
            sCnt = 0;
            for (jj = 0; jj < nSubSamples; jj++)
            {
               if (YY[iL*nSubSamples+jj] != PSUADE_UNDEFINED)
               {
                  means[iL*ntimes+nn] += YY[iL*nSubSamples+jj];
                  sCnt++;
               }
            }
            bins[iL*ntimes+nn] = sCnt;
            if (sCnt < nSubSamples/10 && printLevel >= 5)
               printOutTS(PL_INFO,"RSMSobol1 WARNING: subsample size = %d\n",
                          sCnt);
            if (sCnt >= 1) means[iL*ntimes+nn] /= (double) sCnt;
            else           means[iL*ntimes+nn] = PSUADE_UNDEFINED;
            vars[iL*ntimes+nn] = 0.0;
            if (sCnt > 1)
            {
               for (jj = 0; jj < nSubSamples; jj++)
                  if (YY[iL*nSubSamples+jj] != PSUADE_UNDEFINED)
                     vars[iL*ntimes+nn] += pow(YY[iL*nSubSamples+jj]-
                                               means[iL*ntimes+nn],2.0);
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
            printOutTS(PL_ERROR, "RSMSobol1 ERROR: no feasible region.\n");
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
               ecvs[ii*ntimes+nn] += vars[iL*ntimes+nn]*bins[iL*ntimes+nn]/
                                     totalCnt;
         }
      }
   }
   delete [] bsX;
   delete [] bsY;
   delete [] bsFlags;
   delete [] mSamplePts;
   delete [] samplePts1D;
   delete [] samplePtsND;
   delete [] YY;
   delete [] bins;
   delete faPtr;
   if (constrPtr != NULL) delete constrPtr;

   double *vceMedu = new double[nInputs];
   checkAllocate(vceMedu, "vceMedu in RSMSobol1::analyze3");
   for (ii = 0; ii < nInputs; ii++)
   {
      vceMaxs[ii] = -PSUADE_UNDEFINED;
      vceMins[ii] =  PSUADE_UNDEFINED;
      vceMeds[ii] =  0.0;
      vceMedu[ii] =  0.0;
      for (nn = 0; nn < ntimes; nn++)
      {
         if (vces[ii*ntimes+nn] > vceMaxs[ii])
            vceMaxs[ii] = vces[ii*ntimes+nn];
         if (vces[ii*ntimes+nn] < vceMins[ii]) 
            vceMins[ii] = vces[ii*ntimes+nn];
         vceMeds[ii] += vces[ii*ntimes+nn];
         vceMedu[ii] += (vces[ii*ntimes+nn]-ecvs[ii*ntimes+nn]/nSubSamples);
      }
      vceMeds[ii] /= (double) ntimes;
      vceMedu[ii] /= (double) ntimes;
      if (smean != 0)
      {
         vceMeds[ii] /= (smean * smean);
         vceMedu[ii] /= (smean * smean);
         vceMaxs[ii] /= (smean * smean);
         vceMins[ii] /= (smean * smean);
      }
      else vceMeds[ii] = vceMedu[ii] = vceMaxs[ii] = vceMins[ii] = 0.0;
   }

   vces_ = new double[nInputs_];
   checkAllocate(vces_, "vces_ in RSMSobol1::analyze3");
   if (isScreenDumpModeOn()) printAsterisks(PL_INFO, 0);
   for (ii = 0; ii < nInputs; ii++)
   {
      if (isScreenDumpModeOn())
      {
         printOutTS(PL_INFO,
           "RSMSobol1: Normalized mean VCE for input %3d = %12.4e",
           ii+1, vceMeds[ii]);
         if (ntimes > 1)
              printOutTS(PL_INFO,",bounds = [%12.4e, %12.4e]\n",vceMins[ii],
                         vceMaxs[ii]);
         else printOutTS(PL_INFO, "\n");
      }
      //save vce
      vces_[ii] = vceMeds[ii];
      //cout << vces_[ii] << endl;
   }
   if (isScreenDumpModeOn() && printLevel >= 2)
   {
      ddata = ddata2 = 0.0;
      for (ii = 0; ii < nInputs; ii++)
      {
         printOutTS(PL_INFO,"Unnormalized VCE for input %3d = %12.4e\n",ii+1,
                    vceMeds[ii]*smean*smean);
         printOutTS(PL_INFO,
              "Unnormalized VCE for input %3d = %12.4e (unbiased)\n",
                    ii+1, vceMedu[ii]*smean*smean);
         ddata  += vceMeds[ii]*smean*smean;
         ddata2 += vceMedu[ii]*smean*smean;
      }
      printOutTS(PL_INFO,"Sum of   biased VCEs = %12.4e\n",ddata);
      printOutTS(PL_INFO,"Sum of unbiased VCEs = %12.4e\n",ddata2);
      printOutTS(PL_INFO,"Total variance       = %12.4e\n",smean * smean);
   }

   pData *pObj = NULL;
   if (ioPtr != NULL)
   {
      pObj = ioPtr->getAuxData();
      if (ntimes == 1)
      {
         pObj->nDbles_ = nInputs;
         pObj->dbleArray_ = new double[nInputs];
         for (ii = 0; ii < nInputs; ii++)
            pObj->dbleArray_[ii] = vceMeds[ii] * smean * smean;
         pObj->dbleData_ = smean * smean;
      }
      else
      {
         pObj->nDbles_ = 3*nInputs;
         pObj->dbleArray_ = new double[nInputs*3];
         for (ii = 0; ii < nInputs; ii++)
            pObj->dbleArray_[ii] = vceMeds[ii] * smean * smean;
         for (ii = 0; ii < nInputs; ii++)
            pObj->dbleArray_[nInputs+ii] = vceMins[ii] * smean * smean;
         for (ii = 0; ii < nInputs; ii++)
            pObj->dbleArray_[2*nInputs+ii] = vceMaxs[ii] * smean * smean;
         pObj->dbleData_ = smean * smean;
      }
   }

   if (isScreenDumpModeOn() && printLevel >= 1)
   {
      for (ii = 0; ii < nInputs; ii++) means[ii] = (double) ii;
      printEquals(PL_INFO, 0);
      printOutTS(PL_INFO,"RSMSobol1: ordered normalized VCE : \n");
      printDashes(PL_INFO, 0);
      sortDbleList2(nInputs, vceMeds, means);
      for (ii = nInputs-1; ii >= 0; ii--)
         printOutTS(PL_INFO,"RSMSobol1: Normalized VCE for input %3d = %e\n",
                (int) means[ii]+1,vceMeds[ii]);
      printAsterisks(PL_INFO, 0);
   }

   delete [] vceMaxs;
   delete [] vceMins;
   delete [] vceMeds;
   delete [] vceMedu;
   delete [] means;
   delete [] vars;
   delete [] vces;
   delete [] ecvs;
   return 0.0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
RSMSobol1Analyzer& RSMSobol1Analyzer::operator=(const RSMSobol1Analyzer &)
{
  printOutTS(PL_ERROR, 
       "RSMSobol1 operator= ERROR: operation not allowed.\n");
  exit(1);
  return (*this);
}

// ************************************************************************
// functions for getting results
// ------------------------------------------------------------------------
int RSMSobol1Analyzer::get_nInputs()
{
  return nInputs_;
}
double RSMSobol1Analyzer::get_outputMean()
{
  return outputMean_;
}
double RSMSobol1Analyzer::get_outputStd()
{
  return outputStd_;
}
double RSMSobol1Analyzer::get_vce(int ind)
{
  if (ind < 0 || ind >= nInputs_)
  {
    printf("RSMSobol1 ERROR: get_vce index error.\n");
    return 0.0;
  }
  if (vces_ == NULL)
  {
    printf("RSMSobol1 ERROR: get_vce has not value.\n");
    return 0.0;
  }
  return vces_[ind];
}

