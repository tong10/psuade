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
// Function for the class GradStatAnalyzer  
// Functions for the class FASTAnalyzer
// AUTHOR : CHARLES TONG
// DATE   : 2003
// ************************************************************************
//**/ Analyze how good the previous set of sample is to approximate the
//**/ current set of sample via first order approximations
//**/ (gradient information need to be given)
//**/ - the nInputs entries after the outputID entry of Y should contain
//**/   the gradient information
//**/ - The original outputs and their gradient information are used to
//**/   compute mean and standard deviation
// ************************************************************************
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "GradStatAnalyzer.h"
#include "Sampling.h"
#include "PDFBase.h"
#include "PDFNormal.h"
#include "PsuadeUtil.h"
#include "sysdef.h"
#include "PrintingTS.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
GradStatAnalyzer::GradStatAnalyzer() : Analyzer()
{
  setName("GRADSTAT");
  sampler_     = NULL;
  inputPDFs_   = NULL;
  inputMeans_  = NULL;
  inputStdevs_ = NULL;
  threshold_   = 0.0;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
GradStatAnalyzer::~GradStatAnalyzer()
{
  if (sampler_ != NULL) delete sampler_;
}

// ************************************************************************
// perform analysis
// ------------------------------------------------------------------------
double GradStatAnalyzer::analyze(aData &adata)
{
  int     whichOutput, totalNSamples, prevNSamples, ss, ss2, inputID;
  int     minIndex, nLHSamples, iOne=1, nReps, randomize;
  double  curThresh, distance, minDistance, interpolatedOutput, xtemp;
  double  totalVol, ytemp;
  PDFBase **normalPDFs;

  //**/ ---------------------------------------------------------------
  //**/ extract data
  //**/ ---------------------------------------------------------------
  int nInputs  = adata.nInputs_;
  int nOutputs = adata.nOutputs_;
  int nSamples = adata.nSamples_;
  double *inputLBounds = adata.iLowerB_;
  double *inputUBounds = adata.iUpperB_;
  double *X = adata.sampleInputs_;
  double *Y = adata.sampleOutputs_;
  int outputID = adata.outputID_;
  int nLevels  = adata.currRefineLevel_;
  int *levelSeps = adata.refineSeparators_;
  threshold_ = adata.analysisThreshold_;
  inputPDFs_ = adata.inputPDFs_;
  inputMeans_  = adata.inputMeans_;
  inputStdevs_ = adata.inputStdevs_;
  if (inputPDFs_ != NULL)
  {
    printOutTS(PL_INFO,
         "GradStatAnalyzer INFO: non-uniform distributions\n");
    printOutTS(PL_INFO,
         "    detected. Analysis will use these distributions.\n");
  }

  //**/ ---------------------------------------------------------------
  //**/ error checking
  //**/ ---------------------------------------------------------------
  if (nInputs <= 0 || nOutputs <= 0 || nSamples <= 0)
  {
    printOutTS(PL_ERROR,"GradStatAnalyzer ERROR: invalid arguments.\n");
    printOutTS(PL_ERROR,"    nInputs  = %d\n", nInputs);
    printOutTS(PL_ERROR,"    nOutputs = %d\n", nOutputs);
    printOutTS(PL_ERROR,"    nSamples = %d\n", nSamples);
    return PSUADE_UNDEFINED;
  } 
  whichOutput = outputID;
  if (whichOutput >= nOutputs || whichOutput < 0) whichOutput = 0;
  if ((whichOutput + nInputs) > nOutputs)
  {
    printOutTS(PL_ERROR,"GradStatAnalyzer ERROR: not enough outputs.\n");
    printOutTS(PL_ERROR,
         "                        nOutpus should be nInputs+1.\n");
    return PSUADE_UNDEFINED;
  } 
   
  //**/ ---------------------------------------------------------------
  //**/ first find the mean of the current set of samples
  //**/ ---------------------------------------------------------------
  printAsterisks(PL_INFO, 0);
  printOutTS(PL_INFO,
     "            Gradient-based Uncertainty Analysis\n");
  printDashes(PL_INFO, 0);
  printOutTS(PL_INFO,
     " This analysis is for models that produce derviative\n");
  printOutTS(PL_INFO,
     " information in addition to typical outputs. Thus, the\n");
  printOutTS(PL_INFO,
     " number of outputs is expected to be nInputs+1. The\n");
  printOutTS(PL_INFO,
     " outputs and derivatives are used to construct a respone\n");
  printOutTS(PL_INFO,
     " surface which will be probed to compute basic statistics\n");
  printOutTS(PL_INFO," based on the given distributions.\n");
  printEquals(PL_INFO, 0);
  printOutTS(PL_INFO,
     "GradStatAnalyzer: non-Gradient-based analysis\n");
  psVector VecY;
  VecY.setLength(nSamples);
  for (ss = 0; ss < nSamples; ss++) 
    VecY[ss] = Y[nOutputs*ss+whichOutput];
  computeMeanVariance(VecY);
  printDashes(PL_INFO, 0);

  //**/ ---------------------------------------------------------------
  //**/ perform analysis
  //**/ ---------------------------------------------------------------
  totalNSamples = levelSeps[nLevels];
  if (nLevels > 0)
  {
    prevNSamples  = levelSeps[nLevels-1];
    curThresh = 0.0;
    totalVol = 0.0;
    for (ss = prevNSamples; ss < totalNSamples; ss++)
    {
      minDistance = 1.0e30;
      /* choose the point with the least of the max delta y */
      for (ss2 = 0; ss2 < prevNSamples; ss2++)
      {
        distance = 0.0;
        minIndex = 0;
        for (inputID = 0; inputID < nInputs; inputID++)
        {
          xtemp = X[ss*nInputs+inputID]-X[ss2*nInputs+inputID];
          xtemp *= Y[ss2*nOutputs+whichOutput+inputID+1];
          if (PABS(xtemp) > distance) distance = PABS(xtemp);
        }
        if (distance < minDistance) 
        {
          minIndex = ss2;
          minDistance = distance;
        }
      }
      interpolatedOutput = Y[minIndex*nOutputs+whichOutput];
      for (inputID = 0; inputID < nInputs; inputID++)
      {
        xtemp = X[ss*nInputs+inputID]-X[minIndex*nInputs+inputID];
        interpolatedOutput += 
            xtemp*Y[minIndex*nOutputs+whichOutput+inputID+1];
      }
      ytemp = Y[ss*nOutputs+whichOutput] - interpolatedOutput;

      curThresh += PABS(ytemp);
      totalVol += PABS(Y[ss*nOutputs+whichOutput]);
    }
    curThresh /= (double) (totalNSamples - prevNSamples);
    totalVol /= (double) (totalNSamples - prevNSamples);
    curThresh /= totalVol;
    printOutTS(PL_INFO,  "GradStatAnalyzer: convergence check %e < %e ?\n",
           curThresh, threshold_);
  }
  else curThresh = 1.0;

  //**/ ---------------------------------------------------------------
  //**/ generate LHS samples and perform analysis
  //**/ ---------------------------------------------------------------
  randomize = 1;
  nReps = 1;
  sampler_ = SamplingCreateFromID(PSUADE_SAMP_LHS);
  sampler_->setInputBounds(nInputs, inputLBounds, inputUBounds);
  sampler_->setOutputParams(nOutputs);
  sampler_->setSamplingParams(nSamples, nReps, randomize);
  sampler_->initialize(0);

  //**/ refine local version of sampler
  nLHSamples = nSamples;
  while (nLHSamples < 8000)
  {
    sampler_->refine(2, 1, 1.0, 0, NULL);
    nLHSamples = sampler_->getNumSamples(); 
  }
  if (nLHSamples == nSamples)
  {
    printOutTS(PL_INFO,  "GradStatAnalyzer: nSamples >= 8000 - stop.\n");
    curThresh = 0.0;
  }

  psVector  VecLHSInps, VecLHSOuts;
  psIVector VecLHSStas;
  VecLHSInps.setLength(nSamples*nInputs);
  VecLHSOuts.setLength(nSamples*nOutputs);
  VecLHSStas.setLength(nSamples);
  sampler_->getSamples(nSamples,nInputs,nOutputs,VecLHSInps.getDVector(),
                       VecLHSOuts.getDVector(), VecLHSStas.getIVector());

  //**/ transform if not uniform distribution
  psVector VecT1, VecT2;
  VecT1.setLength(nLHSamples);
  VecT2.setLength(nLHSamples);
  normalPDFs = new PDFBase*[nInputs];
  for (inputID = 0; inputID < nInputs; inputID++)
  {
    if (inputPDFs_ != NULL && inputPDFs_[inputID] == 1)
      normalPDFs[inputID] = (PDFBase *) new PDFNormal(inputMeans_[inputID],
                                            inputStdevs_[inputID]);
    else normalPDFs[inputID] = NULL;
  }
  for (inputID = 0; inputID < nInputs; inputID++)
  {
    if (normalPDFs[inputID] != NULL)
    {
      for (ss = 0; ss < nLHSamples; ss++)
        VecT1[ss] = VecLHSInps[ss*nInputs+inputID];
      normalPDFs[inputID]->invCDF(nLHSamples,VecT1.getDVector(),
                                  VecT2.getDVector());
      for (ss = 0; ss < nLHSamples; ss++)
      {
        VecLHSInps[ss*nInputs+inputID] = VecT2[ss];
        if (VecT2[ss] < inputLBounds[inputID])
          VecLHSInps[ss*nInputs+inputID] = inputLBounds[inputID];
        if (VecT2[ss] > inputUBounds[inputID])
          VecLHSInps[ss*nInputs+inputID] = inputUBounds[inputID];
      }
    }
  }
  for (inputID = 0; inputID < nInputs; inputID++)
    if (normalPDFs[inputID] != NULL) delete normalPDFs[inputID];
  delete [] normalPDFs;

  //**/ finally analyze 
  for (ss = 0; ss < nLHSamples; ss++)
  {
    minDistance = 1.0e30;
    for (ss2 = 0; ss2 < totalNSamples; ss2++)
    {
      distance = 0.0;
      for (inputID = 0; inputID < nInputs; inputID++)
      {
        xtemp = VecLHSInps[ss*nInputs+inputID] -
                 X[ss2*nInputs+inputID];
        xtemp *= Y[ss2*nOutputs+whichOutput+inputID+1];
        if (PABS(xtemp) > distance) distance = PABS(xtemp);
      }
      if (distance < minDistance) 
      {
        minIndex = ss2;
        minDistance = distance;
      }
    }
    interpolatedOutput = Y[minIndex*nOutputs+whichOutput];
    for (inputID = 0; inputID < nInputs; inputID++)
    {
      xtemp = VecLHSInps[ss*nInputs+inputID] - 
                X[minIndex*nInputs+inputID];
      interpolatedOutput += xtemp*Y[minIndex*nOutputs+whichOutput+inputID+1];
    }
    VecLHSOuts[ss] = interpolatedOutput;
  }
  printOutTS(PL_INFO,  "GradStatAnalyzer: Gradient-based analysis\n");
  computeMeanVariance(VecLHSOuts);
  printAsterisks(PL_INFO, 0);
  return curThresh;
}

// *************************************************************************
// Compute agglomerated mean and variance
// -------------------------------------------------------------------------
int GradStatAnalyzer::computeMeanVariance(psVector VecY) 
{
  int    ss, nSamples;
  double mean=0, variance=0;

  nSamples = VecY.length();
  for (ss = 0; ss < nSamples; ss++) mean += VecY[ss];
  mean /= (double) nSamples;
  for (ss = 0; ss < nSamples; ss++) 
    variance += ((VecY[ss] - mean) * (VecY[ss] - mean));
  variance /= (double) (nSamples - 1);
  printOutTS(PL_INFO, "GradStatAnalyzer: mean     = %24.16e\n", mean);
  printOutTS(PL_INFO,
       "GradStatAnalyzer: std dev  = %24.16e\n", sqrt(variance));
  return 0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
GradStatAnalyzer& GradStatAnalyzer::operator=(const GradStatAnalyzer &)
{
  printOutTS(PL_ERROR,
       "GradStatAnalyzer operator= ERROR: operation not allowed.\n");
  exit(1);
  return (*this);
}

