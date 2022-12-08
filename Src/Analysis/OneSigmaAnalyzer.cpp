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
// Functions for the class OneSigmaAnalyzer  
//**/ This analyzer uses GP, MARSB, or Kriging to assess the standard 
//**/ deviation of the approximations at the original set of sample points. 
//**/ This analyzer should be modified to use factorial grid points as test 
//**/ points.
// AUTHOR : CHARLES TONG
// DATE   : 2005
// ************************************************************************

// ************************************************************************
// ************************************************************************
// Author : Charles Tong
// Date   : 2005
// ************************************************************************
// ************************************************************************
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "FuncApprox.h"
#include "OneSigmaAnalyzer.h"
#include "PrintingTS.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
OneSigmaAnalyzer::OneSigmaAnalyzer() : Analyzer()
{
  setName("OneSigma");
  //**/ default is MARSBag
  rsType_ = PSUADE_RS_MARSB; 
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
OneSigmaAnalyzer::~OneSigmaAnalyzer() 
{ 
} 

// ************************************************************************ 
// perform analysis 
// ------------------------------------------------------------------------
double OneSigmaAnalyzer::analyze(aData &adata)
{
  int    ss, ii, nPtsPerDim=64, count, nPts, n1d, itmp, jtmp, iOne=1;
  double Ymax, step, stdev, maxStdev, ddata, avgStdev;
  psVector vecWT, vec1Sam;

  //**/ ---------------------------------------------------------------
  //**/ extract data
  //**/ ---------------------------------------------------------------
  int printLevel = adata.printLevel_;
  int nInputs    = adata.nInputs_;
  int nOutputs   = adata.nOutputs_;
  int nSamples   = adata.nSamples_;
  double *iLowerB = adata.iLowerB_;
  double *iUpperB = adata.iUpperB_;
  double *X = adata.sampleInputs_;
  double *Y = adata.sampleOutputs_;
  int    *S = adata.sampleStates_;
  int outputID = adata.outputID_;
  int wgtID    = adata.regWgtID_;
  if (adata.inputPDFs_ != NULL)
  {
    count = 0;
    for (ii = 0; ii < nInputs; ii++) count += adata.inputPDFs_[ii];
    if (count > 0)
    {
      printOutTS(PL_INFO, 
         "OneSigmaTest INFO: some inputs have non-uniform PDFs,\n");
      printOutTS(PL_INFO, 
         "        but they are not relevant in this analysis.\n");
    }
  }
  int status = 0;
  for (ss = 0; ss < nSamples; ss++)
    if (Y[nOutputs*ss+outputID] > 0.9*PSUADE_UNDEFINED) status = 1;
  if (status == 1)
  {
    printOutTS(PL_ERROR, 
          "OneSigmaTest ERROR: Some outputs are undefined. Prune\n");
    printOutTS(PL_ERROR, 
          "                    the undefined sample points first.\n");
    return PSUADE_UNDEFINED;
  }

  //**/ ---------------------------------------------------------------
  //**/ not implemented for large dimension 
  //**/ ---------------------------------------------------------------
  if (nInputs > 12)
  {
    printOutTS(PL_ERROR, 
         "OneSigmaTest ERROR: does not support nInputs > 12.\n");
    return PSUADE_UNDEFINED;
  }
  if (nInputs == 1 ) n1d = nSamples;
  if (nInputs == 2 ) n1d = 512;
  if (nInputs == 3 ) n1d = 64;
  if (nInputs == 4 ) n1d = 23;
  if (nInputs == 5 ) n1d = 12;
  if (nInputs == 6 ) n1d = 8;
  if (nInputs == 7 ) n1d = 6;
  if (nInputs == 8 ) n1d = 5;
  if (nInputs == 9 ) n1d = 4;
  if (nInputs == 10) n1d = 3;
  if (nInputs >= 11) n1d = 3;

  //**/ ---------------------------------------------------------------
  //**/ error checking and diagnostics
  //**/ ---------------------------------------------------------------
  if (nInputs <= 0 || nOutputs <= 0 || nSamples <= 0)
  {
    printOutTS(PL_ERROR, "OneSigmaTest ERROR: invalid arguments.\n");
    printOutTS(PL_ERROR, "    nInputs  = %d\n", nInputs);
    printOutTS(PL_ERROR, "    nOutputs = %d\n", nOutputs);
    printOutTS(PL_ERROR, "    nSamples = %d\n", nSamples);
    return PSUADE_UNDEFINED;
  } 
  if (printLevel > 0)
  {
    printAsterisks(PL_INFO, 0);
    printOutTS(PL_INFO, "               OneSigma Analysis\n");
    printDashes(PL_INFO, 0);
    printOutTS(PL_INFO, 
       " This analysis computes standard deviations when interpolated\n");
    printOutTS(PL_INFO, " at the sample points.\n");
    printEquals(PL_INFO, 0);
    if (rsType_ == PSUADE_RS_GP3) 
       printOutTS(PL_INFO, "OneSigmaTest: Gaussian process (local).\n");
    if (rsType_ == PSUADE_RS_GP1) 
       printOutTS(PL_INFO, "OneSigmaTest: Gaussian process model.\n");
    if (rsType_ == PSUADE_RS_KR)  
       printOutTS(PL_INFO, "OneSigmaTest: Kriging model.\n");
    if (rsType_ == PSUADE_RS_MARSB) 
       printOutTS(PL_INFO, "OneSigmaTest: MARS with bagging.\n");
    printEquals(PL_INFO, 0);
  }
   
  //**/ ---------------------------------------------------------------
  //**/ see how many data points are good to go
  //**/ ---------------------------------------------------------------
  int nGood = 0;
  for (ss = 0; ss < nSamples; ss++) if (S[ss] == 1) nGood++;
  if (nGood <= 0) return 0;

  // ---------------------------------------------------------------
  // copy the selected sample data to a new array and find max
  // ---------------------------------------------------------------
  psVector vecXT, vecYT;
  vecXT.setLength(nInputs*nGood);
  vecYT.setLength(nGood);
  nGood = 0;
  for (ss = 0; ss < nSamples; ss++)
  {
    if (S[ss] == 1)
    {
      for (ii = 0; ii < nInputs; ii++)
        vecXT[nGood*nInputs+ii] = X[ss*nInputs+ii];
      vecYT[nGood++] = Y[ss*nOutputs+outputID];
    }
  }
  Ymax = 0.0;
  for (ss = 0; ss < nGood; ss++)
    if (PABS(vecYT[ss]) > Ymax) Ymax = PABS(vecYT[ss]);
  printOutTS(PL_INFO,
       "OneSigmaTest: nSamples = %d, of which %d are valid.\n",
       nSamples, nGood);
  printOutTS(PL_INFO,"OneSigmaTest: Output Maximum Norm = %e\n",Ymax);

  //**/ ---------------------------------------------------------------
  //**/ set up the function approximators
  //**/ ---------------------------------------------------------------
  FuncApprox *fa = genFA(rsType_, nInputs, iOne, nGood);
  if (fa == NULL)
  {
    printOutTS(PL_ERROR, 
       "OneSigmaTest ERROR: failed to create function approximator.\n");
    return 1.0e12;
  }
  fa->setNPtsPerDim(nPtsPerDim);
  fa->setBounds(iLowerB, iUpperB);
  fa->setOutputLevel(printLevel);
  if (wgtID >= 0 && wgtID < nOutputs)
  {
    vecWT.setLength(nGood);
    nGood = 0;
    for (ss = 0; ss < nSamples; ss++)
      if (S[ss] == 1) vecWT[nGood++] = Y[ss*nOutputs+wgtID];
    fa->loadWeights(nGood, vecWT.getDVector());
  }

  //**/ ---------------------------------------------------------------
  //**/ train the function approximators
  //**/ ---------------------------------------------------------------
  status = fa->initialize(vecXT.getDVector(), vecYT.getDVector());
  if (status != 0)
  {
    delete fa;
    return 1.0e12;
  }

  //**/ ---------------------------------------------------------------
  //**/ check RS data against the original data
  //**/ ---------------------------------------------------------------
  printOutTS(PL_INFO, 
     "OneSigmaTest: Step 1 - assess uncertainties of evaluations at\n");
  printOutTS(PL_INFO, "                       the original data set.\n");
  adata.sampleErrors_ = new double[nSamples];
  for (ss = 0; ss < nSamples; ss++) adata.sampleErrors_[ss] = 0.0;

  maxStdev = avgStdev = 0.0;
  nPts = 0;
  for (ss = 0; ss < nSamples; ss++) if (S[ss] != 1) nPts++;
  count = 0;
  for (ss = 0; ss < nSamples; ss++)
  {
    if (S[ss] != 1)
    {
      ddata = fa->evaluatePointFuzzy(&X[ss*nInputs], stdev);
      printOutTS(PL_INFO, "Sample Point %5d (of %5d): \n",count+1,nPts);
      for (ii = 0; ii < nInputs; ii++)
        printOutTS(PL_INFO, " X %4d   %12.4e\n", ii+1, X[ss*nInputs+ii]);
      printOutTS(PL_INFO, 
           "    Sample Mean/Error = %12.4e %12.4e\n", ddata, stdev);

      adata.sampleErrors_[ss] = stdev;
      if (PABS(stdev) > maxStdev) maxStdev = PABS(stdev);
      avgStdev += PABS(stdev);
      count++;
    }
  }
  printOutTS(PL_INFO,"OneSigmaTest: max std dev (%d) = %9.2e\n",
             count,maxStdev);
  if (count > 0)
  {
    avgStdev /= (double) count;
    printOutTS(PL_INFO, 
        "OneSigmaTest: (1) data - avg std dev = %9.2e (%d points)\n",
        avgStdev, count);
    printOutTS(PL_INFO, 
        "OneSigmaTest: (1) data - max std dev = %9.2e\n", maxStdev);
  }

  // ----------------------------------------------------------------
  // check RS data against grid data
  // ----------------------------------------------------------------
  printOutTS(PL_INFO, 
     "OneSigmaTest: Step 2 - assess uncertainties of evaluations at\n");
  printOutTS(PL_INFO, 
     "                   fixed lattice points in the parameter space.\n");
  nPts = 1;
  for (ii = 0; ii < nInputs; ii++) nPts *= n1d;
  maxStdev = avgStdev = 0.0;
  count = 0;

  vec1Sam.setLength(nInputs);
  for (ss = 0; ss < nPts; ss++)
  {
    itmp = ss;
    if (ss % (nPts/10) == 0)
      printOutTS(PL_INFO, 
           "OneSigmaTest: processing grid sample %d (out of %d)\n",
           ss+1, nPts);
    for (ii = 0; ii < nInputs; ii++)
    {
      jtmp = itmp % n1d;
      itmp = itmp / n1d;
      step = (iUpperB[ii] - iLowerB[ii]) / (n1d - 1);
      vec1Sam[ii] = jtmp * step + iLowerB[ii];
      ddata = fa->evaluatePointFuzzy(vec1Sam.getDVector(), stdev);
      if (stdev > maxStdev) maxStdev = stdev;
      avgStdev += stdev;
      count++;
    }
  }
  if (count > 0)
  {
    avgStdev /= (double) count;
    printOutTS(PL_INFO, 
         "OneSigmaTest: (2) grid - avg std dev = %9.2e (%d points)\n",
         avgStdev, count);
    printOutTS(PL_INFO, "OneSigmaTest: (2) grid - max std dev = %9.2e\n",
         maxStdev);
  }
  printAsterisks(PL_INFO, 0);

  //**/ ---------------------------------------------------------------
  //**/ clean up 
  //**/ ---------------------------------------------------------------
  delete fa;
  return maxStdev;
}

// ************************************************************************ 
// set parameters
// ------------------------------------------------------------------------
int OneSigmaAnalyzer::setParams(int argc, char **argv)
{
   char  *request = (char *) argv[0];
   if (!strcmp(request, "rstype"))
   {
      if (argc != 2) 
         printOutTS(PL_WARN, "OneSigmaTest WARNING: setParams.\n");
      rsType_ = *(int *) argv[1];
      if (rsType_ != 5 && rsType_ != 7) rsType_ = 7;
   }
   else
   {
      printOutTS(PL_ERROR,"OneSigmaTest ERROR: setParams - not valid.\n");
      exit(1);
   }
   return 0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
OneSigmaAnalyzer& OneSigmaAnalyzer::operator=(const OneSigmaAnalyzer &)
{
   printOutTS(PL_ERROR, 
        "OneSigmaTest operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

