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
// Functions for the class FORMAnalyzer  
// AUTHOR : CHARLES TONG
// DATE   : 2005
// ************************************************************************
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "FuncApprox.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "FORMAnalyzer.h"
#include "PrintingTS.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
FORMAnalyzer::FORMAnalyzer() : Analyzer()
{
  setName("FORM");
  rsType_ = PSUADE_RS_MARS; /* default is MARS */
  printOutTS(PL_ERROR, "FORMAnalyzer currently not supported.\n");
  exit(1);
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
FORMAnalyzer::~FORMAnalyzer() 
{ 
} 

// ************************************************************************ 
// perform analysis 
// ------------------------------------------------------------------------
double FORMAnalyzer::analyze(aData &adata)
{
  int    sInd, ii, status, nPtsPerDim=64, iter, iOne=1;
  double beta, ddata, gdata, Ldata, oldbeta;
  FuncApprox *fa;

  //**/ ---------------------------------------------------------------
  //**/ extract data
  //**/ ---------------------------------------------------------------
  int printLevel = adata.printLevel_;
  int nInputs    = adata.nInputs_;
  int nOutputs   = adata.nOutputs_;
  int nSamples   = adata.nSamples_;
  double *xLower = adata.iLowerB_;
  double *xUpper = adata.iUpperB_;
  double *X      = adata.sampleInputs_;
  double *Y      = adata.sampleOutputs_;
  int outputID   = adata.outputID_;
  int wgtID      = adata.regWgtID_;
  double thresh  = adata.analysisThreshold_;
  if (adata.inputPDFs_ != NULL)
  {
    printOutTS(PL_WARN, 
         "FORM INFO: some inputs have non-uniform PDFs, but\n");
    printOutTS(PL_WARN, 
         "           they will not be relevant in this analysis.\n");
  }

  //**/ ---------------------------------------------------------------
  //**/ error checking and diagnostics
  //**/ ---------------------------------------------------------------
  if (nInputs <= 0 || nOutputs <= 0 || nSamples <= 0)
  {
    printOutTS(PL_ERROR, "FORMAnalyzer ERROR: invalid arguments.\n");
    printOutTS(PL_ERROR, "    nInputs  = %d\n", nInputs);
    printOutTS(PL_ERROR, "    nOutputs = %d\n", nOutputs);
    printOutTS(PL_ERROR, "    nSamples = %d\n", nSamples);
    return PSUADE_UNDEFINED;
  } 
  if (printLevel > 0)
  {
    printAsterisks(PL_INFO, 0);
    printThisFA(rsType_);
    printEquals(PL_INFO, 0);
  }
   
  //**/ ---------------------------------------------------------------
  //**/ copy the selected output data to a new array 
  //**/ ---------------------------------------------------------------
  psVector vecYY, vecWgts;;
  vecYY.setLength(nSamples);
  for (sInd = 0; sInd < nSamples; sInd++)
    vecYY[sInd] = Y[sInd*nOutputs+outputID];
  if (wgtID >= 0)
  {
    vecWgts.setLength(nSamples);
    for (sInd = 0; sInd < nSamples; sInd++)
      vecWgts[sInd] = Y[sInd*nOutputs+wgtID];
  }

  //**/ ---------------------------------------------------------------
  //**/ generate response surface 
  //**/ ---------------------------------------------------------------
  fa = genFA(rsType_, nInputs, iOne=1, nSamples);
  if (fa == NULL)
  {
    printOutTS(PL_INFO,"FORMAnalyzer ERROR: FuncApprox returns NULL.\n");
    delete fa;
    return 1.0e12;
  }
  fa->setNPtsPerDim(nPtsPerDim);
  fa->setBounds(xLower, xUpper);
  fa->setOutputLevel(0);
  if (wgtID >= 0) fa->loadWeights(nSamples, vecWgts.getDVector());

  status = fa->initialize(X, vecYY.getDVector());
  if (status != 0)
  {
    printOutTS(PL_INFO, 
         "FORMAnalyzer ERROR: FuncApprox returns error.\n");
    delete fa;
    return 1.0e12;
  }

  //**/ ---------------------------------------------------------------
  //**/ skip the transformation for now
  //**/ ---------------------------------------------------------------

  //**/ ---------------------------------------------------------------
  //**/ set initial guess and initialize
  //**/ ---------------------------------------------------------------
  psVector vecCurrZ, vecAlphas;
  vecCurrZ.setLength(nInputs);
  for (ii = 0; ii < nInputs; ii++)
    vecCurrZ[ii] = 0.5 * (xUpper[ii] + xLower[ii]);
  vecAlphas.setLength(nInputs);
  beta = 0.0;
  for (ii = 0; ii < nInputs; ii++) beta += vecCurrZ[ii] * vecCurrZ[ii];
  beta = sqrt(beta);

  //**/ ---------------------------------------------------------------
  //**/ iterate (beta = sqrt(Z'*Z), iter = 1) 
  //**/    compute alpha = partial(F(Z))/partial(Z)
  //**/    compute F(Z) 
  //**/    compute Z = - alpha (beta + F(Z)/L)
  //**/    update beta
  //**/    check whether beta has stabilized
  //**/ ---------------------------------------------------------------
  iter = 1;
  while (iter < 1000)
  {
    gdata = fa->evaluatePoint(vecCurrZ.getDVector()) - thresh;
    for (ii = 0; ii < nInputs; ii++)
    {
      vecCurrZ[ii] += 1.0e-5;
      ddata = fa->evaluatePoint(vecCurrZ.getDVector()) - thresh;
      vecCurrZ[ii] -= 1.0e-5;
      vecAlphas[ii] = (ddata - gdata) * 1.0e5;
    }
    Ldata = 0.0;
    for (ii = 0; ii < nInputs; ii++) 
      Ldata += vecAlphas[ii] * vecAlphas[ii];
    Ldata = sqrt(Ldata);
    for (ii = 0; ii < nInputs; ii++) vecAlphas[ii] /= Ldata;
    if (Ldata < 1.0e-12) break;
    for (ii = 0; ii < nInputs; ii++) 
      vecCurrZ[ii] = - vecAlphas[ii] * (beta + gdata / Ldata);
    oldbeta = beta;
    beta = 0.0;
    for (ii = 0; ii < nInputs; ii++) beta += vecCurrZ[ii] * vecCurrZ[ii];
    beta = sqrt(beta);
    if (printLevel > 1)
      printOutTS(PL_INFO, 
           "FORM analyze: beta at iter %4d = %12.4e(%12.4e)\n",iter,
           beta, oldbeta);
    if (PABS((beta-oldbeta)/(beta+oldbeta)) < 1.0e-4) break;
    iter++;
  }
  if (printLevel > 1)
  {
    printOutTS(PL_INFO, "FORMAnalyzer: optimal point = \n");
    for (ii = 0; ii < nInputs; ii++) 
      printOutTS(PL_INFO, "\t X[%3d] = %12.4e\n",ii+1,vecCurrZ[ii]);
    printOutTS(PL_INFO, "optimal distance = %12.4e\n",beta);
    printOutTS(PL_INFO, "number of iterations = %d\n",iter);
  }

  //**/ ---------------------------------------------------------------
  //**/ clean up 
  //**/ ---------------------------------------------------------------
  delete fa;
  return 0.0;
}

// ************************************************************************ 
// set parameters
// ------------------------------------------------------------------------
int FORMAnalyzer::setParams(int argc, char **argv)
{
  char  *request = (char *) argv[0];
  Analyzer::setParams(argc, argv);
  if (!strcmp(request, "rstype"))
  {
    if (argc != 2) printOutTS(PL_INFO,"FORMAnalyzer WARNING: setParams.\n");
    rsType_ = *(int *) argv[1];
    if (rsType_ < 0 || rsType_ >= PSUADE_NUM_RS) rsType_ = 0;
  }
  else
  {
    printOutTS(PL_ERROR, "FORMAnalyzer ERROR: setParams not valid.\n");
    exit(1);
  }
  return 0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
FORMAnalyzer& FORMAnalyzer::operator=(const FORMAnalyzer &)
{
  printOutTS(PL_ERROR, 
       "FORMAnalyzer operator= ERROR: operation not allowed.\n");
  exit(1);
   return (*this);
}

