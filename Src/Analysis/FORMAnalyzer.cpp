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
#include "FuncApprox/FuncApprox.h"
#include "Util/sysdef.h"
#include "Util/PsuadeUtil.h"
#include "FORMAnalyzer.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
FORMAnalyzer::FORMAnalyzer() : Analyzer()
{
   setName("FORM");
   rsType_ = PSUADE_RS_MARS; /* default is MARS */
   printf("FORMAnalyzer currently not supported.\n");
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
   int        nInputs, nOutputs, nSamples, outputID, sInd, ii, status;
   int        nPtsPerDim=64, length, wgtID, iter, printLevel;
   double     *X, *Y, *xLower, *xUpper, *YY, *wgts=NULL, *currZ, beta;
   double     ddata, gdata, Ldata, oldbeta, *alphas, thresh;
   FuncApprox *fa;

   printLevel = adata.printLevel_;
   nInputs    = adata.nInputs_;
   nOutputs   = adata.nOutputs_;
   nSamples   = adata.nSamples_;
   xLower     = adata.iLowerB_;
   xUpper     = adata.iUpperB_;
   X          = adata.sampleInputs_;
   Y          = adata.sampleOutputs_;
   outputID   = adata.outputID_;
   wgtID      = adata.regWgtID_;
   thresh     = adata.analysisThreshold_;
   if (adata.inputPDFs_ != NULL)
   {
      printf("FORM INFO: some inputs have non-uniform PDFs, but\n");
      printf("           they will not be relevant in this analysis.\n");
   }

   if (nInputs <= 0 || nOutputs <= 0 || nSamples <= 0)
   {
      printf("FORMAnalyzer ERROR: invalid arguments.\n");
      printf("    nInputs  = %d\n", nInputs);
      printf("    nOutputs = %d\n", nOutputs);
      printf("    nSamples = %d\n", nSamples);
      return PSUADE_UNDEFINED;
   } 
   if (printLevel > 0)
   {
      printAsterisks(0);
      if (rsType_ == 0) printf("FORM: MARS model.\n");
      if (rsType_ == 1) printf("FORM: linear regression model.\n");
      if (rsType_ == 2) printf("FORM: quadratic regression model.\n");
      if (rsType_ == 3) printf("FORM: cubic regression model.\n");
      if (rsType_ == 4) printf("FORM: quartic regression model.\n");
      if (rsType_ == 5) printf("FORM: artificial neural network model.\n");
      if (rsType_ == 6) printf("FORM: user-specified regression model.\n");
      if (rsType_ == 7) printf("FORM: Gaussian process model.\n");
      printEquals(0);
   }
   
   YY = new double[nSamples];
   for (sInd = 0; sInd < nSamples; sInd++)
      YY[sInd] = Y[sInd*nOutputs+outputID];
   if (wgtID >= 0)
   {
      wgts = new double[nSamples];
      for (sInd = 0; sInd < nSamples; sInd++)
         wgts[sInd] = Y[sInd*nOutputs+wgtID];
   }

   fa = genFA(rsType_, nInputs, nSamples);
   if (fa == NULL)
   {
      printf("FORMAnalyzer ERROR: FuncApprox returns NULL.\n");
      delete [] YY;
      delete fa;
      if (wgts != NULL) delete [] wgts;
      return 1.0e12;
   }
   fa->setNPtsPerDim(nPtsPerDim);
   fa->setBounds(xLower, xUpper);
   fa->setOutputLevel(0);
   if (wgtID >= 0)
   {
      fa->loadWeights(nSamples, wgts);
      delete [] wgts;
   }
   length = -999;
   status = fa->genNDGridData(X, YY, &length, NULL, NULL);
   if (status != 0)
   {
      printf("FORMAnalyzer ERROR: FuncApprox returns error.\n");
      delete [] YY;
      delete fa;
      return 1.0e12;
   }


   currZ = new double[nInputs];
   for (ii = 0; ii < nInputs; ii++)
      currZ[ii] = 0.5 * (xUpper[ii] + xLower[ii]);
   alphas = new double[nInputs];
   beta = 0.0;
   for (ii = 0; ii < nInputs; ii++) beta += currZ[ii] * currZ[ii];
   beta = sqrt(beta);

   iter = 1;
   while (iter < 1000)
   {
      gdata = fa->evaluatePoint(currZ) - thresh;
      for (ii = 0; ii < nInputs; ii++)
      {
         currZ[ii] += 1.0e-5;
         ddata = fa->evaluatePoint(currZ) - thresh;
         currZ[ii] -= 1.0e-5;
         alphas[ii] = (ddata - gdata) * 1.0e5;
      }
      Ldata = 0.0;
      for (ii = 0; ii < nInputs; ii++) Ldata += alphas[ii] * alphas[ii];
      Ldata = sqrt(Ldata);
      for (ii = 0; ii < nInputs; ii++) alphas[ii] /= Ldata;
      if (Ldata < 1.0e-12) break;
      for (ii = 0; ii < nInputs; ii++) 
         currZ[ii] = - alphas[ii] * (beta + gdata / Ldata);
      oldbeta = beta;
      beta = 0.0;
      for (ii = 0; ii < nInputs; ii++) beta += currZ[ii] * currZ[ii];
      beta = sqrt(beta);
      if (printLevel > 1)
         printf("FORM analyze: beta at iter %4d = %12.4e(%12.4e)\n",iter,
                beta, oldbeta);
      if (PABS((beta-oldbeta)/(beta+oldbeta)) < 1.0e-4) break;
      iter++;
   }
   if (printLevel > 1)
   {
      printf("FORMAnalyzer: optimal point = \n");
      for (ii = 0; ii < nInputs; ii++) 
         printf("\t X[%3d] = %12.4e\n",ii+1,currZ[ii]);
      printf("optimal distance = %12.4e\n",beta);
      printf("number of iterations = %d\n",iter);
   }

   delete fa;
   delete [] currZ;
   delete [] alphas;
   delete [] YY;
   return 0.0;
}

// ************************************************************************ 
// set parameters
// ------------------------------------------------------------------------
int FORMAnalyzer::setParams(int argc, char **argv)
{
   char  *request = (char *) argv[0];
   if (!strcmp(request, "rstype"))
   {
      if (argc != 2) printf("FORMAnalyzer WARNING: setParams.\n");
      rsType_ = *(int *) argv[1];
      if (rsType_ < 0 || rsType_ > 7) rsType_ = 0;
   }
   else
   {
      printf("FORMAnalyzer ERROR: setParams not valid.\n");
      exit(1);
   }
   return 0;
}

