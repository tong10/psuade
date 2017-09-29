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
// Functions for the class IntegrationAnalyzer  
// AUTHOR : CHARLES TONG
// DATE   : 2005
// ************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include "IntegrationAnalyzer.h"
#include "Util/PsuadeUtil.h"
#include "Util/sysdef.h"
#define PABS(x) (((x) > 0) ? x : -(x))

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
IntegrationAnalyzer::IntegrationAnalyzer() : Analyzer()
{
   setName("INTEGRATION");
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
IntegrationAnalyzer::~IntegrationAnalyzer()
{
}

// ************************************************************************
// perform analysis
// ------------------------------------------------------------------------
double IntegrationAnalyzer::analyze(aData &adata)
{
   int    nInputs, nOutputs, nSamples, ss, ii, whichOutput;
   int    outputID, nLevels, *levelSeps;
   double result, last, error, *X, *Y, *iLower, *iUpper;

   nInputs   = adata.nInputs_;
   nOutputs  = adata.nOutputs_;
   nSamples  = adata.nSamples_;
   iLower    = adata.iLowerB_;
   iUpper    = adata.iUpperB_;
   X         = adata.sampleInputs_;
   Y         = adata.sampleOutputs_;
   outputID  = adata.outputID_;
   nLevels   = adata.currRefineLevel_;
   levelSeps = adata.refineSeparators_;
   if (adata.inputPDFs_ != NULL)
   {
      printf("IntegrationAnalyzer INFO: some inputs have non-uniform PDFs.\n");
      printf("          However, they will not be relevant in this analysis\n");
   }

   if (nInputs <= 0)
   {
      printf("IntegrationAnalyzer ERROR: invalid nInputs.\n");
      printf("    nInputs  = %d\n", nInputs);
      return PSUADE_UNDEFINED;
   } 
   if (nOutputs <= 0)
   {
      printf("IntegrationAnalyzer ERROR: invalid nOutputs.\n");
      printf("    nOutputs = %d\n", nOutputs);
      return PSUADE_UNDEFINED;
   } 
   if (nSamples <= 0)
   {
      printf("IntegrationAnalyzer ERROR: invalid nSamples.\n");
      printf("    nSamples = %d\n", nSamples);
      return PSUADE_UNDEFINED;
   } 
   whichOutput = outputID;
   if (whichOutput < 0 || whichOutput >= nOutputs)
   {
      printf("IntegrationAnalyzer ERROR: invalid outputID (%d).\n",
             whichOutput+1);
      return PSUADE_UNDEFINED;
   }
   for (ss = 0; ss < nSamples; ss++)
   {
      if (Y[ss] == PSUADE_UNDEFINED)
      {
         printf("IntegrationAnalyzer ERROR: some outputs are undefined.\n");
         printf("                           Prune them before analyze.\n");
         return PSUADE_UNDEFINED;
      }
   }

   result = 0.0;
   for (ss = 0; ss < nSamples; ss++) result += Y[ss];
   result /= (double) nSamples;
   for (ii = 0; ii < nInputs; ii++) result *= (iUpper[ii] - iLower[ii]);
   printAsterisks(0);
   printf("IntegrationAnalyzer: numerical integral = %14.4e\n", result);
   printAsterisks(0);

   if (nLevels <= 0) return result;

   last = 0.0;
   for (ss = 0; ss < levelSeps[nLevels-1]; ss++) last += Y[ss];
   last /= (double) levelSeps[nLevels-1];
   for (ii = 0; ii < nInputs; ii++) last *= (iUpper[ii] - iLower[ii]);

   if (result == 0.0) error = result;
   else               error = PABS((last-result)/result); 
   if (adata.printLevel_ > 0)
      printf("IntegrationAnalyzer: numerical error    = %14.4e\n", error);

   return error;
}

