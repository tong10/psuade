// ************************************************************************
// Copyright (c) 2007   Lawrence Livermore National Security, LLC. 
// Produced at the Lawrence Livermore National Laboratory.
// Written by the PSUADE team.
// All rights reserved.
//
// Please see the COPYRIGHT_and_LICENSE file for the copyright notice,
// disclaimer, contact information and the GNU Lesser General Public License.
//
// PSUADE is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License (as published by the Free Software
// Foundation) version 2.1 dated February 1999.
//
// PSUADE is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the terms and conditions of the GNU General
// Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
// ************************************************************************
// Functions for the TSIAnalyzer class 
// AUTHOR : CHARLES TONG
// DATE   : 2020
// ************************************************************************
#include <stdio.h>
#include <assert.h>
#include <sstream>
using namespace std;

#include "Main/Psuade.h"
#include "TSIAnalyzer.h"
#include "Util/PsuadeUtil.h"

#ifdef HAVE_METIS
extern "C" 
{
void METIS_PartGraphRecursive(int *, int *, int *, int *, int *,
                              int *, int *, int *, int *, int *, int *);
}
#endif

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
TSIAnalyzer::TSIAnalyzer() : Analyzer()
{
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
TSIAnalyzer::~TSIAnalyzer()
{
}

// ************************************************************************
// initialize the sampler data
// ------------------------------------------------------------------------
double TSIAnalyzer::analyze(aData &adata)
{
   int    ss, ii, jj, nInputs, nOutputs, nSamples, outputID, inputID, nnz;
   int    *incrs, itmp, jtmp, *sample2Aggr, *cellsOccupied, n1d, nAggrs;
   int    options[10], *aggrCnts;
#ifdef HAVE_METIS
   int    wgtflag=0, numflag=0, edgeCut=0;
#endif
   int    graphN, *graphI, *graphJ, printLevel, index, totalCnt;
   double *lbounds, *ubounds, *ranges, *X, *Y, *YY, *aggrMean, dvar, dmean;
   double variance, dtmp;
   char   pString[500];

   printLevel = adata.printLevel_;
   nInputs  = adata.nInputs_;
   nOutputs = adata.nOutputs_;
   nSamples = adata.nSamples_;
   X        = adata.sampleInputs_;
   YY       = adata.sampleOutputs_;
   outputID = adata.outputID_;
   lbounds  = adata.iLowerB_;
   ubounds  = adata.iUpperB_;

   if (nInputs <= 0 || nOutputs <= 0)
   {
      printf("TSIAnalyzer ERROR: invalid nInputs or nOutputs.\n");
      printf("   nInputs  = %d\n", nInputs);
      printf("   nOutputs = %d\n", nOutputs);
      return -1;
   }  
   if (nSamples <= 1)
   {
      printf("TSIAnalyzer ERROR: nSamples should be > 1.\n");
      printf("   nSamples = %d\n", nSamples);
      return -1;
   }  
   Y = new double[nSamples];
   for (ss = 0; ss < nSamples; ss++) Y[ss] = YY[ss*nOutputs+outputID];

   ranges  = new double[nInputs];
   for (ii = 0; ii < nInputs; ii++)
   {
      ranges[ii] = ubounds[ii] - lbounds[ii];
      if (ranges[ii] <= 0.0)
      {
         printf("TSIAnalyzer::analyze ERROR - lbound/ubound mismatch.\n");
         exit(1);
      }
   }

   dmean = 0.0;
   for (ss = 0; ss < nSamples; ss++) dmean += Y[ss];
   dmean /= (double) nSamples;
   variance = 0.0;
   for (ss = 0; ss < nSamples; ss++)
   {
      variance += ((Y[ss] - dmean) * (Y[ss] - dmean));
   }
   variance /= (double) (nSamples - 1);
   printf("TSIAnalyzer: output mean     = %e\n", dmean);
   printf("TSIAnalyzer: output variance = %e\n", variance);

   if (nInputs > 21)
   {
      printf("TSIAnalyzer ERROR : nInputs > 21 currently not supported.\n");
      exit(1);
   }
   if (nInputs == 1 ) n1d = nSamples*10;
   if (nInputs == 2 ) n1d = 1024;
   if (nInputs == 3 ) n1d = 100;
   if (nInputs == 4 ) n1d = 36;
   if (nInputs == 5 ) n1d = 16;
   if (nInputs == 6 ) n1d = 11;
   if (nInputs == 7 ) n1d = 8;
   if (nInputs == 8 ) n1d = 6;
   if (nInputs == 9 ) n1d = 5;
   if (nInputs == 10) n1d = 4;
   if (nInputs == 11) n1d = 3;
   if (nInputs == 12) n1d = 3;
   if (nInputs == 13) n1d = 3;
   if (nInputs >= 14) n1d = 2;

   printAsterisks(0);
   printf("           Crude Total Sensitivity Indices\n");
   printEquals(0);
   printf("TSIAnalyzer: number of subdomains = %d\n", nSamples/10);
   printf("This may need to be adjusted for higher accuracy.\n");
   printf("Recommendation: Try different numbers of subdomains.\n");
   printf("                Generally, it is better to have a\n");
   printf("                large number of subdomains.\n");
   printf("Turn on analysis expert mode to change it.\n");
   if (psAnaExpertMode_ != 0)
   {
      strcpy(pString, "Enter the number of subdomains (> 5): ");
      nAggrs = getInt(5, nSamples, pString);
   }
   else nAggrs = nSamples / 10;
   if (printLevel > 1)
      printf("TSIAnalyzer: number of subdomains = %d\n", nAggrs);

   incrs  = new int[nInputs];
   graphN = 1;
   incrs[0] = graphN;
   for (jj = 1; jj < nInputs; jj++)
   {
      graphN *= n1d;
      incrs[jj] = graphN;
   }
   graphI = new int[graphN+1];
   graphJ = new int[graphN*(nInputs-1)*2+1];
   nnz = 0;
   graphI[0] = nnz;
   for (ii = 0; ii < graphN; ii++)
   {
      itmp = ii;
      for (jj = 0; jj < nInputs-1; jj++)
      {
         jtmp = itmp % n1d;
         itmp = itmp / n1d;
         if (jtmp > 0     ) graphJ[nnz++] = ii - incrs[jj];
         if (jtmp < n1d-1) graphJ[nnz++] = ii + incrs[jj];
      }
      graphI[ii+1] = nnz;
   }
   cellsOccupied = new int[graphN];
   sample2Aggr = new int[nSamples];
   aggrMean = new double[nAggrs];
   aggrCnts = new int[nAggrs];

   options[0] = 0;
#ifdef HAVE_METIS
   METIS_PartGraphRecursive(&graphN, graphI, graphJ, NULL, NULL,
             &wgtflag,&numflag,&nAggrs,options,&edgeCut,cellsOccupied);
#else
   printf("TSIAnalyzer ERROR : METIS not installed.\n");
   nInputs = 0;
#endif

   for (inputID = 0; inputID < nInputs; inputID++)
   {
      for (ss = 0; ss < nSamples; ss++)
      {
         itmp = 0;
         for (ii = 0; ii < nInputs; ii++)
         {
            if (ii != inputID)
            {
               itmp = itmp * n1d;
               dtmp = X[ss*nInputs+ii];
               dtmp = (dtmp - lbounds[ii]) / ranges[ii];
               jtmp = (int) (dtmp * n1d);
               itmp += jtmp;
            }
         }
         sample2Aggr[ss] = cellsOccupied[itmp];
      }

      for (ii = 0; ii < nAggrs; ii++)
      {
         aggrMean[ii] = 0.0;
         aggrCnts[ii] = 0;
      }
      for (ss = 0; ss < nSamples; ss++)
      {
          index = sample2Aggr[ss];
          aggrMean[index] += Y[ss];
          aggrCnts[index]++;
      }
      totalCnt = 0;
      for (ii = 0; ii < nAggrs; ii++) totalCnt += aggrCnts[ii];
      for (ii = 0; ii < nAggrs; ii++)
         if (aggrCnts[ii] > 0) aggrMean[ii] /= (double) aggrCnts[ii];
      dmean = 0.0;
      for (ii = 0; ii < nAggrs; ii++)
          dmean += aggrMean[ii] * aggrCnts[ii];
      dmean /= (double) nSamples;
      dvar = 0.0;
      for (ii = 0; ii < nAggrs; ii++)
         dvar += pow(aggrMean[ii] - dmean, 2.0);
      dvar /= (double) (nAggrs - 1.0);

      if (dvar < variance)
         printf("Input %4d : Approximate total sensitivity index = %e\n",
                inputID+1, 1.0-dvar/variance);
      else
         printf("Input %4d : Approximate total sensitivity index = 0.\n",
                inputID+1);
   }

   delete [] ranges;
   delete [] Y;
   delete [] aggrMean;
   delete [] aggrCnts;
   delete [] incrs;
   delete [] sample2Aggr;
   if (cellsOccupied != NULL) delete [] cellsOccupied;
   if (graphI != NULL) delete [] graphI;
   if (graphJ != NULL) delete [] graphJ;
   return 0;
}

