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
// Functions for the class DeltaAnalyzer (Delta test)
//**/ Reference: "Variable selection for Financial Modeling"
//**/            by Q. Yu, E. Severin and A. Lendasse 
// ------------------------------------------------------------------------
// AUTHOR : CHARLES TONG/Michael Snow
// DATE   : 2009
// ************************************************************************
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <time.h>
#include "Psuade.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "DeltaAnalyzer.h"
#include "PrintingTS.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
DeltaAnalyzer::DeltaAnalyzer(): Analyzer(),nBins_(0),nInputs_(0)
{
  setName("DELTATEST");
  mode_ = 0;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
DeltaAnalyzer::~DeltaAnalyzer()
{
  VecMinDeltas_.clean();
  MatDeltaBins_.clean();
  VecOrders_.clean();
  VecRanks_.clean();
}

// ************************************************************************
// perform delta test
// ------------------------------------------------------------------------
double DeltaAnalyzer::analyze(aData &adata)
{
  int    ss, ss2, ii, jj, kk, ll, nBins=1000, *iPtr, converged=0, nIndex=3;
  int    nSelected=0, info, place, count, uniqueFlag=0, reverseCnt=0;
  int    iter=100, stagnate=100;
  double distance, delta, minDist, dtemp, minDelta, deltaSave, auxMin;
  double temperature=.0001, oldDelta=PSUADE_UNDEFINED;
  double bestDelta=PSUADE_UNDEFINED, alpha=0.98, r, ddata, accum;
  char   pString[500];
  FILE   *fp;

  //**/ ---------------------------------------------------------------
  //**/ extract data
  //**/ ---------------------------------------------------------------
  int printLevel = adata.printLevel_;
  double *XX     = adata.sampleInputs_;
  double *YY     = adata.sampleOutputs_;
  int nInputs    = adata.nInputs_;
  nInputs_       = nInputs;
  int nOutputs   = adata.nOutputs_;
  int nSamples   = adata.nSamples_;
  int outputID   = adata.outputID_;
  double *iLowerB = adata.iLowerB_;
  double *iUpperB = adata.iUpperB_;
  nBins_   = nBins;

  if (adata.inputPDFs_ != NULL)
  {
    count = 0;
    for (ii = 0; ii < nInputs; ii++) count += adata.inputPDFs_[ii];
    if (count > 0)
    {
      printOutTS(PL_WARN, 
           "DeltaTest INFO: some inputs have non-uniform PDFs, but\n");
      printOutTS(PL_WARN,
           "          they are not relevant in this analysis.\n");
    }
  }

  //**/ ---------------------------------------------------------------
  //**/ error checking
  //**/ ---------------------------------------------------------------
  if (nSamples <= 1)
  {
    printOutTS(PL_ERROR, 
         "DeltaTest INFO: test not meaningful for nSamples <= 1.\n");
    return PSUADE_UNDEFINED;
  }
  if (XX == NULL || YY == NULL)
  {
    printOutTS(PL_ERROR, "DeltaTest ERROR: no data.\n");
    return PSUADE_UNDEFINED;
  }
  info = 0;
  for (ii = 0; ii < nSamples; ii++)
    if (YY[nOutputs*ii+outputID] == PSUADE_UNDEFINED) info = 1;
  if (info == 1)
  {
    printOutTS(PL_ERROR, 
         "DeltaTest ERROR: Some outputs are undefined.\n");
    printOutTS(PL_ERROR, 
         "                 Prune the undefined's first.\n");
    return PSUADE_UNDEFINED;
  }

  //**/ ---------------------------------------------------------------
  //**/ clean up
  //**/ ---------------------------------------------------------------
  VecMinDeltas_.clean();
  MatDeltaBins_.clean();
  VecOrders_.clean();
  VecRanks_.clean();
 
  //**/ ---------------------------------------------------------------
  //**/ prepare to run
  //**/ ---------------------------------------------------------------
  printAsterisks(PL_INFO, 0);
  printOutTS(PL_INFO,"          Delta Test for variable selection\n");
  printDashes(PL_INFO, 0);
  printOutTS(PL_INFO,"This test has the characteristics that ");
  printOutTS(PL_INFO,"the more important a parameter\n");
  printOutTS(PL_INFO,"is relative to the others, the smaller the ");
  printOutTS(PL_INFO,"subset is at the end of the\n");
  printOutTS(PL_INFO,"the test (sharp zoom into the ");
  printOutTS(PL_INFO,"most important subset).\n");
  printOutTS(PL_INFO,"Thus, the purpose of this test is to ");
  printOutTS(PL_INFO,"identify a subset of important\n");
  printOutTS(PL_INFO,"parameters.\n");
  printOutTS(PL_INFO,"NOTE: If both nInputs and nSamples ");
  printOutTS(PL_INFO,"are large, this test may take a\n");
  printOutTS(PL_INFO,"      long time to run. So, be patient.)\n");
  printEquals(PL_INFO, 0);

  psIVector vecAuxBins, vecInpBins;
  vecAuxBins.setLength(nInputs);
  vecInpBins.setLength(nInputs);

  psIVector vecRanks;
  vecRanks.setLength(nInputs); 
  int *ranks = vecRanks.getIVector();
  VecRanks_.setLength(nInputs_);

  psVector vecOrders;
  vecOrders.setLength(nInputs);
  double *dOrder = vecOrders.getDVector(); 
  VecOrders_.setLength(nInputs_);

  psVector vecDistPairs, vecRangesInv2;
  vecDistPairs.setLength(nSamples*(nSamples-1)/2);
  vecRangesInv2.setLength(nInputs);

  psVector vecYT;
  vecYT.setLength(nSamples);
  for (ii = 0; ii < nSamples; ii++) vecYT[ii] = YY[ii*nOutputs+outputID];

  //**/ ---------------------------------------------------------------
  //**/ get user information
  //**/ ---------------------------------------------------------------
  if (psConfig_.AnaExpertModeIsOn())
  {
    printOutTS(PL_INFO,
         "DeltaTest Option: set the number of neighbors K.\n");
    printOutTS(PL_INFO, 
         "The larger K is, the larger the distinguishing power is.\n");
    sprintf(pString, "What is K (>= 1, <= 20, default=3)? ");
    nIndex = getInt(1, 20, pString);
    sprintf(pString,"How many inputs to select FOR SURE? (0 if not sure) ");
    nSelected = getInt(0, nInputs-1, pString);
    for (ii = 0; ii < nInputs; ii++) vecAuxBins[ii] = 0;
    for (ii = 0; ii < nSelected; ii++)
    {
      sprintf(pString,"Enter the %d-th input to be selected : ", ii+1);
      kk = getInt(1, nInputs, pString);
      vecAuxBins[kk-1] = 1;
    }
    sprintf(pString,"How many iterations for optimization? (> 100) ");
    iter = getInt(1, 100000, pString);
    printEquals(PL_INFO, 0);
  }

  //**/ ---------------------------------------------------------------
  //**/ initialize internal variables
  //**/ ---------------------------------------------------------------
  for (ii = 0; ii < nInputs; ii++) 
  {
    if ((iUpperB[ii] - iLowerB[ii]) > 0)
      vecRangesInv2[ii] = 1.0 / (iUpperB[ii] - iLowerB[ii]) /
                          (iUpperB[ii] - iLowerB[ii]);
    else
    {
      printOutTS(PL_ERROR, "DeltaTest ERROR: problem with input range.\n");
      exit(1);
    }
  }
  for (ii = 0; ii < nInputs; ii++) vecInpBins[ii] = 0;

  psIMatrix matLDeltaBins;
  MatDeltaBins_.setFormat(PS_MAT2D);
  MatDeltaBins_.setDim(nBins_, nInputs_);
  matLDeltaBins = MatDeltaBins_;
  int **deltaBins = matLDeltaBins.getIMatrix2D();

  VecMinDeltas_.setLength(nBins_);
  psVector vecLMinDeltas;
  vecLMinDeltas.setLength(nBins_);
  double *minDeltas = vecLMinDeltas.getDVector();
  for (ii = 0; ii < nBins; ii++) minDeltas[ii] = PSUADE_UNDEFINED;

  psIVector vecMinInds;
  vecMinInds.setLength(nIndex);

  //**/ ---------------------------------------------------------------
  //**/ generate initial configuration and distances
  //**/ ---------------------------------------------------------------
  if (nSelected == 0)
    for (ii = 0; ii < nInputs;ii++) vecInpBins[ii] = PSUADE_rand()%2;
  else
    for (ii = 0; ii < nInputs;ii++) vecInpBins[ii] = vecAuxBins[ii];

  for (ss = 1; ss < nSamples; ss++)
  {
    for (ss2 = 0; ss2 < ss; ss2++)
    {
      distance = 0.0;
      for (ii = 0; ii < nInputs; ii++)
      {
        if (vecInpBins[ii] == 1)
        {
          dtemp = XX[ss*nInputs+ii] - XX[ss2*nInputs+ii];
          distance += dtemp * dtemp * vecRangesInv2[ii];
        }
      }
      vecDistPairs[ss*(ss-1)/2+ss2] = distance;
    }
  }

  //**/ search for min distance 
  delta = 0.0;
  for (ss = 0; ss < nSamples; ss++)
  {
    ddata = 0.0;
    for (jj = 0; jj < nIndex; jj++)
    {
      minDist = PSUADE_UNDEFINED;
      vecMinInds[jj] = -1;
      for (ss2 = 0; ss2 < ss; ss2++)
      {
        kk = ss * (ss - 1) / 2 + ss2;
        if (vecDistPairs[kk] < minDist)
        {
          for (ll = 0; ll < jj; ll++)
            if (ss2 == vecMinInds[ll]) break;
          if (jj == 0 || ll == jj)
          {
            minDist = vecDistPairs[kk];
            vecMinInds[jj] = ss2;
          }
        }
      }
      for (ss2 = ss+1; ss2 < nSamples; ss2++)
      {
        kk = ss2 * (ss2 - 1) / 2 + ss;
        if (vecDistPairs[kk] < minDist)
        {
          for (ll = 0; ll < jj; ll++)
            if (ss2 == vecMinInds[ll]) break;
          if (jj == 0 || ll == jj)
          {
            minDist = vecDistPairs[kk];
            vecMinInds[jj] = ss2;
          }
        }
      }
      if (vecMinInds[jj] == -1)
      {
        printOutTS(PL_ERROR, "DeltaTest ERROR (1).\n");
        exit(1);
      }
      ddata += pow(vecYT[ss] - vecYT[vecMinInds[jj]], 2.0);
    }
    delta += ddata / (double) nIndex;
  }
  printOutTS(PL_INFO,"Current best solution for output %d:\n",outputID+1);
  printOutTS(PL_INFO,
     "To stop the search, create a psuade_stop file in local directory.\n");
  printDashes(PL_INFO, 0);
  delta /= (2.0 * nSamples);
  for (ii = 0; ii < nInputs; ii++) 
    printOutTS(PL_INFO, "%d ", vecInpBins[ii]);
  printOutTS(PL_INFO, " = %e\n", delta);
 
  //**/ ---------------------------------------------------------------
  //**/ going through all combinations
  //**/ ---------------------------------------------------------------
  count = 1;
  auxMin = - PSUADE_UNDEFINED;
  while (count <= 3*iter*nInputs)
  {
    fflush(stdout);
    count++;

    //**/ select or de-select the next input
    if (reverseCnt >= 4*nInputs)
    {	
      //**/printOutTS(PL_WARN, "local minima: kickstarted");
      //**/ re-initialize some of the inputs
      temperature *= nInputs * nInputs;
      for (ii = 0;ii <= nInputs/5;ii++) 
        vecInpBins[PSUADE_rand()%nInputs] ^=1;
      for (ss = 1; ss < nSamples; ss++)
      {
        for (ss2 = 0; ss2 < ss; ss2++)
        {
          distance = 0.0;
          for (ii = 0; ii < nInputs; ii++)
          {
            if (vecInpBins[ii] == 1)
            {
              dtemp = XX[ss*nInputs+ii] - XX[ss2*nInputs+ii];
              distance += dtemp * dtemp * vecRangesInv2[ii];
            }
          }
          vecDistPairs[ss*(ss-1)/2+ss2] = distance;
        }
      }
      reverseCnt = 0;
      place = PSUADE_rand()%(nInputs);
    }
    else 
    {
      if (reverseCnt >= 3*nInputs)
      {
        //printOutTS(PL_WARN, "suspected local minima, checking");
        place = reverseCnt - 3 * nInputs;
      }
      else 
      {
        place = PSUADE_rand()%(nInputs);
      }
    }
    temperature *= alpha;
    vecInpBins[place] ^= 1;

    //**/ update the distance
    delta = 0.0;
    for (ss = 1; ss < nSamples; ss++)
    {
      for (ss2 = 0; ss2 < ss; ss2++)
      {
        kk = ss * (ss - 1) / 2 + ss2;
        dtemp = XX[ss*nInputs+place] - XX[ss2*nInputs+place];
        if (vecInpBins[place] == 1)
             vecDistPairs[kk] += dtemp * dtemp * vecRangesInv2[place];
        else vecDistPairs[kk] -= dtemp * dtemp * vecRangesInv2[place];
      }
    }

    //**/ search for min distance 
    delta = 0.0;
    for (ss = 0; ss < nSamples; ss++)
    {
      ddata = 0.0;
      for (jj = 0; jj < nIndex; jj++)
      {
        minDist = PSUADE_UNDEFINED;
        vecMinInds[jj] = -1;
        for (ss2 = 0; ss2 < ss; ss2++)
        {
          kk = ss * (ss - 1) / 2 + ss2;
          if (vecDistPairs[kk] < minDist)
          {
            for (ll = 0; ll < jj; ll++)
              if (ss2 == vecMinInds[ll]) break;
            if (jj == 0 || ll == jj)
            {
              minDist = vecDistPairs[kk];
              vecMinInds[jj] = ss2;
            }
          }
        }
        for (ss2 = ss+1; ss2 < nSamples; ss2++)
        {
          kk = ss2 * (ss2 - 1) / 2 + ss;
          if (vecDistPairs[kk] < minDist)
          {
            for (ll = 0; ll < jj; ll++)
              if (ss2 == vecMinInds[ll]) break;
            if (jj == 0 || ll == jj)
            {
              minDist = vecDistPairs[kk];
              vecMinInds[jj] = ss2;
            }
          }
        }
        if (vecMinInds[jj] == -1)
        {
          printOutTS(PL_ERROR, "DeltaTest ERROR (1).\n");
          exit(1);
        }
        ddata += pow(vecYT[ss] - vecYT[vecMinInds[jj]], 2.0);
      }
      delta += ddata / (double) nIndex;
    }
    delta /= (2.0 * nSamples);

    //**/ store away the minimum configurations
    if (delta < minDeltas[0])
    {
      //**/ don't find unique, it will mess up ranking
      uniqueFlag = 1;
      //for (ii = 0; ii < nBins; ii++)
      //{
      //  if (minDeltas[ii] != PSUADE_UNDEFINED)
      //  {
      //    for (jj = 0; jj < nInputs; jj++)
      //      if (vecInpBins[jj] != deltaBins[ii][jj]) break;
      //    if (jj == nInputs) {uniqueFlag = 0; break;}
      //  }
      //}
      if (uniqueFlag == 1)
      {
        //**/ put the new minimum at the first position
        minDeltas[0] = delta;
        for (ii = 0; ii < nInputs; ii++) 
          deltaBins[0][ii] = vecInpBins[ii];

        //**/ now sort (compare and swap)
        for (ii = 1; ii < nBins; ii++)
        {
          if (minDeltas[ii] > minDeltas[ii-1])
          {
            dtemp = minDeltas[ii];
            minDeltas[ii] = minDeltas[ii-1];
            minDeltas[ii-1] = dtemp;
            iPtr = deltaBins[ii];
            deltaBins[ii] = deltaBins[ii-1];
            deltaBins[ii-1] = iPtr;
          }
        }
      }
    }

    for (ii = 0; ii < nInputs; ii++) 
      printOutTS(PL_INFO, "%d ", vecInpBins[ii]);
    printOutTS(PL_INFO," = %e (%d of %d)\n",delta,
               count-1,3*iter*nInputs);

    //**/ check convergence
    if (minDeltas[nBins-1] == auxMin) converged++;
    else
    {
      converged = 0;
      auxMin = minDeltas[nBins-1];
    }
    if (converged > stagnate*3*nInputs)
    {
      printOutTS(PL_INFO, "DeltaTest: stagnate for %d iterations, ", 
                 stagnate);
      printOutTS(PL_INFO, "considered converged.\n");
      break;
    }

    //**/ if current minimum is larger than the previous one,
    //**/ with a certain probability, reverse the selection
    if (delta >= oldDelta) 
    {
      r = PSUADE_rand()%100000;
      r /= 100000;
      //**/ printOutTS(PL_WARN, "\n dice roll r %e and comparison %e and
      //**/ temperature %e \n",r,
      //**/ exp(10000*(oldDelta-delta)/temperature),temperature);

      if (r>=exp(-.1*(delta-oldDelta)/(temperature)))
      {
        vecInpBins[place] ^=1;
        reverseCnt++;
        for (ss = 1; ss < nSamples; ss++)
        {
          for (ss2 = 0; ss2 < ss; ss2++)
          {
            kk = ss * (ss - 1) / 2 + ss2;
            dtemp = XX[ss*nInputs+place] - XX[ss2*nInputs+place];
            if (vecInpBins[place] == 1)
                 vecDistPairs[kk] += dtemp * dtemp * vecRangesInv2[place];
            else vecDistPairs[kk] -= dtemp * dtemp * vecRangesInv2[place];
          }
        }
      }
      else
      {
        oldDelta = delta;
        reverseCnt = 0;
      }
    }
    else 
    {
      oldDelta = delta;
      reverseCnt = 0;
    }
    if (oldDelta <= bestDelta) bestDelta = oldDelta;
    fp = fopen("psuade_stop","r");
    if (fp != NULL)
    {
      printOutTS(PL_INFO, "psuade_stop file found ==> terminate.\n");
      printOutTS(PL_INFO, "To restart, delete psuade_stop first.\n");
      fclose(fp);
      break;
    }
  }

  //**/ ---------------------------------------------------------------
  //**/ output the results
  //**/ ---------------------------------------------------------------
  printAsterisks(PL_INFO, 0);
  printOutTS(PL_INFO, 
       "Final Selections (based on %d neighbors) = \n", nIndex);
  printOutTS(PL_INFO, 
       "(The top is the best and then the next few best.)\n");
  printDashes(PL_INFO, 0);

  //save minDeltas and deltaBins
  for (ii=0; ii < nBins_; ii++)
  {
    VecMinDeltas_[ii] = minDeltas[ii];
    for (kk = 0; kk < nInputs_; kk++)
      MatDeltaBins_.setEntry(ii, kk, deltaBins[ii][kk]);
  }

  kk = count = 0;
  while (count < 10)
  {
    uniqueFlag = 1;
    for (ii = nBins-1; ii >= nBins-kk; ii--)
    {
      for (jj = 0; jj < nInputs; jj++)
        if (deltaBins[nBins-kk-1][jj] != deltaBins[ii][jj]) break;
      if (jj == nInputs) {uniqueFlag = 0; break;}
    } 
    if (uniqueFlag ==1 && minDeltas[nBins-kk-1] < 0.99 * PSUADE_UNDEFINED)
    {
      printOutTS(PL_INFO, "Rank %2d => ", count+1);
      for (ii = 0; ii < nInputs; ii++) 
        printOutTS(PL_INFO, "%d ", deltaBins[nBins-kk-1][ii]);
      printOutTS(PL_INFO, ": delta = %11.4e\n", minDeltas[nBins-kk-1]);
      count++;
    }
    kk++;
  }
  printDashes(PL_INFO, 0);
  count = 0;
  for (ii = 0; ii < nInputs; ii++)
  {
    ddata = 0;
    //**/for (kk = 0; kk < nBins; kk++) 
    //**/  count += deltaBins[nBins-kk-1][ii];
    //**/ Oct 2009: smallest data weighted more heavily
    accum = 0.0;
    for (kk = 0; kk < nBins; kk++)
    {
      if (minDeltas[nBins-kk-1] != PSUADE_UNDEFINED)
      {
        ddata += (minDeltas[nBins-1]*deltaBins[nBins-kk-1][ii]/
                  minDeltas[nBins-kk-1]);
        accum += (minDeltas[nBins-1]/minDeltas[nBins-kk-1]);
      }
    }
    //**/ printOutTS(PL_WARN, "%2d ", count);
    ranks[ii] = (int) (ddata / accum * 100);
  }

  //**/ output to a matlab file
  if (plotScilab()) fp = fopen("scilabdelta.sci", "w");
  else              fp = fopen("matlabdelta.m", "w");
  if (fp == NULL)
  {
    printOutTS(PL_INFO,"Delta test ERROR: cannot open graphics files.\n");
    printOutTS(PL_INFO,"                  ==> graphics not generated.\n");
  }
  else
  {
    fwritePlotCLF(fp);
    fprintf(fp, "A = [\n");
    for (ii = 0; ii < nInputs; ii++)
      fprintf(fp, "%e\n", 0.01 * ranks[ii]);
    fprintf(fp, "];\n");
    fprintf(fp, "bar(A, 0.8);\n");
    fwritePlotAxes(fp);
    fwritePlotTitle(fp, "Delta Test Rankings");
    fwritePlotXLabel(fp, "Input parameters");
    fwritePlotYLabel(fp, "Delta Metric (normalized)");
    fclose(fp);
    if (plotScilab()) 
      printOutTS(PL_INFO,
         "Delta test ranking is now in scilabdelta.sci.\n");
    else 
      printOutTS(PL_INFO,"Delta test ranking is now in matlabdelta.m.\n");
  }

  //**/ output the ranks in order
  for (ii = 0; ii < nInputs; ii++) dOrder[ii] = 1.0 * ii;
  sortIntList2a(nInputs, ranks, dOrder);
  printOutTS(PL_INFO,
    "Order of importance (based on frequencies of appearance in search)\n");
  for (ii = 0; ii < nInputs; ii++)
    printOutTS(PL_INFO, "(Delta) Rank %4d : input %4d (score = %d )\n", ii+1, 
               (int) dOrder[nInputs-ii-1]+1, ranks[nInputs-ii-1]);
  printAsterisks(PL_INFO, 0);

  //save dOrder and ranks
  for (ii = 0; ii < nInputs_; ii++)
  {
    VecOrders_[ii] = dOrder[ii];
    VecRanks_[ii]  = ranks[ii];
  }
  //**/ ---------------------------------------------------------------
  //**/ test delta values with most important parameters 
  //**/ ---------------------------------------------------------------
  printOutTS(PL_INFO, 
    "Final test adding the most important parameters incrementally:\n");
  printOutTS(PL_INFO,"You should see the rightmost values ");
  printOutTS(PL_INFO,"decreasing and then increasing.\n");
  printOutTS(PL_INFO,"The lowest point can be used as a separator ");
  printOutTS(PL_INFO,"for classifying important\n");
  printOutTS(PL_INFO,"and less important parameters.\n");
  printDashes(PL_INFO, 0);
  for (ii = 0; ii < nInputs; ii++) vecInpBins[ii] = 0;
  //**/ since it is not possible to compute the delta reduction due to
  //**/ the first input, the following code does the 2nd one first 
  for (ii = 1; ii >= 0; ii--)
  {
    vecInpBins[(int) dOrder[nInputs-ii-1]] = 1;
    delta = 0.0;
    for (ss = 0; ss < nSamples; ss++)
    {
      ddata = 0.0;
      for (jj = 0; jj < nIndex; jj++)
      {
        minDist = PSUADE_UNDEFINED;
        vecMinInds[jj] = -1;
        for (ss2 = 0; ss2 < ss; ss2++)
        {
          distance = 0.0;
          for (kk = 0; kk < nInputs; kk++)
          {
            if (vecInpBins[kk] == 1)
            {
              dtemp = XX[ss*nInputs+kk] - XX[ss2*nInputs+kk];
              distance += dtemp * dtemp * vecRangesInv2[kk];
            }
          }
          if (distance < minDist)
          {
            for (ll = 0; ll < jj; ll++)
              if (ss2 == vecMinInds[ll]) break;
            if (jj == 0 || ll == jj)
            {
              minDist = distance;
              vecMinInds[jj] = ss2;
            }
          }
        }
        for (ss2 = ss+1; ss2 < nSamples; ss2++)
        {
          distance = 0.0;
          for (kk = 0; kk < nInputs; kk++)
          {
            if (vecInpBins[kk] == 1)
            {
              dtemp = XX[ss*nInputs+kk] - XX[ss2*nInputs+kk];
              distance += dtemp * dtemp * vecRangesInv2[kk];
            }
          }
          if (distance < minDist)
          {
            for (ll = 0; ll < jj; ll++)
              if (ss2 == vecMinInds[ll]) break;
            if (jj == 0 || ll == jj)
            {
              minDist = distance;
              vecMinInds[jj] = ss2;
            }
          }
        }
        ddata += pow(vecYT[ss] - vecYT[vecMinInds[jj]], 2.0);
      }
      delta += ddata / (double) nIndex;
    }
    delta /= (2.0 * nSamples);
    minDeltas[ii] = delta;
  }
  deltaSave = minDeltas[1] - minDeltas[0];
  vecInpBins[(int) dOrder[nInputs-2]] = 0;
  vecInpBins[(int) dOrder[nInputs-1]] = 0;
  //**/ now start from the beginning
  minDelta = 1.0e35;
  for (ii = 0; ii < nInputs; ii++)
  {
    vecInpBins[(int) dOrder[nInputs-ii-1]] = 1;
    delta = 0.0;
    for (ss = 0; ss < nSamples; ss++)
    {
      ddata = 0.0;
      for (jj = 0; jj < nIndex; jj++)
      {
        minDist = PSUADE_UNDEFINED;
        vecMinInds[jj] = -1;
        for (ss2 = 0; ss2 < ss; ss2++)
        {
          distance = 0.0;
          for (kk = 0; kk < nInputs; kk++)
          {
            if (vecInpBins[kk] == 1)
            {
              dtemp = XX[ss*nInputs+kk] - XX[ss2*nInputs+kk];
              distance += dtemp * dtemp * vecRangesInv2[kk];
            }
          }
          if (distance < minDist)
          {
            for (ll = 0; ll < jj; ll++)
              if (ss2 == vecMinInds[ll]) break;
            if (jj == 0 || ll == jj)
            {
              minDist = distance;
              vecMinInds[jj] = ss2;
            }
          }
        }
        for (ss2 = ss+1; ss2 < nSamples; ss2++)
        {
          distance = 0.0;
          for (kk = 0; kk < nInputs; kk++)
          {
            if (vecInpBins[kk] == 1)
            {
              dtemp = XX[ss*nInputs+kk] - XX[ss2*nInputs+kk];
              distance += dtemp * dtemp * vecRangesInv2[kk];
            }
          }
          if (distance < minDist)
          {
            for (ll = 0; ll < jj; ll++)
              if (ss2 == vecMinInds[ll]) break;
            if (jj == 0 || ll == jj)
            {
              minDist = distance;
              vecMinInds[jj] = ss2;
            }
          }
        }
        ddata += pow(vecYT[ss] - vecYT[vecMinInds[jj]], 2.0);
      }
      delta += ddata / (double) nIndex;
    }
    delta /= (2.0 * nSamples);
    if (ii == 0)
    {
      deltaSave += delta;
      for (kk = 0; kk < nInputs; kk++) printOutTS(PL_INFO, "0 ");
      printOutTS(PL_INFO, " = %e\n", deltaSave);
    }
    for (kk = 0; kk < nInputs; kk++) 
      printOutTS(PL_INFO, "%d ", vecInpBins[kk]);
    printOutTS(PL_INFO, " = %e\n", delta);
    minDeltas[ii] = delta;
    if (delta < minDelta) minDelta = delta;
  }
  printAsterisks(PL_INFO, 0);
  return minDelta;
}

// ************************************************************************ 
// set parameters
// ------------------------------------------------------------------------
int DeltaAnalyzer::setParams(int argc, char **argv)
{
  char *request = (char *) argv[0];
  Analyzer::setParams(argc, argv);
  if (!strcmp(request, "gdelta")) mode_ = 1;
  return 0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
DeltaAnalyzer& DeltaAnalyzer::operator=(const DeltaAnalyzer &)
{
  printOutTS(PL_ERROR,"DeltaTest operator= ERROR: operation not allowed.\n");
  exit(1);
  return (*this);
}

// ************************************************************************
// functions for getting results
// ------------------------------------------------------------------------
int DeltaAnalyzer::get_mode()
{
   return mode_;
}
int DeltaAnalyzer::get_nBins()
{
  return nBins_;
}
int DeltaAnalyzer::get_nInputs()
{
  return nInputs_;
}
int DeltaAnalyzer::get_nConfig()
{
  return nConfig_;
}
double *DeltaAnalyzer::get_minDeltas()
{
  psVector vecDT;
  vecDT = VecMinDeltas_;
  return vecDT.getDVector();
}
int **DeltaAnalyzer::get_deltaBins()
{
  psIMatrix matIT;
  matIT = MatDeltaBins_;
  return matIT.takeIMatrix2D();
}
double *DeltaAnalyzer::get_dOrder()
{
  psVector vecDT;
  vecDT = VecOrders_;
  return vecDT.getDVector();
}
int *DeltaAnalyzer::get_ranks()
{
  psIVector vecIT;
  vecIT = VecRanks_;
  return vecIT.getIVector();
}

