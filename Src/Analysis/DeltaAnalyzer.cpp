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
// ------------------------------------------------------------------------
// AUTHOR : CHARLES TONG/Michael Snow
// DATE   : 2009
// ************************************************************************
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "Main/Psuade.h"
#include "Util/sysdef.h"
#include "Util/PsuadeUtil.h"
#include "DeltaAnalyzer.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
DeltaAnalyzer::DeltaAnalyzer() : Analyzer()
{
   setName("DELTATEST");
   mode_ = 0;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
DeltaAnalyzer::~DeltaAnalyzer()
{
}

// ************************************************************************
// perform delta test
// ------------------------------------------------------------------------
double DeltaAnalyzer::analyze(aData &adata)
{
   int    printLevel, nSamples, nInputs, nOutputs, outputID, ss, ss2;
   int    *inputBins, *auxBins, nSelected=0, *minIndices, nIndex=3, info;
   int    ii, jj, kk, ll, **deltaBins, nBins=1000, *iPtr, converged=0;
   int    place, count=1, uniqueFlag=0, reverseCnt=0, *ranks;
   int    nConfig=20, iter=100, stagnate=100;
   double *X, *Y, distance, delta, minDist, *minDeltas, *iLowerB, *iUpperB;
   double dtemp, temperature=.0001, oldDelta=PSUADE_UNDEFINED, minDelta;
   double bestDelta=PSUADE_UNDEFINED, *dOrder, *rangesInv2, *distPairs;
   double *YY, alpha=0.98, r, ddata, accum, auxMin, deltaSave;
   char   pString[500];
   FILE   *fp;

   printLevel = adata.printLevel_;
   nSamples   = adata.nSamples_;
   X          = adata.sampleInputs_;
   YY         = adata.sampleOutputs_;
   nInputs    = adata.nInputs_;
   nOutputs   = adata.nOutputs_;
   nSamples   = adata.nSamples_;
   outputID   = adata.outputID_;
   iLowerB    = adata.iLowerB_;
   iUpperB    = adata.iUpperB_;
   if (adata.inputPDFs_ != NULL)
   {
      printf("DeltaAnalyzer INFO: some inputs have non-uniform PDFs, but\n");
      printf("                    they will not be relevant in this analysis.\n");
   }

   if (nSamples <= 1)
   {
      printf("DeltaAnalyzer INFO: not meaningful to do this");
      printf("                    test when nSamples <= 1.\n");
      return PSUADE_UNDEFINED;
   }
   if (X == NULL || YY == NULL)
   {
      printf("DeltaAnalyzer ERROR: no data.\n");
      return PSUADE_UNDEFINED;
   }
   info = 0;
   for (ii = 0; ii < nSamples; ii++)
      if (YY[nOutputs*ii+outputID] == PSUADE_UNDEFINED) info = 1;
   if (info == 1)
   {
      printf("DeltaAnalyzer ERROR: Some outputs are undefined.\n");
      printf("                     Prune the undefined's first.\n");
      return PSUADE_UNDEFINED;
   }

   printAsterisks(0);
   printf("DeltaAnalyzer for variable selection\n");
   printf("This test has the characteristics that the more important\n");
   printf("a parameter is relative to the others, the smaller the \n");
   printf("subset is at the end of the test (sharp zoom into the most\n");
   printf("important subset).\n");
   printf("Note: If both nInputs and nSamples are large, this test\n");
   printf("      may take a long time to run. So, be patient.)\n");
   printEquals(0);
   auxBins = new int[nInputs];
   inputBins = new int[nInputs];
   ranks = new int[nInputs];
   dOrder = new double[nInputs];
   rangesInv2 = new double[nInputs];
   distPairs = new double[nSamples*(nSamples-1)/2];
   Y = new double[nSamples];
   for (ii = 0; ii < nSamples; ii++) Y[ii] = YY[ii*nOutputs+outputID];

   if (psAnaExpertMode_ == 1)
   {
      printf("DeltaAnalyzer Option: set the number of neighbors K.\n");
      printf("The larger K is, the larger the distinguishing power is.\n");
      sprintf(pString, "What is K (>= 1, <= 20, default=3)? ");
      nIndex = getInt(1, 20, pString);
      sprintf(pString,"How many inputs to select FOR SURE? (0 if not sure) ");
      nSelected = getInt(0, nInputs-1, pString);
      for (ii = 0; ii < nInputs; ii++) auxBins[ii] = 0;
      for (ii = 0; ii < nSelected; ii++)
      {
         sprintf(pString,"Enter the %d-th input to be selected : ", ii+1);
         kk = getInt(1, nInputs, pString);
         auxBins[kk-1] = 1;
      }
      sprintf(pString,"How many iterations for optimization? (> 100) ");
      iter = getInt(1, 100000, pString);
      //sprintf(pString,"How many configurations for ranking? (1 if not sure) ");
      //nConfig = getInt(1, nBins, pString);
      printEquals(0);
   }

   for (ii = 0; ii < nInputs; ii++) 
   {
      if ((iUpperB[ii] - iLowerB[ii]) > 0)
         rangesInv2[ii] = 1.0 / (iUpperB[ii] - iLowerB[ii]) /
                                (iUpperB[ii] - iLowerB[ii]);
      else
      {
         printf("DeltaAnalyzer ERROR: problem with input range.\n");
         exit(1);
      }
   }
   for (ii = 0; ii < nInputs; ii++) inputBins[ii] = 0;
   deltaBins = new int*[nBins];
   for (ii = 0; ii < nBins; ii++)
   {
      deltaBins[ii] = new int[nInputs];
      for (jj = 0; jj < nInputs; jj++) deltaBins[ii][jj] = 0;
   }
   minDeltas = new double[nBins];
   for (ii = 0; ii < nBins; ii++) minDeltas[ii] = PSUADE_UNDEFINED;
   minIndices = new int[nIndex];

   srand(time(NULL));  
   if (nSelected == 0)
      for (ii = 0; ii < nInputs;ii++) inputBins[ii]=random()%2;
   else
      for (ii = 0; ii < nInputs;ii++) inputBins[ii]=auxBins[ii];

   for (ss = 1; ss < nSamples; ss++)
   {
      for (ss2 = 0; ss2 < ss; ss2++)
      {
         distance = 0.0;
         for (ii = 0; ii < nInputs; ii++)
         {
            if (inputBins[ii] == 1)
            {
               dtemp = X[ss*nInputs+ii] - X[ss2*nInputs+ii];
               distance += dtemp * dtemp * rangesInv2[ii];
            }
         }
         distPairs[ss*(ss-1)/2+ss2] = distance;
      }
   }

   delta = 0.0;
   for (ss = 0; ss < nSamples; ss++)
   {
      ddata = 0.0;
      for (jj = 0; jj < nIndex; jj++)
      {
         minDist = PSUADE_UNDEFINED;
         minIndices[jj] = -1;
         for (ss2 = 0; ss2 < ss; ss2++)
         {
            kk = ss * (ss - 1) / 2 + ss2;
            if (distPairs[kk] < minDist)
            {
               for (ll = 0; ll < jj; ll++)
                  if (ss2 == minIndices[ll]) break;
               if (jj == 0 || ll == jj)
               {
                  minDist = distPairs[kk];
                  minIndices[jj] = ss2;
               }
            }
         }
         for (ss2 = ss+1; ss2 < nSamples; ss2++)
         {
            kk = ss2 * (ss2 - 1) / 2 + ss;
            if (distPairs[kk] < minDist)
            {
               for (ll = 0; ll < jj; ll++)
                  if (ss2 == minIndices[ll]) break;
               if (jj == 0 || ll == jj)
               {
                  minDist = distPairs[kk];
                  minIndices[jj] = ss2;
               }
            }
         }
         if (minIndices[jj] == -1)
         {
            printf("DeltaAnalyzer ERROR (1).\n");
            exit(1);
         }
         ddata += pow(Y[ss] - Y[minIndices[jj]], 2.0);
      }
      delta += ddata / (double) nIndex;
   }
   printf("Current best solution for output %d:\n",outputID+1);
   printDashes(0);
   delta /= (2.0 * nSamples);
   for (ii = 0; ii < nInputs; ii++) printf("%d ", inputBins[ii]);
   printf(" = %e\n", delta);
 
   count = 1;
   auxMin = - PSUADE_UNDEFINED;
   while (count <= 3*iter*nInputs)
   {
      //if ((count % (3*nInputs) == 0)) printf("%d%% ", count/(3*nInputs));
      fflush(stdout);
      count++;

      if (reverseCnt >= 4*nInputs)
      {	
         temperature*=nInputs*nInputs;
         for (ii = 0;ii <= nInputs/5;ii++) inputBins[random()%nInputs] ^=1;
         for (ss = 1; ss < nSamples; ss++)
         {
            for (ss2 = 0; ss2 < ss; ss2++)
            {
               distance = 0.0;
               for (ii = 0; ii < nInputs; ii++)
               {
                  if (inputBins[ii] == 1)
                  {
                     dtemp = X[ss*nInputs+ii] - X[ss2*nInputs+ii];
                     distance += dtemp * dtemp * rangesInv2[ii];
                  }
               }
               distPairs[ss*(ss-1)/2+ss2] = distance;
            }
         }
         reverseCnt = 0;
         place = random()%(nInputs);
      }
      else 
      {
         if (reverseCnt >= 3*nInputs)
         {
            //printf("suspected local minima, checking");
            place = reverseCnt - 3 * nInputs;
         }
         else 
         {
            place = random()%(nInputs);
         }
      }
      temperature *= alpha;
      inputBins[place] ^= 1;

      delta = 0.0;
      for (ss = 1; ss < nSamples; ss++)
      {
         for (ss2 = 0; ss2 < ss; ss2++)
         {
            kk = ss * (ss - 1) / 2 + ss2;
            dtemp = X[ss*nInputs+place] - X[ss2*nInputs+place];
            if (inputBins[place] == 1)
                 distPairs[kk] += dtemp * dtemp * rangesInv2[place];
            else distPairs[kk] -= dtemp * dtemp * rangesInv2[place];
         }
      }

      delta = 0.0;
      for (ss = 0; ss < nSamples; ss++)
      {
         ddata = 0.0;
         for (jj = 0; jj < nIndex; jj++)
         {
            minDist = PSUADE_UNDEFINED;
            minIndices[jj] = -1;
            for (ss2 = 0; ss2 < ss; ss2++)
            {
               kk = ss * (ss - 1) / 2 + ss2;
               if (distPairs[kk] < minDist)
               {
                  for (ll = 0; ll < jj; ll++)
                     if (ss2 == minIndices[ll]) break;
                  if (jj == 0 || ll == jj)
                  {
                     minDist = distPairs[kk];
                     minIndices[jj] = ss2;
                  }
               }
            }
            for (ss2 = ss+1; ss2 < nSamples; ss2++)
            {
               kk = ss2 * (ss2 - 1) / 2 + ss;
               if (distPairs[kk] < minDist)
               {
                  for (ll = 0; ll < jj; ll++)
                     if (ss2 == minIndices[ll]) break;
                  if (jj == 0 || ll == jj)
                  {
                     minDist = distPairs[kk];
                     minIndices[jj] = ss2;
                  }
               }
            }
            if (minIndices[jj] == -1)
            {
               printf("DeltaAnalyzer ERROR (1).\n");
               exit(1);
            }
            ddata += pow(Y[ss] - Y[minIndices[jj]], 2.0);
         }
         delta += ddata / (double) nIndex;
      }
      delta /= (2.0 * nSamples);
      if ((count % (3*nInputs) == 0))
      {
         for (ii = 0; ii < nInputs; ii++) printf("%d ", deltaBins[nBins-1][ii]);
         printf(" = %e (%d of %d)\n", bestDelta, count/(3*nInputs), iter);
      }

      if (delta < minDeltas[0])
      {
         uniqueFlag = 1;
         for (ii = 0; ii < nBins; ii++)
         {
            if (minDeltas[ii] != PSUADE_UNDEFINED)
            {
               for (jj = 0; jj < nInputs; jj++)
                 if (inputBins[jj] != deltaBins[ii][jj]) break;
               if (jj == nInputs) {uniqueFlag = 0; break;}
            }
         }
         if (uniqueFlag == 1)
         {
            minDeltas[0] = delta;
            for (ii = 0; ii < nInputs; ii++) deltaBins[0][ii] = inputBins[ii];
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

      if (minDeltas[nBins-1] == auxMin) converged++;
      else
      {
         converged = 0;
         auxMin = minDeltas[nBins-1];
      }
      if (converged > stagnate*3*nInputs)
      {
         printf("DeltaAnalyzer: stagnate for %d iterations, ", stagnate); 
         printf("considered converged.\n"); 
         break;
      }

      if (delta >= oldDelta) 
      {
         r = random()%100000;
         r /= 100000;

         if (r>=exp(-.1*(delta-oldDelta)/(temperature)))
         {
            inputBins[place] ^=1;
            reverseCnt++;
            for (ss = 1; ss < nSamples; ss++)
            {
               for (ss2 = 0; ss2 < ss; ss2++)
               {
                  kk = ss * (ss - 1) / 2 + ss2;
                  dtemp = X[ss*nInputs+place] - X[ss2*nInputs+place];
                  if (inputBins[place] == 1)
                       distPairs[kk] += dtemp * dtemp * rangesInv2[place];
                  else distPairs[kk] -= dtemp * dtemp * rangesInv2[place];
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
   }

   printAsterisks(0);
   printf("Final Selections (based on %d neighbors) = \n", nIndex);
   for (kk = 0; kk < 10; kk++)
   {
      printf("Rank %2d => ", kk+1);
      for (ii = 0; ii < nInputs; ii++) printf("%d ", deltaBins[nBins-kk-1][ii]);
      printf(": delta = %11.4e\n", minDeltas[nBins-kk-1]);
   }
   printDashes(0);
   count = 0;
   for (ii = 0; ii < nInputs; ii++)
   {
      ddata = 0;
      accum = 0.0;
      for (kk = 0; kk < nConfig; kk++)
      {
         if (minDeltas[nBins-kk-1] != PSUADE_UNDEFINED)
         {
            ddata += (minDeltas[nBins-1]*deltaBins[nBins-kk-1][ii]/
                      minDeltas[nBins-kk-1]);
            accum += (minDeltas[nBins-1]/minDeltas[nBins-kk-1]);
         }
      }
      ranks[ii] = (int) (ddata / accum * 100);
   }

   fp = fopen("matlabdelta.m", "w");
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
   printf("Delta test ranking is now in matlabdelta.m.\n");

   for (ii = 0; ii < nInputs; ii++) dOrder[ii] = 1.0 * ii;
   sortIntList2a(nInputs, ranks, dOrder);
   printf("Order of importance (based on %d best configurations): \n",
          nConfig);
   for (ii = 0; ii < nInputs; ii++)
      printf("(D)Rank %4d : input %4d (score = %d )\n", ii+1, (int) 
             dOrder[nInputs-ii-1]+1, ranks[nInputs-ii-1]);
   printAsterisks(0);

   printf("Final test using the most important parameters incrementally:\n");
   printDashes(0);
   for (ii = 0; ii < nInputs; ii++) inputBins[ii] = 0;
   for (ii = 1; ii >= 0; ii--)
   {
      inputBins[(int) dOrder[nInputs-ii-1]] = 1;
      delta = 0.0;
      for (ss = 0; ss < nSamples; ss++)
      {
         ddata = 0.0;
         for (jj = 0; jj < nIndex; jj++)
         {
            minDist = PSUADE_UNDEFINED;
            minIndices[jj] = -1;
            for (ss2 = 0; ss2 < ss; ss2++)
            {
               distance = 0.0;
               for (kk = 0; kk < nInputs; kk++)
               {
                  if (inputBins[kk] == 1)
                  {
                     dtemp = X[ss*nInputs+kk] - X[ss2*nInputs+kk];
                     distance += dtemp * dtemp * rangesInv2[kk];
                  }
               }
               if (distance < minDist)
               {
                  for (ll = 0; ll < jj; ll++)
                     if (ss2 == minIndices[ll]) break;
                  if (jj == 0 || ll == jj)
                  {
                     minDist = distance;
                     minIndices[jj] = ss2;
                  }
               }
            }
            for (ss2 = ss+1; ss2 < nSamples; ss2++)
            {
               distance = 0.0;
               for (kk = 0; kk < nInputs; kk++)
               {
                  if (inputBins[kk] == 1)
                  {
                     dtemp = X[ss*nInputs+kk] - X[ss2*nInputs+kk];
                     distance += dtemp * dtemp * rangesInv2[kk];
                  }
               }
               if (distance < minDist)
               {
                  for (ll = 0; ll < jj; ll++)
                     if (ss2 == minIndices[ll]) break;
                  if (jj == 0 || ll == jj)
                  {
                     minDist = distance;
                     minIndices[jj] = ss2;
                  }
               }
            }
            ddata += pow(Y[ss] - Y[minIndices[jj]], 2.0);
         }
         delta += ddata / (double) nIndex;
      }
      delta /= (2.0 * nSamples);
      minDeltas[ii] = delta;
   }
   deltaSave = minDeltas[1] - minDeltas[0];
   inputBins[(int) dOrder[nInputs-2]] = 0;
   inputBins[(int) dOrder[nInputs-1]] = 0;
   minDelta = 1.0e35;
   for (ii = 0; ii < nInputs; ii++)
   {
      inputBins[(int) dOrder[nInputs-ii-1]] = 1;
      delta = 0.0;
      for (ss = 0; ss < nSamples; ss++)
      {
         ddata = 0.0;
         for (jj = 0; jj < nIndex; jj++)
         {
            minDist = PSUADE_UNDEFINED;
            minIndices[jj] = -1;
            for (ss2 = 0; ss2 < ss; ss2++)
            {
               distance = 0.0;
               for (kk = 0; kk < nInputs; kk++)
               {
                  if (inputBins[kk] == 1)
                  {
                     dtemp = X[ss*nInputs+kk] - X[ss2*nInputs+kk];
                     distance += dtemp * dtemp * rangesInv2[kk];
                  }
               }
               if (distance < minDist)
               {
                  for (ll = 0; ll < jj; ll++)
                     if (ss2 == minIndices[ll]) break;
                  if (jj == 0 || ll == jj)
                  {
                     minDist = distance;
                     minIndices[jj] = ss2;
                  }
               }
            }
            for (ss2 = ss+1; ss2 < nSamples; ss2++)
            {
               distance = 0.0;
               for (kk = 0; kk < nInputs; kk++)
               {
                  if (inputBins[kk] == 1)
                  {
                     dtemp = X[ss*nInputs+kk] - X[ss2*nInputs+kk];
                     distance += dtemp * dtemp * rangesInv2[kk];
                  }
               }
               if (distance < minDist)
               {
                  for (ll = 0; ll < jj; ll++)
                     if (ss2 == minIndices[ll]) break;
                  if (jj == 0 || ll == jj)
                  {
                     minDist = distance;
                     minIndices[jj] = ss2;
                  }
               }
            }
            ddata += pow(Y[ss] - Y[minIndices[jj]], 2.0);
         }
         delta += ddata / (double) nIndex;
      }
      delta /= (2.0 * nSamples);
      if (ii == 0)
      {
         deltaSave += delta;
         for (kk = 0; kk < nInputs; kk++) printf("0 ");
         printf(" = %e\n", deltaSave);
      }
      for (kk = 0; kk < nInputs; kk++) printf("%d ", inputBins[kk]);
      printf(" = %e\n", delta);
      minDeltas[ii] = delta;
      if (delta < minDelta) minDelta = delta;
   }
   printAsterisks(0);

   delete [] auxBins;
   delete [] inputBins;
   for (ii = 0; ii < nBins; ii++) delete [] deltaBins[ii];
   delete [] deltaBins;
   delete [] minDeltas;
   delete [] ranks;
   delete [] dOrder;
   delete [] rangesInv2;
   delete [] distPairs;
   delete [] minIndices;
   delete [] Y;
   return minDelta;
}

// ************************************************************************ 
// set parameters
// ------------------------------------------------------------------------
int DeltaAnalyzer::setParams(int argc, char **argv)
{
   char *request = (char *) argv[0];
   if (!strcmp(request, "gdelta")) mode_ = 1;
   return 0;
}

