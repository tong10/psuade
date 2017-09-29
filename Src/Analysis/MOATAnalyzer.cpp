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
// Functions for the class MOATAnalyzer  
// AUTHOR : CHARLES TONG
// DATE   : 2004
// ************************************************************************
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Samplings/MOATConstraints.h"
#include "Util/PsuadeUtil.h"
#include "Util/sysdef.h"
#include "MOATAnalyzer.h"
#include "BinomialAnalyzer.h"
#include "BootstrapAnalyzer.h"
#include "Main/Psuade.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
MOATAnalyzer::MOATAnalyzer() : Analyzer()
{
   setName("MOAT");
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
MOATAnalyzer::~MOATAnalyzer()
{
}

// ************************************************************************
// perform analysis
// ------------------------------------------------------------------------
double MOATAnalyzer::analyze(aData &adata)
{
   int    nInputs, nOutputs, nSamples, outputID, *S, ii, jj, diffIndex;
   int    iD, *indexTrack, index, n1, n2, flag, actualPaths, sigFlag;
   int    diffCnt, *counts, itemp, iaCnt, iD2, nPaths, ip, printLevel;
   double *X, *Y, *YY, *YG, *XG, *Xbase;
   double *means, *modifiedMeans, *stds, *modifiedStds, xtemp1, xtemp2;
   double ytemp1, ytemp2, scale, *xLower, *xUpper, dtemp, dtemp2, thresh;
   double *XSort, *YSort, dsum, dstdev, *YB, binMax=-1.0e35, binMin=1.0e35;
   char   winput[500];
   char   pString[500];
   PsuadeData *ioPtr;
   BinomialAnalyzer  binAnalyzer;
   BootstrapAnalyzer bsAnalyzer;
   aData             bData;
   pData             qData;
   MOATConstraints   *constrPtr;

   printLevel = adata.printLevel_;
   nInputs    = adata.nInputs_;
   nOutputs   = adata.nOutputs_;
   nSamples   = adata.nSamples_;
   xLower     = adata.iLowerB_;
   xUpper     = adata.iUpperB_;
   X          = adata.sampleInputs_;
   Y          = adata.sampleOutputs_;
   S          = adata.sampleStates_;
   outputID   = adata.outputID_;
   ioPtr      = adata.ioPtr_;
   if (adata.inputPDFs_ != NULL)
   {
      printf("MOATAnalyzer INFO: some inputs have non-uniform PDFs.\n");
      printf("          However, they will not be relevant in this analysis\n");
      printf("          (since the sample should have been generated with\n");
      printf("           the desired distributions.)\n");
   }
   if (ioPtr != NULL) ioPtr->getParameter("input_names", qData);

   if (nInputs <= 0 || nSamples <= 0 || nOutputs <= 0 || 
       outputID < 0 || outputID >= nOutputs)
   {
      printf("MOATAnalyzer ERROR: invalid arguments.\n");
      printf("    nInputs  = %d\n", nInputs);
      printf("    nOutputs = %d\n", nOutputs);
      printf("    nSamples = %d\n", nSamples);
      printf("    outputID = %d\n", outputID+1);
      exit(1);
   } 
   n1 = 0;
   for (iD = 0; iD < nSamples; iD++)
      if (S[iD] != 1 || Y[iD*nOutputs+outputID] == PSUADE_UNDEFINED) n1++; 
   if (n1 > 0) printf("MOATAnalyzer INFO: %d invalid data points.\n",n1);

   XSort = new double[nSamples];
   for (ii = 0; ii < nInputs; ii++)
   {
      for (iD = 0; iD < nSamples; iD++)
         XSort[iD] = X[iD*nInputs+ii];
      sortDbleList(nSamples, XSort);
      dtemp = XSort[nSamples-1] - XSort[0];
      dtemp2 = xUpper[ii] - xLower[ii];
      if (PABS(dtemp-dtemp2) > 1.0e-6)
      {
         printf("MOATAnalyzer WARNING: input and data range mismatch.\n");
         printf("    (This can be the result of applying constraints.)\n");
         printf("    diagnostics: \n");
         printf("    Input %3d: original vs new ranges = %e %e\n", ii+1,
                dtemp2, dtemp);
      }
   }
   delete [] XSort;

   printAsterisks(0);
   for (ii = 0; ii < 23; ii++) printf("*");
   printf(" MOAT Analysis ");
   for (ii = 0; ii < 23; ii++) printf("*");
   printf("\n");
   printDashes(0);

   constrPtr = new MOATConstraints();
   constrPtr->initialize(ioPtr);

   YY = new double[nSamples];
   YG = new double[nSamples];
   XG = new double[nSamples];
   for (iD = 0; iD < nSamples; iD++) YY[iD] = Y[nOutputs*iD+outputID];
   counts = new int[nInputs];
   for (ii = 0; ii < nInputs; ii++) counts[ii] = 0;
   means = new double[nInputs];
   modifiedMeans = new double[nInputs];
   stds = new double[nInputs];
   modifiedStds = new double[nInputs];
   for (ii = 0; ii < nInputs; ii++)
      means[ii] = modifiedMeans[ii] = stds[ii] = modifiedStds[ii] = 0.0;
   indexTrack = new int[nSamples];
   for (iD = 0; iD < nSamples; iD++) indexTrack[iD] = -1;
   Xbase = new double[nSamples];
   for (iD = 0; iD < nSamples; iD++) Xbase[iD] = 0.0;

   indexTrack[0] = -1;
   for (iD = 1; iD < nSamples; iD++)
   {
      if (printLevel > 3 && (iD % 10 == 0))
         printf("MOATAnalyzer: processing sample %d\n", iD+1);

      Xbase[iD] = 0.0;
      ytemp1 = YY[iD-1]; 
      ytemp2 = YY[iD]; 
      diffCnt = 0;
      for (ii = 0; ii < nInputs; ii++)
      {
         xtemp1 = X[(iD-1)*nInputs+ii]; 
         xtemp2 = X[iD*nInputs+ii]; 
         if (xtemp1 != xtemp2 && ytemp1 !=  PSUADE_UNDEFINED &&
             ytemp2 != PSUADE_UNDEFINED && S[iD-1] != 0 &&
             S[iD] != 0)
         {
            diffCnt++;
            diffIndex = ii;
         }
      }
      if (diffCnt == 1)
      {
         indexTrack[iD] = diffIndex;
         xtemp1 = X[(iD-1)*nInputs+diffIndex]; 
         xtemp2 = X[iD*nInputs+diffIndex]; 
         scale = constrPtr->getScale(&X[iD*nInputs],diffIndex,flag);
         if (flag == 1) scale = xUpper[diffIndex] - xLower[diffIndex];
         else if (printLevel > 3)
         {
            printf("MOATAnalyzer: sample group %d, ", iD+1);
            printf("input %d, scale = %e\n", diffIndex, scale);
         }
         YG[iD] = (ytemp2-ytemp1)/(xtemp2-xtemp1)*scale;
         if (xtemp2 > xtemp1) XG[iD] = xtemp2;
         else                 XG[iD] = xtemp1;
         counts[diffIndex]++;
         if (xtemp2 > xtemp1) Xbase[iD] = xtemp1;
         else                 Xbase[iD] = xtemp2;
      }
      else
      {
         YG[iD] = PSUADE_UNDEFINED;
         indexTrack[iD] = -1;
      }
   }

   if (nSamples / (nInputs+1) * (nInputs+1) == nSamples)
   {
      for (iD = 0; iD < nSamples; iD+=(nInputs+1))
         indexTrack[iD] = -1;
   }

   for (iD = 0; iD < nSamples; iD++)
   {
      if (YG[iD] != PSUADE_UNDEFINED)
      {
         index = indexTrack[iD];
         if (index >= 0)
         {
            means[index] += YG[iD];
            modifiedMeans[index] += PABS(YG[iD]);
         }
      }
   }
   for (ii = 0; ii < nInputs; ii++)
   {
      if (counts[ii] > 0)
      {
         means[ii] /= (double) counts[ii];
         modifiedMeans[ii] /= (double) counts[ii];
      }
      else 
      {
         printf("MOATAnalyzer analyze: zero data points for input %d\n",
                ii+1);
         means[ii] = 0.0;
         modifiedMeans[ii] = 0.0;
      }
   }
   for (iD = 0; iD < nSamples; iD++)
   {
      if (YG[iD] != PSUADE_UNDEFINED)
      {
         index = indexTrack[iD];
         if (index >= 0)
         {
            stds[index] += (YG[iD] - means[index])*(YG[iD] - means[index]);
            modifiedStds[index] += (PABS(YG[iD]) - modifiedMeans[index])*
                                   (PABS(YG[iD]) - modifiedMeans[index]);
         }
      }
   }
   for (ii = 0; ii < nInputs; ii++)
   {
      if (counts[ii] > 1)
      {
         stds[ii] /= (double) (counts[ii] - 1);
         modifiedStds[ii] /= (double) (counts[ii] - 1);
      }
      else
      {
         printf("MOATAnalyzer analyze: %d data points for input %d\n",
                counts[ii], ii+1);
         stds[ii] = 0.0;
         modifiedStds[ii] = 0.0;
      }
      if (stds[ii] < 0.0) stds[ii] = -sqrt(-stds[ii]);
      else                stds[ii] = sqrt(stds[ii]);
      if (modifiedStds[ii] < 0.0) modifiedStds[ii] = -sqrt(-modifiedStds[ii]);
      else                        modifiedStds[ii] = sqrt(modifiedStds[ii]);
   }
   n1 = 0;
   for (ii = 0; ii < nInputs; ii++) n1 += counts[ii];
   if (n1 <= 0)
   {
      delete [] counts;
      delete constrPtr;
      delete [] YY;
      delete [] YG;
      delete [] XG;
      delete [] means;
      delete [] modifiedMeans;
      delete [] modifiedStds;
      delete [] stds;
      delete [] indexTrack;
      delete [] Xbase;
      printf("MOATAnalyzer analyze INFO: probably not an MOAT sample.\n");
      return 1;
   }

   printDashes(0);
   for (ii = 0; ii < nInputs; ii++)
      printf("Input %3d (unmod. mean & std) = %12.4e %12.4e \n",
             ii+1, means[ii], stds[ii]);

   printAsterisks(0);
   for (ii = 0; ii < 23; ii++) printf("*");
   printf(" MOAT Analysis ");
   for (ii = 0; ii < 23; ii++) printf("*");
   printf("\n");
   printDashes(0);
   for (ii = 0; ii < nInputs; ii++)
      printf("Input %3d (mod. mean & std) = %12.4e %12.4e \n",
             ii+1, modifiedMeans[ii], stds[ii]);

   if (psAnalysisInteractive_ == 1 || psAnaExpertMode_ == 1)
      createScreenDiagramFile(nSamples, nInputs, YG, indexTrack, 
                              modifiedMeans,stds,outputID,qData.strArray_);
   if (psAnalysisInteractive_ == 1 || psAnaExpertMode_ == 1)
      createScatterFile(nSamples, nInputs,YG,Xbase,indexTrack,qData.strArray_);

   if (psAnalysisInteractive_ == 1 || psAnaExpertMode_ == 1)
      createBootstrapFile(nSamples,nInputs,YG,Xbase,indexTrack,qData.strArray_);

   for (ii = 0; ii < nInputs; ii++) means[ii] = ii;
   sortDbleList2(nInputs, modifiedMeans, means);
   printAsterisks(0);
   for (ii = 0; ii < 18; ii++) printf("*");
   printf(" MOAT Analysis (ordered) ");
   for (ii = 0; ii < 18; ii++) printf("*");
   printf("\n");
   printDashes(0);
   for (ii = nInputs-1; ii >= 0; ii--)
   {
      iD = (int) means[ii];
      itemp = counts[iD] - 1;
      printf("Input %3d (mu*, sigma, dof) = %12.4e %12.4e %d\n",
             iD+1, modifiedMeans[ii], stds[iD], itemp);
   }
   printAsterisks(0);
   if (printLevel > 1)
   {
      for (ii = 0; ii < 11; ii++) printf("-");
      printf(" MOAT Analysis (ordered) : +- 1 sigma  ");
      for (ii = 0; ii < 11; ii++) printf("-");
      printf("\n");
      printDashes(0);
      for (ii = nInputs-1; ii >= 0; ii--)
      {
         iD = (int) means[ii];
         if (counts[iD] > 1) dtemp = sqrt((double) counts[iD] - 1.0);
         else                dtemp = 1.0;
         printf("Input %3d bounds = %12.4e %12.4e\n",
                iD+1, modifiedMeans[ii]-modifiedStds[iD]/dtemp,
                modifiedMeans[ii]+modifiedStds[iD]/dtemp);
         if (printLevel > 2 && ii > 0)
         {
            iD2 = (int) means[ii-1];
            if (counts[iD2-1] > 1) dtemp2 = sqrt((double) counts[iD2-1] - 1.0);
            else                   dtemp2 = 1.0;
            if ((modifiedMeans[ii]-modifiedStds[iD]/dtemp) >
                (modifiedMeans[ii-1]+modifiedStds[iD2]/dtemp2))
               printf("=============> Input %3d is different from input %3d.\n",
                  iD+1, iD2+1);
         }
      }
      printAsterisks(0);
      for (ii = 0; ii < 11; ii++) printf("-");
      printf(" MOAT Analysis (ordered) : +- 2 sigma  ");
      for (ii = 0; ii < 11; ii++) printf("-");
      printf("\n");
      printDashes(0);
      for (ii = nInputs-1; ii >= 0; ii--)
      {
         iD = (int) means[ii];
         if (counts[iD] > 1) dtemp = sqrt((double) counts[iD] - 1.0);
         else                dtemp = 1.0;
         printf("Input %3d bounds = %12.4e %12.4e\n",
                iD+1, modifiedMeans[ii]-2.0*modifiedStds[iD]/dtemp,
                modifiedMeans[ii]+2.0*modifiedStds[iD]/dtemp);
         if (printLevel > 2 && ii > 0)
         {
            iD2 = (int) means[ii-1];
            if (counts[iD2-1] > 1) dtemp2 = sqrt((double) counts[iD2-1] - 1.0);
            else                   dtemp2 = 1.0;
            if ((modifiedMeans[ii]-2.0*modifiedStds[iD]/dtemp) >
                (modifiedMeans[ii-1]+2.0*modifiedStds[iD2]/dtemp2))
               printf("=============> Input %3d is different from input %3d.\n",
                  iD+1, iD2+1);
         }
      }
   }

   if (psAnaExpertMode_ == 1)
   {
      printAsterisks(0);
      sprintf(pString, "Perform MOAT interaction study ? (y or n) ");
      getString(pString, winput);
      if (winput[0] == 'y')
      {
         printAsterisks(0);
         for (ii = 0; ii < 17; ii++) printf("*");
         printf(" MOAT Interaction Analysis ");
         for (ii = 0; ii < 17; ii++) printf("*");
         printf("\n");
         printDashes(0);
         printf("** No data for an input means no interaction info.\n");
         printf("** Small stdev means (crudely) little interaction.\n");
         nPaths = nSamples / (nInputs + 1);
         XSort = new double[nPaths];
         YSort = new double[nPaths];
         for (ii = 0; ii < nInputs; ii++)
         {
            actualPaths = 0;
            for (iD = 0; iD < nSamples; iD++)
            {
               if (indexTrack[iD] == ii && YG[iD] != PSUADE_UNDEFINED)
               { 
                  XSort[actualPaths] = XG[iD]; 
                  YSort[actualPaths] = PABS(YG[iD]); 
                  actualPaths++;
               }
            }
            sortDbleList2(actualPaths, XSort, YSort);
            ip = 0; 
            iaCnt = 0;
            dstdev = 0.0;
            while (ip < actualPaths)
            {
               n1 = ip;
               for (ip = n1+1; ip < actualPaths; ip++)
                  if (XSort[ip] != XSort[ip-1]) break;
               n2 = ip;
               if ((n2 - n1) >= 2)
               {
                  dsum = 0.0;
                  for (iD = n1; iD < n2; iD++) dsum += YSort[iD];
                  dsum /= (double) (n2-n1);
                  dtemp = 0.0;
                  for (iD = n1; iD < n2; iD++) 
                     dtemp += (YSort[iD] - dsum) * (YSort[iD] - dsum);
                  dtemp /= (double) (n2-n1-1);
                  dtemp = sqrt(dtemp);
                  dstdev += dtemp;
                  iaCnt++;
                  if (printLevel > 3)
                     printf("Input %3d stdev at %12.4e = %12.4e (%d) \n",
                            ii+1, XSort[n1], dtemp, n2-n1);
               }
               else
                  if (printLevel > 3)
                     printf("Input %3d stdev at %12.4e = not available (%d).\n",
                            ii+1, XSort[n1], n2-n1);
            }
            if (iaCnt > 0)
               printf("Input %3d average interaction measure (std dev) = %11.3e\n",
                      ii+1, dstdev/iaCnt);
         }
         delete [] XSort;
         delete [] YSort;
      }
   }

   if (psAnaExpertMode_ == 1)
   {
      printAsterisks(0);
      sprintf(pString, "Perform hypothesis tests ? (y or n) ");
      getString(pString, winput);
      if (winput[0] == 'y')
      {
         printAsterisks(0);
         for (ii = 0; ii < 19; ii++) printf("*");
         printf("MOAT Hypothesis Tests");
         for (ii = 0; ii < 19; ii++) printf("*");
         printf("\n");
         printDashes(0);
         printf("* This consists of 2 tests: \n");
         printf("* (1) bootstrap confidence interval test (95 %%) \n");
         printf("* (2) binomial test (90 %%) \n");
         printf("* If either test indicates significance, the corresponding\n");
         printf("*   input will be declared significant.\n");
         printDashes(0);
         for (ii = 0; ii < nInputs; ii++)
         {
            if (modifiedMeans[ii] > binMax) binMax = modifiedMeans[ii];
            if (modifiedMeans[ii] < binMin) binMin = modifiedMeans[ii];
         }
         sprintf(pString,"Enter the significance threshold (min,max = %e %e): ",
                 binMin, binMax);
         thresh = binMin - 1.0;
         while (thresh < binMin || thresh > binMax)
            thresh = getDouble(pString);
         YB = new double[nSamples/(nInputs+1)*100];
         bData.nOutputs_ = 1;
         bData.outputID_ = 0;
         bData.sampleOutputs_ = YB;
         bData.analysisThreshold_ = thresh;
         bData.printLevel_ = 0;
         for (ii = 0; ii < nInputs; ii++)
         {
            printf(">>>>>> MOAT hypothesis tests for INPUT %3d : \n", ii+1);
            sigFlag = 0;
            actualPaths = 0;
            for (iD = 0; iD < nSamples; iD++)
            {
               if (YG[iD] != PSUADE_UNDEFINED)
               {
                  index = indexTrack[iD];
                  if (index == ii) YB[actualPaths++] = PABS(YG[iD]);
               }
            }
            bData.nSamples_ = actualPaths;
            dtemp = bsAnalyzer.analyze(bData);
            if (dtemp > thresh) sigFlag = 1; 
            if (printLevel > 1)
            {
               printf("   ---- Bootstrap confidence interval = [0, %e] (>? %e)",
                      dtemp, thresh);
               if (sigFlag == 1) printf(" **");
               printf("\n");
            } 
            jj = actualPaths;
            for (iD = 0; iD < actualPaths; iD++)
            {
               if (YB[iD] > 1.5*thresh)   
               {
                  dtemp = YB[iD];
                  while (dtemp > 1.5*thresh)
                  {
                     YB[jj++] = 1.5 * thresh;
                     dtemp -= (1.5 * thresh);
                     if (jj >= (nSamples/nInputs+1)*100)
                     {
                        printf("MOATAnalyzer::analyze ERROR- binomial test.\n");
                        exit(1);
                     }
                  }
               }
            }
            actualPaths = jj;
            bData.nSamples_ = actualPaths;
            dtemp = binAnalyzer.analyze(bData);
            if (dtemp > 0.1) sigFlag++;
            if (printLevel > 1)
               printf("   ---- Binomial test ERROR = %4.1f %% if declared insignificant.\n",
                      100.0*dtemp);
            if (sigFlag > 0) 
                 printf("<<<<<< INPUT %3d IS SIGNIFICANT (%d) *****\n", ii+1,sigFlag); 
            else printf("<<<<<< INPUT %3d IS NOT SIGNIFICANT\n", ii+1); 
            
         }
         bData.sampleOutputs_ = NULL;
      }
   }

   delete [] counts;
   delete constrPtr;
   delete [] YY;
   delete [] YG;
   delete [] XG;
   delete [] means;
   delete [] modifiedMeans;
   delete [] modifiedStds;
   delete [] stds;
   delete [] indexTrack;
   delete [] Xbase;
   return 0.0;
}

// ************************************************************************
// create Morris diagram matlab/scilab file
// ------------------------------------------------------------------------
int MOATAnalyzer::createScreenDiagramFile(int nSamples, int nInputs, 
                        double *Y, int *indices, double *modifiedMeans,
                        double *stds, int outputID, char **iNames)
{
   int  iD, ii, index, cnt;
   char winput[500], moatFile[500], pString[500];
   FILE *fp;

   printAsterisks(0);
   sprintf(pString,"Create screening diagram? (y or n) ");
   getString(pString, winput);
   if (winput[0] != 'y') return 0;
   sprintf(pString, "matlab/scilab screening diagram file name (no extension): ");
   getString(pString, moatFile);
   cnt = strlen(moatFile);
   if (cnt > 500)
   {
      printf("ERROR: file name too long.\n");
      exit(1);
   }
   moatFile[cnt-1] = '.';
   if (psPlotTool_ == 1)
   {
      moatFile[cnt] = 's';
      moatFile[cnt+1] = 'c';
      moatFile[cnt+2] = 'i';
      moatFile[cnt+3] = '\0';
   }
   else
   {
      moatFile[cnt] = 'm';
      moatFile[cnt+1] = '\0';
   }
   fp = fopen(moatFile, "w");

   if (fp == NULL)
   {
      printf("MOATAnalyzer: cannot open MOAT plot file %s\n",moatFile);
      return 0;
   }

   if (psPlotTool_ == 1)
   {
      fprintf(fp, "// Morris one-at-a-time screening plots\n");
      fprintf(fp, "// Z contains the gradient data\n");
      fprintf(fp, "// C contains the number of data for each input\n");
      fprintf(fp, "Z = zeros(%d,%d);\n",nInputs,nSamples/(nInputs+1));
      fprintf(fp, "C = zeros(%d,1);\n",nInputs);
      for (iD = 0; iD < nSamples; iD++)
      {
         if (Y[iD] != PSUADE_UNDEFINED)
         {
            index = indices[iD];
            if (index >= 0)
            {
               fprintf(fp, "C(%d) = C(%d) + 1;\n", index+1, index+1);
               fprintf(fp, "Z(%d,C(%d)) = %24.16e;\n",index+1,index+1,Y[iD]);
            }
         }
      }
      fprintf(fp, "// compute max and min for axis \n");
      fprintf(fp, "// XM : for computing mean \n");
      fprintf(fp, "// XX : for computing modified mean \n");
      fprintf(fp, "// VV : counts for each input \n");
      fprintf(fp, "// YY : standard deviation for each input \n");
      fprintf(fp, "XX = zeros(%d,1);\n",nInputs);
      fprintf(fp, "XM = zeros(%d,1);\n",nInputs);
      fprintf(fp, "YY = zeros(%d,1);\n",nInputs);
      fprintf(fp, "VV = zeros(%d,1);\n",nInputs);
      fprintf(fp, "for jj = 1 : %d\n", nInputs);
      fprintf(fp, "   for kk = 1 : C(jj)\n");
      fprintf(fp, "      if (VV(jj) <= C(jj))\n");
      fprintf(fp, "         XM(jj) = XM(jj) + Z(jj,kk);\n");
      fprintf(fp, "         XX(jj) = XX(jj) + abs(Z(jj,kk));\n");
      fprintf(fp, "         VV(jj) = VV(jj) + 1;\n");
      fprintf(fp, "      end;\n");
      fprintf(fp, "   end;\n");
      fprintf(fp, "   if (VV(jj) > 0)\n");
      fprintf(fp, "      XM(jj) = XM(jj) / VV(jj);\n");
      fprintf(fp, "      XX(jj) = XX(jj) / VV(jj);\n");
      fprintf(fp, "   end;\n");
      fprintf(fp, "   for kk = 1 : C(jj)\n");
      fprintf(fp, "      if (VV(jj) <= C(jj))\n");
      fprintf(fp, "         YY(jj) = YY(jj) + (Z(jj,kk)-XM(jj))^2;\n");
      fprintf(fp, "      end;\n");
      fprintf(fp, "   end;\n");
      fprintf(fp, "   if (VV(jj) > 1)\n");
      fprintf(fp, "      YY(jj) = sqrt(YY(jj) / (VV(jj)-1));\n");
      fprintf(fp, "   end;\n");
      fprintf(fp, "end;\n");
      fprintf(fp, "// plot sequence of Morris plot\n");
      fprintf(fp, "first = floor(max(C) / 4)\n");
      fprintf(fp, "last  = max(C);\n");
      fprintf(fp, "if (first < 2)\n");
      fprintf(fp, "   first = 2;\n");
      fprintf(fp, "end;\n");
      fprintf(fp, "if (last < first)\n");
      fprintf(fp, "   last = first;\n");
      fprintf(fp, "end;\n");
      fprintf(fp, "inc = floor((last - first + 1)/3);\n");
      fprintf(fp, "if (inc < 1)\n");
      fprintf(fp, "   inc = 1;\n");
      fprintf(fp, "end;\n");
      fprintf(fp, "list = [first : inc : last];\n");
      fprintf(fp, "list = unique(list);\n");
      fprintf(fp, "count = length(list);\n");

      fprintf(fp, "clf\n");
      fprintf(fp, "for mm = first : last\n");
      fprintf(fp, "   ii = list(mm);\n");
      fprintf(fp, "   XX = zeros(%d,1);\n",nInputs);
      fprintf(fp, "   XM = zeros(%d,1);\n",nInputs);
      fprintf(fp, "   YY = zeros(%d,1);\n",nInputs);
      fprintf(fp, "   VV = zeros(%d,1);\n",nInputs);
      fprintf(fp, "   for jj = 1 : %d\n", nInputs);
      fprintf(fp, "      for kk = 1 : ii\n");
      fprintf(fp, "         if (VV(jj) <= C(jj))\n");
      fprintf(fp, "            XM(jj) = XX(jj) + Z(jj,kk);\n");
      fprintf(fp, "            XX(jj) = XX(jj) + abs(Z(jj,kk));\n");
      fprintf(fp, "            VV(jj) = VV(jj) + 1;\n");
      fprintf(fp, "         end;\n");
      fprintf(fp, "      end;\n");
      fprintf(fp, "      if (VV(jj) > 0)\n");
      fprintf(fp, "         XM(jj) = XM(jj) / VV(jj);\n");
      fprintf(fp, "         XX(jj) = XX(jj) / VV(jj);\n");
      fprintf(fp, "      end;\n");
      fprintf(fp, "      for kk = 1 : ii\n");
      fprintf(fp, "         if (VV(jj) <= C(jj))\n");
      fprintf(fp, "            YY(jj) = YY(jj) + (Z(jj,kk)-XM(jj))^2;\n");
      fprintf(fp, "         end;\n");
      fprintf(fp, "      end;\n");
      fprintf(fp, "      if (VV(jj) > 1)\n");
      fprintf(fp, "         YY(jj) = sqrt(YY(jj) / (VV(jj)-1));\n");
      fprintf(fp, "      end;\n");
      fprintf(fp, "   end;\n");
      fprintf(fp, "scf(1);\n");
      fprintf(fp, "   subplot(2,2,ii-first+1);\n");
      fprintf(fp, "   bar(XX, 0.8);\n");
      fwritePlotAxesNoGrid(fp);
      if (iNames != NULL)
      {
         fprintf(fp, "set(gca,'XTick',[1:%d],'XTickLabel',{", nInputs);
         for (ii = 0; ii < nInputs-1; ii++) fprintf(fp,"'%s',", iNames[ii]);
         fprintf(fp,"'%s'})\n", iNames[nInputs-1]);
      }
      else fwritePlotXLabel(fp, "Input Numbers");
      fwritePlotYLabel(fp, "Modified Means (of gradients)");
      fprintf(fp,"pStr=sprintf('Modified Means with %%d replications',ii);\n");
      fprintf(fp,"title(pStr);\n");
      fprintf(fp, "a = gca();\n");
      fprintf(fp, "a.title.text = pStr;\n");
      fprintf(fp, "a.title.font_size = 3;\n");
      fprintf(fp, "a.title.font_style = 4;\n");
      fprintf(fp, "scf(2);\n");
      fprintf(fp, "   subplot(2,2,ii-first+1);\n");
      fprintf(fp, "   bar(YY, 0.8);\n");
      fwritePlotAxesNoGrid(fp);
      fwritePlotXLabel(fp, "Input Numbers");
      fwritePlotYLabel(fp, "Std Devs");
      fprintf(fp,"pStr=sprintf('Std Devs with %%d replications',ii);\n");
      fprintf(fp,"title(pStr);\n");
      fprintf(fp, "a = gca();\n");
      fprintf(fp, "a.title.text = pStr;\n");
      fprintf(fp, "a.title.font_size = 3;\n");
      fprintf(fp, "a.title.font_style = 4;\n");
      fprintf(fp, "end;\n");

      fprintf(fp, "Y = [\n");
      for (ii = 0; ii < nInputs; ii++) fprintf(fp, "%24.16e\n", stds[ii]);
      fprintf(fp, "];\n");
      fprintf(fp, "X = [\n");
      for (ii = 0; ii < nInputs; ii++)
         fprintf(fp, "%24.16e\n", modifiedMeans[ii]);
      fprintf(fp, "];\n");
      fprintf(fp, "scf(3);\n");
      fprintf(fp, "plot(X,Y,'*');\n");
      fprintf(fp, "xstring(X*1.02,Y,{");
      if (iNames != NULL)
      {
         for (ii = 0; ii < nInputs-1; ii++) fprintf(fp, "'%s',",iNames[ii]);
         fprintf(fp, "'%s'},'FontWeight','bold','FontSize',12)\n",
                 iNames[nInputs-1]);
      }
      else
      {
         for (ii = 0; ii < nInputs-1; ii++) fprintf(fp, "'X%d',",ii+1);
         fprintf(fp, "'X%d'},'FontWeight','bold','FontSize',12)\n",nInputs);
      }
      fprintf(fp, "Xmin = min(X) - 0.1 * (max(X) - min(X));\n");
      fprintf(fp, "Ymin = min(Y) - 0.1 * (max(Y) - min(Y));\n");
      fprintf(fp, "Xmax = max(X) + 0.1 * (max(X) - min(X));\n");
      fprintf(fp, "Ymax = max(Y) + 0.1 * (max(Y) - min(Y));\n");
      fprintf(fp, "a=gca();\n");
      fprintf(fp, "a.data_bounds=[Xmin, Ymin; Xmax, Ymax];\n");
      fwritePlotXLabel(fp, "Modified Means (of gradients)");
      fwritePlotYLabel(fp, "Std Deviations (of gradients)");
      fwritePlotAxes(fp);
      sprintf(pString, "Modified Morris Diagram for Output %d", outputID+1);
      fwritePlotTitle(fp, pString);
      printf("MOAT screening diagram scilab file = %s\n", moatFile);
   }
   else
   {
      fprintf(fp, "%% Morris one-at-a-time screening plots\n");
      fprintf(fp, "%% Z contains the gradient data\n");
      fprintf(fp, "%% C contains the number of data for each input\n");
      fprintf(fp, "%% turn plotAll on or off\n");
      fprintf(fp, "plotAll = 1;\n");
      fprintf(fp, "Z = zeros(%d,%d);\n",nInputs,nSamples/(nInputs+1));
      fprintf(fp, "C = zeros(%d,1);\n",nInputs);
      for (iD = 0; iD < nSamples; iD++)
      {
         if (Y[iD] != PSUADE_UNDEFINED)
         {
            index = indices[iD];
            if (index >= 0)
            {
               fprintf(fp, "C(%d) = C(%d) + 1;\n", index+1, index+1);
               fprintf(fp, "Z(%d,C(%d)) = %24.16e;\n",index+1,index+1,Y[iD]);
            }
         }
      }
      fprintf(fp, "%% compute max and min for axis \n");
      fprintf(fp, "%% XM : for computing mean \n");
      fprintf(fp, "%% XX : for computing modified mean \n");
      fprintf(fp, "%% VV : counts for each input \n");
      fprintf(fp, "%% YY : standard deviation for each input \n");
      fprintf(fp, "XX = zeros(%d,1);\n",nInputs);
      fprintf(fp, "XM = zeros(%d,1);\n",nInputs);
      fprintf(fp, "YY = zeros(%d,1);\n",nInputs);
      fprintf(fp, "for jj = 1 : %d\n", nInputs);
      fprintf(fp, "   VV = zeros(%d,1);\n",nInputs);
      fprintf(fp, "   for kk = 1 : C(jj)\n");
      fprintf(fp, "      if (VV(jj) <= C(jj))\n");
      fprintf(fp, "         XM(jj) = XM(jj) + Z(jj,kk);\n");
      fprintf(fp, "         XX(jj) = XX(jj) + abs(Z(jj,kk));\n");
      fprintf(fp, "         VV(jj) = VV(jj) + 1;\n");
      fprintf(fp, "      end;\n");
      fprintf(fp, "   end;\n");
      fprintf(fp, "   if (VV(jj) > 0)\n");
      fprintf(fp, "      XM(jj) = XM(jj) / VV(jj);\n");
      fprintf(fp, "      XX(jj) = XX(jj) / VV(jj);\n");
      fprintf(fp, "   end;\n");
      fprintf(fp, "   VV = zeros(%d,1);\n",nInputs);
      fprintf(fp, "   for kk = 1 : C(jj)\n");
      fprintf(fp, "      if (VV(jj) <= C(jj))\n");
      fprintf(fp, "         YY(jj) = YY(jj) + (Z(jj,kk)-XM(jj))^2;\n");
      fprintf(fp, "         VV(jj) = VV(jj) + 1;\n");
      fprintf(fp, "      end;\n");
      fprintf(fp, "   end;\n");
      fprintf(fp, "   if (VV(jj) > 1)\n");
      fprintf(fp, "      YY(jj) = sqrt(YY(jj) / (VV(jj)-1));\n");
      fprintf(fp, "   end;\n");
      fprintf(fp, "end;\n");
      fprintf(fp, "%% plot sequence of Morris plot\n");
      fprintf(fp, "first = floor(max(C) / 4)\n");
      fprintf(fp, "last  = max(C);\n");
      fprintf(fp, "if (first < 2)\n");
      fprintf(fp, "   first = 2;\n");
      fprintf(fp, "end;\n");
      fprintf(fp, "if (last < first)\n");
      fprintf(fp, "   last = first;\n");
      fprintf(fp, "end;\n");
      fprintf(fp, "inc = floor((last - first + 1)/3);\n");
      fprintf(fp, "if (inc < 1)\n");
      fprintf(fp, "   inc = 1;\n");
      fprintf(fp, "end;\n");
      fprintf(fp, "list = [first : inc : last];\n");
      fprintf(fp, "list = unique(list);\n");
      fprintf(fp, "count = length(list);\n");
      fprintf(fp, "clf\n");
      fprintf(fp, "for mm = 1 : count\n");
      fprintf(fp, "   ii = list(mm);\n");
      fprintf(fp, "   XX = zeros(%d,1);\n",nInputs);
      fprintf(fp, "   XM = zeros(%d,1);\n",nInputs);
      fprintf(fp, "   YY = zeros(%d,1);\n",nInputs);
      fprintf(fp, "   VV = zeros(%d,1);\n",nInputs);
      fprintf(fp, "   for jj = 1 : %d\n", nInputs);
      fprintf(fp, "      cnt = 0;\n");
      fprintf(fp, "      for kk = 1 : ii\n");
      fprintf(fp, "         if (cnt <= C(jj))\n");
      fprintf(fp, "            XM(jj) = XM(jj) + Z(jj,kk);\n");
      fprintf(fp, "            XX(jj) = XX(jj) + abs(Z(jj,kk));\n");
      fprintf(fp, "            cnt = cnt + 1;\n");
      fprintf(fp, "         end;\n");
      fprintf(fp, "      end;\n");
      fprintf(fp, "      if (cnt > 0)\n");
      fprintf(fp, "         XM(jj) = XM(jj) / cnt;\n");
      fprintf(fp, "         XX(jj) = XX(jj) / cnt;\n");
      fprintf(fp, "      end;\n");
      fprintf(fp, "      cnt = 0;\n");
      fprintf(fp, "      for kk = 1 : ii\n");
      fprintf(fp, "         if (cnt <= C(jj))\n");
      fprintf(fp, "            YY(jj) = YY(jj) + (Z(jj,kk)-XM(jj))^2;\n");
      fprintf(fp, "            cnt = cnt + 1;\n");
      fprintf(fp, "         end;\n");
      fprintf(fp, "      end;\n");
      fprintf(fp, "      if (cnt > 1)\n");
      fprintf(fp, "         YY(jj) = sqrt(YY(jj) / (cnt-1));\n");
      fprintf(fp, "      end;\n");
      fprintf(fp, "   end;\n");
      fprintf(fp, "   figure(1)\n");
      fprintf(fp, "   subplot(2,2,mm) \n");
      fprintf(fp, "   bar(XX,0.8)\n");
      fwritePlotAxes(fp);
      if (iNames != NULL)
      {
         fprintf(fp, "set(gca,'XTick',[1:%d],'XTickLabel',{", nInputs);
         for (ii = 0; ii < nInputs-1; ii++) fprintf(fp,"'%s',", iNames[ii]);
         fprintf(fp,"'%s'})\n", iNames[nInputs-1]);
      }
      else fwritePlotXLabel(fp, "Input Numbers");
      fwritePlotYLabel(fp, "Modified Means (of gradients)");
      fprintf(fp,"    pStr=sprintf('%%d replications',list(mm));\n");
      fprintf(fp,"    title(pStr);\n");
      fprintf(fp,"    Xmin = min(XX) - (max(XX) - min(XX)) * 0.1;\n");
      fprintf(fp,"    Xmax = max(XX) + (max(XX) - min(XX)) * 0.1;\n");
      fprintf(fp,"    axis([0 %d Xmin Xmax])\n", nInputs+1);
      fprintf(fp, "   if (plotAll == 1)\n");
      fprintf(fp, "   figure(2)\n");
      fprintf(fp, "   subplot(2,2,mm) \n");
      fprintf(fp, "   bar(YY,0.8)\n");
      fprintf(fp,"    Ymin = min(YY) - (max(YY) - min(YY)) * 0.1;\n");
      fprintf(fp,"    Ymax = max(YY) + (max(YY) - min(YY)) * 0.1;\n");
      fprintf(fp,"    axis([0 %d Ymin Ymax])\n", nInputs+1);
      fwritePlotAxes(fp);
      if (iNames != NULL)
      {
         fprintf(fp, "set(gca,'XTick',[1:%d],'XTickLabel',{", nInputs);
         for (ii = 0; ii < nInputs-1; ii++) fprintf(fp,"'%s',", iNames[ii]);
         fprintf(fp,"'%s'})\n", iNames[nInputs-1]);
      }
      else fwritePlotXLabel(fp, "Input Numbers");
      fwritePlotYLabel(fp, "Std. Devs.");
      fprintf(fp,"    pStr=sprintf('%%d replications',list(mm));\n");
      fprintf(fp,"    title(pStr);\n");
      fprintf(fp, "   end;\n");
      fprintf(fp, "end;\n");

      fprintf(fp, "if (plotAll == 1)\n");
      fprintf(fp, "figure(3)\n");
      fprintf(fp, "Y = [\n");
      for (ii = 0; ii < nInputs; ii++) fprintf(fp, "%24.16e\n", stds[ii]);
      fprintf(fp, "];\n");
      fprintf(fp, "X = [\n");
      for (ii = 0; ii < nInputs; ii++)
         fprintf(fp, "%24.16e\n", modifiedMeans[ii]);
      fprintf(fp, "];\n");
      fprintf(fp, "plot(X,Y,'*','MarkerSize',12)\n");
      fprintf(fp, "text(X*1.01,Y,{");
      if (iNames != NULL)
      {
         for (ii = 0; ii < nInputs-1; ii++) fprintf(fp, "'%s',",iNames[ii]);
         fprintf(fp, "'%s'},'FontWeight','bold','FontSize',12)\n",
                 iNames[nInputs-1]);
      }
      else
      {
         for (ii = 0; ii < nInputs-1; ii++) fprintf(fp, "'X%d',",ii+1);
         fprintf(fp, "'X%d'},'FontWeight','bold','FontSize',12)\n",nInputs);
      }
      fprintf(fp, "Xmin = min(X) - (max(X) - min(X)) * 0.1;\n");
      fprintf(fp, "Ymin = min(Y) - (max(Y) - min(Y)) * 0.1;\n");
      fprintf(fp, "Xmax = max(X) + (max(X) - min(X)) * 0.1;\n");
      fprintf(fp, "Ymax = max(Y) + (max(Y) - min(Y)) * 0.1;\n");
      fprintf(fp, "axis([Xmin Xmax Ymin Ymax])\n");
      fwritePlotAxes(fp);
      fwritePlotXLabel(fp, "Modified Means (of gradients)");
      fwritePlotYLabel(fp, "Std Deviations (of gradients)");
      sprintf(pString, "Modified Morris Diagram for Output %d", outputID+1);
      fwritePlotTitle(fp, pString);
      fprintf(fp, "end;\n");
      printf("MOAT screening diagram matlab file = %s\n", moatFile);
   }
   fclose(fp); 
   printAsterisks(0);
   return 0;
}

// ************************************************************************
// create scatter plot matlab file
// ------------------------------------------------------------------------
int MOATAnalyzer::createScatterFile(int nSamples, int nInputs, double *Y, 
                                    double *X, int *indices, char **iNames)
{
   int  iD, cnt, ii;
   char winput[500], scatterFile[500], pString[500];
   FILE *fp;

   if (psPlotTool_ == 1)
   {
      printf("MOAT: scatter plot not availabe with scilab.\n");
      return 0;
   }
   printAsterisks(0);
   sprintf(pString, "Create scatter plot ? (y or n) ");
   getString(pString, winput);
   if (winput[0] != 'y')
   {
      return 0;
   }
   sprintf(pString,"Enter matlab/scilab scatter plot file name (no extension): ");
   getString(pString, scatterFile);
   cnt = strlen(scatterFile);
   if (cnt > 500)
   {
      printf("ERROR: file name too long.\n");
      exit(1);
   }
   scatterFile[cnt-1] = '.';
   if (psPlotTool_ == 1)
   {
      scatterFile[cnt] = 's';
      scatterFile[cnt+1] = 'c';
      scatterFile[cnt+2] = 'i';
      scatterFile[cnt+3] = '\0';
   }
   else
   {
      scatterFile[cnt] = 'm';
      scatterFile[cnt+1] = '\0';
   }

   fp = fopen(scatterFile, "w");
   if (fp == NULL)
   {
      printf("MOATAnalyzer: cannot open scatterplot file %s\n",scatterFile);
      return 0;
   }

   if (psPlotTool_ == 1)
   {
      fprintf(fp, "// This file contains individual gradient info.\n");
      fprintf(fp, "// The gradients are normalized for the range.\n");
   }
   else
   {
      fprintf(fp, "%% This file contains individual gradient info.\n");
      fprintf(fp, "%% The gradients are normalized for the range.\n");
   }
   fprintf(fp, "clf\n");
   fprintf(fp, "A = [\n");
   for (iD = 1; iD < nSamples; iD++)
   {
      if (Y[iD] != PSUADE_UNDEFINED)
         fprintf(fp,"%4d %24.16e %24.16e %d\n",indices[iD]+1,Y[iD],X[iD],iD+1);
   }
   fprintf(fp, "];\n");

   fprintf(fp, "hold off\n");
   fprintf(fp, "clf\n");
   fprintf(fp, "hold on\n");
   fprintf(fp, "ymin = min(A(:,2));\n");
   fprintf(fp, "ymax = max(A(:,2));\n");
   fprintf(fp, "n = %d;\n", nInputs);
   fprintf(fp, "for ii = 1 : n\n");
   fprintf(fp, "  inds = find(A(:,1)==ii);\n");
   fprintf(fp, "  leng = length(inds);\n");
   fprintf(fp, "  if (leng > 0)\n");
   fprintf(fp, "    [AA,II] = sort(A(inds,3));\n");
   fprintf(fp, "    xx = A(inds(II(1)),1);\n");
   fprintf(fp, "    yy = A(inds(II(1)),2);\n");
   fprintf(fp, "    plot(xx,yy,'ro','MarkerSize',11,'MarkerFaceColor','r',");
   fprintf(fp, "'MarkerEdgeColor','k')\n");
   fprintf(fp, "    for jj = 2 : leng\n");
   fprintf(fp, "      x1 = A(inds(II(jj)),3);\n");
   fprintf(fp, "      x2 = A(inds(II(jj-1)),3);\n");
   fprintf(fp, "      if x1 ~= x2\n");
   fprintf(fp, "        break;\n");
   fprintf(fp, "      else\n");
   fprintf(fp, "        xx = A(inds(II(jj)),1);\n");
   fprintf(fp, "        yy = A(inds(II(jj)),2);\n");
   fprintf(fp, "        plot(xx,yy,'ro','MarkerSize',11,'MarkerFaceColor','r',");
   fprintf(fp, "'MarkerEdgeColor','k')\n");
   fprintf(fp, "        if (jj == leng)\n");
   fprintf(fp, "          jj = jj + 1;\n");
   fprintf(fp, "        end;\n");
   fprintf(fp, "      end;\n");
   fprintf(fp, "    end;\n");
   fprintf(fp, "    if (jj <= leng)\n");
   fprintf(fp, "      xx = A(inds(II(jj)),1);\n");
   fprintf(fp, "      yy = A(inds(II(jj)),2);\n");
   fprintf(fp, "      plot(xx+0.1,yy,'go','MarkerSize',11,'MarkerFaceColor','g',");
   fprintf(fp, "'MarkerEdgeColor','k')\n");
   fprintf(fp, "    end;\n");
   fprintf(fp, "    for kk = jj+1 : leng\n");
   fprintf(fp, "      x1 = A(inds(II(kk)),3);\n");
   fprintf(fp, "      x2 = A(inds(II(kk-1)),3);\n");
   fprintf(fp, "      if x1 ~= x2\n");
   fprintf(fp, "        break;\n");
   fprintf(fp, "      else\n");
   fprintf(fp, "        xx = A(inds(II(kk)),1);\n");
   fprintf(fp, "        yy = A(inds(II(kk)),2);\n");
   fprintf(fp, "        plot(xx+0.1,yy,'go','MarkerSize',11,'MarkerFaceColor','g',");
   fprintf(fp, "'MarkerEdgeColor','k')\n");
   fprintf(fp, "        if (kk == leng)\n");
   fprintf(fp, "          kk = kk + 1;\n");
   fprintf(fp, "        end;\n");
   fprintf(fp, "      end;\n");
   fprintf(fp, "    end;\n");
   fprintf(fp, "    if (kk <= leng)\n");
   fprintf(fp, "      xx = A(inds(II(kk)),1);\n");
   fprintf(fp, "      yy = A(inds(II(kk)),2);\n");
   fprintf(fp, "      plot(xx+0.2,yy,'bo','MarkerSize',11,'MarkerFaceColor','b',");
   fprintf(fp, "'MarkerEdgeColor','k')\n");
   fprintf(fp, "    end;\n");
   fprintf(fp, "    for jj = kk+1 : leng\n");
   fprintf(fp, "      x1 = A(inds(II(jj)),3);\n");
   fprintf(fp, "      x2 = A(inds(II(jj-1)),3);\n");
   fprintf(fp, "      if x1 ~= x2\n");
   fprintf(fp, "        break;\n");
   fprintf(fp, "      else\n");
   fprintf(fp, "        xx = A(inds(II(jj)),1);\n");
   fprintf(fp, "        yy = A(inds(II(jj)),2);\n");
   fprintf(fp, "        plot(xx+0.2,yy,'bo','MarkerSize',11,'MarkerFaceColor','b',");
   fprintf(fp, "'MarkerEdgeColor','k')\n");
   fprintf(fp, "        if (jj == leng)\n");
   fprintf(fp, "          jj = jj + 1;\n");
   fprintf(fp, "        end;\n");
   fprintf(fp, "      end;\n");
   fprintf(fp, "    end;\n");
   fprintf(fp, "    if (jj <= leng)\n");
   fprintf(fp, "      xx = A(inds(II(jj)),1);\n");
   fprintf(fp, "      yy = A(inds(II(jj)),2);\n");
   fprintf(fp, "      plot(xx+0.3,yy,'mo','MarkerSize',11,'MarkerFaceColor','m',");
   fprintf(fp, "'MarkerEdgeColor','k')\n");
   fprintf(fp, "    end;\n");
   fprintf(fp, "    for kk = jj+1 : leng\n");
   fprintf(fp, "      x1 = A(inds(II(kk)),3);\n");
   fprintf(fp, "      x2 = A(inds(II(kk-1)),3);\n");
   fprintf(fp, "      if x1 ~= x2\n");
   fprintf(fp, "        break;\n");
   fprintf(fp, "      else\n");
   fprintf(fp, "        xx = A(inds(II(kk)),1);\n");
   fprintf(fp, "        yy = A(inds(II(kk)),2);\n");
   fprintf(fp, "        plot(xx+0.3,yy,'mo','MarkerSize',11,'MarkerFaceColor','m',");
   fprintf(fp, "'MarkerEdgeColor','k')\n");
   fprintf(fp, "        if (kk == leng)\n");
   fprintf(fp, "          kk = kk + 1;\n");
   fprintf(fp, "        end;\n");
   fprintf(fp, "      end;\n");
   fprintf(fp, "    end;\n");
   fprintf(fp, "    if (kk <= leng)\n");
   fprintf(fp, "      for jj = kk : leng\n");
   fprintf(fp, "        xx = A(inds(II(jj)),1);\n");
   fprintf(fp, "        yy = A(inds(II(jj)),2);\n");
   fprintf(fp, "        plot(xx+0.4,yy,'co','MarkerSize',11,'MarkerFaceColor','c',");
   fprintf(fp, "'MarkerEdgeColor','k')\n");
   fprintf(fp, "      end;\n");
   fprintf(fp, "    end;\n");
   fprintf(fp, "  end;\n");
   fprintf(fp, "end;\n");
   fprintf(fp, "for ii = 1 : n\n");
   fprintf(fp, "   inds = find(A(:,1)==ii);\n");
   fprintf(fp, "   if length(inds) > 0\n");
   fprintf(fp, "      asum = sum(A(inds,2))/length(inds);\n");
   fprintf(fp, "      plot(ii,asum,'k*','MarkerSize',13);\n");
   fprintf(fp, "   end;\n");
   fprintf(fp, "end;\n");
   fprintf(fp, "axis([0  n+1 ymin ymax])\n");
   if (iNames != NULL)
   {
      fprintf(fp, "set(gca,'XTick',[1:%d],'XTickLabel',{", nInputs);
      for (ii = 0; ii < nInputs-1; ii++) fprintf(fp,"'%s',", iNames[ii]);
      fprintf(fp,"'%s'})\n", iNames[nInputs-1]);
   }
   else fwritePlotXLabel(fp, "Input Numbers");
   fwritePlotYLabel(fp, "Individual Gradients");
   fwritePlotTitle(fp, "Scatter Plots for the Gradients (colored)");
   fwritePlotAxes(fp);
   fprintf(fp, "disp('from lo to hi : red,green,blue,magenta,cyan')\n");
   fprintf(fp, "hold off\n");
   fprintf(fp, "disp('Press enter to display the modified gradients')\n");
   fprintf(fp, "hold off\n");
   fprintf(fp, "pause\n");
   fprintf(fp, "clf\n");
   // scatter plots for modified gradients (multi-color)
   fprintf(fp, "hold on\n");
   fprintf(fp, "ymax = max(abs(A(:,2)));\n");
   fprintf(fp, "ymin = min(abs(A(:,2)));\n");
   fprintf(fp, "n = %d;\n", nInputs);
   fprintf(fp, "for ii = 1 : n\n");
   fprintf(fp, "  inds = find(A(:,1)==ii);\n");
   fprintf(fp, "  leng = length(inds);\n");
   fprintf(fp, "  if (leng > 0)\n");
   fprintf(fp, "    [AA,II] = sort(A(inds,3));\n");
   fprintf(fp, "    xx = A(inds(II(1)),1);\n");
   fprintf(fp, "    yy = abs(A(inds(II(1)),2));\n");
   fprintf(fp, "    plot(xx,yy,'ro','MarkerSize',11,'MarkerFaceColor',");
   fprintf(fp, "'r','MarkerEdgeColor','k');\n");
   fprintf(fp, "    for jj = 2 : leng\n");
   fprintf(fp, "      x1 = A(inds(II(jj)),3);\n");
   fprintf(fp, "      x2 = A(inds(II(jj-1)),3);\n");
   fprintf(fp, "      if x1 ~= x2\n");
   fprintf(fp, "        break;\n");
   fprintf(fp, "      else\n");
   fprintf(fp, "        xx = A(inds(II(jj)),1);\n");
   fprintf(fp, "        yy = abs(A(inds(II(jj)),2));\n");
   fprintf(fp, "        plot(xx,yy,'ro','MarkerSize',11,'MarkerFaceColor',");
   fprintf(fp, "'r','MarkerEdgeColor','k');\n");
   fprintf(fp, "        if (jj == leng)\n");
   fprintf(fp, "          jj = jj + 1;\n");
   fprintf(fp, "        end;\n");
   fprintf(fp, "      end;\n");
   fprintf(fp, "    end;\n");
   fprintf(fp, "    if (jj <= leng)\n");
   fprintf(fp, "      xx = A(inds(II(jj)),1);\n");
   fprintf(fp, "      yy = abs(A(inds(II(jj)),2));\n");
   fprintf(fp, "      plot(xx+0.1,yy,'go','MarkerSize',11,'MarkerFaceColor',");
   fprintf(fp, "'g','MarkerEdgeColor','k');\n");
   fprintf(fp, "    end;\n");
   fprintf(fp, "    for kk = jj+1 : leng\n");
   fprintf(fp, "      x1 = A(inds(II(kk)),3);\n");
   fprintf(fp, "      x2 = A(inds(II(kk-1)),3);\n");
   fprintf(fp, "      if x1 ~= x2\n");
   fprintf(fp, "        break;\n");
   fprintf(fp, "      else\n");
   fprintf(fp, "        xx = A(inds(II(kk)),1);\n");
   fprintf(fp, "        yy = abs(A(inds(II(kk)),2));\n");
   fprintf(fp, "        plot(xx+0.1,yy,'go','MarkerSize',11,'MarkerFaceColor',");
   fprintf(fp, "'g','MarkerEdgeColor','k');\n");
   fprintf(fp, "        if (kk == leng)\n");
   fprintf(fp, "          kk = kk + 1;\n");
   fprintf(fp, "        end;\n");
   fprintf(fp, "      end;\n");
   fprintf(fp, "    end;\n");
   fprintf(fp, "    if (kk <= leng)\n");
   fprintf(fp, "      xx = A(inds(II(kk)),1);\n");
   fprintf(fp, "      yy = abs(A(inds(II(kk)),2));\n");
   fprintf(fp, "      plot(xx+0.2,yy,'bo','MarkerSize',11,'MarkerFaceColor',");
   fprintf(fp, "'b','MarkerEdgeColor','k');\n");
   fprintf(fp, "    end;\n");
   fprintf(fp, "    for jj = kk+1 : leng\n");
   fprintf(fp, "      x1 = A(inds(II(jj)),3);\n");
   fprintf(fp, "      x2 = A(inds(II(jj-1)),3);\n");
   fprintf(fp, "      if x1 ~= x2\n");
   fprintf(fp, "        break;\n");
   fprintf(fp, "      else\n");
   fprintf(fp, "        xx = A(inds(II(jj)),1);\n");
   fprintf(fp, "        yy = abs(A(inds(II(jj)),2));\n");
   fprintf(fp, "        plot(xx+0.2,yy,'bo','MarkerSize',11,'MarkerFaceColor',");
   fprintf(fp, "'b','MarkerEdgeColor','k');\n");
   fprintf(fp, "        if (jj == leng)\n");
   fprintf(fp, "          jj = jj + 1;\n");
   fprintf(fp, "        end;\n");
   fprintf(fp, "      end;\n");
   fprintf(fp, "    end;\n");
   fprintf(fp, "    if (jj <= leng)\n");
   fprintf(fp, "      xx = A(inds(II(jj)),1);\n");
   fprintf(fp, "      yy = abs(A(inds(II(jj)),2));\n");
   fprintf(fp, "      plot(xx+0.3,yy,'mo','MarkerSize',11,'MarkerFaceColor',");
   fprintf(fp, "'m','MarkerEdgeColor','k');\n");
   fprintf(fp, "    end;\n");
   fprintf(fp, "    for kk = jj+1 : leng\n");
   fprintf(fp, "      x1 = A(inds(II(kk)),3);\n");
   fprintf(fp, "      x2 = A(inds(II(kk-1)),3);\n");
   fprintf(fp, "      if x1 ~= x2\n");
   fprintf(fp, "        break;\n");
   fprintf(fp, "      else\n");
   fprintf(fp, "        xx = A(inds(II(kk)),1);\n");
   fprintf(fp, "        yy = abs(A(inds(II(kk)),2));\n");
   fprintf(fp, "        plot(xx+0.3,yy,'mo','MarkerSize',11,'MarkerFaceColor',");
   fprintf(fp, "'m','MarkerEdgeColor','k');\n");
   fprintf(fp, "        if (kk == leng)\n");
   fprintf(fp, "          kk = kk + 1;\n");
   fprintf(fp, "        end;\n");
   fprintf(fp, "      end;\n");
   fprintf(fp, "    end;\n");
   fprintf(fp, "    if (kk <= leng)\n");
   fprintf(fp, "      for jj = kk : leng\n");
   fprintf(fp, "        xx = A(inds(II(jj)),1);\n");
   fprintf(fp, "        yy = abs(A(inds(II(jj)),2));\n");
   fprintf(fp, "        plot(xx+0.4,yy,'co','MarkerSize',11,'MarkerFaceColor',");
   fprintf(fp, "'c','MarkerEdgeColor','k');\n");
   fprintf(fp, "      end;\n");
   fprintf(fp, "    end;\n");
   fprintf(fp, "  end;\n");
   fprintf(fp, "end;\n");
   fprintf(fp, "for ii = 1 : n\n");
   fprintf(fp, "   inds = find(A(:,1)==ii);\n");
   fprintf(fp, "   if length(inds) > 0\n");
   fprintf(fp, "      asum = sum(abs(A(inds,2)))/length(inds);\n");
   fprintf(fp, "      plot(ii,asum,'k*','MarkerSize',13);\n");
   fprintf(fp, "   end;\n");
   fprintf(fp, "end;\n");
   fprintf(fp, "axis([0 n+1 ymin ymax])\n");
   fwritePlotTitle(fp, "Scatter Plots for the Modified Gradients (colored)");
   if (iNames != NULL)
   {
      fprintf(fp, "set(gca,'XTick',[1:%d],'XTickLabel',{", nInputs);
      for (ii = 0; ii < nInputs-1; ii++) fprintf(fp,"'%s',", iNames[ii]);
      fprintf(fp,"'%s'})\n", iNames[nInputs-1]);
   }
   else fwritePlotXLabel(fp, "Input Numbers");
   fwritePlotAxes(fp);
   fwritePlotYLabel(fp, "Individual Modified Gradients");
   fprintf(fp, "disp('from lo to hi : red,green,blue,magenta,cyan')\n");
   fprintf(fp, "hold off\n");
   fclose(fp);
   printf("MOAT scatter plot matlab file = %s.\n", scatterFile);
   printAsterisks(0);
   return 0;
}

// ************************************************************************
// create bootstrap mean plot matlab file
// ------------------------------------------------------------------------
int MOATAnalyzer::createBootstrapFile(int nSamples, int nInputs, double *Y, 
                                      double *X, int *indices, char **iNames)
{
   int    nReps, *validCnts, ii, is, jj, ib, input, count, index;
   double **YGs, *means, *stds, *bArray;
   char   winput[500], bootstrapFile[500], pString[500];
   FILE   *fp;

   if (psPlotTool_ == 1)
   {
      printf("MOAT: bootstrap plot not availabe with scilab.\n");
      return 0;
   }
   printAsterisks(0);
   sprintf(pString, "Create bootstrap modified mean plot ? (y or n) ");
   getString(pString, winput);
   if (winput[0] != 'y')
   {
      return 0;
   }
   sprintf(pString, "Enter matlab/scilab bootstrap file name (no extension): ");
   getString(pString, bootstrapFile);
   count = strlen(bootstrapFile);
   if (count > 500)
   {
      printf("ERROR: file name too long.\n");
      exit(1);
   }
   bootstrapFile[count-1] = '.';
   if (psPlotTool_ == 1)
   {
      bootstrapFile[count] = 's';
      bootstrapFile[count+1] = 'c';
      bootstrapFile[count+2] = 'i';
      bootstrapFile[count+3] = '\0';
   }
   else
   {
      bootstrapFile[count] = 'm';
      bootstrapFile[count+1] = '\0';
   }

   fp = fopen(bootstrapFile, "w");
   if (fp == NULL)
   {
      printf("MOATAnalyzer: cannot write to plot file %s\n",bootstrapFile);
      return 0;
   }

   nReps = nSamples / (nInputs + 1);
   validCnts = new int[nInputs];
   for (ii = 0; ii < nInputs; ii++) validCnts[ii] = 0;
   for (is = 0; is < nSamples; is++)
   {
      if (Y[is] != PSUADE_UNDEFINED && indices[is] >= 0)
         validCnts[indices[is]]++;
   }
   for (ii = 0; ii < nInputs; ii++)
      if (validCnts[ii] > nReps) nReps = validCnts[ii];
   YGs = new double*[nInputs]; 
   for (ii = 0; ii < nInputs; ii++)
   {
      if (validCnts[ii] > 0) YGs[ii] = new double[validCnts[ii]];
      else                   YGs[ii] = NULL;
   }
   for (ii = 0; ii < nInputs; ii++) validCnts[ii] = 0;
   for (is = 0; is < nSamples; is++)
   {
      if (Y[is] != PSUADE_UNDEFINED && indices[is] >= 0)
      {
         input = indices[is];
         count = validCnts[input]++;
         YGs[input][count] = PABS(Y[is]);
      }
   }

   bArray = new double[1000];
   means = new double[nInputs];
   stds = new double[nInputs];
   for (ii = 0; ii < nInputs; ii++)
   {
      if (validCnts[ii] > 5)
      {
         for (ib = 0; ib < 1000; ib++)
         {
            bArray[ib] = 0.0;
            for (jj = 0; jj < validCnts[ii]; jj++)
            {
               index = PSUADE_rand() % validCnts[ii];
               bArray[ib] += YGs[ii][index];
            }
            bArray[ib] /= validCnts[ii];
         }
         means[ii] = 0.0;
         for (ib = 0; ib < 1000; ib++) means[ii] += bArray[ib];
         means[ii] /= 1000.0;
         stds[ii] = 0.0;
         for (ib = 0; ib < 1000; ib++)
            stds[ii] += pow(bArray[ib] - means[ii], 2.0);
         stds[ii] /= (1000.0 - 1.0);
         stds[ii] = sqrt(stds[ii]);
      }
      else
      {
         printf("MOATAnalyzer WARNING: input %d needs > 5 replications\n",ii+1);
         printf("                      do perform bootstrapping.\n");
         stds[ii] = 0.0;
         means[ii] = 0.0;
         for (jj = 0; jj < validCnts[ii]; jj++) means[ii] += YGs[ii][jj];
         if (validCnts[ii] > 0) means[ii] /= (double) validCnts[ii];
      }
   }

   fprintf(fp, "%% This file contains modified means of gradients\n");
   fprintf(fp, "%% and also their spreads based on bootstraping.\n");
   fprintf(fp, "Means = [\n");
   for (ii = 0; ii < nInputs; ii++) fprintf(fp,"%24.16e\n", means[ii]);
   fprintf(fp, "];\n");
   fprintf(fp, "Stds = [\n");
   for (ii = 0; ii < nInputs; ii++) fprintf(fp,"%24.16e\n", stds[ii]);
   fprintf(fp, "];\n");
   fprintf(fp, "clf\n");
   fprintf(fp, "n = %d;\n", nInputs);
   fprintf(fp, "%%n = 20;\n");
   fprintf(fp, "%%[Means, I2] = sort(Means,'descend');\n");
   fprintf(fp, "%%Stds = Stds(I2);\n");
   fprintf(fp, "%%Means = Means(1:n);\n");
   fprintf(fp, "%%Stds = Stds(1:n);\n");
   fprintf(fp, "%%Str = {");
   for (ii = 0; ii < nInputs-1; ii++) fprintf(fp,"'X%d',",ii+1);
   fprintf(fp,"'X%d'};\n",nInputs);
   fprintf(fp, "%%Str = Str(I2);\n");
   fprintf(fp, "ymin = min(Means-2*Stds);\n");
   fprintf(fp, "ymax = max(Means+2*Stds);\n");
   fprintf(fp, "h2 = 0.05 * (ymax - ymin);\n");
   fprintf(fp, "for ii = 1:n\n");
   fprintf(fp, "   h = plot(ii,Means(ii),'r*','MarkerSize',13);\n");
   fprintf(fp, "   if (ii == 1)\n");
   fprintf(fp, "      hold on\n");
   fprintf(fp, "   end;\n");
   fprintf(fp, "   XX = [ii ii];\n");
   fprintf(fp, "   YY = [Means(ii)-Stds(ii) Means(ii)+Stds(ii)];\n");
   fprintf(fp, "%%   text(ii-0.2,Means(ii)+Stds(ii)+h2,Str(ii),");
   fprintf(fp, "'FontSize',12,'FontWeight','bold')\n");
   fprintf(fp, "   plot(XX,YY,'-ko','LineWidth',2.0,'MarkerEdgeColor',");
   fprintf(fp, "'k','MarkerFaceColor','g','MarkerSize',11)\n");
   fprintf(fp, "end;\n");
   fprintf(fp, "axis([0 n+1 ymin ymax])\n");
   fwritePlotAxes(fp);
   fwritePlotTitle(fp,"Modified Means Plot (bootstrap)");
   fwritePlotYLabel(fp, "Modified Means (of gradients)");
   if (iNames != NULL)
   {
      fprintf(fp, "set(gca,'XTick',[1:%d],'XTickLabel',{", nInputs);
      for (ii = 0; ii < nInputs-1; ii++) fprintf(fp,"'%s',", iNames[ii]);
      fprintf(fp,"'%s'})\n", iNames[nInputs-1]);
   }
   else fwritePlotXLabel(fp, "Input Numbers");
   fprintf(fp, "hold off\n");
   fclose(fp);
   printf("MOAT bootstrap plot matlab file = %s.\n", bootstrapFile);
   printAsterisks(0);

   for (ii = 0; ii < nInputs; ii++) if (YGs[ii] != NULL) delete [] YGs[ii];
   delete [] YGs;
   delete [] validCnts;
   delete [] means;
   delete [] stds;
   delete [] bArray;
   return 0;
}

