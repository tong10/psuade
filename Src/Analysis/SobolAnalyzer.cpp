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
// Functions for the class SobolAnalyzer (TS, S, interactions)  
// AUTHOR : CHARLES TONG
// DATE   : 2004
// ************************************************************************
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "SobolAnalyzer.h"
#include "Psuade.h"
#include "sysdef.h"
#include "PsuadeUtil.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
SobolAnalyzer::SobolAnalyzer() : Analyzer()
{
   setName("SOBOL");
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
SobolAnalyzer::~SobolAnalyzer()
{
}

// ************************************************************************
// perform analysis
// ------------------------------------------------------------------------
double SobolAnalyzer::analyze(aData &adata)
{
   int     nInputs, nOutputs, nSamples, outputID, count, errFlag;
   int     repID, iD, ii, nReps, ss, errCount, *sIn;
   double  *Y, xtemp1, xtemp2, *means, *modifiedMeans, *stds;
   double  *S, *ST, *PE, tau, sMean, sVar, sMean2, sVar2;
   double  *xIn, *yIn, *xLower, *xUpper, dtemp;

   nInputs  = adata.nInputs_;
   nOutputs = adata.nOutputs_;
   nSamples = adata.nSamples_;
   xLower   = adata.iLowerB_;
   xUpper   = adata.iUpperB_;
   xIn      = adata.sampleInputs_;
   yIn      = adata.sampleOutputs_;
   sIn      = adata.sampleStates_;
   outputID = adata.outputID_;
   if (adata.inputPDFs_ != NULL)
   {
      count = 0;
      for (ii = 0; ii < nInputs; ii++) count += adata.inputPDFs_[ii];
      if (count > 0)
      {
         printf("SobolAnalyzer INFO: some inputs have non-uniform PDFs.\n");
         printf("     However, they are not relevant in this analysis\n");
         printf("     (since the sample should have been generated with\n");
         printf("     the desired distributions.)\n");
      }
   }

   if (nInputs <= 0 || nSamples <= 0 || nOutputs <= 0 || 
       outputID < 0 || outputID >= nOutputs)
   {
      printf("SobolAnalyzer ERROR: invalid arguments.\n");
      printf("    nInputs  = %d\n", nInputs);
      printf("    nOutputs = %d\n", nOutputs);
      printf("    nSamples = %d\n", nSamples);
      printf("    outputID = %d\n", outputID+1);
      return PSUADE_UNDEFINED;
   } 

   count = 0;
   for (ss = 0; ss < nSamples; ss++)
   {
      errFlag = 0;
      if (sIn[ss] != 1) errFlag = 1;
      for (ii = 0; ii < nOutputs; ii++)
         if (yIn[ss*nOutputs+ii] < 0.9*PSUADE_UNDEFINED) errFlag = 1;
      if (errFlag == 0) count++;
   }
   printf("SobolAnalyzer INFO: there are a total of %d sample points.\n",
          nSamples);
   printf("SobolAnalyzer INFO: there are %d valid sample points.\n",
          count);

   nReps = nSamples / (nInputs + 2);
   if ((nReps * (nInputs+2)) == nSamples)
   {
      for (ss = 0; ss < nSamples; ss+=(nInputs+2))
      {
         errCount = 0;
         for (iD = 1; iD <= nInputs; iD++)
         {
            for (ii = 0; ii < nInputs; ii++)
            {
               if (ii == (iD-1))
                    xtemp1 = xIn[(ss+nInputs+1)*nInputs+ii]; 
               else xtemp1 = xIn[ss*nInputs+ii]; 
               xtemp2 = xIn[(ss+iD)*nInputs+ii]; 
               if (xtemp1 != xtemp2) errCount++;
            }
         }
         if (errCount > 0)
         {
            printf("SobolAnalyzer ERROR: invalid sample (%d,%d).\n",
                   ss, errCount);
            printf("SobolAnalyzer requires Sobol samples.\n");
            return PSUADE_UNDEFINED;
         }
      }
   }
   else
   {
      printf("SobolAnalyzer ERROR: invalid sample size.\n");
      printf("SobolAnalyzer requires Sobol samples.\n");
      return PSUADE_UNDEFINED;
   }
   
   Y = new double[nSamples];
   for (ss = 0; ss < nSamples; ss++) Y[ss] = yIn[nOutputs*ss+outputID];
   means = new double[nInputs];
   modifiedMeans = new double[nInputs];
   stds = new double[nInputs];
   for (ii = 0; ii < nInputs; ii++)
      means[ii] = modifiedMeans[ii] = stds[ii] = 0.0;

   MOATAnalyze(nInputs,nSamples,xIn,Y,xLower,xUpper,means,
               modifiedMeans,stds);

   printf("Sobol-OAT Analysis : \n");
   for (ii = 0; ii < nInputs; ii++)
      printf("Input %3d (mod. mean & std) = %12.4e %12.4e\n",
             ii+1, modifiedMeans[ii], stds[ii]);


   for (ss = 0; ss < nSamples; ss++) Y[ss] = yIn[nOutputs*ss+outputID];

   sMean = 0.0;
   count = 0;
   for (repID = 0;  repID < nReps; repID++)
   {
      if (Y[repID*(nInputs+2)+nInputs+1] < 0.9*PSUADE_UNDEFINED) 
      {
         sMean += Y[repID*(nInputs+2)+nInputs+1]; 
         count++;
      }
   }
   if (count <= 1)
   {
      printf("SobolAnalyzer ERROR: too few valid sample points.\n");
      exit(1);
   }
   sMean /= ((double) (count));
   sVar = 0.0;
   for (repID = 0;  repID < nReps; repID++)
      if (Y[repID*(nInputs+2)+nInputs+1] < 0.9*PSUADE_UNDEFINED) 
         sVar += ((Y[repID*(nInputs+2)+nInputs+1] - sMean) * 
                  (Y[repID*(nInputs+2)+nInputs+1] - sMean)); 
   sVar = sVar / (double) (count-1.0);
   if (sVar == 0)
   {
      printf("SobolAnalyzer ERROR: sample variance = 0.0.\n");
      exit(1);
   }

   S  = new double[nInputs];
   ST = new double[nInputs];
   PE = new double[nInputs];
   for (ii = 0; ii < nInputs; ii++)
   {
      tau = 0.0;
      count = 0;
      for (repID = 0;  repID < nReps; repID++)
      {
         if ((Y[repID*(nInputs+2)+ii+1] < 0.9*PSUADE_UNDEFINED) &&
             (Y[repID*(nInputs+2)] < 0.9*PSUADE_UNDEFINED))
         {
            tau += ((Y[repID*(nInputs+2)] - sMean) *
                    (Y[repID*(nInputs+2)+ii+1] - sMean)); 
            count++;
         }
      }
      if (count <= 0)
      {
         printf("SobolAnalyzer ERROR: too few valid sample points for TSI.\n");
         exit(1);
      }
      tau /= ((double) count);
      ST[ii] = 1.0 - tau / sVar; 

      tau = 0.0;
      count = 0;
      for (repID = 0;  repID < nReps; repID++)
      {
         if ((Y[repID*(nInputs+2)+ii+1] < 0.9*PSUADE_UNDEFINED) &&
             (Y[repID*(nInputs+2)+nInputs+1] < 0.9*PSUADE_UNDEFINED))
         {
            tau += ((Y[repID*(nInputs+2)+nInputs+1] - sMean) *
                    (Y[repID*(nInputs+2)+ii+1] - sMean)); 
            count++;
         }
      }
      if (count <= 0)
      {
         printf("SobolAnalyzer ERROR: too few valid sample points for VCE.\n");
         exit(1);
      }
      tau /= ((double) count);
      S[ii] = tau / sVar; 

      sMean2 = 0.0;
      count  = 0;
      for (repID = 0;  repID < nReps; repID++)
      {
         if ((Y[repID*(nInputs+2)+ii+1] < 0.9*PSUADE_UNDEFINED) &&
             (Y[repID*(nInputs+2)+nInputs+1] < 0.9*PSUADE_UNDEFINED))
         {
            sMean2 += (Y[repID*(nInputs+2)+nInputs+1]*Y[repID*(nInputs+2)+ii+1]);
            count++;
         }
      }
      sMean2 = sMean2 / (double) count;
      sVar2 = 0.0;
      for (repID = 0;  repID < nReps; repID++)
      {
         if ((Y[repID*(nInputs+2)+ii+1] < 0.9*PSUADE_UNDEFINED) &&
             (Y[repID*(nInputs+2)+nInputs+1] < 0.9*PSUADE_UNDEFINED))
         {
            dtemp = Y[repID*(nInputs+2)+ii+1] * Y[repID*(nInputs+2)+nInputs+1];
            sVar2 += pow(dtemp - sMean2, 2.0);
         }
      }
      sVar2 = sVar2 / count;
      PE[ii] = 0.6745 * sqrt(sVar2) / sqrt(1.0 * count);
   }

   printf("Sobol Analysis (ST: total sensitivity, PE: probable error):\n");
   for (ii = 0; ii < nInputs; ii++)
      printf("Input %3d (S, ST, PE) = %12.4e %12.4e %12.4e\n",
             ii+1, S[ii], ST[ii], PE[ii]);

   delete [] Y;
   delete [] means;
   delete [] modifiedMeans;
   delete [] stds;
   delete [] ST;
   delete [] S;
   delete [] PE;
   return 0.0;
}

// ************************************************************************
// perform analysis similar to MOAT analysis
// ------------------------------------------------------------------------
int SobolAnalyzer::MOATAnalyze(int nInputs, int nSamples, double *xIn,
                         double *yIn, double *xLower, double *xUpper,
                         double *means, double *modifiedMeans, double *stds)
{
   int    ss, ii, *counts;
   double xtemp1, xtemp2, ytemp1, ytemp2, scale;
   FILE   *fp;

   for (ss = 0; ss < nSamples; ss+=(nInputs+2))
   {
      for (ii = 1; ii <= nInputs; ii++)
      {
         ytemp1 = yIn[ss+ii]; 
         ytemp2 = yIn[ss]; 
         xtemp1 = xIn[(ss+ii)*nInputs+ii-1]; 
         xtemp2 = xIn[ss*nInputs+ii-1]; 
         scale  = xUpper[ii-1] - xLower[ii-1];
         if (xtemp1 != xtemp2)
         {
            yIn[ss+ii] = (ytemp2-ytemp1)/(xtemp2-xtemp1)*scale;
         }
         else
         {
            printf("SobolAnalyzer ERROR: divide by 0.\n");
            printf("     Check sample (Is this Sobol?) \n");
            exit(1);
         }
      }
   }

   counts = new int[nInputs];
   for (ii = 0; ii < nInputs; ii++) counts[ii] = 0;
   for (ss = 0; ss < nSamples; ss+=(nInputs+2))
   {
      for (ii = 1; ii <= nInputs; ii++)
      {
         if (yIn[ss+ii] < 0.9*PSUADE_UNDEFINED)
         {
            means[ii-1] += yIn[ss+ii];
            modifiedMeans[ii-1] += PABS(yIn[ss+ii]);
            counts[ii-1]++;
         }
      }
   }
   for (ii = 0; ii < nInputs; ii++)
   {
      if (counts[ii] > 0)
      {
         means[ii] /= (double) (counts[ii]);
         modifiedMeans[ii] /= (double) (counts[ii]);
      }
   }
   for (ss = 0; ss < nSamples; ss+=(nInputs+2))
   {
      for (ii = 1; ii <= nInputs; ii++)
      {
         if (yIn[ss+ii] < 0.9*PSUADE_UNDEFINED)
            stds[ii-1] += (yIn[ss+ii] - means[ii-1]) *
                          (yIn[ss+ii] - means[ii-1]);
      }
   }
   for (ii = 0; ii < nInputs; ii++)
      if (counts[ii] > 0)
         stds[ii] /= (double) (counts[ii]);
   for (ii = 0; ii < nInputs; ii++) stds[ii] = sqrt(stds[ii]);
   delete [] counts;

   if (psPlotTool_ == 1) fp = fopen("scilabsmp.sci", "w");
   else                  fp = fopen("matlabsmp.m", "w");
   if (fp == NULL)
   {
      printf("ERROR: cannot open file to write statistics.\n");
      return 0;
   }
   fprintf(fp, "Y = [\n");
   for (ii = 0; ii < nInputs; ii++) fprintf(fp, "%24.16e\n", stds[ii]);
   fprintf(fp, "];\n");
   fprintf(fp, "X = [\n");
   for (ii = 0; ii < nInputs; ii++) 
   fprintf(fp, "%24.16e\n",modifiedMeans[ii]);
   fprintf(fp, "];\n");
   fprintf(fp, "plot(X,Y,'*')\n");
   fwritePlotAxes(fp);
   fwritePlotXLabel(fp, "Modified Means");
   fwritePlotYLabel(fp, "Std Devs");
   fprintf(fp, "title('Std Devs vs Modified mean Plot')\n");
   fclose(fp);
   if (psPlotTool_ == 1) printf("scilab file ready: scilabsmp.sci");
   else                  printf("matlab file ready: matlabsmp.m\n");
   return 0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
SobolAnalyzer& SobolAnalyzer::operator=(const SobolAnalyzer &)
{
   printf("SobolAnalyzer operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

