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
// functions for the class LSAnalyzer (Local sensitivity analysis)  
// AUTHOR : CHARLES TONG
// DATE   : 2012
// ************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "LSAnalyzer.h"
#include "Util/sysdef.h"
#include "Util/PsuadeUtil.h"
#include "Main/Psuade.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
LSAnalyzer::LSAnalyzer() : Analyzer()
{
   setName("LSA");
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
LSAnalyzer::~LSAnalyzer()
{
}

// ************************************************************************
// perform analysis
// ------------------------------------------------------------------------
double LSAnalyzer::analyze(aData &adata)
{
   int    nInputs, nOutputs, nSamples, outputID;
   int    ii, ss, whichOutput, ncount, index;
   double *X, *Y, *mEffect, dtemp;
   char   winput1[500], meFileName[500], pString[500];
   FILE   *fp=NULL;

   nInputs  = adata.nInputs_;
   nOutputs = adata.nOutputs_;
   nSamples = adata.nSamples_;
   outputID = adata.outputID_;
   X        = adata.sampleInputs_;
   Y        = adata.sampleOutputs_;
   if (adata.inputPDFs_ != NULL)
   {
      printf("Local SA INFO: some inputs have non-uniform PDFs, but\n");
      printf("               they will not be relevant in this analysis.\n");
   }
   whichOutput = outputID;
   if (whichOutput >= nOutputs || whichOutput < 0) whichOutput = 0;

   if (nInputs <= 0 || nOutputs <= 0 || nSamples <= 0)
   {
      printf("LSAnalyzer ERROR: invalid arguments.\n");
      printf("   nInputs  = %d\n", nInputs);
      printf("   nOutputs = %d\n", nOutputs); 
      printf("   nSamples = %d\n", nSamples);
      return PSUADE_UNDEFINED;
   } 
   if (nSamples != nInputs + 1)
   {
      printf("LSAnalyzer ERROR: nSamples should be equal to nInputs+1.\n");
      printf("   nSamples = %d\n", nSamples);
      return PSUADE_UNDEFINED;
   } 

   if (psPlotTool_ == 0)
   {
      if (psAnaExpertMode_ == 1)
      {
         sprintf(pString,"Create sensitivity plot ? (y or n) ");
         getString(pString, winput1);
         if (winput1[0] == 'y')
         {
            sprintf(pString,"Enter sensitivity plot matlab file name : ");
            getString(pString, meFileName);
            if (strlen(meFileName) < 500)
               meFileName[strlen(meFileName)-1] = '\0';
            else
            {
               printf("ERROR: file name too long.\n");
               exit(1);
            }
         }
         fp = fopen(meFileName, "w");
         if (fp != NULL)
         {
            printf("LSAnalyzer: sensitivity effect matlab file = %s\n",
                   meFileName);
         }
      }
   }

   printf("\n");
   printAsterisks(0);
   printAsterisks(0);
   printf("* Local Sensitivity Analysis\n");
   printDashes(0);
   printf("* total number of samples = %d\n",nSamples);
   printf("* number of Inputs        = %d\n",nInputs);
   printf("* Output number           = %d\n", whichOutput+1);
   printDashes(0);
   if (fp != NULL)
   {
      fprintf(fp,"%% ********************************************** \n");
      fprintf(fp,"%% ********************************************** \n");
      fprintf(fp,"%% * Local Sensitivity Analysis                ** \n");
      fprintf(fp,"%% *-------------------------------------------** \n");
      fprintf(fp,"%% Output %d\n", whichOutput+1);
      fprintf(fp,"%% *-------------------------------------------** \n");
   }

   mEffect = new double[nInputs];
   for (ii = 0; ii < nInputs; ii++) mEffect[ii] = PSUADE_UNDEFINED;

   for (ss = 1; ss < nSamples; ss++)
   {
      ncount = 0;
      for (ii = 0; ii < nInputs; ii++)
      {
         if (PABS(X[nInputs*ss+ii]-X[ii]) > 1.0e-15) 
         {
            ncount++;
            index = ii;
         }
      }
      if (ncount == 1)
      {
         mEffect[index] = Y[nOutputs*ss+whichOutput] - Y[whichOutput];
         dtemp = X[nInputs*ss+ii] - X[ii];
         dtemp /= (PABS(X[nInputs*ss+ii]-X[ii]));
         mEffect[index] *= dtemp; 
      }
   } 
   for (ii = 0; ii < nInputs; ii++)
   {
      if (mEffect[ii] == PSUADE_UNDEFINED)
      {
         printf("LSAnalyzer ERROR: sample not suitable for LSA.\n");
         delete [] mEffect;
         return -1;
      }
   }
   for (ii = 0; ii < nInputs; ii++)
      printf("* Input %3d (importance measure = %12.4e)\n", ii+1, 
             mEffect[ii]);

   if (fp != NULL)
   {
      fprintf(fp, "Y = [\n");
      for (ii = 0; ii < nInputs; ii++)
         fprintf(fp, "%24.16e \n", mEffect[ii]);
      fprintf(fp, "]; \n");
      fwritePlotCLF(fp);
      fprintf(fp, "plot(Y,'*','MarkerSize',12); \n");
      fwritePlotAxes(fp);
      sprintf(pString, "Sensitivity plot for Output %d\n",whichOutput+1);
      fwritePlotTitle(fp, pString);
      sprintf(pString, "Input Numbers\n");
      fwritePlotXLabel(fp, pString);
      sprintf(pString, "Sensitivity Measure\n");
      fwritePlotYLabel(fp, pString);
      fclose(fp);
   }
   printAsterisks(0);

   delete [] mEffect;
   return 0.0;
}

