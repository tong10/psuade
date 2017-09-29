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
// functions for the class FFAnalyzer (Fractional Factorial analysis)  
// AUTHOR : CHARLES TONG
// DATE   : 2005
// ************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "FFAnalyzer.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "Psuade.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
FFAnalyzer::FFAnalyzer() : Analyzer()
{
   setName("FF");
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
FFAnalyzer::~FFAnalyzer()
{
}

// ************************************************************************
// perform analysis
// ------------------------------------------------------------------------
double FFAnalyzer::analyze(aData &adata)
{
   int    nInputs, nOutputs, nSamples, outputID, checkSample, printLevel;
   int    ii, ii2, rr, ss, whichOutput, ncount, nplots, *iArray, nReps, offset;
   double *X, *Y, *txArray, *twArray, *tyArray, **mEffects, accum, ***iEffects;
   double ***lolos, ***lohis, ***hilos, ***hihis, stdev, mean;
   char   pString[500], winput[500];
   FILE   *fp=NULL;

   printLevel = adata.printLevel_;
   nInputs    = adata.nInputs_;
   nOutputs   = adata.nOutputs_;
   nSamples   = adata.nSamples_;
   outputID   = adata.outputID_;
   X          = adata.sampleInputs_;
   Y          = adata.sampleOutputs_;
   if (adata.inputPDFs_ != NULL)
   {
      ncount = 0;
      for (ii = 0; ii < nInputs; ii++) ncount += adata.inputPDFs_[ii];
      if (ncount > 0)
      {
         printf("FFAnalysis INFO: some inputs have non-uniform PDFs, but\n");
         printf("           they are not relevant in this analysis.\n");
      }
   }
   whichOutput = outputID;
   if (whichOutput >= nOutputs || whichOutput < 0) whichOutput = 0;

   if (nInputs <= 0 || nOutputs <= 0 || nSamples <= 0)
   {
      printf("FFAnalysis ERROR: invalid arguments.\n");
      printf("   nInputs  = %d\n", nInputs);
      printf("   nOutputs = %d\n", nOutputs); 
      printf("   nSamples = %d\n", nSamples);
      return PSUADE_UNDEFINED;
   } 
   if (nSamples / 2 * 2 != nSamples)
   {
      printf("FFAnalysis ERROR: nSamples has to be even.\n");
      printf("   nSamples = %d\n", nSamples);
      return PSUADE_UNDEFINED;
   } 

   txArray = new double[nSamples];
   tyArray = new double[nSamples];
   twArray = new double[nSamples];

   printf("\n");
   printAsterisks(0);
   printAsterisks(0);
   printf("* Fractional Factorial Main Effect Analysis\n");
   printf("* This analysis works for FF4, FF5 and PBD designs\n");
   printDashes(0);
   printf("* total number of samples = %10d \n",nSamples);
   printf("* number of Inputs        = %10d \n",nInputs);
   printf("* Output number           = %d\n", whichOutput+1);
   printDashes(0);
   fp = NULL;
   if (psAnaExpertMode_ == 1 || psAnalysisInteractive_ == 1)
   {
      sprintf(pString,"Create main effect plot ? (y or n) ");
      getString(pString, winput);
      if (winput[0] == 'y')
      {
         if (psPlotTool_ == 0) fp = fopen("matlabmeff.m", "w");
         else                  fp = fopen("scilabmeff.sci", "w");
      }
   }
   if (fp != NULL)
   {
      if (psPlotTool_ == 1)
      {
         fprintf(fp,"// ********************************************** \n");
         fprintf(fp,"// ********************************************** \n");
         fprintf(fp,"// * Fractional Factorial Main Effect Analysis ** \n");
         fprintf(fp,"// * (2 level fractional factorial only)       ** \n");
         fprintf(fp,"// *-------------------------------------------** \n");
         fprintf(fp,"// Output %d\n", whichOutput+1);
         fprintf(fp,"// *-------------------------------------------** \n");
      }
      else
      {
         fprintf(fp,"%% ********************************************** \n");
         fprintf(fp,"%% ********************************************** \n");
         fprintf(fp,"%% * Fractional Factorial Main Effect Analysis ** \n");
         fprintf(fp,"%% * (2 level fractional factorial only)       ** \n");
         fprintf(fp,"%% *-------------------------------------------** \n");
         fprintf(fp,"%% Output %d\n", whichOutput+1);
         fprintf(fp,"%% *-------------------------------------------** \n");
      }
   }

   nReps = 1;
   for (ii = 0; ii < nInputs; ii++)
   {
      for (ss = 0; ss < nSamples; ss++)
      {
         txArray[ss] = X[nInputs*ss+ii];
         tyArray[ss] = Y[nOutputs*ss+whichOutput];
      }
      sortDbleList2(nSamples, txArray, tyArray);
      checkSample = 1;
      for (ss = 1; ss < nSamples/2; ss++)
         if (txArray[ss] != txArray[0]) checkSample = 0;
      for (ss = nSamples/2+1; ss < nSamples; ss++)
         if (txArray[ss] != txArray[nSamples/2]) checkSample = 0;
      if (checkSample == 0)
      {
         printf("FFAnalysis ERROR: sample not fractional factorial.\n");
         printf("If you are using replicated Fractional Factorial\n");
         printf("enter the number of replications.\n");
         sprintf(pString, "Number of replications = (2 - %d) ", nSamples/2);
         nReps = getInt(2, nSamples/2, pString);
         break;
      }
   }
   if (nReps > 1)
   {
      for (ii = 0; ii < nInputs; ii++)
      {
         for (ii2 = 0; ii2 < nReps; ii2++)
         {
            for (ss = 0; ss < nSamples/nReps; ss++)
            {
               txArray[ss] = X[(ii2*nSamples/nReps+ss)*nInputs+ii];
               tyArray[ss] = Y[(ii2*nSamples/nReps+ss)*nOutputs+whichOutput];
            }
            sortDbleList2(nSamples/nReps, txArray, tyArray); 
            checkSample = 1;
            for (ss = 1; ss < nSamples/nReps/2; ss++)
               if (txArray[ss] != txArray[0]) checkSample = 0;
            for (ss = nSamples/nReps/2+1; ss < nSamples/nReps; ss++)
               if (txArray[ss] != txArray[nSamples/nReps-1])
                  checkSample = 0;
            if (checkSample == 0)
            {
               printf("FFAnalysis ERROR: sample not fractional factorial.\n");
               printf("                  nor replicated fractional factorial.\n");
               delete [] twArray;
               delete [] txArray;
               delete [] tyArray;
               if(fp != NULL) fclose(fp);
               return 0.0;
            }
         }
      }
   }

   mEffects = new double*[nInputs];
   for (ii = 0; ii < nInputs; ii++) mEffects[ii] = new double[nReps+1];
   iArray  = new int[nInputs];
   for (ii = 0; ii < nInputs; ii++) iArray[ii] = ii;

   for (rr = 0; rr < nReps; rr++)
   {
      offset = rr * nSamples / nReps;
      for (ii = 0; ii < nInputs; ii++)
      {
         for (ss = 0; ss < nSamples/nReps; ss++)
         {
            txArray[ss] = X[(offset+ss)*nInputs+ii];
            tyArray[ss] = Y[(offset+ss)*nOutputs+whichOutput];
         }
         sortDbleList2(nSamples/nReps, txArray, tyArray); 
         accum = 0.0;
         for (ss = nSamples/nReps/2; ss < nSamples/nReps; ss++)
            accum += tyArray[ss];
         for (ss = 0; ss < nSamples/nReps/2; ss++) accum -= tyArray[ss];
         mEffects[ii][rr+1] = accum * 2 / (double) (nSamples/nReps);
      }
   }
   printf("* Fractional Factorial Main Effect (normalized)\n");
   printf("* Note: std err is the standard error or mean.\n");
   printDashes(0);
   for (ii = 0; ii < nInputs; ii++)
   {
      mean = 0.0;
      for (rr = 0; rr < nReps; rr++) mean += mEffects[ii][rr+1];
      mean /= (double) nReps;
      mEffects[ii][0] = mean;
      stdev = 0.0;
      if (nReps > 1)
      {
         for (rr = 0; rr < nReps; rr++)
            stdev += pow(mEffects[ii][rr+1]-mean, 2.0);
         stdev = sqrt(stdev / (double) (nReps - 1.0));
         printf("* Input %3d  =  %12.4e (std err = %12.4e)\n", ii+1,mean,
                  stdev/sqrt(1.0*nReps));
      }  
      else printf("* Input %3d  =  %12.4e\n", ii+1, mean);
   }
   printDashes(0);
   printf("* Fractional Factorial Main Effect (ordered)\n");
   printDashes(0);
   for (ii = 0; ii < nInputs; ii++) twArray[ii] = PABS(mEffects[ii][0]);
   sortDbleList2a(nInputs, twArray, iArray);
   for (ii = nInputs-1; ii >= 0; ii--)
      printf("* Rank %3d : Input %3d (measure = %12.4e)\n", nInputs-ii, 
             iArray[ii]+1, twArray[ii]);

   if (fp != NULL)
   {
      fprintf(fp, "Y = [\n");
      for (ii = 0; ii < nInputs; ii++)
         fprintf(fp, "%24.16e \n", mEffects[ii][0]);
      fprintf(fp, "]; \n");
      fwritePlotCLF(fp);
      fprintf(fp, "plot(Y,'*'); \n");
      fwritePlotAxes(fp);
      fwritePlotTitle(fp, "Main Effect Plot for Output");
      fprintf(fp, "disp('Press enter to continue to the next plot')\n");
      fprintf(fp, "pause\n");
   }

   if (nInputs < 2 || adata.samplingMethod_ == PSUADE_SAMP_PBD)
   {
      delete [] txArray;
      delete [] twArray;
      delete [] tyArray;
      delete [] iArray;
      for(int i = 0; i < nInputs; i++) delete [] mEffects[i];
      delete [] mEffects;
      if (fp != NULL) fclose(fp);
      return 0.0;
   }

   iEffects = new double**[nInputs];
   hihis = new double**[nInputs];
   lolos = new double**[nInputs];
   lohis = new double**[nInputs];
   hilos = new double**[nInputs];
   for (ii = 0; ii < nInputs; ii++)
   {
      iEffects[ii] = new double*[nInputs];
      hihis[ii] = new double*[nInputs];
      lolos[ii] = new double*[nInputs];
      lohis[ii] = new double*[nInputs];
      hilos[ii] = new double*[nInputs];
      for (ii2 = 0; ii2 < nInputs; ii2++)
      {
         hihis[ii][ii2] = new double[nReps+1];
         hilos[ii][ii2] = new double[nReps+1];
         lohis[ii][ii2] = new double[nReps+1];
         lolos[ii][ii2] = new double[nReps+1];
         iEffects[ii][ii2] = new double[nReps+2];
         for (rr = 0; rr <= nReps; rr++)
         {
            iEffects[ii][ii2][rr] = 0.0;
            hihis[ii][ii2][rr] = 0.0;
            lolos[ii][ii2][rr] = 0.0;
            hilos[ii][ii2][rr] = 0.0;
            lohis[ii][ii2][rr] = 0.0;
         }
      }
   }

   printAsterisks(0);
   printf("* Fractional Factorial (2-level) Interaction Analysis\n");
   if (adata.samplingMethod_ == PSUADE_SAMP_FF4 ||
       adata.samplingMethod_ == PSUADE_SAMP_RFF4)
   {
      printf("* Note: Since Fractional Factorial Resolution 4 is used,\n");
      printf("* the first and second order effects are confounded.\n");
   }
   printDashes(0);

   for (rr = 0; rr < nReps; rr++)
   {
      offset = rr * nSamples / nReps;
      for (ii = 0; ii < nInputs; ii++)
      {
         for (ii2 = ii+1; ii2 < nInputs; ii2++)
         {
   
            for (ss = 0; ss < nSamples/nReps; ss++)
            {
               txArray[ss] = X[(offset+ss)*nInputs+ii];
               twArray[ss] = X[(offset+ss)*nInputs+ii2];
               tyArray[ss] = Y[(offset+ss)*nOutputs+whichOutput];
            }


            sortDbleList3(nSamples/nReps, txArray, twArray, tyArray);


            for (ss = 0; ss < nSamples/nReps; ss+=nSamples/nReps/2)
               sortDbleList2(nSamples/nReps/2,&twArray[ss],&tyArray[ss]);

            accum = 0.0;
            for (ss = 0; ss < nSamples/nReps/4; ss++) accum += tyArray[ss];
            lolos[ii][ii2][rr+1] += accum;
            accum = 0.0;
            for (ss = nSamples*3/nReps/4; ss < nSamples/nReps; ss++) accum += tyArray[ss];
            hihis[ii][ii2][rr+1] += accum;
            accum = 0.0;
            for (ss = nSamples/nReps/4; ss < nSamples/nReps/2; ss++) accum += tyArray[ss];
            lohis[ii][ii2][rr+1] += accum;
            accum = 0.0;
            for (ss = nSamples/nReps/2; ss < nSamples*3/nReps/4; ss++)
               accum += tyArray[ss];
            hilos[ii][ii2][rr+1] += accum;
            iEffects[ii][ii2][rr+1] = 0.5 * (lolos[ii][ii2][rr+1] +
                  hihis[ii][ii2][rr+1] - lohis[ii][ii2][rr+1] - hilos[ii][ii2][rr+1]);
         }
      }
      accum = 0.0;
      for (ii = 0; ii < nInputs; ii++)
         for (ii2 = ii+1; ii2 < nInputs; ii2++) accum += PABS(iEffects[ii][ii2][rr+1]);
      if (accum == 0.0) accum = 1.0;
      for (ii = 0; ii < nInputs; ii++)
         for (ii2 = ii+1; ii2 < nInputs; ii2++) 
            iEffects[ii][ii2][rr+1] = iEffects[ii][ii2][rr+1]/accum;
   }
   for (ii = 0; ii < nInputs; ii++)
   {
      for (ii2 = ii+1; ii2 < nInputs; ii2++)
      {
         mean = 0.0;
         for (rr = 0; rr < nReps; rr++) mean += iEffects[ii][ii2][rr+1];
         mean = mean / (double) nReps;
         iEffects[ii][ii2][0] = mean;
         stdev = 0.0;
         if (nReps > 1)
         {
            for (rr = 0; rr < nReps; rr++)
               stdev += pow(iEffects[ii][ii2][rr+1]-mean, 2.0);
            stdev = sqrt(stdev / (double) (nReps - 1.0));
         }
         iEffects[ii][ii2][nReps+1] = stdev/sqrt(1.0*nReps);
         if (nReps == 1)
            printf("* Input %3d %3d =  %12.4e\n",ii+1,ii2+1, mean);
         else
            printf("* Input %3d %3d =  %12.4e (std err = %12.4e)\n",ii+1, 
                   ii2+1, mean, stdev/sqrt(1.0*nReps));
      }
   }
   printAsterisks(0);

   if (fp != NULL)
   {
      ncount = (nInputs - 1) * nInputs / 2;
      if      (ncount == 1) nplots = 1; 
      else if (ncount == 2) nplots = 2; 
      else if (ncount <= 4) nplots = 4; 
      else if (ncount <= 8) nplots = 8; 
      else if (ncount <= 16) nplots = 16; 
      else if (ncount <= 32) nplots = 32; 
      else if (ncount <= 48) nplots = 48; 
      else                   nplots = 48;
      ncount = 0;
      if (psPlotTool_ == 1)
           fprintf(fp, "// matrix of interaction effect: ind1, ind2, effect\n");
      else fprintf(fp, "%% matrix of interaction effect: ind1, ind2, effect\n");
      fprintf(fp, "A2 = [\n");
      for (ii = 0; ii < nInputs; ii++)
         for (ii2 = ii+1; ii2 < nInputs; ii2++)
            fprintf(fp, "%3d %3d %12.4e %12.4e\n",ii+1,ii2+1, 
                   iEffects[ii][ii2][0], iEffects[ii][ii2][nReps+1]);
      fprintf(fp, "];\n");
      for (ii = 0; ii < nInputs; ii++)
      {
         for (ii2 = ii+1; ii2 < nInputs; ii2++)
         {
            fprintf(fp, "Y1 = [\n");
            for (rr = 0; rr < nReps; rr++)
            {
               fprintf(fp, "   %24.16e\n", lolos[ii][ii2][rr+1]);
               fprintf(fp, "   %24.16e\n", hilos[ii][ii2][rr+1]);
            }
            fprintf(fp, "];\n");
            fprintf(fp, "Y2 = [\n");
            for (rr = 0; rr < nReps; rr++)
            {
               fprintf(fp, "   %24.16e\n", lohis[ii][ii2][rr+1]);
               fprintf(fp, "   %24.16e\n", hihis[ii][ii2][rr+1]);
            }
            fprintf(fp, "];\n");
            if (nplots == 2)
            {
               fprintf(fp,"subplot(1,2,%d),",ncount+1);
               ncount++;
            }
            if (nplots == 4)
            {
               fprintf(fp,"subplot(2,2,%d),",ncount+1);
               ncount++;
            }
            if (nplots == 8)
            {
               fprintf(fp,"subplot(2,4,%d),",ncount+1);
               ncount++;
            }
            if (nplots == 16)
            {
               fprintf(fp,"subplot(4,4,%d),",ncount+1);
               ncount++;
            }
            if (nplots == 32 || nplots == 48)
            {
               fprintf(fp,"subplot(4,4,%d),",ncount+1);
               ncount++;
               if (ncount >= 16) ncount = 0;
            }
            fprintf(fp, "plot(Y1,'k')\n");
            if (psPlotTool_ == 1)
            {
               fprintf(fp, "set(gca(),\"auto_clear\",\"off\")\n");
               fprintf(fp, "plot(Y2,'b')\n");
               sprintf(pString, "Interaction(%d,%d)",ii+1,ii2+1);
               fwritePlotTitle(fp, pString);
               fprintf(fp, "set(gca(),\"auto_clear\",\"on\")\n");
            }
            else
            {
               fprintf(fp, "hold on\n");
               fprintf(fp, "plot(Y2,'b')\n");
               fprintf(fp, "hold off\n");
               fprintf(fp, "title('Interaction(%d,%d)')\n",ii+1,ii2+1);
               fprintf(fp, "text(0.2,0.2,'black: P1 (lo to hi), P2 = lo','sc')\n");
               fprintf(fp, "text(0.2,0.3,'blue:  P1 (lo to hi), P2 = hi','sc')\n");
            }
         }
      }
      fclose(fp);
      if (psPlotTool_ == 1)
         printf("The main effect plot has been generated in scilabffme.sci.\n");
      else
         printf("The main effect plot has been generated in matlabffme.m.\n");
   }
   printAsterisks(0);

   delete [] txArray;
   delete [] twArray;
   delete [] tyArray;
   for (ii = 0; ii < nInputs; ii++)
   {
      delete [] mEffects[ii];
      for (ii2 = 0; ii2 < nInputs; ii2++)
      {
         delete [] iEffects[ii][ii2];
         delete [] lolos[ii][ii2];
         delete [] lohis[ii][ii2];
         delete [] hilos[ii][ii2];
         delete [] hihis[ii][ii2];
      }
      delete [] iEffects[ii];
      delete [] lolos[ii];
      delete [] lohis[ii];
      delete [] hilos[ii];
      delete [] hihis[ii];
   }
   delete [] iArray;
   delete [] iEffects;
   delete [] lolos;
   delete [] lohis;
   delete [] hilos;
   delete [] hihis;
   delete [] mEffects;
   return 0.0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
FFAnalyzer& FFAnalyzer::operator=(const FFAnalyzer &)
{
   printf("FFAnalysis operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

