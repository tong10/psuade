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
#include "Util/sysdef.h"
#include "Util/PsuadeUtil.h"
#include "Main/Psuade.h"

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
   int    nInputs, nOutputs, nSamples, outputID, checkSample;
   int    ii, ii2, ss, whichOutput, ncount, nplots, *iArray;
   double *X, *Y, *txArray, *twArray, *tyArray, *mEffect, accum, **iEffect;
   double **lolo, **lohi, **hilo, **hihi;
   char   *cString, winput1[500], winput2[500], meFileName[500], pString[500];
   FILE   *fp=NULL;

   nInputs  = adata.nInputs_;
   nOutputs = adata.nOutputs_;
   nSamples = adata.nSamples_;
   outputID = adata.outputID_;
   X        = adata.sampleInputs_;
   Y        = adata.sampleOutputs_;
   if (adata.inputPDFs_ != NULL)
   {
      printf("FFAnalyzer INFO: some inputs have non-uniform PDFs, but they\n");
      printf("                 will not be relevant in this analysis.\n");
   }
   whichOutput = outputID;
   if (whichOutput >= nOutputs || whichOutput < 0) whichOutput = 0;

   if (nInputs <= 0 || nOutputs <= 0 || nSamples <= 0)
   {
      printf("FFAnalyzer ERROR: invalid arguments.\n");
      printf("   nInputs  = %d\n", nInputs);
      printf("   nOutputs = %d\n", nOutputs); 
      printf("   nSamples = %d\n", nSamples);
      return PSUADE_UNDEFINED;
   } 
   if (nSamples / 2 * 2 != nSamples)
   {
      printf("FFAnalyzer ERROR: nSamples has to be even.\n");
      printf("   nSamples = %d\n", nSamples);
      return PSUADE_UNDEFINED;
   } 

   if (psPlotTool_ == 0)
   {
      if (psAnaExpertMode_ == 1)
      {
         sprintf(pString,"Create FF main effect plot ? (y or n) ");
         getString(pString, winput1);
         if (winput1[0] == 'y')
         {
            sprintf(pString,"Enter scatter plot matlab file name : ");
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
            printf("FFAnalyzer: main effect matlab file = %s\n",
                   meFileName);
         }
      }
      else
      {
         if (psConfig_ != NULL)
         {
            cString = psConfig_->getParameter("FF_matlab_file");
            if (cString != NULL)
            {
               sscanf(cString, "%s %s %s",winput1,winput2,meFileName);
               fp = fopen(meFileName, "w");
               if (fp != NULL)
               {
                  printf("FFAnalyzer: main effect matlab file = %s\n",
                         meFileName);
               }
            }
         }
      }
   }

   txArray = new double[nSamples];
   tyArray = new double[nSamples];
   twArray = new double[nSamples];
   iArray  = new int[nInputs];
   mEffect = new double[nInputs];
   for (ii = 0; ii < nInputs; ii++) mEffect[ii] = 0.0;
   for (ii = 0; ii < nInputs; ii++) iArray[ii] = ii;

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
   if (fp != NULL)
   {
      fprintf(fp,"%% ********************************************** \n");
      fprintf(fp,"%% ********************************************** \n");
      fprintf(fp,"%% * Fractional Factorial Main Effect Analysis ** \n");
      fprintf(fp,"%% * (2 level fractional factorial only)       ** \n");
      fprintf(fp,"%% *-------------------------------------------** \n");
      fprintf(fp,"%% Output %d\n", whichOutput+1);
      fprintf(fp,"%% *-------------------------------------------** \n");
   }

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
         printf("FFAnalyzer ERROR: sample not fractional factorial.\n");
         delete [] iArray;
         delete [] mEffect;
         delete [] twArray;
         delete [] txArray;
         delete [] tyArray;
         return 0.0;
      }
      accum = 0.0;
      for (ss = nSamples/2; ss < nSamples; ss++) accum += tyArray[ss];
      for (ss = 0; ss < nSamples/2; ss++) accum -= tyArray[ss];
      mEffect[ii] = accum * 2 / (double) nSamples;
   }
   accum = 0.0;
   for (ii = 0; ii < nInputs; ii++) accum += PABS(mEffect[ii]);
   if (accum > 0.0)
      for (ii = 0; ii < nInputs; ii++) mEffect[ii] /= accum;
   printf("* Fractional Factorial Main Effect (normalized)\n");
   printDashes(0);
   for (ii = 0; ii < nInputs; ii++)
      printf("* Input %3d  =  %12.4e\n", ii+1, mEffect[ii]);
   printDashes(0);
   printf("* Fractional Factorial Main Effect (ordered)\n");
   printDashes(0);
   for (ii = 0; ii < nInputs; ii++) twArray[ii] = PABS(mEffect[ii]);
   sortDbleList2a(nInputs, twArray, iArray);
   for (ii = nInputs-1; ii >= 0; ii--)
      printf("* Rank %3d : Input %3d (measure = %12.4e)\n", nInputs-ii, 
             iArray[ii]+1, twArray[ii]);

   if (fp != NULL)
   {
      fprintf(fp, "Y = [\n");
      for (ii = 0; ii < nInputs; ii++)
         fprintf(fp, "%24.16e \n", mEffect[ii]);
      fprintf(fp, "]; \n");
      fprintf(fp, "clf\n");
      fprintf(fp, "hold off\n");
      fprintf(fp, "plot(Y,'*'); \n");
      fprintf(fp, "grid\n");
      fprintf(fp, "title('Main Effect Plot for Output %d')\n",whichOutput+1);
      fprintf(fp, "disp('Press enter to continue')\n");
      fprintf(fp, "pause\n");
   }

   if (nInputs < 2 || adata.samplingMethod_ == PSUADE_SAMP_PBD) return 0.0;

   iEffect = new double*[nInputs];
   hihi = new double*[nInputs];
   lolo = new double*[nInputs];
   lohi = new double*[nInputs];
   hilo = new double*[nInputs];
   for (ii = 0; ii < nInputs; ii++)
   {
      iEffect[ii] = new double[nInputs];
      hihi[ii] = new double[nInputs];
      lolo[ii] = new double[nInputs];
      lohi[ii] = new double[nInputs];
      hilo[ii] = new double[nInputs];
   }

   printAsterisks(0);
   printf("* Fractional Factorial (2-level) Interaction Analysis\n");
   printf("* The measures are normalized so that the sum = 1.\n");
   if (adata.samplingMethod_ == PSUADE_SAMP_FF4)
   {
      printf("* Note: Since Fractional Factorial Resolution 4 is used,\n");
      printf("* the first and second order effects are confounded.\n");
   }
   printDashes(0);

   for (ii = 0; ii < nInputs; ii++)
   {
      for (ii2 = ii+1; ii2 < nInputs; ii2++)
      {

         for (ss = 0; ss < nSamples; ss++)
         {
            txArray[ss] = X[nInputs*ss+ii];
            twArray[ss] = X[nInputs*ss+ii2];
            tyArray[ss] = Y[nOutputs*ss+whichOutput];
         }


         sortDbleList3(nSamples, txArray, twArray, tyArray);


         for (ss = 0; ss < nSamples; ss+=nSamples/2)
            sortDbleList2(nSamples/2,&twArray[ss],&tyArray[ss]);

         accum = 0.0;
         for (ss = 0; ss < nSamples/4; ss++) accum += tyArray[ss];
         lolo[ii][ii2] = accum;
         accum = 0.0;
         for (ss = nSamples*3/4; ss < nSamples; ss++) accum += tyArray[ss];
         hihi[ii][ii2] = accum;
         accum = 0.0;
         for (ss = nSamples/4; ss < nSamples/2; ss++) accum += tyArray[ss];
         lohi[ii][ii2] = accum;
         accum = 0.0;
         for (ss = nSamples/2; ss < nSamples*3/4; ss++)
            accum += tyArray[ss];
         hilo[ii][ii2] = accum;
         iEffect[ii][ii2] = 0.5 * (lolo[ii][ii2] +
                  hihi[ii][ii2] - lohi[ii][ii2] - hilo[ii][ii2]);
      }
   }
   accum = 0.0;
   for (ii = 0; ii < nInputs; ii++)
      for (ii2 = ii+1; ii2 < nInputs; ii2++) accum += PABS(iEffect[ii][ii2]);
   if (accum > 0.0)
   {
      for (ii = 0; ii < nInputs; ii++)
      {
         for (ii2 = ii+1; ii2 < nInputs; ii2++)
         {
            printf("* Input %3d %3d =  %12.4e\n",ii+1,ii2+1, 
                   PABS(iEffect[ii][ii2]/accum));
         }
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
      for (ii = 0; ii < nInputs; ii++)
      {
         for (ii2 = ii+1; ii2 < nInputs; ii2++)
         {
            fprintf(fp, "Y1 = [\n");
            fprintf(fp, "   %24.16e\n", lolo[ii][ii2]);
            fprintf(fp, "   %24.16e\n", hilo[ii][ii2]);
            fprintf(fp, "];\n");
            fprintf(fp, "Y2 = [\n");
            fprintf(fp, "   %24.16e\n", lohi[ii][ii2]);
            fprintf(fp, "   %24.16e\n", hihi[ii][ii2]);
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
               if (ncount >= 16)
               {
                  ncount = 0;
                  fprintf(fp, "hold off\n");
               }
            }
            fprintf(fp, "plot(Y1,'k')\n");
            fprintf(fp, "hold on\n");
            fprintf(fp, "plot(Y2,'b')\n");
            fprintf(fp, "hold off\n");
            fprintf(fp, "title('Interaction(%d,%d)')\n",ii+1,ii2+1);
            fprintf(fp, "text(0.2,0.2,'black: P1 (lo to hi), P2 = lo','sc')\n");
            fprintf(fp, "text(0.2,0.3,'blue:  P1 (lo to hi), P2 = hi','sc')\n");
            fprintf(fp, "disp('Press enter to continue') \n");
            fprintf(fp, "pause \n");
         }
      }
      fclose(fp);
      printf("The %s main effect plot has been generated.\n",meFileName);
   }
   printAsterisks(0);
   printAsterisks(0);

   delete [] txArray;
   delete [] twArray;
   delete [] tyArray;
   delete [] mEffect;
   for (ii = 0; ii < nInputs; ii++)
   {
      delete [] iEffect[ii];
      delete [] lolo[ii];
      delete [] lohi[ii];
      delete [] hilo[ii];
      delete [] hihi[ii];
   }
   delete [] iArray;
   delete [] iEffect;
   delete [] lolo;
   delete [] lohi;
   delete [] hilo;
   delete [] hihi;
   return 0.0;
}

