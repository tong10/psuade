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
// Functions for the class RSFuncApproxAnalyzer  
// AUTHOR : CHARLES TONG
// DATE   : 2003
// ************************************************************************
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "Psuade.h"
#include "FuncApprox.h"
#include "Ann.h"
#include "Regression.h"
#include "UserRegression.h"
#include "RSFuncApproxAnalyzer.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
RSFuncApproxAnalyzer::RSFuncApproxAnalyzer() : Analyzer()
{
   setName("RSFA");
   rsType_ = PSUADE_RS_MARS;
   useCV_ = 0;
   numCVGroups_ = 10;
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
RSFuncApproxAnalyzer::~RSFuncApproxAnalyzer() 
{ 
} 

// ************************************************************************ 
// perform analysis 
// ------------------------------------------------------------------------
double RSFuncApproxAnalyzer::analyze(aData &adata)
{
   int        nInputs, nOutputs, nSamples, outputID, nLevels, *levelSeps;
   int        iD, iL, nLast, nPtsPerDim=64, length, status, rsState;
   int        count, iD2, iI, nSubSamples, wgtID, printLevel;
   int        *iArray, *iArray2, testFlag, iOne=1;
   double     ddata, ymax, ymin, *YLocal, *X, *Y, *X2, *Y2, retdata=0;
   double     *lower, *upper, *eArray, *YT, *WW, *wgts, sdata, *XX, *YY;
   double     cvErr1, cvErr1s, cvErr2, cvErr2s,cvMax, cvMaxs, cvMaxBase;
   double     cvMaxBases, CVMaxBases, maxBase, maxBases, *sArray, *sigmas;
   double     CVErr2, CVErr2s, CVErr1, CVErr1s, CVMax, CVMaxs, CVMaxBase;
   double     sumErr1, sumErr1s, sumErr2, sumErr2s, maxErr, maxErrs;
   double     sumErr11, sumErr11s, ymean, yvar, ydiff, *S2, ssum;
   char       *cString;
   char       winput1[500], winput2[500], dataFile[500], errFile[500];
   char       pString[500];
   FILE       *fpData, *fpErr;
   FuncApprox *fa;

   printLevel = adata.printLevel_;
   nInputs   = adata.nInputs_;
   nOutputs  = adata.nOutputs_;
   nSamples  = adata.nSamples_;
   lower     = adata.iLowerB_;
   upper     = adata.iUpperB_;
   X         = adata.sampleInputs_;
   Y         = adata.sampleOutputs_;
   outputID  = adata.outputID_;
   nLevels   = adata.currRefineLevel_;
   levelSeps = adata.refineSeparators_;
   wgtID     = adata.regWgtID_;
   if (adata.inputPDFs_ != NULL)
   {
      count = 0;
      for (iI = 0; iI < nInputs; iI++) count += adata.inputPDFs_[iI];
      if (count > 0)
      {
         printf("RSAnalysis INFO: some inputs have non-uniform PDFs, but\n");
         printf("                they are not relevant in this analysis.\n");
      }
   }

   if (nInputs <= 0 || nOutputs <= 0 || nSamples <= 0)
   {
      printf("RSAnalysis ERROR: invalid arguments.\n");
      printf("   nInputs  = %d\n", nInputs);
      printf("   nOutputs = %d\n", nOutputs);
      printf("   nSamples = %d\n", nSamples);
      return PSUADE_UNDEFINED;
   } 
   status = 0;
   for (iD = 0; iD < nSamples; iD++)
      if (Y[nOutputs*iD+outputID] > 0.9*PSUADE_UNDEFINED) status = 1;
   if (status == 1)
   {
      printf("RSAnalysis ERROR: Some outputs are undefined. Prune the\n");
      printf("                  undefined sample points first.\n");
      return PSUADE_UNDEFINED;
   }
   if (printLevel > 0)
   {
      printAsterisks(0);
      printThisFA(rsType_);
      printEquals(0);
   }
   
   YLocal = new double[nSamples];
   for (iD = 0; iD < nSamples; iD++) YLocal[iD] = Y[iD*nOutputs+outputID];
   ymax = 0.0;
   ymin = PSUADE_UNDEFINED;
   for (iD = 0; iD < nSamples; iD++)
   {
      if (PABS(YLocal[iD]) > ymax) ymax = PABS(YLocal[iD]);
      if (PABS(YLocal[iD]) < ymin) ymin = PABS(YLocal[iD]);
   }
   printf("RSA: Output ID = %d\n", outputID+1);
   printf("RSA: Output Maximum/Minimum = %14.6e %14.6e\n",ymax,ymin);
   printf("INFO: Turn on analysis expert mode to perform cross validation.\n");
   printf("INFO: Set printlevel higher (1-4) to display more information.\n");
   printf("INFO: Set print level to 4 for interpolation error graphics file.\n");
   printf("INFO: Set analysis expert mode to perform cross validation.\n");
   if (ymax == PSUADE_UNDEFINED)
   {
      printf("RSAnalysis ERROR: some outputs are undefined.\n");
      printf("                  Prune them first before analyze.\n");
      delete [] YLocal;
      return PSUADE_UNDEFINED;
   }

   YT = new double[nSamples];
   for (iL = 0; iL < nLevels; iL++)
   {
      nLast = nSamples;
      if (levelSeps != NULL) nLast = levelSeps[iL];

      fa = genFA(rsType_, nInputs, iOne, nLast);
      if (fa == NULL)
      {
         printf("RSFAnalysis INFO: cannot create function approximator.\n");
         delete [] YLocal;
         delete fa;
         delete [] YT;
         return 1.0;
      }
      fa->setNPtsPerDim(nPtsPerDim);
      fa->setBounds(lower, upper);
      if (iL == nLevels-1) fa->setOutputLevel(printLevel);
      if (wgtID >= 0 && wgtID < nOutputs)
      {
         wgts = new double[nLast];
         for (iD = 0; iD < nLast; iD++) wgts[iD] = Y[iD*nOutputs+wgtID];
         fa->loadWeights(nLast, wgts);
         delete [] wgts;
      }
      length = -999;
      status = fa->genNDGridData(X, YLocal, &length, NULL, NULL);
      if (status != 0)
      {
         printf("RSAnalysis ERROR: something wrong in genNDGridData.\n");
         delete [] YLocal;
         delete fa;
         delete [] YT;
         return PSUADE_UNDEFINED;
      }
      if (rsType_ == PSUADE_RS_GP2)
      {
         printf("RSAnalysis INFO: GP2 currently for checking only.\n");
         delete [] YLocal;
         delete [] YT;
         delete fa;
         return 1.0;
      }

      sumErr1  = sumErr2 = maxErr = sumErr1s = sumErr2s = maxErrs = 0.0;
      sumErr11 = sumErr11s = ydiff = 0.0;
      maxBase = maxBases = ymean = yvar = 0.0;
      fa->evaluatePoint(nLast, X, YT);
      for (iD = 0; iD < nLast; iD++)
      {
         ddata = PABS(YT[iD] - YLocal[iD]);
         if (YLocal[iD] != 0.0) sdata = ddata / PABS(YLocal[iD]);
         else                   sdata = ddata;
         sumErr1   += ddata;
         sumErr1s  += sdata;
         sumErr11  += (YT[iD] - YLocal[iD]);
         if (YLocal[iD] != 0.0)
            sumErr11s += (YT[iD] - YLocal[iD]) / PABS(YLocal[iD]);
         else
            sumErr11s += (YT[iD] - YLocal[iD]);
         if (ddata > maxErr ) {maxErr = ddata;  maxBase = PABS(YLocal[iD]);}
         if (sdata > maxErrs) {maxErrs = sdata; maxBases = PABS(YLocal[iD]);}
         sumErr2  += (ddata * ddata);
         sumErr2s += (sdata * sdata);
         ymean += YLocal[iD]; 
      }
      ymean /= (double) nLast;
      for (iD = 0; iD < nLast; iD++)
         yvar += (YLocal[iD] - ymean) * (YLocal[iD] - ymean);
      sumErr1   = sumErr1 / (double) nLast;
      sumErr1s  = sumErr1s / (double) nLast;
      sumErr11  = sumErr11 / (double) nLast;
      sumErr11s = sumErr11s / (double) nLast;
      sumErr2  = sqrt(sumErr2 / (double) nLast);
      sumErr2s = sqrt(sumErr2s / (double) nLast);
      if (printLevel > 2 || iL == nLevels-1)
      {
         printf("RSAnalysis: L %2d: interpolation error on training set \n", iL);
         printf("             avg error far from 0 ==> systematic bias.\n");
         printf("             rms error large      ==> average   error large.\n");
         printf("             max error large      ==> pointwise error large.\n");
         printf("             R-square may not always be a reliable measure.\n");
         printf("  avg error   = %11.3e (unscaled)\n", sumErr11);
         printf("  avg error   = %11.3e (scaled)\n", sumErr11s);
         printf("  rms error   = %11.3e (unscaled)\n", sumErr2);
         printf("  rms error   = %11.3e (scaled)\n", sumErr2s);
         printf("  max error   = %11.3e (unscaled, BASE=%9.3e)\n",
                maxErr, maxBase);
         printf("  max error   = %11.3e (  scaled, BASE=%9.3e)\n",
                maxErrs, maxBases);
         printf("  R-square    = %16.8e\n",1.0 - sumErr2*sumErr2*nLast / yvar);
         printf("Based on %d training points (total=%d).\n",nLast,nSamples);
      }

      if (nLevels > 1 && iL < nLevels-1)
      {
         sumErr1 = sumErr2 = maxErr = sumErr1s = sumErr2s = maxErrs = 0.0;
         sumErr11 = sumErr11s;
         fa->evaluatePoint(nSamples, X, YT);
         for (iD = nLast; iD < nSamples; iD++)
         {
            ddata = PABS(YT[iD] - YLocal[iD]);
            if (YLocal[iD] != 0.0) sdata = ddata / PABS(YLocal[iD]);
            else                   sdata = ddata;
            sumErr1  += ddata;
            sumErr1s += sdata;
            sumErr11  += (YT[iD] - YLocal[iD]);
            if (YLocal[iD] != 0.0)
               sumErr11s += (YT[iD] - YLocal[iD]) / PABS(YLocal[iD]);
            else
               sumErr11s += (YT[iD] - YLocal[iD]);
            if (ddata > maxErr ) {maxErr  = ddata; maxBase = PABS(YLocal[iD]);}
            if (sdata > maxErrs) {maxErrs = sdata; maxBases = PABS(YLocal[iD]);}
            sumErr2  += (ddata * ddata);
            sumErr2s += (sdata * sdata);
         }
         sumErr1   = sumErr1 / (double) (nSamples-nLast);
         sumErr1s  = sumErr1s / (double) (nSamples-nLast);
         sumErr11  = sumErr11 / (double) (nSamples-nLast);
         sumErr11s = sumErr11s / (double) (nSamples-nLast);
         sumErr2   = sqrt(sumErr2 / (double) (nSamples-nLast));
         sumErr2s  = sqrt(sumErr2s / (double) (nSamples-nLast));
         if (printLevel > 1) 
         {
            printf("RSAnalysis: L %2d: Prediction error on the remaining data:\n",
                   iL);
            printf("  avg error = %11.3e (unscaled)\n", sumErr11);
            printf("  avg error = %11.3e (scaled)\n", sumErr11s);
            printf("  rms error = %11.3e (unscaled)\n", sumErr2);
            printf("  rms error = %11.3e (scaled)\n", sumErr2s);
            printf("  max error = %11.3e (unscaled, BASE=%11.3e)\n", maxErr, 
                   maxBase);
            printf("  max error = %11.3e (  scaled, BASE=%11.3e)\n", maxErrs, 
                   maxBases);
            printf("Based on %d training points (rest=%d).\n",nLast,
                   nSamples-nLast);
         }
      }

      if (printLevel > 3 && iL == nLevels-1) 
      {
         if (psPlotTool_ == 1)
         {
            fpErr = fopen("psuade_rsfa_err.sci", "w");
            if (fpErr != NULL)
            {
               fprintf(fpErr, "// Surface Fitting Error Histogram \n");
               fprintf(fpErr, "// Interpolation errors on all points\n");
               fprintf(fpErr, "// col 1: interpolated data\n");
               fprintf(fpErr, "// col 2: training data\n");
               fprintf(fpErr, "// col 3: col 1 - col2\n");
               fprintf(fpErr, "// col 4 on: input data\n");
            }
            else printf("INFO: cannot open file psuade_rsfa_err.sci.\n");
         }
         else
         {
            fpErr = fopen("psuade_rsfa_err.m", "w");
            if (fpErr != NULL)
            {
               fprintf(fpErr, "%% Surface Fitting Error Histogram \n");
               fprintf(fpErr, "%% Interpolation errors on all points\n");
               fprintf(fpErr, "%% col 1: interpolated data\n");
               fprintf(fpErr, "%% col 2: training data\n");
               fprintf(fpErr, "%% col 3: col 1 - col2\n");
               fprintf(fpErr, "%% col 4 on: input data\n");
            }
            else printf("INFO: cannot open file psuade_rsfa_err.m.\n");
         }
         if (fpErr != NULL) fprintf(fpErr, "E = [\n");
         sumErr1 = sumErr2 = maxErr = sumErr1s = sumErr2s = maxErrs = 0.0;
         maxBase = sumErr11 = sumErr11s = 0.0;
         fa->evaluatePoint(nSamples, X, YT);
         for (iD = 0; iD < nSamples; iD++)
         {
            ddata = YT[iD];
            if (fpErr != NULL)
            {
               fprintf(fpErr, "%24.16e %24.16e %24.16e ", ddata, 
                       YLocal[iD], ddata-YLocal[iD]);
               for (iI = 0; iI < nInputs; iI++)
                  fprintf(fpErr, "%24.16e ", X[iD*nInputs+iI]); 
               fprintf(fpErr, "\n");
            }
            ddata = ddata - YLocal[iD];
            sumErr11 += ddata;
            if (YLocal[iD] != 0.0) sumErr11s += ddata / PABS(YLocal[iD]);
            else                   sumErr11s += ddata;
            ddata = PABS(ddata);
            if (YLocal[iD] != 0.0) sdata = ddata / PABS(YLocal[iD]);
            else                   sdata = ddata;
            sumErr1  += ddata;
            sumErr1s += sdata;
            if (ddata > maxErr ) {maxErr  = ddata; maxBase = PABS(YLocal[iD]);}
            if (sdata > maxErrs) {maxErrs = sdata; maxBases = PABS(YLocal[iD]);}
            sumErr2  += (ddata * ddata);
            sumErr2s += (sdata * sdata);
         }
         if (fpErr != NULL)
         {
            fprintf(fpErr, "];\n");
            fwritePlotCLF(fpErr);
            fwritePlotFigure(fpErr, 1);
            fprintf(fpErr, "subplot(1,2,1)\n");
            fprintf(fpErr, "plot(E(:,3),'x')\n");
            fwritePlotAxes(fpErr);
            fwritePlotTitle(fpErr, "Interpolation Error Plot");
            fwritePlotXLabel(fpErr, "Sample Number");
            fwritePlotYLabel(fpErr, "Interpolation Error");
            fprintf(fpErr, "subplot(1,2,2)\n");
            fprintf(fpErr, "xmax = max(E(:,2));\n");
            fprintf(fpErr, "xmin = min(E(:,2));\n");
            fprintf(fpErr, "ymax = max(E(:,1));\n");
            fprintf(fpErr, "ymin = min(E(:,1));\n");
            fprintf(fpErr, "xmin = min(xmin, ymin);\n");
            fprintf(fpErr, "xmax = max(xmax, ymax);\n");
            fprintf(fpErr, "XX   = xmin : xmax-xmin : xmax;\n");
            fprintf(fpErr, "plot(E(:,2), E(:,1),'x', XX, XX)\n");
            if (psPlotTool_ == 1)
            {
               fprintf(fpErr, "a = get(\"current_axes\");\n");
               fprintf(fpErr, "a.data_bounds=[xmin,xmin;xmax,xmax];\n");
            }
            else
            {
               fprintf(fpErr, "axis([xmin xmax xmin xmax])\n");
            }
            fwritePlotAxes(fpErr);
            fwritePlotTitle(fpErr, "Interpolated versus actual data");
            fwritePlotXLabel(fpErr, "Actual data");
            fwritePlotYLabel(fpErr, "Interpolated data");

            fwritePlotFigure(fpErr, 2);
            fwritePlotCLF(fpErr);
            fprintf(fpErr, "subplot(1,2,1)\n");
            if (psPlotTool_ == 1)
            {
               fprintf(fpErr, "ymin = min(E(:,3));\n");
               fprintf(fpErr, "ymax = max(E(:,3));\n");
               fprintf(fpErr, "ywid = 0.1 * (ymax - ymin);\n");
               fprintf(fpErr, "if (ywid < 1.0e-12)\n");
               fprintf(fpErr, "   disp('range too small.')\n");
               fprintf(fpErr, "   halt\n");
               fprintf(fpErr, "end;\n");
               fprintf(fpErr, "histplot(10, E(:,3), style=2);\n");
               fprintf(fpErr, "a = gce();\n");
               fprintf(fpErr, "a.children.fill_mode = \"on\";\n");
               fprintf(fpErr, "a.children.thickness = 2;\n");
               fprintf(fpErr, "a.children.foreground = 0;\n");
               fprintf(fpErr, "a.children.background = 2;\n");
            }
            else
            {
               fprintf(fpErr, "[nk,xk]=hist(E(:,3),10);\n");
               fprintf(fpErr, "bar(xk,nk/%d,1.0)\n",nSamples);
            }
            fwritePlotAxes(fpErr);
            fwritePlotTitle(fpErr, "Interpolation Errors Histogram");
            fwritePlotXLabel(fpErr, "Error");
            fwritePlotYLabel(fpErr, "Probabilities");
            fprintf(fpErr, "subplot(1,2,2)\n");
            if (psPlotTool_ == 1)
            {
               fprintf(fpErr, "ymin = min(E(:,3)/E(:,2));\n");
               fprintf(fpErr, "ymax = max(E(:,3)/E(:,2));\n");
               fprintf(fpErr, "ywid = 0.1 * (ymax - ymin);\n");
               fprintf(fpErr, "if (ywid < 1.0e-12)\n");
               fprintf(fpErr, "   disp('range too small.')\n");
               fprintf(fpErr, "   halt\n");
               fprintf(fpErr, "end;\n");
               fprintf(fpErr, "histplot(10, E(:,3)/E(:,2), style=2);\n");
               fprintf(fpErr, "a = gce();\n");
               fprintf(fpErr, "a.children.fill_mode = \"on\";\n");
               fprintf(fpErr, "a.children.thickness = 2;\n");
               fprintf(fpErr, "a.children.foreground = 0;\n");
               fprintf(fpErr, "a.children.background = 2;\n");
            }
            else
            {
               fprintf(fpErr, "[nk,xk]=hist(E(:,3)./E(:,2),10);\n");
               fprintf(fpErr, "bar(xk,nk/%d,1.0)\n",nSamples);
            }
            fwritePlotAxes(fpErr);
            fwritePlotTitle(fpErr, "Interpolation Errors Histogram (normalized)");
            fwritePlotXLabel(fpErr, "Error");
            fwritePlotYLabel(fpErr, "Probabilities");

            fwritePlotFigure(fpErr, 3);
            fprintf(fpErr, "nn = %d;\n", nInputs);
            iL = (int) pow(1.0*nInputs, 0.5);
            if (iL*iL < nInputs) iL++;
            if (psPlotTool_ == 1) fprintf(fpErr, "drawlater\n");
            fprintf(fpErr, "for ii = 1 : nn\n");
            fprintf(fpErr, "   subplot(%d,%d,ii)\n", iL, iL);
            fprintf(fpErr, "   plot(E(:,ii+3),E(:,3),'x')\n");
            fwritePlotAxes(fpErr);
            if (psPlotTool_ == 1)
            {
               fwritePlotTitle(fpErr,"Error Plot for Input");
               sprintf(winput1, 
                    "a.title.text = \"Error Plot for Input\" + string(ii);\n");
               fprintf(fpErr, "%s", winput1);
               fwritePlotXLabel(fpErr, "Input Values");
            }
            else
            {
               fprintf(fpErr, "title(['Error Plot for Input ',int2str(ii)])\n");
               fwritePlotXLabel(fpErr, "Input Values");
            }
            fwritePlotYLabel(fpErr, "Interpolation Error");
            if (psPlotTool_ == 1) fprintf(fpErr, "drawnow\n");
            fprintf(fpErr, "end\n");

            if (nInputs > 2)
            {
               fwritePlotFigure(fpErr, 4);
               if (psPlotTool_ == 1)
               {
                  fprintf(fpErr, "f = gcf();\n");
                  fprintf(fpErr, "f.color_map = jetcolormap(3);\n");
                  fprintf(fpErr, "drawlater\n");
                  fprintf(fpErr, "param3d1([E(:,4)' ; E(:,4)'],");
                  fprintf(fpErr, "[E(:,5)' ; E(:,5)'],[E(:,2)' ; E(:,2)'])\n");
                  fprintf(fpErr, "e = gce();\n");
                  fprintf(fpErr, "e.children.mark_mode = \"on\";\n");
                  fprintf(fpErr, "e.children.mark_size_unit = \"point\";\n");
                  fprintf(fpErr, "e.children.mark_style = 10;\n");
                  fprintf(fpErr, "e.children.mark_size = 6;\n");
                  fprintf(fpErr, "for i = 1:length(e.children)\n");
                  fprintf(fpErr, "   e.children(i).mark_foreground = 1;\n");
                  fprintf(fpErr, "end\n");
                  fprintf(fpErr, "set(gca(),\"auto_clear\",\"off\")\n");
                  fprintf(fpErr, "drawnow\n");
               }
               else
               {
                  fprintf(fpErr, "   plot3(E(:,4),E(:,5),E(:,2),'bp')\n");
               }
               fwritePlotXLabel(fpErr, "Input 1");
               fwritePlotYLabel(fpErr, "Input 2");
               fwritePlotTitle(fpErr, "Output Error Scatter Plot");
               fwritePlotAxes(fpErr);
            }
            fclose(fpErr);
            if (psPlotTool_ == 1)
               printf("Interpolation error info are in psuade_rsfa_err.sci\n");
            else
               printf("Interpolation error info are in psuade_rsfa_err.m\n");
         }
         sumErr1   = sumErr1 / (double) nSamples;
         sumErr1s  = sumErr1s / (double) nSamples;
         sumErr11  = sumErr11 / (double) nSamples;
         sumErr11s = sumErr11s / (double) nSamples;
         sumErr2   = sqrt(sumErr2);
         sumErr2s  = sqrt(sumErr2s);
         retdata = sumErr1 / ymax;
      }
      delete fa;
   }
   delete [] YT;

   nSubSamples = nSamples / numCVGroups_;
#if 0
   if (psConfig_ != NULL)
   {
      cString = psConfig_->getParameter("RSFA_cross_validation");
      if (cString != NULL)
      {
         useCV_ = 1;
         printf("RSFA: turn on cross validation (from config file)\n");
      }
      cString = psConfig_->getParameter("RSFA_cv_ngroups");
      if (cString != NULL)
      {
         sscanf(cString, "%s %s %d",winput1,winput2,&iD);
         if (iD > 0)
         {
            printf("RSFA: number of CV groups = %d\n",iD);
            numCVGroups_ = iD;
            nSubSamples = nSamples / numCVGroups_;
         }
         else
         {
            printf("RSFA: invalid number of CV groups = %d\n",iD);
            useCV_ = 0;
         }
      }
   }
   else
#endif
   {
      printAsterisks(0);
      printf("Next you will be asked whether to do cross validation or not.\n");
      printf("Since cross validation iterates as many times as the number\n");
      printf("of groups. The rs_expert mode will be turned off. To change\n");
      printf("the default parameters for different response surface, you\n");
      printf("will need to exit, create a config file (use genconfigfile\n");
      printf("in command line mode), and set config option in your data file.\n");
      printDashes(0);
      sprintf(pString, "Perform cross validation ? (y or n) ");
      getString(pString, winput1);
      if (winput1[0] == 'y')
      {
         useCV_ = 1;
         sprintf(pString, "Enter the number of groups to validate : (2 - %d) ",
                 nSamples);
         iD = getInt(1, nSamples, pString);
         printf("RSFA: number of CV groups = %d\n",iD);
         numCVGroups_ = iD;
         nSubSamples = nSamples / numCVGroups_;
      }
   }

   if (useCV_ == 1)
   {
      printf("RSAnalysis: L %2d:cross validation (CV) begins...\n",nLevels);

      adata.sampleErrors_ = new double[nSamples];
      for (iD = 0; iD < nSamples; iD++) adata.sampleErrors_[iD] = 0.0;

      fa = genFA(rsType_, nInputs, iOne, nSamples-nSubSamples);
      if(fa == NULL)
      {
	 printf("RSAnalysis: genFA returned NULL in file %s line %d\n", 
                __FILE__, __LINE__ );
         exit(1);
      }
      fa->setNPtsPerDim(nPtsPerDim);
      fa->setBounds(lower, upper);
      fa->setOutputLevel(printLevel);
      rsState = psRSExpertMode_;
      psRSExpertMode_ = 0;

      XX = new double[nSamples*nInputs];
      YY = new double[nSamples];
      YT = new double[nSamples];
      S2 = new double[nSamples];
      WW = new double[nSamples];
      iArray  = new int[nSamples];
      iArray2 = new int[nSamples];
      eArray  = new double[nSamples];
      sArray  = new double[nSamples];
      sigmas  = new double[nSamples];
      sprintf(pString, "Random selection of leave-out groups ? (y or n) ");
      getString(pString, winput1);
      if (winput1[0] == 'y')
      {
         generateRandomIvector(nSamples, iArray);
      }
      else
      {
         for (iD = 0; iD < nSamples; iD++) iArray[iD] = iD;
      }
      for (iI = 0; iI < nInputs; iI++)
      {
         for (iD = 0; iD < nSamples; iD++)
            XX[iArray[iD]*nInputs+iI] = X[iD*nInputs+iI];
      }
      for (iD = 0; iD < nSamples; iD++) YY[iArray[iD]] = YLocal[iD];
      for (iD = 0; iD < nSamples; iD++) iArray2[iArray[iD]] = iD;
      if (wgtID >= 0 && wgtID < nOutputs)
      {
         for (iD = 0; iD < nSamples; iD++)
            WW[iArray[iD]] = Y[iD*nOutputs+wgtID];
      }

      X2 = new double[nSamples*nInputs];
      Y2 = new double[nSamples];
      wgts = new double[nSamples];
      CVErr1 = CVErr1s = CVErr2 = CVErr2s = CVMax = CVMaxs = 0.0;
      cvMaxBase = cvMaxBases = 0.0;
      for (iD = 0; iD < nSamples; iD+=nSubSamples)
      {
         printf("RSAnalysis:: L %2d:CV processes %d out of %d.\n",nLevels,
                iD/nSubSamples+1, nSamples/nSubSamples);

         count = 0;
         for (iD2 = 0; iD2 < nSamples; iD2++)
         {
            if (iD2 < iD || iD2 >= (iD+nSubSamples))
            {
               for (iI = 0; iI < nInputs; iI++)
                  X2[count*nInputs+iI] = XX[iD2*nInputs+iI];
               if (wgtID >= 0 && wgtID < nOutputs)
                  wgts[count] = WW[iD2*nOutputs+wgtID];
               Y2[count++] = YY[iD2];
            }
         }
         if (wgtID >= 0 && wgtID < nOutputs)
            fa->loadWeights(nSamples-nSubSamples, wgts);

         length = -999;
         status = fa->genNDGridData(X2, Y2, &length, NULL, NULL);
         if (status == -1) break;
         count = nSubSamples;
         if ((iD + nSubSamples) > nSamples) count = nSamples - iD;
         fa->evaluatePointFuzzy(count, &(XX[iD*nInputs]), YT, S2);

         cvErr1 = cvErr2 = cvErr1s = cvErr2s = cvMax = cvMaxs = 0.0;
         for (iD2 = 0; iD2 < count; iD2++)
         {
            ddata = YT[iD2] - YY[iD+iD2];
            eArray[iArray2[iD+iD2]] = ddata;
            sArray[iArray2[iD+iD2]] = YT[iD2];
            sigmas[iArray2[iD+iD2]] = S2[iD2];
            cvErr1  += ddata;
            cvErr2  += (ddata * ddata);
            if (PABS(ddata) > cvMax)
            {
               cvMax  = PABS(ddata);
               cvMaxBase = PABS(YY[iD+iD2]);
            }
            if (YY[iD+iD2] != 0.0) sdata = ddata / PABS(YY[iD+iD2]);
            else                   sdata = ddata;
            cvErr1s += sdata;
            cvErr2s += (sdata * sdata);
            if (PABS(sdata) > cvMaxs)
            {
               cvMaxs = PABS(sdata);
               cvMaxBases = PABS(YY[iD+iD2]);
            }
            if (printLevel > 4) 
               printf("Sample %6d: predicted =  %e, actual =  %e\n", 
                      iArray2[iArray[iD+iD2]], YT[iD2], YY[iD+iD2]);
         }
         adata.sampleErrors_[iD/nSubSamples] = cvErr2s;
         CVErr1  += cvErr1;
         CVErr1s += cvErr1s;
         CVErr2  += cvErr2;
         CVErr2s += cvErr2s;
         cvErr1  = cvErr1 / (double) nSubSamples;
         cvErr1s = cvErr1s / (double) nSubSamples;
         cvErr2  = sqrt(cvErr2 / nSubSamples);
         cvErr2s = sqrt(cvErr2s / nSubSamples);
         if (cvMax > CVMax ) {CVMax = cvMax; CVMaxBase = cvMaxBase;}
         if (cvMaxs > CVMaxs ) {CVMaxs = cvMaxs; CVMaxBases = cvMaxBases;}

         printf("RSA: first member of sample group %5d = %d)\n",
                iD/nSubSamples+1, iArray[iD]+1);
         printf("RSA: CV error for sample group %5d = %11.3e (avg unscaled)\n",
               iD/nSubSamples+1, cvErr1);
         printf("RSA: CV error for sample group %5d = %11.3e (avg scaled)\n",
               iD/nSubSamples+1, cvErr1s);
         printf("RSA: CV error for sample group %5d = %11.3e (rms unscaled)\n",
               iD/nSubSamples+1, cvErr2);
         printf("RSA: CV error for sample group %5d = %11.3e (rms scaled)\n",
               iD/nSubSamples+1, cvErr2s);
         printf("RSA: CV error for sample group %5d = %11.3e (max",
               iD/nSubSamples+1, cvMax);
         printf(" unscaled,BASE=%9.3e)\n", cvMaxBase);
         printf("RSA: CV error for sample group %5d = %11.3e (max",
               iD/nSubSamples+1, cvMaxs);
         printf("   scaled,BASE=%9.3e)\n", cvMaxBases);
      }
      if (status >= 0)
      {
         CVErr1  = CVErr1 / (double) nSamples;
         CVErr1s = CVErr1s / (double) nSamples;
         CVErr2  = sqrt(CVErr2 / nSamples);
         CVErr2s = sqrt(CVErr2s / nSamples);
         printf("RSA: final CV error  = %11.3e (avg unscaled)\n",CVErr1);
         printf("RSA: final CV error  = %11.3e (avg   scaled)\n",CVErr1s);
         printf("RSA: final CV error  = %11.3e (rms unscaled)\n",CVErr2);
         printf("RSA: final CV error  = %11.3e (rms   scaled)\n",CVErr2s);
         printf("RSA: final CV error  = %11.3e (max unscaled, BASE=%9.3e)\n",
                CVMax, CVMaxBase);
         printf("RSFA: final CV error  = %11.3e (max   scaled, BASE=%9.3e)\n",
                CVMaxs, CVMaxBases);
         printf("RSFA: L %2d:cross validation (CV) completed.\n",nLevels);
         if (psPlotTool_ == 1)
         {
            fpData = fopen("RSFA_CV_err.sci", "w");
            if (fpData != NULL)
            {
               fprintf(fpData,"// This file stores CV error for each point\n");
               fprintf(fpData,"// Column 1: error (true - predicted)\n");
               fprintf(fpData,"// Column 2: true\n");
               fprintf(fpData,"// Column 3: predicted\n");
               fprintf(fpData,"// Column 4: standard deviations\n"); 
               fprintf(fpData,"// Set morePlots=1 for normalized residual error plot.\n");
            }
            else printf("INFO: cannot open file RSFA_CV_err.sci.\n");
         }
         else
         {
            fpData = fopen("RSFA_CV_err.m", "w");
            if (fpData != NULL)
            {
               fprintf(fpData,"%% This file stores CV error for each point\n");
               fprintf(fpData,"%% Column 1: error (true - predicted)\n");
               fprintf(fpData,"%% Column 2: true\n");
               fprintf(fpData,"%% Column 3: predicted\n");
               fprintf(fpData,"%% Column 4: standard deviations\n"); 
               fprintf(fpData,"%% Set morePlots=1 for normalized residual error plot.\n");
            }
            else printf("INFO: cannot open file RSFA_CV_err.m.\n");
         }
         if (fpData != NULL)
         {
            fprintf(fpData,"disp('Set morePlots=1 for normalized residual error plot')\n");
            ssum = 0.0;
            fprintf(fpData, "morePlots = 0;\n");
            fprintf(fpData, "A = [\n");
            for (iD = 0; iD < nSamples; iD++)
            {
               fprintf(fpData, "  %e %e %e %e\n", eArray[iD], YLocal[iD],
                       sArray[iD], sigmas[iD]);
               ssum += PABS(sigmas[iD]);
            }
            fprintf(fpData, "];\n");
            fwritePlotCLF(fpData);
            fwritePlotFigure(fpData, 1);
            fprintf(fpData, "subplot(1,2,1)\n");
            if (psPlotTool_ == 1)
            {
               fprintf(fpData, "ymin = min(A(:,1));\n");
               fprintf(fpData, "ymax = max(A(:,1));\n");
               fprintf(fpData, "ywid = 0.1 * (ymax - ymin);\n");
               fprintf(fpData, "if (ywid < 1.0e-12)\n");
               fprintf(fpData, "   disp('range too small.')\n");
               fprintf(fpData, "   halt\n");
               fprintf(fpData, "end;\n");
               fprintf(fpData, "histplot(10, A(:,1), style=2);\n");
               fprintf(fpData, "a = gce();\n");
               fprintf(fpData, "a.children.fill_mode = \"on\";\n");
               fprintf(fpData, "a.children.thickness = 2;\n");
               fprintf(fpData, "a.children.foreground = 0;\n");
               fprintf(fpData, "a.children.background = 2;\n");
            }
            else
            {
               fprintf(fpData, "[nk, xk] = hist(A(:,1), 10);\n");
               fprintf(fpData, "bar(xk,nk/sum(nk))\n");
            }
            fwritePlotAxes(fpData);
            fwritePlotTitle(fpData, "CV Error Histogram");
            fwritePlotXLabel(fpData, "Error");
            fwritePlotYLabel(fpData, "Probabilities");
            if (psPlotTool_ != 1)
            {
               fprintf(fpData,"disp(['Error Mean  = ' num2str(mean(A(:,1)))])\n");
               fprintf(fpData,"disp(['Error stdev = ' num2str(std(A(:,1)))])\n");
            }
            fprintf(fpData, "subplot(1,2,2)\n");
            fprintf(fpData, "xmax = max(A(:,2));\n");
            fprintf(fpData, "xmin = min(A(:,2));\n");
            fprintf(fpData, "if (xmax == xmin) \n");
            fprintf(fpData, "   xmin = 0.9 * xmin; \n");
            fprintf(fpData, "   xmax = 1.1 * xmax; \n");
            fprintf(fpData, "end;\n");
            fprintf(fpData, "if (xmax == xmin) \n");
            fprintf(fpData, "   xmin = -0.1; \n");
            fprintf(fpData, "   xmax = 0.1; \n");
            fprintf(fpData, "end;\n");
            fprintf(fpData, "ymax = max(A(:,3)+A(:,4));\n");
            fprintf(fpData, "ymin = min(A(:,3)-A(:,4));\n");
            fprintf(fpData, "if (ymax == ymin) \n");
            fprintf(fpData, "   ymin = 0.9 * ymin; \n");
            fprintf(fpData, "   ymax = 1.1 * ymax; \n");
            fprintf(fpData, "end;\n");
            fprintf(fpData, "if (ymax == ymin) \n");
            fprintf(fpData, "   ymin = -0.1; \n");
            fprintf(fpData, "   ymax = 0.1; \n");
            fprintf(fpData, "end;\n");
            fprintf(fpData, "xmin = min(xmin, ymin);\n");
            fprintf(fpData, "xmax = max(xmax, ymax);\n");
            fprintf(fpData, "XX = xmin : xmax-xmin : xmax;\n");
            fprintf(fpData, "plot(A(:,2), A(:,3),'*','MarkerSize',15)\n");
            if (psPlotTool_ == 1)
                 fprintf(fpData, "set(gca(),\"auto_clear\",\"off\")\n");
            else fprintf(fpData, "hold on\n");
            fprintf(fpData, "plot(XX, XX)\n");
            if (psPlotTool_ == 1)
            {
               fprintf(fpData, "a = get(\"current_axes\");\n");
               fprintf(fpData, "a.data_bounds=[xmin,xmin;xmax,xmax];\n");
            }
            else
            {
               fprintf(fpData, "axis([xmin xmax xmin xmax])\n");
            }
            if (ssum > 0.0)
            {
               if (psPlotTool_ == 1)
               {
                  fprintf(fpData,"drawlater\n");
                  fprintf(fpData,"for ii = 1 : %d\n", nSamples);
                  fprintf(fpData,"   xx = [A(ii,2) A(ii,2)];\n");
                  fprintf(fpData,"   yy = [A(ii,3)-A(ii,4) A(ii,3)+A(ii,4)];\n");
                  fprintf(fpData,"   plot(xx, yy, 'r-', 'lineWidth', 1);\n");
                  fprintf(fpData,"   plot(xx(1), yy(1),'rv');\n");
                  fprintf(fpData,"   plot(xx(2), yy(2),'r^');\n");
                  fprintf(fpData,"end;\n");
                  fprintf(fpData,"drawnow\n");
                  fprintf(fpData, "set(gca(),\"auto_clear\",\"on\")\n");
               }
               else
               {
                  fprintf(fpData,"for ii = 1 : %d\n", nSamples);
                  fprintf(fpData,"   xx = [A(ii,2) A(ii,2)];\n");
                  fprintf(fpData,"   yy = [A(ii,3)-A(ii,4) A(ii,3)+A(ii,4)];\n");
                  fprintf(fpData,"   plot(xx, yy, 'r-', 'lineWidth', 1);\n");
                  fprintf(fpData,"   plot(xx(1), yy(1),'rv','markerSize',10);\n");
                  fprintf(fpData,"   plot(xx(2), yy(2),'^','markerSize',10);\n");
                  fprintf(fpData,"end;\n");
                  fprintf(fpData,"text(0.1,0.9,'red: +/- 1 sigma','sc',");
                  fprintf(fpData,"'fontSize',11,'fontweight','bold')\n");
               }
            }
            if (psPlotTool_ == 1)
                 fprintf(fpData, "set(gca(),\"auto_clear\",\"on\")\n");
            else fprintf(fpData, "hold off\n");
            fwritePlotAxes(fpData);
            fwritePlotTitle(fpData, "Actual vs predicted data with CV");
            fwritePlotXLabel(fpData, "actual data");
            fwritePlotYLabel(fpData, "predicted data");
            fprintf(fpData,"if morePlots == 1\n");
            fwritePlotFigure(fpData, 2);
            fprintf(fpData, "B = [\n");
            for (iD = 0; iD < nSamples; iD++)
            {
               if (YLocal[iD] == 0) fprintf(fpData, " %e 0 ",YLocal[iD]);
               else fprintf(fpData," %e %e ",YLocal[iD],eArray[iD]/YLocal[iD]);
               for (iD2 = 0; iD2 < nInputs; iD2++)
                  fprintf(fpData," %e ",XX[iD*nInputs+iD2]);
               fprintf(fpData,"\n");
            }
            fprintf(fpData, "];\n");
            if (psPlotTool_ == 1)
                 fprintf(fpData, "[AA,II] = gsort(B(:,1),'g','i');\n");
            else fprintf(fpData, "[AA,II] = sort(B(:,1));\n");
            fprintf(fpData, "BB = B(II,2);\n");
            fprintf(fpData, "plot(AA, BB, '*','markerSize',12);\n");
            fwritePlotAxes(fpData);
            fwritePlotTitle(fpData, "Normalized Resdual Analysis");
            fwritePlotXLabel(fpData, "Actual Data");
            fwritePlotYLabel(fpData, "Normalized Error");
            fprintf(fpData,"end;\n");
            fclose(fpData);
            if (psPlotTool_ == 1)
                 printf("CV error file is RSFA_CV_err.sci\n");
            else printf("CV error file is RSFA_CV_err.m\n");
         }
      }
      psRSExpertMode_ = rsState;
      delete fa;
      delete [] XX;
      delete [] YY;
      delete [] YT;
      delete [] WW;
      delete [] X2;
      delete [] Y2;
      delete [] S2;
      delete [] wgts;
      delete [] iArray;
      delete [] iArray2;
      delete [] eArray;
      delete [] sArray;
      delete [] sigmas;
   }
   delete [] YLocal;

   testFlag = 0;
   if (psConfig_ != NULL)
   {
      cString = psConfig_->getParameter("RSFA_test_datafile");
      if (cString != NULL)
      {
         sscanf(cString, "%s %s %s",winput1,winput2,dataFile);
         fpData = fopen(dataFile, "r");
         if (fpData != NULL)
         {
            testFlag = 1;
            printf("RSFA:: test data file = %s\n", dataFile);
            fclose(fpData);
            fpData = NULL;
         }
      }
      if (testFlag == 1)
      {
         cString = psConfig_->getParameter("RSFA_test_errorfile");
         if (cString != NULL)
         {
            sscanf(cString, "%s %s %s",winput1,winput2,errFile);
            fpErr = fopen(errFile, "w");
            if (fpErr != NULL)
            {
               printf("RSFA:: error file = %s\n", errFile);
               testFlag += 2;
               fclose(fpErr);
            }
         }
      }
   }
   if (testFlag == 1) retdata = validate(adata, dataFile, NULL);
   if (testFlag == 3) retdata = validate(adata, dataFile, errFile);
   printAsterisks(0);

   return retdata;
}

// ************************************************************************ 
// set parameters
// ------------------------------------------------------------------------
int RSFuncApproxAnalyzer::setParams(int argc, char **argv)
{
   int   idata;
   char  *request, *dataFile, *errFile;
   aData *adata;

   request = (char *) argv[0]; 
   if (!strcmp(request, "validate"))
   {
      if (argc != 4) printf("RSAnalysis WARNING: setParams.\n");
      adata    = (aData *) argv[1];
      dataFile = (char *) argv[2];
      errFile  = (char *) argv[3];
      validate(*adata, dataFile, errFile);
   }
   else if (!strcmp(request, "rstype"))
   {
      if (argc != 2) printf("RSAnalysis WARNING: setParams.\n");
      rsType_ = *(int *) argv[1];
      if (rsType_ < 0 || rsType_ > PSUADE_NUM_RS)
      {
         printf("RSAnalysis ERROR: INVALID rstype, set to REGR2.\n");
         rsType_ = PSUADE_RS_REGR2;
      }
   }
   else if (!strcmp(request, "usecv"))
   {
      useCV_ = 1;
   }
   else if (!strcmp(request, "numcvgroups"))
   {
      idata = *(int *) argv[1];
      if (idata > 1) numCVGroups_ = idata;
   }
   else
   {
      printf("RSAnalysis ERROR: setParams - not valid.\n");
      exit(1);
   }
   return 0;
}

// ************************************************************************ 
// validate
// ------------------------------------------------------------------------
double RSFuncApproxAnalyzer::validate(aData &adata, char *dataFile, 
                                      char *errFile)
{
   int        nInputs, nOutputs, outputID, wgtID, length, ii, iOne=1,fatype;
   int        nSamples, status, nTestOut, nTestIn, nPtsPerDim=64, nTestSam;
   double     *X, *Y, *lower, *upper, *XX, *YY, *YLocal, *wgts, *stdevs;
   double     sumErr1, sumErr2, maxErr, sumErr1s, sumErr2s, maxErrs;
   double     dataMax, dataMin, *YT, maxBase, maxBases, sdata, ddata, ssum;
   pData      pPtr, pInputs, pOutputs;
   PsuadeData *ioPtr;
   FuncApprox *fa;
   FILE       *fpErr;

   printf("RSAnalysis: validating against a test set ...\n");

   nInputs   = adata.nInputs_;
   nOutputs  = adata.nOutputs_;
   nSamples  = adata.nSamples_;
   lower     = adata.iLowerB_;
   upper     = adata.iUpperB_;
   X         = adata.sampleInputs_;
   Y         = adata.sampleOutputs_;
   outputID  = adata.outputID_;
   wgtID     = adata.regWgtID_;
   fatype    = adata.faType_;
   YLocal    = new double[nSamples];
   for (ii = 0; ii < nSamples; ii++) YLocal[ii] = Y[ii*nOutputs+outputID];

   ioPtr = new PsuadeData();
   status = ioPtr->readPsuadeFile(dataFile);
   ioPtr->getParameter("output_noutputs", pPtr);
   nTestOut = pPtr.intData_;
   ioPtr->getParameter("method_nsamples", pPtr);
   nTestSam = pPtr.intData_;
   ioPtr->getParameter("input_ninputs", pPtr);
   nTestIn = pPtr.intData_;
   if (nTestSam <= 0)
   {
      printf("RSAnalysis: file %s has no data.\n", dataFile);
      delete ioPtr;
      delete [] YLocal;
      return 1.0e12;
   }
   if (nTestIn != nInputs)
   {
      printf("RSAnalysis: test data has different number of inputs (%d)\n",
              nTestIn);
      printf("            than that of the sample (%d)\n", nInputs);
      delete ioPtr;
      delete [] YLocal;
      return 1.0e12;
   }
   ioPtr->getParameter("input_sample",pInputs);
   XX = pInputs.dbleArray_;
   ioPtr->getParameter("output_sample",pOutputs);
   YY = pOutputs.dbleArray_;

   fa = genFA(fatype, nInputs, iOne, nSamples);
   if (fa == NULL)
   {
      printf("RSAnalysis INFO: cannot create function approximator.\n");
      delete [] YLocal;
      delete ioPtr;
      return 1.0e12;
   }
   fa->setNPtsPerDim(nPtsPerDim);
   fa->setBounds(adata.iLowerB_, adata.iUpperB_);
   fa->setOutputLevel(adata.printLevel_);
   if (wgtID >= 0 && wgtID < nOutputs)
   {
      wgts = new double[nSamples];
      for (ii = 0; ii < nSamples; ii++) wgts[ii] = Y[ii*nOutputs+wgtID];
      fa->loadWeights(nSamples, wgts);
      delete [] wgts;
   }
   length = -999;
   status = fa->genNDGridData(X, YLocal, &length, NULL, NULL);

   if (errFile != NULL) 
   {
      fpErr = fopen(errFile, "w");
      if (fpErr != NULL)
      {
         if (psPlotTool_ == 1)
         {
            fprintf(fpErr, "// Surface Fitting Error Histogram \n");
            fprintf(fpErr, "// Interpolation errors on all points\n");
         }
         else
         {
            fprintf(fpErr, "%% Surface Fitting Error Histogram \n");
            fprintf(fpErr, "%% Interpolation errors on all points\n");
            fprintf(fpErr, "%% Column 1: predicted data \n");
            fprintf(fpErr, "%% Column 2: actual data \n");
            fprintf(fpErr, "%% Column 3: error = predicted - actual data \n");
            fprintf(fpErr, "%% Column 4: normalized error\n");
            fprintf(fpErr, "%% Column 5: predicted standard deviation\n");
         }
         fprintf(fpErr, "E = [\n");
      }
   }
   else fpErr = NULL;

   dataMax = -PSUADE_UNDEFINED;
   dataMin = PSUADE_UNDEFINED;
   for (ii = 0; ii < nTestSam; ii++)
   {
      ddata = PABS(YY[ii*nTestOut+outputID]);
      if (ddata > dataMax) dataMax = ddata;
      if (ddata < dataMin) dataMin = ddata;
   }
   YT = new double[nTestSam];
   stdevs = new double[nTestSam];
   fa->evaluatePointFuzzy(nTestSam, XX, YT, stdevs);
   sumErr1 = sumErr1s = sumErr2 = sumErr2s = maxErr = maxErrs = 0.0;
   for (ii = 0; ii < nTestSam; ii++)
   {
      ddata = YT[ii];
      if (fpErr != NULL) 
      {
         if (YY[ii*nTestOut+outputID] != 0)
            fprintf(fpErr, "%16.8e %16.8e %16.8e %16.8e %16.8e\n", ddata, 
              YY[ii*nTestOut+outputID],ddata-YY[ii*nTestOut+outputID],
              PABS((ddata-YY[ii*nTestOut+outputID])/YY[ii*nTestOut+outputID]),
              stdevs[ii]);
         else
            fprintf(fpErr, "%16.8e %16.8e %16.8e %16.8e %16.8e\n", ddata, 
              YY[ii*nTestOut+outputID],ddata-YY[ii*nTestOut+outputID],
              PABS((ddata-YY[ii*nTestOut+outputID])), stdevs[ii]);
      }
      ddata = ddata - YY[ii*nTestOut+outputID];
      if (YY[ii*nTestOut+outputID] != 0.0)
           sdata = ddata / PABS(YY[ii*nTestOut+outputID]);
      else sdata = ddata;
      sumErr1  += ddata;
      sumErr1s += sdata;
      sumErr2  += (ddata * ddata);
      sumErr2s += (sdata * sdata);
      ddata = PABS(ddata);
      sdata = PABS(sdata);
      if (ddata > maxErr ) 
      {
         maxErr = ddata;
         maxBase = PABS(YY[ii*nTestOut+outputID]);
      }
      if (sdata > maxErrs)
      {
         maxErrs  = sdata;
         maxBases = PABS(YY[ii*nTestOut+outputID]);
      }
   }
   delete [] YT;

   sumErr1  = sumErr1 / (double) nTestSam;
   sumErr1s = sumErr1s / (double) nTestSam;
   sumErr2  = sqrt(sumErr2 / (double) nTestSam);
   sumErr2s = sqrt(sumErr2s / (double) nTestSam);
   printf("RSA: Test data maximum/minimum      = %9.3e %9.3e\n",
          dataMax, dataMin);
   printf("RSA: Prediction errors = %11.3e (max unscaled, BASE=%9.3e)\n",
          maxErr, maxBase);
   printf("RSA: Prediction errors = %11.3e (max   scaled, BASE=%9.3e)\n",
          maxErrs, maxBases);
   printf("RSA: Prediction errors = %11.3e (rms unscaled)\n",sumErr2);
   printf("RSA: Prediction errors = %11.3e (rms   scaled)\n",sumErr2s);
   printf("RSA: Prediction errors = %11.3e (avg unscaled)\n",sumErr1);
   printf("RSA: Prediction errors = %11.3e (avg   scaled)\n",sumErr1s);
   ssum = 0.0;
   for (ii = 0; ii < nTestSam; ii++) ssum += stdevs[ii];
   if (ssum != 0)
      printf("RSA: Average std. dev. = %11.3e (sum of all points)\n",ssum/nSamples);

   if (fpErr != NULL) 
   {
#if 0
      fprintf(fpErr, "];\n");
      fprintf(fpErr, "Emax = max(E(:,3));\n");
      fprintf(fpErr, "Emin = min(E(:,3));\n");
      fprintf(fpErr, "Estp = 0.1 * (Emax - Emin);\n");
      fprintf(fpErr, "X = zeros(11,1);\n");
      fprintf(fpErr, "Y = zeros(10,1);\n");
      fprintf(fpErr, "for i = 1 : 11\n");
      fprintf(fpErr, "   X(i) = Emin + (i - 1) * Estp;\n");
      fprintf(fpErr, "end;\n");
      fprintf(fpErr, "for i = 1 : 10\n");
      fprintf(fpErr, "   inds = find(E(:,3) >= X(i));\n");
      fprintf(fpErr, "   vals = E(inds,3);\n");
      fprintf(fpErr, "   inds = find(vals < X(i+1));\n");
      fprintf(fpErr, "   Y(i) = length(inds);\n");
      fprintf(fpErr, "end;\n");
      fprintf(fpErr, "Y(10) = Y(10) + 1;\n");
      fprintf(fpErr, "bar(X(1:10),Y)\n");
      fprintf("Error file in %s.\n", errFile);
#else
      fprintf(fpErr, "];\n");
      fwritePlotCLF(fpErr);
      fwritePlotFigure(fpErr, 1);
      fprintf(fpErr, "subplot(1,2,1)\n");
      fprintf(fpErr, "plot(E(:,3),'x')\n");
      fwritePlotAxes(fpErr);
      fwritePlotTitle(fpErr, "Error Plot");
      fwritePlotXLabel(fpErr, "Sample Number");
      fwritePlotYLabel(fpErr, "Error");
      fprintf(fpErr, "subplot(1,2,2)\n");
      if (psPlotTool_ == 1)
      {
         fprintf(fpErr, "ymin = min(E(:,3));\n");
         fprintf(fpErr, "ymax = max(E(:,3));\n");
         fprintf(fpErr, "ywid = 0.1 * (ymax - ymin);\n");
         fprintf(fpErr, "if (ywid < 1.0e-12)\n");
         fprintf(fpErr, "   disp('range too small.')\n");
         fprintf(fpErr, "   halt\n");
         fprintf(fpErr, "end;\n");
         fprintf(fpErr, "histplot(10, E(:,3), style=2);\n");
         fprintf(fpErr, "a = gce();\n");
         fprintf(fpErr, "a.children.fill_mode = \"on\";\n");
         fprintf(fpErr, "a.children.thickness = 2;\n");
         fprintf(fpErr, "a.children.foreground = 0;\n");
         fprintf(fpErr, "a.children.background = 2;\n");
      }
      else
      {
         fprintf(fpErr, "[nk,xk]=hist(E(:,3),10);\n");
         fprintf(fpErr, "bar(xk,nk/%d,1.0)\n",nTestSam);
      }
      fwritePlotAxes(fpErr);
      fwritePlotTitle(fpErr, "Histogram of Errors");
      fwritePlotXLabel(fpErr, "Error");
      fwritePlotYLabel(fpErr, "Probabilities");
      fwritePlotFigure(fpErr, 2);
      fprintf(fpErr, "xmax = max(E(:,2));\n");
      fprintf(fpErr, "xmin = min(E(:,2));\n");
      fprintf(fpErr, "if (xmax == xmin) \n");
      fprintf(fpErr, "   xmin = 0.9 * xmin; \n");
      fprintf(fpErr, "   xmax = 1.1 * xmax; \n");
      fprintf(fpErr, "end;\n");
      fprintf(fpErr, "if (xmax == xmin) \n");
      fprintf(fpErr, "   xmin = -0.1; \n");
      fprintf(fpErr, "   xmax = 0.1; \n");
      fprintf(fpErr, "end;\n");
      fprintf(fpErr, "ymax = max(E(:,1)+E(:,5));\n");
      fprintf(fpErr, "ymin = min(E(:,1)-E(:,5));\n");
      fprintf(fpErr, "if (ymax == ymin) \n");
      fprintf(fpErr, "   ymin = 0.9 * ymin; \n");
      fprintf(fpErr, "   ymax = 1.1 * ymax; \n");
      fprintf(fpErr, "end;\n");
      fprintf(fpErr, "if (ymax == ymin) \n");
      fprintf(fpErr, "   ymin = -0.1; \n");
      fprintf(fpErr, "   ymax = 0.1; \n");
      fprintf(fpErr, "end;\n");
      fprintf(fpErr, "xmin = min(xmin, ymin);\n");
      fprintf(fpErr, "xmax = max(xmax, ymax);\n");
      fprintf(fpErr, "plot(E(:,2), E(:,1),'*','MarkerSize',12)\n");
      if (psPlotTool_ == 1)
           fprintf(fpErr, "set(gca(),\"auto_clear\",\"off\")\n");
      else fprintf(fpErr, "hold on\n");
      fprintf(fpErr, "XX = xmin : xmax-xmin : xmax;\n");
      fprintf(fpErr, "plot(XX, XX)\n");
      if (psPlotTool_ == 1)
      {
         fprintf(fpErr, "a = get(\"current_axes\");\n");
         fprintf(fpErr, "a.data_bounds=[xmin,xmin;xmax,xmax];\n");
      }
      else
      {
         fprintf(fpErr, "axis([xmin xmax xmin xmax])\n");
      }
      if (ssum > 0.0)
      {
         fprintf(fpErr,"for ii = 1 : %d\n", nTestSam);
         fprintf(fpErr,"   xx = [E(ii,2) E(ii,2)];\n");
         fprintf(fpErr,"   yy = [E(ii,1)-E(ii,5) E(ii,1)+E(ii,5)];\n");
         fprintf(fpErr,"   plot(xx, yy, 'r-', 'lineWidth', 3);\n");
         fprintf(fpErr,"   plot(xx(1), yy(1), 'rv', 'markerSize', 10);\n");
         fprintf(fpErr,"   plot(xx(2), yy(2), '^', 'markerSize', 10);\n");
         fprintf(fpErr,"end;\n");
      }
      if (psPlotTool_ == 1)
           fprintf(fpErr, "set(gca(),\"auto_clear\",\"on\")\n");
      else 
      {
         fprintf(fpErr, "hold off\n");
         fprintf(fpErr,"text(0.1,0.9,'red: +/- 1 sigma','sc',");
         fprintf(fpErr,"'fontSize',11,'fontweight','bold')\n");
      }
      fwritePlotAxes(fpErr);
      fwritePlotTitle(fpErr, "Interpolated versus actual data");
      fwritePlotXLabel(fpErr, "Actual data");
      fwritePlotYLabel(fpErr, "Interpolated data");
#endif
      fclose(fpErr);
      fpErr = NULL;
      printf("RSAnalysis: individual prediction errors can be found in %s.\n",
             errFile);
   }

   delete [] YLocal;
   delete [] stdevs;
   delete fa;
   delete ioPtr;
   return sumErr1;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
RSFuncApproxAnalyzer& RSFuncApproxAnalyzer::operator=(const RSFuncApproxAnalyzer &)
{
   printf("RSAnalysis operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

