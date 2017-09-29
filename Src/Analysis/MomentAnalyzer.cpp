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
// Functions for the class MomentAnalyzer  
// Reference: from any introductory statistics book
// AUTHOR : CHARLES TONG
// DATE   : 2004
// ************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include "MomentAnalyzer.h"
#include "Main/Psuade.h"
#include "Util/sysdef.h"
#include "Util/PsuadeUtil.h"
#include "DataIO/pData.h"
#include "Samplings/RSConstraints.h"

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
MomentAnalyzer::MomentAnalyzer() : Analyzer()
{
   setName("MOMENT");
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
MomentAnalyzer::~MomentAnalyzer()
{
}

// ************************************************************************
// perform analysis
// ------------------------------------------------------------------------
double MomentAnalyzer::analyze(aData &adata)
{
   int    nInputs, nOutputs, nSamples, outputID, nSubSamples, status, nValid;
   int    nGroups, groupID, whichOutput, ss, printLevel, nConstr, cnt, outID;
   int    errFlag, ii;
   double *gMeans, *gVariances, mean, variance, skewness, kurtosis, *X;
   double *Y2, *YY, *Y=NULL, variance2;
   char   winput[500], filename[500], pString[500];
   FILE   *fp=NULL;
   RSConstraints *constrPtr;
   PsuadeData    *ioPtr;
   pData         pPtr;

   printLevel = adata.printLevel_;
   nInputs  = adata.nInputs_;
   nOutputs = adata.nOutputs_;
   nSamples = adata.nSamples_;
   X        = adata.sampleInputs_;
   YY       = adata.sampleOutputs_;
   outputID = adata.outputID_;
   nSubSamples = adata.nSubSamples_;
   ioPtr       = adata.ioPtr_;
   if (adata.inputPDFs_ != NULL)
   {
      printf("MomentAnalyzer INFO: non-uniform probability distributions\n");
      printf("                     have been defined in the data file, but\n");
      printf("                     they will not be used in this analysis.\n");
   }
   if (ioPtr != NULL) ioPtr->getParameter("output_names", pPtr);

   if (nInputs <= 0 || nOutputs <= 0)
   {
      printf("MomentAnalyzer ERROR: invalid nInputs or nOutputs.\n");
      printf("   nInputs  = %d\n", nInputs);
      printf("   nOutputs = %d\n", nOutputs); 
      return PSUADE_UNDEFINED;
   } 
   if (nSamples <= 1)
   {
      printf("MomentAnalyzer ERROR: nSamples should be > 1.\n");
      printf("   nSamples = %d\n", nSamples);
      return PSUADE_UNDEFINED;
   } 
   if (nSubSamples <= 0)
   {
      printf("MomentAnalyzer ERROR: nSubSamples should be > 0.\n");
      printf("   nSubSamples = %d\n", nSubSamples);
      return PSUADE_UNDEFINED;
   } 
   if (nSamples/nSubSamples*nSubSamples != nSamples)
   {
      printf("MomentAnalyzer ERROR : nSamples != k*nSubSamples.\n");
      printf("   nSamples    = %d\n", nSamples);
      printf("   nSubSamples = %d\n", nSubSamples);
      return PSUADE_UNDEFINED;
   } 
   whichOutput = outputID;
   if (whichOutput < -1 || whichOutput >= nOutputs)
   {
      printf("MomentAnalyzer ERROR: invalid outputID (%d).\n",whichOutput+1);
      printf("   outputID = %d\n", outputID+1);
      return PSUADE_UNDEFINED;
   } 
   if (whichOutput < 0) whichOutput = 0;
   errFlag = 0;
   for (ss = 0; ss < nSamples; ss++)
   {
      for (ii = 0; ii < nOutputs; ii++)
      {
         if (YY[ss*nOutputs+ii] == PSUADE_UNDEFINED)
         {
            errFlag++;
            break;
         }
      }
   }
   if (errFlag == 0) Y = YY;
   else
   {
      printf("MomentAnalyzer INFO: Some outputs are undefined,\n");
      printf("               which are purged before processing.\n");
      Y = new double[nSamples*nOutputs];
      cnt = 0;
      for (ss = 0; ss < nSamples; ss++)
      {
         for (ii = 0; ii < nOutputs; ii++)
            if (YY[ss*nOutputs+ii] == PSUADE_UNDEFINED) break;
         if (ii == nOutputs)
         {
            for (ii = 0; ii < nOutputs; ii++)
               Y[cnt*nOutputs+ii] = YY[ss*nOutputs+ii];
            cnt++;
         } 
      }
      if (nSubSamples == nSamples) nSubSamples = cnt;
      printf("MomentAnalyzer INFO: original sample size = %d\n",nSamples);
      nSamples = cnt;
      printf("MomentAnalyzer INFO: purged   sample size = %d\n",nSamples);
   }
   
   nConstr = 0;
   constrPtr = NULL;
   if (ioPtr != NULL)
   {
      constrPtr = new RSConstraints();
      constrPtr->genConstraints(ioPtr);
      nConstr = constrPtr->genConstraints(ioPtr);
   }
   if (nConstr > 0)
   {
      if (psAnalysisInteractive_ == 1 || psAnaExpertMode_ == 1)
      {
         sprintf(pString,"Generate output distribution plot? (y or n) ");
         getString(pString, winput);
         if (winput[0] == 'y')
         {
            if (outputID < 0 || outputID >= nOutputs)
            {
               sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
               outID = getInt(1, nOutputs, pString);
               outID--;
            }
            else outID = outputID;
            sprintf(pString, "Type in a matlab/scilab file name (no extension): ");
            getString(pString, filename);
            cnt = strlen(filename);
            filename[cnt-1] = '.';
            if (psPlotTool_ == 1)
            {
               filename[cnt] = 's';
               filename[cnt+1] = 'c';
               filename[cnt+2] = 'i';
               filename[cnt+3] = '\0';
            }
            else
            {
               filename[cnt] = 'm';
               filename[cnt+1] = '\0';
            }
            fp = fopen(filename, "w");
            if (fp == NULL)
               printf("ERROR: cannot open file %s.\n", filename);
         }
      }
      Y2 = new double[nSamples];
      nValid = 0;
      if (fp != NULL) fprintf(fp, "Y = [\n");
      for (ss = 0; ss < nSamples; ss++)
      {
         constrPtr->evaluate(&X[ss*nInputs],Y[ss*nOutputs+outputID],status);
         if (status == 1)
         {
            Y2[nValid++] = Y[ss*nOutputs+outputID];
            if (fp != NULL) fprintf(fp, "  %24.16e \n",Y2[nValid-1]);
         }
      }
      if (nValid == 0)
      {
         printf("MomemtAnalyzer INFO: no valid data points.\n");
         delete [] Y2;
         delete constrPtr;
         if (fp != NULL) fclose(fp);
         if (errFlag != 0 && Y != NULL) delete [] Y;
         return 0.0;
      }
      if (fp != NULL)
      {
         fprintf(fp, "];\n");
         if (psPlotTool_ == 1)
         {
            fwritePlotCLF(fp);
            fprintf(fp, "ymin = min(Y);\n");
            fprintf(fp, "ymax = max(Y);\n");
            fprintf(fp, "ywid = 0.1 * (ymax - ymin);\n");
            fprintf(fp, "if (ywid < 1.0e-12)\n");
            fprintf(fp, "   disp('range too small.')\n");
            fprintf(fp, "   halt\n");
            fprintf(fp, "end;\n");
            fprintf(fp, "histplot(10, Y/ywid, style=2);\n");
            fprintf(fp, "a = gce();\n");
            fprintf(fp, "a.children.fill_mode = \"on\";\n");
            fprintf(fp, "a.children.thickness = 2;\n");
            fprintf(fp, "a.children.foreground = 0;\n");
            fprintf(fp, "a.children.background = 2;\n");
            fwritePlotAxes(fp);
            fwritePlotTitle(fp, "Output Distribution");
            if (ioPtr != NULL)
                 fwritePlotXLabel(fp, pPtr.strArray_[outID]);
            else fwritePlotXLabel(fp, "Output Value");
            fwritePlotYLabel(fp, "Probabilities");
         }
         else
         {
            fwritePlotCLF(fp);
            fprintf(fp, "[nk,xk]=hist(Y(:,1),10);\n");
            fprintf(fp, "bar(xk,nk/%d,1.0)\n",nValid);
            fwritePlotAxes(fp);
            fwritePlotTitle(fp, "Output Distribution");
            if (ioPtr != NULL)
                 fwritePlotXLabel(fp, pPtr.strArray_[outID]);
            else fwritePlotXLabel(fp, "Output Value");
            fwritePlotYLabel(fp, "Probabilities");
         }
         printf("Output distribution plot is now in %s.\n",filename);
         fclose(fp);
      }
      computeMeanVariance(nInputs,1,nValid,Y2,&mean,&variance,0);
      if (printLevel >= 0)
      {
         printAsterisks(0);
         printf("*       Sample mean          = %12.4e\n",mean);
         printf("*       Sample std dev       = %12.4e\n",sqrt(variance));
         printf("*       Based on nSamples    = %d\n",nValid);
         printf("*       Output number        = %d\n", outID+1);
         printAsterisks(0);
      }
      delete [] Y2;
      delete constrPtr;
      if (errFlag != 0 && Y != NULL) delete [] Y;
      return sqrt(variance/nValid);
   }

   nGroups    = nSamples / nSubSamples;
   gMeans     = new double[nGroups];
   gVariances = new double[nGroups];

   printf("\n");
   if (printLevel >= 0)
   {
      printAsterisks(0);
      printf("*             Standard error of mean calculation\n");
      printEquals(0);
      printf("* nSamples = %10d\n",nSamples);
      printf("* nGroups  = %10d\n",nGroups);
      printf("* nInputs  = %10d\n",nInputs);
      printf("* nOutputs = %10d\n",nOutputs);
      printDashes(0);
   }
   if (outputID == -1)
   {
      for (ss = 0; ss < nOutputs; ss++)
      {
         if (printLevel >= 0) printf("* outputID = %10d\n", ss+1);
         computeMeanVariance(nInputs,nOutputs,nSubSamples, Y,
                             gMeans,gVariances,ss);
         if (printLevel >= 0)
         {
            printf("*       Sample mean          = %12.4e\n",gMeans[0]);
            printf("*       Sample std dev       = %12.4e\n",sqrt(gVariances[0]));
         }
         if (gVariances[0] > 0)
         {
            computeSkewnessKurtosis(nInputs,nOutputs,nSubSamples, Y, 
                         sqrt(gVariances[0]), &skewness,&kurtosis,ss);
            if (printLevel >= 0)
            {
               printf("*       Sample skewness      = %12.4e\n", skewness);
               printf("*       Sample kurtosis      = %12.4e\n", kurtosis);
            }
         }
         else
         {
            if (printLevel >= 0)
            {
               printf("*       Std dev = 0, skeweness and kurtosis not computed.\n");
            }
            skewness = kurtosis = 0.0;
         }
         if (ss < nOutputs-1 && printLevel >= 0) printEquals(0);
      }
      variance = gVariances[0];
   }
   else
   {
      if (printLevel >= 0) printf("* outputID = %10d\n", outputID+1);
      computeMeanVariance(nInputs,nOutputs,nSubSamples,Y,gMeans,gVariances,outputID);
      if (printLevel >= 0)
      {
         printf("*       Sample mean          = %12.4e\n",gMeans[0]);
         printf("*       Sample std dev       = %12.4e\n",sqrt(gVariances[0]));
      }
      if (gVariances[0] > 0)
      {
         computeSkewnessKurtosis(nInputs,nOutputs,nSubSamples, Y, 
                      sqrt(gVariances[0]), &skewness,&kurtosis,outputID);
         if (printLevel >= 0)
         {
            printf("*       Sample skewness      = %12.4e\n", skewness);
            printf("*       Sample kurtosis      = %12.4e\n", kurtosis);
         }
      }
      else
      {
         if (printLevel >= 0)
         {
            printf("*       Std dev = 0, skeweness and kurtosis not computed.\n");
         }
         skewness = kurtosis = 0.0;
      }
      variance = gVariances[0];
   }
   if (printLevel >= 0)
   {
      printDashes(0);
      printAsterisks(0);
   }

#ifdef HAVE_PYTHON
   PyObject *temp, *Ylist;

   PyDict_SetItemString( AnalysisDataDict, "nSamples",
			 temp=PyInt_FromLong(nSamples) ); Py_DECREF(temp);
   PyDict_SetItemString( AnalysisDataDict, "nGroups",
			 temp=PyInt_FromLong(nGroups) ); Py_DECREF(temp);
   PyDict_SetItemString( AnalysisDataDict, "nInputs",
			 temp=PyInt_FromLong(nInputs) ); Py_DECREF(temp);
   PyDict_SetItemString( AnalysisDataDict, "nOutputs",
			 temp=PyInt_FromLong(nOutputs) ); Py_DECREF(temp);
   PyDict_SetItemString( AnalysisDataDict, "outputID",
			 temp=PyInt_FromLong(whichOutput) ); Py_DECREF(temp);

   PyDict_SetItemString( AnalysisDataDict, "mean",
			 temp=PyFloat_FromDouble(gMeans[0]) ); Py_DECREF(temp);
   PyDict_SetItemString( AnalysisDataDict, "std dev",
			 temp=PyFloat_FromDouble(sqrt(gVariances[0])) );
   Py_DECREF(temp);
   PyDict_SetItemString( AnalysisDataDict, "skewness",
			 temp=PyFloat_FromDouble(skewness) ); Py_DECREF(temp);
   PyDict_SetItemString( AnalysisDataDict, "kurtosis",
			 temp=PyFloat_FromDouble(kurtosis) ); Py_DECREF(temp);

   Ylist = PyList_New(0);
   for (ss = 0; ss < nSamples; ss++)
      //if (YG[ii] != PSUADE_UNDEFINED) {
         PyList_Append( Ylist,
			temp=PyFloat_FromDouble(Y[nOutputs*ss+outputID]) );
	 Py_DECREF(temp);
      //}
   PyDict_SetItemString( AnalysisDataDict, "Y",	Ylist ); Py_DECREF(Ylist);
#endif

   if (psAnalysisInteractive_ ==1 || psAnaExpertMode_ == 1)
   {
      sprintf(pString,"Generate output distribution plot? (y or n) ");
      getString(pString, winput);
      if (winput[0] == 'y')
      {
         if (outputID < 0 || outputID >= nOutputs)
         {
            sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
            outID = getInt(1, nOutputs, pString);
            outID--;
         }
         else outID = outputID;
         sprintf(pString, "Type in a matlab/scilab file name : ");
         getString(pString, filename);
         cnt = strlen(filename);
         filename[cnt-1] = '.';
         if (psPlotTool_ == 1)
         {
            filename[cnt] = 's';
            filename[cnt+1] = 'c';
            filename[cnt+2] = 'i';
            filename[cnt+3] = '\0';
         }
         else
         {
            filename[cnt] = 'm';
            filename[cnt+1] = '\0';
         }
         fp = fopen(filename, "w");
         if (fp != NULL)
         {
            fprintf(fp, "Y = [\n");
            for (ss = 0; ss < nSamples; ss++)
               fprintf(fp, "  %24.16e\n", Y[nOutputs*ss+outID]);  
            fprintf(fp, "];\n");
            if (psPlotTool_ == 1)
            {
               fwritePlotCLF(fp);
               fprintf(fp, "ymin = min(Y);\n");
               fprintf(fp, "ymax = max(Y);\n");
               fprintf(fp, "ywid = 0.1 * (ymax - ymin);\n");
               fprintf(fp, "if (ywid < 1.0e-12)\n");
               fprintf(fp, "   disp('range too small.')\n");
               fprintf(fp, "   halt\n");
               fprintf(fp, "end;\n");
               fprintf(fp, "histplot(10, Y/ywid, style=2);\n");
               fprintf(fp, "a = gce();\n");
               fprintf(fp, "a.children.fill_mode = \"on\";\n");
               fprintf(fp, "a.children.thickness = 2;\n");
               fprintf(fp, "a.children.foreground = 0;\n");
               fprintf(fp, "a.children.background = 2;\n");
               fwritePlotAxes(fp);
               fwritePlotTitle(fp, "Output Distribution");
               if (ioPtr != NULL)
                    fwritePlotXLabel(fp, pPtr.strArray_[outID]);
               else fwritePlotXLabel(fp, "Output Values");
               fwritePlotYLabel(fp, "Probabilities");
            }
            else
            {
               fwritePlotCLF(fp);
               fprintf(fp, "twoPlots = 0;\n");
               fprintf(fp, "if (twoPlots == 1)\n");
               fprintf(fp, "subplot(1,2,1)\n");
               fprintf(fp, "end;\n");
               if (nSamples > 500) fprintf(fp, "[nk,xk]=hist(Y(:,1),20);\n");
               else                fprintf(fp, "[nk,xk]=hist(Y(:,1),10);\n");
               fprintf(fp, "bar(xk,nk/%d,1.0)\n",nSamples);
               fwritePlotAxes(fp);
               fwritePlotTitle(fp, "Probability Distribution");
               if (ioPtr != NULL)
                    fwritePlotXLabel(fp, pPtr.strArray_[outID]);
               else fwritePlotXLabel(fp, "Output Value");
               fwritePlotYLabel(fp, "Probabilities");
               fprintf(fp, "nk = nk / %d;\n", nSamples);
               fprintf(fp, "nk = nk / %d;\n", nSamples);
               if (nSamples > 500) fprintf(fp, "for ii = 2 : 20\n");
               else                fprintf(fp, "for ii = 2 : 10\n");
               fprintf(fp, "   nk(ii) = nk(ii) + nk(ii-1);\n");
               fprintf(fp, "end;\n");
               fprintf(fp, "Yk = sort(Y(:,1));\n");
               fprintf(fp, "Xk = 1 : %d;\n", nSamples);
               fprintf(fp, "Xk = Xk / %d;\n", nSamples);
               fprintf(fp, "if (twoPlots == 1)\n");
               fprintf(fp, "subplot(1,2,2)\n");
               fprintf(fp, "plot(Yk, Xk, 'lineWidth',3)\n");
               fwritePlotAxes(fp);
               fwritePlotTitle(fp, "Cumulative Distribution");
               if (ioPtr != NULL)
                    fwritePlotXLabel(fp, pPtr.strArray_[outID]);
               else fwritePlotXLabel(fp, "Output Values");
               fwritePlotYLabel(fp, "Probabilities");
               fprintf(fp, "end;\n");
            }
            printf("Output distribution plot is now in %s.\n",filename);
            fclose(fp);
         }
         else printf("ERROR: cannot open file %s.\n", filename);
      }
   }

   if (outputID >= 0)
   {
      if (nGroups == 1)
      {
         printAsterisks(0);
         variance2 = gVariances[0];
         delete [] gMeans;
         delete [] gVariances;
         if (printLevel >= 0)
         {
            printf("*       std error of mean       = %16.8e\n", 
                   sqrt(variance2/(double) nSamples));
            printAsterisks(0);
         }
         if (errFlag != 0 && Y != NULL) delete [] Y;
         return sqrt(variance2/(double) nSamples);
      }
      for (groupID = 0; groupID < nGroups; groupID++)
      {
         computeMeanVariance(nInputs,nOutputs,nSubSamples,
                &Y[nSubSamples*groupID*nOutputs],
                &gMeans[groupID],&gVariances[groupID],whichOutput);
      }
      if (printLevel >= 0) printf("Output %d\n", whichOutput+1);

      mean = 0.0;
      for (groupID = 0; groupID < nGroups; groupID++)
         mean += gVariances[groupID];
      mean = mean / (double) nGroups;
      variance2 = 0.0;
      for (groupID = 0; groupID < nGroups; groupID++)
         variance2 += ((gVariances[groupID] - mean) *
                      (gVariances[groupID] - mean));
      variance2 = variance2 / (double) nGroups;
      if (printLevel >= 0)
      {
         printf("*       Mean of group variances = %16.8e\n", mean);
         printf("*       Std error of variance   = %16.8e\n", variance2);
      }
      mean = 0.0;
      for (groupID = 0; groupID < nGroups; groupID++)
         mean += gMeans[groupID];
      mean = mean / (double) nGroups;
      variance2 = 0.0;
      for (groupID = 0; groupID < nGroups; groupID++)
         variance2 += ((gMeans[groupID] - mean) * (gMeans[groupID] - mean));
      variance2 = variance2 / (double) nGroups;
      if (printLevel >= 0)
      {
         printf("*       Mean of group means     = %16.8e\n", mean);
         printf("*       Std error of mean       = %16.8e\n", variance2);
         printAsterisks(0);
      }
   }

   if (printLevel > 1 && outputID >= 0)
      analyzeMore(nInputs, nOutputs, nSamples, nSubSamples, Y, outputID);

   delete [] gMeans;
   delete [] gVariances;
   if (errFlag != 0 && Y != NULL) delete [] Y;
   return sqrt(variance/nValid);
}

// *************************************************************************
// Compute agglomerated mean and variance
// -------------------------------------------------------------------------
int MomentAnalyzer::computeMeanVariance(int nInputs, int nOutputs, 
                        int nSamples, double *Y, double *mean, 
                        double *variance, int outputID)
{
   int    ss;
   double meanTmp, varTmp;

   meanTmp = 0.0;
   for (ss = 0; ss < nSamples; ss++) meanTmp += Y[nOutputs*ss+outputID];
   meanTmp /= (double) nSamples;
   varTmp = 0.0;
   for (ss = 0; ss < nSamples; ss++)
   {
      varTmp += ((Y[nOutputs*ss+outputID] - meanTmp) *
                 (Y[nOutputs*ss+outputID] - meanTmp));
   }
   varTmp /= (double) (nSamples - 1);
   (*mean)     = meanTmp;
   (*variance) = varTmp;
   return 0;
}

// *************************************************************************
// Compute skewness and kurtosis 
// kurtois of a standard normal distribution is 3.
// kurtois-3 measures excess kurtosis
// -------------------------------------------------------------------------
int MomentAnalyzer::computeSkewnessKurtosis(int nInputs, int nOutputs, 
                        int nSamples, double *Y, double stdev,
                        double *skewness, double *kurtosis, int outputID)
{
   int    ss;
   double meanTmp, skewTmp, kurTmp;

   meanTmp = 0.0;
   for (ss = 0; ss < nSamples; ss++) meanTmp += Y[nOutputs*ss+outputID];
   meanTmp /= (double) nSamples;
   skewTmp = kurTmp = 0.0;
   for (ss = 0; ss < nSamples; ss++)
   {
      skewTmp += ((Y[nOutputs*ss+outputID] - meanTmp) *
                  (Y[nOutputs*ss+outputID] - meanTmp) *
                  (Y[nOutputs*ss+outputID] - meanTmp));
      kurTmp += ((Y[nOutputs*ss+outputID] - meanTmp) *
                 (Y[nOutputs*ss+outputID] - meanTmp) *
                 (Y[nOutputs*ss+outputID] - meanTmp) *
                 (Y[nOutputs*ss+outputID] - meanTmp));
   }
   skewTmp /= (double) (nSamples - 1);
   kurTmp  /= (double) (nSamples - 1);
   if (stdev >= 0.0) (*skewness) = skewTmp / (stdev * stdev * stdev);
   else              (*skewness) = 0.0;
   if (stdev >= 0.0) (*kurtosis) = kurTmp / (stdev * stdev * stdev * stdev);
   else              (*kurtosis) = 0.0;
   return 0;
}

// ************************************************************************
// perform more analysis
// ------------------------------------------------------------------------
int MomentAnalyzer::analyzeMore(int nInputs, int nOutputs, int nSamples,
                              int nSubSamples, double *Y, int outputID)
{
   int    nGroups, groupID, myNSubSamples, nLevels=10, level, multiple;
   double *gMeans, *gVariances, mean, variance;

   for (level = 1; level <= nLevels; level++)
   {
      multiple = 1 << level;
      nGroups = nSamples / nSubSamples;
      if ((nGroups % multiple) != 0) break;
      myNSubSamples = nSubSamples * multiple;
      nGroups    = nSamples / myNSubSamples;
      gMeans     = new double[nGroups];
      gVariances = new double[nGroups];


      for (groupID = 0; groupID < nGroups; groupID++)
      {
         computeMeanVariance(nInputs,nOutputs,myNSubSamples,
                &Y[myNSubSamples*groupID*nOutputs],
                &gMeans[groupID],&gVariances[groupID],outputID);
      }
      mean = 0.0;
      for (groupID = 0; groupID < nGroups; groupID++)
         mean += gMeans[groupID];
      mean = mean / (double) nGroups;
      variance = 0.0;
      for (groupID = 0; groupID < nGroups; groupID++)
         variance += ((gMeans[groupID] - mean) * (gMeans[groupID] - mean));
      variance = variance / (double) nGroups;
      printEquals(0);
      printf("*       Agglomerate %d groups into 1 \n", multiple);
      printf("*       Mean of group means     = %16.8e\n", mean);
      printf("*       Std error of mean       = %16.8e\n", variance);
      printEquals(0);
      delete [] gMeans;
      delete [] gVariances;
   }
   return 0;
}

