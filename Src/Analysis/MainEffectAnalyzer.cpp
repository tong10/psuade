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
// functions for the class MainEffectAnalyzer  
// AUTHOR : CHARLES TONG
// DATE   : 2003
// ************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "MainEffectAnalyzer.h"
#include "Util/sysdef.h"
#include "Util/PsuadeUtil.h"
#include "Main/Psuade.h"
#include "Util/Matrix.h"
#include "DataIO/pData.h"
#include "DataIO/PsuadeData.h"
#include "Samplings/RSConstraints.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
MainEffectAnalyzer::MainEffectAnalyzer() : Analyzer()
{
   setName("MainEffect");
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
MainEffectAnalyzer::~MainEffectAnalyzer()
{
}

// ************************************************************************
// perform VCE (McKay) analysis (nSubSamples : number of LHS points)
// ------------------------------------------------------------------------
double MainEffectAnalyzer::analyze(aData &adata)
{
   int    nInputs, nOutputs, nSamples, outputID, nSubSamples, status;
   int    ii, jj, ir, ss, nReplications, ncount, index;
   int    whichOutput, nSubs, iZero=0, printLevel;
   double *txArray, *tyArray, *vce, *meanVCEVar, *varVCEMean;
   double aMean, aVariance, *varVCEVar, totalVCE;
   double *sampleInputs, *sampleOutputs, *X, *Y, ddata;
   double *XX, *YY, **bsVCEs, *iLowerB, *iUpperB;
   char   *cString, winput1[500], winput2[500], meFileName[500];
   char   pString[500];
   FILE   *fp=NULL, *fp1=NULL;
   PsuadeData    *ioPtr;
   RSConstraints *constrPtr=NULL;

   nInputs       = adata.nInputs_;
   nOutputs      = adata.nOutputs_;
   nSamples      = adata.nSamples_;
   outputID      = adata.outputID_;
   sampleInputs  = adata.sampleInputs_;
   sampleOutputs = adata.sampleOutputs_;
   nSubSamples   = adata.nSubSamples_;
   iLowerB       = adata.iLowerB_;
   iUpperB       = adata.iUpperB_;
   printLevel    = adata.printLevel_;
   whichOutput   = outputID;
   ioPtr         = adata.ioPtr_;
   if (adata.inputPDFs_ != NULL)
   {
      printf("MainEffect INFO: some inputs have non-uniform PDFs, but\n");
      printf("                 they will not be relevant in this analysis\n");
      printf("                 (since the sample should have been generated\n");
      printf("                 with the desired distributions.)\n");
   }

   if (nInputs <= 0 || nOutputs <= 0 || nSamples <= 0)
   {
      printf("MainEffectAnalyzer ERROR: invalid arguments.\n");
      printf("    nInputs  = %d\n", nInputs);
      printf("    nOutputs = %d\n", nOutputs);
      printf("    nSamples = %d\n", nSamples);
      exit(1);
   } 
   if (nSamples/nSubSamples*nSubSamples != nSamples)
   {
      printf("MainEffectAnalyzer ERROR: nSamples != k*nLevels.\n");
      printf("    nSamples = %d\n", nSamples);
      printf("    nLevels  = %d\n", nSubSamples);
      exit(1);
   } 
   if (whichOutput >= nOutputs || whichOutput < 0)
   {
      printf("MainEffectAnalyzer ERROR: invalid outputID.\n");
      printf("    nOutputs = %d\n", nOutputs);
      printf("    outputID = %d\n", whichOutput+1);
      exit(1);
   }
   if (ioPtr == NULL)
   {
      printf("MainEffectAnalyzer ERROR: no data (PsuadeData).\n");
      return PSUADE_UNDEFINED;
   }

   if (ioPtr != NULL)
   {
      constrPtr = new RSConstraints();
      constrPtr->genConstraints(ioPtr);
   }
   
   X = sampleInputs;
   Y = new double[nSamples];
   ncount = 0;
   for (ss = 0; ss < nSamples; ss++)
   {
      Y[ss] = sampleOutputs[nOutputs*ss+whichOutput];
      ddata = constrPtr->evaluate(&(X[ss*nInputs]), Y[ss], status);
      if (status == 0) Y[ss] = PSUADE_UNDEFINED; 
      else             ncount++;
   }
   if (ncount == 0)
   {
      printf("MainEffectAnalyzer ERROR: no valid sample point.\n");
      printf("    nSamples before filtering = %d\n", nSamples);
      printf("    nSamples after  filtering = %d\n", ncount);
      printf("    INFO: check your data file for undefined's (1e35)\n");
      return 1.0;
   } 
   if (ncount != nSamples)
   {
      printf("MainEffectAnalyzer: nSamples before filtering = %d\n", nSamples);
      printf("MainEffectAnalyzer: nSamples after filtering  = %d\n", ncount);
   }

   computeMeanVariance(nInputs,1,nSamples,Y,&aMean,&aVariance,0);
   if (printLevel > 0)
   {
      printAsterisks(0);
      printf("* Main Effect Analysis\n");
      printDashes(0);
      printf("* total number of samples = %10d                     **\n",
             nSamples);
      printf("* number of Inputs        = %10d                     **\n",nInputs);
      printDashes(0);
      printf("Output %d\n", whichOutput+1);
      printf("=====> MainEffectAnalyzer: mean               = %12.4e\n",aMean);
      printf("=====> MainEffectAnalyzer: standard deviation = %12.4e\n",
             sqrt(aVariance));
   }

   nSubs = nSubSamples;
   nReplications = nSamples / nSubs;
#if 0
   if (nReplications <= 1)
   {
      printf("MainEffectAnalyzer INFO: go no further since nReps = 1.\n");
      return 1.0;
   }
#endif
   vce           = new double[nInputs];
   meanVCEVar    = new double[nInputs];
   varVCEMean    = new double[nInputs];
   varVCEVar     = new double[nInputs];
   txArray       = new double[nSamples];
   tyArray       = new double[nSamples];

   for (ss = 0; ss < nSamples; ss++) tyArray[ss] = Y[ss];
   for (ii = 0; ii < nInputs; ii++)
   {
      for (ss = 0; ss < nSamples; ss++)
         txArray[ss] = X[nInputs*ss+ii];
      sortDbleList2(nSamples, txArray, tyArray);
      nReplications = 1;
      for (ss = 1; ss < nSamples; ss++)
      {
         if (txArray[ss] == txArray[0]) nReplications++;
         else                           break;
      }
      if (nReplications <= 1)
      {
         printf("analyzer INFO: nReps = 1 for input %d.\n",ii+1);
         printf("         Are you using replicated Latin hypercube?\n");
         printf("         If so, you need to use > 1 replications.\n");
         printf("Since replicated Latin hypercube is not used, a crude\n");
         printf("main effect analysis will be done below.\n");
         computeVCECrude(nInputs, nSamples, X, Y, iLowerB, iUpperB, 
                         aVariance, vce);
         delete [] vce;
         delete [] meanVCEVar;
         delete [] varVCEMean;
         delete [] varVCEVar;
         delete [] Y;
         delete [] txArray;
         delete [] tyArray;
         return 1.0;
      }
   }

   fp = NULL;
   if (psAnalysisInteractive_ == 1)
   {
      sprintf(pString,"Create main effect scatter plot ? (y or n) ");
      getString(pString, winput1);
      if (winput1[0] == 'y')
      {
         if (psPlotTool_ == 1)
            sprintf(pString,"Enter scatter plot file name (ends with .sci): ");
         else
            sprintf(pString,"Enter scatter plot file name (ends with .m): ");
         getString(pString, meFileName);
         meFileName[strlen(meFileName)-1] = '\0';
         fp = fopen(meFileName, "w");
         if (fp != NULL)
         {
            printf("MainEffectAnalyzer: main effect file = %s\n",
                   meFileName);
         }
      }
   }
   else
   {
      if (psConfig_ != NULL)
      {
         cString = psConfig_->getParameter("ME_matlab_file");
         if (cString != NULL)
         {
            sscanf(cString, "%s %s %s",winput1,winput2,meFileName);
            fp = fopen(meFileName, "w");
            if (fp != NULL)
            {
               printf("MainEffectAnalyzer: main effect file = %s\n",
                      meFileName);
            }
         }
      }
   }

   if (fp != NULL)
   {
      if (psPlotTool_ == 1)
      {
         fprintf(fp,"// ********************************************** \n");
         fprintf(fp,"// ********************************************** \n");
         fprintf(fp,"// *  Importance Analysis                      ** \n");
         fprintf(fp,"// *-------------------------------------------** \n");
         fprintf(fp,"// file for main effect plots \n");
      }
      else
      {
         fprintf(fp,"%% ********************************************** \n");
         fprintf(fp,"%% ********************************************** \n");
         fprintf(fp,"%% *  Importance Analysis                      ** \n");
         fprintf(fp,"%% *-------------------------------------------** \n");
         fprintf(fp,"%% file for main effect plots \n");
      }
   }

   if (fp != NULL)
   {
      for (ii = 0; ii < nInputs; ii++)
      {
         for (ss = 0; ss < nSamples; ss++)
            txArray[ss] = X[nInputs*ss+ii];
         sortDbleList2(nSamples, txArray, tyArray);
         fprintf(fp, "X%d = [ \n", ii);
         fprintf(fp,"%24.16e\n",txArray[0]);
         for (ss = 1; ss < nSamples; ss++)
            if (txArray[ss] != txArray[ss-1])
               fprintf(fp,"%24.16e\n",txArray[ss]);
         fprintf(fp, "];\n");
      }
   }

   status = computeVCE(nInputs, nSamples, nSubSamples, X, Y, whichOutput, fp, 
                       varVCEMean, meanVCEVar, varVCEVar, vce);
   if (status == -1)
   {
      computeVCECrude(nInputs, nSamples, X, Y, iLowerB, iUpperB, 
                      aVariance, vce);
      delete [] Y;
      delete [] vce;
      delete [] varVCEMean;
      delete [] meanVCEVar;
      delete [] varVCEVar;
      delete [] txArray;
      delete [] tyArray;
      if (constrPtr != NULL) delete constrPtr;
      return 0.0;
   }


   if (printLevel >= 0)
   {
      printAsterisks(0);
      printf("            McKay's correlation ratio\n");
      printDashes(0);
      printf(" - Turn on higher print levels to display more information\n");
      printf(" - Turn on analysis expert mode for more analysis\n");
      printEquals(0);
      totalVCE = 0.0;
      for (ii = 0; ii < nInputs; ii++)
      {
         totalVCE += vce[ii] / aVariance;
         printf("    INPUT %2d = %9.2e (raw = %9.2e)\n", ii+1, 
                vce[ii]/aVariance, vce[ii]);
      }
      printf("    Total VCE = %9.2e\n", totalVCE);
   }

#if 1
   if (printLevel > 2)
   {
      printAsterisks(0);
      printf("     McKay's biased correlation ratio (stdVCEMean)\n");
      printDashes(0);
      totalVCE = 0.0;
      for (ii = 0; ii < nInputs; ii++)
      {
         printf("    Input %2d = %9.2e >>? %9.2e\n", ii+1,
                varVCEMean[ii]/aVariance,1.0/(double)(nSubs*nSubs));
         totalVCE += varVCEMean[ii] / aVariance;
      }
      printf("    Total VCE = %9.2e\n", totalVCE);
   }
#endif

   if (printLevel > 2)
   {
      printAsterisks(0);
      printf("           Strength of interaction (varVCEVar)\n");
      printDashes(0);
      for (ii = 0; ii < nInputs; ii++)
         printf("    Input %2d = %9.2e (meanVCEVar = %9.2e)\n", ii+1, 
                varVCEVar[ii], meanVCEVar[ii]);
      printAsterisks(0);
      printf("            Hora and Iman sensitivity index\n");
      printf("   (may not be valid in the presence of constraints)\n");
      printDashes(0);
      totalVCE = 0.0;
      for (ii = 0; ii < nInputs; ii++) 
         totalVCE += sqrt(aVariance-meanVCEVar[ii]);
      for (ii = 0; ii < nInputs; ii++)
         printf("    Input %2d = %9.2e\n", ii+1, 
                sqrt(aVariance-meanVCEVar[ii])/totalVCE);
      printAsterisks(0);
   }

#ifdef HAVE_PYTHON
   PyObject *temp, *Xlist, *XIlist, *Ylist, *VCElist, *vVCEMlist;

   PyDict_SetItemString( AnalysisDataDict, "nInputs",
			 temp=PyInt_FromLong(nInputs) ); Py_DECREF(temp);
   PyDict_SetItemString( AnalysisDataDict, "outputID",
			 temp=PyInt_FromLong(whichOutput) ); Py_DECREF(temp);
   PyDict_SetItemString( AnalysisDataDict, "nSamples",
			 temp=PyInt_FromLong(nSamples) ); Py_DECREF(temp);
   PyDict_SetItemString( AnalysisDataDict, "nSubs",
			 temp=PyInt_FromLong(nSubs) ); Py_DECREF(temp);
   PyDict_SetItemString( AnalysisDataDict, "aVariance",
			 temp=PyFloat_FromDouble(aVariance) ); Py_DECREF(temp);

   Xlist     = PyList_New(0);
   VCElist   = PyList_New(0);
   vVCEMlist = PyList_New(0);
   for (ii = 0; ii < nInputs; ii++) {
      PyList_Append( VCElist, temp=PyFloat_FromDouble(vce[ii]) );
      Py_DECREF(temp);
      PyList_Append( vVCEMlist, temp=PyFloat_FromDouble(varVCEMean[ii]) );
      Py_DECREF(temp);

      XIlist = PyList_New(0);
      for (ss = 0; ss < nSamples; ss++) {
         PyList_Append( XIlist, temp=PyFloat_FromDouble(X[nInputs*ss+ii]) );
         Py_DECREF(temp);
      }
      PyList_Append( Xlist, XIlist ); Py_DECREF(XIlist);
   }
   PyDict_SetItemString( AnalysisDataDict, "X", Xlist ); Py_DECREF(Xlist);
   PyDict_SetItemString( AnalysisDataDict, "vce", VCElist ); Py_DECREF(VCElist);
   PyDict_SetItemString( AnalysisDataDict, "varVCEMean", vVCEMlist );
   Py_DECREF(vVCEMlist);

   Ylist = PyList_New(0);
   for (ss = 0; ss < nSamples; ss++) {
      PyList_Append( Ylist, temp=PyFloat_FromDouble(Y[ss]) ); Py_DECREF(temp);
   }
   PyDict_SetItemString( AnalysisDataDict, "Y",	Ylist ); Py_DECREF(Ylist);
#endif

   if (fp != NULL)
   {
      fwritePlotCLF(fp);
      fprintf(fp, "nx = ceil(sqrt(%d));", nInputs);
      for (ii = 0; ii < nInputs; ii++)
      {
         fprintf(fp, "subplot(nx, nx, %d)\n", ii+1);
         if (psPlotTool_ == 1)
            fprintf(fp, "Y = matrix(Y%d_%d,nReps%d,nSymbols%d);\n",
                    whichOutput,ii,ii,ii);
         else
            fprintf(fp, "Y = reshape(Y%d_%d,nReps%d,nSymbols%d);\n",
                    whichOutput,ii,ii,ii);
         fprintf(fp, "ymax = max(Y%d_%d);\n", whichOutput, ii);
         fprintf(fp, "ymin = min(Y%d_%d);\n", whichOutput, ii);
         fprintf(fp, "xmax = max(X%d);\n", ii);
         fprintf(fp, "xmin = min(X%d);\n", ii);
         if (psPlotTool_ == 1)
         {
            fprintf(fp, "a = get(\"current_axes\");\n");
            fprintf(fp, "a.data_bounds=[xmin,ymin;xmax,ymax];\n");
         }
         else
         {
            fprintf(fp, "axis([xmin xmax ymin ymax])\n");
         }
         fprintf(fp, "for k = 1 : nSymbols%d \n",ii);
         fprintf(fp, "   yt = Y(:, k);\n");
         fprintf(fp, "   mean = sum(yt) / nReps%d;\n",ii);
         fprintf(fp, "   x = X%d(k) * ones(nReps%d, 1);\n",ii,ii);
         fprintf(fp, "   plot(x, yt, '.', 'MarkerSize', 8)\n");
         if (psPlotTool_ == 0)
         {
            fprintf(fp, "   if ( k == 1 )\n");
            fprintf(fp, "      hold on\n");
            fprintf(fp, "   end;\n");
         }
         fprintf(fp, "   plot(x(1), mean, 'r*', 'MarkerSize', 12)\n");
         fprintf(fp, "   meanArray(k) = mean;\n");
         fprintf(fp, "   tmaxArray(k) = max(yt);\n");
         fprintf(fp, "   tminArray(k) = min(yt);\n");
         fprintf(fp, "end;\n");
         fwritePlotAxes(fp);
         fwritePlotXLabel(fp, "Input");
         fwritePlotYLabel(fp, "Output");
         fprintf(fp, "meanMean = sum(meanArray)/nSymbols%d;\n",ii);
         fprintf(fp, "adjMean  = (meanArray - meanMean);\n");
         fprintf(fp, "varMean = sum(adjMean.^2)/nSymbols%d\n",ii);
         fprintf(fp, "[xx,I] = sort(X%d);\n", ii);
         fprintf(fp, "yy1 = meanArray(I);\n");
         fprintf(fp, "yy2 = tmaxArray(I);\n");
         fprintf(fp, "yy3 = tminArray(I);\n");
         fprintf(fp, "plot(xx, yy1, 'r-', 'lineWidth', 2.0)\n");
         fprintf(fp, "plot(xx, yy2, 'r-', 'lineWidth', 2.0)\n");
         fprintf(fp, "plot(xx, yy3, 'r-', 'lineWidth', 2.0)\n");
         fprintf(fp, "title('Output %d vs Input %d')\n",
                 whichOutput+1,ii+1);
         fflush(fp);
      }
      fclose(fp);
      printf("Main effect matlab plot %s has been generated.\n", meFileName);
   }

   if (psAnaExpertMode_ == 1)
   {
      printf("Bootstrap analysis takes the sample, replicates it n times,\n");
      printf("and assess whether the sensitivity indices have converged.\n");
      printf("If you are performed iterative analysis with refinements,\n");
      printf("you will need to enter 'no index reuse' below at the first\n");
      printf("iteration and 'yes' afterward until the final refinement.\n");
      sprintf(pString,"Perform bootstrap main effect analysis? (y or n) ");
      getString(pString, winput1);
      ncount = 0;
      if (winput1[0] == 'y')
      {
         sprintf(pString,"Number of bootstrap samples to use (>=100): ");
         ncount = getInt(100, 2000, pString);
         XX = new double[nInputs*nSamples];
         YY = new double[nSamples];
         bsVCEs = new double*[ncount];
         for (ii = 0; ii < ncount; ii++) bsVCEs[ii] = new double[nInputs];
         nReplications = nSamples / nSubSamples;
         winput1[0] = 'n';
         fp1 = fopen(".ME_bootstrap_indset", "r");
         if (fp1 != NULL)
         {
            printf(".ME_bootstrap_indset file found.\n");
            sprintf(pString,"Re-use file? (y or n) ");
            getString(pString, winput1);
            if (winput1[0] == 'y')
            {
               fscanf(fp1, "%d", &ii);
               if (ii != nReplications*ncount)
               {
                  printf("ERROR: expect the first line to be %d.\n",
                         nReplications*ncount);
                  printf("       Instead found the first line to be %d.n",
                         ii);
                  exit(1);
               }
            }
            else
            {
               fclose(fp1);
               fp1 = fopen(".ME_bootstrap_indset", "w");
               if (fp1 == NULL)
               {
                  printf("ERROR: cannot open ME_bootstrap_indset file.\n");
                  exit(1);
               }
               fprintf(fp1, "%d\n", nReplications*ncount);
            }
         }
         else
         {
            fp1 = fopen(".ME_bootstrap_indset", "w");
            if (fp1 == NULL)
            {
               printf("ERROR: cannot open .ME_bootstrap_indset file.\n");
               exit(1);
            }
            fprintf(fp1, "%d\n", nReplications*ncount);
         }
         fp = NULL;
         for (ii = 0; ii < ncount; ii++)
         {
            for (ir = 0; ir < nReplications; ir++)
            {
               if (fp1 != NULL && winput1[0] == 'y')
               {
                  fscanf(fp1, "%d\n", &index);
                  if (index < 0 || index >= nReplications)
                  {
                     printf("ERROR: reading index from file .ME_bootstrap_indset\n");
                     printf("       index read = %d\n", index);
                     printf("       expected   = [0,%d]\n", nReplications-1);
                  }
               }
               else 
                  index = PSUADE_rand() % nReplications;
               if (fp1 != NULL && winput1[0] != 'y')
                  fprintf(fp1, "%d\n", index);
               for (ss = 0; ss < nSubSamples; ss++)
               { 
                  for (jj = 0; jj < nInputs; jj++)
                     XX[(ir*nSubSamples+ss)*nInputs+jj] =
                        X[(index*nSubSamples+ss)*nInputs+jj];
                  YY[ir*nSubSamples+ss] = Y[index*nSubSamples+ss];
               }
            }
            computeMeanVariance(nInputs,1,nSamples,YY,&aMean,
                                &aVariance,0);
            computeVCE(nInputs, nSamples, nSubSamples, XX, YY, 
                       iZero, fp, varVCEMean, meanVCEVar, 
                       varVCEVar, vce);
            for (jj = 0; jj < nInputs; jj++) 
               bsVCEs[ii][jj] = vce[jj]/aVariance;
         }
         if (fp1 != NULL) fclose(fp1);
         printf("Enter name of matlab file to store bootstrap info: ");
         scanf("%s", winput1);
         fgets(winput2,500,stdin);
         fp = fopen(winput1, "w");
         if (fp != NULL)
         {
            fprintf(fp, "%%bootstrap sample of main effects\n");
            fprintf(fp, "%%VCE(i,j) = VCE for input i, bs sample j\n");
            fprintf(fp, "VCE = zeros(%d,%d);\n",nInputs,ncount);
            for (ii = 0; ii < ncount; ii++)
               for (jj = 0; jj < nInputs; jj++)
                  fprintf(fp, "VCE(%d,%d) = %e;\n",jj+1,ii+1,
                          bsVCEs[ii][jj]);
            fclose(fp);
         }
         else
         {
            printf("ERROR: cannot open file %s\n", winput1);
         }
         for (ii = 0; ii < ncount; ii++) delete [] bsVCEs[ii];
         delete [] bsVCEs;
         delete [] XX;
         delete [] YY;
      }
   }

   pData *pPtr = ioPtr->getAuxData();
   pPtr->nDbles_ = nInputs;
   pPtr->dbleArray_ = new double[nInputs];
   for (ii = 0; ii < nInputs; ii++) pPtr->dbleArray_[ii] = vce[ii];
   pPtr->dbleData_ = aVariance;

   delete [] Y;
   delete [] vce;
   delete [] varVCEMean;
   delete [] meanVCEVar;
   delete [] varVCEVar;
   delete [] txArray;
   delete [] tyArray;
   if (constrPtr != NULL) delete constrPtr;

   // return 1.0 to facilitate continuous refinement
   return 1.0;
}

// *************************************************************************
// Compute agglomerated mean and variance
// -------------------------------------------------------------------------
int MainEffectAnalyzer::computeMeanVariance(int nInputs, int nOutputs, 
                        int nSamples, double *Y, double *aMean, 
                        double *aVariance, int outputID)
{
   int    ss, count;
   double mean, variance;

   count = 0;
   mean = 0.0;
   for (ss = 0; ss < nSamples; ss++)
   {
      if (Y[nOutputs*ss+outputID] < 0.9*PSUADE_UNDEFINED)
      {
         mean += Y[nOutputs*ss+outputID];
         count++;
      }
   }
   if (count <= 0)
   {
      printf("MainEffectAnalyzer ERROR: no valid data.\n");
      exit(1);
   }
   mean /= (double) count;
   variance = 0.0;
   for (ss = 0; ss < nSamples; ss++) 
   {
      if (Y[nOutputs*ss+outputID] < 0.9*PSUADE_UNDEFINED)
      {
         variance += ((Y[nOutputs*ss+outputID] - mean) *
                      (Y[nOutputs*ss+outputID] - mean));
      }
   }
   variance /= (double) count;
   (*aMean) = mean;
   (*aVariance) = variance;
   return 0;
}

// *************************************************************************
// Compute VCE 
// -------------------------------------------------------------------------
int MainEffectAnalyzer::computeVCE(int nInputs, int nSamples, int nSubSamples, 
                                   double *X, double *Y, int whichOutput, 
                                   FILE *fp, double *varVCEMean,
                                   double *meanVCEVar, double *varVCEVar,
                                   double *vce)
{
   int    ii, ss, nReplications, nSubs, ncount, symIndex, repID, subID;
   int    validnReps, validnSubs, totalCnt, *bins;
   double *txArray, *tyArray, *vceMean, *vceVariance, fixedInput, aMean;

   bins    = new int[nSubSamples];
   txArray = new double[nSamples];
   tyArray = new double[nSamples];
   vceMean     = new double[nSubSamples];
   vceVariance = new double[nSubSamples];

   for (ii = 0; ii < nInputs; ii++)
   {
      if (fp != NULL) fprintf(fp, "Y%d_%d = [\n",whichOutput,ii);

      for (ss = 0; ss < nSamples; ss++)
      {
         txArray[ss] = X[nInputs*ss+ii];
         tyArray[ss] = Y[ss];
      }
      sortDbleList2(nSamples, txArray, tyArray);
      nReplications = 1;
      for (ss = 1; ss < nSamples; ss++)
      {
         if (txArray[ss] == txArray[0]) nReplications++;
         else                           break;
      }
      nSubs = nSamples / nReplications;

      ncount = 0;
      for (ss = 0; ss < nSamples; ss+=nReplications)
      {
         symIndex = ss / nReplications;
         fixedInput = txArray[ss];
         vceMean[symIndex] = 0.0;
         validnReps = 0;
         for (repID = 0; repID < nReplications; repID++)
         {
            if (PABS(txArray[ss+repID] - fixedInput) < 1.0E-8)
            {
               ncount++;
               if (tyArray[ss+repID] < 0.9*PSUADE_UNDEFINED)
               {
                  validnReps++;
                  vceMean[symIndex] += tyArray[ss+repID];
                  if (fp != NULL)
                     fprintf(fp,"%24.16e\n",tyArray[ss+repID]);
               }
            }
         }
         bins[symIndex] = validnReps;
         if (validnReps > 0) vceMean[symIndex] /= (double) validnReps;
         else                vceMean[symIndex] = PSUADE_UNDEFINED;
         vceVariance[symIndex] = 0.0;
         for (repID = 0; repID < nReplications; repID++)
         {
            if (vceMean[symIndex] != PSUADE_UNDEFINED &&
                tyArray[ss+repID] < 0.9*PSUADE_UNDEFINED)
            {
               vceVariance[symIndex] += 
                    ((tyArray[ss+repID]-vceMean[symIndex])*
                     (tyArray[ss+repID]-vceMean[symIndex]));
            }
         }
         if (validnReps > 0)
              vceVariance[symIndex] /= (double) validnReps;
         else vceVariance[symIndex] = PSUADE_UNDEFINED;
      }
      if (fp != NULL)
      { 
         fprintf(fp,"];\n");
         fprintf(fp, "nSymbols%d = %d;\n", ii, nSubs);
         fprintf(fp, "nReps%d    = %d;\n", ii, nReplications);
      }

      if (ncount != nSamples)
      {
         printf("MainEffectAnalyzer ERROR: not rLHS, input %d\n",
                ii+1);
         printf("       error data = %d (%d)\n",ncount, nSamples);
         printf("       Did you use rLHS ?\n");
         delete [] bins;
         delete [] txArray;
         delete [] tyArray;
         delete [] vceMean;
         delete [] vceVariance;
         return -1;
      }

      totalCnt = 0;
      for (subID = 0; subID < nSubs; subID++) totalCnt += bins[subID];
      varVCEMean[ii] = 0.0;
      meanVCEVar[ii] = 0.0;
      aMean = 0.0;
      validnSubs = 0;
      for (subID = 0; subID < nSubs; subID++)
      {
         if (vceVariance[subID] != PSUADE_UNDEFINED)
         {
            aMean += (vceMean[subID] / totalCnt * bins[subID]);
            validnSubs++;
         }
      }
      for (subID = 0; subID < nSubs; subID++)
      {
         if (vceVariance[subID] != PSUADE_UNDEFINED)
         {
            varVCEMean[ii] += ((vceMean[subID] - aMean) *
                               (vceMean[subID] - aMean) * bins[subID]/totalCnt);
            meanVCEVar[ii] += vceVariance[subID] * bins[subID] / totalCnt;
         }
      }
      varVCEVar[ii] = 0.0;
      for (subID = 0; subID < nSubs; subID++)
      {
         if (vceVariance[subID] != PSUADE_UNDEFINED)
            varVCEVar[ii] += ((vceVariance[subID] - meanVCEVar[ii]) *
                              (vceVariance[subID] - meanVCEVar[ii]) * 
                              bins[subID] / totalCnt);
      }
      if (validnSubs > 0)
         vce[ii] = varVCEMean[ii] - meanVCEVar[ii]/((double) validnSubs);
      else
         vce[ii] = 0.0;
   }
   delete [] txArray;
   delete [] tyArray;
   delete [] vceMean;
   delete [] vceVariance;
   delete [] bins;
   return 0;
}

// *************************************************************************
// Compute VCE 
// -------------------------------------------------------------------------
int MainEffectAnalyzer::computeVCECrude(int nInputs, int nSamples, 
                          double *X, double *Y, double *iLowerB, 
                          double *iUpperB, double aVariance, double *vce)
{
   int    ii, ss, nSize, index;
   int    totalCnt, *bins, nIntervals, *tags;
   double *vceMean, *vceVariance, aMean, ddata, hstep;
   char   pString[500];
   FILE   *fp;

   nSize = 50;
   nIntervals = nSamples / nSize;
   printAsterisks(0);
   printf("                 Crude Main Effect\n");
   printDashes(0);
   printf("MainEffectAnalyzer: number of intervals  = %d\n", nIntervals);
   printf("MainEffectAnalyzer: sample size/interval = %d\n", nSize);
   printf("These may need to be adjusted for higher accuracy.\n");
   printf("Turn on analysis expert mode to change them.\n");
   if (psAnaExpertMode_ == 1)
   {
      sprintf(pString,"number of levels (>5, default = %d): ", nIntervals);
      nIntervals = getInt(5, nSamples/2, pString);
      nSize = nSamples / nIntervals;
   }
   printEquals(0);
   bins = new int[nIntervals];
   tags = new int[nSamples];
   vceMean     = new double[nIntervals];
   vceVariance = new double[nIntervals];

   for (ii = 0; ii < nInputs; ii++)
   {
      hstep = (iUpperB[ii] - iLowerB[ii]) / nIntervals;
      for (ss = 0; ss < nIntervals; ss++)
      {
         vceMean[ss] = 0.0;
         vceVariance[ss] = 0.0;
         bins[ss] = 0;
      }
      for (ss = 0; ss < nSamples; ss++)
      {
         ddata = (X[nInputs*ss+ii] - iLowerB[ii]) / hstep;
         index = (int) ddata;
         if (index >= nIntervals) index = nIntervals - 1;
         tags[ss] = -1;
         if (Y[ss] < 0.9*PSUADE_UNDEFINED)
         {
            vceMean[index] += Y[ss];
            bins[index]++;
            tags[ss] = index;
         }
      }
      for (ss = 0; ss < nIntervals; ss++)
         if (bins[ss] > 0) vceMean[ss] /= (double) bins[ss];
      for (ss = 0; ss < nSamples; ss++)
      {
         index = tags[ss];
         if (Y[ss] < 0.9*PSUADE_UNDEFINED)
         {
            vceVariance[index] += ((Y[ss]-vceMean[index])*
                                   (Y[ss]-vceMean[index]));
         }
      }
      for (ss = 0; ss < nIntervals; ss++)
      {
         if (bins[ss] > 0)
              vceVariance[ss] /= (double) bins[ss];
         else vceVariance[ss] = PSUADE_UNDEFINED;
      }
         
      totalCnt = 0;
      for (ss = 0; ss < nIntervals; ss++) totalCnt += bins[ss];
      vce[ii] = 0.0;
      aMean = 0.0;
      for (ss = 0; ss < nIntervals; ss++)
      {
         if (vceVariance[ss] != PSUADE_UNDEFINED)
            aMean += (vceMean[ss] / totalCnt * bins[ss]);
      }
      for (ss = 0; ss < nIntervals; ss++)
      {
         if (vceVariance[ss] != PSUADE_UNDEFINED)
         {
            vce[ii] += ((vceMean[ss] - aMean) *
                        (vceMean[ss] - aMean) * bins[ss]/totalCnt);
         }
      }
   }

   ddata = 0.0;
   for (ii = 0; ii < nInputs; ii++)
   {
      ddata += vce[ii] / aVariance;
      printf("    INPUT %2d = %9.2e (raw = %9.2e)\n", ii+1, 
             vce[ii]/aVariance, vce[ii]);
   }
   printf("    Total VCE = %9.2e\n", ddata);
   printAsterisks(0);

   fp = fopen("matlabmec.m", "w");
   fwritePlotCLF(fp);
   fprintf(fp, "A = [\n");
   for (ii = 0; ii < nInputs; ii++)
      fprintf(fp, "%e\n", vce[ii]/aVariance);
   fprintf(fp, "];\n");
   fprintf(fp, "bar(A, 0.8);\n");
   fwritePlotAxes(fp);
   fwritePlotTitle(fp, "Approximate VCE Rankings");
   fwritePlotXLabel(fp, "Input parameters");
   fwritePlotYLabel(fp, "approximate VCE");
   fclose(fp);
   printf("Main Effect approximate ranking is now in matlabmec.m.\n");

   delete [] vceMean;
   delete [] vceVariance;
   delete [] bins;
   delete [] tags;
   return 0;
}

