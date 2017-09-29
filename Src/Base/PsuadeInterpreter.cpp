// ************************************************************************
// Copyright (c) 2007   Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory.
// Written by the PSUADE team.
// All rights reserved.
//
// Please see the COPYRIGHT_and_LICENSE file for the copyright notice,
// disclaimer, contact information and the GNU Lesser General Public License.
//
// PSUADE is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License (as published by the Free Software
// Foundation) version 2.1 dated February 1999.
//
// PSUADE is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the terms and conditions of the GNU General
// Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
// ************************************************************************
// Functions for the class PsuadeBase
// AUTHOR : CHARLES TONG
// DATE   : 2005
// ************************************************************************

// ------------------------------------------------------------------------
// system includes
// ------------------------------------------------------------------------
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <PsuadeCmakeConfig.h>

#ifdef WINDOWS
#include <windows.h>
#endif

// ------------------------------------------------------------------------
// local includes : class definition and utilities
// ------------------------------------------------------------------------
#include "Psuade.h"
#include "PsuadeBase.h"
#include "dtype.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "PsuadeConfig.h"
#include "Matrix.h"
#include "Vector.h"

// ------------------------------------------------------------------------
// local includes : function approximator and others
// ------------------------------------------------------------------------
#include "FuncApprox.h"
#include "AnalysisManager.h"
#include "TSIAnalyzer.h"
#include "SobolAnalyzer.h"
#include "PDFManager.h"
#include "Sampling.h"
#include "FunctionInterface.h"
#include "PsuadeData.h"
#include "Optimizer.h"
#include "PsuadeSession.h"
#include "PrintingTS.h"
#include "PDFHistogram.h"

// ------------------------------------------------------------------------
// local defines 
// ------------------------------------------------------------------------
#define PABS(x)  ((x) > 0 ? x : -(x))

// ************************************************************************
// interpret command from interactive session
// ------------------------------------------------------------------------
int PsuadeBase::interpretInteractive()
{
   int    status, nReps, iInd, outputID, ind, ind2, iInd1, iInd2, count;
   int    it, nPtsPerDim=64, iplot1, iplot2, jplot=-1, faLeng, sInd, flag;
   int    ss, ii, jj, kk, ll, iplot3, *states, rsiNOutputs, *rsiSet, count2;
   int    **rsiMatrix, analysisMethod, samplingMethod, faType, setCompare;
   int    otrans, *xsforms, nOut2, oInd, *tags, faFlag, iOne=1;
   int    nPaths, nTrials=100, trial, *indSet, nSmooth, currP, faLimit;
   int    nParts, nPlots, *plotIndices=NULL, faID, nFiles, iplot4;
   int    inputID, oplot1, oplot2, scriptMode=0, *tempI;
   int    method, saveDiag, saveMode1, saveMode2, saveMethod, usePDFs=0;
   int    nParams, commandCnt, iSave;
   int    *SPDFIndices;
   long   nSamplesLong;
   double *tempX, *tempY, *tempW, Ymax, Ymin, Xmin, Xmax, width;
   double *tempInds=NULL, *tempT=NULL, *inputSettings=NULL, *faYIn=NULL;
   double *faXOut=NULL, *faYOut=NULL;
   double thresh, threshL, threshU, GYmax, GYmin, dtemp;
   double *tempV, **moatSample=NULL, currX1, currX2, currY1, gamma;
   double currY2, ddata, filterRange, sLo, sHi;
#ifdef HAVE_SVM
   double sumErr, maxErr, tolerance;
#endif
   double *threshLs, *threshUs, diagMax, diagMin, currX1F, currX2F;
   double ***rsi3Matrix,aVal,bVal,minDist;
   char   command[5001],dataFile[5001],**names,**inputPDFFiles;
   char   winput[50001], lineIn[50001], pString[50001];
   char   dirName[5001],errFile[5001],*targv[8],scriptName[5001], cString[5001];
   char   subdirName[5001], **strptr;
   FILE   *fp, *fpOut, *fErr, *scriptFp=NULL;
   FuncApprox *faPtr=NULL, **faPtrs=NULL, **faPtrsRsEval=NULL;
   pData  pPtr, pINames, pLower, pUpper, pONames, pInputs, pOutputs, pStates;
   pData  *pdata, pPDFs, pMeans, pStds, pPDFSIndices;
   aData  aPtr;
   PsuadeData *ioPtr=NULL;
   Sampling   *sampPtr=NULL, *sampAux=NULL;
   FunctionInterface *funcIO=NULL;
   AnalysisManager   *anaManager=NULL;
   PDFManager *pdfman=NULL;
   Matrix     corMat;
   Vector     vecIn, vecOut, vecUpper, vecLower;
   PsuadeSession *currSession=NULL, *newSession=NULL;

   // loop on the command interpreter
   printf("PSUADE - A Problem Solving environment for \n");
   printf("         Uncertainty Analysis and Design Exploration (%d.%d.%da)\n",
          psuade_VERSION_MAJOR, psuade_VERSION_MINOR, psuade_VERSION_PATCH);
   printf("(for help, enter <help>)\n");
   printEquals(PL_INFO, 0);
   commandCnt = 0;
   currSession = new PsuadeSession();
   while (1)
   {
      if (scriptMode == 1 && scriptFp != NULL)
      {
         fgets(lineIn, 5000, scriptFp);
         if (feof(scriptFp) != 0)
         {
            fclose(scriptFp);
            scriptMode = 0;
         }
         for (ii = 0; ii < 100; ii++) command[ii] = '\0';
         sscanf(lineIn, "%s", command);
         printf("script> %s\n", command);
      }
      else
      {
         printf("psuade> ");
         for (ii = 0; ii < 500; ii++) lineIn[ii] = '\0';
         fgets(lineIn,5000,stdin); 
         winput[0] = '\0';
         pString[0] = '\0';
         command[0] = '\0';
         sscanf(lineIn, "%s", command);
      }
      if (!strcmp(command, "\0"))
      {
         commandCnt++;
         if (commandCnt >= 10)
         {
            printf("Enter carriage return > 10 times ==> terminate.\n");
            return 0;
         }
      }
      else commandCnt = 0;
      if (!strcmp(command, "help") || !strcmp(command, "h"))
      {
         strcpy(winput, "\0");
         sscanf(lineIn, "%s %s %s", command, winput, pString);
         if (!strcmp(winput, "info"))
         {
            printf("Useful information for using PSUADE :\n");
            printf("\tI.    Uncertainty analysis: \n");
            printf("\t\t 1. Sampling: MC, LPTAU, METIS, LH, OA, or OALH\n");
            printf("\t\t    Analyzer: use 'ua' in command line mode\n");
            printf("\t\t 2. Sampling: first construct response surface (rs)\n");
            printf("\t\t    Analyzer: use rsua\n");
            printf("\t\t 3. Sampling: first construct response surface (rs)\n");
            printf("\t\t    Analyzer: use rsuab (rsua with bootstrapping)\n");
            printf("\tII.   Parameter Screening: \n");
            printf("\t\t 1. Sampling: MOAT, GMOAT\n");
            printf("\t\t    Analyzer: moat, moatmo\n");
            printf("\t\t 2. Sampling: LH, LPTAU\n");
            printf("\t\t    Analyzer: mars_sa/gp_sa (MARS/GP screening)\n");
            printf("\t\t 3. Sampling: LH, LPTAU\n");
            printf("\t\t    Analyzer: delta_test (large samples, low dimension)\n");
            printf("\t\t 4. Sampling: LH, MC\n");
            printf("\t\t    Analyzer: sot_sa (sum-of-trees screening)\n");
            printf("\t\t 5. Sampling: FF4, FF5\n");
            printf("\t\t    Analyzer: ff (Fractional Factorial for 2-way analysis)\n");
            printf("\t\t 6. Sampling: LSA\n");
            printf("\t\t    Analyzer: lsa (local sensitivity analysis - linear)\n");
            printf("\tIII.  Classical regression/sensitivity analysis: \n");
            printf("\t\t 1. Sampling: MC, LPTAU, METIS, LH, OA, or OALH\n");
            printf("\t\t    Regression-based correlation analysis (use ca)\n");
            printf("\t\t 2. Sampling: MC, LPTAU, METIS, LH, OA, or OALH\n");
            printf("\t\t    Polynomial regression-based (use rscheck/examine SRCs)\n");
            printf("\tIV.   Response surface analysis: \n");
            printf("\t\t 1. Sampling: LPTAU, METIS, LH, OA, or OALH\n");
            printf("\t\t    Response surface validation: a few options\n");
            printf("\t\t    (a) examine R-squared (may not be reliable)\n");
            printf("\t\t    (b) test it on the training set (use rstest_ts)\n");
            printf("\t\t    (c) test it on another hold-out set (use rstest_ts)\n");
            printf("\t\t    (d) perform cross validation (use rstest_cv)\n");
            printf("\t\t    (e) perform generalization test (rstest_gt)\n");
            printf("\tV.   Global sensitivity analysis (first order): \n");
            printf("\t\t 1. Sampling: replicated Latin hypercube\n");
            printf("\t\t    Analyzer: use 'me'\n");
            printf("\t\t 2. Sampling: any space-filling (large enough) sample\n");
            printf("\t\t    Analyzer: use 'me'\n");
            printf("\t\t 3. Sampling: Sobol' sampling method (large)\n");
            printf("\t\t    Analyzer: use 'sobol'\n");
            printf("\t\t 4. Sampling: any space-filling sample\n");
            printf("\t\t    Analyzer: use rscheck with Legendre and select scaling\n");
            printf("\t\t 5. Sampling: first construct response surface\n");
            printf("\t\t    Analyzer: use rssobol1b or rsmeb (faster)\n");
            printf("\tVI.  Global sensitivity analysis (second order): \n");
            printf("\t\t 1. Sampling: replicated OA\n");
            printf("\t\t    Analyzer: use 'ie' in command line mode\n");
            printf("\t\t 2. Sampling: any space-filling (large enough) sample\n");
            printf("\t\t    Analyzer: use 'ie'\n");
            printf("\t\t 3. Sampling: first construct response surface\n");
            printf("\t\t    Analyzer: use rssobol2b or rsieb (less robust)\n");
            printf("\tVII.  Global sensitivity analysis (total order): \n");
            printf("\t\t 1. Sampling: use any space-filling (very large) sample\n");
            printf("\t\t    Analyzer: use tsi (coarse analysis for low dimensions)\n");
            printf("\t\t 2. Sampling: Sobol' sampling method (large)\n");
            printf("\t\t    Analyzer: use 'sobol'\n");
            printf("\t\t 3. Sampling: first construct response surface\n");
            printf("\t\t    Analyzer: use rssoboltsib in command line mode, or\n");
            printf("\t\t 4. Sampling: use FAST sampling\n");
            printf("\t\t    Analyzer: use FAST analyzer\n");
            printf("\tVIII. Global sensitivity analysis (group): \n");
            printf("\t\t 1. Sampling: first construct response surface\n");
            printf("\t\t    Analyzer: use rssobolg in command line mode\n");
            printf("\tIX.   Hypothesis testing: \n");
            printf("\t\t 1. Sampling: any of your choice \n");
            printf("\t\t    Analyzer: 1test or 2test\n");
            printf("\tX.    Principal component analysis: \n");
            printf("\t\t 1. Sampling: any of your choice \n");
            printf("\t\t    Analyzer: pca in command line mode\n");
            printf("\tXI.   Bayesian inverse UQ: (response surface-based)\n");
            printf("\t\t 1. Sampling: LPTAU, LH, OA, METIS for response surface\n");
            printf("\t\t    Analyzer: rsmcmc (in command line mode), or\n");
            printf("\t\t 2. Sampling: LPTAU, LH, OA, METIS for response surface\n");
            printf("\t\t    Analyzer: mcmc (need re-compilation for efficiency)\n");
            printf("\tXII.  Advanced features: \n");
            printf("\t\t 1. Impose constraints in sampling/analysis (e.g. moatgen)\n");
            printf("\t\t 2. Plot PDF of std dev. of response surfaces (rssd_ua)\n");
            printf("\t\t 3. Intersection or Bayes-like rules (e.g. rsi2, rsvol)\n");
            printf("\t\t 4. Multi-objective optimization (mo_opt)\n");
            printf("\t\t 5. Mixed aleatoric-epistemic uncertainty analysis (aeua)\n");
            printf("\t\t 6. 2nd order analysis (soua) - uncertainty in input PDFs\n");
            printf("\t\t 7. Tools for setting up user application with PSUADE\n");
            printf("\t\t 8. Sample refinement (uniform/adaptive: refine, a_refine)\n");
         }
         else if (!strcmp(winput, "io"))
         {
            printf("Commands for reading/write/updating data to/from files:\n");
            printf("(to see details of each command, use '-h' option)\n");
            printf("\tload         <file> (load a data file in PSUADE format) \n");
            printf("\tloadmore     <file> (add data to current data set)\n");
            printf("\twrite        <file> (write to file in PSUADE format)\n");
            printf("\tnwrite       <file> (write unevaluated points only)\n");
            printf("\tiread        <file> (read from file with only inputs)\n");
            printf("\tiwrite       <file> (write to file with only inputs)\n");
            printf("\towrite       <file> (write to file with only outputs)\n");
            printf("\tread_std     <file> (read data in standard format)\n");
            printf("\twrite_std    <file> (write to file in standard format)\n");
            printf("\tread_csv     <file> (read data in special CSV format)\n");
            printf("\twrite_matlab <file> (write to file in matlab format)\n");
            printf("\twrite_ultra  <file> (write to file in ultra format)\n");
            printf("\tread_xls     <file> (read data in Excel format)\n");
            printf("\twrite_xls    <file> (write to file in Excel format)\n");
            printf("\tupdate       <file> (update sample OUTPUTS from a data file)\n");
            printf("\tiadd         <file> (add more inputs from a data file)\n");
            printf("\toadd         <file> (add more outputs from a data file)\n");
            printf("\tiadd1        <file> (add 1 input with random data in [0,1])\n");
            printf("\toadd1        <file> (add 1 output with random data in [0,1])\n");
            printf("\tireplace     <file> (replace one input from <file>)\n");
            printf("\toreplace     <file> (replace all outputs from <file>)\n");
            printf("\tsplitsample         (split sample into 2 files)\n");
         }
         else if (!strcmp(winput, "stat"))
         {
            printf("Commands for basic statistic on raw sample data:\n");
            printf("(to see details of each command, use '-h' option)\n");
            printf("\tua         (uncertainty analysis on any sample)\n");
            printf("\tca         (correlation analysis on any sample)\n");
            printf("\tme         (main effect analysis on any sample)\n");
            printf("\tie         (2-way interaction effect analysis on any sample)\n");
            printf("\ttsi        (total sensitivity analysis on any sample)\n");
            printf("\tsobol      (Sobol sensitivity analysis on Sobol sample)\n");
            printf("\tfast       (total sensitivity analysis using FAST sample)\n");
            printf("\tanova      (analyis of variation with MARS response surface)\n");
            printf("\t1stest     (1-sample test (Chi-squared, dist fit))\n");
            printf("\t2stest     (2-sample test (T-test,K-S,Mann-Whitney))\n");
            printf("\tgendist    (create a sample given an PDF (for single input))\n");
            printf("\tpdfconvert (convert a sample using selected distribution)\n");
            printf("\trand_draw  (draw a bootstrap from the loaded sample)\n");
            printf("\trand_draw2 (draw a bootstrap from 2 samples of 2 input sets)\n");
            printf("\tgensample  (create a sample using PDF info from the loaded file)\n");
            printf("\tcdf_lookup (look up cumulative probability given a value)\n");
         }
         else if (!strcmp(winput, "screen"))
         {
            printf("Commands for parameter screening:\n");
            printf("(to see details of each command, use '-h' option)\n");
            printf("\tlsa        (parameter screening with local SA method)\n");
            printf("\tmoat       (parameter screening with the Morris method)\n");
            printf("\tmoatmo     (Morris screening for multiple outputs)\n");
            printf("\tff         (parameter screening with fractional factorial)\n");
            printf("\tmars_sa    (parameter screening with MARS)\n");
            printf("\tgp_sa      (parameter screening with Gaussian process)\n");
            printf("\tdelta_test (parameter screening with Delta test)\n");
            printf("\tsot_sa     (parameter screening with sum-of-trees method)\n");
            printf("\tpca        (principal component analysis: on outputs)\n");
         }
         else if (!strcmp(winput, "rs"))
         { 
            printf("Commands for response surface analysis:\n");
            printf("(to see details of each command, use '-h' option)\n");
            printf("\tsvmfind     (search for good parameters for SVM)\n");
            printf("\trscheck     (check quality of RS with CV and rich graphics)\n");
            printf("\trstest_hs   (check quality of RS with another test set)\n");
            printf("\trstest_ts   (check quality of RS with the training set)\n");
            printf("\trstest_cv   (check quality of RS with cross validation)\n");
            printf("\trstest_gt   (check quality of RS with generalization test)\n");
            printf("\trscreate    (create response surface to be used by rseval)\n");
            printf("\tivec_create (create an input vector in local memory)\n");
            printf("\tivec_modify (modify an entry of local input vector)\n");
            printf("\tivec_show   (display values of the local input vector)\n");
            printf("\trseval      (evaluate RS with local vector or from file)\n");
            printf("\trseval_m    (use PSUADE as response surface server)\n");
            printf("\trs_splot    (create scatter plots on response surface)\n");
            printf("\trstgen      (generate a sample (FF/FACT) for rstest_hs)\n");
            printf("\trsvol       (compute volume in constrained region)\n");
            printf("\trsint       (compute volume under response surface)\n");
         }
         else if (!strcmp(winput, "qsa") && !strcmp(pString, "long"))
         {
            printf("Commands for RS-based uncertainty/sensitivity analysis:\n");
            printf("(to see details of each command, use '-h' option)\n");
            printf("\trsua       (RS-based UA on fuzzy response surface)\n");
            printf("\trsuab      (RS-based UA on response surface with bootstrap)\n");
            printf("\trsmeb      (RS-based McKay main effect with bootstrap)\n");
            printf("\trsieb      (RS-based McKay pairwise effect with bootstrap)\n");
            printf("\trssobol1   (RS-based Sobol' main effect)\n");
            printf("\trssobol2   (RS-based Sobol' interaction effect)\n");
            printf("\trssobolg   (RS-based Sobol' group main effect)\n");
            printf("\trssoboltsi (RS-based Sobol' total effect)\n");
            printf("\trssobol1b  (RS-based Sobol' main effect with bootstrap)\n");
            printf("\trssobol2b  (RS-based Sobol' 2-way analysis with bootstrap)\n");
            printf("\trssoboltsib(RS-based Sobol' total effect with bootstrap)\n");
            printf("\tsoua       (RS-based 2nd order analysis: uncertainty in PDF)\n");
            printf("\taeua       (RS-based mixed aleatoric-epistemic analysis)\n");
         }
         else if (!strcmp(winput, "qsa"))
         {
            printf("Commands for RS-based uncertainty/sensitivity analysis:\n");
            printf("(to see more qsa commands, use 'help qsa long')\n");
            printf("(to see details of each command, use '-h' option)\n");
            printf("\trsua       (RS-based UA on fuzzy response surface)\n");
            printf("\trsuab      (RS-based UA on response surface with bootstrap)\n");
            printf("\trsmeb      (RS-based McKay main effect with bootstrapping)\n");
            printf("\trsieb      (RS-based McKay pairwise effect with bootstrap)\n");
            printf("\trssobol1b  (RS-based Sobol' main effect with bootstrap)\n");
            printf("\trssobol2b  (RS-based Sobol' pairwise effect with bootstrap)\n");
            printf("\trssoboltsib(RS-based Sobol' total effect with bootstrap)\n");
            printf("\tsoua       (RS-based 2nd order analysis: uncertainty in PDF)\n");
            printf("\taeua       (RS-based mixed aleatoric-epistemic analysis)\n");
         }
         else if (!strcmp(winput, "calibration"))
         {
            printf("Commands for optimization/calibration:\n");
            printf("(to see details of each command, use '-h' option)\n");
            printf("\trsmcmc     (RS-based Bayesian inversion using MCMC)\n");
            printf("\tmo_opt     (RS-based multi-objective optimization)\n");
            printf("\tmcmc       (Simulation-based MCMC: advanced feature)\n");
         }
         else if (!strcmp(winput, "plot") && !strcmp(pString, "long"))
         { 
            printf("Commands for plotting sample data:\n");
            printf("(to see details of each command, use '-h' option)\n");
            printf("\tsplot      (scatter plot in matlab/scilab)\n");
            printf("\tsplot2     (2-input/1 output scatter plot in matlab/scilab)\n");
            printf("\tsplot3     (3-input/1 output scatter plot in matlab/scilab)\n");
            printf("\trs1        (1-input response surface in matlab/scilab)\n");
            printf("\trs2        (2-input response surface in matlab/scilab)\n");
            printf("\trs3        (3-input response surface in matlab)\n");
            printf("\trs3m       (3-input RS in matlab (movie), output in z-axis)\n");
            printf("\trs4        (4-input response surface in matlab (movie))\n");
            printf("\trssd       (std dev response surface in matlab)\n");
            printf("\trssd_ua    (create pdf of std dev for response surface)\n");
            printf("\trsi2       (2-input RS intersections in matlab)\n");
            printf("\trsi3       (3-input RS intersections in matlab)\n");
            printf("\trsi3m      (3-input RS intersections in matlab, 2D movie)\n");
            printf("\trawi2      (2-input intersections with raw data in matlab)\n");
            printf("\trawi3      (3-input intersections with raw data in matlab)\n");
            printf("\trspairs    (2-input response surface pairs in matlab)\n");
            printf("\trsipairs   (2-input intersecting response surfaces in matlab)\n");
            printf("\tiplot1     (1-input (input only) scatter plot)\n");
            printf("\tiplot2     (2-input (input only) scatter plot)\n");
            printf("\tiplot3     (3-input (input only) scatter plot)\n");
            printf("\tiplot4m    (4-input (input only) scatter plot): movie\n");
            printf("\tiplot2_all (all pairs 2-input plots)\n");
            //printf("\tiplot_pdf  (plot the sample PDF of selected inputs)\n");
            printf("\toplot2     (2-output scatter plot)\n");
            printf("\tihist      (create an input histogram matlab file)\n");
            printf("\tihist2     (create two-input histogram matlab file)\n");
            printf("\tohist      (create an output histogram matlab file)\n");
            printf("\tohist2     (create two-output histogram matlab file)\n");
            printf("\tiotrace    (plot inputs/outputs of each run)\n");
            //printf("\tmeplot  (plot main effects with Pgplot)\n");
            //printf("\tmeplot2 (plot main effect/interaction with Pgplot)\n");
            //printf("\trsplot  (plot response surface with Pgplot)\n");
         }
         else if (!strcmp(winput, "plot"))
         { 
            printf("Commands for plotting sample data:\n");
            printf("(to see more plot commands, use 'help plot long')\n");
            printf("\tsplot      (sample scatter plot versus each input in matlab)\n");
            printf("\tsplot2     (2-input/1 output scatter plot in matlab/scilab)\n");
            printf("\tsplot3     (3-input/1 output scatter plot in matlab/scilab)\n");
            printf("\trs1        (1-input response surface in matlab/scilab)\n");
            printf("\trs2        (2-input response surface in matlab/scilab)\n");
            printf("\trs3        (3-input response surface in matlab)\n");
            printf("\trs3m       (3-input RS in matlab (movie), output in z-axis)\n");
            printf("\trs4        (4-input response surface in matlab (movie))\n");
            printf("\trssd       (std dev response surface in matlab)\n");
            printf("\trssd_ua    (create pdf of std dev for response surface)\n");
            printf("\trsi2       (2-input RS intersections in matlab)\n");
            printf("\trsi3       (3-input RS intersections in matlab)\n");
            printf("\trsi3m      (3-input RS intersections in matlab, 2D movie)\n");
            printf("\trspairs    (all 2-input response surfaces in matlab)\n");
            printf("\tiplot1     (1-input (input only) scatter plot)\n");
            printf("\tiplot2     (2-input (input only) scatter plot)\n");
            printf("\tiplot3     (3-input (input only) scatter plot)\n");
            printf("\tiplot4m    (4-input (input only) scatter plot): movie\n");
            printf("\tiplot2_all (all pairs 2-input plots)\n");
            printf("\toplot2     (2-output scatter plot)\n");
            printf("\tihist      (create an input histogram matlab file)\n");
            printf("\tihist2     (create two-input histogram matlab file)\n");
            printf("\tohist      (create an output histogram matlab file)\n");
            printf("\tohist2     (create two-output histogram matlab file)\n");
         }
         else if (!strcmp(winput, "setup"))
         {
            printf("Commands for setting up/monitoring work flow:\n");
            printf("\tsetupguide     (instructions on how to set up application)\n");
            printf("\tgeninputfile   (create an input file for psuade)\n");
            printf("\tgenbatchfile   (create a LLNL-specific batch file)\n");
            printf("\tgendriver      (create an application driver)\n");
            printf("\tgenexample     (create a simple example for demonstration)\n");
            printf("\tchkjobs        (check job status and create a report)\n");
         }
         else if (!strcmp(winput, "edit"))
         {
            printf("Commands for manipulating/displaying data in local memory:\n");
            printf("\tvalidate      (validate certain sample outputs)\n");
            printf("\tinvalidate    (invalidate selected sample points)\n");
            printf("\tsrandomize    (randomize sample point orders)\n");
            printf("\tifilter       (take out points outside input bounds)\n");
            printf("\tofilter       (take out points outside output bounds)\n");
            printf("\tidelete       (delete one input from data)\n");
            printf("\todelete       (delete one output from data)\n");
            printf("\tsdelete       (delete one sample point)\n");
            printf("\tspurge        (take out invalid sample points)\n");
            printf("\tishuffle      (re-order inputs)\n");
            printf("\tiselect_index (select/re-order inputs based on input number)\n");
            printf("\tiselect_name  (select/re-order inputs based on input names)\n");
            printf("\tsshow         (display one sample point)\n");
            printf("\tsinfo         (display information on the loaded sample)\n");
            printf("\tlist1         (list 1 input/1 output data pair)\n");
            printf("\tlist2         (list 2 inputs/1 output data pair)\n");
            printf("\tlistall       (list all inputs/all outputs data)\n");
            printf("\tmax           (find the maximum (output) sample point)\n");
            printf("\tmin           (find the minimum (output) sample point)\n");
            printf("\tirerange      (change input range in data)\n");
            printf("\tireset        (reset a selected input to certain values)\n");
            printf("\toreset        (reset a selected output to certain values)\n");
            printf("\tifloor        (truncate certain input to integer)\n");
            printf("\ticeil         (round certain input to integer)\n");
            printf("\tiround        (round certain input to the nearest integer)\n");
            printf("\titran         (transform certain input (log or power))\n");
            printf("\totran         (transform certain output (log or power))\n");
            printf("\titag          (tag sample points based on input values)\n");
            printf("\totag          (tag sample points based on output values)\n");
            printf("\trm_dup        (take out duplicate sample points)\n");
            printf("\toop           (replace output1 = a * out2 + b * out3)\n");
            printf("\tsetdriver <s> (set the driver field to <s>)\n");
         }
         else if (!strcmp(winput, "misc"))
         {
            printf("Miscellaneous commands:\n");
            printf("\trename <file> <file2> (rename a file) \n");
            printf("\trun <file>     (run a psuade input script) \n");
            printf("\tquit (exit)    (terminate command line session)\n");
            printf("\trsmax <d>      (set maximum number of data points for RS)\n");
            printf("\tsys <command>  (execute a system command)\n");
            printf("\tprintlevel <d> (set print level)\n");
            printf("\tinteractive    (turn on/off interactive mode: obsolete)\n");
            printf("\toutput_file    (set the default output file name\n");
            printf("\tsqc            (sample quality check: distance metric)\n");
            printf("\tnna            (nearest neighbor analysis: for outliers\n");
            printf("\tshowformat     (show file formats for RS/MOAT constraints)\n");
            printf("\tscilab         (turn on/off scilab (limited support))\n");
            printf("\tstart_matlab   (start matlab in PSUADE interactive mode)\n");
            printf("\tcheckformat    (check format of various files used by PSUADE)\n");
         }
         else if (!strcmp(winput, "advanced"))
         {
            printf("Advanced analysis and control commands:\n");
            printf("\tio_expert       (turn on/off IO expert mode)\n");
            printf("\trs_expert       (turn on/off response surface expert mode)\n");
            printf("\trs_codegen      (turn on/off response surface code generator)\n");
            printf("\tana_expert      (turn on/off analysis expert mode)\n");
            printf("\tsam_expert      (turn on/off sampling expert mode)\n");
            printf("\topt_expert      (turn on/off optimization expert mode)\n");
            printf("\tgenhistogram    (create a histogram from the loaded sample)\n");
            printf("\tgenhistogram2   (similar to genhistogram with a target nbins)\n");
            printf("\tgenconfigfile   (create a template PSUADE config file)\n");
            printf("\tuse_configfile <file> (use config file to select parameters)\n");
            printf("\trefine          (refine a sample: higher mesh resolution)\n");
            printf("\ta_refine        (adaptively refine a sample)\n");
            printf("\ta_refine_metis  (a_refine but with initial sample = METIS)\n");
            printf("\tinterface_track (track the interface separating Y=0 vs Y=1.)\n");
            printf("\tmoat_adjust  <file> (modify MOAT sample from <file>)\n");
            printf("\tgmoat_adjust <file> (modify GMOAT sample from <file>)\n");
            printf("\tmoatgen             (create MOAT adjust file: 1 constr)\n");
            printf("\tmoatgen2            (moatgen with multiple constraints)\n");
            printf("\tmoat_concat  <file> (combine 2 MOAT samples (2 INPUT SETS))\n");
         }
         else if (!strcmp(winput, "future"))
         {
            printf("\tgen_discrete (generate a sample of discrete variables)\n");
            printf("\tgp_sa2       (parameter screening via layered GP, not ready)\n");
            printf("\tsot_sa2      (parameter screening via sum-of-trees+boosting)\n");
            printf("\tgd_test      (Gower/Mahalanobis extrapolation analysis)\n");
            printf("\tpdfcheck     (internal self-check for accuracy of PDFs)\n");
         }
         else
         {
            printf("Help topics:\n");
            printf("\tinfo         (information about the use of PSUADE)\n");
            printf("\tio           (file read/write commands)\n");
            printf("\tstat         (basic statistics)\n");
            printf("\tscreen       (parameter screening commands)\n");
            printf("\trs           (response surface analysis commands)\n");
            printf("\tqsa          (quantitative sensitivity analysis commands)\n");
            printf("\tcalibration  (Bayesian calibration/optimization commands)\n");
            printf("\tplot         (commands for generating visualization plots)\n");
            printf("\tsetup        (commands to set up PSUADE work flow)\n");
            printf("\tedit         (commands to edit sample data in loal memory)\n");
            printf("\tmisc         (miscellaneous commands)\n");
            printf("\tadvanced     (advanced analysis and control commands)\n");
            printf("\t<command -h> (help for a specific command)\n");
         }
      }

      // +++ rename a file
      else if (!strcmp(command, "rename"))
      {
         sscanf(lineIn,"%s %s %s",command,winput,cString);
         if (!strcmp(winput, "-h"))
         {
            printf("syntax: rename <file1> <file2>\n");
            printf("where <file1> and <file2> are source and destination files.\n");
            printf("\nThis command allows you to change file name without\n");
            printf("     leaving the PSUADE command line mode.\n");
            continue;
         }
         fp = fopen(winput,"r");
         if (fp == NULL) printf("ERROR: file %s not found.\n", winput);
         else
         {
            fclose(fp);
            if (rename(winput, cString) == 0) 
                 printf("%s has been renamed to %s\n", winput, cString);
            else printf("ERROR: in renaming %s to %s\n", winput, cString);
         }
      }

      // +++ run
      else if (!strcmp(command, "run"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("syntax: run <file>\n");
            printf("where <file> is a PSUADE input file.\n");
            printf("\nThis command allows you to run PSUADE batch file in\n");
            printf("     command line mode.\n");
            continue;
         }
         strcpy(dataFile, "\0");
         sscanf(lineIn,"%s %s",command,dataFile);
         if ((fp=fopen(dataFile,"r")) == NULL)
         {
            printf("file %s not found.\n", dataFile);
            printf("syntax: run <file>.\n");
            printf("where <file> is a PSUADE input file.\n");
            continue;
         }
         else fclose(fp);
         status = getInputFromFile(dataFile);
         if (status == 0) run();
         else
         {
            printf("ERROR: batch file not valid : check file format.\n");
         }
      }

      // Input/output commands
      // +++ load and loadp
      else if (!strcmp(command, "load") || !strcmp(command, "loadp"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         { 
            printf("load: load PSUADE data from a file to local memory\n");
            printf("syntax: load <filename>.\n");
            printf("where <filename> is a PSUADE data file.\n");
            continue;
         }
         strcpy(dataFile, "\0");
         if (!strcmp(command,"load")) sscanf(lineIn,"%s %s",command,dataFile);
         else                         strcpy(dataFile, "psuadeIO");
         if (!strcmp(dataFile,psOutputFilename_))
         {
            printf("WARNING: you are loading a file with the same\n");
            printf("name as the default output file %s", psOutputFilename_);
            printf("This file may be overwritten during your session.\n");
            printf("You can change the default ouptut filename with\n");
            printf("the 'output_file' command.\n");
         }
         if ((fp=fopen(dataFile,"r")) == NULL)
         {
            printf("file %s not found.\n", dataFile);
            printf("syntax: load <file>.\n");
            printf("where <file> is a PSUADE data file.\n");
            continue;
         }
         else fclose(fp);
         cleanUp();
         if (psConfig_ != NULL) delete psConfig_;
         psConfig_ = new PsuadeConfig();
         psuadeIO_ = new PsuadeData();
         psuadeIO_->setOutputLevel(0);
         status = psuadeIO_->readPsuadeFile(dataFile);
         if (status == 0)
         {
            psuadeIO_->getParameter("input_ninputs", pPtr);
            nInputs_ = pPtr.intData_;
            pINames.clean();
            psuadeIO_->getParameter("input_names", pINames);
            names = pINames.strArray_;
            pLower.clean();
            psuadeIO_->getParameter("input_lbounds", pLower);
            iLowerB_ = pLower.dbleArray_;
            pLower.dbleArray_ = NULL;
            pUpper.clean();
            psuadeIO_->getParameter("input_ubounds", pUpper);
            iUpperB_ = pUpper.dbleArray_;
            pUpper.dbleArray_ = NULL;
            inputNames_ = new char*[nInputs_+1];
            for (ii = 0; ii < nInputs_; ii++)
            {
               inputNames_[ii] = new char[200]; 
               strcpy(inputNames_[ii], names[ii]);
            }
            psuadeIO_->getParameter("input_pdfs", pPDFs);
            psuadeIO_->getParameter("input_means", pMeans);
            psuadeIO_->getParameter("input_stdevs", pStds);
            inputPDFs_  = pPDFs.intArray_;
            pPDFs.intArray_ = NULL;
            inputMeans_ = pMeans.dbleArray_;
            pMeans.dbleArray_ = NULL;
            inputStds_  = pStds.dbleArray_;
            pStds.dbleArray_ = NULL;
            psuadeIO_->getParameter("input_sample_files", pPDFs);
            inputPDFFiles = pPDFs.strArray_;
            psuadeIO_->getParameter("input_sample_indices", pPDFSIndices);
            SPDFIndices = pPDFSIndices.intArray_;
            psuadeIO_->getParameter("input_cor_matrix", pPtr);
            inputCMat_ = new Matrix();
            Matrix *tmpMat = (Matrix *) pPtr.psObject_;
            inputCMat_->load(*tmpMat);

            psuadeIO_->getParameter("output_noutputs", pPtr);
            nOutputs_ = pPtr.intData_;
            pONames.clean();
            psuadeIO_->getParameter("output_names", pONames);
            names = pONames.strArray_;
            outputNames_ = new char*[nOutputs_+1];
            for (ii = 0; ii < nOutputs_; ii++)
            {
               outputNames_[ii] = new char[200]; 
               strcpy(outputNames_[ii], names[ii]);
            }
            psuadeIO_->getParameter("method_sampling", pPtr);
            samplingMethod = pPtr.intData_;
            psuadeIO_->getParameter("method_nsamples", pPtr);
            nSamples_ = pPtr.intData_;
            psuadeIO_->getParameter("method_nreplications",pPtr);
            nReps = pPtr.intData_;
            psuadeIO_->getParameter("input_sample", pPtr);
            sampleInputs_ = pPtr.dbleArray_;
            pPtr.dbleArray_ = NULL;
            psuadeIO_->getParameter("output_sample", pPtr);
            sampleOutputs_  = pPtr.dbleArray_;
            pPtr.dbleArray_ = NULL;
            psuadeIO_->getParameter("output_states", pPtr);
            sampleStates_  = pPtr.intArray_;
            pPtr.intArray_ = NULL;
            pINames.clean();
            pONames.clean();
            if (sampleInputs_ == NULL || sampleOutputs_ == NULL)
            {
               printf("WARNING: no sample matrix or output found.\n");
               nSamples_ = 0;
            }
            else
            {
               printf("load complete : nSamples = %d\n", nSamples_);
               printf("                nInputs  = %d\n", nInputs_);
               printf("                nOutputs = %d\n", nOutputs_);
               if (currSession != NULL) delete currSession;
               currSession = new PsuadeSession();
               psuadeIO_->getSession(currSession);
            }
            currSession->psuadeIO_ = psuadeIO_;
         }
         else
         {
            printf("ERROR: file %s either not found or in wrong format.\n",
                   dataFile);
            cleanUp();
         }
         fflush(stdout);
      }

      // +++ loadmore 
      else if (!strcmp(command, "loadmore")) 
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         { 
            printf("loadmore: add another PSUADE sample to the loaded sample\n");
            printf("          so it is equivalent to sample concatenation.\n");
            printf("syntax: loadmore <filename>.\n");
            printf("where <filename> is a PSUADE data file.\n");
            continue;
         }
         strcpy(dataFile, "\0");
         sscanf(lineIn, "%s %s", command, dataFile);
         if (!strcmp(dataFile,psOutputFilename_))
         {
            printf("WARNING: you are loading a file with the same\n");
            printf("name as the default output file %s", psOutputFilename_);
            printf("This file may be overwritten during your session.\n");
            printf("You can change the default output filename with\n");
            printf("the 'output_file' command.\n");
         }
         if ((fp=fopen(dataFile,"r")) == NULL)
         {
            printf("file %s not found.\n", dataFile);
            printf("syntax: loadmore <filename>.\n");
            printf("where <filename> is a PSUADE data file.\n");
            continue;
         }
         psuadeIO_ = new PsuadeData();
         psuadeIO_->setOutputLevel(0);
         status = psuadeIO_->readPsuadeFile(dataFile);
         if (status == 0)
         {
            flag = 0;
            if (sampleInputs_ != NULL)
            {
               psuadeIO_->getParameter("input_ninputs", pPtr);
               ind = pPtr.intData_;
               pINames.clean();
               psuadeIO_->getParameter("input_names", pINames);
               names = pINames.strArray_;
               pLower.clean();
               psuadeIO_->getParameter("input_lbounds", pLower);
               tempX = pLower.dbleArray_;
               pUpper.clean();
               psuadeIO_->getParameter("input_ubounds", pUpper);
               tempY = pUpper.dbleArray_;
               flag = 1;
               if (ind != nInputs_)
               {
                  printf("ERROR: nInputs are different.\n");
                  printf("       incoming nInputs = %d\n", ind);
                  printf("       expected nInputs = %d\n", nInputs_);
                  flag = 0;
               }
               if (flag == 1)
               {
                  psuadeIO_->getParameter("output_noutputs",pPtr);
                  ind = pPtr.intData_;
                  pONames.clean();
                  psuadeIO_->getParameter("output_names",pONames);
                  names = pONames.strArray_;
                  if (ind != nOutputs_)
                  {
                     printf("ERROR: nOutputs are different.\n");
                     printf("       incoming nOutputs = %d\n", ind);
                     printf("       expected nOutputs = %d\n", nOutputs_);
                     printf("INFO: local data set not changed.\n");
                     flag = 0;
                  }
               }
               if (flag == 1)
               {
                  for (ii = 0; ii < ind; ii++)
                  {
                     if (strcmp(outputNames_[ii],names[ii]))
                     {
                        sprintf(pString,
                                "Output names are different. Override? (y/n) ");
                        getString(pString, winput);
                        if (winput[0] != 'y')
                        {
                           printf("INFO: No override, data set not changed.\n");
                           flag = 0;
                           break;
                        }
                     }
                  }
               }
               if (flag == 1)
               {
                  for (ii = 0; ii < nInputs_; ii++)
                  {
                     if (tempX[ii] < iLowerB_[ii]) iLowerB_[ii] = tempX[ii];
                     if (tempY[ii] > iUpperB_[ii]) iUpperB_[ii] = tempY[ii];
                  }
               }
               tempX = NULL;
               tempY = NULL;
               names = NULL;
               pLower.clean();
               pUpper.clean();
               pINames.clean();
               pONames.clean();
            }
            else
            {
               psuadeIO_->getParameter("input_ninputs", pPtr);
               nInputs_ = pPtr.intData_;
               pINames.clean();
               psuadeIO_->getParameter("input_names", pINames);
               names = pINames.strArray_;
               pLower.clean();
               psuadeIO_->getParameter("input_lbounds", pLower);
               iLowerB_ = pLower.dbleArray_;
               pLower.dbleArray_ = NULL;
               pUpper.clean();
               psuadeIO_->getParameter("input_ubounds", pUpper);
               iUpperB_ = pUpper.dbleArray_;
               pUpper.dbleArray_ = NULL;
               inputNames_ = new char*[nInputs_+1];
               for (ii = 0; ii < nInputs_; ii++)
               {
                  inputNames_[ii] = new char[200]; 
                  strcpy(inputNames_[ii], names[ii]);
               }
               pINames.clean();
               psuadeIO_->getParameter("output_noutputs", pPtr);
               nOutputs_ = pPtr.intData_;
               pONames.clean();
               psuadeIO_->getParameter("output_names", pONames);
               names = pONames.strArray_;
               outputNames_ = new char*[nOutputs_+1];
               for (ii = 0; ii < nOutputs_; ii++)
               {
                  outputNames_[ii] = new char[200]; 
                  strcpy(outputNames_[ii], names[ii]);
               }
               pONames.clean();
               psuadeIO_->getParameter("method_sampling", pPtr);
               samplingMethod = pPtr.intData_;
               psuadeIO_->getParameter("method_nsamples", pPtr);
               nSamples_ = pPtr.intData_;
               psuadeIO_->getParameter("method_nreplications",pPtr);
               nReps = pPtr.intData_;
               psuadeIO_->getParameter("input_sample", pPtr);
               sampleInputs_  = pPtr.dbleArray_;
               pPtr.dbleArray_ = NULL;
               psuadeIO_->getParameter("output_sample", pPtr);
               sampleOutputs_  = pPtr.dbleArray_;
               pPtr.dbleArray_ = NULL;
               psuadeIO_->getParameter("output_states", pPtr);
               sampleStates_  = pPtr.intArray_;
               pPtr.intArray_ = NULL;
               if (sampleOutputs_ == NULL) nSamples_ = 0;
               if (currSession != NULL) delete currSession;
               currSession = new PsuadeSession();
               psuadeIO_->getSession(currSession);
               printf("loadmore complete : nSamples = %d\n", nSamples_);
               printf("                    nInputs  = %d\n", nInputs_);
               printf("                    nOutputs = %d\n", nOutputs_);
            }
            if (flag == 1)
            {
               psuadeIO_->getParameter("method_nsamples", pPtr);
               ind = pPtr.intData_;
               tempX  = sampleInputs_;
               tempY  = sampleOutputs_;
               states = sampleStates_;
               sampleInputs_  = new double[(nSamples_+ind)*nInputs_];
               sampleOutputs_ = new double[(nSamples_+ind)*nOutputs_];
               sampleStates_  = new int[nSamples_+ind];
               for (ii = 0; ii < nSamples_*nInputs_; ii++)
                  sampleInputs_[ii] = tempX[ii];
               for (ii = 0; ii < nSamples_*nOutputs_; ii++)
                  sampleOutputs_[ii] = tempY[ii];
               for (ii = 0; ii < nSamples_; ii++)
                  sampleStates_[ii] = states[ii];
               delete [] tempX;
               delete [] tempY;
               delete [] states;
               psuadeIO_->getParameter("input_sample", pPtr);
               tempX  = pPtr.dbleArray_;
               pPtr.dbleArray_ = NULL;
               psuadeIO_->getParameter("output_sample", pPtr);
               tempY  = pPtr.dbleArray_;
               psuadeIO_->getParameter("output_states", pPtr);
               states  = pPtr.intArray_;
               pPtr.intArray_ = NULL;
               pPtr.dbleArray_ = NULL;
               for (ii = 0; ii < ind*nInputs_; ii++)
                  sampleInputs_[nSamples_*nInputs_+ii] = tempX[ii];
               for (ii = 0; ii < ind*nOutputs_; ii++)
                  sampleOutputs_[nSamples_*nOutputs_+ii] = tempY[ii];
               for (ii = 0; ii < ind; ii++)
                  sampleStates_[nSamples_+ii] = states[ii];
               nSamples_ += ind;
               psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,
                              NULL,NULL,sampleInputs_,NULL,NULL,
                              NULL,NULL,NULL); 
               psuadeIO_->updateOutputSection(nSamples_,nOutputs_,
                                sampleOutputs_, sampleStates_, NULL); 
               if (currSession != NULL) delete currSession;
               currSession = new PsuadeSession();
               psuadeIO_->getSession(currSession);
               printf("loadmore complete : nSamples = %d\n", nSamples_);
               printf("                    nInputs  = %d\n", nInputs_);
               printf("                    nOutputs = %d\n", nOutputs_);
            }
         }
         else
         {
            printf("ERROR: file %s either not found or in wrong format.\n",
                   dataFile);
         }
         fflush(stdout);
         if (dataReg_ != NULL) delete [] dataReg_;
         dataReg_ = NULL;
      }

      // +++ write 
      else if (!strcmp(command, "write"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         { 
            printf("write: save the resident sample to another PSUADE data\n");
            printf("       file.\n");
            printf("syntax: write <filename>.\n");
            printf("where <filename> is a PSUADE data file.\n");
            continue;
         }
         strcpy(dataFile, "\0");
         if (nInputs_ <= 0 || psuadeIO_ == NULL)
         {
            printf("ERROR: data not loaded - nothing to write\n");
            continue;
         }
         strcpy(dataFile, psOutputFilename_);
         sscanf(lineIn, "%s %s", command, dataFile);
         if (!strcmp(dataFile,psInputFilename_))
         {
           printf("WARNING : output file name should not be the same as\n"); 
           printf(" the input file name %s.\n", psInputFilename_);
           printf(" Command not executed.\n");
         }
         else if (!strcmp(dataFile, "\0"))
         {
            printf("ERROR: no file to write to.\n");
            printf("syntax: write <filename>.\n");
            printf("where <filename> is a PSUADE data file.\n");
         }
         else
         {
            if (nOutputs_ > 1)
            {
               sprintf(pString,
                       "Save one output only (all otherwise)? (y or n) "); 
               getString(pString, winput);
            }
            else
            {
               winput[0] = 'n';
               outputID = 0;
            }
            if (winput[0] == 'y')
            {
               sprintf(pString,"Enter output number (1 - %d) : ",nOutputs_);
               outputID = getInt(1, nOutputs_, pString);
               outputID--;
               for (sInd = 0; sInd < nSamples_; sInd++)
                  sampleOutputs_[sInd] = 
                        sampleOutputs_[sInd*nOutputs_+outputID];
               pONames.clean();
               psuadeIO_->getParameter("output_names", pONames);
               names = pONames.strArray_;
               strcpy(names[0], names[outputID]);
               nOutputs_ = 1;
               psuadeIO_->updateAnalysisSection(0, 0, 0, 0, 0, 0);
            }
            else names = NULL;
            psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,iLowerB_,
                              iUpperB_,sampleInputs_,inputNames_,inputPDFs_,
                              inputMeans_,inputStds_,inputCMat_); 
            psuadeIO_->updateOutputSection(nSamples_,nOutputs_,
                             sampleOutputs_, sampleStates_, names); 
            if (winput[0] == 'y') pONames.clean();
            psuadeIO_->updateMethodSection(samplingMethod,nSamples_,nReps,
                                           -1,-1);
            psuadeIO_->writePsuadeFile(dataFile,0);
         }
      }

      // +++ nwrite 
      else if (!strcmp(command, "nwrite"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("nwrite: save the unevaluated points to a PSUADE data file\n");
            printf("  Note: sample points with undefined outputs or status=0\n");
            printf("        are considered unevaluated points.\n");
            printf("syntax: nwrite <filename>.\n");
            printf("where <filename> is a PSUADE data file.\n");
            continue;
         }

         if (nInputs_ <= 0 || psuadeIO_ == NULL)
         {
            printf("ERROR: data not loaded - nothing to write\n");
            continue;
         }
         strcpy(dataFile, psOutputFilename_);
         sscanf(lineIn, "%s %s", command, dataFile);
         if (!strcmp(dataFile,psInputFilename_))
         {
            printf("WARNING : output filename should not be the same as\n"); 
            printf(" the input file %s.\n", psInputFilename_);
            printf(" Command not executed.\n");
         }
         else if (!strcmp(dataFile, "\0"))
         {
            printf("ERROR: need to specify a file to write to.\n");
            printf("syntax: nwrite <filename>.\n");
         }
         else
         {
            tempX  = new double[nInputs_*nSamples_];
            tempY  = new double[nOutputs_*nSamples_];
            states = new int[nSamples_];
            count = 0;
            for (ii = 0; ii < nSamples_; ii++)
            {
               for (jj = 0; jj < nOutputs_; jj++)
                  if (sampleOutputs_[ii*nOutputs_+jj] == PSUADE_UNDEFINED)
                     sampleStates_[ii] = 0;
               if (sampleStates_[ii] == 0)
               {
                  for (jj = 0; jj < nInputs_; jj++)
                     tempX[count*nInputs_+jj] = 
                        sampleInputs_[ii*nInputs_+jj];
                  for (jj = 0; jj < nOutputs_; jj++)
                     tempY[count*nOutputs_+jj] = 
                        sampleOutputs_[ii*nOutputs_+jj]; 
                  states[count++] = sampleStates_[ii]; 
               }
            }
            if (count == 0)
            {
               printf("INFO: no unevaluated sample points\n");
               printf("      ==> no file generated.\n");
            }
            else
            {
               ioPtr = new PsuadeData();
               ioPtr->updateInputSection(count, nInputs_, NULL, iLowerB_,
                                   iUpperB_, tempX, inputNames_,
                                   NULL, NULL,NULL,NULL); 
               ioPtr->updateOutputSection(count, nOutputs_, tempY, states,
                                          outputNames_);
               ioPtr->updateMethodSection(PSUADE_SAMP_MC, count, 1, -1, -1);
               ioPtr->writePsuadeFile(dataFile, 0);
               delete ioPtr;
               ioPtr = NULL;
               printf("Data written to file %s.\n", dataFile);
            }
            delete [] tempX;
            delete [] tempY;
            delete [] states;
         }
      }

      // +++ iread 
      else if (!strcmp(command, "iread"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("iread: read sample inputs from a file written using\n");
            printf("       iwrite or a file with the following format:\n");
            printf("line 1: PSUADE_BEGIN\n");
            printf("line 2: <nSamples> <nInputs>\n");
            printf("line 3: (optional) '#' as first character, then input names\n");
            printf("line 4: 1 <sample point 1 inputs> \n");
            printf("line 5: 2 <sample point 2 inputs> \n");
            printf(".......\n");
            printf("line n: PSUADE_END\n\n");
            printf("syntax: iread <filename>.\n");
            printf("where <filename> is the name of the data file:\n");
            printf("\n Note: iread can be used to read a posterior sample\n");
            printf("         generated by rsmcmc.\n");
            continue;
         }
         strcpy(dataFile, "\0");
         sscanf(lineIn, "%s %s", command, dataFile);
         if ((fp = fopen(dataFile,"r")) == NULL)
         {
            printf("file %s not found.\n", dataFile);
            printf("syntax: iread <filename>.\n");
            printf("where <filename> is the name of the data file:\n");
            continue;
         }
         fclose(fp);
         cleanUp();
         psuadeIO_ = new PsuadeData();
         psuadeIO_->setOutputLevel(0);
         fp = fopen(dataFile,"r");
         fscanf(fp, "%s", winput);
         if (strcmp(winput, "PSUADE_BEGIN"))
         {
            fclose(fp);
            fp = fopen(dataFile,"r");
         }
         else fgets(lineIn, 5000, fp);
         fscanf(fp, "%d %d", &nSamples_, &nInputs_);
         nOutputs_ = 1;
         if (nSamples_ <= 0 || nInputs_ <= 0)
         {
            printf("iread ERROR: some sample parameters <= 0.\n");
            printf("      nSamples read = %d\n", nSamples_);
            printf("      nInputs  read = %d\n", nInputs_);
            printf("Note: the format for the second line should be: \n");
            printf("  <nSamples> <nInputs> \n");
            delete psuadeIO_;
            psuadeIO_ = NULL;
            fclose(fp);
            continue;
         }
         fgets(pString, 5000, fp);
         while (1)
         {
            kk = getc(fp);
            if (kk == '#')
            {
               inputNames_ = new char*[nInputs_];
               for (ii = 0; ii < nInputs_; ii++)
               {
                  inputNames_[ii] = new char[1001];
                  fscanf(fp,"%s", inputNames_[ii]);
               }
               fgets(pString, 5000, fp);
            }
            else
            {
               ungetc(kk, fp);
               break;
            }
         }
         sampleInputs_  = new double[nSamples_*nInputs_];
         sampleOutputs_ = new double[nSamples_*nOutputs_];
         sampleStates_  = new int[nSamples_];
         for (ii = 0; ii < nSamples_; ii++) sampleStates_[ii] = 1;
         for (ii = 0; ii < nSamples_; ii++)
         {
            fscanf(fp, "%d", &jj);
            if ((ii+1) != jj)
            {
               printf("iread ERROR: sample index mismatch.\n");
               printf("         sample number read     = %d\n", jj);
               printf("         sample number expected = %d\n", ii+1);
               printf("INFO: the first field must be the sample number.\n");
               cleanUp();
               break;
            }
            for (jj = 0; jj < nInputs_; jj++)
               fscanf(fp,"%lg", &sampleInputs_[ii*nInputs_+jj]);
            sampleOutputs_[ii] = PSUADE_UNDEFINED;
            sampleStates_[ii] = 0;
         }
         fscanf(fp, "%s", winput);
         fclose(fp);
         if (ii != nSamples_) continue;
         if (inputNames_ == NULL)
         {
            inputNames_ = new char*[nInputs_];
            for (ii = 0; ii < nInputs_; ii++)
            {
               inputNames_[ii] = new char[100];
               sprintf(inputNames_[ii], "X%d", ii+1);
            }
         }
         iLowerB_ = new double[nInputs_];
         iUpperB_ = new double[nInputs_];
         tempW = new double[nSamples_];
         double ireadMax, ireadMin;
         for (ii = 0; ii < nInputs_; ii++)
         {
            for (jj = 0; jj < nSamples_; jj++)
               tempW[jj] = sampleInputs_[jj*nInputs_+ii];
            ireadMax = ireadMin = tempW[0];
            for (jj = 1; jj < nSamples_; jj++)
            {
               if (tempW[jj] > ireadMax) ireadMax = tempW[jj];
               if (tempW[jj] < ireadMin) ireadMin = tempW[jj];
            } 
            iLowerB_[ii] = ireadMin;
            iUpperB_[ii] = ireadMax;
            if (iLowerB_[ii] == iUpperB_[ii]) iUpperB_[ii] += 1.0e-12;
         }
         delete [] tempW;
         tempW = NULL;
         inputPDFs_  = new int[nInputs_];
         inputMeans_ = new double[nInputs_];
         inputStds_  = new double[nInputs_];
         inputCMat_  = new Matrix();
         inputCMat_->setDim(nInputs_, nInputs_);
         for (ii = 0; ii < nInputs_; ii++)
         {
            inputPDFs_[ii] = 0;
            inputMeans_[ii] = 0;
            inputStds_[ii] = 0;
            inputCMat_->setEntry(ii,ii,1.0);
         }
         psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,iLowerB_,iUpperB_,
                                   sampleInputs_,inputNames_,inputPDFs_,inputMeans_,
                                   inputStds_,inputCMat_); 

         outputNames_ = new char*[nOutputs_];
         for (ii = 0; ii < nOutputs_; ii++)
         {
            outputNames_[ii] = new char[100];
            sprintf(outputNames_[ii], "Y%d", ii+1);
         }
         psuadeIO_->updateOutputSection(nSamples_, nOutputs_, sampleOutputs_,
                                        sampleStates_, outputNames_);
         psuadeIO_->updateMethodSection(PSUADE_SAMP_MC, nSamples_, 1, 0, 0);
         psuadeIO_->updateAnalysisSection(0, 0, 0, 0, 0, 0);
         nReps = 1;
         if (currSession != NULL) delete currSession;
         currSession = new PsuadeSession();
         psuadeIO_->getSession(currSession);
         printf("iread: data have been read.\n");
         printf("nSamples = %d\n", nSamples_);
         printf("nInputs  = %d\n", nInputs_);
         printf("nOutputs = %d (set to 1)\n", nOutputs_);
      }

      // +++ iwrite 
      else if (!strcmp(command, "iwrite"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("iwrite: save only the sample inputs to an ASCII file\n");
            printf("syntax: iwrite <filename>.\n");
            printf("where <filename> is the name of the target data file.\n");
            printf("\nThe target file will have the following format: \n");
            printf("line 1: PSUADE_BEGIN\n");
            printf("line 2: <nSamples> <nInputs>\n");
            printf("line 3: # input parameter names\n");
            printf("line 4: 1 <sample point 1 inputs> \n");
            printf("line 5: 2 <sample point 2 inputs> \n");
            printf(".......\n");
            printf("line n: PSUADE_END\n\n");
            continue;
         }
         if (psuadeIO_ == NULL || sampleOutputs_ == NULL)
         {
            printf("ERROR: no sample output data.\n");
            printf("       Use load first to put data into local memory.\n");
            continue;
         }
         strcpy(dataFile, psOutputFilename_);
         sscanf(lineIn, "%s %s", command, dataFile);
         if (!strcmp(dataFile,psInputFilename_))
         {
           printf("WARNING : output file name should not be the same as\n"); 
           printf(" the input file name %s. Try the output_file command.\n", 
                  psInputFilename_);
           printf(" Command not executed.\n");
         }
         else if (!strcmp(dataFile, "\0"))
         {
            printf("ERROR: need to specify a file to write to.\n");
         }
         else
         {
            if (sampleInputs_ == NULL)
               printf("ERROR: no input to write.\n");
            else
            {
               fp = fopen(dataFile, "w");
               if (fp == NULL)
               {
                  printf("ERROR: cannot open file %s.\n", dataFile);
                  continue;
               }
               fprintf(fp, "PSUADE_BEGIN\n");
               fprintf(fp, "%d %d\n", nSamples_, nInputs_);
               fprintf(fp, "# ");
               for (ii = 0; ii < nInputs_; ii++)
                  fprintf(fp, "%s ", inputNames_[ii]);
               fprintf(fp, "\n");
               for (ii = 0; ii < nSamples_; ii++)
               {
                  fprintf(fp, "%d ", ii+1);
                  for (jj = 0; jj < nInputs_; jj++)
                     fprintf(fp,"%24.16e ",sampleInputs_[ii*nInputs_+jj]); 
                  fprintf(fp, "\n");
               }
               fprintf(fp, "PSUADE_END\n");
               fclose(fp);
               printf("iwrite: sample inputs written to %s.\n",dataFile);
            }
         }
      }
      
      // +++ owrite 
      else if (!strcmp(command, "owrite"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("owrite: save only the sample outputs to an ASCII file\n");
            printf("syntax: owrite <filename>.\n");
            printf("where <filename> is the name of the target data file.\n");
            printf("\nThe target file will have the following format: \n");
            printf("line 1: PSUADE_BEGIN\n");
            printf("line 2: <nSamples> <nOutputs>\n");
            printf("line 3: 1 <sample point 1 outputs> \n");
            printf("line 4: 2 <sample point 2 outputs> \n");
            printf(".......\n");
            printf("line n: PSUADE_END\n\n");
            continue;
         }
         if (psuadeIO_ == NULL || sampleOutputs_ == NULL)
         {
            printf("ERROR: no sample output data.\n");
            printf("       Use load first to put data into local memory.\n");
            continue;
         }
         strcpy(dataFile, "\0");
         sscanf(lineIn, "%s %s", command, dataFile);
         strcpy(dataFile, psOutputFilename_);
         sscanf(lineIn, "%s %s", command, dataFile);
         if (!strcmp(dataFile,psInputFilename_))
         {
            printf("WARNING : output file name should not be the same as\n"); 
            printf(" the input file name %s. Try the output_file command.\n", 
                   psInputFilename_);
            printf(" Command not executed.\n");
         }
         else if (!strcmp(dataFile, "\0"))
         {
            printf("ERROR: need to specify a file to write to.\n");
         }
         else
         {
            if (sampleOutputs_ == NULL)
               printf("ERROR: no output to write.\n");
            else
            {
               if (nOutputs_ == 1) outputID = 0;
               else
               {
                  sprintf(pString, "Output number (1 - %d, %d for all): ",
                          nOutputs_, nOutputs_+1);
                  outputID = getInt(1, nOutputs_+1, pString);
                  outputID--;
               }
               fp = fopen(dataFile, "w");
               if (fp == NULL)
               {
                  printf("ERROR: cannot open file %s.\n", dataFile);
                  continue;
               }
               fprintf(fp, "PSUADE_BEGIN\n");
               if (outputID == nOutputs_)
                    fprintf(fp, "%d %d\n", nSamples_, nOutputs_);
               else fprintf(fp, "%d 1\n", nSamples_);
               for (ii = 0; ii < nSamples_; ii++)
               {
                  fprintf(fp, "%d ", ii+1);
                  if (outputID == nOutputs_)
                  {
                     for (jj = 0; jj < nOutputs_; jj++)
                        fprintf(fp,"%24.16e ",sampleOutputs_[ii*nOutputs_+jj]); 
                     fprintf(fp, "\n");
                  }
                  else
                     fprintf(fp,"%24.16e\n",
                             sampleOutputs_[ii*nOutputs_+outputID]);
               }
               fprintf(fp, "PSUADE_END\n");
               fclose(fp);
               printf("owrite: sample outputs written to %s.\n",dataFile);
            }
         }
      }

      // +++ read_std 
      else if (!strcmp(command, "read_std"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("read_std: read a data file in the following format:\n");
            printf("line 1: <nSamples> <nInputs> <nOutputs>\n");
            printf("line 2: <sample 1 inputs> < sample 1 outputs>\n");
            printf("line 3: <sample 2 inputs> < sample 2 outputs>\n");
            printf(".......\n\n");
            printf("syntax: read_std <filename>.\n");
            printf("where <filename> is a data file in standard format.\n");
            continue;
         }
         strcpy(dataFile, "\0");
         sscanf(lineIn, "%s %s", command, dataFile);
         if (!strcmp(dataFile, "\0"))
         {
            printf("ERROR: please specify a file to read from.\n");
            printf("syntax: read_std <filename>.\n");
            printf("where <filename> is a data file in standard format.\n");
            continue;
         }
         if ((fp=fopen(dataFile,"r")) == NULL)
         {
            printf("file %s not found.\n", dataFile);
            printf("syntax: read_std <filename>.\n");
            printf("where <filename> is a data file in standard format.\n");
            continue;
         }
         cleanUp();
         fscanf(fp, "%d %d %d", &nSamples_, &nInputs_, &nOutputs_);
         if (nSamples_ <= 0) 
         {
            printf("ERROR: nSamples <= 0 (%d)\n", nSamples_);
            continue;
         }
         if (nInputs_ <= 0) 
         {
            printf("ERROR: nInputs <= 0 (%d)\n", nInputs_);
            continue;
         }
         flag = 0;
         if (nOutputs_ <= 0) 
         {
            printf("INFO: nOutputs = 0 ==> create a dummy output\n");
            nOutputs_ = 1;
            flag = 1;
         }
         printf("nSamples = %d\n", nSamples_);
         printf("nInputs  = %d\n", nInputs_);
         printf("nOutputs = %d\n", nOutputs_);
         sampleInputs_  = new double[nSamples_*nInputs_];
         sampleOutputs_ = new double[nSamples_*nOutputs_];
         sampleStates_  = new int[nSamples_];
         for (ii = 0; ii < nSamples_; ii++)
         {
            for (jj = 0; jj < nInputs_; jj++)
               fscanf(fp, "%lg", &sampleInputs_[ii*nInputs_+jj]); 
            if (flag == 0)
            {
               for (jj = 0; jj < nOutputs_; jj++)
                  fscanf(fp,"%lg", &sampleOutputs_[ii*nOutputs_+jj]);
            }
            else sampleOutputs_[ii] = PSUADE_UNDEFINED;
            sampleStates_[ii] = 1;
            for (jj = 0; jj < nOutputs_; jj++)
               if (sampleOutputs_[ii*nOutputs_+jj] > 0.5*PSUADE_UNDEFINED)
                  sampleStates_[ii] = 0;
         }
         fclose(fp);
         inputNames_ = new char*[nInputs_+1];
         for (ii = 0; ii < nInputs_; ii++)
         {
            inputNames_[ii] = new char[200]; 
            sprintf(inputNames_[ii], "X%d", ii+1);
         }
         outputNames_ = new char*[nOutputs_+1];
         for (ii = 0; ii < nOutputs_; ii++)
         {
            outputNames_[ii] = new char[200]; 
            sprintf(outputNames_[ii], "Y%d", ii+1);
         }
         if (iLowerB_ != NULL) delete [] iLowerB_;
         iLowerB_ = new double[nInputs_];
         if (iUpperB_ != NULL) delete [] iUpperB_;
         iUpperB_ = new double[nInputs_];
         for (jj = 0; jj < nInputs_; jj++)
         {
            iLowerB_[jj] = sampleInputs_[jj];
            iUpperB_[jj] = sampleInputs_[jj];
            for (ii = 1; ii < nSamples_; ii++)
            {
               if (sampleInputs_[ii*nInputs_+jj] < iLowerB_[jj])
                  iLowerB_[jj] = sampleInputs_[ii*nInputs_+jj];
               if (sampleInputs_[ii*nInputs_+jj] > iUpperB_[jj])
                  iUpperB_[jj] = sampleInputs_[ii*nInputs_+jj];
            }
         }
         samplingMethod = PSUADE_SAMP_MC;
         inputPDFs_  = new int[nInputs_];
         inputMeans_ = new double[nInputs_];
         inputStds_  = new double[nInputs_];
         inputCMat_  = new Matrix();
         inputCMat_->setDim(nInputs_, nInputs_);
         for (ii = 0; ii < nInputs_; ii++)
         {
            inputPDFs_[ii] = 0;
            inputMeans_[ii] = 0;
            inputStds_[ii] = 0;
            inputCMat_->setEntry(ii,ii,1.0);
         }
         psuadeIO_ = new PsuadeData();
         psuadeIO_->setOutputLevel(0);
         psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,iLowerB_,
                               iUpperB_,sampleInputs_,inputNames_,inputPDFs_,
                               inputMeans_,inputStds_,inputCMat_); 
         psuadeIO_->updateOutputSection(nSamples_,nOutputs_,sampleOutputs_, 
                                        sampleStates_, outputNames_); 
         nReps = 1;
         psuadeIO_->updateMethodSection(samplingMethod,nSamples_,nReps,0,0);
         if (currSession != NULL) delete currSession;
         currSession = new PsuadeSession();
         psuadeIO_->getSession(currSession);
      }

      // +++ write_std 
      else if (!strcmp(command, "write_std"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("write_std: write a data file in the standard format.\n");
            printf("           Use <read_std -h> to see the standard\n");
            printf("           format.\n");
            printf("syntax: write_std <filename>.\n");
            printf("where <filename> is the name of the target data file.\n");
            continue;
         }
         if (psuadeIO_ == NULL || sampleOutputs_ == NULL)
         {
            printf("ERROR: no sample output data.\n");
            printf("       Use load to put data into local memory first.\n");
            continue;
         }
         strcpy(dataFile, "\0");
         sscanf(lineIn, "%s %s", command, dataFile);
         if (!strcmp(dataFile,psInputFilename_))
         {
            printf("WARNING : output file name should not be the same as\n"); 
            printf(" the input file name %s.\n", psInputFilename_);
            printf(" Command not executed.\n");
         }
         else if (!strcmp(dataFile, "\0"))
         {
            printf("ERROR: need to specify a file to write to.\n");
            printf("syntax: write_std <filename>.\n");
         }
         else
         {
            if (sampleInputs_ == NULL) printf("ERROR: no inputs to write.\n");
            else
            {
               if (sampleOutputs_ == NULL) outputID = -9;
               else
               {
                  if (nOutputs_ == 1) outputID = 0;
                  else
                  {
                     sprintf(pString,"Enter output ID (1 - %d, 0 for all) : ",
                             nOutputs_);
                     outputID = getInt(0, nOutputs_, pString);
                     outputID--;
                  }
               }
               fp = fopen(dataFile, "w");
               if (fp == NULL)
               {
                  printf("ERROR: cannot open file %s.\n", dataFile);
                  continue;
               }
               if (outputID == -9) 
                    fprintf(fp, "%d %d\n", nSamples_, nInputs_);
               else if (outputID == -1) 
                    fprintf(fp, "%d %d %d\n", nSamples_, nInputs_, nOutputs_);
               else fprintf(fp, "%d %d 1\n", nSamples_, nInputs_);
               for (ii = 0; ii < nSamples_; ii++)
               {
                  for (jj = 0; jj < nInputs_; jj++)
                     fprintf(fp, "%24.16e ", sampleInputs_[ii*nInputs_+jj]); 
                  if (outputID == -9) fprintf(fp, "\n");
                  else if (outputID == -1)
                  {
                     for (jj = 0; jj < nOutputs_; jj++)
                        fprintf(fp,"%24.16e ",sampleOutputs_[ii*nOutputs_+jj]);
                     fprintf(fp,"\n");
                  }
                  else fprintf(fp,"%24.16e\n",
                               sampleOutputs_[ii*nOutputs_+outputID]);
               }
               fclose(fp);
            }
         }
      }

      // +++ read_csv 
      else if (!strcmp(command, "read_csv"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("read_csv: read a data file in CSV format.\n");
            printf("syntax: read_csv <filename>.\n");
            printf("where <filename> is in the following CSV format.\n");
            printf("line 1 : (optional) INPUT,INPUT,...,OUTPUT, ...)\n");
            printf("line 2 : (optional) variable name list)\n");
            printf("line 3 : val1, val2, ... valn\n");
            printf(".....\n");
            printf("line m : val1, val2, ... valn\n\n");
            printf("That is, line 1 uses 'INPUT' and 'OUTPUT' to indicate\n");
            printf("    which columns are inputs and which are outputs.\n");
            printf("    Line 2 specifies input and output variable names.\n");
            printf("    Line 3 has the input/output values of sample point 1.\n");
            printf("    Line 4 has the input/output values of sample point 2...\n");
            continue;
         }
         strcpy(dataFile, "\0");
         sscanf(lineIn, "%s %s", command, dataFile);
         if (!strcmp(dataFile, "\0"))
         {
            printf("ERROR: please specify a file to read from.\n");
            printf("syntax: read_csv <filename>.\n");
            printf("where <filename> is in the following CSV format.\n");
            printf("line 1 : (optional) INPUT,INPUT,...,OUTPUT, ...)\n");
            printf("line 2 : (optional) variable name list)\n");
            printf("line 3 : val1, val2, ... valn\n");
            printf(".....\n");
            printf("line m : val1, val2, ... valn\n\n");
            printf("That is, line 1 uses 'INPUT' and 'OUTPUT' to indicate\n");
            printf("    which columns are inputs and which are outputs.\n");
            printf("    Line 2 specifies input and output variable names.\n");
            printf("    Line 3 has the input/output values of sample point 1.\n");
            printf("    Line 4 has the input/output values of sample point 2...\n");
            continue;
         }
         if ((fp=fopen(dataFile,"r")) == NULL)
         {
            printf("file %s not found.\n", dataFile);
            printf("syntax: read_csv <filename>.\n");
            printf("where <filename> is in the following CSV format.\n");
            printf("line 1 : (optional) INPUT,INPUT,...,OUTPUT, ...)\n");
            printf("line 2 : (optional) variable name list)\n");
            printf("line 3 : val1, val2, ... valn\n");
            printf(".....\n");
            printf("line m : val1, val2, ... valn\n\n");
            printf("That is, line 1 uses 'INPUT' and 'OUTPUT' to indicate\n");
            printf("    which columns are inputs and which are outputs.\n");
            printf("    Line 2 specifies input and output variable names.\n");
            printf("    Line 3 has the input/output values of sample point 1.\n");
            printf("    Line 4 has the input/output values of sample point 2...\n");
            continue;
         }
         cleanUp();

         int  nInOuts;
         char **inames, **onames;
         status = read_csv(dataFile, &nSamples_, &nInputs_, &tempX, &nOutputs_, 
                           &inames, &onames);
         if (status < 0 || nSamples_ <= 0 || nInputs_ <= 0) 
         {
            printf("ERROR: not a valid csv file (%d,%d,%d).\n",status,
                   nSamples_,nInputs_);
            continue;
         }
         nInOuts = nInputs_ + nOutputs_;
         flag = 0;
         if (nOutputs_ == 0)
         {
            sprintf(pString,"Enter number of input variables (1 - %d) : ",nInputs_);
            kk = getInt(1, nInputs_, pString);
            if (nInputs_ == kk)
            {
               printf("INFO: nOutputs = 0 ==> create a dummy output\n");
               nOutputs_ = 1;
               flag = 1;
               nInOuts = nInputs_;
            } 
            else
            {
               nInOuts = nInputs_;
               nOutputs_ = nInputs_ - kk;
               nInputs_  = kk;
            }
         } 
         printf("nSamples = %d\n", nSamples_);
         printf("nInputs  = %d\n", nInputs_);
         printf("nOutputs = %d\n", nOutputs_);
         sampleInputs_  = new double[nSamples_*nInputs_];
         sampleOutputs_ = new double[nSamples_*nOutputs_];
         sampleStates_  = new int[nSamples_];
         for (ii = 0; ii < nSamples_; ii++)
         {
            for (jj = 0; jj < nInputs_; jj++)
               sampleInputs_[ii*nInputs_+jj] = tempX[ii*nInOuts+jj]; 
            if (flag == 0)
            {
               for (jj = 0; jj < nOutputs_; jj++)
                  sampleOutputs_[ii*nOutputs_+jj] = tempX[ii*nInOuts+nInputs_+jj];
            }
            else sampleOutputs_[ii] = PSUADE_UNDEFINED;
            sampleStates_[ii] = 1;
            for (jj = 0; jj < nOutputs_; jj++)
               if (sampleOutputs_[ii*nOutputs_+jj] > 0.5*PSUADE_UNDEFINED)
                  sampleStates_[ii] = 0;
         }
         free(tempX);
         inputNames_ = new char*[nInputs_+1];
         for (ii = 0; ii < nInputs_; ii++)
         {
            inputNames_[ii] = new char[1000]; 
            if      (inames == NULL) sprintf(inputNames_[ii], "X%d", ii+1);
            else if (inames[ii] == NULL) sprintf(inputNames_[ii], "X%d", ii+1);
            else
            {
               strcpy(inputNames_[ii], inames[ii]);
               free(inames[ii]);
            }
         }
         outputNames_ = new char*[nOutputs_+1];
         for (ii = 0; ii < nOutputs_; ii++)
         {
            outputNames_[ii] = new char[200]; 
            if      (onames == NULL) sprintf(outputNames_[ii], "Y%d", ii+1);
            else if (onames[ii] == NULL) sprintf(outputNames_[ii], "Y%d", ii+1);
            else
            {
               strcpy(outputNames_[ii], onames[ii]);
               free(onames[ii]);
            }
         }
         if (inames != NULL) free(inames);
         if (onames != NULL) free(onames);
         if (iLowerB_ != NULL) delete [] iLowerB_;
         iLowerB_ = new double[nInputs_];
         if (iUpperB_ != NULL) delete [] iUpperB_;
         iUpperB_ = new double[nInputs_];
         for (jj = 0; jj < nInputs_; jj++)
         {
            iLowerB_[jj] = sampleInputs_[jj];
            iUpperB_[jj] = sampleInputs_[jj];
            for (ii = 1; ii < nSamples_; ii++)
            {
               if (sampleInputs_[ii*nInputs_+jj] < iLowerB_[jj])
                  iLowerB_[jj] = sampleInputs_[ii*nInputs_+jj];
               if (sampleInputs_[ii*nInputs_+jj] > iUpperB_[jj])
                  iUpperB_[jj] = sampleInputs_[ii*nInputs_+jj];
            }
         }
         samplingMethod = PSUADE_SAMP_MC;
         inputPDFs_  = new int[nInputs_];
         inputMeans_ = new double[nInputs_];
         inputStds_  = new double[nInputs_];
         inputCMat_  = new Matrix();
         inputCMat_->setDim(nInputs_, nInputs_);
         for (ii = 0; ii < nInputs_; ii++)
         {
            inputPDFs_[ii] = 0;
            inputMeans_[ii] = 0;
            inputStds_[ii] = 0;
            inputCMat_->setEntry(ii,ii,1.0);
         }
         psuadeIO_ = new PsuadeData();
         psuadeIO_->setOutputLevel(0);
         psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,iLowerB_,
                               iUpperB_,sampleInputs_,inputNames_,inputPDFs_,
                               inputMeans_,inputStds_,inputCMat_); 
         psuadeIO_->updateOutputSection(nSamples_,nOutputs_,sampleOutputs_, 
                                        sampleStates_, outputNames_); 
         nReps = 1;
         psuadeIO_->updateMethodSection(samplingMethod,nSamples_,nReps,0,0);
         if (currSession != NULL) delete currSession;
         currSession = new PsuadeSession();
         psuadeIO_->getSession(currSession);
      }

      // +++ write_matlab 
      else if (!strcmp(command, "write_matlab"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("write_matlab: write a data file in the matlab format.\n");
            printf("syntax: write_matlab <filename>.\n");
            printf("where <filename> is the name of the target data file.\n\n");
            printf("Info: the target file will have 2 matrices. Matrix X has\n");
            printf("      the sample inputs and Y has the sample outputs.\n");
            continue;
         }
         if (psuadeIO_ == NULL || sampleOutputs_ == NULL)
         {
            printf("ERROR: no sample output data.\n");
            printf("       Use load to put data into local memory first.\n");
            continue;
         }
         strcpy(dataFile, "\0");
         sscanf(lineIn, "%s %s", command, dataFile);
         if (!strcmp(dataFile,psInputFilename_))
         {
            printf("WARNING : output file name should not be the same as\n"); 
            printf(" the input file name %s. \n", psInputFilename_);
            printf(" Command not executed.\n");
         }
         else if (!strcmp(dataFile, "\0"))
         {
            printf("ERROR: need to specify a file to write to.\n");
            printf("syntax: write_matlab <filename>.\n");
         }
         else
         {
            if (sampleInputs_ == NULL) printf("ERROR: no inputs to write.\n");
            else
            {
               if (sampleOutputs_ == NULL)
                  printf("INFO: no outputs to write (only inputs).\n");
               fp = fopen(dataFile, "w");
               if (fp == NULL)
               {
                  printf("ERROR: cannot open file %s.\n", dataFile);
                  continue;
               }
               fprintf(fp, "%% X - input matrix\n");
               fprintf(fp, "%% Y - output matrix\n");
               fprintf(fp, "X = [\n");
               for (ii = 0; ii < nSamples_; ii++)
               {
                  for (jj = 0; jj < nInputs_; jj++)
                     fprintf(fp, "%24.16e ", sampleInputs_[ii*nInputs_+jj]); 
                  fprintf(fp,"\n");
               }
               fprintf(fp,"];\n");
               if (sampleOutputs_ != NULL)
               {
                  fprintf(fp, "Y = [\n");
                  for (ii = 0; ii < nSamples_; ii++)
                  {
                     for (jj = 0; jj < nOutputs_; jj++)
                        fprintf(fp,"%24.16e ",sampleOutputs_[ii*nOutputs_+jj]);
                     fprintf(fp,"\n");
                  }
                  fprintf(fp,"];\n");
               }
               fclose(fp);
            }
         }
      }

      // +++ read_xls  
      else if (!strcmp(command, "read_xls"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("read_xls: read a data file in Excel format.\n");
            printf("syntax: read_xls <filename>.\n");
            printf("where <filename> is a data file in xls format:\n");
            printf("line 1: <nSamples> <nInputs> <nOutputs>\n");
            printf("line 2: (optional) input and output names\n");
            printf("line 3: 1 <sample 1 inputs> < sample 1 outputs>\n");
            printf("line 4: 2 <sample 2 inputs> < sample 2 outputs>\n");
            printf(".......\n");
            continue;
         }
         strcpy(dataFile, "\0");
         sscanf(lineIn,"%s %s",command,dataFile);
         if ((fp=fopen(dataFile,"r")) == NULL)
         {
            printf("file %s not found.\n", dataFile);
            printf("syntax: read_xls <filename>.\n");
            printf("where <filename> is a data file in xls format.\n");
            continue;
         }
         else fclose(fp);
         cleanUp();

         psuadeIO_ = new PsuadeData();
         psuadeIO_->setOutputLevel(0);
         fp = fopen(dataFile,"r");
         fscanf(fp, "%d %d %d", &nSamples_, &nInputs_, &nOutputs_);
         printf("nSamples = %d\n", nSamples_);
         printf("nInputs  = %d\n", nInputs_);
         printf("nOutputs = %d\n", nOutputs_);
         if (nSamples_ <= 0 || nInputs_ <= 0 || nOutputs_ < 0)
         {
            printf("read_xls ERROR: some sample parameters <= 0.\n");
            printf("Note: the first line should be: \n");
            printf("  <nSamples> <nInputs> <nOutputs>\n");
            fclose(fp);
            nSamples_ = nInputs_ = nOutputs_ = 0;
            continue;
         }
         flag = 0;
         if (nOutputs_ == 0)
         {
            flag = 1;
            nOutputs_ = 1;
         }
         sampleInputs_  = new double[nSamples_*nInputs_];
         sampleOutputs_ = new double[nSamples_*nOutputs_];
         sampleStates_  = new int[nSamples_];
         for (ii = 0; ii < nSamples_; ii++) sampleStates_[ii] = 1;
         iLowerB_ = new double[nInputs_];
         iUpperB_ = new double[nInputs_];
         fgets(pString, 50000, fp);
         fgets(pString, 50000, fp);
         int hasComment = 0;
         if (pString[0] == '#') hasComment = 1;
         else
         {
            ii = 0;
            while ((pString[ii] == ' '  || pString[ii] == '\t') && 
                    pString[ii] != '\r' && pString[ii] != '\n' && 
                   ii < 5000) ii++;
            if (pString[ii] == '-') ii++;
            kk = (int) pString[ii] - '0';
            if (kk != 1) hasComment = 1;
            if (kk == 1 && (pString[ii+1] == ' ' || pString[ii+1] == '\t'))
                 hasComment = 0;
            else hasComment = 1;
         }
         fclose(fp);
         fp = fopen(dataFile,"r");
         fscanf(fp, "%d %d %d", &nSamples_, &nInputs_, &nOutputs_);
         if (hasComment == 1)
         {
            inputNames_ = new char*[nInputs_];
            for (ii = 0; ii < nInputs_; ii++)
            {
               inputNames_[ii] = new char[1001];
               fscanf(fp,"%s", inputNames_[ii]);
            }
            outputNames_ = new char*[nOutputs_];
            if (flag == 1) strcpy(outputNames_[0], "Y");
            else
            {
               for (ii = 0; ii < nOutputs_; ii++)
               {
                  outputNames_[ii] = new char[1001];
                  fscanf(fp,"%s", outputNames_[ii]);
               }
            }
            fgets(pString, 5000, fp);
         }
         for (ii = 0; ii < nSamples_; ii++)
         {
            fscanf(fp, "%d", &jj);
            if ((ii+1) != jj)
            {
               printf("read_xls ERROR: sample index mismatch.\n");
               printf("         sample number read     = %d\n", jj);
               printf("         sample number expected = %d\n", ii+1);
               printf("INFO: the first field must be the sample number.\n");
               break;
            }
            for (jj = 0; jj < nInputs_; jj++)
               fscanf(fp,"%lg", &sampleInputs_[ii*nInputs_+jj]); 
            if (flag == 0)
            {
               for (jj = 0; jj < nOutputs_; jj++)
                  fscanf(fp,"%lg", &sampleOutputs_[ii*nOutputs_+jj]); 
            }
            else
            {
               sampleOutputs_[ii] = PSUADE_UNDEFINED;
               sampleStates_[ii] = 0;
            }
         }
         fclose(fp);
         if (ii != nSamples_)
         {
            delete [] sampleInputs_;
            delete [] sampleOutputs_;
            delete [] sampleStates_;
            delete [] iLowerB_;
            delete [] iUpperB_;
            sampleInputs_ = NULL;
            sampleOutputs_ = NULL;
            sampleStates_ = NULL;
            iLowerB_ = NULL;
            iUpperB_ = NULL;
            if (inputNames_ != NULL)
            {
               for (ii = 0; ii < nInputs_; ii++) 
                  if (inputNames_[ii] != NULL) delete [] inputNames_[ii];
               inputNames_ = NULL;
            }
            if (outputNames_ != NULL)
            {
               for (ii = 0; ii < nOutputs_; ii++) 
                  if (outputNames_[ii] != NULL) delete [] outputNames_[ii];
               outputNames_ = NULL;
            }
            nInputs_ = nSamples_ = nOutputs_ = 0;
            continue;
         }
         tempW = new double[nSamples_];
         double xlsMin, xlsMax;
         for (ii = 0; ii < nInputs_; ii++)
         {
            for (jj = 0; jj < nSamples_; jj++) 
               tempW[jj] = sampleInputs_[jj*nInputs_+ii];
            xlsMin = xlsMax = tempW[0];
            for (jj = 1; jj < nSamples_; jj++) 
            {
               if (tempW[jj] > xlsMax) xlsMax = tempW[jj];
               if (tempW[jj] < xlsMin) xlsMin = tempW[jj];
            }
            iLowerB_[ii] = xlsMin;
            iUpperB_[ii] = xlsMax;
            if (iLowerB_[ii] == iUpperB_[ii]) iUpperB_[ii] += 1.0e-12;
         }
         delete [] tempW;
         inputPDFs_  = new int[nInputs_];
         inputMeans_ = new double[nInputs_];
         inputStds_  = new double[nInputs_];
         inputCMat_  = new Matrix();
         inputCMat_->setDim(nInputs_, nInputs_);
         for (ii = 0; ii < nInputs_; ii++)
         {
            inputPDFs_[ii] = 0;
            inputMeans_[ii] = 0;
            inputStds_[ii] = 0;
            inputCMat_->setEntry(ii,ii,1.0);
         }
         if (inputNames_ == NULL)
         {
            inputNames_ = new char*[nInputs_];
            for (ii = 0; ii < nInputs_; ii++)
            {
               inputNames_[ii] = new char[100];
               sprintf(inputNames_[ii], "X%d", ii+1);
            }
         } 
         if (outputNames_ == NULL)
         {
            outputNames_ = new char*[nOutputs_];
            for (ii = 0; ii < nOutputs_; ii++)
            {
               outputNames_[ii] = new char[100];
               sprintf(outputNames_[ii], "Y%d", ii+1);
            }
         } 
         psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,iLowerB_,iUpperB_,
                               sampleInputs_,inputNames_,inputPDFs_,inputMeans_,
                               inputStds_,inputCMat_); 
         psuadeIO_->updateOutputSection(nSamples_, nOutputs_, sampleOutputs_, 
                                        sampleStates_, outputNames_);
         psuadeIO_->updateMethodSection(PSUADE_SAMP_MC, nSamples_, 1, 0, 0);
         psuadeIO_->updateAnalysisSection(0, 0, 0, 0, 0, 0);
         nReps = 1;
         if (currSession != NULL) delete currSession;
         currSession = new PsuadeSession();
         psuadeIO_->getSession(currSession);
         printf("Excel data have been read.\n");
         printf("Use 'write' to store the data in PSUADE format.\n");
      }

      // +++ write_xls 
      else if (!strcmp(command, "write_xls"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("write_xls: write to a data file in Excel format\n");
            printf("           Use <read_xls -h> to see the Excel\n");
            printf("           format.\n");
            printf("syntax: write_xls <filename>.\n");
            printf("where <filename> is the name of the target data file.\n");
            printf("\nThe target file will have the following format:\n");
            printf("line 1: <nSamples> <nInputs> <nOutputs>\n");
            printf("line 2: # input and output names\n");
            printf("line 3: 1 <sample 1 inputs> < sample 1 outputs>\n");
            printf("line 4: 2 <sample 2 inputs> < sample 2 outputs>\n");
            printf(".......\n");
            continue;
         }
         if (psuadeIO_ == NULL || sampleOutputs_ == NULL)
         {
            printf("ERROR: no sample output data to write.\n");
            printf("       Use load to put data into local memory first.\n");
            continue;
         }
         strcpy(dataFile, "\0");
         sscanf(lineIn, "%s %s", command, dataFile);
         if (!strcmp(dataFile,psInputFilename_))
         {
            printf("WARNING : output file name should not be the same as\n"); 
            printf(" the input file name %s.\n", psInputFilename_);
            printf(" Command not executed.\n");
         }
         else if (!strcmp(dataFile, "\0"))
         {
            printf("ERROR: need to specify a file to write to.\n");
            printf("syntax: write_xls <filename>.\n");
         }
         else
         {
            if (sampleInputs_ == NULL || sampleOutputs_ == NULL)
               printf("ERROR: no inputs or outputs to write.\n");
            else
            {
               fp = fopen(dataFile, "w");
               if (fp == NULL)
               {
                  printf("ERROR: cannot open file %s.\n", dataFile);
                  continue;
               }
               fprintf(fp, "%d \t %d \t %d\n", nSamples_, nInputs_, nOutputs_);
               fprintf(fp, "#Sample \t ");
               for (ii = 0; ii < nInputs_; ii++)
                  fprintf(fp, "%s \t ", inputNames_[ii]);
               for (ii = 0; ii < nOutputs_-1; ii++)
                  fprintf(fp, "%s \t ", outputNames_[ii]);
               fprintf(fp, "%s \n", outputNames_[nOutputs_-1]);
               for (ii = 0; ii < nSamples_; ii++)
               {
                  fprintf(fp, "%d", ii+1);
                  for (jj = 0; jj < nInputs_; jj++)
                     fprintf(fp,"\t%24.16e", sampleInputs_[ii*nInputs_+jj]); 
                  for (jj = 0; jj < nOutputs_; jj++)
                     fprintf(fp,"\t%24.16e", sampleOutputs_[ii*nOutputs_+jj]);
                  fprintf(fp, "\n");
               }
               fclose(fp);
            }
         }
      }

      // +++ write_ultra 
      else if (!strcmp(command, "write_ultra"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("write_ultra: write to a data file in the ULTRA data\n");
            printf("             format (suitable for 2 input only).\n");
            printf("syntax: write_ultra <filename>.\n");
            printf("where <filename> is the name of the target data file.\n");
            continue;
         }
         if (psuadeIO_ == NULL || sampleOutputs_ == NULL)
         {
            printf("ERROR: no sample output data.\n");
            printf("       Use load to put data into local memory first.\n");
            continue;
         }
         strcpy(dataFile, "\0");
         sscanf(lineIn, "%s %s", command, dataFile);
         if (!strcmp(dataFile,psInputFilename_))
         {
            printf("WARNING : output file name should not be the same as\n"); 
            printf(" the input file name %s. \n", psInputFilename_);
            printf(" Command not executed.\n");
         }
         else if (!strcmp(dataFile, "\0"))
         {
            printf("ERROR: need to specify a file to write to.\n");
            printf("syntax: write_ultra <filename>.\n");
            printf("where <filename> will be a data file in ultra format.\n");
            continue;
         }
         else
         {
            if (sampleInputs_ == NULL || sampleOutputs_ == NULL)
               printf("ERROR: no inputs or outputs to save.\n");
            else if (nInputs_ != 2)
               printf("ERROR: command available only for nInputs=2.\n");
            else
            {
               if (nOutputs_ == 1) outputID = 0;
               else
               {
                  sprintf(pString, "Enter output number (1 - %d) : ",
                          nOutputs_);
                  outputID = getInt(1, nOutputs_, pString);
                  outputID--;
               }
               tempX = new double[nSamples_];
               tempW = new double[nSamples_];
               tempY = new double[nSamples_];
               tempT = new double[nSamples_];
               tempInds = new double[nSamples_];
               for (ii = 0; ii < nSamples_; ii++)
               {
                  tempX[ii] = sampleInputs_[ii*nInputs_]; 
                  tempW[ii] = sampleInputs_[ii*nInputs_+1]; 
                  tempY[ii] = sampleOutputs_[ii*nOutputs_+outputID];
                  tempInds[ii] = (double) ii;
               }
               sortDbleList2(nSamples_, tempW, tempInds);
               for (ii = 0; ii < nSamples_; ii++) tempT[ii] = tempX[ii];
               for (ii = 0; ii < nSamples_; ii++)
               {
                  ind2 = (int) tempInds[ii];
                  tempX[ii] = tempT[ind2];
               }
               for (ii = 0; ii < nSamples_; ii++) tempT[ii] = tempY[ii];
               for (ii = 0; ii < nSamples_; ii++)
               {
                  ind2 = (int) tempInds[ii];
                  tempY[ii] = tempT[ind2];
               }
               ind = 0;
               for (ii = 1; ii < nSamples_; ii++)
               {
                  if (tempW[ii] != tempW[ii-1])
                  {
                     sortDbleList2(ii-ind, &tempX[ind], &tempY[ind]);
                     ind = ii;
                  }
               } 
               sortDbleList2(nSamples_-ind, &tempX[ind], &tempY[ind]);
               fp = fopen(dataFile, "w");
               if (fp == NULL)
               {
                  printf("ERROR: cannot open output file %s\n", dataFile);
                  exit(1);
               }
               ind = 0;
               ind2 = 1;
               for (ii = 1; ii < nSamples_; ii++)
               {
                  if (tempW[ii] != tempW[ii-1])
                  {
                     fprintf(fp, "#var2_%d\n", ind2);
                     for (jj = ind; jj < ii; jj++)
                        fprintf(fp, "%24.16e %24.16e\n",tempX[jj],tempY[jj]);
                     ind = ii;
                     ind2++;
                  }
               } 
               fprintf(fp, "#var2_%d\n", ind2);
               for (jj = ind; jj < nSamples_; jj++)
                  fprintf(fp, "%14.6e %24.16e\n",tempX[jj],tempY[jj]);
               fclose(fp);
               delete [] tempInds;
               delete [] tempX;
               delete [] tempW;
               delete [] tempY;
               delete [] tempT;
               printf("ultra file created in %s\n", dataFile);
            }
         }
      }

      // +++ update 
      else if (!strcmp(command, "update"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("update: update the current sample outputs from another\n");
            printf("        PSUADE data file (i.e. the current unevaluated\n");
            printf("        sample points will be replaced with evaluated\n");
            printf("        outputs from the provided data file.) New points\n");
            printf("        in the data file will be ignored.\n");
            printf("syntax: update <filename> (merge from <filename>).\n");
            printf("where <filename> is a data file in PSUADE data format.\n");
            continue;
         }
         if (psuadeIO_ == NULL || sampleOutputs_ == NULL)
         {
            printf("ERROR: no sample output data to update.\n");
            printf("       Use load to put data into local memory first.\n");
            continue;
         }
         strcpy(dataFile, "\0");
         sscanf(lineIn, "%s %s", command, dataFile);
         status = 0;
         fp = fopen(dataFile, "r");
         if (fp == NULL)
         {
            printf("ERROR: File %s not found.\n", dataFile);
            printf("syntax: update <filename>).\n");
            printf("where <filename> is a data file in PSUADE data format.\n");
            status = 1;
         }
         else fclose(fp);
         ioPtr = NULL;
         if (status == 0)
         {
            ioPtr = new PsuadeData();
            status = ioPtr->readPsuadeFile(dataFile);
            if (status != 0)
               printf("ERROR: file %s either not found or in wrong format.\n",
                      dataFile);
         }
         if (status == 0)
         {
            ioPtr->getParameter("input_ninputs", pPtr);
            ind = pPtr.intData_;
            flag = 1;
            if (ind != nInputs_) flag = 0;
            if (flag == 1)
            {
               pINames.clean();
               pLower.clean();
               pUpper.clean();
               ioPtr->getParameter("input_names", pINames);
               names = pINames.strArray_;
               ioPtr->getParameter("input_lbounds", pLower);
               tempX = pLower.dbleArray_;
               ioPtr->getParameter("input_ubounds", pUpper);
               tempY = pUpper.dbleArray_;
               tempX = NULL;
               tempY = NULL;
               names = NULL;
               pINames.clean();
               pLower.clean();
               pUpper.clean();
            }
            if (flag == 1)
            {
               ioPtr->getParameter("output_noutputs", pPtr);
               ind = pPtr.intData_;
               if (ind != nOutputs_) flag = 0;
            }
            if (flag == 1)
            {
               pONames.clean();
               ioPtr->getParameter("output_names", pONames);
               names = pONames.strArray_;
               for (ii = 0; ii < ind; ii++)
                  if (strcmp(outputNames_[ii],names[ii])) flag = 0;
               pONames.clean();
            }
            if (flag == 0) printf("ERROR: invalid data file.\n");
            else
            {
               pInputs.clean();
               pOutputs.clean();
               pStates.clean();
               ioPtr->getParameter("input_sample",pInputs);
               tempX = pInputs.dbleArray_;
               ioPtr->getParameter("output_sample",pOutputs);
               tempY = pOutputs.dbleArray_;
               ioPtr->getParameter("output_states",pStates);
               states = pStates.intArray_;
               ioPtr->getParameter("method_nsamples", pPtr);
               ind = pPtr.intData_;
               count = 0;
               for (sInd = 0; sInd < nSamples_; sInd++)
               {
                  if (sampleStates_[sInd] == 1) continue;
                  kk = sInd * nInputs_;
                  for (ii = 0; ii < ind; ii++)
                  {
                     if (states[ii] != 1) continue;
                     jj = ii * nInputs_;
                     for (iInd = 0; iInd < nInputs_; iInd++)
                        if (sampleInputs_[kk+iInd] != tempX[jj+iInd])
                           break;
                     if (iInd == nInputs_)
                     {
                        kk = sInd * nOutputs_;
                        jj = ii * nOutputs_;
                        for (oInd = 0; oInd < nOutputs_; oInd++)
                           sampleOutputs_[kk+oInd] = tempY[jj+oInd];
                        sampleStates_[sInd] = 1;
                        printf("   Matched: sample %d <-- sample %d (%d)\n",
                               sInd+1,ii+1,ind);
                        count++;
                        break;
                     }
                  }
               }
               printf("Number of sample points updated = %d\n", count);
               tempX = NULL;
               tempY = NULL;
               states = NULL;
               pInputs.clean();
               pOutputs.clean();
               pStates.clean();
            }
         }
         psuadeIO_->updateOutputSection(nSamples_,nOutputs_,sampleOutputs_,
                                        sampleStates_, NULL); 
         if (currSession != NULL) delete currSession;
         currSession = new PsuadeSession();
         psuadeIO_->getSession(currSession);
         fflush(stdout);
         if (ioPtr != NULL) delete ioPtr;
         ioPtr = NULL;
      }

      // +++ moat_adjust 
      else if (!strcmp(command, "moat_adjust"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("moat_adjust: adjust a MOAT sample using new input values\n");
            printf("syntax: moat_adjust <filename>.\n");
            printf("where <filename> is a MOAT adjust file (created by moatgen)\n");
            continue;
         }
         if (psuadeIO_ == NULL || sampleOutputs_ == NULL)
         {
            printf("ERROR: no Morris sample to adjust.\n");
            printf("       Use load to put data into local memory first.\n");
            continue;
         }
         strcpy(dataFile, "\0");
         sscanf(lineIn, "%s %s", command, dataFile);
         status = 0;
         fp = fopen(dataFile, "r");
         if (fp == NULL)
         {
            printf("ERROR: File %s not found.\n", dataFile);
            printf("syntax: moat_adjust <filename>.\n");
            printf("where <filename> is a MOAT adjust file (created by moatgen)\n");
            status = 1;
         }
         else fclose(fp);

         if (samplingMethod != PSUADE_SAMP_MOAT)
         {
            printf("Cannot adjust: the sampling method is not MOAT.\n");
            status = 1;
         }
         if (status == 0)
         {
            sampPtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_MOAT);
            sampPtr->setPrintLevel(PL_INTERACTIVE);
            sampPtr->setInputBounds(nInputs_, iLowerB_, iUpperB_);
            sampPtr->setOutputParams(nOutputs_);
            sampPtr->setSamplingParams(nSamples_,nSamples_/(nInputs_+1),0);
            sampPtr->initialize(1);
            sampPtr->loadSamples(nSamples_, nInputs_, nOutputs_, sampleInputs_,
                                 sampleOutputs_, sampleStates_);
            sampPtr->repair(dataFile, 0);
            sampPtr->getSamples(nSamples_, nInputs_, nOutputs_, sampleInputs_,
                                sampleOutputs_, sampleStates_);
            delete sampPtr;
            sampPtr = NULL;
         }
         fflush(stdout);
      }

      // +++ gmoat_adjust 
      else if (!strcmp(command, "gmoat_adjust"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("gmoat_adjust: adjust a MOAT sample using new input values\n");
            printf("Note: gmoat_adjust differs from moat_adjust in that it\n");
            printf("      uses GMOAT instead of MOAT to create the initial\n");
            printf("      sample.\n");
            printf("syntax: gmoat_adjust <filename>.\n");
            printf("where <filename> is a MOAT adjust file (created by moatgen)\n");
            continue;
         }
         if (psuadeIO_ == NULL || sampleOutputs_ == NULL)
         {
            printf("ERROR: no Morris sample to adjust.\n");
            printf("       Use load to put data into local memory first.\n");
            continue;
         }
         strcpy(dataFile, "\0");
         sscanf(lineIn, "%s %s", command, dataFile);
         status = 0;
         fp = fopen(dataFile, "r");
         if (fp == NULL)
         {
            printf("ERROR: File %s not found.\n", dataFile);
            printf("syntax: gmoat_adjust <file>.\n");
            printf("where <file> is a MOAT adjust file (created by moatgen)\n");
            status = 1;
         }
         else fclose(fp);

         if (samplingMethod != PSUADE_SAMP_GMOAT)
         {
            printf("Cannot adjust: the sampling method is not GMOAT.\n");
            status = 1;
         }
         if (status == 0)
         {
            sampPtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_GMOAT);
            sampPtr->setPrintLevel(PL_INTERACTIVE);
            sampPtr->setInputBounds(nInputs_, iLowerB_, iUpperB_);
            sampPtr->setOutputParams(nOutputs_);
            sampPtr->setSamplingParams(nSamples_,nSamples_/(nInputs_+1),0);
            sampPtr->initialize(1);
            sampPtr->loadSamples(nSamples_, nInputs_, nOutputs_, sampleInputs_,
                                 sampleOutputs_, sampleStates_);
            sampPtr->repair(dataFile, 0);
            sampPtr->getSamples(nSamples_, nInputs_, nOutputs_, sampleInputs_,
                                sampleOutputs_, sampleStates_);
            delete sampPtr;
            sampPtr = NULL;
         }
         fflush(stdout);
      }

      // +++ iadd 
      else if (!strcmp(command, "iadd"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("iadd: add more inputs to the existing sample from\n");
            printf("      another PSUADE data file\n");
            printf("       syntax: iadd <filename>\n");
            printf("       where <filename> should be a PSUADE data file\n");
            printf("                        containing additional inputs.\n");
            printf("Note: the sample size in the data file should be the\n");
            printf("      same as that in the resident memory.\n");
            continue;
         }
         if (psuadeIO_ == NULL || sampleOutputs_ == NULL)
         {
            printf("ERROR: no sample output data.\n");
            printf("       Use load to put data into local memory first.\n");
            continue;
         }
         strcpy(dataFile, "\0");
         sscanf(lineIn, "%s %s", command, dataFile);
         status = 0;
         fp = fopen(dataFile, "r");
         if (fp == NULL)
         {
            printf("ERROR: File %s not found.\n", dataFile);
            printf("       syntax: iadd <file>\n");
            printf("       where <file> should be a PSUADE data file.\n");
            status = 1;
         }
         else fclose(fp);
         ioPtr = NULL;
         if (status == 0)
         {
            ioPtr = new PsuadeData();
            status = ioPtr->readPsuadeFile(dataFile);
            if (status != 0)
               printf("ERROR: file %s either not found or in wrong format.\n",
                      dataFile);
         }
         if (status == 0)
         {
            ioPtr->getParameter("method_nsamples", pPtr);
            ind = pPtr.intData_;
            if (ind != nSamples_)
            {
               printf("ERROR: nSamples not the same in both sets of data.\n");
               printf("       nSamples in local memory = %d\n", nSamples_);
               printf("       nSamples from file       = %d\n", ind);
            }
            else
            {
               ioPtr->getParameter("input_ninputs", pPtr);
               ind = pPtr.intData_;

               tempW = iLowerB_;
               ioPtr->getParameter("input_lbounds", pLower);
               tempT = pLower.dbleArray_;
               iLowerB_ = new double[nInputs_+ind];
               for (ii = 0; ii < nInputs_; ii++) iLowerB_[ii] = tempW[ii];
               for (ii = 0; ii < ind; ii++) iLowerB_[nInputs_+ii] = tempT[ii];
               pLower.clean();
               delete [] tempW;
               
               tempW = iUpperB_;
               ioPtr->getParameter("input_ubounds", pUpper);
               tempT = pUpper.dbleArray_;
               iUpperB_ = new double[nInputs_+ind];
               for (ii = 0; ii < nInputs_; ii++) iUpperB_[ii] = tempW[ii];
               for (ii = 0; ii < ind; ii++) iUpperB_[nInputs_+ii] = tempT[ii];
               pUpper.clean();
               delete [] tempW;

               tempI = inputPDFs_;
               ioPtr->getParameter("input_pdfs", pPDFs);
               inputPDFs_ = new int[nInputs_+ind];
               for (jj = 0; jj < nInputs_; jj++) inputPDFs_[jj] = tempI[jj];
               delete [] tempI;
               tempI = pPDFs.intArray_;
               for (jj = nInputs_; jj < nInputs_+ind; jj++)
                  inputPDFs_[jj] = tempI[jj-nInputs_];
               pPDFs.clean();

               ioPtr->getParameter("input_means", pMeans);
               tempW = inputMeans_;
               inputMeans_ = new double[nInputs_+ind];
               for (jj = 0; jj < nInputs_; jj++) inputMeans_[jj] = tempW[jj];
               delete [] tempW;
               tempW = pMeans.dbleArray_;
               for (jj = nInputs_; jj < nInputs_+ind; jj++)
                  inputMeans_[jj] = tempW[jj-nInputs_];
               pMeans.clean();

               ioPtr->getParameter("input_stdevs", pStds);
               tempW = inputStds_;
               inputStds_ = new double[nInputs_+ind];
               for (jj = 0; jj < nInputs_; jj++) inputStds_[jj] = tempW[jj];
               delete [] tempW;
               tempW = pStds.dbleArray_;
               for (jj = nInputs_; jj < nInputs_+ind; jj++)
                  inputStds_[jj] = tempW[jj-nInputs_];
               pStds.clean();

               ioPtr->getParameter("input_cor_matrix", pPtr);
               Matrix *tmpMat1 = (Matrix *) pPtr.psObject_;
               Matrix *tmpMat2 = new Matrix();
               tmpMat2->setDim(nInputs_+ind,nInputs_+ind);
               for (ii = 0; ii < nInputs_; ii++)
               {
                  for (jj = 0; jj < nInputs_; jj++)
                  {
                     ddata = inputCMat_->getEntry(ii,jj);
                     tmpMat2->setEntry(ii,jj,ddata);
                  }
               }
               for (ii = nInputs_; ii < nInputs_+ind; ii++)
               {
                  for (jj = nInputs_; jj < nInputs_+ind; jj++)
                  {
                     ddata = tmpMat1->getEntry(ii-nInputs_,jj-nInputs_);
                     tmpMat2->setEntry(ii,jj,ddata);
                  }
               }
               delete inputCMat_;
               inputCMat_ = tmpMat2;

               names = inputNames_;
               inputNames_ = new char*[nInputs_+ind];
               for (ii = 0; ii < nInputs_; ii++) inputNames_[ii] = names[ii];
               if (names != NULL) delete [] names;
               pINames.clean();
               ioPtr->getParameter("input_names", pINames);
               names = pINames.strArray_;
               for (ii = 0; ii < ind; ii++)
               {
                  inputNames_[nInputs_+ii] = names[ii];
                  pINames.strArray_[ii] = NULL;
               }
               pINames.clean();
               ioPtr->getParameter("input_sample", pInputs);
               tempW = pInputs.dbleArray_;
               tempX = sampleInputs_;
               sampleInputs_ = new double[nSamples_*(nInputs_+ind)];
               kk = 0;
               for (ii = 0; ii < nSamples_; ii++)
               {
                  for (jj = 0; jj < nInputs_; jj++)
                     sampleInputs_[kk++] = tempX[ii*nInputs_+jj];
                  for (jj = 0; jj < ind; jj++)
                     sampleInputs_[kk++] = tempW[ii*ind+jj];
               }
               delete [] tempX;
               nInputs_ += ind;
               pInputs.clean();
               psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,
                              iLowerB_,iUpperB_,sampleInputs_,inputNames_, 
                              inputPDFs_, inputMeans_,inputStds_,inputCMat_); 
            }
            if (currSession != NULL) delete currSession;
            currSession = new PsuadeSession();
            psuadeIO_->getSession(currSession);
         }
         fflush(stdout);
         if (ioPtr != NULL) delete ioPtr;
         ioPtr = NULL;
      }

      // +++ oadd 
      else if (!strcmp(command, "oadd"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("oadd: add more outputs to the existing sample from\n");
            printf("      another PSUADE data file\n");
            printf("      syntax: oadd <filename>\n");
            printf("      where <filename> should be a PSUADE data file\n");
            printf("                       containing additional outputs.\n\n");
            printf("Note: the sample size in the data file should be the\n");
            printf("      same as that in the resident memory.\n");
            continue;
         }
         if (psuadeIO_ == NULL || sampleOutputs_ == NULL)
         {
            printf("ERROR: no data to add to.\n");
            printf("       Use load to put data into local memory first.\n");
            continue;
         }
         strcpy(dataFile, "\0");
         sscanf(lineIn, "%s %s", command, dataFile);
         status = 0;
         fp = fopen(dataFile, "r");
         if (fp == NULL)
         {
            printf("ERROR: File %s not found.\n", dataFile);
            printf("       syntax: oadd <file>\n");
            printf("       where <file> should be a PSUADE data file.\n");
            status = 1;
         }
         else fclose(fp);
         ioPtr = NULL;
         if (status == 0)
         {
            ioPtr = new PsuadeData();
            status = ioPtr->readPsuadeFile(dataFile);
            if (status != 0)
               printf("ERROR: file %s either not found or in wrong format.\n",
                      dataFile);
         }
         if (status == 0)
         {
            ioPtr->getParameter("method_nsamples", pPtr);
            ind = pPtr.intData_;
            if (ind != nSamples_)
            {
               printf("ERROR: nSamples not the same in both sets of data.\n");
               printf("       nSamples in local memory = %d\n", nSamples_);
               printf("       nSamples from file       = %d\n", ind);
            }
            else
            {
               printf("INFO: no check to see if both sets of inputs match.\n");
               ioPtr->getParameter("output_noutputs", pPtr);
               ind = pPtr.intData_;

               names = outputNames_;
               outputNames_ = new char*[nOutputs_+ind];
               for (ii = 0; ii < nOutputs_; ii++) outputNames_[ii] = names[ii];
               if (names != NULL) delete [] names;
               pONames.clean();
               ioPtr->getParameter("output_names", pONames);
               names = pONames.strArray_;
               for (ii = 0; ii < ind; ii++)
               {
                  outputNames_[nOutputs_+ii] = names[ii];
                  names[ii] = NULL;
               }
               pONames.clean();
               ioPtr->getParameter("output_sample", pOutputs);
               tempW = pOutputs.dbleArray_;
               tempX = sampleOutputs_;
               sampleOutputs_ = new double[nSamples_*(nOutputs_+ind)];
               kk = 0;
               for (ii = 0; ii < nSamples_; ii++)
               {
                  for (jj = 0; jj < nOutputs_; jj++)
                     sampleOutputs_[kk++] = tempX[ii*nOutputs_+jj];
                  for (jj = 0; jj < ind; jj++)
                     sampleOutputs_[kk++] = tempW[ii*ind+jj];
               }
               delete [] tempX;
               nOutputs_ += ind;
               pOutputs.clean();
               psuadeIO_->updateOutputSection(nSamples_,nOutputs_,
                                sampleOutputs_, sampleStates_, outputNames_); 
               if (currSession != NULL) delete currSession;
               currSession = new PsuadeSession();
               psuadeIO_->getSession(currSession);
            }
         }
         fflush(stdout);
         if (ioPtr != NULL) delete ioPtr;
         ioPtr = NULL;
      }

      // +++ iadd1 
      else if (!strcmp(command, "iadd1"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("iadd: add one input to the existing sample (at the end).\n");
            printf("      The values of the new input are drawn randomly\n");
            printf("      between 0 and 1.\n");
            printf("syntax: iadd1 \n");
            continue;
         }
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no data to add to.\n");
            printf("       Use load to put data into local memory first.\n");
            continue;
         }
         tempW = iLowerB_;
         iLowerB_ = new double[nInputs_+1];
         for (ii = 0; ii < nInputs_; ii++) iLowerB_[ii] = tempW[ii];
         iLowerB_[nInputs_] = 0;
         delete [] tempW;
         tempW = iUpperB_;
         iUpperB_ = new double[nInputs_+1];
         for (ii = 0; ii < nInputs_; ii++) iUpperB_[ii] = tempW[ii];
         iUpperB_[nInputs_] = 1;
         delete [] tempW;

         if (inputPDFs_ != NULL)
         {
            tempI = inputPDFs_;
            inputPDFs_ = new int[nInputs_+1];
            for (jj = 0; jj < nInputs_; jj++) inputPDFs_[jj] = tempI[jj];
            inputPDFs_[nInputs_] = 0;
            delete [] tempI;
         }
         if (inputMeans_ != NULL)
         {
            tempW = inputMeans_;
            inputMeans_ = new double[nInputs_+1];
            for (jj = 0; jj < nInputs_; jj++) inputMeans_[jj] = tempW[jj];
            inputMeans_[nInputs_] = 0.0;
            delete [] tempW;
         }
         if (inputStds_ != NULL)
         {
            tempW = inputStds_;
            inputStds_ = new double[nInputs_+1];
            for (jj = 0; jj < nInputs_; jj++) inputStds_[jj] = tempW[jj];
            inputStds_[nInputs_] = 0.0;
            delete [] tempW;
         }
         if (inputCMat_ != NULL)
         {
            Matrix *tmpMat2 = new Matrix();
            tmpMat2->setDim(nInputs_+1,nInputs_+1);
            for (ii = 0; ii < nInputs_; ii++)
            {
               for (jj = 0; jj < nInputs_; jj++)
               {
                  ddata = inputCMat_->getEntry(ii,jj);
                  tmpMat2->setEntry(ii,jj,ddata);
               }
            }
            ddata = 1.0;
            tmpMat2->setEntry(nInputs_,nInputs_,ddata);
            delete inputCMat_;
            inputCMat_ = tmpMat2;
         }
         if (inputNames_ != NULL)
         {
            names = inputNames_;
            inputNames_ = new char*[nInputs_+1];
            for (ii = 0; ii < nInputs_; ii++) inputNames_[ii] = names[ii];
            delete [] names;
            inputNames_[nInputs_] = new char[100];
            sprintf(inputNames_[nInputs_], "X%d", nInputs_+1);
         }
         tempX = sampleInputs_;
         sampleInputs_ = new double[nSamples_*(nInputs_+1)];
         for (ii = 0; ii < nSamples_; ii++)
         {
            for (jj = 0; jj < nInputs_; jj++)
               sampleInputs_[ii*(nInputs_+1)+jj] = tempX[ii*nInputs_+jj];
            sampleInputs_[ii*(nInputs_+1)+nInputs_] = PSUADE_drand();
         }
         nInputs_++;
         psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,iLowerB_,
                          iUpperB_,sampleInputs_,inputNames_,inputPDFs_, 
                          inputMeans_,inputStds_,inputCMat_); 
         if (currSession != NULL) delete currSession;
         currSession = new PsuadeSession();
         psuadeIO_->getSession(currSession);
         fflush(stdout);
         printf("iadd1 : one input added (use sinfo to see updated sample).\n");
      }

      // +++ oadd1 
      else if (!strcmp(command, "oadd1"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("iadd: add one output to the existing sample.\n");
            printf("      The values of the new output are drawn randomly\n");
            printf("      between 0 and 1.\n");
            printf("syntax: oadd1 \n");
            continue;
         }
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no data to add to.\n");
            printf("       Use load to put data into local memory first.\n");
            continue;
         }
         if (sampleOutputs_ != NULL)
         {
            if (outputNames_ != NULL)
            {
               names = outputNames_;
               outputNames_ = new char*[nOutputs_+1];
               for (ii = 0; ii < nOutputs_; ii++) outputNames_[ii] = names[ii];
               delete [] names;
               outputNames_[nOutputs_] = new char[100];
               sprintf(outputNames_[nOutputs_], "Y%d", nOutputs_+1);
            }
            tempX = sampleOutputs_;
            sampleOutputs_ = new double[nSamples_*(nOutputs_+1)];
            for (ii = 0; ii < nSamples_; ii++)
            {
               for (jj = 0; jj < nOutputs_; jj++)
                  sampleOutputs_[ii*(nOutputs_+1)+jj] = tempX[ii*nOutputs_+jj];
               sampleOutputs_[ii*(nOutputs_+1)+nOutputs_] = PSUADE_drand();
            }
            nOutputs_++;
            psuadeIO_->updateOutputSection(nSamples_,nOutputs_,sampleOutputs_, 
                             sampleStates_, outputNames_); 
         }
         if (currSession != NULL) delete currSession;
         currSession = new PsuadeSession();
         psuadeIO_->getSession(currSession);
         printf("oadd1 : one output added (use sinfo to see updated sample).\n");
         fflush(stdout);
      }

      // +++ ireplace 
      else if (!strcmp(command, "ireplace"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("ireplace: replace an input in the existing sample from\n");
            printf("          another sample file in ASCII format.\n");
            printf("syntax: ireplace <filename>\n");
            printf("where <filename> should be in the format given below: \n");
            printf("line 1: PSUADE_BEGIN\n");
            printf("line 2: <nSamples>\n");
            printf("line 3: input value for sample point 1\n");
            printf("line 4: input value for sample point 2\n");
            printf(".......\n");
            printf("last line: PSUADE_END\n");
            continue;
         }
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no data to ireplace.\n");
            printf("       Use load to put data into local memory first.\n");
            continue;
         }
         strcpy(dataFile, "\0");
         sscanf(lineIn, "%s %s", command, dataFile);
         status = 0;
         fp = fopen(dataFile, "r");
         if (fp == NULL)
         {
            printf("ERROR: File %s not found.\n", dataFile);
            printf("Use <ireplace -h> to see command syntax\n");
            status = 1;
         }
         else fclose(fp);
         if (status == 0)
         {
            fp = fopen(dataFile, "r");
            fscanf(fp, "%s", winput);
            if (strcmp(winput, "PSUADE_BEGIN"))
            {
               printf("ERROR: file must begin with PSUADE_BEGIN\n");
               fclose(fp);
               continue;
            }
            fscanf(fp, "%d", &kk);
            if (kk != nSamples_)
            {
               fclose(fp);
               printf("ERROR: File and local data do not match.\n");
               printf("     size of sample in local memory = %d\n",nSamples_);
               printf("     nSamples from external file    = %d\n",kk);
               status = 1;
            }
         }
         if (status == 0)
         {
            inputID = 0;
            if (nInputs_ > 1)
            {
               sprintf(pString,"Enter input number to replace (1 - %d) : ",
                       nInputs_);
               inputID = getInt(1, nInputs_, pString);
               inputID--;
            }
            for (ii = 0; ii < nSamples_; ii++)
               fscanf(fp, "%lg", &(sampleInputs_[ii*nInputs_+inputID]));
            Xmin = PSUADE_UNDEFINED;
            Xmax = - Xmin;
            for (ii = 0; ii < nSamples_; ii++)
            {
               if (sampleInputs_[ii*nInputs_+inputID] < Xmin)
                  Xmin = sampleInputs_[ii*nInputs_+inputID];
               else if (sampleInputs_[ii*nInputs_+inputID] > Xmax)
                  Xmax = sampleInputs_[ii*nInputs_+inputID];
            }
            iLowerB_[inputID] = Xmin; 
            iUpperB_[inputID] = Xmax; 
            inputPDFs_[inputID] = 0;
            inputMeans_[inputID] = 0;
            inputStds_[inputID] = 0;
            for (ii = 0; ii < nInputs_; ii++)
               inputCMat_->setEntry(inputID,ii,0);
            for (ii = 0; ii < nInputs_; ii++)
               inputCMat_->setEntry(ii,inputID,0);
            inputCMat_->setEntry(inputID,inputID,1);
            sprintf(inputNames_[inputID], "X%d", inputID+1);
            psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,iLowerB_,
                          iUpperB_,sampleInputs_,inputNames_, 
                          inputPDFs_, inputMeans_,inputStds_,inputCMat_); 
            fscanf(fp, "%s", winput);
            if (strcmp(winput, "PSUADE_END"))
               printf("WARNING: file should end with PSUADE_END\n");
            fclose(fp);
            if (currSession != NULL) delete currSession;
            currSession = new PsuadeSession();
            psuadeIO_->getSession(currSession);
         }
         fflush(stdout);
      }

      // +++ oreplace 
      else if (!strcmp(command, "oreplace"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("oreplace: replace all outputs in the existing sample\n");
            printf("          from another sample file in ASCII format.\n");
            printf("syntax: oreplace <file>\n");
            printf("where <file> should be in the format given below: \n");
            printf("line 1: PSUADE_BEGIN\n");
            printf("line 2: <nSamples> <nOutputs>\n");
            printf("line 3: output value for sample point 1\n");
            printf("line 4: output value for sample point 2\n");
            printf(".......\n");
            printf("last line: PSUADE_END\n");
            continue;
         }
         if (psuadeIO_ == NULL || sampleOutputs_ == NULL)
         {
            printf("ERROR: no data to oreplace.\n");
            printf("       Use load to put data into local memory first.\n");
            continue;
         }
         strcpy(dataFile, "\0");
         sscanf(lineIn, "%s %s", command, dataFile);
         status = 0;
         fp = fopen(dataFile, "r");
         if (fp == NULL)
         {
            printf("ERROR: File %s not found.\n", dataFile);
            printf("Use <oreplace -h> to see command syntax\n");
            status = 1;
         }
         else fclose(fp);
         if (status == 0)
         {
            fp = fopen(dataFile, "r");
            fscanf(fp, "%s", winput);
            if (strcmp(winput, "PSUADE_BEGIN"))
            {
               printf("ERROR: file must begin with PSUADE_BEGIN\n");
               fclose(fp);
               continue;
            }
            fscanf(fp, "%d %d", &kk, &nOut2);
            if (kk != nSamples_ || nOut2 <= 0)
            {
               printf("ERROR: File and local parameters do not match.\n");
               printf("       nSamples (local) = %d\n", nSamples_);
               printf("       nSamples (file)  = %d\n", kk);
               printf("       nOutputs (file)  = %d\n", nOut2);
               fclose(fp);
               status = 1;
            }
         }
         if (status == 0)
         {
            printf("This command only checks that nSamples match.\n");
            printf("So make sure the data are correctly ordered.\n");
            printf("       nSamples (local) = %d\n", nSamples_);
            printf("       nOutputs (file)  = %d\n", nOut2);
            if (nOutputs_ != nOut2)
            {
               if (outputNames_ != NULL) 
               {
                  for (ii = 0; ii < nOutputs_; ii++)
                     delete [] outputNames_[ii];
                  delete [] outputNames_;
               }
               nOutputs_ = nOut2;
               outputNames_ = new char*[nOutputs_];
               for (ii = 0; ii < nOutputs_; ii++)
               {
                  outputNames_[ii] = new char[200];
                  sprintf(outputNames_[ii], "Y%d", ii+1);
               }
               if (sampleOutputs_ != NULL) delete [] sampleOutputs_;
               sampleOutputs_ = NULL;
            }
            if (sampleOutputs_ == NULL)
               sampleOutputs_ = new double[nSamples_*nOutputs_]; 
            for (ii = 0; ii < nSamples_*nOutputs_; ii++)
               fscanf(fp, "%lg", &(sampleOutputs_[ii]));
            for (ii = 0; ii < nSamples_; ii++) sampleStates_[ii] = 1;
            psuadeIO_->updateOutputSection(nSamples_,nOutputs_,
                          sampleOutputs_,sampleStates_,outputNames_); 
            fscanf(fp, "%s", winput);
            if (strcmp(winput, "PSUADE_END"))
               printf("WARNING: file should end with PSUADE_END\n");
            if (currSession != NULL) delete currSession;
            currSession = new PsuadeSession();
            psuadeIO_->getSession(currSession);
         }
         fflush(stdout);
      }

      // +++ moatgen 
      else if (!strcmp(command, "moatgen"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("moatgen: create a Morris adjust file with one constraint\n");
            printf("syntax: moatgen\n");
            printf("Note: a PSUADE datafile should have been loaded before\n");
            printf("      using this command. The data will be used as\n");
            printf("      response surface to find feasible region in\n");
            printf("      creating a MOAT adjust sample.\n");
            printf("Info: when the parameter space is not a hyper-rectangle,\n");
            printf("      the standard MOAT method will have difficulties.\n");
            printf("      This command circumvents the difficulties by creating\n");
            printf("      MOAT sample points that live inside the parameter boundary.\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data to moatgen.\n");
            printf("       Use load to put data into local memory first.\n");
            continue;
         }
         faFlag = 3;
         faPtr = genFAInteractive(psuadeIO_, faFlag);
         if (faPtr == NULL) {printf("ERROR detected in RS.\n"); continue;}
         faPtr->setOutputLevel(outputLevel_);
         psuadeIO_->getParameter("ana_outputid", pPtr);
         outputID = pPtr.intData_;
         Ymax = - 1.0e35;
         Ymin =   1.0e35;
         for (sInd = 0; sInd < nSamples_; sInd++)
         {
            if (sampleOutputs_[sInd*nOutputs_+outputID] > Ymax)
               Ymax = sampleOutputs_[sInd*nOutputs_+outputID];
            if (sampleOutputs_[sInd*nOutputs_+outputID] < Ymin)
               Ymin = sampleOutputs_[sInd*nOutputs_+outputID];
         }
         sprintf(pString,
                 "Enter the lower bound constraint (Ymin=%e) : ",Ymin);
         threshL = getDouble(pString);
         sprintf(pString,
                 "Enter the upper bound constraint (Ymax=%e) : ",Ymax);
         threshU = getDouble(pString);
         if (threshL >= threshU)
         {
            printf("ERROR: lower bound >= upper bound.\n");
            continue;
         }
         nPaths = 5000;
         sprintf(pString, "Enter P (resolution: try 4-10) : ");
         currP = getInt(4, 10, pString);
         sprintf(pString, "Enter the number of trials (> 100) : ");
         nTrials = getInt(101, 10000000, pString);
         moatSample = new double*[nPaths*(nInputs_+1)];
         for (ii = 0; ii < nPaths*(nInputs_+1); ii++)
            moatSample[ii] = new double[nInputs_];
         tempW = new double[nInputs_];
         indSet = new int[nInputs_];
         count = 0;
         filterRange = threshU - threshL;
         for (ii = 0; ii < nPaths; ii++)
         {
            trial = 0; 
            while (trial < nTrials)
            {
               iInd = count;
               trial++;
               for (jj = 0; jj < nInputs_; jj++)
               {
                  ind = PSUADE_rand() % currP;
                  dtemp = ind * (iUpperB_[jj] - iLowerB_[jj]) / (currP - 1.0);
                  moatSample[iInd][jj] = iLowerB_[jj] + dtemp;
               }
               dtemp = faPtr->evaluatePoint(moatSample[iInd]);
               if (dtemp < threshL || dtemp > threshU) continue;

               generateRandomIvector(nInputs_, indSet);

               for (jj = 0; jj < nInputs_; jj++)
               {
                  iInd++;
                  for (kk = 0; kk < nInputs_; kk++)
                     moatSample[iInd][kk] = moatSample[iInd-1][kk];

                  ind2 = indSet[jj];

                  ddata = moatSample[iInd][ind2]; 
                  moatSample[iInd][ind2] = iLowerB_[ind2];
                  currY1 = faPtr->evaluatePoint(moatSample[iInd]);
                  moatSample[iInd][ind2] = iUpperB_[ind2];
                  currY2 = faPtr->evaluatePoint(moatSample[iInd]);
                  currX1 = iLowerB_[ind2];
                  currX2 = iUpperB_[ind2];
                  moatSample[iInd][ind2] = ddata;

                  if (currY2 >= threshU && currY1 >= threshU)
                     currX1 = currX2 = 0.0;
                  else if (currY2 <= threshL && currY1 <= threshL)
                     currX1 = currX2 = 0.0;
                  else if (currY2 > currY1)
                  {
                     if (currY2 <= threshU) currX2 = iUpperB_[ind2];
                     else
                     {
                        sLo = iLowerB_[ind2];
                        sHi = iUpperB_[ind2];
                        while (PABS((currY2-threshU)/filterRange)>0.0001)
                        {
                           moatSample[iInd][ind2] = 0.5 * (sLo + sHi);
                           currY2 = faPtr->evaluatePoint(moatSample[iInd]);
                           if (currY2 > threshU) sHi = 0.5 * (sLo + sHi);
                           else                  sLo = 0.5 * (sLo + sHi);
                        }
                        currX2 = moatSample[iInd][ind2];
                     }
                     if (currY1 >= threshL) currX1 = iLowerB_[ind2];
                     else
                     {
                        sLo = iLowerB_[ind2];
                        sHi = iUpperB_[ind2];
                        while (PABS((currY1-threshL)/filterRange)>0.0001)
                        {
                           moatSample[iInd][ind2] = 0.5 * (sLo + sHi);
                           currY1 = faPtr->evaluatePoint(moatSample[iInd]);
                           if (currY1 < threshL) sLo = 0.5 * (sLo + sHi);
                           else                  sHi = 0.5 * (sLo + sHi);
                        }
                        currX1 = moatSample[iInd][ind2];
                     }
                  }
                  else
                  {
                     if (currY1 <= threshU) currX1 = iLowerB_[ind2];
                     else
                     {
                        sLo = iLowerB_[ind2];
                        sHi = iUpperB_[ind2];
                        while (PABS((currY1-threshU)/filterRange)>0.0001)
                        {
                           moatSample[iInd][ind2] = 0.5 * (sLo + sHi);
                           currY1 = faPtr->evaluatePoint(moatSample[iInd]);
                           if (currY1 > threshU) sLo = 0.5 * (sLo + sHi);
                           else                  sHi = 0.5 * (sLo + sHi);
                        }
                        currX1 = moatSample[iInd][ind2];
                     }
                     if (currY2 >= threshL) currX2 = iUpperB_[ind2];
                     else
                     {
                        sLo = iLowerB_[ind2];
                        sHi = iUpperB_[ind2];
                        while (PABS((currY2-threshL)/filterRange)>0.0001)
                        {
                           moatSample[iInd][ind2] = 0.5 * (sLo + sHi);
                           currY2 = faPtr->evaluatePoint(moatSample[iInd]);
                           if (currY2 < threshL) sHi = 0.5 * (sLo + sHi);
                           else                  sLo = 0.5 * (sLo + sHi);
                        }
                        currX2 = moatSample[iInd][ind2];
                     }
                  }
                  if (PABS(currX2-currX1)<0.1*(iUpperB_[ind2]-iLowerB_[ind2])) 
                     break;
                  moatSample[iInd][ind2] = ddata;
                  tempW[ind2] = PABS(currX2 - currX1) / (currP-1.0);
                  moatSample[iInd][ind2] += tempW[ind2];
                  if (moatSample[iInd][ind2] > iUpperB_[ind2])
                       dtemp = threshL - 1.0;
                  else dtemp = faPtr->evaluatePoint(moatSample[iInd]);
                  if (dtemp < threshL || dtemp > threshU) 
                  {
                     moatSample[iInd][ind2] -= 2.0 * tempW[ind2];
                     if (moatSample[iInd][ind2] < iLowerB_[ind2])
                        break;
                     dtemp = faPtr->evaluatePoint(moatSample[iInd]);
                     if (dtemp < threshL || dtemp > threshU) 
                        break;
                  }
               }
               if (jj == nInputs_) 
               {
                  count += (nInputs_ + 1);
                  if (outputLevel_ > 2)
                     printf("moatgen: path %d (out of %d) found.\n", ii+1,
                            nPaths);
                  break; 
               }
               else
               {
                  if (outputLevel_ > 2)
                     printf("Current path fails (%d out of max %d).\n",
                            trial, nTrials); 
               }
            }
            if (trial >= nTrials)
            {
               printf("moatgen fails to find all possible paths.\n");
               printf("Suggestion: try again with a larger P than %d.\n", currP);
               break;
            }
         }
         for (ii = 0; ii < nPaths; ii++)
         {
            dtemp = faPtr->evaluatePoint(moatSample[ii]);
            if (dtemp < threshL || dtemp > threshU)
               printf("moatgen: sample %d fails final test (%e <? %e <? %e).\n",
                      ii, threshL, dtemp, threshU);
         }
         delete [] tempW;
         delete [] indSet;
         if (trial >= nTrials)
         {
            for (ii = 0; ii < nPaths*(nInputs_+1); ii++)
               delete [] moatSample[ii];
            delete [] moatSample;
            continue;
         }
         fp = fopen("MOAT_adjust_file", "w");
         if (fp == NULL)
         {
            printf("ERROR: cannot open file MOAT_adjust_file.\n");
            continue;
         }
         fprintf(fp, "BEGIN\n");
         fprintf(fp, "%d %d\n", nPaths*(nInputs_+1), nInputs_);
         for (ii = 0; ii < nInputs_; ii++)
            fprintf(fp, "%d ", ii+1);
         fprintf(fp, "\n");
         for (ii = 0; ii < nPaths*(nInputs_+1); ii++)
         {
            for (jj = 0; jj < nInputs_; jj++)
               fprintf(fp, "%e ", moatSample[ii][jj]);
            fprintf(fp, "\n");
         }
         fprintf(fp, "END\n");
         fclose(fp);
         count = nPaths * (nInputs_ + 1); 
         tempW = new double[count*nInputs_];
         states = new int[count];
         for (ii = 0; ii < count; ii++) states[ii] = 1;
         for (ii = 0; ii < count; ii++)
            for (jj = 0; jj < nInputs_; jj++)
               tempW[ii*nInputs_+jj] = moatSample[ii][jj];
         printf("moatgen: check for repeated sample points.\n");
         for (ii = 0; ii < count; ii++)
         {
            status = compareSamples(ii,count,nInputs_, tempW, states);
            if (status >= 0)
               printf("moatgen check: sample %d and %d are identical.\n",
                      ii+1,status+1);
         }
         printf("moatgen: adjust file created in MOAT_adjust_file.\n");
         printf("         Make sure to change the input indices to match\n");
         printf("         the indices in the original MOAT file when used\n");
         printf("         with gmoat_adjust.\n");
         for (ii = 0; ii < nPaths*(nInputs_+1); ii++) delete [] moatSample[ii];
         delete [] moatSample;
         delete [] tempW;
         delete [] states;
         delete faPtr;
         faPtr = NULL;
      }

      // +++ moatgen2 
      else if (!strcmp(command, "moatgen2"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("moatgen2: create a MOAT adjust file with multiple constraints\n");
            printf("syntax: moatgen\n");
            printf("Note: a PSUADE datafile should have been loaded before\n");
            printf("      using this command. The data will be used as a\n");
            printf("      response surface to find feasible region in\n");
            printf("      creating a MOAT adjust sample.\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL)
         {
            printf("ERROR: no input file loaded.\n");
            continue;
         }
         status = 0;
         sprintf(pString,"How many constraint data files are there (1-10)? ");
         nFiles = getInt(1, 10, pString);
         faFlag = 2;
         faPtrs = new FuncApprox*[nFiles];
         threshLs = new double[nFiles];
         threshUs = new double[nFiles];
         for (kk = 0; kk < nFiles; kk++)
         {
            sprintf(pString,"Enter name of constraint file #%d : ", kk+1);
            getString(pString, winput);
            kk = strlen(winput);
            winput[kk-1] = '\0';
            ioPtr = new PsuadeData;
            status = ioPtr->readPsuadeFile(winput);
            if (status != 0)
            {
               printf("moatgen2 FILE READ ERROR: file = %s\n", winput);
               continue;
            }
            ioPtr->getParameter("input_ninputs", pPtr);
            jj = pPtr.intData_;
            if (jj != nInputs_)
            {
               printf("moatgen2 ERROR: nInputs mismatch.\n");
               printf("         local nInputs = %d.\n", nInputs_);
               printf("         file  nInputs = %d.\n", jj);
               exit(1);
            }
            pLower.clean();
            ioPtr->getParameter("input_lbounds", pLower);
            for (ii = 0; ii < nInputs_; ii++)
            {
               if (iLowerB_[ii] != pLower.dbleArray_[ii])
               {
                  printf("moatgen2 ERROR: lower bound mismatch (input %d)\n", 
                         ii+1);
                  exit(1);
               }
            }
            pUpper.clean();
            ioPtr->getParameter("input_ubounds", pUpper);
            for (ii = 0; ii < nInputs_; ii++)
            {
               if (iUpperB_[ii] != pUpper.dbleArray_[ii])
               {
                  printf("moatgen2 ERROR: upper bound mismatch (input %d)\n",
                         ii+1);
                  exit(1);
               }
            }
            faPtrs[kk] = genFAInteractive(ioPtr, faFlag);
            if (faPtrs[kk] == NULL) {printf("ERROR detected in RS.\n"); exit(1);}
            faPtrs[kk]->setOutputLevel(outputLevel_);
            sprintf(pString,"Constraint %d lower bound : ",kk+1);
            threshLs[kk] = getDouble(pString);
            sprintf(pString,"Constraint %d upper bound : ",kk+1);
            threshUs[kk] = getDouble(pString);
            if (threshLs[kk] >= threshUs[kk])
            {
               printf("ERROR: lower bound >= upper bound.\n");
               exit(1);;
            }
            delete ioPtr;
         }
         if (status != 0) continue;
         nPaths = 5000;
         sprintf(pString, "Enter P (resolution: try 4-10) : ");
         currP = getInt(4, 10, pString);
         sprintf(pString, "Enter the number of trials (> 100) : ");
         nTrials = getInt(101, 10000000, pString);
         moatSample = new double*[nPaths*(nInputs_+1)];
         for (ii = 0; ii < nPaths*(nInputs_+1); ii++)
            moatSample[ii] = new double[nInputs_];
         tempW = new double[nInputs_];
         indSet = new int[nInputs_];
         count = 0;
         for (ii = 0; ii < nPaths; ii++)
         {
            trial = 0; 
            while (trial < nTrials)
            {
               iInd = count;
               trial++;
               for (jj = 0; jj < nInputs_; jj++)
               {
                  ind = PSUADE_rand() % currP;
                  dtemp = ind * (iUpperB_[jj] - iLowerB_[jj]) / (currP - 1.0);
                  moatSample[iInd][jj] = iLowerB_[jj] + dtemp;
               }
               for (kk = 0; kk < nFiles; kk++)
               {
                  dtemp = faPtrs[kk]->evaluatePoint(moatSample[iInd]);
                  if (dtemp < threshLs[ii] || dtemp > threshUs[ii]) break;
               }
               if (kk != nFiles) continue;

               generateRandomIvector(nInputs_, indSet);

               for (jj = 0; jj < nInputs_; jj++)
               {
                  iInd++;
                  for (kk = 0; kk < nInputs_; kk++)
                     moatSample[iInd][kk] = moatSample[iInd-1][kk];

                  ind2 = indSet[jj];

                  ddata = moatSample[iInd][ind2]; 
                  currX1F = - PSUADE_UNDEFINED;
                  currX2F =   PSUADE_UNDEFINED;
                  for (kk = 0; kk < nFiles; kk++)
                  {
                     filterRange = threshUs[kk] - threshLs[kk];
                     moatSample[iInd][ind2] = iLowerB_[ind2];
                     currY1 = faPtrs[kk]->evaluatePoint(moatSample[iInd]);
                     moatSample[iInd][ind2] = iUpperB_[ind2];
                     currY2 = faPtrs[kk]->evaluatePoint(moatSample[iInd]);
                     currX1 = iLowerB_[ind2];
                     currX2 = iUpperB_[ind2];

                     if (currY2 >= threshUs[kk] && currY1 >= threshUs[kk])
                        currX1 = currX2 = 0.0;
                     else if (currY2 <= threshLs[kk] && currY1 <= threshLs[kk])
                        currX1 = currX2 = 0.0;
                     else if (currY2 > currY1)
                     {
                        if (currY2 <= threshUs[kk]) currX2 = iUpperB_[ind2];
                        else
                        {
                           sLo = iLowerB_[ind2];
                           sHi = iUpperB_[ind2];
                           while (PABS((currY2-threshUs[kk])/filterRange)>1e-4)
                           {
                              moatSample[iInd][ind2] = 0.5 * (sLo + sHi);
                              currY2=faPtrs[kk]->evaluatePoint(moatSample[iInd]);
                              if (currY2 > threshUs[kk]) sHi = 0.5 * (sLo+sHi);
                              else                       sLo = 0.5 * (sLo+sHi);
                           }
                           currX2 = moatSample[iInd][ind2];
                        }
                        if (currY1 >= threshLs[kk]) currX1 = iLowerB_[ind2];
                        else
                        {
                           sLo = iLowerB_[ind2];
                           sHi = iUpperB_[ind2];
                           while (PABS((currY1-threshLs[kk])/filterRange)>1e-4)
                           {
                              moatSample[iInd][ind2] = 0.5 * (sLo + sHi);
                              currY1=faPtrs[kk]->evaluatePoint(moatSample[iInd]);
                              if (currY1 < threshLs[kk]) sLo = 0.5 * (sLo+sHi);
                              else                       sHi = 0.5 * (sLo+sHi);
                           }
                           currX1 = moatSample[iInd][ind2];
                        }
                     }
                     else
                     {
                        if (currY1 <= threshUs[kk]) currX1 = iLowerB_[ind2];
                        else
                        {
                           sLo = iLowerB_[ind2];
                           sHi = iUpperB_[ind2];
                           while (PABS((currY1-threshUs[kk])/filterRange)>1e-4)
                           {
                              moatSample[iInd][ind2] = 0.5 * (sLo + sHi);
                              currY1=faPtrs[kk]->evaluatePoint(moatSample[iInd]);
                              if (currY1 > threshUs[kk]) sLo = 0.5 * (sLo+sHi);
                              else                       sHi = 0.5 * (sLo+sHi);
                           }
                           currX1 = moatSample[iInd][ind2];
                        }
                        if (currY2 >= threshLs[kk]) currX2 = iUpperB_[ind2];
                        else
                        {
                           sLo = iLowerB_[ind2];
                           sHi = iUpperB_[ind2];
                           while (PABS((currY2-threshLs[kk])/filterRange)>1e-4)
                           {
                              moatSample[iInd][ind2] = 0.5 * (sLo + sHi);
                              currY2=faPtrs[kk]->evaluatePoint(moatSample[iInd]);
                              if (currY2 < threshLs[kk]) sHi = 0.5 * (sLo+sHi);
                              else                       sLo = 0.5 * (sLo+sHi);
                           }
                           currX2 = moatSample[iInd][ind2];
                        }
                     }
                     if (PABS(currX2-currX1)<0.1*(iUpperB_[ind2]-iLowerB_[ind2])) 
                        break;
                     if (currX1 > currX1F) currX1F = currX1;
                     if (currX2 < currX2F) currX2F = currX2;
                  }
                  moatSample[iInd][ind2] = ddata;
                  if (kk != nFiles) break;
                  tempW[ind2] = PABS(currX2F - currX1F) / (currP-1.0);
                  moatSample[iInd][ind2] += tempW[ind2];
                  if (moatSample[iInd][ind2] > iUpperB_[ind2])
                       dtemp = threshLs[kk] - 1.0;
                  else
                  {
                     for (kk = 0; kk < nFiles; kk++)
                     {
                        moatSample[iInd][ind2] = ddata;
                        moatSample[iInd][ind2] += tempW[ind2];
                        dtemp = faPtrs[kk]->evaluatePoint(moatSample[iInd]);
                        if (dtemp < threshLs[kk] || dtemp > threshUs[kk]) 
                        {
                           moatSample[iInd][ind2] -= 2.0 * tempW[ind2];
                           if (moatSample[iInd][ind2] < iLowerB_[ind2])
                              break;
                           dtemp = faPtrs[kk]->evaluatePoint(moatSample[iInd]);
                           if (dtemp < threshLs[kk] || dtemp > threshUs[kk]) 
                              break;
                        }
                     }
                     moatSample[iInd][ind2] = ddata;
                     if (kk != nFiles) break;
                  }
               }
               if (jj == nInputs_) 
               {
                  count += (nInputs_ + 1);
                  if (outputLevel_ > 2)
                     printf("moatgen2: path %d (out of %d) found.\n", ii+1,
                            nPaths);
                  break; 
               }
               else
               {
                  if (outputLevel_ > 2)
                     printf("Current path fails (%d out of max %d).\n",
                            trial, nTrials); 
               }
            }
            if (trial >= nTrials)
            {
               printf("moatgen2 fails to find all possible paths.\n");
               printf("Suggestion: try again with a larger P than %d.\n", currP);
               break;
            }
         }
         for (ii = 0; ii < nPaths; ii++)
         {
            for (kk = 0; kk < nFiles; kk++)
            {
               dtemp = faPtrs[kk]->evaluatePoint(moatSample[ii]);
               if (dtemp < threshLs[kk] || dtemp > threshUs[kk])
               printf("moatgen2: sample %d fails final test (%e <? %e <? %e).\n",
                      ii, threshLs[kk], dtemp, threshUs[kk]);
            }
         }
         delete [] tempW;
         delete [] indSet;
         delete [] threshLs;
         delete [] threshUs;
         for (kk = 0; kk < nFiles; kk++) delete faPtrs[kk];
         delete [] faPtrs;
         faPtrs = NULL;
         if (trial >= nTrials)
         {
            for (ii = 0; ii < nPaths*(nInputs_+1); ii++)
               delete [] moatSample[ii];
            delete [] moatSample;
            continue;
         }
         fp = fopen("MOAT_adjust_file", "w");
         if (fp == NULL)
         {
            printf("ERROR: cannot open file MOAT_adjust_file.\n");
            continue;
         }
         fprintf(fp, "BEGIN\n");
         fprintf(fp, "%d %d\n", nPaths*(nInputs_+1), nInputs_);
         for (ii = 0; ii < nInputs_; ii++)
            fprintf(fp, "%d ", ii+1);
         fprintf(fp, "\n");
         for (ii = 0; ii < nPaths*(nInputs_+1); ii++)
         {
            for (jj = 0; jj < nInputs_; jj++)
               fprintf(fp, "%e ", moatSample[ii][jj]);
            fprintf(fp, "\n");
         }
         fprintf(fp, "END\n");
         fclose(fp);
         count = nPaths * (nInputs_ + 1); 
         tempW = new double[count*nInputs_];
         states = new int[count];
         for (ii = 0; ii < count; ii++) states[ii] = 1;
         for (ii = 0; ii < count; ii++)
            for (jj = 0; jj < nInputs_; jj++)
               tempW[ii*nInputs_+jj] = moatSample[ii][jj];
         printf("moatgen2: check for repeated sample points.\n");
         for (ii = 0; ii < count; ii++)
         {
            status = compareSamples(ii,count,nInputs_, tempW, states);
            if (status >= 0)
               printf("moatgen2 check: sample %d and %d are identical.\n",
                      ii+1,status+1);
         }
         printf("moatgen2: adjust file created in MOAT_adjust_file.\n");
         printf("          Make sure to change the input indices to match\n");
         printf("          the indices in the original MOAT file when used\n");
         printf("          with gmoat_adjust.\n");
         for (ii = 0; ii < nPaths*(nInputs_+1); ii++) delete [] moatSample[ii];
         delete [] moatSample;
         delete [] tempW;
         delete [] states;
      }

      // +++ moat_concatenate 
      else if (!strcmp(command, "moat_concat"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("moat_concat: concatenate 2 MOAT samples (different inputs)\n");
            printf("syntax: moat_concat <file>\n");
            printf("Note: a PSUADE MOAT datafile should have been loaded\n");
            printf("      before using this command. \n");
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL)
         {
            printf("ERROR: no input file loaded.\n");
            continue;
         }
         strcpy(dataFile, "\0");
         sscanf(lineIn,"%s %s",command,dataFile);
         if ((fp=fopen(dataFile,"r")) == NULL)
         {
            printf("file %s not found.\n", dataFile);
            printf("syntax: moat_concat <file>.\n");
            printf("where <file> is a PSUADE data file.\n");
            continue;
         }
         psuadeIO_->getParameter("method_sampling", pPtr);
         if (pPtr.intData_ != PSUADE_SAMP_MOAT)
         {
            printf("ERROR: local sample is not MOAT. \n");
            continue;
         }
         ioPtr = new PsuadeData();
         ioPtr->setOutputLevel(0);
         status = ioPtr->readPsuadeFile(dataFile);
         if (status == 0)
         {
            ioPtr->getParameter("method_sampling", pPtr);
            if (pPtr.intData_ != PSUADE_SAMP_MOAT)
            {
               printf("ERROR: second sample is not MOAT. \n");
               continue;
            }
            nReps = nSamples_ / (nInputs_ + 1);
            ioPtr->getParameter("input_ninputs", pPtr);
            ind = pPtr.intData_;
            ioPtr->getParameter("method_nsamples", pPtr);
            count = pPtr.intData_;
            ll = count / (ind + 1);
            if (nReps != ll)
            {
               printf("ERROR: different number of replications.\n");
               printf("       num_replcations for sample 1 = %d\n",nReps);
               printf("       num_replcations for sample 2 = %d\n",ll);
               continue;
            }
            names = inputNames_;
            inputNames_ = new char*[nInputs_+ind];
            for (ii = 0; ii < nInputs_; ii++) inputNames_[ii] = names[ii];
            pINames.clean();
            psuadeIO_->getParameter("input_names", pINames);
            names = pINames.strArray_;
            for (ii = nInputs_; ii < nInputs_+ind; ii++)
            {
               inputNames_[ii] = new char[200];
               strcpy(inputNames_[ii], names[ii-nInputs_]);
            }
            tempX = sampleInputs_;
            if (sampleOutputs_ != NULL) delete [] sampleOutputs_;
            if (sampleStates_  != NULL) delete [] sampleStates_;
            kk = nReps * (nInputs_ + ind + 1) * (nInputs_ + ind);
            sampleInputs_ = new double[kk];
            kk = nReps * (nInputs_ + ind + 1) * nOutputs_;
            sampleOutputs_ = new double[kk];
            for (ii = 0; ii < kk; ii++) sampleOutputs_[ii] = PSUADE_UNDEFINED;
            kk = nReps * (nInputs_ + ind + 1);
            sampleStates_ = new int[kk];
            for (ii = 0; ii < kk; ii++) sampleStates_[ii] = 0;
            for (ii = 0; ii < nReps; ii++)
            {
               for (jj = 0; jj <= nInputs_; jj++)
               {
                  for (kk = 0; kk < nInputs_; kk++)
                  {
                     ind2 = ii * (nInputs_ + ind + 1) * (nInputs_ + ind);
                     sampleInputs_[ind2+jj*(nInputs_+ind)+kk] = 
                            tempX[ii*(nInputs_+1)*nInputs_+jj*nInputs_+kk];
                  }
               }
               for (jj = nInputs_+1; jj < nInputs_+ind+1; jj++)
               {
                  for (kk = 0; kk < nInputs_; kk++)
                  {
                     ind2 = ii * (nInputs_ + ind + 1) * (nInputs_ + ind);
                     sampleInputs_[ind2+jj*(nInputs_+ind)+kk] = 
                           sampleInputs_[ind2+nInputs_*(nInputs_+ind)+kk]; 
                  }
               }
            }
            if (tempX != NULL) delete [] tempX;
            ioPtr->getParameter("input_sample", pPtr);
            tempX = pPtr.dbleArray_;
            pPtr.dbleArray_ = NULL;
            for (ii = 0; ii < nReps; ii++)
            {
               for (jj = 0; jj <= nInputs_; jj++)
               {
                  for (kk = nInputs_; kk < nInputs_+ind; kk++)
                  {
                     ind2 = ii * (nInputs_ + ind + 1) * (nInputs_ + ind);
                     sampleInputs_[ind2+jj*(nInputs_+ind)+kk] = 
                            tempX[ii*(ind+1)*ind+kk-nInputs_];
                  }
               }
               for (jj = nInputs_+1; jj < nInputs_+ind+1; jj++)
               {
                  for (kk = nInputs_; kk < nInputs_+ind; kk++)
                  {
                     ind2 = ii * (nInputs_ + ind + 1) * (nInputs_ + ind);
                     sampleInputs_[ind2+jj*(nInputs_+ind)+kk] = 
                            tempX[ii*(ind+1)*ind+(jj-nInputs_)*ind+kk-nInputs_];
                  }
               }
            }
            delete [] tempX;
            pLower.clean();
            ioPtr->getParameter("input_lbounds", pLower);
            tempW = pLower.dbleArray_;
            tempV = iLowerB_;
            iLowerB_ = new double[nInputs_+ind];
            for (ii = 0; ii < nInputs_; ii++) iLowerB_[ii] = tempV[ii];
            for (ii = nInputs_; ii < nInputs_+ind; ii++)
               iLowerB_[ii] = tempW[ii-nInputs_];
            if (tempV != NULL) delete [] tempV;
            pLower.clean();
            pUpper.clean();
            ioPtr->getParameter("input_ubounds", pUpper);
            tempW = pUpper.dbleArray_;
            tempV = iUpperB_;
            iUpperB_ = new double[nInputs_+ind];
            for (ii = 0; ii < nInputs_; ii++) iUpperB_[ii] = tempV[ii];
            for (ii = nInputs_; ii < nInputs_+ind; ii++)
               iUpperB_[ii] = tempW[ii-nInputs_];
            if (tempV != NULL) delete [] tempV;
            pUpper.clean();
            delete ioPtr;
            nSamples_ = (nInputs_ + ind + 1) * nReps;;
            tempI = inputPDFs_;
            inputPDFs_ = new int[nInputs_];
            for (ii = 0; ii < nInputs_; ii++) inputPDFs_[ii] = tempI[ii];
            for (ii = nInputs_; ii < nInputs_+ind; ii++) inputPDFs_[ii] = 0;
            tempW = inputMeans_;
            inputMeans_ = new double[nInputs_];
            for (ii = 0; ii < nInputs_; ii++) inputMeans_[ii] = tempW[ii];
            for (ii = nInputs_; ii < nInputs_+ind; ii++) inputMeans_[ii] = 0;
            tempW = inputStds_;
            inputStds_ = new double[nInputs_];
            for (ii = 0; ii < nInputs_; ii++) inputStds_[ii] = tempW[ii];
            for (ii = nInputs_; ii < nInputs_+ind; ii++) inputStds_[ii] = 0;
            Matrix *tmpMat = new Matrix();
            tmpMat->setDim(nInputs_+ind,nInputs_+ind);
            for (ii = 0; ii < nInputs_; ii++)
            {
               for (jj = 0; jj < nInputs_; jj++)
               {
                  ddata = inputCMat_->getEntry(ii,jj);
                  tmpMat->setEntry(ii,jj,ddata);
               }
            }
            for (ii = nInputs_; ii < nInputs_+ind; ii++)
               tmpMat->setEntry(ii,ii,1.0);
            delete inputCMat_;
            inputCMat_ = tmpMat;
            nInputs_ += ind;

            psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,iLowerB_,
                            iUpperB_,sampleInputs_,inputNames_, 
                            inputPDFs_, inputMeans_,inputStds_,inputCMat_); 
            psuadeIO_->updateOutputSection(nSamples_,nOutputs_,
                                sampleOutputs_, sampleStates_, NULL); 
            psuadeIO_->updateMethodSection(-1, nSamples_, 1, -1, -1);
            if (currSession != NULL) delete currSession;
            currSession = new PsuadeSession();
            psuadeIO_->getSession(currSession);
            printf("The two samples have been concatenated.\n");
            printf("The new sample has nInputs = %d\n", nInputs_);
            printf("                  nSamples = %d\n", nSamples_);
            printf("                  nOutputs = %d\n", nOutputs_);
            printf("Use 'write' to write the expanded sample to a file.\n");
         }
         else
         {
            printf("ERROR reading second sample file.\n");
         }
      }

      // +++ splitsample 
      else if (!strcmp(command, "splitsample"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("splitsample: split the existing data set into two subsets.\n");
            printf("Note: the two sets of data will be stored in psuadeSplit1\n");
            printf("      and psuadeSplit2.\n");
            continue;
         }
         if (psuadeIO_ == NULL || sampleOutputs_ == NULL)
         {
            printf("ERROR: no sample output data to split.\n");
            printf("       Use load to put data into local memory first.\n");
         }
         else if (nSamples_ <= 0) printf("Reload data file.\n");
         else
         {
            printf("The current sample size is %d.\n", nSamples_);
            sprintf(pString,
                    "Sample size of the first set? (1 - %d) ",nSamples_-1);
            kk = getInt(1, nSamples_, pString);
            sprintf(pString, "Random draw from original sample ? (y or n) ");
            getString(pString, winput);
            tempX  = new double[kk*nInputs_];
            tempY  = new double[kk*nOutputs_];
            states = new int[kk];
            tags   = new int[nSamples_];
            for (ii = 0; ii < nSamples_; ii++) tags[ii] = 0;
            for (ii = 0; ii < kk; ii++)
            {
               if (winput[0] == 'y')
               {
                  ind = PSUADE_rand() % nSamples_;
                  ind2 = 0;
                  while (tags[ind2] == 1 && ind2 < 1000)
                  {
                     ind = PSUADE_rand() % nSamples_;
                     ind2++;
                  }
                  if (tags[ind] == 1)
                     for (ind = 0; ind < nSamples_; ind++)
                        if (tags[ind] == 0) break;
                  if (tags[ind] == 1)
                  {
                     printf("ERROR: cannot split data. \n");
                     continue;
                  }
               } else ind = ii;
               for (jj = 0; jj < nInputs_; jj++)
                  tempX[ii*nInputs_+jj] = sampleInputs_[ind*nInputs_+jj];
               for (jj = 0; jj < nOutputs_; jj++)
                  tempY[ii*nOutputs_+jj] = sampleOutputs_[ind*nOutputs_+jj];
               tags[ind] = 1;
               states[ii] = sampleStates_[ind];
            }
            ind = 0;
            for (ii = 0; ii < nSamples_; ii++)
            {
               if (tags[ii] == 0)
               {
                  for (jj = 0; jj < nInputs_; jj++)
                     sampleInputs_[ind*nInputs_+jj]=sampleInputs_[ii*nInputs_+jj];
                  for (jj = 0; jj < nOutputs_; jj++)
                     sampleOutputs_[ind*nOutputs_+jj] = 
                                   sampleOutputs_[ii*nOutputs_+jj];
                  sampleStates_[ind] = sampleStates_[ii]; 
                  ind++;
               }
            }
            psuadeIO_->updateInputSection(kk,nInputs_,NULL,NULL,NULL,
                                    tempX,NULL, NULL, NULL, NULL, NULL);
            psuadeIO_->updateOutputSection(kk,nOutputs_,
                                           tempY,states,outputNames_); 
            psuadeIO_->updateMethodSection(-1,kk,-1,-1,-1);
            psuadeIO_->writePsuadeFile("psuadeSplit1",0);

            psuadeIO_->updateInputSection(nSamples_-kk,nInputs_,NULL,NULL,
                                      NULL,sampleInputs_,NULL, 
                                      NULL, NULL, NULL, NULL);
            psuadeIO_->updateOutputSection(nSamples_-kk,nOutputs_,
                             sampleOutputs_,sampleStates_,outputNames_); 
            psuadeIO_->updateMethodSection(-1,nSamples_-kk,-1,-1,-1);
            psuadeIO_->writePsuadeFile("psuadeSplit2",0);
            if (currSession != NULL) delete currSession;
            currSession = NULL;

            printf("The 2 data files are in psuadeSplit1 and psuadeSplit2.\n");
            printf("The loaded data in local memory have been destroyed.\n");
            fflush(stdout);
            delete [] tempX;
            delete [] tempY;
            delete [] states;
            if (inputNames_ != NULL)
            {
               for (ii = 0; ii < nInputs_; ii++) delete [] inputNames_[ii];
               delete [] inputNames_;
            }
            if (outputNames_ != NULL)
            {
               for (ii = 0; ii < nOutputs_; ii++) delete [] outputNames_[ii];
               delete [] outputNames_;
            }
            inputNames_  = NULL;
            outputNames_ = NULL;
            if (sampleInputs_  != NULL) delete [] sampleInputs_;
            if (sampleOutputs_ != NULL) delete [] sampleOutputs_;
            if (sampleStates_  != NULL) delete [] sampleStates_;
            if (iLowerB_       != NULL) delete [] iLowerB_;
            if (iUpperB_       != NULL) delete [] iUpperB_;
            if (inputPDFs_     != NULL) delete [] inputPDFs_;
            if (inputMeans_    != NULL) delete [] inputMeans_;
            if (inputStds_     != NULL) delete [] inputStds_;
            sampleInputs_  = NULL;
            sampleOutputs_ = NULL;
            sampleStates_  = NULL;
            iLowerB_       = NULL;
            iUpperB_       = NULL;
            inputPDFs_     = NULL;
            inputMeans_    = NULL;
            inputStds_     = NULL;
            nSamples_ = 0;
            nInputs_ = 0;
            nOutputs_ = 0;
         }
      }

      // UQ/SA commands
      // +++ uq 
      else if (!strcmp(command, "ua"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("ua: uncertainty analysis (compute moments).\n");
            printf("syntax: ua (after data have been loaded)\n");
            continue;
         }
         if (psuadeIO_ == NULL || sampleOutputs_ == NULL)
         {
            printf("ERROR: no data to analyze (load data first).\n");
            continue;
         }
         outputID = 0;
         sprintf(pString, "Enter output number (1 - %d) : ", nOutputs_);
         outputID = getInt(0, nOutputs_, pString);
         outputID--;
         analysisMethod = PSUADE_ANA_MOMENT;
         anaManager = new AnalysisManager();
         anaManager->setup(analysisMethod, 0);
         psuadeIO_->getParameter("ana_diagnostics",pPtr);
         ii = pPtr.intData_;
         psuadeIO_->updateAnalysisSection(-1,-1,-1,outputLevel_,-1, -1);
         anaManager->analyze(psuadeIO_, 0, NULL, outputID);
         psuadeIO_->updateAnalysisSection(-1,-1,-1,ii,-1, -1);
         delete anaManager;
      }

      // +++ ca 
      else if (!strcmp(command, "ca"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("ca: correlation analysis\n");
            printf("syntax: ca (after data have been loaded)\n");
            continue;
         }
         if (psuadeIO_ == NULL || sampleOutputs_ == NULL)
         {
            printf("ERROR: no data to analyze (load data first).\n");
            continue;
         }
         outputID = 0;
         sprintf(pString, "Enter output number (1 - %d) : ", nOutputs_);
         outputID = getInt(1, nOutputs_, pString);
         outputID--;
         analysisMethod = PSUADE_ANA_CORRELATION;
         anaManager = new AnalysisManager();
         anaManager->setup(analysisMethod, 0);
         psuadeIO_->getParameter("ana_diagnostics",pPtr);
         ii = pPtr.intData_;
         psuadeIO_->updateAnalysisSection(-1,-1,-1,outputLevel_,-1, -1);
         anaManager->analyze(psuadeIO_, 0, NULL, outputID);
         psuadeIO_->updateAnalysisSection(-1,-1,-1,ii,-1, -1);
         delete anaManager;
      }

      // +++ anova 
      else if (!strcmp(command, "anova"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("anova: analysis of variation\n");
            printf("syntax: anova (after data have been loaded)\n");
            continue;
         }
         if (psuadeIO_ == NULL || sampleOutputs_ == NULL)
         {
            printf("ERROR: no data to analyze (load data first).\n");
            continue;
         }
         outputID = 0;
         sprintf(pString, "Enter output number (1 - %d) : ", nOutputs_);
         outputID = getInt(1, nOutputs_, pString);
         outputID--;
         analysisMethod = PSUADE_ANA_ANOVA;
         anaManager = new AnalysisManager();
         anaManager->setup(analysisMethod, 0);
         psuadeIO_->getParameter("ana_diagnostics",pPtr);
         ii = pPtr.intData_;
         psuadeIO_->updateAnalysisSection(-1,-1,-1,outputLevel_,-1, -1);
         anaManager->analyze(psuadeIO_, 0, NULL, outputID);
         psuadeIO_->updateAnalysisSection(-1,-1,-1,ii,-1, -1);
         delete anaManager;
      }

      // +++ moat 
      else if (!strcmp(command, "moat"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("moat: Morris screening analysis\n");
            printf("syntax: moat (after data have been loaded)\n");
            continue;
         }
         if (psuadeIO_ == NULL || sampleOutputs_ == NULL)
         {
            printf("ERROR: no data to analyze (load data first).\n");
            continue;
         }
         outputID = 0;
         sprintf(pString, "Enter output number (1 - %d) = ", nOutputs_);
         outputID = getInt(1, nOutputs_, pString);
         outputID--;
         analysisMethod = PSUADE_ANA_MOAT;
         anaManager = new AnalysisManager();
         anaManager->setup(analysisMethod, 0);
         psuadeIO_->getParameter("ana_diagnostics",pPtr);
         ii = pPtr.intData_;
         psuadeIO_->updateAnalysisSection(-1,-1,-1,outputLevel_,-1, -1);
         anaManager->analyze(psuadeIO_, 0, NULL, outputID);
         psuadeIO_->updateAnalysisSection(-1,-1,-1,ii,-1, -1);
         delete anaManager;
      }

      // +++ moatmo 
      else if (!strcmp(command, "moatmo"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("moatmo: Morris screening analysis for multiple\n");
            printf("        outputs simultaneously\n");
            printf("syntax: moatmo (after data have been loaded)\n");
            continue;
         }
         if (psuadeIO_ == NULL || sampleOutputs_ == NULL)
         {
            printf("ERROR: no data to analyze (load data first).\n");
            continue;
         }
         if (nOutputs_ <= 1)
         {
            printf("INFO: only one output -> use moat instead.\n");
            continue;
         } 
         psuadeIO_->getParameter("ana_diagnostics",pPtr);
         int diagSave = pPtr.intData_;
         psuadeIO_->updateAnalysisSection(-1,-1,-1,outputLevel_,-1, -1);
         analysisMethod = PSUADE_ANA_MOAT;
         anaManager = new AnalysisManager();
         anaManager->setup(analysisMethod, 0);
         double *tempW = new double[nInputs_*nOutputs_];
         int anaExpertSave = psAnaExpertMode_;
         if (psAnaExpertMode_ == 1)
         {
            printf("INFO: Analysis expert mode will be turned off.\n");
            printf("      This command will only create modified means\n");
            printf("      plots. If you want to create screen and scatter\n");
            printf("      plots, you will have to use 'moat' for each\n");
            printf("      output.\n");
            psAnaExpertMode_ = 0;
         }
         for (ii = 0; ii < nOutputs_; ii++)
         {
            anaManager->analyze(psuadeIO_, 0, NULL, ii);
            pdata = psuadeIO_->getAuxData();
            if (pdata->nDbles_ == nInputs_)
            {
               for (jj = 0; jj < nInputs_; jj++) 
                  tempW[ii*nInputs_+jj] = pdata->dbleArray_[jj];
               pdata->clean();
            }
         }
         psuadeIO_->updateAnalysisSection(-1,-1,-1,diagSave,-1, -1);
         psAnaExpertMode_ = anaExpertSave;
         delete anaManager;
         fp = NULL;
         if (psPlotTool_ == 0) fp = fopen("matlabmoatmo.m", "w");
         if (fp != NULL)
         {
            sprintf(pString,"This file contains Morris modified Means");
            fwriteComment(fp, pString);
            fprintf(fp, "A = [\n");
            for (ii = 0; ii < nInputs_; ii++)
            {
               for (jj = 0; jj < nOutputs_; jj++)
                  fprintf(fp,"%16.8e ", tempW[jj*nInputs_+ii]);
               fprintf(fp, "\n");
            }
            fprintf(fp, "];\n");
            fprintf(fp, "A2 = A * inv(diag(max(A)));\n");
            if (inputNames_ == NULL)
            {
               fprintf(fp, "  Str = {");
               for (ii = 0; ii < nInputs_-1; ii++) fprintf(fp,"'X%d',",ii+1);
               fprintf(fp,"'X%d'};\n",nInputs_);
            }
            else
            {
               fprintf(fp, "  XStr = {");
               for (ii = 0; ii < nInputs_-1; ii++)
               {
                  if (inputNames_[ii] != NULL)
                       fprintf(fp,"'%s',",inputNames_[ii]);
                  else fprintf(fp,"'X%d',",ii+1);
               }
               if (inputNames_[nInputs_-1] != NULL)
                    fprintf(fp,"'%s'};\n",inputNames_[nInputs_-1]);
               else fprintf(fp,"'X%d'};\n",nInputs_);
            }
            fprintf(fp, "  YStr = {");
            for (ii = 0; ii < nOutputs_-1; ii++) fprintf(fp,"'%d',",ii+1);
            fprintf(fp,"'%d'};\n",nOutputs_);
            fwriteHold(fp, 0);
            fprintf(fp,"nn = %d;\n", nInputs_);
            fprintf(fp,"mm = %d;\n", nOutputs_);
            fprintf(fp,"X = 0.5 : 1 : nn-0.5;\n");
            fprintf(fp,"Y = 0.5 : 1 : mm-0.5;\n");
            fprintf(fp, "imagesc(X,Y,A2');\n");
            fprintf(fp, "axis([0 nn 0 mm])\n");
            fwriteHold(fp, 1);
            fprintf(fp,"XX = 0 : nn;\n");
            fprintf(fp,"for ii = 1 : mm-1\n");
            fprintf(fp,"  plot(XX,ii*ones(nn+1,1),'linewidth',2)\n");
            fprintf(fp,"end;\n");
            fprintf(fp,"YY = 0 : mm;\n");
            fprintf(fp,"for ii = 1 : nn-1\n");
            fprintf(fp,"  plot(ii*ones(mm+1,1),YY,'linewidth',2)\n");
            fprintf(fp,"end;\n");
            fprintf(fp,"set(gca,'XTickLabel',[]);\n");
            fprintf(fp,"set(gca,'YTickLabel',[]);\n");
            fprintf(fp,"th=text((1:nn)-0.5,repmat(mm+0.05*mm,nn,1),XStr,");
            fprintf(fp,"'HorizontalAlignment','left','rotation',90);\n");
            fprintf(fp,"set(th, 'fontsize', 12)\n");
            fprintf(fp,"set(th, 'fontweight', 'bold')\n");
            fprintf(fp,"th=text(repmat(-0.05,mm,1),(1:mm)-0.5,YStr,");
            fprintf(fp,"'HorizontalAlignment','left','rotation',90);\n");
            fprintf(fp,"set(th, 'fontsize', 12)\n");
            fprintf(fp,"set(th, 'fontweight', 'bold')\n");
            fprintf(fp,"set(gca, 'fontsize', 12)\n");
            fprintf(fp,"set(gca, 'fontweight', 'bold')\n");
            fprintf(fp,"set(gca, 'linewidth',2)\n");
            fprintf(fp,"box on\n");
            fprintf(fp,"colorbar\n");
            fwritePlotTitle(fp,"Morris Relative Importance Measure");
            fwritePlotXLabel(fp, "Inputs");
            fwritePlotYLabel(fp, "Outputs");
            fprintf(fp,"disp('The colors denote relative magnitudes')\n");
            fwriteHold(fp, 0);
            fclose(fp);
            printf("Morris plot file = matlabmoatmo.m\n");
         }
         delete [] tempW;
      }

      // +++ ff 
      else if (!strcmp(command, "ff"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("ff: Fractional factorial screening analysis\n");
            printf("syntax: ff (after data have been loaded)\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data to analyze (load data first).\n");
            continue;
         }
         outputID = 0;
         sprintf(pString, "Enter output number (1 - %d) = ", nOutputs_);
         outputID = getInt(1, nOutputs_, pString);
         outputID--;
         analysisMethod = PSUADE_ANA_FF;
         anaManager = new AnalysisManager();
         anaManager->setup(analysisMethod, 0);
         psuadeIO_->getParameter("ana_diagnostics",pPtr);
         ii = pPtr.intData_;
         psuadeIO_->updateAnalysisSection(-1,-1,-1,outputLevel_,-1, -1);
         anaManager->analyze(psuadeIO_, 0, NULL, outputID);
         psuadeIO_->updateAnalysisSection(-1,-1,-1,ii,-1, -1);
         delete anaManager;
      }

      // +++ lsa 
      else if (!strcmp(command, "lsa"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("lsa: local sensitivity analysis\n");
            printf("syntax: lsa (after data have been loaded)\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data to analyze (load data first).\n");
            continue;
         }
         outputID = 0;
         sprintf(pString, "Enter output number (1 - %d) = ", nOutputs_);
         outputID = getInt(1, nOutputs_, pString);
         outputID--;
         analysisMethod = PSUADE_ANA_LSA;
         anaManager = new AnalysisManager();
         anaManager->setup(analysisMethod, 0);
         psuadeIO_->getParameter("ana_diagnostics",pPtr);
         ii = pPtr.intData_;
         psuadeIO_->updateAnalysisSection(-1,-1,-1,outputLevel_,-1, -1);
         anaManager->analyze(psuadeIO_, 0, NULL, outputID);
         psuadeIO_->updateAnalysisSection(-1,-1,-1,ii,-1, -1);
         delete anaManager;
      }

      // +++ mars_sa 
      else if (!strcmp(command, "mars_sa"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("mars_sa: MARS-based sensitivity analysis\n");
            printf("syntax: mars_sa (after data have been loaded)\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data to analyze (load data first).\n");
            continue;
         }
         outputID = 0;
         sprintf(pString, "Enter output number (1 - %d) = ", nOutputs_);
         outputID = getInt(1, nOutputs_, pString);
         outputID--;
         sprintf(pString,"MARS (0) or MARS with bagging (1) ? ");
         kk = getInt(0, 1, pString);
         if (kk == 0) faType = PSUADE_RS_MARS;
         else         faType = PSUADE_RS_MARSB;
         faPtr = genFA(faType, nInputs_, iOne, nSamples_);
         if (faPtr != NULL)
         {
            faPtr->setBounds(iLowerB_, iUpperB_);
            faPtr->setOutputLevel(outputLevel_);
            if (faType == PSUADE_RS_MARSB) faPtr->setOutputLevel(4);
            tempY = new double[nSamples_];
            for (ii = 0; ii < nSamples_; ii++)
               tempY[ii] = sampleOutputs_[ii*nOutputs_+outputID];
            status = faPtr->initialize(sampleInputs_,tempY);
            if (faType == PSUADE_RS_MARS)
            {
               strcpy(pString, "rank");
               targv[0] = (char *) pString;
               faPtr->setParams(1, targv);
            }
            delete faPtr;
            faPtr = NULL;
            delete [] tempY;
            tempY = NULL;
         }
      }

      // +++ gp_sa 
      else if (!strcmp(command, "gp_sa"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("gp_sa: Gaussian Process-based sensitivity analysis\n");
            printf("syntax: gp_sa (after data have been loaded)\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data to analyze (load data first).\n");
            continue;
         }
         outputID = 0;
         sprintf(pString, "Enter output number (1 - %d) = ", nOutputs_);
         outputID = getInt(1, nOutputs_, pString);
         outputID--;
         faType = -1;
         printf("Which Gaussian process ? \n");
#ifdef HAVE_TPROS
         printf("1. MacKay's Tpros\n");
#endif
#ifdef HAVE_GPMC
         printf("2. Rasmussen's gp-mc\n");
#endif
         printf("3. Kriging\n");
         sprintf(pString, "Enter number (1, 2, or 3) = ");
         faType = getInt(1, 3, pString);
#ifdef HAVE_TPROS
         if (faType == 1) faType = PSUADE_RS_GP1;
#endif
#ifdef HAVE_GPMC
         if (faType == 2) faType = PSUADE_RS_GP2;
#endif
         if (faType == 3) faType = PSUADE_RS_KR;
         faPtr = genFA(faType, nInputs_, iOne, nSamples_);
         if (faPtr != NULL)
         {
            faPtr->setBounds(iLowerB_, iUpperB_);
            faPtr->setOutputLevel(outputLevel_);
            int rsKeep = psRSExpertMode_;
            psRSExpertMode_ = 0;
            tempY = new double[nSamples_];
            for (ii = 0; ii < nSamples_; ii++)
               tempY[ii] = sampleOutputs_[ii*nOutputs_+outputID];
            status = faPtr->initialize(sampleInputs_,tempY);
            strcpy(pString, "rank");
            targv[0] = (char *) pString;
            faPtr->setParams(1, targv);
            psRSExpertMode_ = rsKeep;
            delete faPtr;
            delete [] tempY;
            faPtr = NULL;
            tempY = NULL;
         }
      }

      // +++ sot_sa 
      else if (!strcmp(command, "sot_sa"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("sot_sa: Sum-of-trees-based sensitivity analysis\n");
            printf("syntax: sot_sa (after data have been loaded)\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data to analyze (load data first).\n");
            continue;
         }
         outputID = 0;
         sprintf(pString, "Enter output number (1 - %d) = ", nOutputs_);
         outputID = getInt(1, nOutputs_, pString);
         outputID--;
         faType = PSUADE_RS_SOTS;
         faPtr  = genFA(faType, nInputs_, iOne, nSamples_);
         if (faPtr == NULL)
         {
            printf("ERROR: cannot create response surface.\n");
            continue;
         }
         faPtr->setBounds(iLowerB_, iUpperB_);
         faPtr->setOutputLevel(outputLevel_);
         tempY = new double[nSamples_];
         for (ii = 0; ii < nSamples_; ii++)
            tempY[ii] = sampleOutputs_[ii*nOutputs_+outputID];
         status = faPtr->initialize(sampleInputs_,tempY);
         strcpy(pString, "mode0");
         targv[0] = (char *) pString;
         faPtr->setParams(1, targv);
         strcpy(pString, "rank");
         targv[0] = (char *) pString;
         ddata = faPtr->setParams(1, targv);
         delete faPtr;
         faPtr = NULL;
         delete [] tempY;
         tempY = NULL;
         printf("sot_sa score (sum of all std dev) = %e\n", ddata);
      }

      // +++ me 
      else if (!strcmp(command, "me"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("me: main effect analysis (variance-based)\n");
            printf("syntax: me (after data have been loaded)\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data to analyze (load data first).\n");
            continue;
         }
         printf("INFO: This command operates on the raw sample data.\n");
         printf("      It is intended for few parameters (<10) with a\n");
         printf("      large sample size (thousands to tens of thousands).\n");
         printf("      The alternative is to use rssobol1 or rssobol1b.\n");
         outputID = 0;
         sprintf(pString, "Enter output number (1 - %d) : ", nOutputs_);
         outputID = getInt(1, nOutputs_, pString);
         outputID--;
         analysisMethod = PSUADE_ANA_ME;
         anaManager = new AnalysisManager();
         anaManager->setup(analysisMethod, 0);
         psuadeIO_->getParameter("ana_diagnostics",pPtr);
         ii = pPtr.intData_;
         psuadeIO_->updateAnalysisSection(-1,-1,-1,outputLevel_,-1, -1);
         anaManager->analyze(psuadeIO_, 0, NULL, outputID);
         psuadeIO_->updateAnalysisSection(-1,-1,-1,ii,-1, -1);
         pdata = psuadeIO_->getAuxData(); 
         pdata->clean();
         delete anaManager;
      }

      // +++ ie 
      else if (!strcmp(command, "ie"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("ie: 2-way interaction effect analysis (variance-based)\n");
            printf("syntax: ie (after data have been loaded)\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data to analyze (load data first).\n");
            continue;
         }
         if (nInputs_ <= 2)
         {
            printf("INFO: there is no point doing this for nInputs <= 2\n");
            printf("      For nInputs=2, interaction effect = total variance.\n");
            continue;
         }
         outputID = 0;
         printf("INFO: This command operates on the raw sample data.\n");
         printf("      It is intended for few parameters (<<10) with a\n");
         printf("      large sample size (thousands to tens of thousands).\n");
         printf("      The alternative is to use rssobol2 or rssobol2b.\n");
         sprintf(pString, "Enter output number (1 - %d) : ", nOutputs_);
         outputID = getInt(1, nOutputs_, pString);
         outputID--;
         analysisMethod = PSUADE_ANA_IE;
         anaManager = new AnalysisManager();
         anaManager->setup(analysisMethod, 0);
         psuadeIO_->getParameter("ana_diagnostics",pPtr);
         ii = pPtr.intData_;
         psuadeIO_->updateAnalysisSection(-1,-1,-1,outputLevel_,-1,-1);
         anaManager->analyze(psuadeIO_, 0, NULL, outputID);
         psuadeIO_->updateAnalysisSection(-1,-1,-1,ii,-1,-1);
         delete anaManager;
         pdata = psuadeIO_->getAuxData();
         if (pdata->nDbles_ >= nInputs_)
         {
            printEquals(PL_INFO, 0);
            printf("Pairwise Interaction (+Main) Effect Statistics: \n");
            if (pdata->dbleData_ > 0)
            {
               for (ii = 0; ii < nInputs_; ii++)
                  for (jj = ii+1; jj < nInputs_; jj++)
                     printf("Inputs %4d %4d: sensitivity index = %11.4e ",
                            ii+1,jj+1,pdata->dbleArray_[ii*nInputs_+jj]);
                     printf("(normalized=%11.4e)\n",
                            pdata->dbleArray_[ii*nInputs_+jj]/pdata->dbleData_);
               if (psPlotTool_ == 1) fp = fopen("scilabaie.sci", "w");
               else                  fp = fopen("matlabaie.m", "w");
               if (fp != NULL)
               {
                  sprintf(pString,
                          " This file contains Sobol' second order indices");
                  fwriteComment(fp, pString);
                  sprintf(pString,
                          " set sortFlag = 1 and set nn to be the number");
                  fwriteComment(fp, pString);
                  sprintf(pString," of inputs to display.");
                  fwriteComment(fp, pString);
                  fprintf(fp, "sortFlag = 0;\n");
                  fprintf(fp, "nn = %d;\n", nInputs_);
                  fprintf(fp, "Mids = [\n");
                  for (ii = 0; ii < nInputs_*nInputs_; ii++) 
                     fprintf(fp,"%24.16e\n", pdata->dbleArray_[ii]);
                  fprintf(fp, "];\n");
                  if (inputNames_ == NULL)
                  {
                     fprintf(fp, "  Str = {");
                     for (ii = 0; ii < nInputs_-1; ii++) 
                        fprintf(fp,"'X%d',",ii+1);
                     fprintf(fp,"'X%d'};\n",nInputs_);
                  }
                  else
                  {
                     fprintf(fp, "  Str = {");
                     for (ii = 0; ii < nInputs_-1; ii++)
                     {
                        if (inputNames_[ii] != NULL) 
                             fprintf(fp,"'%s',",inputNames_[ii]);
                        else fprintf(fp,"'X%d',",ii+1);
                     }
                     if (inputNames_[nInputs_-1] != NULL)
                          fprintf(fp,"'%s'};\n",inputNames_[nInputs_-1]);
                     else fprintf(fp,"'X%d'};\n",nInputs_);
                  }
                  fwritePlotCLF(fp);
                  fprintf(fp, "ymin = min(Mids);\n");
                  fprintf(fp, "ymax = max(Mids);\n");
                  fprintf(fp, "h2 = 0.05 * (ymax - ymin);\n");
                  if (psPlotTool_ == 1)
                       fprintf(fp, "Mids = matrix(Mids, %d, %d);\n",
                               nInputs_,nInputs_);
                  else fprintf(fp, "Mids = reshape(Mids, %d, %d);\n",
                               nInputs_,nInputs_);
                  fprintf(fp, "Mids = Mids';\n");
                  if (psPlotTool_ == 1)
                  {
                     fprintf(fp, "drawlater\n");
                     fprintf(fp, "hist3d(Mids);\n");
                     fprintf(fp, "a=gca();\n");
                     fprintf(fp, "a.data_bounds=[0, 0, 0; %d+1, %d+1, ymax];\n",
                             nInputs_, nInputs_);
                     fprintf(fp, "newtick = a.x_ticks;\n");
                     fprintf(fp, "newtick(2) = [1:nn]';\n");
                     fprintf(fp, "newtick(3) = Str';\n");
                     fprintf(fp, "a.x_ticks = newtick;\n");
                     fprintf(fp, "a.x_label.font_size = 3;\n");
                     fprintf(fp, "a.x_label.font_style = 4;\n");
                     fprintf(fp, "a.y_ticks = newtick;\n");
                     fprintf(fp, "a.y_label.font_size = 3;\n");
                     fprintf(fp, "a.y_label.font_style = 4;\n");
                     fprintf(fp, "drawnow\n");
                  }
                  else
                  {
                     fprintf(fp, "bar3(Mids,0.8);\n");
                     fprintf(fp, "axis([0 %d+1 0 %d+1 0 ymax])\n", 
                             nInputs_, nInputs_);
                     fprintf(fp, "set(gca,'XTickLabel',Str);\n");
                     fprintf(fp, "set(gca,'YTickLabel',Str);\n");
                     fprintf(fp, "set(gca, 'fontsize', 12)\n");
                     fprintf(fp, "set(gca, 'fontweight', 'bold')\n");
                  }
                  fwritePlotAxes(fp);
                  fwritePlotTitle(fp,"Sobol 2nd Order Indices (+ 1st order)");
                  fwritePlotZLabel(fp, "Sobol Indices");
                  fwritePlotXLabel(fp, "Inputs");
                  fwritePlotYLabel(fp, "Inputs");
                  fclose(fp);
                  if (psPlotTool_ == 1)
                       printf("ie plot file = scilabaie.sci\n");
                  else printf("ie plot file = matlabaie.m\n");
               }
            }
            else
            {
               printf("Total variance = 0 ==> no interaction effect plot.\n");
            }
            pdata->clean();
         }
      }

      // +++ tsi 
      else if (!strcmp(command, "tsi"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("tsi: total sensitivity analysis (variance-based)\n");
            printf("     (suitable for raw data)\n");
            printf("syntax: tsi (after data have been loaded)\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data to analyze (load data first).\n");
            continue;
         }
         if (nInputs_ >= 20 || nSamples_ < 50*nInputs_) 
         {
            printf("This command is not recommended for small sample or\n");
            printf("large number of inputs.\n");
            printf("Use at most 20 inputs.\n");
            printf("Need at least %d sample points.\n",50*nInputs_);
            continue;
         }
         outputID = 0;
         sprintf(pString, "Enter output number (1 - %d) : ", nOutputs_);
         outputID = getInt(1, nOutputs_, pString);
         outputID--;
         TSIAnalyzer *tsiAnalyzer = new TSIAnalyzer();
         aPtr.printLevel_ = outputLevel_;
         aPtr.nSamples_ = nSamples_;
         aPtr.nInputs_ = nInputs_;
         aPtr.nOutputs_ = nOutputs_;
         aPtr.sampleInputs_ = sampleInputs_;
         aPtr.sampleOutputs_ = sampleOutputs_;
         aPtr.iLowerB_ = iLowerB_;
         aPtr.iUpperB_ = iUpperB_;
         aPtr.outputID_ = outputID;
         aPtr.sampleStates_ = sampleStates_;
         aPtr.ioPtr_ = psuadeIO_;
         tsiAnalyzer->analyze(aPtr);
         delete tsiAnalyzer;
      }

      // +++ sobol 
      else if (!strcmp(command, "sobol"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("tsi: Sobol sensitivity analysis (variance-based)\n");
            printf("     (suitable for raw data)\n");
            printf("syntax: sobol (after data have been loaded)\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data to analyze (load data first).\n");
            continue;
         }
         if (nSamples_ < 50*nInputs_) 
         {
            printf("This command is not recommended for small sample\n");
            printf("Need at least %d sample points.\n",50*nInputs_);
            continue;
         }
         outputID = 0;
         sprintf(pString, "Enter output number (1 - %d) : ", nOutputs_);
         outputID = getInt(1, nOutputs_, pString);
         outputID--;
         SobolAnalyzer *sobolAnalyzer = new SobolAnalyzer();
         aPtr.printLevel_ = outputLevel_;
         aPtr.nSamples_ = nSamples_;
         aPtr.nInputs_ = nInputs_;
         aPtr.nOutputs_ = nOutputs_;
         aPtr.sampleInputs_ = sampleInputs_;
         aPtr.sampleOutputs_ = sampleOutputs_;
         aPtr.iLowerB_ = iLowerB_;
         aPtr.iUpperB_ = iUpperB_;
         aPtr.outputID_ = outputID;
         aPtr.sampleStates_ = sampleStates_;
         aPtr.ioPtr_ = psuadeIO_;
         sobolAnalyzer->analyze(aPtr);
         delete sobolAnalyzer;
      }

      // +++ fast 
      else if (!strcmp(command, "fast"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("fast: Fourier Amplitude Sensitivity Test to\n");
            printf("      compute first order sensitivity indices.\n");
            printf("syntax: fast (after data have been loaded)\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data to analyze (load data first).\n");
            continue;
         }
         outputID = 0;
         sprintf(pString, "Enter output number (1 - %d) = ", nOutputs_);
         outputID = getInt(1, nOutputs_, pString);
         outputID--;
         analysisMethod = PSUADE_ANA_FAST;
         anaManager = new AnalysisManager();
         anaManager->setup(analysisMethod, 0);
         psuadeIO_->getParameter("ana_diagnostics",pPtr);
         ii = pPtr.intData_;
         psuadeIO_->updateAnalysisSection(-1,-1,-1,outputLevel_,-1, -1);
         anaManager->analyze(psuadeIO_, 0, NULL, outputID);
         psuadeIO_->updateAnalysisSection(-1,-1,-1,ii,-1, -1);
         delete anaManager;
      }

      // +++ meplot 
      else if (!strcmp(command, "meplot"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("meplot: main effect plot using Pgplot (obsolete)\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data to analyze (load data first).\n");
            continue;
         }
         iInd = 0;
         if (nInputs_ > 1)
         {
            sprintf(pString, "Enter input number (1 - %d) : ", nInputs_);
            iInd = getInt(1, nInputs_, pString);
            iInd--;
         }
         outputID = 0;
         sprintf(pString, "Enter output number (1 - %d) : ", nOutputs_);
         outputID = getInt(1, nOutputs_, pString);
         outputID--;
         Xmin = iLowerB_[iInd];
         Xmax = iUpperB_[iInd];
         Ymin =   1.0e35;
         Ymax = - 1.0e35;
         for (sInd = 0; sInd < nSamples_; sInd++)
         {
            if (sampleOutputs_[sInd*nOutputs_+outputID] > Ymax)
               Ymax = sampleOutputs_[sInd*nOutputs_+outputID];
            if (sampleOutputs_[sInd*nOutputs_+outputID] < Ymin)
               Ymin = sampleOutputs_[sInd*nOutputs_+outputID];
         }
         width = Xmax - Xmin;
         if (width == 0.0) Xmax = Xmin = 1.0e-3;
         Xmax = Xmax + width * 0.01;
         Xmin = Xmin - width * 0.01;
         width = Ymax - Ymin;
         if (width == 0.0) Ymax = Ymin = 1.0e-3;
         Ymax = Ymax + width * 0.01;
         Ymin = Ymin - width * 0.01;
         Plotbegin(Xmin, Xmax, Ymin, Ymax);
         nSamplesLong = (long) nSamples_;
         tempX = new double[nSamples_];
         tempY = new double[nSamples_];
         for (sInd = 0; sInd < nSamples_; sInd++)
         {
            tempX[sInd] = sampleInputs_[sInd*nInputs_+iInd];
            tempY[sInd] = sampleOutputs_[sInd*nOutputs_+outputID];
         }
         sortDbleList2(nSamples_, tempX, tempY);
         PlotScatter2D(nSamplesLong, tempX, tempY);
         printf("Enter any character and return to continue\n");
         scanf("%s", command);
         Plotend();
         fflush(stdout);
         printf("\n");
         fflush(stdout);
         fgets(lineIn,5000,stdin); 
         delete [] tempX;
         delete [] tempY;
      }

      // +++ meplot2
      else if (!strcmp(command, "meplot2"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("meplot2: 2-way interaction plot using Pgplot (obsolete)\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data to analyze (load data first).\n");
            continue;
         }
         if (nInputs_ < 2)
         {
            printf("ERROR: meplot2 requires 2 or more inputs.\n");
            continue;
         }
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("meplot2: wo-way interaction plot using Pgplot (obsolete)\n");
            continue;
         }
         iInd1 = iInd2 = -1;
         if (nInputs_ == 2)
         {
            iInd1 = 0;
            iInd2 = 1;
         }
         else
         {
            sprintf(pString, "Enter first input number (1 - %d) : ",nInputs_);
            iInd1 = getInt(1, nInputs_, pString);
            iInd1--;
            sprintf(pString, "Enter second input number (1 - %d) : ",nInputs_);
            iInd2 = getInt(1, nInputs_, pString);
            iInd2--;
         }

         outputID = 0;
         sprintf(pString, "Enter output number (1 - %d) : ", nOutputs_);
         outputID = getInt(1, nOutputs_, pString);
         outputID--;
         Xmin = iLowerB_[iInd1];
         Xmax = iUpperB_[iInd1];
         Ymin =   1.0e35;
         Ymax = - 1.0e35;
         for (sInd = 0; sInd < nSamples_; sInd++)
         {
            if (sampleOutputs_[sInd*nOutputs_+outputID] > Ymax)
               Ymax = sampleOutputs_[sInd*nOutputs_+outputID];
            if (sampleOutputs_[sInd*nOutputs_+outputID] < Ymin)
               Ymin = sampleOutputs_[sInd*nOutputs_+outputID];
         }
         width = Xmax - Xmin;
         if (width == 0.0) Xmax = Xmin = 1.0e-3;
         Xmax = Xmax + width * 0.01;
         Xmin = Xmin - width * 0.01;
         width = Ymax - Ymin;
         if (width == 0.0) Ymax = Ymin = 1.0e-3;
         Ymax = Ymax + width * 0.01;
         Ymin = Ymin - width * 0.01;
         Plotbegin(Xmin, Xmax, Ymin, Ymax);
         nSamplesLong = (long) nSamples_;
         tempInds = new double[nSamples_];
         tempX = new double[nSamples_];
         tempW = new double[nSamples_];
         tempY = new double[nSamples_];
         tempT = new double[nSamples_];
         for (sInd = 0; sInd < nSamples_; sInd++)
            tempInds[sInd] = (double) sInd;
         for (sInd = 0; sInd < nSamples_; sInd++)
         {
            tempX[sInd] = sampleInputs_[sInd*nInputs_+iInd1];
            tempW[sInd] = sampleInputs_[sInd*nInputs_+iInd2];
            tempY[sInd] = sampleOutputs_[sInd*nOutputs_+outputID];
         }
         sortDbleList2(nSamples_, tempX, tempInds);
         for (sInd = 0; sInd < nSamples_; sInd++)
         {
            ind2 = (int) tempInds[sInd];
            tempT[sInd] = tempW[ind2];
         }
         for (sInd = 0; sInd < nSamples_; sInd++)
            tempW[sInd] = tempT[sInd];
         for (sInd = 0; sInd < nSamples_; sInd++)
         {
            ind2 = (int) tempInds[sInd];
            tempT[sInd] = tempY[ind2];
         }
         for (sInd = 0; sInd < nSamples_; sInd++)
            tempY[sInd] = tempT[sInd];
         ind2 = 0;
         for (sInd = 1; sInd < nSamples_; sInd++)
         {
            if (PABS(tempX[sInd]-tempX[sInd-1]) > 1.0e-12)
            {
               sortDbleList2(sInd-ind2, &tempW[ind2], &tempY[ind2]);
               ind2 = sInd;
            }
         }
         PlotScatterM2D(nSamplesLong, tempX, tempW, tempY);
         printf("Enter any character and return to continue\n");
         scanf("%s", command);
         Plotend();
         fflush(stdout);
         printf("\n");
         fflush(stdout);
         fgets(lineIn,5000,stdin); 
         delete [] tempX;
         delete [] tempW;
         delete [] tempY;
         delete [] tempT;
         delete [] tempInds;
      }

      // +++ rsplot
      else if (!strcmp(command, "rsplot"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rsplot: two-input RS plot using Pgplot (obsolete)\n");
            printf("syntax: splot (no argument needed)\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data to analyze (load data first).\n");
            continue;
         }
         if (nInputs_ < 2)
         {
            printf("ERROR: rsplot requires 2 or more inputs.\n");
            continue;
         }
         pgPlotResponseSurface();
      }

      // +++ splot 
      else if (!strcmp(command, "splot"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("splot: create scatter plots (against each parameter).\n");
            printf("syntax: splot (no argument needed)\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data to analyze (load data first).\n");
            continue;
         }
         outputID = 0;
         sprintf(pString, "Enter output number (1 - %d) : ", nOutputs_);
         outputID = getInt(1, nOutputs_, pString);
         outputID--;

         if (psPlotTool_ == 1)
         {
            fp = fopen("scilabsp.sci", "w");
            if (fp == NULL)
            {
               printf("ERROR: cannot open file scilabsp.sci.\n");
               continue;
            }
         }
         else
         {
            fp = fopen("matlabsp.m", "w");
            if (fp == NULL)
            {
               printf("ERROR: cannot open file matlabsp.m.\n");
               continue;
            }
         }
         sprintf(pString," plotMode=0  : plot all in a single plot");
         fwriteComment(fp, pString);
         sprintf(pString," plotMode!=0 : plot one at a time");
         fwriteComment(fp, pString);
         fprintf(fp, "plotMode=0;\n");
         fprintf(fp, "Y = [\n");
         for (sInd = 0; sInd < nSamples_; sInd++)
            fprintf(fp, "%24.16e\n",sampleOutputs_[sInd*nOutputs_+outputID]);
         fprintf(fp, "];\n");
         for (iInd = 0; iInd < nInputs_; iInd++)
         {
            fprintf(fp, "X%d = [\n", iInd+1);
            for (sInd = 0; sInd < nSamples_; sInd++)
               fprintf(fp, "%24.16e\n",sampleInputs_[sInd*nInputs_+iInd]);
            fprintf(fp, "];\n");
         }
         fprintf(fp, "S = [\n");
         for (sInd = 0; sInd < nSamples_; sInd++)
            fprintf(fp, "%d\n",sampleStates_[sInd]);
         fprintf(fp, "];\n");
         if (psPlotTool_ == 1)
         {
            for (iInd = 0; iInd < nInputs_; iInd++)
            {
               fwritePlotCLF(fp);
               fprintf(fp, "drawlater\n");
               fprintf(fp, "plot(X%d,Y,'.','markersize',10)\n",iInd+1);
               fprintf(fp, "a = gca();\n");
               fprintf(fp, "a.children.children.mark_foreground = 5;\n");
               sprintf(winput, "%s vs %s", outputNames_[outputID],
                       inputNames_[iInd]);
               fwritePlotTitle(fp, winput);
               fwritePlotAxes(fp);
               fwritePlotXLabel(fp, inputNames_[iInd]);
               fwritePlotYLabel(fp, outputNames_[outputID]);
               fprintf(fp, "drawnow\n");
               if (iInd < nInputs_-1) 
               {
                  fprintf(fp, "disp(\'Press enter to advance,\')\n");
                  fprintf(fp, "halt;\n");
               }
            }
            printf("scilabsp.sci is now available for scatter plots.\n");
         }
         else
         {
            fprintf(fp, "fs=6;\n");
            fwritePlotCLF(fp);
            kk = ((int) pow(1.0*nInputs_-0.1, 0.5)) + 1;
            ll = kk;
            while ((ll - 1) * kk >= nInputs_) ll--; 
            for (iInd = 0; iInd < nInputs_; iInd++)
            {
               fprintf(fp,"if plotMode == 0\n");
               fprintf(fp,"subplot(%d,%d,%d)\n",kk,ll,iInd+1);
               fprintf(fp,"else\n");
               if (iInd > 0)
               {
                  fprintf(fp,"pause\n");
                  fprintf(fp,"disp('Press enter to continue')\n");
               }
               fwritePlotCLF(fp);
               fprintf(fp,"end;\n");
               fprintf(fp,"iset = find(S == 0);\n");
               fprintf(fp,"plot(X%d(iset),Y(iset),'rX','markersize',fs)\n",
                       iInd+1);
               fprintf(fp,"hold on\n");
               fprintf(fp,"iset = find(S == 1);\n");
               fprintf(fp,"plot(X%d(iset),Y(iset),'b*','markersize',fs)\n",
                       iInd+1);
               fprintf(fp,"hold off\n");
               fprintf(fp,"axis([%24.16e %24.16e min(Y) max(Y)])\n",
                       iLowerB_[iInd], iUpperB_[iInd]);
               sprintf(winput, "%s vs %s", outputNames_[outputID],
                       inputNames_[iInd]);
               fwritePlotTitle(fp, winput);
               fwritePlotAxes(fp);
               fwritePlotXLabel(fp, inputNames_[iInd]);
               fwritePlotYLabel(fp, outputNames_[outputID]);
            }
            printf("matlabsp.m is now available for scatter plots.\n");
         }
         fclose(fp);    
      }

      // +++ rawi2
      else if (!strcmp(command, "rawi2"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rawi2: similar to rsi2 except with sample (not RS) data\n");
            printf("syntax: rawi2 (no argument needed).\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data to analyze (load data first).\n");
            continue;
         }
         if (nInputs_ < 2)
         {
            printf("ERROR: rawi2 requires 2 inputs.\n");
            continue;
         }
         if (nOutputs_ < 2)
         {
            printf("ERROR: rawi2 requires 2 or more outputs.\n");
            continue;
         }
         if (psPlotTool_ == 1)
         {
            printf("INFO: rawi2 is currently not available for scilab.\n");
            continue;
         }
         dtemp = pow(1.0*nSamples_, 0.5) + 1.0e-12;
         nPtsPerDim = (int) dtemp;
         if (nPtsPerDim * nPtsPerDim != nSamples_)
         {
            printf("rawi2 error: nSamples must be a square.\n");
         }
         else
         { 
            rsiNOutputs = 2;
            sprintf(pString,"How many outputs to use ? (2 - %d) ",nOutputs_);
            rsiNOutputs = getInt(2, nOutputs_, pString);

            rsiSet = new int[rsiNOutputs];
            if (rsiNOutputs == nOutputs_)
            {
               for (ii = 0; ii < rsiNOutputs; ii++) rsiSet[ii] = ii;
            }
            else
            {
               for (ii = 0; ii < rsiNOutputs; ii++)
               {
                  sprintf(pString,"Enter the %d-th output index (1 - %d) : ",
                          ii+1, nOutputs_);
                  rsiSet[ii] = getInt(1, nOutputs_, pString);
                  rsiSet[ii]--;
               }
            }
            rsiMatrix = new int*[nPtsPerDim];
            for (ii = 0; ii < nPtsPerDim; ii++)
            {
               rsiMatrix[ii] = new int[nPtsPerDim];
               for (jj = 0; jj < nPtsPerDim; jj++)
                  rsiMatrix[ii][jj] = rsiNOutputs;
            }
            fp = fopen("matlabrawi2.m", "w");
            if (fp == NULL)
            {
               printf("ERROR: cannot open file matlabrawi2.m.\n");
               continue;
            }
            fwritePlotCLF(fp);
            faYOut = new double[nSamples_];

            for (ii = 0; ii < rsiNOutputs; ii++)
            {
               jplot = rsiSet[ii];
               for (sInd = 0; sInd < nSamples_; sInd++)
                  faYOut[sInd] = sampleOutputs_[sInd*nOutputs_+jplot];

               Ymin = faYOut[0];
               for (sInd = 1; sInd < nSamples_; sInd++)
                  if (faYOut[sInd] < Ymin) Ymin = faYOut[sInd];
               Ymax = faYOut[0];
               for (sInd = 1; sInd < nSamples_; sInd++)
                  if (faYOut[sInd] > Ymax) Ymax = faYOut[sInd];

               printf("Ymin and Ymax = %e %e\n", Ymin, Ymax);
               sprintf(pString,
                    "Enter the lower threshold for output %d (min = %16.8e) : ",
                    jplot, Ymin);
               threshL = getDouble(pString);
               sprintf(pString,
                    "Enter the upper threshold for output %d (max = %16.8e) : ",
                    jplot, Ymax);
               threshU = getDouble(pString);

               for (sInd = 0; sInd < nSamples_; sInd++)
               {
                  ind  = sInd % nPtsPerDim;
                  ind2 = sInd / nPtsPerDim;
                  if (faYOut[sInd] < threshL) rsiMatrix[ind][ind2]--;
                  if (faYOut[sInd] > threshU) rsiMatrix[ind][ind2]--;
               }
            }
            delete [] faYOut;

            fprintf(fp, "X = [\n");
            for (sInd = 0; sInd < nSamples_; sInd+=nPtsPerDim)
               fprintf(fp, "%e\n", sampleInputs_[sInd*2]);
            fprintf(fp, "];\n");
            fprintf(fp, "Y = [\n");
            for (sInd = 0; sInd < nPtsPerDim; sInd++)
               fprintf(fp, "%e\n", sampleInputs_[sInd*2+1]);
            fprintf(fp, "];\n");
            fprintf(fp, "A = [\n");
            count = 0;
            for (ii = 0;  ii < nPtsPerDim; ii++)
               for (jj = 0;  jj < nPtsPerDim; jj++)
                  if (rsiMatrix[jj][ii] == 0) count++;
            if (count == nPtsPerDim*nPtsPerDim)
            {
               for (ii = 0;  ii < nPtsPerDim; ii++)
                  for (jj = 0;  jj < nPtsPerDim; jj++) fprintf(fp, "0\n");
            }
            else
            {
               for (ii = 0;  ii < nPtsPerDim; ii++)
               {
                  for (jj = 0;  jj < nPtsPerDim; jj++)
                     if (rsiMatrix[jj][ii] == 0) fprintf(fp, "NaN\n");
                     else fprintf(fp, "%d\n", rsiMatrix[jj][ii]);
               }
            }
            fprintf(fp, "];\n");
            fprintf(fp, "A = reshape(A,%d,%d);\n",nPtsPerDim, nPtsPerDim);
            fprintf(fp, "A(%d,%d) = %e;\n", nPtsPerDim, nPtsPerDim, 
                    (double) rsiNOutputs);
            fprintf(fp, "contourf(X,Y,A)\n");
            fprintf(fp, "axis([%e %e %e %e])\n",iLowerB_[0],
                    iUpperB_[0],iLowerB_[1],iUpperB_[1]);
            fwritePlotAxes(fp);
            fwritePlotXLabel(fp, inputNames_[0]);
            fwritePlotYLabel(fp, inputNames_[1]);
            fwritePlotTitle(fp, "Intersection Contour");
            fprintf(fp, "colorbar\n");
            fprintf(fp, "colormap(cool)\n");
            fclose(fp);
            printf("matlabrawi2.m is now available for plotting.\n");

            delete [] rsiSet;
            for (ii = 0; ii < nPtsPerDim; ii++) delete [] rsiMatrix[ii];
            delete [] rsiMatrix;
         }
      }

      // +++ rawi3 
      else if (!strcmp(command, "rawi3"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rawi3: similar to rsi3 except with sample (not RS) data\n");
            printf("syntax: rawi3 (no argument needed).\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data to analyze (load data first).\n");
            continue;
         }
         if (nInputs_ < 3)
         {
            printf("ERROR: rawi3 requires 3 or more inputs.\n");
            continue;
         }
         if (nOutputs_ < 2)
         {
            printf("ERROR: rawi3 requires 2 or more outputs.\n");
            continue;
         }
         if (psPlotTool_ == 1)
         {
            printf("INFO: rawi3 is currently not available for scilab.\n");
            continue;
         }
         dtemp = pow(1.0*nSamples_, 0.333333) + 0.1;
         if (nPtsPerDim*nPtsPerDim*nPtsPerDim != nSamples_)
         {
            printf("rawi3 error: nSamples must be an integer 3rd power.\n");
         }
         iplot1 = 0; iplot2 = 1; iplot3 = 2;

         rsiNOutputs = 2;
         sprintf(pString,"How many outputs to use ? (2 - %d) ",nOutputs_);
         rsiNOutputs = getInt(2, nOutputs_, pString);

         rsiSet = new int[rsiNOutputs];
         if (rsiNOutputs == nOutputs_)
         {
            for (ii = 0; ii < rsiNOutputs; ii++) rsiSet[ii] = ii;
         }
         else
         {
            for (ii = 0; ii < rsiNOutputs; ii++)
            {
               sprintf(pString,"Enter the %d-th output index (1 - %d) : ",
                       ii+1, nOutputs_);
               rsiSet[ii] = getInt(1, nOutputs_, pString);
               rsiSet[ii]--;
            }
         }

         rsi3Matrix = new double**[nPtsPerDim];
         for (ii = 0; ii < nPtsPerDim; ii++)
         {
            rsi3Matrix[ii] = new double*[nPtsPerDim];
            for (jj = 0; jj < nPtsPerDim; jj++)
            {
               rsi3Matrix[ii][jj] = new double[nPtsPerDim];
               for (kk = 0; kk < nPtsPerDim; kk++)
                  rsi3Matrix[ii][jj][kk] = rsiNOutputs;
            }
         }

         faXOut = sampleInputs_;
         faYOut = new double[nSamples_];
         for (ii = 0; ii < rsiNOutputs; ii++)
         {
            jplot = rsiSet[ii];
            for (sInd = 0; sInd < nSamples_; sInd++)
               faYOut[sInd] = sampleOutputs_[sInd*nOutputs_+jplot];

            Ymax = faYOut[0];
            Ymin = faYOut[0];
            for (sInd = 1; sInd < nSamples_; sInd++)
               if (faYOut[sInd] < Ymin) Ymin = faYOut[sInd];
            for (sInd = 1; sInd < nSamples_; sInd++)
               if (faYOut[sInd] > Ymax) Ymax = faYOut[sInd];

            printf("\nOutput %d : Ymin and Ymax found = %e %e.\n", jplot,
                   Ymin, Ymax);
            sprintf(pString,"Enter the lower threshold (min = %e) : ", Ymin);
            threshL = getDouble(pString);
            sprintf(pString,"Enter the upper threshold (max = %e) : ", Ymax);
            threshU = getDouble(pString);

            for (sInd = 0; sInd < nSamples_; sInd++)
            {
               ind  = (sInd % (nPtsPerDim * nPtsPerDim)) % nPtsPerDim;
               ind2 = (sInd % (nPtsPerDim * nPtsPerDim)) / nPtsPerDim;
               kk   = sInd / (nPtsPerDim * nPtsPerDim);
               if (faYOut[sInd] < threshL) rsi3Matrix[ind][ind2][kk]--;
               if (faYOut[sInd] > threshU) rsi3Matrix[ind][ind2][kk]--;
            }
         }
         for (ii = 0; ii < nPtsPerDim; ii++)
         {
            for (jj = 0; jj < nPtsPerDim; jj++)
            {
               if (rsi3Matrix[ii][jj][kk] == 0.0)
                  rsi3Matrix[ii][jj][kk] = 0.5;
            }
         }
         delete [] faYOut;

         fp = fopen("matlabrawi3.m", "w");
         if (fp == NULL)
         {
            printf("ERROR: cannot open file matlabrawi3.m.\n");
            continue;
         }
         fwritePlotCLF(fp);
         fprintf(fp, "xlo = %e; \n", iLowerB_[1]);
         fprintf(fp, "xhi = %e; \n", iUpperB_[1]);
         fprintf(fp, "ylo = %e; \n", iLowerB_[0]);
         fprintf(fp, "yhi = %e; \n", iUpperB_[0]);
         fprintf(fp, "zlo = %e; \n", iLowerB_[2]);
         fprintf(fp, "zhi = %e; \n", iUpperB_[2]);
         fprintf(fp, "X=zeros(%d,%d,%d);\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);
         fprintf(fp, "Y=zeros(%d,%d,%d);\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);
         fprintf(fp, "Z=zeros(%d,%d,%d);\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);
         fprintf(fp, "V=zeros(%d,%d,%d);\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);
         for (jj = 0; jj < nPtsPerDim; jj++)
         {
            fprintf(fp, "Y(:,:,%d) = [\n", jj + 1);
            for (sInd = 0; sInd < nPtsPerDim; sInd++)
            {
               for (ii = 0; ii < nPtsPerDim; ii++)
               {
                  ind = sInd*nPtsPerDim*nPtsPerDim+ii*nPtsPerDim+jj;
                  fprintf(fp, "%e ", faXOut[ind*3]);
               }
               fprintf(fp, "\n");
            }
            fprintf(fp, "];\n");
            fprintf(fp, "X(:,:,%d) = [\n", jj + 1);
            for (sInd = 0; sInd < nPtsPerDim; sInd++)
            {
               for (ii = 0; ii < nPtsPerDim; ii++)
               {
                  ind = sInd*nPtsPerDim*nPtsPerDim+ii*nPtsPerDim+jj;
                  fprintf(fp, "%e ", faXOut[ind*3+1]);
               }
               fprintf(fp, "\n");
            }
            fprintf(fp, "];\n");
            fprintf(fp, "Z(:,:,%d) = [\n", jj + 1);
            for (sInd = 0; sInd < nPtsPerDim; sInd++)
            {
               for (ii = 0; ii < nPtsPerDim; ii++)
               {
                  ind = sInd*nPtsPerDim*nPtsPerDim+ii*nPtsPerDim+jj;
                  fprintf(fp, "%e ", faXOut[ind*3+2]);
               }
               fprintf(fp, "\n");
            }
            fprintf(fp, "];\n");
         }
         count = 0;
         for (jj = 0; jj < nPtsPerDim; jj++)
         {
            fprintf(fp, "V(:,:,%d) = [\n", jj + 1);
            for (sInd = 0; sInd < nPtsPerDim; sInd++)
            {
               for (ii = 0; ii < nPtsPerDim; ii++)
               {
                  ind = sInd*nPtsPerDim*nPtsPerDim+ii*nPtsPerDim+jj;
                  fprintf(fp, "%e ", rsi3Matrix[jj][ii][sInd]);
                  if (rsi3Matrix[jj][ii][sInd] == 0.5) count++;
               }
               fprintf(fp, "\n");
            }
            fprintf(fp, "];\n");
         }
         if (count == nPtsPerDim*nPtsPerDim*nPtsPerDim)
         {
            fprintf(fp, "V(1,1,1)=0;\n");
            fprintf(fp, "V(%d,%d,%d)=1;\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);
         }

         fprintf(fp, "xt = [%e:%e:%e];\n", iLowerB_[1],
                 (iUpperB_[1]-iLowerB_[1])*0.01, iUpperB_[1]);
         fprintf(fp, "yt = [%e:%e:%e];\n", iLowerB_[iplot1],
                 (iUpperB_[0]-iLowerB_[0])*0.01, iUpperB_[0]);
         fprintf(fp, "zt = [%e:%e:%e];\n", iLowerB_[2],
                 (iUpperB_[2]-iLowerB_[2])*0.01, iUpperB_[2]);
         fwritePlotCLF(fp);
         fprintf(fp, "isoval = 0.5;\n");
         fprintf(fp, "h = patch(isosurface(X,Y,Z,V,isoval),... \n");
         fprintf(fp, "          'FaceColor', 'blue', ... \n");
         fprintf(fp, "          'EdgeColor', 'none', ... \n");
         fprintf(fp, "          'AmbientStrength', 0.2, ... \n");
         fprintf(fp, "          'SpecularStrength', 0.7, ... \n");
         fprintf(fp, "          'DiffuseStrength', 0.4);\n");
         fprintf(fp, "isonormals(X,Y,Z,V,h);\n");
         fprintf(fp, "patch(isocaps(X,Y,Z,V,isoval), ...\n");
         fprintf(fp, "      'FaceColor', 'interp', ... \n");
         fprintf(fp, "      'EdgeColor', 'none'); \n");
         fprintf(fp, "axis([xlo xhi ylo yhi zlo zhi])\n");
         fprintf(fp, "daspect([%e,%e,%e])\n",iUpperB_[1]-iLowerB_[1],
                 iUpperB_[0]-iLowerB_[0], iUpperB_[2]-iLowerB_[2]);
         fprintf(fp, "colormap('default'); colorbar\n");
         fprintf(fp, "%%axis tight\n");
         fprintf(fp, "view(3) \n");
         fprintf(fp, "set(gcf,'Renderer','zbuffer')\n");
         fprintf(fp, "box on\n");
         fprintf(fp, "grid on\n");
         fprintf(fp, "lighting phong\n");
         fwritePlotAxes(fp);
         fprintf(fp, "xlabel('%s','FontSize',12,'FontWeight','bold')\n",
                 inputNames_[1]);
         fprintf(fp, "ylabel('%s','Fontsize',12,'FontWeight','bold')\n",
                 inputNames_[0]);
         fprintf(fp, "zlabel('%s','Fontsize',12,'FontWeight','bold')\n",
                 inputNames_[2]);
         fprintf(fp, "title('Intersection Contour','FontWeight',");
         fprintf(fp, "'bold','FontSize',12)\n");
         fprintf(fp, "colorbar\n");
         fprintf(fp, "colormap(cool)\n");
         fclose(fp);
         printf("matlabrawi3.m is now available for plotting.\n");
         delete [] rsiSet;
         for (ii = 0; ii < nPtsPerDim; ii++) 
         {
            for (jj = 0; jj < nPtsPerDim; jj++) 
               delete [] rsi3Matrix[ii][jj];
            delete [] rsi3Matrix[ii];
         }
         delete [] rsi3Matrix;
      }

      // +++ rspairs
      else if (!strcmp(command, "rspairs"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rspairs: generate RS of all 2-input pairs\n");
            printf("syntax: rspairs (no argument needed).\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data to analyze (load data first).\n");
            continue;
         }
         if (nInputs_ < 2)
         {
            printf("ERROR: rspairs requires 2 or more inputs.\n");
            continue;
         }
         if (psPlotTool_ == 1)
         {
            printf("INFO: rspairs is currently not available for scilab.\n");
            continue;
         }
         nPtsPerDim = 64;
         sprintf(pString, "Grid resolution ? (32 - 128) ");
         nPtsPerDim = getInt(32, 128, pString);
         faFlag = 1;
         faPtr = genFAInteractive(psuadeIO_, faFlag);
         if (faPtr == NULL) {printf("ERROR detected.\n"); continue;}
         faPtr->setNPtsPerDim(nPtsPerDim);
         faPtr->setBounds(iLowerB_, iUpperB_);
         faPtr->setOutputLevel(outputLevel_);

         nPlots = 0;
         sprintf(pString, "Enter the number of inputs to plot (1 - %d) : ",
                 nInputs_);
         nPlots = getInt(1, nInputs_, pString);
         plotIndices = new int[nInputs_];
         for (ii = 0; ii < nInputs_; ii++) plotIndices[ii] = ii;
         if (nPlots < nInputs_)
         {
            for (ii = 0; ii < nPlots; ii++)
            {
               sprintf(pString, "Enter the %d-th input (1 - %d) : ", ii+1,
                       nInputs_);
               plotIndices[ii] = getInt(1, nInputs_, pString);
               plotIndices[ii]--;
            }
         }
         inputSettings = new double[nInputs_];
         if (nInputs_ > 2)
         {
            sprintf(pString,
                    "Set other nominal values at mid point ? (y or n) ");
            getString(pString, command);
         }
         if (command[0] == 'y')
         {
            for (iInd1 = 0; iInd1 < nInputs_; iInd1++)
               inputSettings[iInd1] = 0.5*(iLowerB_[iInd1]+iUpperB_[iInd1]);
         }
         else
         {
            printf("Enter data file for nominal values. Format: \n");
            printf("PSUADE_BEGIN\n");
            printf("<numInputs>\n");
            printf("1   <data>\n");
            printf("2   <data>\n");
            printf("..  <data>\n");
            printf("PSUADE_END\n");
            printf("Data file name : ");
            scanf("%s", dataFile);
            fgets(lineIn,5000,stdin); 
            fp = fopen(dataFile, "r");
            if (fp == NULL)
            {
               printf("ERROR: data file %s not found.\n", dataFile);
               delete [] inputSettings;
               continue;
            }
            else
            {
               fscanf(fp, "%s", winput);
               if (strcmp(winput, "PSUADE_BEGIN"))
               {
                  printf("ERROR: file must begin with PSUADE_BEGIN\n");
                  fclose(fp);
                  delete [] inputSettings;
                  continue;
               }
               else
               {
                  fscanf(fp, "%d", &kk);
                  if (kk != nInputs_)
                  {
                     printf("ERROR: input size does not match nInputs.\n");
                     fclose(fp);
                     delete [] inputSettings;
                     continue;
                  }
                  for (ii = 0; ii < nInputs_; ii++)
                  {
                     fscanf(fp, "%d", &ind);
                     if (ind != (ii+1))
                     {
                        printf("ERROR: input index mismatch (%d,%d)\n",
                               ii+1,ind);
                        break;
                     }
                     fscanf(fp, "%lg", &inputSettings[ii]);
                  }
                  if (ii != nInputs_)
                  {
                     delete [] inputSettings;
                     fclose(fp);
                     continue;
                  }
                  fscanf(fp, "%s", winput);
                  fscanf(fp, "%s", winput);
                  fclose(fp);
                  if (strcmp(winput, "PSUADE_END"))
                  {
                     printf("ERROR: file must begin with PSUADE_END\n");
                     delete [] inputSettings;
                     continue;
                  }
               }
            }
         }

         fp = fopen("matlabrspairs.m", "w");
         if (fp == NULL)
         {
            printf("ERROR: cannot open file matlabrspairs.m.\n");
            continue;
         }
         fwritePlotCLF(fp);

         jplot = 0;
         sprintf(pString, "Enter output number (1 - %d) : ", nOutputs_);
         jplot = getInt(1, nOutputs_, pString);
         jplot--;
         Ymax = Ymin = sampleOutputs_[0*nOutputs_+jplot];
         for (sInd = 1; sInd < nSamples_; sInd++)
         {
            if (sampleOutputs_[sInd*nOutputs_+jplot] > Ymax) 
               Ymax = sampleOutputs_[sInd*nOutputs_+jplot];
            if (sampleOutputs_[sInd*nOutputs_+jplot] < Ymin) 
               Ymin = sampleOutputs_[sInd*nOutputs_+jplot];
         }
         printf("Ymin and Ymax found = %e %e.\n", Ymin, Ymax);
         sprintf(pString,"Set lower threshold ? (y or n) ");
         getString(pString, command);
         if (command[0] == 'y')
         {
            sprintf(pString,"Enter the lower threshold (min = %e) : ", 
                    Ymin);
            thresh = getDouble(pString);
            fprintf(fp, "Ymin = %e;\n", thresh);
         }
         else fprintf(fp, "Ymin = -1.0e35;\n");
         sprintf(pString,"Set upper threshold ? (y or n) : ");
         getString(pString, command);
         if (command[0] == 'y')
         {
            sprintf(pString,"Enter the upper threshold (max = %e) : ", 
                    Ymax);
            thresh = getDouble(pString);
            fprintf(fp, "Ymax = %e;\n", thresh);
         }
         else fprintf(fp, "Ymax = 1.0e35;\n");

         faYIn = new double[nSamples_];
         diagMax = - PSUADE_UNDEFINED;
         diagMin =   PSUADE_UNDEFINED;
         for (ii = 0; ii < nPlots; ii++)
         {
            iplot2 = plotIndices[ii];
            for (jj = 0; jj <= ii; jj++)
            {
               iplot1 = plotIndices[jj];
               if (iplot1 == iplot2)
               {
                  for (sInd = 0; sInd < nSamples_; sInd++)
                     faYIn[sInd] = sampleOutputs_[sInd*nOutputs_+jplot];
                  faPtr->gen1DGridData(sampleInputs_,faYIn,iplot2,
                             inputSettings, &faLeng, &faXOut,&faYOut);
                  fprintf(fp, "A = [\n");
                  for (sInd = 0; sInd < faLeng; sInd++)
                  {
                     fprintf(fp, "%e\n", faYOut[sInd]);
                     if (faYOut[sInd] > diagMax) diagMax = faYOut[sInd];
                     if (faYOut[sInd] < diagMin) diagMin = faYOut[sInd];
                  }
                  fprintf(fp, "];\n");
                  fprintf(fp, "X = [\n");
                  for (sInd = 0; sInd < faLeng; sInd++)
                     fprintf(fp, "%e\n", faXOut[sInd]);
                  fprintf(fp, "];\n");
                  fprintf(fp, "B = A;\n");
                  fprintf(fp, "[ia,aa] = find(A<Ymin);\n");
                  fprintf(fp, "for ii = 1 : length(ia)\n");
                  fprintf(fp, "   B(ia(ii)) = NaN;\n");
                  fprintf(fp, "end;\n");
                  fprintf(fp, "n1 = length(ia);\n");
                  fprintf(fp, "[ia,aa] = find(A>Ymax);\n");
                  fprintf(fp, "for ii = 1 : length(ia)\n");
                  fprintf(fp, "   B(ia(ii)) = NaN;\n");
                  fprintf(fp, "end;\n");
                  fprintf(fp, "n2 = length(ia);\n");
                  fprintf(fp, "if (n1 + n2 == %d)\n",faLeng);
                  fprintf(fp, "   B(1) = 0;\n");
                  fprintf(fp, "   B(%d) = 1;\n",faLeng);
                  fprintf(fp, "end;\n");
                  fprintf(fp, "subplot(%d,%d,%d), ",
                          nPlots,nPlots,ii*nPlots+jj+1);
                  fprintf(fp, "plot(X,B,'LineWidth',2.0)\n");
                  fwritePlotAxes(fp);
                  fwritePlotXLabel(fp, inputNames_[iplot1]);
                  fwritePlotYLabel(fp, outputNames_[jplot]);
                  delete [] faXOut;
                  delete [] faYOut;
                  continue;
               }

               for (sInd = 0; sInd < nSamples_; sInd++)
                  faYIn[sInd] = sampleOutputs_[sInd*nOutputs_+jplot];
               faPtr->gen2DGridData(sampleInputs_,faYIn,iplot1,iplot2,
                             inputSettings, &faLeng, &faXOut,&faYOut);

               fprintf(fp, "A = [\n");
               for (sInd = 0; sInd < faLeng; sInd++)
                  fprintf(fp, "%e\n", faYOut[sInd]);
               fprintf(fp, "];\n");
               fprintf(fp, "A = reshape(A,%d,%d);\n", nPtsPerDim, nPtsPerDim);
               fprintf(fp, "X = [\n");
               for (sInd = 0; sInd < faLeng; sInd+=nPtsPerDim)
                  fprintf(fp, "%e\n", faXOut[sInd*2]);
               fprintf(fp, "];\n");
               fprintf(fp, "Y = [\n");
               for (sInd = 0; sInd < nPtsPerDim; sInd++)
                  fprintf(fp, "%e\n", faXOut[sInd*2+1]);
               fprintf(fp, "];\n");
               fprintf(fp, "B = A;\n");
               fprintf(fp, "[ia,ja,aa] = find(A<Ymin);\n");
               fprintf(fp, "for ii = 1 : length(ia)\n");
               fprintf(fp, "   B(ia(ii),ja(ii)) = NaN;\n");
               fprintf(fp, "end;\n");
               fprintf(fp, "n1 = length(ia);\n");
               fprintf(fp, "[ia,ja,aa] = find(A>Ymax);\n");
               fprintf(fp, "for ii = 1 : length(ia)\n");
               fprintf(fp, "   B(ia(ii),ja(ii)) = NaN;\n");
               fprintf(fp, "end;\n");
               fprintf(fp, "n2 = length(ia);\n");
               fprintf(fp, "nA = size(A,1);\n");
               fprintf(fp, "if (n1 + n2 == nA*nA)\n");
               fprintf(fp, "   B(1,1) = 0;\n");
               fprintf(fp, "   B(%d,%d) = 1;\n",nPtsPerDim,nPtsPerDim);
               fprintf(fp, "end;\n");
               fprintf(fp, "subplot(%d,%d,%d), ",nPlots,nPlots,ii*nPlots+jj+1);
               fprintf(fp, "contourf(X,Y,B)\n");
               fwritePlotAxes(fp);
               fwritePlotXLabel(fp, inputNames_[iplot1]);
               fwritePlotYLabel(fp, inputNames_[iplot2]);
               fprintf(fp, "axis([%e %e %e %e])\n",iLowerB_[iplot1],
                       iUpperB_[iplot1],iLowerB_[iplot2],iUpperB_[iplot2]);
               fprintf(fp, "subplot(%d,%d,%d), ",nPlots,nPlots,jj*nPlots+ii+1);
               fprintf(fp, "surf(X,Y,B)\n");
               fwritePlotAxes(fp);
               fwritePlotXLabel(fp, inputNames_[iplot1]);
               fwritePlotYLabel(fp, inputNames_[iplot2]);
               delete [] faXOut;
               delete [] faYOut;
            }
         }
         if (diagMax - diagMin == 0) diagMax += 0.1;
         for (ii = 0; ii < nPlots; ii++)
         {
            iplot1 = plotIndices[ii];
            fprintf(fp, "subplot(%d,%d,%d), ",nPlots,nPlots,ii*nPlots+ii+1);
            fprintf(fp, "axis([%e %e %e %e])\n",iLowerB_[iplot1],
                    iUpperB_[iplot1],diagMin, diagMax);
         }
         fclose(fp);
         printf("matlabrspairs.m is now available for contour plots.\n");
         delete [] faYIn;
         delete [] inputSettings;
         delete [] plotIndices;
         plotIndices = NULL;
         delete faPtr;
         faPtr = NULL;
      }

      // +++ rsipairs
      else if (!strcmp(command, "rsipairs"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rsipairs: generate RS intersections for all input pairs\n");
            printf("syntax: rsipairs (no argument needed).\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data to analyze (load data first).\n");
            continue;
         }
         if (nInputs_ < 2)
         {
            printf("ERROR: rsipairs requires 2 or more inputs.\n");
            continue;
         }
         if (psPlotTool_ == 1)
         {
            printf("INFO: rsipairs is currently not available for scilab.\n");
            continue;
         }
         nPtsPerDim = 64;
         sprintf(pString, "Grid resolution ? (32 - 128) ");
         nPtsPerDim = getInt(32, 128, pString);
         faFlag = 1;
         faPtr = genFAInteractive(psuadeIO_, faFlag);
         if (faPtr == NULL) {printf("ERROR detected.\n"); continue;}
         faPtr->setNPtsPerDim(nPtsPerDim);
         faPtr->setBounds(iLowerB_, iUpperB_);
         faPtr->setOutputLevel(outputLevel_);

         nPlots = 0;
         sprintf(pString, "Enter the number of inputs to plot (1 - %d) : ",
                 nInputs_);
         nPlots = getInt(1, nInputs_, pString);
         plotIndices = new int[nInputs_];
         for (ii = 0; ii < nInputs_; ii++) plotIndices[ii] = ii;
         if (nPlots < nInputs_)
         {
            for (ii = 0; ii < nPlots; ii++)
            {
               sprintf(pString, "Enter the %d-th input (1 - %d) : ", ii+1,
                       nInputs_);
               plotIndices[ii] = getInt(1, nInputs_, pString);
               plotIndices[ii]--;
            }
         }
    
         printf("rsipairs: will use all outputs for constraining.\n");
         rsiNOutputs = nOutputs_;
         rsiSet = new int[rsiNOutputs];
         for (ii = 0; ii < rsiNOutputs; ii++) rsiSet[ii] = ii;
         threshLs = new double[rsiNOutputs];
         threshUs = new double[rsiNOutputs];
         faYIn = new double[nSamples_];
         for (ii = 0; ii < rsiNOutputs; ii++)
         {
            jplot = rsiSet[ii];
            Ymax = Ymin = sampleOutputs_[0*nOutputs_+jplot];
            for (sInd = 1; sInd < nSamples_; sInd++)
            {
               if (sampleOutputs_[sInd*nOutputs_+jplot] > Ymax) 
                  Ymax = sampleOutputs_[sInd*nOutputs_+jplot];
               if (sampleOutputs_[sInd*nOutputs_+jplot] < Ymin) 
                  Ymin = sampleOutputs_[sInd*nOutputs_+jplot];
            }
            printf("Output %d: Ymin and Ymax = %e %e\n",jplot+1,Ymin,Ymax);
            sprintf(pString,
                    "Enter the lower threshold (min = %16.8e) : ", Ymin);
            threshLs[ii] = getDouble(pString);
            sprintf(pString,
                    "Enter the upper threshold (max = %16.8e) : ", Ymax);
            threshUs[ii] = getDouble(pString);
         }
         for (ii = 0; ii < rsiNOutputs; ii++)
         {
            printf("Lower and upper thresholds for output %d = %e %e\n",
                   rsiSet[ii]+1, threshLs[ii], threshUs[ii]);
         }

         faYIn = new double[nSamples_];
         rsiMatrix = new int*[nPtsPerDim];
         for (ii = 0; ii < nPtsPerDim; ii++)
            rsiMatrix[ii] = new int[nPtsPerDim];
         fp = fopen("matlabrsipairs.m", "w");
         if (fp == NULL)
         {
            printf("ERROR: cannot open file matlabrsipairs.m.\n");
            continue;
         }
         fprintf(fp, "hold off\n");
         fwritePlotCLF(fp);

         inputSettings = new double[nInputs_];
         for (ii = 0; ii < nPlots; ii++)
         {
            iplot2 = plotIndices[ii];
            for (jj = 0; jj <= ii; jj++)
            {
               iplot1 = plotIndices[jj];
               for (iInd1 = 0; iInd1 < nInputs_; iInd1++)
               {
                  if (iInd1 != iplot1 && iInd1 != iplot2)
                       inputSettings[iInd1]=0.5*(iLowerB_[iInd1]+iUpperB_[iInd1]);
                  else inputSettings[iInd1]=1.0;
               }
               for (iInd1 = 0; iInd1 < nPtsPerDim; iInd1++)
                  for (iInd2 = 0; iInd2 < nPtsPerDim; iInd2++)
                     rsiMatrix[iInd1][iInd2] = rsiNOutputs;
               for (kk = 0; kk < rsiNOutputs; kk++)
               {
                  jplot = rsiSet[kk];
                  for (sInd = 0; sInd < nSamples_; sInd++)
                     faYIn[sInd] = sampleOutputs_[sInd*nOutputs_+jplot];

                  if (iplot1 == iplot2)
                  {
                     faPtr->gen1DGridData(sampleInputs_,faYIn,iplot1,
                               inputSettings, &faLeng, &faXOut,&faYOut);
                     if (kk == 0)
                     {
                        fprintf(fp, "X = [\n");
                        for (sInd = 0; sInd < faLeng; sInd++)
                           fprintf(fp, "%e\n", faXOut[sInd]);
                        fprintf(fp, "];\n");
                     }
                     for (sInd = 0; sInd < faLeng; sInd++)
                     {
                        if (faYOut[sInd]<threshLs[kk]) rsiMatrix[sInd][sInd]--;
                        if (faYOut[sInd]>threshUs[kk]) rsiMatrix[sInd][sInd]--;
                     }
                  }
                  else
                  {
                     faPtr->gen2DGridData(sampleInputs_,faYIn, iplot1, iplot2, 
                                 inputSettings, &faLeng, &faXOut,&faYOut);
                     if (kk == 0)
                     {
                        fprintf(fp, "X = [\n");
                        for (sInd = 0; sInd < faLeng; sInd+=nPtsPerDim)
                           fprintf(fp, "%e\n", faXOut[sInd*2]);
                        fprintf(fp, "];\n");
                        fprintf(fp, "Y = [\n");
                        for (sInd = 0; sInd < nPtsPerDim; sInd++)
                           fprintf(fp, "%e\n", faXOut[sInd*2+1]);
                        fprintf(fp, "];\n");
                     }
                     for (sInd = 0; sInd < faLeng; sInd++)
                     {
                        ind  = sInd % nPtsPerDim;
                        ind2 = sInd / nPtsPerDim;
                        if (faYOut[sInd] < threshLs[kk]) rsiMatrix[ind][ind2]--;
                        if (faYOut[sInd] > threshUs[kk]) rsiMatrix[ind][ind2]--;
                     }
                  }
                  delete [] faXOut;
                  delete [] faYOut;
               }
               if (iplot1 == iplot2)
               {
                  fprintf(fp, "A = [\n");
                  for (sInd = 0; sInd < nPtsPerDim; sInd++)
                     fprintf(fp, "%d\n", rsiMatrix[sInd][sInd]);
                  fprintf(fp, "];\n");
                  fprintf(fp, "subplot(%d,%d,%d), ",
                          nPlots,nPlots,ii*nPlots+jj+1);
                  fprintf(fp, "plot(X,A,'LineWidth',2.0)\n");
                  fprintf(fp, "axis([%e %e 0 %d])\n",iLowerB_[iplot1],
                          iUpperB_[iplot1],rsiNOutputs+1);
                  fwritePlotAxes(fp);
                  fwritePlotXLabel(fp, inputNames_[iplot1]);
               }
               else
               {
                  fprintf(fp, "A = [\n");
                  for (iInd1 = 0; iInd1 < nPtsPerDim; iInd1++)
                     for (iInd2 = 0; iInd2 < nPtsPerDim; iInd2++)
                        fprintf(fp, "%d\n", rsiMatrix[iInd2][iInd1]);
                  fprintf(fp, "];\n");
                  fprintf(fp, "A = reshape(A,%d,%d);\n",nPtsPerDim,
                          nPtsPerDim);
                  fprintf(fp, "A(1,1) = 0;\n");
                  fprintf(fp, "A(%d,%d) = %d;\n",nPtsPerDim,nPtsPerDim,
                          rsiNOutputs); 
                  fprintf(fp, "subplot(%d,%d,%d)\n",
                          nPlots,nPlots,ii*nPlots+jj+1);
                  fprintf(fp, "contourf(X,Y,A)\n");
                  fprintf(fp, "axis([%e %e %e %e])\n",iLowerB_[iplot1],
                          iUpperB_[iplot1],iLowerB_[iplot2],iUpperB_[iplot2]);
                  fwritePlotAxes(fp);
                  fwritePlotXLabel(fp, inputNames_[iplot1]);
                  fwritePlotYLabel(fp, inputNames_[iplot2]);
                  fprintf(fp, "subplot(%d,%d,%d)\n",
                          nPlots,nPlots,jj*nPlots+ii+1);
                  fprintf(fp, "surf(X,Y,A)\n");
                  fprintf(fp, "axis([%e %e %e %e 0 %d])\n",iLowerB_[iplot1],
                      iUpperB_[iplot1],iLowerB_[iplot2],iUpperB_[iplot2],
                      rsiNOutputs);
                  fwritePlotAxes(fp);
                  fwritePlotXLabel(fp, inputNames_[iplot1]);
                  fwritePlotYLabel(fp, inputNames_[iplot2]);
               }
            }
         }
         fclose(fp);
         printf("matlabrsipairs.m is now available for plotting.\n");

         delete [] inputSettings;
         delete [] plotIndices;
         delete [] threshLs;
         delete [] threshUs;
         delete [] rsiSet;
         delete [] faYIn;
         delete faPtr;
         for (ii = 0; ii < nPtsPerDim; ii++) delete [] rsiMatrix[ii];
         delete [] rsiMatrix;
         faPtr = NULL;
         plotIndices = NULL;
      }

      // +++ rscheck 
      else if (!strcmp(command, "rscheck"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rscheck: check the quality of the RS (training errors\n");
            printf("         and cross validation errors)\n");
            printf("syntax: rscheck (no argument needed)\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data to analyze (load data first).\n");
            continue;
         }
         faType = -1;
         sprintf(pString, "Enter your choice ? ");
         while (faType < 0 || faType > PSUADE_NUM_RS)
         {
            writeFAInfo(-1);
            faType = getFAType(pString);
         }
         outputID = 0;
         sprintf(pString, "Enter output number (1 - %d) = ", nOutputs_);
         outputID = getInt(1, nOutputs_, pString);
         outputID--;
         if (psRSExpertMode_ == 1 && psInteractive_ == 1)
         {
            printf("Available input/output transformations :\n");
            printf("0. no transformation.\n");
            printf("1. log transformation on all the inputs.\n");
            printf("2. log transformation on all the outputs.\n");
            printf("3. log transformation on all inputs and outputs.\n");
            sprintf(pString, "Enter your choice ? ");
            otrans = getInt(0, 3, pString);
            xsforms = new int[2];
            xsforms[0] = otrans & 1;
            xsforms[1] = otrans & 2;
         }
         else
         {
            xsforms = new int[2];
            xsforms[0] = 0;
            xsforms[1] = 0;
         }
         analysisMethod = PSUADE_ANA_RSFA;
         anaManager = new AnalysisManager();
         anaManager->setup(analysisMethod, faType);
         anaManager->loadLogXsformFlags(2, xsforms);
         psuadeIO_->getParameter("ana_diagnostics",pPtr);
         ii = pPtr.intData_;
         psuadeIO_->updateAnalysisSection(-1,-1,-1,outputLevel_,-1,-1);
         anaManager->analyze(psuadeIO_, 1, NULL, outputID);
         psuadeIO_->updateAnalysisSection(-1,-1,-1,ii,-1,-1);
         delete [] xsforms;
         delete anaManager;
      }

      // +++ rstest_ts 
      else if (!strcmp(command, "rstest_ts"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rstest_ts: response surface test on the training set.\n");
            printf("syntax: rstest_ts (no argument needed)\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data to analyze (load data first).\n");
            continue;
         }
         faType = -1;
         sprintf(pString, "Enter your choice ? ");
         while (faType < 0 || faType > PSUADE_NUM_RS)
         {
            writeFAInfo(-1);
            faType = getFAType(pString);
         }
         outputID = 0;
         sprintf(pString, "Enter output number (1 - %d) = ", nOutputs_);
         outputID = getInt(1, nOutputs_, pString);
         outputID--;
         faPtr = genFA(faType, nInputs_, iOne, nSamples_);
         faPtr->setBounds(iLowerB_, iUpperB_);
         faPtr->setOutputLevel(outputLevel_);
         tempY = new double[nSamples_];
         for (ss = 0; ss < nSamples_; ss++)
            tempY[ss] = sampleOutputs_[ss*nOutputs_+outputID];
         status = faPtr->initialize(sampleInputs_,tempY);
         tempV = new double[nSamples_];
         tempW = new double[nSamples_];
         faPtr->evaluatePointFuzzy(nSamples_, sampleInputs_, tempV, tempW);
         if (psPlotTool_ == 1) fp = fopen("RSTest_ts.sci", "w");
         else                  fp = fopen("RSTest_ts.m", "w");
         sprintf(pString," col 1: simulation data, col 2: rs data");
         fwriteComment(fp, pString);
         sprintf(pString," col 3: std dev");
         fwriteComment(fp, pString);
         fprintf(fp, "A = [\n");
         for (ss = 0; ss < nSamples_; ss++)
            fprintf(fp, "%e %e %e\n",tempY[ss],tempV[ss],tempW[ss]);
         fprintf(fp, "];\n");
         fwritePlotCLF(fp);
         fwritePlotFigure(fp, 1);
         fprintf(fp, "subplot(1,2,1)\n");
         if (psPlotTool_ == 1)
         {
            fprintf(fp, "histplot(10, A(:,1)-A(:,2), style=2);\n");
            fprintf(fp, "a = gce();\n");
            fprintf(fp, "a.children.fill_mode = \"on\";\n");
            fprintf(fp, "a.children.thickness = 2;\n");
            fprintf(fp, "a.children.foreground = 0;\n");
            fprintf(fp, "a.children.background = 2;\n");
         }
         else
         {
            fprintf(fp, "[nk,xk] = hist(A(:,1)-A(:,2),10);\n");
            fprintf(fp, "bar(xk,nk/%d,1.0)\n",nSamples_);
         }
         fwritePlotAxes(fp);
         fwritePlotTitle(fp, "Error Plot (unscaled)");
         fwritePlotXLabel(fp, "Error");
         fwritePlotYLabel(fp, "Probabilities");
         fprintf(fp, "subplot(1,2,2)\n");
         fprintf(fp, "xmax = max(A(:,1));\n");
         fprintf(fp, "xmin = min(A(:,1));\n");
         fprintf(fp, "ymax = max(A(:,2));\n");
         fprintf(fp, "ymin = min(A(:,2));\n");
         fprintf(fp, "xmin = min(xmin, ymin);\n");
         fprintf(fp, "xmax = max(xmax, ymax);\n");
         fprintf(fp, "XX   = xmin : xmax-xmin : xmax;\n");
         fprintf(fp, "plot(A(:,1), A(:,2),'x', XX, XX)\n");
         if (psPlotTool_ == 1)
         {
            fprintf(fp, "a = gca();\n");
            fprintf(fp, "a.data_bounds=[xmin,xmin;xmax,xmax];\n");
         }
         else fprintf(fp, "axis([xmin xmax xmin xmax])\n");
         fwritePlotAxes(fp);
         fwritePlotTitle(fp, "Interpolated versus actual data");
         fwritePlotXLabel(fp, "Actual data");
         fwritePlotYLabel(fp, "Interpolated data");
         fclose(fp);
         ddata = 0.0;
         for (ss = 0; ss < nSamples_; ss++)
            ddata += (tempY[ss] - tempV[ss]);
         ddata = ddata / nSamples_;
         printf("Training set avg error (unscaled) = %e\n", ddata);
         ddata = 0.0;
         for (ss = 0; ss < nSamples_; ss++)
         {
            ddata += (tempY[ss] - tempV[ss]);
            if (tempY[ss] != 0.0)
            {
               ddata /= PABS(tempY[ss]);
               dtemp = tempY[ss];
            }
         }
         ddata = ddata / nSamples_;
         printf("Training set avg error (  scaled) = %e (base=%e)\n",
                ddata,dtemp);
         ddata = 0.0;
         for (ss = 0; ss < nSamples_; ss++)
            ddata += pow(tempY[ss] - tempV[ss], 2.0);
         ddata = sqrt(ddata / nSamples_);
         printf("Training set rms error (unscaled) = %e\n", ddata);
         ddata = 0.0;
         for (ss = 0; ss < nSamples_; ss++)
         {
            ddata += pow(tempY[ss] - tempV[ss], 2.0);
            if (tempY[ss] != 0.0)
            {
               ddata /= pow(tempY[ss],2.0);
               dtemp = tempY[ss];
            }
         }
         ddata = sqrt(ddata / nSamples_);
         printf("Training set rms error (  scaled) = %e (base=%e)\n",
                ddata,dtemp);
         ddata = 0.0;
         for (ss = 0; ss < nSamples_; ss++)
            if (PABS(tempV[ss]-tempY[ss]) > ddata)
               ddata = PABS(tempY[ss] - tempV[ss]);
         printf("Training set max error (unscaled) = %e\n", ddata);
         ddata = 0.0;
         for (ss = 0; ss < nSamples_; ss++)
         {
            if (tempY[ss] != 0 && PABS((tempV[ss]-tempY[ss])/tempY[ss])>ddata)
            {
               dtemp = tempY[ss];
               ddata = PABS((tempY[ss]-tempV[ss])/tempY[ss]);
            }
         }
         printf("Training set max error (  scaled) = %e (base=%e)\n",
                ddata,dtemp);
         ddata = 0.0;
         for (ss = 0; ss < nSamples_; ss++)
            ddata += (tempY[ss] - tempV[ss]);
         ddata = ddata / nSamples_;
         printf("Training set error mean = %e\n", ddata);
         dtemp = 0.0;
         for (ss = 0; ss < nSamples_; ss++)
            dtemp += pow(tempY[ss] - tempV[ss] - ddata, 2.0);
         dtemp = sqrt(dtemp / (nSamples_ - 1.0));
         printf("Training set error std  = %e\n", dtemp);
         if (psPlotTool_ == 1)
              printf("rstest_ts plot file RSTest_ts.sci has been created\n");
         else printf("rstest_ts plot file RSTest_ts.m has been created\n");
         delete faPtr;
         delete [] tempY;
         delete [] tempV;
         delete [] tempW;
         faPtr = NULL;
         tempV = NULL;
         tempY = NULL;
         tempW = NULL;
      }

      // +++ rstest_cv 
      else if (!strcmp(command, "rstest_cv"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rstest_cv: cross validation response surface test.\n");
            printf("syntax: rstest_cv (no argument needed)\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data to analyze (load data first).\n");
            continue;
         }
         faType = -1;
         sprintf(pString, "Enter your choice ? ");
         while (faType < 0 || faType > PSUADE_NUM_RS)
         {
            writeFAInfo(-1);
            faType = getFAType(pString);
         }
         outputID = 0;
         sprintf(pString, "Enter output number (1 - %d) = ", nOutputs_);
         outputID = getInt(1, nOutputs_, pString);
         outputID--;
         sprintf(pString, "How many groups (2 - %d)? ", nSamples_);
         nParts = getInt(2, nSamples_, pString);
         count = nSamples_ / nParts;
         tempX = new double[nSamples_*nInputs_];
         tempY = new double[nSamples_];
         tempV = new double[nSamples_];
         for (ss = 0; ss < nSamples_; ss+=count)
         { 
            for (kk = 0; kk < ss*nInputs_; kk++) tempX[kk] = sampleInputs_[kk];
            for (kk = (ss+count)*nInputs_; kk < nSamples_*nInputs_; kk++)
               tempX[kk-count*nInputs_] = sampleInputs_[kk];
            for (kk = 0; kk < ss; kk++)
               tempY[kk] = sampleOutputs_[kk*nOutputs_+outputID];
            for (kk = ss+count; kk < nSamples_; kk++)
               tempY[kk-count] = sampleOutputs_[kk*nOutputs_+outputID];
            count2 = ss;
            if (ss+count < nSamples_) count2 += (nSamples_ - ss - count);
            faPtr = genFA(faType, nInputs_, iOne, count2);
            faPtr->setBounds(iLowerB_, iUpperB_);
            faPtr->setOutputLevel(outputLevel_);
            status = faPtr->initialize(tempX,tempY);
            count2 = nSamples_ - count2;
            faPtr->evaluatePoint(count2, &sampleInputs_[ss*nInputs_], &tempV[ss]);
            delete faPtr;
         }
         if (psPlotTool_ == 1) fp = fopen("RSTest_cv.sci", "w");
         else                  fp = fopen("RSTest_cv.m", "w");
         sprintf(pString," col 1: simulation data, col 2: rs predicted data");
         fwriteComment(fp, pString);
         fprintf(fp, "A = [\n");
         for (ss = 0; ss < nSamples_; ss++)
            fprintf(fp,"%e %e\n",sampleOutputs_[ss*nOutputs_+outputID],tempV[ss]);
         fprintf(fp, "];\n");
         fwritePlotCLF(fp);
         fwritePlotFigure(fp, 1);
         fprintf(fp, "subplot(1,2,1)\n");
         if (psPlotTool_ == 1)
         {
            fprintf(fp, "histplot(10, A(:,1)-A(:,2), style=2);\n");
            fprintf(fp, "a = gce();\n");
            fprintf(fp, "a.children.fill_mode = \"on\";\n");
            fprintf(fp, "a.children.thickness = 2;\n");
            fprintf(fp, "a.children.foreground = 0;\n");
            fprintf(fp, "a.children.background = 2;\n");
         }
         else
         {
            fprintf(fp, "[nk,xk] = hist(A(:,1)-A(:,2),10);\n");
            fprintf(fp, "bar(xk,nk/%d,1.0)\n",nSamples_);
         }
         fwritePlotAxes(fp);
         fwritePlotTitle(fp, "CV Error Plot (unscaled)");
         fwritePlotXLabel(fp, "Error");
         fwritePlotYLabel(fp, "Probabilities");
         fprintf(fp, "subplot(1,2,2)\n");
         fprintf(fp, "xmax = max(A(:,1));\n");
         fprintf(fp, "xmin = min(A(:,1));\n");
         fprintf(fp, "ymax = max(A(:,2));\n");
         fprintf(fp, "ymin = min(A(:,2));\n");
         fprintf(fp, "xmin = min(xmin, ymin);\n");
         fprintf(fp, "xmax = max(xmax, ymax);\n");
         fprintf(fp, "XX   = xmin : xmax-xmin : xmax;\n");
         fprintf(fp, "plot(A(:,1), A(:,2),'x', XX, XX)\n");
         if (psPlotTool_ == 1)
         {
            fprintf(fp, "a = gca();\n");
            fprintf(fp, "a.data_bounds=[xmin,xmin;xmax,xmax];\n");
         }
         else fprintf(fp, "axis([xmin xmax xmin xmax])\n");
         fwritePlotAxes(fp);
         fwritePlotTitle(fp, "Predicted vs actual data");
         fwritePlotXLabel(fp, "Actual data");
         fwritePlotYLabel(fp, "Predicted data");
         fclose(fp);
         for (ss = 0; ss < nSamples_; ss++)
            tempY[ss] = sampleOutputs_[ss*nOutputs_+outputID];
         ddata = 0.0;
         for (ss = 0; ss < nSamples_; ss++)
            ddata += (tempY[ss] - tempV[ss]);
         ddata = ddata / nSamples_;
         printf("CV avg error (unscaled) = %e\n", ddata);
         ddata = 0.0;
         for (ss = 0; ss < nSamples_; ss++)
         {
            ddata += (tempY[ss] - tempV[ss]);
            if (tempY[ss] != 0.0)
            {
               ddata /= PABS(tempY[ss]);
               dtemp = tempY[ss];
            }
         }
         ddata = ddata / nSamples_;
         printf("CV avg error (  scaled) = %e (base=%e)\n",ddata,dtemp);
         ddata = 0.0;
         for (ss = 0; ss < nSamples_; ss++)
            ddata += pow(tempY[ss] - tempV[ss], 2.0);
         ddata = sqrt(ddata / nSamples_);
         printf("CV rms error (unscaled) = %e\n", ddata);
         ddata = 0.0;
         for (ss = 0; ss < nSamples_; ss++)
         {
            ddata += pow(tempY[ss] - tempV[ss], 2.0);
            if (tempY[ss] != 0.0)
            {
               ddata /= pow(tempY[ss],2.0);
               dtemp = tempY[ss];
            }
         }
         ddata = sqrt(ddata / nSamples_);
         printf("CV rms error (  scaled) = %e (base=%e)\n",ddata,dtemp);
         ddata = 0.0;
         for (ss = 0; ss < nSamples_; ss++)
            if (PABS(tempV[ss]-tempY[ss]) > ddata)
               ddata = PABS(tempY[ss] - tempV[ss]);
         printf("CV max error (unscaled) = %e\n", ddata);
         ddata = 0.0;
         for (ss = 0; ss < nSamples_; ss++)
         {
            if (tempY[ss] != 0 && PABS((tempV[ss]-tempY[ss])/tempY[ss])>ddata)
            {
               dtemp = tempY[ss];
               ddata = PABS((tempY[ss]-tempV[ss])/tempY[ss]);
            }
         }
         printf("CV max error (  scaled) = %e (base=%e)\n",ddata,dtemp);
         ddata = 0.0;
         for (ss = 0; ss < nSamples_; ss++)
            ddata += (tempY[ss] - tempV[ss]);
         ddata = ddata / nSamples_;
         printf("CV error mean = %e\n", ddata);
         dtemp = 0.0;
         for (ss = 0; ss < nSamples_; ss++)
            dtemp += pow(tempY[ss] - tempV[ss] - ddata, 2.0);
         dtemp = sqrt(dtemp / (nSamples_ - 1.0));
         printf("CV error std  = %e\n", dtemp);
         if (psPlotTool_ == 1)
              printf("rstest_cv plot file RSTest_cv.sci has been created\n");
         else printf("rstest_cv plot file RSTest_cv.m has been created\n");
         delete [] tempX;
         delete [] tempY;
         delete [] tempV;
         faPtr = NULL;
         tempV = NULL;
         tempY = NULL;
         tempX = NULL;
      }

      // +++ rstest_gt 
      else if (!strcmp(command, "rstest_gt"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rstest_gt: response surface generalization test.\n");
            printf("syntax: rstest_gt (no argument needed)\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data to analyze (load data first).\n");
            continue;
         }
         printf("ERROR: no data to analyze (load data first).\n");
         faType = -1;
         sprintf(pString, "Enter your choice ? ");
         while (faType < 0 || faType > PSUADE_NUM_RS)
         {
            writeFAInfo(-1);
            faType = getFAType(pString);
         }
         outputID = 0;
         sprintf(pString, "Enter output number (1 - %d) = ", nOutputs_);
         outputID = getInt(1, nOutputs_, pString);
         outputID--;
         thresh = -1.0;
         while (thresh < 0.5 || thresh >= 1.0)
         {
            printf("This test partitions each input into two parts and\n");
            printf("uses one part to predict the other part. You can\n");
            printf("decide the percentage of the sample points for \n");
            printf("training (>= 50 percent).\n");
            sprintf(pString,"Enter this percentage in fraction [0.5 - 1) = ");
            thresh = getDouble(pString);
         }
         if (psPlotTool_ == 1) fp = fopen("RSTest_gt.sci", "w");
         else                  fp = fopen("RSTest_gt.m", "w");
         sprintf(pString," col 1-m: inputs, col m+1: simulation output");
         fwriteComment(fp, pString);
         sprintf(pString," col m+2: predicted output");
         fwriteComment(fp, pString);
         fwritePlotCLF(fp);

         int    *acnts;
         double *ameans, *avars;
         tempX = new double[nSamples_*nInputs_];
         tempV = new double[nSamples_*nInputs_];
         tempY = new double[nSamples_];
         tempW = new double[nSamples_];
         tempT = new double[nSamples_];
         ameans = new double[nInputs_*2];
         avars  = new double[nInputs_*2];
         acnts  = new int[nInputs_*2];
         printf("Generalization test results:\n");
         for (ii = 0; ii < nInputs_; ii++)
         {
            fprintf(fp, "A%d = [\n", ii+1);
            ddata = thresh * (iUpperB_[ii] - iLowerB_[ii]) + iLowerB_[ii];
            count = count2 = 0;
            for (ss = 0; ss < nSamples_; ss++)
            { 
               if (sampleInputs_[ss*nInputs_+ii] <= ddata)
               {
                  for (kk = 0; kk < nInputs_; kk++)
                     tempX[count*nInputs_+kk] = sampleInputs_[ss*nInputs_+kk];
                  tempY[count] = sampleOutputs_[ss*nOutputs_+outputID];
                  count++;
               }
               else
               {
                  for (kk = 0; kk < nInputs_; kk++)
                     tempV[count2*nInputs_+kk] = sampleInputs_[ss*nInputs_+kk];
                  tempW[count2] = sampleOutputs_[ss*nOutputs_+outputID];
                  count2++;
               }
            }
            if (count2 == 0 || count == 0)
            {
               printf("ERROR: for input %d - no training or test point.\n", 
                      ii+1);
               continue;
            }
            faPtr = genFA(faType, nInputs_, iOne, count);
            faPtr->setBounds(iLowerB_, iUpperB_);
            faPtr->setOutputLevel(outputLevel_);
            status = faPtr->initialize(tempX,tempY);
            faPtr->evaluatePoint(count2, tempV, tempT);
            delete faPtr;
            for (ss = 0; ss < count2; ss++)
            {
               for (kk = 0; kk < nInputs_; kk++)
                  fprintf(fp, "%e ", tempV[ss*nInputs_+kk]);
               fprintf(fp, "%e ", tempW[ss]);
               fprintf(fp, "%e\n", tempT[ss]);
            }
            ameans[2*ii] = 0.0;
            for (ss = 0; ss < count2; ss++) 
               ameans[2*ii] += (tempW[ss] - tempT[ss]);
            ameans[2*ii] /= (double) count2;
            avars[2*ii] = 0.0;
            for (ss = 0; ss < count2; ss++) 
               avars[2*ii] += pow(tempW[ss]-tempT[ss]-ameans[2*ii],2.0);
            avars[2*ii] /= (double) count2;
            acnts[2*ii] = count2;

            ddata = (1.0 - thresh) * (iUpperB_[ii] - iLowerB_[ii]) + iLowerB_[ii];
            count = count2 = 0;
            for (ss = 0; ss < nSamples_; ss++)
            { 
               if (sampleInputs_[ss*nInputs_+ii] > ddata)
               {
                  for (kk = 0; kk < nInputs_; kk++)
                     tempX[count*nInputs_+kk] = sampleInputs_[ss*nInputs_+kk];
                  tempY[count] = sampleOutputs_[ss*nOutputs_+outputID];
                  count++;
               }
               else
               {
                  for (kk = 0; kk < nInputs_; kk++)
                     tempV[count2*nInputs_+kk] = sampleInputs_[ss*nInputs_+kk];
                  tempW[count2] = sampleOutputs_[ss*nOutputs_+outputID];
                  count2++;
               }
            }
            if (count2 == 0 || count == 0)
            {
               printf("ERROR: for input %d - no training or test point.\n",
                      ii+1);
               continue;
            }
            faPtr = genFA(faType, nInputs_, iOne, count);
            faPtr->setBounds(iLowerB_, iUpperB_);
            faPtr->setOutputLevel(outputLevel_);
            status = faPtr->initialize(tempX,tempY);
            faPtr->evaluatePoint(count2, tempV, tempT);
            delete faPtr;
            for (ss = 0; ss < count2; ss++)
            {
               for (kk = 0; kk < nInputs_; kk++)
                  fprintf(fp, "%e ", tempV[ss*nInputs_+kk]);
               fprintf(fp, "%e ", tempW[ss]);
               fprintf(fp, "%e\n", tempT[ss]);
            }
            fprintf(fp, "];\n");
            ameans[2*ii+1] = 0.0;
            for (ss = 0; ss < count2; ss++) 
               ameans[2*ii+1] += (tempW[ss] - tempT[ss]);
            ameans[2*ii+1] /= (double) count2;
            avars[2*ii+1] = 0.0;
            for (ss = 0; ss < count2; ss++) 
               avars[2*ii+1] += pow(tempW[ss]-tempT[ss]-ameans[2*ii+1],2.0);
            avars[2*ii+1] /= (double) count2;
            acnts[2*ii+1] = count2;
            fprintf(fp, "xmax = max(A%d(:,%d+1));\n", ii+1, nInputs_);
            fprintf(fp, "xmin = min(A%d(:,%d+1));\n", ii+1, nInputs_);
            fprintf(fp, "ymax = max(A%d(:,%d+2));\n", ii+1, nInputs_);
            fprintf(fp, "ymin = min(A%d(:,%d+2));\n", ii+1, nInputs_);
            fprintf(fp, "xmin = min(xmin, ymin);\n");
            fprintf(fp, "xmax = max(xmax, ymax);\n");
            fprintf(fp, "XX   = xmin : xmax-xmin : xmax;\n");
            fprintf(fp, "plot(A%d(:,%d+1), A%d(:,%d+2),'*', XX, XX)\n",
                    ii+1, nInputs_, ii+1, nInputs_);
            if (psPlotTool_ == 1)
            {
               fprintf(fp, "a = gca();\n");
               fprintf(fp, "a.data_bounds=[xmin,xmin;xmax,xmax];\n");
            }
            else fprintf(fp, "axis([xmin xmax xmin xmax])\n");
            fwritePlotAxes(fp);
            sprintf(pString, "Extrapolated vs actual data (input = %d)",ii+1);
            fwritePlotTitle(fp, pString);
            fwritePlotXLabel(fp, "Actual data");
            fwritePlotYLabel(fp, "Predicted data");
            fprintf(fp, "disp('Press enter to continue to the next input')\n");
            if (ii < nInputs_-1) fprintf(fp, "pause\n");
         }
         for (ii = 0; ii < nInputs_; ii++)
         {
            printf("Input %4d: partition 1 error mean/var = %12.4e %12.4e (%d)\n",
                   ii+1, ameans[2*ii], avars[2*ii], acnts[2*ii]); 
            printf("Input %4d: partition 2 error mean/var = %12.4e %12.4e (%d)\n",
                   ii+1, ameans[2*ii+1], avars[2*ii+1], acnts[2*ii]+1); 
         }
         fclose(fp);
         if (psPlotTool_ == 1)
              printf("rstest_gt plot file RSTest_gt.sci has been created\n");
         else printf("rstest_gt plot file RSTest_gt.m has been created\n");
         delete [] tempX;
         delete [] tempW;
         delete [] tempT;
         delete [] tempY;
         delete [] tempV;
         delete [] acnts;
         delete [] ameans;
         delete [] avars;
         tempV = NULL;
         tempY = NULL;
         tempX = NULL;
         tempW = NULL;
         tempT = NULL;
         faPtr = NULL;
      }

      // +++ svmfind 
      else if (!strcmp(command, "svmfind"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("svmfind: find optimal parameters for SVM surface fit\n");
            printf("syntax: svmfind (no argument needed)\n");
            printf("Note: This command is for SVM kern option RBF.\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data to analyze (load data first).\n");
            continue;
         }
#ifndef HAVE_SVM
         printf("SVM not installed.\n");
#else
         outputID = 0;
         sprintf(pString, "Enter output number (1 - %d) = ", nOutputs_);
         outputID = getInt(1, nOutputs_, pString);
         outputID--;
         faType = PSUADE_RS_SVM;
         int saveMode1 = psAnaExpertMode_;
         psAnaExpertMode_ = 0;
         int saveMode2 = psRSExpertMode_;
         psRSExpertMode_ = 0;
         faPtr = genFA(faType, nInputs_, iOne, nSamples_);
         psAnaExpertMode_ = saveMode1;
         psRSExpertMode_ = saveMode2;
         faPtr->setBounds(iLowerB_, iUpperB_);
         tempY = new double[nSamples_];
         for (ii = 0; ii < nSamples_; ii++)
            tempY[ii] = sampleOutputs_[nOutputs_*ii+outputID];
         kk = 2;
         printf("Search for best gamma for RBF kernel.\n");
         printf("SVM tolerance can be adjusted for more accurate result.\n");
         printf("However, small tolerance takes more time to run.\n");
         sprintf(pString, "Enter the desired tolerance (e.g. 1e-4): ");
         tolerance = getDouble(pString);
         double gamma1=1e-10;
         sprintf(pString,"Enter lower bound of gamma to search (1e-6 - 1e-6): ");
         while (gamma1 < 1e-6 || gamma1 > 1e6)
            gamma1 = getDouble(pString);
         double gamma2=1e-10;
         sprintf(pString,"Enter upper bound of gamma to search (1e-6 - 1e-6): ");
         while (gamma2 < 1e-6 || gamma2 > 1e6 || gamma2 < gamma1)
            gamma2 = getDouble(pString);
         Ymin  = 1.0e35;
         gamma = gamma1;
         while (gamma <= gamma2)
         {
            targv[0] = (char *) &gamma;
            targv[1] = (char *) &tolerance;
            targv[2] = (char *) &kk;
            faPtr->setParams(3, targv);
            faPtr->initialize(sampleInputs_, tempY);
            sumErr = maxErr = 0.0;
            for (ii = 0; ii < nSamples_; ii++)
            {
               ddata = sampleOutputs_[ii*nOutputs_+outputID];
               dtemp = faPtr->evaluatePoint(&sampleInputs_[ii*nInputs_]);
               dtemp = PABS(dtemp - ddata);
               if (dtemp > maxErr) maxErr = dtemp;
               sumErr += (dtemp * dtemp);
            }
            sumErr = sqrt(sumErr);
            printf("svmfind: RBF_gamma = %e, error (L2n,max) = %e %e\n",
                   gamma, sumErr, maxErr);
            if (gamma == 1.0e-3 || sumErr < Ymin)
            {
               Ymin = sumErr;
               Xmin = gamma;
            }
            gamma *= 10.0;
         }
         printf("svmfind: best RBF_gamma = %e\n", Xmin);
         targv[0] = (char *) &gamma;
         gamma = Xmin;
         faPtr->setParams(1, targv);
         delete [] tempY;
         delete faPtr;
         faPtr = NULL;
#endif
      }

      // +++ rstest  or rstest_hs
      else if (!strcmp(command, "rstest") || !strcmp(command, "rstest_hs"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rstest_hs: check quality of RS with a holdout set\n");
            printf("syntax: rstest_hs (no argument needed)\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data to analyze (load data first).\n");
            continue;
         }
         strcpy(dataFile, "\0");
         sprintf(pString, "Enter test data file name : ");
         getString(pString, dataFile);
         dataFile[strlen(dataFile)-1] = '\0';
         status = 0;
         fp = fopen(dataFile, "r");
         if (fp == NULL)
         {
            printf("ERROR: Test data file %s not found.\n", dataFile);
            status = 1;
         }
         else fclose(fp);
         ioPtr = new PsuadeData();
         status = ioPtr->readPsuadeFile(dataFile);
         if (status != 0)
         {
            printf("ERROR: file %s either not found or in wrong format.\n",
                   dataFile);
            continue;
         }
         if (status == 0)
         {
            ioPtr->getParameter("output_noutputs", pPtr);
            kk = pPtr.intData_;
            if (kk > 1)
            {
               printf("ERROR: your test data file should only have 1 output\n");
               printf("       per sample point. Fix and do this again.\n");
               status = 1;
            }
            delete ioPtr;
            ioPtr = NULL;
         }
         if (status == 0)
         {
            outputID = 0;
            sprintf(pString,"Enter output number (1 - %d) = ",nOutputs_);
            outputID = getInt(1, nOutputs_, pString);
            outputID--;
            faType = -1;
            faLimit = 9;
            while (faType < 0 || faType >= PSUADE_NUM_RS)
            {
               printf("Enter a response surface method : \n");
               faLimit = writeFAInfo(-1);
               sprintf(pString, "Enter your choice ? ");
               faType = getInt(0, faLimit, pString);
            }
            if (psRSExpertMode_ == 1)
            {
               printf("Available input/output transformations :\n");
               printf("0. no transformation.\n");
               printf("1. log transformation on all the inputs.\n");
               printf("2. log transformation on all the outputs.\n");
               printf("3. log transformation on all inputs and outputs.\n");
               sprintf(pString, "Enter your choice ? ");
               otrans = getInt(0, 3, pString);
               xsforms = new int[2];
               xsforms[0] = otrans & 1;
               xsforms[1] = otrans & 2;
            }
            else
            {
               xsforms = new int[2];
               xsforms[0] = 0;
               xsforms[1] = 0;
            }
            int discFile=1, nInps, nOuts;
            printf("Name of the discrepancy model file? (enter n if none): ");
            scanf("%s", winput);
            fgets(lineIn,5000,stdin); 
            if (!strcmp(winput, "NONE") || winput[0] == 'n') discFile = 0;
            else
            {
               ioPtr = new PsuadeData();
               status = ioPtr->readPsuadeFile(winput);
               if (status == 0)
               {
                  ioPtr->getParameter("input_ninputs", pPtr);
                  nInps = pPtr.intData_;
                  if (nInps < nInputs_)
                  {
                     printf("ERROR: your sample has %d inputs but\n", nInputs_);
                     printf("       your discrepancy model has %d inputs.\n",
                            nInps);
                     delete ioPtr;
                     ioPtr = NULL;
                     delete [] xsforms;
                     continue;
                  }
                  ioPtr->getParameter("output_noutputs", pPtr);
                  nOuts = pPtr.intData_;
                  if (nOuts > 1)
                  {
                     printf("The discrepancy model has nOutputs > 1.\n");
                     printf("This is currently not supported.\n");
                     delete ioPtr;
                     ioPtr = NULL;
                     delete [] xsforms;
                     continue;
                  }
                  faPtrsRsEval = new FuncApprox*[1];
                  faPtrsRsEval[0] = genFAInteractive(ioPtr, 3);
                  delete ioPtr;
                  ioPtr = NULL;
               }
               else
               {
                  printf("ERROR: in reading the discrepancy model file %s.\n",
                         winput);
                  discFile = 0;
                  delete ioPtr;
                  ioPtr = NULL;
                  continue;
               }
            }
            double *samOuts = new double[nSamples_];
            for (ii = 0; ii < nSamples_; ii++) 
               samOuts[ii] = sampleOutputs_[ii*nOutputs_+outputID];
            if (discFile == 1)
            {
               if (psMasterMode_ == 1)
               {
                  fp = fopen("rstest_hs_data", "w");
                  if (fp != NULL)
                  {
                     fprintf(fp,"%% Col 2: RS estimate, \n");
                     fprintf(fp,"%% Col 3: model correction.\n");
                  }
               }
               else fp = NULL; 
                  
               for (ii = 0; ii < nSamples_; ii++)
               {
                  dtemp = 
                     faPtrsRsEval[0]->evaluatePoint(&sampleInputs_[ii*nInputs_]);
                  if (fp != NULL) 
                  {
                     fprintf(fp, "%6d ",ii+1);
                     for (jj = 0; jj < nInputs_; jj++)
                        fprintf(fp, "%12.4e ",sampleInputs_[ii*nInputs_+jj]);
                     fprintf(fp, " %12.4e %12.4e\n",samOuts[ii],dtemp);
                  }
                  samOuts[ii] += dtemp;
               }
               if (fp != NULL)
               {
                  printf("rstest_hs_data has simulations and discrepancies.\n");
                  fclose(fp);
               }
            }
            analysisMethod = PSUADE_ANA_RSFA;
            anaManager = new AnalysisManager();
            anaManager->setup(analysisMethod, faType);
            anaManager->loadLogXsformFlags(2, xsforms);
            aPtr.printLevel_ = outputLevel_;
            aPtr.sampleInputs_ = sampleInputs_;
            aPtr.sampleOutputs_ = samOuts;
            aPtr.nSamples_ = nSamples_;
            aPtr.nInputs_ = nInputs_;
            aPtr.nOutputs_ = 1;
            aPtr.iLowerB_ = iLowerB_;
            aPtr.iUpperB_ = iUpperB_;
            aPtr.outputID_ = 0;
            aPtr.faType_ = faType;
            aPtr.regWgtID_ = -1;
            aPtr.ioPtr_ = NULL;
            aPtr.sampleStates_ = sampleStates_;
            strcpy(pString, "validate");
            targv[0] = (char *) pString;
            targv[1] = (char *) &aPtr;
            targv[2] = (char *) dataFile;
            if (psPlotTool_ == 1) strcpy(errFile, "RSTest_hs.sci");
            else                  strcpy(errFile, "RSTest_hs.m");
            targv[3] = (char *) errFile;
            kk = 4;
            anaManager->specialRequest(PSUADE_ANA_RSFA, kk, targv);
            delete [] xsforms;
            delete anaManager;
            delete [] samOuts;
            if (discFile == 1)
            {
               delete faPtrsRsEval[0];
               delete [] faPtrsRsEval;
               faPtrsRsEval = NULL;
            }
            if (psPlotTool_ == 1)
                 printf("rstest plot file RSTest_hs.sci has been created\n");
            else printf("rstest plot file RSTest_hs.m has been created\n");
         }
      }

      // +++ rstgen 
      else if (!strcmp(command, "rstgen"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("This command creates a sample for response surface\n");
            printf("tests. Possible samples are factorial or fractional\n");
            printf("factorial to cover the corners (uniform distribution).\n");
            printf("syntax: rstgen (no argument needed)\n");
            continue;
         }
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: data not loaded yet.\n");
            printf("       Need to load the response surface sample for\n");
            printf("       extracting input information.\n");
            continue;
         }
         printEquals(PL_INFO, 0);
         printf("You have the following options:\n");
         printf("(1) If the number of inputs m is small, you can use\n");
         printf("    factorial sampling where nSamples = 2^m.\n");
         printf("(2) You can use fractional factorial sampling with\n");
         printf("    resolution V which needs nSamples=128 for m up to 11.\n");
         printf("(3) If the number of inputs m is moderate, you can use\n");
         printf("    fractional factorial sampling with resolution IV\n");
         printf("    which needs nSamples = 32 for m up to 15 and 64 for\n");
         printf("    m up to 32.\n");
         sprintf(pString, "Which option? (1 - 3) ");
         ind = getInt(1, 3, pString);
         if (ind == 1)
         {
            sampPtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_FACT);
            count = (int) pow(2.0, 1.0*nInputs_);
            if (count > 10000)
               printf("nSamples = %d may be too large.\n", nSamples_);
         }
         else if (ind == 2)
         {
            sampPtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_FF5);
            count = 128;
            if (nInputs_ > 11)
            {
               printf("nInputs = %d too large for FF5.\n", nInputs_);
               continue;
            }
         }
         else if (ind == 3)
         {
            sampPtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_FF4);
            count = 64;
            if (nInputs_ > 32)
            {
               printf("nInputs = %d too large for FF4.\n", nInputs_);
               continue;
            }
         }
         sampPtr->setPrintLevel(PL_INTERACTIVE);
         sampPtr->setInputBounds(nInputs_, iLowerB_, iUpperB_);
         sampPtr->setOutputParams(nOutputs_);
         sampPtr->setSamplingParams(count, -1, 0);
         sampPtr->initialize(0);
         count = sampPtr->getNumSamples();
         tempX  = new double[count * nInputs_];
         tempY  = new double[count * nOutputs_];
         states = new int[count];
         sampPtr->getSamples(count,nInputs_,nOutputs_,tempX,tempY,states);
         delete sampPtr;
         sampPtr = NULL;
         int    *iPDFs = new int[nInputs_];
         double *iMeans = new double[nInputs_];
         double *iStds = new double[nInputs_];
         Matrix *iCMat = new Matrix();
         iCMat->setDim(nInputs_, nInputs_);
         for (ii = 0; ii < nInputs_; ii++)
         {
            iPDFs[ii] = 0;
            iMeans[ii] = 0;
            iStds[ii] = 0;
            iCMat->setEntry(ii,ii,1.0);
         }
         ioPtr = new PsuadeData();
         ioPtr->updateInputSection(count, nInputs_, NULL, iLowerB_,
                                   iUpperB_, tempX, inputNames_, iPDFs, 
                                   iMeans, iStds, iCMat);
         delete [] iPDFs;
         delete [] iMeans;
         delete [] iStds;
         delete iCMat;
         ioPtr->updateOutputSection(count, nOutputs_, tempY, states, 
                                    outputNames_); 
         ioPtr->updateMethodSection(PSUADE_SAMP_MC, count, 1, -1, -1);
         ioPtr->writePsuadeFile("psuade_rs_test.data",0);
         printf("The test sample is in file psuade_rs_test.data.\n"); 
         delete [] tempX;
         delete [] tempY;
         delete [] states;
         delete ioPtr;
      }

      // +++ rsint 
      else if (!strcmp(command, "rsint"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rsint: compute volume under sample surface)\n");
            printf("syntax: int (no argument needed)\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL)
         {
            printf("ERROR: data not loaded yet.\n");
            continue;
         }
         outputID = 0;
         sprintf(pString, "Enter output number (1 - %d) : ", nOutputs_);
         outputID = getInt(1, nOutputs_, pString);
         outputID--;
         analysisMethod = PSUADE_ANA_INTEGRATION;
         anaManager = new AnalysisManager();
         anaManager->setup(analysisMethod, 0);
         psuadeIO_->getParameter("ana_diagnostics",pPtr);
         ii = pPtr.intData_;
         psuadeIO_->updateAnalysisSection(-1,-1,-1,outputLevel_,-1,-1);
         anaManager->analyze(psuadeIO_, 0, NULL, outputID);
         psuadeIO_->updateAnalysisSection(-1,-1,-1,ii,-1,-1);
         delete anaManager;
      }

      // +++ rsvol 
      else if (!strcmp(command, "rsvol"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rsvol: calculate volume of the constrained space\n");
            printf("syntax: rsvol (no argument needed)\n");
            printf("This command is useful for estimating how much of\n");
            printf("parameter space has been cut out due to constraints.\n");
            printf("It will report the percentage of the parameter space\n");
            printf("that is feasible in view of constraints.\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL)
         {
            printf("ERROR: data not loaded yet.\n");
            continue;
         }
         outputID = 0;
         sprintf(pString, "Enter output number (1 - %d) : ", nOutputs_);
         outputID = getInt(1, nOutputs_, pString);
         outputID--;
         Ymax = - 1.0e35;
         Ymin =   1.0e35;
         for (sInd = 0; sInd < nSamples_; sInd++)
         {
            if (sampleOutputs_[sInd*nOutputs_+outputID] > Ymax)
               Ymax = sampleOutputs_[sInd*nOutputs_+outputID];
            if (sampleOutputs_[sInd*nOutputs_+outputID] < Ymin)
               Ymin = sampleOutputs_[sInd*nOutputs_+outputID];
         }
         sprintf(pString,"Enter the lower bound constraint (Ymin=%e) : ",
                 Ymin);
         threshL = getDouble(pString);
         sprintf(pString,"Enter the upper bound constraint (Ymax=%e) : ",
                 Ymax);
         threshU = getDouble(pString);
         if (threshL >= threshU)
         {
            printf("ERROR: lower bound >= upper bound.\n");
            printf("       lower bound = %e\n", threshL);
            printf("       upper bound = %e\n", threshU);
            continue;
         }
         faFlag = 3;
         psuadeIO_->updateAnalysisSection(-1, -1, -1, -1, outputID,-1); 
         faPtr = genFAInteractive(psuadeIO_, faFlag);
         faPtr->setOutputLevel(outputLevel_);
         if (nInputs_ <= 51) 
              sampPtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_LPTAU);
         else sampPtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_MC);
         sampPtr->setPrintLevel(PL_INTERACTIVE);
         sampPtr->setInputBounds(nInputs_, iLowerB_, iUpperB_);
         sampPtr->setOutputParams(1);
         kk = 1000000;
         sampPtr->setSamplingParams(kk, 1, 0);
         sampPtr->initialize(0);
         tempX = new double[kk*nInputs_];
         tempY = new double[kk];
         states = new int[kk];
         sampPtr->getSamples(kk, nInputs_, 1, tempX, tempY, states);
         count = 0;
         for (ii = 0; ii < kk; ii++)
         {
            dtemp = faPtr->evaluatePoint(&tempX[ii*nInputs_]);
            if (dtemp >= threshL && dtemp <= threshU) count++; 
         } 
         printf("rsvol: percentage inside the constrained region = %5.2f%%.\n",
                100.0 * count / kk);
         delete faPtr;
         faPtr = NULL;
         delete sampPtr;
         sampPtr = NULL;
         delete [] tempX;
         delete [] tempY;
         delete [] states;
         tempX = NULL;
         tempY = NULL;
         states = NULL;
      }

      // +++ rs_ua, rsua, rsb_ua, rsuab and others
      else if (!strcmp(command, "rs_ua")  || !strcmp(command, "rsua") ||
               !strcmp(command, "rsb_ua") || !strcmp(command, "rs_ua2") ||
               !strcmp(command, "rs_uab") || !strcmp(command, "rsuab") ||
               !strcmp(command, "rs_uap") ||
               !strcmp(command,"rssobol1") || !strcmp(command,"rssobol2") ||
               !strcmp(command,"rssoboltsi") || !strcmp(command,"rssobolg") ||
               !strcmp(command,"rs_qsa") || !strcmp(command,"rs_qsa2") ||
               !strcmp(command, "rs1") || !strcmp(command, "rs2") ||
               !strcmp(command, "rs3") || !strcmp(command, "rs4") ||
               !strcmp(command, "rs3m") || !strcmp(command, "rssd") ||
               !strcmp(command, "rssd_ua") ||
               !strcmp(command, "rsi2") || !strcmp(command, "rsi3") ||
               !strcmp(command, "rsi3m") || !strcmp(command, "rsadd"))

      {
         newSession = new PsuadeSession();
         newSession->outputLevel_ = outputLevel_;
         newSession->nInputs_ = nInputs_;
         newSession->nOutputs_ = nOutputs_;
         newSession->nSamples_ = nSamples_;
         newSession->sampleInputs_ = sampleInputs_;
         newSession->sampleOutputs_ = sampleOutputs_;
         newSession->sampleStates_ = sampleStates_;
         newSession->inputLBounds_ = iLowerB_;
         newSession->inputUBounds_ = iUpperB_;
         newSession->inputNames_ = inputNames_;
         newSession->outputNames_ = outputNames_;
         newSession->psuadeIO_ = psuadeIO_;
         RSBasedAnalysis(lineIn, newSession);
         delete newSession;
         newSession = NULL;
      }

      // +++ rssobol1b and others
      else if (!strcmp(command,"rssobol1b") || !strcmp(command,"rssobol2b") ||
               !strcmp(command,"rssoboltsib") || !strcmp(command,"rsmeb") ||
               !strcmp(command,"rsieb"))
      {
         newSession = new PsuadeSession();
         newSession->outputLevel_ = outputLevel_;
         newSession->nInputs_ = nInputs_;
         newSession->nOutputs_ = nOutputs_;
         newSession->nSamples_ = nSamples_;
         newSession->sampleInputs_ = sampleInputs_;
         newSession->sampleOutputs_ = sampleOutputs_;
         newSession->sampleStates_ = sampleStates_;
         newSession->inputLBounds_ = iLowerB_;
         newSession->inputUBounds_ = iUpperB_;
         newSession->inputNames_ = inputNames_;
         newSession->outputNames_ = outputNames_;
         newSession->psuadeIO_ = psuadeIO_;
         RSBasedAnalysis(lineIn, newSession);
         delete newSession;
         newSession = NULL;
         if (inputNames_ != NULL)
         {
            for (ii = 0; ii < nInputs_; ii++)
               if (inputNames_[ii] != NULL) delete [] inputNames_[ii];
            delete [] inputNames_;
         }
         inputNames_ = NULL;
         if (outputNames_ != NULL)
         {
            for (ii = 0; ii < nOutputs_; ii++)
               if (outputNames_[ii] != NULL) delete [] outputNames_[ii];
            delete [] outputNames_;
         }
         outputNames_ = NULL;
         if (sampleInputs_  != NULL) delete [] sampleInputs_;
         if (sampleOutputs_ != NULL) delete [] sampleOutputs_;
         if (sampleStates_  != NULL) delete [] sampleStates_;
         if (iLowerB_       != NULL) delete [] iLowerB_;
         if (iUpperB_       != NULL) delete [] iUpperB_;
         if (inputPDFs_     != NULL) delete [] inputPDFs_;
         if (inputMeans_    != NULL) delete [] inputMeans_;
         if (inputStds_     != NULL) delete [] inputStds_;
         if (inputCMat_     != NULL) delete inputCMat_;
         sampleInputs_  = NULL;
         sampleOutputs_ = NULL;
         sampleStates_  = NULL;
         iLowerB_       = NULL;
         iUpperB_       = NULL;
         inputPDFs_     = NULL;
         inputMeans_    = NULL;
         inputStds_     = NULL;
         inputCMat_     = NULL;
         nInputs_ = nOutputs_ = nSamples_ = 0;
         printf("LOADED DATA IN LOCAL MEMORY HAVE BEEN ERASED.\n");
         printf("RE-LOAD THE ORIGINAL DATA FILE FOR MORE ANALYSIS.\n");
      }

      // +++ pca 
      else if (!strcmp(command, "pca"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("pca: principal component analysis on the outputs\n");
            printf("syntax: pca (no argument needed)\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL)
         {
            printf("ERROR: data not loaded yet.\n");
            continue;
         }
         aPtr.nSamples_ = nSamples_;
         aPtr.nInputs_  = nInputs_;
         aPtr.nOutputs_ = nOutputs_;
         aPtr.iLowerB_  = iLowerB_;
         aPtr.iUpperB_  = iUpperB_;
         aPtr.sampleInputs_ = sampleInputs_;
         aPtr.sampleOutputs_ = sampleOutputs_;
         aPtr.printLevel_ = outputLevel_;
         anaManager = new AnalysisManager();
         analysisMethod = PSUADE_ANA_PCA;
         anaManager->setup(analysisMethod, 0);
         psuadeIO_->getParameter("ana_diagnostics",pPtr);
         ii = pPtr.intData_;
         psuadeIO_->updateAnalysisSection(-1,-1,-1,outputLevel_,-1,-1);
         anaManager->analyze(psuadeIO_, 0, NULL, outputID);
         psuadeIO_->updateAnalysisSection(-1,-1,-1,ii,-1,-1);
         delete anaManager;
         psuadeIO_->getParameter("output_noutputs", pPtr);
         kk = pPtr.intData_;
         if (kk != nOutputs_) 
         {
            if (outputNames_ != NULL) 
            {
               for (ii = 0; ii < nOutputs_; ii++)
                  delete [] outputNames_[ii];
               delete [] outputNames_;
            }
            nOutputs_ = kk;
            pONames.clean();
            psuadeIO_->getParameter("output_names", pONames);
            names = pONames.strArray_;
            outputNames_ = new char*[nOutputs_+1];
            for (ii = 0; ii < nOutputs_; ii++)
            {
               outputNames_[ii] = new char[200]; 
               strcpy(outputNames_[ii], names[ii]);
            }
            if (sampleOutputs_ != NULL) delete [] sampleOutputs_;
            psuadeIO_->getParameter("output_sample", pPtr);
            sampleOutputs_   = pPtr.dbleArray_;
            pPtr.dbleArray_ = NULL;
         }
      }

      // +++ 1stest 
      else if (!strcmp(command, "1stest"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("1stest: perform one-sample tests (chi-squared, etc.)\n");
            printf("syntax: lstest (no argument needed)\n");
            continue;
         }
         analysisMethod = PSUADE_ANA_1SAMPLE;
         anaManager = new AnalysisManager();
         anaManager->setup(analysisMethod, 0);
         anaManager->analyze(analysisMethod);
         delete anaManager;
      }

      // +++ 2stest 
      else if (!strcmp(command, "2stest"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("2stest: two-sample tests (Kolmogorov Smirnov, etc.)\n");
            printf("syntax: 2stest (no argument needed)\n");
            continue;
         }
         analysisMethod = PSUADE_ANA_2SAMPLE;
         anaManager = new AnalysisManager();
         anaManager->setup(analysisMethod, 0);
         anaManager->analyze(analysisMethod);
         delete anaManager;
         fgets(lineIn,5000,stdin); 
      }

      // +++ iplot1
      else if (!strcmp(command, "iplot1"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("iplot1: plot the sample points in one parameter space\n");
            printf("syntax: iplot1 (no argument needed)\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL)
         {
            printf("ERROR: data not loaded yet.\n");
            continue;
         }
         iplot1 = -1;
         sprintf(pString, "Enter the input number (1 - %d) : ", nInputs_);
         iplot1 = getInt(1, nInputs_, pString);
         iplot1--;

         if (psPlotTool_ == 1)
         {
            fp = fopen("scilabiplt1.sci", "w");
            if (fp == NULL)
            {
               printf("ERROR: cannot open file scilabiplt1.sci.\n");
               continue;
            }
            fwritePlotCLF(fp);
            fprintf(fp, "X = [\n");
            for (sInd = 0; sInd < nSamples_; sInd++)
            {
               if (sampleOutputs_[sInd*nOutputs_] > 0.5*PSUADE_UNDEFINED ||
                   sampleStates_[sInd] != 1)
                  fprintf(fp, "%e 0\n",sampleInputs_[sInd*nInputs_+iplot1]);
               else
                  fprintf(fp, "%e 1\n",sampleInputs_[sInd*nInputs_+iplot1]);
            }
            fprintf(fp, "];\n");
            fprintf(fp, "n = %d;\n", nSamples_);
            fprintf(fp, "ia = find(X(:,2) == 0);\n");
            fprintf(fp, "if (length(ia) > 0)\n");
            fprintf(fp, "   plot(ia, X(ia,1),'rX','markerSize',13)\n");
            fprintf(fp, "   set(gca(),\"auto_clear\",\"off\")\n");
            fprintf(fp, "end;\n");
            fprintf(fp, "ia = find(X(:,2) == 1);\n");
            fprintf(fp, "if (length(ia) > 0)\n");
            fprintf(fp, "   plot(ia, X(ia,1),'b*','markerSize',13)\n");
            fprintf(fp, "end;\n");
            fprintf(fp, "set(gca(),\"auto_clear\",\"on\")\n");
            fprintf(fp, "minX = min(X(:,1)); maxX = max(X(:,1));\n");
            fprintf(fp, "a = gca();\n");
            fprintf(fp, "a.data_bounds=[0,minX;%d,maxX];\n",nSamples_);
            fwritePlotYLabel(fp, inputNames_[iplot1]);
            fwritePlotXLabel(fp, "Sample Number");
            fwritePlotAxes(fp);
            fwritePlotTitle(fp, "1D Input Data Plot");
            fclose(fp);    
            printf("scilabiplt1.sci is now ready for input scatter plots.\n");
         }
         else
         {
            fp = fopen("matlabiplt1.m", "w");
            if (fp == NULL)
            {
               printf("ERROR: cannot open file matlabiplt1.m.\n");
               continue;
            }
            fwritePlotCLF(fp);
            fprintf(fp, "X = [\n");
            for (sInd = 0; sInd < nSamples_; sInd++)
            {
               if (sampleOutputs_[sInd*nOutputs_] > 0.5*PSUADE_UNDEFINED ||
                   sampleStates_[sInd] != 1)
                  fprintf(fp, "%e 0\n",sampleInputs_[sInd*nInputs_+iplot1]);
               else
                  fprintf(fp, "%e 1\n",sampleInputs_[sInd*nInputs_+iplot1]);
            }
            fprintf(fp, "];\n");
            fprintf(fp, "iset = find(X(:,2)==0);\n");
            fprintf(fp, "plot(iset,X(iset,1),'rX','MarkerSize',13)\n");
            fprintf(fp, "hold on\n");
            fprintf(fp, "iset = find(X(:,2)==1);\n");
            fprintf(fp, "plot(iset,X(iset,1),'b*','MarkerSize',13)\n");
            fprintf(fp, "hold off\n");
            fprintf(fp, "axis([0 %d min(X(:,1)) max(X(:,1))])\n", nSamples_);
            fwritePlotYLabel(fp, inputNames_[iplot1]);
            fwritePlotXLabel(fp, "Sample Number");
            fwritePlotAxes(fp);
            fwritePlotTitle(fp, "1D Input Data Plot");
            fclose(fp);    
            printf("matlabiplt1.m is now ready for input scatter plots.\n");
         }
      }

      // +++ iplot2 
      else if (!strcmp(command, "iplot2"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("iplot2: plot sample points in two-parameter space\n");
            printf("syntax: iplot2 (no argument needed)\n");
            printf("This command is useful for examining where in the.\n");
            printf("parameter space failed runs are.\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL)
         {
            printf("ERROR: data not loaded yet.\n");
            continue;
         }
         if (nInputs_ < 2)
         {
            printf("ERROR: iplot2 requires 2 or more inputs.\n");
            continue;
         }
         iplot1 = iplot2 = -1;
         sprintf(pString, "Enter the input for x axis (1 - %d) : ", nInputs_);
         iplot1 = getInt(1, nInputs_, pString);
         iplot1--;
         if (nInputs_ == 2)
         {
            if (iplot1 == 0) iplot2 = 1;
            else             iplot2 = 0;
         }
         while (iplot2 < 0 || iplot2 >= nInputs_ || iplot1 == iplot2)
         {
            sprintf(pString, "Enter the input for y axis (1 - %d), not %d : ",
                    nInputs_, iplot1+1);
            iplot2 = getInt(1, nInputs_, pString);
            iplot2--;
            if (iplot2 == iplot1)
            {
               printf("ERROR: same index for x and y axes.\n");
               iplot2 = -1;
            }
         }

         if (psPlotTool_ == 1)
         {
            fp = fopen("scilabiplt2.sci", "w");
            if (fp == NULL)
            {
               printf("ERROR: cannot open file scilabiplt2.sci.\n");
               continue;
            }
            fwritePlotCLF(fp);
            fprintf(fp, "X = [\n");
            for (sInd = 0; sInd < nSamples_; sInd++)
            {
               if (sampleOutputs_[sInd*nOutputs_] > 0.5*PSUADE_UNDEFINED)
                  fprintf(fp, "%e %e 0\n",sampleInputs_[sInd*nInputs_+iplot1],
                          sampleInputs_[sInd*nInputs_+iplot2]);
               else
                  fprintf(fp, "%e %e 1\n",sampleInputs_[sInd*nInputs_+iplot1],
                          sampleInputs_[sInd*nInputs_+iplot2]);
            }
            fprintf(fp, "];\n");
            fprintf(fp, "ia = find(X(:,3) == 0);\n");
            fprintf(fp, "if (length(ia) > 0)\n");
            fprintf(fp, "   plot(X(ia,1),X(ia,2),'rX')\n");
            fprintf(fp, "   set(gca(),\"auto_clear\",\"off\")\n");
            fprintf(fp, "end;\n");
            fprintf(fp, "ia = find(X(:,3) == 1);\n");
            fprintf(fp, "if (length(ia) > 0)\n");
            fprintf(fp, "   plot(X(ia,1),X(ia,2),'b*')\n");
            fprintf(fp, "end;\n");
            fprintf(fp, "set(gca(),\"auto_clear\",\"on\")\n");
            fprintf(fp, "a=gca();\n");
            fprintf(fp, "a.data_bounds=[%e,%e;%e,%e];\n",
                    iLowerB_[iplot1], iLowerB_[iplot2],
                    iUpperB_[iplot1], iUpperB_[iplot2]);
            fwritePlotXLabel(fp, inputNames_[iplot1]);
            fwritePlotYLabel(fp, inputNames_[iplot2]);
            fwritePlotTitle(fp, "2D Input Data Plot");
            fwritePlotAxes(fp);
            fprintf(fp,"disp('red X: failed runs.')\n");
            fclose(fp);    
            printf("scilabiplt2.sci now has input scatter plots.\n");
         }
         else
         {
            fp = fopen("matlabiplt2.m", "w");
            if (fp == NULL)
            {
               printf("ERROR: cannot open file matlabiplt2.m.\n");
               continue;
            }
            fwritePlotCLF(fp);
            fprintf(fp, "ranflag = 0;\n");
            fprintf(fp, "X = [\n");
            for (sInd = 0; sInd < nSamples_; sInd++)
            {
               if (sampleOutputs_[sInd*nOutputs_] > 0.5*PSUADE_UNDEFINED ||
                   sampleStates_[sInd] != 1)
                  fprintf(fp, "%24.16e %24.16e 0\n",
                          sampleInputs_[sInd*nInputs_+iplot1],
                          sampleInputs_[sInd*nInputs_+iplot2]);
               else
                  fprintf(fp, "%24.16e %24.16e 1\n",
                          sampleInputs_[sInd*nInputs_+iplot1],
                          sampleInputs_[sInd*nInputs_+iplot2]);
            }
            fprintf(fp,"];\n");
            fprintf(fp,"iset = find(X(:,3)==0);\n");
            fprintf(fp,"plot(X(iset,1).*(1+ranflag*rand(size(iset,1),1)/100),");
            fprintf(fp,"X(iset,2).*(1+ranflag*rand(size(iset,1),1)/100),'rX',");
            fprintf(fp,"'markerSize',13)\n");
            fprintf(fp,"hold on\n");
            fprintf(fp,"iset = find(X(:,3)==1);\n");
            fprintf(fp,"plot(X(iset,1).*(1+ranflag*rand(size(iset,1),1)/100),");
            fprintf(fp,"X(iset,2).*(1+ranflag*rand(size(iset,1),1)/100),'b*',");
            fprintf(fp,"'markerSize',13)\n");
            fprintf(fp,"hold off\n");
            fprintf(fp,"axis([%e %e %e %e])\n", iLowerB_[iplot1], iUpperB_[iplot1],
                    iLowerB_[iplot2], iUpperB_[iplot2]);
            fwritePlotXLabel(fp, inputNames_[iplot1]);
            fwritePlotYLabel(fp, inputNames_[iplot2]);
            fwritePlotAxes(fp);
            fwritePlotTitle(fp, "2D Input Data Plot");
            fprintf(fp,"disp('red X: failed runs.')\n");
            fclose(fp);    
            printf("matlabiplt2.m is now ready for input scatter plots.\n");
         }
      }

      // +++ splot2
      else if (!strcmp(command, "splot2"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("splot2: scatter plot 2-inputs/1-output in 3D\n");
            printf("syntax: splot2 (no argument needed)\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL)
         {
            printf("ERROR: data not loaded yet.\n");
            continue;
         }
         if (nInputs_ < 2)
         {
            printf("ERROR: splot2 requires 2 or more inputs.\n");
            continue;
         }
         if (nInputs_ == 2)
         {
            iplot1 = 0;
            iplot2 = 1;
         }
         else
         {
            iplot1 = iplot2 = -1;
            sprintf(pString, "X-axis input ? (1 - %d) ",
                    nInputs_);
            iplot1 = getInt(1, nInputs_, pString);
            iplot1--;
            iplot2 = iplot1;
            while (iplot1 == iplot2)
            {
               sprintf(pString, "Y-axis input ? (1 - %d, not %d) ",
                       nInputs_, iplot1+1);
               iplot2 = getInt(1, nInputs_, pString);
               iplot2--;
               if (iplot1 == iplot2)
                  printf("ERROR: duplicate input number %d.\n",iplot2+1);
            }
         }
         if (nOutputs_ == 1) oplot1 = 0;
         else
         {
            sprintf(pString, "Z-axis output ? (1 - %d) : ",nOutputs_);
            oplot1 = getInt(1, nOutputs_, pString);
            oplot1--;
         }

         if (psPlotTool_ == 1)
         {
            fp = fopen("scilabsp2.sci", "w");
            if (fp == NULL)
            {
               printf("ERROR: cannot open file scilabsp2.sci.\n");
               continue;
            }
            fwritePlotCLF(fp);
            fprintf(fp, "X = [\n");
            for (sInd = 0; sInd < nSamples_; sInd++)
            {
               dtemp = sampleOutputs_[sInd*nOutputs_+oplot1];
               if (dtemp > 0.5*PSUADE_UNDEFINED || sampleStates_[sInd] != 1)
                  fprintf(fp, "%e %e %e 0\n",
                          sampleInputs_[sInd*nInputs_+iplot1],
                          sampleInputs_[sInd*nInputs_+iplot2],
                          sampleOutputs_[sInd*nOutputs_+oplot1]);
               else
                  fprintf(fp, "%e %e %e 1\n",
                          sampleInputs_[sInd*nInputs_+iplot1],
                          sampleInputs_[sInd*nInputs_+iplot2],
                          sampleOutputs_[sInd*nOutputs_+oplot1]);
            }
            fprintf(fp, "];\n");
            fprintf(fp, "f = gcf();\n");
            fprintf(fp, "f.color_map = jetcolormap(3);\n");
            fprintf(fp, "drawlater\n");
            fprintf(fp, "ia1 = find(X(:,4) == 0);\n");
            fprintf(fp, "if (length(ia1) > 0)\n");
            fprintf(fp, "   param3d1([X(ia,1)' ; X(ia,1)'],");
            fprintf(fp, "[X(ia,2)' ; X(ia,2)'],[X(ia,3)' ; X(ia,3)'])\n");
            fprintf(fp, "   e = gce();\n");
            fprintf(fp, "   e.children.mark_mode = \"on\";\n");
            fprintf(fp, "   e.children.mark_size_unit = \"point\";\n");
            fprintf(fp, "   e.children.mark_style = 10;\n");
            fprintf(fp, "   e.children.mark_size = 6;\n");
            fprintf(fp, "   for i = 1:length(e.children)\n");
            fprintf(fp, "      e.children(i).mark_foreground = 3;\n");
            fprintf(fp, "   end\n");
            fprintf(fp, "   set(gca(),\"auto_clear\",\"off\")\n");
            fprintf(fp, "end;\n");
            fprintf(fp, "ia2 = find(X(:,4) == 1);\n");
            fprintf(fp, "if (length(ia2) > 0)\n");
            fprintf(fp, "   param3d1([X(ia2,1)';X(ia2,1)'],[X(ia2,2)';");
            fprintf(fp, "X(ia2,2)'],[X(ia2,3)';X(ia2,3)'])\n");
            fprintf(fp, "   e = gce();\n");
            fprintf(fp, "   e.children.mark_mode = \"on\";\n");
            fprintf(fp, "   e.children.mark_size_unit = \"point\";\n");
            fprintf(fp, "   e.children.mark_style = 10;\n");
            fprintf(fp, "   e.children.mark_size = 6;\n");
            fprintf(fp, "   for i = 1:length(e.children)\n");
            fprintf(fp, "      e.children(i).mark_foreground = 1;\n");
            fprintf(fp, "   end\n");
            fprintf(fp, "end\n");
            fprintf(fp, "drawnow\n");
            fprintf(fp, "set(gca(),\"auto_clear\",\"on\")\n");
            fwritePlotXLabel(fp, inputNames_[iplot1]);
            fwritePlotYLabel(fp, inputNames_[iplot2]);
            fwritePlotZLabel(fp, outputNames_[oplot1]);
            fwritePlotTitle(fp, "3D 2-Input/1-Output Data Plot");
            fwritePlotAxes(fp);
            fprintf(fp,"disp('red  *: failed runs.')\n");
            fprintf(fp,"disp('blue X: good   runs.')\n");
            fclose(fp);
            printf("scilabsp2.sci is now ready for input scatter plots.\n");
         }
         else
         {
            fp = fopen("matlabsp2.m", "w");
            if (fp == NULL)
            {
               printf("ERROR: cannot open file matlabsp2.m.\n");
               continue;
            }
            fwritePlotCLF(fp);
            fprintf(fp, "X = [\n");
            for (sInd = 0; sInd < nSamples_; sInd++)
            {
               dtemp = sampleOutputs_[sInd*nOutputs_+oplot1];
               if (dtemp > 0.5*PSUADE_UNDEFINED || sampleStates_[sInd] != 1)
                  fprintf(fp, "%e %e %e 0\n",
                          sampleInputs_[sInd*nInputs_+iplot1],
                          sampleInputs_[sInd*nInputs_+iplot2],
                          sampleOutputs_[sInd*nOutputs_+oplot1]);
               else
                  fprintf(fp, "%e %e %e 1\n",
                          sampleInputs_[sInd*nInputs_+iplot1],
                          sampleInputs_[sInd*nInputs_+iplot2],
                          sampleOutputs_[sInd*nOutputs_+oplot1]);
            }
            fprintf(fp, "];\n");
            fprintf(fp, "iset = find(X(:,4)==0);\n");
            fprintf(fp, "if size(iset) > 0\n");
            fprintf(fp, "plot3(X(iset,1),X(iset,2),X(iset,3),'rX',");
            fprintf(fp, "'MarkerSize',13)\n");
            fprintf(fp, "hold on\n");
            fprintf(fp, "end;\n");
            fprintf(fp, "iset = find(X(:,4)==1);\n");
            fprintf(fp, "if size(iset) > 0\n");
            fprintf(fp, "plot3(X(iset,1),X(iset,2),X(iset,3),'b*',");
            fprintf(fp, "'MarkerSize',13)\n");
            fprintf(fp, "end;\n");
            fprintf(fp, "hold off\n");
            fwritePlotXLabel(fp, inputNames_[iplot1]);
            fwritePlotYLabel(fp, inputNames_[iplot2]);
            fwritePlotZLabel(fp, outputNames_[oplot1]);
            fwritePlotAxes(fp);
            fwritePlotTitle(fp, "3D 2-Input/1-Output Data Plot");
            fprintf(fp,"disp('red  *: failed runs.')\n");
            fprintf(fp,"disp('blue X: good   runs.')\n"); 
            fclose(fp);    
            printf("matlabsp2.m is now ready for input scatter plots.\n");
         }
      }

      // +++ splot3
      else if (!strcmp(command, "splot3"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("splot3: scatter plot 3-inputs/1-output in 3D\n");
            printf("syntax: splot3 (no argument needed)\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL)
         {
            printf("ERROR: data not loaded yet.\n");
            continue;
         }
         if (nInputs_ < 3)
         {
            printf("ERROR: splot2 requires 3 or more inputs.\n");
            continue;
         }
         if (nInputs_ == 3)
         {
            iplot1 = 0;
            iplot2 = 1;
            iplot3 = 2;
         }
         else
         {
            sprintf(pString, "X-axis input ? (1 - %d) ", nInputs_);
            iplot1 = getInt(1, nInputs_, pString);
            iplot1--;
            sprintf(pString, "Y-axis input ? (1 - %d) ", nInputs_);
            iplot2 = getInt(1, nInputs_, pString);
            iplot2--;
            sprintf(pString, "Z-axis input ? (1 - %d) ", nInputs_);
            iplot3 = getInt(1, nInputs_, pString);
            iplot3--;
         }
         if (nOutputs_ == 1) oplot1 = 0;
         else
         {
            sprintf(pString, "Which output ? (1 - %d) : ",nOutputs_);
            oplot1 = getInt(1, nOutputs_, pString);
            oplot1--;
         }

         if (psPlotTool_ == 1)
         {
            printf("ERROR: this command is currently not supported on scilab.\n");
         }
         else
         {
            fp = fopen("matlabsp3.m", "w");
            if (fp == NULL)
            {
               printf("ERROR: cannot open file matlabsp3.m.\n");
               continue;
            }
            fwritePlotCLF(fp);
            fprintf(fp, "scale = 1;\n");
            fprintf(fp, "X = [\n");
            for (sInd = 0; sInd < nSamples_; sInd++)
            {
               dtemp = sampleOutputs_[sInd*nOutputs_+oplot1];
               if (dtemp > 0.5*PSUADE_UNDEFINED || sampleStates_[sInd] != 1)
                  fprintf(fp, "%e %e %e %e 0\n",
                          sampleInputs_[sInd*nInputs_+iplot1],
                          sampleInputs_[sInd*nInputs_+iplot2],
                          sampleInputs_[sInd*nInputs_+iplot3],
                          sampleOutputs_[sInd*nOutputs_+oplot1]);
               else
                  fprintf(fp, "%e %e %e %e 1\n",
                          sampleInputs_[sInd*nInputs_+iplot1],
                          sampleInputs_[sInd*nInputs_+iplot2],
                          sampleInputs_[sInd*nInputs_+iplot3],
                          sampleOutputs_[sInd*nOutputs_+oplot1]);
            }
            fprintf(fp, "];\n");
            fprintf(fp, "ymin = min(X(:,4));\n");
            fprintf(fp, "ywid = max(X(:,4)) - ymin;\n");
            fprintf(fp, "if (ywid == 0)\n");
            fprintf(fp, "  ywid = 1;\n");
            fprintf(fp, "end;\n");
            fprintf(fp, "ywid = 1 / ywid;\n");
            fprintf(fp, "for ii = 1 : %d\n", nSamples_);
            fprintf(fp, "  if X(ii,5) == 0\n");
            fprintf(fp, "    dsize = (X(ii,4) - ymin) * ywid;\n");
            fprintf(fp, "    isize = ceil((dsize * 100 + 10)*scale);\n");
            fprintf(fp, "    plot3(X(ii,1),X(ii,2),X(ii,3),'r.',");
            fprintf(fp, "'MarkerSize',isize)\n");
            fprintf(fp, "  else\n");
            fprintf(fp, "    dsize = (X(ii,4) - ymin) * ywid;\n");
            fprintf(fp, "    isize = ceil((dsize * 100 + 10)*scale);\n");
            fprintf(fp, "    plot3(X(ii,1),X(ii,2),X(ii,3),'b.',");
            fprintf(fp, "'MarkerSize',isize)\n");
            fprintf(fp, "  end;\n");
            fprintf(fp, "  if ii == 1\n");
            fprintf(fp, "    hold on\n");
            fprintf(fp, "  end;\n");
            fprintf(fp, "end;\n");
            fprintf(fp, "hold off\n");
            fwritePlotXLabel(fp, inputNames_[iplot1]);
            fwritePlotYLabel(fp, inputNames_[iplot2]);
            fwritePlotZLabel(fp, outputNames_[iplot3]);
            fwritePlotAxes(fp);
            fwritePlotTitle(fp, "3D 3-Input/1-Output Data Plot");
            fprintf(fp,"disp('red  *: failed runs.')\n");
            fprintf(fp,"disp('blue X: good   runs.')\n"); 
            fprintf(fp,"disp('Note: sizes of the dots are output magnitudes.')\n");
            fprintf(fp,"disp('Note: use scale to adjust relative sizes of dots')\n");
            fclose(fp);    
            printf("matlabsp3.m is now ready for input scatter plots.\n");
         }
      }

      // +++ iplot3
      else if (!strcmp(command, "iplot3"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("iplot3: plot the sample points in three-parameter space\n");
            printf("syntax: iplot3 (no argument needed)\n");
            printf("This command is useful for examining where in the.\n");
            printf("parameter space failed runs are.\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL)
         {
            printf("ERROR: data not loaded yet.\n");
            continue;
         }
         if (nInputs_ < 3)
         {
            printf("ERROR: iplot3 requires 3 or more inputs.\n");
            continue;
         }
         iplot1 = iplot2 = iplot3 = -1;
         sprintf(pString, "Enter the input for x axis (1 - %d) : ", nInputs_);
         iplot1 = getInt(1, nInputs_, pString);
         iplot1--;
         iplot2 = iplot1;
         while (iplot1 == iplot2)
         {
            sprintf(pString,"Enter the input for y axis (1 - %d), not %d : ",
                    nInputs_, iplot1+1);
            iplot2 = getInt(1, nInputs_, pString);
            iplot2--;
            if (iplot1 == iplot2)
               printf("ERROR: duplicate input number %d.\n",iplot2+1);
         }
         if (nInputs_ == 3) iplot3 = 3 - iplot1 - iplot2;
         while (iplot3 < 0 || iplot3 == iplot1 || iplot3 == iplot2)
         {
            sprintf(pString,
                    "Enter the input for z axis (1 - %d), not %d nor %d: ",
                    nInputs_, iplot1+1, iplot2+1);
            iplot3 = getInt(1, nInputs_, pString);
            iplot3--;
            if (iplot3 == iplot1 || iplot3 == iplot2)
               printf("ERROR: duplicate input number %d.\n",iplot3+1);
         }

         if (psPlotTool_ == 1)
         {
            fp = fopen("scilabiplt3.sci", "w");
            if (fp == NULL)
            {
               printf("ERROR: cannot open file scilabiplt3.sci.\n");
               continue;
            }
            fwritePlotCLF(fp);
            fprintf(fp, "X = [\n");
            for (sInd = 0; sInd < nSamples_; sInd++)
            {
               dtemp = sampleOutputs_[sInd*nOutputs_];
               if (dtemp > 0.5*PSUADE_UNDEFINED || sampleStates_[sInd] != 1)
                  fprintf(fp, "%e %e %e 0\n",
                          sampleInputs_[sInd*nInputs_+iplot1],
                          sampleInputs_[sInd*nInputs_+iplot2],
                          sampleInputs_[sInd*nInputs_+iplot3]);
               else
                  fprintf(fp, "%e %e %e 1\n",
                          sampleInputs_[sInd*nInputs_+iplot1],
                          sampleInputs_[sInd*nInputs_+iplot2],
                          sampleInputs_[sInd*nInputs_+iplot3]);
            }
            fprintf(fp, "];\n");
            fprintf(fp, "f = gcf();\n");
            fprintf(fp, "f.color_map = jetcolormap(3);\n");
            fprintf(fp, "drawlater\n");
            fprintf(fp, "ia1 = find(X(:,4) == 0);\n");
            fprintf(fp, "if (length(ia1) > 0)\n");
            fprintf(fp, "   param3d1([X(ia,1)' ; X(ia,1)'],");
            fprintf(fp, "[X(ia,2)' ; X(ia,2)'],[X(ia,3)' ; X(ia,3)'])\n");
            fprintf(fp, "   e = gce();\n");
            fprintf(fp, "   e.children.mark_mode = \"on\";\n");
            fprintf(fp, "   e.children.mark_size_unit = \"point\";\n");
            fprintf(fp, "   e.children.mark_style = 10;\n");
            fprintf(fp, "   e.children.mark_size = 6;\n");
            fprintf(fp, "   for i = 1:length(e.children)\n");
            fprintf(fp, "      e.children(i).mark_foreground = 3;\n");
            fprintf(fp, "   end\n");
            fprintf(fp, "   set(gca(),\"auto_clear\",\"off\")\n");
            fprintf(fp, "end;\n");
            fprintf(fp, "ia2 = find(X(:,4) == 1);\n");
            fprintf(fp, "if (length(ia2) > 0)\n");
            fprintf(fp, "   param3d1([X(ia2,1)';X(ia2,1)'],[X(ia2,2)';");
            fprintf(fp, "X(ia2,2)'],[X(ia2,3)';X(ia2,3)'])\n");
            fprintf(fp, "   e = gce();\n");
            fprintf(fp, "   e.children.mark_mode = \"on\";\n");
            fprintf(fp, "   e.children.mark_size_unit = \"point\";\n");
            fprintf(fp, "   e.children.mark_style = 10;\n");
            fprintf(fp, "   e.children.mark_size = 6;\n");
            fprintf(fp, "   for i = 1:length(e.children)\n");
            fprintf(fp, "      e.children(i).mark_foreground = 1;\n");
            fprintf(fp, "   end\n");
            fprintf(fp, "end\n");
            fprintf(fp, "drawnow\n");
            fprintf(fp, "set(gca(),\"auto_clear\",\"on\")\n");
            fwritePlotXLabel(fp, inputNames_[iplot1]);
            fwritePlotYLabel(fp, inputNames_[iplot2]);
            fwritePlotZLabel(fp, inputNames_[iplot3]);
            fwritePlotTitle(fp, "3D Input Data Plot");
            fwritePlotAxes(fp);
            fclose(fp);    
            printf("scilabiplt3.sci now has input scatter plots.\n");
         }
         else
         {
            fp = fopen("matlabiplt3.m", "w");
            if (fp == NULL)
            {
               printf("ERROR: cannot open file matlabiplt3.m.\n");
               continue;
            }
            sprintf(pString," ranflag: to distinguish identical points");
            fwriteComment(fp, pString);
            sprintf(pString,"     by adding a small perturbation (when on)");
            fwriteComment(fp, pString);
            fprintf(fp, "ranflag  = 0;\n");
            sprintf(pString," set cvFlag to 1 to use convex hull");
            fwriteComment(fp, pString);
            fprintf(fp, "cvFlag = 0;\n");
            fwritePlotCLF(fp);
            fprintf(fp, "X = [\n");
            for (sInd = 0; sInd < nSamples_; sInd++)
            {
               dtemp = sampleOutputs_[sInd*nOutputs_];
               if (dtemp > 0.5*PSUADE_UNDEFINED || sampleStates_[sInd] != 1)
                  fprintf(fp, "%e %e %e 0\n",
                          sampleInputs_[sInd*nInputs_+iplot1],
                          sampleInputs_[sInd*nInputs_+iplot2],
                          sampleInputs_[sInd*nInputs_+iplot3]);
               else
                  fprintf(fp, "%e %e %e 1\n",
                          sampleInputs_[sInd*nInputs_+iplot1],
                          sampleInputs_[sInd*nInputs_+iplot2],
                          sampleInputs_[sInd*nInputs_+iplot3]);
            }
            fprintf(fp, "];\n");
            fprintf(fp, "iset = find(X(:,4)==1);\n");
            fprintf(fp, "plot3(X(iset,1).*(1+ranflag*rand(size(iset,1),1)/100),");
            fprintf(fp, "X(iset,2).*(1+ranflag*rand(size(iset,1),1)/100),");
            fprintf(fp, "X(iset,3).*(1+ranflag*rand(size(iset,1),1)/100),'b*',");
            fprintf(fp, "'MarkerSize',13)\n");
            fprintf(fp, "hold on\n");
            fprintf(fp, "iset = find(X(:,4)==0);\n");
            fprintf(fp, "plot3(X(iset,1).*(1+ranflag*rand(size(iset,1),1)/100),");
            fprintf(fp, "X(iset,2).*(1+ranflag*rand(size(iset,1),1)/100),");
            fprintf(fp, "X(iset,3).*(1+ranflag*rand(size(iset,1),1)/100),'rX',");
            fprintf(fp, "'MarkerSize',13)\n");
            fprintf(fp, "hold off\n");
            fwritePlotXLabel(fp, inputNames_[iplot1]);
            fwritePlotYLabel(fp, inputNames_[iplot2]);
            fwritePlotZLabel(fp, inputNames_[iplot3]);
            fwritePlotAxes(fp);
            fwritePlotTitle(fp, "3D Input Scatter Plot");
            fprintf(fp, "if cvFlag == 1\n");
            fprintf(fp, "   figure\n");
            fprintf(fp, "   iset = find(X(:,4)==0);\n");
            fprintf(fp, "   m = size(iset,1);\n");
            fprintf(fp, "   if (m > 3)\n");
            fprintf(fp, "      XX = X(iset,1:3);\n");
            fprintf(fp, "      KK = convhulln(XX);\n");
            fprintf(fp, "      trisurf(KK,XX(:,1),XX(:,2),XX(:,3))\n");
            fwritePlotXLabel(fp, inputNames_[iplot1]);
            fwritePlotYLabel(fp, inputNames_[iplot2]);
            fwritePlotZLabel(fp, inputNames_[iplot3]);
            fwritePlotAxes(fp);
            fwritePlotTitle(fp,"3D Input Convex Hull for Good Points");
            fprintf(fp, "   else\n");
            fprintf(fp, "      disp('too few points to display')\n");
            fprintf(fp, "   end;\n");
            fprintf(fp, "end;\n");
            fclose(fp);    
            printf("matlabiplt3.m is now ready for input scatter plots.\n");
            printf("Note: see inside the matlab file to see options to \n");
            printf("      distinguish between valid and invalid runs.\n");
         }
      }

      // +++ iplot4m
      else if (!strcmp(command, "iplot4m"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("iplot4: plot the sample points in four-parameter space\n");
            printf("syntax: iplot4 (no argument needed)\n");
            printf("This command is useful for examining where in the\n");
            printf("parameter space failed runs are. Fourth dimension is\n");
            printf("in time.\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL)
         {
            printf("ERROR: data not loaded yet.\n");
            continue;
         }
         if (nInputs_ < 4)
         {
            printf("ERROR: iplot4m requires 4 or more inputs.\n");
            continue;
         }
         if (psPlotTool_ == 1)
         {
            printf("INFO: iplot4m is currently not available for scilab.\n");
            continue;
         }
         iplot1 = iplot2 = iplot3 = iplot4 = -1;
         sprintf(pString, "Enter the input for x axis (1 - %d) : ", nInputs_);
         iplot1 = getInt(1, nInputs_, pString);
         iplot1--;
         iplot2 = iplot1;
         while (iplot1 == iplot2)
         {
            sprintf(pString,"Enter the input for y axis (1 - %d), not %d : ",
                    nInputs_, iplot1+1);
            iplot2 = getInt(1, nInputs_, pString);
            iplot2--;
            if (iplot1 == iplot2)
               printf("ERROR: duplicate input number %d.\n",iplot2+1);
         }
         iplot3 = -1;
         while (iplot3 < 0 || iplot3 == iplot1 || iplot3 == iplot2)
         {
            sprintf(pString,
                    "Enter the input for z axis (1 - %d), not %d nor %d: ",
                    nInputs_, iplot1+1, iplot2+1);
            iplot3 = getInt(1, nInputs_, pString);
            iplot3--;
            if (iplot3 == iplot1 || iplot3 == iplot2)
               printf("ERROR: duplicate input number %d.\n",iplot3+1);
         }
         iplot4 = -1;
         while (iplot4 < 0 || iplot4 == iplot1 || iplot4 == iplot2 || 
                iplot4 == iplot3)
         {
            sprintf(pString,
               "Enter the input for t axis (1 - %d), not %d nor %d nor %d: ",
                    nInputs_, iplot1+1, iplot2+1, iplot3+1);
            iplot4 = getInt(1, nInputs_, pString);
            iplot4--;
            if (iplot4 == iplot1 || iplot4 == iplot2 || iplot4 == iplot3)
               printf("ERROR: duplicate input number %d.\n",iplot4+1);
         }

         fp = fopen("matlabiplt4m.m", "w");
         if (fp == NULL)
         {
            printf("ERROR: cannot open file matlabiplt4m.m.\n");
            continue;
         }
         sprintf(pString," set cvFlag to 1 to use convex hull");
         fwriteComment(fp, pString);
         fprintf(fp, "cvFlag = 0;\n");
         sprintf(pString," change npart to change resolution of the ");
         fwriteComment(fp, pString);
         sprintf(pString," 4th dimension");
         fwriteComment(fp, pString);
         fwritePlotCLF(fp);
         fprintf(fp, "X = [\n");
         for (sInd = 0; sInd < nSamples_; sInd++)
         {
            dtemp = sampleOutputs_[sInd*nOutputs_];
            if (dtemp > 0.5*PSUADE_UNDEFINED || sampleStates_[sInd] != 1)
               fprintf(fp, "%e %e %e %e 0\n",
                       sampleInputs_[sInd*nInputs_+iplot1],
                       sampleInputs_[sInd*nInputs_+iplot2],
                       sampleInputs_[sInd*nInputs_+iplot3],
                       sampleInputs_[sInd*nInputs_+iplot4]);
            else
               fprintf(fp, "%e %e %e %e 1\n",
                       sampleInputs_[sInd*nInputs_+iplot1],
                       sampleInputs_[sInd*nInputs_+iplot2],
                       sampleInputs_[sInd*nInputs_+iplot3],
                       sampleInputs_[sInd*nInputs_+iplot4]);
         }
         fprintf(fp, "];\n");
         fprintf(fp, "iset = find(X(:,5)==0);\n");
         fprintf(fp, "m = size(iset,1);\n");
         fprintf(fp, "if m == 0\n");
         fprintf(fp, "   disp('No valid point.')\n");
         fprintf(fp, "else\n");
         fprintf(fp, "   if m > 100\n");
         fprintf(fp, "      npart = 10;\n");
         fprintf(fp, "   else\n");
         fprintf(fp, "      npart = m / 10;\n");
         fprintf(fp, "   end;\n");
         fprintf(fp, "   if npart == 0\n");
         fprintf(fp, "      npart = 1;\n");
         fprintf(fp, "   end;\n");
         fprintf(fp, "   xmax = max(X(:,4));\n");
         fprintf(fp, "   xmin = min(X(:,4));\n");
         fprintf(fp, "   for ii = 1 : npart\n");
         fprintf(fp, "      iset=find(X(:,4)<xmin+(xmax-xmin)/npart*ii);\n");
         fprintf(fp, "      XT = X(iset,:);\n");
         fprintf(fp, "      iset=find(XT(:,4)>=xmin+(xmax-xmin)/npart*(ii-1));\n");
         fprintf(fp, "      XV = XT(iset,:);\n");
         fprintf(fp, "      if cvFlag == 0\n");
         fprintf(fp, "        iset = find(XV(:,5)==0);\n");
         fprintf(fp, "        plot3(XV(iset,1),XV(iset,2),XV(iset,3),'r*',");
         fprintf(fp, "'MarkerSize',13)\n");
         fprintf(fp, "        hold on\n");
         fprintf(fp, "        iset = find(XV(:,5)~=0);\n");
         fprintf(fp, "        plot3(XV(iset,1),XV(iset,2),XV(iset,3),'bX',");
         fprintf(fp, "'MarkerSize',13)\n");
         fprintf(fp, "        hold off\n");
         fprintf(fp, "      else\n");
         fprintf(fp, "        iset = find(XV(:,5)==1);\n");
         fprintf(fp, "        mm = size(iset,1);\n");
         fprintf(fp, "        if (mm > 3)\n");
         fprintf(fp, "          XW = XV(iset,1:3);\n");
         fprintf(fp, "          KK = convhulln(XW);\n");
         fprintf(fp, "          trisurf(KK,XW(:,1),XW(:,2),XW(:,3))\n");
         fprintf(fp, "        else\n");
         fprintf(fp, "          disp('too few points to display')\n");
         fprintf(fp, "        end;\n");
         fprintf(fp, "      end;\n");
         fwritePlotXLabel(fp, inputNames_[iplot1]);
         fwritePlotYLabel(fp, inputNames_[iplot2]);
         fwritePlotZLabel(fp, inputNames_[iplot3]);
         fwritePlotAxes(fp);
         sprintf(pString,"4D Input Scatter Plot for %s)",inputNames_[iplot4]);
         fwritePlotTitle(fp, pString);
         fprintf(fp, "      disp('Press ENTER to advance')\n");
         fprintf(fp, "      pause\n");
         fprintf(fp, "   end;\n");
         fprintf(fp, "end;\n");
         fclose(fp);    
         printf("matlabiplt4m.m is now ready for input scatter plots.\n");
         printf("Note: see inside the matlab file to see options to \n");
         printf("      distinguish between valid and invalid runs.\n");
      }

      // +++ oplot2
      else if (!strcmp(command, "oplot2"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("oplot2: plot the sample outputs in two-parameter space.\n");
            printf("syntax: oplot2 (no argument needed)\n");
            printf("This command is useful for examining output correlation\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL)
         {
            printf("ERROR: data not loaded yet.\n");
            continue;
         }
         if (nOutputs_ < 2)
         {
            printf("ERROR: oplot2 requires 2 or more outputs.\n");
            continue;
         }
         oplot1 = oplot2 = -1;
         sprintf(pString,"Enter the output for x axis (1 - %d) : ",nOutputs_);
         oplot1 = getInt(1, nOutputs_, pString);
         oplot1--;
         if (nOutputs_ == 2)
         {
            if (oplot1 == 0) oplot2 = 1;
            else             oplot2 = 0;
         }
         while (oplot2 < 0 || oplot2 >= nOutputs_ || oplot1 == oplot2)
         {
            sprintf(pString,"Enter the output for y axis (1 - %d, not %d): ",
                    nOutputs_, oplot1+1);
            oplot2 = getInt(1, nOutputs_, pString);
            oplot2--;
            if (oplot2 == oplot1)
            {
               printf("ERROR: same index for x and y axes.\n");
               oplot2 = -1;
            }
         }

         if (psPlotTool_ == 1)
         {
            fp = fopen("scilaboplt2.sci", "w");
            if (fp == NULL)
            {
               printf("ERROR: cannot open file scilaboplt2.sci.\n");
               continue;
            }
            fwritePlotCLF(fp);
            fprintf(fp, "Y = [\n");
            for (sInd = 0; sInd < nSamples_; sInd++)
            {
               if ((sampleOutputs_[sInd*nOutputs_+oplot1]>0.5*PSUADE_UNDEFINED)||
                   (sampleOutputs_[sInd*nOutputs_+oplot2]>0.5*PSUADE_UNDEFINED))
                  fprintf(fp, "%e %e 0\n",
                          sampleOutputs_[sInd*nOutputs_+oplot1],
                          sampleOutputs_[sInd*nOutputs_+oplot2]);
               else
                  fprintf(fp, "%e %e 1\n",
                          sampleOutputs_[sInd*nOutputs_+oplot1],
                          sampleOutputs_[sInd*nOutputs_+oplot2]);
            }
            fprintf(fp, "];\n");
            fprintf(fp, "ia1 = find(X(:,4) == 0);\n");
            fprintf(fp, "ia2 = find(X(:,4) == 1);\n");
            fprintf(fp, "drawlater\n");
            fprintf(fp, "if (length(ia1) > 0)\n");
            fprintf(fp, "   plot(Y(ia1,1),Y(ia1,2),'r*')\n");
            fprintf(fp, "   set(gca(),\"auto_clear\",\"off\")\n");
            fprintf(fp, "end;\n");
            fprintf(fp, "if (length(ia2) > 0)\n");
            fprintf(fp, "   plot(Y(ia2,1),Y(ia2,2),'bX')\n");
            fprintf(fp, "end\n");
            fprintf(fp, "set(gca(),\"auto_clear\",\"on\")\n");
            fprintf(fp, "drawnow\n");
            fwritePlotXLabel(fp, outputNames_[oplot1]);
            fwritePlotYLabel(fp, outputNames_[oplot2]);
            fwritePlotAxes(fp);
            fwritePlotTitle(fp, "2D Output Scatter Plot");
            fclose(fp);    
            printf("scilaboplt2.sci is now ready for scatter plots.\n");
         }
         else
         {
            fp = fopen("matlaboplt2.m", "w");
            if (fp == NULL)
            {
               printf("ERROR: cannot open file matlaboplt2.m.\n");
               continue;
            }
            fwritePlotCLF(fp);
            fprintf(fp, "Y = [\n");
            for (sInd = 0; sInd < nSamples_; sInd++)
            {
               if ((sampleOutputs_[sInd*nOutputs_+oplot1]>0.5*PSUADE_UNDEFINED)||
                   (sampleOutputs_[sInd*nOutputs_+oplot2]>0.5*PSUADE_UNDEFINED))
                  fprintf(fp, "%e %e 0\n",
                          sampleOutputs_[sInd*nOutputs_+oplot1],
                          sampleOutputs_[sInd*nOutputs_+oplot2]);
               else
                  fprintf(fp, "%e %e 1\n",
                          sampleOutputs_[sInd*nOutputs_+oplot1],
                          sampleOutputs_[sInd*nOutputs_+oplot2]);
            }
            fprintf(fp, "];\n");
            fprintf(fp, "iset = find(Y(:,3) == 1);\n");
            fprintf(fp, "plot(Y(iset,1),Y(iset,2),\'bX\')\n");
            fprintf(fp, "hold on\n");
            fprintf(fp, "iset = find(Y(:,3) == 0);\n");
            fprintf(fp, "plot(Y(iset,1),Y(iset,2),\'r*\')\n");
            fprintf(fp, "hold off\n");
            fwritePlotXLabel(fp, outputNames_[oplot1]);
            fwritePlotYLabel(fp, outputNames_[oplot2]);
            fwritePlotAxes(fp);
            fwritePlotTitle(fp, "2D Output Scatter Plot");
            fclose(fp);    
            printf("matlaboplt2.m is now ready for scatter plots.\n");
         }
      }

      // +++ rsmcmc
      else if (!strcmp(command, "rsmcmc"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rsmcmc: perform MCMC on response surfaces\n");
            printf("syntax: rsmcmc (no argument needed)\n");
            printf("The response surface is constructed from the sample\n");
            printf("already loaded into local memory.\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL)
         {
            printf("ERROR: data not loaded yet.\n");
            printf("       Need input/output specification file.\n");
            continue;
         }
         analysisMethod = PSUADE_ANA_MCMC;
         anaManager = new AnalysisManager();
         anaManager->setup(analysisMethod, 0);
         psuadeIO_->getParameter("ana_diagnostics",pPtr);
         ii = pPtr.intData_;
         psuadeIO_->updateAnalysisSection(-1,-1,-1,outputLevel_,-1,-1);
         anaManager->analyze(psuadeIO_, 1, NULL, -1);
         psuadeIO_->updateAnalysisSection(-1,-1,-1,ii,-1,-1);
         delete anaManager;
      }

      // +++ mcmc
      else if (!strcmp(command, "mcmc"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("mcmc: perform MCMC on actual simulator\n");
            printf("syntax: mcmc (no argument needed)\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL)
         {
            printf("ERROR: data not loaded yet.\n");
            continue;
         }
         analysisMethod = PSUADE_ANA_MCMC;
         anaManager = new AnalysisManager();
         anaManager->setup(analysisMethod, 0);
         strcpy(winput, "setsim");
         targv[0] = winput;
         anaManager->specialRequest(analysisMethod, 1, targv);
         psuadeIO_->getParameter("ana_diagnostics",pPtr);
         ii = pPtr.intData_;
         psuadeIO_->updateAnalysisSection(-1,-1,-1,outputLevel_,-1,-1);
         anaManager->analyze(psuadeIO_, 1, NULL, outputID);
         psuadeIO_->updateAnalysisSection(-1,-1,-1,ii,-1,-1);
         delete anaManager;
      }

      // miscellaneous commands
      // +++ refine
      else if (!strcmp(command, "refine"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("refine: refine a sample (loaded previously using 'load')\n");
            printf("syntax: refine (no argument needed)\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data (load data first).\n");
            continue;
         }
         psuadeIO_->getParameter("method_sampling", pPtr);
         samplingMethod = pPtr.intData_;
         kk = psSamExpertMode_;
         psSamExpertMode_ = 0;
         sampPtr = (Sampling *) SamplingCreateFromID(samplingMethod);
         sampPtr->setPrintLevel(PL_INTERACTIVE);
         sampPtr->setInputBounds(nInputs_, iLowerB_, iUpperB_);
         sampPtr->setInputParams(nInputs_, NULL, NULL, NULL);
         sampPtr->setOutputParams(nOutputs_);
         psuadeIO_->getParameter("method_nreplications",pPtr);
         nReps = pPtr.intData_;
         sampPtr->setSamplingParams(nSamples_, nReps, -1);
         sampPtr->initialize(1);
         sampPtr->loadSamples(nSamples_, nInputs_, nOutputs_, sampleInputs_,
                              sampleOutputs_, sampleStates_);
         status = sampPtr->refine(2, 0, 0.0, nSamples_, NULL);

         if (status == 0)
         {
            if (sampleInputs_  != NULL) delete [] sampleInputs_;
            if (sampleOutputs_ != NULL) delete [] sampleOutputs_;
            if (sampleStates_  != NULL) delete [] sampleStates_;
            sampleInputs_  = NULL;
            sampleOutputs_ = NULL;
            sampleStates_  = NULL;
            nSamples_ = sampPtr->getNumSamples();
            sampleInputs_  = new double[nInputs_*nSamples_];
            sampleOutputs_ = new double[nOutputs_*nSamples_];
            sampleStates_  = new int[nSamples_];
            sampPtr->getSamples(nSamples_, nInputs_, nOutputs_, sampleInputs_,
                                sampleOutputs_, sampleStates_);
            psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,NULL,NULL,
                                     sampleInputs_,NULL,NULL,NULL,NULL,NULL); 
            psuadeIO_->updateOutputSection(nSamples_,nOutputs_,sampleOutputs_, 
                                           sampleStates_,NULL); 
            psuadeIO_->updateMethodSection(-1, nSamples_, -1, -1, -1);
            printf("refine successful: use write to store data set.\n");
            if (currSession != NULL) delete currSession;
            currSession = new PsuadeSession();
            psuadeIO_->getSession(currSession);
         }
         else
         {
            printf("ERROR: refine not successful.\n");
         }
         delete sampPtr;
         sampPtr = NULL;
         psSamExpertMode_ = kk;
      }

      // +++ a_refine_metis
      else if (!strcmp(command, "a_refine_metis"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("a_refine_metis: adaptively refine a loaded sample\n");
            printf("                when original sample is METIS.\n");
            printf("syntax: a_refine_metis (no argument needed)\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data (load data first).\n");
            continue;
         }
         outputID = 0;
         sprintf(pString,"Select output to use for a_refine (1 - %d) ",nOutputs_);
         outputID = getInt(1, nOutputs_, pString);
         psuadeIO_->getParameter("method_sampling", pPtr);
         samplingMethod = pPtr.intData_;
         if (samplingMethod != PSUADE_SAMP_METIS)
         {
            printf("ERROR: adaptive refinement requires METIS sample.\n");
            continue;
         }
         flag = 1;
         for (ss = 0; ss < nSamples_*nOutputs_; ss++)
         {
            if (sampleOutputs_[ss] != 1 && sampleOutputs_[ss] != 0)
            {
               flag = 0;
               break;
            }
         }
         if (flag == 0)
         {
            printf("This function uses the current loaded sample together\n");
            printf("with previous refinement information (stored in the file\n");
            printf("named 'psuadeMetisInfo') to perform adaptive refinement.\n");
            printf("To work properly, a few pieces of information is needed: \n");
            printf("(1) the original sample size (before any refinement)\n");
            printf("(2) the previous refinement file (psuadeMetisInfo) which\n");
            printf("    needs to be in the current directory.\n");
            sprintf(pString,"What is the original sample size ? ");
            int origNSamples = getInt(1, nSamples_, pString);
            sprintf(pString,"How many sample points to add? (1 - %d) ",nSamples_);
            int refineSize = getInt(1, nSamples_, pString);
            int tempSamExpert = psSamExpertMode_;
            psSamExpertMode_ = 0;
            sampPtr = (Sampling *) SamplingCreateFromID(samplingMethod);
            sampPtr->setPrintLevel(PL_INTERACTIVE);
            sampPtr->setInputBounds(nInputs_, iLowerB_, iUpperB_);
            sampPtr->setInputParams(nInputs_, NULL, NULL, NULL);
            sampPtr->setOutputParams(nOutputs_);
            sampPtr->setSamplingParams(nSamples_, -1, -1);
            sampPtr->initialize(1);
            sampPtr->loadSamples(nSamples_, nInputs_, nOutputs_, sampleInputs_,
                                 sampleOutputs_, sampleStates_);
            sampAux = (Sampling *) SamplingCreateFromID(samplingMethod);
            sampAux->setPrintLevel(PL_INTERACTIVE);
            sampAux->setInputBounds(nInputs_, iLowerB_, iUpperB_);
            sampAux->setInputParams(nInputs_, NULL, NULL, NULL);
            sampAux->setOutputParams(nOutputs_);
            sampAux->setSamplingParams(nSamples_, -1, -1);
            sampAux->initialize(1);
            if (nOutputs_ > 1)
            {
               tempY = new double[nSamples_];
               for (ss = 0; ss < nSamples_; ss++)
                  tempY[ss] = sampleOutputs_[ss*nOutputs_+outputID];
               sampAux->loadSamples(nSamples_, nInputs_, 1, sampleInputs_,
                                    tempY, sampleStates_);
            }
            else
            {
               sampAux->loadSamples(nSamples_, nInputs_, nOutputs_, sampleInputs_,
                                    sampleOutputs_, sampleStates_);
            }
            strcpy(pString, "setUniformRefinement");
            sampAux->setParam(pString);
            status = sampAux->refine(2, 1, 0.0, nSamples_, NULL);
            int nSamples2 = sampAux->getNumSamples();
            if (nSamples2 != 2 * nSamples_)
            {
               printf("a_refine_metis: Catastrophic ERROR: consult developers.\n");
               printf("         refined sample size != 2 * original size\n");
               printf("         May be due to too many levels of refinement.\n");
               delete sampAux;;
               sampAux = NULL;
               delete sampPtr;;
               sampPtr = NULL;
               if (nOutputs_ > 1) delete [] tempY;
               tempY = NULL;
               continue;
            }
            double *samInputs2  = new double[nSamples2*nInputs_];
            double *samOutputs2 = new double[nSamples2];
            int    *samStates2  = new int[nSamples2];
            double *samStds2    = new double[nSamples2];
            sampAux->getSamples(nSamples2, nInputs_, 1, samInputs2, 
                                samOutputs2, samStates2);
            faType = PSUADE_RS_MARSB;
            faPtr = genFA(faType, nInputs_, iOne, nSamples_);
            faPtr->setNPtsPerDim(32);
            faPtr->setBounds(iLowerB_, iUpperB_);
            faPtr->setOutputLevel(0);
            int numMars = 100, ivar1;
            double **marsX = new double*[numMars];
            double **marsY = new double*[numMars];
            for (ii = 0; ii < numMars; ii++)
            {
               marsX[ii] = new double[nInputs_*nSamples_];
               marsY[ii] = new double[nSamples_];
               for (ss = 0; ss < nSamples_; ss++)
               {
                  if (ss < origNSamples)
                       ivar1 = PSUADE_rand() % origNSamples;
                  else ivar1 = ss;
                  for (jj = 0; jj < nInputs_; jj++)
                     marsX[ii][ss*nInputs_+jj] =
                         sampleInputs_[ivar1*nInputs_+jj];
                  marsY[ii][ss] = sampleOutputs_[ivar1*nOutputs_+outputID];
               }
            }
            strcpy(cString, "mars_params");
            int ivar2 = 2 * nInputs_ / 3 + 1;
            int ivar3 = nSamples_;
            if (ivar3 > 300) ivar3 = 300;
            targv[0] = (char *) cString;
            targv[1] = (char *) &ivar3;
            targv[2] = (char *) &ivar2;
            faPtr->setParams(3, targv);
            strcpy(cString, "num_mars");
            targv[0] = (char *) cString;
            targv[1] = (char *) &numMars;
            faPtr->setParams(2, targv);
            strcpy(cString, "mars_sample");
            targv[0] = (char *) cString;
            targv[2] = (char *) &nSamples_;
            for (ii = 0; ii < numMars; ii++)
            {
               targv[1] = (char *) &ii;
               targv[3] = (char *) marsX[ii];
               targv[4] = (char *) marsY[ii];
               faPtr->setParams(5, targv);
            }
            if (nOutputs_ > 1)
                 faPtr->initialize(sampleInputs_,tempY);
            else faPtr->initialize(sampleInputs_,sampleOutputs_);

            faPtr->evaluatePointFuzzy(nSamples2-nSamples_, 
                            &samInputs2[nInputs_*nSamples_], 
                            &samOutputs2[nSamples_], &samStds2[nSamples_]);
            for (ss = 0; ss < nSamples_; ss++) 
               samStds2[ss] = PABS(samStds2[ss+nSamples_]);
            if (outputLevel_ > 4)
            {
               printf("Standard deviations to be used to select refinements.\n");
               for (ss = 0; ss < nSamples_; ss++) 
                  printf("Sample point %7d: stdev = %e\n", ss+1, samStds2[ss]);
            }

            strcpy(pString, "setAdaptiveRefinementBasedOnErrors");
            sampPtr->setParam(pString);
            sprintf(pString, "setRefineSize %d", refineSize);
            sampPtr->setParam(pString);
            sampPtr->refine(2,1,1.0e-6,nSamples_,samStds2);
            if (sampleInputs_  != NULL) delete [] sampleInputs_;
            if (sampleOutputs_ != NULL) delete [] sampleOutputs_;
            if (sampleStates_  != NULL) delete [] sampleStates_;

            nSamples_ = sampPtr->getNumSamples();
            sampleInputs_  = new double[nInputs_*nSamples_];
            sampleOutputs_ = new double[nOutputs_*nSamples_];
            sampleStates_  = new int[nSamples_];
            sampPtr->getSamples(nSamples_, nInputs_, nOutputs_, sampleInputs_,
                                sampleOutputs_, sampleStates_);
            psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,NULL,NULL,
                                       sampleInputs_,NULL,NULL,NULL,NULL,NULL); 
            psuadeIO_->updateOutputSection(nSamples_,nOutputs_,sampleOutputs_, 
                                           sampleStates_,NULL); 
            psuadeIO_->updateMethodSection(-1, nSamples_, -1, -1, -1);
            if (currSession != NULL) delete currSession;
            currSession = new PsuadeSession();
            psuadeIO_->getSession(currSession);
            printf("refine successful: use write to store the refined sample.\n");
            printf("                   Then run the newly created sample points\n");
            printf("                   in this refined sample.\n");
            delete sampPtr;
            delete sampAux;
            sampPtr = sampAux = NULL;
            delete faPtr;
            faPtr = NULL;
            psSamExpertMode_ = tempSamExpert;;
            delete [] samStds2;
            delete [] samInputs2;
            delete [] samOutputs2;
            delete [] samStates2;
            for (ii = 0; ii < numMars; ii++)
            {
               delete [] marsX[ii];
               delete [] marsY[ii];
            }
            delete [] marsX;
            delete [] marsY;
            if (nOutputs_ > 1) delete [] tempY;
         }
         else
         {
            printf("This function uses the current loaded sample together\n");
            printf("with previous refinement information (stored in the file\n");
            printf("named 'psuadeMetisInfo') to perform adaptive refinement.\n");
            printf("The sample has been detected to have 0/1 outputs.\n");
            printf("Adaptive refinement is applied to the 0/1 interfaces.\n");
            int tempSamExpert = psSamExpertMode_;
            psSamExpertMode_ = 0;
            sampPtr = (Sampling *) SamplingCreateFromID(samplingMethod);
            sampPtr->setPrintLevel(PL_INTERACTIVE);
            sampPtr->setInputBounds(nInputs_, iLowerB_, iUpperB_);
            sampPtr->setInputParams(nInputs_, NULL, NULL, NULL);
            sampPtr->setOutputParams(nOutputs_);
            sampPtr->setSamplingParams(nSamples_, -1, -1);
            sampPtr->initialize(1);
            if (nOutputs_ > 1)
            {
               tempY = new double[nSamples_];
               for (ss = 0; ss < nSamples_; ss++)
                  tempY[ss] = sampleOutputs_[ss*nOutputs_+outputID];
               sampPtr->loadSamples(nSamples_, nInputs_, 1, sampleInputs_,
                                    tempY, sampleStates_);
            }
            else
            {
               sampPtr->loadSamples(nSamples_, nInputs_, nOutputs_, sampleInputs_,
                                    sampleOutputs_, sampleStates_);
            }
            strcpy(pString, "setAdaptiveRefinementBasedOnOutputs");
            sampPtr->setParam(pString);
            sprintf(pString, "setRefineSize %d", nSamples_);
            sampPtr->setParam(pString);
            sampPtr->refine(2,1,1.0e-6,0,NULL);
            if (sampleInputs_  != NULL) delete [] sampleInputs_;
            if (sampleStates_  != NULL) delete [] sampleStates_;
            tempV = sampleOutputs_;

            kk = nSamples_;
            nSamples_ = sampPtr->getNumSamples();
            sampleInputs_  = new double[nInputs_*nSamples_];
            sampleOutputs_ = new double[nOutputs_*nSamples_];
            sampleStates_  = new int[nSamples_];
            if (nOutputs_ > 1)
            {
               sampPtr->getSamples(nSamples_, nInputs_, 1, sampleInputs_,
                                   tempY, sampleStates_);
               for (ss = 0; ss < kk*nOutputs_; ii++)
                  sampleOutputs_[ss] = tempV[ss];
               for (ss = kk*nOutputs_; ss < nSamples_*nOutputs_; ii++)
                  sampleOutputs_[ss] = PSUADE_UNDEFINED;
            }
            else
            {
               sampPtr->getSamples(nSamples_, nInputs_, nOutputs_, sampleInputs_,
                                   sampleOutputs_, sampleStates_);
            }
            psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,NULL,NULL,
                                          sampleInputs_,NULL,NULL,NULL,NULL,NULL); 
            psuadeIO_->updateOutputSection(nSamples_,nOutputs_,sampleOutputs_, 
                                           sampleStates_,NULL); 
            psuadeIO_->updateMethodSection(-1, nSamples_, -1, -1, -1);
            if (currSession != NULL) delete currSession;
            currSession = new PsuadeSession();
            psuadeIO_->getSession(currSession);
            printf("refine successful: use write to store the refined sample.\n");
            printf("                   Then run the newly created sample points\n");
            printf("                   in this refined sample.\n");
            delete tempV;
            delete sampPtr;
            sampPtr = NULL;
            if (nOutputs_ > 1) delete [] tempY;
            if (sampleOutputs_ != NULL) delete [] sampleOutputs_;
            psSamExpertMode_ = tempSamExpert;;
         }
      }

      // +++ a_refine
      else if (!strcmp(command, "a_refine"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("a_refine: adaptively refine a loaded sample\n");
            printf("syntax: a_refine (no argument needed)\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data (load data first).\n");
            continue;
         }
         outputID = 0;
         sprintf(pString,"Select output to use for a_refine (1 - %d) ",nOutputs_);
         outputID = getInt(1, nOutputs_, pString);
         outputID--;
         for (ss = 0; ss < nSamples_; ss++)
         {
            if (sampleOutputs_[ss*nInputs_+outputID] == PSUADE_UNDEFINED ||
                sampleStates_[ss] != 1)
            {
               printf("a_refine ERROR: some of the sample outputs are not ready.\n");
               break;
            }
         }
         if (ss != nSamples_) continue;
 
         printf("Please provide the initial sample size (before any refinement).\n");
         printf("below. This information is useful to prevent clustering of\n");
         printf("new sample points (all non-original sample points will be used\n");
         printf("in all bootstrapped samples for prediction error estimation).\n");
         printf("If you enter 0, the current sample size will be used as the\n");
         printf("initial sample size (i.e. no favorites).\n");
         sprintf(pString,"Enter the initial sample size (<= %d): ",nSamples_);
         int origNSamples = getInt(0, nSamples_, pString);
         if (origNSamples == 0) origNSamples = nSamples_;
         psuadeIO_->getParameter("method_sampling", pPtr);
         samplingMethod = PSUADE_SAMP_GMETIS;
         sprintf(pString,"How many sample points to add? (1 - %d) ",nSamples_);
         int refineSize = getInt(1, nSamples_, pString);
         int tempSamExpert = psSamExpertMode_;
         psSamExpertMode_ = 0;
         sampPtr = (Sampling *) SamplingCreateFromID(samplingMethod);
         sampPtr->setPrintLevel(PL_INTERACTIVE);
         sampPtr->setInputBounds(nInputs_, iLowerB_, iUpperB_);
         sampPtr->setOutputParams(nOutputs_);
         sampPtr->setSamplingParams(nSamples_, -1, -1);
         sampPtr->initialize(1);
         sampPtr->loadSamples(nSamples_, nInputs_, nOutputs_, sampleInputs_,
                              sampleOutputs_, sampleStates_);
         sampAux = (Sampling *) SamplingCreateFromID(samplingMethod);
         sampAux->setPrintLevel(PL_INTERACTIVE);
         sampAux->setInputBounds(nInputs_, iLowerB_, iUpperB_);
         sampAux->setOutputParams(1);
         sampAux->setSamplingParams(nSamples_, -1, -1);
         sampAux->initialize(1);
         if (nOutputs_ > 1)
         {
            tempY = new double[nSamples_];
            for (ss = 0; ss < nSamples_; ss++)
               tempY[ss] = sampleOutputs_[ss*nOutputs_+outputID];
            sampAux->loadSamples(nSamples_, nInputs_, 1, sampleInputs_,
                                 tempY, sampleStates_);
         }
         else
         {
            sampAux->loadSamples(nSamples_, nInputs_, nOutputs_, sampleInputs_,
                                 sampleOutputs_, sampleStates_);
         }
         strcpy(pString, "changeInfoName");
         sampAux->setParam(pString);
         strcpy(pString, "setUniformRefinement");
         sampAux->setParam(pString);
         status = sampAux->refine(2, 1, 0.0, nSamples_, NULL);
         int nSamples2 = sampAux->getNumSamples();
         double *samInputs2  = new double[nSamples2*nInputs_];
         double *samOutputs2 = new double[nSamples2];
         int    *samStates2  = new int[nSamples2];
         double *samStds2    = new double[nSamples2];
         sampAux->getSamples(nSamples2, nInputs_, 1, samInputs2, 
                             samOutputs2, samStates2);
         faType = PSUADE_RS_MARSB;
         if (psRSExpertMode_ == 1)
         {
            faType = -1;
            printf("Select a stochastic response surface method : ");
            printf("1. MARS with bootstrapped aggregation (MARSB)\n");
#ifdef HAVE_TPROS
            printf("2. Gaussian process (GP)\n");
#endif
            printf("3. Kriging\n");
            printf("4. Sum of trees method\n");
            sprintf(pString,"Make your choice (default: MARSB) : ");
            faType = getInt(1,4,pString);
            switch (faType)
            {
               case 1: faType = PSUADE_RS_MARSB; break;
#ifdef HAVE_TPROS
               case 2: faType = PSUADE_RS_GP1; break;
#else
               case 2: printf("GP not available. Defaulted to MARSB.\n"); 
                       faType = PSUADE_RS_MARSB; 
                       break;
#endif
               case 3: faType = PSUADE_RS_KR; break;
               case 4: faType = PSUADE_RS_SOTS; break;
            }
         }
         faPtr = genFA(faType, nInputs_, iOne, nSamples_);
         faPtr->setNPtsPerDim(32);
         faPtr->setBounds(iLowerB_, iUpperB_);
         faPtr->setOutputLevel(0);
         int numMars = 100, ivar1;
         double **marsX = new double*[numMars];
         double **marsY = new double*[numMars];
         for (ii = 0; ii < numMars; ii++)
         {
            marsX[ii] = new double[nInputs_*nSamples_];
            marsY[ii] = new double[nSamples_];
            for (ss = 0; ss < nSamples_; ss++)
            {
               if (ss < origNSamples)
                    ivar1 = PSUADE_rand() % origNSamples;
               else ivar1 = ss;
               for (jj = 0; jj < nInputs_; jj++)
                  marsX[ii][ss*nInputs_+jj] =
                      sampleInputs_[ivar1*nInputs_+jj];
               marsY[ii][ss] = sampleOutputs_[ivar1*nOutputs_+outputID];
            }
         }
         strcpy(cString, "mars_params");
         int ivar2 = 2 * nInputs_ / 3 + 1;
         int ivar3 = nSamples_;
         if (ivar3 > 300) ivar3 = 300;
         targv[0] = (char *) cString;
         targv[1] = (char *) &ivar3;
         targv[2] = (char *) &ivar2;
         faPtr->setParams(3, targv);
         strcpy(cString, "num_mars");
         targv[0] = (char *) cString;
         targv[1] = (char *) &numMars;
         faPtr->setParams(2, targv);
         strcpy(cString, "mars_sample");
         targv[0] = (char *) cString;
         targv[2] = (char *) &nSamples_;
         for (ii = 0; ii < numMars; ii++)
         {
            targv[1] = (char *) &ii;
            targv[3] = (char *) marsX[ii];
            targv[4] = (char *) marsY[ii];
            faPtr->setParams(5, targv);
         }
         if (nOutputs_ > 1)
              faPtr->initialize(sampleInputs_,tempY);
         else faPtr->initialize(sampleInputs_,sampleOutputs_);

         faPtr->evaluatePointFuzzy(nSamples2-nSamples_, 
                         &samInputs2[nInputs_*nSamples_], 
                         &samOutputs2[nSamples_], &samStds2[nSamples_]);

         for (ss = 0; ss < nSamples2-nSamples_; ss++) 
            samStds2[ss] = PABS(samStds2[ss+nSamples_]);
         if (outputLevel_ > 4)
         {
            printf("Standard deviations to be used to select refinements.\n");
               for (ss = 0; ss < nSamples2-nSamples_; ss++) 
                  printf("Sample point %7d: stdev = %e\n", ss+1, samStds2[ss]);
         }

         strcpy(pString, "setAdaptiveRefinementBasedOnErrors");
         sampPtr->setParam(pString);
         sprintf(pString, "setRefineSize %d", refineSize);
         sampPtr->setParam(pString);
         sampPtr->refine(2,1,1.0e-6,nSamples_,samStds2);
         if (sampleInputs_  != NULL) delete [] sampleInputs_;
         if (sampleOutputs_ != NULL) delete [] sampleOutputs_;
         if (sampleStates_  != NULL) delete [] sampleStates_;

         nSamples_ = sampPtr->getNumSamples();
         sampleInputs_  = new double[nInputs_*nSamples_];
         sampleOutputs_ = new double[nOutputs_*nSamples_];
         sampleStates_  = new int[nSamples_];
         sampPtr->getSamples(nSamples_, nInputs_, nOutputs_, sampleInputs_,
                             sampleOutputs_, sampleStates_);
         psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,NULL,NULL,
                                    sampleInputs_,NULL,NULL,NULL,NULL,NULL); 
         psuadeIO_->updateOutputSection(nSamples_,nOutputs_,sampleOutputs_, 
                                        sampleStates_,NULL); 
         psuadeIO_->updateMethodSection(-1, nSamples_, -1, -1, -1);
         if (currSession != NULL) delete currSession;
         currSession = new PsuadeSession();
         psuadeIO_->getSession(currSession);
         printf("refine successful: use write to store the refined sample.\n");
         printf("                   Then run the newly created sample points\n");
         printf("                   in this refined sample.\n");
         delete sampPtr;
         delete sampAux;
         sampPtr = sampAux = NULL;
         delete faPtr;
         faPtr = NULL;
         psSamExpertMode_ = tempSamExpert;;
         delete [] samStds2;
         delete [] samInputs2;
         delete [] samOutputs2;
         delete [] samStates2;
         for (ii = 0; ii < numMars; ii++)
         {
            delete [] marsX[ii];
            delete [] marsY[ii];
         }
         delete [] marsX;
         delete [] marsY;
         if (nOutputs_ > 1)
         {
            delete [] tempY;
            tempY = NULL;
         }
      }

      // +++ sqc
      else if (!strcmp(command, "sqc"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("sqc: check sample quality\n");
            printf("syntax: sqc (no argument needed)\n");
            printf("The sample quality is measured by max-min distance\n");
            printf("between sample points.\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data (load data first).\n");
         }
         else
         {
            tempV = new double[2];
            SamplingQuality(nSamples_,nInputs_,sampleInputs_,iLowerB_,
                            iUpperB_,tempV);
            printf("Max-min distance between sample points = %e (normalized)\n", 
                   tempV[0]);
            printf("Min-min distance between sample points = %e (normalized)\n", 
                   tempV[1]);
         }
      }

      // +++ validate 
      else if (!strcmp(command, "validate"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("validate: validate a subset of sample points (outputs)\n");
            printf("syntax: validate (no argument needed)\n");
            printf("This command just set the ready flag of each sample\n");
            printf("point to be 'ready' (NO OUTPUT VALUE IS MODIFIED).\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data (validate assumes outputs exist.)\n");
            continue;
         }
         sprintf(pString,"Enter first sample number (1 - %d) = ", nSamples_);
         ind = getInt(1, nSamples_, pString);
         ind--;
         sprintf(pString, "Enter last sample number (%d - %d) = ", 
                 ind+1, nSamples_);
         ind2 = getInt(ind+1, nSamples_, pString);
         for (ii = ind; ii < ind2; ii++) sampleStates_[ii] = 1;
         psuadeIO_->updateOutputSection(nSamples_,nOutputs_,sampleOutputs_, 
                                        sampleStates_,NULL); 
         psuadeIO_->writePsuadeFile(dataFile,0);
      }

      // +++ invalidate 
      else if (!strcmp(command, "invalidate"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("invalidate: select sample points to be unevaluated\n");
            printf("            based on a certain sample output value\n");
            printf("            (or deselect all sample points).\n");
            printf("syntax: invalidate (no argument needed)\n");
            printf("This command invalidates sample points based on the value\n");
            printf("of an output.  If most sample points are to be invalidated,\n");
            printf("first use 'invalidate' for all sample points, and then use\n");
            printf("'validate' to restore the range of desired sample points.\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data (load data first).\n");
            continue;
         }
         printf("Invalidate sample points based on the sample output values.\n");
         printf("If all sample points are to be invalidated, enter 0 below.\n");
         sprintf(pString, "Enter output number (1 - %d, 0 for all points) : ", 
                 nOutputs_);
         outputID = getInt(0, nOutputs_, pString);
         outputID--;
         if (outputID == -1)
         {
            for (ii = 0; ii < nSamples_; ii++) sampleStates_[ii] = 0;
            psuadeIO_->updateOutputSection(nSamples_,nOutputs_,sampleOutputs_, 
                                           sampleStates_,NULL); 
            printf("Use the write command to store the filtered sample.\n");
         }
         else
         {
            Ymax = - PSUADE_UNDEFINED;
            Ymin =   PSUADE_UNDEFINED;
            for (sInd = 0; sInd < nSamples_; sInd++)
            {
               if (sampleOutputs_[sInd*nOutputs_+outputID] > Ymax)
                  Ymax = sampleOutputs_[sInd*nOutputs_+outputID];
               if (sampleOutputs_[sInd*nOutputs_+outputID] < Ymin)
                  Ymin = sampleOutputs_[sInd*nOutputs_+outputID];
            }
            printf("INFO: Values outside the bounds are invalidated.\n");
            sprintf(pString, "Enter the filter's lower bound (Ymin=%e) : ",
                    Ymin);
            threshL = getDouble(pString);
            sprintf(pString, "Enter the filter's upper bound (Ymax=%e) : ",
                    Ymax);
            threshU = getDouble(pString);
            if (threshL >= threshU)
            {
               printf("ERROR: lower bound >= upper bound.\n");
               continue;
            }
            for (sInd = 0; sInd < nSamples_; sInd++)
            {
               if (sampleOutputs_[sInd*nOutputs_+outputID] < threshL ||
                   sampleOutputs_[sInd*nOutputs_+outputID] > threshU)
                  sampleStates_[sInd] = 0;
            }
            psuadeIO_->updateOutputSection(nSamples_, nOutputs_, sampleOutputs_, 
                                           sampleStates_, NULL);
            if (currSession != NULL) delete currSession;
            currSession = new PsuadeSession();
            psuadeIO_->getSession(currSession);
            printf("INFO: use 'purge' to take out the invalid points.\n");
            printf("INFO: then use 'write' to store the filtered sample.\n");
         }
      }

      // +++ srandomize 
      else if (!strcmp(command, "srandomize"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("srandomize: randomize the orders of the resident sample.\n");
            printf("syntax: srandomize <n>\n");
            printf("That command can be used, for example, with rstest_cv.\n");
            continue;
         }
         if (psuadeIO_ == NULL || nInputs_ <= 0)
         {
            printf("ERROR: data not loaded yet.\n");
            continue;
         }
         tempX = new double[nInputs_ * nSamples_];
         tempY = new double[nOutputs_ * nSamples_];
         tempI = new int[nSamples_];;
         generateRandomIvector(nSamples_, tempI);
         for (ii = 0; ii < nSamples_; ii++)
         {
            ind = tempI[ii];
            for (jj = 0; jj < nInputs_; jj++)
               tempX[ii*nInputs_+jj] = sampleInputs_[ind*nInputs_+jj];
            for (jj = 0; jj < nOutputs_; jj++)
               tempY[ii*nOutputs_+jj] = sampleOutputs_[ind*nOutputs_+jj];
            tempI[ii] = sampleStates_[ind];
         }
         delete [] sampleInputs_;
         delete [] sampleOutputs_;
         delete [] sampleStates_;
         sampleInputs_  = tempX;
         sampleOutputs_ = tempY;
         sampleStates_  = tempI;
      }

      // +++ spurge 
      else if (!strcmp(command, "spurge"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("spurge: remove bad sample points (Output=UNDEFINED or\n");
            printf("        status != 1). There is also an option to remove\n");
            printf("        repeated sample points.)\n");
            printf("syntax: spurge (no argument needed)\n");
            continue;
         }
         if (psuadeIO_ == NULL || nInputs_ <= 0)
         {
            printf("ERROR: data not loaded yet.\n");
            continue;
         }
         sprintf(pString,"Do you want duplicate sample points purged too? (y or n)");
         getString(pString, winput);
         setCompare = 0;
         if (winput[0] == 'y') setCompare = 1;
         kk = 0;
         for (sInd = 0; sInd < nSamples_; sInd++)
         {
            for (oInd = 0; oInd < nOutputs_; oInd++)
               if (sampleOutputs_[sInd*nOutputs_+oInd] == 
                   PSUADE_UNDEFINED || sampleStates_[sInd] == 0) break; 
            if (oInd == nOutputs_)
            {
               if (setCompare == 1)
                  jj = compareSamples(sInd,nSamples_,nInputs_,
                                      sampleInputs_, sampleStates_);
               else jj = -1;
               if (jj < 0 || jj > sInd)
               {
                  for (iInd = 0; iInd < nInputs_; iInd++)
                     sampleInputs_[kk*nInputs_+iInd] = 
                        sampleInputs_[sInd*nInputs_+iInd];
                  for (oInd = 0; oInd < nOutputs_; oInd++)
                     sampleOutputs_[kk*nOutputs_+oInd] = 
                        sampleOutputs_[sInd*nOutputs_+oInd];
                  sampleStates_[kk] = sampleStates_[sInd]; 
                  kk++;
               }
            }
         }
         nSamples_ = kk;
         samplingMethod = PSUADE_SAMP_MC;
         nReps = 1;
         psuadeIO_->updateMethodSection(samplingMethod,nSamples_,nReps,0,-1);
         psuadeIO_->updateInputSection(nSamples_, nInputs_, NULL, iLowerB_,
                          iUpperB_,sampleInputs_,NULL,NULL,NULL,NULL,NULL);
         psuadeIO_->updateOutputSection(nSamples_, nOutputs_, sampleOutputs_, 
                                    sampleStates_, outputNames_);
         if (currSession != NULL) delete currSession;
         currSession = new PsuadeSession();
         psuadeIO_->getSession(currSession);
         printf("Number of sample points after spurge = %d\n", nSamples_);
         printf("spurge completed. Use write to store the reduced sample.\n");
      }

      // +++ rm_dup
      else if (!strcmp(command, "rm_dup"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rm_dup: take out duplicate sample points\n");
            printf("syntax: rm_dup (no argument needed)\n");
            continue;
         }
         if (psuadeIO_ == NULL || nInputs_ <= 0)
         {
            printf("ERROR: data not loaded yet.\n");
            continue;
         }
         kk = 0;
         for (sInd = 0; sInd < nSamples_; sInd++)
         {
            jj = compareSamples(sInd,nSamples_,nInputs_,sampleInputs_,
                                sampleStates_);
            if (jj < 0 || jj > sInd)
            {
               for (iInd = 0; iInd < nInputs_; iInd++)
                  sampleInputs_[kk*nInputs_+iInd] = 
                        sampleInputs_[sInd*nInputs_+iInd];
               for (oInd = 0; oInd < nOutputs_; oInd++)
                  sampleOutputs_[kk*nOutputs_+oInd] = 
                        sampleOutputs_[sInd*nOutputs_+oInd];
               sampleStates_[kk] = sampleStates_[sInd]; 
               kk++;
            }
         }
         nSamples_ = kk;
         samplingMethod = PSUADE_SAMP_MC;
         nReps = 1;
         psuadeIO_->updateMethodSection(samplingMethod,nSamples_,nReps,0,-1);
         psuadeIO_->updateInputSection(nSamples_, nInputs_, NULL, iLowerB_,
                        iUpperB_, sampleInputs_, NULL,NULL,NULL,NULL,NULL);
         psuadeIO_->updateOutputSection(nSamples_, nOutputs_, sampleOutputs_, 
                                    sampleStates_, outputNames_);
         if (currSession != NULL) delete currSession;
         currSession = new PsuadeSession();
         psuadeIO_->getSession(currSession);
         printf("Number of sample points after rm_repeats = %d\n", nSamples_);
         printf("rm_dup completed. Use write to store the reduced sample.\n");
      }

      // +++ ifilter 
      else if (!strcmp(command, "ifilter"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("ifilter: take out sample points based on input constraints\n");
            printf("syntax: ifilter (no argument needed)\n");
            continue;
         }
         if (psuadeIO_ == NULL || nInputs_ <= 0)
         {
            printf("ERROR: data not loaded yet.\n");
            continue;
         }
         inputID = 0;
         if (nInputs_ > 1)
         {
            sprintf(pString, "Enter input number (1 - %d) : ", nInputs_);
            inputID = getInt(1, nInputs_, pString);
            inputID--;
         }
         Xmin = sampleInputs_[inputID];
         for (sInd = 1; sInd < nSamples_; sInd++)
            if (sampleInputs_[sInd*nInputs_+inputID] < Xmin) 
               Xmin = sampleInputs_[sInd*nInputs_+inputID];
         Xmax = sampleInputs_[inputID];
         for (sInd = 1; sInd < nSamples_; sInd++)
            if (sampleInputs_[sInd*nInputs_+inputID] > Xmax) 
               Xmax = sampleInputs_[sInd*nInputs_+inputID];
         printf("Xmin and Xmax found = %e %e.\n", Xmin, Xmax);
         sprintf(pString,"Enter the lower threshold (Xmin = %e) : ",Xmin);
         threshL = getDouble(pString);
         sprintf(pString,"Enter the upper threshold (Xmax = %e) : ",Xmax);
         threshU = getDouble(pString);
         if (threshL >= threshU)
         {
            printf("ERROR: Lower bound should be < upper bound.\n");
            continue;
         }
         kk = 0;
         for (sInd = 0; sInd < nSamples_; sInd++)
         {
            if (sampleInputs_[sInd*nInputs_+inputID] >= threshL && 
                sampleInputs_[sInd*nInputs_+inputID] <= threshU) 
            {
               for (iInd = 0; iInd < nInputs_; iInd++)
                  sampleInputs_[kk*nInputs_+iInd] = 
                     sampleInputs_[sInd*nInputs_+iInd];
               for (oInd = 0; oInd < nOutputs_; oInd++)
                  sampleOutputs_[kk*nOutputs_+oInd] = 
                     sampleOutputs_[sInd*nOutputs_+oInd];
               sampleStates_[kk] = sampleStates_[sInd]; 
               kk++;
            }
         }
         nSamples_ = kk;
         nReps = 1;
         psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,NULL,NULL,
                          sampleInputs_,inputNames_,NULL,NULL,NULL,NULL); 
         psuadeIO_->updateOutputSection(nSamples_, nOutputs_, sampleOutputs_, 
                                        sampleStates_, outputNames_);
         psuadeIO_->updateMethodSection(-1,nSamples_,nReps,0,-1);
         if (currSession != NULL) delete currSession;
         currSession = new PsuadeSession();
         psuadeIO_->getSession(currSession);
         printf("ifilter completed. Use write and load again to continue.\n");
      }

      // +++ ofilter 
      else if (!strcmp(command, "ofilter"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("ofilter: take out sample points based on output constraints\n");
            printf("syntax: ofilter (no argument needed)\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data (load data first).\n");
            continue;
         }
         if (nOutputs_ == 1) outputID = 0;
         else
         {
            sprintf(pString, "Enter output number (1 - %d) : ", nOutputs_);
            outputID = getInt(1, nOutputs_, pString);
            outputID--;
         }
         Ymax = - PSUADE_UNDEFINED;
         Ymin =   PSUADE_UNDEFINED;
         for (sInd = 0; sInd < nSamples_; sInd++)
         {
            if (sampleOutputs_[sInd*nOutputs_+outputID] > Ymax)
               Ymax = sampleOutputs_[sInd*nOutputs_+outputID];
            if (sampleOutputs_[sInd*nOutputs_+outputID] < Ymin)
               Ymin = sampleOutputs_[sInd*nOutputs_+outputID];
         }
         sprintf(pString, "Enter the lower constraint (Ymin=%e) : ",Ymin);
         threshL = getDouble(pString);
         sprintf(pString, "Enter the upper constraint (Ymax=%e) : ",Ymax);
         threshU = getDouble(pString);
         if (threshL >= threshU)
         {
            printf("ERROR: lower bound >= upper bound.\n");
            continue;
         }
         kk = 0;
         for (sInd = 0; sInd < nSamples_; sInd++)
         {
            if (sampleOutputs_[sInd*nOutputs_+outputID] >= threshL &&
                sampleOutputs_[sInd*nOutputs_+outputID] <= threshU)
            {
               for (iInd = 0; iInd < nInputs_; iInd++)
                  sampleInputs_[kk*nInputs_+iInd] = 
                        sampleInputs_[sInd*nInputs_+iInd];
               for (oInd = 0; oInd < nOutputs_; oInd++)
                  sampleOutputs_[kk*nOutputs_+oInd] = 
                        sampleOutputs_[sInd*nOutputs_+oInd];
               sampleStates_[kk] = sampleStates_[sInd]; 
               kk++;
            }
         }
         nSamples_ = kk;
         nReps = 1;
         psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,NULL,NULL,
                          sampleInputs_,inputNames_,NULL,NULL,NULL,NULL); 
         psuadeIO_->updateOutputSection(nSamples_, nOutputs_, sampleOutputs_, 
                                        sampleStates_, outputNames_);
         psuadeIO_->updateMethodSection(-1,nSamples_,nReps,0,-1);
         if (currSession != NULL) delete currSession;
         currSession = new PsuadeSession();
         psuadeIO_->getSession(currSession);
         printf("ofilter completed. Use write and load again to continue.\n");
      }

      // +++ oop 
      else if (!strcmp(command, "oop"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("oop: form linear combinations of output values for each point.\n");
            printf("syntax: oop (no argument needed)\n");
            printf("Sometimes one may want to combine some outputs to form a\n");
            printf("new output (e.g. the difference of two outputs).\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data (load data first).\n");
            continue;
         }
         printf("Set output <1> = <a> * output <2> + <b> * output <3>\n");
         sprintf(pString, "Enter output <1> (1 - %d) : ", nOutputs_);
         ii = getInt(1, nOutputs_, pString);
         ii--;
         sprintf(pString, "Enter the value <a>: ");
         aVal = getDouble(pString);
         sprintf(pString, "Enter output <2> (1 - %d) : ", nOutputs_);
         jj = getInt(1, nOutputs_, pString);
         jj--;
         sprintf(pString, "Enter the value <b>: ");
         bVal = getDouble(pString);
         sprintf(pString, "Enter output <3> (1 - %d) : ", nOutputs_);
         kk = getInt(1, nOutputs_, pString);
         kk--;
         for (sInd = 0; sInd < nSamples_; sInd++)
         {
            sampleOutputs_[sInd*nOutputs_+ii] = 
                        aVal * sampleOutputs_[sInd*nOutputs_+jj] + 
                        bVal * sampleOutputs_[sInd*nOutputs_+kk]; 
         }
         psuadeIO_->updateOutputSection(nSamples_, nOutputs_, sampleOutputs_, 
                                        sampleStates_, outputNames_);
         if (currSession != NULL) delete currSession;
         currSession = new PsuadeSession();
         psuadeIO_->getSession(currSession);
         printf("oop completed. Use write and load again to continue.\n");
      }

      // +++ nna
      else if (!strcmp(command, "nna"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("nna: nearest neighbor analysis (for detecting outliers)\n");
            printf("syntax: nna (no argument needed)\n");
            printf("For each sample point, the output will be compared to\n");
            printf("its nearest neighbor. If the sample point is an outlier,\n");
            printf("it will be shown to have large gradient.\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data (load data first).\n");
            continue;
         }
         outputID = 0;
         printf("nna: for each data point, find nearest neighbor and");
         printf(" plot changes in the outputs.\n");
         sprintf(pString, "Enter output number (1 - %d) = ", nOutputs_);
         outputID = getInt(1, nOutputs_, pString);
         outputID--;
         if (psPlotTool_ == 1)
         {
            fp = fopen("scilabnna.sci", "w");
            if (fp == NULL)
            {
               printf("ERROR: cannot open file scilabnna.sci.\n");
               continue;
            }
         }
         else
         {
            fp = fopen("matlabnna.m", "w");
            if (fp == NULL)
            {
               printf("ERROR: cannot open file matlabnna.m.\n");
               continue;
            }
         }
         sprintf(pString," nearest neighbor analysis");
         fwriteComment(fp, pString);
         sprintf(pString," for detecting outliers.");
         fwriteComment(fp, pString);
         sprintf(pString," The following plot is : ");
         fwriteComment(fp, pString);
         sprintf(pString," Y-axis: delta output / distance");
         fwriteComment(fp, pString);
         sprintf(pString," X-axis: sample number");
         fwriteComment(fp, pString);
         fprintf(fp, "A = [\n");
         for (ii = 0; ii < nSamples_; ii++)
         {
            minDist = PSUADE_UNDEFINED;
            ind = -1;
            for (kk = 0; kk < nSamples_; kk++)
            {
               if (ii != kk)
               {
                  ddata = 0.0;
                  for (jj = 0; jj < nInputs_; jj++)
                  {
                     dtemp = sampleInputs_[ii*nInputs_+jj] - 
                             sampleInputs_[kk*nInputs_+jj];
                     dtemp /= (iUpperB_[jj] - iLowerB_[jj]);
                     ddata += pow(dtemp, 2.0);
                  }
                  if (ddata < minDist && ddata > 0.0)
                  {
                     minDist = ddata;
                     ind = kk;
                  }
               }
            }
            if (minDist > 0)
               dtemp = (sampleOutputs_[ii*nOutputs_+outputID] -
                        sampleOutputs_[ind*nOutputs_+outputID])/minDist;
            else dtemp = 0.0;
            fprintf(fp, "%7d %24.16e %d %e\n", ii+1, dtemp, ind+1, minDist);
         }
         fprintf(fp, "];\n");
         if (psPlotTool_ == 1)
              fprintf(fp, "plot(A(:,1),A(:,2),'x');\n");
         else fprintf(fp, "plot(A(:,1),A(:,2),'x','MarkerSize',12)\n");
         fwritePlotXLabel(fp, "Sample number");
         fwritePlotYLabel(fp, "Delta output / distance");
         fwritePlotTitle(fp, "Nearest Neighbor Analysis");
         fwritePlotAxes(fp);
         fclose(fp);
         if (psPlotTool_ == 1)
              printf("Nearest neighbors analysis is now in file scilabnna.sci\n");
         else printf("Nearest neighbors analysis is not in file matlabnna.m.\n");
      }

      // +++ interface_track
      else if (!strcmp(command, "interface_track"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("interface_track: interface (0/1) tracking \n");
            printf("syntax: interface_track (no argument needed)\n");
            printf("This command tracks the interface between the region\n");
            printf("where the output is '1' versus the region where the\n");
            printf("output is '0' using a polynomial function.\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data (load data first).\n");
            continue;
         }
         outputID = 0;
         sprintf(pString, "Which output to use (1 - %d) = ", nOutputs_);
         outputID = getInt(1, nOutputs_, pString);
         outputID--;
         int nIndex = 3, ss2;
         if (psAnaExpertMode_ == 1)
         {
            printf("In order to track the interface, the outputs of the\n");
            printf("nearest neighbors of each sample point are compared\n");
            printf("against themselves. You need to select the number of\n");
            printf("neighbors (K) to examine. The more neighbors are used,\n");
            printf("the more fuzzy the interface will be. On the other hand,\n");
            printf("the less neighbors are used, the sample size may be too\n");
            printf("small to be useful. Try different ones. The default is 3.\n");
            sprintf(pString, "What is K (>= 1, <= 10, default=3)? ");
            nIndex = getInt(1, 10, pString);
         }
         double *distPairs  = new double[nSamples_ * (nSamples_ - 1) / 2];
         int    *minIndices = new int[nIndex];;
         for (ss = 0; ss < nSamples_; ss++)
         {
            for (ss2 = 0; ss2 < ss; ss2++)
            {
               kk = ss * (ss - 1) / 2 + ss2;
               distPairs[kk] = 0.0;
               for (jj = 0; jj < nInputs_; jj++)
               {
                  ddata = sampleInputs_[ss*nInputs_+jj] - sampleInputs_[ss2*nInputs_+jj];
                  ddata = ddata / (iUpperB_[jj] - iLowerB_[jj]);
                  distPairs[kk] += pow(ddata, 2.0);
               }
            }
         } 
         count = 0;
         tempX = new double[nSamples_*nInputs_];
         for (ss = 0; ss < nSamples_; ss++)
         {
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
               for (ss2 = ss+1; ss2 < nSamples_; ss2++)
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
            }
            flag = 0;
            for (jj = 0; jj < nIndex; jj++)
            {
               if (minIndices[jj] != -1)
               {
                  if (sampleOutputs_[ss*nOutputs_+outputID] != 
                      sampleOutputs_[minIndices[jj]*nOutputs_+outputID]) 
                      flag = 1;
               }
            }
            if (flag == 1) 
            {
               for (jj = 0; jj < nInputs_; jj++) 
                  tempX[count*nInputs_+jj] = sampleInputs_[ss*nInputs_+jj];
               count++;
            }
         }
         delete [] minIndices;
         delete [] distPairs;
         if (count < nInputs_+1)
         {
            printf("Too few valid points (%d).\n", count);
         }
         else
         {
            if (outputLevel_ > 3)
            {
               printf("Listed below are the sample points on the interface.\n");
               for (ss = 0; ss < count; ss++) 
               {
                  printf("%5d  ", ss+1);
                  for (ii = 0; ii < nInputs_; ii++) 
                     printf("%e ", tempX[ss*nInputs_+ii]);
                  printf("\n");
               }
            }
            ioPtr = new PsuadeData();
            printf("Since the interface is at most m-1 dimension where\n");
            printf("m is the number of inputs. Please select one input\n");
            printf("to be the dependent variable in the polynomial equation\n");
            printf("describing the interface in the form of :\n");
            printf("   X(k) = a polynomial in terms of all other X(i), i != k.\n");
            sprintf(pString, "Select k (1 - %d) : ", nInputs_);
            kk = getInt(1, nInputs_, pString);
            kk--;
            double *jLowerB = new double[nInputs_-1];
            double *jUpperB = new double[nInputs_-1];
            jj = 0;
            for (ii = 0; ii < nInputs_; ii++)
            {
               if (ii != kk)
               {
                  jLowerB[jj] = iLowerB_[ii];
                  jUpperB[jj] = iUpperB_[ii];
                  jj++;
               }
            }
            tempY  = new double[count];
            states = new int[count];
            for (ss = 0; ss < count; ss++) 
            {
               tempY[ss] = tempX[ss*nInputs_+kk];
               jj = 0;
               for (ii = 0; ii < nInputs_; ii++)
               {
                  if (ii != kk)
                  {
                     tempX[ss*(nInputs_-1)+jj] = tempX[ss*nInputs_+ii];
                     jj++;
                  }
               }
            }
            int    *iPDFs = new int[nInputs_];
            double *iMeans = new double[nInputs_];
            double *iStds = new double[nInputs_];
            Matrix *iCMat = new Matrix();
            iCMat->setDim(nInputs_, nInputs_);
            for (ii = 0; ii < nInputs_; ii++)
            {
               iPDFs[ii] = 0;
               iMeans[ii] = 0;
               iStds[ii] = 0;
               iCMat->setEntry(ii,ii,1.0);
            }
            ioPtr->updateInputSection(count, nInputs_-1, NULL, jLowerB,
                          jUpperB, tempX, inputNames_,NULL,NULL,NULL,NULL); 
            delete [] iPDFs;
            delete [] iMeans;
            delete [] iStds;
            delete iCMat;
            for (ss = 0; ss < count; ss++) states[ss] = 1;
            ioPtr->updateOutputSection(count, 1, tempY, states,outputNames_);
            ioPtr->updateMethodSection(PSUADE_SAMP_MC, count, 1, -1, -1);
            analysisMethod = PSUADE_ANA_RSFA;
            anaManager = new AnalysisManager();
            printf("Select which polynomial regression to use: \n");
            printf("(1) linear \n");
            printf("(2) quadratic \n");
            printf("(3) cubic \n");
            printf("(4) quartic \n");
            sprintf(pString, "Please select (1 - 4) : ");
            ind = getInt(1, 4, pString);
            faType = PSUADE_RS_REGR1 + ind - 1;
            anaManager->setup(analysisMethod, faType);
            ioPtr->updateAnalysisSection(-1,-1,-1,1,-1,-1);
            anaManager->analyze(ioPtr, 1, NULL, 0);
            delete anaManager;
            delete ioPtr;
            ioPtr = NULL;
            delete [] tempY;
            delete [] states;
            delete [] jLowerB;
            delete [] jUpperB;
         }
         delete [] tempX;
      }

      // +++ pdfconvert 
      else if (!strcmp(command, "pdfconvert"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("pdfconvert: convert a sample based on its pdfs\n");
            printf("syntax: pdfconvert (no argument needed)\n");
            printf("To use this, first create a sample using uniform\n");
            printf("distribution. Then load the sample (make sure that\n");
            printf("before you load the sample, you have modified the\n");
            printf("input section of this sample file to reflect the\n");
            printf("desired distribution. Then use this command to\n");
            printf("convert the sample to the desired distributions.\n");
            printf("Use write to store the converted sample.\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data (load data first).\n");
            continue;
         }
         pdfman = new PDFManager();
         psuadeIO_->getParameter("ana_diagnostics",pPtr);
         kk = pPtr.intData_;
         psuadeIO_->updateAnalysisSection(-1,-1,-1,outputLevel_,-1, -1);
         psuadeIO_->updateMethodSection(PSUADE_SAMP_MC, -1, -1, -1, -1);
         pdfman->initialize(psuadeIO_);
         vecIn.load(nSamples_*nInputs_, sampleInputs_);
         vecOut.setLength(nSamples_*nInputs_);
         vecUpper.load(nInputs_, iUpperB_);
         vecLower.load(nInputs_, iLowerB_);
         pdfman->invCDF(nSamples_, vecIn, vecOut, vecLower, vecUpper);
         psuadeIO_->updateAnalysisSection(-1,-1,-1,ii,-1, -1);
         for (ii = 0; ii < nSamples_*nInputs_; ii++)
            sampleInputs_[ii] = vecOut[ii];
         int    *iPDFs = new int[nInputs_];
         double *iMeans = new double[nInputs_];
         double *iStds = new double[nInputs_];
         Matrix *iCMat = new Matrix();
         iCMat->setDim(nInputs_, nInputs_);
         for (ii = 0; ii < nInputs_; ii++)
         {
            iPDFs[ii] = 0;
            iMeans[ii] = 0;
            iStds[ii] = 0;
            iCMat->setEntry(ii,ii,1.0);
         }
         psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,
                                       NULL,NULL,sampleInputs_,NULL, 
                                       iPDFs,iMeans,iStds,iCMat);
         delete [] iPDFs;
         delete [] iMeans;
         delete [] iStds;
         delete iCMat;
         psuadeIO_->updateAnalysisSection(-1, -1, -1, kk, -1, 0);
         if (currSession != NULL) delete currSession;
         currSession = new PsuadeSession();
         psuadeIO_->getSession(currSession);
         delete pdfman;
         printf("The sample in local memory has been converted based on\n");
         printf("the PDF information in the INPUT section.\n");
         printf("You can now store your new sample using 'write'.\n");
      }

      // +++ rand_draw 
      else if (!strcmp(command, "rand_draw"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rand_draw: draw a random sample from the loaded sample\n");
            printf("syntax: rand_draw <n>\n");
            printf("That is, this command generates a bootstrapped sample\n");
            printf("from the sample which has been loaded to local memory.\n");
            printf("This command is used if your want to draw a sample\n");
            printf("from the posterior sample after Bayesian analysis.\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data (load data first).\n");
            continue;
         }
         sprintf(pString,"Size of the sample to be drawn : (1-1000000) ");
         count = getInt(1, 1000000, pString);
         tempX = new double[nInputs_ * count];
         tempY = new double[nOutputs_ * count];
         tempI = new int[count];;
         for (ii = 0; ii < count; ii++)
         {
            ind = PSUADE_rand() % nSamples_;
            for (jj = 0; jj < nInputs_; jj++)
               tempX[ii*nInputs_+jj] = sampleInputs_[ind*nInputs_+jj];
            for (jj = 0; jj < nOutputs_; jj++)
               tempY[ii*nOutputs_+jj] = sampleOutputs_[ind*nOutputs_+jj];
            tempI[ii] = sampleStates_[ind];
         }
         int    *iPDFs = new int[nInputs_];
         double *iMeans = new double[nInputs_];
         double *iStds = new double[nInputs_];
         Matrix *iCMat = new Matrix();
         iCMat->setDim(nInputs_, nInputs_);
         for (ii = 0; ii < nInputs_; ii++)
         {
            iPDFs[ii] = 0;
            iMeans[ii] = 0;
            iStds[ii] = 0;
            iCMat->setEntry(ii,ii,1.0);
         }
         ioPtr = new PsuadeData();
         ioPtr->updateInputSection(count, nInputs_, NULL, iLowerB_,
                                   iUpperB_, tempX, inputNames_,
                                   iPDFs,iMeans,iStds,iCMat);
         delete [] iPDFs;
         delete [] iMeans;
         delete [] iStds;
         delete iCMat;
         ioPtr->updateOutputSection(count, nOutputs_, tempY, tempI,
                                    outputNames_);
         ioPtr->updateMethodSection(PSUADE_SAMP_MC, count, 1, -1, -1);
         printf("Store random sample to : (filename) ");
         scanf("%s", dataFile);
         fgets(lineIn,5000,stdin); 
         if ((fp = fopen(dataFile, "w")) == NULL)
         {
            printf("ERROR: cannot open file %s\n", dataFile);
         }
         else
         {
            ioPtr->writePsuadeFile(dataFile,0);
            printf("rand_draw completed. Sample has been saved.\n");
         }    
         delete [] tempX;
         delete [] tempY;
         delete [] tempI;
         delete ioPtr;
      }

      // +++ rand_draw2 
      else if (!strcmp(command, "rand_draw2"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rand_draw2: draw a sample randomly from 2 samples\n");
            printf("syntax: rand_draw2 <n>\n");
            printf("That is, this command generates a bootstrapped sample\n");
            printf("from two PSUADE-formatted sample files. This command\n");
            printf("is used if you have 2 posterior samples from MCMC for\n");
            printf("2 sets of inputs, and you would like to propagate them\n");
            printf("through the coupled model which has these two sets of\n");
            printf("inputs.\n");
            continue;
         }
         if (psuadeIO_ != NULL) delete psuadeIO_;

         sprintf(pString,"Enter name of the first file : ");
         getString(pString, winput);
         kk = strlen(winput);
         winput[kk-1] = '\0';
         psuadeIO_ = new PsuadeData;
         status = psuadeIO_->readPsuadeFile(winput);
         if (status != 0)
         {
            printf("rand_draw2 FILE READ ERROR: file = %s\n", winput);
            continue;
         }
         sprintf(pString,"Enter name of the second file : ");
         getString(pString, winput);
         kk = strlen(winput);
         winput[kk-1] = '\0';
         ioPtr = new PsuadeData();
         status = ioPtr->readPsuadeFile(winput);
         if (status != 0)
         {
            printf("rand_draw2 FILE READ ERROR: file = %s\n", winput);
            continue;
         }

         psuadeIO_->getParameter("input_ninputs", pPtr);
         nInputs_ = pPtr.intData_;
         psuadeIO_->getParameter("method_nsamples", pPtr);
         nSamples_ = pPtr.intData_;
         psuadeIO_->getParameter("input_sample", pPtr);
         if (sampleInputs_ != NULL) delete [] sampleInputs_;
         sampleInputs_ = pPtr.dbleArray_;
         pPtr.dbleArray_ = NULL;

         ioPtr->getParameter("input_ninputs", pPtr);
         int nInputs2 = pPtr.intData_;
         ioPtr->getParameter("method_nsamples", pPtr);
         int nSamples2 = pPtr.intData_;
         ioPtr->getParameter("input_sample", pPtr);
         tempX = pPtr.dbleArray_;
         pPtr.dbleArray_ = NULL;

         iLowerB_ = new double[nInputs_+nInputs2];
         iUpperB_ = new double[nInputs_+nInputs2];
         pLower.clean();
         psuadeIO_->getParameter("input_lbounds", pLower);
         tempT = pLower.dbleArray_;
         pLower.dbleArray_ = NULL;
         for (ii = 0; ii < nInputs_; ii++) iLowerB_[ii] = tempT[ii];
         delete [] tempT;
         ioPtr->getParameter("input_lbounds", pLower);
         tempT = pLower.dbleArray_;
         pLower.dbleArray_ = NULL;
         for (ii = 0; ii < nInputs2; ii++) iLowerB_[nInputs_+ii] = tempT[ii];
         delete [] tempT;

         pUpper.clean();
         psuadeIO_->getParameter("input_ubounds", pUpper);
         tempT = pUpper.dbleArray_;
         pUpper.dbleArray_ = NULL;
         for (ii = 0; ii < nInputs_; ii++) iUpperB_[ii] = tempT[ii];
         delete [] tempT;
         ioPtr->getParameter("input_ubounds", pUpper);
         tempT = pUpper.dbleArray_;
         pUpper.dbleArray_ = NULL;
         for (ii = 0; ii < nInputs2; ii++) iUpperB_[nInputs_+ii] = tempT[ii];
         delete [] tempT;

         sprintf(pString,"Size of the sample to be drawn : (1-2000000) ");
         count = getInt(1, 200000, pString);
         tempW = new double[count*(nInputs_+nInputs2)];
         tempY = new double[count];
         tempI = new int[count]; 
         for (ii = 0; ii < count; ii++)
         {
            ind  = PSUADE_rand() % nSamples_;
            ind2 = PSUADE_rand() % nSamples2;
            for (jj = 0; jj < nInputs_; jj++)
               tempW[ii*(nInputs_+nInputs2)+jj] = 
                  sampleInputs_[ind*nInputs_+jj];
            for (jj = 0; jj < nInputs2; jj++)
               tempW[ii*(nInputs_+nInputs2)+nInputs_+jj] = 
                  tempX[ind2*nInputs2+jj];
            tempY[ii] = PSUADE_UNDEFINED;
            tempI[ii] = 0;
         } 

         inputNames_ = new char*[nInputs_+nInputs2];
         psuadeIO_->getParameter("input_names", pINames);
         for (ii = 0; ii < nInputs_; ii++)
         {
            inputNames_[ii] = new char[1001]; 
            strcpy(inputNames_[ii], pINames.strArray_[ii]);
         }
         pINames.clean();
         ioPtr->getParameter("input_names", pINames);
         for (ii = 0; ii < nInputs2; ii++)
         {
            inputNames_[nInputs_+ii] = new char[1001]; 
            strcpy(inputNames_[nInputs_+ii], pINames.strArray_[ii]);
         }
         pINames.clean();
         delete psuadeIO_;
         delete ioPtr;
         psuadeIO_ = NULL;

         ioPtr = new PsuadeData();
         int    *iPDFs = new int[nInputs_+nInputs2];
         double *iMeans = new double[nInputs_+nInputs2];
         double *iStds = new double[nInputs_+nInputs2];
         Matrix *iCMat = new Matrix();
         iCMat->setDim(nInputs_+nInputs2, nInputs_+nInputs2);
         for (ii = 0; ii < nInputs_+nInputs2; ii++)
         {
            iPDFs[ii] = 0;
            iMeans[ii] = 0;
            iStds[ii] = 0;
            iCMat->setEntry(ii,ii,1.0);
         }
         ioPtr->updateInputSection(count, nInputs_+nInputs2, NULL, 
                              iLowerB_, iUpperB_, tempW, inputNames_,
                              iPDFs, iMeans, iStds, iCMat);
         delete [] iPDFs;
         delete [] iMeans;
         delete [] iStds;
         delete iCMat;
 
         nOutputs_ = 1;
         outputNames_ = new char*[1];
         outputNames_[0] = new char[200]; 
         sprintf(outputNames_[0], "Y");
         ioPtr->updateOutputSection(count, nOutputs_, tempY, tempI,
                                    outputNames_);
         ioPtr->updateMethodSection(PSUADE_SAMP_MC, count, 1, -1, -1);

         printf("Store random sample to : (filename) ");
         scanf("%s", dataFile);
         fgets(lineIn,5000,stdin); 
         if ((fp = fopen(dataFile, "w")) == NULL)
         {
            printf("ERROR: cannot open file %s\n", dataFile);
         }
         else
         {
            ioPtr->writePsuadeFile(dataFile,0);
            printf("rand_draw2 completed. Sample has been saved.\n");
         }    
         delete [] tempW;
         delete [] tempX;
         delete [] tempY;
         delete [] tempI;
         for (ii = 0; ii < nInputs_+nInputs2; ii++) 
            delete [] inputNames_[ii];
         delete [] inputNames_;
         delete [] outputNames_[0];
         delete [] outputNames_;
         inputNames_  = NULL;
         outputNames_ = NULL;
         if (sampleInputs_  != NULL) delete [] sampleInputs_;
         if (sampleOutputs_ != NULL) delete [] sampleOutputs_;
         if (sampleStates_  != NULL) delete [] sampleStates_;
         if (iLowerB_       != NULL) delete [] iLowerB_;
         if (iUpperB_       != NULL) delete [] iUpperB_;
         if (inputPDFs_     != NULL) delete [] inputPDFs_;
         if (inputMeans_    != NULL) delete [] inputMeans_;
         if (inputStds_     != NULL) delete [] inputStds_;
         sampleInputs_  = NULL;
         sampleOutputs_ = NULL;
         sampleStates_  = NULL;
         iLowerB_       = NULL;
         iUpperB_       = NULL;
         inputPDFs_     = NULL;
         inputMeans_    = NULL;
         inputStds_     = NULL;
         nSamples_ = 0;
         nInputs_ = 0;
         nOutputs_ = 0;
         inputNames_ = NULL;
         outputNames_ = NULL;
         delete ioPtr;
         ioPtr = NULL;
      }

      // +++ gensample 
      else if (!strcmp(command, "gensample"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("gensample: create a sample from the loaded INPUT PDFs.\n");
            continue;
         }
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no sample input information.\n");
            printf("       Use load to put data into local memory first.\n");
            continue;
         }
         sprintf(pString,"Sample size? (>= 2, <=10000000) ");
         kk = getInt(1,10000000,pString);
         pdfman = new PDFManager();
         pdfman->initialize(nInputs_,inputPDFs_,inputMeans_,inputStds_,
                            *inputCMat_,inputPDFFiles,SPDFIndices);
         vecLower.load(nInputs_, iLowerB_);
         vecUpper.load(nInputs_, iUpperB_);
         vecIn.setLength(kk*nInputs_);
         pdfman->genSample(kk, vecIn, vecLower, vecUpper);
         delete pdfman;
         ioPtr = new PsuadeData();
         tempX = vecIn.getDVector();
         tempY = new double[kk];
         for (ii = 0; ii < kk; ii++) tempY[ii] = PSUADE_UNDEFINED; 
         states = new int[kk];
         for (ii = 0; ii < kk; ii++) states[ii] = 0; 
         ioPtr->updateInputSection(kk,nInputs_,NULL,iLowerB_,iUpperB_,tempX, 
                  inputNames_,inputPDFs_,inputMeans_,inputStds_,inputCMat_);
         if (sampleOutputs_ == NULL)
         {
            names = new char*[1];
            names[0] = new char[100];
            strcpy(names[0], "Y");
         }
         else names = outputNames_;
         ioPtr->updateOutputSection(kk, iOne, tempY, states, names);
         if (sampleOutputs_ == NULL)
         {
            delete [] names[0];
            delete [] names;
            names = NULL;
         }
         ioPtr->updateMethodSection(PSUADE_SAMP_MC, kk, 1, -1, -1);
         ioPtr->writePsuadeFile("psuade_sample",0);
         printf("The sample has been written into the 'psuade_sample' file.\n");
         delete ioPtr;
         delete [] tempY;
         delete [] states;
         tempY = NULL;
         states = NULL;
      }

      // +++ cdf_lookup 
      else if (!strcmp(command, "cdf_lookup"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("cdf_lookup: look up the cumulative probability given a value.\n");
            printf("syntax: cdf_lookup <n>\n");
            printf("Given a distribution type and its parameters, this command\n");
            printf("returns the cumulative probability when provided with a\n");
            printf("lookup value.\n");
            continue;
         }
         int    gtype;
         double slbound, subound, smean, sstdev;
         sprintf(pString,
           "PDF type? (1) N (2) L (3) T (4) Beta (5) Weibull (6) Gamma (7) Exp (8) F : ");
         gtype = getInt(1, 8, pString);
         if (gtype == 1)
         {
            gtype = PSUADE_PDF_NORMAL;
            sprintf(pString, "PDF mean = ");
            smean = getDouble(pString);
            sprintf(pString, "PDF std. dev. = ");
            sstdev = getDouble(pString);
         }
         else if (gtype == 2) 
         {
            gtype = PSUADE_PDF_LOGNORMAL;
            sprintf(pString, "PDF log(mean) = ");
            smean = getDouble(pString);
            sprintf(pString, "PDF std. dev. = ");
            sstdev = getDouble(pString);
         }
         else if (gtype == 3)
         {
            gtype = PSUADE_PDF_TRIANGLE;
            sprintf(pString, "PDF center = ");
            smean = getDouble(pString);
            sprintf(pString, "PDF half base width (assumed isosceles) = ");
            sstdev = getDouble(pString);
         }
         else if (gtype == 4)
         {
            gtype = PSUADE_PDF_BETA;
            sprintf(pString, "PDF alpha = ");
            smean = getDouble(pString);
            sprintf(pString, "PDF beta = ");
            sstdev = getDouble(pString);
         }
         else if (gtype == 5)
         {
            gtype = PSUADE_PDF_WEIBULL;
            sprintf(pString, "PDF lambda = ");
            smean = getDouble(pString);
            sprintf(pString, "PDF k = ");
            sstdev = getDouble(pString);
         }
         else if (gtype == 6)
         {
            gtype = PSUADE_PDF_GAMMA;
            sprintf(pString, "PDF alpha = ");
            smean = getDouble(pString);
            sprintf(pString, "PDF beta = ");
            sstdev = getDouble(pString);
         }
         else if (gtype == 7)
         {
            gtype = PSUADE_PDF_EXPONENTIAL;
            sprintf(pString, "PDF lambda = ");
            smean = getDouble(pString);
            sstdev = 0.0;
         }
         else if (gtype == 8)
         {
            gtype = PSUADE_PDF_F;
            sprintf(pString, "PDF d1 = ");
            smean = getDouble(pString);
            sprintf(pString, "PDF d2 = ");
            sstdev = getDouble(pString);
         }
         corMat.setDim(1,1);
         corMat.setEntry(0,0, 1.0e0);
         pdfman = new PDFManager();
         pdfman->initialize(1, &gtype, &smean, &sstdev, corMat, NULL, NULL);
         vecIn.setLength(iOne);
         vecOut.setLength(iOne);
         sprintf(pString,"Enter parameter value to fetch cumulative probability: ");
         ddata = getDouble(pString);
         vecIn.load(1, &ddata);
         slbound = 0;
         subound = 1;
         vecUpper.load(1, &subound);
         vecLower.load(1, &slbound);
         pdfman->getCDF(1, vecIn, vecOut, vecLower, vecUpper);
         printf("Cumulative probability = %e\n", vecOut[0]);
         delete pdfman;
         pdfman = NULL;
      }

      // +++ output_file 
      else if (!strcmp(command, "output_file"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("output_file: changes the default output filename.\n");
            printf("syntax: output_file <filename>\n");
            continue;
         }
         psOutputFilename_ = PSUADE_strdup(winput);
         if (psInputFilename_ == psOutputFilename_)
         {
            printf("WARNING: The output filename is the same as the input\n");
            printf("         filename. If you save the file it will overwrite\n");
            printf("         the input file.\n");
         }
      }

      // +++ setupguide 
      else if (!strcmp(command, "setupguide"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("setupGuide: show a short guide to help set up an application\n");
            continue;
         }
         setupGuide();
      }

      // +++ genworkflow 
      else if (!strcmp(command, "genworkflow"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("genworkflow: generate a UQ workfow for an application\n");
            continue;
         }
         printf("genworkflow not implemented yet.\n");
      }

      // +++ geninputfile 
      else if (!strcmp(command, "geninputfile"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("geninputfile: create a PSUADE input file\n");
            printf("syntax: geninputfile (no argument needed)\n");
            continue;
         }
         status = genSetup(0, dataFile);
         if (sampleInputs_  != NULL) delete [] sampleInputs_;
         if (sampleOutputs_ != NULL) delete [] sampleOutputs_;
         if (iLowerB_       != NULL) delete [] iLowerB_;
         if (iUpperB_       != NULL) delete [] iUpperB_;
         if (sampleStates_  != NULL) delete [] sampleStates_;
         if (inputNames_ != NULL)
         {
            for (ii = 0; ii < nInputs_; ii++)
               delete [] inputNames_[ii];
            delete [] inputNames_;
         }
         if (outputNames_ != NULL)
         {
            for (ii = 0; ii < nOutputs_; ii++)
               delete [] outputNames_[ii];
            delete [] outputNames_;
         }
         sampleInputs_ = sampleOutputs_ = NULL;
         sampleStates_ = NULL;
         iLowerB_ = iUpperB_ = NULL;
         inputNames_ = outputNames_ = NULL;

         if (status == 0)
         {
            printf("PSUADE can create the sample input data file for you.\n");
            sprintf(pString,"Create the file ? (y or n) ");
            getString(pString, winput);
            if (winput[0] == 'y')
            {
               getInputFromFile(dataFile);
               cleanUp();
               psuadeIO_ = new PsuadeData();
               psuadeIO_->setOutputLevel(0);
               strcpy(winput, psInputFilename_);
               status = psuadeIO_->readPsuadeFile(winput);
               if (status != 0) 
               {
                  printf("ERROR: cannot read file %s\n",winput);
                  continue;
               }
               psuadeIO_->getParameter("input_ninputs", pPtr);
               nInputs_ = pPtr.intData_;
               pINames.clean();
               psuadeIO_->getParameter("input_names", pINames);
               names = pINames.strArray_;
               pLower.clean();
               psuadeIO_->getParameter("input_lbounds", pLower);
               iLowerB_ = pLower.dbleArray_;
               pLower.dbleArray_ = NULL;
               pUpper.clean();
               psuadeIO_->getParameter("input_ubounds", pUpper);
               iUpperB_ = pUpper.dbleArray_;
               pUpper.dbleArray_ = NULL;
               inputNames_ = new char*[nInputs_+1];
               for (ii = 0; ii < nInputs_; ii++)
               {
                  inputNames_[ii] = new char[200]; 
                  strcpy(inputNames_[ii], names[ii]);
               }
               psuadeIO_->getParameter("output_noutputs", pPtr);
               nOutputs_ = pPtr.intData_;
               pONames.clean();
               psuadeIO_->getParameter("output_names", pONames);
               names = pONames.strArray_;
               outputNames_ = new char*[nOutputs_+1];
               for (ii = 0; ii < nOutputs_; ii++)
               {
                  outputNames_[ii] = new char[200]; 
                  strcpy(outputNames_[ii], names[ii]);
               }
               psuadeIO_->getParameter("method_sampling", pPtr);
               samplingMethod = pPtr.intData_;
               psuadeIO_->getParameter("method_nsamples", pPtr);
               nSamples_ = pPtr.intData_;
               psuadeIO_->getParameter("method_nreplications",pPtr);
               nReps = pPtr.intData_;
               psuadeIO_->getParameter("input_sample", pPtr);
               sampleInputs_  = pPtr.dbleArray_;
               psuadeIO_->getParameter("output_sample", pPtr);
               sampleOutputs_  = pPtr.dbleArray_;
               psuadeIO_->getParameter("output_states", pPtr);
               sampleStates_  = pPtr.intArray_;
               pPtr.intArray_ = NULL;
               pPtr.dbleArray_ = NULL;
               pINames.clean();
               pONames.clean();
               printf("==================================================\n");
               printf("The sample matrix is now stored in %s file.\n", 
                      psOutputFilename_);
               printf("You can also use genmars command to convert the\n");
               printf("sample matrix to a row-column format.\n");
               printf("==================================================\n");
            }
         }
      }

      // +++ genbatchfile 
      else if (!strcmp(command, "genbatchfile"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("genbatchfile: create a batch file for ensemble runs\n");
            printf("syntax: genbatchfile (no argument needed)\n");
            continue;
         }
         genBatchFile(0);
      }

      // +++ gendriver 
      else if (!strcmp(command, "gendriver"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("gendriver: create a driver file for runing PSUADE\n");
            printf("syntax: gendriver (no argument needed)\n");
            continue;
         }
         genDriver(0);
      }

      // +++ gendist 
      else if (!strcmp(command, "gendist"))
      {
         int    ns, gtype;
         double *sData, slbound, subound, smean, sstdev;

         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("gendist: create a sample using selected PDFs\n");
            printf("syntax: gendist (no argument needed)\n");
            continue;
         }
         sprintf(pString, "nSamples = ");
         ns = getInt(10, 10000000, pString);
         sprintf(pString, 
           "PDF type = (1) N (2) L (3) T (4) Beta (5) Exp (6) Weibull: ");
         gtype = getInt(1, 6, pString);
         if      (gtype == 1) gtype = PSUADE_PDF_NORMAL;
         else if (gtype == 2) gtype = PSUADE_PDF_LOGNORMAL;
         else if (gtype == 3) gtype = PSUADE_PDF_TRIANGLE;
         else if (gtype == 4) gtype = PSUADE_PDF_BETA;
         else if (gtype == 5) gtype = PSUADE_PDF_EXPONENTIAL;
         else if (gtype == 6) gtype = PSUADE_PDF_WEIBULL;
         sprintf(pString, "PDF parameter 1 (e.g. mean for N) = ");
         smean = getDouble(pString);
         if (gtype != PSUADE_PDF_EXPONENTIAL) 
         {
            sprintf(pString, "PDF parameter 2 (e.g. std dev for N) = ");
            sstdev = getDouble(pString);
         }
         else sstdev = 0.0;
         corMat.setDim(1,1);
         corMat.setEntry(0,0, 1.0e0);
         pdfman = new PDFManager();
         pdfman->initialize(1, &gtype, &smean, &sstdev, corMat, NULL, NULL);
         vecOut.setLength(ns);
         subound =  PSUADE_UNDEFINED;
         slbound = -PSUADE_UNDEFINED;
         vecUpper.load(1, &subound);
         vecLower.load(1, &slbound);
         pdfman->genSample(ns, vecOut, vecLower, vecUpper);
         delete pdfman;
         pdfman = NULL;

         sData = vecOut.getDVector();
         fp = fopen("sample1D", "w");
         if (fp == NULL)
         {
            printf("ERROR: cannot write to file sample1D.\n");
            continue;
         }
         fprintf(fp, "%d\n", ns);
         for (ii = 0; ii < ns; ii++) fprintf(fp, "%e\n", sData[ii]);
         fclose(fp);
         printf("data file created in sample1D.\n");
      }

      // +++ genexample 
      else if (!strcmp(command, "genexample"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("genexample: create a PSUADE example (driver+input file)\n");
            printf("syntax: genexample (no argument needed)\n");
            continue;
         }
         genDriver(1);
         genSetup(1, winput);
         printf("Now use: psuade psuade.in to run the example.\n");
      }

      // +++ genconfigfile 
      else if (!strcmp(command, "genconfigfile"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("genconfigfile: create a PSUADE config template file\n");
            printf("syntax: genconfigure (no argument needed)\n");
            printf("Configure files are used to modify settings.\n");
            printf("Configure files may be used to replace interactive\n");
            printf("queries from PSUADE.\n");
            continue;
         }
         sprintf(pString,"Enter the name of the configure file to write to: ");
         getString(pString, dataFile);
         dataFile[strlen(dataFile)-1] = '\0';
         if (!strcmp(dataFile, "\0")) 
         {
            printf("ERROR: invalid file name.\n");
            continue;
         }
         status = genConfigFileTemplate(dataFile);
         if (status == 0) 
            printf("genconfigfile completed - file name is %s.\n", dataFile);
         else 
            printf("ERROR: cannot write to file %s or no filename given.\n",
                   dataFile);
      }

      // +++ chkjobs 
      else if (!strcmp(command, "chkjobs"))
      {
         int  nPatterns, choice, nJobs1, nJobs2;
         char **patterns, checkFile[200];

         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("chkjobs: monitor the status of ensemble runs \n");
            printf("syntax: chkjobs (no argument needed)\n");
            printf("Use this command with '-h' to see more details.\n");
            continue;
         }
         printf("PSUADE will check for job status in the current directory.\n");
         printf("So first make sure you are in the right directory.\n");
         printf("You should have subdirectories each of which is a job.\n");
         printf("The subdirectory names should be something like: \n");
         printf("<dir_prefix>xxx where xxx is a number from 1 to njobs.\n");
         sprintf(pString,"Enter the subdirectory prefix now : ");
         getString(pString, winput);
         sscanf(winput, "%s", dirName);
         sprintf(pString,"Enter another level of subdirectory, if any : (or NONE) ");
         getString(pString, winput);
         sscanf(winput, "%s", subdirName);
         printf("What is considered to be failed runs?\n");
         printf("(1) that a certain file does not exist\n");
         printf("(2) both (1) and that a pattern does not exist in a file\n");
         printf("(3) that certain pattern(s)s exist in a given file\n");
         sprintf(pString,"Make your selection : ");
         choice = getInt(1, 3, pString);
         nPatterns = 0;
         if (choice == 1)
         {
            sprintf(pString,"Enter the name of file that has to exist : ");
            getString(pString, winput);
            sscanf(winput, "%s", checkFile);
         }
         else
         {
            sprintf(pString,"Enter the name of file to check patterns: ");
            getString(pString, winput);
            sscanf(winput, "%s", checkFile);
            if (choice == 2) nPatterns = 1;
            else
            {
               sprintf(pString,"How many patterns (1 - 10)? ");
               nPatterns = getInt(1, 10, pString);
            }
            patterns = new char*[nPatterns];
            for (ii = 0; ii < nPatterns; ii++)
            {
               patterns[ii] = new char[510]; 
               printf("Enter pattern %d : ",ii+1);
               fgets(patterns[ii], 500, stdin); 
            }
         }
         sprintf(pString,"Enter the first job number to be probed (1 - ?): ");
         nJobs1 = getInt(1, 1000000, pString);
         sprintf(pString,"Enter the last  job number to be probed (%d - ?): ",
                 nJobs1+1);
         nJobs2 = getInt(nJobs1+1, 1000000, pString);
         fErr = fopen("relaunchJobs.py", "w");
         if (fErr == NULL)
            printf("ERROR: cannot open file to store job status info.\n");

         if (fErr != NULL)
         {
            fprintf(fErr, "import os\n");
            fprintf(fErr, "import sys\n\n");
            fprintf(fErr, "jobs = [");
         }
         count = 0;
         for (ii = nJobs1; ii <= nJobs2; ii++)
         {
            if (outputLevel_ > 0)
               printf("Processing job %d\n", ii);
            if (choice == 1)
            {
               if (strncmp(subdirName, "NONE", 4) == 0)
                  sprintf(winput, "%s%d/%s", dirName, ii, checkFile);
               else
                  sprintf(winput,"%s%d/%s/%s",dirName,ii,subdirName,checkFile);
               fp = fopen(winput, "r");
               if (fp == NULL) 
               {
                  if (fErr != NULL)
                  {
                     if (count == 0) fprintf(fErr, "%d", ii);
                     else            fprintf(fErr, ",%d", ii);
                  }
                  else
                  {
                     printf("%s%d/%s fails (1): file does not exist.\n",
                            dirName, ii, checkFile);
                  } 
                  count++;
               }
               else fclose(fp);
            }
            if (choice == 2)
            {
               if (strncmp(subdirName, "NONE", 4) == 0)
                  sprintf(winput, "%s%d/%s", dirName, ii, checkFile);
               else
                  sprintf(winput, "%s%d/%s/%s",dirName,ii,subdirName,checkFile);
               fp = fopen(winput, "r");
               if (fp == NULL) 
               {
                  if (fErr != NULL)
                  {
                     if (count == 0) fprintf(fErr, "%d", ii);
                     else            fprintf(fErr, ",%d", ii);
                  }
                  else
                  {
                     printf("%s%d fails (2a): file does not exist.\n",
                            dirName, ii);
                  }
                  count++;
               }
               else 
               {
                  fclose(fp);
                  sprintf(command, "grep \"%s\" %s > /dev/null", 
                          patterns[0], winput);
                  status = system(command);
                  if (status != 0)
                  {
                     if (fErr != NULL)
                     {
                        if (count == 0) fprintf(fErr, "%d", ii);
                        else            fprintf(fErr, ",%d", ii);
                     }
                     else
                     {
                        printf("%s%d fails (2b): file does not exist.\n",
                               dirName, ii);
                     }
                     count++;
                  }
               }
            }
            if (choice == 3)
            {
               if (strncmp(subdirName, "NONE", 4) == 0)
                  sprintf(winput, "%s%d/%s", dirName, ii, checkFile);
               else
                  sprintf(winput, "%s%d/%s/%s", dirName, ii, subdirName, 
                          checkFile);
               fp = fopen(winput, "r");
               if (fp == NULL) 
               {
                  if (fErr != NULL)
                  {
                     if (count == 0) fprintf(fErr, "%d", ii);
                     else            fprintf(fErr, ",%d", ii);
                  }
                  else
                  {
                     printf("%s%d fails (3a): file does not exist.\n",
                            dirName, ii);
                  }
                  count++;
               }
               else 
               {
                  fclose(fp);
                  for (jj = 0; jj < nPatterns; jj++)
                  {
                     sprintf(command, "grep \"%s\" %s > /dev/null", 
                             patterns[jj], winput);
                     status = system(command);
                     if (status == 0)
                     {
                        if (fErr != NULL)
                        {
                           if (count == 0) fprintf(fErr, "%d", ii);
                           else            fprintf(fErr, ",%d", ii);
                        }
                        else
                        {
                           printf("%s%d fails (3b): file does not exist.\n",
                                  dirName, ii);
                        }
                        count++;
                        break;
                     }
                  }
               }
            }
         }
         if (fErr != NULL)
         {
            fprintf(fErr,"]");
            fprintf(fErr,"\n");
            fprintf(fErr,"for index in jobs:\n");
            fprintf(fErr,"#  insert your clean-up procedure\n");
            fprintf(fErr,"## e.g. cmd = \"/bin/rm -r workdir.\" + str(index)\n");
            fprintf(fErr,"##      os.system(cmd)\n\n");
            fprintf(fErr,"#  insert your re-start procedure\n");
            fprintf(fErr,"## e.g. cmd = \"driver.py psuadeApps_ct.in.\" + str(index)\n");
            fprintf(fErr,"psuadeApps_ct.out.\" + str(index)\n");
            fprintf(fErr,"##      os.system(cmd)\n\n");
            fprintf(fErr,"#  insert your job submission procedure\n");
            fprintf(fErr,"## cmd = \"cd %s.\" + str(index) + \"; \"\n",
                    dirName);
            fprintf(fErr,"## cmd = cmd + \"/usr/bin/psub batchFile.\" +str(index)\n");
            fprintf(fErr,"## os.system(cmd)\n\n");
            fclose(fErr);
            printf("The failed jobs are in the file relaunchJobs.py.\n");
         }
         if (nPatterns > 0)
         {
            for (ii = 0; ii < nPatterns; ii++) delete [] patterns[ii];
            delete [] patterns;
         }
      }

      // +++ list1 
      else if (!strcmp(command, "list1"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("list1: list all points with 1 selected input and 1 output\n");
            printf("syntax: list1 (no argument needed)\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data (load data first).\n");
         }
         else
         {
            sprintf(pString, "Enter input number (1 - %d) : ", nInputs_);
            iInd = getInt(1, nInputs_, pString);
            iInd--;
            if (nOutputs_ == 1) outputID = 0;
            else
            {
               sprintf(pString,"Enter output number (1 - %d) : ", nOutputs_);
               outputID = getInt(1, nOutputs_, pString);
               outputID--;
            }
            sprintf(pString,"Sort the input (y or n) ? ");
            getString(pString, winput);
            if (winput[0] == 'n') 
            {
               sprintf(pString,"Sort the output (y or n) ? ");
               getString(pString, winput);
               if (winput[0] == 'n') winput[0] = 'N'; 
               if (winput[0] == 'y') winput[0] = 'Y'; 
            }
            tempX = new double[nSamples_];
            tempY = new double[nSamples_];
            tempW = new double[nSamples_];
            for (sInd = 0; sInd < nSamples_; sInd++)
            {
               tempX[sInd] = sampleInputs_[sInd*nInputs_+iInd];
               tempY[sInd] = sampleOutputs_[sInd*nOutputs_+outputID];
               tempW[sInd] = (double) sInd + 1;
            }
            if (winput[0] == 'y') 
               sortDbleList3(nSamples_, tempX, tempY, tempW);
            else if (winput[0] == 'Y') 
               sortDbleList3(nSamples_, tempY, tempX, tempW);
            for (sInd = 0; sInd < nSamples_; sInd++)
               printf("%6d: Sample %7d : input = %16.8e, output = %16.8e\n",
                      sInd+1, (int) tempW[sInd], tempX[sInd], tempY[sInd]);
            delete [] tempX;
            delete [] tempY;
            delete [] tempW;
            tempX = tempY = tempW = NULL;
         }
      }

      // +++ list2 
      else if (!strcmp(command, "list2"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("list2: list all points with 2 selected inputs and 1 output\n");
            printf("syntax: list2 (no argument needed)\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data (load data first).\n");
         }
         else
         {
            sprintf(pString, "Enter input number 1 (1 - %d) : ", nInputs_);
            iInd1 = getInt(1, nInputs_, pString);
            iInd1--;
            sprintf(pString, "Enter input number 2 (1 - %d) : ", nInputs_);
            iInd2 = getInt(1, nInputs_, pString);
            iInd2--;
            if (nOutputs_ == 1) outputID = 0;
            else
            {
               sprintf(pString,"Enter output number (1 - %d) : ", nOutputs_);
               outputID = getInt(1, nOutputs_, pString);
               outputID--;
            }
            sprintf(pString,"Sort the inputs (y or n) ? ");
            getString(pString, winput);
            if (winput[0] == 'n') 
            {
               sprintf(pString,"Sort the output (y or n) ? ");
               getString(pString, winput);
               if (winput[0] == 'n') winput[0] = 'N'; 
               if (winput[0] == 'y') winput[0] = 'Y'; 
            }
            tempX = new double[nSamples_];
            tempY = new double[nSamples_];
            tempW = new double[nSamples_];
            tempV = new double[nSamples_];
            for (sInd = 0; sInd < nSamples_; sInd++)
            {
               tempX[sInd] = sampleInputs_[sInd*nInputs_+iInd1];
               tempV[sInd] = sampleInputs_[sInd*nInputs_+iInd2];
               tempY[sInd] = sampleOutputs_[sInd*nOutputs_+outputID];
               tempW[sInd] = (double) sInd;
            }
            if (winput[0] == 'y') 
            {
               sortDbleList4(nSamples_, tempX, tempV, tempY, tempW);
               sInd = 0;
               count = 1; 
               while (sInd < nSamples_) 
               {
                  sInd++;
                  if (tempX[sInd] == tempX[sInd-1]) count++;
                  else
                  {
                     if (count > 1)
                     {
                        sortDbleList4(count, &(tempV[sInd-count]), 
                              &(tempX[sInd-count]), &(tempW[sInd-count]),
                              &(tempY[sInd-count]));
                     }
                     count = 1;
                  }
               }
               if (count > 1)
               {
                  sortDbleList4(count, &(tempV[sInd-count]),
                         &(tempX[sInd-count]), &(tempW[sInd-count]),
                         &(tempY[sInd-count]));
               }
            }
            else if (winput[0] == 'Y') 
            {
               sortDbleList4(nSamples_, tempY, tempX, tempV, tempW);
            }
            for (sInd = 0; sInd < nSamples_; sInd++)
            {
               printf("%6d: Sample %7d : ", sInd+1, ((int) tempW[sInd])+1);
               printf("inputs = (%12.4e, %12.4e), output = %12.4e\n",
                      tempX[sInd], tempV[sInd], tempY[sInd]);
            }
            delete [] tempX;
            delete [] tempY;
            delete [] tempW;
            delete [] tempV;
            tempX = tempY = tempW = tempV = NULL;
         }
      }

      // +++ listall 
      else if (!strcmp(command, "listall"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("listall: list all points with all inputs and 1 output\n");
            printf("syntax: listall (no argument needed)\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data (load data first).\n");
         }
         else
         {
            if (nOutputs_ == 1) outputID = 0;
            else
            {
               sprintf(pString,"Enter output number (1 - %d) : ", nOutputs_);
               outputID = getInt(1, nOutputs_, pString);
               outputID--;
            }
            sprintf(pString,"Sort the output (y or n) ? ");
            getString(pString, winput);
            if (winput[0] == 'y')
            {
               int    *tmpInds = new int[nSamples_];
               double *tmpOuts = new double[nSamples_];
               for (sInd = 0; sInd < nSamples_; sInd++)
               {
                  tmpInds[sInd] = sInd;
                  tmpOuts[sInd] = sampleOutputs_[sInd*nOutputs_+outputID];
               }
               sortDbleList2a(nSamples_, tmpOuts, tmpInds);
               for (sInd = 0; sInd < nSamples_; sInd++)
               {
                  kk = tmpInds[sInd];
                  printf("%6d: Sample %7d : ", sInd+1, kk+1);
                  for (ii = 0; ii < nInputs_; ii++)
                     printf("%12.4e ", sampleInputs_[kk*nInputs_+ii]);
                  printf("= %12.4e\n", tmpOuts[ss]);
               }
               delete [] tmpOuts;
               delete [] tmpInds;
            }
            else
            {
               for (ss = 0; ss < nSamples_; ss++)
               {
                  printf("Sample %7d : ", ss+1);
                  for (ii = 0; ii < nInputs_; ii++)
                     printf("%12.4e ", sampleInputs_[ss*nInputs_+ii]);
                  printf("= %12.4e\n", sampleOutputs_[ss]);
               }
            }
         }
      }

      // +++ disp_sample or sshow
      else if (!strcmp(command, "disp_sample") || !strcmp(command, "sshow"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("disp_sample (or sshow): display one sample point\n");
            printf("syntax: disp_sample or sshow (no argument needed)\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data (load data first).\n");
         }
         else
         {
            sprintf(pString, "Enter sample number (1 - %d) : ", nSamples_);
            ind = getInt(1, nSamples_, pString);
            ind--;
            printf("Sample %7d : \n", ind+1);
            for (ii = 0; ii < nInputs_; ii++)
               printf("   input  %3d = %16.8e\n", ii+1, 
                      sampleInputs_[ind*nInputs_+ii]);
            for (ii = 0; ii < nOutputs_; ii++)
               printf("   output %3d = %16.8e\n", sInd+1, 
                      sampleOutputs_[ind*nOutputs_+ii]);
         }
      }

      // +++ max 
      else if (!strcmp(command, "max"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("max: find the sample point with maximum output value\n");
            printf("syntax: max (no argument needed)\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data (load data first).\n");
         }
         else
         {
            if (nOutputs_ == 1) outputID = 0;
            else
            {
               sprintf(pString,"Enter output number (1 - %d) : ", nOutputs_);
               outputID = getInt(1, nOutputs_, pString);
               outputID--;
            }
            double maxY = sampleOutputs_[outputID];
            int    maxI = 0;
            for (sInd = 1; sInd < nSamples_; sInd++)
            {
               if (sampleOutputs_[sInd*nOutputs_+outputID] > maxY)
               {
                  maxI  = sInd;
                  maxY = sampleOutputs_[sInd*nOutputs_+outputID];
               }
            }
            printf("Sample %d gives maximum for output %d\n",
                    maxI+1,outputID+1);
            for (iInd = 0; iInd < nInputs_; iInd++)
               printf("  input %4d = %16.8e\n", iInd+1,
                      sampleInputs_[maxI*nInputs_+iInd]);
            printf("  ====> output %4d = %16.8e\n", outputID+1,
                   sampleOutputs_[maxI*nOutputs_+outputID]);
         }
      }

      // +++ min 
      else if (!strcmp(command, "min"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("min: find the sample point with minimum output value\n");
            printf("syntax: min (no argument needed)\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data (load data first).\n");
         }
         else
         {
            if (nOutputs_ == 1) outputID = 0;
            else
            {
               sprintf(pString,"Enter output number (1 - %d) : ", nOutputs_);
               outputID = getInt(1, nOutputs_, pString);
               outputID--;
            }
            double minY = sampleOutputs_[outputID];
            int    minI = 0;
            for (sInd = 1; sInd < nSamples_; sInd++)
            {
               if (sampleOutputs_[sInd*nOutputs_+outputID] < minY)
               {
                  minI = sInd;
                  minY = sampleOutputs_[sInd*nOutputs_+outputID];
               }
            }
            printf("Sample %d gives minimum for output %d\n",
                    minI+1,outputID+1);
            for (iInd = 0; iInd < nInputs_; iInd++)
               printf("  input  %4d = %16.8e\n", iInd+1,
                      sampleInputs_[minI*nInputs_+iInd]);
            printf("  ====> output %4d = %16.8e\n", outputID+1,
                   sampleOutputs_[minI*nOutputs_+outputID]);
         }
      }

      // +++ idelete 
      else if (!strcmp(command, "idelete"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("idelete: delete a sample input \n");
            printf("syntax: idelete (no argument needed)\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data (load data first).\n");
            continue;
         }
         if (nInputs_ == 1)
         {
            printf("You have only one input left -> no deletion.\n");
            continue;
         }
         for (ii = 0; ii < nInputs_; ii++)
            printf("Input %3d = %s\n", ii+1, inputNames_[ii]);
            
         int *indexSet = new int[nInputs_];
         int inpCnt    = 0;
         for (ii = 0; ii < nInputs_; ii++) indexSet[ii] = 1;
         sprintf(pString,"How many inputs to remove? (1-%d) ",nInputs_-1);
         inpCnt = getInt(1, nInputs_-1, pString);
         sprintf(pString,"Enter input number (1 - %d) : ", nInputs_);
         for (ii = 0; ii < inpCnt; ii++)
         {
            iInd = getInt(1, nInputs_, pString);
            printf("You are removing input %d (%s)\n",iInd,inputNames_[iInd-1]);
            indexSet[iInd-1] = 0;
         }
         inpCnt = 0;
         for (ii = 0; ii < nInputs_; ii++)
         {
            if (indexSet[ii] == 1)
            {
               iLowerB_[inpCnt] = iLowerB_[ii];
               iUpperB_[inpCnt] = iUpperB_[ii];
               inpCnt++;
            }
         }
         for (sInd = 0; sInd < nSamples_; sInd++)
         {
            kk = 0;
            for (ii = 0; ii < nInputs_; ii++)
            {
               if (indexSet[ii] == 1)
               {
                  sampleInputs_[sInd*inpCnt+kk] = 
                        sampleInputs_[sInd*nInputs_+ii];
                  kk++;
               }
            }
         }
         kk = 0;
         for (ii = 0; ii < nInputs_; ii++)
         {
            if (indexSet[ii] == 1)
            {
               strcpy(inputNames_[kk], inputNames_[ii]);
               kk++;
            }
         }
         kk = 0;
         for (ii = 0; ii < nInputs_; ii++)
         {
            if (indexSet[ii] == 1)
            {
               inputPDFs_[kk]  = inputPDFs_[ii];
               inputMeans_[kk] = inputMeans_[ii];
               inputStds_[kk]  = inputStds_[ii];
               for (jj = 0; jj < nInputs_; jj++)
               {
                  ddata = inputCMat_->getEntry(ii,jj);
                  inputCMat_->setEntry(kk,jj,ddata);
               }
               for (jj = 0; jj < nInputs_; jj++)
               {
                  ddata = inputCMat_->getEntry(jj,ii);
                  inputCMat_->setEntry(jj,kk,ddata);
               }
               kk++;
            }
         }
         Matrix *tmpMat = new Matrix();
         tmpMat->setDim(inpCnt,inpCnt);
         for (ii = 0; ii < inpCnt; ii++)
         {
            for (jj = 0; jj < inpCnt; jj++)
            {
               ddata = inputCMat_->getEntry(ii,jj);
               tmpMat->setEntry(ii,jj,ddata);
            }
         }
         delete [] indexSet;
         delete inputCMat_;
         inputCMat_ = tmpMat;
         nInputs_ = inpCnt;
         psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,iLowerB_,
                                       iUpperB_,sampleInputs_,inputNames_,
                                       inputPDFs_,inputMeans_,inputStds_,
                                       inputCMat_); 
         if (currSession != NULL) delete currSession;
         currSession = new PsuadeSession();
         psuadeIO_->getSession(currSession);
         printf("idelete completed. Use 'write' to store.\n");
      }

      // +++ odelete 
      else if (!strcmp(command, "odelete"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("odelete: delete a sample output \n");
            printf("syntax: odelete (no argument needed)\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data (load data first).\n");
            continue;
         }
         if (nOutputs_ == 1)
         {
            printf("You have only one output left -> no deletion.\n");
            continue;
         }
         if (outputNames_ != NULL) 
         {
            for (ii = 0; ii < nOutputs_; ii++)
               printf("Output %3d = %s\n", ii+1, outputNames_[ii]);
         }
         sprintf(pString,"Enter output number (1 - %d) : ", nOutputs_);
         outputID = getInt(1, nOutputs_, pString);
         outputID--;
         if (outputNames_ != NULL) 
         {
            for (ii = outputID+1; ii < nOutputs_; ii++)
               strcpy(outputNames_[ii-1], outputNames_[ii]);
            delete [] outputNames_[nOutputs_-1];
            outputNames_[nOutputs_-1] = NULL;
         }
         for (ii = 0; ii < nSamples_; ii++)
         {
            for (jj = 0; jj < outputID; jj++)
                  sampleOutputs_[ii*(nOutputs_-1)+jj] = 
                                sampleOutputs_[ii*nOutputs_+jj];
            for (jj = outputID+1; jj < nOutputs_; jj++)
               sampleOutputs_[ii*(nOutputs_-1)+jj-1] = 
                             sampleOutputs_[ii*nOutputs_+jj];
         }
         nOutputs_--;
         psuadeIO_->updateOutputSection(nSamples_,nOutputs_,
                         sampleOutputs_,sampleStates_,outputNames_); 
         if (currSession != NULL) delete currSession;
         currSession = new PsuadeSession();
         psuadeIO_->getSession(currSession);
         printf("odelete completed. Use 'write' to store.\n");
      }

      // +++ sdelete 
      else if (!strcmp(command, "sdelete"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("sdelete: delete a sample point \n");
            printf("syntax: sdelete (no argument needed)\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data (load data first).\n");
            continue;
         }
         if (nSamples_ == 1)
         {
            printf("You have only one sample point left -> no deletion.\n");
            continue;
         }
         sprintf(pString,"Enter sample point to delete (1 - %d) : ",nSamples_);
         ss = getInt(1, nSamples_, pString);
         ss--;
         for (ii = ss; ii < nSamples_; ii++)
         {
            for (jj = 0; jj < nInputs_; jj++)
               sampleInputs_[ii*nInputs_+jj] = 
                             sampleInputs_[(ii+1)*nInputs_+jj];
            for (jj = 0; jj < nOutputs_; jj++)
               sampleOutputs_[ii*nOutputs_+jj] = 
                             sampleOutputs_[(ii+1)*nOutputs_+jj];
         }
         nSamples_--;
         psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,
                                       NULL,NULL,sampleInputs_,NULL, 
                                       NULL, NULL, NULL, NULL);
         psuadeIO_->updateOutputSection(nSamples_,nOutputs_,
                       sampleOutputs_,sampleStates_,outputNames_); 
         psuadeIO_->updateMethodSection(-1,nSamples_,-1,-1,-1);
         if (currSession != NULL) delete currSession;
         currSession = new PsuadeSession();
         psuadeIO_->getSession(currSession);
         printf("sdelete completed. Use 'write' to store.\n");
      }

      // +++ ishuffle 
      else if (!strcmp(command, "ishuffle"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("ishuffle: re-arrange the orderings of input parameters.\n");
            printf("syntax: ishuffle (no argument needed)\n");
            printf("NOTE: this command does not change the data file.\n");
            printf("      until a 'write' is issued.\n");
            continue;
         }
         if (psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data (load data first).\n");
            continue;
         }
         int *indexSet = new int[nInputs_];
         int indexCnt  = 0;
         while (indexCnt < nInputs_)
         {
            sprintf(pString, "Enter the %d-th input (1 - %d) : ",ii+1,nInputs_);
            indexSet[indexCnt] = getInt(1, nInputs_, pString);
            indexCnt++;
         }
         double *samInps   = new double[nSamples_*nInputs_];
         char   **inpNames = new char*[nInputs_];
         for (ii = 0; ii < nInputs_; ii++)
         {
            kk = indexSet[ii] - 1;
            for (jj = 0; jj < nSamples_; jj++)
               samInps[jj*nInputs_+ii] = sampleInputs_[jj*nInputs_+kk];
            if (inputNames_ != NULL) inpNames[ii] = inputNames_[kk];
         }
         delete [] sampleInputs_;
         sampleInputs_ = samInps;
         if (inputNames_ != NULL) 
         {
            for (ii = 0; ii < nInputs_; ii++) 
            {
               inputNames_[ii] = inpNames[ii];
               inpNames[ii] = NULL;
            }
            delete [] inpNames;
         }
         double *tmpLBs = new double[nInputs_];
         double *tmpUBs = new double[nInputs_];
         for (ii = 0; ii < nInputs_; ii++)
         {
            kk = indexSet[ii]-1;
            tmpLBs[ii] = iLowerB_[kk];
            tmpUBs[ii] = iUpperB_[kk];
         }
         delete [] iLowerB_;
         delete [] iUpperB_;
         iLowerB_ = tmpLBs;
         iUpperB_ = tmpUBs;
         int    *tmpPDFs  = new int[nInputs_];
         double *tmpMeans = new double[nInputs_];
         double *tmpStds  = new double[nInputs_];
         tempW = new double[nInputs_];
         for (ii = 0; ii < nInputs_; ii++)
         {
            kk = indexSet[ii] - 1;
            tmpPDFs[ii]  = inputPDFs_[kk];
            tmpMeans[ii] = inputMeans_[kk];
            tmpStds[ii]  = inputStds_[kk];
         }
         delete [] inputPDFs_;
         delete [] inputMeans_;
         delete [] inputStds_;
         inputPDFs_ = tmpPDFs;
         inputMeans_ = tmpMeans;
         inputStds_ = tmpStds;
         Matrix *tmpMat = new Matrix();
         tmpMat->setDim(nInputs_, nInputs_);
         for (ii = 0; ii < nInputs_; ii++)
         {
            kk = indexSet[ii] - 1;
            for (jj = 0; jj < nInputs_; jj++)
            {
               ll = indexSet[jj]-1;
               ddata = inputCMat_->getEntry(kk,ll);
               tmpMat->setEntry(ii,jj,ddata);
            }
         }
         delete [] indexSet;
         delete inputCMat_;
         inputCMat_ = tmpMat;
         psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,iLowerB_,
                              iUpperB_,sampleInputs_,inputNames_,inputPDFs_,
                              inputMeans_,inputStds_,inputCMat_); 
         if (currSession != NULL) delete currSession;
         currSession = new PsuadeSession();
         psuadeIO_->getSession(currSession);
         printf("ishuffle completed. Use 'write' to store.\n");
      }

      // +++ iselect_index 
      else if (!strcmp(command, "iselect_index"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("iselect_index: select and re-order a subset of inputs.\n");
            printf("syntax: iselect_index (no argument needed)\n");
            printf("This command can be used in conjunction with iread,\n");
            printf("iadd and iwrite to convert a posterior sample (from\n");
            printf("rsmcmc) into a sample file that can be understood by\n");
            printf("the S-type PDF for use in multi-stage UQ analysis.\n");
            continue;
         }
         if (psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data (load data first).\n");
            continue;
         }
         int *indexSet = new int[nInputs_+1];
         int indexCnt  = 0;
         while (1)
         {
            sprintf(pString, "Enter the %d-th input (1 - %d, or 0 if done) : ",
                    indexCnt+1,nInputs_);
            indexSet[indexCnt] = getInt(0, nInputs_, pString);
            if (indexSet[indexCnt] == 0) break;
            indexCnt++;
            if (indexCnt > nInputs_)
            {
               printf("WARNING: Your last input should be '0'.\n");
               printf("         The last index is ignored.\n");
               indexCnt = nInputs_;
               break;
            }
         }
         if (indexCnt == 0)
         {
            printf("ERROR: you have not entered any inputs.\n");
            continue;
         }
         double *samInps = new double[nSamples_*indexCnt];
         char   **inpNames = new char*[indexCnt];
         for (ii = 0; ii < indexCnt; ii++)
         {
            kk = indexSet[ii] - 1;
            for (jj = 0; jj < nSamples_; jj++)
               samInps[jj*indexCnt+ii] = sampleInputs_[jj*nInputs_+kk];
            if (inputNames_ != NULL) inpNames[ii] = inputNames_[kk];
         }
         delete [] sampleInputs_;
         sampleInputs_ = samInps;
         for (ii = 0; ii < indexCnt; ii++) 
         {
            if (inputNames_ != NULL) inputNames_[ii] = inpNames[ii];
            inpNames[ii] = NULL;
         }
         delete [] inpNames;
         double *tmpLBs = new double[indexCnt];
         double *tmpUBs = new double[indexCnt];
         for (ii = 0; ii < indexCnt; ii++)
         {
            kk = indexSet[ii] - 1;
            tmpLBs[ii] = iLowerB_[kk];
            tmpUBs[ii] = iUpperB_[kk];
         }
         delete [] iLowerB_;
         delete [] iUpperB_;
         iLowerB_ = tmpLBs;
         iUpperB_ = tmpUBs;
         int    *tmpPDFs = new int[indexCnt];
         double *tmpMeans = new double[indexCnt];
         double *tmpStds  = new double[indexCnt];
         for (ii = 0; ii < indexCnt; ii++)
         {
            kk = indexSet[ii] - 1;
            tmpPDFs[ii]  = inputPDFs_[kk];
            tmpMeans[ii] = inputMeans_[kk];
            tmpStds[ii]  = inputStds_[kk];
         }
         delete [] inputPDFs_;
         delete [] inputMeans_;
         delete [] inputStds_;
         inputPDFs_  = tmpPDFs;
         inputMeans_ = tmpMeans;
         inputStds_  = tmpStds;
         Matrix *tmpMat = new Matrix();
         tmpMat->setDim(indexCnt, indexCnt);
         for (ii = 0; ii < indexCnt; ii++)
         {
            kk = indexSet[ii] - 1;
            for (jj = 0; jj < indexCnt; jj++)
            {
               ll = indexSet[jj]-1;
               ddata = inputCMat_->getEntry(kk,ll);
               tmpMat->setEntry(ii,jj,ddata);
            }
         }
         delete [] indexSet;
         delete inputCMat_;
         inputCMat_ = tmpMat;
         nInputs_ = indexCnt;
         psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,iLowerB_,
                              iUpperB_,sampleInputs_,inputNames_,inputPDFs_,
                              inputMeans_,inputStds_,inputCMat_); 
         samplingMethod = PSUADE_SAMP_MC;
         nReps = 1;
         psuadeIO_->updateMethodSection(samplingMethod,nSamples_,nReps,-1,-1);
         if (currSession != NULL) delete currSession;
         currSession = new PsuadeSession();
         psuadeIO_->getSession(currSession);
         printf("iselect_index completed. Use 'write' to store.\n");
      }

      // +++ iselect_name 
      else if (!strcmp(command, "iselect_name"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("iselect_name: select and re-order a subset of inputs.\n");
            printf("syntax: iselect_index (no argument needed)\n");
            printf("This command can be used in conjunction with iread, \n");
            printf("iadd and iwrite to convert a posterior sample (from\n");
            printf("rsmcmc) into a sample file that can be understood by\n");
            printf("the S-type PDF for multi-stage UQ analysis.\n");
            continue;
         }
         if (psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data (load data first).\n");
            continue;
         }
         if (inputNames_ == NULL)
         {
            printf("ERROR: input name list does not exist.\n");
            continue;
         }
         int *indexSet = new int[nInputs_];
         int indexCnt = 0;
         printf("Current input names:\n");
         for (ii = 0; ii < nInputs_; ii++)
            printf("Input %5d : name = %s\n",ii+1,inputNames_[ii]);
         while (1)
         {
            sprintf(pString,"Enter the %d-th input name (enter 'DONE' if done) : ",
                    indexCnt+1);
            getString(pString, cString);
            kk = strlen(cString);
            cString[kk-1] = '\0';
            if (!strcmp(cString, "DONE")) break;
            for (ii = 0; ii < nInputs_; ii++)
            {
               if (!strcmp(cString, inputNames_[ii])) 
               {
                  indexSet[indexCnt] = ii;
                  break;
               }
            }
            if (ii == nInputs_) printf("ERROR: input name %s not found.\n",cString);
            else               indexCnt++;
            if (indexCnt > nInputs_)
            {
               printf("WARNING: Your last input should be a 'DONE'.\n");
               printf("         The last name is ignored.\n");
               indexCnt = nInputs_;
               break;
            }
         }
         if (indexCnt == 0)
         {
            printf("ERROR: you have not entered any inputs.\n");
            continue;
         }
         double *samInps = new double[nSamples_*indexCnt];
         char   **inpNames  = new char*[indexCnt];
         for (ii = 0; ii < indexCnt; ii++)
         {
            kk = indexSet[ii];
            for (jj = 0; jj < nSamples_; jj++)
               samInps[jj*indexCnt+ii] = sampleInputs_[jj*nInputs_+kk];
            if (inputNames_ != NULL) inpNames[ii] = inputNames_[kk];
         }
         delete [] sampleInputs_;
         sampleInputs_ = samInps;
         for (ii = 0; ii < indexCnt; ii++) 
         {
            if (inputNames_ != NULL) inputNames_[ii] = inpNames[ii];
            inpNames[ii] = NULL;
         }
         delete [] inpNames;
         double *tmpLBs = new double[indexCnt];
         double *tmpUBs = new double[indexCnt];
         for (ii = 0; ii < indexCnt; ii++)
         {
            kk = indexSet[ii];
            tmpLBs[ii] = iLowerB_[kk];
            tmpUBs[ii] = iUpperB_[kk];
         }
         delete [] iLowerB_;
         delete [] iUpperB_;
         iLowerB_ = tmpLBs;
         iUpperB_ = tmpUBs;
         int    *tmpPDFs  = new int[indexCnt];
         double *tmpMeans = new double[indexCnt];
         double *tmpStds  = new double[indexCnt];
         for (ii = 0; ii < indexCnt; ii++)
         {
            kk = indexSet[ii];
            tmpPDFs[ii] = inputPDFs_[kk];
            tmpMeans[ii] = inputMeans_[kk];
            tmpStds[ii] = inputStds_[kk];
         }
         if (inputPDFs_  != NULL) delete [] inputPDFs_;
         if (inputMeans_ != NULL) delete [] inputMeans_;
         if (inputStds_  != NULL) delete [] inputStds_;
         inputPDFs_  = tmpPDFs;
         inputMeans_ = tmpMeans;
         inputStds_  = tmpStds;
         Matrix *tmpMat = new Matrix();
         tmpMat->setDim(indexCnt, indexCnt);
         for (ii = 0; ii < indexCnt; ii++)
         {
            kk = indexSet[ii];
            for (jj = 0; jj < indexCnt; jj++)
            {
               ll = indexSet[jj];
               ddata = inputCMat_->getEntry(kk,ll);
               tmpMat->setEntry(ii,jj,ddata);
            }
         }
         delete [] indexSet;
         delete inputCMat_;
         inputCMat_ = tmpMat;
         nInputs_ = indexCnt;
         psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,iLowerB_,
                              iUpperB_,sampleInputs_,inputNames_,inputPDFs_,
                              inputMeans_,inputStds_,inputCMat_); 
         samplingMethod = PSUADE_SAMP_MC;
         nReps = 1;
         psuadeIO_->updateMethodSection(samplingMethod,nSamples_,nReps,-1,-1);
         if (currSession != NULL) delete currSession;
         currSession = new PsuadeSession();
         psuadeIO_->getSession(currSession);
         printf("iselect_name completed. Use 'write' to store.\n");
      }

      // +++ itag 
      else if (!strcmp(command, "itag"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("itag: tag sample points based on input values\n");
            printf("syntax: itag (no argument needed)\n");
            printf("This command is useful for extracting indices of\n");
            printf("sample points with input values falling in a given\n");
            printf("range. This command can be called multiple times to\n");
            printf("impose multiple filters (AND operation) on the sample.\n");
            printf("This tagged list can be reset to empty by a 'load'.\n");
            printf("NOTE: this command does not change the data file.\n");
            printf("      Only display sample index information after\n");
            printf("      tagging.\n");
            continue;
         }
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no data (load data first).\n");
            continue;
         }
         if (tagArray_ == NULL && nSamples_ > 0) 
         {
            tagArray_ = new int[nSamples_];
            for (ii = 0; ii < nSamples_; ii++) tagArray_[ii] = 1;
         } 
         inputID = 0;
         if (nInputs_ > 1)
         {
            sprintf(pString, "Enter input number (1 - %d) : ", nInputs_);
            inputID = getInt(1, nInputs_, pString);
            inputID--;
         }
         Xmin = sampleInputs_[inputID];
         for (sInd = 1; sInd < nSamples_; sInd++)
            if (sampleInputs_[sInd*nInputs_+inputID] < Xmin)
               Xmin = sampleInputs_[sInd*nInputs_+inputID];
         Xmax = sampleInputs_[inputID];
         for (sInd = 1; sInd < nSamples_; sInd++)
            if (sampleInputs_[sInd*nInputs_+inputID] > Xmax)
               Xmax = sampleInputs_[sInd*nInputs_+inputID];
         printf("Xmin and Xmax found = %e %e.\n", Xmin, Xmax);
         sprintf(pString,"Enter the lower threshold (Xmin = %e) : ",Xmin);
         threshL = getDouble(pString);
         sprintf(pString,"Enter the upper threshold (Xmax = %e) : ",Xmax);
         threshU = getDouble(pString);
         if (threshL >= threshU)
         {
            printf("ERROR: Lower bound should be < upper bound.\n");
            continue;
         }
         for (sInd = 0; sInd < nSamples_; sInd++)
         {
            if (sampleInputs_[sInd*nInputs_+inputID] < threshL ||
                sampleInputs_[sInd*nInputs_+inputID] > threshU)
            {
               tagArray_[sInd] = 0;
            }
         }
         printf("Sample points left after thresholding: [\n");
         for (sInd = 0; sInd < nSamples_; sInd++)
            if (tagArray_[sInd] == 1) printf("%d ", sInd+1);
         printf("]\n");
         printf("itag completed.\n");
      }

      // +++ otag 
      else if (!strcmp(command, "otag"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("otag: tag sample points based on output values\n");
            printf("syntax: otag (no argument needed)\n");
            printf("This command is useful for extracting indices of\n");
            printf("sample points with input values falling in a given\n");
            printf("range. This command can be called multiple times to\n");
            printf("impose multiple filters (AND operation) on the sample.\n");
            printf("This tagged list can be reset to empty by a 'load'.\n");
            printf("Note: this command does not change the data file.\n");
            printf("      Only display sample index information after\n");
            printf("      tagging.\n");
            continue;
         }
         if (psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data (load data first).\n");
            continue;
         }
         if (tagArray_ == NULL && nSamples_ > 0) 
         {
            tagArray_ = new int[nSamples_];
            for (ii = 0; ii < nSamples_; ii++) tagArray_[ii] = 1;
         } 
         outputID = 0;
         sprintf(pString, "Enter output number (1 - %d) : ", nOutputs_);
         outputID = getInt(1, nOutputs_, pString);
         outputID--;
         Ymin = sampleOutputs_[outputID];
         for (sInd = 1; sInd < nSamples_; sInd++)
            if (sampleOutputs_[sInd*nOutputs_+outputID] < Ymin)
               Ymin = sampleOutputs_[sInd*nOutputs_+outputID];
         Ymax = sampleOutputs_[outputID];
         for (sInd = 1; sInd < nSamples_; sInd++)
            if (sampleOutputs_[sInd*nOutputs_+outputID] > Ymax)
               Ymax = sampleOutputs_[sInd*nOutputs_+outputID];
         printf("Ymin and Ymax found = %e %e.\n", Ymin, Ymax);
         sprintf(pString,"Enter the lower threshold (Ymin = %e) : ",Ymin);
         threshL = getDouble(pString);
         sprintf(pString,"Enter the upper threshold (Ymax = %e) : ",Ymax);
         threshU = getDouble(pString);
         if (threshL >= threshU)
         {
            printf("ERROR: Lower bound should be < upper bound.\n");
            continue;
         }
         for (sInd = 0; sInd < nSamples_; sInd++)
         {
            if (sampleOutputs_[sInd*nOutputs_+outputID] < threshL ||
                sampleOutputs_[sInd*nOutputs_+outputID] > threshU)
            {
               tagArray_[sInd] = 0;
            }
         }
         printf("Sample points left after thresholding: [\n");
         for (sInd = 0; sInd < nSamples_; sInd++)
            if (tagArray_[sInd] == 1) printf("%d ", sInd+1);
         printf("]\n");
         printf("otag completed.\n");
      }

      // +++ oreset 
      else if (!strcmp(command, "oreset"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("oreset: reset certain output values with a new one.\n");
            printf("syntax: oreset (no argument needed)\n");
            printf("This command is useful when some sample outputs are to be\n");
            printf("set to a different value (e.g. set undefined to 0).\n");
            continue;
         }
         if (psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data (load data first).\n");
            continue;
         }
         sprintf(pString,"Enter output number (1 - %d) : ", nOutputs_);
         oInd = getInt(1, nInputs_, pString);
         oInd--;
         Ymin = PSUADE_UNDEFINED;
         Ymax = -PSUADE_UNDEFINED;
         for (ii = 0; ii < nSamples_; ii++)
         {
            if (sampleOutputs_[ii*nOutputs_+oInd] < Ymin)
               Ymin = sampleOutputs_[ii*nOutputs_+oInd];
            if (sampleOutputs_[ii*nOutputs_+oInd] > Ymax)
               Ymax = sampleOutputs_[ii*nOutputs_+oInd];
         }
         printf("Lower and upper values for this output : %e %e\n",Ymin, Ymax);
         printf("Now specify the range of current values to be reset: \n");
         sprintf(pString,"Enter the lower bound (inclusive) of this range : ");
         Ymin = getDouble(pString);
         sprintf(pString,"Enter the upper bound (inclusive) of this range : ");
         Ymax = getDouble(pString);
         sprintf(pString,"Enter the desired output value to be set to : ");
         ddata = getDouble(pString);
         for (ii = 0; ii < nSamples_; ii++)
         {
            if (sampleOutputs_[ii*nOutputs_+oInd] >= Ymin &&
                sampleOutputs_[ii*nOutputs_+oInd] <= Ymax)
               sampleOutputs_[ii*nOutputs_+oInd] = ddata;
         }
         psuadeIO_->updateOutputSection(nSamples_,nOutputs_, sampleOutputs_,
                                        sampleStates_, NULL); 
         printf("oreset completed: use 'write' to store.\n");
      }

      // +++ irerange
      else if (!strcmp(command, "irerange"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("irerange: modify a sample input with a different range.\n");
            printf("syntax: irerange (no argument needed)\n");
            printf("Note: Input values will be re-evaluated to the new\n");
            printf("      ranges (that is, they will be re-scaled).\n");
            continue;
         }
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no data (load data first).\n");
            continue;
         }
         sprintf(pString,"Enter input number (1 - %d) : ", nInputs_);
         iInd = getInt(1, nInputs_, pString);
         iInd--;
         printf("Current lower bound for input %d = %24.16e\n",
                iInd+1,iLowerB_[iInd]);
         sprintf(pString,"Enter new lower bound for input %d : ",iInd+1);
         Xmin = getDouble(pString);
         printf("Current upper bound for input %d = %24.16e\n",
                iInd+1,iUpperB_[iInd]);
         sprintf(pString,"Enter new upper bound for input %d : ",iInd+1);
         Xmax = getDouble(pString);
         if (Xmin >= Xmax)
         {
            printf("ERROR: lower bound >= upper bound.\n");
            continue;
         }
         for (ii = 0; ii < nSamples_; ii++)
         {
            dtemp = sampleInputs_[ii*nInputs_+iInd];
            dtemp = (dtemp-iLowerB_[iInd]) / (iUpperB_[iInd]-iLowerB_[iInd]);
            sampleInputs_[ii*nInputs_+iInd] = dtemp * (Xmax - Xmin) + Xmin;
         }
         iLowerB_[iInd] = Xmin;
         iUpperB_[iInd] = Xmax;
         psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,iLowerB_,iUpperB_,
                                       sampleInputs_,NULL,NULL,NULL,NULL,NULL); 
         if (currSession != NULL) delete currSession;
         currSession = new PsuadeSession();
         psuadeIO_->getSession(currSession);
         printf("irerange completed: use 'write' to store.\n");
      }

      // +++ ireset 
      else if (!strcmp(command, "ireset"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("ireset: map an input parameter to some distinct value\n");
            printf("syntax: ireset (no argument needed)\n");
            printf("This command is useful when a sample input is to be\n");
            printf("re-mapped from intervals to distinct values.\n");
            continue;
         }
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no data (load data first).\n");
            continue;
         }
         sprintf(pString,"Enter input number (1 - %d) : ", nInputs_);
         iInd = getInt(1, nInputs_, pString);
         iInd--;
         Xmin = PSUADE_UNDEFINED;
         Xmax = -PSUADE_UNDEFINED;
         for (ii = 0; ii < nSamples_; ii++)
         {
            if (sampleInputs_[ii*nInputs_+iInd] < Xmin)
               Xmin = sampleInputs_[ii*nInputs_+iInd];
            if (sampleInputs_[ii*nInputs_+iInd] > Xmax)
               Xmax = sampleInputs_[ii*nInputs_+iInd];
         }
         printf("Lower and upper values for this input : %e %e\n",Xmin, Xmax);
         printf("Now specify the range of current values to be reset.\n");
         sprintf(pString,"Enter the lower bound (inclusive) of this range : ");
         Xmin = getDouble(pString);
         sprintf(pString,"Enter the upper bound (inclusive) of this range : ");
         Xmax = getDouble(pString);
         sprintf(pString,"Enter the desired input value to be set to : ");
         ddata = getDouble(pString);
         for (jj = 0; jj < nSamples_; jj++)
         {
            if (sampleInputs_[jj*nInputs_+iInd] >= Xmin &&
                sampleInputs_[jj*nInputs_+iInd] <= Xmax)
               sampleInputs_[jj*nInputs_+iInd] = ddata;
         }
         psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,NULL,NULL,
                                       sampleInputs_,NULL,NULL,NULL,NULL,NULL); 
         if (currSession != NULL) delete currSession;
         currSession = new PsuadeSession();
         psuadeIO_->getSession(currSession);
         printf("ireset completed: use 'write' to store.\n");
      }

      // +++ ifloor 
      else if (!strcmp(command, "ifloor"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("ifloor: map an input to nearest smaller integer\n");
            printf("syntax: ifloor (no argument needed)\n");
            printf("This command is useful when some of the inputs are\n");
            printf("discrete. PSUADE creates samples based on continuous\n");
            printf("variables. This command helps to modify the samples\n");
            printf("to accommodate discrete variables.\n");
            continue;
         }
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no data (load data first).\n");
            continue;
         }
         sprintf(pString,"Enter input number (1 - %d) : ", nInputs_);
         iInd = getInt(1, nInputs_, pString);
         iInd--;
         for (jj = 0; jj < nSamples_; jj++)
         {
            kk = (int) sampleInputs_[jj*nInputs_+iInd];
            sampleInputs_[jj*nInputs_+iInd] = (double) kk;
         }
         psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,NULL,NULL,
                                       sampleInputs_,NULL,NULL,NULL,NULL,NULL); 
         if (currSession != NULL) delete currSession;
         currSession = new PsuadeSession();
         psuadeIO_->getSession(currSession);
         printf("ifloor completed: use 'write' to store.\n");
      }

      // +++ iceil 
      else if (!strcmp(command, "iceil"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("iceil: map an input to nearest larger integer\n");
            printf("syntax: iceil (no argument needed)\n");
            printf("This command is useful when some of the inputs are\n");
            printf("discrete. PSUADE creates samples based on continuous\n");
            printf("variables. This command helps to modify the samples\n");
            printf("to accommodate discrete variables.\n");
            continue;
         }
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no data (load data first).\n");
            continue;
         }
         sprintf(pString,"Enter input number (1 - %d) : ", nInputs_);
         iInd = getInt(1, nInputs_, pString);
         iInd--;
         for (jj = 0; jj < nSamples_; jj++)
         {
            ddata = sampleInputs_[jj*nInputs_+iInd];
            kk = (int) ddata;
            if (kk != ddata) kk++;
            sampleInputs_[jj*nInputs_+iInd] = (double) kk;
         }
         psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,NULL,NULL,
                                       sampleInputs_,NULL,NULL,NULL,NULL,NULL); 
         if (currSession != NULL) delete currSession;
         currSession = new PsuadeSession();
         psuadeIO_->getSession(currSession);
         printf("iceil completed: use 'write' to store.\n");
      }

      // +++ iround 
      else if (!strcmp(command, "iround"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("iround: map an input to nearest integer\n");
            printf("syntax: iround (no argument needed)\n");
            printf("This command is useful when some of the inputs are\n");
            printf("discrete. PSUADE creates samples based on continuous\n");
            printf("variables. This command helps to modify the samples\n");
            printf("to accommodate discrete variables.\n");
            continue;
         }
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no data (load data first).\n");
            continue;
         }
         sprintf(pString,"Enter input number (1 - %d) : ", nInputs_);
         iInd = getInt(1, nInputs_, pString);
         iInd--;
         for (jj = 0; jj < nSamples_; jj++)
         {
            ddata = sampleInputs_[jj*nInputs_+iInd];
            kk = (int) ddata;
            if (kk != ddata) kk = (int) (ddata + 0.5);
            sampleInputs_[jj*nInputs_+iInd] = (double) kk;
         }
         psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,NULL,NULL,
                                    sampleInputs_,NULL,NULL,NULL,NULL,NULL); 
         if (currSession != NULL) delete currSession;
         currSession = new PsuadeSession();
         psuadeIO_->getSession(currSession);
         printf("iceil completed: use 'write' to store.\n");
      }

      // +++ itran 
      else if (!strcmp(command, "itran"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("itran: transform an input to log10 or power of 10\n");
            printf("syntax: itran (no argument needed)\n");
            continue;
         }
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no data (load data first).\n");
            continue;
         }
         sprintf(pString,"Enter input number (1 - %d, 0 - all) : ", nInputs_);
         iInd = getInt(0, nInputs_, pString);
         iInd--;
         sprintf(pString,"Enter type: (1) log10(input), (2) 10^(input) : ");
         kk = getInt(1, 2, pString);
         if (iInd >= 0)
         {
            if (kk == 1)
            {
               for (ss = 0; ss < nSamples_; ss++)
                  sampleInputs_[ss*nInputs_+iInd]=log10(sampleInputs_[ss*nInputs_+iInd]);
            }
            else
            {
               for (ss = 0; ss < nSamples_; ss++)
                  sampleInputs_[ss*nInputs_+iInd] = 
                         pow(10.0,sampleInputs_[ss*nInputs_+iInd]);
            }
            if (iLowerB_ != NULL)
            {
               iLowerB_[iInd] = PSUADE_UNDEFINED;
               iUpperB_[iInd] = -PSUADE_UNDEFINED;
               inputPDFs_[iInd] = 0;
               inputMeans_[iInd] = 0;
               inputStds_[iInd] = 0;
               for (ss = 0; ss < nSamples_; ss++)
               {
                  if (sampleInputs_[ss*nInputs_+iInd] < iLowerB_[iInd])
                     iLowerB_[iInd] = sampleInputs_[ss*nInputs_+iInd];
                  if (sampleInputs_[ss*nInputs_+iInd] > iUpperB_[iInd])
                     iUpperB_[iInd] = sampleInputs_[ss*nInputs_+iInd];
               }
            }
            for (ii = 0; ii < nInputs_; ii++)
            {
               inputCMat_->setEntry(iInd,ii,0.0);
               inputCMat_->setEntry(ii,iInd,0.0);
            }
            inputCMat_->setEntry(iInd,iInd,1.0);
         }
         else
         {
            if (kk == 1)
            {
               for (ss = 0; ss < nSamples_*nInputs_; ss++)
                  sampleInputs_[ss] = log10(sampleInputs_[ss]);
            }
            else
            {
               for (ss = 0; ss < nSamples_*nInputs_; ss++)
                  sampleInputs_[ss] = pow(10.0,sampleInputs_[ss]);
            }
            if (iLowerB_ != NULL)
            {
               for (ii = 0; ii < nInputs_; ii++)
               {
                  inputPDFs_[ii] = 0;
                  inputMeans_[ii] = 0;
                  inputStds_[ii] = 0;
                  iLowerB_[ii] = PSUADE_UNDEFINED;
                  iUpperB_[ii] = -PSUADE_UNDEFINED;
                  for (ss = 0; ss < nSamples_; ss++)
                  {
                     if (sampleInputs_[ss*nInputs_+ii] < iLowerB_[ii])
                        iLowerB_[ii] = sampleInputs_[ss*nInputs_+ii];
                     if (sampleInputs_[ss*nInputs_+ii] > iUpperB_[ii])
                        iUpperB_[ii] = sampleInputs_[ss*nInputs_+ii];
                  }
               }
            }
            for (ii = 0; ii < nInputs_; ii++)
            {
               for (jj = 0; jj < nInputs_; jj++)
               {
                  if (ii == jj) inputCMat_->setEntry(ii,jj,1.0);
                  else          inputCMat_->setEntry(ii,jj,0.0);
               }
            }
         }
         psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,iLowerB_,iUpperB_,
                             sampleInputs_,NULL,inputPDFs_,inputMeans_,
                             inputStds_,inputCMat_); 
         if (currSession != NULL) delete currSession;
         currSession = new PsuadeSession();
         psuadeIO_->getSession(currSession);
         printf("Input transformation completed: use 'write' to store.\n");
      }

      // +++ otran 
      else if (!strcmp(command, "otran"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("otran: transform an output to log10 or power of 10\n");
            printf("syntax: otran (no argument needed)\n");
            continue;
         }
         if (psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data (load data first).\n");
            continue;
         }
         sprintf(pString,"Enter output number (1 - %d) : ", nOutputs_);
         oInd = getInt(1, nOutputs_, pString);
         oInd--;
         sprintf(pString,"Enter type: (1) log10(output), (2) 10^(output) : ");
         kk = getInt(1, 2, pString);
         if (kk == 1)
         {
            for (ss = 0; ss < nSamples_; ss++)
               sampleOutputs_[ss*nOutputs_+oInd] = 
                        log10(sampleOutputs_[ss*nOutputs_+oInd]);
         }
         else
         {
            for (ss = 0; ss < nSamples_; ss++)
               sampleOutputs_[ss*nOutputs_+oInd] = 
                        pow(10.0,sampleOutputs_[ss*nOutputs_+oInd]);
         }
         psuadeIO_->updateOutputSection(nSamples_,nOutputs_,
                                        sampleOutputs_, sampleStates_, NULL); 
         if (currSession != NULL) delete currSession;
         currSession = new PsuadeSession();
         psuadeIO_->getSession(currSession);
         printf("Output transformation completed: use 'write' to store.\n");
      }

      // +++ iotrace 
      else if (!strcmp(command, "iotrace"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("iotrace: trace input/output values for each run\n");
            printf("syntax: iotrace (no argument needed)\n");
            printf("This command is useful for tracing which values of\n");
            printf("inputs correspond to which values of the outputs.\n");
            printf("Note: one line for each sample point connecting all\n");
            printf("      inputs and outputs.\n");
            continue;
         }
         if (psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data (load data first).\n");
            continue;
         }
         sprintf(pString,"Number of inputs to plot (1 - %d) : ", nInputs_);
         kk = getInt(1, nInputs_, pString);
         if (kk < nInputs_)
         {
            indSet = new int[kk];
            for (ii = 0; ii < kk; ii++)
            {
               sprintf(pString,"Enter the %d-th input (1 - %d): ",
                       ii+1, nInputs_);
               indSet[ii] = getInt(1,nInputs_,pString);
               indSet[ii]--;
            }
         }
         else
         {
            kk = nInputs_;
            indSet = new int[nInputs_];
            for (ii = 0; ii < nInputs_; ii++) indSet[ii] = ii;
         }
         sprintf(pString,"Number of outputs to plot (1 - %d) : ",nOutputs_);
         ll = getInt(1, nOutputs_, pString);
         if (ll < nOutputs_)
         {
            tempI = new int[ll];
            for (ii = 0; ii < ll; ii++)
            {
               sprintf(pString,"Enter the %d-th output (1 - %d): ",
                       ii+1, nOutputs_);
               tempI[ii] = getInt(1,nOutputs_,pString);
               tempI[ii]--;
            }
         }
         else
         {
            ll = nOutputs_;
            tempI = new int[nOutputs_];
            for (ii = 0; ii < nOutputs_; ii++) tempI[ii] = ii;
         }
         if (psPlotTool_ == 1)
         {
            fp = fopen("scilabiotrace.sci", "w");
            if (fp == NULL)
            {
               printf("ERROR: cannot open file scilabiotrace.sci\n");
               continue;
            }
         }
         else
         {
            fp = fopen("matlabiotrace.m","w");
            if (fp == NULL)
            {
               printf("ERROR: cannot open file matlabiotrace.m.\n");
               continue;
            }
         }
         tempV = new double[nOutputs_];
         tempW = new double[nOutputs_];
         for (ii = 0; ii < nOutputs_; ii++)
         {
            tempV[ii] =  PSUADE_UNDEFINED;
            tempW[ii] = -PSUADE_UNDEFINED;
            for (ss = 0; ss < nSamples_; ss++)
            {
               ddata = sampleOutputs_[ss*nOutputs_+ii];
               tempV[ii] = (ddata < tempV[ii]) ? ddata : tempV[ii];
               tempW[ii] = (ddata > tempW[ii]) ? ddata : tempW[ii];
            }
         }
         fprintf(fp, "Y = 1:%d;\n", kk+ll); 
         fprintf(fp, "X = [\n"); 
         for (ss = 0; ss < nSamples_; ss++)
         {
            for (ii = 0; ii < kk; ii++)
            {
               ddata = (sampleInputs_[ss*nInputs_+indSet[ii]] - 
                        iLowerB_[indSet[ii]]) /
                       (iUpperB_[indSet[ii]] - iLowerB_[indSet[ii]]);
               fprintf(fp, "%e ", ddata);
            }
            for (ii = 0; ii < ll; ii++)
            {
               ddata = (sampleOutputs_[ss*nOutputs_+tempI[ii]] - tempV[tempI[ii]])/
                       (tempW[tempI[ii]] - tempV[tempI[ii]]);
               fprintf(fp, "%e ", ddata);
            }
            fprintf(fp, "\n"); 
         }
         fprintf(fp, "];\n"); 
         fprintf(fp,"plot(X',Y)\n");
         sprintf(pString,"Input/Output Values");
         fwritePlotXLabel(fp, pString);
         sprintf(pString,"Inputs 1:%d, Outputs %d:%d",kk,kk+1,kk+ll);
         fwritePlotYLabel(fp, pString);
         sprintf(pString,"Trace Inputs/Outputs per sample point");
         fwritePlotTitle(fp, pString);
         fwritePlotAxes(fp);
         if (psPlotTool_ == 1)
              printf("Scilab iotrace is now in scilabiotrace.sci\n");
         else printf("Matlab iotrace is now in matlabiotrace.m\n");
         fclose(fp);
         delete [] indSet;
         delete [] tempI;
         delete [] tempV;
         delete [] tempW;
      }

      // +++ iplot2all 
      else if (!strcmp(command, "iplot2_all"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("iplot2_all: plot all pairs of inputs\n");
            printf("syntax: iplot2_all (no argument needed)\n");
            printf("This command is useful to find out which parameters\n");
            printf("are responsible for the failed runs.\n");
            continue;
         }
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no data (load data first).\n");
            continue;
         }
         sprintf(pString,"Number of inputs to plot (1 - %d) : ", nInputs_);
         kk = getInt(1, nInputs_, pString);
         if (kk < nInputs_)
         {
            indSet = new int[kk];
            for (ii = 0; ii < kk; ii++)
            {
               sprintf(pString,"Enter the %d-th input (1 - %d): ",
                       ii+1, nInputs_);
               indSet[ii] = getInt(1,nInputs_,pString);
               indSet[ii]--;
            }
         }
         else
         {
            kk = nInputs_;
            indSet = new int[nInputs_];
            for (ii = 0; ii < nInputs_; ii++) indSet[ii] = ii;
         }

         if (psPlotTool_ == 1)
         {
            fp = fopen("scilabiplt2all.sci", "w");
            if (fp == NULL)
            {
               printf("ERROR: cannot open file scilabiplt2all.sci\n");
               continue;
            }
         }
         else
         {
            fp = fopen("matlabiplt2all.m","w");
            if (fp == NULL)
            {
               printf("ERROR: cannot open file matlabiplt2all.m.\n");
               continue;
            }
         }
         sprintf(pString," plotMode  = 0 : plot all in a single plot");
         fwriteComment(fp, pString);
         sprintf(pString," plotMode ~= 0 : plot one at a time");
         fwriteComment(fp, pString);
         sprintf(pString," ranflag: to distinguish between identical points");
         fwriteComment(fp, pString);
         sprintf(pString,"         by adding a small perturbation (when on)");
         fwriteComment(fp, pString);
         fprintf(fp,"plotMode = 0;\n");
         fprintf(fp,"ranflag  = 0;\n");
         for (ii = 0; ii < kk; ii++)
         {
            fprintf(fp, "X%d = [\n", ii+1); 
            for (jj = 0; jj < nSamples_; jj++)
               fprintf(fp, "%e\n", sampleInputs_[jj*nInputs_+indSet[ii]]);
            fprintf(fp, "];\n"); 
         }
         fprintf(fp, "S = [\n"); 
         for (jj = 0; jj < nSamples_; jj++)
            fprintf(fp, "%d\n", sampleStates_[jj]);
         fprintf(fp, "];\n"); 
         fprintf(fp, "X = [\n"); 
         for (ii = 0; ii < kk; ii++) fprintf(fp, "X%d ", ii+1); 
         fprintf(fp, "S];\n"); 
         for (ii = 0; ii < kk; ii++)
         {
            for (jj = ii; jj < kk; jj++)
            {
               fprintf(fp,"if plotMode == 0\n");
               fprintf(fp,"subplot(%d,%d,%d)\n",kk,kk,ii*kk+jj+1);
               fprintf(fp,"else\n");
               if (ii + jj > 0)
               {
                  fprintf(fp,"pause\n");
                  fprintf(fp,"disp('Press enter to continue')\n");
               }
               fwritePlotCLF(fp);
               fprintf(fp,"end;\n");
               if (psPlotTool_ == 1)
               {
                  fprintf(fp,"iset = find(S==0);\n");
                  fprintf(fp,"plot(X(iset,%d),X(iset,%d),'rX','MarkerSize',13)\n",
                          jj+1,ii+1);
                  fprintf(fp,"iset = find(S~=0);\n");
                  fprintf(fp,"plot(X(iset,%d),X(iset,%d),'bX','MarkerSize',13)\n",
                          jj+1,ii+1);
               }
               else
               {
                  fprintf(fp, "iset = find(S==0);\n");
                  fprintf(fp, "plot(X(iset,%d).*(1+ranflag*", jj+1);
                  fprintf(fp, "rand(size(iset,1),1)/100),X(iset,%d)",
                          ii+1);
                  fprintf(fp, ".*(1+ranflag*rand(size(iset,1),1)/100),'rX',");
                  fprintf(fp, "'MarkerSize',13)\n");
                  if (psPlotTool_ == 1)
                       fprintf(fp, "set(gca(),\"auto_clear\",\"off\")\n");
                  else fprintf(fp, "hold on\n");
                  fprintf(fp, "iset = find(S~=0);\n");
                  fprintf(fp, "plot(X(iset,%d).*(1+ranflag*",jj+1);
                  fprintf(fp, "rand(size(iset,1),1)/100),X(iset,%d)",
                          ii+1);
                  fprintf(fp, ".*(1+ranflag*rand(size(iset,1),1)/100),'b*',");
                  fprintf(fp, "'MarkerSize',13)\n");
               }
               fwritePlotXLabel(fp, inputNames_[indSet[jj]]);
               fwritePlotYLabel(fp, inputNames_[indSet[ii]]);
               fwritePlotAxes(fp);
               if (psPlotTool_ == 1)
               {
                  fprintf(fp, "set(gca(),\"auto_clear\",\"on\")\n");
                  fprintf(fp, "a=gca();\n");
                  fprintf(fp, "a.data_bounds=[%e,%e;%e,%e];\n",
                          iLowerB_[indSet[jj]], iLowerB_[indSet[ii]],
                          iUpperB_[indSet[jj]], iUpperB_[indSet[ii]]);
               }
               else
               {
                  fprintf(fp, "hold off\n");
                  fprintf(fp,"axis([%e %e %e %e])\n",
                          iLowerB_[indSet[jj]], iUpperB_[indSet[jj]],
                          iLowerB_[indSet[ii]], iUpperB_[indSet[ii]]);
               }
               fprintf(fp,"disp('red X: failed runs.')\n");
            }
         }
         if (psPlotTool_ == 1)
              printf("The Scilab file is in scilabiplt2all.sci\n");
         else printf("The Matlab file is in matlabiplt2all.m\n");
         fclose(fp);
         delete [] indSet;
      }

      // +++ ihist 
      else if (!strcmp(command, "ihist"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("ihist: plot histogram for an input\n");
            printf("syntax: ihist (no argument needed)\n");
            printf("This command is useful to examine input distribution\n");
            printf("from a posterior sample, e.g. MCMCPostSample (loaded\n");
            printf("using 'iread').\n");
            continue;
         }
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no data (load data first).\n");
            continue;
         }
         sprintf(pString,"Which input to generate histogram (1 - %d) : ",
                 nInputs_);
         kk = getInt(1, nInputs_, pString);
         kk--;
         if (psPlotTool_ == 1) fp = fopen("scilabihist.sci", "w");
         else                  fp = fopen("matlabihist.m", "w");
         fprintf(fp, "X = [\n");
         for (ss = 0; ss < nSamples_; ss++)
            fprintf(fp, "%e \n",sampleInputs_[ss*nInputs_+kk]);
         fprintf(fp, "];\n");
         fwritePlotCLF(fp);
         fwritePlotFigure(fp, 1);
         if (psPlotTool_ == 1)
         {
            fprintf(fp, "histplot(10, X, style=2);\n");
            fprintf(fp, "a = gce();\n");
            fprintf(fp, "a.children.fill_mode = \"on\";\n");
            fprintf(fp, "a.children.thickness = 2;\n");
            fprintf(fp, "a.children.foreground = 0;\n");
            fprintf(fp, "a.children.background = 2;\n");
         }
         else
         {
            fprintf(fp, "[nk,xk] = hist(X,10);\n");
            fprintf(fp, "bar(xk,nk/%d,1.0)\n",nSamples_);
         }
         fwritePlotAxes(fp);
         sprintf(pString,"Sample Histogram for %s",inputNames_[kk]);
         fwritePlotTitle(fp, pString);
         fwritePlotXLabel(fp, "Input Value");
         fwritePlotYLabel(fp, "Probabilities");
         fclose(fp);
         if (psPlotTool_ == 1) printf("Histogram is available in scilabihist.sci\n");
         else                  printf("Histogram is available in matlabihist.m\n");
      }
         
      // +++ ihist2 
      else if (!strcmp(command, "ihist2"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("ihist2: plot histogram for two inputs\n");
            printf("syntax: ihist2 (no argument needed)\n");
            printf("This command is useful to examine input distributions\n");
            printf("from a posterior sample, e.g. MCMCPostSample (loaded\n");
            printf("using 'iread').\n");
            continue;
         }
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no data (load data first).\n");
            continue;
         }
         if (nInputs_ < 2)
         {
            printf("ERROR: this command is for nInputs_ >= 2.\n");
            continue;
         }
         sprintf(pString,"Select the first input (1 - %d) : ", nInputs_);
         jj = getInt(1, nInputs_, pString);
         jj--;
         sprintf(pString,"Select the second input (1 - %d) : ", nInputs_);
         kk = getInt(1, nInputs_, pString);
         kk--;
         if (psPlotTool_ == 1)
         {
            printf("ERROR: scilab plot not currently available.\n");
            continue;
         }
         fp = fopen("matlabihist2.m", "w");
         fprintf(fp, "X = [\n");
         for (ss = 0; ss < nSamples_; ss++)
            fprintf(fp, "%e %e\n",sampleInputs_[ss*nInputs_+kk],
                    sampleInputs_[ss*nInputs_+kk]);
         fprintf(fp, "];\n");
         fprintf(fp, "X1 = X(:,1);\n");
         fprintf(fp, "X2 = X(:,2);\n");
         fprintf(fp, "nb = 20;\n");
         fprintf(fp, "bins = zeros(nb,nb);\n");
         fprintf(fp, "x1min = %e;\n",iLowerB_[jj]);
         fprintf(fp, "x1max = %e;\n",iUpperB_[jj]);
         fprintf(fp, "if (x1min == x1max)\n");
         fprintf(fp, "   if (x1min == 0)\n");
         fprintf(fp, "      x1min = -0.1;\n");
         fprintf(fp, "      x1max = 0.1;\n");
         fprintf(fp, "   else\n");
         fprintf(fp, "      x1wid = x1max;\n");
         fprintf(fp, "      x1min = x1min - x1wid * 0.1;\n");
         fprintf(fp, "      x1max = x1max + x1wid * 0.1;\n");
         fprintf(fp, "   end\n");
         fprintf(fp, "end\n");
         fprintf(fp, "x1wid = (x1max - x1min) / (nb - 1);\n");
         fprintf(fp, "x2min = %e;\n",iLowerB_[kk]);
         fprintf(fp, "x2max = %e;\n",iUpperB_[kk]);
         fprintf(fp, "if (x2min == x2max)\n");
         fprintf(fp, "   if (x2min == 0)\n");
         fprintf(fp, "      x2min = -0.1;\n");
         fprintf(fp, "      x2max = 0.1;\n");
         fprintf(fp, "   else\n");
         fprintf(fp, "      x2wid = x2max - x2min;\n");
         fprintf(fp, "      x2min = x2min - x2wid * 0.1;\n");
         fprintf(fp, "      x2max = x2max + x2wid * 0.1;\n");
         fprintf(fp, "   end\n");
         fprintf(fp, "end\n");
         fprintf(fp, "x2wid = (x2max - x2min) / (nb - 1);\n");
         fprintf(fp, "for ii = 1 : nb\n");
         fprintf(fp, "   x1l = x1min + (ii-1) * x1wid;\n");
         fprintf(fp, "   x1u = x1min + ii * x1wid;\n");
         fprintf(fp, "   if ii == nb\n");
         fprintf(fp, "      iset = find(X1 >= x1l & X1 <= x1u);\n");
         fprintf(fp, "   else\n");
         fprintf(fp, "      iset = find(X1 >= x1l & X1 < x1u);\n");
         fprintf(fp, "   end\n");
         fprintf(fp, "   if (length(iset) > 0)\n");
         fprintf(fp, "      XX = X2(iset);\n");
         fprintf(fp, "      for jj = 1 : nb\n");
         fprintf(fp, "         x2l = x2min + (jj-1) * x2wid;\n");
         fprintf(fp, "         x2u = x2min + jj * x2wid;\n");
         fprintf(fp, "         if jj == nb\n");
         fprintf(fp, "            jset = find(XX >= x2l & XX <= x2u);\n");
         fprintf(fp, "         else\n");
         fprintf(fp, "            jset = find(XX >= x2l & XX < x2u);\n");
         fprintf(fp, "         end\n");
         fprintf(fp, "         bins(ii,jj) = length(jset);\n");
         fprintf(fp, "      end\n");
         fprintf(fp, "   end\n");
         fprintf(fp, "end\n");
         fwritePlotCLF(fp);
         fwritePlotFigure(fp, 1);
         fprintf(fp,"imagesc(bins)\n");
         fprintf(fp,"xtick = x1min:(x1max-x1min)/4:x1max;\n");
         fprintf(fp,"set(gca,'XTick',0:nb/4:nb);\n");
         fprintf(fp,"set(gca,'XTickLabel', xtick);\n");
         fprintf(fp,"ytick = x2min:(x2max-x2min)/4:x2max;\n");
         fprintf(fp,"set(gca,'YTick',0:nb/4:nb);\n");
         fprintf(fp,"set(gca,'YTickLabel', ytick);\n");
         fprintf(fp,"set(gca,'YDir', 'normal');\n");
         fprintf(fp,"xlabel('%s','FontWeight','bold','FontSize',12)\n",
                 inputNames_[jj]);
         fprintf(fp,"ylabel('%s','FontWeight','bold','FontSize',12)\n",
                 inputNames_[kk]);
         fwritePlotAxesNoGrid(fp);
         sprintf(pString,"Sample Input Histogram");
         fwritePlotTitle(fp, pString);
         fclose(fp);
         printf("Histogram is available in matlabihist2.m\n");
      }

      // +++ ohist 
      else if (!strcmp(command, "ohist"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("ohist: plot histogram for an output\n");
            printf("syntax: ohist (no argument needed)\n");
            printf("This command is useful to examine output distribution\n");
            printf("without using the 'ua' command.\n");
            continue;
         }
         if (psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data (load data first).\n");
            continue;
         }
         sprintf(pString,"Which output to generate histogram (1 - %d) : ",
                 nOutputs_);
         kk = getInt(1, nOutputs_, pString);
         kk--;
         if (psPlotTool_ == 1) fp = fopen("scilabohist.sci", "w");
         else                  fp = fopen("matlabohist.m", "w");
         fprintf(fp, "Y = [\n");
         for (ss = 0; ss < nSamples_; ss++)
            fprintf(fp, "%e \n",sampleOutputs_[ss*nOutputs_+kk]);
         fprintf(fp, "];\n");
         fwritePlotCLF(fp);
         fwritePlotFigure(fp, 1);
         if (psPlotTool_ == 1)
         {
            fprintf(fp, "histplot(10, Y, style=2);\n");
            fprintf(fp, "a = gce();\n");
            fprintf(fp, "a.children.fill_mode = \"on\";\n");
            fprintf(fp, "a.children.thickness = 2;\n");
            fprintf(fp, "a.children.foreground = 0;\n");
            fprintf(fp, "a.children.background = 2;\n");
         }
         else
         {
            fprintf(fp, "[nk,yk] = hist(Y,10);\n");
            fprintf(fp, "bar(yk,nk/%d,1.0)\n",nSamples_);
         }
         fwritePlotAxes(fp);
         sprintf(pString,"Sample Output Histogram for %s",outputNames_[kk]);
         fwritePlotTitle(fp, pString);
         fwritePlotXLabel(fp, "Output Value");
         fwritePlotYLabel(fp, "Probabilities");
         fclose(fp);
         if (psPlotTool_ == 1) printf("Histogram is available in scilabohist.sci\n");
         else                  printf("Histogram is available in matlabohist.m\n");
      }
         
      // +++ ohist2 
      else if (!strcmp(command, "ohist2"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("ohist2: plot histogram for two outputs\n");
            printf("syntax: ohist2 (no argument needed)\n");
            continue;
         }
         if (psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data (load data first).\n");
            continue;
         }
         if (nOutputs_ < 2)
         {
            printf("ERROR: this command is for nOutputs >= 2.\n");
            continue;
         }
         sprintf(pString,"Select the first input (1 - %d) : ", nOutputs_);
         jj = getInt(1, nOutputs_, pString);
         jj--;
         sprintf(pString,"Select the second input (1 - %d) : ", nOutputs_);
         kk = getInt(1, nOutputs_, pString);
         kk--;
         if (psPlotTool_ == 1)
         {
            printf("ERROR: scilab plot not currently available.\n");
            continue;
         }
         fp = fopen("matlabohist2.m", "w");
         fprintf(fp, "Y = [\n");
         for (ss = 0; ss < nSamples_; ss++)
            fprintf(fp, "%e %e\n",sampleInputs_[ss*nInputs_+kk],
                    sampleInputs_[ss*nInputs_+kk]);
         fprintf(fp, "];\n");
         fprintf(fp, "Y1 = Y(:,1);\n");
         fprintf(fp, "Y2 = Y(:,2);\n");
         fprintf(fp, "nb = 20;\n");
         fprintf(fp, "bins = zeros(nb,nb);\n");
         fprintf(fp, "y1min = min(Y1);\n");
         fprintf(fp, "y1max = max(Y1);\n");
         fprintf(fp, "if (y1min == y1max)\n");
         fprintf(fp, "   if (y1min == 0)\n");
         fprintf(fp, "      y1min = -0.1;\n");
         fprintf(fp, "      y1max = 0.1;\n");
         fprintf(fp, "   else\n");
         fprintf(fp, "      y1wid = y1max;\n");
         fprintf(fp, "      y1min = y1min - y1wid * 0.1;\n");
         fprintf(fp, "      y1max = y1max + y1wid * 0.1;\n");
         fprintf(fp, "   end\n");
         fprintf(fp, "end\n");
         fprintf(fp, "y1wid = (y1max - y1min) / (nb - 1);\n");
         fprintf(fp, "y2min = min(Y2);\n");
         fprintf(fp, "y2max = max(Y2);\n");
         fprintf(fp, "if (y2min == y2max)\n");
         fprintf(fp, "   if (y2min == 0)\n");
         fprintf(fp, "      y2min = -0.1;\n");
         fprintf(fp, "      y2max = 0.1;\n");
         fprintf(fp, "   else\n");
         fprintf(fp, "      y2wid = y2max - y2min;\n");
         fprintf(fp, "      y2min = y2min - y2wid * 0.1;\n");
         fprintf(fp, "      y2max = y2max + y2wid * 0.1;\n");
         fprintf(fp, "   end\n");
         fprintf(fp, "end\n");
         fprintf(fp, "y2wid = (y2max - y2min) / (nb - 1);\n");
         fprintf(fp, "for ii = 1 : nb\n");
         fprintf(fp, "   y1l = y1min + (ii-1) * y1wid;\n");
         fprintf(fp, "   y1u = y1min + ii * y1wid;\n");
         fprintf(fp, "   if ii == nb\n");
         fprintf(fp, "      iset = find(Y1 >= y1l & Y1 <= y1u);\n");
         fprintf(fp, "   else\n");
         fprintf(fp, "      iset = find(Y1 >= y1l & Y1 < y1u);\n");
         fprintf(fp, "   end\n");
         fprintf(fp, "   if (length(iset) > 0)\n");
         fprintf(fp, "      YY = Y2(iset);\n");
         fprintf(fp, "      for jj = 1 : nb\n");
         fprintf(fp, "         y2l = y2min + (jj-1) * y2wid;\n");
         fprintf(fp, "         y2u = y2min + jj * y2wid;\n");
         fprintf(fp, "         if jj == nb\n");
         fprintf(fp, "            jset = find(YY >= y2l & YY <= y2u);\n");
         fprintf(fp, "         else\n");
         fprintf(fp, "            jset = find(YY >= y2l & YY < y2u);\n");
         fprintf(fp, "         end\n");
         fprintf(fp, "         bins(ii,jj) = length(jset);\n");
         fprintf(fp, "      end\n");
         fprintf(fp, "   end\n");
         fprintf(fp, "end\n");
         fwritePlotCLF(fp);
         fwritePlotFigure(fp, 1);
         fprintf(fp,"imagesc(bins)\n");
         fprintf(fp,"xtick = y1min:(y1max-y1min)/4:y1max;\n");
         fprintf(fp,"set(gca,'XTick',0:nb/4:nb);\n");
         fprintf(fp,"set(gca,'XTickLabel', xtick);\n");
         fprintf(fp,"ytick = y2min:(y2max-y2min)/4:y2max;\n");
         fprintf(fp,"set(gca,'YTick',0:nb/4:nb);\n");
         fprintf(fp,"set(gca,'YTickLabel', ytick);\n");
         fprintf(fp,"set(gca,'YDir', 'normal');\n");
         fprintf(fp,"xlabel('%s','FontWeight','bold','FontSize',12)\n",
                 outputNames_[jj]);
         fprintf(fp,"ylabel('%s','FontWeight','bold','FontSize',12)\n",
                 outputNames_[kk]);
         fwritePlotAxesNoGrid(fp);
         fwritePlotAxes(fp);
         sprintf(pString,"Sample Output Histogram");
         fwritePlotTitle(fp, pString);
         fclose(fp);
         printf("Histogram is available in matlabohist2.m\n");
      }

      // +++ iplot_pdf 
      else if (!strcmp(command, "iplot_pdf"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("iplot_pdf: plot pdf of selected inputs\n");
            printf("syntax: iplot_pdf (no argument needed)\n");
            continue;
         }
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no data (load data first).\n");
            continue;
         }
         sprintf(pString,"Number of inputs to plot (1 - %d) : ", nInputs_);
         kk = getInt(1, nInputs_, pString);
         if (kk < nInputs_)
         {
            indSet = new int[kk];
            for (ii = 0; ii < kk; ii++)
            {
               sprintf(pString,"Enter the %d-th input (1 - %d): ",
                       ii+1, nInputs_);
               indSet[ii] = getInt(1,nInputs_,pString);
               indSet[ii]--;
            }
         }
         else
         {
            kk = nInputs_;
            indSet = new int[nInputs_];
            for (ii = 0; ii < nInputs_; ii++) indSet[ii] = ii;
         }

         if (psPlotTool_ == 1)
         {
            fp = fopen("scilabipltpdf.sci", "w");
            if (fp == NULL)
            {
               printf("ERROR: cannot open file scilabipltpdf.sci\n");
               continue;
            }
         }
         else
         {
            fp = fopen("matlabipltpdf.m","w");
            if (fp == NULL)
            {
               printf("ERROR: cannot open file matlabipltpdf.m.\n");
               continue;
            }
         }
         fprintf(fp, "XX = [\n"); 
         for (ii = 0; ii < nSamples_; ii++)
         {
            for (jj = 0; jj < kk; jj++)
               fprintf(fp, "%e ", sampleInputs_[ii*nInputs_+indSet[jj]]);
            fprintf(fp, "\n");
         }
         fprintf(fp, "];\n"); 
         for (ii = 0; ii < kk; ii++)
            fprintf(fp, "X%d = XX(:,%d);\n", ii+1, ii+1); 
         ll = (int) sqrt(1.0*kk);
         if (ll * ll != kk) ll++; 
         fwritePlotCLF(fp);
         for (ii = 0; ii < kk; ii++)
         {
            if (kk > 1)
               fprintf(fp,"subplot(%d,%d,%d)\n",ll,ll,ii+1);
            if (psPlotTool_ == 1)
            {
               fprintf(fp, "X = X%d;\n",ii+1);
               fprintf(fp, "ymin = min(X);\n");
               fprintf(fp, "ymax = max(X);\n");
               fprintf(fp, "ywid = 0.1 * (ymax - ymin);\n");
               fprintf(fp, "if (ywid < 1.0e-12)\n");
               fprintf(fp, "   disp('range too small.')\n");
               fprintf(fp, "   halt\n");
               fprintf(fp, "end;\n");
               fprintf(fp, "histplot(10, X, style=2);\n");
               fprintf(fp, "a = gce();\n");
               fprintf(fp, "a.children.fill_mode = \"on\";\n");
               fprintf(fp, "a.children.thickness = 2;\n");
               fprintf(fp, "a.children.foreground = 0;\n");
               fprintf(fp, "a.children.background = 2;\n");
               fwritePlotAxes(fp);
               fwritePlotTitle(fp, "Input Distribution");
               if (inputNames_ != NULL)
                    fwritePlotXLabel(fp, inputNames_[indSet[ii]]);
               else fwritePlotXLabel(fp, "Input Value");
               fwritePlotYLabel(fp, "Probabilities");
            }
            else
            {
               fprintf(fp, "X = X%d;\n",ii+1);
               fprintf(fp, "[nk,xk]=hist(X,10);\n");
               fprintf(fp, "bar(xk,nk/%d,1.0)\n",nSamples_);
               fwritePlotAxes(fp);
               fwritePlotTitle(fp, "Input Distribution");
               if (inputNames_ != NULL)
                  fwritePlotXLabel(fp, inputNames_[indSet[ii]]);
               else fwritePlotXLabel(fp, "Input Value");
               fprintf(fp, "axis tight\n");
               fwritePlotYLabel(fp, "Probabilities");
            }
         }
         if (psPlotTool_ == 1)
              printf("Plot file is scilabipltpdf.sci\n");
         else printf("Plot file is matlabipltpdf.m\n");
         delete [] indSet;
         indSet = NULL;
      }

      // +++ printlevel 
      else if (!strcmp(command, "printlevel"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("printlevel: set printlevel (0 - 5)\n");
            continue;
         }
         sscanf(lineIn, "%s %d", command, &ii);
         if (ii >= PL_MIN && ii <= PL_MAX)
         {
            outputLevel_ = ii;
            printf("psuade printlevel set to %d\n", ii);
         }
         else
         {
            printf("printlevel out of range; set to %d\n", PL_INTERACTIVE);
            outputLevel_ = PL_INTERACTIVE;
         }
         setPrintLevelTS(outputLevel_);
      }

      // +++ turn on printlevel=4 and interactive
      else if (!strcmp(command, "on"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("on: set printlevel to 4 and turn on interactive mode\n");
            continue;
         }
         outputLevel_ = 4;
         psRSExpertMode_ = 1;
         psAnaExpertMode_ = 1;
         psSamExpertMode_ = 1;
      }

      // +++ gp_sa2 
      else if (!strcmp(command, "gp_sa2"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("gp_sa2: GP SA on sub-divisions of each input\n");
            printf("        (not verified yet)\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data (load data first).\n");
            continue;
         }
         outputID = 0;
         sprintf(pString, "Enter output number (1 - %d) = ", nOutputs_);
         outputID = getInt(1, nOutputs_, pString);
         outputID--;
         faType = PSUADE_RS_GP1;
         printf("This command may be very time-consuming to execute.\n");
         count = nSamples_ / 100;
         if (count > 4) count = 4;
         sprintf(pString, "How many partitions (1 - %d) = ", count);
         nParts = getInt(1, count, pString);
         count = nSamples_ / nParts + 1;
         tempX = new double[count*nInputs_];
         tempY = new double[count];
         tempW = new double[nInputs_];
         for (ii = 0; ii < nInputs_; ii++)
         {
            tempW[ii] = 0.0;
            ddata = (iUpperB_[ii] - iLowerB_[ii]) / (double) nParts;
            for (jj = 0; jj < nParts; jj++)
            {
               count = 0;
               for (kk = 0; kk < nSamples_; kk++)
               {
                  dtemp = sampleInputs_[kk*nInputs_+ii];
                  if (dtemp >= ddata*jj+iLowerB_[ii] &&
                      dtemp <= ddata*(jj+1)+iLowerB_[ii])
                  {
                     for (ind = 0; ind < nInputs_; ind++)
                        tempX[count*nInputs_+ind] = 
                              sampleInputs_[kk*nInputs_+ii];
                     tempY[count++] = sampleOutputs_[kk*nOutputs_+outputID];
                  }
               }
               if (count == 0)
               {
                  printf("ERROR encountered.\n");
                  exit(1);
               }
               if (outputLevel_ > 1)
                 printf("Input %3d Part %3d has size %d\n",ii+1,jj+1,count);
               faPtr = genFA(faType, nInputs_, iOne, count);
               if (faPtr != NULL)
               {
                  faPtr->setBounds(iLowerB_, iUpperB_);
                  faPtr->setOutputLevel(outputLevel_);
                  faPtr->initialize(tempX,tempY);
                  strcpy(pString, "rank");
                  targv[0] = (char *) pString;
                  targv[1] = (char *) &ii;
                  dtemp = faPtr->setParams(2, targv);
                  tempW[ii] += dtemp;
                  delete faPtr;
                  faPtr = NULL;
               }
            }
         }
      }

      // +++ sot_sa2 
      else if (!strcmp(command, "sot_sa2"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("sot_sa2: sum-of-trees screening using boosting\n");
            printf("syntax: sot_sa2 (no argument needed)\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data (load data first).\n");
            continue;
         }
         sprintf(pString, "Enter output number (1 - %d) = ", nOutputs_);
         outputID = getInt(1, nOutputs_, pString);
         outputID--;
         faType = PSUADE_RS_SOTS;
         faPtr  = genFA(faType, nInputs_, iOne, nSamples_);
         if (faPtr != NULL)
         {
            faPtr->setBounds(iLowerB_, iUpperB_);
            faPtr->setOutputLevel(outputLevel_);
            tempY = new double[nSamples_];
            for (ii = 0; ii < nSamples_; ii++)
               tempY[ii] = sampleOutputs_[ii*nOutputs_+outputID];
            status = faPtr->initialize(sampleInputs_,tempY);
            strcpy(pString, "mode1");
            targv[0] = (char *) pString;
            faPtr->setParams(1, targv);
            strcpy(pString, "rank");
            targv[0] = (char *) pString;
            faPtr->setParams(1, targv);
            delete faPtr;
            faPtr = NULL;
            delete [] tempY;
            tempY = NULL;
         }
      }

      // +++ delta_test 
      else if (!strcmp(command, "delta_test"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("delta_test: perform Delta test\n");
            printf("syntax: delta_test (no argument needed)\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data (load data first).\n");
            continue;
         }
         outputID = 0;
         sprintf(pString, "Enter output number (1 - %d) = ", nOutputs_);
         outputID = getInt(1, nOutputs_, pString);
         outputID--;
         analysisMethod = PSUADE_ANA_DTEST;
         anaManager = new AnalysisManager();
         anaManager->setup(analysisMethod, 0);
         psuadeIO_->getParameter("ana_diagnostics",pPtr);
         ii = pPtr.intData_;
         psuadeIO_->updateAnalysisSection(-1,-1,-1,outputLevel_,-1,-1);
         anaManager->analyze(psuadeIO_, 0, NULL, outputID);
         psuadeIO_->updateAnalysisSection(-1,-1,-1,ii,-1,-1);
         delete anaManager;
      }

      // +++ eta_test 
      else if (!strcmp(command, "eta_test"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("eta_test: perform Eta test (a different form of Delta test)\n");
            printf("This test is a derivative of delta_test.\n");
            printf("The power of this test has not been verified.\n");
            printf("syntax: eta_test (no argument needed)\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data (load data first).\n");
            continue;
         }
         outputID = 0;
         sprintf(pString, "Enter output number (1 - %d) = ", nOutputs_);
         outputID = getInt(1, nOutputs_, pString);
         outputID--;
         analysisMethod = PSUADE_ANA_ETEST;
         anaManager = new AnalysisManager();
         anaManager->setup(analysisMethod, 0);
         psuadeIO_->getParameter("ana_diagnostics",pPtr);
         ii = pPtr.intData_;
         psuadeIO_->updateAnalysisSection(-1,-1,-1,outputLevel_,-1,-1);
         anaManager->analyze(psuadeIO_, 0, NULL, outputID);
         psuadeIO_->updateAnalysisSection(-1,-1,-1,ii,-1,-1);
         delete anaManager;
      }

      // +++ gd_test 
      else if (!strcmp(command, "gd_test"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("gd_test: perform Gower distance analysis\n");
            printf("The power of this test has not been verified.\n");
            printf("This test is not ready yet.\n");
            printf("syntax: gd_test (no argument needed)\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data (load data first).\n");
            continue;
         }
         outputID = 0;
         sprintf(pString, "Enter output number (1 - %d) : ", nOutputs_);
         outputID = getInt(1, nOutputs_, pString);
         outputID--;
         analysisMethod = PSUADE_ANA_GOWER;
         anaManager = new AnalysisManager();
         anaManager->setup(analysisMethod, 0);
         psuadeIO_->getParameter("ana_diagnostics",pPtr);
         ii = pPtr.intData_;
         psuadeIO_->updateAnalysisSection(-1,-1,-1,outputLevel_,-1,-1);
         anaManager->analyze(psuadeIO_, 0, NULL, outputID);
         psuadeIO_->updateAnalysisSection(-1,-1,-1,ii,-1,-1);
         delete anaManager;
      }

      // +++ rscreate 
      else if (!strcmp(command, "rscreate"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rscreate: create a response surface with data from 'load'\n");
            printf("syntax: rscreate (no argument needed)\n");
            printf("This command is useful for creating a response surface\n");
            printf("on the fly, and then use rseval new sample points.\n");
            printf("The sample to be used for rscreate should have already\n");
            printf("been loaded to local memory using 'load'.\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data (load data first).\n");
            continue;
         }
         faFlag = 3;
         sprintf(pString, "Which output to use (1 - %d, 0 - all) ? ",nOutputs_);
         outputID = getInt(0, nOutputs_, pString);
         psuadeIO_->getParameter("ana_outputid", pPtr);
         kk = pPtr.intData_;
         faPtrsRsEval = new FuncApprox*[nOutputs_];
         if (outputID != 0)
         {
            psuadeIO_->updateAnalysisSection(-1, -1, -1, outputLevel_, outputID-1, -1);
            faPtrsRsEval[0] = genFAInteractive(psuadeIO_, faFlag);
            for (ii = 1; ii < nOutputs_; ii++) faPtrsRsEval[ii] = NULL;
         }
         else
         {
            for (ii = 0; ii < nOutputs_; ii++)
            {
               printf("Creating response surface for output %d\n", ii+1);
               psuadeIO_->updateAnalysisSection(-1, -1, -1, outputLevel_, ii, -1);
               faPtrsRsEval[ii] = genFAInteractive(psuadeIO_, faFlag);
            }
         }
         psuadeIO_->updateAnalysisSection(-1, -1, -1, -1, kk, -1);
         printf("Now you can use the following to evaluate a new point: \n");
         printf("(1) use rseval with ivec_create/ivec_modify/ivec_show\n");
         printf("(2) use rseval and a file containing the new point\n");
         printf("(3) use rseval_m and a file containing new points\n");
      }

      // +++ rseval 
      else if (!strcmp(command, "rseval"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rseval: evaluate a data point after rscreate\n");
            printf("syntax: rseval (no argument needed)\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data (load data first).\n");
            continue;
         }
         if (faPtrsRsEval == NULL) 
         {
            printf("ERROR: response surface has not been created yet.\n");
            printf("       Please load data and use rscreate first.\n");
            continue;
         }
         count = 0;
         sprintf(pString, "Data taken from a file (n - from register)? (y or n) ");
         getString(pString, winput);
         if (winput[0] == 'y')
         {
            printf("File format: \n");
            printf("PSUADE_BEGIN \n");
            printf("<nPts> <nInputs> \n");
            printf("1   <data> <data> ... \n");
            printf("2   <data> <data> ... \n");
            printf("... <data> <data> ... \n");
            printf("PSUADE_END \n");
            printf("Enter file for evaluation data : ");
            scanf("%s", dataFile);
            fgets(lineIn,5000,stdin); 
            fp = fopen(dataFile, "r");
            if (fp == NULL)
            {
               printf("ERROR: data file %s not found.\n", dataFile);
               delete [] inputSettings;
               continue;
            }
            else
            {
               fscanf(fp, "%s", winput);
               if (strcmp(winput, "PSUADE_BEGIN"))
               {
                  printf("ERROR: file must begin with PSUADE_BEGIN\n");
                  fclose(fp);
                  continue;
               }
               else
               {
                  fscanf(fp, "%d %d", &count, &kk);
                  if (count <= 0)
                  {
                     printf("ERROR: invalid sample size\n");
                     fclose(fp);
                     continue;
                  }
                  if (kk != nInputs_)
                  {
                     printf("ERROR: input size does not match nInputs.\n");
                     printf(":      input size in local memory = %d.\n", 
                            nInputs_);
                     printf(":      input size from file       = %d.\n", kk);
                     count = 0;
                     fclose(fp);
                     continue;
                  }
                  printf("Number of points to be evaluated = %d\n", count);
                  fgets(pString, 5000, fp);
                  while (1)
                  {
                     kk = getc(fp);
                     if (kk != '#')
                     {
                        ungetc(kk, fp);
                        break;
                     }
                     else
                     {
                        fgets(pString, 5000, fp);
                     }
                  }
                  inputSettings = new double[count*nInputs_];
                  for (jj = 0; jj < count; jj++)
                  {
                     fscanf(fp, "%d", &ind);
                     if (ind != (jj+1))
                     {
                        printf("ERROR: input index mismatch (%d,%d)\n",
                               jj+1,ind);
                        printf("       read     index = %d\n", ind);
                        printf("       expected index = %d\n", jj+1);
                        count = 0;
                        delete [] inputSettings;
                        inputSettings = NULL;
                        break;
                     }
                     for (ii = 0; ii < nInputs_; ii++)
                        fscanf(fp, "%lg", &(inputSettings[jj*nInputs_+ii]));
                  }
                  if (count > 0)
                  {
                     fscanf(fp, "%s", winput);
                     fscanf(fp, "%s", winput);
                     if (strcmp(winput, "PSUADE_END"))
                     {
                        fclose(fp);
                        printf("ERROR: file must end with PSUADE_END\n");
                        delete [] inputSettings;
                        count = 0;
                        continue;
                     }
                  }
               }
            }
         }
         else
         {
            inputSettings = new double[nInputs_];
            if (dataReg_ == NULL)
            {
               printf("ERROR: data register has not been created.\n");
               continue;
            }
            count = 1;
            for (ii = 0; ii < nInputs_; ii++) 
               inputSettings[ii] = dataReg_[ii];
         }
         if (count > 0)
         {
            flag = 0;
            printf("Fuzzy evaluation ? (y or n) ");
            fgets(winput,10,stdin); 
            if (winput[0] == 'y') flag = 1;
            printf("Evaluated data written to a file ? (y or n) ");
            fgets(winput,10,stdin); 
            if (winput[0] == 'y')
            {
               fp = fopen("eval_sample","w");
               sprintf(pString," Inputs Out1 <Std1> Out2 <Std2> ...");
               fwriteComment(fp, pString);
               if (nOutputs_ == 1 || (nOutputs_ > 1 && faPtrsRsEval[1] == NULL))
                  fprintf(fp, "%d %d 2\n", count, nInputs_);
               else
                  fprintf(fp, "%d %d %d\n", count, nInputs_, 2*nOutputs_);
            }
            else fp = NULL;
            tempY = new double[nOutputs_*count];
            tempW = new double[nOutputs_*count];
            for (ii = 0 ; ii < nOutputs_; ii++)
            {
               if (faPtrsRsEval[ii] != NULL)
               {
                  if (flag == 1)
                     dtemp = faPtrsRsEval[ii]->evaluatePointFuzzy(count,
                                 inputSettings,&(tempY[ii*count]),&(tempW[ii*count]));
                  else
                     dtemp = faPtrsRsEval[ii]->evaluatePoint(count,
                                 inputSettings,&(tempY[ii*count]));
               }
            }
            for (kk = 0 ; kk < count; kk++)
            {
               if (fp != NULL)
                  for (ii = 0 ; ii < nInputs_; ii++)
                     fprintf(fp, "%e ", inputSettings[kk*nInputs_+ii]);
               printf("Interpolated Point %d: ", kk+1);
               for (ii = 0 ; ii < nOutputs_; ii++)
               {
                  if (faPtrsRsEval[ii] != NULL)
                  {
                     if (flag == 1)
                     {
                        printf("output %d = %e (stdev = %e) ",
                               ii+1,tempY[ii*count+kk], tempW[ii*count+kk]);
                        if (fp != NULL)
                           fprintf(fp,"%e %e ",tempY[ii*count+kk],tempW[ii*count+kk]);
                     }
                     else
                     {
                        printf("output %d = %e ", ii+1,tempY[ii*count+kk]);
                        if (fp != NULL) fprintf(fp,"%e ",tempY[ii*count+kk]);
                     }
                  }
               }
               printf("\n");
               if (fp != NULL) fprintf(fp, "\n");
            }
            if (fp != NULL) printf("Evaluated data file in 'eval_sample'\n");
            if (fp != NULL) fclose(fp);
            fp = NULL;
            delete [] inputSettings;
            delete [] tempY;
            delete [] tempW;
         }
      }

      // +++ rseval_m 
      else if (!strcmp(command, "rseval_m"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rseval_m: evaluate points repeated after rscreate.\n");
            printf("          Data taken from a file with format:\n");
            printf("PSUADE_BEGIN \n");
            printf("<number of points> <nInputs> \n");
            printf("1  <input data 1> <input data 2> .. \n");
            printf("2  <input data 1> <input data 2> .. \n");
            printf("PSUADE_END \n");
            printf("syntax: rseval_m <filename>)\n");
            printf("Note: This command enables PSUADE to be a response\n");
            printf("      surface server. So run this command in the \n");
            printf("      background. Run your program which should put\n");
            printf("      new points to be evaluated in a file called\n");
            printf("      rsevalDataIn.1 in the same directory you run\n");
            printf("      PSUADE, which takes the sample points, evaluatesn");
            printf("      them, and puts the results in rsEvalDataOut.1.\n");
            printf("      You can continue with new points in rsEvalDataIn.2\n");
            printf("      and so on.\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data (load data first).\n");
            continue;
         }
         if (faPtrsRsEval == NULL) 
         {
            printf("ERROR: response surface has not been created yet.\n");
            printf("       Please load data and use rscreate first.\n");
            continue;
         }
         sprintf(dataFile, "psComplete");
         if ((fp=fopen(dataFile,"r")) != NULL)
         {
            printf("Please remove psComplete file and re-do.\n");
            continue;
         }
         inputSettings = new double[nInputs_];
         count = 1;
         while ((fp=fopen(dataFile,"r")) == NULL)
         {
            sprintf(winput, "rsevalDataIn.%d", count);
            printf("Put data in a file called %s.\n", winput);
            while ((fp=fopen(winput,"r")) == NULL);
            fscanf(fp, "%s", winput);
            if (strcmp(winput, "PSUADE_BEGIN"))
            {
               printf("ERROR: file must begin with PSUADE_BEGIN\n");
               fclose(fp);
               delete [] inputSettings;
               continue;
            }
            else
            {
               fscanf(fp, "%d %d", &count2, &kk);
               if (count2 <= 0 || kk != nInputs_)
               {
                  printf("ERROR: invalid data in file.\n");
                  printf("       nData = %d, nInputs = %d\n",count2,kk);
                  fclose(fp);
                  delete [] inputSettings;
                  continue;
               }
               sprintf(winput, "rsevalDataOut.%d", count);
               if ((fpOut=fopen(winput,"w")) == NULL)
               {
                  printf("ERROR: cannot open output file.\n");
                  fclose(fp);
                  delete [] inputSettings;
                  continue;
               }
               for (ss = 0; ss < count2; ss++)
               {
                  fscanf(fp, "%d", &ind);
                  if (ind != ss+1)
                  {
                     printf("ERROR: invalid data in file.\n");
                     printf("       sample number = %d\n", ss+1);
                     fclose(fp);
                     continue;
                  }
                  for (ii = 0; ii < nInputs_; ii++)
                     fscanf(fp, "%lg", &inputSettings[ii]);
                  for (ii = 0; ii < nOutputs_; ii++)
                  {
                     if (faPtrsRsEval[ii] != NULL)
                     {
                        dtemp = faPtrsRsEval[ii]->evaluatePoint(inputSettings);
                        fprintf(fpOut, "%e\n", dtemp);
                     } 
                  } 
               } 
               fclose(fpOut);
               fscanf(fp, "%s", winput);
               fscanf(fp, "%s", winput);
               if (strcmp(winput, "PSUADE_END"))
               {
                  fclose(fp);
                  printf("ERROR: file must end with PSUADE_END\n");
                  continue;
               }
            }
            printf("Evaluated data in file rsevalDataOut.%d\n",count);
            printf("If ready for the next set of evaluations, delete\n");
            printf("both the rsevalDataIn.%d and rsevalDataOut.%d ",
                   count, count);
            printf("files.\n");
            printf("If done, create an empty psComplete file to signal");
            printf(" completion.\n");
            sprintf(winput, "rsevalDataIn.%d", count);
            while ((fp=fopen(winput,"r")) != NULL)
            {
               fclose(fp);
#ifdef WINDOWS
               Sleep(1000);
#else
               sleep(1);
#endif
            }
            sprintf(winput, "rsevalDataOut.%d", count);
            while ((fp=fopen(winput,"r")) != NULL)
            {
               fclose(fp);
#ifdef WINDOWS
               Sleep(1000);
#else
               sleep(1);
#endif
            }
            count++;
         }
         printf("File psComplete detected ==> terminate rseval_m.\n");
         delete [] inputSettings;
      }

      // +++ rs_splot 
      else if (!strcmp(command, "rs_splot"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rs_splot: create a RS-based output scatter plot\n");
            printf("syntax: rs_splot (no argument needed)\n");
            printf("This command is useful when the sample size is small\n");
            printf("and you would like to see the trend with scatter plots.\n");
            printf("The idea is to create a large sample for scatter plots\n");
            printf("using response surfaces from small data set.\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data (load data first).\n");
            continue;
         }
         outputID = 0;
         sprintf(pString, "Enter output number (1 - %d) : ", nOutputs_);
         outputID = getInt(1, nOutputs_, pString);
         outputID--;
         Ymax = - 1.0e35;
         Ymin =   1.0e35;
         for (sInd = 0; sInd < nSamples_; sInd++)
         {
            if (sampleOutputs_[sInd*nOutputs_+outputID] > Ymax)
               Ymax = sampleOutputs_[sInd*nOutputs_+outputID];
            if (sampleOutputs_[sInd*nOutputs_+outputID] < Ymin)
               Ymin = sampleOutputs_[sInd*nOutputs_+outputID];
         }
         sprintf(pString,"Enter the lower bound constraint (Ymin=%e) : ",
                 Ymin);
         threshL = getDouble(pString);
         sprintf(pString,"Enter the upper bound constraint (Ymax=%e) : ",
                 Ymax);
         threshU = getDouble(pString);
         if (threshL >= threshU)
         {
            printf("ERROR: lower bound >= upper bound.\n");
            continue;
         }
         faFlag = 3;
         psuadeIO_->updateAnalysisSection(-1, -1, -1, -1, outputID,-1); 
         faPtr = genFAInteractive(psuadeIO_, faFlag);
         faPtr->setOutputLevel(outputLevel_);
         sprintf(pString,
               "Enter Monte Carlo sample size for probing (max=100000): ");
         kk = getInt(100, 100000, pString);
         sampPtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_MC);
         sampPtr->setPrintLevel(PL_INTERACTIVE);
         sampPtr->setInputBounds(nInputs_, iLowerB_, iUpperB_);
         sampPtr->setOutputParams(1);
         sampPtr->setSamplingParams(kk, 1, 0);
         sampPtr->initialize(0);
         tempX = new double[kk*nInputs_];
         tempY = new double[kk];
         states = new int[kk];
         for (ii = 0; ii < kk; ii++) states[ii] = 1;
         sampPtr->getSamples(kk, nInputs_, 1, tempX, tempY, states);
         count = 0;
         for (ii = 0; ii < kk; ii++)
         {
            dtemp = faPtr->evaluatePoint(&tempX[ii*nInputs_]);
            tempY[ii] = dtemp;
            states[ii] = 1;
            if (dtemp >= threshL && dtemp <= threshU) count++; 
            else                                      states[ii] = 0;
         } 
         winput[0] = '1';
         if (winput[0] == '1')
         {
            if (psPlotTool_ == 1) fp = fopen("scilabrssp1.sci", "w");
            else                  fp = fopen("matlabrssp1.m", "w");
            if (fp == NULL)
            {
               printf("ERROR: cannot open file matlabrssp1.m.\n");
               continue;
            }
            fprintf(fp, "Y = [\n");
            for (sInd = 0; sInd < kk; sInd++)
               if (states[sInd] == 1) fprintf(fp, "%e\n",tempY[sInd]);
            fprintf(fp, "];\n");
            for (iInd = 0; iInd < nInputs_; iInd++)
            {
               fprintf(fp, "X%d = [\n", iInd+1);
               for (sInd = 0; sInd < kk; sInd++)
                  if (states[sInd] == 1)
                     fprintf(fp, "%e\n",tempX[sInd*nInputs_+iInd]);
               fprintf(fp, "];\n");
            }
            for (iInd = 0; iInd < nInputs_; iInd++)
            {
               fwritePlotCLF(fp);
               fprintf(fp, "plot(X%d,Y,'X','MarkerSize',13)\n",iInd+1);
               sprintf(pString, "%s vs %s",outputNames_[outputID],inputNames_[iInd]);
               fwritePlotTitle(fp, pString);
               fwritePlotXLabel(fp, inputNames_[iInd]);
               fwritePlotYLabel(fp, outputNames_[outputID]);
               fwritePlotAxes(fp);
               if (iInd < nInputs_-1) 
               {
                  fprintf(fp, "disp(\'Press enter to go to the next plot\')\n");
                  fprintf(fp, "pause\n");
               }
            }
            fclose(fp);    
            if (psPlotTool_ == 1)
                 printf("scilabrssp1.sci is now ready for scatter plots.\n");
            else printf("matlabrssp1.m is now ready for scatter plots.\n");
         }
         delete faPtr;
         faPtr = NULL;
         delete sampPtr;
         sampPtr = NULL;
         delete [] tempX;
         delete [] tempY;
         delete [] states;
         tempX = NULL;
         tempY = NULL;
         states = NULL;
      }

      // +++ interactive 
      else if (!strcmp(command, "interactive"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("This command turns on/off interactive mode.\n");
            printf("Interactive mode, when turned on, will trigger many\n");
            printf("prompts for additonal information during the course\n");
            printf("of setup and analysis.\n");
            printf("THIS COMMAND HAS BEEN SUPERSEDED BY ana_expert,\n");
            printf("    rs_expert, sam_expert, opt_expert AND io_expert.\n");
            continue;
         }
         if (!strcmp(winput, "on"))  psAnaExpertMode_ = 0;
         if (!strcmp(winput, "off")) psAnaExpertMode_ = 1;
         if (psAnaExpertMode_ == 0)
         {
            psAnaExpertMode_ = 1;
            printf("analysis expert mode on\n");
         }
         else
         {
            psAnaExpertMode_ = 0;
            printf("analysis expert mode off\n");
         }
      }

      // +++ dontask 
      else if (!strcmp(command, "dontask"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("This command turns on/off interactive mode.\n");
            printf("Interactive mode, when turned on, will trigger many\n");
            printf("prompts for additonal information during the course\n");
            printf("of setup and analysis.\n");
            continue;
         }
         psInteractive_ = 0;
      }

      // +++ io_expert 
      else if (!strcmp(command, "io_expert"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("This command turns on/off I/O expert mode.\n");
            printf("This mode, when turned on, will trigger prompts for\n");
            printf("additonal information regarding data input/output.\n");
            printf("When this mode is off, default setting will be used.\n");
            continue;
         }
         if (!strcmp(winput, "on"))  psIOExpertMode_ = 0;
         if (!strcmp(winput, "off")) psIOExpertMode_ = 1;
         if (psIOExpertMode_ == 0)
         {
            psIOExpertMode_ = 1;
            printf("IO expert mode on\n");
         }
         else
         {
            psIOExpertMode_ = 0;
            printf("IO expert mode off\n");
         }
      }

      // +++ rs_expert 
      else if (!strcmp(command, "rs_expert"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("This command turns on/off response surface expert mode.\n");
            printf("This mode, when turned on, will trigger prompts for\n");
            printf("additional information regarding creation of response\n");
            printf("surfaces.\n");
            printf("When this mode is off, default setting will be used.\n");
            continue;
         }
         if (!strcmp(winput, "on"))  psRSExpertMode_ = 0;
         if (!strcmp(winput, "off")) psRSExpertMode_ = 1;
         if (psRSExpertMode_ == 0)
         {
            psRSExpertMode_ = 1;
            printf("response surface expert mode on\n");
         }
         else
         {
            psRSExpertMode_ = 0;
            printf("response surface expert mode off\n");
         }
      }

      // +++ rs_codegen 
      else if (!strcmp(command, "rs_codegen"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("This command turns on/off response surface code\n");
            printf("generation mode. This mode, when turned on, will\n");
            printf("generate a psuade_rs.info file whenever a response\n");
            printf("surface check ('rscheck') is performed. This file\n");
            printf("contains information for constructing the stand-alone\n");
            printf("response surface interpolators apart from PSUADE.\n");
            printf("In addition, a psuade_rs.py file will also be created\n");
            printf("as a stand-alone Python interpolator.\n");
            continue;
         }
         if (!strcmp(winput, "on"))  psRSCodeGen_ = 0;
         if (!strcmp(winput, "off")) psRSCodeGen_ = 2;
         if (psRSCodeGen_ == 0)
         {
            psRSCodeGen_ = 2;
            printf("response surface code generation mode on\n");
         }
         else
         {
            psRSCodeGen_ = 0;
            printf("response surface code generation mode off\n");
         }
      }

      // +++ ana_expert 
      else if (!strcmp(command, "ana_expert"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("This command turns on/off data analysis expert mode.\n");
            printf("This mode, when turned on, will trigger prompts for\n");
            printf("additional information regarding data analysis.\n");
            printf("When this mode is off, default setting will be used.\n");
            continue;
         }
         if (!strcmp(winput, "on"))  psAnaExpertMode_ = 0;
         if (!strcmp(winput, "off")) psAnaExpertMode_ = 1;
         if (psAnaExpertMode_ == 0)
         {
            psAnaExpertMode_ = 1;
            printf("analysis expert mode on\n");
         }
         else
         {
            psAnaExpertMode_ = 0;
            printf("analysis expert mode off\n");
         }
      }

      // +++ sam_expert 
      else if (!strcmp(command, "sam_expert"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("This command turns on/off sample design expert mode.\n");
            printf("This mode, when turned on, will trigger prompts for\n");
            printf("additional information regarding creation of samples.\n");
            printf("When this mode is off, default setting will be used.\n");
            continue;
         }
         if (!strcmp(winput, "on"))  psSamExpertMode_ = 0;
         if (!strcmp(winput, "off")) psSamExpertMode_ = 1;
         if (psSamExpertMode_ == 0)
         {
            psSamExpertMode_ = 1;
            printf("sampling expert mode on\n");
         }
         else
         {
            psSamExpertMode_ = 0;
            printf("sampling expert mode off\n");
         }
      }

      // +++ opt_expert 
      else if (!strcmp(command, "opt_expert"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("This command turns on/off optimization expert mode.\n");
            printf("This mode, when turned on, will trigger prompts for\n");
            printf("additonal information regarding numerical optimization.\n");
            printf("When this mode is off, default setting will be used.\n");
            continue;
         }
         if (!strcmp(winput, "on"))  psOptExpertMode_ = 0;
         if (!strcmp(winput, "off")) psOptExpertMode_ = 1;
         if (psOptExpertMode_ == 0)
         {
            psOptExpertMode_ = 1;
            printf("optimization expert mode on\n");
         }
         else
         {
            psOptExpertMode_ = 0;
            printf("optimization expert mode off\n");
         }
      }

      // +++ genhistogram 
      else if (!strcmp(command, "genhistogram"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("This command uses the loaded sample and generates\n");
            printf("a sample of scenarios each associated with a\n");
            printf("probability.\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data (load sample data first).\n");
            continue;
         }
         int *incrs = new int[nInputs_];
         printf("Please enter the number of bins per input dimension.\n");
         for (ii = 0; ii < nInputs_; ii++)
         {
            sprintf(pString, "Number of histogram bins for input %d : ", ii+1);
            incrs[ii] = getInt(2, 1000, pString);
         }
         PDFHistogram *pdfhist = new PDFHistogram(nSamples_, nInputs_, 
                                                  sampleInputs_, incrs, 1);
         delete pdfhist;
         delete [] incrs;
         printf("Histogram file has been created in psuade_pdfhist_sample.\n");
      }

      // +++ genhistogram2 
      else if (!strcmp(command, "genhistogram2"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("This command uses the loaded sample and generates\n");
            printf("a sample of scenarios each associated with a\n");
            printf("probability. The number of scenarios is obtained\n");
            printf("from the user.\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data (load sample data first).\n");
            continue;
         }
         int *incrs = new int[nInputs_];
         sprintf(pString, "Approximate number of scenarios (2-%d): ",nSamples_);
         int nScen = getInt(2, nSamples_, pString);
         int nbins = (int) pow(10.0 * nScen, 1.0/nInputs_);
         int kdiff = 1;
         int ksign = 0;
         PDFHistogram *pdfhist=NULL;
         while (1)
         {
            for (ii = 0; ii < nInputs_; ii++) incrs[ii] = nbins;
            pdfhist = new PDFHistogram(nSamples_,nInputs_,sampleInputs_,incrs,1);
            delete pdfhist;
            fp = fopen("psuade_pdfhist_sample", "r");
            if (fp == NULL)
            {
               printf("ERROR: psuade_pdfhist_sample file not found.\n");
               continue;
            }
            else
            {
               fscanf(fp, "%d", &kk);
               fclose(fp);
               printf("Current iteration: nbins = %d, scenario size = %d\n",
                      nbins,kk);
               kdiff = kk - nScen;
               if (kdiff == 0) break;
               if (kdiff < 0 && ksign == 1) break;
               if (kdiff > 0 && ksign == -1) break;
               if (kdiff < 0) nbins++;
               if (kdiff > 0) nbins--;
               if (kdiff < 0) ksign = -1;
               if (kdiff > 0) ksign = 1;
            }
         }
         delete [] incrs;
         printf("Histogram file has been created in psuade_pdfhist_sample.\n");
      }

      // +++ master 
      else if (!strcmp(command, "master"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("This command turns on/off master mode.\n");
            printf("This mode, when turned on, will trigger prompts\n");
            printf("for additional information in certain analysis.\n");
            continue;
         }
         if (!strcmp(winput, "on"))  psMasterMode_ = 0;
         if (!strcmp(winput, "off")) psMasterMode_ = 1;
         if (psMasterMode_ == 0)
         {
            psMasterMode_ = 1;
            printf("Master mode on\n");
         }
         else
         {
            psMasterMode_ = 0;
            printf("Master mode off\n");
         }
      }

      // +++ gm 
      else if (!strcmp(command, "gm"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("This command turns on/off gm mode.\n");
            printf("This mode, when turned on, will trigger prompts\n");
            printf("for additional information in certain analysis.\n");
            continue;
         }
         if (!strcmp(winput, "on"))  psGMMode_ = 0;
         if (!strcmp(winput, "off")) psGMMode_ = 1;
         if (psGMMode_ == 0)
         {
            psGMMode_ = 1;
            printf("gm mode on\n");
         }
         else
         {
            psGMMode_ = 0;
            printf("gm mode off\n");
         }
      }

      // +++ use_configfile 
      else if (!strcmp(command, "use_configfile"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("This command turns on the use_configFile mode.\n");
            printf("syntax: use_configfile <file>  or\n");
            printf("        use_configfile off     to turn off this mode.\n");
            printf("This mode, when used, will cause some functions\n");
            printf("to read their parameters from the configure file \n");
            printf("instead of prompting users for parameters.\n");
            printf("Note: use genconfigfile to see what are available.\n");
            continue;
         }
         if (!strcmp(winput, "off"))
         {
            printf("Turn off config file usage.\n");
            if (psConfig_ != NULL) delete psConfig_;
            psConfig_ = NULL;
         }
         fp = fopen(winput, "r");
         if (fp == NULL)
         {
            printf("ERROR: config file not found.\n");
            continue;
         }
         if (psConfigFileName_ == NULL) psConfigFileName_ = new char[500];
         strcpy(psConfigFileName_, winput);
         if (psConfig_ != NULL) delete psConfig_;
         psConfig_ = new PsuadeConfig(psConfigFileName_,1);
      }

      // +++ scilab 
      else if (!strcmp(command, "scilab"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("This command switches between matlab/scilab mode.\n");
            printf("When PSUADE creates output graphics, this mode will\n");
            printf("direct PSUADE to create matlab or scilab graphics.\n");
            printf("The default is matlab.\n");
            continue;
         }
         if (psPlotTool_ == 0)
         {
            psPlotTool_ = 1;
            printf("scilab mode on\n");
         }
         else
         {
            psPlotTool_ = 0;
            printf("scilab mode off\n");
         }
      }

      // +++ rsmax 
      else if (!strcmp(command, "rsmax"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rsmax: set maximum sample size for response surface\n"); 
            printf("       (used to guard against large sample sizes that\n");
            printf("        make response surface generation very expensive.)\n");
            printf("       The default max is 5000.\n");
            continue;
         }
         sscanf(lineIn, "%s %d", command, &psFAMaxDataPts_);
         if (psFAMaxDataPts_ <= 1 || psFAMaxDataPts_ > 100000)
            psFAMaxDataPts_ = 100000;
         printf("psuade max data size for response surface set to %d\n", 
                psFAMaxDataPts_);
      }

      // +++ rsca 
      else if (!strcmp(command, "rsca"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rsca: RS-based brute-force calibration (use mesh data)\n");
            printf("When the input dimension is low, an alternative to MCMC\n");
            printf("is to do a brute force search on an interpolated mesh.\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data (load data first).\n");
            continue;
         }
         if (nInputs_ > 4)
         {
            printf("INFO: nInputs_ must be <= 4.\n");
            continue;
         }
         if (psPlotTool_ != 0)
         {
            printf("INFO: rsca is currently not available for scilab.\n");
            continue;
         }
         if (faPtr != NULL) delete faPtr;

         inputSettings = new double[nInputs_];
         iplot1 = 0;
         iplot2 = 1;
         iplot3 = 2;
         kk = nInputs_;
         if (nInputs_ > 3)
         {
            printf("Currently, you can only select at most 3 inputs.\n");
            sprintf(pString, "Enter the number of inputs to study (1 - 3) : ");
            kk = getInt(1, 3, pString);
         }
         if (kk == 1) iplot2 = iplot3 = -2;
         if (kk == 2) iplot3 = -2;
         if (nInputs_ > kk)
         {
            sprintf(pString, "Enter the first input (1 - %d) : ", nInputs_);
            iplot1 = getInt(1, nInputs_, pString);
            iplot1--;
            if (kk > 1)
            {
               iplot2 = iplot1;
               while (iplot1 == iplot2)
               {
                  sprintf(pString, "Enter the second input (1 - %d), not %d : ",
                          nInputs_, iplot1+1);
                  iplot2 = getInt(1, nInputs_, pString);
                  iplot2--;
                  if (iplot1 == iplot2)
                     printf("ERROR: duplicate input number %d.\n",iplot2+1);
               }
            }
            if (kk > 2)
            {
               iplot3 = 3 - iplot1 - iplot2;
               while (iplot3 < 0 || iplot3 == iplot1 || iplot3 == iplot2)
               {
                  sprintf(pString,
                    "Enter the input for t axis (1 - %d), not %d nor %d: ",
                          nInputs_, iplot1+1, iplot2+1);
                  iplot3 = getInt(1, nInputs_, pString);
                  iplot3--;
                  if (iplot3 == iplot1 || iplot3 == iplot2)
                     printf("ERROR: duplicate input number %d.\n",iplot3+1);
               }
            }
         }

         if      (kk == 1) nPtsPerDim = 1024;
         else if (kk == 2) nPtsPerDim = 64;
         else if (kk == 3) nPtsPerDim = 32;
         faFlag = 3;
         faPtr = genFAInteractive(psuadeIO_, faFlag);
         if (faPtr == NULL) {printf("Error detected.\n"); continue;}
         faPtr->setOutputLevel(outputLevel_);
         faPtr->setNPtsPerDim(nPtsPerDim);
         faPtr->setBounds(iLowerB_, iUpperB_);
         if (nInputs_ > kk)
            printf("The other inputs will be set to their nominal values.\n");
         for (iInd1 = 0; iInd1 < nInputs_; iInd1++)
         {
            if (iInd1 != iplot1 && iInd1 != iplot2 && iInd1 != iplot3)
               inputSettings[iInd1] = 0.5*(iLowerB_[iInd1]+iUpperB_[iInd1]);
            else inputSettings[iInd1] = 1.0;
         }
         jplot = 0;
         faYIn = new double[nSamples_];
         for (sInd = 0; sInd < nSamples_; sInd++)
            faYIn[sInd] = sampleOutputs_[sInd*nOutputs_+jplot];
         if (kk == 3)
         {
            printf("Please wait while generating data \n");
            fp = fopen("matlabrsbca3.m", "w");
            if (fp == NULL)
            {
               printf("ERROR: cannot open file matlabrsbca3.m.\n");
               continue;
            }
            fprintf(fp, "twoPlots = 1;\n");
            fwritePlotCLF(fp);
            GYmax = -1.0e35;
            GYmin =  1.0e35;
            for (ii = 0; ii < nPtsPerDim; ii++)
            {
               printf(".");
               fflush(stdout);
               inputSettings[iplot3] = (iUpperB_[iplot3] - iLowerB_[iplot3]) *
                                   ii / (nPtsPerDim - 1.0) + iLowerB_[iplot3];
               faPtr->gen2DGridData(sampleInputs_,faYIn, iplot1, iplot2,
                                   inputSettings, &faLeng, &faXOut,&faYOut);
               if (ii == 0)
               {
                  fprintf(fp, "X = [\n");
                  for (sInd = 0; sInd < faLeng; sInd+=nPtsPerDim)
                     fprintf(fp, "%e ", faXOut[sInd*2]);
                  fprintf(fp, "];\n");
                  fprintf(fp, "Y(:,:) = [\n");
                  for (sInd = 0; sInd < nPtsPerDim; sInd++)
                     fprintf(fp, "%e ", faXOut[sInd*2+1]);
                  fprintf(fp, "];\n");
               }
               for (sInd = 0; sInd < faLeng; sInd++)
                  faYOut[sInd] = exp(-faYOut[sInd]);

               fprintf(fp, "A%d = [\n", ii + 1);
               for (sInd = 0; sInd < nPtsPerDim; sInd++)
               {
                  for (jj = 0; jj < nPtsPerDim; jj++)
                  {
                     dtemp = faYOut[sInd*nPtsPerDim+jj];
                     count = 1;
                     if (sInd > 0)
                     {
                        dtemp += faYOut[(sInd-1)*nPtsPerDim+jj];
                        count++;
                     }
                     if (sInd < nPtsPerDim-1)
                     {
                        dtemp += faYOut[(sInd+1)*nPtsPerDim+jj];
                        count++;
                     }
                     if (jj > 0)
                     {
                        dtemp += faYOut[sInd*nPtsPerDim+jj-1];
                        count++;
                     }
                     if (jj < nPtsPerDim-1)
                     {
                        dtemp += faYOut[sInd*nPtsPerDim+jj+1];
                        count++;
                     }
                     if (sInd > 0 && jj > 0)
                     {
                        dtemp += faYOut[(sInd-1)*nPtsPerDim+jj-1];
                        count++;
                     }
                     if (sInd > 0 && jj < nPtsPerDim-1)
                     {
                        dtemp += faYOut[(sInd-1)*nPtsPerDim+jj+1];
                        count++;
                     }
                     if (sInd < (nPtsPerDim-1) && jj > 0)
                     {
                        dtemp += faYOut[(sInd+1)*nPtsPerDim+jj-1];
                        count++;
                     }
                     if (sInd < (nPtsPerDim-1) && jj < (nPtsPerDim-1))
                     {
                        dtemp += faYOut[(sInd+1)*nPtsPerDim+jj+1];
                        count++;
                     }
                     dtemp /= ((double) count);
                     fprintf(fp, "%e\n", dtemp);
                     if (dtemp > Ymax) Ymax = dtemp;
                     if (dtemp < Ymin) Ymin = dtemp;
                  }
               }
               fprintf(fp, "];\n");
               fprintf(fp, "A%d = reshape(A%d,%d,%d);\n", ii+1, ii+1,
                       nPtsPerDim, nPtsPerDim);
                                                                                
               if (Ymax > GYmax) GYmax = Ymax;
               if (Ymin < GYmin) GYmin = Ymin;
               delete [] faXOut;
               delete [] faYOut;
            }

            fprintf(fp, "disp(\'Please wait while loading data.\')\n");
            fprintf(fp, "hold off\n");
            fwritePlotCLF(fp);
            for (ii = 0; ii < nPtsPerDim; ii++)
            {
               dtemp = (iUpperB_[iplot3] - iLowerB_[iplot3]) *
                       ii / (nPtsPerDim - 1.0) + iLowerB_[iplot3];
               fprintf(fp, "disp(\'Plotting frame %d of %d\')\n",
                       ii+1,nPtsPerDim);
               fprintf(fp, "if twoPlots == 1\n");
               fprintf(fp, "subplot(1,2,1), mesh(X,Y,A%d)\n", ii+1);
               fprintf(fp, "axis([%e %e %e %e %e %e])\n",iLowerB_[iplot1],
                       iUpperB_[iplot1],iLowerB_[iplot2],iUpperB_[iplot2],
                       GYmin, GYmax);
               fwritePlotXLabel(fp, inputNames_[iplot1]);
               fwritePlotYLabel(fp, inputNames_[iplot2]);
               fwritePlotZLabel(fp, outputNames_[jplot]);
               fwritePlotAxes(fp);
               fprintf(fp, "colorbar\n");
               fprintf(fp, "title(\'%s Mesh plot, val(3) = %14.7e\')\n",
                       outputNames_[jplot], dtemp);
               fprintf(fp, "subplot(1,2,2)\n");
               fprintf(fp, "end\n");
               fprintf(fp, "contourf(X,Y,A%d)\n",ii+1);
               fprintf(fp, "axis([%e %e %e %e])\n",iLowerB_[iplot1],
                       iUpperB_[iplot1],iLowerB_[iplot2],iUpperB_[iplot2]);
               fwritePlotXLabel(fp, inputNames_[iplot1]);
               fwritePlotYLabel(fp, inputNames_[iplot2]);
               fwritePlotAxes(fp);
               fprintf(fp, "colorbar\n");
               fprintf(fp, "colormap(jet)\n");
               fprintf(fp, "title(\'%s contour plot, val(3) = %14.7e\')\n",
                       outputNames_[jplot], dtemp);
               fprintf(fp,"pause(1)\n");
            }
            fprintf(fp, "rotate3d on\n");
            fclose(fp);
            printf("\nmatlabrsbca3.m is now available.\n");
            printf("You can identify the max and min from the plots.\n");
            delete [] faYIn;
            delete [] inputSettings;
            delete faPtr;
            faPtr = NULL;
         }

         if (kk == 2)
         {
            printf("Please wait while generating data ");
            fp = fopen("matlabrsbca2.m", "w");
            if (fp == NULL)
            {
               printf("ERROR: cannot open file matlabrsbca2.m.\n");
               continue;
            }
            fprintf(fp, "twoPlots = 1;\n");
            fwritePlotCLF(fp);
            faPtr->gen2DGridData(sampleInputs_,faYIn,iplot1,iplot2,
                          inputSettings, &faLeng, &faXOut,&faYOut);
                                                                                
            for (jj = 0; jj < faLeng; jj++) faYOut[jj] = exp(-faYOut[jj]);
            Ymax = -1.0e35;
            Ymin =  1.0e35;
            fprintf(fp, "A = [\n");
            for (sInd = 0; sInd < nPtsPerDim; sInd++)
            {
               for (jj = 0; jj < nPtsPerDim; jj++)
               {
                  dtemp = faYOut[sInd*nPtsPerDim+jj];
                  count = 1;
                  if (sInd > 0)
                  {
                     dtemp += faYOut[(sInd-1)*nPtsPerDim+jj];
                     count++;
                  }
                  if (sInd < nPtsPerDim-1)
                  {
                     dtemp += faYOut[(sInd+1)*nPtsPerDim+jj];
                     count++;
                  }
                  if (jj > 0)
                  {
                     dtemp += faYOut[sInd*nPtsPerDim+jj-1];
                     count++;
                  }
                  if (jj < nPtsPerDim-1)
                  {
                     dtemp += faYOut[sInd*nPtsPerDim+jj+1];
                     count++;
                  }
                  if (sInd > 0 && jj > 0)
                  {
                     dtemp += faYOut[(sInd-1)*nPtsPerDim+jj-1];
                     count++;
                  }
                  if (sInd > 0 && jj < nPtsPerDim-1)
                  {
                     dtemp += faYOut[(sInd-1)*nPtsPerDim+jj+1];
                     count++;
                  }
                  if (sInd < (nPtsPerDim-1) && jj > 0)
                  {
                     dtemp += faYOut[(sInd+1)*nPtsPerDim+jj-1];
                     count++;
                  }
                  if (sInd < (nPtsPerDim-1) && jj < (nPtsPerDim-1))
                  {
                     dtemp += faYOut[(sInd+1)*nPtsPerDim+jj+1];
                     count++;
                  }
                  dtemp /= ((double) count);
                  fprintf(fp, "%e\n", dtemp);
                  if (dtemp > Ymax) Ymax = dtemp;
                  if (dtemp < Ymin) Ymin = dtemp;
               }
            }
            fprintf(fp, "];\n");
            fprintf(fp, "A = reshape(A,%d,%d);\n", nPtsPerDim, nPtsPerDim);
            fprintf(fp, "X = [\n");
            for (sInd = 0; sInd < faLeng; sInd+=nPtsPerDim)
               fprintf(fp, "%e\n", faXOut[sInd*2]);
            fprintf(fp, "];\n");
            fprintf(fp, "Y = [\n");
            for (sInd = 0; sInd < nPtsPerDim; sInd++)
               fprintf(fp, "%e\n", faXOut[sInd*2+1]);
            fprintf(fp, "];\n");
            fwritePlotCLF(fp);
            fprintf(fp, "if twoPlots == 1\n");
            fprintf(fp, "subplot(1,2,1), mesh(X,Y,A)\n");
            fwritePlotXLabel(fp, inputNames_[iplot1]);
            fwritePlotYLabel(fp, inputNames_[iplot2]);
            fwritePlotZLabel(fp, outputNames_[jplot]);
            fwritePlotAxes(fp);
            fprintf(fp, "title(\'Mesh Plot for %s\')\n",outputNames_[jplot]);
            fprintf(fp, "colorbar\n");
            fprintf(fp, "subplot(1,2,2)\n");
            fprintf(fp, "end\n");
            fprintf(fp, "contourf(X,Y,A)\n");
            fwritePlotXLabel(fp, inputNames_[iplot1]);
            fwritePlotYLabel(fp, inputNames_[iplot2]);
            fwritePlotAxes(fp);
            fprintf(fp, "colorbar\n");
            fprintf(fp, "colormap(jet)\n");
            fprintf(fp, "title(\'Contour Plot for %s\')\n",outputNames_[jplot]);
            fclose(fp);
            printf("matlabrsbca2.m is now available.\n");
            printf("You can identify the max and min from the plots.\n");
            delete [] faXOut;
            delete [] faYIn;
            delete [] faYOut;
            delete [] inputSettings;
            delete faPtr;
            faPtr = NULL;
         }

         if (kk == 1)
         {
            printf("Please wait while generating data \n");
            fp = fopen("matlabrsbca1.m", "w");
            if (fp == NULL)
            {
               printf("ERROR: cannot open file matlabrsbca1.m.\n");
               continue;
            }
            fwritePlotCLF(fp);
            faPtr->gen1DGridData(sampleInputs_,faYIn, iplot1,
                                 inputSettings, &faLeng, &faXOut,&faYOut);
            for (jj = 0; jj < faLeng; jj++)
               faYOut[jj] = exp(-faYOut[jj]);
            Ymax = -1.0e35;
            for (sInd = 0; sInd < faLeng; sInd++)
            {
               if (PABS(faYOut[sInd]) > Ymax) Ymax = PABS(faYOut[sInd]);
               if (PABS(faYOut[sInd]) < Ymin) Ymin = PABS(faYOut[sInd]);
            }
            if (Ymax == 0.0) Ymax = 1.0;
            printf("Ymax = %e\n", Ymax);
            for (jj = 0; jj < faLeng; jj++) faYOut[jj] /= Ymax;
            fprintf(fp, "A = [\n");
            tempV = new double[faLeng];
            sprintf(pString,"How many smoothing step : (0 - 100) : ");
            nSmooth = getInt(0, 1000, pString);
            for (jj = 0; jj < nSmooth; jj++)
            {
               for (sInd = 0; sInd < faLeng; sInd++)
                  tempV[sInd] = faYOut[sInd];
     
               for (sInd = 0; sInd < faLeng; sInd++)
               {
                  dtemp = tempV[sInd];
                  count = 1;
                  if (sInd > 0)
                  {
                     dtemp += tempV[sInd-1];
                     count++;
                  }
                  if (sInd < nPtsPerDim-1)
                  {
                     dtemp += tempV[sInd+1];
                     count++;
                  }
                  dtemp /= ((double) count);
                  faYOut[sInd] = dtemp;
               }
            }
            Ymax = -1.0e35;
            Ymin =  1.0e35;
            for (sInd = 0; sInd < faLeng; sInd++)
            {
               fprintf(fp, "%e\n", faYOut[sInd]);
               if (dtemp > Ymax) Ymax = dtemp;
               if (dtemp < Ymin) Ymin = dtemp;
            }
            fprintf(fp, "];\n");
            fprintf(fp, "X = [\n");
            for (sInd = 0; sInd < faLeng; sInd++)
               fprintf(fp, "%e\n", faXOut[sInd]);
            fprintf(fp, "];\n");
            fprintf(fp, "plot(X,A)\n");
            fwritePlotXLabel(fp, inputNames_[iplot1]);
            fwritePlotYLabel(fp, outputNames_[jplot]);
            fwritePlotAxes(fp);
            fprintf(fp, "title(\'Likelihood Plot for %s\')\n",
                    outputNames_[jplot]);
            fclose(fp);
            printf("matlabrsbca1.m is now available.\n");
            printf("You can identify the max and min from the plots.\n");
            delete [] faXOut;
            delete [] faYIn;
            delete [] faYOut;
            delete [] tempV;
            delete [] inputSettings;
            delete faPtr;
            faPtr = NULL;
         }
      }

      // +++ pdfcheck 
      else if (!strcmp(command, "pdfcheck"))
      {
         int    nSam, *PDFs;
         double *sOutputs, *slbounds, *subounds, *smeans, *sstdevs;
         Vector vIn, vOut;

         PDFs     = new int[4];
         smeans   = new double[4];
         sstdevs  = new double[4];
         slbounds = new double[4];
         subounds = new double[4];
         nSam = 100000;
         sOutputs = new double[nSam];

         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("pdfcheck: 1-sample and 2-sample tests for distributions.\n");
            printf("syntax: pdfcheck (no argument needed)\n");
            printf("This is a demonstration command only. It creates and\n");
            printf("reconstruct different probability distributions.\n");
            continue;
         }
         printAsterisks(PL_INFO, 0);
         printAsterisks(PL_INFO, 0);
         printEquals(PL_INFO, 0);
         printf("##### One-input test (Normal(1,4)): \n");
         printf("Sample mean     should be = 1\n");
         printf("Sample std dev  should be = 2\n");
         printf("Sample skewness should be = 0\n");
         printf("Sample Kurtosis should be = 3\n");
         printDashes(PL_INFO, 0);
         smeans[0] = 1.0;
         sstdevs[0] = 2.0;
         slbounds[0] = -6.0;
         subounds[0] = 8.0;
         pdfman = new PDFManager();
         PDFs[0] = PSUADE_PDF_NORMAL;
         ddata = 1.0;
         corMat.setDim(1,1);
         corMat.setEntry(0,0,ddata);
         pdfman->initialize(1, PDFs, smeans, sstdevs, corMat, NULL, NULL);
         vecLower.load(1, slbounds);
         vecUpper.load(1, subounds);
         vIn.setLength(nSam); 
         vOut.setLength(nSam); 
         pdfman->genSample(nSam, vOut, vecLower, vecUpper);
         analysisMethod = PSUADE_ANA_MOMENT;
         anaManager = new AnalysisManager();
         anaManager->setup(analysisMethod, 0);
         anaManager->analyze(analysisMethod,nSam,vecLower,vecUpper,vIn,vOut,
                             outputLevel_);
         delete anaManager;
         delete pdfman;

         pdfman = new PDFManager();
         vOut.setLength(nSam); 
         printAsterisks(PL_INFO, 0);
         printAsterisks(PL_INFO, 0);
         printEquals(PL_INFO, 0);
         printf("##### One-input test (LogNormal(1,2)): \n");
         printf("Sample mean     should be = 4.5\n");
         printf("Sample std dev  should be = 5.8\n");
         printDashes(PL_INFO, 0);
         smeans[0] = 1.0;
         sstdevs[0] = 1.0;
         slbounds[0] = 0.0;
         subounds[0] = 200.0;
         PDFs[0] = PSUADE_PDF_LOGNORMAL;
         ddata = 1.0;
         corMat.setDim(1,1);
         corMat.setEntry(0,0,ddata);
         pdfman->initialize(1, PDFs, smeans, sstdevs, corMat, NULL, NULL);
         vecLower.load(1, slbounds);
         vecUpper.load(1, subounds);
         vOut.setLength(nSam); 
         pdfman->genSample(nSam, vOut, vecLower, vecUpper);
         analysisMethod = PSUADE_ANA_MOMENT;
         anaManager = new AnalysisManager();
         anaManager->setup(analysisMethod, 0);
         vIn.setLength(nSam); 
         anaManager->analyze(analysisMethod,nSam,vecLower,vecUpper,vIn,vOut,
                             outputLevel_);
         delete anaManager;
         delete pdfman;

         pdfman = new PDFManager();
         vOut.setLength(nSam); 
         printAsterisks(PL_INFO, 0);
         printAsterisks(PL_INFO, 0);
         printEquals(PL_INFO, 0);
         printf("##### One-input test (Weibull(1,1)): \n");
         printf("Sample mean     should be = 1\n");
         printf("Sample std dev  should be = 1\n");
         printf("Sample skewness should be = 2\n");
         printDashes(PL_INFO, 0);
         smeans[0] = 1.0;
         sstdevs[0] = 1.0;
         slbounds[0] = 0.0;
         subounds[0] = 10.0;
         PDFs[0] = PSUADE_PDF_WEIBULL;
         ddata = 1.0;
         corMat.setDim(1,1);
         corMat.setEntry(0,0,ddata);
         pdfman->initialize(1, PDFs, smeans, sstdevs, corMat, NULL, NULL);
         vecLower.load(1, slbounds);
         vecUpper.load(1, subounds);
         vOut.setLength(nSam); 
         pdfman->genSample(nSam, vOut, vecLower, vecUpper);
         analysisMethod = PSUADE_ANA_MOMENT;
         anaManager = new AnalysisManager();
         anaManager->setup(analysisMethod, 0);
         vIn.setLength(nSam); 
         anaManager->analyze(analysisMethod,nSam,vecLower,vecUpper,vIn,vOut,
                             outputLevel_);
         delete anaManager;
         delete pdfman;

         printAsterisks(PL_INFO, 0);
         printAsterisks(PL_INFO, 0);
         printEquals(PL_INFO, 0);
         printf("##### Two-input test (Normal(0,1)): \n");
         printf("Sample mean     should be = 0\n");
         printf("Sample std dev  should be = 1.414\n");
         printf("Sample skewness should be = 0\n");
         printf("Sample Kurtosis should be = 3\n");
         printDashes(PL_INFO, 0);
         pdfman = new PDFManager();
         PDFs[0] = PSUADE_PDF_NORMAL;
         PDFs[1] = PSUADE_PDF_NORMAL;
         smeans[0] = smeans[1] = 0.0;
         sstdevs[0] = sstdevs[1] = 1.0;
         slbounds[0] = slbounds[1] = -3.0;
         subounds[0] = subounds[1] =  3.0;
         ddata = 1.0;
         corMat.setDim(2,2);
         corMat.setEntry(0,0,ddata);
         corMat.setEntry(1,1,ddata);
         pdfman->initialize(2, PDFs, smeans, sstdevs, corMat, NULL, NULL);
         vecLower.load(2, slbounds);
         vecUpper.load(2, subounds);
         vIn.setLength(2*nSam); 
         pdfman->genSample(nSam, vIn, vecLower, vecUpper);
         vOut.setLength(nSam); 
         for (ii = 0; ii < nSam; ii++) vOut[ii] = vIn[ii*2] + vIn[ii*2+1];
         analysisMethod = PSUADE_ANA_MOMENT;
         anaManager = new AnalysisManager();
         anaManager->setup(analysisMethod, 0);
         anaManager->analyze(analysisMethod,nSam,vecLower,vecUpper,vIn,vOut,
                             outputLevel_);
         delete anaManager;
         delete pdfman;

         pdfman = new PDFManager();
         vOut.setLength(nSam); 
         printAsterisks(PL_INFO, 0);
         printAsterisks(PL_INFO, 0);
         printEquals(PL_INFO, 0);
         printf("##### Two-input test (LogNormal(1,1)): \n");
         printf("Sample mean     should be = 8.9634\n");
         printf("Sample std dev  should be = 8.3081\n");
         printDashes(PL_INFO, 0);
         smeans[0] = smeans[1] = 1.0;
         sstdevs[0] = sstdevs[1] = 1.0;
         slbounds[0] = slbounds[1] =  0.0;
         subounds[0] = subounds[1] =  200.0;
         PDFs[0] = PSUADE_PDF_LOGNORMAL;
         PDFs[1] = PSUADE_PDF_LOGNORMAL;
         ddata = 1.0;
         corMat.setDim(2,2);
         corMat.setEntry(0,0,ddata);
         corMat.setEntry(1,1,ddata);
         pdfman->initialize(2, PDFs, smeans, sstdevs, corMat, NULL, NULL);
         vecLower.load(2, slbounds);
         vecUpper.load(2, subounds);
         vIn.setLength(2*nSam); 
         pdfman->genSample(nSam, vIn, vecLower, vecUpper);
         vOut.setLength(nSam); 
         for (ii = 0; ii < nSam; ii++) vOut[ii] = vIn[ii*2] + vIn[ii*2+1];
         analysisMethod = PSUADE_ANA_MOMENT;
         anaManager = new AnalysisManager();
         anaManager->setup(analysisMethod, 0);
         anaManager->analyze(analysisMethod,nSam,vecLower,vecUpper,vIn,vOut,
                             outputLevel_);
         delete anaManager;
         delete pdfman;

         printAsterisks(PL_INFO, 0);
         printAsterisks(PL_INFO, 0);
         printEquals(PL_INFO, 0);
         printf("##### Two-input test (Normal(0,1)) with cor = 0.5: \n");
         printf("Sample mean     should be = 0\n");
         printf("Sample std dev  should be = 1.73\n");
         printf("Sample skewness should be = 0\n");
         printf("Sample Kurtosis should be = 3\n");
         printDashes(PL_INFO, 0);
         pdfman = new PDFManager();
         PDFs[0] = PSUADE_PDF_NORMAL;
         PDFs[1] = PSUADE_PDF_NORMAL;
         smeans[0] = smeans[1] = 0.0;
         sstdevs[0] = sstdevs[1] = 1.0;
         slbounds[0] = slbounds[1] = -5.0;
         subounds[0] = subounds[1] =  5.0;
         ddata = 1.0;
         corMat.setDim(2,2);
         corMat.setEntry(0,0,ddata);
         corMat.setEntry(1,1,ddata);
         ddata = 0.5;
         corMat.setEntry(0,1,ddata);
         corMat.setEntry(1,0,ddata);
         pdfman->initialize(2, PDFs, smeans, sstdevs, corMat, NULL, NULL);
         vecLower.load(2, slbounds);
         vecUpper.load(2, subounds);
         vIn.setLength(2*nSam); 
         pdfman->genSample(nSam, vIn, vecLower, vecUpper);
         vOut.setLength(nSam); 
         for (ii = 0; ii < nSam; ii++) vOut[ii] = vIn[ii*2] + vIn[ii*2+1];
         analysisMethod = PSUADE_ANA_MOMENT;
         anaManager = new AnalysisManager();
         anaManager->setup(analysisMethod, 0);
         anaManager->analyze(analysisMethod,nSam,vecLower,vecUpper, vIn, vOut,
                             outputLevel_);
         delete anaManager;
         delete pdfman;

         pdfman = new PDFManager();
         vOut.setLength(nSam); 
         printAsterisks(PL_INFO, 0);
         printAsterisks(PL_INFO, 0);
         printEquals(PL_INFO, 0);
         printf("##### Two-input test (LogNormal(1,1)) with cor = 0.5: \n");
         printf("Sample mean     should be = 7.4\n");
         printf("Sample mean     should be = 1.0\n");
         printDashes(PL_INFO, 0);
         smeans[0] = smeans[1] = 1.0;
         sstdevs[0] = sstdevs[1] = 1.0;
         slbounds[0] = slbounds[1] =  0.0;
         subounds[0] = subounds[1] =  100.0;
         PDFs[0] = PSUADE_PDF_LOGNORMAL;
         PDFs[1] = PSUADE_PDF_LOGNORMAL;
         ddata = 1.0;
         corMat.setDim(2,2);
         corMat.setEntry(0,0,ddata);
         corMat.setEntry(1,1,ddata);
         ddata = 0.5;
         corMat.setEntry(0,1,ddata);
         corMat.setEntry(1,0,ddata);
         pdfman->initialize(2, PDFs, smeans, sstdevs, corMat, NULL, NULL);
         vecLower.load(2, slbounds);
         vecUpper.load(2, subounds);
         vIn.setLength(2*nSam); 
         pdfman->genSample(nSam, vIn, vecLower, vecUpper);
         vOut.setLength(nSam); 
         for (ii = 0; ii < nSam; ii++) vOut[ii] = vIn[ii*2] + vIn[ii*2+1];
         analysisMethod = PSUADE_ANA_MOMENT;
         anaManager = new AnalysisManager();
         anaManager->setup(analysisMethod, 0);
         anaManager->analyze(analysisMethod,nSam,vecLower,vecUpper,vIn,vOut,
                             outputLevel_);
         delete anaManager;
         delete pdfman;

         printEquals(PL_INFO, 0);
         printEquals(PL_INFO, 0);
         printf("4-input test (Normal(0,1), cor = 0.5, LogNormal(1,1), no cor: \n");
         printf("Sample mean     should be = 8.96\n");
         printf("Sample std dev  should be = 8.54\n");
         printDashes(PL_INFO, 0);
         pdfman = new PDFManager();
         PDFs[0] = PSUADE_PDF_NORMAL;
         PDFs[1] = PSUADE_PDF_NORMAL;
         PDFs[2] = PSUADE_PDF_LOGNORMAL;
         PDFs[3] = PSUADE_PDF_LOGNORMAL;
         smeans[0] = smeans[1] = 0.0;
         smeans[2] = smeans[3] = 1.0;
         sstdevs[0] = sstdevs[1] = sstdevs[2] = sstdevs[3] = 1.0;
         slbounds[0] = slbounds[1] = -10.0;
         slbounds[2] = slbounds[3] =  0.0;
         subounds[0] = subounds[1] =  10.0;
         subounds[2] = subounds[3] =  100.0;
         corMat.setDim(4,4);
         ddata = 1.0;
         corMat.setEntry(0,0,ddata);
         corMat.setEntry(1,1,ddata);
         corMat.setEntry(2,2,ddata);
         corMat.setEntry(3,3,ddata);
         ddata = 0.5;
         corMat.setEntry(0,1,ddata);
         corMat.setEntry(1,0,ddata);
         pdfman->initialize(4, PDFs, smeans, sstdevs, corMat, NULL, NULL);
         vecLower.load(4, slbounds);
         vecUpper.load(4, subounds);
         vIn.setLength(4*nSam); 
         pdfman->genSample(nSam, vIn, vecLower, vecUpper);
         vOut.setLength(nSam); 
         for (ii = 0; ii < nSam; ii++)
            vOut[ii] = vIn[ii*4] + vIn[ii*4+1] + vIn[ii*4+2] + vIn[ii*4+3];
         analysisMethod = PSUADE_ANA_MOMENT;
         anaManager = new AnalysisManager();
         anaManager->setup(analysisMethod, 0);
         anaManager->analyze(analysisMethod,nSam,vecLower,vecUpper, vIn, vOut,
                             outputLevel_);
         delete anaManager;
         delete pdfman;

         delete [] slbounds;
         delete [] subounds;
         delete [] sOutputs;
         delete [] smeans;
         delete [] sstdevs;
         delete [] PDFs;
      }

      // +++ gen_discrete 
      else if (!strcmp(command, "gen_discrete"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("gen_discrete: generate a sample of discrete variables.\n");
            printf("syntax: gen_discrete (no argument needed)\n");
            printf("You will be prompted for a parameter file describing\n");
            printf("the variables, their ranges, and their probabilities.\n");
            continue;
         }
         sampPtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_DISCRETE);
         sampPtr->setPrintLevel(outputLevel_);
         sampPtr->setOutputParams(iOne);
         sampPtr->setSamplingParams(iOne, -1, 0);
         sampPtr->initialize(0);
         count = sampPtr->getNumSamples();
         kk    = sampPtr->getNumInputs();
         tempX  = new double[count * kk];
         tempY  = new double[count * iOne];
         states = new int[count];
         sampPtr->getSamples(count,kk,iOne,tempX,tempY,states);
         delete sampPtr;
         sampPtr = NULL;
         tempT = new double[kk];
         tempW = new double[kk];
         for (ii = 0; ii < kk; ii++) tempT[ii] = PSUADE_UNDEFINED;
         for (ii = 0; ii < kk; ii++) tempW[ii] = -PSUADE_UNDEFINED;
         for (ii = 0; ii < kk; ii++)
         {
            for (jj = 0; jj < count; jj++)
            {
               if (tempX[jj*kk+ii] < tempT[ii]) tempT[ii] = tempX[jj*kk+ii];
               if (tempX[jj*kk+ii] > tempW[ii]) tempW[ii] = tempX[jj*kk+ii];
            }
         }
         names = new char*[kk];
         for (ii = 0; ii < kk; ii++)
         {
            names[ii] = new char[100];
            sprintf(names[ii], "X%d", ii+1);
         }
         strptr = new char*[1];
         strptr[0] = new char[100];
         strcpy(strptr[0], "Y");
         int    *iPDFs = new int[nInputs_];
         double *iMeans = new double[nInputs_];
         double *iStds = new double[nInputs_];
         Matrix *iCMat = new Matrix();
         iCMat->setDim(nInputs_, nInputs_);
         for (ii = 0; ii < nInputs_; ii++)
         {
            iPDFs[ii] = 0;
            iMeans[ii] = 0;
            iStds[ii] = 0;
            iCMat->setEntry(ii,ii,1.0);
         }
         ioPtr = new PsuadeData();
         ioPtr->updateInputSection(count, kk, NULL, tempT, tempW, tempX, 
                                   names, iPDFs, iMeans, iStds, iCMat); 
         delete [] iPDFs;
         delete [] iMeans;
         delete [] iStds;
         delete iCMat;
         ioPtr->updateOutputSection(count, iOne, tempY, states, strptr); 
         ioPtr->updateMethodSection(PSUADE_SAMP_MC, count, 1, -1, -1);
         ioPtr->writePsuadeFile("psuade_discrete_sample",0);
         printf("The test sample is in file psuade_discrete_sample.\n"); 
         delete [] tempX;
         delete [] tempY;
         delete [] tempT;
         delete [] tempW;
         delete [] states;
         delete [] strptr[0];
         delete [] strptr;
         delete ioPtr;
         ioPtr = NULL;
         tempX = tempY = tempT = tempW = NULL;
      }

      // +++ mo_opt
      else if (!strcmp(command, "mo_opt"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("mo_opt: multi-objective optimization\n");
            printf("syntax: mo_opt\n");
            printf("Note: no need to load data file first\n");
            continue;
         }
         sprintf(pString,"Enter the name of a PSUADE input file: ");
         getString(pString, dataFile);
         kk = strlen(dataFile);
         dataFile[kk-1] = '\0';
         printf("PSUADE input file = %s\n", dataFile);
         fp = fopen(dataFile, "r");
         if (fp == NULL)
         {
            printf("ERROR: file %s not found.\n", dataFile);
            continue;
         }
         fclose(fp);
         cleanUp();
         ioPtr = new PsuadeData();
         status = ioPtr->readPsuadeFile(dataFile);
         if (status != 0) 
         {
            printf("ERROR: cannot read file %s.\n", dataFile);
            continue;
         }
         ioPtr->getParameter("input_ninputs", pPtr);
         nInputs_ = pPtr.intData_;
         ioPtr->getParameter("output_noutputs", pPtr);
         nOutputs_ = pPtr.intData_;
         ioPtr->getParameter("input_lbounds", pLower);
         double *iLowerB = pLower.dbleArray_;
         ioPtr->getParameter("input_ubounds", pUpper);
         double *iUpperB = pUpper.dbleArray_;
         printf("Enter the name of the optimization driver. It can\n");
         printf("be an executable or a PSUADE data file to be used\n");
         printf("to create a response surface. In the latter case,\n");
         printf("make sure the number of inputs and outputs are the\n");
         printf("same as the PSUADE input file (%s).\n", dataFile);
         sprintf(pString,"Enter name of optimization driver file : ");
         getString(pString, winput);
         kk = strlen(winput);
         winput[kk-1] = '\0';
         printf("optimization driver (RS data) file = %s\n", winput);
         fp = fopen(winput, "r");
         if (fp == NULL)
         {
            printf("ERROR: file %s not found.\n", winput);
            continue;
         }
         ioPtr->updateApplicationSection(winput,winput,NULL,NULL,NULL,-1);
         ioPtr->updateOptimizationSection(10, 1, 1.0e-6);

         int    optCnt=0;
         double *optData[4]; 
         double *optInitX, *optInitY, *optX, *optY;
         optInitX = new double[nInputs_];
         for (ii = 0; ii < nInputs_; ii++)
            optInitX[ii] = 0.5 * (iLowerB[ii] + iUpperB[ii]);
         optInitY = new double[nOutputs_];
         for (ii = 0; ii < nOutputs_; ii++) optInitY[ii] = 0.0;
         optX = new double[nInputs_];
         for (ii = 0; ii < nInputs_; ii++) optX[ii] = 0.0;
         optY = new double[nOutputs_];
         for (ii = 0; ii < nOutputs_; ii++) optY[ii] = 0.0;
         optData[0] = optInitX;
         optData[1] = optInitY;
         optData[2] = optX;
         optData[3] = optY;
         funcIO = createFunctionInterfaceSimplified(ioPtr);
         status = OptimizerSearch(ioPtr, funcIO, optData, &optCnt);
         delete [] optInitX;
         delete [] optInitY;
         delete [] optX;
         delete [] optY;
         delete funcIO;
         delete ioPtr;
         pLower.clean();
         pUpper.clean();
      }

      // +++ so_ua or soua  
      else if ((!strcmp(command, "so_ua")) || (!strcmp(command, "soua")))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("soua: UQ for second order uncertainty analysis\n");
            printf("      (uncertainty in input probability distributions)\n");
            printf("syntax: soua (no argument needed).\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data to analyze (load data first).\n");
            continue;
         }
         if (nInputs_ == 1)
         {
            printf("ERROR: nInputs_ must be 2 or more.\n");
            continue;
         }
         if (faPtr != NULL) delete faPtr;
         nPtsPerDim = 64;
         sprintf(pString, "Grid resolution ? (32 - 256) ");
         nPtsPerDim = getInt(32, 256, pString);
         faFlag = 3;
         faID = -1;
         while (faID < 0)
         {
            faPtr = genFAInteractive(psuadeIO_, faFlag);
            if (faPtr == NULL) {printf("ERROR detected.\n"); continue;}
            faID = faPtr->getID();
         }
         if (faPtr == NULL) {printf("ERROR detected.\n"); continue;}
         faPtr->setNPtsPerDim(nPtsPerDim);
         faPtr->setBounds(iLowerB_, iUpperB_);
         faPtr->setOutputLevel(outputLevel_);
         psuadeIO_->getParameter("ana_outputid", pPtr);
         outputID = pPtr.intData_;
         printf("INFO: outputID = %d\n", outputID+1);
         printf("INFO: outputID obtained from the loaded sample file.\n");
         printf("      If a different output should be selected, modify\n");
         printf("      the sample file, re-load, and re-run this.\n");

         printf("To perform 2nd order uncertainty analysis, uncertainties\n");
         printf("about the input distributions are to be specified.\n");
         printf("For example, if an input has normal distribution, then\n");
         printf("uncertainties (in upper and lower bounds) of the mean\n");
         printf("and std dev are to be specified. In the following, such\n");
         printf("lower and upper bound information are to be collected.\n");
         printf("In this early implementation, no joint distributions\n");
         printf("are included in the analysis. If no 2nd order uncertainty\n");
         printf("is needed for an input, simply set the lower and upper\n");
         printf("bounds the same.\n");
         
         nParams = 2 * nInputs_;
         tempV = new double[nParams];
         tempW = new double[nParams];
         for (ii = 0; ii < nInputs_; ii++)
         {
            printf("Enter uncertainties for input %d:\n", ii+1); 
            if (inputPDFs_ == NULL || inputPDFs_[ii] == 0)
            {
               printf("Lower and upper bounds = %e %e\n",
                      iLowerB_[ii],iUpperB_[ii]);
               sprintf(pString, "Enter lower bound for input lower bound: ");
               tempV[ii*2] = getDouble(pString); 
               sprintf(pString, "Enter upper bound for input lower bound: ");
               tempW[ii*2] = getDouble(pString); 
               sprintf(pString, "Enter lower bound for input upper bound: ");
               tempV[ii*2+1] = getDouble(pString); 
               sprintf(pString, "Enter upper bound for input upper bound: ");
               tempW[ii*2+1] = getDouble(pString); 
            }
            else if (inputPDFs_[ii] == PSUADE_PDF_NORMAL)
            {
               printf("Normal distribution mean and std dev = %e %e\n",
                      inputMeans_[ii],inputStds_[ii]);
               sprintf(pString, "Enter lower bound for input mean : ");
               tempV[ii*2] = getDouble(pString); 
               sprintf(pString, "Enter upper bound for input mean : ");
               tempW[ii*2] = getDouble(pString); 
               sprintf(pString, "Enter lower bound for input std dev : ");
               tempV[ii*2+1] = getDouble(pString); 
               sprintf(pString, "Enter upper bound for input std dev : ");
               tempW[ii*2+1] = getDouble(pString); 
            }
            else if (inputPDFs_[ii] == PSUADE_PDF_LOGNORMAL)
            {
               printf("LogNormal distribution mean and std dev = %e %e\n",
                      inputMeans_[ii],inputStds_[ii]);
               sprintf(pString, "Enter lower bound for input mean : ");
               tempV[ii*2] = getDouble(pString); 
               sprintf(pString, "Enter upper bound for input mean : ");
               tempW[ii*2] = getDouble(pString); 
               sprintf(pString, "Enter lower bound for input std dev : ");
               tempV[ii*2+1] = getDouble(pString); 
               sprintf(pString, "Enter upper bound for input std dev : ");
               tempW[ii*2+1] = getDouble(pString); 
            }
            else if (inputPDFs_[ii] == PSUADE_PDF_BETA)
            {
               printf("Beta distribution alpha and beta = %e %e\n",
                      inputMeans_[ii],inputStds_[ii]);
               sprintf(pString, "Enter lower bound for input alpha : ");
               tempV[ii*2] = getDouble(pString); 
               sprintf(pString, "Enter upper bound for input alpha : ");
               tempW[ii*2] = getDouble(pString); 
               sprintf(pString, "Enter lower bound for input beta : ");
               tempV[ii*2+1] = getDouble(pString); 
               sprintf(pString, "Enter upper bound for input beta : ");
               tempW[ii*2+1] = getDouble(pString); 
            }
            else if (inputPDFs_[ii] == PSUADE_PDF_TRIANGLE)
            {
               printf("Triangle distribution mean and half width = %e %e\n",
                      inputMeans_[ii],inputStds_[ii]);
               sprintf(pString, "Enter lower bound for input mean : ");
               tempV[ii*2] = getDouble(pString); 
               sprintf(pString, "Enter upper bound for input mean : ");
               tempW[ii*2] = getDouble(pString); 
               sprintf(pString, "Enter lower bound for input half width : ");
               tempV[ii*2+1] = getDouble(pString); 
               sprintf(pString, "Enter upper bound for input half width : ");
               tempW[ii*2+1] = getDouble(pString); 
            }
            else 
            {
               printf("ERROR: Other distributions currently not supported.\n");
               continue;
               delete [] tempV;
               delete [] tempW;
            }
         }
 
         int    nValid=0, *validFlags;
         validFlags = new int[nInputs_*2];
         for (ii = 0; ii < nInputs_*2; ii++)
         {
            if (tempV[ii] != tempW[ii])
            {
               nValid++; 
               validFlags[ii] = 1;
            }
            else
            {
               validFlags[ii] = 0;
               for (jj = ii+1; jj < nInputs_*2; jj++)
               {
                  tempV[nValid+jj-ii-1] = tempV[jj];
                  tempW[nValid+jj-ii-1] = tempW[jj];
               }
            }
         }
         if (nValid == 0) 
         {
            printf("INFO: nothing to perturb.\n");
            delete [] tempV;
            delete [] tempW;
            delete [] validFlags;
            continue;
         }
         int    nSams=100, nSams2=1000;
         int    *samStates = new int[nSams];
         double *samIns  = new double[nSams*nParams];
         double *samOuts = new double[nSams];
         sampPtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_LPTAU);
         sampPtr->setPrintLevel(PL_INTERACTIVE);
         sampPtr->setInputBounds(nValid, tempV, tempW);
         sampPtr->setSamplingParams(nSams,-1,0);
         ii = 1;
         sampPtr->setOutputParams(ii);
         sampPtr->initialize(0);
         sampPtr->getSamples(nSams, nValid, ii, samIns, samOuts, samStates);
         delete sampPtr;
         sampPtr = NULL;

         int    *sPDFs  = new int[nInputs_];
         double *sMeans = new double[nInputs_];
         double *sStds  = new double[nInputs_];
         double *oneSample = new double[nInputs_];
         if (inputMeans_ == NULL)
            for (ii = 0; ii < nInputs_; ii++) sPDFs[ii] = 0;
         else
            for (ii = 0; ii < nInputs_; ii++) sPDFs[ii] = inputPDFs_[ii];
         corMat.setDim(nInputs_,nInputs_);
         ddata = 1.0;
         for (ii = 0; ii < nInputs_; ii++) corMat.setEntry(ii,ii,ddata);
         if (psPlotTool_ == 1) fp = fopen("scilabsoua.sci", "w");
         else                  fp = fopen("matlabsoua.m", "w");
         fwritePlotCLF(fp);
         for (ss = 0; ss < nSams; ss++)
         { 
            kk = 0;
            for (ii = 0; ii < nInputs_; ii++)
            {
               if (validFlags[2*ii] == 1) 
               {
                  sMeans[ii] = samIns[ss*nParams+kk];
                  kk++;
               }
               else sMeans[ii] = inputMeans_[ii];
               if (validFlags[2*ii+1] == 1) 
               {
                  sStds[ii] = samIns[ss*nParams+kk];
                  kk++;
               }
               else sMeans[ii] = inputStds_[ii];
            }
            pdfman = new PDFManager();
            pdfman->initialize(nInputs_, sPDFs, sMeans, sStds, corMat, NULL, NULL);
            vecLower.load(nInputs_, iLowerB_);
            vecUpper.load(nInputs_, iUpperB_);
            vecIn.setLength(nSams2*nInputs_);
            pdfman->genSample(nSams2, vecIn, vecLower, vecUpper);
            fprintf(fp, "Y = [\n");
            count = 0;
            for (kk = 0; kk < nSams2; kk++)
            {
               for (ii = 0; ii < nInputs_; ii++)
                  oneSample[ii] = vecIn[kk*nInputs_+ii];
               dtemp = faPtr->evaluatePoint(oneSample);
               fprintf(fp, "%e\n", dtemp);
               count++;
            }
            fprintf(fp, "];\n");
            if (psPlotTool_ == 1) fprintf(fp, "Y = gsort(Y,'g','i');\n");
            else                  fprintf(fp, "Y = sort(Y);\n");
            fprintf(fp, "X = 1:%d;\n", count);
            fprintf(fp, "X = 0.001 * X;\n");
            fprintf(fp, "plot(Y,X)\n");
            sprintf(winput, "Cumulative Distributions");
            fwritePlotTitle(fp, winput);
            fwritePlotAxes(fp);
            if (outputNames_ != NULL) sprintf(winput, "%s", outputNames_[outputID]);
            else                     sprintf(winput, "Output Values");
            fwritePlotXLabel(fp, winput);
            sprintf(winput, "Probabilities");
            fwritePlotYLabel(fp, winput);
            if (ss == 0)
            {
               if (psPlotTool_ == 1)
                    fprintf(fp, "set(gca(),\"auto_clear\",\"off\")\n");
               else fprintf(fp, "hold on\n");
            }
            delete pdfman;
         }
         fclose(fp);
         printf("Plot file for 2nd order uncertainty analysis is now");
         if (psPlotTool_ == 1) printf("in scilabsoua.sci file.\n");
         else                  printf("in matlabsoua.m file.\n");
         delete [] sPDFs;
         delete [] sMeans;
         delete [] sStds;
         delete [] oneSample;
         delete [] samIns;
         delete [] samOuts;
         delete [] samStates;
         delete [] tempV;
         delete [] tempW;
         delete [] validFlags;
         delete faPtr;
         faPtr = NULL;
      }

      // +++ ae_ua or aeua   
      else if ((!strcmp(command, "ae_ua")) || (!strcmp(command, "aeua")))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("aeua: UQ for aleatoric-epistemic uncertainty analysis\n");
            printf("syntax: aeua (no argument needed).\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
         {
            printf("ERROR: no data to analyze (load data first).\n");
            continue;
         }
         if (nInputs_ == 1)
         {
            printf("ERROR: nInputs_ must be 2 or more.\n");
            continue;
         }
         if (faPtr != NULL) delete faPtr;
         psuadeIO_->getParameter("ana_outputid", pPtr);
         iSave = pPtr.intData_;
         sprintf(pString,"Enter output number (1 - %d) : ",nOutputs_);
         outputID = getInt(1, nOutputs_, pString);
         outputID--;
         psuadeIO_->updateAnalysisSection(-1,-1,-1,-1,outputID,-1);
         psuadeIO_->getParameter("ana_rstype",pPtr);
         faType = pPtr.intData_;
         printf("The response surface model is taken from your data file:\n");
         printThisFA(faType);
         faPtr = genFA(faType, nInputs_, outputLevel_, nSamples_);
         if (faPtr == NULL) 
         {
            printf("ERROR : cannot construct response surface.\n"); 
            continue;
         }
         printf("Also, the aleatoric uncertainty information will be taken\n");
         printf("      from the INPUT section of your data file.\n");
         nPtsPerDim = 64;
         faPtr->setNPtsPerDim(nPtsPerDim);
         faPtr->setBounds(iLowerB_, iUpperB_);
         faPtr->setOutputLevel(outputLevel_);
         tempY = new double[nSamples_];
         for (ii = 0; ii < nSamples_; ii++) 
            tempY[ii] = sampleOutputs_[ii*nOutputs_+outputID];
         status = faPtr->initialize(sampleInputs_,tempY);
         delete [] tempY;
         tempY = NULL;
         printf("INFO: outputID = %d\n", outputID+1);

         printf("Step 1: select aleatoric and epistemic parameters\n");
         int nEpistemic=0;
         int *uTypes = new int[nInputs_];
         for (ii = 0; ii < nInputs_; ii++) uTypes[ii] = 0;
         kk = 1;
         while (kk > 0)
         {
            sprintf(pString,
                    "Select epistemic parameters (1 - %d, 0 if done) : ",
                    nInputs_); 
            kk = getInt(0, nInputs_, pString);
            if (kk > 0)
            {
               uTypes[kk-1] = 1;
               nEpistemic++;
            }
         } 
         if (nEpistemic == 0 || nEpistemic == nInputs_)
         {
            printf("At least 1 and at most %d epistemic parameters are\n",
                   nInputs_-1);
            printf("required for this command.\n");
            continue;
         }
         printf("You have specified %d epistemic parameters.\n",nEpistemic);

         int    nSams=20000, nSams2=2000;
         int    *samStates = new int[nSams];
         double *samIns  = new double[nSams*nEpistemic];
         double *samOuts = new double[nSams];
         double *lbs = new double[nEpistemic];
         double *ubs = new double[nEpistemic];
         kk = 0;
         for (ii = 0; ii < nInputs_; ii++)
         {
            if (uTypes[ii] != 0)
            {
               lbs[kk] = iLowerB_[ii];
               ubs[kk] = iUpperB_[ii];
               kk++;
            }
         }
         sampPtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_LHS);
         sampPtr->setPrintLevel(PL_INTERACTIVE);
         sampPtr->setInputBounds(nEpistemic, lbs, ubs);
         kk = 1;
         sampPtr->setOutputParams(kk);
         sampPtr->setSamplingParams(nSams,-1,-1);
         sampPtr->initialize(0);
         sampPtr->getSamples(nSams,nEpistemic,kk,samIns,samOuts,samStates);
         delete sampPtr;
         sampPtr = NULL;
         delete [] lbs;
         delete [] ubs;

         int    nAleatoric = nInputs_ - nEpistemic;
         int    *iPDFs2 = new int[nAleatoric];
         double *iMeans2 = new double[nAleatoric];
         double *iStds2 = new double[nAleatoric];
         double *oneSample = new double[nInputs_];
         lbs = new double[nAleatoric];
         ubs = new double[nAleatoric];
         ddata = 1.0;
         kk = 0;
         corMat.setDim(nAleatoric,nAleatoric);
         for (ii = 0; ii < nInputs_; ii++)
         {
            if (uTypes[ii] != 1)
            {
               corMat.setEntry(kk,kk,ddata);
               iPDFs2[kk] = inputPDFs_[ii]; 
               iMeans2[kk] = inputMeans_[ii]; 
               iStds2[kk] = inputStds_[ii]; 
               lbs[kk] = iLowerB_[ii];
               ubs[kk] = iUpperB_[ii];
               kk++;
            }
         }
         pdfman = new PDFManager();
         pdfman->initialize(nAleatoric,iPDFs2,iMeans2,iStds2,corMat,
                            inputPDFFiles,SPDFIndices);
         vecLower.load(nAleatoric, lbs);
         vecUpper.load(nAleatoric, ubs);
         vecIn.setLength(nSams2*nAleatoric);
         pdfman->genSample(nSams2, vecIn, vecLower, vecUpper);
         delete pdfman;
         delete [] iMeans2;
         delete [] iStds2;

         if (psPlotTool_ == 1) fp = fopen("scilabaeua.sci", "w");
         else                  fp = fopen("matlabaeua.m", "w");
         if (fp == NULL)
         {
            printf("aeua ERROR: cannot open plot file.\n");
            continue;
         }
         fwritePlotCLF(fp);
         double *Ymaxs = new double[nSams2];
         double *Ymins = new double[nSams2];
         for (ss = 0; ss < nSams2; ss++)
         {
            Ymaxs[ss] = -PSUADE_UNDEFINED;
            Ymins[ss] =  PSUADE_UNDEFINED;
         }
         ss = 0;
         int converged = 0;
         int upperCnt, lowerCnt, convergedCnt=0;
         double upperAcc, lowerAcc;
         while (ss < nSams && (converged == 0 || ss < 1000))
         {
            if (outputLevel_ > 1) printf("Epistemic sample %d\n",ss+1);
            if (ss < 50) 
            {
               fprintf(fp, "Y%d = [\n",ss+1);
               printf(".");
               fflush(stdout);
            }
            count2 = 0;
            for (ii = 0; ii < nInputs_; ii++)
            {
               if (uTypes[ii] == 1)
               {
                  oneSample[ii] = samIns[ss*nEpistemic+count2];
                  count2++;
               }
            }
            count = lowerCnt = upperCnt = 0;
            upperAcc = lowerAcc = 0;
            for (kk = 0; kk < nSams2; kk++)
            {
               flag = 0;
               for (ii = 0; ii < nAleatoric; ii++)
               {
                  if (vecIn[kk*nAleatoric+ii] < lbs[ii] ||
                      vecIn[kk*nAleatoric+ii] > ubs[ii]) flag++;
               }
               if (flag == 0)
               {
                  count2 = 0;
                  for (ii = 0; ii < nInputs_; ii++)
                  {
                     if (uTypes[ii] == 0)
                     {
                        oneSample[ii] = vecIn[kk*nAleatoric+count2];
                        count2++;
                     }
                  }
                  dtemp = faPtr->evaluatePoint(oneSample);
                  if (ss < 50) fprintf(fp, "%e\n", dtemp);
                  if (dtemp > Ymaxs[kk]) 
                  {
                     if (Ymaxs[kk] != 0) upperAcc += PABS((Ymaxs[kk]-dtemp)/Ymaxs[kk]);
                     Ymaxs[kk] = dtemp;
                     upperCnt++;
                  }
                  if (dtemp < Ymins[kk])
                  {
                     if (Ymins[kk] != 0) lowerAcc += PABS((Ymins[kk]-dtemp)/Ymins[kk]);
                     Ymins[kk] = dtemp;
                     lowerCnt++;
                  }
                  count++;
               }
            }
            ddata = 100.0 * upperAcc / nSams;
            if (outputLevel_ > 2 && ddata > 0.0 && ss > 0) 
               printf("  Upper lifted  %7.3f percent on average from last time\n",ddata);
            ddata = 100.0 * lowerAcc / nSams;
            if (outputLevel_ > 2 && ddata > 0.0 && ss > 0) 
               printf("  Lower dropped %7.3f percent on average from last time\n",ddata);
            if (upperCnt+lowerCnt == 0) convergedCnt++;
            else                        convergedCnt = converged = 0;
            if (outputLevel_ > 1) 
               printf("  Convergence indicator = %5d (0 20 times => converged)\n",
                      upperCnt+lowerCnt);
            if (convergedCnt >= 20) converged = 1;
            if (count < 0.5 * nSams2)
            {
               printf("WARNING: < half of the points are within bounds.\n");
               printf("         Input ranges may need to be widened.\n");
            }
            if (count == 0)
            {
               printf("ERROR: none of the sample points are within bounds.\n");
               continue;
            }
            if (ss < 50) fprintf(fp, "];\n");
            if (ss < 50) 
            {
               if (psPlotTool_ == 1) 
                    fprintf(fp,"Y%d = gsort(Y%d,'g','i');\n",ss+1,ss+1);
               else fprintf(fp,"Y%d = sort(Y%d);\n",ss+1,ss+1);
               fprintf(fp, "X = 1:%d;\n",count);
               fprintf(fp, "X = X / %d;\n", count);
               fprintf(fp, "plot(Y%d,X)\n",ss+1);
               fprintf(fp, "drawnow\n");
               if (ss == 0)
               {
                  sprintf(winput, "Cumulative Distributions");
                  fwritePlotTitle(fp, winput);
                  fwritePlotAxes(fp);
                  if (outputNames_ != NULL) 
                       sprintf(winput,"%s",outputNames_[outputID]);
                  else sprintf(winput,"Output Values");
                  fwritePlotXLabel(fp, winput);
                  sprintf(winput, "Probabilities");
                  fwritePlotYLabel(fp, winput);
                  if (psPlotTool_ == 1)
                       fprintf(fp, "set(gca(),\"auto_clear\",\"off\")\n");
                  else fprintf(fp, "hold on\n");
               }
            }
            ss++;
         }
         printf("\n");
         count = 0;
         fprintf(fp, "YU = [\n");
         for (ss = 0; ss < nSams2; ss++)
         {
            if (Ymaxs[ss] > -PSUADE_UNDEFINED) 
            {
               fprintf(fp, "%e\n", Ymaxs[ss]);
               count++;
            }
         }
         fprintf(fp, "];\n");
         if (psPlotTool_ == 1) 
              fprintf(fp,"YU = gsort(YU,'g','i');\n");
         else fprintf(fp,"YU = sort(YU);\n");
         fprintf(fp, "X = 1:%d;\n",count);
         fprintf(fp, "X = X / %d;\n", count);
         fprintf(fp, "plot(YU,X,'r-','lineWidth',3)\n");

         count = 0;
         fprintf(fp, "YL = [\n");
         for (ss = 0; ss < nSams2; ss++)
         {
            if (Ymins[ss] < PSUADE_UNDEFINED) 
            {
               fprintf(fp, "%e\n", Ymins[ss]);
               count++;
            }
         }
         fprintf(fp, "];\n");
         if (psPlotTool_ == 1) 
              fprintf(fp,"YL = gsort(YL,'g','i');\n");
         else fprintf(fp,"YL = sort(YL);\n");
         fprintf(fp, "X = 1:%d;\n",count);
         fprintf(fp, "X = X / %d;\n", count);
         fprintf(fp, "plot(YL,X,'r-','lineWidth',3)\n");
         fclose(fp);
         printf("Plot file for aleatoric-epistemic analysis is now ");
         if (psPlotTool_ == 1) printf("in scilabaeua.sci.\n");
         else                  printf("in matlabaeua.m.\n");
         psuadeIO_->updateAnalysisSection(-1,-1,-1,-1,iSave,-1);
         delete [] oneSample;
         delete [] samIns;
         delete [] samOuts;
         delete [] samStates;
         delete [] uTypes;
         delete [] lbs;
         delete [] ubs;
         delete [] Ymaxs;
         delete [] Ymins;
         delete faPtr;
         faPtr = NULL;
         oneSample = NULL;
         samIns = NULL;
         samOuts = NULL;
         samStates = NULL;
      }

      // +++ sys
      else if (!strcmp(command, "sys"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("sys: execute a system command \n");
            printf("syntax: sys <command>\n");
            continue;
         }
         count = 0;
         for (ii = 0; ii < 498; ii++)
            if (lineIn[ii] == 's' && lineIn[ii+1] == 'y' &&
                lineIn[ii+2] == 's') break;
         if (ii != 498) system(&lineIn[ii+4]);
      }

      // +++ ivec_create
      else if (!strcmp(command,"ivec_create") || !strcmp(command,"vcreate"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("ivec_create: create an internal vector (for rseval)\n");
            printf("syntax: ivec_create or vcreate (no argument needed)\n");
            printf("This command is used together with rscreate and rseval\n");
            printf("to create a response surface on the fly and evaluate\n");
            printf("a new sample point placed into the local register.\n");
            continue;
         }
         if (nInputs_ <= 0 || psuadeIO_ == NULL)
         {
            printf("ERROR: no data to analyze (load data first).\n");
            continue;
         }
         if (dataReg_ != NULL) delete [] dataReg_;
         dataReg_ = new double[nInputs_];
         for (ii = 0; ii < nInputs_; ii++) 
            dataReg_[ii] = iLowerB_[ii] + PSUADE_drand() *
                           (iUpperB_[ii] - iLowerB_[ii]);; 
         printf("Internal sample vector created and the values have been set");
         printf(" to be the mid points.\n");
      }

      // +++ ivec_modify
      else if (!strcmp(command,"ivec_modify") || !strcmp(command,"vmodify"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("ivec_modify: modify an entry in the internal vector\n");
            printf("syntax: ivec_modify or vmodify (no argument needed)\n");
            printf("This command is used after ivec_create to modify\n");
            printf("individual entries in the internal register.\n");
            continue;
         }
         if (dataReg_ == NULL)
         {
           printf("ERROR: use ivec_create before this command.\n");
           continue;
         }
         sprintf(pString, "Which entry to modify ? (1 - %d) ", nInputs_);
         ind = getInt(1, nInputs_, pString);
         ind--;
         sprintf(pString, "What is the new value ? (%e - %e) ",
                 iLowerB_[ind], iUpperB_[ind]);
         ddata = getDouble(pString);
         if (ddata < iLowerB_[ind] || ddata > iUpperB_[ind])
            printf("WARNING: data out of range (extrapolation).\n");
         dataReg_[ind] = ddata;
      }

      // +++ ivec_show
      else if (!strcmp(command,"ivec_show") || !strcmp(command, "vshow"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("ivec_show: display the internal vector\n");
            printf("syntax: ivec_show or vshow (no argument needed)\n");
            printf("This command is used together with ivec_create and\n");
            printf("ivec_modify.\n");
            continue;
         }
         if (dataReg_ == NULL)
         {
           printf("ERROR: use ivec_create before this command.\n");
           continue;
         }
         for (ii = 0; ii < nInputs_; ii++)
            printf("Input %3d = %e\n", ii+1, dataReg_[ii]);
      }

      // +++ showformat
      else if (!strcmp(command, "showformat"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("showformat: display different PSUADE file formats \n");
            printf("syntax: showformat (no argument needed)\n");
            printf("User-provided information are often needed in design\n");
            printf("and analysis. These information are provided by users\n");
            printf("to PSUADE at various stages. This command lists many\n");
            printf("such file formats.\n");
            continue;
         }
         printDashes(PL_INFO, 0);
         printf("PSUADE data file format (can be read by load): \n");
         printf("  run 'geninputfile' command to create an example.\n");
         printDashes(PL_INFO, 0);
         printf("Standard file format (can be read by read_std): \n");
         printf("  line 1: nSamples nInputs nOutputs\n");
         printf("  line 2: <sample 1 inputs> < sample 1 outputs>\n");
         printf("  line 3: <sample 2 inputs> < sample 2 outputs>\n");
         printf("  .......\n");
         printDashes(PL_INFO, 0);
         printf("Xls file format (can be read by read_xls): \n");
         printf("  line 1: nSamples nInputs nOutputs\n");
         printf("  line 2: (optional) # sample input and output names\n");
         printf("  line 3: 1 <sample 1 inputs> < sample 1 outputs>\n");
         printf("  line 4: 2 <sample 2 inputs> < sample 2 outputs>\n");
         printf("  .......\n");
         printDashes(PL_INFO, 0);
         printf("MCMC experimental data file format (O1 = output 1): \n");
         printf("  line 1: PSUADE_BEGIN\n");
         printf("  line 2: nExperiments(p) nOutputs(n) nDesignInputs designInputList\n");
         printf("  line 3: 1 <designInput values...> <O1 mean> <O1 std dev> ... <On std dev> \n");
         printf("  line 4: 2 <designInput values...> <O1 mean> <O1 std dev> ... <On std dev> \n");
         printf("  ...\n");
         printf("  line  : p <designInput values...> <O1 mean> <O1 std dev> ... <On std dev> \n");
         printf("  line  : PSUADE_END\n");
         printDashes(PL_INFO, 0);
         printf("Sample Input Only file format (used in PDF S type, rs_uab, iread): \n");
         printf("  line 1: PSUADE_BEGIN\n");
         printf("  line 2: <number of sample points> <number of inputs>\n");
         printf("  line 3: (optional) : '#' followed by input names\n");
         printf("  line 4: 1 sample point 1 inputs \n");
         printf("  line 5: 2 sample point 2 inputs \n");
         printf("  line 6: 3 sample point 3 inputs \n");
         printf("  ...\n");
         printf("  line n: PSUADE_END\n");
         printDashes(PL_INFO, 0);
         printf("RSConstraints file format (used in Analysis: rs_constraint): \n");
         printf("  line 1: nInputs\n ");
         printf("  line 2: <input (or 0)> <value (nominal val if 0)> \n ");
         printf("  line 3: <input (or 0)> <value (nominal val if 0)> \n ");
         printf("  ... \n");
         printDashes(PL_INFO, 0);
         printf("RS index file format (used in Analysis: rs_index_file): \n");
         printf("  line 1: nInputs in rs data (driver) file\n");
         printf("  line 2: 1 <num> <default if num == 0>\n");
         printf("  line 3: 2 <num> <0 if num != 0>\n");
         printf("  line 4: 3 <num> <default if num == 0>\n");
         printf("  line 5: 4 <num> <0 if num != 0>\n");
         printf("  ...\n");
         printDashes(PL_INFO, 0);
         printf("MOATConstraints file format (used in Analysis: moat_constraint): \n");
         printf("  line 1: nInputs \n");
         printf("  line 2: <input (or 0)> <value (nominal val if 0)> \n");
         printf("  line 3: <input (or 0)> <value (nominal val if 0)> \n");
         printf("  ... \n");
         printDashes(PL_INFO, 0);
      }

      // +++ setdriver 
      else if (!strcmp(command, "setdriver"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("setdriver: set the driver field in the PSUADE file\n");
            printf("syntax: setdriver\n");
            continue;
         }
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: data not loaded\n");
            continue;
         }
         printf("Which driver do you want to set: \n");
         printf("1. simulation driver (driver) \n");
         printf("2. optimization driver (opt_driver) \n");
         printf("3. auxiliary optimization driver (aux_opt_driver) \n");
         printf("4. ensemble simulation driver (ensemble_driver) \n");
         printf("5. ensemble optimization driver (ensemble_opt_driver) \n");
         sprintf(pString, "Enter (1 - 5): ");
         kk = getInt(1,5,pString);
         printf("Enter name of the driver : ");
         fgets(dataFile, 5000, stdin);
         dataFile[strlen(dataFile)-1] = '\0';
         for (ii = 0; ii < strlen(dataFile)-1; ii++)
            if (dataFile[ii] == ' ') break;
         if (ii == strlen(dataFile)-1 && (fp=fopen(dataFile,"r")) == NULL)
         {
            printf("WARNING: file %s not found.\n", dataFile);
         }
         if (kk == 1)
            psuadeIO_->updateApplicationSection(dataFile,NULL,NULL,NULL,NULL,-1);
         else if (kk == 2)
            psuadeIO_->updateApplicationSection(NULL,dataFile,NULL,NULL,NULL,-1);
         else if (kk == 3)
            psuadeIO_->updateApplicationSection(NULL,NULL,dataFile,NULL,NULL,-1);
         else if (kk == 4)
            psuadeIO_->updateApplicationSection(NULL,NULL,NULL,dataFile,NULL,-1);
         else if (kk == 5)
            psuadeIO_->updateApplicationSection(NULL,NULL,NULL,NULL,dataFile,-1);
         printf("Use 'write' to update your PSUADE input file.\n");
      }

      // +++ start_matlab 
      else if ((!strcmp(command, "start_matlab")))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("start_matlab: start matlab from within PSUADE\n");
            printf("              interactive mode.\n");
            printf("syntax: start_matlab (no argument needed)\n");
            continue;
         }
         strcpy(command, "which matlab");
         status = system(command);
         if (status != 0)
            printf("Matlab not found (have you set the path?)\n");
         else
         {   
            printf("Two modes to start Matlab: \n");
            printf("1. desktop mode (create a matlab window)\n");
            printf("2. nodesktop mode (use current window for matlab)\n");
            sprintf(pString,"Start Matlab in desktop mode (y or n)? ");
            getString(pString, winput);
            if (winput[0] == 'y') strcpy(command, "matlab");
            else                  strcpy(command, "matlab -nodesktop");
            status = system(command);
         }
      }

      // +++ checkformat 
      else if ((!strcmp(command, "checkformat")))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("checkformat: check format of files used by PSUADE.\n");
            printf("syntax: checkformat (no argument needed)\n");
            continue;
         }
         printf("Please select from the following file formats to check: \n");
         printf("1. OUU Z3 sample file.\n");
         printf("2. OUU Z4 sample file.\n");
         printf("3. psuade input file\n");
         printf("4. MCMC specification file\n");
         printf("5. MCMC posterior file\n");
         printf("6. S-type PDF sample file\n");
         sprintf(pString, "Enter your selection : ");
         int option = getInt(1,6,pString);
         printf("Please enter your sample file: ");
         scanf("%s", cString);
         fgets(winput, 1000, stdin);
         if (option == 1 || option == 2)
         {
            sprintf(pString, "How many inputs are in this sample file? ");
            kk = getInt(1,10000000,pString);
            status = checkOUUFileFormat(cString, option, kk, outputLevel_);
         }
         else if (option == 3)
         {
            ioPtr = new PsuadeData();
            status = ioPtr->readPsuadeFile(cString);
            if (status != 0) 
            {
               printf("PSUADE data file format is not valid.\n");
               status = -1;
            }
            ioPtr->getParameter("input_ninputs", pPtr);
            if (pPtr.intData_ <= 0) 
            {
               printf("PSUADE data file has no inputs.\n");
               status = -1;
            }
            ioPtr->getParameter("input_sample", pPtr);
            if (pPtr.dbleArray_ == NULL) 
                 printf("PSUADE data file has no sample.\n");
            else delete [] pPtr.dbleArray_;
            pPtr.dbleArray_ = NULL;

         }
         else if (option == 4)
         {
            status = checkMCMCFileFormat(cString, 0, outputLevel_);
         }
         else if (option == 5)
         {
            status = checkMCMCFileFormat(cString, 1, outputLevel_);
         }
         else if (option == 6)
         {
            status = checkSPDFFileFormat(cString, outputLevel_);
         }
         if (status == 0) printf("PASSED: file format validated.\n");
         else             printf("FAILED: invalid file format.\n");
      }

      // +++ sample_info or sinfo
      else if (!strcmp(command, "sample_info") || !strcmp(command, "sinfo"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("sample_info: (or sinfo) display loaded sample information.\n");
            printf("syntax: sample_info or sinfo (no argument needed)\n");
            continue;
         }
         printf("Sample in memory: nSamples = %d\n", nSamples_);
         printf("                  nInputs  = %d\n", nInputs_);
         printf("                  nOutputs = %d\n", nOutputs_);
         if (inputNames_ != NULL)
         {
            printf("Input names: \n");
            for (ii = 0; ii < nInputs_; ii++)
               if (inputNames_[ii] != NULL)
                  printf("  Input %4d: %s\n",ii+1,inputNames_[ii]);
         }
         if (outputNames_ != NULL)
         {
            printf("Output names: \n");
            for (ii = 0; ii < nOutputs_; ii++)
               if (outputNames_[ii] != NULL)
                  printf("  Output %4d: %s\n",ii+1,outputNames_[ii]);
         }
         count = 0;
         for (ss = 0; ss < nSamples_; ss++) if (sampleStates_[ss] == 1) count++;
         printf("Number of valid sample points (state = 1)  = %d\n",count);
         count = 0;
         for (ss = 0; ss < nSamples_; ss++)
         {
            flag = 1;
            for (ii = 0; ii < nOutputs_; ii++)
               if (sampleOutputs_[ss*nOutputs_+ii] == PSUADE_UNDEFINED) flag = 0;
            if (flag == 1) count++;
         }
         printf("Number of sample points with valid outputs = %d\n",count);
      }

      // +++ quit 
      else if ((!strcmp(command, "quit")) || (!strcmp(command, "q")))
      {
         break;
      }
      else if ((!strcmp(command, "exit")))
      {
         break;
      }
      else if (command[0] == '#')
      {
         printf("#\n");
      }
      else
      {
         sprintf(scriptName, "%s.psu", command);
         fp = fopen(scriptName, "r");
         if (fp != NULL && scriptMode == 0)
         {
            scriptFp = fp;
            printf("Script file %s found.\n", scriptName);
            scriptMode = 1;
         }
         else if (fp != NULL && scriptMode == 1)
         {
            printf("ERROR: only one level of script interpretation allowed.\n");
            fclose(fp);
         }
         else
         {
            printf("command %s not recognized\n", command);
            fflush(stdout);
         }
      }
   }

   // quit psuade interactive session
   if (faPtr != NULL) delete faPtr;
   if (faPtrsRsEval != NULL)
   {
      for (ii = 0; ii < nOutputs_; ii++)
         if (faPtrsRsEval[ii] != NULL) delete faPtrsRsEval[ii];
      delete faPtrsRsEval;
   }
   delete currSession;
   faPtr = NULL;
   faPtrsRsEval = NULL;
   currSession = NULL;
   return 0;
}

