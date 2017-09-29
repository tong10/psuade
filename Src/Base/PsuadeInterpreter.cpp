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
#include "PDFBase.h"
#include "PDFManager.h"
#include "PDFNormal.h"
#include "Matrix.h"
#include "Vector.h"
#include "Sampling.h"
#include "FunctionInterface.h"
#include "PsuadeData.h"
#include "Optimizer.h"

// ------------------------------------------------------------------------
// local defines 
// ------------------------------------------------------------------------
#define PABS(x)  ((x) > 0 ? x : -(x))

// ************************************************************************
// interpret command from interactive session
// ------------------------------------------------------------------------
int PsuadeBase::interpretInteractive()
{
   int    nInputs=0, nSamples, nOutputs, status, *sampleStates=NULL, nReps;
   int    iInd, outputID, ind, ind2, iInd1, iInd2, count;
   int    it, nPtsPerDim=64, iplot1, iplot2, jplot=-1, faLeng, sInd, flag;
   int    ss, ii, jj, kk, ll, iplot3, *states, rsiNOutputs, *rsiSet, count2;
   int    **rsiMatrix, analysisMethod, samplingMethod, faType, setCompare;
   int    otrans, *xsforms, nOut2, oInd, *tags, faFlag, iOne=1;
   int    nPaths, nTrials=100, trial, *indSet, nSmooth, currP, faLimit;
   int    nParts, nPlots, *plotIndices=NULL, faID, nFiles, iplot4;
   int    inputID, oplot1, oplot2, scriptMode=0, *tempI, *tagArray=NULL;
   int    method, saveDiag, saveMode1, saveMode2, saveMethod, usePDFs=0;
   int    nParams, *iPDFs, nInputs2, nSamples2, commandCnt, iSave;
   long   nSamplesLong;
   double *sampleInputs=NULL, *sampleOutputs=NULL, *tempX, *tempY, *tempW;
   double *iLowerB=NULL, *iUpperB=NULL, Ymax, Ymin, Xmin, Xmax, width;
   double *tempInds, *tempT, *inputSettings, *faYIn, *faXOut, *faYOut;
   double thresh, threshL, threshU, GYmax, GYmin, dtemp;
   double *tempV, **moatSample, currX1, currX2, currY1, gamma;
   double currY2, ddata, filterRange, sLo, sHi, *dataReg=NULL;
#ifdef HAVE_SVM
   double sumErr, maxErr, tolerance;
#endif
   double *threshLs, *threshUs, diagMax, diagMin, currX1F, currX2F;
   double ***rsi3Matrix, aVal, bVal, minDist, *tempM, *iMeans, *iStds;
   string sparam;
   char   command[100], dataFile[100], **inputNames=NULL, **names, *names2;
   char   **outputNames=NULL, winput[501], lineIn[500], pString[500];
   char   dirName[500],errFile[100],*targv[8],scriptName[500], cString[500];
   char   subdirName[500];
   FILE   *fp, *fpOut, *fErr, *scriptFp=NULL;
   FuncApprox *faPtr=NULL, **faPtrs=NULL, **faPtrsRsEval=NULL;
   pData  pPtr, pINames, pLower, pUpper, pONames, pInputs, pOutputs, pStates;
   pData  *pdata, pPDFs, pMeans, pStds;
   aData  aPtr;
   PsuadeData *ioPtr=NULL;
   Sampling   *sampPtr=NULL, *sampAux;;
   FunctionInterface *funcIO=NULL;
   AnalysisManager   *anaManager;
   PDFManager *pdfman;
   Matrix     corMat;
   Vector     vecIn, vecOut, vecUpper, vecLower;

   printf("PSUADE - A Problem Solving environment for \n");
   printf("         Uncertainty Analysis and Design Exploration (%d.%d.%d)\n",
          psuade_VERSION_MAJOR, psuade_VERSION_MINOR, psuade_VERSION_PATCH);
   printf("(for help, enter <help>)\n");
   printEquals(0);
   commandCnt = 0;
   while (1)
   {
      if (scriptMode == 1 && scriptFp != NULL)
      {
         fgets(lineIn, 500, scriptFp);
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
         fgets(lineIn,500,stdin); 
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
            printf("\t\t    Analyzer: use rs_ua (with system-generated sample from rs)\n");
            printf("\t\t 3. Sampling: first construct response surface (rs)\n");
            printf("\t\t    Analyzer: use rs_uap (user provides a 'posterior' sample)\n");
            printf("\tII.   Parameter Screening: \n");
            printf("\t\t 1. Sampling: MOAT, GMOAT\n");
            printf("\t\t    Analyzer: moat\n");
            printf("\t\t 2. Sampling: LH, LPTAU\n");
            printf("\t\t    Analyzer: mars_sa (MARS importance screening)\n");
            printf("\t\t 3. Sampling: LH, LPTAU\n");
            printf("\t\t    Analyzer: delta_test\n");
            printf("\t\t 4. Sampling: LH, MC\n");
            printf("\t\t    Analyzer: sot_sa (sum-of-trees screening)\n");
            printf("\t\t 5. Sampling: FF4, FF5\n");
            printf("\t\t    Analyzer: ff (Fractional Factorial screening)\n");
            printf("\t\t 6. Sampling: LSA\n");
            printf("\t\t    Analyzer: lsa (local sensitivity analysis)\n");
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
            printf("\t\t    (e) perform generalization test (rstest_gt, most rigorous)\n");
            printf("\tV.   Global sensitivity analysis (first order): \n");
            printf("\t\t 1. Sampling: replicated Latin hypercube\n");
            printf("\t\t    Analyzer: use 'me'\n");
            printf("\t\t 2. Sampling: any space-filling (large enough) sample\n");
            printf("\t\t    Analyzer: use 'me'\n");
            printf("\t\t 3. Sampling: any space-filling sample\n");
            printf("\t\t    Analyzer: use rscheck with Legendre and select scaling\n");
            printf("\t\t 4. Sampling: first construct response surface\n");
            printf("\t\t    Analyzer: use rssobol1b in command line mode\n");
            printf("\tVI.  Global sensitivity analysis (second order): \n");
            printf("\t\t 1. Sampling: replicated OA\n");
            printf("\t\t    Analyzer: use 'ie' in command line mode\n");
            printf("\t\t 2. Sampling: any space-filling (large enough) sample\n");
            printf("\t\t    Analyzer: use 'ie'\n");
            printf("\t\t 3. Sampling: first construct response surface\n");
            printf("\t\t    Analyzer: use rssobol2 in command line mode\n");
            printf("\tVII.  Global sensitivity analysis (total order): \n");
            printf("\t\t 1. Sampling: use any space-filling (very large) sample\n");
            printf("\t\t    Analyzer: use tsi command (approximate analysis)\n");
            printf("\t\t 2. Sampling: first construct response surface\n");
            printf("\t\t    Analyzer: use rssoboltsib in command line mode, or\n");
            printf("\t\t 3. Sampling: use FAST sampling\n");
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
            printf("\t\t 1. Sampling: LPTAU, LH, OA, METIS to create response surface\n");
            printf("\t\t    Analyzer: rsmcmc (in command line mode), or\n");
            printf("\t\t 2. Sampling: LPTAU, LH, OA, METIS to create response surface\n");
            printf("\t\t    Analyzer: mcmc (need special re-compilation for efficiency)\n");
            printf("\tXII.  Advanced features: \n");
            printf("\t\t 1. Impose constraints in sampling and analysis (e.g. moatgen)\n");
            printf("\t\t 2. Plot PDF of std dev. (error) of response surfaces (rssd_ua)\n");
            printf("\t\t 3. Intersection or Bayes-like rules (e.g. rsi2, rsvol)\n");
            printf("\t\t 4. Multi-objective optimization (mo_opt)\n");
            printf("\t\t 5. Mixed aleatoric-epistemic uncertainty analysis (ae_ua)\n");
            printf("\t\t 6. Second order analysis (so_ua) - uncertainty in distribution\n");
            printf("\t\t 7. Tools for setting up user application with PSUADE\n");
            printf("\t\t 8. Sampling refinements (uniform/adaptive: refine, a_refine)\n");
         }
         else if (!strcmp(winput, "io"))
         {
            printf("Commands for reading/write/updating data to/from files:\n");
            printf("(to see more io commands, use 'help io long')\n");
            printf("\tload         <file> (load a data file) \n");
            printf("\tloadmore     <file> (add data to current data set)\n");
            printf("\twrite        <file> (write to file in PSUADE format)\n");
            printf("\tnwrite       <file> (write unevaluated points only)\n");
            printf("\tiwrite       <file> (write to file with only inputs)\n");
            printf("\towrite       <file> (write to file with only outputs)\n");
            printf("\tread_std     <file> (read data in standard format)\n");
            printf("\twrite_std    <file> (write to file in standard format)\n");
            printf("\twrite_matlab <file> (write to file in matlab format)\n");
            printf("\tread_xls     <file> (read data in Excel format)\n");
            printf("\twrite_xls    <file> (write to file in Excel format)\n");
            printf("\tupdate       <file> (update sample OUTPUTS from a PSUADE file)\n");
            printf("\tiadd         <file> (add more inputs from a PSUADE file)\n");
            printf("\toadd         <file> (add more outputs from a PSUADE file)\n");
            printf("\tireplace     <file> (replace one input from <file>)\n");
            printf("\toreplace     <file> (replace all outputs from <file>)\n");
            printf("\tsplitdata           (split into 2 data files (in PSUADE format)\n");
         }
         else if (!strcmp(winput, "stat"))
         {
            printf("Commands for basic statistic on raw data:\n");
            printf("\tua         (uncertainty analysis + matlab/scilab plot)\n");
            printf("\trs_ua      (uncertainty analysis on fuzzy response surface)\n");
            printf("\trs_uab     (RS-based UA with bagging and optional discrepancy)\n");
            printf("\trs_uap     (RS-based UA with posterior sample from mcmc)\n");
            printf("\tca         (correlation analysis)\n");
            printf("\tanova      (analyis of variation)\n");
            printf("\t1stest     (1-sample test (Chi-squared, dist fit))\n");
            printf("\t2stest     (2-sample test (T-test,K-S,Mann-Whitney))\n");
            printf("\tgendist    (create a sample given an PDF (for 1 input))\n");
            printf("\tpdfconvert (convert a sample using selected distribution)\n");
            printf("\trand_draw  (draw a bootstrapped sample from the loaded sample)\n");
            printf("\trand_draw2 (draw a sample from 2 samples with 2 set of inputs)\n");
         }
         else if (!strcmp(winput, "screen"))
         {
            printf("Commands for parameter screening:\n");
            printf("\tmoat       (parameter screening with the Morris method)\n");
            printf("\tdelta_test (parameter screening with Delta test)\n");
            printf("\tsot_sa     (parameter screening with sum-of-trees method)\n");
            printf("\tff         (parameter screening with Fractional factorial)\n");
            printf("\tlsa        (parameter screening with local SA method)\n");
            printf("\tmars_sa    (parameter screening with MARS)\n");
            printf("\tgp_sa      (parameter screening with Gaussian process)\n");
            printf("\tpca        (principal component analysis: output screening)\n");
         }
         else if (!strcmp(winput, "rs"))
         { 
            printf("Commands for response surface analysis:\n");
            printf("\tsvmfind     (search for good parameters for SVM)\n");
            printf("\trscheck     (check quality of RS with CV and rich graphics)\n");
            printf("\trstest_hs   (check quality of RS with another test set)\n");
            printf("\trstest_ts   (check quality of RS with the training set)\n");
            printf("\trstest_cv   (check quality of RS with cross validation)\n");
            printf("\trstest_gt   (check quality of RS with generalization test)\n");
            printf("\trscreate    (create response surface to be used by rseval)\n");
            printf("\tivec_create (create random input vector in local memory)\n");
            printf("\tivec_modify (modify 1 entry of input vector in local memory)\n");
            printf("\tivec_show   (display the input vector in local memory)\n");
            printf("\trseval      (evaluate RS at given points)\n");
            printf("\trseval_m    (use PSUADE as response surface server)\n");
            printf("\trs_splot    (create scatter plots on response surface)\n");
            printf("\trstgen      (generate a sample (FF/FACT) for rstest_hs)\n");
            printf("\trsvol       (compute percentage volume of constrained region)\n");
            printf("\tint         (integration: find volume under response surface)\n");
         }
         else if (!strcmp(winput, "qsa") && !strcmp(pString, "long"))
         {
            printf("Commands for quantitative uncertainty and sensitivity analysis:\n");
            printf("\tme         (main effect analysis + matlab plot on any sample)\n");
            printf("\tie         (2-way interaction effect analysis on any sample)\n");
            printf("\ttsi        (total sensitivity analysis on any sample)\n");
            printf("\trssobol1   (RS-based Sobol' main effect)\n");
            printf("\trssobol2   (RS-based Sobol' interaction effect)\n");
            printf("\trssoboltsi (RS-based Sobol' total effect)\n");
            printf("\trssobolg   (RS-based Sobol' group main effect)\n");
            printf("\trssobol1b  (RS-based Sobol' main effect with bootstrapping)\n");
            printf("\trssoboltsib(RS-based Sobol' total effect with bootstrapping)\n");
            printf("\tso_ua      (RS-based 2nd order analysis - uncertainty in PDFs)\n");
            printf("\tae_ua      (RS-based mixed aleatoric-epistemic analysis)\n");
         }
         else if (!strcmp(winput, "qsa"))
         {
            printf("Commands for quantitative uncertainty and sensitivity analysis:\n");
            printf("(to see more qsa commands, use 'help qsa long')\n");
            printf("\tme         (main effect analysis + matlab plot on any sample)\n");
            printf("\tie         (2-way interaction effect analysis on any sample)\n");
            printf("\ttsi        (total sensitivity analysis on any sample)\n");
            printf("\trssobol1b  (RS-based Sobol' main effect with bootstrapping)\n");
            printf("\trssobol2   (RS-based Sobol' interaction effect)\n");
            printf("\trssobolg   (RS-based Sobol' group main effect)\n");
            printf("\trssoboltsib(RS-based Sobol' total effect with bootstrapping)\n");
            printf("\tso_ua      (RS-based 2nd order analysis - uncertainty in PDFs)\n");
            printf("\tae_ua      (RS-based mixed aleatoric-epistemic analysis)\n");
         }
         else if (!strcmp(winput, "calibration"))
         {
            printf("Commands for optimization/calibration:\n");
            printf("\trsmcmc     (RS-based Bayesian inversion using MCMC)\n");
            printf("\tmo_opt     (RS-based multi-objective optimization)\n");
            printf("\tmcmc       (Simulation-based MCMC: advanced feature)\n");
         }
         else if (!strcmp(winput, "plot") && !strcmp(pString, "long"))
         { 
            printf("Commands for plotting sample data:\n");
            printf("\tsplot      (scatter plot in matlab/scilab)\n");
            printf("\tsplot2     (2-input/1 output scatter plot in matlab/scilab)\n");
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
            printf("\tiplot_pdf  (plot the sample PDF of selected inputs)\n");
            printf("\toplot2     (2-output scatter plot)\n");
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
            printf("\tishuffle      (re-arrange the orders of the input parameters)\n");
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
            printf("\titran         (transform certain input (log or power)\n");
            printf("\totran         (transform certain output (log or power)\n");
            printf("\titag          (tag sample points based on input values)\n");
            printf("\totag          (tag sample points based on output values)\n");
            printf("\tpurge         (take out invalid sample points)\n");
            printf("\trm_duplicates (take out duplicate sample points)\n");
            printf("\toop           (replace output1 = a * out2 + b * out3)\n");
            printf("\tsetdriver <s> (set the driver field to <s>)\n");
         }
         else if (!strcmp(winput, "misc"))
         {
            printf("Miscellaneous commands:\n");
            printf("\trun <file>      (run a psuade input script) \n");
            printf("\tquit (exit)     (terminate command line session)\n");
            printf("\trsmax <d>       (set maximum number of data points for RS)\n");
            printf("\tsys <command>   (execute a system command)\n");
            printf("\tprintlevel <d>  (set print level)\n");
            printf("\tinteractive     (turn on/off advanced interactive mode: obsolete)\n");
            printf("\toutput_file     (set the default output file name\n");
            printf("\tsqc             (sample quality check: distance metric)\n");
            printf("\tnna             (nearest neighbor analysis: for outliers\n");
            printf("\tshow_format     (display file formats for RS/MOAT constraints)\n");
            printf("\tscilab          (turn on/off scilab (limited support))\n");
            printf("\tstart_matlab    (start matlab in PSUADE interactive mode)\n");
         }
         else if (!strcmp(winput, "advanced"))
         {
            printf("Advanced analysis and control commands:\n");
            printf("\tio_expert       (turn on/off IO expert mode)\n");
            printf("\trs_expert       (turn on/off response surface expert mode)\n");
            printf("\tana_expert      (turn on/off analysis expert mode)\n");
            printf("\tsam_expert      (turn on/off sampling expert mode)\n");
            printf("\topt_expert      (turn on/off optimization expert mode)\n");
            printf("\tgenconfigfile   (create a template PSUADE config file)\n");
            printf("\tuse_configfile <file> (use config file to select parameters)\n");
            printf("\trefine          (refine a sample: higher mesh resolution)\n");
            printf("\ta_refine        (adaptively refine a sample based on MarsBag)\n");
            printf("\tinterface_track (track the interface separating Y=0 vs Y=1.)\n");
            printf("\tmoatrepair   <file>  (repair MOAT sample from <file>)\n");
            printf("\tgmoatrepair  <file>  (repair GMOAT sample from <file>)\n");
            printf("\tmoatgen              (create MOAT repair file: 1 constr)\n");
            printf("\tmoatgen2             (moatgen with multiple constraints)\n");
            printf("\tmoat_concat  <file>  (concatenate 2 MOAT samples (2 INPUT SETS))\n");
         }
         else if (!strcmp(winput, "future"))
         {
            printf("\tgen_discrete (generate a sample of discrete variables)\n");
            printf("\tgp_sa2       (parameter screening with layered GP, not ready)\n");
            printf("\tsot_sa2      (parameter screening with sum-of-trees+boosting)\n");
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

      // ================================================================
      // run
      // ================================================================
      else if (!strcmp(command, "run"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("syntax: run <file>.\n");
            printf("where <file> is a PSUADE input file.\n");
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
         getInputFromFile(dataFile);
         run();
      }

      // ================================================================
      // Input/output commands
      // ================================================================
      // load data from a psuade file. 
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
            printf("This file may be overwritten during your interactive session.\n");
            printf("You can change the default ouptut filename with output_file.\n");
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
         if (psuadeIO_ != NULL) delete psuadeIO_;
         if (sampler_ != NULL) SamplingDestroy(sampler_);
         sampler_ = NULL;
         if (tagArray != NULL) delete [] tagArray;
         tagArray = NULL;
         psuadeIO_ = new PsuadeData();
         psuadeIO_->setOutputLevel(0);
         status = psuadeIO_->readPsuadeFile(dataFile);
         if (status == 0)
         {
            if (inputNames != NULL)
            {
               for (ii = 0; ii < nInputs; ii++)
                  if (inputNames[ii] != NULL) delete [] inputNames[ii];
               delete [] inputNames;
            }
            if (outputNames != NULL)
            {
               for (ii = 0; ii < nOutputs; ii++)
                  if (outputNames[ii] != NULL) delete [] outputNames[ii];
               delete [] outputNames;
            }
            psuadeIO_->getParameter("input_ninputs", pPtr);
            nInputs = pPtr.intData_;
            pINames.clean();
            psuadeIO_->getParameter("input_names", pINames);
            names = pINames.strArray_;
            pLower.clean();
            psuadeIO_->getParameter("input_lbounds", pLower);
            iLowerB = pLower.dbleArray_;
            pLower.dbleArray_ = NULL;
            pUpper.clean();
            psuadeIO_->getParameter("input_ubounds", pUpper);
            iUpperB = pUpper.dbleArray_;
            pUpper.dbleArray_ = NULL;
            inputNames = new char*[nInputs+1];
            for (ii = 0; ii < nInputs; ii++)
            {
               inputNames[ii] = new char[200]; 
               strcpy(inputNames[ii], names[ii]);
            }
            psuadeIO_->getParameter("output_noutputs", pPtr);
            nOutputs = pPtr.intData_;
            pONames.clean();
            psuadeIO_->getParameter("output_names", pONames);
            names = pONames.strArray_;
            outputNames = new char*[nOutputs+1];
            for (ii = 0; ii < nOutputs; ii++)
            {
               outputNames[ii] = new char[200]; 
               strcpy(outputNames[ii], names[ii]);
            }
            psuadeIO_->getParameter("method_sampling", pPtr);
            samplingMethod = pPtr.intData_;
            psuadeIO_->getParameter("method_nsamples", pPtr);
            nSamples = pPtr.intData_;
            psuadeIO_->getParameter("method_nreplications",pPtr);
            nReps = pPtr.intData_;
            psuadeIO_->getParameter("input_sample", pPtr);
            if (sampleInputs != NULL) delete [] sampleInputs;
            sampleInputs = pPtr.dbleArray_;
            pPtr.dbleArray_ = NULL;
            psuadeIO_->getParameter("output_sample", pPtr);
            if (sampleOutputs != NULL) delete [] sampleOutputs;
            sampleOutputs  = pPtr.dbleArray_;
            pPtr.dbleArray_ = NULL;
            psuadeIO_->getParameter("output_states", pPtr);
            if (sampleStates != NULL) delete [] sampleStates;
            sampleStates  = pPtr.intArray_;
            pPtr.intArray_ = NULL;
            pINames.clean();
            pONames.clean();
            if (sampleInputs == NULL || sampleOutputs == NULL)
            {
               printf("WARNING: no sample matrix nor output found.\n");
               nSamples = 0;
            }
            else
            {
               printf("load complete : nSamples = %d\n", nSamples);
               printf("                nInputs  = %d\n", nInputs);
               printf("                nOutputs = %d\n", nOutputs);
            }
         }
         else
         {
            cleanUp();
            if (psuadeIO_ != NULL) delete psuadeIO_;
            if (sampler_ != NULL) SamplingDestroy(sampler_);
            sampler_ = NULL;
            psuadeIO_ = NULL;
         }
         fflush(stdout);
         if (dataReg != NULL) delete [] dataReg;
         dataReg = NULL;
      }

      // load data from a psuade file and add to the existing data set
      // ----------------------------------------------------------------
      else if (!strcmp(command, "loadmore")) 
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         { 
            printf("loadmore: add another PSUADE sample to the existing sample\n");
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
            printf("This file may be overwritten during your interactive session.\n");
            printf("You can change the default ouptut filename with output_file.\n");
         }
         cleanUp();
         if (psuadeIO_ != NULL) delete psuadeIO_;
         if (sampler_ != NULL) SamplingDestroy(sampler_);
         sampler_ = NULL;
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
            if (sampleInputs != NULL)
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
               if (ind != nInputs)
               {
                  printf("ERROR: nInputs are different.\n");
                  printf("       incoming nInputs = %d\n", ind);
                  printf("       expected nInputs = %d\n", nInputs);
                  flag = 0;
               }
               if (flag == 1)
               {
                  psuadeIO_->getParameter("output_noutputs",pPtr);
                  ind = pPtr.intData_;
                  pONames.clean();
                  psuadeIO_->getParameter("output_names",pONames);
                  names = pONames.strArray_;
                  if (ind != nOutputs)
                  {
                     printf("ERROR: nOutputs are different.\n");
                     printf("       incoming nOutputs = %d\n", ind);
                     printf("       expected nOutputs = %d\n", nOutputs);
                     flag = 0;
                  }
               }
               if (flag == 1)
               {
                  for (ii = 0; ii < ind; ii++)
                  {
                     if (strcmp(outputNames[ii],names[ii]))
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
               if (flag == 0)
               {
                  if (sampleInputs  != NULL) delete [] sampleInputs;
                  if (sampleOutputs != NULL) delete [] sampleOutputs;
                  if (sampleStates  != NULL) delete [] sampleStates;
                  if (iLowerB       != NULL) delete [] iLowerB;
                  if (iUpperB       != NULL) delete [] iUpperB;
                  sampleInputs  = NULL;
                  sampleOutputs = NULL;
                  sampleStates  = NULL;
                  iLowerB       = NULL;
                  iUpperB       = NULL;
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
               if (inputNames != NULL)
               {
                  for (ii = 0; ii < nInputs; ii++)
                     if (inputNames[ii] != NULL) delete [] inputNames[ii];
                  delete [] inputNames;
               }
               if (outputNames != NULL)
               {
                  for (ii = 0; ii < nOutputs; ii++)
                     if (outputNames[ii] != NULL) delete [] outputNames[ii];
                  delete [] outputNames;
               }
               psuadeIO_->getParameter("input_ninputs", pPtr);
               nInputs = pPtr.intData_;
               pINames.clean();
               psuadeIO_->getParameter("input_names", pINames);
               names = pINames.strArray_;
               pLower.clean();
               psuadeIO_->getParameter("input_lbounds", pLower);
               iLowerB = pLower.dbleArray_;
               pLower.dbleArray_ = NULL;
               pUpper.clean();
               psuadeIO_->getParameter("input_ubounds", pUpper);
               iUpperB = pUpper.dbleArray_;
               pUpper.dbleArray_ = NULL;
               inputNames = new char*[nInputs+1];
               for (ii = 0; ii < nInputs; ii++)
               {
                  inputNames[ii] = new char[200]; 
                  strcpy(inputNames[ii], names[ii]);
               }
               psuadeIO_->getParameter("output_noutputs", pPtr);
               nOutputs = pPtr.intData_;
               pONames.clean();
               psuadeIO_->getParameter("output_names", pONames);
               names = pONames.strArray_;
               outputNames = new char*[nOutputs+1];
               for (ii = 0; ii < nOutputs; ii++)
               {
                  outputNames[ii] = new char[200]; 
                  strcpy(outputNames[ii], names[ii]);
               }
               psuadeIO_->getParameter("method_sampling", pPtr);
               samplingMethod = pPtr.intData_;
               psuadeIO_->getParameter("method_nsamples", pPtr);
               nSamples = pPtr.intData_;
               psuadeIO_->getParameter("method_nreplications",pPtr);
               nReps = pPtr.intData_;
               psuadeIO_->getParameter("input_sample", pPtr);
               sampleInputs  = pPtr.dbleArray_;
               pPtr.dbleArray_ = NULL;
               psuadeIO_->getParameter("output_sample", pPtr);
               sampleOutputs  = pPtr.dbleArray_;
               pPtr.dbleArray_ = NULL;
               psuadeIO_->getParameter("output_states", pPtr);
               sampleStates  = pPtr.intArray_;
               pPtr.intArray_ = NULL;
               pINames.clean();
               pONames.clean();
               printf("loadmore completed: nSamples = %d\n", nSamples);
            }
            if (flag == 1)
            {
               psuadeIO_->getParameter("method_nsamples", pPtr);
               ind = pPtr.intData_;
               tempX  = sampleInputs;
               tempY  = sampleOutputs;
               states = sampleStates;
               sampleInputs  = new double[(nSamples+ind)*nInputs];
               sampleOutputs = new double[(nSamples+ind)*nOutputs];
               sampleStates  = new int[nSamples+ind];
               for (ii = 0; ii < nSamples*nInputs; ii++)
                  sampleInputs[ii] = tempX[ii];
               for (ii = 0; ii < nSamples*nOutputs; ii++)
                  sampleOutputs[ii] = tempY[ii];
               for (ii = 0; ii < nSamples; ii++)
                  sampleStates[ii] = states[ii];
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
               for (ii = 0; ii < ind*nInputs; ii++)
                  sampleInputs[nSamples*nInputs+ii] = tempX[ii];
               for (ii = 0; ii < ind*nOutputs; ii++)
                  sampleOutputs[nSamples*nOutputs+ii] = tempY[ii];
               for (ii = 0; ii < ind; ii++)
                  sampleStates[nSamples+ii] = states[ii];
               nSamples += ind;
               psuadeIO_->updateInputSection(nSamples,nInputs,NULL,
                                             NULL,NULL,sampleInputs,NULL); 
               psuadeIO_->updateOutputSection(nSamples,nOutputs,
                                sampleOutputs, sampleStates, NULL); 
               printf("loadmore completed: nSamples = %d\n", nSamples);
            }
         }
         fflush(stdout);
         if (dataReg != NULL) delete [] dataReg;
         dataReg = NULL;
      }

      // write the data to an output file - this is needed if PDF ~=U
      // ----------------------------------------------------------------
      else if (!strcmp(command, "write"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         { 
            printf("write: save the current sample to a PSUADE data file\n");
            printf("syntax: write <filename>.\n");
            printf("where <filename> is a PSUADE data file.\n");
            continue;
         }
         strcpy(dataFile, "\0");
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: data not loaded - nothing to write\n");
         }
         else
         {
            strcpy(dataFile, psOutputFilename_);
            sscanf(lineIn, "%s %s", command, dataFile);
            if (!strcmp(dataFile,psInputFilename_))
            {
              printf("WARNING : output filename should not have the same name as the\n"); 
              printf(" input file %s. Try the output_file command.\n", psInputFilename_);
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
               if (nOutputs > 1)
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
                  sprintf(pString,"Enter output number (1 - %d) : ",nOutputs);
                  outputID = getInt(1, nOutputs, pString);
                  outputID--;
                  for (sInd = 0; sInd < nSamples; sInd++)
                     sampleOutputs[sInd] = 
                           sampleOutputs[sInd*nOutputs+outputID];
                  pONames.clean();
                  psuadeIO_->getParameter("output_names", pONames);
                  names = pONames.strArray_;
                  strcpy(names[0], names[outputID]);
                  nOutputs = 1;
                  psuadeIO_->updateAnalysisSection(0, 0, 0, 0, 0, 0);
               }
               else names = NULL;
               psuadeIO_->updateInputSection(nSamples,nInputs,NULL,NULL,
                                             NULL,sampleInputs,NULL); 
               psuadeIO_->updateOutputSection(nSamples,nOutputs,
                                sampleOutputs, sampleStates, names); 
               if (winput[0] == 'y') pONames.clean();
               psuadeIO_->updateMethodSection(samplingMethod,nSamples,nReps,
                                              -1,-1);
               psuadeIO_->writePsuadeFile(dataFile,0);
            }
         }
      }

      // write the unevaluated sample points data to an output file 
      else if (!strcmp(command, "nwrite"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("nwrite: save the unevaluated points to a PSUADE file\n");
            printf("syntax: nwrite <filename>.\n");
            printf("where <filename> is a PSUADE data file.\n");
            continue;
         }

         if (psuadeIO_ == NULL)
         {
            printf("ERROR: data not loaded - nothing to write\n");
         }
         else
         {
            strcpy(dataFile, psOutputFilename_);
            sscanf(lineIn, "%s %s", command, dataFile);
            if (!strcmp(dataFile,psInputFilename_))
            {
              printf("WARNING : output filename should not have the same name as the\n"); 
              printf(" input file %s. Try the output_file command.\n", psInputFilename_);
              printf(" Command not executed.\n");
            }
            else if (!strcmp(dataFile, "\0"))
            {
               printf("ERROR: need to specify a file to write to.\n");
               printf("syntax: nwrite <filename>.\n");
            }
            else
            {
               tempX  = new double[nInputs*nSamples];
               tempY  = new double[nOutputs*nSamples];
               states = new int[nSamples];
               count = 0;
               for (ii = 0; ii < nSamples; ii++)
               {
                  for (jj = 0; jj < nOutputs; jj++)
                     if (sampleOutputs[ii*nOutputs+jj] == PSUADE_UNDEFINED)
                        sampleStates[ii] = 0;
                  if (sampleStates[ii] == 0)
                  {
                     for (jj = 0; jj < nInputs; jj++)
                        tempX[count*nInputs+jj] = 
                           sampleInputs[ii*nInputs+jj];
                     for (jj = 0; jj < nOutputs; jj++)
                        tempY[count*nOutputs+jj] = 
                           sampleOutputs[ii*nOutputs+jj]; 
                     states[count++] = sampleStates[ii]; 
                  }
               }
               if (count == 0)
                  printf("INFO: no unevaluated sample points (no file generated).\n");
               else
               {
                  ioPtr = new PsuadeData();
                  ioPtr->updateInputSection(count, nInputs, NULL, iLowerB,
                                            iUpperB, tempX, inputNames);
                  ioPtr->updateOutputSection(count, nOutputs, tempY, states,
                                             outputNames);
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
      }

      // write the data to an output file in input only format 
      // ----------------------------------------------------------------
      else if (!strcmp(command, "iwrite"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("iwrite: save only the sample inputs to an ASCII file\n");
            printf("syntax: iwrite <filename>.\n");
            printf("where <filename> is the name of the target data file.\n");
            continue;
         }
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no data.\n");
            printf("       Use load first to put data into local memory.\n");
         }
         else
         {
            strcpy(dataFile, psOutputFilename_);
            sscanf(lineIn, "%s %s", command, dataFile);
            if (!strcmp(dataFile,psInputFilename_))
            {
              printf("WARNING : output filename should not have the same name as the\n"); 
              printf(" input file %s. Try the output_file command.\n", psInputFilename_);
              printf(" Command not executed.\n");
            }
            else if (!strcmp(dataFile, "\0"))
            {
               printf("ERROR: need to specify a file to write to.\n");
            }
            else
            {
               if (sampleInputs == NULL)
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
                  fprintf(fp, "%d %d\n", nSamples, nInputs);
                  for (ii = 0; ii < nSamples; ii++)
                  {
                     fprintf(fp, "%d ", ii+1);
                     for (jj = 0; jj < nInputs; jj++)
                        fprintf(fp,"%24.16e ",sampleInputs[ii*nInputs+jj]); 
                     fprintf(fp, "\n");
                  }
                  fprintf(fp, "PSUADE_END\n");
                  fclose(fp);
                  printf("iwrite: sample inputs written to %s.\n",dataFile);
               }
            }
         }
      }

      // write the data to an output file in output only format 
      // ----------------------------------------------------------------
      else if (!strcmp(command, "owrite"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("owrite: save only the sample outputs to an ASCII file\n");
            printf("syntax: owrite <filename>.\n");
            printf("where <filename> is the name of the target data file.\n");
            continue;
         }
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no data.\n");
            printf("       Use load first to put data into local memory.\n");
         }
         else
         {
            strcpy(dataFile, "\0");
            sscanf(lineIn, "%s %s", command, dataFile);
            strcpy(dataFile, psOutputFilename_);
            sscanf(lineIn, "%s %s", command, dataFile);
            if (!strcmp(dataFile,psInputFilename_))
            {
              printf("WARNING : output filename should not have the same name as the\n"); 
              printf(" input file %s. Try the output_file command.\n", psInputFilename_);
              printf(" Command not executed.\n");
            }
            else if (!strcmp(dataFile, "\0"))
            {
               printf("ERROR: need to specify a file to write to.\n");
            }
            else
            {
               if (sampleOutputs == NULL)
                  printf("ERROR: no output to write.\n");
               else
               {
                  if (nOutputs == 1) outputID = 0;
                  else
                  {
                     sprintf(pString, "Output number (1 - %d, %d for all): ",
                             nOutputs, nOutputs+1);
                     outputID = getInt(1, nOutputs+1, pString);
                     outputID--;
                  }
                  fp = fopen(dataFile, "w");
                  if (fp == NULL)
                  {
                     printf("ERROR: cannot open file %s.\n", dataFile);
                     continue;
                  }
                  fprintf(fp, "PSUADE_BEGIN\n");
                  if (outputID == nOutputs)
                       fprintf(fp, "%d %d\n", nSamples, nOutputs);
                  else fprintf(fp, "%d 1\n", nSamples);
                  for (ii = 0; ii < nSamples; ii++)
                  {
                     fprintf(fp, "%d ", ii+1);
                     if (outputID == nOutputs)
                     {
                        for (jj = 0; jj < nOutputs; jj++)
                           fprintf(fp,"%24.16e ",sampleOutputs[ii*nOutputs+jj]); 
                        fprintf(fp, "\n");
                     }
                     else
                        fprintf(fp,"%24.16e\n",
                                sampleOutputs[ii*nOutputs+outputID]);
                  }
                  fprintf(fp, "PSUADE_END\n");
                  fclose(fp);
                  printf("owrite: sample outputs written to %s.\n",dataFile);
               }
            }
         }
      }

      // read data in standard format 
      // ----------------------------------------------------------------
      else if (!strcmp(command, "read_std"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("read_std: read a data file in the standard format.\n");
            printf("line 1: nSamples nInputs nOutputs\n");
            printf("line 2: <sample point 1 inputs> < sample point 1 outputs>\n");
            printf("line 3: <sample point 2 inputs> < sample point 2 outputs>\n");
            printf(".......\n");
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
         fscanf(fp, "%d %d %d", &nSamples, &nInputs, &nOutputs);
         if (nSamples <= 0) 
         {
            printf("ERROR: nSamples <= 0 (%d)\n", nSamples);
            continue;
         }
         if (nInputs <= 0) 
         {
            printf("ERROR: nInputs <= 0 (%d)\n", nInputs);
            continue;
         }
         flag = 0;
         if (nOutputs <= 0) 
         {
            printf("INFP: nOutputs = 0 ==> create a dummy output\n");
            nOutputs = 1;
            flag = 1;
         }
         printf("nSamples = %d\n", nSamples);
         printf("nInputs  = %d\n", nInputs);
         printf("nOutputs = %d\n", nOutputs);
         sampleInputs = new double[nSamples*nInputs];
         sampleOutputs = new double[nSamples*nOutputs];
         sampleStates = new int[nSamples];
         if (sampler_ != NULL) SamplingDestroy(sampler_);
         sampler_ = NULL;
         for (ii = 0; ii < nSamples; ii++)
         {
            for (jj = 0; jj < nInputs; jj++)
               fscanf(fp, "%lg", &sampleInputs[ii*nInputs+jj]); 
            if (flag == 0)
            {
               for (jj = 0; jj < nOutputs; jj++)
                  fscanf(fp,"%lg", &sampleOutputs[ii*nOutputs+jj]);
            }
            else sampleOutputs[ii] = PSUADE_UNDEFINED;
            sampleStates[ii] = 1;
            for (jj = 0; jj < nOutputs; jj++)
               if (sampleOutputs[ii*nOutputs+jj] > 0.5*PSUADE_UNDEFINED)
                  sampleStates[ii] = 0;
         }
         fclose(fp);
         inputNames = new char*[nInputs+1];
         for (ii = 0; ii < nInputs; ii++)
         {
            inputNames[ii] = new char[200]; 
            sprintf(inputNames[ii], "X%d", ii+1);
         }
         outputNames = new char*[nOutputs+1];
         for (ii = 0; ii < nOutputs; ii++)
         {
            outputNames[ii] = new char[200]; 
            sprintf(outputNames[ii], "Y%d", ii+1);
         }
         if (iLowerB != NULL) delete [] iLowerB;
         iLowerB = new double[nInputs];
         if (iUpperB != NULL) delete [] iUpperB;
         iUpperB = new double[nInputs];
         for (jj = 0; jj < nInputs; jj++)
         {
            iLowerB[jj] = sampleInputs[jj];
            iUpperB[jj] = sampleInputs[jj];
            for (ii = 1; ii < nSamples; ii++)
            {
               if (sampleInputs[ii*nInputs+jj] < iLowerB[jj])
                  iLowerB[jj] = sampleInputs[ii*nInputs+jj];
               if (sampleInputs[ii*nInputs+jj] > iUpperB[jj])
                  iUpperB[jj] = sampleInputs[ii*nInputs+jj];
            }
         }
         samplingMethod = PSUADE_SAMP_MC;
         if (psuadeIO_ != NULL) delete psuadeIO_;
         psuadeIO_ = new PsuadeData();
         psuadeIO_->setOutputLevel(0);
         psuadeIO_->updateInputSection(nSamples,nInputs,NULL,iLowerB,
                                       iUpperB,sampleInputs,inputNames); 
         psuadeIO_->updateOutputSection(nSamples,nOutputs,sampleOutputs, 
                                        sampleStates, outputNames); 
         nReps = 1;
         psuadeIO_->updateMethodSection(samplingMethod,nSamples,nReps,0,0);
      }

      // write the data to an output file in standard format 
      // ----------------------------------------------------------------
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
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no data.\n");
            printf("       Use load to put data into local memory first.\n");
         }
         else
         {
            strcpy(dataFile, "\0");
            sscanf(lineIn, "%s %s", command, dataFile);
            if (!strcmp(dataFile,psInputFilename_))
            {
              printf("WARNING : output filename should not have the same name as the\n"); 
              printf(" input file %s.\n", psInputFilename_);
              printf(" Command not executed.\n");
            }
            else if (!strcmp(dataFile, "\0"))
            {
               printf("ERROR: need to specify a file to write to.\n");
               printf("syntax: write_std <filename>.\n");
            }
            else
            {
               if (sampleInputs == NULL)
                  printf("ERROR: no inputs to write.\n");
               else
               {
                  if (sampleOutputs == NULL) outputID = -9;
                  else
                  {
                     if (nOutputs == 1) outputID = 0;
                     else
                     {
                        sprintf(pString, "Enter output number (1 - %d, 0 for all) : ",
                                nOutputs);
                        outputID = getInt(0, nOutputs, pString);
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
                       fprintf(fp, "%d %d\n", nSamples, nInputs);
                  else if (outputID == -1) 
                       fprintf(fp, "%d %d %d\n", nSamples, nInputs, nOutputs);
                  else fprintf(fp, "%d %d 1\n", nSamples, nInputs);
                  for (ii = 0; ii < nSamples; ii++)
                  {
                     for (jj = 0; jj < nInputs; jj++)
                        fprintf(fp, "%24.16e ", sampleInputs[ii*nInputs+jj]); 
                     if (outputID == -9) fprintf(fp, "\n");
                     else if (outputID == -1)
                     {
                        for (jj = 0; jj < nOutputs; jj++)
                           fprintf(fp,"%24.16e ",sampleOutputs[ii*nOutputs+jj]);
                        fprintf(fp,"\n");
                     }
                     else fprintf(fp,"%24.16e\n",
                                  sampleOutputs[ii*nOutputs+outputID]);
                  }
                  fclose(fp);
               }
            }
         }
      }

      // write the data to an output file in matlab format 
      // ----------------------------------------------------------------
      else if (!strcmp(command, "write_matlab"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("write_matlab: write a data file in the matlab format.\n");
            printf("syntax: write_matlab <filename>.\n");
            printf("where <filename> is the name of the target data file.\n");
            continue;
         }
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no data.\n");
            printf("       Use load to put data into local memory first.\n");
         }
         else
         {
            strcpy(dataFile, "\0");
            sscanf(lineIn, "%s %s", command, dataFile);
            if (!strcmp(dataFile,psInputFilename_))
            {
              printf("WARNING : output filename should not have the same name as the\n"); 
              printf(" input file %s. \n", psInputFilename_);
              printf(" Command not executed.\n");
            }
            else if (!strcmp(dataFile, "\0"))
            {
               printf("ERROR: need to specify a file to write to.\n");
               printf("syntax: write_matlab <filename>.\n");
            }
            else
            {
               if (sampleInputs == NULL)
                  printf("ERROR: no inputs to write.\n");
               else
               {
                  if (sampleOutputs == NULL)
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
                  for (ii = 0; ii < nSamples; ii++)
                  {
                     for (jj = 0; jj < nInputs; jj++)
                        fprintf(fp, "%24.16e ", sampleInputs[ii*nInputs+jj]); 
                     fprintf(fp,"\n");
                  }
                  fprintf(fp,"];\n");
                  if (sampleOutputs != NULL)
                  {
                     fprintf(fp, "Y = [\n");
                     for (ii = 0; ii < nSamples; ii++)
                     {
                        for (jj = 0; jj < nOutputs; jj++)
                           fprintf(fp,"%24.16e ",sampleOutputs[ii*nOutputs+jj]);
                        fprintf(fp,"\n");
                     }
                     fprintf(fp,"];\n");
                  }
                  fclose(fp);
               }
            }
         }
      }

      // read the data in Excel format and convert to psuade form
      // ----------------------------------------------------------------
      else if (!strcmp(command, "read_xls"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("read_xls: read a data file in Excel format.\n");
            printf("syntax: read_xls <filename>.\n");
            printf("where <filename> is a data file in xls format:\n");
            printf("line 1: nSamples nInputs nOutputs\n");
            printf("line 2: 1 <sample point 1 inputs> < sample point 1 outputs>\n");
            printf("line 3: 2 <sample point 2 inputs> < sample point 2 outputs>\n");
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
         if (psuadeIO_ != NULL) delete psuadeIO_;
         if (sampler_  != NULL) SamplingDestroy(sampler_);
         sampler_ = NULL;
         psuadeIO_ = new PsuadeData();
         psuadeIO_->setOutputLevel(0);
         fp = fopen(dataFile,"r");
         fscanf(fp, "%d %d %d", &nSamples, &nInputs, &nOutputs);
         printf("nSamples = %d\n", nSamples);
         printf("nInputs  = %d\n", nInputs);
         printf("nOutputs = %d\n", nOutputs);
         if (nSamples <= 0 || nInputs <= 0 || nOutputs < 0)
         {
            printf("read_xls ERROR: some sample parameters <= 0.\n");
            printf("Note: the first line should be: \n");
            printf("  <nSamples> <nInputs> <nOutputs>\n");
            fclose(fp);
            continue;
         }
         flag = 0;
         if (nOutputs == 0)
         {
            flag = 1;
            nOutputs = 1;
         }
         sampleInputs = new double[nSamples*nInputs];
         sampleOutputs = new double[nSamples*nOutputs];
         sampleStates = new int[nSamples];
         for (ii = 0; ii < nSamples; ii++) sampleStates[ii] = 1;
         if (iLowerB != NULL) delete [] iLowerB;
         iLowerB = new double[nInputs];
         if (iUpperB != NULL) delete [] iUpperB;
         iUpperB = new double[nInputs];
         for (ii = 0; ii < nSamples; ii++)
         {
            fscanf(fp, "%d", &jj);
            if ((ii+1) != jj)
            {
               printf("read_xls ERROR: sample index mismatch.\n");
               printf("         sample number read     = %d\n", jj);
               printf("         sample number expected = %d\n", ii+1);
               printf("INFO: the first field in each line must be the sample number.\n");
               break;
            }
            for (jj = 0; jj < nInputs; jj++)
               fscanf(fp,"%lg", &sampleInputs[ii*nInputs+jj]); 
            if (flag == 0)
            {
               for (jj = 0; jj < nOutputs; jj++)
                  fscanf(fp,"%lg", &sampleOutputs[ii*nOutputs+jj]); 
            }
            else
            {
               sampleOutputs[ii] = PSUADE_UNDEFINED;
               sampleStates[ii] = 0;
            }
         }
         fclose(fp);
         if (ii != nSamples) continue;
         inputNames = new char*[nInputs];
         outputNames = new char*[nOutputs];
         tempW = new double[nSamples];
         for (ii = 0; ii < nInputs; ii++)
         {
            inputNames[ii] = new char[100]; 
            sprintf(inputNames[ii], "X%d", ii+1);
            for (jj = 0; jj < nSamples; jj++) 
               tempW[jj] = sampleInputs[jj*nInputs+ii];
            sortDbleList(nSamples, tempW);
            iLowerB[ii] = tempW[0];
            iUpperB[ii] = tempW[nSamples-1];
         }
         delete [] tempW;
         for (ii = 0; ii < nOutputs; ii++)
         {
            outputNames[ii] = new char[100]; 
            sprintf(outputNames[ii], "Y%d", ii+1);
         }
         psuadeIO_->updateInputSection(nSamples,nInputs,NULL,iLowerB,iUpperB,
                                       sampleInputs,inputNames); 
         psuadeIO_->updateOutputSection(nSamples, nOutputs, sampleOutputs, 
                                        sampleStates, outputNames);
         psuadeIO_->updateMethodSection(PSUADE_SAMP_MC, nSamples, 1, 0, 0);
         psuadeIO_->updateAnalysisSection(0, 0, 0, 0, 0, 0);
         nReps = 1;
         printf("Excel data have been read.\n");
         printf("Use 'write' to store the data in PSUADE format.\n");
      }

      // write the data to an output file in Excel format 
      // ----------------------------------------------------------------
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
            continue;
         }

         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no data to write.\n");
            printf("       Use load to put data into local memory first.\n");
         }
         else
         {
            strcpy(dataFile, "\0");
            sscanf(lineIn, "%s %s", command, dataFile);
            if (!strcmp(dataFile,psInputFilename_))
            {
              printf("WARNING : output filename should not have the same name as the\n"); 
              printf(" input file %s.\n", psInputFilename_);
              printf(" Command not executed.\n");
            }
            else if (!strcmp(dataFile, "\0"))
            {
               printf("ERROR: need to specify a file to write to.\n");
               printf("syntax: write_xls <filename>.\n");
            }
            else
            {
               if (sampleInputs == NULL || sampleOutputs == NULL)
                  printf("ERROR: no inputs or outputs to write.\n");
               else
               {
                  fp = fopen(dataFile, "w");
                  if (fp == NULL)
                  {
                     printf("ERROR: cannot open file %s.\n", dataFile);
                     continue;
                  }
                  fprintf(fp, "%d %d %d\n", nSamples, nInputs, nOutputs);
                  for (ii = 0; ii < nSamples; ii++)
                  {
                     fprintf(fp, "%d", ii+1);
                     for (jj = 0; jj < nInputs; jj++)
                        fprintf(fp,"\t%24.16e", sampleInputs[ii*nInputs+jj]); 
                     for (jj = 0; jj < nOutputs; jj++)
                        fprintf(fp,"\t%24.16e", sampleOutputs[ii*nOutputs+jj]);
                     fprintf(fp, "\n");
                  }
                  fclose(fp);
               }
            }
         }
      }

      // write the data to an output file in ULTRA format 
      // ----------------------------------------------------------------
      else if (!strcmp(command, "write_ultra"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("write_ultra: write to a data file in the ULTRA data format\n");
            printf("             (suitable for 2 input only).\n");
            printf("syntax: write_ultra <filename>.\n");
            printf("where <filename> is the name of the target data file.\n");
            continue;
         }
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no data.\n");
            printf("       Use load to put data into local memory first.\n");
         }
         else
         {
            strcpy(dataFile, "\0");
            sscanf(lineIn, "%s %s", command, dataFile);
            if (!strcmp(dataFile,psInputFilename_))
            {
              printf("WARNING : output filename should not have the same name as the\n"); 
              printf(" input file %s. \n", psInputFilename_);
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
               if (sampleInputs == NULL || sampleOutputs == NULL)
                  printf("ERROR: no inputs or outputs to save.\n");
               else if (nInputs != 2)
                  printf("ERROR: command available only for nInputs=2.\n");
               else
               {
                  if (nOutputs == 1) outputID = 0;
                  else
                  {
                     sprintf(pString, "Enter output number (1 - %d) : ",
                             nOutputs);
                     outputID = getInt(1, nOutputs, pString);
                     outputID--;
                  }
                  tempX = new double[nSamples];
                  tempW = new double[nSamples];
                  tempY = new double[nSamples];
                  tempT = new double[nSamples];
                  tempInds = new double[nSamples];
                  for (ii = 0; ii < nSamples; ii++)
                  {
                     tempX[ii] = sampleInputs[ii*nInputs]; 
                     tempW[ii] = sampleInputs[ii*nInputs+1]; 
                     tempY[ii] = sampleOutputs[ii*nOutputs+outputID];
                     tempInds[ii] = (double) ii;
                  }
                  sortDbleList2(nSamples, tempW, tempInds);
                  for (ii = 0; ii < nSamples; ii++) tempT[ii] = tempX[ii];
                  for (ii = 0; ii < nSamples; ii++)
                  {
                     ind2 = (int) tempInds[ii];
                     tempX[ii] = tempT[ind2];
                  }
                  for (ii = 0; ii < nSamples; ii++) tempT[ii] = tempY[ii];
                  for (ii = 0; ii < nSamples; ii++)
                  {
                     ind2 = (int) tempInds[ii];
                     tempY[ii] = tempT[ind2];
                  }
                  ind = 0;
                  for (ii = 1; ii < nSamples; ii++)
                  {
                     if (tempW[ii] != tempW[ii-1])
                     {
                        sortDbleList2(ii-ind, &tempX[ind], &tempY[ind]);
                        ind = ii;
                     }
                  } 
                  sortDbleList2(nSamples-ind, &tempX[ind], &tempY[ind]);
                  fp = fopen(dataFile, "w");
                  if (fp == NULL)
                  {
                     printf("ERROR: cannot open output file %s\n", dataFile);
                     exit(1);
                  }
                  ind = 0;
                  ind2 = 1;
                  for (ii = 1; ii < nSamples; ii++)
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
                  for (jj = ind; jj < nSamples; jj++)
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
      }

      // merge in the output of identical sample point from another file
      // ----------------------------------------------------------------
      else if (!strcmp(command, "update"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("update: update the current sample outputs from another\n");
            printf("        PSUADE data file (i.e. the current unevaluated\n");
            printf("        sample points will be replaced with evalauted\n");
            printf("        outputs from the provided data file.) New points\n");
            printf("        in the PSUADE data file will be ignored.\n");
            printf("syntax: update <filename> (merge from <filename>).\n");
            printf("where <filename> is a data file in PSUADE data format.\n");
            continue;
         }
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no data to update.\n");
            printf("       Use load to put data into local memory first.\n");
         }
         else
         {
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
            }
            if (status == 0)
            {
               ioPtr->getParameter("input_ninputs", pPtr);
               ind = pPtr.intData_;
               flag = 1;
               if (ind != nInputs) flag = 0;
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
                  if (ind != nOutputs) flag = 0;
               }
               if (flag == 1)
               {
                  pONames.clean();
                  ioPtr->getParameter("output_names", pONames);
                  names = pONames.strArray_;
                  for (ii = 0; ii < ind; ii++)
                     if (strcmp(outputNames[ii],names[ii])) flag = 0;
                  pONames.clean();
               }
               if (flag == 0)
                  printf("ERROR: invalid data file.\n");
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
                  for (sInd = 0; sInd < nSamples; sInd++)
                  {
                     if (sampleStates[sInd] == 1) continue;
                     kk = sInd * nInputs;
                     for (ii = 0; ii < ind; ii++)
                     {
                        if (states[ii] != 1) continue;
                        jj = ii * nInputs;
                        for (iInd = 0; iInd < nInputs; iInd++)
                           if (sampleInputs[kk+iInd] != tempX[jj+iInd])
                              break;
                        if (iInd == nInputs)
                        {
                           kk = sInd * nOutputs;
                           jj = ii * nOutputs;
                           for (oInd = 0; oInd < nOutputs; oInd++)
                              sampleOutputs[kk+oInd] = tempY[jj+oInd];
                           sampleStates[sInd] = 1;
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
            fflush(stdout);
            if (ioPtr != NULL) delete ioPtr;
            ioPtr = NULL;
         }
      }

      // repair an MOAT sample file
      // ----------------------------------------------------------------
      else if (!strcmp(command, "moatrepair"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("moatrepair: repair a Morris sample using new input values\n");
            printf("syntax: moatrepair <filename>.\n");
            printf("where <filename> is a MOAT repair file (created by moatgen).\n");
            continue;
         }
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no Morris sample to repair.\n");
            printf("       Use load to put data into local memory first.\n");
         }
         else
         {
            strcpy(dataFile, "\0");
            sscanf(lineIn, "%s %s", command, dataFile);
            status = 0;
            fp = fopen(dataFile, "r");
            if (fp == NULL)
            {
               printf("ERROR: File %s not found.\n", dataFile);
               printf("syntax: moatrepair <file>.\n");
               printf("where <file> is a MOAT repair file (created by moatgen).\n");
               status = 1;
            }
            else fclose(fp);

            if (samplingMethod != PSUADE_SAMP_MOAT)
            {
               printf("Cannot repair: the sampling method is not MOAT.\n");
               status = 1;
            }
            if (status == 0)
            {
               sampPtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_MOAT);
               sampPtr->setPrintLevel(0);
               sampPtr->setInputBounds(nInputs, iLowerB, iUpperB);
               sampPtr->setOutputParams(nOutputs);
               sampPtr->setSamplingParams(nSamples,nSamples/(nInputs+1),0);
               sampPtr->initialize(1);
               sampPtr->loadSamples(nSamples, nInputs, nOutputs, sampleInputs,
                                    sampleOutputs, sampleStates);
               sampPtr->repair(dataFile, 0);
               sampPtr->getSamples(nSamples, nInputs, nOutputs, sampleInputs,
                                   sampleOutputs, sampleStates);
               delete sampPtr;
               sampPtr = NULL;
            }
            fflush(stdout);
         }
      }

      // repair an GMOAT sample file
      // ----------------------------------------------------------------
      else if (!strcmp(command, "gmoatrepair"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("gmoatrepair: repair a Morris sample using new input values\n");
            printf("Note: gmoatrepair differs from moatrepair in that it uses\n");
            printf("      GMOAT instead of MOAT to create the initial sample.\n");
            printf("syntax: gmoatrepair <filename>.\n");
            printf("where <filename> is a MOAT repair file (created by moatgen).\n");
            continue;
         }
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no Morris sample to repair.\n");
            printf("       Use load to put data into local memory first.\n");
         }
         else
         {
            strcpy(dataFile, "\0");
            sscanf(lineIn, "%s %s", command, dataFile);
            status = 0;
            fp = fopen(dataFile, "r");
            if (fp == NULL)
            {
               printf("ERROR: File %s not found.\n", dataFile);
               printf("syntax: gmoatrepair <file>.\n");
               printf("where <file> is a MOAT repair file (created by moatgen).\n");
               status = 1;
            }
            else fclose(fp);

            if (samplingMethod != PSUADE_SAMP_GMOAT)
            {
               printf("Cannot repair: the sampling method is not GMOAT.\n");
               status = 1;
            }
            if (status == 0)
            {
               sampPtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_GMOAT);
               sampPtr->setPrintLevel(0);
               sampPtr->setInputBounds(nInputs, iLowerB, iUpperB);
               sampPtr->setOutputParams(nOutputs);
               sampPtr->setSamplingParams(nSamples,nSamples/(nInputs+1),0);
               sampPtr->initialize(1);
               sampPtr->loadSamples(nSamples, nInputs, nOutputs, sampleInputs,
                                    sampleOutputs, sampleStates);
               sampPtr->repair(dataFile, 0);
               sampPtr->getSamples(nSamples, nInputs, nOutputs, sampleInputs,
                                   sampleOutputs, sampleStates);
               delete sampPtr;
               sampPtr = NULL;
            }
            fflush(stdout);
         }
      }

      // put additional inputs from another file which has the same
      // number of sample points
      // ----------------------------------------------------------------
      else if (!strcmp(command, "iadd"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("iadd: add more inputs to the existing sample from another\n");
            printf("      PSUADE data file\n");
            printf("       syntax: iadd <filename> (<filename> contains new inputs)\n");
            printf("       where <filename> should be a PSUADE data file.\n");
            continue;
         }
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no data to add to.\n");
            printf("       Use load to put data into local memory first.\n");
         }
         else
         {
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
            }
            if (status == 0)
            {
               ioPtr->getParameter("method_nsamples", pPtr);
               ind = pPtr.intData_;
               if (ind != nSamples)
               {
                  printf("ERROR: nSamples not the same in both sets of data.\n");
                  printf("       nSamples in local memory = %d\n", nSamples);
                  printf("       nSamples from file       = %d\n", ind);
               }
               else
               {
                  ioPtr->getParameter("input_ninputs", pPtr);
                  ind = pPtr.intData_;

                  tempW = iLowerB;
                  ioPtr->getParameter("input_lbounds", pLower);
                  tempT = pLower.dbleArray_;
                  iLowerB = new double[nInputs+ind];
                  for (ii = 0; ii < nInputs; ii++)
                     iLowerB[ii] = tempW[ii];
                  for (ii = 0; ii < ind; ii++)
                     iLowerB[nInputs+ii] = tempT[ii];
                  pLower.clean();
                  delete [] tempW;
                  
                  tempW = iUpperB;
                  ioPtr->getParameter("input_ubounds", pUpper);
                  tempT = pUpper.dbleArray_;
                  iUpperB = new double[nInputs+ind];
                  for (ii = 0; ii < nInputs; ii++)
                     iUpperB[ii] = tempW[ii];
                  for (ii = 0; ii < ind; ii++)
                     iUpperB[nInputs+ii] = tempT[ii];
                  pUpper.clean();
                  delete [] tempW;

                  names = inputNames;
                  inputNames = new char*[nInputs+ind];
                  for (ii = 0; ii < nInputs; ii++)
                     inputNames[ii] = names[ii];
                  if (names != NULL) delete [] names;
                  pINames.clean();
                  ioPtr->getParameter("input_names", pINames);
                  names = pINames.strArray_;
                  for (ii = 0; ii < ind; ii++)
                  {
                     inputNames[nInputs+ii] = names[ii];
                     names[ii] = NULL;
                  }
                  pINames.clean();
                  ioPtr->getParameter("input_sample", pInputs);
                  tempW = pInputs.dbleArray_;
                  tempX = sampleInputs;
                  sampleInputs = new double[nSamples*(nInputs+ind)];
                  kk = 0;
                  for (ii = 0; ii < nSamples; ii++)
                  {
                     for (jj = 0; jj < nInputs; jj++)
                        sampleInputs[kk++] = tempX[ii*nInputs+jj];
                     for (jj = 0; jj < ind; jj++)
                        sampleInputs[kk++] = tempW[ii*ind+jj];
                  }
                  delete [] tempX;
                  nInputs += ind;
                  pInputs.clean();
                  psuadeIO_->updateInputSection(nSamples,nInputs,NULL,
                                 iLowerB,iUpperB,sampleInputs,inputNames); 
               }
            }
            fflush(stdout);
            if (ioPtr != NULL) delete ioPtr;
            ioPtr = NULL;
         }
      }

      // put additional outputs from another file which has the same
      // sample points
      // ----------------------------------------------------------------
      else if (!strcmp(command, "oadd"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("oadd: add more outputs to the existing sample from another\n");
            printf("      PSUADE data file\n");
            printf("      syntax: oadd <filename> (<filename> contains new outputs)\n");
            printf("      where <filename> should be a PSUADE data file.\n");
            continue;
         }
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no data to add to.\n");
            printf("       Use load to put data into local memory first.\n");
         }
         else
         {
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
            }
            if (status == 0)
            {
               ioPtr->getParameter("method_nsamples", pPtr);
               ind = pPtr.intData_;
               if (ind != nSamples)
               {
                  printf("ERROR: nSamples not the same in both sets of data.\n");
                  printf("       nSamples in local memory = %d\n", nSamples);
                  printf("       nSamples from file       = %d\n", ind);
               }
               else
               {
                  printf("INFO: no checking for whether both sets of inputs match.\n");
                  ioPtr->getParameter("output_noutputs", pPtr);
                  ind = pPtr.intData_;

                  names = outputNames;
                  outputNames = new char*[nOutputs+ind];
                  for (ii = 0; ii < nOutputs; ii++)
                     outputNames[ii] = names[ii];
                  if (names != NULL) delete [] names;
                  pONames.clean();
                  ioPtr->getParameter("output_names", pONames);
                  names = pONames.strArray_;
                  for (ii = 0; ii < ind; ii++)
                  {
                     outputNames[nOutputs+ii] = names[ii];
                     names[ii] = NULL;
                  }
                  pONames.clean();
                  ioPtr->getParameter("output_sample", pOutputs);
                  tempW = pOutputs.dbleArray_;
                  tempX = sampleOutputs;
                  sampleOutputs = new double[nSamples*(nOutputs+ind)];
                  kk = 0;
                  for (ii = 0; ii < nSamples; ii++)
                  {
                     for (jj = 0; jj < nOutputs; jj++)
                        sampleOutputs[kk++] = tempX[ii*nOutputs+jj];
                     for (jj = 0; jj < ind; jj++)
                        sampleOutputs[kk++] = tempW[ii*ind+jj];
                  }
                  delete [] tempX;
                  nOutputs += ind;
                  pOutputs.clean();
                  psuadeIO_->updateOutputSection(nSamples,nOutputs,
                                   sampleOutputs, sampleStates, outputNames); 
               }
            }
            fflush(stdout);
            if (ioPtr != NULL) delete ioPtr;
            ioPtr = NULL;
         }
      }

      // replace one existing input with an input from a given file
      // ----------------------------------------------------------------
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
            printf("line 2: nSamples\n");
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
         }
         else
         {
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
               if (kk != nSamples)
               {
                  fclose(fp);
                  printf("ERROR: File and local data do not match.\n");
                  printf("       nSamples for sample in local memory = %d\n",nSamples);
                  printf("       nSamples from external file         = %d\n",kk);
                  status = 1;
               }
            }
            if (status == 0)
            {
               inputID = 0;
               if (nInputs > 1)
               {
                  sprintf(pString,"Enter input number to replace (1 - %d) : ",nInputs);
                  inputID = getInt(1, nInputs, pString);
                  inputID--;
               }
               for (ii = 0; ii < nSamples; ii++)
                  fscanf(fp, "%lg", &(sampleInputs[ii*nInputs+inputID]));
               Xmin = PSUADE_UNDEFINED;
               Xmax = - Xmin;
               for (ii = 0; ii < nSamples; ii++)
               {
                  if (sampleInputs[ii*nInputs+inputID] < Xmin)
                     Xmin = sampleInputs[ii*nInputs+inputID];
                  else if (sampleInputs[ii*nInputs+inputID] > Xmax)
                     Xmax = sampleInputs[ii*nInputs+inputID];
               }
               iLowerB[inputID] = Xmin; 
               iUpperB[inputID] = Xmax; 
               sprintf(inputNames[inputID], "X%d", inputID+1);
               psuadeIO_->updateInputSection(nSamples,nInputs,NULL,iLowerB,
                             iUpperB,sampleInputs,inputNames); 
               fscanf(fp, "%s", winput);
               if (strcmp(winput, "PSUADE_END"))
                  printf("WARNING: file should end with PSUADE_END\n");
               fclose(fp);
            }
            fflush(stdout);
         }
      }

      // replace all existing output with outputs from a given file
      // ----------------------------------------------------------------
      else if (!strcmp(command, "oreplace"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("oreplace: replace all outputs in the existing sample from\n");
            printf("          another sample file in ASCII format.\n");
            printf("syntax: oreplace <file>\n");
            printf("where <file> should be in the format given below: \n");
            printf("line 1: PSUADE_BEGIN\n");
            printf("line 2: nSamples nOutputs\n");
            printf("line 3: output value for sample point 1\n");
            printf("line 4: output value for sample point 2\n");
            printf(".......\n");
            printf("last line: PSUADE_END\n");
            continue;
         }
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no data to oreplace.\n");
            printf("       Use load to put data into local memory first.\n");
         }
         else
         {
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
               if (kk != nSamples || nOut2 <= 0)
               {
                  printf("ERROR: File and local parameters do not match.\n");
                  printf("       nSamples (local) = %d\n", nSamples);
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
               printf("       nSamples (local) = %d\n", nSamples);
               printf("       nOutputs (file)  = %d\n", nOut2);
               if (nOutputs != nOut2)
               {
                  if (outputNames != NULL) 
                  {
                     for (ii = 0; ii < nOutputs; ii++)
                        delete [] outputNames[ii];
                     delete [] outputNames;
                  }
                  nOutputs = nOut2;
                  outputNames = new char*[nOutputs];
                  for (ii = 0; ii < nOutputs; ii++)
                  {
                     outputNames[ii] = new char[200];
                     sprintf(outputNames[ii], "Y%d", ii+1);
                  }
                  if (sampleOutputs != NULL) delete [] sampleOutputs;
                  sampleOutputs = NULL;
               }
               if (sampleOutputs == NULL)
                  sampleOutputs = new double[nSamples*nOutputs]; 
               for (ii = 0; ii < nSamples*nOutputs; ii++)
                  fscanf(fp, "%lg", &(sampleOutputs[ii]));
               for (ii = 0; ii < nSamples; ii++) sampleStates[ii] = 1;
               psuadeIO_->updateOutputSection(nSamples,nOutputs,
                             sampleOutputs,sampleStates,outputNames); 
               fscanf(fp, "%s", winput);
               if (strcmp(winput, "PSUADE_END"))
                  printf("WARNING: file should end with PSUADE_END\n");
            }
            fflush(stdout);
         }
      }

      // generate MOAT repair
      // ----------------------------------------------------------------
      else if (!strcmp(command, "moatgen") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("moatgen: create a Morris repair file with one constraint\n");
            printf("syntax: moatgen\n");
            printf("Note: a PSUADE datafile should have been loaded before\n");
            printf("      using this command. The data will be used as response\n");
            printf("      surface to find feasible region in creating a MOAT\n");
            printf("      repair sample.\n");
         }
      }
      else if (!strcmp(command, "moatgen") && nInputs >= 1)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("moatgen: create a Morris repair file with one constraint\n");
            printf("syntax: moatgen\n");
            printf("Note: a PSUADE datafile should have been loaded before\n");
            printf("      using this command. The data will be used as response\n");
            printf("      surface to find feasible region in creating a MOAT\n");
            printf("      repair sample.\n");
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
         for (sInd = 0; sInd < nSamples; sInd++)
         {
            if (sampleOutputs[sInd*nOutputs+outputID] > Ymax)
               Ymax = sampleOutputs[sInd*nOutputs+outputID];
            if (sampleOutputs[sInd*nOutputs+outputID] < Ymin)
               Ymin = sampleOutputs[sInd*nOutputs+outputID];
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
         moatSample = new double*[nPaths*(nInputs+1)];
         for (ii = 0; ii < nPaths*(nInputs+1); ii++)
            moatSample[ii] = new double[nInputs];
         tempW = new double[nInputs];
         indSet = new int[nInputs];
         count = 0;
         filterRange = threshU - threshL;
         for (ii = 0; ii < nPaths; ii++)
         {
            trial = 0; 
            while (trial < nTrials)
            {
               iInd = count;
               trial++;
               for (jj = 0; jj < nInputs; jj++)
               {
                  ind = PSUADE_rand() % currP;
                  dtemp = ind * (iUpperB[jj] - iLowerB[jj]) / (currP - 1.0);
                  moatSample[iInd][jj] = iLowerB[jj] + dtemp;
               }
               dtemp = faPtr->evaluatePoint(moatSample[iInd]);
               if (dtemp < threshL || dtemp > threshU) continue;

               generateRandomIvector(nInputs, indSet);

               for (jj = 0; jj < nInputs; jj++)
               {
                  iInd++;
                  for (kk = 0; kk < nInputs; kk++)
                     moatSample[iInd][kk] = moatSample[iInd-1][kk];

                  ind2 = indSet[jj];

                  ddata = moatSample[iInd][ind2]; 
                  moatSample[iInd][ind2] = iLowerB[ind2];
                  currY1 = faPtr->evaluatePoint(moatSample[iInd]);
                  moatSample[iInd][ind2] = iUpperB[ind2];
                  currY2 = faPtr->evaluatePoint(moatSample[iInd]);
                  currX1 = iLowerB[ind2];
                  currX2 = iUpperB[ind2];
                  moatSample[iInd][ind2] = ddata;

                  if (currY2 >= threshU && currY1 >= threshU)
                     currX1 = currX2 = 0.0;
                  else if (currY2 <= threshL && currY1 <= threshL)
                     currX1 = currX2 = 0.0;
                  else if (currY2 > currY1)
                  {
                     if (currY2 <= threshU) currX2 = iUpperB[ind2];
                     else
                     {
                        sLo = iLowerB[ind2];
                        sHi = iUpperB[ind2];
                        while (PABS((currY2-threshU)/filterRange)>0.0001)
                        {
                           moatSample[iInd][ind2] = 0.5 * (sLo + sHi);
                           currY2 = faPtr->evaluatePoint(moatSample[iInd]);
                           if (currY2 > threshU) sHi = 0.5 * (sLo + sHi);
                           else                  sLo = 0.5 * (sLo + sHi);
                        }
                        currX2 = moatSample[iInd][ind2];
                     }
                     if (currY1 >= threshL) currX1 = iLowerB[ind2];
                     else
                     {
                        sLo = iLowerB[ind2];
                        sHi = iUpperB[ind2];
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
                     if (currY1 <= threshU) currX1 = iLowerB[ind2];
                     else
                     {
                        sLo = iLowerB[ind2];
                        sHi = iUpperB[ind2];
                        while (PABS((currY1-threshU)/filterRange)>0.0001)
                        {
                           moatSample[iInd][ind2] = 0.5 * (sLo + sHi);
                           currY1 = faPtr->evaluatePoint(moatSample[iInd]);
                           if (currY1 > threshU) sLo = 0.5 * (sLo + sHi);
                           else                  sHi = 0.5 * (sLo + sHi);
                        }
                        currX1 = moatSample[iInd][ind2];
                     }
                     if (currY2 >= threshL) currX2 = iUpperB[ind2];
                     else
                     {
                        sLo = iLowerB[ind2];
                        sHi = iUpperB[ind2];
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
                  if (PABS(currX2-currX1)<0.1*(iUpperB[ind2]-iLowerB[ind2])) 
                     break;
                  moatSample[iInd][ind2] = ddata;
                  tempW[ind2] = PABS(currX2 - currX1) / (currP-1.0);
                  moatSample[iInd][ind2] += tempW[ind2];
                  if (moatSample[iInd][ind2] > iUpperB[ind2])
                       dtemp = threshL - 1.0;
                  else dtemp = faPtr->evaluatePoint(moatSample[iInd]);
                  if (dtemp < threshL || dtemp > threshU) 
                  {
                     moatSample[iInd][ind2] -= 2.0 * tempW[ind2];
                     if (moatSample[iInd][ind2] < iLowerB[ind2])
                        break;
                     dtemp = faPtr->evaluatePoint(moatSample[iInd]);
                     if (dtemp < threshL || dtemp > threshU) 
                        break;
                  }
               }
               if (jj == nInputs) 
               {
                  count += (nInputs + 1);
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
            for (ii = 0; ii < nPaths*(nInputs+1); ii++)
               delete [] moatSample[ii];
            delete [] moatSample;
            continue;
         }
         fp = fopen("MOAT_repair_file", "w");
         if (fp == NULL)
         {
            printf("ERROR: cannot open file MOAT_repair_file.\n");
            continue;
         }
         fprintf(fp, "BEGIN\n");
         fprintf(fp, "%d %d\n", nPaths*(nInputs+1), nInputs);
         for (ii = 0; ii < nInputs; ii++)
            fprintf(fp, "%d ", ii+1);
         fprintf(fp, "\n");
         for (ii = 0; ii < nPaths*(nInputs+1); ii++)
         {
            for (jj = 0; jj < nInputs; jj++)
               fprintf(fp, "%e ", moatSample[ii][jj]);
            fprintf(fp, "\n");
         }
         fprintf(fp, "END\n");
         fclose(fp);
         count = nPaths * (nInputs + 1); 
         tempW = new double[count*nInputs];
         states = new int[count];
         for (ii = 0; ii < count; ii++) states[ii] = 1;
         for (ii = 0; ii < count; ii++)
            for (jj = 0; jj < nInputs; jj++)
               tempW[ii*nInputs+jj] = moatSample[ii][jj];
         printf("moatgen: check for repeated sample points.\n");
         for (ii = 0; ii < count; ii++)
         {
            status = compareSamples(ii,count,nInputs, tempW, states);
            if (status >= 0)
               printf("moatgen check: sample %d and %d are identical.\n",
                      ii+1,status+1);
         }
         printf("moatgen: repair file created in MOAT_repair_file.\n");
         printf("         Make sure to change the input indices to match\n");
         printf("         the indices in the original MOAT file when used\n");
         printf("         with gmoatrepair.\n");
         for (ii = 0; ii < nPaths*(nInputs+1); ii++) delete [] moatSample[ii];
         delete [] moatSample;
         delete [] tempW;
         delete [] states;
         delete faPtr;
         faPtr = NULL;
      }

      // generate MOAT repair file with multiple constraints
      // ----------------------------------------------------------------
      else if (!strcmp(command, "moatgen2") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("moatgen2: create a Morris repair file with multiple constraints\n");
            printf("syntax: moatgen\n");
            printf("Note: a PSUADE datafile should have been loaded before\n");
            printf("      using this command. The data will be used as response\n");
            printf("      surface to find feasible region in creating a MOAT\n");
            printf("      repair sample.\n");
            continue;
         }
         printf("ERROR: no input file loaded.\n");
      }
      else if (!strcmp(command, "moatgen2"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("moatgen2: create a Morris repair file with multiple constraints\n");
            printf("syntax: moatgen\n");
            printf("Note: a PSUADE datafile should have been loaded before\n");
            printf("      using this command. The data will be used as response\n");
            printf("      surface to find feasible region in creating a MOAT\n");
            printf("      repair sample.\n");
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
               exit(1);
            }
            ioPtr->getParameter("input_ninputs", pPtr);
            jj = pPtr.intData_;
            if (jj != nInputs)
            {
               printf("moatgen2 ERROR: nInputs mismatch.\n");
               printf("         local nInputs = %d.\n", nInputs);
               printf("         file  nInputs = %d.\n", jj);
               exit(1);
            }
            pLower.clean();
            ioPtr->getParameter("input_lbounds", pLower);
            for (ii = 0; ii < nInputs; ii++)
            {
               if (iLowerB[ii] != pLower.dbleArray_[ii])
               {
                  printf("moatgen2 ERROR: lower bound mismatch (input %d)\n", ii+1);
                  exit(1);
               }
            }
            pUpper.clean();
            ioPtr->getParameter("input_ubounds", pUpper);
            for (ii = 0; ii < nInputs; ii++)
            {
               if (iUpperB[ii] != pUpper.dbleArray_[ii])
               {
                  printf("moatgen2 ERROR: upper bound mismatch (input %d)\n", ii+1);
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
         nPaths = 5000;
         sprintf(pString, "Enter P (resolution: try 4-10) : ");
         currP = getInt(4, 10, pString);
         sprintf(pString, "Enter the number of trials (> 100) : ");
         nTrials = getInt(101, 10000000, pString);
         moatSample = new double*[nPaths*(nInputs+1)];
         for (ii = 0; ii < nPaths*(nInputs+1); ii++)
            moatSample[ii] = new double[nInputs];
         tempW = new double[nInputs];
         indSet = new int[nInputs];
         count = 0;
         for (ii = 0; ii < nPaths; ii++)
         {
            trial = 0; 
            while (trial < nTrials)
            {
               iInd = count;
               trial++;
               for (jj = 0; jj < nInputs; jj++)
               {
                  ind = PSUADE_rand() % currP;
                  dtemp = ind * (iUpperB[jj] - iLowerB[jj]) / (currP - 1.0);
                  moatSample[iInd][jj] = iLowerB[jj] + dtemp;
               }
               for (kk = 0; kk < nFiles; kk++)
               {
                  dtemp = faPtrs[kk]->evaluatePoint(moatSample[iInd]);
                  if (dtemp < threshLs[ii] || dtemp > threshUs[ii]) break;
               }
               if (kk != nFiles) continue;

               generateRandomIvector(nInputs, indSet);

               for (jj = 0; jj < nInputs; jj++)
               {
                  iInd++;
                  for (kk = 0; kk < nInputs; kk++)
                     moatSample[iInd][kk] = moatSample[iInd-1][kk];

                  ind2 = indSet[jj];

                  ddata = moatSample[iInd][ind2]; 
                  currX1F = - PSUADE_UNDEFINED;
                  currX2F =   PSUADE_UNDEFINED;
                  for (kk = 0; kk < nFiles; kk++)
                  {
                     filterRange = threshUs[kk] - threshLs[kk];
                     moatSample[iInd][ind2] = iLowerB[ind2];
                     currY1 = faPtrs[kk]->evaluatePoint(moatSample[iInd]);
                     moatSample[iInd][ind2] = iUpperB[ind2];
                     currY2 = faPtrs[kk]->evaluatePoint(moatSample[iInd]);
                     currX1 = iLowerB[ind2];
                     currX2 = iUpperB[ind2];

                     if (currY2 >= threshUs[kk] && currY1 >= threshUs[kk])
                        currX1 = currX2 = 0.0;
                     else if (currY2 <= threshLs[kk] && currY1 <= threshLs[kk])
                        currX1 = currX2 = 0.0;
                     else if (currY2 > currY1)
                     {
                        if (currY2 <= threshUs[kk]) currX2 = iUpperB[ind2];
                        else
                        {
                           sLo = iLowerB[ind2];
                           sHi = iUpperB[ind2];
                           while (PABS((currY2-threshUs[kk])/filterRange)>1e-4)
                           {
                              moatSample[iInd][ind2] = 0.5 * (sLo + sHi);
                              currY2=faPtrs[kk]->evaluatePoint(moatSample[iInd]);
                              if (currY2 > threshUs[kk]) sHi = 0.5 * (sLo+sHi);
                              else                       sLo = 0.5 * (sLo+sHi);
                           }
                           currX2 = moatSample[iInd][ind2];
                        }
                        if (currY1 >= threshLs[kk]) currX1 = iLowerB[ind2];
                        else
                        {
                           sLo = iLowerB[ind2];
                           sHi = iUpperB[ind2];
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
                        if (currY1 <= threshUs[kk]) currX1 = iLowerB[ind2];
                        else
                        {
                           sLo = iLowerB[ind2];
                           sHi = iUpperB[ind2];
                           while (PABS((currY1-threshUs[kk])/filterRange)>1e-4)
                           {
                              moatSample[iInd][ind2] = 0.5 * (sLo + sHi);
                              currY1=faPtrs[kk]->evaluatePoint(moatSample[iInd]);
                              if (currY1 > threshUs[kk]) sLo = 0.5 * (sLo+sHi);
                              else                       sHi = 0.5 * (sLo+sHi);
                           }
                           currX1 = moatSample[iInd][ind2];
                        }
                        if (currY2 >= threshLs[kk]) currX2 = iUpperB[ind2];
                        else
                        {
                           sLo = iLowerB[ind2];
                           sHi = iUpperB[ind2];
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
                     if (PABS(currX2-currX1)<0.1*(iUpperB[ind2]-iLowerB[ind2])) 
                        break;
                     if (currX1 > currX1F) currX1F = currX1;
                     if (currX2 < currX2F) currX2F = currX2;
                  }
                  moatSample[iInd][ind2] = ddata;
                  if (kk != nFiles) break;
                  tempW[ind2] = PABS(currX2F - currX1F) / (currP-1.0);
                  moatSample[iInd][ind2] += tempW[ind2];
                  if (moatSample[iInd][ind2] > iUpperB[ind2])
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
                           if (moatSample[iInd][ind2] < iLowerB[ind2])
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
               if (jj == nInputs) 
               {
                  count += (nInputs + 1);
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
            for (ii = 0; ii < nPaths*(nInputs+1); ii++)
               delete [] moatSample[ii];
            delete [] moatSample;
            continue;
         }
         fp = fopen("MOAT_repair_file", "w");
         if (fp == NULL)
         {
            printf("ERROR: cannot open file MOAT_repair_file.\n");
            continue;
         }
         fprintf(fp, "BEGIN\n");
         fprintf(fp, "%d %d\n", nPaths*(nInputs+1), nInputs);
         for (ii = 0; ii < nInputs; ii++)
            fprintf(fp, "%d ", ii+1);
         fprintf(fp, "\n");
         for (ii = 0; ii < nPaths*(nInputs+1); ii++)
         {
            for (jj = 0; jj < nInputs; jj++)
               fprintf(fp, "%e ", moatSample[ii][jj]);
            fprintf(fp, "\n");
         }
         fprintf(fp, "END\n");
         fclose(fp);
         count = nPaths * (nInputs + 1); 
         tempW = new double[count*nInputs];
         states = new int[count];
         for (ii = 0; ii < count; ii++) states[ii] = 1;
         for (ii = 0; ii < count; ii++)
            for (jj = 0; jj < nInputs; jj++)
               tempW[ii*nInputs+jj] = moatSample[ii][jj];
         printf("moatgen2: check for repeated sample points.\n");
         for (ii = 0; ii < count; ii++)
         {
            status = compareSamples(ii,count,nInputs, tempW, states);
            if (status >= 0)
               printf("moatgen2 check: sample %d and %d are identical.\n",
                      ii+1,status+1);
         }
         printf("moatgen2: repair file created in MOAT_repair_file.\n");
         printf("          Make sure to change the input indices to match\n");
         printf("          the indices in the original MOAT file when used\n");
         printf("          with gmoatrepair.\n");
         for (ii = 0; ii < nPaths*(nInputs+1); ii++) delete [] moatSample[ii];
         delete [] moatSample;
         delete [] tempW;
         delete [] states;
      }

      // concatenate 2 MOAT samples with 2 different input sets
      // ----------------------------------------------------------------
      else if (!strcmp(command, "moat_concat") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("moat_concat: concatenate 2 MOAT samples (different inputs)\n");
            printf("syntax: moat_concat <file>\n");
            printf("Note: a PSUADE MOAT datafile should have been loaded before\n");
            printf("      using this command. \n");
         }
         printf("ERROR: no datafile loaded (use load command).\n");
      }
      else if (!strcmp(command, "moat_concat") && nInputs >= 1)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("moat_concat: concatenate 2 MOAT samples (different inputs)\n");
            printf("syntax: moat_concat <file> \n");
            printf("Note: a PSUADE MOAT datafile should have been loaded before\n");
            printf("      using this command. \n");
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
            nReps = nSamples / (nInputs + 1);
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
            names = inputNames;
            inputNames = new char*[nInputs+ind];
            for (ii = 0; ii < nInputs; ii++) inputNames[ii] = names[ii];
            pINames.clean();
            psuadeIO_->getParameter("input_names", pINames);
            names = pINames.strArray_;
            for (ii = nInputs; ii < nInputs+ind; ii++)
            {
               inputNames[ii] = new char[200];
               strcpy(inputNames[ii], names[ii-nInputs]);
            }
            tempX = sampleInputs;
            if (sampleOutputs != NULL) delete [] sampleOutputs;
            if (sampleStates  != NULL) delete [] sampleStates;
            kk = nReps * (nInputs + ind + 1) * (nInputs + ind);
            sampleInputs = new double[kk];
            kk = nReps * (nInputs + ind + 1) * nOutputs;
            sampleOutputs = new double[kk];
            for (ii = 0; ii < kk; ii++) sampleOutputs[ii] = PSUADE_UNDEFINED;
            kk = nReps * (nInputs + ind + 1);
            sampleStates = new int[kk];
            for (ii = 0; ii < kk; ii++) sampleStates[ii] = 0;
            for (ii = 0; ii < nReps; ii++)
            {
               for (jj = 0; jj <= nInputs; jj++)
               {
                  for (kk = 0; kk < nInputs; kk++)
                  {
                     ind2 = ii * (nInputs + ind + 1) * (nInputs + ind);
                     sampleInputs[ind2+jj*(nInputs+ind)+kk] = 
                            tempX[ii*(nInputs+1)*nInputs+jj*nInputs+kk];
                  }
               }
               for (jj = nInputs+1; jj < nInputs+ind+1; jj++)
               {
                  for (kk = 0; kk < nInputs; kk++)
                  {
                     ind2 = ii * (nInputs + ind + 1) * (nInputs + ind);
                     sampleInputs[ind2+jj*(nInputs+ind)+kk] = 
                           sampleInputs[ind2+nInputs*(nInputs+ind)+kk]; 
                  }
               }
            }
            if (tempX != NULL) delete [] tempX;
            ioPtr->getParameter("input_sample", pPtr);
            tempX = pPtr.dbleArray_;
            pPtr.dbleArray_ = NULL;
            for (ii = 0; ii < nReps; ii++)
            {
               for (jj = 0; jj <= nInputs; jj++)
               {
                  for (kk = nInputs; kk < nInputs+ind; kk++)
                  {
                     ind2 = ii * (nInputs + ind + 1) * (nInputs + ind);
                     sampleInputs[ind2+jj*(nInputs+ind)+kk] = 
                            tempX[ii*(ind+1)*ind+kk-nInputs];
                  }
               }
               for (jj = nInputs+1; jj < nInputs+ind+1; jj++)
               {
                  for (kk = nInputs; kk < nInputs+ind; kk++)
                  {
                     ind2 = ii * (nInputs + ind + 1) * (nInputs + ind);
                     sampleInputs[ind2+jj*(nInputs+ind)+kk] = 
                            tempX[ii*(ind+1)*ind+(jj-nInputs)*ind+kk-nInputs];
                  }
               }
            }
            delete [] tempX;
            pLower.clean();
            ioPtr->getParameter("input_lbounds", pLower);
            tempW = pLower.dbleArray_;
            tempV = iLowerB;
            iLowerB = new double[nInputs+ind];
            for (ii = 0; ii < nInputs; ii++) iLowerB[ii] = tempV[ii];
            for (ii = nInputs; ii < nInputs+ind; ii++)
               iLowerB[ii] = tempW[ii-nInputs];
            if (tempV != NULL) delete [] tempV;
            pLower.clean();
            pUpper.clean();
            ioPtr->getParameter("input_ubounds", pUpper);
            tempW = pUpper.dbleArray_;
            tempV = iUpperB;
            iUpperB = new double[nInputs+ind];
            for (ii = 0; ii < nInputs; ii++) iUpperB[ii] = tempV[ii];
            for (ii = nInputs; ii < nInputs+ind; ii++)
               iUpperB[ii] = tempW[ii-nInputs];
            if (tempV != NULL) delete [] tempV;
            pUpper.clean();
            delete ioPtr;
            nSamples = (nInputs + ind + 1) * nReps;;
            nInputs += ind;
            psuadeIO_->updateInputSection(nSamples,nInputs,NULL,iLowerB,
                                     iUpperB,sampleInputs,inputNames); 
            psuadeIO_->updateOutputSection(nSamples,nOutputs,
                                sampleOutputs, sampleStates, NULL); 
            psuadeIO_->updateMethodSection(-1, nSamples, 1, -1, -1);
            printf("The two samples have been concatenated.\n");
            printf("The new sample has nInputs = %d\n", nInputs);
            printf("                  nSamples = %d\n", nSamples);
            printf("                  nOutputs = %d\n", nOutputs);
            printf("Use 'write' to write the expanded sample to a file.\n");
         }
         else
         {
            printf("ERROR reading second sample file.\n");
         }
      }

      // split the file into 2 separate data file
      // ----------------------------------------------------------------
      else if (!strcmp(command, "splitdata"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("splitdata: split the existing data set into two subsets.\n");
            printf("Note: the two sets of data will be stored in psuadeSplit1\n");
            printf("      and psuadeSplit2.\n");
            continue;
         }
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no data to split.\n");
            printf("       Use load to put data into local memory first.\n");
         }
         else if (nSamples <= 0) printf("Reload data file.\n");
         else
         {
            printf("The current sample size is %d.\n", nSamples);
            sprintf(pString,
                    "Sample size of the first set? (1 - %d) ",nSamples-1);
            kk = getInt(1, nSamples, pString);
            sprintf(pString, "Random draw from original sample ? (y or n) ");
            getString(pString, winput);
            tempX  = new double[kk*nInputs];
            tempY  = new double[kk*nOutputs];
            states = new int[kk];
            tags   = new int[nSamples];
            for (ii = 0; ii < nSamples; ii++) tags[ii] = 0;
            for (ii = 0; ii < kk; ii++)
            {
               if (winput[0] == 'y')
               {
                  ind = PSUADE_rand() % nSamples;
                  ind2 = 0;
                  while (tags[ind2] == 1 && ind2 < 1000)
                  {
                     ind = PSUADE_rand() % nSamples;
                     ind2++;
                  }
                  if (tags[ind] == 1)
                     for (ind = 0; ind < nSamples; ind++)
                        if (tags[ind] == 0) break;
                  if (tags[ind] == 1)
                  {
                     printf("ERROR: cannot split data. \n");
                     continue;
                  }
               } else ind = ii;
               for (jj = 0; jj < nInputs; jj++)
                  tempX[ii*nInputs+jj] = sampleInputs[ind*nInputs+jj];
               for (jj = 0; jj < nOutputs; jj++)
                  tempY[ii*nOutputs+jj] = sampleOutputs[ind*nOutputs+jj];
               tags[ind] = 1;
               states[ii] = sampleStates[ind];
            }
            ind = 0;
            for (ii = 0; ii < nSamples; ii++)
            {
               if (tags[ii] == 0)
               {
                  for (jj = 0; jj < nInputs; jj++)
                     sampleInputs[ind*nInputs+jj]=sampleInputs[ii*nInputs+jj];
                  for (jj = 0; jj < nOutputs; jj++)
                     sampleOutputs[ind*nOutputs+jj] = 
                                   sampleOutputs[ii*nOutputs+jj];
                  sampleStates[ind] = sampleStates[ii]; 
                  ind++;
               }
            }
            psuadeIO_->updateInputSection(kk,nInputs,NULL,NULL,NULL,
                                          tempX,NULL); 
            psuadeIO_->updateOutputSection(kk,nOutputs,
                                           tempY,states,outputNames); 
            psuadeIO_->updateMethodSection(-1,kk,-1,-1,-1);
            psuadeIO_->writePsuadeFile("psuadeSplit1",0);

            psuadeIO_->updateInputSection(nSamples-kk,nInputs,NULL,NULL,
                                          NULL,sampleInputs,NULL); 
            psuadeIO_->updateOutputSection(nSamples-kk,nOutputs,
                             sampleOutputs,sampleStates,outputNames); 
            psuadeIO_->updateMethodSection(-1,nSamples-kk,-1,-1,-1);
            psuadeIO_->writePsuadeFile("psuadeSplit2",0);

            printf("The 2 data files are in psuadeSplit1 and psuadeSplit2.\n");
            printf("The loaded data in local memory have been destroyed.\n");
            fflush(stdout);
            delete [] tempX;
            delete [] tempY;
            delete [] states;
            delete [] sampleInputs;
            delete [] sampleOutputs;
            delete [] sampleStates;
            if (inputNames != NULL)
            {
               for (ii = 0; ii < nInputs; ii++) delete [] inputNames[ii];
               delete [] inputNames;
            }
            if (outputNames != NULL)
            {
               for (ii = 0; ii < nOutputs; ii++) delete [] outputNames[ii];
               delete [] outputNames;
            }
            inputNames  = NULL;
            outputNames = NULL;
            sampleInputs = NULL;
            sampleOutputs = NULL;
            sampleStates = NULL;
            nSamples = 0;
            nInputs = 0;
            nOutputs = 0;
         }
      }

      // ================================================================
      // UQ/SA commands
      // ================================================================
      // uncertainty analysis + output to a matlab file
      else if (!strcmp(command, "ua") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("ua: uncertainty analysis (compute moments).\n");
            printf("syntax: ua (after data have been loaded)\n");
            continue;
         }
         printf("ERROR: no data to analyze (load data first).\n");
      }
      else if (!strcmp(command, "ua") && nInputs > 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("ua: uncertainty analysis (compute moments).\n");
            printf("syntax: ua (after data have been loaded)\n");
            continue;
         }
         if (psuadeIO_ == NULL) printf("ERROR: data not loaded yet.\n");
         else
         {
            outputID = 0;
            if (nOutputs > 1)
            {
               sprintf(pString, "Enter output number (1 - %d) : ", 
                       nOutputs);
               outputID = getInt(0, nOutputs, pString);
               outputID--;
            }
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
      }

      // correlation analysis 
      // ----------------------------------------------------------------
      else if (!strcmp(command, "ca") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("ca: correlation analysis\n");
            printf("syntax: ca (after data have been loaded)\n");
            continue;
         }
         printf("ERROR: no data to analyze (load data first).\n");
      }
      else if (!strcmp(command, "ca") && nInputs > 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("ca: correlation analysis\n");
            printf("syntax: ca (after data have been loaded)\n");
            continue;
         }
         if (psuadeIO_ == NULL) printf("ERROR: data not loaded yet.\n");
         else
         {
            outputID = 0;
            if (nOutputs > 1)
            {
               sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
               outputID = getInt(1, nOutputs, pString);
               outputID--;
            }
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
      }

      // analysis of variation
      // ----------------------------------------------------------------
      else if (!strcmp(command, "anova") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("anova: analysis of variation\n");
            printf("syntax: anova (after data have been loaded)\n");
            continue;
         }
         printf("ERROR: no data to analyze (load data first).\n");
      }
      else if (!strcmp(command, "anova") && nInputs > 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("anova: analysis of variation\n");
            printf("syntax: anova (after data have been loaded)\n");
            continue;
         }
         if (psuadeIO_ == NULL) printf("ERROR: data not loaded yet.\n");
         else
         {
            outputID = 0;
            if (nOutputs > 1)
            {
               sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
               outputID = getInt(1, nOutputs, pString);
               outputID--;
            }
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
      }

      // Morris' std deviation versus revised mean plot
      else if (!strcmp(command, "moat") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("moat: Morris screening analysis\n");
            printf("syntax: moat (after data have been loaded)\n");
            continue;
         }
         printf("ERROR: no data to analyze (load data first).\n");
      }
      else if (!strcmp(command, "moat") && nInputs > 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("moat: Morris screening analysis\n");
            printf("syntax: moat (after data have been loaded)\n");
            continue;
         }
         if (psuadeIO_ == NULL) printf("ERROR: data not loaded yet.\n");
         else
         {
            outputID = 0;
            if (nOutputs > 1)
            {
               sprintf(pString, "Enter output number (1 - %d) = ", nOutputs);
               outputID = getInt(1, nOutputs, pString);
               outputID--;
            }
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
      }

      // fractional factorial first and second order effects
      // ----------------------------------------------------------------
      else if (!strcmp(command, "ff") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("ff: Fractional factorial screening analysis\n");
            printf("syntax: ff (after data have been loaded)\n");
            continue;
         }
         printf("ERROR: no data to analyze (load data first).\n");
      }
      else if (!strcmp(command, "ff") && nInputs > 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("ff: Fractional factorial screening analysis\n");
            printf("syntax: ff (after data have been loaded)\n");
            continue;
         }
         if (psuadeIO_ == NULL) printf("ERROR: data not loaded yet.\n");
         else
         {
            outputID = 0;
            if (nOutputs > 1)
            {
               sprintf(pString, "Enter output number (1 - %d) = ", nOutputs);
               outputID = getInt(1, nOutputs, pString);
               outputID--;
            }
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
      }

      // local sensitivity analysis fof first order effects
      // ----------------------------------------------------------------
      else if (!strcmp(command, "lsa") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("lsa: local sensitivity analysis\n");
            printf("syntax: lsa (after data have been loaded)\n");
            continue;
         }
         printf("ERROR: no data to analyze (load data first).\n");
      }
      else if (!strcmp(command, "lsa") && nInputs > 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("lsa: local sensitivity analysis\n");
            printf("syntax: lsa (after data have been loaded)\n");
            continue;
         }
         if (psuadeIO_ == NULL) printf("ERROR: data not loaded yet.\n");
         else
         {
            outputID = 0;
            if (nOutputs > 1)
            {
               sprintf(pString, "Enter output number (1 - %d) = ", nOutputs);
               outputID = getInt(1, nOutputs, pString);
               outputID--;
            }
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
      }

      // MARS screening of parameters
      // ----------------------------------------------------------------
      else if (!strcmp(command, "mars_sa") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("mars_sa: MARS-based sensitivity analysis\n");
            printf("syntax: mars_sa (after data have been loaded)\n");
            continue;
         }
         printf("ERROR: no data to analyze (load data first).\n");
      }
      else if (!strcmp(command, "mars_sa") && nInputs > 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("mars_sa: MARS-based sensitivity analysis\n");
            printf("syntax: mars_sa (after data have been loaded)\n");
            continue;
         }
         if (psuadeIO_ == NULL) printf("ERROR: data not loaded yet.\n");
         else
         {
            outputID = 0;
            if (nOutputs > 1)
            {
               sprintf(pString, "Enter output number (1 - %d) = ", nOutputs);
               outputID = getInt(1, nOutputs, pString);
               outputID--;
            }
            sprintf(pString,"MARS (0) or MARS with bagging (1) ? ");
            kk = getInt(0, 1, pString);
            if (kk == 0) faType = PSUADE_RS_MARS;
            else         faType = PSUADE_RS_MARSB;
            faPtr = genFA(faType, nInputs, iOne, nSamples);
            if (faPtr != NULL)
            {
               faPtr->setBounds(iLowerB, iUpperB);
               faPtr->setOutputLevel(outputLevel_);
               if (faType == PSUADE_RS_MARSB) faPtr->setOutputLevel(4);
               tempY = new double[nSamples];
               for (ii = 0; ii < nSamples; ii++)
                  tempY[ii] = sampleOutputs[ii*nOutputs+outputID];
               kk = -999;
               status = faPtr->genNDGridData(sampleInputs,tempY,&kk,NULL,NULL);
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
      }

      // GP screening of parameters
      // ----------------------------------------------------------------
      else if (!strcmp(command, "gp_sa") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("gp_sa: Gaussian Process-based sensitivity analysis\n");
            printf("syntax: gp_sa (after data have been loaded)\n");
            continue;
         }
         printf("ERROR: no data to analyze (load data first).\n");
      }
      else if (!strcmp(command, "gp_sa") && nInputs > 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("gp_sa: Gaussian Process-based sensitivity analysis\n");
            printf("syntax: gp_sa (after data have been loaded)\n");
            continue;
         }
         if (psuadeIO_ == NULL) printf("ERROR: data not loaded yet.\n");
         else
         {
            outputID = 0;
            if (nOutputs > 1)
            {
               sprintf(pString, "Enter output number (1 - %d) = ", nOutputs);
               outputID = getInt(1, nOutputs, pString);
               outputID--;
            }
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
            faPtr = genFA(faType, nInputs, iOne, nSamples);
            if (faPtr != NULL)
            {
               faPtr->setBounds(iLowerB, iUpperB);
               faPtr->setOutputLevel(outputLevel_);
               if (faType == PSUADE_RS_KR)
               {
                  strcpy(pString, "setMode3");
                  targv[0] = (char *) pString;
                  faPtr->setParams(1, targv);
               }
               int rsKeep = psRSExpertMode_;
               psRSExpertMode_ = 0;
               tempY = new double[nSamples];
               for (ii = 0; ii < nSamples; ii++)
                  tempY[ii] = sampleOutputs[ii*nOutputs+outputID];
               kk = -999;
               status = faPtr->genNDGridData(sampleInputs,tempY,&kk,NULL,NULL);
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
      }

      // sum-of-trees screening of parameters
      // ----------------------------------------------------------------
      else if (!strcmp(command, "sot_sa") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("sot_sa: Sum-of-trees-based sensitivity analysis\n");
            printf("syntax: sot_sa (after data have been loaded)\n");
            continue;
         }
         printf("ERROR: no data to analyze (load data first).\n");
      }
      else if (!strcmp(command, "sot_sa") && nInputs > 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("sot_sa: Sum-of-trees-based sensitivity analysis\n");
            printf("syntax: sot_sa (after data have been loaded)\n");
            continue;
         }
         if (psuadeIO_ == NULL) printf("ERROR: data not loaded yet.\n");
         else
         {
            outputID = 0;
            if (nOutputs > 1)
            {
               sprintf(pString, "Enter output number (1 - %d) = ", nOutputs);
               outputID = getInt(1, nOutputs, pString);
               outputID--;
            }
            faType = PSUADE_RS_SOTS;
            faPtr  = genFA(faType, nInputs, iOne, nSamples);
            if (faPtr != NULL)
            {
               faPtr->setBounds(iLowerB, iUpperB);
               faPtr->setOutputLevel(outputLevel_);
               tempY = new double[nSamples];
               for (ii = 0; ii < nSamples; ii++)
                  tempY[ii] = sampleOutputs[ii*nOutputs+outputID];
               kk = -999;
               status = faPtr->genNDGridData(sampleInputs,tempY,&kk,NULL,NULL);
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
         }
      }

      // quantitative main effect analysis + output to a matlab file
      // ----------------------------------------------------------------
      else if (!strcmp(command, "me") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("me: main effect analysis (variance-based)\n");
            printf("syntax: me (after data have been loaded)\n");
            continue;
         }
         printf("ERROR: no data to analyze (load data first).\n");
      }
      else if (!strcmp(command, "me") && nInputs > 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("me: main effect analysis (variance-based)\n");
            printf("syntax: me (after data have been loaded)\n");
            continue;
         }
         if (psuadeIO_ == NULL) printf("ERROR: data not loaded yet.\n");
         else
         {
            outputID = 0;
            if (nOutputs > 1)
            {
               sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
               outputID = getInt(1, nOutputs, pString);
               outputID--;
            }
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
      }

      // quantitative interaction effect analysis 
      // ----------------------------------------------------------------
      else if (!strcmp(command, "ie") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("ie: 2-way interaction effect analysis (variance-based)\n");
            printf("syntax: ie (after data have been loaded)\n");
            continue;
         }
         printf("ERROR: no data to analyze.\n");
      }
      else if (!strcmp(command, "ie") && nInputs > 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("ie: 2-way interaction effect analysis (variance-based)\n");
            printf("syntax: ie (after data have been loaded)\n");
            continue;
         }
         if (psuadeIO_ == NULL) printf("ERROR: data not loaded yet.\n");
         else
         {
            outputID = 0;
            if (nOutputs > 1)
            {
               sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
               outputID = getInt(1, nOutputs, pString);
               outputID--;
            }
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
            if (pdata->nDbles_ >= nInputs)
            {
               printEquals(0);
               printf("Two-way Interaction Effect Statistics: \n");
               if (pdata->dbleData_ > 0)
               {
                  for (ii = 0; ii < nInputs; ii++)
                     for (jj = ii+1; jj < nInputs; jj++)
                        printf("Inputs %4d %4d: Sobol' interaction effect = %12.4e\n",
                               ii+1,jj+1,pdata->dbleArray_[ii*nInputs+jj]);
                  if (psPlotTool_ == 1) fp = fopen("scilabaie.sci", "w");
                  else                  fp = fopen("matlabaie.m", "w");
                  if (fp != NULL)
                  {
                     if (psPlotTool_ == 1)
                     {
                        fprintf(fp,"// This file contains Sobol' second order indices\n");
                        fprintf(fp,"// set sortFlag = 1 and set nn to be the number\n");
                        fprintf(fp,"// of inputs to display.\n");
                     }
                     else
                     {
                        fprintf(fp,"%% This file contains Sobol' second order indices\n");
                        fprintf(fp,"%% set sortFlag = 1 and set nn to be the number\n");
                        fprintf(fp,"%% of inputs to display.\n");
                     }
                     fprintf(fp, "sortFlag = 0;\n");
                     fprintf(fp, "nn = %d;\n", nInputs);
                     fprintf(fp, "Mids = [\n");
                     for (ii = 0; ii < nInputs*nInputs; ii++) 
                        fprintf(fp,"%24.16e\n", pdata->dbleArray_[ii]);
                     fprintf(fp, "];\n");
                     if (inputNames == NULL)
                     {
                        fprintf(fp, "  Str = {");
                        for (ii = 0; ii < nInputs-1; ii++) fprintf(fp,"'X%d',",ii+1);
                        fprintf(fp,"'X%d'};\n",nInputs);
                     }
                     else
                     {
                        fprintf(fp, "  Str = {");
                        for (ii = 0; ii < nInputs-1; ii++)
                        {
                           if (inputNames[ii] != NULL) 
                                fprintf(fp,"'%s',",inputNames[ii]);
                           else fprintf(fp,"'X%d',",ii+1);
                        }
                        if (inputNames[nInputs-1] != NULL)
                             fprintf(fp,"'%s'};\n",inputNames[nInputs-1]);
                        else fprintf(fp,"'X%d'};\n",nInputs);
                     }
                     fwritePlotCLF(fp);
                     fprintf(fp, "ymin = min(Mids);\n");
                     fprintf(fp, "ymax = max(Mids);\n");
                     fprintf(fp, "h2 = 0.05 * (ymax - ymin);\n");
                     if (psPlotTool_ == 1)
                          fprintf(fp, "Mids = matrix(Mids, %d, %d);\n",nInputs,nInputs);
                     else fprintf(fp, "Mids = reshape(Mids, %d, %d);\n",nInputs,nInputs);
                     fprintf(fp, "Mids = Mids';\n");
                     if (psPlotTool_ == 1)
                     {
                        fprintf(fp, "drawlater\n");
                        fprintf(fp, "hist3d(Mids);\n");
                        fprintf(fp, "a=gca();\n");
                        fprintf(fp, "a.data_bounds=[0, 0, 0; %d+1, %d+1, ymax];\n",nInputs,
                                nInputs);
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
                        fprintf(fp, "axis([0 %d+1 0 %d+1 0 ymax])\n", nInputs, nInputs);
                        fprintf(fp, "set(gca,'XTickLabel',Str);\n");
                        fprintf(fp, "set(gca,'YTickLabel',Str);\n");
                        fprintf(fp, "set(gca, 'fontsize', 12)\n");
                        fprintf(fp, "set(gca, 'fontweight', 'bold')\n");
                     }
                     fwritePlotAxes(fp);
                     fwritePlotTitle(fp,"Sobol Second Order Indices (+ first order)");
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
                  printf("Total variance = 0. Hence, no interaction effect plot.\n");
               }
               pdata->clean();
            }
         }
      }

      // quantitative total sensitivity analysis 
      // ----------------------------------------------------------------
      else if (!strcmp(command, "tsi") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("tsi: total sensitivity analysis (variance-based)\n");
            printf("     (suitable for raw data)\n");
            printf("syntax: tsi (after data have been loaded)\n");
            continue;
         }
         printf("ERROR: no data to analyze.\n");
      }
      else if (!strcmp(command, "tsi") && nInputs > 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("tsi: total sensitivity analysis (variance-based)\n");
            printf("     (suitable for raw data)\n");
            printf("syntax: tsi (after data have been loaded)\n");
            continue;
         }
         if (psuadeIO_ == NULL) printf("ERROR: data not loaded yet.\n");
         else
         {
            if (nInputs > 10 || nSamples < 50*nInputs) 
            {
               printf("This command is not recommended for small sample or.\n");
               printf("large number of inputs.\n");
               printf("Need at most 10 inputs.\n");
               printf("Need at least %d sample points.\n",50*nInputs);
               continue;
            }
            outputID = 0;
            if (nOutputs > 1)
            {
               sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
               outputID = getInt(1, nOutputs, pString);
               outputID--;
            }
            TSIAnalyzer *tsiAnalyzer = new TSIAnalyzer();
            aPtr.printLevel_ = outputLevel_;
            aPtr.nSamples_ = nSamples;
            aPtr.nInputs_ = nInputs;
            aPtr.nOutputs_ = nOutputs;
            aPtr.sampleInputs_ = sampleInputs;
            aPtr.sampleOutputs_ = sampleOutputs;
            aPtr.iLowerB_ = iLowerB;
            aPtr.iUpperB_ = iUpperB;
            aPtr.outputID_ = outputID;
            aPtr.sampleStates_ = sampleStates;
            aPtr.ioPtr_ = psuadeIO_;
            tsiAnalyzer->analyze(aPtr);
            delete tsiAnalyzer;
         }
      }

      // scatterplot of output against any one input using Pgplot
      // ----------------------------------------------------------------
      else if (!strcmp(command, "meplot") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("meplot: main effect plot using Pgplot (obsolete)\n");
            continue;
         }
         printf("ERROR: no data to analyze.\n");
      }
      else if (!strcmp(command, "meplot") && nInputs > 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("meplot: main effect plot using Pgplot (obsolete)\n");
            continue;
         }
         iInd = 0;
         if (nInputs > 1)
         {
            sprintf(pString, "Enter input number (1 - %d) : ", nInputs);
            iInd = getInt(1, nInputs, pString);
            iInd--;
         }
         outputID = 0;
         if (nOutputs > 1)
         {
            sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
            outputID = getInt(1, nOutputs, pString);
            outputID--;
         }
         Xmin = iLowerB[iInd];
         Xmax = iUpperB[iInd];
         Ymin =   1.0e35;
         Ymax = - 1.0e35;
         for (sInd = 0; sInd < nSamples; sInd++)
         {
            if (sampleOutputs[sInd*nOutputs+outputID] > Ymax)
               Ymax = sampleOutputs[sInd*nOutputs+outputID];
            if (sampleOutputs[sInd*nOutputs+outputID] < Ymin)
               Ymin = sampleOutputs[sInd*nOutputs+outputID];
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
         nSamplesLong = (long) nSamples;
         tempX = new double[nSamples];
         tempY = new double[nSamples];
         for (sInd = 0; sInd < nSamples; sInd++)
         {
            tempX[sInd] = sampleInputs[sInd*nInputs+iInd];
            tempY[sInd] = sampleOutputs[sInd*nOutputs+outputID];
         }
         sortDbleList2(nSamples, tempX, tempY);
         PlotScatter2D(nSamplesLong, tempX, tempY);
         printf("Enter any character and return to continue\n");
         scanf("%s", command);
         Plotend();
         fflush(stdout);
         printf("\n");
         fflush(stdout);
         fgets(lineIn,500,stdin); 
         delete [] tempX;
         delete [] tempY;
      }

      // ----------------------------------------------------------------
      else if (!strcmp(command, "meplot2") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("meplot2: two-way interaction plot using Pgplot (obsolete)\n");
            continue;
         }
         printf("ERROR: no data to analyze.\n");
      }
      else if (!strcmp(command, "meplot2") && nInputs < 2)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("meplot2: two-way interaction plot using Pgplot (obsolete)\n");
            continue;
         }
         printf("ERROR: meplot2 requires 2 or more inputs.\n");
      }
      else if (!strcmp(command, "meplot2") && nInputs >= 2)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("meplot2: two-way interaction plot using Pgplot (obsolete)\n");
            continue;
         }
         iInd1 = iInd2 = -1;
         if (nInputs == 2)
         {
            iInd1 = 0;
            iInd2 = 1;
         }
         else
         {
            sprintf(pString, "Enter first input number (1 - %d) : ",nInputs);
            iInd1 = getInt(1, nInputs, pString);
            iInd1--;
            sprintf(pString, "Enter second input number (1 - %d) : ",nInputs);
            iInd2 = getInt(1, nInputs, pString);
            iInd2--;
         }

         outputID = 0;
         if (nOutputs > 1)
         {
            sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
            outputID = getInt(1, nOutputs, pString);
            outputID--;
         }
         Xmin = iLowerB[iInd1];
         Xmax = iUpperB[iInd1];
         Ymin =   1.0e35;
         Ymax = - 1.0e35;
         for (sInd = 0; sInd < nSamples; sInd++)
         {
            if (sampleOutputs[sInd*nOutputs+outputID] > Ymax)
               Ymax = sampleOutputs[sInd*nOutputs+outputID];
            if (sampleOutputs[sInd*nOutputs+outputID] < Ymin)
               Ymin = sampleOutputs[sInd*nOutputs+outputID];
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
         nSamplesLong = (long) nSamples;
         tempInds = new double[nSamples];
         tempX = new double[nSamples];
         tempW = new double[nSamples];
         tempY = new double[nSamples];
         tempT = new double[nSamples];
         for (sInd = 0; sInd < nSamples; sInd++)
            tempInds[sInd] = (double) sInd;
         for (sInd = 0; sInd < nSamples; sInd++)
         {
            tempX[sInd] = sampleInputs[sInd*nInputs+iInd1];
            tempW[sInd] = sampleInputs[sInd*nInputs+iInd2];
            tempY[sInd] = sampleOutputs[sInd*nOutputs+outputID];
         }
         sortDbleList2(nSamples, tempX, tempInds);
         for (sInd = 0; sInd < nSamples; sInd++)
         {
            ind2 = (int) tempInds[sInd];
            tempT[sInd] = tempW[ind2];
         }
         for (sInd = 0; sInd < nSamples; sInd++)
            tempW[sInd] = tempT[sInd];
         for (sInd = 0; sInd < nSamples; sInd++)
         {
            ind2 = (int) tempInds[sInd];
            tempT[sInd] = tempY[ind2];
         }
         for (sInd = 0; sInd < nSamples; sInd++)
            tempY[sInd] = tempT[sInd];
         ind2 = 0;
         for (sInd = 1; sInd < nSamples; sInd++)
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
         fgets(lineIn,500,stdin); 
         delete [] tempX;
         delete [] tempW;
         delete [] tempY;
         delete [] tempT;
         delete [] tempInds;
      }

      // ----------------------------------------------------------------
      else if (!strcmp(command, "rsplot") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rsplot: two-input RS plot using Pgplot (obsolete)\n");
            printf("syntax: splot (no argument needed)\n");
            continue;
         }
         printf("ERROR: no data to analyze.\n");
      }
      else if (!strcmp(command, "rsplot") && nInputs < 2)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rsplot: two-input RS plot using Pgplot (obsolete)\n");
            printf("syntax: splot (no argument needed)\n");
            continue;
         }
         printf("ERROR: rsplot requires 2 or more inputs.\n");
      }
      else if (!strcmp(command, "rsplot") && nInputs >= 2)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rsplot: two-input RS plot using Pgplot (obsolete)\n");
            printf("syntax: splot (no argument needed)\n");
            continue;
         }
         pgPlotResponseSurface();
      }

      // generate matlab scatter plot file
      // ----------------------------------------------------------------
      else if (!strcmp(command, "splot") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("splot: create scatter plots (against each parameter).\n");
            printf("syntax: splot (no argument needed)\n");
            continue;
         }
         printf("ERROR: no data to analyze (load data first).\n");
      }
      else if (!strcmp(command, "splot") && nInputs > 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("splot: create scatter plots (against each parameter).\n");
            printf("syntax: splot (no argument needed)\n");
            continue;
         }
         if (psuadeIO_ == NULL) printf("ERROR: data not loaded yet.\n");
         else
         {
            outputID = 0;
            if (nOutputs > 1)
            {
               sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
               outputID = getInt(1, nOutputs, pString);
               outputID--;
            }

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
            if (psPlotTool_ == 1)
            {
               fprintf(fp, "// plotMode=0  : plot all in a single plot\n");
               fprintf(fp, "// plotMode!=0 : plot one at a time\n");
            }
            else
            {
               fprintf(fp, "%% plotMode=0  : plot all in a single plot\n");
               fprintf(fp, "%% plotMode!=0 : plot one at a time\n");
            }
            fprintf(fp, "plotMode=0;\n");
            fprintf(fp, "Y = [\n");
            for (sInd = 0; sInd < nSamples; sInd++)
               fprintf(fp, "%24.16e\n",sampleOutputs[sInd*nOutputs+outputID]);
            fprintf(fp, "];\n");
            for (iInd = 0; iInd < nInputs; iInd++)
            {
               fprintf(fp, "X%d = [\n", iInd+1);
               for (sInd = 0; sInd < nSamples; sInd++)
                  fprintf(fp, "%24.16e\n",sampleInputs[sInd*nInputs+iInd]);
               fprintf(fp, "];\n");
            }
            fprintf(fp, "S = [\n");
            for (sInd = 0; sInd < nSamples; sInd++)
               fprintf(fp, "%d\n",sampleStates[sInd]);
            fprintf(fp, "];\n");
            if (psPlotTool_ == 1)
            {
               for (iInd = 0; iInd < nInputs; iInd++)
               {
                  fwritePlotCLF(fp);
                  fprintf(fp, "drawlater\n");
                  fprintf(fp, "plot(X%d,Y,'.','markersize',10)\n",iInd+1);
                  fprintf(fp, "a = gca();\n");
                  fprintf(fp, "a.children.children.mark_foreground = 5;\n");
                  sprintf(winput, "%s vs %s", outputNames[outputID],
                          inputNames[iInd]);
                  fwritePlotTitle(fp, winput);
                  fwritePlotAxes(fp);
                  fwritePlotXLabel(fp, inputNames[iInd]);
                  fwritePlotYLabel(fp, outputNames[outputID]);
                  fprintf(fp, "drawnow\n");
                  if (iInd < nInputs-1) 
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
               kk = ((int) pow(1.0*nInputs-0.1, 0.5)) + 1;
               ll = kk;
               while ((ll - 1) * kk >= nInputs) ll--; 
               for (iInd = 0; iInd < nInputs; iInd++)
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
                  fprintf(fp,"plot(X%d(iset),Y(iset),'rX','markersize',fs)\n",iInd+1);
                  fprintf(fp,"hold on\n");
                  fprintf(fp,"iset = find(S == 1);\n");
                  fprintf(fp,"plot(X%d(iset),Y(iset),'b*','markersize',fs)\n",iInd+1);
                  fprintf(fp,"hold off\n");
                  fprintf(fp,"axis([%24.16e %24.16e min(Y) max(Y)])\n",iLowerB[iInd],
                          iUpperB[iInd]);
                  sprintf(winput, "%s vs %s", outputNames[outputID],
                          inputNames[iInd]);
                  fwritePlotTitle(fp, winput);
                  fwritePlotAxes(fp);
                  fwritePlotXLabel(fp, inputNames[iInd]);
                  fwritePlotYLabel(fp, outputNames[outputID]);
               }
               printf("matlabsp.m is now available for scatter plots.\n");
            }
            fclose(fp);    
         }
      }

      // generate response surface of any one inputs and write the
      // grid data to file for display with matlab
      // ----------------------------------------------------------------
      else if (!strcmp(command, "rs1") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rs1: response surface plot in one parameter\n");
            printf("syntax: rs1 (no argument needed)\n");
            continue;
         }
         printf("ERROR: no data to analyze (load data first).\n");
      }
      else if (!strcmp(command, "rs1") && nInputs >= 1)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rs1: response surface plot in one parameter\n");
            printf("syntax: rs1 (no argument needed)\n");
            continue;
         }
         nPtsPerDim = 256;
         faFlag = 1;
         faPtr = genFAInteractive(psuadeIO_, faFlag);
         if (faPtr == NULL) {printf("ERROR detected.\n"); continue;}
         faPtr->setNPtsPerDim(nPtsPerDim);
         faPtr->setBounds(iLowerB, iUpperB);
         faPtr->setOutputLevel(outputLevel_);
         inputSettings = new double[nInputs];
         if (nInputs > 1)
         {
            sprintf(pString, "Enter the input for x axis (1 - %d) : ", nInputs);
            iplot1 = getInt(1, nInputs, pString);
            iplot1--;
         }
         else iplot1 = 0;
         if (nInputs > 1)
         {
            sprintf(pString,
                    "Set other nominal values automatically ? (y or n) ");
            getString(pString, command);
         }
         if (command[0] == 'y')
         {
            for (iInd1 = 0; iInd1 < nInputs; iInd1++)
            {
               if (iInd1 != iplot1)
                    inputSettings[iInd1] = 0.5*(iLowerB[iInd1]+iUpperB[iInd1]);
               else inputSettings[iInd1] = 1.0;
            }
         }
         else
         {
            for (iInd1 = 0; iInd1 < nInputs; iInd1++)
            {
               if (iInd1 != iplot1)
               {
                  sprintf(pString,
                          "Enter nominal value for input %d (%e - %e): ", 
                          iInd1+1, iLowerB[iInd1], iUpperB[iInd1]);
                  inputSettings[iInd1] = getDouble(pString);
               }
               else inputSettings[iInd1] = 1.0;
            }
         }
         jplot = 0;
         if (nOutputs > 1)
         {
            sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
            jplot = getInt(1, nOutputs, pString);
            jplot--;
         }

         faYIn = new double[nSamples];
         for (sInd = 0; sInd < nSamples; sInd++)
            faYIn[sInd] = sampleOutputs[sInd*nOutputs+jplot];

         faPtr->gen1DGridData(sampleInputs,faYIn,iplot1,
                       inputSettings, &faLeng, &faXOut,&faYOut);

         if (psPlotTool_ == 1)
         {
            fp = fopen("scilabrs1.sci", "w");
            if (fp == NULL)
            {
               printf("ERROR: cannot open file scilabrs1.sci.\n");
               continue;
            }
            fwritePlotCLF(fp);
            fprintf(fp, "A = [\n");
            for (sInd = 0; sInd < faLeng; sInd++)
               fprintf(fp, "%e\n", faYOut[sInd]);
            fprintf(fp, "];\n");
            fprintf(fp, "X = [\n");
            for (sInd = 0; sInd < faLeng; sInd++)
               fprintf(fp, "%e\n", faXOut[sInd]);
            fprintf(fp, "];\n");
            fprintf(fp, "plot(X,A,'-')\n");
            fprintf(fp, "a = gca();\n");
            fprintf(fp, "a.children.children.thickness = 4;\n");
            fwritePlotAxes(fp);
            fwritePlotXLabel(fp, inputNames[iplot1]);
            fwritePlotYLabel(fp, outputNames[jplot]);
            sprintf(winput, "Plot for %s", outputNames[jplot]);
            fwritePlotTitle(fp, winput);
            fclose(fp);
            printf("scilabrs1.sci is now available.\n");
         }
         else
         {
            fp = fopen("matlabrs1.m", "w");
            if (fp == NULL)
            {
               printf("ERROR: cannot open file matlabrs1.m.\n");
               continue;
            }
            fwritePlotCLF(fp);
            fprintf(fp, "clf\n");
            fprintf(fp, "A = [\n");
            for (sInd = 0; sInd < faLeng; sInd++)
               fprintf(fp, "%e\n", faYOut[sInd]);
            fprintf(fp, "];\n");
            fprintf(fp, "X = [\n");
            for (sInd = 0; sInd < faLeng; sInd++)
               fprintf(fp, "%e\n", faXOut[sInd]);
            fprintf(fp, "];\n");
            fprintf(fp, "plot(X,A,'-','lineWidth',4)\n");
            fprintf(fp, "hold on\n");
            Ymin = faYOut[0];
            for (sInd = 1; sInd < faLeng; sInd++)
               if (faYOut[sInd] < Ymin) Ymin = faYOut[sInd];
            Ymax = faYOut[0];
            for (sInd = 1; sInd < faLeng; sInd++)
               if (faYOut[sInd] > Ymax) Ymax = faYOut[sInd];
            printf("Ymin and Ymax found = %e %e.\n", Ymin, Ymax);
            sprintf(pString,"Set lower threshold ? (y or n) : ");
            getString(pString, command);
            fprintf(fp, "yminFlag = 0;\n");
            if (command[0] == 'y')
            {
               sprintf(pString,"Enter the lower threshold (min = %e) : ", 
                       Ymin);
               thresh = getDouble(pString);
               fprintf(fp, "ymin = %e;\n", thresh);
               fprintf(fp, "plot(X,ones(%d,1)*ymin,'r-')\n",faLeng);
            }
            sprintf(pString,"Set upper threshold ? (y or n) : ");
            getString(pString, command);
            if (command[0] == 'y')
            {
               sprintf(pString,"Enter the upper threshold (max = %e) : ", 
                       Ymax);
               thresh = getDouble(pString);
               fprintf(fp, "ymax = %e;\n", thresh);
               fprintf(fp, "plot(X,ones(%d,1)*ymax,'r-')\n",faLeng);
            }
            fwritePlotAxes(fp);
            fwritePlotXLabel(fp, inputNames[iplot1]);
            fwritePlotYLabel(fp, outputNames[jplot]);
            sprintf(winput, "Plot for %s", outputNames[jplot]);
            fwritePlotTitle(fp, winput);
            fclose(fp);
            printf("matlabrs1.m is now available.\n");
         }
         delete [] faXOut;
         delete [] faYIn;
         delete [] faYOut;
         delete [] inputSettings;
         delete faPtr;
         faPtr = NULL;
      }

      // generate response surface of any two inputs and write the
      // grid data to file for display with matlab
      // ----------------------------------------------------------------
      else if (!strcmp(command, "rs2") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rs2: response surface plot in two parameters\n");
            printf("syntax: rs2 (no argument needed)\n");
            continue;
         }
         printf("ERROR: no data to analyze (load data first).\n");
      }
      else if (!strcmp(command, "rs2") && nInputs < 2)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rs2: response surface plot in two parameters\n");
            printf("syntax: rs2 (no argument needed)\n");
            continue;
         }
         printf("ERROR: rs requires 2 or more inputs.\n");
      }
      else if (!strcmp(command, "rs2") && nInputs >= 2)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rs2: response surface plot in two parameters\n");
            printf("syntax: rs2 (no argument needed)\n");
            continue;
         }
         nPtsPerDim = 64;
         sprintf(pString, "Grid resolution ? (32 - 256) ");
         nPtsPerDim = getInt(32, 256, pString);
         faFlag = 1;
         faPtr = genFAInteractive(psuadeIO_, faFlag);
         if (faPtr == NULL) {printf("ERROR detected.\n"); continue;}
         faPtr->setNPtsPerDim(nPtsPerDim);
         faPtr->setBounds(iLowerB, iUpperB);
         faPtr->setOutputLevel(outputLevel_);
         inputSettings = new double[nInputs];
         iplot1 = iplot2 = -1;
         sprintf(pString, "Enter the input for x axis (1 - %d) : ", nInputs);
         iplot1 = getInt(1, nInputs, pString);
         iplot1--;
         if (nInputs == 2)
         {
            if (iplot1 == 0) iplot2 = 1;
            else             iplot2 = 0;
         }
         while (iplot2 < 0 || iplot2 >= nInputs || iplot1 == iplot2)
         {
            sprintf(pString, "Enter the input for y axis (1 - %d), not %d : ",
                    nInputs, iplot1+1);
            iplot2 = getInt(1, nInputs, pString);
            iplot2--;
            if (iplot2 == iplot1)
            {
               printf("ERROR: cannot have same index for x and y axes.\n");
               iplot2 = -1;
            }
         }
         if (nInputs > 2)
         {
            sprintf(pString,
                    "Set other nominal values automatically ? (y or n) ");
            getString(pString, command);
         }
         if (command[0] == 'y')
         {
            for (iInd1 = 0; iInd1 < nInputs; iInd1++)
            {
               if (iInd1 != iplot1 && iInd1 != iplot2)
                    inputSettings[iInd1] = 0.5*(iLowerB[iInd1]+iUpperB[iInd1]);
               else inputSettings[iInd1] = 1.0;
            }
         }
         else
         {
            for (iInd1 = 0; iInd1 < nInputs; iInd1++)
            {
               if (iInd1 != iplot1 && iInd1 != iplot2)
               {
                  sprintf(pString,
                          "Enter nominal value for input %d (%e - %e): ", 
                          iInd1+1, iLowerB[iInd1], iUpperB[iInd1]);
                  inputSettings[iInd1] = getDouble(pString);
               }
               else inputSettings[iInd1] = 1.0;
            }
         }
         jplot = 0;
         if (nOutputs > 1)
         {
            sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
            jplot = getInt(1, nOutputs, pString);
            jplot--;
         }

         faYIn = new double[nSamples];
         for (sInd = 0; sInd < nSamples; sInd++)
            faYIn[sInd] = sampleOutputs[sInd*nOutputs+jplot];

         faPtr->gen2DGridData(sampleInputs,faYIn,iplot1,iplot2,
                       inputSettings, &faLeng, &faXOut,&faYOut);

         if (psPlotTool_ == 1)
         {
            fp = fopen("scilabrs2.sci", "w");
            if (fp == NULL)
            {
               printf("ERROR: cannot open file scilabrs2.sci.\n");
               continue;
            }
            fprintf(fp, "twoPlots = 1;\n");
            fprintf(fp, "A = [\n");
            for (sInd = 0; sInd < faLeng; sInd++)
               fprintf(fp, "%e\n", faYOut[sInd]);
            fprintf(fp, "];\n");
            fprintf(fp, "A = matrix(A,%d,%d);\n", nPtsPerDim, nPtsPerDim);
            fprintf(fp, "x = [\n");
            for (sInd = 0; sInd < faLeng; sInd+=nPtsPerDim)
               fprintf(fp, "%e\n", faXOut[sInd*2]);
            fprintf(fp, "];\n");
            fprintf(fp, "y = [\n");
            for (sInd = 0; sInd < nPtsPerDim; sInd++)
               fprintf(fp, "%e\n", faXOut[sInd*2+1]);
            fprintf(fp, "];\n");
            fwritePlotCLF(fp);
            fprintf(fp, "if twoPlots == 1\n");
            fprintf(fp, "drawlater\n");
            fprintf(fp, "subplot(1,2,1)\n");
            fprintf(fp, "mesh(x,y,A)\n");
            fprintf(fp, "h = get(\"hdl\");\n");
            fprintf(fp, "h.color_flag=1;\n");
            fprintf(fp, "h.color_mode=-2;\n");
            fprintf(fp, "bmin = min(min(A)); bmax = max(max(A));\n");
            fprintf(fp, "xset(\"colormap\",jetcolormap(64));\n");
            fprintf(fp, "colorbar(bmin,bmax);\n");
            fwritePlotAxes(fp);
            fwritePlotXLabel(fp, inputNames[iplot1]);
            fwritePlotYLabel(fp, inputNames[iplot2]);
            fwritePlotZLabel(fp, outputNames[jplot]);
            sprintf(winput, "Mesh Plot for %s", outputNames[jplot]);
            fwritePlotTitle(fp, winput);
            fprintf(fp, "a=gca();\n");
            fprintf(fp, "a.data_bounds=[%e,%e;%e,%e];\n",iLowerB[iplot1],
                    iLowerB[iplot2], iUpperB[iplot1], iUpperB[iplot2]);
            fprintf(fp, "a.axes_visible=\"on\";\n");
            fprintf(fp, "drawnow\n");
            fprintf(fp, "subplot(1,2,2)\n");
            fprintf(fp, "end;\n");
            fprintf(fp, "drawlater\n");
            fprintf(fp, "B = A;\n");
            fprintf(fp, "nX = length(x);\n");
            fprintf(fp, "nY = length(y);\n");
            fprintf(fp, "for ii = 1 : nX\n");
            fprintf(fp, "for jj = 1 : nY\n");
            fprintf(fp, "B(ii,jj) = A(nX-ii+1,jj);\n");
            fprintf(fp, "end;\n");
            fprintf(fp, "end;\n");
            fprintf(fp, "a=gca();\n");
            fprintf(fp, "a.data_bounds=[%e,%e;%e,%e];\n",iLowerB[iplot1],
                    iLowerB[iplot2], iUpperB[iplot1], iUpperB[iplot2]);
            fprintf(fp, "bmin = min(min(B)); bmax = max(max(B));\n");
            fprintf(fp, "Matplot1((B-bmin)/(bmax-bmin)*64,[%e,%e,%e,%e])\n",
                    iLowerB[iplot1],iLowerB[iplot2],iUpperB[iplot1], 
                    iUpperB[iplot2]);
            fprintf(fp, "set(gca(),\"auto_clear\",\"off\")\n");
            fprintf(fp, "contour2d(x,y,flipdim(B',1),6);\n");
            fprintf(fp, "xset(\"colormap\",jetcolormap(64));\n");
            fprintf(fp, "colorbar(bmin,bmax);\n");
            fwritePlotAxes(fp);
            fwritePlotXLabel(fp, inputNames[iplot1]);
            fwritePlotYLabel(fp, inputNames[iplot2]);
            sprintf(winput, "Contour Plot for %s", outputNames[jplot]);
            fwritePlotTitle(fp, winput);
            fprintf(fp, "drawnow\n");
            fclose(fp);
            printf("scilabrs2.sci is now available for response surface and ");
            printf("contour plots\n");
         }
         else
         {
            fp = fopen("matlabrs2.m", "w");
            if (fp == NULL)
            {
               printf("ERROR: cannot open file matlabrs2.m.\n");
               continue;
            }
            fprintf(fp, "twoPlots = 1;\n");
            fprintf(fp, "A = [\n");
            for (sInd = 0; sInd < faLeng; sInd++)
               fprintf(fp, "%e\n", faYOut[sInd]);
            fprintf(fp, "];\n");
            fprintf(fp, "A = reshape(A,%d,%d);\n", nPtsPerDim, nPtsPerDim);
            fprintf(fp, "x = [\n");
            for (sInd = 0; sInd < faLeng; sInd+=nPtsPerDim)
               fprintf(fp, "%e\n", faXOut[sInd*2]);
            fprintf(fp, "];\n");
            fprintf(fp, "y = [\n");
            for (sInd = 0; sInd < nPtsPerDim; sInd++)
               fprintf(fp, "%e\n", faXOut[sInd*2+1]);
            fprintf(fp, "];\n");
            fprintf(fp, "B = A;\n");
            Ymin = faYOut[0];
            for (sInd = 1; sInd < faLeng; sInd++)
               if (faYOut[sInd] < Ymin) Ymin = faYOut[sInd];
            Ymax = faYOut[0];
            for (sInd = 1; sInd < faLeng; sInd++)
               if (faYOut[sInd] > Ymax) Ymax = faYOut[sInd];
            printf("Ymin and Ymax found = %e %e.\n", Ymin, Ymax);
            sprintf(pString,"Set lower threshold ? (y or n) : ");
            getString(pString, command);
            fprintf(fp, "n1 = 0;\n");
            fprintf(fp, "n2 = 0;\n");
            if (command[0] == 'y')
            {
               sprintf(pString,"Enter the lower threshold (min = %e) : ",Ymin);
               thresh = getDouble(pString);
               fprintf(fp, "ymin = %e;\n", thresh);
               fprintf(fp, "[ia,ja,aa] = find(A<ymin);\n");
               fprintf(fp, "for ii = 1 : length(ia)\n");
               //fprintf(fp, "   B(ia(ii),ja(ii)) = %e;\n",Ymin-PABS(Ymin)*0.9);
               fprintf(fp, "   B(ia(ii),ja(ii)) = NaN;\n");
               fprintf(fp, "end;\n");
               fprintf(fp, "n1 = length(ia);\n");
            }
            sprintf(pString,"Set upper threshold ? (y or n) : ");
            getString(pString, command);
            if (command[0] == 'y')
            {
               sprintf(pString,"Enter the upper threshold (max = %e) : ",Ymax);
               thresh = getDouble(pString);
               fprintf(fp, "ymax = %e;\n", thresh);
               fprintf(fp, "[ia,ja,aa] = find(A>ymax);\n");
               fprintf(fp, "for ii = 1 : length(ia)\n");
               //fprintf(fp, "   B(ia(ii),ja(ii)) = %e;\n",Ymin-PABS(Ymin)*0.9);
               fprintf(fp, "   B(ia(ii),ja(ii)) = NaN;\n");
               fprintf(fp, "end;\n");
               fprintf(fp, "n2 = length(ia);\n");
            }
            fprintf(fp, "nB = size(B,1);\n");
            fprintf(fp, "if (n1 + n2 == nB * nB)\n");
            fprintf(fp, "   B(1,1) = 0;\n");
            fprintf(fp, "   B(%d,%d) = 1;\n",nPtsPerDim,nPtsPerDim);
            fprintf(fp, "end\n");
            fprintf(fp, "clf\n");
            fprintf(fp, "if twoPlots == 1\n");
            fprintf(fp, "subplot(1,2,1), mesh(x,y,A)\n");
            fwritePlotAxes(fp);
            fwritePlotXLabel(fp, inputNames[iplot1]);
            fwritePlotYLabel(fp, inputNames[iplot2]);
            fwritePlotZLabel(fp, outputNames[jplot]);
            sprintf(winput, "Mesh Plot for %s", outputNames[jplot]);
            fwritePlotTitle(fp, winput);
            fprintf(fp, "colorbar\n");
            fprintf(fp, "subplot(1,2,2)\n");
            fprintf(fp, "end\n");
            fprintf(fp, "contourf(x,y,B)\n");
            fwritePlotAxes(fp);
            fwritePlotXLabel(fp, inputNames[iplot1]);
            fwritePlotYLabel(fp, inputNames[iplot2]);
            fprintf(fp, "colorbar\n");
            fprintf(fp, "colormap(jet)\n");
            sprintf(winput, "Contour Plot for %s", outputNames[jplot]);
            fwritePlotTitle(fp, winput);
            fclose(fp);
            printf("matlabrs2.m is now available for response surface and ");
            printf("contour plots\n");
         }
         delete [] faXOut;
         delete [] faYIn;
         delete [] faYOut;
         delete [] inputSettings;
         delete faPtr;
         faPtr = NULL;
      }

      // generate response surface of any 3 inputs and write the
      // grid data to file for display with matlab
      // ----------------------------------------------------------------
      else if (!strcmp(command, "rs3") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rs3: response surface plot in three parameters\n");
            printf("syntax: rs3 (no argument needed)\n");
            continue;
         }
         printf("ERROR: no data to analyze (load data first).\n");
      }
      else if (!strcmp(command, "rs3") && nInputs < 3)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rs3: response surface plot in three parameters\n");
            printf("syntax: rs3 (no argument needed)\n");
            continue;
         }
         printf("ERROR: rs3 requires 3 or more inputs.\n");
      }
      else if (!strcmp(command, "rs3") && nInputs >= 3)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rs3: response surface plot in three parameters\n");
            printf("syntax: rs3 (no argument needed)\n");
            continue;
         }
         if (psPlotTool_ == 1)
         {
            printf("INFO: rs3 is currently not available in scilab.\n");
            continue;
         }
         nPtsPerDim = 16;
         sprintf(pString, "Grid resolution ? (16 - 32) ");
         nPtsPerDim = getInt(16, 32, pString);
         faFlag = 1;
         faPtr = genFAInteractive(psuadeIO_, faFlag);
         if (faPtr == NULL) {printf("ERROR detected.\n"); continue;}
         faPtr->setNPtsPerDim(nPtsPerDim);
         faPtr->setBounds(iLowerB, iUpperB);
         faPtr->setOutputLevel(outputLevel_);

         inputSettings = new double[nInputs];
         iplot1 = iplot2 = iplot3 = -1;
         sprintf(pString, "Enter the input for x axis (1 - %d) : ", nInputs);
         iplot1 = getInt(1, nInputs, pString);
         iplot1--;
         iplot2 = iplot1;
         while (iplot1 == iplot2)
         {
            sprintf(pString, "Enter the input for y axis (1 - %d), not %d : ",
                    nInputs, iplot1+1);
            iplot2 = getInt(1, nInputs, pString);
            iplot2--;
            if (iplot1 == iplot2)
               printf("ERROR: duplicate input number %d.\n",iplot2+1);
         }
         if (nInputs == 3) iplot3 = 3 - iplot1 - iplot2;
         while (iplot3 < 0 || iplot3 == iplot1 || iplot3 == iplot2)
         {
            sprintf(pString,
                    "Enter the input for z axis (1 - %d), not %d nor %d: ",
                    nInputs, iplot1+1, iplot2+1);
            iplot3 = getInt(1, nInputs, pString);
            iplot3--;
            if (iplot3 == iplot1 || iplot3 == iplot2)
               printf("ERROR: duplicate input number %d.\n",iplot3+1);
         }
         strcpy(command, "y\n");
         if (nInputs > 3)
         {
            sprintf(pString,"Set other nominal values automatically? (y/n) ");
            getString(pString, command);
         }
         if (command[0] == 'y')
         {
            for (iInd1 = 0; iInd1 < nInputs; iInd1++)
            {
               if (iInd1 != iplot1 && iInd1 != iplot2 && iInd1 != iplot3)
                    inputSettings[iInd1] = 0.5*(iLowerB[iInd1]+iUpperB[iInd1]);
               else inputSettings[iInd1] = 1.0;
            }
         }
         else
         {
            for (iInd1 = 0; iInd1 < nInputs; iInd1++)
            {
               if (iInd1 != iplot1 && iInd1 != iplot2 && iInd1 != iplot3)
               {
                  inputSettings[iInd1] = iLowerB[iInd1] - 1.0;
                  sprintf(pString,
                          "Enter nominal value for input %d (%e - %e): ", 
                          iInd1+1, iLowerB[iInd1], iUpperB[iInd1]);
                  while (inputSettings[iInd1] < iLowerB[iInd1] ||
                         inputSettings[iInd1] > iUpperB[iInd1])
                     inputSettings[iInd1] = getDouble(pString);
               }
               else inputSettings[iInd1] = 1.0;
            }
         }
         jplot = 0;
         if (nOutputs > 1)
         {
            sprintf(pString, "Enter the output number (1 - %d) : ", nOutputs);
            jplot = getInt(1, nOutputs, pString);
            jplot--;
         }
         faYIn = new double[nSamples];
         for (sInd = 0; sInd < nSamples; sInd++)
            faYIn[sInd] = sampleOutputs[sInd*nOutputs+jplot];

         printf("Please wait while generating the RS data \n");
         faPtr->gen3DGridData(sampleInputs,faYIn, iplot1, iplot2, iplot3, 
                              inputSettings, &faLeng, &faXOut,&faYOut);

         GYmin = faYOut[0];
         for (sInd = 1; sInd < faLeng; sInd++)
            if (faYOut[sInd] < GYmin) GYmin = faYOut[sInd];
         GYmax = faYOut[0];
         for (sInd = 1; sInd < faLeng; sInd++)
            if (faYOut[sInd] > GYmax) GYmax = faYOut[sInd];
         printf("\nYmin and Ymax found = %e %e.\n", GYmin, GYmax);
         threshL = GYmin - 0.2 * PABS(GYmax-GYmin);
         gamma = threshL;
         sprintf(pString,"Set lower threshold? (y or n) ");
         getString(pString, command);
         if (command[0] == 'y')
         {
            sprintf(pString,"Enter the lower threshold (min = %e): ",GYmin);
            threshL = getDouble(pString);
            if (threshL < GYmin)
            {
               threshL = GYmin;
               printf("rs3 INFO: lower threshold set to %e.\n", threshL);
            }
         }
         threshU = GYmax + 0.2 * PABS(GYmax-GYmin);
         sprintf(pString,"Set upper threshold? (y or n) ");
         getString(pString, command);
         if (command[0] == 'y')
         {
            sprintf(pString,"Enter the upper threshold (max = %e): ",GYmax);
            threshU = getDouble(pString);
            if (threshU > GYmax)
            {
               threshU = GYmax;
               printf("rs3 INFO: upper threshold set to %e.\n", threshU);
            }
         }
         if (threshL >= threshU)
         {
            printf("rs3 ERROR: lower threshold (%e) >= upper threshold (%e)\n",
                   threshL, threshU);
            continue;
         }
         fp = fopen("matlabrs3.m", "w");
         if (fp == NULL)
         {
            printf("ERROR: cannot open file matlabrs3.m.\n");
            continue;
         }
         fwritePlotCLF(fp);
         fprintf(fp, "xlo = %e; \n", iLowerB[iplot2]);
         fprintf(fp, "xhi = %e; \n", iUpperB[iplot2]);
         fprintf(fp, "ylo = %e; \n", iLowerB[iplot1]);
         fprintf(fp, "yhi = %e; \n", iUpperB[iplot1]);
         fprintf(fp, "zlo = %e; \n", iLowerB[iplot3]);
         fprintf(fp, "zhi = %e; \n", iUpperB[iplot3]);
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
                  if (faYOut[ind] < threshL)
                  {
                     fprintf(fp, "%e ", gamma);
                     count++;
                  }
                  else if (faYOut[ind] > threshU)
                  {
                     fprintf(fp, "%e ", gamma);
                     count++;
                  }
                  else fprintf(fp, "%e ", faYOut[ind]);
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
         fprintf(fp, "xt = [%e:%e:%e];\n", iLowerB[iplot2],
                 (iUpperB[iplot2]-iLowerB[iplot2])*0.01, iUpperB[iplot2]);
         fprintf(fp, "yt = [%e:%e:%e];\n", iLowerB[iplot1],
                 (iUpperB[iplot1]-iLowerB[iplot1])*0.01, iUpperB[iplot1]);
         fprintf(fp, "zt = [%e:%e:%e];\n", iLowerB[iplot3],
                 (iUpperB[iplot3]-iLowerB[iplot3])*0.01, iUpperB[iplot3]);
         fprintf(fp, "isoval = %e;\n", gamma);
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
         fprintf(fp, "daspect([%e,%e,%e])\n",iUpperB[iplot2]-iLowerB[iplot2],
                 iUpperB[iplot1]-iLowerB[iplot1],
                 iUpperB[iplot3]-iLowerB[iplot3]);
         fprintf(fp, "   xlabel('%s','FontSize',12,'FontWeight','bold')\n",
                 inputNames[iplot2]);
         fprintf(fp, "   ylabel('%s','Fontsize',12,'FontWeight','bold')\n",
                 inputNames[iplot1]);
         fprintf(fp, "   zlabel('%s','Fontsize',12,'FontWeight','bold')\n",
                 inputNames[iplot3]);
         fprintf(fp, "   title('%s','Fontsize',12,'FontWeight','bold')\n", 
                 outputNames[jplot]);
         fwritePlotAxes(fp);
         fprintf(fp, "colormap('default'); colorbar\n");
         fprintf(fp, "%%axis tight\n");
         fprintf(fp, "view(3) \n");
         fprintf(fp, "set(gcf,'Renderer','zbuffer')\n");
         fprintf(fp, "lighting phong\n");
         fprintf(fp, "cin = input('generate slices ? (y or n) ','s');\n");
         fprintf(fp, "if (cin == 'y')\n");
         fprintf(fp, "xin = input('axis to slide through ? (x,y,z) ','s');\n");
         fprintf(fp, "for i = 1 : 101\n");
         fprintf(fp, "   display(['displaying ' int2str(i) ' of 100'])\n");
         fprintf(fp, "   if (xin == 'y')\n");
         fprintf(fp, "      h = slice(X,Y,Z,V,xt(i),[],[]);\n");
         fprintf(fp, "   elseif (xin == 'x')\n");
         fprintf(fp, "      h = slice(X,Y,Z,V,[],yt(i),[]);\n");
         fprintf(fp, "   elseif (xin == 'z')\n");
         fprintf(fp, "      h = slice(X,Y,Z,V,[],[],zt(i));\n");
         fprintf(fp, "   end\n");
         fprintf(fp, "   axis([%11.4e %11.4e %11.4e %11.4e %11.4e %11.4e ",
                 iLowerB[iplot2], iUpperB[iplot2], iLowerB[iplot1],
                 iUpperB[iplot1], iLowerB[iplot3], iUpperB[iplot3]);
         fprintf(fp, "%11.4e %11.4e])\n",
                 threshL-0.2*(threshU-threshL),threshU+0.2*(threshU-threshL));
         fwritePlotAxes(fp);
         fprintf(fp, "   xlabel('%s','FontSize',12,'FontWeight','bold')\n",
                 inputNames[iplot2]);
         fprintf(fp, "   ylabel('%s','Fontsize',12,'FontWeight','bold')\n",
                 inputNames[iplot1]);
         fprintf(fp, "   zlabel('%s','Fontsize',12,'FontWeight','bold')\n",
                 inputNames[iplot3]);
         fprintf(fp, "   title('3D Contour Plot',");
         fprintf(fp, "'FontWeight','bold','FontSize',12)\n");
         fprintf(fp, "   view(3)\n");
         fprintf(fp, "   colorbar\n");
         fprintf(fp, "   pause(1)\n");
         fprintf(fp, "   if (i < 101)\n");
         fprintf(fp, "      clf\n");
         fprintf(fp, "   end\n");
         fprintf(fp, "end\n");
         fprintf(fp, "end\n");
         fclose(fp);
         printf("\nmatlabrs3.m is now available.\n");
         delete [] faYIn;
         delete faPtr;
         faPtr = NULL;
         delete [] inputSettings;
         delete [] faXOut;
         delete [] faYOut;
      }

      // generate response surface of any 3 inputs and write the
      // grid data to file for display with matlab (movie)
      // ----------------------------------------------------------------
      else if (!strcmp(command, "rs3m") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rs3m: response surface plot in 3 parameters using movie mode\n");
            printf("syntax: rs3m (no argument needed)\n");
            continue;
         }
         printf("ERROR: no data to analyze (load data first).\n");
      }
      else if (!strcmp(command, "rs3m") && nInputs < 3)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rs3m: response surface plot in 3 parameters using movie mode\n");
            printf("syntax: rs3m (no argument needed)\n");
            continue;
         }
         printf("ERROR: rs3m requires 3 or more inputs.\n");
      }
      else if (!strcmp(command, "rs3m") && nInputs >= 3)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rs3m: response surface plot in 3 parameters using movie mode\n");
            printf("syntax: rs3m (no argument needed)\n");
            continue;
         }
         if (psPlotTool_ == 1)
         {
            printf("INFO: rs3m is currently not available in scilab.\n");
            continue;
         }
         nPtsPerDim = 16;
         sprintf(pString, "Grid resolution ? (16 - 32) ");
         nPtsPerDim = getInt(16, 32, pString);
         faFlag = 1;
         faPtr = genFAInteractive(psuadeIO_, faFlag);
         if (faPtr == NULL) {printf("ERROR detected.\n"); continue;}
         faPtr->setNPtsPerDim(nPtsPerDim);
         faPtr->setBounds(iLowerB, iUpperB);
         faPtr->setOutputLevel(outputLevel_);

         inputSettings = new double[nInputs];
         iplot1 = iplot2 = iplot3 = -1;
         sprintf(pString, "Enter the input for x axis (1 - %d) : ", nInputs);
         iplot1 = getInt(1, nInputs, pString);
         iplot1--;
         iplot2 = iplot1;
         while (iplot1 == iplot2)
         {
            sprintf(pString, "Enter the input for y axis (1 - %d), not %d : ",
                    nInputs, iplot1+1);
            iplot2 = getInt(1, nInputs, pString);
            iplot2--;
            if (iplot1 == iplot2)
               printf("ERROR: duplicate input number.\n");
         }
         if (nInputs == 3) iplot3 = 3 - iplot1 - iplot2;
         while (iplot3 < 0 || iplot3 == iplot1 || iplot3 == iplot2)
         {
            sprintf(pString,
                    "Enter the input for t axis (1 - %d), not %d nor %d: ",
                    nInputs, iplot1+1, iplot2+1);
            iplot3 = getInt(1, nInputs, pString);
            iplot3--;
            if (iplot3 == iplot1 || iplot3 == iplot2)
               printf("ERROR: duplicate input number %d.\n", iplot3+1);
         }
         if (nInputs > 3)
         {
            sprintf(pString,
                    "Set other nominal values at mid point ? (y or n) ");
            getString(pString, command);
            if (command[0] == 'y')
            {
               for (iInd1 = 0; iInd1 < nInputs; iInd1++)
               {
                  if (iInd1 != iplot1 && iInd1 != iplot2 && iInd1 != iplot3)
                     inputSettings[iInd1] = 0.5*(iLowerB[iInd1]+iUpperB[iInd1]);
                  else inputSettings[iInd1] = 1.0;
               }
            }
            else
            {
               for (iInd1 = 0; iInd1 < nInputs; iInd1++)
               {
                  if (iInd1 != iplot1 && iInd1 != iplot2 && iInd1 != iplot3)
                  {
                     sprintf(pString,"Enter setting for input %d (%e - %e): ", 
                             iInd1+1, iLowerB[iInd1], iUpperB[iInd1]);
                     inputSettings[iInd1] = getDouble(pString);
                  }
                  else inputSettings[iInd1] = 1.0;
               }
            }
         }
         jplot = 0;
         if (nOutputs > 1)
         {
            sprintf(pString,"Enter the output number (1 - %d) : ", nOutputs);
            jplot = getInt(1, nOutputs, pString);
            jplot--;
         }
         faYIn = new double[nSamples];
         for (sInd = 0; sInd < nSamples; sInd++)
            faYIn[sInd] = sampleOutputs[sInd*nOutputs+jplot];

         fp = fopen("matlabrs3m.m", "w");
         if (fp == NULL)
         {
            printf("ERROR: cannot open file matlabrs3m.m.\n");
            continue;
         }

         faPtr->gen3DGridData(sampleInputs,faYIn, iplot1, iplot2, iplot3, 
                              inputSettings, &faLeng, &faXOut,&faYOut);
         GYmin = faYOut[0];
         for (sInd = 1; sInd < faLeng; sInd++)
            if (faYOut[sInd] < GYmin) GYmin = faYOut[sInd];
         Ymax = faYOut[0];
         for (sInd = 1; sInd < faLeng; sInd++)
            if (faYOut[sInd] > GYmax) GYmax = faYOut[sInd];
         printf("\nYmin and Ymax found = %e %e.\n", GYmin, GYmax);
         threshL = GYmin - 0.2 * PABS(GYmin);
         sprintf(pString, "Set lower threshold? (y or n) ");
         getString(pString, command);
         if (command[0] == 'y')
         {
            sprintf(pString,"Enter the lower threshold (min = %e): ",GYmin);
            threshL = getDouble(pString);
         }
         threshU = GYmax + 0.2 * PABS(GYmax);
         sprintf(pString, "Set upper threshold? (y or n) ");
         getString(pString, command);
         if (command[0] == 'y')
         {
            sprintf(pString,"Enter the upper threshold (min = %e): ",GYmax);
            threshU = getDouble(pString);
         }
         fprintf(fp, "twoPlots = 1;\n");
         fprintf(fp, "disp(\'Please wait while loading data.\')\n");
         fprintf(fp, "hold off\n");
         fprintf(fp, "clf\n");

         for (ii = 0; ii < nPtsPerDim; ii++)
         {
            inputSettings[iplot3] = (iUpperB[iplot3] - iLowerB[iplot3]) / 
                                    (nPtsPerDim - 1.0) * ii + iLowerB[iplot3];
            fprintf(fp, "x = [\n");
            for (sInd = 0; sInd < faLeng; sInd+=nPtsPerDim*nPtsPerDim)
               fprintf(fp, "%e\n", faXOut[sInd*3]);
            fprintf(fp, "];\n");
            fprintf(fp, "y = [\n");
            for (sInd = 0; sInd < nPtsPerDim*nPtsPerDim; sInd+=nPtsPerDim)
               fprintf(fp, "%e\n", faXOut[sInd*3+1]);
            fprintf(fp, "];\n");

            fprintf(fp, "A%d = [\n", ii + 1);
            for (sInd = 0; sInd < faLeng; sInd+=nPtsPerDim)
               fprintf(fp, "%e\n", faYOut[sInd+ii]);
            fprintf(fp, "];\n");
            fprintf(fp, "A%d = reshape(A%d,%d,%d);\n", ii+1, ii+1,
                    nPtsPerDim, nPtsPerDim);
            fprintf(fp, "disp(\'Plotting frame %d of %d\')\n",ii+1,nPtsPerDim);
            fprintf(fp, "B%d = A%d;\n", ii+1, ii+1);
            fprintf(fp, "yLo = %e;\n", threshL);
            fprintf(fp, "yHi = %e;\n", threshU);
            fprintf(fp, "nA  = size(A%d,1);\n", ii+1);
            fprintf(fp, "[ia,ja,aa] = find(A%d<yLo);\n", ii+1);
            fprintf(fp, "for ii = 1 : length(ia)\n");
            fprintf(fp, "   B%d(ia(ii),ja(ii)) = NaN;\n", ii+1);
            fprintf(fp, "end;\n");
            fprintf(fp, "n1 = length(ia);\n");
            fprintf(fp, "[ia,ja,aa] = find(A%d>yHi);\n", ii+1);
            fprintf(fp, "for ii = 1 : length(ia)\n");
            fprintf(fp, "   B%d(ia(ii),ja(ii)) = NaN;\n", ii+1);
            fprintf(fp, "end;\n");
            fprintf(fp, "n2 = length(ia);\n");
            fprintf(fp, "if (n1 + n2 == nA*nA)\n");
            fprintf(fp, "   B%d(1,1) = 0;\n",ii+1);
            fprintf(fp, "   B%d(%d,%d) = 1;\n",ii+1,nPtsPerDim,
                    nPtsPerDim);
            fprintf(fp, "end;\n");
            fprintf(fp, "if twoPlots == 1\n");
            fprintf(fp, "subplot(1,2,1), surf(x,y,A%d)\n", ii+1);
            fprintf(fp, "axis([%e %e %e %e %e %e])\n",iLowerB[iplot1],
                    iUpperB[iplot1],iLowerB[iplot2],iUpperB[iplot2],
                    GYmin, GYmax); 
            fwritePlotAxes(fp);
            fprintf(fp, "xlabel('%s','FontSize',12,'FontWeight','bold')\n",
                    inputNames[iplot1]);
            fprintf(fp, "ylabel('%s','Fontsize',12,'FontWeight','bold')\n",
                    inputNames[iplot2]);
            fprintf(fp, "zlabel('%s','Fontsize',12,'FontWeight','bold')\n",
                    outputNames[jplot]);
            fprintf(fp, "colorbar\n");
            fprintf(fp, "title(\'%s Mesh plot, val(3) = %14.7e\',",
                    outputNames[jplot], inputSettings[iplot3]);
            fprintf(fp, "'FontWeight','bold','FontSize',12)\n");
            fprintf(fp, "subplot(1,2,2)\n");
            fprintf(fp, "end\n");
            fprintf(fp, "contourf(x,y,B%d)\n",ii+1);
            fprintf(fp, "axis([%e %e %e %e])\n",iLowerB[iplot1],
                    iUpperB[iplot1],iLowerB[iplot2],iUpperB[iplot2]);
            fwritePlotAxes(fp);
            fprintf(fp, "colorbar\n");
            fprintf(fp, "colormap(jet)\n");
            fprintf(fp, "title(\'%s contour plot, val(3) = %14.7e\',",
                        outputNames[jplot], inputSettings[iplot3]);
            fprintf(fp, "'FontWeight','bold','FontSize',12)\n");
            fprintf(fp,"pause(1)\n");
         }
         delete [] faXOut;
         delete [] faYOut;
         fprintf(fp, "rotate3d on\n");
         fclose(fp);
         printf("matlabrs3m.m is now available for response surface and ");
         printf("contour plots\n");
         delete [] faYIn;
         delete faPtr;
         faPtr = NULL;
         delete [] inputSettings;
      }

      // generate response surface of any 4 inputs and write the
      // grid data to file for display with matlab
      // ----------------------------------------------------------------
      else if (!strcmp(command, "rs4") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rs4: response surface plot in 4 parameters\n");
            printf("syntax: rs4 (no argument needed)\n");
            continue;
         }
         printf("ERROR: no data to analyze.\n");
      }
      else if (!strcmp(command, "rs4") && nInputs < 4)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rs4: response surface plot in 4 parameters\n");
            printf("syntax: rs4 (no argument needed)\n");
            continue;
         }
         printf("ERROR: rs4 requires 4 or more inputs.\n");
      }
      else if (!strcmp(command, "rs4") && nInputs >= 4)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rs4: response surface plot in 4 parameters\n");
            printf("syntax: rs4 (no argument needed)\n");
            continue;
         }
         if (psPlotTool_ == 1)
         {
            printf("INFO: rs4 is currently not available in scilab.\n");
            continue;
         }
         nPtsPerDim = 16;
         printf("Note: if matlab crashes, it may be due to high grid resolution.\n");
         sprintf(pString, "Grid resolution ? (16 - 32) ");
         nPtsPerDim = getInt(16, 32, pString);
         faFlag = 1;
         faPtr = genFAInteractive(psuadeIO_, faFlag);
         if (faPtr == NULL) {printf("ERROR detected.\n"); continue;}
         faPtr->setNPtsPerDim(nPtsPerDim);
         faPtr->setBounds(iLowerB, iUpperB);
         faPtr->setOutputLevel(outputLevel_);

         inputSettings = new double[nInputs];
         iplot1 = iplot2 = iplot3 = iplot4 = -1;
         sprintf(pString, "Enter the input for x axis (1 - %d) : ", nInputs);
         iplot1 = getInt(1, nInputs, pString);
         iplot1--;
         iplot2 = iplot1;
         while (iplot1 == iplot2)
         {
            sprintf(pString, "Enter the input for y axis (1 - %d), not %d : ",
                    nInputs, iplot1+1);
            iplot2 = getInt(1, nInputs, pString);
            iplot2--;
            if (iplot1 == iplot2)
               printf("ERROR: duplicate input number %d.\n",iplot2+1);
         }
         iplot3 = iplot1;
         while (iplot3 == iplot1 || iplot3 == iplot2)
         {
            sprintf(pString, "Enter the input for z axis (1 - %d), not %d,%d: ",
                    nInputs, iplot1+1, iplot2+1);
            iplot3 = getInt(1, nInputs, pString);
            iplot3--;
            if (iplot3 == iplot1 || iplot3 == iplot2)
               printf("ERROR: duplicate input number %d.\n",iplot3+1);
         }
         if (nInputs == 4) iplot4 = 6 - iplot1 - iplot2 - iplot3;
         while (iplot4 < 0 || iplot4 == iplot1 || iplot4 == iplot2 || 
                iplot4 == iplot3)
         {
            sprintf(pString,
                    "Enter the input for t axis (1 - %d), not %d nor %d,%d: ",
                    nInputs, iplot1+1, iplot2+1, iplot3+1);
            iplot4 = getInt(1, nInputs, pString);
            iplot4--;
            if (iplot4 == iplot1 || iplot4 == iplot2 || iplot4 == iplot3)
               printf("ERROR: duplicate input number %d.\n",iplot4+1);
         }
         strcpy(command, "y\n");
         if (nInputs > 4)
         {
            sprintf(pString,"Set other nominal values at mid point ? (y/n) ");
            getString(pString, command);
         }
         if (command[0] == 'y')
         {
            for (iInd1 = 0; iInd1 < nInputs; iInd1++)
            {
               if (iInd1 != iplot1 && iInd1 != iplot2 && iInd1 != iplot3 &&
                   iInd1 != iplot4)
                    inputSettings[iInd1] = 0.5*(iLowerB[iInd1]+iUpperB[iInd1]);
               else inputSettings[iInd1] = 1.0;
            }
         }
         else
         {
            for (iInd1 = 0; iInd1 < nInputs; iInd1++)
            {
               if (iInd1 != iplot1 && iInd1 != iplot2 && iInd1 != iplot3 &&
                   iInd1 != iplot4)
               {
                  inputSettings[iInd1] = iLowerB[iInd1] - 1.0;
                  sprintf(pString,
                          "Enter nominal value for input %d (%e - %e): ", 
                          iInd1+1, iLowerB[iInd1], iUpperB[iInd1]);
                  while (inputSettings[iInd1] < iLowerB[iInd1] ||
                         inputSettings[iInd1] > iUpperB[iInd1])
                     inputSettings[iInd1] = getDouble(pString);
               }
               else inputSettings[iInd1] = 1.0;
            }
         }
         jplot = 0;
         if (nOutputs > 1)
         {
            sprintf(pString, "Enter the output number (1 - %d) : ", nOutputs);
            jplot = getInt(1, nOutputs, pString);
            jplot--;
         }
         faYIn = new double[nSamples];
         for (sInd = 0; sInd < nSamples; sInd++)
            faYIn[sInd] = sampleOutputs[sInd*nOutputs+jplot];

         faPtr->gen4DGridData(sampleInputs,faYIn, iplot1, iplot2, iplot3, 
                              iplot4, inputSettings, &faLeng, &faXOut,&faYOut);
         GYmin =   PSUADE_UNDEFINED;
         GYmax = - PSUADE_UNDEFINED;
         for (sInd = 0; sInd < faLeng; sInd++)
            if (faYOut[sInd] < GYmin) GYmin = faYOut[sInd];
         for (sInd = 0; sInd < faLeng; sInd++)
            if (faYOut[sInd] > GYmax) GYmax = faYOut[sInd];
         printf("\nYmin and Ymax found = %e %e.\n", GYmin, GYmax);
         threshL = GYmin - 0.2 * PABS(GYmax - GYmin);
         sprintf(pString,"Set lower threshold? (y or n) ");
         gamma = threshL;
         getString(pString, command);
         if (command[0] == 'y')
         {
            sprintf(pString,"Enter the lower threshold (min = %e): ",GYmin);
            threshL = getDouble(pString);
         }
         threshU = GYmax + 0.2 * PABS(GYmax - GYmin);
         sprintf(pString,"Set upper threshold? (y or n) ");
         getString(pString, command);
         if (command[0] == 'y')
         {
            sprintf(pString,"Enter the upper threshold (max = %e): ",GYmax);
            threshU = getDouble(pString);
         }

         fp = fopen("matlabrs4.m", "w");
         if (fp == NULL)
         {
            printf("ERROR: cannot open file matlabrs4.m.\n");
            continue;
         }

         fprintf(fp, "%% user adjustable parameter section begins *****n");
         fprintf(fp, "%% use nSubplots, nSubNx and nSubNy to spread \n");
         fprintf(fp, "%% the movie frames into a number of subplots.\n");
         fprintf(fp, "nSubplots = 1;\n");
         fprintf(fp, "nSubNx = 1;\n");
         fprintf(fp, "nSubNy = 1;\n");
         fprintf(fp, "%% user adjustable parameter section ends *****\n");
         fprintf(fp, "clf\n");
         fprintf(fp, "nFrames = %d;\n", nPtsPerDim);
         fprintf(fp, "nSubCnt = 0;\n");
         fprintf(fp, "isoval = %e;\n", threshL);
         fprintf(fp, "xlo = %e; \n", iLowerB[iplot2]);
         fprintf(fp, "xhi = %e; \n", iUpperB[iplot2]);
         fprintf(fp, "ylo = %e; \n", iLowerB[iplot1]);
         fprintf(fp, "yhi = %e; \n", iUpperB[iplot1]);
         fprintf(fp, "zlo = %e; \n", iLowerB[iplot3]);
         fprintf(fp, "zhi = %e; \n", iUpperB[iplot3]);
         fprintf(fp, "X=zeros(%d,%d,%d);\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);
         fprintf(fp, "Y=zeros(%d,%d,%d);\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);
         fprintf(fp, "Z=zeros(%d,%d,%d);\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);
         fprintf(fp, "V=zeros(%d,%d,%d);\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);
         for (ii = 0; ii < nPtsPerDim; ii++)
         {
            fprintf(fp, "Y(:,:,%d) = [\n", ii + 1);
            for (sInd = 0; sInd < nPtsPerDim; sInd++)
            {
               for (jj = 0; jj < nPtsPerDim; jj++)
               {
                  ind = (sInd*nPtsPerDim*nPtsPerDim+jj*nPtsPerDim+ii)*nPtsPerDim;
                  fprintf(fp, "%e ", faXOut[ind*4]);
               }
               fprintf(fp, "\n");
            }
            fprintf(fp, "];\n");
            fprintf(fp, "X(:,:,%d) = [\n", ii + 1);
            for (sInd = 0; sInd < nPtsPerDim; sInd++)
            {
               for (jj = 0; jj < nPtsPerDim; jj++)
               {
                  ind = (sInd*nPtsPerDim*nPtsPerDim+jj*nPtsPerDim+ii)*nPtsPerDim;
                  fprintf(fp, "%e ", faXOut[ind*4+1]);
               }
               fprintf(fp, "\n");
            }
            fprintf(fp, "];\n");
            fprintf(fp, "Z(:,:,%d) = [\n", ii + 1);
            for (sInd = 0; sInd < nPtsPerDim; sInd++)
            {
               for (jj = 0; jj < nPtsPerDim; jj++)
               {
                  ind = (sInd*nPtsPerDim*nPtsPerDim+jj*nPtsPerDim+ii)*nPtsPerDim;
                  fprintf(fp, "%e ", faXOut[ind*4+2]);
               }
               fprintf(fp, "\n");
            }
            fprintf(fp, "];\n");
         }
         for (ll = 0; ll < nPtsPerDim; ll++)
         {
            for (ii = 0; ii < nPtsPerDim; ii++)
            {
               count = 0;
               fprintf(fp, "V(:,:,%d) = [\n", ii + 1);
               for (sInd = 0; sInd < nPtsPerDim; sInd++)
               {
                  for (jj = 0; jj < nPtsPerDim; jj++)
                  {
                     ind = ((sInd*nPtsPerDim+jj)*nPtsPerDim+ii)*nPtsPerDim+ll;
                     if (faYOut[ind] < threshL)
                     {
                        fprintf(fp, "%e ", threshL);
                        count++;
                     }
                     else if (faYOut[ind] > threshU)
                     {
                        fprintf(fp, "%e ", threshL);
                        count++;
                     }
                     else fprintf(fp, "%e ", faYOut[ind]);
                  }
                  fprintf(fp, "\n");
               }
               fprintf(fp, "];\n");
               if (count == nPtsPerDim*nPtsPerDim)
               {
                  if (threshL-0.2*(threshU-threshL) > gamma)
                     fprintf(fp, "V(:,:,%d) = %e * ones(%d,%d);\n",ii+1,gamma,
                             nPtsPerDim, nPtsPerDim);
                  else
                     fprintf(fp, "V(:,:,%d) = %e * ones(%d,%d);\n",ii+1,
                             threshL-0.2*(threshU-threshL),
                             nPtsPerDim,nPtsPerDim);
                  printf("Frame %d, slice %d nonfeasible -> set to ground.\n",
                         ll+1, ii+1);
               }
            }
            fprintf(fp, "frame = %d;\n", ll+1);
            fprintf(fp, "if nSubplots > 1\n");
            fprintf(fp, "   if frame <= 2\n");
            fprintf(fp, "      nSubCnt = nSubCnt + 1;\n");
            fprintf(fp, "      subplot(nSubNx, nSubNy, nSubCnt)\n");
            fprintf(fp, "   elseif frame == nFrames\n");
            fprintf(fp, "      subplot(nSubNx, nSubNy, nSubplots)\n");
            fprintf(fp, "   else\n");
            fprintf(fp, "      ft1 = (nFrames-1) / (nSubplots-1);\n");
            fprintf(fp, "      ft2 = round(ft1 * (nSubCnt-1)) + 2;\n");
            fprintf(fp, "      if frame == ft2\n");
            fprintf(fp, "         nSubCnt = nSubCnt + 1;\n");
            fprintf(fp, "         subplot(nSubNx, nSubNy, nSubCnt)\n");
            fprintf(fp, "      end\n");
            fprintf(fp, "   end\n");
            fprintf(fp, "else\n");
            fprintf(fp, "   clf\n");
            fprintf(fp, "end\n");
            fprintf(fp, "disp('Frame %d of %d')\n", ll+1, nPtsPerDim);
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
            fprintf(fp, "daspect([xhi-xlo, yhi-ylo, zhi-zlo])\n");
            fprintf(fp, "colormap('default')\n");
            fprintf(fp, "if nSubplots > 1\n");
            fprintf(fp, "   colorbar\n");
            fprintf(fp, "elseif nSubCnt == 2\n");
            fprintf(fp, "   colorbar\n");
            fprintf(fp, "end\n");
            fprintf(fp, "%%axis tight\n");
            fprintf(fp, "view(3) \n");
            fprintf(fp, "set(gcf,'Renderer','zbuffer')\n");
            fprintf(fp, "box on\n");
            fprintf(fp, "grid on\n");
            fprintf(fp, "lighting phong\n");
            fwritePlotAxes(fp);
            fprintf(fp, "xlabel('%s','FontSize',12,'FontWeight','bold')\n",
                    inputNames[iplot2]);
            fprintf(fp, "ylabel('%s','Fontsize',12,'FontWeight','bold')\n",
                    inputNames[iplot1]);
            fprintf(fp, "zlabel('%s','Fontsize',12,'FontWeight','bold')\n",
                    inputNames[iplot3]);
            fprintf(fp, "title('3D Plot at %s=%e',",
                    inputNames[iplot4],faXOut[ll*4+3]);
            fprintf(fp, "'FontWeight','bold','FontSize',12)\n");
            fprintf(fp, "pause(1)\n");
         }
         delete [] faXOut;
         delete [] faYOut;
         fclose(fp);
         printf("\nmatlabrs4.m is now available.\n");
         delete [] faYIn;
         delete faPtr;
         faPtr = NULL;
         delete [] inputSettings;
      }

      // generate standard deviation response surface and write the
      // grid data to file for display with matlab
      // ----------------------------------------------------------------
      else if (!strcmp(command, "rssd") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rssd: response surface plots for the std. deviations.\n");
            printf("INFO: rssd not available for >2 inputs for scilab.\n");
            printf("syntax: rssd (no argument needed.\n");
            continue;
         }
         printf("ERROR: no data to analyze (load data first).\n");
      }
      else if (!strcmp(command, "rssd"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rssd: response surface plots for the std. deviations.\n");
            printf("INFO: rssd not available for >2 inputs for scilab.\n");
            printf("syntax: rssd (no argument needed.\n");
            continue;
         }
         inputSettings = new double[nInputs];
         iplot1 = iplot2 = iplot3 = iplot4 = -1;
         sprintf(pString, "Enter the input for x axis (1 - %d) : ", nInputs);
         iplot1 = getInt(1, nInputs, pString);
         iplot1--;
         count = 1;
         if (nInputs == 2) iplot2 = nInputs - iplot1 - 1;
         else if (nInputs > 2)
         {
            iplot2 = iplot1;
            while (iplot1 == iplot2)
            {
               sprintf(pString, "Y-axis input ? (1-%d, 0 if not used, not %d) ",
                       nInputs, iplot1+1);
               iplot2 = getInt(0, nInputs, pString);
               iplot2--;
               if (iplot2 == -1) break;
               if (iplot1 == iplot2)
                  printf("ERROR: duplicate input number %d.\n",iplot2+1);
            }
         }
         if (iplot2 != -1) count++;
         if (psPlotTool_ == 0 && iplot2 != -1)
         {
            if (nInputs == 3) iplot3 = nInputs - iplot1 - iplot2;
            else if (nInputs > 3)
            {
               iplot3 = iplot1;
               while (iplot3 == iplot1 || iplot3 == iplot2)
               {
                  sprintf(pString,
                     "Z axis input ? (1-%d, 0 if not used, not %d nor %d) ",
                     nInputs, iplot1+1, iplot2+1);
                  iplot3 = getInt(0, nInputs, pString);
                  iplot3--;
                  if (iplot3 == -1) break;
                  if (iplot3 == iplot1 || iplot3 == iplot2)
                     printf("ERROR: duplicate input number %d.\n",iplot3+1);
               }
            }
            if (iplot3 != -1) count++;
            if (nInputs >= 4 && iplot3 != -1)
            {
               while (iplot4 < 0 || iplot4 == iplot1 || iplot4 == iplot2 || 
                      iplot4 == iplot3)
               {
                  sprintf(pString,
                     "Enter the input for t axis (1 - %d), not %d nor %d,%d: ",
                          nInputs, iplot1+1, iplot2+1, iplot3+1);
                  iplot4 = getInt(1, nInputs, pString);
                  iplot4--;
                  if (iplot4 == iplot1 || iplot4 == iplot2 || iplot4 == iplot3)
                     printf("ERROR: duplicate input number %d.\n",iplot4+1);
               }
            }
            if (iplot4 != -1) count++;
         }
         strcpy(command, "y\n");
         if (nInputs > count)
         {
            sprintf(pString,"Set other nominal values at mid point ? (y/n) ");
            getString(pString, command);
         }
         if (command[0] == 'y')
         {
            for (iInd1 = 0; iInd1 < nInputs; iInd1++)
            {
               if (iInd1 != iplot1 && iInd1 != iplot2 && iInd1 != iplot3 &&
                   iInd1 != iplot4)
                    inputSettings[iInd1] = 0.5*(iLowerB[iInd1]+iUpperB[iInd1]);
               else inputSettings[iInd1] = 1.0;
            }
         }
         else
         {
            for (iInd1 = 0; iInd1 < nInputs; iInd1++)
            {
               if (iInd1 != iplot1 && iInd1 != iplot2 && iInd1 != iplot3 &&
                   iInd1 != iplot4)
               {
                  inputSettings[iInd1] = iLowerB[iInd1] - 1.0;
                  sprintf(pString,
                          "Enter nominal value for input %d (%e - %e): ", 
                          iInd1+1, iLowerB[iInd1], iUpperB[iInd1]);
                  while (inputSettings[iInd1] < iLowerB[iInd1] ||
                         inputSettings[iInd1] > iUpperB[iInd1])
                     inputSettings[iInd1] = getDouble(pString);
               }
               else inputSettings[iInd1] = 1.0;
            }
         }
         if      (iplot2 == -1) nPtsPerDim = 1024;
         else if (iplot3 == -1) nPtsPerDim = 128;
         else if (iplot4 == -1) nPtsPerDim = 24;
         else                   nPtsPerDim = 10;
         printf("This command works with the following response surfaces:\n");
         printf("1. Linear    regression\n");
         printf("2. Quadratic regression\n");
         printf("3. cubic     regression\n");
         printf("4. quartic   regression\n");
         printf("5. GP1\n");
         printf("6. GP2\n");
         printf("7. MarsBagg\n");
         printf("8. Tree GP\n");
         printf("9. Kriging\n");
         sprintf(pString, "Enter your choice: (1, 2, ..., 9) ");
         faType = getInt(1, 9, pString);
         if      (faType == 1) faType = PSUADE_RS_REGR1;
         else if (faType == 2) faType = PSUADE_RS_REGR2;
         else if (faType == 3) faType = PSUADE_RS_REGR3;
         else if (faType == 4) faType = PSUADE_RS_REGR4;
         else if (faType == 5) faType = PSUADE_RS_GP1;
         else if (faType == 6) faType = PSUADE_RS_GP2;
         else if (faType == 7) faType = PSUADE_RS_MARSB;
         else if (faType == 8) faType = PSUADE_RS_TGP;
         else if (faType == 9) faType = PSUADE_RS_KR;
         faFlag = 1;
         faPtr = genFA(faType, nInputs, iOne, nSamples);
         if (faPtr == NULL) {printf("ERROR detected.\n"); continue;}
         faPtr->setNPtsPerDim(nPtsPerDim);
         faPtr->setBounds(iLowerB, iUpperB);
         faPtr->setOutputLevel(outputLevel_);

         jplot = 0;
         if (nOutputs > 1)
         {
            sprintf(pString, "Enter the output number (1 - %d) : ",nOutputs);
            jplot = getInt(1, nOutputs, pString);
            jplot--;
         }
         faYIn = new double[nSamples];
         for (sInd = 0; sInd < nSamples; sInd++)
            faYIn[sInd] = sampleOutputs[sInd*nOutputs+jplot];

         faLeng = -999;
         if (iplot2 == -1)
            faPtr->gen1DGridData(sampleInputs,faYIn,iplot1,
                                 inputSettings, &faLeng, &faXOut,&faYOut);
         else if (iplot3 == -1)
            faPtr->gen2DGridData(sampleInputs,faYIn,iplot1,iplot2,
                                inputSettings, &faLeng, &faXOut,&faYOut);
         else if (iplot4 == -1)
            faPtr->gen3DGridData(sampleInputs,faYIn, iplot1, iplot2, iplot3, 
                                 inputSettings, &faLeng, &faXOut,&faYOut);
         else
            faPtr->gen4DGridData(sampleInputs,faYIn, iplot1, iplot2, 
                                 iplot3, iplot4, inputSettings, &faLeng, 
                                 &faXOut,&faYOut);

         tempW = new double[faLeng];
         tempX = new double[faLeng*nInputs];
         for (sInd = 0; sInd < faLeng; sInd++)
            for (jj = 0; jj < nInputs; jj++)
               tempX[sInd*nInputs+jj] = inputSettings[jj];
         for (sInd = 0; sInd < faLeng; sInd++)
         {
            tempX[sInd*nInputs+iplot1] = faXOut[sInd*count];
            if (iplot2 != -1)
               tempX[sInd*nInputs+iplot2] = faXOut[sInd*count+1];
            if (iplot3 != -1)
               tempX[sInd*nInputs+iplot3] = faXOut[sInd*count+2];
            if (iplot4 != -1)
               tempX[sInd*nInputs+iplot4] = faXOut[sInd*count+3];
         }
         faPtr->evaluatePointFuzzy(faLeng, tempX, faYOut, tempW);
         gamma = PSUADE_UNDEFINED;
         for (sInd = 0; sInd < faLeng; sInd++)
            if (tempW[sInd] < gamma) gamma = tempW[sInd];
         
         if (psPlotTool_ == 1)
         {
            fp = fopen("scilabrssd.sci", "w");
            if (fp == NULL)
            {
               printf("ERROR: cannot open file scilabrssd.sci.\n");
               continue;
            }
         }
         else
         {
            fp = fopen("matlabrssd.m", "w");
            if (fp == NULL)
            {
               printf("ERROR: cannot open file matlabrssd.m.\n");
               continue;
            }
         }
         fwritePlotCLF(fp);
         if (count == 1)
         {
            fprintf(fp, "A = [\n");
            for (sInd = 0; sInd < faLeng; sInd++)
               fprintf(fp, "%e\n", tempW[sInd]);
            fprintf(fp, "];\n");
            fprintf(fp, "X = [\n");
            for (sInd = 0; sInd < faLeng; sInd++)
               fprintf(fp, "%e\n", faXOut[sInd]);
            fprintf(fp, "];\n");
            if (psPlotTool_ == 1)
            {
               fprintf(fp, "plot(X,A);");
               fprintf(fp, "a = gca();\n");
               fprintf(fp, "a.children.children.thickness = 4;\n");
               fprintf(fp, "set(gca(),\"auto_clear\",\"off\")\n");
            }
            else 
            {
               fprintf(fp, "plot(X,A,'lineWidth',4)\n");
               fprintf(fp, "hold on\n");
            }
            fwritePlotAxes(fp);
            fwritePlotXLabel(fp, inputNames[iplot1]);
            fwritePlotYLabel(fp, outputNames[jplot]);
            sprintf(winput, "Std. Dev. Plot for %s", outputNames[jplot]);
            fwritePlotTitle(fp, winput);
         }
         else if (count == 2)
         {
            if (psPlotTool_ == 0) fprintf(fp, "twoPlots = 1;\n");
            fprintf(fp, "A = [\n");
            for (sInd = 0; sInd < faLeng; sInd++)
               fprintf(fp, "%e\n", tempW[sInd]);
            fprintf(fp, "];\n");
            if (psPlotTool_ == 1) 
               fprintf(fp, "A = matrix(A,%d,%d);\n", nPtsPerDim, nPtsPerDim);
            else
               fprintf(fp, "A = reshape(A,%d,%d);\n", nPtsPerDim, nPtsPerDim);
            fprintf(fp, "X = [\n");
            for (sInd = 0; sInd < faLeng; sInd+=nPtsPerDim)
               fprintf(fp, "%e\n", faXOut[sInd*2]);
            fprintf(fp, "];\n");
            fprintf(fp, "Y = [\n");
            for (sInd = 0; sInd < nPtsPerDim; sInd++)
               fprintf(fp, "%e\n", faXOut[sInd*2+1]);
            fprintf(fp, "];\n");
            if (psPlotTool_ == 1)
            {
               fprintf(fp, "mesh(X,Y,A)\n");
               fprintf(fp, "h = get(\"hdl\");\n");
               fprintf(fp, "h.color_flag=1;\n");
               fprintf(fp, "h.color_mode=-2;\n");
               fprintf(fp, "bmin = min(min(A)); bmax = max(max(A));\n");
               fprintf(fp, "xset(\"colormap\",jetcolormap(64));\n");
               fprintf(fp, "colorbar(bmin,bmax);\n");
               fwritePlotAxes(fp);
               fwritePlotXLabel(fp, inputNames[iplot1]);
               fwritePlotYLabel(fp, inputNames[iplot2]);
               fwritePlotZLabel(fp, outputNames[jplot]);
               sprintf(winput, "Std. Dev. Plot for %s", outputNames[jplot]);
               fwritePlotTitle(fp, winput);
               fprintf(fp, "scf(2);\n");
               fprintf(fp, "a=gca();\n");
               fprintf(fp, "a.data_bounds=[%e,%e;%e,%e];\n",
                       iLowerB[iplot1], iLowerB[iplot2],
                       iUpperB[iplot1], iUpperB[iplot2]);
               fprintf(fp, "a.axes_visible=\"on\";\n");
               fprintf(fp, "B = A;\n");
               fprintf(fp, "nX = length(X);\n");
               fprintf(fp, "nY = length(Y);\n");
               fprintf(fp, "for ii = 1 : nX\n");
               fprintf(fp, "for jj = 1 : nY\n");
               fprintf(fp, "B(ii,jj) = A(nX-ii+1,jj);\n");
               fprintf(fp, "end;\n");
               fprintf(fp, "end;\n");
               fprintf(fp, "Matplot1((B-bmin)/(bmax-bmin)*64,[%e,%e,%e,%e])\n",
                       iLowerB[iplot1],iLowerB[iplot2],
                       iUpperB[iplot1], iUpperB[iplot2]);
               fprintf(fp, "xset(\"colormap\",jetcolormap(64));\n");
               fprintf(fp, "colorbar(bmin,bmax);\n");
               fprintf(fp, "a.thickness = 2;\n");
               fprintf(fp, "a.font_size = 3;\n");
               fprintf(fp, "a.font_style = 4;\n");
               fprintf(fp, "a.box = \"on\";\n");
               fprintf(fp, "a.grid = [1 1];\n");
               fwritePlotAxes(fp);
               fwritePlotXLabel(fp, inputNames[iplot1]);
               fwritePlotYLabel(fp, inputNames[iplot2]);
               sprintf(winput, "Std. Dev. Plot for %s", outputNames[jplot]);
               fwritePlotTitle(fp, winput);
            }
            else
            { 
               fprintf(fp, "if twoPlots == 1\n");
               fprintf(fp, "subplot(1,2,1), surf(X,Y,A)\n");
               fwritePlotAxes(fp);
               fwritePlotXLabel(fp, inputNames[iplot1]);
               fwritePlotYLabel(fp, inputNames[iplot2]);
               fwritePlotZLabel(fp, outputNames[jplot]);
               sprintf(winput, "Std. Dev. Plot for %s", outputNames[jplot]);
               fwritePlotTitle(fp, winput);
               fprintf(fp, "colorbar\n");
               fprintf(fp, "subplot(1,2,2)\n");
               fprintf(fp, "end\n");
               fprintf(fp, "contourf(X,Y,A)\n");
               fwritePlotAxes(fp);
               fwritePlotXLabel(fp, inputNames[iplot1]);
               fwritePlotYLabel(fp, inputNames[iplot2]);
               sprintf(winput, "Std. Dev. Plot for %s", outputNames[jplot]);
               fwritePlotTitle(fp, winput);
               fprintf(fp, "colorbar\n");
               fprintf(fp, "colormap(jet)\n");
            }
         }
         else if (count == 3)
         {
            fprintf(fp,"xlo = %e; \n", iLowerB[iplot2]);
            fprintf(fp,"xhi = %e; \n", iUpperB[iplot2]);
            fprintf(fp,"ylo = %e; \n", iLowerB[iplot1]);
            fprintf(fp,"yhi = %e; \n", iUpperB[iplot1]);
            fprintf(fp,"zlo = %e; \n", iLowerB[iplot3]);
            fprintf(fp,"zhi = %e; \n", iUpperB[iplot3]);
            fprintf(fp,"X=zeros(%d,%d,%d);\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);
            fprintf(fp,"Y=zeros(%d,%d,%d);\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);
            fprintf(fp,"Z=zeros(%d,%d,%d);\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);
            fprintf(fp,"V=zeros(%d,%d,%d);\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);
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
            GYmax = - PSUADE_UNDEFINED;
            GYmin =   PSUADE_UNDEFINED;
            for (jj = 0; jj < nPtsPerDim; jj++)
            {
               fprintf(fp, "V(:,:,%d) = [\n", jj + 1);
               for (sInd = 0; sInd < nPtsPerDim; sInd++)
               {
                  for (ii = 0; ii < nPtsPerDim; ii++)
                  {
                     ind = sInd*nPtsPerDim*nPtsPerDim+ii*nPtsPerDim+jj;
                     fprintf(fp, "%e ", tempW[ind]);
                     if (tempW[ind] > GYmax) GYmax = tempW[ind];
                     if (tempW[ind] < GYmin) GYmin = tempW[ind];
                  }
                  fprintf(fp, "\n");
               }
               fprintf(fp, "];\n");
            }
            fprintf(fp, "xt = [%e:%e:%e];\n", iLowerB[iplot2],
                    (iUpperB[iplot2]-iLowerB[iplot2])*0.01, iUpperB[iplot2]);
            fprintf(fp, "yt = [%e:%e:%e];\n", iLowerB[iplot1],
                    (iUpperB[iplot1]-iLowerB[iplot1])*0.01, iUpperB[iplot1]);
            fprintf(fp, "zt = [%e:%e:%e];\n", iLowerB[iplot3],
                    (iUpperB[iplot3]-iLowerB[iplot3])*0.01, iUpperB[iplot3]);
            fprintf(fp, "clf\n");
            fprintf(fp, "isoval = %e;\n", gamma);
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
            fprintf(fp, "daspect([xhi-xlo, yhi-ylo, zhi-zlo])\n");
            fprintf(fp, "colormap('default'); colorbar\n");
            fprintf(fp, "%%axis tight\n");
            fprintf(fp, "view(3) \n");
            fprintf(fp, "set(gcf,'Renderer','zbuffer')\n");
            fprintf(fp, "box on\n");
            fprintf(fp, "grid on\n");
            fprintf(fp, "lighting phong\n");
            fwritePlotAxes(fp);
            fwritePlotXLabel(fp, inputNames[iplot2]);
            fwritePlotYLabel(fp, inputNames[iplot1]);
            fwritePlotZLabel(fp, inputNames[iplot3]);
            sprintf(winput, "Std. Dev. Plot for %s", outputNames[jplot]);
            fwritePlotTitle(fp, winput);
            fprintf(fp,"cin = input('generate slices ? (y or n) ','s');\n");
            fprintf(fp,"if (cin == 'y')\n");
            fprintf(fp,"xin = input('axis to slide through? (x,y,z) ','s');\n");
            fprintf(fp,"N = 101;\n");
            fprintf(fp,"for i = 1 : N\n");
            fprintf(fp,"   display(['displaying ' int2str(i) ' of 101'])\n");
            fprintf(fp,"   if (xin == 'y')\n");
            fprintf(fp,"      h = slice(X,Y,Z,V,xt(i),[],[]);\n");
            fprintf(fp,"   elseif (xin == 'x')\n");
            fprintf(fp,"      h = slice(X,Y,Z,V,[],yt(i),[]);\n");
            fprintf(fp,"   elseif (xin == 'z')\n");
            fprintf(fp,"      h = slice(X,Y,Z,V,[],[],zt(i));\n");
            fprintf(fp,"   end\n");
            fprintf(fp,"   axis([%11.4e %11.4e %11.4e %11.4e %11.4e %11.4e ",
                    iLowerB[iplot2], iUpperB[iplot2], iLowerB[iplot1],
                    iUpperB[iplot1], iLowerB[iplot3], iUpperB[iplot3]);
            fprintf(fp, "%11.4e %11.4e])\n",
                    GYmin-0.1*(GYmax-GYmin),GYmax+0.1*(GYmax-GYmin));
            fwritePlotAxes(fp);
            fwritePlotXLabel(fp, inputNames[iplot2]);
            fwritePlotYLabel(fp, inputNames[iplot1]);
            fwritePlotZLabel(fp, inputNames[iplot3]);
            sprintf(winput, "Std. Dev. Slice Plot for %s", outputNames[jplot]);
            fwritePlotTitle(fp, winput);
            fprintf(fp, "   view(3)\n");
            fprintf(fp, "   colorbar\n");
            fprintf(fp, "   pause(1)\n");
            fprintf(fp, "   if (i < 101)\n");
            fprintf(fp, "      clf\n");
            fprintf(fp, "   end\n");
            fprintf(fp, "end\n");
            fprintf(fp, "end\n");
         }
         else if (count == 4)
         {
            fprintf(fp,"xlo = %e; \n", iLowerB[iplot2]);
            fprintf(fp,"xhi = %e; \n", iUpperB[iplot2]);
            fprintf(fp,"ylo = %e; \n", iLowerB[iplot1]);
            fprintf(fp,"yhi = %e; \n", iUpperB[iplot1]);
            fprintf(fp,"zlo = %e; \n", iLowerB[iplot3]);
            fprintf(fp,"zhi = %e; \n", iUpperB[iplot3]);
            fprintf(fp,"X=zeros(%d,%d,%d);\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);
            fprintf(fp,"Y=zeros(%d,%d,%d);\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);
            fprintf(fp,"Z=zeros(%d,%d,%d);\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);
            fprintf(fp,"V=zeros(%d,%d,%d);\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);
            for (ii = 0; ii < nPtsPerDim; ii++)
            {
               fprintf(fp, "Y(:,:,%d) = [\n", ii + 1);
               for (sInd = 0; sInd < nPtsPerDim; sInd++)
               {
                  for (jj = 0; jj < nPtsPerDim; jj++)
                  {
                     ind = (sInd*nPtsPerDim*nPtsPerDim+jj*nPtsPerDim+ii)*
                           nPtsPerDim;
                     fprintf(fp, "%e ", faXOut[ind*4]);
                  }
                  fprintf(fp, "\n");
               }
               fprintf(fp, "];\n");
               fprintf(fp, "X(:,:,%d) = [\n", ii + 1);
               for (sInd = 0; sInd < nPtsPerDim; sInd++)
               {
                  for (jj = 0; jj < nPtsPerDim; jj++)
                  {
                     ind = (sInd*nPtsPerDim*nPtsPerDim+jj*nPtsPerDim+ii)*
                           nPtsPerDim;
                     fprintf(fp, "%e ", faXOut[ind*4+1]);
                  }
                  fprintf(fp, "\n");
               }
               fprintf(fp, "];\n");
               fprintf(fp, "Z(:,:,%d) = [\n", ii + 1);
               for (sInd = 0; sInd < nPtsPerDim; sInd++)
               {
                  for (jj = 0; jj < nPtsPerDim; jj++)
                  {
                     ind = (sInd*nPtsPerDim*nPtsPerDim+jj*nPtsPerDim+ii)*
                           nPtsPerDim;
                     fprintf(fp, "%e ", faXOut[ind*4+2]);
                  }
                  fprintf(fp, "\n");
               }
               fprintf(fp, "];\n");
            }
            fprintf(fp, "xt = [%e:%e:%e];\n", iLowerB[iplot2],
                    (iUpperB[iplot2]-iLowerB[iplot2])*0.05, iUpperB[iplot2]);
            fprintf(fp, "yt = [%e:%e:%e];\n", iLowerB[iplot1],
                    (iUpperB[iplot1]-iLowerB[iplot1])*0.05, iUpperB[iplot1]);
            fprintf(fp, "zt = [%e:%e:%e];\n", iLowerB[iplot3],
                    (iUpperB[iplot3]-iLowerB[iplot3])*0.05, iUpperB[iplot3]);
            for (ll = 0; ll < nPtsPerDim; ll++)
            {
               for (ii = 0; ii < nPtsPerDim; ii++)
               {
                  fprintf(fp, "V(:,:,%d) = [\n", ii + 1);
                  for (sInd = 0; sInd < nPtsPerDim; sInd++)
                  {
                     for (jj = 0; jj < nPtsPerDim; jj++)
                     {
                        ind=((sInd*nPtsPerDim+jj)*nPtsPerDim+ii)*nPtsPerDim+ll;
                        fprintf(fp, "%e ", tempW[ind]);
                     }
                     fprintf(fp, "\n");
                  }
                  fprintf(fp, "];\n");
               }
               fprintf(fp, "disp('Frame %d of %d')\n", ll+1, nPtsPerDim);
               fprintf(fp, "clf\n");
               fprintf(fp, "isoval = %e;\n", gamma);
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
               fprintf(fp, "daspect([xhi-xlo, yhi-ylo, zhi-zlo])\n");
               fprintf(fp, "colormap('default'); colorbar\n");
               fprintf(fp, "%%axis tight\n");
               fprintf(fp, "view(3) \n");
               fprintf(fp, "set(gcf,'Renderer','zbuffer')\n");
               fprintf(fp, "lighting phong\n");
               fwritePlotAxes(fp);
               fwritePlotXLabel(fp, inputNames[iplot2]);
               fwritePlotYLabel(fp, inputNames[iplot1]);
               fwritePlotZLabel(fp, inputNames[iplot3]);
               fprintf(fp, "title('3D Std Dev Isosurface Plot at %s=%e',",
                       inputNames[iplot4],faXOut[ll*4+3]);
               fprintf(fp, "'FontWeight','bold','FontSize',12)\n");
               fprintf(fp, "pause(1)\n");
            }
         }
         fclose(fp);
         if (psPlotTool_ == 1)
              printf("\nscilabrssd.sci is now available.\n");
         else printf("\nmatlabrssd.m is now available.\n");
         delete [] faXOut;
         delete [] faYOut;
         delete [] faYIn;
         delete [] tempW;
         delete [] tempX;
         delete faPtr;
         faPtr = NULL;
         delete [] inputSettings;
      }

      // uncertainty analysis of standard deviations from RS fit 
      // ----------------------------------------------------------------
      else if (!strcmp(command, "rssd_ua") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rssd_ua: generate pdf for std. deviations from RS fit\n");
            printf("syntax: rssd_ua (no argument needed).\n");
            continue;
         }
         printf("ERROR: no data to analyze (load data first).\n");
      }
      else if (!strcmp(command, "rssd_ua") && nInputs > 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rssd_ua: generate pdf for std. deviations from RS fit\n");
            printf("syntax: rssd_ua (no argument needed).\n");
            continue;
         }
         printf("This command works with the following response surfaces:\n");
         printf("1. Linear    regression\n");
         printf("2. Quadratic regression\n");
         printf("3. cubic     regression\n");
         printf("4. quartic   regression\n");
         printf("5. GP1\n");
         printf("6. GP2\n");
         printf("7. MarsBagg\n");
         printf("8. Tree GP\n");
         printf("9. Kriging\n");
         sprintf(pString, "Enter your choice: (1, 2, ..., 9) ");
         faType = getInt(1, 9, pString);
         if      (faType == 1) faType = PSUADE_RS_REGR1;
         else if (faType == 2) faType = PSUADE_RS_REGR2;
         else if (faType == 3) faType = PSUADE_RS_REGR3;
         else if (faType == 4) faType = PSUADE_RS_REGR4;
         else if (faType == 5) faType = PSUADE_RS_GP1;
         else if (faType == 6) faType = PSUADE_RS_GP2;
         else if (faType == 7) faType = PSUADE_RS_MARSB;
         else if (faType == 8) faType = PSUADE_RS_TGP;
         else if (faType == 9) faType = PSUADE_RS_KR;

         jplot = 0;
         if (nOutputs > 1)
         {
            sprintf(pString, "Enter the output number (1 - %d) : ",nOutputs);
            jplot = getInt(1, nOutputs, pString);
            jplot--;
         }

         printf("rssd_ua: setting up function approximator\n");
         faPtr = genFA(faType, nInputs, iOne, nSamples);
         if (faPtr == NULL) {printf("ERROR detected in RS.\n"); continue;}
         faYIn = new double[nSamples];
         for (sInd = 0; sInd < nSamples; sInd++)
            faYIn[sInd] = sampleOutputs[sInd*nOutputs+jplot];
         faPtr->setBounds(iLowerB, iUpperB);
         faPtr->setOutputLevel(outputLevel_);
         faLeng = -999;
         faPtr->genNDGridData(sampleInputs,faYIn,&faLeng,NULL,NULL);

         printf("rssd_ua: creating a large sample for constructing PDF\n");
         sampPtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_LPTAU);
         sampPtr->setPrintLevel(0);
         sampPtr->setInputBounds(nInputs, iLowerB, iUpperB);
         sampPtr->setOutputParams(1);
         count = 100000;
         sampPtr->setSamplingParams(count, -1, 1);
         sampPtr->initialize(0);
         tempX  = new double[count*nInputs];
         tempY  = new double[count];
         tempW  = new double[count];
         states = new int[count];
         sampPtr->getSamples(count, nInputs, 1, tempX, tempY, states);
         faPtr->evaluatePointFuzzy(count, tempX, tempY, tempW);

         analysisMethod = PSUADE_ANA_MOMENT;
         anaManager = new AnalysisManager();
         anaManager->setup(analysisMethod, 0);
         vecUpper.load(nInputs, iUpperB);
         vecLower.load(nInputs, iLowerB);
         vecIn.load(count*nInputs, tempX);
         vecOut.load(count, tempW);
         printf("rssd_ua: analyzing the distribution of std deviation\n");
         anaManager->analyze(analysisMethod, count, vecLower, vecUpper,
                             vecIn, vecOut, 0);
         flag = 1;
         for (ii = 0; ii < count; ii++)
         {
            if (tempY[ii] == 0.0) flag = 0; 
            else                  tempW[ii] /= tempY[ii];
         }
         if (flag == 1)
         { 
            sprintf(pString,"analyze std dev with normalized output (y or n)? ");
            getString(pString, winput);
            if (winput[0] == 'y' ) 
            {
               vecOut.load(count, tempW);
               anaManager->analyze(analysisMethod, count, vecLower, vecUpper,
                                   vecIn, vecOut, 0);
            }
         }
         delete anaManager;
         delete sampPtr;
         delete faPtr;
         delete [] tempX;
         delete [] tempY;
         delete [] tempW;
         delete [] states;
         delete [] faYIn;
         sampPtr = NULL;
         tempX  = NULL;
         tempY  = NULL;
         tempW  = NULL;
         states = NULL;
         faYIn = NULL;
         faPtr = NULL;
      }

      // generate intersection surfaces for multiple outputs for 
      // display with matlab
      // ----------------------------------------------------------------
      else if (!strcmp(command, "rsi2") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rsi2: generate intersection surfaces for >1 outputs.\n");
            printf("syntax: rsi2 (no argument needed).\n");
            continue;
         }
         printf("ERROR: no data to analyze (load data first).\n");
      }
      else if (!strcmp(command, "rsi2") && nInputs < 2)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rsi2: generate intersection surfaces for >1 outputs.\n");
            printf("syntax: rsi2 (no argument needed).\n");
            continue;
         }
         printf("ERROR: rsi2 requires 2 or more inputs.\n");
      }
      else if (!strcmp(command, "rsi2") && nOutputs < 2)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rsi2: generate intersection surfaces for >1 outputs.\n");
            printf("syntax: rsi2 (no argument needed).\n");
            continue;
         }
         printf("ERROR: rsi2 requires 2 or more outputs.\n");
      }
      else if (!strcmp(command, "rsi2") && nInputs >= 2)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rsi2: generate intersection surfaces for >1 outputs.\n");
            printf("syntax: rsi2 (no argument needed).\n");
            continue;
         }
         if (psPlotTool_ == 1)
         {
            printf("INFO: rsi2 is currently not available for scilab.\n");
            continue;
         }
         nPtsPerDim = 32;
         sprintf(pString, "Grid resolution ? (32 - 256) ");
         nPtsPerDim = getInt(32, 256, pString);
         faFlag = 1;
         faPtr = genFAInteractive(psuadeIO_, faFlag);
         if (faPtr == NULL) {printf("ERROR detected.\n"); continue;}
         faPtr->setNPtsPerDim(nPtsPerDim);
         faPtr->setBounds(iLowerB, iUpperB);
         faPtr->setOutputLevel(outputLevel_);

         inputSettings = new double[nInputs];
         iplot1 = iplot2 = -1;
         sprintf(pString, "Enter the input for x axis (1 - %d) : ", nInputs);
         iplot1 = getInt(1, nInputs, pString);
         iplot1--;
         iplot2 = iplot1;
         while (iplot1 == iplot2)
         {
            sprintf(pString,"Enter the input for y axis (1 - %d), not %d : ",
                   nInputs, iplot1+1);
            iplot2 = getInt(1, nInputs, pString);
            iplot2--;
            if (iplot1 == iplot2)
               printf("ERROR: duplicate input number %d.\n", iplot2+1);
         }
         sprintf(pString,"Set other nominal values automatically ? (y or n) ");
         getString(pString, winput);
         if (winput[0] == 'y')
         {
            for (iInd1 = 0; iInd1 < nInputs; iInd1++)
            {
               if (iInd1 != iplot1 && iInd1 != iplot2)
                    inputSettings[iInd1] = 0.5*(iLowerB[iInd1]+iUpperB[iInd1]);
               else inputSettings[iInd1] = 1.0;
            }
         }
         else
         {
            for (iInd1 = 0; iInd1 < nInputs; iInd1++)
            {
               if (iInd1 != iplot1 && iInd1 != iplot2)
               {
                  
                  inputSettings[iInd1] = iLowerB[iInd1] - 1.0;
                  while (inputSettings[iInd1] < iLowerB[iInd1] ||
                         inputSettings[iInd1] > iUpperB[iInd1])
                  {
                     sprintf(pString,
                             "Enter nominal value for input %d (%e - %e):",
                             iInd1+1, iLowerB[iInd1], iUpperB[iInd1]);
                     inputSettings[iInd1] = getDouble(pString);
                  }
               }
               else inputSettings[iInd1] = 1.0;
            }
         }

         rsiNOutputs = 2;
         if (nOutputs > 2)
         {
            sprintf(pString,"How many outputs to use ? (2 - %d) ",nOutputs);
            rsiNOutputs = getInt(2, nOutputs, pString);
         }

         rsiSet = new int[rsiNOutputs];
         if (rsiNOutputs == nOutputs)
         {
            for (ii = 0; ii < rsiNOutputs; ii++) rsiSet[ii] = ii;
         }
         else
         {
            for (ii = 0; ii < rsiNOutputs; ii++)
            {
               sprintf(pString,"Enter the %d-th output index (1 - %d) : ",
                       ii+1, nOutputs);
               rsiSet[ii] = getInt(1, nOutputs, pString);
               rsiSet[ii]--;
            }
         }
          
         faYIn = new double[nSamples];
         rsiMatrix = new int*[nPtsPerDim];
         for (ii = 0; ii < nPtsPerDim; ii++)
         {
            rsiMatrix[ii] = new int[nPtsPerDim];
            for (jj = 0; jj < nPtsPerDim; jj++)
               rsiMatrix[ii][jj] = rsiNOutputs;
         }
         fp = fopen("matlabrsi2.m", "w");
         if (fp == NULL)
         {
            printf("ERROR: cannot open file matlabrsi2.m.\n");
            continue;
         }
         fprintf(fp, "twoPlots = 1;\n");
         fprintf(fp, "clf\n");

         for (ii = 0; ii < rsiNOutputs; ii++)
         {
            jplot = rsiSet[ii];
            for (sInd = 0; sInd < nSamples; sInd++)
               faYIn[sInd] = sampleOutputs[sInd*nOutputs+jplot];

            faPtr->gen2DGridData(sampleInputs,faYIn, iplot1, iplot2, 
                                 inputSettings, &faLeng, &faXOut,&faYOut);

            Ymin = faYOut[0];
            for (sInd = 1; sInd < faLeng; sInd++)
               if (faYOut[sInd] < Ymin) Ymin = faYOut[sInd];
            Ymax = faYOut[0];
            for (sInd = 1; sInd < faLeng; sInd++)
               if (faYOut[sInd] > Ymax) Ymax = faYOut[sInd];

            printf("Ymin and Ymax = %e %e\n", Ymin, Ymax);
            sprintf(pString,
                    "Enter the lower threshold for output %d (min = %16.8e) : ",
                    jplot+1, Ymin);
            threshL = getDouble(pString);
            sprintf(pString,
                    "Enter the upper threshold for output %d (max = %16.8e) : ",
                    jplot+1, Ymax);
            threshU = getDouble(pString);

            if (ii == 0)
            {
               fprintf(fp, "x = [\n");
               for (sInd = 0; sInd < faLeng; sInd+=nPtsPerDim)
                  fprintf(fp, "%e\n", faXOut[sInd*2]);
               fprintf(fp, "];\n");
               fprintf(fp, "y = [\n");
               for (sInd = 0; sInd < nPtsPerDim; sInd++)
                  fprintf(fp, "%e\n", faXOut[sInd*2+1]);
               fprintf(fp, "];\n");
            }
            if (ii < 4)
            {
               fprintf(fp, "A%d = [\n", ii+1);
               for (sInd = 0; sInd < faLeng; sInd++)
                  fprintf(fp, "%e\n", faYOut[sInd]);
               fprintf(fp, "];\n");
               fprintf(fp, "A%d = reshape(A%d,%d,%d);\n",ii+1,ii+1,
                       nPtsPerDim,nPtsPerDim);
               fprintf(fp, "yLo = %e;\n", threshL);
               fprintf(fp, "yHi = %e;\n", threshU);
               fprintf(fp, "nA  = size(A%d,1);\n", ii+1);
               fprintf(fp, "[ia,ja,aa] = find(A%d<yLo);\n", ii+1);
               fprintf(fp, "for ii = 1 : length(ia)\n");
               fprintf(fp, "   A%d(ia(ii),ja(ii)) = NaN;\n", ii+1); 
               fprintf(fp, "end;\n");
               fprintf(fp, "n1 = length(ia);\n");
               fprintf(fp, "[ia,ja,aa] = find(A%d>yHi);\n", ii+1);
               fprintf(fp, "for ii = 1 : length(ia)\n");
               fprintf(fp, "   A%d(ia(ii),ja(ii)) = NaN;\n", ii+1); 
               fprintf(fp, "end;\n");
               fprintf(fp, "n2 = length(ia);\n");
               fprintf(fp, "if (n1 + n2 == nA*nA)\n");
               fprintf(fp, "   A%d(1,1) = 0;\n",ii+1);
               fprintf(fp, "   A%d(%d,%d) = 1;\n",ii+1,nPtsPerDim,
                       nPtsPerDim);
               fprintf(fp, "end;\n");
               if (ii == 0) fprintf(fp, "subplot(2,3,1)\n");
               if (ii == 1) fprintf(fp, "subplot(2,3,2)\n");
               if (ii == 2) fprintf(fp, "subplot(2,3,3)\n");
               if (ii == 3) fprintf(fp, "subplot(2,3,4)\n");
               fprintf(fp, "contourf(x,y,A%d)\n", ii+1);
               fprintf(fp, "axis([%e %e %e %e])\n",iLowerB[iplot1],
                       iUpperB[iplot1],iLowerB[iplot2],iUpperB[iplot2]);
               fwritePlotAxes(fp);
               fprintf(fp, "xlabel('%s','FontSize',12,'FontWeight','bold')\n",
                       inputNames[iplot1]);
               fprintf(fp, "ylabel('%s','Fontsize',12,'FontWeight','bold')\n",
                       inputNames[iplot2]);
               fprintf(fp, "title('Plot for %s',",outputNames[jplot]);
               fprintf(fp, "'FontWeight','bold','FontSize',12)\n");
               fprintf(fp, "colorbar\n");
            }

            for (sInd = 0; sInd < faLeng; sInd++)
            {
               ind  = sInd % nPtsPerDim;
               ind2 = sInd / nPtsPerDim;
               if (faYOut[sInd] < threshL) rsiMatrix[ind][ind2]--;
               if (faYOut[sInd] > threshU) rsiMatrix[ind][ind2]--;
            }
            delete [] faXOut;
            delete [] faYOut;
         }


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
         fprintf(fp, "if twoPlots == 1\n");
         fprintf(fp, "subplot(2,3,5), mesh(y,x,A)\n");
         fprintf(fp, "axis([%e %e %e %e])\n",iLowerB[iplot1],
                 iUpperB[iplot1],iLowerB[iplot2],iUpperB[iplot2]);
         fwritePlotAxes(fp);
         fprintf(fp, "xlabel('%s','FontSize',12,'FontWeight','bold')\n",
                 inputNames[iplot1]);
         fprintf(fp, "ylabel('%s','Fontsize',12,'FontWeight','bold')\n",
                 inputNames[iplot2]);
         fprintf(fp, "title('Intersection Mesh','FontWeight',");
         fprintf(fp, "'bold','FontSize',12)\n");
         fprintf(fp, "colorbar\n");
         fprintf(fp, "colormap(cool)\n");
         fprintf(fp, "end\n");
         fprintf(fp, "subplot(2,3,6), contourf(x,y,A)\n");
         fprintf(fp, "axis([%e %e %e %e])\n",iLowerB[iplot1],
                 iUpperB[iplot1],iLowerB[iplot2],iUpperB[iplot2]);
         fwritePlotAxes(fp);
         fprintf(fp, "xlabel('%s','FontSize',12,'FontWeight','bold')\n",
                 inputNames[iplot1]);
         fprintf(fp, "ylabel('%s','Fontsize',12,'FontWeight','bold')\n",
                 inputNames[iplot2]);
         fprintf(fp, "title('Intersection Contour','FontWeight',");
         fprintf(fp, "'bold','FontSize',12)\n");
         fprintf(fp, "colorbar\n");
         fprintf(fp, "colormap(cool)\n");
         fclose(fp);
         printf("matlabrsi2.m is now available for plotting.\n");

         delete [] inputSettings;
         delete [] rsiSet;
         delete [] faYIn;
         delete faPtr;
         faPtr = NULL;
         for (ii = 0; ii < nPtsPerDim; ii++) delete [] rsiMatrix[ii];
         delete [] rsiMatrix;
      }

      // generate 3D response surface and write the grid data to file
      // for display with matlab (movie mode)
      // ----------------------------------------------------------------
      else if (!strcmp(command, "rsi3m") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rsi3m: generate intersection surfaces for >1 outputs\n");
            printf("syntax: rsi3m (no argument needed).\n");
            continue;
         }
         printf("ERROR: no data to analyze (load data first).\n");
      }
      else if (!strcmp(command, "rsi3m") && nInputs < 3)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rsi3m: generate intersection surfaces for >1 outputs\n");
            printf("syntax: rsi3m (no argument needed).\n");
            continue;
         }
         printf("ERROR: rsi3m requires 3 or more inputs.\n");
      }
      else if (!strcmp(command, "rsi3m") && nOutputs < 2)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rsi3m: generate intersection surfaces for >1 outputs\n");
            printf("syntax: rsi3m (no argument needed).\n");
            continue;
         }
         printf("ERROR: rsi3m requires 2 or more outputs.\n");
      }
      else if (!strcmp(command, "rsi3m") && nInputs >= 3)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rsi3m: generate intersection surfaces for >1 outputs\n");
            printf("syntax: rsi3m (no argument needed).\n");
            continue;
         }
         if (psPlotTool_ == 1)
         {
            printf("INFO: rsi3m is currently not available for scilab.\n");
            continue;
         }
         nPtsPerDim = 24;
         sprintf(pString, "Grid resolution ? (16 - 32) ");
         nPtsPerDim = getInt(16, 32, pString);
         faFlag = 1;
         faPtr = genFAInteractive(psuadeIO_, faFlag);
         if (faPtr == NULL) {printf("ERROR detected.\n"); continue;}
         faPtr->setNPtsPerDim(nPtsPerDim);
         faPtr->setBounds(iLowerB, iUpperB);
         faPtr->setOutputLevel(outputLevel_);

         inputSettings = new double[nInputs];
         iplot1 = iplot2 = iplot3 = -1;
         sprintf(pString, "Enter the input for x axis (1 - %d) : ", nInputs);
         iplot1 = getInt(1, nInputs, pString);
         iplot1--;
         iplot2 = iplot1;
         while (iplot1 == iplot2)
         {
            sprintf(pString,"Enter the input for y axis (1 - %d), not %d : ",
                    nInputs, iplot1+1);
            iplot2 = getInt(1, nInputs, pString);
            iplot2--;
            if (iplot1 == iplot2)
               printf("ERROR: duplicate input number %d.\n",iplot2+1);
         }
         if (nInputs == 3) iplot3 = 3 - iplot1 - iplot2;
         while (iplot3 < 0 || iplot3 == iplot1 || iplot3 == iplot2)
         {
            sprintf(pString,
                    "Enter the input for t axis (1 - %d), not %d nor %d: ",
                    nInputs, iplot1+1, iplot2+1);
            iplot3 = getInt(1, nInputs, pString);
            iplot3--;
            if (iplot3 == iplot1 || iplot3 == iplot2)
               printf("ERROR: duplicate input number %d.\n",iplot3+1);
         }
         sprintf(pString,"Set other nominal values automatically ? (y or n) ");
         getString(pString, winput);
         if (winput[0] == 'y')
         {
            for (iInd1 = 0; iInd1 < nInputs; iInd1++)
            {
               if (iInd1 != iplot1 && iInd1 != iplot2 && iInd1 != iplot3)
                    inputSettings[iInd1] = 0.5*(iLowerB[iInd1]+iUpperB[iInd1]);
               else inputSettings[iInd1] = 1.0;
            }
         }
         else
         {
            for (iInd1 = 0; iInd1 < nInputs; iInd1++)
            {
               if (iInd1 != iplot1 && iInd1 != iplot2 && iInd1 != iplot3)
               {
                  inputSettings[iInd1] = iLowerB[iInd1] - 1.0;
                  sprintf(pString,
                          "Enter nominal value for input %d (%e - %e): ", 
                          iInd1+1, iLowerB[iInd1], iUpperB[iInd1]);
                  while (inputSettings[iInd1] < iLowerB[iInd1] ||
                         inputSettings[iInd1] > iUpperB[iInd1])
                     inputSettings[iInd1] = getDouble(pString);
               }
               else inputSettings[iInd1] = 1.0;
            }
         }

         rsiNOutputs = 2;
         if (nOutputs > 2)
         {
            sprintf(pString,"How many outputs to use ? (2 - %d) ",nOutputs);
            rsiNOutputs = getInt(2, nOutputs, pString);
         }

         rsiSet = new int[rsiNOutputs];
         if (rsiNOutputs == nOutputs)
         {
            for (ii = 0; ii < rsiNOutputs; ii++) rsiSet[ii] = ii;
         }
         else
         {
            for (ii = 0; ii < rsiNOutputs; ii++)
            {
               sprintf(pString,"Enter the %d-th output index (1 - %d) : ",
                       ii+1, nOutputs);
               rsiSet[ii] = getInt(1, nOutputs, pString);
               rsiSet[ii]--;
            }
         }
         faYIn = new double[nSamples];

         printf("Please wait while generating the RS data \n");
         fp = fopen("matlabrsi3m.m", "w");
         if (fp == NULL)
         {
            printf("ERROR: cannot open file matlabrsi3m.m.\n");
            continue;
         }
         fprintf(fp, "hold off\n");
         fprintf(fp, "clf\n");
         fprintf(fp, "disp(\'Please wait while loading.\')\n");
         fprintf(fp, "pause(1)\n");
         for (jj = 0; jj < nPtsPerDim; jj++)
            fprintf(fp, "M%d = %e * ones(%d);\n", jj+1, 1.0*rsiNOutputs, 
                    nPtsPerDim);
         for (ii = 0; ii < rsiNOutputs; ii++)
         {
            jplot = rsiSet[ii];
            for (sInd = 0; sInd < nSamples; sInd++)
               faYIn[sInd] = sampleOutputs[sInd*nOutputs+jplot];

            faPtr->gen3DGridData(sampleInputs,faYIn, iplot1, iplot2, iplot3, 
                                 inputSettings, &faLeng, &faXOut,&faYOut);
            GYmax = -1.0e35;
            GYmin =  1.0e35;
            GYmin = faYOut[0];
            for (sInd = 1; sInd < faLeng; sInd++)
               if (faYOut[sInd] < GYmin) GYmin = faYOut[sInd];
            GYmax = faYOut[0];
            for (sInd = 1; sInd < faLeng; sInd++)
               if (faYOut[sInd] > GYmax) GYmax = faYOut[sInd];

            for (jj = 0; jj < nPtsPerDim; jj++)
            {
               printf(".");
               fflush(stdout);

               fprintf(fp, "A%d_%d = [\n", ii+1, jj+1);
               for (sInd = 0; sInd < faLeng; sInd+=nPtsPerDim)
                  fprintf(fp, "%e\n", faYOut[sInd+jj]);
               fprintf(fp, "];\n");
               fprintf(fp, "A%d_%d = reshape(A%d_%d,%d,%d);\n", ii+1, jj+1,
                       ii+1, jj+1, nPtsPerDim, nPtsPerDim);

               if (ii == 0 && jj == 0)
               {
                  fprintf(fp, "x = [\n");
                  for (sInd = 0; sInd < faLeng; sInd+=nPtsPerDim*nPtsPerDim)
                     fprintf(fp, "%e\n", faXOut[sInd*3]);
                  fprintf(fp, "];\n");
                  fprintf(fp, "y = [\n");
                  for (sInd = 0; sInd < nPtsPerDim*nPtsPerDim; sInd+=nPtsPerDim)
                     fprintf(fp, "%e\n", faXOut[sInd*3+1]);
                  fprintf(fp, "];\n");
               }
            }
            delete [] faXOut;
            delete [] faYOut;
            printf("\nOutput %d : Ymin and Ymax found = %e %e.\n", jplot+1,
                   GYmin, GYmax);
            sprintf(pString,"Enter the lower threshold (min = %e) : ", GYmin);
            threshL = getDouble(pString);
            sprintf(pString,"Enter the upper threshold (max = %e) : ", GYmax);
            threshU = getDouble(pString);

            for (jj = 0; jj < nPtsPerDim; jj++)
            {
               fprintf(fp, "B%d_%d = A%d_%d;\n",ii+1,jj+1,ii+1,jj+1);
               fprintf(fp, "nA  = size(A%d_%d,1);\n", ii+1, jj+1);
               fprintf(fp, "n1 = 0;\n");
               fprintf(fp, "n2 = 0;\n");
               if (threshL > GYmin)
               { 
                  fprintf(fp, "yLo = %e;\n", threshL);
                  fprintf(fp, "[ia,ja,aa] = find(A%d_%d<yLo);\n",ii+1,jj+1);
                  fprintf(fp, "for ii = 1 : length(ia)\n");
                  fprintf(fp, "   B%d_%d(ia(ii),ja(ii))=NaN;\n",ii+1,jj+1);
                  fprintf(fp, "   M%d(ia(ii),ja(ii))=M%d(ia(ii),ja(ii))-1;\n", 
                          jj+1,jj+1);
                  fprintf(fp, "end;\n");
                  fprintf(fp, "n1 = length(ia);\n");
               }
               if (threshU < GYmax)
               { 
                  fprintf(fp, "yHi = %e;\n", threshU);
                  fprintf(fp, "[ia,ja,aa] = find(A%d_%d>yHi);\n",ii+1,jj+1);
                  fprintf(fp, "for ii = 1 : length(ia)\n");
                  fprintf(fp, "   B%d_%d(ia(ii),ja(ii))=NaN;\n",ii+1,jj+1);
                  fprintf(fp, "   M%d(ia(ii),ja(ii))=M%d(ia(ii),ja(ii))-1;\n", 
                          jj+1,jj+1);
                  fprintf(fp, "end;\n");
                  fprintf(fp, "n1 = length(ia);\n");
               }
               fprintf(fp, "if (n1+n2 == nA*nA)\n");
               fprintf(fp, "   B%d_%d(1,1)=0;\n",ii+1,jj+1);
               fprintf(fp, "   B%d_%d(%d,%d)=1;\n",ii+1,jj+1,
                       nPtsPerDim,nPtsPerDim);
               fprintf(fp, "end;\n");
            }
         }
         for (jj  = 0; jj < nPtsPerDim; jj++)
         {
            fprintf(fp, "[ia,ja,aa] = find(M%d == 0);\n", jj+1);
            fprintf(fp, "nM  = size(M%d,1);\n", jj+1);
            fprintf(fp, "for ii = 1 : length(ia)\n");
            fprintf(fp, "   M%d(ia(ii),ja(ii)) = NaN;\n", jj+1);
            fprintf(fp, "end;\n");
            fprintf(fp, "if (length(ia) == nM*nM)\n");
            fprintf(fp, "   M%d(1,1) = 0;\n", jj+1);
            fprintf(fp, "   M%d(nM,nM) = %e;\n", jj+1, 1.0*rsiNOutputs);
            fprintf(fp, "end;\n");
            fprintf(fp, "Mmax = max(max(M%d));\n", jj+1);
            fprintf(fp, "if (Mmax ~= %d)\n", rsiNOutputs);
            fprintf(fp, "   M%d(%d,%d) = %d;\n", jj+1, nPtsPerDim,
                    nPtsPerDim, rsiNOutputs);
            fprintf(fp, "end;\n");
            fprintf(fp, "Mmin = min(min(M%d));\n", jj+1);
            fprintf(fp, "if (Mmin ~= 0)\n");
            fprintf(fp, "   M%d(1,1) = 0;\n", jj+1);
            fprintf(fp, "end;\n");
         }

         for (jj = 0; jj < nPtsPerDim; jj++)
         {
            inputSettings[iplot3] = (iUpperB[iplot3] - iLowerB[iplot3]) *
                                    jj / (nPtsPerDim - 1.0) + iLowerB[iplot3];
            fprintf(fp, "disp(\'Plotting frame %d of %d\')\n",jj+1,nPtsPerDim);
            fprintf(fp, "subplot(2,3,1), contourf(x,y,B1_%d)\n", jj+1);
            fwritePlotAxes(fp);
            fprintf(fp,"title(\'Contour Plot for %s\',",outputNames[rsiSet[0]]);
            fprintf(fp, "'FontSize',12,'FontWeight','bold')\n"); 
            fprintf(fp, "xlabel('%s','FontSize',12,'FontWeight','bold')\n",
                    inputNames[iplot1]);
            fprintf(fp, "ylabel('%s','FontSize',12,'FontWeight','bold');\n",
                    inputNames[iplot2]);
            fprintf(fp, "subplot(2,3,2), contourf(x,y,B2_%d)\n", jj+1);
            fwritePlotAxes(fp);
            fprintf(fp,"title(\'Contour Plot for %s\',",outputNames[rsiSet[1]]);
            fprintf(fp, "'FontSize',12,'FontWeight','bold')\n"); 
            fprintf(fp, "xlabel('%s','FontSize',12,'FontWeight','bold')\n",
                    inputNames[iplot1]);
            fprintf(fp, "ylabel('%s','FontSize',12,'FontWeight','bold');\n",
                    inputNames[iplot2]);
            if (rsiNOutputs > 2)
            {
               fprintf(fp, "subplot(2,3,3), contourf(x,y,B3_%d)\n", jj+1);
               fwritePlotAxes(fp);
               fprintf(fp,"title(\'Contour Plot for %s\',",
                       outputNames[rsiSet[2]]);
               fprintf(fp, "'FontSize',12,'FontWeight','bold')\n"); 
               fprintf(fp, "xlabel('%s','FontSize',12,'FontWeight','bold')\n",
                       inputNames[iplot1]);
               fprintf(fp, "ylabel('%s','FontSize',12,'FontWeight','bold');\n",
                       inputNames[iplot2]);
            }
            if (rsiNOutputs > 3)
            {
               fprintf(fp, "subplot(2,3,4), contourf(x,y,B4_%d)\n", jj+1);
               fwritePlotAxes(fp);
               fprintf(fp, "title(\'Contour Plot for ");
               fprintf(fp, "%s\','FontSize',12,'FontWeight','bold')\n", 
                       outputNames[rsiSet[3]]);
               fprintf(fp, "xlabel('%s','FontSize',12,'FontWeight','bold')\n",
                       inputNames[iplot1]);
               fprintf(fp, "ylabel('%s','FontSize',12,'FontWeight','bold');\n",
                       inputNames[iplot2]);
            }
            if (rsiNOutputs > 4)
            {
               fprintf(fp, "subplot(2,3,5), contourf(x,y,B5_%d)\n", jj+1);
               fwritePlotAxes(fp);
               fprintf(fp, "title(\'Contour Plot for ");
               fprintf(fp, "%s\','FontSize',12,'FontWeight','bold')\n", 
                       outputNames[rsiSet[4]]);
               fprintf(fp, "xlabel('%s','FontSize',12,'FontWeight','bold')\n",
                       inputNames[iplot1]);
               fprintf(fp, "ylabel('%s','FontSize',12,'FontWeight','bold');\n",
                       inputNames[iplot2]);
            }
            fprintf(fp, "subplot(2,3,6), contourf(x,y,M%d)\n",jj+1);
            fwritePlotAxes(fp);
            fwritePlotXLabel(fp, inputNames[iplot1]);
            fwritePlotYLabel(fp, inputNames[iplot2]);
            fprintf(fp, "title('Intersection: Input %s = %11.4e',",
                    inputNames[iplot3], inputSettings[iplot3]);
            fprintf(fp, "'FontSize',12,'FontWeight','bold')\n");
            fprintf(fp, "colorbar\n");
            fprintf(fp, "colormap(jet)\n");
            fprintf(fp,"pause(1)\n");
         }
         fclose(fp);
         printf("matlabrsi3m.m is now available for response surface and ");
         printf("contour plots\n");

         delete [] rsiSet;
         delete [] faYIn;
         delete faPtr;
         faPtr = NULL;
         delete [] inputSettings;
      }

      // generate 3D response surface and write the grid data to file
      // for display with matlab
      // ----------------------------------------------------------------
      else if (!strcmp(command, "rsi3") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rsi3: generate intersection surfaces for >1 outputs\n");
            printf("syntax: rsi3 (no argument needed).\n");
            continue;
         }
         printf("ERROR: no data to analyze (load data first).\n");
      }
      else if (!strcmp(command, "rsi3") && nInputs < 3)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rsi3: generate intersection surfaces for >1 outputs\n");
            printf("syntax: rsi3 (no argument needed).\n");
            continue;
         }
         printf("ERROR: rsi3 requires 3 or more inputs.\n");
      }
      else if (!strcmp(command, "rsi3") && nOutputs < 2)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rsi3: generate intersection surfaces for >1 outputs\n");
            printf("syntax: rsi3 (no argument needed).\n");
            continue;
         }
         printf("ERROR: rsi3 requires 2 or more outputs.\n");
      }
      else if (!strcmp(command, "rsi3") && nInputs >= 3)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rsi3: generate intersection surfaces for >1 outputs\n");
            printf("syntax: rsi3 (no argument needed).\n");
            continue;
         }
         if (psPlotTool_ == 1)
         {
            printf("INFO: rsi3 is currently not available for scilab.\n");
            continue;
         }
         nPtsPerDim = 24;
         sprintf(pString, "Grid resolution ? (16 - 32) ");
         nPtsPerDim = getInt(16, 32, pString);
         faFlag = 1;
         faPtr = genFAInteractive(psuadeIO_, faFlag);
         if (faPtr == NULL) {printf("ERROR detected.\n"); continue;}
         faPtr->setNPtsPerDim(nPtsPerDim);
         faPtr->setBounds(iLowerB, iUpperB);
         faPtr->setOutputLevel(outputLevel_);

         inputSettings = new double[nInputs];
         iplot1 = iplot2 = iplot3 = -1;
         sprintf(pString, "Enter the input for x axis (1 - %d) : ", nInputs);
         iplot1 = getInt(1, nInputs, pString);
         iplot1--;
         iplot2 = iplot1;
         while (iplot1 == iplot2)
         {
            sprintf(pString,"Enter the input for y axis (1 - %d), not %d : ",
                    nInputs, iplot1+1);
            iplot2 = getInt(1, nInputs, pString);
            iplot2--;
            if (iplot1 == iplot2)
               printf("ERROR: duplicate input number %d.\n",iplot2+1);
         }
         if (nInputs == 3) iplot3 = 3 - iplot1 - iplot2;
         while (iplot3 < 0 || iplot3 == iplot1 || iplot3 == iplot2)
         {
            sprintf(pString,
                    "Enter the input for z axis (1 - %d), not %d nor %d: ",
                    nInputs, iplot1+1, iplot2+1);
            iplot3 = getInt(1, nInputs, pString);
            iplot3--;
            if (iplot3 == iplot1 || iplot3 == iplot2)
               printf("ERROR: duplicate input number %d.\n",iplot3+1);
         }
         sprintf(pString,"Set other nominal values automatically ? (y or n) ");
         getString(pString, winput);
         if (winput[0] == 'y')
         {
            for (iInd1 = 0; iInd1 < nInputs; iInd1++)
            {
               if (iInd1 != iplot1 && iInd1 != iplot2 && iInd1 != iplot3)
                    inputSettings[iInd1] = 0.5*(iLowerB[iInd1]+iUpperB[iInd1]);
               else inputSettings[iInd1] = 1.0;
            }
         }
         else
         {
            for (iInd1 = 0; iInd1 < nInputs; iInd1++)
            {
               if (iInd1 != iplot1 && iInd1 != iplot2 && iInd1 != iplot3)
               {
                  inputSettings[iInd1] = iLowerB[iInd1] - 1.0;
                  sprintf(pString,
                          "Enter nominal value for input %d (%e - %e): ", 
                          iInd1+1, iLowerB[iInd1], iUpperB[iInd1]);
                  while (inputSettings[iInd1] < iLowerB[iInd1] ||
                         inputSettings[iInd1] > iUpperB[iInd1])
                     inputSettings[iInd1] = getDouble(pString);
               }
               else inputSettings[iInd1] = 1.0;
            }
         }

         rsiNOutputs = 1;
         if (nOutputs > 1)
         {
            sprintf(pString,"How many outputs to use ? (2 - %d) ",nOutputs);
            rsiNOutputs = getInt(2, nOutputs, pString);
         }
         threshLs = new double[rsiNOutputs];
         threshUs = new double[rsiNOutputs];

         rsiSet = new int[rsiNOutputs];
         if (rsiNOutputs == nOutputs)
         {
            for (ii = 0; ii < rsiNOutputs; ii++) rsiSet[ii] = ii;
         }
         else
         {
            for (ii = 0; ii < rsiNOutputs; ii++)
            {
               sprintf(pString,"Enter the %d-th output index (1 - %d) : ",
                       ii+1, nOutputs);
               rsiSet[ii] = getInt(1, nOutputs, pString);
               rsiSet[ii]--;
            }
         }
         faYIn = new double[nSamples];

         printf("Please wait while generating the RS data \n");
         fp = fopen("matlabrsi3.m", "w");
         if (fp == NULL)
         {
            printf("ERROR: cannot open file matlabrsi3.m.\n");
            continue;
         }
         fwritePlotCLF(fp);
         fprintf(fp, "xlo = %e; \n", iLowerB[iplot2]);
         fprintf(fp, "xhi = %e; \n", iUpperB[iplot2]);
         fprintf(fp, "ylo = %e; \n", iLowerB[iplot1]);
         fprintf(fp, "yhi = %e; \n", iUpperB[iplot1]);
         fprintf(fp, "zlo = %e; \n", iLowerB[iplot3]);
         fprintf(fp, "zhi = %e; \n", iUpperB[iplot3]);
         fprintf(fp, "X=zeros(%d,%d,%d);\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);
         fprintf(fp, "Y=zeros(%d,%d,%d);\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);
         fprintf(fp, "Z=zeros(%d,%d,%d);\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);
         fprintf(fp, "V=zeros(%d,%d,%d);\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);

         for (ii = 0; ii < rsiNOutputs; ii++)
         {
            jplot = rsiSet[ii];
            for (sInd = 0; sInd < nSamples; sInd++)
               faYIn[sInd] = sampleOutputs[sInd*nOutputs+jplot];

            faPtr->gen3DGridData(sampleInputs,faYIn, iplot1, iplot2, iplot3, 
                                 inputSettings, &faLeng, &faXOut,&faYOut);
            GYmax = -1.0e35;
            GYmin =  1.0e35;
            GYmin = faYOut[0];
            for (sInd = 1; sInd < faLeng; sInd++)
               if (faYOut[sInd] < GYmin) GYmin = faYOut[sInd];
            GYmax = faYOut[0];
            for (sInd = 1; sInd < faLeng; sInd++)
               if (faYOut[sInd] > GYmax) GYmax = faYOut[sInd];

            printf("\nOutput %d : Ymin and Ymax found = %e %e.\n", jplot+1,
                   GYmin, GYmax);
            sprintf(pString,"Enter the lower threshold (min = %e) : ", GYmin);
            threshLs[ii] = getDouble(pString);
            sprintf(pString,"Enter the upper threshold (max = %e) : ", GYmax);
            threshUs[ii] = getDouble(pString);
            if (ii == 0) gamma = threshLs[ii];  
            else         gamma = (gamma < threshLs[ii]) ? gamma : threshLs[ii];

            if (ii == 0)
            {
               for (jj = 0; jj < nPtsPerDim; jj++)
               {
                  fprintf(fp, "Y(:,:,%d) = [\n", jj + 1);
                  for (sInd = 0; sInd < nPtsPerDim; sInd++)
                  {
                     for (kk = 0; kk < nPtsPerDim; kk++)
                     {
                        ind = sInd*nPtsPerDim*nPtsPerDim+kk*nPtsPerDim+jj;
                        fprintf(fp, "%e ", faXOut[ind*3]);
                     }
                     fprintf(fp, "\n");
                  }
                  fprintf(fp, "];\n");
                  fprintf(fp, "X(:,:,%d) = [\n", jj + 1);
                  for (sInd = 0; sInd < nPtsPerDim; sInd++)
                  {
                     for (kk = 0; kk < nPtsPerDim; kk++)
                     {
                        ind = sInd*nPtsPerDim*nPtsPerDim+kk*nPtsPerDim+jj;
                        fprintf(fp, "%e ", faXOut[ind*3+1]);
                     }
                     fprintf(fp, "\n");
                  }
                  fprintf(fp, "];\n");
                  fprintf(fp, "Z(:,:,%d) = [\n", jj + 1);
                  for (sInd = 0; sInd < nPtsPerDim; sInd++)
                  {
                     for (kk = 0; kk < nPtsPerDim; kk++)
                     {
                        ind = sInd*nPtsPerDim*nPtsPerDim+kk*nPtsPerDim+jj;
                        fprintf(fp, "%e ", faXOut[ind*3+2]);
                     }
                     fprintf(fp, "\n");
                  }
                  fprintf(fp, "];\n");
               }
            }
            delete [] faXOut;
            delete [] faYOut;
         }

         for (ii = 0; ii < rsiNOutputs; ii++)
         {
            jplot = rsiSet[ii];
            for (sInd = 0; sInd < nSamples; sInd++)
               faYIn[sInd] = sampleOutputs[sInd*nOutputs+jplot];

            faPtr->gen3DGridData(sampleInputs,faYIn, iplot1, iplot2, iplot3, 
                                 inputSettings, &faLeng, &faXOut,&faYOut);
            for (jj = 0; jj < nPtsPerDim; jj++)
            {
               fprintf(fp, "V%d(:,:,%d) = [\n", ii+1, jj+1);
               for (sInd = 0; sInd < nPtsPerDim; sInd++)
               {
                  for (kk = 0; kk < nPtsPerDim; kk++)
                  {
                     ind = sInd*nPtsPerDim*nPtsPerDim+kk*nPtsPerDim+jj;
                     if (faYOut[ind] < threshLs[ii])
                     {
                        fprintf(fp, "%e ", gamma);
                        count++;
                     }
                     else if (faYOut[ind] > threshUs[ii])
                     {
                        fprintf(fp, "%e ", gamma);
                        count++;
                     }
                     else fprintf(fp, "%e ", faYOut[ind]);
                  }
                  fprintf(fp, "\n");
               }
               fprintf(fp, "];\n");
            }
            delete [] faXOut;
            delete [] faYOut;
            if (ii == 0) fprintf(fp, "V = V%d;\n", ii+1);
            else         fprintf(fp, "V = min(V, V%d);\n", ii+1);
         }

         fprintf(fp, "xt = [%e:%e:%e];\n", iLowerB[iplot2],
                 (iUpperB[iplot2]-iLowerB[iplot2])*0.01, iUpperB[iplot2]);
         fprintf(fp, "yt = [%e:%e:%e];\n", iLowerB[iplot1],
                 (iUpperB[iplot1]-iLowerB[iplot1])*0.01, iUpperB[iplot1]);
         fprintf(fp, "zt = [%e:%e:%e];\n", iLowerB[iplot3],
                 (iUpperB[iplot3]-iLowerB[iplot3])*0.01, iUpperB[iplot3]);
         fprintf(fp, "isoval = %e;\n", gamma);
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
         fprintf(fp, "daspect([%e,%e,%e])\n",iUpperB[iplot2]-iLowerB[iplot2],
                 iUpperB[iplot1]-iLowerB[iplot1],
                 iUpperB[iplot3]-iLowerB[iplot3]);
         fprintf(fp, "   xlabel('%s','FontSize',12,'FontWeight','bold')\n",
                 inputNames[iplot2]);
         fprintf(fp, "   ylabel('%s','Fontsize',12,'FontWeight','bold')\n",
                 inputNames[iplot1]);
         fprintf(fp, "   zlabel('%s','Fontsize',12,'FontWeight','bold')\n",
                 inputNames[iplot3]);
         fwritePlotAxes(fp);
         fprintf(fp, "%%colormap('default'); colorbar\n");
         fprintf(fp, "%%axis tight\n");
         fprintf(fp, "view(3) \n");
         fprintf(fp, "set(gcf,'Renderer','zbuffer')\n");
         fprintf(fp, "lighting phong\n");
         fprintf(fp, "cin = input('generate slices ? (y or n) ','s');\n");
         fprintf(fp, "if (cin == 'y')\n");
         fprintf(fp, "xin = input('axis to slide through ? (x,y,z) ','s');\n");
         fprintf(fp, "for i = 1 : 101\n");
         fprintf(fp, "   if (xin == 'y')\n");
         fprintf(fp, "      h = contourslice(X,Y,Z,V,xt(i),[],[],101);\n");
         fprintf(fp, "   elseif (xin == 'x')\n");
         fprintf(fp, "      h = contourslice(X,Y,Z,V,[],yt(i),[],101);\n");
         fprintf(fp, "   elseif (xin == 'z')\n");
         fprintf(fp, "      h = contourslice(X,Y,Z,V,[],[],zt(i),101);\n");
         fprintf(fp, "   end\n");
         fprintf(fp, "   axis([%11.4e %11.4e %11.4e %11.4e %11.4e %11.4e ",
                 iLowerB[iplot2], iUpperB[iplot2], iLowerB[iplot1],
                 iUpperB[iplot1], iLowerB[iplot3], iUpperB[iplot3]);
         fprintf(fp, "%11.4e %11.4e])\n",
                 threshL-0.2*(threshU-threshL),threshU+0.2*(threshU-threshL));
         fwritePlotAxes(fp);
         fprintf(fp, "   xlabel('%s','FontSize',12,'FontWeight','bold')\n",
                 inputNames[iplot2]);
         fprintf(fp, "   ylabel('%s','Fontsize',12,'FontWeight','bold')\n",
                 inputNames[iplot1]);
         fprintf(fp, "   zlabel('%s','Fontsize',12,'FontWeight','bold')\n",
                 inputNames[iplot3]);
         fwritePlotAxes(fp);
         fprintf(fp, "colormap('default'); colorbar\n");
         fprintf(fp, "%%axis tight\n");
         fprintf(fp, "view(3) \n");
         fprintf(fp, "set(gcf,'Renderer','zbuffer')\n");
         fprintf(fp, "lighting phong\n");
         fprintf(fp, "end\n");
         fprintf(fp, "end\n");
         fclose(fp);
         printf("matlabrsi3.m is now available for response surface and ");
         printf("contour plots\n");

         delete [] rsiSet;
         delete [] faYIn;
         delete faPtr;
         faPtr = NULL;
         delete [] inputSettings;
      }

      // generate intersection surfaces for multiple outputs for 
      // display with matlab (using raw instead of RS data)
      // ----------------------------------------------------------------
      else if (!strcmp(command, "rawi2") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rawi2: similar to rsi2 except with sample (not RS) data\n");
            printf("syntax: rawi2 (no argument needed).\n");
            continue;
         }
         printf("ERROR: no data to analyze (load data first).\n");
      }
      else if (!strcmp(command, "rawi2") && nInputs < 2)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rawi2: similar to rsi2 except with sample (not RS) data\n");
            printf("syntax: rawi2 (no argument needed).\n");
            continue;
         }
         printf("ERROR: rawi2 requires 2 inputs.\n");
      }
      else if (!strcmp(command, "rawi2") && nOutputs < 2)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rawi2: similar to rsi2 except with sample (not RS) data\n");
            printf("syntax: rawi2 (no argument needed).\n");
            continue;
         }
         printf("ERROR: rawi2 requires 2 or more outputs.\n");
      }
      else if (!strcmp(command, "rawi2") && nInputs >= 2)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rawi2: similar to rsi2 except with sample (not RS) data\n");
            printf("syntax: rawi2 (no argument needed).\n");
            continue;
         }
         if (psPlotTool_ == 1)
         {
            printf("INFO: rawi2 is currently not available for scilab.\n");
            continue;
         }
         dtemp = pow(1.0*nSamples, 0.5) + 1.0e-12;
         nPtsPerDim = (int) dtemp;
         if (nPtsPerDim * nPtsPerDim != nSamples)
         {
            printf("rawi2 error: nSamples must be a square.\n");
         }
         else
         { 
            rsiNOutputs = 2;
            if (nOutputs > 2)
            {
               sprintf(pString,"How many outputs to use ? (2 - %d) ",nOutputs);
               rsiNOutputs = getInt(2, nOutputs, pString);
            }

            rsiSet = new int[rsiNOutputs];
            if (rsiNOutputs == nOutputs)
            {
               for (ii = 0; ii < rsiNOutputs; ii++) rsiSet[ii] = ii;
            }
            else
            {
               for (ii = 0; ii < rsiNOutputs; ii++)
               {
                  sprintf(pString,"Enter the %d-th output index (1 - %d) : ",
                          ii+1, nOutputs);
                  rsiSet[ii] = getInt(1, nOutputs, pString);
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
            fprintf(fp, "clf\n");
            faYOut = new double[nSamples];

            for (ii = 0; ii < rsiNOutputs; ii++)
            {
               jplot = rsiSet[ii];
               for (sInd = 0; sInd < nSamples; sInd++)
                  faYOut[sInd] = sampleOutputs[sInd*nOutputs+jplot];

               Ymin = faYOut[0];
               for (sInd = 1; sInd < nSamples; sInd++)
                  if (faYOut[sInd] < Ymin) Ymin = faYOut[sInd];
               Ymax = faYOut[0];
               for (sInd = 1; sInd < nSamples; sInd++)
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

               for (sInd = 0; sInd < nSamples; sInd++)
               {
                  ind  = sInd % nPtsPerDim;
                  ind2 = sInd / nPtsPerDim;
                  if (faYOut[sInd] < threshL) rsiMatrix[ind][ind2]--;
                  if (faYOut[sInd] > threshU) rsiMatrix[ind][ind2]--;
               }
            }
            delete [] faYOut;

            fprintf(fp, "X = [\n");
            for (sInd = 0; sInd < nSamples; sInd+=nPtsPerDim)
               fprintf(fp, "%e\n", sampleInputs[sInd*2]);
            fprintf(fp, "];\n");
            fprintf(fp, "Y = [\n");
            for (sInd = 0; sInd < nPtsPerDim; sInd++)
               fprintf(fp, "%e\n", sampleInputs[sInd*2+1]);
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
            fprintf(fp, "axis([%e %e %e %e])\n",iLowerB[0],
                    iUpperB[0],iLowerB[1],iUpperB[1]);
            fwritePlotAxes(fp);
            fwritePlotXLabel(fp, inputNames[0]);
            fwritePlotYLabel(fp, inputNames[1]);
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

      // generate 3D response surface and write the grid data to file
      // for display with matlab using raw data
      // ----------------------------------------------------------------
      else if (!strcmp(command, "rawi3") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rawi3: similar to rsi3 except with sample (not RS) data\n");
            printf("syntax: rawi3 (no argument needed).\n");
            continue;
         }
         printf("ERROR: no data to analyze (load data first).\n");
      }
      else if (!strcmp(command, "rawi3") && nInputs < 3)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rawi3: similar to rsi3 except with sample (not RS) data\n");
            printf("syntax: rawi3 (no argument needed).\n");
            continue;
         }
         printf("ERROR: rawi3 requires 3 or more inputs.\n");
      }
      else if (!strcmp(command, "rawi3") && nOutputs < 2)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rawi3: similar to rsi3 except with sample (not RS) data\n");
            printf("syntax: rawi3 (no argument needed).\n");
            continue;
         }
         printf("ERROR: rawi3 requires 2 or more outputs.\n");
      }
      else if (!strcmp(command, "rawi3"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rawi3: similar to rsi3 except with sample (not RS) data\n");
            printf("syntax: rawi3 (no argument needed).\n");
            continue;
         }
         if (psPlotTool_ == 1)
         {
            printf("INFO: rawi3 is currently not available for scilab.\n");
            continue;
         }
         dtemp = pow(1.0*nSamples, 0.333333) + 0.1;
         if (nPtsPerDim*nPtsPerDim*nPtsPerDim != nSamples)
         {
            printf("rawi3 error: nSamples must be an integer 3rd power.\n");
         }
         iplot1 = 0; iplot2 = 1; iplot3 = 2;

         rsiNOutputs = 2;
         if (nOutputs > 2)
         {
            sprintf(pString,"How many outputs to use ? (2 - %d) ",nOutputs);
            rsiNOutputs = getInt(2, nOutputs, pString);
         }

         rsiSet = new int[rsiNOutputs];
         if (rsiNOutputs == nOutputs)
         {
            for (ii = 0; ii < rsiNOutputs; ii++) rsiSet[ii] = ii;
         }
         else
         {
            for (ii = 0; ii < rsiNOutputs; ii++)
            {
               sprintf(pString,"Enter the %d-th output index (1 - %d) : ",
                       ii+1, nOutputs);
               rsiSet[ii] = getInt(1, nOutputs, pString);
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

         faXOut = sampleInputs;
         faYOut = new double[nSamples];
         for (ii = 0; ii < rsiNOutputs; ii++)
         {
            jplot = rsiSet[ii];
            for (sInd = 0; sInd < nSamples; sInd++)
               faYOut[sInd] = sampleOutputs[sInd*nOutputs+jplot];

            Ymax = faYOut[0];
            Ymin = faYOut[0];
            for (sInd = 1; sInd < nSamples; sInd++)
               if (faYOut[sInd] < Ymin) Ymin = faYOut[sInd];
            for (sInd = 1; sInd < nSamples; sInd++)
               if (faYOut[sInd] > Ymax) Ymax = faYOut[sInd];

            printf("\nOutput %d : Ymin and Ymax found = %e %e.\n", jplot,
                   Ymin, Ymax);
            sprintf(pString,"Enter the lower threshold (min = %e) : ", Ymin);
            threshL = getDouble(pString);
            sprintf(pString,"Enter the upper threshold (max = %e) : ", Ymax);
            threshU = getDouble(pString);

            for (sInd = 0; sInd < nSamples; sInd++)
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
         fprintf(fp, "clf\n");
         fprintf(fp, "xlo = %e; \n", iLowerB[1]);
         fprintf(fp, "xhi = %e; \n", iUpperB[1]);
         fprintf(fp, "ylo = %e; \n", iLowerB[0]);
         fprintf(fp, "yhi = %e; \n", iUpperB[0]);
         fprintf(fp, "zlo = %e; \n", iLowerB[2]);
         fprintf(fp, "zhi = %e; \n", iUpperB[2]);
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

         fprintf(fp, "xt = [%e:%e:%e];\n", iLowerB[1],
                 (iUpperB[1]-iLowerB[1])*0.01, iUpperB[1]);
         fprintf(fp, "yt = [%e:%e:%e];\n", iLowerB[iplot1],
                 (iUpperB[0]-iLowerB[0])*0.01, iUpperB[0]);
         fprintf(fp, "zt = [%e:%e:%e];\n", iLowerB[2],
                 (iUpperB[2]-iLowerB[2])*0.01, iUpperB[2]);
         fprintf(fp, "clf\n");
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
         fprintf(fp, "daspect([%e,%e,%e])\n",iUpperB[1]-iLowerB[1],
                 iUpperB[0]-iLowerB[0], iUpperB[2]-iLowerB[2]);
         fprintf(fp, "colormap('default'); colorbar\n");
         fprintf(fp, "%%axis tight\n");
         fprintf(fp, "view(3) \n");
         fprintf(fp, "set(gcf,'Renderer','zbuffer')\n");
         fprintf(fp, "box on\n");
         fprintf(fp, "grid on\n");
         fprintf(fp, "lighting phong\n");
         fwritePlotAxes(fp);
         fprintf(fp, "xlabel('%s','FontSize',12,'FontWeight','bold')\n",
                 inputNames[1]);
         fprintf(fp, "ylabel('%s','Fontsize',12,'FontWeight','bold')\n",
                 inputNames[0]);
         fprintf(fp, "zlabel('%s','Fontsize',12,'FontWeight','bold')\n",
                 inputNames[2]);
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

      // generate response surface of all 2-input pairs and write the
      // grid data to file for display with matlab
      // ----------------------------------------------------------------
      else if (!strcmp(command, "rspairs") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rspairs: generate RS of all 2-input pairs\n");
            printf("syntax: rspairs (no argument needed).\n");
            continue;
         }
         printf("ERROR: no data to analyze (load data first).\n");
      }
      else if (!strcmp(command, "rspairs") && nInputs < 2)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rspairs: generate RS of all 2-input pairs\n");
            printf("syntax: rspairs (no argument needed).\n");
            continue;
         }
         printf("ERROR: rspairs requires 2 or more inputs.\n");
      }
      else if (!strcmp(command, "rspairs") && nInputs >= 2)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rspairs: generate RS of all 2-input pairs\n");
            printf("syntax: rspairs (no argument needed).\n");
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
         faPtr->setBounds(iLowerB, iUpperB);
         faPtr->setOutputLevel(outputLevel_);

         nPlots = 0;
         sprintf(pString, "Enter the number of inputs to plot (1 - %d) : ",
                nInputs);
         nPlots = getInt(1, nInputs, pString);
         plotIndices = new int[nInputs];
         for (ii = 0; ii < nInputs; ii++) plotIndices[ii] = ii;
         if (nPlots < nInputs)
         {
            for (ii = 0; ii < nPlots; ii++)
            {
               sprintf(pString, "Enter the %d-th input (1 - %d) : ", ii+1,
                       nInputs);
               plotIndices[ii] = getInt(1, nInputs, pString);
               plotIndices[ii]--;
            }
         }
         inputSettings = new double[nInputs];
         if (nInputs > 2)
         {
            sprintf(pString,
                    "Set other nominal values at mid point ? (y or n) ");
            getString(pString, command);
         }
         if (command[0] == 'y')
         {
            for (iInd1 = 0; iInd1 < nInputs; iInd1++)
               inputSettings[iInd1] = 0.5*(iLowerB[iInd1]+iUpperB[iInd1]);
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
            fgets(lineIn,500,stdin); 
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
                  if (kk != nInputs)
                  {
                     printf("ERROR: input size does not match nInputs.\n");
                     fclose(fp);
                     delete [] inputSettings;
                     continue;
                  }
                  for (ii = 0; ii < nInputs; ii++)
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
                  if (ii != nInputs)
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
         fprintf(fp, "clf\n");

         jplot = 0;
         if (nOutputs > 1)
         {
            sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
            jplot = getInt(1, nOutputs, pString);
            jplot--;
         }
         Ymax = Ymin = sampleOutputs[0*nOutputs+jplot];
         for (sInd = 1; sInd < nSamples; sInd++)
         {
            if (sampleOutputs[sInd*nOutputs+jplot] > Ymax) 
               Ymax = sampleOutputs[sInd*nOutputs+jplot];
            if (sampleOutputs[sInd*nOutputs+jplot] < Ymin) 
               Ymin = sampleOutputs[sInd*nOutputs+jplot];
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

         faYIn = new double[nSamples];
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
                  for (sInd = 0; sInd < nSamples; sInd++)
                     faYIn[sInd] = sampleOutputs[sInd*nOutputs+jplot];
                  faPtr->gen1DGridData(sampleInputs,faYIn,iplot2,
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
                  fwritePlotXLabel(fp, inputNames[iplot1]);
                  fwritePlotYLabel(fp, outputNames[jplot]);
                  delete [] faXOut;
                  delete [] faYOut;
                  continue;
               }

               for (sInd = 0; sInd < nSamples; sInd++)
                  faYIn[sInd] = sampleOutputs[sInd*nOutputs+jplot];
               faPtr->gen2DGridData(sampleInputs,faYIn,iplot1,iplot2,
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
               fwritePlotXLabel(fp, inputNames[iplot1]);
               fwritePlotYLabel(fp, inputNames[iplot2]);
               fprintf(fp, "axis([%e %e %e %e])\n",iLowerB[iplot1],
                       iUpperB[iplot1],iLowerB[iplot2],iUpperB[iplot2]);
               fprintf(fp, "subplot(%d,%d,%d), ",nPlots,nPlots,jj*nPlots+ii+1);
               fprintf(fp, "surf(X,Y,B)\n");
               fwritePlotAxes(fp);
               fwritePlotXLabel(fp, inputNames[iplot1]);
               fwritePlotYLabel(fp, inputNames[iplot2]);
               delete [] faXOut;
               delete [] faYOut;
            }
         }
         if (diagMax - diagMin == 0) diagMax += 0.1;
         for (ii = 0; ii < nPlots; ii++)
         {
            iplot1 = plotIndices[ii];
            fprintf(fp, "subplot(%d,%d,%d), ",nPlots,nPlots,ii*nPlots+ii+1);
            fprintf(fp, "axis([%e %e %e %e])\n",iLowerB[iplot1],
                    iUpperB[iplot1],diagMin, diagMax);
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

      // generate intersection surfaces for multiple outputs for 
      // display with matlab
      // ----------------------------------------------------------------
      else if (!strcmp(command, "rsipairs") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rsipairs: generate RS intersections for all input pairs\n");
            printf("syntax: rsipairs (no argument needed).\n");
            continue;
         }
         printf("ERROR: no data to analyze (load data first).\n");
      }
      else if (!strcmp(command, "rsipairs") && nInputs < 2)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rsipairs: generate RS intersections for all input pairs\n");
            printf("syntax: rsipairs (no argument needed).\n");
            continue;
         }
         printf("ERROR: rsipairs requires 2 or more inputs.\n");
      }
      else if (!strcmp(command, "rsipairs") && nInputs >= 2)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rsipairs: generate RS intersections for all input pairs\n");
            printf("syntax: rsipairs (no argument needed).\n");
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
         faPtr->setBounds(iLowerB, iUpperB);
         faPtr->setOutputLevel(outputLevel_);

         nPlots = 0;
         sprintf(pString, "Enter the number of inputs to plot (1 - %d) : ",
                nInputs);
         nPlots = getInt(1, nInputs, pString);
         plotIndices = new int[nInputs];
         for (ii = 0; ii < nInputs; ii++) plotIndices[ii] = ii;
         if (nPlots < nInputs)
         {
            for (ii = 0; ii < nPlots; ii++)
            {
               sprintf(pString, "Enter the %d-th input (1 - %d) : ", ii+1,
                       nInputs);
               plotIndices[ii] = getInt(1, nInputs, pString);
               plotIndices[ii]--;
            }
         }
    
         printf("rsipairs: will use all outputs for constraining.\n");
         rsiNOutputs = nOutputs;
         rsiSet = new int[rsiNOutputs];
         for (ii = 0; ii < rsiNOutputs; ii++) rsiSet[ii] = ii;
         threshLs = new double[rsiNOutputs];
         threshUs = new double[rsiNOutputs];
         faYIn = new double[nSamples];
         for (ii = 0; ii < rsiNOutputs; ii++)
         {
            jplot = rsiSet[ii];
            Ymax = Ymin = sampleOutputs[0*nOutputs+jplot];
            for (sInd = 1; sInd < nSamples; sInd++)
            {
               if (sampleOutputs[sInd*nOutputs+jplot] > Ymax) 
                  Ymax = sampleOutputs[sInd*nOutputs+jplot];
               if (sampleOutputs[sInd*nOutputs+jplot] < Ymin) 
                  Ymin = sampleOutputs[sInd*nOutputs+jplot];
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

         faYIn = new double[nSamples];
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
         fprintf(fp, "clf\n");

         inputSettings = new double[nInputs];
         for (ii = 0; ii < nPlots; ii++)
         {
            iplot2 = plotIndices[ii];
            for (jj = 0; jj <= ii; jj++)
            {
               iplot1 = plotIndices[jj];
               for (iInd1 = 0; iInd1 < nInputs; iInd1++)
               {
                  if (iInd1 != iplot1 && iInd1 != iplot2)
                       inputSettings[iInd1]=0.5*(iLowerB[iInd1]+iUpperB[iInd1]);
                  else inputSettings[iInd1]=1.0;
               }
               for (iInd1 = 0; iInd1 < nPtsPerDim; iInd1++)
                  for (iInd2 = 0; iInd2 < nPtsPerDim; iInd2++)
                     rsiMatrix[iInd1][iInd2] = rsiNOutputs;
               for (kk = 0; kk < rsiNOutputs; kk++)
               {
                  jplot = rsiSet[kk];
                  for (sInd = 0; sInd < nSamples; sInd++)
                     faYIn[sInd] = sampleOutputs[sInd*nOutputs+jplot];

                  if (iplot1 == iplot2)
                  {
                     faPtr->gen1DGridData(sampleInputs,faYIn,iplot1,
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
                     faPtr->gen2DGridData(sampleInputs,faYIn, iplot1, iplot2, 
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
                  fprintf(fp, "axis([%e %e 0 %d])\n",iLowerB[iplot1],
                          iUpperB[iplot1],rsiNOutputs+1);
                  fwritePlotAxes(fp);
                  fwritePlotXLabel(fp, inputNames[iplot1]);
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
                  fprintf(fp, "axis([%e %e %e %e])\n",iLowerB[iplot1],
                          iUpperB[iplot1],iLowerB[iplot2],iUpperB[iplot2]);
                  fwritePlotAxes(fp);
                  fwritePlotXLabel(fp, inputNames[iplot1]);
                  fwritePlotYLabel(fp, inputNames[iplot2]);
                  fprintf(fp, "subplot(%d,%d,%d)\n",
                          nPlots,nPlots,jj*nPlots+ii+1);
                  fprintf(fp, "surf(X,Y,A)\n");
                  fprintf(fp, "axis([%e %e %e %e 0 %d])\n",iLowerB[iplot1],
                      iUpperB[iplot1],iLowerB[iplot2],iUpperB[iplot2],
                      rsiNOutputs);
                  fwritePlotAxes(fp);
                  fwritePlotXLabel(fp, inputNames[iplot1]);
                  fwritePlotYLabel(fp, inputNames[iplot2]);
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

      // check quality of the response surface 
      // ----------------------------------------------------------------
      else if (!strcmp(command, "rscheck") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rscheck: check the quality of the RS (training errors\n");
            printf("         and cross validation errors)\n");
            printf("syntax: rscheck (no argument needed)\n");
            continue;
         }
         printf("ERROR: no data to analyze (load data first).\n");
      }
      else if (!strcmp(command, "rscheck") && nInputs > 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rscheck: check the quality of the RS (training errors\n");
            printf("         and cross validation errors)\n");
            printf("syntax: rscheck (no argument needed)\n");
            continue;
         }
         if (psuadeIO_ == NULL) printf("ERROR: data not loaded yet.\n");
         else
         {
            faType = -1;
            sprintf(pString, "Enter your choice ? ");
            while (faType < 0 || faType > PSUADE_NUM_RS)
            {
               writeFAInfo();
               faType = getFAType(pString);
            }
            outputID = 0;
            if (nOutputs > 1)
            {
               sprintf(pString, "Enter output number (1 - %d) = ", nOutputs);
               outputID = getInt(1, nOutputs, pString);
               outputID--;
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
      }

      // response surface test on training set
      // ----------------------------------------------------------------
      else if (!strcmp(command, "rstest_ts") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rstest_ts: response surface test on the training set.\n");
            printf("syntax: rstest_ts (no argument needed)\n");
            continue;
         }
         printf("ERROR: no data to analyze (load data first).\n");
      }
      else if (!strcmp(command, "rstest_ts") && nInputs > 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rstest_ts: response surface test on the training set.\n");
            printf("syntax: rstest_ts (no argument needed)\n");
            continue;
         }
         if (psuadeIO_ == NULL) printf("ERROR: data not loaded yet.\n");
         else
         {
            faType = -1;
            sprintf(pString, "Enter your choice ? ");
            while (faType < 0 || faType > PSUADE_NUM_RS)
            {
               writeFAInfo();
               faType = getFAType(pString);
            }
            outputID = 0;
            if (nOutputs > 1)
            {
               sprintf(pString, "Enter output number (1 - %d) = ", nOutputs);
               outputID = getInt(1, nOutputs, pString);
               outputID--;
            }
            faPtr = genFA(faType, nInputs, iOne, nSamples);
            faPtr->setBounds(iLowerB, iUpperB);
            faPtr->setOutputLevel(outputLevel_);
            faLeng = -999;
            tempY = new double[nSamples];
            for (ss = 0; ss < nSamples; ss++)
               tempY[ss] = sampleOutputs[ss*nOutputs+outputID];
            status = faPtr->genNDGridData(sampleInputs,tempY,&faLeng,NULL,NULL);
            tempV = new double[nSamples];
            tempW = new double[nSamples];
            faPtr->evaluatePointFuzzy(nSamples, sampleInputs, tempV, tempW);
            if (psPlotTool_ == 1)
            {
               fp = fopen("RSTest_ts.sci", "w");
               fprintf(fp, "// col 1: simulation data, col 2: rs data\n");
               fprintf(fp, "// col 3: std dev\n");
            }
            else
            {
               fp = fopen("RSTest_ts.m", "w");
               fprintf(fp, "%% col 1: simulation data, col 2: rs data\n");
               fprintf(fp, "%% col 3: std dev\n");
            }
            fprintf(fp, "A = [\n");
            for (ss = 0; ss < nSamples; ss++)
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
               fprintf(fp, "bar(xk,nk/%d,1.0)\n",nSamples);
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
            for (ss = 0; ss < nSamples; ss++)
               ddata += (tempY[ss] - tempV[ss]);
            ddata = ddata / nSamples;
            printf("Training set avg error (unscaled) = %e\n", ddata);
            ddata = 0.0;
            for (ss = 0; ss < nSamples; ss++)
            {
               ddata += (tempY[ss] - tempV[ss]);
               if (tempY[ss] != 0.0)
               {
                  ddata /= PABS(tempY[ss]);
                  dtemp = tempY[ss];
               }
            }
            ddata = ddata / nSamples;
            printf("Training set avg error (  scaled) = %e (base=%e)\n",ddata,dtemp);
            ddata = 0.0;
            for (ss = 0; ss < nSamples; ss++)
               ddata += pow(tempY[ss] - tempV[ss], 2.0);
            ddata = sqrt(ddata / nSamples);
            printf("Training set rms error (unscaled) = %e\n", ddata);
            ddata = 0.0;
            for (ss = 0; ss < nSamples; ss++)
            {
               ddata += pow(tempY[ss] - tempV[ss], 2.0);
               if (tempY[ss] != 0.0)
               {
                  ddata /= pow(tempY[ss],2.0);
                  dtemp = tempY[ss];
               }
            }
            ddata = sqrt(ddata / nSamples);
            printf("Training set rms error (  scaled) = %e (base=%e)\n",ddata,dtemp);
            ddata = 0.0;
            for (ss = 0; ss < nSamples; ss++)
               if (PABS(tempV[ss]-tempY[ss]) > ddata)
                  ddata = PABS(tempY[ss] - tempV[ss]);
            printf("Training set max error (unscaled) = %e\n", ddata);
            ddata = 0.0;
            for (ss = 0; ss < nSamples; ss++)
            {
               if (tempY[ss] != 0 && PABS((tempV[ss]-tempY[ss])/tempY[ss])>ddata)
               {
                  dtemp = tempY[ss];
                  ddata = PABS((tempY[ss]-tempV[ss])/tempY[ss]);
               }
            }
            printf("Training set max error (  scaled) = %e (base=%e)\n",ddata,dtemp);
            ddata = 0.0;
            for (ss = 0; ss < nSamples; ss++)
               ddata += (tempY[ss] - tempV[ss]);
            ddata = ddata / nSamples;
            printf("Training set error mean = %e\n", ddata);
            dtemp = 0.0;
            for (ss = 0; ss < nSamples; ss++)
               dtemp += pow(tempY[ss] - tempV[ss] - ddata, 2.0);
            dtemp = sqrt(dtemp / (nSamples - 1.0));
            printf("Training set error std  = %e\n", dtemp);
            if (psPlotTool_ == 1)
                 printf("rstest_ts plot file in RSTest_ts.sci\n");
            else printf("rstest_ts plot file in RSTest_ts.m\n");
            delete faPtr;
            delete [] tempY;
            delete [] tempV;
            delete [] tempW;
            faPtr = NULL;
            tempV = NULL;
            tempY = NULL;
            tempW = NULL;
         }
      }

      // response surface test - cross validation
      // ----------------------------------------------------------------
      else if (!strcmp(command, "rstest_cv") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rstest_cv: cross validation response surface test.\n");
            printf("syntax: rstest_cv (no argument needed)\n");
            continue;
         }
         printf("ERROR: no data to analyze (load data first).\n");
      }
      else if (!strcmp(command, "rstest_cv") && nInputs > 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rstest_cv: cross validation response surface test.\n");
            printf("syntax: rstest_cv (no argument needed)\n");
            continue;
         }
         if (psuadeIO_ == NULL) printf("ERROR: data not loaded yet.\n");
         else
         {
            faType = -1;
            sprintf(pString, "Enter your choice ? ");
            while (faType < 0 || faType > PSUADE_NUM_RS)
            {
               writeFAInfo();
               faType = getFAType(pString);
            }
            outputID = 0;
            if (nOutputs > 1)
            {
               sprintf(pString, "Enter output number (1 - %d) = ", nOutputs);
               outputID = getInt(1, nOutputs, pString);
               outputID--;
            }
            sprintf(pString, "How many groups (2 - %d)? ", nSamples);
            nParts = getInt(2, nSamples, pString);
            count = nSamples / nParts;
            tempX = new double[nSamples*nInputs];
            tempY = new double[nSamples];
            tempV = new double[nSamples];
            for (ss = 0; ss < nSamples; ss+=count)
            { 
               for (kk = 0; kk < ss*nInputs; kk++) tempX[kk] = sampleInputs[kk];
               for (kk = (ss+count)*nInputs; kk < nSamples*nInputs; kk++)
                  tempX[kk-count*nInputs] = sampleInputs[kk];
               for (kk = 0; kk < ss; kk++)
                  tempY[kk] = sampleOutputs[kk*nOutputs+outputID];
               for (kk = ss+count; kk < nSamples; kk++)
                  tempY[kk-count] = sampleOutputs[kk*nOutputs+outputID];
               count2 = ss;
               if (ss+count < nSamples) count2 += (nSamples - ss - count);
               faPtr = genFA(faType, nInputs, iOne, count2);
               faPtr->setBounds(iLowerB, iUpperB);
               faPtr->setOutputLevel(outputLevel_);
               faLeng = -999;
               status = faPtr->genNDGridData(tempX,tempY,&faLeng,NULL,NULL);
               count2 = nSamples - count2;
               faPtr->evaluatePoint(count2, &sampleInputs[ss*nInputs], &tempV[ss]);
               delete faPtr;
            }
            if (psPlotTool_ == 1)
            {
               fp = fopen("RSTest_cv.sci", "w");
               fprintf(fp, "// col 1: simulation data, col 2: rs data\n");
            }
            else
            {
               fp = fopen("RSTest_cv.m", "w");
               fprintf(fp, "%% col 1: simulation data, col 2: rs predicted data\n");
            }
            fprintf(fp, "A = [\n");
            for (ss = 0; ss < nSamples; ss++)
               fprintf(fp, "%e %e\n",sampleOutputs[ss*nOutputs+outputID],tempV[ss]);
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
               fprintf(fp, "bar(xk,nk/%d,1.0)\n",nSamples);
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
            for (ss = 0; ss < nSamples; ss++)
               tempY[ss] = sampleOutputs[ss*nOutputs+outputID];
            ddata = 0.0;
            for (ss = 0; ss < nSamples; ss++)
               ddata += (tempY[ss] - tempV[ss]);
            ddata = ddata / nSamples;
            printf("CV avg error (unscaled) = %e\n", ddata);
            ddata = 0.0;
            for (ss = 0; ss < nSamples; ss++)
            {
               ddata += (tempY[ss] - tempV[ss]);
               if (tempY[ss] != 0.0)
               {
                  ddata /= PABS(tempY[ss]);
                  dtemp = tempY[ss];
               }
            }
            ddata = ddata / nSamples;
            printf("CV avg error (  scaled) = %e (base=%e)\n",ddata,dtemp);
            ddata = 0.0;
            for (ss = 0; ss < nSamples; ss++)
               ddata += pow(tempY[ss] - tempV[ss], 2.0);
            ddata = sqrt(ddata / nSamples);
            printf("CV rms error (unscaled) = %e\n", ddata);
            ddata = 0.0;
            for (ss = 0; ss < nSamples; ss++)
            {
               ddata += pow(tempY[ss] - tempV[ss], 2.0);
               if (tempY[ss] != 0.0)
               {
                  ddata /= pow(tempY[ss],2.0);
                  dtemp = tempY[ss];
               }
            }
            ddata = sqrt(ddata / nSamples);
            printf("CV rms error (  scaled) = %e (base=%e)\n",ddata,dtemp);
            ddata = 0.0;
            for (ss = 0; ss < nSamples; ss++)
               if (PABS(tempV[ss]-tempY[ss]) > ddata)
                  ddata = PABS(tempY[ss] - tempV[ss]);
            printf("CV max error (unscaled) = %e\n", ddata);
            ddata = 0.0;
            for (ss = 0; ss < nSamples; ss++)
            {
               if (tempY[ss] != 0 && PABS((tempV[ss]-tempY[ss])/tempY[ss])>ddata)
               {
                  dtemp = tempY[ss];
                  ddata = PABS((tempY[ss]-tempV[ss])/tempY[ss]);
               }
            }
            printf("CV max error (  scaled) = %e (base=%e)\n",ddata,dtemp);
            ddata = 0.0;
            for (ss = 0; ss < nSamples; ss++)
               ddata += (tempY[ss] - tempV[ss]);
            ddata = ddata / nSamples;
            printf("CV error mean = %e\n", ddata);
            dtemp = 0.0;
            for (ss = 0; ss < nSamples; ss++)
               dtemp += pow(tempY[ss] - tempV[ss] - ddata, 2.0);
            dtemp = sqrt(dtemp / (nSamples - 1.0));
            printf("CV error std  = %e\n", dtemp);
            if (psPlotTool_ == 1)
                 printf("rstest_cv plot file in RSTest_cv.sci\n");
            else printf("rstest_cv plot file in RSTest_cv.m\n");
            delete [] tempX;
            delete [] tempY;
            delete [] tempV;
            faPtr = NULL;
            tempV = NULL;
            tempY = NULL;
            tempX = NULL;
         }
      }

      // response surface test - generalization test
      // ----------------------------------------------------------------
      else if (!strcmp(command, "rstest_gt") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rstest_gt: response surface generalization test.\n");
            printf("syntax: rstest_gt (no argument needed)\n");
            continue;
         }
         printf("ERROR: no data to analyze (load data first).\n");
      }
      else if (!strcmp(command, "rstest_gt") && nInputs > 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rstest_gt: response surface generalization test.\n");
            printf("syntax: rstest_gt (no argument needed)\n");
            continue;
         }
         if (psuadeIO_ == NULL) printf("ERROR: data not loaded yet.\n");
         else
         {
            faType = -1;
            sprintf(pString, "Enter your choice ? ");
            while (faType < 0 || faType > PSUADE_NUM_RS)
            {
               writeFAInfo();
               faType = getFAType(pString);
            }
            outputID = 0;
            if (nOutputs > 1)
            {
               sprintf(pString, "Enter output number (1 - %d) = ", nOutputs);
               outputID = getInt(1, nOutputs, pString);
               outputID--;
            }
            thresh = -1.0;
            while (thresh < 0.5 || thresh >= 1.0)
            {
               printf("This test partitions each input into two parts and uses\n");
               printf("one part to predict the other part. You can decide the\n");
               printf("percentage of the sample points for training (>= 50 percent).\n");
               sprintf(pString, "Enter this percentage [0.5 - 1) = ");
               thresh = getDouble(pString);
            }
            if (psPlotTool_ == 1)
            {
               fp = fopen("RSTest_gt.sci", "w");
               fprintf(fp, "// col 1-m: inputs, col m+1: simulation output\n");
               fprintf(fp, "// col m+2: predicted output \n");
            }
            else
            {
               fp = fopen("RSTest_gt.m", "w");
               fprintf(fp, "%% col 1-m: inputs, col m+1: simulation output\n");
               fprintf(fp, "%% col m+2: predicted output \n");
            }
            fwritePlotCLF(fp);

            int    *acnts;
            double *ameans, *avars;
            tempX = new double[nSamples*nInputs];
            tempV = new double[nSamples*nInputs];
            tempY = new double[nSamples];
            tempW = new double[nSamples];
            tempT = new double[nSamples];
            ameans = new double[nInputs*2];
            avars  = new double[nInputs*2];
            acnts  = new int[nInputs*2];
            printf("Generalization test results:\n");
            for (ii = 0; ii < nInputs; ii++)
            {
               fprintf(fp, "A%d = [\n", ii+1);
               ddata = thresh * (iUpperB[ii] - iLowerB[ii]) + iLowerB[ii];
               count = count2 = 0;
               for (ss = 0; ss < nSamples; ss++)
               { 
                  if (sampleInputs[ss*nInputs+ii] <= ddata)
                  {
                     for (kk = 0; kk < nInputs; kk++)
                        tempX[count*nInputs+kk] = sampleInputs[ss*nInputs+kk];
                     tempY[count] = sampleOutputs[ss*nOutputs+outputID];
                     count++;
                  }
                  else
                  {
                     for (kk = 0; kk < nInputs; kk++)
                        tempV[count2*nInputs+kk] = sampleInputs[ss*nInputs+kk];
                     tempW[count2] = sampleOutputs[ss*nOutputs+outputID];
                     count2++;
                  }
               }
               if (count2 == 0 || count == 0)
               {
                  printf("ERROR: for input %d - no training or test point.\n", ii+1);
                  continue;
               }
               faPtr = genFA(faType, nInputs, iOne, count);
               faPtr->setBounds(iLowerB, iUpperB);
               faPtr->setOutputLevel(outputLevel_);
               faLeng = -999;
               status = faPtr->genNDGridData(tempX,tempY,&faLeng,NULL,NULL);
               faPtr->evaluatePoint(count2, tempV, tempT);
               delete faPtr;
               for (ss = 0; ss < count2; ss++)
               {
                  for (kk = 0; kk < nInputs; kk++)
                     fprintf(fp, "%e ", tempV[ss*nInputs+kk]);
                  fprintf(fp, "%e ", tempW[ss]);
                  fprintf(fp, "%e\n", tempT[ss]);
               }
               ameans[2*ii] = 0.0;
               for (ss = 0; ss < count2; ss++) ameans[2*ii] += (tempW[ss] - tempT[ss]);
               ameans[2*ii] /= (double) count2;
               avars[2*ii] = 0.0;
               for (ss = 0; ss < count2; ss++) 
                  avars[2*ii] += pow(tempW[ss]-tempT[ss]-ameans[2*ii],2.0);
               avars[2*ii] /= (double) count2;
               acnts[2*ii] = count2;

               ddata = (1.0 - thresh) * (iUpperB[ii] - iLowerB[ii]) + iLowerB[ii];
               count = count2 = 0;
               for (ss = 0; ss < nSamples; ss++)
               { 
                  if (sampleInputs[ss*nInputs+ii] > ddata)
                  {
                     for (kk = 0; kk < nInputs; kk++)
                        tempX[count*nInputs+kk] = sampleInputs[ss*nInputs+kk];
                     tempY[count] = sampleOutputs[ss*nOutputs+outputID];
                     count++;
                  }
                  else
                  {
                     for (kk = 0; kk < nInputs; kk++)
                        tempV[count2*nInputs+kk] = sampleInputs[ss*nInputs+kk];
                     tempW[count2] = sampleOutputs[ss*nOutputs+outputID];
                     count2++;
                  }
               }
               if (count2 == 0 || count == 0)
               {
                  printf("ERROR: for input %d - no training or test point.\n",ii+1);
                  continue;
               }
               faPtr = genFA(faType, nInputs, iOne, count);
               faPtr->setBounds(iLowerB, iUpperB);
               faPtr->setOutputLevel(outputLevel_);
               faLeng = -999;
               status = faPtr->genNDGridData(tempX,tempY,&faLeng,NULL,NULL);
               faPtr->evaluatePoint(count2, tempV, tempT);
               delete faPtr;
               for (ss = 0; ss < count2; ss++)
               {
                  for (kk = 0; kk < nInputs; kk++)
                     fprintf(fp, "%e ", tempV[ss*nInputs+kk]);
                  fprintf(fp, "%e ", tempW[ss]);
                  fprintf(fp, "%e\n", tempT[ss]);
               }
               fprintf(fp, "];\n");
               ameans[2*ii+1] = 0.0;
               for (ss = 0; ss < count2; ss++) ameans[2*ii+1] += (tempW[ss] - tempT[ss]);
               ameans[2*ii+1] /= (double) count2;
               avars[2*ii+1] = 0.0;
               for (ss = 0; ss < count2; ss++) 
                  avars[2*ii+1] += pow(tempW[ss]-tempT[ss]-ameans[2*ii+1],2.0);
               avars[2*ii+1] /= (double) count2;
               acnts[2*ii+1] = count2;
               fprintf(fp, "xmax = max(A%d(:,%d+1));\n", ii+1, nInputs);
               fprintf(fp, "xmin = min(A%d(:,%d+1));\n", ii+1, nInputs);
               fprintf(fp, "ymax = max(A%d(:,%d+2));\n", ii+1, nInputs);
               fprintf(fp, "ymin = min(A%d(:,%d+2));\n", ii+1, nInputs);
               fprintf(fp, "xmin = min(xmin, ymin);\n");
               fprintf(fp, "xmax = max(xmax, ymax);\n");
               fprintf(fp, "XX   = xmin : xmax-xmin : xmax;\n");
               fprintf(fp, "plot(A%d(:,%d+1), A%d(:,%d+2),'*', XX, XX)\n",
                       ii+1, nInputs, ii+1, nInputs);
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
               if (ii < nInputs-1) fprintf(fp, "pause\n");
            }
            for (ii = 0; ii < nInputs; ii++)
            {
               printf("Input %4d: partition 1 error mean/variance = %12.4e %12.4e (%d)\n",
                      ii+1, ameans[2*ii], avars[2*ii], acnts[2*ii]); 
               printf("Input %4d: partition 2 error mean/variance = %12.4e %12.4e (%d)\n",
                      ii+1, ameans[2*ii+1], avars[2*ii+1], acnts[2*ii]+1); 
            }
            fclose(fp);
            if (psPlotTool_ == 1)
                 printf("rstest_gt plot file in RSTest_gt.sci\n");
            else printf("rstest_gt plot file in RSTest_gt.m\n");
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
      }

      // find best RS parameters for SVM
      // ----------------------------------------------------------------
      else if (!strcmp(command, "svmfind") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("svmfind: find optimal parameters for SVM surface fit\n");
            printf("syntax: svmfind (no argument needed)\n");
            continue;
         }
         printf("ERROR: no data to analyze.\n");
      }
      else if (!strcmp(command, "svmfind") && nInputs > 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("svmfind: find optimal parameters for SVM surface fit\n");
            printf("syntax: svmfind (no argument needed)\n");
            continue;
         }
         if (psuadeIO_ == NULL) printf("ERROR: data not loaded yet.\n");
         else
         {
#ifndef HAVE_SVM
            printf("SVM not installed.\n");
#else
            outputID = 0;
            if (nOutputs > 1)
            {
               sprintf(pString, "Enter output number (1 - %d) = ", nOutputs);
               outputID = getInt(1, nOutputs, pString);
               outputID--;
            }
            faType = PSUADE_RS_SVM;
            kk = psAnalysisInteractive_;
            psAnalysisInteractive_ = 0;
            jj = psRSExpertMode_;
            psRSExpertMode_ = 0;
            faPtr = genFA(faType, nInputs, iOne, nSamples);
            psAnalysisInteractive_ = kk;
            psRSExpertMode_ = jj;
            faPtr->setBounds(iLowerB, iUpperB);
            tempY = new double[nSamples];
            for (ii = 0; ii < nSamples; ii++)
               tempY[ii] = sampleOutputs[nOutputs*ii+outputID];
            status = -999;
            gamma = 1.0e-6;
            kk = 3;
            if (psRSExpertMode_ == 1)
            {
               printf("SVM kernel options: \n");
               printf("1. linear\n");
               printf("2. third order polynomial\n");
               printf("3. radial basis function\n");
               printf("4. sigmoid function\n");
               sprintf(pString, "kernel selection (1-4): ");
               kk = getInt(1, 4, pString);
            }
            kk--;
            printf("SVM tolerance can be adjusted for more accurate result.\n");
            printf("However, small tolerance takes more time to run.\n");
            tolerance = 0.0;
            sprintf(pString, "Enter the desired tolerance (0.1 - 1e-6): ");
            while (tolerance > 0.1 || tolerance < 1.0e-6)
               tolerance = getDouble(pString);
            while (gamma <= 1.0e6)
            {
               targv[0] = (char *) &gamma;
               targv[1] = (char *) &tolerance;
               targv[2] = (char *) &kk;
               faPtr->setParams(3, targv);
               faPtr->genNDGridData(sampleInputs, tempY, &status, NULL, NULL);
               sumErr = maxErr = 0.0;
               for (ii = 0; ii < nSamples; ii++)
               {
                  ddata = sampleOutputs[ii*nOutputs+outputID];
                  dtemp = faPtr->evaluatePoint(&sampleInputs[ii*nInputs]);
                  dtemp = PABS(dtemp - ddata);
                  if (dtemp > maxErr) maxErr = dtemp;
                  sumErr += (dtemp * dtemp);
               }
               sumErr = sqrt(sumErr);
               printf("svmfind: RBF_gamma = %e, error (L2n,max) = %e %e\n",
                      gamma, sumErr, maxErr);
               if (gamma == 1.0e-6 || sumErr < Ymin)
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
      }

      // check quality of the response surface using another data set 
      // ----------------------------------------------------------------
      else if ((!strcmp(command, "rstest") ||
                !strcmp(command, "rstest_hs")) && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rstest_hs: check quality of RS with a holdout set\n");
            printf("syntax: rstest_hs (no argument needed)\n");
            continue;
         }
         printf("ERROR: no data to build response surface.\n");
      }
      else if ((!strcmp(command, "rstest") ||
                !strcmp(command, "rstest_hs")) && nInputs > 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rstest_hs: check quality of RS with a holdout set\n");
            printf("syntax: rstest_hs (no argument needed)\n");
            continue;
         }
         if (psuadeIO_ == NULL) printf("ERROR: data not loaded yet.\n");
         else
         {
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
            if (status == 0)
            {
               ioPtr->getParameter("output_noutputs", pPtr);
               kk = pPtr.intData_;
               if (kk > 1)
               {
                  printf("ERROR: your test data file should only have one output\n");
                  printf("       per sample point. Fix this and do this again.\n");
                  status = 1;
               }
               delete ioPtr;
               ioPtr = NULL;
            }
            if (status == 0)
            {
               outputID = 0;
               if (nOutputs > 1)
               {
                  sprintf(pString,"Enter output number (1 - %d) = ",nOutputs);
                  outputID = getInt(1, nOutputs, pString);
                  outputID--;
               }
               faType = -1;
               faLimit = 9;
               while (faType < 0 || faType >= PSUADE_NUM_RS)
               {
                  printf("Enter a response surface method : \n");
                  faLimit = writeFAInfo();
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
               printf("Have a discrepancy model file? If so, enter name (enter n otherwise): ");
               scanf("%s", winput);
               fgets(lineIn,500,stdin); 
               if (!strcmp(winput, "NONE") || winput[0] == 'n') discFile = 0;
               else
               {
                  ioPtr = new PsuadeData();
                  status = ioPtr->readPsuadeFile(winput);
                  if (status == 0)
                  {
                     ioPtr->getParameter("input_ninputs", pPtr);
                     nInps = pPtr.intData_;
                     if (nInps < nInputs)
                     {
                        printf("ERROR: your sample has %d inputs but\n", nInputs);
                        printf("       your discrepancy model has %d inputs.\n", nInps);
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
                     printf("ERROR: in reading the discrepancy model file %s.\n",winput);
                     discFile = 0;
                     delete ioPtr;
                     ioPtr = NULL;
                     continue;
                  }
               }
               double *samOuts = new double[nSamples];
               for (ii = 0; ii < nSamples; ii++) 
                  samOuts[ii] = sampleOutputs[ii*nOutputs+outputID];
               if (discFile == 1)
               {
                  if (outputLevel_ > 4)
                  {
                     fp = fopen("rstest_hs_data", "w");
                     if (fp != NULL)
                        fprintf(fp,"%% Col 2: RS estimate, Col 3: model correction.\n");
                  }
                  else fp = NULL; 
                     
                  for (ii = 0; ii < nSamples; ii++)
                  {
                     dtemp = faPtrsRsEval[0]->evaluatePoint(&sampleInputs[ii*nInputs]);
                     if (fp != NULL) 
                     {
                        fprintf(fp, "%6d ",ii+1);
                        for (jj = 0; jj < nInputs; jj++)
                           fprintf(fp, "%12.4e ",sampleInputs[ii*nInputs+jj]);
                        fprintf(fp, " %12.4e %12.4e\n",samOuts[ii],dtemp);
                     }
                     samOuts[ii] += dtemp;
                  }
                  if (fp != NULL)
                  {
                     printf("rstest_hs_data contains simulations plus discrepancies.\n");
                     fclose(fp);
                  }
               }
               analysisMethod = PSUADE_ANA_RSFA;
               anaManager = new AnalysisManager();
               anaManager->setup(analysisMethod, faType);
               anaManager->loadLogXsformFlags(2, xsforms);
               aPtr.printLevel_ = outputLevel_;
               aPtr.sampleInputs_ = sampleInputs;
               aPtr.sampleOutputs_ = samOuts;
               aPtr.nSamples_ = nSamples;
               aPtr.nInputs_ = nInputs;
               aPtr.nOutputs_ = 1;
               aPtr.iLowerB_ = iLowerB;
               aPtr.iUpperB_ = iUpperB;
               aPtr.outputID_ = 0;
               aPtr.faType_ = faType;
               aPtr.regWgtID_ = -1;
               aPtr.ioPtr_ = NULL;
               aPtr.sampleStates_ = sampleStates;
               strcpy(pString, "validate");
               targv[0] = (char *) pString;
               targv[1] = (char *) &aPtr;
               targv[2] = (char *) dataFile;
               if (psPlotTool_ == 1) strcpy(errFile, "psuade_error_file.sci");
               else                  strcpy(errFile, "psuade_error_file.m");
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
            }
         }
      }

      // generate sample for response surface test
      // ----------------------------------------------------------------
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
         printEquals(0);
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
            count = (int) pow(2.0, 1.0*nInputs);
            if (count > 10000)
               printf("nSamples = %d may be too large.\n", nSamples);
         }
         else if (ind == 2)
         {
            sampPtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_FF5);
            count = 128;
            if (nInputs > 11)
            {
               printf("nInputs = %d too large for FF5.\n", nInputs);
               continue;
            }
         }
         else if (ind == 3)
         {
            sampPtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_FF4);
            count = 64;
            if (nInputs > 32)
            {
               printf("nInputs = %d too large for FF4.\n", nInputs);
               continue;
            }
         }
         sampPtr->setPrintLevel(0);
         sampPtr->setInputBounds(nInputs, iLowerB, iUpperB);
         sampPtr->setOutputParams(nOutputs);
         sampPtr->setSamplingParams(count, -1, 0);
         sampPtr->initialize(0);
         count = sampPtr->getNumSamples();
         tempX  = new double[count * nInputs];
         tempY  = new double[count * nOutputs];
         states = new int[count];
         sampPtr->getSamples(count,nInputs,nOutputs,tempX,tempY,states);
         delete sampPtr;
         sampPtr = NULL;
         ioPtr = new PsuadeData();
         ioPtr->updateInputSection(count, nInputs, NULL, iLowerB,
                                   iUpperB, tempX, inputNames); 
         ioPtr->updateOutputSection(count, nOutputs, tempY, states, 
                                    outputNames); 
         ioPtr->updateMethodSection(PSUADE_SAMP_MC, count, 1, -1, -1);
         ioPtr->writePsuadeFile("psuade_rs_test.data",0);
         printf("The test sample is in file psuade_rs_test.data.\n"); 
         delete [] tempX;
         delete [] tempY;
         delete [] states;
         delete ioPtr;
      }

      // numerical integration
      // ----------------------------------------------------------------
      else if (!strcmp(command, "int") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("int: numerical integration(volume under sample surface)\n");
            printf("syntax: int (no argument needed)\n");
            continue;
         }
         printf("ERROR: no data to integrate (load data first).\n");
      }
      else if (!strcmp(command, "int") && nInputs > 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("int: numerical integration(volume under sample surface)\n");
            printf("syntax: int (no argument needed)\n");
            continue;
         }
         if (psuadeIO_ == NULL) printf("ERROR: data not loaded yet.\n");
         else
         {
            outputID = 0;
            if (nOutputs > 1)
            {
               sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
               outputID = getInt(1, nOutputs, pString);
               outputID--;
            }
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
      }

      // calculate volume of constrained space
      // ----------------------------------------------------------------
      else if (!strcmp(command, "rsvol") && nInputs <= 0)
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
         printf("ERROR: no data (load data first).\n");
      }
      else if (!strcmp(command, "rsvol") && nInputs > 0)
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
         if (psuadeIO_ == NULL) printf("ERROR: data not loaded yet.\n");
         else
         {
            outputID = 0;
            if (nOutputs > 1)
            {
               sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
               outputID = getInt(1, nOutputs, pString);
               outputID--;
            }
            Ymax = - 1.0e35;
            Ymin =   1.0e35;
            for (sInd = 0; sInd < nSamples; sInd++)
            {
               if (sampleOutputs[sInd*nOutputs+outputID] > Ymax)
                  Ymax = sampleOutputs[sInd*nOutputs+outputID];
               if (sampleOutputs[sInd*nOutputs+outputID] < Ymin)
                  Ymin = sampleOutputs[sInd*nOutputs+outputID];
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
            if (nInputs <= 51) 
                 sampPtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_LPTAU);
            else sampPtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_MC);
            sampPtr->setPrintLevel(0);
            sampPtr->setInputBounds(nInputs, iLowerB, iUpperB);
            sampPtr->setOutputParams(1);
            kk = 1000000;
            sampPtr->setSamplingParams(kk, 1, 0);
            sampPtr->initialize(0);
            tempX = new double[kk*nInputs];
            tempY = new double[kk];
            states = new int[kk];
            sampPtr->getSamples(kk, nInputs, 1, tempX, tempY, states);
            count = 0;
            for (ii = 0; ii < kk; ii++)
            {
               dtemp = faPtr->evaluatePoint(&tempX[ii*nInputs]);
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
      }

      // uncertainty analysis on fuzzy response surface
      // ----------------------------------------------------------------
      else if (!strcmp(command, "rs_ua") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rs_ua: uncertainty analysis on fuzzy response surface\n");
            printf("syntax: rs_ua (no argument needed)\n");
            printf("This command perform uncertainty analysis on the response\n");
            printf("surface built from the loaded sample. It uses 100\n");
            printf("instantiations and computes CDF for each one of them\n");
            printf("for plotting (graphics output only).\n");
            continue;
         }
         printf("ERROR: no data (load data first).\n");
      }
      else if (!strcmp(command, "rs_ua") && nInputs > 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rs_ua: uncertainty analysis on fuzzy response surface\n");
            printf("syntax: rs_ua (no argument needed)\n");
            printf("This command perform uncertainty analysis on the response\n");
            printf("surface built from the loaded sample. It uses 100\n");
            printf("instantiations and computes CDF for each one of them\n");
            printf("for plotting (graphics output only).\n");
         }
         else
         {
            outputID = 0;
            sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
            outputID = getInt(1, nOutputs, pString);
            outputID--;
            int nSamp = 1000;
            sprintf(pString, "MC Sample size ? (10000 - 100000) ");
            if (nSamp < 10000 || nSamp > 100000)
               nSamp = getInt(10000, 100000, pString);
            flag = 0;
            printf("Would you like to save the sample ? (y or n) ");
            fgets(winput,10,stdin); 
            if (winput[0] == 'y') flag = 1; 

            printf("Phase 1 of 4: create response surface\n");
            faType = -1;
            faPtr = genFA(faType, nInputs, iOne, nSamples);
            faPtr->setNPtsPerDim(32);
            faPtr->setBounds(iLowerB, iUpperB);
            faPtr->setOutputLevel(0);
            if (nOutputs > 1)
            {
               tempY = new double[nSamples];
               for (ss = 0; ss < nSamples; ss++) 
                  tempY[ss] = sampleOutputs[ss*nOutputs+outputID];
            }
            else tempY = sampleOutputs;
            count = -999;
            faPtr->genNDGridData(sampleInputs,tempY,&count,NULL,NULL);
            if (nOutputs > 1) delete [] tempY;
            tempY = NULL;
            
            printEquals(0);
            printf("Phase 2 of 4: create MC sample\n");
            int    *samStates  = new int[nSamp];
            double *samInputs  = new double[nInputs*nSamp];
            double *samOutputs = new double[nSamp];
            double *samStds    = new double[nSamp];
            sampPtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_MC);
            sampPtr->setPrintLevel(0);
            sampPtr->setInputBounds(nInputs, iLowerB, iUpperB);
            sampPtr->setOutputParams(1);
            sampPtr->setSamplingParams(nSamp, -1, -1);
            sampPtr->initialize(0);
            sampPtr->getSamples(nSamp,nInputs,1,samInputs,samOutputs,samStates);

            printf("Phase 3 of 4: evaluate sample\n");
            faPtr->evaluatePointFuzzy(nSamp,samInputs,samOutputs,samStds); 
            delete [] samStates;
            delete faPtr;
            delete sampPtr;
            samStates = NULL;
            faPtr = NULL;
            sampPtr = NULL;

            int nbins = 20, ntimes=10;
            int **Fcounts = new int*[ntimes];
            double Fmax=-PSUADE_UNDEFINED;
            double Fmin=PSUADE_UNDEFINED;
            PDFNormal *rsPDF;

            printf("Phase 4 of 4: binning\n");
            for (ss = 0; ss < nSamp; ss++)
            {
               if (samOutputs[ss] > Fmax) Fmax = samOutputs[ss];
               if (samOutputs[ss] < Fmin) Fmin = samOutputs[ss];
               if (samOutputs[ss]+3*samStds[ss] > Fmax) 
                  Fmax = samOutputs[ss] + 3 * samStds[ss];
               if (samOutputs[ss]-3*samStds[ss] < Fmin) 
                  Fmin = samOutputs[ss] - 3 * samStds[ss];
            }
            Fmax = Fmax + 0.1 * (Fmax - Fmin);
            Fmin = Fmin + 0.1 * (Fmax - Fmin);
            if (Fmax == Fmin)
            {
               Fmax = Fmax + 0.1 * PABS(Fmax);
               Fmin = Fmin - 0.1 * PABS(Fmin);
            }
            for (ii = 0; ii < ntimes; ii++)
            {
               Fcounts[ii] = new int[nbins];
               for (kk = 0; kk < nbins; kk++) Fcounts[ii][kk] = 0;
            }
            if (flag == 1)
            {
               fp = fopen("rsua_sample","w");
               fprintf(fp, "%% inputs, output, output-3 sigma, output+3sigma\n");
               fprintf(fp, "%d %d 3\n", nSamp, nInputs);
            }
            double mean=0, stdev=0;
            double *samFuzzy = new double[ntimes*nInputs];
            for (ss = 0; ss < nSamp; ss++)
            {
               mean += samOutputs[ss];
               if (samStds[ss] != 0)
               {
                  rsPDF = new PDFNormal(samOutputs[ss],samStds[ss]);
                  rsPDF->genSample(ntimes,samFuzzy,
                                   samOutputs[ss]-3*samStds[ss],
                                   samOutputs[ss]+3*samStds[ss]);
               }
               else
               {
                  for (ii = 0; ii < ntimes; ii++) 
                     samFuzzy[ii] = samOutputs[ss];
               }
               if (flag == 1)
               {
                  for (ii = 0; ii < nInputs; ii++) 
                     fprintf(fp, "%e ", samInputs[ss*nInputs+ii]);
                  fprintf(fp, "%e ", samOutputs[ss]);
                  fprintf(fp, "%e ", samOutputs[ss]-3*samStds[ss]);
                  fprintf(fp, "%e\n", samOutputs[ss]+3*samStds[ss]);
               }
               for (ii = 0; ii < ntimes; ii++) 
               {
                  ddata = samFuzzy[ii] - Fmin;
                  if (Fmax > Fmin)
                     ddata = ddata / ((Fmax - Fmin) / nbins);
                  else ddata = nbins / 2;
                  kk = (int) ddata;
                  if (kk < 0)      kk = 0;
                  if (kk >= nbins) kk = nbins - 1;
                  Fcounts[ii][kk]++;
               }
               if (ss % (nSamp / 8) == 0)
               {
                  printf(".");
                  fflush(stdout);
               }
            }
            mean /= (double) nSamp;
            for (ss = 0; ss < nSamp; ss++)
               stdev += pow(samOutputs[ss]-mean, 2.0);
            stdev = sqrt(stdev/(double) nSamp);
            printf("\n");
            printAsterisks(0);
            printf("Sample mean    = %e\n", mean);
            printf("Sample std dev = %e\n", stdev);
            printAsterisks(0);
            if (flag == 1)
            {
               fclose(fp);
               printf("A MC sample has been written to the file 'rsua_sample'.\n");
            }
            delete [] samStds;
            delete [] samInputs;
            delete [] samOutputs;
            delete [] samFuzzy;

            if (psPlotTool_ == 1)
            {
               fp = fopen("scilabrsua.sci", "w");
               if (fp == NULL)
               {
                  printf("rs_ua ERROR: cannot open scilab file.\n");
                  continue;
               }
            }
            else
            {
               fp = fopen("matlabrsua.m", "w");
               if (fp == NULL)
               {
                  printf("rs_ua ERROR: cannot open matlab file.\n");
                  continue;
               }
            }
            fwritePlotCLF(fp);
            fprintf(fp, "twoPlots = 1;\n");
            fprintf(fp, "if (twoPlots == 1)\n");
            fprintf(fp, "subplot(1,2,1)\n");
            fprintf(fp, "end;\n");
            fprintf(fp, "X = [\n");
            for (kk = 0; kk < nbins; kk++)
               fprintf(fp, "%e\n", (Fmax-Fmin)/nbins*(0.5+kk)+Fmin);
            fprintf(fp, "];\n");
            for (ii = 0; ii < ntimes; ii++)
            {
               fprintf(fp, "N%d = [\n", ii+1);
               for (kk = 0; kk < nbins; kk++)
                  fprintf(fp, "%d\n",  Fcounts[ii][kk]);
               fprintf(fp, "];\n");
            }
            fprintf(fp, "N = [");
            for (ii = 0; ii < ntimes; ii++)
               fprintf(fp, "N%d/sum(N%d) ", ii+1, ii+1);
            fprintf(fp, "];\n");
            fprintf(fp, "bar(X,N)\n");
            fwritePlotAxes(fp);
            fwritePlotTitle(fp, "Probability Distributions");
            fwritePlotXLabel(fp, "Output Value");
            fwritePlotYLabel(fp, "Probabilities");
            fprintf(fp, "if (twoPlots == 1)\n");
            for (ii = 0; ii < ntimes; ii++)
            {
               fprintf(fp, "for ii = 2 : %d\n", nbins);
               fprintf(fp, "   N%d(ii) = N%d(ii) + N%d(ii-1);\n",ii+1,ii+1,ii+1);
               fprintf(fp, "end;\n");
            }
            fprintf(fp, "N = [");
            for (ii = 0; ii < ntimes; ii++)
               fprintf(fp, "N%d/N%d(%d) ", ii+1, ii+1, nbins);
            fprintf(fp, "];\n");
            fprintf(fp, "subplot(1,2,2)\n");
            fprintf(fp, "bar(X,N)\n");
            fwritePlotAxes(fp);
            fwritePlotTitle(fp, "Cumulative Distribution");
            fwritePlotXLabel(fp, "Output Values");
            fwritePlotYLabel(fp, "Probabilities");
            fprintf(fp, "end;\n");
            fclose(fp);
            if (psPlotTool_ == 1)
                 printf("Output distribution plots is in scilabrsua.sci.\n");
            else printf("Output distribution plots is in matlabrsua.m.\n");
            printf("Note: The different colors in a given bar correspond\n");
            printf("      to probabilities with different bootstrapped sample.\n");
            for (ii = 0; ii < ntimes; ii++) delete [] Fcounts[ii];
            delete [] Fcounts;
         }
      }

      // principal component analysis 
      // ----------------------------------------------------------------
      else if (!strcmp(command, "pca") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("pca: principal component analysis on the outputs\n");
            printf("syntax: pca (no argument needed)\n");
            continue;
         }
         printf("ERROR: no data to analyze (load data first).\n");
      }
      else if (!strcmp(command, "pca") && nInputs > 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("pca: principal component analysis on the outputs\n");
            printf("syntax: pca (no argument needed)\n");
            continue;
         }
         if (psuadeIO_ == NULL) printf("ERROR: data not loaded yet.\n");
         else
         {
            aPtr.nSamples_ = nSamples;
            aPtr.nInputs_ = nInputs;
            aPtr.nOutputs_ = nOutputs;
            aPtr.iLowerB_ = iLowerB;
            aPtr.iUpperB_ = iUpperB;
            aPtr.sampleInputs_ = sampleInputs;
            aPtr.sampleOutputs_ = sampleOutputs;
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
            if (kk != nOutputs) 
            {
               if (outputNames != NULL) 
               {
                  for (ii = 0; ii < nOutputs; ii++)
                     delete [] outputNames[ii];
                  delete [] outputNames;
               }
               nOutputs = kk;
               pONames.clean();
               psuadeIO_->getParameter("output_names", pONames);
               names = pONames.strArray_;
               outputNames = new char*[nOutputs+1];
               for (ii = 0; ii < nOutputs; ii++)
               {
                  outputNames[ii] = new char[200]; 
                  strcpy(outputNames[ii], names[ii]);
               }
               if (sampleOutputs != NULL) delete [] sampleOutputs;
               psuadeIO_->getParameter("output_sample", pPtr);
               sampleOutputs   = pPtr.dbleArray_;
               pPtr.dbleArray_ = NULL;
            }
         }
      }

      // One sample test (chi-squred, distribution fitting) 
      // ----------------------------------------------------------------
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

      // 2 sample test (Kolmogorov Smirnov, Mann-Whitney, T-test) 
      // ----------------------------------------------------------------
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
         fgets(lineIn,500,stdin); 
      }

      // 1D input scatter plot
      // ----------------------------------------------------------------
      else if (!strcmp(command, "iplot1") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("iplot1: plot the sample points in one parameter space\n");
            printf("syntax: iplot1 (no argument needed)\n");
            continue;
         }
         printf("ERROR: no data to analyze (load data first).\n");
      }
      else if (!strcmp(command, "iplot1") && nInputs >= 1)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("iplot1: plot the sample points in one parameter space\n");
            printf("syntax: iplot1 (no argument needed)\n");
            continue;
         }
         if (psuadeIO_ == NULL) printf("ERROR: data not loaded yet.\n");
         else
         {
            iplot1 = -1;
            sprintf(pString, "Enter the input number (1 - %d) : ", nInputs);
            iplot1 = getInt(1, nInputs, pString);
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
               for (sInd = 0; sInd < nSamples; sInd++)
               {
                  if (sampleOutputs[sInd*nOutputs] > 0.5*PSUADE_UNDEFINED ||
		      sampleStates[sInd] != 1)
                     fprintf(fp, "%e 0\n",sampleInputs[sInd*nInputs+iplot1]);
                  else
                     fprintf(fp, "%e 1\n",sampleInputs[sInd*nInputs+iplot1]);
               }
               fprintf(fp, "];\n");
               fprintf(fp, "n = %d;\n", nSamples);
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
               fprintf(fp, "a.data_bounds=[0,minX;%d,maxX];\n",nSamples);
               fwritePlotYLabel(fp, inputNames[iplot1]);
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
               for (sInd = 0; sInd < nSamples; sInd++)
               {
                  if (sampleOutputs[sInd*nOutputs] > 0.5*PSUADE_UNDEFINED ||
		      sampleStates[sInd] != 1)
                     fprintf(fp, "%e 0\n",sampleInputs[sInd*nInputs+iplot1]);
                  else
                     fprintf(fp, "%e 1\n",sampleInputs[sInd*nInputs+iplot1]);
               }
               fprintf(fp, "];\n");
               fprintf(fp, "iset = find(X(:,2)==0);\n");
               fprintf(fp, "plot(iset,X(iset,1),'rX','MarkerSize',13)\n");
               fprintf(fp, "hold on\n");
               fprintf(fp, "iset = find(X(:,2)==1);\n");
               fprintf(fp, "plot(iset,X(iset,1),'b*','MarkerSize',13)\n");
               fprintf(fp, "hold off\n");
               fprintf(fp, "axis([0 %d min(X(:,1)) max(X(:,1))])\n", nSamples);
               fwritePlotYLabel(fp, inputNames[iplot1]);
               fwritePlotXLabel(fp, "Sample Number");
               fwritePlotAxes(fp);
               fwritePlotTitle(fp, "1D Input Data Plot");
               fclose(fp);    
               printf("matlabiplt1.m is now ready for input scatter plots.\n");
            }
         }
      }

      // 2D input scatter plot
      // ----------------------------------------------------------------
      else if (!strcmp(command, "iplot2") && nInputs <= 0)
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
         printf("ERROR: no data to analyze (load data first).\n");
      }
      else if (!strcmp(command, "iplot2") && nInputs < 2)
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
         printf("ERROR: iplot2 requires 2 or more inputs.\n");
      }
      else if (!strcmp(command, "iplot2") && nInputs >= 2)
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
         if (psuadeIO_ == NULL) printf("ERROR: data not loaded yet.\n");
         else
         {
            iplot1 = iplot2 = -1;
            sprintf(pString, "Enter the input for x axis (1 - %d) : ", nInputs);
            iplot1 = getInt(1, nInputs, pString);
            iplot1--;
            if (nInputs == 2)
            {
               if (iplot1 == 0) iplot2 = 1;
               else             iplot2 = 0;
            }
            while (iplot2 < 0 || iplot2 >= nInputs || iplot1 == iplot2)
            {
               sprintf(pString, "Enter the input for y axis (1 - %d), not %d : ",
                       nInputs, iplot1+1);
               iplot2 = getInt(1, nInputs, pString);
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
               for (sInd = 0; sInd < nSamples; sInd++)
               {
                  if (sampleOutputs[sInd*nOutputs] > 0.5*PSUADE_UNDEFINED)
                     fprintf(fp, "%e %e 0\n",sampleInputs[sInd*nInputs+iplot1],
                             sampleInputs[sInd*nInputs+iplot2]);
                  else
                     fprintf(fp, "%e %e 1\n",sampleInputs[sInd*nInputs+iplot1],
                             sampleInputs[sInd*nInputs+iplot2]);
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
                       iLowerB[iplot1], iLowerB[iplot2],
                       iUpperB[iplot1], iUpperB[iplot2]);
               fwritePlotXLabel(fp, inputNames[iplot1]);
               fwritePlotYLabel(fp, inputNames[iplot2]);
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
               for (sInd = 0; sInd < nSamples; sInd++)
               {
                  if (sampleOutputs[sInd*nOutputs] > 0.5*PSUADE_UNDEFINED ||
                      sampleStates[sInd] != 1)
                     fprintf(fp, "%24.16e %24.16e 0\n",
                             sampleInputs[sInd*nInputs+iplot1],
                             sampleInputs[sInd*nInputs+iplot2]);
                  else
                     fprintf(fp, "%24.16e %24.16e 1\n",
                             sampleInputs[sInd*nInputs+iplot1],
                             sampleInputs[sInd*nInputs+iplot2]);
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
               fprintf(fp,"axis([%e %e %e %e])\n", iLowerB[iplot1], iUpperB[iplot1],
                       iLowerB[iplot2], iUpperB[iplot2]);
               fwritePlotXLabel(fp, inputNames[iplot1]);
               fwritePlotYLabel(fp, inputNames[iplot2]);
               fwritePlotAxes(fp);
               fwritePlotTitle(fp, "2D Input Data Plot");
               fprintf(fp,"disp('red X: failed runs.')\n");
               fclose(fp);    
               printf("matlabiplt2.m is now ready for input scatter plots.\n");
            }
         }
      }

      // 3D input/output scatter plot
      // ----------------------------------------------------------------
      else if (!strcmp(command, "splot2") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("splot2: scatter plot 2-inputs/1-output in 3D\n");
            printf("syntax: splot2 (no argument needed)\n");
            continue;
         }
         printf("ERROR: no data to analyze (load data first).\n");
      }
      else if (!strcmp(command, "splot2") && nInputs < 2)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("splot2: scatter plot 2-inputs/1-output in 3D\n");
            printf("syntax: splot2 (no argument needed)\n");
            continue;
         }
         printf("ERROR: splot2 requires 2 or more inputs.\n");
      }
      else if (!strcmp(command, "splot2") && nInputs >= 2)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("splot2: scatter plot 2-inputs/1-output in 3D\n");
            printf("syntax: splot2 (no argument needed)\n");
            continue;
         }
         if (psuadeIO_ == NULL) printf("ERROR: data not loaded yet.\n");
         else
         {
            if (nInputs == 2)
            {
               iplot1 = 0;
               iplot2 = 1;
            }
            else
            {
               iplot1 = iplot2 = -1;
               sprintf(pString, "X-axis input ? (1 - %d) ",
                       nInputs);
               iplot1 = getInt(1, nInputs, pString);
               iplot1--;
               iplot2 = iplot1;
               while (iplot1 == iplot2)
               {
                  sprintf(pString, "Y-axis input ? (1 - %d, not %d) ",
                          nInputs, iplot1+1);
                  iplot2 = getInt(1, nInputs, pString);
                  iplot2--;
                  if (iplot1 == iplot2)
                     printf("ERROR: duplicate input number %d.\n",iplot2+1);
               }
            }
            if (nOutputs == 1) oplot1 = 0;
            else
            {
               sprintf(pString, "Z-axis output ? (1 - %d) : ",nOutputs);
               oplot1 = getInt(1, nOutputs, pString);
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
               for (sInd = 0; sInd < nSamples; sInd++)
               {
                  dtemp = sampleOutputs[sInd*nOutputs+oplot1];
                  if (dtemp > 0.5*PSUADE_UNDEFINED || sampleStates[sInd] != 1)
                     fprintf(fp, "%e %e %e 0\n",
                             sampleInputs[sInd*nInputs+iplot1],
                             sampleInputs[sInd*nInputs+iplot2],
                             sampleOutputs[sInd*nOutputs+oplot1]);
                  else
                     fprintf(fp, "%e %e %e 1\n",
                             sampleInputs[sInd*nInputs+iplot1],
                             sampleInputs[sInd*nInputs+iplot2],
                             sampleOutputs[sInd*nOutputs+oplot1]);
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
               fwritePlotXLabel(fp, inputNames[iplot1]);
               fwritePlotYLabel(fp, inputNames[iplot2]);
               fwritePlotZLabel(fp, outputNames[oplot1]);
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
               for (sInd = 0; sInd < nSamples; sInd++)
               {
                  dtemp = sampleOutputs[sInd*nOutputs+oplot1];
                  if (dtemp > 0.5*PSUADE_UNDEFINED || sampleStates[sInd] != 1)
                     fprintf(fp, "%e %e %e 0\n",
                             sampleInputs[sInd*nInputs+iplot1],
                             sampleInputs[sInd*nInputs+iplot2],
                             sampleOutputs[sInd*nOutputs+oplot1]);
                  else
                     fprintf(fp, "%e %e %e 1\n",
                             sampleInputs[sInd*nInputs+iplot1],
                             sampleInputs[sInd*nInputs+iplot2],
                             sampleOutputs[sInd*nOutputs+oplot1]);
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
               fwritePlotXLabel(fp, inputNames[iplot1]);
               fwritePlotYLabel(fp, inputNames[iplot2]);
               fwritePlotZLabel(fp, outputNames[oplot1]);
               fwritePlotAxes(fp);
               fwritePlotTitle(fp, "3D 2-Input/1-Output Data Plot");
               fprintf(fp,"disp('red  *: failed runs.')\n");
               fprintf(fp,"disp('blue X: good   runs.')\n"); 
               fclose(fp);    
               printf("matlabsp2.m is now ready for input scatter plots.\n");
            }
         }
      }

      // 3D input scatter plot
      // ----------------------------------------------------------------
      else if (!strcmp(command, "iplot3") && nInputs <= 0)
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
         printf("ERROR: no data to analyze (load data first).\n");
      }
      else if (!strcmp(command, "iplot3") && nInputs < 3)
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
         printf("ERROR: iplot3 requires 3 or more inputs.\n");
      }
      else if (!strcmp(command, "iplot3") && nInputs >= 3)
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
         if (psuadeIO_ == NULL) printf("ERROR: data not loaded yet.\n");
         else
         {
            iplot1 = iplot2 = iplot3 = -1;
            sprintf(pString, "Enter the input for x axis (1 - %d) : ", nInputs);
            iplot1 = getInt(1, nInputs, pString);
            iplot1--;
            iplot2 = iplot1;
            while (iplot1 == iplot2)
            {
               sprintf(pString,"Enter the input for y axis (1 - %d), not %d : ",
                       nInputs, iplot1+1);
               iplot2 = getInt(1, nInputs, pString);
               iplot2--;
               if (iplot1 == iplot2)
                  printf("ERROR: duplicate input number %d.\n",iplot2+1);
            }
            if (nInputs == 3) iplot3 = 3 - iplot1 - iplot2;
            while (iplot3 < 0 || iplot3 == iplot1 || iplot3 == iplot2)
            {
               sprintf(pString,
                       "Enter the input for z axis (1 - %d), not %d nor %d: ",
                       nInputs, iplot1+1, iplot2+1);
               iplot3 = getInt(1, nInputs, pString);
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
               for (sInd = 0; sInd < nSamples; sInd++)
               {
                  dtemp = sampleOutputs[sInd*nOutputs];
                  if (dtemp > 0.5*PSUADE_UNDEFINED || sampleStates[sInd] != 1)
                     fprintf(fp, "%e %e %e 0\n",
                             sampleInputs[sInd*nInputs+iplot1],
                             sampleInputs[sInd*nInputs+iplot2],
                             sampleInputs[sInd*nInputs+iplot3]);
                  else
                     fprintf(fp, "%e %e %e 1\n",
                             sampleInputs[sInd*nInputs+iplot1],
                             sampleInputs[sInd*nInputs+iplot2],
                             sampleInputs[sInd*nInputs+iplot3]);
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
               fwritePlotXLabel(fp, inputNames[iplot1]);
               fwritePlotYLabel(fp, inputNames[iplot2]);
               fwritePlotZLabel(fp, inputNames[iplot3]);
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
               fprintf(fp, "%% ranflag: to distinguish identical points\n");
               fprintf(fp, "%%     by adding a small perturbation (when on)\n");
               fprintf(fp, "ranflag  = 0;\n");
               fprintf(fp, "%% set cvFlag to 1 to use convex hull\n");
               fprintf(fp, "cvFlag = 0;\n");
               fwritePlotCLF(fp);
               fprintf(fp, "X = [\n");
               for (sInd = 0; sInd < nSamples; sInd++)
               {
                  dtemp = sampleOutputs[sInd*nOutputs];
                  if (dtemp > 0.5*PSUADE_UNDEFINED || sampleStates[sInd] != 1)
                     fprintf(fp, "%e %e %e 0\n",
                             sampleInputs[sInd*nInputs+iplot1],
                             sampleInputs[sInd*nInputs+iplot2],
                             sampleInputs[sInd*nInputs+iplot3]);
                  else
                     fprintf(fp, "%e %e %e 1\n",
                             sampleInputs[sInd*nInputs+iplot1],
                             sampleInputs[sInd*nInputs+iplot2],
                             sampleInputs[sInd*nInputs+iplot3]);
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
               fwritePlotXLabel(fp, inputNames[iplot1]);
               fwritePlotYLabel(fp, inputNames[iplot2]);
               fwritePlotZLabel(fp, inputNames[iplot3]);
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
               fwritePlotXLabel(fp, inputNames[iplot1]);
               fwritePlotYLabel(fp, inputNames[iplot2]);
               fwritePlotZLabel(fp, inputNames[iplot3]);
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
      }

      // 4D input scatter plot
      // ----------------------------------------------------------------
      else if (!strcmp(command, "iplot4m") && nInputs <= 0)
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
         printf("ERROR: no data to analyze (load data first).\n");
      }
      else if (!strcmp(command, "iplot4m") && nInputs < 4)
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
         printf("ERROR: iplot4m requires 4 or more inputs.\n");
      }
      else if (!strcmp(command, "iplot4m") && nInputs >= 4)
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
         if (psuadeIO_ == NULL) printf("ERROR: data not loaded yet.\n");
         else
         {
            if (psPlotTool_ == 1)
            {
               printf("INFO: iplot4m is currently not available for scilab.\n");
               continue;
            }
            iplot1 = iplot2 = iplot3 = iplot4 = -1;
            sprintf(pString, "Enter the input for x axis (1 - %d) : ", nInputs);
            iplot1 = getInt(1, nInputs, pString);
            iplot1--;
            iplot2 = iplot1;
            while (iplot1 == iplot2)
            {
               sprintf(pString,"Enter the input for y axis (1 - %d), not %d : ",
                       nInputs, iplot1+1);
               iplot2 = getInt(1, nInputs, pString);
               iplot2--;
               if (iplot1 == iplot2)
                  printf("ERROR: duplicate input number %d.\n",iplot2+1);
            }
            iplot3 = -1;
            while (iplot3 < 0 || iplot3 == iplot1 || iplot3 == iplot2)
            {
               sprintf(pString,
                       "Enter the input for z axis (1 - %d), not %d nor %d: ",
                       nInputs, iplot1+1, iplot2+1);
               iplot3 = getInt(1, nInputs, pString);
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
                       nInputs, iplot1+1, iplot2+1, iplot3+1);
               iplot4 = getInt(1, nInputs, pString);
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
            fprintf(fp, "%% set cvFlag to 1 to use convex hull\n");
            fprintf(fp, "cvFlag = 0;\n");
            fprintf(fp, "%% change npart to change resolution of the \n");
            fprintf(fp, "%% 4th dimension\n");
            fwritePlotCLF(fp);
            fprintf(fp, "X = [\n");
            for (sInd = 0; sInd < nSamples; sInd++)
            {
               dtemp = sampleOutputs[sInd*nOutputs];
               if (dtemp > 0.5*PSUADE_UNDEFINED || sampleStates[sInd] != 1)
                  fprintf(fp, "%e %e %e %e 0\n",
                          sampleInputs[sInd*nInputs+iplot1],
                          sampleInputs[sInd*nInputs+iplot2],
                          sampleInputs[sInd*nInputs+iplot3],
                          sampleInputs[sInd*nInputs+iplot4]);
               else
                  fprintf(fp, "%e %e %e %e 1\n",
                          sampleInputs[sInd*nInputs+iplot1],
                          sampleInputs[sInd*nInputs+iplot2],
                          sampleInputs[sInd*nInputs+iplot3],
                          sampleInputs[sInd*nInputs+iplot4]);
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
            fwritePlotXLabel(fp, inputNames[iplot1]);
            fwritePlotYLabel(fp, inputNames[iplot2]);
            fwritePlotZLabel(fp, inputNames[iplot3]);
            fwritePlotAxes(fp);
            sprintf(pString,"4D Input Scatter Plot for %s)",inputNames[iplot4]);
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
      }

      // 2D output scatter plot
      // ----------------------------------------------------------------
      else if (!strcmp(command, "oplot2") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("oplot2: plot the sample outputs in two-parameter space.\n");
            printf("syntax: oplot2 (no argument needed)\n");
            printf("This command is useful for examining output correlation\n");
            continue;
         }
         printf("ERROR: no data to analyze (load data first).\n");
      }
      else if (!strcmp(command, "oplot2") && nOutputs < 2)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("oplot2: plot the sample outputs in two-parameter space.\n");
            printf("syntax: oplot2 (no argument needed)\n");
            printf("This command is useful for examining output correlation\n");
            continue;
         }
         printf("ERROR: iplot2 requires 2 or more outputs.\n");
      }
      else if (!strcmp(command, "oplot2") && nOutputs >= 2)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("oplot2: plot the sample outputs in two-parameter space.\n");
            printf("syntax: oplot2 (no argument needed)\n");
            printf("This command is useful for examining output correlation\n");
            continue;
         }
         if (psuadeIO_ == NULL) printf("ERROR: data not loaded yet.\n");
         else
         {
            oplot1 = oplot2 = -1;
            sprintf(pString,"Enter the output for x axis (1 - %d) : ",nOutputs);
            oplot1 = getInt(1, nOutputs, pString);
            oplot1--;
            if (nOutputs == 2)
            {
               if (oplot1 == 0) oplot2 = 1;
               else             oplot2 = 0;
            }
            while (oplot2 < 0 || oplot2 >= nOutputs || oplot1 == oplot2)
            {
               sprintf(pString,"Enter the output for y axis (1 - %d, not %d): ",
                       nOutputs, oplot1+1);
               oplot2 = getInt(1, nOutputs, pString);
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
               for (sInd = 0; sInd < nSamples; sInd++)
               {
                  if ((sampleOutputs[sInd*nOutputs+oplot1]>0.5*PSUADE_UNDEFINED)||
                      (sampleOutputs[sInd*nOutputs+oplot2]>0.5*PSUADE_UNDEFINED))
                     fprintf(fp, "%e %e 0\n",
                             sampleOutputs[sInd*nOutputs+oplot1],
                             sampleOutputs[sInd*nOutputs+oplot2]);
                  else
                     fprintf(fp, "%e %e 1\n",
                             sampleOutputs[sInd*nOutputs+oplot1],
                             sampleOutputs[sInd*nOutputs+oplot2]);
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
               fwritePlotXLabel(fp, outputNames[oplot1]);
               fwritePlotYLabel(fp, outputNames[oplot2]);
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
               for (sInd = 0; sInd < nSamples; sInd++)
               {
                  if ((sampleOutputs[sInd*nOutputs+oplot1]>0.5*PSUADE_UNDEFINED)||
                      (sampleOutputs[sInd*nOutputs+oplot2]>0.5*PSUADE_UNDEFINED))
                     fprintf(fp, "%e %e 0\n",
                             sampleOutputs[sInd*nOutputs+oplot1],
                             sampleOutputs[sInd*nOutputs+oplot2]);
                  else
                     fprintf(fp, "%e %e 1\n",
                             sampleOutputs[sInd*nOutputs+oplot1],
                             sampleOutputs[sInd*nOutputs+oplot2]);
               }
               fprintf(fp, "];\n");
               fprintf(fp, "iset = find(Y(:,3) == 1);\n");
               fprintf(fp, "plot(Y(iset,1),Y(iset,2),\'bX\')\n");
               fprintf(fp, "hold on\n");
               fprintf(fp, "iset = find(Y(:,3) == 0);\n");
               fprintf(fp, "plot(Y(iset,1),Y(iset,2),\'r*\')\n");
               fprintf(fp, "hold off\n");
               fwritePlotXLabel(fp, outputNames[oplot1]);
               fwritePlotYLabel(fp, outputNames[oplot2]);
               fwritePlotAxes(fp);
               fwritePlotTitle(fp, "2D Output Scatter Plot");
               fclose(fp);    
               printf("matlaboplt2.m is now ready for scatter plots.\n");
            }
         }
      }

      // response surface based MCMC
      // ----------------------------------------------------------------
      else if (!strcmp(command, "rsmcmc") && nInputs <= 0)
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
         printf("ERROR: no data to analyze (load data first).\n");
      }
      else if (!strcmp(command, "rsmcmc") && nInputs > 0)
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

      // simulation based MCMC (not done yet)
      else if (!strcmp(command, "mcmc") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("mcmc: perform MCMC on actual simulator\n");
            printf("syntax: mcmc (no argument needed)\n");
            continue;
         }
         printf("ERROR: no data to analyze (load data first).\n");
      }
      else if (!strcmp(command, "mcmc") && nInputs > 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("mcmc: perform MCMC on actual simulator\n");
            printf("syntax: mcmc (no argument needed)\n");
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

      // ================================================================
      // miscellaneous commands
      // ================================================================
      // sample refinement
      else if (!strcmp(command, "refine"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("refine: refine a sample (loaded previously using 'load')\n");
            printf("syntax: refine (no argument needed)\n");
            continue;
         }
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no data (load data first).\n");
         }
         else
         {
            psuadeIO_->getParameter("method_sampling", pPtr);
            samplingMethod = pPtr.intData_;
            kk = psSamExpertMode_;
            psSamExpertMode_ = 0;
            sampPtr = (Sampling *) SamplingCreateFromID(samplingMethod);
            sampPtr->setPrintLevel(0);
            sampPtr->setInputBounds(nInputs, iLowerB, iUpperB);
            sampPtr->setInputParams(nInputs, NULL, NULL, NULL);
            sampPtr->setOutputParams(nOutputs);
            psuadeIO_->getParameter("method_nreplications",pPtr);
            nReps = pPtr.intData_;
            sampPtr->setSamplingParams(nSamples, nReps, -1);
            sampPtr->initialize(1);
            sampPtr->loadSamples(nSamples, nInputs, nOutputs, sampleInputs,
                                 sampleOutputs, sampleStates);
            status = sampPtr->refine(2, 0, 0.0, nSamples, NULL);

            if (status == 0)
            {
               if (sampleInputs  != NULL) delete [] sampleInputs;
               if (sampleOutputs != NULL) delete [] sampleOutputs;
               if (sampleStates  != NULL) delete [] sampleStates;
               nSamples = sampPtr->getNumSamples();
               sampleInputs  = new double[nInputs*nSamples];
               sampleOutputs = new double[nOutputs*nSamples];
               sampleStates  = new int[nSamples];
               sampPtr->getSamples(nSamples, nInputs, nOutputs, sampleInputs,
                                   sampleOutputs, sampleStates);
               psuadeIO_->updateInputSection(nSamples,nInputs,NULL,NULL,NULL,
                                             sampleInputs,NULL); 
               psuadeIO_->updateOutputSection(nSamples,nOutputs,sampleOutputs, 
                                              sampleStates,NULL); 
               psuadeIO_->updateMethodSection(-1, nSamples, -1, -1, -1);
               printf("refine successful: use write to store data set.\n");
            }
            else
            {
               printf("ERROR: refine not successful.\n");
            }
            delete sampPtr;
            sampPtr = NULL;
            psSamExpertMode_ = kk;
         }
      }

      // adaptive sample refinement
      // ----------------------------------------------------------------
      else if (!strcmp(command, "a_refine"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("a_refine: adaptively refine a loaded sample\n");
            printf("syntax: a_refine (no argument needed)\n");
            continue;
         }
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no data (load data first).\n");
         }
         outputID = 0;
         if (nOutputs > 1)
         {
            sprintf(pString,"Select output to use for a_refine (1 - %d) ",nOutputs);
            outputID = getInt(1, nOutputs, pString);
         }
         psuadeIO_->getParameter("method_sampling", pPtr);
         samplingMethod = pPtr.intData_;
         if (samplingMethod != PSUADE_SAMP_METIS)
         {
            printf("ERROR: adaptive refinement requires METIS sample.\n");
            continue;
         }
         flag = 1;
         for (ss = 0; ss < nSamples*nOutputs; ss++)
         {
            if (sampleOutputs[ss] != 1 && sampleOutputs[ss] != 0)
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
            int origNSamples = getInt(1, nSamples, pString);
            sprintf(pString,"How many sample points to add? (1 - %d) ",nSamples);
            int refineSize = getInt(1, nSamples, pString);
            int tempSamExpert = psSamExpertMode_;
            psSamExpertMode_ = 0;
            sampPtr = (Sampling *) SamplingCreateFromID(samplingMethod);
            sampPtr->setPrintLevel(0);
            sampPtr->setInputBounds(nInputs, iLowerB, iUpperB);
            sampPtr->setInputParams(nInputs, NULL, NULL, NULL);
            sampPtr->setOutputParams(nOutputs);
            sampPtr->setSamplingParams(nSamples, -1, -1);
            sampPtr->initialize(1);
            sampPtr->loadSamples(nSamples, nInputs, nOutputs, sampleInputs,
                                 sampleOutputs, sampleStates);
            sampAux = (Sampling *) SamplingCreateFromID(samplingMethod);
            sampAux->setPrintLevel(0);
            sampAux->setInputBounds(nInputs, iLowerB, iUpperB);
            sampAux->setInputParams(nInputs, NULL, NULL, NULL);
            sampAux->setOutputParams(nOutputs);
            sampAux->setSamplingParams(nSamples, -1, -1);
            sampAux->initialize(1);
            if (nOutputs > 1)
            {
               tempY = new double[nSamples];
               for (ss = 0; ss < nSamples; ss++)
                  tempY[ss] = sampleOutputs[ss*nOutputs+outputID];
               sampAux->loadSamples(nSamples, nInputs, 1, sampleInputs,
                                    tempY, sampleStates);
            }
            else
            {
               sampAux->loadSamples(nSamples, nInputs, nOutputs, sampleInputs,
                                    sampleOutputs, sampleStates);
            }
            sparam.clear();
            sparam.append("setUniformRefinement");
            sampAux->setParam(sparam);
            status = sampAux->refine(2, 1, 0.0, nSamples, NULL);
            int nSamples2 = sampAux->getNumSamples();
            if (nSamples2 != 2 * nSamples)
            {
               printf("a_refine ERROR: something is wrong. Consult developers.\n");
               printf("                refined sample size != 2 * original size\n");
               printf("                May be due to too many levels of refinement.\n");
               delete sampAux;;
               sampAux = NULL;
               delete sampPtr;;
               sampPtr = NULL;
               if (nOutputs > 1) delete [] tempY;
               tempY = NULL;
               continue;
            }
            double *samInputs2  = new double[nSamples2*nInputs];
            double *samOutputs2 = new double[nSamples2];
            int    *samStates2  = new int[nSamples2];
            double *samStds2    = new double[nSamples2];
            sampAux->getSamples(nSamples2, nInputs, 1, samInputs2, 
                                samOutputs2, samStates2);
            faType = PSUADE_RS_MARSB;
            faPtr = genFA(faType, nInputs, iOne, nSamples);
            faPtr->setNPtsPerDim(32);
            faPtr->setBounds(iLowerB, iUpperB);
            faPtr->setOutputLevel(0);
            int numMars = 100, ivar1;
            double **marsX = new double*[numMars];
            double **marsY = new double*[numMars];
            for (ii = 0; ii < numMars; ii++)
            {
               marsX[ii] = new double[nInputs*nSamples];
               marsY[ii] = new double[nSamples];
               for (ss = 0; ss < nSamples; ss++)
               {
                  if (ss < origNSamples)
                       ivar1 = PSUADE_rand() % origNSamples;
                  else ivar1 = ss;
                  for (jj = 0; jj < nInputs; jj++)
                     marsX[ii][ss*nInputs+jj] =
                         sampleInputs[ivar1*nInputs+jj];
                  marsY[ii][ss] = sampleOutputs[ivar1*nOutputs+outputID];
               }
            }
            strcpy(cString, "mars_params");
            int ivar2 = 2 * nInputs / 3 + 1;
            targv[0] = (char *) cString;
            targv[1] = (char *) &nSamples;
            targv[2] = (char *) &ivar2;
            faPtr->setParams(3, targv);
            strcpy(cString, "num_mars");
            targv[0] = (char *) cString;
            targv[1] = (char *) &numMars;
            faPtr->setParams(2, targv);
            strcpy(cString, "mars_sample");
            targv[0] = (char *) cString;
            targv[2] = (char *) &nSamples;
            for (ii = 0; ii < numMars; ii++)
            {
               targv[1] = (char *) &ii;
               targv[3] = (char *) marsX[ii];
               targv[4] = (char *) marsY[ii];
               faPtr->setParams(5, targv);
            }
            count = -999;
            if (nOutputs > 1)
                 faPtr->genNDGridData(sampleInputs,sampleOutputs,&count,NULL,NULL);
            else faPtr->genNDGridData(sampleInputs,tempY,&count,NULL,NULL);

            faPtr->evaluatePointFuzzy(nSamples2-nSamples, 
                            &samInputs2[nInputs*nSamples], 
                            &samOutputs2[nSamples], &samStds2[nSamples]);
            for (ss = 0; ss < nSamples; ss++) 
               samStds2[ss] = PABS(samStds2[ss+nSamples]);
            if (outputLevel_ > 4)
            {
               printf("Standard deviations to be used to select refinements.\n");
               for (ss = 0; ss < nSamples; ss++) 
                  printf("Sample point %7d: stdev = %e\n", ss+1, samStds2[ss]);
            }

            sparam.clear();
            sparam.append("setAdaptiveRefinementBasedOnErrors");
            sampPtr->setParam(sparam);
            sparam.clear();
            sprintf(cString, "setRefineSize %d", refineSize);
            sparam.append(cString);
            sampPtr->setParam(sparam);
            sampPtr->refine(2,1,1.0e-6,nSamples,samStds2);
            if (sampleInputs  != NULL) delete [] sampleInputs;
            if (sampleOutputs != NULL) delete [] sampleOutputs;
            if (sampleStates  != NULL) delete [] sampleStates;

            nSamples = sampPtr->getNumSamples();
            sampleInputs  = new double[nInputs*nSamples];
            sampleOutputs = new double[nOutputs*nSamples];
            sampleStates  = new int[nSamples];
            sampPtr->getSamples(nSamples, nInputs, nOutputs, sampleInputs,
                                sampleOutputs, sampleStates);
            psuadeIO_->updateInputSection(nSamples,nInputs,NULL,NULL,NULL,
                                          sampleInputs,NULL); 
            psuadeIO_->updateOutputSection(nSamples,nOutputs,sampleOutputs, 
                                           sampleStates,NULL); 
            psuadeIO_->updateMethodSection(-1, nSamples, -1, -1, -1);
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
            sampPtr->setPrintLevel(0);
            sampPtr->setInputBounds(nInputs, iLowerB, iUpperB);
            sampPtr->setInputParams(nInputs, NULL, NULL, NULL);
            sampPtr->setOutputParams(nOutputs);
            sampPtr->setSamplingParams(nSamples, -1, -1);
            sampPtr->initialize(1);
            if (nOutputs > 1)
            {
               tempY = new double[nSamples];
               for (ss = 0; ss < nSamples; ss++)
                  tempY[ss] = sampleOutputs[ss*nOutputs+outputID];
               sampPtr->loadSamples(nSamples, nInputs, 1, sampleInputs,
                                    tempY, sampleStates);
            }
            else
            {
               sampPtr->loadSamples(nSamples, nInputs, nOutputs, sampleInputs,
                                    sampleOutputs, sampleStates);
            }
            sparam.clear();
            sparam.append("setAdaptiveRefinementBasedOnOutputs");
            sampPtr->setParam(sparam);
            sparam.clear();
            sprintf(cString, "setRefineSize %d", nSamples);
            sparam.append(cString);
            sampPtr->setParam(sparam);
            sampPtr->refine(2,1,1.0e-6,0,NULL);
            if (sampleInputs  != NULL) delete [] sampleInputs;
            if (sampleStates  != NULL) delete [] sampleStates;
            tempV = sampleOutputs;

            kk = nSamples;
            nSamples = sampPtr->getNumSamples();
            sampleInputs  = new double[nInputs*nSamples];
            sampleOutputs = new double[nOutputs*nSamples];
            sampleStates  = new int[nSamples];
            if (nOutputs > 1)
            {
               sampPtr->getSamples(nSamples, nInputs, 1, sampleInputs,
                                   tempY, sampleStates);
               for (ss = 0; ss < kk*nOutputs; ii++)
                  sampleOutputs[ss] = tempV[ss];
               for (ss = kk*nOutputs; ss < nSamples*nOutputs; ii++)
                  sampleOutputs[ss] = PSUADE_UNDEFINED;
            }
            else
            {
               sampPtr->getSamples(nSamples, nInputs, nOutputs, sampleInputs,
                                   sampleOutputs, sampleStates);
            }
            psuadeIO_->updateInputSection(nSamples,nInputs,NULL,NULL,NULL,
                                          sampleInputs,NULL); 
            psuadeIO_->updateOutputSection(nSamples,nOutputs,sampleOutputs, 
                                           sampleStates,NULL); 
            psuadeIO_->updateMethodSection(-1, nSamples, -1, -1, -1);
            printf("refine successful: use write to store the refined sample.\n");
            printf("                   Then run the newly created sample points\n");
            printf("                   in this refined sample.\n");
            delete tempV;
            delete sampPtr;
            sampPtr = NULL;
            if (sampleOutputs != NULL) delete [] sampleOutputs;
            psSamExpertMode_ = tempSamExpert;;
         }
      }

      // sample quality check 
      // ----------------------------------------------------------------
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
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no data (load data first).\n");
         }
         else
         {
            dtemp = SamplingQuality(nSamples, nInputs, sampleInputs);
            printf("Max-min distance between sample points = %e\n", dtemp);
         }
      }

      // validate certain sample outputs (a range)
      // ----------------------------------------------------------------
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
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no data (validate assumes outputs exist.)\n");
         }
         else
         {
            sprintf(pString,"Enter first sample number (1 - %d) = ", 
                    nSamples);
            ind = getInt(1, nSamples, pString);
            ind--;
            sprintf(pString, "Enter last sample number (%d - %d) = ", 
                    ind+1, nSamples);
            ind2 = getInt(ind+1, nSamples, pString);
            psuadeIO_->updateInputSection(nSamples,nInputs,NULL,NULL,NULL,
                                          sampleInputs,NULL); 
            for (ii = ind; ii < ind2; ii++) sampleStates[ii] = 1;
            psuadeIO_->updateOutputSection(nSamples,nOutputs,sampleOutputs, 
                                           sampleStates,NULL); 
            psuadeIO_->updateMethodSection(samplingMethod,nSamples,nReps,
                                           -1,-1);
            psuadeIO_->writePsuadeFile(dataFile,0);
         }
      }

      // invalidate all sample outputs (set to undefined)
      // ----------------------------------------------------------------
      else if (!strcmp(command, "invalidate"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("invalidate: select sample points to be unevaluated\n");
            printf("            based on a certain sample output value\n");
            printf("            (or deselect all sample points).\n");
            printf("syntax: invalidate (no argument needed)\n");
            printf("This command just set the ready flag of selected sample\n");
            printf("points to be 'not ready' based on some an output value.\n");
            printf("If a range of sample points is to be unevaluated irrespective\n");
            printf("of output values, first use 'invalidate' for all sample points,\n");
            printf("and then use 'validate' to restore the range of desired \n");
            printf("sample points.\n");
            continue;
         }
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no data (load data first).\n");
         }
         else
         {
            printf("Invalidate sample points based on the value of an output.\n");
            printf("If a range of sample points is to be unevaluated, first\n");
            printf("use 'invalidate' for all sample points, and then use\n");
            printf("'validate' to restore the range of desired sample points.\n");
            printf("If all sample points are to be invalidated, enter 0.\n");
            sprintf(pString, "Enter output number (1 - %d, 0 for all points) : ", 
                    nOutputs);
            outputID = getInt(0, nOutputs, pString);
            outputID--;
            if (outputID == -1)
            {
               for (ii = 0; ii < nSamples; ii++) sampleStates[ii] = 0;
               psuadeIO_->updateInputSection(nSamples,nInputs,NULL,NULL,NULL,
                                             sampleInputs,NULL); 
               psuadeIO_->updateOutputSection(nSamples,nOutputs,sampleOutputs, 
                                              sampleStates,NULL); 
               psuadeIO_->updateMethodSection(samplingMethod,nSamples,nReps,
                                              -1,-1);
               printf("Use the write command to store the filtered sample.\n");
            }
            else
            {
               Ymax = - PSUADE_UNDEFINED;
               Ymin =   PSUADE_UNDEFINED;
               for (sInd = 0; sInd < nSamples; sInd++)
               {
                  if (sampleOutputs[sInd*nOutputs+outputID] > Ymax)
                     Ymax = sampleOutputs[sInd*nOutputs+outputID];
                  if (sampleOutputs[sInd*nOutputs+outputID] < Ymin)
                     Ymin = sampleOutputs[sInd*nOutputs+outputID];
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
               for (sInd = 0; sInd < nSamples; sInd++)
               {
                  if (sampleOutputs[sInd*nOutputs+outputID] < threshL ||
                      sampleOutputs[sInd*nOutputs+outputID] > threshU)
                     sampleStates[sInd] = 0;
               }
               psuadeIO_->updateInputSection(nSamples,nInputs,NULL,NULL,NULL,
                                             sampleInputs, inputNames); 
               psuadeIO_->updateOutputSection(nSamples, nOutputs, sampleOutputs, 
                                              sampleStates, outputNames);
               printf("INFO: use 'purge' to take out the invalid points.\n");
               printf("INFO: then use 'write' to store the filtered sample.\n");
            }
         }
      }

      // randomize the orders of the resident sample
      // ----------------------------------------------------------------
      else if (!strcmp(command, "srandomize") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("srandomize: randomize the orders of the resident sample.\n");
            printf("syntax: srandomize <n>\n");
            printf("That command can be used, for example, with rstest_cv.\n");
            continue;
         }
         printf("ERROR: no data (load data first).\n");
      }
      else if (!strcmp(command, "srandomize") && nInputs > 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("srandomize: randomize the orders of the resident sample.\n");
            printf("syntax: srandomize <n>\n");
            continue;
         }
         if (psuadeIO_ == NULL) printf("ERROR: data not loaded yet.\n");
         else
         {
            tempX = new double[nInputs * nSamples];
            tempY = new double[nOutputs * nSamples];
            tempI = new int[nSamples];;
            generateRandomIvector(nSamples, tempI);
            for (ii = 0; ii < nSamples; ii++)
            {
               ind = tempI[ii];
               for (jj = 0; jj < nInputs; jj++)
                  tempX[ii*nInputs+jj] = sampleInputs[ind*nInputs+jj];
               for (jj = 0; jj < nOutputs; jj++)
                  tempY[ii*nOutputs+jj] = sampleOutputs[ind*nOutputs+jj];
               tempI[ii] = sampleStates[ind];
            }
            delete [] sampleInputs;
            delete [] sampleOutputs;
            delete [] sampleStates;
            sampleInputs = tempX;
            sampleOutputs = tempY;
            sampleStates = tempI;
         }
      }

      // take out failed sample points
      // ----------------------------------------------------------------
      else if (!strcmp(command, "purge") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("purge: take out bad sample points (Output=UNDEFINED or status!=1)\n");
            printf("       (It also has an option to take out repeated sample points.)\n");
            printf("syntax: purge (no argument needed)\n");
            continue;
         }
         printf("ERROR: no data (load data first).\n");
      }
      else if (!strcmp(command, "purge") && nInputs > 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("purge: take out bad sample points (Output=UNDEFINED or status!=1)\n");
            printf("       (It also has an option to take out repeated sample points.)\n");
            printf("syntax: purge (no argument needed)\n");
            continue;
         }
         if (psuadeIO_ == NULL) printf("ERROR: data not loaded yet.\n");
         else
         {
            sprintf(pString,
                    "Do you want repeated sample points purged too? (y or n)");
            getString(pString, winput);
            setCompare = 0;
            if (winput[0] == 'y') setCompare = 1;
            kk = 0;
            for (sInd = 0; sInd < nSamples; sInd++)
            {
               for (oInd = 0; oInd < nOutputs; oInd++)
                  if (sampleOutputs[sInd*nOutputs+oInd] == 
                      PSUADE_UNDEFINED || sampleStates[sInd] == 0) break; 
               if (oInd == nOutputs)
               {
                  if (setCompare == 1)
                     jj = compareSamples(sInd,nSamples,nInputs,
                                         sampleInputs, sampleStates);
                  else jj = -1;
                  if (jj < 0 || jj > sInd)
                  {
                     for (iInd = 0; iInd < nInputs; iInd++)
                        sampleInputs[kk*nInputs+iInd] = 
                           sampleInputs[sInd*nInputs+iInd];
                     for (oInd = 0; oInd < nOutputs; oInd++)
                        sampleOutputs[kk*nOutputs+oInd] = 
                           sampleOutputs[sInd*nOutputs+oInd];
                     sampleStates[kk] = sampleStates[sInd]; 
                     kk++;
                  }
               }
            }
            nSamples = kk;
            samplingMethod = PSUADE_SAMP_MC;
            nReps = 1;
            psuadeIO_->updateMethodSection(samplingMethod,nSamples,nReps,
                                           0,-1);
            psuadeIO_->updateInputSection(nSamples, nInputs, NULL, iLowerB,
                                      iUpperB, sampleInputs, NULL);
            psuadeIO_->updateOutputSection(nSamples, nOutputs, sampleOutputs, 
                                       sampleStates, outputNames);
            printf("Number of sample points after purge = %d\n", nSamples);
            printf("purge completed. Use write to store the reduced sample.\n");
         }
      }

      // remove duplicate points 
      // ----------------------------------------------------------------
      else if (!strcmp(command, "rm_duplicates") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rm_duplicates: take out duplicate sample points\n");
            printf("syntax: rm_duplicates (no argument needed)\n");
            continue;
         }
         printf("ERROR: no data (load data first).\n");
      }
      else if (!strcmp(command, "rm_duplicates") && nInputs > 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rm_duplicates: take out duplicate sample points\n");
            printf("syntax: rm_duplicates (no argument needed)\n");
            continue;
         }
         if (psuadeIO_ == NULL) printf("ERROR: data not loaded yet.\n");
         else
         {
            kk = 0;
            for (sInd = 0; sInd < nSamples; sInd++)
            {
               jj = compareSamples(sInd,nSamples,nInputs,sampleInputs,
                                   sampleStates);
               if (jj < 0 || jj > sInd)
               {
                  for (iInd = 0; iInd < nInputs; iInd++)
                     sampleInputs[kk*nInputs+iInd] = 
                           sampleInputs[sInd*nInputs+iInd];
                  for (oInd = 0; oInd < nOutputs; oInd++)
                     sampleOutputs[kk*nOutputs+oInd] = 
                           sampleOutputs[sInd*nOutputs+oInd];
                  sampleStates[kk] = sampleStates[sInd]; 
                  kk++;
               }
            }
            nSamples = kk;
            samplingMethod = PSUADE_SAMP_MC;
            nReps = 1;
            psuadeIO_->updateMethodSection(samplingMethod,nSamples,nReps,
                                           0,-1);
            psuadeIO_->updateInputSection(nSamples, nInputs, NULL, iLowerB,
                                      iUpperB, sampleInputs, NULL);
            psuadeIO_->updateOutputSection(nSamples, nOutputs, sampleOutputs, 
                                       sampleStates, outputNames);
            printf("Number of sample points after rm_repeats = %d\n", nSamples);
            printf("rm_duplicates completed. Use write to store the reduced sample.\n");
         }
      }

      // take out infeasible sample points
      // ----------------------------------------------------------------
      else if (!strcmp(command, "ifilter") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("ifilter: take out sample points based on input constraints\n");
            printf("syntax: ifilter (no argument needed)\n");
            continue;
         }
         printf("ERROR: no data (load data first).\n");
      }
      else if (!strcmp(command, "ifilter") && nInputs > 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("ifilter: take out sample points based on input constraints\n");
            printf("syntax: ifilter (no argument needed)\n");
            continue;
         }
         if (psuadeIO_ == NULL) printf("ERROR: data not loaded yet.\n");
         else
         {
            inputID = 0;
            if (nInputs > 1)
            {
               sprintf(pString, "Enter input number (1 - %d) : ", nInputs);
               inputID = getInt(1, nInputs, pString);
               inputID--;
            }
            Xmin = sampleInputs[inputID];
            for (sInd = 1; sInd < nSamples; sInd++)
               if (sampleInputs[sInd*nInputs+inputID] < Xmin) 
                  Xmin = sampleInputs[sInd*nInputs+inputID];
            Xmax = sampleInputs[inputID];
            for (sInd = 1; sInd < nSamples; sInd++)
               if (sampleInputs[sInd*nInputs+inputID] > Xmax) 
                  Xmax = sampleInputs[sInd*nInputs+inputID];
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
            for (sInd = 0; sInd < nSamples; sInd++)
            {
               if (sampleInputs[sInd*nInputs+inputID] >= threshL && 
                   sampleInputs[sInd*nInputs+inputID] <= threshU) 
               {
                  for (iInd = 0; iInd < nInputs; iInd++)
                     sampleInputs[kk*nInputs+iInd] = 
                        sampleInputs[sInd*nInputs+iInd];
                  for (oInd = 0; oInd < nOutputs; oInd++)
                     sampleOutputs[kk*nOutputs+oInd] = 
                        sampleOutputs[sInd*nOutputs+oInd];
                  sampleStates[kk] = sampleStates[sInd]; 
                  kk++;
               }
            }
            nSamples = kk;
            samplingMethod = PSUADE_SAMP_MC;
            nReps = 1;
            psuadeIO_->updateInputSection(nSamples,nInputs,NULL,NULL,NULL,
                                          sampleInputs,inputNames); 
            psuadeIO_->updateOutputSection(nSamples, nOutputs, sampleOutputs, 
                                           sampleStates, outputNames);
            psuadeIO_->updateMethodSection(PSUADE_SAMP_MC,nSamples,nReps,0,-1);
            printf("ifilter completed. Use write and load again to continue.\n");
         }
      }

      // take out sample points outside bounds
      // ----------------------------------------------------------------
      else if (!strcmp(command, "ofilter") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("ofilter: take out sample points based on output constraints\n");
            printf("syntax: ofilter (no argument needed)\n");
            continue;
         }
         printf("ERROR: no data (load data first).\n");
      }
      else if (!strcmp(command, "ofilter") && nInputs > 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("ofilter: take out sample points based on output constraints\n");
            printf("syntax: ofilter (no argument needed)\n");
            continue;
         }
         if (psuadeIO_ == NULL) printf("ERROR: data not loaded yet.\n");
         else
         {
            if (nOutputs == 1) outputID = 0;
            else
            {
               sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
               outputID = getInt(1, nOutputs, pString);
               outputID--;
            }
            Ymax = - PSUADE_UNDEFINED;
            Ymin =   PSUADE_UNDEFINED;
            for (sInd = 0; sInd < nSamples; sInd++)
            {
               if (sampleOutputs[sInd*nOutputs+outputID] > Ymax)
                  Ymax = sampleOutputs[sInd*nOutputs+outputID];
               if (sampleOutputs[sInd*nOutputs+outputID] < Ymin)
                  Ymin = sampleOutputs[sInd*nOutputs+outputID];
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
            for (sInd = 0; sInd < nSamples; sInd++)
            {
               if (sampleOutputs[sInd*nOutputs+outputID] >= threshL &&
                   sampleOutputs[sInd*nOutputs+outputID] <= threshU)
               {
                  for (iInd = 0; iInd < nInputs; iInd++)
                     sampleInputs[kk*nInputs+iInd] = 
                           sampleInputs[sInd*nInputs+iInd];
                  for (oInd = 0; oInd < nOutputs; oInd++)
                     sampleOutputs[kk*nOutputs+oInd] = 
                           sampleOutputs[sInd*nOutputs+oInd];
                  sampleStates[kk] = sampleStates[sInd]; 
                  kk++;
               }
            }
            nSamples = kk;
            samplingMethod = PSUADE_SAMP_MC;
            nReps = 1;
            psuadeIO_->updateInputSection(nSamples,nInputs,NULL,NULL,NULL,
                                          sampleInputs,inputNames); 
            psuadeIO_->updateOutputSection(nSamples, nOutputs, sampleOutputs, 
                                           sampleStates, outputNames);
            psuadeIO_->updateMethodSection(PSUADE_SAMP_MC,nSamples,nReps,0,-1);
            printf("ofilter completed. Use write and load again to continue.\n");
         }
      }

      // manipulate outputs
      // ----------------------------------------------------------------
      else if (!strcmp(command, "oop") && nInputs <= 0)
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
         printf("ERROR: no data (load data first).\n");
      }
      else if (!strcmp(command, "oop") && nInputs > 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("oop: form linear combinations of output values for each point.\n");
            printf("syntax: oop (no argument needed)\n");
            continue;
         }
         if (psuadeIO_ == NULL) printf("ERROR: data not loaded yet.\n");
         else
         {
            printf("Set output <1> = <a> * output <2> + <b> * output <3>\n");
            sprintf(pString, "Enter output <1> (1 - %d) : ", nOutputs);
            ii = getInt(1, nOutputs, pString);
            ii--;
            sprintf(pString, "Enter the value <a>: ");
            aVal = getDouble(pString);
            sprintf(pString, "Enter output <2> (1 - %d) : ", nOutputs);
            jj = getInt(1, nOutputs, pString);
            jj--;
            sprintf(pString, "Enter the value <b>: ");
            bVal = getDouble(pString);
            sprintf(pString, "Enter output <3> (1 - %d) : ", nOutputs);
            kk = getInt(1, nOutputs, pString);
            kk--;
            for (sInd = 0; sInd < nSamples; sInd++)
            {
               sampleOutputs[sInd*nOutputs+ii] = 
                           aVal * sampleOutputs[sInd*nOutputs+jj] + 
                           bVal * sampleOutputs[sInd*nOutputs+kk]; 
            }
            psuadeIO_->updateOutputSection(nSamples, nOutputs, sampleOutputs, 
                                           sampleStates, outputNames);
            printf("oop completed. Use write and load again to continue.\n");
         }
      }

      // nearest neighbor analysis (for detecting outliers) 
      // ----------------------------------------------------------------
      else if (!strcmp(command, "nna") && nInputs <= 0)
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
         printf("ERROR: no data to analyze (load data first).\n");
      }
      else if (!strcmp(command, "nna") && nInputs > 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("nna: nearest neighbor analysis (for detecting outliers)\n");
            printf("syntax: nna (no argument needed)\n");
            continue;
         }
         if (psuadeIO_ == NULL) printf("ERROR: data not loaded yet.\n");
         else
         {
            outputID = 0;
            printf("nna: for each data point, find nearest neighbor and");
            printf(" plot changes in the outputs.\n");
            if (nOutputs > 1)
            {
               sprintf(pString, "Enter output number (1 - %d) = ", nOutputs);
               outputID = getInt(1, nOutputs, pString);
               outputID--;
            }
            if (psPlotTool_ == 1)
            {
               fp = fopen("scilabnna.sci", "w");
               if (fp == NULL)
               {
                  printf("ERROR: cannot open file scilabnna.sci.\n");
                  continue;
               }
               fprintf(fp, "// nearest neighbor analysis\n");
               fprintf(fp, "// for detecting outliers.\n");
               fprintf(fp, "// The following plot is : \n");
               fprintf(fp, "// Y-axis: delta output / distance\n");
               fprintf(fp, "// X-axis: sample number\n");
            }
            else
            {
               fp = fopen("matlabnna.m", "w");
               if (fp == NULL)
               {
                  printf("ERROR: cannot open file matlabnna.m.\n");
                  continue;
               }
               fprintf(fp, "%% nearest neighbor analysis\n");
               fprintf(fp, "%% for detecting outliers.\n");
               fprintf(fp, "%% The following plot is : \n");
               fprintf(fp, "%% Y-axis: delta output / distance\n");
               fprintf(fp, "%% X-axis: sample number\n");
            }
            fprintf(fp, "A = [\n");
            for (ii = 0; ii < nSamples; ii++)
            {
               minDist = PSUADE_UNDEFINED;
               ind = -1;
               for (kk = 0; kk < nSamples; kk++)
               {
                  if (ii != kk)
                  {
                     dtemp = 0.0;
                     for (jj = 0; jj < nInputs; jj++)
                     {
                        dtemp = sampleInputs[ii*nInputs+jj] - 
                                sampleInputs[kk*nInputs+jj];
                        dtemp /= (iUpperB[ii] - iLowerB[ii]);
                        dtemp += pow(dtemp, 2.0);
                     }
                     if (dtemp < minDist && dtemp > 0.0)
                     {
                        minDist = dtemp;
                        ind = kk;
                     }
                  }
               }
               dtemp = (sampleOutputs[ii*nOutputs+outputID] -
                        sampleOutputs[ind*nOutputs+outputID])/dtemp;
               fprintf(fp, "%7d %24.16e\n", ii+1, dtemp);
            }
            fprintf(fp, "];\n");
            if (psPlotTool_ == 1)
                 fprintf(fp, "plot(A(:,1),A(:,2),'x');\n");
            else fprintf(fp, "plot(A(:,1),A(:,2),'x','MarkerWeight','bold')\n");
            fwritePlotXLabel(fp, "Sample number");
            fwritePlotYLabel(fp, "Delta output / distance");
            fwritePlotTitle(fp, "Delta test Statistics");
            fwritePlotAxes(fp);
            fclose(fp);
            if (psPlotTool_ == 1)
                 printf("Nearest neighbors analysis is now in file scilabnna.sci\n");
            else printf("Nearest neighbors analysis is not in file matlabnna.m.\n");
         }
      }

      // polynomial regression for interface
      // ----------------------------------------------------------------
      else if (!strcmp(command, "interface_track") && nInputs <= 0)
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
         printf("ERROR: no data to analyze (load data first).\n");
      }
      else if (!strcmp(command, "interface_track") && nInputs > 0)
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
         if (psuadeIO_ == NULL) printf("ERROR: data not loaded yet.\n");
         else
         {
            outputID = 0;
            if (nOutputs > 1)
            {
               sprintf(pString, "Which output to use (1 - %d) = ", nOutputs);
               outputID = getInt(1, nOutputs, pString);
               outputID--;
            }
            int nIndex = 3, ss2;
            if (psAnaExpertMode_ == 1)
            {
               printf("In order to track the interface, the outputs of the nearest\n");
               printf("neighbors of each sample point are compared against themselves.\n");
               printf("You need to select the number of neighbors (K) to examine. The more\n");
               printf("neighbors are used, the more fuzzy the interface will be. On the\n");
               printf("other hand, the less neighbors are used, the sample size may be\n");
               printf("too small to be useful. Try different ones. The default is 3.\n");
               sprintf(pString, "What is K (>= 1, <= 10, default=3)? ");
               nIndex = getInt(1, 10, pString);
            }
            double *distPairs  = new double[nSamples * (nSamples - 1) / 2];
            int    *minIndices = new int[nIndex];;
            for (ss = 0; ss < nSamples; ss++)
            {
               for (ss2 = 0; ss2 < ss; ss2++)
               {
                  kk = ss * (ss - 1) / 2 + ss2;
                  distPairs[kk] = 0.0;
                  for (jj = 0; jj < nInputs; jj++)
                  {
                     ddata = sampleInputs[ss*nInputs+jj] - sampleInputs[ss2*nInputs+jj];
                     ddata = ddata / (iUpperB[jj] - iLowerB[jj]);
                     distPairs[kk] += pow(ddata, 2.0);
                  }
               }
            } 
            count = 0;
            tempX = new double[nSamples*nInputs];
            for (ss = 0; ss < nSamples; ss++)
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
               }
               flag = 0;
               for (jj = 0; jj < nIndex; jj++)
               {
                  if (minIndices[jj] != -1)
                  {
                     if (sampleOutputs[ss*nOutputs+outputID] != 
                         sampleOutputs[minIndices[jj]*nOutputs+outputID]) 
                         flag = 1;
                  }
               }
               if (flag == 1) 
               {
                  for (jj = 0; jj < nInputs; jj++) 
                     tempX[count*nInputs+jj] = sampleInputs[ss*nInputs+jj];
                  count++;
               }
            }
            delete [] minIndices;
            delete [] distPairs;
            if (count < nInputs+1)
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
                     for (ii = 0; ii < nInputs; ii++) 
                        printf("%e ", tempX[ss*nInputs+ii]);
                     printf("\n");
                  }
               }
               ioPtr = new PsuadeData();
               printf("Since the interface is at most m-1 dimension where\n");
               printf("m is the number of inputs. Please select one input\n");
               printf("to be the dependent variable in the polynomial equation\n");
               printf("describing the interface in the form of :\n");
               printf("   X(k) = a polynomial in terms of all other X(i), i != k.\n");
               sprintf(pString, "Select k (1 - %d) : ", nInputs);
               kk = getInt(1, nInputs, pString);
               kk--;
               double *jLowerB = new double[nInputs-1];
               double *jUpperB = new double[nInputs-1];
               jj = 0;
               for (ii = 0; ii < nInputs; ii++)
               {
                  if (ii != kk)
                  {
                     jLowerB[jj] = iLowerB[ii];
                     jUpperB[jj] = iUpperB[ii];
                     jj++;
                  }
               }
               tempY  = new double[count];
               states = new int[count];
               for (ss = 0; ss < count; ss++) 
               {
                  tempY[ss] = tempX[ss*nInputs+kk];
                  jj = 0;
                  for (ii = 0; ii < nInputs; ii++)
                  {
                     if (ii != kk)
                     {
                        tempX[ss*(nInputs-1)+jj] = tempX[ss*nInputs+ii];
                        jj++;
                     }
                  }
               }
               ioPtr->updateInputSection(count, nInputs-1, NULL, jLowerB,
                                         jUpperB, tempX, inputNames);
               for (ss = 0; ss < count; ss++) states[ss] = 1;
               ioPtr->updateOutputSection(count, 1, tempY, states,outputNames);
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
      }

      // RS-based UA with posterior sample
      // ----------------------------------------------------------------
      else if (!strcmp(command, "rs_uap") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rs_uap: rs_ua when the response surface is appended\n");
            printf("        with a discrepancy model and the sample is.\n");
            printf("        provided by users.\n");
            printf("syntax: rs_uap (no argument needed)\n");
            continue;
         }
         printf("ERROR: no data (load data first).\n");
      }
      else if (!strcmp(command, "rs_uap") && nInputs > 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rs_uap: rs_ua when the response surface is appended\n");
            printf("        with a discrepancy model and the sample is.\n");
            printf("        provided by users.\n");
            printf("syntax: rs_uap (no argument needed)\n");
            continue;
         }
         if (psuadeIO_ == NULL) printf("ERROR: data not loaded yet.\n");
         else if (nOutputs > 1)
         {
            printf("Currently this command does not support nOutputs > 1.\n");
            printf("Use 'write' to choose one output for processing.\n");
            continue;
         }
         else
         {
            int discFile=1, nInps, nOuts;
            faFlag = 3;
            psuadeIO_->getParameter("ana_outputid", pPtr);
            kk = pPtr.intData_;
            outputID = 0;
            faPtrsRsEval = new FuncApprox*[2];
            psuadeIO_->updateAnalysisSection(-1, -1, -1, -1, outputID, -1);
            faPtrsRsEval[0] = genFAInteractive(psuadeIO_, faFlag);
            faPtrsRsEval[1] = NULL;
            psuadeIO_->updateAnalysisSection(-1, -1, -1, -1, kk, -1);
            printf("Enter discrepancy model PSUADE file (if none, enter NONE): ");
            scanf("%s", winput);
            fgets(lineIn,500,stdin); 
            if (!strcmp(winput, "NONE")) discFile = 0; 
            else
            {
               ioPtr = new PsuadeData();
               status = ioPtr->readPsuadeFile(winput);
               if (status == 0)
               {
                  ioPtr->getParameter("input_ninputs", pPtr);
                  nInps = pPtr.intData_;
                  if (nInps < nInputs)
                  {
                     printf("Discrepancy model has %d inputs.\n", nInps);
                     printf("So the first %d inputs in the model file ",nInps);
                     printf("are assumed to associate with the inputs of\n");
                     printf("the discrepancy model.\n");
                  }
                  ioPtr->getParameter("output_noutputs", pPtr);
                  nOuts = pPtr.intData_;
                  if (nOuts > 1)
                  {
                     printf("The discrepancy model has nOutputs > 1.\n");
                     printf("This is currently not supported.\n");
                     continue;
                  }
                  faPtrsRsEval[1] = genFAInteractive(ioPtr, 3);
                  delete ioPtr;
               }
               else
               {
                  printf("ERROR: in reading the discrepancy model file %s.\n",winput);
                  discFile = 0;
                  continue;
               }
               ioPtr = NULL;
            }

            int dnInps, dnSamp;
            double *inputVals=NULL, *outVals;
            printf("Enter sample file: ");
            scanf("%s", dataFile);
            fgets(lineIn,500,stdin); 
            fp = fopen(dataFile, "r");
            if (fp == NULL)
            {
               printf("ERROR: sample data file %s not found.\n", dataFile);
               continue;
            }
            else
            {
               fscanf(fp, "%s", winput);
               if (strcmp(winput, "PSUADE_BEGIN"))
               {
                  printf("ERROR: file must begin with PSUADE_BEGIN\n");
                  fclose(fp);
                  printf("File format: \n");
                  printf("PSUADE_BEGIN \n");
                  printf("<nPts> <nInputs> \n");
                  printf("1 <input 1> <input 2> ... <output>\n");
                  printf("2 <input 1> <input 2> ... <output>\n");
                  printf("...... \n");
                  printf("PSUADE_END \n");
                  continue;
               }
               else
               {
                  fscanf(fp, "%d %d", &dnSamp, &dnInps);
                  if (dnSamp <= 0)
                  {
                     printf("ERROR: invalid sample size\n");
                     fclose(fp);
                     continue;
                  }
                  printf("Sample size read = %d\n", dnSamp);
                  if (dnInps != nInputs)
                  {
                     printf("ERROR: input size does not match nInputs.\n");
                     printf(":      input size in local memory = %d.\n", 
                            nInputs);
                     printf(":      input size from file       = %d.\n",dnInps);
                     fclose(fp);
                     continue;
                  }
                  inputVals = new double[dnSamp*dnInps];
                  outVals = new double[dnSamp];
                  for (jj = 0; jj < dnSamp; jj++)
                  {
                     fscanf(fp, "%d", &ind);
                     if (ind != (jj+1))
                     {
                        printf("ERROR: input index mismatch (%d,%d)\n",
                               jj+1,ind);
                        printf("       read     index = %d\n", ind);
                        printf("       expected index = %d\n", jj+1);
                        delete [] inputVals;
                     delete [] outVals;
                        break;
                     }
                     for (ii = 0; ii < nInputs; ii++)
                        fscanf(fp, "%lg", &(inputVals[jj*dnInps+ii]));
                  }
                  fscanf(fp, "%s", winput);
                  fscanf(fp, "%s", winput);
                  if (strcmp(winput, "PSUADE_END"))
                  {
                     fclose(fp);
                     printf("ERROR: file must end with PSUADE_END\n");
                     delete [] inputVals;
                     delete [] outVals;
                     continue;
                  }
               }
               fclose(fp);
            }
            faPtrsRsEval[0]->evaluatePoint(dnSamp, inputVals, outVals);
            if (discFile == 1)
            {
               for (jj = 0; jj < dnSamp; jj++)
               {
                  dtemp = faPtrsRsEval[1]->evaluatePoint(&inputVals[jj*nInputs]);
                  outVals[jj] += dtemp;
               }
            }
            double mean=0.0, stdev=0.0;
            for (jj = 0; jj < dnSamp; jj++) mean += outVals[jj];
            mean /= (double) dnSamp;
            for (jj = 0; jj < dnSamp; jj++)
               stdev += pow(outVals[jj] - mean, 2.0);
            stdev = sqrt(stdev / dnSamp);
            if (psPlotTool_ == 1)
            {
               fp = fopen("scilabrsuap.sci", "w");
               if (fp == NULL)
               {
                  printf("rs_uap ERROR: cannot open scilabrsupa.sci file.\n");
                  continue;
               }
            }
            else
            {
               fp = fopen("matlabrsuap.m", "w");
               if (fp == NULL)
               {
                  printf("rs_uap ERROR: cannot open matlabrsupa.m file.\n");
                  continue;
               }
            }
            fprintf(fp, "Y = [ \n");
            for (jj = 0; jj < dnSamp; jj++)
               fprintf(fp, "%16.8e\n", outVals[jj]);
            fprintf(fp, "];\n");
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
               fprintf(fp, "[nk,xk]=hist(Y,10);\n");
               fprintf(fp, "bar(xk,nk/%d,1.0)\n",dnSamp);
            }
            fwritePlotAxes(fp);
            fwritePlotTitle(fp, "Probability Distribution");
            fwritePlotXLabel(fp, "Output Value");
            fwritePlotYLabel(fp, "Probabilities");
            fclose(fp);
            if (psPlotTool_ == 1)
                 printf("rs_uap: distribution in scilabrsuap.scin");
            else printf("rs_uap: distribution in matlabrsuap.m.\n");
            printf("Sample mean  = %e\n", mean);
            printf("Sample stdev = %e\n", stdev);
            delete [] inputVals;
            delete [] outVals;
            delete faPtrsRsEval[0];
            if (discFile == 1) delete faPtrsRsEval[1];
            delete [] faPtrsRsEval;
            faPtrsRsEval = NULL;
         }
      }

      // RS-based UA with bootstrap and can be with posterior
      // ----------------------------------------------------------------
      else if (!strcmp(command, "rs_uab") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rs_uab: this is a generic RS-based command that can\n");
            printf("        accommodate a discrepancy model, a pre-generated\n");
            printf("        sample (in a PSUADE file), and bootstrapping.\n");
            printf("syntax: rs_uab (no argument needed)\n");
            continue;
         }
         printf("ERROR: no data (load data first).\n");
      }
      else if (!strcmp(command, "rs_uab") && nInputs > 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rs_uab: this is a generic RS-based command that can\n");
            printf("        accommodate a discrepancy model, a pre-generated\n");
            printf("        sample (in a PSUADE file), and bootstrapping.\n");
            printf("syntax: rs_uab (no argument needed)\n");
            continue;
         }
         if (psuadeIO_ == NULL) printf("ERROR: data not loaded yet.\n");
         else if (nOutputs > 1)
         {
            printf("Currently this command does not support nOutputs > 1.\n");
            printf("Use 'write' to generate a one-output data file first.\n");
            continue;
         }
         else
         {
            int discFile=1, nInps, nOuts;
            outputID = 0;
            faPtrsRsEval = new FuncApprox*[2];
            faPtrsRsEval[1] = NULL;

            printf("Enter discrepancy model PSUADE file (if none, just 'n'): ");
            scanf("%s", winput);
            fgets(lineIn,500,stdin); 
            if (winput[0] == 'n') discFile = 0; 
            else
            {
               ioPtr = new PsuadeData();
               status = ioPtr->readPsuadeFile(winput);
               if (status == 0)
               {
                  ioPtr->getParameter("input_ninputs", pPtr);
                  nInps = pPtr.intData_;
                  if (nInps < nInputs)
                  {
                     printf("Discrepancy model has %d inputs.\n", nInps);
                     printf("So the first %d inputs in the model file ",nInps);
                     printf("are assumed to associate with the inputs of\n");
                     printf("the discrepancy model.\n");
                  }
                  ioPtr->getParameter("output_noutputs", pPtr);
                  nOuts = pPtr.intData_;
                  if (nOuts > 1)
                  {
                     printf("The discrepancy model has nOutputs > 1.\n");
                     printf("This is currently not supported.\n");
                     continue;
                  }
                  faPtrsRsEval[1] = genFAInteractive(ioPtr, 3);
                  delete ioPtr;
               }
               else
               {
                  printf("ERROR: in reading the discrepancy model file %s.\n",winput);
                  discFile = 0;
                  delete [] faPtrsRsEval;
                  faPtrsRsEval = NULL;
                  continue;
               }
               ioPtr = NULL;
            }

            int    dnSamp;
            double *inputVals=NULL;
            printf("Enter sample file: ");
            scanf("%s", dataFile);
            fgets(lineIn,500,stdin); 
            fp = fopen(dataFile, "r");
            if (fp == NULL)
            {
               printf("ERROR: sample data file %s not found.\n", dataFile);
               if (discFile == 1) delete faPtrsRsEval[1];
               delete [] faPtrsRsEval;
               faPtrsRsEval = NULL;
               if (ioPtr != NULL) delete ioPtr;
               ioPtr = NULL;
               continue;
            }
            else
            {
               fscanf(fp, "%s", winput);
               if (strcmp(winput, "PSUADE_BEGIN"))
               {
                  printf("ERROR: file must begin with PSUADE_BEGIN\n");
                  fclose(fp);
                  printf("File format: \n");
                  printf("PSUADE_BEGIN \n");
                  printf("<nPts> <nInputs> \n");
                  printf("1 <input 1> <input 2> ... <output>\n");
                  printf("2 <input 1> <input 2> ... <output>\n");
                  printf("...... \n");
                  printf("PSUADE_END \n");
                  delete [] faPtrsRsEval;
                  faPtrsRsEval = NULL;
                  if (ioPtr != NULL) delete ioPtr;
                  ioPtr = NULL;
                  continue;
               }
               else
               {
                  fscanf(fp, "%d %d", &dnSamp, &kk);
                  if (dnSamp <= 0)
                  {
                     printf("ERROR: invalid sample size\n");
                     fclose(fp);
                     delete [] faPtrsRsEval;
                     faPtrsRsEval = NULL;
                     if (ioPtr != NULL) delete ioPtr;
                     ioPtr = NULL;
                     continue;
                  }
                  if (kk != nInputs)
                  {
                     printf("ERROR: input size does not match nInputs.\n");
                     printf(":      input size in local memory = %d.\n", 
                            nInputs);
                     printf(":      input size from file       = %d.\n",kk);
                     fclose(fp);
                     delete [] faPtrsRsEval;
                     faPtrsRsEval = NULL;
                     if (ioPtr != NULL) delete ioPtr;
                     ioPtr = NULL;
                     continue;
                  }
                  inputVals = new double[dnSamp*nInputs];
                  for (jj = 0; jj < dnSamp; jj++)
                  {
                     fscanf(fp, "%d", &ind);
                     if (ind != (jj+1))
                     {
                        printf("ERROR: input index mismatch (%d,%d)\n",
                               jj+1,ind);
                        printf("       read     index = %d\n", ind);
                        printf("       expected index = %d\n", jj+1);
                        delete [] inputVals;
                        delete [] faPtrsRsEval;
                        faPtrsRsEval = NULL;
                        if (ioPtr != NULL) delete ioPtr;
                        ioPtr = NULL;
                        fclose(fp);
                        break;
                     }
                     for (ii = 0; ii < nInputs; ii++)
                        fscanf(fp, "%lg", &(inputVals[jj*nInputs+ii]));
                  }
                  if (jj != dnSamp) continue;
                  fscanf(fp, "%s", winput);
                  fscanf(fp, "%s", winput);
                  if (strcmp(winput, "PSUADE_END"))
                  {
                     fclose(fp);
                     printf("ERROR: file must end with PSUADE_END\n");
                     delete [] inputVals;
                     delete [] faPtrsRsEval;
                     faPtrsRsEval = NULL;
                     if (ioPtr != NULL) delete ioPtr;
                     ioPtr = NULL;
                     continue;
                  }
               }
               fclose(fp);
            }
            int nbs;
            sprintf(pString, "How many bootstrapped samples to use (10 - 1000) : ");
            nbs = getInt(1, 1000, pString);
            printf("Write the CDFs to a matlab/scilab file? (y or n) ");
            scanf("%s", winput);
            fgets(lineIn,500,stdin); 
            flag = 0;
            if (winput[0] == 'y')
            {
               if (dnSamp > 50000)
               {
                  printf("INFO: sample size %d too large (>50000) for matlab plot.\n",
                         dnSamp);
                  printf("      CDF plots not to be generated.\n");
               }
               else
               {
                  flag = 1; 
                  if (psPlotTool_ == 1)
                  {
                     fp = fopen("scilabrsuab_cdf.sci", "w");
                     if (fp == NULL)
                     {
                        printf("ERROR: cannot open file.\n");
                        flag = 0;
                     }
                     else
                     {
                        fprintf(fp, "// CDFs for rs_uab\n");
                        fwritePlotCLF(fp);
                     }
                  }
                  else
                  {
                     fp = fopen("matlabrsuab_cdf.m", "w");
                     if (fp == NULL)
                     {
                        printf("ERROR: cannot open file.\n");
                        flag = 0;
                     }
                     else
                     {
                        fprintf(fp, "%% CDFs for rs_uab\n");
                        fwritePlotCLF(fp);
                     }
                  }
               }
            }
            double mean=0.0, stdev=0.0;
            double *outVals;
            if (nbs == 1) nSamples2 = nSamples;
            else
            {
               nSamples2 = (int) (0.9 * nSamples);
               if ((double) nSamples2 / (double) nSamples < 0.9) nSamples2++;
            }
            faPtrsRsEval[0] = genFA(-1, nInputs, -1, nSamples2);
            if (faPtrsRsEval[0] == NULL)
            {
               printf("ERROR: cannot generate response surface.\n");
               delete [] faPtrsRsEval;
               faPtrsRsEval = NULL;
               if (inputVals != NULL) delete [] inputVals;
               inputVals = NULL;
               if (ioPtr != NULL) delete ioPtr;
               ioPtr = NULL;
               continue;
            }
            faPtrsRsEval[0]->setBounds(iLowerB, iUpperB);
            faPtrsRsEval[0]->setOutputLevel(0);
            outVals = new double[dnSamp];
            tempX = new double[nSamples*nInputs];
            tempY = new double[nSamples];
            tempI = new int[nSamples];
            tempV = new double[nbs];
            tempW = new double[nbs];
            for (it = 0; it < nbs; it++)
            {
               printf("rs_uab: iteration %d\n", it+1);
               if (nbs == 1)
               {
                  for (jj = 0; jj < nSamples*nInputs; jj++)
                     tempX[jj] = sampleInputs[jj];
                  for (jj = 0; jj < nSamples; jj++)
                     tempY[jj] = sampleOutputs[jj*nOutputs+outputID];
               }
               else
               {   
                  for (jj = 0; jj < nSamples; jj++) tempI[jj] = 0;
                  kk = 0;
                  while (kk < nSamples2)
                  {
                     ind = PSUADE_rand() % nSamples;
                     if (tempI[ind] == 0)
                     {
                        for (ii = 0; ii < nInputs; ii++)
                           tempX[kk*nInputs+ii] = sampleInputs[ind*nInputs+ii];
                        tempY[kk] = sampleOutputs[ind*nOutputs+outputID];
                        tempI[ind] = 1;
                        kk++;
                     }
                  }
               }
               if (discFile == 1)
               {
                  for (jj = 0; jj < dnSamp; jj++)
                  {
                     dtemp = faPtrsRsEval[1]->evaluatePoint(&tempX[jj*nInputs]);
                     tempY[jj] += dtemp;
                  }
               }
               kk = -999;
               status = faPtrsRsEval[0]->genNDGridData(tempX,tempY,&kk,NULL,NULL);
               faPtrsRsEval[0]->evaluatePoint(dnSamp, inputVals, outVals);
               mean = stdev = 0.0;
               for (jj = 0; jj < dnSamp; jj++) mean += outVals[jj];
               mean /= (double) dnSamp;
               for (jj = 0; jj < dnSamp; jj++)
                  stdev += pow(outVals[jj] - mean, 2.0);
               stdev = sqrt(stdev / dnSamp);
               tempV[it] = mean;
               tempW[it] = stdev;
               if (fp != NULL && flag == 1)
               {
                  fprintf(fp, "Y = [\n");
                  for (jj = 0; jj < dnSamp; jj++) fprintf(fp,"%e\n",outVals[jj]);
                  fprintf(fp, "];\n");
                  fprintf(fp, "Y%d = sort(Y);\n",it+1);
                  fprintf(fp, "X%d = (1 : %d)';\n", it+1, dnSamp);
                  fprintf(fp, "X%d = X%d / %d;\n", it+1, it+1, dnSamp);
                  if (it == 0)
                  {
                     fprintf(fp, "YY = Y%d;\n", it+1);
                     fprintf(fp, "XX = X%d;\n", it+1);
                  }
                  else
                  {
                     fprintf(fp, "YY = [YY Y%d];\n", it+1);
                     fprintf(fp, "XX = [XX X%d];\n", it+1);
                  }
               }
            }
            if (fp != NULL)
            {
               fprintf(fp, "plot(YY, XX, 'lineWidth',3)\n");
               fwritePlotTitle(fp, "Cumulative Distribution");
               fwritePlotAxes(fp);
               fwritePlotXLabel(fp, "Output Value");
               fwritePlotYLabel(fp, "Probabilities");
               fclose(fp);
               if (psPlotTool_ == 1)
                    printf("rs_uab: scilabrsuab_cdf.sci has the CDF plots.\n");
               else printf("rs_uab: matlabrsuab_cdf.m has the CDF plots.\n");
            }
            delete faPtrsRsEval[0];
            faPtrsRsEval[0] = NULL;
            if (discFile == 1) delete faPtrsRsEval[1];
            delete [] faPtrsRsEval;
            faPtrsRsEval = NULL;
            delete [] inputVals;
            delete [] outVals;
            delete [] tempX;
            delete [] tempY;
            delete [] tempI;
            if (it == nbs && nbs > 1)
            {
               mean = 0.0;
               for (jj = 0; jj < nbs; jj++) mean += tempV[jj];
               mean /= (double) nbs;
               for (jj = 0; jj < nbs; jj++)
                  stdev += pow(tempV[jj] - mean, 2.0);
               stdev = sqrt(stdev / nbs);
               printf("rs_uab: Sample Mean  = %e (%e)\n", mean, stdev);
               mean = 0.0;
               for (jj = 0; jj < nbs; jj++) mean += tempW[jj];
               mean /= (double) nbs;
               for (jj = 0; jj < nbs; jj++)
                  stdev += pow(tempW[jj] - mean, 2.0);
               stdev = sqrt(stdev / nbs);
               printf("rs_uab: Sample Stdev = %e (%e)\n", mean, stdev);
            }
            else if (kk == nbs && nbs == 1)
            {
               printf("rs_uab: Sample Mean  = %e\n", tempV[0]);
               printf("rs_uab: Sample Stdev = %e\n", tempW[0]);
            }
            delete [] tempW;
            delete [] tempV;
         }
      }

      // convert data set based on distribution
      // ----------------------------------------------------------------
      else if (!strcmp(command, "pdfconvert") && nInputs <= 0)
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
         printf("ERROR: no data (load data first).\n");
      }
      else if (!strcmp(command, "pdfconvert") && nInputs > 0)
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
         if (psuadeIO_ == NULL) printf("ERROR: data not loaded yet.\n");
         else
         {
            pdfman = new PDFManager();
            psuadeIO_->getParameter("ana_diagnostics",pPtr);
            kk = pPtr.intData_;
            psuadeIO_->updateAnalysisSection(-1,-1,-1,outputLevel_,-1, -1);
            psuadeIO_->updateMethodSection(PSUADE_SAMP_MC, -1, -1, -1, -1);
            pdfman->initialize(psuadeIO_);
            vecIn.load(nSamples*nInputs, sampleInputs);
            vecOut.setLength(nSamples*nInputs);
            vecUpper.load(nInputs, iUpperB);
            vecLower.load(nInputs, iLowerB);
            pdfman->invCDF(nSamples, vecIn, vecOut, vecLower, vecUpper);
            psuadeIO_->updateAnalysisSection(-1,-1,-1,ii,-1, -1);
            for (ii = 0; ii < nSamples*nInputs; ii++)
               sampleInputs[ii] = vecOut[ii];
            psuadeIO_->updateInputSection(nSamples,nInputs,NULL,
                                          NULL,NULL,sampleInputs,NULL); 
            psuadeIO_->updateAnalysisSection(-1, -1, -1, kk, -1, 0);
            delete pdfman;
            printf("The sample in local memory has been converted based on\n");
            printf("the PDF information in the INPUT section.\n");
            printf("You can now store your new sample using 'write'.\n");
         }
      }

      // draw a random sample from the resident sample
      // ----------------------------------------------------------------
      else if (!strcmp(command, "rand_draw") && nInputs <= 0)
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
         printf("ERROR: no data (load data first).\n");
      }
      else if (!strcmp(command, "rand_draw") && nInputs > 0)
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
         if (psuadeIO_ == NULL) printf("ERROR: data not loaded yet.\n");
         else
         {
            sprintf(pString,"Size of the sample to be drawn : (1-1000000) ");
            count = getInt(1, 1000000, pString);
            tempX = new double[nInputs * count];
            tempY = new double[nOutputs * count];
            tempI = new int[count];;
            for (ii = 0; ii < count; ii++)
            {
               ind = PSUADE_rand() % nSamples;
               for (jj = 0; jj < nInputs; jj++)
                  tempX[ii*nInputs+jj] = sampleInputs[ind*nInputs+jj];
               for (jj = 0; jj < nOutputs; jj++)
                  tempY[ii*nOutputs+jj] = sampleOutputs[ind*nOutputs+jj];
               tempI[ii] = sampleStates[ind];
            }
            ioPtr = new PsuadeData();
            ioPtr->updateInputSection(count, nInputs, NULL, iLowerB,
                                      iUpperB, tempX, inputNames);
            ioPtr->updateOutputSection(count, nOutputs, tempY, tempI,
                                       outputNames);
            ioPtr->updateMethodSection(PSUADE_SAMP_MC, count, 1, -1, -1);
            printf("Store random sample to : (filename) ");
            scanf("%s", dataFile);
            fgets(lineIn,500,stdin); 
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
      }

      // draw a random sample from the 2 samples
      // ----------------------------------------------------------------
      else if (!strcmp(command, "rand_draw2"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rand_draw2: draw a sample randomly from 2 samples\n");
            printf("syntax: rand_draw2 <n>\n");
            printf("That is, this command generates a bootstrapped sample\n");
            printf("from the sample which has been loaded to local memory\n");
            printf("(mix and match). This command is used if you have 2\n");
            printf("samples with posterior distributions for 2 sets of\n");
            printf("inputs, and you would like to propagate these 2 sets\n");
            printf("of distributions through the coupled model which\n");
            printf("comprises the two sets of inputs.\n");
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
         ioPtr = new PsuadeData;
         status = ioPtr->readPsuadeFile(winput);
         if (status != 0)
         {
            printf("rand_draw2 FILE READ ERROR: file = %s\n", winput);
            continue;
         }

         psuadeIO_->getParameter("input_ninputs", pPtr);
         nInputs = pPtr.intData_;
         psuadeIO_->getParameter("method_nsamples", pPtr);
         nSamples = pPtr.intData_;
         psuadeIO_->getParameter("input_sample", pPtr);
         if (sampleInputs != NULL) delete [] sampleInputs;
         sampleInputs = pPtr.dbleArray_;
         pPtr.dbleArray_ = NULL;

         ioPtr->getParameter("input_ninputs", pPtr);
         nInputs2 = pPtr.intData_;
         ioPtr->getParameter("method_nsamples", pPtr);
         nSamples2 = pPtr.intData_;
         ioPtr->getParameter("input_sample", pPtr);
         tempX = pPtr.dbleArray_;
         pPtr.dbleArray_ = NULL;

         iLowerB = new double[nInputs+nInputs2];
         iUpperB = new double[nInputs+nInputs2];
         pLower.clean();
         psuadeIO_->getParameter("input_lbounds", pLower);
         tempT = pLower.dbleArray_;
         pLower.dbleArray_ = NULL;
         for (ii = 0; ii < nInputs; ii++) iLowerB[ii] = tempT[ii];
         delete [] tempT;
         ioPtr->getParameter("input_lbounds", pLower);
         tempT = pLower.dbleArray_;
         pLower.dbleArray_ = NULL;
         for (ii = 0; ii < nInputs2; ii++) iLowerB[nInputs+ii] = tempT[ii];
         delete [] tempT;

         pUpper.clean();
         psuadeIO_->getParameter("input_ubounds", pUpper);
         tempT = pUpper.dbleArray_;
         pUpper.dbleArray_ = NULL;
         for (ii = 0; ii < nInputs; ii++) iUpperB[ii] = tempT[ii];
         delete [] tempT;
         ioPtr->getParameter("input_ubounds", pUpper);
         tempT = pUpper.dbleArray_;
         pUpper.dbleArray_ = NULL;
         for (ii = 0; ii < nInputs2; ii++) iUpperB[nInputs+ii] = tempT[ii];
         delete [] tempT;

         sprintf(pString,"Size of the sample to be drawn : (1-2000000) ");
         count = getInt(1, 200000, pString);
         tempW = new double[count*(nInputs+nInputs2)];
         tempY = new double[count];
         tempI = new int[count]; 
         for (ii = 0; ii < count; ii++)
         {
            ind  = PSUADE_rand() % nSamples;
            ind2 = PSUADE_rand() % nSamples2;
            for (jj = 0; jj < nInputs; jj++)
               tempW[ii*(nInputs+nInputs2)+jj] = 
                  sampleInputs[ind*nInputs+jj];
            for (jj = 0; jj < nInputs2; jj++)
               tempW[ii*(nInputs+nInputs2)+nInputs+jj] = 
                  tempX[ind2*nInputs2+jj];
            tempY[ii] = PSUADE_UNDEFINED;
            tempI[ii] = 0;
         } 
         delete psuadeIO_;
         delete ioPtr;
         psuadeIO_ = NULL;
         ioPtr = new PsuadeData();
         inputNames = new char*[nInputs+nInputs2];
         for (ii = 0; ii < nInputs+nInputs2; ii++)
         {
            inputNames[ii] = new char[200]; 
            sprintf(inputNames[ii], "X%d", ii+1);
         }
         ioPtr->updateInputSection(count, nInputs+nInputs2, NULL, 
                              iLowerB, iUpperB, tempW, inputNames);
         nOutputs = 1;
         outputNames = new char*[1];
         outputNames[0] = new char[200]; 
         sprintf(outputNames[0], "Y");
         ioPtr->updateOutputSection(count, nOutputs, tempY, tempI,
                                    outputNames);
         ioPtr->updateMethodSection(PSUADE_SAMP_MC, count, 1, -1, -1);
         printf("Store random sample to : (filename) ");
         scanf("%s", dataFile);
         fgets(lineIn,500,stdin); 
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
         delete [] sampleInputs;
         delete [] iLowerB;
         delete [] iUpperB;
         sampleInputs = NULL;
         iLowerB = NULL;
         iUpperB = NULL;
         for (ii = 0; ii < nInputs+nInputs2; ii++)
            delete [] inputNames[ii];
         delete [] inputNames;
         inputNames = NULL;
         delete [] outputNames[0];
         delete [] outputNames;
         outputNames = NULL;
         delete ioPtr;
      }

      // Changes the default name of the outputFile 
      // ----------------------------------------------------------------
      else if (!strcmp(command, "output_file"))
      {
        sscanf(lineIn,"%s %s",command,winput);
        if (!strcmp(winput, "-h"))
        {
          printf("output_file: changes the default output filename.\n");
          printf("syntax: output_file (no argument needed)\n");
          continue;
        }
        psOutputFilename_ = PSUADE_strdup(winput);
        if (psInputFilename_ == psOutputFilename_)
        {
          printf("WARNING: The output filename is the same as the input filename\n");
          printf("         If you save the file it will overwrite the input file.\n");
        }
        
      }

      // a short guide to set up
      // ----------------------------------------------------------------
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

      // generate an application workflow
      // ----------------------------------------------------------------
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

      // generate an input file for psuade 
      // ----------------------------------------------------------------
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
         if (sampleInputs  != NULL) delete [] sampleInputs;
         if (sampleOutputs != NULL) delete [] sampleOutputs;
         if (iLowerB       != NULL) delete [] iLowerB;
         if (iUpperB       != NULL) delete [] iUpperB;
         if (sampleStates  != NULL) delete [] sampleStates;
         if (inputNames != NULL)
         {
            for (ii = 0; ii < nInputs; ii++)
               delete [] inputNames[ii];
            delete [] inputNames;
         }
         if (outputNames != NULL)
         {
            for (ii = 0; ii < nOutputs; ii++)
               delete [] outputNames[ii];
            delete [] outputNames;
         }
         sampleInputs = sampleOutputs = NULL;
         sampleStates = NULL;
         iLowerB = iUpperB = NULL;
         inputNames = outputNames = NULL;

         if (status == 0)
         {
            printf("PSUADE can create the sample input data file for you.\n");
            sprintf(pString,"Create the file ? (y or n) ");
            getString(pString, winput);
            if (winput[0] == 'y')
            {
               getInputFromFile(dataFile);
               cleanUp();
               if (psuadeIO_ != NULL) delete psuadeIO_;
               if (sampler_ != NULL) SamplingDestroy(sampler_);
               sampler_ = NULL;
               if (inputNames != NULL) 
               {
                  for (ii = 0; ii < nInputs; ii++) delete [] inputNames[ii];
                  delete [] inputNames;
               }
               if (outputNames != NULL) 
               {
                  for (ii = 0; ii < nOutputs; ii++) delete [] outputNames[ii];
                  delete [] outputNames;
               }
               psuadeIO_ = new PsuadeData();
               psuadeIO_->setOutputLevel(0);
               strcpy(winput, psInputFilename_);
               status = psuadeIO_->readPsuadeFile(winput);
               psuadeIO_->getParameter("input_ninputs", pPtr);
               nInputs = pPtr.intData_;
               pINames.clean();
               psuadeIO_->getParameter("input_names", pINames);
               names = pINames.strArray_;
               pLower.clean();
               psuadeIO_->getParameter("input_lbounds", pLower);
               iLowerB = pLower.dbleArray_;
               pLower.dbleArray_ = NULL;
               pUpper.clean();
               psuadeIO_->getParameter("input_ubounds", pUpper);
               iUpperB = pUpper.dbleArray_;
               pUpper.dbleArray_ = NULL;
               inputNames = new char*[nInputs+1];
               for (ii = 0; ii < nInputs; ii++)
               {
                  inputNames[ii] = new char[200]; 
                  strcpy(inputNames[ii], names[ii]);
               }
               psuadeIO_->getParameter("output_noutputs", pPtr);
               nOutputs = pPtr.intData_;
               pONames.clean();
               psuadeIO_->getParameter("output_names", pONames);
               names = pONames.strArray_;
               outputNames = new char*[nOutputs+1];
               for (ii = 0; ii < nOutputs; ii++)
               {
                  outputNames[ii] = new char[200]; 
                  strcpy(outputNames[ii], names[ii]);
               }
               psuadeIO_->getParameter("method_sampling", pPtr);
               samplingMethod = pPtr.intData_;
               psuadeIO_->getParameter("method_nsamples", pPtr);
               nSamples = pPtr.intData_;
               psuadeIO_->getParameter("method_nreplications",pPtr);
               nReps = pPtr.intData_;
               psuadeIO_->getParameter("input_sample", pPtr);
               sampleInputs  = pPtr.dbleArray_;
               psuadeIO_->getParameter("output_sample", pPtr);
               sampleOutputs  = pPtr.dbleArray_;
               psuadeIO_->getParameter("output_states", pPtr);
               sampleStates  = pPtr.intArray_;
               pPtr.intArray_ = NULL;
               pPtr.dbleArray_ = NULL;
               pINames.clean();
               pONames.clean();
               printf("==================================================\n");
               printf("The sample matrix is now stored in %s file.\n", psOutputFilename_);
               printf("You can also use genmars command to convert the\n");
               printf("sample matrix to a row-column format.\n");
               printf("==================================================\n");
            }
         }
      }

      // generate a LLNL-specific batch file
      // ----------------------------------------------------------------
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

      // generate an input file for psuade 
      // ----------------------------------------------------------------
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

      // generate data based on some distribution
      // ----------------------------------------------------------------
      else if (!strcmp(command, "gendist"))
      {
         int    ns, gtype;
         double *sData, slbound, subound, smean, sstdev;

         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("gendist: create a sample using some probability distributions\n");
            printf("syntax: gendist (no argument needed)\n");
            continue;
         }
         sprintf(pString, "nSamples = ");
         ns = getInt(10, 10000000, pString);
         sprintf(pString, "distribution type = (1) N (2) L (3) T (4) Beta : ");
         gtype = getInt(1, 3, pString);
         if      (gtype == 1) gtype = PSUADE_PDF_NORMAL;
         else if (gtype == 2) gtype = PSUADE_PDF_LOGNORMAL;
         else if (gtype == 3) gtype = PSUADE_PDF_TRIANGLE;
         else if (gtype == 4) gtype = PSUADE_PDF_BETA;
         sprintf(pString, "sample mean = ");
         smean = getDouble(pString);
         sprintf(pString, "sample standard deviation = ");
         sstdev = getDouble(pString);
         corMat.setDim(1,1);
         corMat.setEntry(0,0, 1.0e0);
         pdfman = new PDFManager();
         pdfman->initialize(1, &gtype, &smean, &sstdev, corMat);
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

      // generate an example
      // ----------------------------------------------------------------
      else if (!strcmp(command, "genexample"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("genexample: create a PSUADE example (driver + psuade input file)\n");
            printf("syntax: genexample (no argument needed)\n");
            continue;
         }
         genDriver(1);
         genSetup(1, winput);
         printf("Now use: psuade psuade.in to run the example.\n");
      }

      // generate an example Psuade config file
      // ----------------------------------------------------------------
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

      // check job status and generate a report
      // ----------------------------------------------------------------
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

      // list one input and one output pair
      // ----------------------------------------------------------------
      else if (!strcmp(command, "list1"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("list1: list all points with 1 selected input and 1 output\n");
            printf("syntax: list1 (no argument needed)\n");
            continue;
         }
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no data (load data first).\n");
         }
         else
         {
            sprintf(pString, "Enter input number (1 - %d) : ", nInputs);
            iInd = getInt(1, nInputs, pString);
            iInd--;
            if (nOutputs == 1) outputID = 0;
            else
            {
               sprintf(pString,"Enter output number (1 - %d) : ", nOutputs);
               outputID = getInt(1, nOutputs, pString);
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
            tempX = new double[nSamples];
            tempY = new double[nSamples];
            tempW = new double[nSamples];
            for (sInd = 0; sInd < nSamples; sInd++)
            {
               tempX[sInd] = sampleInputs[sInd*nInputs+iInd];
               tempY[sInd] = sampleOutputs[sInd*nOutputs+outputID];
               tempW[sInd] = (double) sInd + 1;
            }
            if (winput[0] == 'y') 
               sortDbleList3(nSamples, tempX, tempY, tempW);
            else if (winput[0] == 'Y') 
               sortDbleList3(nSamples, tempY, tempX, tempW);
            for (sInd = 0; sInd < nSamples; sInd++)
               printf("%6d: Sample %7d : input = %16.8e, output = %16.8e\n",
                      sInd+1, (int) tempW[sInd], tempX[sInd], tempY[sInd]);
            delete [] tempX;
            delete [] tempY;
            delete [] tempW;
            tempX = tempY = tempW = NULL;
         }
      }

      // list two input and one output pair
      // ----------------------------------------------------------------
      else if (!strcmp(command, "list2"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("list2: list all points with 2 selected inputs and 1 output\n");
            printf("syntax: list2 (no argument needed)\n");
            continue;
         }
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no data (load data first).\n");
         }
         else
         {
            sprintf(pString, "Enter input number 1 (1 - %d) : ", nInputs);
            iInd1 = getInt(1, nInputs, pString);
            iInd1--;
            sprintf(pString, "Enter input number 2 (1 - %d) : ", nInputs);
            iInd2 = getInt(1, nInputs, pString);
            iInd2--;
            if (nOutputs == 1) outputID = 0;
            else
            {
               sprintf(pString,"Enter output number (1 - %d) : ", nOutputs);
               outputID = getInt(1, nOutputs, pString);
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
            tempX = new double[nSamples];
            tempY = new double[nSamples];
            tempW = new double[nSamples];
            tempV = new double[nSamples];
            for (sInd = 0; sInd < nSamples; sInd++)
            {
               tempX[sInd] = sampleInputs[sInd*nInputs+iInd1];
               tempV[sInd] = sampleInputs[sInd*nInputs+iInd2];
               tempY[sInd] = sampleOutputs[sInd*nOutputs+outputID];
               tempW[sInd] = (double) sInd;
            }
            if (winput[0] == 'y') 
            {
               sortDbleList4(nSamples, tempX, tempV, tempY, tempW);
               sInd = 0;
               count = 1; 
               while (sInd < nSamples) 
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
               sortDbleList4(nSamples, tempY, tempX, tempV, tempW);
            }
            for (sInd = 0; sInd < nSamples; sInd++)
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

      // list all inputs and one output 
      // ----------------------------------------------------------------
      else if (!strcmp(command, "listall"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("listall: list all points with all inputs and 1 output\n");
            printf("syntax: listall (no argument needed)\n");
            continue;
         }
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no data (load data first).\n");
         }
         else
         {
            if (nOutputs == 1) outputID = 0;
            else
            {
               sprintf(pString,"Enter output number (1 - %d) : ", nOutputs);
               outputID = getInt(1, nOutputs, pString);
               outputID--;
            }
            sprintf(pString,"Sort the output (y or n) ? ");
            getString(pString, winput);
            if (winput[0] == 'y')
            {
               tempI = new int[nSamples];
               tempY = new double[nSamples];
               for (sInd = 0; sInd < nSamples; sInd++)
               {
                  tempI[sInd] = sInd;
                  tempY[sInd] = sampleOutputs[sInd*nOutputs+outputID];
               }
               sortDbleList2a(nSamples, tempY, tempI);
               for (sInd = 0; sInd < nSamples; sInd++)
               {
                  kk = tempI[sInd];
                  printf("%6d: Sample %7d : ", sInd+1, kk+1);
                  for (ii = 0; ii < nInputs; ii++)
                     printf("%12.4e ", sampleInputs[kk*nInputs+ii]);
                  printf("= %12.4e\n", tempY[sInd]);
               }
               delete [] tempY;
               delete [] tempI;
               tempY = NULL;
               tempI = NULL;
            }
            else
            {
               for (sInd = 0; sInd < nSamples; sInd++)
               {
                  printf("Sample %7d : ", sInd+1);
                  for (ii = 0; ii < nInputs; ii++)
                     printf("%12.4e ", sampleInputs[sInd*nInputs+ii]);
                  printf("= %12.4e\n", sampleOutputs[sInd]);
               }
            }
         }
      }

      // list one sample point 
      // ----------------------------------------------------------------
      else if (!strcmp(command, "disp_sample") || !strcmp(command, "sshow"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("disp_sample (or sshow): display one sample point\n");
            printf("syntax: disp_sample or sshow (no argument needed)\n");
            continue;
         }
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no data (load data first).\n");
         }
         else
         {
            sprintf(pString, "Enter sample number (1 - %d) : ", nSamples);
            iInd = getInt(1, nSamples, pString);
            iInd--;
            printf("Sample %7d : \n", iInd+1);
            for (sInd = 0; sInd < nInputs; sInd++)
               printf("   input  %3d = %16.8e\n", sInd+1, 
                      sampleInputs[iInd*nInputs+sInd]);
            for (sInd = 0; sInd < nOutputs; sInd++)
               printf("   output %3d = %16.8e\n", sInd+1, 
                      sampleOutputs[iInd*nOutputs+sInd]);
         }
      }

      // find the maximum sample point
      // ----------------------------------------------------------------
      else if (!strcmp(command, "max"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("max: find the sample point with maximum output value\n");
            printf("syntax: max (no argument needed)\n");
            continue;
         }
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no data (load data first).\n");
         }
         else
         {
            if (nOutputs == 1) outputID = 0;
            else
            {
               sprintf(pString,"Enter output number (1 - %d) : ", nOutputs);
               outputID = getInt(1, nOutputs, pString);
               outputID--;
            }
            Ymax = sampleOutputs[outputID];
            ind  = 0;
            for (sInd = 1; sInd < nSamples; sInd++)
            {
               if (sampleOutputs[sInd*nOutputs+outputID] > Ymax)
               {
                  ind  = sInd;
                  Ymax = sampleOutputs[sInd*nOutputs+outputID];
               }
            }
            printf("Sample %d gives maximum for output %d\n",
                    ind+1,outputID+1);
            for (iInd = 0; iInd < nInputs; iInd++)
               printf("  input %4d = %16.8e\n", iInd+1,
                      sampleInputs[ind*nInputs+iInd]);
            printf("  ====> output %4d = %16.8e\n", outputID+1,
                   sampleOutputs[ind*nOutputs+outputID]);
         }
      }

      // find the minimum sample point
      // ----------------------------------------------------------------
      else if (!strcmp(command, "min"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("min: find the sample point with minimum output value\n");
            printf("syntax: min (no argument needed)\n");
            continue;
         }
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no data (load data first).\n");
         }
         else
         {
            if (nOutputs == 1) outputID = 0;
            else
            {
               sprintf(pString,"Enter output number (1 - %d) : ", nOutputs);
               outputID = getInt(1, nOutputs, pString);
               outputID--;
            }
            Ymin = sampleOutputs[outputID];
            ind  = 0;
            for (sInd = 1; sInd < nSamples; sInd++)
            {
               if (sampleOutputs[sInd*nOutputs+outputID] < Ymin)
               {
                  ind  = sInd;
                  Ymin = sampleOutputs[sInd*nOutputs+outputID];
               }
            }
            printf("Sample %d gives minimum for output %d\n",
                    ind+1,outputID+1);
            for (iInd = 0; iInd < nInputs; iInd++)
               printf("  input  %4d = %16.8e\n", iInd+1,
                      sampleInputs[ind*nInputs+iInd]);
            printf("  ====> output %4d = %16.8e\n", outputID+1,
                   sampleOutputs[ind*nOutputs+outputID]);
         }
      }

      // remove a certain input
      // ----------------------------------------------------------------
      else if (!strcmp(command, "idelete"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("idelete: delete a sample input \n");
            printf("syntax: idelete (no argument needed)\n");
            continue;
         }
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no data (load data first).\n");
         }
         else
         {
            if (nInputs == 1)
            {
               printf("You have only one input left -> no deletion.\n");
               continue;
            }
            else
            {
               pINames.clean();
               psuadeIO_->getParameter("input_names", pINames);
               names = pINames.strArray_;
               for (ii = 0; ii < nInputs; ii++)
                  printf("Input %3d = %s\n", ii+1, names[ii]);
                  
               indSet = new int[nInputs];
               for (ii = 0; ii < nInputs; ii++) indSet[ii] = 1;

               sprintf(pString,"How many inputs to remove? (1-%d) ",nInputs-1);
               count = getInt(1, nInputs-1, pString);
               sprintf(pString,"Enter input number (1 - %d) : ", nInputs);
               for (ii = 0; ii < count; ii++)
               {
                  iInd = getInt(1, nInputs, pString);
                  printf("You are removing input %d (%s)\n",iInd,names[iInd-1]);
                  indSet[iInd-1] = 0;
               }
               kk = 0;
               for (ii = 0; ii < nInputs; ii++)
               {
                  if (indSet[ii] == 1)
                  {
                     iLowerB[kk] = iLowerB[ii];
                     iUpperB[kk] = iUpperB[ii];
                     kk++;
                  }
               }
               count = kk;
               for (sInd = 0; sInd < nSamples; sInd++)
               {
                  kk = 0;
                  for (ii = 0; ii < nInputs; ii++)
                  {
                     if (indSet[ii] == 1)
                     {
                        sampleInputs[sInd*count+kk] = 
                              sampleInputs[sInd*nInputs+ii];
                        kk++;
                     }
                  }
               }
               kk = 0;
               for (ii = 0; ii < nInputs; ii++)
                  if (indSet[ii] == 1)
                     strcpy(names[kk++], names[ii]);
               nInputs = count;
               psuadeIO_->updateInputSection(nSamples,nInputs,NULL,iLowerB,
                                             iUpperB,sampleInputs,names); 
               delete [] indSet;
               printf("idelete completed. Use 'write' to store.\n");
            }
         }
      }

      // remove a selected output
      // ----------------------------------------------------------------
      else if (!strcmp(command, "odelete"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("odelete: delete a sample output \n");
            printf("syntax: odelete (no argument needed)\n");
            continue;
         }
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no data (load data first).\n");
         }
         else
         {
            if (nOutputs == 1)
            {
               printf("You have only one output left -> no deletion.\n");
               continue;
            }
            else if (nOutputs > 1)
            {
               if (outputNames != NULL) 
               {
                  for (ii = 0; ii < nOutputs; ii++)
                     printf("Output %3d = %s\n", ii+1, outputNames[ii]);
               }
               sprintf(pString,"Enter output number (1 - %d) : ", nOutputs);
               outputID = getInt(1, nOutputs, pString);
               outputID--;
               if (outputNames != NULL) 
               {
                  for (ii = outputID+1; ii < nOutputs; ii++)
                     strcpy(outputNames[ii-1], outputNames[ii]);
                  delete [] outputNames[nOutputs-1];
                  outputNames[nOutputs-1] = NULL;
               }
               for (ii = 0; ii < nSamples; ii++)
               {
                  for (jj = 0; jj < outputID; jj++)
                     sampleOutputs[ii*(nOutputs-1)+jj] = 
                                   sampleOutputs[ii*nOutputs+jj];
                  for (jj = outputID+1; jj < nOutputs; jj++)
                     sampleOutputs[ii*(nOutputs-1)+jj-1] = 
                                   sampleOutputs[ii*nOutputs+jj];
               }
               nOutputs--;
               psuadeIO_->updateOutputSection(nSamples,nOutputs,
                             sampleOutputs,sampleStates,outputNames); 
               printf("odelete completed. Use 'write' to store.\n");
            }
         }
      }

      // remove a selected sample point 
      // ----------------------------------------------------------------
      else if (!strcmp(command, "sdelete"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("sdelete: delete a sample point \n");
            printf("syntax: sdelete (no argument needed)\n");
            continue;
         }
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no data (load data first).\n");
         }
         else
         {
            if (nSamples == 1)
            {
               printf("You have only one sample point left -> no deletion.\n");
               continue;
            }
            else if (nSamples > 1)
            {
               sprintf(pString,"Enter sample point to delete (1 - %d) : ", nSamples);
               ss = getInt(1, nSamples, pString);
               ss--;
               for (ii = ss; ii < nSamples; ii++)
               {
                  for (jj = 0; jj < nInputs; jj++)
                     sampleInputs[ii*nInputs+jj] = 
                                   sampleInputs[(ii+1)*nInputs+jj];
                  for (jj = 0; jj < nOutputs; jj++)
                     sampleOutputs[ii*nOutputs+jj] = 
                                   sampleOutputs[(ii+1)*nOutputs+jj];
               }
               nSamples--;
               psuadeIO_->updateInputSection(nSamples,nInputs,NULL,
                                             NULL,NULL,sampleInputs,NULL); 
               psuadeIO_->updateOutputSection(nSamples,nOutputs,
                             sampleOutputs,sampleStates,outputNames); 
               psuadeIO_->updateMethodSection(-1,nSamples,-1,-1,-1);
               printf("sdelete completed. Use 'write' to store.\n");
            }
         }
      }

      // change the order of the input parameters
      // ----------------------------------------------------------------
      else if (!strcmp(command, "ishuffle"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("ishuffle: re-arrange orders of the input parameters.\n");
            printf("syntax: ishuffle (no argument needed)\n");
            printf("NOTE: this command does not change the data file.\n");
            printf("      until a 'write' is issued.\n");
            continue;
         }
         if (psuadeIO_ == NULL || nSamples <= 0)
         {
            printf("ERROR: no data (load data first).\n");
         }
         else
         {
            indSet = new int[nInputs];
            ii = 0;
            while (ii < nInputs)
            {
               sprintf(pString, "Enter the %d-th input (1 - %d) : ",ii+1,nInputs);
               indSet[ii] = getInt(1, nInputs, pString);
               ii++;
            }
            tempX = new double[nSamples*nInputs];
            if (inputNames != NULL) names  = new char*[nInputs];
            for (ii = 0; ii < nInputs; ii++)
            {
               kk = indSet[ii]-1;
               for (jj = 0; jj < nSamples; jj++)
                  tempX[jj*nInputs+ii] = sampleInputs[jj*nInputs+kk];
               if (inputNames != NULL) names[ii] = inputNames[kk];
            }
            if (inputNames != NULL) 
            {
               for (ii = 0; ii < nInputs; ii++) 
               {
                  inputNames[ii] = names[ii];
                  names[ii] = NULL;
               }
            }
            delete [] sampleInputs;
            delete [] names;
            sampleInputs = tempX;
            tempX = NULL;
            delete [] indSet;
            printf("ishuffle completed. Use 'write' to store.\n");
         }
      }

      // tag sample points based on input values 
      // ----------------------------------------------------------------
      else if (!strcmp(command, "itag"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("itag: tag sample points based on input values\n");
            printf("syntax: itag (no argument needed)\n");
            printf("This command is useful for extracting indices of\n");
            printf("sample points with input values falling in a given\n");
            printf("range. This can be done multiple times to impose\n");
            printf("multiple filters (AND operation). Load the sample\n");
            printf("again to reset previously-set filters.\n");
            printf("NOTE: this command does not change the data file.\n");
            printf("      Only display sample index information after\n");
            printf("      tagging.\n");
            continue;
         }
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no data (load data first).\n");
         }
         else
         {
            if (tagArray == NULL && nSamples > 0) 
            {
               tagArray = new int[nSamples];
               for (ii = 0; ii < nSamples; ii++) tagArray[ii] = 1;
            } 
            inputID = 0;
            if (nInputs > 1)
            {
               sprintf(pString, "Enter input number (1 - %d) : ", nInputs);
               inputID = getInt(1, nInputs, pString);
               inputID--;
            }
            Xmin = sampleInputs[inputID];
            for (sInd = 1; sInd < nSamples; sInd++)
               if (sampleInputs[sInd*nInputs+inputID] < Xmin)
                  Xmin = sampleInputs[sInd*nInputs+inputID];
            Xmax = sampleInputs[inputID];
            for (sInd = 1; sInd < nSamples; sInd++)
               if (sampleInputs[sInd*nInputs+inputID] > Xmax)
                  Xmax = sampleInputs[sInd*nInputs+inputID];
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
            for (sInd = 0; sInd < nSamples; sInd++)
            {
               if (sampleInputs[sInd*nInputs+inputID] < threshL ||
                   sampleInputs[sInd*nInputs+inputID] > threshU)
               {
                  tagArray[sInd] = 0;
               }
            }
            printf("Sample points left after thresholding: [\n");
            for (sInd = 0; sInd < nSamples; sInd++)
               if (tagArray[sInd] == 1) printf("%d ", sInd+1);
            printf("]\n");
            printf("itag completed.\n");
         }
      }

      // tag sample points based on output values 
      // ----------------------------------------------------------------
      else if (!strcmp(command, "otag"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("otag: tag sample points based on output values\n");
            printf("syntax: otag (no argument needed)\n");
            printf("This command is useful for extracting indices of\n");
            printf("sample points with output values falling in a given\n");
            printf("range. This can be done multiple times to impose\n");
            printf("multiple filters (AND operation). Load the sample\n");
            printf("again to reset previously-set filters.\n");
            printf("Note: this command does not change the data file.\n");
            printf("      Only display sample index information after\n");
            printf("      tagging.\n");
            continue;
         }
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no data (load data first).\n");
         }
         else
         {
            if (tagArray == NULL && nSamples > 0) 
            {
               tagArray = new int[nSamples];
               for (ii = 0; ii < nSamples; ii++) tagArray[ii] = 1;
            } 
            outputID = 0;
            if (nOutputs > 1)
            {
               sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
               outputID = getInt(1, nOutputs, pString);
               outputID--;
            }
            Ymin = sampleOutputs[outputID];
            for (sInd = 1; sInd < nSamples; sInd++)
               if (sampleOutputs[sInd*nOutputs+outputID] < Ymin)
                  Ymin = sampleOutputs[sInd*nOutputs+outputID];
            Ymax = sampleOutputs[outputID];
            for (sInd = 1; sInd < nSamples; sInd++)
               if (sampleOutputs[sInd*nOutputs+outputID] > Ymax)
                  Ymax = sampleOutputs[sInd*nOutputs+outputID];
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
            for (sInd = 0; sInd < nSamples; sInd++)
            {
               if (sampleOutputs[sInd*nOutputs+outputID] < threshL ||
                   sampleOutputs[sInd*nOutputs+outputID] > threshU)
               {
                  tagArray[sInd] = 0;
               }
            }
            printf("Sample points left after thresholding: [\n");
            for (sInd = 0; sInd < nSamples; sInd++)
               if (tagArray[sInd] == 1) printf("%d ", sInd+1);
            printf("]\n");
            printf("otag completed.\n");
         }
      }

      // reset an output to certain value
      // ----------------------------------------------------------------
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
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no data (load data first).\n");
         }
         else
         {
            sprintf(pString,"Enter output number (1 - %d) : ", nOutputs);
            oInd = getInt(1, nInputs, pString);
            oInd--;
            Ymin = PSUADE_UNDEFINED;
            Ymax = -PSUADE_UNDEFINED;
            for (ii = 0; ii < nSamples; ii++)
            {
               if (sampleOutputs[ii*nOutputs+oInd] < Ymin)
                  Ymin = sampleOutputs[ii*nOutputs+oInd];
               if (sampleOutputs[ii*nOutputs+oInd] > Ymax)
                  Ymax = sampleOutputs[ii*nOutputs+oInd];
            }
            printf("Lower and upper values for this output : %e %e\n",Ymin, Ymax);
            printf("Now specify the range of current values to be reset: \n");
            sprintf(pString,"Enter the lower bound (inclusive) of this range : ");
            Ymin = getDouble(pString);
            sprintf(pString,"Enter the upper bound (inclusive) of this range : ");
            Ymax = getDouble(pString);
            sprintf(pString,"Enter the desired output value to be set to : ");
            ddata = getDouble(pString);
            for (ii = 0; ii < nSamples; ii++)
            {
               if (sampleOutputs[ii*nOutputs+oInd] >= Ymin &&
                   sampleOutputs[ii*nOutputs+oInd] <= Ymax)
                  sampleOutputs[ii*nOutputs+oInd] = ddata;
            }
            psuadeIO_->updateOutputSection(nSamples,nOutputs, sampleOutputs,
                                           sampleStates, NULL); 
            printf("oreset completed: use 'write' to store.\n");
         }
      }

      // re-range an input in the data file
      // ----------------------------------------------------------------
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
         }
         else
         {
            sprintf(pString,"Enter input number (1 - %d) : ", nInputs);
            iInd = getInt(1, nInputs, pString);
            iInd--;
            printf("Current lower bound for input %d = %24.16e\n",
                   iInd+1,iLowerB[iInd]);
            sprintf(pString,"Enter new lower bound for input %d : ",iInd+1);
            Xmin = getDouble(pString);
            printf("Current upper bound for input %d = %24.16e\n",
                   iInd+1,iUpperB[iInd]);
            sprintf(pString,"Enter new upper bound for input %d : ",iInd+1);
            Xmax = getDouble(pString);
            if (Xmin >= Xmax)
            {
               printf("ERROR: lower bound >= upper bound.\n");
               continue;
            }
            for (ii = 0; ii < nSamples; ii++)
            {
               dtemp = sampleInputs[ii*nInputs+iInd];
               dtemp = (dtemp-iLowerB[iInd]) / (iUpperB[iInd]-iLowerB[iInd]);
               sampleInputs[ii*nInputs+iInd] = dtemp * (Xmax - Xmin) + Xmin;
            }
            iLowerB[iInd] = Xmin;
            iUpperB[iInd] = Xmax;
            psuadeIO_->updateInputSection(nSamples,nInputs,NULL,iLowerB,iUpperB,
                                          sampleInputs,NULL); 
            printf("irerange completed: use 'write' to store.\n");
         }
      }

      // reset an input to certain values
      // ----------------------------------------------------------------
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
         }
         else
         {
            sprintf(pString,"Enter input number (1 - %d) : ", nInputs);
            iInd = getInt(1, nInputs, pString);
            iInd--;
            Xmin = PSUADE_UNDEFINED;
            Xmax = -PSUADE_UNDEFINED;
            for (ii = 0; ii < nSamples; ii++)
            {
               if (sampleInputs[ii*nInputs+iInd] < Xmin)
                  Xmin = sampleInputs[ii*nInputs+iInd];
               if (sampleInputs[ii*nInputs+iInd] > Xmax)
                  Xmax = sampleInputs[ii*nInputs+iInd];
            }
            printf("Lower and upper values for this input : %e %e\n",Xmin, Xmax);
            printf("Now specify the range of current values to be reset.\n");
            sprintf(pString,"Enter the lower bound (inclusive) of this range : ");
            Xmin = getDouble(pString);
            sprintf(pString,"Enter the upper bound (inclusive) of this range : ");
            Xmax = getDouble(pString);
            sprintf(pString,"Enter the desired input value to be set to : ");
            ddata = getDouble(pString);
            for (jj = 0; jj < nSamples; jj++)
            {
               if (sampleInputs[jj*nInputs+iInd] >= Xmin &&
                   sampleInputs[jj*nInputs+iInd] <= Xmax)
                  sampleInputs[jj*nInputs+iInd] = ddata;
            }
            psuadeIO_->updateInputSection(nSamples,nInputs,NULL,NULL,NULL,
                                          sampleInputs,NULL); 
            printf("ireset completed: use 'write' to store.\n");
         }
      }

      // truncate an input to integer
      // ----------------------------------------------------------------
      else if (!strcmp(command, "ifloor"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("ifloor: map an input to nearest smaller integers\n");
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
         }
         else
         {
            sprintf(pString,"Enter input number (1 - %d) : ", nInputs);
            iInd = getInt(1, nInputs, pString);
            iInd--;
            for (jj = 0; jj < nSamples; jj++)
            {
               kk = (int) sampleInputs[jj*nInputs+iInd];
               sampleInputs[jj*nInputs+iInd] = (double) kk;
            }
            psuadeIO_->updateInputSection(nSamples,nInputs,NULL,NULL,NULL,
                                          sampleInputs,NULL); 
            printf("ifloor completed: use 'write' to store.\n");
         }
      }

      // round an input to integer
      // ----------------------------------------------------------------
      else if (!strcmp(command, "iceil"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("iceil: map an input to nearest larger integers\n");
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
         }
         else
         {
            sprintf(pString,"Enter input number (1 - %d) : ", nInputs);
            iInd = getInt(1, nInputs, pString);
            iInd--;
            for (jj = 0; jj < nSamples; jj++)
            {
               kk = (int) (sampleInputs[jj*nInputs+iInd] + 0.5);
               sampleInputs[jj*nInputs+iInd] = (double) kk;
            }
            psuadeIO_->updateInputSection(nSamples,nInputs,NULL,NULL,NULL,
                                          sampleInputs,NULL); 
            printf("iceil completed: use 'write' to store.\n");
         }
      }

      // transform an input to log or power
      // ----------------------------------------------------------------
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
         }
         else
         {
            sprintf(pString,"Enter input number (1 - %d) : ", nInputs);
            iInd = getInt(1, nInputs, pString);
            iInd--;
            sprintf(pString,"Enter type: (1) log10(input), (2) 10^(input) : ");
            kk = getInt(1, 2, pString);
            if (kk == 1)
            {
               for (ss = 0; ss < nSamples; ss++)
                  sampleInputs[ss*nInputs+iInd] = log10(sampleInputs[ss*nInputs+iInd]);
            }
            else
            {
               for (ss = 0; ss < nSamples; ss++)
                  sampleInputs[ss*nInputs+iInd] = pow(10.0,sampleInputs[ss*nInputs+iInd]);
            }
            if (iLowerB != NULL)
            {
               iLowerB[iInd] = PSUADE_UNDEFINED;
               iUpperB[iInd] = -PSUADE_UNDEFINED;
               for (ss = 0; ss < nSamples; ss++)
               {
                  if (sampleInputs[ss*nInputs+iInd] < iLowerB[iInd])
                     iLowerB[iInd] = sampleInputs[ss*nInputs+iInd];
                  if (sampleInputs[ss*nInputs+iInd] > iUpperB[iInd])
                     iUpperB[iInd] = sampleInputs[ss*nInputs+iInd];
               }
            }
            psuadeIO_->updateInputSection(nSamples,nInputs,NULL,NULL,NULL,
                                          sampleInputs,NULL); 
            printf("Input transformation completed: use 'write' to store.\n");
         }
      }

      // transform an output to log or power
      // ----------------------------------------------------------------
      else if (!strcmp(command, "otran"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("otran: transform an output to log10 or power of 10\n");
            printf("syntax: otran (no argument needed)\n");
            continue;
         }
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no data (load data first).\n");
         }
         else
         {
            sprintf(pString,"Enter output number (1 - %d) : ", nOutputs);
            oInd = getInt(1, nOutputs, pString);
            oInd--;
            sprintf(pString,"Enter type: (1) log10(output), (2) 10^(output) : ");
            kk = getInt(1, 2, pString);
            if (kk == 1)
            {
               for (ss = 0; ss < nSamples; ss++)
                  sampleOutputs[ss*nOutputs+oInd] = 
                           log10(sampleOutputs[ss*nOutputs+oInd]);
            }
            else
            {
               for (ss = 0; ss < nSamples; ss++)
                  sampleOutputs[ss*nOutputs+oInd] = 
                           pow(10.0,sampleOutputs[ss*nOutputs+oInd]);
            }
            psuadeIO_->updateOutputSection(nSamples,nOutputs,
                                           sampleOutputs, sampleStates, NULL); 
            printf("Output transformation completed: use 'write' to store.\n");
         }
      }

      // trace input and output values for each run
      // ----------------------------------------------------------------
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
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: no data (load data first).\n");
         }
         else
         {
            sprintf(pString,"Number of inputs to plot (1 - %d) : ", nInputs);
            kk = getInt(1, nInputs, pString);
            if (kk < nInputs)
            {
               indSet = new int[kk];
               for (ii = 0; ii < kk; ii++)
               {
                  sprintf(pString,"Enter the %d-th input (1 - %d): ",
                          ii+1, nInputs);
                  indSet[ii] = getInt(1,nInputs,pString);
                  indSet[ii]--;
               }
            }
            else
            {
               kk = nInputs;
               indSet = new int[nInputs];
               for (ii = 0; ii < nInputs; ii++) indSet[ii] = ii;
            }
            sprintf(pString,"Number of outputs to plot (1 - %d) : ",nOutputs);
            ll = getInt(1, nOutputs, pString);
            if (ll < nOutputs)
            {
               tempI = new int[ll];
               for (ii = 0; ii < ll; ii++)
               {
                  sprintf(pString,"Enter the %d-th output (1 - %d): ",
                          ii+1, nOutputs);
                  tempI[ii] = getInt(1,nOutputs,pString);
                  tempI[ii]--;
               }
            }
            else
            {
               ll = nOutputs;
               tempI = new int[nOutputs];
               for (ii = 0; ii < nOutputs; ii++) tempI[ii] = ii;
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
            tempV = new double[nOutputs];
            tempW = new double[nOutputs];
            for (ii = 0; ii < nOutputs; ii++)
            {
               tempV[ii] =  PSUADE_UNDEFINED;
               tempW[ii] = -PSUADE_UNDEFINED;
               for (ss = 0; ss < nSamples; ss++)
               {
                  ddata = sampleOutputs[ss*nOutputs+ii];
                  tempV[ii] = (ddata < tempV[ii]) ? ddata : tempV[ii];
                  tempW[ii] = (ddata > tempW[ii]) ? ddata : tempW[ii];
               }
            }
            fprintf(fp, "Y = 1:%d;\n", kk+ll); 
            fprintf(fp, "X = [\n"); 
            for (ss = 0; ss < nSamples; ss++)
            {
               for (ii = 0; ii < kk; ii++)
               {
                  ddata = (sampleInputs[ss*nInputs+indSet[ii]] - iLowerB[indSet[ii]]) /
                          (iUpperB[indSet[ii]] - iLowerB[indSet[ii]]);
                  fprintf(fp, "%e ", ddata);
               }
               for (ii = 0; ii < ll; ii++)
               {
                  ddata = (sampleOutputs[ss*nOutputs+tempI[ii]] - tempV[tempI[ii]])/
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
      }

      // plot all pairs of inputs
      // ----------------------------------------------------------------
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
         }
         else
         {
            sprintf(pString,"Number of inputs to plot (1 - %d) : ", nInputs);
            kk = getInt(1, nInputs, pString);
            if (kk < nInputs)
            {
               indSet = new int[kk];
               for (ii = 0; ii < kk; ii++)
               {
                  sprintf(pString,"Enter the %d-th input (1 - %d): ",
                          ii+1, nInputs);
                  indSet[ii] = getInt(1,nInputs,pString);
                  indSet[ii]--;
               }
            }
            else
            {
               kk = nInputs;
               indSet = new int[nInputs];
               for (ii = 0; ii < nInputs; ii++) indSet[ii] = ii;
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
            if (psPlotTool_ == 1)
            {
               fprintf(fp,"// plotMode  = 0 : plot all in a single plot\n");
               fprintf(fp,"// plotMode ~= 0 : plot one at a time\n");
               fprintf(fp,"// ranflag: to distinguish between identical points\n");
               fprintf(fp,"//         by adding a small perturbation (when on)\n");
            }
            else 
            {
               fprintf(fp,"%% plotMode  = 0 : plot all in a single plot\n");
               fprintf(fp,"%% plotMode ~= 0 : plot one at a time\n");
               fprintf(fp,"%% ranflag: to distinguish between identical points\n");
               fprintf(fp,"%%         by adding a small perturbation (when on)\n");
            }
            fprintf(fp,"plotMode = 0;\n");
            fprintf(fp,"ranflag  = 0;\n");
            for (ii = 0; ii < kk; ii++)
            {
               fprintf(fp, "X%d = [\n", ii+1); 
               for (jj = 0; jj < nSamples; jj++)
                  fprintf(fp, "%e\n", sampleInputs[jj*nInputs+indSet[ii]]);
               fprintf(fp, "];\n"); 
            }
            fprintf(fp, "S = [\n"); 
            for (jj = 0; jj < nSamples; jj++)
               fprintf(fp, "%d\n", sampleStates[jj]);
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
                     fprintf(fp,"plot(X(:,%d),X(:,%d),'x');", jj+1, ii+1);
                  else
                  {
                     fprintf(fp, "iset = find(X(:,%d)==0);\n",kk+1);
                     fprintf(fp, "plot(X(iset,%d).*(1+ranflag*", jj+1);
                     fprintf(fp, "rand(size(iset,1),1)/100),X(iset,%d)",
                             ii+1);
                     fprintf(fp, ".*(1+ranflag*rand(size(iset,1),1)/100),'rX',");
                     fprintf(fp, "'MarkerSize',13)\n");
                     if (psPlotTool_ == 1)
                          fprintf(fp, "set(gca(),\"auto_clear\",\"off\")\n");
                     else fprintf(fp, "hold on\n");
                     fprintf(fp, "iset = find(X(:,%d)~=0);\n",kk+1);
                     fprintf(fp, "plot(X(iset,%d).*(1+ranflag*",jj+1);
                     fprintf(fp, "rand(size(iset,1),1)/100),X(iset,%d)",
                             ii+1);
                     fprintf(fp, ".*(1+ranflag*rand(size(iset,1),1)/100),'b*',");
                     fprintf(fp, "'MarkerSize',13)\n");
                  }
                  fwritePlotXLabel(fp, inputNames[indSet[jj]]);
                  fwritePlotYLabel(fp, inputNames[indSet[ii]]);
                  fwritePlotAxes(fp);
                  if (psPlotTool_ == 1)
                  {
                     fprintf(fp, "set(gca(),\"auto_clear\",\"on\")\n");
                     fprintf(fp, "a=gca();\n");
                     fprintf(fp, "a.data_bounds=[%e,%e;%e,%e];\n",
                             iLowerB[indSet[jj]], iLowerB[indSet[ii]],
                             iUpperB[indSet[jj]], iUpperB[indSet[ii]]);
                  }
                  else
                  {
                     fprintf(fp, "hold off\n");
                     fprintf(fp,"axis([%e %e %e %e])\n",
                             iLowerB[indSet[jj]], iUpperB[indSet[jj]],
                             iLowerB[indSet[ii]], iUpperB[indSet[ii]]);
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
      }

      // plot input pdf from sample
      // ----------------------------------------------------------------
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
         }
         else
         {
            sprintf(pString,"Number of inputs to plot (1 - %d) : ", nInputs);
            kk = getInt(1, nInputs, pString);
            if (kk < nInputs)
            {
               indSet = new int[kk];
               for (ii = 0; ii < kk; ii++)
               {
                  sprintf(pString,"Enter the %d-th input (1 - %d): ",
                          ii+1, nInputs);
                  indSet[ii] = getInt(1,nInputs,pString);
                  indSet[ii]--;
               }
            }
            else
            {
               kk = nInputs;
               indSet = new int[nInputs];
               for (ii = 0; ii < nInputs; ii++) indSet[ii] = ii;
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
            for (ii = 0; ii < nSamples; ii++)
            {
               for (jj = 0; jj < kk; jj++)
                  fprintf(fp, "%e ", sampleInputs[ii*nInputs+indSet[jj]]);
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
                  if (inputNames != NULL)
                       fwritePlotXLabel(fp, inputNames[indSet[ii]]);
                  else fwritePlotXLabel(fp, "Input Value");
                  fwritePlotYLabel(fp, "Probabilities");
               }
               else
               {
                  fprintf(fp, "X = X%d;\n",ii+1);
                  fprintf(fp, "[nk,xk]=hist(X,10);\n");
                  fprintf(fp, "bar(xk,nk/%d,1.0)\n",nSamples);
                  fwritePlotAxes(fp);
                  fwritePlotTitle(fp, "Input Distribution");
                  if (inputNames != NULL)
                     fwritePlotXLabel(fp, inputNames[indSet[ii]]);
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
      }

      // set diagnostics output level
      // ----------------------------------------------------------------
      else if (!strcmp(command, "printlevel"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("printlevel: set printlevel (0 - 5)\n");
            continue;
         }
         sscanf(lineIn, "%s %d", command, &ii);
         if (ii >= -1 && ii <= 5) outputLevel_ = ii;
         printf("psuade printlevel set to %d\n", ii);
      }

      // set diagnostics output level
      // ----------------------------------------------------------------
      else if (!strcmp(command, "on"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("on: set printlevel to 4 and turn on interactive mode\n");
            continue;
         }
         outputLevel_ = 4;
         psAnalysisInteractive_ = 1;
      }

      // screening of parameters with layered GP
      // ----------------------------------------------------------------
      else if (!strcmp(command, "gp_sa2") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("gp_sa2: GP SA on sub-divisions of each input\n");
            printf("        (not verified yet)\n");
            continue;
         }
         printf("ERROR: no data (load data first).\n");
      }
      else if (!strcmp(command, "gp_sa2") && nInputs > 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("gp_sa2: GP SA on sub-divisions of each input\n");
            printf("        (not verified yet)\n");
            continue;
         }
         if (psuadeIO_ == NULL) printf("ERROR: data not loaded yet.\n");
         else
         {
            outputID = 0;
            if (nOutputs > 1)
            {
               sprintf(pString, "Enter output number (1 - %d) = ", nOutputs);
               outputID = getInt(1, nOutputs, pString);
               outputID--;
            }
            faType = PSUADE_RS_GP1;
            printf("This command may be very time-consuming to execute.\n");
            count = nSamples / 100;
            if (count > 4) count = 4;
            sprintf(pString, "How many partitions (1 - %d) = ", count);
            nParts = getInt(1, count, pString);
            count = nSamples / nParts + 1;
            tempX = new double[count*nInputs];
            tempY = new double[count];
            tempW = new double[nInputs];
            for (ii = 0; ii < nInputs; ii++)
            {
               tempW[ii] = 0.0;
               ddata = (iUpperB[ii] - iLowerB[ii]) / (double) nParts;
               for (jj = 0; jj < nParts; jj++)
               {
                  count = 0;
                  for (kk = 0; kk < nSamples; kk++)
                  {
                     dtemp = sampleInputs[kk*nInputs+ii];
                     if (dtemp >= ddata*jj+iLowerB[ii] &&
                         dtemp <= ddata*(jj+1)+iLowerB[ii])
                     {
                        for (ind = 0; ind < nInputs; ind++)
                           tempX[count*nInputs+ind] = 
                                 sampleInputs[kk*nInputs+ii];
                        tempY[count++] = sampleOutputs[kk*nOutputs+outputID];
                     }
                  }
                  if (count == 0)
                  {
                     printf("ERROR encountered.\n");
                     exit(1);
                  }
                  if (outputLevel_ > 1)
                    printf("Input %3d Part %3d has size %d\n",ii+1,jj+1,count);
                  faPtr = genFA(faType, nInputs, iOne, count);
                  if (faPtr != NULL)
                  {
                     faPtr->setBounds(iLowerB, iUpperB);
                     faPtr->setOutputLevel(outputLevel_);
                     iInd = -999;
                     faPtr->genNDGridData(tempX,tempY,&iInd,NULL,NULL);
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
      }

      // sum-of-trees screening of parameters using boosting
      // ----------------------------------------------------------------
      else if (!strcmp(command, "sot_sa2") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("sot_sa2: sum-of-trees screening using boosting\n");
            printf("syntax: sot_sa2 (no argument needed)\n");
            continue;
         }
         printf("ERROR: no data (load data first).\n");
      }
      else if (!strcmp(command, "sot_sa2") && nInputs > 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("sot_sa2: sum-of-trees screening using boosting\n");
            printf("syntax: sot_sa2 (no argument needed)\n");
            continue;
         }
         if (psuadeIO_ == NULL) printf("ERROR: data not loaded yet.\n");
         else
         {
            outputID = 0;
            if (nOutputs > 1)
            {
               sprintf(pString, "Enter output number (1 - %d) = ", nOutputs);
               outputID = getInt(1, nOutputs, pString);
               outputID--;
            }
            faType = PSUADE_RS_SOTS;
            faPtr  = genFA(faType, nInputs, iOne, nSamples);
            if (faPtr != NULL)
            {
               faPtr->setBounds(iLowerB, iUpperB);
               faPtr->setOutputLevel(outputLevel_);
               tempY = new double[nSamples];
               for (ii = 0; ii < nSamples; ii++)
                  tempY[ii] = sampleOutputs[ii*nOutputs+outputID];
               kk = -999;
               status = faPtr->genNDGridData(sampleInputs,tempY,&kk,NULL,NULL);
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
      }

      // Delta test for variable selection
      // ----------------------------------------------------------------
      else if (!strcmp(command, "delta_test") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("delta_test: perform Delta test\n");
            printf("syntax: delta_test (no argument needed)\n");
            continue;
         }
         printf("ERROR: no data (load data first).\n");
      }
      else if (!strcmp(command, "delta_test") && nInputs > 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("delta_test: perform Delta test\n");
            printf("syntax: delta_test (no argument needed)\n");
            continue;
         }
         if (psuadeIO_ == NULL) printf("ERROR: data file not loaded.\n");
         else
         {
            outputID = 0;
            if (nOutputs > 1)
            {
               sprintf(pString, "Enter output number (1 - %d) = ", nOutputs);
               outputID = getInt(1, nOutputs, pString);
               outputID--;
            }
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
      }

      // Eta test for variable selection
      // ----------------------------------------------------------------
      else if (!strcmp(command, "eta_test") && nInputs <= 0)
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
         printf("ERROR: no data (load data first).\n");
      }
      else if (!strcmp(command, "eta_test") && nInputs > 0)
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
         if (psuadeIO_ == NULL) printf("ERROR: data file not loaded.\n");
         else
         {
            outputID = 0;
            if (nOutputs > 1)
            {
               sprintf(pString, "Enter output number (1 - %d) = ", nOutputs);
               outputID = getInt(1, nOutputs, pString);
               outputID--;
            }
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
      }

      // Gower distance analysis
      // ----------------------------------------------------------------
      else if (!strcmp(command, "gd_test") && nInputs <= 0)
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
         printf("ERROR: no data (load data first).\n");
      }
      else if (!strcmp(command, "gd_test") && nInputs > 0)
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
         if (psuadeIO_ == NULL) printf("ERROR: data file not loaded.\n");
         else
         {
            outputID = 0;
            if (nOutputs > 1)
            {
               sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
               outputID = getInt(1, nOutputs, pString);
               outputID--;
               fgets(lineIn,500,stdin); 
            }
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
      }

      // Sobol main effect
      // ----------------------------------------------------------------
      else if (!strcmp(command, "rssobol1") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rssobol1: RS-based Sobol' first order indices\n");
            printf("syntax: rssobol1 (no argument needed)\n");
            printf("note: to compute error bars for indices, turn on\n");
            printf("      ana_expert mode and set ntimes to 100.\n");
            continue;
         }
         printf("ERROR: no data (load data first).\n");
      }
      else if (!strcmp(command, "rssobol1") && nInputs > 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rssobol1: RS-based Sobol' first order indices\n");
            printf("syntax: rssobol1 (no argument needed)\n");
            printf("note: to compute error bars for indices, turn on\n");
            printf("      ana_expert mode and set ntimes to 100.\n");
            continue;
         }
         if (psuadeIO_ == NULL) printf("ERROR: data file not loaded.\n");
         else
         {
            outputID = 0;
            if (nOutputs > 1)
            {
               sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
               outputID = getInt(1, nOutputs, pString);
               outputID--;
            }
            analysisMethod = PSUADE_ANA_RSSOBOL1;
            anaManager = new AnalysisManager();
            anaManager->setup(analysisMethod, 0);
            psuadeIO_->getParameter("ana_diagnostics",pPtr);
            ii = pPtr.intData_;
            psuadeIO_->updateAnalysisSection(-1,-1,-1,-2,-1,-1);
            anaManager->analyze(psuadeIO_, 0, NULL, outputID);
            psuadeIO_->updateAnalysisSection(-1,-1,-1,ii,-1,-1);
            pdata = psuadeIO_->getAuxData();
            if (pdata->nDbles_ >= nInputs)
            {
               printEquals(0);
               printf("Main Effect Statistics: \n");
               if (pdata->dbleData_ > 0)
               {
                  for (ii = 0; ii < nInputs; ii++)
                  {
                     printf("Input %4d: Sobol' main effect = %12.4e",
                            ii+1,pdata->dbleArray_[ii]);
                     if (pdata->nDbles_ == 3*nInputs)
                          printf(", bounds = [%12.4e, %12.4e]\n",
                                 pdata->dbleArray_[ii+nInputs],
                                 pdata->dbleArray_[ii+2*nInputs]);
                     else printf(" (unnormalized)\n");
                  } 
                  if (psPlotTool_ == 1)
                  {
                     printf("rssobol1: scilab graphics is not currently available.\n");
                     fp = NULL;
                  }
                  else fp = fopen("matlabrssobol1.m", "w");
                  if (fp != NULL)
                  {
                     fprintf(fp, "%% This file contains Sobol' indices\n");
                     fprintf(fp, "%% set sortFlag = 1 and set nn to be the number\n");
                     fprintf(fp, "%% of inputs to display.\n");
                     fprintf(fp, "sortFlag = 0;\n");
                     fprintf(fp, "nn = %d;\n", nInputs);
                     fprintf(fp, "Mids = [\n");
                     for (ii = 0; ii < nInputs; ii++) 
                        fprintf(fp,"%24.16e\n", 
                                pdata->dbleArray_[ii]);
                     fprintf(fp, "];\n");
                     if (pdata->nDbles_ == 3*nInputs)
                     {
                        fprintf(fp, "Mins = [\n");
                        for (ii = 0; ii < nInputs; ii++) 
                           fprintf(fp,"%24.16e\n", 
                                   pdata->dbleArray_[nInputs+ii]);
                        fprintf(fp, "];\n");
                        fprintf(fp, "Maxs = [\n");
                        for (ii = 0; ii < nInputs; ii++) 
                           fprintf(fp,"%24.16e\n", 
                                   pdata->dbleArray_[2*nInputs+ii]);
                        fprintf(fp, "];\n");
                     }
                     if (inputNames == NULL)
                     {
                        fprintf(fp, "  Str = {");
                        for (ii = 0; ii < nInputs-1; ii++) fprintf(fp,"'X%d',",ii+1);
                        fprintf(fp,"'X%d'};\n",nInputs);
                     }
                     else
                     {
                        fprintf(fp, "  Str = {");
                        for (ii = 0; ii < nInputs-1; ii++)
                        {
                           if (inputNames[ii] != NULL) 
                                fprintf(fp,"'%s',",inputNames[ii]);
                           else fprintf(fp,"'X%d',",ii+1);
                        }
                        if (inputNames[nInputs-1] != NULL)
                             fprintf(fp,"'%s'};\n",inputNames[nInputs-1]);
                        else fprintf(fp,"'X%d'};\n",nInputs);
                     }
                     fprintf(fp, "clf\n");
                     fprintf(fp, "if (sortFlag == 1)\n");
                     fprintf(fp, "  [Mids, I2] = sort(Mids,'descend');\n");
                     if (pdata->nDbles_ == 3*nInputs)
                     {
                        fprintf(fp, "  Maxs = Maxs(I2);\n");
                        fprintf(fp, "  Mins = Mins(I2);\n");
                     }
                     fprintf(fp, "  Str  = Str(I2);\n");
                     fprintf(fp, "  I2 = I2(1:nn);\n");
                     fprintf(fp, "  Mids = Mids(1:nn);\n");
                     if (pdata->nDbles_ == 3*nInputs)
                     {
                        fprintf(fp, "  Maxs = Maxs(1:nn);\n");
                        fprintf(fp, "  Mins = Mins(1:nn);\n");
                     }
                     fprintf(fp, "  Str  = Str(1:nn);\n");
                     fprintf(fp, "end\n");
                     if (pdata->nDbles_ == 3*nInputs)
                     {
                        fprintf(fp, "ymin = min(Mins);\n");
                        fprintf(fp, "ymax = max(Maxs);\n");
                     }
                     else
                     {
                        fprintf(fp, "ymin = min(Mids);\n");
                        fprintf(fp, "ymax = max(Mids);\n");
                     }
                     fprintf(fp, "h2 = 0.05 * (ymax - ymin);\n");
                     fprintf(fp, "bar(Mids,0.8);\n");
                     if (pdata->nDbles_ == 3*nInputs)
                     {
                        fprintf(fp, "for ii = 1:nn\n");
                        fprintf(fp, "%% h = plot(ii,Means(ii),'r*','MarkerSize',13);\n");
                        fprintf(fp, "   if (ii == 1)\n");
                        fprintf(fp, "      hold on\n");
                        fprintf(fp, "   end;\n");
                        fprintf(fp, "   XX = [ii ii];\n");
                        fprintf(fp, "   YY = [Mins(ii) Maxs(ii)];\n");
                        fprintf(fp, "   plot(XX,YY,'-ko','LineWidth',3.0,'MarkerEdgeColor',");
                        fprintf(fp, "'k','MarkerFaceColor','g','MarkerSize',13)\n");
                        fprintf(fp, "end;\n");
                     }
                     fprintf(fp, "ymin=0;\n");
                     fprintf(fp, "axis([0 nn+1 ymin ymax])\n");
                     fwritePlotAxes(fp);
                     fprintf(fp, "set(gca,'XTickLabel',[]);\n");
                     fprintf(fp, "th=text(1:nn, repmat(ymin-0.07*(ymax-ymin),nn,1),Str,");
                     fprintf(fp, "'HorizontalAlignment','left','rotation',90);\n");
                     fprintf(fp, "set(th, 'fontsize', 12)\n");
                     fprintf(fp, "set(th, 'fontweight', 'bold')\n");
                     fwritePlotTitle(fp,"Sobol First Order Indices");
                     fwritePlotYLabel(fp, "Sobol Indices");
                     fprintf(fp, "hold off\n");
                     fclose(fp);
                     printf("rssobol1 plot matlab file = matlabrssobol1.m\n");
                  }
               }
               else
               {
                  printf("Total variance = 0. Hence, no main effect plot.\n");
               }
               pdata->clean();
            }
            delete anaManager;
         }
      }

      // rssobol1 with bootstrap
      // ----------------------------------------------------------------
      else if (!strcmp(command, "rs_qsa"))
      {
         printf("INFO: this command has been splitted into rssobol1b\n");
         printf("      and rssoboltsib.\n");
      }
      else if (!strcmp(command, "rssobol1b") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rssobol1b: RS-based Sobol' sensitivity analysis\n");
            printf("syntax: rssobol1b (no argument needed)\n");
            printf("Note: This command differs from rssobol1 and 'me'\n");
            printf("      in that it uses bootstrapped samples multiple\n");
            printf("      times to get the errors in Sobol' indices due\n");
            printf("      response surface errors.\n");
            continue;
         }
         printf("ERROR: no data (load data first).\n");
      }
      else if (!strcmp(command, "rssobol1b") && nInputs > 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rssobol1b: RS-based Sobol' sensitivity analysis\n");
            printf("syntax: rssobol1b (no argument needed)\n");
            printf("Note: This command differs from rssobol1 and 'me'\n");
            printf("      in that it uses bootstrapped samples multiple\n");
            printf("      times to get the errors in Sobol' indices due\n");
            printf("      response surface errors.\n");
         }
         if (psuadeIO_ == NULL) printf("ERROR: data file not loaded.\n");
         else
         {
            if (nSamples < 5)
            {
               printf("INFO: This command is not suitable for small samples.\n");
               printf("      nSamples needs to be at least 5.\n");
               continue;
            }
            printf("INFO: MAKE SURE YOU SET THE RESPONSE SURFACE TYPE IN YOUR\n");
            printf("      DATA FILE.\n");
            outputID = 0;
            if (nOutputs > 1)
            {
               sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
               outputID = getInt(1, nOutputs, pString);
               outputID--;
            }

            analysisMethod = PSUADE_ANA_RSSOBOL1;
            anaManager = new AnalysisManager();
            anaManager->setup(analysisMethod, 0);
 
            psuadeIO_->getParameter("ana_diagnostics",pPtr);
            saveMode1 = psAnaExpertMode_;
            psAnaExpertMode_ = 0;
            saveMode2 = psAnalysisInteractive_;
            psAnalysisInteractive_ = 0;

            tempX = new double[nSamples*nInputs];
            tempY = new double[nSamples];
            tempI = new int[nSamples];
            states = new int[nSamples];

            sprintf(pString, "How many bootstrapped samples to use (10 - 1000) : ");
            count = getInt(10, 1000, pString);
            tempT = new double[count*nInputs];

            printf("NOTE: the response surface type will be obtained from you data\n");
            printf("      file. If you would like to fine-tune your response surface,\n");
            printf("      do not turn on rs_expert mode, set the parameters in the\n");
            printf("      config file.\n");
            for (kk = 0; kk < count; kk++)
            {
               printf("rssobol1b: iteration %d\n", kk+1);
               for (jj = 0; jj < nSamples; jj++) tempI[jj] = 0;
               ss = nSamples2 = 0;
               while (ss < nSamples)
               {
                  ind = PSUADE_rand() % nSamples;
                  if (tempI[ind] == 0)
                  {
                     for (ii = 0; ii < nInputs; ii++)
                        tempX[nSamples2*nInputs+ii] = sampleInputs[ind*nInputs+ii];
                     tempY[nSamples2] = sampleOutputs[ind*nOutputs+outputID];
                     states[nSamples2] = sampleStates[ind];
                     tempI[ind] = 1;
                     nSamples2++;
                  }
                  ss++;
               }
               psuadeIO_->updateInputSection(nSamples2,nInputs,NULL,NULL,
                                             NULL,tempX,NULL); 
               psuadeIO_->updateOutputSection(nSamples2,1,tempY,states,
                                              outputNames);
               psuadeIO_->updateMethodSection(PSUADE_SAMP_MC,nSamples2,
                                              -1,-1,-1);

               anaManager->analyze(psuadeIO_, 0, NULL, 0);
               pdata = psuadeIO_->getAuxData(); 
               if (pdata->nDbles_ != nInputs)
               {
                  printf("ERROR: nInputs do not match (%d, %d).\n",
                         pdata->nDbles_, nInputs);
                  printf("       Consult PSUADE developers.\n");
                  delete [] tempT;
                  delete [] tempX;
                  delete [] tempY;
                  delete [] tempI;
                  delete [] states;
                  continue;
               }

               if (pdata->dbleData_ > 0)
                  for (ii = 0; ii < nInputs; ii++)
                     tempT[kk*nInputs+ii] = 
                          pdata->dbleArray_[ii]/pdata->dbleData_;
               else
                  for (ii = 0; ii < nInputs; ii++)
                     tempT[kk*nInputs+ii] = pdata->dbleArray_[ii];

               pdata->clean();
            }
            delete [] tempX;
            delete [] tempY;
            delete [] tempI;
            delete [] states;
            tempM = new double[nInputs];
            for (ii = 0; ii < nInputs; ii++)
            {
               tempM[ii] = tempT[ii];
               for (jj = 1; jj < count; jj++) tempM[ii] += tempT[jj*nInputs+ii];
               tempM[ii] /= (double) count;
            }
            tempV = new double[nInputs];
            for (ii = 0; ii < nInputs; ii++)
            {
               tempV[ii] = pow(tempT[ii]-tempM[ii], 2.0);
               for (jj = 1; jj < count; jj++) 
                  tempV[ii] += pow(tempT[jj*nInputs+ii]-tempM[ii],2.0);
               tempV[ii] /= (double) (count - 1);
               tempV[ii] = sqrt(tempV[ii]);
            }
            printEquals(0);
            printf("Sobol1 Statistics (based on %d replications): \n",
                   count);
            for (ii = 0; ii < nInputs; ii++)
               printf("Input %4d: mean = %16.8e, std = %16.8e\n",ii+1,
                      tempM[ii],tempV[ii]);

            if (psPlotTool_ == 1)
            {
               fp = fopen("scilabrssobol1b.sci","w");
               fprintf(fp, "// This file contains first order Sobol' indices\n");
               fprintf(fp, "// with error bars coming from bootstrapping.\n");
               fprintf(fp, "// to select the most important ones to display,\n");
               fprintf(fp, "// set sortFlag = 1 and set nn to be the number\n");
               fprintf(fp, "// of inputs to display.\n");
            }
            else
            {
               fp = fopen("matlabrssobol1b.m","w");
               fprintf(fp, "%% This file contains first order Sobol' indices\n");
               fprintf(fp, "%% with error bars coming from bootstrapping.\n");
               fprintf(fp, "%% to select the most important ones to display,\n");
               fprintf(fp, "%% set sortFlag = 1 and set nn to be the number\n");
               fprintf(fp, "%% of inputs to display.\n");
            }
            fprintf(fp, "sortFlag = 0;\n");
            fprintf(fp, "nn = %d;\n", nInputs);
            fprintf(fp, "Means = [\n");
            for (ii = 0; ii < nInputs; ii++) fprintf(fp,"%24.16e\n",tempM[ii]);
            fprintf(fp, "];\n");
            fprintf(fp, "Stds = [\n");
            for (ii = 0; ii < nInputs; ii++) fprintf(fp,"%24.16e\n",tempV[ii]);
            fprintf(fp, "];\n");
            if (inputNames == NULL)
            {
               fprintf(fp, "  Str = {");
               for (ii = 0; ii < nInputs-1; ii++) fprintf(fp,"'X%d',",ii+1);
               fprintf(fp,"'X%d'};\n",nInputs);
            }
            else
            {
               fprintf(fp, "  Str = {");
               for (ii = 0; ii < nInputs-1; ii++)
               {
                  if (inputNames[ii] != NULL) 
                       fprintf(fp,"'%s',",inputNames[ii]);
                  else fprintf(fp,"'X%d',",ii+1);
               }
               if (inputNames[nInputs-1] != NULL) 
                    fprintf(fp,"'%s'};\n",inputNames[nInputs-1]);
               else fprintf(fp,"'X%d'};\n",nInputs);
            }
            fwritePlotCLF(fp);
            fprintf(fp, "if (sortFlag == 1)\n");
            if (psPlotTool_ == 1)
                 fprintf(fp, "  [Means, I2] = gsort(Means);\n");
            else fprintf(fp, "  [Means, I2] = sort(Means,'descend');\n");
            fprintf(fp, "  Stds = Stds(I2);\n");
            fprintf(fp, "  I2 = I2(1:nn);\n");
            fprintf(fp, "  Means = Means(1:nn);\n");
            fprintf(fp, "  Stds = Stds(1:nn);\n");
            fprintf(fp, "  Str  = Str(I2);\n");
            fprintf(fp, "end\n");
            fprintf(fp, "ymin = min(Means-Stds);\n");
            fprintf(fp, "ymax = max(Means+Stds);\n");
            fprintf(fp, "h2 = 0.05 * (ymax - ymin);\n");
            if (psPlotTool_ == 1) fprintf(fp, "drawlater\n");
            fprintf(fp, "bar(Means,0.8);\n");
            fprintf(fp, "for ii = 1:nn\n");
            fprintf(fp, "   if (ii == 1)\n");
            if (psPlotTool_ == 1)
                 fprintf(fp, "      set(gca(),\"auto_clear\",\"off\")\n");
            else fprintf(fp, "      hold on\n");
            fprintf(fp, "   end;\n");
            fprintf(fp, "   XX = [ii ii];\n");
            fprintf(fp, "   YY = [Means(ii)-Stds(ii) Means(ii)+Stds(ii)];\n");
            fprintf(fp, "   plot(XX,YY,'-ko','LineWidth',3.0,'MarkerEdgeColor',");
            fprintf(fp, "'k','MarkerFaceColor','g','MarkerSize',13)\n");
            fprintf(fp, "end;\n");
            fwritePlotAxes(fp);
            if (psPlotTool_ == 1)
            {
               fprintf(fp, "a=gca();\n");
               fprintf(fp, "a.data_bounds=[0, ymin; nn+1, ymax];\n");
               fprintf(fp, "newtick = a.x_ticks;\n");
               fprintf(fp, "newtick(2) = [1:nn]';\n");
               fprintf(fp, "newtick(3) = Str';\n");
               fprintf(fp, "a.x_ticks = newtick;\n");
               fprintf(fp, "a.x_label.font_size = 3;\n");
               fprintf(fp, "a.x_label.font_style = 4;\n");
            }
            else
            {
               fprintf(fp, "axis([0  nn+1 ymin ymax])\n");
               fprintf(fp, "set(gca,'XTickLabel',[]);\n");
               fprintf(fp, "th=text(1:nn, repmat(ymin-0.05*(ymax-ymin),nn,1),Str,");
               fprintf(fp, "'HorizontalAlignment','left','rotation',90);\n");
               fprintf(fp, "set(th, 'fontsize', 12)\n");
               fprintf(fp, "set(th, 'fontweight', 'bold')\n");
            }
            fwritePlotTitle(fp,"First Order Sobol Indices (bootstrap)");
            fwritePlotYLabel(fp, "First Order Sobol Index");
            if (psPlotTool_ == 1)
            {
               fprintf(fp, "drawnow\n");
               printf("Scilab plot for main effect is in scilabrssobol1b.sci.\n");
            }
            else printf("Matlab plot for main effect is in matlabrssobol1b.m.\n");
            fclose(fp);

            delete anaManager;
            delete [] tempT;
            delete [] tempM;
            delete [] tempV;
            tempX = tempY = tempM = tempT = tempV = NULL;
            tempI = NULL;

            psAnaExpertMode_ = saveMode1;
            psAnalysisInteractive_ = saveMode2;
            cleanUp();
            if (inputNames != NULL)
            {
               for (ii = 0; ii < nInputs; ii++)
                  if (inputNames[ii] != NULL) delete [] inputNames[ii];
               delete [] inputNames;
            }
            inputNames = NULL;
            if (outputNames != NULL)
            {
               for (ii = 0; ii < nOutputs; ii++)
                  if (outputNames[ii] != NULL) delete [] outputNames[ii];
               delete [] outputNames;
            }
            outputNames = NULL;
            if (sampleInputs != NULL) delete [] sampleInputs;
            if (sampleOutputs != NULL) delete [] sampleOutputs;
            if (sampleStates != NULL) delete [] sampleStates;
            sampleInputs = NULL;
            sampleOutputs = NULL;
            sampleStates = NULL;
            nInputs = 0;
            nOutputs = 0;
            printf("LOADED DATA IN LOCAL MEMORY HAVE BEEN ERASED.\n");
            printf("RE-LOAD THE ORIGINAL DATA FILE FOR MORE ANALYSIS.\n");
         }
      }

      // rssoboltsi with bootstrap
      // ----------------------------------------------------------------
      else if (!strcmp(command, "rssoboltsib") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rssoboltsib: RS-based Sobol' total sensitivity analysis\n");
            printf("syntax: rssoboltsib (no argument needed)\n");
            printf("Note: This command differs from rssoboltsi 'tsi'\n");
            printf("      in that it uses bootstrapped samples multiple\n");
            printf("      times to get the errors in Sobol' indices due\n");
            printf("      response surface errors.\n");
            continue;
         }
         printf("ERROR: no data (load data first).\n");
      }
      else if (!strcmp(command, "rssoboltsib") && nInputs > 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rssoboltsib: RS-based Sobol' total sensitivity analysis\n");
            printf("syntax: rssoboltsib (no argument needed)\n");
            printf("Note: This command differs from rssoboltsi 'tsi'\n");
            printf("      in that it uses bootstrapped samples multiple\n");
            printf("      times to get the errors in Sobol' indices due\n");
            printf("      response surface errors.\n");
         }
         if (psuadeIO_ == NULL) printf("ERROR: data file not loaded.\n");
         else
         {
            if (nSamples < 5)
            {
               printf("INFO: This command is not suitable for small samples.\n");
               printf("      nSamples needs to be at least 5.\n");
               continue;
            }
            printf("INFO: MAKE SURE YOU SET THE RESPONSE SURFACE TYPE IN YOUR\n");
            printf("      DATA FILE.\n");
            outputID = 0;
            if (nOutputs > 1)
            {
               sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
               outputID = getInt(1, nOutputs, pString);
               outputID--;
            }

            analysisMethod = PSUADE_ANA_RSSOBOLTSI;
            anaManager = new AnalysisManager();
            anaManager->setup(analysisMethod, 0);
 
            saveMode1 = psAnaExpertMode_;
            psAnaExpertMode_ = 0;
            saveMode2 = psAnalysisInteractive_;
            psAnalysisInteractive_ = 0;

            tempX = new double[nSamples*nInputs];
            tempY = new double[nSamples];
            tempI = new int[nSamples];
            states = new int[nSamples];

            sprintf(pString, "How many bootstrapped samples to use (10 - 1000) : ");
            count = getInt(10, 1000, pString);
            tempT = new double[count*nInputs];

            printf("NOTE: the response surface type will be obtained from you data\n");
            printf("      file. If you would like to fine-tune your response surface,\n");
            printf("      do not turn on rs_expert mode, set the parameters in the\n");
            printf("      config file.\n");
            for (kk = 0; kk < count; kk++)
            {
               printf("rssoboltsib: iteration %d\n", kk+1);
               for (jj = 0; jj < nSamples; jj++) tempI[jj] = 0;
               ss = nSamples2 = 0;
               while (ss < nSamples)
               {
                  ind = PSUADE_rand() % nSamples;
                  if (tempI[ind] == 0)
                  {
                     for (ii = 0; ii < nInputs; ii++)
                        tempX[nSamples2*nInputs+ii] = sampleInputs[ind*nInputs+ii];
                     tempY[nSamples2] = sampleOutputs[ind*nOutputs+outputID];
                     states[nSamples2] = sampleStates[ind];
                     tempI[ind] = 1;
                     nSamples2++;
                  }
                  ss++;
               }
               psuadeIO_->updateInputSection(nSamples2,nInputs,NULL,NULL,
                                             NULL,tempX,NULL); 
               psuadeIO_->updateOutputSection(nSamples2,1,tempY,states,
                                              outputNames);
               psuadeIO_->updateMethodSection(PSUADE_SAMP_MC,nSamples2,
                                              -1,-1,-1);

               anaManager->analyze(psuadeIO_, 0, NULL, 0);
               pdata = psuadeIO_->getAuxData(); 
               if (pdata->nDbles_ != nInputs)
               {
                  printf("ERROR: nInputs do not match (%d, %d).\n",
                         pdata->nDbles_, nInputs);
                  printf("       Consult PSUADE developers.\n");
                  delete [] tempX;
                  delete [] tempY;
                  delete [] tempI;
                  delete [] tempT;
                  delete [] states;
                  continue;
               }

               if (pdata->dbleData_ > 0)
                  for (ii = 0; ii < nInputs; ii++)
                     tempT[kk*nInputs+ii] = 
                          pdata->dbleArray_[ii]/pdata->dbleData_;
               else
                  for (ii = 0; ii < nInputs; ii++)
                     tempT[kk*nInputs+ii] = pdata->dbleArray_[ii];

               pdata->clean();
            }
            delete [] tempX;
            delete [] tempY;
            delete [] tempI;
            delete [] states;
            tempM = new double[nInputs];
            for (ii = 0; ii < nInputs; ii++)
            {
               tempM[ii] = tempT[ii];
               for (jj = 1; jj < count; jj++) tempM[ii] += tempT[jj*nInputs+ii];
               tempM[ii] /= (double) count;
            }
            tempV = new double[nInputs];
            for (ii = 0; ii < nInputs; ii++)
            {
               tempV[ii] = pow(tempT[ii]-tempM[ii], 2.0);
               for (jj = 1; jj < count; jj++) 
                  tempV[ii] += pow(tempT[jj*nInputs+ii]-tempM[ii],2.0);
               tempV[ii] /= (double) (count - 1);
               tempV[ii] = sqrt(tempV[ii]);
            }
            printEquals(0);
            printf("Sobol' Total Senstivity Statistics (based on %d replications): \n",
                   count);
            for (ii = 0; ii < nInputs; ii++)
               printf("Input %4d: mean = %16.8e, std = %16.8e\n",ii+1,
                      tempM[ii],tempV[ii]);

            if (psPlotTool_ == 1)
            {
               fp = fopen("scilabrssoboltsib.sci","w");
               fprintf(fp, "// This file contains total order Sobol' indices\n");
               fprintf(fp, "// with error bars coming from bootstrapping.\n");
               fprintf(fp, "// to select the most important ones to display,\n");
               fprintf(fp, "// set sortFlag = 1 and set nn to be the number\n");
               fprintf(fp, "// of inputs to display.\n");
            }
            else
            {
               fp = fopen("matlabrssoboltsib.m","w");
               fprintf(fp, "%% This file contains total order Sobol' indices\n");
               fprintf(fp, "%% with error bars coming from bootstrapping.\n");
               fprintf(fp, "%% to select the most important ones to display,\n");
               fprintf(fp, "%% set sortFlag = 1 and set nn to be the number\n");
               fprintf(fp, "%% of inputs to display.\n");
            }
            fprintf(fp, "sortFlag = 0;\n");
            fprintf(fp, "nn = %d;\n", nInputs);
            fprintf(fp, "Means = [\n");
            for (ii = 0; ii < nInputs; ii++) fprintf(fp,"%24.16e\n",tempM[ii]);
            fprintf(fp, "];\n");
            fprintf(fp, "Stds = [\n");
            for (ii = 0; ii < nInputs; ii++) fprintf(fp,"%24.16e\n",tempV[ii]);
            fprintf(fp, "];\n");
            if (inputNames == NULL)
            {
               fprintf(fp, "  Str = {");
               for (ii = 0; ii < nInputs-1; ii++) fprintf(fp,"'X%d',",ii+1);
               fprintf(fp,"'X%d'};\n",nInputs);
            }
            else
            {
               fprintf(fp, "  Str = {");
               for (ii = 0; ii < nInputs-1; ii++)
               {
                  if (inputNames[ii] != NULL) 
                       fprintf(fp,"'%s',",inputNames[ii]);
                  else fprintf(fp,"'X%d',",ii+1);
               }
               if (inputNames[nInputs-1] != NULL) 
                    fprintf(fp,"'%s'};\n",inputNames[nInputs-1]);
               else fprintf(fp,"'X%d'};\n",nInputs);
            }
            fprintf(fp, "clf\n");
            fprintf(fp, "if (sortFlag == 1)\n");
            if (psPlotTool_ == 1)
                 fprintf(fp, "  [Means, I2] = gsort(Means);\n");
            else fprintf(fp, "  [Means, I2] = sort(Means,'descend');\n");
            fprintf(fp, "  Stds = Stds(I2);\n");
            fprintf(fp, "  I2 = I2(1:nn);\n");
            fprintf(fp, "  Means = Means(1:nn);\n");
            fprintf(fp, "  Stds = Stds(1:nn);\n");
            fprintf(fp, "  Str  = Str(I2);\n");
            fprintf(fp, "end\n");
            fprintf(fp, "ymin = min(Means-Stds);\n");
            fprintf(fp, "ymax = max(Means+Stds);\n");
            fprintf(fp, "h2 = 0.05 * (ymax - ymin);\n");
            if (psPlotTool_ == 1) fprintf(fp, "drawlater\n");
            fprintf(fp, "bar(Means,0.8);\n");
            fprintf(fp, "for ii = 1:nn\n");
            fprintf(fp, "   if (ii == 1)\n");
            if (psPlotTool_ == 1)
                 fprintf(fp, "      set(gca(),\"auto_clear\",\"off\")\n");
            else fprintf(fp, "      hold on\n");
            fprintf(fp, "   end;\n");
            fprintf(fp, "   XX = [ii ii];\n");
            fprintf(fp, "   YY = [Means(ii)-Stds(ii) Means(ii)+Stds(ii)];\n");
            fprintf(fp, "   plot(XX,YY,'-ko','LineWidth',3.0,'MarkerEdgeColor',");
            fprintf(fp, "'k','MarkerFaceColor','g','MarkerSize',12)\n");
            fprintf(fp, "end;\n");
            fwritePlotAxes(fp);
            if (psPlotTool_ == 1)
            {
               fprintf(fp, "a=gca();\n");
               fprintf(fp, "a.data_bounds=[0, ymin; nn+1, ymax];\n");
               fprintf(fp, "newtick = a.x_ticks;\n");
               fprintf(fp, "newtick(2) = [1:nn]';\n");
               fprintf(fp, "newtick(3) = Str';\n");
               fprintf(fp, "a.x_ticks = newtick;\n");
               fprintf(fp, "a.x_label.font_size = 3;\n");
               fprintf(fp, "a.x_label.font_style = 4;\n");
            }
            else
            {
               fprintf(fp, "axis([0  nn+1 ymin ymax])\n");
               fprintf(fp, "set(gca,'XTickLabel',[]);\n");
               fprintf(fp, "th=text(1:nn, repmat(ymin-0.05*(ymax-ymin),nn,1),Str,");
               fprintf(fp, "'HorizontalAlignment','left','rotation',90);\n");
               fprintf(fp, "set(th, 'fontsize', 12)\n");
               fprintf(fp, "set(th, 'fontweight', 'bold')\n");
            }
            fwritePlotTitle(fp,"Total Order Sobol Indices (bootstrap)");
            fwritePlotYLabel(fp, "Total Order Sobol Index");
            if (psPlotTool_ == 1)
            {
               fprintf(fp, "drawnow\n");
               printf("Scilab plot for total effect is in scilabrssobol1tsi.sci.\n");
            }
            else printf("Matlab plot for total effect is in matlabrssobol1tsi.m.\n");
            fclose(fp);

            delete anaManager;
            delete [] tempT;
            delete [] tempM;
            delete [] tempV;
            tempX = tempY = tempM = tempT = tempV = NULL;
            tempI = NULL;

            psAnaExpertMode_ = saveMode1;
            psAnalysisInteractive_ = saveMode2;
            cleanUp();
            if (psuadeIO_ != NULL) delete psuadeIO_;
            psuadeIO_ = NULL;
            if (inputNames != NULL)
            {
               for (ii = 0; ii < nInputs; ii++)
                  if (inputNames[ii] != NULL) delete [] inputNames[ii];
               delete [] inputNames;
            }
            inputNames = NULL;
            if (outputNames != NULL)
            {
               for (ii = 0; ii < nOutputs; ii++)
                  if (outputNames[ii] != NULL) delete [] outputNames[ii];
               delete [] outputNames;
            }
            outputNames = NULL;
            if (sampleInputs != NULL) delete [] sampleInputs;
            if (sampleOutputs != NULL) delete [] sampleOutputs;
            if (sampleStates != NULL) delete [] sampleStates;
            sampleInputs = NULL;
            sampleOutputs = NULL;
            sampleStates = NULL;
            nInputs = 0;
            nOutputs = 0;
            printf("LOADED DATA IN LOCAL MEMORY HAVE BEEN ERASED.\n");
            printf("RE-LOAD THE ORIGINAL DATA FILE FOR MORE ANALYSIS.\n");
         }
      }

      // Sobol interaction effect
      // ----------------------------------------------------------------
      else if (!strcmp(command, "rssobol2") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rssobol2: RS-based 2-input Sobol' decomposition\n");
            printf("syntax: rssobol2 (no argument needed)\n");
            continue;
         }
         printf("ERROR: no data (load data first).\n");
      }
      else if (!strcmp(command, "rssobol2") && nInputs > 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rssobol2: RS-based 2-input Sobol' decomposition\n");
            printf("syntax: rssobol2 (no argument needed)\n");
            continue;
         }
         if (psuadeIO_ == NULL) printf("ERROR: data file not loaded.\n");
         else
         {
            outputID = 0;
            if (nOutputs > 1)
            {
               sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
               outputID = getInt(1, nOutputs, pString);
               outputID--;
            }
            analysisMethod = PSUADE_ANA_RSSOBOL2;
            anaManager = new AnalysisManager();
            anaManager->setup(analysisMethod, 0);
            psuadeIO_->getParameter("ana_diagnostics",pPtr);
            ii = pPtr.intData_;
            psuadeIO_->updateAnalysisSection(-1,-1,-1,outputLevel_,-1,-1);
            anaManager->analyze(psuadeIO_, 0, NULL, outputID);
            psuadeIO_->updateAnalysisSection(-1,-1,-1,ii,-1,-1);
            pdata = psuadeIO_->getAuxData();
            if (pdata->nDbles_ >= nInputs)
            {
               printEquals(0);
               printf("Two-way Interaction Effect Statistics: \n");
               if (pdata->dbleData_ > 0)
               {
                  for (ii = 0; ii < nInputs; ii++)
                  {
                     for (jj = ii+1; jj < nInputs; jj++)
                     {
                        printf("Inputs %4d %4d: Sobol' interaction effect = ",ii+1,jj+1);
                        printf("%12.4e (unnormalized)\n",pdata->dbleArray_[ii*nInputs+jj]);
                     }
                  }
                  if (psPlotTool_ == 1)
                  {
                     fp = fopen("scilabrssobol2.sci", "w");
                     if (fp == NULL)
                     {
                        printf("ERROR opening file.\n");
                        continue;
                     }
                     fprintf(fp, "// This file contains Sobol' second order indices\n");
                     fprintf(fp, "// set sortFlag = 1 and set nn to be the number\n");
                     fprintf(fp, "// of inputs to display.\n");
                  }
                  else
                  {
                     fp = fopen("matlabrssobol2.m", "w");
                     if (fp == NULL)
                     {
                        printf("ERROR opening file.\n");
                        continue;
                     }
                     fprintf(fp, "%% This file contains Sobol' second order indices\n");
                     fprintf(fp, "%% set sortFlag = 1 and set nn to be the number\n");
                     fprintf(fp, "%% of inputs to display.\n");
                  }
                  fprintf(fp, "sortFlag = 0;\n");
                  fprintf(fp, "nn = %d;\n", nInputs);
                  fprintf(fp, "Mids = [\n");
                  for (ii = 0; ii < nInputs*nInputs; ii++) 
                     fprintf(fp,"%24.16e\n", pdata->dbleArray_[ii]);
                  fprintf(fp, "];\n");
                  if (inputNames == NULL)
                  {
                     fprintf(fp, "Str = {");
                     for (ii = 0; ii < nInputs-1; ii++) fprintf(fp,"'X%d',",ii+1);
                     fprintf(fp,"'X%d'};\n",nInputs);
                  }
                  else
                  {
                     fprintf(fp, "Str = {");
                     for (ii = 0; ii < nInputs-1; ii++)
                     {
                        if (inputNames[ii] != NULL) 
                             fprintf(fp,"'%s',",inputNames[ii]);
                        else fprintf(fp,"'X%d',",ii+1);
                     }
                     if (inputNames[nInputs-1] != NULL)
                          fprintf(fp,"'%s'};\n",inputNames[nInputs-1]);
                     else fprintf(fp,"'X%d'};\n",nInputs);
                  }
                  fwritePlotCLF(fp);
                  fprintf(fp, "ymin = min(Mids);\n");
                  fprintf(fp, "ymax = max(Mids);\n");
                  fprintf(fp, "h2 = 0.05 * (ymax - ymin);\n");
                  if (psPlotTool_ == 1)
                  {
                     fprintf(fp, "Mids = matrix(Mids, %d, %d);\n",nInputs,nInputs);
                     fprintf(fp, "Mids = Mids';\n");
                     fprintf(fp, "drawlater\n");
                     fprintf(fp, "hist3d(Mids);\n");
                     fprintf(fp, "a=gca();\n");
                     fprintf(fp, "a.data_bounds=[0, 0, 0; %d+1, %d+1, ymax];\n",
                             nInputs,nInputs);
                     fprintf(fp, "newtick = a.x_ticks;\n");
                     fprintf(fp, "newtick(2) = [1:%d]';\n",nInputs);
                     fprintf(fp, "newtick(3) = Str';\n");
                     fprintf(fp, "a.x_ticks = newtick;\n");
                     fprintf(fp, "a.x_label.font_size = 3;\n");
                     fprintf(fp, "a.x_label.font_style = 4;\n");
                     fprintf(fp, "a.y_ticks = newtick;\n");
                     fprintf(fp, "a.y_label.font_size = 3;\n");
                     fprintf(fp, "a.y_label.font_style = 4;\n");
                     fprintf(fp, "a.rotation_angles = [5 -70];\n");
                     fprintf(fp, "drawnow\n");
                  }
                  else
                  {
                     fprintf(fp, "Mids = reshape(Mids, %d, %d);\n",nInputs,nInputs);
                     fprintf(fp, "Mids = Mids';\n");
                     fprintf(fp, "bar3(Mids,0.8);\n");
                     fprintf(fp, "axis([0 %d+1 0 %d+1 0 ymax])\n",nInputs,nInputs);
                     fprintf(fp, "set(gca,'XTickLabel',Str);\n");
                     fprintf(fp, "set(gca,'YTickLabel',Str);\n");
                     fprintf(fp, "set(gca, 'fontsize', 12)\n");
                     fprintf(fp, "set(gca, 'fontweight', 'bold')\n");
                  }
                  fwritePlotAxes(fp);
                  fwritePlotTitle(fp,"Sobol Second Order Indices (+ first order)");
                  fwritePlotZLabel(fp, "Sobol Indices");
                  fwritePlotXLabel(fp, "Inputs");
                  fwritePlotYLabel(fp, "Inputs");
                  fclose(fp);
                  if (psPlotTool_ == 1)
                       printf("rssobol2 plot file = scilabrssobol2.sci\n");
                  else printf("rssobol2 plot file = matlabrssobol2.m\n");
               }
               else
               {
                  printf("Total variance = 0. Hence, no interaction effect plot.\n");
               }
               pdata->clean();
            }
            delete anaManager;
         }
      }

      // several qsa methods 
      // ----------------------------------------------------------------
      else if (!strcmp(command, "rs_qsa2") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rs_qsa: RS-based quantitative sensitivity analysis\n");
            printf("syntax: rs_qsa (no argument needed)\n");
            printf("Note: to facilitate processing, all expert modes have\n");
            printf("      been suppressed.\n");
            printf("Note: This command differs from rssobol1, rssoboltsi,\n");
            printf("      and the command 'me' in that it uses bootstrapped\n");
            printf("      samples multiple times to get the errors in Sobol'\n");
            printf("      indices due to response surface errors.\n");
            continue;
         }
         printf("ERROR: no data (load data first).\n");
      }
      else if (!strcmp(command, "rs_qsa") && nInputs > 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rs_qsa: RS-based quantitative sensitivity analysis\n");
            printf("syntax: rs_qsa (no argument needed)\n");
            printf("Note: to facilitate processing, all expert modes have\n");
            printf("      been suppressed.\n");
            printf("Note: This command differs from rssobol1, rssoboltsi,\n");
            printf("      and the command 'me' in that it uses bootstrapped\n");
            printf("      samples multiple times to get the errors in Sobol'\n");
            printf("      indices due to response surface errors.\n");
         }
         if (psuadeIO_ == NULL) printf("ERROR: data file not loaded.\n");
         else
         {
            outputID = 0;
            if (nOutputs > 1)
            {
               sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
               outputID = getInt(1, nOutputs, pString);
               outputID--;
            }

            printf("Which main or total effect analyzer ? \n");
            printf("1. Sobol' main effect using McKay's method (replicated LH).\n");
            printf("2. Sobol' main effect method using on numerical integration. \n");
            printf("3. Sobol' total sensitivity using numerical integration.\n");
            sprintf(pString, "Which  method (1, 2, or 3) ? ");
            method = getInt(1, 3, pString);
            if      (method == 1) analysisMethod = PSUADE_ANA_ME;
            else if (method == 2) analysisMethod = PSUADE_ANA_RSSOBOL1;
            else if (method == 3) analysisMethod = PSUADE_ANA_RSSOBOLTSI;
            anaManager = new AnalysisManager();
            anaManager->setup(analysisMethod, 0);
 
            if (method == 1)
            {
               faType = -1;
               faPtr = genFA(faType, nInputs, iOne, nSamples);
               faPtr->setBounds(iLowerB, iUpperB);
               faPtr->setOutputLevel(0);
            }

            psuadeIO_->getParameter("ana_diagnostics",pPtr);
            saveDiag = pPtr.intData_;
            psuadeIO_->updateAnalysisSection(-1,-1,-1,-1,-1,-1);
            saveMode1 = psAnaExpertMode_;
            psAnaExpertMode_ = 0;
            saveMode2 = psAnalysisInteractive_;
            psAnalysisInteractive_ = 0;
            psuadeIO_->getParameter("method_sampling", pPtr);
            saveMethod = pPtr.intData_;

            tempM = new double[nSamples*nInputs];
            tempV = new double[nSamples];
            states = new int[nSamples];
            count2 = 100000;
            nReps  = 200;
            tempX  = new double[count2*nInputs];
            for (ii = 0; ii < count2*nInputs; ii++) tempX[ii] = 0.0;
            tempY  = new double[count2];
            for (ii = 0; ii < count2; ii++) tempY[ii] = 0.0;
            tempI  = new int[count2];
            for (ii = 0; ii < count2; ii++) tempI[ii] = 1;

            sprintf(pString, "How many times to run it (10 - 1000) : ");
            count = getInt(10, 1000, pString);
            tempT = new double[count*nInputs];
            psuadeIO_->getParameter("ana_use_input_pdfs", pPtr);
            usePDFs = pPtr.intData_;
            if (usePDFs == 1)
            {
               pdfman = new PDFManager();
               pdfman->initialize(psuadeIO_);
               vecOut.setLength(count2*nInputs);
               vecUpper.load(nInputs, iUpperB);
               vecLower.load(nInputs, iLowerB);
            }

            for (kk = 0; kk < count; kk++)
            {
               printf("rsme: iteration %d\n", kk+1);
               for (jj = 0; jj < nSamples; jj++)
               {
                  ind = PSUADE_rand() % nSamples;
                  for (ii = 0; ii < nInputs; ii++)
                     tempM[jj*nInputs+ii] = sampleInputs[ind*nInputs+ii];
                  tempV[jj] = sampleOutputs[ind*nOutputs+outputID];
                  states[jj] = sampleStates[jj];
               }
               if (method == 1)
               {
                  jj = -999;
                  status = faPtr->genNDGridData(tempM,tempV,&jj,NULL,NULL);

                  sampPtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_LHS);
                  sampPtr->setPrintLevel(0);
                  sampPtr->setInputBounds(nInputs, iLowerB, iUpperB);
                  sampPtr->setOutputParams(1);
                  sampPtr->setSamplingParams(count2, nReps, 1);
                  sampPtr->initialize(0);
                  sampPtr->getSamples(count2, nInputs, 1, tempX, tempY, tempI);
                  if (usePDFs == 1)
                  {
                     vecIn.load(count2*nInputs, tempX);
                     pdfman->invCDF(count2, vecIn, vecOut, vecLower, vecUpper);
                     for (ii = 0; ii < count2*nInputs; ii++)
                     tempX[ii] = vecOut[ii];
                  }

                  faPtr->evaluatePoint(count2, tempX, tempY);

                  psuadeIO_->updateInputSection(count2,nInputs,NULL,NULL,NULL,
                                                tempX,NULL);
                  psuadeIO_->updateOutputSection(count2,1,tempY,tempI,NULL);
                  psuadeIO_->updateMethodSection(PSUADE_SAMP_LHS,count2,
                                                 nReps,-1,-1);
               }
               else
               {
                  psuadeIO_->updateInputSection(nSamples,nInputs,NULL,NULL,
                                                NULL,tempM,NULL);
                  psuadeIO_->updateOutputSection(nSamples,1,tempV,states,
                                                 outputNames);
                  psuadeIO_->updateMethodSection(PSUADE_SAMP_MC,nSamples,
                                                 -1,-1,-1);
               }
               anaManager->analyze(psuadeIO_, 0, NULL, 0);
               pdata = psuadeIO_->getAuxData();
               if (pdata->nDbles_ != nInputs)
               {
                  printf("ERROR: nInputs do not match (%d, %d).\n",
                         pdata->nDbles_, nInputs);
                  printf("       Consult PSUADE developers.\n");
                  delete [] tempX;
                  delete [] tempY;
                  delete [] tempI;
                  delete [] tempT;
                  delete [] tempM;
                  delete [] tempV;
                  delete [] states;
                  if (method == 1) delete sampPtr;
                  continue;
               }

               if (pdata->dbleData_ > 0)
                  for (ii = 0; ii < nInputs; ii++)
                     tempT[kk*nInputs+ii] =
                          pdata->dbleArray_[ii]/pdata->dbleData_;
               else
                  for (ii = 0; ii < nInputs; ii++)
                     tempT[kk*nInputs+ii] = pdata->dbleArray_[ii];

               pdata->clean();
               if (method == 1) delete sampPtr;
            }
            if (usePDFs == 1) delete pdfman;
            delete [] tempM;
            delete [] tempV;
            delete [] states;
            tempM = new double[nInputs];
            for (ii = 0; ii < nInputs; ii++)
            {
               tempM[ii] = tempT[ii];
               for (jj = 1; jj < count; jj++) tempM[ii] += tempT[jj*nInputs+ii];
               tempM[ii] /= (double) count;
            }
            tempV = new double[nInputs];
            for (ii = 0; ii < nInputs; ii++)
            {
               tempV[ii] = pow(tempT[ii]-tempM[ii], 2.0);
               for (jj = 1; jj < count; jj++)
                  tempV[ii] += pow(tempT[jj*nInputs+ii]-tempM[ii],2.0);
               tempV[ii] /= (double) (count - 1);
               tempV[ii] = sqrt(tempV[ii]);
            }
            printEquals(0);
            printf("Statistics (based on %d replications): \n",
                   count);
            for (ii = 0; ii < nInputs; ii++)
               printf("Input %4d: mean = %16.8e, std = %16.8e\n",ii+1,
                      tempM[ii],tempV[ii]);
            delete anaManager;
            delete [] tempX;
            delete [] tempY;
            delete [] tempI;
            delete [] tempT;
            delete [] tempM;
            delete [] tempV;
            if (faPtr != NULL) delete faPtr;
            faPtr = NULL;
            tempX = tempY = tempM = tempT = tempV = NULL;
            tempI = NULL;

            psuadeIO_->updateInputSection(nSamples,nInputs,NULL,NULL,NULL,
                                          sampleInputs,NULL);
            psuadeIO_->updateOutputSection(nSamples,nOutputs,sampleOutputs,
                                           sampleStates,outputNames);
            psuadeIO_->updateMethodSection(saveMethod,nSamples,-1,-1,-1);
            psuadeIO_->updateAnalysisSection(-1,-1,-1,saveDiag,-1,-1);
            psAnaExpertMode_ = saveMode1;
            psAnalysisInteractive_ = saveMode2;
         }
      }

      // Sobol total sensitivity effect
      // ----------------------------------------------------------------
      else if (!strcmp(command, "rssoboltsi") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rssoboltsi: RS-based Sobol' total sensitivity indices\n");
            printf("syntax: rssoboltsi (no argument needed)\n");
            printf("note: to compute error bars for indices, turn on\n");
            printf("      ana_expert mode and set ntimes to 100.\n");
            continue;
         }
         printf("ERROR: no data (load data first).\n");
      }
      else if (!strcmp(command, "rssoboltsi") && nInputs > 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rssoboltsi: RS-based Sobol' total sensitivity indices\n");
            printf("syntax: rssoboltsi (no argument needed)\n");
            printf("note: to compute error bars for indices, turn on\n");
            printf("      ana_expert mode and set ntimes to 100.\n");
            continue;
         }
         if (psuadeIO_ == NULL) printf("ERROR: data file not loaded.\n");
         else
         {
            outputID = 0;
            if (nOutputs > 1)
            {
               sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
               outputID = getInt(1, nOutputs, pString);
               outputID--;
            }
            analysisMethod = PSUADE_ANA_RSSOBOLTSI;
            anaManager = new AnalysisManager();
            anaManager->setup(analysisMethod, 0);
            psuadeIO_->getParameter("ana_diagnostics",pPtr);
            ii = pPtr.intData_;
            psuadeIO_->updateAnalysisSection(-1,-1,-1,-2,-1,-1);
            anaManager->analyze(psuadeIO_, 0, NULL, outputID);
            psuadeIO_->updateAnalysisSection(-1,-1,-1,ii,-1,-1);
            delete anaManager;
         }
      }

      // Sobol group main effect
      else if (!strcmp(command, "rssobolg") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rssobolg: RS-based group Sobol' decomposition\n");
            printf("syntax: rssobolg (no argument needed)\n");
            continue;
         }
         printf("ERROR: no data (load data first).\n");
      }
      else if (!strcmp(command, "rssobolg") && nInputs > 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rssobolg: RS-based group Sobol' decomposition\n");
            printf("syntax: rssobolg (no argument needed)\n");
            continue;
         }
         if (psuadeIO_ == NULL) printf("ERROR: data file not loaded.\n");
         else
         {
            outputID = 0;
            if (nOutputs > 1)
            {
               sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
               outputID = getInt(1, nOutputs, pString);
               outputID--;
            }
            analysisMethod = PSUADE_ANA_RSSOBOLG;
            anaManager = new AnalysisManager();
            anaManager->setup(analysisMethod, 0);
            psuadeIO_->getParameter("ana_diagnostics",pPtr);
            ii = pPtr.intData_;
            psuadeIO_->updateAnalysisSection(-1,-1,-1,outputLevel_,-1,-1);
            anaManager->analyze(psuadeIO_, 0, NULL, outputID);
            psuadeIO_->updateAnalysisSection(-1,-1,-1,ii,-1,-1);
            fgets(lineIn,500,stdin); 
            delete anaManager;
         }
      }

      // create response surface 
      else if (!strcmp(command, "rscreate") && nInputs <= 0)
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
         printf("ERROR: no data (load data first).\n");
      }
      else if (!strcmp(command, "rscreate") && nInputs >= 1)
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
         faFlag = 3;
         sprintf(pString, "Which output to use (1 - %d, 0 - all) ? ",nOutputs);
         outputID = getInt(0, nOutputs, pString);
         psuadeIO_->getParameter("ana_outputid", pPtr);
         kk = pPtr.intData_;
         faPtrsRsEval = new FuncApprox*[nOutputs];
         if (outputID != 0)
         {
            psuadeIO_->updateAnalysisSection(-1, -1, -1, outputLevel_, outputID-1, -1);
            faPtrsRsEval[0] = genFAInteractive(psuadeIO_, faFlag);
            for (ii = 1; ii < nOutputs; ii++) faPtrsRsEval[ii] = NULL;
         }
         else
         {
            for (ii = 0; ii < nOutputs; ii++)
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

      // evaluate the response surface at a given point
      // ----------------------------------------------------------------
      else if (!strcmp(command, "rseval") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rseval: evaluate a data point after rscreate\n");
            printf("syntax: rseval (no argument needed)\n");
            continue;
         }
         printf("ERROR: response surface has not been created yet.\n");
         printf("       Please load data and use rscreate first.\n");
      }
      else if (!strcmp(command, "rseval") && nInputs >= 1)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rseval: evaluate a data point after rscreate\n");
            printf("syntax: rseval (no argument needed)\n");
            continue;
         }
         if (faPtrsRsEval == NULL) 
            printf("ERROR: no response surface to evaluate.\n");
         else 
         {
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
               fgets(lineIn,500,stdin); 
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
                     if (kk != nInputs)
                     {
                        printf("ERROR: input size does not match nInputs.\n");
                        printf(":      input size in local memory = %d.\n", 
                               nInputs);
                        printf(":      input size from file       = %d.\n", kk);
                        count = 0;
                        fclose(fp);
                        continue;
                     }
                     printf("Number of points to be evaluated = %d\n", count);
                     inputSettings = new double[count*nInputs];
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
                           break;
                        }
                        for (ii = 0; ii < nInputs; ii++)
                           fscanf(fp, "%lg", &(inputSettings[jj*nInputs+ii]));
                     }
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
            else
            {
               inputSettings = new double[nInputs];
               if (dataReg == NULL)
               {
                  printf("ERROR: data register has not been created.\n");
                  continue;
               }
               count = 1;
               for (ii = 0; ii < nInputs; ii++) 
                  inputSettings[ii] = dataReg[ii];
            }
            if (count > 0)
            {
               printf("Evaluated data written to a file ? (y or n) ");
               fgets(winput,10,stdin); 
               if (winput[0] == 'y')
               {
                  fp = fopen("eval_sample","w");
                  fprintf(fp, "%% Inputs Out1 Std1 Out2 Std2 ...\n");
                  if (nOutputs == 1 || (nOutputs > 1 && faPtrsRsEval[1] == NULL))
                     fprintf(fp, "%d %d 2\n", count, nInputs);
                  else
                     fprintf(fp, "%d %d %d\n", count, nInputs, 2*nOutputs);
               }
               else fp = NULL;
               tempY = new double[nOutputs*count];
               tempW = new double[nOutputs*count];
               for (ii = 0 ; ii < nOutputs; ii++)
               {
                  if (faPtrsRsEval[ii] != NULL)
                     dtemp = faPtrsRsEval[ii]->evaluatePointFuzzy(count,
                                    inputSettings,&(tempY[ii*count]),&(tempW[ii*count]));
               }
               for (kk = 0 ; kk < count; kk++)
               {
                  if (fp != NULL)
                     for (ii = 0 ; ii < nInputs; ii++)
                        fprintf(fp, "%e ", inputSettings[kk*nInputs+ii]);
                  printf("Interpolated Point %d: ", kk+1);
                  for (ii = 0 ; ii < nOutputs; ii++)
                  {
                     if (faPtrsRsEval[ii] != NULL)
                     {
                        printf("output %d = %e (stdev = %e) ",
                               ii+1,tempY[ii*count+kk], tempW[ii*count+kk]);
                        if (fp != NULL)
                           fprintf(fp,"%e %e ",tempY[ii*count+kk],tempW[ii*count+kk]);
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
      }

      // evaluate the response surface at given points repeatedly
      // ----------------------------------------------------------------
      else if (!strcmp(command, "rseval_m") && nInputs <= 0)
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
         printf("ERROR: response surface has not been created yet.\n");
         printf("       Please load data and use rscreate first.\n");
      }
      else if (!strcmp(command, "rseval_m") && nInputs >= 1)
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
         if (faPtrsRsEval == NULL) 
            printf("ERROR: no response surface to evaluate.\n");
         else 
         {
            sprintf(dataFile, "psComplete");
            if ((fp=fopen(dataFile,"r")) != NULL)
            {
               printf("Please remove psComplete file and re-do.\n");
               continue;
            }
            inputSettings = new double[nInputs];
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
                  if (count2 <= 0 || kk != nInputs)
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
                     for (ii = 0; ii < nInputs; ii++)
                        fscanf(fp, "%lg", &inputSettings[ii]);
                     for (ii = 0; ii < nOutputs; ii++)
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
      }

      // probe response surface and create scatter plots
      // ----------------------------------------------------------------
      else if (!strcmp(command, "rs_splot") && nInputs <= 0)
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
         printf("ERROR: no data (load data first).\n");
      }
      else if (!strcmp(command, "rs_splot") && nInputs > 0)
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
         if (psuadeIO_ == NULL) printf("ERROR: data not loaded yet.\n");
         else
         {
            outputID = 0;
            if (nOutputs > 1)
            {
               sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
               outputID = getInt(1, nOutputs, pString);
               outputID--;
            }
            Ymax = - 1.0e35;
            Ymin =   1.0e35;
            for (sInd = 0; sInd < nSamples; sInd++)
            {
               if (sampleOutputs[sInd*nOutputs+outputID] > Ymax)
                  Ymax = sampleOutputs[sInd*nOutputs+outputID];
               if (sampleOutputs[sInd*nOutputs+outputID] < Ymin)
                  Ymin = sampleOutputs[sInd*nOutputs+outputID];
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
            sampPtr->setPrintLevel(0);
            sampPtr->setInputBounds(nInputs, iLowerB, iUpperB);
            sampPtr->setOutputParams(1);
            sampPtr->setSamplingParams(kk, 1, 0);
            sampPtr->initialize(0);
            tempX = new double[kk*nInputs];
            tempY = new double[kk];
            states = new int[kk];
            for (ii = 0; ii < kk; ii++) states[ii] = 1;
            sampPtr->getSamples(kk, nInputs, 1, tempX, tempY, states);
            count = 0;
            for (ii = 0; ii < kk; ii++)
            {
               dtemp = faPtr->evaluatePoint(&tempX[ii*nInputs]);
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
               for (iInd = 0; iInd < nInputs; iInd++)
               {
                  fprintf(fp, "X%d = [\n", iInd+1);
                  for (sInd = 0; sInd < kk; sInd++)
                     if (states[sInd] == 1)
                        fprintf(fp, "%e\n",tempX[sInd*nInputs+iInd]);
                  fprintf(fp, "];\n");
               }
               for (iInd = 0; iInd < nInputs; iInd++)
               {
                  fwritePlotCLF(fp);
                  fprintf(fp, "plot(X%d,Y,'X','MarkerSize',13)\n",iInd+1);
                  sprintf(pString, "%s vs %s",outputNames[outputID],inputNames[iInd]);
                  fwritePlotTitle(fp, pString);
                  fwritePlotXLabel(fp, inputNames[iInd]);
                  fwritePlotYLabel(fp, outputNames[outputID]);
                  fwritePlotAxes(fp);
                  if (iInd < nInputs-1) 
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
      }
      // set interactive mode
      // ----------------------------------------------------------------
      else if (!strcmp(command, "interactive"))
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
         if (psAnalysisInteractive_ == 0)
         {
            psAnalysisInteractive_ = 1;
            printf("interactive mode on\n");
         }
         else
         {
            psAnalysisInteractive_ = 0;
            printf("interactive mode off\n");
         }
      }

      // set IO expert mode
      // ----------------------------------------------------------------
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

      // set response surface expert mode
      // ----------------------------------------------------------------
      else if (!strcmp(command, "rs_expert"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("This command turns on/off response surface expert mode.\n");
            printf("This mode, when turned on, will trigger prompts for\n");
            printf("additional information regarding creation of response surfaces.\n");
            printf("When this mode is off, default setting will be used.\n");
            continue;
         }
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

      // set analysis expert mode
      // ----------------------------------------------------------------
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

      // set sampling expert mode
      // ----------------------------------------------------------------
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

      // set optimization expert mode
      // ----------------------------------------------------------------
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

      // set grand master mode
      // ----------------------------------------------------------------
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
         if (psGMMode_ == 0)
         {
            psGMMode_ = 1;
            printf("GM mode on\n");
         }
         else
         {
            psGMMode_ = 0;
            printf("GM mode off\n");
         }
      }

      // use configure file instead of interactive 
      // ----------------------------------------------------------------
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

      // set scilab as plot tool
      // ----------------------------------------------------------------
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

      // set max data size 
      // ----------------------------------------------------------------
      else if (!strcmp(command, "rsmax"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rsmax: set maximum sample size for response surface\n"); 
            printf("       (used to guard against large sample sizes that\n");
            printf("        make response surface generation very expensive.)\n");
            printf("       The default max is 4000.\n");
            continue;
         }
         sscanf(lineIn, "%s %d", command, &psFAMaxDataPts_);
         if (psFAMaxDataPts_ <= 1 || psFAMaxDataPts_ > 100000)
            psFAMaxDataPts_ = 100000;
         printf("psuade max data size for response surface set to %d\n", 
                psFAMaxDataPts_);
      }

      // RS-based brute force calibration 
      // ----------------------------------------------------------------
      else if (!strcmp(command, "rsca") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rsca: RS-based brute-force calibration (use mesh data)\n");
            printf("When the input dimension is low, an alternative to MCMC\n");
            printf("is to do a brute force search on an interpolated mesh.\n");
            continue;
         }
         printf("ERROR: no data (load data first).\n");
      }
      else if (!strcmp(command, "rsca") && nInputs > 4)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rsca: RS-based brute-force calibration (use mesh data)\n");
            printf("When the input dimension is low, an alternative to MCMC\n");
            printf("is to do a brute force search on an interpolated mesh.\n");
            continue;
         }
      }
      else if (!strcmp(command, "rsca") && nInputs >= 1)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("rsca: RS-based brute-force calibration (use mesh data)\n");
            printf("When the input dimension is low, an alternative to MCMC\n");
            printf("is to do a brute force search on an interpolated mesh.\n");
            continue;
         }
         if (psPlotTool_ != 0)
         {
            printf("INFO: rsca is currently not available for scilab.\n");
            continue;
         }
         if (faPtr != NULL) delete faPtr;

         inputSettings = new double[nInputs];
         iplot1 = 0;
         iplot2 = 1;
         iplot3 = 2;
         kk = nInputs;
         if (nInputs > 3)
         {
            printf("Currently, you can only select at most 3 inputs.\n");
            sprintf(pString, "Enter the number of inputs to study (1 - 3) : ");
            kk = getInt(1, 3, pString);
         }
         if (kk == 1) iplot2 = iplot3 = -2;
         if (kk == 2) iplot3 = -2;
         if (nInputs > kk)
         {
            sprintf(pString, "Enter the first input (1 - %d) : ", nInputs);
            iplot1 = getInt(1, nInputs, pString);
            iplot1--;
            if (kk > 1)
            {
               iplot2 = iplot1;
               while (iplot1 == iplot2)
               {
                  sprintf(pString, "Enter the second input (1 - %d), not %d : ",
                          nInputs, iplot1+1);
                  iplot2 = getInt(1, nInputs, pString);
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
                          nInputs, iplot1+1, iplot2+1);
                  iplot3 = getInt(1, nInputs, pString);
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
         faPtr->setBounds(iLowerB, iUpperB);
         if (nInputs > kk)
            printf("The other inputs will be set to their nominal values.\n");
         for (iInd1 = 0; iInd1 < nInputs; iInd1++)
         {
            if (iInd1 != iplot1 && iInd1 != iplot2 && iInd1 != iplot3)
               inputSettings[iInd1] = 0.5*(iLowerB[iInd1]+iUpperB[iInd1]);
            else inputSettings[iInd1] = 1.0;
         }
         jplot = 0;
         faYIn = new double[nSamples];
         for (sInd = 0; sInd < nSamples; sInd++)
            faYIn[sInd] = sampleOutputs[sInd*nOutputs+jplot];
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
               inputSettings[iplot3] = (iUpperB[iplot3] - iLowerB[iplot3]) *
                                   ii / (nPtsPerDim - 1.0) + iLowerB[iplot3];
               faPtr->gen2DGridData(sampleInputs,faYIn, iplot1, iplot2,
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
            fprintf(fp, "clf\n");
            for (ii = 0; ii < nPtsPerDim; ii++)
            {
               dtemp = (iUpperB[iplot3] - iLowerB[iplot3]) *
                       ii / (nPtsPerDim - 1.0) + iLowerB[iplot3];
               fprintf(fp, "disp(\'Plotting frame %d of %d\')\n",
                       ii+1,nPtsPerDim);
               fprintf(fp, "if twoPlots == 1\n");
               fprintf(fp, "subplot(1,2,1), mesh(X,Y,A%d)\n", ii+1);
               fprintf(fp, "axis([%e %e %e %e %e %e])\n",iLowerB[iplot1],
                       iUpperB[iplot1],iLowerB[iplot2],iUpperB[iplot2],
                       GYmin, GYmax);
               fwritePlotXLabel(fp, inputNames[iplot1]);
               fwritePlotYLabel(fp, inputNames[iplot2]);
               fwritePlotZLabel(fp, outputNames[jplot]);
               fwritePlotAxes(fp);
               fprintf(fp, "colorbar\n");
               fprintf(fp, "title(\'%s Mesh plot, val(3) = %14.7e\')\n",
                       outputNames[jplot], dtemp);
               fprintf(fp, "subplot(1,2,2)\n");
               fprintf(fp, "end\n");
               fprintf(fp, "contourf(X,Y,A%d)\n",ii+1);
               fprintf(fp, "axis([%e %e %e %e])\n",iLowerB[iplot1],
                       iUpperB[iplot1],iLowerB[iplot2],iUpperB[iplot2]);
               fwritePlotXLabel(fp, inputNames[iplot1]);
               fwritePlotYLabel(fp, inputNames[iplot2]);
               fwritePlotAxes(fp);
               fprintf(fp, "colorbar\n");
               fprintf(fp, "colormap(jet)\n");
               fprintf(fp, "title(\'%s contour plot, val(3) = %14.7e\')\n",
                       outputNames[jplot], dtemp);
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
            fprintf(fp, "clf\n");
            faPtr->gen2DGridData(sampleInputs,faYIn,iplot1,iplot2,
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
            fprintf(fp, "clf\n");
            fprintf(fp, "if twoPlots == 1\n");
            fprintf(fp, "subplot(1,2,1), mesh(X,Y,A)\n");
            fwritePlotXLabel(fp, inputNames[iplot1]);
            fwritePlotYLabel(fp, inputNames[iplot2]);
            fwritePlotZLabel(fp, outputNames[jplot]);
            fwritePlotAxes(fp);
            fprintf(fp, "title(\'Mesh Plot for %s\')\n",outputNames[jplot]);
            fprintf(fp, "colorbar\n");
            fprintf(fp, "subplot(1,2,2)\n");
            fprintf(fp, "end\n");
            fprintf(fp, "contourf(X,Y,A)\n");
            fwritePlotXLabel(fp, inputNames[iplot1]);
            fwritePlotYLabel(fp, inputNames[iplot2]);
            fwritePlotAxes(fp);
            fprintf(fp, "colorbar\n");
            fprintf(fp, "colormap(jet)\n");
            fprintf(fp, "title(\'Contour Plot for %s\')\n",outputNames[jplot]);
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
            fprintf(fp, "clf\n");
            faPtr->gen1DGridData(sampleInputs,faYIn, iplot1,
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
            fwritePlotXLabel(fp, inputNames[iplot1]);
            fwritePlotYLabel(fp, outputNames[jplot]);
            fwritePlotAxes(fp);
            fprintf(fp, "title(\'Likelihood Plot for %s\')\n",
                    outputNames[jplot]);
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

      // 1-sample and 2-sample tests
      // ----------------------------------------------------------------
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
         printAsterisks(0);
         printAsterisks(0);
         printEquals(0);
         printf("##### One-input test (Normal(1,4)): \n");
         printf("Sample mean     should be = 1\n");
         printf("Sample std dev  should be = 2\n");
         printf("Sample skewness should be = 0\n");
         printf("Sample Kurtosis should be = 3\n");
         printDashes(0);
         smeans[0] = 1.0;
         sstdevs[0] = 2.0;
         slbounds[0] = -6.0;
         subounds[0] = 8.0;
         pdfman = new PDFManager();
         PDFs[0] = PSUADE_PDF_NORMAL;
         ddata = 1.0;
         corMat.setDim(1,1);
         corMat.setEntry(0,0,ddata);
         pdfman->initialize(1, PDFs, smeans, sstdevs, corMat);
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
         printAsterisks(0);
         printAsterisks(0);
         printEquals(0);
         printf("##### One-input test (LogNormal(1,2)): \n");
         printf("Sample mean     should be = 4.5\n");
         printf("Sample std dev  should be = 5.8\n");
         printDashes(0);
         smeans[0] = 1.0;
         sstdevs[0] = 1.0;
         slbounds[0] = 0.0;
         subounds[0] = 200.0;
         PDFs[0] = PSUADE_PDF_LOGNORMAL;
         ddata = 1.0;
         corMat.setDim(1,1);
         corMat.setEntry(0,0,ddata);
         pdfman->initialize(1, PDFs, smeans, sstdevs, corMat);
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
         printAsterisks(0);
         printAsterisks(0);
         printEquals(0);
         printf("##### One-input test (Weibull(1,1)): \n");
         printf("Sample mean     should be = 1\n");
         printf("Sample std dev  should be = 1\n");
         printf("Sample skewness should be = 2\n");
         printDashes(0);
         smeans[0] = 1.0;
         sstdevs[0] = 1.0;
         slbounds[0] = 0.0;
         subounds[0] = 10.0;
         PDFs[0] = PSUADE_PDF_WEIBULL;
         ddata = 1.0;
         corMat.setDim(1,1);
         corMat.setEntry(0,0,ddata);
         pdfman->initialize(1, PDFs, smeans, sstdevs, corMat);
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

         printAsterisks(0);
         printAsterisks(0);
         printEquals(0);
         printf("##### Two-input test (Normal(0,1)): \n");
         printf("Sample mean     should be = 0\n");
         printf("Sample std dev  should be = 1.414\n");
         printf("Sample skewness should be = 0\n");
         printf("Sample Kurtosis should be = 3\n");
         printDashes(0);
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
         pdfman->initialize(2, PDFs, smeans, sstdevs, corMat);
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
         printAsterisks(0);
         printAsterisks(0);
         printEquals(0);
         printf("##### Two-input test (LogNormal(1,1)): \n");
         printf("Sample mean     should be = 8.9634\n");
         printf("Sample std dev  should be = 8.3081\n");
         printDashes(0);
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
         pdfman->initialize(2, PDFs, smeans, sstdevs, corMat);
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

         printAsterisks(0);
         printAsterisks(0);
         printEquals(0);
         printf("##### Two-input test (Normal(0,1)) with cor = 0.5: \n");
         printf("Sample mean     should be = 0\n");
         printf("Sample std dev  should be = 1.73\n");
         printf("Sample skewness should be = 0\n");
         printf("Sample Kurtosis should be = 3\n");
         printDashes(0);
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
         pdfman->initialize(2, PDFs, smeans, sstdevs, corMat);
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
         printAsterisks(0);
         printAsterisks(0);
         printEquals(0);
         printf("##### Two-input test (LogNormal(1,1)) with cor = 0.5: \n");
         printf("Sample mean     should be = 7.4\n");
         printf("Sample mean     should be = 1.0\n");
         printDashes(0);
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
         pdfman->initialize(2, PDFs, smeans, sstdevs, corMat);
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

         printEquals(0);
         printEquals(0);
         printf("4-input test (Normal(0,1), cor = 0.5, LogNormal(1,1), no cor: \n");
         printf("Sample mean     should be = 8.96\n");
         printf("Sample std dev  should be = 8.54\n");
         printDashes(0);
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
         pdfman->initialize(4, PDFs, smeans, sstdevs, corMat);
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

      // generate discrete value sample
      // ----------------------------------------------------------------
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
         names2 = new char[100];
         strcpy(names2, "Y");
         ioPtr = new PsuadeData();
         ioPtr->updateInputSection(count, kk, NULL, tempT,
                                   tempW, tempX, names); 
         ioPtr->updateOutputSection(count, iOne, tempY, states, &names2); 
         ioPtr->updateMethodSection(PSUADE_SAMP_MC, count, 1, -1, -1);
         ioPtr->writePsuadeFile("psuade_discrete_sample",0);
         printf("The test sample is in file psuade_discrete_sample.\n"); 
         delete [] tempX;
         delete [] tempY;
         delete [] tempT;
         delete [] tempW;
         delete [] states;
         delete ioPtr;
         ioPtr = NULL;
         tempX = tempY = tempT = tempW = NULL;
      }

      // response-surface based multi-objective optimization 
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
         sprintf(pString,"Enter a PSUADE input file (specify inputs/outputs): ");
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
         if (psuadeIO_ != NULL) delete psuadeIO_;
         if (sampler_ != NULL) SamplingDestroy(sampler_);
         if (tagArray != NULL) delete [] tagArray;
         tagArray = NULL;
         if (inputNames != NULL)
         {
            for (ii = 0; ii < nInputs; ii++)
               if (inputNames[ii] != NULL) delete [] inputNames[ii];
            delete [] inputNames;
         }
         if (outputNames != NULL)
         {
            for (ii = 0; ii < nOutputs; ii++)
               if (outputNames[ii] != NULL) delete [] outputNames[ii];
            delete [] outputNames;
         }
         ioPtr = new PsuadeData();
         ioPtr->readPsuadeFile(dataFile);
         ioPtr->getParameter("input_ninputs", pPtr);
         nInputs = pPtr.intData_;
         ioPtr->getParameter("output_noutputs", pPtr);
         nOutputs = pPtr.intData_;
         ioPtr->getParameter("input_lbounds", pLower);
         iLowerB = pLower.dbleArray_;
         ioPtr->getParameter("input_ubounds", pUpper);
         iUpperB = pUpper.dbleArray_;
         printf("Next, enter the name of the application driver. It can be\n");
         printf("an executable or a PSUADE data file to be used to create\n");
         printf("a response surface. In the latter case, make sure the\n");
         printf("number of inputs and outputs are the same as the PSUADE\n");
         printf("input file (%s).\n", dataFile);
         sprintf(pString, "Enter the application driver file now: ");
         getString(pString, winput);
         kk = strlen(winput);
         winput[kk-1] = '\0';
         printf("application driver (RS data) file = %s\n", winput);
         fp = fopen(winput, "r");
         if (fp == NULL)
         {
            printf("ERROR: file %s not found.\n", winput);
            continue;
         }
         ioPtr->updateApplicationSection(winput, winput,-1);
         ioPtr->updateOptimizationSection(10, 1.0e-6);

         int    optCnt=0;
         double *optData[4]; 
         double *optInitX, *optInitY, *optX, *optY;
         optInitX = new double[nInputs];
         for (ii = 0; ii < nInputs; ii++)
            optInitX[ii] = 0.5 * (iLowerB[ii] + iUpperB[ii]);
         optInitY = new double[nOutputs];
         for (ii = 0; ii < nOutputs; ii++) optInitY[ii] = 0.0;
         optX = new double[nInputs];
         for (ii = 0; ii < nInputs; ii++) optX[ii] = 0.0;
         optY = new double[nOutputs];
         for (ii = 0; ii < nOutputs; ii++) optY[ii] = 0.0;
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
         iLowerB = NULL;
         iUpperB = NULL;
      }

      // UQ for aleatoric and epistemic UQ using inner-outer iteration
      // ----------------------------------------------------------------
      else if (!strcmp(command, "so_ua") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("so_ua: UQ for second order uncertainty analysis\n");
            printf("       (uncertainty in parameter probability distributions)\n");
            printf("syntax: so_ua (no argument needed).\n");
            continue;
         }
         printf("ERROR: no data to analyze (load data first).\n");
      }
      else if (!strcmp(command, "so_ua") && nInputs == 1)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("so_ua: UQ for second order uncertainty analysis\n");
            printf("       (uncertainty in parameter probability distributions)\n");
            printf("syntax: so_ua (no argument needed).\n");
            continue;
         }
         printf("ERROR: nInputs must be 2 or more.\n");
      }
      else if (!strcmp(command, "so_ua") && nInputs > 1)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("so_ua: UQ for second order uncertainty analysis\n");
            printf("       (uncertainty in parameter probability distributions)\n");
            printf("syntax: so_ua (no argument needed).\n");
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
         faPtr->setBounds(iLowerB, iUpperB);
         faPtr->setOutputLevel(outputLevel_);
         psuadeIO_->getParameter("ana_outputid", pPtr);
         outputID = pPtr.intData_;
         printf("INFO: outputID = %d\n", outputID+1);
         printf("INFO: outputID obtained from the loaded sample file.\n");
         printf("      If a different output should be selected, modify\n");
         printf("      the sample file, re-load, and re-run this.\n");

         psuadeIO_->getParameter("input_pdfs", pPDFs);
         psuadeIO_->getParameter("input_means", pMeans);
         psuadeIO_->getParameter("input_stdevs", pStds);
         iPDFs  = pPDFs.intArray_;
         iMeans = pMeans.dbleArray_;
         iStds  = pStds.dbleArray_;
         if (pPDFs.nInts_ == 0)
         {
            iPDFs = NULL;
            iMeans = NULL;
            iStds = NULL;
         }
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
         
         nParams = 2 * nInputs;
         tempV = new double[nParams];
         tempW = new double[nParams];
         for (ii = 0; ii < nInputs; ii++)
         {
            printf("Enter uncertainties for input %d:\n", ii+1); 
            if (iPDFs == NULL || iPDFs[ii] == 0)
            {
               printf("Lower and upper bounds = %e %e\n",
                      iLowerB[ii],iUpperB[ii]);
               sprintf(pString, "Enter lower bound for input lower bound: ");
               tempV[ii*2] = getDouble(pString); 
               sprintf(pString, "Enter upper bound for input lower bound: ");
               tempW[ii*2] = getDouble(pString); 
               sprintf(pString, "Enter lower bound for input upper bound: ");
               tempV[ii*2+1] = getDouble(pString); 
               sprintf(pString, "Enter upper bound for input upper bound: ");
               tempW[ii*2+1] = getDouble(pString); 
            }
            else if (iPDFs[ii] == PSUADE_PDF_NORMAL)
            {
               printf("Normal distribution mean and std dev = %e %e\n",
                      iMeans[ii],iStds[ii]);
               sprintf(pString, "Enter lower bound for input mean : ");
               tempV[ii*2] = getDouble(pString); 
               sprintf(pString, "Enter upper bound for input mean : ");
               tempW[ii*2] = getDouble(pString); 
               sprintf(pString, "Enter lower bound for input std dev : ");
               tempV[ii*2+1] = getDouble(pString); 
               sprintf(pString, "Enter upper bound for input std dev : ");
               tempW[ii*2+1] = getDouble(pString); 
            }
            else if (iPDFs[ii] == PSUADE_PDF_LOGNORMAL)
            {
               printf("LogNormal distribution mean and std dev = %e %e\n",
                      iMeans[ii],iStds[ii]);
               sprintf(pString, "Enter lower bound for input mean : ");
               tempV[ii*2] = getDouble(pString); 
               sprintf(pString, "Enter upper bound for input mean : ");
               tempW[ii*2] = getDouble(pString); 
               sprintf(pString, "Enter lower bound for input std dev : ");
               tempV[ii*2+1] = getDouble(pString); 
               sprintf(pString, "Enter upper bound for input std dev : ");
               tempW[ii*2+1] = getDouble(pString); 
            }
            else if (iPDFs[ii] == PSUADE_PDF_BETA)
            {
               printf("Beta distribution alpha and beta = %e %e\n",
                      iMeans[ii],iStds[ii]);
               sprintf(pString, "Enter lower bound for input alpha : ");
               tempV[ii*2] = getDouble(pString); 
               sprintf(pString, "Enter upper bound for input alpha : ");
               tempW[ii*2] = getDouble(pString); 
               sprintf(pString, "Enter lower bound for input beta : ");
               tempV[ii*2+1] = getDouble(pString); 
               sprintf(pString, "Enter upper bound for input beta : ");
               tempW[ii*2+1] = getDouble(pString); 
            }
            else if (iPDFs[ii] == PSUADE_PDF_TRIANGLE)
            {
               printf("Triangle distribution mean and half width = %e %e\n",
                      iMeans[ii],iStds[ii]);
               sprintf(pString, "Enter lower bound for input mean : ");
               tempV[ii*2] = getDouble(pString); 
               sprintf(pString, "Enter upper bound for input mean : ");
               tempW[ii*2] = getDouble(pString); 
               sprintf(pString, "Enter lower bound for input half width : ");
               tempV[ii*2+1] = getDouble(pString); 
               sprintf(pString, "Enter upper bound for input half width : ");
               tempW[ii*2+1] = getDouble(pString); 
            }
         }
 
         int    nSams=100, nSams2=1000;
         int    *samStates = new int[nSams];
         double *samIns  = new double[nSams*nParams];
         double *samOuts = new double[nSams];
         sampPtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_LPTAU);
         sampPtr->setPrintLevel(0);
         sampPtr->setInputBounds(nParams, tempV, tempW);
         sampPtr->setSamplingParams(nSams,-1,0);
         ii = 1;
         sampPtr->setOutputParams(ii);
         sampPtr->initialize(0);
         sampPtr->getSamples(nSams, nParams, ii, samIns, samOuts, samStates);
         delete sampPtr;
         sampPtr = NULL;

         int    *sPDFs  = new int[nInputs];
         double *sMeans = new double[nInputs];
         double *sStds  = new double[nInputs];
         double *oneSample = new double[nInputs];
         if (iMeans == NULL)
         {
            for (ii = 0; ii < nInputs; ii++) sPDFs[ii] = 0;
         }
         else
         {
            for (ii = 0; ii < nInputs; ii++) sPDFs[ii] = iPDFs[ii];
         }  
         corMat.setDim(nInputs,nInputs);
         ddata = 1.0;
         for (ii = 0; ii < nInputs; ii++) corMat.setEntry(ii,ii,ddata);
         if (psPlotTool_ == 1) fp = fopen("scilabsoua.sci", "w");
         else                  fp = fopen("matlabsoua.m", "w");
         fwritePlotCLF(fp);
         for (ss = 0; ss < nSams; ss++)
         { 
            for (ii = 0; ii < nInputs; ii++)
            {
               sMeans[ii] = samIns[ss*nParams+2*ii];
               sStds[ii] = samIns[ss*nParams+2*ii+1];
            }
            pdfman = new PDFManager();
            pdfman->initialize(nInputs, sPDFs, sMeans, sStds, corMat);
            vecLower.load(nInputs, iLowerB);
            vecUpper.load(nInputs, iUpperB);
            vecIn.setLength(nSams2*nInputs);
            pdfman->genSample(nSams2, vecIn, vecLower, vecUpper);
            fprintf(fp, "Y = [\n");
            count = 0;
            for (kk = 0; kk < nSams2; kk++)
            {
               flag = 0;
               for (ii = 0; ii < nInputs; ii++)
                  if (vecIn[kk*nInputs+ii] < iLowerB[ii] ||
                      vecIn[kk*nInputs+ii] > iUpperB[ii]) flag++;
               if (flag == 0)
               {
                  for (ii = 0; ii < nInputs; ii++)
                      oneSample[ii] = vecIn[kk*nInputs+ii];
                  dtemp = faPtr->evaluatePoint(oneSample);
                  fprintf(fp, "%e\n", dtemp);
                  count++;
               }
            }
            if (count < 0.5 * nSams2)
            {
               printf("WARNING: < half of the points are within bounds.\n");
               printf("         Input ranges needs to be adjusted.\n");
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
            if (outputNames != NULL) sprintf(winput, "%s", outputNames[outputID]);
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
         delete faPtr;
         faPtr = NULL;
      }

      // UQ for aleatoric and epistemic UQ using inner-outer iteration
      // ----------------------------------------------------------------
      else if (!strcmp(command, "ae_ua") && nInputs <= 0)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("ae_ua: UQ for aleatoric-epistemic uncertainty analysis\n");
            printf("syntax: ae_ua (no argument needed).\n");
            continue;
         }
         printf("ERROR: no data to analyze (load data first).\n");
      }
      else if (!strcmp(command, "ae_ua") && nInputs == 1)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("ae_ua: UQ for aleatoric-epistemic uncertainty analysis\n");
            printf("syntax: ae_ua (no argument needed).\n");
            continue;
         }
         printf("ERROR: nInputs must be 2 or more.\n");
      }
      else if (!strcmp(command, "ae_ua") && nInputs > 1)
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("ae_ua: UQ for aleatoric-epistemic uncertainty analysis\n");
            printf("syntax: ae_ua (no argument needed).\n");
            continue;
         }
         if (faPtr != NULL) delete faPtr;
         psuadeIO_->getParameter("ana_outputid", pPtr);
         iSave = pPtr.intData_;
         sprintf(pString,"Enter output number (1 - %d) : ",nOutputs);
         outputID = getInt(1, nOutputs, pString);
         outputID--;
         psuadeIO_->updateAnalysisSection(-1,-1,-1,-1,outputID,-1);
         nPtsPerDim = 64;
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
         faPtr->setBounds(iLowerB, iUpperB);
         faPtr->setOutputLevel(outputLevel_);
         printf("INFO: outputID = %d\n", outputID+1);
         printf("INFO: outputID obtained from the loaded sample file.\n");
         printf("      If a different output should be selected, modify\n");
         printf("      the sample file, re-load, and re-run this.\n");

         psuadeIO_->getParameter("input_pdfs", pPDFs);
         psuadeIO_->getParameter("input_means", pMeans);
         psuadeIO_->getParameter("input_stdevs", pStds);
         iPDFs  = pPDFs.intArray_;
         iMeans = pMeans.dbleArray_;
         iStds  = pStds.dbleArray_;
         if (iPDFs == NULL)
         {
            printf("ERROR: cannot perform because all variables are epistemic.\n");
            continue;
         }

         printf("First divide your parameters into aleatoric and epistemic sets.\n");
         int nEpistemic=0;
         int *uTypes = new int[nInputs];
         for (ii = 0; ii < nInputs; ii++) uTypes[ii] = 0;
         kk = 1;
         while (kk > 0)
         {
            sprintf(pString, "Select epistemic parameters (1 - %d, 0 if done) : ",
                    nInputs); 
            kk = getInt(0, nInputs, pString);
            if (kk > 0)
            {
               uTypes[kk-1] = 1;
               nEpistemic++;
            }
         } 
         if (nEpistemic == 0 || nEpistemic == nInputs)
         {
            printf("At least 1 and at most %d epistemic parameters are\n",nInputs-1);
            printf("needed for this command.\n");
            continue;
         }
         printf("You have specified %d epistemic parameters.\n", nEpistemic);

         int    nSams=100, nSams2=1000;
         int    *samStates = new int[nSams];
         double *samIns  = new double[nSams*nEpistemic];
         double *samOuts = new double[nSams];
         double *lbs = new double[nEpistemic];
         double *ubs = new double[nEpistemic];
         kk = 0;
         for (ii = 0; ii < nInputs; ii++)
         {
            if (uTypes[ii] != 0)
            {
               lbs[kk] = iLowerB[ii];
               ubs[kk] = iUpperB[ii];
               kk++;
            }
         }
         sampPtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_LPTAU);
         sampPtr->setPrintLevel(0);
         sampPtr->setInputBounds(nEpistemic, lbs, ubs);
         kk = 1;
         sampPtr->setOutputParams(kk);
         sampPtr->setSamplingParams(nSams,-1,-1);
         sampPtr->initialize(0);
         sampPtr->getSamples(nSams, nEpistemic, kk, samIns, samOuts, samStates);
         delete sampPtr;
         sampPtr = NULL;
         delete [] lbs;
         delete [] ubs;

         int    nAleatoric = nInputs - nEpistemic;
         int    *iPDFs2 = new int[nAleatoric];
         double *iMeans2 = new double[nAleatoric];
         double *iStds2 = new double[nAleatoric];
         double *oneSample = new double[nAleatoric];
         lbs = new double[nAleatoric];
         ubs = new double[nAleatoric];
         ddata = 1.0;
         kk = 0;
         corMat.setDim(nAleatoric,nAleatoric);
         for (ii = 0; ii < nInputs; ii++)
         {
            if (uTypes[ii] != 1)
            {
               corMat.setEntry(kk,kk,ddata);
               iPDFs2[kk] = iPDFs[ii]; 
               iMeans2[kk] = iMeans[ii]; 
               iStds2[kk] = iStds[ii]; 
               lbs[kk] = iLowerB[ii];
               ubs[kk] = iUpperB[ii];
               kk++;
            }
         }
         pdfman = new PDFManager();
         pdfman->initialize(nAleatoric, iPDFs2, iMeans2, iStds2, corMat);
         vecLower.load(nAleatoric, lbs);
         vecUpper.load(nAleatoric, ubs);
         vecIn.setLength(nSams2*nAleatoric);
         pdfman->genSample(nSams2, vecIn, vecLower, vecUpper);
         delete pdfman;
         delete [] iMeans2;
         delete [] iStds2;

         if (psPlotTool_ == 1) fp = fopen("scilabaeua.sci", "w");
         else                  fp = fopen("matlabaeua.m", "w");
         fwritePlotCLF(fp);
         for (ss = 0; ss < nSams; ss++)
         {
            fprintf(fp, "Y = [\n");
            count2 = 0;
            for (ii = 0; ii < nInputs; ii++)
            {
               if (uTypes[ii] == 1)
               {
                  oneSample[ii] = samIns[ss*nEpistemic+count2];
                  count2++;
               }
            }
            count = 0;
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
                  for (ii = 0; ii < nInputs; ii++)
                  {
                     if (uTypes[ii] == 0)
                     {
                        oneSample[ii] = vecIn[kk*nAleatoric+count2];
                        count2++;
                     }
                  }
                  dtemp = faPtr->evaluatePoint(oneSample);
                  fprintf(fp, "%e\n", dtemp);
                  count++;
               }
            }
            if (count < 0.5 * nSams2)
            {
               printf("WARNING: less than half of the points are within bounds.\n");
               printf("         Input ranges may need to be adjusted.\n");
            }
            if (count == 0)
            {
               printf("ERROR: none of the sample points are within bounds.\n");
               continue;
            }
            fprintf(fp, "];\n");
            if (psPlotTool_ == 1) fprintf(fp, "Y = gsort(Y,'g','i');\n");
            else                  fprintf(fp, "Y = sort(Y);\n");
            fprintf(fp, "X = 1:%d;\n",count);
            fprintf(fp, "X = X / %d;\n", count);
            fprintf(fp, "plot(Y,X)\n");
            if (ss == 0)
            {
               sprintf(winput, "Cumulative Distributions");
               fwritePlotTitle(fp, winput);
               fwritePlotAxes(fp);
               if (outputNames != NULL) sprintf(winput, "%s", outputNames[outputID]);
               else                     sprintf(winput, "Output Values");
               fwritePlotXLabel(fp, winput);
               sprintf(winput, "Probabilities");
               fwritePlotYLabel(fp, winput);
               if (psPlotTool_ == 1)
                    fprintf(fp, "set(gca(),\"auto_clear\",\"off\")\n");
               else fprintf(fp, "hold on\n");
            }
         }
         fclose(fp);
         printf("Plot file for aleatoric-epistemic uncertainty analysis is now");
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
         delete faPtr;
         faPtr = NULL;
         oneSample = NULL;
         samIns = NULL;
         samOuts = NULL;
         samStates = NULL;
      }

      // system command
      // ----------------------------------------------------------------
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

      // create random input vector
      // ----------------------------------------------------------------
      else if ((!strcmp(command,"ivec_create") || !strcmp(command,"vcreate")) 
               && (nInputs <= 0))
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
         printf("ERROR: no data (load data first).\n");
      }
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
         if (dataReg != NULL) delete [] dataReg;
         dataReg = new double[nInputs];
         for (ii = 0; ii < nInputs; ii++) 
            dataReg[ii] = iLowerB[ii] + PSUADE_drand() *
                          (iUpperB[ii] - iLowerB[ii]);; 
         printf("Internal sample vector created and the values have been set");
         printf(" to be the mid points.\n");
      }

      // modify an entry in the input vector
      // ----------------------------------------------------------------
      else if ((!strcmp(command,"ivec_modify") || !strcmp(command,"vmodify"))
               && dataReg == NULL)
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
         printf("ERROR: use ivec_create before this command.\n");
      }
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
         sprintf(pString, "Which entry to modify ? (1 - %d) ", nInputs);
         ind = getInt(1, nInputs, pString);
         ind--;
         sprintf(pString, "What is the new value ? (%e - %e) ",
                 iLowerB[ind], iUpperB[ind]);
         ddata = getDouble(pString);
         if (ddata < iLowerB[ind] || ddata > iUpperB[ind])
            printf("WARNING: data out of range (extrapolation).\n");
         dataReg[ind] = ddata;
      }

      // show the input vector
      // ----------------------------------------------------------------
      else if ((!strcmp(command,"ivec_show") || !strcmp(command, "vshow"))
               && dataReg == NULL)
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
         printf("ERROR: use ivec_create before this command.\n");
      }
      else if (!strcmp(command, "ivec_show") || !strcmp(command, "vshow"))
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
         for (ii = 0; ii < nInputs; ii++)
            printf("Input %3d = %e\n", ii+1, dataReg[ii]);
      }

      // show auxillary file formats
      // ----------------------------------------------------------------
      else if (!strcmp(command, "show_format"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("show_format: display different PSUADE file formats \n");
            printf("syntax: show_format (no argument needed)\n");
            printf("User-provided information are often needed in design\n");
            printf("and analysis. These information are provided by users\n");
            printf("to PSUADE at various stages. This command lists many\n");
            printf("such file formats.\n");
            continue;
         }
         printf("RSConstraints file format: \n");
         printf("  line 1: nInputs\n ");
         printf("  line 2: <input (or 0)> <value (nominal val if 0)> \n ");
         printf("  line 3: <input (or 0)> <value (nominal val if 0)> \n ");
         printf("  ... \n");
         printf("RS index file format: \n");
         printf("  line 1: nInputs in rs data (driver) file\n");
         printf("  line 2: 1 <num> <default if num == 0>\n");
         printf("  line 3: 2 <num> <0 if num != 0>\n");
         printf("  line 4: 3 <num> <default if num == 0>\n");
         printf("  line 5: 4 <num> <0 if num != 0>\n");
         printf("  ...\n");
         printf("MOATConstraints file format: \n");
         printf("  line 1: nInputs \n");
         printf("  line 2: <input (or 0)> <value (nominal val if 0)> \n");
         printf("  line 3: <input (or 0)> <value (nominal val if 0)> \n");
         printf("  ... \n");
      }

      // set the driver field in the data file
      // ----------------------------------------------------------------
      else if (!strcmp(command, "setdriver"))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("setdriver: set the driver field in the PSUADE file\n");
            printf("syntax: setdrive <driver executable name>\n");
            continue;
         }
         if (psuadeIO_ == NULL)
         {
            printf("ERROR: data not loaded\n");
         }
         else
         {
            strcpy(dataFile, "\0");
            sscanf(lineIn, "%s %s", command, dataFile);
            if ((fp=fopen(dataFile,"r")) == NULL)
            {
               printf("file %s not found.\n", dataFile);
               printf("syntax: setdriver <file>.\n");
               printf("where <file> is an executable driver file.\n");
               continue;
            }
            psuadeIO_->updateApplicationSection(dataFile,NULL,-1);
            printf("Use 'write' to update your PSUADE input file.\n");
         }
      }

      // start matlab
      // ----------------------------------------------------------------
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

      // outputs current sample in memory
      // ----------------------------------------------------------------
      else if ((!strcmp(command, "curr_sample")) || (!strcmp(command, "sinfo")))
      {
         sscanf(lineIn,"%s %s",command,winput);
         if (!strcmp(winput, "-h"))
         {
            printf("curr_sample: (or sinfo) display loaded sample information.\n");
            printf("syntax: curr_sample or sinfo (no argument needed)\n");
            continue;
         }
         printf("Sample in memory: nSamples = %d\n", nSamples);
         printf("                  nInputs  = %d\n", nInputs);
         printf("                  nOutputs = %d\n", nOutputs);
         if (inputNames != NULL)
         {
            printf("Input names: \n");
            for (ii = 0; ii < nInputs; ii++)
               if (inputNames[ii] != NULL)
                  printf("  Input %4d: %s\n",ii+1,inputNames[ii]);
         }
         if (outputNames != NULL)
         {
            printf("Output names: \n");
            for (ii = 0; ii < nOutputs; ii++)
               if (outputNames[ii] != NULL)
                  printf("  Output %4d: %s\n",ii+1,outputNames[ii]);
         }
         count = 0;
         for (ss = 0; ss < nSamples; ss++) if (sampleStates[ss] == 1) count++;
         printf("Number of valid sample points (state = 1)  = %d\n",count);
         count = 0;
         for (ss = 0; ss < nSamples; ss++)
         {
            flag = 1;
            for (ii = 0; ii < nOutputs; ii++)
               if (sampleOutputs[ss*nOutputs+ii] == PSUADE_UNDEFINED) flag = 0;
            if (flag == 1) count++;
         }
         printf("Number of sample points with valid outputs = %d\n",count);
      }

      // quit psuade interactive session
      // ----------------------------------------------------------------
      else if ((!strcmp(command, "quit")) || (!strcmp(command, "q")))
      {
         break;
      }
      else if ((!strcmp(command, "exit")))
      {
         break;
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

   // -------------------------------------------------------------------
   // quit psuade interactive session
   // -------------------------------------------------------------------
   if (sampleInputs  != NULL) delete [] sampleInputs;
   if (sampleOutputs != NULL) delete [] sampleOutputs;
   if (sampleStates  != NULL) delete [] sampleStates;
   if (iLowerB       != NULL) delete [] iLowerB;
   if (iUpperB       != NULL) delete [] iUpperB;
   if (dataReg       != NULL) delete [] dataReg;
   if (inputNames != NULL)
   {
      for (ii = 0; ii < nInputs; ii++)
         delete [] inputNames[ii];
      delete [] inputNames;
   }
   if (outputNames != NULL)
   {
      for (ii = 0; ii < nOutputs; ii++)
         delete [] outputNames[ii];
      delete [] outputNames;
   }
   if (faPtr != NULL) delete faPtr;
   if (tagArray != NULL) delete [] tagArray;
   if (faPtrsRsEval != NULL)
   {
      for (ii = 0; ii < nOutputs; ii++)
         if (faPtrsRsEval[ii] != NULL) delete faPtrsRsEval[ii];
      delete faPtrsRsEval;
   }
   return 0;
}

