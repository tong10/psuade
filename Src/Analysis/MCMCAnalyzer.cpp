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
// Functions for the class MCMCAnalyzer
// ------------------------------------------------------------------------
// AUTHOR : CHARLES TONG
// DATE   : 2007
// Latest revision: May 2013
// Add model inadequacy function
// ************************************************************************
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

#include "MCMCAnalyzer.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "FunctionInterface.h"
#include "pData.h"
#include "PDFBase.h"
#include "PDFNormal.h"
#include "PDFLogNormal.h"
#include "PDFTriangle.h"
#include "PDFBeta.h"
#include "PDFWeibull.h"
#include "PDFGamma.h"
#include "FuncApprox.h"
#include "FuncApprox/Mars.h"
#include "Psuade.h"
#include "Sampling.h"
#include "RSConstraints.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
MCMCAnalyzer::MCMCAnalyzer() : Analyzer()
{
   setName("MCMC");
   mode_ = 0; // RS mode (mode = 1 : simulator-based mode)
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
MCMCAnalyzer::~MCMCAnalyzer()
{
}

// ************************************************************************
// perform MCMC analysis 
// ------------------------------------------------------------------------
double MCMCAnalyzer::analyze(aData &adata)
{
   int    nInputs, nOutputs, nSamples, its, length, status, iOne=1;
   int    ii, ii2, jj, jj2, index, maxPts=257, ngrid, sumBins, ione=1;
   int    nbins, **bins, *scores, iMax, increment, fail=0, NTotal, index2;
   int    maxSamples, burnInSamples, *pdfFlags, printLevel, ****bins2;
   int    nPlots, *plotIndices=NULL, kk, kk2, faType, setCutOff=0, *Ivec;
   int    *rsIndices=NULL, nFail, cnt, freq=1, scheme=0, modelFormFlag=0;
   int    dnSamples=0, dnInputs=0, modelFormScheme=0, rsErrFlag=0;
   double *dSamInputs=NULL, *dSamMeans=NULL, *dSamStdevs=NULL;
   double *XRange, *XGuess, *XGuessD, *XDist, *means, *sigmas, Xtemp, *YY;
   double *X, *lower, *upper, ddata, ddata2, Ymax, *Xmax, *Y;
   double *inputMeans, *inputStdevs, Ytemp, Ytemp2, stdev, stdev2;
   double cutoffHi, cutoffLo, *rsValues=NULL;
   char   lineIn[500], charString[500], cfname[501], *rsFile;
   FILE   *fp, *fp2=NULL;
   string sfname, iline;
   size_t compFlag;
   ifstream   ifile;
   pData      pPtr, qData, pInputs, pOutputs;
   FuncApprox **faPtrs=NULL, **faPtrs1=NULL;
   PDFBase    **inputPDFs;
   FunctionInterface *funcIO=NULL;
   PsuadeData *ioPtr=NULL;
   RSConstraints *constrPtr;

   printAsterisks(0);
   printf("*                     MCMC Optimizer\n");
   printEquals(0);
   printf("TO GAIN ACCESS TO DIFFERENT OPTIONS: TURN ON\n\n");
   printf(" * ana_expert to finetune MCMC parameters, \n");
   printf("   (e.g. sample size for burn-in can be adjusted).\n");
   printf(" * rs_expert to finetune response surface used in MCMC,\n");
   printf(" * printlevel 3 to display more diagnostic information.\n");
   printf(" * printlevel 4 to display more iteration information.\n");
   printf(" * printlevel >=5 reserved only for expert diagnostics only.\n");
   printDashes(0);
   printf("FEATURES AVAILABLE IN THE CURRENT VERSION OF MCMC:\n");
   printf(" * Support different prior distributions (default: uniform)\n");
   printf("   - set ana_expert to use other than uniform distributions.\n");
   printf(" * Support multiple outputs (likelihood from multiple outputs)\n");
   if (psGMMode_ == 1)
   {
      printf(" * Support 2 types of likelihood function\n");
      printf(" * 1. Gaussian likelihood function (default)\n");
      printf(" * 2. top-hat likelihood functions (for nOutputs=1 only)\n");
      printf("   - use ana_expert to select\n");
   }
   printf(" * For some response surfaces such as polynomial and legendre\n");
   printf("   polynomials, bootstrapped MARS, and Gaussian process (GP1),\n");
   printf("   the fitting errors (response surface errors) may be used in\n");
   printf("   the likelihood functions. There are two types of such errors:\n");
   printf("   1. fitting errors (errors incurred in interpolation), and\n");
   printf("   2. model form errors (in this case a discrepancy model is used.)\n");
   printf("      - use ana_expert to select these options\n");
   printf(" * Some input parameters may be disabled (set to default values)\n");
   printf("   - in case these parameters are not to be calibrated\n");
   printf("   - use rs_index_file in PSUADE data file's ANALYSIS section\n");
   printf("     to select them and to set them to default values.\n");
   printf(" * A sample may be generated from the posteriors.\n");
   printf("   - use ana_expert to select this option\n");
   printf(" * This analysis can be terminated gracefully by creating a file\n");
   printf("   named 'psuade_stop' in the same directory during execution.\n");
   printf("   (in case PHASE 2 takes too long.\n");
   printEquals(0);

   printLevel  = adata.printLevel_;
   nInputs     = adata.nInputs_;
   nOutputs    = adata.nOutputs_;
   nSamples    = adata.nSamples_;
   X           = adata.sampleInputs_;
   Y           = adata.sampleOutputs_;
   lower       = adata.iLowerB_;
   upper       = adata.iUpperB_;
   ioPtr       = adata.ioPtr_;
   pdfFlags    = adata.inputPDFs_;
   inputMeans  = adata.inputMeans_;
   inputStdevs = adata.inputStdevs_;
   if (ioPtr != NULL) ioPtr->getParameter("input_names", qData);

   if (psAnaExpertMode_ == 0)
   {
      if (pdfFlags != NULL)
      {
         ii2 = 0;
         for (ii = 0; ii < nInputs; ii++) ii2 += pdfFlags[ii];
         if (ii2 > 0)
         {
            printf("MCMC ERROR: non-uniform prior distribution currently\n");
            printf("            not supported in non-expert mode.\n");
            printf("            Continue with uniform priors.\n");
            for (ii = 0; ii < nInputs; ii++) pdfFlags[ii] = 0;
         }
      }
   }

   // error checking
   if (nInputs <= 0 || nOutputs <= 0)
   {
      printf("MCMC ERROR: invalid nInputs or nOutputs.\n");
      printf("    nInputs  = %d\n", nInputs);
      printf("    nOutputs = %d\n", nOutputs);
      return PSUADE_UNDEFINED;
   }
   if (mode_ == 1 && ioPtr == NULL)
   {
      printf("MCMC ERROR : Direct simulation mode is requested but no\n");
      printf("             information on the simulator is given (missing\n");
      printf("             driver program in your PSUADE data file).\n");
      return PSUADE_UNDEFINED;
   }
   if (nOutputs > 1 && mode_ == 1 && ioPtr != NULL)
   {
      printf("MCMC ERROR : Direct simulation mode does not support nOutputs>1.\n");
      printf("             Only response surface mode supports nOutputs > 1.\n");
      printf("             (nOutputs > 1 can be handled in your program by\n");
      printf("             computing the likelihood function yourself and\n");
      printf("             then use mean of zero and std dev of 1).\n");
      return PSUADE_UNDEFINED;
   }
   status = 0;
   for (ii = 0; ii < nSamples*nOutputs; ii++)
      if (Y[ii] > 0.9*PSUADE_UNDEFINED) status = 1;
   if (status == 1)
   {
      printf("MCMC ERROR: Some outputs are undefined (having value larger\n");
      printf("            than 1e35). Prune the undefined sample points\n");
      printf("            first and re-run.\n");
      return PSUADE_UNDEFINED;
   }

   if (psGMMode_ == 1 && nOutputs == 1)
   {
      printf("*** OPTION TO USE TOP-HAT LIKELIHOOD FUNCTION:\n\n");
      printf("You may use a top-hat (uniform) likelihood function instead\n");
      printf("of the default Gaussian likelihood function. This option may\n");
      printf("sometimes be desirable for quantifying epistemic uncertainties\n");
      printf("with the purpose of identifying the FEASIBLE region (NOTE: using\n");
      printf("MCMC for epistemic uncertainties may not be needed. Alternate\n");
      printf("ways include using 'ofilter' in command line mode). With this\n");
      printf("option the lower and upper cutoffs for the simulator (or response\n");
      printf("surface) output will have to be provided to run the optimization.\n");
      printf("Also, with this option no data is requested for building the\n");
      printf("likelihood function. Thus, this works if you just want to find\n");
      printf("the part of the parameter space where the simulator output is\n");
      printf("bounded.\n");
      printf("Only nOutputs=1 is allowed for this option to simplify things.\n");
      printf("CAUTION: prescribing cutoff carelessly (too narrow) may\n");
      printf("         result in zero posteriors and failure.\n");
      printf("NOTE: IF YOU DON'T KNOW WHAT THIS FEATURE IS, JUST RESPOND 'n'.\n");
      printf("===> Select top-hat likelihood function option ? (y or n) ");
      scanf("%s", charString);
      fgets(lineIn,500,stdin);
      if (charString[0] != 'y') setCutOff = 0;
      else
      {
         sprintf(charString,"Enter lower cutoff : ");
         cutoffLo = getDouble(charString);
         sprintf(charString,"Enter upper cutoff ( > %e) : ", cutoffLo);
         cutoffHi = cutoffLo;
         while (cutoffHi <= cutoffLo)
         {
            cutoffHi = getDouble(charString);
            if (cutoffHi <= cutoffLo)
               printf("Upper cutoff must be > %e.\n", cutoffLo);
         }
         setCutOff = 1;
      }
      printEquals(0);
   }

   // get output means and standard deviations from a file
   if (setCutOff == 0)
   {
      printf("*** NEED DATA TO CREATE GAUSSIAN LIKELIHOOD FUNCTION: \n\n");
      printf("MCMC is going to create a Gaussian likelihood function.\n");
      printf("Please provide a data file containing design parameter\n");
      printf("values, mean, and standard deviation of the observation\n");
      printf("data for each output.\n");
      printf("NOTE: the design parameters are useful only if you wish to\n");
      printf("      include a discrepancy function (if you don't know what\n");
      printf("      this is, simply do not specify any design parameters).\n");
      printf("IMPORTANT: IF m DESIGN PARAMETERS ARE SPECIFIED, THE FIRST\n");
      printf("      m PARAMETERS IN THE SIMULATOR (OR RESPONSE SURFACE)\n");
      printf("      WILL BE TAKEN AS THE m DESIGN (AND NOT CALIBRATION)\n");
      printf("      PARAMETERS.\n");
      printDashes(0);
      printf("*** THE FORMAT OF THIS DATA FILE IS: (O1 means output 1) \n");
      printf("PSUADE_BEGIN\n");
      printf("<no. of experiments p> <no. of design parameters m> <nOutputs n>\n");
      printf("1 <design values...> <O1 mean> <O1 std dev> ... <On std dev> \n");
      printf("2 <design values...> <O1 mean> <O1 std dev> ... <On std dev> \n");
      printf("...\n");
      printf("p <design values...> <O1 mean> <O1 std dev> ... <On std dev> \n");
      printf("PSUADE_END\n");
      printDashes(0);
      printf("The likelihood function to be constructed is in the form of:\n");
      printf("      C exp(-0.5*S) \n");
      printf("where C is some internal normalization constant and\n");
      printf("      S = 1/ p sum_{k=1}^p sum_{i=1)^n (Y_ki - m_ki)^2/sd_ki^2\n");
      printf("where n is the number of outputs and m_ki and\n");
      printf("      sd_ki are the mean and standard deviation of output i of\n");
      printf("      group k.\n");
      printDashes(0);
      printf("NOTE: Alternately, your simulator (or response surface) output\n");
      printf("      may instead be some error measure from comparison of all\n");
      printf("      model outputs with observational data. In this case, set\n");
      printf("      nOutputs=1, mean = 0, and standard deviation = 1 in this\n");
      printf("      specification file (that is, your simulation output is S\n");
      printf("      above and MCMC will generate likelihood C exp(-0.5 S).\n");
      printf("===> Enter the spec file for building the likelihood function : ");
      cin >> sfname;
      fgets(lineIn, 500, stdin); 
      kk = sfname.size();
      if (kk < 500)
      {
         cfname[kk] = '\0';
         sfname.copy(cfname, kk, 0);
         ifile.open(cfname);
         if (! ifile.is_open())
         {
            printf("MCMC ERROR : cannot open spec file %s.\n", cfname);
            return PSUADE_UNDEFINED;
         }
      }
      else
      {
         printf("MCMC ERROR: file name too long.\n");
         return PSUADE_UNDEFINED;
      }
      getline (ifile, iline);
      compFlag = iline.compare(0,12,"PSUADE_BEGIN");
      if (compFlag == 0)
      {
         ifile >> dnSamples;
         if (dnSamples <= 0)
         {
            printf("MCMC ERROR: number of experimental data points <= 0.\n");
            ifile.close();
            return PSUADE_UNDEFINED;
         }
         printf("SPEC FILE: Number of observation set    = %d\n", dnSamples);
         ifile >> dnInputs;
         if (dnInputs < 0)
         {
            printf("MCMC ERROR: number of experimental design variables < 0.\n");
            ifile.close();
            return PSUADE_UNDEFINED;
         }
         if (dnInputs > nInputs)
         {
            printf("MCMC ERROR: number of experimental design variables %d\n",
                   dnInputs);
            printf("            cannot be larger than number of inputs to\n");
            printf("            the simulator %d.\n", nInputs);
            ifile.close();
            return PSUADE_UNDEFINED;
         }
         printf("SPEC FILE: Number of design parameters  = %d\n", dnInputs);
         ifile >> kk;
         if (kk != nOutputs)
         {
            printf("MCMC ERROR: number of outputs for each experiment does not\n");
            printf("            match the simulator nOutputs.\n");
            ifile.close();
            return PSUADE_UNDEFINED;
         }
         printf("SPEC FILE: Number of simulation outputs = %d\n", nOutputs);
         if (dnInputs > 0) dSamInputs = new double[dnSamples*dnInputs];
         dSamMeans = new double[dnSamples*nOutputs];
         dSamStdevs = new double[dnSamples*nOutputs];
         for (ii = 0; ii < dnSamples; ii++)
         {
            ifile >> kk;
            if (kk != ii+1)
            {
               printf("MCMC ERROR: invalid index %d at line %d in spec file.\n",
                      kk, ii+2);
               printf("            (Expecting %d).\n", ii+1);
               printf("==> check line %d\n", ii+3);
               ifile.close();
               if (dSamInputs != NULL) delete [] dSamInputs;
               delete [] dSamMeans;
               delete [] dSamStdevs;
               return PSUADE_UNDEFINED;
            }
            for (jj = 0; jj < dnInputs; jj++) ifile >> dSamInputs[ii*dnInputs+jj];
            for (jj = 0; jj < nOutputs; jj++)
            {
               ifile >> dSamMeans[ii*nOutputs+jj];
               ifile >> dSamStdevs[ii*nOutputs+jj];
               if (printLevel >= 3)
                  printf("Data mean/stdev = %16.8e %16.8e\n",dSamMeans[ii*nOutputs+jj],
                         dSamStdevs[ii*nOutputs+jj]);
               if (dSamStdevs[ii*nOutputs+jj] < 0.0)
               {
                  ifile.close();
                  printf("MCMC ERROR: std dev in spec file should be > 0.\n");
                  printf("==> check the last entry in line %d\n", ii+3);
                  if (dSamInputs != NULL) delete [] dSamInputs;
                  delete [] dSamMeans;
                  delete [] dSamStdevs;
                  return PSUADE_UNDEFINED;
               }
            }
         }
      }
      else
      {
         printf("MCMC ERROR: PSUADE_BEGIN missing at the beginning of the spec file.\n");
         ifile.close();
         return PSUADE_UNDEFINED;
      }
      getline (ifile, iline);
      getline (ifile, iline);
      compFlag = iline.compare(0,10,"PSUADE_END");
      if (compFlag != 0)
      {
         printf("MCMC ERROR: PSUADE_END missing at the end of the spec file.\n");
         ifile.close();
         delete [] dSamMeans;
         delete [] dSamStdevs;
         return PSUADE_UNDEFINED;
      }
      ifile.close();
      printEquals(0);

      if (psAnaExpertMode_ == 1 && mode_ == 0)
      {
         printf("*** OPTION TO INCLUDE RESPONSE SURFACE UNCERTAINTIES: \n\n"); 
         printf("To incorporate the response surface errors into the\n");
         printf("likelihood function, make sure that either Gaussian\n");
         printf("process/Kriging, polynomial regression, or bootstrapped\n");
         printf("MARS response surface is selected in the simulation\n");
         printf("data file. Otherwise, no RS uncertainties will be included.\n\n");
         printf("NOTE: if you don't know what this feature is, just say no.\n");
         printf("===> Include response surface uncertainties? (y or n) ");
         scanf("%s", charString);
         fgets(lineIn,500,stdin);
         if (charString[0] == 'y') rsErrFlag = 1;
         printEquals(0);
      }
   }

   //    There is an option to specify a rs index file in the data file.
   //    This option allows one to disable a certain input in the MCMC
   //    optimization ==> faPtrs.
   if (ioPtr != NULL)
   {
      ioPtr->getParameter("ana_rstype", pPtr);
      faType = pPtr.intData_;
      ioPtr->getParameter("ana_rsindexfile", pPtr);
      rsFile = pPtr.strArray_[0];
      if (strcmp(rsFile, "NONE"))
      {
         printf("A response surface index file has been specified.\n");
         fp = fopen(rsFile, "r");
         if (fp == NULL)
         {
            printf("MCMC ERROR: rs_index_file %s not found.\n",rsFile);
            exit(1);
         }
         else
         {
            printf("INFO: rs_index_file %s found.\n",rsFile);
            fscanf(fp,"%d", &kk);
            if (kk != nInputs)
            {
               printf("MCMC ERROR: invalid nInputs in rs_index_file.\n");
               printf("  Data format should be: \n");
               printf("  line 1: nInputs in rs data (driver) file\n");
               printf("  line 2: 1 <1 or 0> <default value if first number==0>\n");
               printf("  line 3: 2 <2 or 0> <0 if first number != 0>\n");
               printf("  line 4: 3 <3 or 0> <default value if first number==0>\n");
               printf("  line 5: 4 <4 or 0> <0 if first number != 0>\n");
               printf("  ...\n");
               exit(1);
            }
            rsIndices = new int[nInputs];
            rsValues = new double[nInputs];
            for (ii = 0; ii < nInputs; ii++) rsIndices[ii] = 0;
            for (ii = 0; ii < nInputs; ii++)
            {
               fscanf(fp, "%d", &kk);
               if (kk != ii+1)
               {
                  printf("MCMC ERROR: 1st index in indexFile = %d (must be %d]).\n",
                         rsIndices[ii], ii+1);
                  printf("  Data format should be: \n");
                  printf("  line 1: nInputs in rs data (driver) file\n");
                  printf("  line 2: 1 <1 or 0> <default value if first number==0>\n");
                  printf("  line 3: 2 <2 or 0> <0 if first number != 0>\n");
                  printf("  line 4: 3 <3 or 0> <default value if first number==0>\n");
                  printf("  line 5: 4 <4 or 0> <0 if first number != 0>\n");
                  printf("  ...\n");
                  exit(1);
               }
               fscanf(fp, "%d", &rsIndices[ii]);
               if (rsIndices[ii] == 0)
                  printf("MCMC INFO: input %3d inactive\n",ii+1);

               if (rsIndices[ii] < 0 || rsIndices[ii] > nInputs)
               {
                  printf("MCMC INFO: input %3d = %d invalid\n",ii+1,rsIndices[ii]);
                  exit(1);
               }
               rsIndices[ii]--;
               fscanf(fp, "%lg", &rsValues[ii]);
               printf("RS_index %5d = %5d %16.8e\n",ii+1,rsIndices[ii]+1,
                         rsValues[ii]);
            }
            fclose(fp);
            printf("Response surface index information: \n");
            for (ii = 0; ii < nInputs; ii++)
               printf("Input %4d: index = %4d, default = %e\n",
                      ii+1, rsIndices[ii]+1, rsValues[ii]);
         }
      }
   }
   else if (mode_ == 0)
   {
      printf("MCMC INFO: since ioPtr=NULL, assume MARS as reponse surface.\n");
      faType = PSUADE_RS_MARS;
   }
   if (nSamples > 0)
   {
      printf("MCMC INFO: creating response surfaces for all outputs.\n");
      faPtrs = new FuncApprox*[nOutputs];
      YY = new double[nSamples];
      for (ii = 0; ii < nOutputs; ii++)
      {
         faType = -1;
         printf("MCMC INFO: creating response surface for output %d.\n",ii+1);
         faPtrs[ii] = genFA(faType, nInputs, iOne, nSamples);
         faPtrs[ii]->setNPtsPerDim(16);
         faPtrs[ii]->setBounds(lower, upper);
         faPtrs[ii]->setOutputLevel(0);
         length = -999;
         for (kk = 0; kk < nSamples; kk++) YY[kk] = Y[kk*nOutputs+ii];
         status = faPtrs[ii]->genNDGridData(X, YY, &length, NULL, NULL);
         if (status != 0)
         {
            printf("MCMC ERROR: Something wrong in creating response surface.\n");
            printf("            Consult PSUADE developers.\n");
            for (kk = 0; kk <= ii; kk++) delete [] faPtrs[kk];
            delete [] faPtrs;
            delete [] YY;
            return PSUADE_UNDEFINED;
         }
      }
      delete [] YY;
   }
   else if (mode_ == 0)
   {
      printf("MCMC ERROR: no sample given for generating likelihood function.\n");
      return PSUADE_UNDEFINED;
   }

   //    ==> burnInSamples, maxSamples, nbins, plotIndices, nPlots
   //    ==> modelFormFlag, modelFormScheme, fp2 (posterior sample)
   burnInSamples = 20000;
   maxSamples = 20000;
   nbins = 20;
   printEquals(0);
   printf("*** CURRENT SETTINGS OF MCMC PARAMETERS: \n\n");
   printf("MCMC Burn-in sample size      (default) = %d\n", burnInSamples);
   printf("MCMC sample increment         (default) = %d\n", maxSamples);
   printf("MCMC no. of bins in histogram (default) = %d\n", nbins);
   printf("NOTE: sample increment - sample size to run before convergence check\n");
   printf("NOTE: histogram nBins  - define granularity of histogram bar graph\n");
   printf("Turn on ana_expert mode to change these default settings.\n\n");
   if (psAnaExpertMode_ == 1)
   {
      sprintf(charString,"Enter sample size for burn in (5000 - 1000000): ");
      burnInSamples = getInt(5000, 1000000, charString);
      sprintf(charString,"Enter sample increment (5000 - 1000000): ");
      maxSamples = getInt(5000, 1000000, charString);
      sprintf(charString,"Enter the number of histogram bins (10 - 25) : ");
      nbins = getInt(10, 25, charString);
   }
   if (psAnaExpertMode_ == 1)
   {
      printEquals(0);
      printf("*** OPTION TO CREATE POSTERIOR PLOTS FOR SELECTED INPUTS:\n\n"); 
      printf("MCMC will create MATLAB files for the posterior distributions.\n");
      printf("You can choose to generate posterior plots for all inputs, or \n");
      printf("just a selected few (in case there are too many inputs).\n");
      printf("Select inputs for which posterior plots are to be generated.\n");
      sprintf(charString,"Enter input number (%d - %d, -1 for all, 0 to terminate) : ",
              dnInputs, nInputs);
      kk = 1;
      plotIndices = new int[nInputs];
      nPlots = 0;
      while (kk != 0)
      {
         kk = getInt(-1, nInputs, charString);
         if (kk == -1)
         {
            for (ii = dnInputs; ii < nInputs; ii++)
            {
               if (rsIndices == NULL || rsIndices[ii] >= 0) 
                  plotIndices[nPlots++] = ii;
            }
            break;
         }
         if (kk != 0)
         {
            if (kk <= dnInputs) 
               printf("Input %d is not a calibration (but a design) parameter.\n", kk);
            if (rsIndices != NULL && rsIndices[kk-1] < 0) 
               printf("Input %d has been fixed by the rs index file.\n",kk+1);
            else 
               plotIndices[nPlots++] = kk - 1;
         }
      }
      printf("MCMC Plot summary: input number to be plotted are (%d):\n", nPlots);
      for (ii = 0; ii < nPlots; ii++)
         printf("   Input %4d\n", plotIndices[ii]+1);
   }
   else
   {
      nPlots = nInputs - dnInputs;
      plotIndices = new int[nPlots];
      for (ii = dnInputs; ii < nInputs; ii++) plotIndices[ii-dnInputs] = ii;
   }
   if (psAnaExpertMode_ == 1 && setCutOff == 0)
   {
      printEquals(0);
      printf("*** OPTION TO ADD A DISCREPANCY FUNCTION:\n\n"); 
      printf("To use this feature, first make sure that the observation\n");
      printf("data file specified earlier has design parameters specified.\n");
      printf("If you do not know what this is, simply say no.\n");
      printf("REQUIREMENT: THE m DESIGN PARAMETERS IN THE OBSERVATION DATA\n");
      printf("             FILE NEED TO CORRESPOND TO THE FIRST m INPUTS\n");
      printf("             IN THE SIMULATOR DATA FILE.\n");
      printf("NOTE: if you don't know what this feature is, just say no.\n");
      printf("===> Select this option ? (y or n) ");
      scanf("%s", charString);
      fgets(lineIn,500,stdin);
      if (charString[0] == 'y')
      {
         modelFormFlag = 1;
#if 1
         printf("\n");
         printf("*** Select from the following discrepancy function forms:\n");
         printf("1. d(X), where X are design parameters specified in the\n");
         printf("   likelihood data file. In this case, the first %d parameters\n",
                dnInputs);
         printf("   are design variables and the rest are calibration parameters\n");
         printf("   (except the ones fixed by the rs_index list.)\n");
         printf("2. d(X,theta), where theta in addition are the remaining %d \n",
                nInputs-dnInputs);
         printf("   calibration parameters (except the ones fixed by the rs_index\n");
         printf("   list.)\n");
         sprintf(charString,"===> Make your selection: (1 or 2) ");
         modelFormScheme = getInt(1,2,charString);
         modelFormScheme--;
#endif
         if (rsIndices != NULL)
         {
            printf("INFO: the response surface indexing capability is disabled with\n");
            printf("      the use of discrepancy function capability.\n");
            delete [] rsIndices;
            if (rsValues != NULL) delete [] rsValues;
            rsIndices = NULL;
            rsValues = NULL;
         }
      }

      printEquals(0);
      printf("*** OPTION TO CREATE A SAMPLE FROM THE POSTERIOR DISTRIBUTIONS:\n\n");
      printf("In addition to generating the posterior distributions, you can\n");
      printf("also draw a sample from these posteriors. The posterior sample\n");
      printf("can be used as prior sample for another simulator/emulator.\n");
      printf("NOTE: if you don't know what this feature is, just say no.\n");
      printf("===> Select this option ? (y or n) ");
      scanf("%s", charString);
      fgets(lineIn,500,stdin);
      if (charString[0] == 'y')
      {
         fp2 = fopen("MCMCPostSample", "w");
         if (fp2 == NULL)
         {
            printf("ERROR: cannot write to file 'MCMCPostSample'.\n");
            exit(1);
         }
         fprintf(fp2, "PSUADE_BEGIN\n");
         fprintf(fp2, "<nSamples: to be filled in>  %d\n",nInputs);
      }
      printAsterisks(0);
   } 
   else
   {
      fp2 = NULL;
   }
   printEquals(0);

   //    setup input PDF, if there is any
   inputPDFs = new PDFBase*[nInputs];
   for (ii = 0; ii < nInputs; ii++)
   {
      if (pdfFlags != NULL && pdfFlags[ii] == PSUADE_PDF_NORMAL)
      {
         inputPDFs[ii] = (PDFBase *) new PDFNormal(inputMeans[ii], inputStdevs[ii]);
         if (printLevel > 2) 
            printf("Parameter %3d has normal prior distribution (%e,%e).\n",
                   ii+1, inputMeans[ii], inputStdevs[ii]);
      }
      else if (pdfFlags != NULL && pdfFlags[ii] == PSUADE_PDF_LOGNORMAL)
      {
         inputPDFs[ii] = (PDFBase *) new PDFLogNormal(inputMeans[ii],inputStdevs[ii]);
         if (printLevel > 2)
            printf("Parameter %3d has lognormal prior distribution.\n",ii+1);
      }
      else if (pdfFlags != NULL && pdfFlags[ii] == PSUADE_PDF_TRIANGLE)
      {
         inputPDFs[ii] = (PDFBase *) new PDFTriangle(inputMeans[ii],inputStdevs[ii]);
         if (printLevel > 2)
            printf("Parameter %3d has triangle prior distribution.\n",ii+1);
      }
      else if (pdfFlags != NULL && pdfFlags[ii] == PSUADE_PDF_BETA)
      {
         inputPDFs[ii] = (PDFBase *) new PDFBeta(inputMeans[ii],inputStdevs[ii]);
         if (printLevel > 2)
            printf("Parameter %3d has beta prior distribution.\n",ii+1);
      }
      else if (pdfFlags != NULL && pdfFlags[ii] == PSUADE_PDF_WEIBULL)
      {
         inputPDFs[ii] = (PDFBase *) new PDFWeibull(inputMeans[ii],inputStdevs[ii]);
         if (printLevel > 2)
            printf("Parameter %3d has Weibull prior distribution.\n",ii+1);
      }
      else if (pdfFlags != NULL && pdfFlags[ii] == PSUADE_PDF_GAMMA)
      {
         inputPDFs[ii] = (PDFBase *) new PDFGamma(inputMeans[ii],inputStdevs[ii]);
         if (printLevel > 2)
            printf("Parameter %3d has gamma prior distribution.\n",ii+1);
      }
      else if (pdfFlags != NULL && pdfFlags[ii] == 1000+PSUADE_PDF_NORMAL)
      {
         inputPDFs[ii] = NULL;
         if (printLevel > 2)
         {
            printf("Parameter %3d: multi-parameter normal distribution.\n",ii);
            printf("               curently not supported.\n");
            return -1.0;
         }
      }
      else if (pdfFlags != NULL && pdfFlags[ii] == 1000+PSUADE_PDF_LOGNORMAL)
      {
         inputPDFs[ii] = NULL;
         if (printLevel > 2)
         {
            printf("Parameter %3d: multi-parameter lognormal distribution.\n",ii+1);
            printf("               curently not supported.\n");
            return -1.0;
         }
      }
      else if (pdfFlags == NULL || pdfFlags[ii] == PSUADE_PDF_UNIFORM)
      {
         inputPDFs[ii] = NULL;
         if (printLevel > 2)
            printf("Parameter %3d has uniform prior distribution.\n",ii);
      }
      else if (pdfFlags != NULL && pdfFlags[ii] == PSUADE_PDF_USER)
      {
         inputPDFs[ii] = NULL;
         if (printLevel > 2)
         {
            printf("Parameter %3d: user-provided distribution currently not\n",ii+1);
            printf("               supported.\n");
            return -1.0;
         }
      }
   }

   //    ==> funcIO, freq (how often to do simulation and emulation)
   maxPts = nbins * 5;
   if (psAnaExpertMode_ == 1)
   {
      printf("Since MCMC uses many function evaluations to construct\n");
      printf("the proposal distributions, you have the option to set\n");
      printf("how many points are used to construct it in order to\n");
      printf("reduce the inference cost. The default is %d.\n", maxPts); 
      sprintf(charString, "Enter the sample size (%d - %d): ",nbins*2,nbins*10);
      maxPts = getInt(nbins*2, nbins*10, charString);
      maxPts = maxPts / nbins * nbins;
   }
   if (mode_ == 1 && ioPtr != NULL) 
   {
      funcIO = createFunctionInterface(ioPtr);
      if (nSamples == 0)
         printf("MCMC: DIRECT SIMULATION has been set up.\n");
      else
         printf("MCMC: DIRECT SIMULATION PLUS RESPONSE SURFACE have been set up.\n");
      printf("MCMC INFO: make sure simulation output is in the right form\n");
      printf("     which should not have been translated (mean) nor scaled (std)\n");
      printf("     unless you compute the error measure yourself, in which case\n");
      printf("     you should have use nOutputs=1, mean=0, and std dev=1 in the\n");
      printf("     spec file.\n");
      if (psAnaExpertMode_ == 1 && nSamples > 0)
      {
         printf("Since MCMC uses many function evaluations, you have the option\n");
         printf("to set how many times the simulator is invoked when creating the \n");
         printf("proposal distribution during each MCMC iteration. A frequency of\n");
         printf("f means that each MCMC step uses f simulator runs. These simulator\n");
         printf("runs will be supplemented with evaluations from the given response\n");
         printf("surface. The default is f=10 (if you do not know what this is, enter 10).\n"); 
         sprintf(charString, "Enter frequency f (1 - %d): ", maxPts);
         kk = getInt(1, maxPts, charString);
         freq = maxPts * nInputs / kk;
         if (freq * kk != maxPts) freq++;
      }
      else if (nSamples > 0) freq = maxPts / 10;
      else                   freq = 1;
      if (nSamples > 0)
         printf("Frequency of invoking the simulator has been set to %d\n", freq);
   }
   if (funcIO == NULL && faPtrs == NULL)
   {
      printf("MCMC ERROR: missing simulator and sample data - cannot proceed.\n");
      return PSUADE_UNDEFINED;
   }

   //    ==> faPtrs1
   if (modelFormFlag == 1)
   {
      int dleng=-999, *LHSamStates;
      int LHSnSamples, dfaType;
      int dnPerDim=16;
      double *dOneSample, expdata, simdata;
      double *LHSamOutputs, *LHSamInputs, *tSamInputs, *settings;
      Sampling *sampler;

#if 0
      double *tSamOutputs;
      FuncApprox **efaPtrs=NULL;
      LHSnSamples = nSamples;
      if (nInputs > 50) sampler = SamplingCreateFromID(PSUADE_SAMP_LHS);
      else              sampler = SamplingCreateFromID(PSUADE_SAMP_LPTAU);
      sampler->setInputBounds(nInputs, lower, upper);
      sampler->setOutputParams(1);
      sampler->setSamplingParams(LHSnSamples, 1, 1);
      sampler->initialize(0);
      LHSamInputs  = new double[LHSnSamples*nInputs];
      LHSamOutputs = new double[LHSnSamples];
      LHSamStates  = new int[LHSnSamples];
      settings = new double[nInputs];
      for (ii2 = 0; ii2 < nInputs; ii2++)
         settings[ii2] = 0.5*(lower[ii2] + upper[ii2]);
      sampler->getSamples(LHSnSamples, nInputs, 1, LHSamInputs, LHSamOutputs, 
                          LHSamStates);
      for (kk = 0; kk < LHSnSamples; kk++)
      {
         for (ii2 = 0; ii2 < dnInputs; ii2++)
            LHSamInputs[kk*nInputs+ii2] = X[kk*nInputs+ii2];
         for (ii2 = 0; ii2 < nInputs; ii2++)
            if (rsIndices != NULL && rsIndices[ii2] < 0)
               LHSamInputs[kk*nInputs+ii2] = rsValues[ii2];
      }
      delete [] LHSamStates;
      delete sampler;

      dleng = -999;
      tSamOutputs = new double[dnSamples];
      efaPtrs = new FuncApprox*[nOutputs];
      for (ii = 0; ii < nOutputs; ii++)
      {
         efaPtrs[ii] = genFA(faType, dnInputs, dnSamples);
         efaPtrs[ii]->setNPtsPerDim(dnPerDim);
         efaPtrs[ii]->setBounds(lower, upper);
         efaPtrs[ii]->setOutputLevel(0);
         for (kk = 0; kk < dnSamples; kk++)
            tSamOutputs[kk] = dSamMeans[kk*nOutputs+ii];
         efaPtrs[ii]->genNDGridData(dSamInputs,tSamOutputs,&dleng,NULL,NULL);
      }
      delete [] tSamOutputs;

      dfaType = -1;
      printf("Select response surface type for the discrepancy function.\n");
      while (dfaType < 0 || dfaType >= PSUADE_NUM_RS)
      {
         writeFAInfo();
         sprintf(charString, "Please enter your choice ? ");
         dfaType = getInt(0, PSUADE_NUM_RS-1, charString);
      }

      faPtrs1 = new FuncApprox*[nOutputs];
      dOneSample = new double[nInputs];
      tSamInputs = new double[LHSnSamples*nInputs];
      for (ii = 0; ii < nOutputs; ii++)
      {
         for (kk = 0; kk < LHSnSamples; kk++)
         {
            for (ii2 = 0; ii2 < nInputs; ii2++)
               dOneSample[ii2] = LHSamInputs[kk*nInputs+ii2];

            expdata = efaPtrs[ii]->evaluatePoint(dOneSample);
            if (modelFormScheme == 1)
            {
               if (funcIO != NULL)
                  funcIO->evaluate(kk+1,nInputs,dOneSample,1,&simdata,0);
               else
                  simdata = faPtrs[ii]->evaluatePoint(dOneSample);
            }
            else
            {
               simdata = 0.0;
               for (kk2 = 0; kk2 < LHSnSamples; kk2++)
               {
                  for (ii2 = dnInputs; ii2 < nInputs; ii2++)
                     dOneSample[ii2] = LHSamInputs[kk2*nInputs+ii2];
                  if (funcIO != NULL)
                     funcIO->evaluate(kk2+1,nInputs,dOneSample,1,&Ytemp,0);
                  else
                     Ytemp = faPtrs[ii]->evaluatePoint(dOneSample);
                  simdata += Ytemp;
               }
               simdata /= (double) LHSnSamples;
            }

            LHSamOutputs[kk] = expdata - simdata;
         }
         printf("Creating discrepancy response surface for output %d\n", ii+1);
         if (modelFormScheme == 0) faPtrs1[ii] = genFA(dfaType, dnInputs, LHSnSamples);
         else                      faPtrs1[ii] = genFA(dfaType, nInputs, LHSnSamples);
         if (faPtrs1[ii] == NULL)
         {
            printf("MCMC ERROR: cannot create discrepancy function for output %d.\n",ii+1);
            return -1.0;
         }
         faPtrs1[ii]->setNPtsPerDim(dnPerDim);
         faPtrs1[ii]->setBounds(lower, upper);
         faPtrs1[ii]->setOutputLevel(0);
         if (modelFormScheme == 0) 
         {
            for (kk = 0; kk < LHSnSamples; kk++)
               for (ii2 = 0; ii2 < dnInputs; ii2++)
                  tSamInputs[kk*dnInputs+ii2] = LHSamInputs[kk*nInputs+ii2];
            dleng = -999;
            faPtrs1[ii]->genNDGridData(tSamInputs,LHSamOutputs,&dleng,NULL,NULL);
         }
         else
         {
            dleng = -999;
            faPtrs1[ii]->genNDGridData(LHSamInputs,LHSamOutputs,&dleng,NULL,NULL);
         }
         delete efaPtrs[ii];
      }
      delete efaPtrs;
#else
      LHSnSamples = dnSamples;
      if (nInputs > 50) sampler = SamplingCreateFromID(PSUADE_SAMP_LHS);
      else              sampler = SamplingCreateFromID(PSUADE_SAMP_LPTAU);
      sampler->setInputBounds(nInputs, lower, upper);
      sampler->setOutputParams(1);
      sampler->setSamplingParams(LHSnSamples, 1, 1);
      sampler->initialize(0);
      LHSamInputs  = new double[LHSnSamples*nInputs];
      LHSamOutputs = new double[LHSnSamples];
      LHSamStates  = new int[LHSnSamples];
      settings = new double[nInputs];
      for (ii2 = 0; ii2 < nInputs; ii2++)
      {
         if (inputPDFs == NULL || (inputPDFs != NULL && inputPDFs[ii2] == NULL))
            settings[ii2] = 0.5*(lower[ii2] + upper[ii2]);
         else
            settings[ii2] = inputPDFs[ii2]->getMean();
      }
      sampler->getSamples(LHSnSamples, nInputs, 1, LHSamInputs, LHSamOutputs, 
                          LHSamStates);
      for (kk = 0; kk < LHSnSamples; kk++)
      {
         for (ii2 = 0; ii2 < dnInputs; ii2++)
            LHSamInputs[kk*nInputs+ii2] = dSamInputs[kk*dnInputs+ii2];
      }
      delete [] LHSamStates;
      delete sampler;

      dfaType = -1;
      printf("*** SELECT RESPONSE SURFACE TYPE FOR DISCREPANCY FUNCTION:\n");
      while (dfaType < 0 || dfaType >= PSUADE_NUM_RS)
      {
         writeFAInfo();
         sprintf(charString, "===> Enter your choice : ");
         dfaType = getInt(0, PSUADE_NUM_RS-1, charString);
      }

      faPtrs1 = new FuncApprox*[nOutputs];
      dOneSample = new double[nInputs];
      tSamInputs = new double[LHSnSamples*nInputs];
      int askFlag = 0;
      double *usrNominals = new double[nInputs];
      for (ii = 0; ii < nOutputs; ii++)
      {
         for (kk = 0; kk < dnSamples; kk++)
         {
            for (ii2 = 0; ii2 < nInputs; ii2++)
               dOneSample[ii2] = LHSamInputs[kk*nInputs+ii2];

            expdata = dSamMeans[kk*nOutputs+ii];
            if (modelFormScheme == 1)
            {
               if (funcIO != NULL)
                    funcIO->evaluate(kk+1,nInputs,dOneSample,1,&simdata,0);
               else simdata = faPtrs[ii]->evaluatePoint(dOneSample);
            }
            else
            {
               simdata = 0.0;
               if (LHSnSamples < 100000000)
               {
                  if (psAnaExpertMode_ == 1 && askFlag == 0)
                  {
                     printf("To create discrepancy functions, the calibration\n");
                     printf("parameters need to be set to some nominal values.\n");
                     printf("You can choose the nominal values, or it will be\n");
                     printf("set to the mid points of the ranges.\n");
                     printf("Set nomininal values yourself ? (y or n) ");
                     scanf("%s", charString);
                     fgets(lineIn,500,stdin);
                     if (charString[0] == 'y')
                     {
                        for (ii2 = dnInputs; ii2 < nInputs; ii2++)
                        {
                           printf("Input %d has lower and upper bounds = %e %e\n",
                                  ii2+1, lower[ii2], upper[ii2]);
                           sprintf(charString, "Nominal value for input %d : ",ii2+1);
                           dOneSample[ii2] = getDouble(charString);
                           usrNominals[ii2] = dOneSample[ii2];
                        }
                     }
                     else
                     {
                        for (ii2 = dnInputs; ii2 < nInputs; ii2++)
                        {
                           usrNominals[ii2] = settings[ii2];
                           printf("Discrepancy function calibration input %d default = %e \n",
                                  usrNominals[ii2]);
                        }
                     }
                     askFlag = 1;
                  }
                  if (askFlag == 0)
                  {
                     for (ii2 = dnInputs; ii2 < nInputs; ii2++)
                        dOneSample[ii2] = settings[ii2];
                  }
                  else
                  {
                     for (ii2 = dnInputs; ii2 < nInputs; ii2++)
                        dOneSample[ii2] = usrNominals[ii2];
                  }
                  if (funcIO != NULL)
                       funcIO->evaluate(kk+1,nInputs,dOneSample,1,&simdata,0);
                  else simdata = faPtrs[ii]->evaluatePoint(dOneSample);
                  if (printLevel > 4)
                  {
                     printf("Experiment %4d : ", kk+1);
                     for (ii2 = 0; ii2 < nInputs; ii2++)
                        printf("%12.4e ", dOneSample[ii2]);
                     printf(" = %12.4e (simulation = %12.4e)\n", expdata, simdata);
                  }
               }
               else
               {
                  for (kk2 = 0; kk2 < LHSnSamples; kk2++)
                  {
                     for (ii2 = dnInputs; ii2 < nInputs; ii2++)
                        dOneSample[ii2] = LHSamInputs[kk2*nInputs+ii2];
                     if (funcIO != NULL)
                        funcIO->evaluate(kk+1,nInputs,dOneSample,1,&Ytemp,0);
                     else
                        Ytemp = faPtrs[ii]->evaluatePoint(dOneSample);
                     simdata += Ytemp;
                  }
                  simdata /= (double) LHSnSamples;
               }
            }

            LHSamOutputs[kk] = expdata - simdata;
         }
         printf("Creating discrepancy response surface for output %d\n", ii+1);
         if (modelFormScheme == 0) faPtrs1[ii] = genFA(dfaType,dnInputs,iOne,LHSnSamples);
         else                      faPtrs1[ii] = genFA(dfaType,nInputs,iOne,LHSnSamples);
         if (faPtrs1[ii] == NULL)
         {
            printf("MCMC ERROR: cannot create discrepancy function for output %d.\n",ii+1);
            return -1.0;
         }
         faPtrs1[ii]->setNPtsPerDim(dnPerDim);
         faPtrs1[ii]->setBounds(lower, upper);
         faPtrs1[ii]->setOutputLevel(0);
         if (modelFormScheme == 0) 
         {
            for (kk = 0; kk < LHSnSamples; kk++)
               for (ii2 = 0; ii2 < dnInputs; ii2++)
                  tSamInputs[kk*dnInputs+ii2] = LHSamInputs[kk*nInputs+ii2];
            dleng = -999;
            faPtrs1[ii]->genNDGridData(tSamInputs,LHSamOutputs,&dleng,NULL,NULL);

            PsuadeData *dataPtr = new PsuadeData(); 
            char       **iNames;
            int        iOne=1, *states;
            if (qData.strArray_ == NULL)
            {
               iNames = new char*[dnInputs];
               for (ii2 = 0; ii2 < dnInputs; ii2++)
               {
                  iNames[ii2] = new char[100];
                  sprintf(iNames[ii2], "X%d", ii2+1);
               }
            }
            else iNames = qData.strArray_;
            dataPtr->updateInputSection(LHSnSamples, dnInputs, NULL, lower, upper,
                                        tSamInputs, iNames);
            if (qData.strArray_ == NULL)
            { 
               for (ii2 = 0; ii2 < dnInputs; ii2++) delete [] iNames[ii2];
               delete [] iNames;
            }
            states = new int[LHSnSamples];
            for (kk = 0; kk < LHSnSamples; kk++) states[kk] = 1;
            iNames = new char*[1];
            iNames[0] = new char[100];
            sprintf(iNames[0], "Y%d", ii+1);
            dataPtr->updateOutputSection(LHSnSamples, iOne, LHSamOutputs, states,
                                         iNames);
            delete [] states;
            delete [] iNames[0];
            delete [] iNames;
            dataPtr->updateMethodSection(PSUADE_SAMP_MC, LHSnSamples, 1, -1, -1);
            sprintf(charString, "DiscrepancyModel%d", ii+1);
            dataPtr->writePsuadeFile(charString, 0);
            printf("MCMC INFO: a sample (inputs and outputs) for the discrepancy\n");
            printf("           model is now in DiscrepancyModel%d.\n",ii+1);
            delete dataPtr;
         }
         else
         {
            dleng = -999;
            faPtrs1[ii]->genNDGridData(LHSamInputs,LHSamOutputs,&dleng,NULL,NULL);
            PsuadeData *dataPtr = new PsuadeData(); 
            char       **iNames;
            int        iOne=1, *states;
            if (qData.strArray_ == NULL)
            {
               iNames = new char*[nInputs];
               for (ii2 = 0; ii2 < nInputs; ii2++)
               {
                  iNames[ii2] = new char[100];
                  sprintf(iNames[ii2], "X%d", ii2+1);
               }
            }
            else iNames = qData.strArray_;

            dataPtr->updateInputSection(LHSnSamples, nInputs, NULL, lower, upper,
                                        LHSamInputs, iNames);
            if (qData.strArray_ == NULL)
            { 
               for (ii2 = 0; ii2 < nInputs; ii2++) delete [] iNames[ii2];
               delete [] iNames;
            }
            states = new int[LHSnSamples];
            for (kk = 0; kk < LHSnSamples; kk++) states[kk] = 1;
            sprintf(charString, "Y%d", ii+1);
            iNames = new char*[1];
            iNames[0] = new char[100];
            sprintf(iNames[0], "Y%d", ii+1);
            dataPtr->updateOutputSection(LHSnSamples, iOne, LHSamOutputs, states,
                                         iNames);
            delete [] states;
            delete [] iNames[0];
            delete [] iNames;
            dataPtr->updateMethodSection(PSUADE_SAMP_MC, LHSnSamples, 1, -1, -1);
            sprintf(charString, "DiscrepancyModel%d", ii+1);
            dataPtr->writePsuadeFile(charString, 0);
            printf("MCMC INFO: a sample (inputs and outputs) for the discrepancy\n");
            printf("           model is now in DiscrepancyModel%d.\n",ii+1);
            delete dataPtr;
         }
#if 0
         for (kk = 0; kk < dnSamples; kk++)
         {
            for (ii = 0; ii < dnInputs; ii++)
               dOneSample[ii] = dSamInputs[kk*dnInputs+ii];
            for (ii = dnInputs; ii < nInputs; ii++)
               dOneSample[ii] = settings[ii];
            simdata = faPtrs[0]->evaluatePoint(dOneSample);
            simdata += faPtrs1[0]->evaluatePoint(dOneSample);
            expdata = dSamMeans[kk*dnInputs];
            printf("Data sample %4d: sim vs exp data = %16.8e %16.8e\n", kk+1,
                   simdata, expdata);
         }
#endif
      }
#endif
      delete [] usrNominals;
      delete [] LHSamInputs;
      delete [] dOneSample;
      delete [] LHSamOutputs;
      delete [] settings;
      delete [] tSamInputs;
   }

   //    set up constraint filters, if any
   printEquals(0);
   printf("MCMC INFO : creating constraints, if there is any.\n");
   printf("            Constraints remove infeasible regions from the priors.\n");
   printf("            Constraints can be specified by RS constraint files.\n");
   constrPtr = new RSConstraints();
   constrPtr->genConstraints(ioPtr);
   printEquals(0);

   //    Set initial point = mid points of all inputs
   XRange  = new double[nInputs];
   XGuess  = new double[nInputs];
   XGuessD = new double[nInputs];
   for (ii = 0; ii < nInputs; ii++)
   {
      XRange[ii] = upper[ii] - lower[ii]; 
      if (inputPDFs[ii] == NULL)
           XGuess[ii] = 0.5 * (lower[ii] + upper[ii]);
      else XGuess[ii] = inputPDFs[ii]->getMean();
   }
   if (rsIndices != NULL)
   {
      for (ii = 0; ii < nInputs; ii++)
         if (rsIndices[ii] < 0) XGuess[ii] = rsValues[ii];
   }

   //    Run the algorithm (first perform burn in)
   Xmax = new double[nInputs];
   for (ii = 0; ii < nInputs; ii++) Xmax[ii] = 0.0;
   Ymax = 0.0;
   Ivec = new int[nInputs];
   Ivec[nInputs-1] = -1;
   XDist = new double[maxPts+1];
   printAsterisks(0);
   printf("MCMC PHASE 1: BURN In \n");
   fflush(stdout);
   increment = burnInSamples / nInputs;
   for (its = 0; its < increment; its++)
   {
      if ((its+1) % (increment/10) == 0)
      {
         printf("%3.0f%% ", 100.0 * ((its+1.0) / increment));
         fflush(stdout);
      }
      if (printLevel >= 4 && (its % 20) == 0)
         printf("Iteration %d of %d\n", its+1, increment);
      jj = Ivec[nInputs-1];
      generateRandomIvector(nInputs, Ivec);
      if (Ivec[0] == jj && nInputs > 1)
      {
         Ivec[0] = Ivec[nInputs-1];
         Ivec[nInputs-1] = jj;
      }
      for (kk = 0; kk < nInputs; kk++)
      {
         ii = Ivec[kk];
         if ((rsIndices == NULL ||
             (rsIndices != NULL && rsIndices[ii] >=0)) && ii >= dnInputs)
         {
            Xtemp = XGuess[ii];
            if (psGMMode_ == 1 && nOutputs > 1)
            {
               fp = fopen("proposalDistribution", "w");
               if (fp == NULL)
               {
                  printf("MCMC ERROR: cannot write to proposalDistribution file.\n");
               }
               else
               {
                  fprintf(fp,"// This file contains individual terms in the exponent\n");
                  fprintf(fp,"// of the likelihood function, i.e.\n");
                  fprintf(fp,"// E_ki = (Y_ki - m_ki)^2/sd_ki^2 in\n");
                  fprintf(fp,"// S = 1/p sum_{k=1}^p sum_{i=1)^n (Y_ki - m_ki)^2/sd_ki^2\n");
                  fprintf(fp,"// Column 1 : active input variable\n");
                  fprintf(fp,"// Column 2 : value of the active variable\n");
                  fprintf(fp,"// Column 3 : k\n");
                  fprintf(fp,"// Column 4 : E_ki for i = 1, .., n\n");
               }
            }
            nFail = 0;
            cnt = 0;
            for (jj = 0; jj <= maxPts; jj++)
            {
               XGuess[ii] = lower[ii]+jj*XRange[ii]/maxPts;

               if (nOutputs > 1)
               {
                  XDist[jj] = 0.0;
                  for (ii2 = 0; ii2 < dnInputs; ii2++) XGuessD[ii2] = XGuess[ii2];
                  for (kk2 = 0; kk2 < dnSamples; kk2++) 
                  {
                     for (ii2 = 0; ii2 < dnInputs; ii2++) 
                     {
                        XGuess[ii2] = dSamInputs[kk2*dnInputs+ii2];
                        if (XGuess[ii2] < lower[ii2] || XGuess[ii2] > upper[ii2])
                        {
                           printf("MCMC ERROR: design parameter out of bound.\n");
                           printf("Consult PSUADE developers for help.\n");
                           exit(1);
                        }
                     }
                     if (psGMMode_ == 1 && fp != NULL) 
                        fprintf(fp, "%6d %16.8e %6d ",ii+1,XGuess[ii],kk2+1);
                     for (ii2 = 0; ii2 < nOutputs; ii2++) 
                     {
                        Ytemp = faPtrs[ii2]->evaluatePoint(XGuess);
                        if (faPtrs1 != NULL && faPtrs1[ii2] != NULL)
                           Ytemp += faPtrs1[ii2]->evaluatePoint(XGuess);
                        Ytemp = pow((Ytemp-dSamMeans[kk2*nOutputs+ii2])/
                                    dSamStdevs[kk2*nOutputs+ii2],2.0);
                        XDist[jj] += Ytemp;
                        if (psGMMode_ == 1) fprintf(fp, "%16.8e ", Ytemp);
                     }
                     if (psGMMode_ == 1 && fp != NULL) fprintf(fp, "\n");
                  }
                  for (ii2 = 0; ii2 < dnInputs; ii2++) XGuess[ii2] = XGuessD[ii2];
                  XDist[jj] = exp(-0.5 * XDist[jj] / dnSamples);
               }
               else
               {
                  cnt++;
                  if (cnt >= freq) cnt = 0; 
                  XDist[jj] = 0.0;
                  for (ii2 = 0; ii2 < dnInputs; ii2++) XGuessD[ii2] = XGuess[ii2];
                  if (setCutOff == 0)
                  {
                     for (kk2 = 0; kk2 < dnSamples; kk2++) 
                     {
                        for (ii2 = 0; ii2 < dnInputs; ii2++) 
                        {
                           XGuess[ii2] = dSamInputs[kk2*dnInputs+ii2];
                           if (XGuess[ii2] < lower[ii2] || XGuess[ii2] > upper[ii2])
                           {
                              printf("MCMC ERROR: design parameter out of bound.\n");
                              printf("Consult PSUADE developers for help.\n");
                              exit(1);
                           }
                        }
                        if (funcIO == NULL)
                        {
                           Ytemp = faPtrs[0]->evaluatePoint(XGuess);
                        }
                        else
                        {
                           if (cnt == 0)
                                funcIO->evaluate(jj,nInputs,XGuess,1,&Ytemp,0);
                           else Ytemp = faPtrs[0]->evaluatePoint(XGuess);
                        }
                        if (faPtrs1 != NULL && faPtrs1[0] != NULL)
                           Ytemp += faPtrs1[0]->evaluatePoint(XGuess);

                        Ytemp = pow((Ytemp-dSamMeans[kk2])/dSamStdevs[kk2],2.0);
                        XDist[jj] += Ytemp;
                     }
                     XDist[jj] = exp(-0.5 * XDist[jj] / dnSamples);
                  }
                  else
                  {
                     if (cnt == 0)
                          funcIO->evaluate(jj,nInputs,XGuess,1,&Ytemp,0);
                     else Ytemp = faPtrs[0]->evaluatePoint(XGuess);
                     if (Ytemp < cutoffLo || Ytemp > cutoffHi) XDist[jj] = 0;
                     else XDist[jj] = 1;
                  }
                  for (ii2 = 0; ii2 < dnInputs; ii2++) XGuess[ii2] = XGuessD[ii2];
               }

               Ytemp = constrPtr->evaluate(XGuess,XDist[jj],status);
               if (status == 0)
               {
                  XDist[jj] = 0.0;
                  nFail++;
               }
               ddata = 1.0;
               if (inputPDFs != NULL)
               {
                  for (ii2 = dnInputs; ii2 < nInputs; ii2++)
                  {
                     if (inputPDFs[ii2] != NULL)
                     {
                        inputPDFs[ii2]->getPDF(ione,&XGuess[ii2],&ddata2);
                        ddata *= ddata2;
                     }
                  }
               }
               XDist[jj] *= ddata; 
               if (XDist[jj] > Ymax)
               {
                  Ymax = XDist[jj];
                  for (ii2 = 0; ii2 < nInputs; ii2++) Xmax[ii2] = XGuess[ii2];
               }
               if (jj > 0) XDist[jj] += XDist[jj-1];
            }
            if (psGMMode_ == 1 && nOutputs > 1 && fp != NULL) fclose(fp);
            
            if (printLevel >= 4)
               printf("proposal distribution max = %e\n", XDist[maxPts]);
            if (XDist[maxPts] - XDist[0] > 1.0e-15)
            {
               for (jj = 1; jj <= maxPts; jj++)
                  XDist[jj] = (XDist[jj] - XDist[0]) / (XDist[maxPts] - XDist[0]);
               XDist[0] = 0;
               XGuess[ii] = PSUADE_drand();
               index = binarySearchDble(XGuess[ii], XDist, maxPts+1);
               if (index == maxPts) index = maxPts - 1;
               if (index >= 0) ddata = (double) index;
               else
               {
                  index = - index - 1;
                  if (PABS(XDist[index]-XDist[index+1]) > 1.0e-16)
                     ddata=index +
                        (XGuess[ii]-XDist[index])/(XDist[index+1]-XDist[index]);
                  else ddata = (double) index;
               }
               XGuess[ii] = lower[ii]+ddata*XRange[ii]/maxPts;
            }
            else
            {
               XGuess[ii] = Xtemp;
               if (nFail == maxPts+1) 
                  printf("INFO : Constraints have resulted in zero distribution.\n");
               printf("MCMC iteration %7d : cannot find good value for input %d\n",
                      its+1,ii+1);
               printf("accumulative distribution = %e\n\n", XDist[maxPts]);
               printf("Recommendations: see if doing the following helps:\n");
               printf("1. shift the data mean and/or widen its std. dev. in the spec file.\n");
               printf("2. use different parameter ranges (re-baseline).\n");
               exit(1);
            }
         }
      }
   }
   printf("MCMC PHASE 1 completed\n");

   // Continue to run the Gibbs algorithm
   bins = new int*[nbins];
   for (ii = 0; ii < nbins; ii++)
   {
      bins[ii] = new int[nInputs];
      for (jj = 0; jj < nInputs; jj++) bins[ii][jj] = 0;
   }
   bins2 = new int***[nbins];
   for (jj = 0; jj < nbins; jj++)
   {
      bins2[jj] = new int**[nbins];
      for (jj2 = 0; jj2 < nbins; jj2++)
      {
         bins2[jj][jj2] = new int*[nInputs];
         for (ii = 0; ii < nInputs; ii++)
         {
            bins2[jj][jj2][ii] = new int[nInputs];
            for (ii2 = 0; ii2 < nInputs; ii2++)
               bins2[jj][jj2][ii][ii2] = 0;
         }
      }
   }
   scores = new int[nInputs];
   for (ii = 0; ii < nInputs; ii++) scores[ii] = -1;
   means = new double[nInputs];
   sigmas = new double[nInputs];
   for (ii = 0; ii < nInputs; ii++) means[ii] = sigmas[ii] = 0.0;
   NTotal = 0;
   its = 0;
   printf("MCMC PHASE 2: CREATE POSTERIORS\n");
   fflush(stdout);
   increment = maxSamples / nInputs;
   while (1)
   {
      if (((its+1) % (increment/10)) == 0)
      {
         printf("%3.0f%% ",100.0*(((its % increment)+ 1.0) / increment));
         fflush(stdout);
      }
      if (printLevel >= 4 && (its % 20) == 0)
         printf("Iteration %d of %d\n", its+1, increment);

      if (psGMMode_ == 1 && printLevel > 4) 
      {
         if (psPlotTool_ == 1) fp = fopen("psTrack.sci", "w");
         else                  fp = fopen("psTrack.m", "w");
         printf("MCMC INFO: since you have turned on printlevel = %d,\n",
                printLevel);
         printf("           the proposal distributions at every iteration\n");
         printf("           will be generated. This may create a very large\n");
         printf("           file.\n");
      }
      if (scheme == 0)
      {
         jj = Ivec[nInputs-1];
         generateRandomIvector(nInputs, Ivec);
         if (Ivec[0] == jj && nInputs > 1)
         {
            Ivec[0] = Ivec[nInputs-1];
            Ivec[nInputs-1] = jj;
         }
         for (kk = 0; kk < nInputs; kk++)
         {
            ii = Ivec[kk];
            if ((rsIndices == NULL ||
                (rsIndices != NULL && rsIndices[ii] >=0)) && ii >= dnInputs)
            {
               Xtemp = XGuess[ii];
               nFail = 0;
               cnt = 0;
               for (jj = 0; jj <= maxPts; jj++)
               {
                  XGuess[ii] = lower[ii]+jj*XRange[ii]/maxPts;
   
                  if (nOutputs > 1)
                  {
                     XDist[jj] = 0.0;
                     for (ii2 = 0; ii2 < dnInputs; ii2++) XGuessD[ii2] = XGuess[ii2];
                     for (kk2 = 0; kk2 < dnSamples; kk2++) 
                     {
                        for (ii2 = 0; ii2 < dnInputs; ii2++) 
                        {
                           XGuess[ii2] = dSamInputs[kk2*dnInputs+ii2];
                           if (XGuess[ii2] < lower[ii2] || XGuess[ii2] > upper[ii2])
                           {
                              printf("MCMC ERROR: design parameter out of bound.\n");
                              printf("Consult PSUADE developers for help.\n");
                              exit(1);
                           }
                        }
                        for (ii2 = 0; ii2 < nOutputs; ii2++) 
                        {
                           Ytemp2 = stdev = stdev2 = 0.0;
                           if (rsErrFlag == 1)
                           {
                              Ytemp = faPtrs[ii2]->evaluatePointFuzzy(XGuess,stdev);
                              if (faPtrs1 != NULL && faPtrs1[ii2] != NULL)
                                 Ytemp2 = faPtrs1[ii2]->evaluatePointFuzzy(XGuess,stdev2);
                           }
                           else
                           {
                              Ytemp = faPtrs[ii2]->evaluatePoint(XGuess);
                              if (faPtrs1 != NULL && faPtrs1[ii2] != NULL)
                                 Ytemp2 = faPtrs1[ii2]->evaluatePoint(XGuess);
                           }
                           Ytemp += Ytemp2;
                           Ytemp = pow((Ytemp-dSamMeans[kk2*nOutputs+ii2]),2.0) /
                                   (pow(dSamStdevs[kk2*nOutputs+ii2],2.0) + 
                                    stdev*stdev + stdev2 * stdev2);
                           XDist[jj] += Ytemp;
                        }
                     }
                     for (ii2 = 0; ii2 < dnInputs; ii2++) XGuess[ii2] = XGuessD[ii2];
                     XDist[jj] = exp(-0.5 * XDist[jj] / dnSamples);
                  }
                  else
                  {
                     cnt++;
                     if (cnt >= freq) cnt = 0; 
                     XDist[jj] = 0.0;
                     for (ii2 = 0; ii2 < dnInputs; ii2++) XGuessD[ii2] = XGuess[ii2];
                     if (setCutOff == 0)
                     {
                        for (kk2 = 0; kk2 < dnSamples; kk2++) 
                        {
                           for (ii2 = 0; ii2 < dnInputs; ii2++) 
                           {
                              XGuess[ii2] = dSamInputs[kk2*dnInputs+ii2];
                              if (XGuess[ii2] < lower[ii2] || XGuess[ii2] > upper[ii2])
                              {
                                 printf("MCMC ERROR: design parameter out of bound.\n");
                                 printf("Consult PSUADE developers for help.\n");
                                 exit(1);
                              }
                           }
                           stdev = Ytemp2 = stdev2 = 0.0;
                           if (rsErrFlag == 1)
                           {
                              if (funcIO == NULL || cnt > 0)
                                   Ytemp = faPtrs[0]->evaluatePointFuzzy(XGuess,stdev);
                              else funcIO->evaluate(jj,nInputs,XGuess,1,&Ytemp,0);
                              if (faPtrs1 != NULL && faPtrs1[0] != NULL)
                                 Ytemp2 = faPtrs1[0]->evaluatePointFuzzy(XGuess,stdev2);
                           }
                           else
                           {
                              if (funcIO == NULL || cnt > 0)
                                   Ytemp = faPtrs[0]->evaluatePoint(XGuess);
                              else funcIO->evaluate(jj,nInputs,XGuess,1,&Ytemp,0);
                              if (faPtrs1 != NULL && faPtrs1[0] != NULL)
                                 Ytemp2 = faPtrs1[0]->evaluatePoint(XGuess);
                           }
                           Ytemp += Ytemp2;
                           if (setCutOff == 1)
                           {
                              if (Ytemp < cutoffLo || Ytemp > cutoffHi)
                                 Ytemp = PSUADE_UNDEFINED;
                              if (Ytemp != PSUADE_UNDEFINED) Ytemp = 1.0;
                              else                           Ytemp = 0.0;
                                 XDist[jj] += Ytemp;
                           }
                           else 
                           {
                              Ytemp = pow((Ytemp-dSamMeans[kk2]),2.0) /
                                      (pow(dSamStdevs[kk2],2.0)+stdev*stdev+stdev2*stdev2);
                              XDist[jj] += Ytemp;
                           }
                        }
                        XDist[jj] = exp(-0.5 * XDist[jj] / dnSamples);
                     }
                     else
                     {
                        stdev = Ytemp2 = stdev2 = 0.0;
                        if (rsErrFlag == 1)
                        {
                           if (funcIO == NULL || cnt > 0)
                                Ytemp = faPtrs[0]->evaluatePointFuzzy(XGuess,stdev);
                           else funcIO->evaluate(jj,nInputs,XGuess,1,&Ytemp,0);
                        }
                        else
                        {
                           if (funcIO == NULL || cnt > 0)
                                Ytemp = faPtrs[0]->evaluatePoint(XGuess);
                           else funcIO->evaluate(jj,nInputs,XGuess,1,&Ytemp,0);
                        }
                        if (Ytemp < cutoffLo || Ytemp > cutoffHi) XDist[jj] = 0;
                        else XDist[jj] = 1;
                     }
                     for (ii2 = 0; ii2 < dnInputs; ii2++) XGuess[ii2] = XGuessD[ii2];
                  }

                  Ytemp = constrPtr->evaluate(XGuess,XDist[jj],status);
                  if (status == 0)
                  {
                     XDist[jj] = 0.0;
                     nFail++;
                  }

                  ddata = 1.0;
                  if (inputPDFs != NULL)
                  {
                     for (ii2 = dnInputs; ii2 < nInputs; ii2++)
                     {
                        if (inputPDFs[ii2] != NULL)
                        {
                           inputPDFs[ii2]->getPDF(ione,&XGuess[ii2],&ddata2);
                           ddata *= ddata2;
                        }
                     }
                  }
                  XDist[jj] *= ddata; 

                  if (XDist[jj] > Ymax)
                  {
                     Ymax = XDist[jj];
                     for (ii2 = 0; ii2 < nInputs; ii2++) Xmax[ii2] = XGuess[ii2];
                  }
                  if (jj > 0) XDist[jj] += XDist[jj-1];
               }

               if (XDist[maxPts] - XDist[0] > 1.0e-15)
               {
                  for (jj = 1; jj <= maxPts; jj++)
                     XDist[jj] = (XDist[jj] - XDist[0]) / (XDist[maxPts] - XDist[0]);
                  XDist[0] = 0;
                  if (psGMMode_ == 1 && printLevel > 4) 
                  {
                     fprintf(fp, "A%d = [\n", ii+1);
                     for (jj = 0; jj <= maxPts; jj++) fprintf(fp, "%e\n",XDist[jj]);
                     fprintf(fp, "];\n");
                  }
                  XGuess[ii] = PSUADE_drand();
                  index = binarySearchDble(XGuess[ii], XDist, maxPts+1);
                  if (index == maxPts) index = maxPts - 1;
                  if (index >= 0) ddata = (double) index;
                  else
                  {
                     index = - index - 1;
                     if (PABS(XDist[index]-XDist[index+1]) > 1.0e-16)
                        ddata=index +
                           (XGuess[ii]-XDist[index])/(XDist[index+1]-XDist[index]);
                     else ddata = (double) index;
                  }
                  XGuess[ii] = lower[ii]+ddata*XRange[ii]/maxPts;
               }
               else 
               {
                  XGuess[ii] = Xtemp;
                  if (printLevel > 3)
                     printf("MCMC iteration %7d : no modification in input %d\n",
                            its+1,ii+1);
                  if (nFail == maxPts+1) 
                     printf("INFO : Constraints have resulted in zero distribution.\n");
               }
            }
            if (psGMMode_ == 1 && printLevel > 4) 
            {
               if (psPlotTool_ == 1) fprintf(fp, "clf();\n");
               else                  fprintf(fp, "clf\n");
               fprintf(fp, "X%d = [\n",kk+1);
               for (jj = 0; jj <= maxPts; jj++) 
                  fprintf(fp, "%e\n", lower[kk]+jj*XRange[kk]/maxPts);
               fprintf(fp, "];\n");
               fprintf(fp, "A%d = [\n", kk+1);
               for (jj = 0; jj <= maxPts; jj++) fprintf(fp, "%e\n", XDist[jj]);
               fprintf(fp, "];\n");
               fprintf(fp, "plot(X%d,A%d)\n", kk+1, kk+1);
               fwritePlotAxes(fp);
               sprintf(charString, "Proposal Distribution for Input %d",kk+1);
               fwritePlotTitle(fp, charString);
               fprintf(fp, "disp('Press enter to advance')\n");
               fprintf(fp, "pause\n");
               printf("MCMC INFO: proposal distribution stored at iteration %d.\n",
                      its);
            }
         }
      }
      else
      {
         printf("MCMC: Metropolis Hasting not implemented yet.\n");
         exit(1);
      }

      for (ii = 0; ii < nInputs; ii++) 
      {
         ddata = (XGuess[ii] - lower[ii]) / (XRange[ii] / nbins);
         index = (int) ddata;
         if (index >= nbins) index--;
         bins[index][ii]++;
         means[ii] += XGuess[ii];
         sigmas[ii] += XGuess[ii] * XGuess[ii];
      }
      for (ii = 0; ii < nInputs; ii++) 
      {
         for (ii2 = 0; ii2 < nInputs; ii2++) 
         {
            ddata = (XGuess[ii] - lower[ii]) / (XRange[ii] / nbins);
            index = (int) ddata;
            if (index >= nbins) index--;
            ddata = (XGuess[ii2] - lower[ii2]) / (XRange[ii2] / nbins);
            index2 = (int) ddata;
            if (index2 >= nbins) index2--;
            bins2[index][index2][ii][ii2]++;
         }
      }
      NTotal++;

      if (fp2 != NULL)
      {
         fprintf(fp2, "%d ", NTotal);
         for (jj = 0; jj < nInputs; jj++)
         {
            if (rsIndices != NULL && rsIndices[jj] >= 0)
               fprintf(fp2, "%e ", XGuess[jj]);
            else if (rsIndices == NULL)
               fprintf(fp2, "%e ", XGuess[jj]);
         }
         fprintf(fp2, "\n");
      }

      if (psGMMode_ == 1 && printLevel > 4)
      {
         if (psPlotTool_ == 1) 
            printf("MCMCAnalyzer: diagnostics file is in psTrack.sci\n");
         else
            printf("MCMCAnalyzer: diagnostics file is in psTrack.m\n");
         fclose(fp);
         printf("Enter 'q' and return to continue (no more stop).\n");
         printf("Enter any other character and return to advance one step.\n");
         scanf("%s", charString);
         if (charString[0] == 'q') printLevel = 4;
         fgets(lineIn,500,stdin);
      }
      its++;
      if (its*nInputs >= maxSamples*10) break;
      if (its % increment == 0)
      {
         printf("\n\n");
         fail = 0;
         for (ii = 0; ii < nInputs; ii++) 
         {
            if ((rsIndices == NULL || (rsIndices != NULL && rsIndices[ii] >=0)) &&
                ii >= dnInputs)
            {
               index = 0;
               iMax = bins[0][ii];
               for (jj = 0; jj < nbins; jj++) 
               {
                  if (bins[jj][ii] > iMax)
                  {
                     index = jj;
                     iMax = bins[jj][ii];
                  }
               }
               if (index != scores[ii])
               {
                  fail = 1;
                  printf("MCMC: input %3d (bin %2d vs %2d) not converged yet\n",
                         ii+1, index, scores[ii]);
               }
               scores[ii] = index;
               if (NTotal == 0)
               {
                  printf("MCMC: input %3d statistics (mean) unavailable\n",ii+1);
                  break;
               }
               else
                  printf("MCMC: input %3d statistics (mean) = %e\n",
                         ii+1, means[ii]/NTotal); 
            }
         }
         if (ii != nInputs) break;
         for (ii = 0; ii < nInputs; ii++) 
            if ((rsIndices == NULL || (rsIndices != NULL && rsIndices[ii] >=0)) && 
                ii >= dnInputs)
               printf("MCMC: input %3d  at peak of likelihood = %e\n", ii+1,
                      Xmax[ii]);
         for (ii2 = 0; ii2 < nOutputs; ii2++) 
         {
            if (funcIO != NULL)
               funcIO->evaluate(maxPts+1,nInputs,XGuess,1,&Ytemp,0);
            else
               Ytemp = faPtrs[ii2]->evaluatePoint(XGuess);
            if (faPtrs1 != NULL && faPtrs1[ii2] != NULL)
               Ytemp += faPtrs1[ii2]->evaluatePoint(XGuess);
            printf("MCMC: output %3d at peak of likelihood = %e\n", ii2+1, Ytemp);
         }  
         printf("MCMC: likelihood at peak = %e\n", Ymax);
         if (fail == 0) break;
         printf("MCMC: need another %d runs.\n", increment*nInputs);
         if (psPlotTool_ == 1) fp = fopen("scilabmcmc.sci", "w");
         else                  fp = fopen("matlabmcmc.m", "w");
         ngrid = (int) (pow(1.0*nPlots, 0.5001) + 1.0); 
         fwritePlotCLF(fp);
         for (kk = 0; kk < nPlots; kk++)
         {
            ii = plotIndices[kk];
            fprintf(fp, "subplot(%d,%d,%d)\n", ngrid, ngrid, kk+1);
            fprintf(fp, "X%d = [\n", ii+1);
            for (jj = 0; jj < nbins; jj++)
               fprintf(fp, "%e ",XRange[ii]/nbins*(jj+0.5)+lower[ii]);
            fprintf(fp, "];\n");
            fprintf(fp, "N%d = [\n", ii+1);
            sumBins = 0;
            for (jj = 0; jj < nbins; jj++) sumBins += bins[jj][ii];
            if (sumBins == 0) sumBins = 1;
            for (jj = 0; jj < nbins; jj++)
               fprintf(fp, "%e ", (double) bins[jj][ii]/(double) sumBins);
            fprintf(fp, "];\n");
            if (psPlotTool_ == 1)
            {
               fprintf(fp, "bar(X%d, N%d, 1.0);\n", ii+1, ii+1);
               fprintf(fp, "xmin = min(X%d);\n", ii+1);
               fprintf(fp, "xmax = max(X%d);\n", ii+1);
               fprintf(fp, "xwid = xmax - xmin;\n");
               fprintf(fp, "xmin = xmin - 0.5 * xwid / %d;\n", nbins);
               fprintf(fp, "xmax = xmax + 0.5 * xwid / %d;\n", nbins);
               fprintf(fp, "ymax = max(N%d);\n", ii+1);
               fprintf(fp, "e = gce();\n");
               fprintf(fp, "e.children.thickness = 2;\n");
               fprintf(fp, "e.children.foreground = 0;\n");
               fprintf(fp, "e.children.background = 2;\n");
               fprintf(fp, "a = gca();\n");
               fprintf(fp, "a.data_bounds=[xmin,0;xmax,ymax];\n");
            }
            else
            {
               fprintf(fp, "bar(X%d,N%d,1.0)\n", ii+1, ii+1);
               fprintf(fp, "xmin = min(X%d);\n", ii+1);
               fprintf(fp, "xmax = max(X%d);\n", ii+1);
               fprintf(fp, "xwid = xmax - xmin;\n");
               fprintf(fp, "xmin = xmin - 0.5 * xwid / %d;\n", nbins);
               fprintf(fp, "xmax = xmax + 0.5 * xwid / %d;\n", nbins);
               fprintf(fp, "ymax = max(N%d);\n", ii+1);
               fprintf(fp, "axis([xmin xmax 0 ymax])\n");
            }
            fwritePlotXLabel(fp, "Input Numbers");
            fwritePlotYLabel(fp, "Probabilities");
            fwritePlotAxesNoGrid(fp);
         }
         fclose(fp);
         if (psPlotTool_ == 1)
            printf("MCMC: scilabmcmc.sci file has been created.\n");
         else
            printf("MCMC: matlabmcmc.m file has been created.\n");
      }
      fp = fopen("psuade_stop", "r");
      if (fp != NULL)
      {
         fclose(fp);
         break;
      }
   }
   printf("MCMC PHASE 2 completed\n");
   if (NTotal == 0)
   {
      printf("MCMC FAILS: close to zero probability distribution detected.\n");
      printf("Recommendations: see the following helps:\n");
      printf("1. widen the data standard deviation in the spec file.\n");
      printf("2. use different parameter ranges (re-baseline).\n");
   }
   else  
   {
      ngrid = (int) (pow(1.0*nPlots, 0.5001) + 1.0); 
      if (psPlotTool_ == 1) fp = fopen("scilabmcmc.sci", "w");
      else                  fp = fopen("matlabmcmc.m", "w");
      fwritePlotCLF(fp);
      for (kk = 0; kk < nPlots; kk++)
      {
         ii = plotIndices[kk];
         fprintf(fp, "subplot(%d,%d,%d)\n", ngrid, ngrid, kk+1);
         fprintf(fp, "X%d = [\n", ii+1);
         for (jj = 0; jj < nbins; jj++)
            fprintf(fp, "%e ", XRange[ii]/nbins*(jj+0.5)+lower[ii]);
         fprintf(fp, "];\n");
         fprintf(fp, "N%d = [\n", ii+1);
         sumBins = 0;
         for (jj = 0; jj < nbins; jj++) sumBins += bins[jj][ii];
         if (sumBins == 0) sumBins = 1;
         for (jj = 0; jj < nbins; jj++)
            fprintf(fp, "%e ", (double) bins[jj][ii]/(double) sumBins);
         fprintf(fp, "];\n");
         if (psPlotTool_ == 1)
         {
            fprintf(fp, "bar(X%d, N%d, 1.0);\n", ii+1, ii+1);
            fprintf(fp, "xmin = min(X%d);\n", ii+1);
            fprintf(fp, "xmax = max(X%d);\n", ii+1);
            fprintf(fp, "xwid = xmax - xmin;\n");
            fprintf(fp, "xmin = xmin - 0.5 * xwid / %d;\n", nbins);
            fprintf(fp, "xmax = xmax + 0.5 * xwid / %d;\n", nbins);
            fprintf(fp, "ymax = max(N%d);\n", ii+1);
            fprintf(fp, "e = gce();\n");
            fprintf(fp, "e.children.thickness = 2;\n");
            fprintf(fp, "e.children.foreground = 0;\n");
            fprintf(fp, "e.children.background = 2;\n");
            fprintf(fp, "a = gca();\n");
            fprintf(fp, "a.data_bounds=[xmin,0;xmax,ymax];\n");
         }
         else
         {
            fprintf(fp, "bar(X%d,N%d,1.0)\n", ii+1, ii+1);
            fprintf(fp, "xmin = min(X%d);\n", ii+1);
            fprintf(fp, "xmax = max(X%d);\n", ii+1);
            fprintf(fp, "xwid = xmax - xmin;\n");
            fprintf(fp, "xmin = xmin - 0.5 * xwid / %d;\n", nbins);
            fprintf(fp, "xmax = xmax + 0.5 * xwid / %d;\n", nbins);
            fprintf(fp, "ymax = max(N%d);\n", ii+1);
            fprintf(fp, "axis([xmin xmax 0 ymax])\n");
         }
         if (qData.strArray_ != NULL)
              sprintf(charString, "%s", qData.strArray_[ii]);
         else sprintf(charString, "Input %d", ii+1);
         fwritePlotXLabel(fp, charString);
         fwritePlotYLabel(fp, "Probabilities");
         fwritePlotAxesNoGrid(fp);
      }
      fclose(fp);
      printAsterisks(0);
      if (psPlotTool_ == 1)
         printf("MCMC: final scilabmcmc.sci file has been created.\n");
      else
         printf("MCMC: final matlabmcmc.m file has been created.\n");
      printEquals(0);

      if (psPlotTool_ == 1) fp = fopen("scilabmcmc2.sci", "w");
      else                  fp = fopen("matlabmcmc2.m", "w");
      fwritePlotCLF(fp);
      fprintf(fp, "ns = 0;\n");
      for (kk = 0; kk < nPlots; kk++)
      {
         ii = plotIndices[kk];
         for (kk2 = 0; kk2 <= kk; kk2++)
         {
            ii2 = plotIndices[kk2];
            fprintf(fp, "subplot(%d,%d,%d)\n",nPlots,nPlots,kk2*nPlots+kk+1);
            if (ii == ii2)
            {
               fprintf(fp, "X%d = [\n", ii+1);
               for (jj = 0; jj < nbins; jj++)
                  fprintf(fp, "%e ", XRange[ii]/nbins*(jj+0.5)+lower[ii]);
               fprintf(fp, "];\n");
               fprintf(fp, "N%d = [\n", ii+1);
               sumBins = 0;
               for (jj = 0; jj < nbins; jj++) sumBins += bins[jj][ii];
               if (sumBins == 0) sumBins = 1;
               for (jj = 0; jj < nbins; jj++)
                  fprintf(fp, "%e ", (double) bins[jj][ii]/(double) sumBins);
               fprintf(fp, "];\n");
            
               if (psPlotTool_ == 1)
               {
                  fprintf(fp, "bar(X%d, N%d, 1.0);\n",ii+1,ii+1);
                  fprintf(fp, "xmin = min(X%d);\n", ii+1);
                  fprintf(fp, "xmax = max(X%d);\n", ii+1);
                  fprintf(fp, "xwid = xmax - xmin;\n");
                  fprintf(fp, "xmin = xmin - 0.5 * xwid / %d;\n", nbins);
                  fprintf(fp, "xmax = xmax + 0.5 * xwid / %d;\n", nbins);
                  fprintf(fp, "ymax = max(N%d);\n", ii+1);
                  fprintf(fp, "e = gce();\n");
                  fprintf(fp, "e.children.thickness = 2;\n");
                  fprintf(fp, "e.children.foreground = 0;\n");
                  fprintf(fp, "e.children.background = 2;\n");
                  fprintf(fp, "a = gca();\n");
                  fprintf(fp, "a.data_bounds=[xmin,0;xmax,ymax];\n");
               }
               else
               {
                  fprintf(fp, "bar(X%d,N%d,1.0)\n", ii+1, ii+1);
                  fprintf(fp, "xmin = min(X%d);\n", ii+1);
                  fprintf(fp, "xmax = max(X%d);\n", ii+1);
                  fprintf(fp, "xwid = xmax - xmin;\n");
                  fprintf(fp, "xmin = xmin - 0.5 * xwid / %d;\n", nbins);
                  fprintf(fp, "xmax = xmax + 0.5 * xwid / %d;\n", nbins);
                  fprintf(fp, "ymax = max(N%d);\n", ii+1);
                  fprintf(fp, "axis([xmin xmax 0 ymax])\n");
               }
               if (qData.strArray_ != NULL)
                    sprintf(charString, "%s", qData.strArray_[ii]);
               else sprintf(charString, "Input %d", ii+1);
               fwritePlotXLabel(fp, charString);
               fwritePlotYLabel(fp, "Probabilities");
               fwritePlotAxes(fp);
            }
            else
            {
               fprintf(fp, "Y%d = [\n", ii+1);
               for (jj = 0; jj < nbins; jj++)
                  fprintf(fp, "%e ", XRange[ii]/nbins*(jj+0.5)+lower[ii]);
               fprintf(fp, "];\n");
               fprintf(fp, "X%d = [\n", ii2+1);
               for (jj = 0; jj < nbins; jj++)
                  fprintf(fp, "%e ", XRange[ii2]/nbins*(jj+0.5)+lower[ii2]);
               fprintf(fp, "];\n");
               fprintf(fp, "N%d_%d = [\n", ii+1, ii2+1);
               for (jj = 0; jj < nbins; jj++)
               {
                  for (jj2 = 0; jj2 < nbins; jj2++)
                     fprintf(fp, "%d ", bins2[jj][jj2][ii][ii2]);
                  fprintf(fp, "\n");
               }
               fprintf(fp, "];\n");
               fprintf(fp, "n = length(X%d);\n", ii2+1);
               fprintf(fp, "XT = X%d;\n", ii2+1);
               fprintf(fp, "YT = Y%d;\n", ii+1);
               fprintf(fp, "HX = (XT(n) - XT(1)) / (n-1);\n");
               fprintf(fp, "HY = (YT(n) - YT(1)) / (n-1);\n");
               fprintf(fp, "ZZ = N%d_%d;\n", ii+1, ii2+1);
               fprintf(fp, "for kk = 1 : ns\n");
               fprintf(fp, "ZZ1 = ZZ;\n");
               fprintf(fp, "for ii = 2 : n-1\n");
               fprintf(fp, "   for jj = 2 : n-1\n");
               fprintf(fp, "      ZZ(ii,jj) = ZZ(ii,jj) + ZZ1(ii+1,jj);\n");
               fprintf(fp, "      ZZ(ii,jj) = ZZ(ii,jj) + ZZ1(ii-1,jj);\n");
               fprintf(fp, "      ZZ(ii,jj) = ZZ(ii,jj) + ZZ1(ii,jj+1);\n");
               fprintf(fp, "      ZZ(ii,jj) = ZZ(ii,jj) + ZZ1(ii,jj-1);\n");
               fprintf(fp, "      ZZ(ii,jj) = ZZ(ii,jj) + ZZ1(ii+1,jj+1);\n");
               fprintf(fp, "      ZZ(ii,jj) = ZZ(ii,jj) + ZZ1(ii-1,jj-1);\n");
               fprintf(fp, "      ZZ(ii,jj) = ZZ(ii,jj) + ZZ1(ii-1,jj+1);\n");
               fprintf(fp, "      ZZ(ii,jj) = ZZ(ii,jj) + ZZ1(ii+1,jj-1);\n");
               fprintf(fp, "      ZZ(ii,jj) = ZZ(ii,jj) / 9;\n");
               fprintf(fp, "   end;\n");
               fprintf(fp, "end;\n");
               fprintf(fp, "end;\n");
               if (psPlotTool_ == 1)
               {
                  fprintf(fp, "XX = [XT(1):HX:XT(n)];\n");
                  fprintf(fp, "YY = [YT(1):HY:YT(n)];\n");
                  fprintf(fp, "DD = splin2d(XX,YY,ZZ);\n");
                  fprintf(fp, "HX = 0.01 * (XT(n) - XT(1));\n");
                  fprintf(fp, "HY = 0.01 * (YT(n) - YT(1));\n");
                  fprintf(fp, "X2 = [XT(1):HX:XT(n)];\n");
                  fprintf(fp, "Y2 = [YT(1):HY:YT(n)];\n");
                  fprintf(fp, "[XI, YI] = ndgrid(X2, Y2);\n");
                  fprintf(fp, "disp('interpolation')\n");
                  fprintf(fp, "ZI =interp2d(XI, YI, XX, YY, DD, \"natural\");\n");
                  fprintf(fp, "disp('interpolation done')\n");
                  fprintf(fp, "ZB = ZI;\n");
                  fprintf(fp, "nX = length(X2);\n");
                  fprintf(fp, "nY = length(Y2);\n");
                  fprintf(fp, "for ii = 1 : nX\n");
                  fprintf(fp, "for jj = 1 : nY\n");
                  fprintf(fp, "ZI(ii,jj) = ZB(ii,nY-jj+1);\n");
                  fprintf(fp, "end;\n");
                  fprintf(fp, "end;\n");
                  fprintf(fp,"zmax = max(max(ZI));\n");
                  fprintf(fp,"zmin = min(min(ZI)) / zmax;\n");
                  fprintf(fp,"ZI   = ZI / zmax;\n");
                  fprintf(fp,"zmax = 1;\n");
                  fprintf(fp,"Matplot1((ZI-zmin)/(zmax-zmin)*64,[%e,%e,%e,%e])\n",
                          lower[ii2], lower[ii], upper[ii2], upper[ii]);
                  fprintf(fp, "xset(\"colormap\",jetcolormap(64));\n");
                  fprintf(fp, "colorbar(zmin,zmax);\n");
                  fwritePlotAxesNoGrid(fp);
                  fprintf(fp, "contour2d(X2,Y2,ZB,5,rect=[%e,%e,%e,%e]);\n",
                          lower[ii2], lower[ii], upper[ii2], upper[ii]);
                  fprintf(fp, "xset(\"fpf\",\" \");\n");
               }
               else
               {
                  fprintf(fp,"[YY,XX]=meshgrid(XT(1):HX:XT(n),YT(1):HY:YT(n));\n");
                  fprintf(fp, "HX = 0.01 * (XT(n) - XT(1));\n");
                  fprintf(fp, "HY = 0.01 * (YT(n) - YT(1));\n");
                  fprintf(fp,"[YI,XI]=meshgrid(XT(1):HX:XT(n),YT(1):HY:YT(n));\n");
                  fprintf(fp, "ZI=interp2(YY, XX, ZZ, YI, XI, 'spline');\n");
                  fprintf(fp, "pcolor(XI,YI,ZI)\n");
                  fprintf(fp, "shading interp\n");
                  fprintf(fp, "hold on\n");
                  fprintf(fp, "contour(XI,YI,ZI,5,'k')\n");
                  fwritePlotAxesNoGrid(fp);
               }
               if (qData.strArray_ != NULL)
                    sprintf(charString, "%s", qData.strArray_[ii]);
               else sprintf(charString, "Input %d", ii+1);
               fwritePlotXLabel(fp, charString);
               if (qData.strArray_ != NULL)
                    sprintf(charString, "%s", qData.strArray_[ii2]);
               else sprintf(charString, "Input %d", ii2+1);
               fwritePlotYLabel(fp, charString);
            }
         }
      }
      fclose(fp);
      if (psPlotTool_ == 1)
         printf("MCMC: scilabmcmc2.sci file (2-input analysis) is ready.\n");
      else
         printf("MCMC: matlabmcmc2.m file (2-input analysis) is ready.\n");
      printAsterisks(0);
   }

   //    clean up
   if (inputPDFs != NULL)
   {
      for (ii = 0; ii < nInputs; ii++)
         if (inputPDFs[ii] != NULL) delete inputPDFs[ii];
      delete [] inputPDFs;
   }
   delete [] scores;
   delete [] Xmax;
   delete [] means;
   delete [] sigmas;
   delete [] XDist;
   delete [] XGuess;
   delete [] XGuessD;
   delete [] XRange;
   if (plotIndices != NULL) delete [] plotIndices;
   if (dSamMeans != NULL) delete [] dSamMeans;
   if (dSamStdevs != NULL) delete [] dSamStdevs;
   if (faPtrs != NULL)
   {
      for (ii = 0; ii < nOutputs; ii++) 
         if (faPtrs[ii] != NULL) delete faPtrs[ii];
      delete [] faPtrs;
   }
   if (faPtrs1 != NULL)
   {
      for (ii = 0; ii < nOutputs; ii++) 
         if (faPtrs1[ii] != NULL) delete faPtrs1[ii];
      delete [] faPtrs1;
   }
   for (ii = 0; ii < nbins; ii++) delete [] bins[ii];
   delete [] bins;
   for (jj = 0; jj < nbins; jj++)
   {
      for (jj2 = 0; jj2 < nbins; jj2++)
      {
         for (ii = 0; ii < nInputs; ii++) delete [] bins2[jj][jj2][ii];
         delete [] bins2[jj][jj2];
      }
      delete [] bins2[jj];
   }
   delete [] bins2;
   if (dSamInputs != NULL) delete [] dSamInputs;
   if (rsIndices != NULL) delete [] rsIndices;
   if (rsValues  != NULL) delete [] rsValues;
   if (fp2 != NULL)
   {
      fprintf(fp2, "PSUADE_END\n");
      fclose(fp2);
      printf("MCMC: check the 'MCMCPostSample' file for a posterior sample.\n");
   }
   delete constrPtr;
   delete [] Ivec;
   return 0.0;
}

// ************************************************************************
// set parameters
// ------------------------------------------------------------------------
int MCMCAnalyzer::setParams(int argc, char **argv)
{
   char  *request = (char *) argv[0];
   if (!strcmp(request, "setsim")) mode_ = 1;
   else
   {
      printf("MCMCAnalyzer ERROR: setParams - not valid.\n");
      exit(1);
   }
   return 0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
MCMCAnalyzer& MCMCAnalyzer::operator=(const MCMCAnalyzer &)
{
   printf("MCMCAnalyzer operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

