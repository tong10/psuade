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
#include <algorithm>
#include <assert.h>

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
#include "PrintingTS.h"
#include "TwoSampleAnalyzer.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
MCMCAnalyzer::MCMCAnalyzer() : Analyzer(), nInputs_(0), nOutputs_(0), means_(0),
                        sigmas_(0), mostLikelyInput_(0), mostLikelyOutput_(0)
{
   setName("MCMC");
   mode_ = 0; // RS mode (mode = 1 : simulator-based mode)
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
MCMCAnalyzer::~MCMCAnalyzer()
{
   if(means_) delete[] means_;
   if(sigmas_) delete[] sigmas_;
   if(mostLikelyInput_) delete[] mostLikelyInput_;
   if(mostLikelyOutput_) delete[] mostLikelyOutput_;
}

// ************************************************************************
// perform MCMC analysis 
// ------------------------------------------------------------------------
double MCMCAnalyzer::analyze(aData &adata)
{
   int    ii, ii2, jj, kk, status, cnt, iOne=1, iZero=0;
   int    nInputs, nOutputs, nSamples, maxPts=257, nbins, printLevel, faType;
   int    maxSamples, burnInSamples, modelFormFlag=0, numChains=3, nChainGood=0;
   int    nPlots, *plotIndices=NULL, *pdfFlags, *designParams=NULL;
   int    *rsIndices=NULL, freq=1, dnSamples=0, dnInputs=0, rsErrFlag=0;
   double *dSamInputs=NULL, *dSamMeans=NULL, *dSamStdevs=NULL;
   double *X=NULL, *Y=NULL, *YY=NULL, *lower=NULL, *upper=NULL;
   double *inputMeans=NULL, *inputStdevs=NULL, *rsValues=NULL;
   double *discFuncConstantMeans=NULL, *discFuncConstantStds=NULL;
   double psrfThreshold=1.05, genPosteriors=0;
   char   lineIn[1001], charString[1001], cfname[501], *rsFile=NULL;
   FILE   *fp=NULL;
   string sfname, iline;
   size_t compFlag;
   ifstream   ifile;
   pData      pPtr, qData, pOutputs;
   FuncApprox **faPtrs=NULL, **faPtrs1=NULL;
   PDFBase    **inputPDFs;
   FunctionInterface *funcIO=NULL;
   PsuadeData *ioPtr=NULL;
   RSConstraints *constrPtr;

   printLevel  = adata.printLevel_;
   printAsterisks(PL_INFO, 0);
   printOutTS(PL_INFO,"*                     MCMC Optimizer\n");
   printEquals(PL_INFO, 0);
   if (printLevel > 0)
   {
   printOutTS(PL_INFO,"TO GAIN ACCESS TO DIFFERENT OPTIONS: TURN ON\n\n");
   printOutTS(PL_INFO," * ana_expert to finetune MCMC parameters, \n");
   printOutTS(PL_INFO,"   (e.g. sample size for burn-in can be adjusted).\n");
   printOutTS(PL_INFO," * rs_expert to finetune response surface used in MCMC,\n");
   printOutTS(PL_INFO," * printlevel 3 to display more diagnostic information.\n");
   printOutTS(PL_INFO," * printlevel 4 to display more iteration information.\n");
   printOutTS(PL_INFO," * printlevel >=5 reserved only for expert diagnostics only.\n");
   printDashes(PL_INFO,0);
   printOutTS(PL_INFO,"FEATURES AVAILABLE IN THE CURRENT VERSION OF MCMC:\n");
   printOutTS(PL_INFO," * Support different prior distributions (default: uniform)\n");
   printOutTS(PL_INFO,"   - set ana_expert to use other than uniform distributions.\n");
   printOutTS(PL_INFO," * Support multiple outputs (likelihood from multiple outputs)\n");
   printOutTS(PL_INFO," * For some response surfaces such as polynomial and Legendre\n");
   printOutTS(PL_INFO,"   polynomials, bootstrapped MARS, and Gaussian process (GP1),\n");
   printOutTS(PL_INFO,"   the fitting errors (response surface errors) may be used in\n");
   printOutTS(PL_INFO,"   the likelihood functions. There are two types of such errors:\n");
   printOutTS(PL_INFO,"   1. fitting errors (errors incurred in interpolation), and\n");
   printOutTS(PL_INFO,"   2. model form errors (in this case a discrepancy model is used.)\n");
   printOutTS(PL_INFO,"      - use ana_expert to select these options\n");
   printOutTS(PL_INFO," * Some input parameters may be designated as design parameters\n");
   printOutTS(PL_INFO,"   - they will have to be specified in the data spec file\n");
   printOutTS(PL_INFO," * Some input parameters may be disabled (set to default values)\n");
   printOutTS(PL_INFO,"   - in case these parameters are not to be calibrated\n");
   printOutTS(PL_INFO,"   - use rs_index_file in PSUADE data file's ANALYSIS section\n");
   printOutTS(PL_INFO,"     to select them and to set them to default values.\n");
   printOutTS(PL_INFO,"   - this feature is not available with discrepancy modeling.\n");
   printOutTS(PL_INFO," * A sample may be generated from the posteriors.\n");
   printOutTS(PL_INFO,"   - use ana_expert to be able to select this option\n");
   printOutTS(PL_INFO," * This analysis can be terminated gracefully by creating a file\n");
   printOutTS(PL_INFO,"   named 'psuade_stop' in the same directory during execution.\n");
   printOutTS(PL_INFO,"   (in case it takes too long).\n");
   if (psAnaExpertMode_ == 1)
   {
      printOutTS(PL_INFO," * For multi-modal posteriors, it may be necessary to use a large\n");
      printOutTS(PL_INFO,"   number of chains. Turn on ana_expert mode to enable this.\n");
      printOutTS(PL_INFO," * You can turn on/off master mode on the fly by creating the \n");
      printOutTS(PL_INFO,"   'psuade_master' or 'psuade_nomaster' file in the run directory.\n");
   }
   if (psGMMode_ == 1)
   {
      printOutTS(PL_INFO," * In GM mode, MCMCDistTrack.m will be created at each point\n");
      printOutTS(PL_INFO,"   to give a snapshot of the proposal distribution scaning\n");
      printOutTS(PL_INFO,"   the selected input.\n");
      printOutTS(PL_INFO," * You can turn on/off GM mode on the fly by creating the \n");
      printOutTS(PL_INFO,"   'psuade_gm' or 'psuade_nogm' file in the run directory.\n");
   }
   if (psMasterMode_ == 1)
   {
      printOutTS(PL_INFO," * In Master mode, MCMCChainHistogram.m will be created to give\n");
      printOutTS(PL_INFO,"   the current histograms of input values.\n");
   }
   printEquals(PL_INFO, 0);
   }

   nInputs     = adata.nInputs_;
   nInputs_    = nInputs;
   nOutputs    = adata.nOutputs_;
   nOutputs_   = nOutputs;
   nSamples    = adata.nSamples_;
   X           = adata.sampleInputs_;
   Y           = adata.sampleOutputs_;
   lower       = adata.iLowerB_;
   upper       = adata.iUpperB_;
   pdfFlags    = adata.inputPDFs_;
   inputMeans  = adata.inputMeans_;
   inputStdevs = adata.inputStdevs_;
   ioPtr       = adata.ioPtr_;
   if (ioPtr != NULL) ioPtr->getParameter("input_names", qData);

   if (psAnaExpertMode_ == 0)
   {
      if (pdfFlags != NULL)
      {
         ii2 = 0;
         for (ii = 0; ii < nInputs; ii++) ii2 += pdfFlags[ii];
         if (ii2 > 0)
         {
            printOutTS(PL_ERROR,"MCMC ERROR: non-uniform prior distribution currently\n");
            printOutTS(PL_ERROR,"     supported only in ana_expert mode.\n");
            printOutTS(PL_ERROR,"     CONTINUE with uniform priors.\n");
            for (ii = 0; ii < nInputs; ii++) pdfFlags[ii] = 0;
         }
      }
   }

   // error checking
   if (nInputs <= 0 || nOutputs <= 0)
   {
      printOutTS(PL_ERROR,"MCMC ERROR: invalid nInputs or nOutputs.\n");
      printOutTS(PL_ERROR,"    nInputs  = %d\n", nInputs);
      printOutTS(PL_ERROR,"    nOutputs = %d\n", nOutputs);
      return PSUADE_UNDEFINED;
   }
   if (mode_ == 1 && ioPtr == NULL)
   {
      printOutTS(PL_ERROR,"MCMC ERROR: Direct simulation mode is requested but no\n");
      printOutTS(PL_ERROR,"     information on the simulator is given (missing\n");
      printOutTS(PL_ERROR,"     driver program in your PSUADE data file).\n");
      return PSUADE_UNDEFINED;
   }
   if (nOutputs > 1 && mode_ == 1 && ioPtr != NULL)
   {
      printOutTS(PL_ERROR,"MCMC ERROR: Direct simulation mode does not support nOutputs>1.\n");
      printOutTS(PL_ERROR,"     Only response surface mode supports nOutputs > 1.\n");
      printOutTS(PL_ERROR,"     (nOutputs > 1 can be handled in your program by\n");
      printOutTS(PL_ERROR,"     computing the likelihood function yourself and\n");
      printOutTS(PL_ERROR,"     then use mean of zero and std dev of 1).\n");
      return PSUADE_UNDEFINED;
   }
   status = 0;
   for (ii = 0; ii < nSamples*nOutputs; ii++)
      if (Y[ii] > 0.9*PSUADE_UNDEFINED) status = 1;
   if (status == 1)
   {
      printOutTS(PL_ERROR,"MCMC ERROR: Some outputs are undefined (having value larger\n");
      printOutTS(PL_ERROR,"     than or equal to 1e35). Prune the undefined sample points\n");
      printOutTS(PL_ERROR,"     first and re-run.\n");
      return PSUADE_UNDEFINED;
   }
   fp = fopen("psuade_stop", "r");
   if (fp != NULL)
   {
      fclose(fp);
      fp = NULL;
      printOutTS(PL_ERROR,"MCMC INFO: psuade_stop FILE FOUND. WILL BE REMOVED\n");
      strcpy(charString, "psuade_stop");
      unlink(charString);
   }
   fp = fopen("psuade_master", "r");
   if (fp != NULL)
   {
      fclose(fp);
      fp = NULL;
      printOutTS(PL_ERROR,"MCMC INFO: psuade_master FILE FOUND. WILL BE REMOVED\n");
      strcpy(charString, "psuade_master");
      unlink(charString);
   }
   fp = fopen("psuade_nomaster", "r");
   if (fp != NULL)
   {
      fclose(fp);
      fp = NULL;
      printOutTS(PL_ERROR,"MCMC INFO: psuade_nomaster FILE FOUND. WILL BE REMOVED\n");
      strcpy(charString, "psuade_nomaster");
      unlink(charString);
   }
   fp = fopen("psuade_gm", "r");
   if (fp != NULL)
   {
      fclose(fp);
      fp = NULL;
      printOutTS(PL_ERROR,"MCMC INFO: psuade_gm FILE FOUND. WILL BE REMOVED\n");
      strcpy(charString, "psuade_gm");
      unlink(charString);
   }
   fp = fopen("psuade_nogm", "r");
   if (fp != NULL)
   {
      fclose(fp);
      fp = NULL;
      printOutTS(PL_ERROR,"MCMC INFO: psuade_nogm FILE FOUND. WILL BE REMOVED\n");
      strcpy(charString, "psuade_nogm");
      unlink(charString);
   }

   // clean up
   if (means_) delete[] means_;
   if (sigmas_) delete[] sigmas_;
   if (mostLikelyInput_) delete[] mostLikelyInput_;
   if (mostLikelyOutput_) delete[] mostLikelyOutput_;
   means_ = NULL;
   sigmas_ = NULL;
   mostLikelyInput_ = NULL;
   mostLikelyOutput_ = NULL;

   // get experimental data information from the spec file
   if (printLevel > 0)
   {
   printOutTS(PL_INFO,"*** NEED DATA TO CREATE GAUSSIAN LIKELIHOOD FUNCTION: \n\n");
   printOutTS(PL_INFO,"MCMC is going to create a Gaussian likelihood function.\n");
   printOutTS(PL_INFO,"Please provide a data file containing design parameter\n");
   printOutTS(PL_INFO,"values, mean, and standard deviation of the observation\n");
   printOutTS(PL_INFO,"data for each output.\n");
   printOutTS(PL_INFO,"NOTE: the design parameters are useful only if you wish to\n");
   printOutTS(PL_INFO,"      include a discrepancy function (if you don't know what\n");
   printOutTS(PL_INFO,"      this is, simply do not specify any design parameters).\n");
   printOutTS(PL_INFO,"IMPORTANT: IF m DESIGN PARAMETERS ARE SPECIFIED, YOU NEED\n");
   printOutTS(PL_INFO,"      TO SPECIFY WHICH ONES (IN THE INPUT SPACE) THEY ARE.\n");
   printOutTS(PL_INFO,"      THESE DESIGN PARAMETERS WILL BE EXCLUDED FROM BEING\n");
   printOutTS(PL_INFO,"      CALIBRATION PARAMETERS.\n");
   printDashes(PL_INFO, 0);
   printOutTS(PL_INFO,"*** THE FORMAT OF THIS DATA FILE IS: (O1 means output 1, \n");
   printOutTS(PL_INFO,"    m - no. of design parameters, p - no. of experiments): \n");
   printOutTS(PL_INFO,"PSUADE_BEGIN\n");
   printOutTS(PL_INFO,"<p> <nOutputs> <m> <design parameter identifiers>\n");
   printOutTS(PL_INFO,"1 <design values...> <O1 mean> <O1 std dev> ... <On std dev> \n");
   printOutTS(PL_INFO,"2 <design values...> <O1 mean> <O1 std dev> ... <On std dev> \n");
   printOutTS(PL_INFO,"...\n");
   printOutTS(PL_INFO,"p <design values...> <O1 mean> <O1 std dev> ... <On std dev> \n");
   printOutTS(PL_INFO,"PSUADE_END\n");
   printDashes(PL_INFO, 0);
   printOutTS(PL_INFO,"The likelihood function to be constructed is in the form of:\n");
   printOutTS(PL_INFO,"      C exp(-0.5*S) \n");
   printOutTS(PL_INFO,"where C is some internal normalization constant and\n");
   printOutTS(PL_INFO,"      S = 1/ p sum_{k=1}^p sum_{i=1)^n (Y_ki - m_ki)^2/sd_ki^2\n");
   printOutTS(PL_INFO,"where n is the number of outputs and m_ki and\n");
   printOutTS(PL_INFO,"      sd_ki are the mean and standard deviation of output i of\n");
   printOutTS(PL_INFO,"      group k.\n");
   printDashes(PL_INFO, 0);
   printOutTS(PL_INFO,"NOTE: Alternately, your simulator (or response surface) output\n");
   printOutTS(PL_INFO,"      may instead be some error measure from comparison of all\n");
   printOutTS(PL_INFO,"      model outputs with observational data. In this case, set\n");
   printOutTS(PL_INFO,"      nOutputs=1, mean = 0, and standard deviation = 1 in this\n");
   printOutTS(PL_INFO,"      specification file (that is, your simulation output is S\n");
   printOutTS(PL_INFO,"      above and MCMC will generate likelihood C exp(-0.5 S).\n");
   }

   printf("===> Enter the spec file for building the likelihood function : ");
   cin >> sfname;
   fgets(lineIn, 1000, stdin); 
   kk = sfname.size();
   if (kk <= 1000)
   {
      cfname[kk] = '\0';
      sfname.copy(cfname, kk, 0);
      ifile.open(cfname);
      if (! ifile.is_open())
      {
         printOutTS(PL_ERROR,"MCMC ERROR : cannot open spec file %s.\n", cfname);
         return PSUADE_UNDEFINED;
      }
   }
   else
   {
      printOutTS(PL_ERROR, "MCMC ERROR: file name too long.\n");
      return PSUADE_UNDEFINED;
   }
   getline (ifile, iline);
   compFlag = iline.compare(0,12,"PSUADE_BEGIN");
   if (compFlag == 0)
   {
      ifile >> dnSamples;
      if (dnSamples <= 0)
      {
         printOutTS(PL_ERROR,"MCMC ERROR: number of experimental data points <= 0.\n");
         ifile.close();
         return PSUADE_UNDEFINED;
      }
      printOutTS(PL_INFO,"SPEC FILE: Number of observation set    = %d\n", dnSamples);

      ifile >> kk;
      if (kk != nOutputs)
      {
         printOutTS(PL_ERROR,"MCMC ERROR: number of outputs for each experiment does\n");
         printOutTS(PL_ERROR,"     not match the simulator nOutputs.\n");
         printOutTS(PL_ERROR,"     %d versus %d\n", kk, nOutputs);
         ifile.close();
         return PSUADE_UNDEFINED;
      }
      printOutTS(PL_INFO,"SPEC FILE: Number of simulation/emulation outputs = %d\n",nOutputs);

      ifile >> dnInputs;
      if (dnInputs < 0)
      {
         printOutTS(PL_ERROR,"MCMC ERROR: number of experimental design variables < 0.\n");
         ifile.close();
         return PSUADE_UNDEFINED;
      }
      if (dnInputs > nInputs)
      {
         printOutTS(PL_ERROR,"MCMC ERROR: number of experimental design variables %d\n",
                    dnInputs);
         printOutTS(PL_ERROR,"     cannot be larger than number of inputs to the\n");
         printOutTS(PL_ERROR,"     simulator %d.\n", nInputs);
         ifile.close();
         return PSUADE_UNDEFINED;
      }
      printOutTS(PL_INFO,"SPEC FILE: Number of design parameters  = %d\n", dnInputs);
      if (dnInputs > 0)
      {
         designParams = new int[nInputs];
         for (ii = 0; ii < nInputs; ii++) designParams[ii] = 0;
         cnt = 0;
         for (ii = 0; ii < dnInputs; ii++)
         {
            ifile >> kk;
            if (kk <= 0 || kk > nInputs)
            {
               printOutTS(PL_ERROR,"MCMC ERROR: invalid design parameter %d.\n",kk);
               ifile.close();
               delete [] designParams;
               return PSUADE_UNDEFINED;
            }
            if (kk <= cnt)
            {
               printOutTS(PL_ERROR,"MCMC ERROR: design parameters should be in ascending order.\n");
               ifile.close();
               delete [] designParams;
               return PSUADE_UNDEFINED;
            }
            designParams[kk-1] = 1;
            printOutTS(PL_INFO,"SPEC FILE: input %d is a design parameter\n", kk);
            cnt = kk;
         }
         dSamInputs = new double[dnSamples*dnInputs];
      }
      dSamMeans = new double[dnSamples*nOutputs];
      dSamStdevs = new double[dnSamples*nOutputs];
      for (ii = 0; ii < dnSamples; ii++)
      {
         ifile >> kk;
         if (kk != ii+1)
         {
            printOutTS(PL_ERROR,"MCMC ERROR: invalid index %d at line %d in spec file.\n",
                   kk, ii+2);
            printOutTS(PL_ERROR,"            (Expecting %d).\n", ii+1);
            printOutTS(PL_ERROR,"==> check line %d\n", ii+3);
            ifile.close();
            if (dSamInputs != NULL) delete [] dSamInputs;
            delete [] dSamMeans;
            delete [] dSamStdevs;
            return PSUADE_UNDEFINED;
         }
         if (printLevel > 0)
            printOutTS(PL_INFO,"Calibration Data Set %d\n", kk);
         for (jj = 0; jj < dnInputs; jj++)
         {
            ifile >> dSamInputs[ii*dnInputs+jj];
            if (printLevel > 0)
               printOutTS(PL_INFO,"   Design parameter %d = %e\n", jj+1, 
                          dSamInputs[ii*dnInputs+jj]);
         }
         for (jj = 0; jj < nOutputs; jj++)
         {
            ifile >> dSamMeans[ii*nOutputs+jj];
            ifile >> dSamStdevs[ii*nOutputs+jj];
            if (printLevel > 0)
               printOutTS(PL_INFO,"      Data mean/stdev = %16.8e %16.8e\n",
                          dSamMeans[ii*nOutputs+jj],dSamStdevs[ii*nOutputs+jj]);
            if (dSamStdevs[ii*nOutputs+jj] < 0.0)
            {
               ifile.close();
               printOutTS(PL_ERROR,"MCMC ERROR: std dev in spec file should be > 0.\n");
               printOutTS(PL_ERROR,"==> check the last entry in line %d\n", ii+3);
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
      printOutTS(PL_ERROR, 
           "MCMC ERROR: PSUADE_BEGIN missing at the beginning of the spec file.\n");
      ifile.close();
      return PSUADE_UNDEFINED;
   }
   getline (ifile, iline);
   getline (ifile, iline);
   compFlag = iline.compare(0,10,"PSUADE_END");
   if (compFlag != 0)
   {
      printOutTS(PL_ERROR,"MCMC ERROR: PSUADE_END missing at the end of the spec file.\n");
      ifile.close();
      delete [] dSamMeans;
      delete [] dSamStdevs;
      return PSUADE_UNDEFINED;
   }
   ifile.close();
   printEquals(PL_INFO, 0);

   if (psAnaExpertMode_ == 1 && mode_ == 0)
   {
      printOutTS(PL_INFO,"*** OPTION TO INCLUDE RESPONSE SURFACE UNCERTAINTIES: \n\n");
      printOutTS(PL_INFO,"To incorporate the response surface errors into the\n");
      printOutTS(PL_INFO,"likelihood function, make sure that either Gaussian\n");
      printOutTS(PL_INFO,"process/Kriging, polynomial regression, or bootstrapped\n");
      printOutTS(PL_INFO,"MARS response surface is selected in the simulation\n");
      printOutTS(PL_INFO,"data file. Otherwise, no RS uncertainties will be included.\n\n");
      printOutTS(PL_INFO,"NOTE: if you don't know what this feature is, just say no.\n");
      printf( "===> Include response surface uncertainties? (y or n) ");
      scanf("%s", charString);
      fgets(lineIn,1000,stdin);
      if (charString[0] == 'y') rsErrFlag = 1;
      printEquals(PL_INFO, 0);
   }

   if (ioPtr != NULL)
   {
      ioPtr->getParameter("ana_rstype", pPtr);
      faType = pPtr.intData_;
      ioPtr->getParameter("ana_rsindexfile", pPtr);
      rsFile = pPtr.strArray_[0];
      if (strcmp(rsFile, "NONE"))
      {
         printOutTS(PL_INFO,"A response surface index file has been specified.\n");
         fp = fopen(rsFile, "r");
         if (fp == NULL)
         {
            printOutTS(PL_ERROR,"MCMC ERROR: rs_index_file %s not found.\n",rsFile);
            return PSUADE_UNDEFINED;
         }
         else
         {
            printOutTS(PL_INFO,"INFO: rs_index_file %s found.\n",rsFile);
            fscanf(fp,"%d", &kk);
            if (kk != nInputs)
            {
               printOutTS(PL_ERROR,"MCMC ERROR: invalid nInputs in rs_index_file.\n");
               printOutTS(PL_ERROR,"  Data format should be: \n");
               printOutTS(PL_ERROR,"  line 1: nInputs in rs data (driver) file\n");
               printOutTS(PL_ERROR,"  line 2: 1 <1 or 0> <default value if first number==0>\n");
               printOutTS(PL_ERROR,"  line 3: 2 <2 or 0> <0 if first number != 0>\n");
               printOutTS(PL_ERROR,"  line 4: 3 <3 or 0> <default value if first number==0>\n");
               printOutTS(PL_ERROR,"  line 5: 4 <4 or 0> <0 if first number != 0>\n");
               printOutTS(PL_ERROR,"  ...\n");
               fclose(fp);
               return PSUADE_UNDEFINED;
            }
            rsIndices = new int[nInputs];
            rsValues = new double[nInputs];
            for (ii = 0; ii < nInputs; ii++) rsIndices[ii] = 0;
            for (ii = 0; ii < nInputs; ii++)
            {
               fscanf(fp, "%d", &kk);
               if (kk != ii+1)
               {
                  printOutTS(PL_ERROR,"MCMC ERROR: 1st index in indexFile = %d (must be %d]).\n",
                             kk, ii+1);
                  printOutTS(PL_ERROR,"  Data format should be: \n");
                  printOutTS(PL_ERROR,"  line 1: nInputs in rs data (driver) file\n");
                  printOutTS(PL_ERROR,"  line 2: 1 <1 or 0> <default value if first number==0>\n");
                  printOutTS(PL_ERROR,"  line 3: 2 <2 or 0> <0 if first number != 0>\n");
                  printOutTS(PL_ERROR,"  line 4: 3 <3 or 0> <default value if first number==0>\n");
                  printOutTS(PL_ERROR,"  line 5: 4 <4 or 0> <0 if first number != 0>\n");
                  printOutTS(PL_ERROR,"  ...\n");
                  fclose(fp);
                  return PSUADE_UNDEFINED;
               }
               fscanf(fp, "%d", &rsIndices[ii]);
               if (rsIndices[ii] == 0)
                  printOutTS(PL_INFO,"MCMC INFO: input %3d inactive\n",ii+1);

               if (rsIndices[ii] == 0 && designParams != NULL && designParams[ii] == 1)
               {
                  printOutTS(PL_ERROR,"MCMC ERROR: inactive input %d cannot be design parameter.\n",
                             ii+1);
                  fclose(fp);
                  return PSUADE_UNDEFINED;
               }

               if (rsIndices[ii] < 0 || rsIndices[ii] > nInputs)
               {
                  printOutTS(PL_ERROR, "MCMC INFO: input %3d = %d invalid\n",ii+1,rsIndices[ii]);
                  fclose(fp);
                  return PSUADE_UNDEFINED;
               }
               rsIndices[ii]--;
               fscanf(fp, "%lg", &rsValues[ii]);
            }
            fclose(fp);
            printOutTS(PL_INFO, "Response surface index information: \n");
            for (ii = 0; ii < nInputs; ii++)
               printOutTS(PL_INFO, "Input %4d: index = %4d, default = %e\n",
                      ii+1, rsIndices[ii]+1, rsValues[ii]);
         }
      }
   }
   else if (mode_ == 0)
   {
      printOutTS(PL_INFO,"MCMC INFO: since ioPtr=NULL, assume MARS as reponse surface.\n");
      faType = PSUADE_RS_MARS;
   }
   if (nSamples > 0)
   {
      printOutTS(PL_INFO,"MCMC INFO: CREATING RESPONSE SURFACES FOR ALL OUTPUTS.\n");
      faPtrs = new FuncApprox*[nOutputs];
      YY = new double[nSamples];
      for (ii = 0; ii < nOutputs; ii++)
      {
         faType = -1;
         printOutTS(PL_INFO,"MCMC INFO: CREATING RESPONSE SURFACE FOR OUTPUT %d.\n",ii+1);
         faPtrs[ii] = genFA(faType, nInputs, iZero, nSamples);
         faPtrs[ii]->setNPtsPerDim(16);
         faPtrs[ii]->setBounds(lower, upper);
         faPtrs[ii]->setOutputLevel(0);
         for (kk = 0; kk < nSamples; kk++) YY[kk] = Y[kk*nOutputs+ii];

         status = faPtrs[ii]->initialize(X, YY);
         if (status != 0)
         {
            printOutTS(PL_ERROR,"MCMC ERROR: Unable to create response surface.\n");
            printOutTS(PL_ERROR,"            Consult PSUADE developers.\n");
            for (kk = 0; kk < ii; kk++) delete [] faPtrs[kk];
            delete [] faPtrs;
            delete [] YY;
            return PSUADE_UNDEFINED;
         }
      }
      delete [] YY;
   }
   else if (mode_ == 0)
   {
      printOutTS(PL_ERROR,"MCMC ERROR: no sample given for generating likelihood function.\n");
      return PSUADE_UNDEFINED;
   }

   //    ==> burnInSamples, maxSamples, nbins, plotIndices, nPlots
   maxSamples = 10000;
   burnInSamples = maxSamples / 2;
   nbins = 20;
   printEquals(PL_INFO, 0);
   printOutTS(PL_INFO,"*** CURRENT SETTINGS OF MCMC PARAMETERS: \n\n");
   printOutTS(PL_INFO,"MCMC Burn-in sample size      (default) = %d\n", burnInSamples);
   printOutTS(PL_INFO,"MCMC sample increment         (default) = %d\n", maxSamples);
   printOutTS(PL_INFO,"MCMC no. of bins in histogram (default) = %d\n", nbins);
   printOutTS(PL_INFO,"NOTE: sample increment - sample size to run before convergence check\n");
   printOutTS(PL_INFO, "NOTE: histogram nBins  - define granularity of histogram bar graph\n");
   printOutTS(PL_INFO, "Turn on ana_expert mode to change these default settings.\n\n");
   if (psAnaExpertMode_ == 1)
   {
      sprintf(charString,"Enter sample increment (5000 - 50000): ");
      maxSamples = getInt(5000, 200000, charString);
      burnInSamples = maxSamples / 2;
      sprintf(charString,"Enter the number of histogram bins (10 - 25) : ");
      nbins = getInt(10, 50, charString);
   }
   if (psAnaExpertMode_ == 1)
   {
      printEquals(PL_INFO, 0);
      printOutTS(PL_INFO,"*** OPTION TO CREATE POSTERIOR PLOTS FOR SELECTED INPUTS:\n\n");
      printOutTS(PL_INFO,"MCMC will create MATLAB files for the posterior distributions.\n");
      printOutTS(PL_INFO,"You can choose to generate posterior plots for all inputs, or \n");
      printOutTS(PL_INFO,"just a selected few (in case there are too many inputs).\n");
      printf("Select inputs for which posterior plots are to be generated.\n");
      sprintf(charString,"Enter input number (-1 for all, 0 to terminate) : ");
      kk = 1;
      plotIndices = new int[nInputs];
      nPlots = 0;
      while (kk != 0 || nPlots < 1)
      {
         kk = getInt(-1, nInputs, charString);
         if (kk == -1)
         {
            for (ii = 0; ii < nInputs; ii++)
            {
               if (rsIndices == NULL || rsIndices[ii] >= 0)
                  if (designParams == NULL || designParams[ii] == 0) 
                     plotIndices[nPlots++] = ii;
            }
            break;
         }
         if (kk != 0)
         {
            if (rsIndices != NULL && rsIndices[kk-1] < 0) 
               printOutTS(PL_ERROR,
                          "Input %d has been fixed by the rs index file (no plot).\n",
                          kk+1);
            else if (designParams != NULL && designParams[kk-1] == 1)
               printOutTS(PL_ERROR,"Input %d is a design parameter (no plot).\n",kk);
            else 
               plotIndices[nPlots++] = kk - 1;
         }
         if (kk == 0 && nPlots == 0)
            printOutTS(PL_ERROR,
                 "You need to set at least 1 input for plotting posteriors.\n");
      }
      if (nPlots > 1) sortIntList(nPlots, plotIndices);
   }
   else
   {
      plotIndices = new int[nInputs];
      nPlots = 0;
      for (ii = 0; ii < nInputs; ii++) 
         if (rsIndices == NULL || rsIndices[ii] >= 0)
            if (designParams == NULL || designParams[ii] == 0) 
               plotIndices[nPlots++] = ii;
   }
   printOutTS(PL_INFO, 
              "MCMC Plot summary: input number to be plotted are (%d):\n", nPlots);
   for (ii = 0; ii < nPlots; ii++)
      printOutTS(PL_INFO, "   Input %4d\n", plotIndices[ii]+1);

   // option to add discrepancy function and a posterior sample
   // ==> modelFormFlag, genPosteriors
   if (psAnaExpertMode_ == 1)
   {
      printEquals(PL_INFO, 0);
      printOutTS(PL_INFO,"*** OPTION TO ADD A DISCREPANCY FUNCTION:\n\n");
      printOutTS(PL_INFO,"To use this feature, first make sure that the observation\n");
      printOutTS(PL_INFO,"data file specified earlier has design parameters specified\n");
      printOutTS(PL_INFO,"since the discrepancy function is to be a function of these\n");
      printOutTS(PL_INFO,"design parameters (if not, a constant discrepancy function\n");
      printOutTS(PL_INFO,"is to be created).\n");
      printOutTS(PL_INFO,"NOTE: if you don't know what this feature is, just say NO.\n");
      printf("===> Add discrepancy function ? (y or n) ");
      scanf("%s", charString);
      fgets(lineIn,1000,stdin);
      if (charString[0] == 'y') modelFormFlag = 1;

      printEquals(PL_INFO, 0);
      printOutTS(PL_INFO,"*** OPTION TO CREATE A SAMPLE FROM THE POSTERIOR DISTRIBUTIONS:\n\n");
      printOutTS(PL_INFO,"In addition to generating the posterior distributions, you can\n");
      printOutTS(PL_INFO,"also draw a sample from these posteriors. The posterior sample\n");
      printOutTS(PL_INFO,"can be used as prior sample for another simulator/emulator.\n");
      printOutTS(PL_INFO,"NOTE: if you don't know what this feature is, just say no.\n");
      printf("===> Generate posterior samples for the calibration parameters? (y or n) ");
      scanf("%s", charString);
      fgets(lineIn,1000,stdin);
      if (charString[0] == 'y') genPosteriors = 1;
   } 
   printEquals(PL_INFO, 0);

   // setup input PDF, if there is any
   if (printLevel > 2) 
      printOutTS(PL_INFO,"*** INFORMATION ON PARAMETER PRIOR DISTRIBUTIONS\n");
   inputPDFs = new PDFBase*[nInputs];
   for (ii = 0; ii < nInputs; ii++)
   {
      if (pdfFlags != NULL && pdfFlags[ii] == PSUADE_PDF_NORMAL)
      {
         inputPDFs[ii] = (PDFBase *) new PDFNormal(inputMeans[ii], inputStdevs[ii]);
         if (printLevel > 2) 
            printOutTS(PL_INFO,"Parameter %3d has normal prior distribution (%e,%e).\n",
                   ii+1, inputMeans[ii], inputStdevs[ii]);
      }
      else if (pdfFlags != NULL && pdfFlags[ii] == PSUADE_PDF_LOGNORMAL)
      {
         inputPDFs[ii] = (PDFBase *) new PDFLogNormal(inputMeans[ii],inputStdevs[ii]);
         if (printLevel > 2)
            printOutTS(PL_INFO,"Parameter %3d has lognormal prior distribution.\n",ii+1);
      }
      else if (pdfFlags != NULL && pdfFlags[ii] == PSUADE_PDF_TRIANGLE)
      {
         inputPDFs[ii] = (PDFBase *) new PDFTriangle(inputMeans[ii],inputStdevs[ii]);
         if (printLevel > 2)
            printOutTS(PL_INFO,"Parameter %3d has triangle prior distribution.\n",ii+1);
      }
      else if (pdfFlags != NULL && pdfFlags[ii] == PSUADE_PDF_BETA)
      {
         inputPDFs[ii] = (PDFBase *) new PDFBeta(inputMeans[ii],inputStdevs[ii]);
         if (printLevel > 2)
            printOutTS(PL_INFO,"Parameter %3d has beta prior distribution.\n",ii+1);
      }
      else if (pdfFlags != NULL && pdfFlags[ii] == PSUADE_PDF_WEIBULL)
      {
         inputPDFs[ii] = (PDFBase *) new PDFWeibull(inputMeans[ii],inputStdevs[ii]);
         if (printLevel > 2)
            printOutTS(PL_INFO,"Parameter %3d has Weibull prior distribution.\n",ii+1);
      }
      else if (pdfFlags != NULL && pdfFlags[ii] == PSUADE_PDF_GAMMA)
      {
         inputPDFs[ii] = (PDFBase *) new PDFGamma(inputMeans[ii],inputStdevs[ii]);
         if (printLevel > 2)
            printOutTS(PL_INFO,"Parameter %3d has gamma prior distribution.\n",ii+1);
      }
      else if (pdfFlags != NULL && pdfFlags[ii] == 1000+PSUADE_PDF_NORMAL)
      {
         inputPDFs[ii] = NULL;
         if (printLevel > 2)
         {
            printOutTS(PL_INFO,"Parameter %3d: multi-parameter normal distribution.\n",ii+1);
            printOutTS(PL_INFO,"               curently not supported.\n");
            return -1.0;
         }
      }
      else if (pdfFlags != NULL && pdfFlags[ii] == 1000+PSUADE_PDF_LOGNORMAL)
      {
         inputPDFs[ii] = NULL;
         if (printLevel > 2)
         {
            printOutTS(PL_INFO,"Parameter %3d: multi-parameter lognormal distribution.\n",ii+1);
            printOutTS(PL_INFO,"               curently not supported.\n");
            return -1.0;
         }
      }
      else if (pdfFlags == NULL || pdfFlags[ii] == PSUADE_PDF_UNIFORM)
      {
         inputPDFs[ii] = NULL;
         if (printLevel > 2)
            printOutTS(PL_INFO,"Parameter %3d has uniform prior distribution.\n",ii+1);
      }
      else if (pdfFlags != NULL && pdfFlags[ii] == PSUADE_PDF_USER)
      {
         inputPDFs[ii] = NULL;
         if (printLevel > 2)
         {
            printOutTS(PL_INFO,"Parameter %3d: user-provided distribution currently not\n",ii+1);
            printOutTS(PL_INFO,"               supported.\n");
            return -1.0;
         }
      }
   }
   if (printLevel > 2) printEquals(PL_INFO, 0);

   //    ==> funcIO, freq (how often to do simulation and emulation)
   maxPts = nbins * 5;
   if (psAnaExpertMode_ == 1)
   {
      printOutTS(PL_INFO,"*** SETTING PROPOSAL DISTRIBUTION RESOLUTION\n");
      printOutTS(PL_INFO,"Since MCMC uses many function evaluations to construct\n");
      printOutTS(PL_INFO,"the proposal distributions, you have the option to set\n");
      printOutTS(PL_INFO,"how many points are used to construct it in order to\n");
      printOutTS(PL_INFO,"keep the inference cost reasonable.\n");
      printf("Sample size to construct proposal distributions. Default is %d.\n",maxPts);
      sprintf(charString,"Enter new sample size (%d - %d): ",nbins*3,nbins*10);
      maxPts = getInt(nbins*3, nbins*10, charString);
      maxPts = maxPts / nbins * nbins;
      printOutTS(PL_INFO,"Proposal distribution sample size = %d.\n", maxPts);
   }
   if (mode_ == 1 && ioPtr != NULL) 
   {
      funcIO = createFunctionInterface(ioPtr);
      if (nSamples == 0)
         printOutTS(PL_INFO, "MCMC: DIRECT SIMULATION has been set up.\n");
      else
         printOutTS(PL_INFO,"MCMC: DIRECT SIMULATION PLUS RESPONSE SURFACE have been set up.\n");
      printOutTS(PL_INFO,"MCMC INFO: make sure simulation output is in the right form\n");
      printOutTS(PL_INFO,"     which should not have been translated (mean) nor scaled (std)\n");
      printOutTS(PL_INFO,"     unless you compute the error measure yourself, in which case\n");
      printOutTS(PL_INFO,"     you should have use nOutputs=1, mean=0, and std dev=1 in the\n");
      printOutTS(PL_INFO,"     spec file.\n");
      if (psAnaExpertMode_ == 1 && nSamples > 0)
      {
         printOutTS(PL_INFO,"Since MCMC uses many function evaluations, you have the option\n");
         printOutTS(PL_INFO,"to set how frequent the simulator is invoked (response surface\n");
         printOutTS(PL_INFO,"will be used otherwise) in constructing the proposal distribution\n");
         printOutTS(PL_INFO,"during each MCMC iteration. A frequency of f means that each MCMC\n");
         printOutTS(PL_INFO,"step uses f simulator runs. These f simulator runs will then be\n");
         printOutTS(PL_INFO,"supplemented with evaluations from the given response surface.\n");
         printOutTS(PL_INFO,"The default is f=10 (if you do not know what this is, enter 10).\n");
         sprintf(charString,"Max. number of simulator runs per MCMC step (1 - %d, default=10)? ",
                 maxPts);
         kk = getInt(1, maxPts, charString);
         freq = maxPts * nInputs / kk;
         if (freq * kk != maxPts) freq++;
      }
      else if (nSamples > 0) freq = maxPts / 10;
      else                   freq = 1;
      if (nSamples > 0)
         printOutTS(PL_INFO,"Frequency of invoking the simulator has been set to %d\n", freq);
   }
   if (funcIO == NULL && faPtrs == NULL)
   {
      printOutTS(PL_ERROR,"MCMC ERROR: missing simulator and sample data - cannot proceed.\n");
      return PSUADE_UNDEFINED;
   }

   if (modelFormFlag == 1)
   {
      int    *ExpSamStates, ind, ExpNSamples, dfaType, dnPerDim=16;
      double *dOneSample, expdata, simdata;
      double *ExpSamOutputs, *ExpSamInputs, *tSamInputs, *settings;

      ExpNSamples   = dnSamples;
      ExpSamInputs  = new double[ExpNSamples*nInputs];
      ExpSamOutputs = new double[ExpNSamples];
      ExpSamStates  = new int[ExpNSamples];
      discFuncConstantMeans = new double[nOutputs];
      discFuncConstantStds  = new double[nOutputs];
      for (ii2 = 0; ii2 < nOutputs; ii2++)
         discFuncConstantMeans[ii2] = discFuncConstantStds[ii2] = PSUADE_UNDEFINED;

      printOutTS(PL_INFO,"*** SELECT RESPONSE SURFACE TYPE FOR DISCREPANCY FUNCTION:\n");
      dfaType = -1;
      while (dfaType < 0 || dfaType >= PSUADE_NUM_RS)
      {
         writeFAInfo(-1);
         sprintf(charString, "===> Enter your choice : ");
         dfaType = getInt(0, PSUADE_NUM_RS-1, charString);
      }

      settings = new double[nInputs];
      for (ii2 = 0; ii2 < nInputs; ii2++)
      {
         if (inputPDFs == NULL || (inputPDFs != NULL && inputPDFs[ii2] == NULL))
            settings[ii2] = 0.5*(lower[ii2] + upper[ii2]);
         else
            settings[ii2] = inputPDFs[ii2]->getMean();
         if (rsIndices != NULL && rsIndices[ii2] < 0)
            settings[ii2] = rsValues[ii2];
      }

      faPtrs1 = new FuncApprox*[nOutputs];
      dOneSample = new double[nInputs];
      tSamInputs = new double[ExpNSamples*nInputs];
      int        askFlag = 0, *states=NULL;
      double     *tLowers = new double[nInputs];
      double     *tUppers = new double[nInputs];
      PsuadeData *dataPtr = new PsuadeData(); 
      char       **iNames;
      for (ii = 0; ii < nOutputs; ii++)
      {
         for (kk = 0; kk < dnSamples; kk++)
         {
            cnt = 0;
            for (ii2 = 0; ii2 < nInputs; ii2++)
            {
               if (designParams != NULL && designParams[ii2] == 1)
               {
                  dOneSample[ii2] = dSamInputs[kk*dnInputs+cnt];
                  cnt++;
               }
               else dOneSample[ii2] = settings[ii2];
            }

            simdata = 0.0;
            if (psAnaExpertMode_ == 1 && askFlag == 0)
            {
               printOutTS(PL_INFO,"To create discrepancy functions, the calibration\n");
               printOutTS(PL_INFO,"parameters need to be set to some nominal values.\n");
               printOutTS(PL_INFO,"You can choose the nominal values, or it will be\n");
               printOutTS(PL_INFO,"set to the input means or mid points of the ranges.\n");
               printf( "Set nomininal values yourself ? (y or n) ");
               scanf("%s", charString);
               fgets(lineIn,1000,stdin);
               if (charString[0] == 'y')
               {
                  for (ii2 = 0; ii2 < nInputs; ii2++)
                  {
                     if ((rsIndices == NULL || rsIndices[ii2] >= 0) &&
                         (designParams == NULL || designParams[ii2] == 0))
                     {
                        printOutTS(PL_INFO, "Input %d has lower and upper bounds = %e %e\n",
                                   ii2+1, lower[ii2], upper[ii2]);
                        sprintf(charString, "Nominal value for input %d : ",ii2+1);
                        dOneSample[ii2] = getDouble(charString);
                        settings[ii2]   = dOneSample[ii2];
                     }
                  }
               }
               askFlag = 1;
            }

            if (funcIO != NULL)
                 funcIO->evaluate(kk+1,nInputs,dOneSample,1,&simdata,0);
            else simdata = faPtrs[ii]->evaluatePoint(dOneSample);
            expdata = dSamMeans[kk*nOutputs+ii];

            if (printLevel >= 4)
            {
               printOutTS(PL_INFO, "Experiment %4d (out of %d) : ",kk+1,dnSamples);
               for (ii2 = 0; ii2 < nInputs; ii2++)
                     printOutTS(PL_INFO,"Input %7d = %12.4e ",ii2+1,dOneSample[ii2]);
               printOutTS(PL_INFO, "simuation, experimental data = %12.4e %12.4e\n",
                          simdata, expdata);
            }

            ExpSamOutputs[kk] = expdata - simdata;
         }

         if (dnInputs > 0) 
         {
            for (kk = 0; kk < ExpNSamples*dnInputs; kk++)
               tSamInputs[kk] = dSamInputs[kk];
            iNames = new char*[dnInputs];
            cnt = 0;
            for (ii2 = 0; ii2 < nInputs; ii2++)
            {
               if (designParams[ii2] == 1)
               {
                  iNames[cnt] = new char[100];
                  tLowers[cnt] = lower[ii2];
                  tUppers[cnt] = upper[ii2];
                  if (qData.strArray_ == NULL)
                       sprintf(iNames[cnt], "X%d", ii2+1);
                  else strcpy(iNames[cnt], qData.strArray_[ii2]);
                  cnt++;
               }
            }
            dataPtr->updateInputSection(ExpNSamples, dnInputs, NULL, tLowers, tUppers,
                                        dSamInputs, iNames, NULL,NULL,NULL,NULL);
            for (ii2 = 0; ii2 < dnInputs; ii2++) delete [] iNames[ii2];
            delete [] iNames;
         }
         else
         {
            iNames = new char*[1];
            iNames[0] = new char[100];
            sprintf(iNames[0], "X0");
            for (ii2 = 0; ii2 < ExpNSamples; ii2++) tSamInputs[ii2] = 0.5;
            tLowers[0] = 0.0;
            tUppers[0] = 1.0;
            dataPtr->updateInputSection(ExpNSamples, iOne, NULL, tLowers, tUppers,
                                        tSamInputs, iNames, NULL,NULL,NULL,NULL);
            delete [] iNames[0];
            delete [] iNames;
         }

         states = new int[ExpNSamples];
         for (kk = 0; kk < ExpNSamples; kk++) states[kk] = 1;
         iNames = new char*[1];
         iNames[0] = new char[100];
         sprintf(iNames[0], "Y%d", ii+1);
         dataPtr->updateOutputSection(ExpNSamples,iOne,ExpSamOutputs,states,iNames);
         delete [] states;
         delete [] iNames[0];
         delete [] iNames;
         dataPtr->updateMethodSection(PSUADE_SAMP_MC, ExpNSamples, 1, -1, -1);
         sprintf(charString, "psDiscrepancyModel%d", ii+1);
         dataPtr->writePsuadeFile(charString, 0);
         delete dataPtr;

         printOutTS(PL_INFO,"Creating discrepancy response surface for output %d\n",ii+1);
         faPtrs1[ii] = NULL;
         if (dnInputs > 0 && dnSamples > 1)
         {
            faPtrs1[ii] = genFA(dfaType,dnInputs,iOne,ExpNSamples);
            if (faPtrs1[ii] == NULL)
            {
               printOutTS(PL_ERROR,
                    "MCMC ERROR: cannot create discrepancy func for output %d.\n",ii+1);
               return -1.0;
            }
         }
         if (faPtrs1[ii] != NULL)
         {
            faPtrs1[ii]->setNPtsPerDim(dnPerDim);
            faPtrs1[ii]->setBounds(lower, upper);
            faPtrs1[ii]->setOutputLevel(0);
            faPtrs1[ii]->initialize(tSamInputs,ExpSamOutputs);
         }
         else
         {
            discFuncConstantMeans[ii] = 0.0;
            for (kk = 0; kk < ExpNSamples; kk++)
               discFuncConstantMeans[ii] += ExpSamOutputs[kk];
            discFuncConstantMeans[ii] /= (double) ExpNSamples;
            discFuncConstantStds[ii] = 0.0;
            for (kk = 0; kk < ExpNSamples; kk++)
               discFuncConstantStds[ii] += pow(ExpSamOutputs[kk]-discFuncConstantMeans[ii],2.0);
            discFuncConstantStds[ii] = sqrt(discFuncConstantStds[ii]/ExpNSamples);
         }
      }
      delete [] ExpSamInputs;
      delete [] dOneSample;
      delete [] ExpSamOutputs;
      delete [] settings;
      delete [] tSamInputs;
      delete [] tLowers;
      delete [] tUppers;
   }

   //    set up constraint filters, if any
   printEquals(PL_INFO, 0);
   printOutTS(PL_INFO,"MCMC INFO: creating constraints, if there is any.\n");
   printOutTS(PL_INFO,"     Constraints remove infeasible regions from the priors.\n");
   printOutTS(PL_INFO,"     Constraints can be specified by RS constraint files.\n");
   constrPtr = new RSConstraints();
   constrPtr->genConstraints(ioPtr);
   printEquals(PL_INFO, 0);
   if (psAnaExpertMode_ == 1)
   {
      sprintf(charString, "How many MCMC chains? (2-20, default=3) : ");
      numChains = getInt(2,20,charString);
      sprintf(charString, "PSRF threshold? (1.0 - 1.2, default = 1.05) : ");
      psrfThreshold = getDouble(charString);
      if (psrfThreshold < 1.0 || psrfThreshold > 1.2)
      {
         printf("MCMC : invalid PSRF threshold ==> reset to 1.05.\n");
         psrfThreshold = 1.05;
      }
   }

   //    set up for MCMC iterations
   int    *Ivec, **bins, ****bins2, globalIts, countTrack, dcnt;
   int    mcmcFail=0, sumBins, index2, nFail;
   int    ii3, jj2, kk2, index, length, count, iChain, chainCnt, mcmcIts;
   int    maxGlobalIts=20, chainCntSave, *chainStatus;
   double *XRange=NULL, *XGuess=NULL, *XDist=NULL, *XDesignS, *YDesignS;
   double *YDesignStds=NULL, *XGuessS=NULL, *YGuessS=NULL, *YGuessStds=NULL;
   double Xtemp, Ytemp, Ytemp2, *Xmax, Ymax, *s2Vec, *SDist;
   double ***XChains=NULL, stdev, stdev2, ddata, ddata2, WStat, BStat;
   double *chainMeans=NULL, *chainStdevs=NULL, *psrfs=NULL;
   TwoSampleAnalyzer *s2Analyzer=NULL;
   XRange  = new double[nInputs];
   for (ii = 0; ii < nInputs; ii++) XRange[ii] = upper[ii] - lower[ii]; 
   XDist   = new double[maxPts+1];
   SDist   = new double[maxPts+1];
   XGuess  = new double[nInputs];
   XGuessS = new double[dnSamples*nInputs*(maxPts+1)];
   YGuessS = new double[dnSamples*nOutputs*(maxPts+1)];
   YGuessStds = new double[dnSamples*nOutputs*(maxPts+1)];
   XDesignS = new double[dnSamples*nInputs*(maxPts+1)];
   YDesignS = new double[dnSamples*nOutputs*(maxPts+1)];
   YDesignStds = new double[dnSamples*nOutputs*(maxPts+1)];
   means_ = new double[nInputs_];
   sigmas_ = new double[nInputs_];
   for (ii = 0; ii < nInputs; ii++) means_[ii] = sigmas_[ii] = 0.0;
   Xmax = new double[nInputs];
   for (ii = 0; ii < nInputs; ii++) Xmax[ii] = 0;
   mostLikelyInput_ = new double[nInputs_];
   for (ii = 0; ii < nInputs_; ii++) mostLikelyInput_[ii] = 0;
   mostLikelyOutput_ = new double[nOutputs_];
   for (ii = 0; ii < nOutputs_; ii++) mostLikelyOutput_[ii] = 0;
   Ymax = -PSUADE_UNDEFINED;
   Ivec = new int[nInputs];
   Ivec[nInputs-1] = -1;
   XChains = new double**[numChains];
   for (ii = 0; ii < numChains; ii++)
   {
      XChains[ii] = new double*[maxGlobalIts*maxSamples];
      for (jj = 0; jj < maxGlobalIts*maxSamples; jj++)
         XChains[ii][jj] = new double[nInputs+1];
   }
   chainMeans = new double[numChains];
   chainStdevs = new double[numChains];
   chainStatus  = new int[numChains];
   assert(chainStatus);
   for (ii = 0; ii < numChains; ii++) chainMeans[ii] = chainStdevs[ii] = 0.0;
   for (ii = 0; ii < numChains; ii++) chainStatus[ii] = 0;
   psrfs = new double[nInputs];
   for (ii = 0; ii < nInputs; ii++) psrfs[ii] = 0.0;
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
   s2Vec = new double[maxGlobalIts*maxSamples];
   if (printLevel > 3) s2Analyzer = new TwoSampleAnalyzer();

   Sampling *sampler;
   if (nInputs > 50) sampler = SamplingCreateFromID(PSUADE_SAMP_LHS);
   else              sampler = SamplingCreateFromID(PSUADE_SAMP_LPTAU);
   sampler->setInputBounds(nInputs, lower, upper);
   sampler->setOutputParams(1);
   sampler->setSamplingParams(numChains, 1, 1);
   sampler->initialize(0);
   double *mcmcSeeds = new double[numChains*nInputs];
   double *tmpOuts = new double[numChains];
   int    *tmpStates = new int[numChains];
   sampler->getSamples(numChains,nInputs,1,mcmcSeeds,tmpOuts,tmpStates);
   delete [] tmpOuts;
   delete [] tmpStates;
   delete sampler;
   for (iChain = 0; iChain < numChains; iChain++)
   {
      for (ii = 0; ii < nInputs; ii++)
      {
         ddata = mcmcSeeds[iChain*nInputs+ii];
         ddata = (ddata - lower[ii]) / XRange[ii];
         mcmcSeeds[iChain*nInputs+ii] = ddata;
      }
   }
   
   // run the Gibbs algorithm
   printAsterisks(PL_INFO, 0);
   printOutTS(PL_INFO, "MCMC begins ... \n");
   fflush(stdout);
   fp = NULL;
   globalIts = chainCnt = 0;
   while (globalIts < maxGlobalIts)
   {
      for (iChain = 0; iChain < numChains; iChain++)
      {
         printOutTS(PL_INFO,"MCMC : chain %d, iteration = %d\n",iChain+1,
                    globalIts+1);
         if (iChain == 0) chainCntSave = chainCnt;
         else             chainCnt     = chainCntSave;
         if (chainCnt == 0)
         {
            for (ii = 0; ii < nInputs; ii++)
            {
               if (designParams == NULL || designParams[ii] == 0)
                    XGuess[ii] = mcmcSeeds[iChain*nInputs+ii];
               else XGuess[ii] = 0.5;
               XChains[iChain][0][ii] = XGuess[ii];
               XGuess[ii] = XGuess[ii] * (upper[ii] - lower[ii]) + lower[ii];
            }
            if (rsIndices != NULL)
            {
               for (ii = 0; ii < nInputs; ii++)
                  if (rsIndices[ii] < 0) XGuess[ii] = rsValues[ii];
            }
            XChains[iChain][0][nInputs] = -1;
         }
         else
         {
            for (ii = 0; ii < nInputs; ii++)
            {
               ddata = XChains[iChain][chainCnt-1][ii];
               XGuess[ii] = ddata * (upper[ii] - lower[ii]) + lower[ii];
            }
         }
         
         mcmcIts = countTrack = 0;
         while (mcmcIts < maxSamples)
         {
            count = (mcmcIts+1) / (maxSamples/10);
            if (count != countTrack)
            { 
               countTrack++;
               printOutTS(PL_INFO, "%3.0f%% ",10.0*countTrack );
               fflush(stdout);
            }
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
                   (rsIndices != NULL && rsIndices[ii] >=0)) && 
                   (designParams == NULL || designParams[ii] == 0))
               {
                  Xtemp = XGuess[ii];
                  if (psGMMode_ == 1)
                  {
                     fp = fopen("MCMCDistTrack.m", "w");
                     if (fp == NULL)
                     {
                        printOutTS(PL_WARN,"MCMC ERROR: cannot write to MCMCDistTrack.m file.\n");
                     }
                     else
                     {
                        fprintf(fp,"%% This file contains individual terms in the exponent\n");
                        fprintf(fp,"%% of the likelihood function, i.e. in each row \n");
                        fprintf(fp,"%% first column : \n");
                        fprintf(fp,"%% S = 1/(p*n) sum_{k=1}^p sum_{i=1)^n (Y_ki - m_ki)^2/sd_ki^2\n");
                        fprintf(fp,"%% next columns: log-likelihood terms\n");
                        fprintf(fp,"nOuts = %d;\n", nOutputs);
                        fprintf(fp,"nObs  = %d;\n", dnSamples);
                        fprintf(fp,"nPts  = %d;\n", maxPts+1);
                        fprintf(fp,"A = [\n");
                     }
                  }
                  if (nOutputs == 1)
                  {
                     cnt = 0;
                     for (jj = 0; jj <= maxPts; jj++)
                     {
                        XGuess[ii] = lower[ii]+jj*XRange[ii]/maxPts;
                     
                        cnt++;
                        if (cnt >= freq) cnt = 0; 
                        index = jj * dnSamples;
                        for (kk2 = 0; kk2 < dnSamples; kk2++) 
                        {
                           dcnt = 0;
                           for (ii2 = 0; ii2 < nInputs; ii2++) 
                           {
                              XGuessS[(index+kk2)*nInputs+ii2] = XGuess[ii2];
                              if (designParams != NULL && designParams[ii2] == 1)
                              {
                                 XGuessS[(index+kk2)*nInputs+ii2] = dSamInputs[kk2*dnInputs+dcnt];
                                 XDesignS[(index+kk2)*dnInputs+dcnt] = dSamInputs[kk2*dnInputs+dcnt];
                                 if (XGuessS[(index+kk2)*nInputs+ii2] < lower[ii2] || 
                                     XGuessS[(index+kk2)*nInputs+ii2] > upper[ii2])
                                 {
                                    printOutTS(PL_INFO,"WARNING: design parameter value out of bound.\n");
                                    printOutTS(PL_INFO,"   Input = %d\n", ii2+1);
                                    printOutTS(PL_INFO,"   Bounds = [%24.16e, %24.16e].\n", 
                                               lower[ii2],upper[ii2]);
                                    printOutTS(PL_INFO,"   Design value = %e\n", 
                                               XGuessS[(index+kk2)*nInputs+ii2]);
                                    printOutTS(PL_INFO,"   Design sample = %d\n", kk2+1);
                                 }
                                 dcnt++;
                              }
                           }
                        }
                        for (kk2 = 0; kk2 < dnSamples; kk2++) 
                        {
                           YDesignS[index+kk2] = YDesignStds[index+kk2] = 0.0;
                           YGuessS[index+kk2] = YGuessStds[index+kk2] = 0.0;
                        }
                     }

                     if (rsErrFlag == 1)
                     {
                        if (funcIO == NULL || cnt > 0)
                           faPtrs[0]->evaluatePointFuzzy((maxPts+1)*dnSamples,
                                                         XGuessS,YGuessS,YGuessStds);
                        else
                        {
                           for (kk2 = 0; kk2 < (maxPts+1)*dnSamples; kk2++) 
                           {
                              funcIO->evaluate(kk2+1,nInputs,&XGuessS[kk2*nInputs],1,
                                               &YGuessS[kk2],0);
                              YGuessStds[kk2] = 0.0;
                           }
                        }
                        if (faPtrs1 != NULL && faPtrs1[0] != NULL)
                        {
                           faPtrs1[0]->evaluatePointFuzzy((maxPts+1)*dnSamples,XDesignS,
                                                          YDesignS,YDesignStds);
                        }
                        else if (discFuncConstantMeans != NULL &&
                                 discFuncConstantMeans[0] != PSUADE_UNDEFINED)
                        {
                           for (kk2 = 0; kk2 < (maxPts+1)*dnSamples; kk2++) 
                              YDesignS[kk2] = discFuncConstantMeans[0];
                        }
                     }
                     else
                     {
                        if (funcIO == NULL || cnt > 0)
                           faPtrs[0]->evaluatePoint((maxPts+1)*dnSamples,XGuessS,YGuessS);
                        else
                        {
                           for (kk2 = 0; kk2 < (maxPts+1)*dnSamples; kk2++) 
                              funcIO->evaluate(kk2+1,nInputs,&XGuessS[kk2*nInputs],1,&YGuessS[kk2],0);
                        }
                        if (faPtrs1 != NULL && faPtrs1[0] != NULL)
                        {
                           faPtrs1[0]->evaluatePoint((maxPts+1)*dnSamples,XDesignS,YDesignS);
                        }
                        else if (discFuncConstantMeans != NULL &&
                                 discFuncConstantMeans[0] != PSUADE_UNDEFINED)
                        {
                           for (kk2 = 0; kk2 < (maxPts+1)*dnSamples; kk2++) 
                              YDesignS[kk2] = discFuncConstantMeans[0];
                        }
                     }

                     for (jj = 0; jj <= maxPts; jj++)
                     {
                        index = jj * dnSamples;
                        for (kk2 = 0; kk2 < dnSamples; kk2++) YGuessS[index+kk2] += YDesignS[index+kk2];

                        XDist[jj] = 0.0;
                        for (kk2 = 0; kk2 < dnSamples; kk2++)
                        {
                           Ytemp = YGuessS[index+kk2];
                           stdev = YGuessStds[index+kk2];
                           stdev2 = YDesignStds[index+kk2];
                           Ytemp2 = pow((Ytemp-dSamMeans[kk2]),2.0) /
                                  (pow(dSamStdevs[kk2],2.0)+stdev*stdev+stdev2*stdev2);
                           XDist[jj] += Ytemp2;
                        }
                        XDist[jj] = XDist[jj] / dnSamples;
                        if (psGMMode_ == 1 && fp != NULL) 
                        {
                           fprintf(fp, "%e ", XDist[jj]);
                           for (kk2 = 0; kk2 < dnSamples; kk2++)
                           {
                              Ytemp = YGuessS[index+kk2];
                              stdev = YGuessStds[index+kk2];
                              stdev2 = YDesignStds[index+kk2];
                              Ytemp2 = pow((Ytemp-dSamMeans[kk2]),2.0) /
                                     (pow(dSamStdevs[kk2],2.0)+stdev*stdev+stdev2*stdev2);
                              fprintf(fp, "%e ", Ytemp2);
                           }
                           fprintf(fp, "\n");
                        }
                     }
                  }
                  else
                  {
                     for (jj = 0; jj <= maxPts; jj++)
                     {
                        XGuess[ii] = lower[ii]+jj*XRange[ii]/maxPts;
                        index = jj * dnSamples;
                        for (kk2 = 0; kk2 < dnSamples; kk2++) 
                        {
                           dcnt = 0;
                           for (ii2 = 0; ii2 < nInputs; ii2++) 
                           {
                              XGuessS[(index+kk2)*nInputs+ii2] = XGuess[ii2];
                              if (designParams != NULL && designParams[ii2] == 1)
                              {
                                 XGuessS[(index+kk2)*nInputs+ii2] = dSamInputs[kk2*dnInputs+dcnt];
                                 XDesignS[(index+kk2)*dnInputs+dcnt] = dSamInputs[kk2*dnInputs+dcnt];
                                 if (XGuessS[(index+kk2)*nInputs+ii2] < lower[ii2] || 
                                     XGuessS[(index+kk2)*nInputs+ii2] > upper[ii2])
                                 {
                                    printOutTS(PL_INFO,"WARNING: design parameter value out of bound.\n");
                                    printOutTS(PL_INFO,"   Input = %d\n", ii2+1);
                                    printOutTS(PL_INFO,"   Bounds = [%24.16e, %24.16e].\n", 
                                               lower[ii2],upper[ii2]);
                                    printOutTS(PL_INFO,"   Design value = %e\n", 
                                               XGuessS[(index+kk2)*nInputs+ii2]);
                                    printOutTS(PL_INFO,"   Design sample = %d\n", kk2+1);
                                 }
                                 dcnt++;
                              }
                           }
                        }
                        for (ii2 = 0; ii2 < dnSamples*nOutputs; ii2++) 
                        {
                           YDesignS[index*nOutputs+ii2] = YDesignStds[index*nOutputs+ii2] = 0.0;
                           YGuessS[index*nOutputs+ii2] = YGuessStds[index*nOutputs+ii2] = 0.0;
                        }
                     }
                     for (ii2 = 0; ii2 < nOutputs; ii2++) 
                     {
                        if (rsErrFlag == 1)
                        {
                           faPtrs[ii2]->evaluatePointFuzzy((maxPts+1)*dnSamples,XGuessS,
                                           &YGuessS[ii2*dnSamples*(maxPts+1)],
                                           &YGuessStds[ii2*dnSamples*(maxPts+1)]);
                           if (faPtrs1 != NULL && faPtrs1[ii2] != NULL)
                           {
                              faPtrs1[ii2]->evaluatePointFuzzy((maxPts+1)*dnSamples,XDesignS,
                                               &YDesignS[ii2*dnSamples*(maxPts+1)],
                                               &YDesignStds[ii2*dnSamples*(maxPts+1)]);
                           }
                           else if (discFuncConstantMeans != NULL &&
                                    discFuncConstantMeans[ii2] != PSUADE_UNDEFINED)
                           {
                              for (kk2 = 0; kk2 < dnSamples*(maxPts+1); kk2++)
                              {
                                 YDesignS[ii2*dnSamples*(maxPts+1)+kk2] = discFuncConstantMeans[ii2];
                                 YDesignStds[ii2*dnSamples*(maxPts+1)+kk2] = 0.0;
                              }
                           }
                        }
                        else
                        {
                           faPtrs[ii2]->evaluatePoint((maxPts+1)*dnSamples,XGuessS,
                                               &YGuessS[ii2*dnSamples*(maxPts+1)]);
                           if (faPtrs1 != NULL && faPtrs1[ii2] != NULL)
                           {
                              faPtrs1[ii2]->evaluatePoint(dnSamples*(maxPts+1),XDesignS,
                                                          &YDesignS[ii2*dnSamples*(maxPts+1)]);
                           }
                           else if (discFuncConstantMeans != NULL &&
                                    discFuncConstantMeans[ii2] != PSUADE_UNDEFINED)
                           {
                              for (kk2 = 0; kk2 < dnSamples*(maxPts+1); kk2++)
                                 YDesignS[ii2*dnSamples*(maxPts+1)+kk2] = discFuncConstantMeans[ii2];
                           }
                           for (kk2 = 0; kk2 < dnSamples*(maxPts+1); kk2++)
                              YGuessStds[ii2*dnSamples*(maxPts+1)+kk2] = 
                                    YDesignStds[ii2*dnSamples*(maxPts+1)+kk2] = 0.0;
                        }
                     }
                     for (jj = 0; jj <= maxPts; jj++)
                     {
                        XDist[jj] = 0.0;
                        index = jj * dnSamples;
                        for (ii2 = 0; ii2 < nOutputs; ii2++) 
                        {
                           for (kk2 = 0; kk2 < dnSamples; kk2++)
                           {
                              Ytemp = YGuessS[ii2*dnSamples*(maxPts+1)+index+kk2] + 
                                      YDesignS[ii2*dnSamples*(maxPts+1)+index+kk2];
                              stdev = YGuessStds[ii2*dnSamples*(maxPts+1)+index+kk2];
                              stdev2 = YDesignStds[ii2*dnSamples*(maxPts+1)+index+kk2];
                              Ytemp2 = pow((Ytemp-dSamMeans[kk2*nOutputs+ii2]),2.0) /
                                       (pow(dSamStdevs[kk2*nOutputs+ii2],2.0) + 
                                        stdev*stdev + stdev2*stdev2);
                              XDist[jj] += Ytemp2;
                           }
                        }
                        XDist[jj] = XDist[jj] / (dnSamples*nOutputs);
                        if (psGMMode_ == 1 && fp != NULL) 
                        {
                           fprintf(fp, "%e ", XDist[jj]);
                           for (ii2 = 0; ii2 < nOutputs; ii2++) 
                           {
                              for (kk2 = 0; kk2 < dnSamples; kk2++)
                              {
                                 Ytemp = YGuessS[ii2*dnSamples*(maxPts+1)+index+kk2] + 
                                         YDesignS[ii2*dnSamples*(maxPts+1)+index+kk2];
                                 stdev = YGuessStds[ii2*dnSamples*(maxPts+1)+index+kk2];
                                 stdev2 = YDesignStds[ii2*dnSamples*(maxPts+1)+index+kk2];
                                 Ytemp2 = pow((Ytemp-dSamMeans[kk2*nOutputs+ii2]),2.0) /
                                          (pow(dSamStdevs[kk2*nOutputs+ii2],2.0) + 
                                           stdev*stdev + stdev2*stdev2);
                                 fprintf(fp, "%e ", Ytemp2);
                              }
                           }
                           fprintf(fp, "\n");
                        }
                     }
                  }

                  nFail = 0;
                  for (jj = 0; jj <= maxPts; jj++)
                  {
                     XGuess[ii] = lower[ii]+jj*XRange[ii]/maxPts;
                     Ytemp = constrPtr->evaluate(XGuess,XDist[jj],status);
                     if (status == 0)
                     {
                        XDist[jj] = 0.0;
                        nFail++;
                     }

                     ddata = 1.0;
                     if (inputPDFs != NULL)
                     {
                        for (ii2 = 0; ii2 < nInputs; ii2++)
                        {
                           if ((designParams == NULL || designParams[ii2] == 0) && 
                               inputPDFs[ii2] != NULL &&
                               (rsIndices == NULL || rsIndices[ii2] >= 0))
                           {
                              inputPDFs[ii2]->getPDF(iOne,&XGuess[ii2],&ddata2);
                              ddata *= ddata2;
                           }
                        }
                     }
                     SDist[jj] = ddata; 
                  }

                  ddata = XDist[0];
                  for (jj = 1; jj <= maxPts; jj++) 
                  {
                     ddata2 = XDist[jj];
                     if (ddata2 < ddata) ddata = ddata2;
                  }
                  XChains[iChain][chainCnt][nInputs] = ddata;
                  for (jj = 0; jj <= maxPts; jj++)
                  {
                     XDist[jj] = SDist[jj] * exp(-0.5 * (XDist[jj] - ddata));

                     ddata2 = XDist[jj] * exp(-0.5*ddata);
                     if (ddata2 > Ymax)
                     {
                        Ymax = ddata2;
                        for (ii2 = 0; ii2 < nInputs; ii2++) Xmax[ii2] = XGuess[ii2];
                        Xmax[ii] = lower[ii] + (upper[ii] - lower[ii]) * jj / maxPts;
                     }
                     if (jj > 0) XDist[jj] += XDist[jj-1];
                  }
                  if (psGMMode_ == 1 && fp != NULL)
                  {
                     fprintf(fp,"];\n");
                     fprintf(fp,"S = [\n");
                     for (jj = 0; jj <= maxPts; jj++) fprintf(fp, "%e\n", SDist[jj]);
                     fprintf(fp,"];\n");
                     fprintf(fp,"P = A(:,1) .* S;\n");
                     fprintf(fp,"subplot(1,2,1)\n");
                     fprintf(fp,"X = %e : %e : %e;\n",lower[ii],(upper[ii]-lower[ii])/maxPts-1e-8,
                             upper[ii]);
                     fprintf(fp,"if length(X) > nPts\n");
                     fprintf(fp,"   X = X(1:nPts);\n");
                     fprintf(fp,"end;\n");
                     fprintf(fp,"plot(X,P,'lineWidth',2.0)\n");
                     fprintf(fp,"set(gca,'linewidth',2)\n");
                     fprintf(fp,"set(gca,'fontweight','bold')\n");
                     fprintf(fp,"set(gca,'fontsize',12)\n");
                     fprintf(fp,"ylabel('log likelihood','FontWeight','bold','FontSize',12)\n");
                     fprintf(fp,"title('Input %d','FontWeight','bold','FontSize',12)\n",ii+1);
                     fprintf(fp,"grid on\n");
                     fprintf(fp,"box on\n");
                     fprintf(fp,"subplot(1,2,2)\n");
                     fprintf(fp,"P2 = exp(-0.5*P);\n");
                     fprintf(fp,"for ii = 2 : %d\n",maxPts+1);
                     fprintf(fp,"   P2(ii) = P2(ii) + P2(ii-1);\n");
                     fprintf(fp,"end;\n");
                     fprintf(fp,"P2 = P2 / P2(nPts);\n");
                     fprintf(fp,"plot(X,P2,'lineWidth',2.0)\n");
                     fprintf(fp,"set(gca,'linewidth',2)\n");
                     fprintf(fp,"set(gca,'fontweight','bold')\n");
                     fprintf(fp,"set(gca,'fontsize',12)\n");
                     fprintf(fp,"title('Input %d','FontWeight','bold','FontSize',12)\n",ii+1);
                     fprintf(fp,"ylabel('likelihood CDF','FontWeight','bold','FontSize',12)\n");
                     fprintf(fp,"grid on\n");
                     fprintf(fp,"box on\n");
                     for (ii2 = 0; ii2 < nInputs; ii2++)
                        printf("Current guess %d = %e\n",ii2+1, XGuess[ii2]);
                  }

                  if (printLevel > 3)
                     printOutTS(PL_INFO,"proposal distribution max = %e\n", XDist[maxPts]);
                  if (XDist[maxPts] - XDist[0] > 0.0e-24)
                  {
                     for (jj = 1; jj <= maxPts; jj++)
                        XDist[jj] = (XDist[jj] - XDist[0]) / (XDist[maxPts] - XDist[0]);
                     XDist[0] = 0;
                     XGuess[ii] = PSUADE_drand();
                     index = binarySearchDble(XGuess[ii], XDist, maxPts+1);
                     //if (index == maxPts) index = maxPts - 1;
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
                     if (printLevel > 2)
                        printOutTS(PL_INFO,"MCMC iteration %7d : no modification in input %d\n",
                                   mcmcIts+1,ii+1);
                     if (nFail == maxPts+1) 
                     {
                        printOutTS(PL_ERROR,"ERROR: Constraints have resulted in zero distribution.\n");
                        exit(1);
                     }
                  }
                  if (psGMMode_ == 1 && fp != NULL)
                  {
                     fclose(fp);
                     fp = NULL; 
                     printf("INFO: Next selection for input %d = %e (%d)\n", ii+1, XGuess[ii],index);
                     printf("Now examine MCMCDistTrack.m for diagnostics.\n");
                     printf("When you are done, type 1 to continue or 0 to exit.\n");
                     scanf("%d", &ii2);
                     if (ii2 == 0) exit(1);
                  }

                  if (mcmcIts >= maxSamples/2 || globalIts > 0)
                  {
                     for (ii2 = 0; ii2 < nInputs; ii2++) 
                     {
                        ddata = (XGuess[ii2] - lower[ii2]) / XRange[ii2];
                        XChains[iChain][chainCnt][ii2] = ddata;
                        index = (int) (ddata * nbins);
                        if (index >= nbins) index = nbins - 1;
                        bins[index][ii2]++;
                     }
                     for (ii2 = 0; ii2 < nInputs; ii2++) 
                     {
                        ddata = (XGuess[ii2] - lower[ii2]) / XRange[ii2];
                        index = (int) (ddata * nbins);
                        if (index >= nbins) index = nbins - 1;
                        for (ii3 = 0; ii3 < nInputs; ii3++) 
                        {
                           ddata2 = (XGuess[ii3] - lower[ii3]) / XRange[ii3];
                           index2 = (int) (ddata2 * nbins);
                           if (index2 >= nbins) index2 = nbins - 1;
                           bins2[index][index2][ii2][ii3]++;
                        }
                     }
                     chainCnt++;
                  }
                  mcmcIts++;
               }
               fp = fopen("psuade_stop", "r");
               if (fp != NULL)
               {
                  printOutTS(PL_ERROR, "MCMC INFO: psuade_stop FILE FOUND - TERMINATE MCMC.\n");
                  fclose(fp);
                  return 0.0;
               }
               if (mcmcIts >= maxSamples) break;
            }
         }
         if (countTrack <= 10) printOutTS(PL_INFO,"100%%\n");
         else                  printOutTS(PL_INFO,"\n");
      }
 
      globalIts++;
      printOutTS(PL_INFO, "\nIteration %d summary: \n", globalIts);
      mcmcFail = nInputs - dnInputs;
      if (rsIndices != NULL)
      {
         for (ii = 0; ii < nInputs; ii++)
            if (rsIndices[ii] < 0) mcmcFail--; 
      }
      for (ii = 0; ii < nInputs; ii++)
      {
         if ((rsIndices == NULL || (rsIndices != NULL && rsIndices[ii] >=0)) && 
             (designParams == NULL || designParams[ii] == 0))
         {
            if (printLevel > 2) printOutTS(PL_INFO, "Input = %d\n", ii+1);
          
            for (iChain = 0; iChain < numChains; iChain++)
            {
               ddata = 0.0;
               for (jj = 0; jj < chainCnt; jj++) 
                  ddata += XChains[iChain][jj][ii];
               ddata /= chainCnt;
               ddata2 = 0.0;
               for (jj = 0; jj < chainCnt; jj++) 
                   ddata2 += pow(XChains[iChain][jj][ii]-ddata,2.0);
               ddata2 /= (double) (chainCnt - 1);
               chainMeans[iChain] = ddata;
               chainStdevs[iChain] = ddata2;
               if (globalIts > 2 && chainStdevs[iChain] < 1.0e-20) 
               {
                  printf("MCMC INFO: chain %d disabled.\n",iChain+1);
                  chainStatus[iChain] = 1;
               }
            }
            nChainGood = 0;
            for (iChain = 0; iChain < numChains; iChain++)
            {
               if (chainStatus[iChain] == 0) nChainGood++;
            }
            if (nChainGood <= 1)
            {
               printf("MCMC ERROR: too few chains <= 1.\n");
               printf("Suggestion: you may want to relax on the experimental data\n");
               printf("            uncertainties (make them larger).\n");
               printf("            To see if this is the problem, turn on printlevel\n");
               printf("            to 3 and run again. If the variance of the chains\n");
               printf("            are small, small data uncertainties is probably the\n");
               printf("            problem.\n");
               exit(1);
            }
            WStat = 0.0;
            for (iChain = 0; iChain < numChains; iChain++)
            {
               if (chainStatus[iChain] == 0) WStat += chainStdevs[iChain];
            }
            WStat /= (double) nChainGood;
            if (WStat < 0) WStat = PSUADE_UNDEFINED;
            if (printLevel > 2) printf("  Within  chain variance W = %e\n", WStat);
            ddata = 0.0;
            for (iChain = 0; iChain < numChains; iChain++)
            {
               if (chainStatus[iChain] == 0)
                  ddata += chainMeans[iChain];
            }
            ddata /= (double) nChainGood;
            BStat = 0.0;
            for (iChain = 0; iChain < numChains; iChain++)
            {
               if (chainStatus[iChain] == 0)
                  BStat += pow(chainMeans[iChain]-ddata,2.0);
            }
            BStat = BStat / (nChainGood - 1.0) * chainCnt;
            if (printLevel > 2) printf("  Between chain variance B = %e\n", BStat/chainCnt);
            ddata = (1 - 1.0/chainCnt) * WStat + BStat / chainCnt;
            ddata = ddata / WStat * (numChains + 1) / numChains - 
                      (chainCnt - 1.0) / (double) (chainCnt * numChains); 
            if (ddata < 0) ddata2 = PSUADE_UNDEFINED;
            else           ddata2 = sqrt(ddata);
            if (printLevel > 2)
            {
               for (iChain = 0; iChain < numChains; iChain++)
                  printf("  Chain %d mean, var = %e %e\n", iChain+1,
                         chainMeans[iChain]*XRange[ii]+lower[ii],
                         chainStdevs[iChain]*XRange[ii]*XRange[ii]);
               printf("  Chain length             = %d\n", chainCnt);
               printf("  Weighted average of B, W = %e\n", ddata);
            }
            printf("  Input %d PSRF = %e\n", ii+1, ddata2);
            psrfs[ii] = ddata2;
            if (ddata2 < psrfThreshold)
            {
               printOutTS(PL_INFO,"MCMC INFO : PSRF < %e ==> converged.\n",
                          psrfThreshold);
               mcmcFail--;
            }
            ddata = 0.0;
            for (iChain = 0; iChain < numChains; iChain++)
            {
               if (chainStatus[iChain] == 0)
                  for (jj = 0; jj < chainCnt; jj++)
                     ddata += XChains[iChain][jj][ii];
            }
            ddata /= (double) (nChainGood * chainCnt);
            means_[ii] = ddata;
            ddata2 = 0.0;
            for (iChain = 0; iChain < numChains; iChain++)
            {
               if (chainStatus[iChain] == 0)
                  for (jj = 0; jj < chainCnt; jj++)
                     ddata2 += pow(XChains[iChain][jj][ii]-ddata,2.0);
            }
            ddata2 /= (double) (chainCnt*nChainGood-1);
            sigmas_[ii] = sqrt(ddata2);
         }
      }
      if (mcmcFail == 0 && printLevel > 3 && s2Analyzer != NULL)
      {
         if (rsIndices != NULL)
         {
            for (ii = 0; ii < nInputs; ii++)
               if (rsIndices[ii] < 0) mcmcFail--; 
         }
         for (ii = 0; ii < nInputs; ii++)
         {
            if ((rsIndices == NULL || (rsIndices != NULL && rsIndices[ii] >=0)) && 
                (designParams == NULL || designParams[ii] == 0))
            {
               printOutTS(PL_INFO, "Geweke Input = %d\n", ii+1);
          
               cnt = chainCnt / 2;
               ddata2 = (double) numChains;
               for (iChain = 0; iChain < numChains; iChain++)
               {
                  for (jj = 0; jj < 2*cnt; jj++) 
                     s2Vec[jj] = XChains[iChain][chainCnt-2*cnt+jj][ii];
                  ddata = s2Analyzer->TAnalyze(cnt,s2Vec,cnt,&s2Vec[cnt],1);
               }
            }
         }
      }

      for (ii = 0; ii < nInputs; ii++) 
      {
         if ((rsIndices == NULL || (rsIndices != NULL && rsIndices[ii] >=0)) && 
             (designParams == NULL || designParams[ii] == 0))
         {
            printOutTS(PL_INFO,"MCMC: input %3d value at peak of likelihood = %e\n",
                       ii+1, Xmax[ii]);
            ddata = means_[ii]*(upper[ii]-lower[ii])+lower[ii];
            printOutTS(PL_INFO,"MCMC: input %3d mean    = %e\n", ii+1, ddata);
            ddata = sigmas_[ii]*(upper[ii]-lower[ii]);
            printOutTS(PL_INFO,"MCMC: input %3d std dev = %e\n", ii+1, ddata);
            mostLikelyInput_[ii] = Xmax[ii];
         }
      }

      if (psMasterMode_ == 1)
      {
         fp = fopen("MCMCChainHistogram.m", "w");
         if (fp != NULL)
         {
            for (ii = 0; ii < nInputs; ii++)
            {
               if ((rsIndices == NULL || rsIndices[ii] >= 0) &&
                   (designParams == NULL || designParams[ii] == 0)) 
               {
                  fprintf(fp,"X%d = [\n",ii+1);
                  for (jj = 0; jj < chainCnt; jj++)
                  {
                     for (iChain = 0; iChain < numChains; iChain++)
                        fprintf(fp,"%e ", XChains[iChain][jj][ii]);
                     fprintf(fp,"\n");
                  }
                  fprintf(fp,"];\n");
                  fprintf(fp,"[nk, nx] = hist(X%d);\n", ii+1);
                  fprintf(fp,"nk = nk / %d;\n", chainCnt);
                  fprintf(fp,"bar(nx, nk)\n");
                  fprintf(fp,"%%plot(X%d)\n", ii+1);
                  sprintf(charString, "Input %d", ii+1);
                  fwritePlotTitle(fp, charString);
                  fwritePlotXLabel(fp, "Input Value");
                  fwritePlotYLabel(fp, "Input Count");
                  fwritePlotAxes(fp);
                  fprintf(fp,"text(0.05,0.9,'colors are for different chains','sc')\n");
                  fprintf(fp,"disp('Colors in the histograms are for different chains.\n");
                  fprintf(fp,"disp('Press enter to continue to next input')\n");
                  fprintf(fp,"pause\n");
               }
            }
            fclose(fp);
            fp = NULL;
            printf("MCMC: A file called MCMCChainHistogram.m has been created.\n");
            printf("      Use Matlab to view the histograms of all chains for all\n");
            printf("      inputs to assess convergence.\n");
            sprintf(charString, "Enter 1 to continue or 0 to terminate : ");
            ii = getInt(0, 10, charString);
            if (ii == 0) break;
         }
         else
         {
            printOutTS(PL_INFO,"MCMC INFO: cannot create MCMCChainHistogram.m file.\n");
         } 
      }   
 
      if (mcmcFail == 0) break;

      genMatlabFile(nInputs,lower,upper,XRange,nPlots,plotIndices,nbins,
                    bins,bins2,qData,numChains,chainCnt,XChains,chainStatus);
 
      fp = fopen("psuade_stop", "r");
      if (fp != NULL)
      {
         printOutTS(PL_INFO,"MCMC INFO: psuade_stop FILE FOUND - TERMINATE MCMC.\n");
         fclose(fp);
         fp = NULL;
         strcpy(charString, "psuade_stop");
         unlink(charString);
         break;
      }
      fp = fopen("psuade_nogm", "r");
      if (fp != NULL)
      {
         printOutTS(PL_INFO,"MCMC INFO: psuade_nogm FILE FOUND. GM mode is now off.\n");
         fclose(fp);
         fp = NULL;
         psGMMode_ = 0;
         strcpy(charString, "psuade_nogm");
         unlink(charString);
      }
      fp = fopen("psuade_gm", "r");
      if (fp != NULL)
      {
         printOutTS(PL_INFO,"MCMC INFO: psuade_gm FILE FOUND. GM mode is now on.\n");
         fclose(fp);
         fp = NULL;
         psGMMode_ = 1;
         strcpy(charString, "psuade_gm");
         unlink(charString);
      }
      fp = fopen("psuade_master", "r");
      if (fp != NULL)
      {
         printOutTS(PL_INFO,"MCMC INFO: psuade_master FILE FOUND. Master mode is now on.\n");
         fclose(fp);
         fp = NULL;
         psMasterMode_ = 1;
         strcpy(charString, "psuade_master");
         unlink(charString);
      }
      fp = fopen("psuade_nomaster", "r");
      if (fp != NULL)
      {
         printOutTS(PL_INFO,"MCMC INFO: psuade_nomaster FILE FOUND. Master mode is now off.\n");
         fclose(fp);
         fp = NULL;
         psMasterMode_ = 0;
         strcpy(charString, "psuade_nomaster");
         unlink(charString);
      }
      fp = fopen("psuade_print", "r");
      if (fp != NULL)
      {
         printOutTS(PL_INFO,"MCMC INFO: psuade_print FILE FOUND. Print level is set to 3.\n");
         fclose(fp);
         fp = NULL;
         printLevel = 3;
         strcpy(charString, "psuade_print");
         unlink(charString);
      }
   }
   if (globalIts >= maxGlobalIts)
   {
      mcmcFail = 0;
      for (ii = 0; ii < nInputs; ii++) 
      {
         if ((rsIndices == NULL || (rsIndices != NULL && rsIndices[ii] >=0)) && 
             (designParams == NULL || designParams[ii] == 0))
            if (psrfs[ii] > psrfThreshold) mcmcFail = 1;
      }
      if (mcmcFail == 1) 
         printOutTS(PL_INFO,"MCMC maximum iterations exceeded but no convergence.\n");
   }
   else printOutTS(PL_INFO, "MCMC iterations completed\n");

   for (ii = 0; ii < nInputs; ii++) 
      for (jj = 0; jj < nbins; jj++) bins[jj][ii] = 0;
   for (ii = 0; ii < nInputs; ii++) 
      for (ii2 = 0; ii2 < nInputs; ii2++) 
         for (jj = 0; jj < nbins; jj++)
            for (jj2 = 0; jj2 < nbins; jj2++)
               bins2[jj][jj2][ii][ii2] = 0;
   for (iChain = 0; iChain < numChains; iChain++) 
   { 
      if (chainStatus[iChain] == 0)
      {
         for (jj = 0; jj < chainCnt; jj++) 
         {
            for (ii2 = 0; ii2 < nInputs; ii2++) 
            {
               ddata = XChains[iChain][jj][ii2];
               index = (int) (ddata * nbins);
               if (index >= nbins) index = nbins - 1;
               bins[index][ii2]++;
            }
            for (ii2 = 0; ii2 < nInputs; ii2++) 
            {
               ddata = XChains[iChain][jj][ii2];
               index = (int) (ddata * nbins);
               if (index >= nbins) index = nbins - 1;
               for (ii3 = 0; ii3 < nInputs; ii3++) 
               {
                  ddata2 = XChains[iChain][jj][ii3];
                  index2 = (int) (ddata2 * nbins);
                  if (index2 >= nbins) index2 = nbins - 1;
                  bins2[index][index2][ii2][ii3]++;
               }
            }
         }
      }
   }
  
   genMatlabFile(nInputs,lower,upper,XRange,nPlots,plotIndices,nbins,
                 bins,bins2,qData,numChains,chainCnt,XChains,chainStatus);

   if (genPosteriors == 1)
   {
      fp = fopen("MCMCPostSample", "w");
      if (fp != NULL)
      {
         fprintf(fp, "PSUADE_BEGIN\n");
         cnt = nChainGood * chainCnt;
         if (cnt > 200000) cnt = 200000;
         cnt /= nChainGood;
         fprintf(fp, "%d %d\n", cnt*nChainGood,nInputs);
         if (qData.strArray_ != NULL)
         {
            fprintf(fp, "# ");
            for (jj = 0; jj < nInputs; jj++)
               fprintf(fp,"%s ", qData.strArray_[jj]);
            fprintf(fp, "\n");
         }
         ii2 = 0;
         for (iChain = 0; iChain < numChains; iChain++)
         { 
            if (chainStatus[iChain] == 0)
            {
               for (ii = chainCnt-cnt; ii < chainCnt; ii++)
               {
                  fprintf(fp, "%d ", ii2+1);
                  for (jj = 0; jj < nInputs; jj++)
                  {
                     if ((rsIndices == NULL || rsIndices[jj] >= 0) &&
                         (designParams == NULL || designParams[jj] == 0)) 
                     {
                        ddata = XChains[iChain][ii][jj] * XRange[jj] + lower[jj];
                        fprintf(fp, "%e ", ddata);
                     }
                     else if (rsIndices != NULL && rsIndices[jj] == 0)
                        fprintf(fp, "%e ", rsValues[jj]);
                     else if (designParams != NULL && designParams[jj] != 0) 
                        fprintf(fp, "%e ", 0.5 * (upper[jj] + lower[jj]));
                  }
                  fprintf(fp, "\n");
                  ii2++;
               }
            }
         }
         fprintf(fp, "PSUADE_END\n");
         fprintf(fp, "#N=%d;\n",nChainGood*cnt);
         fprintf(fp, "#m=%d;\n",cnt);
         for (iChain = 0; iChain < numChains; iChain++)
         {
            if (chainStatus[iChain] == 0)
               fprintf(fp, "#A%d = A(%d*m+1:%d*m,:);\n",iChain,iChain,iChain+1);
         }
         fprintf(fp, "#for ii = 2 : %d\n", nInputs+1);
         for (iChain = 0; iChain < numChains; iChain++)
         {
            if (chainStatus[iChain] == 0)
            {
               fprintf(fp, "#subplot(*,*,%d)\n",iChain+1);
               fprintf(fp, "#hist(A%d(:,ii))\n",iChain+1);
            }
            fprintf(fp, "#ii-1\n");
            fprintf(fp, "#pause;\n");
         }
         fprintf(fp, "#end;\n");
         fclose(fp);
      }
      printOutTS(PL_INFO,"MCMC: 'MCMCPostSample' file has a posterior sample.\n");
   }

   int    nInps, nOuts, nSams, *states;
   double *allOuts;
   char   **oNames;
   PsuadeData *filePtr1, *filePtr2;
   if (modelFormFlag == 1)
   {
      sprintf(charString, "psDiscrepancyModel1");
      filePtr1 = new PsuadeData();
      status = filePtr1->readPsuadeFile(charString);
      if (status < 0) exit(1);
      if (status > 0)
      {
         printf("MCMC ERROR: cannot read file %s in PSUADE format.\n",charString);
         exit(1);
      } 
   }
   if (modelFormFlag == 1 && status == 0)
   {
      filePtr1->getParameter("input_ninputs", pPtr);
      nInps = pPtr.intData_;
      filePtr1->getParameter("output_noutputs", pPtr);
      nOuts = pPtr.intData_;
      filePtr1->getParameter("method_nsamples", pPtr);
      nSams = pPtr.intData_;
      filePtr1->getParameter("output_sample", pOutputs);
      unlink(charString);
      allOuts = new double[nOutputs * nSams];
      for (jj = 0; jj < nSams; jj++) 
         allOuts[jj*nOutputs] = pOutputs.dbleArray_[jj];
      pOutputs.clean();
      oNames = new char*[nOutputs];
      oNames[0] = new char[100];
      sprintf(oNames[0], "Y1");
      for (ii = 1; ii < nOutputs; ii++)
      {
         filePtr2 = new PsuadeData();
         sprintf(charString, "psDiscrepancyModel%d", ii+1);
         status = filePtr2->readPsuadeFile(charString);
         if (status != 0) break;
         filePtr2->getParameter("input_ninputs", pPtr);
         if (pPtr.intData_ != nInps) break;
         filePtr2->getParameter("output_noutputs", pPtr);
         if (pPtr.intData_ != nOuts) break;
         filePtr2->getParameter("method_nsamples", pPtr);
         if (pPtr.intData_ != nSams) break;
         filePtr2->getParameter("output_sample", pOutputs);
         delete filePtr2;
         unlink(charString);
         for (jj = 0; jj < nSams; jj++) 
            allOuts[jj*nOutputs+ii] = pOutputs.dbleArray_[jj];
         pOutputs.clean();
         oNames[ii] = new char[100];
         sprintf(oNames[ii], "Y%d", ii+1);
      }
      if (nOutputs == 1)
      { 
         sprintf(charString, "psDiscrepancyModel");
         filePtr1->writePsuadeFile(charString, 0);
      }
      else if (ii == nOutputs)
      {
         states = new int[nSams];
         for (jj = 0; jj < nSams; jj++) states[jj] = 1;
         filePtr1->updateOutputSection(nSams,nOuts,allOuts,states,oNames);
         sprintf(charString, "psDiscrepancyModel");
         filePtr2->writePsuadeFile(charString, 0);
         printOutTS(PL_INFO,"MCMC INFO: a sample (inputs/outputs) the discrepancy model\n");
         printOutTS(PL_INFO,"           is now in psDiscrepancyModel.\n");
         delete [] states;
         for (ii = 1; ii < nOutputs; ii++) delete [] oNames[ii];
      }
      else
      {
         printOutTS(PL_INFO,"MCMC INFO: unsuccessful creation of discrepancy sample file\n");
      }
      delete [] oNames[0];
      delete [] oNames;
      delete filePtr1;
      delete [] allOuts;
   }
  
   // clean up
   for (ii = 0; ii < numChains; ii++)
   {
      for (jj = 0; jj < maxGlobalIts*maxSamples; jj++)
         delete [] XChains[ii][jj];
      delete [] XChains[ii];
   }
   delete [] XChains;
   if (inputPDFs != NULL)
   {
      for (ii = 0; ii < nInputs; ii++)
         if (inputPDFs[ii] != NULL) delete inputPDFs[ii];
      delete [] inputPDFs;
   }
   delete [] Xmax;
   delete [] XDist;
   delete [] SDist;
   delete [] XGuess;
   delete [] XGuessS;
   delete [] YGuessS;
   delete [] YGuessStds;
   delete [] XDesignS;
   delete [] YDesignS;
   delete [] YDesignStds;
   delete [] XRange;
   delete [] psrfs;
   delete [] s2Vec;
   delete [] mcmcSeeds;
   if (s2Analyzer != NULL) delete s2Analyzer;
   if (discFuncConstantMeans != NULL) delete [] discFuncConstantMeans;
   if (discFuncConstantStds  != NULL) delete [] discFuncConstantStds;
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
   delete constrPtr;
   delete [] Ivec;
   return 0.0;
}

// ************************************************************************
// write to matlab file 
// ------------------------------------------------------------------------
double MCMCAnalyzer::genMatlabFile(int nInputs, double *lower, double *upper,
                                   double *XRange, int nPlots, int *plotIndices,
                                   int nbins, int **bins, int ****bins2, 
                                   pData &qData, int nChains, int chainCnt,
                                   double ***XChains, int *chainStatus)
{
   int    kk, kk2, ii2, jj, jj2, sumBins;
   double ddata;
   char   cfname[1001], charString[1001];;
   FILE   *fp;

   if (psPlotTool_ == 1) strcpy(cfname, "scilabmcmc2.sci");
   else                  strcpy(cfname, "matlabmcmc2.m");
   fp = fopen(cfname, "w");
   if (fp == NULL)
   {
      printOutTS(PL_ERROR, "ERROR: cannot open %s file.\n", cfname);
      return 0;
   }
   sprintf(charString,"This file shows posteriors plots");
   fwriteComment(fp, charString);
   sprintf(charString,"ns  - set to 1 for 1-step smoothing of 2D contours");
   fwriteComment(fp, charString);
   sprintf(charString,"ns1 - set to 1 for 1-step smoothing of 1D histgrams");
   fwriteComment(fp, charString);
   fprintf(fp, "ns  = 0;\n");
   fprintf(fp, "ns1 = 0;\n");
   fwritePlotCLF(fp);
   fprintf(fp, "active = [\n");
   for (kk = 0; kk < nInputs; kk++)
   {
      ii2 = binarySearchInt(kk, plotIndices, nPlots);
      if (ii2 < 0) fprintf(fp, "0\n");
      else         fprintf(fp, "1\n");
   }
   fprintf(fp, "];\n");
   fprintf(fp, "L = [\n");
   for (kk = 0; kk < nInputs; kk++) fprintf(fp, "%e ",lower[kk]);
   fprintf(fp, "];\n");
   fprintf(fp, "U = [\n");
   for (kk = 0; kk < nInputs; kk++) fprintf(fp, "%e ",upper[kk]);
   fprintf(fp, "];\n");
   fprintf(fp, "iStr = {\n");
   for (kk = 0; kk < nInputs-1; kk++)
   {
      if (qData.strArray_ != NULL)
           fprintf(fp, "'%s',", qData.strArray_[kk]);
      else fprintf(fp, "'Input %d',", kk+1);
   }
   if (qData.strArray_ != NULL)
        fprintf(fp, "'%s'};\n", qData.strArray_[nInputs-1]);
   else fprintf(fp, "'Input %d'};\n", nInputs);
   fprintf(fp, "X = zeros(%d,%d);\n", nInputs, nbins);
   fprintf(fp, "D = zeros(%d,%d);\n", nInputs, nbins);
   fprintf(fp, "NC = zeros(%d,%d,%d,%d);\n",nInputs,nInputs,nbins,nbins);
   for (kk = 0; kk < nInputs; kk++)
   {
      for (kk2 = 0; kk2 < nInputs; kk2++)
      {
         if (kk == kk2)
         {
            fprintf(fp, "X(%d,:) = [\n", kk+1);
            for (jj = 0; jj < nbins; jj++)
               fprintf(fp, "%e ", XRange[kk]/nbins*(jj+0.5)+lower[kk]);
            fprintf(fp, "];\n");
            fprintf(fp, "D(%d,:) = [\n", kk+1);
            sumBins = 0;
            for (jj = 0; jj < nbins; jj++) sumBins += bins[jj][kk];
            if (sumBins == 0) sumBins = 1;
            for (jj = 0; jj < nbins; jj++)
               fprintf(fp, "%e ", (double) bins[jj][kk]/(double) sumBins);
            fprintf(fp, "];\n");
         }
         else
         {
            fprintf(fp, "NC(%d,%d,:,:) = [\n", kk+1, kk2+1);
            for (jj = 0; jj < nbins; jj++)
            {
               for (jj2 = 0; jj2 < nbins; jj2++)
                  fprintf(fp, "%d ", bins2[jj][jj2][kk][kk2]);
               fprintf(fp, "\n");
            }
            fprintf(fp, "]';\n");
         }
      }
   }
   fprintf(fp, "nInps  = length(active);\n");
   fprintf(fp, "nPlots = 0;\n");
   fprintf(fp, "for ii = 1 : nInps\n");
   fprintf(fp, "   if (active(ii) == 1)\n");
   fprintf(fp, "      nPlots = nPlots + 1;\n");
   fprintf(fp, "      active(ii) = nPlots;\n");
   fprintf(fp, "   end;\n");
   fprintf(fp, "end;\n");
   fprintf(fp, "for ii = 1 : nInps\n");
   fprintf(fp, "  for jj = ii : nInps\n");
   fprintf(fp, "    if (active(ii) ~= 0 & active(jj) ~= 0)\n");
   fprintf(fp, "      index = (active(ii)-1) * nPlots + active(jj);\n");
   fprintf(fp, "      subplot(nPlots,nPlots,index)\n");
   fprintf(fp, "      if (ii == jj)\n");
   fprintf(fp, "        n = length(D(ii,:));\n");
   fprintf(fp, "        DN = D(ii,:);\n");
   fprintf(fp, "        for kk = 1 : ns1\n");
   fprintf(fp, "          DN1 = DN;\n");
   fprintf(fp, "          for ll = 2 : n-1\n");
   fprintf(fp, "            DN(ll) = DN(ll) + DN1(ll+1);\n");
   fprintf(fp, "            DN(ll) = DN(ll) + DN1(ll-1);\n");
   fprintf(fp, "            DN(ll) = DN(ll) / 3;\n");
   fprintf(fp, "          end;\n");
   fprintf(fp, "        end;\n");
   fprintf(fp, "        bar(X(ii,:), DN, 1.0);\n");
   fprintf(fp, "        xmin = min(X(ii,:));\n");
   fprintf(fp, "        xmax = max(X(ii,:));\n");
   fprintf(fp, "        xwid = xmax - xmin;\n");
   fprintf(fp, "        xmin = xmin - 0.5 * xwid / %d;\n", nbins);
   fprintf(fp, "        xmax = xmax + 0.5 * xwid / %d;\n", nbins);
   fprintf(fp, "        ymax = max(DN);\n");
   if (psPlotTool_ == 1)
   {
      fprintf(fp, "        e = gce();\n");
      fprintf(fp, "        e.children.thickness = 2;\n");
      fprintf(fp, "        e.children.foreground = 0;\n");
      fprintf(fp, "        e.children.background = 2;\n");
      fprintf(fp, "        a = gca();\n");
      fprintf(fp, "        a.data_bounds=[xmin,0;xmax,ymax];\n");
      fprintf(fp, "        a.x_label.text = iStr(ii);\n");
      fprintf(fp, "        a.x_label.font_size = 3;\n");
      fprintf(fp, "        a.x_label.font_style = 4;\n");
      fprintf(fp, "        a.grid = [1 1];\n");
      fprintf(fp, "        a.y_label.text = Str(jj);\n");
      fprintf(fp, "        a.y_label.font_size = 3;\n");
      fprintf(fp, "        a.y_label.font_style = 4;\n");
      fprintf(fp, "        a.thickness = 2;\n");
      fprintf(fp, "        a.font_size = 3;\n");
      fprintf(fp, "        a.font_style = 4;\n");
      fprintf(fp, "        a.box = \"on\";\n");
   }
   else
   {
      fprintf(fp, "        axis([xmin xmax 0 ymax])\n");
      fprintf(fp, "        set(gca,'linewidth',2)\n");
      fprintf(fp, "        set(gca,'fontweight','bold')\n");
      fprintf(fp, "        set(gca,'fontsize',12)\n");
      fprintf(fp, "        xlabel(iStr(ii),'FontWeight','bold','FontSize',12)\n");
      fprintf(fp, "        ylabel('Probabilities','FontWeight','bold','FontSize',12)\n");
      fprintf(fp, "        grid on\n");
      fprintf(fp, "        box on\n");
   }
   fprintf(fp, "      else\n");
   fprintf(fp, "        n = length(X(jj,:));\n");
   fprintf(fp, "        XT = X(jj,:);\n");
   fprintf(fp, "        YT = X(ii,:);\n");
   fprintf(fp, "        HX = (XT(n) - XT(1)) / (n-1);\n");
   fprintf(fp, "        HY = (YT(n) - YT(1)) / (n-1);\n");
   fprintf(fp, "        ZZ = squeeze(NC(ii,jj,:,:));\n");
   fprintf(fp, "        for kk = 1 : ns\n");
   fprintf(fp, "          ZZ1 = ZZ;\n");
   fprintf(fp, "          for ll = 2 : n-1\n");
   fprintf(fp, "            for mm = 2 : n-1\n");
   fprintf(fp, "              ZZ(ll,mm) = ZZ(ll,mm) + ZZ1(ll+1,mm);\n");
   fprintf(fp, "              ZZ(ll,mm) = ZZ(ll,mm) + ZZ1(ll-1,mm);\n");
   fprintf(fp, "              ZZ(ll,mm) = ZZ(ll,mm) + ZZ1(ll,mm+1);\n");
   fprintf(fp, "              ZZ(ll,mm) = ZZ(ll,mm) + ZZ1(ll,mm-1);\n");
   fprintf(fp, "              ZZ(ll,mm) = ZZ(ll,mm) + ZZ1(ll+1,mm+1);\n");
   fprintf(fp, "              ZZ(ll,mm) = ZZ(ll,mm) + ZZ1(ll-1,mm-1);\n");
   fprintf(fp, "              ZZ(ll,mm) = ZZ(ll,mm) + ZZ1(ll-1,mm+1);\n");
   fprintf(fp, "              ZZ(ll,mm) = ZZ(ll,mm) + ZZ1(ll+1,mm-1);\n");
   fprintf(fp, "              ZZ(ll,mm) = ZZ(ll,mm) / 9;\n");
   fprintf(fp, "            end;\n");
   fprintf(fp, "          end;\n");
   fprintf(fp, "        end;\n");
   fprintf(fp, "        ZZ = ZZ / (sum(sum(ZZ)));\n");
   if (psPlotTool_ == 1)
   {
      fprintf(fp, "        XX = [XT(1):HX:XT(n)];\n");
      fprintf(fp, "        YY = [YT(1):HY:YT(n)];\n");
      fprintf(fp, "        DD = splin2d(XX,YY,ZZ);\n");
      fprintf(fp, "        HX = 0.01 * (XT(n) - XT(1));\n");
      fprintf(fp, "        HY = 0.01 * (YT(n) - YT(1));\n");
      fprintf(fp, "        X2 = [XT(1):HX:XT(n)];\n");
      fprintf(fp, "        Y2 = [YT(1):HY:YT(n)];\n");
      fprintf(fp, "        [XI, YI] = ndgrid(X2, Y2);\n");
      fprintf(fp, "        disp('interpolation')\n");
      fprintf(fp, "        ZI =interp2d(XI, YI, XX, YY, DD, \"natural\");\n");
      fprintf(fp, "        disp('interpolation done')\n");
      fprintf(fp, "        ZB = ZI;\n");
      fprintf(fp, "        nX = length(X2);\n");
      fprintf(fp, "        nY = length(Y2);\n");
      fprintf(fp, "        for ii = 1 : nX\n");
      fprintf(fp, "          for jj = 1 : nY\n");
      fprintf(fp, "            ZI(ii,jj) = ZB(ii,nY-jj+1);\n");
      fprintf(fp, "          end;\n");
      fprintf(fp, "        end;\n");
      fprintf(fp, "        zmax = max(max(ZI));\n");
      fprintf(fp, "        zmin = min(min(ZI)) / zmax;\n");
      fprintf(fp, "        ZI   = ZI / zmax;\n");
      fprintf(fp, "        zmax = 1;\n");
      fprintf(fp, "        Matplot1((ZI-zmin)/(zmax-zmin)*64,[L(ii2),L(ii),U(ii2),U(ii)]);\n");
      fprintf(fp, "        xset(\"colormap\",jetcolormap(64));\n");
      fprintf(fp, "        colorbar(zmin,zmax);\n");
      fprintf(fp, "        contour2d(X2,Y2,ZB,5,rect=[L(ii2),L(ii),U(ii2),U(ii)]);\n");
      fprintf(fp, "        xset(\"fpf\",\" \");\n");
      fprintf(fp, "        a = gca();\n");
      fprintf(fp, "        a.x_label.text = iStr(jj);\n");
      fprintf(fp, "        a.x_label.font_size = 3;\n");
      fprintf(fp, "        a.x_label.font_style = 4;\n");
      fprintf(fp, "        a.y_label.text = iStr(ii);\n");
      fprintf(fp, "        a.y_label.font_size = 3;\n");
      fprintf(fp, "        a.y_label.font_style = 4;\n");
      fwritePlotAxesNoGrid(fp);
   }
   else
   {
      fprintf(fp,"%%      [YY,XX]=meshgrid(XT(1):HX:XT(n),YT(1):HY:YT(n));\n");
      fprintf(fp,"%%      HX = 0.01 * (XT(n) - XT(1));\n");
      fprintf(fp,"%%      HY = 0.01 * (YT(n) - YT(1));\n");
      fprintf(fp,"%%      [YI,XI]=meshgrid(XT(1):HX:XT(n),YT(1):HY:YT(n));\n");
      fprintf(fp,"%%      ZI=interp2(YY, XX, ZZ, YI, XI, 'spline');\n");
      fprintf(fp,"%%      pcolor(XI,YI,ZI)\n");
      fprintf(fp,"%%      shading interp\n");
      fprintf(fp,"%%      hold on\n");
      fprintf(fp,"%%      contour(XI,YI,ZI,5,'k')\n");
      fprintf(fp,"        imagesc(ZZ')\n");
      fprintf(fp,"        xtick = L(ii):(U(ii)-L(ii))/4:U(ii);\n");
      fprintf(fp,"        set(gca,'XTick',0:n/4:n);\n");
      fprintf(fp,"        set(gca,'XTickLabel', xtick);\n");
      fprintf(fp,"        ytick = L(jj):(U(jj)-L(jj))/4:U(jj);\n");
      fprintf(fp,"        set(gca,'YTick',0:n/4:n);\n");
      fprintf(fp,"        set(gca,'YTickLabel', ytick);\n");
      fprintf(fp,"        set(gca,'YDir', 'normal');\n");
      fprintf(fp,"        xlabel(iStr(jj),'FontWeight','bold','FontSize',12)\n");
      fprintf(fp,"        ylabel(iStr(ii),'FontWeight','bold','FontSize',12)\n");
      fwritePlotAxesNoGrid(fp);
   }
   fprintf(fp, "      end;\n");
   fprintf(fp, "    end;\n");
   fprintf(fp, "  end;\n");
   fprintf(fp, "end;\n");
   double minData=PSUADE_UNDEFINED;
   for (int iChain = 0; iChain < nChains; iChain++) 
   { 
      if (chainStatus[iChain] == 0)
      {
         for (jj = 0; jj < chainCnt; jj++) 
         {
            ddata = XChains[iChain][jj][nInputs];
            if (ddata < minData) minData = ddata;
         }
      }
   }
   fprintf(fp, "      subplot(nPlots,nPlots,1)\n");
   if (psPlotTool_ == 0)
   {
      fprintf(fp,"set(gcf,'NextPlot','add');\n");
      fprintf(fp,"axes;\n");
      fprintf(fp,"h=title('MCMC Posterior Distributions, best -log(likelihood=%e)',",
              minData);
      fprintf(fp,"'fontSize',12,'fontWeight','bold');\n");
      fprintf(fp,"set(gca,'Visible','off');\n");
      fprintf(fp,"set(h,'Visible','on');\n");
   }
   fprintf(fp,"negll = %e;\n", minData);
   fclose(fp);
   if (psPlotTool_ == 1)
        printOutTS(PL_INFO, "MCMC: scilabmcmc2.sci file has been created.\n");
   else printOutTS(PL_INFO, "MCMC: matlabmcmc2.m file has been created.\n");
   return 0;
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
      printOutTS(PL_ERROR, "MCMCAnalyzer ERROR: setParams - not valid.\n");
      exit(1);
   }
   return 0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
MCMCAnalyzer& MCMCAnalyzer::operator=(const MCMCAnalyzer &)
{
   printOutTS(PL_ERROR, "MCMCAnalyzer operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

// ************************************************************************
// functions for getting results
// ------------------------------------------------------------------------
int MCMCAnalyzer::get_nInputs()
{
   return nInputs_;
}
int MCMCAnalyzer::get_nOutputs()
{
   return nOutputs_;
}
double *MCMCAnalyzer::get_means()
{
   double* retVal = NULL;
   if (means_)
   {
      retVal = new double[nInputs_];
      std::copy(means_, means_+nInputs_, retVal);
   }
   return retVal;
}
double *MCMCAnalyzer::get_mostLikelyInput()
{
   double* retVal = NULL;
   if (mostLikelyInput_)
   {
      retVal = new double[nInputs_];
      std::copy(mostLikelyInput_, mostLikelyInput_+nInputs_, retVal);
   }
   return retVal;
}
double *MCMCAnalyzer::get_mostLikelyOutput()
{
   double* retVal = NULL;
   if (mostLikelyOutput_)
   {
      retVal = new double[nOutputs_];
      std::copy(mostLikelyOutput_, mostLikelyOutput_+nOutputs_, retVal);
   }
   return retVal;
}
