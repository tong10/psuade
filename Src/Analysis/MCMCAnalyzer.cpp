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
// ************************************************************************
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

#include "MCMCAnalyzer.h"
#include "Util/sysdef.h"
#include "Util/PsuadeUtil.h"
#include "DataIO/FunctionInterface.h"
#include "DataIO/pData.h"
#include "PDFLib/PDFBase.h"
#include "PDFLib/PDFNormal.h"
#include "PDFLib/PDFLogNormal.h"
#include "PDFLib/PDFTriangle.h"
#include "PDFLib/PDFBeta.h"
#include "PDFLib/PDFWeibull.h"
#include "PDFLib/PDFGamma.h"
#include "FuncApprox/FuncApprox.h"
#include "FuncApprox/Mars.h"
#include "Main/Psuade.h"
#include "Samplings/Sampling.h"
#include "Samplings/RSConstraints.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
MCMCAnalyzer::MCMCAnalyzer() : Analyzer()
{
   setName("MCMC");
   mode_ = 0; // RS mode (mode = 1 : simulation-based mode)
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
   int    nInputs, nOutputs, nSamples, its, length, status;
   int    ii, ii2, jj, jj2, index, maxPts=100, ngrid, sumBins, ione=1;
   int    nbins, **bins, *scores, iMax, increment, fail=0, NTotal, index2;
   int    maxSamples, burnInSamples, *pdfFlags, printLevel, ****bins2;
   int    nPlots, *plotIndices=NULL, kk, kk2, faType, setCutOff=0, *Ivec;
   int    *rsIndices=NULL, nFail, cnt, freq=10, scheme=0, modelFormFlag=0;
   int    dnSamples=0, dnInputs=0, modelFormScheme=0, rsErrFlag=0;
   double *dSamInputs=NULL, *dSamMeans=NULL, *dSamStdevs=NULL;
   double *XRange, *XGuess, *XDist, *means, *sigmas, Xtemp, *YY;
   double *X, *lower, *upper, ddata, Ymax, *Xmax, *Y;
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
   printDashes(0);
   printf("To gain access to different options, turn on \n");
   printf(" * ana_expert to finetune MCMC parameters, \n");
   printf("   (e.g. sample size for burn-in can be adjusted).\n");
   printf(" * rs_expert to finetune response surface used in MCMC,\n");
   printf(" * printlevel 3 to display more diagnostic information.\n");
   printf(" * printlevel >=5 reserved only for expert diagnostics only.\n");
   printf("Features available in this current version of MCMC:\n");
   printf(" * Support different distributions \n");
   printf("   - turn on ana_expert mode to use other than uniform distributions.\n");
   printf(" * Support multiple outputs (likelihood from multiple outputs)\n");
   printf(" * Support Gaussian likelihood function (default)\n");
   printf(" * Support flat-top likelihood functions (for nOutputs=1 only)\n");
   printf("   - use ana_expert to select\n");
   printf(" * For some response surfaces such as polynomial and legendre\n");
   printf("   polynomials, bootstrapped MARS, and Gaussian process (GP1),\n");
   printf("   the fitting errors will be used to create likelihood functions.\n");
   printf("   - select simulation response surface in the simulation data file.\n");
   printf("   - select discrepancy response surface on the fly.\n");
   printf(" * Some input parameters may be disabled (set to default values)\n");
   printf("   - in case these parameters are not to be calibrated\n");
   printf("   - use rs_index_file in PSUADE data file for response surface\n");
   printf(" * A sample may be generated from the posteriors.\n");
   printf("   - use ana_expert to select this optionn");
   printf(" * Posteriors can be propagated through another output.\n");
   printf("   - use ana_expert to select this option\n");
   printf(" * A discrepancy model can be used for model form uncertainty.\n");
   printf("   - use ana_expert to select this option\n");
   printf(" * Create a file named 'psuade_stop' to gracefully terminate\n");
   printf("   before convergence is achieved.\n");
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
            printf("            not supported.\n");
            printf("            Continue with uniform priors.\n");
            for (ii = 0; ii < nInputs; ii++) pdfFlags[ii] = 0;
         }
      }
   }

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
      printf("MCMC ERROR : Direct simulation mode does not support nOutputs > 1.\n");
      printf("             Only response surface mode supports nOutputs > 1.\n");
      printf("             (nOutputs > 1 can be handled in your program by\n");
      printf("              computing the likelihood function yourself and\n");
      printf("              then use mean of zero and std dev of 1).\n");
      return PSUADE_UNDEFINED;
   }

   if (psAnaExpertMode_ == 1 && nOutputs == 1)
   {
      printEquals(0);
      printf("You have the option to use a flat top (uniform) likelihood\n");
      printf("function (non-informative posterior) instead of the default\n");
      printf("Gaussian likelihood with a nonzero tail. This option may\n");
      printf("sometimes be desirable for handling epistemic uncertainties.\n");
      printf("With this option the lower and upper cutoffs for the simulator\n");
      printf("or response surface output will have to be provided to run the\n");
      printf("optimization. Only nOutputs=1 is here to simplify things.\n");
      printf("NOTE: prescribing cutoff carelessly may result in zero\n");
      printf("      posteriors.\n");
      printf("Select the flat top likelihood function option ? (y or n) ");
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

   if (setCutOff == 0)
   {
      printEquals(0);
      printf("MCMC is going to create the likelihood function for you. You will\n");
      printf("need to provide a data file containing the design parameter values,\n");
      printf("mean and standard deviation of the target data for each output.\n");
      printf("NOTE: the design parameters are useful only if you wish to include\n");
      printf("      a discrepancy function (if you don't know what this is, simply\n");
      printf("      ignore the design parameters (set it to 0).\n");
      printf("IMPORTANT: The p design parameters will be considered the same as the\n");
      printf("           first p parameters in the simulation file.\n");
      printDashes(0);
      printf("The format of the data file is: \n");
      printf("PSUADE_BEGIN\n");
      printf("<number of experiments p> <number of design parameters> <nOutputs n>\n");
      printf("1 <design values...> <Output1 mean> <Output 1 std dev> ... <Output n std dev> \n");
      printf("2 <design values...> <Output1 mean> <Output 1 std dev> ... <Output n std dev> \n");
      printf("...\n");
      printf("p <design values...> <Outputn mean> <Output 1 std dev> ... <Output n std dev> \n");
      printf("PSUADE_END\n");
      printDashes(0);
      printf("NOTE: Alternately, you can compute the chi-squared yourself and then\n");
      printf("      specify nOutputs=1, mean = 0, and standard deviation = 1 in this\n"); 
      printf("      specification file.\n\n");
      printf("Now enter the spec file for building the likelihood function : ");
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
         ifile >> dnInputs;
         if (dnInputs < 0)
         {
            printf("MCMC ERROR: number of experimental design variables < 0.\n");
            ifile.close();
            return PSUADE_UNDEFINED;
         }
         if (dnInputs > nInputs)
         {
            printf("MCMC ERROR: number of experimental design variables %d cannot be\n",
                   dnInputs);
            printf("            larger than number of inputs to the simulator %d.\n",
                   nInputs);
            ifile.close();
            return PSUADE_UNDEFINED;
         }
         ifile >> kk;
         if (kk != nOutputs)
         {
            printf("MCMC ERROR: number of outputs for each experiment does not\n");
            printf("            match the simulator nOutputs.\n");
            ifile.close();
            return PSUADE_UNDEFINED;
         }
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
         printf("MCMC ERROR: PSUADE_BEGIN not found in the first line of the spec file.\n");
         ifile.close();
         return PSUADE_UNDEFINED;
      }
      getline (ifile, iline);
      getline (ifile, iline);
      compFlag = iline.compare(0,10,"PSUADE_END");
      if (compFlag != 0)
      {
         printf("MCMC ERROR: PSUADE_END not found in spec file.\n");
         ifile.close();
         delete [] dSamMeans;
         delete [] dSamStdevs;
         return PSUADE_UNDEFINED;
      }
      ifile.close();

      if (psAnaExpertMode_ == 1 && mode_ == 0)
      {
         printf("You can choose to incorporate the response surface errors into\n");
         printf("the likelihood function. This will only work for response\n");
         printf("surfaces that provide also response surface errors such as\n");
         printf("GP, polynomial regression, bootstrapped MARS.\n");
         printf("Select this option ? (y or n) ");
         scanf("%s", charString);
         fgets(lineIn,500,stdin);
         if (charString[0] == 'y') rsErrFlag = 1;
      }
   }

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
               printf("  data format should be: \n");
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
                  printf("  data format should be: \n");
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
         }
         if (printLevel >= 4)
         {
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
         printf("MCMC INFO: creating response surface for output %d.\n",ii+1);
         faPtrs[ii] = genFA(faType, nInputs, nSamples);
         faPtrs[ii]->setBounds(lower, upper);
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

   burnInSamples = 10000;
   maxSamples = 10000;
   nbins = 20;
   printEquals(0);
   printf("MCMC Burn-in sample size      (default) = %d\n", burnInSamples);
   printf("MCMC sample increment         (default) = %d\n", maxSamples);
   printf("MCMC no. of bins in histogram (default) = %d\n", nbins);
   printf("Note: sample increment - sample size to run before convergence check\n");
   printf("Note: histogram nBins  - define granularity of histogram bar graph\n");
   printf("Turn on ana_expert mode to change these default settings.\n");
   printEquals(0);
   if (psAnaExpertMode_ == 1)
   {
      sprintf(charString,"Enter sample size for burn in (5000 - 1000000): ");
      burnInSamples = getInt(5000, 1000000, charString);
      sprintf(charString,"Enter sample increment (5000 - 1000000): ");
      maxSamples = getInt(5000, 1000000, charString);
      sprintf(charString,"Enter the number of histogram bins (10 - 25) : ");
      nbins = getInt(10, 25, charString);
      printf("MCMC will create MATLAB files for the posterior distributions.\n");
      printf("You can choose to generate posterior plots of all inputs, or \n");
      printf("just a selected few (in case there are too many inputs).\n");
      printf("Select inputs for which posterior plots are to be generated.\n");
      sprintf(charString,"Enter the input number (%d - %d, 0 to terminate) : ",
              dnInputs, nInputs);
      kk = 1;
      plotIndices = new int[nInputs];
      nPlots = 0;
      while (kk != 0)
      {
         kk = getInt(0, nInputs, charString);
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
      printf("MCMC Plot summary: input number to be plotted are:\n");
      for (ii = 0; ii < nPlots; ii++)
         printf("   Input %4d\n", plotIndices[ii]+1);

      printEquals(0);
      printf("You have the option to add a discrepancy function to either\n");
      printf("the simulator or the response surface used for generating the\n");
      printf("likelihood function for Bayesian inference. If this option is\n");
      printf("selected, you will need to make sure the data file you provided\n");
      printf("earlier has design parameters specified.\n");
      printf("ONE REQUIREMENT: the m design parameters in the likelihood file\n");
      printf("      need to correspond to the first m inputs in the simulator.\n");
      printf("Select this option ? (y or n) ");
      scanf("%s", charString);
      fgets(lineIn,500,stdin);
      if (charString[0] == 'y')
      {
         modelFormFlag = 1;
#if 1
         printf("Please select among from the following discrepancy function forms:\n");
         printf("1. d(X), where X are design parameters specified in likelihood data file.\n");
         printf("   In this case, the first %d parameters are design variables.\n",dnInputs);
         printf("   and the rest are calibration parameters (except the ones fixed by\n");
         printf("   the rs_index list.)\n");
         printf("2. d(X,theta), where theta in addition are the remaining %d calibration\n",
                nInputs-dnInputs);
         printf("   parameters (except the ones fixed by the rs_index list.)\n");
         sprintf(charString,"Make your selection: (1 or 2) ");
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
      printf("In addition to generating the posterior distributions, you can\n");
      printf("also draw a sample from these posteriors. The posterior sample\n");
      printf("can be used as prior sample for another simulator/emulator.\n");
      printf("Select this option ? (y or n) ");
      scanf("%s", charString);
      fgets(lineIn,500,stdin);
      if (charString[0] == 'y')
      {
         fp2 = fopen("MCMCPostSample", "w");
         if (fp2 == NULL)
         {
            printf("MCMC ERROR: cannot open MCMCPostSample file.\n");
            printf("            Check permission on local directory.\n");
            exit(1);
         }
         fprintf(fp2,"%% first %d columns - input, last column - dummpy output.\n",
                 nInputs);
      }
      printAsterisks(0);
   } 
   else
   {
      nPlots = nInputs - dnInputs;
      plotIndices = new int[nPlots];
      for (ii = dnInputs; ii < nInputs; ii++) plotIndices[ii-dnInputs] = ii;
      fp2 = NULL;
   }

   inputPDFs = new PDFBase*[nInputs];
   for (ii = 0; ii < nInputs; ii++)
   {
      if (pdfFlags != NULL && pdfFlags[ii] == PSUADE_PDF_NORMAL)
      {
         inputPDFs[ii] = (PDFBase *) new PDFNormal(inputMeans[ii], inputStdevs[ii]);
         if (printLevel > 2) 
            printf("Parameter %3d has normal prior distribution (%e,%e).\n",
                   ii, inputMeans[ii], inputStdevs[ii]);
      }
      else if (pdfFlags != NULL && pdfFlags[ii] == PSUADE_PDF_LOGNORMAL)
      {
         inputPDFs[ii] = (PDFBase *) new PDFLogNormal(inputMeans[ii],inputStdevs[ii]);
         if (printLevel > 2)
            printf("Parameter %3d has lognormal prior distribution.\n",ii);
      }
      else if (pdfFlags != NULL && pdfFlags[ii] == PSUADE_PDF_TRIANGLE)
      {
         inputPDFs[ii] = (PDFBase *) new PDFTriangle(inputMeans[ii],inputStdevs[ii]);
         if (printLevel > 2)
            printf("Parameter %3d has triangle prior distribution.\n",ii);
      }
      else if (pdfFlags != NULL && pdfFlags[ii] == PSUADE_PDF_BETA)
      {
         inputPDFs[ii] = (PDFBase *) new PDFBeta(inputMeans[ii],inputStdevs[ii]);
         if (printLevel > 2)
            printf("Parameter %3d has beta prior distribution.\n",ii);
      }
      else if (pdfFlags != NULL && pdfFlags[ii] == PSUADE_PDF_WEIBULL)
      {
         inputPDFs[ii] = (PDFBase *) new PDFWeibull(inputMeans[ii],inputStdevs[ii]);
         if (printLevel > 2)
            printf("Parameter %3d has Weibull prior distribution.\n",ii);
      }
      else if (pdfFlags != NULL && pdfFlags[ii] == PSUADE_PDF_GAMMA)
      {
         inputPDFs[ii] = (PDFBase *) new PDFGamma(inputMeans[ii],inputStdevs[ii]);
         if (printLevel > 2)
            printf("Parameter %3d has gamma prior distribution.\n",ii);
      }
      else
      {
         inputPDFs[ii] = NULL;
         if (printLevel > 2)
            printf("Parameter %3d has uniform prior distribution.\n",ii);
      }
   }

   if (mode_ == 1 && ioPtr != NULL)
   {
      funcIO = createFunctionInterface(ioPtr);
      printf("MCMC: direct simulation (not response surface) has been set up.\n");
      printf("MCMC INFO: make sure simulation output is in the right form\n");
      printf("     which should not have been translated (mean) nor scaled (std)\n");
      printf("     unless you compute the chi-square yourself in which case\n");
      printf("     you should have use nOutputs=1, mean=0, and std dev=1 in the\n");
      printf("     spec file.\n");
      if (psAnaExpertMode_ == 1 && nSamples > 0)
      {
         printf("Since MCMC uses many function evaluations, you have the option\n");
         printf("of setting how frequent the simulator (rather than the response\n");
         printf("surface) is used in creating the proposal distribution. This is\n");
         printf("sometimes needed if the response surface may not be accurate enough\n");
         printf("but the simulator is expensive to run. A frequency of f means that\n");
         printf("each MCMC step uses f simulator runs. Default is f=10.\n"); 
         sprintf(charString, "Enter frequency f (1 - %d): ", maxPts);
         kk = getInt(1, maxPts, charString);
         freq = maxPts * nInputs / kk;
         if (freq * kk != maxPts) freq++;
      }
      else if (nSamples > 0) freq = maxPts / 10;
      else                   freq = maxPts;
   }
   if (funcIO == NULL && faPtrs == NULL)
   {
      printf("MCMC ERROR: missing simulator and sample data - cannot proceed.\n");
      return PSUADE_UNDEFINED;
   }

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
      printf("Select response surface type for the discrepancy function.\n");
      while (dfaType < 0 || dfaType >= PSUADE_NUM_RS)
      {
         writeFAInfo();
         sprintf(charString, "Please enter your choice ? ");
         dfaType = getInt(0, PSUADE_NUM_RS-1, charString);
      }

      faPtrs1 = new FuncApprox*[nOutputs];
      dOneSample = new double[nInputs];
      tSamInputs = new double[nSamples*nInputs];
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
                  for (ii2 = dnInputs; ii2 < nInputs; ii2++)
                     dOneSample[ii2] = settings[ii2];
                  if (funcIO != NULL)
                       funcIO->evaluate(kk+1,nInputs,dOneSample,1,&simdata,0);
                  else simdata = faPtrs[ii]->evaluatePoint(dOneSample);
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
            for (kk = 0; kk < nSamples; kk++)
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
            sprintf(charString, "DiscrepancySample%d", ii+1);
            dataPtr->writePsuadeFile(charString, 0);
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
            sprintf(charString, "DiscrepancySample%d", ii+1);
            dataPtr->writePsuadeFile(charString, 0);
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
      delete [] LHSamInputs;
      delete [] dOneSample;
      delete [] LHSamOutputs;
      delete [] settings;
      delete [] tSamInputs;
   }

   printf("MCMC INFO : creating constraints, if there is any.\n");
   printf("            Constraints remove infeasible regions from the priors.\n");
   printf("            Constraints can be specified by RS constraint files.\n");
   constrPtr = new RSConstraints();
   constrPtr->genConstraints(ioPtr);

   XRange = new double[nInputs];
   XGuess = new double[nInputs];
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

   Xmax = new double[nInputs];
   Ivec = new int[nInputs];
   Ivec[nInputs-1] = -1;
   for (ii = 0; ii < nInputs; ii++) Xmax[ii] = 0.0;
   Ymax = 0.0;
   XDist = new double[maxPts];
   printf("MCMC Phase 1: burn in \n");
   fflush(stdout);
   increment = burnInSamples / nInputs;
   for (its = 0; its < increment; its++)
   {
      if ((its+1) % (increment/10) == 0)
      {
         printf("%3.0f%% ", 100.0 * ((its+1.0) / increment));
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
             (rsIndices != NULL && rsIndices[ii] >=0)) && ii >= dnInputs)
         {
            Xtemp = XGuess[ii];
            nFail = 0;
            cnt = 0;
            for (jj = 0; jj < maxPts; jj++)
            {
               XGuess[ii] = lower[ii]+jj*XRange[ii]/(maxPts-1.0);

               if (nOutputs > 1)
               {
                  XDist[jj] = 0.0;
                  for (kk2 = 0; kk2 < dnSamples; kk2++) 
                  {
                     for (ii2 = 0; ii2 < dnInputs; ii2++) 
                        XGuess[ii2] = dSamInputs[kk2*dnInputs+ii2];
                     for (ii2 = 0; ii2 < nOutputs; ii2++) 
                     {
                        Ytemp = faPtrs[ii2]->evaluatePoint(XGuess);
                        if (faPtrs1 != NULL && faPtrs1[ii2] != NULL)
                           Ytemp += faPtrs1[ii2]->evaluatePoint(XGuess);
                        Ytemp = pow((Ytemp-dSamMeans[kk2*nOutputs+ii2])/
                                    dSamStdevs[kk2*nOutputs+ii2],2.0);
                        XDist[jj] += Ytemp;
                     }
                  }
                  XDist[jj] = exp(-0.5 * XDist[jj] / dnSamples);
               }
               else
               {
                  cnt++;
                  if (cnt > freq) cnt = 0; 
                  XDist[jj] = 0.0;
                  for (kk2 = 0; kk2 < dnSamples; kk2++) 
                  {
                     for (ii2 = 0; ii2 < dnInputs; ii2++) 
                        XGuess[ii2] = dSamInputs[kk2*dnInputs+ii2];
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
                        Ytemp = pow((Ytemp-dSamMeans[kk2])/dSamStdevs[kk2],2.0);
                        XDist[jj] += Ytemp;
                     }
                  }
                  XDist[jj] = exp(-0.5 * XDist[jj] / dnSamples);
               }

               Ytemp = constrPtr->evaluate(XGuess,XDist[jj],status);
               if (status == 0)
               {
                  XDist[jj] = 0.0;
                  nFail++;
               }
               if (inputPDFs != NULL && inputPDFs[ii] != NULL)
                    inputPDFs[ii]->getPDF(ione,&XGuess[ii],&ddata);
               else ddata = 1.0;
               XDist[jj] *= ddata; 
               if (XDist[jj] > Ymax)
               {
                  Ymax = XDist[jj];
                  for (ii2 = 0; ii2 < nInputs; ii2++) Xmax[ii2] = XGuess[ii2];
               }
               if (jj > 0) XDist[jj] += XDist[jj-1];
            }
            if (XDist[maxPts-1] > 0.0)
            {
               for (jj = 0; jj < maxPts; jj++) XDist[jj] /= XDist[maxPts-1];
               XGuess[ii] = PSUADE_drand();
               index = binarySearchDble(XGuess[ii], XDist, maxPts);
               if (index < 0) index = - index - 1;
               if (index == maxPts-1) ddata = (double) index;
               else
                  ddata=index +
                     (XGuess[ii]-XDist[index])/(XDist[index+1]-XDist[index]);
               XGuess[ii] = lower[ii]+ddata*XRange[ii]/(maxPts-1.0);
            }
            else
            {
               XGuess[ii] = Xtemp;
               if (printLevel > 4) 
                  printf("MCMC iteration %7d : no modification in input %d\n",
                         its+1,ii+1);
               if (nFail == maxPts) 
                  printf("INFO : Constraints have resulted in zero distribution.\n");
            }
         }
      }
   }
   printf("\n");

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
   printf("MCMC Phase 2: create posterior (10.)\n");
   fflush(stdout);
   increment = maxSamples / nInputs;
   while (1)
   {
      if (((its+1) % (increment/10)) == 0)
      {
         printf("%3.0f%% ",100.0*(((its % increment)+ 1.0) / increment));
         fflush(stdout);
      }

      if (printLevel > 4) 
      {
         if (psPlotTool_ == 1) fp = fopen("psTrack.sci", "w");
         else                  fp = fopen("psTrack.m", "w");
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
               for (jj = 0; jj < maxPts; jj++)
               {
                  XGuess[ii] = lower[ii]+jj*XRange[ii]/(maxPts-1.0);
   
                  if (nOutputs > 1)
                  {
                     XDist[jj] = 0.0;
                     for (kk2 = 0; kk2 < dnSamples; kk2++) 
                     {
                        for (ii2 = 0; ii2 < dnInputs; ii2++) 
                           XGuess[ii2] = dSamInputs[kk2*dnInputs+ii2];
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
                     XDist[jj] = exp(-0.5 * XDist[jj] / dnSamples);
                  }
                  else
                  {
                     cnt++;
                     if (cnt > freq) cnt = 0; 
                     XDist[jj] = 0.0;
                     for (kk2 = 0; kk2 < dnSamples; kk2++) 
                     {
                        for (ii2 = 0; ii2 < dnInputs; ii2++) 
                           XGuess[ii2] = dSamInputs[kk2*dnInputs+ii2];
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

                  Ytemp = constrPtr->evaluate(XGuess,XDist[jj],status);
                  if (status == 0)
                  {
                     XDist[jj] = 0.0;
                     nFail++;
                  }

                  if (inputPDFs != NULL && inputPDFs[ii] != NULL)
                       inputPDFs[ii]->getPDF(ione,&XGuess[ii],&ddata);
                  else ddata = 1.0;
                  XDist[jj] *= ddata; 

                  if (XDist[jj] > Ymax)
                  {
                     Ymax = XDist[jj];
                     for (ii2 = 0; ii2 < nInputs; ii2++) Xmax[ii2] = XGuess[ii2];
                  }
                  if (jj > 0) XDist[jj] += XDist[jj-1];
               }

               if (XDist[maxPts-1] > 0.0)
               {
                  for (jj = 0; jj < maxPts; jj++) XDist[jj] /= XDist[maxPts-1];
                  if (printLevel > 4)
                  {
                     fprintf(fp, "A%d = [\n", ii+1);
                     for (jj = 0; jj < maxPts; jj++) fprintf(fp, "%e\n",XDist[jj]);
                     fprintf(fp, "];\n");
                  }
                  XGuess[ii] = PSUADE_drand();
                  index = binarySearchDble(XGuess[ii], XDist, maxPts);
                  if (index < 0) index = - index - 1;
                  if (index == maxPts-1) ddata = (double) index;
                  else
                     ddata=index+
                        (XGuess[ii]-XDist[index])/(XDist[index+1]-XDist[index]);
                  XGuess[ii] = lower[ii]+ddata*XRange[ii]/(maxPts-1.0);
               }
               else 
               {
                  XGuess[ii] = Xtemp;
                  if (printLevel > 4) 
                     printf("MCMC iteration %7d : no modification in input %d\n",
                            its+1,ii+1);
                  if (nFail == maxPts) 
                     printf("INFO : Constraints have resulted in zero distribution.\n");
               }
            }
            if (printLevel > 4)
            {
               ngrid = (int) (pow(1.0*nInputs, 0.5001) + 1.0); 
               if (psPlotTool_ == 1) fprintf(fp, "clf();\n");
               else                  fprintf(fp, "clf\n");
               for (ii = 0; ii < nInputs; ii++)
               {
                  fprintf(fp, "subplot(%d,%d,%d)\n", ngrid, ngrid, ii+1);
                  fprintf(fp, "X%d = [\n", ii+1);
                  for (jj = 0; jj < maxPts; jj++) 
                     fprintf(fp, "%e\n",
                             lower[ii]+jj*XRange[ii]/(maxPts-1.0));
                  fprintf(fp, "];\n");
                  fprintf(fp, "A%d = [\n", ii+1);
                  for (jj = 0; jj < maxPts; jj++) fprintf(fp, "%e\n", XDist[jj]);
                  fprintf(fp, "];\n");
                  fprintf(fp, "plot(X%d,A%d)\n", ii+1, ii+1);
                  fwritePlotAxes(fp);
               }
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
         for (jj = dnInputs; jj < nInputs; jj++)
            if (rsIndices != NULL && rsIndices[jj] >= 0)
               fprintf(fp, "%e ", XGuess[jj]);
         fprintf(fp2, "1.0e35\n");
      }

      if (printLevel > 4)
      {
         if (psPlotTool_ == 1) 
            printf("MCMCAnalyzer: diagnostics file is in psTrack.sci\n");
         else
            printf("MCMCAnalyzer: diagnostics file is in psTrack.m\n");
         fclose(fp);
         printf("Enter any character and return to continue.\n");
         scanf("%s", charString);
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
                  printf("MCMC: input %3d statistics (mean,stdev) unavailable\n",
                         ii+1);
                  break;
               }
               else
                  printf("MCMC: input %3d statistics (mean,stdev) = %e %e\n",
                         ii+1, means[ii]/NTotal, 
                         sqrt(sigmas[ii]/NTotal-pow(means[ii]/NTotal,2.0)));
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
               fprintf(fp, "xmin = min(X%d);\n", ii+1);
               fprintf(fp, "xmax = max(X%d);\n", ii+1);
               fprintf(fp, "bar(X%d, N%d, 1.0);\n", ii+1, ii+1);
               fprintf(fp, "a = gce();\n");
               fprintf(fp, "a.children.thickness = 2;\n");
               fprintf(fp, "a.children.foreground = 0;\n");
               fprintf(fp, "a.children.background = 2;\n");
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
   if (NTotal == 0)
   {
      printf("MCMC FAILS: probability distribution may be too narrow.\n");
      printf("            Try again with a larger standard deviation.\n");
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
            fprintf(fp, "xmin = min(X%d);\n", ii+1);
            fprintf(fp, "xmax = max(X%d);\n", ii+1);
            fprintf(fp, "bar(X%d, N%d, 1.0);\n", ii+1, ii+1);
            fprintf(fp, "a = gce();\n");
            fprintf(fp, "a.children.thickness = 2;\n");
            fprintf(fp, "a.children.foreground = 0;\n");
            fprintf(fp, "a.children.background = 2;\n");
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
      if (psPlotTool_ == 1)
         printf("MCMC: final scilabmcmc.sci file has been created.\n");
      else
         printf("MCMC: final matlabmcmc.m file has been created.\n");

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
                  fprintf(fp, "xmin = min(X%d);\n", ii+1);
                  fprintf(fp, "xmax = max(X%d);\n", ii+1);
                  fprintf(fp, "bar(X%d, N%d, 1.0);\n",ii+1,ii+1);
                  fprintf(fp, "a = gce();\n");
                  fprintf(fp, "a.children.thickness = 2;\n");
                  fprintf(fp, "a.children.foreground = 0;\n");
                  fprintf(fp, "a.children.background = 2;\n");
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
      printf("\n");
   }

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
      fclose(fp2);
      printf("MCMC: check the MCMCPostSample file for a posterior sample.\n");
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

