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
// Functions for the class TwoParamAnalyzer 
// AUTHOR : CHARLES TONG
// DATE   : updated in 2006
// ************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "TwoParamAnalyzer.h"
#include "Psuade.h"
#include "PsuadeUtil.h"
#include "sysdef.h"
#include "Matrix.h"
#include "PsuadeData.h"
#include "pData.h"
#include "RSConstraints.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
TwoParamAnalyzer::TwoParamAnalyzer() : Analyzer()
{
   setName("TwoParam");
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
TwoParamAnalyzer::~TwoParamAnalyzer()
{
}

// ************************************************************************
// perform VCE analysis for two parameters
// ------------------------------------------------------------------------
double TwoParamAnalyzer::analyze(aData &adata)
{
   int    nInputs, nOutputs, nSamples, outputID, nSubSamples;
   int    ii, jj, ss, repID, subID, nReps1, ncount, index, status;
   int    symIndex, whichOutput, ii2, ss2, nReps, ir, *bins;
   int    errflag, nSubs, printLevel, iZero=0, iOne=1;
   int    totalCnt, validnReps, *pdfFlags;
   double *vce, *meanVCEVar, *vceMean, *vceVariance, *varVCEMean;
   double fixedInput, aMean, aVariance, **bsVCEs, ddata;
   double *txArray, *tyArray, *twArray, fixedInput2;
   double *varVCEVar, *mainEffects, *STI, *xLower, *xUpper;
   double *XIn, *YIn, *XX, *YY, *X, *Y;
   char   pString[500], winput[501];
   FILE   *fp, *fp1=NULL;
   PsuadeData    *ioPtr=NULL;
   RSConstraints *constrPtr=NULL;
   pData         pCorMat;
   Matrix        *corMatp, corMat;

   nInputs  = adata.nInputs_;
   nOutputs = adata.nOutputs_;
   nSamples = adata.nSamples_;
   xLower   = adata.iLowerB_;
   xUpper   = adata.iUpperB_;
   XIn      = adata.sampleInputs_;
   YIn      = adata.sampleOutputs_;
   outputID    = adata.outputID_;
   printLevel  = adata.printLevel_;
   nSubSamples = adata.nSubSamples_;
   ioPtr       = adata.ioPtr_;
   pdfFlags    = adata.inputPDFs_;

   if (nInputs <= 0 || nOutputs <= 0 || nSamples <= 0)
   {
      printf("Two-way Effect ERROR: invalid arguments.\n");
      printf("    nInputs  = %d\n", nInputs);
      printf("    nOutputs = %d\n", nOutputs);
      printf("    nSamples = %d\n", nSamples);
      return PSUADE_UNDEFINED;
   } 
   if (nSamples/nSubSamples*nSubSamples != nSamples)
   {
      printf("Two-way Effect ERROR: nSamples != k*nSubSamples.\n");
      printf("    nSamples    = %d\n", nSamples);
      printf("    nSubSamples = %d\n", nSubSamples);
      return PSUADE_UNDEFINED;
   } 
   whichOutput = outputID;
   if (whichOutput >= nOutputs || whichOutput < 0)
   {
      printf("Two-way Effect ERROR: invalid outputID.\n");
      printf("    nOutputs = %d\n", nOutputs);
      printf("    outputID = %d\n", whichOutput+1);
      printf("    outputID reset to 1\n");
      whichOutput = 0;
   }
   if (ioPtr == NULL)
   {
      printf("Two-way Effect ERROR: no data (PsuadeData).\n");
      return PSUADE_UNDEFINED;
   }
   ioPtr->getParameter("input_cor_matrix", pCorMat);
   corMatp = (Matrix *) pCorMat.psObject_;
   for (ii = 0; ii < nInputs; ii++)
   {
      for (jj = 0; jj < ii; jj++)
      {
         if (corMatp->getEntry(ii,jj) != 0.0)
         {
            printf("* Two-way Effect INFO: this method should not be used\n");
            printf("*                  if inputs are correlated with joint\n");
            printf("*                  PDFs. Use group variance-based method.\n");
            return PSUADE_UNDEFINED;
         }     
      }
   }
   status = 0;
   for (ss = 0; ss < nSamples; ss++)
      if (YIn[nOutputs*ss+whichOutput] > 0.9*PSUADE_UNDEFINED) status = 1;
   if (status == 1)
   {
      printf("Two-way Effect ERROR: Some outputs are undefined. Prune\n");
      printf("                           the undefined sample points first.\n");
      return PSUADE_UNDEFINED;
   }

   for (ii = 0; ii < nInputs; ii++)
   {
      if (pdfFlags != NULL && pdfFlags[ii] != PSUADE_PDF_UNIFORM)
      {
         printf("* Two-way Effect INFO: some inputs have non-uniform PDFs.\n");
         printf("*           However, they will not be relevant in this analysis\n");
         printf("*           (since the sample should have been generated with\n");
         printf("*            the desired distributions.)\n");
         break;
      }
   }

   if (ioPtr != NULL)
   {
      constrPtr = new RSConstraints();
      constrPtr->genConstraints(ioPtr);
   }
   X = XIn;
   Y = new double[nSamples];
   ncount = 0;
   for (ss = 0; ss < nSamples; ss++)
   {
      Y[ss] = YIn[nOutputs*ss+whichOutput];
      ddata = constrPtr->evaluate(&(X[ss*nInputs]), Y[ss], status);
      if (status == 0) Y[ss] = PSUADE_UNDEFINED;
      else             ncount++;
   }
   if (ncount == 0)
   {
      printf("Two-way Effect ERROR: no valid sample point.\n");
      printf("    nSamples before filtering = %d\n", nSamples);
      printf("    nSamples after  filtering = %d\n", ncount);
      printf("    INFO: check your data file for undefined's (1e35)\n");
      return 1.0;
   }
   if (ncount != nSamples)
   {
      printf("* Two-way Effect: nSamples before filtering = %d\n", nSamples);
      printf("* Two-way Effect: nSamples after filtering  = %d\n", ncount);
   }

   computeMeanVariance(nInputs,1,nSamples,Y,&aMean,&aVariance,0);
   printf("* =====> Two-way Effect: mean     = %12.4e\n", aMean);
   printf("* =====> Two-way Effect: variance = %12.4e (sd =%12.4e)\n", 
          aVariance, sqrt(aVariance));
   if (PABS(aVariance) < 1.0e-15)
   {
      printf("Two-way Effect INFO: std dev too small ==> terminate.\n");
      return PSUADE_UNDEFINED;
   }
   if (psAnaExpertMode_ == 1 && printLevel > 4)
   {
      printf("Perform crude analysis? (y or n) ");
      fgets(pString,400,stdin);
      if (pString[0] == 'y')
      {
         vce = new double[nInputs*nInputs];
         status = computeVCECrude(nInputs, nSamples, X, Y, aVariance,
                                  xLower, xUpper, vce);
         if (status == 0)
         {
            pData *pPtr = ioPtr->getAuxData();
            pPtr->nDbles_ = nInputs * nInputs;
            pPtr->dbleArray_ = new double[nInputs * nInputs];
            for (ii = 0; ii < nInputs*nInputs; ii++)
               pPtr->dbleArray_[ii] = vce[ii];
            pPtr->dbleData_ = aVariance;
         }
         delete [] Y;
         delete [] vce;
         return 1.0;
      }
   }

   vceMean       = new double[nSubSamples];
   vceVariance   = new double[nSubSamples];
   bins          = new int[nSubSamples];
   vce           = new double[nInputs*nInputs];
   varVCEMean    = new double[nInputs*nInputs];
   varVCEVar     = new double[nInputs*nInputs];
   meanVCEVar    = new double[nInputs*nInputs];
   mainEffects   = new double[nInputs];
   STI           = new double[nInputs];
   errflag       = 0;

   printf("\n");
   printAsterisks(0);
   printAsterisks(0);
   printf("* Second order Sensitivities\n");
   printf("* Note: This method can only be accurate for large sample size.\n");
   printf("* Note: For random samples, a crude analysis will be performed.\n");
   printf("*       (use analysis expert mode to set number of bins).\n");
   printDashes(0);
   printf("* total number of samples = %10d\n",nSamples);
   printf("* number of nInputs       = %10d\n",nInputs);
   printDashes(0);
   printf("Output %d\n", whichOutput+1);

   txArray  = new double[nSamples];
   twArray  = new double[nSamples];
   tyArray  = new double[nSamples];
   if (tyArray == NULL)
   {
      printf("Two-way Effect ERROR:: memory allocation problem.\n");
      printf("                       Consult PSUADE developers.\n");
   }

   for (ii = 0; ii < nInputs; ii++)
   {
      for (ii2 = ii+1; ii2 < nInputs; ii2++)
      {
         if (printLevel > 4 || nSamples > 100000)
            printf("Two-way Effect:: input pairs = %d %d\n", ii+1, ii2+1);
         for (ss = 0; ss < nSamples; ss++)
         {
            txArray[ss] = X[nInputs*ss+ii];
            twArray[ss] = X[nInputs*ss+ii2];
            tyArray[ss] = Y[ss];
         }

         sortDbleList3(nSamples, txArray, twArray, tyArray);

         for (ss = 1; ss < nSamples; ss++)
            if (PABS(txArray[ss]-txArray[ss-1]) > 1.0E-10)
               break;
         nReps1 = ss;
         if (nReps1 <= 1)
         {
            printf("Two-way Effect INFO: nReps < 1 for input %d.\n",ii+1);
            printf("        ==> not replicated orthogonal array nor Factorial\n");
            printf("        ==> crude 2-way interaction analysis.\n");
            status = computeVCECrude(nInputs, nSamples, X, Y, aVariance,
                                     xLower, xUpper, vce);
            if (status == 0)
            {
               pData *pPtr = ioPtr->getAuxData();
               pPtr->nDbles_ = nInputs * nInputs;
               pPtr->dbleArray_ = new double[nInputs * nInputs];
               for (ii = 0; ii < nInputs*nInputs; ii++) 
                  pPtr->dbleArray_[ii] = vce[ii];
               pPtr->dbleData_ = aVariance;
            }
            delete [] vce;
            delete [] vceMean;
            delete [] vceVariance;
            delete [] varVCEMean;
            delete [] varVCEVar;
            delete [] meanVCEVar;
            delete [] mainEffects;
            delete [] bins;
            delete [] STI;
            delete [] txArray;
            delete [] tyArray;
            delete [] twArray;
            delete [] Y;
            return 1.0;
         }

         for (ss = 0; ss < nSamples; ss+=nReps1)
         {
            fixedInput = txArray[ss];
            ncount = 0;
            for (ss2 = 0; ss2 < nReps1; ss2++)
               if (PABS(txArray[ss+ss2]-fixedInput) < 1.0E-10)
                  ncount++;
            if (ncount != nReps1) errflag = 1;
         }
         if (errflag != 0)
         {
            printf("Two-way Effect WARNING(2): replications not satisified.\n");
            printf("    Are you using replicated orthogonal array or Factorial?\n");
            printf("    If so, you need to use > 1 replications. A crude\n");
            printf("    2-way interaction analysis will be done instead.\n");
            status = computeVCECrude(nInputs, nSamples, X, Y, aVariance,
                                     xLower, xUpper, vce);
            if (status == 0)
            {
               pData *pPtr = ioPtr->getAuxData();
               pPtr->nDbles_ = nInputs * nInputs;
               pPtr->dbleArray_ = new double[nInputs * nInputs];
               for (ii = 0; ii < nInputs*nInputs; ii++) 
                  pPtr->dbleArray_[ii] = vce[ii];
               pPtr->dbleData_ = aVariance;
            }
            delete [] vce;
            delete [] vceMean;
            delete [] vceVariance;
            delete [] varVCEMean;
            delete [] varVCEVar;
            delete [] meanVCEVar;
            delete [] mainEffects;
            delete [] bins;
            delete [] STI;
            delete [] txArray;
            delete [] tyArray;
            delete [] twArray;
            delete [] Y;
            return 1.0;
         }

         for (ss = 0; ss < nSamples; ss+=nReps1)
            sortDbleList2(nReps1,&twArray[ss],&tyArray[ss]);

         for (ss = 1; ss < nReps1; ss++)
            if (PABS(twArray[ss]-twArray[ss-1]) > 1.0E-10)
               break;
         nReps = ss;
         if (nReps <= 1 || (nReps1/nReps*nReps != nReps1))
         {
            printf("Two-way Effect WARNING(3): replications not satisified.\n");
            printf("    Are you using replicated orthogonal array or Factorial?\n");
            printf("    If so, you need to use > 1 replications.\n");
            printf("A Crude 2-way interaction analysis will be done instead.\n");
            status = computeVCECrude(nInputs, nSamples, X, Y, aVariance,
                                     xLower, xUpper, vce);
            if (status == 0)
            {
               pData *pPtr = ioPtr->getAuxData();
               pPtr->nDbles_ = nInputs * nInputs;
               pPtr->dbleArray_ = new double[nInputs * nInputs];
               for (ii = 0; ii < nInputs*nInputs; ii++)   
                  pPtr->dbleArray_[ii] = vce[ii];
               pPtr->dbleData_ = aVariance;
            }
            delete [] vce;
            delete [] vceMean;
            delete [] vceVariance;
            delete [] varVCEMean;
            delete [] varVCEVar;
            delete [] meanVCEVar;
            delete [] mainEffects;
            delete [] bins;
            delete [] STI;
            delete [] txArray;
            delete [] tyArray;
            delete [] twArray;
            delete [] Y;
            return 1.0;
         }

         for (ss = 0; ss < nSamples; ss+=nReps1)
         {
            for (ss2 = 0; ss2 < nReps1; ss2+=nReps)
            {
               ncount = 1;
               fixedInput = twArray[ss+ss2];
               for (repID = 1; repID < nReps; repID++)
               {
                  if (PABS(twArray[ss+ss2+repID]-fixedInput)<1.0E-10)
                  ncount++;
               }
               if (ncount != nReps)
               {
                  printf("Two-way Effect WARNING(3): sample not rOA/FACT.\n");
                  errflag = 1;
               }
            }
            if (errflag != 0) break;
         }
         if (errflag != 0)
         {
            printf("Two-way Effect WARNING(4): replications not satisified.\n");
            printf("    Are you using replicated orthogonal array or Factorial?\n");
            printf("    If so, you need to use > 1 replications.\n");
            printf("A Crude 2-way interaction analysis will be done instead.\n");
            status = computeVCECrude(nInputs, nSamples, X, Y, aVariance,
                                     xLower, xUpper, vce);
            if (status == 0)
            {
               pData *pPtr = ioPtr->getAuxData();
               pPtr->nDbles_ = nInputs * nInputs;
               pPtr->dbleArray_ = new double[nInputs * nInputs];
               for (ii = 0; ii < nInputs*nInputs; ii++)
                  pPtr->dbleArray_[ii] = vce[ii];
               pPtr->dbleData_ = aVariance;
            }
            delete [] vce;
            delete [] vceMean;
            delete [] vceVariance;
            delete [] varVCEMean;
            delete [] varVCEVar;
            delete [] meanVCEVar;
            delete [] mainEffects;
            delete [] bins;
            delete [] STI;
            delete [] txArray;
            delete [] tyArray;
            delete [] twArray;
            delete [] Y;
            return 1.0;
         }

         for (ss = 0; ss < nSamples; ss+=nReps)
         {
            symIndex = ss / nReps;
            fixedInput  = txArray[ss];
            fixedInput2 = twArray[ss];
            vceMean[symIndex] = 0.0;
            validnReps = 0;
            for (repID = 0; repID < nReps; repID++)
            {
               if (PABS(txArray[ss+repID] - fixedInput) < 1.0E-10 &&
                   PABS(twArray[ss+repID] - fixedInput2) < 1.0E-10)
               {
                  if (tyArray[ss+repID] < 0.9*PSUADE_UNDEFINED)
                  {
                     validnReps++;
                     vceMean[symIndex] += tyArray[ss+repID];
                  }
               }
               else
               {
                  printf("Two-way Effect ERROR: consult developers.\n");
                  exit(1);
               }
            }
            bins[symIndex] = validnReps;
            if (validnReps > 0) vceMean[symIndex] /= (double) validnReps;
            vceVariance[symIndex] = 0.0;
            for (repID = 0; repID < nReps; repID++)
            {
               if (PABS(txArray[ss+repID] - fixedInput) < 1.0E-10 &&
                   PABS(twArray[ss+repID] - fixedInput2) < 1.0E-10)
               {
                  if (tyArray[ss+repID] < 0.9*PSUADE_UNDEFINED)
                  {
                     vceVariance[symIndex] += 
                       ((tyArray[ss+repID]-vceMean[symIndex])*
                        (tyArray[ss+repID]-vceMean[symIndex]));
                  }
               }
            }
            if (validnReps > 0)
                 vceVariance[symIndex] /= ((double) validnReps);
            else vceVariance[symIndex] = PSUADE_UNDEFINED;
         }
         nSubs = nSamples / nReps;

         totalCnt = 0;
         for (subID = 0; subID < nSubs; subID++) totalCnt += bins[subID];

         varVCEMean[ii*nInputs+ii2] = 0.0;
         meanVCEVar[ii*nInputs+ii2] = 0.0;
         aMean = 0.0;
         for (subID = 0; subID < nSubs; subID++)
         {
            if (vceVariance[subID] != PSUADE_UNDEFINED)
            {
               aMean += (vceMean[subID] / totalCnt * bins[subID]);
            }
         }
         for (subID = 0; subID < nSubs; subID++)
         {
            if (vceVariance[subID] != PSUADE_UNDEFINED)
            {
               varVCEMean[ii*nInputs+ii2] += ((vceMean[subID] - aMean)*
                                              (vceMean[subID] - aMean) *
                                              bins[subID] / totalCnt);
               meanVCEVar[ii*nInputs+ii2] += vceVariance[subID] *
                                             bins[subID] / totalCnt;
            }
         }
         vce[ii*nInputs+ii2] = varVCEMean[ii*nInputs+ii2];

         varVCEVar[ii*nInputs+ii2] = 0.0; 
         for (subID = 0; subID < nSubs; subID++)
            if (vceVariance[subID] != PSUADE_UNDEFINED)
               varVCEVar[ii*nInputs+ii2] += 
                  ((vceVariance[subID]-meanVCEVar[ii*nInputs+ii2])*
                   (vceVariance[subID]-meanVCEVar[ii*nInputs+ii2]) * 
                   bins[subID] / totalCnt);
         if (printLevel > 4 || nSamples > 100000) 
            printf("Interaction %3d %d = %e\n",ii+1,ii2+1,
                   vce[ii*nInputs+ii2]/aVariance);
      }
      if (errflag != 0) break;
   }

   for (ii = 0; ii < nInputs; ii++)
   {
      mainEffects[ii] = analyze1D(nInputs,nOutputs,nSamples,X,Y,
                         ii,whichOutput,nReps1,aMean,aVariance,&STI[ii]);
      printf("Main effect (Inputs %3d) = %12.4e\n", ii+1, 
             mainEffects[ii]/aVariance);
   }

   pData *pPtr = ioPtr->getAuxData();
   pPtr->nDbles_ = nInputs * nInputs;
   pPtr->dbleArray_ = new double[nInputs * nInputs];
   for (ii = 0; ii < nInputs*nInputs; ii++) pPtr->dbleArray_[ii] = vce[ii];
   pPtr->dbleData_ = aVariance;

#if 0
   for (ii = 0; ii < nInputs; ii++)
   {
      for (ii2 = ii+1; ii2 < nInputs; ii2++)
      {
         //varVCEMean[ii*nInputs+ii2] -= 
         //               (mainEffects[ii] + mainEffects[ii2]);
         varVCEMean[ii*nInputs+ii2] -= (STI[ii] + STI[ii2]);
         //if (varVCEMean[ii*nInputs+ii2] < 0.0) 
         //   varVCEMean[ii*nInputs+ii2] = 0.0;
      }
   }
#endif

   if (errflag == 0)
   {
      printEquals(0);
      printf("* First order sensitivity indices\n");
      printDashes(0);
      printf("* The first number is a measure of variances due to\n");
      printf("* varying the corresponding input alone.\n");
      printf("* The measure in parentheis is one minus the mean of\n");
      printf("* variances indicating the contribution by terms with\n");
      printf("* this input alone relative to the interaction terms\n");
      printf("* involving this input with other inputs.\n");
      printDashes(0);
      for (ii = 0; ii < nInputs; ii++)
      {
         if (mainEffects[ii] < 1.0e-10 * aVariance)
            mainEffects[ii] = 0.0;
         if (STI[ii] < 0.0) STI[ii] = 0.0;
         printf("* Sensitivity index (Inputs %3d) = %12.4e (%12.4e)\n",
                ii+1,mainEffects[ii]/aVariance,STI[ii]);
      }
      printEquals(0);
      printf("* (first + second) order sensitivity index\n");
      printDashes(0);
      for (ii = 0; ii < nInputs; ii++)
      {
         for (ii2 = ii+1; ii2 < nInputs; ii2++)
         {
            printf("*   SENSITIVITY INDEX  ( Inputs %2d %2d ) = %12.4e\n",
                   ii+1,ii2+1,vce[ii*nInputs+ii2]/aVariance);
         }
      }
      if (constrPtr == NULL)
      {
         printEquals(0);
         printf("* Just the second order sensitivity index\n");
         printf("* (by subtracting mainEffects i, and just from VCE(i,j))\n");
         printf("* Valid only for orthogonal (uncorrelated) inputs.\n");
         printDashes(0);
         for (ii = 0; ii < nInputs; ii++)
         {
            for (ii2 = ii+1; ii2 < nInputs; ii2++)
            {
               ddata = vce[ii*nInputs+ii2];
               ddata -= (mainEffects[ii] + mainEffects[ii2]);
               ddata /= aVariance;
               printf("*2nd order sensitivity index (Inputs %2d %2d ) = %12.4e\n",
                      ii+1,ii2+1, ddata);
            }
         }
      }
   }

#if 0
   if (errflag == 0)
   {
      printf("*=======================================================**\n");
      printf("* Total variance due to the complementary set           **\n");
      printf("*-------------------------------------------------------**\n");
      for (ii = 0; ii < nInputs; ii++)
      {
         for (ii2 = ii+1; ii2 < nInputs; ii2++)
         {
            printf("*   Mean(Var) ratio   (Inputs %2d,%2d) = %12.4e\n",
                   ii+1,ii2+1, meanVCEVar[ii*nInputs+ii2]/aVariance);
         }
      }
      printf("*=======================================================**\n");
      printf("* Total second order sensitivity index                  **\n");
      printf("*-------------------------------------------------------**\n");
      for (ii = 0; ii < nInputs; ii++)
      {
         for (ii2 = ii+1; ii2 < nInputs; ii2++)
         {
            printf("*   Sensitivity index (Inputs %2d,%2d) = %12.4e\n",
                   ii+1,ii2+1, 1.0-meanVCEVar[ii*nInputs+ii2]/aVariance);
         }
      }
      printf("*=======================================================**\n");
      printf("* Strength of higher order interactions                 **\n");
      printf("*-------------------------------------------------------**\n");
      for (ii = 0; ii < nInputs; ii++)
      {
         for (ii2 = ii+1; ii2 < nInputs; ii2++)
         {
            printf("*   var VCE variance   (Inputs %2d,%2d) = %12.4e\n",
                   ii+1,ii2+1,varVCEVar[ii*nInputs+ii2]);
         }
      }
   }
#endif
   printAsterisks(0);
   printAsterisks(0);

   if (psAnaExpertMode_ == 1)
   {
      printf("Bootstrap analysis takes the sample, replicates it n times,\n");
      printf("and assess whether the sensitivity indices have converged.\n");
      printf("If you are performed iterative analysis with refinements,\n");
      printf("you will need to enter 'no index reuse' below at the first\n");
      printf("iteration and 'yes' afterward until the final refinement.\n");
      sprintf(pString,"Perform bootstrap interaction analysis? (y or n) ");
      getString(pString, winput);

      ncount = 0;
      if (winput[0] == 'y')
      {
         sprintf(pString,"Enter number of bootstrap samples (>=100): ");
         ncount = getInt(100, 2000, pString);
         XX = new double[nInputs*nSamples];
         YY = new double[nSamples];
         bsVCEs = new double*[ncount];
         for (ii = 0; ii < ncount; ii++)
            bsVCEs[ii] = new double[nInputs*nInputs];
         nReps = nSamples / nSubSamples;

         fp1 = fopen(".IE_bootstrap_indset", "r");
         if (fp1 != NULL)
         {
            printf(".IE_bootstrap_indset file found.\n");
            sprintf(pString,"Re-use file? (y or n) ");
            getString(pString, winput);
            if (winput[0] == 'y')
            {
               fscanf(fp1, "%d", &ii);
               if (ii != nReps*ncount)
               {
                  printf("ERROR: expect the first line to be %d.\n",
                         nReps*ncount);
                  printf("       Instead found the first line to be %d\n",
                         ii);
                  exit(1);
               }
            }
            else
            {
               fclose(fp1);
               fp1 = fopen(".IE_bootstrap_indset", "w");
               if (fp1 == NULL)
               {
                  printf("ERROR: cannot open ME_bootstrap_indset file.\n");
                  exit(1);
               }
               fprintf(fp1, "%d\n", nReps*ncount);
            }
         }
         else
         {
            fp1 = fopen(".IE_bootstrap_indset", "w");
            if (fp1 == NULL)
            {
               printf("ERROR: cannot open .IE_bootstrap_indset file.\n");
               exit(1);
            }
            fprintf(fp1, "%d\n", nReps*ncount);
         }
         fp = NULL;

         for (ii = 0; ii < ncount; ii++)
         {
            if (ii % (ncount / 10) == 0)
               printf("Processing %d(a) of %d\n", ii+1, ncount);
            for (ir = 0; ir < nReps; ir++)
            {
               if (fp1 != NULL && winput[0] == 'y')
                  fscanf(fp1, "%d\n", &index);
               else
                  index = PSUADE_rand() % nReps;
               if (fp1 != NULL && winput[0] != 'y')
                  fprintf(fp1, "%d\n", index);

               for (ss = 0; ss < nSubSamples; ss++)
               {
                  for (ii2 = 0; ii2 < nInputs; ii2++)
                     XX[(ir*nSubSamples+ss)*nInputs+ii2] =
                        XIn[(index*nSubSamples+ss)*nInputs+ii2];
                  YY[ir*nSubSamples+ss] = 
                     YIn[(index*nSubSamples+ss)*nOutputs+whichOutput];
               }
            }
            fp = NULL;
            if (ii % (ncount / 10) == 0)
               printf("Processing %d(b) of %d\n", ii+1, ncount);
            computeMeanVariance(nInputs,iOne,nSamples,YY,&aMean,&aVariance,
                                iZero);
            if (PABS(aVariance) < 1.0e-15)
            {
               printf("INFO: variance too small ==> terminate.\n");
               continue;
            }
            computeVCE2(nInputs, nSamples, nSubSamples, XX, YY, fp, aMean, 
                        varVCEMean, meanVCEVar, varVCEVar, vce);
            for (ii2 = 0; ii2 < nInputs; ii2++)
               mainEffects[ii2] = analyze1D(nInputs,iOne,nSamples,XX,YY,
                                   ii2,iZero,nReps1,aMean,aVariance,&STI[ii2]);
            for (ii2 = 0; ii2 < nInputs; ii2++)
               for (ss = ii2+1; ss < nInputs; ss++)
                  bsVCEs[ii][ii2*nInputs+ss] = (vce[ii2*nInputs+ss] - 
                          mainEffects[ii2] - mainEffects[ss]) / aVariance;
         }
         if (fp1 != NULL) fclose(fp1);
         if (psPlotTool_ == 0)
         {
            printf("Enter name of matlab/scilab file to store bootstrap info: ");
            // winput is a char array of size 501 so just for defensive programming 
            // Bill Oliver added a width specifier
            scanf("%500s", winput);
            fgets(pString,500,stdin);
            fp = fopen(winput, "w");
         }
         else fp = NULL;
         if (fp != NULL)
         {
            if (psPlotTool_ == 1)
            {
               fprintf(fp, "//bootstrap sample of interaction effects\n");
               fprintf(fp, "//VCE((i-1)*n+j,k) = VCE(i,j), bs sample k\n");
            }
            else
            {
               fprintf(fp, "%%bootstrap sample of interaction effects\n");
               fprintf(fp, "%%VCE((i-1)*n+j,k) = VCE(i,j), bs sample k\n");
            }
            fprintf(fp, "VCE = zeros(%d,%d);\n",nInputs*nInputs,ncount);
            for (ss = 0; ss < ncount; ss++)
            {
               for (ii = 0; ii < nInputs; ii++)
               {
                  for (ii2 = ii+1; ii2 < nInputs; ii2++)
                  {
                     if (ss == 0)
                     {
                        if (psPlotTool_ == 1)
                             fprintf(fp,"//Interaction(%d,%d) \n",ii+1,ii2+1);
                        else fprintf(fp,"%%Interaction(%d,%d) \n",ii+1,ii2+1);
                     }
                     fprintf(fp, "VCE(%d,%d) = %e;\n",ii*nInputs+ii2+1,ss+1,
                             bsVCEs[ss][ii*nInputs+ii2]);
                  }
               }
            }
            fclose(fp);
         }
         else
         {
            printf("ERROR: cannot open file %s\n", winput);
         }
         for (ii = 0; ii < ncount; ii++) delete [] bsVCEs[ii];
         delete [] bsVCEs;
         delete [] XX;
         delete [] YY;
      }
   }

   delete [] vce;
   delete [] varVCEMean;
   delete [] meanVCEVar;
   delete [] varVCEVar;
   delete [] vceMean;
   delete [] vceVariance;
   delete [] mainEffects;
   delete [] STI;
   delete [] txArray;
   delete [] twArray;
   delete [] tyArray;
   delete [] Y;
   delete [] bins;
   if (constrPtr != NULL) delete constrPtr;
   // return 1.0 to facilitate continuous refinement
   return 1.0;
}

// *************************************************************************
// Compute agglomerated mean and variance
// -------------------------------------------------------------------------
int TwoParamAnalyzer::computeMeanVariance(int nInputs, int nOutputs, 
                                     int nSamples, double *y, double *aMean, 
                                     double *aVariance, int outputID)
{
   int    ss, count;
   double mean, variance;

   mean = 0.0;
   count = 0;
   for (ss = 0; ss < nSamples; ss++)
   {
      if (y[nOutputs*ss+outputID] < 0.9*PSUADE_UNDEFINED)
      {
         mean += y[nOutputs*ss+outputID];
         count++;
      }

   }
   if (count <= 0)
   {
      printf("Two-way Effect ERROR: no valid data.\n");
      exit(1);
   }
   mean /= (double) count;
   variance = 0.0;
   for (ss = 0; ss < nSamples; ss++) 
   {
      if (y[nOutputs*ss+outputID] < 0.9*PSUADE_UNDEFINED)
      {
         variance += ((y[nOutputs*ss+outputID] - mean) *
                      (y[nOutputs*ss+outputID] - mean));
      }
   }
   variance /= (double) count;
   (*aMean) = mean;
   (*aVariance) = variance;
   return 0;
}

// ************************************************************************
// perform VCE (McKay) analysis
// ------------------------------------------------------------------------
double TwoParamAnalyzer::analyze1D(int nInputs, int nOutputs, int nSamples,
                           double *x, double *y, int inputID, int outputID,
                           int nReps, double aMean, double aVariance,
                           double *retData)
{
   int    ss, repID, subID, ncount, symIndex, nSubs;
   int    whichOutput, validnReps, totalCnt, *bins;
   double meanVCEVar, vce, *vceMean, *vceVariance, varVCEMean;
   double fixedInput, *txArray, *tyArray, tmean;

   whichOutput = outputID;
   if (whichOutput >= nOutputs || whichOutput < 0) whichOutput = 0;

   if (nReps <= 1) return 1.0;
   nSubs         = nSamples / nReps;
   vceMean       = new double[nSubs];
   vceVariance   = new double[nSubs];
   txArray       = new double[nSamples];
   tyArray       = new double[nSamples];
   bins          = new int[nSubs];

   for (ss = 0; ss < nSamples; ss++)
   {
      txArray[ss] = x[nInputs*ss+inputID];
      tyArray[ss] = y[nOutputs*ss+whichOutput];
   }
   sortDbleList2(nSamples, txArray, tyArray);
   ncount = 0;
   for (ss = 0; ss < nSamples; ss+=nReps)
   {
      symIndex = ss / nReps;
      fixedInput = txArray[ss];
      vceMean[symIndex] = 0.0;
      validnReps = 0;
      for (repID = 0; repID < nReps; repID++)
      {
         if (PABS(txArray[ss+repID] - fixedInput) < 1.0E-10)
         {
            ncount++;
            if (tyArray[ss+repID] < 0.9*PSUADE_UNDEFINED)
            {
               vceMean[symIndex] += tyArray[ss+repID];
               validnReps++;
            }
         }
      }
      if (validnReps > 0) vceMean[symIndex] /= (double) validnReps;
      vceVariance[symIndex] = 0.0;
      tmean = vceMean[symIndex];
      for (repID = 0; repID < nReps; repID++)
      {
         if (PABS(txArray[ss+repID] - fixedInput) < 1.0E-10)
         {
            if (tyArray[ss+repID] < 0.9*PSUADE_UNDEFINED)
               vceVariance[symIndex] += ((tyArray[ss+repID]-tmean)*
                                         (tyArray[ss+repID]-tmean));
         }
      }
      if (validnReps > 0) vceVariance[symIndex] /= validnReps;
      else                vceVariance[symIndex] = PSUADE_UNDEFINED;
      bins[symIndex] = validnReps;
   }
   if (ncount != nSamples)
   {
      printf("Two-way Effect ERROR: sample not valid, input %d\n",
             inputID);
      printf("       error data = %d (%d)\n",ncount, nSamples);
      //for (ss = 0; ss < nSamples; ss++)
      //   printf("%5d : %16.8e \n", ss, txArray[ss]);
      exit(1);
   }

   varVCEMean = 0.0;
   meanVCEVar = 0.0;
   totalCnt = 0;
   for (subID = 0; subID < nSubs; subID++) totalCnt += bins[subID];
   aMean = 0.0;
   for (subID = 0; subID < nSubs; subID++)
   {
      if (vceVariance[subID] != PSUADE_UNDEFINED)
         aMean += (vceMean[subID] / totalCnt * bins[subID]);
   }
   for (subID = 0; subID < nSubs; subID++)
   {
      if (vceVariance[subID] != PSUADE_UNDEFINED)
      {
         varVCEMean += ((vceMean[subID]-aMean) * (vceMean[subID]-aMean) *
                        bins[subID] / totalCnt);
         meanVCEVar += vceVariance[subID] * bins[subID] / totalCnt;
      }
   }
   vce = varVCEMean;

   delete [] txArray;
   delete [] tyArray;
   delete [] vceMean;
   delete [] vceVariance;
   delete [] bins;
   (*retData) = 1.0 - meanVCEVar / aVariance;
   return vce;
}

// *************************************************************************
// Compute two-way VCE 
// -------------------------------------------------------------------------
int TwoParamAnalyzer::computeVCE2(int nInputs, int nSamples, int nSubSamples,
                                 double *X, double *Y, FILE *fp, 
                                 double aMean, double *varVCEMean,
                                 double *meanVCEVar, double *varVCEVar,
                                 double *vce)
{
   int    ii, ii2, ss, nReps1, nReps, nSubs, symIndex, repID, subID;
   int    totalCnt, *bins, validnReps;
   double *txArray, *twArray, *tyArray, fixedInput, fixedInput2, *vceMean;
   double *vceVariance;


   txArray = new double[nSamples];
   twArray = new double[nSamples];
   tyArray = new double[nSamples];
   vceMean     = new double[nSubSamples];
   vceVariance = new double[nSubSamples];
   bins        = new int[nSubSamples];

   for (ii = 0; ii < nInputs; ii++)
   {
      for (ii2 = ii+1; ii2 < nInputs; ii2++)
      {
         for (ss = 0; ss < nSamples; ss++)
         {
            txArray[ss] = X[nInputs*ss+ii];
            twArray[ss] = X[nInputs*ss+ii2];
            tyArray[ss] = Y[ss];
         }

         sortDbleList3(nSamples, txArray, twArray, tyArray);

         for (ss = 1; ss < nSamples; ss++)
            if (PABS(txArray[ss]-txArray[ss-1]) > 1.0E-10)
               break;
         nReps1 = ss;

         for (ss = 0; ss < nSamples; ss+=nReps1)
            sortDbleList2(nReps1,&twArray[ss],&tyArray[ss]);

         for (ss = 1; ss < nReps1; ss++)
            if (PABS(twArray[ss]-twArray[ss-1]) > 1.0E-10)
               break;
         nReps = ss;

         for (ss = 0; ss < nSamples; ss+=nReps)
         {
            symIndex = ss / nReps;
            fixedInput  = txArray[ss];
            fixedInput2 = twArray[ss];
            vceMean[symIndex] = 0.0;
            validnReps = 0;
            for (repID = 0; repID < nReps; repID++)
            {
               if (PABS(txArray[ss+repID] - fixedInput) < 1.0E-10 &&
                   PABS(twArray[ss+repID] - fixedInput2) < 1.0E-10)
               {
                  if (tyArray[ss+repID] < 0.9*PSUADE_UNDEFINED)
                  {
                     vceMean[symIndex] += tyArray[ss+repID];
                     validnReps++;
                  }
               }
            }
            if (validnReps > 0) vceMean[symIndex] /= (double) validnReps;
            vceVariance[symIndex] = 0.0;
            for (repID = 0; repID < nReps; repID++)
            {
               if (tyArray[ss+repID] < 0.9*PSUADE_UNDEFINED)
                  vceVariance[symIndex] += 
                       ((tyArray[ss+repID]-vceMean[symIndex])*
                        (tyArray[ss+repID]-vceMean[symIndex]));
            }
            if (validnReps > 0) 
                 vceVariance[symIndex] /= ((double) validnReps);
            else vceVariance[symIndex] = PSUADE_UNDEFINED;
         }
         nSubs = nSamples / nReps;

         varVCEMean[ii*nInputs+ii2] = 0.0;
         meanVCEVar[ii*nInputs+ii2] = 0.0;
         totalCnt = 0;
         for (subID = 0; subID < nSubs; subID++) totalCnt += bins[subID];
         for (subID = 0; subID < nSubs; subID++)
         {
            if (vceVariance[subID] != PSUADE_UNDEFINED)
            {
               varVCEMean[ii*nInputs+ii2] += ((vceMean[subID] - aMean)*
                                              (vceMean[subID] - aMean) *
                                              bins[subID] / totalCnt);
               meanVCEVar[ii*nInputs+ii2] += vceVariance[subID] * bins[subID] /
                                             totalCnt;
            }
         }
         vce[ii*nInputs+ii2] = varVCEMean[ii*nInputs+ii2];

         varVCEVar[ii*nInputs+ii2] = 0.0; 
         for (subID = 0; subID < nSubs; subID++)
            if (vceVariance[subID] != PSUADE_UNDEFINED)
               varVCEVar[ii*nInputs+ii2] += 
                  ((vceVariance[subID]-meanVCEVar[ii*nInputs+ii2])*
                   (vceVariance[subID]-meanVCEVar[ii*nInputs+ii2]) *
                   bins[subID] / totalCnt);
      }
   }
   delete [] txArray;
   delete [] twArray;
   delete [] tyArray;
   delete [] vceMean;
   delete [] bins;
   delete [] vceVariance;
   return 0;
}

// *************************************************************************
// Compute VCE 
// -------------------------------------------------------------------------
int TwoParamAnalyzer::computeVCECrude(int nInputs, int nSamples, 
                             double *X, double *Y,double aVariance,
                             double *iLowerB, double *iUpperB, double *vce)
{
   int    ii, ii2, ss, nSize, nIntervals, ind1, ind2, ind;
   int    totalCnt, *bins, *tags, nIntervals1, nSize1;
   double *vceMean, *vceVariance, ddata, hstep1, hstep2, aMean, *vce1;
   char   pString[500];

   nSize = 10;
   ddata = 1.0 * nSamples / nSize;
   ddata = pow(ddata, 0.5);
   nIntervals = (int) ddata;
   nIntervals1 = nIntervals * nIntervals;
   if (nIntervals < 4)
   {
      printf("ERROR: sample size too small.\n");
      printf("       Need sample size >= 160 to give meaningful results.\n");
      return -1;
   }
   printAsterisks(0);
   printf("           Crude 2-way Interaction Effect\n");
   printEquals(0);
   printf("Two-way Effect: default number of levels in each input  = %d\n", 
          nIntervals);
   printf("Two-way Effect: default number of level for main effect = %d\n", 
          nIntervals1);
   printf("Two-way Effect: sample size for each 2-dimensional box  = %d\n",
          nSize);
   printDashes(0);
   printf("* These may need to be adjusted for higher accuracy.\n");
   printf("* Recommend sample size per 2D box to be at least 10 to\n");
   printf("* ensure sufficiency accuracy.\n");
   printf("* Turn on analysis expert mode to try different number of levels.\n");
   if (psAnaExpertMode_ == 1)
   {
      printEquals(0);
      sprintf(pString,"Number of levels for 2-way analysis (>= 4, default = %d): ",
              nIntervals);
      nIntervals = getInt(4, nSamples, pString);
      nSize = nSamples / nIntervals / nIntervals;
      if (nSize < 10)
      {
         printf("* Using number of levels = %d for 2-way effect will give\n", 
                nIntervals);
         printf("  each 2D box less than 10 sample points. That is too small.\n");
         nSize = 10;
         ddata = 1.0 * nSamples / nSize;
         ddata = pow(ddata, 0.5);
         nIntervals = (int) ddata;
         printf("* Number of levels for 2-way analysis defaulted to %d\n", nIntervals);
      } 
      sprintf(pString,"Number of levels for main effect (>= %d): ",nIntervals);
      nIntervals1 = getInt(nIntervals, nSamples, pString);
      nSize1 = nSamples / nIntervals1;
      if (nSize1 < 10)
      {
         printf("Sample size per level for main effect %d too small (should be >= 10).\n",
                nSize1);
         nSize1 = 10;
         nIntervals1 = nSamples / nSize1;
         printf("Default number of levels for main effect to %d\n", nIntervals1);
      } 
   }
   printEquals(0);
   nSize = nIntervals * nIntervals;
   if (nSize < nIntervals1) nSize = nIntervals1;
   vceMean     = new double[nSize];
   vceVariance = new double[nSize];
   bins        = new int[nSize];
   tags        = new int[nSamples];
   vce1        = new double[nInputs];

   for (ii = 0; ii < nInputs; ii++)
   {
      hstep1 = (iUpperB[ii] - iLowerB[ii]) / nIntervals1;
      for (ss = 0; ss < nIntervals1; ss++)
      {
         vceMean[ss] = 0.0;
         vceVariance[ss] = 0.0;
         bins[ss] = 0;
         tags[ss] = 0;
      }
      for (ss = 0; ss < nSamples; ss++)
      {
         ddata = (X[nInputs*ss+ii] - iLowerB[ii]) / hstep1;
         ind1 = (int) ddata;
         if (ind1 >= nIntervals1) ind1 = nIntervals1 - 1;
         tags[ind1] = -1;
         if (Y[ss] < 0.9*PSUADE_UNDEFINED)
         {
            vceMean[ind1] += Y[ss];
            bins[ind1]++;
            tags[ss] = ind1;
         }
      }
      for (ss = 0; ss < nIntervals1; ss++)
         if (bins[ss] > 0) vceMean[ss] /= (double) bins[ss];

      for (ss = 0; ss < nSamples; ss++)
      {
         ind = tags[ss];
         if (ind >= 0 && Y[ss] < 0.9*PSUADE_UNDEFINED)
         {
            vceVariance[ind] += ((Y[ss]-vceMean[ind])*
                                 (Y[ss]-vceMean[ind]));
         }
      }
      for (ss = 0; ss < nIntervals1; ss++)
      {
         if (bins[ss] > 0)
            vceVariance[ss] /= (double) bins[ss];
         else vceVariance[ss] = PSUADE_UNDEFINED;
      }

      totalCnt = 0;
      for (ss = 0; ss < nIntervals1; ss++) totalCnt += bins[ss];
      aMean = 0.0;
      for (ss = 0; ss < nIntervals1; ss++)
      {
         if (vceVariance[ss] != PSUADE_UNDEFINED)
            aMean += (vceMean[ss] / totalCnt * bins[ss]);
      }
      vce1[ii] = 0.0;
      for (ss = 0; ss < nIntervals1; ss++)
      {
         if (vceVariance[ss] != PSUADE_UNDEFINED)
         {
            vce1[ii] += ((vceMean[ss] - aMean) *
                         (vceMean[ss] - aMean) * bins[ss]/totalCnt);
         }
      }
      printf("*   SENSITIVITY INDEX  ( Inputs %2d ) = %12.4e",
             ii+1,vce1[ii]/aVariance);
      printf(" (unnormalized = %12.4e)\n",vce1[ii]);
   }
   for (ii = 0; ii < nInputs; ii++)
   {
      for (ii2 = ii+1; ii2 < nInputs; ii2++)
      {
         hstep1 = (iUpperB[ii] - iLowerB[ii]) / nIntervals;
         hstep2 = (iUpperB[ii2] - iLowerB[ii2]) / nIntervals;
         for (ss = 0; ss < nIntervals*nIntervals; ss++)
         {
            vceMean[ss] = 0.0;
            vceVariance[ss] = 0.0;
            bins[ss] = 0;
            tags[ss] = 0;
         }
         for (ss = 0; ss < nSamples; ss++)
         {
            ddata = (X[nInputs*ss+ii] - iLowerB[ii]) / hstep1;
            ind1 = (int) ddata;
            if (ind1 >= nIntervals) ind1 = nIntervals - 1;
            ddata = (X[nInputs*ss+ii2] - iLowerB[ii2]) / hstep2;
            ind2 = (int) ddata;
            if (ind2 >= nIntervals) ind2 = nIntervals - 1;
            ind = ind1 * nIntervals + ind2;
            tags[ind] = -1;
            if (Y[ss] < 0.9*PSUADE_UNDEFINED)
            {
               vceMean[ind] += Y[ss];
               bins[ind]++;
               tags[ss] = ind;
            }
         }

         for (ss = 0; ss < nIntervals*nIntervals; ss++)
            if (bins[ss] > 0) vceMean[ss] /= (double) bins[ss];

         for (ss = 0; ss < nSamples; ss++)
         {
            ind = tags[ss];
            if (ind >= 0 && Y[ss] < 0.9*PSUADE_UNDEFINED)
            {
               vceVariance[ind] += ((Y[ss]-vceMean[ind])*
                                    (Y[ss]-vceMean[ind]));
            }
         }
         for (ss = 0; ss < nIntervals*nIntervals; ss++)
         {
            if (bins[ss] > 0)
                 vceVariance[ss] /= (double) bins[ss];
            else vceVariance[ss] = PSUADE_UNDEFINED;
         }

         totalCnt = 0;
         for (ss = 0; ss < nIntervals*nIntervals; ss++)
            totalCnt += bins[ss];
         aMean = 0.0;
         for (ss = 0; ss < nIntervals*nIntervals; ss++)
         {
            if (vceVariance[ss] != PSUADE_UNDEFINED)
               aMean += (vceMean[ss] / totalCnt * bins[ss]);
         }
         vce[ii*nInputs+ii2] = 0.0;
         for (ss = 0; ss < nIntervals*nIntervals; ss++)
         {
            if (vceVariance[ss] != PSUADE_UNDEFINED)
            {
               vce[ii*nInputs+ii2] += ((vceMean[ss] - aMean) *
                             (vceMean[ss] - aMean) * bins[ss]/totalCnt);
            }
         }
         if ( vce[ii*nInputs+ii2] < 0.0) vce[ii*nInputs+ii2] = 0.0;
         printf("*   SENSITIVITY INDEX  (Inputs %2d %2d) = %12.4e",
                ii+1,ii2+1,vce[ii*nInputs+ii2]/aVariance);
         printf(" (unnormalized = %12.4e) : 1st + 2nd order\n",
                vce[ii*nInputs+ii2]);
      }
   }
   printAsterisks(0);
   delete [] vceMean;
   delete [] vceVariance;
   delete [] bins;
   delete [] tags;
   delete [] vce1;
   return 0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
TwoParamAnalyzer& TwoParamAnalyzer::operator=(const TwoParamAnalyzer &)
{
   printf("Two-way Effect operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

