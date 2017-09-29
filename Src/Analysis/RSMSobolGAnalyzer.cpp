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
// Functions for the class RSMSobolGAnalyzer  
// Sobol' group main effect analysis (with response surface)
// AUTHOR : CHARLES TONG
// DATE   : 2007
// ************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

#include "Util/PsuadeUtil.h"
#include "Util/sysdef.h"
#include "Util/Vector.h"
#include "Util/Matrix.h"
#include "Main/Psuade.h"
#include "FuncApprox/FuncApprox.h"
#include "Samplings/RSConstraints.h"
#include "Samplings/Sampling.h"
#include "PDFLib/PDFManager.h"
#include "DataIO/PsuadeData.h"
#include "DataIO/pData.h"
#include "RSMSobolGAnalyzer.h"

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
RSMSobolGAnalyzer::RSMSobolGAnalyzer() : Analyzer()
{
   setName("RSMSOBOLG");
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
RSMSobolGAnalyzer::~RSMSobolGAnalyzer()
{
}

// ************************************************************************
// perform analysis
// ------------------------------------------------------------------------
double RSMSobolGAnalyzer::analyze(aData &adata)
{
   int        nInputs, nOutputs, nSamples, *S, fileValid, ii, ii2, jj, kk;
   int        nGroups, groupID, **groupMembers, length, status, sCnt;
   int        ir, index, currNSamples, nSubSamplesG, nInputsG, nInputsN;
   int        nSubSamplesN, nSamp, *pdfFlags, outputID, noPDF=1, *SS;
   int        printLevel, *bins, totalCnt, *iArray, *pdfFlagsG, *pdfFlagsN;
   size_t     compFlag;
   double     *xLower, *xUpper, *X, *Y, *Y2, *XX, *YY, *XXG, *XXN, ddata;
   double     *oneSamplePt, variance, *cLower, *cUpper, *means, *vars;
   double     *inputMeans, *inputStdevs, vce, dmean, ecv, *mSamplePts;
   double     *inputMeansG, *inputMeansN, *inputStdevsG, *inputStdevsN;
   PsuadeData *ioPtr;
   ifstream   ifile;
   char       cfname[255], pString[500];
   string     sfname, iline;
   RSConstraints *constrPtr;
   FuncApprox    *faPtr;
   Vector        vecIn, vecOut, vecUB, vecLB;
   pData         pCorMat;
   Matrix        *corMatp, corMat;
   PDFManager    *pdfman, *pdfmanN, *pdfmanG;
   Sampling      *sampler;

   printLevel = adata.printLevel_;
   nInputs  = adata.nInputs_;
   nOutputs = adata.nOutputs_;
   nSamples = adata.nSamples_;
   xLower   = adata.iLowerB_;
   xUpper   = adata.iUpperB_;
   X        = adata.sampleInputs_;
   Y2       = adata.sampleOutputs_;
   S        = adata.sampleStates_;
   outputID = adata.outputID_;
   ioPtr    = adata.ioPtr_;
   pdfFlags    = adata.inputPDFs_;
   inputMeans  = adata.inputMeans_;
   inputStdevs = adata.inputStdevs_;
   if (pdfFlags != NULL)
   {
      for (ii = 0; ii < nInputs; ii++)
         if (pdfFlags[ii] != 0) noPDF = 0;
   }
   if (noPDF == 1) printf("RSMSobolG INFO: all uniform distributions.\n");
   else
   {
      printf("RSMSobolG INFO: non-uniform distributions detected which\n");
      printf("                will be used in this analysis.\n");
   }

   if (nInputs <= 1 || nSamples <= 0 || nOutputs <= 0)
   {
      printf("RSMSobolG ERROR: invalid arguments.\n");
      printf("   nInputs  = %d\n", nInputs);
      printf("   nOutputs = %d\n", nOutputs);
      printf("   nSamples = %d\n", nSamples);
      return PSUADE_UNDEFINED;
   } 
   if (nInputs <= 2)
   {
      printf("RSMSobolG INFO: nInputs == 2.\n");
      printf("   You do not need to perform this when nInputs = 2.\n");
      return PSUADE_UNDEFINED;
   }
   if (outputID >= nOutputs || outputID < 0)
   {
      printf("RSMSobolG ERROR: invalid output ID (%d).\n", outputID);
      return PSUADE_UNDEFINED;
   }
   if (ioPtr == NULL)
   {
      printf("RSMSobolG ERROR: no data.\n");
      return PSUADE_UNDEFINED;
   } 
   Y = new double[nSamples];
   for (ii = 0; ii < nSamples; ii++) Y[ii] = Y2[ii*nOutputs+outputID];

   printAsterisks(0);
   printf("To use this function, you need to provide a file\n");
   printf("specifying group information, in the form of : \n");
   printf("line 1: PSUADE_BEGIN\n");
   printf("line 2: <d> specifying the number of groups\n");
   printf("line 3 to line <d>+2: group number, size, input numbers\n"); 
   printf("last line: PSUADE_END\n");
   fileValid = 500;
   while (fileValid > 250)
   {
      printf("Enter the group file : ");
      cin >> sfname;
      fileValid = sfname.size();

      if (fileValid < 250)
      {
         sfname.copy(cfname, fileValid, 0);
         cfname[fileValid] = '\0';
         ifile.open(cfname);
         if (ifile.is_open()) break;
         fileValid = 500;
      }
      else
      {
         printf("Invalid file name (either nonexistent or name too long).\n");
      }
   }
   if (ifile.is_open())
   {
      getline (ifile, iline);
      compFlag = iline.compare("PSUADE_BEGIN");
      if (compFlag == 0)
      {
         ifile >> nGroups;
         if (nGroups <= 0)
         {
            printf("RSMSobolG ERROR: nGroups <= 0.\n");
            exit(1);
         }
         groupMembers = new int*[nGroups];
         for (ii = 0; ii < nGroups; ii++)
         {
            ifile >> groupID;
            if (groupID != ii+1)
            {
               printf("RSMSobolG ERROR: invalid groupID %d",groupID);
               printf(" should be %d\n", ii+1);
               exit(1);
            }
            ifile >> length;
            if (length <= 0 || length >= nInputs)
            {
               printf("RSMSobolG ERROR: invalid group length.\n");
               exit(1);
            }
            groupMembers[ii] = new int[nInputs];
            for (jj = 0; jj < nInputs; jj++) groupMembers[ii][jj] = 0;
            sCnt = 1;
            for (jj = 0; jj < length; jj++)
            {
               ifile >> index;
               if (index <= 0 || index > nInputs)
               {
                  printf("RSMSobolG ERROR: invalid group member.\n");
                  exit(1);
               }
               groupMembers[ii][index-1] = sCnt++;
            }
         }
         getline(ifile, iline);
         getline(ifile, iline);
         compFlag = iline.compare("PSUADE_END");
         if (compFlag != 0)
         {
            printf("RSMSobolG ERROR: PSUADE_END not found.\n");
            exit(1);
         }
      }
      else
      {
         printf("RSMSobolG ERROR: PSUADE_BEGIN not found.\n");
         exit(1);
      }
      ifile.close();
   }

   ioPtr->getParameter("input_cor_matrix", pCorMat);
   corMatp = (Matrix *) pCorMat.psObject_;
   for (ii = 0; ii < nGroups; ii++)
   {
      for (ii2 = 0; ii2 < nInputs; ii2++)
      {
         for (jj = 0; jj < nInputs; jj++)
         {
            if (((groupMembers[ii][ii2] != 0 && groupMembers[ii][jj] == 0) ||
                 (groupMembers[ii][ii2] == 0 && groupMembers[ii][jj] != 0)) &&
                corMatp->getEntry(ii2,jj) != 0.0)
            {
               printf("RSMSobolG INFO: this method cannot handle correlated\n");
               printf("         inputs (joint PDF) across different groups.\n");
               for (ir = 0; ir < nGroups; ir++) delete [] groupMembers[ir];
               delete [] groupMembers;
               return PSUADE_UNDEFINED;
           }
         }
      }
   }

   constrPtr = new RSConstraints();
   constrPtr->genConstraints(ioPtr);

   faPtr = genFAInteractive(ioPtr, 0);
   length = -999;
   status = faPtr->genNDGridData(X, Y, &length, NULL, NULL);

   printAsterisks(0);
   if (psAnaExpertMode_ == 1)
   {
      printf("* RSMSobolGAnalyzer creates a sample of size M1 for each\n");
      printf("* subgroup of inputs and M2 for the rest of the inputs. The\n");
      printf("* total sample size is thus:\n");
      printf("* N = M1 * M2 * nGroups.\n");
      printf("* As a user, please decide on M1 and M2.\n");
      printEquals(0);

      sprintf(pString,"Enter M1 (suggestion: 100 - 1000) : ");
      nSubSamplesG = getInt(100, 100000, pString);
      if (nSubSamplesG > 10000)
         printf("An nSubSamples of %d may take very long time.\n",nSubSamplesG);

      sprintf(pString, "Enter M2 (suggestion: > 1000) : ");
      nSubSamplesN = getInt(1000, 100000, pString);
      if (nSubSamplesN > 10000)
      {
         printf("An nLevels of %d may take very long time.\n",nSubSamplesN);
      }
      printAsterisks(0);
   }
   else
   {
      nSubSamplesG = 500;
      nSubSamplesN = 2000;
      printf("RSMSobolG: default M1 = %d.\n", nSubSamplesG);
      printf("RSMSobolG: default M2 = %d.\n", nSubSamplesN);
      printf("To change these settings, turn on ana_expert mode and rerun.\n");
   }
   printEquals(0);

   nSamp = 50000;
   if (printLevel > 1)
   {
      printf("RSMSobolG INFO: creating a sample for basic statistics.\n");
      printf("                sample size = %d\n", nSamp);
   }

   XX = new double[nSamp*8*nInputs];
   YY = new double[nSamp*8];

   for (ir = 0; ir < 3; ir++)
   {
      printf("RSMSobolG: compute statistics, sample size = %d\n",nSamp);
       
      if (noPDF == 0)
      {
         pdfman = new PDFManager();
         pdfman->initialize(nInputs,pdfFlags,inputMeans,
                            inputStdevs,*corMatp);
         vecLB.load(nInputs, xLower);
         vecUB.load(nInputs, xUpper);
         vecOut.setLength(nSamp*nInputs);
         pdfman->genSample(nSamp, vecOut, vecLB, vecUB);
         for (ii = 0; ii < nSamp*nInputs; ii++) XX[ii] = vecOut[ii];
         delete pdfman;
      }
      else
      {
         if (nInputs > 51)
              sampler = SamplingCreateFromID(PSUADE_SAMP_LHS);
         else sampler = SamplingCreateFromID(PSUADE_SAMP_LPTAU);
         sampler->setInputBounds(nInputs, xLower, xUpper);
         sampler->setOutputParams(1);
         sampler->setSamplingParams(nSamp, 1, 1);
         sampler->initialize(0);
         SS = new int[nSamp];
         sampler->getSamples(nSamp, nInputs, 1, XX, YY, SS);
         delete [] SS;
         delete sampler;
      }

      if (ir < 3 && printLevel > 1)
      {
         printf("RSMSobolG: running the sample with response surface...\n");
      }
      faPtr->evaluatePoint(nSamp, XX, YY);
      if (ir < 3 && printLevel > 1)
      {
         printf("RSMSobolG: done running the sample with response surface.\n");
      }

      for (ii = 0; ii < nSamp; ii++)
      {
         oneSamplePt = &(XX[ii*nInputs]);
         ddata = constrPtr->evaluate(oneSamplePt,YY[ii],status);
         if (status == 0) YY[ii] = PSUADE_UNDEFINED;
      }

      dmean = 0.0;
      sCnt = 0;
      for (ii = 0; ii < nSamp; ii++)
      {
         if (YY[ii] != PSUADE_UNDEFINED)
         {
            dmean += YY[ii];
            sCnt++;
         }
      }
      if (sCnt > 1) dmean /= (double) sCnt;
      else
      {
         printf("RSMSobolG ERROR: too few samples that satisify\n"); 
         printf("the constraints (%d out of %d)\n", sCnt, nSamp);
         delete [] XX;
         delete [] YY;
         delete faPtr;
         return PSUADE_UNDEFINED;
      }
      variance = 0.0;
      for (ii = 0; ii < nSamp; ii++)
      {
         if (YY[ii] != PSUADE_UNDEFINED)
            variance += (YY[ii] - dmean) * (YY[ii] - dmean) ;
      }
      variance /= (double) sCnt;
      if (printLevel > 3 || ir == 2)
      {
         printf("RSMSobolG: sample mean    (based on N = %d) = %10.3e\n",
                sCnt, dmean);
         printf("RSMSobolG: sample std dev (based on N = %d) = %10.3e\n",
                sCnt, sqrt(variance));
      }
      nSamp *= 2;
   }
   if (variance == 0.0) variance = 1.0;
   delete [] XX;
   delete [] YY;

   cLower = new double[nInputs];
   cUpper = new double[nInputs];
   XXG    = new double[nSubSamplesG*nInputs];
   XXN    = new double[nSubSamplesN*nInputs];
   YY     = new double[nSubSamplesN+nSubSamplesG];
   means  = new double[nSubSamplesG];
   vars   = new double[nSubSamplesG];
   bins   = new int[nSubSamplesG];
   iArray = new int[nInputs];
   mSamplePts   = new double[nInputs*nSubSamplesN];
   inputMeansG  = new double[nInputs];
   inputMeansN  = new double[nInputs];
   inputStdevsG = new double[nInputs];
   inputStdevsN = new double[nInputs];
   pdfFlagsG    = new int[nInputs];
   pdfFlagsN    = new int[nInputs];

   printAsterisks(0);
   for (ii = 0; ii < nGroups; ii++)
   {
      if (printLevel > 1)
      {
         printf("RSMSobolG: processing group %d\n", ii+1);
         printf("           group members: ");
         for (jj = 0; jj < nInputs; jj++)
         {
            if (groupMembers[ii][jj] != 0)
               printf("%d ", jj+1);
         }
         printf("\n");
      }
      nInputsN = 0;
      for (jj = 0; jj < nInputs; jj++)
      {
         if (groupMembers[ii][jj] == 0)
         {
            cLower[nInputsN] = xLower[jj];
            cUpper[nInputsN] = xUpper[jj];
            iArray[nInputsN] = jj;
            nInputsN++;
         }
      }
      if (noPDF == 0) 
      {
         nInputsN = 0;
         for (jj = 0; jj < nInputs; jj++)
         {
            if (groupMembers[ii][jj] == 0)
            {
               pdfFlagsN[nInputsN] = pdfFlags[jj];
               inputMeansN[nInputsN] = inputMeans[jj];
               inputStdevsN[nInputsN] = inputStdevs[jj];
               nInputsN++;
            }
         }
         corMat.setDim(nInputsN, nInputsN);
         for (jj = 0; jj < nInputsN; jj++)
            for (ii2 = 0; ii2 < nInputsN; ii2++)
               corMat.setEntry(jj,ii2,corMatp->getEntry(iArray[jj],
                                                        iArray[ii2])); 
         pdfmanN = new PDFManager();
         pdfmanN->initialize(nInputsN, pdfFlagsN, inputMeansN,
                             inputStdevsN, corMat);
         vecLB.load(nInputsN, cLower);
         vecUB.load(nInputsN, cUpper);
         vecOut.setLength(nSubSamplesN*nInputsN);
         pdfmanN->genSample(nSubSamplesN, vecOut, vecLB, vecUB);
         for (jj = 0; jj < nSubSamplesN*nInputsN; jj++)
            XXN[jj] = vecOut[jj];
         delete pdfmanN;
      }
      else
      {
         if (nInputsN > 51)
              sampler = SamplingCreateFromID(PSUADE_SAMP_LHS);
         else sampler = SamplingCreateFromID(PSUADE_SAMP_LPTAU);
         sampler->setInputBounds(nInputsN, cLower, cUpper);
         sampler->setOutputParams(1);
         sampler->setSamplingParams(nSubSamplesN, 1, 1);
         sampler->initialize(0);
         SS = new int[nSubSamplesN];
         sampler->getSamples(nSubSamplesN,nInputsN,1,XXN,YY,SS);
         delete [] SS;
         delete sampler;
      }

      currNSamples = nSubSamplesG / 8;
      nInputsG = 0;
      for (jj = 0; jj < nInputs; jj++)
      {
         if (groupMembers[ii][jj] != 0)
         {
            cLower[nInputsG] = xLower[jj];
            cUpper[nInputsG] = xUpper[jj];
            iArray[nInputsG] = jj;
            nInputsG++;
         }
      }
      if (noPDF == 0) 
      {
         nInputsG = 0;
         for (jj = 0; jj < nInputs; jj++)
         {
            if (groupMembers[ii][jj] != 0)
            {
               pdfFlagsG[nInputsG] = pdfFlags[jj];
               inputMeansG[nInputsG] = inputMeans[jj];
               inputStdevsG[nInputsG] = inputStdevs[jj];
               nInputsG++;
            }
         }
         corMat.setDim(nInputsG, nInputsG);
         for (jj = 0; jj < nInputsG; jj++)
            for (ii2 = 0; ii2 < nInputsG; ii2++)
               corMat.setEntry(jj,ii2,
                     corMatp->getEntry(iArray[jj],iArray[ii2])); 
         vecLB.load(nInputsG, cLower);
         vecUB.load(nInputsG, cUpper);
         vecOut.setLength(nSubSamplesG*nInputsG);
      }

      for (ir = 0; ir < 4; ir++)
      {
         if (printLevel > 2)
             printf("RSMSobolG: processing refinement %d\n",ir+1);
         if (printLevel > 3)
             printf("nSamplesG = %d, nSamplesN = %d\n",currNSamples, 
                    nSubSamplesN);
         if (noPDF == 0) 
         {
            pdfmanG = new PDFManager();
            pdfmanG->initialize(nInputsG, pdfFlagsG, inputMeansG,
                                inputStdevsG, corMat);
            pdfmanG->genSample(nSubSamplesG, vecOut, vecLB, vecUB);
            for (jj = 0; jj < nSubSamplesG*nInputsG; jj++)
               XXG[jj] = vecOut[jj];
            delete pdfmanG;
         }
         else
         {
            if (nInputsN > 51)
                 sampler = SamplingCreateFromID(PSUADE_SAMP_LHS);
            else sampler = SamplingCreateFromID(PSUADE_SAMP_LPTAU);
            sampler->setInputBounds(nInputsG, cLower, cUpper);
            sampler->setOutputParams(1);
            sampler->setSamplingParams(nSubSamplesG, 1, 1);
            sampler->initialize(0);
            SS = new int[nSubSamplesG];
            sampler->getSamples(nSubSamplesG,nInputsG,1,XXG,YY,SS);
            delete [] SS;
            delete sampler;
         }

         for (ii2 = 0; ii2 < nSubSamplesG; ii2++)
         {
            for (jj = 0; jj < nSubSamplesN; jj++)
            {
               sCnt = 0;
               for (kk = 0; kk < nInputs; kk++)
               {
                  if (groupMembers[ii][kk] != 0)
                  {
                     mSamplePts[jj*nInputs+kk] = XXG[ii2*nInputsG+sCnt];
                     sCnt++;
                  }
               }
            }

            for (jj = 0; jj < nSubSamplesN; jj++)
            {
               sCnt = 0;
               for (kk = 0; kk < nInputs; kk++)
               {
                  if (groupMembers[ii][kk] == 0)
                  {
                     mSamplePts[jj*nInputs+kk] = XXN[jj*nInputsN+sCnt];
                     sCnt++;
                  }
               }
            }

            faPtr->evaluatePoint(nSubSamplesN, mSamplePts, YY);

            for (jj = 0; jj < nSubSamplesN; jj++)
            {
               oneSamplePt = &mSamplePts[jj*nInputs];
               ddata = constrPtr->evaluate(oneSamplePt,YY[jj],status);
               if (status == 0) YY[jj] = PSUADE_UNDEFINED;
            }

            means[ii2] = 0.0;
            sCnt = 0;
            for (jj = 0; jj < nSubSamplesN; jj++)
            {
               if (YY[jj] != PSUADE_UNDEFINED)
               {
                  means[ii2] += YY[jj];
                  sCnt++;
               }
            }
            bins[ii2] = sCnt;
            if (sCnt < 1 && printLevel >= 5)
               printf("RSMSobolG WARNING: subsample size = 0.\n"); 
            if (sCnt < 1) means[ii2] = PSUADE_UNDEFINED;
            else          means[ii2] /= (double) sCnt;

            vars[ii2] = 0.0;
            ddata = means[ii2];
            for (jj = 0; jj < nSubSamplesN; jj++)
            {
               if (YY[jj] != PSUADE_UNDEFINED)
                  vars[ii2] += (YY[jj] - ddata) * (YY[jj] - ddata);
            }
            if (sCnt < 1) vars[ii2] = PSUADE_UNDEFINED;
            else          vars[ii2] /= (double) sCnt;

            if (printLevel > 4)
            {
               printf("RSMSobolG: Group %d\n", ii+1);
               printf("  refinement = %d, size = %d (%d),", ir, sCnt,
                      nSubSamplesN);
               printf(" mean = %12.4e, var = %12.4e\n", means[ii2], vars[ii2]);
            }
         }

         totalCnt = 0;
         for (ii2 = 0; ii2 < nSubSamplesG; ii2++) totalCnt += bins[ii2];
         if (totalCnt == 0)
         {
            printf("RSMSobolG ERROR: empty constrained space.\n");
            exit(1);
         }
 
         dmean = 0.0;
         for (ii2 = 0; ii2 < nSubSamplesG; ii2++)
         {
            if (means[ii2] != PSUADE_UNDEFINED)
               dmean += means[ii2] * bins[ii2] / totalCnt;
         }

         vce = 0.0;
         for (ii2 = 0; ii2 < nSubSamplesG; ii2++)
            if (means[ii2] != PSUADE_UNDEFINED)
               vce += (means[ii2] - dmean) * (means[ii2] - dmean) *
                      bins[ii2] / totalCnt;

         ecv = 0.0;
         for (ii2 = 0; ii2 < nSubSamplesG; ii2++)
         {
            if (vars[ii2] != PSUADE_UNDEFINED)
               ecv += vars[ii2] * bins[ii2] / totalCnt;
         }

         if (printLevel > 2)
         {
            printf("Unnormalized VCE (refinement=%3d) ", ir);
            printf("for input group %3d = %12.4e\n", ii+1, vce);
         }
         if (printLevel > 3)
         {
            printf("Unnormalized ECV (refinement=%3d) ", ir); 
            printf("for input group %3d = %12.4e\n", ii+1, ecv);
         }
         if (ir == 3)
            printf("** Normalized VCE for input group %3d = %12.4e\n",
                   ii+1, vce/variance);
         currNSamples *= 2;
      }
   }
   printAsterisks(0);
    
   delete constrPtr;
   delete faPtr;
   delete [] cLower;
   delete [] cUpper;
   delete [] XXG;
   delete [] XXN;
   delete [] inputMeansG;
   delete [] inputMeansN;
   delete [] inputStdevsG;
   delete [] inputStdevsN;
   delete [] pdfFlagsG;
   delete [] pdfFlagsN;
   delete [] mSamplePts;
   delete [] YY;
   delete [] Y;
   delete [] means;
   delete [] vars;
   delete [] bins;
   delete [] iArray;
   return 0.0;
}

