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
// Functions for the class FASTAnalyzer  
// AUTHOR : CHARLES TONG
// DATE   : 2005
// ************************************************************************
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "FASTAnalyzer.h"
#include "sysdef.h"
#include "PsuadeUtil.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

#define  PSUADE_FAST_MaxDimension  50
                                                                                
static unsigned long
PSUADE_FAST_OMEGA[PSUADE_FAST_MaxDimension] =
{
      1,    3,    1,    5,   11,    1,   17,   23,   19,   25,
     41,   31,   23,   87,   67,   73,   85,  143,  149,   99,
    119,  237,  267,  283,  151,  385,  157,  215,  449,  163,
    337,  253,  375,  441,  673,  773,  875,  873,  587,  849,
    623,  637,  891,  943, 1171, 1225, 1335, 1725, 1663, 2019
};

static unsigned long
PSUADE_FAST_DELTA[PSUADE_FAST_MaxDimension] =
{
      4,    8,    6,   10,   20,   22,   32,   40,   38,   26,
     56,   62,   46,   76,   96,   60,   86,  126,  134,  112,
     92,  128,  154,  196,   34,  416,  106,  208,  328,  198,
    382,   88,  348,  186,  140,  170,  284,  568,  302,  438,
    410,  248,  448,  388,  596,  216,  100,  488,  166,    0
};
                                                                                
// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
FASTAnalyzer::FASTAnalyzer() : Analyzer()
{
   setName("FAST");
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
FASTAnalyzer::~FASTAnalyzer() 
{ 
} 

// ************************************************************************ 
// perform analysis 
// ------------------------------------------------------------------------
double FASTAnalyzer::analyze(aData &adata)
{
   int    nInputs, nOutputs, nSamples, outputID, M, ii, N2, maxInd;
   int    ss, printLevel, count;
   double *fourierCoefs, *Y, *YY, *fourierCoefs2, retdata;
   double maxData, fsum;

   printLevel = adata.printLevel_;
   nInputs    = adata.nInputs_;
   nOutputs   = adata.nOutputs_;
   nSamples   = adata.nSamples_;
   Y          = adata.sampleOutputs_;
   outputID   = adata.outputID_;
   if (adata.inputPDFs_ != NULL)
   {
      count = 0;
      for (ii = 0; ii < nInputs; ii++) count += adata.inputPDFs_[ii];
      if (count > 0)
      {
         printf("FAST INFO: some inputs have non-uniform PDFs, but\n");
         printf("     they are not relevant in this analysis.\n");
         printf("     (To perform this analysis with desired distributions,\n");
         printf("     first create a FAST sample, then prescribe PDFs and\n");
         printf("     run the sample through 'pdfconvert' before running\n");
         printf("     the simulations.\n");
      }
   }

   if (nInputs <= 0 || nOutputs <= 0 || nSamples <= 0)
   {
      printf("FAST ERROR: invalid arguments.\n");
      printf("    nInputs  = %d\n", nInputs);
      printf("    nOutputs = %d\n", nOutputs);
      printf("    nSamples = %d\n", nSamples);
      return PSUADE_UNDEFINED;
   } 
   if (outputID < 0 || outputID >= nOutputs)
   {
      printf("FAST ERROR: invalid outputID.\n");
      printf("    outputID = %d\n", outputID+1);
      return PSUADE_UNDEFINED;
   } 
   
   YY = new double[nSamples];
   for (ss = 0; ss < nSamples; ss++) YY[ss] = Y[ss*nOutputs+outputID];
   M = computeCoefficents(nSamples, nInputs, YY, &fourierCoefs,
                          printLevel);
   
   printEquals(0);
   printf("* Fourier Amplitude Sampling Test (FAST) coefficients\n");
   printDashes(0);
   printf("* M = %d\n", M);
   fsum = 0.0;
   for (ii = 0; ii < nInputs; ii++)
   {
      printf("* Input %4d = %14.6e\n", ii+1, fourierCoefs[ii]);
      fsum += fourierCoefs[ii];
   }
   printf("* Sum of FAST coefficients = %14.6e\n", fsum);
   printf("* FAST variance            = %14.6e\n", fourierCoefs[nInputs]);

   if (M % 2 == 0)
   {
      for (ss = 0; ss < nSamples; ss+=2) YY[ss/2] = Y[ss*nOutputs+outputID];
      N2 = (nSamples + 1) / 2;
      M  = computeCoefficents(N2, nInputs, YY, &fourierCoefs2, 0);
      if (printLevel >= 2)
      {
         printDashes(0);
         printf("* Fourier Amplitude Sampling Test (FAST) coarse coefficients\n");
         printDashes(0);
         for (ii = 0; ii < nInputs; ii++)
            printf("* Input %4d = %14.6e\n", ii+1, fourierCoefs2[ii]);
         printf("* FAST variance            = %14.6e\n", fourierCoefs2[nInputs]);
      }
   }
   printEquals(0);

   retdata = 1.0;
   if (M % 2 == 0)
   {
      maxInd = 0;
      maxData = fourierCoefs[0];
      for (ii = 1; ii < nInputs; ii++)
      {
         if (fourierCoefs[0] > maxData)
         {
            maxData = fourierCoefs[0];
            maxInd = ii;
         }
      }
      retdata = PABS(1.0 - fourierCoefs2[maxInd]/fourierCoefs[maxInd]);
      delete [] fourierCoefs2;
   }
   delete [] YY;
   delete [] fourierCoefs;
   if (printLevel >= 2)
      printf("* FAST convergence rate = %12.4e\n",retdata);
   return retdata;
}

// ************************************************************************
// calculate frequencies
// ------------------------------------------------------------------------
int FASTAnalyzer::calculateOmegas(int nInputs, int nSamples, int *omegas)
{
   int ii;

   omegas[0] = PSUADE_FAST_OMEGA[nInputs-1];
   for (ii = 1; ii < nInputs; ii++)
      omegas[ii] = omegas[ii-1] + PSUADE_FAST_DELTA[nInputs-1-ii];
   return 0;
}
                                                                                
// ************************************************************************ 
// compute Fourier coefficients
// ------------------------------------------------------------------------
int FASTAnalyzer::computeCoefficents(int nSamples, int nInputs, double *Y,
                                     double **fourierCoefs, int flag)
{
   int    *omegas, M, ii, jj, ss, N;
   double *fourierReal, *fourierImag, ps_pi=3.14159, ds, freq;
   double fastReal, fastImag, fastCoef, dataReal, dataImag;

   omegas = new int[nInputs];
   calculateOmegas(nInputs, nSamples, omegas);
   M = (nSamples - 1) / (2 * omegas[nInputs-1]);
   if ((2 * M * omegas[nInputs-1] + 1) != nSamples)
   {
      printf("FAST ERROR: not FAST samples ?\n");
      delete [] omegas;
      exit(1);
   } 
   
   fourierReal = new double[M*nInputs];
   fourierImag = new double[M*nInputs];
   for (ii = 0; ii < M*nInputs; ii++)
      fourierReal[ii] = fourierImag[ii] = 0.0;
   ds = ps_pi / (double) nSamples;
  
   for (ii = 0; ii < M; ii++)
   {
      for (jj = 0; jj < nInputs; jj++)
      {
         for (ss = 0; ss < nSamples; ss++)
         {
            freq = - ps_pi / 2.0 + ds * 0.5 * (2 * ss + 1);
            fourierReal[ii*nInputs+jj] += 
               Y[ss]*cos((ii+1)*omegas[jj]*freq)*ds;
            fourierImag[ii*nInputs+jj] += 
               Y[ss]*sin((ii+1)*omegas[jj]*freq)*ds;
         }
         for (ss = 0; ss < (nSamples-1)/2; ss++)
         {
            freq = - ps_pi + ds * (ss + 1);
            fourierReal[ii*nInputs+jj] += 
               Y[(nSamples+1)/2-ss-2]*cos((ii+1)*omegas[jj]*freq)*ds;
            fourierImag[ii*nInputs+jj] += 
               Y[(nSamples+1)/2-ss-2]*sin((ii+1)*omegas[jj]*freq)*ds;
         }
         for (ss = 0; ss < (nSamples-1)/2; ss++)
         {
            freq = ps_pi - ds * (ss + 1);
            fourierReal[ii*nInputs+jj] += 
               Y[(nSamples+1)/2+ss]*cos((ii+1)*omegas[jj]*freq)*ds;
            fourierImag[ii*nInputs+jj] += 
               Y[(nSamples+1)/2+ss]*sin((ii+1)*omegas[jj]*freq)*ds;
         }
      }
   }
   for (ii = 0; ii < M*nInputs; ii++)
   {
      fourierReal[ii] /= (2.0 * ps_pi);
      fourierImag[ii] /= (2.0 * ps_pi);
   }
   (*fourierCoefs) = new double[nInputs+1];
   for (jj = 0; jj < nInputs; jj++)
   {
      (*fourierCoefs)[jj] = 0.0;
      for (ii = 0; ii < M; ii++)
         (*fourierCoefs)[jj] += 
             (fourierReal[ii*nInputs+jj]*fourierReal[ii*nInputs+jj]);
      for (ii = 0; ii < M; ii++)
         (*fourierCoefs)[jj] += 
             (fourierImag[ii*nInputs+jj]*fourierImag[ii*nInputs+jj]);
      (*fourierCoefs)[jj] *= 2.0;
   }

   N = M * omegas[nInputs-1];
   fastReal = fastImag = 0.0;
   if (flag >= 3)
   {
      for (jj = 0; jj < nInputs; jj++)
         printf("FAST: input %4d fundamental frequency = %d\n",jj+1,
                omegas[jj]);
   }
   for (ii = 0; ii < N; ii++)
   {
      dataReal = dataImag = 0.0;
      for (ss = 0; ss < nSamples; ss++)
      {
         freq = - ps_pi / 2.0 + ds * 0.5 * (2 * ss + 1);
         dataReal += Y[ss]*cos((ii+1)*freq)*ds;
         dataImag += Y[ss]*sin((ii+1)*freq)*ds;
      }
      for (ss = 0; ss < (nSamples-1)/2; ss++)
      {
         freq = - ps_pi + ds * (ss + 1);
         dataReal += Y[(nSamples+1)/2-ss-2]*cos((ii+1)*freq)*ds;
         dataImag += Y[(nSamples+1)/2-ss-2]*sin((ii+1)*freq)*ds;
      }
      for (ss = 0; ss < (nSamples-1)/2; ss++)
      {
         freq = ps_pi - ds * (ss + 1);
         dataReal += Y[(nSamples+1)/2+ss]*cos((ii+1)*freq)*ds;
         dataImag += Y[(nSamples+1)/2+ss]*sin((ii+1)*freq)*ds;
      }
      dataReal /= (2.0 * ps_pi);
      dataImag /= (2.0 * ps_pi);
      fastReal += dataReal * dataReal;
      fastImag += dataImag * dataImag;
      if (flag >= 3)
      {
         printf("FAST: frequency %5d : data = %9.1e (%9.1e %9.1e) ",
                ii+1,dataReal*dataReal+dataImag*dataImag,dataReal,dataImag);
         for (jj = 0; jj < nInputs; jj++)
            if ((ii + 1) / omegas[jj] * omegas[jj] == (ii + 1))
               printf("(%4d) ", jj+1);
         printf("\n");
      }
   }
   fastCoef = 2.0 * (fastReal + fastImag);

   for (jj = 0; jj < nInputs; jj++) (*fourierCoefs)[jj] /= fastCoef;
   (*fourierCoefs)[nInputs] = fastCoef;
   delete [] fourierReal;
   delete [] fourierImag;
   delete [] omegas;
   return M;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
FASTAnalyzer& FASTAnalyzer::operator=(const FASTAnalyzer &)
{
   printf("FAST operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

