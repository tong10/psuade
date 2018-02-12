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
// Functions for Saltelli's Fourier Amplitude Sampling Test (FAST) class 
// AUTHOR : CHARLES TONG
// DATE   : 2006
// ************************************************************************
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "SFASTSampling.h"

#define PABS(x) ((x) > 0 ? x : -(x))

// ************************************************************************
// Constructor 
// ------------------------------------------------------------------------
SFASTSampling::SFASTSampling() : Sampling()
{
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
SFASTSampling::~SFASTSampling()
{
}

// ************************************************************************
// initialize the sampling data (default M = 4)
// ------------------------------------------------------------------------
int SFASTSampling::initialize(int initLevel)
{
   int    M=4, NS, *omegas, ii, jj, jj2;
   double *ranges, ds, ss, ps_pi=3.14159, ddata;

   if (nSamples_ == 0)
   {
      printf("SFASTSampling::initialize ERROR - nSamples = 0.\n");
      exit(1);
   }
   if (nInputs_ == 0 || lowerBounds_ == NULL || upperBounds_ == NULL)
   {
      printf("SFASTSampling::initialize ERROR - input not set up.\n");
      exit(1);
   }

   deleteSampleData();
   if (initLevel != 0) return 0;

   NS = nSamples_ / nInputs_;
   if (((NS-1)/(2*M)*(2*M)+1) != NS) 
   {
      NS = 2 * M * 8 + 1;
      nSamples_ = NS * nInputs_; 
      printf("SFAST::initialize - nSamples reset to %d\n", nSamples_);
   }
   allocSampleData();
   omegas = new int[nInputs_];
   NS = nSamples_ / nInputs_;
   ds = ps_pi / (double) (2 * NS);
   for (jj = 0; jj < nInputs_; jj++)
   {
      omegas[0] = (NS - 1) / (2 * M);
      calculateOmegas(nInputs_, omegas, jj);
      for (ii = 0; ii < NS; ii++)
      {
         ss = - 0.5 * ps_pi + ds * (2 * ii + 1);
         for (jj2 = 0; jj2 < nInputs_; jj2++)
         {
            ddata = ss * omegas[jj2];
            sampleMatrix_[jj*NS+ii][jj2] = 0.5 + asin(sin(ddata))/ps_pi;
         }
      }
   }

   ranges = new double[nInputs_];
   for (jj = 0; jj < nInputs_; jj++) 
      ranges[jj] = upperBounds_[jj] - lowerBounds_[jj];

   for (ii = 0; ii < nSamples_; ii++) 
   {
      for (jj = 0; jj < nInputs_; jj++)
      {
         ddata = sampleMatrix_[ii][jj];
         sampleMatrix_[ii][jj] = ddata * ranges[jj] + lowerBounds_[jj];
      }
   } 

   if (printLevel_ > 4)
   {
      printf("SFASTSampling::initialize: nSamples = %d\n", nSamples_);
      printf("SFASTSampling::initialize: nInputs  = %d\n", nInputs_);
      printf("SFASTSampling::initialize: nOutputs = %d\n", nOutputs_);
      for (jj = 0; jj < nInputs_; jj++)
         printf("    SFASTSampling input %3d = [%e %e]\n", jj+1,
                lowerBounds_[jj], upperBounds_[jj]);
   }

   delete [] ranges;
   delete [] omegas;
   return 0;
}

// ************************************************************************
// refine the sample space
// ------------------------------------------------------------------------
int SFASTSampling::refine(int refineRatio, int randomize, double thresh,
                          int nSamples, double *sampleErrors)
{
   int    *omegas, NS, omi, M=4, nLevels, *oldStates, oldNS;
   int    ii, jj, jj2, fail;
   double *ranges, **oldSamples, *oldOutputs, ds, ss, ps_pi=3.14159, ddata;

   (void) randomize;
   (void) thresh;
   (void) nSamples;
   (void) sampleErrors;

   omegas = new int[nInputs_];
   NS  = nSamples_ / nInputs_;
   omi = (NS - 1) / (2 * M);
   if (((omi*2*M+1)*NS) != nSamples_)
   {
      printf("SFASTSampling refine ERROR : invalid sample size.\n");
      delete [] omegas;
      exit(1);
   } 
   
   nLevels     = refineRatio;
   oldSamples  = sampleMatrix_;
   oldOutputs  = sampleOutput_;
   oldStates   = sampleStates_;
   oldNS       = NS;
   sampleMatrix_ = NULL;
   sampleOutput_ = NULL;
   sampleStates_ = NULL;
   deleteSampleData();
   NS = (NS - 1) * 2 + 1;
   nSamples_ = NS * nInputs_;
   allocSampleData();

   ds = ps_pi / (double) (2 * NS);
   for (jj = 0; jj < nInputs_; jj++)
   {
      omegas[0] = (NS - 1) / (2 * M);
      calculateOmegas(nInputs_, omegas, jj);
      for (ii = 0; ii < NS; ii++)
      {
         ss = - 0.5 * ps_pi + ds * (2 * ii + 1);
         for (jj2 = 0; jj2 < nInputs_; jj2++)
         {
            ddata = ss * omegas[jj2];
            sampleMatrix_[jj*NS+ii][jj2] = 0.5 + asin(sin(ddata))/ps_pi;
         }
      }
   }

   ranges = new double[nInputs_];
   for (jj = 0; jj < nInputs_; jj++) 
      ranges[jj] = upperBounds_[jj] - lowerBounds_[jj];

   for (ii = 0; ii < nSamples_; ii++) 
   {
      for (jj = 0; jj < nInputs_; jj++)
      {
         ddata = sampleMatrix_[ii][jj];
         sampleMatrix_[ii][jj] = ddata * ranges[jj] + lowerBounds_[jj];
      }
   } 

   fail = 0;
   for (jj = 0; jj < nInputs_; jj++)
   {
      for (ii = 0; ii < oldNS; ii++)
      {
         for (jj2 = 0; jj2 < nInputs_; jj2++)
         {
            ddata = oldSamples[jj*oldNS+ii][jj2] - 
                    sampleMatrix_[jj*NS+ii*2][jj2];
            if (PABS(ddata) > 1.0e-13) fail = 1;
         }
         for (jj2 = 0; jj2 < nOutputs_; jj2++)
            sampleOutput_[(jj*NS+ii*2)*nOutputs_+jj2] = 
               oldOutputs[(jj*oldNS+ii)*nOutputs_+jj2];
         sampleStates_[jj*NS+ii*2] = oldStates[jj*oldNS+ii];
      }
   }
   if (fail == 1)
   {
      printf("SFASTSampling fails checking.\n");
      exit(1);
   }

   if (printLevel_ > 4)
   {
      printf("SFASTSampling::refine: nSamples = %d\n", nSamples_);
      printf("SFASTSampling::refine: nInputs  = %d\n", nInputs_);
      printf("SFASTSampling::refine: nOutputs = %d\n", nOutputs_);
      for (jj = 0; jj < nInputs_; jj++)
         printf("    SFASTSampling input %3d = [%e %e]\n", jj+1,
                lowerBounds_[jj], upperBounds_[jj]);
   }

   delete [] omegas;
   delete [] ranges;
   return 0;
}

// ************************************************************************
// calculate frequencies
// ------------------------------------------------------------------------
int SFASTSampling::calculateOmegas(int nInputs, int *omegas, int index)
{
   int om1, step, offset, ii;

   om1    = omegas[0];
   step   = om1 / 16;
   offset = 1;
   for (ii = 0; ii < nInputs; ii++) 
   {
      if (ii != index)
      {
         omegas[ii] = offset;
         offset += step;
      }
   }
   omegas[index] = om1;
   return 0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
SFASTSampling& SFASTSampling::operator=(const SFASTSampling &)
{
   printf("SFASTSampling operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

