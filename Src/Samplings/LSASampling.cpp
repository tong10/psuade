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
// Functions for the Local Sensitivity Analysis sampling class
// AUTHOR : CHARLES TONG
// DATE   : 2012
// ************************************************************************
using namespace std;
#include "Util/sysdef.h"
#include "Util/PsuadeUtil.h"
#include "LSASampling.h"

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
LSASampling::LSASampling() : Sampling()
{
   samplingID_ = PSUADE_SAMP_LSA;
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
LSASampling::~LSASampling()
{
}

// ************************************************************************
// initialization 
// ------------------------------------------------------------------------
int LSASampling::initialize(int initLevel)
{
   int    ii, ss;
   double ddata;

   if (nInputs_ == 0 || lowerBounds_ == NULL || upperBounds_ == NULL)
   {
      printf("ss::initialize ERROR - input not set up.\n");
      exit(1);
   }

   if (printLevel_ > 4)
   {
      printf("LSASampling::initialize: nSamples = %d\n", nSamples_);
      printf("LSASampling::initialize: nInputs  = %d\n", nInputs_);
      printf("LSASampling::initialize: nOutputs = %d\n", nOutputs_);
      for (ii = 0; ii < nInputs_; ii++)
         printf("    LSASampling input %3d = [%e %e]\n", ii+1,
                lowerBounds_[ii], upperBounds_[ii]);
   }

   deleteSampleData();
   if (initLevel != 0) return 0;

   nSamples_ = nInputs_ + 1;
   allocSampleData();
   for (ss = 0; ss < nSamples_; ss++)
   { 
      for (ii = 0; ii < nInputs_; ii++) 
         sampleMatrix_[ss][ii] = 0.5*(upperBounds_[ii]+lowerBounds_[ii]);
      if (ss > 0)
      {
         ddata = PSUADE_drand();
         if (ddata >= 0.5) ddata = 1.0;
         else              ddata = 0.0;
         sampleMatrix_[ss][ss-1] = lowerBounds_[ss-1] + ddata *
                (upperBounds_[ss-1] - lowerBounds_[ss-1]);
      }
   }
   return 0;
}

// ************************************************************************
// refine 
// ------------------------------------------------------------------------
int LSASampling::refine(int, int, double, int, double *)
{
   printf("LSASampling ERROR - refine not available.\n");
   return -1;
}

