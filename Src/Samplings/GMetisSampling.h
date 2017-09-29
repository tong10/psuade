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
// Definition for the sampling class based on the GMetis partitioning
// AUTHOR : CHARLES TONG
// DATE   : 2010
// ************************************************************************
#ifndef __GMETISSAMPLINGH__
#define __GMETISSAMPLINGH__

#include <iostream>
#include <string.h>
using    namespace std;
#include "Sampling.h"

/**
 * @name GMetis (domain decomposition) samping method
 *
 **/

// ************************************************************************
// class definition
// ************************************************************************
class GMetisSampling : public Sampling 
{
   int refineType_;
   int refineSize_;
   int n1d_;
   int nAggrs_;
   int *aggrCnts_;
   int **aggrLabels_;
   int graphN_;
   int *graphI_;
   int *graphJ_;
   int *cellsOccupied_;

public:

   GMetisSampling();
   ~GMetisSampling();

   int initialize(int);
   int refine(int, int, double, int, double *);
   int setParam(string sparam);
};

#endif // __GMETISSAMPLINGH__

