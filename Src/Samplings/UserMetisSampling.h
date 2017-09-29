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
// Definition for the sampling class based on the Metis partitioning
// AUTHOR : CHARLES TONG
// DATE   : 2006
// ************************************************************************
#ifndef __USERMETISSAMPLINGH__
#define __USERMETISSAMPLINGH__

#include <iostream>
#include <string.h>
using namespace std;
#include "Sampling.h"

/**
 * @name Metis (domain decomposition) samping method
 *
 **/

// ************************************************************************
// class definition
// ************************************************************************
                                                                                
class UserMetisSampling : public Sampling 
{
   int n1d_;
   int nAggrs_;
   int *aggrCnts_;
   int **aggrLabels_;
   int *cellsOccupied_;

public:

   UserMetisSampling();
   ~UserMetisSampling();

   int initialize(int);
   int refine(int, int, double, int, double *);
   int setParam(string);
};

#endif // __USERMETISSAMPLINGH__

