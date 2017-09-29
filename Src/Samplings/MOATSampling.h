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
// Definition for the Morris one-at-a-time sampling class
// AUTHOR : CHARLES TONG
// DATE   : 2004
// ************************************************************************
#ifndef __MOATSAMPLINGH_
#define __MOATSAMPLINGH_

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Sampling.h"
#include "../DataIO/PsuadeData.h"

/**
 * @name Morris one-at-a-time samping method
 *
 **/

// ************************************************************************
// class definition
// ************************************************************************
class MOATSampling: public Sampling 
{
   int P_;
   int *inputSubset_;

public:

   MOATSampling();
   ~MOATSampling();

   int initialize(int);
   int refine(int, int, double, int, double *);
   int repair(char *, int);
   int genRepair(int, double *, double *);
   int genRepair(PsuadeData *);

private:

   int initializeHighDimension();
   int generate(double **);
   int merge();
   int checkSample(int, int, double **);
   int checkSample2(int, int, double *);
};

#endif // __MOATSAMPLINGH_

