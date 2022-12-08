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
// Definition for the GPBasedSampling class
// AUTHOR : CHARLES TONG
// DATE   : 2017
// ************************************************************************
#ifndef __GPBASEDSAMPLINGH__
#define __GPBASEDSAMPLINGH__

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Sampling.h"
#include "psVector.h"
#include "psMatrix.h"

/**
 * @name sequential samping method
 *
 **/

// ************************************************************************
// class definition
// ************************************************************************
class GPBasedSampling: public Sampling 
{
   psVector VecGPParams_;
   int MMVorMAV_;

public:

   /** constructor */
   GPBasedSampling();

   /** destructor */
   ~GPBasedSampling();

   /** initialization 
       @param flag: flag to signal how far to initialize
    */
   int initialize(int flag);

   /** set MMV  (minimize maximum variances)
    */
   void setMMV();

   /** set MAV (minimize average variances) 
    */
   void setMAV();

   /** set MAV (minimize average variances) 
    */
   double evaluateCandidateSet(psMatrix &, psMatrix &);

   /** This function overloads the assignment operator
       @param obj : Sampling object
    */
   GPBasedSampling& operator=(const GPBasedSampling &);

private:
   void constructCMatrix(psMatrix &, psMatrix &, psVector);
   int  genDesigns(psVector&,psMatrix&,int,psMatrix &,psVector &);
   int  genDesigns2(psVector&,psMatrix&,int,psMatrix &,psVector &);
   int  genDesignsUltimate(psVector&,psMatrix&,int,psMatrix &);
   int  genInitialDesigns(psVector,psMatrix&,int,psMatrix &,psVector &);
   int  genInitialDesigns2(psVector,psMatrix&,int,psMatrix &,psVector &);
   int  evaluateUltimate(int, psIVector &, int, int, int, psVector &,
                         psMatrix &, psIVector &, double &, long &);
   int  evaluateUltimateOpt(int, int, psVector &, psMatrix &, 
                            psIVector &, double &);
   int  printDesigns(psMatrix &, FILE *);
};

#endif // __GPBASEDSAMPLINGH__

