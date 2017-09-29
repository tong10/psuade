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
// Definitions for the class OUUOptimizer (2 stage, 4 variable type, OUU)
// AUTHOR : CHARLES TONG
// DATE   : 2015
// ************************************************************************
#ifndef __OUUOPTIMIZER__
#define __OUUOPTIMIZER__

#include "Optimizer.h"
#include "oData.h"

// ************************************************************************
// class definition
// ************************************************************************
class OUUOptimizer : public Optimizer
{
   int optCode_;

public:

   /** constructor */
   OUUOptimizer();

   /** destructor */
   ~OUUOptimizer();

   /** run optimization 
     @param odata - an object that contains all needed data
     */
   void optimize(oData *odata);
   /** select which local optimizer to use 
     @param optcode - local optimizer code
     */
   void setLocalOptimizer(int optcode);

   /** assign operator override */ 
   OUUOptimizer& operator=(const OUUOptimizer &);
};

// ************************************************************************
// friend function
// ------------------------------------------------------------------------
extern "C"
{
   void validate(int,double *samInputs,double *samOutputs,double *errors);
}

#endif // __OUUOPTIMIZER__

