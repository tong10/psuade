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
// Definition for the abstract class MOATConstraints
// AUTHOR : CHARLES TONG
// DATE   : 2009
// ************************************************************************
                                                                                             
#ifndef __MOATCONSTRAINTS__
#define __MOATCONSTRAINTS__

#include "FuncApprox/FuncApprox.h"
#include "DataIO/PsuadeData.h"

/**
 * @name Defining output constraints
 *
 **/

// ************************************************************************
// class definition
// ************************************************************************
                                                                                
class MOATConstraints
{

   int        nConstraints_;
   int        *constraintNInputs_;
   int        **constraintInputIndices_;
   double     **constraintInputValues_;
   int        nInputs_;
   double     *XLBounds_;
   double     *XUBounds_;
   FuncApprox **constraintFAs_;
   double     *YLBounds_;
   double           *YUBounds_;

public:

   MOATConstraints();
   ~MOATConstraints();

   int    initialize(PsuadeData *psuadeIO);
   double getScale(double *, int, int &);
};

#endif // __MOATCONSTRAINTS__

