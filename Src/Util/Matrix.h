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
// Matrix functions
// AUTHOR : CHARLES TONG
// DATE   : 2008
// ************************************************************************

#ifndef __MATRIXH__
#define __MATRIXH__

/**
 * @name Matrix class
 *
 **/
/*@{*/

#include "Vector.h"

// ************************************************************************
// class definition
// ************************************************************************

class Matrix
{
   int    nRows_, nCols_;
   double **Mat_;
   int    status_;
   double determinant_;

public:

   Matrix();
   // Copy constructor by Bill Oliver
   Matrix(const Matrix & ma);
   // overload = operator
   Matrix & operator=(const Matrix & ma);
   ~Matrix();
   int    nrows();
   int    ncols();
   int    load(Matrix &);
   int    setDim(int, int);
   void   setEntry(const int, const int, const double);
   double getEntry(const int, const int);
   double getDeterminant();
   int    submatrix(Matrix &, const int, const int *);
   int    CholDecompose();
   void   CholMatvec(Vector &);
   void   CholSolve(Vector &);
   void   CholTSolve(Vector &);
   void   print();

private:
   double computeDeterminant(int, double **);
};

/*@}*/

#endif /* __MATRIXH__ */

