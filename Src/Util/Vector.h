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
// psVector functions
// AUTHOR : CHARLES TONG
// DATE   : 2008
// ************************************************************************

#ifndef __VECTORH__
#define __VECTORH__

/**
 * @name psVector class
 *
 **/
/*@{*/

// ************************************************************************
// class definition
// ************************************************************************

class psVector
{

   int    length_;
   double *Vec_;

public:

   psVector();
   // Copy Constructor by Bill Oliver
   psVector(const psVector & v);

   // Overload operator= by Bill Oliver
   psVector & operator=(const psVector & v);

   ~psVector();
   int     length();
   int     load(psVector &);
   int     setLength(int);
   int     load(int, double *);
   int     add(int, double *);
   double& operator[](int ind);
   double  *getDVector();
};

/*@}*/

#endif /* __VECTORH__ */

