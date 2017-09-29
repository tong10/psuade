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
// Vector functions
// AUTHOR : CHARLES TONG
// DATE   : 2008
// ************************************************************************

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "Util/Vector.h"

//#define PS_DEBUG

// ************************************************************************
// Constructor
// ------------------------------------------------------------------------
Vector::Vector()
{
#ifdef PS_DEBUG
   printf("Vector constructor\n");
#endif
   length_ = 0;
   Vec_ = NULL;
#ifdef PS_DEBUG
   printf("Vector constructor ends\n");
#endif
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
Vector::~Vector()
{
#ifdef PS_DEBUG
   printf("Vector destructor\n");
#endif
   if (Vec_ != NULL) delete [] Vec_;
   Vec_ = NULL;
   length_ = 0;
#ifdef PS_DEBUG
   printf("Vector destructor ends\n");
#endif
}

// ************************************************************************
// get length 
// ------------------------------------------------------------------------
int Vector::length() 
{
   return length_;
}

// ************************************************************************
// load vector
// ------------------------------------------------------------------------
int Vector::load(Vector &inVec)
{
#ifdef PS_DEBUG
   printf("Vector load\n");
#endif
   assert(this != &inVec);
   if (Vec_ != NULL) delete [] Vec_;
   Vec_ = NULL;
   length_ = inVec.length();
   if (length_ <= 0) return -1;
   Vec_ = new double[length_];
   assert(Vec_ != NULL);
   for (int ii = 0; ii < length_; ii++) Vec_[ii] = inVec[ii];
#ifdef PS_DEBUG
   printf("Vector load ends\n");
#endif
   return 0;
}

// ************************************************************************
// set dimension
// ------------------------------------------------------------------------
int Vector::setLength(int leng)
{
#ifdef PS_DEBUG
   printf("Vector setLength\n");
#endif
   if (Vec_ != NULL) delete [] Vec_;
   Vec_ = NULL;
   assert(leng > 0);
   length_ = leng;
   Vec_ = new double[leng];
   assert(Vec_ != NULL);
   for (int ii = 0; ii < leng; ii++) Vec_[ii] = 0.0;
#ifdef PS_DEBUG
   printf("Vector setLength ends\n");
#endif
   return 0;
}

// ************************************************************************
// load vector
// ------------------------------------------------------------------------
int Vector::load(int leng, double *data)
{
#ifdef PS_DEBUG
   printf("Vector load, length = %d\n", leng);
#endif
   if (Vec_ != NULL) delete [] Vec_;
   Vec_ = NULL;
   assert(leng > 0);
   assert(data != NULL);
   length_ = leng;
   Vec_ = new double[leng];
   assert(Vec_ != NULL);
   for (int ii = 0; ii < leng; ii++) Vec_[ii] = data[ii];
#ifdef PS_DEBUG
   printf("Vector load ends\n");
#endif
   return 0;
}

// ************************************************************************
// get entry
// ------------------------------------------------------------------------
double& Vector::operator[](int ind) 
{
   if (ind < 0 || ind >= length_)
   {
      printf("Vector operator[] ERROR: index = %d (0, %d)\n",ind,length_=1);
      exit(1);
   }
   return Vec_[ind];
}

// ************************************************************************
// get vector
// ------------------------------------------------------------------------
double *Vector::getDVector() 
{
   return Vec_;
}

// ************************************************************************
// add to vector
// ------------------------------------------------------------------------
int Vector::add(int leng, double *data)
{
#ifdef PS_DEBUG
   printf("Vector add, length = %d\n", leng);
#endif
   int    ii;
   double *tmpVec = Vec_;
   Vec_ = new double[leng+length_];
   for (ii = 0; ii < length_; ii++) Vec_[ii] = tmpVec[ii];
   if (data == NULL)
        for (ii = 0; ii < leng; ii++) Vec_[length_+ii] = 0.0;
   else for (ii = 0; ii < leng; ii++) Vec_[length_+ii] = data[ii];
   delete [] tmpVec;
   length_ += leng;
#ifdef PS_DEBUG
   printf("Vector add ends\n");
#endif
   return 0;
}

