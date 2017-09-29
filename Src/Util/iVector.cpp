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
// iVector functions
// AUTHOR : CHARLES TONG
// DATE   : 2010
// ************************************************************************

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "iVector.h"

//#define PS_DEBUG

// ************************************************************************
// Constructor
// ------------------------------------------------------------------------
iVector::iVector()
{
#ifdef PS_DEBUG
   printf("iVector constructor\n");
#endif
   length_ = 0;
   Vec_ = NULL;
#ifdef PS_DEBUG
   printf("iVector constructor ends\n");
#endif
}

// ************************************************************************
// Copy Constructor by Bill Oliver
// ------------------------------------------------------------------------
iVector::iVector(const iVector & v){

  length_ = v.length_;
  Vec_ = new int[length_];
  for(int i = 0; i < length_; i++){

    Vec_[i] = v.Vec_[i];
  }

}

// ************************************************************************
// overload operator= by Bill Oliver
// ------------------------------------------------------------------------
iVector & iVector::operator=(const iVector & v){

  if(this == &v)  return *this;
  length_ = v.length_;
  delete [] Vec_;
  Vec_ = new int[length_];
  for(int i = 0; i < length_; i++)
    Vec_[i] = v.Vec_[i];

  return *this;

}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
iVector::~iVector()
{
#ifdef PS_DEBUG
   printf("iVector destructor\n");
#endif
   if (Vec_ != NULL) delete [] Vec_;
   Vec_ = NULL;
   length_ = 0;
#ifdef PS_DEBUG
   printf("iVector destructor ends\n");
#endif
}

// ************************************************************************
// get length 
// ------------------------------------------------------------------------
int iVector::length() 
{
   return length_;
}

// ************************************************************************
// load vector
// ------------------------------------------------------------------------
int iVector::load(iVector &inVec)
{
#ifdef PS_DEBUG
   printf("iVector load\n");
#endif
   assert(this != &inVec);
   if (Vec_ != NULL) delete [] Vec_;
   Vec_ = NULL;
   length_ = inVec.length();
   if (length_ <= 0) return -1;
   Vec_ = new int[length_];
   assert(Vec_ != NULL);
   for (int ii = 0; ii < length_; ii++) Vec_[ii] = inVec[ii];
#ifdef PS_DEBUG
   printf("iVector load ends\n");
#endif
   return 0;
}

// ************************************************************************
// set dimension
// ------------------------------------------------------------------------
int iVector::setLength(int leng)
{
#ifdef PS_DEBUG
   printf("iVector setLength\n");
#endif
   if (Vec_ != NULL) delete [] Vec_;
   Vec_ = NULL;
   assert(leng > 0);
   length_ = leng;
   Vec_ = new int[leng];
   assert(Vec_ != NULL);
   for (int ii = 0; ii < leng; ii++) Vec_[ii] = 0;
#ifdef PS_DEBUG
   printf("iVector setLength ends\n");
#endif
   return 0;
}

// ************************************************************************
// load vector
// ------------------------------------------------------------------------
int iVector::load(int leng, int *data)
{
#ifdef PS_DEBUG
   printf("iVector load, length = %d\n", leng);
#endif
   if (Vec_ != NULL) delete [] Vec_;
   Vec_ = NULL;
   assert(leng > 0);
   assert(data != NULL);
   length_ = leng;
   Vec_ = new int[leng];
   assert(Vec_ != NULL);
   for (int ii = 0; ii < leng; ii++) Vec_[ii] = data[ii];
#ifdef PS_DEBUG
   printf("iVector load ends\n");
#endif
   return 0;
}

// ************************************************************************
// get entry
// ------------------------------------------------------------------------
int& iVector::operator[](int ind) 
{
   if (ind < 0 || ind >= length_)
   {
      printf("iVector operator[] ERROR: index = %d (0, %d)\n",ind,length_=1);
      exit(1);
   }
   return Vec_[ind];
}

// ************************************************************************
// get vector
// ------------------------------------------------------------------------
int *iVector::getIVector() 
{
   return Vec_;
}


