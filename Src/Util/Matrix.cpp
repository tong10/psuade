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
// psMatrix functions (modified from Pelikan's (University of Cincinati) code
// AUTHOR : CHARLES TONG
// DATE   : 2008
// ************************************************************************
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "Matrix.h"

//#define PS_DEBUG 1
// ************************************************************************
// ************************************************************************
// external functions
// ------------------------------------------------------------------------
extern "C" {
  void dgetri_(int *, double *, int *, int *, double *, int *, int *);
  void dgels_(char *, int *, int *, int *, double *, int *,
              double *, int *, double *, int *, int *);
  void dsyev_(char *, char *, int *, double *, int *, double *,
              double *, int *, int *);
  void dpotrf_(char *, int *, double *, int *, int *);
  void dpotrs_(char *, int *, int *, double *, int *, double *,int *,int *);
  void dgetrf_(int *, int *, double *, int *, int *, int *);
  void dgetrs_(char *,int *,int *,double *,int *,int *,double *,int *,int *);
}

// ************************************************************************
// Constructor
// ------------------------------------------------------------------------
psMatrix::psMatrix()
{
#ifdef PS_DEBUG
   printf("psMatrix constructor\n");
#endif
   nRows_ = 0;
   nCols_ = 0;
   Mat_ = NULL;
   status_ = 0;
   determinant_ = 0.0;
   pivots_ = 0;
#ifdef PS_DEBUG
   printf("psMatrix constructor ends\n");
#endif
}

// ************************************************************************
// Copy Constructor by Bill Oliver
// ------------------------------------------------------------------------
psMatrix::psMatrix(const psMatrix & ma)
{
   int ii, jj;
   nRows_ = ma.nRows_;
   nCols_ = ma.nCols_;
   status_ = ma.status_;
   if (nRows_ > 0 && nCols_ > 0)
   {
      Mat_ = new double*[nRows_];
      for (ii = 0; ii < nRows_; ii++)
      {
         Mat_[ii] = new double[nCols_];
         for(jj = 0; jj < nCols_; jj++)
            Mat_[ii][jj] = ma.Mat_[ii][jj];
      }
   }
   determinant_ = ma.determinant_;
}

// ************************************************************************
// operator=  by Bill Oliver
// ------------------------------------------------------------------------
psMatrix & psMatrix::operator=(const psMatrix & ma)
{
   int ii, jj;
   if (this == &ma) return *this;

   if (Mat_ != NULL)
   {
      for(ii = 0; ii < nRows_; ii++) delete [] Mat_[ii];
      delete [] Mat_;
   }
   Mat_ = NULL;

   nRows_ = ma.nRows_;
   nCols_ = ma.nCols_;
   status_ = ma.status_;
   determinant_ = ma.determinant_;
   if (nRows_ > 0 && nCols_ > 0)
   {
      Mat_ = new double*[nRows_];
      for(ii = 0; ii < nRows_; ii++)
      {
         Mat_[ii] = new double[nCols_];
         for(jj = 0; jj < nCols_; jj++) Mat_[ii][jj] = ma.Mat_[ii][jj];
      }
   }
   return *this;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
psMatrix::~psMatrix()
{
#ifdef PS_DEBUG
   printf("psMatrix destructor\n");
#endif
   if (Mat_ != NULL)
   {
      for (int ii = 0; ii < nRows_; ii++)
         if (Mat_[ii] != NULL) delete [] Mat_[ii];
      delete [] Mat_;
   }
   nRows_ = nCols_ = 0;
   Mat_ = NULL;
   if (pivots_ != NULL) delete [] pivots_;
#ifdef PS_DEBUG
   printf("psMatrix destructor ends\n");
#endif
}

// ************************************************************************
// get number of rows 
// ------------------------------------------------------------------------
int psMatrix::nrows()
{
   return nRows_;
}

// ************************************************************************
// get number of columns 
// ------------------------------------------------------------------------
int psMatrix::ncols()
{
   return nCols_;
}

// ************************************************************************
// load matrix from another matrix
// ------------------------------------------------------------------------
int psMatrix::load(psMatrix &inMat)
{
   int ii, jj;

#ifdef PS_DEBUG
   printf("psMatrix load\n");
#endif
   if (Mat_ != NULL)
   {
      for (ii = 0; ii < nRows_; ii++)
         if (Mat_[ii] != NULL) delete [] Mat_[ii];
      delete [] Mat_;
   }
   Mat_ = NULL;

   assert(this != &inMat);
   nRows_ = inMat.nrows();
   nCols_ = inMat.ncols();
   if (nRows_ > 0 && nCols_ > 0)
   {
      Mat_ = new double*[nRows_];
      assert(Mat_ != NULL);
      for (ii = 0; ii < nRows_; ii++)
      {
         Mat_[ii] = new double[nCols_];
         assert(Mat_[ii] != NULL);
         for (jj = 0; jj < nCols_; jj++) 
            Mat_[ii][jj] = inMat.getEntry(ii,jj);
      }
   }
   status_ = 0;
   determinant_ = inMat.determinant_;
#ifdef PS_DEBUG
   printf("psMatrix load ends\n");
#endif
   return 0;
}

// ************************************************************************
// load matrix from doubles
// ------------------------------------------------------------------------
int psMatrix::load(int nrows, int ncols, double **mat)
{
  int ii;
#ifdef PS_DEBUG
  printf("psMatrix load\n");
#endif
  if (Mat_ != NULL)
  {
    for (ii = 0; ii < nRows_; ii++)
      if (Mat_[ii] != NULL) delete [] Mat_[ii];
    delete [] Mat_;
  }
  Mat_ = NULL;

  assert(nrows);
  assert(ncols);
  assert(mat);
  nRows_ = nrows;
  nCols_ = ncols;
  Mat_   = mat;
  mat    = NULL;
  status_ = 0;
  determinant_ = 0.0;
#ifdef PS_DEBUG
  printf("psMatrix load ends\n");
#endif
  return 0;
}

// ************************************************************************
// load matrix from doubles
// ------------------------------------------------------------------------
int psMatrix::load(int nrows, int ncols, double *mat)
{
   int ii, jj;

#ifdef PS_DEBUG
   printf("psMatrix load\n");
#endif
   if (Mat_ != NULL)
   {
      for (ii = 0; ii < nRows_; ii++)
         if (Mat_[ii] != NULL) delete [] Mat_[ii];
      delete [] Mat_;
   }
   Mat_ = NULL;

   assert(nrows);
   assert(ncols);
   assert(mat);
   nRows_ = nrows;
   nCols_ = ncols;
   Mat_   = new double*[nrows];
   for (ii = 0; ii < nrows; ii++)
   {
     Mat_[ii] = new double[ncols];
     for (jj = 0; jj < ncols; jj++)
       Mat_[ii][jj] = mat[ii+jj*nrows];
   }
   status_ = 0;
   determinant_ = 0.0;
#ifdef PS_DEBUG
   printf("psMatrix load ends\n");
#endif
   return 0;
}

// ************************************************************************
// set matrix dimension
// ------------------------------------------------------------------------
int psMatrix::setDim(int nrows, int ncols)
{
   int ii, jj;

   if (Mat_ != NULL)
   {
      for (ii = 0; ii < nRows_; ii++)
         if (Mat_[ii] != NULL) delete [] Mat_[ii];
      delete [] Mat_;
   }
   Mat_ = NULL;

   nRows_ = nrows;
   nCols_ = ncols;
   if (nRows_ <= 0 || nCols_ <= 0) return -1;
   Mat_ = new double*[nRows_];
   assert(Mat_ != NULL);
   for (ii = 0; ii < nRows_; ii++)
   {
      Mat_[ii] = new double[nCols_];
      assert(Mat_[ii] != NULL);
      for (jj = 0; jj < nCols_; jj++) Mat_[ii][jj] = 0.0;
   }
   status_ = 0;
   return 0;
}

// ************************************************************************
// set entry
// ------------------------------------------------------------------------
void psMatrix::setEntry(const int row, const int col, const double ddata)
{
   if (row < 0 || row >= nRows_ || col < 0 || col >= nCols_)
   {
      printf("psMatrix setEntry ERROR: index (%d,%d) out of range (%d,%d)\n",
             row, col, nRows_, nCols_);
      exit(1);
   }
   Mat_[row][col] = ddata;
   status_ = 0;
}

// ************************************************************************
// get entry
// ------------------------------------------------------------------------
double psMatrix::getEntry(const int row, const int col)
{
   assert(row >= 0 && row < nRows_);
   assert(col >= 0 && col < nCols_);
#ifdef PS_DEBUG
   printf("psMatrix getEntry (%d,%d) : Mat(%d,%d) = %e\n",
          nRow_, nCols_, row, col, Mat_[row][col]);
#endif
   return Mat_[row][col];
}

// ************************************************************************
// get determinant 
// ------------------------------------------------------------------------
double psMatrix::getDeterminant()
{
   return determinant_;
}

// ************************************************************************
// get matrix in 1 dimensional array
// ------------------------------------------------------------------------
void psMatrix::getMatrix1D(psVector &mat)
{
   int ii, jj;
#ifdef PS_DEBUG
   printf("psMatrix getMatrix1D\n");
#endif
   assert(nRows_ >= 0);
   assert(nCols_ >= 0);
   mat.setLength(nRows_ * nCols_);
   for (ii = 0; ii < nRows_; ii++)
     for (jj = 0; jj < nCols_; jj++)
       mat[ii+jj*nRows_] = Mat_[ii][jj];
#ifdef PS_DEBUG
   printf("psMatrix getMatrix1D ends\n");
#endif
   return;
}

// ************************************************************************
// get matrix 
// ------------------------------------------------------------------------
double **psMatrix::getMatrix2D()
{
  return Mat_;
}

// ************************************************************************
// extract submatrix
// ------------------------------------------------------------------------
int psMatrix::submatrix(psMatrix &inMat, const int num, const int *indices)
{
   int nrows, ncols, ii, jj, row, col;

#ifdef PS_DEBUG
   printf("psMatrix submatrix\n");
#endif
   nrows = inMat.nrows();
   ncols = inMat.ncols();
   if (nrows != ncols)
   {
      printf("psMatrix::submatrix ERROR : incoming matrix is rectangular.\n");
      exit(1);
   }
   if (num <= 0 || indices == NULL)
   {
      printf("psMatrix::submatrix ERROR : no incoming indices.\n");
      exit(1);
   }
   for (ii = 0; ii < num; ii++)
   {
      row = indices[ii];
      if (row < 0 || row >= nrows)
      {
         printf("psMatrix::submatrix ERROR : index out of bound (%d)\n",row); 
         exit(1);
      }
   }
   setDim(num, num);
   for (ii = 0; ii < num; ii++)
   {
      row = indices[ii];
      for (jj = 0; jj < num; jj++)
      {
         col = indices[jj];
         Mat_[ii][jj] = inMat.getEntry(row, col);
      }
   }
#ifdef PS_DEBUG
   printf("psMatrix::submatrix: incoming matrix\n");
   for (ii = 0; ii < nrows; ii++)
   {
      for (jj = 0; jj < nrows; jj++) printf("%e ",inMat.getEntry(ii,jj));;
      printf("\n");
   }
   printf("psMatrix::submatrix: outgoing matrix\n");
   for (ii = 0; ii < num; ii++)
   {
      for (jj = 0; jj < num; jj++) printf("%e ",Mat_[ii][jj]);
      printf("\n");
   }
#endif
#ifdef PS_DEBUG
   printf("psMatrix submatrix ends\n");
#endif
   return 0;
}

#if 1
// ************************************************************************
// Cholesky decomposition (A = L L^T)
// ------------------------------------------------------------------------
int psMatrix::CholDecompose()
{
  int    ii, jj, status=0;
  double ddata, *mat;
  char   uplo='L';

#ifdef PS_DEBUG
  printf("psMatrix CholDecompose\n");
#endif
  if (status_ != 0) 
  {
     printf("psMatrix ERROR : matrix has been decomposed.\n");
     return -1;
  }
  assert(nRows_ == nCols_);

  mat = new double[nRows_*nRows_];
  for (ii = 0; ii < nRows_; ii++)
  {
    for (jj = 0; jj < nCols_; jj++)
      mat[ii+jj*nRows_] = Mat_[ii][jj];
  } 
  dpotrf_(&uplo, &nRows_, mat, &nRows_, &status);
  for (ii = 0; ii < nRows_; ii++)
  {
    for (jj = 0; jj < nCols_; jj++)
      Mat_[ii][jj] = mat[ii+jj*nRows_]; 
  } 
  if (status != 0)
  {
    printf("psMatrix ERROR (1): failed in Cholesky factorization.\n");
    //for (ii = 0; ii < nRows_; ii++)
    //{
    //  printf("Matrix: \n");
    //  for (jj = 0; jj < nCols_; jj++)
    //    printf("%12.4e ", Mat_[ii][jj]);
    //  printf("\n");
    //} 
    delete [] mat;
    return status;
  }
  delete [] mat;
  status_ = 1;
#ifdef PS_DEBUG
  printf("psMatrix CholDecompose ends\n");
#endif
  return 0;
}
#else
// ************************************************************************
// Cholesky decomposition (A = L L^T)
// ------------------------------------------------------------------------
int psMatrix::CholDecompose()
{
   int     ii, jj, kk;
   double  ddata;

#ifdef PS_DEBUG
   printf("psMatrix CholDecompose\n");
#endif
   if (status_ != 0) 
   {
      printf("psMatrix ERROR : matrix has been decomposed.\n");
      exit(1);
   }
   assert(nRows_ == nCols_);

#ifdef PS_DEBUG
   determinant_ = computeDeterminant(nRows_, Mat_); 
   printf("psMatrix determinant = %e\n", determinant_);
   for (ii = 0; ii < nRows_; ii++)
      for (jj = ii; jj < nCols_; jj++)
         printf("psMatrix (%d,%d) = %e\n", ii+1, jj+1, Mat_[ii][jj]);
#endif
   for (ii = 0; ii < nRows_; ii++)
   {
      for (jj = 0; jj <= ii; jj++)
      {
         ddata = Mat_[ii][jj];
         for (kk = 0; kk < jj; kk++) ddata -= Mat_[ii][kk] * Mat_[jj][kk];
         if (ii == jj)
         {
            if (ddata <= 0.0)
            {
               printf("CholDecompose : matrix not positive definite.\n");
               printf("dim = (%d,%d) : %e\n", ii+1, jj+1, ddata);
               return -1;
            }
            Mat_[ii][ii] = sqrt(ddata);
         }
         else Mat_[ii][jj] = ddata / Mat_[jj][jj];
#ifdef PS_DEBUG
         printf("psMatrix Chol (%d,%d) = %e\n", ii+1, jj+1, Mat_[jj][ii]);
#endif
      }
   } 
   for (ii = 0; ii < nRows_; ii++)
      for (jj = ii+1; jj < nCols_; jj++) Mat_[ii][jj] = Mat_[jj][ii];
   status_ = 1;
#ifdef PS_DEBUG
   printf("psMatrix CholDecompose ends\n");
#endif
   return 0;
}
#endif

// ************************************************************************
// matrix vector multiply (by the L factor)  
// ------------------------------------------------------------------------
int psMatrix::CholLMatvec(psVector &ivec, psVector &ovec)
{
  int    ii, jj, status=0;
  double ddata;

#ifdef PS_DEBUG
  printf("psMatrix CholMatvec\n");
#endif
  assert(ivec.length() == nCols_);
  if (status_ == 0) status = CholDecompose();
  if (status != 0)
  {
    printf("psMatrix ERROR (2): failed in Cholesky factorization.\n");
    return status;
  }
  ovec.setLength(nRows_);
  for (ii = nRows_-1; ii >= 0; ii--)
  {
    ddata = 0.0;
    for (jj = 0; jj <= ii; jj++) ddata += Mat_[ii][jj] * ivec[jj];
    ovec[ii] = ddata;
  }
#ifdef PS_DEBUG
  printf("psMatrix CholMatvec ends\n");
#endif
  return 0;
}

// ************************************************************************
// Cholesky L-solve 
// ------------------------------------------------------------------------
int psMatrix::CholSolve(psVector &ivec, psVector &ovec)
{
  int    ii, jj, iOne=1, status=0;
  double *vec, *mat;
  char   uplo='L';

#ifdef PS_DEBUG
  printf("psMatrix CholSolve\n");
#endif
  assert(ivec.length() == nCols_);
  if (status_ == 0) status = CholDecompose();
  if (status != 0)
  {
    printf("psMatrix ERROR (3): failed in Cholesky factorization.\n");
    return status;
  }
  ovec.setLength(nRows_);
  for (ii = 0; ii < nRows_; ii++) ovec[ii] = ivec[ii];
  vec = ovec.getDVector();
  mat = new double[nRows_*nRows_];
  for (ii = 0; ii < nRows_; ii++)
  {
    for (jj = 0; jj < nCols_; jj++)
      mat[ii+jj*nRows_] = Mat_[ii][jj];
  } 
  dpotrs_(&uplo,&nRows_,&iOne, mat, &nRows_,vec, &nRows_, &status);
  if (status != 0)
  {
    printf("psMatrix ERROR (1): failed in Cholesky solve.\n");
    return status;
  }
  delete [] mat;
#ifdef PS_DEBUG
   printf("psMatrix CholSolve ends\n");
#endif
  return 0;
}

// ************************************************************************
// Cholesky L-solve 
// ------------------------------------------------------------------------
int psMatrix::CholLSolve(psVector &ivec, psVector &ovec)
{
  int    ii, jj, status=0;
  double ddata;

#ifdef PS_DEBUG
  printf("psMatrix CholLSolve\n");
#endif
  assert(ivec.length() == nCols_);
  if (status_ == 0) status = CholDecompose();
  if (status != 0)
  {
    printf("psMatrix ERROR (4): failed in Cholesky factorization.\n");
    return status;
  }
  ovec.setLength(nRows_);
  for (ii = 0; ii < nRows_; ii++)
  {
    ddata = ivec[ii];
    for (jj = 0; jj < ii; jj++) ddata -= Mat_[ii][jj] * ovec[jj];
    ovec[ii] = ddata / Mat_[ii][ii];
  }
#ifdef PS_DEBUG
  printf("psMatrix CholLSolve ends\n");
#endif
  return 0;
}

// ************************************************************************
// Cholesky LT-solve 
// ------------------------------------------------------------------------
int psMatrix::CholLTSolve(psVector &ivec, psVector &ovec)
{
  int    ii, jj, status=0;
  double ddata;

#ifdef PS_DEBUG
  printf("psMatrix CholTSolve (transpose)\n");
#endif
  assert(ivec.length() == nCols_);
  if (status_ == 0) status = CholDecompose();
  if (status != 0)
  {
    printf("psMatrix ERROR (1): failed in Cholesky factorization.\n");
    return status;
  }
  ovec.setLength(nRows_);
  for (ii = nRows_-1; ii > 0; ii--)
  {
    ddata = ivec[ii];
    for (jj = ii+1; jj < nRows_; jj++) ddata -= Mat_[jj][ii] * ovec[jj];
    ovec[ii] = ddata / Mat_[ii][ii];
  }
#ifdef PS_DEBUG
  printf("psMatrix CholTSolve ends\n");
#endif
  return 0;
}

// ************************************************************************
// Compute LU factorization 
// ------------------------------------------------------------------------
int psMatrix::LUDecompose()
{
  int    ii, jj, lwork, status=0;
  double *localMatrix, *work;

  assert(nRows_ == nCols_);
  if (pivots_ != NULL) pivots_ = new int[nRows_];
  pivots_ = new int[nRows_];
  localMatrix = new double[nRows_ * nRows_];
  work = new double[nRows_ * nRows_];
  for (ii = 0; ii < nRows_; ii++)
    for (jj = 0; jj < nRows_; jj++)
      localMatrix[ii*nRows_+jj] = Mat_[jj][ii];
  lwork = nRows_ * nRows_;
  dgetrf_(&nRows_, &nRows_, localMatrix, &nRows_, pivots_, &status);
  for (ii = 0; ii < nRows_; ii++)
    for (jj = 0; jj < nRows_; jj++)
      Mat_[ii][jj] = localMatrix[jj*nRows_+ii];
  delete [] localMatrix;
  delete [] work;
  if (status != 0)
  {
    printf("psMatrix computeInverse ERROR: failed in LU factorization\n");
    delete [] pivots_;
    pivots_ = NULL;
  }
  else status_ = 2;
  return status;
}

// ************************************************************************
// LU solve 
// ------------------------------------------------------------------------
int psMatrix::LUSolve(psVector &ivec, psVector &ovec)
{
  int    ii, jj, iOne=1, status=0;
  double *vec, *mat;
  char   trans='N';

#ifdef PS_DEBUG
  printf("psMatrix LUSolve\n");
#endif
  assert(ivec.length() == nCols_);
  if (status_ != 2)
  {
    printf("psMatrix ERROR (3): LU factorization has not been called.\n");
    status = -1;
    return status;
  }
  ovec.setLength(nRows_);
  for (ii = 0; ii < nRows_; ii++) ovec[ii] = ivec[ii];
  vec = ovec.getDVector();
  mat = new double[nRows_*nRows_];
  for (ii = 0; ii < nRows_; ii++)
  {
    for (jj = 0; jj < nCols_; jj++)
      mat[ii+jj*nRows_] = Mat_[ii][jj];
  } 
  dgetrs_(&trans,&nRows_,&iOne,mat,&nRows_,pivots_,vec,&nRows_,&status);
  if (status != 0)
  {
    printf("psMatrix ERROR (1): failed in LU solve.\n");
    delete [] mat;
    return status;
  }
  delete [] mat;
#ifdef PS_DEBUG
   printf("psMatrix LUSolve ends\n");
#endif
  return 0;
}

// ************************************************************************
// print matrix
// ------------------------------------------------------------------------
void psMatrix::print()
{
   int ii, jj;
   printf("psMatrix print: \n");
   for (ii = 0; ii < nRows_; ii++)
   {
      for (jj = 0; jj < nCols_; jj++) printf("%e ", Mat_[ii][jj]);
      printf("\n");
   }
}

// ************************************************************************
// Compute determinant (by Bourke)
// ------------------------------------------------------------------------
double psMatrix::computeDeterminant()
{
   int    ii, jj, kk, ind;
   double result = 0.0;
   double **localMat = NULL;

   assert(nRows_ == nCols_);
   if (status_ == 1)
   {
     result = 1.0;
     for (ii = 0; ii < nRows_; ii++) result *= Mat_[ii][ii];
     return result;
   }
   if (nRows_ == 1)
   {
      result = Mat_[0][0];
   }
   else if (nRows_ == 2)
   {
      result = Mat_[0][0] * Mat_[1][1] - Mat_[1][0] * Mat_[0][1];
   }
   else
   {
      result = 0.0;
      for (ii = 0; ii < nRows_; ii++)
      {
         localMat = new double*[nRows_-1];
         for (kk = 0; kk < nRows_-1; kk++)
            localMat[kk] = new double[nRows_-1];
         for (kk = 1; kk < nRows_; kk++)
         {
            ind = 0;
            for (jj = 0; jj < nRows_; jj++)
            {
               if (jj == ii) continue;
               localMat[kk-1][ind] = Mat_[kk][jj];
               ind++;
            }
         }
         result += pow(-1.0,1.0+ii+1.0) * Mat_[0][ii] * 
                   computeDeterminant(nRows_-1, localMat);
         for (kk = 0; kk < nRows_-1; kk++) delete [] localMat[kk];
         delete [] localMat;
      }
   }
   return result;
}

// ************************************************************************
// Compute determinant (by Bourke)
// ------------------------------------------------------------------------
double psMatrix::computeDeterminant(int ndim, double **mat)
{
   int    ii, jj, kk, ind;
   double result = 0.0;
   double **localMat = NULL;

   if (ndim == 1)
   {
      result = mat[0][0];
   }
   else if (ndim == 2)
   {
      result = mat[0][0] * mat[1][1] - mat[1][0] * mat[0][1];
   }
   else
   {
      result = 0.0;
      for (ii = 0; ii < ndim; ii++)
      {
         localMat = new double*[ndim-1];
         for (kk = 0; kk < ndim-1; kk++)
            localMat[kk] = new double[ndim-1];
         for (kk = 1; kk < ndim; kk++)
         {
            ind = 0;
            for (jj = 0; jj < ndim; jj++)
            {
               if (jj == ii) continue;
               localMat[kk-1][ind] = mat[kk][jj];
               ind++;
            }
         }
         result += pow(-1.0,1.0+ii+1.0) * mat[0][ii] * 
                   computeDeterminant(ndim-1, localMat);
         for (kk = 0; kk < ndim-1; kk++) delete [] localMat[kk];
         delete [] localMat;
      }
   }
   return(result);
}

// ************************************************************************
// Compute inverse 
// ------------------------------------------------------------------------
int psMatrix::computeInverse(psMatrix &inverse)
{
  int    ii, jj, *ipiv, lwork, status;
  double *localMatrix, *work;

  assert(nRows_ == nCols_);
  ipiv = new int[nRows_];
  localMatrix = new double[nRows_ * nRows_];
  work = new double[nRows_ * nRows_];
  for (ii = 0; ii < nRows_; ii++)
    for (jj = 0; jj < nRows_; jj++)
      localMatrix[ii*nRows_+jj] = Mat_[jj][ii];
  lwork = nRows_ * nRows_;
  dgetrf_(&nRows_, &nRows_, localMatrix, &nRows_, ipiv, &status);
  if (status != 0)
  {
    printf("psMatrix computeInverse ERROR: failed in LU factorization\n");
    return status;
  }
  dgetri_(&nRows_, localMatrix, &nRows_, ipiv, work, &lwork, &status);
  if (status != 0)
  {
    printf("psMatrix computeInverse ERROR: failed in matrix inverse\n");
    return status;
  }
  inverse.setDim(nRows_, nRows_);
  for (ii = 0; ii < nRows_; ii++)
    for (jj = 0; jj < nRows_; jj++)
      inverse.setEntry(ii, jj, localMatrix[jj*nRows_+ii]);
  delete [] work;
  delete [] ipiv;
  delete [] localMatrix;
  return status;
}

// ************************************************************************
// matrix vector multiply
// ------------------------------------------------------------------------
int psMatrix::matvec(psVector &inVec, psVector &outVec, int transp)
{
  int    ii, jj;
  double ddata, *vdata;

  if (transp == 0)
  {
    assert(inVec.length() == nCols_);
    outVec.setLength(nRows_);
    vdata = inVec.getDVector();
    for (ii = 0; ii < nRows_; ii++)
    {
      ddata = 0.0;
      for (jj = 0; jj < nCols_; jj++)
        ddata += Mat_[ii][jj] * vdata[jj];
      outVec[ii] = ddata;
    }
  }
  else
  {
    assert(inVec.length() == nRows_);
    outVec.setLength(nCols_);
    vdata = inVec.getDVector();
    for (ii = 0; ii < nCols_; ii++)
    {
      ddata = 0.0;
      for (jj = 0; jj < nRows_; jj++)
        ddata += Mat_[jj][ii] * vdata[jj];
      outVec[ii] = ddata;
    }
  }
  return 0;
}

// ************************************************************************
// matrix matrix multiply
// ------------------------------------------------------------------------
void psMatrix::matmult(psMatrix &inMat, psMatrix &outMat)
{
  int      ii, jj, kk, ncols;
  double   ddata, *vdata, **idata, **odata;
  psVector colVec;

  assert(inMat.nrows() == nCols_);
  outMat.setDim(nRows_, inMat.ncols());
  idata = inMat.getMatrix2D();
  odata = outMat.getMatrix2D();
  ncols = inMat.ncols();
  for (jj = 0; jj < ncols; jj++)
  {
    for (ii = 0; ii < nRows_; ii++)
    {
      ddata = 0.0;
      for (kk = 0; kk < nCols_; kk++)
        ddata += Mat_[kk][ii] * idata[jj][kk];
      odata[ii][jj] = ddata;
    }
  }
}

// ************************************************************************
// matrix transpose 
// ------------------------------------------------------------------------
void psMatrix::transpose()
{
  int    ii, jj;
  double **tmpMat;

  assert(nCols_ > 0 && nRows_ > 0);
  tmpMat = new double*[nRows_];
  for (ii = 0; ii < nRows_; ii++)
  {
     tmpMat[ii] = new double[nCols_];
     for (jj = 0; jj < nCols_; jj++) tmpMat[ii][jj] = Mat_[jj][ii];
  }
  for (ii = 0; ii < nCols_; ii++) delete [] Mat_[ii];
  delete [] Mat_;
  Mat_ = tmpMat;
  ii = nRows_;
  nRows_ = nCols_;
  nCols_ = ii;
}

// ************************************************************************
// eigen solve 
// ------------------------------------------------------------------------
void psMatrix::eigenSolve(psMatrix &eigMat, psVector &eigVals, int flag)
{
   int    ii, jj, lwork, N, info;
   double *work, *eigs, *mat;
   char   jobz='V', uplo='U';

   if (flag == 1) jobz = 'N';
   N     = nCols_;
   lwork = 3 * N;
   work  = new double[lwork];
   eigs  = new double[N];
   mat   = new double[N*N];
   for (ii = 0; ii < N; ii++)
      for (jj = 0; jj < N; jj++)
         mat[ii*N+jj] = Mat_[ii][jj];
   dsyev_(&jobz,&uplo,&N,mat,&N,eigs,work,&lwork,&info);
   if (info != 0)
   {
      printf("ERROR: dsyev returns a nonzero (%d).\n", info);
      delete [] mat;
      delete [] eigs;
      delete [] work;
      exit(1);
   }

   eigMat.setDim(N,N);
   for (ii = 0; ii < N; ii++)
      for (jj = 0; jj < N; jj++)
         eigMat.setEntry(jj,ii,mat[ii*N+jj]);
   eigVals.setLength(N);
   for (ii = 0; ii < N; ii++) eigVals[ii] = eigs[ii];
   delete [] mat;
   delete [] eigs;
   delete [] work;
}

// ************************************************************************
// matrix solve using QR 
// ------------------------------------------------------------------------
void psMatrix::matSolve(psVector &invec, psVector &outvec)
{
   int    iOne=1, info, lwork, nn, ii, jj;
   double *work, *b, *x, *dmat;
   char   trans[1];

   nn = invec.length();
   (*trans) = 'N';
   b = invec.getDVector();
   outvec.setLength(nn);
   x = outvec.getDVector();
   for (ii = 0; ii < nn; ii++) x[ii] = b[ii];
   lwork = 2 * nn * nn;
   work = new double[lwork];
   dmat = new double[nn*nn];
   for (ii = 0; ii < nn; ii++)
      for (jj = 0; jj < nn; jj++)
         dmat[ii*nn+jj] = Mat_[ii][jj];
   dgels_(trans, &nn, &nn, &iOne, dmat, &nn, x, &nn, work, &lwork, &info);
   if (info != 0)
   {
      printf("psMatrix matSolve ERROR: dgels returns error %d.\n",info);
      exit(1);
   }
   delete [] dmat;
   delete [] work;
   return;
}

