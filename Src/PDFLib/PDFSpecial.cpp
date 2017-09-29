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
// Functions for the user-defined distribution
// AUTHOR : CHARLES TONG
// DATE   : 2013
// ************************************************************************
#include <stdio.h>
#include <math.h>
#include "PsuadeUtil.h"
#include "PDFSpecial.h"
#define PABS(x) (((x) >= 0) ? x : -(x))

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
PDFSpecial::PDFSpecial(int scount)
{
   int  ii, jj, nn;
   char pString[1000], filename[1001];
   FILE *fp=NULL;

   printf("PDFSpecial constructor: expecting a sample file.\n");
   printf("                        having the following format: \n");
   printf("line 1: <number of sample points> <number of inputs>\n");
   printf("line 2: 1 sample point 1 inputs \n");
   printf("line 3: 2 sample point 2 inputs \n");
   printf("line 4: 3 sample point 3 inputs \n");
   printf("...\n");
   sprintf(pString,"Enter name of sample file : ");
   getString(pString, filename);
   nn = strlen(filename);
   if (nn > 1000)
   {
      printf("PDFSpecial constructor ERROR: file name too long.\n");
      exit(1);
   }
   filename[nn-1] = '\0';
   nSamples_ = 0;
   samples_ = NULL;
   fp = fopen(filename, "r");
   if (fp == NULL)
   {
      printf("PDFSpecial constructor ERROR: cannot open sample file %s.\n",
             filename);
      exit(1);
   }
   fscanf(fp, "%d %d", &nSamples_, &nInputs_);
   if (nSamples_ < 1 || nInputs_ < 1)
   {
      printf("PDFSpecial constructor ERROR: nSamples or nInputs <= 0.\n");
      exit(1);
   }
   if (nInputs_ != scount)
   {
      printf("PDFSpecial constructor ERROR: nInputs does not match.\n");
      printf("               nInputs in your sample file    = %d\n", nInputs_);
      printf("               nInputs from psuade input file = %d\n", scount);
      exit(1);
   }
   samples_ = new double[nSamples_*nInputs_];
   for (ii = 0; ii < nSamples_; ii++)
   {
      fscanf(fp, "%d", &nn);
      if (nn != (ii+1))
      {
         printf("PDFSpecial constructor ERROR: invalid sample number.\n");
         printf("                    Expected: %d\n", ii+1);
         printf("                    Read:     %d\n", nn);
         printf("Advice: check your line number %d.\n\n",ii+2);
         printf("Correct Format: \n");
         printf("line 1: <number of sample points> <number of inputs>\n");
         printf("line 2: 1 sample point 1 inputs \n");
         printf("line 3: 2 sample point 2 inputs \n");
         printf("line 4: 3 sample point 3 inputs \n");
         printf("...\n");
         fclose(fp);
         exit(1);
      } 
      for (jj = 0; jj < nInputs_; jj++)
         fscanf(fp, "%lg", &samples_[ii*nInputs_+jj]);
   }
   fclose(fp);
   printf("PDFSpecial constructor INFO: sample file read.\n");
   printf("           sample size   = %d\n", nSamples_);
   printf("           no. of inputs = %d\n", nInputs_);
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
PDFSpecial::~PDFSpecial()
{
   if (samples_ != NULL) delete [] samples_;
}

// ************************************************************************
// forward transformation to range
// ------------------------------------------------------------------------
int PDFSpecial::getPDF(int length, double *inData, double *outData)
{
   printf("PDFSpecial::getPDF not available.\n");
   for (int ii = 0; ii < length; ii++) outData[ii] = 0;
   return -1;
}

// ************************************************************************
// look up cumulative density
// ------------------------------------------------------------------------
int PDFSpecial::getCDF(int length, double *inData, double *outData)
{
   printf("PDFSpecial::getCDF not available.\n");
   for (int ii = 0; ii < length; ii++) outData[ii] = 0;
   return -1;
}

// ************************************************************************
// transformation to range
// ------------------------------------------------------------------------
int PDFSpecial::invCDF(int length, double *inData, double *outData,
                       double lower, double upper)
{
   printf("PDFSpecial::invCDF not available.\n");
   for (int ii = 0; ii < length; ii++) outData[ii] = 0;
   return -1;
}

// ************************************************************************
// generate a sample
// ------------------------------------------------------------------------
int PDFSpecial::genSample(int length,double *outData,double,double)
{
   int    ii, jj, ind;

   for (ii = 0; ii < length; ii++)
   {
      ind = PSUADE_rand() % nSamples_;
      for (jj = 0; jj < nInputs_; jj++)
         outData[ii*nInputs_+jj] = samples_[ind*nInputs_+jj];
   }
   return 0;
}

// ************************************************************************
// get mean
// ------------------------------------------------------------------------
double PDFSpecial::getMean()
{
   printf("PDFSpecial::getMean not available for this distribution.\n");
   return 0.0;
}

