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
PDFSpecial::PDFSpecial(int scount, char *fname, int *indices)
{
   int    ii, jj, nn, nInps;
   double *oneSample;
   char   pString[1001], filename[1001];
   FILE   *fp=NULL;

   if (fname == NULL || !strcmp(fname, "NONE"))
   {
      printf("PDFSpecial constructor: expecting a sample file.\n");
      printf("                        having the following format: \n");
      printf("line 1: (optional) PSUADE_BEGIN\n");
      printf("line 2: <number of sample points> <number of inputs>\n");
      printf("line 3: (optional) : '#' followed by input names\n");
      printf("line 4: 1 sample point 1 inputs \n");
      printf("line 5: 2 sample point 2 inputs \n");
      printf("line 6: 3 sample point 3 inputs \n");
      printf("...\n");
      printf("last n: (optional) PSUADE_END\n");
      sprintf(pString,"Enter name of sample file : ");
      getString(pString, filename);
      nn = strlen(filename);
      if (nn > 1000)
      {
         printf("PDFSpecial constructor ERROR: file name too long.\n");
         exit(1);
      }
      filename[nn-1] = '\0';
   }
   else strcpy(filename, fname);

   nSamples_ = 0;
   samples_ = NULL;
   nInputs_ = scount;
   fp = fopen(filename, "r");
   if (fp == NULL)
   {
      printf("PDFSample ERROR: cannot open sample file %s.\n",filename);
      exit(1);
   }
   fscanf(fp, "%s", pString);
   if (strcmp(pString, "PSUADE_BEGIN"))
   {
      fclose(fp);
      fp = fopen(filename, "r");
   } 
   fscanf(fp, "%d %d", &nSamples_, &nInps);
   if (nSamples_ < 1)
   {
      printf("PDFSample ERROR: sample file has nSamples <= 0.\n");
      exit(1);
   }
   if (nInps < 1)
   {
      printf("PDFSample ERROR: sample file has nInputs <= 0.\n");
      exit(1);
   }
   if (nInputs_ != nInps && indices == NULL)
   {
      printf("PDFSample ERROR: nInputs does not match.\n");
      printf("          nInputs in your sample file    = %d\n",nInps);
      printf("          nInputs from psuade input file = %d\n",nInputs_);
      exit(1);
   }
   if (indices != NULL)
   {
      for (ii = 0; ii < scount; ii++)
      {
         if (indices[ii] < 0 || indices[ii] >= nInps)
         {
            printf("PDFSample ERROR: sample index > nInputs in sample file.\n");
            printf("          sample index requested         = %d\n",
                   indices[ii]+1);
            printf("          nInputs in your sample file    = %d\n",nInps);
            exit(1);
         } 
      }
   }
   fgets(pString, 1000, fp);
   while (1)
   {
      nn = getc(fp);
      if (nn == '#') fgets(pString, 1000, fp);
      else
      {
         ungetc(nn, fp);
         break;
      }
   }
   samples_  = new double[nSamples_*nInputs_];
   oneSample = new double[nInps];
   for (ii = 0; ii < nSamples_; ii++)
   {
      fscanf(fp, "%d", &nn);
      if (nn != (ii+1))
      {
         printf("PDFSample ERROR: invalid sample number.\n");
         printf("          Expected: %d\n", ii+1);
         printf("          Read:     %d\n", nn);
         printf("Advice: check your data format and line number %d.\n\n",ii+2);
         printf("Correct Format: \n");
         printf("line 1: (optional) PSUADE_BEGIN\n");
         printf("line 2: <number of sample points> <number of inputs>\n");
         printf("line 3: (optional) : '#' followed by input names\n");
         printf("line 4: 1 sample point 1 inputs \n");
         printf("line 5: 2 sample point 2 inputs \n");
         printf("line 6: 3 sample point 3 inputs \n");
         printf("...\n");
         printf("last n: (optional) PSUADE_END\n");
         fclose(fp);
         exit(1);
      } 
      for (jj = 0; jj < nInps; jj++)
         fscanf(fp, "%lg", &oneSample[jj]);
      for (jj = 0; jj < nInputs_; jj++)
      {
         nn = indices[jj];
         samples_[ii*nInputs_+jj] = oneSample[nn] ;
      }
      fgets(pString, 1000, fp);
   }
   fclose(fp);
   delete [] oneSample;
   printf("PDFSample INFO: sample file '%s' has been read.\n", fname);
   printf("   sample size   = %d\n", nSamples_);
   printf("   no. of inputs = %d\n", nInputs_);
   if (indices != NULL)
   {
      for (ii = 0; ii < nInputs_; ii++)
         printf("   Input %d has column %d in the sample file.\n",
                ii+1, indices[ii]+1);
   }
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
   printf("PDFSample::getPDF not available.\n");
   for (int ii = 0; ii < length; ii++) outData[ii] = 0;
   return -1;
}

// ************************************************************************
// look up cumulative density
// ------------------------------------------------------------------------
int PDFSpecial::getCDF(int length, double *inData, double *outData)
{
   printf("PDFSample::getCDF not available.\n");
   for (int ii = 0; ii < length; ii++) outData[ii] = 0;
   return -1;
}

// ************************************************************************
// transformation to range
// ------------------------------------------------------------------------
int PDFSpecial::invCDF(int length, double *inData, double *outData,
                       double lower, double upper)
{
   printf("PDFSample::invCDF not available.\n");
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
   printf("PDFSample::getMean not available for this distribution.\n");
   return 0.0;
}

// ************************************************************************
// get mean
// ------------------------------------------------------------------------
int PDFSpecial::setParam(char *sparam)
{
   char winput[1001];
   sscanf("%s", winput);
   if (!strcmp(winput, "setIput"))
   {
      sscanf("%s %d", winput, &whichInput_);
      if (whichInput_ < 0) whichInput_ = 0;
   }
   return 0;
}

