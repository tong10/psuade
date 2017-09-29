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
// Functions for the histogram-based distribution
// AUTHOR : CHARLES TONG
// DATE   : 2015
// ************************************************************************
#include <stdio.h>
#include <math.h>
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "PDFHistogram.h"
#include "PrintingTS.h"
#define PABS(x) (((x) >= 0) ? x : -(x))

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
PDFHistogram::PDFHistogram(char *fname, int nInps, int *indices, int *incrs)
{
   int    ii, jj, ss, nn;
   double *oneSample, ddata, dmin, dmax;
   char   pString[1001], filename[1001];
   FILE   *fp=NULL;

   if (fname == NULL || !strcmp(fname, "NONE"))
   {
      printf("PDFHistogram constructor: expecting a sample file.\n");
      printf("                          having the following format: \n");
      printf("line 1: (optional) PSUADE_BEGIN\n");
      printf("line 2: <number of sample points> <number of inputs>\n");
      printf("line 3: (optional) : '#' followed by input names\n");
      printf("line 4: 1 sample point 1 inputs \n");
      printf("line 5: 2 sample point 2 inputs \n");
      printf("line 6: 3 sample point 3 inputs \n");
      printf("...\n");
      printf("line n: (optional) PSUADE_END\n");
      sprintf(pString,"Enter name of sample file : ");
      getString(pString, filename);
      nn = strlen(filename);
      if (nn > 1000)
      {
         printf("PDFHistogram constructor ERROR: file name too long.\n");
         exit(1);
      }
      filename[nn-1] = '\0';
   }
   else strcpy(filename, fname);

   nSamples_ = 0;
   samples_  = NULL;
   nInputs_  = nInps;
   fp = fopen(filename, "r");
   if (fp == NULL)
   {
      printf("PDFHistogram ERROR: cannot open sample file %s.\n",
             filename);
      exit(1);
   }
   fscanf(fp, "%s", pString);
   if (strcmp(pString, "PSUADE_BEGIN"))
   {
      fclose(fp);
      fp = fopen(filename, "r");
   } 
   fscanf(fp, "%d %d", &nSamples_, &nn);
   if (nSamples_ < 10000)
   {
      printf("PDFHistogram ERROR: sample file needs to be >= 10000.\n");
      fclose(fp);
      exit(1);
   }
   if (nn < 1 || nn > 15)
   {
      printf("PDFHistogram ERROR: sample file has nInputs <= 0 or > 15.\n");
      fclose(fp);
      exit(1);
   }
   if (nInputs_ != nn)
   {
      printf("PDFHistogram ERROR: nInputs does not match.\n");
      printf("          nInputs in your sample file    = %d\n",nn);
      printf("          nInputs from psuade input file = %d\n",nInputs_);
      fclose(fp);
      exit(1);
   }
   if (indices != NULL)
   {
      for (ii = 0; ii < nInputs_; ii++)
      {
         if (indices[ii] < 0 || indices[ii] >= nInps)
         {
            printf("PDFHistogram ERROR: wrong sample index.\n");
            printf("             sample index requested         = %d\n",
                   indices[ii]+1);
            fclose(fp);
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
         printf("PDFHistogram ERROR: invalid sample number.\n");
         printf("             Expected: %d\n", ii+1);
         printf("             Read:     %d\n", nn);
         printf("Advice: check your data format and line number %d.\n\n",ii+2);
         printf("Correct Format: \n");
         printf("line 1: (optional) PSUADE_BEGIN\n");
         printf("line 2: <number of sample points> <number of inputs>\n");
         printf("line 3: (optional) : '#' followed by input names\n");
         printf("line 4: 1 sample point 1 inputs \n");
         printf("line 5: 2 sample point 2 inputs \n");
         printf("line 6: 3 sample point 3 inputs \n");
         printf("...\n");
         printf("line n: (optional) PSUADE_END\n");
         fclose(fp);
         exit(1);
      } 
      for (jj = 0; jj < nInps; jj++)
         fscanf(fp, "%lg", &oneSample[jj]);
      for (jj = 0; jj < nInputs_; jj++)
      {
         if (indices != NULL) nn = indices[jj];
         else                 nn = jj;
         samples_[ii*nInputs_+jj] = oneSample[nn] ;
      }
      fgets(pString, 1000, fp);
   }
   fscanf(fp, "%s", pString);
   fclose(fp);
   delete [] oneSample;
   printOutTS(PL_INFO,"PDFHistogram INFO: sample file '%s' has been read.\n", 
              fname);
   printOutTS(PL_INFO,"   Sample size   = %d\n", nSamples_);
   printOutTS(PL_INFO,"   No. of inputs = %d\n", nInputs_);
   if (indices != NULL)
   {
      for (ii = 0; ii < nInputs_; ii++)
         printOutTS(PL_INFO,"   Input %d has column %d in the sample file.\n",
                    ii+1, indices[ii]+1);
   }

   lowerBs_ = new double[nInputs_];
   upperBs_ = new double[nInputs_];
   for (ii = 0; ii < nInputs_; ii++)
   {
      dmin = dmax = samples_[ii];
      for (nn = 1; nn < nSamples_; nn++)
      {
         ddata = samples_[nn*nInputs_+ii];
         if (ddata < dmin) dmin = ddata;
         if (ddata > dmax) dmax = ddata;
      }
      lowerBs_[ii] = dmin - 0.01 * (dmax - dmin);
      upperBs_[ii] = dmax + 0.01 * (dmax - dmin);
      if (lowerBs_[ii] == upperBs_[ii])
      {
         printf("PDFHistogram ERROR: upper bound=lower bound for input %d.\n",
                ii+1);
         exit(1);
      }
   }

   incrs_ = new int[nInps];
   for (ii = 0; ii < nInputs_; ii++)
   {
       if (incrs[ii] <= 1)
       {
          printf("PDFHistogram ERROR: invalid partition.\n");
          exit(1);
       }
       incrs_[ii] = incrs[ii];
   }
   
   initHistogram();
   for (ss = 0; ss < nSamples_; ss++)
   {
      for (ii = 0; ii < nInputs_; ii++) 
      {
         ddata = samples_[ss*nInputs_+ii];
         ddata = (ddata - lowerBs_[ii]) / (upperBs_[ii]-lowerBs_[ii]);
         if (ddata == 1.0) jj = incrs_[ii] - 1;
         else              jj = (int) (ddata * incrs_[ii]);
         indexSet_[ii] = jj;
      }
      mergeHistogram(indexSet_,ss);
   }
   finalizeHistogram();
}

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
PDFHistogram::PDFHistogram(int nSamp,int nInps,double *samInputs,
                           int *incrs, int flag)
{
   int    ii, jj, ss;
   double ddata, dmin, dmax;
   char   pString[1001];

   if (nSamp < 10000)
   {
      printf("PDFHistogram ERROR: sample file needs to be >= 10000.\n");
      exit(1);
   }
   if (nInps < 1 || nInps > 15)
   {
      printf("PDFHistogram ERROR: sample file has nInputs <= 0 or > 15.\n");
      exit(1);
   }
   nSamples_ = nSamp;
   nInputs_  = nInps;
   samples_  = new double[nSamples_*nInputs_];
   for (ii = 0; ii < nSamp*nInps; ii++) samples_[ii] = samInputs[ii];
   printOutTS(PL_INFO,"PDFHistogram  Sample size   = %d\n", nSamples_);
   printOutTS(PL_INFO,"PDFHistogram  No. of inputs = %d\n", nInputs_);

   lowerBs_ = new double[nInputs_];
   upperBs_ = new double[nInputs_];
   for (ii = 0; ii < nInputs_; ii++)
   {
      dmin = dmax = samples_[ii];
      for (ss = 1; ss < nSamples_; ss++)
      {
         ddata = samples_[ss*nInputs_+ii];
         if (ddata < dmin) dmin = ddata;
         if (ddata > dmax) dmax = ddata;
      }
      lowerBs_[ii] = dmin - 0.01 * (dmax - dmin);
      upperBs_[ii] = dmax + 0.01 * (dmax - dmin);
      if (lowerBs_[ii] == upperBs_[ii])
      {
         printf("PDFHistogram ERROR: upper bound=lower bound for input %d.\n",
                ii+1);
         exit(1);
      }
   }

   incrs_ = new int[nInps];
   for (ii = 0; ii < nInputs_; ii++)
   {
       if (incrs[ii] <= 1)
       {
          printf("PDFHistogram ERROR: invalid partition.\n");
          exit(1);
       }
       incrs_[ii] = incrs[ii];
   }
   
   initHistogram();
   for (ss = 0; ss < nSamples_; ss++)
   {
      for (ii = 0; ii < nInputs_; ii++) 
      {
         ddata = samples_[ss*nInputs_+ii];
         ddata = (ddata - lowerBs_[ii]) / (upperBs_[ii]-lowerBs_[ii]);
         if (ddata == 1.0) jj = incrs_[ii] - 1;
         else              jj = (int) (ddata * incrs_[ii]);
         indexSet_[ii] = jj;
      }
      mergeHistogram(indexSet_,ss);
   }
   finalizeHistogram();
   if (flag == 1)
   {
      FILE *fp = fopen("psuade_pdfhist_sample","w");
      if (fp == NULL)
      {
         printf("PDFHistogram ERROR: no histogram file created.\n");
         return;
      }
      fprintf(fp, "%d %d\n", nHist_, nInputs_);
      for (ii = 0; ii < nHist_; ii++)
      {
         for (jj = 0; jj < nInputs_; jj++)
            fprintf(fp, "%16.8e ", histMedians_[ii*nInputs_+jj]);
         fprintf(fp, "%16.8e\n", 1.0*histCnts_[ii]/nSamples_);
      }      
      fclose(fp);
   }
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
PDFHistogram::~PDFHistogram()
{
   if (samples_   != NULL) delete [] samples_;
   if (lowerBs_   != NULL) delete [] lowerBs_;
   if (upperBs_   != NULL) delete [] upperBs_;
   if (incrs_     != NULL) delete [] incrs_;
   if (histCells_ != NULL) delete [] histCells_;
   if (indexSet_  != NULL) delete [] indexSet_;
   if (histMedians_   != NULL) delete [] histMedians_;
   if (histCDF_   != NULL) delete [] histCDF_;
   if (histCnts_  != NULL) delete [] histCnts_;
   if (histMap_   != NULL)
   {
      for (int ii = 0; ii < nHist_; ii++) delete [] histMap_[ii];
      delete [] histMap_;
   }
}

// ************************************************************************
// forward transformation to range
// ------------------------------------------------------------------------
int PDFHistogram::getPDF(int length, double *inData, double *outData)
{
   int    ss, ii, jj;
   double ddata;
   for (ss = 0; ss < length; ss++)
   {
      for (ii = 0; ii < nInputs_; ii++) 
      {
         ddata = inData[ss*nInputs_+ii];
         ddata = (ddata - lowerBs_[ii]) / (upperBs_[ii]-lowerBs_[ii]);
         if (ddata == 1.0) jj = incrs_[ii] - 1;
         else              jj = (int) (ddata * incrs_[ii]);
         indexSet_[ii] = jj;
      }
      outData[ss] = findProbability(indexSet_);
   }
   return 0;
}

// ************************************************************************
// look up cumulative density
// ------------------------------------------------------------------------
int PDFHistogram::getCDF(int length, double *inData, double *outData)
{
   printf("PDFHistogram::getCDF not available.\n");
   for (int ii = 0; ii < length; ii++) outData[ii] = 0;
   return -1;
}

// ************************************************************************
// transformation to range
// ------------------------------------------------------------------------
int PDFHistogram::invCDF(int length, double *inData, double *outData,
                       double lower, double upper)
{
   printf("PDFHistogram::invCDF not available.\n");
   for (int ii = 0; ii < length; ii++) outData[ii] = 0;
   return -1;
}

// ************************************************************************
// generate a sample
// ------------------------------------------------------------------------
int PDFHistogram::genSample(int length,double *outData,double, double)
{
   int    ii, jj, kk, ind, count;
   double ddata;

   count = 0;
   while (count < length)
   {
      ddata = PSUADE_drand();
      ind = searchHistogram(ddata);
      for (ii = 0; ii < histCnts_[ind]; ii++)
      kk = PSUADE_rand() % histCnts_[ind];
      jj = histMap_[ind][kk];
      for (ii = 0; ii < nInputs_; ii++)
         outData[count*nInputs_+ii] = samples_[jj*nInputs_+ii];
      count++;
   }
   return 0;
}

// ************************************************************************
// get mean
// ------------------------------------------------------------------------
double PDFHistogram::getMean()
{
   printf("PDFHistogram::getMean not available for this distribution.\n");
   return 0.0;
}

// ************************************************************************
// search the histogram CDF using a random number, return cell index
// ------------------------------------------------------------------------
int PDFHistogram::searchHistogram(double prob)
{
   int    ii;
   double ddata;
   if (prob < histCDF_[0]) return 0; 
   if (prob >= histCDF_[nHist_-1]) return (nHist_-1); 
   for (ii = 1; ii < nHist_; ii++)
   {
      ddata = 0.5 * (histCDF_[ii-1] + histCDF_[ii]);
      if (prob <= ddata) return (ii-1);
      else if (prob < histCDF_[ii]) return ii;
   }
   return 0;
}

// ************************************************************************
// initialize histogram
// ------------------------------------------------------------------------
void PDFHistogram::initHistogram()
{
   nHist_     = 0;
   histCells_ = NULL;
   indexSet_  = new int[nInputs_];
   histMedians_ = NULL;
   histCDF_   = NULL;
   histCnts_  = NULL;
   histMap_   = NULL;
}

// ************************************************************************
// merge histogram
// ------------------------------------------------------------------------
void PDFHistogram::mergeHistogram(int *indexSet, int samNum)
{
   int    ii, jj, index, *tmpInt, **tmpInt2;
   double *tmpDble;
   ii = 0;
   for (ii = 0; ii < nHist_; ii++)
   {
      for (jj = 0; jj < nInputs_; jj++)
         if (histCells_[ii*nInputs_+jj] != indexSet[jj]) break;
      if (jj == nInputs_) break;
   }
   index = ii;
   //*** new entry
   if (index == nHist_)
   {
      if (nHist_ % 1000 == 0)
      {
         tmpInt = histCells_;
         histCells_ = new int[(nHist_+1000)*nInputs_];
         for (ii = 0; ii < nHist_*nInputs_; ii++) 
            histCells_[ii] = tmpInt[ii];
         if (tmpInt != NULL) delete [] tmpInt;
         tmpInt = histCnts_;
         histCnts_ = new int[nHist_+1000];
         for (ii = 0; ii < nHist_; ii++) histCnts_[ii] = tmpInt[ii];
         for (ii = nHist_; ii < nHist_+1000; ii++) histCnts_[ii] = 0;
         if (tmpInt != NULL) delete [] tmpInt;
         tmpInt2 = histMap_;
         histMap_ = new int*[nHist_+1000];
         for (ii = 0; ii < nHist_; ii++) histMap_[ii] = tmpInt2[ii];
         for (ii = nHist_; ii < nHist_+1000; ii++) histMap_[ii] = NULL;
         if (tmpInt2 != NULL) delete [] tmpInt2;
         tmpDble = histMedians_;
         histMedians_ = new double[(nHist_+1000)*nInputs_];
         for (ii = 0; ii < nHist_*nInputs_; ii++) 
            histMedians_[ii] = tmpDble[ii];
         for (ii = nHist_*nInputs_; ii < (nHist_+1000)*nInputs_; ii++) 
            histMedians_[ii] = 0;
      }
      histMap_[nHist_] = new int[100]; 
      histMap_[nHist_][0] = samNum;
      histCnts_[nHist_]++;
      for (ii = 0; ii < nInputs_; ii++) 
      {
         histCells_[nHist_*nInputs_+ii] = indexSet[ii]; 
         histMedians_[nHist_*nInputs_+ii] = samples_[samNum*nInputs_+ii]; 
      }
      nHist_++;
   }
   else
   {
      if (histCnts_[index] % 100 == 0)
      {
         tmpInt = histMap_[index];
         histMap_[index] = new int[histCnts_[index]+100]; 
         for (ii = 0; ii < histCnts_[index]; ii++) 
            histMap_[index][ii] = tmpInt[ii];
      }
      histMap_[index][histCnts_[index]] = samNum;
      histCnts_[index]++;
      for (ii = 0; ii < nInputs_; ii++) 
         histMedians_[index*nInputs_+ii] += samples_[samNum*nInputs_+ii]; 
   }
}

// ************************************************************************
// finalize histogram
// ------------------------------------------------------------------------
void PDFHistogram::finalizeHistogram()
{
   int ss, ii;
   for (ss = 0; ss < nHist_; ss++) 
      for (ii = 0; ii < nInputs_; ii++) 
         histMedians_[ss*nInputs_+ii] /= (double) histCnts_[ss];
   histCDF_ = new double[nHist_]; 
   for (ss = 0; ss < nHist_; ss++)
      histCDF_[ss] = 1.0 * histCnts_[ss] / nSamples_; 
   for (ss = 1; ss < nHist_; ss++)
      histCDF_[ss] += histCDF_[ss-1];
}

// ************************************************************************
// find probability
// ------------------------------------------------------------------------
double PDFHistogram::findProbability(int *indexSet)
{
   int    ss, ii, jj;
   double ddata;
   for (ss = 0; ss < nHist_; ss++)
   {
      for (ii = 0; ii < nInputs_; ii++) 
         if (histCells_[ss*nInputs_+ii] != indexSet[ii]) break;
      if (ii == nInputs_) break;
   }
   return (1.0*histCnts_[ss]/nSamples_);
}

