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
// Functions for the sparse grid sampling class 
// AUTHOR : CHARLES TONG
// DATE   : 2010
// ************************************************************************
#include <math.h>
#include <stdio.h>
#include <sstream>
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "Vector.h"
#include "SparseGridSampling.h"

#define HUQ_nLevels_ 8
#define HUQ_nTerms_  4
static double 
HUQ_GL[HUQ_nLevels_][HUQ_nTerms_] =
{
   {0.5, 0.0, 0.0, 0.0},
   {7.8867513459481287e-1, 0.0, 0.0, 0.0},
   {5.0e-01, 8.8729833462074170e-1, 0.0, 0.0},
   {6.6999052179242813e-1, 9.3056815579702623e-1, 0.0, 0.0},
   {0.5, 7.6923465505284150e-1, 9.5308992296933193e-1},
   {6.1930959304159849e-1, 8.3060469323313235e-1, 9.6623475710157603e-1, 0},
   {0.5, 7.0292257568869854e-1, 8.7076559279969723e-1, 9.7455395617137919e-1},
   {5.9171732124782495e-1, 7.6276620495816450e-1, 8.9833323870681348e-1, 9.8014492824876809e-1}

};
static double HUQ_GLW[HUQ_nLevels_][HUQ_nTerms_] =
{
   {1.0, 0.0, 0.0, 0.0},
   {0.5, 0.0, 0.0, 0.0},
   {4.4444444444444570e-1, 2.7777777777777712e-1, 0, 0},
   {3.2607257743127516e-1, 1.7392742256872484e-1, 0, 0},
   {2.8444444444444655e-1, 2.3931433524968501e-1, 1.1846344252809174e-1, 0},
   {2.3395696728634746e-1, 1.8038078652407072e-1, 8.5662246189581834e-2, 0},
   {2.089795918367362e-1, 1.909150252525609e-1, 1.3985269574463935e-1, 6.4742483084431701e-2},
   {1.8134189168918213e-1, 1.5685332293894469e-1, 1.1119051722668793e-1, 5.0614268145185180e-2}
};
static int HUQ_GLn[HUQ_nLevels_] =
{
   1, 1, 2, 2, 3, 3, 4
};

#define PABS(x)  ((x) > 0 ? x : -(x))

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
SparseGridSampling::SparseGridSampling() : Sampling()
{
   samplingID_ = PSUADE_SAMP_SG;
   sampleWeights_ = NULL;
   pOrder_ = 1; 
}

// ************************************************************************
// copy constructor added by Bill Oliver
// ------------------------------------------------------------------------
SparseGridSampling::SparseGridSampling(const SparseGridSampling &gs) : Sampling()
{
   pOrder_ = gs.pOrder_;
   nSamples_ = gs.nSamples_;
   sampleWeights_ = new double[nSamples_];
   for(int i = 0; i < nSamples_; i++)
      sampleWeights_[i] = gs.sampleWeights_[i];
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
SparseGridSampling::~SparseGridSampling()
{
   if (sampleWeights_ != NULL) delete [] sampleWeights_;
}

//*************************************************************************
//* initialize the sampling data
//*------------------------------------------------------------------------
int SparseGridSampling::initialize(int initLevel)
{
   int    ii, jj, kk, ll, minQ, maxQ, bQ, index, size = 0;
   int    nPerms, **pcePerms, *counts, total, *midx, newLeng, nVecs, numNew;
   int    *keepFlags;
   double **newn, *neww, **ddata, val, *ranges;
   Vector *Vnodes;
   Vector  Vweights;
   FILE   *fp;

   if (nInputs_ == 0 || lowerBounds_ == NULL || upperBounds_ == NULL)
   {
      printf("SparseGridSampling::initialize ERROR - input not set up.\n");
      exit(1);
   }

   deleteSampleData();

   if (initLevel != 0) return 0;

   pOrder_ = 0;
   while (pOrder_ < 2 || pOrder_ > 7)
   {
      printf("SparseGridSampling: enter polynomial order (2 - 7): ");
      scanf("%d", &pOrder_);
   }
   if (pOrder_ > 7)
   {
      printf("SparseGridSampling ERROR: does not support pOrder > 7.\n");
      exit(1);
   }

   minQ = pOrder_ - nInputs_;
   if (minQ < 0) minQ = 0;
   maxQ = pOrder_ - 1;
   midx = new int[nInputs_];
   Vnodes = new Vector[nInputs_];
   nVecs = 0;
   for (ii = minQ; ii <= maxQ; ii++)
   {
      val   = (pow(-1.0, 1.0*(maxQ-ii))) * 
               nChooseK(nInputs_-1,nInputs_+ii-pOrder_);
      bQ = (int) val;
      nPerms = GenSequence(nInputs_, nInputs_+ii, &pcePerms);
      counts = new int[nPerms];
      total = 0;
      for (jj = 0; jj < nPerms; jj++)
      {
         counts[jj] = 1;
         for (kk = 0; kk < nInputs_; kk++)
         {
            index = pcePerms[jj][kk];
            counts[jj] *= HUQ_GLn[index];
         } 
         total += counts[jj];
      }
      for (kk = 0; kk < nInputs_; kk++) Vnodes[kk].add(total, NULL);
      Vweights.add(total,NULL);
      for (jj = 0; jj < nPerms; jj++)
      {
         for (kk = 0; kk < nInputs_; kk++)
         {
            index = pcePerms[jj][kk];
            midx[kk] = index;
         }

         size = KronProd(nInputs_, midx, &newn, &neww);
         // Bill Oliver added check on value of size
         if(size <= 0){
	   printf("size variable is <= 0 in file %s line %d\n", __FILE__, __LINE__);
           exit(1);

         }
         for (ll = nVecs; ll < nVecs+counts[jj]; ll++)
         {
            for (kk = 0; kk < nInputs_; kk++)
               Vnodes[kk][ll] = newn[ll-nVecs][kk];
            Vweights[ll] = bQ * neww[ll-nVecs];
         } 
         nVecs += counts[jj];
	 for(int i = 0; i < size; i++)
	   delete newn[i];
	 delete newn;
	 delete neww;
         
      }
      delete [] counts;
      for (jj = 0; jj < nPerms; jj++) delete [] pcePerms[jj];
      delete [] pcePerms;
      /* need to prune repeated ones */
      ddata = new double*[nVecs];
      for (jj = 0; jj < nVecs; jj++) 
      {
         ddata[jj] = new double[nInputs_+1];
         for (kk = 0; kk < nInputs_; kk++) 
            ddata[jj][kk] = Vnodes[kk][jj];
         ddata[jj][nInputs_] = Vweights[jj];
      }
      newLeng = sortAndDelete(nVecs, nInputs_+1, ddata);
      for (jj = 0; jj < newLeng; jj++) 
      {
         for (kk = 0; kk < nInputs_; kk++) 
            Vnodes[kk][jj] = ddata[jj][kk];
         Vweights[jj] = ddata[jj][nInputs_];
      }
      for (jj = 0; jj < nVecs; jj++) delete [] ddata[jj];
      delete [] ddata;
      nVecs = newLeng;
   }

   val = HUQ_GL[0][0];
   for (ii = 0; ii < nInputs_; ii++)
   {
      keepFlags = new int[nVecs];
      for (jj = 0; jj < nVecs; jj++) keepFlags[jj] = -1;
      numNew = 0;
      for (jj = 0; jj < nVecs; jj++) 
      {
         if (Vnodes[ii][jj] != val)
         {
            keepFlags[numNew] = jj;
            numNew++;
         }
      }
      if (numNew > 0)
      {
         for (jj = 0; jj < nInputs_; jj++) Vnodes[jj].add(numNew, NULL);
         Vweights.add(numNew,NULL);
         for (jj = 0; jj < nInputs_; jj++)
         {
            numNew = 0;
            for (kk = 0; kk < nVecs; kk++) 
            {
               if (keepFlags[kk] != -1)
               {
                  Vnodes[jj][nVecs+numNew] = Vnodes[jj][keepFlags[kk]];
                  Vweights[nVecs+numNew] = Vweights[keepFlags[kk]];
                  numNew++;
               }
            }
         }
         for (kk = nVecs; kk < nVecs+numNew; kk++) 
            Vnodes[ii][kk] = 2 * val - Vnodes[ii][kk];
         nVecs += numNew;
      }
      delete [] keepFlags;
   }
   /* need to prune repeated ones */
   ddata = new double*[nVecs];
   for (jj = 0; jj < nVecs; jj++) 
   {
      ddata[jj] = new double[nInputs_+1];
      for (kk = 0; kk < nInputs_; kk++) 
         ddata[jj][kk] = Vnodes[kk][jj];
      ddata[jj][nInputs_] = Vweights[jj];
   }
   newLeng = sortAndDelete(nVecs, nInputs_+1, ddata);
   for (jj = 0; jj < newLeng; jj++) 
   {
      for (kk = 0; kk < nInputs_; kk++) 
         Vnodes[kk][jj] = ddata[jj][kk];
      Vweights[jj] = ddata[jj][nInputs_];
   }
   for (jj = 0; jj < nVecs; jj++) delete [] ddata[jj];
   delete [] ddata;
   nVecs = newLeng;
   val = 0.0;
   for (jj = 0; jj < nVecs; jj++) val += Vweights[jj];
   for (jj = 0; jj < nVecs; jj++) Vweights[jj] /= val;

   nSamples_ = nVecs;
   allocSampleData();
   ranges = new double[nInputs_];
   for (ii = 0; ii < nInputs_; ii++) 
      ranges[ii] = upperBounds_[ii] - lowerBounds_[ii];
   sampleWeights_ = new double[nSamples_];
   for (kk = 0; kk < nSamples_; kk++)
   {
      for (ii = 0; ii < nInputs_; ii++) 
         sampleMatrix_[kk][ii] = 
                   Vnodes[ii][kk] * ranges[ii] + lowerBounds_[ii];
   }
   delete [] ranges;
   
   fp = fopen("ps_sparse_grid_info", "w");
   fprintf(fp, "%d %d %d\n", nSamples_, nInputs_, pOrder_);
   for (kk = 0; kk < nSamples_; kk++)
   {
      for (ii = 0; ii < nInputs_; ii++) 
         fprintf(fp, "%24.16e ", Vnodes[ii][kk]);
      sampleWeights_[kk] = Vweights[kk]; 
      fprintf(fp, "%24.16e ", sampleWeights_[kk]);
      fprintf(fp, "\n");
   }
   for (ii = 0; ii < nInputs_; ii++) 
      fprintf(fp, "%24.16e %24.16e\n", lowerBounds_[ii], upperBounds_[ii]);
   fclose(fp);
   printf("Sparse grid information has been stored to ps_sparse_grid_info.\n");
   printf("You need this file to build sparse grid response surfaces.\n");

   delete [] midx;
   delete [] Vnodes;
   return 0;
}

// ************************************************************************
// refine the sample space
// ------------------------------------------------------------------------
int SparseGridSampling::refine(int, int, double , int, double *)
{
   printf("SparseGridSampling::refine ERROR - not available.\n");
   exit(1);
   return 0;
}

// ************************************************************************
// Kronecker product 
// ------------------------------------------------------------------------
int SparseGridSampling::KronProd(int nn, int *midx, double ***newn2, 
                                 double **neww2)
{
   int    total, ii, index, *counts, cnt;
   double **newn, *neww;

   /* total number of new rows */
   total = 1;
   for (ii = 0; ii < nn; ii++)
   {
      index = midx[ii];
      total *= HUQ_GLn[index];
   }

   /* create storage */
   newn = new double*[total];
   for (ii = 0; ii < total; ii++)
      newn[ii] = new double[nn];
   neww = new double[total];
   counts = new int[nn];
   for (ii = 0; ii < nn; ii++) counts[ii] = 0;
   for (ii = 0; ii < total; ii++) neww[ii] = 1.0;

   /* create data */
   cnt = 0;
   while (cnt < total)
   {
      for (ii = 0; ii < nn; ii++) 
      {
         newn[cnt][ii] = HUQ_GL[midx[ii]][counts[ii]];
         neww[cnt] *= HUQ_GLW[midx[ii]][counts[ii]];
      }
      counts[nn-1]++;
      ii = nn - 1; 
      while (counts[ii] >= midx[ii] && ii > 0)
      {
         counts[ii-1]++;
         counts[ii] = 0;
         ii--;
      }
      cnt++;
   }
   delete [] counts;
   (*newn2) = newn;    
   (*neww2) = neww;
   // This function is used in a loop so we need to return total    
   return total;
}

// ************************************************************************
// generate sequence
// ------------------------------------------------------------------------
int SparseGridSampling::GenSequence(int nn, int rsum, int ***perms)
{
   int *flags, **outPerms, nPerms, cnt, ii, jj, idata, idata2;
   flags = new int[nn];
   for (ii = 0; ii < nn; ii++) flags[ii] = 0;
   idata = rsum - nn;
   flags[0] = idata;
   if      ((rsum-nn) == 0) nPerms = 1;
   else if ((rsum-nn) == 1) nPerms = nn;
   else                     
      nPerms = computeNumPCEPermutations(nn, idata) - 
               computeNumPCEPermutations(nn, idata-1);
   // Bill Oliver added defensive programming by checking the value of nPerms
   if(nPerms <= 0)
   {
      printf("Problem in SparseGridSampling::GenSequence with calculation\n");
      printf("of nPerms (%d).\n", nPerms);
      exit(1);
   }
   outPerms = new int*[nPerms];
   for (ii = 0; ii < nPerms; ii++) outPerms[ii] = new int[nn];
   for (ii = 0; ii < nn; ii++) outPerms[0][ii] = flags[ii];
   idata2 = 0;
   cnt = 1;
   while (flags[nn-1] < idata)
   {
      if (idata2 == nn-1)
      {
         for (jj = idata2-1; jj >= 0; jj--)
         {
            idata2 = jj;
            if (flags[jj] != 0) break;
         }
      }
      flags[idata2]--;
      idata2++;
      flags[idata2] = idata;
      for (jj = 0; jj < idata2; jj++) flags[idata2] -= flags[jj];
      if (idata2 < nn)
         for (jj = idata2+1; jj < nn; jj++) flags[jj] = 0;
      for (jj = 0; jj < nn; jj++) outPerms[cnt][jj] = flags[jj];
      cnt++;
   }
   (*perms) = outPerms;
   delete [] flags;
   return nPerms;
}

// ************************************************************************
// n choose k
// ------------------------------------------------------------------------
int SparseGridSampling::nChooseK(int n, int k)
{
   int ii, idata=n;

   idata = 1;
   for (ii = n; ii > k; ii--) idata *= ii;
   for (ii = 2; ii <= n-k; ii++) idata /= ii;
   return idata;
}

// ************************************************************************
// set input settings
// ------------------------------------------------------------------------
int SparseGridSampling::setParam(string sparam)
{
   int           pos;
   istringstream buffer;
   string        substr;

   pos = sparam.find("pOrder");
   if (pos >= 0)
   {
      substr = sparam.substr(8);
      buffer.str(substr);
      buffer >> pOrder_;
      if (pOrder_ < 1) pOrder_ = 1;
   }
   return 0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
SparseGridSampling& SparseGridSampling::operator=(const SparseGridSampling & gs)
{
   // Bill Oliver modified the operator= to work
   if(this == &gs) return *this;
   delete [] sampleWeights_;
   pOrder_ = gs.pOrder_;
   nSamples_ = gs.nSamples_;
   sampleWeights_ = new double[nSamples_];
   for(int i = 0; i < nSamples_; i++)
      sampleWeights_[i] = gs.sampleWeights_[i];
   return (*this);
}

