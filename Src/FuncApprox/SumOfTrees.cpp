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
// Functions for the class SumOfTrees
// AUTHOR : CHARLES TONG
// DATE   : 2010
// ************************************************************************
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "SumOfTrees.h"
#include "Util/sysdef.h"
#include "Util/PsuadeUtil.h"
#include "Main/Psuade.h"

#define PABS(x) ((x) > 0 ? (x) : (-(x)))

// ************************************************************************
// Constructor for object class SumOfTrees
// ------------------------------------------------------------------------
SumOfTrees::SumOfTrees(int nInputs,int nSamples) : 
                 FuncApprox(nInputs,nSamples)
{
   char pString[501];

   faID_ = PSUADE_RS_SOTS;

   mode_ = 0;

   numTrees_ = 100;
   treePtrs_ = NULL;
   minPtsPerNode_ = 10;
   shrinkFactor_ = 0.005;
   tolerance_ = 1.0e-6;

   if (psRSExpertMode_ == 1)
   {
      sprintf(pString,"Enter the desired number of trees (>10): ");
      numTrees_ = getInt(10, 10000, pString);
      sprintf(pString,"Enter minimum points per node (>1, default=10): ");
      minPtsPerNode_ = getInt(2, nSamples/5, pString);
   }

   treePtrs_ = new TreeNode*[numTrees_];
   for (int ii = 0; ii < numTrees_; ii++)
      treePtrs_[ii] = new TreeNode();
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
SumOfTrees::~SumOfTrees()
{
   if (treePtrs_ != NULL)
   {
      for (int ii = 0; ii < numTrees_; ii++)
         if (treePtrs_[ii] != NULL) delete treePtrs_[ii];
   }
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int SumOfTrees::genNDGridData(double *X, double *Y, int *N, double **X2, 
                              double **Y2)
{
   int    ii, ss, totPts;
   double *XX, *YY, *XL, *HX;

   initTrees();

   buildTrees(X, Y);

   if ((*N) == -999) return 0;
  
   totPts = nPtsPerDim_;
   for (ii = 1; ii < nInputs_; ii++) totPts = totPts * nPtsPerDim_;
   HX = new double[nInputs_];
   for (ii = 0; ii < nInputs_; ii++) 
      HX[ii] = (upperBounds_[ii]-lowerBounds_[ii]) / (double) (nPtsPerDim_-1); 
 
   (*X2) = new double[totPts*nInputs_];
   XX = (*X2);
   (*Y2) = new double[totPts];
   YY = (*Y2);
   XL = new double[nInputs_];
   (*N) = totPts;
 
   for (ii = 0; ii < nInputs_; ii++) XL[ii] = lowerBounds_[ii];
 
   for (ss = 0; ss < totPts; ss++)
   {
      for (ii = 0; ii < nInputs_; ii++ ) XX[ss*nInputs_+ii] = XL[ii];
      for (ii = 0; ii < nInputs_; ii++ ) 
      {
         XL[ii] += HX[ii];
         if (XL[ii] < upperBounds_[ii] || 
              PABS(XL[ii] - upperBounds_[ii]) < 1.0E-7) break;
         else XL[ii] = lowerBounds_[ii];
      }
   }
 
   evaluatePoint(totPts, XX, YY);

   delete [] XL;
   delete [] HX;
   return 0;
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int SumOfTrees::gen1DGridData(double *X, double *Y, int ind1, 
                              double *settings, int *N, double **X2, 
                              double **Y2)
{
   int    ii, ss, totPts;
   double *XT, *XX, *YY, HX;

   initTrees();

   buildTrees(X, Y);

   totPts = nPtsPerDim_;
   HX = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 

   (*X2) = new double[totPts];
   XX = (*X2);
   (*Y2) = new double[totPts];
   YY = (*Y2);
   (*N) = totPts;

   XT = new double[totPts*nInputs_];
   for (ss = 0; ss < totPts; ss++) 
      for (ii = 0; ii < nInputs_; ii++) XT[ss*nInputs_+ii] = settings[ii]; 
    
   for (ss = 0; ss < totPts; ss++) 
   {
      XT[ss*nInputs_+ind1]  = HX * ss + lowerBounds_[ind1];
      XX[ss] = HX * ss + lowerBounds_[ind1];
      YY[ss] = 0.0;
   }

   evaluatePoint(totPts, XT, YY);

   delete [] XT;
   return 0;
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int SumOfTrees::gen2DGridData(double *X, double *Y, int ind1, int ind2, 
                              double *settings, int *N, double **X2, 
                              double **Y2)
{
   int    ii, ss, jj, index, totPts;
   double *XT, *XX, *YY, *HX;
 
   initTrees();

   buildTrees(X, Y);

   totPts = nPtsPerDim_ * nPtsPerDim_;
   HX = new double[2];
   HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
   HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 

   (*X2) = new double[2*totPts];
   XX = (*X2);
   (*Y2) = new double[totPts];
   YY = (*Y2);
   (*N) = totPts;

   XT = new double[totPts*nInputs_];
   for (ss = 0; ss < totPts; ss++) 
      for (ii = 0; ii < nInputs_; ii++) XT[ss*nInputs_+ii] = settings[ii]; 
    
   for (ii = 0; ii < nPtsPerDim_; ii++) 
   {
      for (jj = 0; jj < nPtsPerDim_; jj++)
      {
         index = ii * nPtsPerDim_ + jj;
         XT[index*nInputs_+ind1] = HX[0] * ii + lowerBounds_[ind1];
         XT[index*nInputs_+ind2] = HX[1] * jj + lowerBounds_[ind2];
         XX[index*2]   = HX[0] * ii + lowerBounds_[ind1];
         XX[index*2+1] = HX[1] * jj + lowerBounds_[ind2];
      }
   }

   evaluatePoint(totPts, XT, YY);

   delete [] XT;
   delete [] HX;
   return 0;
}

// ************************************************************************
// Generate 3D results for display
// ------------------------------------------------------------------------
int SumOfTrees::gen3DGridData(double *X, double *Y, int ind1, int ind2, 
                              int ind3, double *settings, int *N, 
                              double **X2, double **Y2)
{
   int    ii, ss, jj, ll, index, totPts;
   double *XT, *XX, *YY, *HX;

   initTrees();

   buildTrees(X, Y);

   totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
   HX = new double[3];
   HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
   HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 
   HX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_ - 1); 

   (*X2) = new double[3*totPts];
   XX = (*X2);
   (*Y2) = new double[totPts];
   YY = (*Y2);
   (*N) = totPts;

   XT = new double[totPts*nInputs_];
   for (ss = 0; ss < totPts; ss++)
      for (ii = 0; ii < nInputs_; ii++) XT[ss*nInputs_+ii] = settings[ii];

   for (ii = 0; ii < nPtsPerDim_; ii++) 
   {
      for (jj = 0; jj < nPtsPerDim_; jj++)
      {
         for (ll = 0; ll < nPtsPerDim_; ll++)
         {
            index = ii * nPtsPerDim_ * nPtsPerDim_ + jj * nPtsPerDim_ + ll;
            XT[index*nInputs_+ind1]  = HX[0] * ii + lowerBounds_[ind1];
            XT[index*nInputs_+ind2]  = HX[1] * jj + lowerBounds_[ind2];
            XT[index*nInputs_+ind3]  = HX[2] * ll + lowerBounds_[ind3];
            XX[index*3]   = HX[0] * ii + lowerBounds_[ind1];
            XX[index*3+1] = HX[1] * jj + lowerBounds_[ind2];
            XX[index*3+2] = HX[2] * ll + lowerBounds_[ind3];
         }
      }
   }

   evaluatePoint(totPts, XT, YY);

   delete [] XT;
   delete [] HX;
   return 0;
}

// ************************************************************************
// Generate 4D results for display
// ------------------------------------------------------------------------
int SumOfTrees::gen4DGridData(double *X, double *Y, int ind1, int ind2, 
                              int ind3, int ind4, double *settings, 
                              int *N, double **X2, double **Y2)
{
   int    ii, ss, jj, ll, mm, index, totPts;
   double *XT, *XX, *YY, *HX;

   initTrees();

   buildTrees(X, Y);

   totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
   HX = new double[4];
   HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
   HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 
   HX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_ - 1); 
   HX[3] = (upperBounds_[ind4] - lowerBounds_[ind4]) / (nPtsPerDim_ - 1); 

   (*X2) = new double[4*totPts];
   XX = (*X2);
   (*Y2) = new double[totPts];
   YY = (*Y2);
   (*N) = totPts;

   XT = new double[totPts*nInputs_];
   for (ss = 0; ss < totPts; ss++) 
      for (ii = 0; ii < nInputs_; ii++) XT[ss*nInputs_+ii] = settings[ii]; 
    
   for (ii = 0; ii < nPtsPerDim_; ii++) 
   {
      for (jj = 0; jj < nPtsPerDim_; jj++)
      {
         for (ll = 0; ll < nPtsPerDim_; ll++)
         {
            for (mm = 0; mm < nPtsPerDim_; mm++)
            {
               index = ii*nPtsPerDim_*nPtsPerDim_ * nPtsPerDim_ +
                       jj*nPtsPerDim_*nPtsPerDim_ + ll*nPtsPerDim_ + mm;
               XT[index*nInputs_+ind1]  = HX[0] * ii + lowerBounds_[ind1];
               XT[index*nInputs_+ind2]  = HX[1] * jj + lowerBounds_[ind2];
               XT[index*nInputs_+ind3]  = HX[2] * ll + lowerBounds_[ind3];
               XT[index*nInputs_+ind4]  = HX[3] * mm + lowerBounds_[ind4];
               XX[index*4]   = HX[0] * ii + lowerBounds_[ind1];
               XX[index*4+1] = HX[1] * jj + lowerBounds_[ind2];
               XX[index*4+2] = HX[2] * ll + lowerBounds_[ind3];
               XX[index*4+3] = HX[3] * mm + lowerBounds_[ind4];
            }
         }
      }
   }

   evaluatePoint(totPts, XT, YY);

   delete [] XT;
   delete [] HX;
   return 0;
}

// ************************************************************************
// Evaluate a given point 
// ------------------------------------------------------------------------
double SumOfTrees::evaluatePoint(double *X)
{
   int    ii;
   double Y=0.0;

   for (ii = 0; ii < numTrees_; ii++) Y += evaluateTree(treePtrs_[ii], X);
   Y /= (double) numTrees_;
   return Y;
}

// ************************************************************************
// Evaluate a number of points
// ------------------------------------------------------------------------
double SumOfTrees::evaluatePoint(int npts, double *X, double *Y)
{
   int ii, ss;

   for (ss = 0; ss < npts; ss++)
   {
      Y[ss] = 0.0;
      for (ii = 0; ii < numTrees_; ii++)
         Y[ss] += evaluateTree(treePtrs_[ii], &X[ss*nInputs_]);
      if (mode_ == 0) Y[ss] /= (double) numTrees_;
   }
   return 0.0;
}

// ************************************************************************
// Evaluate a given point with standard deviation
// ------------------------------------------------------------------------
double SumOfTrees::evaluatePointFuzzy(double *X, double &std)
{
   double Y=0.0;
   Y = evaluatePoint(X);
   std = 0.0;
   return Y;
}

// ************************************************************************
// Evaluate a number of points with standard deviations
// ------------------------------------------------------------------------
double SumOfTrees::evaluatePointFuzzy(int npts, double *X, double *Y,
                                       double *Ystd)
{
   evaluatePoint(npts, X, Y);
   for (int ss = 0; ss < npts; ss++) Ystd[ss] = 0.0;
   return 0.0;
}

// ************************************************************************
// initialize the trees
// ------------------------------------------------------------------------
int SumOfTrees::initTrees()
{
   int    ii;

   if (treePtrs_ != NULL)
   {
      for (ii = 0; ii < numTrees_; ii++)
         if (treePtrs_[ii] != NULL) delete treePtrs_[ii];
   }
   for (ii = 0; ii < numTrees_; ii++)
      treePtrs_[ii] = new TreeNode();
   return 0;
}   

// ************************************************************************
// build a number of trees 
// ------------------------------------------------------------------------
int SumOfTrees::buildTrees(double *X, double *Y)
{
   int      jj, kk, mm, ss;
   double   *XT, *YT, *YY, mult, checksum, *inputScores;

   //** allocate temporary storage
   XT = new double[nSamples_*nInputs_];
   YT = new double[nSamples_];
   YY = new double[nSamples_];
   inputScores = new double[nInputs_];

   //** initialize local variables (YY keeps intermediate Y)
   for (jj = 0; jj < nInputs_; jj++) inputScores[jj] = 0.0; 
   for (jj = 0; jj < nSamples_; jj++) YY[jj] = Y[jj];

   mult = shrinkFactor_;
   for (mm = 0; mm < numTrees_; mm++)
   {
      if (mode_ == 0)
      {
         for (ss = 0; ss < nSamples_; ss++)
         {
            kk = PSUADE_rand() % nSamples_;
            for (jj = 0; jj < nInputs_; jj++)
               XT[ss*nInputs_+jj] = X[kk*nInputs_+jj];
            YT[ss] = Y[kk];
         }
      }
      else
      {
         for (jj = 0; jj < nSamples_*nInputs_; jj++) XT[jj] = X[jj];
         for (jj = 0; jj < nSamples_; jj++) YT[jj] = YY[jj];
      }
      buildOneTree(treePtrs_[mm], nSamples_, XT, YT, 0);
      if (mode_ == 1)
      {
         for (jj = 0; jj < nSamples_; jj++) 
            YY[jj] -= (mult * evaluateTree(treePtrs_[mm], &X[jj*nInputs_]));
         checksum = 0.0;
         for (jj = 0; jj < nSamples_; jj++) checksum += YY[jj] * YY[jj];
         printf("SumOfTree boosting at iteration %d = %e\n", mm+1, 
                checksum/nSamples_);
         mult = shrinkFactor_ * mult / (shrinkFactor_ - mult * mult);
         if (mm == numTrees_-1) mult = 1.0;
         printf("mult %d = %e\n", mm, mult);
      }
   }

   delete [] XT;
   delete [] YT;
   delete [] YY;
   delete [] inputScores;
   return 0;
}

// ************************************************************************
// build a tree 
// ------------------------------------------------------------------------
int SumOfTrees::buildOneTree(TreeNode *tnode, int leng, double *XT, 
                             double *YT, int level)
{
   int    ii, jj, kk, minInd, repeatCnt, inpImpurity, locImpurity;
   double mean1, mean2, sosSum, sosMin, Ymax, *X1, *Y1;
   double lmean1, lmean2, sos1, sos2, sos, mean, maxImpurity, ddata;

   if (leng < 2*minPtsPerNode_)
   {
      tnode->leftNode_ = NULL;
      tnode->rightNode_ = NULL;
      tnode->whichInput_ = -1;
      mean = 0.0;
      for (jj = 0; jj < leng; jj++) mean += YT[jj];
      if (leng > 0) mean /= (double) leng;
      else          mean = 0.0;
      tnode->nodeValue_ = mean;
      ddata = 0.0;
      for (jj = 0; jj < leng; jj++) ddata += pow(YT[jj]-mean,2.0);
      tnode->nodeStdev_ = sqrt(ddata/leng);
      return -1;
   }

   Ymax = PABS(YT[0]);
   for (jj = 1; jj < leng; jj++) 
      if (PABS(YT[jj]) > Ymax) Ymax = PABS(YT[jj]);
   if (Ymax == 0.0) Ymax = 1.0;
      
   X1 = new double[leng*nInputs_];
   Y1 = new double[leng];
   for (jj = 0; jj < leng; jj++) Y1[jj] = YT[jj] / Ymax;
   mean = sos = 0.0;
   for (jj = 0; jj < leng; jj++) mean += Y1[jj];
   mean /= leng;
   for (jj = 0; jj < leng; jj++) sos += (Y1[jj]-mean)*(Y1[jj]-mean);
   sos /= leng;

   maxImpurity = - PSUADE_UNDEFINED;
   inpImpurity = -1;
   for (ii = 0; ii < nInputs_; ii++)
   {
      for (jj = 0; jj < leng; jj++)
      {
         X1[jj] = XT[jj*nInputs_+ii];
         Y1[jj] = YT[jj] / Ymax;
      }
      sortDbleList2(leng, X1, Y1);
      mean1 = sos1 = 0.0;
      for (jj = 0; jj < 2; jj++) mean1 += Y1[jj];
      mean1 /= 2.0;
      for (jj = 0; jj < 2; jj++) sos1 += (Y1[jj] - mean1) * (Y1[jj] - mean1);
      mean2 = sos2 = 0.0;
      for (jj = 2; jj < leng; jj++) mean2 += Y1[jj];
      mean2 /= (leng - 2.0);
      for (jj = 2; jj < leng; jj++)
         sos2 += (Y1[jj] - mean2) * (Y1[jj] - mean2);
      sosSum = (sos1 + sos2) / leng;
      sosMin = sosSum;
      minInd = 1;
      repeatCnt = 0;
      for (jj = 2; jj < leng - 2; jj++)
      {
         lmean1 = mean1;
         lmean2 = mean2;
         mean1 = (mean1 * jj + Y1[jj]) / (jj + 1.0);
         mean2 = (mean2 * (leng - jj) - Y1[jj]) / (leng-jj-1.0);
         sos1 = sos1 + lmean1 * lmean1 * jj + Y1[jj] * Y1[jj] - 
                (jj + 1.0) * mean1 * mean1;
         sos2 = sos2 + lmean2 * lmean2 * (leng - jj) - Y1[jj] * Y1[jj] - 
                (leng - jj - 1.0) * mean2 * mean2;
#if 0
#endif
         sosSum = (sos1 + sos2) / leng;
         if (sosSum == sosMin) repeatCnt++;
         else                  repeatCnt = 0; 
         if (sosSum < sosMin) 
         {
            sosMin = sosSum;
            minInd = jj;
         }
      }
      //printf("SumOfTree: input = %3d, max decrease in impurity at %d (%d)\n",
      //       ii+1, minInd+repeatCnt/2,leng);
      //printf("SumOfTree: input = %3d, actual decrease in impurity = %e\n",
      //       ii+1, sos - sosMin);
      ddata = sos - sosMin;
      if (ddata > maxImpurity)
      {
         maxImpurity = ddata;
         inpImpurity = ii;
         locImpurity = minInd + repeatCnt / 2;
      }
   }
   if (outputLevel_ >= 2)
      printf("SumOfTree: level = %d, leng = %d, input selected = %d, impurity = %e\n",
             level, leng, inpImpurity+1, maxImpurity);

   if (maxImpurity / sos < tolerance_) inpImpurity = -1;

   if (inpImpurity >= 0)
   {
      for (jj = 0; jj < leng; jj++)
      {
         X1[jj] = XT[jj*nInputs_+inpImpurity];
         Y1[jj] = 1.0 * jj;
      }
      sortDbleList2(leng, X1, Y1);
      for (jj = 0; jj < leng; jj++)
      {
         kk = (int) (Y1[jj] + 0.0001);
         for (ii = 0; ii < nInputs_; ii++)
            X1[jj*nInputs_+ii] = XT[kk*nInputs_+ii];
         Y1[jj] = YT[kk];
      }
      for (jj = 0; jj < leng*nInputs_; jj++) XT[jj] = X1[jj];
      for (jj = 0; jj < leng; jj++) YT[jj] = Y1[jj];

      tnode->cutPoint_ = XT[locImpurity*nInputs_+inpImpurity]; 
      tnode->whichInput_ = inpImpurity; 
      if (outputLevel_ >= 3)
         printf("SumTree Level = %d, Cutting input %d at %e\n", level, 
                inpImpurity+1, tnode->cutPoint_);
      tnode->leftNode_ = new TreeNode();
      buildOneTree(tnode->leftNode_, locImpurity+1, XT, YT, level+1);
      tnode->rightNode_ = new TreeNode();
      buildOneTree(tnode->rightNode_, leng-locImpurity-1, 
                   &XT[(locImpurity+1)*nInputs_],
                   &YT[locImpurity+1], level+1);
   }
   else
   {
      tnode->leftNode_ = NULL;
      tnode->rightNode_ = NULL;
      tnode->whichInput_ = -1;
      mean = 0.0;
      for (jj = 0; jj < leng; jj++) mean += YT[jj];
      if (leng > 0) mean /= (double) leng;
      else          mean = 0.0;
      tnode->nodeValue_ = mean;
      ddata = 0.0;
      for (jj = 0; jj < leng; jj++) ddata += pow(YT[jj]-mean,2.0);
      tnode->nodeStdev_ = sqrt(ddata/leng);
   }
   delete [] X1;
   delete [] Y1;
   return 0;
}

// ************************************************************************
// get the information about splitting 
// ------------------------------------------------------------------------
double SumOfTrees::tabulateSplits(TreeNode *tptr, double *scores, int level)
{
   double accum=0;
   if (tptr->whichInput_ < 0) return tptr->nodeStdev_;
   scores[tptr->whichInput_] += 1.0 / pow(2.0, 1.0*level);
   if (tptr->leftNode_  != NULL)
      accum += tabulateSplits(tptr->leftNode_,scores,level+1);
   if (tptr->rightNode_ != NULL) 
      accum += tabulateSplits(tptr->rightNode_,scores,level+1);
   return accum;
}

// ************************************************************************
// Evaluate a point for a tree
// ------------------------------------------------------------------------
double SumOfTrees::evaluateTree(TreeNode *tptr, double *X)
{
   int    ind;
   double Y=0.0;
   if (tptr->whichInput_ < 0) 
   {
      Y = tptr->nodeValue_;
      //for (int ii = 0; ii < nInputs_; ii++)
      //   printf("X %2d = %e\n", ii+1, X[ii]);
      //printf("returns Y = %e\n", Y);
   }
   else                        
   {
      ind = tptr->whichInput_;
      if (X[ind] <= tptr->cutPoint_) 
      {
         //printf("Go to X %d (%e) < %e\n", ind, X[ind], tptr->cutPoint_);
         Y = evaluateTree(tptr->leftNode_, X); 
         //printf("returns L Y = %e\n", Y);
      }
      else
      {
         //printf("Go to X %d (%e) > %e\n", ind, X[ind], tptr->cutPoint_);
         Y = evaluateTree(tptr->rightNode_, X); 
         //printf("returns R Y = %e\n", Y);
      }
   }
   return Y;
}

// ************************************************************************
// set parameters
// ------------------------------------------------------------------------
double SumOfTrees::setParams(int targc, char **targv)
{
   int    mm, ii, *iArray;
   double *inputScores, mmax;
   FILE   *fp;

   if (targc > 0 && !strcmp(targv[0], "mode0"))
   {
      mode_ = 0;
   }
   else if (targc > 0 && !strcmp(targv[0], "mode1"))
   {
      mode_ = 1;
   }
   else if (targc > 0 && !strcmp(targv[0], "rank"))
   {
      double val=0.0;
      inputScores = new double[nInputs_];
      for (ii = 0; ii < nInputs_; ii++) inputScores[ii] = 0.0;
      for (mm = 0; mm < numTrees_; mm++)
      {
         val += tabulateSplits(treePtrs_[mm], inputScores, 0);
      }
      val /= (double) numTrees_;
      mmax = 0.0;
      for (ii = 0; ii < nInputs_; ii++)
         if (inputScores[ii] > mmax) mmax = inputScores[ii];
      if (mmax > 0)
      {
         for (ii = 0; ii < nInputs_; ii++)
            inputScores[ii] = inputScores[ii] / mmax * 100;
      }

      fp = fopen("matlabsot.m", "w");
      fwritePlotCLF(fp);
      fprintf(fp, "A = [\n");
      for (ii = 0; ii < nInputs_; ii++)
         fprintf(fp, "%e\n", 0.01 * inputScores[ii]);
      fprintf(fp, "];\n");
      fprintf(fp, "bar(A, 0.8);\n");
      fwritePlotAxes(fp);
      fwritePlotTitle(fp, "Sum-of-trees Rankings");
      fwritePlotXLabel(fp, "Input parameters");
      fwritePlotYLabel(fp, "Sum-of-trees Metric (normalized)");
      fclose(fp);
      printf("Sum-of-trees ranking is now in matlabsot.m.\n");

      iArray = new int[nInputs_];
      for (ii = 0; ii < nInputs_; ii++) iArray[ii] = ii;
      sortDbleList2a(nInputs_, inputScores, iArray);
      if (targc == 1)
      {
         printAsterisks(0);
         printf("* SumOfTrees screening rankings ");
         if (mode_ == 0) printf("(with bootstrapping)\n");
         else            printf("(with boosting)\n");
         printf("* Minimum points per node = %d\n", minPtsPerNode_);
         printAsterisks(0);
         for (ii = nInputs_-1; ii >= 0; ii--)
            printf("*  Rank %3d : Input = %3d (score = %4.1f)\n",
                   nInputs_-ii, iArray[ii]+1, inputScores[ii]);
         printAsterisks(0);
      }
      delete [] iArray;
      delete [] inputScores;
      return val;
   }
   return 0.0;
}


