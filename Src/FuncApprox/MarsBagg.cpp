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
// Functions for the class Mars
// AUTHOR : CHARLES TONG
// DATE   : 2008
// ************************************************************************
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "MarsBagg.h"
#include "Util/sysdef.h"
#include "Util/PsuadeUtil.h"
#include "Main/Psuade.h"

#define PABS(x) ((x) > 0 ? (x) : (-(x)))

// ************************************************************************
// Constructor for object class MarsBagg
// ------------------------------------------------------------------------
MarsBagg::MarsBagg(int nInputs,int nSamples) : FuncApprox(nInputs,nSamples)
{
#ifdef HAVE_MARS
   int  ii, itmp;
   char pString[500], *targv[3];

   faID_ = PSUADE_RS_MARSB;

   numMars_ = 51;

   mode_ = 0;

   printAsterisks(0);
   maxBasis_ = 51;
   if (maxBasis_ > nSamples-10) maxBasis_ = nSamples - 10;
   varPerBasis_ = nInputs;
   if (outputLevel_ > 1)
   {
      printf("MARS with bagging: default mode = mean (options: mean, median).\n");
      printf("                   number of instantiations   = 51\n");
      printf("                   no. of basis functions     = %d\n",maxBasis_);
      printf("                   no. of variables per basis = %d\n",varPerBasis_);
   }
   if (psRSExpertMode_ == 1)
   {
      sprintf(pString, "MARS with bagging: mean (0) or median (1) mode ? ");
      mode_ = getInt(0, 1, pString);
      sprintf(pString, "How many instantiation of MARS (10-5000) ? ");
      numMars_ = getInt(10, 5000, pString);
      sprintf(pString, "How many basis functions in MARS (< %d) ? ",nSamples);
      maxBasis_ = getInt(1, nSamples, pString);
      sprintf(pString, "How many variables per basis (<= %d) ? ",nInputs);
      varPerBasis_ = getInt(1, nInputs, pString);
   }
   else if (outputLevel_ > 1)
   {
      printf("MARS with bagging: default mode = mean.\n");
      printf("                   To change mode, turn on interactive.\n");
      numMars_  = 51;
      printf("MARS with bagging: default number of instantiation = 51.\n");
   }

   strcpy(pString, "mars_params");
   targv[0] = (char *) pString;
   targv[1] = (char *) &maxBasis_;
   targv[2] = (char *) &varPerBasis_;
   itmp = psRSExpertMode_;
   psRSExpertMode_ = 0;
   marsObjs_ = new Mars*[numMars_];
   for (ii = 0; ii < numMars_; ii++) 
   {
      marsObjs_[ii] = new Mars(nInputs_, nSamples_);
      marsObjs_[ii]->setParams(3, targv);
   }
   psRSExpertMode_ = itmp;

   dataSetX_ = NULL;
   dataSetY_ = NULL;

#else
   printf("PSUADE ERROR : MARSBAGG not installed.\n");
   exit(1);
#endif
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
MarsBagg::~MarsBagg()
{
   int ii;

   if (marsObjs_ != NULL) 
   {
      for (ii = 0; ii < numMars_; ii++) delete marsObjs_[ii];
      delete [] marsObjs_;
   }
}

// ************************************************************************
// Set bounds for object class FuncApprox
// ------------------------------------------------------------------------
int MarsBagg::setBounds( double *lower, double *upper )
{
   for (int ii=0 ; ii<numMars_; ii++) marsObjs_[ii]->setBounds(lower, upper);
   return 0;
}

// ************************************************************************
// load output weights
// ------------------------------------------------------------------------
int MarsBagg::loadWeights(int n, double *wgts)
{
   for (int ii=0 ; ii<numMars_; ii++) marsObjs_[ii]->loadWeights(n, wgts);
   return 0;
}

// ************************************************************************
// set number of points to generate in each dimension
// ------------------------------------------------------------------------
void MarsBagg::setNPtsPerDim(int npoints)
{
   nPtsPerDim_ = npoints;
   for (int ii=0 ; ii<numMars_; ii++) marsObjs_[ii]->setNPtsPerDim(npoints);
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int MarsBagg::genNDGridData(double *XX, double *Y, int *N, double **XX2, 
                            double **Y2)
{
#ifdef HAVE_MARS
   int    totPts, ii, ss, jj, index, SAFlag, *iArray, *iCnts;
   double *XXt, *Yt, *XB, *YB, **YM, **SAIndices, *means, *stdevs;

   XB = new double[nInputs_ * nSamples_];
   YB = new double[nSamples_];
   if (outputLevel_ >= 4)
   {
      SAIndices = new double*[numMars_];
      iArray    = new int[nInputs_];
      for (ii = 0; ii < numMars_; ii++)
      {
         SAIndices[ii] = new double[nInputs_];
         for (jj = 0; jj < nInputs_; jj++) SAIndices[ii][jj] = 0.0;
      }
   }
   else
   {
      SAIndices = NULL;
      iArray = NULL;
   }
   SAFlag = 0;
   if ((*N) == -999)
   {
      iCnts = new int[nSamples_];
      for (ii = 0; ii < numMars_; ii++)
      {
         if (outputLevel_ >= 4)
            printf("MarsBagg::genNDGridData : creating Mars #%d (of %d)\n",
                   ii+1, numMars_);
         if (dataSetX_ == NULL)
         {
            for (ss = 0; ss < nSamples_; ss++) iCnts[ss] = 2;
            for (ss = 0; ss < nSamples_; ss++)
            {
               index = -1;
               while (index == -1)
               {
                  index = PSUADE_rand() % nSamples_;
                  if (iCnts[index] <= 0) index = -1; 
               }
               for (jj = 0; jj < nInputs_; jj++)
                  XB[ss*nInputs_+jj] = XX[index*nInputs_+jj]; 
               YB[ss] = Y[index]; 
               iCnts[index]--;
            }
         }
         else
         {
            for (ss = 0; ss < nSamples_; ss++)
            {
               for (jj = 0; jj < nInputs_; jj++)
                  XB[ss*nInputs_+jj] = dataSetX_[ii][ss*nInputs_+jj]; 
               YB[ss] = dataSetY_[ii][ss]; 
            }
         }
         marsObjs_[ii]->setOutputLevel(outputLevel_);
         marsObjs_[ii]->genNDGridData(XB, YB, N, NULL, NULL);
         if (outputLevel_ > 4) 
            SAFlag += getImportance(nInputs_, SAIndices[ii]);
      }
      delete [] iCnts;
      if (SAFlag == 0 && SAIndices != NULL)
      {
         means  = new double[nInputs_];
         stdevs = new double[nInputs_];
         for (jj = 0; jj < nInputs_; jj++)
         {
            means[jj] = stdevs[jj] = 0.0;
            for (ii = 0; ii < numMars_; ii++) means[jj] += SAIndices[ii][jj];
            means[jj]/= (double) numMars_;
            for (ii = 0; ii < numMars_; ii++)
               stdevs[jj] += pow(SAIndices[ii][jj] - means[jj], 2.0e0);
            stdevs[jj] /= (double) numMars_;
            stdevs[jj] = sqrt(stdevs[jj]);
         }     
         for (ii = 0; ii < nInputs_; ii++) iArray[ii] = ii;
         sortDbleList2a(nInputs_, means, iArray);
         printf("* ========== MARS screening rankings =========== *\n");
         for (ii = nInputs_-1; ii >= 0; ii--)
            printf("*  Rank %3d : Input = %3d (measure = %9.3e, stdev = %9.3e)\n",
                   nInputs_-ii, iArray[ii]+1, means[ii], stdevs[ii]);
         printf("* ============================================== *\n");
         delete [] means;
         delete [] stdevs;
      }
   }
   else
   {
      totPts = nPtsPerDim_;
      for (ii = 1; ii < nInputs_; ii++) totPts = totPts * nPtsPerDim_;
      (*XX2) = new double[nInputs_ * totPts];
      (*Y2)  = new double[totPts];
      for (ii = 0; ii < totPts; ii++) (*Y2)[ii] = 0.0;

      if (mode_ != 0)
      {
         YM = new double*[totPts];
         for (ss = 0; ss < totPts; ss++) YM[ss] = new double[numMars_];
      }

      for (ii = 0; ii < numMars_; ii++)
      {
         for (ss = 0; ss < nSamples_; ss++)
         {
            index = PSUADE_rand() % nSamples_;
            for (jj = 0; jj < nInputs_; jj++)
               XB[ss*nInputs_+jj] = XX[index*nInputs_+jj]; 
            YB[ss] = Y[index]; 
         }
         marsObjs_[ii]->genNDGridData(XB, YB, N, &XXt, &Yt);

         if (outputLevel_ >= 4) 
            SAFlag += getImportance(nInputs_, SAIndices[ii]);

         if (mode_ != 0)
         {
            for (ss = 0; ss < totPts; ss++) YM[ss][ii] = Yt[ss];
         }
         else
         {
            for (ss = 0; ss < totPts; ss++) (*Y2)[ss] += Yt[ss];
         }

         if (ii == numMars_-1)
            for (ss = 0; ss < totPts*nInputs_; ss++) (*XX2)[ss] = XXt[ss];
         delete [] XXt;
         delete [] Yt;
      }

      if (SAFlag == 0 && SAIndices != NULL)
      {
         means  = new double[nInputs_];
         stdevs = new double[nInputs_];
         for (jj = 0; jj < nInputs_; jj++)
         {
            means[jj] = stdevs[jj] = 0.0;
            for (ii = 0; ii < numMars_; ii++) means[jj] += SAIndices[ii][jj];
            means[jj]/= (double) numMars_;
            for (ii = 0; ii < numMars_; ii++)
               stdevs[jj] += pow(SAIndices[ii][jj] - means[jj], 2.0e0);
            stdevs[jj] /= (double) numMars_;
            stdevs[jj] = sqrt(stdevs[jj]);
         }     
         for (ii = 0; ii < nInputs_; ii++) iArray[ii] = ii;
         sortDbleList2a(nInputs_, means, iArray);
         printf("* ========= MARSBAG screening rankings ========= *\n");
         for (ii = nInputs_-1; ii >= 0; ii--)
            printf("*  Rank %3d : Input = %3d (measure = %9.3e, stdev = %9.3e)\n",
                   nInputs_-ii, iArray[ii]+1, means[ii], stdevs[ii]);
         printf("* ============================================== *\n");
         delete [] means;
         delete [] stdevs;
      }

      (*N) = totPts;

      if (mode_ != 0)
      {
         for (ss = 0; ss < totPts; ss++)
         {
            sortDbleList(numMars_, YM[ss]); 
            (*Y2)[ss] = YM[ss][numMars_/2];
            delete [] YM[ss];
         }
         delete [] YM;
      }
      else
      { 
         for (ss = 0; ss < totPts; ss++) (*Y2)[ss] /= (double) numMars_;
      }
   }
   if (SAIndices != NULL)
   {
      for (ii = 0; ii < numMars_; ii++) delete [] SAIndices[ii];
      delete [] SAIndices;
      delete [] iArray;
   }
   delete [] XB;
   delete [] YB;
#else
   printf("PSUADE ERROR : MARS not used.\n");
   return -1;
#endif
   return 0;
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int MarsBagg::gen1DGridData(double *X, double *Y, int ind1, 
                            double *settings, int *N, double **XX, 
                            double **YY)
{
#ifdef HAVE_MARS
   int    totPts, ii, ss, jj, index;
   double *XB, *YB, *XXt, *Yt, **YM;

   XB = new double[nInputs_ * nSamples_];
   YB = new double[nSamples_];
   totPts = nPtsPerDim_;
   (*N)   = totPts;
   (*XX)  = new double[nInputs_ * totPts];
   (*YY)  = new double[totPts];
   for (ii = 0; ii < totPts; ii++) (*YY)[ii] = 0.0;

   if (mode_ != 0)
   {
      YM = new double*[totPts];
      for (ss = 0; ss < totPts; ss++) YM[ss] = new double[numMars_];
   }

   for (ii = 0; ii < numMars_; ii++)
   {
      if (dataSetX_ == NULL)
      {
         for (ss = 0; ss < nSamples_; ss++)
         {
            index = PSUADE_rand() % nSamples_;
            for (jj = 0; jj < nInputs_; jj++)
               XB[ss*nInputs_+jj] = X[index*nInputs_+jj]; 
            YB[ss] = Y[index]; 
         }
      }
      else
      {
         for (ss = 0; ss < nSamples_; ss++)
         {
            for (jj = 0; jj < nInputs_; jj++)
               XB[ss*nInputs_+jj] = dataSetX_[ii][ss*nInputs_+jj]; 
            YB[ss] = dataSetY_[ii][ss]; 
         }
      }
      marsObjs_[ii]->gen1DGridData(XB, YB, ind1, settings, N, &XXt, &Yt);

      if (mode_ != 0)
      {
         for (ss = 0; ss < totPts; ss++) YM[ss][ii] = Yt[ss];
      }
      else
      {
         for (ss = 0; ss < totPts; ss++) (*YY)[ss] += Yt[ss];
      }
      if (ii == numMars_-1)
         for (ss = 0; ss < totPts*nInputs_; ss++) (*XX)[ss] = XXt[ss];
      delete [] XXt;
      delete [] Yt;
   }

   if (mode_ != 0)
   {
      for (ss = 0; ss < totPts; ss++)
      {
         sortDbleList(numMars_, YM[ss]); 
         (*YY)[ss] = YM[ss][numMars_/2];
         delete [] YM[ss];
      }
      delete [] YM;
   }
   else
   { 
      for (ss = 0; ss < totPts; ss++) (*YY)[ss] /= (double) numMars_;
   }
   delete [] XB;
   delete [] YB;
#else
   printf("PSUADE ERROR : MARS not used.\n");
   return -1;
#endif
   return 0;
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int MarsBagg::gen2DGridData(double *X, double *Y, int ind1, int ind2, 
                            double *settings, int *N, double **XX, 
                            double **YY)
{
#ifdef HAVE_MARS
   int    totPts, ii, ss, jj, index;
   double *XB, *YB, *XXt, *Yt, **YM;

   XB = new double[nInputs_ * nSamples_];
   YB = new double[nSamples_];
   totPts = nPtsPerDim_ * nPtsPerDim_;
   (*N)   = totPts;
   (*XX)  = new double[nInputs_ * totPts];
   (*YY)  = new double[totPts];
   for (ii = 0; ii < totPts; ii++) (*YY)[ii] = 0.0;

   if (mode_ != 0)
   {
      YM = new double*[totPts];
      for (ss = 0; ss < totPts; ss++) YM[ss] = new double[numMars_];
   }

   for (ii = 0; ii < numMars_; ii++)
   {
      if (outputLevel_ >= 1)
         printf("MarsBagg::gen2DGridData : creating Mars #%d (of %d)\n",
                ii+1, numMars_);
      if (dataSetX_ == NULL)
      {
         for (ss = 0; ss < nSamples_; ss++)
         {
            index = PSUADE_rand() % nSamples_;
            for (jj = 0; jj < nInputs_; jj++)
               XB[ss*nInputs_+jj] = X[index*nInputs_+jj]; 
            YB[ss] = Y[index]; 
         }
      }
      else
      {
         for (ss = 0; ss < nSamples_; ss++)
         {
            for (jj = 0; jj < nInputs_; jj++)
               XB[ss*nInputs_+jj] = dataSetX_[ii][ss*nInputs_+jj]; 
            YB[ss] = dataSetY_[ii][ss]; 
         }
      }
      marsObjs_[ii]->gen2DGridData(XB, YB, ind1, ind2, settings, N, 
                                   &XXt, &Yt);

      if (mode_ != 0)
      {
         for (ss = 0; ss < totPts; ss++) YM[ss][ii] = Yt[ss];
      }
      else
      {
         for (ss = 0; ss < totPts; ss++) (*YY)[ss] += Yt[ss];
      }

      if (ii == numMars_-1)
         for (ss = 0; ss < totPts*nInputs_; ss++) (*XX)[ss] = XXt[ss];
      delete [] XXt;
      delete [] Yt;
   }

   if (mode_ != 0)
   {
      for (ss = 0; ss < totPts; ss++)
      {
         sortDbleList(numMars_, YM[ss]); 
         (*YY)[ss] = YM[ss][numMars_/2];
         delete [] YM[ss];
      }
      delete [] YM;
   }
   else
   { 
      for (ss = 0; ss < totPts; ss++) (*YY)[ss] /= (double) numMars_;
   }
   delete [] XB;
   delete [] YB;
#else
   printf("PSUADE ERROR : MARS not used.\n");
   return -1;
#endif
   return 0;
}

// ************************************************************************
// Generate 3D results for display
// ------------------------------------------------------------------------
int MarsBagg::gen3DGridData(double *X, double *Y, int ind1, int ind2, 
                            int ind3, double *settings, int *N, double **XX, 
                            double **YY)
{
#ifdef HAVE_MARS
   int    totPts, ii, ss, jj, index;
   double *XB, *YB, *XXt, *Yt, **YM;

   XB = new double[nInputs_ * nSamples_];
   YB = new double[nSamples_];
   totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
   (*N)   = totPts;
   (*XX)  = new double[nInputs_ * totPts];
   (*YY)  = new double[totPts];
   for (ii = 0; ii < totPts; ii++) (*YY)[ii] = 0.0;

   if (mode_ != 0)
   {
      YM = new double*[totPts];
      for (ss = 0; ss < totPts; ss++) YM[ss] = new double[numMars_];
   }

   for (ii = 0; ii < numMars_; ii++)
   {
      if (dataSetX_ == NULL)
      {
         for (ss = 0; ss < nSamples_; ss++)
         {
            index = PSUADE_rand() % nSamples_;
            for (jj = 0; jj < nInputs_; jj++)
               XB[ss*nInputs_+jj] = X[index*nInputs_+jj]; 
            YB[ss] = Y[index]; 
         }
      }
      else
      {
         for (ss = 0; ss < nSamples_; ss++)
         {
            for (jj = 0; jj < nInputs_; jj++)
               XB[ss*nInputs_+jj] = dataSetX_[ii][ss*nInputs_+jj]; 
            YB[ss] = dataSetY_[ii][ss]; 
         }
      }
      marsObjs_[ii]->gen3DGridData(XB, YB, ind1, ind2, ind3, settings, 
                                   N, &XXt, &Yt);

      if (mode_ != 0)
      {
         for (ss = 0; ss < totPts; ss++) YM[ss][ii] = Yt[ss];
      }
      else
      {
         for (ss = 0; ss < totPts; ss++) (*YY)[ss] += Yt[ss];
      }

      if (ii == numMars_-1)
         for (ss = 0; ss < totPts*nInputs_; ss++) (*XX)[ss] = XXt[ss];
      delete [] XXt;
      delete [] Yt;
   }

   if (mode_ != 0)
   {
      for (ss = 0; ss < totPts; ss++)
      {
         sortDbleList(numMars_, YM[ss]); 
         (*YY)[ss] = YM[ss][numMars_/2];
         delete [] YM[ss];
      }
      delete [] YM;
   }
   else
   { 
      for (ss = 0; ss < totPts; ss++) (*YY)[ss] /= (double) numMars_;
   }
   delete [] XB;
   delete [] YB;
#else
   printf("PSUADE ERROR : MARS not used.\n");
   return -1;
#endif
   return 0;
}

// ************************************************************************
// Generate 4D results for display
// ------------------------------------------------------------------------
int MarsBagg::gen4DGridData(double *X, double *Y, int ind1, int ind2, 
                            int ind3, int ind4, double *settings, int *N, 
                            double **XX, double **YY)
{
#ifdef HAVE_MARS
   int    totPts, ii, ss, jj, index;
   double *XB, *YB, *XXt, *Yt, **YM;

   XB = new double[nInputs_ * nSamples_];
   YB = new double[nSamples_];
   totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
   (*N)   = totPts;
   (*XX)  = new double[nInputs_ * totPts];
   (*YY)  = new double[totPts];
   for (ii = 0; ii < totPts; ii++) (*YY)[ii] = 0.0;

   if (mode_ != 0)
   {
      YM = new double*[totPts];
      for (ss = 0; ss < totPts; ss++) YM[ss] = new double[numMars_];
   }

   for (ii = 0; ii < numMars_; ii++)
   {
      if (dataSetX_ == NULL)
      {
         for (ss = 0; ss < nSamples_; ss++)
         {
            index = PSUADE_rand() % nSamples_;
            for (jj = 0; jj < nInputs_; jj++)
               XB[ss*nInputs_+jj] = X[index*nInputs_+jj]; 
            YB[ss] = Y[index]; 
         }
      }
      else
      {
         for (ss = 0; ss < nSamples_; ss++)
         {
            for (jj = 0; jj < nInputs_; jj++)
               XB[ss*nInputs_+jj] = dataSetX_[ii][ss*nInputs_+jj]; 
            YB[ss] = dataSetY_[ii][ss]; 
         }
      }
      marsObjs_[ii]->gen4DGridData(XB, YB, ind1, ind2, ind3, ind4,
                                   settings, N, &XXt, &Yt);

      if (mode_ != 0)
      {
         for (ss = 0; ss < totPts; ss++) YM[ss][ii] = Yt[ss];
      }
      else
      {
         for (ss = 0; ss < totPts; ss++) (*YY)[ss] += Yt[ss];
      }

      if (ii == numMars_-1)
         for (ss = 0; ss < totPts*nInputs_; ss++) (*XX)[ss] = XXt[ss];
      delete [] XXt;
      delete [] Yt;
   }

   if (mode_ != 0)
   {
      for (ss = 0; ss < totPts; ss++)
      {
         sortDbleList(numMars_, YM[ss]); 
         (*YY)[ss] = YM[ss][numMars_/2];
         delete [] YM[ss];
      }
      delete [] YM;
   }
   else
   { 
      for (ss = 0; ss < totPts; ss++) (*YY)[ss] /= (double) numMars_;
   }
   delete [] XB;
   delete [] YB;
#else
   printf("PSUADE ERROR : MARS not used.\n");
   return -1;
#endif
   return 0;
}

// ************************************************************************
// Evaluate a given point 
// ------------------------------------------------------------------------
double MarsBagg::evaluatePoint(double *X)
{
   double Y=0.0;
#ifdef HAVE_MARS
   int    ii;
   double Yt, *YM;

   if (mode_ != 0) YM = new double[numMars_];

   for (ii = 0; ii < numMars_; ii++) 
   {
      Yt = marsObjs_[ii]->evaluatePoint(X);
      if (mode_ != 0) YM[ii] = Yt;
      else            Y += Yt;
   }
printf("\n");
   if (mode_ != 0) Y = YM[numMars_/2];
   else            Y /= (double) numMars_;
#else
   printf("PSUADE ERROR : MARS not used.\n");
#endif
   return Y;
}

// ************************************************************************
// Evaluate a number of points
// ------------------------------------------------------------------------
double MarsBagg::evaluatePoint(int npts, double *X, double *Y)
{
#ifdef HAVE_MARS
   int    in, ii;
   double YY, Yt, *YM;
FILE *fp = fopen("marsfile", "a");

   if (mode_ != 0) YM = new double[numMars_];

   for (in = 0; in < npts; in++) 
   {
      YY = 0.0;
      for (ii = 0; ii < numMars_; ii++) 
      {
         Yt = marsObjs_[ii]->evaluatePoint(&(X[in*nInputs_]));
         if (mode_ != 0) YM[ii] = Yt;
         else            YY += Yt;
fprintf(fp, "%14.6e ", Yt);
      }
fprintf(fp, "\n");
      if (mode_ != 0) Y[in] = YM[numMars_/2];
      else            Y[in] = YY / (double) numMars_;
   }
fclose(fp);
#else
   printf("PSUADE ERROR : MARS not used.\n");
#endif
   return 0.0;
}

// ************************************************************************
// Evaluate a given point with standard deviation
// ------------------------------------------------------------------------
double MarsBagg::evaluatePointFuzzy(double *X, double &std)
{
   double Ymean=0.0;
#ifdef HAVE_MARS
   int    ii;
   double *Y;

   Y = new double[numMars_];
   for (ii = 0; ii < numMars_; ii++) 
   {
      Y[ii] = marsObjs_[ii]->evaluatePoint(X);
      Ymean += Y[ii];
   }
   Ymean /= (double) numMars_;
   std = 0.0;
   for (ii = 0; ii < numMars_; ii++) 
      std += (Y[ii] - Ymean) * (Y[ii] - Ymean);
   std = sqrt(std / (double) numMars_);
   delete [] Y;
#else
   printf("PSUADE ERROR : MARS not used.\n");
#endif
   return Ymean;
}

// ************************************************************************
// Evaluate a number of points with standard deviations
// ------------------------------------------------------------------------
double MarsBagg::evaluatePointFuzzy(int npts, double *X, double *Y,
                                    double *Ystd)
{
#ifdef HAVE_MARS
   int    in, ii;
   double *YY;

   YY = new double[numMars_];
   for (in = 0; in < npts; in++) 
   {
      for (ii = 0; ii < numMars_; ii++) 
         YY[ii] = marsObjs_[ii]->evaluatePoint(&(X[in*nInputs_]));
      Y[in] = 0.0;
      for (ii = 0; ii < numMars_; ii++) Y[in] += YY[ii];
      Y[in] /= (double) numMars_;
      Ystd[in] = 0.0;
      for (ii = 0; ii < numMars_; ii++) 
         Ystd[in] += (YY[ii] - Y[in]) * (YY[ii] - Y[in]);
      Ystd[in] = sqrt(Ystd[in] / (double) numMars_);
   }
   delete [] YY;
#else
   printf("PSUADE ERROR : MARS not used.\n");
#endif
   return 0.0;
}

// ************************************************************************
// get the importance indicators
// ------------------------------------------------------------------------
int MarsBagg::getImportance(int nInputs, double *indicators)
{
   int    lineLeng=500, ii, jj, nCount;
   double dmax;
   char   line[500], word1[500], word2[500], word3[500];
   FILE   *fp;

   fp = fopen(".psuade_mars", "r");
   if (fp == NULL) return -1;
   else
   {
      nCount = 0;
      strcpy(word1, "none");
      while ((fgets(line, lineLeng, fp) != NULL) && (feof(fp) == 0))
      {
         sscanf(line,"%s %s %s", word1, word2, word3);
         if (!strcmp(word1,"relative") && !strcmp(word2,"variable") &&
             !strcmp(word3,"importance:"))
         {
            fgets(line, lineLeng, fp);
            fgets(line, lineLeng, fp);
            if (feof(fp) == 0)
            {
               for (ii = 0; ii < nInputs_; ii+=6)
               {
                  for (jj = 0; jj < 6; jj++)
                  {
                     if (ii+jj < nInputs_)
                     {
                        fscanf(fp,"%lg", &(indicators[ii+jj]));
                        nCount++;
                     }
                  }
                  fgets(line, lineLeng, fp);
                  fgets(line, lineLeng, fp);
                  fgets(line, lineLeng, fp);
               }
            }
         }
      }
      if (nCount != nInputs_)
      {
         fclose(fp);
         return -1;
      }
      dmax = indicators[0];
      for (ii = 1; ii < nInputs_; ii++)
         if (indicators[ii] > dmax) dmax = indicators[ii];
      for (ii = 0; ii < nInputs_; ii++) indicators[ii] /= dmax;
      fclose(fp);
   }
   return 0;
}

// ************************************************************************
// set parameters
// ------------------------------------------------------------------------
double MarsBagg::setParams(int targc, char **targv)
{
   int    ii, itmp, leng;
   double *Xdata, *Ydata;
   char   cString[500], *argv[3];

   if (targc == 1 && !strcmp(targv[0], "median"))
   {
      mode_ = 1;
   }
   else if (targc == 2 && !strcmp(targv[0], "num_mars"))
   {
      if (marsObjs_ != NULL) 
      {
         for (ii = 0; ii < numMars_; ii++) delete marsObjs_[ii];
         delete [] marsObjs_;
      }
      numMars_ = *(int *) targv[1];
      if (numMars_ < 2) numMars_ = 2;
      printf("MARS with bagging: no. of MARs set to = %d.\n", numMars_);
      strcpy(cString, "mars_params");
      argv[0] = (char *) cString;
      argv[1] = (char *) &maxBasis_;
      argv[2] = (char *) &varPerBasis_;
      itmp = psRSExpertMode_;
      psRSExpertMode_ = 0;
      marsObjs_ = new Mars*[numMars_];
      for (ii = 0; ii < numMars_; ii++) 
      {
         marsObjs_[ii] = new Mars(nInputs_, nSamples_);
         marsObjs_[ii]->setParams(3, argv);
      }
      psRSExpertMode_ = itmp;
   }
   else if (targc == 5 && !strcmp(targv[0], "mars_sample"))
   {
      itmp = *(int *) targv[1];
      if (itmp < 0 || itmp >= numMars_)
      {
         printf("MarsBag ERROR: in loading sample - invalid index.\n");
         exit(1);
      }
      leng = *(int *) targv[2];
      if (leng != nSamples_)
      {
         printf("MarsBag ERROR: in loading sample - nSamples mismatch.\n");
         exit(1);
      }
      Xdata = (double *) targv[3];
      Ydata = (double *) targv[4];
      if (dataSetX_ == NULL)
      {
         dataSetX_ = new double*[numMars_];
         for (ii = 0; ii < numMars_; ii++) dataSetX_[ii] = NULL;
         dataSetY_ = new double*[numMars_];
         for (ii = 0; ii < numMars_; ii++) dataSetY_[ii] = NULL;
      }
      if (dataSetX_[itmp] == NULL)
      {
         dataSetX_[itmp] = new double[leng*nInputs_];
         dataSetY_[itmp] = new double[leng];
      }
      for (ii = 0; ii < leng*nInputs_; ii++) dataSetX_[itmp][ii] = Xdata[ii];
      for (ii = 0; ii < leng; ii++) dataSetY_[itmp][ii] = Ydata[ii];
   }
   else if (targc == 3 && !strcmp(targv[0], "mars_params"))
   {
      maxBasis_ = *(int *) targv[1];
      varPerBasis_ = *(int *) targv[2];
      printf("MARS with bagging: numBasis    set to = %d.\n", maxBasis_);
      printf("MARS with bagging: varPerBasis set to = %d.\n", varPerBasis_);
   }
   else
   {
      printf("MarsBagg setParams ERROR: invalid command %s.\n", targv[0]);
   }
   return 0.0;
}

