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
// Functions for the class Kriging
// AUTHOR : CHARLES TONG
// DATE   : 2012
// ************************************************************************
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "Kriging.h"
#include "sysdef.h"
#include "Psuade.h"
#include "PsuadeUtil.h"
#include "PsuadeConfig.h"
#include "Sampling.h"

#define PABS(x) ((x) > 0 ? (x) : (-(x)))

// ************************************************************************
// ************************************************************************
// external functions
// ------------------------------------------------------------------------
extern "C" {
   void dpotrf_(char *, int *, double *, int *, int *);
   void dgetrf_(int *, int *, double *, int *, int *, int *);
   void dpotrs_(char *, int *, int *, double *, int *, double *,
                int *, int *);
   void dgetrs_(char *, int *, int *, double *Cmatrix_, int *nRows, int *LUPivots_,
                double *, int *, int *);
   void bobyqa_(int *,int *, double *, double *, double *, double *,
                double *, int *, int *, double*);
}

// ************************************************************************
// ************************************************************************
// internal global variables
// ------------------------------------------------------------------------
int    KRI_outputLevel=0;
int    KRI_iter=-1;
int    KRI_nInputs=-1;
int    KRI_nSamples=-1;
int    KRI_pOrder=-1;
double *KRI_XDists=NULL;
double *KRI_SMatrix=NULL;
double *KRI_FMatrix=NULL;
double *KRI_FMatTmp=NULL;
double *KRI_MMatrix=NULL;
double *KRI_X=NULL;
double *KRI_Y=NULL;
double *KRI_Ytmp=NULL;
double KRI_OptY=0.0;
double *KRI_OptThetas=NULL;
double *KRI_dataStdDevs=NULL;
double KRI_YStd=1.0;

// ************************************************************************
// ************************************************************************
// resident function to perform evaluation 
// ------------------------------------------------------------------------
#ifdef __cplusplus
extern "C"
{
#endif
   void *kribobyqaevalfunc_(int *nInps, double *XValues, double *YValue)
   {
      int    nBasis, count, Cleng, ii, jj, kk, status;
      double dist, ddata, ddata2;
      char   uplo='L';
      FILE   *fp=NULL;

      fp = fopen("ps_terminate", "r");
      if (fp != NULL)
      {
         fclose(fp);
         printf("Kriging: terminating ....\n");
         (*YValue) = KRI_OptY;
         return NULL;
      }
      KRI_iter++;
      nBasis = 1;
      if (KRI_pOrder == 1) nBasis = KRI_nInputs + 1;
      Cleng = KRI_nSamples + nBasis;

      count = 0;
      for (jj = 0; jj < KRI_nSamples; jj++)
      {
         KRI_SMatrix[jj*KRI_nSamples+jj] = 1.0 + 1e-15 * KRI_nSamples;
         if (KRI_dataStdDevs != NULL) 
            KRI_SMatrix[jj*KRI_nSamples+jj] += pow(KRI_dataStdDevs[jj]/KRI_YStd,2.0);
         for (kk = jj+1; kk < KRI_nSamples; kk++)
         {
            dist = 0.0;
            for (ii = 0; ii < KRI_nInputs; ii++)
               dist += pow(KRI_XDists[count*KRI_nInputs+ii]/XValues[ii], 2.0);
            dist = exp(-dist);
            if (dist < 1.0e-16) dist = 0;
            KRI_SMatrix[jj*KRI_nSamples+kk] = dist;
            KRI_SMatrix[kk*KRI_nSamples+jj] = dist;
            count++;
         }
      }

      for (jj = 0; jj < KRI_nSamples; jj++)
      {
         KRI_FMatrix[jj] = 1.0;
         KRI_FMatTmp[jj] = 1.0;
         if (KRI_pOrder == 1)
         {
            for (kk = 1; kk < KRI_nInputs+1; kk++)
            {
               ddata = KRI_X[jj*KRI_nInputs+kk-1];
               KRI_FMatTmp[kk*KRI_nSamples+jj] = ddata;
               KRI_FMatrix[kk*KRI_nSamples+jj] = ddata;
            }
         }
      }
      dpotrf_(&uplo, &KRI_nSamples, KRI_SMatrix, &KRI_nSamples, &status);
      kk = nBasis;
      dpotrs_(&uplo, &KRI_nSamples, &kk, KRI_SMatrix, &KRI_nSamples,
              KRI_FMatTmp, &KRI_nSamples, &status);
      for (jj = 0; jj < nBasis; jj++)
      {
         for (kk = 0; kk < nBasis; kk++)
         {
            ddata = 0.0;
            for (ii = 0; ii < KRI_nSamples; ii++)
               ddata += KRI_FMatrix[ii+jj*KRI_nSamples]*
                        KRI_FMatTmp[ii+kk*KRI_nSamples];
            KRI_MMatrix[jj+kk*nBasis] = ddata;
         }
      }
      if (nBasis > 1) dpotrf_(&uplo, &nBasis, KRI_MMatrix, &nBasis, &status);
      for (jj = 0; jj < KRI_nSamples; jj++) KRI_Ytmp[jj] = KRI_Y[jj];
      kk = 1;
      dpotrs_(&uplo,&KRI_nSamples,&kk,KRI_SMatrix,&KRI_nSamples,KRI_Ytmp,
              &KRI_nSamples, &status);
      for (jj = 0; jj < nBasis; jj++)
      {
         ddata = 0.0;
         for (ii = 0; ii < KRI_nSamples; ii++)
            ddata += KRI_FMatrix[ii+jj*KRI_nSamples]*KRI_Ytmp[ii];
         KRI_Ytmp[KRI_nSamples+jj] = ddata;
      }
      if (nBasis == 1)
      {
         if (KRI_MMatrix[0] == 0)
         {
            printf("Kriging ERROR: divide by 0.\n");
            exit(1);
         }
         KRI_Ytmp[KRI_nSamples] /= KRI_MMatrix[0];
      }
      else
      {
         kk = 1;
         dpotrs_(&uplo,&kk,&kk,KRI_MMatrix,&nBasis,&KRI_Ytmp[KRI_nSamples],
                 &kk, &status);
      }
      for (jj = 0; jj < KRI_nSamples; jj++)
      {
         ddata = 0.0;
         for (ii = 0; ii < nBasis; ii++)
            ddata += KRI_FMatrix[jj+ii*KRI_nSamples]*KRI_Ytmp[ii+KRI_nSamples];
         KRI_Ytmp[jj] = KRI_Y[jj] - ddata;
      }
      kk = 1;
      dpotrs_(&uplo,&KRI_nSamples,&kk,KRI_SMatrix,&KRI_nSamples,KRI_Ytmp,
              &KRI_nSamples, &status);

      ddata = 0.0;
      for (jj = 0; jj < KRI_nSamples; jj++)
         ddata += KRI_Ytmp[jj] * KRI_Y[jj];
      ddata /= (double) KRI_nSamples;

      for (jj = 0; jj < KRI_nSamples; jj++) 
      {
         ddata2 = KRI_SMatrix[KRI_nSamples*jj+jj];
         if (ddata2 < 0.0) ddata2 = - ddata2;
         ddata *= pow(ddata2, 2.0/(double) KRI_nSamples);
      }

      (*YValue) = ddata;

      if (ddata < KRI_OptY)
      {
         KRI_OptY = ddata;
         for (ii = 0; ii < KRI_nInputs; ii++)
            KRI_OptThetas[ii] = XValues[ii];
         if (KRI_outputLevel > 1)
         {
            printf("\t Kriging (2): iteration %d\n", KRI_iter);
            for (ii = 0; ii < KRI_nInputs; ii++)
               printf("\t    Current best theta %d = %e\n",ii+1,XValues[ii]);
            printf("\t    Current best objective value so far = %e\n",ddata);
         }
      }
      if (psRSExpertMode_ == 1 && KRI_outputLevel > 3)
      {
         for (ii = 0; ii < KRI_nInputs; ii++)
            printf("\t    theta %d = %e\n", ii+1, XValues[ii]);
         printf("\t    Current Kriging objective value = %e\n", ddata);
         printf("\t* To terminate optimization early, create a file\n");
         printf("\t* called 'ps_terminate' in your local directory.\n");
         fp = fopen("ps_terminate", "r");
         if (fp != NULL)
         {
            printf("\t*** ps_terminate file found.\n");
            (*YValue) = 0.0;
            fclose(fp);
         }
      }
      return NULL;
   }
#ifdef __cplusplus
}
#endif

// ************************************************************************
// Constructor for object class Kriging
// ------------------------------------------------------------------------
Kriging::Kriging(int nInputs,int nSamples) : FuncApprox(nInputs,nSamples)
{
   int    ii, jj;
   double ddata;
   char   pString[500], winput[500], winput2[500], fname[500], *strPtr;
   FILE   *fp;

   faID_ = PSUADE_RS_KR;
   XNormalized_ = NULL;
   YNormalized_ = NULL;
   Cmatrix_ = NULL;
   pOrder_  = 0;
   initFlag_ = 0;
   workLength_ = 0;
   workArray_ = NULL;
   workX_ = NULL;
   dataStdDevs_ = NULL;
   optTolerance_ = 1.0e-6;
   fastMode_ = 2;
   Thetas_ = new double[nInputs_];
   for (ii = 0; ii < nInputs_; ii++) Thetas_[ii] = 0.1;

   printAsterisks(0);
   printf("*                Kriging Analysis\n");
   printf("* Set printlevel to 1-4 to see Kriging details.\n");
   printf("* Turn on rs_expert mode to set slow or fast mode.\n");
   printf("* Fast mode: no optimization of hyperparameters.\n");
   printf("*      - turn on rs_expert to set hyperparameters.\n");
   printf("*      - default values = 1.0\n");
   printf("* Slow mode (default): hyperparameters are optimized.\n");
   printf("*      - to change optimization parameters, turn\n");
   printf("*        rs_expert mode.\n");
   printf("* Very slow mode: use multi-start global optimization.\n");
   printf("*      - to change optimization parameters, turn\n");
   printf("*        rs_expert mode.\n");
   printEquals(0); 

   // read configure file 
   if (psRSExpertMode_ == 0 && psConfig_ != NULL)
   {
      strPtr = psConfig_->getParameter("KRI_MODE");
      if (strPtr != NULL)
      {
         sscanf(strPtr, "%s %s %d", winput, winput2, &ii);
         if (ii < 1 || ii > 3)
         {
            printf("Kriging INFO: mode from config file not valid.\n");
            printf("              mode kept at %d.\n", fastMode_);
         }
         else
         {
            fastMode_ = ii;
            printf("Kriging INFO: mode from config file = %d.\n",
                   fastMode_);
         }
      }
      strPtr = psConfig_->getParameter("KRI_TOL");
      if (strPtr != NULL)
      {
         sscanf(strPtr, "%s %s %lg", winput, winput2, &optTolerance_);
         if (optTolerance_ < 0.0 || optTolerance_ >= 1.0)
         {
            optTolerance_ = 1.0e-6;
            printf("Kriging INFO: tolerance from config file not valid.\n");
            printf("              tolerance kept at %d.\n", optTolerance_);
         }
         else
         {
            printf("Kriging INFO: tolerance from config file = %e.\n",
                   optTolerance_);
         }
      }
      strPtr = psConfig_->getParameter("KRI_LENG_SCALE");
      if (strPtr != NULL)
      {
         sscanf(strPtr, "%s %d %s %lg", winput, &ii, winput2, &ddata);
         if (ii < 1 || ii > nInputs_)
         {
            optTolerance_ = 1.0e-6;
            printf("Kriging INFO: invalid input number for length scale.\n");
            printf("              Input number read = %d.\n", ii);
         }
         else
         {
            Thetas_[ii-1] = ddata;
            printf("Kriging INFO: length scale for input %d set to %e.\n",
                   ii, ddata);
         }
      }
      strPtr = psConfig_->getParameter("KRI_DATA_STDEV_FILE");
      if (strPtr != NULL)
      {
         sscanf(strPtr, "%s %s %s", winput, winput2, fname);
         fp = fopen(fname, "r");
         if (fp == NULL)
         {
            printf("Kriging INFO: data variance file not found.\n");
         }
         else
         {
            fscanf(fp, "%d", &ii); 
            if (ii != nSamples_)
            {
               printf("Kriging ERROR: std. dev. file should have %d entries.\n",
                      nSamples_);
               fclose(fp);
            }
            else
            {
               dataStdDevs_ = new double[nSamples_];
               for (ii = 0; ii < nSamples_; ii++)
               {
                  fscanf(fp, "%d %lg", &jj, &dataStdDevs_[ii]); 
                  if (jj != ii+1)
                  {
                     printf("Kriging ERROR: std. dev. file has problem in line %d.\n",
                            jj+1);
                     delete [] dataStdDevs_;
                     dataStdDevs_ = NULL;
                     break;
                  }
               } 
               fclose(fp);
               if (dataStdDevs_ != NULL)
               {
                  printf("Kriging INFO: std. dev. file has been read.\n");
                  KRI_dataStdDevs = dataStdDevs_; 
               }
            }
         }
      }
   }
   
   // if configure file is not used, ask for parameters 
   if (psRSExpertMode_ == 1)
   {
      printf("There are three modes available: \n");
      printf("(1) fast mode with pre-specified thetas\n");
      printf("(2) slow mode with optimization on the thetas\n");
      printf("(3) very slow mode with multi-start optimization\n");
      sprintf(pString, "Please select mode (1, 2, or 3) : ");
      fastMode_ = getInt(1,3,pString);
      if (fastMode_ == 1)
      {
         printf("Kriging: Length scales are correlation lengths\n");
         printf("         in the random parameter space.\n");
         printf("Current initial length scales are:\n");
         for (ii = 0; ii < nInputs_; ii++)
            printf("Input %d: %e\n", ii+1, Thetas_[ii]);
         sprintf(pString, "Change length scales (thetas)? (y or n) ");
         getString(pString, winput);
         if (winput[0] == 'y')
         {
            for (ii = 0; ii < nInputs_; ii++)
            {
               sprintf(pString,"Enter theta for input %d (>0): ", ii+1);
               Thetas_[ii] = getDouble(pString);
               if (Thetas_[ii] <= 0.0)
               {
                  printf("ERROR: theta <= 0 not valid.\n");
                  exit(1);
               }
               if (ii == 0)
               {
                  sprintf(pString,"Use %e for all other thetas? (y or n) ",
                          Thetas_[0]);
                  getString(pString, winput);
                  if (winput[0] == 'y')
                  {
                     for (jj = 1; jj < nInputs_; jj++) Thetas_[jj] = Thetas_[0];
                     break;
                  }
               }
            }
         }
      }
      else
      {
         sprintf(pString, "Enter optimization tolerance (default = 1e-6): ");
         optTolerance_ = getDouble(pString);
         if (optTolerance_ <= 0 || optTolerance_ > 0.5)
         {
            printf("Kriging INFO: optimization tolerance should be in (0,0.5]).\n");
            printf("              Tolerance set to default = 1.0e-6.\n");
            optTolerance_ = 1.0e-6;
         }
         if (fastMode_ == 2)
         {
            printf("Kriging: Current initial length scales (thetas) are:\n");
            for (ii = 0; ii < nInputs_; ii++)
               printf("Input %d: %e\n", ii+1, Thetas_[ii]);
            printf("If some knowledge is available about the relative\n");
            printf("importance of some parameters, different initial\n");
            printf("thetas can be entered to reflect this knowledge \n");
            printf("(sensitive parameters have larger thetas.)\n");
            sprintf(pString, "Change initial thetas? (y or n) ");
            getString(pString, winput);
            if (winput[0] == 'y')
            {
               for (ii = 0; ii < nInputs_; ii++)
               {
                  sprintf(pString,"Enter theta for input %d : ", ii+1);
                  Thetas_[ii] = getDouble(pString);
                  if (Thetas_[ii] <= 0.0)
                     printf("Warning: theta < 0 not recommended.\n");
                  if (ii == 0)
                  {
                     sprintf(pString,"Use %e for all other thetas? (y or n) ",
                             Thetas_[0]);
                     getString(pString, winput);
                     if (winput[0] == 'y')
                     {
                        for (jj = 1; jj < nInputs_; jj++) Thetas_[jj] = Thetas_[0];
                        break;
                     }
                  }
               }
            }
         }
      }
   }
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
Kriging::~Kriging()
{
   if (Thetas_ != NULL) delete [] Thetas_;
   if (XNormalized_ != NULL) delete [] XNormalized_;
   if (YNormalized_ != NULL) delete [] YNormalized_;
   if (Cmatrix_ != NULL) delete [] Cmatrix_;
   if (workArray_ != NULL) delete [] workArray_;
   if (workX_ != NULL) delete [] workX_;
   if (dataStdDevs_ != NULL) delete [] dataStdDevs_;
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int Kriging::genNDGridData(double *X, double *Y, int *N2, double **X2,
                           double **Y2)
{
   int totPts;

   if (outputLevel_ >= 1) printf("Kriging training begins....\n");
   train(X,Y);
   if (outputLevel_ >= 1) printf("Kriging training completed.\n");
   if ((*N2) == -999 || X2 == NULL || Y2 == NULL) return 0;
  
   genNDGrid(N2, X2);
   if ((*N2) == 0) return 0;
   totPts = (*N2);

   (*Y2) = new double[totPts];
   if (outputLevel_ >= 1) printf("Kriging interpolation begins....\n");
   predict(totPts, *X2, *Y2, NULL);
   if (outputLevel_ >= 1) printf("Kriging interpolation completed.\n");
   return 0;
}

// ************************************************************************
// Generate 1D results for display
// ------------------------------------------------------------------------
int Kriging::gen1DGridData(double *X, double *Y, int ind1, double *settings, 
                           int *n, double **X2, double **Y2)
{
   int    ii, kk, totPts;
   double HX, *XX, *YY;

   if (outputLevel_ >= 1) printf("Kriging training begins....\n");
   train(X,Y);
   if (outputLevel_ >= 1) printf("Kriging training completed.\n");
  
   totPts = nPtsPerDim_;
   HX = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 

   (*X2) = new double[totPts];
   XX = new double[totPts*nInputs_];
   for (ii = 0; ii < nPtsPerDim_; ii++) 
   {
      for (kk = 0; kk < nInputs_; kk++) 
         XX[ii*nInputs_+kk] = settings[kk]; 
      XX[ii*nInputs_+ind1] = HX * ii + lowerBounds_[ind1];
      (*X2)[ii] = HX * ii + lowerBounds_[ind1];
   }
    
   YY = new double[totPts];
   if (outputLevel_ >= 1) printf("Kriging interpolation begins....\n");
   predict(totPts, XX, YY, NULL);
   if (outputLevel_ >= 1) printf("Kriging interpolation completed.\n");
   (*n) = totPts;
   (*Y2) = YY;
   delete [] XX;
   return 0;
}

// ************************************************************************
// Generate 2D results for display
// ------------------------------------------------------------------------
int Kriging::gen2DGridData(double *X, double *Y, int ind1, int ind2, 
                           double *settings, int *n, double **X2, double **Y2)
{
   int    ii, jj, kk, totPts, index;
   double *HX, *XX, *YY;

   if (outputLevel_ >= 1) printf("Kriging training begins....\n");
   train(X, Y);
   if (outputLevel_ >= 1) printf("Kriging training completed.\n");
  
   totPts = nPtsPerDim_ * nPtsPerDim_;
   HX    = new double[2];
   HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
   HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 

   XX = new double[totPts*nInputs_];
   (*X2) = new double[2*totPts];
   for (ii = 0; ii < nPtsPerDim_; ii++) 
   {
      for (jj = 0; jj < nPtsPerDim_; jj++) 
      {
         index = ii * nPtsPerDim_ + jj;
         for (kk = 0; kk < nInputs_; kk++) 
            XX[index*nInputs_+kk] = settings[kk]; 
         XX[index*nInputs_+ind1]  = HX[0] * ii + lowerBounds_[ind1];
         XX[index*nInputs_+ind2]  = HX[1] * jj + lowerBounds_[ind2];
         (*X2)[index*2]   = HX[0] * ii + lowerBounds_[ind1];
         (*X2)[index*2+1] = HX[1] * jj + lowerBounds_[ind2];
      }
   }
    
   YY = new double[totPts];
   if (outputLevel_ >= 1) printf("Kriging interpolation begins....\n");
   predict(totPts, XX, YY, NULL);
   if (outputLevel_ >= 1) printf("Kriging interpolation completed.\n");
   (*n) = totPts;
   (*Y2) = YY;
   delete [] XX;
   delete [] HX;
   return 0;
}

// ************************************************************************
// Generate 3D results for display
// ------------------------------------------------------------------------
int Kriging::gen3DGridData(double *X, double *Y, int ind1, int ind2, 
                           int ind3, double *settings, int *n, double **X2, 
                           double **Y2)
{
   int    ii, jj, ll, kk, totPts, index;
   double *HX, *XX, *YY;

   if (outputLevel_ >= 1) printf("Kriging training begins....\n");
   train(X, Y);
   if (outputLevel_ >= 1) printf("Kriging training completed.\n");
  
   totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
   HX    = new double[3];
   HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
   HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 
   HX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_ - 1); 

   XX = new double[totPts*nInputs_];
   (*X2) = new double[3*totPts];
   for (ii = 0; ii < nPtsPerDim_; ii++) 
   {
      for (jj = 0; jj < nPtsPerDim_; jj++) 
      {
         for (ll = 0; ll < nPtsPerDim_; ll++) 
         {
            index = ii * nPtsPerDim_ * nPtsPerDim_ + jj * nPtsPerDim_ + ll;
            for (kk = 0; kk < nInputs_; kk++) 
               XX[index*nInputs_+kk] = settings[kk]; 
            XX[index*nInputs_+ind1]  = HX[0] * ii + lowerBounds_[ind1];
            XX[index*nInputs_+ind2]  = HX[1] * jj + lowerBounds_[ind2];
            XX[index*nInputs_+ind3]  = HX[2] * ll + lowerBounds_[ind3];
            (*X2)[index*3]   = HX[0] * ii + lowerBounds_[ind1];
            (*X2)[index*3+1] = HX[1] * jj + lowerBounds_[ind2];
            (*X2)[index*3+2] = HX[2] * ll + lowerBounds_[ind3];
         }
      }
   }
    
   YY = new double[totPts];
   if (outputLevel_ >= 1) printf("Kriging interpolation begins....\n");
   predict(totPts, XX, YY, NULL);
   if (outputLevel_ >= 1) printf("Kriging interpolation completed.\n");
   (*n) = totPts;
   (*Y2) = YY;
   delete [] XX;
   delete [] HX;
   return 0;
}

// ************************************************************************
// Generate 4D results for display
// ------------------------------------------------------------------------
int Kriging::gen4DGridData(double *X, double *Y, int ind1, int ind2, 
                           int ind3, int ind4, double *settings, int *n, 
                           double **X2, double **Y2)
{
   int    ii, jj, ll, mm, kk, totPts, index;
   double *HX, *XX, *YY;

   if (outputLevel_ >= 1) printf("Kriging training begins....\n");
   train(X, Y);
   if (outputLevel_ >= 1) printf("Kriging training completed.\n");
  
   totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
   HX    = new double[4];
   HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
   HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 
   HX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_ - 1); 
   HX[3] = (upperBounds_[ind4] - lowerBounds_[ind4]) / (nPtsPerDim_ - 1); 

   XX = new double[totPts*nInputs_];
   (*X2) = new double[4*totPts];
   for (ii = 0; ii < nPtsPerDim_; ii++) 
   {
      for (jj = 0; jj < nPtsPerDim_; jj++) 
      {
         for (ll = 0; ll < nPtsPerDim_; ll++) 
         {
            for (mm = 0; mm < nPtsPerDim_; mm++) 
            {
               index = ii*nPtsPerDim_*nPtsPerDim_*nPtsPerDim_ +
                       jj*nPtsPerDim_*nPtsPerDim_ + ll*nPtsPerDim_ + mm;
               for (kk = 0; kk < nInputs_; kk++) 
                  XX[index*nInputs_+kk] = settings[kk]; 
               XX[index*nInputs_+ind1]  = HX[0] * ii + lowerBounds_[ind1];
               XX[index*nInputs_+ind2]  = HX[1] * jj + lowerBounds_[ind2];
               XX[index*nInputs_+ind3]  = HX[2] * ll + lowerBounds_[ind3];
               XX[index*nInputs_+ind4]  = HX[3] * mm + lowerBounds_[ind4];
               (*X2)[index*4]   = HX[0] * ii + lowerBounds_[ind1];
               (*X2)[index*4+1] = HX[1] * jj + lowerBounds_[ind2];
               (*X2)[index*4+2] = HX[2] * ll + lowerBounds_[ind3];
               (*X2)[index*4+3] = HX[3] * mm + lowerBounds_[ind4];
            }
         }
      }
   }
    
   YY = new double[totPts];
   if (outputLevel_ >= 1) printf("Kriging interpolation begins....\n");
   predict(totPts, XX, YY, NULL);
   if (outputLevel_ >= 1) printf("Kriging interpolation completed.\n");
   (*n) = totPts;
   (*Y2) = YY;
   delete [] XX;
   delete [] HX;
   return 0;
}

// ************************************************************************
// Evaluate a given point 
// ------------------------------------------------------------------------
double Kriging::evaluatePoint(double *X)
{
   int    iOne=1;
   double Y=0.0;
   predict(iOne, X, &Y, NULL);
   return Y;
}

// ************************************************************************
// Evaluate a number of points 
// ------------------------------------------------------------------------
double Kriging::evaluatePoint(int npts, double *X, double *Y)
{
   predict(npts, X, Y, NULL);
   return 0.0;
}

// ************************************************************************
// Evaluate a given point, return also the standard deviation 
// ------------------------------------------------------------------------
double Kriging::evaluatePointFuzzy(double *X, double &Ystd)
{
   int    iOne=1;
   double Y=0.0;
   predict(iOne, X, &Y, &Ystd);
   return Y;
}

// ************************************************************************
// Evaluate a number of points, return also the standard deviation 
// ------------------------------------------------------------------------
double Kriging::evaluatePointFuzzy(int npts, double *X, double *Y, 
                                   double *Ystds)
{
   predict(npts, X, Y, Ystds);
   return 0.0;
}

// ************************************************************************
// training 
// ------------------------------------------------------------------------
double Kriging::train(double *X, double *Y)
{
   int    ii, jj, kk, nBasis, count, status, nDists, Cleng;
   double *XDists, dist, ddata; 
   char   pString[500];

   // clean up 
   if (XNormalized_ != NULL) delete [] XNormalized_;
   if (YNormalized_ != NULL) delete [] YNormalized_;
   if (Cmatrix_ != NULL) delete [] Cmatrix_;
   XNormalized_ = NULL;
   YNormalized_ = NULL;
   Cmatrix_ = NULL;

   // normalize the input and outputs 
   XNormalized_ = new double[nSamples_*nInputs_];
   initInputScaling(X, XNormalized_, 1);
   YNormalized_ = new double[nSamples_];
   initOutputScaling(Y, YNormalized_);
   KRI_YStd = YStd_;

   // compute distances between all pairs of inputs (nDists, XDists)
   computeDistances(&XDists, &nDists);

   // fast mode: compute correlation ==> Cmatrix_
   if (fastMode_ == 1)
   {
      printf("Kriging training (1) begins.... (order = %d)\n",pOrder_);
      if (nSamples_ <= nInputs_ + 1) pOrder_ = 0;
      if (pOrder_ == 1) nBasis = nInputs_ + 1;
      else              nBasis = 1;

      Cleng = nSamples_ + nBasis;
      Cmatrix_ = new double[Cleng*Cleng];
      count = 0;
      for (jj = 0; jj < nSamples_; jj++)
      {
         Cmatrix_[jj*Cleng+jj] = 1.0 + 1e-15 * nSamples_;
         if (dataStdDevs_ != NULL) 
            Cmatrix_[jj*Cleng+jj] += pow(dataStdDevs_[jj]/YStd_, 2.0);
         for (kk = jj+1; kk < nSamples_; kk++)
         {
            dist = 0.0;
            for (ii = 0; ii < nInputs_; ii++)
               dist += pow(XDists[count*nInputs_+ii]/Thetas_[ii], 2.0);
            Cmatrix_[jj*Cleng+kk] = Cmatrix_[kk*Cleng+jj] = exp(-dist);
            count++;
         }
      }
      for (jj = 0; jj < nSamples_; jj++)
      {
         Cmatrix_[jj*Cleng+nSamples_] = 1.0;
         Cmatrix_[nSamples_*Cleng+jj] = 1.0; 
      }
      Cmatrix_[nSamples_*Cleng+nSamples_] = 1.0; 
      if (pOrder_ == 1) 
      {
         for (jj = 0; jj < nSamples_; jj++)
         {
            for (kk = nSamples_+1; kk < nSamples_+nInputs_+1; kk++)
            {
               ddata = XNormalized_[jj*nInputs_+kk-nSamples_-1];
               Cmatrix_[kk*Cleng+jj] = ddata;
               Cmatrix_[jj*Cleng+kk] = ddata;
            }
         }
      }
      for (jj = nSamples_; jj < Cleng; jj++)
         for (kk = nSamples_; kk < Cleng; kk++)
            Cmatrix_[jj*Cleng+kk] = 0.0;

      LUPivots_ = new int[Cleng];
      dgetrf_(&Cleng, &Cleng, Cmatrix_, &Cleng, LUPivots_, &status);

      if (status != 0) 
      {
         printf("Kriging ERROR: LU decomposition not successful.\n");
         exit(1);
      }
      printf("Kriging training (1) ends.\n");
   }
   else
   // slow mode: optimize
   {
      int    maxfun, pLevel, nPts, iOne=1, iZero=0, nSamOpt;
      int    *samStates, mode;
      double *TValues, *TUppers, *TLowers, rhobeg=1.0, rhoend=1.0e-4;
      double *work, *samInputs, *samOutputs, *optThetas, *MMatrix=NULL;
      double optY=PSUADE_UNDEFINED, *SMatrix=NULL, *FMatrix=NULL;
      double *FMatTmp;
      FILE   *fp=NULL;
      Sampling *sampler;

      fp = fopen("ps_terminate", "r");
      if (fp != NULL)
      {
         printf("Kriging ERROR: remove the 'ps_terminate' file\n");
         printf("               first and re-do.\n");
         fclose(fp);
         exit(1);
      }

      if (nSamples_ <= nInputs_ + 1) pOrder_ = 0;
      if (pOrder_ == 1) nBasis = nInputs_ + 1;
      else              nBasis = 1;
      Cleng = nSamples_ + nBasis;

      TUppers = new double[nInputs_+1];
      TLowers = new double[nInputs_+1];
      for (ii = 0; ii < nInputs_; ii++)
      {
         if (XMeans_[ii] == 0 && XStds_[ii] == 1.0)
         {
            TUppers[ii] = 30.0 * (upperBounds_[ii]-lowerBounds_[ii]);;
            TLowers[ii] = 0.01 * (upperBounds_[ii]-lowerBounds_[ii]);;
         }
         else
         {
            TUppers[ii] = 30.0*2;
            TLowers[ii] = 0.01*2;
         }
         if (psGMMode_ == 1)
         {
            sprintf(pString,
                    "Kriging: Enter optimization lower bound for input %d : ",ii+1);
            TLowers[ii] = getDouble(pString);
            if (TLowers[ii] <= 0.0)
            {
               printf("Kriging ERROR: lower bound <= 0\n");
               exit(1);
            }
            sprintf(pString,
                    "Kriging: Enter optimization upper bound for input %d : ",ii+1);
            TUppers[ii] = getDouble(pString);
            if (TLowers[ii] > TUppers[ii])
            {
               printf("Kriging ERROR: lower bound >= upper bound\n");
               exit(1);
            }
         }
      }
      rhobeg = TUppers[0] - TLowers[0];
      for (ii = 1; ii < nInputs_; ii++)
      {
         ddata = TUppers[ii] - TLowers[ii];
         if (ddata < rhobeg) rhobeg = ddata;
      }
      rhobeg *= 0.5;
      rhoend = rhobeg * optTolerance_;
      TValues = new double[nInputs_+1];
      nPts = (nInputs_ + 1) * (nInputs_ + 2) / 2;
      work = new double[(nPts+5)*(nPts+nInputs_)+3*nInputs_*(nInputs_+5)/2+1];
      printEquals(0);
      printf("* Kriging optimization tolerance = %e\n", rhoend);

      if (fastMode_ == 2)
      {
         printf("Kriging training (2) begins.... (order = %d)\n",pOrder_);
         nSamOpt = 1;
         samInputs  = new double[nSamOpt * nInputs_];
         for (ii = 0; ii < nInputs_; ii++) samInputs[ii] = Thetas_[ii];
      }
      else
      {
         printf("Kriging training (3) begins.... (order = %d)\n",pOrder_);
         mode = 1;
         nSamOpt = 10;
         if (psGMMode_ == 1)
         {
            printf("Kriging: slow mode with multi-start optimization.\n");
            printf("Choose sampling method to generate multi-start.\n");
            sprintf(pString, "Sampling method (1-QMC, 2-LHS, 3-FF) : ");
            mode = getInt(1,3,pString);
            if (mode != 3)
            {
               sprintf(pString, "Enter sample size (1 - 100) : ");
               nSamOpt = getInt(1,100,pString);
            } 
         }
         if (mode == 1)
            sampler = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_LPTAU);
         else if (mode == 2)
            sampler = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_LHS);
         else
            sampler = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_FF4);
 
         sampler->setPrintLevel(0);
         sampler->setInputBounds(nInputs_, TLowers, TUppers);
         sampler->setOutputParams(iOne);
         sampler->setSamplingParams(nSamOpt, iOne, iZero);
         sampler->initialize(0);
         nSamOpt = sampler->getNumSamples();
         samInputs  = new double[nSamOpt * nInputs_];
         samOutputs = new double[nSamOpt];
         samStates  = new int[nSamOpt];
         sampler->getSamples(nSamOpt, nInputs_, iOne, samInputs,
                             samOutputs, samStates);
         delete [] samOutputs;
         delete [] samStates;
         delete sampler;
      }

      KRI_XDists = XDists;
      KRI_nSamples = nSamples_;
      KRI_nInputs = nInputs_;
      KRI_pOrder = pOrder_;
      KRI_X = XNormalized_;
      KRI_Y = YNormalized_;
      KRI_Ytmp = new double[Cleng];
      KRI_OptThetas = Thetas_;
      maxfun = 10000;
      KRI_OptY = PSUADE_UNDEFINED;
      KRI_outputLevel = outputLevel_;
      optThetas = new double[nInputs_];
      SMatrix = new double[nSamples_*nSamples_];
      FMatrix = new double[nSamples_*nBasis];
      FMatTmp = new double[nSamples_*nBasis];
      MMatrix = new double[nBasis*nBasis];
      KRI_SMatrix = SMatrix;
      KRI_FMatrix = FMatrix;
      KRI_FMatTmp = FMatTmp;
      KRI_MMatrix = MMatrix;

      for (kk = 0; kk < nSamOpt; kk++)
      {
         if (outputLevel_ >= 1) 
            printf("Kriging multi-start optimization: start = %d (%d)\n",
                   kk+1, nSamOpt);
         KRI_iter = 0;
         pLevel = 8888;
         for (ii = 0; ii < nInputs_; ii++) 
            TValues[ii] = samInputs[kk*nInputs_+ii];
         if (outputLevel_ >= 1) 
         {
            for (ii = 0; ii < nInputs_; ii++) 
               printf("Kriging: Input %4d initial length scale = %e\n",
                      ii+1,TValues[ii]);
         }
#ifdef HAVE_BOBYQA
         bobyqa_(&nInputs_,&nPts,TValues,TLowers,TUppers,&rhobeg,&rhoend,
                 &pLevel, &maxfun, work);
         if (outputLevel_ >= 1) 
         {
            for (ii = 0; ii < nInputs_; ii++) 
               printf("Kriging: Input %4d final length scale = %e\n",
                      ii+1,TValues[ii]);
            printf("Kriging final objective value = %e\n", KRI_OptY);
         }
         if (KRI_OptY < optY)
         {
            optY = KRI_OptY;
            for (ii = 0; ii < nInputs_; ii++)
               optThetas[ii] = KRI_OptThetas[ii];
         }
         fp = fopen("ps_terminate", "r");
         if (fp != NULL)
         {
            fclose(fp);
            printf("Kriging: terminating ....\n");
            break;
         }
         fp = fopen("ps_print", "r");
         if (fp != NULL)
         {
            fclose(fp);
            printf("Kriging: output level set to 2.\n");
            outputLevel_ = 2;
            KRI_outputLevel = outputLevel_;
         }
         fp = fopen("ps_rs_expert", "r");
         if (fp != NULL)
         {
            fclose(fp);
            printf("Kriging: turn on rs_expert mode.\n");
            psRSExpertMode_ = 1;
         }
#else
         printf("ERROR: Bobyqa optimizer not installed.\n");
         exit(1);
#endif
      }
      for (ii = 0; ii < nInputs_; ii++) Thetas_[ii] = optThetas[ii];
      if (outputLevel_ >= 1) 
      {
         for (ii = 0; ii < nInputs_; ii++) 
            printf("Kriging: Input %4d optimal length scale = %e\n",
                   ii+1,Thetas_[ii]);
      }
      delete [] work;
      delete [] TValues;
      delete [] TUppers;
      delete [] TLowers;
      delete [] KRI_Ytmp;
      delete [] samInputs;
      delete [] optThetas;
      KRI_OptThetas = NULL;
      KRI_XDists = NULL;
      KRI_X = NULL;
      KRI_Y = NULL;
      KRI_Ytmp = NULL;
      if (SMatrix != NULL) delete [] SMatrix;
      if (FMatrix != NULL) delete [] FMatrix;
      if (FMatTmp != NULL) delete [] FMatTmp;
      if (MMatrix != NULL) delete [] MMatrix;
      KRI_SMatrix = NULL;
      KRI_FMatrix = NULL;
      KRI_FMatTmp = NULL;
      KRI_MMatrix = NULL;

      count = 0;
      Cmatrix_ = new double[Cleng*Cleng];
      for (jj = 0; jj < nSamples_; jj++)
      {
         Cmatrix_[jj*Cleng+jj] = 1.0 + 1e-15 * nSamples_;
         if (dataStdDevs_ != NULL) 
            Cmatrix_[jj*Cleng+jj] += pow(dataStdDevs_[jj]/YStd_, 2.0);
         for (kk = jj+1; kk < nSamples_; kk++)
         {
            dist = 0.0;
            for (ii = 0; ii < nInputs_; ii++)
               dist += pow(XDists[count*nInputs_+ii]/Thetas_[ii], 2.0);
            dist = exp(-dist);
            if (dist < 1.0e-16) dist = 0.0;
            Cmatrix_[jj*Cleng+kk] = dist;
            Cmatrix_[kk*Cleng+jj] = dist;
            count++;
         }
      }
      for (jj = 0; jj < nSamples_; jj++)
      {
         Cmatrix_[jj*Cleng+nSamples_] = 1.0;
         Cmatrix_[nSamples_*Cleng+jj] = 1.0; 
      }
      Cmatrix_[nSamples_*Cleng+nSamples_] = 0.0;
      if (pOrder_ == 1) 
      {
         for (jj = 0; jj < nSamples_; jj++)
         {
            for (kk = nSamples_+1; kk < nSamples_+nInputs_+1; kk++)
            {
               ddata = XNormalized_[jj*nInputs_+kk-nSamples_-1];
               Cmatrix_[kk*Cleng+jj] = ddata;
               Cmatrix_[jj*Cleng+kk] = ddata;
            }
         }
      }
      for (jj = nSamples_; jj < Cleng; jj++)
         for (kk = nSamples_; kk < Cleng; kk++)
            Cmatrix_[jj*Cleng+kk] = 0.0;

      LUPivots_ = new int[Cleng];
      dgetrf_(&Cleng, &Cleng, Cmatrix_, &Cleng, LUPivots_, &status);

      if (status != 0) 
      {
         printf("Kriging ERROR: LU decomposition not successful.\n");
         exit(1);
      }
      if (fastMode_ == 2)
           printf("Kriging training (2) ends.\n");
      else printf("Kriging training (3) ends.\n");
   }
   delete [] XDists;
   return 0.0;
}

// ************************************************************************
// predict 
// ------------------------------------------------------------------------
double Kriging::predict(int length, double *X, double *Y, double *YStds)
{
   int    ii, jj, kk, status, nRows;
   double ddata, dist, mean, stdev;
   char   trans='N';

   nRows = nSamples_ + 1;
   if (pOrder_ == 1) nRows = nRows + nInputs_;
   if (workArray_ == NULL)
   {
      workArray_ = new double[2*nRows*length];
      workX_ = new double[length*nInputs_];
      workLength_ = length;
   }
   else if (workArray_ != NULL && length > workLength_)
   {
      delete [] workArray_;
      delete [] workX_;
      workLength_ = length;
      workArray_ = new double[2*nRows*length];
      workX_ = new double[length*nInputs_];
   }
   for (ii = 0; ii < nInputs_; ii++)
   {
      mean  = XMeans_[ii];
      stdev = XStds_[ii];
      for (jj = 0; jj < length; jj++)
         workX_[jj*nInputs_+ii] = (X[jj*nInputs_+ii] - mean) / stdev;
   }

#if 1
   for (kk = 0; kk < length; kk++)
   {
      for (jj = 0; jj < nSamples_; jj++)
      {
         dist = 0.0;
         for (ii = 0; ii < nInputs_; ii++)
         {
            ddata = XNormalized_[jj*nInputs_+ii] - workX_[kk*nInputs_+ii];
            dist += ddata * ddata / (Thetas_[ii] * Thetas_[ii]);
         }
         ddata = exp(-dist);
         if (ddata < 1.0e-16) ddata = 0.0;
         workArray_[kk*nRows+jj] = ddata;
         workArray_[nRows*length+kk*nRows+jj] = ddata;
      }
      workArray_[kk*nRows+nSamples_] = 1.0;
      workArray_[nRows*length+kk*nRows+nSamples_] = 1.0;
      if (pOrder_ == 1)
      {
         for (jj = 0; jj < nInputs_; jj++)
            workArray_[kk*nRows+nSamples_+jj+1] = 
               workArray_[nRows*length+kk*nRows+nSamples_+jj+1] = 
                                                  workX_[kk*nInputs_+jj];
      }
   }
   dgetrs_(&trans,&nRows,&length,Cmatrix_,&nRows,LUPivots_,workArray_,
           &nRows,&status);
   for (kk = 0; kk < length; kk++)
   {
      ddata = 0.0;
      for (jj = 0; jj < nSamples_; jj++)
         ddata += YNormalized_[jj] * workArray_[kk*nRows+jj];
      Y[kk] = ddata * YStd_ + YMean_;
      if (YStds != NULL)
      {
         ddata = 0.0;
         for (jj = 0; jj < nRows; jj++)
            ddata += workArray_[kk*nRows+jj] * 
                     workArray_[nRows*length+kk*nRows+jj];
         ddata = 1.0 - ddata;
         ddata = ddata * YStd_ * YStd_;
         if (ddata < 0.0)
         {
            printf("Kriging WARNING: prediction variance < 0\n");
            ddata = 0.0;
         }
         YStds[kk] = sqrt(ddata);
      }
   }
#else
   for (kk = 0; kk < length; kk++)
   {
      for (jj = 0; jj < nSamples_; jj++)
      {
         dist = 0.0;
         for (ii = 0; ii < nInputs_; ii++)
         {
            ddata = XNormalized_[jj*nInputs_+ii] - workX_[kk*nInputs_+ii];
            dist += pow(ddata/Thetas_[ii], 2.0);
         }
         ddata = exp(-dist);
         workArray_[jj] = workArray_[jj+nRows] = ddata;
      }
      workArray_[nSamples_] = workArray_[nSamples_+nRows] = 1.0;
      if (pOrder_ == 1)
      {
         for (jj = 0; jj < nInputs_; jj++)
         {
            workArray_[nSamples_+jj+1] = workX_[kk*nInputs_+jj];
            workArray_[nRows+nSamples_+jj+1] = workX_[kk*nInputs_+jj];
         }
      }
      dgetrs_(&trans,&nRows,&iOne,Cmatrix_,&nRows,LUPivots_,workArray_,
              &nRows,&status);
      ddata = 0.0;
      for (jj = 0; jj < nSamples_; jj++)
         ddata += YNormalized_[jj] * workArray_[jj];
      Y[kk] = ddata * YStd_ + YMean_;
      if (YStds != NULL)
      {
         ddata = 0.0;
         for (jj = 0; jj < nRows; jj++)
            ddata += workArray_[jj] * workArray_[nRows+jj];
         ddata = 1.0 - ddata;
         ddata = ddata * YStd_ * YStd_;
         if (ddata < 0.0)
         {
            printf("Kriging WARNING: prediction variance < 0\n");
            ddata = 0.0;
         }
         YStds[kk] = sqrt(ddata);
      }
   }
#endif
   return 0.0;
}

// ************************************************************************
// compute distances between all pairs of inputs
// ------------------------------------------------------------------------
int Kriging::computeDistances(double **XDists, int *length)
{
   int    ii, jj, kk, count;
   double *LDists, dist;

   LDists = new double[(nSamples_*(nSamples_-1)/2)*nInputs_];
   if (LDists == NULL) 
   {
      printf("Kriging ERROR: allocation problem.\n");
      exit(1);
   }
   count = 0;
   for (jj = 0; jj < nSamples_; jj++)
   {
      for (kk = jj+1; kk < nSamples_; kk++)
      {
         dist = 0.0;
         for (ii = 0; ii < nInputs_; ii++)
         {
            LDists[count*nInputs_+ii] = XNormalized_[jj*nInputs_+ii] - 
                                XNormalized_[kk*nInputs_+ii];
            if (LDists[count*nInputs_+ii] < 0) 
               LDists[count*nInputs_+ii] = - LDists[count*nInputs_+ii]; 
            dist += pow(LDists[count*nInputs_+ii], 2.0);
         }
         if (dist == 0.0)
         {
            printf("Kriging ERROR: repeated sample points.\n");
            printf("               Prune repeated points and re-run.\n");
            printf("Sample %d : (with sample point %d)\n", kk+1, jj+1);
            for (ii = 0; ii < nInputs_; ii++)
               printf("   Input %d : %e\n",ii+1,
                      XNormalized_[kk*nInputs_+ii]*XStds_[ii] +XMeans_[ii]);
            exit(1);
         }
         count++;
      }
   }
   (*length) = count;
   (*XDists) = LDists;
   return 0;
}

// ************************************************************************
// set parameters
// ------------------------------------------------------------------------
double Kriging::setParams(int targc, char **targv)
{
   int    ii, *iArray;
   double *lengthScales, mmax, range;
   char   pString[500];
   FILE   *fp=NULL;

   if (targc > 0 && !strcmp(targv[0], "setMode1"))
      fastMode_ = 1;
   else if (targc > 0 && !strcmp(targv[0], "setMode2"))
      fastMode_ = 2;
   else if (targc > 0 && !strcmp(targv[0], "setMode3"))
      fastMode_ = 3;
   else if (targc > 0 && !strcmp(targv[0], "rank"))
   {
      lengthScales = new double[nInputs_];
      mmax = 0.0;
      for (ii = 0; ii < nInputs_; ii++)
      {
         lengthScales[ii] = 1.0 / Thetas_[ii];
         if (XMeans_[ii] == 0 && XStds_[ii] == 1)
         {
            range = upperBounds_[ii] - lowerBounds_[ii];
            lengthScales[ii] *= range;
         }
         if (lengthScales[ii] > mmax) mmax = lengthScales[ii];
      }
      for (ii = 0; ii < nInputs_; ii++)
         lengthScales[ii] = lengthScales[ii] / mmax * 100.0;

      if (psPlotTool_ == 1)
           fp = fopen("scilabkrisa.sci", "w");
      else fp = fopen("matlabkrisa.m", "w");
      if (fp == NULL)
      {
         printf("Kriging ERROR: something wrong with opening a write file.\n");
      }
      else
      {
         fprintf(fp, "n = %d;\n", nInputs_);
         fprintf(fp, "Y = [\n");
         for (ii = 0; ii < nInputs_; ii++)
            fprintf(fp, "%24.16e \n", PABS(lengthScales[ii]));
         fprintf(fp, "]; \n");
         fprintf(fp, "ymax = max(Y);\n");
         fprintf(fp, "ymin = 0;\n");
         fprintf(fp, "if (ymax == ymin)\n");
         fprintf(fp, "   ymax = ymax * 0.1;\n");
         fprintf(fp, "end;\n");
         fwritePlotCLF(fp);
         fprintf(fp, "bar(Y,0.8);\n");
         fwritePlotAxes(fp);
         sprintf(pString, "Kriging Ranking");
         fwritePlotTitle(fp, pString);
         sprintf(pString, "Input Numbers");
         fwritePlotXLabel(fp, pString);
         sprintf(pString, "Kriging Measure");
         fwritePlotYLabel(fp, pString);
         if (psPlotTool_ == 1)
         {
            fprintf(fp, "a=gca();\n");
            fprintf(fp, "a.data_bounds=[0, ymin; n+1, ymax+0.01*(ymax-ymin)];\n");
         }
         else
         {
            fprintf(fp,"axis([0 n+1 ymin ymax+0.01*(ymax-ymin)])\n");
         }
         fclose(fp);
         if (psPlotTool_ == 1)
              printf("Kriging ranking in file scilabkrisa.sci\n");
         else printf("Kriging ranking in file matlabkrisa.m\n");
      }
      iArray = new int[nInputs_];
      for (ii = 0; ii < nInputs_; ii++) iArray[ii] = ii;
      sortDbleList2a(nInputs_, lengthScales, iArray);
      printAsterisks(0);
      printf("* Kriging screening rankings\n");
      printAsterisks(0);
      for (ii = nInputs_-1; ii >= 0; ii--)
         printf("*  Rank %3d : Input = %4d (score = %5.1f) (ref = %e)\n",
                nInputs_-ii, iArray[ii]+1, lengthScales[ii], 
                0.01*lengthScales[ii]*mmax);
      printAsterisks(0);
      delete [] lengthScales;
   }
   return 0.0;
}

