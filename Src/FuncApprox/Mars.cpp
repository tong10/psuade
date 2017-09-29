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
// DATE   : 2003
// ************************************************************************
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "Mars.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "Psuade.h"
#include "PsuadeConfig.h"

#define PABS(x) ((x) > 0 ? (x) : (-(x)))

#ifdef HAVE_MARS
extern "C" 
{
  void mars_process(int,int,double**,double*,float*,int,int,int*,
                 float*,int*);
  void mars_fmod(int,int,double**,double*,float*,int*);
}
#endif

// ************************************************************************
// Constructor for object class Mars
// ------------------------------------------------------------------------
Mars::Mars(int nInputs,int nSamples) : FuncApprox(nInputs,nSamples)
{
#ifdef HAVE_MARS
   int  ss, ii;
   char pString[501], *cString, winput1[500], winput2[500];

   faID_ = PSUADE_RS_MARS;

   nBasisFcns_ = 100;
   if (nBasisFcns_ > nSamples) nBasisFcns_ = nSamples - 1;

   if (nInputs >= 8) maxVarPerBasis_ = 8;
   else              maxVarPerBasis_ = nInputs;

   varFlags_ = new int[nInputs];
   for (ii = 0; ii < nInputs; ii++) varFlags_[ii] = 1;
   for (ii = 0; ii < nSamples_; ii++) weights_[ii] = 1000.0;
   chooseWght_ = 0;

   if (psRSExpertMode_ != 1 && psConfig_ != NULL)
   {
      cString = psConfig_->getParameter("MARS_num_bases");
      if (cString != NULL)
      {
         sscanf(cString, "%s %s %d", winput1, winput2, &ii);
         if (ii < 10 || ii > nSamples)
         {
            printf("Mars INFO: nbasis from config file not valid.\n");
            printf("           nbasis kept at %d.\n", nBasisFcns_);
         }
         else 
         {
            nBasisFcns_ = ii;
            printf("Mars INFO: number of basis set to %d (config).\n",nBasisFcns_);
         }
      }
      cString = psConfig_->getParameter("MARS_interaction");
      if (cString != NULL)
      {
         sscanf(cString, "%s %s %d", winput1, winput2, &ii);
         if (ii > nInputs || ii < 1)
         {
            printf("Mars INFO: interaction from config file not valid.\n");
            printf("           interaction kept at %d.\n", maxVarPerBasis_);
         }
         else 
         {
            maxVarPerBasis_ = ii;
            printf("Mars INFO: interaction set to %d (config).\n", maxVarPerBasis_);
         }
      }
   }

   if (psRSExpertMode_ == 1 && psInteractive_ == 1)
   {
      printf("Mars: Current number of basis functions = %d\n", nBasisFcns_);
      if (nSamples > 10)
      {
         sprintf(pString,"Enter the number of basis functions (>=10, <= %d): ",
                 nSamples);
         nBasisFcns_ = getInt(10, nSamples, pString);
      }
      else
      {
         sprintf(pString,"Enter the number of basis functions (>=%d, <= %d): ",
                 nSamples, nSamples);
         nBasisFcns_ = getInt(nSamples, nSamples, pString);
      }
      printf("Mars: Current degree of interactions    = %d\n",maxVarPerBasis_);
      sprintf(pString, "Enter the degree of interactions (<=%d) : ", nInputs);
      maxVarPerBasis_ = getInt(1, nInputs, pString);
   }

   wgts_ = new float[nSamples];
   for (ss = 0; ss < nSamples; ss++) wgts_[ss] = 1.0;
   fm_ = NULL;
   im_ = NULL;
   noGen_ = 0;
#else
   printf("PSUADE ERROR : Mars not installed.\n");
   exit(1);
#endif
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
Mars::~Mars()
{
   if (wgts_     != NULL) delete [] wgts_;
   if (varFlags_ != NULL) delete [] varFlags_;
   if (fm_       != NULL) delete [] fm_;
   if (im_       != NULL) delete [] im_;
}

// ************************************************************************
// initialize
// ------------------------------------------------------------------------
int Mars::initialize(double *XIn, double *YIn)
{
#ifdef HAVE_MARS
   int    length, ss, ii, nsm, nin, nfm, nim, errflag, count;
   double **X, *XX, *Y, dmean, dstd;
   char   *cStr, pString[500], response[500], line[1001], *cString;
   char   word[501], equal[101];
   FILE   *fp=NULL;

   response[0] = 'n';
   if (psRSExpertMode_ != 1 && psConfig_ != NULL)
   {
      cString = psConfig_->getParameter("normalize_outputs");
      if (cString != NULL) 
      {
         printf("Mars INFO: output to be normalized.\n");
         response[0] = 'y';
      }
   }
   if (psRSExpertMode_ == 1 && psInteractive_ == 1)
   {
      sprintf(pString, "Mars: normalize output? (y or n) ");
      getString(pString, response);
   }

   XX = new double[nSamples_*nInputs_];
   initInputScaling(XIn, XX, 1);
   Y = new double[nSamples_];
   if (response[0] == 'y')
   {
      initOutputScaling(YIn, Y);
   }
   else
   {
      for (ii = 0; ii < nSamples_; ii++) Y[ii] = YIn[ii];
      YMean_ = 0.0;
      YStd_ = 1.0;
   }
   
   if (fm_ != NULL) delete [] fm_;
   if (im_ != NULL) delete [] im_;
   length = 3 + nBasisFcns_ * (5 * maxVarPerBasis_ + nSamples_+6) +
            2 * nInputs_ + nSamples_;
   fm_    = new float[length];
   for (ii = 0; ii < length; ii++) fm_[ii] = PSUADE_UNDEFINED;
   length = 21 + nBasisFcns_ * (3 * maxVarPerBasis_ + 8);
   im_    = new int[length];
   for (ii = 0; ii < length; ii++) im_[ii] = -9999;
   if (chooseWght_ == 1)
   {
      for (ss = 0; ss < nSamples_; ss++)
      {
         if (PABS(Y[ss]) < 1.0e-12) wgts_[ss] = 1.0;
         else                       wgts_[ss] = (float) PABS(Y[ss]);
      }
   }
   else
   {
      for (ss = 0; ss < nSamples_; ss++) wgts_[ss] = (float) weights_[ss];
   }

   X = new double*[nSamples_];
   for (ss = 0; ss < nSamples_; ss++) X[ss] = new double[nInputs_];
   for (ss = 0; ss < nSamples_; ss++) 
      for (ii = 0; ii < nInputs_; ii++) X[ss][ii] = XX[ss*nInputs_+ii];
   delete [] XX;

   if (outputLevel_ >= 2) 
      printf("Mars: nBasis, maxVarPerBasis = %d %d\n", 
             nBasisFcns_, maxVarPerBasis_);
   if (outputLevel_ >= 2 || psMasterMode_ == 1)
   {
      printf("Entering Mars (process)\n");
      printf("If it crashes here, it is mars_process problem.\n");
      printf("One way to solve the problem is to use different nSamples.\n");
   }
   mars_process(nSamples_, nInputs_, X, Y, wgts_, nBasisFcns_,
                maxVarPerBasis_, varFlags_, fm_, im_);
   if (outputLevel_ >= 2 || psMasterMode_ == 1) 
      printf("Returned from Mars (process).\n");
   for (ss = 0; ss < nSamples_; ss++) delete [] X[ss];
   delete [] X;

   if (noGen_ == 1 || psRSCodeGen_ == 0) return 0;
   fp = fopen("psuade_rs.info", "w");
   if (fp == NULL)
   {
      printf("ERROR: Cannot open file psuade_rs.info.\n");
      return 0;
   }
   fprintf(fp,"/* This file contains information to re-construct MARS\n");
   fprintf(fp,"   response surface offline. Follow the steps below:\n");
   fprintf(fp,"   1. Rename this file to, say, main.c\n");
   fprintf(fp,"   2. Compile main.c (cc -o main main.c -lm) and run \n");
   fprintf(fp,"   3. Run: main input output (input file has number of \n");
   fprintf(fp,"      inputs followed by the input values. Upon termination,\n");
   fprintf(fp,"      the result will be stored in 'output') */\n");
   fprintf(fp,"/* *************************************/\n");
   fprintf(fp,"/* MARS interpolator from PSUADE.      */\n");
   fprintf(fp,"/* ====================================*/\n");
   fprintf(fp,"#include <math.h>\n");
   fprintf(fp,"#include <stdlib.h>\n");
   fprintf(fp,"#include <stdio.h>\n");
   fprintf(fp,"int interpolate(int,double*,double*);\n");
   fprintf(fp,"main(int argc, char **argv) {\n");
   fprintf(fp,"  int    i, iOne=1, nInps;\n");
   fprintf(fp,"  double X[%d], Y, Std;\n",nInputs_);
   fprintf(fp,"  FILE   *fIn=NULL, *fOut=NULL;\n");
   fprintf(fp,"  if (argc != 3) {\n");
   fprintf(fp,"     printf(\"ERROR: not enough argument.\\n\");\n");
   fprintf(fp,"     exit(1);\n");
   fprintf(fp,"  }\n");
   fprintf(fp,"  fIn = fopen(argv[1], \"r\");\n");
   fprintf(fp,"  if (fIn == NULL) {\n");
   fprintf(fp,"     printf(\"ERROR: cannot open input file.\\n\");\n");
   fprintf(fp,"     exit(1);\n");
   fprintf(fp,"  }\n");
   fprintf(fp,"  fscanf(fIn, \"%%d\", &nInps);\n");
   fprintf(fp,"  if (nInps != %d) {\n", nInputs_);
   fprintf(fp,"    printf(\"ERROR - wrong nInputs.\\n\");\n");
   fprintf(fp,"    exit(1);\n");
   fprintf(fp,"  }\n");
   fprintf(fp,"  for (i=0; i<%d; i++) fscanf(fIn, \"%%lg\", &X[i]);\n",nInputs_);
   fprintf(fp,"  fclose(fIn);\n");
   fprintf(fp,"  interpolate(iOne,X,&Y);\n");
   fprintf(fp,"  printf(\"Y = %%e\\n\", Y);\n");
   fprintf(fp,"  fOut = fopen(argv[2], \"w\");\n");
   fprintf(fp,"  if (fOut == NULL) {\n");
   fprintf(fp,"     printf(\"ERROR: cannot open output file.\\n\");\n");
   fprintf(fp,"     exit(1);\n");
   fprintf(fp,"  }\n");
   fprintf(fp,"  fprintf(fOut,\" %%e\\n\", Y);\n");
   fprintf(fp,"  fclose(fOut);\n");
   fprintf(fp,"}\n\n");
   fprintf(fp,"/* *************************************/\n");
   fprintf(fp,"/* **** MARS interpolation function ****/\n");
   fprintf(fp,"/* X[0], X[1],   .. X[m-1]   - first point\n");
   fprintf(fp," * X[m], X[m+1], .. X[2*m-1] - second point\n");
   fprintf(fp," * ... */\n");
   fprintf(fp,"/* ====================================*/\n");
   fprintf(fp,"double Ymean = %24.16e;\n",YMean_);
   fprintf(fp,"double Ystd  = %24.16e;\n",YStd_);
   fprintf(fp,"int icat(double, int, double *);\n");
   fprintf(fp,"int    *getIm();\n");
   fprintf(fp,"double *getFm();\n");
   fprintf(fp,"double *getXM();\n");
   fprintf(fp,"double *getXS();\n");
   fprintf(fp,"int interpolate(int npts, double *X,double *Y){\n");
   fprintf(fp,"  int    k, nk, ss, nn, ip, ind;\n"); 
   fprintf(fp,"  int    nInps=%d, *im=NULL;\n",nInputs_);
   fprintf(fp,"  double Yt, az, *tb, *cm, phi, t, u, v, *XX;\n");
   fprintf(fp,"  double *fm=NULL, *XM=NULL, *XS=NULL, dt;\n");
   fprintf(fp,"  XX = (double *) malloc(nInps*sizeof(double));\n");
   fprintf(fp,"  fm = getFm();\n");
   fprintf(fp,"  im = getIm();\n");
   fprintf(fp,"  XM = getXM();\n");
   fprintf(fp,"  XS = getXS();\n");
   fprintf(fp,"  nk  = im[4];\n");
   fprintf(fp,"  ind = im[10] - 1; az  = fm[ind];\n");
   fprintf(fp,"  ind = im[11] - 1; tb  = &fm[ind];\n"); 
   fprintf(fp,"  ind = im[14] - 1; cm  = &fm[ind];\n");
   fprintf(fp,"  for (ss = 0; ss < npts; ss++) {\n");
   fprintf(fp,"    Yt = az;\n"); 
   fprintf(fp,"    for (nn = 0; nn < nk; nn++) {\n");
   fprintf(fp,"      if (tb[nn*5] != 0.0) {\n");
   fprintf(fp,"        phi = 1.0;\n");
   fprintf(fp,"        ip = nn;\n");
   fprintf(fp,"        while (ip > -1) {\n");
   fprintf(fp,"          t = tb[ip*5+1];\n");
   fprintf(fp,"          v = t;\n");
   fprintf(fp,"          if (v < 0) v = - v;\n");
   fprintf(fp,"          ind = floor(v+0.1) - 1;\n"); 
   fprintf(fp,"          if (cm[2*ind] <= 0.0) {\n");
   fprintf(fp,"            v = 1.0;\n");
   fprintf(fp,"            if (t < 0.0) v = -1.0;\n");
   fprintf(fp,"            dt = (X[ind+nInps*ss] - XM[ind])/XS[ind];\n");
   fprintf(fp,"            u = v * (dt-tb[ip*5+2]); \n");
   fprintf(fp,"            if (u < 0.0) u = 0.0;\n");
   fprintf(fp,"          }\n");
   fprintf(fp,"          else {\n");
   fprintf(fp,"            dt = (X[ind+nInps*ss] - XM[ind])/XS[ind];\n");
   fprintf(fp,"            k = icat(dt, ind, cm);\n");
   fprintf(fp,"            if (k != 0) {\n");
   fprintf(fp,"              ind = floor(tb[ip*5+2]+0.1) - 1;\n");
   fprintf(fp,"              u = cm[k+ind];\n");
   fprintf(fp,"            }\n");
   fprintf(fp,"            else u = 0.0;\n");
   fprintf(fp,"            if (t < 0.0) {\n");
   fprintf(fp,"              if (u == 0.0) u = 1.0;\n");
   fprintf(fp,"              else          u = 0.0;\n");
   fprintf(fp,"            }\n");
   fprintf(fp,"          } \n");
   fprintf(fp,"          if (u == 0.0) {\n");
   fprintf(fp,"            phi = 0.0;\n");
   fprintf(fp,"            break;\n");
   fprintf(fp,"          }\n");
   fprintf(fp,"          else {\n");
   fprintf(fp,"            phi *= u;\n");
   fprintf(fp,"            ip = floor(tb[ip*5+3] + 0.1) - 1;\n");
   fprintf(fp,"          }\n");
   fprintf(fp,"        }\n");
   fprintf(fp,"        Yt += tb[nn*5] * phi;\n");
   fprintf(fp,"      }\n");
   fprintf(fp,"    }\n");
   fprintf(fp,"    Y[ss] = Yt * Ystd + Ymean;\n");
   fprintf(fp,"  }\n");
   fprintf(fp,"  free(XX);\n");
   fprintf(fp,"  return 0;\n");
   fprintf(fp,"}\n");
   fprintf(fp,"/* ====================================*/\n");
   fprintf(fp,"int icat(double X, int input, double *cm) {\n");
   fprintf(fp,"  int j0, j1, j2, k, rdata;\n");
   fprintf(fp,"  j0 = floor(cm[2*input] + 0.1);\n");
   fprintf(fp,"  j1 = j0; j2 = floor(cm[2*input+1] + 0.1);\n");
   fprintf(fp,"  while (j2 != (j1+1)) {\n");
   fprintf(fp,"    k = floor(0.5*(j1+j2)) - 1;\n");
   fprintf(fp,"    if (cm[k] == X) {\n");
   fprintf(fp,"      rdata = k - j0; break;\n");
   fprintf(fp,"    }\n");
   fprintf(fp,"    else if (cm[k] < X) j1 = k;\n");
   fprintf(fp,"    else                j2 = k;\n");
   fprintf(fp,"    if (X == cm[j1-1]) rdata = j2 - j0;\n");
   fprintf(fp,"    else {\n");
   fprintf(fp,"      if (X == cm[j2-1]) rdata = j2 - j0;\n");
   fprintf(fp,"      else               rdata = 0;\n");
   fprintf(fp,"    }\n");
   fprintf(fp,"  }\n");
   fprintf(fp,"  return rdata;\n");
   fprintf(fp,"}\n");
   fprintf(fp,"/* ====================================*/\n");
   length = 3 + nBasisFcns_ * (5 * maxVarPerBasis_ + nSamples_ + 6) + 
            2 * nInputs_ + nSamples_;
   count = length;
   for (ii = length-1; ii >= 0; ii--)
   {
      if (fm_[ii] < 0.9*PSUADE_UNDEFINED) break;
      count--;
   }
   if (count > 0)
   {
      fprintf(fp,"/* \n");
      fprintf(fp,"FM %d \n",count);
      fprintf(fp,"*/ \n");
      fprintf(fp,"static double\n");
      fprintf(fp,"FM[%d] = \n", count);
      fprintf(fp,"{\n");
      for (ii = 0; ii < count-1; ii++) fprintf(fp," %24.16e ,\n", fm_[ii]);
      fprintf(fp," %24.16e };\n", fm_[count-1]);
   }
   length = 21 + nBasisFcns_ * (3 * maxVarPerBasis_ + 8);
   count = length;
   for (ii = length-1; ii >= 0; ii--)
   {
      if (im_[ii] != -9999) break;
      count--;
   }
   if (count > 0)
   {
      fprintf(fp,"/* \n");
      fprintf(fp,"IM %d \n",count);
      fprintf(fp,"*/ \n");
      fprintf(fp,"static int\n");
      fprintf(fp,"IM[%d] = \n", count);
      fprintf(fp,"{\n");
      for (ii = 0; ii < count-1; ii++) fprintf(fp," %d ,\n", im_[ii]);
      fprintf(fp," %d };\n", im_[count-1]);
   }
   fprintf(fp,"/* \n");
   fprintf(fp,"XM %d \n",count);
   fprintf(fp,"*/ \n");
   fprintf(fp,"static double\n");
   fprintf(fp,"XM[%d] = \n", nInputs_);
   fprintf(fp,"{\n");
   for (ii = 0; ii < nInputs_; ii++)
      fprintf(fp, " %24.16e , \n",XMeans_[ii]);
   fprintf(fp,"};\n");
   fprintf(fp,"/* \n");
   fprintf(fp,"XS %d \n",count);
   fprintf(fp,"*/ \n");
   fprintf(fp,"static double\n");
   fprintf(fp,"XS[%d] = \n", nInputs_);
   fprintf(fp,"{\n");
   for (ii = 0; ii < nInputs_; ii++)
      fprintf(fp, " %24.16e , \n",XStds_[ii]);
   fprintf(fp,"};\n");
   fprintf(fp,"/* ====================================*/\n");
   fprintf(fp,"int    *getIm() {return IM;}\n");
   fprintf(fp,"double *getFm() {return FM;}\n");
   fprintf(fp,"double *getXS() {return XS;}\n");
   fprintf(fp,"double *getXM() {return XM;}\n");
   fprintf(fp,"/* ====================================*/\n");
   fclose(fp);
   if (outputLevel_ > 0)
      printf("MARS parameters are now stored in psuade_rs.info\n");

   fp = fopen("psuade_rs.py", "w");
   if (fp == NULL)
   {
      printf("ERROR: Cannot open file psuade_rs.py.\n");
      return 0;
   }
   fwriteRSPythonHeader(fp);
   fprintf(fp,"#==================================================\n");
   fprintf(fp,"# MARS interpolation\n");
   fprintf(fp,"#==================================================\n");
   fwriteRSPythonCommon(fp);
   fprintf(fp,"nInputs = %d\n",nInputs_);
   length = 3 + nBasisFcns_ * (5 * maxVarPerBasis_ + nSamples_ + 6) +
            2 * nInputs_ + nSamples_;
   count = length;
   for (ii = length-1; ii >= 0; ii--)
   {
      if (fm_[ii] < 0.9*PSUADE_UNDEFINED) break;
      count--;
   }
   fprintf(fp,"FM = [\n"); 
   for (ii = 0; ii < count; ii++)
      fprintf(fp, "%16.8e ,\n", fm_[ii]);
   fprintf(fp,"]\n"); 
   length = 21 + nBasisFcns_ * (3 * maxVarPerBasis_ + 8);
   count = length;
   for (ii = length-1; ii >= 0; ii--)
   {
      if (im_[ii] != -9999) break;
      count--;
   }
   count = 15;
   fprintf(fp,"IM = [\n"); 
   for (ii = 0; ii < count; ii++)
      fprintf(fp, "%d ,\n", im_[ii]);
   fprintf(fp,"]\n"); 
   fprintf(fp,"XM = [\n");
   for (ii = 0; ii < nInputs_; ii++)
      fprintf(fp, " %24.16e , \n",XMeans_[ii]);
   fprintf(fp,"]\n");
   fprintf(fp,"XS = [\n");
   for (ii = 0; ii < nInputs_; ii++)
      fprintf(fp, " %24.16e , \n",XStds_[ii]);
   fprintf(fp,"]\n");
   fprintf(fp,"###################################################\n");
   fprintf(fp,"def icat(X, input, fm, indcm) :\n");
   fprintf(fp,"  j0 = int(fm[indcm+2*input] + 0.1)\n");
   fprintf(fp,"  j1 = j0 \n");
   fprintf(fp,"  j2 = int(fm[indcm+2*input+1] + 0.1)\n");
   fprintf(fp,"  while (j2 != (j1+1)) :\n");
   fprintf(fp,"    k = int(0.5*(j1+j2)+1.0e-8) - 1\n");
   fprintf(fp,"    if (fm[indcm+k] == X) :\n");
   fprintf(fp,"      rdata = k - j0;\n");
   fprintf(fp,"      break\n");
   fprintf(fp,"    elif (fm[indcm+k] < X) :\n");
   fprintf(fp,"      j1 = k\n");
   fprintf(fp,"    else :\n");
   fprintf(fp,"      j2 = k\n");
   fprintf(fp,"    if (X == fm[indcm+j1-1]) : \n");
   fprintf(fp,"      rdata = j2 - j0\n");
   fprintf(fp,"    else :\n");
   fprintf(fp,"      if (X == fm[indcm+j2-1]) : \n");
   fprintf(fp,"        rdata = j2 - j0\n");
   fprintf(fp,"      else :\n");
   fprintf(fp,"        rdata = 0\n");
   fprintf(fp,"  return rdata\n");
   fprintf(fp,"###################################################\n");
   fprintf(fp,"# Interpolation function  \n");
   fprintf(fp,"# X[0], X[1],   .. X[m-1]   - first point\n");
   fprintf(fp,"# X[m], X[m+1], .. X[2*m-1] - second point\n");
   fprintf(fp,"# ... \n");
   fprintf(fp,"#==================================================\n");
   fprintf(fp,"def interpolate(XX): \n");
   fprintf(fp,"  nSamp = int(len(XX) / nInputs + 1.0e-8)\n");
   fprintf(fp,"  X = nInputs * [0.0]\n");
   fprintf(fp,"  Ys = 2 * nSamp * [0.0]\n");
   fprintf(fp,"  for ss in range(nSamp) : \n");
   fprintf(fp,"    for ii in range(nInputs) : \n");
   fprintf(fp,"      X[ii] = XX[ss*nInputs+ii]\n");
   fprintf(fp,"    nk  = IM[4]\n");
   fprintf(fp,"    ind = IM[10] - 1; az  = FM[ind];\n");
   fprintf(fp,"    indtb = IM[11] - 1\n");
   fprintf(fp,"    indcm = IM[14] - 1\n");
   fprintf(fp,"    Yt = az;\n");
   fprintf(fp,"    for nn in range(nk) : \n");
   fprintf(fp,"      if (FM[indtb+nn*5] != 0.0) :\n");
   fprintf(fp,"        phi = 1.0\n");
   fprintf(fp,"        ip = nn\n");
   fprintf(fp,"        while (ip > -1) :\n");
   fprintf(fp,"          ind = int(indtb+ip*5)+1\n");
   fprintf(fp,"          t = FM[ind]\n");
   fprintf(fp,"          v = t\n");
   fprintf(fp,"          if (v < 0): \n");
   fprintf(fp,"            v = - v\n");
   fprintf(fp,"          ind = int(v+0.1) - 1\n");
   fprintf(fp,"          ind1 = indcm+2*ind\n");
   fprintf(fp,"          if (FM[ind1] <= 0.0) :\n");
   fprintf(fp,"            v = 1.0\n");
   fprintf(fp,"            if (t < 0.0) : \n");
   fprintf(fp,"              v = -1.0\n");
   fprintf(fp,"            dt = (X[ind] - XM[ind])/XS[ind]\n");
   fprintf(fp,"            ind2 = int(indtb+ip*5+2)\n");
   fprintf(fp,"            u = v * (dt-FM[ind2]) \n");
   fprintf(fp,"            if (u < 0.0) : \n");
   fprintf(fp,"              u = 0.0\n");
   fprintf(fp,"          else :\n");
   fprintf(fp,"            dt = (X[ind] - XM[ind])/XS[ind]\n");
   fprintf(fp,"            k = icat(dt, ind, FM, indcm)\n");
   fprintf(fp,"            if (k != 0) :\n");
   fprintf(fp,"              ind = int(FM[indtb+ip*5+2]+0.1) - 1\n");
   fprintf(fp,"              u = FM[indcm+k+ind]\n");
   fprintf(fp,"            else : \n");
   fprintf(fp,"              u = 0.0\n");
   fprintf(fp,"            if (t < 0.0) :\n");
   fprintf(fp,"              if (u == 0.0) : \n");
   fprintf(fp,"                u = 1.0\n");
   fprintf(fp,"              else : \n");
   fprintf(fp,"                u = 0.0\n");
   fprintf(fp,"          if (u == 0.0) :\n");
   fprintf(fp,"            phi = 0.0\n");
   fprintf(fp,"            break;\n");
   fprintf(fp,"          else :\n");
   fprintf(fp,"            phi *= u;\n");
   fprintf(fp,"            ind3 = int(indtb+ip*5+3)\n");
   fprintf(fp,"            ip = int(FM[ind3] + 0.1) - 1\n");
   fprintf(fp,"        Yt += FM[indtb+nn*5] * phi;\n");
   fprintf(fp,"    Ys[ss*2] = Yt * %24.16e + %24.16e\n",YStd_,YMean_);
   fprintf(fp,"  return Ys\n");
   fprintf(fp,"###################################################\n");
   fprintf(fp,"# main program\n");
   fprintf(fp,"#==================================================\n");
   fprintf(fp,"infileName  = sys.argv[1]\n");
   fprintf(fp,"outfileName = sys.argv[2]\n");
   fprintf(fp,"inputs = getInputData(infileName)\n");
   fprintf(fp,"outputs = interpolate(inputs)\n");
   fprintf(fp,"genOutputFile(outfileName, outputs)\n");
   fprintf(fp,"###################################################\n");
   fclose(fp);
   if (outputLevel_ > 0)
      printf("FILE psuade_rs.py contains the final MARS interpolator\n");
   return 0;
#else
   printf("PSUADE ERROR : Mars not installed.\n");
   return -1;
#endif
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int Mars::genNDGridData(double *XIn, double *YIn, int *NN, double **XX, 
                        double **YY)
{
#ifdef HAVE_MARS
   int    totPts, ss, ii;
   double *HX, *Xloc, **X2;

   initialize(XIn, YIn);
   if ((*NN) == -999) return 0;
  
   totPts = nPtsPerDim_;
   for (ii = 1; ii < nInputs_; ii++) totPts = totPts * nPtsPerDim_;
   HX = new double[nInputs_];
   for (ii = 0; ii < nInputs_; ii++) 
      HX[ii] = (upperBounds_[ii]-lowerBounds_[ii]) / (double) (nPtsPerDim_-1); 
 
   (*XX) = new double[nInputs_ * totPts];
   (*YY) = new double[totPts];
   (*NN) = totPts;
   X2 = new double*[totPts];
   for (ss = 0; ss < totPts; ss++) X2[ss] = new double[nInputs_];
   Xloc  = new double[nInputs_];
 
   for (ii = 0; ii < nInputs_; ii++) Xloc[ii] = lowerBounds_[ii];
   for (ss = 0; ss < totPts; ss++)
   {
      for (ii = 0; ii < nInputs_; ii++ ) X2[ss][ii] = Xloc[ii];
      for (ii = 0; ii < nInputs_; ii++ ) 
      {
         Xloc[ii] += HX[ii];
         if (Xloc[ii] < upperBounds_[ii] || 
              PABS(Xloc[ii] - upperBounds_[ii]) < 1.0E-7) break;
         else Xloc[ii] = lowerBounds_[ii];
      }
   }
 
   for (ii = 0; ii < nInputs_; ii++)
   {
      for (ss = 0; ss < totPts; ss++)
         X2[ss][ii] = (X2[ss][ii] - XMeans_[ii]) / XStds_[ii];
   } 
   if (outputLevel_ >= 2 || psMasterMode_ == 1)
   {
      printf("Entering Mars (fmod)\n");
      printf("If it crashes here, it is mars_fmod problem.\n");
   }
   mars_fmod(totPts, nInputs_, X2, *YY, fm_, im_);
   if (outputLevel_ >= 2 || psMasterMode_ == 1) 
      printf("Returned from Mars (fmod).\n");

   for (ii = 0; ii < totPts; ii++)
      (*YY)[ii] = ((*YY)[ii] * YStd_) + YMean_;

   delete [] Xloc;
   for (ss = 0; ss < totPts; ss++ ) delete [] X2[ss];
   delete [] X2;
   delete [] HX;
#else
   printf("PSUADE ERROR : Mars not installed.\n");
   return -1;
#endif
   return 0;
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int Mars::gen1DGridData(double *XIn, double *YIn, int ind1, 
                        double *settings, int *NN, double **XX, double **YY)
{
#ifdef HAVE_MARS
   int    ss, ii, totPts;
   double HX, **X2;

   initialize(XIn, YIn);
   if ((*NN) == -999) return 0;
  
   totPts = nPtsPerDim_;
   HX = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 

   X2 = new double*[totPts];
   for (ss = 0; ss < totPts; ss++) 
   {
      X2[ss] = new double[nInputs_];
      for (ii = 0; ii < nInputs_; ii++)
      {
         X2[ss][ii] = settings[ii]; 
         X2[ss][ii] = (X2[ss][ii] - XMeans_[ii]) / XStds_[ii];
      }
   }
   (*XX) = new double[totPts];
   (*YY) = new double[totPts];
   (*NN) = totPts;

   for (ss = 0; ss < totPts; ss++) 
   {
      X2[ss][ind1]  = HX * ss + lowerBounds_[ind1];
      X2[ss][ind1] = (X2[ss][ind1] - XMeans_[ind1]) / XStds_[ind1];
      (*XX)[ss]     = HX * ss + lowerBounds_[ind1];
      (*YY)[ss]     = 0.0;
   }

   if (outputLevel_ >= 2 || psMasterMode_ == 1)
   {
      printf("Entering Mars (fmod)\n");
      printf("If it crashes here, it is mars_fmod problem.\n");
   }
   mars_fmod(totPts, nInputs_, X2, *YY, fm_, im_);
   if (outputLevel_ >= 2 || psMasterMode_ == 1) 
      printf("Returned from Mars (fmod).\n");
 
   for (ii = 0; ii < totPts; ii++)
      (*YY)[ii] = ((*YY)[ii] * YStd_) + YMean_;

   for (ss = 0; ss < totPts; ss++) delete [] X2[ss];
   delete [] X2;
#else
   printf("PSUADE ERROR : Mars not installed.\n");
   return -1;
#endif
   return 0;
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int Mars::gen2DGridData(double *XIn, double *YIn, int ind1, int ind2, 
                        double *settings, int *NN, double **XX, double **YY)
{
#ifdef HAVE_MARS
   int    ss, ii, jj, totPts, index;
   double *HX, **X2;

   initialize(XIn, YIn);
   if ((*NN) == -999) return 0;
  
   totPts = nPtsPerDim_ * nPtsPerDim_;
   HX = new double[2];
   HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
   HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 

   X2 = new double*[totPts];
   for (ss = 0; ss < totPts; ss++) 
   {
      X2[ss] = new double[nInputs_];
      for (ii = 0; ii < nInputs_; ii++)
      {
         X2[ss][ii] = settings[ii]; 
         X2[ss][ii] = (X2[ss][ii] - XMeans_[ii]) / XStds_[ii];
      }
   }
   (*XX) = new double[2*totPts];
   (*YY) = new double[totPts];
   (*NN) = totPts;

   for (ii = 0; ii < nPtsPerDim_; ii++) 
   {
      for (jj = 0; jj < nPtsPerDim_; jj++)
      {
         index = ii * nPtsPerDim_ + jj;
         X2[index][ind1]  = HX[0] * ii + lowerBounds_[ind1];
         X2[index][ind1] = (X2[index][ind1] - XMeans_[ind1]) / XStds_[ind1];
         X2[index][ind2]  = HX[1] * jj + lowerBounds_[ind2];
         X2[index][ind2] = (X2[index][ind2] - XMeans_[ind2]) / XStds_[ind2];
         (*XX)[index*2]   = HX[0] * ii + lowerBounds_[ind1];
         (*XX)[index*2+1] = HX[1] * jj + lowerBounds_[ind2];
      }
   }

   if (outputLevel_ >= 2 || psMasterMode_ == 1) 
   {
      printf("Entering Mars (fmod)\n");
      printf("If it crashes here, it is mars_fmod problem.\n");
   }
   mars_fmod(totPts, nInputs_, X2, *YY, fm_, im_);
   if (outputLevel_ >= 2 || psMasterMode_ == 1)
      printf("Returned from Mars (fmod).\n");
   
   for (ii = 0; ii < totPts; ii++)
      (*YY)[ii] = ((*YY)[ii] * YStd_) + YMean_;

   for (ss = 0; ss < totPts; ss++) delete [] X2[ss];
   delete [] X2;
   delete [] HX;
#else
   printf("PSUADE ERROR : Mars not installed.\n");
   return -1;
#endif
   return 0;
}

// ************************************************************************
// Generate 3D results for display
// ------------------------------------------------------------------------
int Mars::gen3DGridData(double *XIn, double *YIn, int ind1, int ind2, int ind3,
                        double *settings, int *NN, double **XX, double **YY)
{
#ifdef HAVE_MARS
   int    ss, ii, jj, ll, totPts, index;
   double *HX, **X2;

   initialize(XIn, YIn);
   if ((*NN) == -999) return 0;
  
   totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
   HX = new double[3];
   HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
   HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 
   HX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_ - 1); 

   X2 = new double*[totPts];
   for (ss = 0; ss < totPts; ss++) 
   {
      X2[ss] = new double[nInputs_];
      for (ii = 0; ii < nInputs_; ii++) X2[ss][ii] = settings[ii]; 
   }
   (*XX) = new double[3*totPts];
   (*YY) = new double[totPts];
   (*NN) = totPts;

   for (ii = 0; ii < nPtsPerDim_; ii++) 
   {
      for (jj = 0; jj < nPtsPerDim_; jj++)
      {
         for (ll = 0; ll < nPtsPerDim_; ll++)
         {
            index = ii * nPtsPerDim_ * nPtsPerDim_ + jj * nPtsPerDim_ + ll;
            X2[index][ind1]  = HX[0] * ii + lowerBounds_[ind1];
            X2[index][ind1] = (X2[index][ind1] - XMeans_[ind1]) / XStds_[ind1];
            X2[index][ind2]  = HX[1] * jj + lowerBounds_[ind2];
            X2[index][ind2] = (X2[index][ind2] - XMeans_[ind2]) / XStds_[ind2];
            X2[index][ind3]  = HX[2] * ll + lowerBounds_[ind3];
            X2[index][ind3] = (X2[index][ind3] - XMeans_[ind3]) / XStds_[ind3];
            (*XX)[index*3]   = HX[0] * ii + lowerBounds_[ind1];
            (*XX)[index*3+1] = HX[1] * jj + lowerBounds_[ind2];
            (*XX)[index*3+2] = HX[2] * ll + lowerBounds_[ind3];
         }
      }
   }

   if (outputLevel_ >= 2 || psMasterMode_ == 1)
   {
      printf("Entering Mars (fmod)\n");
      printf("If it crashes here, it is mars_fmod problem.\n");
   }
   mars_fmod(totPts, nInputs_, X2, *YY, fm_, im_);
   if (outputLevel_ >= 2 || psMasterMode_ == 1)
      printf("Returned from Mars (fmod).\n");

   for (ii = 0; ii < totPts; ii++)
      (*YY)[ii] = ((*YY)[ii] * YStd_) + YMean_;

   for (ss = 0; ss < totPts; ss++) delete [] X2[ss];
   delete [] X2;
   delete [] HX;
#else
   printf("PSUADE ERROR : Mars not installed.\n");
   return -1;
#endif
   return 0;
}

// ************************************************************************
// Generate 4D results for display
// ------------------------------------------------------------------------
int Mars::gen4DGridData(double *XIn, double *YIn, int ind1, int ind2, int ind3,
                        int ind4, double *settings, int *NN, double **XX, 
                        double **YY)
{
#ifdef HAVE_MARS
   int    ss, ii, jj, ll, mm, totPts, index;
   double *HX, **X2;

   initialize(XIn, YIn);
   if ((*NN) == -999) 
   {
      (*XX) = NULL;
      (*YY) = NULL;
      return 0;
   }
  
   totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
   HX = new double[4];
   HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
   HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 
   HX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_ - 1); 
   HX[3] = (upperBounds_[ind4] - lowerBounds_[ind4]) / (nPtsPerDim_ - 1); 

   X2 = new double*[totPts];
   for (ss = 0; ss < totPts; ss++) 
   {
      X2[ss] = new double[nInputs_];
      for (ii = 0; ii < nInputs_; ii++) X2[ss][ii] = settings[ii]; 
   }
   (*XX) = new double[4*totPts];
   (*YY) = new double[totPts];
   (*NN) = totPts;

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
               X2[index][ind1]  = HX[0] * ii + lowerBounds_[ind1];
               X2[index][ind1]  = (X2[index][ind1]-XMeans_[ind1])/XStds_[ind1];
               X2[index][ind2]  = HX[1] * jj + lowerBounds_[ind2];
               X2[index][ind2]  = (X2[index][ind2]-XMeans_[ind2])/XStds_[ind2];
               X2[index][ind3]  = HX[2] * ll + lowerBounds_[ind3];
               X2[index][ind3]  = (X2[index][ind3]-XMeans_[ind3])/XStds_[ind3];
               X2[index][ind4]  = HX[3] * mm + lowerBounds_[ind4];
               X2[index][ind4]  = (X2[index][ind4]-XMeans_[ind4])/XStds_[ind4];
               (*XX)[index*4]   = HX[0] * ii + lowerBounds_[ind1];
               (*XX)[index*4+1] = HX[1] * jj + lowerBounds_[ind2];
               (*XX)[index*4+2] = HX[2] * ll + lowerBounds_[ind3];
               (*XX)[index*4+3] = HX[3] * mm + lowerBounds_[ind4];
            }
         }
      }
   }
   if (outputLevel_ >= 2 || psMasterMode_ == 1) 
   {
      printf("Entering Mars (fmod)\n");
      printf("If it crashes here, it is mars_fmod problem.\n");
   }
   mars_fmod(totPts, nInputs_, X2, *YY, fm_, im_);
   if (outputLevel_ >= 2 || psMasterMode_ == 1)
      printf("Returned from Mars (fmod).\n");

   for (ii = 0; ii < totPts; ii++)
      (*YY)[ii] = ((*YY)[ii] * YStd_) + YMean_;

   for (ss = 0; ss < totPts; ss++) delete [] X2[ss];
   delete [] X2;
   delete [] HX;
#else
   printf("PSUADE ERROR : Mars not installed.\n");
   return -1;
#endif
   return 0;
}

// ************************************************************************
// Generate results to a file called psuade_grid_data
// ------------------------------------------------------------------------
int Mars::writeToFileGrid2DData(double *XX, double *Y, int ind1,
                                int ind2, double *settings)
{
#ifdef HAVE_MARS
   int    ss, ii, jj, totPts, index;
   double *HX, **X2, dtmp, *YY, **X;
   FILE   *fp;

   X = new double*[nSamples_];
   for (ss = 0; ss < nSamples_; ss++) X[ss] = new double[nInputs_];
   for (ss = 0; ss < nSamples_; ss++) 
      for (ii = 0; ii < nInputs_; ii++) X[ss][ii] = XX[nInputs_*ss+ii];

   if (chooseWght_ == 1)
   {
      for (ss = 0; ss < nSamples_; ss++)
      {
         if (PABS(Y[ss]) < 1.0e-12) wgts_[ss] = 1.0;
         else                       wgts_[ss] = (float) PABS(Y[ss]);
      }
   }
   else
   {
      for (ss = 0; ss < nSamples_; ss++) wgts_[ss] = (float) weights_[ss];
   }

   if (outputLevel_ >= 2 || psMasterMode_ == 1)
   {
      printf("Entering Mars (process)\n");
      printf("If it crashes here, it is mars_process problem.\n");
      printf("One way to solve the problem is to use different nSamples.\n");
   }
   mars_process(nSamples_, nInputs_, X, Y, wgts_, nBasisFcns_,
                maxVarPerBasis_, varFlags_, fm_, im_);
   if (outputLevel_ >= 2 || psMasterMode_ == 1)
      printf("Returned from Mars (process)\n");
   for (ss = 0; ss < nSamples_; ss++) delete [] X[ss];
   delete [] X;

   totPts = nPtsPerDim_ * nPtsPerDim_;
   HX = new double[2];
   HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
   HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 

   X2 = new double*[totPts];
   for (ss = 0; ss < totPts; ss++) 
   {
      X2[ss] = new double[nInputs_];
      for (ii = 0; ii < nInputs_; ii++) X2[ss][ii] = settings[ii]; 
   }
   YY = new double[2 * totPts];

   for (ii = 0; ii < nPtsPerDim_; ii++) 
   {
      for (jj = 0; jj < nPtsPerDim_; jj++) 
      {
        index = ii * nPtsPerDim_ + jj;
        X2[index][ind1] = HX[0] * ii + lowerBounds_[ind1];
        X2[index][ind2] = HX[1] * jj + lowerBounds_[ind2];
      }
   }

   if (outputLevel_ >= 2 || psMasterMode_ == 1)
   {
      printf("Entering Mars (fmod)\n");
      printf("If it crashes here, it is mars_fmod problem.\n");
   }
   mars_fmod(totPts, nInputs_, X2, YY, fm_, im_);
   if (outputLevel_ >= 2 || psMasterMode_ == 1)
      printf("Returned from Mars (fmod)\n");

   fp = fopen("psuade_grid_data", "w");
   if(fp == NULL)
   {
      printf("fopen returned NULL in file %s line %d, exiting\n", 
             __FILE__, __LINE__);
      exit(1);
   }
   fprintf(fp, "%d\n",nPtsPerDim_);
   dtmp = lowerBounds_[ind1];
   for (ii = 0; ii < nPtsPerDim_; ii++) 
   {
      fprintf(fp, "%e\n", dtmp);
      dtmp += HX[0];
   }
   fprintf(fp, "%d\n",nPtsPerDim_);
   dtmp = lowerBounds_[ind2];
   for (ii = 0; ii < nPtsPerDim_; ii++) 
   {
      fprintf(fp, "%e\n", dtmp);
      dtmp += HX[1];
   }
   for (ii = 0; ii < nPtsPerDim_; ii++)
      for (jj = 0; jj < nPtsPerDim_; jj++)
         fprintf(fp, "%e\n", YY[ii*nPtsPerDim_+jj]);
   fclose(fp);

   for (ii = 0; ii < totPts; ii++) delete [] X2[ii];
   delete [] X2;
   delete [] YY;
   delete [] HX;
#else
   printf("PSUADE ERROR : Mars not installed.\n");
#endif
   return 0;
}

// ************************************************************************
// Evaluate a given point 
// ------------------------------------------------------------------------
double Mars::evaluatePoint(double *X)
{
   double Y=0.0;
#ifdef HAVE_MARS
   int    ii;
   double **X2;

   if (XMeans_ == NULL)
   {
      printf("PSUADE ERROR : Mars not initialized.\n");
      exit(1);
   }
   X2    = new double*[1];
   X2[0] = new double[nInputs_];
   for (ii = 0; ii < nInputs_; ii++)
      X2[0][ii] = (X[ii] - XMeans_[ii]) / XStds_[ii]; 
   mars_fmod(1, nInputs_, X2, &Y, fm_, im_);
   Y = Y * YStd_ + YMean_;
   delete [] X2[0];
   delete [] X2;
#else
   printf("PSUADE ERROR : Mars not installed.\n");
#endif
   return Y;
}

// ************************************************************************
// Evaluate a number of points
// ------------------------------------------------------------------------
double Mars::evaluatePoint(int npts, double *X, double *Y)
{
#ifdef HAVE_MARS
   int    ii, kk;
   double **X2;

   if (XMeans_ == NULL)
   {
      printf("PSUADE ERROR : Mars not initialized.\n");
      exit(1);
   }
   X2 = new double*[npts];
   for (kk = 0; kk < npts; kk++)
   {
      X2[kk] = new double[nInputs_];
      for (ii = 0; ii < nInputs_; ii++) 
         X2[kk][ii] = (X[kk*nInputs_+ii] - XMeans_[ii]) / XStds_[ii]; 
   }
   mars_fmod(npts, nInputs_, X2, Y, fm_, im_);
   for (kk = 0; kk < npts; kk++) Y[kk] = Y[kk] * YStd_ + YMean_;

   for (kk = 0; kk < npts; kk++) delete [] X2[kk];
   delete [] X2;
#else
   printf("PSUADE ERROR : Mars not installed.\n");
#endif
   return 0.0;
}

// ************************************************************************
// Evaluate a given point with standard deviation
// ------------------------------------------------------------------------
double Mars::evaluatePointFuzzy(double *X, double &std)
{
   double Y=0.0;
#ifdef HAVE_MARS
   int    ii;
   double **X2;

   if (XMeans_ == NULL)
   {
      printf("PSUADE ERROR : Mars not initialized.\n");
      exit(1);
   }
   X2    = new double*[1];
   X2[0] = new double[nInputs_];
   for (ii = 0; ii < nInputs_; ii++)
      X2[0][ii] = (X[ii] - XMeans_[ii]) / XStds_[ii]; 
   mars_fmod(1, nInputs_, X2, &Y, fm_, im_);
   Y = Y * YStd_ + YMean_;
   delete [] X2[0];
   delete [] X2;
#else
   printf("PSUADE ERROR : Mars not installed.\n");
#endif
   std = 0.0;
   return Y;
}

// ************************************************************************
// Evaluate a number of points with standard deviations
// ------------------------------------------------------------------------
double Mars::evaluatePointFuzzy(int npts,double *X, double *Y,double *Ystd)
{
#ifdef HAVE_MARS
   int    ii, kk;
   double **X2;

   if (XMeans_ == NULL)
   {
      printf("PSUADE ERROR : Mars not initialized.\n");
      exit(1);
   }
   X2 = new double*[npts];
   for (kk = 0; kk < npts; kk++)
   {
      X2[kk] = new double[nInputs_];
      for (ii = 0; ii < nInputs_; ii++)
         X2[kk][ii] = (X[kk*nInputs_+ii] - XMeans_[ii]) / XStds_[ii]; 
      Ystd[kk] = 0.0;
   }
   mars_fmod(npts, nInputs_, X2, Y, fm_, im_);
   for (kk = 0; kk < npts; kk++) Y[kk] = Y[kk] * YStd_ + YMean_;

   for (kk = 0; kk < npts; kk++) delete [] X2[kk];
   delete [] X2;
#else
   printf("PSUADE ERROR : Mars not installed.\n");
#endif
   return 0.0;
}

// ************************************************************************
// print MARS output information
// ------------------------------------------------------------------------
void Mars::printStat()
{
   int    lineLeng=500, i1, i2, nk=0, ind1, ind2;
   double *coefs;
   char   line[500], word1[100], word2[100], word3[100];
   FILE   *fp;
#if 0
   int    nGroups, *groupMembers, ivars, parent;
   double *vars, *knots, dtmp1, dtmp2, dtmp3, dtmp4, dtmp5, dtmp6;
#endif

   // first get the nk parameter
   fp = fopen(".psuade_mars", "r");
   if (fp == NULL) 
   {
      printf("Mars importance information not available.\n");
      return;
   }
   while ((fgets(line, lineLeng, fp) != NULL) && (feof(fp) == 0))
   {
      sscanf(line,"%s %s", word1, word2);
      if (!strcmp(word1,"input") && !strcmp(word2,"parameters"))
      {
         fgets(line, lineLeng, fp);
         fgets(line, lineLeng, fp);
         sscanf(line,"%d %d %d", &i1, &i2, &nk);
         break;
      }
   }
   fclose(fp);
   if (nk <= 0) 
   {
      printf("Mars::printStat ERROR - no. of basis functions = %d <= 0.\n",nk);
      return;
   }

   coefs = new double[nk+1];
   for (i1 = 0; i1 <= nk; i1++) coefs[i1] = -9999.0;
   while ((fgets(line, lineLeng, fp) != NULL) && (feof(fp) == 0))
   {
      sscanf(line,"%s %s", word1, word2);
      if (!strcmp(word1,"final") && !strcmp(word2,"model"))
      {
         for (i1 = 0; i1 <= nk; i1+=6)
         {
            fgets(line, lineLeng, fp);
            fscanf(fp, "%s", word3);
            if (strcmp(word3, "bsfn:")) 
            {
               nk = i1 - 1;
               break;
            }
            else
            {
               for (i2 = 0; i2 < 6; i2++)
               {
                  fscanf(fp, "%d", &ind1);
                  if (ind1 != (i1 + i2))
                  {
                     printf("Mars : nk reset from %d to %d\n", nk, i1+i2-1);
                     nk = i1 + i2 - 1; 
                     break;
                  }
               }
               if (i2 == 6) fgets(line, lineLeng, fp);
            } 
            fscanf(fp, "%s", word3);
            if (strcmp(word3, "coef:")) 
            {
               printf("Mars::printStat WARNING reading coefficients (%s).\n",
                      word3);
               printf("%s\n", line);
            }
            for (i2 = 0; i2 < 6; i2++) 
               if (i1+i2 <= nk) fscanf(fp, "%lg", &coefs[i1+i2]);
            fgets(line, lineLeng, fp);
         }
         break;
      }
   }
   fclose(fp);
   if (coefs[nk] == -9999.0)
   {
      printf("Mars::printStat ERROR - invalid Mars coefficients.\n");
      delete [] coefs;
      return;
   }
 
#if 0
   fp = fopen(".psuade_mars", "r");
   while ((fgets(line, lineLeng, fp) != NULL) && (feof(fp) == 0))
   {
      sscanf(line,"%s %s %s", word1, word2, word3);
      if (!strcmp(word1,"forward") && !strcmp(word2,"stepwise"))
      {
         nGroups = nk;
         groupMembers = new int*[nk];
         for (i1 = 0; i1 < nk; i1++) 
         {
            groupMembers[i1] = new int[2];
            for (i2 = 0; i2 < 2; i2++) groupMembers[i1][i2] = -1; 
         }
         knots = new double[nk+1];
         vars = new double[nk+1];
         ivars = new int[nk+1];
         parent = new int[nk+1];
         fgets(line, lineLeng, fp);
         fgets(line, lineLeng, fp);
         ind1 = 0;
         nGroups = 0;
         while (ind1 < nk)
         {
            fgets(line, lineLeng, fp);
            sscanf(line, "%d %lg", &i1, &dtmp1);
            i2 = (int) dtmp1;
            if ((dtmp1 - (double) i2) != 0) i2 = -1;
            if ((i1 - i2) != 1)
            {
               if (i1 == 0)
               {
                  sscanf(line, "%d %lg %lg %lg",&groupMembers[nGroups][0], 
                         &dtmp1, &dtmp2, &dtmp3); 
                  ivars[i1] = 0;
                  knots[i1] = 0.0;
                  parent[i1] = 0;
               }
               else
               {
                  sscanf(line, "%d %lg %lg %lg %lg %lg %lg", 
                      &groupMembers[nGroups][0], &dtmp1, &dtmp2, &dtmp3, 
                      &dtmp4, &dtmp5, &dtmp6);
                  ivars[i1] = (int) dtmp4;
                  knots[i1] = dtmp5;
                  parent[i1] = dtmp6;
               }
               nGroups++;
            }
            else
            {   
               sscanf(line, "%d %d %lg %lg %lg %lg %lg %lg", 
                      &groupMembers[nGroups][0], &groupMembers[nGroups][1],
                      &dtmp1, &dtmp2, &dtmp3, &dtmp4, &dtmp5, &dtmp6);
               ivars[i1] = ivars[ind1+1] = (int) dtmp4;
               knots[i1] = knots[ind1+1] = dtmp5;
               parent[i1] = parent[ind1+1] = (int) dtmp6;
               nGroups++;
            }
            ind1 = groupMembers[nGroups-1][0];
         }

         // analyze and print out

         for (i1 = 0; i1 < nGroups; i1++)
         {
            ind1 = groupMembers[i1][0];
            ind2 = groupMembers[i1][1];
            if (parent[ind1] == 0)
            {
               if (ind2 > 0 && coefs[ind2] != 0.0)
               {
                  printf("BF%4d : variable %3d > %9.3e has coefficient = %9.3e\n",
                         ind2, ivars[ind2], knots[ind2], coefs[ind2]);
               }
               if (ind1 > 0 && ind2 > 0 && coefs[ind1] != 0.0)
               {
                  printf("BF%4d : variable %3d < %9.3e has coefficient = %9.3e\n",
                         ind1, ivars[ind1], knots[ind1], coefs[ind1]);
               }
               else if (ind1 > 0 && ind2 < 0 && coefs[ind1] != 0.0)
               {
                  printf("BF%4d : variable %3d > %9.3e has coefficient = %9.3e\n",
                         ind1, ivars[ind1], knots[ind1], coefs[ind1]);
               }
            }
         }
         for (i1 = 0; i1 < nGroups; i1++)
         {
            ind1 = groupMembers[i1][0];
            if (parent[ind1] != 0)
            {
               if (coefs[ind1] != 0.0)
                  printf("Interaction: variables (%3d, %3d) with strength = %9.3e\n",
                         ivars[ind1], ivars[parent[ind1]], PABS(coefs[ind1]));
               ind1 = groupMembers[i1][1];
               if (ind1 > 0 && coefs[ind1] != 0.0)
                  printf("Interaction: variables (%3d, %3d) with strength = %9.3e\n",
                         ivars[ind1], ivars[parent[ind1]], PABS(coefs[ind1]));
            }
         }
         for (i1 = 0; i1 < nk; i1++) delete [] groupMembers[i1];
         delete [] groupMembers;
         delete [] knots;
         delete [] vars;
         delete [] ivars;
         delete [] parent;
         break;
      }
   }
   delete [] coefs;
   fclose(fp);
#endif

   fp = fopen(".psuade_mars", "r");
   if(fp == NULL)
   {
      printf("fopen returned NULL in file %s line %d, exiting\n",
             __FILE__, __LINE__);
      exit(1);
   }
   while ((fgets(line, lineLeng, fp) != NULL) && (feof(fp) == 0))
   {
      sscanf(line,"%s %s", word1, word2);
      if (!strcmp(word1,"anova") && !strcmp(word2,"decomposition"))
      {
         printf("\n Mars ANOVA Decomposition\n");
         fgets(line, lineLeng, fp);
         printf("%s", line);
         fgets(line, lineLeng, fp);
         ind2 = 1;
         sscanf(line, "%d", &ind1);
         while (ind2 == ind1)
         {
            printf("%s", line);
            fgets(line, lineLeng, fp);
            ind1 = -1;
            sscanf(line, "%d", &ind1);
            ind2++;
         }
         break;
      }
   }
   fclose(fp);
}

// ************************************************************************
// set parameters
// ------------------------------------------------------------------------
double Mars::setParams(int targc, char **targv)
{
   int    lineLeng=500, *iArray, ii, jj, nCount;
   double *dArray, dmax;
   char   line[500], word1[500], word2[500], word3[500];
   FILE   *fp, *fp2;

   if (targc > 0 && !strcmp(targv[0], "rank"))
   {
      fp = fopen(".psuade_mars", "r");
      if (fp != NULL)
      {
         nCount = 0;
         strcpy(word1, "none");
         dArray = new double[nInputs_];
         iArray = new int[nInputs_];
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
                           fscanf(fp,"%lg", &(dArray[ii+jj]));
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
             printf("* Mars gives insufficient importance measures.\n");
             printf("* To inspect Mars output, examine .psuade_mars.\n");
         }
         else
         {
            dmax = dArray[0];
            for (ii = 1; ii < nInputs_; ii++) 
               if (dArray[ii] > dmax) dmax = dArray[ii];
            for (ii = 0; ii < nInputs_; ii++) dArray[ii] /= (dmax / 100.0);

            if (psPlotTool_ == 1) fp2 = fopen("scilabmarsa.sci", "w");
            else                  fp2 = fopen("matlabmarsa.m", "w");
            if (fp2 == NULL)
            {
               printf("Mars ERROR: cannot open matlab/scilab file.\n");
               delete [] dArray;
               delete [] iArray;
               return -1.0;
            }
            fwritePlotCLF(fp2);
            fprintf(fp2, "n = %d;\n", nInputs_);
            fprintf(fp2, "Y = [\n");
            for (ii = 0; ii < nInputs_; ii++)
               fprintf(fp2, "%24.16e\n", dArray[ii]);
            fprintf(fp2, "];\n");
            fprintf(fp2, "ymax = max(Y);\n");
            fprintf(fp2, "ymin = min(Y);\n");
            fprintf(fp2, "if (ymax == ymin)\n");
            fprintf(fp2, "   ymax = ymax * 1.01;\n");
            fprintf(fp2, "   ymin = ymin * 0.99;\n");
            fprintf(fp2, "end;\n");
            fprintf(fp2, "if (ymax == ymin)\n");
            fprintf(fp2, "   ymax = ymax + 0.01;\n");
            fprintf(fp2, "   ymin = ymin - 0.01;\n");
            fprintf(fp2,"end;\n");
            fprintf(fp2, "bar(Y,0.8);\n");
            fwritePlotAxes(fp2);
            fwritePlotTitle(fp2, "MARS Rankings");
            fwritePlotXLabel(fp2, "Input parameters");
            fwritePlotYLabel(fp2, "MARS ranks");
            if (psPlotTool_ == 1)
            {
               fprintf(fp, "a=gca();\n");
               fprintf(fp, "a.data_bounds=[0, ymin-0.01*(ymax-ymin); ");
               fprintf(fp, "n+1, ymax+0.01*(ymax-ymin)];\n");
            }
            else
            {
               fprintf(fp,"axis([0 n+1 ymin-0.01*(ymax-ymin) ");
               fprintf(fp,"ymax+0.01*(ymax-ymin)])\n");
            }
            fclose(fp2);
            if (psPlotTool_ == 1)
                 printf("Mars ranking is now in scilabmarsa.sci.\n");
            else printf("Mars ranking is now in matlabmarsa.m.\n");

            for (ii = 0; ii < nInputs_; ii++) iArray[ii] = ii;
            sortDbleList2a(nInputs_, dArray, iArray);
            printAsterisks(PL_INFO, 0);
            printf("* Mars screening rankings \n");
            printAsterisks(PL_INFO, 0);
            for (ii = nInputs_-1; ii >= 0; ii--)
               printf("*  Rank %3d : Input = %3d (score = %4.1f)\n", 
                      nInputs_-ii, iArray[ii]+1, dArray[ii]);
            printAsterisks(PL_INFO, 0);
         }
         delete [] dArray;
         delete [] iArray;
         fclose(fp);
      }
      else printf("Mars: cannot rank - .psuade_mars file not found.\n");
   }
   else if (targc == 3 && !strcmp(targv[0], "mars_params"))
   {
      ii = *(int *) targv[1];
      if (ii > 0)
      {
         nBasisFcns_ = ii;
         if (nBasisFcns_ < 20) nBasisFcns_ = 20;
         if (outputLevel_ >= 2)
            printf("Mars: numBasis    set to = %d.\n", nBasisFcns_);
      }
      jj = *(int *) targv[2];
      if (jj > 0)
      {
         maxVarPerBasis_ = jj;
         if (maxVarPerBasis_ < 1) maxVarPerBasis_ = 1;
         if (outputLevel_ >= 2)
         printf("Mars: varPerBasis set to = %d.\n", maxVarPerBasis_);
      }
   }
   else if (targc == 1 && !strcmp(targv[0], "no_gen"))
   {
      noGen_ = 1;
   }
   else 
   {
      printf("Mars setParams ERROR: invalid command.\n");
   }
   return 0.0;
}

