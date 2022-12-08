// ************************************************************************
// Copyright (c) 2007   Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory.
// Written by the PSUADE team. 
// All rights reserved.
//
// Please see the COPYRIGHT and LICENSE file for the copyright notice,
// disclaimer, contact information and the GNU Lesser General Public License.
//
// PSUADE is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free 
// Software Foundation) version 2.1 dated February 1999.
//
// PSUADE is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the terms and conditions of the GNU Lesser
// General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
// ************************************************************************
// Functions for the class GowerAnalyzer (Gower's distance analysis) 
// ************************************************************************
// AUTHOR : CHARLES TONG
// DATE   : 2009
// ************************************************************************
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "GowerAnalyzer.h"
#include "MainEffectAnalyzer.h"
#include "PsuadeUtil.h"
#include "sysdef.h"
#include "Psuade.h"
#include "FuncApprox.h"
#include "Sampling.h"
#include "PrintingTS.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
GowerAnalyzer::GowerAnalyzer() : Analyzer()
{
  setName("GOWER");
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
GowerAnalyzer::~GowerAnalyzer()
{
}

// ************************************************************************
// perform mean/variance analysis
// ------------------------------------------------------------------------
double GowerAnalyzer::analyze(aData &adata)
{
  int    ss, ss2, ii, iZero=0, nLHSample=100000, nLHSSub=500, iOne=1;
  double gower, dmax, dmin, ddata, vsum;
  char   dataFile[500], lineIn[500];
  FILE   *fp;
  pData  pPtr;

  //**/ ---------------------------------------------------------------
  //**/ extract data
  //**/ ---------------------------------------------------------------
  int printLevel = adata.printLevel_;
  int nSamples   = adata.nSamples_;
  int nInputs    = adata.nInputs_;
  int nOutputs   = adata.nOutputs_;
  int outputID   = adata.outputID_;
  double *XIn    = adata.sampleInputs_;
  double *YIn    = adata.sampleOutputs_;
  double *lowerB = adata.iLowerB_;
  double *upperB = adata.iUpperB_;

  psVector  vecX2, vecY2;
  psIVector vecS2;
  vecY2.setLength(nSamples);
  for (ii = 0; ii < nSamples; ii++) 
    vecY2[ii] = YIn[ii*nOutputs+outputID];

  //**/ ---------------------------------------------------------------
  //**/ generate sample and evaluate with response surface 
  //**/ ---------------------------------------------------------------
  //FuncApprox *faPtr = genFA(PSUADE_RS_MARS, nInputs, iOne, nSamples);
  FuncApprox *faPtr = genFA(-1, nInputs, iOne, nSamples);
  int status = faPtr->initialize(XIn, vecY2.getDVector());

  Sampling *sampPtr=(Sampling *) SamplingCreateFromID(PSUADE_SAMP_LHS);
  sampPtr->setInputBounds(nInputs, lowerB, upperB);
  sampPtr->setOutputParams(1);
  sampPtr->setSamplingParams(nLHSample, 200, 0);
  sampPtr->initialize(0);
  vecX2.setLength(nLHSample*nInputs);
  vecY2.setLength(nLHSample);
  vecS2.setLength(nLHSample);
  sampPtr->getSamples(nLHSample, nInputs, 1, vecX2.getDVector(), 
                      vecY2.getDVector(), vecS2.getIVector());
  faPtr->evaluatePoint(nLHSample, vecX2.getDVector(), 
                       vecY2.getDVector());

  //**/ ---------------------------------------------------------------
  //**/ create ME and study main effect
  //**/ ---------------------------------------------------------------
  MainEffectAnalyzer *mePtr = new MainEffectAnalyzer();
  
  psVector vecVCE, vecVVM, vecVVV, vecMVV;
  vecVCE.setLength(nInputs);
  vecVVM.setLength(nInputs);
  vecVVV.setLength(nInputs);
  vecMVV.setLength(nInputs);
  mePtr->computeVCE(nInputs, nLHSample, nLHSSub, vecX2.getDVector(), 
         vecY2.getDVector(), iZero, NULL, vecVVM.getDVector(), 
         vecMVV.getDVector(),vecVVV.getDVector(),vecVCE.getDVector());

  vsum = 0.0;
  for (ii = 0; ii < nInputs; ii++) vsum += vecVCE[ii];
  if (vsum <= 0.0)
  {
    printOutTS(PL_INFO,"GowerAnalyzer INFO: vce sum <= 0.0.\n");
    printOutTS(PL_INFO,"                    vce not used for scaling.");
    for (ii = 0; ii < nInputs; ii++) vecVCE[ii] = 1.0;
  }
  else
  {
    printOutTS(PL_INFO,"GowerAnalyzer INFO: vce used for scaling.\n");
    for (ii = 0; ii < nInputs; ii++) vecVCE[ii] /= vsum;
    for (ii = 0; ii < nInputs; ii++)
      printOutTS(PL_INFO,  "VCE %4d = %12.4e\n", ii+1, vecVCE[ii]);
  } 
  delete faPtr;
  delete sampPtr;
  delete mePtr;

  //**/ ---------------------------------------------------------------
  //**/ get data to be predicted
  //**/ ---------------------------------------------------------------
  printf("Next enter a sample for prediction (in PSUADE format).\n");
  printf("Enter file name for unobserved data (outside convex hull): ");
  scanf("%s", dataFile);
  fp = fopen(dataFile,"r");
  if (fp == NULL)
  {
    printOutTS(PL_ERROR,"GowerAnalyzer ERROR: cannot open file %s.\n", 
               dataFile);
    return PSUADE_UNDEFINED;
  }
  fclose(fp);
  fgets(lineIn,500,stdin);

  PsuadeData *pIO = new PsuadeData();
  pIO->setOutputLevel(0);
  status = pIO->readPsuadeFile(dataFile);
  if (status != 0)
  {
    printf("ERROR: cannot read file %s in PSUADE format.\n",dataFile);
    exit(1);
  }
  pIO->getParameter("input_ninputs", pPtr);
  int nInputs2 = pPtr.intData_;
  pIO->getParameter("method_nsamples", pPtr);
  int nSamples2 = pPtr.intData_;
  pIO->getParameter("input_sample", pPtr);
  double *X3 = pPtr.dbleArray_;
  if (nInputs != nInputs2)
  {
    printOutTS(PL_ERROR,
         "GowerAnalyzer ERROR: different input dimensions %d %d\n",
         nInputs, nInputs2);
    delete pIO;
    return PSUADE_UNDEFINED;
  }

  psVector vecRanges;
  vecRanges.setLength(nInputs);
  for (ii = 0; ii < nInputs; ii++)
  {
    dmax = - PSUADE_UNDEFINED;
    dmin =   PSUADE_UNDEFINED;
    for (ss = 0; ss < nSamples; ss++)
    {
      if (XIn[ss*nInputs+ii] < dmin) dmin = XIn[ss*nInputs+ii];
      if (XIn[ss*nInputs+ii] > dmax) dmax = XIn[ss*nInputs+ii];
    }
    vecRanges[ii] = dmax - dmin;
    if (vecRanges[ii] == 0.0)
    {
      printOutTS(PL_ERROR,
         "GowerAnalyzer ERROR: some input range = 0 (%d).\n",ii+1);
      delete pIO;
      return PSUADE_UNDEFINED;
    }
  }

  //**/ ---------------------------------------------------------------
  //**/ Gower distance processing 
  //**/ ---------------------------------------------------------------
  if (plotScilab())
  {
    fp = fopen("psuade_gower_data.sci", "w");
    if (fp != NULL) fprintf(fp, "// Gower distance statistics\n");
  }
  else
  {
    fp = fopen("psuade_gower_data.m", "w");
    if (fp != NULL) fprintf(fp, "%% Gower distance statistics\n");
  }
  if (fp == NULL)
  {
    printOutTS(PL_ERROR,  
       "GowerAnalyzer ERROR: cannot write to psuade_gower_data.m\n");
    delete pIO; 
    return 1.0;
  }

  fprintf(fp, "G = [\n");
  for (ss = 0; ss < nSamples; ss++)
  {
    for (ss2 = 0; ss2 < nSamples2; ss2++)
    {
      gower = 0.0;
      for (ii = 0; ii < nInputs; ii++)
        gower += PABS(XIn[ss*nInputs+ii]-X3[ss2*nInputs+ii])/
                 vecRanges[ii];
      gower /= (double) nInputs;
      fprintf(fp, "%e ", gower);
    }
    fprintf(fp, "\n");
  }
  fprintf(fp,"];\n");
  fwritePlotCLF(fp);
  fprintf(fp,"for ii = 1 : %d\n", nSamples2);
  fprintf(fp,"   X = G(:,ii);\n");
  if (plotScilab())
       fprintf(fp,"   X = gsort(X,'g','i');\n");
  else fprintf(fp,"   X = sort(X);\n");
  fprintf(fp,"   Y = [1:%d]' / %d;\n", nSamples, nSamples);
  fprintf(fp,"   plot(X,Y)\n");
  if (plotMatlab())
  {
    fprintf(fp,"   if (ii == 1)\n");
    fprintf(fp,"      hold on\n");
    fprintf(fp,"   end;\n");
  }
  fprintf(fp,"end;\n");
  fwritePlotAxes(fp);
  fwritePlotXLabel(fp, "Gower distance");
  fwritePlotYLabel(fp, "Probabilities");

  //**/ ---------------------------------------------------------------
  //**/ scaled Gower distance 
  //**/ ---------------------------------------------------------------
  fprintf(fp, "\n");
  if (plotMatlab())
  {
    fprintf(fp, "figure(2);\n");
    fprintf(fp, "%% Scaled Gower distance statistics\n");
  }
  else
  {
    fprintf(fp, "scf(2);\n");
    fprintf(fp, "// Scaled Gower distance statistics\n");
  }
  fprintf(fp, "G2 = [\n");
  for (ss = 0; ss < nSamples; ss++)
  {
    for (ss2 = 0; ss2 < nSamples2; ss2++)
    {
      gower = 0.0;
      for (ii = 0; ii < nInputs; ii++)
        gower += PABS(XIn[ss*nInputs+ii]-X3[ss2*nInputs+ii]) / 
                 vecRanges[ii] * vecVCE[ii];
      gower /= (double) nInputs;
      fprintf(fp, "%e ", gower);
    }
    fprintf(fp, "\n");
  }
  fprintf(fp,"];\n");
  fwritePlotCLF(fp);
  fprintf(fp,"for ii = 1 : %d\n", nSamples2);
  fprintf(fp,"   X = G2(:,ii);\n");
  if (plotScilab())
       fprintf(fp,"   X = gsort(X,'g','i');\n");
  else fprintf(fp,"   X = sort(X);\n");
  fprintf(fp,"   Y = [1:%d]' / %d;\n", nSamples, nSamples);
  fprintf(fp,"   plot(X,Y)\n");
  if (plotMatlab())
  {
    fprintf(fp,"   if (ii == 1)\n");
    fprintf(fp,"      hold on\n");
    fprintf(fp,"   end;\n");
  }
  fprintf(fp,"end;\n");
  fwritePlotAxes(fp);
  fwritePlotXLabel(fp, "Scaled Gower distance");
  fwritePlotYLabel(fp, "Probabilities");

  //**/ ---------------------------------------------------------------
  //**/ Mahalanobis distance processing 
  //**/ ---------------------------------------------------------------
  fprintf(fp, "\n");
  if (plotMatlab())
  {
    fprintf(fp,"figure(3);\n");
    fprintf(fp,"%% Mahalanobis distance statistics\n");
    fprintf(fp,"%% showing distance from centers of clusters/std dev\n");
  }
  else
  {
    fprintf(fp,"scf(3);\n");
    fprintf(fp,"// Mahalanobis distance statistics\n");
    fprintf(fp,"// showing distance from centers of clusters/std dev\n");
  }
  fprintf(fp, "M = [\n");

  psVector vecDMeans, vecDVars;
  vecDMeans.setLength(nInputs);
  vecDVars.setLength(nInputs);

  for (ii = 0; ii < nInputs; ii++)
  {
    vecDMeans[ii] = 0.0;
    for (ss = 0; ss < nSamples; ss++) vecDMeans[ii] += XIn[ss*nInputs+ii];
    vecDMeans[ii] /= nSamples;
    vecDVars[ii] = 0.0;
    for (ss = 0; ss < nSamples; ss++) 
      vecDVars[ii] += pow(XIn[ss*nInputs+ii] - vecDMeans[ii], 2.0);
    vecDVars[ii] /= nSamples;
    if (vecDVars[ii] == 0.0)
      printf("GowerAnalyzer WARNING: input %d has zero variance.\n",ii);
  }
  for (ss2 = 0; ss2 < nSamples2; ss2++)
  {
    ddata = 0.0;
    for (ii = 0; ii < nInputs; ii++)
    {
      if (vecDVars[ii] != 0.0)
        ddata += pow(X3[ss2*nInputs+ii]-vecDMeans[ii], 2.0)/vecDVars[ii];
      else
      {
        ddata = PSUADE_UNDEFINED;
        break;
      }
    }
    if (ddata != PSUADE_UNDEFINED) ddata = sqrt(ddata);
    fprintf(fp, "%e\n", ddata);
  }
  fprintf(fp,"];\n");
  fwritePlotCLF(fp);
  fprintf(fp,"X = [1:%d]';\n", ss2);
  fprintf(fp,"plot(X,M,'b*')\n");
  fwritePlotAxes(fp);
  fwritePlotXLabel(fp, "Sample Number");
  fwritePlotYLabel(fp, "Mahalanobis distance");
  fclose(fp);
  if (plotScilab())
     printOutTS(PL_INFO,
          "GowerAnalyzer: psuade_gower_data.sci file created.\n");
  else
     printOutTS(PL_INFO,
           "GowerAnalyzer: psuade_gower_data.m file created.\n");

  //**/ ---------------------------------------------------------------
  //**/ clean up
  //**/ ---------------------------------------------------------------
  delete pIO;
  return 0.0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
GowerAnalyzer& GowerAnalyzer::operator=(const GowerAnalyzer &)
{
   printOutTS(PL_ERROR,
              "GowerAnalyzer operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

