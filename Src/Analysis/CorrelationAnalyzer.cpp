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
// Functions for the class CorrelationAnalyzer  
// AUTHOR : CHARLES TONG
// DATE   : 2005
// ************************************************************************
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

#include "CorrelationAnalyzer.h"
#include "sysdef.h"
#include "Psuade.h"
#include "PsuadeUtil.h"
#include "FuncApprox.h"
#include "PrintingTS.h"
#include "psVector.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
CorrelationAnalyzer::CorrelationAnalyzer() : Analyzer(), nInputs_(0), 
                                             nOutputs_(0) 
{
  setName("CORRELATION");
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
CorrelationAnalyzer::~CorrelationAnalyzer()
{
  cleanUp();
}

// ************************************************************************
// perform mean/variance analysis
// ------------------------------------------------------------------------
double CorrelationAnalyzer::analyze(aData &adata)
{
  int     ii, jj, ss, info, idata, count, iOne=1;
  double  xmean, xvar, ymean, yvar, yval, ddata, dmean, dvar;
  FILE    *fp;
  FuncApprox *faPtr=NULL;
  PsuadeData *ioPtr=NULL;
  pData   pData;
  psVector vecXMeans, vecYMeans, vecXVars, vecYVars, vecXVals, vecYVals; 
  psVector vecXLocal, vecYLocal, vecWLocal, vecXX, vecYY;

  //**/ ---------------------------------------------------------------
  // clean up
  //**/ ---------------------------------------------------------------
  cleanUp();

  //**/ ---------------------------------------------------------------
  // extract data
  //**/ ---------------------------------------------------------------
  nInputs_ = adata.nInputs_;
  int nInputs  = nInputs_;
  nOutputs_ = adata.nOutputs_;
  int nOutputs = nOutputs_;
  int nSamples = adata.nSamples_;
  double *X = adata.sampleInputs_;
  double *Y = adata.sampleOutputs_;
  int outputID = adata.outputID_;
  int printLevel = adata.printLevel_;
  if (adata.inputPDFs_ != NULL)
  {
    count = 0;
    for (ii = 0; ii < nInputs; ii++) count += adata.inputPDFs_[ii];
    if (count > 0)
    {
      printOutTS(PL_INFO, 
         "Correlation INFO: non-uniform probability distributions\n");
      printOutTS(PL_INFO, 
         "           have been defined in the data file, but they\n");
      printOutTS(PL_INFO,
         "           will not be used in this analysis.\n");
    }
  }
  ioPtr = adata.ioPtr_;
  if (ioPtr != NULL) ioPtr->getParameter("input_names", pData);

  //**/ ---------------------------------------------------------------
  // error checking
  //**/ ---------------------------------------------------------------
  if (nInputs <= 0 || nOutputs <= 0 || nSamples <= 0)
  {
    printOutTS(PL_ERROR, "Correlation ERROR: invalid arguments.\n");
    printOutTS(PL_ERROR, "    nInputs  = %d\n", nInputs);
    printOutTS(PL_ERROR, "    nOutputs = %d\n", nOutputs);
    printOutTS(PL_ERROR, "    nSamples = %d\n", nSamples);
    return PSUADE_UNDEFINED;
  } 
  if (outputID < 0 || outputID >= nOutputs)
  {
    printOutTS(PL_ERROR, "Correlation ERROR: invalid outputID.\n");
    printOutTS(PL_ERROR, "    nOutputs = %d\n", nOutputs);
    printOutTS(PL_ERROR, "    outputID = %d\n", outputID+1);
    return PSUADE_UNDEFINED;
  } 
  if (nSamples == 1)
  {
    printOutTS(PL_ERROR, 
         "Correlation INFO: analysis not meaningful for nSamples=1\n");
    return PSUADE_UNDEFINED;
  } 
  info = 0; 
  for (ss = 0; ss < nSamples; ss++)
    if (Y[nOutputs*ss+outputID] == PSUADE_UNDEFINED) info = 1;
  if (info == 1)
  {
    printOutTS(PL_ERROR,
         "Correlation ERROR: Some outputs are undefined.\n");
    printOutTS(PL_ERROR,
         "                   Prune the undefined's first.\n");
    return PSUADE_UNDEFINED;
  } 
   
  //**/ ---------------------------------------------------------------
  // first find the mean of the current set of samples
  //**/ ---------------------------------------------------------------
  printAsterisks(PL_INFO, 0);
  printOutTS(PL_INFO, "*                   Correlation Analysis\n");
  printEquals(PL_INFO, 0);
  printOutTS(PL_INFO, "*  Basic Statistics\n");
  printDashes(PL_INFO, 0);
  printOutTS(PL_INFO, "* Output of interest = %d\n", outputID+1);
  printDashes(PL_INFO, 0);
  computeMeanVariance(nSamples,nOutputs,Y,&ymean,&yvar,outputID,1);
  outputMean_ = ymean;
  outputVar_ = yvar;

  //**/ ---------------------------------------------------------------
  // compute the Pearson product moment correlation coefficient (PEAR)
  //**/ ---------------------------------------------------------------
  printEquals(PL_INFO, 0);
  printOutTS(PL_INFO,
       "*  Pearson correlation coefficients (PEAR) - linear -\n");
  printOutTS(PL_INFO,
       "*  which gives a measure of relationship between X_i's & Y.\n");
  printDashes(PL_INFO, 0);
  vecXMeans.setLength(nInputs);
  vecXVars.setLength(nInputs);
  vecXVals.setLength(nInputs);
  VecInpMeans_.setLength(nInputs_);
  VecInpVars_.setLength(nInputs_);
  VecInpPearsonCoef_.setLength(nInputs_);

  for (ii = 0; ii < nInputs; ii++)
  {
    computeMeanVariance(nSamples, nInputs, X, &dmean, &dvar, ii, 0);
    VecInpMeans_[ii] = vecXMeans[ii] = dmean;
    VecInpVars_[ii] = vecXVars[ii] = dvar;
  }
  computeCovariance(nSamples,nInputs,X,nOutputs,Y,vecXMeans.getDVector(),
           vecXVars.getDVector(),ymean,yvar,outputID,vecXVals.getDVector());
  for (ii = 0; ii < nInputs; ii++)
  {
    printOutTS(PL_INFO,
         "* Pearson Correlation coeff.  (Input %3d) = %e\n", 
         ii+1, vecXVals[ii]);
    VecInpPearsonCoef_[ii] = vecXVals[ii];
  }

  //**/ ---------------------------------------------------------------
  // now write these information to a plot file
  //**/ ---------------------------------------------------------------
  if (plotScilab())
  {
    fp = fopen("scilabca.sci","w");
    if (fp == NULL)
      printOutTS(PL_INFO,
            "CorrelationAnalysis: cannot write to scilab file.\n");
    else
      fprintf(fp,"// This file contains correlation coefficients.\n");
  }
  else
  {
    fp = fopen("matlabca.m","w");
    if (fp == NULL)
      printOutTS(PL_WARN, 
           "CorrelationAnalysis: cannot write to matlab file.\n");
    else
      fprintf(fp,"%% This file contains correlation coefficients.\n");
  }
  if (fp != NULL)
  {
    fprintf(fp, "sortFlag = 0;\n");
    fprintf(fp, "nn = %d;\n", nInputs);
    fprintf(fp, "PCC = [\n");
    for (ii = 0; ii < nInputs; ii++) fprintf(fp,"%24.16e\n", vecXVals[ii]);
    fprintf(fp, "];\n");
    if (pData.strArray_ != NULL)
    {
      if (plotScilab()) fprintf(fp, "  Str = [");
      else              fprintf(fp, "  Str = {");
      for (ii = 0; ii < nInputs-1; ii++) fprintf(fp,"'X%d',",ii+1);
      fprintf(fp,"'X%d'];\n",nInputs);
    }
    else
    {
      if (plotScilab()) fprintf(fp, "  Str = [");
      else              fprintf(fp, "  Str = {");
      for (ii = 0; ii < nInputs-1; ii++)
      {
        if (pData.strArray_[ii] != NULL) 
             fprintf(fp,"'%s',",pData.strArray_[ii]);
        else fprintf(fp,"'X%d',",ii+1);
      }
      if (plotScilab()) 
      {
        if (pData.strArray_[nInputs-1] != NULL) 
             fprintf(fp,"'%s'];\n",pData.strArray_[nInputs-1]);
        else fprintf(fp,"'X%d'];\n",nInputs);
      }
      else
      {
        if (pData.strArray_[nInputs-1] != NULL) 
             fprintf(fp,"'%s'};\n",pData.strArray_[nInputs-1]);
        else fprintf(fp,"'X%d'};\n",nInputs);
      }
    }
    fwritePlotCLF(fp);
    fprintf(fp, "if (sortFlag == 1)\n");
    if (plotScilab())
         fprintf(fp, "  [PCC, II] = gsort(PCC,'g','d');\n");
    else fprintf(fp, "  [PCC, II] = sort(PCC,'descend');\n");
    fprintf(fp, "  II   = II(1:nn);\n");
    fprintf(fp, "  PCC  = PCC(1:nn);\n");
    fprintf(fp, "  Str1 = Str(II);\n");
    fprintf(fp, "else\n");
    fprintf(fp, "  Str1 = Str;\n");
    fprintf(fp, "end\n");
    fprintf(fp, "ymin = min(PCC);\n");
    fprintf(fp, "ymax = max(PCC);\n");
    fprintf(fp, "h2 = 0.05 * (ymax - ymin);\n");
    fprintf(fp, "subplot(1,2,1)\n");
    fprintf(fp, "bar(PCC,0.8);\n");
    fwritePlotAxes(fp);
    fwritePlotTitle(fp,"Pearson Correlation Coefficients");
    fwritePlotYLabel(fp, "Correlation Coefficient");
    if (plotScilab())
    {
      fprintf(fp,"a=gca();\n");
      fprintf(fp,"a.data_bounds=[0, ymin; nn+1, ymax];\n");
      fprintf(fp,"a.x_ticks(2) = [1:nn]';\n");
      fprintf(fp,"a.x_ticks(3) = Str1';\n");
      fprintf(fp,"a.x_label.font_size = 3;\n");
      fprintf(fp,"a.x_label.font_style = 4;\n");
    }
    else
    {
      fprintf(fp,"axis([0 nn+1 ymin ymax])\n");
      fprintf(fp,"set(gca,'XTickLabel',[]);\n");
      fprintf(fp,"th=text(1:nn, repmat(ymin-0.05*(ymax-ymin),nn,1),");
      fprintf(fp,"Str1,'HorizontalAlignment','left','rotation',90);\n");
      fprintf(fp,"set(th, 'fontsize', 12)\n");
      fprintf(fp,"set(th, 'fontweight', 'bold')\n");
    }
  }

  //**/ ---------------------------------------------------------------
  // compute write these information to a plot file
  //**/ ---------------------------------------------------------------
  vecYMeans.setLength(nOutputs);
  vecYVars.setLength(nOutputs);
  vecYVals.setLength(nOutputs);

  VecOutMeans_.setLength(nOutputs_);
  VecOutVars_.setLength(nOutputs_);
  VecOutPearsonCoef_.setLength(nOutputs_);

  if (nOutputs > 1)
  {
    printEquals(PL_INFO, 0);
    printOutTS(PL_INFO, "*  PEAR (linear) for Y_i's versus Y.\n");
    printDashes(PL_INFO, 0);
    for (ii = 0; ii < nOutputs; ii++)
    {
      computeMeanVariance(nSamples, nOutputs, Y, &dmean, &dvar, ii, 0);
      VecOutMeans_[ii] = vecYMeans[ii] = dmean;
      VecOutVars_[ii] = vecYVars[ii] = dvar;
    }
    computeCovariance(nSamples,nOutputs,Y,nOutputs,Y,vecYMeans.getDVector(),
           vecYVars.getDVector(),ymean,yvar,outputID,vecYVals.getDVector());
    for (ii = 0; ii < nOutputs; ii++)
    {
      VecOutPearsonCoef_[ii] = vecYVals[ii];
      if (ii != outputID)
        printOutTS(PL_INFO, 
             "* Pearson Correlation coeff. (Output %2d) = %e\n", ii+1,
             vecYVals[ii]);
    }
  }
  printEquals(PL_INFO, 0);

  //**/ why did I do this?
  //**/YY = new double[nSamples];
  //**/for (ss = 0; ss < nSamples; ss++) YY[ss] = Y[ss*nOutputs+outputID];
  //**/faPtr = genFA(PSUADE_RS_REGR1, nInputs, nSamples);
  //**/faPtr->setOutputLevel(0);
  //**/faPtr->initialize(X, YY);
  //**/delete faPtr;
  //**/delete [] YY;

  //**/ ---------------------------------------------------------------
  //**/ compute the Spearman coefficient (SPEA)
  //**/ ---------------------------------------------------------------
  printOutTS(PL_INFO, 
       "*  Spearman coefficients (SPEA) - nonlinear relationship -  *\n");
  printOutTS(PL_INFO, 
       "*  which gives a measure of relationship between X_i's & Y. *\n");
  printOutTS(PL_INFO, 
       "*  (Idea: use the ranks for the inputs instead)             *\n");
  printDashes(PL_INFO, 0);
  vecXLocal.setLength(nSamples);
  vecYY.setLength(nSamples);
  VecInpSpearmanCoef_.setLength(nInputs_);
  for (ii = 0; ii < nInputs; ii++)
  {
    for (ss = 0; ss < nSamples; ss++)
    {
      vecXLocal[ss] = X[ss*nInputs+ii];
      vecYY[ss] = Y[ss*nOutputs+outputID];
    }
    sortDbleList2(nSamples, vecXLocal.getDVector(), vecYY.getDVector());
    for (ss = 0; ss < nSamples; ss++) vecXLocal[ss] = (double) ss;
    sortDbleList2(nSamples, vecYY.getDVector(), vecXLocal.getDVector());
    for (ss = 0; ss < nSamples; ss++) vecYY[ss] = (double) ss;
    computeMeanVariance(nSamples,1,vecXLocal.getDVector(),&xmean,&xvar,0,0);
    computeMeanVariance(nSamples,1,vecYY.getDVector(),&ymean,&yvar,0,0);
    computeCovariance(nSamples,1,vecXLocal.getDVector(),1,
           vecYY.getDVector(),&xmean,&xvar,ymean,yvar,0,&ddata);
    VecInpSpearmanCoef_[ii] = vecXVals[ii] = ddata;
    printOutTS(PL_INFO, 
         "* Spearman coefficient         (Input %3d ) = %e\n", ii+1,
         vecXVals[ii]);
  }
  if (fp != NULL)
  {
    fprintf(fp, "SPEA = [\n");
    for (ii = 0; ii < nInputs; ii++) fprintf(fp,"%24.16e\n", vecXVals[ii]);
    fprintf(fp, "];\n");
    fprintf(fp, "if (sortFlag == 1)\n");
    if (plotScilab())
         fprintf(fp, "  [SPEA, II] = gsort(SPEA,'g','d');\n");
    else fprintf(fp, "  [SPEA, II] = sort(SPEA,'descend');\n");
    fprintf(fp, "  II   = II(1:nn);\n");
    fprintf(fp, "  SPEA = SPEA(1:nn);\n");
    fprintf(fp, "  Str2 = Str(II);\n");
    fprintf(fp, "else\n");
    fprintf(fp, "  Str2 = Str;\n");
    fprintf(fp, "end\n");
    fprintf(fp, "ymin = min(SPEA);\n");
    fprintf(fp, "ymax = max(SPEA);\n");
    fprintf(fp, "h2 = 0.05 * (ymax - ymin);\n");
    fprintf(fp, "subplot(1,2,2)\n");
    fprintf(fp, "bar(SPEA,0.8);\n");
    fwritePlotAxes(fp);
    fwritePlotTitle(fp,"Spearman Correlation Coefficients");
    fwritePlotYLabel(fp, "Correlation Coefficient");
    if (plotScilab())
    {
      fprintf(fp, "a=gca();\n");
      fprintf(fp, "a.data_bounds=[0, ymin; nn+1, ymax];\n");
      fprintf(fp, "a.x_ticks(2) = [1:nn]';\n");
      fprintf(fp, "a.x_ticks(3) = Str2';\n");
      fprintf(fp, "a.x_label.font_size = 3;\n");
      fprintf(fp, "a.x_label.font_style = 4;\n");
      printEquals(PL_INFO, 0);
      printOutTS(PL_INFO, 
           " Correlation analysis plot file = scilabca.sci.\n");
    }
    else
    {
      fprintf(fp,"axis([0 nn+1 ymin ymax])\n");
      fprintf(fp,"set(gca,'XTickLabel',[]);\n");
      fprintf(fp,"th=text(1:nn,repmat(ymin-0.05*(ymax-ymin),nn,1),");
      fprintf(fp,"Str2,'HorizontalAlignment','left','rotation',90);\n");
      fprintf(fp,"set(th, 'fontsize', 12)\n");
      fprintf(fp,"set(th, 'fontweight', 'bold')\n");
      printEquals(PL_INFO, 0);
      printOutTS(PL_INFO,"Correlation analysis plot file = matlabca.m\n");
    }
    fclose(fp);
    fp = NULL;
  }
  printEquals(PL_INFO, 0);
  for (ii = 0; ii < nInputs; ii++) vecXMeans[ii] = (double) ii;
  for (ii = 0; ii < nInputs; ii++) vecXVals[ii] = PABS(vecXVals[ii]);
  sortDbleList2(nInputs, vecXVals.getDVector(), vecXMeans.getDVector());
  for (ii = nInputs-1; ii >= 0; ii--)
    printOutTS(PL_INFO,
         "* Spearman coefficient(ordered) (Input %3d ) = %e\n",
         (int) (vecXMeans[ii]+1), vecXVals[ii]);

  vecYLocal.setLength(nSamples);
  VecOutSpearmanCoef_.setLength(nOutputs_);
  if (nOutputs > 1)
  {
    printEquals(PL_INFO, 0);
    printOutTS(PL_INFO,
       "*  SPEA (nonlinear) for Y_i's versus Y.                     *\n");
    printDashes(PL_INFO, 0);
    for (ii = 0; ii < nOutputs; ii++)
    {
      if (ii != outputID)
      {
        for (ss = 0; ss < nSamples; ss++)
        {
          vecYLocal[ss] = Y[ss*nOutputs+ii];
          vecYY[ss] = Y[ss*nOutputs+outputID];
        }
        sortDbleList2(nSamples,vecYLocal.getDVector(),vecYY.getDVector());
        for (ss = 0; ss < nSamples; ss++) vecYLocal[ss] = (double) ss;
        sortDbleList2(nSamples,vecYY.getDVector(),vecYLocal.getDVector());
        for (ss = 0; ss < nSamples; ss++) vecYY[ss] = (double) ss;
        computeMeanVariance(nSamples,1,vecYLocal.getDVector(),&xmean,
                            &xvar,0,0);
        computeMeanVariance(nSamples,1,vecYY.getDVector(),&ymean,&yvar,0,0);
        computeCovariance(nSamples,1,vecYLocal.getDVector(),1,
               vecYY.getDVector(),&xmean,&xvar,ymean,yvar,0,&yval);
        VecOutSpearmanCoef_[ii] = vecYVals[ii] = yval;
        printOutTS(PL_INFO,
           "* Spearman coefficient        (Input %3d ) = %e\n", ii+1,
           vecYVals[ii]);
      }
    }
  }

  //**/ ---------------------------------------------------------------
  //**/ run Kendall tau rank correlation test
  //**/ ---------------------------------------------------------------
  VecInpKendallCoef_.setLength(nInputs_);
  if (printLevel > 1)
  {
    int nc=0, nd=0;
    printEquals(PL_INFO, 0);
    printOutTS(PL_INFO, "*  Kendall coefficients of concordance\n");
    printOutTS(PL_INFO, 
       "*  (Idea: use the ranks for both inputs and outputs)\n");
    printOutTS(PL_INFO, 
       "* Kendall coefficients (in the scale of -1 to 1) measure the strength\n");
    printOutTS(PL_INFO, 
       "* of relationship between an input and the output. The direction of\n");
    printOutTS(PL_INFO, 
       "* the relationship is indicated by the sign of the coefficient.\n");
    printOutTS(PL_INFO, 
       "* This nonparametric test is an alternative to Pearson's correlation\n");
    printOutTS(PL_INFO, 
       "* when the model input-output has a nonlinear relationship, and this\n");
    printOutTS(PL_INFO, 
       "* method is an alternative to the Spearman's correlation for small\n");
    printOutTS(PL_INFO, 
       "* sample size and there are many tied ranks.\n");
    printDashes(PL_INFO, 0);
    for (ii = 0; ii < nInputs; ii++)
    {
      nc = nd = 0;
      for (ss = 0; ss < nSamples; ss++)
      {
        vecXLocal.getDVector()[ss] = X[ss*nInputs+ii];
        vecYY[ss] = Y[ss*nOutputs+outputID];
      }
      //**/ reorder input ii based on magnitude of output
      sortDbleList2(nSamples, vecYY.getDVector(), vecXLocal.getDVector());
      for (ss = 0; ss < nSamples; ss++) vecYY[ss] = (double) ss;
      //**/ get the orderings from above
      sortDbleList2(nSamples, vecXLocal.getDVector(), vecYY.getDVector());
      for (ss = 0; ss < nSamples; ss++) vecXLocal[ss] = (double) ss;
      for (ss = 0; ss < nSamples; ss++)
      {
        for (jj = ss+1; jj < nSamples; jj++)
        {
          ddata = vecXLocal[ss] - vecXLocal[jj];
          if (ddata != 0.0)
          {
            ddata = (vecYY[ss] - vecYY[jj]) / ddata;
            if (ddata >= 0.0) nc++;
            else              nd++;
          }
          else
          {
            ddata = vecYY[ss] - vecYY[jj];
            if (ddata >= 0.0) nc++;
            else              nd++;
          }
        }
      }
      //printOutTS(PL_INFO, 
      //     "* Kendall coefficient         (Input %3d ) = %10.2e \n",
      //     ii+1, 2.0*(nc-nd)/(nSamples*(nSamples-1)));
      vecXVals[ii] = 2.0 * (nc - nd) / (nSamples * (nSamples - 1));
      VecInpKendallCoef_[ii] = vecXVals[ii];
    }
    printDashes(PL_INFO, 0);
    for (ii = 0; ii < nInputs; ii++) vecXMeans[ii] = (double) ii;
    for (ii = 0; ii < nInputs; ii++) vecXVals[ii] = PABS(vecXVals[ii]);
    sortDbleList2(nInputs, vecXVals.getDVector(), vecXMeans.getDVector());
    for (ii = nInputs-1; ii >= 0; ii--)
      printOutTS(PL_INFO,
           "* Kendall coefficient(ordered) (Input %3d ) = %e\n",
           (int) (vecXMeans[ii]+1), vecXVals[ii]);
  }

#if 0
  //**/ ---------------------------------------------------------------
  //**/ run linear regression analysis on the rank-ordered data
  //**/ ---------------------------------------------------------------
  if (printLevel > 2)
  {
    printEquals(PL_INFO, 0);
    printOutTS(PL_INFO, 
      "*  Regression analysis on rank-ordered inputs/outputs       *\n");
    vecXX.setLength(nSamples*nInputs);
    vecWLocal.setLength(nSamples);
    for (ss = 0; ss < nSamples*nInputs; ss++) vecXX[ss] = X[ss];
    for (ss = 0; ss < nSamples; ss++)
    {
      vecYY[ss] = Y[ss*nOutputs+outputID];
      vecWLocal[ss] = (double) ss;
    }
    sortDbleList2(nSamples, vecYY.getDVector(), vecWLocal.getDVector());
    for (ss = 0; ss < nSamples; ss++)
    {
      idata = (int) vecWLocal[ss];
      vecYY[idata] = (double) (ss + 1);
    }
    for (ii = 0; ii < nInputs; ii++)
    {
      for (ss = 0; ss < nSamples; ss++)
      {
        vecWLocal[ss] = (double) ss;
        vecXLocal[ss] = vecXX[nSamples*ii+ss];
      }
      sortDbleList2(nSamples,vecXLocal.getDVector(),vecWLocal.getDVector());
      for (ss = 0; ss < nSamples; ss++)
      {
        idata = (int) vecWLocal[ss];
        vecXX[ii*nSamples+idata] = (double) (ss + 1);
      }
    }
    faPtr = genFA(PSUADE_RS_REGR1, nInputs, iOne, nSamples);
    faPtr->setOutputLevel(-1);
    faPtr->initialize(vecXX.getDVector(), vecYY.getDVector());
    delete faPtr;
  }
#endif
  printAsterisks(PL_INFO, 0);
  return 0.0;
}

// *************************************************************************
// Compute mean and variance
// -------------------------------------------------------------------------
int CorrelationAnalyzer::computeMeanVariance(int nSamples, int xDim, 
              double *X, double *xmean, double *xvar, int xID, int flag)
{
  int    ss;
  double mean, variance;

  mean = 0.0;
  for (ss = 0; ss < nSamples; ss++) mean += X[xDim*ss+xID];
  mean /= (double) nSamples;
  variance = 0.0;
  for (ss = 0; ss < nSamples; ss++) 
    variance += ((X[xDim*ss+xID] - mean) * (X[xDim*ss+xID] - mean));
  variance /= (double) (nSamples - 1);
  (*xmean) = mean;
  (*xvar)  = variance;
  if (flag == 1)
  {
    printOutTS(PL_INFO, "Correlation: mean     = %e\n", mean);
    printOutTS(PL_INFO, "Correlation: variance = %e\n", variance);
  }
  return 0;
}

// *************************************************************************
// Compute agglomerated covariances
// -------------------------------------------------------------------------
int CorrelationAnalyzer::computeCovariance(int nSamples,int nX,double *X,
             int nY, double *Y, double *xmeans, double *xvars, double ymean,
             double yvar, int yID, double *Rvalues)
{
  int    ii, ss;
  double denom, numer;

  for (ii = 0; ii < nX; ii++)
  {
    numer = 0.0;
    for (ss = 0; ss < nSamples; ss++)
      numer += ((X[ss*nX+ii] - xmeans[ii]) * (Y[ss*nY+yID] - ymean));
    numer /= (double) (nSamples - 1);
    denom = sqrt(xvars[ii] * yvar);
    if (denom == 0.0)
    {
      printOutTS(PL_INFO,"Correlation ERROR: denom=0 for input %d\n",ii+1);
      printOutTS(PL_INFO, 
           "denom = xvar * yvar : xvar = %e, yvar = %e\n",xvars[ii],yvar);
      Rvalues[ii] = 0.0;
    }
    else Rvalues[ii] = numer / denom;
  }
  return 0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
CorrelationAnalyzer& CorrelationAnalyzer::operator=(const CorrelationAnalyzer &)
{
  printOutTS(PL_ERROR,"Correlation operator= ERROR: operation not allowed.\n");
  exit(1);
  return (*this);
}

// ************************************************************************
// functions for getting results
// ------------------------------------------------------------------------
int CorrelationAnalyzer::get_nInputs()
{
  return nInputs_;
}
int CorrelationAnalyzer::get_nOutputs()
{
  return nOutputs_;
}
double CorrelationAnalyzer::get_outputMean()
{
  return outputMean_;
}
double CorrelationAnalyzer::get_outputVar()
{
  return outputVar_;
}
double *CorrelationAnalyzer::get_inputMeans()
{
  psVector vecT;
  vecT = VecInpMeans_;
  return vecT.takeDVector();
}
double *CorrelationAnalyzer::get_inputVars()
{
  psVector vecT;
  vecT = VecInpVars_;
  return vecT.takeDVector();
}
double *CorrelationAnalyzer::get_outputMeans()
{
  psVector vecT;
  vecT = VecOutMeans_;
  return vecT.takeDVector();
}
double *CorrelationAnalyzer::get_outputVars()
{
  psVector vecT;
  vecT = VecOutVars_;
  return vecT.takeDVector();
}
double *CorrelationAnalyzer::get_inputPearsonCoef()
{
  psVector vecT;
  vecT = VecInpPearsonCoef_;
  return vecT.takeDVector();
}
double *CorrelationAnalyzer::get_outputPearsonCoef()
{
  psVector vecT;
  vecT = VecOutPearsonCoef_;
  return vecT.takeDVector();
}
double *CorrelationAnalyzer::get_inputSpearmanCoef()
{
  psVector vecT;
  vecT = VecInpSpearmanCoef_;
  return vecT.takeDVector();
}
double *CorrelationAnalyzer::get_outputSpearmanCoef()
{
  psVector vecT;
  vecT = VecOutSpearmanCoef_;
  return vecT.takeDVector();
}
double *CorrelationAnalyzer::get_inputKendallCoef()
{
  psVector vecT;
  vecT = VecInpKendallCoef_;
  return vecT.takeDVector();
}
int CorrelationAnalyzer::cleanUp()
{
  nInputs_ = 0;
  nOutputs_ = 0;
  outputMean_ = 0.0;
  outputVar_ = 0.0;
  VecInpMeans_.clean();
  VecInpVars_.clean();
  VecOutMeans_.clean();
  VecOutVars_.clean();
  VecInpPearsonCoef_.clean();
  VecOutPearsonCoef_.clean();
  VecInpSpearmanCoef_.clean();
  VecOutSpearmanCoef_.clean();
  VecInpKendallCoef_.clean();
  return 0;
}

