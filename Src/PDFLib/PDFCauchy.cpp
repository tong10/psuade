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
// Functions for the Cauchy distribution
// AUTHOR : CHARLES TONG
// DATE   : 2016
// ************************************************************************
#include <stdio.h>
#include <math.h>
#include "Psuade.h"
#include "PsuadeUtil.h"
#include "PDFCauchy.h"
#define PABS(x) (((x) >= 0) ? x : -(x))
#define PS_PI 3.14159265358

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
PDFCauchy::PDFCauchy(double x0, double gamma)
{
  X0_    = x0;
  gamma_ = gamma;
  if (gamma <= 0)
  {
    printf("PDFCauchy: gamma needs to be > 0.\n");
    exit(1);
  }
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
PDFCauchy::~PDFCauchy()
{
}

// ************************************************************************
// forward transformation 
// ------------------------------------------------------------------------
int PDFCauchy::getPDF(int length, double *inData, double *outData)
{
  int    ii;
  double mult, xdata;

  if (psConfig_.PDFDiagnosticsIsOn())
    printf("PDFCauchy: getPDF begins (length = %d)\n",length);
  mult = 1.0 / (PS_PI * gamma_);
  for (ii = 0; ii < length; ii++)
  {
    xdata = (inData[ii] - X0_) / gamma_;
    outData[ii] = mult / (1.0 + xdata * xdata);
  }
  if (psConfig_.PDFDiagnosticsIsOn()) printf("PDFCauchy: getPDF ends.\n");
  return 0;
}

// ************************************************************************
// look up cumulative density
// ------------------------------------------------------------------------
int PDFCauchy::getCDF(int length, double *inData, double *outData)
{
  int    ii;
  double ddata;

  if (psConfig_.PDFDiagnosticsIsOn())
    printf("PDFCauchy: getCDF begins (length = %d)\n",length);
  for (ii = 0; ii < length; ii++)
  {
    ddata = (inData[ii] - X0_) / gamma_;
    outData[ii] = atan(ddata) / PS_PI + 0.5;
  }
  if (psConfig_.PDFDiagnosticsIsOn()) printf("PDFCauchy: getCDF ends.\n");
  return 0;
}

// ************************************************************************
// transformation to range
// ------------------------------------------------------------------------
int PDFCauchy::invCDF(int length, double *inData, double *outData)
{
  int    ii;
  double ddata, dtm, xlo, xmi, xhi, ylo, ymi, yhi;

  //**/ -------------------------------------------------------------
  //**/ map the input data onto the CDF
  //**/ -------------------------------------------------------------
  if (psConfig_.PDFDiagnosticsIsOn())
    printf("PDFCauchy: invCDF begins (length = %d)\n",length);

  //**/ -------------------------------------------------------------
  //**/ search for xlo and xhi
  //**/ -------------------------------------------------------------
  double XLO = 0.0;
  double YLO = 1;
  while (YLO > 1.0e-8)
  {
    XLO -= 1.0;
    dtm = (XLO - X0_) / gamma_;
    YLO = atan(dtm) / PS_PI + 0.5;
  }
  double XHI = 0.0;
  double YHI = 0.0;
  while (YHI < 0.9999999)
  {
    XHI += 1.0;
    dtm = (XHI - X0_) / gamma_;
    YHI = atan(dtm) / PS_PI + 0.5;
  }

  //**/ -------------------------------------------------------------
  //**/ now do the search
  //**/ -------------------------------------------------------------
  for (ii = 0; ii < length; ii++)
  {
    ddata = inData[ii];
    if (ddata <= 0.0 || ddata >= 1)
    {
      printf("PDFCauchy invCDF ERROR - CDF value %e not in (0,1).\n",
             ddata);
      exit(1);
    }
    xlo = XLO;
    ylo = YLO;
    xhi = XHI;
    yhi = YHI;
    if      (ddata <= ylo) outData[ii] = xlo;
    else if (ddata >= yhi) outData[ii] = xhi;
    else
    {
      while (PABS(ddata-ylo) > 1.0e-12 || PABS(ddata-yhi) > 1.0e-12)
      {
        xmi = 0.5 * (xhi + xlo);
        dtm = (xmi - X0_) / gamma_;
        ymi = atan(dtm) / PS_PI + 0.5;
        if (ddata > ymi) 
        {
          xlo = xmi;
          ylo = ymi;
        }
        else
        {
          xhi = xmi;
          yhi = ymi;
        }
      }
      if (PABS(ddata-ylo) < PABS(ddata-yhi)) outData[ii] = xlo;
      else                                   outData[ii] = xhi;
    }
  }
  if (psConfig_.PDFDiagnosticsIsOn()) printf("PDFCauchy: invCDF ends.\n");
  return 0;
}

// ************************************************************************
// generate a sample
// ------------------------------------------------------------------------
int PDFCauchy::genSample(int length, double *outData, double *lowers,
                         double *uppers)
{
  int    ii;
  double UU, xhi, xlo, xmi, yhi, ylo, ymi, mult, lower, upper;

  //**/ -------------------------------------------------------------
  //**/ upper and lower bounds has to be in (0,1), and upper > lower
  //**/ -------------------------------------------------------------
  if (lowers == NULL || uppers == NULL)
  {
    printf("PDFCauchy genSample ERROR - lower/upper bound unavailable.\n"); 
    exit(1);
  }
  lower = lowers[0];
  upper = uppers[0];
  if (upper <= lower)
  {
    printf("PDFCauchy genSample ERROR - lower bound >= upper bound.\n");
    exit(1);
  }

  //**/ -------------------------------------------------------------
  //**/ generate sample
  //**/ -------------------------------------------------------------
  if (psConfig_.PDFDiagnosticsIsOn())
    printf("PDFCauchy: genSample begins (length = %d)\n",length);
  for (ii = 0; ii < length; ii++)
  {
    UU = PSUADE_drand();
    xlo = lower;
    mult = (xlo - X0_) / gamma_;
    ylo = atan(mult) / PS_PI + 0.5;
    xhi = upper;
    mult = (xhi - X0_) / gamma_;
    yhi = atan(mult) / PS_PI + 0.5;
    UU = UU * (yhi - ylo) + ylo;
    if      (UU <= ylo) outData[ii] = xlo;
    else if (UU >= yhi) outData[ii] = xhi;
    else
    {
      while (PABS(UU-ylo) > 1.0e-12 || PABS(UU-yhi) > 1.0e-12)
      {
        xmi = 0.5 * (xhi + xlo);
        mult = (xmi - X0_) / gamma_;
        ymi = atan(mult) / PS_PI + 0.5;
        if (UU > ymi) 
        {
          xlo = xmi;
          ylo = ymi;
        }
        else
        {
          xhi = xmi;
          yhi = ymi;
        }
      }
      if (PABS(UU-ylo) < PABS(UU-yhi)) outData[ii] = xlo;
      else                             outData[ii] = xhi;
    }
  }
  if (psConfig_.PDFDiagnosticsIsOn()) 
    printf("PDFCauchy: genSample ends.\n");
  return 0;
}

// ************************************************************************
// get mean
// ------------------------------------------------------------------------
double PDFCauchy::getMean()
{
  printf("PDFCauchy ERROR: mean is not defined.\n");
  return 0.0;
}

