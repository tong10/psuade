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
// Functions for the class PsuadeSession  
// AUTHOR : CHARLES TONG
// DATE   : 2014
// ************************************************************************
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "PsuadeSession.h"
#include "PsuadeUtil.h"
#include "sysdef.h"
#include "PDFBase.h"
#include "pData.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
PsuadeSession::PsuadeSession()
{
  owned_ = 0;
  outputLevel_ = 0;
  nSamples_ = 0;
  nInputs_ = 0;
  nOutputs_ = 0;
  psuadeIO_ = NULL;
}

// ************************************************************************
// Copy constructor by Bill Oliver 
// ------------------------------------------------------------------------
PsuadeSession::PsuadeSession(const PsuadeSession &ps)
{
  int ii;

  outputLevel_ = ps.outputLevel_;
  nSamples_ = ps.nSamples_;
  nInputs_  = ps.nInputs_;
  nOutputs_ = ps.nOutputs_;
  rsType_   = ps.rsType_;
  owned_    = 1;
  psuadeIO_ = ps.psuadeIO_;

  vecInpLBounds_ = ps.vecInpLBounds_;
  vecInpUBounds_ = ps.vecInpUBounds_;
  vecSamInputs_  = ps.vecSamInputs_;
  vecSamOutputs_ = ps.vecSamOutputs_;
  vecSamStates_  = ps.vecSamStates_;
  vecInpPDFs_    = ps.vecInpPDFs_;
  vecInpMeans_   = ps.vecInpMeans_;
  vecInpStds_    = ps.vecInpStds_;
  inputNames_    = ps.inputNames_;
  outputNames_   = ps.outputNames_;
  corMatrix_.load((psMatrix &) ps.corMatrix_); 
}

// ************************************************************************
// operator= 
// ------------------------------------------------------------------------
PsuadeSession & PsuadeSession::operator=(const PsuadeSession & ps)
{
  int ii;

  if (this == &ps) return *this;

  outputLevel_ = ps.outputLevel_;
  nSamples_ = ps.nSamples_;
  nInputs_  = ps.nInputs_;
  nOutputs_ = ps.nOutputs_;
  rsType_   = ps.rsType_;
  owned_    = 1;
  psuadeIO_ = ps.psuadeIO_;

  vecInpLBounds_ = ps.vecInpLBounds_;
  vecInpUBounds_ = ps.vecInpUBounds_;
  vecSamInputs_  = ps.vecSamInputs_;
  vecSamOutputs_ = ps.vecSamOutputs_;
  vecSamStates_  = ps.vecSamStates_;
  vecInpPDFs_    = ps.vecInpPDFs_;
  vecInpMeans_   = ps.vecInpMeans_;
  vecInpStds_    = ps.vecInpStds_;
  inputNames_    = ps.inputNames_;
  outputNames_   = ps.outputNames_;
  corMatrix_.load((psMatrix &) ps.corMatrix_); 
  return *this;
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
PsuadeSession::~PsuadeSession()
{
  if (owned_)
  {
    inputNames_.clean();
    outputNames_.clean();
    vecInpLBounds_.clean();
    vecInpUBounds_.clean();
    vecSamInputs_.clean();
    vecSamOutputs_.clean();
    vecSamStates_.clean();
    vecInpPDFs_.clean();
    vecInpMeans_.clean();
    vecInpStds_.clean();
  }
  psuadeIO_ = NULL;
}

