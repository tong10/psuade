// ************************************************************************
// Copyright (c) 2007   Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory.
// Written by the Charles Tong <tong10@llnl.gov>.
//
// CODE UCRL-CODE-235523 All rights reserved.
//
// This file is part of PSUADE.
//
// Please see the GPL.pdf file for the copyright notice,
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License (as published by the
// Free Software Foundation) version 2.1 dated February 1991.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the terms and conditions of the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
// ************************************************************************
// Functions for the class PDFManager
// AUTHOR : CHARLES TONG
// DATE   : 2008
// ************************************************************************

#ifndef __PDFMANAGERH__
#define __PDFMANAGERH__

#include "Util/Matrix.h"
#include "Util/Vector.h"
#include "DataIO/PsuadeData.h"
#include "PDFLib/PDFBase.h"
#include "PDFLib/PDFMVNormal.h"
#include "PDFLib/PDFMVLogNormal.h"

/**
 * @name class for probability density functions
 *
 **/

// ************************************************************************
// class definition
// ************************************************************************
class PDFManager 
{

   int nInputs_;
   int *pdfMap_;
   int nGNormal_;
   int *gnormalInputs_;
   int nGLognormal_;
   int *glognormalInputs_;
   int printLevel_;
   PDFBase        **PDFptrs_;
   PDFMVNormal    *PDFMVNormalPtr_;
   PDFMVLogNormal *PDFMVLogNormalPtr_;
   PsuadeData     *PSIOptr_;
   Vector         means_;
   Matrix         corMat_;

public:

   PDFManager();    
   ~PDFManager();
   int initialize(PsuadeData *);
   int initialize(int, int *, double *, double *, Matrix &);
   int getPDF(int, Vector &, Vector &, Vector &, Vector &);
   int genSample(int, Vector &, Vector &, Vector &);
   int genSample();
   int invCDF(int, Vector &, Vector &, Vector &, Vector &);
   int getCDF(int, Vector &, Vector &, Vector &, Vector &);
};

#endif // __PDFMANAGERH__

