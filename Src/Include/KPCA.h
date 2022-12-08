// ************************************************************************
// Copyright (c) 2007   Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory.
// Written by the PSUADE team.
// All rights reserved.
//
// Please see the COPYRIGHT_and_LICENSE file for the copyright notice,
// disclaimer, contact information and the GNU Lesser General Public License.
//
// DASSI is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License (as published by the Free Software
// Foundation) version 2.1 dated February 1999.
//
// DASSI is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the terms and conditions of the GNU General
// Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
// ************************************************************************
// KPCA functions
// DATE : 2015
// ************************************************************************
#ifndef __KPCAH__
#define __KPCAH__

/**
 * @name KPCA class
 *
 **/
/*@{*/

#include "psMatrix.h"
#include "psVector.h"

#define DA_GAUSSIAN   0
#define DA_LINEAR     1
#define DA_QUADRATIC  2
#define DA_CUBIC      3
#define DA_P4         4
#define DA_P5         5
 
// ************************************************************************
// class definition
// ************************************************************************
class KPCA
{
protected:
   int      kernelType_;
   int      rDim_;
   int      initialized_;
   int      outputLevel_;
   double   convergenceTol_;
   double   gaussScale_;
   double   trueScale_;

   psMatrix snapMat_;
   psMatrix featureMat_;
   psMatrix kernelMat_;
   psMatrix Umat_;
   psVector Svec_;
   psMatrix Vmat_;
   psMatrix eigenMatrix_;
   psVector eigenValues_;
   psVector kernelMeans_;

public:

   KPCA();
   ~KPCA();
   int setKernel(const int kernel);
   int setgaussScale(const float scale);
   int setReducedDimension(const int dim);
   int setConvergenceTol(double);
   int setOutputLevel(int);
   int loadSnapshots(psMatrix);
   int initialize(psMatrix &);
   int validateKPCA(psMatrix &, int);
   int genBetaFromStdRV(psVector &, psVector &);
   int getSnapshotDim();
   int getNumSnapshots();
   int getReducedDimension();
   int getSnapshot(int, psVector&);
   int getFeatureStdRVValues(psMatrix&);
   int computePreImage(psVector &, psVector &);
   int computePreImage(psVector &, psVector &, int);
   int genSnapshotFromStdRV(psVector &, psVector &);
   int projectFeatureToStdRV(psVector &, psVector &);
   int saveKPCA(char *);
   int retrieveKPCA(char *);
   int genKPCAModelFile(psMatrix &, int, int, char *);
   int inverseKPCA(char *, psVector, psVector &);
   int readSnapshotsFile(char *, psMatrix &); 
   int getKernel();

private:
   int iterate(psVector, psVector, psVector &);
   double computeInnerProduct(psVector &v1, psVector &v2);
   double computeInnerProductPrime(psVector &v1, psVector &v2);
   int getEigenMatrix(psMatrix &mat);
   int getEigenValues(psVector &vec);
   int genNormalSample(psVector &vec);
   int computeKernelMatrix();
};

/*@}*/

#endif /* __KPCAH__ */

