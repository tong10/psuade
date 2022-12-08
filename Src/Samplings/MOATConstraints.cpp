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
// Functions for the class MOATConstraint
// AUTHOR : CHARLES TONG
// DATE   : 2009
// ************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include "PsuadeData.h"
#include "MOATConstraints.h"
#include "PsuadeUtil.h"
#define PABS(X) (((X) > 0)? X : -(X))

// ************************************************************************
// Constructor 
// ------------------------------------------------------------------------
MOATConstraints::MOATConstraints()
{
  nConstraints_ = 0;
  constraintFAs_ = NULL;
  constraintInputIndices_ = NULL;
  constraintInputValues_ = NULL;
  nInputs_ = 0;
  vecConstraintNInputs_.clean();
  vecXLBounds_.clean();
  vecXUBounds_.clean();
  vecYLBounds_.clean();
  vecYUBounds_.clean();
}

// ************************************************************************
// Copy Constructor by Bill Oliver
// ------------------------------------------------------------------------
MOATConstraints::MOATConstraints(const MOATConstraints & mc)
{
  nConstraints_ = mc.nConstraints_;
  constraintFAs_ = new FuncApprox*[nConstraints_];
  constraintInputIndices_ = new int*[nConstraints_];
  constraintInputValues_ = new double*[nConstraints_];
  nInputs_ = mc.nInputs_;
  sizeXLBounds_ = mc.sizeXLBounds_;
  sizeXUBounds_ = mc.sizeXUBounds_;
  vecXLBounds_ = mc.vecXLBounds_;
  vecXUBounds_ = mc.vecXUBounds_;
  vecYLBounds_ = mc.vecYLBounds_;
  vecYUBounds_ = mc.vecYUBounds_;
  vecConstraintNInputs_ = mc.vecConstraintNInputs_;
  for (int ii = 0; ii < nConstraints_; ii++)
  {
    constraintInputIndices_[ii] = mc.constraintInputIndices_[ii];
    constraintInputValues_[ii] = mc.constraintInputValues_[ii];
    constraintFAs_[ii] = mc.constraintFAs_[ii];
    for (int jj = 0; jj < vecConstraintNInputs_[ii]; jj++) 
    {
      constraintFAs_[ii][jj] = mc.constraintFAs_[ii][jj];
      constraintInputIndices_[ii][jj] = mc.constraintInputIndices_[ii][jj];
      constraintInputValues_[ii][jj] = mc.constraintInputValues_[ii][jj];
    }
  }
} 
   
// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
MOATConstraints::~MOATConstraints()
{
  int ii;

  if (constraintFAs_ != NULL)
  {
    for (ii = 0; ii < nConstraints_; ii++)
      if (constraintFAs_[ii] != NULL) delete constraintFAs_[ii];
    delete [] constraintFAs_;
  }
  if (constraintInputIndices_ != NULL)
  {
    for (ii = 0; ii < nConstraints_; ii++)
      if (constraintInputIndices_[ii] != NULL)
        delete [] constraintInputIndices_[ii];
    delete [] constraintInputIndices_;
  }
  if (constraintInputValues_ != NULL)
  {
    for (ii = 0; ii < nConstraints_; ii++)
      if (constraintInputValues_[ii] != NULL)
        delete [] constraintInputValues_[ii];
    delete [] constraintInputValues_;
  }
  vecConstraintNInputs_.clean();
  vecXLBounds_.clean();
  vecXUBounds_.clean();
  vecYLBounds_.clean();
  vecYUBounds_.clean();
}

// ************************************************************************
// generate all MOAT  constraints
// ------------------------------------------------------------------------
int MOATConstraints::initialize(PsuadeData *psuadeIO)
{
  int        printLevel, ii, jj, status, nInputsChk;
  double     *filterLBounds, *filterUBounds;
  char       **filterIndexFiles, **filterDataFiles;
  FILE       *fp;
  pData      pPtr, flbPtr, fubPtr, fdfPtr, fifPtr, iLPtr, iUPtr;
  PsuadeData *pIO;
  psIVector  vecSortArray;

  psuadeIO->getParameter("input_ninputs", pPtr);
  nInputs_ = pPtr.intData_;
  psuadeIO->getParameter("input_lbounds", iLPtr);
  sizeXLBounds_ = iLPtr.nDbles_;
  vecXLBounds_.load(sizeXLBounds_, iLPtr.dbleArray_);
  iLPtr.dbleArray_ = NULL;
  psuadeIO->getParameter("input_ubounds", iUPtr);
  sizeXUBounds_ = iUPtr.nDbles_;
  vecXUBounds_.load(sizeXUBounds_, iUPtr.dbleArray_);
  iUPtr.dbleArray_ = NULL;
  psuadeIO->getParameter("ana_diagnostics", pPtr);
  printLevel = pPtr.intData_;
  psuadeIO->getParameter("ana_num_moatfilters", pPtr);
  nConstraints_ = pPtr.intData_;
  if (nConstraints_ > 0)
  {
    if (printLevel > 0)
      printf("MOATConstraints: number of filters = %d\n",nConstraints_);
    psuadeIO->getParameter("ana_moatfilterlbounds", flbPtr);
    filterLBounds = flbPtr.dbleArray_;
    vecYLBounds_.setLength(nConstraints_);
    for (ii = 0; ii < nConstraints_; ii++) vecYLBounds_[ii] = filterLBounds[ii];
    psuadeIO->getParameter("ana_moatfilterubounds", fubPtr);
    filterUBounds = fubPtr.dbleArray_;
    vecYUBounds_.setLength(nConstraints_);
    for (ii = 0; ii < nConstraints_; ii++) vecYUBounds_[ii] = filterUBounds[ii];
    psuadeIO->getParameter("ana_moatfilterdatafile", fdfPtr);
    filterDataFiles = fdfPtr.strArray_;
    psuadeIO->getParameter("ana_moatfilterindexfile", fifPtr);
    filterIndexFiles = fifPtr.strArray_;
    for (ii = 0; ii < nConstraints_; ii++)
    {
      if (strcmp(filterDataFiles[ii], "NULL"))
      {
        fp = fopen(filterDataFiles[ii], "r");
        if (fp == NULL)
        {
          printf("MOATConstraints ERROR: filter data file %s not found.\n",
                 filterDataFiles[ii]);
          exit(1);
        }
        else fclose(fp);
      }
    }

    constraintFAs_ = new FuncApprox*[nConstraints_];
    vecConstraintNInputs_.setLength(nConstraints_);
    constraintInputIndices_ = new int*[nConstraints_];
    constraintInputValues_ = new double*[nConstraints_];
    for (ii = 0; ii < nConstraints_; ii++)
    {
      if (printLevel > 0)
        printf("MOATConstraints: initializing filter %d\n",ii+1);
      constraintFAs_[ii] = NULL;
      vecConstraintNInputs_[ii] = 0;
      constraintInputIndices_[ii] = NULL;
      constraintInputValues_[ii] = NULL;
      if (strcmp(filterDataFiles[ii], "NULL"))
      {
        pIO = new PsuadeData();
        status = pIO->readPsuadeFile(filterDataFiles[ii]);
        if (status != 0) 
        {
          printf("MOATConstraints ERROR: cannot read filter file %s.\n",
                 filterDataFiles[ii]);
          exit(1);
        }
        constraintFAs_[ii] = genFAInteractive(pIO, 2);
        pIO->getParameter("input_ninputs", pPtr);
        nInputsChk = pPtr.intData_;
        delete pIO;

        if (!strcmp(filterIndexFiles[ii],"NULL"))
        {
          printf("MOATConstraints ERROR: filter index file not given.\n");
          exit(1);
        }
        if (printLevel > 0)
          printf("MOATConstraints: filter %d has index file = %s\n",
                 ii+1, filterIndexFiles[ii]);
        fp = fopen(filterIndexFiles[ii],"r");
        if (fp != NULL)
        {
          fscanf(fp, "%d", &jj);
          vecConstraintNInputs_[ii] = jj;
          if (vecConstraintNInputs_[ii] != nInputsChk)
          {
            printf("MOATConstraints: filter %d must have %d(%d) inputs\n",
                   ii+1, nInputsChk, vecConstraintNInputs_[ii]);
            exit(1);
          }
          constraintInputIndices_[ii] = new int[vecConstraintNInputs_[ii]];
          constraintInputValues_[ii] = new double[vecConstraintNInputs_[ii]];
          for (jj = 0; jj < vecConstraintNInputs_[ii]; jj++)
          {
            fscanf(fp,"%d %lg",&constraintInputIndices_[ii][jj],
                         &(constraintInputValues_[ii][jj]));
            if (printLevel > 0)
              printf("MOATConstraints: filter %d has input %d = %d\n",
                              ii+1, jj+1, constraintInputIndices_[ii][jj]);
            if (constraintInputIndices_[ii][jj] >= 0)
              constraintInputIndices_[ii][jj]--;
            if (constraintInputIndices_[ii][jj] >= nInputs_ ||
                constraintInputIndices_[ii][jj] < -1)
            {
              printf("MOATConstraints ERROR: invalid filter input %d.\n",
                     constraintInputIndices_[ii][jj]);
              printf("MOATConstraints input format: ");
              printf("nInputs \n ");
              printf("<input (or 0)> <value (nominal val if 0)> \n ");
              printf("<input (or 0)> <value (nominal val if 0)> \n ");
              printf("... \n ");
              exit(1);
            }
          }
          vecSortArray.setLength(vecConstraintNInputs_[ii]);
          for (jj = 0; jj < vecConstraintNInputs_[ii]; jj++)
            vecSortArray[jj] = constraintInputIndices_[ii][jj];
          sortIntList(vecConstraintNInputs_[ii], vecSortArray.getIVector());
          for (jj = 1; jj < vecConstraintNInputs_[ii]; jj++)
          {
            if (vecSortArray[jj] >= 0 && vecSortArray[jj]==vecSortArray[jj-1])
            {
              printf("MOATConstraints ERROR: filter has duplicate");
              printf(" input %d.\n", vecSortArray[jj]);
              printf("MOATConstraints input format: ");
              printf("nInputs \n ");
              printf("<input (or 0)> <value (nominal val if 0)> \n ");
              printf("<input (or 0)> <value (nominal val if 0)> \n ");
              printf("... \n ");
              exit(1);
            }
          }
          fclose(fp);
        }
        else
        {
          printf("MOATConstraints ERROR : filter index file %s not found.\n",
                 filterIndexFiles[ii]);
          exit(1);
        }
      }
    }
  }
  return 0;
}

// ************************************************************************
// compute the MOAT scale
// ------------------------------------------------------------------------
double MOATConstraints::getScale(double *sampleInputs, int diffIndex, 
                                 int &flag)
{
  int    ii, jj, searchIndex;
  double filterRange, XL, XU, scale;
  double currX1, currX2, currY1, currY2, sLow, sHi;
  psVector vecSamplePt, vecXLs, vecXUs;

  if (nConstraints_ <= 0) 
  {
    flag = 1;
    return 0.0;
  }
  vecXLs.setLength(nConstraints_);
  vecXUs.setLength(nConstraints_);
  for (ii = 0; ii < nConstraints_; ii++)
  {
    vecSamplePt.setLength(vecConstraintNInputs_[ii]);
    for (jj = 0; jj < vecConstraintNInputs_[ii]; jj++)
      if (constraintInputIndices_[ii][jj] == -1)
        vecSamplePt[jj] = constraintInputValues_[ii][jj];
    filterRange = vecYUBounds_[ii] - vecYLBounds_[ii];
    searchIndex = -1;
    for (jj = 0; jj < vecConstraintNInputs_[ii]; jj++)
    {
      if (constraintInputIndices_[ii][jj] == diffIndex)
      {
        searchIndex = jj;
        break;
      }
    }
    if (searchIndex < 0) 
    {
      vecXLs[ii] = vecXLBounds_[diffIndex];
      vecXUs[ii] = vecXUBounds_[diffIndex];
    }
    else
    {
      for (jj = 0; jj < vecConstraintNInputs_[ii]; jj++)
        if (constraintInputIndices_[ii][jj] != -1)
          vecSamplePt[jj]=sampleInputs[constraintInputIndices_[ii][jj]];
      vecSamplePt[searchIndex] = vecXLBounds_[diffIndex];
      currY1 = constraintFAs_[ii]->evaluatePoint(vecSamplePt.getDVector());
      vecSamplePt[searchIndex] = vecXUBounds_[diffIndex];
      currY2 = constraintFAs_[ii]->evaluatePoint(vecSamplePt.getDVector());
      currX1 = vecXLBounds_[diffIndex];
      currX2 = vecXUBounds_[diffIndex];
      if (currY2 >= vecYUBounds_[ii] && currY1 >= vecYUBounds_[ii])
        currX1 = currX2 = 0.0;
      else if (currY2 <= vecYLBounds_[ii] && currY1 <= vecYLBounds_[ii])
        currX1 = currX2 = 0.0;
      else if (currY2 > currY1)
      {
        if (currY2 <= vecYUBounds_[ii]) currX2 = vecXUBounds_[diffIndex];
        else
        {
          sLow = vecXLBounds_[diffIndex];
          sHi  = vecXUBounds_[diffIndex];
          while (PABS((currY2-vecYUBounds_[ii])/filterRange)>0.001)
          {
            vecSamplePt[searchIndex] = 0.5 * (sLow + sHi);
            currY2 = constraintFAs_[ii]->evaluatePoint(vecSamplePt.getDVector());
            if (currY2 > vecYUBounds_[ii]) sHi  = 0.5 * (sLow + sHi);
            else                           sLow = 0.5 * (sLow + sHi);
          }
          currX2 = vecSamplePt[searchIndex];
        }
        if (currY1 >= vecYLBounds_[ii]) currX1 = vecXLBounds_[diffIndex];
        else
        {
          sLow = vecXLBounds_[diffIndex];
          sHi  = vecXUBounds_[diffIndex];
          while (PABS((currY1-vecYLBounds_[ii])/filterRange)>0.001)
          {
            vecSamplePt[searchIndex] = 0.5 * (sLow + sHi);
            currY1 = constraintFAs_[ii]->evaluatePoint(vecSamplePt.getDVector());
            if (currY1 < vecYLBounds_[ii]) sLow = 0.5 * (sLow + sHi);
            else                           sHi  = 0.5 * (sLow + sHi);
          }
          currX1 = vecSamplePt[searchIndex];
        }
      }
      else
      {
        if (currY1 <= vecYUBounds_[ii]) currX1 = vecXLBounds_[diffIndex];
        else
        {
          sLow = vecXLBounds_[diffIndex];
          sHi  = vecXUBounds_[diffIndex];
          while (PABS((currY1-vecYUBounds_[ii])/filterRange)>0.001)
          {
            vecSamplePt[searchIndex] = 0.5 * (sLow + sHi);
            currY1 = constraintFAs_[ii]->evaluatePoint(vecSamplePt.getDVector());
            if (currY1 > vecYUBounds_[ii]) sLow = 0.5 * (sLow + sHi);
            else                           sHi  = 0.5 * (sLow + sHi);
          }
          currX1 = vecSamplePt[searchIndex];
        }
        if (currY2 >= vecYLBounds_[ii]) currX2 = vecXUBounds_[diffIndex];
        else
        {
          sLow = vecXLBounds_[diffIndex];
          sHi  = vecXUBounds_[diffIndex];
          while (PABS((currY2-vecYLBounds_[ii])/filterRange)>0.001)
          {
            vecSamplePt[searchIndex] = 0.5 * (sLow + sHi);
            currY2 = constraintFAs_[ii]->evaluatePoint(vecSamplePt.getDVector());
            if (currY2 < vecYLBounds_[ii]) sHi  = 0.5 * (sLow + sHi);
            else                           sLow = 0.5 * (sLow + sHi);
          }   
          currX2 = vecSamplePt[searchIndex];
        }
      }
      if (currX1 < currX2)
      {
        vecXLs[ii] = currX1;
        vecXUs[ii] = currX2;
      }
      else
      {
        vecXLs[ii] = currX2;
        vecXUs[ii] = currX1;
      }
    }
  }
  XU = vecXUs[0];
  for (ii = 1; ii < nConstraints_; ii++) if (vecXUs[ii] < XU) XU = vecXUs[ii];
  XL = vecXLs[0];
  for (ii = 1; ii < nConstraints_; ii++) if (vecXLs[ii] > XL) XL = vecXLs[ii];
  scale = XU - XL;
  flag = 0;
  if (scale < 0.0) scale = 0.0;
  return scale;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
MOATConstraints& MOATConstraints::operator=(const MOATConstraints &mc)
{
  int ii, jj;
  if (this == &mc) return (*this);

  if (constraintFAs_ != NULL)
  {
    for (ii = 0; ii < nConstraints_; ii++)
      if (constraintFAs_[ii] != NULL) delete constraintFAs_[ii];
    delete [] constraintFAs_;
  }
  if (constraintInputIndices_ != NULL)
  {
    for (ii = 0; ii < nConstraints_; ii++)
      if (constraintInputIndices_[ii] != NULL)
        delete [] constraintInputIndices_[ii];
    delete [] constraintInputIndices_;
  }
  if (constraintInputValues_ != NULL)
  {
    for (ii = 0; ii < nConstraints_; ii++)
      if (constraintInputValues_[ii] != NULL)
        delete [] constraintInputValues_[ii];
    delete [] constraintInputValues_;
  }

  nConstraints_ = mc.nConstraints_;
  constraintFAs_ = new FuncApprox*[nConstraints_];
  constraintInputIndices_ = new int*[nConstraints_];
  constraintInputValues_ = new double*[nConstraints_];
  nInputs_ = mc.nInputs_;
  vecConstraintNInputs_ = mc.vecConstraintNInputs_;
  sizeXLBounds_ = mc.sizeXLBounds_;
  sizeXUBounds_ = mc.sizeXUBounds_;
  vecXLBounds_ = mc.vecXLBounds_;
  vecXUBounds_ = mc.vecXUBounds_;
  vecYLBounds_ = mc.vecYLBounds_;
  vecYUBounds_ = mc.vecYUBounds_;
  for (ii = 0; ii < nConstraints_; ii++) 
  {
    constraintInputIndices_[ii] = mc.constraintInputIndices_[ii];
    constraintInputValues_[ii] = mc.constraintInputValues_[ii];
    constraintFAs_[ii] = mc.constraintFAs_[ii];
    for (jj = 0; jj < vecConstraintNInputs_[ii]; jj++)
    {
      constraintFAs_[ii][jj] = mc.constraintFAs_[ii][jj];
      constraintInputIndices_[ii][jj] = mc.constraintInputIndices_[ii][jj];
      constraintInputValues_[ii][jj] = mc.constraintInputValues_[ii][jj];
    }
  }
  return (*this);
}

