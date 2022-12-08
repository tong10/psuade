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
// Functions for the class RSConstraints
// AUTHOR : CHARLES TONG
// DATE   : 2006
// ************************************************************************

#include <stdio.h>
#include <stdlib.h>
#include "PsuadeData.h"
#include "RSConstraints.h"
#include "PsuadeUtil.h"
#include "Globals.h"

// ************************************************************************
// Constructor 
// ------------------------------------------------------------------------
RSConstraints::RSConstraints()
{
  nConstraints_ = 0;
  constraintFAs_ = NULL;
  constraintNInputs_ = NULL;
  constraintInputIndices_ = NULL;
  constraintInputValues_ = NULL;
  nInputs_ = 0;
  vecLBounds_.clean();
  vecUBounds_.clean();
}
   
// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
RSConstraints::~RSConstraints()
{
  int ii;

  if (constraintFAs_ != NULL)
  {
    for (ii = 0; ii < nConstraints_; ii++)
      if (constraintFAs_[ii] != NULL) delete constraintFAs_[ii];
    delete [] constraintFAs_;
  }
  if (constraintNInputs_ != NULL) delete [] constraintNInputs_;
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
  vecLBounds_.clean();
  vecUBounds_.clean();
}

// ************************************************************************
// return number of constraints
// ------------------------------------------------------------------------
int RSConstraints::getNumConstraints()
{
  return nConstraints_;
}

// ************************************************************************
// generate all output constraints
// ------------------------------------------------------------------------
int RSConstraints::genConstraints(PsuadeData *psuadeIO)
{
  int        printLevel, ii, jj, status, nInputsChk;
  char       **filterIndexFiles, **filterDataFiles;
  FILE       *fp;
  pData      pPtr, flbPtr, fubPtr, fdfPtr, fifPtr;
  PsuadeData *pIO;
  psIVector  vecSortArray;

  psuadeIO->getParameter("input_ninputs", pPtr);
  nInputs_ = pPtr.intData_;
  psuadeIO->getParameter("ana_diagnostics", pPtr);
  printLevel = pPtr.intData_;
  psuadeIO->getParameter("ana_num_rsfilters", pPtr);
  nConstraints_ = pPtr.intData_;
  if (nConstraints_ > 0)
  {
    if (printLevel > 0)
      printf("RSConstraints: number of filters = %d\n", nConstraints_);

    psuadeIO->getParameter("ana_rsfilterlbounds", flbPtr);
    vecLBounds_.load(nConstraints_, flbPtr.dbleArray_);
    psuadeIO->getParameter("ana_rsfilterubounds", fubPtr);
    vecUBounds_.load(nConstraints_, fubPtr.dbleArray_);

    psuadeIO->getParameter("ana_rsfilterdatafile", fdfPtr);
    filterDataFiles = fdfPtr.strArray_;
    psuadeIO->getParameter("ana_rsfilterindexfile", fifPtr);
    filterIndexFiles = fifPtr.strArray_;
    for (ii = 0; ii < nConstraints_; ii++)
    {
      if (strcmp(filterDataFiles[ii], "NULL"))
      {
        fp = fopen(filterDataFiles[ii], "r");
        if (fp == NULL)
        {
          printf("RSConstraint ERROR: filter data file %s not found.\n",
                 filterDataFiles[ii]);
          exit(1);
        }
        else fclose(fp);
      }
    }

    constraintFAs_ = new FuncApproxFilter*[nConstraints_];
    constraintNInputs_ = new int[nConstraints_];
    constraintInputIndices_ = new int*[nConstraints_];
    constraintInputValues_ = new double*[nConstraints_];
    for (ii = 0; ii < nConstraints_; ii++)
    {
      if (printLevel > 0)
        printf("RSConstraints: initializing filter %d\n",ii+1);
      constraintFAs_[ii] = NULL;
      constraintNInputs_[ii] = 0;
      constraintInputIndices_[ii] = NULL;
      constraintInputValues_[ii] = NULL;
      if (strcmp(filterDataFiles[ii], "NULL"))
      {
        pIO = new PsuadeData();
        status = pIO->readPsuadeFile(filterDataFiles[ii]);
        if (status != 0) 
        {
          printf("RSConstraints ERROR: cannot read filter data file %s.\n",
                 filterDataFiles[ii]);
          exit(1);
        }
        if (printLevel > 0)
          printf("RSConstraints: filter data file = %s\n",
                 filterDataFiles[ii]);
        pIO->getParameter("input_ninputs", pPtr);
        nInputsChk = pPtr.intData_;
        delete pIO;

        if (!strcmp(filterIndexFiles[ii],"NULL"))
        {
          printf("RSConstraints ERROR: filter index file not given.\n");
          exit(1);
        }
        if (printLevel > 0)
          printf("RSConstraints: filter %d has index file = %s\n",
                 ii+1, filterIndexFiles[ii]);
        fp = fopen(filterIndexFiles[ii],"r");
        if (fp != NULL)
        {
          fscanf(fp, "%d", &constraintNInputs_[ii]);
          if (constraintNInputs_[ii] != nInputsChk)
          {
            printf("RSConstraints: filter %d must have %d(%d) inputs\n",
                   ii+1, nInputsChk, constraintNInputs_[ii]);
            exit(1);
          }
          constraintInputIndices_[ii] = new int[constraintNInputs_[ii]];
          constraintInputValues_[ii] = new double[constraintNInputs_[ii]];
          for (jj = 0; jj < constraintNInputs_[ii]; jj++)
          {
            fscanf(fp,"%d %lg",&constraintInputIndices_[ii][jj],
                   &(constraintInputValues_[ii][jj]));
            printf("RSConstraints: filter %d has input %d = %d\n",
                   ii+1, jj+1, constraintInputIndices_[ii][jj]);
            if (constraintInputIndices_[ii][jj] >= 0)
              constraintInputIndices_[ii][jj]--;
            if (constraintInputIndices_[ii][jj] >= nInputs_ ||
                constraintInputIndices_[ii][jj] < -1)
            {
              printf("RSConstraints ERROR: invalid filter input %d.\n",
                     constraintInputIndices_[ii][jj]+1);
              printf("RSConstraints input format: ");
              printf("nInputs \n ");
              printf("<input (or 0)> <value (nominal val if 0)> \n ");
              printf("<input (or 0)> <value (nominal val if 0)> \n ");
              exit(1);
            }
          }
          vecSortArray.setLength(constraintNInputs_[ii]);
          for (jj = 0; jj < constraintNInputs_[ii]; jj++)
            vecSortArray[jj] = constraintInputIndices_[ii][jj];
          sortIntList(constraintNInputs_[ii], vecSortArray.getIVector());
          for (jj = 1; jj < constraintNInputs_[ii]; jj++)
          {
            if (vecSortArray[jj] >= 0 && vecSortArray[jj]==vecSortArray[jj-1])
            {
              printf("RSConstraints WARNING: filter has duplicate");
              printf(" input %d.\n", vecSortArray[jj]+1);
            }
          }
          fclose(fp);
        }
        else
        {
          printf("RSConstraints ERROR : filter index file %s not found.\n",
                 filterIndexFiles[ii]);
          exit(1);
        }
        printf("RSConstraints INFO: creating filter response surface %d (%d)\n",
               ii+1, nConstraints_);
        constraintFAs_[ii] = new FuncApproxFilter(filterDataFiles[ii]);
        constraintFAs_[ii]->setYBounds(flbPtr.dbleArray_[ii],
                                       fubPtr.dbleArray_[ii]);
      }
    }
  }
  return nConstraints_;
}

// ************************************************************************
// evaluate
// ------------------------------------------------------------------------
double RSConstraints::evaluate(double *sampleInputs, double sampleOutput, 
                               int &flag)
{
  int    ii, jj, status=1, numFailed;
  double dtemp=0.0;
  psVector vecSamplePt;

  if (nConstraints_ <= 0) 
  {
    flag = 1;
    return 0.0;
  }
  numFailed = 0;
  for (ii = 0; ii < nConstraints_; ii++)
  {
    if (constraintFAs_[ii] != NULL)
    {
      vecSamplePt.setLength(constraintNInputs_[ii]);
      for (jj = 0; jj < constraintNInputs_[ii]; jj++)
      {
        if (constraintInputIndices_[ii][jj] == -1)
          vecSamplePt[jj] = constraintInputValues_[ii][jj];
        if (constraintInputIndices_[ii][jj] != -1)
          vecSamplePt[jj]=sampleInputs[constraintInputIndices_[ii][jj]];
      }
      dtemp = constraintFAs_[ii]->evaluatePoint(vecSamplePt.getDVector(),status);
      if (status == 0) numFailed++;
    }
    else
    {
      if (sampleOutput < vecLBounds_[ii] || sampleOutput > vecUBounds_[ii])
        numFailed++;
    }
  }
  flag = 1;
  //**/ set intersection 
  if (psConfig_.RSConstraintSetOp_ == 0 && numFailed > 0) flag = 0;
  //**/ set union 
  if (psConfig_.RSConstraintSetOp_ == 1 && numFailed == nConstraints_) 
    flag = 0;
  return dtemp;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
RSConstraints& RSConstraints::operator=(const RSConstraints &)
{
  printf("RSConstraints operator= ERROR: operation not allowed.\n");
  exit(1);
  return (*this);
}

