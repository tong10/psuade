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
// Functions for the communication manager class
// AUTHOR : CHARLES TONG
// DATE   : 2003
// ************************************************************************
                                                                                
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "CommManager.h"

// ************************************************************************
// select communication protocol
// ------------------------------------------------------------------------
#if defined(HAVE_MPICH) 
#include "CommMPICH.h"
#endif

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
CommManager::CommManager(int argc, void **argv)
{
#if defined(HAVE_MPICH)
  int mypid;

  // -----------------------------------------------------------------
  // If MPICH communication protocol is requested
  // -----------------------------------------------------------------
  gComm_    = new CommMPICH(argc,argv);
  numProcs_ = gComm_->getNumProcs();
  mypid     = gComm_->getPID();
#if 0
  if (mypid == 0)
  {
     printf("CommManager : MPICH interface activated.\n");
     printf("CommManager : number of processors = %d\n",numProcs_);
  }
#endif

#else

  // -----------------------------------------------------------------
  // if no communication protocol is defined, use single processor 
  // -----------------------------------------------------------------
  gComm_    = NULL;
  numProcs_ = 1;
#if 0
  printf("CommManager : Single processor configuration.\n");
#endif

#endif
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
CommManager::~CommManager()
{
   int mypid=0;

   if (numProcs_ > 1) 
   {
      mypid = gComm_->getPID();
      delete gComm_;
   }
#if 0
   printf("Proc %4d: CommManager shutting down...\n", mypid);
#endif
}

