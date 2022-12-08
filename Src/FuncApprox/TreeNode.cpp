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
// Functions for the class TreeNode
// AUTHOR : CHARLES TONG
// DATE   : 2010
// ************************************************************************
#include <stdlib.h>
#include "TreeNode.h"

// ************************************************************************
// Constructor for object class SumOfTrees
// ------------------------------------------------------------------------
TreeNode::TreeNode()
{
  leftNode_   = NULL;
  rightNode_  = NULL;
  cutPoint_   = 0.0;
  whichInput_ = -1;
  nodeValue_  = 0.0;
  nodeStdev_  = 0.0;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
TreeNode::~TreeNode()
{
  if (leftNode_  != NULL) delete leftNode_;
  if (rightNode_ != NULL) delete rightNode_;
}

