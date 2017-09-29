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
// Definitions for the communication manager
// AUTHOR : CHARLES TONG
// DATE   : 2003
// ************************************************************************
#ifndef __COMMMANAGERH__
#define __COMMMANAGERH__

#include <stdio.h>
#include "sysdef.h"
#include "BaseComm.h"

class CommManager 
{
   int      numProcs_;
   BaseComm *gComm_;

public:

   CommManager(int, void **);
   ~CommManager();

   inline int getNumProcs() {return numProcs_;}
   inline int getPID()      {if (numProcs_==1) return 0; 
                             else return gComm_->getPID();}

   inline void synchronize() 
               {if (numProcs_ != 1) gComm_->synchronize();}

   inline int  send(void *msg,int leng,int dtyp,int msgid,int dest)
               {if (gComm_ != NULL) 
                     return gComm_->send(msg,leng,dtyp,msgid,dest);
                else return -1;}

   inline int  recv(void *msg,int leng,int dtyp,int msgid,int src)
               {if (gComm_ != NULL)
                     return gComm_->recv(msg,leng,dtyp,msgid,src);
                else return -1;}

   inline int  iRecv(void *msg,int leng,int dtyp,int msgid,int src)
               {if (gComm_ != NULL) 
                     return gComm_->iRecv(msg,leng,dtyp,msgid,src);
                else return -1;}

   inline int  disableIrecv(int proc)
               {if (gComm_ != NULL) return gComm_->disableIrecv(proc);
                else                return -1;}

   inline int  iProbe(int src, int msgID)
               {if (gComm_ != NULL) return gComm_->iProbe(src, msgID);
                else                return -1;}

   inline void wait(int handler)
               {if (gComm_ != NULL) gComm_->wait(handler);}

   inline int  waitAny()
               {if (gComm_ != NULL) return gComm_->waitAny();
                else                return -1;}

   inline void bcast(void *msg, int leng, int dtyp, int src)
               {if (gComm_ != NULL) gComm_->bcast(msg,leng,dtyp,src);}

   inline void allReduce(void *msg, int leng, int dtyp, char op)
               {if (gComm_ != NULL) gComm_->allReduce(msg,leng,dtyp,op);}

   inline void shutdown() {if (gComm_ != NULL) gComm_->shutdown();}
};

#endif // __COMMMANAGERH__

