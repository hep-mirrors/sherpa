#include "AMEGIC++/Amplitude/Amplitude_Manipulator.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;

void Amplitude_Manipulator::SetPrev(Point* p)
{
  if (p->left==0) return;
  
  p->left->prev  = p;
  p->right->prev = p;
  if (p->middle) p->middle->prev = p;

  SetPrev(p->left);
  if (p->middle) SetPrev(p->middle);
  SetPrev(p->right);
}

void Amplitude_Manipulator::FixSign(Single_Amplitude* first_amp)
{
  int fermnumber = 0;
  for (short int i=0;i<N;i++) {
    if (fl[i].IsFermion()) fermnumber++;
  }

  int* perm           = new int[fermnumber];
  Single_Amplitude* f = first_amp;
  
  Point* p;

  while (f) { 
    p = f->GetPointlist();
    p[0].prev = 0;
    SetPrev(p);
    f->sign = 1;
    GetPerm(perm,f,f->sign);
    f->sign *= Permutation(perm,fermnumber);
    f = f->Next;
  }
  
  delete[] perm;
}

void Amplitude_Manipulator::GetPerm(int* perm,Single_Amplitude* f,int& sign)
{
  Point* p = f->GetPointlist();
  int depth = 2*N-3+dec;

  for (short int i=0;i<depth;i++) p[i].m = 0;
  
  Point* pnext;

  int pnumb = 0;

  do {
    pnext = FindNext(p);
    
    if (pnext) {
      Point* pb;
      Point* pe;
      GetFermionLine(pnext,pb,pe);
      perm[pnumb]   = pb->number;
      perm[pnumb+1] = pe->number;
      int swap = SetFermionNumberFlow(pb,pe);
      sign *= SetPropOrientation(pb,pe);
      if (swap) f->AddSpinorDirection(pe->number,pb->number);      
           else f->AddSpinorDirection(pb->number,pe->number);      
      pnumb+=2;
    }
  } 
  while(pnext);  
}

Point* Amplitude_Manipulator::FindNext(Point* p)
{
  if (p==0) return 0;
  if (p->fl.IsFermion() && p->m==0) return p;
  
  Point* ptmp = FindNext(p->left);
  if (ptmp!=0) return ptmp;
  if (p->middle) {
    ptmp = FindNext(p->middle);
    if (ptmp!=0) return ptmp;
  }
  return FindNext(p->right);
}
 
void Amplitude_Manipulator::GetFermionLine(Point* pcurr,Point*& pbegin,Point*& pend)
{
  pbegin = pend = 0;

  pbegin = BackwardLine(pcurr);
  pend   = ForwardLine(pcurr);
  
  if (b[pbegin->number]==-1 ||
      b[pend->number]  ==-1) {
    //pure Majorana cases
    if (pbegin->fl.Majorana() && pend->fl.Majorana()) {
      if (pbegin->number<pend->number) {
	Point* h = pbegin;
	pbegin   = pend;
	pend     = h;
	return;
      }
      return;
    }
    //Majorana propagator case
    if (b[pbegin->number]==-1  && b[pend->number]==-1 && 
	!pbegin->fl.IsAnti() && !pend->fl.IsAnti()) {
      if (pbegin->number<pend->number) {
	Point* h = pbegin;
	pbegin   = pend;
	pend     = h;
	return;
      }
      return;
    }
    if (b[pbegin->number]==-1  && b[pend->number]==-1 && 
	pbegin->fl.IsAnti() && pend->fl.IsAnti()) {
      if (pbegin->number>pend->number) {
	Point* h = pbegin;
	pbegin   = pend;
	pend     = h;
	return;
      }
      return;
    }
    //special cases
    if (b[pbegin->number]==-1 && b[pend->number]==-1 && 
	pbegin->fl.Majorana() && !pend->fl.IsAnti()) return;
    if (b[pbegin->number]==-1 && b[pend->number]==1 && 
	pbegin->fl.Majorana() && pend->fl.IsAnti())  return;
    if (b[pbegin->number]==1 && b[pend->number]==-1 && 
	pbegin->fl.IsAnti()  && pend->fl.Majorana()) {
      Point* h = pbegin;
      pbegin   = pend;
      pend     = h;
      return;
    }
    //normal cases
    if (!pbegin->fl.IsAnti() && b[pbegin->number]==-1) {
      Point* h = pbegin;
      pbegin   = pend;
      pend     = h;
      return;
    }
    if (pend->fl.IsAnti() && b[pend->number]==-1) {
      Point* h = pbegin;
      pbegin   = pend;
      pend     = h;
      return;
    }
    return;
  }     

  if (fl[pbegin->number].IsAnti() && fl[pbegin->number].IsChargino() && fl[pend->number].Majorana()) {
    Point* h = pbegin;
    pbegin   = pend;
    pend     = h;
    return;
  }
  
  if (!fl[pend->number].IsAnti() && fl[pend->number].IsChargino() && fl[pbegin->number].Majorana()) {
    Point* h = pbegin;
    pbegin   = pend;
    pend     = h;
    return;
    }
  
  if (fl[pend->number].IsAnti() && fl[pbegin->number].Majorana() && 
      fl[pend->number].IsChargino()) return;
  
  if (!fl[pbegin->number].IsAnti() && fl[pbegin->number].IsChargino() &&
      fl[pend->number].IsAnti() && fl[pend->number].IsChargino()) return;
  
  if (fl[pbegin->number].IsAnti() && fl[pbegin->number].IsChargino() &&
      !fl[pend->number].IsAnti() && fl[pend->number].IsChargino()) {
    Point* h = pbegin;
    pbegin   = pend;
    pend     = h;
    return;
  }
  
  //standard final fermion line
  if (!fl[pbegin->number].IsAnti() && !fl[pbegin->number].Majorana()) return;
  if (fl[pend->number].IsAnti() && !fl[pend->number].Majorana()) return;
  if ( (fl[pbegin->number].IsAnti()  && !fl[pbegin->number].Majorana()) ||
       (!fl[pend->number].IsAnti()   && !fl[pend->number].Majorana() && 
	!fl[pend->number].IsChargino()) ) {
    Point* h = pbegin;
    pbegin   = pend;
    pend     = h;
    return;
  }
 
  if (fl[pbegin->number].Majorana() && fl[pend->number].Majorana()) return;
  
  msg_Error()<<"ERROR in Amplitude_Manipulator::GetFermionLine(). Continue run."<<endl;
  return;
}

Point* Amplitude_Manipulator::ForwardLine(Point* p)
{
  p->m = 1;
  
  if (p->left==0) return p;
  
  if (p->left->fl.IsFermion())  return ForwardLine(p->left);
  if (p->middle) {
    if (p->middle->fl.IsFermion())  return ForwardLine(p->middle);
  }


  if (p->right->fl.IsFermion()) return ForwardLine(p->right);

  msg_Error()<<"ERROR in Amplitude_Manipulator::ForwardLine :"<<std::endl
	     <<"   Dead fermion line in Amplitude_Manipulator::ForwardLine. Continue run."<<endl;
  return 0;
}

Point* Amplitude_Manipulator::BackwardLine(Point* p)
{  
  p->m = 1;

  if (p->prev==0) return p;  

  if (p->prev->fl.IsFermion()) return BackwardLine(p->prev);
  if (p->prev->left==p){
    if(p->prev->right->fl.IsFermion()) return ForwardLine(p->prev->right);
    return ForwardLine(p->prev->middle);
  }
  if (p->prev->middle==p){
    if(p->prev->right->fl.IsFermion()) return ForwardLine(p->prev->right);
    return ForwardLine(p->prev->left);
  }
  if (p->prev->right==p){
    if(p->prev->left->fl.IsFermion()) return ForwardLine(p->prev->left);
    return ForwardLine(p->prev->middle);
  }

  msg_Error()<<"ERROR in Amplitude_Manipulator::BackwardLine :"<<std::endl
	     <<"   Dead fermion line in Amplitude_Manipulator::BackwardLine. Continue run."<<endl;
  return 0;
}

int Amplitude_Manipulator::SetPropOrientation(Point* pb,Point* pe)
{
  int sign = 1;
  
  if (pb->prev==0) ForwardLineOrientation(pb,sign);
  else BackwardLineOrientation(pb,sign);
  return sign;
}

void Amplitude_Manipulator::ForwardLineOrientation(Point* p,int& sign)
{
  if (p->prev==0) {
    if (b[p->number]==-1) {
      // ----<---O Orientation
      // msg_Debugging()<<p->number<<" is vbar"<<endl;
    }
  }
  if (p->left==0) {
    if (b[p->number]==-1) {
      // ----<---O Orientation
      // msg_Debugging()<<p->number<<" is u"<<endl;
    }
    if (b[p->number]==1) {
      // O----<--- Orientation
      // msg_Debugging()<<p->number<<" is v"<<endl;
    }
    return;
  }

  int minus = 1;

  if (p->number>99 && p->m==1 && !p->fl.Majorana()) {
    // ====>===== Fermion number flow
    // ----<----- Orientation
    // ---->----- Momentum Flow
    minus = -1;
  }
  if (p->number>99 && p->m==-1 && !p->fl.Majorana()) {
    // ====<===== Fermion number flow
    // ----<----- Orientation
    // ---->----- Momentum Flow
    minus = -1;
  }

  if (p->m==-1) {
    // ====>===== Fermion number flow
    // ----<----- Orientation
    // ---->----- Momentum Flow

    //Gamma'
    int ferm = 0;
    int vect = 0;
    int majo = 0;
    
    if (p->fl.IsFermion()) ferm++;
    if (p->fl.IsVector())  vect++;
    if (p->fl.Majorana())  majo++;
    if (p->left->fl.IsFermion()) ferm++;
    if (p->left->fl.IsVector())  vect++;
    if (p->left->fl.Majorana())  majo++;
    if (p->right->fl.IsFermion()) ferm++;
    if (p->right->fl.IsVector())  vect++;
    if (p->right->fl.Majorana())  majo++;

    if (vect==1 && ferm==2 && majo!=2) {
      Complex h = p->cpl[0];
      p->cpl[0] = -p->cpl[1];
      p->cpl[1] = -h;
    }   
  }

  if (minus==-1) {
    sign *= -1;
    //msg_Debugging()<<"FL Flavour(-1) opposite to spin flow: "<<p->fl<<";"<<p->t<<endl;
  }
  else {
    if (p->number>99) {
      //msg_Debugging()<<"FL Flavour in spin flow: "<<p->fl<<";"<<p->t<<endl;
    }
  }

  if (p->left->fl.IsFermion())     {ForwardLineOrientation(p->left,sign); return;}
  if (p->middle) {
    if (p->middle->fl.IsFermion()) {ForwardLineOrientation(p->middle,sign); return;}
  }

  if (p->right->fl.IsFermion())    {ForwardLineOrientation(p->right,sign);return;}

  msg_Error()<<"ERROR in Amplitude_Manipulator::ForwardLineOrientation :"<<std::endl 
	     <<"   Dead fermion line. Continue run."<<endl;
}

void Amplitude_Manipulator::BackwardLineOrientation(Point* p,int& sign)
{  
  if (p->left==0) {  
    if (b[p->number]==-1) {
      // ----<---O Orientation
      // msg_Debugging()<<p->number<<" is vbar"<<endl;
    }
    if (b[p->number]==1) {
      // O----<--- Orientation
      // msg_Debugging()<<p->number<<" is ubar"<<endl;
    }    
  }
  if (p->prev==0) {  
    if (b[p->number]==-1) {
      // ----<---O Orientation
      // msg_Debugging()<<p->number<<" is u"<<endl;
    }
    return;
  }

  int minus = 1;

  if (p->number>99 && p->m==-1) {
    // ====>===== Fermion number flow
    // ----<----- Orientation
    // ----<----- Momentum Flow
    //Okay
  }
  if (p->number>99 && p->m==1) {
    // ====<===== Fermion number flow
    // ----<----- Orientation
    // ----<----- Momentum Flow
    //Spinorflow != Momentumflow
    //Okay
  }
  
  if (p->m==-1) {
    // ====>===== Fermion number flow
    // ----<----- Orientation
    // ----<----- Momentum Flow
    //Gamma'
    int ferm = 0;
    int vect = 0;
    int majo = 0;
    
    if ((p->prev)->fl.IsFermion()) ferm++;
    if ((p->prev)->fl.IsVector())  vect++;
    if ((p->prev)->fl.Majorana())  majo++;
    if ((p->prev)->left->fl.IsFermion()) ferm++;
    if ((p->prev)->left->fl.IsVector())  vect++;
    if ((p->prev)->left->fl.Majorana())  majo++;
    if ((p->prev)->right->fl.IsFermion()) ferm++;
    if ((p->prev)->right->fl.IsVector())  vect++;
    if ((p->prev)->right->fl.Majorana())  majo++;

    if (vect==1 && ferm==2 && majo!=2) {
      Complex h         = (p->prev)->cpl[0];
      (p->prev)->cpl[0] = -(p->prev)->cpl[1];
      (p->prev)->cpl[1] = -h;
      }
  }

  if (minus==-1) {
    sign *= -1;
    // msg_Debugging()<<"BL Flavour(-1) opposite to spin flow: "<<p->fl<<";"<<p->t<<endl;
  }
  else {
    if (p->number>99) {
      // msg_Debugging()<<"BL Flavour in spin flow             : "<<p->fl<<";"<<p->t<<endl;
    }
  }

  if (p->prev->fl.IsFermion()) { BackwardLineOrientation(p->prev,sign); return; }
  if (p->prev->left==p){
    if (p->prev->right->fl.IsFermion()) { ForwardLineOrientation(p->prev->right,sign); return; }
    ForwardLineOrientation(p->prev->middle,sign); return;
  }
  if (p->prev->middle==p){
    if(p->prev->right->fl.IsFermion()) { ForwardLineOrientation(p->prev->right,sign); return;}
    ForwardLineOrientation(p->prev->left,sign); return;
  }
  if (p->prev->right==p){
    if (p->prev->left->fl.IsFermion()) { ForwardLineOrientation(p->prev->left,sign); return; }
    ForwardLineOrientation(p->prev->middle,sign); return;
  }

  msg_Error()<<"ERROR in Amplitude_Manipulator::BackwardLineOrientation :"<<std::endl 
	     <<"   Dead fermion line. Continue run."<<endl;
}

int Amplitude_Manipulator::SetFermionNumberFlow(Point* pb,Point* pe)
{
  int okay = 0;
  if (okay==0 && (b[pb->number]==-1 || (b[pe->number]==-1))) {
    //Initial Line
    //special cases
    if (b[pb->number]==-1 && b[pe->number]==1 &&
	pb->fl.Majorana() && pe->fl.IsAnti()) okay = 2; 
    if (b[pe->number]==-1 && b[pb->number]==1 &&
	pe->fl.Majorana() && pb->fl.IsAnti() && okay==0)  okay = 1; 
    if (b[pe->number]==-1 && b[pb->number]==-1 &&
	pb->fl.Majorana() && !pe->fl.IsAnti() && okay==0) okay = 2;    
    if (b[pb->number]==-1 && b[pe->number]==-1 &&
	pe->fl.Majorana() && !pb->fl.IsAnti() && okay==0) okay = 1;    
    //normal cases
    if (b[pb->number]==-1 && pb->fl.IsAnti()  && okay==0) okay = 2;
    if (b[pe->number]==-1 && !pe->fl.IsAnti() && okay==0) okay = 2;
  }
     
  if (okay==0 && b[pb->number]==1 && b[pe->number]==1) {
    if (!pb->fl.IsAnti())            okay = 2;
    if (pe->fl.IsAnti() && okay==0)  okay = 2;
  }

  if (pb->fl.Majorana() && pe->fl.Majorana()) okay = 0; 
  
  if (okay==2) {
    Point* h = pb;
    pb       = pe;
    pe       = h;
  }

  int majoflag = 0;
  int chinum = 0;
  if (!pb->fl.Majorana() && !pe->fl.Majorana()) {
    if (pb->fl.IsChargino()) chinum++;
    if (pe->fl.IsChargino()) chinum++;
    if (chinum!=1) {
    if (b[pb->number]==-1 && b[pe->number]==-1 &&
	((!pb->fl.IsAnti() && !pe->fl.IsAnti()) ||
	 (pb->fl.IsAnti()  && pe->fl.IsAnti()))) majoflag = 1;
    if (b[pb->number]==-1 && b[pe->number]==1  &&
	((pb->fl.IsAnti() && !pe->fl.IsAnti()) ||
	 (pe->fl.IsAnti() && !pb->fl.IsAnti()))) majoflag = 1;
    if (b[pe->number]==-1 && b[pb->number]==1  &&
	((pb->fl.IsAnti() && !pe->fl.IsAnti()) ||
	 (pe->fl.IsAnti() && !pb->fl.IsAnti()))) majoflag = 1;
  
    }
    else {
    if (b[pb->number]==-1 && b[pe->number]==-1 &&
	((!pb->fl.IsAnti() && pe->fl.IsAnti()) ||
	 (pb->fl.IsAnti()  && !pe->fl.IsAnti()))) majoflag = 1;
    if (b[pb->number]==-1 && b[pe->number]==1  &&
	((pb->fl.IsAnti() && pe->fl.IsAnti()) ||
	 (!pe->fl.IsAnti() && !pb->fl.IsAnti()))) majoflag = 1;
    if (b[pe->number]==-1 && b[pb->number]==1  &&
	((pb->fl.IsAnti() && pe->fl.IsAnti()) ||
	 (!pe->fl.IsAnti() && !pb->fl.IsAnti()))) majoflag = 1;
    }
  }
  
  int fermflag = 0;
  
  //new
  if (pb->fl.Majorana() && pe->fl.Majorana()) fermflag = 2;
  
  if (majoflag) {
    if (!pb->fl.IsAnti() && b[pb->number]==-1) majoflag=1;
    if (pb->fl.IsAnti()  && b[pb->number]==-1) majoflag=-1;
    if (!pb->fl.IsAnti() && b[pb->number]==1)  majoflag=-1;
    if (pb->fl.IsAnti()  && b[pb->number]==1)  majoflag=1;

    if (pb->prev==0) SetForwardFNFlow(pb,majoflag,fermflag);
    else SetBackwardFNFlow(pb,majoflag,fermflag);

    if (!pe->fl.IsAnti() && b[pe->number]==-1) majoflag=1;
    if (pe->fl.IsAnti()  && b[pe->number]==-1) majoflag=-1;
    if (!pe->fl.IsAnti() && b[pe->number]==1)  majoflag=-1;
    if (pe->fl.IsAnti()  && b[pe->number]==1)  majoflag=1;
    
    if (pe->prev==0) SetForwardFNFlow(pe,majoflag,fermflag);
    else SetBackwardFNFlow(pe,majoflag,fermflag);
  }
  if (fermflag) {
    if (pb->prev==0) SetForwardFNFlow(pb,0,fermflag);
    else SetBackwardFNFlow(pb,0,fermflag);
    if (fermflag==1) return 1;
  }
  if (!fermflag && !majoflag) {   
    if (pb->prev==0) SetForwardFNFlow(pb,0,fermflag);
    else SetBackwardFNFlow(pb,0,fermflag);
  }
  return 0;
}

void Amplitude_Manipulator::SetForwardFNFlow(Point* p,int majoflag,int& fermflag)
{
  //nothing......
  //p->m = 1;
  
  if (majoflag==-1) p->m = -1;
  
  if (p->fl.Majorana() && majoflag) return;
  
  //new version
  if (fermflag) {
    if (fermflag==2) {
      if (!p->fl.Majorana()) {
	if (p->fl.IsAnti()) fermflag=-1;
	else fermflag=1;
	
	int done = 0;
	//set m's of all majos before
	if (p->prev->fl.IsFermion() && p->prev->fl.Majorana()) { SetMajoFlowBackward(p->prev,fermflag); done = 1;}
	if (p->prev->left==p && done==0) {
	  if (p->prev->right->fl.IsFermion() && p->prev->right->fl.Majorana()) { SetMajoFlowForward(p->prev->right,fermflag); done = 1; }
	  if (done==0) SetMajoFlowForward(p->prev->middle,fermflag); done = 1;
	}
	if (p->prev->middle==p && done==0){
	  if(p->prev->right->fl.IsFermion() && p->prev->right->fl.Majorana()) { SetMajoFlowForward(p->prev->right,fermflag); done = 1;}
	if (done==0) SetMajoFlowForward(p->prev->left,fermflag); done = 1;
	}
	if (p->prev->right==p && done==0){
	  if (p->prev->left->fl.IsFermion() && p->prev->left->fl.Majorana()) { SetMajoFlowForward(p->prev->left,fermflag); done = 1; }
	  if (done==0) SetMajoFlowForward(p->prev->middle,fermflag);
	}
      }
    }
    if (fermflag==-1) p->m = -1;
    if (fermflag==1) p->m = 1;
  }    

  //new
  if (p->left==0) return;

  if (p->left->fl.IsFermion())     { SetForwardFNFlow(p->left,majoflag,fermflag);   return; }
  if (p->middle) {
    if (p->middle->fl.IsFermion()) { SetForwardFNFlow(p->middle,majoflag,fermflag); return; }
  }
  if (p->right->fl.IsFermion())    { SetForwardFNFlow(p->right,majoflag,fermflag);  return; }

  msg_Error()<<"ERROR in Amplitude_Manipulator::SetForwardFNFlow : Dead fermion line, continue run."<<endl;
}

void Amplitude_Manipulator::SetMajoFlowForward(Point* p,int fermflag)
{
  if (p==0) return;

  if (fermflag==-1) p->m = 1;
  if (fermflag==1) p->m = -1;

  if (p->left==0) return;

  if (p->left->fl.IsFermion() && p->left->fl.Majorana()) { 
    SetMajoFlowForward(p->left,fermflag);   
    return;
  }
  if (p->middle) {
    if (p->middle->fl.IsFermion() && p->middle->fl.Majorana()) { 
      SetMajoFlowForward(p->middle,fermflag); return; 
    }
  }
  if (p->right->fl.IsFermion() && p->right->fl.Majorana()) { 
    SetMajoFlowForward(p->right,fermflag); return;
  }
}
void Amplitude_Manipulator::SetBackwardFNFlow(Point* p,int majoflag,int& fermflag)
{  
  if (p->fl.Majorana() && majoflag) return;
  if (majoflag==-1) p->m = 1;
  else p->m = -1;
  
  //new version
  if (fermflag) {
    if (fermflag==2) {
      if (!p->fl.Majorana()) {
	if (!p->fl.IsAnti()) fermflag=-1;
	else fermflag=1;
	
	int done = 0;
	//set m's of all majos before
	if (p->left->fl.IsFermion() && p->left->fl.Majorana()) { 
	  SetMajoFlowForward(p->left,fermflag);   
	  done = 1; 
	}
	if (p->middle) {
	  if (p->middle->fl.IsFermion() && p->middle->fl.Majorana() && done==0) { 
	    SetMajoFlowForward(p->middle,fermflag); done = 1; 
	}
	}
	if (p->right->fl.IsFermion() && p->right->fl.Majorana() && done==0) { 
	  SetMajoFlowForward(p->right,fermflag);
	}
      }
    }
    if (fermflag==-1) p->m = 1;
    if (fermflag==1) p->m = -1;
  }   
    
  //new
  if (p->prev==0) return;  

  if (p->prev->fl.IsFermion()) { SetBackwardFNFlow(p->prev,majoflag,fermflag); return; }
  if (p->prev->left==p){
    if (p->prev->right->fl.IsFermion()) { SetForwardFNFlow(p->prev->right,majoflag,fermflag); return; }
    SetForwardFNFlow(p->prev->middle,majoflag,fermflag); return;
  }
  if (p->prev->middle==p){
    if(p->prev->right->fl.IsFermion()) { SetForwardFNFlow(p->prev->right,majoflag,fermflag); return;}
    SetForwardFNFlow(p->prev->left,majoflag,fermflag); return;
  }
  if (p->prev->right==p){
    if (p->prev->left->fl.IsFermion()) { SetForwardFNFlow(p->prev->left,majoflag,fermflag); return; }
    SetForwardFNFlow(p->prev->middle,majoflag,fermflag); return;
  }
  msg_Error()<<"ERROR in Amplitude_Manipulator::SetBackwardFNFlow : Dead fermion line, continue run."<<endl;
}

void Amplitude_Manipulator::SetMajoFlowBackward(Point* p,int fermflag)
{
  if (p->prev==0) return;

  if (fermflag==-1) p->m = -1;
  if (fermflag==1) p->m = 1;
  
  if (p->prev->fl.IsFermion() && p->prev->fl.Majorana()) { SetMajoFlowBackward(p->prev,fermflag); return;}
  if (p->prev->left==p) {
    if (p->prev->right->fl.IsFermion() && p->prev->right->fl.Majorana()) { SetMajoFlowForward(p->prev->right,fermflag); return; }
    SetMajoFlowForward(p->prev->middle,fermflag); return;
  }
      if (p->prev->middle==p){
	if(p->prev->right->fl.IsFermion() && p->prev->right->fl.Majorana()) { SetMajoFlowForward(p->prev->right,fermflag); return;}
	SetMajoFlowForward(p->prev->left,fermflag); return;
      }
      if (p->prev->right==p){
	if (p->prev->left->fl.IsFermion() && p->prev->left->fl.Majorana()) { SetMajoFlowForward(p->prev->left,fermflag); return; }
	SetMajoFlowForward(p->prev->middle,fermflag); return;
      }
}

int Amplitude_Manipulator::Permutation(int* perm,int fermnumber)
{
  int steps = 0;
  for (short int i=0;i<fermnumber;i++) {
    for (short int j=i+1;j<fermnumber;j++) {
      if (perm[i]>perm[j]) {
	int h   = perm[i];
	perm[i] = perm[j];
	perm[j] = h;
	steps++;
      }
    }
  }

  if (!(steps%2)) return 1;
  return -1;
}

