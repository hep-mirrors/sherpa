#include "Amplitude_Manipulator.H"
#include "Run_Parameter.H"
#include "Message.H"

using namespace AMEGIC;
using namespace AORGTOOLS;
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
    if (fl[i].isfermion()) fermnumber++;
  }

  int* perm           = new int[fermnumber];
  Single_Amplitude* f = first_amp;
  
  Point* p;

  int count = 1;

  while (f) { 
    count++;
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
  int depth = 2*N-3;

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
      SetFermionNumberFlow(pb,pe);
      sign *= SetPropOrientation(pb,pe);
      f->AddSpinorDirection(pb->number,pe->number);
      pnumb+=2;
    }
  } 
  while(pnext);  
  
}

Point* Amplitude_Manipulator::FindNext(Point* p)
{
  if (p==0) return 0;
  if (p->fl.isfermion() && p->m==0) return p;
  
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
	!pbegin->fl.isanti() && !pend->fl.isanti()) {
      if (pbegin->number<pend->number) {
	Point* h = pbegin;
	pbegin   = pend;
	pend     = h;
	return;
      }
      return;
    }
    if (b[pbegin->number]==-1  && b[pend->number]==-1 && 
	pbegin->fl.isanti() && pend->fl.isanti()) {
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
	pbegin->fl.Majorana() && !pend->fl.isanti()) {
      return;
    }
    if (b[pbegin->number]==-1 && b[pend->number]==1 && 
	pbegin->fl.Majorana() && pend->fl.isanti()) {
      return;
    }
    if (b[pbegin->number]==1 && b[pend->number]==-1 && 
	pbegin->fl.isanti()  && pend->fl.Majorana()) {
      Point* h = pbegin;
      pbegin   = pend;
      pend     = h;
      return;
    }
    //normal cases
    if (!pbegin->fl.isanti() && b[pbegin->number]==-1) {
      Point* h = pbegin;
      pbegin   = pend;
      pend     = h;
      return;
    }
    if (pend->fl.isanti() && b[pend->number]==-1) {
      Point* h = pbegin;
      pbegin   = pend;
      pend     = h;
      return;
    }    
    return;
  }     
  //final fermion line
  if (!fl[pbegin->number].isanti() && !fl[pbegin->number].Majorana()) return;
  if (fl[pend->number].isanti()    && !fl[pend->number].Majorana()) return;
  if ( (fl[pbegin->number].isanti()  && !fl[pbegin->number].Majorana()) ||
       (!fl[pend->number].isanti()   && !fl[pend->number].Majorana()) ) {
    Point* h = pbegin;
    pbegin   = pend;
    pend     = h;
    return;
  }
  
  if (fl[pbegin->number].Majorana() && fl[pend->number].Majorana()) {
    if (pbegin->number<pend->number) {
      Point* h = pbegin;
      pbegin   = pend;
      pend     = h;
      return;
    }
    return;
  }

  cerr<<"Error in Amplitude_Manipulator::GetFermionLine()"<<endl;
  return;
}

Point* Amplitude_Manipulator::ForwardLine(Point* p)
{
  p->m = 1;
  if (p->left==0) return p;
  
  if (p->left->fl.isfermion())  return ForwardLine(p->left);
  if (p->middle) {
    if (p->middle->fl.isfermion())  return ForwardLine(p->middle);
  }


  if (p->right->fl.isfermion()) return ForwardLine(p->right);

  cerr<<"Dead fermion line!!!"<<endl;
}

Point* Amplitude_Manipulator::BackwardLine(Point* p)
{  
  p->m = 1;
  if (p->prev==0) return p;  

  if (p->prev->fl.isfermion()) return BackwardLine(p->prev);
  if (p->prev->left==p)        return ForwardLine(p->prev->right);
  if (p->prev->middle==p)      return ForwardLine(p->prev->middle);
  
  if (p->prev->right==p)       return ForwardLine(p->prev->left);

  cerr<<"Dead fermion line!!!"<<endl;
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
      AORGTOOLS::msg.Debugging()<<p->number<<" is vbar"<<endl;
    }
  }
  if (p->left==0) {
    if (b[p->number]==-1) {
      // ----<---O Orientation
      AORGTOOLS::msg.Debugging()<<p->number<<" is u"<<endl;
    }
    if (b[p->number]==1) {
      // O----<--- Orientation
      AORGTOOLS::msg.Debugging()<<p->number<<" is v"<<endl;
    }
    return;
  }

  int minus = 1;

  if (p->number>99 && p->m==1) {
    // ====>===== Fermion number flow
    // ----<----- Orientation
    // ---->----- Momentum Flow
    minus = -1;
  }
  if (p->number>99 && p->m==-1) {
    // ====<===== Fermion number flow
    // ----<----- Orientation
    // ---->----- Momentum Flow

    minus = -1;
  }

  if (p->m==1) {
    // ====>===== Fermion number flow
    // ----<----- Orientation
    // ---->----- Momentum Flow

    //Gamma'
    int ferm = 0;
    int vect = 0;
    int majo = 0;
    
    if (p->fl.isfermion()) ferm++;
    if (p->fl.isvector())  vect++;
    if (p->fl.Majorana())  majo++;
    if (p->left->fl.isfermion()) ferm++;
    if (p->left->fl.isvector())  vect++;
    if (p->left->fl.Majorana())  majo++;
    if (p->right->fl.isfermion()) ferm++;
    if (p->right->fl.isvector())  vect++;
    if (p->right->fl.Majorana())  majo++;

    if (vect==1 && ferm==2 && majo!=2) {
      Complex h = p->cpl[0];
      p->cpl[0] = -p->cpl[1];
      p->cpl[1] = -h;
    }    
  }

  if (minus==-1) {
    sign *= -1;
    //if (!p->fl.isanti()) p->fl = p->fl.bar();
    AORGTOOLS::msg.Debugging()<<"FL Flavour(-1) opposite to spin flow: "
			      <<p->fl<<";"<<p->t<<endl;
  }
  else {
    if (p->number>99) {
      //if (p->fl.isanti()) p->fl = p->fl.bar();
      AORGTOOLS::msg.Debugging()<<"FL Flavour in spin flow: "<<p->fl<<";"<<p->t<<endl;
    }
  }

  if (p->left->fl.isfermion())     {ForwardLineOrientation(p->left,sign); return;}
  if (p->middle) {
    if (p->middle->fl.isfermion()) {ForwardLineOrientation(p->middle,sign); return;}
  }

  if (p->right->fl.isfermion())    {ForwardLineOrientation(p->right,sign);return;}

  cerr<<"Dead fermion line!!!"<<endl;
}

void Amplitude_Manipulator::BackwardLineOrientation(Point* p,int& sign)
{  
  if (p->left==0) {  
    if (b[p->number]==-1) {
      // ----<---O Orientation
      AORGTOOLS::msg.Debugging()<<p->number<<" is vbar"<<endl;
    }
    if (b[p->number]==1) {
      // O----<--- Orientation
      AORGTOOLS::msg.Debugging()<<p->number<<" is ubar"<<endl;
    }    
  }
  if (p->prev==0) {  
    if (b[p->number]==-1) {
      // ----<---O Orientation
      AORGTOOLS::msg.Debugging()<<p->number<<" is u"<<endl;
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
    
    if ((p->prev)->fl.isfermion()) ferm++;
    if ((p->prev)->fl.isvector())  vect++;
    if ((p->prev)->fl.Majorana())  majo++;
    if ((p->prev)->left->fl.isfermion()) ferm++;
    if ((p->prev)->left->fl.isvector())  vect++;
    if ((p->prev)->left->fl.Majorana())  majo++;
    if ((p->prev)->right->fl.isfermion()) ferm++;
    if ((p->prev)->right->fl.isvector())  vect++;
    if ((p->prev)->right->fl.Majorana())  majo++;

    if (vect==1 && ferm==2 && majo!=2) {
      Complex h         = (p->prev)->cpl[0];
      (p->prev)->cpl[0] = -(p->prev)->cpl[1];
      (p->prev)->cpl[1] = -h;
    }
  }

  if (minus==-1) {
    sign *= -1;
    //if (!p->fl.isanti()) p->fl = p->fl.bar();
    AORGTOOLS::msg.Debugging()<<"BL Flavour(-1) opposite to spin flow: "
			      <<p->fl<<";"<<p->t<<endl;
  }
  else {
    if (p->number>99) {
      //if (p->fl.isanti()) p->fl = p->fl.bar();
      AORGTOOLS::msg.Debugging()<<"BL Flavour in spin flow             : "
				<<p->fl<<";"<<p->t<<endl;
    }
  }

  if (p->prev->fl.isfermion()) { BackwardLineOrientation(p->prev,sign);        return; }
  if (p->prev->left==p)        { ForwardLineOrientation(p->prev->right,sign);  return; }
  if (p->prev->middle==p)      { ForwardLineOrientation(p->prev->middle,sign); return; }
  if (p->prev->right==p)       { ForwardLineOrientation(p->prev->left,sign);   return; }

  cerr<<"Dead fermion line!!!"<<endl;
}

void Amplitude_Manipulator::SetFermionNumberFlow(Point* pb,Point* pe)
{
  int okay = 0;
  if (okay==0 && (b[pb->number]==-1 || (b[pe->number]==-1))) {
    //Initial Line
    //special cases
    if (b[pb->number]==-1 && b[pe->number]==1 &&
	pb->fl.Majorana() && pe->fl.isanti()) okay = 2; 
    if (b[pe->number]==-1 && b[pb->number]==1 &&
	pe->fl.Majorana() && pb->fl.isanti() && okay==0)  okay = 1; 
    if (b[pe->number]==-1 && b[pb->number]==-1 &&
	pb->fl.Majorana() && !pe->fl.isanti() && okay==0) okay = 2;    
    if (b[pb->number]==-1 && b[pe->number]==-1 &&
	pe->fl.Majorana() && !pb->fl.isanti() && okay==0) okay = 1;    
    //normal cases
    if (b[pb->number]==-1 && pb->fl.isanti()  && okay==0) okay = 2;
    if (b[pe->number]==-1 && !pe->fl.isanti() && okay==0) okay = 2;
  }
  
  if (okay==0 && b[pb->number]==1 && b[pe->number]==1) {
    if (!pb->fl.isanti())            okay = 2;
    if (pe->fl.isanti() && okay==0)  okay = 2;
  }

  // X0 -------O--->---O------X0
  //             e+/e- ?????????????

  if (okay==2) {
    Point* h = pb;
    pb       = pe;
    pe       = h;
  }

  int majoflag = 0;

  if (!pb->fl.Majorana() && !pe->fl.Majorana()) {
    if (b[pb->number]==-1 && b[pe->number]==-1 &&
	(  (!pb->fl.isanti() && !pe->fl.isanti()) ||
	   (pb->fl.isanti()  && pe->fl.isanti()) ) ) majoflag = 1;
    if (b[pb->number]==-1 && b[pe->number]==1  &&
	(  (pb->fl.isanti() && !pe->fl.isanti()) ||
	   (pe->fl.isanti() && !pb->fl.isanti()) ) ) majoflag = 1;
    if (b[pe->number]==-1 && b[pb->number]==1  &&
	(  (pb->fl.isanti() && !pe->fl.isanti()) ||
	   (pe->fl.isanti() && !pb->fl.isanti()) ) ) majoflag = 1;
  }
  
  if (majoflag) {
    AORGTOOLS::msg.Debugging()<<"Fermion Flow(Majo): "
			      <<pb->number<<" <-> "<<pe->number<<endl;
    if (!pb->fl.isanti() && b[pb->number]==-1) majoflag=1;
    if (pb->fl.isanti()  && b[pb->number]==-1) majoflag=-1;
    if (!pb->fl.isanti() && b[pb->number]==1)  majoflag=-1;
    if (pb->fl.isanti()  && b[pb->number]==1)  majoflag=1;

    if (pb->prev==0) SetForwardFNFlow(pb,majoflag);
    else SetBackwardFNFlow(pb,majoflag);

    if (!pe->fl.isanti() && b[pe->number]==-1) majoflag=1;
    if (pe->fl.isanti()  && b[pe->number]==-1) majoflag=-1;
    if (!pe->fl.isanti() && b[pe->number]==1)  majoflag=-1;
    if (pe->fl.isanti()  && b[pe->number]==1)  majoflag=1;
    
    if (pe->prev==0) SetForwardFNFlow(pe,majoflag);
    else SetBackwardFNFlow(pe,majoflag);
  }
  else {   
    AORGTOOLS::msg.Debugging()<<"Fermion Flow: "
			      <<pb->number<<" -> "<<pe->number<<endl;
    if (pb->prev==0) SetForwardFNFlow(pb,0);
    else SetBackwardFNFlow(pb,0);
  }
}

void Amplitude_Manipulator::SetForwardFNFlow(Point* p,int majoflag)
{
  if (p->left==0) return;
  
  //nothing......
  //p->m = 1;
  
  if (majoflag==-1) p->m = -1;
  
  if (p->fl.Majorana() && majoflag) return;
  
  if (p->left->fl.isfermion())     { SetForwardFNFlow(p->left,majoflag);   return; }
  if (p->middle) {
    if (p->middle->fl.isfermion()) { SetForwardFNFlow(p->middle,majoflag); return; }
  }
  if (p->right->fl.isfermion())    { SetForwardFNFlow(p->right,majoflag);  return; }

  cerr<<"Dead fermion line!!!"<<endl;
}

void Amplitude_Manipulator::SetBackwardFNFlow(Point* p,int majoflag)
{  
  if (p->prev==0) return;  

  if (p->fl.Majorana() && majoflag) return;
  if (majoflag==-1) p->m = 1;
  else p->m = -1;

  if (p->prev->fl.isfermion()) { SetBackwardFNFlow(p->prev,majoflag);        return; }
  if (p->prev->left==p)        { SetForwardFNFlow(p->prev->right,majoflag);  return; }
  if (p->prev->middle==p)      { SetForwardFNFlow(p->prev->middle,majoflag); return; }
  if (p->prev->right==p)       { SetForwardFNFlow(p->prev->left,majoflag);   return; }

  cerr<<"Dead fermion line!!!"<<endl;
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

