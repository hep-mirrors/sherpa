//#include <iomanip>
#include "Amplitude_Generator.H"
#include "Amplitude_Manipulator.H"
#include "Run_Parameter.H"
#include "Message.H"
#include "Model.H"

using namespace AMEGIC;
using namespace AORGTOOLS;
using namespace APHYTOOLS;
using namespace std;

Amplitude_Generator::Amplitude_Generator(int _no,Flavour* _fl,int* _b,Topology* _top,
					 Basic_Sfuncs* _BS,String_Handler* _shand) 
  : N(_no), fl(_fl), b(_b), top(_top), BS(_BS), shand(_shand)
{
  single_top = top->Get(N-2);
  
  // 2 incoming
  prenum  = 1000;
  prea    = new Pre_Amplitude[prenum];
  short int i;
  for(i=0;i<prenum;i++) prea[i].p = new Point[single_top->depth];
}

Amplitude_Generator::~Amplitude_Generator() 
{
  for(short int i=0;i<prenum;i++) delete[] prea[i].p;
  delete[] prea;
}

void Amplitude_Generator::Set_End(Point* p,int* &perm,int& pnum)
{
  p->b     = 0;
  p->fl    = Flavour(kf::none);
  if ((p->left==0) && (p->right==0)) {
    p->number = *perm;
    p->fl = fl[*perm];
    p->b  = b[*perm];
    if (p->fl.IsBoson()) {
      p->Lorentz->type = lf::Pol;
      p->Lorentz->SetParticleArg(0);
    }
    else {
      p->Lorentz->type = lf::None;
      p->Lorentz->SetParticleArg();
    }

    perm++;
    return;
  }
  p->number = pnum;
  pnum++;
  Set_End(p->left,perm,pnum);
  Set_End(p->right,perm,pnum);
}

void Amplitude_Generator::Next_P(Point* p,Point* &hit)
{
  if (hit) return;
  if (p==0) return;
  if ((p->left!=0) && (p->right!=0)) {
    if ((p->left->fl==Flavour(kf::none)) || (p->right->fl==Flavour(kf::none))) {
      hit = p;
      return;
    }
  }
  Next_P(p->left,hit); 
  Next_P(p->right,hit);
}


void Amplitude_Generator::Print_P(Point* p)
{
  if (!(AORGTOOLS::rpa.gen.Debugging())) return;
  if ((p->left==0) && (p->right==0)) {
    AORGTOOLS::msg.Out()<<"EndPoint : "<<p->fl<<"("<<p->b<<")"<<endl;
    return;
  }
  AORGTOOLS::msg.Out()<<"left : ";
  Print_P(p->left);
  AORGTOOLS::msg.Out()<<"right : ";
  Print_P(p->right);
}

int Amplitude_Generator::Match_Vertex(Single_Vertex* v,Flavour* flav,Complex* cpl)
{
  if (flav[0] == v->in[0]) {
    short int hit = 1;
    if (flav[1] != Flavour(kf::none)) {if (flav[1] != v->in[1]) hit = 0;}
    else {flav[1] = v->in[1];}
    if (flav[2] != Flavour(kf::none)) {if (flav[2] != v->in[2]) hit = 0;}
    else {flav[2] = v->in[2];}
    if (hit==1) {
      for (short int j=0;j<4;j++) cpl[j] = v->cpl[j];
      return 1;
    }
  }
  return 0;
}

int Amplitude_Generator::Check_End(Point* p,Flavour infl,Vertex* v) 
{
  if (p==0) return 1;
  if (p->left==0) return 1;
  if (((p->left->fl)!=Flavour(kf::none)) && ((p->right->fl)!=Flavour(kf::none))) { 
    short int j,k,hit;
    Flavour flav[3];
    Complex cpl[4];
    for (j=0;j<4;j++) cpl[j] = Complex(0.,0.);

    for (j=0;j<v->Max_Number();j++) {
      if ((*v)[j]->on) { 
	hit = 0; 
	flav[0] = infl;
	flav[1] = p->left->fl;
	flav[2] = p->right->fl;

	if (infl.Majorana()) {
	  if (p->left->fl.IsFermion()) {
	    if (p->b*p->left->b==1) flav[1] = flav[1].Bar(); 
	    if (p->right->b == -1)  flav[2] = flav[2].Bar();
	  }
	  else {
	    if (p->b*p->right->b==1) flav[2] = flav[2].Bar(); 
	    if (p->left->b == -1)    flav[1] = flav[1].Bar();
	  }	  
	}
	else {
	  if (infl.IsBoson()) {
	    if (p->left->b == -1) flav[1] = flav[1].Bar();
	    if (p->right->b == -1) flav[2] = flav[2].Bar();
	  }
	  else {
	    if (infl.IsAnti()) {
	      if (p->b*p->left->b==1) flav[1] = flav[1].Bar(); 
	      if (p->right->b == -1)  flav[2] = flav[2].Bar();
	    }
	    else {
	      if (p->b*p->right->b==1) flav[2] = flav[2].Bar(); 
	      if (p->left->b == -1)    flav[1] = flav[1].Bar();
	    }
	  }
	}
	if (Match_Vertex((*v)[j],flav,cpl)) {
	  for (k=0;k<4;k++) p->cpl[k] = cpl[k];
	  p->v = (*v)[j];
	  *(p->Color)   = *((*v)[j]->Color);
	  *(p->Lorentz) = *((*v)[j]->Lorentz);
	  
	  return 1;
	}
      }
    }
  }
  else return 1;
  return 0;
}

void Amplitude_Generator::Set_Props(Point* pl,int dep,Single_Amplitude* &first,int* perm)
{
  int ap = 0;
  int lanz = 1;

  Point* preah;
  preah = new Point[dep];

  Point* p;

  short int i,j,k;
  int help = 0;
  top->Copy(pl,prea[0].p,help);

  prea[0].on = 1;

  int sw1;
  Flavour flav[3];
  Complex cpl[4];

  Vertex* v = mo->Get_Vertex();
    
  int first_try = 1;
  
  for (;;) {    
    sw1 = 1;
    if (prea[ap].on) {
      p = 0;
      Next_P(&prea[ap].p[0],p);
      if (p==0) {
	if (first_try==1) p = &prea[ap].p[0];
	else sw1 = 0;
      }
      first_try = 0;
    }
    else sw1 = 0;
    if (sw1) {
      for (i=0;i<v->Max_Number();i++) {
	if ((*v)[i]->on) {
	  flav[0] = p->fl;
	  flav[1] = p->left->fl;
	  flav[2] = p->right->fl;

	  if (p->left->fl  == Flavour(kf::none)) p->left->b  = 0;
	  if (p->right->fl == Flavour(kf::none)) p->right->b = 0;

	  if (flav[0].Majorana()) {
	    if (p->left->fl != Flavour(kf::none)) {
	      if (p->left->fl.IsFermion()) {
		if (p->b*p->left->b == 1)  flav[1] = flav[1].Bar();
	      }
	    }
	    if (p->right->fl != Flavour(kf::none)) {
	      if (p->right->fl.IsFermion()) {
		if (p->b*p->right->b == 1)  flav[2] = flav[2].Bar();
	      }
	    }
	  }
	  else {
	    if (flav[0].IsBoson()) {
	      if (p->left->b   == -1) flav[1] = flav[1].Bar();
	      if (p->right->b  == -1) flav[2] = flav[2].Bar();
	      if (p->left->fl  == Flavour(kf::none)) p->left->b  = -1;
	      if (p->right->fl == Flavour(kf::none)) p->right->b = -1;
	    }
	    else {
	      if (flav[0].IsAnti()) {
		if (p->b*p->left->b == 1)  flav[1] = flav[1].Bar();
		if (p->right->b     ==-1)  flav[2] = flav[2].Bar();
		if (p->left->fl     == Flavour(kf::none)) p->left->b = p->b;
	      }
	      else {
		if (p->b*p->right->b == 1) flav[2] = flav[2].Bar();
		if (p->left->b       ==-1) flav[1] = flav[1].Bar();
		if (p->right->fl     == Flavour(kf::none)) p->right->b = p->b;
	      }
	    }
	  }
	  sw1 = 0;
	  if (Match_Vertex((*v)[i],flav,cpl)) {
	    if (flav[0].Majorana()) {    
	      if (flav[1].IsFermion() && p->left->b==0)  p->left->b  = p->b;
	      if (flav[2].IsFermion() && p->right->b==0) p->right->b = p->b;
	    }
	    sw1 = Check_End(p->left,flav[1],v);	    
	    if (sw1) {
	      sw1 = Check_End(p->right,flav[2],v);
	    }
	  }
	  if (sw1) {
	    int ll = 0;
	    top->Copy(prea[ap].p,preah,ll);
	    if (p->left->fl==Flavour(kf::none))  p->left->fl  = flav[1];
	    if (p->right->fl==Flavour(kf::none)) p->right->fl = flav[2];
	    p->v          = (*v)[i];
	    *(p->Color)   = *((*v)[i]->Color);
	    *(p->Lorentz) = *((*v)[i]->Lorentz);
	  
	    for (k=0;k<4;k++) p->cpl[k] = cpl[k];
	    ll = 0;top->Copy(prea[ap].p,prea[lanz].p,ll);
	    ll = 0;top->Copy(preah,prea[ap].p,ll);
	    prea[lanz].on = 1;
	    lanz++;
	  }
	}
      }
      prea[ap].on = 0;
    }
    if (ap==lanz-1) break;
    ap++;
  }

  int count=0;
  Single_Amplitude* n;
  n = first;
  if (n) while (n->Next) n = n->Next;
  Single_Amplitude* gra;

  
  for (i=0;i<lanz;i++) {
    if (prea[i].on) {
      int sw1 = 1;
      if (AORGTOOLS::rpa.me.Model()==AORGTOOLS::Model_Type::QCD ||
	  AORGTOOLS::rpa.me.Model()==AORGTOOLS::Model_Type::pure_QCD) {
	sw1 = 0;
	for (j=0;j<dep;j++) {
	  if (((prea[i].p[j].fl).IsBoson()) && 
	      (prea[i].p[j].fl!=Flavour(kf::gluon))) sw1++;
	}
	if (AORGTOOLS::rpa.me.Model()==AORGTOOLS::Model_Type::pure_QCD) {
	  if (sw1>0) sw1 = 0;
	  else  sw1 = 1;
	}
	else {
	  if (sw1>1) sw1 = 0;
	  else sw1 = 1;

	}   
      }
      // test if 3-Vertex
      for (j=0;j<dep;j++) {
	if (prea[i].p[j].left!=0) {
	  if ((prea[i].p[j].v)->in[1]==(prea[i].p[j].v)->in[2]) {
	    if (prea[i].p[j].left->number>prea[i].p[j].right->number) {
	      sw1 = 0;
	      break;
	    }
	  }
	}
      }
      //test if 4-vertex
      if (sw1) {
	for (j=0;j<dep;j++) {
	  if (prea[i].p[j].left!=0) {
	    if ((prea[i].p[j].v)->in[1]==(prea[i].p[j].v)->in[2]) {
	      if (prea[i].p[j].left->number>prea[i].p[j].right->number) {
		sw1 = 0;
		break;
	      }
	    }
	  }
	}
      }

      // test if 5-vertex - Starform
      // smallest number must be on the upper side
      if (sw1) {
	int numbers[4];
	short int k;
	for (j=0;j<dep;j++) {
	  if (prea[i].p[j].left!=0) {
	    if ((prea[i].p[j].v)->in[1]==(prea[i].p[j].v)->in[2]) {
	      for (k=0;k<4;k++) numbers[k] = -1;
	      //first gluon
	      if (prea[i].p[j].left->left!=0) {
		if ((prea[i].p[j].left->v)->in[1]==(prea[i].p[j].left->v)->in[2]) {
		  numbers[0] = prea[i].p[j].left->left->number;
		  numbers[1] = prea[i].p[j].left->right->number;
		}
	      }
	      //second
	      if (prea[i].p[j].right->left!=0) {
		if ((prea[i].p[j].right->v)->in[1]==(prea[i].p[j].right->v)->in[2]) {
		  numbers[2] = prea[i].p[j].right->left->number;
		  numbers[3] = prea[i].p[j].right->right->number;
		}
	      }
	      // check if match
	      int sw3 = 1;
	      for (k=0;k<4;k++) {
		if (numbers[k]== -1) {
		  sw3 = 0;
		  break;
		}
	      }
	      if (sw3) {
		if (numbers[0]>numbers[2] || numbers[0]>numbers[3]) sw1 = 0;
	      }
	    }
	  }
	}
      }
      if (sw1) {
	short int k,l;
	for (k=0;k<dep;k++) {
	  if ((prea[i].p[k].number>99) && ((prea[i].p[k].fl).IsBoson())) 
	    prea[i].p[k].number += 100;
	}
	gra = new Single_Amplitude(prea[i].p,b,dep,N,top,BS,fl,shand);
	AORGTOOLS::msg.Tracking()<<"*";AORGTOOLS::msg.Tracking().flush();
	if (first) n->Next = gra;
	else first   = gra; 
	n = gra;
	count++;
      }
    }
  }
  delete[] preah;
}

void Amplitude_Generator::Unite(Point* p,Point* pdel)
{
  int depth = single_top->depth;
  Point psave;
  for (short int i=0;i<depth;i++) {
    if (p[i].number==pdel[i].number) {
      if (p[i].fl!=pdel[i].fl) {
	// Double Counting !!!!!!!!!!!!!!
	int sw1 = 1;
	for (short int j=0;j<p[i].nextra;j++) {
	  if (p[i].extrafl[j]==pdel[i].fl) {
	    sw1 = 0;
	    break;
	  }
	}
	if (sw1==1) {
	  psave = p[i];
	
	  if (p[i].nextra>0) delete[] p[i].extrafl;
	  delete[] p[i].cpl;
	  
	  int nfl  = 1+pdel[i].nextra+p[i].nextra;
	  int ncpl = pdel[i].ncpl+p[i].ncpl;
	  
	  p[i].extrafl = new Flavour[nfl];
	  p[i].cpl     = new Complex[ncpl];
	  
	  //Flavour
	  int count = 0;
	  for (short int j=0;j<psave.nextra;j++)
	    p[i].extrafl[j] = psave.extrafl[j];   
	  count += psave.nextra;
	  
	  p[i].extrafl[count] = pdel[i].fl;
	  count++;
	  for (short int j=0;j<pdel[i].nextra;j++)
	    p[i].extrafl[count+j] = pdel[i].extrafl[j];			
	  p[i].nextra = nfl;
	  
	  //Couplings
	  count = 0;
	  for (short int j=0;j<psave.ncpl;j++)
	    p[i].cpl[j] = psave.cpl[j];   
	  count += psave.ncpl;
	  for (short int j=0;j<pdel[i].ncpl;j++)
	    p[i].cpl[count+j] = pdel[i].cpl[j];			
	  p[i].ncpl = ncpl;
	  
	  //previous couplings too
	  int hit = -1;
	  for (short int j=0;j<depth;j++) {
	    if (p[j].left==&p[i] || p[j].right==&p[i]) {
	      hit = j;
	      break;
	    }
	  }
	  if (hit!=-1) {
	    psave = p[hit];
	    int ncpl = pdel[hit].ncpl+p[hit].ncpl;	  
	    delete[] p[hit].cpl;	  
	    p[hit].cpl = new Complex[ncpl];
	    
	    //Couplings
	    count = 0;
	    for (short int j=0;j<psave.ncpl;j++)
	      p[hit].cpl[j] = psave.cpl[j];   
	    count += psave.ncpl;
	    for (short int j=0;j<pdel[hit].ncpl;j++)
	      p[hit].cpl[count+j] = pdel[hit].cpl[j];			
	    p[hit].ncpl = ncpl;	  
	  }
	  else AORGTOOLS::msg.Error()<<"Error in Amplitude_Generator"<<endl;
	}
      }
    }
  } 
}

int Amplitude_Generator::Single_Compare(Point* p1, Point* p2)
{
  //zero check
  if (p1==0) {
    if (p2==0) return 1;
    else return 0;
  }
  else {
    if (p2==0) return 0;
  }
  
  //Flavour equal....
  if (p1->fl!=p2->fl) return 0;

  //outgoing number equal
  if ((p1->left==0) && (p2->left==0)) {
    if (p1->number!=p2->number) return 0;
                           else return 1;
  }

  //Check extended Color_Functions
  if (p1->Color->String()==p2->Color->String()) {
    if (p1->Color->Next && p2->Color->Next) {
      Color_Function* ctmp1 = p1->Color;
      Color_Function* ctmp2 = p2->Color;
      while(ctmp1->Next) {
	if((ctmp1->Next)->String()!=(ctmp2->Next)->String()) return 0;
	ctmp1 = ctmp1->Next;
	ctmp2 = ctmp2->Next;
      }
    }
    else {
      if (p1->Color->Next!=0 || p2->Color->Next!=0) return 0;
    }
  }
  else return 0;
  
  if (Single_Compare(p1->middle,p2->middle)) {
    int sw1 = Single_Compare(p1->left,p2->left);
    if (sw1) sw1 = Single_Compare(p1->right,p2->right);
    if (sw1==0) {
      sw1 = Single_Compare(p1->left,p2->right);
      if (sw1) sw1 = Single_Compare(p1->right,p2->left);
    }
    return sw1;
  }

  if (Single_Compare(p1->middle,p2->left)) {
    int sw1 = Single_Compare(p1->left,p2->middle);
    if (sw1) sw1 = Single_Compare(p1->right,p2->right);
    if (sw1==0) {
      sw1 = Single_Compare(p1->left,p2->right);
      if (sw1) sw1 = Single_Compare(p1->right,p2->middle);
    }
    return sw1;
  }

  if (Single_Compare(p1->middle,p2->right)) {
    int sw1 = Single_Compare(p1->right,p2->middle);
    if (sw1) sw1 = Single_Compare(p1->left,p2->left);
    if (sw1==0) {
      sw1 = Single_Compare(p1->right,p2->left);
      if (sw1) sw1 = Single_Compare(p1->left,p2->middle);
    }
    return sw1;
  }

  return 0;

  /*
  if (sw1==0) {
    int gluoncheck = 1;
    if ((p1->left!=0) && (p1->right!=0) && (p1->middle==0)) {
      if (p1->left->fl!=p1->right->fl) gluoncheck = 0;
    }
    //new
    if ((p1->left!=0) && (p1->right!=0) && (p1->middle!=0)) {
      //this means no two identical flavours available
      if ((p1->left->fl!=p1->right->fl) && (p1->left->fl!=p1->middle->fl)) gluoncheck = 0;
    }
    if ((p2->left!=0) && (p2->right!=0) && (p2->middle==0)) {
      if (p2->left->fl!=p2->right->fl) gluoncheck = 0;
    }
    //new
    if ((p2->left!=0) && (p2->right!=0) && (p2->middle!=0)) {
      //this means no two identical flavours available
      if ((p2->left->fl!=p2->right->fl) && (p2->left->fl!=p2->middle->fl)) gluoncheck = 0;
    }
    if (gluoncheck==1) { 
      sw1 = Single_Compare(p1->left,p2->right);
      if (sw1) Single_Compare(p1->right,p2->left);
      if ((p1->middle!=0 || p2->middle!=0) && sw1) 
	sw1 = Single_Compare(p1->middle,p2->middle);
    }
    else {
      //two outgoing majorana's
      int majoranacheck = 1;
      if ((p1->left!=0) && (p1->right!=0) && (p1->middle==0)) {
	if (!p1->left->fl.Majorana()) majoranacheck = 0;
	if (!p1->right->fl.Majorana()) majoranacheck = 0;
      }
      //new
      if ((p1->left!=0) && (p1->right!=0) && (p1->middle!=0)) {
	//this means no two Majoranas available
	if ((!p1->left->fl.Majorana()) || (!p1->right->fl.Majorana())) {
	  if (!p1->middle->fl.Majorana()) majoranacheck = 0;
	}
      }
      if ((p2->left!=0) && (p2->right!=0) && (p2->middle==0)) {
	if (!p2->left->fl.Majorana()) majoranacheck = 0;
	if (!p2->right->fl.Majorana()) majoranacheck = 0;
      }
      //new
      if ((p2->left!=0) && (p2->right!=0) && (p2->middle!=0)) {
	//this means no two Majoranas available
	if ((!p2->left->fl.Majorana()) || (!p2->right->fl.Majorana())) {
	  if (!p2->middle->fl.Majorana()) majoranacheck = 0;
	}
      }
      if (majoranacheck == 0 && p1->fl.Majorana()) {
	majoranacheck = 1;
	if ((p1->left!=0) && (p1->right!=0) && (p1->middle==0)) {
	  if (!p1->left->fl.Majorana() && !p1->right->fl.Majorana()) majoranacheck = 0;
	}
	//new
	if ((p1->left!=0) && (p1->right!=0) && (p1->middle!=0)) {
	  if ((!p1->left->fl.Majorana())  && 
	      (!p1->right->fl.Majorana()) && 
	      (!p1->middle->fl.Majorana())  ) majoranacheck = 0;
	}
	if ((p2->left!=0) && (p2->right!=0) && (p2->middle==0)) {
	  if (!p2->left->fl.Majorana() && !p2->right->fl.Majorana()) majoranacheck = 0;
	}
	//new
	if ((p2->left!=0) && (p2->right!=0) && (p2->middle!=0)) {
	  if ((!p2->left->fl.Majorana())  && 
	      (!p2->right->fl.Majorana()) &&
	      (!p2->middle->fl.Majorana())   ) majoranacheck = 0;
	}
      }
      if (majoranacheck==1) { 
	sw1 = Single_Compare(p1->left,p2->right);
	if (sw1) Single_Compare(p1->right,p2->left);
	if ((p1->middle!=0 || p2->middle!=0) && sw1) 
	  sw1 = Single_Compare(p1->middle,p2->middle);
      }
    }
  }
  return sw1;
  */

}

void Amplitude_Generator::Kill_Off(Single_Amplitude* &first)
{
  Single_Amplitude* last;
  last = first;
  Single_Amplitude* f1 = first;
  Single_Amplitude* f2;
  while (f1) {
    if (f1->on==0) {
      if (f1==first) {
	first = f1->Next;
	f1 = f1->Next;
	delete last;
	last = first;
      }
      else {
	last->Next = f1->Next;
	f2 = f1;
	f1 = f1->Next;
	delete f2;
      }
    }
    else {
      last = f1;
      f1 = f1->Next;
    }
  }
}

void Amplitude_Generator::Compare(Single_Amplitude* &first)
{
  Single_Amplitude* f1;
  Single_Amplitude* f2;
  Point* p1;
  Point* p2;  
  int count1 = 1;
  int count2 = 1;
  
  f1 = first;
  while (f1) { 
    p1 = f1->GetPointlist();
    f2 = f1->Next;
    count2 = count1+1;
    while (f2) {
      p2 = f2->GetPointlist();
      int sw1 = Single_Compare(p1,p2);
      if (sw1==1) {
	f2->on = 0;
	AORGTOOLS::msg.Debugging()<<count1<<" = "<<count2<<endl;
      }
      count2++;
      f2 = f2->Next;
    }
    count1++;
    f1 = f1->Next;
  }
  Kill_Off(first);
}

//==================================================================

Point* Amplitude_Generator::FindNext(Point* p)
{
  if (p==0) return 0;
  if (p->left->m==0) return p;
  if (p->right->m==0) return p;
  if ((p->middle!=0) && (p->middle->m==0)) return p;
  
  FindNext(p->left);
  FindNext(p->right);
  FindNext(p->middle);
}

int Amplitude_Generator::ShrinkProps(Point*& p,Point*& pnext, Point*& pcopy, Point*& beg_pcopy,
				     vector<Point*>& pcollist)
{
  if (p->left==0 || pnext->left==0) return 0;
  
  if (p->v->nleg==4 || pnext->v->nleg==4) return 0;
  
  if (pnext->m==1) return 0;
    
  int hit = 0;
  int in[4] = {0};
  Point* ptmp;
  
  Flavour flav[4];
  
  flav[0] = p->fl;
  if (p->number<99 && p->b==1) in[0] = 1;

  if (p->left->number==pnext->number) {
    flav[1] = pnext->left->fl;
    if (pnext->left->number<99 && pnext->left->b==-1)   in[1] = 1;
    flav[2] = pnext->right->fl;
    if (pnext->right->number<99 && pnext->right->b==-1) in[2] = 1;
    flav[3] = p->right->fl;
    if (p->right->number<99 && p->right->b==-1)         in[3] = 1;
  }
  if (p->right->number==pnext->number) {
    flav[1] = p->left->fl;
    if (p->left->number<99 && p->left->b==-1)           in[1] = 1;
    flav[2] = pnext->left->fl;
    if (pnext->left->number<99 && pnext->left->b==-1)   in[2] = 1;
    flav[3] = pnext->right->fl;
    if (pnext->right->number<99 && pnext->right->b==-1) in[3] = 1;
  }  
 
  //barflags 
  
  Vertex* v = mo->Get_Vertex();
  
  for (short int i=0;i<v->Max_Number4();i++) {
    if ((*v)(i)->on) {
      
      Single_Vertex test;
      
      test = *((*v))(i);
      
      if (in[0]) test.in[0] = test.in[0].Bar();
      if (in[1]) test.in[1] = test.in[1].Bar();
      if (in[2]) test.in[3] = test.in[3].Bar();
      if (in[3]) test.in[2] = test.in[2].Bar();

      if (flav[0]==test.in[0] && 
	  flav[1]==test.in[1] &&
	  flav[2]==test.in[3] &&
	  flav[3]==test.in[2]) 
	{
	  AORGTOOLS::msg.Debugging()<<"4 leg Vertex found !!!"<<endl;
	  
	  hit = 1;
	  
	  pcopy->v      = (*v)(i);
	  for (short int k=0;k<4;k++) pcopy->cpl[k] = (*v)(i)->cpl[k];
	  //set pcopy legs
	  
	  if (p->left->number==pnext->number) {	    
	    pcopy->middle = pcopy->left->right;
	    pcopy->left   = pcopy->left->left;
	    pcopy->right  = pcopy->right;
	  }
	  if (p->right->number==pnext->number) {	    
	    pcopy->left   = pcopy->left;
	    pcopy->middle = pcopy->right->left;
	    pcopy->right  = pcopy->right->right;
	  }
	  
	  //setting the contraction flags
	  
	  //if (p->right->number == pnext->number) p->right->m = 1;
	  //if (p->left->number == pnext->number) p->left->m = 1;
	  pnext->m  = 1;
	  
	  if ((*v)(i)->ncf==1) {
	    
	    *(pcopy->Color)   = *((*v)(i)->Color);
	    *(pcopy->Lorentz) = *((*v)(i)->Lorentz);
	    
	    break;
	  }
	  else {
	    
	    for (short int k=0;k<(*v)(i)->ncf;k++) {

	      *(pcopy->Color)   = ((*v)(i)->Color)[k];

	      if (((*v)(i)->Color[k]).Next!=0) {
		Color_Function* cforig = ((*v)(i)->Color[k]).Next;
		Color_Function* cfcopy = pcopy->Color;
		while (cforig) {
		  cfcopy->Next = new Color_Function(*cforig);
		  cfcopy = cfcopy->Next;
		  cforig = cforig->Next;
		}
	      }	      
	      
	      *(pcopy->Lorentz) = ((*v)(i)->Lorentz)[k];
	      
	      Color_Function* cfmemo = pcopy->Color;
	      
	      //set the contraction indices: 4 -> (prop->number) 
	      while (cfmemo) {
		for (short int i=0;i<3;i++) {
		  if (cfmemo->type==cf::D && i==2) break;
		  if (cfmemo->partarg[i]==4) {
		    cfmemo->partarg[i] = pnext->number;
		  }
		  }
		cfmemo = cfmemo->Next;
	      }

	      //fill the vector pcollist
	      int ll = 0;
	      Point* ptmp = new Point[single_top->depth];
	      top->Copy(beg_pcopy,ptmp,ll);
	      pcollist.push_back(ptmp);
	    }
	    break;
	  }
	}
      else hit = 0;
    }
  }
  return hit;
}

int Amplitude_Generator::EvalPointlist(Point*& porig, Point*& pcopy,Point*& beg_pcopy,
				       vector<Point*>& pcollist)
{
  if(porig==0) return 0;
  
  if (ShrinkProps(porig,porig->right,pcopy,beg_pcopy,pcollist))    return 1;
  
  if (ShrinkProps(porig,porig->left,pcopy,beg_pcopy,pcollist))     return 1;
  
  if (EvalPointlist(porig->left,pcopy->left,beg_pcopy,pcollist))   return 1;
  if (EvalPointlist(porig->right,pcopy->right,beg_pcopy,pcollist)) return 1;
  if (EvalPointlist(porig->middle,pcopy->middle,beg_pcopy,pcollist)) return 1;

  return 0;
}

void Amplitude_Generator::CheckFor4Vertices(Single_Amplitude* &first)
{
  Single_Amplitude* f1 = first;
  Single_Amplitude* extraAmpl;
  Single_Amplitude* f2;
  Point* p;
  int dep       = single_top->depth;
  Point* pcopy  = new Point[dep];
  vector<Point*> pcollist;//vector with different color structures
  int counter   = 0;
  int amplcount = 0;
  
  AORGTOOLS::msg.Debugging()<<"=============== CheckFor4Vertices ======================="<<endl;  
  
  while (f1) { 
    amplcount++;
    
    p = f1->GetPointlist();
    
    for (int i=0; i<dep; i++) p[i].m = 0; 
    
    while (p) {
      //p initial Pointlist with m flags set, pcopy is the contracted Pointlist !!!
      int ll = 0;
      top->Copy(p,pcopy,ll);
      if (EvalPointlist(p,pcopy,pcopy,pcollist)) {
	if (pcollist.size()==0) {
	  counter++;
	  extraAmpl = new Single_Amplitude(pcopy,b,dep,N,top,BS,fl,shand);
	  extraAmpl->Next = 0;
	  
	  f2 = first;
	  
	  while (f2) {
	    if (f2->Next==0) {
	      f2->Next = extraAmpl;
	      break;
	    }
	    f2 = f2->Next;
	  }
	}
	else { 
	  for (short int j=0;j<pcollist.size();j++) {
	    counter++;
	    extraAmpl = new Single_Amplitude(pcollist[j],b,dep,N,top,BS,fl,shand);
	    extraAmpl->Next = 0;
	    
	    f2 = first;
	    //muss ueberdacht werden, wohin zeigt der Zeiger ???
	    while (f2) {
	      if (f2->Next==0) {
		f2->Next = extraAmpl;
		break;
	      }
	      f2 = f2->Next;
	    }
	  }
	  //erase all elements of pcollist
	  pcollist.clear();
	}
	p = FindNext(p);
      }
      else break;
    }
    f1 = f1->Next;
  }
  AORGTOOLS::msg.Debugging()<<"=========================================================="<<endl;
}

//=====================================================================


Single_Amplitude* Amplitude_Generator::Matching()
{
  int nloop = N-1;
  short int i,j;
  int* ii = new int[nloop];
  for (i=0;i<nloop;i++) ii[i] = 0;
  int* perm = new int[N];
  int sw1;
  int qsum,l1sum,l2sum,l3sum; 
  int chsum,neusum;           
  Single_Amplitude* first_amp;
  int over = 0;
  perm[0]  = 0;
  first_amp = 0;

  long int count = 1;

  for(;;) {
    sw1 = 1;
    for (j=0;j<nloop;j++) perm[j+1] = ii[j];
    
    for(i=0;i<N;i++) {
      for (j=i+1;j<N;j++) 
	if (perm[i]==perm[j]) {sw1 = 0;break;}
    }

    if (sw1) count++;

    if (sw1) { 
      if (fl[0].IsQuark() && !fl[0].IsAnti()) {
	if (!(fl[perm[N-1]].IsQuark() && 
	      (b[perm[N-1]]==-1 && (fl[perm[N-1]].IsAnti() || fl[perm[N-1]].Majorana())) ||
	      (b[perm[N-1]]==1  && (!(fl[perm[N-1]].IsAnti()) || fl[perm[N-1]].Majorana())) ))  sw1 = 0;
      }
      if (fl[0].IsQuark() && fl[0].IsAnti()) {
	if (!(fl[perm[1]].IsQuark() && 
	      (b[perm[1]]==-1 && 
	       (!fl[perm[1]].IsAnti() || fl[perm[1]].Majorana())) ||
	      (b[perm[1]]==1 &&  
	       ((fl[perm[1]].IsAnti()) || fl[perm[1]].Majorana())) ))  sw1 = 0;
      }
      //Anti-Leptons too!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (fl[0].IsLepton()) {
	if (!(fl[perm[N-1]].IsLepton() //|| fl[perm[N-1]].IsSlepton() 
	      && 
	      (b[perm[N-1]]==-1 && (fl[perm[N-1]].IsAnti() || fl[perm[N-1]].Majorana())) ||
	      (b[perm[N-1]]==1  && (!(fl[perm[N-1]].IsAnti()) || fl[perm[N-1]].Majorana())) )) 
	  sw1 = 0;
      }
      if (fl[0]==Flavour(kf::photon)) {
	//changed to any place 
	//if (fl[perm[1]]==Flavour(kf::photon)) sw1 = 0;
	//if (fl[perm[N-1]]==Flavour(kf::photon)) sw1 = 0;
      }
    }
    if (sw1) {
      chsum=neusum=qsum=l1sum=l2sum=l3sum=0;
      for (j=0;j<N-1;j++) {
	/*
	if (fl[perm[j]].IsChargino()) {
	  if (fl[perm[j]].IsAnti()) {
	    if (b[perm[j]]==1) chsum--;
	                  else chsum++;
	  }
	  else {
	    if (b[perm[j]]==1) chsum++;
	                  else chsum--;
	  }
	  if (!fl[0].IsAnti()) {
	    if ( b[perm[N-1]]==-1 && fl[perm[N-1]].IsAnti() && chsum>0 ) 
	      {sw1=0;break;} 
	    if ( b[perm[N-1]]==1 && !(fl[perm[N-1]].IsAnti()) && chsum>0 )
	      {sw1=0;break;} 
	  }
	  else {
	    if ( b[perm[1]]==-1 && fl[perm[1]].IsAnti() && chsum<0 ) 
	      {sw1=0;break;}
	    if ( b[perm[1]]==1 && !(fl[perm[1]].IsAnti()) && chsum<0 )
	      {sw1=0;break;}
	  }
	}
	*/
	if (fl[perm[j]].IsNeutralino()) {
	  // no rules!!!!
	}
	if (fl[perm[j]].IsQuark()) {
	  if (fl[perm[j]].IsAnti()) {
	    if (b[perm[j]]==1) qsum--;
	    else qsum++;
	  }
	  else {
	    if (b[perm[j]]==1) qsum++;
	    else qsum--;
	  }
	  if (!fl[0].IsAnti()) {
	    if ( b[perm[N-1]]==-1 && fl[perm[N-1]].IsAnti() && qsum>0 ) 
	      {sw1=0;break;} // s-channel
	    if ( b[perm[N-1]]==1 && !(fl[perm[N-1]].IsAnti()) && qsum>0 )
	      {sw1=0;break;} // t-channel,
	  }
	  else {
	    if ( b[perm[1]]==-1 && fl[perm[1]].IsAnti() && qsum<0 ) 
	      {sw1=0;break;} // s-channel
	    if ( b[perm[1]]==1 && !(fl[perm[1]].IsAnti()) && qsum<0 )
	      {sw1=0;break;} // t-channel,	    
	  }
	}
	if (fl[perm[j]].IsLepton()) {
	  if ( (fl[perm[j]]==Flavour(kf::e)) || 
	       (fl[perm[j]]==Flavour(kf::nue)) ) {
	    if (b[perm[j]]==1) l1sum++;
	    else l1sum--;
	  }
	  if ( (fl[perm[j]]==Flavour(kf::e).Bar()) || 
	       (fl[perm[j]]==Flavour(kf::nue).Bar()) ) {
	    if (b[perm[j]]==1) l1sum--;
	    else l1sum++;
	  }
	  if ( (fl[perm[j]]==Flavour(kf::mu)) || 
	       (fl[perm[j]]==Flavour(kf::numu)) ) {
	    if (b[perm[j]]==1) l2sum++;
	    else l2sum--;
	  }
	  if ( (fl[perm[j]]==Flavour(kf::mu).Bar()) || 
	       (fl[perm[j]]==Flavour(kf::numu).Bar()) )  {
	    if (b[perm[j]]==1) l2sum--;
	    else l2sum++;
	  }
	  if ( (fl[perm[j]]==Flavour(kf::tau)) || 
	       (fl[perm[j]]==Flavour(kf::nutau)) ) {
	    if (b[perm[j]]==1) l3sum++;
	    else l3sum--;
	  }
	  if ( (fl[perm[j]]==Flavour(kf::tau).Bar()) || 
	       (fl[perm[j]]==Flavour(kf::nutau).Bar()) )  {
	    if (b[perm[j]]==1) l3sum--;
	    else l3sum++;
	  }
	}
	if ( (b[perm[N-1]]==-1) && (fl[perm[N-1]].IsAnti()) && 
	     ( (l1sum>0) || (l2sum>0) || (l3sum>0))  ) {sw1=0;break;} // s-channel, e+e-.
	if ( (b[perm[N-1]]==1) && !(fl[perm[N-1]].IsAnti()) &&
	     ( (l1sum>0) || (l2sum>0) || (l3sum>0)) )  {sw1=0;break;} // t-channel, e+e-.
      }
    }
    if (AORGTOOLS::rpa.gen.Tracking()) {
      if (!(count%1000000)) {long int fak =1;
      for (long int k=1;k<=N;k++) fak *= k;
      //AORGTOOLS::msg.Out()<<"Count: "<<count<<"/"<<fak<<endl;//AORGTOOLS::msg.Out()<<"|";AORGTOOLS::msg.Out().flush();}
      }
    } 
    if (sw1) {
      long int fak =1;
      for (long int k=1;k<=N;k++) fak *= k;
      for (j=0;j<single_top->number;j++) {
	perm++;
	int pnum = 100;
	Set_End(&single_top->p[j][0],perm,pnum);
	perm -= N;
	single_top->p[j][0].number = *perm;
	single_top->p[j][0].fl     = fl[*perm];
	single_top->p[j][0].b      = b[*perm];
	Set_Props(single_top->p[j],2*N-3,first_amp,perm);	
	//Print_P(&single_top->p[j][0]);
      }
    }     
    for (j=nloop-1;j>=0;j--) {
      if ((ii[j]+1)<N) {
	ii[j]++;            
	break;
      }
      else {
	ii[j] = 0;
	if (j==0) over = 1;
      }
    }
    if (over) break;
  }
  delete[] ii;
  delete[] perm;
  
  CheckFor4Vertices(first_amp);
  
  Compare(first_amp);
  
  Amplitude_Manipulator(N,fl,b).FixSign(first_amp);
  
  return first_amp;
}





















