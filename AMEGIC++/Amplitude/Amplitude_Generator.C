//#include <iomanip>
#include "AMEGIC++/Amplitude/Amplitude_Generator.H"
#include "AMEGIC++/Amplitude/Amplitude_Manipulator.H"
#include "AMEGIC++/Main/Tools.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Permutation.H"
#include <cassert>

namespace AMEGIC {
  class Compare_Pre_Amplitudes {
  public:
    int operator()(const AMEGIC::Pre_Amplitude & a, const AMEGIC::Pre_Amplitude & b) {
      if (a.perm<b.perm) return 1;
      return 0;
    }
  };
}

using namespace AMEGIC;
using namespace ATOOLS;
using namespace MODEL;
using namespace std;

Amplitude_Generator::Amplitude_Generator(int _no,Flavour* _fl,int* _b,
					 Amegic_Model * _model,Topology * _top,
					 std::vector<double> _order,int _ntchan_min,
					 Basic_Sfuncs* _BS,String_Handler* _shand, bool create_4V) 
  : fl(_fl), b(_b), p_model(_model), top(_top), 
    N(_no), order(_order), ntchan_min(_ntchan_min),
    BS(_BS), shand(_shand), m_create_4V(create_4V)
{
  single_top = top->Get(N-2);
  
  // 2 incoming
  prea.push_back(Pre_Amplitude());
  prea.back().p = new Point[single_top->depth];

  // fill hash table
  
  Vertex * v = p_model->p_vertex;
 
  for (int i=0;i<v->MaxNumber();++i) {
    if ((*v)[i]->on) {
      v_table[(*v)[i]->in[0]].push_back((*v)[i]);
    }
  }
}

Amplitude_Generator::~Amplitude_Generator() 
{
  for(size_t i=0;i<prea.size();++i) delete [] prea[i].p;

  for(unsigned int i=0;i<prea_table.size();i++) delete[] prea_table[i].p;
  prea_table.clear();
}

void Amplitude_Generator::Set_End(Point* p,int* &perm,int& pnum)
{
  assert(p != NULL);
  p->b     = 0;
  p->fl    = Flavour(kf_none);
  if ((p->left==0) && (p->right==0)) {
    p->number = *perm;
    p->fl = fl[*perm];
    p->b  = b[*perm];
    if (p->Lorentz) delete p->Lorentz;
    if (p->fl.IsBoson()) {
      p->Lorentz = LF_Getter::GetObject("Pol",LF_Key());
      p->Lorentz->SetParticleArg(0);
    }
    else {
      p->Lorentz = LF_Getter::GetObject("None",LF_Key());
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
    if ((p->left->fl==Flavour(kf_none)) || (p->right->fl==Flavour(kf_none))) {
      hit = p;
      return;
    }
  }
  Next_P(p->left,hit); 
  Next_P(p->right,hit);
}


void Amplitude_Generator::Print_P(Point* p)
{
  if (!(msg_LevelIsDebugging())) return;
  if ((p->left==0) && (p->right==0)) {
    msg_Out()<<"  "<<p->fl<<"("<<p->b<<")"<<endl;
    return;
  }
  if (p->cpl.size()>1) msg_Out()<<"cpl: "<<p->cpl[0]<<" "<<p->cpl[1]<<"\n";
  msg_Indent();
  msg_Out()<<"left : \n";
  Print_P(p->left);
  msg_Out()<<"right : \n";
  Print_P(p->right);
  if(p->middle){
    msg_Out()<<" middle : \n";
    Print_P(p->middle);
  }
}

int Amplitude_Generator::MatchVertex(AMEGIC::Single_Vertex* v,Flavour* flav,vector<Complex>& cpl)
{
  if (v->dec>0) return false;
  if (flav[0] == v->in[0]) {
    int hit = 1;
    if (flav[1] != Flavour(kf_none)) {if (flav[1] != v->in[1]) hit = 0;}
    else { flav[1] = v->in[1];}
    if (flav[2] != Flavour(kf_none)) {if (flav[2] != v->in[2]) hit = 0;}
    else { flav[2] = v->in[2];}
    if (hit==1) {
      cpl.clear();
      for (size_t j=0;j<v->cpl.size();j++) cpl.push_back(v->Coupling(j));
      return 1;
    }
  }
  return 0;
}

int Amplitude_Generator::CheckEnd(Point* p,Flavour infl) 
{
  if (p==0) return 1;
  if (p->left==0) return 1;
  if (((p->left->fl)!=Flavour(kf_none)) && ((p->right->fl)!=Flavour(kf_none))) { 
    Flavour flav[3];
    Flavour s_flav[3];
    vector <Complex> cpl;
    cpl.clear();

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

    s_flav[0]=flav[0];
    s_flav[1]=flav[1];
    s_flav[2]=flav[2];

    // look for flav[0]
    Vertex_List & vl= v_table[flav[0]];

    for (size_t j=0;j<vl.size();j++) {
      if (MatchVertex(vl[j],flav,cpl)) {
	p->cpl.clear();
	for (size_t k=0;k<cpl.size();k++) p->cpl.push_back(cpl[k]);
	p->v = vl[j];
	*p->Color = vl[j]->Color.back();
	if (p->Lorentz) delete p->Lorentz;
	p->Lorentz = vl[j]->Lorentz.front()->GetCopy();
	p->t = vl[j]->t;
	return 1;
      }
    }
  }
  else return 1;
  return 0;
}

void Amplitude_Generator::SetProps(Point* pl,int dep,Single_Amplitude* &first,int* perm,int topcount, int permcount)
{
  int ap = 0;
  int lanz = 1;

  Point* preah;
  preah = new Point[dep];

  Point* p;

  int help = 0;

  top->Copy(pl,prea[0].p,help);

  prea[0].on = 1;

  int sw1;
  Flavour flav[3];
  Flavour s_flav[3];
  vector<Complex> cpl;

  int first_try = 1;
  
  for (;;) {    
    for (int k(prea.size());k<=ap;++k) {
      prea.push_back(Pre_Amplitude());
      prea.back().p = new Point[single_top->depth];
    }
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
      flav[0] = p->fl;
      flav[1] = p->left->fl;
      flav[2] = p->right->fl;

      if (p->left->fl  == Flavour(kf_none)) p->left->b  = 0;
      if (p->right->fl == Flavour(kf_none)) p->right->b = 0;

      if (flav[0].Majorana()) {
	if (p->left->fl != Flavour(kf_none)) {
	  if (p->left->fl.IsFermion()) {
	    if (p->b*p->left->b == 1)  flav[1] = flav[1].Bar();
	  }
	  else if (p->left->b   == -1) flav[1] = flav[1].Bar();
	}
	if (p->left->fl==Flavour(kf_none)) p->left->b = p->b;
	
	if (p->right->fl != Flavour(kf_none)) {
	  if (p->right->fl.IsFermion()) {
	    if (p->b*p->right->b == 1)  flav[2] = flav[2].Bar();
	}
	else if (p->right->b  == -1) flav[2] = flav[2].Bar();  
	}
	if (p->right->fl==Flavour(kf_none)) p->right->b = p->b;
      }
      else {
	if (flav[0].IsBoson()) {
	  if (p->left->b   == -1) flav[1] = flav[1].Bar();
	  if (p->right->b  == -1) flav[2] = flav[2].Bar();
	  if (p->left->fl  == Flavour(kf_none)) p->left->b  = -1;
	  if (p->right->fl == Flavour(kf_none)) p->right->b = -1;
	}
	else {
	  if (flav[0].IsAnti()) {
	    if (p->b*p->left->b == 1)  flav[1] = flav[1].Bar();
	    if (p->right->b     ==-1)  flav[2] = flav[2].Bar();
	    if (p->left->fl     == Flavour(kf_none)) p->left->b = p->b;
	  }
	  else {
	    if (p->b*p->right->b == 1) flav[2] = flav[2].Bar();
	    if (p->left->b       ==-1) flav[1] = flav[1].Bar();
	    if (p->right->fl     == Flavour(kf_none)) p->right->b = p->b;
	  }
	}
      }
      s_flav[0]=flav[0];
      s_flav[1]=flav[1];
      s_flav[2]=flav[2];

      // look for flav[0]
      Vertex_List & vl= v_table[flav[0]];

      for (size_t i=0;i<vl.size();i++) {
	sw1 = 0;
	flav[0]=s_flav[0];
	flav[1]=s_flav[1];
	flav[2]=s_flav[2];

	if (MatchVertex(vl[i],flav,cpl)) {
	  if (flav[0].Majorana()) {    
	    if (flav[1].IsFermion() && p->left->b==0)  p->left->b  = p->b;
	    if (flav[2].IsFermion() && p->right->b==0) p->right->b = p->b;
	  }
	  sw1 = CheckEnd(p->left,flav[1]);
	  if (sw1) {
	    sw1 = CheckEnd(p->right,flav[2]);
	  }
	}
	if (sw1) {
	  for (int k(prea.size());k<=lanz;++k) {
	    prea.push_back(Pre_Amplitude());
	    prea.back().p = new Point[single_top->depth];
	  }
	  //match
	  int ll = 0;
	  top->Copy(prea[ap].p,preah,ll);
	  if (p->left->fl==Flavour(kf_none))  p->left->fl  = flav[1];
	  if (p->right->fl==Flavour(kf_none)) p->right->fl = flav[2];
	  p->v          = vl[i];
	  *p->Color = vl[i]->Color.back();
	  if (p->Lorentz) delete p->Lorentz;
	  p->Lorentz = vl[i]->Lorentz.front()->GetCopy();
	  p->t = vl[i]->t;
	  
	  p->cpl.clear();
	  for (size_t k=0;k<cpl.size();k++) p->cpl.push_back(cpl[k]);
	  ll = 0;top->Copy(prea[ap].p,prea[lanz].p,ll);
	  ll = 0;top->Copy(preah,prea[ap].p,ll);
	  prea[lanz].on = 1;
	  lanz++;
	}
      }
      prea[ap].on = 0;
    }
    if (ap==lanz-1) break;
    ap++;
  }

  for (int i=0;i<lanz;i++) {
    if (prea[i].on) {
      int ll = 0;
      Point* ph;
      ph = new Point[dep];
      top->Copy(prea[i].p,ph,ll);
      prea_table.push_back(Pre_Amplitude(ph,topcount,permcount));
      //      counter_table.push_back(counters);
    }
  }

  delete[] preah;
}


void Amplitude_Generator::CreateSingleAmplitudes(Single_Amplitude * & first) {

  int count=0;
  Single_Amplitude* n;
  int dep=single_top->depth;

  n = first;
  if (n) while (n->Next) n = n->Next;
  Single_Amplitude* gra;
  for (size_t i=0;i<prea_table.size();i++) {
    int sw1 = 1;
    // test if 3-Vertex
    if (sw1) {
      for (int j=0;j<dep;j++) {
	if (prea_table[i].p[j].left!=0) {
	  if ((prea_table[i].p[j].v)->in[1]==(prea_table[i].p[j].v)->in[2]) {
	    if (prea_table[i].p[j].left->number>prea_table[i].p[j].right->number) {
	      sw1 = 0;
	      break;
	    }
	  }
	}
      }
    }
    //test if 4-vertex
    if (sw1) {
      for (int j=0;j<dep;j++) {
	if (prea_table[i].p[j].left!=0) {
	  if ((prea_table[i].p[j].v)->in[1]==(prea_table[i].p[j].v)->in[2]) {
	    if (prea_table[i].p[j].left->number>prea_table[i].p[j].right->number) {
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
      for (int j=0;j<dep;j++) {
	if (prea_table[i].p[j].left!=0) {
	  if ((prea_table[i].p[j].v)->in[1]==(prea_table[i].p[j].v)->in[2]) {
	    for (int k=0;k<4;k++) numbers[k] = -1;
	    //first gluon
	    if (prea_table[i].p[j].left->left!=0) {
	      if ((prea_table[i].p[j].left->v)->in[1]==(prea_table[i].p[j].left->v)->in[2]) {
		numbers[0] = prea_table[i].p[j].left->left->number;
		numbers[1] = prea_table[i].p[j].left->right->number;
	      }
	    }
	    //second
	    if (prea_table[i].p[j].right->left!=0) {
	      if ((prea_table[i].p[j].right->v)->in[1]==(prea_table[i].p[j].right->v)->in[2]) {
		numbers[2] = prea_table[i].p[j].right->left->number;
		numbers[3] = prea_table[i].p[j].right->right->number;
	      }
	    }
	    // check if match
	    int sw3 = 1;
	    for (int k=0;k<4;k++) {
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
      for (int k=0;k<dep;++k) {
	if ((prea_table[i].p[k].number>99) && ((prea_table[i].p[k].fl).IsBoson())) 
	  prea_table[i].p[k].number += 100;
      }
      if (CheckOrders(prea_table[i].p) && 
	  CheckTChannels(prea_table[i].p)) {
	gra = new Single_Amplitude(prea_table[i].p,prea_table[i].top,prea_table[i].perm,b,dep,N,top,BS,fl,shand);
	count++;
	if (first) n->Next = gra;
	else first   = gra; 
	n = gra;
      }
    }
  }
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
	  p[i].cpl.clear();
	  
	  int nfl  = 1+pdel[i].nextra+p[i].nextra;
	  
	  p[i].extrafl = new Flavour[nfl];
	  
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
	  p[i].cpl.clear();
	  for (size_t j=0;j<psave.Ncpl();j++) p[i].cpl.push_back(psave.cpl[j]);   
	  for (size_t j=0;j<pdel[i].Ncpl();j++) p[i].cpl.push_back(pdel[i].cpl[j]);			
	  
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
	    p[hit].cpl.clear();
	    
	    //Couplings
	    for (size_t j=0;j<psave.Ncpl();j++) p[hit].cpl.push_back(psave.cpl[j]);   
	    for (size_t j=0;j<pdel[hit].Ncpl();j++) p[hit].cpl.push_back(pdel[hit].cpl[j]);			
	  }
	  else 
	    msg_Error()<<"ERROR in Amplitude_Generator"<<endl
			       <<"   Continue and hope for the best ..."<<std::endl;
	}
      }
    }
  } 
}

int Amplitude_Generator::CompareColors(Point* p1,Point* p2)
{
  Color_Function * c1 = p1->Color;
  Color_Function * c2 = p2->Color;

  if (c1->Next()==0 && c2->Next()==0) return 1;

  if ((c1->Next()!=0 && c2->Next()==0) || (c1->Next()==0 && c2->Next()!=0)) return 0; 

  int l1[3],l2[3],l1n[3],l2n[3];
  
  if (c1->Type()==cf::F) {
    if (c1->Next()->Type()==cf::F) {
      if (c1->String()==c2->String()) {
	for (int i=0;i<3;i++){
	  l1[i]=c1->ParticleArg(i);
	  l2[i]=c2->ParticleArg(i);
	  l1n[i]=c1->Next()->ParticleArg(i);
	  l2n[i]=c2->Next()->ParticleArg(i);
	}
	if (l1[0]!=l2[0] && l1n[0]!=l2n[0]) {
	  if (l1[1]==l2[1] && l1[2]==l2[2] && l1n[1]==l2n[1] && l1n[2]==l2n[2]) return 0;
	}
	else return 1;
      }
      else if (c1->ParticleArg(1)==c2->Next()->ParticleArg(0) &&
	       c1->Next()->ParticleArg(1)==c2->ParticleArg(1) ) {
	return 1;
      }
      else return 0;
    }
    else 
      msg_Error()<<"ERROR in Amplitude_Generator::CompareColors :"<<std::endl
		 <<"   Color structure not supported. Continue and hope for the best. "<<endl;
  }

  for (int i=0;i<3;i++){
    l1[i]=c1->ParticleArg(i);
    l2[i]=c2->ParticleArg(i);
    l1n[i]=c1->Next()->ParticleArg(i);
    l2n[i]=c2->Next()->ParticleArg(i);
  }
  if (c1->Next()->Type()!=cf::T)
    msg_Error()<<"ERROR in Amplitude_Generator::CompareColors :"<<std::endl
	       <<"   Unexpected sequence in color structure. Continue and hope for the best. "<<endl;
  
  if (l1[0]!=l2[0] && l1n[0]!=l2n[0]) {
    if (l1[1]==l2[1] && l1[2]==l2[2] && l1n[1]==l2n[1] && l1n[2]==l2n[2]) return 0;
  }
  return 1;
}

int Amplitude_Generator::SingleCompare(Point* p1,Point* p2)
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

  if ((p1->left==0) && (p2->left!=0)) return 0;
  if ((p1->left!=0) && (p2->left==0)) return 0;

  //Check extended Color_Functions
  if (p1->Color->Type()!=p2->Color->Type()) return 0;
  
  if (!CompareColors(p1,p2)) return 0;

  
  // return 1 if equal and 0 if different
  
  if (p1->Lorentz->Type()=="C4GS" && p2->Lorentz->Type()=="C4GS") 
    return Compare5Vertex(p1,p2);

  if (SingleCompare(p1->middle,p2->middle)) {
    int sw1 = SingleCompare(p1->left,p2->left);
    if (sw1) sw1 = SingleCompare(p1->right,p2->right);
    if (sw1==0) {
      sw1 = SingleCompare(p1->left,p2->right);
      if (sw1) sw1 = SingleCompare(p1->right,p2->left);
    }
    return sw1;
  }

  if (SingleCompare(p1->middle,p2->left)) {
    int sw1 = SingleCompare(p1->left,p2->middle);
    if (sw1) sw1 = SingleCompare(p1->right,p2->right);
    if (sw1==0) {
      sw1 = SingleCompare(p1->left,p2->right);
      if (sw1) sw1 = SingleCompare(p1->right,p2->middle);
    }
    return sw1;
  }

  if (SingleCompare(p1->middle,p2->right)) {
    int sw1 = SingleCompare(p1->right,p2->middle);
    if (sw1) sw1 = SingleCompare(p1->left,p2->left);
    if (sw1==0) {
      sw1 = SingleCompare(p1->right,p2->left);
      if (sw1) sw1 = SingleCompare(p1->left,p2->middle);
    }
    return sw1;
  }  
  return 0;
}

int Amplitude_Generator::Compare5Vertex(Point* p1,Point* p2)
{
  Point** pts1=new Point*[4];
  Point** pts2=new Point*[4];
  Point *p41,*p42;
  if (p1->left->fl.Kfcode()==kf_shgluon) {
    pts1[0]=p1->left->left;
    pts1[1]=p1->left->middle;
    pts1[2]=p1->left->right;
    pts1[3]=p1->right;
    p41=p1->left;
  }
  else {
    pts1[0]=p1->left;
    pts1[1]=p1->right->left;
    pts1[2]=p1->right->middle;
    pts1[3]=p1->right->right;
    p41=p1->right;
  }
  if (p2->left->fl.Kfcode()==kf_shgluon) {
    pts2[0]=p2->left->left;
    pts2[1]=p2->left->middle;
    pts2[2]=p2->left->right;
    pts2[3]=p2->right;
    p42=p2->left;
  }
  else {
    pts2[0]=p2->left;
    pts2[1]=p2->right->left;
    pts2[2]=p2->right->middle;
    pts2[3]=p2->right->right;
    p42=p2->right;
  }

  if (!CompareColors(p41,p42)) {
    delete[] pts1;
    delete[] pts2;
    return 0;
  }
  
  int hit = 0;
  Permutation perm(4);
  for (int i=0;i<perm.MaxNumber()&&!hit;i++) {
    int* pp=perm.Get(i);
    hit = 1;
    for (int j=0;j<4&&hit;j++) hit = SingleCompare(pts1[j],pts2[pp[j]]);
  }

  delete[] pts1;
  delete[] pts2;
  return hit;
}

int Amplitude_Generator::Kill_Off(Single_Amplitude* &first)
{
  Single_Amplitude* last;
  last = first;
  Single_Amplitude* f1 = first;
  Single_Amplitude* f2;
  int count=0;
  while (f1) {
    if (f1->on==0) {
      ++count;
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
  return count;
}				     


void Amplitude_Generator::CountOrders(Single_Amplitude * & first)
{
  Single_Amplitude* last;
  last = first;
  Single_Amplitude* f1 = first;
  Single_Amplitude* f2;
  int count=0;
  std::vector<double> hitmax;
  while (f1) {
    std::vector<int> hit;
    f1->GetPointlist()->FindOrder(hit);
    bool valid(true);
    for (size_t i(0);i<Min(order.size(),hit.size());++i)
      if (hit[i]>order[i]) valid=false;
    msg_Debugging()<<"Order check: "<<hit<<" vs. "<<order<<" -> "<<valid<<"\n";
    if (hit.size()>hitmax.size())
      hitmax.resize(hit.size(),0);
    for (size_t i(0);i<Min(hit.size(),hitmax.size());++i)
      if (order.size()<=i || hit[i]<=order[i])
	hitmax[i]=Max(hitmax[i],(double)hit[i]);
    if (!valid ||
	!CheckTChannels(f1->GetPointlist())) {
      ++count;
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
  order=hitmax;
  msg_Tracking()<<"Kicked number of diagrams (Amplitude_Generator::CountOrders()) "<<count<<endl;
}

bool Amplitude_Generator::CheckTChannels(Point * p) {
  int ntchan(-1);
  if (p->fl.Strong()) ntchan++;
  p->CountT(ntchan,0);
  //to treat decays correctly
  if (ntchan==-1) ntchan++;
  //
  msg_Debugging()<<METHOD<<" yields "<<ntchan<<" t-channel props, "
		 <<"("<<ntchan_min<<"), start = "<<p->fl<<"."<<std::endl;
  if (ntchan>=ntchan_min) return true;
  return false;
}
 
bool Amplitude_Generator::CheckOrders(Point * p)
{
  std::vector<int> hit;
  p->FindOrder(hit);
  bool valid(true);
  for (size_t i(0);i<Min(order.size(),hit.size());++i)
    if (hit[i]>order[i]) valid=false;
  msg_Debugging()<<"Order check: "<<hit<<" vs. "<<order<<" -> "<<valid<<"\n";
  return valid;
}
 
void Amplitude_Generator::Compare(Single_Amplitude* &first)
{
  Single_Amplitude* f1;
  Single_Amplitude* f2;

  // count amplitudes
  int nampl=0;
  f1 = first;
  while (f1) {
    ++nampl;
    f1=f1->Next;
  } 

  // determine new start_amplitude and 
  f1 = first;
  // start comparison
  Point* p1;
  Point* p2;  

  while (f1) { 
    if (f1->on) {
      p1 = f1->GetPointlist();
      f2 = f1->Next;
      while (f2) {
	if (f2->on) {
	  p2 = f2->GetPointlist();
	  //++ncomps;
	  if (SingleCompare(p1,p2)) {
	    //++noffs;
	    f2->on = 0;
	  }
	}
	f2 = f2->Next;
      }
    }
    f1 = f1->Next;
  }
  
  Kill_Off(first);

  f1 = first;
  nampl=0;
  while (f1) {
    ++nampl;
    f1=f1->Next;
  } 
}

//==================================================================

Point* Amplitude_Generator::FindNext(Point* p)
{
  if (p==0) return 0;
  if (p->left->m==0) return p;
  if (p->right->m==0) return p;
  if ((p->middle!=0) && (p->middle->m==0)) return p;
  
  if (FindNext(p->left))   return p;
  if (FindNext(p->right))  return p;
  if (FindNext(p->middle)) return p;
  return 0;
}

int Amplitude_Generator::ShrinkProps(Point*& p,Point*& pnext, Point*& pcopy, Point*& beg_pcopy,
				     vector<Point*>& pcollist)
{
  
  if (p->left==0 || pnext->left==0) return 0;
  if (p->v->nleg==4 || pnext->v->nleg==4) return 0;
  
  if (pnext->m==1) return 0;
  
  if (pnext->fl.IsFermion()) return 0;
    
  int hit = 0;
  int in[4] = {0};
    
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
  Vertex* v = p_model->p_vertex;
  
  for (short int i=0;i<v->MaxNumber4();i++) {
      
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
	  
	  hit = 1;
	  
	  pcopy->v      = (*v)(i);
	  pcopy->cpl.clear();
	  for (size_t k=0;k<(*v)(i)->cpl.size();k++) pcopy->cpl.push_back((*v)(i)->Coupling(k));
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
	  pnext->m  = 1;

	  if ((*v)(i)->Color.size()==1) {
	    *pcopy->Color = (*v)(i)->Color.back();
            if (pcopy->Lorentz) delete pcopy->Lorentz;
	    pcopy->Lorentz = (*v)(i)->Lorentz.front()->GetCopy();
	    pcopy->t = (*v)(i)->t;
            break;
          }
          else {
	    for (size_t k=0;k<(*v)(i)->Color.size();k++) {
	      *pcopy->Color = (*v)(i)->Color[k];
              if (pcopy->Lorentz) delete pcopy->Lorentz;
	      pcopy->Lorentz = (*v)(i)->Lorentz[k]->GetCopy();
	      pcopy->t = (*v)(i)->t;

	      Color_Function* cfmemo = pcopy->Color;
             
	      //set the contraction indices: 4 -> (prop->number) 
	      while (cfmemo) {
		cfmemo->Replace(4,pnext->number);
		cfmemo = cfmemo->Next();
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
  
  while (f1) { 
    if (f1->on) {
      p = f1->GetPointlist();
      
      for (int i=0; i<dep; i++) p[i].m = 0; 
    
      while (p) {
	//p initial Pointlist with m flags set, pcopy is the contracted Pointlist !!!
	int ll = 0;
	top->Copy(p,pcopy,ll);
	if (EvalPointlist(p,pcopy,pcopy,pcollist)) {
	  if (pcollist.size()==0) {
	    
	    extraAmpl = new Single_Amplitude(pcopy,f1->topnum,f1->permnum,b,dep,N,top,BS,fl,shand);
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
	    for (size_t j=0;j<pcollist.size();j++) {
	      
	      extraAmpl = new Single_Amplitude(pcollist[j],f1->topnum,f1->permnum,b,dep,N,top,BS,fl,shand);
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
	    //erase all elements of pcollist
	    pcollist.clear();
	  }
	  p = FindNext(p);
	}
	else break;
      }
    }
    f1 = f1->Next;
  }
  delete [] pcopy;
}

void Amplitude_Generator::Kill5VertexArtefacts(Single_Amplitude* first)
{
  Single_Amplitude* amp=first;
  while (amp) {
    int tcnt = 0;
    if (Is5VertexArtefact(amp->GetPointlist(),tcnt)) amp->on=0;      
    amp=amp->Next;
  }
}

int Amplitude_Generator::Is5VertexArtefact(Point* p, int &tcnt)
{
  if (!p) return 0;
  if (!p->left) return 0;
  switch(p->t) {
  case 0:
    break;
  case 1: 
    if (tcnt!=-1 || p->fl.Kfcode()!=kf_shgluon) return 1;
    tcnt++;
    break;
  case -1:
    if (tcnt!=0) return 1;
    tcnt--;
  }
  if (Is5VertexArtefact(p->left,tcnt)) return 1;
  if (Is5VertexArtefact(p->right,tcnt)) return 1;
  if (Is5VertexArtefact(p->middle,tcnt)) return 1;
  return 0;
}

int Amplitude_Generator::Count4G(Point * p) {
  if (!p) return 0;
  int v4=0;
  v4+=Count4G(p->left);
  v4+=Count4G(p->right);
  if (p->middle) {
    v4+=Count4G(p->middle);
    Flavour gluon=Flavour(kf_gluon);
    if ((p->fl==gluon)&&(p->left->fl==gluon)&&(p->middle->fl==gluon)&&(p->right->fl==gluon))
      v4+=1;
  }
  return v4;
}

int Amplitude_Generator::CountRealAmplitudes(Single_Amplitude* first) 
{
  Single_Amplitude * f1=first;
  int counts[4]={0,0,0,0};
  while (f1) {
    Point * p =f1->GetPointlist();
    int v4 = Count4G(p);
    if (v4<4) ++(counts[v4]);
    else {
      cerr<<" to many four verticies in one amplitude "<<endl;
    }
    f1=f1->Next; 
  }
  int ra=counts[0];
  int three=3;
  for (int i=1;i<4;++i) {
    ra+=counts[i]/three;
    three*=3;
  } 
  return ra;
}

Single_Amplitude* Amplitude_Generator::Matching()
{
  int nloop = N-1;
  short int i,j;
  int* ii = new int[nloop];
  for (i=0;i<nloop;i++) ii[i] = 0;
  int* perm = new int[N];
  int sw1;
  int qsum,lsum; 
  Single_Amplitude* first_amp;
  int over = 0;
  perm[0]  = 0;
  first_amp = 0;

  int count = 0;
  int start_top = 0;
  int end_top = single_top->number;

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
	if (!((fl[perm[N-1]].IsQuark() && 
	       (b[perm[N-1]]==-1 && (fl[perm[N-1]].IsAnti() || fl[perm[N-1]].Majorana()))) ||
	      (b[perm[N-1]]==1  && (!(fl[perm[N-1]].IsAnti()) || fl[perm[N-1]].Majorana())) )) sw1 = 0;
      }
      if (fl[0].IsQuark() && fl[0].IsAnti()) {
	if (!((fl[perm[1]].IsQuark() && 
	       (b[perm[1]]==-1 && (!fl[perm[1]].IsAnti() || fl[perm[1]].Majorana()))) ||
	      (b[perm[1]]==1  && ((fl[perm[1]].IsAnti()) || fl[perm[1]].Majorana())) )) sw1 = 0;
      }
      if (fl[0].IsLepton() && !fl[0].IsAnti()) {
	if (!((fl[perm[N-1]].IsLepton() && 
	       (b[perm[N-1]]==-1 && (fl[perm[N-1]].IsAnti() || fl[perm[N-1]].Majorana()))) ||
	      (b[perm[N-1]]==1  && (!(fl[perm[N-1]].IsAnti()) || fl[perm[N-1]].Majorana())) )) sw1 = 0;
      }
      if (fl[0].IsLepton() && fl[0].IsAnti()) {
   	if (!((fl[perm[1]].IsLepton() && 
	       (b[perm[1]]==-1 && (!fl[perm[1]].IsAnti() || fl[perm[1]].Majorana()))) ||
   	      (b[perm[1]]==1  && ((fl[perm[1]].IsAnti()) || fl[perm[1]].Majorana())) )) sw1 = 0;
      }
    }

    if (sw1) {
      qsum=lsum=0;
      for (j=0;j<N-1;j++) {
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
 	    if (b[perm[j]]==1) lsum+=LeptonNumber(fl[perm[j]]);
 	    else lsum-=LeptonNumber(fl[perm[j]]);
	}

	if (!fl[0].IsAnti()) {
	  if ( b[perm[N-1]]==-1 && fl[perm[N-1]].IsAnti() && lsum>0 ) 
	      {sw1=0;break;} // s-channel
	  if ( b[perm[N-1]]==1 && !(fl[perm[N-1]].IsAnti()) && lsum>0 )
	    {sw1=0;break;} // t-channel,
	}
	else {
	  if ( b[perm[1]]==-1 && fl[perm[1]].IsAnti() && lsum<0 ) 
	    {sw1=0;break;} // s-channel
	  if ( b[perm[1]]==1 && !(fl[perm[1]].IsAnti()) && lsum<0 )
	    {sw1=0;break;} // t-channel,	    
	}
      }
    }
    if (sw1) {
      for (j=start_top;j<end_top;j++) {
	perm++;
	int pnum = 100;
	Set_End(&single_top->p[j][0],perm,pnum);
	perm -= N;
	single_top->p[j][0].number = *perm;
	single_top->p[j][0].fl     = fl[*perm];
	single_top->p[j][0].b      = b[*perm];

	SetProps(single_top->p[j],2*N-3,first_amp,perm,j,count);
	
	//msg_Out()<<"-----------------------------------"<<std::endl;
	//msg_Out()<<"  "<<single_top->p[j][0].fl<<"("<<single_top->p[j][0].b<<")"<<endl;
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
  
  CreateSingleAmplitudes(first_amp);
  
  CountOrders(first_amp);
  
  if (m_create_4V) CheckFor4Vertices(first_amp);
  Kill5VertexArtefacts(first_amp);
  Kill_Off(first_amp);
  Compare(first_amp);

  Single_Amplitude*  f1 = first_amp;
  while (f1) { 
    (f1->GetPointlist())->ResetFlag();
    f1 = f1->Next;
  }

  return first_amp;
}



















