#include "AMEGIC++/Amplitude/Zfunctions/Basic_Func.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "AMEGIC++/Amplitude/Pfunc.H"
#include <cassert>

using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;

Basic_Func::~Basic_Func()
{
}

void Basic_Func::SetArgCouplProp(int narg,int* _arg,Complex* _coupl,
				 int _pn,Argument* _ps,Pfunc_List* _pl) 
{
  arg = _arg; coupl = _coupl;ps = _ps;pl = _pl;pn= _pn;
  for (short int i=0;i<narg;i+=2) {
    Map(arg[i]);
  }
  for (short int i=0;i<pn;i++) Map(ps[i].numb,ps[i].maped);
}

void Basic_Func::Map(int& numb) 
{
  if (iabs(numb)>99) {
    Pfunc* p = NULL;
    for (Pfunc_Iterator pit=pl->begin();pit!=pl->end();++pit) {
      p = *pit;
      if (p->arg[0]==ATOOLS::iabs(numb)) break;
    }
    assert(p!=NULL);
    numb = (numb>0) ? p->momnum:-p->momnum;
  }
}

void Basic_Func::Map(int& numb,bool& maped) 
{
  if(maped) return;
  maped=true;

  if (numb<0) {
    msg_Error()<<"Negative Number in Basic_Func::Map() -> numb = "<<numb<<endl;
    abort();
  }
  if (numb>99) {
    Pfunc* p = NULL;
    for (Pfunc_Iterator pit=pl->begin();pit!=pl->end();++pit) {
      p = *pit;
      if (p->arg[0]==numb) break;
    }
    assert(p!=NULL);
    numb = p->momnum;
  }
}

double Basic_Func::GetPMass(int a,int sign)
{
  if (sign!=mt::p_s)return 0.;
  int b;
  for(b=0;b<pn;b++)if(ps[b].numb==ATOOLS::iabs(a))break;
  Pfunc* p1 = NULL;
  int hit = 0;
  for (Pfunc_Iterator pit=pl->begin();pit!=pl->end();++pit) {
    p1 = *pit;
    if (p1->momnum==ATOOLS::iabs(a) && (p1->fl).Kfcode()==ps[b].kfcode) {
      hit = 1;
      break;
    }
  }
  if(hit) return (p1->fl).Mass();
  msg_Error()<<"Basic_Func::GetPMass: Propagator not found! "<<a<<","<<b<<endl
		     <<ps[0].numb<<"."<<ps[1].numb<<"."<<pn<<endl;
  abort();
  return 0.;
}

Kabbala Basic_Func::X(const int,const int,const int) {
  std::cerr<<"Calling Basic_Func::X3"<<std::endl;
  return Kabbala();
}

Kabbala Basic_Func::V(const int a,const int b)
{
  std::cerr<<"Calling Basic_Func::V"<<std::endl;
  return Kabbala();
}

Kabbala Basic_Func::Vcplx(const int a,const int b,const int s)
{
  std::cerr<<"Calling Basic_Func::Vcplx"<<std::endl;
  return Kabbala();
}













