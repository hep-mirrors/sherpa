#include "Basic_Func.H"
#include "Run_Parameter.H"
#include "Message.H"

using namespace AMEGIC;
using namespace AORGTOOLS;
using namespace AMATOOLS;
using namespace std;

void Basic_Func::SetArgCouplProp(int narg,Argument* _arg,Complex* _coupl,
				 int _pn,Argument* _ps,list<Pfunc*>* _pl) 
{
  arg = _arg; coupl = _coupl;ps = _ps;pl = _pl;pn= _pn;
  for (short int i=0;i<narg;i+=2) {
    Map(arg[i].numb,arg[i].maped);
    //changing sign according to spinordirection
    if (arg[i].spinortype==Spinor::v || arg[i].spinortype==Spinor::vbar) {
      arg[i].numb   *= -1; //for mu
      arg[i+1].numb *= -1; //for sign
    }
  }

  for (short int i=0;i<pn;i++) Map(ps[i].numb,ps[i].maped);

/* cout<<"*****SetArgCouplProp nachher: "<<endl;
  for (short int i=0;i<narg;i+=2) {
    cout<<"("<<i<<":"<<arg[i].numb<<"/"<<arg[i].maped<<") ";
  }
  cout<<endl;
  for (short int i=0;i<pn;i++) {
    cout<<"("<<i<<":"<<ps[i].numb<<"/"<<ps[i].maped<<") ";
  }
  cout<<endl;*/
}

void Basic_Func::Map(int& numb,bool& maped) 
{

  if(maped) return;
  maped=true;

  if (numb<0) {
    AORGTOOLS::msg.Error()<<"Negative Number in Basic_Func::Map() -> numb = "<<numb<<endl;
    abort();
  }
  if (numb>99) {
    Pfunc* p;
    for (list<Pfunc*>::iterator pit=pl->begin();pit!=pl->end();++pit) {
      p = *pit;
      if (p->arg[0]==AMATOOLS::iabs(numb)) break;
    }
    numb = (numb>0) ? p->momnum:-p->momnum;
  }
}

double Basic_Func::GetPMass(int a,int sign)
{
  if (sign!=mt::p_s)return 0.;
  int b;
  for(b=0;b<pn;b++)if(ps[b].numb==AMATOOLS::iabs(a))break;
  Pfunc* p1;
  int hit = 0;
  for (list<Pfunc*>::iterator pit=pl->begin();pit!=pl->end();++pit) {
    p1 = *pit;
    if (p1->momnum==AMATOOLS::iabs(a) && (p1->fl).Kfcode()==ps[b].kfcode) {
      hit = 1;
      break;
    }
  }
  //cout<<"Basic_Func::GetPMass: "<<ps[b].numb<<" "<<a<<" Hit:"<<hit<<" "<<p1->fl<<endl;
  if(hit) return (p1->fl).Mass();
  AORGTOOLS::msg.Error()<<"Basic_Func::GetPMass: Propagator not found! "<<a<<","<<b<<endl;
  AORGTOOLS::msg.Error()<<ps[0].numb<<"."<<ps[1].numb<<"."<<pn<<endl;
  abort();
  return 0.;
}














