#include "Basic_Func.H"
#include "Run_Parameter.H"
#include "Message.H"

using namespace AMEGIC;
using namespace AORGTOOLS;
using namespace std;

void Basic_Func::SetArgCouplProp(int narg,Argument* _arg,Complex* _coupl,
				 int pn,Argument* _ps,list<Pfunc*>* _pl) 
{
  arg = _arg; coupl = _coupl;ps = _ps;pl = _pl;
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















