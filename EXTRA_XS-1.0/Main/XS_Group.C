#include "XS_Group.H"
#include "Run_Parameter.H"
#include "Message.H"
#include "Random.H"
#include "MathTools.H"

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>

using namespace EXTRAXS;
using namespace PHASIC;
using namespace AORGTOOLS;
using namespace AMATOOLS;
using namespace APHYTOOLS;
using namespace std;

XS_Group::XS_Group(int _nin,int _nout,string _name) 
{
  Init(_nin,_nout,0);
  name = _name; 
}

XS_Group::~XS_Group()
{
  for(int i=xsecs.size();i>0;i--) {
    if (xsecs[i-1]) delete xsecs[i-1];
  }
  delete broker;
}

void XS_Group::Add(XS_Base * _xsec) 
{
  if (xsecs.size()==0) {
    nin  = _xsec->Nin();
    nout = _xsec->Nout();

    fl  = new Flavour[nin+nout];
    for (short int i=0;i<nin+nout;i++) {
      fl[i]  = _xsec->Flavs()[i];
    }
  }
  else {
    if ( (nin != _xsec->Nin()) || (nout != _xsec->Nout())) {
      AORGTOOLS::msg.Error()<<"Error : Cannot add Process "<<_xsec->Name()
			    <<" to group "<<name<<" ! "<<endl
			    <<"   Inconsistent number of external legs."<<endl 
			    <<"  Before : ("<<nin<<" -> "<<nout<<" )"<<endl
			    <<"  Now    : ("<<_xsec->Nin()<<" -> "<<_xsec->Nout()<<" )"<<endl;
      return;
    }
  }  
  AORGTOOLS::msg.Tracking()<<"Add xs "<<_xsec->Name()<<" to group "<<name<<" ! "<<endl; 
  xsecs.push_back(_xsec);
};

AMEGIC::Process_Base * XS_Group::CreateBroker() {
  broker = new Broker_Group(this);
  for (int i=0; i<xsecs.size(); i++) broker->Add(xsecs[i]);
  return broker;
}
