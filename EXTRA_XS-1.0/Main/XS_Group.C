#include "XS_Group.H"
#include "Run_Parameter.H"
#include "Message.H"
#include "Random.H"
#include "MathTools.H"
#include "Process_Group.H"

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>

using namespace EXTRAXS;
using namespace PHASIC;
using namespace AMEGIC;
using namespace AORGTOOLS;
using namespace AMATOOLS;
using namespace APHYTOOLS;
using namespace std;

XS_Group::XS_Group(int _nin,int _nout, string _name) 
{
  Init(_nin,_nout,0);
  name = _name; 
}

XS_Group::~XS_Group()
{
  for(int i=xsecs.size();i>0;i--) {
    if (xsecs[i-1]) delete xsecs[i-1];
  }
}

void XS_Group::Add(XS_Base * _xsec, bool createbroker) 
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
  if (createbroker) _xsec->MakeBroker(broker->ISR(),broker->Beam(),seldata,broker);
  xsecs.push_back(_xsec);
};

void XS_Group::MakeBroker(ISR::ISR_Handler * isr, BEAM::Beam_Handler * beam,
			  APHYTOOLS::Selector_Data * _seldata, AMEGIC::Process_Group * _broker) 
{
  seldata = _seldata;
  Pol_Info * _plavs = 0;
  Flavour * _fl = new Flavour[Nin()+Nout()];
  for (int i=0; i<(Nin()+Nout()); i++) _fl[i] = Flavs()[i];
  broker = new Process_Group(Nin(),Nout(),_fl,isr,beam,seldata,2,
			     rpa.me.KFactorScheme(),rpa.me.ScaleScheme(),
			     _plavs,AMEGIC::FORCED_MODE);
  _broker->Add(broker);
  delete [] _fl;
}
