#include "Channel_Interface.H"
#include "Run_Parameter.H"
#include "Message.H"

using namespace PHASIC;

Channel_Interface::Channel_Interface(int _nin,int _nout,ATOOLS::Flavour * _fl,ATOOLS::Flavour _res) 
{  
  if (_nout != 2 || _nin!=2) {
    ATOOLS::msg.Error()<<"Tried to initialize Channel_Interface with "<<_nin<<" -> "<<_nout<<std::endl;
    abort();
  }
  nin  = _nin; nout = _nout;
  ms   = new double[nin+nout];
  for (short int i=0;i<nin+nout;i++) ms[i] = ATOOLS::sqr(_fl[i].Mass());
  rannum = 3;
  rans   = new double[rannum];

  s      = smax  = pt2max = ATOOLS::sqr(ATOOLS::rpa.gen.Ecms());
  pt2min = 0.;
  E      = 0.5 * sqrt(s);
  name   = "Channel Interface";

  mass   = width = 0.; 
  type   = 0;
  if (_res!=ATOOLS::Flavour(ATOOLS::kf::none)) {
    mass = _res.Mass(); width = _res.Width(); type = 1;
  }
}

void Channel_Interface::GeneratePoint(ATOOLS::Vec4D * p,double * _ran=0) {
  ATOOLS::msg.Error()<<"Channel_Interface::GeneratePoint(): Virtual method called!"<<std::endl;
}

void Channel_Interface::GenerateWeight(ATOOLS::Vec4D * p) {
  ATOOLS::msg.Error()<<"Channel_Interface::GenerateWeight(): Virtual method called!"<<std::endl;
}

void Channel_Interface::ISRInfo(int & _type,double & _mass,double & _width) {
  _type = type; _mass = mass; _width = width;
}

int Channel_Interface::ChNumber()
{ 
  return chnumber;      
}
