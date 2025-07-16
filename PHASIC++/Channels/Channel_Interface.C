#include "PHASIC++/Channels/Channel_Interface.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace PHASIC;

Channel_Interface::Channel_Interface(int nin,int nout,ATOOLS::Flavour *flavour,ATOOLS::Flavour res):
  Single_Channel(nin,nout,flavour)
{  
  if (nout != 2 || nin!=2) {
    msg_Error()<<"Channel_Interface::Channel_Interface(..): "
		       <<"Cannot handle "<<nin<<" -> "<<nout<<" processes. Abort."<<std::endl;
    exit(169);
  }
  m_rannum = 3;
  delete p_rans;
  p_rans = new double[m_rannum];
  s = smax = pt2max = ATOOLS::sqr(ATOOLS::rpa->gen.Ecms());
  pt2min = 0.;
  E = 0.5*sqrt(s);
  m_name = "Channel Interface";
  mass = width = 0.; 
  type = 0;
  if (res!=ATOOLS::Flavour(kf_none)) {
    mass = res.Mass(); width = res.Width(); type = 1;
  }
}

void Channel_Interface::GeneratePoint(ATOOLS::Vec4D * p,double *ran=0) 
{
  msg_Error()<<"Channel_Interface::GeneratePoint(): Virtual method called!"<<std::endl;
}

void Channel_Interface::GenerateWeight(ATOOLS::Vec4D * p)
{
  msg_Error()<<"Channel_Interface::GenerateWeight(): Virtual method called!"<<std::endl;
}

void Channel_Interface::ISRInfo(int &_type,double &_mass,double &_width) 
{
  _type = type; _mass = mass; _width = width;
}

int Channel_Interface::ChNumber()
{ 
  return chnumber;      
}
