#include "HADRONS++/PS_Library/Two_Body_PSs.H"
#include "PHASIC++/Channels/Channel_Elements.H"
#include "PHASIC++/Channels/Channel_Basics.H"
#include "ATOOLS/Org/Message.H"

using namespace HADRONS; 
using namespace PHASIC; 
using namespace ATOOLS; 
using namespace std; 

Iso2Channel::Iso2Channel(const ATOOLS::Flavour * fl) :
  Single_Channel(1,2,fl),
  m_decvec(Vec4D(fl[0].HadMass(),0.,0.,0.))
{
  for (short int i=0;i<nin+nout;i++) ms[i] = ATOOLS::sqr(fl[i].HadMass());
														// get masses^2
  msg_Tracking()<<"Init Iso2Channel("
	   <<fl[0]<<"->"<<fl[1]<<" "<<fl[2]<<", "
	   <<ms[0]<<"->"<<ms[1]<<" "<<ms[2]<<")"<<endl;
  rannum = 2;
  rans   = new double[rannum];
}


void Iso2Channel::GeneratePoint(ATOOLS::Vec4D * p,PHASIC::Cut_Data *,double * _ran)
{
  CE.Isotropic2Momenta(p[0],ms[1],ms[2],p[1],p[2],_ran[0],_ran[1]);
}


void Iso2Channel::GenerateWeight(ATOOLS::Vec4D * p,PHASIC::Cut_Data *)
{
  weight = 1. / ( CE.Isotropic2Weight(p[1],p[2]) * pow(2.*M_PI,2.*3.-4.) );
}


Iso1Channel::Iso1Channel(const ATOOLS::Flavour * fl) :
  Single_Channel(1,1,fl)
{
  msg_Tracking()<<"Init Iso1Channel("
           <<fl[0]<<"->"<<fl[1]<<endl;
}


void Iso1Channel::GeneratePoint(ATOOLS::Vec4D * p,PHASIC::Cut_Data *,double * _ran)
{
  p[1]=p[0];
}


void Iso1Channel::GenerateWeight(ATOOLS::Vec4D * p,PHASIC::Cut_Data *)
{
  weight = 1.0;
}

