#include "AMISIC++/Tools/Fixed_Variable_Channel.H"

#include "PHASIC++/Channels/Channel_Elements.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Math/Random.H"

using namespace PHASIC;

Fixed_Variable_Channel::
Fixed_Variable_Channel(int nin,int nout,ATOOLS::Flavour *fl,
		       const std::string &variable):
  PHASIC::Channel_Interface(nin,nout,fl),
  p_variable(ATOOLS::Variable_Getter::GetObject(variable,""))
{
  if (ms[2]!=ms[3])
    THROW(not_implemented, "MPI phase space only implemented for m2=m3.");
}

void Fixed_Variable_Channel::
GeneratePoint(ATOOLS::Vec4D *p,double *ran)
{
  if (p_variable->SelectorID()=="PT") {
    Ehat=sqrt((p[0]+p[1]).Abs2());
    pt=m_value;
    double mt2=pt*pt+ms[2];
    if (Ehat/2.0>sqrt(mt2)) {
      weight=1.0/sqrt((Ehat*Ehat)/(4.*pt*pt)-1.0);
      double pz(sqrt(Ehat*Ehat/4.0-mt2));
      if (ATOOLS::ran->Get()<0.5) pz=-pz;
      p[2]=ATOOLS::Vec4D(Ehat/2.0,pt*cos(2.0*M_PI*ran[1]),
                         pt*sin(2.0*M_PI*ran[1]),pz);
      m_trigger=true;
    }
    else {
      weight=ATOOLS::Accu();
      m_trigger=false;
      p[2]=ATOOLS::Vec4D(Ehat/2.0,Ehat/2.0*cos(2.0*M_PI*ran[1]),
			 Ehat/2.0*sin(2.0*M_PI*ran[1]),0.0);
    }
    p[3]=ATOOLS::Vec4D(Ehat/2.0,ATOOLS::Vec3D()-ATOOLS::Vec3D(p[2]));
    return;
  }
  msg_Error()<<"Fixed_Variable_Channel::GeneratePoint(..): "
	     <<"Cannot handle "<<p_variable->Name()
	     <<"! Setting weight to 0."<<std::endl;
  weight=0.0;
}
  
void Fixed_Variable_Channel::GenerateWeight(ATOOLS::Vec4D *_p)
{
  weight/=PHASIC::CE.Isotropic2Weight(_p[2],_p[3])*
    pow(2.0*M_PI,2.0*3.0-4.0);
}

std::string Fixed_Variable_Channel::ChID()
{
  return m_chid;
}
