#include "BEAM/Main/EPA.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"

using namespace BEAM;
using namespace ATOOLS;

EPA::EPA(const Flavour _beam,const double _mass,
	 const double _charge,const double _energy,
	 const double _pol,const int _dir,Data_Reader *const read):
  Beam_Base("EPA",_beam,_energy,_pol,_dir),
  m_mass(_mass), m_charge(_charge)
{ 
  m_bunch=Flavour(kf_photon);
  m_vecout=Vec4D(m_energy,0.,0.,_dir*m_energy);
  std::string num(_dir>0?"1":"2");
  m_q2Max=read->GetValue<double>("EPA_q2Max_"+num,2.0);
  m_pt_min=read->GetValue<double>("EPA_ptMin_"+num,0.0);
  m_aqed=read->GetValue<double>("EPA_AlphaQED",0.0072992701);
  if (m_pt_min>1.0) {
    /* pt_min > 1 - according to approximation of 'qmi' calculation in CalculateWeight */
    THROW(critical_error,"Too big p_T cut ( "+ToString(m_pt_min)+")");
  }
}

EPA::~EPA() 
{
}

double EPA::phi(double x, double qq)
{
  const double a = 7.16;
  const double b = -3.96;
  const double c = .028;
  double y,qq1,f;
  qq1=1+qq;
  y= x*x/(1-x);
  f=(1+a*y)*(-log(qq1/qq)+1/qq1+1/(2*qq1*qq1)+1/(3*qq1*qq1*qq1));
  f+=(1-b)*y/(4*qq*qq1*qq1*qq1);
  f+=c*(1+y/4)*(log((qq1-b)/qq1)+b/qq1+b*b/(2*qq1*qq1)+b*b*b/(3*qq1*qq1*qq1));
  return f;
}

bool EPA::CalculateWeight(double x,double q2)
{
  const double alpha = m_aqed;
  const double qz = 0.71;
  m_x = x; m_Q2 = q2;
  if (x==1) {
    m_weight=0.0;
    return 1;
  }
  double f, qmi, qma;
  qma=m_q2Max/qz;
  // x = omega/E = (E-E')/E  ; E,E' - incoming and outgoing protons energy
  //                           omega = E-E' - energy of emitted photon
  qmi= m_mass*m_mass * x*x /(1-x)/qz;
  qmi+=m_pt_min*m_pt_min /(1-x)/qz;
                                                                                
  f = alpha/M_PI*(phi(x,qma)-phi(x,qmi))*(1-x)/x;
  f *= m_charge*m_charge;
  if (f < 0) f = 0.;
  m_weight = f;
  return 1;
}

Beam_Base *EPA::Copy() 
{
  return new EPA(*this);
}

double EPA::Weight(Flavour fl)                
{ 
  return m_weight; 
}

ATOOLS::Vec4D EPA::OutMomentum()              
{ 
  return m_x*m_vecout; 
}

ATOOLS::Flavour EPA::Remnant()                
{ 
  return m_beam;
}



