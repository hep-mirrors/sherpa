#include "PHASIC++/Channels/Decay_Dalitz.H"
#include "PHASIC++/Channels/Channel_Elements.H"
#include "PHASIC++/Channels/Channel_Basics.H"

using namespace PHASIC; 
using namespace ATOOLS; 
using namespace std; 


Decay_Dalitz::Decay_Dalitz(
	const ATOOLS::Flavour * fl,
        const double& mass, const double& width,
	size_t dir, size_t p1, size_t p2,
        const ATOOLS::Mass_Selector* masssel) :
  Single_Channel(1,3,fl),
  m_decvec(Vec4D(fl[0].HadMass(),0.,0.,0.)),
  m_pmass(mass), m_pwidth(width),
  m_sexp(.5),
  m_p1(p1), m_p2(p2), m_dir(dir), m_mode(0),
  p_masssel(masssel)
{
  for (short int i=0;i<m_nin+m_nout;i++) p_ms[i] = p_masssel->Mass2(fl[i]);
  m_smin = ATOOLS::sqr(p_masssel->Mass(fl[m_p1])+p_masssel->Mass(fl[m_p2]));
  m_smax = ATOOLS::sqr(p_masssel->Mass(fl[0])-p_masssel->Mass(fl[m_dir]));
  if (sqrt(m_smin)<m_pmass*10.) m_mode = 1;

  m_rannum = 5;
  p_rans   = new double[m_rannum];
}


void Decay_Dalitz::GeneratePoint(ATOOLS::Vec4D * p,PHASIC::Cut_Data *,double * _ran)
{
  double sprop;
  if (m_mode==1) sprop = CE.MassivePropMomenta(m_pmass,m_pwidth,1,m_smin,m_smax,_ran[0]);
  else sprop = CE.MasslessPropMomenta(m_sexp,m_smin,m_smax,_ran[0]);     
  CE.Isotropic2Momenta(p[0],p_ms[m_dir],sprop,p[m_dir],m_pvec,_ran[1],_ran[2]);
  CE.Isotropic2Momenta(m_pvec,p_ms[m_p1],p_ms[m_p2],p[m_p1],p[m_p2],_ran[3],_ran[4]);
}


void Decay_Dalitz::GenerateWeight(ATOOLS::Vec4D * p,PHASIC::Cut_Data *)
{
  m_weight = 1.;
  double sprop  = (p[m_p1]+p[m_p2]).Abs2();
  if (m_mode==1) 
    m_weight *= CE.MassivePropWeight(m_pmass,m_pwidth,1,m_smin,m_smax,sprop);
  else 
    m_weight *= CE.MasslessPropWeight(m_sexp,m_smin,m_smax,sprop);     
  m_weight   *= CE.Isotropic2Weight(p[m_dir],p[m_p1]+p[m_p2]);
  m_weight   *= CE.Isotropic2Weight(p[m_p1],p[m_p2]);
  m_weight    =  1./(m_weight * pow(2.*M_PI,3.*3.-4.));  
}
