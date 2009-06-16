#include "HADRONS++/PS_Library/Three_Body_PSs.H"
#include "PHASIC++/Channels/Channel_Elements.H"
#include "PHASIC++/Channels/Channel_Basics.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "HADRONS++/PS_Library/ResonanceFlavour.H"

using namespace HADRONS; 
using namespace PHASIC; 
using namespace ATOOLS; 
using namespace std; 

Dalitz::Dalitz(
	const ATOOLS::Flavour * fl,
	SimpleResonanceFlavour res,
	const int p1, 
	const int p2 ) :
  Single_Channel(1,3,fl),
  m_decvec(Vec4D(fl[0].HadMass(),0.,0.,0.)),
  m_pmass(res.Mass()), m_pwidth(res.Width()), 
  m_p1(p1), m_p2(p2), m_mode(0),
  m_sexp(.5)
{
  name = string("Dalitz_")+res.Name()+string("_")+ToString(m_p1)+ToString(m_p2);
											// generate channel name
  for (short int i=0;i<nin+nout;i++) ms[i] = ATOOLS::sqr(fl[i].HadMass());
											// set masses^2
  msg_Tracking()<<"Init Dalitz("<<name<<" : "
	   <<fl[0]<<"->"<<fl[1]<<" "<<fl[2]<<" "<<fl[3]<<", "
	   <<ms[0]<<"->"<<ms[1]<<" "<<ms[2]<<" "<<ms[3]<<")"<<endl;
  for (int i=1;i<4;i++) {
    if (m_p1!=i && m_p2!=i) { m_dir=i; break; }
  }				// get the one with no resonance
  m_smin = ATOOLS::sqr(fl[m_p1].HadMass()+fl[m_p2].HadMass());
  m_smax = ATOOLS::sqr(fl[0].HadMass()-fl[m_dir].HadMass());
  if (sqrt(m_smin)<m_pmass*10.) m_mode = 1;

  rannum = 5;
  rans   = new double[rannum];
}


void Dalitz::GeneratePoint(ATOOLS::Vec4D * p,PHASIC::Cut_Data *,double * _ran)
{
  /*
  double sprop;
  if (m_mode==1) sprop = CE.MassivePropMomenta(m_pmass,m_pwidth,1,m_smin,m_smax,_ran[0]);
  else sprop = CE.MasslessPropMomenta(m_sexp,m_smin,m_smax,_ran[0]);     
  CE.Isotropic2Momenta(p[0],ms[m_dir],sprop,p[m_dir],m_pvec,_ran[1],_ran[2]);
  CE.Isotropic2Momenta(m_pvec,ms[m_p1],ms[m_p2],p[m_p1],p[m_p2],_ran[3],_ran[4]);
  */
  double sprop;
  if (m_mode==1) sprop = CE.MassivePropMomenta(m_pmass,m_pwidth,1,m_smin,m_smax,_ran[0]);
  else sprop = CE.MasslessPropMomenta(m_sexp,m_smin,m_smax,_ran[0]);     
  CE.Isotropic2Momenta(p[0],ms[m_dir],sprop,p[m_dir],m_pvec,_ran[1],_ran[2]);
  CE.Isotropic2Momenta(m_pvec,ms[m_p1],ms[m_p2],p[m_p1],p[m_p2],_ran[3],_ran[4]);
}


void Dalitz::GeneratePoint(ATOOLS::Vec4D * p,double * _ran)
{
  double sprop;
  if (m_mode==1) sprop = CE.MassivePropMomenta(m_pmass,m_pwidth,1,m_smin,m_smax,_ran[0]);
  else sprop = CE.MasslessPropMomenta(m_sexp,m_smin,m_smax,_ran[0]);     
  CE.Isotropic2Momenta(p[0],ms[m_dir],sprop,p[m_dir],m_pvec,_ran[1],_ran[2]);
  CE.Isotropic2Momenta(m_pvec,ms[m_p1],ms[m_p2],p[m_p1],p[m_p2],_ran[3],_ran[4]);
}


void Dalitz::GenerateWeight(ATOOLS::Vec4D * p,PHASIC::Cut_Data *)
{
  /*
  weight = 1.;
  double sprop  = (p[m_p1]+p[m_p2]).Abs2();
  if (m_mode==1) 
    weight *= CE.MassivePropWeight(m_pmass,m_pwidth,1,m_smin,m_smax,sprop);
  else 
    weight *= CE.MasslessPropWeight(m_sexp,m_smin,m_smax,sprop);     
  weight   *= CE.Isotropic2Weight(p[m_dir],p[m_p1]+p[m_p2]);
  weight   *= CE.Isotropic2Weight(p[m_p1],p[m_p2]);
  weight    =  1./(weight * pow(2.*M_PI,3.*3.-4.));  
  */
  weight = 1.;
  double sprop  = (p[m_p1]+p[m_p2]).Abs2();
  if (m_mode==1) 
    weight *= CE.MassivePropWeight(m_pmass,m_pwidth,1,m_smin,m_smax,sprop);
  else 
    weight *= CE.MasslessPropWeight(m_sexp,m_smin,m_smax,sprop);     
  weight   *= CE.Isotropic2Weight(p[m_dir],p[m_p1]+p[m_p2]);
  weight   *= CE.Isotropic2Weight(p[m_p1],p[m_p2]);
  weight    =  1./(weight * pow(2.*M_PI,3.*3.-4.));  
}


void Dalitz::GenerateWeight(ATOOLS::Vec4D * p)
{
  weight = 1.;
  double sprop  = (p[m_p1]+p[m_p2]).Abs2();
  if (m_mode==1) 
    weight *= CE.MassivePropWeight(m_pmass,m_pwidth,1,m_smin,m_smax,sprop);
  else 
    weight *= CE.MasslessPropWeight(m_sexp,m_smin,m_smax,sprop);     
  weight   *= CE.Isotropic2Weight(p[m_dir],p[m_p1]+p[m_p2]);
  weight   *= CE.Isotropic2Weight(p[m_p1],p[m_p2]);
  weight    =  1./(weight * pow(2.*M_PI,3.*3.-4.));  
}
