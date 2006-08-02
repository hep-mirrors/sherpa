#include "HD_ME_Base.H"
#include "Message.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE HADRONS::HD_ME_Base
#define PARAMETER_TYPE HADRONS::Flavour_Info
#include "Getter_Function.C"

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

HD_ME_Base::HD_ME_Base(int _nout,Flavour * _flavs) :
  m_nout(_nout), p_flavs(_flavs)
{
  p_masses  = new double[1+m_nout];
  p_masses2 = new double[1+m_nout];
  for (int i=0;i<1+m_nout;i++) {
    p_masses[i]  = p_flavs[i].PSMass();
    p_masses2[i] = p_masses[i]*p_masses[i];
  }
  m_metype = "not named yet";
}

HD_ME_Base::~HD_ME_Base()
{
  if (p_masses)  { delete [] p_masses;  p_masses  = NULL; }
  if (p_masses2) { delete [] p_masses2; p_masses2 = NULL; }
}



Isotropic::Isotropic(int _nout,Flavour * _flavs,string _met) :
  HD_ME_Base(_nout,_flavs) 
{
  m_metype = "Isotropic";
  if (msg.LevelIsDebugging()) {
    msg.Tracking()<<"Initialised Isotropic("<<m_nout<<") for "<<endl
	     <<"   "<<p_flavs[0]<<" ->";
    for (int i=1;i<1+m_nout;i++) msg.Tracking()<<" "<<p_flavs[i];
    msg.Tracking()<<endl;
  }
}

DECLARE_GETTER(Isotropic_Getter,"Isotropic",HD_ME_Base,Flavour_Info);
HD_ME_Base* Isotropic_Getter::operator()(const Flavour_Info &parameters) const
{ return new Isotropic(parameters.nout, parameters.flavs); }
void Isotropic_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ str<<"implement me"; }
