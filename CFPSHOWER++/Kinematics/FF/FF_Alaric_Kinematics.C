#include "CFPSHOWER++/Kinematics/Kinematics_Base.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "CFPSHOWER++/Tools/Parton.H"
#include "PHASIC++/Channels/Antenna_Kinematics.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Message.H"

namespace CFPSHOWER {
  class FF_Alaric_Kinematics : public Kinematics_Base {
  private:
    double m_QTilde2, m_x;
    ATOOLS::Vec4D m_KTilde;
    ATOOLS::Vec4D Nvector(Splitting & split, Configuration & config);
    void   FixKTilde(Splitting & split, Configuration & config);
    bool   KinCheck(Splitting & split);
    double CalculateY(Splitting & split);
    double CalculateX(Splitting & split);
    double CalculateKT2(Splitting & split);
    void   FillConfig(Splitting & split, Configuration & config,
		      PHASIC::Ant_Args & kinargs);
  public:
    FF_Alaric_Kinematics(const Kernel_Info & info);
    ~FF_Alaric_Kinematics();
    
    bool Init(Splitting & split, Configuration & config,
	      const ATOOLS::Mass_Selector * msel);
    bool operator()(Splitting & split, Configuration & config);
    void CalculateJacobean(Splitting & split,Configuration & config);
    bool UpdateSystem(Splitting & split, Configuration & config);
  };
}

using namespace CFPSHOWER;
using namespace ATOOLS;

FF_Alaric_Kinematics::FF_Alaric_Kinematics(const Kernel_Info & info)  :
  Kinematics_Base(info)
{
  m_name   = std::string("FF: Alaric Kinematics");
  m_scheme = kin_type::Alaric;
}

bool FF_Alaric_Kinematics::Init(Splitting & split, Configuration & config,
				const Mass_Selector * msel) {
  FixKTilde(split,config);
  double kappa = m_KTilde.Abs2()/m_QTilde2;
  split.SetIntPhi(1.);
  //Max( 1.,
  //		       2./m_QTilde2*split.GetSpectator()->Mom()*
  //		       (m_KTilde+(1.-kappa)*split.GetSplitter()->Mom()) ));
  split.SetZmin(0.);
  split.SetZmax(1.-sqrt(kappa*split.Tcut()/m_QTilde2));
  InitSystem(split,msel);
  return true;
}

void FF_Alaric_Kinematics::FixKTilde(Splitting & split, Configuration & config)
{
  m_KTilde = Vec4D(0.,0.,0.,0.);
  for (Parton_List::iterator pit=config.begin();pit!=config.end();pit++) {
    if ((*pit)->Mom()[0]<0.) m_KTilde -= (*pit)->Mom();
  }
  m_QTilde2 = 2.*split.GetSplitter()->Mom()*m_KTilde;
}

Vec4D FF_Alaric_Kinematics::Nvector(Splitting & split, Configuration & config)
{
  return m_KTilde+(1.-split.Z())*m_psplit;
}

bool FF_Alaric_Kinematics::
operator()(Splitting & split, Configuration & config) {
  if (!KinCheck(split)) return false;
  PHASIC::Ant_Args kinargs(m_y,m_x,split.Phi());
  FillConfig(split,config,kinargs);
  if (ConstructFSAntenna(split.Mass2(0), split.Mass2(1), m_mspect2,
			 m_psplit, m_pspect, kinargs)<0) return false;
  split.SetMom(0,kinargs.m_pi);
  split.SetMom(1,kinargs.m_pj);
  split.SetMom(2,kinargs.m_pk);
  split.SetKinSpect(Nvector(split,config));
  split.SetAllMoms(kinargs.m_p);
  split.SetKT2(CalculateKT2(split));
  //msg_Out()<<METHOD<<"(pi = "<<kinargs.m_pi<<", pj = "<<kinargs.m_pj<<", "
  //	   <<"Kt = "<<m_KTilde<<")\n"
  //	   <<"-------------------------------------------------------------\n";
  return (split.KT2()>split.Tcut());
}

bool FF_Alaric_Kinematics::KinCheck(Splitting & split) {
  if (split.Q2()<sqr(split.Mass(0)+split.Mass(1)+m_mspect)) return false;
  split.SetY(m_y = CalculateY(split));
  m_x = CalculateX(split);
  return (m_y>=0.0 && m_y<=1.0);
}

double FF_Alaric_Kinematics::CalculateX(Splitting & split) {
  return split.Z()/(1.+m_y);
}

double FF_Alaric_Kinematics::CalculateY(Splitting & split) {
  return split.T()/(m_QTilde2*split.Z()*(1.-split.Z()) - split.T());  
}

double FF_Alaric_Kinematics::CalculateKT2(Splitting & split) {
  return 2.*((split.Mom(0)*split.Mom(1)) * (split.Mom(1)*split.Mom(2))/
	     (split.Mom(0)*split.Mom(2)));
}

void FF_Alaric_Kinematics::
FillConfig(Splitting & split, Configuration & config,PHASIC::Ant_Args & kinargs) {
  kinargs.m_p.reserve(config.size());
  kinargs.m_b.reserve(config.size());
  for (Parton_List::iterator pit=config.begin();pit!=config.end();pit++) {
    if ((*pit)==split.GetSplitter())  kinargs.m_ij = kinargs.m_p.size();
    if ((*pit)==split.GetSpectator()) kinargs.m_k  = kinargs.m_p.size();
    kinargs.m_p.push_back((*pit)->Mom());
    kinargs.m_b.push_back(((*pit)->Flav().Strong() ? 1 : 2) |
			  ((*pit)->Mom()[0]<0.     ? 4 : 0));
  }
}

void FF_Alaric_Kinematics::
CalculateJacobean(Splitting & split,Configuration & config) {
  m_weight = 1.+split.Y();
}


bool FF_Alaric_Kinematics::
UpdateSystem(Splitting & split, Configuration & config) {
  Parton_List::iterator pit=config.begin();
  for (size_t i=0;i<split.NumberP();i++,pit++) {
    (*pit)->SetMom(split.P(i));
  }
  split.GetSpectator()->SetMom(split.Mom(2));
  return true;
}

DECLARE_GETTER(FF_Alaric_Kinematics,"Kinematics_FF_Alaric",
	       Kinematics_Base,Kernel_Info);

Kinematics_Base * ATOOLS::Getter<Kinematics_Base,Kernel_Info,FF_Alaric_Kinematics>::
operator()(const Parameter_Type & info) const
{
  if (info.Type()==kernel_type::FF &&
      info.KinType()==kin_type::Alaric) {
    return new FF_Alaric_Kinematics(info);
  }
  return NULL;
}

void ATOOLS::Getter<Kinematics_Base,Kernel_Info,FF_Alaric_Kinematics>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"Alaric Kinematics (FF)";
}


