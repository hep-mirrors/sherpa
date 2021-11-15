#include "CFPSHOWER++/Kinematics/Kinematics_Base.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "ATOOLS/Org/Message.H"

namespace CFPSHOWER {
  class FF_CS_Kinematics : public Kinematics_Base {
  private:    
    bool   KinCheck(Splitting & split);
    double CalculateY(Splitting & split);
  public:
    FF_CS_Kinematics(const Kernel_Info & info);
    ~FF_CS_Kinematics() {}
    bool Init(Splitting & split, Configuration & config,
	      const ATOOLS::Mass_Selector * msel);
    bool operator()(Splitting & split, Configuration & config);
    bool InitLimits(Splitting & split, Configuration & config);
    void CalculateJacobean(Splitting & split, Configuration & config);
    bool UpdateSystem(Splitting & split, Configuration & config);
  };
}

using namespace CFPSHOWER;
using namespace ATOOLS;

FF_CS_Kinematics::FF_CS_Kinematics(const Kernel_Info & info)  :
  Kinematics_Base(info)
{
  m_name   = std::string("FF: CS Kinematics");
  m_scheme = kin_type::CataniSeymour;
}

bool FF_CS_Kinematics::Init(Splitting & split, Configuration & config,
			    const ATOOLS::Mass_Selector * msel) {
  //msg_Out()<<METHOD<<"(tcut = "<<split.Tcut()<<", Q2 = "<<split.Q2()<<")\n";
  double disc = 1.-4.*split.Tcut()/split.Q2();
  if (disc<0.) return false;
  double delta = sqrt(disc);
  split.SetZmin((1.-delta)/2.);
  split.SetZmax((1.+delta)/2.);
  InitSystem(split,msel);
  return true;
}

bool FF_CS_Kinematics::KinCheck(Splitting & split) {
  split.SetY(m_y    = CalculateY(split));
  return (m_y>=0.0 && m_y<=1.0);
}

double FF_CS_Kinematics::CalculateY(Splitting & split) {
  return ( ( split.T() / ((1.-split.Z())*split.Z()) +
	     split.Mass2(0) * split.Z()/(1.-split.Z()) +
	     split.Mass2(1) * (1.-split.Z())/split.Z()) /
	   ( split.Q2()-split.Mass2(0)-split.Mass2(1)-m_mspect2) );
}

bool FF_CS_Kinematics::operator()(Splitting & split, Configuration & config) {
  if (!KinCheck(split)) return false;
  PHASIC::Kin_Args kinargs(split.Y(),split.Z(),split.Phi());
  if (PHASIC::ConstructFFDipole(split.Mass2(0),split.Mass2(1),
				m_msplit2,m_mspect2,
				m_psplit,m_pspect,kinargs) < 0) return false;
  split.SetMom(0,kinargs.m_pi);
  split.SetMom(1,kinargs.m_pj);
  split.SetMom(2,kinargs.m_pk);
  split.SetKinSpect(kinargs.m_pk);
  //msg_Out()<<METHOD
  //	   <<"(z = "<<split.Z()<<" vs "<<m_z[0]<<" from "<<split.Ztest()<<", "
  //	   <<"y = "<<m_y<<", phi = "<<split.Phi()<<"):\n"
  //	   <<m_psplit<<" + "<<m_pspect<<" --> \n"
  //	   <<kinargs.m_pi<<" + "<<kinargs.m_pj<<" + "<<kinargs.m_pk<<"\n";
  return true;
}

void FF_CS_Kinematics::CalculateJacobean(Splitting & split,
					 Configuration & config) {
  m_weight = 1.-CalculateY(split);
}

bool FF_CS_Kinematics::UpdateSystem(Splitting & split,Configuration & config) {
  Vec4D check = (split.GetSplitter()->Mom()+split.GetSpectator()->Mom());
  for (size_t i=0;i<3;i++) check -= split.Mom(i);
  if (dabs(check.Abs2())>1.e-6 || dabs(check[0])>1.e-3) {
    msg_Error()<<"\n"<<METHOD<<" throws momentum conservation error: "
	       <<"check = "<<check<<"\n";
    return false;
  }
  split.GetSpectator()->SetMom(split.Mom(2));
  return true;
}



DECLARE_GETTER(FF_CS_Kinematics,"Kinematics_FF_Catani-Seymour",
	       Kinematics_Base,Kernel_Info);

Kinematics_Base * ATOOLS::Getter<Kinematics_Base,Kernel_Info,FF_CS_Kinematics>::
operator()(const Parameter_Type & info) const
{
  if (info.Type()==kernel_type::FF &&
      info.KinType()==kin_type::CataniSeymour) {
    return new FF_CS_Kinematics(info);
  }
  return NULL;
}

void ATOOLS::Getter<Kinematics_Base,Kernel_Info,FF_CS_Kinematics>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"Catani-Seymour Kinematics (FF)";
}


