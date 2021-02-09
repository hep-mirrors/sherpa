#include "CFPSHOWER++/Kinematics/Kinematics_Base.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "ATOOLS/Org/Message.H"

namespace CFPSHOWER {
  class FF_CS_Kinematics : public Kinematics_Base {
  private:
    ATOOLS::Vec4D m_psplit, m_pspect, m_pboth;
    double        m_msplit, m_mspect, m_msplit2, m_mspect2, m_m[2], m_m2[2];
    double        m_Q2, m_Q, m_z[2], m_y;
    bool          m_ismassive;
    
    void Init(const Splitting & split, Configuration & config,
	      const ATOOLS::Mass_Selector * msel);
    bool KinCheck(Splitting & split);
    double CalculateY(Splitting & split);
  public:
    FF_CS_Kinematics(const Kernel_Info & info);
    bool operator()(Splitting & split, Configuration & config);
    void CalculateInvariants(Splitting & split, Configuration & config);
    void CalculateJacobean(Splitting & split, Configuration & config);
    bool UpdateSystem(Splitting & split, Configuration & config);
  };
}

using namespace CFPSHOWER;
using namespace ATOOLS;

FF_CS_Kinematics::FF_CS_Kinematics(const Kernel_Info & info)  :
  Kinematics_Base(info)
{
  SetName("FF: CS Kinematics");
  SetScheme(kin_type::CS);
}

void FF_CS_Kinematics::Init(const Splitting & split, Configuration & config,
			    const ATOOLS::Mass_Selector * msel) {
  m_psplit = split.GetSplitter()->Mom();
  m_pspect = split.GetSpectator()->Mom();
  m_pboth  = m_psplit + m_pspect;
  m_msplit = msel->Mass(split.GetSplitter()->Flav());  m_msplit2 = sqr(m_msplit);
  m_mspect = msel->Mass(split.GetSpectator()->Flav()); m_mspect2 = sqr(m_mspect);
  m_Q2     = m_pboth.Abs2();  m_Q = sqrt(m_Q2);
  m_ismassive = (m_mspect>0. || m_m[0]>0. || m_m[1]>0.);
}

bool FF_CS_Kinematics::KinCheck(Splitting & split) {
  split.SetY(m_y = CalculateY(split));
  return (m_y>=0.0 && m_y<=1.0);
}

double FF_CS_Kinematics::CalculateY(Splitting & split) {
  //msg_Out()<<METHOD<<": t = "<<split.T()<<", z = "<<split.Z()
  //	   <<"["<<split.Zmin()<<", "<<split.Zmax()<<"], "
  //	   <<"and Q2 = "<<split.Q2()
  //	   <<" yields y = "<<(split.T()/((1.-split.Z())*split.Z()*split.Q2()))<<"\n";
  return split.T()/((1.-split.Z())*split.Z()*split.Q2());
}

bool FF_CS_Kinematics::operator()(Splitting & split, Configuration & config) {
  if (!KinCheck(split)) return false;
  PHASIC::Kin_Args kinargs(m_y,split.Z(),split.Phi());
  for (size_t i=0;i<2;i++) m_m2[i] = sqr(split.GetKernel()->GetFlavs()[i].Mass());
  if (PHASIC::ConstructFFDipole(m_m2[0],m_m2[1],m_msplit2,m_mspect2,
				m_psplit,m_pspect,kinargs) < 0) return false;
  split.SetMom(0,kinargs.m_pi);
  split.SetMom(1,kinargs.m_pj);
  split.SetMom(2,kinargs.m_pk);
  return true;
}

void FF_CS_Kinematics::CalculateInvariants(Splitting & split, Configuration & config) {}

void FF_CS_Kinematics::CalculateJacobean(Splitting & split, Configuration & config) {
  m_weight = 1.-m_y;
}

bool FF_CS_Kinematics::UpdateSystem(Splitting & split,Configuration & config) {
  //msg_Out()<<"   "<<split.GetSplitter()->Mom()<<" + "<<split.GetSpectator()->Mom()<<"\n"
  //	   <<"   "<<split.Mom(0)<<" + "<<split.Mom(1)<<" + "<<split.Mom(2)<<"\n";
  Vec4D check = (split.GetSplitter()->Mom()+split.GetSpectator()->Mom());
  //msg_Out()<<"   sum = "<<check<<" --> "<<(check-split.Mom(0)-split.Mom(1)-split.Mom(2))<<"\n";
  for (size_t i=0;i<3;i++) check -= split.Mom(i);
  if (dabs(check.Abs2())>1.e-6 || dabs(check[0])>1.e-3) {
    msg_Error()<<"\n"<<METHOD<<" throws momentum conservation error: check = "<<check<<"\n";
    return false;
  }
  split.GetSpectator()->SetMom(split.Mom(2));
  return true;
}



DECLARE_GETTER(FF_CS_Kinematics,"Kinematics_FF_Catani-Seymour",Kinematics_Base,Kernel_Info);

Kinematics_Base * ATOOLS::Getter<Kinematics_Base,Kernel_Info,FF_CS_Kinematics>::
operator()(const Parameter_Type & info) const
{
  if ((info.Type()==kernel_type::FF &&
       info.KinType()==kin_type::CS) ||
      (info.Type()==kernel_type::FF &&
       info.KinType()==kin_type::PanGlobal &&
       info.LogType()==log_type::coll)) {
    return new FF_CS_Kinematics(info);
  }
  return NULL;
}

void ATOOLS::Getter<Kinematics_Base,Kernel_Info,FF_CS_Kinematics>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"Catani-Seymour Kinematics (FF)";
}


