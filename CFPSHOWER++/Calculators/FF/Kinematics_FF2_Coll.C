#include "CFPSHOWER++/Calculators/FF/Kinematics_FF2_Coll.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "ATOOLS/Org/Message.H"

using namespace CFPSHOWER;
using namespace ATOOLS;

Kinematics_FF2_Coll::Kinematics_FF2_Coll() :
  Kinematics_Base(log_type::coll, kernel_type::FF)
{
  m_name = "FF2(coll)";
}

Kinematics_FF2_Coll::~Kinematics_FF2_Coll() {}

bool Kinematics_FF2_Coll::operator()(Splitting & split,Configuration & config) {
  Init(split);
  return ConstructSystem(split);
}

bool Kinematics_FF2_Coll::KinCheck(Splitting & split) {
  return (split.y()>=0.0 && split.y()<=1.0);
}

double Kinematics_FF2_Coll::CalculateY(Splitting & split) {
  if (split.GetSplitter()->Flav().IsFermion())
    return split.t(0)/((1.-split.z(0))*split.Q2());
  if (split.GetSplitter()->Flav().IsVector()) {
    if (split.GetKernel()->GetFlavs()[0].IsVector())
      return split.t(0)/((1.-split.z(0))*split.z(0)*split.Q2());
    if (split.GetKernel()->GetFlavs()[0].IsFermion())
      return split.t(0)/split.Q2();
  }
  return split.t(0)/((1.-split.z(0))*split.z(0)*split.Q2());
}

bool Kinematics_FF2_Coll::ConstructSystem(Splitting & split) {
  split.Set_y(CalculateY(split));
  if (!KinCheck(split)) return false;
  PHASIC::Kin_Args kinargs(split.y(),split.z(0),split.phi(0));
  if (PHASIC::ConstructFFDipole(split.m2(0),split.m2(1),
				split.msplit2(),split.mspect2(),
				split.GetSplitter()->Mom(),split.GetSpectator()->Mom(),
				kinargs) < 0) return false;
  split.Set_ztilde(0,split.z(0));
  split.Set_ztilde(1,1.-split.z(0));
  split.SetMomentum(0,m_psplit = kinargs.m_pi);
  split.SetMomentum(1,m_pnew   = kinargs.m_pj);
  split.SetSpectatorMomentum(m_pspect = kinargs.m_pk);
  return true;
}

void Kinematics_FF2_Coll::CalculateWeight(Splitting & split) {
  m_weight = 1.-split.y();
  //
  //double Q2red = split.Q2red(), Q2red_y = Q2red & split.y(); 
  //m_weight = ( Q2red / Lambda(split.Q2(),split.msplit2(),split.mspect2()) *
  //	       Q2red_y /(Q2red_y+split.m2(0)+split.m2(1)-split.msplit2()) );
}

bool Kinematics_FF2_Coll::UpdateSystem(Splitting & split,Configuration & config) {
  split.GetSpectator()->SetMom(split.SpectatorMomentum());
  return true;
}




DECLARE_GETTER(Kinematics_FF2_Coll,"Kinematics_FF2_Coll",Kinematics_Base,Kernel_Info);

Kinematics_Base * ATOOLS::Getter<Kinematics_Base,Kernel_Info,Kinematics_FF2_Coll>::
operator()(const Parameter_Type & info) const
{
  if (info.Type()==kernel_type::FF &&
      info.LogType()==log_type::coll &&
      info.GetFlavs().size()==2) {
    return new Kinematics_FF2_Coll();
  }
  return NULL;
}

void ATOOLS::Getter<Kinematics_Base,Kernel_Info,Kinematics_FF2_Coll>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"Collinear 1->2 kinematics (FF)";
}


