#include "CFPSHOWER++/Calculators/FF/Kinematics_FF2_Soft.H"
#include "CFPSHOWER++/Calculators/FF/Kinematics_FF2_Coll.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "ATOOLS/Org/Message.H"

using namespace CFPSHOWER;
using namespace ATOOLS;

Kinematics_FF2_Soft::Kinematics_FF2_Soft() :
  Kinematics_Base(log_type::soft, kernel_type::FF)
{
  m_name = "FF2(soft)";
}

Kinematics_FF2_Soft::~Kinematics_FF2_Soft() {}

bool Kinematics_FF2_Soft::operator()(Splitting & split,Configuration & config) {
  Init(split);
  if (!KinCheck(split,config)) return false;
  return ConstructSystem(split,config);
  /*
  Init(split);
  m_momsum = SumMomenta(config);
  if (!CalculateTransverseMomentum(split) ||
      !CalculateRescaleFactor() ||
      !KinCheck(split,config)) return false;
  ConstructSystem(split,config);
  CalculateWeight(split);
  */
}

double Kinematics_FF2_Soft::CalculateY(Splitting & split) {
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

bool Kinematics_FF2_Soft::KinCheck(Splitting & split, Configuration & config) {
  return (split.y()>=0.0 && split.y()<=1.0);
}

bool Kinematics_FF2_Soft::ConstructSystem(Splitting & split, Configuration & config) {
  split.Set_y(CalculateY(split));
  if (!KinCheck(split,config)) return false;
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
  /*
  m_moms.resize(config->size());
  RescaleAndBoostMomenta(split,config);
  split.SetMomentum(0,m_psplit);
  split.SetMomentum(1,m_pnew);
  split.SetSpectatorMomentum(m_pspect);
  */
}

void Kinematics_FF2_Soft::CalculateWeight(Splitting & split) {
  // this is only relevant for massive splittings
  //double Q2red = split.Q2red(), Q2red_y = Q2red * split.y(); 
  //m_weight = ( Q2red / Lambda(split.Q2(),split.msplit2(),split.mspect2()) *
  //	       Q2red_y /(Q2red_y+split.m2(0)+split.m2(1)-split.msplit2()) );
  m_weight = 1.-split.y();
}

bool Kinematics_FF2_Soft::UpdateSystem(Splitting & split,Configuration & config) {
  split.GetSpectator()->SetMom(split.SpectatorMomentum());
  return true;
}

bool Kinematics_FF2_Soft::CalculateTransverseMomentum(Splitting & split) {
  Vec4D nT = Vec4D(0., cross(Vec3D(m_psplit), Vec3D(m_pspect)));
  if (nT.PSpat2()<=1.e-8) {
    nT = Vec4D(0., 1., 1., 0.);
    Poincare zrot(m_psplit,Vec4D::ZVEC);
    zrot.RotateBack(nT);
  }
  Vec4D lT = LT(m_psplit, m_pspect, nT);
  m_pnew  = ((exp(split.y())*m_psplit + exp(-split.y())*m_pspect) / m_Q +
	     (nT/nT.PSpat() * cos(split.phi(0)) +
	      lT/sqrt(dabs(lT.Abs2())) * sin(split.phi(0))) ) * sqrt(split.t(0));
  return true;
}

bool Kinematics_FF2_Soft::CalculateRescaleFactor() {
  double A = (m_pboth*(m_momsum+m_pnew))/m_Q2;
  double B = (2.*m_momsum[0]*m_pnew[0])/m_Q2;
  if (A*A<B) return false;
  m_momscale = 1. - A + sqrt(A*A-B);
  return (m_momscale>=0.);
}

void Kinematics_FF2_Soft::RescaleAndBoostMomenta(Splitting & split, Configuration & config) {
  m_pplus     = m_momsum - (1.-m_momscale)*m_pboth;
  m_pminus    = m_pplus + m_pnew;
  m_newsystem = Poincare(m_pminus);
  m_psplit   *= m_momscale;
  m_pspect   *= m_momscale;
  m_pnew     *= m_momscale;
  m_newsystem.Boost(m_psplit);
  m_newsystem.Boost(m_pspect);
  m_newsystem.Boost(m_pnew);  
  for (Parton_List::iterator pit=config.begin();pit!=config.end();pit++) {
    if ((*pit)==split.GetSplitter())
      m_moms.push_back(m_psplit);
    else if ((*pit)==split.GetSpectator())
      m_moms.push_back(m_psplit);
    else {
      Vec4D ppit = m_momscale * (*pit)->Mom();
      m_newsystem.Boost(ppit);
      m_moms.push_back(ppit);
    }
  }
}

DECLARE_GETTER(Kinematics_FF2_Soft,"Kinematics_FF2_Soft",Kinematics_Base,Kernel_Info);

Kinematics_Base * ATOOLS::Getter<Kinematics_Base,Kernel_Info,Kinematics_FF2_Soft>::
operator()(const Parameter_Type & info) const
{
  if (info.Type()==kernel_type::FF &&
      info.LogType()==log_type::soft &&
      info.GetFlavs().size()==2) {
    return new Kinematics_FF2_Soft();
  }
  return NULL;
}

void ATOOLS::Getter<Kinematics_Base,Kernel_Info,Kinematics_FF2_Soft>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"Soft 1->2 kinematics (FF)";
}

