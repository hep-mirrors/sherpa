#include "CFPSHOWER++/Calculators/FF/SF_FF2_Coll.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "ATOOLS/Org/Message.H"

using namespace CFPSHOWER;
using namespace ATOOLS;

bool SF_FF2_Coll::Construct(Splitting & split,Configuration & config) {
  if (!KinCheck(split)) return false;
  return ConstructSystem(split);
}

bool SF_FF2_Coll::KinCheck(Splitting & split) {
  split.SetY(m_y = CalculateY(split));
  return (m_y>=0.0 && m_y<=1.0);
}

void SF_FF2_Coll::CalculateInvariants(Splitting & split) {
  m_Q2 = 0.;
  for (size_t i=0;i<2;i++) {
    m_pp[i][i] = 0.;
    for (size_t j=i;j<3;j++) m_Q2 += m_pp[i][j] = m_pp[j][i] = m_moms[i]*m_moms[j];
  }
  for (size_t i=0;i<2;i++) m_z[i] = m_pp[i][2]/(m_pp[0][2]+m_pp[1][2]);
  m_y = m_pp[0][1]/m_Q2;
}

double SF_FF2_Coll::CalculateY(Splitting & split) {
  if (split.GetSplitter()->Flav().IsFermion())
    return split.T()/((1.-split.Z())*split.Q2());
  if (split.GetSplitter()->Flav().IsVector()) {
    if (split.GetKernel()->GetFlavs()[0].IsVector())
      return split.T()/((1.-split.Z())*split.Z()*split.Q2());
    if (split.GetKernel()->GetFlavs()[0].IsFermion())
      return split.T()/split.Q2();
  }
  return split.T()/((1.-split.Z())*split.Z()*split.Q2());
}

bool SF_FF2_Coll::ConstructSystem(Splitting & split) {
  if (!KinCheck(split)) return false;
  PHASIC::Kin_Args kinargs(m_y,split.Z(),split.Phi());
  if (PHASIC::ConstructFFDipole(m_m2[0], m_m2[1], m_msplit2,m_mspect2,
				m_psplit,m_pspect,kinargs) < 0) return false;
  m_moms[0] = kinargs.m_pi;
  m_moms[1] = kinargs.m_pj;
  m_moms[2] = kinargs.m_pk;
  return true;
}

void SF_FF2_Coll::CalculateJacobean(Splitting & split) {
  m_weight = 1.-m_y;
  if (false && m_ismassive) {
    double Q2red   = m_Q2-m_m2[0]-m_m2[1]-m_mspect2;
    double Q2red_y = Q2red * m_y;
    m_weight *= ( Q2red / Lambda(m_Q2,m_msplit2,m_mspect2) *
		  Q2red_y /(Q2red_y+m_m2[0]+m_m2[1]-m_msplit2) );
  }
}

bool SF_FF2_Coll::UpdateSystem(Splitting & split,Configuration & config) {
  split.GetSpectator()->SetMom(m_moms[2]);
  return true;
}

