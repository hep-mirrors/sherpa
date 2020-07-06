#include "CFPSHOWER++/Calculators/FF/SF_FF2_Soft.H"
#include "CFPSHOWER++/Calculators/FF/SF_FF2_Coll.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "ATOOLS/Org/Message.H"

using namespace CFPSHOWER;
using namespace ATOOLS;

bool SF_FF2_Soft::Construct(Splitting & split,Configuration & config) {
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

void SF_FF2_Soft::CalculateInvariants(Splitting & split) {
  m_Q2 = 0.;
  for (size_t i=0;i<2;i++) {
    m_pp[i][i] = i<2 ? m_m2[i] : m_mspect2;
    for (size_t j=i+1;j<3;j++) {
      m_Q2 += m_pp[i][j] = m_pp[j][i] = m_moms[i]*m_moms[j];
    }
  }
  for (size_t i=0;i<2;i++) m_z[i] = m_pp[i][2]/(m_pp[0][2]+m_pp[1][2]);
  m_y = m_pp[0][1]/m_Q2;
}

double SF_FF2_Soft::CalculateY(Splitting & split) {
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

bool SF_FF2_Soft::KinCheck(Splitting & split, Configuration & config) {
  split.SetY(m_y = CalculateY(split));
  return (m_y>=0.0 && m_y<=1.0);
}

bool SF_FF2_Soft::ConstructSystem(Splitting & split, Configuration & config) {
  if (!KinCheck(split,config)) return false;
  PHASIC::Kin_Args kinargs(m_y,split.Z(),split.Phi());
  if (PHASIC::ConstructFFDipole(m_m2[0], m_m2[1], m_msplit2,m_mspect2,
				m_psplit,m_pspect,kinargs) < 0) return false;
  m_moms[0] = kinargs.m_pi;
  m_moms[1] = kinargs.m_pj;
  m_moms[2] = kinargs.m_pk;
  return true;
}

void SF_FF2_Soft::CalculateJacobean(Splitting & split) {
  m_weight = 1.-m_y;
  if (false && m_ismassive) {
    double Q2red   = m_Q2-m_m2[0]-m_m2[1]-m_mspect2;
    double Q2red_y = Q2red * m_y;
    m_weight *= ( Q2red / Lambda(m_Q2,m_msplit2,m_mspect2) *
		  Q2red_y /(Q2red_y+m_m2[0]+m_m2[1]-m_msplit2) );
  }
}

bool SF_FF2_Soft::UpdateSystem(Splitting & split,Configuration & config) {
  split.GetSpectator()->SetMom(m_moms[2]);
  return true;
}
/*
bool SF_FF2_Soft::CalculateTransverseMomentum(Splitting & split) {
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

bool SF_FF2_Soft::CalculateRescaleFactor() {
  double A = (m_pboth*(m_momsum+m_pnew))/m_Q2;
  double B = (2.*m_momsum[0]*m_pnew[0])/m_Q2;
  if (A*A<B) return false;
  m_momscale = 1. - A + sqrt(A*A-B);
  return (m_momscale>=0.);
}

void SF_FF2_Soft::RescaleAndBoostMomenta(Splitting & split, Configuration & config) {
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
*/
