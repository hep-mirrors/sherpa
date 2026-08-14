#include "ATOOLS/Math/BreitBoost.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace ATOOLS;

BreitBoost::BreitBoost(const Vec4D& Q,
                       const Vec4D& hadin) {
  _init(Q,hadin);
}

BreitBoost::BreitBoost(const Vec4D& lepin, const Vec4D& lepout,
                       const Vec4D& hadin) {
  _init(lepin-lepout,hadin);
}

BreitBoost::BreitBoost(Cluster_Amplitude *const ampl) {
  // identify lepton and hadron incoming legs
  int lep_leg(-1), had_leg(-1);
  for (size_t i(0); i < ampl->NIn(); ++i) {
    if (ampl->Leg(i)->Flav().IsLepton()) lep_leg = i;
    else                                 had_leg = i;
  }
  if (lep_leg < 0 || had_leg < 0)
    THROW(fatal_error, "Cannot identify lepton/hadron legs in BreitBoost");
  const Vec4D lepin = -ampl->Leg(lep_leg)->Mom();
  const Vec4D hadin = rpa->gen.PBeam(had_leg);

  Vec4D lepout(0,0,0,0);
  for (size_t i(ampl->NIn()); i < ampl->Legs().size(); ++i) {
    if (ampl->Leg(i)->Flav().IsLepton()) {
      lepout += ampl->Leg(i)->Mom();
    }
  }
  /// construct from momenta, with
  /// incoming hadron projected to zero mass
  _init(lepin-lepout, {hadin.PSpat(), Vec3D(hadin)});
}



void BreitBoost::Apply(Cluster_Amplitude *const ampl) const {
  for (size_t i(0);i<ampl->Legs().size();++i) {
    ampl->Leg(i)->SetMom((*this)*ampl->Leg(i)->Mom());
  }
}

void BreitBoost::_init(const Vec4D& Q,
                       const Vec4D& hadin) {
  m_Q2 = Q.Abs2();
  m_x = Min(-m_Q2/(2.*hadin*Q),1.);
  const double E = sqrt(-m_Q2)/(2.*m_x);

  const Vec4D p(sqrt(E*E+hadin.Abs2()),0.,0.,-E);
  const Vec4D q(0.,0.,0.,2.*m_x*E);
  reserve(3);
  emplace_back(hadin+Q);                  // boost to cms
  emplace_back(back()*hadin,-Vec4D::ZVEC); // align with z axis
  emplace_back(p+q);                      // boost from breit
  back().Invert();                        // needs to be inverted
}
