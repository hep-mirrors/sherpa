#include "METOOLS/HadronCurrents/FormFactors/FF_N_N.H"
#include "METOOLS/HadronCurrents/Tools.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

FF_N_N::FF_N_N(const FF_Parameters & params) :
  FormFactor_Base(params),
  m_tc(0.), m_t0(0.), m_FA0(1.)
{
  if      (m_name=="F1p") m_which = Which::F1p;
  else if (m_name=="F2p") m_which = Which::F2p;
  else if (m_name=="F1n") m_which = Which::F1n;
  else if (m_name=="F2n") m_which = Which::F2n;
  else if (m_name=="FA")  m_which = Which::FA;
  else THROW(fatal_error,"FF_N_N has no form factor '"+m_name+"'; use "
             "F1p, F2p, F1n, F2n or FA.");

  const Model mod(ModelAs<Model>());
  const bool axial(m_which==Which::FA);
  if (!(mod==Model::None || mod==Model::Point ||
        (!axial && mod==Model::Kelly) ||
        (axial && (mod==Model::Meyer || mod==Model::Lattice))))
    THROW(fatal_error,"FF_N_N '"+m_name+"': unknown model "
          +std::to_string(m_ffmodel)+". Vector: 0, 1, 10 (Kelly). "
          "Axial: 0, 1, 20 (Meyer deuterium), 21 (lattice).");

  GeneralModel empty;
  const GeneralModel & md = p_model ? *p_model : empty;
  m_mN  = 0.5*(m_masses[m_pi[0]]+m_masses[m_pi[1]]);
  m_muP = md("MU_P",2.793);
  m_muN = md("MU_N",-1.913);
  m_L2  = md("LAMBDA2",0.7174);
  m_An  = md("KELLY_AN",1.7);
  m_Bn  = md("KELLY_BN",3.3);
  m_gA  = md("GA",1.2754);

  // z-expansions, positive F_A convention of the tex table.
  if (mod==Model::Meyer) {
    m_tc = 9.*sqr(Flavour(kf_pi_plus).HadMass()); m_t0 = -0.28;
    m_ak = {0.759,-2.30,0.6,3.8,-2.3,-2.16,0.896,1.58,-0.823};
  }
  else if (mod==Model::Lattice) {
    m_tc = 9.*sqr(0.134); m_t0 = -0.50;
    m_ak = {0.72115656,-1.70104640,0.26324902,1.53433681,
            0.01061114,-1.49893610,0.67062898};
  }
  // F_A = g_A Ftilde(Q^2)/Ftilde(0): Ftilde(0)=1.2723 for Meyer.
  if (!m_ak.empty()) m_FA0 = ZSeries(0.);

  // p_model dangles once the owning current's setup returns, as in FF_P_X.
  p_model = NULL;
  msg_Tracking()<<"### FF_N_N '"<<m_name<<"': model "<<m_ffmodel
                <<", F(0) = "<<Value(0.)<<", F(1 GeV^2) = "<<Value(1.)<<"\n";
}

double FF_N_N::Kelly(const double & a1,const double & b1,const double & b2,
                     const double & b3,const double & tau) const
{
  // G = (1 + a1 tau)/(1 + b1 tau + b2 tau^2 + b3 tau^3), a0 = 1.
  return (1.+a1*tau)/(1.+tau*(b1+tau*(b2+tau*b3)));
}

double FF_N_N::Z(const double & Q2) const
{
  const double A(sqrt(Max(0.,m_tc+Q2))), B(sqrt(m_tc-m_t0));
  return (A-B)/(A+B);
}

double FF_N_N::ZSeries(const double & Q2) const
{
  const double z(Z(Q2));
  double zk(1.), sum(0.);
  for (size_t k(0);k<m_ak.size();++k) { sum += m_ak[k]*zk; zk *= z; }
  return sum;
}

double FF_N_N::Value(const double & Q2) const
{
  const Model mod(ModelAs<Model>());
  if (mod==Model::None) return 0.;

  if (m_which==Which::FA) {
    if (mod==Model::Point) return m_gA;
    return m_gA*ZSeries(Q2)/m_FA0;
  } 
  // need the above to allow us to set g_A to whatever we like
  //  while using the axial FF models we want to, which may use different norms.

  const bool p (m_which==Which::F1p || m_which==Which::F2p);
  const bool f1(m_which==Which::F1p || m_which==Which::F1n);
  if (mod==Model::Point) return (p && f1) ? 1. : 0.;

  // Kelly Sachs form factors. G_E^n in closed form, A tau/(1+B tau) G_D,
  // which is identical to the rational form of the tex table.
  const double tau(Q2/(4.*sqr(m_mN)));
  double GE, GM;
  if (p) {
    GE = Kelly(-0.24,10.98,12.82,21.97,tau);
    GM = m_muP*Kelly(0.12,10.97,18.86,6.55,tau);
  }
  else {
    GE = m_An*tau/(1.+m_Bn*tau)/sqr(1.+Q2/m_L2);
    GM = m_muN*Kelly(2.33,14.72,24.2,84.1,tau);
  }
  // G_E = F1 - tau F2, G_M = F1 + F2.
  return f1 ? (GE+tau*GM)/(1.+tau) : (GM-GE)/(1.+tau);
}

Complex FF_N_N::operator()(const ATOOLS::Vec4D_Vector& moms)
{
  const Vec4D q(moms[m_pi[0]]-moms[m_pi[1]]);
  return Complex(Value(-q.Abs2()),0.);
}

DECLARE_FF_GETTER(FF_N_N,"FF_N_N")

FormFactor_Base * ATOOLS::Getter<FormFactor_Base,FF_Parameters,FF_N_N>::
operator()(const METOOLS::FF_Parameters &params) const
{
  if (params.m_pi.size()!=2) return NULL;
  return new FF_N_N(params);
}
