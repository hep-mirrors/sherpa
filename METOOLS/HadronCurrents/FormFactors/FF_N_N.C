#include "METOOLS/HadronCurrents/FormFactors/FF_N_N.H"
#include "METOOLS/HadronCurrents/Tools.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "MODEL/Main/Model_Base.H"

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

FF_N_N::FF_N_N(const FF_Parameters & params) :
  FormFactor_Base(params), m_vertex(FF_N_N_vertex::unknown)
{
  if      (m_name=="F1") m_fftype = FFtype::F1;
  else if (m_name=="F2") m_fftype = FFtype::F2;
  else if (m_name=="FA") m_fftype = FFtype::FA;
  else if (m_name=="FP") m_fftype = FFtype::FP;
  else THROW(fatal_error,"FF_N_N has no form factor '"+m_name+"'; use "
             "F1, F2, FA or FP.");
  ValidateModel();

  GeneralModel empty;
  const GeneralModel & md = p_model ? *p_model : empty;
  FixMode(md);
  m_mN  = 0.5*(m_masses[m_pi[0]]+m_masses[m_pi[1]]);
  m_mpi = Flavour(kf_pi).HadMass();
  m_s2  = std::abs(MODEL::s_model->ComplexConstant(string("csin2_thetaW")));
  m_T3  = (m_vertex==FF_N_N_vertex::Z_p) ? 0.5 : -0.5;
  m_gA  = md("GA",1.2754);

  // p_model dangles once the owning current's setup returns, as in FF_P_X.
  p_model = NULL;
  msg_Tracking()<<"### FF_N_N '"<<m_name<<"': model "<<m_ffmodel
                <<", F(0) = "<<Value(0.)<<", F(1 GeV^2) = "<<Value(1.)<<"\n";
}

void FF_N_N::FixMode(const GeneralModel & md)
{
  const kf_code kf0(m_flavs[m_pi[0]].Kfcode()), kf1(m_flavs[m_pi[1]].Kfcode());
  const bool cc(kf0!=kf1), isp(kf0==kf_p_plus);
  const int  boson(std::abs(int(md("BOSON", cc ? 24. : 22.))));
  if      (boson==22 && !cc) m_vertex = isp ? FF_N_N_vertex::gamma_p : FF_N_N_vertex::gamma_n;
  else if (boson==23 && !cc) m_vertex = isp ? FF_N_N_vertex::Z_p : FF_N_N_vertex::Z_n;
  else if (boson==24 &&  cc) m_vertex = FF_N_N_vertex::W_pn;
  else THROW(fatal_error,"FF_N_N: BOSON "+std::to_string(boson)+
             " does not match the nucleon legs.");
}

////////////////////////////////////////////////////////////////////////////
// Vector models. G_E = F1 - tau F2, G_M = F1 + F2.
////////////////////////////////////////////////////////////////////////////
void FF_N_N::Kelly_p(const double & Q2,double & F1,double & F2) const
{
  // Kelly, PRC 70 (2004) 068202, doi:10.1103/PhysRevC.70.068202
  const double mu_p(2.793);
  const double tau(Q2/(4.*sqr(m_mN)));
  const double GE((1.-0.24*tau)/(1.+tau*(10.98+tau*(12.82+tau*21.97))));
  const double GM(mu_p*(1.+0.12*tau)/(1.+tau*(10.97+tau*(18.86+tau*6.55))));
  F1 = (GE+tau*GM)/(1.+tau);
  F2 = (GM-GE)/(1.+tau);
}

void FF_N_N::Kelly_n(const double & Q2,double & F1,double & F2) const
{
  // Kelly, PRC 70 (2004) 068202; G_E^n Galster-like A tau/(1+B tau) G_D
  const double mu_n(-1.913), A(1.7), B(3.3), L2(0.7174);
  const double tau(Q2/(4.*sqr(m_mN)));
  const double GE(A*tau/(1.+B*tau)/sqr(1.+Q2/L2));
  const double GM(mu_n*(1.+2.33*tau)/(1.+tau*(14.72+tau*(24.2+tau*84.1))));
  F1 = (GE+tau*GM)/(1.+tau);
  F2 = (GM-GE)/(1.+tau);
}

////////////////////////////////////////////////////////////////////////////
// Axial models. F_A = g_A Ftilde(Q^2)/Ftilde(0), positive-F_A convention.
////////////////////////////////////////////////////////////////////////////
double FF_N_N::Meyer(const double & Q2) const
{
  // Meyer, Betancourt, Gran, Hill, PRD 93 (2016) 113015, arXiv:1603.03048
  const double tc(9.*sqr(Flavour(kf_pi_plus).HadMass())), t0(-0.28);
  const std::vector<double> ak{0.759,-2.30,0.6,3.8,-2.3,-2.16,0.896,1.58,-0.823};
  auto F = [&](const double & Qsq) {
    const double A(sqrt(Max(0.,tc+Qsq))), B(sqrt(tc-t0)), z((A-B)/(A+B));
    double zk(1.), sum(0.);
    for (const double & a : ak) { sum += a*zk; zk *= z; }
    return sum;
  };
  return m_gA*F(Q2)/F(0.);
}

double FF_N_N::Lattice(const double & Q2) const
{
  // lattice-QCD average z-expansion
  const double tc(9.*sqr(0.134)), t0(-0.50);
  const std::vector<double> ak{0.72115656,-1.70104640,0.26324902,1.53433681,
                               0.01061114,-1.49893610,0.67062898};
  auto F = [&](const double & Qsq) {
    const double A(sqrt(Max(0.,tc+Qsq))), B(sqrt(tc-t0)), z((A-B)/(A+B));
    double zk(1.), sum(0.);
    for (const double & a : ak) { sum += a*zk; zk *= z; }
    return sum;
  };
  return m_gA*F(Q2)/F(0.);
}

////////////////////////////////////////////////////////////////////////////
// Model switches: the only place a new model number is wired in.
////////////////////////////////////////////////////////////////////////////
void FF_N_N::Vector(const double & Q2,double & F1p,double & F2p,
                    double & F1n,double & F2n) const
{
  F1p = F2p = F1n = F2n = 0.;
  switch (GetModel()) {
    case Model::Point: F1p = 1.; return;
    case Model::Kelly: Kelly_p(Q2,F1p,F2p); Kelly_n(Q2,F1n,F2n); return;
    default: return;
  }
}

double FF_N_N::Axial(const double & Q2) const
{
  switch (GetModel()) {
    case Model::Point:   return m_gA;
    case Model::Meyer:   return Meyer(Q2);
    case Model::Lattice: return Lattice(Q2);
    default:             return 0.;
  }
}

////////////////////////////////////////////////////////////////////////////
// Boson combinations.
////////////////////////////////////////////////////////////////////////////
double FF_N_N::Value(const double & Q2) const
{
  if (m_fftype==FFtype::FA || m_fftype==FFtype::FP) {
    double FA(Axial(Q2));
    if      (m_vertex==FF_N_N_vertex::gamma_p || m_vertex==FF_N_N_vertex::gamma_n) FA = 0.;
    else if (m_vertex==FF_N_N_vertex::Z_p     || m_vertex==FF_N_N_vertex::Z_n)     FA *= m_T3;
    return m_fftype==FFtype::FA ? FA : 4.*sqr(m_mN)/(Q2+sqr(m_mpi))*FA;
  }
  double F1p, F2p, F1n, F2n;
  Vector(Q2,F1p,F2p,F1n,F2n);
  const double Fp(m_fftype==FFtype::F1 ? F1p : F2p), Fn(m_fftype==FFtype::F1 ? F1n : F2n);
  switch (m_vertex) {
    case FF_N_N_vertex::gamma_p: return Fp;
    case FF_N_N_vertex::gamma_n: return Fn;
    case FF_N_N_vertex::Z_p:     return m_T3*(Fp-Fn)-2.*m_s2*Fp;
    case FF_N_N_vertex::Z_n:     return m_T3*(Fp-Fn)-2.*m_s2*Fn;
    case FF_N_N_vertex::W_pn:    return Fp-Fn;
    default:                   return 0.;
  }
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
