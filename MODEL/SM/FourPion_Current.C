#include "MODEL/SM/FourPion_Current.H"
#include "METOOLS/Explicit/Form_Factor.H"
#include "METOOLS/Explicit/Current.H"
#include "METOOLS/Explicit/Vertex.H"
#include "MODEL/Main/Single_Vertex.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Phys/Flavour.H"

using namespace METOOLS;
using namespace ATOOLS;

// ----------------------------------------------------------------------
// parameters: the fit of arXiv:0804.0359; omega -> 3 pi of arXiv:hep-ph/0512180
// ----------------------------------------------------------------------

void FourPion_Parameters::ReadMW(Scoped_Settings s,double &m,double &g) {
  m = s["Mass"].SetDefault(m).Get<double>();
  g = s["Width"].SetDefault(g).Get<double>();
}

void FourPion_Parameters::ReadB(Scoped_Settings s,double b[3]) {
  const std::vector<double> v(s.SetDefault(std::vector<double>{b[0],b[1],b[2]})
                              .GetVector<double>());
  if (v.size()!=3) THROW(fatal_error,"FOUR_PION_FORM_FACTOR: Q^2 tower weights need 3 entries.");
  for (size_t i(0); i<3; ++i) b[i] = v[i];
}

void FourPion_Parameters::Read() {
  Scoped_Settings s{ Settings::GetMainSettings()["FOUR_PION_FORM_FACTOR"] };
  // fixed masses and widths, arXiv:0804.0359 appendix: a1 (A.18), rho, rho', rho''
  // (A.13, PDG), omega (A.25), f0(1370) (A.22)
  a1m = 1.23; a1g = 0.2; rhom = 0.7755; rhog = 0.1494; rho1m = 1.459; rho1g = 0.4;
  rho2m = 1.72; rho2g = 0.25; omm = 0.78265; omg = 0.00849; f0m = 1.35; f0g = 0.2;
  pim = 0.13957018;
  // the fitted couplings and Q^2 towers of arXiv:0804.0359 Table II, as PHOKHARA 10
  // stores them (const_and_model_paramall10.0.dat, its own normalisation of the c's)
  c_a1 = 112.388625; c_f0 = 21.412294; c_rr = 2.46245925; c_om = 1.46987657;
  tm[0] = 1.437; tg[0] = 0.52; tm[1] = 1.73829367; tg[1] = 0.45; tm[2] = 2.12538719; tg[2] = 0.3;
  const double ba1[3] = { -0.0660483408, -0.0214852818, -0.0042749303 };
  const double bf0[3] = { 72691.5547, -2589.85498, 1956.40495 };
  const double bom[3] = { -0.331507263, 0.0131990469, -0.00526863683 };
  for (size_t i(0); i<3; ++i) { b_a1[i] = ba1[i]; b_f0[i] = bf0[i]; b_om[i] = bom[i]; }
  // arXiv:0804.0359: g_omegarhopi, g_rhopipi (A.7), g_rhogamma (below A.10),
  // beta of B_rho (A.15)
  g_omrhopi = 42.3; g_rhopipi = 5.997; g10 = 0.1212; beta_t = -0.145;
  ReadMW(s["a1"],a1m,a1g);
  ReadMW(s["rho"],rhom,rhog);
  ReadMW(s["rho'"],rho1m,rho1g);
  ReadMW(s["rho''"],rho2m,rho2g);
  ReadMW(s["omega"],omm,omg);
  ReadMW(s["f0"],f0m,f0g);
  pim = s["m_pi"].SetDefault(pim).Get<double>();
  c_a1 = s["Coupling_a1"].SetDefault(c_a1).Get<double>();
  c_f0 = s["Coupling_f0"].SetDefault(c_f0).Get<double>();
  c_rr = s["Coupling_rhorho"].SetDefault(c_rr).Get<double>();
  c_om = s["Coupling_omega"].SetDefault(c_om).Get<double>();
  for (size_t i(0); i<3; ++i) ReadMW(s["Tower"][ToString(i+1)],tm[i],tg[i]);
  ReadB(s["Tower_a1"],b_a1);
  ReadB(s["Tower_f0"],b_f0);
  ReadB(s["Tower_omega"],b_om);
  g_omrhopi = s["g_omegarhopi"].SetDefault(g_omrhopi).Get<double>();
  g_rhopipi = s["g_rhopipi"].SetDefault(g_rhopipi).Get<double>();
  g10 = s["g_rhorho"].SetDefault(g10).Get<double>();
  beta_t = s["beta_rho"].SetDefault(beta_t).Get<double>();
}

// ----------------------------------------------------------------------
// the current
// ----------------------------------------------------------------------

// g(Q^2) of arXiv:0804.0359 eq. (A.17), from Kuhn, Santamaria, Z. Phys. C48 (1990)
// 445. The paper switches at (m_a1 + m_pi)^2, PHOKHARA 10 (gfun_4pi) and this
// code at (m_rho + m_pi)^2.
double FourPionCurrent::GA1(double s,const FourPion_Parameters &p) {
  if (s>sqr(p.rhom+p.pim)) return s*1.623+10.38-9.32/s+0.65/sqr(s);
  const double c1(s-9.*sqr(p.pim));
  return 4.1*pow(c1,3)*(1.-3.3*c1+5.8*sqr(c1));
}

const FourPionCurrent &FourPionCurrent::Instance() {
  static const FourPionCurrent cur;
  return cur;
}

FourPionCurrent::FourPionCurrent() {
  m_par.Read();
  m_par.con = m_par.a1g*m_par.a1m/GA1(sqr(m_par.a1m),m_par);
}

// BW3 of arXiv:0804.0359 eq. (A.12)
Complex FourPionCurrent::BWrho(double m,double g,double s) const {
  const double sth(4.*sqr(m_par.pim));
  const double gs(s>0. ? g*sqr(m)/s*pow(std::abs(s-sth)/(sqr(m)-sth),1.5) : 0.);
  return sqr(m)/Complex(sqr(m)-s,-std::sqrt(std::max(s,0.))*gs);
}

// F_rho(Q^2, beta) of arXiv:0804.0359 eq. (A.11)
Complex FourPionCurrent::Tower(const double b[3],double q2) const {
  const FourPion_Parameters &p(m_par);
  const double c2(4.*sqr(p.pim)/q2);
  auto bw = [&](double m,double g) {
    const double c1(sqr(m)/q2);
    const double mg(m*g*std::sqrt(std::abs(c1*pow((1.-c2)/(c1-c2),3))));
    return Complex(sqr(m),0.)/Complex(sqr(m)-q2,-mg);
  };
  Complex res(bw(p.rhom,p.rhog));
  for (size_t k(0); k<3; ++k) res += b[k]*bw(p.tm[k],p.tg[k]);
  return res/(1.+b[0]+b[1]+b[2]);
}

Complex FourPionCurrent::OmegaPi(double q2) const {
  return m_par.c_om*Tower(m_par.b_om,q2);
}

Complex FourPionCurrent::OmegaThreePion(double spm,double sp0,double sm0) const {
  const FourPion_Parameters &p(m_par);
  return p.g_omrhopi*p.g_rhopipi*(BWrho(p.rhom,p.rhog,spm)+BWrho(p.rhom,p.rhog,sp0)
                                  +BWrho(p.rhom,p.rhog,sm0));
}

std::array<Complex,4> FourPionCurrent::Coefficients00(double q2,const Vec4D p[4],
                                                      bool omega) const {
  const FourPion_Parameters &P(m_par);
  const Vec4D &q1(p[0]), &q2v(p[1]), &q3(p[2]), &q4(p[3]);
  const Vec4D qq(q1+q2v+q3+q4);
  std::array<Complex,4> G = {{ 0., 0., 0., 0. }};
  const double q_q1(qq*q1), q_q2(qq*q2v), q_q3(qq*q3), q_q4(qq*q4);
  const double q1_2m4(q1*(q2v-q4)), q1_3m2(q1*(q3-q2v)), q3_2m4(q3*(q2v-q4));
  const double q2_3m1(q2v*(q3-q1)), q2_4m1(q2v*(q4-q1)), q3_4m1(q3*(q4-q1));
  const double q4_3m1(q4*(q3-q1)), q4_3m2(q4*(q3-q2v));
  const double s12((q1+q2v).Abs2()), s13((q1+q3).Abs2()), s14((q1+q4).Abs2());
  const double s23((q2v+q3).Abs2()), s24((q2v+q4).Abs2()), s34((q3+q4).Abs2());
  // the a1 pi and f0 rho terms (arXiv:0804.0359 appendix, as PHOKHARA 10 had3_4pi)
  {
    const double qmq3_2((q1+q2v+q4).Abs2()), qmq4_2((q1+q2v+q3).Abs2());
    // BW_a1, eq. (A.16)
    auto bwa1 = [&](double s) {
      return Complex(sqr(P.a1m),0.)/Complex(sqr(P.a1m)-s,-GA1(s,P)*P.con); };
    // B_rho, eq. (A.14)
    auto bwt = [&](double s) {
      return (BWrho(P.rhom,P.rhog,s)+P.beta_t*BWrho(P.rho1m,P.rho1g,s))/(1.+P.beta_t); };
    const Complex c0(Tower(P.b_a1,q2)), c5(bwa1(qmq3_2)), c6(bwa1(qmq4_2));
    const Complex t124(c5*bwt(s24)*P.c_a1*c0), t214(c5*bwt(s14)*P.c_a1*c0);
    const Complex t231(c6*bwt(s13)*P.c_a1*c0), t123(c6*bwt(s23)*P.c_a1*c0);
    // BW_f0, eq. (A.21), and T_rho with beta-bar = 0.08, -0.0075, eqs. (A.19, A.20)
    const double c2(4.*sqr(P.pim)/s12), c1(sqr(P.f0m)/s12);
    const double gf0(P.f0m*P.f0g*std::sqrt(std::abs(c1*(1.-c2)/(c1-c2))));
    const Complex bwf0(Complex(sqr(P.f0m),0.)/Complex(sqr(P.f0m)-s12,-gf0));
    const Complex bwr((BWrho(P.rhom,P.rhog,s34)+0.08*BWrho(P.rho1m,P.rho1g,s34)
                       -0.0075*BWrho(P.rho2m,P.rho2g,s34))/(1.+0.08-0.0075));
    const Complex ss(bwr*bwf0*P.c_f0*Tower(P.b_f0,q2));
    const Complex cf1(t123*(-1.-q1_3m2/qmq4_2)+t124*(1.-q1_2m4/qmq3_2)
                      +t214*(3.+q2_4m1/qmq3_2)+t231*(-3.-q2_3m1/qmq4_2));
    const Complex cf2(t123*(-3.-q1_3m2/qmq4_2)+t124*(3.-q1_2m4/qmq3_2)
                      +t214*(1.+q2_4m1/qmq3_2)+t231*(-1.-q2_3m1/qmq4_2));
    const Complex cf3(t123*(1.-q1_3m2/qmq4_2)+t124*(1.+q1_2m4/qmq3_2)
                      +t214*(1.-q2_4m1/qmq3_2)+t231*(1.-q2_3m1/qmq4_2)-3.*ss);
    G[0] += cf1; G[1] += cf2; G[2] += cf3; G[3] -= cf3;
  }
  if (omega) {
    // the omega pi0 term, arXiv:0804.0359 eqs. (A.5, A.6), BW_omega (A.24), H_rho
    // (A.23): 2 eps(mu, Q, q1, eps(., q2, q3, q4)) written out in dot products
    const double q12(q1*q2v), q13(q1*q3), q14(q1*q4);
    const double q23(q2v*q3), q24(q2v*q4);
    const Complex c(OmegaPi(q2));
    auto bwom = [&](double s) { return Complex(sqr(P.omm),0.)/Complex(sqr(P.omm)-s,-P.omm*P.omg); };
    const Complex A(c*bwom((q2v+q3+q4).Abs2())*OmegaThreePion(s34,s23,s24));
    const Complex B(c*bwom((q1+q3+q4).Abs2())*OmegaThreePion(s34,s13,s14));
    G[1] += A*2.*(q14*q_q3-q13*q_q4); G[2] += A*2.*(q12*q_q4-q14*q_q2); G[3] += A*2.*(q13*q_q2-q12*q_q3);
    G[0] += B*2.*(q24*q_q3-q23*q_q4); G[2] += B*2.*(q12*q_q4-q24*q_q1); G[3] += B*2.*(q23*q_q1-q12*q_q3);
  }
  // the "rho rho" term, arXiv:0804.0359 eqs. (A.8)-(A.10)
  {
    auto bwr = [&](double s) {
      return BWrho(P.rhom,P.rhog,s)/sqr(P.rhom)-BWrho(P.rho1m,P.rho1g,s)/sqr(P.rho1m); };
    const Complex c0r(P.g10*pow(P.g_rhopipi,3)*bwr(q2)*P.c_rr);
    const Complex c13(bwr(s13)), c14(bwr(s14)), c23(bwr(s23)), c24(bwr(s24));
    const Complex c11(c13*c24), c22(c14*c23);
    G[0] += c0r*(c11*(q_q2-q_q4+2.*q3_2m4)+c22*(-q_q2+q_q3+2.*q4_3m2)-2.*c14+2.*c13);
    G[1] += c0r*(c11*(-q_q1+q_q3+2.*q4_3m1)+c22*(q_q1-q_q4-2.*q3_4m1)-2.*c24+2.*c23);
    G[2] += c0r*(c11*(-q_q2+q_q4-2.*q1_2m4)+c22*(-q_q1+q_q4+2.*q2_4m1)-2.*c13-2.*c23);
    G[3] += c0r*(c11*(q_q1-q_q3-2.*q2_3m1)+c22*(q_q2-q_q3-2.*q1_3m2)+2.*c14+2.*c24);
  }
  return G;
}

std::array<Complex,4> FourPionCurrent::CoefficientsCC(double q2,const Vec4D p[4]) const {
  // the four J_00 terms: argument slots (0, 1 | 2, 3) filled by the indices below
  static const int perm[4][4] = { {0,1,2,3}, {3,1,2,0}, {0,2,1,3}, {3,2,1,0} };
  std::array<Complex,4> G = {{ 0., 0., 0., 0. }};
  for (const auto &pm : perm) {
    const Vec4D q[4] = { p[pm[0]], p[pm[1]], p[pm[2]], p[pm[3]] };
    const std::array<Complex,4> g(Coefficients00(q2,q,false));
    for (size_t i(0); i<4; ++i) G[pm[i]] += g[i];
  }
  return G;
}

// ----------------------------------------------------------------------
// FFFourPion: the coefficients at the vertex's own legs
// ----------------------------------------------------------------------
// FOUR_PION_FORM_FACTOR: the coefficients G_i of FourPionCurrent on the
// gamma* 4 pi contact vertex (VPPPP), at the photon virtuality and the
// momenta of the vertex's own pion legs. The legs are those declared in
// LowEnergy_Model::InitFourPionVertices: 0 the photon, then {pi+, pi-, pi0,
// pi0} or {pi+, pi+, pi-, pi-}; momenta outgoing from the vertex. At Born
// level and with initial-state radiation these are the external pions.
namespace METOOLS {
  class FFFourPion: public Linear_Current_FF {
  private:
    bool m_charged;
    std::vector<Complex> m_c;
  public:
    FFFourPion(const Vertex_Key &key);
    const std::vector<Complex> &Coefficients();
  };
}

FFFourPion::FFFourPion(const Vertex_Key &key):
  Linear_Current_FF("FourPion",key), m_charged(true), m_c(5,Complex(0.,0.))
{
  size_t n0(0);
  for (const Flavour &f : key.p_mv->in) if (f.Kfcode()==kf_pi) ++n0;
  m_charged = (n0==0);
  if (key.p_mv->in.size()!=5)
    THROW(fatal_error,"FourPion form factor on a vertex without five legs.");
}

const std::vector<Complex> &FFFourPion::Coefficients() {
  // momenta by declared leg, outgoing from the vertex
  const std::vector<size_t> &id(p_v->V()->id);
  Vec4D leg[5];
  for (size_t i(0); i+1<id.size(); ++i) leg[id[i]] = p_v->J(i)->P();
  leg[id.back()] = -p_v->JC()->P();
  const double q2(leg[0].Abs2());
  std::array<Complex,4> G;
  if (m_charged) {
    // declared {pi+, pi+, pi-, pi-} = legs 1, 2, 3, 4; J_cc(q1..q4): q2, q3
    // one charge, q1, q4 the other
    const size_t order[4] = { 3, 1, 2, 4 };
    const Vec4D p[4] = { leg[order[0]], leg[order[1]], leg[order[2]], leg[order[3]] };
    G = FourPionCurrent::Instance().CoefficientsCC(q2,p);
    for (size_t i(0); i<4; ++i) m_c[order[i]] = G[i];
  }
  else {
    // declared {pi+, pi-, pi0, pi0} = legs 1, 2, 3, 4; J_00(pi0, pi0, pi+, pi-)
    const size_t order[4] = { 3, 4, 1, 2 };
    const Vec4D p[4] = { leg[order[0]], leg[order[1]], leg[order[2]], leg[order[3]] };
    G = FourPionCurrent::Instance().Coefficients00(q2,p,true);
    for (size_t i(0); i<4; ++i) m_c[order[i]] = G[i];
  }
  return m_c;
}

DECLARE_GETTER(FFFourPion,"FourPion",Form_Factor,Vertex_Key);
Form_Factor *Getter<Form_Factor,Vertex_Key,FFFourPion>::
operator()(const Vertex_Key &args) const { return new FFFourPion(args); }
void ATOOLS::Getter<Form_Factor,Vertex_Key,FFFourPion>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"gamma* -> 4 pi, coefficients at the vertex's pion legs"; }
