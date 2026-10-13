#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/Flavour.H"
#include "YFS/Main/Define_Dipoles.H"
#include "YFS/Main/Coulomb.H"
#include <complex>


using namespace MODEL;
using namespace YFS;

Coulomb::Coulomb(bool enabled) : m_coulomb(enabled) {
  m_s = sqr(rpa->gen.Ecms());
  m_MW = Flavour(kf_Wplus).Mass();
  m_GW = Flavour(kf_Wplus).Width();
  m_cms = rpa->gen.Ecms();
  m_subtract = 0;
  if(m_coulomb) rpa->gen.AddCitation(1,"Coulomb corrections for WW threshold as described in \\cite{Bardin:1993mc,Fadin:1995fp,Jadach:1995sp}");
}


Coulomb::~Coulomb(){
  
}


void Coulomb::Calculate(const ATOOLS::Vec4D p1, const ATOOLS::Vec4D p2){
  m_s1 = p1.Abs2();
  m_s2 = p2.Abs2();
  //Eq 9 in https://arxiv.org/pdf/hep-ph/9507422.pdf
  double E = (m_s-4*sqr(m_MW))/(4*m_MW);
  double sarg = sqrt(sqr(E)+sqr(m_GW));
  double pp1 = sqrt(0.5*m_MW*(sarg-E));
  double pp2 = sqrt(0.5*m_MW*(sarg+E));
  double pp = 1./(4.*m_s)*(sqr(m_s) -2*m_s*(m_s1+m_s2) +sqr(m_s1-m_s2));
  double p  = sqrt(pp);
  m_p = p;
  double absKappa = m_MW*sqrt(sqr(E)+sqr(m_GW));
  double arg = (absKappa-pp)/(2*p*pp1);
  m_weight = 1 + m_alpha*sqrt(m_s)/(4*p)*(M_PI - 2.*atan(arg));
  if(IsBad(m_weight)){
    msg_Error()<<METHOD<<std::endl
                   <<"E = "<< E<<std::endl
                   <<"sqrt(s) = "<< sqrt(m_s)<<std::endl
                   <<"sqrt(s1) = "<< sqrt(m_s1)<<std::endl
                   <<"sqrt(s2) = "<< sqrt(m_s1-m_s2)<<std::endl
                   <<"p1 = "<< pp1<<std::endl
                   <<"p2 = "<< pp2<<std::endl
                   <<"pp = "<< pp<<std::endl
                   <<"k = "<< absKappa<<std::endl
                   <<"arg = "<< arg<<std::endl
                   <<"arctan(arg) = "<< atan(arg)<<std::endl
                   <<"weight = "<< m_weight<<std::endl;
  }
}


void Coulomb::Subtract(){
  m_weight-= sqrt(m_s)*m_alpha/(2*m_p)*M_PI;
  // m_weight-= 0;
}


namespace {
  typedef std::complex<double> cplx;
  //! Gauss-Legendre nodes/weights on [0,1] (n = 48), built once
  struct GL {
    std::vector<double> x, w;
    GL(int n) {
      for (int i(1); i <= n; ++i) {
        double z(cos(M_PI*(i - 0.25)/(n + 0.5))), pp(0.);
        for (int it(0); it < 100; ++it) {
          double p1(1.), p2(0.);
          for (int j(1); j <= n; ++j) { const double p3(p2); p2 = p1; p1 = ((2.*j - 1.)*z*p2 - (j - 1.)*p3)/j; }
          pp = n*(z*p1 - p2)/(z*z - 1.);
          const double dz(p1/pp); z -= dz;
          if (std::abs(dz) < 1e-15) break;
        }
        x.push_back(0.5*(1. - z)); w.push_back(1./((1. - z*z)*pp*pp));
      }
    }
  };
  const GL &Gauss() { static const GL g(48); return g; }

  //! FKMS kinematics: p, kappa = p1 - i p2 (eqs. 3-6)
  bool Kin(double s, double s1, double s2, double MW, double GW, double &p, cplx &kappa)
  {
    const double pp((sqr(s) - 2.*s*(s1 + s2) + sqr(s1 - s2))/(4.*s));
    if (!(pp > 0.) || !(s > 0.)) return false;
    p = sqrt(pp);
    const double E((s - 4.*sqr(MW))/(4.*MW)), r(sqrt(sqr(E) + sqr(GW)));
    kappa = cplx(sqrt(0.5*MW*(r - E)), -sqrt(0.5*MW*(r + E)));
    return true;
  }
}

double Coulomb::FirstOrder(double s, double s1, double s2) const
{
  double p; cplx kappa;
  if (!Kin(s, s1, s2, m_MW, m_GW, p, kappa)) return 0.;
  // eq. (9); Re(kappa) = p1
  return m_alpha*sqrt(s)/(4.*p)*(M_PI - 2.*atan((std::norm(kappa) - sqr(p))/(2.*p*kappa.real())));
}

double Coulomb::HigherOrderWeight(double s, double s1, double s2, int order) const
{
  double p; cplx kappa;
  if (!Kin(s, s1, s2, m_MW, m_GW, p, kappa)) return 1.;
  const cplx I(0., 1.);
  const double rs(sqrt(s)), d1(FirstOrder(s, s1, s2));
  const GL &g(Gauss());
  if (order == 1) {
    // eq. (7): f = 1 + f1 + f2, |f|^2 - 1 - 2 Re f1 = 2 Re f2 + |f1|^2
    const cplx D((kappa + I*p)/(kappa - I*p));
    const cplx f1(m_alpha*rs/(4.*I*p)*std::log(D));
    cplx J(0.);
    for (size_t i(0); i < g.x.size(); ++i) {
      const double x(g.x[i]);
      J += g.w[i]*std::log((1. + x*D)/(1. + x/D))/x;
    }
    const cplx f2(sqr(m_alpha)*s/(16.*I*p*kappa)*J);
    return 1. + 2.*f2.real() + std::norm(f1);
  }
  // eq. (2), all orders: f = 1 + alpha sqrt(s) kappa int_0^1 dx x^(-a) /
  // (kappa^2 (1+x)^2 + p^2 (1-x)^2), a = alpha sqrt(s)/(4 kappa); x = t^2
  // tames the integrable x^(-a) at x -> 0 (Re a > 0, |a| << 1)
  const cplx a(m_alpha*rs/(4.*kappa));
  cplx J(0.);
  for (size_t i(0); i < g.x.size(); ++i) {
    const double t(g.x[i]), x(t*t);
    if (x <= 0.) continue;
    J += g.w[i]*2.*t*std::exp(-a*log(x))/(kappa*kappa*sqr(1. + x) + sqr(p)*sqr(1. - x));
  }
  const cplx f(1. + m_alpha*rs*kappa*J);
  return std::norm(f) - d1;
}
