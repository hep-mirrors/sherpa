#include "AHADIC++/Tools/Strong_Coupling.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "ATOOLS/Math/Random.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Message.H"

using namespace AHADIC;
using namespace ATOOLS;

Strong_Coupling::Strong_Coupling(const asform::code form,
				 const double pt02) :
  m_form(form), m_pt02(pt02),
  m_beta0(12.*M_PI/(33.-2.*4.)), m_gamma(m_beta0/M_PI),
  m_Lambda(0.349), m_Lambda2(sqr(m_Lambda)), 
  m_a(3.008), m_b(1.425), m_c(0.908), m_d(0.84), m_m2(sqr(1.204))
{
  MODEL::Running_AlphaS * as = static_cast<MODEL::Running_AlphaS *>
    (MODEL::s_model->GetScalarFunction(std::string("alpha_S")));
  switch (m_form) {
  case asform::GDH_inspired:
    //std::cout<<"Check this: gamma = "<<m_gamma<<" from beta = "<<m_beta0<<std::endl
    //	     <<"with further params {a,b,c,d} = {"<<m_a<<","<<m_b<<","<<m_c<<","<<m_d<<"}"
    //	     <<" and m, mg, n = "<<sqrt(m_m2)<<", "<<sqrt(mg2(0.))<<", "<<n(0.)<<std::endl;
    m_asmax   = (*this)(0.);
    break;
  case asform::constant:
  case asform::fall_off:
  case asform::fixed:
  case asform::fall_off_times_fix:
  default:
    m_beta0   = 12.*M_PI/(33.-2.*3.);
    m_Lambda2 = m_pt02*exp(-1./(m_beta0*(*as)(m_pt02)));
    m_kappa2  = m_pt02/exp(1.)*m_beta0*log(m_pt02/m_Lambda2);
    m_asmax   = (*this)(m_pt02);
  }
  //   std::ofstream was;
  //   was.open("as_in_ahadic_test.dat");
  //   for (double Q(0.9);Q<10.;Q*=1.001) {
  //     if (Q<1.) {
  //       was<<Q<<" "<<(*this)(sqr(Q))<<" --> "<<((*this)(Q*Q)/(Q*Q))<<"\n";
  //     }
  //     else {
  //       was<<Q<<" "<<(*this)(sqr(Q))<<" "<<(*as)(sqr(Q))<<" --> "<<((*this)(Q*Q)/(Q*Q))<<"\n";
  //     }
  //   }
  //   was.close();
  //   exit(1);
}

const double Strong_Coupling::operator()(double q2) const {
  if (m_form==asform::GDH_inspired) {
    double Q2(dabs(q2)),Q(sqrt(Q2));
    return m_gamma*n(Q2)/(log((Q2+mg2(Q2))/m_Lambda2));
  }
  else if (m_form==asform::constant) return m_asmax;
  else if (m_form==asform::fall_off) {
    if (q2<m_pt02) return q2/m_kappa2*exp(-q2/m_pt02);
  }
  else if (m_form==asform::fixed || 
	   m_form==asform::fall_off_times_fix) {
    static double asfix=hadpars.Get(std::string("asfix"));
    if (m_form==asform::fixed) return asfix;
    if(q2<m_pt02) return q2/m_kappa2*exp(-q2/m_pt02)*asfix;
             else return asfix/(m_beta0*log(q2/m_Lambda2));
  }
  return 1./(m_beta0*log(q2/m_Lambda2));
}

const double Strong_Coupling::n(const double Q) const {
  double crit= m_gamma/((1.+Q/m_Lambda)*log(m_m2/m_Lambda2)-m_gamma);
  //std::cout<<METHOD<<": Q/Lambda = "<<Q<<" / "<<m_Lambda
  //	   <<" and m2/Lambda2 = "<<m_m2<<" / "<<m_Lambda2<<std::endl
  //	   <<"   --> "<<((1.+Q/m_Lambda)*log(m_m2/m_Lambda2))<<" - "<<m_gamma<<std::endl;
  return M_PI*(1.+1./(crit+pow(m_b*Q,m_c)));
}

const double Strong_Coupling::mg2(const double Q) const {
  return m_m2/sqr(1.+pow(m_a*Q,m_d));
}

const double Strong_Coupling::SelectPT(const double pt2max) const {
  double pt2;
  while (true) {
    pt2 = m_pt02*(pow(pt2max/m_pt02+1.,ran.Get())-1.);
    if ((*this)(pt2)/m_asmax>ran.Get()) break;
  }
  return pt2;
}
