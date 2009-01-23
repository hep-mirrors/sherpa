#include "Strong_Coupling.H"
#include "Hadronisation_Parameters.H"
#include "Running_AlphaS.H"
#include "Random.H"
#include "Model_Base.H"
#include "Message.H"

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
  case asform::IR_cutoff:
    m_asmax   = (*this)(0.);
    if (m_asmax<0.) {
      msg_Error()<<"Error in "<<METHOD<<":"<<std::endl
		 <<"   Maximal alphaS too small for pt_0^2 = "<<m_pt02<<": "<<m_asmax<<"."<<std::endl
		 <<"   Will abort the run."<<std::endl;
      abort();
    }
    break;
  case asform::GDH_inspired:
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
  return;

  std::ofstream was;
  was.open("as_in_ahadic_test.dat");
  for (double Q(0.0);Q<0.1;Q+=.001) {
    was<<Q<<" "<<(*this)(sqr(Q))<<"\n";
  }
  for (double Q(0.1);Q<10.;Q*=1.001) {
    if (Q<1.) {
      was<<Q<<" "<<(*this)(sqr(Q))<<"\n";
    }
    else {
      was<<Q<<" "<<(*this)(sqr(Q))<<" "<<(*as)(sqr(Q))<<"\n";
    }
  }
  was.close();
  exit(1);
}

const double Strong_Coupling::operator()(double q2) const {
  double Q2(dabs(q2)), Q(sqrt(Q2));

  switch (m_form) {
  case asform::IR_cutoff:
    return 1./(m_beta0*log((Q2+m_pt02)/m_Lambda2));
  case asform::GDH_inspired:
    return m_gamma*n(Q2)/(log((Q2+mg2(Q2))/m_Lambda2));
  case asform::constant: 
    return m_asmax;
  case asform::fall_off:
    if (Q2<m_pt02) return Q2/m_kappa2*exp(-Q2/m_pt02);
    break;
  case asform::fixed:
  case asform::fall_off_times_fix:
  default:
    static double asfix=hadpars.Get(std::string("asfix"));
    if (m_form==asform::fixed) return asfix;
    if(Q2<m_pt02) return Q2/m_kappa2*exp(-Q2/m_pt02)*asfix;
    else return asfix/(m_beta0*log(q2/m_Lambda2));
  }
  return 1./(m_beta0*log(Q2/m_Lambda2));
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

const double Strong_Coupling::SelectPT(const double pt2max,const bool pole) {
  double pt2;
  while (true) {
    if (pole) pt2 = m_pt02*(pow(pt2max/m_pt02+1.,ran.Get())-1.);
         else pt2 = ran.Get()*pt2max;
    if ((*this)(pt2)/m_asmax>ran.Get()) break;
  }
  return m_lastpt2 = pt2;
}
