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
  m_eta(hadpars.Get(std::string("pt_exponent"))),
  m_a(3.008), m_b(1.425), m_c(0.908), m_d(0.84), m_m2(sqr(1.204))
{
  MODEL::Running_AlphaS * as = static_cast<MODEL::Running_AlphaS *>
    (MODEL::s_model->GetScalarFunction(std::string("alpha_S")));
  switch (m_form) {
  case asform::IRregularised:
    m_beta0   = 12.*M_PI/(33.-2.*3.);
    m_Lambda2 = m_pt02*exp(-m_beta0/(*as)(m_pt02));
  case asform::GDH_inspired:
    m_asmax   = (*this)(0.);
    if (m_asmax<0.) {
      msg_Error()<<"Error in "<<METHOD<<":"<<std::endl
		 <<"   Maximal alphaS too small for pt_0^2 = "
		 <<m_pt02<<": "<<m_asmax<<"."<<std::endl
		 <<"   Will abort the run."<<std::endl;
      abort();
    }
    break;
  case asform::IRregularised_IR0:
    m_beta0   = 12.*M_PI/(33.-2.*3.);
    m_Lambda2 = m_pt02*exp(-m_beta0/(*as)(m_pt02));
    m_asmax   = (*as)(m_pt02); 
    break;
  case asform::constant:
  default:
    m_asmax = hadpars.Get(std::string("asfix"));
    break;
  }
  return;

  std::ofstream was;
  was.open("as_in_ahadic_test.dat");
  was<<"asmax for pt_0^2 = "<<m_pt02<<": "<<m_asmax
     <<" for Lambda^2 = "<<m_Lambda2<<"."<<std::endl;
  for (double Q(0.0);Q<0.1;Q+=.001) {
    was<<Q<<" "<<(*this)(sqr(Q))<<"\n";
  }
  for (double Q(0.1);Q<10.;Q*=1.001) {
    if (Q<1.) {
      was<<Q<<" "<<(*this)(sqr(Q))<<"\n";
    }
    else {
      //m_beta0 = 12.*M_PI/(33.-2.*4.);
      was<<Q<<" "<<(*this)(sqr(Q))<<" "<<(*as)(sqr(Q))<<"\n";
    }
  }
  for (double Q(10.);Q<100.;Q*=1.1) {
    //m_beta0 = 12.*M_PI/(33.-2.*5.);
    was<<Q<<" "<<(*this)(sqr(Q))<<" "<<(*as)(sqr(Q))<<"\n";
  }
  was.close();
  exit(1);
}

const double Strong_Coupling::operator()(double q2,bool reweight) const {
  double Q2(dabs(q2)), Q; 

  switch (m_form) {
  case asform::IRregularised_IR0:
    if (reweight) return m_beta0/log((Q2+m_pt02)/m_Lambda2) * Q2/(Q2+m_pt02);
  case asform::IRregularised:
    return m_beta0/log((Q2+m_pt02)/m_Lambda2);
  case asform::GDH_inspired:
    Q = sqrt(Q2);
    return m_gamma*n(Q)/(log((Q2+mg2(Q))/m_Lambda2));
  case asform::constant: 
    return m_asmax;
  }
  return 1./(m_beta0*log(Q2/m_Lambda2));
}

const double Strong_Coupling::SelectPT(const double & scale2max,const double & scale2min) {
  double pt2(0.);
  double mini(m_pt02+scale2min),maxi(m_pt02+scale2max),expo(m_eta==1.?0.:1./1.-m_eta);;
  bool   runit(true);
  while (runit) {
    switch (m_form) {
    case asform::Exponential:
      pt2 = -m_pt02*log(exp(-scale2min/m_pt02)+ran.Get()*(exp(-scale2max/m_pt02)-exp(-scale2min/m_pt02)));
      runit=false;
      break;
    case asform::IRregularised_IR0: 
      if (scale2max<=m_pt02) {
	pt2 = scale2min+(scale2max-scale2min)*sqrt(ran.Get());
	if ((*this)(pt2,false)/m_asmax* sqr(m_pt02/(m_pt02+pt2))>ran.Get()) runit = false;
      }
      else {
	pt2 = -m_pt02+mini*pow(maxi/mini,ran.Get());
	if ((*this)(pt2,false)/m_asmax>ran.Get()) runit = false;
      }
      break;
    case asform::IRregularised: 
    case asform::GDH_inspired:
    case asform::constant: 
    default:
      if (m_eta==1.) pt2 = -m_pt02+mini*pow(maxi/mini,ran.Get());
      else {
	double rn(ran.Get());
	pt2 = -m_pt02+pow(mini*(1.-rn)+maxi*rn,expo);
      }
      if ((*this)(pt2)/m_asmax>ran.Get()) runit = false;
      break;
    }
  }
  return m_lastpt2 = pt2;
}

const double Strong_Coupling::n(const double Q) const {
  double crit= m_gamma/((1.+Q/m_Lambda)*log(m_m2/m_Lambda2)-m_gamma);
  return M_PI*(1.+1./(crit+pow(m_b*Q,m_c)));
}

const double Strong_Coupling::mg2(const double Q) const {
  return m_m2/sqr(1.+pow(m_a*Q,m_d));
}

