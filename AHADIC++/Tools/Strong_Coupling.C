#include "AHADIC++/Tools/Strong_Coupling.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "ATOOLS/Math/Random.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Message.H"

using namespace AHADIC;
using namespace ATOOLS;


Strong_Coupling::Strong_Coupling(const asform::code form) :
  m_form(form), 
  m_beta0(12.*M_PI/(33.-2.*4.)), m_Lambda(0.349), m_Lambda2(sqr(m_Lambda)),
  m_pt02(dabs(hadpars.Get(std::string("pt02")))), 
  m_pt2max(10000.), m_lastpt2(0.), m_asmax(1.), 
  m_eta(hadpars.Get(std::string("pt_exponent"))),
  m_gamma(m_beta0/M_PI), m_a(3.008), m_b(1.425), m_c(0.908), m_d(0.84), m_m2(sqr(1.204))
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
    was<<Q<<" "<<(*this)(sqr(Q),true)<<"\n";
  }
  for (double Q(0.1);Q<10.;Q*=1.001) {
    if (Q<1.) {
      was<<Q<<" "<<(*this)(sqr(Q),true)<<"\n";
    }
    else {
      //m_beta0 = 12.*M_PI/(33.-2.*4.);
      was<<Q<<" "<<(*this)(sqr(Q),true)<<" "<<(*as)(sqr(Q))<<"\n";
    }
  }
  for (double Q(10.);Q<100.;Q*=1.1) {
    //m_beta0 = 12.*M_PI/(33.-2.*5.);
    was<<Q<<" "<<(*this)(sqr(Q),true)<<" "<<(*as)(sqr(Q))<<"\n";
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
  double pt2(0.), pt2max(Min(scale2max,m_pt2max)), ran1;
  double mini(m_pt02+scale2min),maxi(m_pt02+pt2max),expo(dabs(m_eta-1.)<0.01?0.:1.-m_eta), invexpo(1./expo);
  double mini_pow(pow(mini,expo)), maxi_pow(pow(maxi,expo));
  bool   runit(true);
  while (runit) {
    ran1 = ran.Get();
    switch (m_form) {
    case asform::IRregularised_IR0: 
      if (m_eta==1.) {
	if (pt2max<=m_pt02) {
	  pt2 = scale2min+(pt2max-scale2min)*sqrt(ran1);
	  if ((*this)(pt2,false)/m_asmax* sqr(m_pt02/(m_pt02+pt2))>ran.Get()) runit = false;
	}
	else {
	  pt2 = -m_pt02+mini*pow(maxi/mini,ran1);
	  if ((*this)(pt2,false)/m_asmax * pt2/(m_pt02+pt2) > ran.Get()) runit = false;
	}
      }
      else {
	if (pt2max<=m_pt02) {
	  pt2 = pow(pow(scale2min,1.+m_eta)*(1.-ran1)+pow(pt2max,1.+m_eta)*ran1,1./(1.+m_eta));
	  //std::cout<<"Try this for pt2max = "<<pt2max<<", pt2min = "<<scale2min
	  //	   <<" --> "<<pt2<<"."<<std::endl;
	  if ((*this)(pt2,false)/m_asmax* pow(m_pt02/(m_pt02+pt2),1.+m_eta)>ran.Get()) runit = false;
	}
	else {
	  pt2 = -m_pt02+mini*pow(maxi/mini,ran1);
	  //std::cout<<"Try this for pt2max = "<<pt2max<<", pt2min = "<<scale2min
	  //	   <<" --> "<<pt2<<"."<<std::endl;
	  if ((*this)(pt2,false)/m_asmax * pow(pt2/(m_pt02+pt2),m_eta) > ran.Get()) runit = false;
	}
      }
      break;
    case asform::IRregularised: 
    case asform::GDH_inspired:
    case asform::constant: 
    default:
      if (m_eta==1.) 
	pt2 = -m_pt02+mini*pow(maxi/mini,ran1);
      else 
	pt2 = -m_pt02+pow(ran1*maxi_pow+(1.-ran1)*mini_pow,invexpo);
      if ((*this)(pt2)/m_asmax>ran.Get()) runit = false;
      break;
    }
  }
  //std::cout<<"In "<<METHOD<<" with "<<scale2max<<" --> "<<pt2max
  //	   <<" leads to pt2 = "<<pt2<<"."<<std::endl;
  if (pt2>400.) 
    msg_Tracking()<<"Surprise in "<<METHOD<<": pt = "<<sqrt(pt2)
		  <<" from pt_min,max = "<<sqrt(scale2min)<<", "<<sqrt(scale2max)
		  <<" --> "<<sqrt(m_pt2max)<<std::endl
		  <<"    as(pt^2) = "<<(*this)(pt2,false)<<", asmax = "<<m_asmax
		  <<" ("<<int(m_form)<<", eta = "<<m_eta<<", pt0^2 = "<<m_pt02
		  <<", random = "<<ran1<<")."<<std::endl;
  return m_lastpt2 = pt2;
}

const double Strong_Coupling::n(const double Q) const {
  double crit= m_gamma/((1.+Q/m_Lambda)*log(m_m2/m_Lambda2)-m_gamma);
  return M_PI*(1.+1./(crit+pow(m_b*Q,m_c)));
}

const double Strong_Coupling::mg2(const double Q) const {
  return m_m2/sqr(1.+pow(m_a*Q,m_d));
}

