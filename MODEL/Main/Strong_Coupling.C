#include "MODEL/Main/Strong_Coupling.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"

using namespace MODEL;
using namespace ATOOLS;


Strong_Coupling::Strong_Coupling(Running_AlphaS * as,const asform::code & asf,
				 const double & pt02) :
  m_form(asf), p_as(as), m_pt02(pt02)
{
  switch (m_form) {
  case asform::constant:
  case asform::frozen:
  case asform::smooth:
  case asform::IR0:
    m_asmax   = (*p_as)(m_pt02); 
    break;
  case asform::GDH_inspired:
    m_beta0   = 12.*M_PI/(33.-2.*4.);
    m_Lambda  = 0.349; 
    m_Lambda2 = sqr(m_Lambda);
    m_gamma   = m_beta0/M_PI;
    m_a       = 3.008;
    m_b       = 1.425;
    m_c       = 0.908;
    m_d       = 0.84;
    m_m2      = sqr(1.204);
    m_asmax   = (*this)(0.);
    if (m_asmax<0.) {
      msg_Error()<<"Error in "<<METHOD<<":"<<std::endl
		 <<"   Maximal alphaS too small for pt_0^2 = "
		 <<m_pt02<<": "<<m_asmax<<"."<<std::endl
		 <<"   Will abort the run."<<std::endl;
      abort();
    }
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

double Strong_Coupling::operator()(double q2,bool reweight) const {
  double Q2(dabs(q2)), Q; 

  switch (m_form) {
  case asform::constant:
    return m_asmax;
  case asform::frozen:
    if (Q2<m_pt02) return m_asmax;
    return (*p_as)(Q2);
  case asform::smooth:
    return (*p_as)(Q2+m_pt02);
  case asform::IR0:
    if (Q2<m_pt02) return m_asmax*Q2/m_pt02;
    return (*p_as)(Q2);
  case asform::GDH_inspired:
    Q = sqrt(Q2);
    return m_gamma*n(Q)/(log((Q2+mg2(Q))/m_Lambda2));
  }
  return m_asmax;
}

/*
double Strong_Coupling::
SelectPT(const double & scale2max,const double & scale2min,
	 const bool & expo) {
  double pt2(0.), pt2max(Min(scale2max,m_pt2max)), ran1;
  double mini(m_pt02+scale2min),maxi(m_pt02+pt2max);
  bool   runit(true);
  while (runit) {
    ran1 = ran->Get();
    switch (m_form) {
    case asform::IRregularised_IR0: 
      if (pt2max<=m_pt02) {
	pt2 = scale2min+(pt2max-scale2min)*sqrt(ran1);
	if ((*this)(pt2,false)/m_asmax * sqr(m_pt02/(m_pt02+pt2)) *
	    (expo?exp(-4.*pt2/m_pt02):1.)>ran->Get()) 
	  runit = false;
      }
      else {
	pt2 = -m_pt02+mini*pow(maxi/mini,ran1);
	if ((*this)(pt2,false)/m_asmax * pt2/(m_pt02+pt2) *
	    (expo?exp(-4.*pt2/m_pt02):1.)> ran->Get()) 
	  runit = false;
      }
      break;
    case asform::IRregularised: 
    case asform::GDH_inspired:
    case asform::constant: 
    default:
      pt2 = -m_pt02+mini*pow(maxi/mini,ran1);
      if ((*this)(pt2)/m_asmax>ran->Get()) runit = false;
      break;
    }
  }
  //msg_Out()<<METHOD<<"("<<m_form<<", "<<scale2max<<" -> "
  //	   <<"pt^2_max = "<<pt2max<<", "
  //	   <<"pt0^2 = "<<m_pt02<<") = "<<pt2<<".\n";
  return m_lastpt2 = pt2;
}
*/

double Strong_Coupling::n(const double Q) const {
  double crit= m_gamma/((1.+Q/m_Lambda)*log(m_m2/m_Lambda2)-m_gamma);
  return M_PI*(1.+1./(crit+pow(m_b*Q,m_c)));
}

double Strong_Coupling::mg2(const double Q) const {
  return m_m2/sqr(1.+pow(m_a*Q,m_d));
}

