#include "AMISIC++/Tools/Hadronic_XSec_Calculator.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace AMISIC;
using namespace ATOOLS;
using namespace std;

Hadronic_XSec_Calculator::Hadronic_XSec_Calculator() :
  m_Ecms(rpa->gen.Ecms()), m_s(m_Ecms*m_Ecms),
  m_massp(Flavour(kf_p_plus).Mass()), m_masspi(Flavour(kf_pi).Mass()),
  m_pomeron(0.0808), m_reggeon(-0.4525),m_slope(2.3),
  m_xsecpom(21.70),
  m_xsecregge(rpa->gen.Beam1().IsAnti()^rpa->gen.Beam2().IsAnti()?98.39:56.08)
{}
  
void Hadronic_XSec_Calculator::operator()()
{
  m_xstot = CalculateTotalXSec();
  m_xsel  = CalculateElasticXSec(m_xstot);
  m_xssd  = CalculateSingleDXSec();
  m_xsdd  = CalculateDoubleDXSec();

  m_xsnd  = m_xstot-m_xsel-2.0*m_xssd-m_xsdd;
  msg_Info()<<METHOD<<": Results are {\n"
	    <<"   \\sigma_{tot} = "<<m_xstot<<" mb\n"
	    <<"   \\sigma_{el}  = "<<m_xsel<<" mb\n"
	    <<"   \\sigma_{sd}  = "<<2.0*m_xssd<<" mb\n"
	    <<"   \\sigma_{dd}  = "<<m_xsdd<<" mb\n"
	    <<"   \\sigma_{nd}  = "<<m_xsnd<<" mb = "
	    <<(m_xsnd*1.e9/rpa->Picobarn())<<" GeV^-2\n}"<<std::endl;
  // convert all cross sections to 1/GeV^2
  m_xstot *= 1.e9/rpa->Picobarn();
  m_xsel  *= 1.e9/rpa->Picobarn();
  m_xsnd  *= 1.e9/rpa->Picobarn();
  m_xssd  *= 1.e9/rpa->Picobarn();
  m_xsdd  *= 1.e9/rpa->Picobarn();
}

double Hadronic_XSec_Calculator::CalculateTotalXSec() {
  // standard two-component fit: pomeron + reggeon
  return m_xsecpom*pow(m_s,m_pomeron)+m_xsecregge*pow(m_s,m_reggeon);
}

double Hadronic_XSec_Calculator::CalculateElasticXSec(const double & xstot) {
  // standard two-component fit: pomeron + reggeon
  return 0.0511*xstot*xstot/(4.*(m_slope+pow(m_s,m_pomeron))-4.2);
}


double Hadronic_XSec_Calculator::CalculateSingleDXSec() {
  double ap   = 0.25;
  double mmin = m_massp+2.*m_masspi, mmin2 = sqr(mmin), mmax2 = 0.213*m_s;
  double cres = 2., bax = -0.47+150./m_s, mres2 = 2.;
  double JAX  =
    0.5/ap*log((m_slope+ap*log(m_s/mmin2))/(m_slope+ap*log(m_s/mmax2))) +
    0.5*cres/(m_slope+ap*log(m_s/mmin2)+bax)*log(1.0+sqr(mres2/mmin2));
  return 0.0336*pow(m_xsecpom,1.5)*JAX;
}

double Hadronic_XSec_Calculator::CalculateDoubleDXSec() {
  double ap = 0.25, Del0 = 3.2-9.0/log(m_s)+17.4/sqr(log(m_s)), s0 = 8.;
  double y0 = log(m_s/sqr(m_massp)), ymin = 4.0*log(1.0+2.0*m_masspi/m_massp);

  double mmin1 = m_massp+2.*m_masspi, mmin12 = sqr(mmin1);
  double mmin2 = mmin1, mmin22 = sqr(mmin2);
  double mmax2 = 0.213*m_s, mres1 = 2., mres2 = 2.;
  double cres  = 2.;
  double mmxxx = m_s*(0.07-0.44/log(m_s)+1.36/sqr(log(m_s)));
  double bxx   = -1.05+40./sqrt(m_s)+8000./sqr(m_s);
  double JXX  =
    0.5/ap*((y0-ymin)*(log((y0-ymin)/Del0)-1.0)+Del0) +
    0.5*cres/ap*log(1.0+sqr(mres2/mmin2))*log(log(m_s*s0/(mmin12*mmin2*mres2))/
					      log(m_s*s0/(mmxxx*mres2*mmin2))) +
    0.5*cres/ap*log(1.0+sqr(mres1/mmin1))*log(log(m_s*s0/(mmin22*mmin1*mres1))/
					      log(m_s*s0/(mmxxx*mres1*mmin1))) +
    sqr(cres)/(2.0*ap*log(m_s*s0/(mres1*mres2*mmin1*mmin2))+bxx)*
    log(1.0+sqr(mres1/mmin1))*log(1.0+sqr(mres2/mmin2));
  return 0.0084*m_xsecpom*JXX;
}


