#include "EXTRA_XS/Two2Two/LoopHelpers.H"
#include "PHASIC++/Process/Process_Info.H"
#include "PHASIC++/Process/External_ME_Args.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "MODEL/UFO/UFO_Model.H"
#include "MODEL/Main/Model_Base.H"

#include "EXTRA_XS/Main/ME2_Base.H"

#include <vector>
#include <map>
#include <string>

using namespace EXTRAXS;
using namespace MODEL;
using namespace ATOOLS;
using namespace PHASIC;
using namespace std;


namespace EXTRAXS {
  class yy_yy : public ME2_Base {
  private:
    double m_aqed, m_aqed2, m_m2cut;
    Flavour_Vector m_SMs, m_BSMs;
    double m_s, m_t, m_u;
    
    double  Ampl2();
    Complex Ampl(vector<size_t> & hels);
    Complex ampl_f(vector<size_t> & hels, const double & m2);
    Complex ampl_v(vector<size_t> & hels, const double & m2);
    Complex ampl_s(vector<size_t> & hels, const double & m2);

    Complex pppp_f(const double & s,const double & t,const double & u,const double & m2);
    Complex mppp_f(const double & s,const double & t,const double & u,const double & m2);
    Complex mmpp_f(const double & s,const double & t,const double & u,const double & m2);
    Complex pppp_v(const double & s,const double & t,const double & u,const double & m2);
    Complex mppp_v(const double & s,const double & t,const double & u,const double & m2);
    Complex mmpp_v(const double & s,const double & t,const double & u,const double & m2);

    void Test();
  public:
    yy_yy(const External_ME_Args& args);
    double operator()(const ATOOLS::Vec4D_Vector& mom);
  };
}

DECLARE_TREEME2_GETTER(EXTRAXS::yy_yy,"yy_yy")
Tree_ME2_Base * ATOOLS::Getter<PHASIC::Tree_ME2_Base,
			       PHASIC::External_ME_Args,EXTRAXS::yy_yy>::
operator()(const External_ME_Args &args) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;
  if (MODEL::s_model->Name()!="SM") return NULL;

  const Flavour_Vector fl=args.Flavours();
  if (fl.size()!= 4) return NULL;
  if (fl[0].IsPhoton() && fl[1].IsPhoton() && fl[2].IsPhoton() && fl[3].IsPhoton()) {
    if (args.m_orders[0]==0 && args.m_orders[1]==4) {
      std::cout<<"   initialising ME.\n";
      return new yy_yy(args);
    }
  }
  return NULL;
}



yy_yy::yy_yy(const External_ME_Args& args) :
  ME2_Base(args), m_m2cut(1.e-8) {
  m_aqed = MODEL::s_model->ScalarConstant("alpha_QED"); m_aqed2 = sqr(m_aqed);
  msg_Out()<<m_aqed<<", "<<m_aqed2<<"\n";
  // To check analytic continuation and effect of mass cut Test();
  Flavour flav;
  for (size_t kfc=1;kfc<7;kfc++) {
    flav = Flavour(kfc);
    if (!flav.IsOn() || dabs(flav.Charge())<1.e-12) continue;
    if (kfc<=6) m_SMs.push_back(Flavour(kfc));
    else m_BSMs.push_back(Flavour(kfc));
  }
  for (size_t kfc=11;kfc<18;kfc+=2) {
    flav = Flavour(kfc);
    if (!flav.IsOn() || dabs(flav.Charge())<1.e-12) continue;
    if (kfc<=16) m_SMs.push_back(Flavour(kfc));
    else m_BSMs.push_back(Flavour(kfc));
  }
  m_SMs.push_back(Flavour(kf_Wplus));
  // Settings& settings = Settings::GetMainSettings();
  // will extract other things here ....
  msg_Info()<<"**** Initialising light-by-light scattering ME.\n"
	   <<"**** SM particles:\n";
  for (vector<Flavour>::iterator flit=m_SMs.begin();flit!=m_SMs.end();flit++)
    msg_Info()<<"  ** "<<(*flit)<<": mass = "<<(*flit).Mass()<<", "
	     <<"charge = "<<(*flit).Charge()<<".\n";
  if (!m_BSMs.empty()) {
    msg_Info()<<"**** BSM particles:\n";
    for (vector<Flavour>::iterator flit=m_BSMs.begin();flit!=m_BSMs.end();flit++)
      msg_Info()<<"  ** "<<(*flit)<<": mass = "<<(*flit).Mass()<<", "
	       <<"charge = "<<(*flit).Charge()<<".\n";
  }
}

double yy_yy::operator()(const ATOOLS::Vec4D_Vector& mom) {
  m_s = (mom[0]+mom[1]).Abs2();
  m_t = (mom[0]-mom[2]).Abs2();
  m_u = (mom[0]-mom[3]).Abs2();
  double me2 = Ampl2();
  msg_Out()<<"Check this: s, t, u = "<<m_s<<", "<<m_t<<", "<<m_u<<" --> "<<me2<<"\n";
  exit(1);
  return me2;
}

double yy_yy::Ampl2() {
  vector<size_t> hels; hels.resize(4,0);
  double res = 0., add;
  do {
    hels[0] += 1;
    for (size_t pos=0;pos<3;pos++) {
      if (hels[pos]==2) {
	hels[pos+1]++;
	hels[pos]=0;
      }
    }
    res += add = sqr(abs(Ampl(hels)));
  } while (hels[3]<2);
  return res;
}

Complex yy_yy::Ampl(vector<size_t> & hels) {
  Complex ampl = Complex(0.,0.), add;
  for (Flavour_Vector::iterator fit=m_SMs.begin();fit!=m_SMs.end();fit++) {
    Flavour flav = (*fit);
    if (flav.IntCharge()==0) continue;
    double q4     = pow(flav.Charge(),4);
    double m2     = sqr(flav.Mass(true));
    double strong = flav.Strong() ? dabs(flav.StrongCharge()) : 1.;  
    ampl +=  add =
      ( m_aqed2 * q4 * strong *
	( flav.IsFermion()   ? -8.*ampl_f(hels,m2) :
	  ( flav.IsVector()  ? 12.*ampl_v(hels,m2) : ampl_s(hels,m2) )
	  )
	);
  }
  for (Flavour_Vector::iterator fit=m_BSMs.begin();fit!=m_BSMs.end();fit++) {
    Flavour flav = (*fit);
    if (flav.IntCharge()==0) continue;
    double q4     = pow(flav.Charge(),4);
    double m2     = sqr(flav.Mass(true));
    double strong = flav.Strong() ? dabs(flav.StrongCharge()) : 1.;  
    ampl +=  add =
      ( m_aqed2 * q4 * strong *
	( flav.IsFermion()   ? -8.*ampl_f(hels,m2) :
	  ( flav.IsVector()  ? 12.*ampl_v(hels,m2) : ampl_s(hels,m2) )
	  )
	);
  }
  return ampl;
}

Complex yy_yy::ampl_f(vector<size_t> & hels,const double & m2) {
  // 4 equal helicities: 2 combinations
  if ( (hels[0]==1 && hels[1]==1 && hels[2]==1 && hels[3]==1) ||
       (hels[0]==0 && hels[1]==0 && hels[2]==0 && hels[3]==0) ) return pppp_f(m_s,m_t,m_u,m2);
  // 3 equal helicities: 8 combinations
  if ( (hels[0]==0 && hels[1]==1 && hels[2]==1 && hels[3]==1) ||
       (hels[0]==1 && hels[1]==0 && hels[2]==0 && hels[3]==0) ) return mppp_f(m_s,m_t,m_u,m2);
  if ( (hels[0]==1 && hels[1]==0 && hels[2]==1 && hels[3]==1) ||
       (hels[0]==0 && hels[1]==1 && hels[2]==0 && hels[3]==0) ) return mppp_f(m_s,m_u,m_t,m2);
  if ( (hels[0]==1 && hels[1]==1 && hels[2]==0 && hels[3]==1) ||
       (hels[0]==0 && hels[1]==0 && hels[2]==1 && hels[3]==0) ) return mppp_f(m_u,m_t,m_s,m2);
  if ( (hels[0]==1 && hels[1]==1 && hels[2]==1 && hels[3]==0) ||
       (hels[0]==0 && hels[1]==0 && hels[2]==0 && hels[3]==1) ) return mppp_f(m_t,m_s,m_u,m2);
  // 2 equal helicities: 6 combinations
  if ( (hels[0]==0 && hels[1]==0 && hels[2]==1 && hels[3]==1) ||
       (hels[0]==1 && hels[1]==1 && hels[2]==0 && hels[3]==0) ) return mmpp_f(m_s,m_t,m_u,m2);
  if ( (hels[0]==0 && hels[1]==1 && hels[2]==0 && hels[3]==1) ||
       (hels[0]==1 && hels[1]==0 && hels[2]==1 && hels[3]==0) ) return mmpp_f(m_u,m_t,m_s,m2);
  if ( (hels[0]==0 && hels[1]==1 && hels[2]==1 && hels[3]==0) ||
       (hels[0]==1 && hels[1]==0 && hels[2]==0 && hels[3]==1) ) return mmpp_f(m_t,m_s,m_u,m2);
  return Complex(0.,0.);
}

Complex yy_yy::pppp_f(const double & s,const double & t,const double & u,const double & m2) {
  // technical cut: small mass
  if (m2 < m_m2cut*Max(dabs(s),dabs(t))) return -Complex(1.,0.);
  Complex res = ( -Complex(1.,0.)+2.*m2*(EXTRA_XS::D02(s,t,m2) +
					 EXTRA_XS::D02(s,u,m2) +
					 EXTRA_XS::D02(u,t,m2) ) );
  return res;
}
Complex yy_yy::mppp_f(const double & s,const double & t,const double & u,const double & m2) {
  // technical cut: small mass
  if (m2 < m_m2cut*Max(dabs(s),dabs(t))) return -Complex(1.,0.);
  Complex res = ( -Complex(1.,0.) +
		  m2*(s*s+t*t+u*u) * (EXTRA_XS::C02(s,m2)/(u*t) +
				      EXTRA_XS::C02(t,m2)/(s*u) +
				      EXTRA_XS::C02(u,m2)/(s*t) ) +
		  m2*( (2*m2+s*t/u) * EXTRA_XS::D02(s,t,m2) +
		       (2*m2+s*u/t) * EXTRA_XS::D02(s,u,m2) +
		       (2*m2+t*u/s) * EXTRA_XS::D02(u,t,m2) )
	   );
  return res;
}

Complex yy_yy::mmpp_f(const double & s,const double & t,const double & u,const double & m2) {
  // technical cut: small mass
  if (m2 < m_m2cut*Max(dabs(s),dabs(t))) {
    return ( Complex(1.,0)
	     - (u-t)/s * ( EXTRA_XS::B02(u) - EXTRA_XS::B02(t) )
	     - (t*u/(s*s)-1./2.) * ( 2.*u*EXTRA_XS::C02(u) + 2.*t*EXTRA_XS::C02(t) -
				     t*u* EXTRA_XS::D02(u,t) )
	     );
  }
  Complex res = ( Complex(1.,0.) -
		  (u-t)/s * (EXTRA_XS::B02(u,m2) -
			     EXTRA_XS::B02(t,m2) ) -
		  (4*m2/s + 2*(t*u/(s*s)-1./2.) ) * (u*EXTRA_XS::C02(u,m2) +
						     t*EXTRA_XS::C02(t,m2) -
						     t*u*EXTRA_XS::D02(u,t,m2) ) + 
		  2*m2*s * (m2/s-1./2.) * (EXTRA_XS::D02(s,t,m2)+
					   EXTRA_XS::D02(s,u,m2)+
					   EXTRA_XS::D02(u,t,m2) )
		  );
  return res;
}

Complex yy_yy::ampl_v(vector<size_t> & hels,const double & m2) {
  // 4 equal helicities: 2 combinations
  if ( (hels[0]==1 && hels[1]==1 && hels[2]==1 && hels[3]==1) ||
       (hels[0]==0 && hels[1]==0 && hels[2]==0 && hels[3]==0) ) return pppp_f(m_s,m_t,m_u,m2);
  // 3 equal helicities: 8 combinations
  if ( (hels[0]==0 && hels[1]==1 && hels[2]==1 && hels[3]==1) ||
       (hels[0]==1 && hels[1]==0 && hels[2]==0 && hels[3]==0) ) return mppp_f(m_s,m_t,m_u,m2);
  if ( (hels[0]==1 && hels[1]==0 && hels[2]==1 && hels[3]==1) ||
       (hels[0]==0 && hels[1]==1 && hels[2]==0 && hels[3]==0) ) return mppp_f(m_s,m_u,m_t,m2);
  if ( (hels[0]==1 && hels[1]==1 && hels[2]==0 && hels[3]==1) ||
       (hels[0]==0 && hels[1]==0 && hels[2]==1 && hels[3]==0) ) return mppp_f(m_u,m_t,m_s,m2);
  if ( (hels[0]==1 && hels[1]==1 && hels[2]==1 && hels[3]==0) ||
       (hels[0]==0 && hels[1]==0 && hels[2]==0 && hels[3]==1) ) return mppp_f(m_t,m_s,m_u,m2);
  // 2 equal helicities: 6 combinations
  if ( (hels[0]==0 && hels[1]==0 && hels[2]==1 && hels[3]==1) ||
       (hels[0]==1 && hels[1]==1 && hels[2]==0 && hels[3]==0) ) return mmpp_f(m_s,m_t,m_u,m2);
  if ( (hels[0]==0 && hels[1]==1 && hels[2]==0 && hels[3]==1) ||
       (hels[0]==1 && hels[1]==0 && hels[2]==1 && hels[3]==0) ) return mmpp_f(m_u,m_t,m_s,m2);
  if ( (hels[0]==0 && hels[1]==1 && hels[2]==1 && hels[3]==0) ||
       (hels[0]==1 && hels[1]==0 && hels[2]==0 && hels[3]==1) ) return mmpp_f(m_t,m_s,m_u,m2);
  return Complex(0.,0.);
}

Complex yy_yy::pppp_v(const double & s,const double & t,const double & u,const double & m2) {
  return pppp_f(s,t,u,m2);
}

Complex yy_yy::mppp_v(const double & s,const double & t,const double & u,const double & m2) {
  return mppp_f(s,t,u,m2);
}

Complex yy_yy::mmpp_v(const double & s,const double & t,const double & u,const double & m2) {
  // technical cut: small mass
  if (m2 < m_m2cut*Max(dabs(s),dabs(t))) {
  }
  Complex res = ( Complex(1.,0)
		  - (u-t)/s * ( EXTRA_XS::B02(u,m2) - EXTRA_XS::B02(t,m2) )
		  - (t*u/(s*s)-4./3.+4.*m2/s) * ( 2.*u*EXTRA_XS::C02(u,m2) + 2.*t*EXTRA_XS::C02(t,m2) -
						  t*u*EXTRA_XS::D02(u,t,m2) ) +
		  ( 2.*m2*s*(m2/s-4./3.)+2./3.*s*s ) * (EXTRA_XS::D02(s,t) +
							EXTRA_XS::D02(s,u) +
							EXTRA_XS::D02(u,t) )
		  );
  return res;  
}


Complex yy_yy::ampl_s(vector<size_t> & hels,const double & m2) {
  return Complex(0.,0.);  
}

void yy_yy::Test() {
  map<string, Histogram * > histograms;
  size_t nbins = 80;
  double expmin = -12., expmax = -4.;
  histograms["M2_2_Zero_electron"] = new Histogram(0,expmin,expmax,nbins);
  histograms["M2_2_Zero_muon"]     = new Histogram(0,expmin,expmax,nbins);
  m_s = 100.;
  m_t = -10.;
  m_u = -m_s-m_t;
  double me2_e, me2_mu;
  double binsize = (expmax-expmin)/double(nbins), expo;
  for (size_t i=0;i<nbins;i++) {
    expo    = expmin+(expmax-expmin)*double(i)/double(nbins);
    m_m2cut = pow(10.,expo);
    m_SMs.push_back(Flavour(kf_e));
    me2_e   = Ampl2();
    m_SMs.clear();
    m_SMs.push_back(Flavour(kf_mu));
    me2_mu  = Ampl2();
    m_SMs.clear();
    histograms["M2_2_Zero_electron"]->Insert(expo+binsize/2.,me2_e);
    histograms["M2_2_Zero_muon"]->Insert(expo+binsize/2.,me2_mu);
  }
  Histogram * histogram;
  string name;
  for (map<string,Histogram *>::iterator 
	 hit=histograms.begin();hit!=histograms.end();hit++) {
    histogram = hit->second;
    name  = hit->first+string(".dat");
    histogram->Output(name);
    delete histogram;
  }
  histograms.clear();
  exit(1);
}

