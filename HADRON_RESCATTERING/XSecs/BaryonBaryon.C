#include "HADRON_RESCATTERING/XSecs/BaryonBaryon.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Org/Message.H"

using namespace HADRON_RESCATTERING;
using namespace ATOOLS;
using namespace std;

BaryonBaryon::BaryonBaryon() :
  m_test(true) {
  if (m_test) { Tests(); exit(1); }
}

BaryonBaryon::~BaryonBaryon() {}

double BaryonBaryon::
XStot(const ATOOLS::Flavour & A,const ATOOLS::Flavour & B,
      const double & s) {
  if (!(A.IsBaryon() && B.IsBaryon())) return 0.;
  if ((A.Kfcode()==2212 && B.Kfcode()==2212) ||
      (A.Kfcode()==2112 && B.Kfcode()==2112)) {
    if ((A.IsAnti() && !B.IsAnti()) || (!A.IsAnti() && B.IsAnti()))
      return m_NN.ppbartot(s); 
    return m_NN.pptot(s);
  }
  if ((A.Kfcode()==2212 && B.Kfcode()==2112) ||
      (A.Kfcode()==2112 && B.Kfcode()==2212)) {
    return m_NN.pntot(s);
  }
  m_s   = s;
  m_mA = A.HadMass(); m_mA2 = sqr(m_mA);
  m_mB = B.HadMass(); m_mB2 = sqr(m_mB);
  return 0.;
}

double BaryonBaryon::
XSel(const ATOOLS::Flavour & A,const ATOOLS::Flavour & B,
     const double & s) {
  if (!(A.IsBaryon() && B.IsBaryon())) return 0.;
  if ((A.Kfcode()==2212 && B.Kfcode()==2212) ||
      (A.Kfcode()==2112 && B.Kfcode()==2112)) {
    if ((A.IsAnti() && !B.IsAnti()) || (!A.IsAnti() && B.IsAnti()))
      return m_NN.ppbarel(s); 
    return m_NN.ppel(s);
  }
  if ((A.Kfcode()==2212 && B.Kfcode()==2112) ||
      (A.Kfcode()==2112 && B.Kfcode()==2212)) {
    return m_NN.pnel(s);
  }
  return 0.;
}

double BaryonBaryon::xstot(long int & A, long int & B,const double & plab) {
  return 0.;
}
  
double BaryonBaryon::xsel(long int & A, long int & B, const double & plab) {
  return 0.;
}
void BaryonBaryon::Tests() {
  size_t bins = 10000;
  double pmin = 0., pmax = 10., pinc = (pmax-pmin)/double(bins);
  map<string,Histogram *>  histos;
  histos["pp_total_low"]      = new Histogram(0,pmin,pmax,bins);
  histos["pp_elastic_low"]    = new Histogram(0,pmin,pmax,bins);
  histos["pn_total_low"]      = new Histogram(0,pmin,pmax,bins);
  histos["pn_elastic_low"]    = new Histogram(0,pmin,pmax,bins);
  histos["ppbar_total_low"]   = new Histogram(0,pmin,pmax,bins);
  histos["ppbar_elastic_low"] = new Histogram(0,pmin,pmax,bins);
  histos["ppbar_annihil_low"] = new Histogram(0,pmin,pmax,bins);
  histos["ppbar_CEX_low"]     = new Histogram(0,pmin,pmax,bins);
  histos["pnbar_annihil_low"] = new Histogram(0,pmin,pmax,bins);
  Flavour flA(kf_p_plus), flB = flA;
  double  plab, s;
  for (int i=0;i<bins;i++) {
    plab   = pmin+i*pinc;
    flB    = flA;
    s      = ( sqr(flA.HadMass()) + sqr(flB.HadMass()) +
	       2.*flA.HadMass()*sqrt(sqr(flB.HadMass())+sqr(plab)) );
    histos["pp_total_low"]->Insert(i,   XStot(flA,flB,s));
    histos["pp_elastic_low"]->Insert(i, XSel(flA,flB,s));
    flB    = flA.Bar();
    histos["ppbar_total_low"]->Insert(i,   XStot(flA,flB,s));
    histos["ppbar_elastic_low"]->Insert(i, XSel(flA,flB,s));
    histos["ppbar_annihil_low"]->Insert(i, m_NN.ppbarAnnihil(s));
    histos["ppbar_CEX_low"]->Insert(i, m_NN.ppbarCEX(s));
    flB    = Flavour(kf_n);
    s      = ( sqr(flA.HadMass()) + sqr(flB.HadMass()) +
	       2.*flA.HadMass()*sqrt(sqr(flB.HadMass())+sqr(plab)) );
    histos["pn_total_low"]->Insert(i,   XStot(flA,flB,s));
    histos["pn_elastic_low"]->Insert(i, XSel(flA,flB,s));
    flB    = flB.Bar();
    histos["pnbar_annihil_low"]->Insert(i, m_NN.pnbarAnnihil(s));
  }
  Histogram * histo;
  std::string name;
  for (std::map<std::string,Histogram *>::iterator
	 hit=histos.begin();hit!=histos.end();hit++) {
    histo = hit->second;
    name  = std::string("XSecs/")+hit->first+std::string(".dat");
    histo->Output(name);
    delete histo;
  }
  histos.clear();
}


