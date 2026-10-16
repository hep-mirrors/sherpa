#include "HADRON_RESCATTERING/XSecs/BaryonBaryon.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Poincare.H"
using namespace HADRON_RESCATTERING;
using namespace ATOOLS;
using namespace std;

BaryonBaryon::BaryonBaryon() : ScatteringBase()
 {
  m_test = false;
  if (m_test) 
  { 
    Tests();
  } 
}

BaryonBaryon::~BaryonBaryon() {msg_Out() << "BaryonBaryon::Destructor called" << std::endl;}

double BaryonBaryon::
XStot(const ATOOLS::Flavour & A,const ATOOLS::Flavour & B,
      const double & s) {
  //Scoped_Settings s{Settings::GetMainSettings()["FormFactors"]};
  //m_bb_xstot = s["BBScattering"].SetDefault(scatmodel::off).Get<scatmodel::code>();
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

double BaryonBaryon::pnToDGamma(ATOOLS::Particle * A , ATOOLS::Particle * B) {


  if ((A->Flav().Kfcode()==2212 && B->Flav().Kfcode() ==2112) || (A->Flav().Kfcode()==2112 && B->Flav().Kfcode()==2212)) //pn pair.
  {
    std::vector<double> aVec ={2.30346, -9.366346*10, 2.565390*std::pow(10,3), 
                               -2.5594101*std::pow(10,4), 1.43513109*std::pow(10,5),
                               -5.0357289*std::pow(10,5), 1.14924802*std::pow(10,6), 
                               -1.72368391*std::pow(10,6), 1.67934876*std::pow(10,6),
                               -1.01988855*std::pow(10,6), 3.4984035*std::pow(10,5),
                               -5.1662760*std::pow(10,4)};
    std::vector<double> bVec = {-5.1885, 2.9196};

    Vec4D p1 = A->Momentum(); Vec4D p2 = B->Momentum();
    Vec4D P = p1 + p2;    Poincare boost(P);
    Vec4D p1cm = p1;  boost.Boost(p1cm);
    // k = |p1 - p2| in the c.m. frame, i.e. twice the c.m. momentum of each
    // nucleon (arXiv:1504.07242, after Eq. (1)); the fits are in kappa = k/GeV.
    double p_rel = 2.*sqrt(p1cm[1]*p1cm[1]+p1cm[2]*p1cm[2]+p1cm[3]*p1cm[3]);
    double sigma = 0; 

if (p_rel < 1.28) {
    double kappa = p_rel; // kappa = k/(1 GeV), already dimensionless
    sigma = aVec[0] / kappa;          // a_{-1} * kappa^{-1}
    double kp = 1.0;
    for (size_t n = 1; n < aVec.size(); ++n) {
        sigma += aVec[n] * kp;        // a_{n-1} * kappa^{n-1}
        kp *= kappa;
    }
}
else {
    // note: b1=-5.1885, b2=2.9196, formula is exp(-b1*kappa - b2*kappa^2)
    sigma = std::exp(-bVec[0]*p_rel - bVec[1]*p_rel*p_rel);
}
    return sigma;

  }
  return 0.;
}

double BaryonBaryon::NNToDPi(ATOOLS::Particle * A , ATOOLS::Particle * B) {
  double a = 170., b = 1.34, c= 1.77, d = 0.38, e = 0.096;// a is in micro barns

  Vec4D p1 = A->Momentum(); Vec4D p2 = B->Momentum();
  Vec4D P = p1 + p2;    Poincare boost(P);
  Vec4D p1cm = p1;  boost.Boost(p1cm);
  double p_rel = sqrt(p1cm[1]*p1cm[1]+p1cm[2]*p1cm[2]+p1cm[3]*p1cm[3]);
  double s     = P[0]*P[0]-P[1]*P[1]-P[2]*P[2]-P[3]*P[3]; double sqrtS = std::sqrt(s);
  double sThreshold = std::pow(mDeuteron + mPion, 2);
  if (s < sThreshold) return 0.;
  double kallen = (s - std::pow(mDeuteron + mPion, 2)) 
              * (s - std::pow(mDeuteron - mPion, 2));
  double q = (kallen > 0.) ? std::sqrt(kallen) / (2.0 * sqrtS) : 0.0;

  double etaDefined = (q)/mPion;

  if (q ==0. ) return 0.;

  double denom = (c-std::exp(d*etaDefined));
  double sigma = a * std::pow(etaDefined,b)/(denom*denom + e);
  return sigma;
}

double BaryonBaryon::nnToDPiMinus(ATOOLS::Particle * A , ATOOLS::Particle * B) {
  return NNToDPi(A,B);
}
double BaryonBaryon::pnToDPiZero(ATOOLS::Particle * A , ATOOLS::Particle * B) {
  return 0.5*NNToDPi(A,B);
}
double BaryonBaryon::ppToDPiPlus(ATOOLS::Particle * A , ATOOLS::Particle * B) {
  return NNToDPi(A,B);
}

double BaryonBaryon::ppToPiPlusPiZero(ATOOLS::Particle * A , ATOOLS::Particle * B) {
  std::vector<double> piPlusPiZero = {5.099*std::pow(10,15), 1.656*10, 2.333*std::pow(10,7), 1.133*10, 2.868*std::pow(10,16)};

  Vec4D p1 = A->Momentum(); Vec4D p2 = B->Momentum();
  Vec4D P = p1 + p2;    Poincare boost(P);
  Vec4D p1cm = p1;  boost.Boost(p1cm);
  // k = |p1 - p2| in the c.m. frame, i.e. twice the c.m. momentum of each
  // nucleon (arXiv:1504.07242, after Eq. (1)); the fits are in kappa = k/GeV.
  double p_rel = 2.*sqrt(p1cm[1]*p1cm[1]+p1cm[2]*p1cm[2]+p1cm[3]*p1cm[3]);
  // double s     = P[0]*P[0]-P[1]*P[1]-P[2]*P[2]-P[3]*P[3]; double sqrtS = std::sqrt(s);
  double sigma = 0;
  sigma = piPlusPiZero[0]*std::pow(p_rel,piPlusPiZero[1])/(std::pow( (piPlusPiZero[2] - std::exp(piPlusPiZero[3]*p_rel)),2 ) + piPlusPiZero[4]);
  return sigma;
}

  double BaryonBaryon::nnToPiMinusPiZero(ATOOLS::Particle * A , ATOOLS::Particle * B)
  {
    return ppToPiPlusPiZero(A,B);
  }


double BaryonBaryon::pnToPiZeroPiZero(ATOOLS::Particle * A , ATOOLS::Particle * B) {
  std::vector<double> piZeroPiZero = {2.855*std::pow(10,6), 1.311*10, 2.961*1000, 5.572, 1.461*std::pow(10,6)};

  Vec4D p1 = A->Momentum(); Vec4D p2 = B->Momentum();
  Vec4D P = p1 + p2;    Poincare boost(P);
  Vec4D p1cm = p1;  boost.Boost(p1cm);
  // k = |p1 - p2| in the c.m. frame, i.e. twice the c.m. momentum of each
  // nucleon (arXiv:1504.07242, after Eq. (1)); the fits are in kappa = k/GeV.
  double p_rel = 2.*sqrt(p1cm[1]*p1cm[1]+p1cm[2]*p1cm[2]+p1cm[3]*p1cm[3]);
  // double s     = P[0]*P[0]-P[1]*P[1]-P[2]*P[2]-P[3]*P[3]; double sqrtS = std::sqrt(s);
  double sigma = 0;
  sigma = piZeroPiZero[0]*std::pow(p_rel,piZeroPiZero[1])/(std::pow( (piZeroPiZero[2] - std::exp(piZeroPiZero[3]*p_rel)),2 ) + piZeroPiZero[4]);

  return sigma;
}

double BaryonBaryon::pnToPiPlusPiMinus(ATOOLS::Particle * A , ATOOLS::Particle * B) {
  std::vector<double> piPlusPiMinusOne = {6.465*std::pow(10,6), 1.051*10, 1.979*1000, 5.363, 6.045*std::pow(10,5)};
  std::vector<double> piPlusPiMinusTwo = {2.549*std::pow(10,15), 1.657*10, 2.330*std::pow(10,7), 1.119*10, 2.868*std::pow(10,16)};

  Vec4D p1 = A->Momentum(); Vec4D p2 = B->Momentum();
  Vec4D P = p1 + p2;    Poincare boost(P);
  Vec4D p1cm = p1;  boost.Boost(p1cm);
  // k = |p1 - p2| in the c.m. frame, i.e. twice the c.m. momentum of each
  // nucleon (arXiv:1504.07242, after Eq. (1)); the fits are in kappa = k/GeV.
  double p_rel = 2.*sqrt(p1cm[1]*p1cm[1]+p1cm[2]*p1cm[2]+p1cm[3]*p1cm[3]);
  // double s     = P[0]*P[0]-P[1]*P[1]-P[2]*P[2]-P[3]*P[3]; double sqrtS = std::sqrt(s);
  double sigma = 0;
  sigma = piPlusPiMinusOne[0]*std::pow(p_rel,piPlusPiMinusOne[1])/(std::pow( (piPlusPiMinusOne[2] - std::exp(piPlusPiMinusOne[3]*p_rel)),2 ) + piPlusPiMinusOne[4]) + 
          piPlusPiMinusTwo[0]*std::pow(p_rel,piPlusPiMinusTwo[1])/(std::pow( (piPlusPiMinusTwo[2] - std::exp(piPlusPiMinusTwo[3]*p_rel)),2 ) + piPlusPiMinusTwo[4]);
  return sigma;
}

double BaryonBaryon::xstot(long int & A, long int & B,const double & plab) {
  return 0.;
}
  
double BaryonBaryon::xsel(long int & A, long int & B, const double & plab) {
  return 0.;
}

void BaryonBaryon::Tests() {
    size_t bins = 10000;
    double kmin = 0.05, kmax = 4.0;
    double bin_width = (kmax - kmin) / double(bins);

    map<string, Histogram*> histos;
    histos["pnToDGamma"]        = new Histogram(0, kmin, kmax, bins);

    histos["nnToDPiMinus"] = new Histogram(0, kmin, kmax, bins);
    histos["pnToDPiZero"]  = new Histogram(0, kmin, kmax, bins);
    histos["ppToDPiPlus"]  = new Histogram(0, kmin, kmax, bins);

    histos["ppToPiPlusPiZero"]  = new Histogram(0, kmin, kmax, bins);
    histos["pnToPiZeroPiZero"]  = new Histogram(0, kmin, kmax, bins);
    histos["pnToPiPlusPiMinus"] = new Histogram(0, kmin, kmax, bins);
    histos["nnToPiMinusPiZero"] = new Histogram(0, kmin, kmax, bins);

    Flavour flp(kf_p_plus), fln(kf_n);

    for (int i = 0; i < bins; i++) {
        double k   = kmin + (i + 0.5) * bin_width;  // bin centre, one per bin
        // the particles must be built with the flavours of each pair, the
        // flavours are not picked up if changed after construction
        // k = |p1 - p2| in the c.m. frame, each nucleon has |p| = k/2
        double q   = k/2.;
        Vec4D pp1(sqrt(flp.HadMass()*flp.HadMass() + q*q), 0., 0.,  q);
        Vec4D pp2(sqrt(flp.HadMass()*flp.HadMass() + q*q), 0., 0., -q);
        Vec4D pn1(sqrt(fln.HadMass()*fln.HadMass() + q*q), 0., 0.,  q);
        Vec4D pn2(sqrt(fln.HadMass()*fln.HadMass() + q*q), 0., 0., -q);
        Particle p1(-1, flp, pp1), p2(-1, flp, pp2);
        Particle n1(-1, fln, pn1), n2(-1, fln, pn2);

        histos["pnToDGamma"]       ->Insert(k, pnToDGamma(&p1, &n2));
        histos["pnToDPiZero"]      ->Insert(k, pnToDPiZero(&p1, &n2));
        histos["pnToPiZeroPiZero"] ->Insert(k, pnToPiZeroPiZero(&p1, &n2));
        histos["pnToPiPlusPiMinus"]->Insert(k, pnToPiPlusPiMinus(&p1, &n2));

        histos["nnToPiMinusPiZero"]->Insert(k, nnToPiMinusPiZero(&n1, &n2));
        histos["nnToDPiMinus"]     ->Insert(k, nnToDPiMinus(&n1, &n2));

        histos["ppToDPiPlus"]      ->Insert(k, ppToDPiPlus(&p1, &p2));
        histos["ppToPiPlusPiZero"] ->Insert(k, ppToPiPlusPiZero(&p1, &p2));
    }

    for (auto& [name, histo] : histos) {
        histo->Output("DeuteronFormationFromXSecs/" + name + ".dat");
        delete histo;
    }
    histos.clear();
}