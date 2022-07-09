#include "AMISIC++/Tools/Over_Estimator.H"
#include "AMISIC++/Perturbative/MI_Processes.H"
#include "AMISIC++/Tools/MI_Parameters.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace AMISIC;
using namespace ATOOLS;

Over_Estimator::Over_Estimator() :
  m_muR_fac(1.), m_muF_fac(1.),
  m_nbins(100), m_pref(1.), m_bfac(1.)
{}

Over_Estimator::~Over_Estimator() {}

void Over_Estimator::Initialize(MI_Processes * procs) {
  // The main trick here is to fix the maximal prefactor m_pref such that the
  // approximate cross section is larger than the actual one (given by the sum of
  // cross sections of all included parton-level processes), over all transverse
  // momenta.  We observe that the cross section dsigma/dp_T^2 falls quickly with
  // the transverse momentum squared p_T^2, suggesting a logarithmic binning in
  // p_T^2.  This is achieved in the FixMaximum method.
  m_s      = sqr(procs->Ecms());
  m_pt02   = sqr(procs->PT0());
  m_ptmin2 = sqr(procs->PTmin());
  p_alphaS = procs->AlphaS();
  for (size_t i=0;i<2;i++) {
    p_pdf[i]  = procs->PDF(i);
    m_xmin[i] = Max(1.e-6,p_pdf[i]->XMin());
  }
  FixMaximum(procs);
}

void Over_Estimator::FixMaximum(MI_Processes * procs) {
  m_muR_fac      = (*mipars)("RenScale_Factor"); 
  m_muF_fac      = (*mipars)("FacScale_Factor"); 
  double pt2step = log(m_s/(4.*m_ptmin2))/double(m_nbins);
  double int1(0.),int2(0.);
  for (size_t bin=0;bin<m_nbins;bin++) {
    // The actual value of p_T^2 for the bin, giving the maximal rapidity range
    // through s' = x_1 x_2 s = 4 p_T^2 cosh (Delta y), where Delta y is the (maximal)
    // rapidty distance of the tow outgoing partons, or twice the maximal rapidity of
    // each single parton.  yvol is the very crude and massively overestimated rapidty
    // volume both partons can occupy.
    double pt2    = m_ptmin2*exp(pt2step*bin);
    m_xt          = sqrt(4.*pt2/m_s);
    double yvol   = sqr(2.*log(1./m_xt*(1.+sqrt(1.-m_xt*m_xt))));
    double approx = ApproxME(pt2);    
    double exact  = ExactME(procs,pt2);
    // In both the approximate and the exact ME we factor out an appoximated, regularised
    // (with m_pt02) t-channel propagator and the rapidity volume to define a constant
    // prefactor that ensures that the approximation is always larger than the exact
    // calculation.
    double test   = Max(approx,exact)*yvol*sqr(pt2+m_pt02/4.);
    if (bin>=1) {
      int1 += approx*yvol*m_ptmin2*(exp(pt2step*bin)-exp(pt2step*(bin-1.)));
      int2 += exact *yvol*m_ptmin2*(exp(pt2step*bin)-exp(pt2step*(bin-1.)));
    }
    if (double(bin/10)==bin/10.)
      msg_Tracking()<<METHOD<<"(pt = "<<sqrt(pt2)<<"): "<<approx<<" vs. "<<exact<<".\n";
    
    if (test>m_pref) m_pref = test;
  }
}


double Over_Estimator::ApproxME(const double & pt2) {
  // Approximate differential cross section is given by
  // dsigma = pi/2 * alphaS^2(pT^2+pT0^2)/(pT^2+pT0^2)^2 *
  //          [sum_q q(xT,pT^2) + 9/4*g(xT,pT^2)]^2
  // which assume that s- and u-channel contributions and interference can be neglected,
  // that the t-channel exchange of gluons in elastic scattering dominates, and that the
  // product of PDFs is largest for both x being minimal, i.e. for x = xT = 4pT^2/s.
  double scale = pt2+m_pt02;
  double est   = M_PI/2.*sqr((*p_alphaS)(Max(m_pt02,m_muR_fac * scale/4.))) / sqr(scale);
  for (size_t i=0;i<2;i++) {
    double pdfsum = 0.;
    double Q2     = m_muF_fac*Max(pt2,p_pdf[i]->Q2Min());
    if (m_xt>m_xmin[i]) {
      p_pdf[i]->Calculate(m_xt,Q2);
      for (Flavour_Set::const_iterator fl=p_pdf[i]->Partons().begin();
	   fl!=p_pdf[i]->Partons().end();fl++) {
	// only allow u, d, s, c, b quarks and gluons
	if (fl->Kfcode()>=6 && fl->Kfcode()!=21) continue;
	pdfsum += Max(0.,((*fl).IsGluon()?9./4.:1.)*p_pdf[i]->GetXPDF(*fl));
      }
    }
    est *= pdfsum;
  }
  return est;
}

double Over_Estimator::ExactME(MI_Processes * procs,const double & pt2) {
  // For the exact MEs we assume a "minimal" kinematics with smallest values for
  // s', t' and u' given by multiples of pT^2.  
  double shat  = 4.*pt2, that=-2.*pt2, uhat=-2.*pt2;
  double x1    = m_xt, x2 = m_xt;
  double scale = pt2;
  if (x1<m_xmin[0] || x2<m_xmin[1]) return 0.;
  procs->CalcPDFs(x1,x2,m_muF_fac * scale);
  return (*procs)(shat,that,uhat);
}
  
double Over_Estimator::operator()(const double & pt2) {
  return m_pref/sqr(pt2+m_pt02/4.);
}

double Over_Estimator::TrialPT2(const double & Q2) {
  // Produce an overestimated q2 by solving for q2
  // random * exp[-int_{pt02}^{Q2} dpt2 prefb/(pt2+pt02/4)^2] =
  //          exp[-int_{q2}^{Q2}   dpt2 prefb/(pt2+pt02/4)^2]
  double Q2tilde = Q2+m_pt02/4.;
  double prefb   = m_pref*m_bfac/m_xsnd;
  double pt2     = prefb*Q2tilde/(prefb-Q2tilde*log(ran->Get())) - m_pt02/4.;
  return pt2;
}

void Over_Estimator::Test(const double & Q2,const long int & n) {
  msg_Out()<<METHOD<<" for Q^2 = "<<Q2<<", s = "<<m_s<<".\n";
  Histogram histo(0,0.0,Q2,100);
  for (size_t i=0;i<n;i++) {
    double pt2 = Q2;
    size_t trial(0);
    while (pt2>m_pt02) {
      pt2 = TrialPT2(pt2);
      if (trial++==0) histo.Insert(pt2);
    }
  }
  histo.Finalize();
  histo.Output("Over_PT2");
  msg_Out()<<METHOD<<": finished "<<n<<" dry runs.\n";
}
