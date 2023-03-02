#include "AMISIC++/Tools/Over_Estimator.H"
#include "AMISIC++/Perturbative/MI_Processes.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace AMISIC;
using namespace ATOOLS;

Over_Estimator::Over_Estimator() :
  m_muR_fac(1.), m_muF_fac(1.), m_pref(0.), m_bfac(1.), p_prefs(NULL)
{}

Over_Estimator::~Over_Estimator() { if (p_prefs) delete p_prefs; }

void Over_Estimator::Initialize(MI_Processes * procs) {
  /////////////////////////////////////////////////////////////////////////////////
  // The main trick here is to fix the maximal prefactor m_pref such that the
  // approximate cross section is larger than the actual one (given by the sum of
  // cross sections of all included parton-level processes), over all transverse
  // momenta.  We observe that the cross section dsigma/dp_T^2 falls quickly with
  // the transverse momentum squared p_T^2, suggesting a logarithmic binning in
  // p_T^2.  This is achieved in the FixMaximum method.
  /////////////////////////////////////////////////////////////////////////////////
  p_procs     = procs;
  m_s         = p_procs->S();
  m_pt02      = p_procs->PT02();
  m_ptmin2    = p_procs->PT2Min();
  p_alphaS    = p_procs->AlphaS();
  m_muR_fac   = (*mipars)("RenScale_Factor"); 
  m_muF_fac   = (*mipars)("FacScale_Factor");
  m_pt2bins   = size_t((*mipars)("nPT_bins"));
  m_xsnd      = p_procs->GetXSecs()->XSnd();
  for (size_t i=0;i<2;i++) {
    p_pdf[i]  = p_procs->PDF(i);
    m_xmin[i] = Max(1.e-6,p_pdf[i]->XMin());
  }
  FixMaximum();
}

void Over_Estimator::UpdateS() {
  ///////////////////////////////////////////////////////////////////////////////////
  // Updating to new variable centre-of-mass energy if necessary, and fixing a
  // suitable overestimating prefactor by interpolation in the look-up tables.
  // The non-diffractive cross section will have been calculated during the update
  // of the MI_Processes.
  ///////////////////////////////////////////////////////////////////////////////////
  m_s    = p_procs->S();
  m_xsnd = p_procs->GetXSecs()->XSnd();
  m_pref = (*p_prefs)(m_s);
}

void Over_Estimator::FixMaximum() {
  ///////////////////////////////////////////////////////////////////////////////////
  // Looping over all relevant s values from the Sudakov tables in the MI_Processes
  // to fill a look-up table of maxima.
  ///////////////////////////////////////////////////////////////////////////////////
  axis sbins = p_procs->GetSudakov()->GetSbins();
  p_prefs    = new OneDim_Table(sbins);
  for (size_t sbin=0;sbin<sbins.m_nbins;sbin++) {
    m_s = sbins.x(sbin);   
    double pt2step = log(m_s/(4.*m_ptmin2))/double(m_pt2bins);
    double maxpref = 0.;
    for (size_t bin=0;bin<m_pt2bins;bin++) {
      /////////////////////////////////////////////////////////////////////////////////
      // The actual value of p_T^2 for the bin, giving the maximal rapidity range
      // through s' = x_1 x_2 s = 4 p_T^2 cosh (Delta y), where Delta y is the (maximal)
      // rapidty distance of the tow outgoing partons, or twice the maximal rapidity of
      // each single parton.  yvol is the very crude and massively overestimated rapidity
      // volume both partons can occupy.
      /////////////////////////////////////////////////////////////////////////////////
      double pt2    = m_ptmin2*exp(pt2step*bin);
      double xt     = sqrt(4.*pt2/m_s);
      double yvol   = sqr(2.*log(1./xt*(1.+sqrt(1.-xt*xt))));
      double approx = ApproxME(pt2,xt);    
      double exact  = ExactME(pt2,xt);
      /////////////////////////////////////////////////////////////////////////////////
      // In both the approximate and the exact ME we factor out an appoximated, regularised
      // (with m_pt02) t-channel propagator and the rapidity volume to define a constant
      // prefactor that ensures that the approximation is always larger than the exact
      // calculation.
      /////////////////////////////////////////////////////////////////////////////////
      double test   = Max(approx,exact)*yvol*sqr(pt2+m_pt02/4.);
      if (test > maxpref) maxpref = test;
    }
    p_prefs->Fill(sbin,maxpref);
  }
}

double Over_Estimator::operator()(const double & pt2,const double & yvol) {
  /////////////////////////////////////////////////////////////////////////////////
  // Overestimator for differential cross section.
  // Have to divide out the volume of the rapidity integration which is also
  // not present in the true differential cross section calculated by the
  // MI_Processes.
  /////////////////////////////////////////////////////////////////////////////////
  return m_pref/(yvol * sqr(pt2+m_pt02/4.));
}


double Over_Estimator::ApproxME(const double & pt2,const double & xt) {
  /////////////////////////////////////////////////////////////////////////////////
  // Approximate differential cross section is given by
  // dsigma = pi/2 * alphaS^2(pT^2+pT0^2)/(pT^2+pT0^2)^2 *
  //          [sum_q q(xT,pT^2) + 9/4*g(xT,pT^2)]^2
  // which assumes
  // negligible s- and u-channel contributions and interference,
  // t-channel gluon-exchange dominance, and that the
  // product of PDFs is largest for both x being minimal, i.e. for x = xT = 4pT^2/s.
  /////////////////////////////////////////////////////////////////////////////////
  double reg_pt2 = pt2+m_pt02;
  double est     = M_PI/2.*sqr( (*p_alphaS)(Max(m_pt02,m_muR_fac * reg_pt2/4.)) / reg_pt2);
  for (size_t i=0;i<2;i++) {
    double pdfsum = 0.;
    if (xt>m_xmin[i]) {
      p_pdf[i]->Calculate(xt,Min(Max(m_muF_fac*pt2,p_pdf[i]->Q2Min()),p_pdf[i]->Q2Max()));
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

double Over_Estimator::ExactME(const double & pt2,const double & xt) {
  /////////////////////////////////////////////////////////////////////////////////
  // For the exact MEs we assume a "minimal" kinematics with smallest values for
  // s', t' and u' given by multiples of pT^2.  
  /////////////////////////////////////////////////////////////////////////////////
  if (xt>m_xmin[0] && xt>m_xmin[1]) return (*p_procs)(4.*pt2,-2.*pt2,-2.*pt2,xt,xt);
  return 0.;
}
  
double Over_Estimator::TrialPT2(const double & Q2) {
  /////////////////////////////////////////////////////////////////////////////////
  // Produce an overestimated q2 by solving for q2
  // random * exp[-int_{pt02}^{Q2} dpt2 prefb/(pt2+pt02/4)^2] =
  //          exp[-int_{q2}^{Q2}   dpt2 prefb/(pt2+pt02/4)^2]
  /////////////////////////////////////////////////////////////////////////////////
  double Q2tilde = Q2+m_pt02/4.;
  double prefb   = m_pref*m_bfac/m_xsnd;
  return prefb*Q2tilde/(prefb-Q2tilde*log(ran->Get())) - m_pt02/4.;
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
