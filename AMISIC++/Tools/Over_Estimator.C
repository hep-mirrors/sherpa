#include "AMISIC++/Tools/Over_Estimator.H"
#include "AMISIC++/Perturbative/MI_Processes.H"
#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace AMISIC;
using namespace ATOOLS;

//////////////////////////////////////////////////////////////////////////////
// This encodes a quick'n'dirty overestimator for the Sudakov form factor.
// The main trick is to construct an approximate cross section that is
// larger than the actual one (given by the sum of cross sections of all
// included parton-level processes), over all transverse momenta.
// As they are driven by t-channel exchange, usually as factors 1/|t|^2
// which is normalised in the IR by replacing 1/|t|^2 -> 1/(|t|+pt_0^2)^2
// we assume a form of N/(|t|+pt_0^2/4) and fix N.  We also observe that
// the cross section dsigma/dp_T^2 falls quickly with the transverse
// momentum squared p_T^2, suggesting a  logarithmic binning in p_T^2.
// This is achieved in the FixMaximum method.
/////////////////////////////////////////////////////////////////////////////

Over_Estimator::Over_Estimator() :
  m_muR_fac(1.), m_muF_fac(1.), m_pref(0.), m_npt2bins(1000),
  p_prefs(nullptr)
{}

Over_Estimator::~Over_Estimator() { if (p_prefs) delete p_prefs; }

void Over_Estimator::
Initialize(PDF::ISR_Handler * isr,MI_Processes * procs,axis * sbins) {
  ///////////////////////////////////////////////////////////////////////////
  // Inheriting all relevant inputs from the MI_Processes which we aim to
  // over-estimate.  
  ///////////////////////////////////////////////////////////////////////////
  m_s         = procs->S();
  m_pt02      = procs->PT02();
  m_ptmin2    = procs->PT2Min();
  p_alphaS    = procs->AlphaS();
  m_muR_fac   = (*mipars)("RenScale_Factor");
  m_muF_fac   = (*mipars)("FacScale_Factor");
  for (size_t i=0;i<2;i++) {
    p_pdf[i]  = isr->PDF(i);
    m_xmin[i] = Max(1.e-6,p_pdf[i]->XMin());
    m_xmax[i] = Max(1.-1.e-6,p_pdf[i]->XMax());
  }
  FixMaximum(procs,sbins);
  Output();
}

void Over_Estimator::UpdateS(const double & s,const double & pt02,const double & ptmin2) {
  ////////////////////////////////////////////////////////////////////////////
  // Updating to new variable centre-of-mass energy if necessary, and fixing a
  // suitable overestimating prefactor by interpolation in the look-up tables.
  ////////////////////////////////////////////////////////////////////////////
  m_pref   = (*p_prefs)(m_s=s);
  m_pt02   = pt02;
  m_ptmin2 = ptmin2;
}

void Over_Estimator::FixMaximum(MI_Processes * procs,axis * sbins) {
  ////////////////////////////////////////////////////////////////////////////
  // Looping over all relevant s values from the Sudakov tables in the 
  // MI_Processes to fill a look-up table of maxima.
  // Note: This does not include any effect of the matter overlap
  //       (its weight should be below unity)
  ////////////////////////////////////////////////////////////////////////////
  p_prefs    = new OneDim_Table(*sbins);
  for (size_t sbin=0;sbin<sbins->m_nbins;sbin++) {
    m_s      = sbins->x(sbin);
    m_pt02   = mipars->CalculatePT02(m_s);
    m_ptmin2 = mipars->CalculatePT02(m_s);
    double ratioN  = pow(m_s/m_ptmin2,1./double(m_npt2bins));
    double maxpref = 0.;
    for (size_t i=0;i<m_npt2bins;i++) {
      ///////////////////////////////////////////////////////////////////////
      // The actual value of p_T^2 for the bin, giving the maximal rapidity 
      // range through s' = x_1 x_2 s = 4 p_T^2 cosh (Delta y), where Delta y 
      // is the (maximal) rapidty distance of the two outgoing partons, or 
      // twice the maximal rapidity of each single parton.  yvol is the 
      // very crude and massively overestimated rapidity volume both partons
      // can occupy.
      ///////////////////////////////////////////////////////////////////////
      double pt2    = m_ptmin2 * pow(ratioN,i);
      double xt     = sqrt(4.*pt2/m_s);
      if (xt*xt > m_xmax[0]*m_xmax[1]) continue;
      double yvol   = sqr(2.*log(1./xt*(1.+sqrt(1.-xt*xt))));
      double approx = ApproxME(pt2,xt);
      double exact  = (*procs)(4.*pt2,-2.*pt2,-2.*pt2,xt,xt);
      ///////////////////////////////////////////////////////////////////////
      // In both the approximate and the exact ME we factor out an appoximated,
      // regularised (with m_pt02) t-channel propagator and the rapidity
      // volume to define a constant prefactor that ensures that the
      // approximation is always larger than the exact calculation.
      ///////////////////////////////////////////////////////////////////////
      double test   = procs->PDFnorm()*Max(approx,exact)*yvol*sqr(pt2+m_pt02/4.);
      if (test > maxpref) { maxpref = test; }
    }
    p_prefs->Fill(sbin,maxpref);
    if (sbins->m_nbins==1) m_pref = maxpref;
  }
}

double Over_Estimator::ApproxME(const double & pt2,const double & xt) {
  ///////////////////////////////////////////////////////////////////////////
  // Approximate differential cross section is given by
  // dsigma = pi/2 * alphaS^2[(pT^2+pT0^2)/4]/(pT^2+pT0^2)^2 *
  //          [sum_q q(xT,pT^2) + 9/4*g(xT,pT^2)]^2
  // which assumes
  // - negligible s- and u-channel contributions and interference,
  // - t-channel gluon-exchange dominance, and that
  // - the product of PDFs is largest for both x being minimal, i.e.
  //   for x = xT = 4pT^2/s.
  ///////////////////////////////////////////////////////////////////////////
  double reg_pt2 = pt2+m_pt02;
  double est     = ( M_PI/2. *
		     sqr( (*p_alphaS)(Max(m_pt02,m_muR_fac*reg_pt2/4.)) /
			  reg_pt2) );
  for (size_t i=0;i<2;i++) {
    double pdfsum = 0.;
    if (xt>m_xmin[i]) {
      p_pdf[i]->Calculate(xt,Min(Max(m_muF_fac*pt2,
				     p_pdf[i]->Q2Min()),p_pdf[i]->Q2Max()));
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

double Over_Estimator::operator()(const double & pt2,const double & yvol) {
  ///////////////////////////////////////////////////////////////////////////
  // Overestimator for differential cross section.
  // Have to divide out the volume of the rapidity integration which is also
  // not present in the true differential cross section calculated by the
  // MI_Processes.
  //////////////////////////////////////////////////////////////////////////
  double myvol = sqr(log(m_s/(4.*pt2)*sqr(1.+sqrt(1.-4.*pt2/m_s))));
  return m_pref/(yvol * sqr(pt2+m_pt02/4.));
}

double Over_Estimator::TrialPT2(const double & Q2) {
  ///////////////////////////////////////////////////////////////////////////
  // Produce an overestimated q2 by solving for q2
  // random = exp[-int_{q^2}^{Q2} dpt2 prefb/(pt2+pt02/4)^2]
  ///////////////////////////////////////////////////////////////////////////
  double Q2tilde = Q2+m_pt02/4.;
  return ( m_pref * Q2tilde/(m_pref-Q2tilde*log(Max(1.e-12,ran->Get()))) -
	   m_pt02/4. );
}

void Over_Estimator::Output() {
  msg_Info()<<"   "<<std::string(77,'-')<<"\n"
	    <<"   | Fixed (energy-dependent) maxima for "
	    <<"quick'n'dirty Sudakov evaluation"
	    <<"      |\n"
	    <<"   | Energy [GeV]       | Maximum"
	    <<std::string(46,' ')<<"|\n";
  for (size_t i=0;i<p_prefs->GetAxis().m_nbins;i++)
    msg_Info()<<"   | "
	      <<std::setprecision(6)<<std::setw(18)
	      <<sqrt(p_prefs->GetAxis().x(i))<<" | "
	      <<std::setprecision(8)<<std::setw(16)<<p_prefs->Value(i)
	      <<std::string(37,' ')<<"|\n";
  msg_Info()<<"   "<<std::string(77,'-')<<"\n\n";      
}

void Over_Estimator::Test(const double & Q2,const long int & n) {
  msg_Out()<<METHOD<<" for Q^2 = "<<Q2<<", s = "<<m_s<<".\n";
  Histogram histo(0,0.0,Q2,100);
  for (size_t i=0;i<n;i++) {
    double pt2 = Q2;
    size_t trial(0);
    while (pt2>m_pt02) {
      TrialPT2(pt2);
      if (trial++==0) histo.Insert(pt2);
    }
  }
  histo.Finalize();
  histo.Output("Over_PT2");
  msg_Out()<<METHOD<<": finished "<<n<<" dry runs.\n";
}
