#include "AMISIC++/Tools/Interaction_Probability.H"
#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace AMISIC;
using namespace ATOOLS;

/////////////////////////////////////////////////////////////////////////////
// The impact-parameter integral over the interaction probability should
// yield the non-diffractive cross section - we use this condition to
// fix the (energy-dependent) matter radii of the incident hadronic
// states entering the matter overlap.
//
// All equations in this file refer to
// Sjostrand-van der Zijl, PRD 36 (1987) 2019 or
// Corke-Sjostrand, JHEP 05 (2011) 009.
/////////////////////////////////////////////////////////////////////////////

Interaction_Probability::Interaction_Probability() :
  p_mo(NULL), p_procs(NULL), p_sbins(NULL), p_bbins(NULL),
  p_k(NULL), p_norm(NULL), p_diffxsec(NULL),
  m_pdfnorm(1.), m_bmax(1.e6), m_smin(1.e99), m_xs_hard(0.), m_xs_test(0.),
  m_dynamic(false), m_test(false), m_ana(true) {}

Interaction_Probability::~Interaction_Probability() {
  if (m_ana) FinishAnalysis();
  if (p_k)        delete p_k;
  if (p_norm)     delete p_norm;
  if (p_diffxsec) delete p_diffxsec;
}

void Interaction_Probability::
Initialize(Matter_Overlap * mo,MI_Processes * processes,axis * sbins) {
  if (m_ana) InitAnalysis();
  p_mo       = mo;
  p_procs    = processes;
  m_pdfnorm  = p_procs->PDFnorm();
  p_sbins    = sbins;
  p_k        = new OneDim_Table(*p_sbins);
  p_norm     = new OneDim_Table(*p_sbins);
  m_dynamic  = mo->IsDynamic();
  m_bmax     = mo->Bmax();
  p_bbins    = mo->GetBBins();
  p_diffxsec = new TwoDim_Table(*p_sbins,*p_bbins);
  FixKandSmin();
  OutputTables();
}

double Interaction_Probability::SelectB(const double & s) {
  /////////////////////////////////////////////////////////////////////////////
  // Selecting b according to (normalised)  d^2b P_int(b): 
  // the norm is the non-diffractive cross section.
  // Note: the m_pdfnorm accounts for possible terms alpha in the PDF of
  //       photons, and the XSndNorm takes care of some (external) rescaling.
  /////////////////////////////////////////////////////////////////////////////  
  double disc = 0., diff = 0., b = 0.;
  double xsnd = ( m_pdfnorm * p_procs->GetXSecs()->XSndNorm() *
		  p_procs->GetXSecs()->XSnd(s) );
  do {
    b    = ran->Get()*p_mo->Bmax();
    diff = (*this)(s,b);
    disc = 2.*M_PI*b*diff/xsnd;
  } while (disc<ran->Get());
  return b;
}

void Interaction_Probability::FixKandSmin() {
  /////////////////////////////////////////////////////////////////////////////
  // We iteratively solve for a prefactor k to rescale the radii in the
  // matter overlap such that
  // int d^2b P_int(b) = xs_ND, see above for the norms.
  /////////////////////////////////////////////////////////////////////////////
  msg_Info()<<"   "<<string(77,'-')<<"\n"
	    <<"   | Initializing the Interaction Probability: fixing K's"
	    <<string(22,' ')<<"|\n";
  PInt_Dyn_Integrand integrand(p_diffxsec);
  Gauss_Integrator   gauss(&integrand);
  for (size_t sbin=0;sbin<p_sbins->m_nbins;sbin++) {
    double s     = p_sbins->x(sbin), k=1., xs_test = 0.;
    double xs_nd = (m_pdfnorm * p_procs->GetXSecs()->XSndNorm() *
		    p_procs->GetXSecs()->XSnd(s));
    if (s<=4. || p_procs->GetXSecs()->XSratio(s)<=1.) continue;
    integrand.SetS(s);
    do {
      if (p_procs->GetXSecs()->XSratio(s)<=0.5) { k = 1.; break; }
      p_mo->SetKRadius(k);
      InitializeTable(sbin);
      xs_test = gauss.Integrate(0.,p_mo->Bmax(),1.e-3);
      if (dabs(xs_test/xs_nd)<1.e-3) { k = 0.; break; }
      /////////////////////////////////////////////////////////////////////////
      // we rescale the k with the ratio of implied and exact ND cross section,
      // until we are within 2% of each other.
      /////////////////////////////////////////////////////////////////////////
      k *= Min(5., Max(0.2, sqrt(xs_nd/xs_test)));
    } while (dabs(1.-xs_test/xs_nd)>0.02);
    p_k->Fill(sbin,k);
    msg_Info()<<"   | "
	      <<"E = "<<std::setw(8)<<std::setprecision(6)<<sqrt(s)<<": "
	      <<"xs_ND = "<<std::setw(12)<<std::setprecision(6)<<xs_nd<<", "
	      <<"ratio = "<<std::setw(12)<<std::setprecision(6)
	      <<p_procs->GetXSecs()->XSratio(s)<<" -> "
	      <<"k = "<<std::setw(8)<<std::setprecision(6)<<k<<"  |\n";
    if (s<m_smin) m_smin = s;
  }
  msg_Info()<<"   "<<string(77,'-')<<"\n";
}

void Interaction_Probability::InitializeTables() {
  for (size_t sbin=0;sbin<p_sbins->m_nbins;sbin++) InitializeTable(sbin);
}

void Interaction_Probability::
InitializeTable(const size_t & sbin,const bool & out) {
  /////////////////////////////////////////////////////////////////////////
  // Filling, in impact parameter space b the exponent for the interaction
  // probability, i.e.
  // int [ dx_1 dx_2 f(x_1) f(x_2) dpT^2 O(x_1,x_2, ...) dsigma/dpT^2 ]
  // for a single value of s (given by the s-bin)
  // Ultimately this will result in a 2-dim look-up table in s and b.
  /////////////////////////////////////////////////////////////////////////
  MI_Integrator * integrator = p_procs->GetIntegrator();
  double s  = p_sbins->x(sbin), prev = 0., xsfix = 0.;
  p_procs->UpdateS(s);
  if (!p_mo->IsDynamic()) xsfix = (*integrator)(s,nullptr,0.);
  if (out)
    msg_Info()<<"   | "<<METHOD<<"(E = "
	      <<std::setw(8)<<std::setprecision(6)<<sqrt(s)<<", " 
	      <<"dyn = "<<p_mo->IsDynamic()<<"):          |\n";
  for (size_t bbin=0;bbin<p_bbins->m_nbins;bbin++) {
    double b  = p_bbins->x(bbin);
    double xs = m_pdfnorm * ( p_mo->IsDynamic() ?
			      (*integrator)(s,p_mo,b) :
			      xsfix * (*p_mo)(b) );
    p_diffxsec->Fill(sbin,bbin,xs);
    if (prev!=0. && xs/prev < 1.e-6) {
      for (size_t i=bbin+1;i<p_bbins->m_nbins;i++)
	p_diffxsec->Fill(sbin,i,0.);
      break;
    }
    prev = xs;
  }
}

bool Interaction_Probability::CheckTables() {
  XS_Integrand     integrand(p_diffxsec,0.);
  Gauss_Integrator gauss(&integrand);
  for (size_t sbin=0;sbin<p_sbins->m_nbins;sbin++) {
    double s = p_sbins->x(sbin);
    if (s<=1. || p_procs->GetXSecs()->XSratio(s)<=1.) continue;
    integrand.SetS(s);
    m_xs_test = gauss.Integrate(0.,p_bbins->m_xmax,1.e-6);
    m_xs_hard = p_procs->GetXSecs()->XShard(s);
    if (dabs(1.-m_xs_test/m_xs_hard)>10.e-2) {
      msg_Error()<<"Warning in "<<METHOD<<"(s = "<<s<<"):\n"
		 <<"   Cross sections differ by more than 10 % between "
		 <<"hard xs = "<<m_xs_hard<<" "
		 <<"and b-integrated cross section = "<<m_xs_test<<";\n"
		 <<"   Will continue and hope for the best.\n";
      return false;
    }
  }
  return true;
}

void Interaction_Probability::OutputTables() {
  msg_Info()<<"   "<<std::string(77,'-')<<"\n"
	    <<"   | "<<METHOD<<" look-up tables and values:          |\n"
	    <<"   | "<<std::setw(15)<<"E_{c.m.} [GeV]"<<" | "
	    <<std::setw(6)<<"b [fm]"<<" | "
	    <<std::setw(14)<<"xs_hard/xs_ND"<<" | "
	    <<std::setw(12)<<"k"<<" | "
	    <<std::setw(14)<<"b*dP_int(b)"<<" |\n";
  for (size_t i=0;i<p_sbins->m_nbins;i++) {
    double s       = p_sbins->x(i);
    double xsratio = p_procs->XSratio(s);
    msg_Info()<<"   | "<<std::setprecision(6)<<std::setw(15)<<sqrt(s)<<" | "
	      <<std::string(6,' ')<<" | "
	      <<std::setprecision(6)<<std::setw(14)<<xsratio<<" | "
	      <<std::setprecision(6)<<std::setw(12)<<p_k->Value(i)<<" | "
	      <<std::string(14,' ')<<" |\n";
    for (size_t j=0;j<20;j++) {
      double b = p_mo->Bmax()*double(j)/40.;
      msg_Info()<<"   | "<<std::string(15,' ')<<" | "
		<<std::setprecision(3)<<std::setw(6)
		<<(b*rpa->hBar()*rpa->c()*1.e12)<<" | "
		<<std::setprecision(6)<<std::setw(14)
		<<((*p_diffxsec)(s,b)/xsratio)<<" | "
		<<std::string(12,' ')<<" | "
		<<std::setprecision(6)<<std::setw(14)
		<<(b*(*this)(s,b))<<" |\n";
    }
  }
  msg_Info()<<"   "<<std::string(77,'-')<<"\n"
	    <<"   | hard cross section with and without overlap:"
	    <<std::string(30,' ')<<"|\n"
	    <<"   | xs_hard = "<<std::setprecision(5)<<std::setw(6)
	    <<(m_xs_hard*rpa->Picobarn()/1.e9)<<" mb "
	    <<" (without overlap) vs. "<<std::setprecision(5)<<std::setw(6)
	    <<(m_xs_test*rpa->Picobarn()/1.e9)<<" (with overlap) mb.      |\n"
	    <<"   "<<std::string(77,'-')<<"\n\n\n";
}

void Interaction_Probability::InitAnalysis() {
  double bmax = 20.*rpa->hBar()*rpa->c()*1.e12;
  m_histos[std::string("b_times_P_in")] = new Histogram(0, 0., bmax, 200);
}

void Interaction_Probability::Analyse() {
  double s = sqr(rpa->gen.Ecms());
  for (size_t i=0;i<100;i++) {
    double b = 10./double(100)*(double(i)+0.5)*rpa->hBar()*rpa->c()*1.e12;
    m_histos[std::string("b_times_P_in")]->Insert(b,b*(*this)(s,b));
  }
}

void Interaction_Probability::FinishAnalysis() {
  Histogram * histo;
  std::string name;
  for (std::map<std::string,Histogram *>::iterator
	 hit=m_histos.begin();hit!=m_histos.end();hit++) {
    histo = hit->second;
    name  = std::string("MPI_Analysis/")+hit->first+std::string(".dat");
    histo->Finalize();
    histo->Output(name);
    delete histo;
  }
  m_histos.clear();
}
