#include "AMISIC++/Tools/Interaction_Probability.H"
#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include <algorithm>

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
  p_k(NULL), p_norm(NULL), p_diffxsec(NULL), p_inverted(NULL), 
  m_pdfnorm(1.), m_bmax(1.e6), m_smin(1.e99), m_xs_hard(0.), m_xs_test(0.),
  m_dynamic(false), m_test(false), m_ana(false), m_n_variations(1) {}

Interaction_Probability::~Interaction_Probability() {
  if (m_ana) FinishAnalysis();
  if (p_k)        delete p_k;
  if (p_norm)     delete p_norm;
  if (p_diffxsec) delete p_diffxsec;
  if (p_inverted) delete p_inverted;
  for (auto p : p_k_variations) {
    if (p) delete p;
  }
}

void Interaction_Probability::
Initialize(Matter_Overlap * mo,MI_Processes * processes,axis * sbins) {
  if (m_ana) InitAnalysis();
  p_mo         = mo;
  p_procs      = processes;
  m_pdfnorm    = p_procs->PDFnorm();
  p_sbins      = sbins;
  p_k          = new OneDim_Table(*p_sbins);
  p_norm       = new OneDim_Table(*p_sbins);
  m_dynamic    = mo->IsDynamic();
  m_bmax       = mo->Bmax();
  p_bbins      = mo->GetBBins();
  p_diffxsec   = new TwoDim_Table(*p_sbins,*p_bbins);
  FixKandSmin();
  FillIntegrated();
  OutputTables();
}

double Interaction_Probability::SelectB(const double & s) {
  /////////////////////////////////////////////////////////////////////////////
  // Selecting b according to the (normalised)  d^2b P_int(b),
  // where the norm is the non-diffractive cross section.
  // We have encoded this in a look-up table obtained from inverting the
  // probability distribution.
  /////////////////////////////////////////////////////////////////////////////
  return (*p_inverted)(s,ran->Get());
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

  // Initialize K-factor tables and differential XS tables for variations
  p_k_variations.resize(m_n_variations);
  for (size_t ivar=0; ivar<m_n_variations; ++ivar) {
    p_k_variations[ivar] = new OneDim_Table(*p_sbins);
  }

  PInt_Dyn_Integrand integrand(p_diffxsec);
  Gauss_Integrator   gauss(&integrand);
  for (size_t sbin=0;sbin<p_sbins->m_nbins;sbin++) {
    double s     = p_sbins->x(sbin), k=1., xs_test = 0.;
    double xs_nd = (m_pdfnorm * p_procs->GetXSecs()->XSndNorm() * 
                    p_procs->GetXSecs()->XSnd(s));
    if (s<=4. || p_procs->GetXSecs()->XSratio(s)<=1.) continue;
    integrand.SetS(s);
    p_mo->SetMatterFormVariationIndex(0);
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

    for (size_t ivar=0; ivar<m_n_variations; ++ivar) {
      double xs_nd_var = (m_pdfnorm * m_sigma_nd_variations[ivar] *
                          p_procs->GetXSecs()->XSnd(s));
      double k_var = 1.0;
      double xs_test_var = 0.0;
      double xsratio_var = p_procs->GetXSecs()->XSratio(s, ivar);

      // Compute the variation's own bmax and bmin
      double bmax_var = p_mo->ComputeBmaxForVariation(ivar);
      double bmin_var = p_mo->ComputeBminForVariation(ivar);

      // Create a temporary b-grid for this variation
      size_t nbbins = size_t((*mipars)["nB_bins"]);
      axis temp_bbins(nbbins, bmin_var, bmax_var, axis_mode::log);
      TwoDim_Table * temp_diffxsec = new TwoDim_Table(*p_sbins, temp_bbins);

      PInt_Dyn_Integrand integrand_var(temp_diffxsec);
      Gauss_Integrator gauss_var(&integrand_var);
      integrand_var.SetS(s);

      p_mo->SetMatterFormVariationIndex(ivar);

      // Solve for K-factor for this variation
      do {
        if (xsratio_var<=0.5) { k_var = 1.; break; }
        p_mo->SetKRadius(k_var);
        FillTableOnAxis(sbin, temp_diffxsec, &temp_bbins, ivar);
        xs_test_var = gauss_var.Integrate(0., bmax_var, 1.e-3);
        if (dabs(xs_test_var/xs_nd_var)<1.e-3) { k_var = 0.; break; }
        k_var *= Min(5., Max(0.2, sqrt(xs_nd_var/xs_test_var)));
      } while (dabs(1.-xs_test_var/xs_nd_var)>0.02);

      p_k_variations[ivar]->Fill(sbin, k_var);

      msg_Info()<<"   | Variation "<<ivar<<" : "
                <<"xs_ND = "<<std::setw(12)<<std::setprecision(6)<<xs_nd_var<<", "
                <<"ratio = "<<std::setw(12)<<std::setprecision(6)<<xsratio_var<<" -> "
                <<"k = "<<std::setw(8)<<std::setprecision(6)<<k_var<<"  |\n";
    }

    // Restore nominal K-factor in matter overlap
    p_mo->SetMatterFormVariationIndex(0);
    p_mo->SetKRadius(k);
    InitializeTable(sbin);

    if (s<m_smin) m_smin = s;
  }
  msg_Info()<<"   "<<string(77,'-')<<"\n";
}

void Interaction_Probability::FillIntegrated() {
  TwoDim_Table integrated(*p_sbins,*p_bbins);
  for (size_t sbin=0;sbin<p_sbins->m_nbins;sbin++) {
    double s   = p_sbins->x(sbin);
    double b1  = 0.,            b2 = 0.;
    double xs1 = (*this)(s,b1), xs2 = 0., norm = 0., integ = 0.;
    for (size_t bbin=1;bbin<p_bbins->m_nbins;bbin++) {
      b2    = p_bbins->x(bbin);
      xs2   = (*this)(s,b2);
      norm += (b2-b1)*(b2*xs2+b1*xs1)/2.;
      b1    = b2;
      xs1   = xs2;
    }
    integrated.Fill(sbin,0,0.);
    b1 = 0.; xs1 = (*this)(s,b1), b2 = 0., xs2 = 0.;
    for (size_t bbin=1;bbin<p_bbins->m_nbins;bbin++) {
      b2     = p_bbins->x(bbin);
      xs2    = (*this)(s,b2);
      integ += (b2-b1)*(b2*xs2+b1*xs1)/2./norm;
      integrated.Fill(sbin,bbin,integ);
      b1     = b2;
      xs1    = xs2;
    }
  }
  p_inverted = integrated.Invert(1,1000);
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

void Interaction_Probability::
FillTableOnAxis(const size_t & sbin, TwoDim_Table* diffxsec,
                axis* bbins, size_t ivar) {
  MI_Integrator * integrator = p_procs->GetIntegrator();
  double s  = p_sbins->x(sbin), prev = 0., xsfix = 0.;
  p_procs->UpdateS(s);

  double ptmin2_var = mipars->CalculatePTmin2(s, ivar);
  double pt02_var   = mipars->CalculatePT02(s, ivar);
  integrator->SetPT2min(ptmin2_var);
  p_procs->SetPT02(pt02_var);
  
  if (!p_mo->IsDynamic()) xsfix = (*integrator)(s,nullptr,0.);
  
  for (size_t bbin=0;bbin<bbins->m_nbins;bbin++) {
    double b  = bbins->x(bbin);
    double xs = m_pdfnorm * ( p_mo->IsDynamic() ?
                              (*integrator)(s,p_mo,b) :
                              xsfix * (*p_mo)(b) );
    diffxsec->Fill(sbin,bbin,xs);
    if (prev!=0. && xs/prev < 1.e-6) {
      for (size_t i=bbin+1;i<bbins->m_nbins;i++)
        diffxsec->Fill(sbin,i,0.);
      break;
    }
    prev = xs;
  }

  // Restore nominal values
  integrator->SetPT2min(mipars->CalculatePTmin2(s));
  p_procs->SetPT02(mipars->CalculatePT02(s));
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
            <<"   "<<std::string(77,'-')<<"\n\n";

  if (m_n_variations == 1) return;
  for (size_t ivar=0; ivar<m_n_variations; ++ivar) {
    double bmax_var = p_mo->ComputeBmaxForVariation(ivar);
    double bmin_var = p_mo->ComputeBminForVariation(ivar);

    size_t nbbins = size_t((*mipars)["nB_bins"]);
    axis temp_bbins(nbbins, bmin_var, bmax_var, axis_mode::log);
    TwoDim_Table * temp_diffxsec = new TwoDim_Table(*p_sbins, temp_bbins);

    // Compute tables with variation K-factor
    for (size_t sbin=0;sbin<p_sbins->m_nbins;sbin++) {
      double k_var = p_k_variations[ivar]->Value(sbin);
      p_mo->SetMatterFormVariationIndex(ivar);
      p_mo->SetKRadius(k_var);
      FillTableOnAxis(sbin, temp_diffxsec, &temp_bbins, ivar);
    }

    msg_Info()<<"   "<<std::string(77,'-')<<"\n"
              <<"   | "<<METHOD<<" Variation "<<ivar<<std::setw(25)<<": "<<"|\n"
              <<"   | "<<std::setw(15)<<"E_{c.m.} [GeV]"<<" | "
              <<std::setw(6)<<"b [fm]"<<" | "
              <<std::setw(14)<<"xs_hard/xs_ND"<<" | "
              <<std::setw(12)<<"k"<<" | "
              <<std::setw(14)<<"b*dP_int(b)"<<" |\n";

    for (size_t j=0;j<p_sbins->m_nbins;j++) {
      double s           = p_sbins->x(j);
      double xsratio_var = p_procs->GetXSecs()->XSratio(s, ivar);

      msg_Info()<<"   | "<<std::setprecision(6)<<std::setw(15)<<sqrt(s)<<" | "
                <<std::string(6,' ')<<" | "
                <<std::setprecision(6)<<std::setw(14)<<xsratio_var<<" | "
                <<std::setprecision(6)<<std::setw(12)<<p_k_variations[ivar]->Value(j)<<" | "
                <<std::string(14,' ')<<" |\n";
      for (size_t k=0;k<20;k++) {
        double b = bmax_var*double(k)/40.;
        msg_Info()<<"   | "<<std::string(15,' ')<<" | "
                  <<std::setprecision(3)<<std::setw(6)
                  <<(b*rpa->hBar()*rpa->c()*1.e12)<<" | "
                  <<std::setprecision(6)<<std::setw(14)
                  <<((*temp_diffxsec)(s,b)/xsratio_var)<<" | "
                  <<std::string(12,' ')<<" | "
                  <<std::setprecision(6)<<std::setw(14)
                  <<(b*(1.-exp(-(*temp_diffxsec)(s,b))))<<" |\n";
      }
    }
    delete temp_diffxsec;

    // Restore nominal K-factor
    p_mo->SetMatterFormVariationIndex(0);
    if (p_sbins->m_nbins > 0) {
      double k_nom = p_k->Value(0);
      p_mo->SetKRadius(k_nom);
    }
  }
  msg_Info()<<"   "<<std::string(77,'-')<<"\n\n";
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
