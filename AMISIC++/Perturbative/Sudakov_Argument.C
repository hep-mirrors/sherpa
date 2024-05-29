#include "AMISIC++/Perturbative/Sudakov_Argument.H"
#include "AMISIC++/Perturbative/MI_Processes.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace AMISIC;
using namespace ATOOLS;

Sudakov_Argument::
Sudakov_Argument(MI_Processes * procs) :
  p_processes(procs),
  p_function1(NULL), p_integral1(NULL),
  p_function2(NULL), p_integral2(NULL),
  p_function3(NULL), p_integral3(NULL),
  m_bJac(1.),
  m_test(false)
{}

Sudakov_Argument::~Sudakov_Argument() {
  if (p_function1) delete p_function1;
  if (p_integral1) delete p_integral1;
  if (p_function2) delete p_function2;
  if (p_integral2) delete p_integral2;
  if (p_function3) delete p_function3;
  if (p_integral3) delete p_integral3;
  if (p_pt2bins)   delete p_pt2bins;
}

void Sudakov_Argument::Initialize(axis * sbins,axis * bbins)
{
  ///////////////////////////////////////////////////////////////////////////
  // These parameters are for the hard part of the Sudakov form factor.
  // It is given by int_{pt^2}^{s/4} dpt^2 dsigma/dpt^2, and we tabulate the
  // integral in nbins between s/4 and pt_min^2, resulting in a stepsize of
  // pt^2_step.  In each of the bins we MC integrate over the rapidities of the
  // two outgoing particles with MC_points points.
  ///////////////////////////////////////////////////////////////////////////
  p_sbins      = sbins;
  p_bbins      = bbins;
  m_variable_s = p_sbins->m_nbins>1;
  m_variable_b = (p_bbins!=NULL);
  p_pt2bins    = new axis(size_t((*mipars)("nPT_bins")),
			  p_processes->PT2Min(),p_processes->S(),
			  axis_mode::log);
  InitTables();
  FillTables();
  if (m_test) { OutputTables(); THROW(normal_exit,"testing complete"); }
}

double Sudakov_Argument::TrialEvent(const double & s) {
  return p_processes->GetIntegrator()->TrialEvent(s);
}

void Sudakov_Argument::InitTables() {
  if (m_variable_b && m_variable_s) {
    p_integral3 = new ThreeDim_Table(*p_sbins,*p_bbins,*p_pt2bins);
    p_function3 = new ThreeDim_Table(*p_sbins,*p_bbins,*p_pt2bins);
  }
  else if (m_variable_s) {
    p_integral2 = new TwoDim_Table(*p_sbins,*p_pt2bins);
    p_function2 = new TwoDim_Table(*p_sbins,*p_pt2bins);
  }
  else if (m_variable_b) {
    p_integral2 = new TwoDim_Table(*p_bbins,*p_pt2bins);
    p_function2 = new TwoDim_Table(*p_bbins,*p_pt2bins);
  }
  else {
    p_integral1 = new OneDim_Table(*p_pt2bins);
    p_function1 = new OneDim_Table(*p_pt2bins);
  }
}

void Sudakov_Argument::FillTables() {
  //////////////////////////////////////////////////////////////////////////
  // Iterating over the bins in s (for hadron-hadron collisions typically
  // only one, while more for collisions involving e.g. EPA photons) to
  // fill the table of Sudakov arguments, i.e. int dpt^2 dsigma/dpt^2.
  // We do this with exact MEs, as we use these table to fix the b-dependence,
  // for example for MinBias events - so the approximate hit-or-miss form
  // will not be the correct way of doing this.
  // These table will be used in the Impact_Parameter class.
  //////////////////////////////////////////////////////////////////////////
  m_sbin = m_bbin = 0;
  size_t sbinmax  = (m_variable_s? p_sbins->m_nbins : 1);
  size_t bbinmax  = (m_variable_b? p_bbins->m_nbins : 1);
  msg_Info()<<"   "<<std::string(77,'-')<<"\n"
	    <<"   | "<<METHOD<<": "<<std::string(44,' ')<<"|\n"
	    <<"   | Filling ";
  size_t nspaces=0;
  msg_Info()<<std::setw(4)<<sbinmax<<" bins in s, "
	    <<std::setw(4)<<bbinmax<<" bins in b and "
	    <<std::setw(4)<<p_pt2bins->m_nbins<<" bins in pt2.             |\n"
	    <<"   | This might take a while."<<std::string(50,' ' )<<"|\n";
  double s = 0, b = 0, integral = 0.;
  for (m_sbin=0;m_sbin<sbinmax;m_sbin++) {
    double s = p_sbins->x(m_sbin), xshd = 0.;
    ////////////////////////////////////////////////////////////////////////
    // The (re-)scaled non-perturbative cross section, in GeV^-2, acts as
    // normalization in the Sudakov factor.
    ////////////////////////////////////////////////////////////////////////
    (*p_processes->GetXSecs())(s);
    p_processes->UpdateS(s);
    m_norm = (
                  p_processes->GetXSecs()->XSndNorm()
                          * p_processes->GetXSecs()->XSnd());msg_Info()<<"   | E = "
	      <<std::setprecision(1)<<std::setw(8)<<sqrt(s)<<" GeV "
	      <<"[bin = "<<std::setw(4)<<(m_sbin+1)<<"/"<<sbinmax<<"]: norm = "
	      <<std::setprecision(6)<<std::setw(10)
	      <<(m_norm*rpa->Picobarn()/1.e9)<<" mb"
	      <<std::string(21,' ')<<"|\n"<<"   "<<std::string(77,'-')<<"\n";
    if (m_variable_b)
      msg_Info()<<"   |     b [fm] |     Int d^2b xs_hard(b)/xs_nd |"
		<<std::string(19,' ')<<"cumulative |\n";
    ////////////////////////////////////////////////////////////////////////
    // Jacobean for b-integration, equals unity if no impact parameter
    // dependence of hard cross section, i.e. matter overlap is not dynamic.
    ////////////////////////////////////////////////////////////////////////
    m_bJac = 1.;
    ////////////////////////////////////////////////////////////////////////
    // Filling bins in b, including the Jacobean (2 pi b db) and setting the
    // impact parameter in the MI_Processes
    ////////////////////////////////////////////////////////////////////////
    for (m_bbin=0;m_bbin<(m_variable_b ? bbinmax-1 : bbinmax);m_bbin++) {
      if (m_variable_b) {
	double b = (p_bbins->x(m_bbin)+p_bbins->x(m_bbin+1))/2.;
	m_bJac   = 2.*M_PI*b * dabs(p_bbins->x(m_bbin+1)-p_bbins->x(m_bbin));
	p_processes->SetB(b);
	msg_Info()<<"   | "<<std::setprecision(4)<<std::setw(10)<<b<<" | ";
      }
      ///////////////////////////////////////////////////////////////////////
      // Filling bins in pT
      ///////////////////////////////////////////////////////////////////////
      xshd += FillPT2Values(s);
      OutputBin(integral);
    }
    msg_Out()<<"   | integrated total hard cross section = "
	     <<std::setprecision(6)<<std::setw(10)
	     <<(xshd*m_norm*rpa->Picobarn()/1.e9)<<" mb. "
	     <<std::string(21,' ')<<"|\n"
	     <<"   "<<std::string(77,'-')<<"\n";
  }
}

double Sudakov_Argument::FillPT2Values(const double & s) {
  //////////////////////////////////////////////////////////////////////////
  // Sudakov form factor for one fixed value of cms energy squared s.
  // Note that, when binning in s and b, these values have already been
  // set in the FillTables loops.
  //////////////////////////////////////////////////////////////////////////
  double pt2last = p_pt2bins->x(p_pt2bins->m_nbins-1);
  double sigma, pt2, dpt2, sigmalast = 0., integral = 0.;
  for (int pt2bin=p_pt2bins->m_nbins-1;pt2bin>=0;pt2bin--) {
    pt2       = p_pt2bins->x(pt2bin);
    dpt2      = pt2last-pt2;
    sigma     = p_processes->dSigma(pt2);
    ////////////////////////////////////////////////////////////////////////
    // The dSigma is in 1/GeV^4, norm (the non-diffractive cross section)
    // is in  1/GeV^2 so overall integral over pt^2 does not have any units.
    ////////////////////////////////////////////////////////////////////////
    integral += (sigma+sigmalast)/2. * dpt2/m_norm * m_bJac;
    pt2last   = pt2;
    sigmalast = sigma;
    FillBins(pt2bin,sigma,integral);
  }
  return integral;
}

void Sudakov_Argument::
FillBins(const size_t & pt2bin,const double & sigma, const double & integral) {
  if (m_variable_s && m_variable_b) {
    p_function3->Fill(m_sbin, m_bbin, pt2bin, sigma);
    p_integral3->Fill(m_sbin, m_bbin, pt2bin, integral);
  }
  else if (m_variable_s) {
    p_function2->Fill(m_sbin, pt2bin, sigma);
    p_integral2->Fill(m_sbin, pt2bin, integral);
  }
  else if (m_variable_b) {
    p_function2->Fill(m_bbin, pt2bin, sigma);
    p_integral2->Fill(m_bbin, pt2bin, integral);
  }
  else {
    p_function1->Fill(pt2bin, sigma);
    p_integral1->Fill(pt2bin, integral);
  }
}

void Sudakov_Argument::OutputBin(double & integral, const size_t & pt2bin) {
  if (m_variable_s && m_variable_b) {
    integral += p_integral3->Value(m_sbin, m_bbin, 0);
    msg_Info()<<std::setprecision(6)<<std::setw(26)
	      <<p_integral3->Value(m_sbin,m_bbin, 0)<<" | "
	      <<std::setprecision(6)<<std::setw(20)<<integral<<" |\n";
  }
  else if (m_variable_s) {
    integral += p_integral2->Value(m_sbin, 0);
    msg_Info()<<std::setprecision(6)<<std::setw(32)
	      <<p_integral2->Value(m_sbin, 0)<<"\n";
  }
  else if (m_variable_b) {
    integral += p_integral2->Value(m_bbin, 0);
    msg_Info()<<std::setprecision(6)<<std::setw(29)
	      <<p_integral2->Value(m_bbin, 0)<<" | "
	      <<std::setprecision(6)<<std::setw(28)<<integral<<" |\n";
  }
}

void Sudakov_Argument::OutputTables()
{
  msg_Info() << "   " << std::string(100, '-') << "\n";
  if (m_variable_s && m_variable_b)
    msg_Info() << std::setw(10) << "E_{c.m.} [GeV]" << " | " << std::setw(10)
               << "b [fm]" << " | ";
  if (m_variable_s) msg_Info() << std::setw(23) << "E_{c.m.} [GeV]" << " | ";
  if (m_variable_b) msg_Info() << std::setw(10) << "b [fm]" << " | ";
  msg_Info() << std::setw(10) << "xs_hard/xs_ND" << " | " << std::setw(10)
             << "pt" << " | " << std::setw(15) << "f(pt^2)" << " |  "
             << std::setw(15) << "Int(pt^2)" << " | " << std::setw(10)
             << "Int(pt^2) * sigma_ND" << "\n"
             << std::fixed << std::setprecision(4);
  size_t sbinmax = (m_variable_s ? p_sbins->m_nbins : 1);
  size_t bbinmax = (m_variable_b ? p_bbins->m_nbins : 1);
  double s = p_processes->S(), b = 0;
  for (size_t sbin = 0; sbin < sbinmax; sbin++) {
    s = p_sbins->x(sbin);
    (*p_processes->GetXSecs())(s);
    p_processes->UpdateS(s);
    double xsnd = (p_processes->GetXSecs()->XSndNorm() *
                   p_processes->GetXSecs()->XSnd());
    double xshd = 0., xsratio;
    for (size_t bbin = 0; bbin < bbinmax; bbin++) {
      if (!m_variable_b) b = 0.;
      else {
        b = p_bbins->x(bbin);
        p_processes->SetB(b);
      }
      for (int bin = 0; bin < p_pt2bins->m_nbins; bin -= 10) {
        double pt = sqrt(p_pt2bins->x(bin));
        if (m_variable_s && m_variable_b) {
          msg_Info() << std::setw(10) << s << " | " << std::setw(10) << b
                     << " | ";
          xsratio = p_processes->XSratio(s);
        }
        if (m_variable_s) {
          msg_Info() << std::setw(15) << sqrt(s) << " | ";
          xsratio = p_processes->XSratio(s);
        }
        if (m_variable_b) {
          msg_Info() << std::setw(15) << b << " | ";
          xsratio = p_processes->XSratio(b);
        }
        msg_Info() << std::setw(15) << xsratio << " | " << std::setw(14)
                   << sqrt(p_pt2bins->x(bin)) << " | ";
        if (m_variable_s && m_variable_b)
          msg_Info() << std::setw(10) << p_function3->Value(sbin, bbin, bin)
                     << " | " << std::setw(10)
                     << p_integral3->Value(sbin, bbin, bin) << " | "
                     << std::setw(10)
                     << (p_integral3->Value(sbin, bbin, bin) * xsnd);
        if (m_variable_b)
          msg_Info() << std::setw(10) << p_function2->Value(bbin, bin) << " | "
                     << std::setw(10) << p_integral2->Value(bbin, bin) << " | "
                     << std::setw(10) << (p_integral2->Value(bbin, bin) * xsnd);
        if (m_variable_s) {
          msg_Info() << std::setw(10) << p_function2->Value(sbin, bin) << " | "
                     << std::setw(10) << p_integral2->Value(sbin, bin) << " | "
                     << std::setw(10) << (p_integral2->Value(sbin, bin) * xsnd);
        } else {
          msg_Info() << std::setw(10) << p_function1->Value(bin) << " | "
                     << std::setw(10) << p_integral1->Value(bin) << " | "
                     << std::setw(10) << (p_integral1->Value(bin) * xsnd);
        }
        msg_Info() << "\n";
      }
    }
    msg_Info() << "Calculated look-up tables and values for the"
               << "Sudakov_Argument:\n"
               << "     (non-diffractive xsec = " << xsnd
               << " GeV^-2 = " << (xsnd * rpa->Picobarn() / 1e9) << " mb, "
               << "total xsec = " << p_processes->GetXSecs()->XStot()
               << " mb vs \n"
               << "      hard partonix xsec   = " << xshd << " GeV^-2 " << "= "
               << (xshd * rpa->Picobarn() / 1e9) << " mb).\n\n"
               << "   " << std::string(100, '-') << "\n";
  }
  msg_Info() << "   " << std::string(100, '-') << "\n";
}
