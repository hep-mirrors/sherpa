#include "AMISIC++/Tools/Interaction_Probability.H"
#include "AMISIC++/Perturbative/MI_Processes.H"
#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace AMISIC;
using namespace ATOOLS;


/////////////////////////////////////////////////////////////////////////////
// All equations in this file refer to
// Sjostrand-van der Zijl, PRD 36 (1987) 2019 or
// Corke-Sjostrand, JHEP 05 (2011) 009.
/////////////////////////////////////////////////////////////////////////////

Interaction_Probability::Interaction_Probability() :
  p_k(NULL), p_integral(NULL), p_expO(NULL), p_fc(NULL),
  m_test(false) {}

Interaction_Probability::~Interaction_Probability() {
  if (p_k)        delete p_k;
  if (p_integral) delete p_integral;
  if (p_expO)     delete p_expO;
  if (p_fc)       delete p_fc;
}

void Interaction_Probability::
Initialize(Matter_Overlap * mo,MI_Processes * processes,axis * sbins) {
  p_mo = mo;
  if (p_mo->IsDynamic()) InitializeDynamic(processes,sbins);
  else InitializeStatic(processes,sbins);
}

void Interaction_Probability::
InitializeDynamic(MI_Processes * processes,axis * sbins) {
  ///////////////////////////////////////////////////////////////////////////
  // This fixes the potentially s-dependent (re-)scaling factor(s) for the
  // radius of the dynamic matter overlap and tabulates them in p_k.
  // p_k:        solving for k Eq. (CS, 20):
  //             sigma_ND = int_0^bmax db [ 2.*M_PI*b * (1. - exp(-n(b, k)))],
  //             where k modifies the radius as R -> k*R and where
  //             n(b,k) is given by Eq. (CS, 18) with k entering as above.
  ///////////////////////////////////////////////////////////////////////////
  nbar_DynIntegrand nbar(processes,p_mo->GetBBins());
  Gauss_Integrator  gauss(&nbar);
  p_k = new OneDim_Table(*sbins);
  msg_Info()<<"   "<<std::string(77,'-')<<"\n"
	    <<"   | "<<METHOD<<":"<<std::string(36,' ')<<"|\n"
	    <<"   | fixing dynamic overlap radius factors 1/k "
	    <<"to reproduce sigma_ND           |\n";
  for (size_t sbin=0;sbin<sbins->m_nbins;sbin++) {
    double s    = sbins->x(sbin);
    (*processes->GetXSecs())(s);
    double xsnd = processes->GetXSecs()->XSnd();
    double k    = 1;
    size_t iter = 0;
    msg_Info()<<"   | E = "<<sqrt(s)<<" GeV, physical sigma_ND = "
	      <<std::setprecision(6)<<std::setw(10)
	      <<(xsnd*rpa->Picobarn()/1.e9)
	      <<" mb. "<<std::string(25,' ')<<"|\n"
	      <<"   "<<std::string(77,'-')<<"\n"
	      <<"   |   iter |               k |   bmax [fm] |"
	      <<"       sigma_ND [mb] (from P_int) |\n";
    do {
      p_mo->FixKRadius(k);
      double bmax     = nbar.Initialize();
      double d2b_Pint = gauss.Integrate(0.,Min(bmax,p_mo->Bmax()),1.e-8,1);
      std::setprecision(6);
      msg_Info()<<"   | "<<std::setprecision(3)<<std::setw(6)<<iter<<" | "
		<<std::setprecision(6)<<std::setw(15)<<k<<" | "
		<<std::setprecision(6)<<std::setw(11)<<bmax<<" | "
		<<std::setprecision(6)<<std::setw(32)
		<<(d2b_Pint*rpa->Picobarn()/1.e9)<<" |\n";
      if (dabs(1.-d2b_Pint/xsnd)<0.1) break;
      k  *= xsnd/d2b_Pint;
    } while (iter++<10);
    p_k->Fill(sbin,k);
  }
  msg_Info()<<"   "<<std::string(77,'-')<<"\n";
}

void Interaction_Probability::
InitializeStatic(MI_Processes * processes,axis * sbins) {
  ///////////////////////////////////////////////////////////////////////////
  // This fixes all helper functions relevant for the case of static form
  // factors and, thereby, kinematic-independent matter overlaps.
  // We initialize the following tables here:
  // p_k:        solving for k the equation (SZ, 26)
  //             sigma_hard/sigma_ND =
  //             int_0^bmax db [ 2.*M_PI*b * k*O(b)]/
  //             int_0^bmax db [ 2.*M_PI*b * (1. - exp(-k*O(b)))];
  // p_integral: int_0^bmax db [ 2.*M_PI*b * (1. - exp(-k*O(b)))]
  //             the integral over impact parameter space of the
  //             Interaction_Probaility, P_int (24);
  // p_expO:     int_0^bmax db [ 2.*M_PI*b * O(b)*(1-exp(-k*O(b)))] /
  //             int_0^bmax db [ 2.*M_PI*b * (1. - exp(-k*O(b)))],
  //             the expectation value <O> of the Matter_Overlap, (29);
  // p_fc:       int_0^bmax db [ 2.*M_PI*b * O(b)*(1-exp(-k*O(b)))] /
  //             int_0^bmax db [ 2.*M_PI*b * O(b)]
  //             the f_c of (31).
  // In all cases O(b) is obtained from the Matter_Overlap p_mo and k
  // is taken from the solution of Eq (SZ, 26) for k.
  ///////////////////////////////////////////////////////////////////////////
  p_k        = new OneDim_Table(*sbins);
  p_integral = new OneDim_Table(*sbins);
  p_expO     = new OneDim_Table(*sbins);
  p_fc       = new OneDim_Table(*sbins);
  FixK4Static(processes,sbins);
  FixOExp(sbins);
  msg_Info()<<"   "<<std::string(77,'-')<<"\n"
	    <<"   | "<<METHOD<<" for "<<sbins->m_nbins<<" bins in s: "
	    <<std::string(15,' ')<<"|\n"
	    <<"   | k         = "<<std::setprecision(6)<<std::setw(10)
	    <<p_k->Value(0)<<std::string(52,' ')<<"|\n"
	    <<"   | integral  = "<<std::setprecision(6)<<std::setw(10)
	    <<p_integral->Value(0)<<std::string(52,' ')<<"|\n"
	    <<"   | exp[O]    = "<<std::setprecision(6)<<std::setw(10)
	    <<p_expO->Value(0)<<std::string(52,' ')<<"|\n"
	    <<"   | fc        = "<<std::setprecision(6)<<std::setw(10)
	    <<p_fc->Value(0)<<std::string(52,' ')<<"|\n"
	    <<"   | int d2b O = "<<std::setprecision(6)<<std::setw(10)
	    <<(*p_mo).Integral()
	    <<std::string(52,' ')<<"|\n"
	    <<"   | ratio     = "<<std::setprecision(6)<<std::setw(10)
	    <<(p_k->Value(0)*(*p_mo).Integral()/p_integral->Value(0))
	    <<std::string(52,' ')<<"|\n";
  OutputTables(processes,sbins);
  if (m_test) {
    OutputTables(processes,sbins);
    THROW(normal_exit,"testing complete");
  }
}

void Interaction_Probability::
FixK4Static(MI_Processes * processes,axis * sbins) {
  //////////////////////////////////////////////////////////////////////////
  // Fix the prefactor(s) k in the impact-parameter dependent interaction
  // probability P_int(b) = 1-exp[-k O(b)] with O(b) the matter overlap, by
  // demanding that, cf Eq (SZ, 26),
  //         [k int d^2b O(b)]/[int d^2b P_int(b)] = sigma_hard/sigma_ND
  // and solving iteratively for k with Newton-Raphson.
  // We fill two look-up tables here: the s-dependent k values and the
  // equally s-dependent (through k) int d^2b P_int(b).
  //////////////////////////////////////////////////////////////////////////
  for (size_t bin=0;bin<sbins->m_nbins;bin++) {
    double s       = sbins->x(bin);
    double xsratio = processes->XSratio(s);
    double k       = Max(0., NewtonRaphson(xsratio));
    p_k->Fill(bin,k);
    p_integral->Fill(bin,Integral(k, 0));
  }
}

void Interaction_Probability::FixOExp(axis * sbins) {
  //////////////////////////////////////////////////////////////////////////
  // Filling two more look-up tables: <O> from Eq. (SZ, 29) and fc from
  // Eq. (SZ, 31), both depend on s.
  //////////////////////////////////////////////////////////////////////////
  for (size_t bin=0;bin<sbins->m_nbins;bin++) {
    double k     = p_k->Value(bin);
    double intP  = p_integral->Value(bin);
    double intOP = Integral(k, 2);
    p_expO->Fill(bin, intP>1.e-12 ? intOP/intP : 0.);
    p_fc->Fill(bin,   intP>1.e-12 ? intOP/(*p_mo).Integral() : 0.);
  }
}

double Interaction_Probability::NewtonRaphson(const double & ratio) {
  //////////////////////////////////////////////////////////////////////////
  // Newton-Raphson method to find the solution for k in Eq. (26)
  //////////////////////////////////////////////////////////////////////////
  double k = 1.0, f0, f1;
  do {
    double intP0 = Integral(k,0); // b-integral of   {1-exp[-k O(b)]}
    double intP1 = Integral(k,1); // b-integral of   O(b) exp[-k O(b)]
    f0 = k*(*p_mo).Integral()/intP0 - ratio;
    f1 = (*p_mo).Integral()*(intP0 - k*intP1)/sqr(intP0);
    k -= f0/f1;
    if (intP0<=1.e-12) return 0.;
  } while (dabs(f0/f1)>1.e-6 && k>0.);
  return k;
}

double Interaction_Probability::Integral(const double & k,const int & diff) {
  //////////////////////////////////////////////////////////////////////////
  // Integrals to be calculated:
  // diff = 0: int d^2b P_int(b),           denominator in Eqs.(26), (32)
  // diff = 1: int d^2b O(b) exp[-k O(b)],  for Newton-Raphson method
  // diff = 2: int d^2b O(b) P_int(b),      numerator in Eq. (31)
  //////////////////////////////////////////////////////////////////////////
  if (diff==0) {
    P_Integrand p(p_mo,k);
    Gauss_Integrator integrator(&p);
    return integrator.Integrate(0.,p_mo->Bmax(),1.e-8,1);
  }
  else if (diff==1) {
    OtimesExp_Integrand oe(p_mo,k);
    Gauss_Integrator integrator(&oe);
    return integrator.Integrate(0.,p_mo->Bmax(),1.e-8,1);
  }
  else if (diff==2) {
    OtimesP_Integrand op(p_mo,k);
    Gauss_Integrator integrator(&op);
    return integrator.Integrate(0.,p_mo->Bmax(),1.e-8,1);
  }
  return 0.;
}

double P_Integrand::operator()(double b) {
  //////////////////////////////////////////////////////////////////////////
  // Integrand for d^2b [1-exp(- k O(b)] where O(b) is the time-integrated
  // matter overlap, being the tricky part of the denominator in Eq.(26).
  //////////////////////////////////////////////////////////////////////////
  return 2.*M_PI*b*(1. - exp(-m_k*(*p_mo)(b)));
}

double OtimesExp_Integrand::operator()(double b) {
  //////////////////////////////////////////////////////////////////////////
  // Integrand for d^2b O(b) exp(- k O(b)] where O(b) is the time-integrated
  // matter overlap, used by the Newton-Raphson method
  //////////////////////////////////////////////////////////////////////////
  return 2.*M_PI*b*(*p_mo)(b)*exp(-m_k*(*p_mo)(b));
}

double OtimesP_Integrand::operator()(double b) {
  //////////////////////////////////////////////////////////////////////////
  // Integrand for d^2b O(b) {1 - exp(- k O(b)]} where O(b) is the
  // time-integrated matter overlap, used in Eqs. (SZ, 29) and (SZ, 31) to
  // calculate <O> and fc.
  //////////////////////////////////////////////////////////////////////////
  return 2.*M_PI*b*(*p_mo)(b)*(1-exp(-m_k*(*p_mo)(b)));
}

nbar_DynIntegrand::nbar_DynIntegrand(MI_Processes * procs,axis * bbins) :
  //////////////////////////////////////////////////////////////////////////
  // nbar(b) to be used for the interaction probability, Eq. (CS, 13):
  // P_int = d^2b {1-exp[-nbar(b)]}.
  // As the matter overlap is dynamic it becomes part of the partonic cross
  // section calculation which enters nbar(b).  We therefore pre-tabulate
  // and interpolate values for different b's.
  //////////////////////////////////////////////////////////////////////////
  p_procs(procs), m_table(OneDim_Table(*bbins)) {}

double nbar_DynIntegrand::Initialize() {
  //////////////////////////////////////////////////////////////////////////
  // Pre-integrating the (dynamic) matter-overlap dependent parton-level
  // cross section, Eq. (CS, 18), in bins of b.  The such created look-up
  // table is used in the integral over impact parameter of Eq. (CS, 20) which
  // we use to fix the radius of the overlap to reproduce the non-diffractive
  // cross section within 10%.
  //////////////////////////////////////////////////////////////////////////
  double xsMO, xsSum = 0., xsmax, xs = p_procs->SigmaTot(false);
  for (size_t bbin=0;bbin<m_table.GetAxis().m_nbins;bbin++) {
    p_procs->SetB(m_table.GetAxis().x(bbin));
    xsSum += xsMO = p_procs->SigmaTot(true);
    m_table.Fill(bbin,xsMO);
    if (bbin==0) xsmax = xsMO;
    else if (xsMO/xsmax<1.e-6 && bbin<m_table.GetAxis().m_nbins-1) {
      for (size_t i=bbin+1;i<m_table.GetAxis().m_nbins;i++) m_table.Fill(i,0);
      return m_table.GetAxis().x(bbin);
    }
  }
  return m_table.GetAxis().m_xmax;
}

double nbar_DynIntegrand::operator()(double b) {
  //////////////////////////////////////////////////////////////////////////
  // The integrand for the b-integration in Eq. (CS, 20).
  // sigma_ND = d^2b P_int(b) = d^2b {1-exp[-nbar(b)]}.
  // As the matter overlap is dynamic it becomes part of the partonic cross
  // section calculation which enters nbar(b).  We therefore pre-tabulate
  // and interpolate values for different b's.
  //////////////////////////////////////////////////////////////////////////
  return 2.*M_PI*b * (1.-exp(-m_table(b)));
}


void Interaction_Probability::
OutputTables(MI_Processes * processes,axis * sbins) {
  msg_Info()<<"   "<<std::string(77,'-')<<"\n"
	    <<"   | "<<METHOD<<" look-up tables and values:          |\n"
	    <<"   | "<<std::setw(15)<<"E_{c.m.} [GeV]"<<" | "
	    <<std::setw(15)<<"xs_hard/xs_ND"<<" | "
	    <<std::setw(10)<<"k"<<" | "
	    <<std::setw(10)<<"<O>"<<" |  "
	    <<std::setw(10)<<"fc"<<" |\n"
	    <<std::fixed<<std::setprecision(4);
  for (size_t bin=0;bin<sbins->m_nbins;bin++) {
    double s       = sbins->x(bin);
    double xsratio = processes->XSratio(s);
    msg_Info()<<"   | "<<std::setprecision(6)<<std::setw(15)<<sqrt(s)<<" | "
	      <<std::setprecision(6)<<std::setw(15)<<xsratio<<" | "
	      <<std::setprecision(6)<<std::setw(10)<<p_k->Value(bin)<<" | "
	      <<std::setprecision(6)<<std::setw(10)<<p_expO->Value(bin)<<" | "
	      <<std::setprecision(6)<<std::setw(11)<<p_fc->Value(bin)<<" | \n";
  }
  msg_Info()<<"   "<<std::string(77,'-')<<"\n\n";
}

