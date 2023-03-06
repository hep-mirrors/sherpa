#include "AMISIC++/Tools/Interaction_Probability.H"
#include "AMISIC++/Perturbative/MI_Processes.H"
#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace AMISIC;
using namespace ATOOLS;


/////////////////////////////////////////////////////////////////////////////////
// All equations in this file refer to 
// Sjostrand-van der Zijl, PRD 36 (1987) 2019.
/////////////////////////////////////////////////////////////////////////////////

Interaction_Probability::Interaction_Probability() :
  m_test(false) {}

Interaction_Probability::~Interaction_Probability() {
  delete p_k;
  delete p_integral;
  delete p_expO;
  delete p_fc;
}

void Interaction_Probability::Initialize(REMNANTS::Remnant_Handler * remnant_handler,
					 MI_Processes * processes) {
  m_mo.Initialize(remnant_handler);
  axis sbins = processes->GetSudakov()->GetSbins();
  p_k        = new OneDim_Table(sbins);
  p_integral = new OneDim_Table(sbins);
  p_expO     = new OneDim_Table(sbins);
  p_fc       = new OneDim_Table(sbins);
  FixK(processes);
  FixOExp();
  if (m_test) OutputTables(processes);
}

void Interaction_Probability::FixK(MI_Processes * processes) {
  /////////////////////////////////////////////////////////////////////////////////
  // Fix the prefactor(s) k in the impact-parameter dependent interaction
  // probability P_int(b) = 1-exp[-k O(b)] with O(b) the matter overlap, by
  // demanding that
  //         [k int d^2b O(b)]/[int d^2b P_int(b)] = sigma_hard/sigma_ND
  // and solving iteratively for k with Newton-Raphson.
  // We fill two look-up tables here: the s-dependent k values and the
  // equally s-dependent (through k) int d^2b P_int(b).
  /////////////////////////////////////////////////////////////////////////////////
  axis sbins = p_k->GetAxis();
  for (size_t bin=0;bin<sbins.m_nbins;bin++) {
    double s       = sbins.x(bin);
    double xsratio = processes->GetSudakov()->XSratio(s);
    double k       = Max(0., NewtonRaphson(xsratio));
    p_k->Fill(bin,k);
    p_integral->Fill(bin,Integral(k, 0));
  }
}

void Interaction_Probability::FixOExp() {
  /////////////////////////////////////////////////////////////////////////////////
  // Filling two more look-up tables: <O> and fc, both depend on s
  /////////////////////////////////////////////////////////////////////////////////
  axis sbins = p_k->GetAxis();
  for (size_t bin=0;bin<sbins.m_nbins;bin++) {
    double k     = p_k->Value(bin);
    double intP  = p_integral->Value(bin);
    double intOP = Integral(k, 2);
    p_expO->Fill(bin, intP>1.e-12 ? intOP/intP : 0.);
    p_fc->Fill(bin,   intP>1.e-12 ? intOP/m_mo.Integral() : 0.);
  }
}

double Interaction_Probability::NewtonRaphson(const double & ratio) {
  /////////////////////////////////////////////////////////////////////////////////
  // Newton-Raphson method to find the solution for k in Eq. (26)
  /////////////////////////////////////////////////////////////////////////////////
  double k = 1.0, f0, f1;
  do {
    double intP0 = Integral(k,0); // b-integral of   {1-exp[-k O(b)]} 
    double intP1 = Integral(k,1); // b-integral of   O(b) exp[-k O(b)] 
    f0 = k*m_mo.Integral()/intP0 - ratio;          
    f1 = m_mo.Integral()*(intP0 - k*intP1)/sqr(intP0); 
    k -= f0/f1;
    if (intP0<=1.e-12) return 0.;
  } while (dabs(f0/f1)>1.e-6 && k>0.);
  return k;
}

double Interaction_Probability::Integral(const double & k,const int & diff) {
  /////////////////////////////////////////////////////////////////////////////////
  // Integrals to be calculated:
  // diff = 0: int d^2b P_int(b),           denominator in Eqs.(26), (32)
  // diff = 1: int d^2b O(b) exp[-k O(b)],  necessary for Newton-Raphson method
  // diff = 2: int d^2b O(b) P_int(b),      numrtator in Eq. (31) 
  /////////////////////////////////////////////////////////////////////////////////
  if (diff==0) {
    P_Integrand p(&m_mo,k);
    Gauss_Integrator integrator(&p);
    return integrator.Integrate(0.,m_mo.Bmax(),1.e-8,1);
  }
  else if (diff==1) {
    OtimesExp_Integrand oe(&m_mo,k);
    Gauss_Integrator integrator(&oe);
    return integrator.Integrate(0.,m_mo.Bmax(),1.e-8,1);
  }
  else if (diff==2) {
    OtimesP_Integrand op(&m_mo,k);
    Gauss_Integrator integrator(&op);
    return integrator.Integrate(0.,m_mo.Bmax(),1.e-8,1);
  }
  return 0.;
}

void Interaction_Probability::OutputTables(MI_Processes * processes) {
  axis sbins = p_k->GetAxis();
  msg_Out()<<"-------------------------------------------------------------------------------\n"
	   <<"Calculated look-up tables and values for the Interaction Probability:\n"
	   <<std::setw(15)<<"E_{c.m.} [GeV]"<<" | "
	   <<std::setw(15)<<"xs_hard/xs_ND"<<" | "
	   <<std::setw(10)<<"k"<<" | "
	   <<std::setw(10)<<"<O>"<<" |  "
	   <<std::setw(10)<<"fc\n"
	   <<std::fixed<<std::setprecision(4);
  for (size_t bin=0;bin<sbins.m_nbins;bin++) {
    double s       = sbins.x(bin);
    double xsratio = processes->GetSudakov()->XSratio(s);
    msg_Out()<<std::setw(15)<<sqrt(s)<<" | "
	     <<std::setw(15)<<xsratio<<" | "
	     <<std::setw(10)<<p_k->Value(bin)<<" | "
	     <<std::setw(10)<<p_expO->Value(bin)<<" | "
	     <<std::setw(10)<<p_fc->Value(bin)<<"\n";
  }
  msg_Out()<<"-------------------------------------------------------------------------------\n";
  exit(1);
}

double P_Integrand::operator()(double b) {
  /////////////////////////////////////////////////////////////////////////////////
  // Integrand for d^2b [1-exp(- k O(b)] where O(b) is the time-integrated
  // matter overlap, being the tricky part of the denominator in Eq.(26).
  /////////////////////////////////////////////////////////////////////////////////
  return 2.*M_PI*b*(1. - exp(-m_k*(*p_mo)(b)));
}

double OtimesExp_Integrand::operator()(double b) {
  /////////////////////////////////////////////////////////////////////////////////
  // Integrand for d^2b O(b) exp(- k O(b)] where O(b) is the time-integrated
  // matter overlap, used by the Newton-Raphson method
  /////////////////////////////////////////////////////////////////////////////////
  return 2.*M_PI*b*(*p_mo)(b)*exp(-m_k*(*p_mo)(b));
}

double OtimesP_Integrand::operator()(double b) {
  /////////////////////////////////////////////////////////////////////////////////
  // Integrand for d^2b O(b) {1 - exp(- k O(b)]} where O(b) is the time-integrated
  // matter overlap, used in Eqs. (29) and (31) to calculate <o> and fc.
  /////////////////////////////////////////////////////////////////////////////////
  return 2.*M_PI*b*(*p_mo)(b)*(1-exp(-m_k*(*p_mo)(b)));
}







