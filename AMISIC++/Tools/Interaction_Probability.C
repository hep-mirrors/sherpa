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
  p_mo(new Matter_Overlap()), m_prefK(1.), m_integral(0.),
  m_output(false)
{}

Interaction_Probability::~Interaction_Probability() {
  delete p_mo;
  delete p_prefKs;
  delete p_oexpvalues;
  delete p_integrals;
}

void Interaction_Probability::Initialize(MI_Processes * processes) {
  p_mo->Initialize();
  FixPrefactor(processes);
  CalculateOExpValue(processes);
}

double Interaction_Probability::operator()(const double & b) {
  /////////////////////////////////////////////////////////////////////////////////
  // Interaction probability, Eq. (24)
  /////////////////////////////////////////////////////////////////////////////////
  return 1.-exp(-m_prefK * (*p_mo)(b));
}

void Interaction_Probability::CalculateIntegral() {
  /////////////////////////////////////////////////////////////////////////////////
  // Integral int d^2b P_int(b), denominator in Eqs.(26), (32)
  /////////////////////////////////////////////////////////////////////////////////
  IP_Integrand ipint(this);
  Gauss_Integrator integrator(&ipint);
  m_integral  = integrator.Integrate(0.,p_mo->Bmax(),1.e-8,1);
  if (m_output) msg_Out()<<"     --> "<<METHOD<<"(bmax = "<<p_mo->Bmax()<<") "
			 <<"yields integral "<<m_integral<<".\n";
}

void Interaction_Probability::CalculateBNorm() {
  IP_Integrand bipint(this,2);
  Gauss_Integrator integrator(&bipint);
  m_integralB = integrator.Integrate(0.,p_mo->Bmax(),1.e-8,1);
  m_bnorm     = m_integralB/m_integral;
  if (m_output) msg_Out()<<"     --> "<<METHOD<<"(bmax = "<<p_mo->Bmax()<<") "
			 <<"yields norm "<<m_bnorm<<".\n";
}

void Interaction_Probability::CalculateOExpValue(MI_Processes * processes) {
  /////////////////////////////////////////////////////////////////////////////////
  // Integral of int d^2b O(b) P_int(b), numerator of f_c in Eq.(31)
  /////////////////////////////////////////////////////////////////////////////////
  O_ExpV_Integrand oexpvint(this);
  Gauss_Integrator integrator(&oexpvint);
  m_oexpvalue = integrator.Integrate(0.,p_mo->Bmax(),1.e-8,1)/m_integral;
  if (m_output) msg_Out()<<"     --> "<<METHOD<<"(bmax = "<<p_mo->Bmax()<<") "
			 <<"yields <value> = "<<m_oexpvalue<<".\n";
}

bool Interaction_Probability::FixPrefactor(MI_Processes * processes) {
  /////////////////////////////////////////////////////////////////////////////////
  // In this method we fix the prefactor(s) k in the impact-parameter dependent
  // interaction probability P_int(b) = 1-exp[-k O(b)] where O(b) is the matter
  // overlap.  This is done by demanding that
  //         [k int d^2b O(b)]/[int d^2b P_int(b)] = sigma_hard/sigma_ND
  // and solving iteratively numerically for k.
  /////////////////////////////////////////////////////////////////////////////////
  axis sbins   = processes->GetSudakov()->GetSbins();
  p_prefKs     = new OneDim_Table(sbins);
  p_oexpvalues = new OneDim_Table(sbins);
  p_integrals  = new OneDim_Table(sbins);
  msg_Out()<<"   * "<<METHOD<<"(for "<<sbins.m_nbins<<" bins "
	   <<"in s in ["<<sbins.m_xmin<<", "<<sbins.m_xmax<<"])\n";
  double s, xsratio, faclow, fachigh, reslow, reshigh, deltafac, deltares;
  double overlap = p_mo->Integral(); 
  for (size_t sbin=0;sbin<sbins.m_nbins;sbin++) {
    s       = sbins.x(sbin);   
    xsratio = processes->XSratio(s);
    msg_Out()<<"     - "<<METHOD<<"(s = "<<s<<", ratio = "<<xsratio<<", "
	     <<"xs_hard = "<<(xsratio*processes->GetXSecs()->XSnd(s)*rpa->Picobarn())<<" pb, "
	     <<"xs_nd = "<<(processes->GetXSecs()->XSnd(s)*rpa->Picobarn())<<" pb).\n";
    if (xsratio==0.) {
      p_prefKs->Fill(s,fachigh);
      p_integrals->Fill(s,m_integral);
    }
    else {
      m_prefK = faclow = 1.;
      CalculateIntegral();
      reslow  = faclow * p_mo->Integral()/m_integral;
      fachigh = 50.;
      do {
	m_prefK = fachigh;
	CalculateIntegral();
	reshigh  = fachigh * overlap/m_integral;
	if (m_output)
	  msg_Out()<<"         . k = ["<<faclow<<", "<<fachigh<<"] --> "
		   <<"res = ["<<reslow<<", "<<reshigh
		   <<"] from integral = "<<m_integral<<".\n";
	deltafac = fachigh-faclow;
	deltares = reshigh-reslow;
	faclow   = fachigh;
	reslow   = reshigh;
	fachigh += deltafac/deltares * (xsratio-reshigh);
      } while (dabs(1.-reshigh/xsratio)>1.e-8);
      m_prefK = fachigh;
      CalculateIntegral();
      p_prefKs->Fill(s,fachigh);
      p_integrals->Fill(s,m_integral);
      msg_Out()<<"       geometric rescaling factor = "<<m_prefK<<" yields "
	       <<"sigma/sigmaND = "<<(m_prefK*overlap/m_integral)<<".\n";
      CalculateBNorm();
    }
  }
  return true; 
}

double IP_Integrand::operator()(double b) {
  /////////////////////////////////////////////////////////////////////////////////
  // Integrand for d^2b O(b) = 2 pi b db O(b), where O(b) is the time-integrated
  // matter overlap, being the tricky part of the numerator in Eq.(xxx).
  /////////////////////////////////////////////////////////////////////////////////
  return 2.*M_PI*b*(*p_ip)(b);
}

double O_ExpV_Integrand::operator()(double b) {
  /////////////////////////////////////////////////////////////////////////////////
  // Integrand for d^2b O(b) P_int(b) = 2 pi b db O(b) P_int(b) , where O(b) is the
  // time-integrated matter overlap and P_int(b) is the interaction probability, given
  // by 1-exp[-k O(b)].  This is the tricky part of the numerator in Eq.(xxx).
  /////////////////////////////////////////////////////////////////////////////////
  return 2.*M_PI*b * (*p_mo)(b) * (*p_ip)(b);
}

