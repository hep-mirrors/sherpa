#include "AMISIC++/Tools/Interaction_Probability.H"
#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Org/Message.H"

using namespace AMISIC;
using namespace ATOOLS;


// All equations in this file refer to 
// Sjostrand-van der Zijl, PRD 36 (1987) 2019.

Interaction_Probability::Interaction_Probability() :
  p_mo(new Matter_Overlap()), m_prefK(1.), m_integral(0.)
{}

void Interaction_Probability::Initialize(const double & xsecratio) {
  p_mo->Initialize();
  m_bmax = p_mo->Bmax();
  FixPrefactor(xsecratio);
  CalculateOExpValue();
}

double Interaction_Probability::operator()(const double & b) {
  return 1.-exp(-m_prefK * (*p_mo)(b));
}

void Interaction_Probability::CalculateIntegral() {
  // Integral int d^2b P_int(b), denominator in Eqs.(26), (32)
  IP_Integrand ipint(this);
  Gauss_Integrator integrator(&ipint);
  m_integral = integrator.Integrate(0.,p_mo->Bmax(),1.e-8,1);
}

void Interaction_Probability::CalculateOExpValue() {
  // Integral of int d^2b O(b) P_int(b), numerator of f_c in Eq.(31)
  O_ExpV_Integrand oexpvint(this);
  Gauss_Integrator integrator(&oexpvint);
  m_oexpvalue = integrator.Integrate(0.,p_mo->Bmax(),1.e-8,1)/m_integral;
  msg_Tracking()<<METHOD<<" yields "<<m_oexpvalue<<".\n";
}

bool Interaction_Probability::FixPrefactor(const double & xsecratio) {
  // In this method we fix the prefactor k in the impact-parameter dependent
  // interaction probability P_int(b) = 1-exp[-k O(b)] where O(b) is the matter
  // overlap.  This is done by demanding that
  //         [int_b O(b)]/[int_b db P_int(b)] = sigma_ND/sigma_tot
  // and solving iteratively numerically for k.
  DEBUG_FUNC(xsecratio);
  if (xsecratio<=1.) return false;
  double faclow  = 1.;
  SetPrefactor(faclow);
  CalculateIntegral();
  double reslow  = faclow * p_mo->Integral()/m_integral;
  double fachigh = 2., reshigh, deltafac, deltares;
  msg_Debugging()<<"Start iteration for int(overlap) = "
		 <<p_mo->Integral()<<" aiming for ratio "<<xsecratio<<"\n";
  do {
    SetPrefactor(fachigh);
    CalculateIntegral();
    reshigh  = fachigh * p_mo->Integral()/m_integral;
    msg_Debugging()<<"k = ["<<faclow<<", "<<fachigh<<"] --> "
		   <<"res = ["<<reslow<<", "<<reshigh
		   <<"] from integral = "<<m_integral<<".\n";
    deltafac = fachigh-faclow;
    deltares = reshigh-reslow;
    faclow   = fachigh;
    reslow   = reshigh;
    fachigh += deltafac/deltares * (xsecratio-reshigh);
  } while (dabs(1.-reshigh/xsecratio)>1.e-8);
  SetPrefactor(fachigh);
  msg_Debugging()<<"==> geometric rescaling factor = "<<m_prefK<<".\n";
  return true; 
}

double IP_Integrand::operator()(double b) {
  // Integrand for d^2b O(b) = 2 pi b db O(b), where O(b) is the time-integrated
  // matter overlap, being the tricky part of the numerator in Eq.(xxx).
  return 2.*M_PI*b*(*p_ip)(b);
}

double O_ExpV_Integrand::operator()(double b) {
  // Integrand for d^2b O(b) P_int(b) = 2 pi b db O(b) P_int(b) , where O(b) is the
  // time-integrated matter overlap and P_int(b) is the interaction probability, given
  // by 1-exp[-k O(b)].  This is the tricky part of the numerator in Eq.(xxx).
  return 2.*M_PI*b * (*p_mo)(b) * (*p_ip)(b);
}

