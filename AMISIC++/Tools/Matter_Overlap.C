#include "AMISIC++/Tools/Matter_Overlap.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

using namespace AMISIC;
using namespace ATOOLS;


// All equations in this file refer to 
// Sjostrand-van der Zijl, PRD 36 (1987) 2019.

Matter_Overlap::Matter_Overlap() :
  m_bstep(0.), m_bmax(0.), m_integral(0.), m_norm(1./M_PI) ///(4.*M_PI*M_PI))
{}

Matter_Overlap::~Matter_Overlap() {}

void Matter_Overlap::Initialize() {
  InitializeFormFactors();
  CalculateIntegral();
  msg_Tracking()<<METHOD<<"(form = "<<m_overlapform<<" --> r = "<<m_radius12<<"), "
		<<"integral = "<<m_integral<<" norm = "<<(m_norm*m_norm1)<<".\n";
}

double Matter_Overlap::operator()(double b) {
  // Matter overlap in two forms available, but only Single_Gaussian fully
  // functional at the moment.
  switch (m_overlapform) {
    case overlap_form::Single_Gaussian:
      return m_norm * m_norm1 * exp(-b*b/m_radius12);
    case overlap_form::Double_Gaussian:
      double b2(b*b);
      return m_norm * (m_norm1 * exp(-b2/m_radius12) +
		       m_norm2 * exp(-b2/m_radius22) +
		       m_norm3 * exp(-b2/m_radius32));
  }
}

double Matter_Overlap::SelectB() const {
  double b, b2, radius;
  switch (m_overlapform) {
    case overlap_form::Single_Gaussian:
      radius = m_radius1;
      break;
    case overlap_form::Double_Gaussian:
      double rand = ran->Get();
      if ((rand-=sqr(m_fraction1))<=0.)        radius = m_radius1;
      else if ((rand-=sqr(1-m_fraction1))<=0.) radius = m_radius2;
      else                                     radius = m_radius3;
      break;
  }
  do {
    auto b = dabs(ran->GetGaussian())*radius;
    if (b>m_bmax) b = dabs(ran->GetGaussian())*radius;
  } while (b>m_bmax);
  return b;
}


void Matter_Overlap::InitializeFormFactors() {
  // Matter overlap in two forms available, but only Single_Gaussian fully
  // functional at the moment.
  m_overlapform = mipars->GetOverlapForm();
  switch (m_overlapform) {
    case overlap_form::Single_Gaussian:
      m_fraction1 = 1.;
      m_radius1   = (*mipars)("Matter_Radius1");
      m_radius12  = 2.*sqr(m_radius1);
      m_norm1     = 1./m_radius12;
      m_bstep     = m_radius1/100.;
      break;
    case overlap_form::Double_Gaussian:
      m_fraction1 = (*mipars)("Matter_Fraction1");
      m_radius1   = (*mipars)("Matter_Radius1");
      m_radius12  = 2.*sqr(m_radius1);
      m_radius2   = (*mipars)("Matter_Radius2");
      m_radius22  = 2.*sqr(m_radius2);
      m_radius32  = (m_radius12+m_radius22)/2.;
      m_radius3   = sqrt(m_radius32);
      m_norm1     = sqr(m_fraction1)/m_radius12;
      m_norm2     = sqr(1.-m_fraction1)/m_radius22;
      m_norm3     = 2.*m_fraction1*(1.-m_fraction1)/m_radius32;
      m_bstep     = Min(m_radius1,m_radius2)/100.;
      break;
  }
}
  
void Matter_Overlap::CalculateIntegral() {
  // Integral int d^2b O(b), numerator Eq.(32)
  MO_Integrand moint(this);
  Gauss_Integrator integrator(&moint);
  double bmin = 0., bstep = m_bstep, previous = 0., previous2 = 0., result = 0.;
  do {
    result  += previous = integrator.Integrate(bmin,bmin+bstep,1.e-8,1);
    bmin    += bstep;
  } while (dabs(previous/result)>1.e-10);
  m_bmax     = bmin;
  m_integral = result;
}

double MO_Integrand::operator()(double b) {
  // Integrand for d^2b O(b) = 2 pi b db O(b), where O(b) is the time-integrated
  // matter overlap, being the tricky part of the numerator in Eq.(32).
  // This does not include the prefactor k.
  return 2.*M_PI*b*(*p_mo)(b);
}

