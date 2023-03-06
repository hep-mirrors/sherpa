#include "AMISIC++/Tools/Matter_Overlap.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

using namespace AMISIC;
using namespace ATOOLS;


/////////////////////////////////////////////////////////////////////////////////
// All equations in this file refer to 
// Sjostrand-van der Zijl, PRD 36 (1987) 2019.
/////////////////////////////////////////////////////////////////////////////////

Matter_Overlap::Matter_Overlap() :
  ///////////////////////////////////////////////////////////////////////////////
  // The norm come from the (implicitly normalised) Gaussian matter distributions
  // of the hadron form factors ~pi^{-3/2} times the pi^2 from the time-integrated
  // overlap when they collide.
  ///////////////////////////////////////////////////////////////////////////////
  m_bstep(0.), m_bmax(0.), m_integral(0.), m_norm(1./M_PI)
{}

Matter_Overlap::~Matter_Overlap() {}

void Matter_Overlap::Initialize(REMNANTS::Remnant_Handler * remnant_handler) {
  InitializeFormFactors();
  CalculateIntegral();
}

double Matter_Overlap::operator()(double b) {
  /////////////////////////////////////////////////////////////////////////////////
  // Matter overlap in two forms available: Single_Gaussian and Double_Gaussian
  /////////////////////////////////////////////////////////////////////////////////
  switch (m_overlapform) {
  case overlap_form::code::Single_Gaussian:
    return m_norm * m_norm1 * exp(-b*b/m_radius12);
  case overlap_form::code::Double_Gaussian:
  default:
    double b2(b*b);
    return m_norm * (m_norm1 * exp(-b2/m_radius12) +
		     m_norm2 * exp(-b2/m_radius22) +
		     m_norm3 * exp(-b2/m_radius32));
  }
}

double Matter_Overlap::SelectB(const bool & mode) const {
  /////////////////////////////////////////////////////////////////////////////////
  // Algorithm:
  // 1. select a radius R according to matter content:
  //    - for single Gaussian, there is no selection to be made
  //    - for double Gaussian, one of the three radii is picked.
  // 2. Select b according to d^2b O(b) = d b^2 exp(-b^2/R^2).
  /////////////////////////////////////////////////////////////////////////////////
  double b, radius;
  switch (m_overlapform) {
    case overlap_form::code::Single_Gaussian:
      radius = m_radius1;
      break;
    case overlap_form::code::Double_Gaussian:
      double rand = ran->Get();
      if ((rand-=sqr(m_fraction1))<=0.)        radius = m_radius1;
      else if ((rand-=sqr(1-m_fraction1))<=0.) radius = m_radius2;
      else                                     radius = m_radius3;
      break;
  }
  /////////////////////////////////////////////////////////////////////////////////
  // b from Matter_Overlap, hence r^2_overlap = 2*r^2_formfactor
  /////////////////////////////////////////////////////////////////////////////////
  if (mode) radius *= sqrt(2.);
  do {
    b = sqrt(-log(Max(1.e-12,ran->Get())))*radius;
  } while (b>m_bmax);
  return b;
}


void Matter_Overlap::InitializeFormFactors() {
  /////////////////////////////////////////////////////////////////////////////////
  // Matter overlap in two forms available
  /////////////////////////////////////////////////////////////////////////////////
  m_overlapform = mipars->GetOverlapForm();
  switch (m_overlapform) {
    case overlap_form::code::Single_Gaussian:
      m_fraction1 = 1.;
      m_radius1   = (*mipars)("Matter_Radius1");
      m_radius12  = sqr(m_radius1);
      m_norm1     = 1./m_radius12;
      m_bstep     = m_radius1/100.;
      break;
    case overlap_form::code::Double_Gaussian:
      m_fraction1 = (*mipars)("Matter_Fraction1");
      m_radius1   = (*mipars)("Matter_Radius1");
      m_radius12  = sqr(m_radius1);
      m_radius2   = (*mipars)("Matter_Radius2");
      m_radius22  = sqr(m_radius2);
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
  /////////////////////////////////////////////////////////////////////////////////
  // Integral int d^2b O(b), numerator Eq.(32)
  /////////////////////////////////////////////////////////////////////////////////
  MO_Integrand moint(this);
  Gauss_Integrator integrator(&moint);
  double bmin = 0., bstep = m_bstep, previous, result = 0.;
  do {
    result  += previous = integrator.Integrate(bmin,bmin+bstep,1.e-8,1);
    bmin    += bstep;
  } while (dabs(previous/result)>1.e-10);
  m_bmax     = bmin;
  m_integral = result;
}

Vec4D Matter_Overlap::SelectPositionForScatter(const double & b) const {
  double b1, b2, cosphi2;
  do {
    b1 = SelectB();
    b2 = SelectB();
    cosphi2 = (b1*b1-b2*b2-b*b)/(2.*b2*b);
  } while (cosphi2>1. || cosphi2<-1.);
  double sinphi2 = (ran->Get()>0.5?-1.:1.)*sqrt(1.-sqr(cosphi2));
  return Vec4D(0.,b/2.+b2*cosphi2,b2*sinphi2,0.);
}

ATOOLS::Vec4D Matter_Overlap::SelectRelativePositionForParton() const {
  double b   = SelectB();
  double phi = 2.*M_PI*ran->Get();
  return Vec4D(0.,b*cos(phi),b*sin(phi),0.);
}

double MO_Integrand::operator()(double b) {
  /////////////////////////////////////////////////////////////////////////////////
  // Integrand for d^2b O(b) = 2 pi b db O(b), where O(b) is the time-integrated
  // matter overlap, being the tricky part of the numerator in Eq.(32).
  // This does not include the prefactor k.
  /////////////////////////////////////////////////////////////////////////////////
  return 2.*M_PI*b*(*p_mo)(b);
}


