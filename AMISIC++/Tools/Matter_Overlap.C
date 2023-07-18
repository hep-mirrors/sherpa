#include "AMISIC++/Tools/Matter_Overlap.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

using namespace AMISIC;
using namespace ATOOLS;


// All equations in this file refer to
// Sjostrand-van der Zijl, PRD 36 (1987) 2019.

Matter_Overlap::Matter_Overlap() :
  // The norm come from the (implicitly normalised) Gaussian matter distributions
  // of the hadron form factors ~pi^{-3/2} times the pi^2 from the time-integrated
  // overlap when they collide.
  m_bstep(0.), m_bmax(0.), m_integral(0.), m_norm(1./M_PI)
{
  for (size_t i=0;i<4;i++) {
    m_radius[i]   = (i==0) ? 1. : 0.;
    m_radius2[i]  = sqr(m_radius[i]);
    m_rnorm[i]    = (i==0) ? 1./m_radius2[i] : 0.;
    m_fraction[i] = (i==0) ? 1. : 0.;
  }
}

Matter_Overlap::~Matter_Overlap() {}

void Matter_Overlap::Initialize(REMNANTS::Remnant_Handler * remnant_handler) {
  InitializeFormFactors(remnant_handler);
  CalculateIntegral();
}

double Matter_Overlap::operator()(double b) {
  // Matter overlap in two forms available: Single_Gaussian and Double_Gaussian
  double result = 0., b2 = b*b;
  for (size_t i=0;i<4;i++) result += (m_rnorm[i]>0) ? m_rnorm[i] * exp(-b2/m_radius2[i]) : 0.;
  return m_norm * result;
}

double Matter_Overlap::SelectB(const bool & mode) const {
  // Algorithm:
  // 1. select a radius R according to matter content:
  // 2. Select b according to d^2b O(b) = d b^2 exp(-b^2/R^2).
  double b, radius, rand = ran->Get();
  for (size_t i=3;i>=0;i--) {
    rand -= m_fraction[i];
    if (rand<=0.) { radius = m_radius[i]; break; }
  }
  // b from Matter_Overlap, hence r^2_overlap = 2*r^2_formfactor (?)
  if (mode) radius *= sqrt(2.);
  do {
    b = sqrt(-log(Max(1.e-12,ran->Get())))*radius;
  } while (b>m_bmax);
  return b;
}


void Matter_Overlap::InitializeFormFactors(REMNANTS::Remnant_Handler * remnant_handler) {
  /////////////////////////////////////////////////////////////////////////////////
  // Initialise matter overlap from the form factors: 
  // could be single- or double Gaussians
  /////////////////////////////////////////////////////////////////////////////////
  double fraction[2], radius[2][2];
  for (size_t i=0;i<2;i++) {
    fraction[i] = remnant_handler->GetRemnant(i)->GetFF()->Fraction1();
    for (size_t j=0;j<2;j++) {
      radius[i][j] = ( j==0 ?
		       remnant_handler->GetRemnant(i)->GetFF()->Radius1() :
		       remnant_handler->GetRemnant(i)->GetFF()->Radius2() );
    }
  }
  double minR = 1.;
  for (size_t i=0;i<2;i++) {
    for (size_t j=0;j<2;j++) {
      m_fraction[2*i+j] = ( (i==0 ? fraction[i] : 1.-fraction[i] ) *
			    (j==0 ? fraction[j] : 1.-fraction[j] ) );
      m_radius2[2*i+j]  = sqr(radius[0][i]) + sqr(radius[1][j]);
      m_rnorm[2*i+j]    = ( radius[0][i] > 0. && radius[1][j]>0. ? m_fraction[2*i+j]/m_radius2[2*i+j] : 0. );
      m_radius[2*i+j]   = sqrt(m_radius2[2*i+j]);
      if (m_radius[2*i+j] < minR && m_radius[2*i+j]>0.) minR = m_radius[2*i+j];
    }
  }
  m_bstep = minR/100.;
}
  
void Matter_Overlap::CalculateIntegral() {
  // Integral int d^2b O(b), numerator Eq.(32)
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
  // Integrand for d^2b O(b) = 2 pi b db O(b), where O(b) is the time-integrated
  // matter overlap, being the tricky part of the numerator in Eq.(32).
  // This does not include the prefactor k.
  return 2.*M_PI*b*(*p_mo)(b);
}


