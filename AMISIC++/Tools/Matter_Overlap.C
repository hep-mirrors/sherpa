#include "AMISIC++/Tools/Matter_Overlap.H"
#include "AMISIC++/Tools/MI_Parameters.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

using namespace AMISIC;
using namespace REMNANTS;
using namespace ATOOLS;


/////////////////////////////////////////////////////////////////////////////
// All equations in this file refer to either
// Sjostrand-van der Zijl, PRD 36 (1987) 2019, denoted as (SZ, ...), or
// Corke-Sjostrand, JHEP 05 (2011) 009, denoted as (CS, ... )
/////////////////////////////////////////////////////////////////////////////


Matter_Overlap::Matter_Overlap() :
  ///////////////////////////////////////////////////////////////////////////
  // Our default are static form factors without any dependence on the parton
  // kinematics, equivalent to the SZ model.  This can be overwritten in the
  // "Initialize" method
  ///////////////////////////////////////////////////////////////////////////
  m_dynamic(false),
  ///////////////////////////////////////////////////////////////////////////
  // Parameters for a grid definition of the matter overlap, plus an integral
  // to fix them, and an axis of values for use if the overlap is dynamic
  ///////////////////////////////////////////////////////////////////////////
  m_bstep(0.), m_bmax(0.), m_integral(0.), p_bbins(nullptr),
  ///////////////////////////////////////////////////////////////////////////
  // A parameter to rescale the effective radius of the matter overlap for
  // dynamic hadron radii - it is fixed such that the hard (perturbative)
  // parton-parton cross section is reproduced when modified with the dynamic
  // hadron radii.
  ///////////////////////////////////////////////////////////////////////////
  m_kradius(1.),
  ///////////////////////////////////////////////////////////////////////////
  // The norm come from the (implicitly normalised) Gaussian-form matter
  // distributions of the hadron form factors ~pi^{-3/2} times the pi^2 from
  // the time-integrated overlap when they collide.
  ///////////////////////////////////////////////////////////////////////////
  m_norm(1./M_PI),
  ///////////////////////////////////////////////////////////////////////////
  m_matter_form_ivar(0)
{
  ///////////////////////////////////////////////////////////////////////////
  // This is a default initialization of default values.  Apart from the
  // actual radius it can serve for both the single_gaussian and the
  // x-dependent from factor forms.
  ///////////////////////////////////////////////////////////////////////////
  p_ffs[0]  = p_ffs[1]  = NULL;
  m_form[0] = m_form[1] = matter_form::single_gaussian;
  for (size_t i=0;i<4;i++) {
    m_radius[i]   = (i==0) ? 1. : 0.;
    m_radius2[i]  = sqr(m_radius[i]);
    m_rnorm[i]    = (i==0) ? 1./m_radius2[i] : 0.;
    m_fraction[i] = (i==0) ? 1. : 0.;
  }
  m_maxradius  = m_radius[0];
}

Matter_Overlap::~Matter_Overlap() { if (p_bbins) delete p_bbins; }

double Matter_Overlap::operator()(const double b) {
  ///////////////////////////////////////////////////////////////////////////
  // Matter overlap in three forms available: Dynamic (i.e. kinematics
  // dependent single Gaussian), single_gaussian and double_gaussian
  ///////////////////////////////////////////////////////////////////////////
  return (*this)(b, m_kradius);
}

double Matter_Overlap::operator()(const double b, const double k) {
  const double b2 = b*b;
  const double k2 = sqr(k);
  if (k2<=0.) return 0.;
  if (m_dynamic) {
    double effradius2 = k2 * m_dynradius2_variations[m_matter_form_ivar];
    if (effradius2<=0.) return 0.;
    return m_norm/effradius2 * exp(-b2/effradius2);
  }
  const std::array<double,4> & r2    = m_r2_variations[m_matter_form_ivar];
  const std::array<double,4> & rnorm = m_rnorm_variations[m_matter_form_ivar];
  double result = 0.;
  for (size_t i=0;i<4;i++) {
    if (rnorm[i]<=0.) continue;
    double denom = k2 * r2[i];
    if (denom<=0.) continue;
    result += rnorm[i] * exp(-b2/denom);
  }
  return m_norm/k2 * result;
}

double Matter_Overlap::MaxValue(const double & b) {
  ///////////////////////////////////////////////////////////////////////////
  // Maximal value for dynamic matter overlap:
  // 1/[Pi*(R_{min,1}^2 + R_{min,2}^2)] * exp[-b^2/R_max^2]
  // with R_max^2 = R_{max,1}^2 + R_{max,2}^2.
  ///////////////////////////////////////////////////////////////////////////
  return MaxValue(b, m_kradius);
}

double Matter_Overlap::MaxValue(const double & b, double k) {
  if (!m_dynamic) return (*this)(b, k);
  return ( m_norm/(sqr(k*m_fixradius_variations[m_matter_form_ivar]))*
	   exp(-sqr(b/(k*m_maxradius_variations[m_matter_form_ivar])) ));
}

double Matter_Overlap::SelectB() const {
  ///////////////////////////////////////////////////////////////////////////
  // Select an impact parameter b in fm from the matter overlap distribution.
  //
  // Step 1 — Choose an effective Gaussian radius for this event:
  //
  //   Dynamic case: use the kinematics-dependent radius sqrt(m_dynradius2).
  //
  //   Double-Gaussian case: the overlap of two double-Gaussian form factors
  //   decomposes into four Gaussian components (indexed 0..3), each with a
  //   combined radius m_radius[i] and a fractional weight m_fraction[i]
  //   (the four weights sum to 1).  We pick one component with probability
  //   proportional to its weight via inverse-CDF sampling: subtract weights
  //   from a uniform random variable until it crosses zero.
  //   If no component is selected due to floating-point rounding, the first
  //   component (i = 0) is used as a fallback.
  //
  //   Single-Gaussian case: only one component exists, so m_radius[0] is
  //   used directly (no random selection needed).
  //
  // Step 2 — Sample b from the 2D Gaussian with the chosen radius R:
  //
  //   d^2b O(b) ~ b db exp(-b^2/R^2)
  //
  //   Setting x = b^2/R^2 gives an exponential distribution in x, which is
  //   sampled by inversion: x = -log(ran), so b = R*sqrt(-log(ran)).
  ///////////////////////////////////////////////////////////////////////////

  // Step 1
  double effradius = m_kradius * m_radius[0];

  if (m_dynamic) {
    effradius = m_kradius * sqrt(m_dynradius2);
  }
  else if (m_form[0] == matter_form::double_gaussian ||
           m_form[1] == matter_form::double_gaussian) {
    double remaining = ran->Get();
    for (int i = 3; i >= 0; i--) {
      remaining -= m_fraction[i];
      if (remaining <= 0.) {
        effradius = m_kradius * m_radius[i];
        break;
      }
    }
  }

  // Step 2
  double b;
  do {
    b = sqrt(-log(Max(1.e-12, ran->Get()))) * effradius;
  } while (b >= m_bmax);
  return b;
}

bool Matter_Overlap::
        SelectPositionForScatter(double B,
                                 const double & x0, const double & Q20,
                                 const double & x1, const double & Q21,
                                 Vec4D & pos) const {
  ///////////////////////////////////////////////////////////////////////////
  // Independently select two impact parameters b0 and b1 w.r.t. the incoming
  // beams, from their respective Q^2-dependent form factors until a combination
  // that is allowed given the overall impact parameter B is found.
  //
  // Coordinate system: Origin at midpoint between beams
  //   - Beam 0 at (+B/2, 0)
  //   - Beam 1 at (-B/2, 0)
  //   - B is the full transverse distance between beam centers
  //
  // Position is returned in millimeters (converted from fm).
  ///////////////////////////////////////////////////////////////////////////
  // converting the impact parameter B from the internal unit fm to millimeters
  B *= 1.e-12;
  double       b0, b1, cosphi1, sinphi1;
  size_t       trials     = 0;
  const size_t max_trials = 10000;

  // Sample b0 and b1 from Q^2-dependent form factors and reject if any of the
  // three triangle inequalities is violated: B < b0+b1, b0 < B+b1, b1 < B+b0.
  do {
    b0 = p_ffs[0]->B(x0, Q20);
    b1 = p_ffs[1]->B(x1, Q21);
  } while ((B >= b0 + b1 || b0 >= B + b1 || b1 >= B + b0) &&
           (++trials) < max_trials);

  if (trials >= max_trials) return false;

  cosphi1 = (sqr(B) + sqr(b1) - sqr(b0)) / (2. * B * b1);
  cosphi1 = std::max(-1., std::min(1., cosphi1));

  sinphi1 = (ran->Get() > 0.5 ? -1. : 1.) * sqrt(1. - sqr(cosphi1));

  // Position relative to beam 1 at (-B/2, 0)
  pos = Vec4D(0., -B / 2. + b1 * cosphi1, b1 * sinphi1, 0.);

  return true;
}

void Matter_Overlap::Initialize(Remnant_Handler * const rh,
				PDF::ISR_Handler * const isr) {
  for (size_t i=0;i<2;i++) {
    p_ffs[i]  = rh->GetFormFactor(i);
    m_form[i] = p_ffs[i]->GetForm();
  }
  InitializeFFParams(isr);
  Output(CalculateIntegral());
  size_t nbins = size_t((*mipars)["nB_bins"]);
  double bmin  = 0.00001*m_radius[0];
  p_bbins = new axis(nbins, bmin, m_bmax, axis_mode::log);
}

size_t Matter_Overlap::MatterFormVariationSize() const {
  size_t n0 = p_ffs[0] ? p_ffs[0]->MatterFormVariationSize() : size_t(1);
  size_t n1 = p_ffs[1] ? p_ffs[1]->MatterFormVariationSize() : size_t(1);
  return std::max(n0, n1);
}

void Matter_Overlap::InitializeFFParams(PDF::ISR_Handler * const isr) {
  ///////////////////////////////////////////////////////////////////////////
  //  Initialize form factor parameter with a suitable integrand for the
  //  check for proper normalisation in CalculateIntegral.
  ///////////////////////////////////////////////////////////////////////////
  MO_Integrand * moint = NULL;
  if ((m_form[0]==matter_form::single_gaussian ||
       m_form[0]==matter_form::double_gaussian) &&
      (m_form[1]==matter_form::single_gaussian ||
       m_form[1]==matter_form::double_gaussian) ) {
    InitializeStaticFFParams();
  }
  else if (m_form[0]==matter_form::x_dependent_gaussian &&
	   m_form[1]==matter_form::x_dependent_gaussian) {
    InitializeDynamicFFParams(isr);
  }
  else {
    msg_Error()<<"Error in "<<METHOD<<": "
	       <<"combination of matter forms not implemented yet.\n"
	       <<"   Will continue with static formfactors.\n";
    THROW(fatal_error,"Inconsistent settings in AMISIC");
  }
}

void Matter_Overlap::InitializeStaticFFParams() {
  ///////////////////////////////////////////////////////////////////////////
  // Initialise matter overlap from the form factors:
  // could be single- or double Gaussians, e.g. from Eq (SZ 19)
  ///////////////////////////////////////////////////////////////////////////
  m_dynamic = false;
  const size_t n_matter_form_variations = MatterFormVariationSize();
  m_r2_variations.assign(n_matter_form_variations, std::array<double,4>{0.,0.,0.,0.});
  m_rnorm_variations.assign(n_matter_form_variations, std::array<double,4>{0.,0.,0.,0.});
  double fraction[2], radius[2][2];
  for (size_t i=0;i<2;i++) {
    fraction[i] = ( (m_form[i]==matter_form::single_gaussian) ?
		    1 : p_ffs[i]->Fraction1() );
    for (size_t j=0;j<2;j++)
      radius[i][j] = ( j==0 ? p_ffs[i]->Radius1() : p_ffs[i]->Radius2() );
  }
  double minR  = 1.e6;
  m_fixradius = 0.;
  for (size_t i=0;i<2;i++) {
    for (size_t j=0;j<2;j++) {
      m_fraction[2*i+j] = ( (i==0 ? fraction[0] : 1.-fraction[0] ) *
			    (j==0 ? fraction[1] : 1.-fraction[1] ) );
      m_radius2[2*i+j]  = sqr(radius[0][i]) + sqr(radius[1][j]);
      m_radius[2*i+j]   = sqrt(m_radius2[2*i+j]);
      m_rnorm[2*i+j]    = ( radius[0][i] > 0. && radius[1][j]>0. ?
			    m_fraction[2*i+j]/m_radius2[2*i+j] : 0. );
      if (m_rnorm[2*i+j]<=0.) continue;
      if (m_radius[2*i+j] < minR && m_radius[2*i+j]>0.) minR = m_radius[2*i+j];
      if (m_fixradius < m_radius[2*i+j]) m_fixradius = m_radius[2*i+j];
    }
  }
  m_dynradius2 = m_radius2[0];
  m_bstep = minR/100.;

  for (size_t ivar=0; ivar<n_matter_form_variations; ++ivar) {
    double frac[2];
    double rad[2][2];
    frac[0] = (p_ffs[0]->GetForm()==matter_form::double_gaussian) ?
      p_ffs[0]->Fraction1At(ivar) : 1.0;
    frac[1] = (p_ffs[1]->GetForm()==matter_form::double_gaussian) ?
      p_ffs[1]->Fraction1At(ivar) : 1.0;
    rad[0][0] = p_ffs[0]->Radius1At(ivar);
    rad[0][1] = p_ffs[0]->Radius2At(ivar);
    rad[1][0] = p_ffs[1]->Radius1At(ivar);
    rad[1][1] = p_ffs[1]->Radius2At(ivar);

    for (size_t i=0;i<2;i++) {
      for (size_t j=0;j<2;j++) {
        const size_t idx = 2*i+j;
        const double frac_pair = ( (i==0 ? frac[0] : 1.-frac[0]) *
                                   (j==0 ? frac[1] : 1.-frac[1]) );
        const double r2 = sqr(rad[0][i]) + sqr(rad[1][j]);
        const double rnorm = (rad[0][i]>0. && rad[1][j]>0.) ? frac_pair/r2 : 0.;
        m_r2_variations[ivar][idx] = r2;
        m_rnorm_variations[ivar][idx] = rnorm;
      }
    }
  }
}

void Matter_Overlap::InitializeDynamicFFParams(PDF::ISR_Handler * const isr) {
  ///////////////////////////////////////////////////////////////////////////
  // Initialise matter overlap from the form factors:
  // this may involving obtaining the maximal radius of the parton cloud based
  // on the parton dynamics.
  // For the time being we only assume the x-dependent form factors with
  // a radius that increases with smaller x.
  ///////////////////////////////////////////////////////////////////////////
  m_dynamic    = true;
  m_radius2[0] = sqr(p_ffs[0]->Radius1())+sqr(p_ffs[1]->Radius1());
  m_radius[0]  = m_fixradius = sqrt(m_radius2[0]);
  m_rnorm[0]   = 1./m_radius2[0];
  m_bstep      = Min(1.,m_radius[0])/100.;
  for (size_t i=0;i<2;i++) {
    m_xmin[i] = isr->XMin(0); m_xmax[i] = isr->XMax(0);
  }
  FixDynamicRadius(m_xmin[0],m_xmin[1]);
  m_maxradius = sqrt(m_dynradius2);

  const size_t n_matter_form_variations = MatterFormVariationSize();
  m_dynradius2_variations.assign(n_matter_form_variations, 0.);
  m_fixradius_variations.assign(n_matter_form_variations, 0.);
  m_maxradius_variations.assign(n_matter_form_variations, 0.);
  if (n_matter_form_variations <= 1) return;
  for (size_t ivar=0; ivar<n_matter_form_variations; ++ivar) {
    const double R0_1 = p_ffs[0]->Radius1At(ivar);
    const double R0_2 = p_ffs[1]->Radius1At(ivar);
    m_fixradius_variations[ivar] = sqrt(sqr(R0_1) + sqr(R0_2));
    const double Rmax_1 = p_ffs[0]->Radius(m_xmin[0], -1., ivar);
    const double Rmax_2 = p_ffs[1]->Radius(m_xmin[1], -1., ivar);
    m_maxradius_variations[ivar] = sqrt(sqr(Rmax_1) + sqr(Rmax_2));
  }
}

void Matter_Overlap::FixDynamicRadius(const double & x1,  const double & x2,
				      const double & Q21, const double & Q22) {
  ///////////////////////////////////////////////////////////////////////////
  // For dynamic form factors we assume that the (Gaussian) radius depends on
  // x and Q2 of the partons, we combine them for the matter profile in the
  // spirit of Eq. (CS, 11).
  ///////////////////////////////////////////////////////////////////////////
  if (!m_dynamic) return;
  m_dynradius2 = (sqr(p_ffs[0]->Radius(Max(m_xmin[0],Min(m_xmax[0],x1)),Q21))+
                  sqr(p_ffs[1]->Radius(Max(m_xmin[1],Min(m_xmax[1],x2)),Q22)) );

  const size_t n_matter_form_variations = m_dynradius2_variations.size();
  if (n_matter_form_variations <= 1) return;
  const double x1c = Max(m_xmin[0],Min(m_xmax[0],x1));
  const double x2c = Max(m_xmin[1],Min(m_xmax[1],x2));
  for (size_t ivar=0; ivar<n_matter_form_variations; ++ivar) {
    m_dynradius2_variations[ivar] =
      ( sqr(p_ffs[0]->Radius(x1c, Q21, ivar)) +
        sqr(p_ffs[1]->Radius(x2c, Q22, ivar)) );
  }
}


double Matter_Overlap::CalculateIntegral() {
  ///////////////////////////////////////////////////////////////////////////
  // Integral int d^2b O(b), numerator Eq.(SZ, 32)
  ///////////////////////////////////////////////////////////////////////////
  MO_Integrand moint(this);
  Gauss_Integrator integrator(&moint);
  double bmin = 0., bstep = m_bstep, previous, result = 0.;
  do {
    result  += previous = integrator.Integrate(bmin,bmin+bstep,1.e-8,1);
    bmin    += bstep;
  } while (dabs(previous/result)>1.e-10);
  m_bmax     = bmin+bstep;
  m_integral = result;
  return m_integral;
}

double MO_Integrand::operator()(double b) {
  ///////////////////////////////////////////////////////////////////////////
  // Integrand for d^2b O(b) = 2 pi b db O(b), where O(b) is the
  // time-integrated matter overlap, part of the numerator in Eq.(SZ, 32).
  // This should integrate to 1.
  ///////////////////////////////////////////////////////////////////////////
  return 2.*M_PI*b*(*p_mo)(b);
}

double Matter_Overlap::ComputeBmaxForVariation(size_t ivar) {
  ///////////////////////////////////////////////////////////////////////////
  // Compute bmax for a specific matter form variation using the same
  // algorithm as CalculateIntegral, but without modifying the permanent
  // m_bmax member. Temporarily switches m_kradius, then restores it.
  ///////////////////////////////////////////////////////////////////////////
  const double saved_kradius = m_kradius;
  SetMatterFormVariationIndex(ivar);
  SetKRadius(1.0);
  double bstep;
  if (m_dynamic) {
    FixDynamicRadius(m_xmin[0], m_xmin[1]);
    bstep = Min(1., m_fixradius_variations[m_matter_form_ivar]) / 100.;
  } else {
    double minR = 1.e6;
    for (size_t i = 0; i < 4; i++) {
      if (m_rnorm_variations[m_matter_form_ivar][i] <= 0.) continue;
      double r = sqrt(m_r2_variations[m_matter_form_ivar][i]);
      if (r > 0. && r < minR) minR = r;
    }
    bstep = minR / 100.;
  }
  MO_Integrand moint(this);
  Gauss_Integrator integrator(&moint);
  double bmin = 0., previous, result = 0.;
  do {
    result += previous = integrator.Integrate(bmin, bmin+bstep, 1.e-8, 1);
    bmin   += bstep;
  } while (dabs(previous/result)>1.e-10);
  double bmax_var = bmin + bstep;

  SetMatterFormVariationIndex(0);
  SetKRadius(saved_kradius);
  return bmax_var;
}

double Matter_Overlap::ComputeBminForVariation(size_t ivar) {
  ///////////////////////////////////////////////////////////////////////////
  // Compute bmin for a specific matter form variation, analogous to how
  // the nominal bmin = 0.00001 * m_radius[0] is computed in Initialize.
  ///////////////////////////////////////////////////////////////////////////
  SetMatterFormVariationIndex(ivar);
  double bmin = 0.;
  if (m_dynamic) {
    const double R0_nom = p_ffs[0]->Radius1();
    const double R0_var = p_ffs[0]->Radius1At(m_matter_form_ivar);
    const double scale  = (R0_nom > 0.) ? R0_var / R0_nom : 1.;
    bmin = 0.00001 * sqrt(m_radius2[0]) * scale;
  } else {
    bmin = 0.00001 * sqrt(m_r2_variations[m_matter_form_ivar][0]);
  }
  SetMatterFormVariationIndex(0);
  return bmin;
}

void Matter_Overlap::Output(const double & check) {
  msg_Info()<<"   "<<std::string(77,'-')<<"\n"
	    <<"   | Matter_Overlap:"<<std::string(59,' ')<<"|\n"
	    <<"   | Integral up to b_max = "
	    <<std::setprecision(6)<<std::setw(8)<<m_bmax<<" fm yields "
	    <<std::setprecision(6)<<std::setw(8)<<check<<"."
	    <<std::string(23,' ')<<"|\n";
  for (size_t i=0;i<2;i++) {
    msg_Info()<<"   | "<<std::setw(20)<<m_form[i]<<", R_1 = "
	      <<std::setprecision(4)<<std::setw(6)
	      <<p_ffs[i]->Radius1()<<" fm";
    if (m_form[i]==matter_form::double_gaussian) {
      msg_Info()<<", f_1 = "<<std::setprecision(4)<<std::setw(6)
		<<p_ffs[i]->Fraction1()<<", "
		<<"R_2 = "<<std::setprecision(4)<<std::setw(6)
		<<p_ffs[i]->Radius2()<<" fm"
		<<std::string(6,' ')<<"|\n";
    }
    else if (m_form[i]==matter_form::x_dependent_gaussian) {
      msg_Info()<<", alpha = "<<std::setprecision(4)<<std::setw(6)
                <<p_ffs[i]->SoftExponent()
                <<std::string(21,' ')<<"|\n";
    }
    else msg_Info()<<std::string(37,' ')<<"|\n";
  }
  if (m_dynamic)
    msg_Info()<<"   | Maximal x-dependent radius: "
	            <<std::setprecision(4)<<std::setw(6)<<m_maxradius<<" fm. "
	            <<std::string(35,' ')<<"|\n";
  for (size_t ivar=1; ivar<MatterFormVariationSize(); ++ivar) Output(ivar);
  msg_Info()<<"   "<<std::string(77,'-')<<"\n\n";
}

void Matter_Overlap::Output(size_t ivar) {
  SetMatterFormVariationIndex(ivar);
  MO_Integrand moint(this);
  Gauss_Integrator integrator(&moint);
  double bmin = 0., bstep = m_bstep, previous, result = 0.;
  do {
    result += previous = integrator.Integrate(bmin, bmin+bstep, 1.e-8, 1);
    bmin   += bstep;
  } while (dabs(previous/result)>1.e-10 && bmin<m_bmax);

  msg_Info()<<"   "<<std::string(77,'-')<<"\n"
            <<"   |          v"<<std::setw(4)<<ivar<<":"<<std::string(59,' ')<<"|\n"
            <<"   | Integral up to b_max = "
            <<std::setprecision(6)<<std::setw(8)<<(m_bmax)<<" fm yields "
            <<std::setprecision(6)<<std::setw(8)<<result<<"."
            <<std::string(23,' ')<<"|\n";
  for (size_t i=0;i<2;i++) {
    msg_Info()<<"   | "<<std::setw(20)<<m_form[i]<<", R_1 = "
              <<std::setprecision(4)<<std::setw(6)
              <<p_ffs[i]->Radius1At(ivar)<<" fm";
    if (m_form[i]==matter_form::double_gaussian) {
      msg_Info()<<", f_1 = "<<std::setprecision(4)<<std::setw(6)
                <<p_ffs[i]->Fraction1At(ivar)<<", "
                <<"R_2 = "<<std::setprecision(4)<<std::setw(6)
                <<p_ffs[i]->Radius2At(ivar)<<" fm"
                <<std::string(6,' ')<<"|\n";
    }
    else if (m_form[i]==matter_form::x_dependent_gaussian) {
      msg_Info()<<", alpha = "<<std::setprecision(4)<<std::setw(6)
                <<p_ffs[i]->SoftExponentAt(ivar)
                <<std::string(21,' ')<<"|\n";
    }
    else msg_Info()<<std::string(37,' ')<<"|\n";
  }
  SetMatterFormVariationIndex(0);
  if (m_dynamic)
    msg_Info()<<"   | Maximal x-dependent radius: "
	            <<std::setprecision(4)<<std::setw(6)<<m_maxradius_variations[ivar]<<" fm. "
              <<std::string(35,' ')<<"|\n";
}


