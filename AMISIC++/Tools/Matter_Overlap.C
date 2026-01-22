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
  m_norm(1./M_PI), m_invGeV2fm(rpa->hBar()*rpa->c()*1.e12),
  m_variation_index(0)
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

double Matter_Overlap::operator()(double b) {
  ///////////////////////////////////////////////////////////////////////////
  // Matter overlap in three forms available: Dynamic (i.e. kinematics
  // dependent single Gaussian), single_gaussian and double_gaussian
  ///////////////////////////////////////////////////////////////////////////
  return EvaluateAt(b, m_kradius);
}
// ue-reweighting
double Matter_Overlap::EvaluateAt(const double & b,const double & k) const {
  const double b2 = b*b;
  const double k2 = sqr(k);
  if (k2<=0.) return 0.;
  if (m_dynamic) {
    double effradius2 = k2 * m_dynradius2;
    if (effradius2<=0.) return 0.;
    return m_norm/effradius2 * exp(-b2/effradius2);
  }
  const size_t idx = m_variation_index;
  const std::array<double,4> & r2 = m_r2_variations[idx];
  const std::array<double,4> & rnorm = m_rnorm_variations[idx];
  double result = 0.;
  for (size_t i=0;i<4;i++) {
    if (rnorm[i]<=0.) continue;
    double denom = k2 * r2[i];
    if (denom<=0.) continue;
    result += rnorm[i] * exp(-b2/denom);
  }
  return m_norm/k2 * result;
}
// ue-reweighting

double Matter_Overlap::MaxValue(const double & b) {
  ///////////////////////////////////////////////////////////////////////////
  // Maximal value for dynamic matter overlap:
  // 1/[Pi*(R_{min,1}^2 + R_{min,2}^2)] * exp[-b^2/R_max^2]
  // with R_max^2 = R_{max,1}^2 + R_{max,2}^2.
  ///////////////////////////////////////////////////////////////////////////
  if (!m_dynamic) return (*this)(b);
  return ( m_norm/(sqr(m_kradius*m_fixradius))*
	   exp(-sqr(b/(m_kradius*m_maxradius)) ));
}

double Matter_Overlap::SelectB() const {
  ///////////////////////////////////////////////////////////////////////////
  // Algorithm:
  // 1. select a radius R according to matter content:
  // 2. Select b according to d^2b O(b) = d b^2 exp(-b^2/R^2).
  ///////////////////////////////////////////////////////////////////////////
  double effradius = m_kradius * m_radius[0], b=0.;
  if (m_dynamic) effradius = m_kradius * sqrt(m_dynradius2);
  else if (m_form[0]==matter_form::double_gaussian ||
	   m_form[1]==matter_form::double_gaussian) {
    double rand = ran->Get();
    for (int i=3;i>=0;i--) {
      rand -= m_fraction[i];
      if (rand<=1.e-6) { effradius = m_kradius * m_radius[i]; break; }
    }
  }
  do { b = sqrt(-log(Max(1.e-12,ran->Get())))*effradius; } while(b>=m_bmax);
  return b;
}

Vec4D Matter_Overlap::
SelectPositionForScatter(const double & B,
			 const double & x0, const double & Q20,
			 const double & x1, const double & Q21) const {
  ///////////////////////////////////////////////////////////////////////////
  // Independently select two impact paraemters b0 and b1 w.r.t.\ the incoming
  // beams, from their respective form factors until a combination that is
  // allowed given the overall impact paramter B is found.
  // Position is given in fm, calculation in 1/GeV.
  ///////////////////////////////////////////////////////////////////////////
  double b0, b1, cosphi2;
  size_t trials = 0;
  do {
    b0      = p_ffs[0]->B(x0,Q20);
    b1      = p_ffs[1]->B(x1,Q21);
    cosphi2 = (sqr(b0)-sqr(b1)-sqr(B))/(2.*b1*B);
  } while ( (cosphi2>1. || cosphi2<-1.) && (++trials)<10000);
  if (trials>=9999) { b1 = B/2.; cosphi2 = -1.; }
  double sinphi2 = (ran->Get()>0.5?-1.:1.)*sqrt(1.-sqr(cosphi2));
  return m_invGeV2fm*Vec4D(0.,B/2.+b1*cosphi2,b1*sinphi2,0.);
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

size_t Matter_Overlap::VariationSize() const {
  size_t n0 = p_ffs[0] ? p_ffs[0]->VariationSize() : size_t(1);
  size_t n1 = p_ffs[1] ? p_ffs[1]->VariationSize() : size_t(1);
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
  const size_t nvar = VariationSize();
  m_r2_variations.assign(nvar, std::array<double,4>{0.,0.,0.,0.});
  m_rnorm_variations.assign(nvar, std::array<double,4>{0.,0.,0.,0.});
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
      m_fraction[2*i+j] = ( (i==0 ? fraction[i] : 1.-fraction[i] ) *
			    (j==0 ? fraction[j] : 1.-fraction[j] ) );
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

  for (size_t ivar=0; ivar<nvar; ++ivar) {
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

void Matter_Overlap::Output(const double & check) {
  msg_Info()<<"   "<<std::string(77,'-')<<"\n"
	    <<"   | Matter_Overlap:"<<std::string(59,' ')<<"|\n"
	    <<"   | Integral up to b_max = "
	    <<std::setprecision(6)<<std::setw(8)<<(m_bmax*m_invGeV2fm)<<" fm yields "
	    <<std::setprecision(6)<<std::setw(8)<<check<<"."
	    <<std::string(23,' ')<<"|\n";
  for (size_t i=0;i<2;i++) {
    msg_Info()<<"   | "<<std::setw(20)<<m_form[i]<<", R_1 = "
	      <<std::setprecision(4)<<std::setw(6)
	      <<(p_ffs[0]->Radius1()*m_invGeV2fm)<<" fm";
    if (m_form[i]==matter_form::double_gaussian) {
      msg_Info()<<", f_1 = "<<std::setprecision(4)<<std::setw(6)
		<<p_ffs[i]->Fraction1()<<", "
		<<"R_2 = "<<std::setprecision(4)<<std::setw(6)
		<<(p_ffs[i]->Radius2()*m_invGeV2fm)<<" fm"
		<<std::string(6,' ')<<"|\n";
    }
    else msg_Info()<<std::string(37,' ')<<"|\n";
  }
  if (m_dynamic)
    msg_Info()<<"   | Maximal x-dependent radius: "
	      <<std::setprecision(4)<<std::setw(6)<<m_maxradius<<" 1/GeV = "
	      <<std::setprecision(4)<<std::setw(6)<<(m_maxradius*m_invGeV2fm)<<" fm. "
	      <<std::string(20,' ')<<"|\n";
  
  // ue-reweighting
  const size_t nvar = VariationSize();
  if (nvar > 1) {
    MO_Integrand moint(this);
    Gauss_Integrator integrator(&moint);
    for (size_t ivar=0; ivar<nvar; ++ivar) {
      SetVariationIndex(ivar);
      const double bmax_var = Bmax();
      double bmin = 0., bstep = m_bstep, previous, result = 0.;
      do {
        result  += previous = integrator.Integrate(bmin,bmin+bstep,1.e-8,1);
        bmin    += bstep;
      } while (dabs(previous/result)>1.e-10 && bmin < bmax_var);
      
      msg_Info()<<"   "<<std::string(77,'-')<<"\n"
                <<"   | Variation "<<ivar<<":"<<std::string(62,' ')<<"|\n"
                <<"   | Integral up to b_max = "
                <<std::setprecision(6)<<std::setw(8)<<(bmax_var*m_invGeV2fm)<<" fm yields "
                <<std::setprecision(6)<<std::setw(8)<<result<<"."
                <<std::string(23,' ')<<"|\n";
      for (size_t i=0;i<2;i++) {
        msg_Info()<<"   | "<<std::setw(20)<<m_form[i]<<", R_1 = "
                  <<std::setprecision(4)<<std::setw(6)
                  <<(p_ffs[i]->Radius1At(ivar)*m_invGeV2fm)<<" fm";
        if (m_form[i]==matter_form::double_gaussian) {
          msg_Info()<<", f_1 = "<<std::setprecision(4)<<std::setw(6)
                    <<p_ffs[i]->Fraction1At(ivar)<<", "
                    <<"R_2 = "<<std::setprecision(4)<<std::setw(6)
                    <<(p_ffs[i]->Radius2At(ivar)*m_invGeV2fm)<<" fm"
                    <<std::string(6,' ')<<"|\n";
        }
        else msg_Info()<<std::string(37,' ')<<"|\n";
      }
    }
    SetVariationIndex(0);
  }
  // ue-reweighting
  
  msg_Info()<<"   "<<std::string(77,'-')<<"\n\n";
}


