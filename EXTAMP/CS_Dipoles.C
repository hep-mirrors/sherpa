#include "EXTAMP/CS_Dipoles.H"
#include "ATOOLS/Org/Exception.H"

#include "assert.h"

using namespace EXTAMP;

///////////////////////////////////////////////////////////////
////////// CalcKinematics METHODS /////////////////////////////
///////////////////////////////////////////////////////////////

void FF_Dipole::CalcKinematics(const ATOOLS::Vec4D_Vector& p)
{
  /* Implementation of hep-ph/9605323v3 eq. (5.3) - (5.6) */
  
  assert(K()>1); assert(Emitter()>1); assert(Emitted()>1);

  const ATOOLS::Vec4D& pi = p[I()];
  const ATOOLS::Vec4D& pj = p[J()];
  const ATOOLS::Vec4D& pk = p[K()];
  
  m_kin.m_y         = pi*pj/(pi*pj+pj*pk+pk*pi);
  m_kin.m_zi        = pi*pk/(pj*pk+pi*pk);
  m_kin.m_zj        = 1.0-m_kin.m_zi;
  m_kin.m_pk_tilde  = 1.0/(1.0-m_kin.m_y)*pk;
  m_kin.m_pij_tilde = pi+pj-m_kin.m_y/(1.0-m_kin.m_y)*pk;
  m_kin.m_pi        = pi;
  m_kin.m_pj        = pj;
  m_kin.m_pk        = pk;
  
  /* Replace emitter momentum with combined momentum of (ij) and
     remove emitted. */
  m_kin.m_born_mom = p;
  m_kin.m_born_mom[Emitter()] = m_kin.m_pij_tilde;
  m_kin.m_born_mom[K()]       = m_kin.m_pk_tilde;
  m_kin.m_born_mom.erase(m_kin.m_born_mom.begin()+Emitted());

}

void FI_Dipole::CalcKinematics(const ATOOLS::Vec4D_Vector& p)
{
  /* Implementation of hep-ph/9605323v3 (5.37) - (5.42) with a=k */

  assert(K()<2); assert(I()>1); assert(J()>1);
  
  const ATOOLS::Vec4D& pi = p[I()];
  const ATOOLS::Vec4D& pj = p[J()];
  const ATOOLS::Vec4D& pa = p[K()];
  
  m_kin.m_x         = (pi*pa + pj*pa - pi*pj)/((pi+pj)*pa);
  m_kin.m_zi        = pi*pa/(pi*pa+pj*pa);
  m_kin.m_zj        = pj*pa/(pi*pa+pj*pa);
  m_kin.m_pa_tilde  = m_kin.m_x*pa;
  m_kin.m_pij_tilde = pi+pj-(1.0-m_kin.m_x)*pa;
  m_kin.m_pi        = pi;
  m_kin.m_pj        = pj;
  m_kin.m_pa        = pa;

  /* Replace emitter momentum with combined momentum of (ij) and
     remove emitted. */
  m_kin.m_born_mom = p;
  m_kin.m_born_mom[Emitter()] = m_kin.m_pij_tilde;
  m_kin.m_born_mom[K()]       = m_kin.m_pa_tilde;
  m_kin.m_born_mom.erase(m_kin.m_born_mom.begin()+Emitted());
  
}

void IF_Dipole::CalcKinematics(const ATOOLS::Vec4D_Vector& p)
{
  /* Implementation of hep-ph/9605323v3 (5.62) - (5.64) */

  assert(K()>1); assert(Emitter()<2); assert(Emitted()>1);

  const ATOOLS::Vec4D& pa = p[Emitter()];
  const ATOOLS::Vec4D& pi = p[Emitted()];
  const ATOOLS::Vec4D& pk = p[K()];

  m_kin.m_x         = (pk*pa + pi*pa - pi*pk)/((pk+pi)*pa);
  m_kin.m_ui        = pi*pa/(pi*pa+pk*pa);

  m_kin.m_pk_tilde  = pk+pi-(1.0-m_kin.m_x)*pa;
  m_kin.m_pai_tilde = m_kin.m_x*pa;
  m_kin.m_pa        = pa;
  m_kin.m_pi        = pi;
  m_kin.m_pk        = pk;

  /* Replace emitter momentum with combined momentum of (ij) and
     remove emitted. */
  m_kin.m_born_mom = p;
  m_kin.m_born_mom[Emitter()] = m_kin.m_pai_tilde;
  m_kin.m_born_mom[K()]       = m_kin.m_pk_tilde;
  m_kin.m_born_mom.erase(m_kin.m_born_mom.begin()+Emitted());

}

void II_Dipole::CalcKinematics(const ATOOLS::Vec4D_Vector& p)
{
  /* Implementation of hep-ph/9605323v3 (5.137) - (5.140) */

  const ATOOLS::Vec4D& pa = p[Emitter()];
  const ATOOLS::Vec4D& pi = p[Emitted()];
  const ATOOLS::Vec4D& pb = p[K()];
  
  m_kin.m_x         = (pa*pb-pi*pa-pi*pb)/(pa*pb);
  m_kin.m_v         = (pa*pi)/(pa*pb);
  m_kin.m_pb_tilde  = pb;
  m_kin.m_pai_tilde = m_kin.m_x*pa;
  m_kin.m_pa        = pa;
  m_kin.m_pi        = pi;
  m_kin.m_pb        = pb;
  m_kin.m_born_mom  = p;

  /* Apply transformation (5.139) */
  ATOOLS::Vec4D Ka      = pa+pb-pi;
  ATOOLS::Vec4D Katilde = m_kin.m_pai_tilde + pb;
  for(size_t n(0); n<p.size(); n++)
    m_kin.m_born_mom[n] = p[n]-2.0*p[n]*(Ka+Katilde)/(Ka+Katilde).Abs2()*(Ka+Katilde)+2.0*(p[n]*Ka)/Ka.Abs2()*Katilde;
  
  /* Replace emitter momentum with combined momentum of (ij) and
     remove emitted. */
  m_kin.m_born_mom[Emitter()] = m_kin.m_pai_tilde;
  m_kin.m_born_mom[K()]       = m_kin.m_pb_tilde;
  m_kin.m_born_mom.erase(m_kin.m_born_mom.begin()+Emitted());

}

///////////////////////////////////////////////////////////////
////////// CalcKinDependentPrefac METHODS /////////////////////
///////////////////////////////////////////////////////////////

double FF_Dipole::CalcKinDependentPrefac() const
{
  const ATOOLS::Vec4D& pi = m_kin.m_pi;
  const ATOOLS::Vec4D& pj = m_kin.m_pj;
  
  /* hep-ph/9605323v3 eq. (5.2) */
  return -1.0/(2.0*pi*pj);
}

double FI_Dipole::CalcKinDependentPrefac() const
{
  const ATOOLS::Vec4D& pi = m_kin.m_pi;
  const ATOOLS::Vec4D& pj = m_kin.m_pj;
  const double& x = m_kin.m_x;
  
  /* hep-ph/9605323v3 eq. (5.36) */
  return -1.0/(2.0*pi*pj*x);
}

double IF_Dipole::CalcKinDependentPrefac() const
{
  const ATOOLS::Vec4D& pi = m_kin.m_pi;
  const ATOOLS::Vec4D& pa = m_kin.m_pa;
  const double& x = m_kin.m_x;
  
  /* hep-ph/9605323v3 eq. (5.61) */
  return -1.0/(2.0*pi*pa*x);
}

double II_Dipole::CalcKinDependentPrefac() const
{
  const ATOOLS::Vec4D& pi = m_kin.m_pi;
  const ATOOLS::Vec4D& pa = m_kin.m_pa;
  const double& x = m_kin.m_x;
  
  /* hep-ph/9605323v3 eq. (5.136) */
  return -1.0/(2.0*pi*pa*x);
}

///////////////////////////////////////////////////////////////
////////// CalcA METHODS  /////////////////////////////////////
///////////////////////////////////////////////////////////////

double FF_Dipole::CalcA() const
{
  double zi = m_kin.m_zi;
  double zj = m_kin.m_zj;
  const double& y(m_kin.m_y);
  
  /* q->qg expression depends on the flavour assignment being
     i=quark, j=gluon. Need to respect that here by swapping z_i,z_j
     if neccessary. Does not affect other splittings, so can be done
     for all cases. */
  if(FlavI().IsGluon()) std::swap(zi,zj);
  
  /* Coefficients of \delta_{ss^\prime} or -g^{\mu\nu} in eq. (5.7)
     - (5.9). */
  if(FlavType()==FlavourType::qtoqg)
    return 2.0/(1.0-zi*(1.0-y)) - (1.+ zi);
  if(FlavType()==FlavourType::gtoqq)
    return 1.0;
  if(FlavType()==FlavourType::gtogg)
    return 1.0/(1.0-zi*(1.0-y)) + 1.0/(1-zj*(1.0-y)) - 2.0;
  
  THROW(fatal_error, "Internal error");
}

double FI_Dipole::CalcA() const
{
  double zi = m_kin.m_zi;
  double zj = m_kin.m_zj;
  const double& x(m_kin.m_x);
  
  /* q->qg expression depends on the flavour assignment being
     i=quark, j=gluon. Need to respect that here by swapping z_i,z_j
     if neccessary. Does not affect other splittings, so can be done
     for all cases. */
  if(FlavI().IsGluon()) std::swap(zi,zj);
  
  /* Coefficients of \delta_{ss^\prime} or -g^{\mu\nu} in eq. (5.39)
     - (5.41). */
  if(FlavType()==FlavourType::qtoqg)
    return 2.0/(1.0-zi+(1.0-x)) - (1.+ zi);
  if(FlavType()==FlavourType::gtoqq)
    return 1.0;
  if(FlavType()==FlavourType::gtogg)
    return 1.0/(1.0-zi+(1.0-x)) + 1.0/(1-zj+(1.0-x)) - 2.0;
  
  THROW(fatal_error, "Internal error");
}

double IF_Dipole::CalcA() const
{
  const double& x  = (m_kin.m_x);
  const double& ui = (m_kin.m_ui);
  
  /* Need this to distinguish (5.65) from (5.66) */
  const ATOOLS::Flavour& flav_a = RealFlavours()[std::min(I(),J())];
  
  /* Coefficients of \delta_{ss^\prime} or -g^{\mu\nu} in eq. (5.65)
     - (5.68). */
  if((FlavType()==FlavourType::qtoqg) && flav_a.IsQuark())
    return 2.0/(1.0-x+ui) - (1.+x);
  if((FlavType()==FlavourType::qtoqg) && flav_a.IsGluon())
    return 1.0-2.0*x*(1.0-x);
  if(FlavType()==FlavourType::gtoqq)
    return x;
  if(FlavType()==FlavourType::gtogg)
    return 1/(1.0-x+ui)-1.0+x*(1.0-x);
    
  THROW(fatal_error, "Internal error");
}

double II_Dipole::CalcA() const
{
  const double& x  = (m_kin.m_x);
  const double& z  = (SubtractionType()==1) ? m_kin.m_x+m_kin.m_v : x;
  
  /* Need this to distinguish (5.145) from (5.147) */
  const ATOOLS::Flavour& flav_a = RealFlavours()[std::min(I(),J())];
  
  /* Coefficients of \delta_{ss^\prime} or -g^{\mu\nu} in eq. (5.145)
     - (5.148). */
  if((FlavType()==FlavourType::qtoqg) && flav_a.IsQuark())
    return 2.0/(1.0-x) - (1.+z);
  
  if((FlavType()==FlavourType::qtoqg) && flav_a.IsGluon())
    return 1.0-2.0*z*(1.0-z);
  
  if(FlavType()==FlavourType::gtoqq)
    return z;
  
  if(FlavType()==FlavourType::gtogg)
    return x/(1.0-x)+z*(1.0-z);
  
  THROW(fatal_error, "Internal error");
}

///////////////////////////////////////////////////////////////
////////// CalcPtilde METHODS  ////////////////////////////////
///////////////////////////////////////////////////////////////

ATOOLS::Vec4D FF_Dipole::CalcPtilde() const
{
  /* \mu-\nu tensor structure in hep-ph/9605323v3 eq. (5.8), (5.9)  */
  return m_kin.m_zi*m_kin.m_pi - m_kin.m_zj*m_kin.m_pj;
}
  
ATOOLS::Vec4D FI_Dipole::CalcPtilde() const
{
  /* \mu-\nu tensor structure in hep-ph/9605323v3 eq. (5.40), (5.41)  */
  return m_kin.m_zi*m_kin.m_pi - m_kin.m_zj*m_kin.m_pj;
}
  
ATOOLS::Vec4D IF_Dipole::CalcPtilde() const
{
  /* \mu-\nu tensor structure in hep-ph/9605323v3 eq. (5.67), (5.68)  */
  return m_kin.m_pi/m_kin.m_ui - m_kin.m_pk/(1.0-m_kin.m_ui);
}

ATOOLS::Vec4D II_Dipole::CalcPtilde() const
{
  /* \mu-\nu tensor structure in hep-ph/9605323v3 eq. (5.147), (5.148)  */
  return m_kin.m_pi - (m_kin.m_pi*m_kin.m_pa)/(m_kin.m_pb*m_kin.m_pa) * m_kin.m_pb;
}

///////////////////////////////////////////////////////////////
////////// CalcB METHODS  /////////////////////////////////////
///////////////////////////////////////////////////////////////

double FF_Dipole::CalcB() const
{
  const double& zi(m_kin.m_zi);
  const double& zj(m_kin.m_zj);

  if(FlavType()==FlavourType::qtoqg)
    return 0.0;
  if(FlavType()==FlavourType::gtoqq)
    return +4.0*zi*zj;
  if(FlavType()==FlavourType::gtogg)
    return -2.0*zi*zj;

  THROW(fatal_error, "Internal error");
}
  
double FI_Dipole::CalcB() const
{
  const double& zi(m_kin.m_zi);
  const double& zj(m_kin.m_zj);

  if(FlavType()==FlavourType::qtoqg)
    return -1.0;
  if(FlavType()==FlavourType::gtoqq)
    return +4.0*zi*zj;
  if(FlavType()==FlavourType::gtogg)
    return -2.0*zi*zj;

  THROW(fatal_error, "Internal error");
}

double IF_Dipole::CalcB() const
{
  const double& x(m_kin.m_x);
    
  if(FlavType()==FlavourType::qtoqg)
    return -1.0;
  if(FlavType()==FlavourType::gtoqq)
    return -4.0*(1.0-x)/x;
  if(FlavType()==FlavourType::gtogg)
    return -2.0*(1.0-x)/x;
    
  THROW(fatal_error, "Internal error");
}

double II_Dipole::CalcB() const
{
  const double& x(m_kin.m_x);
  const double& v(m_kin.m_v);

  if(FlavType()==FlavourType::qtoqg)
    return -1.0;
      
  if(FlavType()==FlavourType::gtoqq)
    {
      switch(SubtractionType())
	{
	case 0:
	  return  -4.0*(1.0-x)/x;
	case 1:
	  return  -4.0*( (x+v)/(ATOOLS::sqr(x+v) +v*(1-x-v)) -1);
	case 2:
	  return  -4.0*(1.0/(x+v)-1.0);
	default:
	  THROW(not_implemented, "Not implemented");
	}
    }
      
  if(FlavType()==FlavourType::gtogg)
    {
      switch(SubtractionType())
	{
	case 0:
	  return -2.0*(1.0-x)/x;
	case 1:
	  return -2.0*( (x+v)/(ATOOLS::sqr(x+v) +v*(1-x-v)) -1);
	case 2:
	  return -2.0*(1.0/(x+v)-1.0);
	default:
	  THROW(not_implemented, "Not implemented");
	}
    }
      
  THROW(fatal_error, "Internal error");
}
