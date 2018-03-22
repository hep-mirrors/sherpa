#include "EXTAMP/CS_Dipoles.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Math/Poincare.H"

#include "assert.h"

using namespace EXTAMP;

///////////////////////////////////////////////////////////////
////////// CalcKinematics METHODS /////////////////////////////
///////////////////////////////////////////////////////////////

void FF_Dipole::CalcKinematics(const ATOOLS::Vec4D_Vector& p)
{
  assert(K()>1); assert(Emitter()>1); assert(Emitted()>1);

  switch(DipCase()){
  case CS: {
    const ATOOLS::Vec4D& pi = p[I()];
    const ATOOLS::Vec4D& pj = p[J()];
    const ATOOLS::Vec4D& pk = p[K()];

    /* Implementation of hep-ph/9605323v3 eq. (5.3) - (5.6) */
    m_kin.m_y         = pi*pj/(pi*pj+pj*pk+pk*pi);
    m_kin.m_alphamin  = m_kin.m_y; 
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
    break;
    }

  case IDa:
  case IDb:{
    const size_t emitter = I();
    const ATOOLS::Vec4D& pi = p[J()];  // emitted (J<I)
    const ATOOLS::Vec4D& pa = p[I()];  // emitter
    const ATOOLS::Vec4D& pb = p[K()];  // colour spectator
    ATOOLS::Vec4D n;
    /* e-__e+__W+__W-__j__b__bb */
    /* 0   1   2   3   4  5  6  */
    switch(DipCase()){
      case IDa:
        if      (emitter==5) n = p[0]+p[1]-pa-p[3]-p[6];   // b-quark (5) emitts
        else if (emitter==6) n = p[0]+p[1]-pa-p[2]-p[5];   // bbar-quark (6) emitts
        else    THROW(fatal_error, "Invalid emitter.");
        break;
      case IDb:
        if      (emitter==5) n = p[0]+p[1]-pa-p[2];
        else if (emitter==6) n = p[0]+p[1]-pa-p[3];
        else    THROW(fatal_error, "Invalid emitter.");
    }

    m_kin.m_vitilde   = pa*pi/(pa*n);
    m_kin.m_alphamin  = m_kin.m_vitilde; 
    m_kin.m_viab      = pa*pb/((pa+pb)*pi);
    m_kin.m_zain      = pa*n/((pa+pi)*n);
    m_kin.m_xain      = (pa-pi)*n/(pa*n);
    m_kin.m_pij_tilde = pa/m_kin.m_zain;
    m_kin.m_pi        = pi;
    m_kin.m_pa        = pa;

    /* setting momenta of underlying Born config. */
    m_kin.m_born_mom = p;         
    /* Apply transformation of spectator(s) */
    ATOOLS::Vec4D Ka      = n-pi;
    ATOOLS::Vec4D Katilde = n-(1-m_kin.m_xain)*pa;
    std::vector<size_t> kj;
    switch(DipCase()){
      case IDa:
        if      (emitter==5) kj = {2}; //  b-quark(5) emitts,    W+(2) takes recoil
        else if (emitter==6) kj = {3}; //  bbar-quark(6) emitts, W-(3) takes recoil
        break;
      case IDb:
        if      (emitter==5) kj = {3,6};
        else if (emitter==6) kj = {2,5};
        break;
    }
    for(auto j: kj){
      m_kin.m_born_mom[j] = p[j]-2.0*p[j]*(Ka+Katilde)/(Ka+Katilde).Abs2()*(Ka+Katilde)
                            +2.0*(p[j]*Ka)/Ka.Abs2()*Katilde;
    }      
    if(DipCase()==IDb){
      if      (emitter==5) m_kin.m_viab = pa*m_kin.m_born_mom[6]/((pa+m_kin.m_born_mom[6])*pi);
      else if (emitter==6) m_kin.m_viab = pa*m_kin.m_born_mom[5]/((pa+m_kin.m_born_mom[5])*pi);
    }

      /* Replace emitter momentum with combined momentum of (ij) and
         remove emitted. */
      m_kin.m_born_mom[Emitter()] = m_kin.m_pij_tilde;
      m_kin.m_born_mom.erase(m_kin.m_born_mom.begin()+Emitted());
    break;
    }

  case IDin:{ // (J<I)
    /* p __p __W+__W-__j__j__j */
    /* 0   1   2   3   4  5  6  */
    ATOOLS::Vec4D pi;  // emitted
    ATOOLS::Vec4D pa;  // emitter
    ATOOLS::Vec4D pb;  // colour spectator
    ATOOLS::Vec4D n;   // aux vector: p_in - sum p_id
    std::vector<size_t> kj;  // kinematic spectator

    if( FlavI().IsbQuark() && FlavJ().IsGluon() ){ /* b emitts g */
    pa = p[I()];
    pi = p[J()];
    pb = p[K()];
    n  = p[0]+p[1]-pa-p[3]-pb;
    kj = {2}; // W+(2) takes recoil
    }
    else if( FlavJ().IsbQuark() && FlavI().IsGluon() ){ /* b emitts g */
    pa = p[J()];
    pi = p[I()];
    pb = p[K()];
    n  = p[0]+p[1]-pa-p[3]-pb;
    kj = {2}; // W+(2) takes recoil
    }
    else if( FlavI().IsbbarQuark() && FlavJ().IsGluon() ){ /* bbar emitts g */
    pa = p[I()];
    pi = p[J()];
    pb = p[K()];
    n  = p[0]+p[1]-pa-p[2]-pb;
    kj = {3}; // W-(3) takes recoil
    }
    else if( FlavJ().IsbbarQuark() && FlavI().IsGluon() ){ /* bbar emitts g */
    pa = p[J()];
    pi = p[I()];
    pb = p[K()];
    n  = p[0]+p[1]-pa-p[2]-pb;
    kj = {3}; // W-(3) takes recoil
    }
//    else if( FlavI().IsbQuark() && FlavJ().IsbbarQuark() ){ /* g emitts b & bbar */
//    pa = p[I()];
//    pi = p[J()];
//    pb = p[K()];
//    // TODO: better assignment?
//    n  = p[0]+p[1]-pa-pb;
//    kj = {2,3}; // W+(2) & W-(3) take recoil
//    }
//    else if( FlavJ().IsbQuark() && FlavI().IsbbarQuark() ){ /* g emitts b & bbar */
//    pa = p[J()];
//    pi = p[I()];
//    pb = p[K()];
//    // TODO: better assignment?
//    n  = p[0]+p[1]-pa-pb;
//    kj = {2,3}; // W+(2) & W-(3) take recoil
//    }
    else {
      DEBUG_VAR(RealFlavours()); 
      DEBUG_VAR(I());
      DEBUG_VAR(J());
      DEBUG_VAR(K());
      THROW(fatal_error, "Invalid assignment.");
    }

    /* set kinematic variables */
    m_kin.m_vitilde   = pa*pi/(pa*n);
    m_kin.m_alphamin  = m_kin.m_vitilde; 
    m_kin.m_viab      = pa*pb/((pa+pb)*pi);
    m_kin.m_zain      = pa*n/((pa+pi)*n);
    m_kin.m_xain      = (pa-pi)*n/(pa*n);
    m_kin.m_pij_tilde = pa/m_kin.m_zain;
    m_kin.m_pi        = pi;
    m_kin.m_pa        = pa;
    m_kin.m_n         = n;

    /* setting momenta of underlying Born config. */
    m_kin.m_born_mom = p;         
    /* Apply transformation of kinematic spectator(s) */
    ATOOLS::Vec4D Ka      = n-pi;
    ATOOLS::Vec4D Katilde = n-(1-m_kin.m_xain)*pa;
    for(auto j: kj){
      m_kin.m_born_mom[j] = p[j]-2.0*p[j]*(Ka+Katilde)/(Ka+Katilde).Abs2()*(Ka+Katilde)
                            +2.0*(p[j]*Ka)/Ka.Abs2()*Katilde;
    }      
    /* Replace emitter momentum with combined momentum of (ij) and
       remove emitted. */
    m_kin.m_born_mom[Emitter()] = m_kin.m_pij_tilde;
    m_kin.m_born_mom.erase(m_kin.m_born_mom.begin()+Emitted());
    break;
    }
  }
}

void FI_Dipole::CalcKinematics(const ATOOLS::Vec4D_Vector& p)
{
  assert(K()<2); assert(I()>1); assert(J()>1);

  switch(DipCase()){
  case IDin:
  case CS: {
    const ATOOLS::Vec4D& pi = p[I()];
    const ATOOLS::Vec4D& pj = p[J()];
    const ATOOLS::Vec4D& pa = p[K()];
    
    /* Implementation of hep-ph/9605323v3 (5.37) - (5.42) with a=k */
    m_kin.m_x         = (pi*pa + pj*pa - pi*pj)/((pi+pj)*pa);
    m_kin.m_alphamin  = 1 - m_kin.m_x;
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
    break;
    }
//  case IDin:{ // (J<I)
//    /* p __p __W+__W-__j__j__j */
//    /* 0   1   2   3   4  5  6  */
//    ATOOLS::Vec4D pi;  // emitted
//    ATOOLS::Vec4D pa;  // emitter
//    ATOOLS::Vec4D pb;  // colour spectator
//    ATOOLS::Vec4D n;   // aux vector: p_in - sum p_id
//    std::vector<size_t> kj;  // kinematic spectator
//
//    if( FlavI().IsbQuark() && FlavJ().IsGluon() ){ /* b emitts g */
//    pa = p[I()];
//    pi = p[J()];
//    pb = p[K()];
//    n  = p[2]+p[J()];
//    kj = {2}; // W+(2) takes recoil
//    }
//    else if( FlavJ().IsbQuark() && FlavI().IsGluon() ){ /* b emitts g */
//    pa = p[J()];
//    pi = p[I()];
//    pb = p[K()];
//    n  = p[2]+p[I()];
//    kj = {2}; // W+(2) takes recoil
//    }
//    else if( FlavI().IsbbarQuark() && FlavJ().IsGluon() ){ /* bbar emitts g */
//    pa = p[I()];
//    pi = p[J()];
//    pb = p[K()];
//    n  = p[3]+p[J()];
//    kj = {3}; // W-(3) takes recoil
//    }
//    else if( FlavJ().IsbbarQuark() && FlavI().IsGluon() ){ /* bbar emitts g */
//    pa = p[J()];
//    pi = p[I()];
//    pb = p[K()];
//    n  = p[3]+p[I()];
//    kj = {3}; // W-(3) takes recoil
//    }
//    else if( FlavI().IsbQuark() && FlavJ().IsbbarQuark() ){ /* g emitts b & bbar */
//    pa = p[I()];
//    pi = p[J()];
//    pb = p[K()];
//    // TODO: better assignment?
//    n  = p[2]+p[3]+pi;
//    kj = {2,3}; // W+(2) & W-(3) take recoil
//    }
//    else if( FlavJ().IsbQuark() && FlavI().IsbbarQuark() ){ /* g emitts b & bbar */
//    pa = p[J()];
//    pi = p[I()];
//    pb = p[K()];
//    // TODO: better assignment?
//    n  = p[2]+p[3]+pi;
//    kj = {2,3}; // W+(2) & W-(3) take recoil
//    }
//    else THROW(fatal_error, "Invalid assignment.");
//
//    /* set kinematic variables */
//    m_kin.m_vitilde   = pa*pi/(pa*n);
//    m_kin.m_alphamin  = m_kin.m_vitilde; 
//    m_kin.m_viab      = pa*pb/((pa+pb)*pi);
//    m_kin.m_zain      = pa*n/((pa+pi)*n);           // zain aka zjin
//    m_kin.m_xain      = (pa-pi)*n/(pa*n);
//    m_kin.m_pij_tilde = pa/m_kin.m_zain;
//    m_kin.m_pi        = pi;                         // needed in CalcB
//    m_kin.m_pa        = pa;                         // needed in CalcB
//    m_kin.m_n         = n;                          // needed in CalcB
//
//    /* setting momenta of underlying Born config. */
//    m_kin.m_born_mom = p;         
//    /* Apply transformation of kinematic spectator(s) */
//    ATOOLS::Vec4D Ka      = n-pi;
//    ATOOLS::Vec4D Katilde = n-(1-m_kin.m_xain)*pa;
//    for(auto j: kj){
//      m_kin.m_born_mom[j] = p[j]-2.0*p[j]*(Ka+Katilde)/(Ka+Katilde).Abs2()*(Ka+Katilde)
//                            +2.0*(p[j]*Ka)/Ka.Abs2()*Katilde;
//    }      
//    /* Replace emitter momentum with combined momentum of (ij) and
//       remove emitted. */
//    m_kin.m_born_mom[Emitter()] = m_kin.m_pij_tilde;
//    m_kin.m_born_mom.erase(m_kin.m_born_mom.begin()+Emitted());
//    break;
//    }
  default:
     THROW(fatal_error, "Invalid dipole case!");
     break;
  }
}

void IF_Dipole::CalcKinematics(const ATOOLS::Vec4D_Vector& p)
{
  /* Implementation of hep-ph/9605323v3 (5.62) - (5.64) */

  assert(K()>1); assert(Emitter()<2); assert(Emitted()>1);

  switch(DipCase()){
  case CS:
  case IDin: {
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
    break;
    }
  }
}

void II_Dipole::CalcKinematics(const ATOOLS::Vec4D_Vector& p)
{
  /* Implementation of hep-ph/9605323v3 (5.137) - (5.140) */

  switch(DipCase()){
  case CS:
  case IDin: {
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
    break;
    }
  }
}

///////////////////////////////////////////////////////////////
////////// CalcKinDependentPrefac METHODS /////////////////////
///////////////////////////////////////////////////////////////

double FF_Dipole::CalcKinDependentPrefac() const
{
  if(SubtractionType()!=0) THROW(not_implemented, "Not implemented");

  switch(DipCase()){
  case CS: {
    const ATOOLS::Vec4D& pi = m_kin.m_pi;
    const ATOOLS::Vec4D& pj = m_kin.m_pj;
  
    /* hep-ph/9605323v3 eq. (5.2) */
    return -1.0/(2.0*pi*pj);
    }
  case IDa:
  case IDb:
  case IDin: {
    const ATOOLS::Vec4D& pi = m_kin.m_pi;
    const ATOOLS::Vec4D& pa = m_kin.m_pa;
  
    return -1.0/(2.0*pi*pa);
    }
  }
}

double FI_Dipole::CalcKinDependentPrefac() const
{
  if(SubtractionType()!=0) THROW(not_implemented, "Not implemented");

  switch(DipCase()){
  case CS: {
    const ATOOLS::Vec4D& pi = m_kin.m_pi;
    const ATOOLS::Vec4D& pj = m_kin.m_pj;
    const double& x = m_kin.m_x;
  
    /* hep-ph/9605323v3 eq. (5.36) */
    return -1.0/(2.0*pi*pj*x);
    }
  case IDin: {
    const ATOOLS::Vec4D& pi = m_kin.m_pi;
    const ATOOLS::Vec4D& pa = m_kin.m_pa;
  
    return -1.0/(2.0*pi*pa);
    }
  }
}

double IF_Dipole::CalcKinDependentPrefac() const
{
  if(SubtractionType()!=0) THROW(not_implemented, "Not implemented");

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
  switch(DipCase()){
  case CS:{ 
    double zi = m_kin.m_zi;
    double zj = m_kin.m_zj;
    const double& y(m_kin.m_y);
    
    /* q->qg expression depends on the flavour assignment being
       i=quark, j=gluon. Need to respect that here by swapping z_i,z_j
       if neccessary. Does not affect other splittings, so can be done
       for all cases. */
    /* pseudo-dipoles: In my test-case, i=gluon and j=quark. Kinematic variables are
               defined appropriatly. */
    if(FlavI().IsGluon()) std::swap(zi,zj);
    
    /* Coefficients of \delta_{ss^\prime} or -g^{\mu\nu} in eq. (5.7)
       - (5.9). */
    switch(FlavType()){
    case FlavourType::qtoqg:
      return 2.0/(1.0-zi*(1.0-y)) - (1.+ zi);
    case FlavourType::gtoqq:
      return 1.0;
    case FlavourType::gtogg:
      return 1.0/(1.0-zi*(1.0-y)) + 1.0/(1-zj*(1.0-y)) - 2.0;
    default:
      THROW(fatal_error, "Internal error");
    }
  }

  case IDa: 
  case IDb: 
  case IDin: 
    switch(FlavType()){
    case FlavourType::qtoqg:
      return 2.*m_kin.m_viab/m_kin.m_zain - (1+m_kin.m_zain);
//    case FlavourType::gtoqq:
//      return 1.;
//    case FlavourType::gtogg:
//      THROW(fatal_error, "g>gg in final-state not possible with psuedo-dipoles");
    }
  }
  
  THROW(fatal_error, "Internal error");
}

double FI_Dipole::CalcA() const
{
  switch(DipCase()){
  case IDin: 
  case CS:{ 
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
    switch(FlavType()){
    case FlavourType::qtoqg:
      return 2.0/(1.0-zi+(1.0-x)) - (1.+ zi);
    case FlavourType::gtoqq:
      return 1.0;
    case FlavourType::gtogg:
      return 1.0/(1.0-zi+(1.0-x)) + 1.0/(1-zj+(1.0-x)) - 2.0;
    }
  }
//  case IDin: 
//    switch(FlavType()){
//    case FlavourType::qtoqg:
//      return 2.*m_kin.m_viab/m_kin.m_zain - (1+m_kin.m_zain);
//    case FlavourType::gtoqq:
//      return 1.;
//    case FlavourType::gtogg:
//      THROW(fatal_error, "g>gg in final-state not possible with psuedo-dipoles");
//      return 0.;
//    }
  }
  
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
  switch(DipCase()){
  case IDa:
  case IDb:
  case IDin:
  case CS: 
    /* \mu-\nu tensor structure in hep-ph/9605323v3 eq. (5.8), (5.9)  */
    return m_kin.m_zi*m_kin.m_pi - m_kin.m_zj*m_kin.m_pj;
//  case IDin:
//    return m_kin.m_n*m_kin.m_pa/(m_kin.m_pi*m_kin.m_pa)*m_kin.m_pi - m_kin.m_n;
  }
}
  
ATOOLS::Vec4D FI_Dipole::CalcPtilde() const
{
  switch(DipCase()){
  case IDin:
  case CS: 
    /* \mu-\nu tensor structure in hep-ph/9605323v3 eq. (5.40), (5.41)  */
    return m_kin.m_zi*m_kin.m_pi - m_kin.m_zj*m_kin.m_pj;
//  case IDin:
//    return m_kin.m_n*m_kin.m_pa/(m_kin.m_pi*m_kin.m_pa)*m_kin.m_pi - m_kin.m_n;
  }
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
  switch(DipCase()){
  case IDa:
  case IDb:
  case IDin: // ID applies only for q>qg case: no diffenrence to CS
  case CS:{ 
    const double& zi(m_kin.m_zi);
    const double& zj(m_kin.m_zj);

    switch(FlavType()){
    case FlavourType::qtoqg:
      return 0.0;
    case FlavourType::gtoqq:
      return +4.0*zi*zj;
    case FlavourType::gtogg:
      return -2.0*zi*zj;
    }
    }
//  case IDin:{
//    const double& zain(m_kin.m_zain);
//    const double  zajn(1.-m_kin.m_zain);
//
//    switch(FlavType()){
//    case FlavourType::qtoqg:
//      return 0.0;
//    case FlavourType::gtoqq:
//      return +4.0*zain*zajn;
//    }
//    }
  }
  THROW(fatal_error, "Internal error");
}
  
double FI_Dipole::CalcB() const
{
  switch(DipCase()){
  case IDin:
  case CS:{ 
    const double& zi(m_kin.m_zi);
    const double& zj(m_kin.m_zj);
  
    switch(FlavType()){
    case FlavourType::qtoqg:
      return -1.0;
    case FlavourType::gtoqq:
      return +4.0*zi*zj;
    case FlavourType::gtogg:
      return -2.0*zi*zj;
    }
  }
//  case IDin:{
//    const double& zain(m_kin.m_zain);
//    const double  zajn(1.-m_kin.m_zain);
//
//    switch(FlavType()){
//    case FlavourType::qtoqg:
//      return 0.0;
//    case FlavourType::gtoqq:
//      return +4.0*zain*zajn;
//    }
//    }
  }

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
