#include "PHASIC++/Process/CS_Dipole.H"
#include "PHASIC++/Process/Subprocess_Info.H"
#include "ATOOLS/Org/Getter_Function.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE PHASIC::CS_Dipole
#define PARAMETER_TYPE PHASIC::Dipole_Info
#include "ATOOLS/Org/Getter_Function.C"

#include "MODEL/Main/Coupling_Data.H"
#include "ATOOLS/Org/Exception.H"

#include "assert.h"

using namespace PHASIC;

std::ostream& PHASIC::operator<<(std::ostream& str,const SplittingType& st)
{
  switch(st)
    {
    case SplittingType::FF: return str<<"FF";
    case SplittingType::IF: return str<<"IF";
    case SplittingType::FI: return str<<"FI";
    case SplittingType::II: return str<<"II";
    default: THROW(fatal_error, "Internal error");
    }
}


std::ostream& PHASIC::operator<<(std::ostream& str,const FlavourType& ft)
{
  switch(ft)
    {
    case FlavourType::gtogg: return str<<"g->gg";
    case FlavourType::gtoqq: return str<<"g->qq";
    case FlavourType::qtoqg: return str<<"q->qg";
    default: THROW(fatal_error, "Internal error");
    }
}


std::ostream& PHASIC::operator<<(std::ostream& str, const Dipole_Info& di)
{
  return str << di.m_real_flavs <<
    " i=" << di.m_real_i <<
    " j=" << di.m_real_j <<
    " k=" << di.m_real_k <<
    " "   << di.m_flav_type <<
    " ["  << di.m_split_type << "]";
}


Dipole_Info::Dipole_Info(const ATOOLS::Flavour_Vector& flavs,
			 const size_t& i, const size_t& j, const size_t& k,
			 const int& subtrtype, const double& alphamin, const double& alphamax)
  : m_real_flavs(flavs), m_subtype(subtrtype), m_alphamin(alphamin), m_alphamax(alphamax)
{

  /* Position of flavours i,j,k in the real emission flavour vector */
  m_real_i = i; m_real_k = k; m_real_j = j;

  bool IS_emitter(i<2||j<2);
  bool IS_spectator(k<2);

  if(IS_emitter && IS_spectator)
    m_split_type = SplittingType::II;
  else if(IS_emitter && !IS_spectator)
    m_split_type = SplittingType::IF;
  else if(!IS_emitter && IS_spectator)
    m_split_type = SplittingType::FI;
  else if (!IS_emitter && !IS_spectator)
    m_split_type = SplittingType::FF;
  else
    THROW(fatal_error, "Internal error");

  /* g -> gg splitting */
  if(m_real_flavs[i].IsGluon() && m_real_flavs[j].IsGluon())
    m_flav_type = FlavourType::gtogg;
  /* g -> qqbar splitting */
  else if(m_real_flavs[i].IsQuark() && m_real_flavs[j].IsQuark())
    m_flav_type = FlavourType::gtoqq;
  /* q -> qg splitting */
  else 
    m_flav_type = FlavourType::qtoqg;

}


CS_Dipole::CS_Dipole(const Dipole_Info& di)
  : m_dip_info(di), p_aqcd(NULL), p_aqed(NULL)
{
  m_born_flavs = ConstructBornFlavours(I(),J(),di.m_real_flavs);
  m_id_vector  = ConstructIDVector    (I(),J(),di.m_real_flavs);
  
  PHASIC::Subprocess_Info ii;  PHASIC::Subprocess_Info fi;
  for (size_t n=0; n<2; n++)
    ii.m_ps.push_back(PHASIC::Subprocess_Info(Flavours()[n]));
  for (size_t n=2; n<Flavours().size(); n++)
    fi.m_ps.push_back(PHASIC::Subprocess_Info(Flavours()[n]));

  m_symfac = ii.ISSymmetryFactor();
  m_symfac*= fi.FSSymmetryFactor();

}

ATOOLS::Flavour_Vector CS_Dipole::ConstructBornFlavours(const size_t& i, const size_t& j,
							const ATOOLS::Flavour_Vector& flavs)
{
  /* Convention: select the smaller inded among i,j and identify it as
     'emitting' particle (important for initial state splittings).
     Construct new flavour order by replacing particle at this index
     with combined flavour (ij) and by removing the 'emitted' parton.
     Important: follow the same convention for momenta! */
  
  size_t emitter = std::min(i,j);
  size_t emitted = std::max(i,j);
  
  const ATOOLS::Flavour& fl_ij = CombinedFlavour(i,j,flavs);
  ATOOLS::Flavour_Vector ret = flavs;

  /* Now assign combined flavour to the spot of the emitter */
  ret[emitter] = fl_ij;

  /* ... and remove the emitted one */
  ret.erase(ret.begin()+emitted);
  
  return ret;
}


std::vector<size_t> CS_Dipole::ConstructIDVector(const size_t& i, const size_t& j,
						 const ATOOLS::Flavour_Vector& flavs)
{
  /* Follow the same convention here as in ConstructBornFlavours,
     since the index vector must correspond to the flavour ordering in
     the born configuration. */
  
  size_t emitter = std::min(i,j);
  size_t emitted = std::max(i,j);
  size_t combined_id = (1<<i)|(1<<j);

  std::vector<size_t> ret(flavs.size(), 0);

  /* Construct the ID vecot for the real emission config */
  for(size_t i(0); i<ret.size(); i++) ret[i] = 1<<i;

  /* Now assign combined ID to the spot of the emitter */
  ret[emitter] = combined_id;

  /* ... and remove the emitted id */
  ret.erase(ret.begin()+emitted);

  return ret;
}


ATOOLS::Flavour CS_Dipole::CombinedFlavour(const size_t& i, const size_t& j,
					   const ATOOLS::Flavour_Vector& flavs)
{
  /* Convert any incoming flavours to outgoing flavours */
  const ATOOLS::Flavour fli =  i < 2 ? flavs[i].Bar() : flavs[i];
  const ATOOLS::Flavour flj =  j < 2 ? flavs[j].Bar() : flavs[j];

  ATOOLS::Flavour flij;
  
  if (fli.IsQuark() && flj.IsQuark())
    flij = ATOOLS::Flavour(kf_gluon);
  else if (fli.IsGluon() && flj.IsGluon())
    flij = ATOOLS::Flavour(kf_gluon);
  else if (fli.IsQuark())
    flij = fli;
  else if (flj.IsQuark())
    flij = flj;
  else
    THROW(fatal_error, "Internal error");

  /* Convert outgoing combined flavour back to incoming flavour for
     initial state splittings */
  return (i<2 || j<2) ? flij.Bar() : flij;
}


void CS_Dipole::SetCouplings(MODEL::Coupling_Data* p_rqcd,
			     MODEL::Coupling_Data* p_rqed)
{
  p_aqcd=p_rqcd;
  p_aqed=p_rqed;
}


bool CS_Dipole::PassesAlphaMin() const
{
  const double& alphamin = Info().m_alphamin;
  return LastKinematics()->PassesAlphaMin(alphamin);
}


bool CS_Dipole::PassesAlphaCuts() const
{
  const double& alphamin = Info().m_alphamin;
  const double& alphamax = Info().m_alphamax;
  return LastKinematics()->PassesAlphaCuts(alphamin, alphamax);
}


///////////////////////////////////////////////////////////////
////////// FINAL FINAL ////////////////////////////////////////
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


///////////////////////////////////////////////////////////////
////////// FINAL INITIAL //////////////////////////////////////
///////////////////////////////////////////////////////////////


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


///////////////////////////////////////////////////////////////
////////// INITIAL FINAL //////////////////////////////////////
///////////////////////////////////////////////////////////////


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


///////////////////////////////////////////////////////////////
////////// INITIAL INITIAL ////////////////////////////////////
///////////////////////////////////////////////////////////////


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


////////////////////////////////////////////////////////////////////////
/////////////////////////* Getter mechanism *///////////////////////////
////////////////////////////////////////////////////////////////////////


typedef ATOOLS::Getter_Function<CS_Dipole,Dipole_Info> CS_Dipole_Getter;


CS_Dipole* CS_Dipole::Get(const Dipole_Info& di)
{
  CS_Dipole_Getter::Getter_List glist(CS_Dipole_Getter::GetGetters());
  for (CS_Dipole_Getter::Getter_List::const_iterator git(glist.begin());
       git!=glist.end();++git) {
    CS_Dipole *dip=(*git)->GetObject(di);
    if (dip) return dip;
  }
  return NULL;
}


CS_Dipole *CS_Dipole::Get(const std::string& tag,
			  const Dipole_Info& pi)
{
  CS_Dipole* dip=CS_Dipole_Getter::GetObject(tag, pi);
  if (!dip)
    THROW(fatal_error, "Did not find CS_Dipole "+tag);
  return dip;
}
