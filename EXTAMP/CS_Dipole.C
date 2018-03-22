#include "EXTAMP/CS_Dipole.H"
#include "PHASIC++/Process/Spin_Color_Correlated_ME2.H"

#include "MODEL/Main/Coupling_Data.H"
#include "ATOOLS/Org/Exception.H"

using namespace EXTAMP;

std::ostream& EXTAMP::operator<<(std::ostream& str,const SplittingType& st)
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


std::ostream& EXTAMP::operator<<(std::ostream& str,const FlavourType& ft)
{
  switch(ft)
    {
    case FlavourType::gtogg: return str<<"g->gg";
    case FlavourType::gtoqq: return str<<"g->qq";
    case FlavourType::qtoqg: return str<<"q->qg";
    default: THROW(fatal_error, "Internal error");
    }
}


std::ostream& EXTAMP::operator<<(std::ostream& str,const DipoleCase& dc)
{
  switch(dc)
    {
    case DipoleCase::CS:   return str<<"CS";
    case DipoleCase::IDa:  return str<<"IDa";
    case DipoleCase::IDb:  return str<<"IDb";
    case DipoleCase::IDin: return str<<"IDin";
    default: THROW(fatal_error, "Internal error");
    }
}


std::ostream& EXTAMP::operator<<(std::ostream& str, const Dipole_Info& di)
{
  return str << di.m_real_flavs <<
    " i=" << di.m_real_i <<
    " j=" << di.m_real_j <<
    " k=" << di.m_real_k <<
    " "   << di.m_flav_type <<
    " ["  << di.m_split_type << "]" <<
    " ("  << di.m_case << ")";
}


Dipole_Info::Dipole_Info(const ATOOLS::Flavour_Vector& flavs,
			 const size_t& i, const size_t& j, const size_t& k,
             const int& subtrtype, const double& alphamin, const double& alphamax, DipoleCase dpc)
  : m_real_flavs(flavs), m_subtype(subtrtype), m_alphamin(alphamin), m_alphamax(alphamax)
{

  /* Position of flavours i,j,k in the real emission flavour vector */
  m_real_i = i; m_real_k = k; m_real_j = j;
  m_case = dpc;

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
  : m_dip_info(di)
{
  m_born_flavs = ConstructBornFlavours(I(),J(),di.m_real_flavs);
  m_id_vector  = ConstructIDVector    (I(),J(),di.m_real_flavs);

  switch(FlavType())
    {
    case EXTAMP::FlavourType::qtoqg:
      m_const_prefac = 8.0*M_PI;
      break;
      
    case EXTAMP::FlavourType::gtoqq:
      m_const_prefac = 8.0*M_PI*m_TR/m_CA;
      break;
      
    case EXTAMP::FlavourType::gtogg:
      m_const_prefac = 16.0*M_PI;
      break;
      
    default:
      THROW(fatal_error, "Internal error");
    }
  

  switch (SubtractionType())
    {
    case 0: break;
    case 1: break;
    case 2: break;
    default: THROW(not_implemented, "Subtraction type "+
		   ATOOLS::ToString(SubtractionType())+
		   " not implemented");
    }

  /* TODO: pass orders correctly!! */
  PHASIC::External_ME_Args args(ATOOLS::Flavour_Vector(Flavours().begin(),Flavours().begin()+2),
				ATOOLS::Flavour_Vector(Flavours().begin()+2, Flavours().end()),
				std::vector<double>());
  p_corr_me = PHASIC::Spin_Color_Correlated_ME2::GetME2(args);
  if(!p_corr_me) THROW(fatal_error, "Could not find correlated ME for this process.");

}


double CS_Dipole::Calc() const
{
  double alphas = p_corr_me->AlphaQCD();

  return alphas * m_const_prefac * CalcKinDependentPrefac() * CalcCorrelator();
}


double CS_Dipole::CalcCorrelator() const
{
  /* <1,...,m;a,b| T_ij T_k |b,a;m,...m1> * T_ij^{-2} */
  double TijTk = p_corr_me->CalcColorCorrelator(Momenta(), BornIJ(), BornK());
  
  /* <1,...,m;a,b| ptilde^\mu T_ij T_k ptilde^\nu |b,a;m,...m1> * T_ij^{-2} * ptilde^{-2} */
  double SC = (FlavType()==FlavourType::qtoqg) ? 0.0 : p_corr_me->CalcSpinCorrelator(Momenta(), CalcPtilde(), BornIJ(), BornK());
  
  return ( CalcA()*TijTk + CalcB()*SC );
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


bool CS_Dipole::PassesAlphaMin(const ATOOLS::Vec4D_Vector& p) const
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

