#include "EXTAMP/Dipole_Wrapper_Process.H"
#include "EXTAMP/External_ME_Interface.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Process/Spin_Color_Correlated_ME2.H"
#include "PDF/Main/NLOMC_Base.H"

#include "ATOOLS/Org/MyStrStream.H"

#include <assert.h>

using namespace EXTAMP;

Dipole_Wrapper_Process::Dipole_Wrapper_Process(const RS_Process& rsproc,
					       CS_Dipole* dipole,
					       BEAM::Beam_Spectra_Handler* beam,
					       PDF::ISR_Handler* isr) : p_dipole(dipole)

{
  /* Follow sherpa convetions here and set i<j */
  m_norm = rsproc.NormFac();



  m_born_procinfo = ConstructBornProcessInfo(rsproc.Info(),I(),J(),
					     Dipole()->FlavIJ());
  PHASIC::Process_Base::SortFlavours(m_born_procinfo);
  m_born_flavs = m_born_procinfo.ExtractFlavours();


  /* This wrapper has to call it's Init method itself, since it's not
     initialized through the framework like other procs. This
     re-orders flavours according to Sherpa conventions (reordering is
     also applied to the Process_Info instance). */
  PHASIC::Process_Base::Init(rsproc.Info(), beam, isr, 0);

  /* Not done in any of the PHASIC base classes */
  m_mincpl = rsproc.Info().m_mincpl;
  m_maxcpl = rsproc.Info().m_maxcpl;

  /* Construct an index maps such that 
     Flavours[i] == Dipole()->Flavours[m_indexmap[i]] 
     Flavours[m_inversemap] == Dipole()->Flavours[i] */
  m_indexmap   = ConstructIndexMapping(Dipole()->Flavours(),
				       BornFlavours(),NIn());
  m_inversemap = InvertIndexMapping(m_indexmap);

  /* Cross-check index mapping */
  for(size_t i(0); i<BornFlavours().size(); i++)
    {
      assert(BornFlavours()[i] == dipole->Flavours()[m_indexmap[i]]);
      assert(BornFlavours()[m_inversemap[i]] == dipole->Flavours()[i]);
    }
  
  /* Set name in accordance with Sherpa conventions so that it can be
     identified by MC@NLO process. Save name of born config as well
     for NLO_subevts. */
  m_born_name = PHASIC::Process_Base::GenerateName(m_born_procinfo.m_ii,
						   m_born_procinfo.m_fi);
  m_name = rsproc.Name()+"_RS"+
    ATOOLS::ToString(I())+"_"+
    ATOOLS::ToString(J())+"_"+
    ATOOLS::ToString(K());

  /* Set couplings: copy pointers from dipole correlated ME */
  m_cpls.insert(std::make_pair("Alpha_QCD",Dipole()->CorrelatedME()->RunningQCD()));
  m_cpls.insert(std::make_pair("Alpha_QED",Dipole()->CorrelatedME()->RunningQED()));

    /* Fill Comninable Info */
  m_cluster_flav_map
    = External_ME_Interface::ConstructCombinableMap(Flavours(),Info(),NIn());

  m_id_vector = ConstructIDVector();

  /* Cross-check index mapping in ID vector */
  for(size_t i(0); i<BornFlavours().size(); i++)
    {
      const ATOOLS::Flavour& born_flav = BornFlavours()[i];

      /* This is a vector of integers representing the particles that
	 have been clusteded into the born index i */
      std::vector<int> id = ATOOLS::ID(IDVector()[i]);

      /* Only check non-clustered flavours to avoid doing the flavour
	 combination here*/
      if (id.size()>1) continue;
      
      /* Now check that the flavour in real configuration has been
	 properly mapped */
      const ATOOLS::Flavour& real_flav = Flavours()[id.front()];
      if(real_flav != born_flav)
	THROW(fatal_error, "Inconsistent flavour mapping");
    }

  /* Cross-check index mapping in ID vector */
  for(size_t i(0); i<BornFlavours().size(); i++)
    {
      const ATOOLS::Flavour& born_flav = BornFlavours()[i];

      /* This is a vector of integers representing the particles that
	 have been clusteded into the born index i */
      std::vector<int> id = ATOOLS::ID(IDVector()[i]);

      /* Only check non-clustered flavours to avoid doing the flavour
	 combination here*/
      if (id.size()>1) continue;

      /* Now check that the flavour in real configuration has been
	 properly mapped */
      const ATOOLS::Flavour& real_flav = Flavours()[id.front()];
      if(real_flav != born_flav)
	THROW(fatal_error, "Inconsistent flavour mapping");
    }

  m_moms.resize(Dipole()->Flavours().size());
}


PHASIC::Process_Info Dipole_Wrapper_Process::ConstructBornProcessInfo
(const PHASIC::Process_Info& rsinfo,
 size_t i, size_t j, const ATOOLS::Flavour& flav_ij)
{
  PHASIC::Process_Info ret(rsinfo);
  ret.m_fi.m_nlotype&=~ATOOLS::nlo_type::real;
  ret.Combine(i,j,flav_ij);
  return ret;
}


std::vector<size_t> Dipole_Wrapper_Process::ConstructIndexMapping
(const ATOOLS::Flavour_Vector& dipole_flavs,
 const ATOOLS::Flavour_Vector& process_flavs,
 size_t nin)
{
  /* Construct Process_Info with same ordering as in dipole_flavs */
  PHASIC::Subprocess_Info in;
  PHASIC::Subprocess_Info fi;
  for(size_t i(0); i<nin; i++)
    in.m_ps.push_back(dipole_flavs[i]);
  for(size_t i(nin); i<dipole_flavs.size(); i++)
    fi.m_ps.push_back(dipole_flavs[i]);
  PHASIC::Process_Info pi(in,fi);
  assert(pi.ExtractFlavours() == dipole_flavs);
  
  /* Tag all flavours */
  for(size_t i(0); i<nin; i++)
    pi.m_ii.m_ps[i].m_tag=i;
  int cnt(nin); pi.m_fi.SetTags(cnt);
  assert(cnt==dipole_flavs.size());

  /* Apply Sherpa's re-ordering conventions */
  PHASIC::Process_Base::SortFlavours(pi);
  assert(pi.ExtractFlavours() == process_flavs);

  /* Extract map from tags and re-ordered flavours */
  std::vector<size_t> ret;
  ret.resize(dipole_flavs.size());
  for (size_t i(0);i<nin;++i)
    ret[i]=pi.m_ii.m_ps[i].m_tag;
  std::vector<int> fi_tags;
  pi.m_fi.GetTags(fi_tags);
  size_t nout = dipole_flavs.size()-nin;
  for (size_t i(0);i<nout;++i) ret[nin+i]=fi_tags[i];
  return ret;
}


std::vector<size_t> Dipole_Wrapper_Process::InvertIndexMapping
(const std::vector<size_t>& map)
{
  std::vector<size_t> ret(map.size());
  for(size_t i(0); i<map.size(); i++)
    ret[map[i]] = i;
  return ret;
}


void Dipole_Wrapper_Process::SetSubEventProperties(ATOOLS::NLO_subevt& sub)
{
  sub.p_fl     = &BornFlavours()[0];
  sub.p_mom    = &Momenta()[0];
  sub.m_n      =  BornFlavours().size();
  sub.m_i      =  I();
  sub.m_j      =  J();
  sub.m_k      =  K();
  sub.m_kt     =  BornK();
  sub.m_ijt    =  BornIJ();
  sub.p_id     = &IDVector()[0];
  sub.m_me     =  0.0;
  sub.m_result =  0.0;
  sub.m_trig   =  false;
  sub.p_proc   =  this;
  sub.m_pname  =  m_born_name;
  sub.m_pname  =  sub.m_pname.substr(0,sub.m_pname.rfind("__"));
}


void Dipole_Wrapper_Process::AssignSubEvent(ATOOLS::NLO_subevt* sub)
{
  p_subevent = sub;
}


std::vector<size_t> Dipole_Wrapper_Process::ConstructIDVector() const
{
  std::vector<size_t> tmp = Dipole()->IDVector();
  std::vector<size_t> ret = tmp;
  for(size_t i(0); i<ret.size(); i++)
    ret[i] = tmp[m_indexmap[i]];
  return ret;
}


void Dipole_Wrapper_Process::FillProcessMap(PHASIC::NLOTypeStringProcessMap_Map *apmap)
{
  /* Copied from AMEGIC::Single_DipoleTerm */
  p_apmap=apmap;
  if (p_apmap->find(ATOOLS::nlo_type::rsub)==p_apmap->end())
    (*p_apmap)[ATOOLS::nlo_type::rsub] = new PHASIC::StringProcess_Map();
  std::string fname(m_name);
  size_t len(m_pinfo.m_addname.length());
  if (len) fname=fname.erase(fname.rfind(m_pinfo.m_addname),len);
  (*(*p_apmap)[ATOOLS::nlo_type::rsub])[fname]=this;
}


bool Dipole_Wrapper_Process::Combinable(const size_t &idi,const size_t &idj)
{
  bool ret = m_cluster_flav_map.find(idi | idj)!=m_cluster_flav_map.end();
  return ret;
}

  
const ATOOLS::Flavour_Vector &Dipole_Wrapper_Process::CombinedFlavour(const size_t &idij)
{
  std::map<size_t, ATOOLS::Flavour_Vector>::const_iterator it = m_cluster_flav_map.find(idij);
  if(it==m_cluster_flav_map.end())
    THROW(fatal_error, "Internal error");
  return it->second;
}


void Dipole_Wrapper_Process::SetNLOMC(PDF::NLOMC_Base *const nlomc)
{
  PHASIC::Process_Base::SetNLOMC(nlomc);
  Dipole()->SetSubtractionType(nlomc->SubtractionType());
}


void Dipole_Wrapper_Process::CalcKinematics(const ATOOLS::Vec4D_Vector& p)
{
  if(!Dipole()) THROW(fatal_error, "Invalid dipole pointer");

  /* Dipole and dipole wrapper share the same flavour and momentum
     ordering in the real emission configuration, pass momentum on
     as-is */
  DEBUG_VAR(Dipole()->Info());
  Dipole()->CalcKinematics(p);

  /* Apply re-mapping of (Born-) momenta, invert incoming momenta for
     NLO_subevts */
  for(size_t i(0); i<NIn(); i++)
    {m_moms[i] = -(Dipole()->Momenta()[m_indexmap[i]]);
    /*DEBUG_VAR(m_moms[i]);*/}
  for(size_t i(NIn()); i<m_moms.size(); i++)
    {m_moms[i] =  (Dipole()->Momenta()[m_indexmap[i]]);
    /*DEBUG_VAR(m_moms[i]);*/}
}


double Dipole_Wrapper_Process::Calc(ATOOLS::NLO_subevt* const evt)
{
  double sign = MCModeSign(evt);
  double dxs  = (sign==0.0) ? 0.0 : sign*Dipole()->Calc()/NormFac();
  
  /* NLO_subevts of RS proc are added up by PHASIC::Single_Process,
     hence add a minus sign here */
  evt->m_me     = -dxs;
  evt->m_mewgt  = -dxs;
  evt->m_result = -dxs;
  return dxs;
}


void Dipole_Wrapper_Process::CalculateScale(const ATOOLS::Vec4D_Vector& real_p,
					    const ATOOLS::Vec4D_Vector& born_p,
					    ATOOLS::NLO_subevt* const evt)
{
  /* The scale setter knows about the RS and all other dipoles via the
     subevtlist of the RS. Scale setting might be based on the real
     configuration so we have to set the momenta of the RS here before
     calling the scale setter. */
  RS_Process* rsproc = ((RS_Process*)(evt->p_real->p_proc));
  rsproc->Integrator()->SetMomenta(real_p);

  /* Locally set flavours to born configuration in order not to
     confuse the scale setter */
  ATOOLS::Flavour_Vector tmp = m_flavs;
  m_flavs = m_born_flavs;
  ScaleSetter()->CalculateScale(born_p);
  m_flavs = tmp;
}


double Dipole_Wrapper_Process::Partonic(const ATOOLS::Vec4D_Vector &p, const int mode)
{
  CalcKinematics(p);

  /* Use the trigger method of the RS process to set the m_trig member
     of this subevent. Also triggers the real and all other dipoles
     but whatever... */
  PHASIC::Process_Base* rsproc = (PHASIC::Process_Base*)(p_subevent->p_real->p_proc);
  ATOOLS::NLO_subevtlist* subs = rsproc->GetRSSubevtList();
  rsproc->Selector()->RSTrigger(subs);

  /* This method is called in MC@NLO matching from outside an
     RS_Process. Hence have to calc scales ourselves */
  CalculateScale(p, m_moms, p_subevent);

  double dxs = Calc(p_subevent);

  return m_lastxs = m_mewgtinfo.m_B = dxs;
}


int Dipole_Wrapper_Process::MCModeSign(ATOOLS::NLO_subevt* const evt) const
{
  /* m_mcmode is nonzero if this process is used to calcualte the
     DA-DS term or the DA term as defined in arXiv:1111.1200. Since
     DA=DS apart from theta functions implementing the PS starting
     scale in DA, we just have to implement a sign \in {-1,0,1}
     here. */
  double kt2    = GetKT2ofSplitting(*Dipole()->LastKinematics());
  double kt2max = GetMaxKT2ForDA();
  evt->m_kt2    = kt2;
  
  int DS(1),DA(kt2<=kt2max*(1.0+1e-6));

  /* For the calculation of DA-DS, another minus sign is added in
     PHASIC::Single_Process. Need to compensate for that here by yet
     another minus sign. */
  if (m_mcmode==0) return + DS;
  if (m_mcmode==1) return -(DA-DS);
  if (m_mcmode==2) return + DA;
  
  THROW(fatal_error, "Unknown MC-mode "+ATOOLS::ToString(m_mcmode));
}


double Dipole_Wrapper_Process::GetKT2ofSplitting(const Dipole_Kinematics& kin) const
{
  /* Get hardness of splitting as defined by shower. Only needed for
     MC@NLO. */
  if(!p_nlomc) return 0.0;
  double x  = kin.ShowerX() ;
  double y  = kin.ShowerY() ;
  double Q2 = kin.ShowerQ2();
  return p_nlomc->KT2(*p_subevent,x,y,Q2);
}


double Dipole_Wrapper_Process::GetMaxKT2ForDA() const
{
  /* Get the maximum hardness for a splitting: PS starting scale */
  if(p_subevent->m_mu2.size()>(ATOOLS::stp::size+ATOOLS::stp::res))
    return p_subevent->m_mu2[(ATOOLS::stp::id(ATOOLS::stp::size+ATOOLS::stp::res))];
  return p_subevent->m_mu2[(ATOOLS::stp::res)];
}
