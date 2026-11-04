#include "PHASIC++/Process/MCatNLO_Process.H"

#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "PHASIC++/Selectors/Jet_Finder.H"
#include "PHASIC++/Main/Event_Reader.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Process/ME_Generators.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "PHASIC++/Process/Single_Process.H"
#include "PHASIC++/Channels/BBar_Multi_Channel.H"
#include "PDF/Main/NLOMC_Base.H"
#include "PDF/Main/Shower_Base.H"
#include "PDF/Main/Jet_Criterion.H"
#include "PDF/Main/Cluster_Definitions_Base.H"
#include "PDF/Main/ISR_Handler.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/Weight_Info.H"

#include <cassert>

using namespace ATOOLS;
using namespace PHASIC;
using namespace PDF;

MCatNLO_Process::MCatNLO_Process
(ME_Generators& gens,NLOTypeStringProcessMap_Map *pmap):
  m_gens(gens), p_bviproc(NULL),
  p_viprocqcd(NULL), p_rsprocqcd(NULL), p_ddprocqcd(NULL),
  p_viprocew(NULL), p_rsprocqed(NULL), p_ddprocqed(NULL),
  p_bproc(NULL), p_rprocqcd(NULL), p_rprocqed(NULL),
  m_qcd(false), m_ew(false),
  p_ampl(NULL)
{
  static bool ref(false);
  p_apmap=pmap;
  if (!ref) {
    ref=true;
    rpa->gen.AddCitation
      (1,"The automation of MCatNLO is published in \\cite{Hoeche:2011fd}.");
  }
  m_wassevent=false;
}

MCatNLO_Process::~MCatNLO_Process()
{
  if (p_ampl) p_ampl->Delete();
  if (p_rprocqcd) delete p_rprocqcd;
  if (p_rprocqed) delete p_rprocqed;
  if (p_bproc) delete p_bproc;
  if (p_ddprocqcd) delete p_ddprocqcd;
  if (p_ddprocqed) delete p_ddprocqed;
  if (p_rsprocqcd) delete p_rsprocqcd;
  if (p_rsprocqed) delete p_rsprocqed;
  if (p_viprocqcd) delete p_viprocqcd;
  if (p_viprocew) delete p_viprocew;
  if (p_bviproc) delete p_bviproc;
}

void MCatNLO_Process::Init(const Process_Info &pi,
			  BEAM::Beam_Spectra_Handler *const beam,
			   PDF::ISR_Handler *const isr,
			   YFS::YFS_Handler *const yfs,
			   const int mode)
{
  RegisterDefaults();
  Scoped_Settings s{ Settings::GetMainSettings()["MC@NLO"] };
  nlo_type::code bvicode = (nlo_type::code) s["PP_BVI_MODE"].Get<int>();
  if (!bvicode&nlo_type::born) THROW(not_implemented, "B required for MC@NLO");

  // build process infos for every process needed, then initialise this process
  // we need: B, VI, DD, and RS for S, DA-DS, and H events, respectively
  // in addition we need B and R for later interface to shower,
  // possible spin correlations etc.

  // 1) build a process info for this MCatNLO_Process,
  //    add the extra parton according to NLO corrections type,
  //    and then init ourselves
  //    -> this will become m_pinfo
  Process_Info cpi(pi);
  // book-keep user settings,
  // if not set, need Amegic for BVI and DADS, Comix for RS, B, R
  // (for historic reasons, still default to Amegic for B)
  std::string usermegen(pi.m_megenerator), userrsmegen(pi.m_rsmegenerator);
  std::string defbvimegen("Amegic"), defddmegen("Amegic");
  std::string defrsmegen("Comix"), defbmegen("Amegic"), defrmegen("Comix");
  cpi.m_fi.SetNLOType(nlo_type::born|nlo_type::loop|
		      nlo_type::vsub|nlo_type::real|nlo_type::rsub);

  if (pi.m_fi.m_nlocpl.size()>2) THROW(fatal_error, "MC@NLO only valid for 1 or 2 sets of NLO orders");
  // identify whether we want QCD, EW or both
  for (size_t i(0);i<pi.m_fi.m_nlocpl.size();++i) {
    if      (pi.m_fi.m_nlocpl[i][0]==1. && pi.m_fi.m_nlocpl[i][1]==0.) m_qcd=true;
    else if (pi.m_fi.m_nlocpl[i][0]==0. && pi.m_fi.m_nlocpl[i][1]==1.) m_ew=true;
    else msg_Error() << "Ignoring invalid NLO order " << m_pinfo.m_fi.m_nlocpl[i] << '\n';
  }
  if (!m_qcd && !m_ew) THROW(fatal_error, "No valid NLO orders found");
  if (p_read && m_ew) THROW(not_implemented, "HDF5 readin is only available for MC@NLO QCD");

  if (m_ew) // MC@NLO (QCD+)EW
    cpi.m_fi.m_ps.push_back(Subprocess_Info(kf_ewjet,"",""));
  else // MC@NLO QCD only
    cpi.m_fi.m_ps.push_back(Subprocess_Info(kf_jet,"",""));

  Process_Base::Init(cpi,beam,isr,yfs);
  m_pinfo.m_fi.SetNLOType(pi.m_fi.NLOType());

  // 2) build a process info for generating the name of this MCatNLO_Process,
  //    remove the extra parton again to avoid confusion in the output
  Process_Info npi(m_pinfo);
  npi.m_fi.m_ps.pop_back();
  m_name=GenerateName(m_pinfo.m_ii,npi.m_fi);

  // sanity check modes are set correctly
  if (pi.Has(nlo_type::real)!=pi.Has(nlo_type::rsub))
    THROW(fatal_error, "R/S can't be initialised separately.");

  // init B proc
  Process_Info spi(pi);
  ++spi.m_fi.m_nmax;
  spi.m_fi.SetNLOType(cpi.m_fi.NLOType());
  Process_Info bpi(spi);
  if (bpi.m_megenerator=="") bpi.m_megenerator=defbmegen;
  p_bproc=InitProcess(bpi,nlo_type::lo,false);
  if (p_bproc==NULL) return;

  // build process infos
  Process_Info bvipi(spi), sqcdpi(spi), sewpi(spi);
  bvipi.m_fi.m_nlocpl={{0.,0.}};
  sqcdpi.m_fi.m_nlocpl={{1.,0.}};
  sqcdpi.m_maxcpl[0]+=1.;
  sqcdpi.m_mincpl[0]+=1.;
  sewpi.m_fi.m_nlocpl={{0.,1.}};
  sewpi.m_maxcpl[1]+=1.;
  sewpi.m_mincpl[1]+=1.;

  // QCD BVI can't be separated for HDF5 so we don't use p_viprocqcd
  if (bvipi.m_megenerator=="") bvipi.m_megenerator=defbvimegen;
  if (p_read) {
    bvipi.m_fi.m_nlocpl={{1.,0.}};
    p_bviproc=InitProcess(bvipi,bvicode,false);
  }
  else p_bviproc=InitProcess(bvipi,nlo_type::lo,false);

  if (m_qcd) {
    // init R proc
    Process_Info rqcdpi(sqcdpi);
    rqcdpi.m_megenerator=rqcdpi.m_rsmegenerator;
    if (rqcdpi.m_megenerator=="") rqcdpi.m_megenerator=defrmegen;
    p_rprocqcd=InitProcess(rqcdpi,nlo_type::lo,true);

    if (!p_read) {
      // init VI proc
      Process_Info viqcdpi(sqcdpi);
      if (viqcdpi.m_megenerator=="") viqcdpi.m_megenerator=defbvimegen;
      p_viprocqcd=InitProcess(viqcdpi, bvicode&(nlo_type::loop|nlo_type::vsub), false);
      p_viprocqcd->FillProcessMap(p_apmap);
    }

    // init DADS proc
    Process_Info ddqcdpi(sqcdpi);
    if (ddqcdpi.m_megenerator=="") ddqcdpi.m_megenerator=defddmegen;
    p_ddprocqcd=InitProcess(ddqcdpi,nlo_type::rsub,1);

    // init RS proc
    Process_Info rsqcdpi(sqcdpi);
    rsqcdpi.m_integrator=rsqcdpi.m_rsintegrator;
    rsqcdpi.m_itmin=rsqcdpi.m_rsitmin;
    rsqcdpi.m_itmax=rsqcdpi.m_rsitmax;
    rsqcdpi.m_megenerator=rsqcdpi.m_rsmegenerator;
    if (rsqcdpi.m_megenerator=="") rsqcdpi.m_megenerator=defrsmegen;
    p_rsprocqcd=InitProcess(rsqcdpi,nlo_type::real|nlo_type::rsub,1|2);

    p_rsprocqcd->FillProcessMap(p_apmap);
    p_ddprocqcd->FillProcessMap(p_apmap);
    p_rprocqcd->FillProcessMap(p_apmap);

    // set flags to book keeping processes
    p_rprocqcd->SetLookUp(false);
    p_rprocqcd->SetParent(this);
  }

  if (m_ew) {
    // init R proc
    Process_Info rqedpi(sewpi);
    rqedpi.m_megenerator=rqedpi.m_rsmegenerator;
    if (rqedpi.m_megenerator=="") rqedpi.m_megenerator=defrmegen;
    p_rprocqed=InitProcess(rqedpi,nlo_type::lo|nlo_type::qed,true);

    // init VI proc
    Process_Info viewpi(sewpi);
    if (viewpi.m_megenerator=="") viewpi.m_megenerator=defbvimegen;
    p_viprocew=InitProcess(viewpi, bvicode&(nlo_type::loop|nlo_type::vsub)|nlo_type::qed, false);

    // init DADS proc
    Process_Info ddqedpi(sewpi);
    if (ddqedpi.m_megenerator=="") ddqedpi.m_megenerator=defddmegen;
    p_ddprocqed=InitProcess(ddqedpi,nlo_type::rsub|nlo_type::qed,1);

    // init RS proc
    Process_Info rsqedpi(sewpi);
    rsqedpi.m_integrator=rsqedpi.m_rsintegrator;
    rsqedpi.m_itmin=rsqedpi.m_rsitmin;
    rsqedpi.m_itmax=rsqedpi.m_rsitmax;
    rsqedpi.m_megenerator=rsqedpi.m_rsmegenerator;
    if (rsqedpi.m_megenerator=="") rsqedpi.m_megenerator=defrsmegen;
    p_rsprocqed=InitProcess(rsqedpi,nlo_type::real|nlo_type::rsub|nlo_type::qed,1|2);

    p_rsprocqed->FillProcessMap(p_apmap);
    p_viprocew->FillProcessMap(p_apmap);
    p_ddprocqed->FillProcessMap(p_apmap);
    p_rprocqed->FillProcessMap(p_apmap);

    // set flags to book keeping processes
    p_rprocqed->SetLookUp(false);
    p_rprocqed->SetParent(this);
  }

  p_bviproc->FillProcessMap(p_apmap);
  p_bproc->FillProcessMap(p_apmap);

  // set flags to book keeping processes
  p_bproc->SetLookUp(false);
  p_bproc->SetParent(this);

  m_psmode   = s["PSMODE"].Get<int>();
  m_hpsmode  = s["HPSMODE"].Get<int>();
  if (m_hpsmode == -1) {
    if (m_pinfo.m_ckkw&1) m_hpsmode=0;
    else m_hpsmode=4;
  }
  m_kfacmode = s["KFACTOR_MODE"].Get<int>();
  m_fomode   = s["FOMODE"].Get<int>();
  m_rsscale  = s["RS_SCALE"].Get<std::string>();
  if (m_fomode!=2) {
    p_bviproc->SetMCMode(1);
    if (m_qcd) {
      if (p_read){
        p_bviproc->SetSProc(p_ddprocqcd);
      }
      else {
        p_viprocqcd->SetSProc(p_ddprocqcd);
        p_viprocqcd->SetMCMode(1);
      }
      p_ddprocqcd->SetMCMode(2);
      p_rsprocqcd->SetMCMode(2);
    }
    if (m_ew) {
      p_viprocew->SetSProc(p_ddprocqed);
      p_viprocew->SetMCMode(1);
      p_ddprocqed->SetMCMode(2);
      p_rsprocqed->SetMCMode(2);
    }
  }

  // sanity checks:
  // whether we have a corresponding B/R process for every BVI/RS process
  if (p_bproc->Size()!=p_bviproc->Size()) {
    for (size_t i(0);i<p_bproc->Size();++i)
      msg_Debugging()<<"["<<i<<"]B:   "<<(*p_bproc)[i]->Name()<<std::endl;
    for (size_t i(0);i<p_bviproc->Size();++i)
      msg_Debugging()<<"["<<i<<"]BVI: "<<(*p_bviproc)[i]->Name()<<std::endl;
    THROW(fatal_error,"B and BVI have different size");
  }
  if (m_qcd) {
    if (p_rsprocqcd->Size()!=p_rprocqcd->Size()) {
      for (size_t i(0);i<p_rprocqcd->Size();++i)
        msg_Debugging()<<"["<<i<<"]R_QCD: "<<(*p_rprocqcd)[i]->Name()<<std::endl;
      for (size_t i(0);i<p_rsprocqcd->Size();++i)
        msg_Debugging()<<"["<<i<<"]RS_QCD:  "<<(*p_rsprocqcd)[i]->Name()<<std::endl;
      THROW(fatal_error,"QCD R and RS have different size");
    }
    if (p_bproc->Size()!=p_viprocqcd->Size()) {
      for (size_t i(0);i<p_bproc->Size();++i)
        msg_Debugging()<<"["<<i<<"]B:   "<<(*p_bproc)[i]->Name()<<std::endl;
      for (size_t i(0);i<p_viprocqcd->Size();++i)
        msg_Debugging()<<"["<<i<<"]VI_QCD: "<<(*p_viprocqcd)[i]->Name()<<std::endl;
      THROW(fatal_error,"B and VI_QCD have different size");
    }
    for (size_t i(0);i<p_rsprocqcd->Size();++i)
      if ((*p_rsprocqcd)[i]->Flavours()!=(*p_rprocqcd)[i]->Flavours())
        THROW(fatal_error,"Ordering differs in QCD R and RS");
    if (!p_read) for (size_t i(0);i<p_viprocqcd->Size();++i)
      if ((*p_viprocqcd)[i]->Flavours()!=(*p_bproc)[i]->Flavours())
        THROW(fatal_error,"Ordering differs in B and VI_QCD");
  }
  if (m_ew) {
    if (p_rsprocqed->Size()!=p_rprocqed->Size()) {
      for (size_t i(0);i<p_rprocqed->Size();++i)
        msg_Debugging()<<"["<<i<<"]R_QED: "<<(*p_rprocqed)[i]->Name()<<std::endl;
      for (size_t i(0);i<p_rsprocqed->Size();++i)
        msg_Debugging()<<"["<<i<<"]RS_QED:  "<<(*p_rsprocqed)[i]->Name()<<std::endl;
      THROW(fatal_error,"QED R and RS have different size");
    }
    if (p_bproc->Size()!=p_viprocew->Size()) {
      for (size_t i(0);i<p_bproc->Size();++i)
        msg_Debugging()<<"["<<i<<"]B:   "<<(*p_bproc)[i]->Name()<<std::endl;
      for (size_t i(0);i<p_viprocqcd->Size();++i)
        msg_Debugging()<<"["<<i<<"]VI_EW: "<<(*p_viprocew)[i]->Name()<<std::endl;
      THROW(fatal_error,"B and VI_EW have different size");
    }
    for (size_t i(0);i<p_rsprocqed->Size();++i)
      if ((*p_rsprocqed)[i]->Flavours()!=(*p_rprocqed)[i]->Flavours())
        THROW(fatal_error,"Ordering differs in QED R and RS");
    for (size_t i(0);i<p_viprocew->Size();++i)
      if ((*p_viprocew)[i]->Flavours()!=(*p_bproc)[i]->Flavours())
        THROW(fatal_error,"Ordering differs in B and VI_EW");
  }

  // set the correct type for the reweighting
  if (m_qcd) for (size_t i(0);i<p_rsprocqcd->Size();++i)
    (*p_rsprocqcd)[i]->GetMEwgtinfo()->m_type=mewgttype::H;
  if (m_ew) for (size_t i(0);i<p_rsprocqed->Size();++i)
    (*p_rsprocqed)[i]->GetMEwgtinfo()->m_type=mewgttype::H;
}

void MCatNLO_Process::RegisterDefaults() const
{
  Scoped_Settings s{ Settings::GetMainSettings()["MC@NLO"] };
  s["HPSMODE"].SetDefault(-1);  // H event shower mode
  s["KFACTOR_MODE"].SetDefault(14);  // K-factor mode
  s["FOMODE"].SetDefault(0);  // fixed order mode
  s["RS_SCALE"].SetDefault("");  // RS scale
  s["PP_BVI_MODE"].SetDefault(7);  // BVI mode
}

Process_Base* MCatNLO_Process::InitProcess
(const Process_Info &pi,nlo_type::code nlotype,const int real)
{
  Process_Info cpi(pi);
  cpi.m_fi.SetNLOType(nlotype);
  if (real) {
    if (m_ew)
      cpi.m_fi.m_ps.push_back(Subprocess_Info(kf_ewjet,"",""));
    else if (m_qcd)
      cpi.m_fi.m_ps.push_back(Subprocess_Info(kf_jet,"",""));
    else THROW(fatal_error, "Internal error.");
  }
  Process_Base* proc = m_gens.InitializeProcess(cpi,false);
  return proc;
}

bool MCatNLO_Process::InitSubtermInfo()
{
  if (m_qcd) for (size_t i(0);i<p_ddprocqcd->Size();++i) {
    NLO_subevtlist *subs((*p_ddprocqcd)[i]->GetSubevtList());
    for (size_t j(0);j<subs->size()-1;++j) {
      NLO_subevt *sub((*subs)[j]);
      for (size_t ij(0);ij<sub->m_n;++ij)
	      for (size_t k(0);k<sub->m_n;++k)
          if (k != ij && sub->p_fl[k] == sub->p_fl[sub->m_kt] &&
              sub->p_fl[ij]==sub->p_fl[sub->m_ijt]) {
	    m_iinfo[sub->m_pname].insert(IDip_ID(ij,k));
	  }
      m_dinfo[subs->back()->m_pname].insert(*sub);
    }
  }
  if (m_ew) for (size_t i(0);i<p_ddprocqed->Size();++i) {
    NLO_subevtlist *subs((*p_ddprocqed)[i]->GetSubevtList());
    for (size_t j(0);j<subs->size()-1;++j) {
      NLO_subevt *sub((*subs)[j]);
      for (size_t ij(0);ij<sub->m_n;++ij)
	      for (size_t k(0);k<sub->m_n;++k)
          if (k != ij && sub->p_fl[k] == sub->p_fl[sub->m_kt] &&
              sub->p_fl[ij]==sub->p_fl[sub->m_ijt]) {
	    m_iinfo[sub->m_pname].insert(IDip_ID(ij,k));
	  }
      m_dinfo[subs->back()->m_pname].insert(*sub);
    }
  }
  return true;
}

Process_Base *MCatNLO_Process::FindProcess
(const Cluster_Amplitude *ampl,const nlo_type::code type,
 const bool error) const
{
  std::string name(Process_Base::GenerateName(ampl));
  StringProcess_Map::const_iterator pit(p_apmap->find(type)->second->find(name));
  if (pit!=p_apmap->find(type)->second->end()) {
    if (type==nlo_type::lo) {
      for (size_t i(0);i<p_bproc->Size();++i)
	      if ((*p_bproc)[i]==pit->second) return pit->second;
      for (size_t i(0);i<p_bviproc->Size();++i)
	      if ((*p_bviproc)[i]==pit->second) return pit->second;
      for (size_t i(0);i<p_rprocqcd->Size();++i)
	      if ((*p_rprocqcd)[i]==pit->second) return pit->second;
    }
    if (type==(nlo_type::lo|nlo_type::qed)) {
      for (size_t i(0);i<p_rprocqed->Size();++i)
	      if ((*p_rprocqed)[i]==pit->second) return pit->second;
    }
    if (type&nlo_type::vsub) {
      if (p_read) for (size_t i(0);i<p_bviproc->Size();++i)
	      if ((*p_bviproc)[i]==pit->second) return pit->second;
      if (!p_read) for (size_t i(0);i<p_viprocqcd->Size();++i)
	      if ((*p_viprocqcd)[i]==pit->second) return pit->second;
      for (size_t i(0);i<p_viprocew->Size();++i)
	      if ((*p_viprocew)[i]==pit->second) return pit->second;
    }
    if (type&nlo_type::rsub) {
      for (size_t i(0);i<p_rsprocqcd->Size();++i)
	      if ((*p_rsprocqcd)[i]==pit->second) return pit->second;
      for (size_t i(0);i<p_rsprocqed->Size();++i)
	      if ((*p_rsprocqed)[i]==pit->second) return pit->second;
    }
  }
  if (error)
    THROW(fatal_error,"Process '"+name+"'("+ToString(type)+") not found");
  return NULL;
}

Cluster_Amplitude *MCatNLO_Process::CreateAmplitude(const NLO_subevt *sub) const
{
  Cluster_Amplitude *ampl = Cluster_Amplitude::New();
  ampl->SetNIn(m_nin);
  ampl->SetMS(m_qcd?p_rsprocqcd->Generator():p_rsprocqed->Generator());
  ampl->SetMuF2(sub->m_mu2[stp::fac]);
  ampl->SetMuR2(sub->m_mu2[stp::ren]);
  Int_Vector ci(sub->m_n,0), cj(sub->m_n,0);
  for (size_t i=0;i<sub->m_n;++i) {
    ampl->CreateLeg(i<m_nin?-sub->p_mom[i]:sub->p_mom[i],
		    i<m_nin?sub->p_fl[i].Bar():sub->p_fl[i],
		    ColorID(ci[i],cj[i]),sub->p_id[i]);
    if (!sub->IsReal() && sub->p_id[i]&(1<<sub->m_i)) {
      if ((sub->p_id[i]&(1<<sub->m_j))==0)
	THROW(fatal_error,"Internal error");
      ampl->Legs().back()->SetK(1<<sub->m_k);
    }
  }
  ampl->Decays()=*sub->p_dec;
  return ampl;
}

Weights_Map MCatNLO_Process::Differential(const Vec4D_Vector &p,
					  Variations_Mode varmode)
{
  THROW(fatal_error,"Invalid function call");
}

Weights_Map MCatNLO_Process::LocalKFactor(Cluster_Amplitude& ampl)
{
  DEBUG_FUNC(Name());

  // validate ampl
  Cluster_Amplitude *rampl(ampl.Prev());
  if (rampl->Legs().size()!=m_nin+m_nout) return 0.0;

  // evaluate RS process
  msg_Debugging()<<*rampl<<"\n";
  Process_Base *rsproc(FindProcess(rampl,nlo_type::rsub,false));
  if (rsproc==NULL) return 0.0;
  rsproc->SetEventReader(NULL);
  msg_Debugging()<<"Found '"<<rsproc->Name()<<"'\n";
  Cluster_Amplitude *crampl(rampl->Copy());
  int mm(p_bproc->Generator()->SetMassMode(0));
  p_bproc->Generator()->ShiftMasses(crampl);
  p_bproc->Generator()->SetMassMode(mm);
  int rmode = rampl->ColorMap().empty() ? 128 : 0;
  Weights_Map rs = rsproc->Differential(*crampl, Variations_Mode::all, rmode);
  Weights_Map r = rsproc->GetSubevtList()->back()->m_results;
  if (rmode && r.Nominal() == 0.0 &&
      rsproc->Differential(*crampl, Variations_Mode::nominal_only, 64).Nominal() != 0.0) {
    for (int i(0); i < 100 && r.Nominal() == 0.0; ++i) {
      rs = rsproc->Differential(*crampl, Variations_Mode::all, rmode);
      r = rsproc->GetSubevtList()->back()->m_results;
    }
  }
  crampl->Delete();
  msg_Debugging() << "H = " << rs.Nominal() << ", R = " << r.Nominal() << " -> "
                  << rs.Nominal() / r.Nominal() << "\n";
  if (r.Nominal() == 0.0) {
    return 0.0;
  }

  // evaluate BVI process
  msg_Debugging()<<ampl<<"\n";
  Process_Base *bviproc(FindProcess(&ampl,nlo_type::vsub,false));
  if (bviproc==NULL) return 0.0;
  bviproc->SetEventReader(NULL);
  msg_Debugging()<<"Found '"<<bviproc->Name()<<"'\n";
  Process_Base *bproc(FindProcess(&ampl));
  Cluster_Amplitude *campl(ampl.Copy());
  p_bproc->Generator()->SetMassMode(0);
  p_bproc->Generator()->ShiftMasses(campl);
  p_bproc->Generator()->SetMassMode(mm);
  bproc->GetMEwgtinfo()->m_type = mewgttype::none;
  int bmode = ampl.ColorMap().empty() ? 128 : 0;
  Weights_Map b = bproc->Differential(*campl, Variations_Mode::all, bmode);
  if (bmode && b.Nominal() == 0.0 &&
      bproc->Differential(*campl, Variations_Mode::nominal_only, 64).Nominal() != 0.0) {
    for (int i(0); i < 100 && b.Nominal() == 0.0; ++i) {
      b = bproc->Differential(*campl, Variations_Mode::all, bmode);
    }
  }
  if (b.Nominal() == 0.0) {
    campl->Delete();
    return 0.0;
  }
  bviproc->BBarMC()->GenerateEmissionPoint(*campl);
  Weights_Map bvi = bviproc->Differential(*campl, Variations_Mode::all, bmode);
  campl->Delete();

  // eventually calculate local K factor
  msg_Debugging() << "Calc'ing LocalKFactor weights.\n";
  const double random(ran->Get());

  // Below, we'll Iterate through all weight tuples, calculating the
  // corresponding K factors. For this, we use the absolute version of the
  // variations, which makes it more straightforward to deal with each
  // variation individually in the loop below.
  bvi.MakeAbsolute();
  b.MakeAbsolute();
  rs.MakeAbsolute();
  r.MakeAbsolute();

  // Copy BVI variations to get the same variations structure; the actual K
  // factor data is set in the loop below; note that via bvi.MakeAbsolute(), we
  // can be certain that global factors such as base_weight and
  // nominal_prefactor are set to 1.0
  Weights_Map kfacs = bvi;

  for (const auto& kv : bvi) {
    const std::string& key = kv.first;
    const Weights& weights = kv.second;
    const auto num_weights = weights.Size();
    for (int i {0}; i < num_weights; i++) {
      msg_Debugging() << weights.Name(i) << '\n';

      const LocalKFactorInfo info = CalculateLocalKFactorInfo(
          bvi.Get(key, i), b.Get(key, i), rs.Get(key, i), r.Get(key, i));
      if (info.s == 0.0 && info.h == 0.0) {
        kfacs[key][i] = 0.0;
        continue;
      }
      // select S or H and set corresponding K factor; update cluster
      // amplitudes if this is the nominal calculation
      const double selectionwgt {dabs(info.s) / (dabs(info.s) + dabs(info.h))};
      if (selectionwgt > random) {
        const double kfac {info.s / selectionwgt};
        msg_Debugging() << "S selected ( w = " << kfac << " )\n";
        if (key == "Main" && i == 0 && m_kfacmode / 10) {
          for (Cluster_Amplitude* campl(ampl.Next()); campl;
              campl = campl->Next()) {
            campl->SetLKF(bvi[key][i] / b[key][i]);
            campl->SetNLO(2);
          }
        }
        kfacs[key][i] = kfac;
      } else {
        const double kfac {info.h / (1.0 - selectionwgt)};
        msg_Debugging() << "H selected ( w = " << kfac << " )\n";
        if (key == "Main" && i == 0 && m_kfacmode / 10)
          ampl.SetNLO(m_hpsmode);
        kfacs[key][i] = kfac;
      }
      for (Cluster_Amplitude *campl(&ampl);campl;
	   campl=campl->Next()) msg_Debugging()<<*campl<<"\n";
    }
  }
  kfacs.MakeRelative();
  return kfacs;
}

MCatNLO_Process::LocalKFactorInfo MCatNLO_Process::CalculateLocalKFactorInfo(
    double bvi, double b, double rs, double r)
{
  LocalKFactorInfo info;
  const double bvib(b ? bvi / b : 0.0), rsr(r ? rs / r : 0.0);
  if (m_kfacmode % 10 == 0) {
    info.s = bvib * (1.0 - rsr);
    info.h = rsr;
  } else if (m_kfacmode % 10 == 1) {
    info.s = bvib * (1.0 - rsr);
    info.h = 0;
  } else if (m_kfacmode % 10 == 2) {
    info.s = 0;
    info.h = rsr;
  } else if (m_kfacmode % 10 == 3) {
    info.s = bvib;
    info.h = 0.;
  } else if (m_kfacmode % 10 == 4) {
    info.s = bvib + (b?rs/b:0.);
    info.h = 0.;
  } else
    THROW(fatal_error, "Unknown Kfactor mode.");
  msg_Debugging() << "BVI = " << bvi << ", B = " << b << " -> S = " << info.s
                  << ", H = " << info.h << "\n";
  return info;
}

Cluster_Amplitude *MCatNLO_Process::GetAmplitude()
{
  if (p_ampl==NULL) return NULL;
  Cluster_Amplitude *ampl(p_ampl->CopyAll());
  ME_Generator_Base *gen((ME_Generator_Base*)ampl->MS());
  int mm(gen->SetMassMode(1));
  int stat(gen->ShiftMasses(ampl));
  if (stat<0) {
    msg_Tracking() << METHOD << "(): Mass shift failed." << std::endl;
    gen->SetMassMode(mm);
    return NULL;
  }
  if (stat==1) {
    stat=p_shower->GetClusterDefinitions()->ReCluster(ampl);
    if (stat!=1) {
      msg_Debugging()<<METHOD<<"(): Reclustering failed."<<std::endl;
    }
  }
  gen->SetMassMode(mm);
  return ampl;
}

Weights_Map MCatNLO_Process::OneHEvent(const int wmode)
{
  msg_Debugging()<<"H event\n";
  m_wassevent=false;
  Process_Base *rsproc(NULL), *_rproc(NULL);
  if (m_evttype==1) {
    p_selected=p_rprocqcd;
    _rproc=p_rprocqcd;
    rsproc=p_rsprocqcd;
  }
  if (m_evttype==2) {
    p_selected=p_rprocqed;
    _rproc=p_rprocqed;
    rsproc=p_rsprocqed;
  }
  if (p_ampl) p_ampl->Delete();
  const size_t selectedindex(_rproc->SynchronizeSelectedIndex(*rsproc));
  assert(selectedindex != std::numeric_limits<size_t>::max());
  Process_Base *rproc((*_rproc)[selectedindex]);
  rproc->SetLastXS(rsproc->Selected()->LastXS());
  rproc->Integrator()->SetMax
    (rsproc->Selected()->Integrator()->Max());
  Vec4D_Vector &p(rsproc->Selected()->Integrator()->Momenta());
  rproc->SetFixedScale(rsproc->Selected()->ScaleSetter(1)->Scales());
  rproc->ScaleSetter(1)->CalculateScale(Vec4D_Vector());
  rproc->SetFixedScale(std::vector<double>());
  rproc->GetMEwgtinfo()->m_mur2=
    rsproc->Selected()->GetMEwgtinfo()->m_mur2;
  rproc->GetMEwgtinfo()->m_muf2=
    rsproc->Selected()->GetMEwgtinfo()->m_muf2;
  rproc->Integrator()->SetMomenta(p);
  Color_Integrator *ci(rproc->Integrator()->ColorIntegrator().get()),
    *rci(rsproc->Selected()->Integrator()->ColorIntegrator().get());
  if (ci && rci) {
    ci->SetI(rci->I());
    ci->SetJ(rci->J());
  }
  p_ampl=rsproc->Selected()->GetSubevtList()->back()->p_ampl;
  if (p_ampl) p_ampl=p_ampl->CopyAll();
  else p_ampl=((Single_Process*)(rproc))->Cluster(p,1);
  p_selected->Selected()->SetMEwgtinfo(*rsproc->Selected()->GetMEwgtinfo());
  if (p_ampl==NULL) {
    msg_Error()<<METHOD<<"(): No valid clustering. Skip event."<<std::endl;
    return {0};
  }
  Scale_Setter_Base *scs(rsproc->Selected()->ScaleSetter(1));
  for (Cluster_Amplitude *ampl(p_ampl);ampl;ampl=ampl->Next()) {
    ampl->SetMuR2(scs->Scale(stp::ren));
    ampl->SetMuF2(scs->Scale(stp::fac));
    ampl->SetMuQ2(scs->Scale(stp::res));
  }
  if (p_ampl->Next()) {
    p_ampl->Next()->SetNLO
      (p_ampl->Next()->NLO()|m_hpsmode);
  }
  Selector_Base *jf=rsproc->Selected()->
    Selector()->GetSelector("Jetfinder");
  Weights_Map wgtmap;
  wgtmap.insert({"QCUT", Weights {Variations_Type::qcut}});
  if (jf && m_nout-1<m_pinfo.m_fi.NMaxExternal()) {
    for (Cluster_Amplitude *ampl(p_ampl);
	 ampl;ampl=ampl->Next()) ampl->SetJF(jf);
    Jet_Finder *cjf(static_cast<Jet_Finder*>(jf));
    HEventVeto_Args hva(cjf,cjf->JC()->Value(p_ampl));
    for (size_t i {0}; i < s_variations->Size(Variations_Type::qcut) + 1;
         ++i) {
      const double fac {
          i != 0 ? s_variations->Qcut_Parameters(i - 1).m_scale_factor : 1.0};
      int stat = hva.m_jcv < sqr(hva.p_jf->Qcut() * fac);
      wgtmap["QCUT"][i] = (stat ? 1.0 : 0.0);
    }
  }
  msg_Debugging()<<"R selected via NLO "<<*p_ampl<<"\n";
  for (Cluster_Amplitude *campl(p_ampl->Next());campl;
       campl=campl->Next()) msg_Debugging()<<*campl<<"\n";
  return wgtmap;
}

Weights_Map MCatNLO_Process::OneSEvent(const int wmode)
{
  msg_Debugging()<<"S event\n";
  m_wassevent=true;
  Process_Base *bviproc(NULL);
  if (m_evttype==0) bviproc=p_bviproc;
  if (m_evttype==1) bviproc=p_viprocqcd;
  if (m_evttype==2) bviproc=p_viprocew;
  if (p_ampl) p_ampl->Delete();
  p_ampl=NULL;
  const size_t selectedindex(p_bproc->SynchronizeSelectedIndex(*bviproc));
  assert(selectedindex != std::numeric_limits<size_t>::max());
  Process_Base *bproc((*p_bproc)[selectedindex]);
  bproc->Integrator()->SetMax
    (bviproc->Selected()->Integrator()->Max());
  bproc->Integrator()->SetEnhanceFactor
    (bviproc->Selected()->Integrator()->EnhanceFactor());
  bproc->SetSubEvt(bviproc->Selected()->SubEvt());
  if (bproc->SubEvt()) bproc->SubEvt()->p_proc=bviproc->Selected();
  bproc->SetLastXS(bviproc->Selected()->LastXS());
  Vec4D_Vector &p(bviproc->Selected()->Integrator()->Momenta());
  p_ampl = dynamic_cast<Single_Process*>
    (bviproc->Selected())->Cluster(p);
  SortFlavours(p_ampl);
  p_ampl->SetProcs(p_apmap);
  p_ampl->SetIInfo(&m_iinfo);
  p_ampl->SetDInfo(&m_dinfo);
  p_ampl->Decays()=m_decins;
  if (p_ampl->JF<void>())
    for (Cluster_Amplitude *sampl(p_ampl);
	 sampl;sampl=sampl->Next()) sampl->SetNLO(2);
  double lkf(bviproc->Selected()->Last()/bviproc->Selected()->LastBValue());
  for (Cluster_Amplitude *ampl(p_ampl);
       ampl;ampl=ampl->Next()) ampl->SetLKF(lkf);
  int stat(1);
  Weights_Map mcwgt(1);
  if (!(m_psmode&2)) {
    p_nlomc->SetShower(p_shower);
    if (!m_fomode) {
      stat=p_nlomc->GeneratePoint(p_ampl);
      mcwgt=p_nlomc->WeightsMap();
    }
  }
  Cluster_Amplitude *next(p_ampl), *ampl(p_ampl->Prev());
  if (ampl) {
    p_ampl=NULL;
    Process_Base *rproc(NULL);
    Flavour_Vector afl(ampl->Legs().size());
    for (size_t i(0);i<afl.size();++i) {
      afl[i]=i<m_nin?ampl->Leg(i)->Flav().Bar():ampl->Leg(i)->Flav();
    }
    if (m_qcd) for (size_t i(0);i<p_rprocqcd->Size();++i)
      if (afl==(*p_rprocqcd)[i]->Flavours()) {
        rproc=(*p_rprocqcd)[i];
        p_selected=p_rprocqcd;
        p_rprocqcd->SetSelected(rproc);
        break;
      }
    if (m_ew && !rproc) for (size_t i(0);i<p_rprocqed->Size();++i)
      if (afl==(*p_rprocqed)[i]->Flavours()) {
        rproc=(*p_rprocqed)[i];
        p_selected=p_rprocqed;
        p_rprocqed->SetSelected(rproc);
        break;
      }
    if (rproc==NULL) THROW(fatal_error,"Invalid splitting");
    rproc->Integrator()->PSHandler()->SetEnhanceWeight(
        bproc->Integrator()->PSHandler()->EnhanceWeight());
    rproc->Integrator()->SetMax(bproc->Integrator()->Max());
    rproc->Integrator()->SetEnhanceFactor(
        bproc->Selected()->Integrator()->EnhanceFactor());
    rproc->Integrator()->SetMomenta(*ampl);
    rproc->GetMEwgtinfo()->m_mur2=bproc->GetMEwgtinfo()->m_mur2;
    Cluster_Leg *lij(NULL);
    for (size_t i(0);i<next->Legs().size();++i) {
      next->Leg(i)->SetStat(1);
      if (next->Leg(i)->K()) {
	lij=next->Leg(i);
      }
    }
    if (!lij) THROW(fatal_error,"Internal error");
    std::vector<int> ids(ID(lij->Id()));
    size_t iid(1<<ids.front()), jid(1<<ids.back()), kid(lij->K());
    ids.push_back(ID(lij->K()).front());
    if (ids.size()!=3) THROW(fatal_error,"Internal error");
    Cluster_Amplitude *campl(ampl->Copy());
    campl->IdSort();
    Cluster_Definitions_Base *clus(p_shower->GetClusterDefinitions());
    Cluster_Param kt2=clus->Cluster
      (Cluster_Config(campl,ids[0],ids[1],ids[2],lij->Flav(),
		      bviproc->Generator(),NULL,1));
    if (kt2.m_kt2<=0.0) {
      kt2=clus->Cluster
	(Cluster_Config(campl,ids[1],ids[0],ids[2],lij->Flav(),
			bviproc->Generator(),NULL,1));
    }
    campl->Delete();
    ampl->SetKT2(kt2.m_kt2);
    ampl->SetMu2(kt2.m_kt2);
    p_ampl=ampl;
    ampl->SetMuF2(next->MuF2());
    ampl->SetMuR2(next->MuR2());
    ampl->SetOrderQCD(next->OrderQCD()+1);
    ampl->Next()->SetNLO(4);
    ampl->SetJF(ampl->Next()->JF<void>());
    ampl->Next()->SetCA(clus);
    next->SetKin(kt2.m_kin);
    while (ampl->Next()) {
      ampl=ampl->Next();
      if (!(ampl->Leg(0)->Id()&p_ampl->Leg(0)->Id()))
	std::swap<Cluster_Leg*>(ampl->Legs()[0],ampl->Legs()[1]);
      if (ampl->IdNew()&iid) ampl->SetIdNew(ampl->IdNew()|jid);
      for (size_t i(0);i<ampl->Legs().size();++i) {
	ampl->Leg(i)->SetNMax(p_ampl->Leg(i)->NMax());
	size_t cid(ampl->Leg(i)->Id());
	if (cid&iid) {
	  for (size_t j(0);j<ampl->Legs().size();++j)
	    if (ampl->Leg(j)->K()==cid)
	      ampl->Leg(j)->SetK(cid|jid);
	  ampl->Leg(i)->SetId(cid|jid);
	  if (ampl->Prev()->Prev()==NULL) {
	    ampl->Leg(i)->SetK(kid);
	  }
	}
      }
    }
    msg_Debugging()<<"R selected via Sudakov "<<*p_ampl
		   <<" ( w = "<<mcwgt.Nominal()<<" )\n";
    p_selected->Selected()->SetMEwgtinfo(*bviproc->Selected()->GetMEwgtinfo());
    return mcwgt;
  }
  p_selected=p_bproc;
  ampl=p_ampl;
  if (!(m_psmode&2)) ampl->SetNLO(4);
  else ampl->SetNLO(0);
  bproc->Integrator()->SetMomenta(*p_ampl);
  msg_Debugging()<<"B selected "<<*p_ampl
		 <<" ( w = "<<mcwgt.Nominal()<<" )\n";
  for (Cluster_Amplitude *campl(p_ampl->Next());campl;
       campl=campl->Next()) msg_Debugging()<<*campl<<"\n";
  bproc->SetMEwgtinfo(*bviproc->Selected()->GetMEwgtinfo());
  if (m_psmode&2) return Weights_Map{1.0};
  return stat ? mcwgt : Weights_Map{0.0};
}

Weight_Info *MCatNLO_Process::OneEvent(const int wmode,
                                       Variations_Mode varmode,
                                       const int mode)
{
  DEBUG_FUNC(m_name<<" "<<p_read);
  double SB(0), SQCD(0), SEW(0), HQCD(0), HQED(0);
  Weight_Info *winfo(NULL);
  if (p_read==NULL) {
    SB=p_bviproc->Integrator()->SelectionWeight(wmode);
    if (m_qcd) {
      SQCD=p_viprocqcd->Integrator()->SelectionWeight(wmode);
      HQCD=p_rsprocqcd->Integrator()->SelectionWeight(wmode);
    }
    if (m_ew) {
      SEW=p_viprocew->Integrator()->SelectionWeight(wmode);
      HQED=p_rsprocqed->Integrator()->SelectionWeight(wmode);
    }
  }
  else {
    msg_Debugging()<<"Reader type is "<<Demangle(typeid(*p_read).name())<<"\n";
    Cluster_Amplitude *ampl(p_read->ReadEvent());
    if (ampl==NULL) return NULL;
    msg_Debugging()<<*ampl<<"\n";
    if (p_read->SubEvt()->m_n) {
      msg_Debugging()<<"S event\n";
      p_read->SetAmpl(ampl);
      p_bviproc->SetEventReader(p_read);
      winfo = p_bviproc->OneEvent(wmode);
      p_bviproc->SetEventReader(NULL);
      SB=1;
    }
    else {
      msg_Debugging()<<"H event\n";
      p_read->SetAmpl(ampl);
      p_rsprocqcd->SetEventReader(p_read);
      winfo = p_rsprocqcd->OneEvent(wmode);
      p_rsprocqcd->SetEventReader(NULL);
      HQCD=1;
    }
    if (winfo==NULL) {
      msg_Debugging()<<"No weight info\n";
      return winfo;
    }
    msg_Debugging()<<"Weight info: "<<*winfo<<"\n";
  }
  // msg_Debugging()<<METHOD<<" qcd: "<<m_qcd<<", ew: "<<m_ew<<", SB: "<<SB<<", "
  const double r(ran->Get() * (SB + SQCD + SEW + HQCD + HQED));
  if (SB > r) {
    msg_Debugging()<<METHOD<<" SB\n";
    p_selected = p_bviproc;
    m_evttype=0;
    if (p_read==NULL) winfo = p_bviproc->OneEvent(wmode, varmode, mode);
    if (winfo && m_fomode != 2) {
      // calculate and apply weight factor
      const Weights_Map Swgts {OneSEvent(wmode)};
      assert(p_ampl);
      msg_Debugging()<<METHOD<<" B:\n";
      const double Bsel(p_bproc->Selected()->Integrator()->SelectionWeight(wmode));
      msg_Debugging()<<METHOD<<" BVI:\n";
      const double Ssel(p_bviproc->Selected()->Integrator()->SelectionWeight(wmode));
      const double selwgtratio(Bsel / Ssel);
      const double wgtfac(selwgtratio);
      winfo->m_weightsmap *= wgtfac;
      winfo->m_weightsmap *= Swgts;
      *(p_selected->Selected()->GetMEwgtinfo()) *= Swgts.Nominal()*wgtfac;
    }
  } else if (SB + SQCD > r) {
    msg_Debugging()<<METHOD<<" SQCD\n";
    p_selected = p_viprocqcd;
    m_evttype=1;
    if (p_read==NULL) winfo = p_viprocqcd->OneEvent(wmode, varmode, mode);
    if (winfo && m_fomode != 2) {
      // calculate and apply weight factor
      const Weights_Map Swgts {OneSEvent(wmode)};
      assert(p_ampl);
      msg_Debugging()<<METHOD<<" B:\n";
      const double Bsel(p_bproc->Selected()->Integrator()->SelectionWeight(wmode));
      msg_Debugging()<<METHOD<<" VI QCD:\n";
      const double Ssel(p_viprocqcd->Selected()->Integrator()->SelectionWeight(wmode));
      const double selwgtratio(Bsel / Ssel);
      const double wgtfac(selwgtratio);
      winfo->m_weightsmap *= wgtfac;
      winfo->m_weightsmap *= Swgts;
      *(p_selected->Selected()->GetMEwgtinfo()) *= Swgts.Nominal()*wgtfac;
    }
  } else if (SB + SQCD + SEW > r) {
    msg_Debugging()<<METHOD<<" SEW\n";
    p_selected = p_viprocew;
    m_evttype=2;
    if (p_read==NULL) winfo = p_viprocew->OneEvent(wmode, varmode, mode);
    if (winfo && m_fomode != 2) {
      // calculate and apply weight factor
      const Weights_Map Swgts {OneSEvent(wmode)};
      assert(p_ampl);
      msg_Debugging()<<METHOD<<" B:\n";
      const double Bsel(p_bproc->Selected()->Integrator()->SelectionWeight(wmode));
      msg_Debugging()<<METHOD<<" VI EW:\n";
      const double Ssel(p_viprocew->Selected()->Integrator()->SelectionWeight(wmode));
      const double selwgtratio(Bsel / Ssel);
      const double wgtfac(selwgtratio);
      winfo->m_weightsmap *= wgtfac;
      winfo->m_weightsmap *= Swgts;
      *(p_selected->Selected()->GetMEwgtinfo()) *= Swgts.Nominal()*wgtfac;
    }
  } else if (SB + SQCD + SEW + HQCD > r) {
    msg_Debugging()<<METHOD<<" HQCD\n";
    p_selected = p_rsprocqcd;
    m_evttype=1;
    if (p_read==NULL) winfo = p_rsprocqcd->OneEvent(wmode, varmode, mode);
    if (winfo && m_fomode != 2) {
      // calculate and apply weight factor
      const Weights_Map Hwgts {OneHEvent(wmode)};
      msg_Debugging()<<METHOD<<" R QCD:\n";
      const double Rsel(p_rprocqcd->Selected()->Integrator()->SelectionWeight(wmode));
      msg_Debugging()<<METHOD<<" RS QCD:\n";
      const double RSsel(p_rsprocqcd->Selected()->Integrator()->SelectionWeight(wmode));
      const double selwgtratio(Rsel / RSsel);
      const double wgtfac(selwgtratio);
      winfo->m_weightsmap *= wgtfac;
      winfo->m_weightsmap *= Hwgts;
      *p_selected->Selected()->GetMEwgtinfo() *= Hwgts.Nominal()*wgtfac;
    }
  } else {
    msg_Debugging()<<METHOD<<" HQED\n";
    p_selected = p_rsprocqed;
    m_evttype=2;
    if (p_read==NULL) winfo = p_rsprocqed->OneEvent(wmode, varmode, mode);
    if (winfo && m_fomode != 2) {
      // calculate and apply weight factor
      const Weights_Map Hwgts {OneHEvent(wmode)};
      const double Rsel(p_rprocqed->Selected()->Integrator()->SelectionWeight(wmode));
      const double RSsel(p_rsprocqed->Selected()->Integrator()->SelectionWeight(wmode));
      const double selwgtratio(Rsel / RSsel);
      const double wgtfac(selwgtratio);
      winfo->m_weightsmap *= wgtfac;
      winfo->m_weightsmap *= Hwgts;
      *p_selected->Selected()->GetMEwgtinfo() *= Hwgts.Nominal()*wgtfac;
    }
  }
  Mass_Selector *ms(Selected()->Generator());
  for (Cluster_Amplitude *ampl(p_ampl);
       ampl;ampl=ampl->Next()) {
    ampl->SetNLO(1|ampl->NLO());
    ampl->SetMS(ms);
  }
  if (winfo && winfo->m_weightsmap.IsZero()) {
    delete winfo;
    winfo=NULL;
  }
  if (winfo==NULL) return winfo;

  if (rpa->gen.HardSC() || (rpa->gen.SoftSC() && !Flavour(kf_tau).IsStable())) {
    msg_Debugging()<<"Calcing Differential for spin correlations using "
                   <<Selected()->Generator()->Name()<<":"<<std::endl;
    ME_Weight_Info mwi(*Selected()->GetMEwgtinfo());
    if (Selected()->Integrator()->ColorIntegrator()!=NULL)
      while (Selected()->Differential(*p_ampl,Variations_Mode::nominal_only,1|2|4|128).Nominal()==0.0);
    else
      Selected()->Differential(*p_ampl,Variations_Mode::nominal_only,1|2|4|128);
    Selected()->SetMEwgtinfo(mwi);
  }
  return winfo;
}

void MCatNLO_Process::InitPSHandler
(const double &maxerror,const std::string eobs,const std::string efunc)
{
  p_bviproc->InitPSHandler(maxerror,eobs,efunc);
  p_bproc->Integrator()->SetPSHandler
    (p_bviproc->Integrator()->PSHandler());
  if (m_qcd) {
    p_viprocqcd->InitPSHandler(maxerror,eobs,efunc);
    p_ddprocqcd->InitPSHandler(maxerror,eobs,efunc);
    p_rsprocqcd->InitPSHandler(maxerror,eobs,efunc);
    p_rsprocqcd->Integrator()->SetEnhanceFactor
      (p_rsprocqcd->Integrator()->EnhanceFactor()*p_int->RSEnhanceFactor());
    p_rprocqcd->Integrator()->SetPSHandler
      (p_rsprocqcd->Integrator()->PSHandler());
  }
  if (m_ew) {
    p_viprocew->InitPSHandler(maxerror,eobs,efunc);
    p_ddprocqed->InitPSHandler(maxerror,eobs,efunc);
    p_rsprocqed->InitPSHandler(maxerror,eobs,efunc);
    p_rsprocqed->Integrator()->SetEnhanceFactor
      (p_rsprocqed->Integrator()->EnhanceFactor()*p_int->RSEnhanceFactor());
    p_rprocqed->Integrator()->SetPSHandler
      (p_rsprocqed->Integrator()->PSHandler());
  }
}

bool MCatNLO_Process::CalculateTotalXSec(const std::string &resultpath,
                                         const bool         create)
{
  Vec4D_Vector p(m_qcd?p_rsprocqcd->NIn()+p_rsprocqcd->NOut():
                       p_rsprocqed->NIn()+p_rsprocqed->NOut());
  Cluster_Amplitude *ampl(Cluster_Amplitude::New());
  for (int i(0);i<p.size();++i) {
    ampl->CreateLeg(Vec4D(),Flavour(kf_jet));
  }
  PHASIC::Phase_Space_Handler *psh(NULL);
  do {
    if (m_qcd) {
      psh = p_ddprocqcd->Integrator()->PSHandler();
      psh->TestPoint(&p.front(),&p_ddprocqcd->Info(),p_ddprocqcd->Generator());
      for (size_t i(0);i<p.size();++i) ampl->Leg(i)->SetMom(p[i]);
      p_ddprocqcd->Differential(*ampl,Variations_Mode::nominal_only,4);
    }
    if (m_ew) {
      psh = p_ddprocqed->Integrator()->PSHandler();
      psh->TestPoint(&p.front(),&p_ddprocqed->Info(),p_ddprocqed->Generator());
      for (size_t i(0);i<p.size();++i) ampl->Leg(i)->SetMom(p[i]);
      p_ddprocqed->Differential(*ampl,Variations_Mode::nominal_only,4);
    }
  } while (!InitSubtermInfo());
  ampl->Delete();

  p_bviproc->SetEventReader(p_read);
  bool res(p_bviproc->CalculateTotalXSec(resultpath,create));
  p_bviproc->SetEventReader(NULL);
  if (m_qcd) {
    if (!p_read) {
      p_viprocqcd->SetEventReader(p_read);
      res &= p_viprocqcd->CalculateTotalXSec(resultpath,create);
      p_viprocqcd->SetEventReader(NULL);
    }
    psh=p_rsprocqcd->Integrator()->PSHandler();
    if (!p_read && psh->AbsError()==0.0)
      psh->SetAbsError(psh->Error()*rpa->Picobarn()*
          dabs(p_bviproc->Integrator()->TotalResult()));
    p_rsprocqcd->SetEventReader(p_read);
    res &= p_rsprocqcd->CalculateTotalXSec(resultpath,create);
    p_rsprocqcd->SetEventReader(NULL);
  }
  if (m_ew) {
    p_viprocew->SetEventReader(p_read);
    res &= p_viprocew->CalculateTotalXSec(resultpath,create);
    p_viprocew->SetEventReader(NULL);
    psh=p_rsprocqed->Integrator()->PSHandler();
    if (!p_read && psh->AbsError()==0.0)
      psh->SetAbsError(psh->Error()*rpa->Picobarn()*
          dabs(p_bviproc->Integrator()->TotalResult()));
    p_rsprocqed->SetEventReader(p_read);
    res &= p_rsprocqed->CalculateTotalXSec(resultpath,create);
    p_rsprocqed->SetEventReader(NULL);
  }

  if (p_read) {
    p_read->WarmUp();
    p_int->SetMax(p_read->UnitWeight()/rpa->Picobarn());
  }
  for (size_t i(0);i<p_bviproc->Size();++i)
    (*p_bproc)[i]->Integrator()->SetMax
      ((*p_bviproc)[i]->Integrator()->Max());
  return res;
}

void MCatNLO_Process::FinalizeEventReader()
{
  // Sub-processes only ever hold borrowed copies of our p_read (see
  // OneEvent / CalculateTotalXSec); clear them before the base class
  // deletes the reader we own.
  if (p_bviproc) p_bviproc->FinalizeEventReader();
  if (p_viprocqcd) p_viprocqcd->FinalizeEventReader();
  if (p_viprocew) p_viprocew->FinalizeEventReader();
  if (p_rsprocqcd) p_rsprocqcd->FinalizeEventReader();
  if (p_rsprocqed) p_rsprocqed->FinalizeEventReader();
  if (p_bproc) p_bproc->FinalizeEventReader();
  if (p_rprocqcd) p_rprocqcd->FinalizeEventReader();
  if (p_rprocqed) p_rprocqed->FinalizeEventReader();
  if (p_ddprocqcd) p_ddprocqcd->FinalizeEventReader();
  if (p_ddprocqed) p_ddprocqed->FinalizeEventReader();
  Process_Base::FinalizeEventReader();
}

void MCatNLO_Process::SetLookUp(const bool lookup)
{
  m_lookup=lookup;
    if (m_qcd) {
    if (!p_read) p_viprocqcd->SetLookUp(lookup);
    p_ddprocqcd->SetLookUp(lookup);
    p_rsprocqcd->SetLookUp(lookup);
  }
  if (m_ew) {
    p_viprocew->SetLookUp(lookup);
    p_ddprocqed->SetLookUp(lookup);
    p_rsprocqed->SetLookUp(lookup);
  }
  p_bviproc->SetLookUp(lookup);
}

bool MCatNLO_Process::InitScale()
{
  bool res(true);
  if (m_qcd) {
    if (!p_read && !p_viprocqcd->InitScale()) res=false;
    if (!p_ddprocqcd->InitScale()) res=false;
    if (!p_rsprocqcd->InitScale()) res=false;
    if (!p_rprocqcd->InitScale()) res=false;
  }
  if (m_ew) {
    if (!p_viprocew->InitScale()) res=false;
    if (!p_ddprocqed->InitScale()) res=false;
    if (!p_rsprocqed->InitScale()) res=false;
    if (!p_rprocqed->InitScale()) res=false;
  }
  if (!p_bviproc->InitScale()) res=false;
  if (!p_bproc->InitScale()) res=false;
  return res;
}

void MCatNLO_Process::SetScale(const Scale_Setter_Arguments &scale)
{
  if (m_qcd) {
    if (!p_read) p_viprocqcd->SetScale(scale);
    p_ddprocqcd->SetScale(scale);
    if (m_rsscale!="") {
      Scale_Setter_Arguments rsscale(scale);
      rsscale.m_scale=m_rsscale;
      p_rsprocqcd->SetScale(rsscale);
      p_rprocqcd->SetScale(rsscale);
    }
    else {
      p_rsprocqcd->SetScale(scale);
      p_rprocqcd->SetScale(scale);
    }
  }
  if (m_ew) {
    p_viprocew->SetScale(scale);
    p_ddprocqed->SetScale(scale);
    if (m_rsscale!="") {
      Scale_Setter_Arguments rsscale(scale);
      rsscale.m_scale=m_rsscale;
      p_rsprocqed->SetScale(rsscale);
      p_rprocqed->SetScale(rsscale);
    }
    else {
      p_rsprocqed->SetScale(scale);
      p_rprocqed->SetScale(scale);
    }
  }
  p_bviproc->SetScale(scale);
  p_bproc->SetScale(scale);
}

void MCatNLO_Process::SetKFactor(const KFactor_Setter_Arguments &args)
{
  if (m_qcd) {
    if (!p_read) p_viprocqcd->SetKFactor(args);
    p_ddprocqcd->SetKFactor(args);
    p_rsprocqcd->SetKFactor(args);
    p_rprocqcd->SetKFactor(args);
  }
  if (m_ew) {
    p_viprocew->SetKFactor(args);
    p_ddprocqed->SetKFactor(args);
    p_rsprocqed->SetKFactor(args);
    p_rprocqed->SetKFactor(args);
  }
  p_bviproc->SetKFactor(args);
  p_bproc->SetKFactor(args);
}

void MCatNLO_Process::InitializeTheReweighting(ATOOLS::Variations_Mode mode)
{
  if (m_qcd) {
    if (!p_read) p_viprocqcd->InitializeTheReweighting(mode);
    p_ddprocqcd->InitializeTheReweighting(mode);
    p_rsprocqcd->InitializeTheReweighting(mode);
    p_rprocqcd->InitializeTheReweighting(mode);  
  }
  if (m_ew) {
    p_viprocew->InitializeTheReweighting(mode);
    p_ddprocqed->InitializeTheReweighting(mode);
    p_rsprocqed->InitializeTheReweighting(mode);
    p_rprocqed->InitializeTheReweighting(mode);
  }
  p_bviproc->InitializeTheReweighting(mode);
  p_bproc->InitializeTheReweighting(mode);
}

void MCatNLO_Process::SetFixedScale(const std::vector<double> &s)
{
  if (m_qcd) {
    if (!p_read) p_viprocqcd->SetFixedScale(s);
    p_ddprocqcd->SetFixedScale(s);
    p_rsprocqcd->SetFixedScale(s);
    p_rprocqcd->SetFixedScale(s);
  }
  if (m_ew) {
    p_viprocew->SetFixedScale(s);
    p_ddprocqed->SetFixedScale(s);
    p_rsprocqed->SetFixedScale(s);
    p_rprocqed->SetFixedScale(s);
  }
  p_bviproc->SetFixedScale(s);
  p_bproc->SetFixedScale(s);
}

bool MCatNLO_Process::IsGroup() const
{
  return true;
}

size_t MCatNLO_Process::Size() const
{
  if (p_read) return 2;
  if (m_qcd&m_ew) return 5;
  return 3;
}

Process_Base *MCatNLO_Process::operator[](const size_t &i)
{
  if (i==0) return p_bviproc;
  if (m_qcd) {
    if (i==1) return p_rsprocqcd;
    if (i==2 && !p_read) return p_viprocqcd;
  }
  if (m_ew) {
    if (i==1+2*m_qcd) return p_rsprocqed;
    if (i==2+2*m_qcd) return p_viprocew;
  }
  THROW(fatal_error,"Invalid MCatNLO_Process index");
}

void MCatNLO_Process::SetSelector(const Selector_Key &key)
{
  if (m_qcd) {
    if (!p_read) p_viprocqcd->SetSelector(key);
    p_ddprocqcd->SetSelector(key);
    p_rsprocqcd->SetSelector(key);
    p_rprocqcd->SetSelector(key);
  }
  if (m_ew) {
    p_viprocew->SetSelector(key);
    p_ddprocqed->SetSelector(key);
    p_rsprocqed->SetSelector(key);
    p_rprocqed->SetSelector(key);
  }
  p_bviproc->SetSelector(key);
  p_bproc->SetSelector(key);
}

void MCatNLO_Process::SetShower(PDF::Shower_Base *const ps)
{
  p_shower=ps;
  if (m_qcd) {
    if (!p_read) p_viprocqcd->SetShower(ps);
    p_ddprocqcd->SetShower(ps);
    p_rsprocqcd->SetShower(ps);
    p_rprocqcd->SetShower(ps);
  }
  if (m_ew) {
    p_viprocew->SetShower(ps);
    p_ddprocqed->SetShower(ps);
    p_rsprocqed->SetShower(ps);
    p_rprocqed->SetShower(ps);
  }
  p_bviproc->SetShower(ps);
  p_bproc->SetShower(ps);
}

void MCatNLO_Process::SetNLOMC(PDF::NLOMC_Base *const mc)
{
  p_nlomc=mc;
  if (m_qcd) {
    if (!p_read) p_viprocqcd->SetNLOMC(mc);
    p_ddprocqcd->SetNLOMC(mc);
    p_rsprocqcd->SetNLOMC(mc);
    p_rprocqcd->SetNLOMC(mc);
  }
  if (m_ew) {
    p_viprocew->SetNLOMC(mc);
    p_ddprocqed->SetNLOMC(mc);
    p_rsprocqed->SetNLOMC(mc);
    p_rprocqed->SetNLOMC(mc);
  }
  p_bviproc->SetNLOMC(mc);
  p_bproc->SetNLOMC(mc);
}
