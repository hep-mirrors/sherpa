#include "PHASIC++/Process/POWHEG_Process.H"

#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "PHASIC++/Selectors/Selector.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Process/ME_Generators.H"
#include "PHASIC++/Channels/POWHEG_Multi_Channel.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "PHASIC++/Process/Single_Process.H"
#include "PDF/Main/POWHEG_Base.H"
#include "PDF/Main/Shower_Base.H"
#include "PDF/Main/Jet_Criterion.H"
#include "PDF/Main/Cluster_Definitions_Base.H"
#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Message.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaS.H"

using namespace ATOOLS;
using namespace PHASIC;
using namespace PDF;

POWHEG_Process::POWHEG_Process
(ME_Generators& gens,const Process_Vector *procs):
  m_gens(gens), p_bproc(NULL), p_viproc(NULL),
  p_rproc(NULL), p_sproc(NULL), p_powheg(NULL),
  m_pmap(this), p_mc(NULL), m_smode(1), p_ampl(NULL)
{
  m_tinfo=0;
  for (size_t i(0);i<procs->size();++i) {
    POWHEG_Process *proc(dynamic_cast<POWHEG_Process*>((*procs)[i]));
    if (proc!=NULL) m_pprocs.push_back(proc);
  }
  p_as=dynamic_cast<MODEL::Running_AlphaS*>
    (gens.Model()->GetScalarFunction("alpha_S"));
  static bool ref(false);
  if (!ref) {
    ref=true;
    rpa->gen.AddCitation
      (1,"The automation of POWHEG is published in \\cite{Hoeche:2010pf}.");
  }
}

POWHEG_Process::~POWHEG_Process()
{
  if (p_ampl) p_ampl->Delete();
  for (std::map<nlo_type::code,Process_Map*>::const_iterator
	 pmit(m_pmap.begin());pmit!=m_pmap.end();++pmit) delete pmit->second;
  if (p_sproc) delete p_sproc;
  if (p_rproc) delete p_rproc;
  if (p_viproc) delete p_viproc;
  if (p_bproc) delete p_bproc;
}

void POWHEG_Process::Init(const Process_Info &pi,
			  BEAM::Beam_Spectra_Handler *const beam,
			  PDF::ISR_Handler *const isr)
{
  Process_Info cpi(pi);
  cpi.m_fi.SetNLOType(nlo_type::born|nlo_type::loop|
		      nlo_type::vsub|nlo_type::real|nlo_type::rsub);
  if (cpi.m_fi.m_nloewtype&nlo_type::real)
    cpi.m_fi.m_ps.push_back(Subprocess_Info(kf_photon,"",""));
  else
    cpi.m_fi.m_ps.push_back(Subprocess_Info(kf_jet,"",""));
  Process_Base::Init(cpi,beam,isr);
  m_pinfo.m_fi.SetNLOType(pi.m_fi.NLOType());
  Process_Info npi(m_pinfo);
  npi.m_fi.m_ps.pop_back();
  m_name=GenerateName(m_pinfo.m_ii,npi.m_fi);
  if (pi.Has(nlo_type::real)!=pi.Has(nlo_type::rsub))
    THROW(fatal_error, "R/S can't be initialised separately.");
  Process_Info spi(pi);
  spi.m_fi.SetNLOType(cpi.m_fi.NLOType());
  p_bproc=InitProcess(spi,nlo_type::lo,false);
  p_viproc=InitProcess(spi,nlo_type::loop|nlo_type::vsub|
		       (pi.m_fi.NLOType()&nlo_type::polecheck),false);
  p_rproc=InitProcess(spi,nlo_type::lo,true);
  p_sproc=InitProcess(spi,nlo_type::real|nlo_type::rsub,true);
  Data_Reader read(" ",";","!","=");
  read.SetInputPath(rpa->GetPath());
  read.SetInputFile(rpa->gen.Variable("INTEGRATION_DATA_FILE"));
  if (!read.ReadFromFile(m_rmode,"PP_RMODE")) m_rmode=1;
  else msg_Info()<<METHOD<<"(): Set RB integration mode "<<m_rmode<<".\n";
  for (size_t i(0);i<m_pprocs.size();++i) {
    if (m_pprocs[i]->Name()==Name()) THROW(fatal_error,"Doubled process");
    for (std::map<nlo_type::code,Process_Map*>::const_iterator
	   pmit(m_pprocs[i]->m_pmap.begin());
	 pmit!=m_pprocs[i]->m_pmap.end();++pmit) {
      m_pmap[pmit->first]->insert(pmit->second->begin(),pmit->second->end());
    }
  }
  if (p_sproc->Size()!=p_rproc->Size())
    THROW(fatal_error,"R and S have different size");
  for (size_t i(0);i<p_sproc->Size();++i)
    if ((*p_sproc)[i]->Flavours()!=(*p_rproc)[i]->Flavours())
      THROW(fatal_error,"Ordering differs in S and R");
}

Process_Base* POWHEG_Process::InitProcess
(const Process_Info &pi,nlo_type::code nlotype,const bool real)
{
  Process_Info cpi(pi);
  cpi.m_fi.SetNLOType(nlotype);
  cpi.m_megenerator="Amegic";
  if (real) {
    if (cpi.m_fi.m_nloqcdtype==nlotype)
      cpi.m_fi.m_ps.push_back(Subprocess_Info(kf_jet,"",""));
    else if (cpi.m_fi.m_nloewtype==nlotype)
      cpi.m_fi.m_ps.push_back(Subprocess_Info(kf_photon,"",""));
    else THROW(fatal_error, "Internal error.");
  }
  Process_Base *proc(m_gens.InitializeProcess(cpi,false));
  if (proc==NULL) {
    msg_Error()<<cpi<<std::endl;
    THROW(not_implemented, "Process not found.");
  }
  proc->SetParent(this);
  for (size_t i=0;i<proc->Size();++i) {
    std::string fname((*proc)[i]->Name());
    size_t pos=fname.find("EW");
    if (pos!=std::string::npos) fname=fname.substr(0,pos-2);
    pos=fname.find("QCD");
    if (pos!=std::string::npos) fname=fname.substr(0,pos-2);
    if (nlotype&nlo_type::vsub) nlotype=nlo_type::vsub;
    if (nlotype&nlo_type::rsub) nlotype=nlo_type::rsub;
    if (m_pmap.find(nlotype)==m_pmap.end()) m_pmap[nlotype] = new Process_Map();
    (*m_pmap[nlotype])[fname]=(*proc)[i];
  }
  return proc;
}

bool POWHEG_Process::InitSubtermInfo()
{
  std::map<std::string,size_t> idmap;
  for (size_t i(0);i<p_sproc->Size();++i) {
    NLO_subevtlist *subs((*p_sproc)[i]->GetSubevtList());
    for (size_t j(0);j<subs->size()-1;++j) {
      NLO_subevt *sub((*subs)[j]);
      std::string tag(sub->IDString());
      std::map<std::string,size_t>::iterator dit(idmap.find(tag));
      if (dit==idmap.end()) {
	idmap[tag]=m_rho.size();
	m_rho.push_back(0.0);
	dit=idmap.find(tag);
      }
      m_dmap[sub->p_id]=dit->second;
      if (sub->m_nco==2) continue;
      for (size_t ij(0);ij<sub->m_n;++ij)
	for (size_t k(0);k<sub->m_n;++k)
	  if (k!=ij && sub->p_fl[k]==sub->p_fl[sub->m_kt] && 
	      sub->p_fl[ij]==sub->p_fl[sub->m_ijt])
	    m_iinfo[sub->m_pname].insert(IDip_ID(ij,k));
      m_dinfo[subs->back()->m_pname].insert(*sub);
      Process_Base *bproc(FindProcess(sub));
      bproc->Integrator()->RBMap()[sub] = new RB_Data();
    }
  }
  return true;
}

Process_Base *POWHEG_Process::FindProcess
(const Cluster_Amplitude *ampl,const nlo_type::code type,
 const bool error) const
{
  std::string name(Process_Base::GenerateName(ampl));
  Process_Map::const_iterator pit(m_pmap.find(type)->second->find(name));
  if (pit!=m_pmap.find(type)->second->end()) return pit->second;
  if (error)
    THROW(fatal_error,"Process '"+name+"'("+ToString(type)+") not found");
  return NULL;
}

Process_Base *POWHEG_Process::FindProcess
(const NLO_subevt *sub,const nlo_type::code type) const
{
  Process_Map::const_iterator pit
    (m_pmap.find(type)->second->find(sub->m_pname));
  if (pit!=m_pmap.find(type)->second->end()) return pit->second;
  THROW(fatal_error,"Process '"+sub->m_pname+"'("+ToString(type)+") not found");
  return NULL;
}

Cluster_Amplitude *POWHEG_Process::CreateAmplitude(const NLO_subevt *sub) const
{
  Cluster_Amplitude *ampl = Cluster_Amplitude::New();
  ampl->SetNIn(m_nin);
  ampl->SetMS(p_sproc->Generator());
  ampl->SetMuF2(sub->m_muf2);
  ampl->SetMuR2(sub->m_mur2);
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

void POWHEG_Process::SetRBMap(Cluster_Amplitude *ampl)
{
  if (ampl->Legs().size()==m_nin+m_nout) {
    Process_Base::SetRBMap(ampl);
    return;
  }
  Process_Base::SortFlavours(ampl);
  Process_Base *proc(FindProcess(ampl,nlo_type::lo,false));
  if (proc==NULL) THROW(fatal_error,"No process for '"
			+Process_Base::GenerateName(ampl)+"'.");
  do {
    ampl->SetRBMap(&proc->Integrator()->RBMap());
    ampl->SetProcs(&m_pmap);
    ampl->SetIInfo(&m_iinfo);
    ampl->SetDInfo(&m_dinfo);
    ampl=ampl->Next();
  } while (ampl);
}

double POWHEG_Process::GetRho(const int mode)
{
  DEBUG_FUNC(mode);
  double asum(0.0);
  size_t aidx(p_mc->ActiveIdx());
  for (size_t i(0);i<m_rho.size();++i) {
    m_rho[i]=0.0;
    asum+=p_mc->SelectionWeight(i);
  }
  std::map<Process_Base*,Rho_Map> rhomap;
  for (size_t i(0);i<p_sproc->Size();++i) {
    NLO_subevtlist *subs((*p_sproc)[i]->GetSubevtList());
    (*m_pmap[nlo_type::lo])[subs->back()->m_pname]->
      SetLast(subs->back()->m_last[mode],mode);
    msg_Debugging()<<i<<": "<<subs->back()->m_pname<<"\n";
    msg_Indent();
    Rho_Map rhops;
    if (m_smode!=1) {
      if ((*p_sproc)[i]->IsMapped()) {
	rhops=rhomap[(*p_sproc)[i]->MapProc()];
	msg_Debugging()<<"map rho from '"<<
	  (*p_sproc)[i]->MapProc()->Name()<<"'\n";
      }
      else {
	Cluster_Amplitude *rampl(CreateAmplitude(subs->back()));
	rampl->SetProcs(&m_pmap);
	rampl->SetDInfo(&m_dinfo);
	bool lookup(m_lookup);
	p_bproc->SetLookUp(false);
	rhops=p_powheg->GetRho(rampl);
	p_bproc->SetLookUp(lookup);
	rampl->Delete();
      }
      rhomap[(*p_sproc)[i]]=rhops;
    }
    double rsum(0.0);
    std::vector<double> rhos(m_rho.size(),0.0);
    for (size_t j(0);j<subs->size()-1;++j) {
      NLO_subevt *sub((*subs)[j]);
      size_t id(m_dmap.find(sub->p_id)->second);
      double wa(p_mc->SelectionWeight(id)/asum);
      if (wa==0.0) continue;
      m_rho[id]+=sub->m_last[mode]/wa;
      msg_Debugging()<<"S '"<<sub->m_pname<<"'[("<<sub->m_i
		     <<","<<sub->m_j<<")("<<sub->m_k<<")|"<<id
		     <<"] -> "<<-sub->m_me;
      if (m_smode==1) {
	rhos[id]+=sub->m_me;
	rsum+=sub->m_me;
      }
      else {
	Rho_Map::const_iterator rit(rhops.find(*sub));
	if (rit==rhops.end()) THROW(fatal_error,"Internal error");
	if (rit->second.m_t>rhops.m_t0) {
	  msg_Debugging()<<" / "<<(rit->second.m_trig?rit->second.m_d:0.0)
			 <<" ("<<rit->second.m_d<<")";
	  rhos[id]+=rit->second.m_d;
	  rsum+=rit->second.m_d;
	}
	else {
	  rhos[id]+=sub->m_me;
	  rsum+=sub->m_me;
	}
      }
      msg_Debugging()<<"\n";
    }
    if (rsum==0.0) {
      m_zh[mode].push_back(ZH_Key(subs->back(),(*p_rproc)[i],
				  subs->back()->m_last[mode]));
      m_zhsum[mode]+=subs->back()->m_last[mode];
      m_rho[aidx]+=subs->back()->m_last[mode];
      continue;
    }
    else {
      for (size_t j(0);j<rhos.size();++j) {
	rhos[j]/=rsum;
	double wa(p_mc->SelectionWeight(j)/asum);
	if (wa==0.0) continue;
	double res(subs->back()->m_last[mode]*rhos[j]/wa);
	m_rho[j]+=res;
	if (res) msg_Debugging()
	  <<"R w/ B '"<<subs->back()->m_pname<<"'["<<j<<"] -> "
	  <<res<<" ("<<subs->back()->m_last[mode]<<")\n";
      }
    }
    if (m_smode==1) continue;
    double qij2(-1.0);
    for (size_t j(0);j<subs->size()-1;++j) {
      NLO_subevt *sub((*subs)[j]);
      size_t id(m_dmap.find(sub->p_id)->second);
      if (id==aidx && subs->back()->m_me && sub->m_nco!=2) {
	Process_Base *bproc((*m_pmap[nlo_type::lo])[sub->m_pname]);
	RB_Data *rbd(bproc->Integrator()->RBMap()[sub]);
	if (rbd->m_ktres==0.0) continue;
	if (qij2<0.0) qij2=PDF::Qij2Min(p_mc->RealMoms(),subs);
	if (qij2<p_powheg->KT2Min()) continue;
        ZH_Pair zh(p_powheg->ZHSplit
                   (bproc->Get<Single_Process>()->LastXS()*
		    bproc->SymFac(),qij2,rbd));
        double res(zh.second/(zh.first+zh.second)
		   *dabs(rhos[id])*subs->back()->m_last[mode]);
        msg_Debugging()<<bproc->Name()<<": B = "
                       <<bproc->Get<Single_Process>()->LastXS()*bproc->SymFac()
                       <<", R = "<<subs->back()->m_last[mode]
		       <<", Z = "<<zh.first<<", H = "<<zh.second
		       <<", \\rho["<<id<<"] = "<<rhos[id]<<" -> "<<res<<"\n";
	m_zh[mode].push_back(ZH_Key(sub,(*p_rproc)[i],res));
	m_zhsum[mode]+=res;
      }
    }
  }
  if (m_smode==1) return 1.0;
  return m_rho[aidx]/p_sproc->Last(mode);
}

double POWHEG_Process::Differential(const Vec4D_Vector &p)
{
  m_last[0]=0.0;
  if (!m_pinfo.Has(nlo_type::real)) p_sproc->MultiplyLast(0.0,0);
  else {
    p_sproc->Differential(p);
    if (m_smode==0) p_sproc->MultiplyLast(GetRho(0),0);
    if ((m_smode&1) && p_powheg && 
	(m_rmode==0 || (*p_bproc)[0]->Integrator()->WeightHisto())) {
      for (size_t i(0);i<p_sproc->Size();++i) {
	NLO_subevtlist *subs((*p_sproc)[i]->GetSubevtList());
	if (subs->empty() || !subs->back()->IsReal())
	  THROW(fatal_error,"Internal error");
        if (subs->back()->m_last[0]==0.0) continue;
	Cluster_Amplitude *rampl(CreateAmplitude(subs->back()));
	rampl->SetProcs(&m_pmap);
	rampl->SetIInfo(&m_iinfo);
	rampl->SetDInfo(&m_dinfo);
	bool lookup(m_lookup);
	SetLookUp(false);
	p_powheg->AddRBPoint(rampl);
	SetLookUp(lookup);
	rampl->Delete();
      }
      p_int->PSHandler()->CalculatePS();
    }
  }
  p_bproc->Differential(*p_mc->BAmpl());
  p_bproc->MultiplyLast(p_mc->BRWeight(),0);
  m_lastb=p_bproc->Last(0);
  if (!m_pinfo.Has(nlo_type::born)) p_bproc->MultiplyLast(0.0,0);
  if (!m_pinfo.Has(nlo_type::vsub)) p_viproc->MultiplyLast(0.0,0);
  else {
    if (p_viproc->Differential(*p_mc->BAmpl())!=0.0)
      p_viproc->MultiplyLast(p_mc->BRWeight(),0);
  }
  m_last[0]=p_sproc->Last(0)+p_bproc->Last(0)+p_viproc->Last(0);
  m_lastvi=p_viproc->Last(0);
  m_lastrs=p_sproc->Last(0);
  return m_last[0];
}

double POWHEG_Process::Differential2()
{
  m_last[1]=0.0;
  if (m_nin!=2 || p_int->ISR()==NULL ||
      p_int->ISR()->On()!=3) return 0.0;
  p_bproc->Differential2();
  p_bproc->MultiplyLast(p_mc->BRWeight(),1);
  if (m_pinfo.Has(nlo_type::real)) {
    p_sproc->Differential2();
    if (m_smode==0) p_sproc->MultiplyLast(GetRho(1),1);
  }
  m_lastb+=p_bproc->Last(1);
  if (!m_pinfo.Has(nlo_type::born)) p_bproc->MultiplyLast(0.0,1);
  if (!m_pinfo.Has(nlo_type::vsub)) p_viproc->MultiplyLast(0.0,1);
  else {
    if (p_viproc->Differential2()!=0.0)
      p_viproc->MultiplyLast(p_mc->BRWeight(),1);
  }
  m_last[1]=p_sproc->Last(1)+p_bproc->Last(1)+p_viproc->Last(1);
  m_lastvi+=p_viproc->Last(1);
  m_lastrs+=p_sproc->Last(1);
  return m_last[1];
}

double POWHEG_Process::LocalKFactor(const Vec4D_Vector &ppb)
{
  DEBUG_FUNC(Name());
  p_int->RestoreInOrder();
  Vec4D_Vector pb(ppb);
  if (pb[0][3]<0.0)
    for (size_t i=0;i<pb.size();++i) pb[i]=Vec4D(pb[i][0],-pb[i]);
  p_mc->GenerateEmissionPoint(pb);
  p_mc->GenerateWeight(NULL,NULL);
  m_smode=0;
  for (size_t i(0);i<2;++i) m_zh[i].clear();
  m_zhsum[1]=m_zhsum[0]=0.0;
  Differential(p_mc->RealMoms());
  Differential2();
  m_smode=1;
  if (IsBad(Last()) || IsBad(p_bproc->Last())) {
    PRINT_INFO("nan encountered, returning 1.0");
    return 1.0;
  }
  if (p_bproc->Last()==0.0) {
    msg_Info()<<METHOD<<"(): Born zero. Returning 1.\n";
    // probably because ismc->Weight()==0.0
    // see this mainly near the mass cut threshold so far
    DEBUG_VAR(Name());
    DEBUG_VAR(ppb);
    DEBUG_VAR(Last());
    DEBUG_VAR((ppb[2]+ppb[3]).Mass());
    return 1.0;
  }
  double xsc(1.0/(1.0+dabs(Last()-(m_zhsum[0]+m_zhsum[1]))/
		  dabs(m_zhsum[0]+m_zhsum[1])));
  if (fabs((1.0-xsc)*Last()/p_bproc->Last())>1000.0) {
    msg_Info()<<METHOD<<"(): Large local kfactor ( K = "
	      <<(1.0-xsc)*Last()/p_bproc->Last()<<" ). Returning 1.\n";
    DEBUG_VAR(Name());
    DEBUG_VAR(ppb);
    DEBUG_VAR(p_bproc->Last());
    DEBUG_VAR(Last()<<" "<<xsc);
    DEBUG_VAR((ppb[2]+ppb[3]).Mass());
    return 1.0;
  }
  DEBUG_INFO("Found kfactor "<<(1.0-xsc)*Last()/p_bproc->Last());
  return (1.0-xsc)*Last()/p_bproc->Last();
}

Cluster_Amplitude *POWHEG_Process::GetAmplitude()
{
  Cluster_Amplitude *ampl(p_ampl);
  p_ampl=NULL;
  return ampl;
}

double POWHEG_Process::SelectProcess()
{
  DEBUG_FUNC("");
  if (p_ampl) p_ampl->Delete();
  p_ampl=NULL;
  int mode(p_int->InSwaped());
  double xsc(1.0/(1.0+dabs(m_last[mode]-m_zhsum[mode])/dabs(m_zhsum[mode])));
  if (xsc>ran->Get()) return SelectZHProcess();
  return SelectBProcess();
}

double POWHEG_Process::SelectRProcess()
{
  int mode(p_int->InSwaped());
  msg_Debugging()<<"R-like event ( RS = "<<p_sproc->Last(mode)<<" )\n";
  p_selected=p_rproc;
  double rsum(0.0), psum(0.0);
  for (size_t i(0);i<p_rproc->Size();++i) rsum+=(*p_rproc)[i]->Last(mode);
  rsum*=ran->Get();
  for (size_t i(0);i<p_rproc->Size();++i)
    if ((psum+=(*p_rproc)[i]->Last(mode))>=rsum) {
      p_rproc->SetSelected((*p_rproc)[i]);
      return (*p_rproc)[i]->Integrator()->SelectionWeight(0)/
	p_int->SelectionWeight(0);
    }
  THROW(fatal_error,"Internal error");
}

double POWHEG_Process::SelectZHProcess()
{
  int mode(p_int->InSwaped());
  msg_Debugging()<<"RB-like event ( ZH = "<<m_zhsum[mode]<<" )\n";
  if (m_zhsum[mode]==0.0) return 0.0*SelectRProcess();
  p_selected=p_rproc;
  double zhsum(0.0), psum(0.0);
  for (size_t i(0);i<m_zh[mode].size();++i)
    zhsum+=dabs(m_zh[mode][i].m_res);
  zhsum*=ran->Get();
  for (size_t i(0);i<m_zh[mode].size();++i) {
    if ((psum+=dabs(m_zh[mode][i].m_res))>=zhsum) {
      Process_Base *rproc(m_zh[mode][i].p_r);
      p_rproc->SetSelected(rproc);
      rproc->Integrator()->SetMomenta(p_mc->RealMoms());
      if (p_int->InSwaped()) rproc->Integrator()->SwapInOrder();
      Selector_Base *jf=(*p_bproc)[0]->
	Selector()->GetSelector("Jetfinder");
      if (jf) {
	Cluster_Amplitude *rampl
	  (CreateAmplitude(m_zh[mode][i].p_sub->p_real));
	rampl->SetJF(jf);
	bool res(p_shower->JetVeto(rampl));
	rampl->Delete();
	if (res) return 0.0;
      }
      return rproc->Integrator()->SelectionWeight(0)/
	p_int->SelectionWeight(0);
    }
  }
  THROW(fatal_error,"Internal error");
}

double POWHEG_Process::SelectBProcess()
{
  int mode(p_int->InSwaped());
  msg_Debugging()<<"B-like event ( B = "<<p_bproc->Last(mode)<<" )\n";
  if (p_bproc->Last(mode)==0.0) return SelectRProcess();
  double rsum(p_bproc->Last(mode)*ran->Get()), psum(0.0);
  Process_Base *bproc(NULL);
  for (size_t i(0);i<p_bproc->Size();++i) {
    bproc=(*p_bproc)[i];
    if ((psum+=bproc->Last(mode))>=rsum) break;
  }
  p_ampl = dynamic_cast<Single_Process*>(bproc)->Cluster(256);
  SortFlavours(p_ampl);
  if (!p_ampl) {
    msg_Error()<<METHOD<<"(): no cluster amplitude available. new event.\n";
    return 0.;
  }
  p_ampl->SetProcs(&m_pmap);
  p_ampl->SetIInfo(&m_iinfo);
  p_ampl->SetDInfo(&m_dinfo);
  p_ampl->SetRBMap(&bproc->Integrator()->RBMap());
  p_ampl->Decays()=m_decins;
  bool lookup(m_lookup);
  SetLookUp(false);
  p_powheg->SetShower(p_shower);
  int stat(p_powheg->GeneratePoint(p_ampl));
  SetLookUp(lookup);
  Cluster_Amplitude *next(p_ampl), *ampl(p_ampl->Prev());
  if (ampl) {
    p_ampl=NULL;
    Process_Base *rproc(FindProcess(ampl));
    if (rproc==NULL) THROW(fatal_error,"Invalid splitting");
    p_selected=p_rproc;
    p_rproc->SetSelected(rproc);
    rproc->Integrator()->SetMomenta(*ampl);
    Cluster_Leg *lij(NULL);
    for (size_t i(0);i<next->Legs().size();++i)
      if (next->Leg(i)->K()) {
	lij=next->Leg(i);
	break;
      }
    std::vector<int> ids(ID(lij->Id()));
    size_t iid(1<<ids.front()), jid(1<<ids.back()), kid(lij->K());
    ids.push_back(ID(lij->K()).front());
    if (ids.size()!=3) THROW(fatal_error,"Internal error");
    Cluster_Amplitude *campl(ampl->Copy());
    campl->IdSort();
    Cluster_Definitions_Base *clus(p_shower->GetClusterDefinitions());
    CParam kt2(clus->KPerp2(*campl,ids[0],ids[1],ids[2],
			    lij->Flav(),p_bproc->Generator(),1));
    Vec4D_Vector pp(clus->Combine(*campl,ids[0],ids[1],ids[2],lij->Flav(),
				  p_bproc->Generator(),kt2.m_kin));
    campl->Delete();
    ampl->SetKT2(kt2.m_kt2);
    ampl->SetMu2(kt2.m_kt2);
    for (size_t i(0);i<pp.size();++i)
      for (size_t j(0);j<next->Legs().size();++j)
	if (next->Leg(j)->Id()&(1<<i)) next->Leg(j)->SetMom(pp[i]);
    if (p_int->InSwaped())
      std::swap<Cluster_Leg*>(next->Legs()[0],next->Legs()[1]);
    bproc->Integrator()->SetMomenta(*next);
    next = dynamic_cast<Single_Process*>(bproc)->Cluster(256);
    (p_ampl=ampl)->SetNext(next);
    ampl->SetMuF2(next->MuF2());
    ampl->SetMuR2(next->MuR2());
    next->SetKin(kt2.m_kin);
    if (!(ampl->Leg(0)->Id()&next->Leg(0)->Id()))
      std::swap<Cluster_Leg*>(ampl->Legs()[0],ampl->Legs()[1]);
    while (ampl->Next()) {
      ampl=ampl->Next();
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
	    ampl->SetIdNew(ampl->Leg(i)->Id());
	  }
	}
      }
    }
    msg_Debugging()<<"R selected via Sudakov "<<*p_ampl<<"\n";
    return rproc->Integrator()->SelectionWeight(0)/
      p_int->SelectionWeight(0);
  }
  p_selected=p_bproc;
  p_bproc->SetSelected(bproc);
  p_bproc->Integrator()->SetMomenta(*p_ampl);
  msg_Debugging()<<"B selected "<<*p_ampl<<"\n";
  return stat?bproc->Integrator()->SelectionWeight(0)/
    p_int->SelectionWeight(0):0.0;
}

Weight_Info *POWHEG_Process::OneEvent(const int wmode,const int mode)
{
  m_smode=0;
  p_selected=this;
  p_rproc->Integrator()->RestoreInOrder();
  for (size_t i(0);i<2;++i) m_zh[i].clear();
  m_zhsum[1]=m_zhsum[0]=0.0;
  Weight_Info *winfo(p_int->PSHandler()->OneEvent(this,mode));
  if (winfo!=NULL) {
    double cf(SelectProcess());
    if (cf) winfo->m_weight*=cf;
    else {
      delete winfo;
      winfo=NULL;
    }
  }
  m_smode=1;
  return winfo;
}

bool POWHEG_Process::CalculateTotalXSec(const std::string &resultpath,
					const bool create) 
{ 
  p_int->Reset();
  SP(Phase_Space_Handler) psh(p_int->PSHandler());
  p_sproc->Integrator()->SetPSHandler(psh);
  p_rproc->Integrator()->SetPSHandler(psh);
  if (p_int->ISR()) {
    if (m_nin==2) {
      if (m_flavs[0].Mass()!=p_int->ISR()->Flav(0).Mass() ||
          m_flavs[1].Mass()!=p_int->ISR()->Flav(1).Mass()) {
        p_int->ISR()->SetPartonMasses(m_flavs);
      }
    }
  }
  psh->InitCuts();
  if (p_int->ISR())
    p_int->ISR()->SetSprimeMin(psh->Cuts()->Smin());
  psh->CreateIntegrators();
  p_int->SetResultPath(resultpath);
  p_int->ReadResults();
  exh->AddTerminatorObject(p_int);
  psh->InitIncoming();
  double var(p_int->TotalVar());
  msg_Info()<<METHOD<<"(): Calculate xs for '"<<m_name<<"' ("
            <<(p_sproc->Generator()?p_sproc->Generator()->Name():"unknown");
  if(m_pinfo.m_fi.NLOType()&nlo_type::loop)
    msg_Info()<<","<<m_pinfo.m_loopgenerator;
  msg_Info()<<")"<<std::endl;
  double totalxs(psh->Integrate()/rpa->Picobarn());
  if (!IsEqual(totalxs,p_int->TotalResult())) {
    msg_Error()<<"Result of PS-Integrator and summation do not coincide!\n"
	       <<"  '"<<m_name<<"': "<<totalxs
	       <<" vs. "<<p_int->TotalResult()<<std::endl;
  }
  if (p_int->TotalXS()>0.0) {
    p_int->SetTotal();
    if (var==p_int->TotalVar()) {
      exh->RemoveTerminatorObject(p_int);
      return 1;
    }
    p_int->StoreResults();
    p_int->StoreBackupResults();
    exh->RemoveTerminatorObject(p_int);
    return 1;
  }
  exh->RemoveTerminatorObject(p_int);
  return 0;
}

void POWHEG_Process::SetLookUp(const bool lookup)
{
  m_lookup=lookup;
  p_bproc->SetLookUp(lookup);
  p_viproc->SetLookUp(lookup);
  p_rproc->SetLookUp(lookup);
  p_sproc->SetLookUp(lookup);
}

void POWHEG_Process::SetScale(const Scale_Setter_Arguments &scale)
{
  p_bproc->SetScale(scale);
  p_viproc->SetScale(scale);
  p_rproc->SetScale(scale);
  p_sproc->SetScale(scale);
}

void POWHEG_Process::SetKFactor(const KFactor_Setter_Arguments &args)
{
  p_bproc->SetKFactor(args);
  p_viproc->SetKFactor(args);
  p_sproc->SetKFactor(args);
  p_rproc->SetKFactor(args);
}

void POWHEG_Process::SetFixedScale(const std::vector<double> &s)
{
  p_bproc->SetFixedScale(s);
  p_viproc->SetFixedScale(s);
  p_rproc->SetFixedScale(s);
  p_sproc->SetFixedScale(s);
}

bool POWHEG_Process::IsGroup() const
{
  return true;
}

size_t POWHEG_Process::Size() const
{
  return 3;
}

Process_Base *POWHEG_Process::operator[](const size_t &i)
{
  if (i==0) return p_sproc;
  if (i==1) return p_viproc;
  return p_bproc;
}

void POWHEG_Process::SetClusterDefinitions
(PDF::Cluster_Definitions_Base *const cluster)
{
  p_bproc->Generator()->SetClusterDefinitions(cluster);
  p_rproc->Generator()->SetClusterDefinitions(cluster);
}

void POWHEG_Process::SetSelector(const Selector_Key &key)
{
  p_bproc->SetSelector(key);
  p_viproc->SetSelector(key);
  p_rproc->SetSelector(key);
  p_sproc->SetSelector(key);
}

void POWHEG_Process::SetShower(PDF::Shower_Base *const ps)
{
  p_shower=ps;
  p_bproc->SetShower(ps);
  p_viproc->SetShower(ps);
  p_rproc->SetShower(ps);
  p_sproc->SetShower(ps);
}

bool POWHEG_Process::Trigger(const ATOOLS::Vec4D_Vector &p)
{
  return p_sproc->Trigger(p);
}

bool POWHEG_Process::FillIntegrator(Phase_Space_Handler *const psh)
{
  return p_bproc->FillIntegrator(psh);
}

bool POWHEG_Process::InitIntegrator(Phase_Space_Handler *const psh)
{
  p_mc = new POWHEG_Multi_Channel(this,p_sproc,psh);
  psh->SetFSRIntegrator(p_mc);
  return true;
}

void POWHEG_Process::UpdateIntegrator(Phase_Space_Handler *const psh)
{
  p_bproc->UpdateIntegrator(psh);
}

void POWHEG_Process::InitCuts(Cut_Data *const cuts)
{
  p_bproc->InitCuts(cuts);
}

void POWHEG_Process::BuildCuts(Cut_Data *const cuts)
{
  p_bproc->BuildCuts(cuts);
}

