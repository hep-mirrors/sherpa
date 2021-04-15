#include "COMIX/Main/Single_Process.H"

#include "COMIX/Main/Process_Group.H"
#include "COMIX/Main/Single_Dipole_Term.H"
#include "PDF/Main/ISR_Handler.H"
#include "COMIX/Phasespace/PS_Generator.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "PHASIC++/Channels/Multi_Channel.H"
#include "PHASIC++/Main/Color_Integrator.H"
#include "PHASIC++/Main/Helicity_Integrator.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Scales/KFactor_Setter_Base.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"

using namespace COMIX;
using namespace PHASIC;
using namespace ATOOLS;

COMIX::Single_Process::Single_Process():
  COMIX::Process_Base(this),
  p_bg(NULL), p_hc(NULL), p_map(NULL),
  p_loop(NULL), p_kpterms(NULL),
  m_checkpoles(false), m_allowmap(true)
{
  Settings& s = Settings::GetMainSettings();
  m_itype = s["NLO_IMODE"].Get<cs_itype::type>();
}

COMIX::Single_Process::~Single_Process()
{
  if (p_kpterms) delete p_kpterms;
  if (p_loop) delete p_loop;
  if (p_map) {
    for (size_t i(0);i<m_subs.size();++i) {
      delete [] m_subs[i]->p_id;
      delete [] m_subs[i]->p_fl;
      if (i<m_subs.size()-1) delete static_cast
	  <Single_Dipole_Term*>(m_subs[i]->p_proc);
      delete m_subs[i];
    }
  }
  else if (p_bg) {
    NLO_subevtlist *subs(GetSubevtList());
    if (subs && subs->size())
      for (size_t i(0);i<subs->size()-1;++i)
	     if ((*subs)[i]->p_proc)
	       delete static_cast
		 <Single_Dipole_Term*>((*subs)[i]->p_proc);
  }
  if (p_bg!=NULL) delete p_bg;
}

bool COMIX::Single_Process::Initialize
(std::map<std::string,std::string> *const pmap,
 std::vector<Single_Process*> *const procs,
 const std::vector<int> &blocks,size_t &nproc)
{
  DEBUG_FUNC("");
  DEBUG_VAR(m_pinfo);
  m_p.resize(m_nin+m_nout);
  if (!COMIX::Process_Base::Initialize
      (pmap,procs,blocks,nproc)) return false;
  if (p_bg!=NULL) delete p_bg;
  p_bg=NULL;
  if (pmap->find(m_name)!=pmap->end()) {
    if ((*pmap)[m_name]!="x") THROW(fatal_error,"Internal error");
    return false;
  }
  std::string mapfile(rpa->gen.Variable("SHERPA_CPP_PATH")
		      +"/Process/Comix/"+m_name+".map");
  m_hasmapfile=FileExists(mapfile,1);
  My_In_File mapstream(mapfile);
  if (m_hasmapfile) {
    msg_Tracking()<<METHOD<<"(): Map file '"<<mapfile<<"' found.\n";
    if (mapstream.Open()) {
      std::string cname, mapname;
      *mapstream>>cname>>mapname;
      if (cname!=m_name)
	THROW(fatal_error,"Corrupted map file '"+mapfile+"'");
      (*pmap)[m_name]=mapname;
      if (mapname!=m_name) {
	msg_Debugging()<<METHOD<<"(): Map '"<<m_name
		       <<"' onto '"<<mapname<<"'\n";
	if (Parent()->Name()==Name() && mapname=="x") return false;
	if (blocks.size()) {
	  for (size_t i(0);i<procs->size();++i)
	    if ((*procs)[i]->Name()==mapname) {
	      msg_Debugging()<<"Keep "<<m_name<<"\n";
	      return true;
	    }
	  return false;
	}
	return true;
      }
      else {
	++nproc;
	if (blocks.size()) {
	  if (std::find(blocks.begin(),blocks.end(),nproc-1)
	      ==blocks.end()) return false;
	  msg_Debugging()<<"Keep "<<nproc<<"\n";
	}
      }
    }
  }
  msg_Debugging()<<"'"<<m_name<<"' not pre-mapped"<<std::endl;
  p_model->GetCouplings(m_cpls);
  p_bg = new Amplitude();
  p_bg->SetCTS(p_cts);
  double isf(m_pinfo.m_ii.ISSymmetryFactor());
  double fsf(m_pinfo.m_fi.FSSymmetryFactor());
  Subprocess_Info info(m_pinfo.m_ii);
  info.Add(m_pinfo.m_fi);
  p_bg->SetDecayInfos(m_decins);
  p_bg->SetNoDecays(m_pinfo.m_nodecs);
  int stype(0), smode(0);
  // stype -> 0 QCD, 1 QED
  // smode -> 0 LO, 1 RS, 2 I, 4 B, 8 Polecheck, 16 V
  if (m_pinfo.m_fi.m_nlocpl.size()) {
    if (m_pinfo.m_fi.m_nlocpl[0]==0. &&
	m_pinfo.m_fi.m_nlocpl[1]==0.) stype=0;
    else if (m_pinfo.m_fi.m_nlocpl[0]==1. &&
	     m_pinfo.m_fi.m_nlocpl[1]==0.) stype=0;
    else if (m_pinfo.m_fi.m_nlocpl[0]==0. &&
	     m_pinfo.m_fi.m_nlocpl[1]==1.) stype=1;
    else THROW(fatal_error,"Cannot do NLO QCD+EW");
    msg_Debugging()<<"Subtraction type: "<<(sbt::subtype)(stype+1)<<"\n";
  }
  m_checkpoles = ToType<size_t>(rpa->gen.Variable("CHECK_POLES"));
  if (m_pinfo.m_fi.NLOType()&nlo_type::rsub) smode=1;
  else {
    if (m_pinfo.m_fi.NLOType()&nlo_type::vsub) smode|=2;
    if (m_pinfo.m_fi.NLOType()&nlo_type::born) smode|=4;
    if (m_pinfo.m_fi.NLOType()&nlo_type::loop) smode|=16;
    if (m_checkpoles) smode|=8;
  }
  msg_Debugging()<<"Subtraction Mode: "<<smode<<std::endl;
  std::vector<int> mincpl(m_pinfo.m_mincpl.size());
  std::vector<int> maxcpl(m_pinfo.m_maxcpl.size());
  std::vector<int> minacpl(m_pinfo.m_minacpl.size());
  std::vector<int> maxacpl(m_pinfo.m_maxacpl.size());
  for (size_t i(0);i<mincpl.size();++i) mincpl[i]=m_pinfo.m_mincpl[i]*2;
  for (size_t i(0);i<maxcpl.size();++i) maxcpl[i]=m_pinfo.m_maxcpl[i]*2;
  for (size_t i(0);i<minacpl.size();++i) minacpl[i]=m_pinfo.m_minacpl[i];
  for (size_t i(0);i<maxacpl.size();++i) maxacpl[i]=m_pinfo.m_maxacpl[i];
  if (smode&18)
    for (int i(0);i<2;++i) {
      maxcpl[i]-=m_pinfo.m_fi.m_nlocpl[i]*2;
      mincpl[i]-=m_pinfo.m_fi.m_nlocpl[i]*2;
      maxacpl[i]-=m_pinfo.m_fi.m_nlocpl[i];
      minacpl[i]-=m_pinfo.m_fi.m_nlocpl[i];
    }
  if (p_bg->Initialize(m_nin,m_nout,m_flavs,isf,fsf,mapstream,
		       &*p_model,&m_cpls,stype,smode,m_itype,
		       maxcpl,mincpl,maxacpl,minacpl,
		       m_pinfo.m_ntchan,m_pinfo.m_mtchan,m_name)) {
    m_mincpl.resize(p_bg->MinCpl().size());
    m_maxcpl.resize(p_bg->MaxCpl().size());
    for (size_t i(0);i<m_maxcpl.size();++i) {
      m_mincpl[i]=p_bg->MinCpl()[i]/2.0;
      m_maxcpl[i]=p_bg->MaxCpl()[i]/2.0;
    }
    if (smode&1) {
      NLO_subevtlist *subs(GetSubevtList());
      for (size_t i(0);i<subs->size()-1;++i)
	(*subs)[i]->p_proc =
	  new Single_Dipole_Term(this,(*subs)[i],(*subs)[i]);
      subs->back()->p_proc=this;
      for (size_t i(0);i<subs->size();++i)
	(*subs)[i]->p_real=subs->back();
    }
    if (smode&2) {
      int massive(0);
      std::vector<size_t> pl;
      for (size_t i(0);i<m_nin+m_nout;++i)
	if (m_flavs[i].Strong()) {
	  if (m_flavs[i].Mass()) massive=1;
	  pl.push_back(i);
	}
      p_bg->DInfo()->SetMassive(massive);
      p_kpterms = new KP_Terms(this, ATOOLS::sbt::qcd, pl);
      p_kpterms->SetIType(m_itype);
      p_kpterms->SetCoupling(&m_cpls);
      m_mewgtinfo.m_type|=
	((m_itype&cs_itype::I)?mewgttype::VI:mewgttype::none);
      m_mewgtinfo.m_type|=
	(((m_itype&cs_itype::K)||(m_itype&cs_itype::P))?
	 mewgttype::KP:mewgttype::none);
      if (smode&4) m_mewgtinfo.m_type|=mewgttype::B;
    }
    if (smode&16) {
      Process_Info cinfo(m_pinfo);
      cinfo.m_fi.m_nlotype=nlo_type::loop;
      for (size_t i(0);i<m_maxcpl.size();++i) {
	cinfo.m_mincpl[i]=p_bg->MinCpl()[i]/2.0+cinfo.m_fi.m_nlocpl[i];
	cinfo.m_maxcpl[i]=p_bg->MaxCpl()[i]/2.0+cinfo.m_fi.m_nlocpl[i];
      }
      p_loop = PHASIC::Virtual_ME2_Base::GetME2(cinfo);
      if (p_loop==NULL) {
	msg_Error()<<METHOD<<"(): "<<cinfo<<"\n";
	THROW(not_implemented,"No virtual ME for this process");
      }
      p_loop->SetCouplings(m_cpls);
      p_loop->SetNorm(1.0/(isf*fsf));
      p_loop->SetSubType((sbt::subtype)(stype+1));
      m_mewgtinfo.m_type|=mewgttype::VI;
      int helpi;
      Settings& s = Settings::GetMainSettings();
      m_checkpoles = s["CHECK_POLES"].Get<size_t>();
      p_loop->SetPoleCheck(m_checkpoles);
    }
    p_bg->SetLoopME(p_loop);
    p_bg->FillCombinations(m_ccombs,m_cflavs);
    if (m_pinfo.m_fi.m_nlotype&(nlo_type::loop|nlo_type::vsub))
      for (size_t i(0);i<m_maxcpl.size();++i) {
	m_mincpl[i]=p_bg->MinCpl()[i]/2.0+m_pinfo.m_fi.m_nlocpl[i];
	m_maxcpl[i]=p_bg->MaxCpl()[i]/2.0+m_pinfo.m_fi.m_nlocpl[i];
      }
    (*pmap)[m_name]=m_name;
    if (m_pinfo.m_cls==cls::sum) ConstructColorMatrix();
    p_hc = new Hard_Matrix();
    return true;
  }
  mapfile=rpa->gen.Variable("SHERPA_CPP_PATH")
    +"/Process/Comix/"+Parent()->Name()+".map";
  std::string str, tmp;
  My_In_File in(mapfile);
  if (in.Open())
    for (getline(*in,tmp);in->good();
	 getline(*in,tmp)) str+=tmp+"\n";
  in.Close();
  My_Out_File out(mapfile);
  if (!out.Open()) THROW(fatal_error,"Cannot open '"+mapfile+"'");
  *out<<str;
  *out<<m_name<<" x\n";
  out.Close();
  (*pmap)[m_name]="x";
  return false;
}

void COMIX::Single_Process::MapSubEvts(const int mode)
{
  m_subs.resize(p_map->p_bg->SubEvts().size());
  const NLO_subevtlist &subs(p_bg->SubEvts());
  const NLO_subevtlist &rsubs(p_map->p_bg->SubEvts());
  for (size_t i(0);i<m_subs.size();++i) {
    m_subs[i] = new NLO_subevt(*rsubs[i]);
    Flavour *fls(new Flavour[m_subs[i]->m_n]);
    size_t *ids(new size_t[m_subs[i]->m_n]);
    m_subs[i]->p_fl = fls;
    m_subs[i]->p_id = ids;
    m_subs[i]->p_mom=rsubs[i]->p_mom;
    m_subs[i]->p_dec=rsubs[i]->p_dec;
    for (size_t j(0);j<m_subs[i]->m_n;++j) {
      fls[j]=ReMap(rsubs[i]->p_fl[j],0);
      ids[j]=rsubs[i]->p_id[j];
    }
    if (i+1<m_subs.size()) {
      if (mode&1)
	delete static_cast<Single_Dipole_Term*>(subs[i]->p_proc);
      m_subs[i]->p_proc =
	new Single_Dipole_Term(this,rsubs[i],m_subs[i]);
    }
    else {
      m_subs[i]->p_proc=this;
    }
    m_subs[i]->m_pname=static_cast<PHASIC::Process_Base*>
      (m_subs[i]->p_proc)->Name();
  }
  for (size_t i(0);i<m_subs.size();++i)
    m_subs[i]->p_real=m_subs.back();
}

bool COMIX::Single_Process::MapProcess()
{
  std::string mapname((*p_pmap)[m_name]);
  if (mapname!=m_name) {
    for (size_t i(0);i<p_umprocs->size();++i)
      if ((*p_umprocs)[i]->Name()==mapname) {
	p_mapproc=p_map=(*p_umprocs)[i];
	m_maxcpl=p_map->m_maxcpl;
	m_mincpl=p_map->m_mincpl;
	m_mewgtinfo.m_type=p_map->m_mewgtinfo.m_type;
	if (p_map->p_kpterms) {
	  std::vector<size_t> pl;
	  for (size_t j(0);j<m_nin+m_nout;++j)
	    if (m_flavs[j].Strong()) pl.push_back(j);
	  p_kpterms = new KP_Terms(p_map, ATOOLS::sbt::qcd, pl);
	  p_kpterms->SetIType(p_map->m_itype);
	  p_kpterms->SetCoupling(&p_map->m_cpls);
	}
	msg_Tracking()<<"Mapped '"<<m_name<<"' -> '"<<mapname<<"'.\n";
	std::string mapfile(rpa->gen.Variable("SHERPA_CPP_PATH")
			    +"/Process/Comix/"+m_name+".map");
	My_In_File map(mapfile);
	if (!map.Open()) {
	  THROW(fatal_error,"Map file '"+mapfile+"' not found");
	}
	else {
	  size_t nfmap;
	  std::string cname, cmapname;
	  *map>>cname>>cmapname>>nfmap;
	  if (cname!=m_name || cmapname!=mapname || map->eof())
	    THROW(fatal_error,"Corrupted map file '"+mapfile+"'");
	  for (size_t j(0);j<nfmap;++j) {
	    long int src, dest;
	    *map>>src>>dest;
	    Flavour ft((kf_code)(std::abs(src)),src<0);
	    Flavour fb((kf_code)(std::abs(dest)),dest<0);
	    m_fmap[ft]=fb;
	    msg_Debugging()<<"  fmap '"<<ft<<"' onto '"<<fb<<"'\n";
	  }
	  *map>>cname;
	  if (cname!="eof")
	    THROW(fatal_error,"Corrupted map file '"+mapfile+"'");
	}
	MapSubEvts(0);
	p_map->p_bg->FillCombinations(m_ccombs,m_cflavs);
	for (CFlavVector_Map::iterator fit(m_cflavs.begin());
	     fit!=m_cflavs.end();++fit)
	  for (size_t i(0);i<fit->second.size();++i)
	    fit->second[i]=ReMap(fit->second[i],0);
	return true;
      }
    THROW(fatal_error,"Map process '"+mapname+"' not found");
  }
  if (!m_hasmapfile && m_allowmap &&
      m_pinfo.m_special.find("MapOff")==std::string::npos) {
  for (size_t i(0);i<p_umprocs->size();++i) {
    msg_Debugging()<<METHOD<<"(): Try mapping '"
		   <<Name()<<"' -> '"<<(*p_umprocs)[i]->Name()<<"'\n";
    if (p_bg->Map(*(*p_umprocs)[i]->p_bg,m_fmap)) {
      p_mapproc=p_map=(*p_umprocs)[i];
      if (p_kpterms) {
	delete p_kpterms;
	std::vector<size_t> pl;
	for (size_t j(0);j<m_nin+m_nout;++j)
	  if (m_flavs[j].Strong()) pl.push_back(j);
	p_kpterms = new KP_Terms(p_map, ATOOLS::sbt::qcd, pl);
	p_kpterms->SetIType(p_map->m_itype);
	p_kpterms->SetCoupling(&p_map->m_cpls);
      }
      mapname=p_map->Name();
      msg_Tracking()<<"Mapped '"<<m_name<<"' -> '"
		    <<mapname<<"'."<<std::endl;
      std::string mapfile(rpa->gen.Variable("SHERPA_CPP_PATH")
			  +"/Process/Comix/"+m_name+".map");
      if (!FileExists(mapfile,1)) {
	My_Out_File map(mapfile);
	if (map.Open()) {
	  *map<<m_name<<" "<<mapname<<"\n"<<m_fmap.size()<<"\n";
	  for (Flavour_Map::const_iterator 
		 fit(m_fmap.begin());fit!=m_fmap.end();++fit) {
	    msg_Debugging()<<"  fmap '"<<fit->first
			   <<"' onto '"<<fit->second<<"'\n";
	    long int src(fit->first), dest(fit->second);
	    *map<<src<<" "<<dest<<"\n";
	  }
	  *map<<"eof\n";
	}
      }
      MapSubEvts(1);
      delete p_bg;
      p_bg=NULL;
      (*p_pmap)[m_name]=mapname;
      return true;
    }
  }
  }
  if (msg_LevelIsTracking()) {
    msg_Tracking()<<ComixLogo()<<" initialized '"<<m_name<<"', ";
    p_bg->PrintStatistics(msg->Tracking(),0);
  }
  if (!m_hasmapfile) p_bg->WriteOutAmpFile(m_name);
  p_umprocs->push_back(this);
  return false;
}

bool COMIX::Single_Process::GeneratePoint()
{
  SetZero();
  if (m_pinfo.m_cls==cls::sum) {
    m_zero=false;
    return true;
  }
  m_zero=true;
  if (p_map!=NULL && m_lookup && p_map->m_lookup) 
    return !(m_zero=p_map->m_zero);
  if (!p_int->ColorIntegrator()->GeneratePoint()) return false;
  if (p_int->HelicityIntegrator()!=NULL && 
      !p_int->HelicityIntegrator()->GeneratePoint()) return false;
  m_zero=false;
  return true;
}

Weights_Map COMIX::Single_Process::Differential(
    const Cluster_Amplitude &ampl,
    Variations_Mode varmode,
    int mode)
{
  DEBUG_FUNC(Name());
  m_zero=false;
  if ((mode&128)==0) p_int->ColorIntegrator()->SetPoint(&ampl);
  return PHASIC::Process_Base::Differential(ampl,varmode,mode);
}

double COMIX::Single_Process::SetZero()
{
  if (m_pinfo.m_fi.NLOType()&nlo_type::rsub) {
    const NLO_subevtlist &rsubs(p_map?m_subs:p_bg->SubEvts());
    for (size_t i(0);i<rsubs.size();++i) rsubs[i]->Reset();
  }
  m_last["ME"]=0.0;
  return m_w=m_dxs=m_lastxs=0.0;
}

double COMIX::Single_Process::Partonic(const Vec4D_Vector &p, int mode)
{
  Single_Process *sp(p_map!=NULL?p_map:this);
  if (mode==1) {
    UpdateKPTerms(mode);
    return m_lastxs=m_dxs+KPTerms(mode,p_int->ISR()->PDF(0),
                                       p_int->ISR()->PDF(1));
  }
  if (m_zero || !Selector()->Result()) return m_lastxs;
  for (size_t i(0);i<m_nin+m_nout;++i) m_p[i]=p[i];
  if (p_map!=NULL && m_lookup && p_map->m_lookup) {
    m_dxs=p_map->m_dxs;
    m_w=p_map->m_w;
    if (m_pinfo.m_fi.NLOType()&nlo_type::rsub) {
      const NLO_subevtlist &rsubs(p_map->p_bg->SubEvts());
      for (size_t i(0);i<rsubs.size();++i) {
	m_subs[i]->CopyXSData(rsubs[i]);
	for (Cluster_Amplitude *campl(m_subs[i]->p_ampl);
	     campl;campl=campl->Next())
	  for (size_t i(0);i<campl->Legs().size();++i)
	    campl->Leg(i)->SetFlav(ReMap(campl->Leg(i)->Flav(),0));
      }
    }
  }
  else {
    if (m_pinfo.m_fi.NLOType()&nlo_type::rsub &&
        !sp->p_bg->RSTrigger(Selector(),m_mcmode))
      return m_lastxs=m_dxs=0.0;
    sp->p_scale->CalculateScale(p);
    if (m_pinfo.m_cls==cls::sample) {
      m_dxs=sp->p_bg->Differential();
      m_w=p_int->ColorIntegrator()->GlobalWeight();
    }
    else {
      m_w=1.0;
      m_dxs=0.0;
      ComputeHardMatrix(2);
      for (size_t i(0); i<m_cols.m_perms.size(); ++i)
	for (size_t j(0); j<m_cols.m_perms.size(); ++j)
	  m_dxs+=((*p_hc)[i][j]*m_cols.m_colfacs[i][j]).real();
    }
    if (p_int->HelicityIntegrator()!=NULL) 
      m_w*=p_int->HelicityIntegrator()->Weight();
    int isb(m_dxs==sp->p_bg->Born());
    double kb(sp->p_bg->Born()?sp->KFactor(1|2):1.0);
    double kf(m_mewgtinfo.m_K=isb?kb:
	      sp->KFactor(sp->p_bg->Born()?0:2));
    m_mewgtinfo.m_B=sp->p_bg->Born()*kb/kf;
    m_dxs+=sp->p_bg->Born()*(kb/kf-1.0);
    m_w*=kf;
    m_dxs*=m_w;
    if (m_pinfo.m_fi.NLOType()&nlo_type::rsub) {
      const NLO_subevtlist &rsubs(sp->p_bg->SubEvts());
      for (size_t i(0);i<rsubs.size()-1;++i) {
	rsubs[i]->Mult(sp->p_bg->KT2Trigger(rsubs[i],m_mcmode));
	rsubs[i]->m_K=sp->KFactorSetter()->KFactor(*rsubs[i]);
	rsubs[i]->Mult(rsubs[i]->m_K/kf);
      }
      rsubs.back()->m_K=kf;
      if (p_map==NULL) p_bg->SubEvts().MultME(m_w);
      else {
	for (size_t i(0);i<rsubs.size();++i)
	  m_subs[i]->CopyXSData(rsubs[i]);
	m_subs.MultME(m_w);
      }
    }
  }
  UpdateKPTerms(mode);
  double kpterms(KPTerms(mode,p_int->ISR()->PDF(0),p_int->ISR()->PDF(1)));
  FillMEWeights(m_mewgtinfo);
  m_mewgtinfo*=m_w;
  m_mewgtinfo.m_KP=kpterms;
  return m_lastxs=m_dxs+kpterms;
}

const Hard_Matrix *COMIX::Single_Process::ComputeHardMatrix
(Cluster_Amplitude *const ampl,const int mode)
{
  if (p_map!=NULL) return p_map->ComputeHardMatrix(ampl,mode);
  DEBUG_FUNC(Name()<<", mode = "<<mode);
  msg_Debugging()<<*ampl<<"\n";
  Vec4D_Vector p(ampl->Legs().size());
  for (size_t i(0);i<p.size();++i)
    p[i]=i<ampl->NIn()?-ampl->Leg(i)->Mom():ampl->Leg(i)->Mom();
  p_bg->SetMomenta(p);
  std::vector<double> s(ScaleSetter(1)->Scales().size(),0.0);
  s[stp::fac]=ampl->MuF2();
  s[stp::ren]=ampl->MuR2();
  s[stp::res]=ampl->MuQ2();
  if (s.size()>stp::size+stp::res) s[stp::size+stp::res]=ampl->KT2();
  SetFixedScale(s);
  ScaleSetter(1)->CalculateScale(p);
  if (!(mode&1)) ComputeHardMatrix(mode);
  else {
    std::vector<int> ci(ampl->Legs().size()), cj(ci);
    for (size_t i(0);i<ampl->Legs().size();++i) {
      ci[i]=ampl->Leg(i)->Col().m_i;
      cj[i]=ampl->Leg(i)->Col().m_j;
    }
    p_hc->resize(1,std::vector<Complex>(1));
    (*p_hc)[0][0]=p_bg->Differential(ci,cj,(mode&2)?-1:1);
    return p_hc;
  }
  SetFixedScale(std::vector<double>());
  return p_hc;
}

void COMIX::Single_Process::ComputeHardMatrix(const int mode)
{
  std::vector<Spin_Amplitudes> hc;
  for (size_t k(0);k<m_cols.m_perms.size();++k) {
    const std::vector<int> &perm(m_cols.m_perms[k]);
#ifdef DEBUG__BG
    msg_Debugging()<<"Permutation "<<perm<<std::endl;
#endif
    PHASIC::Int_Vector ci(m_nin+m_nout,0), cj(m_nin+m_nout,0);
    int idx(0);
    for (size_t j(0);j<perm.size();++j) {
      Flavour fl(m_flavs[perm[j]]);
      if (perm[j]<m_nin) fl=fl.Bar();
      int cc(fl.StrongCharge());
      if (cc<0 || cc==8) cj[perm[j]]=idx;
      if (cc>0 || cc==8) ci[perm[j]]=++idx;
    }
    if (m_flavs[perm[0]].StrongCharge()==8) cj[perm[0]]=idx;
    double me(p_bg->Differential(ci,cj,(mode&2)?-1:1));
    std::vector<Spin_Amplitudes> amps;
    std::vector<std::vector<Complex> > cols;
    FillAmplitudes(amps,cols);
    hc.push_back(amps.front());
#ifdef DEBUG__BG
    msg_Debugging()<<"Add permutation "<<k<<" "<<perm
		   <<" ( me2 = "<<me<<" ) {\n";
    {
      msg_Indent();
      msg_Debugging()<<"ci = "<<ci<<"\n";
      msg_Debugging()<<"cj = "<<cj<<"\n";
      msg_Debugging()<<hc.back();
    }
    msg_Debugging()<<"}\n";
#endif
  }
  double w(p_bg->ISSymmetryFactor()*
	   p_bg->FSSymmetryFactor());
  p_hc->resize(m_cols.m_perms.size(),
	       std::vector<Complex>
	       (m_cols.m_perms.size()));
  size_t np(p_hc->size()), nh(hc.front().size());
  for (size_t i(0);i<np;++i)
    for (size_t j(i);j<np;++j) {
      (*p_hc)[i][j]=Complex(0.0,0.0);
      for (size_t k(0);k<nh;++k)
	(*p_hc)[i][j]+=hc[i][k]*std::conj(hc[j][k]);
      (*p_hc)[j][i]=std::conj((*p_hc)[i][j]/=w);
    }
}

void COMIX::Single_Process::UpdateKPTerms(const int mode)
{
  DEBUG_FUNC("");
  m_x[0]=m_x[1]=1.0;
  if (!(m_pinfo.m_fi.NLOType()&nlo_type::vsub)) return;
  if (!((m_itype&cs_itype::K) || (m_itype&cs_itype::P))) return;
  const Vec4D &p0(p_int->Momenta()[0]), &p1(p_int->Momenta()[1]);
  double eta0(p0[3]>0.0?p0.PPlus()/rpa->gen.PBeam(0).PPlus():
	      p0.PMinus()/rpa->gen.PBeam(1).PMinus());
  double eta1(p1[3]<0.0?p1.PMinus()/rpa->gen.PBeam(1).PMinus():
	      p1.PPlus()/rpa->gen.PBeam(0).PPlus());
  Single_Process *sp(p_map!=NULL?p_map:this);
  double w(1.0);
  bool map(p_map!=NULL && m_lookup && p_map->m_lookup);
  if (p_int->ISR()->PDF(0) && p_int->ISR()->PDF(0)->Contains(m_flavs[0])) {
    m_x[0]=map?p_map->m_x[0]:eta0+p_ismc->ERan("z_1")*(1.0-eta0);
    w*=(1.0-eta0);
  }
  if (p_int->ISR()->PDF(1) && p_int->ISR()->PDF(1)->Contains(m_flavs[1])) {
    m_x[1]=map?p_map->m_x[1]:eta1+p_ismc->ERan("z_2")*(1.-eta1);
    w*=(1.0-eta1);
  }
  p_kpterms->Calculate(p_int->Momenta(),sp->p_bg->DSij(),
		       m_x[0],m_x[1],eta0,eta1,w);
}

double COMIX::Single_Process::KPTerms
(const int mode, PDF::PDF_Base *pdfa, PDF::PDF_Base *pdfb, double sf)
{
  if (!(m_pinfo.m_fi.NLOType()&nlo_type::vsub)) return 0.0;
  const Vec4D &p0(p_int->Momenta()[0]), &p1(p_int->Momenta()[1]);
  double eta0(p0[3]>0.0?p0.PPlus()/rpa->gen.PBeam(0).PPlus():
	      p0.PMinus()/rpa->gen.PBeam(1).PMinus());
  double eta1(p1[3]<0.0?p1.PMinus()/rpa->gen.PBeam(1).PMinus():
	      p1.PPlus()/rpa->gen.PBeam(0).PPlus());
  double muf2(ScaleSetter(1)->Scale(stp::fac,1));
  return m_w*p_kpterms->Get(pdfa,pdfb,m_x[0],m_x[1],eta0,eta1,
			    muf2,muf2,sf,sf,m_flavs[0],m_flavs[1]);
}

void COMIX::Single_Process::FillMEWeights(ME_Weight_Info &wgtinfo) const
{
  wgtinfo.m_swap=m_p[0][3]<m_p[1][3];
  wgtinfo.m_y1=m_x[wgtinfo.m_swap];
  wgtinfo.m_y2=m_x[1-wgtinfo.m_swap];
  (p_map?p_map:this)->p_bg->FillMEWeights(wgtinfo);
  if (p_kpterms) p_kpterms->FillMEwgts(wgtinfo);
}

int COMIX::Single_Process::PerformTests()
{
  return Tests();
}

bool COMIX::Single_Process::Tests()
{
  msg_Debugging()<<METHOD<<"(): Test '"<<m_name<<"'."<<std::endl;
  if (p_map!=NULL) {
    p_int->SetColorIntegrator(p_map->Integrator()->ColorIntegrator());
    p_int->SetHelicityIntegrator(p_map->Integrator()->HelicityIntegrator());
    p_psgen=p_map->p_psgen;
    return true;
  }
  if (p_bg==NULL) {
    msg_Error()<<METHOD<<"(): No amplitude for '"<<Name()<<"'"<<std::endl;
    return false;
  }
  if (m_gpath.length()>0) {
    std::string script("/plot_graphs");
    if (!FileExists(rpa->gen.Variable("SHERPA_CPP_PATH")+script))
      Copy(rpa->gen.Variable("SHERPA_SHARE_PATH")+script,
           rpa->gen.Variable("SHERPA_CPP_PATH")+script);
    m_gpath+=std::string("/Comix");
    MakeDir(m_gpath,448);
    p_bg->WriteOutGraphs(m_gpath+"/"+ShellName(m_name)+".tex");
  }
  if (p_int->HelicityScheme()==hls::sample) {
    p_int->SetHelicityIntegrator(std::make_shared<Helicity_Integrator>());
    p_bg->SetHelicityIntegrator(&*p_int->HelicityIntegrator());
    Flavour_Vector fl(m_nin+m_nout);
    for (size_t i(0);i<fl.size();++i) fl[i]=m_flavs[i];
    if (!p_int->HelicityIntegrator()->Construct(fl)) return false;
  }
  p_int->SetColorIntegrator(std::make_shared<Color_Integrator>());
  p_bg->SetColorIntegrator(&*p_int->ColorIntegrator());
  Idx_Vector ids(m_nin+m_nout,0);
  Int_Vector acts(m_nin+m_nout,0), types(m_nin+m_nout,0);
  for (size_t i(0);i<ids.size();++i) {
    ids[i]=i;
    acts[i]=m_flavs[i].Strong();
    if (acts[i]) {
      if (m_flavs[i].StrongCharge()==8) types[i]=0;
      else if (m_flavs[i].IsAnti()) types[i]=i<m_nin?1:-1;
      else types[i]=i<m_nin?-1:1;
    }
  }
  if (!p_int->ColorIntegrator()->
      ConstructRepresentations(ids,types,acts)) return false;
  const DecayInfo_Vector &dinfos(p_bg->DecayInfos());
  std::vector<size_t> dids(dinfos.size());
  acts.resize(dids.size());
  types.resize(dids.size());
  for (size_t i(0);i<dids.size();++i) {
    dids[i]=dinfos[i]->m_id;
    acts[i]=dinfos[i]->m_fl.Strong();
    if (acts[i]) {
      if (dinfos[i]->m_fl.StrongCharge()==8) types[i]=0;
      else if (dinfos[i]->m_fl.IsAnti()) types[i]=-1;
      else types[i]=1;
    }
  }
  p_int->ColorIntegrator()->SetDecayIds(dids,types,acts);
  Phase_Space_Handler::TestPoint(&m_p.front(),&Info(),Generator(),1);
  bool res(p_bg->GaugeTest(m_p));
  if (!res) {
    msg_Info()<<METHOD<<"(): Gauge test failed for '"
	      <<m_name<<"'."<<std::endl;
  }
  else if (!msg_LevelIsTracking()) msg_Info()<<"."<<std::flush;
  return res;
}

bool COMIX::Single_Process::Trigger(const ATOOLS::Vec4D_Vector &p)
{
  DEBUG_FUNC(m_pinfo.m_fi.NLOType());
  if (m_zero) return false;
  if (p_map!=NULL && m_lookup && p_map->m_lookup)
    return Selector()->Result();
  if (m_pinfo.m_fi.NLOType()&nlo_type::rsub) {
    Amplitude *bg(p_map!=NULL?p_map->p_bg:p_bg);
    Selector()->SetResult(1);
    if (bg->SetMomenta(p)) return true;
    Selector()->SetResult(0);
    return false;
  }
  (p_map!=NULL?p_map->p_bg:p_bg)->SetMomenta(p);
  return Selector()->Trigger(p);
}

void COMIX::Single_Process::InitPSGenerator(const size_t &ismode)
{
  if (p_map!=NULL) {
    p_psgen=p_map->p_psgen;
    if (p_psgen == nullptr)
      p_psgen = std::make_shared<PS_Generator>(p_map);
  } else {
    p_psgen = std::make_shared<PS_Generator>(this);
  }
}

void COMIX::Single_Process::ConstructPSVertices(PS_Generator *ps)
{
  if (m_psset.find(ps)!=m_psset.end()) return;
  m_psset.insert(ps);
  if (p_bg!=NULL) ps->Construct(p_bg,GetSubevtList());
  else p_map->ConstructPSVertices(ps);
}

Amplitude *COMIX::Single_Process::GetAmplitude() const
{
  if (p_map) return p_map->p_bg;
  return p_bg;
}

bool COMIX::Single_Process::FillIntegrator(Phase_Space_Handler *const psh)
{
  bool res(COMIX::Process_Base::FillIntegrator(psh));
  if ((m_pinfo.m_fi.NLOType()&nlo_type::vsub) && p_ismc) {
    p_ismc->AddERan("z_1");
    p_ismc->AddERan("z_2");
  }
  return res;
}

Flavour COMIX::Single_Process::ReMap
(const Flavour &fl,const size_t &id) const
{
  if (p_map==NULL) return fl;
  Flavour_Map::const_iterator fit(m_fmap.find(fl));
  if (fit!=m_fmap.end()) return fit->second;
  fit=m_fmap.find(fl.Bar());
  if (fit!=m_fmap.end()) return fit->second.Bar();
  THROW(fatal_error,"Invalid flavour '"+ToString(fl)+"'");
  return fl;
}

bool COMIX::Single_Process::Combinable
(const size_t &idi,const size_t &idj)
{
  Combination_Set::const_iterator 
    cit(m_ccombs.find(std::pair<size_t,size_t>(idi,idj)));
  return cit!=m_ccombs.end();
}

const Flavour_Vector &COMIX::Single_Process::
CombinedFlavour(const size_t &idij)
{
  CFlavVector_Map::const_iterator fit(m_cflavs.find(idij));
  if (fit==m_cflavs.end()) THROW(fatal_error,"Invalid request");
  return fit->second;
}

void COMIX::Single_Process::FillAmplitudes
(std::vector<Spin_Amplitudes> &amps,
 std::vector<std::vector<Complex> > &cols)
{
  (p_map?p_map->p_bg:p_bg)->FillAmplitudes(amps,cols);
}

NLO_subevtlist *COMIX::Single_Process::GetSubevtList()
{
  if (m_pinfo.m_fi.NLOType()&nlo_type::rsub)
    return &(p_map?m_subs:p_bg->SubEvts());
  return NULL;
}

NLO_subevtlist *COMIX::Single_Process::GetRSSubevtList()
{
  return GetSubevtList();
}

void COMIX::Single_Process::SetScale(const Scale_Setter_Arguments &args)
{
  PHASIC::Single_Process::SetScale(args);
  Scale_Setter_Base *scs(p_map?p_map->p_scale:p_scale);
  NLO_subevtlist *subs(GetSubevtList());
  if (subs) {
    for (size_t i(0);i<subs->size()-1;++i)
      static_cast<Single_Dipole_Term*>
	((*subs)[i]->p_proc)->SetScaleSetter(scs);
  }
}

void COMIX::Single_Process::SetShower(PDF::Shower_Base *const ps)
{
  PHASIC::Single_Process::SetShower(ps);
  NLO_subevtlist *subs(GetSubevtList());
  if (subs) {
    for (size_t i(0);i<subs->size()-1;++i)
      static_cast<Single_Dipole_Term*>
	((*subs)[i]->p_proc)->SetShower(ps);
  }
}

void COMIX::Single_Process::SetNLOMC(PDF::NLOMC_Base *const mc)
{
  PHASIC::Single_Process::SetNLOMC(mc);
  if (p_bg) p_bg->SetNLOMC(mc);
  if (p_kpterms) p_kpterms->SetNLOMC(mc);
}

size_t COMIX::Single_Process::SetMCMode(const size_t mcmode)
{
  size_t cmcmode(m_mcmode);
  m_mcmode=mcmode;
  NLO_subevtlist *subs(GetSubevtList());
  if (subs) {
    for (size_t i(0);i<subs->size()-1;++i)
      static_cast<Single_Dipole_Term*>
	((*subs)[i]->p_proc)->SetMCMode(mcmode);
  }
  return cmcmode;
}

void COMIX::Single_Process::SetLookUp(const bool lookup)
{
  m_lookup=lookup;
  NLO_subevtlist *subs(GetSubevtList());
  if (subs) {
    for (size_t i(0);i<subs->size()-1;++i)
      static_cast<Single_Dipole_Term*>
	((*subs)[i]->p_proc)->SetLookUp(m_lookup);
  }
  if (p_loop && lookup==0) p_loop->SwitchMode(lookup);
}
