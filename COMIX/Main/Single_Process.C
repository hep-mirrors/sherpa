#include "COMIX/Main/Single_Process.H"

#include "COMIX/Main/Process_Group.H"
#include "COMIX/Main/Single_Dipole_Term.H"
#include "PDF/Main/ISR_Handler.H"
#include "COMIX/Phasespace/PS_Generator.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "PHASIC++/Main/Color_Integrator.H"
#include "PHASIC++/Main/Helicity_Integrator.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Scales/KFactor_Setter_Base.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/My_MPI.H"

using namespace COMIX;
using namespace PHASIC;
using namespace ATOOLS;

COMIX::Single_Process::Single_Process():
  COMIX::Process_Base(this),
  p_bg(NULL), p_map(NULL),
  p_loop(NULL), p_kpterms(NULL),
  m_checkpoles(false), m_allowmap(true)
{
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
 std::vector<Single_Process*> *const procs)
{
  m_p.resize(m_nin+m_nout);
  if (!COMIX::Process_Base::Initialize(pmap,procs)) return false;
  if (p_bg!=NULL) delete p_bg;
  p_bg=NULL;
  if (pmap->find(m_name)!=pmap->end()) {
    if ((*pmap)[m_name]!="x") THROW(fatal_error,"Internal error");
    return false;
  }
  std::string mapfile(rpa->gen.Variable("SHERPA_CPP_PATH")
		      +"/Process/Comix/"+m_name+".map");
  if (FileExists(mapfile)) {
    msg_Tracking()<<METHOD<<"(): Map file '"<<mapfile<<"' found.\n";
    My_In_File mapstream(mapfile);
    if (mapstream.Open()) {
      std::string cname, mapname;
      *mapstream>>cname>>mapname;
      if (cname!=m_name)
	THROW(fatal_error,"Corrupted map file '"+mapfile+"'");
      (*pmap)[m_name]=mapname;
      if (mapname!=m_name) {
	msg_Debugging()<<METHOD<<"(): Map '"<<m_name
		       <<"' onto '"<<mapname<<"'\n";
	return true;
      }
    }
  }
  msg_Debugging()<<"'"<<m_name<<"' not pre-mapped"<<std::endl;
  p_model->GetCouplings(m_cpls);
  p_bg = new Amplitude();
  double isf(m_pinfo.m_ii.ISSymmetryFactor());
  double fsf(m_pinfo.m_fi.FSSymmetryFactor());
  Subprocess_Info info(m_pinfo.m_ii);
  info.Add(m_pinfo.m_fi);
  p_bg->SetDecayInfos(m_decins);
  p_bg->SetNoDecays(m_pinfo.m_nodecs);
  std::vector<Flavour> flavs(m_nin+m_nout);
  for (size_t i(0);i<m_nin+m_nout;++i) flavs[i]=m_flavs[i];
  int smode(0);
  if (m_pinfo.m_fi.NLOType()&nlo_type::rsub) smode=1;
  if (m_pinfo.m_fi.NLOType()&nlo_type::vsub) {
    smode=2;
    if (m_pinfo.m_fi.NLOType()&nlo_type::born) smode|=4;
  }
  if (m_pinfo.m_fi.NLOType()&nlo_type::loop) {
    smode|=16;
    if (m_checkpoles) smode|=8;
  }
  std::vector<int> mincpl(m_pinfo.m_mincpl.size());
  std::vector<int> maxcpl(m_pinfo.m_maxcpl.size());
  for (size_t i(0);i<mincpl.size();++i) mincpl[i]=m_pinfo.m_mincpl[i]*2;
  for (size_t i(0);i<maxcpl.size();++i) maxcpl[i]=m_pinfo.m_maxcpl[i]*2;
  if (p_bg->Initialize(m_nin,m_nout,flavs,isf,fsf,&*p_model,
		       &m_cpls,smode,maxcpl,mincpl,
		       m_pinfo.m_ntchan,m_pinfo.m_mtchan,m_name)) {
    nlo_type::code nlot(nlo_type::loop|nlo_type::vsub);
    m_maxcpl[1]=p_bg->MaxCpl()[1]/2.0+((m_pinfo.m_fi.m_nloewtype&nlot)?1:0);
    m_mincpl[1]=p_bg->MinCpl()[1]/2.0+((m_pinfo.m_fi.m_nloewtype&nlot)?1:0);
    m_maxcpl[0]=p_bg->MaxCpl()[0]/2.0+((m_pinfo.m_fi.m_nloqcdtype&nlot)?1:0);
    m_mincpl[0]=p_bg->MinCpl()[0]/2.0+((m_pinfo.m_fi.m_nloqcdtype&nlot)?1:0);
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
      for (size_t i(m_nin);i<m_nin+m_nout;++i)
	if (flavs[i].Mass()) massive=1;
      p_kpterms = new KP_Terms(this,massive);
      p_bg->DInfo()->SetMassive(massive);
      p_kpterms->SetAlpha(p_bg->DInfo()->AMax(0),
			  p_bg->DInfo()->AMax(2),
			  p_bg->DInfo()->AMax(1),
			  p_bg->DInfo()->AMax(3));
      m_mewgtinfo.m_type|=mewgttype::VI|mewgttype::KP;
      if (smode&4) m_mewgtinfo.m_type|=mewgttype::B;
    }
    if (smode&16) {
      Process_Info cinfo(m_pinfo);
      cinfo.m_fi.m_nloqcdtype=nlo_type::loop;
      cinfo.m_maxcpl[1]=p_bg->MaxCpl()[1]/2.0+
	((cinfo.m_fi.m_nloewtype&nlo_type::loop)?1:0);
      cinfo.m_mincpl[1]=p_bg->MinCpl()[1]/2.0+
	((cinfo.m_fi.m_nloewtype&nlo_type::loop)?1:0);
      cinfo.m_maxcpl[0]=p_bg->MaxCpl()[0]/2.0+
	((cinfo.m_fi.m_nloqcdtype&nlo_type::loop)?1:0);
      cinfo.m_mincpl[0]=p_bg->MinCpl()[0]/2.0+
	((cinfo.m_fi.m_nloqcdtype&nlo_type::loop)?1:0);
      p_loop = PHASIC::Virtual_ME2_Base::GetME2(cinfo);
      if (p_loop==NULL) {
	msg_Error()<<METHOD<<"(): "<<cinfo<<"\n";
	THROW(not_implemented,"No virtual ME for this process");
      }
      p_loop->SetCouplings(m_cpls);
      p_loop->SetNorm(1.0/(isf*fsf));
      m_mewgtinfo.m_type|=mewgttype::VI;
      int helpi;
      Data_Reader reader(" ",";","!","=");
      reader.AddComment("#");
      reader.SetInputPath(rpa->GetPath());
      reader.SetInputFile(rpa->gen.Variable("ME_DATA_FILE"));
      if (reader.ReadFromFile(helpi,"CHECK_POLES")) {
	m_checkpoles=helpi;
	msg_Tracking()<<"Set pole check mode "<<m_checkpoles<<".\n";
      }
    }
    p_bg->SetLoopME(p_loop);
    (*pmap)[m_name]=m_name;
    return true;
  }
#ifdef USING__MPI
  if (MPI::COMM_WORLD.Get_rank()==0) {
#endif
  if (Parent()->Name()==Name()) return false;
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
#ifdef USING__MPI
  }
#endif
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
	  p_kpterms = new KP_Terms
	    (p_map,p_map->p_kpterms->MassKern()?1:0);
	  p_kpterms->SetAlpha(p_map->p_bg->DInfo()->AMax(0),
			      p_map->p_bg->DInfo()->AMax(2),
			      p_map->p_bg->DInfo()->AMax(1),
			      p_map->p_bg->DInfo()->AMax(3));
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
	return true;
      }
    THROW(fatal_error,"Map process '"+mapname+"' not found");
  }
  std::string ampfile(rpa->gen.Variable("SHERPA_CPP_PATH")
		      +"/Process/Comix/"+m_name+".map");
  if (!FileExists(ampfile) && m_allowmap &&
      m_pinfo.m_special.find("MapOff")==std::string::npos) {
  for (size_t i(0);i<p_umprocs->size();++i) {
    msg_Debugging()<<METHOD<<"(): Try mapping '"
		   <<Name()<<"' -> '"<<(*p_umprocs)[i]->Name()<<"'\n";
    if (p_bg->Map(*(*p_umprocs)[i]->p_bg,m_fmap)) {
      p_mapproc=p_map=(*p_umprocs)[i];
      if (p_kpterms) {
	delete p_kpterms;
	p_kpterms = new KP_Terms
	  (p_map,p_map->p_kpterms->MassKern()?1:0);
	p_kpterms->SetAlpha(p_map->p_bg->DInfo()->AMax(0),
			    p_map->p_bg->DInfo()->AMax(2),
			    p_map->p_bg->DInfo()->AMax(1),
			    p_map->p_bg->DInfo()->AMax(3));
      }
      mapname=p_map->Name();
      msg_Tracking()<<"Mapped '"<<m_name<<"' -> '"
		    <<mapname<<"'."<<std::endl;
#ifdef USING__MPI
      if (MPI::COMM_WORLD.Get_rank()==0) {
#endif
      std::string mapfile(rpa->gen.Variable("SHERPA_CPP_PATH")
			  +"/Process/Comix/"+m_name+".map");
      if (!FileExists(mapfile)) {
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
#ifdef USING__MPI
      }
#endif
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
  p_bg->WriteOutAmpFile(m_name);
  p_umprocs->push_back(this);
  return false;
}

bool COMIX::Single_Process::GeneratePoint()
{
  SetZero();
  m_zero=true;
  if (p_map!=NULL && m_lookup && p_map->m_lookup) 
    return !(m_zero=p_map->m_zero);
  if (!p_int->ColorIntegrator()->GeneratePoint()) return false;
  if (p_int->HelicityIntegrator()!=NULL && 
      !p_int->HelicityIntegrator()->GeneratePoint()) return false;
  m_zero=false;
  return true;
}

double COMIX::Single_Process::Differential
(const Cluster_Amplitude &ampl,int mode) 
{
  DEBUG_FUNC(Name());
  m_zero=false;
  p_int->ColorIntegrator()->SetPoint(&ampl);
  return PHASIC::Process_Base::Differential(ampl,mode);
}

double COMIX::Single_Process::SetZero()
{
  if (m_pinfo.m_fi.NLOType()&nlo_type::rsub) {
    const NLO_subevtlist &rsubs(p_map?m_subs:p_bg->SubEvts());
    for (size_t i(0);i<rsubs.size();++i) rsubs[i]->Reset();
  }
  return m_w=m_dxs=m_lastxs=m_last=0.0;
}

double COMIX::Single_Process::Partonic
(const Vec4D_Vector &p,const int mode) 
{
  Single_Process *sp(p_map!=NULL?p_map:this);
  if (mode==1)
    return m_lastxs=m_dxs+m_w*GetKPTerms(m_flavs,mode);
  if (m_zero || !Selector()->Result()) return m_lastxs;
  for (size_t i(0);i<m_nin+m_nout;++i) {
    m_p[i]=p[i];
    double psm(m_flavs[i].Mass());
    if (m_p[i][0]<psm) return m_lastxs;
  }
  if (p_map!=NULL && m_lookup && p_map->m_lookup) {
    m_dxs=p_map->m_dxs;
    m_w=p_map->m_w;
    if (m_pinfo.m_fi.NLOType()&nlo_type::rsub) {
      const NLO_subevtlist &rsubs(p_map->p_bg->SubEvts());
      for (size_t i(0);i<rsubs.size();++i)
	m_subs[i]->CopyXSData(rsubs[i]);
    }
  }
  else {
    if (m_pinfo.m_fi.NLOType()&nlo_type::rsub &&
	!sp->p_bg->JetTrigger(Selector(),m_mcmode))
      return m_lastxs=m_dxs=0.0;
    sp->p_scale->CalculateScale(p);
    m_dxs=sp->p_bg->Differential();
    m_w=p_int->ColorIntegrator()->GlobalWeight();
    if (p_int->HelicityIntegrator()!=NULL) 
      m_w*=p_int->HelicityIntegrator()->Weight();
    int isb(m_dxs==sp->p_bg->Born());
    double kb(sp->p_bg->Born()?sp->KFactor(1):1.0);
    double kf(m_mewgtinfo.m_K=isb?kb:sp->KFactor());
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
  double kpterms(m_w*GetKPTerms(m_flavs,mode));
  FillMEWeights(m_mewgtinfo);
  m_mewgtinfo*=m_w;
  m_mewgtinfo.m_KP=kpterms;
  return m_lastxs=m_dxs+kpterms;
}

double COMIX::Single_Process::GetKPTerms
(const Flavour_Vector &fl,const int mode)
{
  m_x[0]=m_x[1]=1.0;
  if (!(m_pinfo.m_fi.NLOType()&nlo_type::vsub)) return 0.0;
  if (mode==0) {
    Single_Process *sp(p_map!=NULL?p_map:this);
    double eta0=1.0, eta1=1.0;
    double w=sp->p_bg->Coupling(0)/(2.0*M_PI);
    bool map(p_map!=NULL && m_lookup && p_map->m_lookup);
    if (m_flavs[0].Strong()) {
      eta0=p_int->ISR()->X1();
      m_x[0]=map?p_map->m_x[0]:eta0+ran->Get()*(1.0-eta0);
      w*=(1.0-eta0);
    }
    if (m_flavs[1].Strong()) {
      eta1=p_int->ISR()->X2();
      m_x[1]=map?p_map->m_x[1]:eta1+ran->Get()*(1.-eta1);
      w*=(1.0-eta1);
    }
    p_kpterms->SetDSij(sp->p_bg->DSij());
    p_kpterms->Calculate
      (p_int->Momenta(),m_x[0],m_x[1],eta0,eta1,w);
  }
  double eta0(0.0), eta1(0.0);
  if (mode==0) {
    eta0=p_int->Momenta()[0].PPlus()/rpa->gen.PBeam(0).PPlus();
    eta1=p_int->Momenta()[1].PMinus()/rpa->gen.PBeam(1).PMinus();
  }
  else {
    eta0=p_int->Momenta()[0].PPlus()/rpa->gen.PBeam(1).PMinus();
    eta1=p_int->Momenta()[1].PMinus()/rpa->gen.PBeam(0).PPlus();
  }
  return p_kpterms->Get(m_x[0],m_x[1],eta0,eta1,fl,mode);
}

void COMIX::Single_Process::FillMEWeights(ME_Weight_Info &wgtinfo) const
{
  wgtinfo.m_y1=m_x[0];
  wgtinfo.m_y2=m_x[1];
  (p_map?p_map:this)->p_bg->FillMEWeights(wgtinfo);
  if (p_kpterms) p_kpterms->FillMEwgts(wgtinfo);
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
      Copy(rpa->gen.Variable("SHERPA_SHARE_PATH")+script+".sh",
           rpa->gen.Variable("SHERPA_CPP_PATH")+script);
    m_gpath+=std::string("/Comix");
    MakeDir(m_gpath,448);
    p_bg->WriteOutGraphs(m_gpath+"/"+ShellName(m_name)+".tex");
  }
  if (p_int->HelicityScheme()==hls::sample) {
    p_int->SetHelicityIntegrator(new Helicity_Integrator());
    p_bg->SetHelicityIntegrator(&*p_int->HelicityIntegrator());
    Flavour_Vector fl(m_nin+m_nout);
    for (size_t i(0);i<fl.size();++i) fl[i]=m_flavs[i];
    if (!p_int->HelicityIntegrator()->Construct(fl)) return false;
  }
  p_int->SetColorIntegrator(new Color_Integrator());
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
  if (m_zero) return false;
  if (p_map!=NULL && m_lookup && p_map->m_lookup)
    return Selector()->Result();
  if (m_pinfo.m_fi.NLOType()&nlo_type::rsub) {
    if (!Selector()->NoJetTrigger(p)) return false;
    Amplitude *bg(p_map!=NULL?p_map->p_bg:p_bg);
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
    if (p_psgen==NULL) p_psgen = new PS_Generator(p_map);
  }
  else {
    p_psgen = new PS_Generator(this);
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
  if (p_bg!=NULL) return p_bg;
  return p_map->p_bg;
}

bool COMIX::Single_Process::FillIntegrator(Phase_Space_Handler *const psh)
{
  return COMIX::Process_Base::FillIntegrator(psh);
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
  if (p_map) return p_map->Combinable(idi,idj);
  return p_bg->Combinable(idi,idj);
}

const Flavour_Vector &COMIX::Single_Process::
CombinedFlavour(const size_t &idij)
{
  if (p_map) {
    CFlavVector_Map::const_iterator fit(m_cfmap.find(idij));
    if (fit!=m_cfmap.end()) return fit->second;
    Flavour_Vector cf(p_map->CombinedFlavour(idij));
    for (size_t i(0);i<cf.size();++i) cf[i]=ReMap(cf[i],0);
    m_cfmap[idij]=cf;
    return m_cfmap[idij];
  }
  return p_bg->CombinedFlavour(idij);
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
}
