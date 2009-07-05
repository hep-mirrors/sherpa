#include "COMIX/Main/Single_Process.H"

#include "COMIX/Main/Process_Group.H"
#include "COMIX/Amplitude/Matrix_Element.H"
#include "PDF/Main/ISR_Handler.H"
#include "COMIX/Phasespace/PS_Generator.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Main/Color_Integrator.H"
#include "PHASIC++/Main/Helicity_Integrator.H"
#include "PHASIC++/Scales/KFactor_Setter_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"

#include <sys/stat.h>

// #define TIME_CDBG_ME

using namespace COMIX;
using namespace PHASIC;
using namespace ATOOLS;

COMIX::Single_Process::Single_Process(): COMIX::Process_Base(this),
  p_bg(NULL), p_map(NULL), m_etime(0.0), m_en(0.0) {}

COMIX::Single_Process::~Single_Process()
{
  if (p_bg!=NULL) delete p_bg;
}

bool COMIX::Single_Process::Initialize
(std::map<std::string,std::string> *const pmap,
 std::vector<Single_Process*> *const procs)
{
  m_p.resize(m_nin+m_nout);
  if (!Process_Base::Initialize(pmap,procs)) return false;
  if (p_bg!=NULL) delete p_bg;
  p_bg=NULL;
  if (pmap->find(m_name)!=pmap->end()) {
    std::string mapname((*pmap)[m_name]);
    msg_Debugging()<<METHOD<<"(): Map '"<<m_name<<"' onto '"
		   <<mapname<<"'"<<std::endl;
    if (mapname=="x") return false;
    if (mapname!=m_name) return true;
  }
  msg_Debugging()<<"'"<<m_name<<"' not pre-mapped"<<std::endl;
  p_bg = new Matrix_Element();
  Subprocess_Info info(m_pinfo.m_ii);
  info.Add(m_pinfo.m_fi);
  p_bg->SetDecayInfos(info.GetDecayInfos());
  std::vector<Flavour> flavs(m_nin+m_nout);
  for (size_t i(0);i<m_nin+m_nout;++i) flavs[i]=m_flavs[i];
  if (p_bg->Initialize(m_nin,m_nout,flavs,&*p_model,
		       m_pinfo.m_oew,m_pinfo.m_oqcd)) {
    m_oew=p_bg->MaxOrderEW();
    m_oqcd=p_bg->MaxOrderQCD();
    (*pmap)[m_name]=m_name;
    return true;
  }
  (*pmap)[m_name]="x";
  return false;
}

bool COMIX::Single_Process::MapProcess()
{
  std::string mapname((*p_pmap)[m_name]);
  if (mapname!=m_name) {
    for (size_t i(0);i<p_umprocs->size();++i)
      if ((*p_umprocs)[i]->Name()==mapname) {
	p_mapproc=p_map=(*p_umprocs)[i];
	m_oew=p_map->p_bg->MaxOrderEW();
	m_oqcd=p_map->p_bg->MaxOrderQCD();
	msg_Tracking()<<"Mapped '"<<m_name<<"' -> '"
		      <<mapname<<"'."<<std::endl;
	std::string mapfile(rpa.gen.Variable("SHERPA_CPP_PATH")
			    +"/Process/Comix");
	MakeDir(mapfile,true);
	mapfile+="/"+m_name+".fmap";
	std::ifstream map(mapfile.c_str());
	if (!map.good()) {
	  THROW(fatal_error,"Flavour map file '"+mapfile+"' not found");
	}
	else {
	  while (!map.eof()) {
	    long int src, dest;
	    map>>src>>dest;
	    Flavour ft((kf_code)(abs(src)),src<0);
	    Flavour fb((kf_code)(abs(dest)),dest<0);
	    m_fmap[ft]=fb;
	    msg_Debugging()<<"  fmap '"<<ft<<"' onto '"<<fb<<"'\n";
	  }
	}
	return true;
      }
    THROW(fatal_error,"Map process '"+mapname+"' not found");
  }
  for (size_t i(0);i<p_umprocs->size();++i) {
    msg_Debugging()<<METHOD<<"(): Try mapping '"
		   <<Name()<<"' -> '"<<(*p_umprocs)[i]->Name()<<"'\n";
    if (p_bg->Map((*p_umprocs)[i]->p_bg,m_fmap)) {
      p_mapproc=p_map=(*p_umprocs)[i];
      mapname=p_map->Name();
      msg_Tracking()<<"Mapped '"<<m_name<<"' -> '"
		    <<mapname<<"'."<<std::endl;
      std::string mapfile(rpa.gen.Variable("SHERPA_CPP_PATH")
			  +"/Process/Comix");
      MakeDir(mapfile,true);
      mapfile+="/"+m_name+".fmap";
      std::ofstream map(mapfile.c_str());
      if (map.good()) {
	for (Flavour_Map::const_iterator 
	       fit(m_fmap.begin());fit!=m_fmap.end();++fit) {
	  msg_Debugging()<<"  fmap '"<<fit->first
			 <<"' onto '"<<fit->second<<"'\n";
	  long int src(fit->first), dest(fit->second);
	  map<<src<<" "<<dest<<"\n";
	}
      }
      delete p_bg;
      p_bg=NULL;
      (*p_pmap)[m_name]=mapname;
      return true;
    }
  }
  if (msg_LevelIsTracking()) {
    msg_Tracking()<<ComixLogo()<<" initialized '"<<m_name<<"', ";
    p_bg->PrintStatistics(msg->Tracking(),0);
  }
  p_umprocs->push_back(this);
  return false;
}

bool COMIX::Single_Process::GeneratePoint()
{
  m_zero=true;
  m_lastxs=m_last=0.0;
  if (p_map!=NULL && m_lookup && p_map->m_lookup) 
    return !(m_zero=p_map->m_zero);
  if (!p_int->ColorIntegrator()->GeneratePoint()) return false;
  if (p_int->HelicityIntegrator()!=NULL && 
      !p_int->HelicityIntegrator()->GeneratePoint()) return false;
  return !(m_zero=false);
}

double COMIX::Single_Process::Differential(const Cluster_Amplitude &ampl) 
{
  m_zero=false;
  p_int->ColorIntegrator()->SetPoint(&ampl);
  return PHASIC::Process_Base::Differential(ampl);
}

double COMIX::Single_Process::Differential(const Vec4D_Vector &p) 
{
  if (p_map!=NULL && m_lookup && p_map->m_lookup)
    SetTrigger(p_map->Trigger());
  if (m_zero || !Trigger()) return m_lastxs=m_last=0.0;
  for (size_t i(0);i<m_nin+m_nout;++i) {
    m_p[i]=p[i];
    double psm(m_flavs[i].Mass());
    if (m_p[i][0]<psm) return m_lastxs=m_last=0.0;
    if (i<m_nin && psm==0.0) m_p[i][0]=dabs(m_p[i][3]);
  }
  if (p_map!=NULL && m_lookup && p_map->m_lookup) {
    m_lastxs=p_map->m_lastxs;
  }
  else {
#ifdef TIME_CDBG_ME
    double stime(rpa.gen.Timer().UserTime());
#endif
    m_lastxs=(p_map!=NULL?p_map->p_bg:p_bg)->Differential(m_p);
#ifdef TIME_CDBG_ME
    m_etime+=rpa.gen.Timer().UserTime()-stime;
    ++m_en;
#endif
    m_lastxs*=p_int->ColorIntegrator()->GlobalWeight();
    if (p_int->HelicityIntegrator()!=NULL) 
      m_lastxs*=p_int->HelicityIntegrator()->Weight();
  }
  if (m_lastxs<=0.) return m_lastxs=m_last=0.;
  p_scale->CalculateScale(p_int->PSHandler()->LabPoint());
  if (p_int->ISR() && m_nin==2) { 
    if (p_int->ISR()->On()) {
      p_int->ISR()->MtxLock();
      if (!p_int->ISR()->CalculateWeight
	  (p_scale->Scale(stp::fac))) {
	p_int->ISR()->MtxUnLock();
	return m_last=m_lastlumi=0.0;
      }
      m_lastlumi=p_int->ISR()->Weight(&m_flavs.front()); 
      p_int->ISR()->MtxUnLock();
    }
    else m_lastlumi=1.;
  }
  else m_lastlumi=1.;
  return m_last=m_lastxs*m_lastlumi*KFactor();
}

double COMIX::Single_Process::Differential2() 
{
  if (m_zero || m_lastxs==0.0 || !Trigger()) return 0.0;
  if (p_int->ISR() && m_nin==2) {
    p_scale->CalculateScale2(p_int->PSHandler()->LabPoint());
    if (m_flavs[0]==m_flavs[1] || p_int->ISR()->On()==0) return 0.;
    p_int->ISR()->MtxLock();
    if (!p_int->ISR()->CalculateWeight2
	(p_scale->Scale(stp::fac))) {
      p_int->ISR()->MtxUnLock();
      return 0.0;
    }
    double tmp=m_lastxs*p_int->ISR()->
      Weight2(&m_flavs.front())*KFactor2(); 
    p_int->ISR()->MtxUnLock();
    m_last+=tmp;
    return tmp;
  }
  return 0;
}

bool COMIX::Single_Process::Tests()
{
  msg_Debugging()<<METHOD<<"(): Test '"<<m_name<<"' <- '"
		 <<p_int->PSHandler()->Process()->
    Process()->Name()<<"'."<<std::endl;
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
  MakeDir(m_gpath,448);
  if (m_gpath.length()>0) p_bg->PrintGraphs(m_gpath+"/"+m_name+".tex");
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
    if (!m_flavs[i].IsFermion()) types[i]=0;
    else if (m_flavs[i].IsAnti()) types[i]=i<m_nin?1:-1;
    else types[i]=i<m_nin?-1:1;
  }
  if (!p_int->ColorIntegrator()->ConstructRepresentations(ids,types,acts)) return false;
  Phase_Space_Handler::TestPoint(&m_p.front(),m_nin,m_nout,m_flavs);
  bool res(p_bg->GaugeTest(m_p));
  if (!res) {
    msg_Info()<<METHOD<<"(): Gauge test failed for '"
	      <<m_name<<"'."<<std::endl;
  }
  else if (!msg_LevelIsTracking()) msg_Info()<<"."<<std::flush;
  return res;
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
  if (p_bg!=NULL) ps->Construct(p_bg->GetAmplitude());
}

Amplitude *COMIX::Single_Process::GetAmplitude() const
{
  if (p_bg!=NULL) return p_bg->GetAmplitude();
  return p_map->p_bg->GetAmplitude();
}

Matrix_Element *COMIX::Single_Process::GetME() const
{
  if (p_bg!=NULL) return p_bg;
  return p_map->p_bg;
}

bool COMIX::Single_Process::FillIntegrator(Phase_Space_Handler *const psh)
{
  return Process_Base::FillIntegrator(psh);
}

void COMIX::Single_Process::UpdateIntegrator(Phase_Space_Handler *const psh)
{
  Process_Base::UpdateIntegrator(psh);
}

Flavour COMIX::Single_Process::ReMap(const Flavour &fl) const
{
  if (p_map==NULL) return fl;
  Flavour_Map::const_iterator fit(m_fmap.find(fl));
  if (fit!=m_fmap.end()) return fit->second;
  fit=m_fmap.find(fl.Bar());
  if (fit!=m_fmap.end()) return fit->second.Bar();
  THROW(fatal_error,"Invalid flavour '"+ToString(fl)+"'");
  return fl;
}

void COMIX::Single_Process::SetKFactorOn(const bool on)
{
  PHASIC::Single_Process::SetKFactorOn(on);
  if (p_map) p_map->SetKFactorOn(on);
}

bool COMIX::Single_Process::Combinable
(const size_t &idi,const size_t &idj)
{
  if (p_map) return p_map->Combinable(idi,idj);
  return p_bg->GetAmplitude()->Combinable(idi,idj);
}

Flavour COMIX::Single_Process::
CombinedFlavour(const size_t &idij)
{
  if (p_map) return ReMap(p_map->CombinedFlavour(idij));
  return p_bg->GetAmplitude()->CombinedFlavour(idij);
}
