#include "COMIX/Phasespace/PS_Channel.H"

#include "COMIX/Main/Process_Base.H"
#include "COMIX/Phasespace/PS_Current.H"
#include "PHASIC++/Main/Color_Integrator.H"
#include "PHASIC++/Channels/Channel_Elements.H"
#include "ATOOLS/Org/STL_Tools.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Data_Writer.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"

using namespace COMIX;
using namespace PHASIC;
using namespace ATOOLS;

const double s_alphamin(1.0e-12);
const double s_pwmin(1.0e-6), s_pmmin(1.0e-6);

PS_Channel::PS_Channel(const size_t &_nin,const size_t &_nout,
		       ATOOLS::Flavour *_fl,Process_Base *const xs):
  p_xs(xs), m_n(_nin+_nout), m_lid(1), m_rid(2), m_nopt(0),
  p_psid(new PSId_Map()), p_cid(new CId_Map())
{
  nin=_nin;
  nout=_nout;
  m_p.resize(1<<(m_n+1));
  ms = new double[m_n];
  for (size_t i(0);i<m_n;++i) {
    ms[i]=sqr(_fl[i].Mass());
  }
  name="CDBG_Channel";
  Data_Reader read(" ",";","!","=");
  read.SetInputFile("Integration.dat");
  if (!read.ReadFromFile(m_zmode,"CDXS_ZMODE")) m_zmode=0;
  else msg_Info()<<METHOD<<"(): Set zero treatment mode "<<m_zmode<<".\n";
  if (!read.ReadFromFile(m_bmode,"CDXS_BMODE")) m_bmode=1;
  else msg_Info()<<METHOD<<"(): Set boundary mode "<<m_bmode<<".\n";
  if (!read.ReadFromFile(m_omode,"CDXS_OMODE")) m_omode=3;
  else msg_Info()<<METHOD<<"(): Set optimization mode "<<m_omode<<".\n";
  if (!read.ReadFromFile(m_vmode,"CDXS_VMODE")) m_vmode=9;
  else msg_Info()<<METHOD<<"(): Set Vegas mode "<<m_vmode<<".\n";
  if (!read.ReadFromFile(m_tmode,"CDXS_TMODE")) m_tmode=1;
  else msg_Info()<<METHOD<<"(): Set t-channel mode "<<m_tmode<<".\n";
  if (!read.ReadFromFile(m_vsopt,"CDXS_VSOPT")) m_vsopt=0;
  else msg_Info()<<METHOD<<"(): Set Vegas opt start "<<m_vsopt<<".\n";
  if (!read.ReadFromFile(m_nvints,"CDXS_VINTS")) m_nvints=8;
  else msg_Info()<<METHOD<<"(): Set Vegas intervals "<<m_nvints<<".\n";
  if (!read.ReadFromFile(m_texp,"CDXS_TEXP")) m_texp=0.9;
  else msg_Info()<<METHOD<<"(): Set t-channel exp "<<m_texp<<".\n";
  if (!read.ReadFromFile(m_sexp,"CDXS_SEXP")) m_sexp=0.5;
  else msg_Info()<<METHOD<<"(): Set s-channel exp "<<m_sexp<<".\n";
  if (!read.ReadFromFile(m_thexp,"CDXS_THEXP")) m_thexp=1.5;
  else msg_Info()<<METHOD<<"(): Set threshold exp "<<m_thexp<<".\n";
  if (!read.ReadFromFile(m_mfac,"CDXS_MFAC")) m_mfac=1.0;
  else msg_Info()<<METHOD<<"(): Set m_{min} factor "<<m_mfac<<".\n";
  if (!(m_vmode&8)) m_nvints=Max(10,Min(m_nvints,500));
  if (m_vsopt>0) (m_vmode&=~1)|=2;
  m_nr=3*nout-4;
  rannum=m_nr+m_n-2+1;
  rans=new double[rannum];
#ifdef USING__Threading
  int helpi(2);
  if (!read.ReadFromFile(helpi,"COMIX_PS_THREADS")) helpi=2;
  else msg_Tracking()<<METHOD<<"(): Set number of threads "<<helpi<<".\n";
  if (nout<=4) helpi=0;
  else if (nout<=6) helpi=Min(helpi,2);
  if (helpi>0) {
    m_cts.resize(helpi);
    for (size_t i(0);i<m_cts.size();++i) {
      CDBG_PS_TID *tid(new CDBG_PS_TID(this));
      m_cts[i] = tid;
      pthread_cond_init(&tid->m_s_cnd,NULL);
      pthread_cond_init(&tid->m_t_cnd,NULL);
      pthread_mutex_init(&tid->m_s_mtx,NULL);
      pthread_mutex_init(&tid->m_t_mtx,NULL);
      pthread_mutex_lock(&tid->m_s_mtx);
      pthread_mutex_lock(&tid->m_t_mtx);
      tid->m_s=1;
      int tec(0);
      if ((tec=pthread_create(&tid->m_id,NULL,&TGenerateWeight,(void*)tid)))
	THROW(fatal_error,"Cannot create thread "+ToString(i));
    }
  }
  pthread_mutex_init(&m_vgs_mtx,NULL);
  pthread_mutex_init(&m_wvgs_mtx,NULL);
#endif
}

PS_Channel::~PS_Channel()
{
#ifdef USING__Threading
  pthread_mutex_destroy(&m_wvgs_mtx);
  pthread_mutex_destroy(&m_vgs_mtx);
  for (size_t i(0);i<m_cts.size();++i) {
    CDBG_PS_TID *tid(m_cts[i]);
    tid->m_s=0;
    pthread_cond_wait(&tid->m_s_cnd,&tid->m_s_mtx);
    int tec(0);
    if ((tec=pthread_join(tid->m_id,NULL)))
      THROW(fatal_error,"Cannot join thread"+ToString(i));
    pthread_mutex_unlock(&tid->m_t_mtx);
    pthread_mutex_unlock(&tid->m_s_mtx);
    pthread_mutex_destroy(&tid->m_t_mtx);
    pthread_mutex_destroy(&tid->m_s_mtx);
    pthread_cond_destroy(&tid->m_t_cnd);
    pthread_cond_destroy(&tid->m_s_cnd);
    delete tid;
  }
#endif
  for (Vegas_Map::const_iterator vit(m_vmap.begin());
       vit!=m_vmap.end();++vit) delete vit->second;
  delete p_psid;
  delete p_cid;
}

#ifdef USING__Threading
CDBG_PS_TID *PS_Channel::GetTId() const
{
  pthread_t tid(pthread_self());
  for (size_t i(0);i<m_cts.size();++i)
    if (pthread_equal(tid,m_cts[i]->m_id)) return m_cts[i];
  return NULL;
}
#endif

const std::string &PS_Channel::GetPSId(const size_t &id)
{
  PSId_Map *psid(p_psid);
#ifdef USING__Threading
  CDBG_PS_TID *tid(GetTId());
  if (tid) psid=tid->p_psid;
#endif
  PSId_Map::const_iterator iit(psid->find(id));
  if (iit!=psid->end()) return iit->second;
  (*psid)[id]=PSId(id);
  return (*psid)[id];
}

const std::vector<int> &PS_Channel::GetCId(const size_t &id)
{
  CId_Map *cid(p_cid);
#ifdef USING__Threading
  CDBG_PS_TID *tid(GetTId());
  if (tid) cid=tid->p_cid;
#endif
  CId_Map::const_iterator iit(cid->find(id));
  if (iit!=cid->end()) return iit->second;
  (*cid)[id]=ID(id);
  return (*cid)[id];
}

size_t PS_Channel::SId(const size_t &id) const
{
  return (id&3)==3?(1<<m_n)-1-id:id;
}

Vegas *PS_Channel::GetVegas(const std::string &tag,int ni)
{
  Vegas_Map::iterator vit(m_vmap.find(tag));
  if (vit!=m_vmap.end()) return vit->second;
  bool ibi(ni>0);
  if (!ibi) ni=m_nvints;
  Vegas *vegas(new Vegas(1,ni,"CDBG_"+tag,0));
  m_vmap[tag] = vegas;
  if (ibi) vegas->InitBinInfo();
#ifndef CHECK_POINT
  vegas->SetCheckMode(0);
#endif
  if (!ibi && (m_vmode&8)) vegas->SetAutoRefine();
  if (!(m_vmode&4)) vegas->SetOutputMode(0);
  if (m_vmap.size()==1)
    msg_Tracking()<<"  Init internal Vegas map ( "
		  <<m_nvints<<" bins )."<<std::endl;
  if (m_vmode&4)
    msg_Tracking()<<"  Init Vegas "<<std::setw(3)<<std::right
		  <<m_vmap.size()<<" ( "<<ni<<" bins ) '"
		  <<std::setw(35)<<std::left<<tag<<std::right<<"'\n";
  return vegas;
}

PHASIC::Vegas *PS_Channel::GetPVegas
(const Current_Base *cur,const size_t &id)
{
  if (cur!=NULL) {
#ifdef USING__Threading
    pthread_mutex_lock(&m_vgs_mtx);
#endif
    Vegas *vgs(NULL);
    CVegas_Map::const_iterator vit(m_pcmap.find(cur));
    if (vit!=m_pcmap.end()) vgs=vit->second;
    else vgs=m_pcmap[cur]=GetVegas("P_"+cur->PSInfo());
#ifdef USING__Threading
    pthread_mutex_unlock(&m_vgs_mtx);
#endif
    return vgs;
  }
#ifdef USING__Threading
  pthread_mutex_lock(&m_vgs_mtx);
#endif
  Vegas *vgs(NULL);
  IVegas_Map::const_iterator vit(m_pimap.find(id));
  if (vit!=m_pimap.end()) vgs=vit->second;
  else vgs=m_pimap[id]=GetVegas("P_"+GetPSId(id));
#ifdef USING__Threading
  pthread_mutex_unlock(&m_vgs_mtx);
#endif
  return vgs;
}

PHASIC::Vegas *PS_Channel::GetSVegas
(const size_t &id,const Current_Base *cur)
{
#ifdef USING__Threading
  pthread_mutex_lock(&m_vgs_mtx);
#endif
  Vegas *vgs(NULL);
  ICVegas_Map::const_iterator vit(m_sicmap.find(id));
  if (vit!=m_sicmap.end()) {
    CVegas_Map::const_iterator it(vit->second.find(cur));
    if (it!=vit->second.end()) vgs=it->second;
  }
  if (vgs==NULL) vgs=m_sicmap[id][cur]=
    GetVegas("S_"+GetPSId(id)+"_"+cur->PSInfo());
#ifdef USING__Threading
  pthread_mutex_unlock(&m_vgs_mtx);
#endif
  return vgs;
}

PHASIC::Vegas *PS_Channel::GetTVegas
(const size_t &id,const Current_Base *cur)
{
#ifdef USING__Threading
  pthread_mutex_lock(&m_vgs_mtx);
#endif
  Vegas *vgs(NULL);
  ICVegas_Map::const_iterator vit(m_ticmap.find(id));
  if (vit!=m_ticmap.end()) {
    CVegas_Map::const_iterator it(vit->second.find(cur));
    if (it!=vit->second.end()) vgs=it->second;
  }
  if (vgs==NULL) vgs=m_ticmap[id][cur]=
    GetVegas("T_"+GetPSId(id)+"_"+cur->PSInfo());
#ifdef USING__Threading
  pthread_mutex_unlock(&m_vgs_mtx);
#endif
  return vgs;
}

bool PS_Channel::Zero(Vertex_Base *const vtx) const
{
  if (m_czmode&1) return vtx->Zero();
  return false;
}

double PS_Channel::SCut(const size_t &id)
{
  if (id&3) return p_cuts->Getscut(GetPSId((1<<m_n)-1-id));
  return p_cuts->Getscut(GetPSId(id));
}

double PS_Channel::PropMomenta(const Current_Base *cur,const size_t &id,
			       const double &smin,const double &smax,
			       const double *rn)
{
  const double *cr(rn);
  if (cur!=NULL && cur->OnShell())
    return sqr(cur->Flav().Mass());
  if (m_vmode&1) {
    m_vgs.push_back(GetPVegas(cur,id));
    cr=m_vgs.back()->GeneratePoint(rn);
    m_rns.push_back(cr[0]);
#ifdef DEBUG__BG
    msg_Debugging()<<"    generate point "<<m_vgs.back()->Name()<<"\n";
#endif
  }
  if (cur!=NULL) {
    if (cur->Width()>s_pwmin)
      return CE.MassivePropMomenta(cur->Mass(),cur->Width(),1,smin,smax,*cr);
    if (cur->Mass()>s_pmmin) 
      return CE.ThresholdMomenta(m_thexp,m_mfac*cur->Mass(),smin,smax,*cr);
    return CE.MasslessPropMomenta(m_sexp,smin,smax,*cr);
  }
  return CE.MasslessPropMomenta(m_sexp,smin,smax,*cr);
}

double PS_Channel::PropWeight(const Current_Base *cur,const size_t &id,
			      const double &smin,const double &smax,
			      const double &s)
{
  double wgt(1.0), rn;
  if (cur!=NULL) {
    if (cur->OnShell()) return (cur->Mass()*cur->Width())/M_PI;
    if (cur->Width()>s_pwmin) 
      wgt=CE.MassivePropWeight(cur->Mass(),cur->Width(),1,smin,smax,s,rn);
    else if (cur->Mass()>s_pmmin) 
      wgt=CE.ThresholdWeight(m_thexp,m_mfac*cur->Mass(),smin,smax,s,rn);
    else wgt=CE.MasslessPropWeight(m_sexp,smin,smax,s,rn);
  }
  else wgt=CE.MasslessPropWeight(m_sexp,smin,smax,s,rn);
  if (m_vmode&3) {
    Vegas *cvgs(GetPVegas(cur,id));
#ifdef USING__Threading
    pthread_mutex_lock(&m_wvgs_mtx);
#endif
    m_wvgs.push_back(cvgs);
    m_wrns.push_back(rn);
#ifdef USING__Threading
    pthread_mutex_unlock(&m_wvgs_mtx);
#endif
    wgt/=cvgs->GenerateWeight(&rn);
#ifdef DEBUG__BG
    msg_Debugging()<<"    generate weight "<<m_wvgs.back()->Name()<<"\n";
#endif
  }
  return wgt;
}

void PS_Channel::TChannelBounds
(const size_t &aid,const size_t &lid,double &ctmin,double &ctmax,
 const ATOOLS::Vec4D &pa,const ATOOLS::Vec4D &pb,
 const double &s1,const double &s2)
{
  if (m_bmode==0) return;
  const Int_Vector &aidi(GetCId(aid));
  if (aidi.front()==aidi.back()) {
    const Int_Vector &aidj(GetCId(lid));
    if (aidj.front()==aidj.back())
      SingleTChannelBounds(aidi.front(),aidj.front(),
			   ctmin,ctmax,pa,pb,s1,s2,0);
    const Int_Vector &aidk(GetCId((1<<m_n)-1-m_rid-aid-lid));
    if (aidk.front()==aidk.back())
      SingleTChannelBounds(GetCId(m_rid).front(),aidk.front(),
			   ctmin,ctmax,pb,pa,s2,s1,1);
  }
#ifdef DEBUG__BG
  msg_Debugging()<<"    ctmin = "<<ctmin<<", ctmax = "<<ctmax<<"\n";
#endif    
}

void PS_Channel::SingleTChannelBounds
(const size_t &a,const size_t &j,double &rctmin,double &rctmax,
 const ATOOLS::Vec4D &pa,const ATOOLS::Vec4D &pb,
 const double &s1,const double &s2,const int mode)
{
  double ctmin=p_cuts->cosmin[a][j];
  double ctmax=p_cuts->cosmax[a][j];
#ifdef DEBUG__BG
  msg_Debugging()<<"    set t_{"<<a<<","<<j<<"} ctmin = "
		 <<ctmin<<", ctmax = "<<ctmax<<" ("<<mode<<")\n";
#endif    
  double s12((pa+pb).Abs2()), Q12(sqrt(s12));
  double E1((s12+s1-s2)/(2.0*Q12)), pcm2(E1*E1-s1);
  double tmax(p_cuts->scut[a][j]);
  if (tmax<0.0) {
    double sa(pa.Abs2()), Ea((s12+sa-pb.Abs2())/(2.0*Q12));
    double tctmax(E1*Ea+(tmax-s1-sa)/2.0);
    ctmax=Min(tctmax/sqrt(pcm2*(Ea*Ea-sa)),ctmax);
  }
  double pt2(sqr(p_cuts->etmin[j])-s1);
  double ct(sqrt(Max(0.0,1.0-pt2/pcm2)));
  ctmin=Max(ctmin,-ct);
  ctmax=Min(ctmax,ct);
#ifdef DEBUG__BG
  msg_Debugging()<<"    reset t_{"<<a<<","<<j<<"} ctmin = "
		 <<ctmin<<", ctmax = "<<ctmax<<" ("<<mode<<")\n";
#endif    
  if (ctmin>=ctmax) {
    ctmin=p_cuts->cosmin[a][j];
    ctmax=p_cuts->cosmax[a][j];
  }
  rctmin=Max(rctmin,ctmin);
  rctmax=Min(rctmax,ctmax);
}

void PS_Channel::TChannelMomenta
(Current_Base *cur,const size_t &id,const size_t &aid,
 const Vec4D &pa,const Vec4D &pb,Vec4D &p1,Vec4D &p2,
 const double &s1,const double &s2,const double *rns)
{
  const double *cr(rns);
  if (m_vmode&1) {
    m_vgs.push_back(GetTVegas(id,cur));
    cr=m_vgs.back()->GeneratePoint(rns);
    m_rns.push_back(cr[0]);
#ifdef DEBUG__BG
    msg_Debugging()<<"    generate point "<<m_vgs.back()->Name()<<"\n";
#endif
  }
  double ctmin(-1.0), ctmax(1.0);
  TChannelBounds(aid,id,ctmin,ctmax,pa,pb,s1,s2);
  CE.TChannelMomenta(pa,pb,p1,p2,s1,s2,cur->Mass(),
		     m_texp,ctmax,ctmin,1.0,0,cr[0],rns[1]);
}

double PS_Channel::TChannelWeight
(Current_Base *cur,const size_t &id,const size_t &aid,
 const Vec4D &pa,const Vec4D &pb,Vec4D &p1,Vec4D &p2)
{
  double ctmin(-1.0), ctmax(1.0), rns[2];
  TChannelBounds(aid,id,ctmin,ctmax,pa,pb,p1.Abs2(),p2.Abs2());
  double wgt(CE.TChannelWeight(pa,pb,p1,p2,cur->Mass(),
			       m_texp,ctmax,ctmin,1.0,0,rns[0],rns[1]));
  if (m_vmode&3) {
    Vegas *cvgs(GetTVegas(id,cur));
#ifdef USING__Threading
    pthread_mutex_lock(&m_wvgs_mtx);
#endif
    m_wvgs.push_back(cvgs);
    m_wrns.push_back(rns[0]);
#ifdef USING__Threading
    pthread_mutex_unlock(&m_wvgs_mtx);
#endif
    wgt/=cvgs->GenerateWeight(rns);
#ifdef DEBUG__BG
    msg_Debugging()<<"    generate weight "<<m_wvgs.back()->Name()<<"\n";
#endif
  }
  return wgt;
}

void PS_Channel::SChannelBounds
(const size_t &id,const size_t &lid,double &ctmin,double &ctmax)
{
  if (m_bmode==0) return;
  const Int_Vector &aid(GetCId((id&lid)==lid?id:(1<<m_n)-1-id));
  if (aid.size()==2) {
    ctmin=p_cuts->cosmin[aid.front()][aid.back()];
    ctmax=p_cuts->cosmax[aid.front()][aid.back()];
#ifdef DEBUG__BG
    msg_Debugging()<<"    set s_{"<<aid.front()<<","<<aid.back()
		   <<"} ctmin = "<<ctmin<<", ctmax = "<<ctmax<<"\n";
#endif    
  }
}

void PS_Channel::SChannelMomenta
(Current_Base *cur,const size_t &id,const Vec4D &pa,Vec4D &p1,Vec4D &p2,
 const double &s1,const double &s2,const double *rns)
{
  const double *cr(rns);
  if (m_vmode&1) {
    m_vgs.push_back(GetSVegas(id,cur));
    cr=m_vgs.back()->GeneratePoint(rns);
    m_rns.push_back(cr[0]);
#ifdef DEBUG__BG
    msg_Debugging()<<"    generate point "<<m_vgs.back()->Name()<<"\n";
#endif
  }
  double ctmin(-1.0), ctmax(1.0);
  SChannelBounds(cur->CId(),id,ctmin,ctmax);  
  CE.Isotropic2Momenta(pa,s1,s2,p1,p2,cr[0],rns[1],ctmin,ctmax);
}

double PS_Channel::SChannelWeight
(Current_Base *cur,const size_t &id,Vec4D &p1,Vec4D &p2)
{
  double ctmin(-1.0), ctmax(1.0), rns[2];
  SChannelBounds(cur->CId(),id,ctmin,ctmax);
  double wgt(CE.Isotropic2Weight(p1,p2,rns[0],rns[1],ctmin,ctmax));
  if (m_vmode&3) {
    Vegas *cvgs(GetSVegas(id,cur));
#ifdef USING__Threading
    pthread_mutex_lock(&m_wvgs_mtx);
#endif
    m_wvgs.push_back(cvgs);
    m_wrns.push_back(rns[0]);
#ifdef USING__Threading
    pthread_mutex_unlock(&m_wvgs_mtx);
#endif
    wgt/=cvgs->GenerateWeight(rns);
#ifdef DEBUG__BG
    msg_Debugging()<<"    generate weight "<<m_wvgs.back()->Name()<<"\n";
#endif
  }
  return wgt;
}

bool PS_Channel::GeneratePoint
(Current_Base *const ja,Current_Base *const jb,
 Current_Base *const jc,size_t &nr)
{
  size_t aid(ja->CId()), bid(jb->CId()), cid(jc->CId());
  if (((cid&m_lid)==m_lid)^((cid&m_rid)==m_rid)) {
    size_t pid(aid-(m_rid+bid));
    double se(SCut(bid)), sp(SCut(pid));
    double rtsmax((m_p[aid]+m_p[m_rid]).Mass());
    if (CIdCount(bid)>1) {
      double smin(se), smax(sqr(rtsmax-sqrt(sp)));
      se=PropMomenta(jb,bid,smin,smax,&rans[nr++]);
    }
    if (CIdCount(pid)>1) {
      double smin(sp), smax(sqr(rtsmax-sqrt(se)));
      sp=PropMomenta(((PS_Current*)jc)->SCC(),pid,
		     smin,smax,&rans[nr++]);
    }
    TChannelMomenta(jc,bid,(1<<m_n)-1-aid,m_p[aid],m_p[m_rid],
		    m_p[bid],m_p[pid],se,sp,&rans[nr]);
    nr+=2;
    m_p[cid]=m_p[aid]-m_p[bid];
#ifdef DEBUG__BG
    msg_Debugging()<<"  t "<<nr<<": ("<<PSId(ja->CId())
		   <<","<<PSId(m_rid)<<")-"<<PSId(jc->CId())
		   <<"->("<<PSId(jb->CId())<<","<<PSId(pid)
		   <<") m_"<<PSId(bid)<<" = "<<sqrt(se)
		   <<", m_"<<PSId(pid)<<" = "<<sqrt(sp)
		   <<" -> "<<PSId(cid)<<"\n";
#endif
  }
  else {
    size_t lid(SId(aid)), rid(SId(bid));
    double rts(m_p[cid].Mass()), sl(SCut(lid)), sr(SCut(rid));
    if (CIdCount(lid)>1) {
      double smin(sl), smax(sqr(rts-sqrt(sr)));
      sl=PropMomenta(ja,lid,smin,smax,&rans[nr++]);
    }
    if (CIdCount(rid)>1) {
      double smin(sr), smax(sqr(rts-sqrt(sl)));
      sr=PropMomenta(jb,rid,smin,smax,&rans[nr++]);
    }
    SChannelMomenta(jc,lid,m_p[cid],m_p[aid],m_p[bid],sl,sr,&rans[nr]);
    nr+=2;
    m_p[(1<<m_n)-1-aid]=m_p[aid];
    m_p[(1<<m_n)-1-bid]=m_p[bid];
#ifdef DEBUG__BG
    msg_Debugging()<<"  s "<<nr<<": ("<<PSId(cid)
		   <<")->("<<PSId(aid)<<","<<PSId(bid)
		   <<") m_"<<PSId(cid)<<" = "<<rts
		   <<", m_"<<PSId(lid)<<" = "<<sqrt(sl)
		   <<", m_"<<PSId(rid)<<" = "<<sqrt(sr)<<"\n";
#endif
  }
  return true;
}

bool PS_Channel::GeneratePoint
(const size_t &id,size_t &nr,Vertex_Vector &v)
{
  for (size_t i(0);i<v.size() && nr<m_nr;++i) {
    if (v[i]==NULL) continue;
    size_t cid(v[i]->JC()->CId());
    size_t aid(v[i]->JA()->CId()), bid(v[i]->JB()->CId());
    if (aid==id || bid==id || cid==id || (1<<m_n)-1-cid==id) {
      Current_Base *ja(v[i]->JA()), *jb(v[i]->JB()), *jc(v[i]->JC());
      if (aid==id) { 
	std::swap<size_t>(aid,cid);
	std::swap<Current_Base*>(ja,jc);
      }
      else if (bid==id) {
	std::swap<size_t>(bid,cid);
	std::swap<Current_Base*>(jb,jc);
      }
      if (!GeneratePoint(ja,jb,jc,nr)) return false;
      v[i]=NULL;
      if (CIdCount(SId(aid))>1) GeneratePoint(aid,nr,v);
      if (CIdCount(SId(bid))>1) GeneratePoint(bid,nr,v);
      break;
    }
  }
  return true;
}

bool PS_Channel::GeneratePoint(Vertex_Vector v)
{
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"(): {\n";
#endif
  size_t nr(0), lid((1<<m_n)-2), rid(m_rid);
  for (size_t n(2);n<=m_n-2;++n) {
    for (size_t i(0);i<v.size() && nr<m_nr;++i) {
      if (v[i]==NULL) continue;
#ifdef DEBUG__BG
      msg_Debugging()<<" "<<lid<<" "<<*v[i]<<"\n";
#endif
      size_t cid(v[i]->JC()->CId());
      size_t aid(v[i]->JA()->CId()), bid(v[i]->JB()->CId());
      if (aid==lid || bid==lid || cid==lid) {
	Current_Base *ja(v[i]->JA()), *jb(v[i]->JB()), *jc(v[i]->JC());
	if (bid==lid) {
	  std::swap<size_t>(aid,bid);
	  std::swap<Current_Base*>(ja,jb);
	}
	else if (cid==lid) {
	  std::swap<size_t>(aid,cid);
	  std::swap<Current_Base*>(ja,jc);
	}
	if ((cid&(lid|rid))==(lid|rid) || (aid&rid && bid&rid)) {
	  std::swap<size_t>(bid,cid);
	  std::swap<Current_Base*>(jb,jc);
	}
	if (cid==rid) {
	  v[i]=NULL;
	  if (bid!=3) m_p[bid]=m_p[aid-cid];
	  if (CIdCount(bid)>1) GeneratePoint(bid,nr,v);
	  break;
	}
	if (!GeneratePoint(ja,jb,jc,nr)) return false;
	v[i]=NULL;
	if (CIdCount(bid)>1) GeneratePoint(bid,nr,v);
	lid=cid;
      }
    }
  }
  if (nr!=m_nr) THROW(fatal_error,"Internal error");
#ifdef DEBUG__BG
  msg_Debugging()<<"} -> "<<nr<<"\n";
#endif
  return true;
}

bool PS_Channel::GenerateChannel
(Current_Base *const cur,Vertex_Vector &v)
{
  if (cur->NIn()==0) return true;
#ifdef DEBUG__BG
  msg_Indent();
  msg_Debugging()<<METHOD<<"("<<PSId(cur->CId())
		 <<","<<cur->J<PS_Info>().size()
		 <<"): n = "<<v.size()<<" {\n";
#endif
  double sum(0.0);
  Double_Vector psum;
  Vertex_Vector vtcs;
  for (size_t i(0);i<cur->In().size();++i)
    if (!Zero(cur->In()[i])) {
      vtcs.push_back(cur->In()[i]);
      psum.push_back(sum+=((PS_Vertex*)vtcs.back())->Alpha());
    }
  Vertex_Base *vtx(NULL);
  for (size_t i(0);i<psum.size();++i)
    if (psum[i]>=rans[m_nr+v.size()]*sum) {
      vtx=vtcs[i];
      break;
    }
  if (vtx==NULL) {
    if (m_czmode==0) THROW(fatal_error,"No vertex in z mode 0");
#ifdef DEBUG__BG
    msg_Debugging()<<METHOD<<"(): No vertex, switching z mode."<<std::endl;
#endif
    m_czmode=0;
    v.clear();
    return GenerateChannel((*p_cur)[m_n-1].back(),v);
  }
  v.push_back(vtx);
#ifdef DEBUG__BG
  msg_Debugging()<<"  "<<*cur<<" <- ("<<vtcs.size()<<") "
		 <<std::flush<<*vtx<<"\n";
#endif
  if (v.size()<m_n-2 && !GenerateChannel(vtx->JA(),v)) return false;
  if (v.size()<m_n-2 && !GenerateChannel(vtx->JB(),v)) return false;
#ifdef DEBUG__BG
  msg_Debugging()<<"}\n";
#endif
  return true;
}

bool PS_Channel::GenerateChannel(Vertex_Vector &v)
{
  m_czmode=m_zmode;
  if (!GenerateChannel((*p_cur)[m_n-1].back(),v)) return false;
  if (v.size()!=m_n-2) THROW(fatal_error,"Internal error");
#ifdef DEBUG__BG
  for (size_t i(0);i<v.size();++i)
    msg_Debugging()<<*v[i]<<"\n";
#endif
  return true;
}

bool PS_Channel::GenerateChannels()
{
  PHASIC::Process_Base *cur(p_xs->Process());
  p_gen=cur->Get<Process_Base>()->PSGenerator();
  if (p_gen==NULL) 
    THROW(fatal_error,"No phasespace generator for "+cur->Name());
  p_gen->SetZMode(m_zmode);
  if (!p_gen->Evaluate()) return false;
  p_cur = (Current_Matrix*)(&p_gen->Graphs());
  return true;
}

void PS_Channel::GeneratePoint
(ATOOLS::Vec4D *p,PHASIC::Cut_Data *cuts,double *rn) 
{
  if (!GenerateChannels()) return;
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"(): {\n";
#endif
  p_cuts=cuts;
  p_gen->SetPrefMasses(p_cuts);
  m_p[(1<<m_n)-1-3]=m_p[3]=
    (m_p[(1<<m_n)-1-1]=m_p[1]=p[0])+
    (m_p[(1<<m_n)-1-2]=m_p[2]=p[1]);
  for (int i(0);i<rannum;++i) rans[i]=rn[i];
  Vertex_Vector v;
  if (!GenerateChannel(v)) return;
  m_vgs.clear();
  m_rns.clear();
  m_wvgs.clear();
  m_wrns.clear();
  if (!GeneratePoint(v)) return;
  Vec4D sum(-p[0]-p[1]);
  for (size_t i(2);i<m_n;++i) sum+=p[i]=m_p[1<<i];
  if (!IsEqual(sum,Vec4D(),sqrt(Accu()))) msg_Error()
    <<METHOD<<"(): Four momentum not conserved. Diff "<<sum<<std::endl;
#ifdef DEBUG__BG
  msg_Debugging()<<"}\n";
#endif
}

double PS_Channel::GenerateWeight
(Current_Base *const ja,Current_Base *const jb,
 Current_Base *const jc,size_t &nr)
{
  double wgt(1.0);
  size_t aid(ja->CId()), bid(jb->CId()), cid(jc->CId());
  if (((cid&m_lid)==m_lid)^((cid&m_rid)==m_rid)) {
    size_t pid(aid-(m_rid+bid));
    aid=(1<<m_n)-1-aid;
    m_p[pid]=-m_p[aid]-m_p[m_rid]-m_p[bid];
    double se(SCut(bid)), sp(SCut(pid));
    double rtsmax((m_p[aid]+m_p[m_rid]).Mass());
    if (CIdCount(bid)>1) {
      double smin(se), smax(sqr(rtsmax-sqrt(sp)));
      wgt*=PropWeight(jb,bid,smin,smax,se=m_p[bid].Abs2());
    }
    if (CIdCount(pid)>1) {
      double smin(sp), smax(sqr(rtsmax-sqrt(se)));
      wgt*=PropWeight(((PS_Current*)jc)->SCC(),pid,
		      smin,smax,sp=m_p[pid].Abs2());
    }
    wgt*=TChannelWeight(jc,bid,aid,-m_p[aid],-m_p[m_rid],
			m_p[bid],m_p[pid]);
    nr+=2;
#ifdef DEBUG__BG
    msg_Debugging()<<"    t "<<nr<<": ("<<PSId(ja->CId())
		   <<","<<PSId(m_rid)<<")-"<<PSId(jc->CId())
		   <<"->("<<PSId(jb->CId())<<","<<PSId(pid)
		   <<") m_"<<PSId(bid)<<" = "<<sqrt(se)
		   <<", m_"<<PSId(pid)<<" = "<<sqrt(sp)
		   <<" -> m_"<<PSId(aid^m_rid)<<" = "
		   <<(m_p[aid]+m_p[m_rid]).Mass()<<"\n";
#endif
  }
  else {
    size_t lid(SId(aid)), rid(SId(bid));
    double rts(m_p[cid].Mass()), sl(SCut(lid)), sr(SCut(rid));
    if (CIdCount(lid)>1) {
      double smin(sl), smax(sqr(rts-sqrt(sr)));
      wgt*=PropWeight(ja,lid,smin,smax,sl=m_p[lid].Abs2());
    }
    if (CIdCount(rid)>1) {
      double smin(sr), smax(sqr(rts-sqrt(sl)));
      wgt*=PropWeight(jb,rid,smin,smax,sr=m_p[rid].Abs2());
    }
    wgt*=SChannelWeight(jc,lid,m_p[lid],m_p[rid]);
    nr+=2;
#ifdef DEBUG__BG
    msg_Debugging()<<"    s "<<nr<<": ("<<PSId(cid)
		   <<")->("<<PSId(aid)<<","<<PSId(bid)
		   <<") m_"<<PSId(SId(cid))<<" = "<<rts
		   <<", m_"<<PSId(lid)<<" = "<<sqrt(sl)
		   <<", m_"<<PSId(rid)<<" = "<<sqrt(sr)<<"\n";
#endif
  }
  return wgt;
}

bool PS_Channel::GenerateWeight(PS_Current *const cur)
{
  double wgt(0.0), asum(0.0);
#ifdef DEBUG__BG
  msg_Debugging()<<"  J_"<<PSId(cur->CId())<<" ("
		 <<((Current_Base*)cur)->J<PS_Info>().size()
		 <<"/"<<cur->Zero()<<"): {\n";
#endif
  for (size_t i(0);i<cur->In().size();++i) {
    PS_Vertex *v((PS_Vertex *)cur->In()[i]);
    if (!Zero(v) && v->Alpha()>0.0) {
      size_t nr(0);
      Current_Base *ja(v->JA()), *jb(v->JB()), *jc(cur);
      size_t aid(ja->CId()), bid(jb->CId()), cid(jc->CId());
      double cw(ja->J<PS_Info>().front()[0]*jb->J<PS_Info>().front()[0]);
      if ((((aid&m_lid)==m_lid)^((aid&m_rid)==m_rid)) ||
	  (((bid&m_lid)==m_lid)^((bid&m_rid)==m_rid))) {
	if ((bid&m_lid)==m_lid) {
	  std::swap<size_t>(aid,bid);
	  std::swap<Current_Base*>(ja,jb);
	}
	else if ((cid&m_lid)!=m_lid) {
	  std::swap<size_t>(aid,cid);
	  std::swap<Current_Base*>(ja,jc);
	}
	if ((cid&(m_lid|m_rid))==(m_lid|m_rid) || 
	    (aid&m_rid && bid&m_rid)) {
	  std::swap<size_t>(bid,cid);
	  std::swap<Current_Base*>(jb,jc);
	}
	if (cid==m_rid) {
#ifdef DEBUG__BG
	  msg_Debugging()<<"    kill "<<PSId(aid)<<" "<<PSId(bid)
			 <<" "<<PSId(cid)<<"\n";
#endif
	  v->SetWeight(cw);
	  wgt+=v->Alpha()*v->Weight();
	  asum+=v->Alpha();
#ifdef DEBUG__BG
	  msg_Debugging()<<"    w = "<<1.0/v->Weight()
			 <<", a = "<<((PS_Vertex*)v)->Alpha()<<"\n";
#endif
	  continue;
	}
      }
      else {
	if (aid&(m_lid+m_rid) && CIdCount(aid)<CIdCount(cid)) { 
	  std::swap<size_t>(aid,cid);
	  std::swap<Current_Base*>(ja,jc);
	}
	else if (bid&(m_lid+m_rid) && CIdCount(bid)<CIdCount(cid)) { 
	  std::swap<size_t>(bid,cid);
	  std::swap<Current_Base*>(jb,jc);
	}
      }
      v->SetWeight(GenerateWeight(ja,jb,jc,nr)*cw);
      wgt+=v->Alpha()*v->Weight();
      asum+=v->Alpha();
#ifdef DEBUG__BG
      msg_Debugging()<<"    w = "<<1.0/v->JA()->J<PS_Info>().front()[0]
		     <<" * "<<1.0/v->JB()->J<PS_Info>().front()[0]
		     <<" * "<<cw/v->Weight()<<" = "<<1.0/v->Weight()
		     <<", a = "<<v->Alpha()<<"\n";
#endif
    }
  }
  wgt/=asum;
  if (m_omode>0)
    for (size_t i(0);i<cur->In().size();++i) {
      PS_Vertex *v((PS_Vertex*)cur->In()[i]);
      if (!Zero(v) && v->Alpha()>0.0) {
#ifdef DEBUG__BG
	msg_Debugging()<<"    V_{"<<PSId(v->JA()->CId())
		       <<","<<PSId(v->JB()->CId())
		       <<"}: set w = "<<v->Weight()/wgt<<"\n";
#endif
	if (wgt>0.0) v->SetWeight(v->Weight()/wgt);
	else v->SetWeight(0.0);
      }
    }
#ifdef DEBUG__BG
  msg_Debugging()<<"  } -> w = "<<1.0/(wgt*asum)
		 <<" * "<<asum<<" = "<<1.0/wgt<<"\n";
#endif
  cur->ResetJ();
  cur->AddJ(PS_Info(0,0,0,wgt));
  return true;
}

#ifdef USING__Threading
void *PS_Channel::TGenerateWeight(void *arg)
{
  CDBG_PS_TID *tid((CDBG_PS_TID*)arg);
  while (true) {
    // wait for channel to signal
    pthread_mutex_lock(&tid->m_s_mtx);
    pthread_mutex_unlock(&tid->m_s_mtx);
    pthread_cond_signal(&tid->m_s_cnd);
    if (tid->m_s==0) return NULL;
    // worker routine
    for (tid->m_i=tid->m_b;tid->m_i<tid->m_e;++tid->m_i) 
      if (!tid->p_psc->GenerateWeight
	  ((PS_Current*)(*tid->p_psc->p_cur)
	   [tid->m_n][tid->m_i]))
	THROW(fatal_error,"Generate weight failed");
    // signal channel to continue
    pthread_cond_wait(&tid->m_t_cnd,&tid->m_t_mtx);
  }
  return NULL;
}
#endif

bool PS_Channel::GenerateWeight()
{
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"(): {\n";
#endif
  for (size_t n(2);n<m_n;++n) {
#ifdef USING__Threading
    if (m_cts.empty()) {
      for (size_t i(0);i<(*p_cur)[n].size();++i) 
	if (!GenerateWeight((PS_Current*)(*p_cur)[n][i])) return 0.0;
    }
    else {
      // start calculator threads
      size_t d((*p_cur)[n].size()/m_cts.size());
      if ((*p_cur)[n].size()%m_cts.size()>0) ++d;
      for (size_t j(0), i(0);j<m_cts.size()&&i<(*p_cur)[n].size();++j) {
	CDBG_PS_TID *tid(m_cts[j]);
	tid->m_n=n;
	tid->m_b=i;
	tid->m_e=Min(i+=d,(*p_cur)[n].size());
	pthread_cond_wait(&tid->m_s_cnd,&tid->m_s_mtx);
      }
      // suspend calculator threads
      for (size_t j(0), i(0);j<m_cts.size()&&i<(*p_cur)[n].size();++j) {
	i+=d;
	CDBG_PS_TID *tid(m_cts[j]);
	pthread_mutex_lock(&tid->m_t_mtx);
	pthread_mutex_unlock(&tid->m_t_mtx);
	pthread_cond_signal(&tid->m_t_cnd);
      }
    }
#else
    for (size_t i(0);i<(*p_cur)[n].size();++i) 
      if (!GenerateWeight((PS_Current*)(*p_cur)[n][i])) return 0.0;
#endif
  }
  weight=1.0/(*p_cur)[m_n-1].back()->J<PS_Info>().front()[0]/
    pow(2.0*M_PI,3.0*nout-4.0);
#ifdef DEBUG__BG
  msg_Debugging()<<"} -> "<<weight<<"\n";
#endif
  return true;
}

void PS_Channel::GenerateWeight(ATOOLS::Vec4D *p,PHASIC::Cut_Data *cuts) 
{
  p_cuts=cuts;
  for (size_t i(0);i<m_n;++i) {
    m_p[1<<i]=i<2?-p[i]:p[i];
#ifdef DEBUG__BG
    msg_Debugging()<<"  p_"<<i<<" = "<<m_p[1<<i]<<"\n";
#endif
  }
  for (size_t n(2);n<p_cur->size();++n)
    for (size_t i(0);i<(*p_cur)[n].size();++i) {
      Current_Base *cur((*p_cur)[n][i]);
      if (cur->In().empty()) THROW(fatal_error,"Internal error");
      m_p[(1<<m_n)-1-cur->CId()]=
	-(m_p[cur->CId()]=m_p[cur->In().front()->JA()->CId()]
	  +m_p[cur->In().front()->JB()->CId()]);
#ifdef DEBUG__BG
	msg_Debugging()<<"  -p_"<<PSId((1<<m_n)-1-cur->CId())
		       <<" = p_"<<PSId(cur->CId())
		       <<" = "<<m_p[cur->CId()]<<"\n";
#endif
    }
  if (!GenerateWeight())
    THROW(fatal_error,"Internal error");
}

void PS_Channel::AddPoint(double value)
{ 
  Single_Channel::AddPoint(value);
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"(): value = "<<value<<" {\n";
#endif
  if (m_omode>0)
    for (size_t n(2);n<m_n;++n) 
      for (size_t i(0);i<(*p_cur)[n].size();++i) {
	for (size_t j(0);j<(*p_cur)[n][i]->In().size();++j) {
	  PS_Vertex *v((PS_Vertex *)(*p_cur)[n][i]->In()[j]);
	  if (!Zero(v)) v->AddPoint(value);
#ifdef DEBUG__BG
	  msg_Debugging()<<"    V_{"<<PSId(v->JA()->CId())
			 <<","<<PSId(v->JB()->CId())<<"}: <w> = "
			 <<v->Mean()<<" +- "<<v->Sigma()
			 <<", max = "<<v->Max()<<"\n";
#endif
	}
      }
  if (m_vmode&3) {
#ifdef CHECK_POINT
    for (int i(m_vgs.size()-1);i>=0;--i) {
#ifdef DEBUG__BG
      msg_Debugging()<<"  add point "<<m_vgs[i]->Name()<<"\n";
#endif
      m_vgs[i]->AddPoint(value,&m_rns[i]);
    }
#endif
    for (int i(m_wvgs.size()-1);i>=0;--i) {
#ifdef DEBUG__BG
      msg_Debugging()<<"  add point "<<m_wvgs[i]->Name()<<"\n";
#endif
      m_wvgs[i]->AddPoint(value,&m_wrns[i]);
    }
  }
#ifdef DEBUG__BG
  msg_Debugging()<<"}\n";
#endif
}

void PS_Channel::Optimize()  
{
  ++m_nopt;
  if (m_omode>0) {
    msg_Tracking()<<METHOD<<"(): mode = "<<m_omode<<" {\n";
    msg_Tracking()<<"  "<<std::string(108,'-')<<"\n";
    for (size_t n(2);n<m_n;++n) {
      for (size_t i(0);i<(*p_cur)[n].size();++i) {
	double csum(0.0), wmean(0.0), nc(0.0);
	for (size_t j(0);j<(*p_cur)[n][i]->In().size();++j) {
	  PS_Vertex *v((PS_Vertex *)(*p_cur)[n][i]->In()[j]);
	  if (v->Alpha()>0.0) {
	    v->SetOldAlpha(v->Alpha());
	    if (m_omode&2) v->SetAlpha(v->Alpha()*sqrt(v->Mean()));
	    csum+=v->Alpha();
	    wmean+=v->Mean();
	    ++nc;
	  }
	}
	csum/=nc;
	wmean/=nc;	
	bool printed(false);
	for (size_t j(0);j<(*p_cur)[n][i]->In().size();++j) {
	  PS_Vertex *v((PS_Vertex *)(*p_cur)[n][i]->In()[j]);
	  if (v->Alpha()>0.0) {
	    v->SetAlpha(v->Alpha()/csum);
	    double ee(Min(1.0,v->Sigma()/v->Mean()));
	    v->SetAlpha(pow(v->OldAlpha(),ee)*pow(v->Alpha(),1.0-ee));
	    if (v->Alpha()<s_alphamin) v->SetAlpha(0.0);
	  }
	  double dev(int((v->Alpha()/v->OldAlpha()-1.0)*10000)/100.0);
	  if (v->Alpha()!=1.0) {
	    printed=true;
	    double re(int(v->Sigma()/v->Mean()*10000)/100.0);
	    msg_Tracking()<<"  V_{"<<std::setw(6)<<PSId(v->JA()->CId())
			  <<","<<std::setw(6)<<PSId(v->JB()->CId())
			  <<"}: w' = "<<std::setw(15)<<std::right
			  <<v->Mean()/wmean<<" +- "<<std::setw(6)<<re
			  <<" %  =>  a = "<<std::setw(15)<<v->OldAlpha()
			  <<" -> "<<std::setw(15)<<v->Alpha();
	    if (v->Alpha()<s_alphamin) {
	      msg_Tracking()<<std::left<<" (      off )\n";
	    }
	    else {
	      msg_Tracking()<<" ( "<<std::setw(6)
			    <<std::right<<dev<<std::left<<" % )\n";
	    }
	  }
	  if (m_omode&1) v->Reset();
	}
	if (printed)
	  msg_Tracking()<<"  "<<std::string(108,'-')<<"\n";
      }
    }
    msg_Tracking()<<"}"<<std::endl;
    if (m_vmode&2 && m_nopt>=m_vsopt) (m_vmode&=~2)|=1;
  }
  if (m_vmode&1) {
    for (Vegas_Map::const_iterator vit(m_vmap.begin());
	 vit!=m_vmap.end();++vit) {
      msg_Indent();
      vit->second->Optimize();
    }
  }
} 

void PS_Channel::EndOptimize()  
{
  m_omode=0;
  if (m_vmode&1) {
    for (Vegas_Map::const_iterator vit(m_vmap.begin());
	 vit!=m_vmap.end();++vit) vit->second->EndOptimize();
  }
} 

void PS_Channel::ISRInfo(int &type,double &mass,double &width) 
{
  type=0; 
  mass=width=0.0;
}

void PS_Channel::ISRInfo
(std::vector<int> &ts,std::vector<double> &ms,std::vector<double> &ws) const
{
  SP(PS_Generator) ps(p_xs->PSGenerator());
  if (ps==NULL) {
    ps=(*p_xs->Process())[0]->Get<Process_Base>()->PSGenerator();
  }
  msg_Debugging()<<METHOD<<"(): Add isr infos {\n";
  bool addth(ps->ThresholdMass()>0.0);
  const Double_Vector &mps(ps->ISRMasses()), &wps(ps->ISRWidths());
  for (size_t i(0);i<mps.size();++i) {
    msg_Debugging()<<"  resonance "<<i<<": "<<mps[i]<<" / "<<wps[i]<<"\n";
    if (IsEqual(mps[i],ps->ThresholdMass(),s_pwmin)) addth=false;
    ts.push_back(1);
    ms.push_back(mps[i]);
    ws.push_back(wps[i]);
  }
  if (addth) {
    msg_Debugging()<<"  threshold  : "<<ps->ThresholdMass()<<"\n";
    ts.push_back(2);
    ms.push_back(ps->ThresholdMass());
    ws.push_back(0.0);
    ts.push_back(2);
    ms.push_back(2.0*ps->ThresholdMass());
    ws.push_back(0.0);
  }
  msg_Debugging()<<"}\n";
}

int PS_Channel::ChNumber()
{ 
  return m_num; 
}

void PS_Channel::SetChNumber(int n) 
{ 
  m_num=n; 
}

std::string PS_Channel::ChID() 
{
  return name;
}

void PS_Channel::WriteOut(std::string pid)
{ 
  {
    Data_Writer writer;
    writer.SetOutFileMode(fom::permanent|fom::nosearch);
    writer.SetOutputPath(pid);
    writer.SetOutputFile("_"+name+"_PS");
    writer.WriteToFile(m_zmode,"m_zmode");
    writer.WriteToFile(m_bmode,"m_bmode");
    writer.WriteToFile(m_omode,"m_omode");
    writer.WriteToFile(m_vmode,"m_vmode");
    writer.WriteToFile(m_tmode,"m_tmode");
    writer.WriteToFile(m_vsopt,"m_vsopt");
    writer.WriteToFile(m_nvints,"m_nvints");
    writer.WriteToFile(m_texp,"m_texp");
    writer.WriteToFile(m_sexp,"m_sexp");
    writer.WriteToFile(m_thexp,"m_thexp");
    writer.WriteToFile(m_mfac,"m_mfac");
    writer.WriteToFile(m_nopt,"m_nopt");
  }
  if (m_vmode>0) {
    std::vector<std::string> vids;
    for (Vegas_Map::const_iterator vit(m_vmap.begin());
	 vit!=m_vmap.end();++vit) {
      msg_Debugging()<<"write out vegas '"<<vit->first<<"'\n";
      vit->second->WriteOut(pid);
      vids.push_back(vit->first);
    }
    Data_Writer writer;
    writer.SetOutFileMode(fom::permanent|fom::nosearch);
    writer.SetOutputPath(pid);
    writer.SetOutputFile("_"+name+"_VI");
    writer.SetVectorType(vtc::vertical);
    writer.VectorToFile(vids);
  }
  p_cur = (Current_Matrix*)
    (&p_xs->Process()->Get<Process_Base>()->PSGenerator()->Graphs());
  std::vector<std::vector<std::string> > pvds;
  for (size_t pc(0), n(2);n<m_n;++n)
    for (size_t i(0);i<(*p_cur)[n].size();++i) {
      pvds.resize(pvds.size()+(*p_cur)[n][i]->In().size(),
		  std::vector<std::string>(7));
      for (size_t j(0);j<(*p_cur)[n][i]->In().size();++j) {
	PS_Vertex *v((PS_Vertex *)(*p_cur)[n][i]->In()[j]);
	pvds[pc][0]=v->VId();
	pvds[pc][1]=ToString(v->Alpha(),12);
	pvds[pc][2]=ToString(v->OldAlpha(),12);
	pvds[pc][3]=ToString(v->N(),12);
	pvds[pc][4]=ToString(v->Sum(),12);
	pvds[pc][5]=ToString(v->Sum2(),12);
	pvds[pc][6]=ToString(v->Max(),12);
	++pc;
      }
    }
  Data_Writer writer;
  writer.SetOutFileMode(fom::permanent|fom::nosearch);
  writer.SetOutputPath(pid);
  writer.SetOutputFile("_"+name+"_PV");
  writer.MatrixToFile(pvds);
}

void PS_Channel::ReadIn(std::string pid)
{
  {
    Data_Reader reader;
    reader.SetInFileMode(fom::permanent|fom::nosearch);
    reader.SetAddCommandLine(false);
    reader.SetInputPath(pid);
    reader.SetInputFile("_"+name+"_PS");
    reader.ReadFromFile(m_zmode,"m_zmode");
    reader.ReadFromFile(m_bmode,"m_bmode");
    reader.ReadFromFile(m_omode,"m_omode");
    reader.ReadFromFile(m_vmode,"m_vmode");
    reader.ReadFromFile(m_tmode,"m_tmode");
    reader.ReadFromFile(m_vsopt,"m_vsopt");
    reader.ReadFromFile(m_nvints,"m_nvints");
    reader.ReadFromFile(m_texp,"m_texp");
    reader.ReadFromFile(m_sexp,"m_sexp");
    reader.ReadFromFile(m_thexp,"m_thexp");
    reader.ReadFromFile(m_mfac,"m_mfac");
    reader.ReadFromFile(m_nopt,"m_nopt");
  }
  p_gen=p_xs->Process()->Get<Process_Base>()->PSGenerator();
  p_gen->SetPrefMasses
    (p_xs->Process()->Integrator()->PSHandler()->Cuts());
  Data_Reader reader;
  reader.SetInFileMode(fom::permanent|fom::nosearch);
  reader.SetAddCommandLine(false);
  reader.SetInputPath(pid);
  reader.SetInputFile("_"+name+"_PV");
  std::vector<std::vector<std::string> > pvds;
  reader.MatrixFromFile(pvds);
  p_cur = (Current_Matrix*)
    (&p_xs->Process()->Get<Process_Base>()->PSGenerator()->Graphs());
  for (size_t pc(0), n(2);n<m_n;++n){
    for (size_t i(0);i<(*p_cur)[n].size();++i)
      for (size_t j(0);j<(*p_cur)[n][i]->In().size();++j) {
	PS_Vertex *v((PS_Vertex *)(*p_cur)[n][i]->In()[j]);
	if (pc>=pvds.size() || pvds[pc][0]!=v->VId()) 
	  THROW(fatal_error,"Corrupted input file");
	v->SetAlpha(ToType<double>(pvds[pc][1]));
	v->SetOldAlpha(ToType<double>(pvds[pc][2]));
	v->SetN(ToType<double>(pvds[pc][3]));
	v->SetSum(ToType<double>(pvds[pc][4]));
	v->SetSum2(ToType<double>(pvds[pc][5]));
	v->SetMax(ToType<double>(pvds[pc][6]));
	++pc;
      }
  }
  if (m_vmode>0) {
    Data_Reader reader;
    reader.SetInFileMode(fom::permanent|fom::nosearch);
    reader.SetAddCommandLine(false);
    reader.SetInputPath(pid);
    reader.SetInputFile("_"+name+"_VI");
    reader.SetVectorType(vtc::vertical);
    std::vector<std::string> vids;
    if (reader.VectorFromFile(vids)) {
      for (size_t i(0);i<vids.size();++i) {
	msg_Debugging()<<"read in vegas '"<<vids[i]<<"'\n";
	Vegas *vegas(NULL);
	if (vids[i].find("C_")!=0) {
	  vegas=GetVegas(vids[i]);
	}
	else {
	  size_t pos(vids[i].rfind('_'));
	  vegas=GetVegas(vids[i],ToType<int>
			 (vids[i].substr(pos+1)));
	}
	vegas->ReadIn(pid);
      }
    }
  }
}
