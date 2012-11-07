#include "COMIX/Phasespace/PS_Generator.H"

#include "COMIX/Phasespace/PS_Current.H"
#include "COMIX/Phasespace/PS_Vertex.H"
#include "COMIX/Phasespace/PS_Channel.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/STL_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Main/Color_Integrator.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/Single_Process.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "COMIX/Main/Process_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Smart_Pointer.H"
#include "ATOOLS/Org/Smart_Pointer.C"

#include <iomanip>

using namespace COMIX;
using namespace ATOOLS;
using namespace PHASIC;

const double s_pwmin(1.0e-6);

namespace ATOOLS { template class SP(PS_Generator); }

PS_Generator::PS_Generator(Process_Base *const xs):
  p_xs(xs), m_n(0), m_zmode(1), m_pmsinit(0),
  m_thmass(0.0), m_chmass(0.0)
{
  Data_Reader read(" ",";","!","=");
  read.SetInputFile("Integration.dat");
  if (!read.ReadFromFile(m_itmin,"CDXS_ITMIN")) m_itmin=5000;
  else msg_Info()<<METHOD<<"(): Set iteration minimum "<<m_itmin<<".\n";
  if (!read.ReadFromFile(m_itmax,"CDXS_ITMAX")) m_itmax=50000;
  else msg_Info()<<METHOD<<"(): Set iteration maximum "<<m_itmax<<".\n";
  if (!read.ReadFromFile(m_ecmode,"CDXS_ECMODE")) m_ecmode=2;
  else msg_Info()<<METHOD<<"(): Set extra channel mode "<<m_ecmode<<".\n";
  if (!read.ReadFromFile(m_chmass,"CDXS_PS_CHTH")) m_chmass=0.01;
  else msg_Info()<<METHOD<<"(): Set channel mass threshold "<<m_chmass<<".\n";
  m_chmass*=rpa->gen.Ecms();
  p_xs->ConstructPSVertices(this);
  AddSTCC();
#ifdef USING__Threading
  int helpi(0);
  if (!read.ReadFromFile(helpi,"COMIX_PG_THREADS")) helpi=0;
  else msg_Tracking()<<METHOD<<"(): Set number of threads "<<helpi<<".\n";
  if (helpi>0) {
    m_cts.resize(helpi);
    for (size_t i(0);i<m_cts.size();++i) {
      CDBG_PG_TID *tid(new CDBG_PG_TID(this));
      m_cts[i] = tid;
      pthread_cond_init(&tid->m_s_cnd,NULL);
      pthread_cond_init(&tid->m_t_cnd,NULL);
      pthread_mutex_init(&tid->m_s_mtx,NULL);
      pthread_mutex_init(&tid->m_t_mtx,NULL);
      pthread_mutex_lock(&tid->m_s_mtx);
      pthread_mutex_lock(&tid->m_t_mtx);
      tid->m_s=1;
      int tec(0);
      if ((tec=pthread_create(&tid->m_id,NULL,&TCalcJL,(void*)tid)))
	THROW(fatal_error,"Cannot create thread "+ToString(i));
    }
  }
#endif
}

PS_Generator::~PS_Generator()
{
#ifdef USING__Threading
  for (size_t i(0);i<m_cts.size();++i) {
    CDBG_PG_TID *tid(m_cts[i]);
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
  CleanUp();
}

void PS_Generator::CleanUp()
{
  for (size_t i(0);i<m_cur.size();++i) 
    for (size_t j(0);j<m_cur[i].size();++j) delete m_cur[i][j]; 
  m_n=0;
  m_cl=Int_Matrix();
  m_cur=Current_Matrix();
  m_ctt=Current_Vector();
  m_swidths=m_smasses=Double_Vector();
  m_cmap=CB_MMap();
  m_cbmap=CB_Map();
  m_thmass=0.0;
}

size_t PS_Generator::NChannels() const
{
  size_t nch(0);
  for (size_t i(0);i<m_cur.size();++i)
    nch+=m_cur[i].size();
  return nch;
}    

void PS_Generator::SetColors(const Int_Vector &rc,
			     const Int_Vector &ac)
{
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"():\n";
#endif
  for (size_t i(0);i<m_n;++i) {
    m_cl[i][0]=rc[i];
    m_cl[i][1]=ac[i];
#ifdef DEBUG__BG
    msg_Debugging()<<"m_cl["<<i<<"][0] = "<<m_cl[i][0]
		   <<", m_cl["<<i<<"][1] = "<<m_cl[i][1]<<"\n";
#endif
  }
}

#ifdef USING__Threading
void *PS_Generator::TCalcJL(void *arg)
{
  CDBG_PG_TID *tid((CDBG_PG_TID*)arg);
  while (true) {
    // wait for generator to signal
    pthread_mutex_lock(&tid->m_s_mtx);
    pthread_mutex_unlock(&tid->m_s_mtx);
    pthread_cond_signal(&tid->m_s_cnd);
    if (tid->m_s==0) return NULL;
    // worker routine
    for (tid->m_i=tid->m_b;tid->m_i<tid->m_e;++tid->m_i) 
      tid->p_psg->m_cur[tid->m_n][tid->m_i]->Evaluate();
    // signal generator to continue
    pthread_cond_wait(&tid->m_t_cnd,&tid->m_t_mtx);
  }
  return NULL;
}
#endif

void PS_Generator::CalcJL()
{
  for (size_t i(0);i<m_cur[1].size();++i) 
    m_cur[1][i]->ConstructJ(Vec4D(),0,m_cl[i][0],m_cl[i][1]);
  if (m_zmode>0) {
    for (size_t n(2);n<m_n;++n) {
#ifdef USING__Threading
      if (m_cts.empty()) {
	for (size_t i(0);i<m_cur[n].size();++i) 
	  m_cur[n][i]->Evaluate();
      }
      else {
	// start calculator threads
	size_t d(m_cur[n].size()/m_cts.size());
	if (m_cur[n].size()%m_cts.size()>0) ++d;
	for (size_t j(0), i(0);j<m_cts.size()&&i<m_cur[n].size();++j) {
	  CDBG_PG_TID *tid(m_cts[j]);
	  tid->m_n=n;
	  tid->m_b=i;
	  tid->m_e=Min(i+=d,m_cur[n].size());
	  pthread_cond_wait(&tid->m_s_cnd,&tid->m_s_mtx);
	}
	// suspend calculator threads
	for (size_t j(0), i(0);j<m_cts.size()&&i<m_cur[n].size();++j) {
	  i+=d;
	  CDBG_PG_TID *tid(m_cts[j]);
	  pthread_mutex_lock(&tid->m_t_mtx);
	  pthread_mutex_unlock(&tid->m_t_mtx);
	  pthread_cond_signal(&tid->m_t_cnd);
	}
      }
#else
      for (size_t i(0);i<m_cur[n].size();++i) 
	m_cur[n][i]->Evaluate();
#endif
    }
    for (size_t n(m_n-2);n>=2;--n)
      for (size_t i(0);i<m_cur[n].size();++i)
	m_cur[n][i]->ResetZero();
  }
}

bool PS_Generator::Evaluate()
{
  if (m_zmode>0) {
    PHASIC::Process_Base *cur(p_xs->Process());
    while ((*cur)[0]!=cur) {
      double sum(0.0);
      Double_Vector sums;
      std::vector<PHASIC::Process_Base*> pss;
      for (size_t i(0);i<cur->Size();++i)
	if (!(*cur)[i]->Get<PHASIC::Single_Process>()->Zero()) {
	  if ((*cur)[i]->Integrator()->SumSqr()==0) sum+=1.0;
	  else sum+=(*cur)[i]->Integrator()->SumSqr();
	  sums.push_back(sum);
	  pss.push_back((*cur)[i]);
	}
      double disc(sum*ran->Get());
      for (size_t i(0);i<pss.size();++i)
	if (disc<=sums[i]) {
	  cur=pss[i];
	  break;
	}
    }
    SP(Color_Integrator) ci(cur->Integrator()->ColorIntegrator());
    if (ci==NULL) 
      THROW(fatal_error,"No color integrator for "+cur->Name());
    SetColors(ci->I(),ci->J());
  }
  CalcJL();
  return true;
}

bool PS_Generator::AddCurrent
(Current *const ref,const Flavour &fl,
 const size_t &n,const int mode,
 const double &m,const double &w,Current *const scc)
{
  Current_Key ckey(fl,NULL,0);
  PS_Current *cur(new PS_Current(ckey));
  if (mode&1) {
    cur->SetMass(m);
    cur->SetWidth(w);
  }
  cur->SetId(ref->Id());
  cur->SetKey(m_cur[n].size());
  cur->SetDirection(ref->Direction());
  cur->SetCut(ref->Cut());
  cur->SetOnShell(ref->OnShell());
  cur->SetOrderEW(ref->OrderEW());
  cur->SetOrderQCD(ref->OrderQCD());
  cur->SetNTChannel(ref->NTChannel());
  cur->SetSCC(scc);
  m_cur[n].push_back(cur);
  if (scc==NULL) m_tccs[cur->CId()].push_back(cur);
  m_cmap.insert(CB_Pair(ref,cur));
  m_cbmap.insert(CB_Pair(cur,ref));
  bool found(false);
  for (size_t i(0);i<m_ctt.size();++i)
    if (m_ctt[i]->PSInfo()==cur->PSInfo()) {
      found=true;
      break;
    }
  if (!found) m_ctt.push_back(cur);
  return found;
}

int PS_Generator::DecayType(const Current *jc,
			    const Current *ja,const Current *jb) const
{
  if (jc->CId()==(1<<m_n)-1-3) return 0;
  if ((jc->CId()&3)==1 || (jc->CId()&3)==2 ||
      jc->Flav().Mass()!=0.0 ||
      ja->Flav().Mass()!=0.0 || jb->Flav().Mass()!=0.0) return 0;
  if (ja->Flav()==jc->Flav() &&
      jb->Flav()==jc->Flav()) return 0;// temporary (->4)
  if (ja->Flav()==jb->Flav().Bar()) return 0;// temporary (->1)
  if (ja->Flav()==jc->Flav()) return 2;
  if (jb->Flav()==jc->Flav()) return 3;
  return 0;
}

class CB_PSSort {
public:

  bool operator()(const Vertex_Key &ka,const Vertex_Key &kb) const
  {
    if (ka.p_a<kb.p_a) return true;
    if (ka.p_a>kb.p_a) return false;
    if (ka.p_b<kb.p_b) return true;
    if (ka.p_b>kb.p_b) return false;
    if (ka.p_c<kb.p_c) return true;
    if (ka.p_c>kb.p_c) return false;
    return false;
  }

};

bool PS_Generator::Construct(Amplitude *const ampl)
{
  m_pmsinit=0;
  Current_Matrix curs(ampl->Currents());
  msg_Debugging()<<METHOD<<"(): '"<<ampl<<"' {\n";
  {
    msg_Indent();
  if (m_n>0) {
    if (m_n!=curs.size()) 
      THROW(fatal_error,"Invalid number of external particles");
  }
  else {
    m_n=curs.size();
    m_cl.resize(m_n,Int_Vector(2));
    m_cur.resize(m_n);
  }
  for (size_t n(1);n<m_n;++n) {
    for (size_t j(0);j<curs[n].size();++j) {
      if (curs[n][j]->Sub() ||
	  curs[n][j]->Flav().IsDummy()) continue;
      bool found(false);
      for (size_t i(0);i<m_cur[n].size();++i)
	if (m_cur[n][i]->Id()==curs[n][j]->Id() &&
	    (n==1 ||
	     (m_cur[n][i]->Flav().Mass()==curs[n][j]->Flav().Mass() &&
	      m_cur[n][i]->Flav().Width()==curs[n][j]->Flav().Width() &&
	      m_cur[n][i]->OrderEW()==curs[n][j]->OrderEW() &&
	      m_cur[n][i]->OrderQCD()==curs[n][j]->OrderQCD() &&
	      m_cur[n][i]->NTChannel()==curs[n][j]->NTChannel()))) {
	  Current *ref(m_cbmap[m_cur[n][i]]);
	  for (CB_MMap::const_iterator cit(m_cmap.lower_bound(ref));
	       cit!=m_cmap.upper_bound(ref);++cit) {
	    m_cmap.insert(CB_Pair(curs[n][j],cit->second));
	  }
	  found=true;
	  break;
	}
      if (found) {
	for (CB_MMap::const_iterator cit(m_cmap.lower_bound(curs[n][j]));
	     cit!=m_cmap.upper_bound(curs[n][j]);++cit) {
	  std::set<Vertex_Key,CB_PSSort> v3;
	  const Vertex_Vector &in(curs[n][j]->In());
	  const Vertex_Vector &rin(cit->second->In());
	  for (size_t i(0);i<in.size();++i) {
	    if (curs[n][j]->In()[i]->JE()) continue;
	    Current *ja(curs[n][j]->In()[i]->JA());
	    Current *jb(curs[n][j]->In()[i]->JB());
	    if (ja->PSInfo()<jb->PSInfo()) std::swap<Current*>(ja,jb);
	    for (CB_MMap::const_iterator ait(m_cmap.lower_bound(ja));
		 ait!=m_cmap.upper_bound(ja);++ait)
	      for (CB_MMap::const_iterator bit(m_cmap.lower_bound(jb));
		   bit!=m_cmap.upper_bound(jb);++bit) {
		Vertex_Key vkey(ait->second,bit->second,NULL,cit->second,NULL);
		int type(DecayType(curs[n][j],ja,jb));
		bool vf(false);
		for (size_t k(0);k<rin.size();++k) 
		  if (vkey.p_a==rin[k]->JA() && vkey.p_b==rin[k]->JB() &&
		      type==((PS_Vertex*)rin[k])->Type()) {
		    vf=true;
		    break;
		  }
		if (vf || v3.find(vkey)!=v3.end()) continue;
		v3.insert(vkey);
		PS_Vertex *vtx(new PS_Vertex(vkey));
		vtx->SetJA(vkey.p_a);
		vtx->SetJB(vkey.p_b);
		vtx->SetJC(vkey.p_c);
		vtx->SetType(type);
	      }
	  }
	}
	continue;
      }
      curs[n][j]->Print();
      if (curs[n][j]->Flav().Width()<s_pwmin &&
	  !curs[n][j]->Cut() && curs[n][j]->Flav().Mass()>0.0 &&
	  curs[n][j]->Flav().Mass()<m_chmass && n>1 && n<m_n-1) 
	AddCurrent(curs[n][j],curs[n][j]->Flav(),n,1);
      else AddCurrent(curs[n][j],curs[n][j]->Flav(),n);
      std::set<Vertex_Key,CB_PSSort> v3;
      const Vertex_Vector &in(curs[n][j]->In());
      for (size_t i(0);i<in.size();++i) {
	if (curs[n][j]->In()[i]->JE()) continue;
	Current *ja(curs[n][j]->In()[i]->JA());
	Current *jb(curs[n][j]->In()[i]->JB());
	if (ja->PSInfo()<jb->PSInfo()) std::swap<Current*>(ja,jb);
	for (CB_MMap::const_iterator ait(m_cmap.lower_bound(ja));
	     ait!=m_cmap.upper_bound(ja);++ait)
	  for (CB_MMap::const_iterator bit(m_cmap.lower_bound(jb));
	       bit!=m_cmap.upper_bound(jb);++bit)
	    for (CB_MMap::const_iterator cit(m_cmap.lower_bound(curs[n][j]));
		 cit!=m_cmap.upper_bound(curs[n][j]);++cit) {
	      Vertex_Key vkey(ait->second,bit->second,NULL,cit->second,NULL);
	      if (v3.find(vkey)!=v3.end()) continue;
	      v3.insert(vkey);
	      PS_Vertex *vtx(new PS_Vertex(vkey));
	      vtx->SetJA(vkey.p_a);
	      vtx->SetJB(vkey.p_b);
	      vtx->SetJC(vkey.p_c);
	      vtx->SetType(DecayType(curs[n][j],ja,jb));
	    }
      }
      for (CB_MMap::const_iterator cit(m_cmap.lower_bound(curs[n][j]));
	   cit!=m_cmap.upper_bound(curs[n][j]);++cit)
	cit->second->Print();
    }
  }
  }
  msg_Debugging()<<"}\n";
  for (size_t j(m_n-2);j>1;--j)
    for (Current_Vector::iterator cit(m_cur[j].begin());
	 cit!=m_cur[j].end();++cit)
      if ((*cit)->Dangling()) {
#ifdef DEBUG__BG
	msg_Debugging()<<"  delete current "<<**cit
		       <<", "<<(*cit)->Dangling()<<"\n";
#endif
	for (Current_Vector::iterator ctit(m_ctt.begin());
	     ctit!=m_ctt.end();++ctit)
	  if (*ctit==*cit) {
	    m_ctt.erase(ctit);
	    break;
	  }
	delete *cit;
	cit=--m_cur[j].erase(cit);
      }
  msg_Debugging()<<METHOD<<"(): Phase space statistics (n="
		 <<m_n<<") {\n  level currents vertices\n"<<std::right;
  double ismass(0.0);
  size_t csum(0), vsum(0), itmin(1);
  std::map<size_t,double> mmin;
  for (size_t i(1);i<m_n;++i) {
    csum+=m_cur[i].size();
    size_t cvsum(0);
    for (size_t j(0);j<m_cur[i].size();++j) {
      size_t cid(m_cur[i][j]->CId());
      double mass(i==1?m_cur[i][j]->Flav().Mass():m_cur[i][j]->Mass());
      if ((cid&3)==1 || (cid&3)==2) {
	if (i==1) ismass+=mass;
	mass=0.0;
      }
      for (size_t k(0);k<m_cur[i][j]->In().size();++k)
	mass=Max(mass,mmin[m_cur[i][j]->In()[k]->JA()->CId()]+
		 mmin[m_cur[i][j]->In()[k]->JB()->CId()]);
      if (mmin.find(cid)==mmin.end()) mmin[cid]=mass;
      else mmin[cid]=Max(mmin[cid],mass);
      if (i==m_n-1) {
	m_thmass=Max(m_thmass,Max(mass,ismass));
      }
#ifdef DEBUG__BG
      msg_Debugging()<<"  set min mass "
		     <<m_cur[i][j]->PSInfo()<<" -> "<<mass<<"\n";
#endif
      if ((cid==3 || cid==(size_t)((1<<m_n)-1-3)) &&
	  m_cur[i][j]->Mass()>0.0) {
	bool found(false);
	for (size_t k(0);k<m_smasses.size();++k)
	  if (m_smasses[k]==m_cur[i][j]->Mass()&&
	      m_swidths[k]==m_cur[i][j]->Width()) {
	    found=true;
	    break;
	  }
	if (!found) {
	  m_smasses.push_back(m_cur[i][j]->Mass());
	  m_swidths.push_back(m_cur[i][j]->Width());
	}
      }
#ifdef DEBUG__BG
      msg_Debugging()<<"  set min mass "
		     <<m_cur[i][j]->PSInfo()<<" -> "<<mass<<"\n";
#endif
      cvsum+=m_cur[i][j]->NIn();
    }
    msg_Debugging()<<"  "<<std::setw(5)<<i<<" "<<std::setw(8)
		   <<m_cur[i].size()<<" "<<std::setw(8)<<cvsum<<"\n";
    vsum+=cvsum;
    if (i>1 && m_cur[i].size()>0) itmin*=cvsum/m_cur[i].size();
  }
  msg_Debugging()<<std::left<<"} -> "<<csum<<" currents, "
		 <<vsum<<" vertices\n";
  itmin=Min(Max(itmin,m_itmin),m_itmax);
  msg_Tracking()<<METHOD<<"(): Set iteration minimum "<<itmin<<".\n";
  p_xs->Process()->Integrator()->SetItMin(itmin);
  msg_Tracking()<<METHOD<<"(): Masses {\n"
		<<"  threshold  : "<<m_thmass<<"\n";
  for (size_t i(0);i<m_smasses.size();++i)
    msg_Tracking()<<"  resonance "<<i<<": "
		  <<m_smasses[i]<<" / "<<m_swidths[i]<<"\n";
  msg_Tracking()<<"}\n";
  return true;
}

void PS_Generator::AddSTCC()
{
#ifdef DEBUG__BG
  DEBUG_FUNC("");
#endif
  for (size_t n(2);n<m_n-2;++n) {
    size_t oldsize(m_cur[n].size());
    for (size_t j(0);j<oldsize;++j) {
      if (((PS_Current*)m_cur[n][j])->SCC()) continue;
      size_t cid(m_cur[n][j]->CId());
      if ((cid&3)==0 || (cid&3)==3) continue;
      std::set<std::string> added;
      size_t pid(~3&cid);
      if (IdCount(pid)>1) {
	TCC_Map::const_iterator it(m_tccs.find(pid));
	if (it==m_tccs.end()) continue;
	for (size_t i(0);i<it->second.size();++i) {
	  if (added.empty()) 
	    ((PS_Current*)m_cur[n][j])->SetSCC(it->second[i]);
	  else if (added.find(it->second[i]->PSInfo())==added.end())
	    AddExtraCurrent(m_cur[n][j],n,m_cur[n][j]->Flav().Mass(),
			    m_cur[n][j]->Flav().Width(),it->second[i]);
	  added.insert(it->second[i]->PSInfo());
	}
      }
    }
  }
}

void PS_Generator::AddExtraCurrent
(Current *const cur,const size_t &n,
 const double &m,const double &w,Current *const scc)
{
  AddCurrent(cur,cur->Flav(),n,1,m,w);
#ifdef DEBUG__BG
  msg_Debugging()<<"  Add "<<m_cur[n].back()->PSInfo()
		 <<(scc?" ("+scc->PSInfo()+") ":"")<<" {\n";
#endif
  const Vertex_Vector &in(cur->In());
  for (size_t i(0);i<in.size();++i) {
    Vertex_Key vkey(in[i]->JA(),in[i]->JB(),NULL,m_cur[n].back(),NULL);
    PS_Vertex *vtx(new PS_Vertex(vkey));
    vtx->SetJA(vkey.p_a);
    vtx->SetJB(vkey.p_b);
    vtx->SetJC(vkey.p_c);
#ifdef DEBUG__BG
    msg_Debugging()<<"    "<<*vtx<<"\n";
#endif
  }
#ifdef DEBUG__BG
  msg_Debugging()<<"  } <-> {\n";
#endif
  const Vertex_Vector &out(cur->Out());
  for (size_t i(0);i<out.size();++i) {
    Current *ja(out[i]->JA()), *jb(out[i]->JB());
    if (ja==cur) ja=m_cur[n].back();
    else jb=m_cur[n].back();
    Vertex_Key vkey(ja,jb,NULL,out[i]->JC(),NULL);
    PS_Vertex *vtx(new PS_Vertex(vkey));
    vtx->SetJA(vkey.p_a);
    vtx->SetJB(vkey.p_b);
    vtx->SetJC(vkey.p_c);
#ifdef DEBUG__BG
    msg_Debugging()<<"    "<<*vtx<<"\n";
#endif
  }
#ifdef DEBUG__BG
  msg_Debugging()<<"  }\n";
#endif
}

void PS_Generator::SetPrefMasses(Cut_Data *const cuts)
{
  if (m_pmsinit) return;
  m_pmsinit=1;
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"() {\n";
#endif
  std::map<size_t,double> mmin;
  for (size_t n(2);n<m_n-2;++n) {
    size_t oldsize(m_cur[n].size());
    for (size_t j(0);j<oldsize;++j) {
      size_t cid(m_cur[n][j]->CId());
      size_t pid((cid&3)?(1<<m_n)-1-cid:cid);
      double psmin(sqrt(cuts->Getscut(PSId(pid))));
      double mass(Max(psmin,m_cur[n][j]->Mass()));
      if (m_cur[n][j]->OnShell()) {
	mmin[cid]=mass;
	continue;
      }
      if ((cid&3)==1 || (cid&3)==2) mass=0.0;
      for (size_t k(0);k<m_cur[n][j]->In().size();++k) {
	mass=Max(mass,mmin[m_cur[n][j]->In()[k]->JA()->CId()]+
		 mmin[m_cur[n][j]->In()[k]->JB()->CId()]);
      }
      if (mmin.find(cid)==mmin.end()) mmin[cid]=mass;
      else mmin[cid]=Max(mmin[cid],mass);
      if (m_cur[n][j]->Cut() || (cid&3)==1 || (cid&3)==2) continue;
      if ((m_ecmode&1) && mass>m_chmass &&
	  mass<rpa->gen.Ecms() && !IsEqual(mass,m_cur[n][j]->Mass()))
	AddExtraCurrent(m_cur[n][j],n,mass,0.0);
      if ((m_ecmode&2) && m_cur[n][j]->Mass()>m_chmass &&
	  m_cur[n][j]->Width()>s_pwmin && IdCount(pid)>2)
	AddExtraCurrent(m_cur[n][j],n,m_cur[n][j]->Mass(),0.0);
    }
  }
#ifdef DEBUG__BG
  msg_Debugging()<<"}\n";
#endif
}
