#include "COMIX/Amplitude/Amplitude.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/STL_Tools.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace COMIX;
using namespace ATOOLS;

static const double invsqrttwo(1.0/sqrt(2.0));

Amplitude::Amplitude():
  p_model(NULL), m_n(0), m_nf(6), m_oew(99), m_oqcd(99), m_ngpl(3),
  m_pmode('D')
{
  Data_Reader read(" ",";","!","=");
  std::string prec;
  if (!read.ReadFromFile(prec,"COMIX_PMODE")) prec="D";
  else msg_Tracking()<<METHOD<<"(): Set precision "<<prec<<".\n";
  if (prec!="D" && prec!="Q") THROW(not_implemented,"Invalid precision mode");
  m_pmode=prec[0];
  int helpi(0);
  if (!read.ReadFromFile(helpi,"COMIX_PG_MODE")) helpi=0;
  else msg_Info()<<METHOD<<"(): Set print graph mode "<<helpi<<".\n";
  m_pgmode=helpi;
  if (!read.ReadFromFile(helpi,"COMIX_VL_MODE")) helpi=0;
  else msg_Info()<<METHOD<<"(): Set vertex label mode "<<helpi<<".\n";
  Vertex_Base::SetVLMode(helpi);
  if (!read.ReadFromFile(helpi,"COMIX_N_GPL")) helpi=3;
  else msg_Info()<<METHOD<<"(): Set graphs per line "<<helpi<<".\n";
  m_ngpl=Max(1,Min(helpi,5));
#ifdef USING__Threading
  if (!read.ReadFromFile(helpi,"COMIX_ME_THREADS")) helpi=0;
  else msg_Tracking()<<METHOD<<"(): Set number of threads "<<helpi<<".\n";
  if (helpi>0) {
    m_cts.resize(helpi);
    for (size_t i(0);i<m_cts.size();++i) {
      CDBG_ME_TID *tid(new CDBG_ME_TID(this));
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

Amplitude::~Amplitude()
{
#ifdef USING__Threading
  for (size_t i(0);i<m_cts.size();++i) {
    CDBG_ME_TID *tid(m_cts[i]);
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
  }
#endif
  CleanUp();
}

size_t Amplitude::MakeId(const Int_Vector &ids,const int t)
{
  size_t id(0);
  if (t>0) {
    for (size_t i(0);i<ids.size();++i) 
      if (ids[i]>0) id+=1<<i;
  }
  else {
    for (size_t i(0);i<ids.size();++i) 
      if (ids[i]<0 || ids[i]==3) id+=1<<i;
  }
  return id;
}

Int_Vector Amplitude::MakeId(const size_t &id,const size_t &n)
{
  size_t ic(id);
  Int_Vector ids(n,0);
  for (size_t i(0);i<ids.size();++i) {
    size_t c(1<<i);
    if (ic&c) {
      ids[i]=1;
      ic-=c;
    }
  }
  if (ic!=0) THROW(fatal_error,"Invalid particle number");
  return ids;
}

void Amplitude::CleanUp()
{
  for (size_t i(0);i<m_cur.size();++i) 
    for (size_t j(0);j<m_cur[i].size();++j) delete m_cur[i][j]; 
  m_n=0;
  m_fl=Flavour_Vector();
  m_p=Vec4D_Vector();
  m_ch=Int_Vector();
  m_cl=Int_Matrix();
  m_cur=Current_Matrix();
  m_chirs=Int_Matrix();
  m_ress=QComplex_Vector();
  m_dirs=Int_Vector();
  m_cchirs=LongUIntMap_Matrix();
  m_combs.clear();
  m_flavs.clear();
}

bool Amplitude::MatchIndices(const Int_Vector &ids,const size_t &n,
				  const size_t &i,const size_t &j,
				  const size_t &k)
{
  for (size_t l(0);l<n;++l) {
    bool found(false), twice(false);
    for (size_t m(0);m<m_cur[i][j]->Id().size();++m) 
      if (m_cur[i][j]->Id()[m]==ids[l]) {
	if (found) twice=true;
	found=true;
      }
    for (size_t m(0);m<m_cur[n-i][k]->Id().size();++m) 
      if (m_cur[n-i][k]->Id()[m]==ids[l]) {
	if (found) twice=true;
	found=true;
      }
    if (!found || twice) return false;
  }
  return true;
}

int Amplitude::CheckDecay(const ATOOLS::Flavour &fl,
			  const Int_Vector &ids) const
{
  size_t cid(0);
  if (m_decid.empty()) return 0;
  for (size_t i(0);i<ids.size();++i) cid+=1<<ids[i];
  for (size_t i(0);i<m_decid.size();++i) {
    size_t did(m_decid[i].m_id);
    if (did&(1<<0)) did=(1<<m_n)-1-did;
    if (did==cid) {
      if (fl==m_decid[i].m_fl || fl.Bar()==m_decid[i].m_fl) return i+1;
//       msg_Debugging()<<"delete prop "<<fl<<" "<<ids<<" "<<cid
//  		     <<", requested "<<m_decid[i].m_fl<<" "<<did<<"\n";
      return -1;
    }
    if (!((did&cid)==0 || (did&cid)==cid || (did&cid)==did)) {
//       msg_Debugging()<<"delete prop "<<fl<<" "<<ids<<" "<<cid
//  		     <<", requested "<<m_decid[i].m_fl<<" "<<did<<"\n";
      return -1;
    }
  }
  return 0;
}

void Amplitude::AddCurrent(const Int_Vector &ids,const size_t &n,
				const Flavour &fl,const int dir)
{
  // add new currents
  size_t oewmax(0), oqcdmax(0);
  int dec(CheckDecay(fl,ids));
  if (dec<0) return;
  std::map<std::string,Current_Base*> curs;
  Current_Key ckey(dir>0?fl.Bar():fl,p_model);
  Current_Base *cur(Current_Getter::GetObject
		    (std::string(1,m_pmode)+ckey.Type(),ckey));
  if (cur==NULL) return;
  cur->SetDirection(dir);
  if (dec!=0) {
    cur->SetCut(m_decid[dec-1].m_nmax);
    cur->SetOnShell(m_decid[dec-1].m_osd);
  }
  std::set<Vertex_Key> v3;
  // compose current from all possible subcurrents
  for (size_t i(1);i<n;++i) {
    for (size_t j(0);j<m_cur[i].size();++j) {
      for (size_t k(0);k<m_cur[n-i].size();++k) {
	if (!MatchIndices(ids,n,i,j,k)) continue;
	Vertex_Key vkey(m_cur[i][j],m_cur[n-i][k],cur,p_model);
	if (v3.find(vkey.SwapAB())!=v3.end()) continue;
	Vertex_Base *v(Vertex_Getter::GetObject
		       (std::string(1,m_pmode)+vkey.ID(),vkey));
	if (v!=NULL) {
	  size_t oew(vkey.p_a->OrderEW()+
		     vkey.p_b->OrderEW()+v->OrderEW());
	  size_t oqcd(vkey.p_a->OrderQCD()+
		      vkey.p_b->OrderQCD()+v->OrderQCD());
	  if (!v->Active() || oew>m_oew || oqcd>m_oqcd) {
#ifdef DEBUG__BG
	    msg_Debugging()<<"delete vertex {"<<vkey.p_a->Flav()<<",("
			   <<vkey.p_a->OrderEW()<<","<<vkey.p_a->OrderQCD()
			   <<")}{"<<vkey.p_b->Flav()<<",("
			   <<vkey.p_b->OrderEW()<<","<<vkey.p_b->OrderQCD()
			   <<")}-"<<v->Tag()<<"("<<v->OrderEW()<<","
			   <<v->OrderQCD()<<")->{"<<cur->Flav()<<"} => ("
			   <<oew<<","<<oqcd<<") vs. max = ("<<m_oew<<","
			   <<m_oqcd<<"), act = "<<v->Active()<<"\n";
#endif
	    delete v;
	    continue;
	  }
	  std::string okey("("+ToString(oew)+","+ToString(oqcd)+")");
	  if (oew!=cur->OrderEW() || oqcd!=cur->OrderQCD()) {
	    std::map<std::string,Current_Base*>::iterator 
	      cit(curs.find(okey));
	    if (cit!=curs.end()) cur=cit->second;
	    else {
	      if (cur->OrderEW()>0 || cur->OrderQCD()>0)
		cur=Current_Getter::GetObject
		  (std::string(1,m_pmode)+ckey.Type(),ckey);
	      if (n<m_n-1) {
		if (dec!=0) {
		  cur->SetCut(m_decid[dec-1].m_nmax);
		  cur->SetOnShell(m_decid[dec-1].m_osd);
		}
		cur->SetOrderEW(oew);
		cur->SetOrderQCD(oqcd);
		curs[okey]=cur;
	      }
	      else {
		oewmax=Max(oewmax,oew); 
		oqcdmax=Max(oqcdmax,oqcd);
	      }
	    }
	  }
	  v->SetJA(vkey.p_a);
	  v->SetJB(vkey.p_b);
	  v->SetJC(cur);
	  v3.insert(Vertex_Key(vkey.p_a,vkey.p_b,cur,p_model));
	}
      }
    }
  }
  if (v3.empty() && n>1) {
    delete cur;
    return;
  }
  if (n==1 || n==m_n-1) curs[""]=cur;
  if (n==m_n-1) {
    m_oew=oewmax;
    m_oqcd=oqcdmax;
  }
  Int_Vector isfs(ids.size());
  for (size_t i(0);i<ids.size();++i)
    isfs[i]=m_fl[ids[i]].IsFermion();
  for (std::map<std::string,Current_Base*>::iterator 
	 cit(curs.begin());cit!=curs.end();++cit) {
    cit->second->SetId(ids);
    cit->second->SetFId(isfs);
    cit->second->FindPermutations();
    cit->second->SetKey(m_cur[n].size());
    m_cur[n].push_back(cit->second);
    cit->second->Print();
  }
}

bool Amplitude::Construct(Flavour_Vector &fls,
			       Int_Vector ids,const size_t &n)
{
  if (ids.size()==n) {
    if (n==m_n-1) {
      if (!m_fl.front().IsOn()) return false;
      AddCurrent(ids,n,m_fl.front().Bar(),m_dirs.front());
    }
    else {
      for (size_t i(0);i<fls.size();++i) {
	AddCurrent(ids,n,fls[i],0);
	if (fls[i].Bar()!=fls[i])
	  AddCurrent(ids,n,fls[i].Bar(),0);
      }
    }
    return true;
  }
  // currents are unordered -> use ordered indexing
  size_t last(ids.empty()?0:ids.back());
  ids.push_back(0);
  if (n==m_n-1) {
    // calculate only one final current
    ids.back()=last+1;
    if (!Construct(fls,ids,n)) return false;
    return m_cur.back().size();
  }
  for (size_t i(last+1);i<m_n;++i) {
    // fill currents 0..n-1 for external partons
    // currents 0..n-2 each index for internal partons
    ids.back()=i;
    if (!Construct(fls,ids,n) && n==1) return false;
  }
  return true;
}

bool Amplitude::Construct(const Flavour_Vector &flavs)
{
  m_fl=flavs;
  m_n=m_fl.size();
  m_p.resize(m_n);
  m_ch.resize(m_n);
  m_cl.resize(m_n,Int_Vector(2));
  m_cur.resize(m_n);
  Int_Vector ids(1);
  Flavour_Vector fls(p_model->IncludedFlavours());
  for (size_t i(0);i<m_n;++i) {
    ids.back()=i;
    if (!m_fl[i].IsOn()) return false;
    AddCurrent(ids,1,m_fl[i],m_dirs[i]);
  }
  ids.clear();
  for (size_t i(2);i<m_n;++i)
    if (!Construct(fls,ids,i)) return false;
  for (size_t j(m_n-2);j>1;--j)
    for (Current_Vector::iterator cit(m_cur[j].begin());
	 cit!=m_cur[j].end();++cit)
      if ((*cit)->Dangling()) {
#ifdef DEBUG__BG
	msg_Debugging()<<"delete current "<<**cit<<", "<<(*cit)->Dangling()
		       <<", O("<<(*cit)->OrderEW()<<","<<(*cit)->OrderQCD()
		       <<") vs. O_{max}("<<m_oew<<","<<m_oqcd<<")\n";
#endif
	delete *cit;
	cit=--m_cur[j].erase(cit);
      }
  FillCombinations();
  msg_Debugging()<<METHOD<<"(): Amplitude statistics (n="
		 <<m_n<<") {\n  level currents vertices\n"<<std::right;
  size_t csum(0), vsum(0);
  for (size_t i(1);i<m_n;++i) {
    csum+=m_cur[i].size();
    size_t cvsum(0);
    for (size_t j(0);j<m_cur[i].size();++j) cvsum+=m_cur[i][j]->NIn();
    msg_Debugging()<<"  "<<std::setw(5)<<i<<" "<<std::setw(8)
		   <<m_cur[i].size()<<" "<<std::setw(8)<<cvsum<<"\n";
    vsum+=cvsum;
  }
  msg_Debugging()<<std::left<<"} -> "<<csum<<" currents, "
		 <<vsum<<" vertices"<<std::endl;
  return true;
}

void Amplitude::FillCombinations()
{
  msg_Debugging()<<METHOD<<"(): {\n";
  msg_Debugging()<<"  flavours {\n";
  for (size_t i(2);i<m_n-1;++i)
    for (size_t j(0);j<m_cur[i].size();++j) {
      size_t id(m_cur[i][j]->CId());
      m_flavs[id].push_back(m_cur[i][j]->Flav());
      m_flavs[(1<<m_n)-1-id].push_back(m_cur[i][j]->Flav().Bar());
      msg_Debugging()<<"    "<<ID(id)<<" / "<<ID((1<<m_n)-1-id)
		     <<" -> "<<m_cur[i][j]->Flav()<<"\n";
    }
  msg_Debugging()<<"  } -> "<<m_flavs.size()<<"\n";
  msg_Debugging()<<"  combinations {\n";
  for (size_t i(2);i<m_n;++i)
    for (size_t j(0);j<m_cur[i].size();++j) {
      Vertex_Vector ins(m_cur[i][j]->In());
      for (size_t k(0);k<ins.size();++k) {
	size_t ida(ins[k]->JA()->CId());
	size_t idb(ins[k]->JB()->CId());
	size_t idc((1<<m_n)-1-ins[k]->JC()->CId());
	msg_Debugging()<<"    "<<ID(ida)
		       <<" "<<ID(idb)<<" "<<ID(idc)<<"\n";
	m_combs.insert(std::pair<size_t,size_t>(ida,idb));
	m_combs.insert(std::pair<size_t,size_t>(idb,ida));
	m_combs.insert(std::pair<size_t,size_t>(idb,idc));
	m_combs.insert(std::pair<size_t,size_t>(idc,idb));
	m_combs.insert(std::pair<size_t,size_t>(idc,ida));
	m_combs.insert(std::pair<size_t,size_t>(ida,idc));
      }
    }
  msg_Debugging()<<"  } -> "<<m_combs.size()<<"\n";
  msg_Debugging()<<"}\n";
}

bool Amplitude::Map(const Amplitude &ampl,Flavour_Map &flmap)
{
  flmap.clear();
  msg_Debugging()<<METHOD<<"(): {\n";
  size_t svlmode(Vertex_Base::VLMode());
  Vertex_Base::SetVLMode(7);
  for (size_t n(1);n<m_n;++n) {
    if (ampl.m_cur[n].size()!=m_cur[n].size()) {
      msg_Debugging()<<"  current count differs\n} no match\n";
      Vertex_Base::SetVLMode(svlmode);
      flmap.clear();
      return false;
    }
    for (size_t i(0);i<m_cur[n].size();++i) {
      msg_Debugging()<<"  check m_cur["<<n<<"]["<<i<<"] {\n";
      if (flmap.find(ampl.m_cur[n][i]->Flav())==flmap.end()) {
	msg_Debugging()<<"    mapped ["<<n<<"]["<<i<<"] "
		       <<ampl.m_cur[n][i]->Flav()
		       <<" -> "<<m_cur[n][i]->Flav()<<"\n";
	if (ampl.m_cur[n][i]->Flav().IsAnti()^
	    m_cur[n][i]->Flav().IsAnti()) {
	  msg_Debugging()<<"    particle type differs\n  }\n";
	  msg_Debugging()<<"} no match\n";
          Vertex_Base::SetVLMode(svlmode);
	  flmap.clear();
	  return false;
	}
	if (ampl.m_cur[n][i]->Flav().StrongCharge()!=
	    m_cur[n][i]->Flav().StrongCharge()) {
	  msg_Debugging()<<"    color structure differs\n  }\n";
	  msg_Debugging()<<"} no match\n";
          Vertex_Base::SetVLMode(svlmode);
	  flmap.clear();
	  return false;
	}
	flmap[ampl.m_cur[n][i]->Flav()]=m_cur[n][i]->Flav();
      }
      else {
	if (ampl.m_cur[n][i]->Flav()!=ampl.m_cur[n][i]->Flav()) {
	  msg_Debugging()<<"    current differs\n  }\n";
	  msg_Debugging()<<"} no match\n";
	  Vertex_Base::SetVLMode(svlmode);
	  flmap.clear();
	  return false;
	}
      }
      Vertex_Vector vin(m_cur[n][i]->In());
      Vertex_Vector avin(ampl.m_cur[n][i]->In());
      if (avin.size()!=vin.size()) {
	msg_Debugging()<<"    vertex count differs\n  }\n";
	msg_Debugging()<<"} no match\n";
	Vertex_Base::SetVLMode(svlmode);
	flmap.clear();
	return false;
      }
      for (size_t j(0);j<vin.size();++j) {
#ifdef DEBUG__BG
	msg_Debugging()<<"    check m_in["<<j<<"] {\n";
#endif
	if (!avin[j]->Map(*vin[j])) {
	  msg_Debugging()<<"    } no match\n  }\n} no match\n";
	  Vertex_Base::SetVLMode(svlmode);
	  flmap.clear();
	  return false;
	}
#ifdef DEBUG__BG
	msg_Debugging()<<"    }\n";
#endif
      }
      msg_Debugging()<<"  }\n";
    }
  }
  msg_Debugging()<<"} matched\n";
  Vertex_Base::SetVLMode(svlmode);
  return true;
}

#ifdef USING__Threading
void *Amplitude::TCalcJL(void *arg)
{
  CDBG_ME_TID *tid((CDBG_ME_TID*)arg);
  while (true) {
    // wait for amplitude to signal
    pthread_mutex_lock(&tid->m_s_mtx);
    pthread_mutex_unlock(&tid->m_s_mtx);
    pthread_cond_signal(&tid->m_s_cnd);
    if (tid->m_s==0) return NULL;
    // worker routine
    for (tid->m_i=tid->m_b;tid->m_i<tid->m_e;++tid->m_i) 
      tid->p_ampl->m_cur[tid->m_n][tid->m_i]->Evaluate();
    // signal amplitude to continue
    pthread_cond_wait(&tid->m_t_cnd,&tid->m_t_mtx);
  }
  return NULL;
}
#endif

void Amplitude::CalcJL()
{
  for (size_t i(0);i<m_cur[1].size();++i) 
    m_cur[1][i]->ConstructJ(m_p[i],m_ch[i],m_cl[i][0],m_cl[i][1]);
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
	CDBG_ME_TID *tid(m_cts[j]);
	tid->m_n=n;
	tid->m_b=i;
	tid->m_e=Min(i+=d,m_cur[n].size());
	pthread_cond_wait(&tid->m_s_cnd,&tid->m_s_mtx);
      }
      // suspend calculator threads
      for (size_t j(0), i(0);j<m_cts.size()&&i<m_cur[n].size();++j) {
	i+=d;
	CDBG_ME_TID *tid(m_cts[j]);
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
}

void Amplitude::ResetZero()
{
  for (size_t n(m_n-2);n>=2;--n) {
    for (size_t i(0);i<m_cur[n].size();++i) 
      m_cur[n][i]->ResetZero();
  }
}

void Amplitude::SetMomenta(const Vec4D_Vector &moms)
{
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"():\n";
#endif
  Vec4D sum;
  for (size_t i(0);i<m_n;++i) {
    m_p[i]=m_dirs[i]>0?-moms[i]:moms[i];
    sum+=m_p[i];
#ifdef DEBUG__BG
    msg_Debugging()<<"set p["<<i<<"] = "<<m_p[i]
		   <<" ("<<sqrt(dabs(m_p[i].Abs2()))<<")\n";
#endif
  }
  static double accu(sqrt(Accu()));
  if (!IsEqual(sum,Vec4D(),accu)) 
    msg_Error()<<METHOD<<"(): Four momentum not conserved. sum = "
	       <<sum<<"."<<std::endl;
}

void Amplitude::SetColors(const Int_Vector &rc,
			       const Int_Vector &ac,const bool set)
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
  if (set) {
    Vertex_Base::SetCIMin(1);
    Vertex_Base::SetCIMax(0);
  }
  else {
    Vertex_Base::SetCIMin(1);
    Vertex_Base::SetCIMax(3);
  }
}

bool Amplitude::Evaluate(const Int_Vector &chirs)
{
  for (size_t j(0);j<m_n;++j) m_ch[j]=chirs[j];
  CalcJL();
  size_t ihp(MakeId(m_ch,1)), ihm(MakeId(m_ch,-1));
  QComplex res;
  if (m_pmode=='D') res=m_cur[1].front()->
    Contract<double>(*m_cur.back()[0],ihm,ihp);
  else if (m_pmode=='Q') res=m_cur[1].front()->
    Contract<long double>(*m_cur.back()[0],ihm,ihp);
  else THROW(not_implemented,"Internal error");
#ifdef DEBUG__BG
  msg_Debugging()<<"A"<<chirs<<" = "<<res<<" "
		 <<std::abs(res)<<" {"<<ihm<<","<<ihp<<"}\n";
#endif
  m_res=(res*std::conj(res)).real();
  return true;
}

bool Amplitude::EvaluateAll()
{
  for (size_t j(0);j<m_n;++j) m_ch[j]=0;
  CalcJL();
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"(): "<<m_ress.size()<<" amplitudes {\n";
#endif
  if (m_pmode=='D') {
    DComplex_Vector ress(m_ress.size(),DComplex(0.0));
    m_cur[1].front()->Contract(*m_cur.back()[0],m_cchirs,ress);
    for (size_t i(0);i<m_ress.size();++i) m_ress[i]=ress[i];
  }
  else if (m_pmode=='Q') {
    for (size_t i(0);i<m_ress.size();++i) m_ress[i]=0.0;
    m_cur[1].front()->Contract(*m_cur.back()[0],m_cchirs,m_ress);
  }
  else {
    THROW(not_implemented,"Internal error");
  }
#ifdef DEBUG__BG
  msg_Debugging()<<"}\n";
#endif
  double csum(0.0);
  for (size_t j(0);j<m_ress.size();++j) {
#ifdef DEBUG__BG
    msg_Debugging()<<"A["<<j<<"]"<<m_chirs[j]
		   <<" = "<<m_ress[j]<<" -> "<<std::abs(m_ress[j])<<"\n";
#endif
    csum+=(m_ress[j]*std::conj(m_ress[j])).real();
  }
  m_res=csum;
  return true;
}

bool Amplitude::CheckChirs(const Int_Vector &chirs)
{
  Int_Vector q(m_nf+1,0);
  size_t p(0), m(0), mp(0);
  if (m_oew>0) return true;
  for (size_t i(0);i<chirs.size();++i) {
    if (m_fl[i].IsMassive()) ++mp;
    if (m_fl[i].IsQuark() && !m_fl[i].IsMassive()) 
      q[m_fl[i].Kfcode()]+=chirs[i];
    if (chirs[i]>0) ++p;
    else if (chirs[i]<0) ++m;
    else THROW(fatal_error,"Invalid helicities");
  }
  for (size_t i(0);i<q.size();++i) 
    if (q[i]!=0) return false;
  return mp>0 || (p>1 && m>1);
}

bool Amplitude::ConstructChirs(Int_Vector chirs,const size_t &i)
{
  if (i==chirs.size()) {
    if (CheckChirs(chirs)) {
#ifdef DEBUG__BG
      msg_Debugging()<<METHOD<<"(): Add configuration "<<chirs<<"\n";
#endif
      m_chirs.push_back(chirs);
      m_ress.push_back(0.0);
    }
    return true;
  }
  if (chirs[i]!=0) {
    ConstructChirs(chirs,i+1);
  }
  else {
    switch (m_fl[i].IntSpin()) {
    case 0:
      ConstructChirs(chirs,i+1);
      break;
    case 1:
      for (int ch(1);ch>=-1;ch-=2) {
	chirs[i]=ch;
	ConstructChirs(chirs,i+1);
      }
      break;
    case 2:
      for (int ch(m_fl[i].IsMassive()?3:1);ch>=-1;ch-=2) {
	chirs[i]=ch;
	ConstructChirs(chirs,i+1);
      }
      break;
    default:
      THROW(not_implemented,"Cannot handle spin "+
	    ToString(m_fl[i].Spin())+" particles");
    }
  }
  return true;
}

bool Amplitude::Construct(const Int_Vector &incs,
			       const Flavour_Vector &flavs,
			       Model *const model)
{
  CleanUp();
  p_model=model;
  m_dirs=incs;
  if (!Construct(flavs)) return false;
  Int_Vector chirs(flavs.size(),0);
  if (!ConstructChirs(chirs,0)) return false;
  for (size_t i(0);i<m_ress.size();++i) {
    size_t ihp(MakeId(m_chirs[i],1)), ihm(MakeId(m_chirs[i],-1));
    m_cchirs[ihm][ihp]=i;
#ifdef DEBUG__BG
    msg_Debugging()<<"map A"<<m_chirs[i]<<" -> {"<<ihm
		   <<MakeId(ihm,m_n)<<","<<ihp<<MakeId(ihp,m_n)<<"}\n";
#endif
  }
  return true;
}

bool Amplitude::Combinable(const size_t &idi,const size_t &idj) const
{
  Combination_Set::const_iterator 
    cit(m_combs.find(std::pair<size_t,size_t>(idi,idj)));
  return cit!=m_combs.end();
}

const ATOOLS::Flavour_Vector &
Amplitude::CombinedFlavour(const size_t &idij) const
{
  CFlavVector_Map::const_iterator fit(m_flavs.find(idij));
  if (fit==m_flavs.end()) THROW(fatal_error,"Invalid request");
  return fit->second;
}

void Amplitude::SetGauge(const size_t &n)
{
  Vec4D k(1.0,0.0,1.0,0.0);
  switch(n) {
  case 1: k=Vec4D(1.0,1.0,0.0,0.0); break;
  case 2: k=Vec4D(1.0,invsqrttwo,invsqrttwo,0.0); break;
  case 3: k=Vec4D(1.0,invsqrttwo,-invsqrttwo,0.0); break;
  case 4: k=Vec4D(1.0,invsqrttwo,0.0,invsqrttwo); break;
  case 5: k=Vec4D(1.0,invsqrttwo,0.0,-invsqrttwo); break;
  case 6: k=Vec4D(1.0,0.0,invsqrttwo,invsqrttwo); break;
  case 7: k=Vec4D(1.0,0.0,invsqrttwo,-invsqrttwo); break;
  }
  for (size_t j(1);j<m_cur.size();++j)
    for (size_t i(0);i<m_cur[j].size();++i) m_cur[j][i]->SetGauge(k);
}

bool Amplitude::GaugeTest(const Vec4D_Vector &moms)
{
  msg_Tracking()<<METHOD<<"(): Performing gauge test ..."<<std::flush;
  msg_Indent();
  SetGauge(0);
  SetMomenta(moms);
  if (!EvaluateAll()) return false;
  QComplex_Vector ress(m_ress);
  SetGauge(1);
  SetMomenta(moms);
  if (!EvaluateAll()) return false;
  double mean(0.0);
  for (size_t i(0);i<m_ress.size();++i) mean+=std::abs(ress[i]);
  mean/=m_ress.size();
  msg_Debugging()<<METHOD<<"(): {\n";
#ifdef USE__Strict_QCD_Gauge_Test
  if (m_oew>0) {
#else
  if (true) {
#endif
    double xs(0.0), rxs(0.0);
    for (size_t i(0);i<m_ress.size();++i) {
      xs+=sqr(std::abs(m_ress[i]));
      rxs+=sqr(std::abs(ress[i]));
      msg_Debugging()<<"A("<<m_chirs[i]<<") = "<<std::abs(m_ress[i])
		     <<" vs. "<<std::abs(ress[i])<<"\n";
    }
    msg_Debugging()<<"\\sigma_{tot} = "<<xs<<" vs. "<<rxs
		   <<" -> dev. "<<xs/rxs-1.0<<"\n";
    if (!IsEqual(xs,rxs)) {
      msg_Error().precision(12);
      msg_Error()<<"\n"<<METHOD<<"(): Large deviation {\n      "
		  <<std::setw(18)<<std::right<<xs<<"\n   vs "
		  <<std::setw(18)<<rxs<<"\n   => "<<std::setw(18)
		  <<(xs/rxs-1.0)<<"\n}"<<std::left<<std::endl;
      msg_Error().precision(6);
      return true;
    }
    if (!IsEqual(xs,rxs,rpa.gen.Accu())) {
      msg_Error().precision(12);
      msg_Error()<<"\n"<<METHOD<<"(): Gauge test failed {\n      "
		  <<std::setw(18)<<std::right<<xs<<"\n   vs "
		  <<std::setw(18)<<rxs<<"\n   => "<<std::setw(18)
		  <<(xs/rxs-1.0)<<"\n}"<<std::left<<std::endl;
      msg_Error().precision(6);
      return false;
    }
  }
  else {
    for (size_t i(0);i<m_ress.size();++i) {
      msg_Debugging()<<"A("<<m_chirs[i]
		     <<") = "<<m_ress[i]<<" vs. "<<ress[i]<<" -> dev. "
		     <<m_ress[i].real()/ress[i].real()-1.0<<" "
		     <<m_ress[i].imag()/ress[i].imag()-1.0<<"\n";
      double accu(sqrt(Accu()));
      if (!IsEqual(m_ress[i].real(),ress[i].real(),accu) ||
	  !IsEqual(m_ress[i].imag(),ress[i].imag(),accu)) {
	double rrat(mean/Max(dabs(m_ress[i].real()),
			     dabs(ress[i].real()))*Accu());
	double irat(mean/Max(dabs(m_ress[i].imag()),
			     dabs(ress[i].imag()))*Accu());
	if ((IsEqual(m_ress[i].real(),ress[i].real(),rrat) ||
	     (m_ress[i].real()==0.0 && ress[i].real()==0.0) ||
	     (IsZero(m_ress[i].real(),rrat) && IsZero(ress[i].real(),rrat)))&&
	    (IsEqual(m_ress[i].imag(),ress[i].imag(),irat) ||
	     (m_ress[i].imag()==0.0 && ress[i].imag()==0.0) ||
	     (IsZero(m_ress[i].imag(),irat) && IsZero(ress[i].imag(),irat)))) {
	  msg_Error().precision(12);
	  msg_Tracking()
	    <<METHOD<<"(): Large deviation for small numbers {\n"
	    <<"      ("<<std::setw(18)<<std::right<<m_ress[i].real()
	    <<","<<std::setw(18)<<m_ress[i].imag()<<")\n   vs ("
	    <<std::setw(18)<<ress[i].real()<<","
	    <<std::setw(18)<<ress[i].imag()<<")\n   => ("
	    <<std::setw(18)<<(m_ress[i].real()/ress[i].real()-1.0)
	    <<","<<std::setw(18)
	    <<(m_ress[i].imag()/ress[i].imag()-1.0)<<")\n  ref {"
	    <<std::setw(18)<<rrat<<","<<std::setw(18)<<irat
	    <<"}\n}"<<std::left<<std::endl;
	  msg_Error().precision(6);
	}
	else {
	  msg_Error().precision(12);
	  msg_Error()
	    <<"\n"<<METHOD<<"(): Gauge test failed {\n"
	    <<"      ("<<std::setw(18)<<std::right<<m_ress[i].real()
	    <<","<<std::setw(18)<<m_ress[i].imag()<<")\n   vs ("
	    <<std::setw(18)<<ress[i].real()<<","
	    <<std::setw(18)<<ress[i].imag()<<")\n   => ("
	    <<std::setw(18)<<(m_ress[i].real()/ress[i].real()-1.0)
	    <<","<<std::setw(18)
	    <<(m_ress[i].imag()/ress[i].imag()-1.0)<<")\n  ref {"
	    <<std::setw(18)<<rrat<<","<<std::setw(18)<<irat
	    <<"}\n}"<<std::left<<std::endl;
	  msg_Error().precision(6);
 	  return false;
	}
      }
    }
  }
  msg_Debugging()<<"}\n";
  msg_Tracking()<<"satisfied."<<std::endl;
  return true;
}

void Amplitude::WriteOutGraph
(std::ostream &str,Graph_Node *graph,size_t &ng,
 std::set<std::string> &cvs) const
{
  if ((*graph)->empty()) {
    size_t fp(0), nf(0);
    for (size_t j(1);j<graph->size();++j)
      if ((*graph)[j].find("%%")==std::string::npos) {
	std::string cl((*graph)[j]);
	size_t bpos(cl.find("F="));
	if (bpos!=std::string::npos) {
	  ++nf;
	  size_t epos(bpos+=2);
	  for (;cl[epos]>=48 && cl[epos]<=57;++epos);
	  fp+=ToType<size_t>(cl.substr(bpos,epos-bpos));
	}
	bpos=cl.find("T=");
	if (bpos!=std::string::npos) {
	  size_t epos(cl.find(')',bpos+=2));
	  cvs.insert(cl.substr(bpos,epos-bpos+1));
	}
      }
    str<<"  \\parbox{"<<(5*m_n+20)<<"mm}{\\begin{center}\n  Graph "<<++ng;
    if (nf>0) str<<", $\\sum \\rm F$="<<fp<<" ("<<(fp%2==0?'+':'-')<<")";
    str<<"\\\\[6mm]\n  \\begin{fmfgraph*}("<<(10*m_n)<<","<<(10*m_n)<<")\n";
    str<<"    \\fmfsurround{"<<graph->front()<<"}\n";
    for (size_t j(0);j<m_n;++j) 
      str<<"    \\fmfv{decor.size=0ex,label=$J_{"
	 <<j<<"}$}{j_"<<(1<<j)<<"}\n";
    for (size_t j(1);j<graph->size();++j)
      if ((*graph)[j].find("%%")==std::string::npos) 
	str<<(*graph)[j]<<"\n";
    str<<"  \\end{fmfgraph*}\\end{center}\\vspace*{5mm}}";
    if (ng>0 && ng%m_ngpl==0) str<<" \\\\\n\n";
    else str<<" &\n\n";
  }
  else {
    for (size_t i(0);i<(*graph)->size();++i)
      WriteOutGraph(str,(*graph)()[i],ng,cvs);
  }
}

void Amplitude::WriteOutGraphs(const std::string &file) const
{
  Graph_Node graphs("j_1",true);
  graphs.push_back("    %% "+graphs.back());
  m_cur.back().front()->CollectGraphs(&graphs);
  std::ofstream str(file.c_str());
  str<<"\\documentclass[a4paper]{article}\n\n";
  str<<"\\usepackage{feynmp}\n";
  str<<"\\usepackage{amsmath}\n";
  str<<"\\usepackage{amssymb}\n";
  str<<"\\usepackage{longtable}\n";
  str<<"\\usepackage{pst-text}\n\n";
  str<<"\\setlength{\\headheight}{-1cm}\n";
  str<<"\\setlength{\\headsep}{0mm}\n";
  str<<"\\setlength{\\oddsidemargin}{-1in}\n";
  str<<"\\setlength{\\evensidemargin}{-1in}\n";
  str<<"\\setlength{\\textheight}{28truecm}\n";
  str<<"\\setlength{\\textwidth}{21truecm}\n\n";
  str<<"\\newcommand{\\p}{+}\n";
  str<<"\\newcommand{\\m}{-}\n\n";
  str<<"\\begin{document}\n";
  std::string texfile(file.substr(0,file.find(".")));
  if (texfile.rfind("/")!=std::string::npos)
    texfile=texfile.substr(texfile.rfind("/")+1);
  str<<"\\begin{fmffile}{"<<texfile<<"_fg}\n\n";
  str<<"  \\fmfset{thick}{1.25thin}\n";
  str<<"  \\fmfset{arrow_len}{2mm}\n";
  str<<"  \\fmfset{curly_len}{1.5mm}\n";
  str<<"  \\fmfset{wiggly_len}{1.5mm}\n";
  str<<"  \\fmfset{dot_len}{1mm}\n";
  str<<"  \\unitlength=0.5mm\n\n";
  str<<"  \\pagestyle{empty}\n\n";
  str<<"  \\begin{longtable}{ccc}\n\n";
  size_t ng(0);
  std::set<std::string> cvs;
  WriteOutGraph(str,&graphs,ng,cvs);
  str<<"\n  \\end{longtable}\n\n";
  if (!cvs.empty() && (m_pgmode&1)) {
    str<<"  \\begin{longtable}{c}\n";
    for (std::set<std::string>::const_iterator vit(cvs.begin());
	 vit!=cvs.end();++vit) str<<"    $"<<*vit<<"$\\\\\n";
    str<<"  \\end{longtable}\n\n";
  }
  str<<"\\end{fmffile}\n";
  str<<"\\end{document}\n";
}

void Amplitude::PrintStatistics
(std::ostream &str,const int mode) const
{
  if (mode&1) 
    str<<"Amplitude statistics (n="
       <<m_n<<") {\n  level currents vertices\n"<<std::right;
  size_t csum(0), vsum(0);
  for (size_t i(1);i<m_n;++i) {
    csum+=m_cur[i].size();
    size_t cvsum(0);
    for (size_t j(0);j<m_cur[i].size();++j) cvsum+=m_cur[i][j]->NIn();
    if (mode&1)
      str<<"  "<<std::setw(5)<<i<<" "<<std::setw(8)
	 <<m_cur[i].size()<<" "<<std::setw(8)<<cvsum<<"\n";
    vsum+=cvsum;
  }
  if (mode&1) str<<std::left<<"} -> ";
  str<<csum<<" currents, "<<vsum<<" vertices"<<std::endl;
}

