#include "COMIX/Amplitude/Amplitude.H"

#include "MODEL/Interaction_Models/Single_Vertex.H"
#include "PHASIC++/Main/Color_Integrator.H"
#include "PHASIC++/Main/Helicity_Integrator.H"
#include "PHASIC++/Process/Process_Base.H"
#include "METOOLS/Explicit/Dipole_Kinematics.H"
#include "METOOLS/Main/Spin_Structure.H"
#include "ATOOLS/Phys/Spinor.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/STL_Tools.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/My_MPI.H"

using namespace COMIX;
using namespace ATOOLS;

static bool csscite(false);

static const double invsqrttwo(1.0/sqrt(2.0));

Amplitude::Amplitude():
  p_model(NULL), m_nin(0), m_nout(0), m_n(0), m_nf(6), m_ngpl(3),
  m_oew(99), m_oqcd(99), m_maxoew(99), m_maxoqcd(99),
  m_minntc(0), m_maxntc(99),
  m_pmode('D'), p_dinfo(new Dipole_Info()), p_colint(NULL), p_helint(NULL),
  m_trig(true), p_loop(NULL)
{
  p_dinfo->SetMassive(0);
  Data_Reader read(" ",";","!","=");
  std::string prec;
  if (!read.ReadFromFile(prec,"COMIX_PMODE")) prec="D";
  else msg_Tracking()<<METHOD<<"(): Set precision "<<prec<<".\n";
  if (prec!="D") THROW(not_implemented,"Invalid precision mode");
  m_pmode=prec[0];
  int helpi(0);
  if (!read.ReadFromFile(helpi,"COMIX_PG_MODE")) helpi=0;
  else msg_Info()<<METHOD<<"(): Set print graph mode "<<helpi<<".\n";
  m_pgmode=helpi;
  if (!read.ReadFromFile(helpi,"COMIX_VL_MODE")) helpi=0;
  else msg_Info()<<METHOD<<"(): Set vertex label mode "<<helpi<<".\n";
  Vertex::SetVLMode(helpi);
  if (!read.ReadFromFile(helpi,"COMIX_N_GPL")) helpi=3;
  else msg_Info()<<METHOD<<"(): Set graphs per line "<<helpi<<".\n";
  m_ngpl=Max(1,Min(helpi,5));
  double helpd;
  if (!read.ReadFromFile(helpd,"DIPOLE_AMIN")) helpd=Max(rpa->gen.Accu(),1.0e-8);
  else msg_Tracking()<<METHOD<<"(): Set dipole \\alpha_{cut} "<<helpd<<".\n";
  p_dinfo->SetAMin(helpd);
  if (!read.ReadFromFile(helpd,"DIPOLE_ALPHA")) helpd=1.0;
  else msg_Tracking()<<METHOD<<"(): Set dipole \\alpha_{max} "<<helpd<<".\n";
  double amax(helpd);
  if (!read.ReadFromFile(helpd,"DIPOLE_ALPHA_FF")) helpd=amax;
  else msg_Tracking()<<METHOD<<"(): Set FF dipole \\alpha_{max} "<<helpd<<".\n";
  p_dinfo->SetAMax(0,helpd);
  if (!read.ReadFromFile(helpd,"DIPOLE_ALPHA_FI")) helpd=amax;
  else msg_Tracking()<<METHOD<<"(): Set FI dipole \\alpha_{max} "<<helpd<<".\n";
  p_dinfo->SetAMax(2,helpd);
  if (!read.ReadFromFile(helpd,"DIPOLE_ALPHA_IF")) helpd=amax;
  else msg_Tracking()<<METHOD<<"(): Set IF dipole \\alpha_{max} "<<helpd<<".\n";
  p_dinfo->SetAMax(1,helpd);
  if (!read.ReadFromFile(helpd,"DIPOLE_ALPHA_II")) helpd=amax;
  else msg_Tracking()<<METHOD<<"(): Set II dipole \\alpha_{max} "<<helpd<<".\n";
  p_dinfo->SetAMax(3,helpd);
  if (!read.ReadFromFile(helpd,"DIPOLE_KAPPA")) helpd=2.0/3.0;
  else msg_Tracking()<<METHOD<<"(): Set dipole \\kappa="<<helpd<<"\n.";
  p_dinfo->SetKappa(helpd);
  if (!read.ReadFromFile(helpi,"DIPOLE_NF_GSPLIT"))
    helpi=Flavour(kf_jet).Size()/2;
  else msg_Tracking()<<METHOD<<"(): Set dipole N_f="<<helpi<<"\n.";
  p_dinfo->SetNf(helpi);
  if (!read.ReadFromFile(helpd,"DIPOLE_KT2MAX")) helpd=sqr(rpa->gen.Ecms());
  else msg_Tracking()<<METHOD<<"(): Set dipole \\k_{T,max}^2 "<<helpd<<".\n";
  p_dinfo->SetKT2Max(helpd);
  p_dinfo->SetDRMode(0);
  m_sccmur=read.GetValue("USR_WGT_MODE",1);
  m_smth=read.GetValue("NLO_SMEAR_THRESHOLD",0.0);
  m_smpow=read.GetValue("NLO_SMEAR_POWER",4.0);
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
  if (p_dinfo) delete p_dinfo;
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
  for (size_t j(0);j<m_scur.size();++j) delete m_scur[j]; 
  m_cur=Current_Matrix();
  m_scur=Current_Vector();
  m_n=0;
  m_fl=Flavour_Vector();
  m_p=Vec4D_Vector();
  m_ch=Int_Vector();
  m_cl=Int_Matrix();
  m_cress=m_rress=m_ress=Spin_Structure<DComplex>();
  m_cchirs=m_dirs=Int_Vector();
  m_combs.clear();
  m_flavs.clear();
  for (size_t i(0);i<m_subs.size();++i) {
    delete [] m_subs[i]->p_id;
    delete [] m_subs[i]->p_fl;
    delete m_subs[i];
  }
  m_subs.clear();
}

Current *Amplitude::CopyCurrent(Current *const c)
{
  Current_Key ckey(c->Flav(),p_model,c->Id().size());
  Current *cur(Current_Getter::GetObject
	       (std::string(1,m_pmode)+ckey.Type(),ckey));
  if (cur==NULL) return NULL;
  cur->SetDirection(c->Direction());
  cur->SetCut(c->Cut());
  cur->SetOnShell(c->OnShell());
  Int_Vector ids(c->Id()), isfs(ids.size()), pols(ids.size());
  for (size_t i(0);i<ids.size();++i) {
    isfs[i]=m_fl[ids[i]].IsFermion();
    switch (m_fl[ids[i]].IntSpin()) {
    case 0: pols[i]=1; break;
    case 1: pols[i]=2; break;
    case 2: pols[i]=m_fl[ids[i]].IsMassive()?3:2; break;
    default:
      THROW(not_implemented,"Cannot handle spin "+
	    ToString(m_fl[i].Spin())+" particles");
    }
  }
  cur->SetId(ids);
  cur->SetFId(isfs);
  cur->FindPermutations();
  cur->InitPols(pols);
  cur->SetOrderEW(c->OrderEW());
  cur->SetOrderQCD(c->OrderQCD());
  return cur;
}

bool Amplitude::AddRSDipoles()
{
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"(): {\n";
#endif
  Current_Vector scur;
  for (size_t j(0);j<m_n;++j) {
    for (size_t i(0);i<m_cur[2].size();++i) {
      if (m_cur[2][i]->CId()&m_cur[1][j]->CId()) continue;
      // begin temporary
      if (!(m_cur[2][i]->In().front()->OrderQCD()==1 &&
	    m_cur[1][j]->Flav().StrongCharge())) continue;
      // end temporary
      if (m_cur[2][i]->In().size()!=1)
	THROW(not_implemented,"Invalid current");
      if (!AddRSDipole(m_cur[2][i],m_cur[1][j],scur)) return false;
    }
  }
  m_cur[2].insert(m_cur[2].end(),scur.begin(),scur.end());
#ifdef DEBUG__BG
  msg_Debugging()<<"} -> "<<m_scur.size()<<" dipoles\n";
#endif
  if (!csscite) {
    rpa->gen.AddCitation
      (1,"Comix subtraction is published in \\cite{Hoeche:2012xx}.");
    csscite=true;
  }
  return true;
}

bool Amplitude::AddRSDipole
(Current *const c,Current *const s,Current_Vector &scur)
{
#ifdef DEBUG__BG
  msg_Indent();
#endif
  size_t isid((1<<m_nin)-1);
  if ((c->CId()&isid)==isid || c->Flav().IsDummy()) return true;
  Vertex *cin(c->In().front());
  if (cin->JA()->Direction()>0 ||
      cin->JB()->Direction()>0) {
    if (cin->JA()->Flav().Mass() ||
	cin->JB()->Flav().Mass()) return true;
  }
  else {
    double mm(Flavour(p_dinfo->Nf()).Mass());
    if (cin->JA()->Flav().Mass()>mm &&
	cin->JB()->Flav().Mass()>mm) return true;
  }
  for (size_t k(0);k<m_scur.size();++k)
    if (m_scur[k]->CId()==s->CId() &&
	m_scur[k]->Sub()->CId()==c->CId()) return true;
  Current *jijt(CopyCurrent(c)), *jkt(CopyCurrent(s));
  jijt->SetSub(jkt);
  jijt->SetKey(m_scur.size());
  scur.push_back(jijt);
  jkt->SetSub(jijt);
  jkt->SetKey(m_scur.size());
  m_scur.push_back(jkt);
  Vertex_Key svkey(cin->JA(),cin->JB(),NULL,jijt,p_model);
  svkey.m_p=std::string(1,m_pmode);
  svkey.p_dinfo=p_dinfo;
  svkey.p_k=s;
  svkey.p_kt=jkt;
  MODEL::VMIterator_Pair vmp(p_model->GetVertex(svkey.ID()));
  for (MODEL::Vertex_Map::const_iterator vit(vmp.first);
       vit!=vmp.second;++vit) {
    svkey.p_mv=vit->second;
    Vertex *v(new Vertex(svkey));
    v->SetJA(svkey.p_a);
    v->SetJB(svkey.p_b);
    v->SetJC(svkey.p_c);
  }
#ifdef DEBUG__BG
  jijt->Print();
  jkt->Print();
#endif
  return true;
}

bool Amplitude::AddVIDipoles()
{
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"(): {\n";
#endif
  Current_Vector scur;
  for (size_t i(0);i<m_n;++i) {
    for (size_t j(i+1);j<m_n;++j) {
      int sc[2]={m_cur[1][i]->Flav().StrongCharge(),
		 m_cur[1][j]->Flav().StrongCharge()};
      // begin temporary
      if (sc[0]==0 || sc[1]==0) continue;
      // end temporary
      if (!AddVIDipole(m_cur[1][i],m_cur[1][j],scur)) return false;
    }
  }
  m_cur[1].insert(m_cur[1].end(),scur.begin(),scur.end());
#ifdef DEBUG__BG
  msg_Debugging()<<"} -> "<<m_scur.size()<<" dipoles\n";
#endif
  if (!csscite) {
    rpa->gen.AddCitation
      (1,"Comix subtraction is published in \\cite{Hoeche:2012xx}.");
    csscite=true;
  }
  return true;
}

bool Amplitude::AddVIDipole
(Current *const c,Current *const s,Current_Vector &scur)
{
#ifdef DEBUG__BG
  msg_Indent();
#endif
  size_t isid((1<<m_nin)-1);
  if ((c->CId()&isid)==isid || c->Flav().IsDummy()) return true;
  Current *jijt(CopyCurrent(c)), *jkt(CopyCurrent(s));
  jijt->SetSub(jkt);
  jijt->SetKey(m_scur.size());
  scur.push_back(jijt);
  jkt->SetSub(jijt);
  jkt->SetKey(m_scur.size());
  m_scur.push_back(jkt);
  Vertex_Key svkey(c,NULL,NULL,jijt,p_model);
  svkey.m_p=std::string(1,m_pmode);
  svkey.p_dinfo=p_dinfo;
  svkey.p_k=s;
  svkey.p_kt=jkt;
  MODEL::VMIterator_Pair vmp(p_model->GetVertex(svkey.ID(1)));
  for (MODEL::Vertex_Map::const_iterator vit(vmp.first);
       vit!=vmp.second;++vit) {
    svkey.p_mv=vit->second;
    Vertex *v(new Vertex(svkey));
    v->SetJA(svkey.p_a);
    v->SetJB(svkey.p_b);
    v->SetJC(svkey.p_c);
  }
  if (svkey.p_mv==NULL) {
    std::swap<Current*>(svkey.p_a,svkey.p_b);
    vmp=p_model->GetVertex(svkey.ID(0));
    for (MODEL::Vertex_Map::const_iterator vit(vmp.first);
	 vit!=vmp.second;++vit) {
      svkey.p_mv=vit->second;
      Vertex *v(new Vertex(svkey));
      v->SetJA(svkey.p_a);
      v->SetJB(svkey.p_b);
      v->SetJC(svkey.p_c);
    }
  }
#ifdef DEBUG__BG
  jijt->Print();
  jkt->Print();
#endif
  return true;
}

bool Amplitude::MatchIndices
(const Int_Vector &ids,const Vertex_Key &vkey) const
{
  size_t n(ids.size());
  for (size_t o(0);o<n;++o) {
    bool found(false);
    for (size_t p(0);p<vkey.p_a->Id().size();++p) 
      if (vkey.p_a->Id()[p]==ids[o]) {
	if (found) return false;
	found=true;
      }
    for (size_t p(0);p<vkey.p_b->Id().size();++p) 
      if (vkey.p_b->Id()[p]==ids[o]) {
	if (found) return false;
	found=true;
      }
    if (vkey.p_e)
      for (size_t p(0);p<vkey.p_e->Id().size();++p) 
	if (vkey.p_e->Id()[p]==ids[o]) {
	  if (found) return false;
	  found=true;
	}
    if (!found) return false;
  }
  return true;
}

int Amplitude::CheckDecay(const ATOOLS::Flavour &fl,
			       const Int_Vector &ids) const
{
  size_t cid(0);
  if (m_decid.empty() && m_ndc.empty()) return 0;
  for (size_t i(0);i<ids.size();++i) cid+=1<<ids[i];
  if ((cid&(1<<m_nin)-1)==0)
    for (size_t i(0);i<m_ndc.size();++i)
      if (m_ndc[i].Includes(fl)) return -1;
  if (m_decid.empty()) return 0;
  for (size_t i(0);i<m_decid.size();++i) {
    size_t did(m_decid[i]->m_id);
    Flavour dfl(m_decid[i]->m_fl);
    if (did&(1<<0)) {
      did=(1<<m_n)-1-did;
      dfl=dfl.Bar();
    }
    if (did==cid) {
      if (dfl.Includes(fl)) return i+1;
      return -1;
    }
    if (!((did&cid)==0 || (did&cid)==cid || (did&cid)==did)) {
      return -1;
    }
  }
  return 0;
}

Vertex *Amplitude::AddCurrent
(const Current_Key &ckey,Vertex_Key &vkey,const size_t &n,
 const int dec,size_t &oewmax,size_t &oqcdmax,
 std::map<std::string,Current*> &curs)
{
  if (vkey.p_mv==NULL || !vkey.p_mv->on || 
      vkey.p_mv->dec<0) return NULL;
  vkey.m_p=std::string(1,m_pmode);
  Vertex *v(new Vertex(vkey));
  size_t oew(vkey.p_a->OrderEW()+
	     vkey.p_b->OrderEW()+v->OrderEW());
  size_t oqcd(vkey.p_a->OrderQCD()+
	      vkey.p_b->OrderQCD()+v->OrderQCD());
  size_t ntc(vkey.p_a->NTChannel()+vkey.p_b->NTChannel());
  bool isa((vkey.p_a->CId()&1)^(vkey.p_a->CId()&2));
  bool isb((vkey.p_b->CId()&1)^(vkey.p_b->CId()&2));
  if (n<m_n-1 && vkey.p_e==NULL) {
    ntc+=(isa^isb)&&!ckey.m_fl.Strong();
  }
  if (vkey.p_e) {
    oew+=vkey.p_e->OrderEW();
    oqcd+=vkey.p_e->OrderQCD();
    ntc+=vkey.p_e->NTChannel();
    if (n<m_n-1) {
      bool ise((vkey.p_e->CId()&1)^(vkey.p_e->CId()&2));
      ntc+=(isa^isb^ise)&&!ckey.m_fl.Strong();
    }
  }
  if (!v->Active() || oew>m_oew || oqcd>m_oqcd ||
      oew>m_maxoew || oqcd>m_maxoqcd ||
      (n==m_n-1 && ((m_oew<99 && oew!=m_oew) ||
		    (m_oqcd<99 && oqcd!=m_oqcd) ||
		    ntc<m_minntc || ntc>m_maxntc))) {
#ifdef DEBUG__BG
    msg_Debugging()<<"delete vertex {"<<vkey.p_a->Flav()<<",("
		   <<vkey.p_a->OrderEW()<<","<<vkey.p_a->OrderQCD()
		   <<")}{"<<vkey.p_b->Flav()<<",("
		   <<vkey.p_b->OrderEW()<<","<<vkey.p_b->OrderQCD()
		   <<")}-"<<vkey.Type()<<"("<<v->OrderEW()<<","
		   <<v->OrderQCD()<<")->{"<<vkey.p_c->Flav()<<"} => ("
		   <<oew<<","<<oqcd<<") vs. max = ("<<m_oew<<","
		   <<m_oqcd<<"), act = "<<v->Active()<<", n t-ch = "
		   <<ntc<<" vs "<<m_minntc<<"/"<<m_maxntc<<"\n";
#endif
    delete v;
    return NULL;
  }
  Current *sub(NULL);
  if (p_dinfo->Mode()) {
    size_t cid(vkey.p_a->CId()|vkey.p_b->CId());
    if (vkey.p_a->Sub()) {
      if (vkey.p_a->Id().size()==1 &&
	  p_dinfo->Mode()==1) return NULL;
      sub=vkey.p_a->Sub();
    }
    if (vkey.p_b->Sub()) {
      if ((vkey.p_b->Id().size()==1 &&
	   p_dinfo->Mode()==1) || sub) return NULL;
      sub=vkey.p_b->Sub();
    }
    if (vkey.p_e) {
      cid|=vkey.p_e->CId();
      if (vkey.p_e->Sub()) {
	if ((vkey.p_e->Id().size()==1 &&
	     p_dinfo->Mode()==1) || sub) return NULL;
	sub=vkey.p_e->Sub();
      }
    }
    if (sub && (sub->CId()&cid)) {
#ifdef DEBUG__BG
      msg_Debugging()<<"delete vertex {"<<sub->Id()
		     <<","<<sub->Sub()->Id()
		     <<"}<->"<<ID(cid)<<"\n";
#endif
      delete v;
      return NULL;
    }
  }
  std::string okey
    ("("+(n<m_n-1?ToString(oew)+","
	  +ToString(oqcd)+";"+ToString(ntc):"")
     +(sub?",S"+ToString(sub->Sub()->Id())+
       ToString(sub->Sub()->Sub()->Id()):"")+")");
  if (oew!=vkey.p_c->OrderEW() || oqcd!=vkey.p_c->OrderQCD() ||
      ntc!=vkey.p_c->NTChannel() || sub!=vkey.p_c->Sub()) {
    std::map<std::string,Current*>::iterator 
      cit(curs.find(okey));
    if (cit!=curs.end()) vkey.p_c=cit->second;
    else {
      if (vkey.p_c->OrderEW()>0 || vkey.p_c->OrderQCD()>0 ||
	  vkey.p_c->NTChannel()>0 || vkey.p_c->Sub()!=sub)
	vkey.p_c=Current_Getter::GetObject
	  (std::string(1,m_pmode)+ckey.Type(),ckey);
      if (n<m_n-1) {
	if (dec!=0) {
	  vkey.p_c->SetCut(m_decid[dec-1]->m_nmax);
	  vkey.p_c->SetOnShell(m_decid[dec-1]->m_osd);
	}
	vkey.p_c->SetOrderEW(oew);
	vkey.p_c->SetOrderQCD(oqcd);
	vkey.p_c->SetNTChannel(ntc);
      }
      else {
	oewmax=Max(oewmax,oew); 
	oqcdmax=Max(oqcdmax,oqcd);
      }
      vkey.p_c->SetSub(sub);
      curs[okey]=vkey.p_c;
    }
  }
  v->SetJA(vkey.p_a);
  v->SetJB(vkey.p_b);
  v->SetJE(vkey.p_e);
  v->SetJC(vkey.p_c);
  return v;
}

void Amplitude::AddCurrent(const Int_Vector &ids,const size_t &n,
			   const Flavour &fl,const int dir)
{
  // add new currents
  size_t oewmax(0), oqcdmax(0);
  int dec(CheckDecay(fl,ids));
  if (dec<0) return;
  std::map<std::string,Current*> curs;
  Current_Key ckey(dir>0?fl.Bar():fl,p_model,ids.size());
  Current *cur(Current_Getter::GetObject
	       (std::string(1,m_pmode)+ckey.Type(),ckey));
  if (cur==NULL) return;
  cur->SetDirection(dir);
  if (dec!=0) {
    cur->SetCut(m_decid[dec-1]->m_nmax);
    cur->SetOnShell(m_decid[dec-1]->m_osd);
  }
  std::set<Vertex_Key> v3;
  // compose current from all possible subcurrents
  for (size_t i(1);i<n;++i) {
    for (size_t j(0);j<m_cur[i].size();++j) {
      for (size_t k(0);k<m_cur[n-i].size();++k) {
	Vertex_Key vkey(m_cur[i][j],m_cur[n-i][k],NULL,cur,p_model);
	if (!MatchIndices(ids,vkey) ||
	    v3.find(vkey.SwapAB())!=v3.end()) continue;
	MODEL::VMIterator_Pair vmp(p_model->GetVertex(vkey.ID()));
	for (MODEL::Vertex_Map::const_iterator vit(vmp.first);
	     vit!=vmp.second;++vit) {
	  vkey.p_mv=vit->second;
	  if (AddCurrent(ckey,vkey,n,dec,oewmax,oqcdmax,curs))
	    v3.insert(Vertex_Key(vkey.p_a,vkey.p_b,NULL,cur,p_model));
	}
      }
    }
  }
  if (p_model->VInfo()&1)
    for (size_t i(1);i<n-1;++i) {
      for (size_t j(1);j<n-i;++j) {
	for (size_t k(0);k<m_cur[i].size();++k) {
	  for (size_t l(0);l<m_cur[j].size();++l) {
	    for (size_t m(0);m<m_cur[n-i-j].size();++m) {
	      Vertex_Key vkey(m_cur[i][k],m_cur[j][l],
			      m_cur[n-i-j][m],cur,p_model);
	      if (!MatchIndices(ids,vkey) ||
		  v3.find(vkey.SwapAB())!=v3.end() ||
		  v3.find(vkey.SwapBE())!=v3.end() ||
		  v3.find(vkey.SwapEA())!=v3.end() ||
		  v3.find(vkey.SwapEA().SwapAB())!=v3.end() ||
		  v3.find(vkey.SwapEA().SwapBE())!=v3.end()) continue;
	      MODEL::VMIterator_Pair vmp(p_model->GetVertex(vkey.ID()));
	      for (MODEL::Vertex_Map::const_iterator vit(vmp.first);
		   vit!=vmp.second;++vit) {
		vkey.p_mv=vit->second;
		if (AddCurrent(ckey,vkey,n,dec,oewmax,oqcdmax,curs))
		  v3.insert(Vertex_Key(vkey.p_a,vkey.p_b,
				       vkey.p_e,cur,p_model));
	      }
	    }
	  }
	}
      }
    }
  if (v3.empty() && n>1) {
    delete cur;
    return;
  }
  if (n==1) curs[""]=cur;
  if (n==m_n-1) {
    m_maxoew=oewmax;
    m_maxoqcd=oqcdmax;
  }
  Int_Vector isfs(ids.size()), pols(ids.size());
  for (size_t i(0);i<ids.size();++i) {
    isfs[i]=m_fl[ids[i]].IsFermion();
    switch (m_fl[ids[i]].IntSpin()) {
    case 0: pols[i]=1; break;
    case 1: pols[i]=2; break;
    case 2: pols[i]=m_fl[ids[i]].IsMassive()?3:2; break;
    default:
      THROW(not_implemented,"Cannot handle spin "+
	    ToString(m_fl[i].Spin())+" particles");
    }
  }
  for (std::map<std::string,Current*>::iterator 
	 cit(curs.begin());cit!=curs.end();++cit) {
    cit->second->SetId(ids);
    cit->second->SetFId(isfs);
    cit->second->FindPermutations();
    cit->second->InitPols(pols);
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
      if (p_dinfo->Mode()==0) {
	if (!m_fl.front().IsOn()) return false;
	AddCurrent(ids,n,m_fl.front().Bar(),m_dirs.front());
      }
      else {
	size_t kid(0);
	for (size_t i(0);i<ids.size();++i) {
	  if (ids[i]>(int)kid) break;
	  else ++kid;
	}
	if (kid>=m_cur[1].size()) THROW(fatal_error,"Internal error");
	if (!m_fl[kid].IsOn()) return false;
	AddCurrent(ids,n,m_fl[kid].Bar(),m_dirs[kid]);
	if (kid==0 && (m_cur.back().size()==0 || 
		       m_cur.back().back()->CId()&1)) return false;
      }
    }
    else {
      if (m_affm.size()) {
	for (size_t i(0);i<m_affm[n].size();++i) {
	  Flavour cfl((kf_code)abs(m_affm[n][i]),m_affm[n][i]<0);
	  AddCurrent(ids,n,cfl,0);
	}
      }
      else {
	for (size_t i(0);i<fls.size();++i)
	  AddCurrent(ids,n,fls[i],0);
      }
    }
    return true;
  }
  size_t last(ids.empty()?0:ids.back());
  ids.push_back(0);
  if (p_dinfo->Mode()==0 && n==m_n-1) {
    ids.back()=last+1;
    if (!Construct(fls,ids,n)) return false;
    return m_cur.back().size();
  }
  int first(p_dinfo->Mode()&&ids.size()==1?last:last+1);
  for (size_t i(first);i<m_n;++i) {
    ids.back()=i;
    if (!Construct(fls,ids,n) && (n==1 || n==m_n-1)) return false;
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
  if ((p_dinfo->Mode()&2) && !AddVIDipoles()) return false;
  if (!Construct(fls,ids,2)) return false;
  if (p_dinfo->Mode()==1 && !AddRSDipoles()) return false;
  for (size_t i(3);i<m_n;++i)
    if (!Construct(fls,ids,i)) return false;
  if (p_dinfo->Mode()) {
    for (Current_Vector::iterator cit(m_cur.back().begin());
	 cit!=m_cur.back().end();++cit)
      if ((*cit)->Sub()==NULL && (*cit)->CId()&1) {
#ifdef DEBUG__BG
	msg_Debugging()<<"delete obsolete current "<<*cit<<"\n";
	(*cit)->Print();
#endif
	delete *cit;
	cit=--m_cur.back().erase(cit);
      }
  }
  for (size_t j(m_n-2);j>1;--j)
    for (Current_Vector::iterator cit(m_cur[j].begin());
	 cit!=m_cur[j].end();++cit)
      if ((*cit)->Dangling()) {
#ifdef DEBUG__BG
	msg_Debugging()<<"delete current "<<**cit<<", "<<(*cit)->Dangling()
		       <<", O("<<(*cit)->OrderEW()<<","<<(*cit)->OrderQCD()
		       <<") vs. O("<<m_oew<<","<<m_oqcd<<") / O_{max}("
		       <<m_maxoew<<","<<m_maxoqcd<<"), n t-ch = "
		       <<(*cit)->NTChannel()<<" vs. "
		       <<m_minntc<<"/"<<m_maxntc<<"\n";
#endif
	if ((*cit)->Sub() && (*cit)->Sub()->Sub()==*cit)
	  (*cit)->Sub()->SetSub(NULL);
	delete *cit;
	cit=--m_cur[j].erase(cit);
      }
  if (p_dinfo->Mode()) {
    for (Current_Vector::iterator cit(m_cur.back().begin());
	 cit!=m_cur.back().end();++cit)
      if ((*cit)->Sub()==NULL) {
	Current *c(*cit);
	m_cur.back().erase(cit);
	m_cur.back().insert(m_cur.back().begin(),c);
	break;
      }
    for (Current_Vector::iterator cit(m_scur.begin());
	 cit!=m_scur.end();++cit)
      if ((*cit)->Sub()==NULL) {
#ifdef DEBUG__BG
	msg_Debugging()<<"delete dipole current "<<**cit<<"\n";
#endif
	delete *cit;
	cit=--m_scur.erase(cit);
      }
  }
  FillCombinations();
  msg_Debugging()<<METHOD<<"(): Amplitude statistics (n="
		 <<m_n<<") {\n  level currents vertices\n"<<std::right;
  size_t csum(0), vsum(0), scsum(0), svsum(0);
  for (size_t i(1);i<m_n;++i) {
    size_t cvsum(0), csvsum(0);
    for (size_t j(0);j<m_cur[i].size();++j) {
      if (m_cur[i][j]->Sub()) {
	++scsum;
	csvsum+=m_cur[i][j]->NIn();
      }
      else {
	++csum;
	cvsum+=m_cur[i][j]->NIn();
      }
    }
    msg_Debugging()<<"  "<<std::setw(5)<<i<<" "<<std::setw(8)
		   <<m_cur[i].size()<<" "<<std::setw(8)<<cvsum<<"\n";
    vsum+=cvsum;
    svsum+=csvsum;
  }
  msg_Debugging()<<std::left<<"} -> "<<csum<<"(+"<<scsum<<") currents, "
		 <<vsum<<"(+"<<svsum<<") vertices"<<std::endl;
  return true;
}

bool Amplitude::ReadInAmpFile(const std::string &name)
{
  std::string ampfile(rpa->gen.Variable("SHERPA_CPP_PATH")
		      +"/Process/Comix/"+name+".map");
  My_In_File amp(ampfile);
  if (!amp.Open()) return false;
  std::string cname, cmname;
  *amp>>cname>>cmname;
  if (cname!=name || cmname!=name || amp->eof())
    THROW(fatal_error,"Corrupted map file '"+ampfile+"'");
  m_affm.resize(m_n);
  for (size_t size, i(2);i<m_n;++i) {
    *amp>>size;
    m_affm[i].resize(size);
    for (size_t j(0);j<size;++j) *amp>>m_affm[i][j];
  }
  *amp>>cname;
  if (cname!="eof")
    THROW(fatal_error,"Corrupted map file '"+ampfile+"'");
  return true;
}

void Amplitude::WriteOutAmpFile(const std::string &name)
{
#ifdef USING__MPI
  if (MPI::COMM_WORLD.Get_rank()) return;
#endif
  std::string ampfile(rpa->gen.Variable("SHERPA_CPP_PATH")
		      +"/Process/Comix/"+name+".map");
  if (FileExists(ampfile)) return;
  std::ofstream amp(ampfile.c_str());
  if (!amp.good()) return;
  amp<<name<<" "<<name<<"\n";
  m_affm.resize(m_n);
  for (size_t i(2);i<m_n;++i) {
    std::set<long int> kfcs;
    m_affm[i].clear();
    for (size_t j(0);j<m_cur[i].size();++j) {
      long int kfc(m_cur[i][j]->Flav());
      if (kfcs.find(kfc)==kfcs.end()) m_affm[i].push_back(kfc);
      kfcs.insert(kfc);
    }
    amp<<m_affm[i].size();
    for (size_t j(0);j<m_affm[i].size();++j) amp<<" "<<m_affm[i][j];
    amp<<"\n";
  }
  amp<<"eof\n";
}

bool Amplitude::Initialize
(const size_t &nin,const size_t &nout,const std::vector<Flavour> &flavs,
 const double &isf,const double &fsf,MODEL::Model_Base *const model,
 MODEL::Coupling_Map *const cpls,const int smode,const size_t &oew,
 const size_t &oqcd,const size_t &maxoew,const size_t &maxoqcd,
 const size_t &minntc,const size_t &maxntc,const std::string &name)
{
  CleanUp();
  m_nin=nin;
  m_nout=nout;
  m_n=m_nin+m_nout;
  m_oew=oew;
  m_oqcd=oqcd;
  m_minntc=minntc;
  m_maxntc=maxntc;
  m_maxoew=Min(oew,maxoew);
  m_maxoqcd=Min(oqcd,maxoqcd);
  p_dinfo->SetMode(smode);
  ReadInAmpFile(name);
  Int_Vector incs(m_nin,1);
  incs.resize(flavs.size(),-1);
  if (!Construct(incs,flavs,model,cpls)) return false;
  std::map<Flavour,size_t> fc;
  for (size_t i(nin);i<flavs.size();++i) {
    std::map<Flavour,size_t>::iterator fit(fc.find(flavs[i]));
    if (fit==fc.end()) {
      fc[flavs[i]]=0;
      fit=fc.find(flavs[i]);
    }
    ++fit->second;
  }
  m_fsf=fsf;
  m_sf=m_fsf*isf;
  return true;
}

void Amplitude::ConstructNLOEvents()
{
  msg_Debugging()<<METHOD<<"(): {\n";
  for (size_t i(0);i<m_scur.size();++i) {
    Dipole_Kinematics *kin(m_scur[i]->Sub()->In().front()->Kin());
    NLO_subevt *sub(new NLO_subevt());
    m_subs.push_back(sub);
    sub->m_idx=i;
    sub->m_n=m_nin+m_nout-1;
    Flavour *fls(new Flavour[sub->m_n]);
    size_t *ids(new size_t[sub->m_n]);
    sub->m_i=kin->JI()->Id().front();
    sub->m_j=kin->JJ()->Id().front();
    sub->m_k=kin->JK()->Id().front();
    Current_Vector cur(sub->m_n,NULL);
    for (size_t k(0), j(0);k<m_nin+m_nout;++k) {
      if (k==sub->m_j) continue;
      if (k==sub->m_i) {
	cur[j]=kin->JIJT();
	sub->m_ijt=j;
      }
      else if (k==sub->m_k) {
	cur[j]=kin->JKT();
	sub->m_kt=j;
      }
      else {
	cur[j]=m_cur[1][k];
      }
      fls[j]=cur[j]->Flav();
      ids[j]=cur[j]->CId();
      ++j;
    }
    kin->SetCurrents(cur);
    sub->p_mom=&kin->Momenta().front();
    for (size_t j(0);j<m_nin;++j) fls[j]=fls[j].Bar();
    sub->p_fl=fls;
    sub->p_id=ids;
    sub->p_dec=&m_decid;
    PHASIC::Process_Info cpi;
    for (size_t j(0);j<m_nin;++j)
      cpi.m_ii.m_ps.push_back(PHASIC::Subprocess_Info(fls[j]));
    for (size_t j(m_nin);j<m_nin+m_nout-1;++j)
      cpi.m_fi.m_ps.push_back(PHASIC::Subprocess_Info(fls[j]));
    PHASIC::Process_Base::SortFlavours(cpi);
    sub->m_pname=PHASIC::Process_Base::GenerateName(cpi.m_ii,cpi.m_fi);
    msg_Indent();
    msg_Debugging()<<*sub<<"\n";
  }
  NLO_subevt *sub(new NLO_subevt());
  m_subs.push_back(sub);
  sub->m_n=m_nin+m_nout;
  Flavour *fls(new Flavour[sub->m_n]);
  size_t *ids(new size_t[sub->m_n]);
  for (size_t i(0);i<m_nin+m_nout;++i) {
    fls[i]=m_fl[i];
    ids[i]=1<<i;
  }
  sub->p_mom=&m_p.front();
  sub->p_fl=fls;
  sub->p_id=ids;
  sub->p_dec=&m_decid;
  PHASIC::Process_Info cpi;
  for (size_t j(0);j<m_nin;++j)
    cpi.m_ii.m_ps.push_back(PHASIC::Subprocess_Info(fls[j]));
  for (size_t j(m_nin);j<m_nin+m_nout;++j)
    cpi.m_fi.m_ps.push_back(PHASIC::Subprocess_Info(fls[j]));
  PHASIC::Process_Base::SortFlavours(cpi);
  sub->m_pname=PHASIC::Process_Base::GenerateName(cpi.m_ii,cpi.m_fi);
  sub->m_i=sub->m_j=sub->m_k=0;
  {
    msg_Indent();
    msg_Debugging()<<*sub<<"\n";
  }
  for (size_t i(0);i<m_subs.size();++i) m_subs[i]->p_real=sub;
  m_sid.resize(1,0);
  for (size_t i(1);i<m_cur.back().size();++i)
    for (size_t j(0);j<m_subs.size();++j)
      if (m_cur.back()[i]->Sub()==m_scur[j]) {
	m_sid.push_back(j);
	break;
      }
  msg_Debugging()<<"}\n";
}

void Amplitude::ConstructDSijMap()
{
  m_dsm.clear();
  m_dsf.clear();
  size_t c(0), ifp(0);
  std::vector<int> plist(m_fl.size());
  for (size_t i(0);i<m_fl.size();++i)
    plist[i]=m_fl[i].Strong()?c++:0;
  m_dsij.resize(c,std::vector<double>(c));
  Flavour ifl(m_cur[1][0]->Flav());
  if (ifl.IsFermion() && !ifl.IsAnti()) {
    for (int i(m_fl.size()-1);i>0;--i)
      if (m_cur[1][i]->Flav().IsFermion()) ++ifp;
  }
  m_dsf.resize(m_cur.back().size(),ifp%2?-1.0:1.0);
  for (size_t j(1);j<m_cur.back().size();++j) {
    m_dsm[j]=std::pair<int,int>
      (plist[m_cur.back()[j]->Sub()->Sub()->Id().front()],
       plist[m_cur.back()[j]->Sub()->Id().front()]);
    Flavour ffl(m_cur.back()[j]->Sub()->Flav());
    if (!ffl.IsFermion()) continue;
    int idx(m_cur.back()[j]->Sub()->Id().front());
    if (ffl.IsAnti()) {
      for (int i(0);i<idx;++i)
	if (m_cur[1][i]->Flav().IsFermion()) m_dsf[j]=-m_dsf[j];
    }
    else {
      for (int i(m_fl.size()-1);i>idx;--i)
	if (m_cur[1][i]->Flav().IsFermion()) m_dsf[j]=-m_dsf[j];
    }
  }
}

void Amplitude::FillCombinations()
{
  msg_Debugging()<<METHOD<<"(): {\n";
  msg_Debugging()<<"  flavours {\n";
  for (size_t i(2);i<m_n-1;++i)
    for (size_t j(0);j<m_cur[i].size();++j) {
      if (m_cur[i][j]->Sub() || 
	  m_cur[i][j]->Flav().IsDummy()) continue;
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
      if (m_cur[i][j]->Sub() ||
	  m_cur[i][j]->Flav().IsDummy()) continue;
      Vertex_Vector ins(m_cur[i][j]->In());
      for (size_t k(0);k<ins.size();++k) {
	if (ins[k]->JE() ||
	    ins[k]->JA()->Flav().IsDummy() ||
	    ins[k]->JA()->Flav().IsDummy()) continue;
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
  size_t svlmode(Vertex::VLMode());
  Vertex::SetVLMode(7);
  for (size_t n(1);n<m_n;++n) {
    if (ampl.m_cur[n].size()!=m_cur[n].size()) {
      msg_Debugging()<<"  current count differs\n} no match\n";
      Vertex::SetVLMode(svlmode);
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
          Vertex::SetVLMode(svlmode);
	  flmap.clear();
	  return false;
	}
	if (ampl.m_cur[n][i]->Flav().Mass()!=
	    m_cur[n][i]->Flav().Mass() ||
	    ampl.m_cur[n][i]->Flav().Width()!=
	    m_cur[n][i]->Flav().Width()) {
	  msg_Debugging()<<"    mass or width differs\n  }\n";
	  msg_Debugging()<<"} no match\n";
          Vertex::SetVLMode(svlmode);
	  flmap.clear();
	  return false;
	}
	if (ampl.m_cur[n][i]->Flav().StrongCharge()!=
	    m_cur[n][i]->Flav().StrongCharge()) {
	  msg_Debugging()<<"    color structure differs\n  }\n";
	  msg_Debugging()<<"} no match\n";
          Vertex::SetVLMode(svlmode);
	  flmap.clear();
	  return false;
	}
	flmap[ampl.m_cur[n][i]->Flav()]=m_cur[n][i]->Flav();
      }
      else {
	if (m_cur[n][i]->Flav()!=
	    flmap[ampl.m_cur[n][i]->Flav()]) {
	  msg_Debugging()<<"    current differs\n  }\n";
	  msg_Debugging()<<"} no match\n";
	  Vertex::SetVLMode(svlmode);
	  flmap.clear();
	  return false;
	}
      }
      Vertex_Vector vin(m_cur[n][i]->In());
      Vertex_Vector avin(ampl.m_cur[n][i]->In());
      if (avin.size()!=vin.size()) {
	msg_Debugging()<<"    vertex count differs\n  }\n";
	msg_Debugging()<<"} no match\n";
	Vertex::SetVLMode(svlmode);
	flmap.clear();
	return false;
      }
      if (n>1)
      for (size_t j(0);j<vin.size();++j) {
#ifdef DEBUG__BG
	msg_Debugging()<<"    check m_in["<<j<<"] {\n";
#endif
	if (!avin[j]->Map(*vin[j])) {
	  msg_Debugging()<<"    } no match\n  }\n} no match\n";
	  Vertex::SetVLMode(svlmode);
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
  Vertex::SetVLMode(svlmode);
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
  SetCouplings();
  for (size_t i(0);i<m_n;++i) 
    m_cur[1][i]->ConstructJ(m_p[i],m_ch[i],m_cl[i][0],m_cl[i][1]);
  for (size_t i(m_n);i<m_cur[1].size();++i) m_cur[1][i]->Evaluate();
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

void Amplitude::SetCouplings() const
{
#ifdef DEBUG__CF
  msg_Debugging()<<METHOD<<"(): {\n";
#endif
  for (size_t i(0);i<m_cpls.size();++i) {
    double fac(1.0);
    Vertex *v(m_cpls[i].p_v);
    size_t oqcd(m_cpls[i].m_oqcd), oew(m_cpls[i].m_oew);
    MODEL::Coupling_Data *aqcd(m_cpls[i].p_aqcd);
    MODEL::Coupling_Data *aqed(m_cpls[i].p_aqed);
    if (aqcd && oqcd) {
#ifdef DEBUG__CF
      msg_Debugging()<<"  qcd: "<<sqrt(aqcd->Factor())<<" ^ "<<oqcd
		     <<" = "<<pow(aqcd->Factor(),oqcd/2.0)<<"\n";
#endif
      fac*=pow(aqcd->Factor(),oqcd/2.0);
    }
    if (aqed && oew) {
#ifdef DEBUG__CF
      msg_Debugging()<<"  qed: "<<sqrt(aqed->Factor())<<" ^ "<<oew
		     <<" = "<<pow(aqed->Factor(),oew/2.0)<<"\n";
#endif
      fac*=pow(aqed->Factor(),oew/2.0);
    }
    v->SetCplFac(fac);
  }
#ifdef DEBUG__CF
  msg_Debugging()<<"}\n";
#endif
}

void Amplitude::ResetZero()
{
  for (size_t n(m_n-2);n>=2;--n) {
    for (size_t i(0);i<m_cur[n].size();++i) 
      m_cur[n][i]->ResetZero();
  }
}

bool Amplitude::SetMomenta(const Vec4D_Vector &moms)
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
  if (m_subs.empty()) return true;
  p_dinfo->SetStat(1);
  for (size_t i(0);i<m_cur[1].size();++i) m_cur[1][i]->SetP(m_p[i]);
  for (size_t i(0);i<m_scur.size();++i)
    m_scur[i]->Sub()->In().front()->Kin()->Evaluate();
  return p_dinfo->Stat();
}

bool Amplitude::JetTrigger
(PHASIC::Combined_Selector *const sel,const int mode)
{
  if (m_subs.empty() || sel==NULL) return true;
  NLO_subevtlist tmp;
  tmp.resize(1,m_subs.back());
  bool trig(m_trig=sel->JetTrigger(m_p,&tmp));
  for (size_t i(0);i<m_scur.size();++i) {
    tmp.back()=m_subs[i];
    Dipole_Kinematics *kin(m_scur[i]->Sub()->In().front()->Kin());
    int ltrig(sel->JetTrigger(kin->Momenta(),&tmp));
    kin->SetF(1.0);
    if (m_smth) {
      double a(m_smth>0.0?kin->KT2():kin->Y());
      if (a>0.0 && a<dabs(m_smth)) {
	kin->SetF(pow(a/dabs(m_smth),m_smpow));
	if (ltrig==0) kin->SetF(-kin->F());
	ltrig=1;
      }
    }
    kin->AddTrig(ltrig);
    if (mode==1) {
      int da(kin->KT2()<m_subs[i]->m_mu2[stp::res]);
      int ds(kin->Y()<p_dinfo->AMax(kin->Type()));
      kin->SetTrig(ltrig&abs(ds-da));
    }
    else if (mode==2) {
      int da(kin->KT2()<m_subs[i]->m_mu2[stp::res]);
      kin->SetTrig(ltrig&da);
    }
    trig|=kin->Trig();
  }
  return trig;
}

double Amplitude::KT2Trigger(NLO_subevt *const sub,const int mode)
{
  if (mode==0) return 1.0;
  Dipole_Kinematics *kin
    (m_scur[sub->m_idx]->Sub()->In().front()->Kin());
  if (mode==1) {
    int da(kin->KT2()<sub->m_mu2[stp::res]);
    int ds(kin->Y()<p_dinfo->AMax(kin->Type()));
    kin->SetTrig(abs(ds-da));
    return ds-da;
  }
  if (mode==2) {
    int da(kin->KT2()<sub->m_mu2[stp::res]);
    kin->SetTrig(da);
    return da;
  }
  THROW(not_implemented,"Invalid call");
  return 0.0;
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
    Color_Calculator::SetCIMin(1);
    Color_Calculator::SetCIMax(0);
  }
  else {
    Color_Calculator::SetCIMin(1);
    Color_Calculator::SetCIMax(3);
  }
}

double Amplitude::EpsSchemeFactor(const Vec4D_Vector &mom) const
{
  if (p_loop) return p_loop->Eps_Scheme_Factor(mom);
  return 4.0*M_PI;
}

bool Amplitude::Evaluate(const Int_Vector &chirs)
{
  THROW(not_implemented,"Helicity sampling currently disabled");
  return true;
}

bool Amplitude::EvaluateAll()
{
  if (p_loop) p_dinfo->SetDRMode(p_loop->DRMode());
  for (size_t i(0);i<m_subs.size();++i) m_subs[i]->Reset();
  for (size_t j(0);j<m_n;++j) m_ch[j]=0;
  MODEL::Coupling_Data *cpl(m_cpls.front().p_aqcd);
  double mu2(cpl?cpl->Scale():-1.0);
  p_dinfo->SetMu2(mu2);
  CalcJL();
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"(): "<<m_ress.size()<<" amplitudes {\n";
#endif
  if (m_pmode=='D') {
    DComplex_Vector ress(m_ress.size(),DComplex(0.0));
    m_cur.back()[0]->Contract(*m_cur[1].front(),m_cchirs,ress);
    for (size_t i(0);i<m_ress.size();++i) m_ress[i]=ress[i];
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
    msg_Debugging()<<"A["<<j<<"]"<<m_ress(j)
		   <<" = "<<m_ress[j]<<" -> "<<std::abs(m_ress[j])<<"\n";
#endif
    csum+=(m_ress[j]*std::conj(m_ress[j])).real();
  }
  m_res=csum/m_sf;
  m_cmur[1]=m_cmur[0]=csum=0.0;
  if (p_dinfo->Mode()) {
    double asf(cpl->Default()*cpl->Factor()/(2.0*M_PI));
    if (p_loop) {
      p_loop->SetRenScale(mu2);
      m_p[0]=-m_p[0];
      m_p[1]=-m_p[1];
      p_loop->Calc(m_p);
      m_p[0]=-m_p[0];
      m_p[1]=-m_p[1];
    }
    if (p_dinfo->Mode()==1) {
      if (!m_trig) m_res=0.0;
      m_subs.back()->m_me=m_subs.back()->m_mewgt=m_res;
    }
    if (p_dinfo->Mode()&2) m_dsij[0][0]=m_res;
    for (size_t j(1);j<m_cur.back().size();++j) {
#ifdef DEBUG__BG
      msg_Debugging()<<"Dipole "<<m_cur.back()[j]->Sub()->Sub()->Id()
		     <<" <-> "<<m_cur.back()[j]->Sub()->Id()<<"\n";
#endif
      if (m_pmode=='D') {
	for (size_t i(0);i<m_rress.size();++i) m_cress[i]=m_rress[i]=0.0;
 	if (p_dinfo->Mode()==1)
	  m_cur.back()[j]->Contract(*m_cur.back()[j]->Sub(),m_cchirs,m_rress);
 	else
	  for (size_t i(0);i<m_rress.size();++i) m_rress[i]=m_dsf[j]*m_ress[i]; 
#ifdef DEBUG__BGS_AMAP
	for (size_t i(0);i<m_rress.size();++i) m_rress[i]=0.0;
	m_cur.back()[j]->Contract(*m_cur.back()[j]->Sub(),m_cchirs,m_rress);
	for (size_t i(0);i<m_rress.size();++i)
	  if (m_ress[i].real() || m_ress[i].imag()) {
	    msg_Debugging()<<"Check amplitude mapping ("<<m_dsf[j]<<"): "
			   <<m_rress[i]<<" vs. "<<m_dsf[j]*m_ress[i]
			   <<" -> "<<m_rress[i]/(m_dsf[j]*m_ress[i])<<"\n";
	    if (!IsEqual(m_rress[i].real(),m_dsf[j]*m_ress[i].real(),1.0e-6)) {
	      msg_Error()<<METHOD<<"(): Mapping error in";
	      for (size_t l(0);l<m_nin;++l) msg_Error()<<" "<<m_fl[l];
	      msg_Error()<<" ->";
	      for (size_t l(m_nin);l<m_nin+m_nout;++l) msg_Error()<<" "<<m_fl[l];
	      msg_Error()<<" ("<<m_cur.back()[j]->Sub()->Sub()->Id().front()<<","
			 <<m_cur.back()[j]->Sub()->Id().front()<<")["<<i<<"]: "
			 <<m_rress[i]<<" / "<<m_dsf[j]*m_ress[i]<<" = "
			 <<m_rress[i]/(m_dsf[j]*m_ress[i])<<" ("<<m_dsf[j]<<")\n";
	    }
	  }
#endif
	m_cur.back()[j]->Contract(*m_cur.back()[j]->Sub(),m_cchirs,m_cress,1);
	m_cur.back()[j]->Contract(*m_cur.back()[j]->Sub(),m_cchirs,m_cress,2);
      }
      else {
	THROW(not_implemented,"Internal error");
      }
      double ccsum(0.0);
      for (size_t i(0);i<m_rress.size();++i) {
#ifdef DEBUG__BG
	msg_Debugging()<<"A["<<i<<"]"<<m_ress(i)<<" = "
		       <<m_rress[i]<<" * "<<m_cress[i]<<" -> "
		       <<m_rress[i]*std::conj(m_cress[i])<<"\n";
#endif
	ccsum+=(m_rress[i]*std::conj(m_cress[i])).real();
      }
#ifdef DEBUG__BG
      msg_Debugging()<<"ccsum = "<<ccsum<<" for "
		     <<m_cur.back()[j]->Sub()->Id()
		     <<m_cur.back()[j]->Sub()->Sub()->Id()<<"\n";
#endif
      ccsum/=m_sf;
      if (p_dinfo->Mode()==1) {
	if (m_smth) {
	  Dipole_Kinematics *kin=m_cur.back()[j]->
	    Sub()->Sub()->In().front()->Kin();
	  if (kin->F()!=1.0) {
	    double x(dabs(kin->F()));
	    if (m_trig) {
	      m_subs.back()->m_me+=(1.0-x)*ccsum;
	      m_subs.back()->m_mewgt+=(1.0-x)*ccsum;
	    }
	    if (kin->F()<0.0) x=0.0;
	    ccsum*=x;
	  }
	}
	m_subs[m_sid[j]]->m_me=m_subs[m_sid[j]]->m_mewgt=ccsum;
      }
      else {
	Dipole_Kinematics *kin=m_cur.back()[j]->
	  Sub()->Sub()->In().front()->Kin();
	double lf(log(2.0*M_PI*mu2/EpsSchemeFactor(m_p)/
		      dabs(kin->JIJT()->P()*kin->JK()->P())));
#ifdef DEBUG__BG
	msg_Debugging()<<"e^2 = "<<kin->Res(2)<<", e = "<<kin->Res(1)
		       <<", f = "<<kin->Res(0)<<", l = "<<lf
		       <<" ( "<<p_loop<<" )\n";
#endif
	m_dsij[m_dsm[j].first][m_dsm[j].second]=ccsum;
	m_dsij[m_dsm[j].second][m_dsm[j].first]=ccsum;
	m_cmur[1]-=ccsum*asf*kin->Res(2);
	m_cmur[0]-=ccsum*asf*(kin->Res(1)+lf*kin->Res(2));
	ccsum*=-asf*(kin->Res(0)+lf*kin->Res(1)+0.5*sqr(lf)*kin->Res(2));
      }
      csum+=ccsum;
    }
    if (p_loop) {
      double cw(p_loop->Mode()?1.0:m_res);
      if (p_loop->Mode() && p_loop->ColMode()==0)
         cw*=3.0/p_colint->GlobalWeight();
      if (p_dinfo->Mode()&8) {
	double e1p(-m_cmur[0]/m_res/asf), e2p(-m_cmur[1]/m_res/asf);
	if (!IsEqual(e2p,p_loop->ME_E2()))
	  msg_Error()<<METHOD<<"(): Double pole does not match. V -> "
		     <<p_loop->ME_E2()<<", I -> "<<e2p<<", rel. diff. "
		     <<(e2p/p_loop->ME_E2()-1.0)<<".\n";
	if (!IsEqual(e1p,p_loop->ME_E1()))
	  msg_Error()<<METHOD<<"(): Single pole does not match. V -> "
		     <<p_loop->ME_E1()<<", I -> "<<e1p<<", rel. diff. "
		     <<(e1p/p_loop->ME_E1()-1.0)<<".\n";
      }
      csum+=cw*asf*p_loop->ME_Finite();
      if (m_sccmur) {
	double b0(11.0/6.0*3.0-2.0/3.0*0.5*Flavour(kf_quark).Size()/2);
	m_cmur[0]+=cw*asf*(p_loop->ME_E1()+m_maxoqcd*b0);
	m_cmur[1]+=cw*asf*p_loop->ME_E2();
      }
      else {
	m_cmur[0]+=cw*asf*p_loop->ScaleDependenceCoefficient(1);
	m_cmur[1]+=cw*asf*p_loop->ScaleDependenceCoefficient(2);
      }
    }
    if ((p_dinfo->Mode()&18) && !(p_dinfo->Mode()&4)) m_res=0.0;
  }
#ifdef DEBUG__BG
  msg_Debugging()<<"m_res = "<<m_res<<", csum = "<<csum
		 <<" -> "<<m_res+csum<<"\n";
#endif
  m_res+=csum;
  return true;
}

double Amplitude::Differential(const std::vector<Vec4D> &momenta,
			       NLO_subevt *const sub)
{
  p_sub=sub;
  return Differential(momenta,p_colint->I(),p_colint->J());
}

double Amplitude::Differential
(const std::vector<Vec4D> &momenta,
 const Int_Vector &ci,const Int_Vector &cj,const bool set)
{
  SetColors(ci,cj,set);
  if (p_helint!=NULL && p_helint->On()) {
    Evaluate(p_helint->Chiralities());
    return m_res;
  }
  EvaluateAll();
  return m_res;
}

bool Amplitude::ConstructChirs()
{
  m_ress=Spin_Structure<DComplex>(m_fl,0.0);
  m_cchirs.resize(m_ress.size());
  m_cur.back().front()->HM().resize(m_ress.size());
  for (size_t i(0);i<m_ress.size();++i) {
    m_cur.back().front()->HM()[i]=i;
    m_cchirs[i]=1;
  }
  if (p_dinfo->Mode()==0) return true;
  m_cress=m_rress=m_ress;
  for (size_t i(0);i<m_scur.size();++i) {
    Dipole_Kinematics *kin(m_scur[i]->Sub()->In().front()->Kin());
    Current *last(NULL);
    for (size_t j(1);j<m_cur.back().size();++j)
      if (m_cur.back()[j]->Sub()==m_scur[i]) last=m_cur.back()[j];
    last->HM().resize(m_ress.size(),0);
    kin->PM().resize(m_ress.size(),0);
    int ii(kin->JI()?kin->JI()->Id().front():-1);
    int ij(kin->JJ()?kin->JJ()->Id().front():-1);
    size_t ik(kin->JK()->Id().front());
    Flavour_Vector cfl(1,m_scur[i]->Flav());
    for (size_t j(0);j<m_n;++j) if (j!=ik) cfl.push_back(m_fl[j]);
    Spin_Structure<int> tch(cfl,0.0);
    for (size_t j(0);j<m_ress.size();++j) {
      Int_Vector id(m_ress.GetSpinCombination(j)), di(1,id[ik]);
      for (size_t k(0);k<m_n;++k)
	if (k!=ik) di.push_back(id[k]);
      last->HM()[tch.GetNumber(di)]=j;
      if (ii==-1 || ij==-1 || id[ii]==id[ij]) continue;
      std::swap<int>(id[ii],id[ij]);
      kin->PM()[j]=m_ress.GetNumber(id);
    }
  }
  return true;
}

bool Amplitude::Construct
(const Int_Vector &incs,const Flavour_Vector &flavs,
 MODEL::Model_Base *const model,MODEL::Coupling_Map *const cpls)
{
  p_model=model;
  m_dirs=incs;
  if (!Construct(flavs)) return false;
  if (!ConstructChirs()) return false;
  if (p_dinfo->Mode()==1) ConstructNLOEvents();
  if (p_dinfo->Mode()&2) ConstructDSijMap();
  MODEL::Coupling_Data *rqcd(cpls->Get("Alpha_QCD"));
  MODEL::Coupling_Data *rqed(cpls->Get("Alpha_QED"));
  for (size_t i(0);i<m_cur[m_n-1].size();++i) {
    for (size_t j(0);j<m_cur[m_n-1][i]->In().size();++j) {
      Vertex *v(m_cur[m_n-1][i]->In()[j]);
      int oqcd(v->OrderQCD()+v->JA()->OrderQCD()+v->JB()->OrderQCD());
      int oew(v->OrderEW()+v->JA()->OrderEW()+v->JB()->OrderEW());
      if (v->JE()) {
	oqcd+=v->JE()->OrderQCD();
	oew+=v->JE()->OrderEW();
      }
      m_cpls.push_back(Coupling_Info(v,oqcd,oew,rqcd,rqed));
    }
  }
  if (p_dinfo->Mode()==1)
  for (size_t i(0);i<m_scur.size();++i) {
    MODEL::Coupling_Data *aqcd(cpls->Get("Alpha_QCD",m_subs[i]));
    if (aqcd==NULL && rqcd) {
      aqcd = new MODEL::Coupling_Data(*rqcd,m_subs[i]);
      cpls->insert(std::make_pair("Alpha_QCD",aqcd));
    }
    MODEL::Coupling_Data *aqed(cpls->Get("Alpha_QED",m_subs[i]));
    if (aqed==NULL && rqed) {
      aqed = new MODEL::Coupling_Data(*rqed,m_subs[i]);
      cpls->insert(std::make_pair("Alpha_QED",aqed));
    }
    Current *last(NULL);
    for (size_t j(0);j<m_cur.back().size();++j)
      if (m_cur.back()[j]->Sub()==m_scur[i]) {
	last=m_cur.back()[j];
	break;
      }
    if (last==NULL) THROW(fatal_error,"Internal error");
    for (size_t j(0);j<last->In().size();++j) {
      Vertex *v(last->In()[j]);
      int oqcd(v->OrderQCD()+v->JA()->OrderQCD()+v->JB()->OrderQCD());
      int oew(v->OrderEW()+v->JA()->OrderEW()+v->JB()->OrderEW());
      if (v->JE()) {
	oqcd+=v->JE()->OrderQCD();
	oew+=v->JE()->OrderEW();
      }
      m_cpls.push_back(Coupling_Info(v,oqcd,oew,aqcd,aqed));
    }
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

void Amplitude::FillAmplitudes
(std::vector<Spin_Amplitudes> &amps,
 std::vector<std::vector<Complex> > &cols)
{
  cols.push_back(std::vector<Complex>(1,1.0));
  amps.push_back(Spin_Amplitudes(m_fl,Complex(0.0,0.0)));
  for (size_t i(0);i<m_ress.size();++i)
    amps.back().Insert(Complex(m_ress[i]),m_ress(i));
}

double Amplitude::Coupling(const int mode) const
{
  MODEL::Coupling_Data *cpl(m_cpls.front().p_aqcd);
  return cpl->Default()*cpl->Factor();
}

void Amplitude::FillMEWeights(ME_wgtinfo &wgtinfo) const
{
  if (wgtinfo.m_nx<2) return;
  for (int i=0;i<2;i++) wgtinfo.p_wx[i]=m_cmur[i];
}

void Amplitude::SetGauge(const size_t &n)
{
  Vec4D k(1.0,0.0,1.0,0.0);
  switch(n) {
  case 1: k=Vec4D(1.0,0.0,invsqrttwo,invsqrttwo); break;
  case 2: k=Vec4D(1.0,invsqrttwo,0.0,invsqrttwo); break;
  case 3: k=Vec4D(1.0,invsqrttwo,invsqrttwo,0.0); break;
  }
  for (size_t j(1);j<m_cur.size();++j)
    for (size_t i(0);i<m_cur[j].size();++i) m_cur[j][i]->SetGauge(k);
  for (size_t i(0);i<m_scur.size();++i) m_scur[i]->SetGauge(k);
}

bool Amplitude::GaugeTest(const Vec4D_Vector &moms,const int mode)
{
  if (mode==0) {
    size_t nt(0);
    bool cnt(true);
    while (cnt) {
      while (!p_colint->GeneratePoint());
      SetColors(p_colint->I(),p_colint->J());
      if (!GaugeTest(moms,1)) return false;
      if (m_res!=0.0) cnt=false;
      if (cnt) {
	if (++nt>100) 
	  msg_Error()<<METHOD<<"(): Zero result. Redo gauge test."<<std::endl;
	while (!p_colint->GeneratePoint());
      }
    }
    return true;
  }
  msg_Tracking()<<METHOD<<"(): Performing gauge test ..."<<std::flush;
  msg_Indent();
  if (m_pmode=='D') {
    int sd(Spinor<double>::DefaultGauge());
    Spinor<double>::SetGauge(sd>0?sd-1:sd+1);
  }
  SetGauge(1);
  SetMomenta(moms);
  if (!EvaluateAll()) return false;
  Spin_Structure<DComplex> ress(m_ress);
  if (m_pmode=='D') Spinor<double>::ResetGauge();
  SetGauge(0);
  SetMomenta(moms);
  if (!EvaluateAll()) return false;
  double mean(0.0);
  for (size_t i(0);i<m_ress.size();++i) mean+=std::abs(ress[i]);
  mean/=m_ress.size();
  msg_Debugging()<<METHOD<<"(): {\n";
#ifdef USE__Strict_QCD_Gauge_Test
  if (m_maxoew>0) {
#else
  if (true) {
#endif
    double xs(0.0), rxs(0.0);
    for (size_t i(0);i<m_ress.size();++i) {
      xs+=sqr(std::abs(m_ress[i]));
      rxs+=sqr(std::abs(ress[i]));
      msg_Debugging()<<"A("<<m_ress(i)<<") = "<<std::abs(m_ress[i])
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
    if (!IsEqual(xs,rxs,rpa->gen.Accu())) {
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
      msg_Debugging()<<"A("<<m_ress(i)
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
  msg_Tracking()<<METHOD<<"(): Write diagrams to '"<<file<<"'.\n";
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
