#include "COMIX/Amplitude/Amplitude.H"

#include "MODEL/Main/Single_Vertex.H"
#include "PHASIC++/Main/Color_Integrator.H"
#include "PHASIC++/Main/Helicity_Integrator.H"
#include "PHASIC++/Process/Process_Base.H"
#include "METOOLS/Explicit/Dipole_Kinematics.H"
#include "METOOLS/Main/Spin_Structure.H"
#include "ATOOLS/Phys/Spinor.H"
#include "ATOOLS/Math/Permutation.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/My_Limits.H"
#include "ATOOLS/Org/STL_Tools.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/My_File.H"
#include "ATOOLS/Org/Scoped_Settings.H"

#include <algorithm>
#include <cassert>

using namespace COMIX;
using namespace ATOOLS;

static bool csscite(false);

static const double invsqrttwo(1.0/sqrt(2.0));

template<typename Type> inline Type
GetParameter(const std::string &name)
{ return ToType<Type>(rpa->gen.Variable(name)); }

Amplitude::Amplitude():
  p_model(NULL), m_nin(0), m_nout(0), m_dec(0), m_n(0), m_wfmode(0), m_ngpl(3),
  m_maxcpl(2,99), m_mincpl(2,0), m_maxacpl(2,99), m_minacpl(2,0),
  m_minntc(0), m_maxntc(99), m_pmode('D'), p_dinfo(new Dipole_Info()),
  p_colint(NULL), p_helint(NULL), m_trig(true), p_loop(NULL)
{
  Settings& s = Settings::GetMainSettings();
  Scoped_Settings comixsettings{ s["COMIX"] };
  m_sccmur = s["USR_WGT_MODE"].Get<bool>();
  m_murcoeffvirt = s["NLO_MUR_COEFFICIENT_FROM_VIRTUAL"].Get<bool>();
  p_dinfo->SetMassive(0);
  m_pmode = comixsettings["PMODE"].Get<std::string>()[0];
  m_wfmode = comixsettings["WF_MODE"].Get<int>();
  m_pgmode = comixsettings["PG_MODE"].Get<int>();
  m_ngpl = Min(1, Max(5, comixsettings["N_GPL"].Get<int>()));
  p_dinfo->SetAMin(s["DIPOLES"]["AMIN"].Get<double>());
  const auto amax = s["DIPOLES"]["ALPHA"].Get<double>();
  double cur{ 0.0 };
  cur = s["DIPOLES"]["ALPHA_FF"].Get<double>();
  p_dinfo->SetAMax(0,cur?cur:amax);
  cur = s["DIPOLES"]["ALPHA_FI"].Get<double>();
  p_dinfo->SetAMax(2,cur?cur:amax);
  cur = s["DIPOLES"]["ALPHA_IF"].Get<double>();
  p_dinfo->SetAMax(1,cur?cur:amax);
  cur = s["DIPOLES"]["ALPHA_II"].Get<double>();
  p_dinfo->SetAMax(3,cur?cur:amax);
  p_dinfo->SetKappa(s["DIPOLES"]["KAPPA"].Get<double>());
  p_dinfo->SetNf(s["DIPOLES"]["NF_GSPLIT"].Get<int>());
  p_dinfo->SetKT2Max(s["DIPOLES"]["KT2MAX"].Get<double>());
  p_dinfo->SetDRMode(0);
  const subscheme::code subtype{ s["DIPOLES"]["SCHEME"].Get<subscheme::code>() };
  p_dinfo->SetSubType(subtype);
  if (subtype==subscheme::Dire) p_dinfo->SetKappa(1.0);
  m_smth=GetParameter<double>("NLO_SMEAR_THRESHOLD");
  m_smpow=GetParameter<double>("NLO_SMEAR_POWER");
}

Amplitude::~Amplitude()
{
  if (p_dinfo) delete p_dinfo;
  CleanUp();
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
  m_cress=m_ress=std::vector<Spin_Structure<DComplex> >();
  m_on=m_son=std::vector<std::pair<size_t,size_t> >();
  m_cchirs=m_dirs=Int_Vector();
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
  cur->SetOrder(c->Order());
  return cur;
}

bool Amplitude::AddRSDipoles()
{
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"(): {\n";
#endif
  Current_Vector scur;
  if (m_stype&1) {
  for (size_t j(0);j<m_n;++j) {
    Flavour flc(m_cur[1][j]->Flav());
    for (size_t i(0);i<m_cur[2].size();++i) {
      if (m_cur[2][i]->CId()&m_cur[1][j]->CId()) continue;
      Vertex *v(m_cur[2][i]->In().front());
      Flavour fla(v->J(0)->Flav()), flb(v->J(1)->Flav());
      bool isa(v->J(0)->CId()&(1<<m_nin)-1);
      bool isb(v->J(1)->CId()&(1<<m_nin)-1);
      if (!(v->Order(0)==1 && flc.StrongCharge())) continue;
      if (m_cur[2][i]->In().size()!=1)
	THROW(not_implemented,"Invalid current");
      if (!AddRSDipole(m_cur[2][i],m_cur[1][j],scur,1)) return false;
    }
  }
  }
  if (m_stype&2) {
  for (size_t j(0);j<m_n;++j) {
    Flavour flc(m_cur[1][j]->Flav());
    for (size_t i(0);i<m_cur[2].size();++i) {
      if (m_cur[2][i]->CId()&m_cur[1][j]->CId()) continue;
      Vertex *v(m_cur[2][i]->In().front());
      Flavour fla(v->J(0)->Flav()), flb(v->J(1)->Flav());
      bool isa(v->J(0)->CId()&(1<<m_nin)-1);
      bool isb(v->J(1)->CId()&(1<<m_nin)-1);
      if (!((/*!isa &&*/ fla.IsPhoton() && flc.Charge()) ||
	    (/*!isb &&*/ flb.IsPhoton() && flc.Charge()))) continue;
      if (m_cur[2][i]->In().size()!=1)
	THROW(not_implemented,"Invalid current");
      if (!AddRSDipole(m_cur[2][i],m_cur[1][j],scur,2)) return false;
    }
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
(Current *const c,Current *const s,Current_Vector &scur,int stype)
{
#ifdef DEBUG__BG
  msg_Indent();
#endif
  size_t isid((1<<m_nin)-1);
  if ((c->CId()&isid)==isid || c->Flav().IsDummy()) return true;
  Vertex *cin(c->In().front());
  if (cin->J(0)->Direction()>0 ||
      cin->J(1)->Direction()>0) {
    if (cin->J(0)->Flav().Mass() ||
	cin->J(1)->Flav().Mass()) return true;
  }
  else {
    double mm(Flavour(p_dinfo->Nf()).Mass());
    if (cin->J(0)->Flav().Mass()>mm &&
	cin->J(1)->Flav().Mass()>mm) return true;
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
  Vertex_Key *svkey(Vertex_Key::New(cin->J(),jijt,p_model));
  svkey->m_p=std::string(1,m_pmode);
  svkey->p_dinfo=p_dinfo;
  svkey->m_stype=stype;
  svkey->p_k=s;
  svkey->p_kt=jkt;
  MODEL::VMIterator_Pair vmp(p_model->GetVertex(svkey->ID()));
  for (MODEL::Vertex_Map::const_iterator vit(vmp.first);
       vit!=vmp.second;++vit) {
    svkey->p_mv=vit->second;
    Vertex *v(new Vertex(*svkey));
    v->AddJ(svkey->m_j);
    v->SetJC(svkey->p_c);
  }
  svkey->Delete();
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
      if (m_stype==2) {
	sc[0]=m_cur[1][i]->Flav().IntCharge();
	sc[1]=m_cur[1][j]->Flav().IntCharge();
      }
      if (sc[0]==0 || sc[1]==0) continue;
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
  Current_Vector j(2,c);
  j[1]=NULL;
  Vertex_Key *svkey(Vertex_Key::New(j,jijt,p_model));
  svkey->m_p=std::string(1,m_pmode);
  svkey->p_dinfo=p_dinfo;
  svkey->p_k=s;
  svkey->p_kt=jkt;
  MODEL::VMIterator_Pair vmp(p_model->GetVertex(svkey->ID()));
  for (MODEL::Vertex_Map::const_iterator vit(vmp.first);
       vit!=vmp.second;++vit) {
    svkey->p_mv=vit->second;
    Vertex *v(new Vertex(*svkey));
    v->AddJ(svkey->m_j);
    v->SetJC(svkey->p_c);
  }
  if (svkey->p_mv==NULL) {
    std::swap<Current*>(svkey->m_j[0],svkey->m_j[1]);
    vmp=p_model->GetVertex(svkey->ID());
    for (MODEL::Vertex_Map::const_iterator vit(vmp.first);
	 vit!=vmp.second;++vit) {
      svkey->p_mv=vit->second;
      Vertex *v(new Vertex(*svkey));
      v->AddJ(svkey->m_j);
      v->SetJC(svkey->p_c);
    }
  }
  svkey->Delete();
#ifdef DEBUG__BG
  jijt->Print();
  jkt->Print();
#endif
  return true;
}

bool Amplitude::MatchDecay(const Vertex_Key &vkey) const
{
  std::vector<size_t> c(vkey.m_j.size());
  for (size_t j(0);j<c.size();++j) {
    size_t jid(vkey.m_j[j]->CId());
    for (size_t i(0);i<m_decid.size();++i) {
      size_t did(m_decid[i]->m_id);
      if ((did&jid) && (did&~jid)) c[j]|=(1<<i);
    }
  }
  for (size_t i(1);i<c.size();++i)
    if (c[i]!=c[0]) return false;
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
 const int dec,std::vector<int> &maxcpl,std::vector<int> &mincpl,
 std::map<std::string,Current*> &curs)
{
  if (vkey.p_mv==NULL || 
      vkey.p_mv->dec<0) return NULL;
  vkey.m_p=std::string(1,m_pmode);
  Vertex *v(new Vertex(vkey));
  size_t ntc(0), ist(0);
  std::vector<int> order(v->Order());
  for (size_t i(0);i<vkey.m_j.size();++i) {
    if (order.size()<vkey.m_j[i]->Order().size())
      order.resize(vkey.m_j[i]->Order().size(),0);
    for (size_t j(0);j<vkey.m_j[i]->Order().size();++j)
      order[j]+=vkey.m_j[i]->Order(j);
    ntc+=vkey.m_j[i]->NTChannel();
    if (bool(vkey.m_j[i]->CId()&1)^bool(vkey.m_j[i]->CId()&2)) ++ist;
  }
  if (ist==1) ntc+=!ckey.m_fl.Strong();
  bool valid(true);
  for (size_t i(0);i<Min(order.size(),m_maxcpl.size());++i)
    if (order[i]>m_maxcpl[i]) valid=false;
  for (size_t i(0);i<Min(order.size(),m_maxacpl.size());++i)
    if (order[i]>m_maxacpl[i]) valid=false;
  if (!v->Active() || !valid ||
      (n==m_n-1 && (ntc<m_minntc || ntc>m_maxntc))) {
#ifdef DEBUG__BG
    msg_Debugging()<<"delete vertex "<<vkey.ID()<<"-"<<vkey.Type()
		   <<v->Order()<<")->{"<<vkey.p_c->Flav()<<"} => "
		   <<order<<" vs. max = "<<m_maxcpl<<"/"<<m_maxacpl
		   <<", act = "<<v->Active()<<", n t-ch = "
		   <<ntc<<" vs "<<m_minntc<<"/"<<m_maxntc<<"\n";
#endif
    delete v;
    return NULL;
  }
  Current *sub(NULL);
  if (p_dinfo->Mode()) {
    size_t cid(0);
    for (size_t i(0);i<vkey.m_j.size();++i) {
      cid|=vkey.m_j[i]->CId();
      if (vkey.m_j[i]->Sub()) {
	if ((vkey.m_j[i]->Id().size()==1 &&
	     p_dinfo->Mode()==1) || sub) return NULL;
	sub=vkey.m_j[i]->Sub();
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
    ("("+ToString(order)+(n<m_n-1?";"+ToString(ntc):"")
     +(sub?",S"+ToString(sub->Sub()->Id())+
       ToString(sub->Sub()->Sub()->Id()):"")+")");
  if (order!=vkey.p_c->Order() ||
      ntc!=vkey.p_c->NTChannel() || sub!=vkey.p_c->Sub()) {
    std::map<std::string,Current*>::iterator 
      cit(curs.find(okey));
    if (cit!=curs.end()) vkey.p_c=cit->second;
    else {
      if (vkey.p_c->Order().size() ||
	  vkey.p_c->NTChannel()>0 || vkey.p_c->Sub()!=sub)
	vkey.p_c=Current_Getter::GetObject
	  (std::string(1,m_pmode)+ckey.Type(),ckey);
      if (n<m_n-1) {
	if (dec!=0) {
	  vkey.p_c->SetCut(m_decid[dec-1]->m_nmax);
	  vkey.p_c->SetOnShell(m_decid[dec-1]->m_osd);
	}
	vkey.p_c->SetNTChannel(ntc);
      }
      else {
	if (maxcpl.size()<order.size()) maxcpl.resize(order.size(),0);
	if (mincpl.size()<order.size()) mincpl.resize(order.size(),99);
	for (size_t i(0);i<order.size();++i) {
	  maxcpl[i]=Max(maxcpl[i],order[i]);
	  mincpl[i]=Min(mincpl[i],order[i]);
	}
      }
      vkey.p_c->SetOrder(order);
      vkey.p_c->SetSub(sub);
      curs[okey]=vkey.p_c;
    }
  }
  v->AddJ(vkey.m_j);
  v->SetJC(vkey.p_c);
  return v;
}

void Amplitude::AddCurrent(const Int_Vector &ids,const size_t &n,
			   const Flavour &fl,const int dir)
{
  // add new currents
  std::vector<int> maxcpl, mincpl;
  int dec(CheckDecay(fl,ids));
  if (dec<0) return;
  size_t cid(0);
  for (size_t i(0);i<ids.size();++i) cid|=1<<ids[i];
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
  bool one(false);
  size_t nmax(Min(n,p_model->MaxLegs(ckey.m_fl.Bar())-1));
  for (size_t nc(1);nc<nmax;++nc) {
    Int_Vector jc(nc+1,1);
    Current_Matrix jj(nc+1,m_cur[1]);
    for (size_t jl(jc.size()-2);jc[0]<=n-nc;++jc[jl]) {
      if (jc[jl]>n-nc) { jc[jl--]=1; continue; }
      jj[jl]=m_cur[jc[jl]];
      if(jl<jc.size()-2) { --jc[++jl]; continue; }
      jc.back()=n;
      bool zero(false);
      for (size_t i(0);i<jc.size()-1;++i) {
	if (jj[i].empty()) zero=true;
	jc.back()-=jc[i];
      }
      if (zero || jc.back()<=0) continue;
      if ((jj.back()=m_cur[jc.back()]).empty()) continue;
      Int_Vector cc(jj.size(),0);
      Current_Vector cj(cc.size(),NULL);
      Vertex_Key *vkey(Vertex_Key::New(cj,cur,p_model));
      for (size_t i(0);i<cj.size();++i) cj[i]=jj[i].front();
      for (size_t cl(cc.size()-1);cc[0]<jj[0].size();++cc[cl]) {
        if (cc[cl]==jj[cl].size()) { cc[cl--]=0; continue; }
	cj[cl]=jj[cl][cc[cl]];
	if (cl<cc.size()-1) { --cc[++cl]; continue; }
	bool ord(true);
	for (size_t i(0);i<cj.size()-1;++i)
	  if (cj[i]->Id().front()>=cj[i+1]->Id().front())
	    { ord=false; break; }
	if (!ord) continue;
	size_t tid(cid);
	for (size_t i(0);i<cj.size();++i) {
	  vkey->m_j[i]=cj[i];
	  size_t cur(cj[i]->CId());
	  if ((tid&cur)!=cur) break;
	  tid&=~cur;
	}
	if (tid) continue;
	if (m_dec && !MatchDecay(*vkey)) continue;
	Permutation perm(cj.size());
	for (int nperm(perm.MaxNumber()), i(0);i<nperm;++i) {
	  int f(0), *p(perm.Get(i));
	  for (size_t i(0);i<cj.size();++i) vkey->m_j[i]=cj[p[i]];
	  MODEL::VMIterator_Pair vmp(p_model->GetVertex(vkey->ID()));
	  for (MODEL::Vertex_Map::const_iterator
		 vit(vmp.first);vit!=vmp.second;++vit) {
	    vkey->p_mv=vit->second;
	    if (AddCurrent(ckey,*vkey,n,dec,maxcpl,mincpl,curs)) f=true;
	  }
	  if (f) { one=true; break; }
	}
      }
      vkey->Delete();
    }
  }
  if (!one && n>1) {
    delete cur;
    return;
  }
  if (n==1) curs[""]=cur;
  else delete cur;
  Int_Vector isfs(ids.size()), pols(ids.size());
  for (size_t i(0);i<ids.size();++i) {
    isfs[i]=m_fl[ids[i]].IsFermion();
    switch (m_fl[ids[i]].IntSpin()) {
    case 0: pols[i]=1; break;
    case 1: pols[i]=2; break;
    case 2: pols[i]=m_fl[ids[i]].IsMassive()?3:2; break;
    case 3: pols[i]=m_fl[ids[i]].IsMassive()?5:2; break;
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
	size_t cid(0);
	for (size_t i(0);i<ids.size();++i) cid|=1<<ids[i];
	const std::vector<long int> &caffm(m_affm[n][cid]);
	for (size_t i(0);i<caffm.size();++i) {
	  Flavour cfl((kf_code)std::abs(caffm[i]),caffm[i]<0);
	  if (cfl.IsOn()) AddCurrent(ids,n,cfl,0);
	}
      }
      else {
	for (size_t i(0);i<fls.size();++i)
	  if (fls[i].IsOn()) AddCurrent(ids,n,fls[i],0);
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
  Prune();
  if (p_dinfo->Mode()) {
    for (size_t i(0);i<m_cur.back().size();++i)
      if (m_cur.back()[i]->Sub()==NULL && 
	  (i && m_cur.back()[i-1]->Sub())) {
	Current *c(m_cur.back()[i]);
	m_cur.back().erase(m_cur.back().begin()+i);
	m_cur.back().insert(m_cur.back().begin(),c);
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

void Amplitude::Prune()
{
  for (size_t j(m_n-2);j>1;--j)
    for (Current_Vector::iterator cit(m_cur[j].begin());
	 cit!=m_cur[j].end();++cit)
      if ((*cit)->Dangling()) {
#ifdef DEBUG__BG
	msg_Debugging()<<"delete current "<<**cit<<", "<<(*cit)->Dangling()
		       <<", O"<<(*cit)->Order()<<" vs. O_{min}"
		       <<m_mincpl<<"/"<<m_minacpl<<" .. O_{max}"
		       <<m_maxcpl<<"/"<<m_maxacpl
		       <<", n t-ch = "<<(*cit)->NTChannel()<<" vs. "
		       <<m_minntc<<"/"<<m_maxntc<<"\n";
#endif
	if ((*cit)->Sub() && (*cit)->Sub()->Sub()==*cit)
	  (*cit)->Sub()->SetSub(NULL);
	delete *cit;
	cit=--m_cur[j].erase(cit);
      }
}

bool Amplitude::ReadInAmpFile(const std::string &name,My_In_File &amp)
{
  if (amp()==NULL) return false;
  amp->seekg(0);
  if (!amp->good()) return false;
  std::string cname, cmname;
  *amp>>cname>>cmname;
  if (cname!=name || cmname!=name || amp->eof())
    THROW(fatal_error,"Corrupted map file '"+name+".map'");
  m_affm.resize(m_n);
  for (size_t size, i(2);i<m_n;++i) {
    while (true) {
      long int kfc;
      *amp>>kfc;
      if (kfc==0) break;
      std::vector<long int> &caffm(m_affm[i][kfc]);
      *amp>>size;
      caffm.resize(size);
      for (size_t j(0);j<size;++j) *amp>>caffm[j];
    }
  }
  *amp>>cname;
  if (cname!="eof")
    THROW(fatal_error,"Corrupted map file '"+name+".map'");
  return true;
}

void Amplitude::WriteOutAmpFile(const std::string &name)
{
  std::string ampfile(rpa->gen.Variable("SHERPA_CPP_PATH")
		      +"/Process/Comix/"+name+".map");
  if (FileExists(ampfile,1)) return;
  My_Out_File amp(ampfile);
  if (!amp.Open()) return;
  *amp<<name<<" "<<name<<"\n";
  m_affm.resize(m_n);
  for (size_t i(2);i<m_n;++i) {
    std::map<size_t,std::set<long int> > kfcs;
    m_affm[i].clear();
    for (size_t j(0);j<m_cur[i].size();++j) {
      std::set<long int> &ckfcs(kfcs[m_cur[i][j]->CId()]);
      std::vector<long int> &caffm(m_affm[i][m_cur[i][j]->CId()]);
      long int kfc(m_cur[i][j]->Flav());
      if (ckfcs.find(kfc)==ckfcs.end()) caffm.push_back(kfc);
      ckfcs.insert(kfc);
    }
    for (size_t j(0);j<m_affm[i].size();++j)
      for (std::map<size_t,std::vector<long int> >::const_iterator
	     it(m_affm[i].begin());it!=m_affm[i].end();++it) {
	*amp<<it->first<<" "<<it->second.size()<<" ";
	for (size_t k(0);k<it->second.size();++k) *amp<<it->second[k]<<" ";
      }
    *amp<<"0\n";
  }
  *amp<<"eof\n";
}

bool Amplitude::Initialize
(const size_t &nin,const size_t &nout,const std::vector<Flavour> &flavs,
 const double &isf,const double &fsf,My_In_File &ampfile,
 MODEL::Model_Base *const model,MODEL::Coupling_Map *const cpls,
 const int stype,const int smode,const cs_itype::type itype,
 const std::vector<int> &maxcpl, const std::vector<int> &mincpl,
 const std::vector<int> &maxacpl, const std::vector<int> &minacpl,
 const size_t &minntc,const size_t &maxntc,const std::string &name)
{
  DEBUG_FUNC(flavs<<", stype="<<(sbt::subtype)(stype+1)<<", smode="<<smode
             <<", itype="<<itype<<", cpls="<<mincpl<<"/"<<minacpl
	     <<".."<<maxcpl<<"/"<<maxacpl);
  CleanUp();
  m_nin=nin;
  m_nout=nout;
  m_n=m_nin+m_nout;
  m_minntc=minntc;
  m_maxntc=maxntc;
  m_stype=stype;
  m_maxcpl=maxcpl;
  m_mincpl=mincpl;
  m_maxacpl=maxacpl;
  m_minacpl=minacpl;
  p_dinfo->SetMode(smode);
  p_dinfo->SetIType(itype);
  ReadInAmpFile(name,ampfile);
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
    bool order(true);
    if (sub->m_i>m_nin && sub->m_j>m_nin &&
	kin->JI()->Flav().IsBoson() &&
        kin->JJ()->Flav().IsFermion()) order=false;
    sub->m_oqcd=m_maxcpl[0]/2;
    sub->m_oew=m_maxcpl[1]/2;
    Current_Vector cur(sub->m_n,NULL);
    for (size_t k(0), j(0);k<m_nin+m_nout;++k) {
      if ((k==sub->m_j && order) ||
          (k==sub->m_i && !order)) continue;
      if ((k==sub->m_i && order) ||
          (k==sub->m_j && !order)) {
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
    kin->SetSubevt(sub);
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
    sub->m_pname=PHASIC::Process_Base::GenerateName(cpi.m_ii,cpi.m_fi);
    sub->m_stype=(sbt::subtype)(m_stype);
    msg_Indent();
    msg_Debugging()<<*sub<<"\n";
  }
  NLO_subevt *sub(new NLO_subevt());
  m_subs.push_back(sub);
  sub->m_idx=m_subs.size()-1;
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
  sub->m_stype=sbt::none;
  sub->m_i=sub->m_j=sub->m_k=0;
  sub->m_oqcd=m_maxcpl[0]/2;
  sub->m_oew=m_maxcpl[1]/2;
  {
    msg_Indent();
    msg_Debugging()<<*sub<<"\n";
  }
  for (size_t i(0);i<m_subs.size();++i) m_subs[i]->p_real=sub;
  for (size_t i(0);i<m_cur.back().size();++i)
    if (m_cur.back()[i]->Sub()==NULL) m_sid[i]=m_subs.size()-1;
    else
      for (size_t j(0);j<m_scur.size();++j)
      if (m_cur.back()[i]->Sub()==m_scur[j]) {
	m_sid[i]=j;
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
  if (m_stype==1) {
    for (size_t i(0);i<m_fl.size();++i)
      plist[i]=m_fl[i].Strong()?c++:0;
  }
  else {
    for (size_t i(0);i<m_fl.size();++i)
      plist[i]=m_fl[i].Charge()?c++:0;
  }
  m_dsij.resize(c,std::vector<double>(c));
  Flavour ifl(m_cur[1][0]->Flav());
  if (ifl.IsFermion() && !ifl.IsAnti()) {
    for (int i(m_fl.size()-1);i>0;--i)
      if (m_cur[1][i]->Flav().IsFermion()) ++ifp;
  }
  m_dsf.resize(m_cur.back().size(),
	       std::pair<size_t,double>(-1,ifp%2?-1.0:1.0));
  for (size_t j(0);j<m_cur.back().size();++j) {
    if (m_cur.back()[j]->Sub()==NULL) continue;
    m_dsm[j]=std::pair<int,int>
      (plist[m_cur.back()[j]->Sub()->Sub()->Id().front()],
       plist[m_cur.back()[j]->Sub()->Id().front()]);
    for (size_t i(0);i<m_cur.back().size();++i)
      if (m_cur.back()[i]->Sub()==NULL &&
	  m_cur.back()[i]->Order()==m_cur.back()[j]->Order()) {
	m_dsf[j].first=i;
	break;
      }
    Flavour ffl(m_cur.back()[j]->Sub()->Flav());
    if (!ffl.IsFermion()) continue;
    int idx(m_cur.back()[j]->Sub()->Id().front());
    if (ffl.IsAnti()) {
      for (int i(0);i<idx;++i)
	if (m_cur[1][i]->Flav().IsFermion())
	  m_dsf[j].second=-m_dsf[j].second;
    }
    else {
      for (int i(m_fl.size()-1);i>idx;--i)
	if (m_cur[1][i]->Flav().IsFermion())
	  m_dsf[j].second=-m_dsf[j].second;
    }
  }
}

size_t Amplitude::BornID(const size_t &id,const NLO_subevt *sub) const
{
  if (sub==NULL) return id;
  Int_Vector mid(ID(id));
  for (size_t n(0);n<mid.size();++n)
    for (size_t m(0);m<sub->m_n;++m)
      if (sub->p_id[m]&(1<<mid[n])) {
	mid[n]=m;
	break;
      }
  size_t cid(0);
  for (Int_Vector::iterator it(mid.begin());it!=mid.end();)
    if (cid&(1<<*it)) it=mid.erase(it);
    else cid|=1<<*it++;
  return cid;
}

void Amplitude::FillCombinations
(Combination_Set &combs,CFlavVector_Map &flavs,
 SizeT_Map *brs,const NLO_subevt *sub) const
{
  size_t idk(sub?1<<sub->m_k:0), idij(sub?(1<<sub->m_i)|(1<<sub->m_j):0);
  msg_Debugging()<<METHOD<<"(sub = "<<ID(idij)<<"<->"<<ID(idk)<<"): {\n";
  msg_Debugging()<<"  flavours {\n";
  for (size_t i(2);i<m_n-1;++i)
    for (size_t j(0);j<m_cur[i].size();++j) {
      if (m_cur[i][j]->Flav().IsDummy()) continue;
      if (sub==NULL) {
	if (m_cur[i][j]->Sub()) continue;
      }
      else {
	if (m_cur[i][j]->CId()&(idij|idk))
	  if (m_cur[i][j]->Sub()==NULL ||
	      m_cur[i][j]->Id().size()==2 ||
	      m_cur[i][j]->Sub()->CId()!=idk ||
	      m_cur[i][j]->Sub()->Sub()->CId()!=idij) continue;
      }
      size_t id(BornID(m_cur[i][j]->CId(),sub));
      size_t cid(BornID((1<<m_n)-1-m_cur[i][j]->CId(),sub));
      Flavour_Vector &fs(flavs[id]), &cfs(flavs[cid]);
      Flavour f(m_cur[i][j]->Flav()), cf(f.Bar());
      if (std::find(fs.begin(),fs.end(),f)==fs.end()) fs.push_back(f);
      if (std::find(cfs.begin(),cfs.end(),cf)==cfs.end()) cfs.push_back(cf);
      msg_Debugging()<<"    "<<ID(id)<<" / "<<ID(cid)
		     <<" -> "<<m_cur[i][j]->Flav()<<"\n";
    }
  msg_Debugging()<<"  } -> "<<flavs.size()<<"\n";
  msg_Debugging()<<"  combinations {\n";
  for (size_t i(2);i<m_n;++i)
    for (size_t j(0);j<m_cur[i].size();++j) {
      if (m_cur[i][j]->Flav().IsDummy()) continue;
      if (sub==NULL) {
	if (m_cur[i][j]->Sub()) continue;
      }
      else {
	if (m_cur[i][j]->CId()&(idij|idk))
	  if (m_cur[i][j]->Sub()==NULL ||
	      m_cur[i][j]->Id().size()==2 ||
	      m_cur[i][j]->Sub()->CId()!=idk ||
	      m_cur[i][j]->Sub()->Sub()->CId()!=idij) continue;
      }
      Vertex_Vector ins(m_cur[i][j]->In());
      for (size_t k(0);k<ins.size();++k) {
	if (ins[k]->J().size()>2 ||
	    ins[k]->J(0)->Flav().IsDummy() ||
	    ins[k]->J(1)->Flav().IsDummy()) continue;
	size_t ida(BornID(ins[k]->J(0)->CId(),sub));
	size_t idb(BornID(ins[k]->J(1)->CId(),sub));
	size_t idc(BornID((1<<m_n)-1-ins[k]->JC()->CId(),sub));
	if (brs) {
	  (*brs)[ida]=ins[k]->J(0)->CId();
	  (*brs)[idb]=ins[k]->J(1)->CId();
	  (*brs)[idc]=(1<<m_n)-1-ins[k]->JC()->CId();
	}
	msg_Debugging()<<"    "<<ID(ida)
		       <<" "<<ID(idb)<<" "<<ID(idc)<<"\n";
	combs.insert(std::pair<size_t,size_t>(ida,idb));
	combs.insert(std::pair<size_t,size_t>(idb,ida));
	combs.insert(std::pair<size_t,size_t>(idb,idc));
	combs.insert(std::pair<size_t,size_t>(idc,idb));
	combs.insert(std::pair<size_t,size_t>(idc,ida));
	combs.insert(std::pair<size_t,size_t>(ida,idc));
      }
    }
  msg_Debugging()<<"  } -> "<<combs.size()<<"\n";
  msg_Debugging()<<"}\n";
}

bool Amplitude::Map(const Amplitude &ampl,Flavour_Map &flmap)
{
  flmap.clear();
  msg_Debugging()<<METHOD<<"(): {\n";
  for (size_t n(1);n<m_n;++n) {
    if (ampl.m_cur[n].size()!=m_cur[n].size()) {
      msg_Debugging()<<"  current count differs\n} no match\n";
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
	  flmap.clear();
	  return false;
	}
	if (ampl.m_cur[n][i]->Flav().Mass()!=
	    m_cur[n][i]->Flav().Mass() ||
	    ampl.m_cur[n][i]->Flav().Width()!=
	    m_cur[n][i]->Flav().Width()) {
	  msg_Debugging()<<"    mass or width differs\n  }\n";
	  msg_Debugging()<<"} no match\n";
	  flmap.clear();
	  return false;
	}
	if (ampl.m_cur[n][i]->Flav().StrongCharge()!=
	    m_cur[n][i]->Flav().StrongCharge()) {
	  msg_Debugging()<<"    color structure differs\n  }\n";
	  msg_Debugging()<<"} no match\n";
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
	  flmap.clear();
	  return false;
	}
      }
      Vertex_Vector vin(m_cur[n][i]->In());
      Vertex_Vector avin(ampl.m_cur[n][i]->In());
      if (avin.size()!=vin.size()) {
	msg_Debugging()<<"    vertex count differs\n  }\n";
	msg_Debugging()<<"} no match\n";
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
  return true;
}

#ifdef USING__Threading
void *Amplitude::TCalcJL(void *arg)
{
  CDBG_ME_TID *tid((CDBG_ME_TID*)arg);
  pthread_mutex_lock(&tid->m_s_mtx);
  while (true) {
    pthread_cond_wait(&tid->m_s_cnd,&tid->m_s_mtx);
    if (tid->m_s==0) return NULL;
    for (tid->m_i=tid->m_b;tid->m_i<tid->m_e;++tid->m_i) 
      tid->p_ampl->m_cur[tid->m_n][tid->m_i]->Evaluate();
    pthread_cond_signal(&tid->m_t_cnd,&tid->m_t_mtx);
  }
  pthread_mutex_unlock(&tid->m_s_mtx);
  return NULL;
}
#endif

void Amplitude::CalcJL()
{
  SetCouplings();
  for (size_t i(0);i<m_n;++i) 
    m_cur[1][i]->ConstructJ(m_p[i],m_ch[i],m_cl[i][0],m_cl[i][1],m_wfmode);
  for (size_t i(m_n);i<m_cur[1].size();++i) m_cur[1][i]->Evaluate();
  for (size_t n(2);n<m_n;++n) {
#ifndef USING__Threading
    for (size_t i(0);i<m_cur[n].size();++i)
      m_cur[n][i]->Evaluate();
#else
    if (p_cts->empty())
      for (size_t i(0);i<m_cur[n].size();++i) 
	m_cur[n][i]->Evaluate();
    else {
      size_t d(m_cur[n].size()/p_cts->size());
      if (m_cur[n].size()%p_cts->size()>0) ++d;
      for (size_t j(0), i(0);j<p_cts->size()&&i<m_cur[n].size();++j) {
	CDBG_ME_TID *tid((*p_cts)[j]);
	tid->p_ampl=this;
	tid->m_n=n;
	tid->m_b=i;
	tid->m_e=Min(i+=d,m_cur[n].size());
	pthread_cond_signal(&tid->m_s_cnd,&tid->m_s_mtx);
      }
      for (size_t j(0), i(0);j<p_cts->size()&&i<m_cur[n].size();++j) {
	i+=d;
	CDBG_ME_TID *tid((*p_cts)[j]);
	pthread_cond_wait(&tid->m_t_cnd,&tid->m_t_mtx);
      }
    }
#endif
  }
}

void Amplitude::SetCouplings() const
{
#ifdef DEBUG__CF
  msg_Debugging()<<METHOD<<"(): {\n";
#endif
  for (size_t i(0);i<m_cpls.size();++i) {
    MODEL::Coupling_Data *aqcd(m_cpls[i].p_aqcd);
    MODEL::Coupling_Data *aqed(m_cpls[i].p_aqed);
    double gsfac(aqcd?sqrt(aqcd->Factor()):1.0);
    double gwfac(aqed?sqrt(aqed->Factor()):1.0);
    double fac(1.0);
    Vertex *v(m_cpls[i].p_v);
    size_t oqcd(m_cpls[i].m_oqcd), oew(m_cpls[i].m_oew);
    if (aqcd && oqcd) {
#ifdef DEBUG__CF
      msg_Debugging()<<"  qcd: "<<aqcd->Factor()<<" ^ "<<oqcd/2.
		     <<" = "<<intpow(gsfac,oqcd)<<"\n";
#endif
      fac*=intpow(gsfac,oqcd);
    }
    if (aqed && oew) {
#ifdef DEBUG__CF
      msg_Debugging()<<"  qed: "<<aqed->Factor()<<" ^ "<<oew/2.
		     <<" = "<<intpow(gwfac,oew)<<"\n";
#endif
      fac*=intpow(gwfac,oew);
    }
    v->SetCplFac(fac);
  }
#ifdef DEBUG__CF
  msg_Debugging()<<"}\n";
#endif
}

void Amplitude::ResetJ()
{
  for (size_t n(1);n<=m_n-1;--n) {
    for (size_t i(0);i<m_cur[n].size();++i) 
      m_cur[n][i]->ResetJ();
  }
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
  Vec4D sum;
#endif
  for (size_t i(0);i<m_n;++i) {
    m_p[i]=m_dirs[i]>0?-moms[i]:moms[i];
#ifdef DEBUG__BG
    sum+=m_p[i];
    msg_Debugging()<<"set p["<<i<<"] = "<<m_p[i]
		   <<" ("<<sqrt(dabs(m_p[i].Abs2()))<<")\n";
#endif
  }
#ifdef DEBUG__BG
  static double accu(sqrt(rpa->gen.Accu()));
  if (!IsEqual(sum,Vec4D(),accu)) 
    msg_Error()<<METHOD<<"(): Four momentum not conserved. sum = "
  	       <<sum<<"."<<std::endl;
#endif
  if (m_subs.empty()) return true;
  p_dinfo->SetStat(1);
  for (size_t i(0);i<m_cur[1].size();++i) m_cur[1][i]->SetP(m_p[i]);
  for (size_t i(0);i<m_scur.size();++i)
    m_scur[i]->Sub()->In().front()->Kin()->Evaluate();
  return p_dinfo->Stat();
}

bool Amplitude::RSTrigger
(PHASIC::Combined_Selector *const sel,const int mode)
{
  if (m_subs.empty() || sel==NULL) return true;
  sel->RSTrigger(&m_subs);
  int trig=m_trig=m_subs.back()->m_trig;
  for (size_t i(0);i<m_scur.size();++i) {
    Dipole_Kinematics *kin(m_scur[i]->Sub()->In().front()->Kin());
    int ltrig(m_subs[i]->m_trig);
    kin->SetF(1.0);
    if (m_smth) {
      double a(m_smth>0.0?m_subs[i]->m_kt2=kin->KT2():kin->Y());
      if (a>0.0 && a<dabs(m_smth)) {
        kin->SetF(pow(a/dabs(m_smth),m_smpow));
        if (ltrig==0) kin->SetF(-kin->F());
        ltrig=1;
      }
    }
    kin->AddTrig(ltrig?1:0);
    m_subs[i]->m_trig=kin->Trig()?ltrig:0;
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
    double kt2j(sub->m_mu2[stp::res]);
    if (sub->m_mu2[stp::size+stp::res])
      kt2j=sub->m_mu2[stp::size+stp::res];
    int da(kin->A()>0.0 && kin->KT2()<kt2j);
    int ds(kin->Y()<p_dinfo->AMax(kin->Type()));
    kin->AddTrig(abs(ds-da));
    return ds-da;
  }
  if (mode==2) {
    double kt2j(sub->m_mu2[stp::res]);
    if (sub->m_mu2[stp::size+stp::res])
      kt2j=sub->m_mu2[stp::size+stp::res];
    int da(kin->A()>0.0 && kin->KT2()<kt2j);
    kin->AddTrig(da);
    return da;
  }
  if (mode==3) {
    int ds(kin->Y()<p_dinfo->AMax(kin->Type()));
    kin->AddTrig(abs(ds-1));
    return ds-1;
  }
  THROW(not_implemented,"Invalid call");
  return 0.0;
}

void Amplitude::SetCTS(void *const cts)
{
#ifdef USING__Threading
  p_cts=(CDBG_ME_TID_Vector*)cts;
#endif
}

void Amplitude::SetColors(const Int_Vector &rc,
			  const Int_Vector &ac,const int set)
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
  if (set==1) {
    Color_Calculator::SetCIMin(1);
    Color_Calculator::SetCIMax(0);
  }
  else if (set==-1) {
    int cimin(std::numeric_limits<int>::max()), cimax(0);
    for (size_t i(0);i<m_n;++i)
      if (m_cl[i][0]) {
	if (m_cl[i][0]<cimin) cimin=m_cl[i][0];
	if (m_cl[i][0]>cimax) cimax=m_cl[i][0];
      }
    msg_Debugging()<<"cimin = "<<cimin
		   <<", cimax = "<<cimax<<"\n";
    Color_Calculator::SetCIMin(cimin);
    Color_Calculator::SetCIMax(cimax);
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

bool Amplitude::EvaluateAll(const bool& mode)
{
  Complex total(Complex(0,0));
  if (p_loop) p_dinfo->SetDRMode(p_loop->DRMode());
  for (size_t i(0);i<m_subs.size();++i) m_subs[i]->Reset(0);
  for (size_t j(0);j<m_n;++j) m_ch[j]=0;
  MODEL::Coupling_Data *cpl(m_cpls.front().p_aqcd);
  // if (m_stype==1) cpl=m_cpls.front().p_aqed;
  double mu2(cpl?cpl->Scale():-1.0);
  p_dinfo->SetMu2(mu2);
  CalcJL();
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"(): "<<m_ress.size()<<" amplitudes {\n";
#endif
  if (m_pmode=='D') {
    for (size_t i(0);i<m_ress.size();++i) {
      if (m_cur.back()[i]->Sub()) continue;
      for (size_t j(0);j<m_ress[i].size();++j) m_ress[i][j]=0.0;
      m_cur.back()[i]->Contract(*m_cur[1].front(),m_cchirs,m_ress[i]);
#ifdef DEBUG__BG
      for (size_t j(0);j<m_ress[i].size();++j)
	msg_Debugging()<<"A["<<i<<"]["<<j<<"]"<<m_ress[i](j)<<" = "
		       <<m_ress[i][j]<<" -> "<<std::abs(m_ress[i][j])<<"\n";
#endif
    }
    for (size_t i(0);i<m_ress.size();++i) {
      if (m_cur.back()[i]->Sub()==NULL) continue;
      for (size_t j(0);j<m_ress[i].size();++j) m_cress[i][j]=m_ress[i][j]=0.0;
      if (p_dinfo->Mode()==1)
	m_cur.back()[i]->Contract(*m_cur.back()[i]->Sub(),m_cchirs,m_ress[i]);
      else {
	for (size_t j(0);j<m_ress[i].size();++j)
	  m_ress[i][j]=m_dsf[i].second*m_ress[m_dsf[i].first][j];
#ifdef DEBUG__BGS_AMAP
	Spin_Structure<DComplex> ress(m_fl,0.0);
	m_cur.back()[i]->Contract(*m_cur.back()[i]->Sub(),m_cchirs,ress);
	for (size_t j(0);j<m_ress[i].size();++j)
	  if (m_ress[i][j].real() || m_ress[i][j].imag()) {
	    msg_Debugging()<<"Check amplitude mapping ("<<m_dsf[i].first<<","
			   <<m_dsf[i].second<<"): "<<m_ress[i][j]<<" vs. "
			   <<ress[j]<<" -> "<<m_ress[i][j]/ress[j]<<"\n";
	    if (!IsEqual(m_ress[i][j].real(),ress[j].real(),1.0e-6)) {
	      msg_Error()<<METHOD<<"(): Mapping error in";
	      for (size_t l(0);l<m_nin;++l) msg_Error()<<" "<<m_fl[l];
	      msg_Error()<<" ->";
	      for (size_t l(m_nin);l<m_nin+m_nout;++l) msg_Error()<<" "<<m_fl[l];
	      msg_Error()<<" ("<<m_cur.back()[i]->Sub()->Sub()->Id().front()<<","
			 <<m_cur.back()[i]->Sub()->Id().front()<<")["<<j<<"]: "
			 <<m_ress[i][j]<<" / "<<ress[j]<<" = "<<m_ress[i][j]/ress[j]
			 <<" ("<<m_dsf[i].first<<","<<m_dsf[i].second<<")\n";
	    }
	  }
#endif
      }
      m_cur.back()[i]->Contract(*m_cur.back()[i]->Sub(),m_cchirs,m_cress[i],1);
      m_cur.back()[i]->Contract(*m_cur.back()[i]->Sub(),m_cchirs,m_cress[i],2);
#ifdef DEBUG__BG
      for (size_t j(0);j<m_ress[i].size();++j)
	msg_Debugging()<<"A["<<i<<"]["<<j<<"]"<<m_ress[i](j)<<" = "
		       <<m_ress[i][j]<<" * "<<m_cress[i][j]<<" -> "
		       <<m_ress[i][j]*std::conj(m_cress[i][j])<<"\n";
#endif
    }
  }
  else {
    THROW(not_implemented,"Internal error");
  }
#ifdef DEBUG__BG
  msg_Debugging()<<"}\n";
#endif
  double csum(0.0);
  for (size_t k(0);k<m_on.size();++k) {
    int i(m_on[k].first), j(m_on[k].second);
#ifdef DEBUG__BG
    std::vector<int> cpls(m_cur.back()[i]->Order());
    if (cpls.size()<m_cur.back()[j]->Order().size())
      cpls.resize(m_cur.back()[j]->Order().size(),0);
    for (size_t l(0);l<m_cur.back()[j]->Order().size();++l)
      cpls[l]+=m_cur.back()[j]->Order()[l];
    msg_Debugging()<<"ME^2 "<<cpls<<" from "
		   <<m_cur.back()[i]->Order()<<"x"
		   <<m_cur.back()[j]->Order()<<" {\n";
#endif
    for (size_t k(0);k<m_ress[i].size();++k)
      if (m_ress[i][k]!=Complex(0.0,0.0) &&
	  m_ress[j][k]!=Complex(0.0,0.0)) {
#ifdef DEBUG__BG
	msg_Debugging()<<"  A"<<m_ress[i](k)<<" = "
		       <<m_ress[i][k]<<" * "<<m_ress[j][k]<<" -> "
		       <<(m_ress[i][k]*std::conj(m_ress[j][k])).real()<<"\n";
#endif
	csum+=(m_ress[i][k]*std::conj(m_ress[j][k])).real();
        std::cout << "Component " << k << " is non-zero: " << m_ress[i][k] << std::endl;
        total += m_ress[i][k] * conj(m_ress[i][k]);
      }
#ifdef DEBUG__BG
    msg_Debugging()<<"}\n";
#endif
  }
  //std::cout << "Total-Sherpa" << total << std::endl;
  // Helizitätsmittelung
  //std::cout << "Total-Sherpa" << total/m_sf << std::endl;
  Complex s = (m_p[0]+m_p[1])*(m_p[0]+m_p[1]);
  Complex beta = sqrt(1.0-((4*m_p[2].Abs2())/s));
  Vec4<Complex> x = (1.0 / (4.0*beta))*((1.0-sqr(beta))*m_p[3]-sqr(1.0-beta)*m_p[2]);
  Vec4<Complex> y = (1.0 / (4.0*beta))*((1.0-sqr(beta))*m_p[2]-sqr(1.0-beta)*m_p[3]);
  Complex m_squared_massive = (1.0/4.0) * (sqr((8.0*Complex(0.30827124008247381, 0)*Complex(0.30827124008247381, 0)) / 3.0))
                              * (s/(m_p[2].Abs2()*m_p[2].Abs2()*m_p[2].Abs2())) * sqr(1.0+beta) * ((1.0-6.0*m_p[2].Abs2()/s+18.0*(sqr(m_p[2].Abs2())/sqr(s))) * (m_p[0] * y) * (m_p[1] * y))
                              + (s /(8.0*sqr(m_p[2].Abs2())))*(1.0-4.0*(m_p[2].Abs2()/s)+10.0*(sqr(m_p[2].Abs2())/sqr(s))) * ((m_p[0] * y)*(m_p[1] * x) + (m_p[0]*x)*(m_p[1]*y));
  if ((m_res - m_squared_massive).real() > 1e-12 || (m_res - m_squared_massive).imag() > 1e-12)
    std::cout << "Large diff between literature " << m_squared_massive << " and Sherpa " << total/m_sf << std::endl;
  else std::cout << "passed" << std::endl;
  m_born=m_res=csum/m_sf;
  m_cmur[1]=m_cmur[0]=csum=0.0;
  if (p_dinfo->Mode()) {
    assert(cpl!=NULL);
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
      if (!m_trig) m_born=m_res=0.0;
      m_subs.back()->m_me=m_subs.back()->m_mewgt=m_res;
    }
    if (p_dinfo->Mode()&2) {
      for (size_t i(0);i<m_dsij.size();++i)
	for (size_t j(0);j<m_dsij[i].size();++j) m_dsij[i][j]=0.0;
      m_dsij[0][0]=m_res;
    }
    for (size_t k(0);k<m_son.size();++k) {
      int i(m_son[k].first), j(m_son[k].second);
#ifdef DEBUG__BG
      std::vector<int> cpls(m_cur.back()[i]->Order());
      if (cpls.size()<m_cur.back()[j]->Order().size())
	cpls.resize(m_cur.back()[j]->Order().size(),0);
      for (size_t l(0);l<m_cur.back()[j]->Order().size();++l)
	cpls[l]+=m_cur.back()[j]->Order()[l];
      msg_Debugging()<<"ME^2 "<<cpls<<" from "
		     <<m_cur.back()[i]->Order()<<"x"
		     <<m_cur.back()[j]->Order()<<", dipole "
		     <<m_cur.back()[i]->Sub()->Sub()->Id()
		     <<m_cur.back()[i]->Sub()->Id()<<" (type="
		     <<m_cur.back()[i]->Sub()->Sub()->In().front()->SType()<<") {\n";
#endif
      double ccsum(0.0);
      for (size_t k(0);k<m_ress[i].size();++k)
	if (m_ress[i][k]!=Complex(0.0,0.0) &&
	    m_ress[j][k]!=Complex(0.0,0.0)) {
#ifdef DEBUG__BG
	  msg_Debugging()<<"  A["<<i<<"/"<<j<<"]"<<m_ress[i](k)<<" = "
			 <<m_ress[i][k]<<" * "<<m_cress[j][k]<<" -> "
			 <<m_ress[i][k]*std::conj(m_cress[j][k])<<"\n";
#endif
	ccsum+=(m_ress[i][k]*std::conj(m_cress[j][k])).real();
      }
#ifdef DEBUG__BG
      msg_Debugging()<<"}\nccsum = "<<ccsum<<"\n";
#endif
      ccsum/=m_sf;
      if (p_dinfo->Mode()==1) {
	if (m_smth) {
	  Dipole_Kinematics *kin=m_cur.back()[i]->
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
	m_subs[m_sid[i]]->m_me=m_subs[m_sid[i]]->m_mewgt=ccsum;
      }
      else {
	Dipole_Kinematics *kin=m_cur.back()[j]->
	  Sub()->Sub()->In().front()->Kin();
	m_p[0]=-m_p[0];
	m_p[1]=-m_p[1];
	//double lf(log(2.0*M_PI*mu2/EpsSchemeFactor(m_p)/
	//	      dabs(kin->JIJT()->P()*kin->JK()->P())));
  
  double lf(0.);
  if (!p_loop || !(p_loop->fixedIRscale())) 
      lf = log(2.0*M_PI*mu2/EpsSchemeFactor(m_p)/dabs(kin->JIJT()->P()*kin->JK()->P()));
  else{
      double irscale=p_loop->IRscale();
      lf = log(2.0*M_PI*sqr(irscale)/EpsSchemeFactor(m_p)/dabs(kin->JIJT()->P()*kin->JK()->P()));
  }

	m_p[0]=-m_p[0];
	m_p[1]=-m_p[1];
#ifdef DEBUG__BG
	msg_Debugging()<<"e^2 = "<<kin->Res(2)<<", e = "<<kin->Res(1)
		       <<", f = "<<kin->Res(0)<<", l = "<<lf
		       <<" ( "<<p_loop<<" )\n";
#endif
	m_dsij[m_dsm[j].first][m_dsm[j].second]+=ccsum;
	m_dsij[m_dsm[j].second][m_dsm[j].first]+=ccsum;
	m_cmur[1]-=ccsum*asf*kin->Res(2);
	m_cmur[0]-=ccsum*asf*(kin->Res(1)+lf*kin->Res(2));
	ccsum*=-asf*(kin->Res(0)+lf*kin->Res(1)+0.5*sqr(lf)*kin->Res(2));
      }
      if (mode || p_dinfo->IType()&cs_itype::I) csum+=ccsum;
    }
    if (p_loop) {
      double cw(p_loop->Mode()?1.0/m_sf:m_res);
      if (p_loop->Mode() && p_loop->ColMode()==0)
         cw*=1.0/p_colint->GlobalWeight();
      if (p_dinfo->Mode()&8) {
	double e1p(-m_cmur[0]/m_res/asf), e2p(-m_cmur[1]/m_res/asf);
	double e1l(p_loop->ME_E1()), e2l(p_loop->ME_E2());
	if (p_loop->Mode()) {
	  e1l/=p_loop->ME_Born()*m_sf;
	  e2l/=p_loop->ME_Born()*m_sf;
	}
	if (!IsEqual(e2p,e2l))
	  msg_Error()<<METHOD<<"(): Double pole does not match. V -> "
		     <<e2l<<", I -> "<<e2p<<", rel. diff. "
		     <<(e2p/e2l-1.0)<<".\n";
	if (!IsEqual(e1p,e1l))
	  msg_Error()<<METHOD<<"(): Single pole does not match. V -> "
		     <<e1l<<", I -> "<<e1p<<", rel. diff. "
		     <<(e1p/e1l-1.0)<<".\n";
      }
      csum+=cw*asf*p_loop->ME_Finite();
      if (m_murcoeffvirt && p_loop->ProvidesPoles()) {
        if (m_sccmur) {
          double b0(11.0/6.0*3.0-2.0/3.0*0.5*Flavour(kf_quark).Size()/2.);
          m_cmur[0]+=cw*asf*(p_loop->ME_E1()+m_maxcpl[0]/2.*b0);
          m_cmur[1]+=cw*asf*p_loop->ME_E2();
        }
        else {
          m_cmur[0]+=cw*asf*p_loop->ScaleDependenceCoefficient(1);
          m_cmur[1]+=cw*asf*p_loop->ScaleDependenceCoefficient(2);
        }
      }
      else {
	double b0(11.0/6.0*3.0-2.0/3.0*0.5*Flavour(kf_quark).Size()/2.);
        m_cmur[0]=cw*asf*m_maxcpl[0]/2.*b0;
        m_cmur[1]=0.;
      }
    }
    if (!(p_dinfo->Mode()&4)) {
      if (p_dinfo->Mode()&2) m_born=m_res=0.0;
      if ((p_dinfo->Mode()&16) && p_loop) m_born=m_res=0.0;
    }
  }
#ifdef DEBUG__BG
  msg_Debugging()<<"m_res = "<<m_res<<", csum = "<<csum
		 <<" -> "<<m_res+csum<<"\n";
#endif
  m_res+=csum;
  ResetJ();
  return true;
}

double Amplitude::Differential(NLO_subevt *const sub)
{
  p_sub=sub;
  return Differential(p_colint->I(),p_colint->J());
}

double Amplitude::Differential
(const Int_Vector &ci,const Int_Vector &cj,const int set)
{
  SetColors(ci,cj,set);
  if (p_helint!=NULL && p_helint->On()) {
    Evaluate(p_helint->Chiralities());
    return m_res;
  }
  Spinor<double> p(1, m_p[0]);
  Spinor<double> q(1, m_p[1]);
  Spinor<double> k(1, m_p[2]);
  Spinor<double> l(1, m_p[3]);
  Spinor<double> g(1, Vec4D(1.0,0.0,1.0,0.0));
  Spinor<double> h(1, Vec4D(1.0,0.0,1.0,0.0));
  /*Spinor<double> p(1, -(m_p[0][1]+Complex(0, -1)*m_p[0][2]) / (sqrt(m_p[0][0]+m_p[0][3])), sqrt(m_p[0][0]+m_p[0][3]));
  Spinor<double> q(1, -(m_p[1][1]+Complex(0, -1)*m_p[1][2]) / (sqrt(m_p[1][0]+m_p[1][3])), sqrt(m_p[1][0]+m_p[1][3]));
  Spinor<double> k(1, -(m_p[2][1]+Complex(0, -1)*m_p[2][2]) / (sqrt(m_p[2][0]+m_p[2][3])), sqrt(m_p[2][0]+m_p[2][3]));
  Spinor<double> l(1, -(m_p[3][1]+Complex(0, -1)*m_p[3][2]) / (sqrt(m_p[3][0]+m_p[3][3])), sqrt(m_p[3][0]+m_p[3][3]));
  Vec4D gauge = Vec4D(1.0,0.0,invsqrttwo,invsqrttwo);
  Spinor<double> g(1, -(gauge[1]+Complex(0, -1)*gauge[2]) / (sqrt(gauge[0]+gauge[3])), sqrt(gauge[0]+gauge[3]));
  Spinor<double> h(1, -(gauge[1]+Complex(0, -1)*gauge[2]) / (sqrt(gauge[0]+gauge[3])), sqrt(gauge[0]+gauge[3]));*/
  Spinor<double> p_up(1, p[1], -p[0]);
  Spinor<double> q_up(1, q[1], -q[0]);
  Spinor<double> k_up(1, k[1], -k[0]);
  Spinor<double> l_up(1, l[1], -l[0]);
  Spinor<double> g_up(1, g[1], -g[0]);
  Spinor<double> h_up(1, h[1], -h[0]);
  //Complex m_ppmm = (2.0*Complex(0, 1)*0.00756234/((m_p[0]+m_p[1])*(m_p[0]+m_p[1])) * (conj(k*p)*(l*q)*conj(h*k)*(l*g)) / (conj(g*l)*(h*k)));
  //Complex m_mmpp = (2.0*Complex(0, 1)*0.00756234/((m_p[0]+m_p[1])*(m_p[0]+m_p[1]))) * ((k*q)*conj(l*p)*(h*k)*conj(l*g)) / ((g*l)*conj(h*k));
  //std::cout << "Literature +,-,++,-- / (-,+,--,++)^dagger matrix element" << m_ppmm << std::endl;
  //std::cout << "Literature +,-,--,++ / (-,+,++,--)^dagger matrix element"<< m_mmpp << std::endl;
  // zwei Komponenten sind gleich modulo dagger, daher Faktor 2, aber 1/4 wegen Helizitätsmittelung, daher nur 1/2
  //std::cout << "Total " << std::complex<double>(0.5,0)*(m_ppmm*conj(m_ppmm)+m_mmpp*conj(m_mmpp))  << std::endl;
  //std::cout << "Total-2 " << ((32.0*sqr(0.30827124008247381)*sqr(0.30827124008247381))/(sqr((m_p[0]+m_p[1])*(m_p[0]+m_p[1])))) * ((m_p[0]*m_p[3])*(m_p[1]*m_p[2])+(m_p[0]*m_p[2])*(m_p[1]*m_p[3])) << std::endl;
  EvaluateAll();
  //std::cout << m_res << std::endl;
  std::cout << METHOD << "Compare literature and Comix result: " << std::endl;
  if (m_p[2].Mass()<1e-6){
    std::cout << "hier komme ich hin" << std::endl;
    Complex m_ppmm = ((2.0*Complex(0, 1)*Complex(0.30827124008247381, 0)*Complex(0.30827124008247381, 0))/((m_p[0]+m_p[1])*(m_p[0]+m_p[1])) * (conj(k_up*p)*(l_up*q)*conj(h_up*k)*(l_up*g)) / (conj(h_up*k)*(l_up*g)));
    Complex m_mmpp = ((2.0*Complex(0, 1)*Complex(0.30827124008247381, 0)*Complex(0.30827124008247381, 0))/((m_p[0]+m_p[1])*(m_p[0]+m_p[1]))) * ((k_up*q)*conj(l_up*p)*(h_up*k)*conj(l_up*g)) / ((h_up*k)*conj(l_up*g));
    if ((m_res - (std::complex<double>(0.5,0)*(m_ppmm*conj(m_ppmm)+m_mmpp*conj(m_mmpp)))).real() > 1e-12 ||
        (m_res - (std::complex<double>(0.5,0)*(m_ppmm*conj(m_ppmm)+m_mmpp*conj(m_mmpp)))).imag() > 1e-12)
      std::cout << "Large diff between literature " << std::complex<double>(0.5,0)*(m_ppmm*conj(m_ppmm)+m_mmpp*conj(m_mmpp))
                << " and Sherpa " << m_res << std::endl;
    else std::cout << "passed" << std::endl;
  }
  else{
    Complex s = (m_p[0]+m_p[1])*(m_p[0]+m_p[1]);
    Complex beta = sqrt(1.0-((4*m_p[2].Abs2())/s));
    Vec4<Complex> x = (1.0 / (4.0*beta))*((1.0-sqr(beta))*m_p[3]-sqr(1.0-beta)*m_p[2]);
    Vec4<Complex> y = (1.0 / (4.0*beta))*((1.0-sqr(beta))*m_p[2]-sqr(1.0-beta)*m_p[3]);
    Complex m_squared_massive = (1.0/4.0) * (sqr((8.0*Complex(0.30827124008247381, 0)*Complex(0.30827124008247381, 0)) / 3.0))
      * (s/(m_p[2].Abs2()*m_p[2].Abs2()*m_p[2].Abs2())) * sqr(1.0+beta) * ((1.0-6.0*m_p[2].Abs2()/s+18.0*(sqr(m_p[2].Abs2())/sqr(s))) * (m_p[0] * y) * (m_p[1] * y))
      + (s /(8.0*sqr(m_p[2].Abs2())))*(1.0-4.0*(m_p[2].Abs2()/s)+10.0*(sqr(m_p[2].Abs2())/sqr(s))) * ((m_p[0] * y)*(m_p[1] * x) + (m_p[0]*x)*(m_p[1]*y));
    if ((m_res - m_squared_massive).real() > 1e-12 || (m_res - m_squared_massive).imag() > 1e-12)
      std::cout << "Large diff between literature " << m_squared_massive << " and Sherpa " << m_res << std::endl;
    else std::cout << "passed" << std::endl;
  }

  return m_res;
}

bool Amplitude::ConstructChirs()
{
  DEBUG_FUNC("");
  for (size_t i(0);i<m_cur.back().size();++i) {
    m_ress.push_back(Spin_Structure<DComplex>(m_fl,0.0));
    m_cur.back()[i]->HM().resize(m_ress.back().size());
    for (size_t j(0);j<m_ress.back().size();++j)
      m_cur.back()[i]->HM()[j]=j;
  }
  m_cchirs.resize(m_ress.back().size());
  for (size_t j(0);j<m_ress.back().size();++j) m_cchirs[j]=1;
  if (p_dinfo->Mode()==0) return true;
  m_cress=m_ress;
  for (size_t i(0);i<m_scur.size();++i) {
    Dipole_Kinematics *kin(m_scur[i]->Sub()->In().front()->Kin());
    Current *last(NULL);
    for (size_t j(0);j<m_cur.back().size();++j) {
      if (m_cur.back()[j]->Sub()!=m_scur[i]) continue;
      last=m_cur.back()[j];
      last->HM().resize(m_ress.front().size(),0);
      kin->PM().resize(m_ress.front().size(),0);
      int ii(kin->JI()?kin->JI()->Id().front():-1);
      int ij(kin->JJ()?kin->JJ()->Id().front():-1);
      size_t ik(kin->JK()->Id().front());
      Flavour_Vector cfl(1,m_scur[i]->Flav());
      for (size_t j(0);j<m_n;++j) if (j!=ik) cfl.push_back(m_fl[j]);
      Spin_Structure<int> tch(cfl,0.0);
      for (size_t j(0);j<m_ress.front().size();++j) {
	Int_Vector id(m_ress.front().GetSpinCombination(j)), di(1,id[ik]);
	for (size_t k(0);k<m_n;++k)
	  if (k!=ik) di.push_back(id[k]);
	last->HM()[tch.GetNumber(di)]=j;
	if (ii==-1 || ij==-1 || id[ii]==id[ij]) continue;
	std::swap<int>(id[ii],id[ij]);
	kin->PM()[j]=m_ress.front().GetNumber(id);
      }
    }
  }
  return true;
}

bool Amplitude::CheckOrders()
{
  std::vector<int> maxcpl(2,0), mincpl(2,99);
  std::vector<int> maxacpl(2,0), minacpl(2,99);
  msg_Debugging()<<m_cur.back().size()<<std::endl;
  for (size_t i(0);i<m_cur.back().size();++i) {
    bool any(true);
    const std::vector<int> &icpls(m_cur.back()[i]->Order());
    if (m_maxacpl.size()<icpls.size()) m_maxacpl.resize(icpls.size(),99);
    if (m_minacpl.size()<icpls.size()) m_minacpl.resize(icpls.size(),0);
    for (size_t k(0);k<icpls.size();++k)
      if (icpls[k]<m_minacpl[k] || icpls[k]>m_maxacpl[k]) any=false;
    if (any) {
      any=false;
      if (maxacpl.size()<icpls.size()) maxacpl.resize(icpls.size(),0);
      if (minacpl.size()<icpls.size()) minacpl.resize(icpls.size(),99);
      for (size_t k(0);k<icpls.size();++k) {
	maxacpl[k]=Max(maxacpl[k],icpls[k]);
	minacpl[k]=Min(minacpl[k],icpls[k]);
      }
    for (size_t j(0);j<m_cur.back().size();++j) {
      const std::vector<int> &jcpls(m_cur.back()[j]->Order());
      if (m_cur.back()[i]->Sub()!=
	  m_cur.back()[j]->Sub()) continue;
      bool on(true);
      std::vector<int> cpls(icpls);
      if (cpls.size()<jcpls.size()) cpls.resize(jcpls.size(),0);
      for (size_t k(0);k<jcpls.size();++k) cpls[k]+=jcpls[k];
      if (m_maxcpl.size()<cpls.size()) m_maxcpl.resize(cpls.size(),99);
      if (m_mincpl.size()<cpls.size()) m_mincpl.resize(cpls.size(),0);
      for (size_t k(0);k<cpls.size();++k)
        if (cpls[k]<m_mincpl[k] || cpls[k]>m_maxcpl[k]) on=false;
      for (size_t k(cpls.size());k<m_mincpl.size();++k)
        if (m_mincpl[k]) on=false;
      if (on) {
	if (maxcpl.size()<cpls.size()) maxcpl.resize(cpls.size(),0);
	if (mincpl.size()<cpls.size()) mincpl.resize(cpls.size(),99);
	for (size_t k(0);k<cpls.size();++k) {
	  maxcpl[k]=Max(maxcpl[k],cpls[k]);
	  mincpl[k]=Min(mincpl[k],cpls[k]);
	}
	if (m_cur.back()[i]->Sub())
	  m_son.push_back(std::pair<size_t,size_t>(i,j));
	else m_on.push_back(std::pair<size_t,size_t>(i,j));
	any=true;
      }
    }
    }
    if (!any) {
#ifdef DEBUG__BG
      msg_Debugging()<<"delete obsolete current: "<<m_cur.back()[i]<<"\n";
      m_cur.back()[i]->Print();
#endif
      delete m_cur.back()[i];
      m_cur.back().erase(m_cur.back().begin()+i);
      for (size_t k(0);k<m_son.size();++k)
	if (m_son[k].second>i) --m_son[k].second;
      for (size_t k(0);k<m_on.size();++k)
	if (m_on[k].second>i) --m_on[k].second;
      --i;
    }
  }
  bool valid(false);
  for (size_t k(0);k<maxcpl.size();++k)
    if (maxcpl[k]) valid=true;
  if (!valid) return false;
  m_maxcpl=maxcpl;
  m_mincpl=mincpl;
  m_maxacpl=maxacpl;
  m_minacpl=minacpl;
  if (m_cur.back().empty()) return false;
  Prune();
  return true;
}

bool Amplitude::ConstructCouplings(MODEL::Coupling_Map *const cpls)
{
  DEBUG_FUNC("");
  MODEL::Coupling_Data *rqcd(cpls->Get("Alpha_QCD"));
  MODEL::Coupling_Data *rqed(cpls->Get("Alpha_QED"));
  for (long int i(0);i<m_cur.back().size();++i) {
    for (size_t j(0);j<m_cur.back()[i]->In().size();++j) {
      Vertex *v(m_cur.back()[i]->In()[j]);
      int oqcd(v->Order(0)), oew(v->Order(1));
      for (size_t i(0);i<v->J().size();++i) {
	oqcd+=v->J(i)->Order(0);
	oew+=v->J(i)->Order(1);
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
      int oqcd(v->Order(0)), oew(v->Order(1));
      for (size_t i(0);i<v->J().size();++i) {
	oqcd+=v->J(i)->Order(0);
	oew+=v->J(i)->Order(1);
      }
      m_cpls.push_back(Coupling_Info(v,oqcd,oew,aqcd,aqed));
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
  if (!CheckOrders()) return false;
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
  if (!ConstructChirs()) return false;
  m_sid.resize(m_cur.back().size(),0);
  if (p_dinfo->Mode()==1) ConstructNLOEvents();
  if (p_dinfo->Mode()&2) ConstructDSijMap();
  if (!ConstructCouplings(cpls)) return false;
  return true;
}

void Amplitude::FillAmplitudes
(std::vector<Spin_Amplitudes> &amps,
 std::vector<std::vector<Complex> > &cols)
{
  cols.push_back(std::vector<Complex>(1,1.0));
  amps.push_back(Spin_Amplitudes(m_fl,Complex(0.0,0.0)));
  for (size_t i(0);i<m_ress.front().size();++i) {
    Complex sum(m_ress.front()[i]);
    for (size_t j(1);j<m_ress.size();++j) sum+=m_ress[j][i];
    amps.back().Insert(sum,m_ress.front()(i));
  }
}

double Amplitude::Coupling(const int mode) const
{
  MODEL::Coupling_Data *cpl(m_cpls.front().p_aqcd);
  return cpl->Default()*cpl->Factor();
}

void Amplitude::FillMEWeights(ME_Weight_Info &wgtinfo) const
{
  if (wgtinfo.m_wren.size()<2) return;
  for (size_t i=0;i<2;i++) wgtinfo.m_wren[i]=m_cmur[i];
  wgtinfo.m_VI=m_res-m_born;
}

void Amplitude::SetNLOMC(PDF::NLOMC_Base *const mc)
{
  for (size_t i(0);i<m_scur.size();++i)
    m_scur[i]->Sub()->In().front()->Kin()->SetNLOMC(mc);
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

//TODO: Include Gauge(1) into calculations to enable gauge test -> Ask Stefan which gamma matrices he has used
//      Fail of gauge test and different cross section expected because only one gamma matrix implemented
bool Amplitude::GaugeTest(const Vec4D_Vector &moms,const int mode)
{
  Spinor<double> p(1, moms[0]);
  Spinor<double> q(1, moms[1]);
  Spinor<double> k(1, moms[2]);
  Spinor<double> l(1, moms[3]);
  Spinor<double> g(1, Vec4D(1.0,0.0,1.0,0.0));
  Spinor<double> h(1, Vec4D(1.0,0.0,1.0,0.0));
  /*Spinor<double> p(1, -(moms[0][1]+Complex(0, -1)*moms[0][2]) / (sqrt(moms[0][0]+moms[0][3])), sqrt(moms[0][0]+moms[0][3]));
  Spinor<double> q(1, -(moms[1][1]+Complex(0, -1)*moms[1][2]) / (sqrt(moms[1][0]+moms[1][3])), sqrt(moms[1][0]+moms[1][3]));
  Spinor<double> k(1, -(moms[2][1]+Complex(0, -1)*moms[2][2]) / (sqrt(moms[2][0]+moms[2][3])), sqrt(moms[2][0]+moms[2][3]));
  Spinor<double> l(1, -(moms[3][1]+Complex(0, -1)*moms[3][2]) / (sqrt(moms[3][0]+moms[3][3])), sqrt(moms[3][0]+moms[3][3]));
  Vec4D gauge = Vec4D(1.0,0.0,invsqrttwo,invsqrttwo);
  Spinor<double> g(1, -(gauge[1]+Complex(0, -1)*gauge[2]) / (sqrt(gauge[0]+gauge[3])), sqrt(gauge[0]+gauge[3]));
  Spinor<double> h(1, -(gauge[1]+Complex(0, -1)*gauge[2]) / (sqrt(gauge[0]+gauge[3])), sqrt(gauge[0]+gauge[3]));*/
  Spinor<double> p_up(1, p[1], -p[0]);
  Spinor<double> q_up(1, q[1], -q[0]);
  Spinor<double> k_up(1, k[1], -k[0]);
  Spinor<double> l_up(1, l[1], -l[0]);
  Spinor<double> g_up(1, g[1], -g[0]);
  Spinor<double> h_up(1, h[1], -h[0]);
  std::cout << p << q << k << l << g << h << std::endl;
  Complex m_ppmm = ((2.0*Complex(0, 1)*Complex(0.30827124008247381, 0)*Complex(0.30827124008247381, 0))/((moms[0]+moms[1])*(moms[0]+moms[1])) * (conj(k_up*p)*(l_up*q)*conj(h_up*k)*(l_up*g)) / (conj(h_up*k)*(l_up*g)));
  Complex m_mmpp = ((2.0*Complex(0, 1)*Complex(0.30827124008247381, 0)*Complex(0.30827124008247381, 0))/((moms[0]+moms[1])*(moms[0]+moms[1]))) * ((k_up*q)*conj(l_up*p)*(h_up*k)*conj(l_up*g)) / ((h_up*k)*conj(l_up*g));
  //Complex m_ppmm = (2.0*Complex(0, 1)*0.00756234/((moms[0]+moms[1])*(moms[0]+moms[1])) * (conj(k*p)*(l*q)*conj(h*k)*(l*g)) / (conj(g*l)*(h*k)));
  //Complex m_mmpp = (2.0*Complex(0, 1)*0.00756234/((moms[0]+moms[1])*(moms[0]+moms[1]))) * ((k*q)*conj(l*p)*(h*k)*conj(l*g)) / ((g*l)*conj(h*k));
  if (moms[2].Mass()<1e-6 && moms[3].Mass()<1e-6){
    std::cout << "Literature +,-,++,-- / (-,+,--,++)^dagger matrix element" << m_ppmm << std::endl;
    std::cout << "Literature +,-,--,++ / (-,+,++,--)^dagger matrix element"<< m_mmpp << std::endl;
    // zwei Komponenten sind gleich modulo dagger, daher Faktor 2, aber 1/4 wegen Helizitätsmittelung, daher nur 1/2
    std::cout << "Total " << std::complex<double>(0.5,0)*(m_ppmm*conj(m_ppmm)+m_mmpp*conj(m_mmpp))  << std::endl;
    std::cout << "Total-2 " << ((32.0*sqr(0.30827124008247381)*sqr(0.30827124008247381))/(sqr((moms[0]+moms[1])*(moms[0]+moms[1])))) * ((moms[0]*moms[3])*(moms[1]*moms[2])+(moms[0]*moms[2])*(moms[1]*moms[3])) << std::endl;
  }
  if (mode==0) {
    size_t nt(0);
    bool cnt(true);
    while (cnt) {
      while (!p_colint->GeneratePoint());
      SetColors(p_colint->I(),p_colint->J());
      if (!GaugeTest(moms,1)) return false;
      double res(m_born?m_born:m_res);
      if (res!=0.0) cnt=false;
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
  PHASIC::Virtual_ME2_Base *loop(p_loop);
  p_loop=NULL;
  SetGauge(1);
  SetMomenta(moms);
  if (!EvaluateAll(true)) {
    p_loop=loop;
    return false;
  }
  double res(m_born?m_born:m_res);
  if (m_pmode=='D') Spinor<double>::ResetGauge();
  SetGauge(0);
  SetMomenta(moms);
  if (!EvaluateAll(true)) {
    p_loop=loop;
    return false;
  }
  p_loop=loop;
  double res2(m_born?m_born:m_res);
  msg_Debugging()<<METHOD<<"(): {\n";
  msg_Debugging()<<"  \\sigma_{tot} = "<<res<<" vs. "<<res2
                 <<" -> dev. "<<res2/res-1.0<<"\n";
  if (!IsEqual(res2,res)) {
    msg_Error().precision(12);
    msg_Error()<<"\n"<<METHOD<<"(): Large deviation {\n      "
               <<std::setw(18)<<std::right<<res2<<"\n   vs "
	       <<std::setw(18)<<res<<"\n   => "<<std::setw(18)
	       <<(res2/res-1.0)<<"\n}"<<std::left<<std::endl;
    msg_Error().precision(6);
    return true;
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
  for (size_t i(0);i<m_cur.back().size();++i)
    if (m_cur.back()[i]->Sub()==NULL) {
      Graph_Node graphs("j_1",true);
      graphs.push_back("    %% "+graphs.back());
      m_cur.back()[i]->CollectGraphs(&graphs);
      WriteOutGraph(str,&graphs,ng,cvs);
    }
  str<<"\n  \\end{longtable}\n\n";
  if (!cvs.empty() && (m_pgmode&1)) {
    str<<"  \\begin{longtable}{c}\n";
    for (std::set<std::string>::const_iterator
	   vit(cvs.begin());vit!=cvs.end();++vit) {
      std::string label(*vit);
      for (size_t pos(label.find(",,"));
	   pos!=std::string::npos;pos=label.find(",,",pos+1))
	label.replace(pos,2,",");
      str<<"    $"<<label<<"$\\\\\n";
    }
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
