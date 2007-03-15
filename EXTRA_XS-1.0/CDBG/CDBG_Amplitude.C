#include "CDBG_Amplitude.H"

#include "Message.H"
#include "Exception.H"
#include "MyStrStream.H"
#include "Vector_Current.H"
#include "Tensor_Current.H"
#include "STL_Tools.H"
#ifdef PROFILE__all
#include "prof.hh"
#else
#define PROFILE_HERE
#define PROFILE_LOCAL(NAME)
#endif

using namespace EXTRAXS;
using namespace ATOOLS;

// #define DEBUG__BG

#ifdef DEBUG__BG
static size_t s_nonzero(0), s_nonzeroall(0);
#endif

static const double sqrttwo(sqrt(2.0));

CDBG_Amplitude::CDBG_Amplitude():
  p_p(NULL), p_ssp(NULL),
  p_fl(NULL),
  p_sp(NULL), p_sm(NULL),
  p_ep(NULL), p_em(NULL), p_ssep(NULL), p_ssem(NULL),
  p_ch(NULL), p_cl(NULL),
  m_n(0),
  m_mode(0), m_nf(3)
{
  Vec4D k(1.0,0.0,1.0,0.0);
  m_kp=Spinor(1,k);
  m_km=Spinor(-1,k);
}

CDBG_Amplitude::~CDBG_Amplitude()
{
#ifdef DEBUG__BG
  msg_Info()<<"CDBG:     J!=0 -> "<<s_nonzero<<"\n";
  msg_Info()<<"CDBG: all J!=0 -> "<<s_nonzeroall<<"\n";
#endif
  CleanUp();
}

void CDBG_Amplitude::CleanUp()
{
  if (p_p!=NULL) delete p_p;
  if (p_sp!=NULL) delete p_sp;
  if (p_fl!=NULL) delete p_fl;
  if (p_sp!=NULL) delete p_sp;
  if (p_sm!=NULL) delete p_sm;
  if (p_ep!=NULL) delete p_ep;
  if (p_em!=NULL) delete p_em;
  if (p_ssep!=NULL) delete p_ssep;
  if (p_ssem!=NULL) delete p_ssem;
  if (p_ch!=NULL) delete p_ch;
  if (p_cl!=NULL) {
    for (size_t i(0);i<m_n;++i) delete p_cl[i];
    delete p_cl;
  }
  m_n=0;
  p_ssp=p_p=NULL;
  p_fl=NULL;
  p_sp=p_sm=NULL;
  p_ssem=p_ssep=p_ep=p_em=NULL;
  p_ch=NULL;
  p_cl=NULL;
  m_cur=Current_Matrix();
  m_chirs=Int_Matrix();
  m_ress.clear();
  m_hmap.clear();
  m_hamps.clear();
  m_maxid=0;
}

typedef std::pair<Current*,Current*> Current_Pair;

bool CDBG_Amplitude::MatchIndices(const Int_Vector &ids,const size_t &n,
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

void CDBG_Amplitude::AddCurrent(const std::vector<int> &ids,const size_t &n,
				const Flavour &fl)
{
  // add new current
  Current *cur;
  if (fl.IsFermion()) { THROW(not_implemented,"Flavour not available"); }
  else if (fl.IsVector()) cur = new Vector_Current(this,fl);
  else if (fl.IsTensor()) cur = new Tensor_Current(this,fl);
  else THROW(fatal_error,"Invalid flavour");
  cur->SetId(ids);
  cur->SetKey(m_cur[n].size());
  std::set<Current_Pair> v3;
  // compose current from all possible subcurrents
  for (size_t i(1);i<n;++i) {
    for (size_t j(0);j<m_cur[i].size();++j) {
      for (size_t k(0);k<m_cur[n-i].size();++k) {
	if (!MatchIndices(ids,n,i,j,k) || v3.find
	    (Current_Pair(m_cur[n-i][k],m_cur[i][j]))!=
	    v3.end()) continue;
	Vertex *v(cur->GetVertex(m_cur[i][j],m_cur[n-i][k]));
	if (v!=NULL) 
	  v3.insert(Current_Pair(m_cur[i][j],m_cur[n-i][k]));
      }
    }
  }
  if (!v3.empty() || n==1) {
    m_cur[n].push_back(cur);
    cur->Print();
  }
  else delete cur;
}

void CDBG_Amplitude::Construct(std::vector<int> ids,const size_t &n)
{
  if (ids.size()==n) {
    if (n==1) {
      AddCurrent(ids,n,p_fl[ids.back()]);
    }
    else {
      AddCurrent(ids,n,kf::gluon);
      // add pseudoparticles for 4-gluon vertex
      if (n>1 && n<m_n-1) AddCurrent(ids,n,kf::gluonqgc);
    }
    return;
  }
  // currents are unordered -> use ordered indexing
  size_t last(ids.empty()?-1:ids.back());
  ids.push_back(0);
  if (n==m_n-1) {
    // calculate only one final current
    ids.back()=last+1;
    Construct(ids,n);
    return;
  }
  for (size_t m(n>1?m_n-1:m_n), i(last+1);i<m;++i) {
    // fill currents 0..n-1 for external partons
    // currents 0..n-2 each index for internal partons
    ids.back()=i;
    Construct(ids,n);
  }
}

bool CDBG_Amplitude::Construct(const std::vector<ATOOLS::Flavour> &flavs)
{
  m_n=flavs.size();
  p_p = new Vec4D[m_n];
  p_ssp = new Vec4D[m_n];
  p_fl = new Flavour[m_n];
  for (size_t i(0);i<m_n;++i) p_fl[i]=flavs[i];
  p_sm = new CSpinor[m_n];
  p_sp = new CSpinor[m_n];
  p_em = new CVec4D[m_n];
  p_ep = new CVec4D[m_n];
  p_ssem = new CVec4D[m_n];
  p_ssep = new CVec4D[m_n];
  p_ch = new int[m_n];
  p_cl = new int*[m_n];
  for (size_t i(0);i<m_n;++i) p_cl[i] = new int[2];
  m_cur.resize(m_n);
  Int_Vector ids;
  size_t tsum(0);
  for (size_t i(1);i<m_n;++i) {
    Construct(ids,i);
    tsum+=m_cur[i].size();
    msg_Debugging()<<"P("<<i<<","<<m_n<<"): "<<m_cur[i].size()<<"\n";
  }
  msg_Debugging()<<"total: cur = "<<tsum<<"\n";
  return true;
}

CVec4D CDBG_Amplitude::EM(const Vec4D &p)
{
  Spinor pp(1,p);
  CVec4D e;
  e[0]=pp.U1()*m_km.U1()+pp.U2()*m_km.U2();
  e[3]=pp.U1()*m_km.U1()-pp.U2()*m_km.U2();
  e[1]=pp.U1()*m_km.U2()+pp.U2()*m_km.U1();
  e[2]=Complex(0.0,1.0)*(pp.U1()*m_km.U2()-pp.U2()*m_km.U1());
  return e/(sqrttwo*std::conj(m_kp*pp));
}

CVec4D CDBG_Amplitude::EP(const Vec4D &p)
{
  Spinor pm(-1,p);
  CVec4D e;
  e[0]=pm.U1()*m_kp.U1()+pm.U2()*m_kp.U2();
  e[3]=pm.U1()*m_kp.U1()-pm.U2()*m_kp.U2();
  e[1]=pm.U1()*m_kp.U2()+pm.U2()*m_kp.U1();
  e[2]=Complex(0.0,-1.0)*(pm.U1()*m_kp.U2()-pm.U2()*m_kp.U1());
  return e/(sqrttwo*std::conj(m_km*pm));
}

namespace ATOOLS {

  bool operator<(const Complex &a,const Complex &b)
  {
    if (a.real()<b.real()) return true;
    else if (a.real()==b.real()) {
      return a.imag()<b.imag();
    }
    return false;
  }

  bool operator<(const CVec4D &a,const CVec4D &b)
  {
    if (a[0]<b[0]) return true;
    else if (a[0]==b[0]) {
      if (a[1]<b[1]) return true;
      else if (a[1]==b[1]) {
	if (a[2]<b[2]) return true;
	else if (a[2]==b[2]) {
	  return a[3]<b[3];
	}
      }
    }
    return false;
  }

  bool operator<(const CAsT4D &a,const CAsT4D &b)
  {
    if (a[0]<b[0]) return true;
    else if (a[0]==b[0]) {
      if (a[1]<b[1]) return true;
      else if (a[1]==b[1]) {
	if (a[2]<b[2]) return true;
	else if (a[2]==b[2]) {
	  if (a[3]<b[3]) return true;
	  else if (a[3]==b[3]) {
	    if (a[4]<b[4]) return true;
	    else if (a[4]==b[4]) {
	      return a[5]<b[5];
	    }
	  }
	}
      }
    }
    return false;
  }

}

void CDBG_Amplitude::CalcJL(const size_t &n,const size_t &id)
{
  if (n==1) {
    m_jv=p_ch[id]>0?p_ep[id]:p_em[id];
    m_jv(0)=p_cl[id][0];
    m_jv(1)=p_cl[id][1];
    m_cur[n][id]->SetP(p_p[id]);
    m_cur[n][id]->SetCurrent();
    return;
  }
  m_cur[n][id]->Evaluate();
#ifdef DEBUG__BG
  Vector_Current *vc(dynamic_cast<Vector_Current*>(m_cur[n][id]));
  if (vc!=NULL) {
    if (!vc->J().empty()) {
      ++s_nonzero; 
      std::set<CVec4D> evs;
      for (size_t i(0);i<vc->J().size();++i) {
	CVec4D ev(vc->J()[i]);
	if (evs.find(ev)==evs.end()) {
	  evs.insert(ev);
	  ++s_nonzeroall;
	}
      }
    }
  }
  else {
    Tensor_Current *tc(dynamic_cast<Tensor_Current*>(m_cur[n][id]));
    if (tc!=NULL) {
      if (!tc->J().empty()) {
	++s_nonzero; 
	std::set<CAsT4D> ets;
	for (size_t i(0);i<tc->J().size();++i) {
	  CAsT4D et(tc->J()[i]);
	  if (ets.find(et)==ets.end()) {
	    ets.insert(et);
	    ++s_nonzeroall;
	  }
	}
      }
    }
  }
#endif
}

void CDBG_Amplitude::CalcJL()
{
  PROFILE_HERE;
  for (size_t n(1);n<m_n;++n) 
    for (size_t i(0);i<m_cur[n].size();++i) 
      CalcJL(n,i);
}

void CDBG_Amplitude::SetMomenta(const std::vector<Vec4D> &moms)
{
  PROFILE_HERE;
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"():\n";
#endif
  Vec4D sum;
  for (size_t i(0);i<m_n;++i) {
    p_ssp[i]=p_p[i]=moms[i];
    p_ssep[i]=p_ep[i]=EP(p_p[i]);
    p_ssem[i]=p_em[i]=p_ep[i].Conj();
    sum+=p_p[i];
#ifdef DEBUG__BG
    msg_Debugging()<<"set p["<<i<<"] = "<<p_p[i]<<"\n";
#endif
  }
  static double accu(sqrt(Accu()));
  CVec4D::SetAccu(accu);
  if (!((CVec4D)sum).IsZero()) 
    msg.Error()<<METHOD<<"(): Four momentum not conserved. sum = "
	       <<sum<<"."<<std::endl;
  CVec4D::SetAccu(Accu());
}

void CDBG_Amplitude::SetColors(const std::vector<int> &rc,
			     const std::vector<int> &ac)
{
  PROFILE_HERE;
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"():\n";
#endif
  for (size_t i(0);i<m_n;++i) {
    if (m_mode>0) {
      p_cl[i][1]=i;
      p_cl[i][0]=i+1;
    }
    else {
      if (p_fl[i].IsGluon()) {
	p_cl[i][0]=rc[i];
	p_cl[i][1]=ac[i];
      }
      else {
	if (p_fl[i].IsAnti()) p_cl[i][1]=ac[i];
	else p_cl[i][0]=rc[i];
      }
    }
#ifdef DEBUG__BG
    msg_Debugging()<<"p_cl["<<i<<"][0] = "<<p_cl[i][0]
		   <<", p_cl["<<i<<"][1] = "<<p_cl[i][1]<<"\n";
#endif
  }
  if (m_mode>0) p_cl[0][1]=m_n;
}

size_t CDBG_Amplitude::MakeId(const Int_Vector &ids) const
{
  if (ids.size()!=m_n) 
    THROW(fatal_error,"Invalid particle number");
  size_t id(0);
  for (size_t i(0);i<ids.size();++i) 
    if (ids[i]>0) id+=1<<i;
#ifdef DEBUG__CDBCF
  msg_Debugging()<<METHOD<<ids<<" -> "<<id<<"\n";
#endif
  return id;
}

Int_Vector CDBG_Amplitude::MakeId(const size_t &id) const
{
  size_t ic(id);
  Int_Vector ids(m_n,-1);
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

Complex CDBG_Amplitude::Evaluate(const std::vector<int> &chirs)
{
  PROFILE_HERE;
  size_t id(MakeId(chirs));
  if (id>m_maxid || m_hamps[id]==0) {
    msg.Error()<<METHOD<<"(): Configuration "
	       <<chirs<<" does not exist."<<id<<std::endl;
    return 0.0;
  }
  for (size_t j(0);j<m_n;++j) p_ch[j]=chirs[j];
  CalcJL();
  Complex res(*m_cur[1].back()**m_cur.back()[0]);
  msg_Debugging()<<"A("<<chirs<<") = "<<res<<" "
		 <<std::abs(res)<<"\n";
  return res;
}

bool CDBG_Amplitude::EvaluateAll()
{
  PROFILE_HERE;
  for (size_t i(0);i<m_ress.size();++i) {
    if (m_hmap[i]>0) {
      m_ress[i]=std::conj(m_ress[m_hmap[i]]);
      continue;
    }
    for (size_t j(0);j<m_n;++j) p_ch[j]=m_chirs[i][j];
    CalcJL();
    m_ress[i]=*m_cur[1].back()**m_cur.back()[0];
#ifdef DEBUG__BG
    msg_Debugging()<<"A("<<m_chirs[i]<<") = "<<m_ress[i]<<" "
		   <<std::abs(m_ress[i])<<"\n";
#endif
  }
  return true;
}

bool CDBG_Amplitude::EvaluateAll(const Int_Vector &perm)
{
  PROFILE_HERE;
  for (size_t i(0);i<m_ress.size();++i) {
    if (m_hmap[i]>=0) {
      m_ress[i]=std::conj(m_ress[m_hmap[i]]);
      continue;
    }
    for (size_t j(0);j<m_n;++j) {
      p_p[j]=p_ssp[perm[j]];
      p_ep[j]=p_ssep[perm[j]];
      p_em[j]=p_ssem[perm[j]];
      p_ch[j]=m_chirs[i][perm[j]];
    }
    CalcJL();
    m_ress[i]=*m_cur[1].back()**m_cur.back()[0];
#ifdef DEBUG__BG
    Int_Vector rc(m_chirs[i].size());
    for (size_t j(0);j<m_chirs[i].size();++j) rc[j]=p_ch[j];
    msg_Debugging()<<"A["<<i<<"]("<<rc<<") = "<<m_ress[i]<<" "
		   <<std::abs(m_ress[i])<<"\n";
#endif
  }
  return true;
}

bool CDBG_Amplitude::CheckChirs(const std::vector<int> &chirs)
{
  size_t p(0), m(0);
  std::vector<int> q(m_nf,0);
  for (size_t i(0);i<chirs.size();++i) {
    if (p_fl[i].IsQuark()) q[p_fl[i].Kfcode()]+=chirs[i];
    if (chirs[i]>0) ++p;
    else if (chirs[i]<0) ++m;
    else THROW(fatal_error,"Invalid helicities");
  }
  for (size_t i(0);i<q.size();++i) 
    if (q[i]!=0) return false;
  return p>1 && m>1;
}

bool CDBG_Amplitude::MapChirs(const std::vector<int> &chirs)
{
  std::vector<int> rchirs(chirs.size());
  for (size_t i(0);i<chirs.size();++i) rchirs[i]=-chirs[i];
  m_hmap.push_back(-1);
  for (size_t i(0);i<m_chirs.size();++i) {
    bool hit(true);
    for (size_t j(0);j<rchirs.size();++j)
      if (m_chirs[i][j]!=rchirs[j]) hit=false;
    if (!hit) continue;
    m_hmap.back()=i;
    msg_Debugging()<<"mapped helicities "<<(m_hmap.size()-1)
		   <<" "<<chirs<<" -> "<<m_chirs[i]<<"\n";
    return true;
  }
  return false;
}

bool CDBG_Amplitude::ConstructChirs(std::vector<int> chirs,const size_t &i)
{
  if (i==chirs.size()) {
    if (CheckChirs(chirs)) {
      msg_Debugging()<<"added helicities "<<chirs<<"\n";
      m_chirs.push_back(chirs);
      m_maxid=Max(m_maxid,MakeId(chirs));
      MapChirs(chirs);
      m_ress.push_back(0.0);
    }
    return true;
  }
  if (chirs[i]!=0) {
    ConstructChirs(chirs,i+1);
  }
  else {
    for (int ch(1);ch>=-1;ch-=2) {
      chirs[i]=ch;
      ConstructChirs(chirs,i+1);
    }
  }
  return true;
}

bool CDBG_Amplitude::Construct(std::vector<ATOOLS::Flavour> flavs,
			       std::vector<int> chirs)
{
  CleanUp();
  if (flavs.size()!=chirs.size()) 
    THROW(fatal_error,"Inconsistent indices");
  if (!Construct(flavs)) return false;
  if (!ConstructChirs(chirs,0)) return false;
  m_hamps.resize(m_maxid+1,0);
  for (size_t i(0);i<m_chirs.size();++i)
    m_hamps[MakeId(m_chirs[i])]=1;
  return true;
}

bool CDBG_Amplitude::GaugeTest(const std::vector<Vec4D> &moms)
{
  msg_Info()<<METHOD<<"(): Performing gauge test ..."<<std::flush;
  msg_Indent();
  Vec4D k(1.0,1.0,0.0,0.0);
  m_kp=Spinor(1,k);
  m_km=Spinor(-1,k);
  SetMomenta(moms);
  if (!EvaluateAll()) return false;
  Complex_Vector ress(m_ress);
  k=Vec4D(1.0,0.0,1.0,0.0);
  m_kp=Spinor(1,k);
  m_km=Spinor(-1,k);
  SetMomenta(moms);
  if (!EvaluateAll()) return false;
  msg_Debugging()<<METHOD<<"():\n";
  for (size_t i(0);i<m_ress.size();++i) {
    msg_Debugging()<<"A("<<m_chirs[i]
		   <<") = "<<m_ress[i]<<" vs. "<<ress[i]<<" -> dev. "
		   <<m_ress[i].real()/ress[i].real()-1.0<<" "
		   <<m_ress[i].imag()/ress[i].imag()-1.0<<"\n";
    double accu(sqrt(Accu()));
    if (!IsEqual(m_ress[i].real(),ress[i].real(),accu) ||
	!IsEqual(m_ress[i].imag(),ress[i].imag(),accu)) {
      msg.Error().precision(12);
      msg.Error()<<METHOD<<"(): Gauge test failed. "
		 <<m_ress[i]<<" vs. "<<ress[i]<<"."<<std::endl;
      msg.Error().precision(6);
      msg_Debugging()<<"}\n";
      return false;
    }
  }
  msg_Debugging()<<"}\n";
  msg_Info()<<"satisfied."<<std::endl;
  return true;
}

bool CDBG_Amplitude::CyclicityTest()
{
  if (m_mode==0) return true;
  msg_Debugging()<<"Cyclicity test:\n";
  std::vector<Vec4D> moms(m_n);
  std::vector<int> rc(m_n), ac(m_n);
  for (size_t k(0);k<m_n;++k) moms[k]=p_p[k];
  for (size_t i(0);i<m_ress.size();++i) {
    for (size_t k(0);k<m_n;++k) {
      msg_Debugging()<<"P(";
      for (size_t j(0);j<m_n;++j) {
	size_t ni(j+k<m_n?j+k:j+k-m_n);
 	p_ch[j]=m_chirs[i][ni];
 	p_p[j]=moms[ni];
 	p_ep[j]=EP(p_p[j]);
 	p_em[j]=p_ep[j].Conj();
	msg_Debugging()<<" "<<ni;
      }
      CalcJL();
      Complex res(*m_cur[1].back()**m_cur.back()[0]);
      msg_Debugging()<<" ): A("<<m_chirs[i]<<") = "<<res<<" vs. "
		     <<m_ress[i]<<" -> dev. "
		     <<m_ress[i].real()/res.real()-1.0<<" "
		     <<m_ress[i].imag()/res.imag()-1.0<<"\n";
      double accu(sqrt(Accu()));
      if (!IsEqual(res.real(),m_ress[i].real(),accu) ||
	  !IsEqual(res.imag(),m_ress[i].imag(),accu)) {
	msg.Error()<<METHOD<<"(): Cyclicity violated. "
		   <<res<<" vs. "<<m_ress[i]<<std::endl;
	return false;
      }
    }
  }
  SetMomenta(moms);
  return true;
}

bool CDBG_Amplitude::ReflectionTest()
{
  if (m_mode==0) return true;
  msg_Debugging()<<"Reflection test:\n";
  std::vector<Vec4D> moms(m_n);
  std::vector<int> rc(m_n), ac(m_n);
  for (size_t k(0);k<m_n;++k) moms[k]=p_p[k];
  for (size_t i(0);i<m_ress.size();++i) {
    msg_Debugging()<<"R(";
    for (size_t j(0);j<m_n;++j) {
      p_ch[j]=m_chirs[i][m_n-j-1];
      p_p[j]=moms[m_n-j-1];
      p_ep[j]=EP(p_p[j]);
      p_em[j]=p_ep[j].Conj();
      msg_Debugging()<<" "<<m_n-j-1;
    }
    CalcJL();
    Complex res(*m_cur[1].back()**m_cur.back()[0]*(m_n%2==0?1.0:-1.0));
    msg_Debugging()<<" ): (-1)^"<<m_n<<"*A("<<m_chirs[i]<<") = "
		   <<res<<" vs. "<<m_ress[i]<<" -> dev. "
		   <<m_ress[i].real()/res.real()-1.0<<" "
		   <<m_ress[i].imag()/res.imag()-1.0<<"\n";
    double accu(sqrt(Accu()));
    if (!IsEqual(res.real(),m_ress[i].real(),accu) ||
	!IsEqual(res.imag(),m_ress[i].imag(),accu)) {
      msg.Error()<<METHOD<<"(): Reflection violated. "
		 <<res<<" vs. "<<m_ress[i]<<std::endl;
      return false;
    }
  }
  SetMomenta(moms);
  return true;
}

bool CDBG_Amplitude::DWITest()
{
  if (m_mode==0) return true;
  msg_Debugging()<<"Dual Ward identity test:\n";
  std::vector<Vec4D> moms(m_n);
  std::vector<int> rc(m_n), ac(m_n), id(m_n);
  for (size_t k(0);k<m_n;++k) moms[k]=p_p[k];
  for (size_t i(0);i<m_ress.size();++i) {
    Complex sum(0.0,0.0);
    for (size_t k(0);k<m_n-1;++k) {
      id[k]=0;
      p_ch[k]=m_chirs[i][0];
      p_p[k]=moms[0];
      p_ep[k]=EP(p_p[k]);
      p_em[k]=p_ep[k].Conj();
      for (size_t j(0);j<m_n;++j) {
	if (j==k) continue;
	id[j]=j<k?j+1:j;
	p_ch[j]=m_chirs[i][j<k?j+1:j];
	p_p[j]=moms[j<k?j+1:j];
	p_ep[j]=EP(p_p[j]);
	p_em[j]=p_ep[j].Conj();
      }
      CalcJL();
      sum+=m_ress[i]=*m_cur[1].back()**m_cur.back()[0];
      msg_Debugging()<<"O(";
      for (size_t j(0);j<m_n;++j) msg_Debugging()<<" "<<id[j];
      msg_Debugging()<<" ): A("<<m_chirs[i]<<") = "<<m_ress[i]<<"\n";
    }
    msg_Debugging()<<"sum = "<<sum<<"\n";
    if (!IsZero(sum)) {
      msg.Error()<<METHOD<<"(): Dual Ward identity violated. "
		 <<sum<<std::endl;
      return false;
    }
  }
  SetMomenta(moms);
  return true;
}
