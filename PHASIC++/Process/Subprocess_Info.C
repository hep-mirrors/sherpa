#include "PHASIC++/Process/Subprocess_Info.H"

#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

using namespace PHASIC;
using namespace ATOOLS;

std::ostream &PHASIC::operator<<(std::ostream &ostr,const Subprocess_Info &info)
{
  info.Print(ostr);
  return ostr;
}

Subprocess_Info::Subprocess_Info
(const ATOOLS::Flavour &fl,const std::string &id,
 const std::string &pol,const std::string &mpl):
  m_fl(fl), m_id(id), m_pol(pol), m_mpl(mpl),
  m_nmax(0), m_nmin(100), m_tag(0), m_osf(0),
  m_nlotype(nlo_type::lo), m_nlocpl(2,0),
  m_asscontribs(asscontrib::none) {}

Subprocess_Info::~Subprocess_Info()
{
  DeleteDecayInfos();
}

void Subprocess_Info::DeleteDecayInfos()
{
  for (size_t i(0);i<m_decins.size();++i)
    if (!m_decins[i]) delete m_decins[i];
}

std::string Subprocess_Info::MultiplicityTag() const
{
  size_t pn(0);
  std::string id;
  for (size_t i(0);i<m_ps.size();++i)
    if (m_ps[i].NExternal()>1) {
      if (pn>0) id+=ToString(pn);
      id+="["+m_ps[i].MultiplicityTag()+"]";
      pn=0;
    }
    else ++pn;
  return id+ToString(pn);
}

size_t Subprocess_Info::NExternal() const
{
  if (m_ps.empty()) return 1;
  size_t n(0);
  for (size_t i(0);i<m_ps.size();++i) n+=m_ps[i].NExternal();
  return n;
}

size_t Subprocess_Info::NTotalExternal() const
{
  if (m_ps.empty()) return m_fl.Size();
  size_t n(0);
  for (size_t i(0);i<m_ps.size();++i) n+=m_ps[i].NTotalExternal();
  return n;
}

void Subprocess_Info::SetExternal(const std::vector<ATOOLS::Flavour> &fl,size_t &n)
{
  if (m_ps.empty()) m_fl=fl[n++];
  else for (size_t i(0);i<m_ps.size();++i) m_ps[i].SetExternal(fl,n);
}

void Subprocess_Info::SetExternal(const std::vector<ATOOLS::Flavour> &fl)
{
  size_t n(0);
  SetExternal(fl,n);
}

bool Subprocess_Info::SetExternal(const ATOOLS::Flavour &fl,
			    const size_t &i,size_t &n)
{
  if (m_ps.empty()) {
    if (n==i) m_fl=fl;
    return i==n++;
  }
  for (size_t j(0);j<m_ps.size();++j) 
    if (m_ps[j].SetExternal(fl,i,n)) return true;
  return false;
}

void Subprocess_Info::SetExternal(const ATOOLS::Flavour &fl,const size_t &i)
{
  size_t n(0);
  SetExternal(fl,i,n);
}

void Subprocess_Info::GetExternal(std::vector<ATOOLS::Flavour> &fl) const
{
  if (m_ps.empty()) fl.push_back(m_fl);
  else for (size_t i(0);i<m_ps.size();++i) m_ps[i].GetExternal(fl);
}

std::vector<ATOOLS::Flavour> Subprocess_Info::GetExternal() const
{
  std::vector<ATOOLS::Flavour> fl;
  GetExternal(fl);
  return fl;
}

bool Subprocess_Info::GetExternal(ATOOLS::Flavour &fl,
			    const size_t &i,size_t &n) const
{
  if (m_ps.empty()) {
    if (n==i) fl=m_fl;
    return i==n++;
  }
  for (size_t j(0);j<m_ps.size();++j) 
    if (m_ps[j].GetExternal(fl,i,n)) return true;
  return false;
}

ATOOLS::Flavour Subprocess_Info::GetExternal(const size_t &i) const
{
  size_t n(0);
  Flavour fl(kf_none);
  GetExternal(fl,i,n);
  return fl;
}

void Subprocess_Info::Add(const Subprocess_Info &info)
{
  m_ps.insert(m_ps.end(),info.m_ps.begin(),info.m_ps.end());
}

bool Subprocess_Info::AddDecay
(const Subprocess_Info &ii,const Subprocess_Info &fi, int osf)
{
  if (m_ps.empty()) {
    if (m_fl==ii.m_ps.front().m_fl &&
	m_id==ii.m_ps.front().m_id) {
      m_ps=fi.m_ps;
      m_nlotype=fi.m_nlotype;
      m_osf=osf;
    }
    return m_ps.size()>0;
  }
  for (size_t i(0);i<m_ps.size();++i)
    if (m_ps[i].AddDecay(ii,fi,osf)) return true;
  return false;
}

size_t Subprocess_Info::GetDecayInfos
(DecayInfo_Vector &ids,size_t &n,bool init)
{
  if (init) ids.clear();
  if (m_ps.empty()) return 1<<n++;
  size_t cont(0);
  DecayInfo_Vector subsequentdecays;
  for (size_t j(0);j<m_ps.size();++j) {
    size_t decinfs(ids.size());
    cont+=m_ps[j].GetDecayInfos(ids,n,false);
    if (decinfs+1==ids.size()) subsequentdecays.push_back(ids.back());
  }
  ids.push_back(new Decay_Info(cont,m_fl,m_nmax,m_osf));
  ids.back()->SetSubsequentDecayInfos(subsequentdecays);
  return cont;
}

void Subprocess_Info::BuildDecayInfos(size_t nin)
{
  DeleteDecayInfos();
  size_t n(nin);
  GetDecayInfos(m_decins,n,true);
  delete m_decins.back();
  m_decins.pop_back();
}

double Subprocess_Info::Factorial(const double &n) const
{
  if (n<=0.0) return 1.0;
  return n*Factorial(n-1.0);
}

double Subprocess_Info::ISSymmetryFactor() const
{
  double sf(1.0);
  msg_Debugging()<<METHOD<<"(): {\n";
  for (size_t i(0);i<m_ps.size();++i) {
    double pols(2.0*m_ps[i].m_fl.Spin()+1.0);
    if (m_ps[i].m_fl.IntSpin()==2 &&
	m_ps[i].m_fl.Mass()==0.0) pols=2.0;
    sf*=pols;
    msg_Debugging()<<"     "<<std::setw(15)
		   <<m_ps[i].m_fl<<" -> "<<pols;
    if (m_ps[i].m_fl.Strong()) {
      sf*=abs(m_ps[i].m_fl.StrongCharge());
      msg_Debugging()<<"*"<<abs(m_ps[i].m_fl.StrongCharge());
    }
    msg_Debugging()<<"\n";
  }
  msg_Debugging()<<"} -> "<<sf<<"\n";
  return sf;
}

double Subprocess_Info::FSSymmetryFactor() const
{
  double sf(1.0);
  msg_Debugging()<<METHOD<<"(): {\n";
  std::map<Flavour,size_t> fc;
  for (size_t i(0);i<m_ps.size();++i) {
    if (m_ps[i].m_ps.size()) {
      msg_Indent();
      sf*=m_ps[i].FSSymmetryFactor();
    }
    std::map<Flavour,size_t>::iterator fit(fc.find(m_ps[i].m_fl));
    if (fit==fc.end()) {
      fc[m_ps[i].m_fl]=0;
      fit=fc.find(m_ps[i].m_fl);
    }
    ++fit->second;
  }
  for (std::map<Flavour,size_t>::const_iterator fit(fc.begin());
       fit!=fc.end();++fit) {
    msg_Debugging()<<"  "<<std::setw(2)<<fit->second<<" "
		   <<std::setw(15)<<fit->first<<" -> "
		   <<std::setw(12)<<Factorial(fit->second)<<"\n";
    sf*=Factorial(fit->second);
  }
  msg_Debugging()<<"} -> "<<sf<<"\n";
  return sf;
}

size_t Subprocess_Info::NMinExternal() const
{
  if (m_ps.empty()) return 1;
  size_t n(m_nmin-m_ps.size());
  for (size_t i(0);i<m_ps.size();++i) n+=m_ps[i].NMinExternal();
  return n;
}

size_t Subprocess_Info::NMaxExternal() const
{
  if (m_ps.empty()) return 1;
  size_t n(m_nmax-m_ps.size());
  for (size_t i(0);i<m_ps.size();++i) n+=m_ps[i].NMaxExternal();
  return n;
}

bool Subprocess_Info::IsGroup() const
{
  if (m_ps.empty()) return m_fl.IsGroup();
  size_t naregroup(0);
  for (size_t i(0);i<m_ps.size();++i) naregroup+=m_ps[i].IsGroup();
  return naregroup;
}

void Subprocess_Info::SetNMax(const Subprocess_Info &ref)
{
  m_nmin=Min(m_ps.size(),ref.m_nmin);
  m_nmax=Max(m_ps.size(),ref.m_nmax);
  size_t lim(Min(m_ps.size(),ref.m_ps.size()));
  for (size_t j(0);j<lim;++j) m_ps[j].SetNMax(ref.m_ps[j]);
}

void Subprocess_Info::GetNMax(const Subprocess_Info &ref)
{
  m_nmin=Min(m_nmin,ref.m_ps.size());
  m_nmax=Max(m_nmax,ref.m_ps.size());
  size_t lim(Min(m_ps.size(),ref.m_ps.size()));
  for (size_t j(lim);j<ref.m_ps.size();++j) 
    m_ps.push_back(Subprocess_Info(ref.m_ps[j].m_fl,ref.m_ps[j].m_id));
  for (size_t j(0);j<ref.m_ps.size();++j) m_ps[j].GetNMax(ref.m_ps[j]);
}

void Subprocess_Info::SetTags(int& start)
{
  if (m_ps.size()==0) {
    m_tag=start;
    ++start;
  }
  else {
    for (size_t i=0; i<m_ps.size(); ++i) {
      m_ps[i].SetTags(start);
    }
  }
}

void Subprocess_Info::SetTags(const std::vector<int>& tags)
{
  int n=0;
  SetTags(tags,n);
}

void Subprocess_Info::SetTags(const std::vector<int>& tags,int &n)
{
  if (m_ps.size()==0) m_tag=tags[n++];
  else for (size_t i=0;i<m_ps.size();++i) m_ps[i].SetTags(tags,n);
}

void Subprocess_Info::GetTags(std::vector<int>& tags) const
{
  if (m_ps.size()==0) {
    tags.push_back(m_tag);
  }
  else {
    for (size_t i=0; i<m_ps.size(); ++i) {
      m_ps[i].GetTags(tags);
    }
  }
}

int Subprocess_Info::Combine
(const size_t &i,const size_t &j,const Flavour &flij,int &cnt)
{
  if (m_ps.size()==0) {
    ++cnt;
    if (cnt-1==i) m_fl=flij;
    if (cnt-1==j) return -1;
    return 0;
  }
  else {
    for (std::vector<Subprocess_Info>::iterator 
	   psit(m_ps.begin());psit!=m_ps.end();++psit) {
      int stat(psit->Combine(i,j,flij,cnt));
      if (stat<0) psit=m_ps.erase(psit)-1;
    }
    return 1;
  }
}

void Subprocess_Info::ExtractMPL(std::vector<Flavour_Vector> &fl) const
{
  if (m_ps.size()) {
    for (size_t i=0;i<m_ps.size();++i) m_ps[i].ExtractMPL(fl);
    return;
  }
  if (m_mpl=="") {
    fl.push_back(Flavour_Vector(1,m_fl));
    return;
  }
  fl.push_back(Flavour_Vector());
  std::string mpl(m_mpl);
  for (size_t pos(mpl.find(','));pos!=std::string::npos;
       mpl=mpl.substr(pos+1),pos=mpl.find(',')) {
    std::string cur(mpl.substr(0,pos));
    fl.back().push_back(Flavour(ToType<long int>(cur)));
  }
  fl.back().push_back(Flavour(ToType<long int>(mpl)));
}

bool Subprocess_Info::operator<(const Subprocess_Info &pi) const
{
  if (m_ps.size()<pi.m_ps.size()) return true;
  if (m_ps.size()>pi.m_ps.size()) return false;
  if (m_ps.empty()) return m_fl<pi.m_fl;
  for (size_t i(0);i<m_ps.size();++i) {
    if (m_ps[i]<pi.m_ps[i]) return true;
    if (!(m_ps[i]==pi.m_ps[i])) return false;
  }
  return false;
}

bool Subprocess_Info::operator==(const Subprocess_Info &pi) const
{
  if (m_ps.size()!=pi.m_ps.size()) return false;
  if (m_ps.empty()) return m_fl==pi.m_fl;
  for (size_t i(0);i<m_ps.size();++i)
    if (!(m_ps[i]==pi.m_ps[i])) return false;
  return true;
}

void Subprocess_Info::Print(std::ostream &ostr,const size_t &ni) const
{
  ostr<<std::string(ni,' ')<<m_fl<<" "<<m_mpl;
  if (m_id!="") ostr<<"["<<m_id<<"]";
  if (m_osf) ostr<<" OS";
  if (m_ps.size()>0) {
    ostr<<" ("<<m_ps.size()<<")";
    ostr<<", NLO{"<<m_nlotype<<"} "<<m_nlocpl;
    ostr<<", ASS{"<<m_asscontribs<<"}";
    if (m_nmax>0) ostr<<" {"<<m_nmin<<","<<m_nmax<<"}";
    ostr <<": {\n";
    for (size_t i(0);i<m_ps.size();++i) m_ps[i].Print(ostr,ni+2);
    ostr<<std::string(ni,' ')<<"}";
  }
  ostr<<"\n";
}

void Subprocess_Info::PrintFlavours(std::ostream &ostr) const
{
  if (m_ps.size()>0) {
    for (size_t i(0);i<m_ps.size();++i) m_ps[i].PrintFlavours(ostr);
  } else {
    ostr<<(long int)m_fl<<" "<<m_mpl;
  }
}
