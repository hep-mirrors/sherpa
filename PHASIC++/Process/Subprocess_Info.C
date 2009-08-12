#include "PHASIC++/Process/Subprocess_Info.H"

#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

using namespace PHASIC;
using namespace ATOOLS;

std::string PHASIC::PSId(const size_t &id)
{
  size_t ic(id);
  std::string idt;
  for (size_t i(0);ic>0;++i) {
    size_t c(1<<i);
    if (ic&c) {
      char nic[3];
      if (sprintf(nic,"%i",(int)i)<=0)
	THROW(fatal_error,"Conversion error");
      idt+=nic;
      ic-=c;
    }
  }
  return idt;
}

std::ostream &PHASIC::operator<<(std::ostream &ostr,const Subprocess_Info &info)
{
  info.Print(ostr);
  return ostr;
}

Subprocess_Info::Subprocess_Info
(const ATOOLS::Flavour &fl,const std::string &id,const std::string &pol):
  m_fl(fl), m_id(id), m_pol(pol), m_nmax(0), m_tag(0), m_osf(0),
  m_nloqcdtype(nlo_type::lo), m_nloewtype(nlo_type::lo) {}

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
  else for (size_t i(0);i<m_ps.size();++i) m_ps[i].SetExternal(fl);
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
      m_nloqcdtype=fi.m_nloqcdtype;
      m_nloewtype=fi.m_nloewtype;
      m_osf=osf;
    }
    return m_ps.size()>0;
  }
  for (size_t i(0);i<m_ps.size();++i)
    if (m_ps[i].AddDecay(ii,fi,osf)) return true;
  return false;
}

size_t Subprocess_Info::GetDecayInfos
(Decay_Info_Vector &ids,size_t &n) const
{
  if (m_ps.empty()) return 1<<n++;
  size_t cont(0);
  for (size_t j(0);j<m_ps.size();++j) 
    cont+=m_ps[j].GetDecayInfos(ids,n);
  ids.push_back(Decay_Info(cont,m_fl,m_nmax,m_osf));
  return cont;
}

Decay_Info_Vector Subprocess_Info::GetDecayInfos() const
{
  size_t n(0);
  std::vector<Decay_Info> ids;
  GetDecayInfos(ids,n);
  ids.pop_back();
  return ids;
}

size_t Subprocess_Info::NMaxExternal() const
{
  if (m_ps.empty()) return 1;
  size_t n(m_nmax-m_ps.size());
  for (size_t i(0);i<m_ps.size();++i) n+=m_ps[i].NMaxExternal();
  return n;
}

void Subprocess_Info::SetNMax(const Subprocess_Info &ref)
{
  m_nmax=Max(m_ps.size(),ref.m_nmax);
  size_t lim(Min(m_ps.size(),ref.m_ps.size()));
  for (size_t j(0);j<lim;++j) m_ps[j].SetNMax(ref.m_ps[j]);
}

void Subprocess_Info::GetNMax(const Subprocess_Info &ref)
{
  m_nmax=ref.m_ps.size();
  size_t lim(Min(m_ps.size(),ref.m_ps.size()));
  for (size_t j(lim);j<ref.m_ps.size();++j) 
    m_ps.push_back(Subprocess_Info(ref.m_ps[j].m_fl,ref.m_ps[j].m_id));
  for (size_t j(0);j<ref.m_ps.size();++j) m_ps[j].GetNMax(ref.m_ps[j]);
}

nlo_type::code Subprocess_Info::NLOType() const
{
  if (m_nloewtype==nlo_type::lo) {
    return m_nloqcdtype;
  }
  else if (m_nloqcdtype==nlo_type::lo) {
    return m_nloewtype;
  }
  else {
    THROW(fatal_error, "Can't handle NLO EW and NLO QCD in one amplitude.");
    return nlo_type::lo;
  }
}

void Subprocess_Info::SetNLOType(nlo_type::code nlotype)
{
  if (m_nloewtype==nlo_type::lo) {
    m_nloqcdtype=nlotype;
  }
  else if (m_nloqcdtype==nlo_type::lo) {
    m_nloewtype=nlotype;
  }
  else {
    THROW(fatal_error, "Tried to set NLOType for non-NLO amplitude.");
  }
}

void Subprocess_Info::Print(std::ostream &ostr,const size_t &ni) const
{
  ostr<<std::string(ni,' ')<<m_fl;
  if (m_id!="") ostr<<"["<<m_id<<"]";
  if (m_osf) ostr<<" OS";
  if (m_ps.size()>0) {
    ostr<<" ("<<m_ps.size()<<")";
    ostr<<", NLO{"<<m_nloqcdtype<<","<<m_nloewtype<<"}";
    if (m_nmax>0) ostr<<"{"<<m_nmax<<"}";
    ostr <<": {\n";
    for (size_t i(0);i<m_ps.size();++i) m_ps[i].Print(ostr,ni+2);
    ostr<<std::string(ni,' ')<<"}";
  }
  ostr<<"\n";
}
