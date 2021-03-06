#include "PHASIC++/Process/Process_Info.H"
#include "ATOOLS/Org/Exception.H"

#include "ATOOLS/Org/Message.H"

using namespace PHASIC;

std::ostream &PHASIC::operator<<(std::ostream &ostr,const Process_Info &info)
{
  ostr<<"("<<&info<<"){\n";
  {
    ostr<<"  cls = "<<info.m_cls<<", hls = "<<info.m_hls<<"\n";
    ostr<<"  mincpl = "<<info.m_mincpl<<"/"<<info.m_minacpl
	<<", maxcpl = "<<info.m_maxcpl<<"/"<<info.m_maxacpl<<"\n";
    ostr<<"  ckkw = "<<info.m_ckkw
	<<", nlo = "<<info.m_nlomode<<", mhv = "<<info.m_amegicmhv<<"\n";
    ostr<<"  scale = '"<<info.m_scale<<"', kfactor = '"<<info.m_kfactor<<"'\n";
    ostr<<"  megenerator = '"<<info.m_megenerator
	<<"',  loopgenerator = '"<<info.m_loopgenerator<<"'\n";
    if (info.m_special!="") ostr<<"  special = "<<info.m_special<<"\n";
    ostr<<"  gpath = '"<<info.m_gpath
	<<"', min t-channels = "<<info.m_ntchan
	<<"', max t-channels = "<<info.m_mtchan<<"\n";
    if (info.m_nodecs.size()) ostr<<"  nodecs = "<<info.m_nodecs<<"\n";
    info.m_ii.Print(ostr,2);
    info.m_fi.Print(ostr,2);
  }
  ostr<<"}";
  return ostr;
}

ATOOLS::Flavour_Vector Process_Info::ExtractFlavours() const
{
  ATOOLS::Flavour_Vector flavs=m_ii.GetExternal();
  ATOOLS::Flavour_Vector fi=m_fi.GetExternal();
  flavs.insert(flavs.end(), fi.begin(), fi.end());
  return flavs;
}

std::vector<ATOOLS::Flavour_Vector> Process_Info::ExtractMPL() const
{
  std::vector<ATOOLS::Flavour_Vector> fl;
  m_ii.ExtractMPL(fl);
  m_fi.ExtractMPL(fl);
  return fl;
}

bool Process_Info::Has(ATOOLS::nlo_type::code nlotype) const
{
  return m_fi.m_nlotype&nlotype;
}

int Process_Info::Combine(const size_t &i,const size_t &j,
			  const ATOOLS::Flavour &flij)
{
  int cnt(0);
  int res(m_ii.Combine(i,j,flij,cnt));
  if (res<0) THROW(fatal_error,"Removed initial state particle");
  res=m_fi.Combine(i,j,flij,cnt);
  return -res;
}

bool Process_Info::operator<(const Process_Info &pi) const
{
  if (m_ii<pi.m_ii) return true;
  if (!(m_ii==pi.m_ii)) return false;
  if (m_fi<pi.m_fi) return true;
  if (!(m_fi==pi.m_fi)) return false;
  return false;
}
