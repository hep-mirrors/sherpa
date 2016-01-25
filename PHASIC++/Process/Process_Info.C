#include "PHASIC++/Process/Process_Info.H"
#include "ATOOLS/Org/Exception.H"

#include "ATOOLS/Org/Message.H"

using namespace PHASIC;

std::ostream &PHASIC::operator<<(std::ostream &ostr,const Process_Info &info)
{
  ostr<<"("<<&info<<"){\n";
  {
    ostr<<"  cls = "<<info.m_cls<<", hls = "<<info.m_hls<<"\n";
    ostr<<"  maxcpl = "<<info.m_maxcpl<<", mincpl = "<<info.m_mincpl<<"\n";
    ostr<<"  ckkw = "<<info.m_ckkw
	<<", nlo = "<<info.m_nlomode<<", mhv = "<<info.m_amegicmhv<<"\n";
    ostr<<"  scale = '"<<info.m_scale<<"', kfactor = '"<<info.m_kfactor<<"'\n";
    ostr<<"  megenerator = '"<<info.m_megenerator
	<<"',  loopgenerator = '"<<info.m_loopgenerator<<"'\n  selectorfile = '"
        <<info.m_selectorfile<<"'\n";
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

bool Process_Info::Has(nlo_type::code nlotype) const
{
  if (m_fi.m_nloewtype==nlo_type::lo) {
    return (m_fi.m_nloqcdtype&nlotype) ? true : false;
  }
  else if (m_fi.m_nloqcdtype==nlo_type::lo) {
    return (m_fi.m_nloewtype&nlotype) ? true : false;
  }
  else {
    THROW(fatal_error, "Can't handle NLO EW and NLO QCD in one amplitude.");
  }
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
