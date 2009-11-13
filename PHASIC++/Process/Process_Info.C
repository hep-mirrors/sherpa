#include "PHASIC++/Process/Process_Info.H"

#include "ATOOLS/Org/Message.H"

using namespace PHASIC;

std::ostream &PHASIC::operator<<(std::ostream &ostr,const Process_Info &info)
{
  ostr<<"("<<&info<<"){\n";
  {
    ostr<<"  cls = "<<info.m_cls<<", hls = "<<info.m_hls<<"\n";
    ostr<<"  oew = "<<info.m_oew<<", oqcd = "<<info.m_oqcd
	<<", maxoew = "<<info.m_maxoew<<", maxoqcd = "<<info.m_maxoqcd<<"\n";
    ostr<<"  psmc = "<<info.m_psmc<<", ckkw = "<<info.m_ckkw
	<<", nlo = "<<info.m_nlomode<<", mhv = "<<info.m_amegicmhv<<"\n";
    ostr<<"  scale = '"<<info.m_scale<<"', kfactor = '"<<info.m_kfactor<<"'\n";
    ostr<<"  gpath = '"<<info.m_gpath<<"'\n";
    ostr<<"  loopgenerator = '"<<info.m_loopgenerator<<"', selectorfile = '"
        <<info.m_selectorfile<<"', mpi process = "<<info.m_mpiprocess<<"\n";
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

