#include "ATOOLS/Phys/Weight_Info.H"

#include "ATOOLS/Phys/Blob.H"

using namespace ATOOLS;

template Weight_Info &Blob_Data_Base::Get<Weight_Info>();
template PDF_Info &Blob_Data_Base::Get<PDF_Info>();
template ME_Weight_Info &Blob_Data_Base::Get<ME_Weight_Info>();

namespace ATOOLS {
  template <> Blob_Data<Weight_Info>::~Blob_Data() {}
  template class Blob_Data<Weight_Info>;

  template <> Blob_Data<PDF_Info>::~Blob_Data() {}
  template class Blob_Data<PDF_Info>;

  template <> Blob_Data<ME_Weight_Info*>::~Blob_Data() {}
  template class Blob_Data<ME_Weight_Info*>;
}

ME_Weight_Info &ME_Weight_Info::operator*=(const double &scal)
{
  m_B*=scal;
  m_VI*=scal;
  m_KP*=scal;
  m_RS*=scal;
  if (m_type&mewgttype::muR)
    for (size_t i(0);i<m_wren.size();++i) m_wren[i]*=scal;
  if (m_type&mewgttype::muF)
    for (size_t i(0);i<m_wfac.size();++i) m_wfac[i]*=scal;
  for (size_t i(0);i<m_dadsinfos.size();++i) m_dadsinfos[i].m_wgt*=scal;
  return *this;
}

std::ostream & ATOOLS::operator<<(std::ostream & s,
                                  const ATOOLS::ME_Weight_Info & mwi)
{
  s<<"type="<<mwi.m_type<<", B="<<mwi.m_B<<", VI="<<mwi.m_VI<<", KP="<<mwi.m_KP
                        <<", RS="<<mwi.m_RS<<std::endl;
  s<<"muR2="<<mwi.m_mur2<<", muF2="<<mwi.m_muf2
   <<", oqcd="<<mwi.m_oqcd<<", oew="<<mwi.m_oew
   <<", fl1="<<mwi.m_fl1<<", fl2="<<mwi.m_fl2
   <<", x1="<<mwi.m_x1<<", x2="<<mwi.m_x2
   <<", x1p="<<mwi.m_y1<<", x2p="<<mwi.m_y2<<std::endl;
  s<<"wren="<<mwi.m_wren<<std::endl;
  s<<"wfac="<<mwi.m_wfac<<std::endl;
  for (size_t i(0);i<mwi.m_dadsinfos.size();++i)
    s<<mwi.m_dadsinfos[i]<<std::endl;
  return s;
}
