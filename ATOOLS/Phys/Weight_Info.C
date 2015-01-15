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
  m_w0*=scal;
  if (m_type&mewgttype::muR)
    for (size_t i(0);i<m_wren.size();++i) m_wren[i]*=scal;
  if (m_type&mewgttype::muF)
    for (size_t i(0);i<m_wfac.size();++i) m_wfac[i]*=scal;
  for (size_t i(0);i<m_dadsinfos.size();++i) m_dadsinfos[i].m_wgt*=scal;
  return *this;
}

