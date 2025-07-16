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

  std::ostream& operator<<(std::ostream& s,
                           const PDF_Info& pi)
  {
    s<<" pdf1 = ("<<pi.m_fl1<<","<<pi.m_x1<<","<<sqrt(pi.m_muf12)
                  <<":"<<pi.m_xf1<<") , "
     <<" pdf2 = ("<<pi.m_fl2<<","<<pi.m_x2<<","<<sqrt(pi.m_muf22)
                  <<":"<<pi.m_xf2<<")";
    return s;
  }

  std::ostream& operator<<(std::ostream& s,
                           const Weight_Info& wi)
  {
    return s<<" w = "<<wi.m_weightsmap.Nominal()<<", dxs = "<<wi.m_dxs
            <<", trials = "<<wi.m_ntrial
            <<", pdfs = { "<<wi.m_pdf<<" }"<<std::endl;
  }
}

bool PDF_Info::operator==(const PDF_Info& rhs) const
{
  if (m_fl1 != rhs.m_fl1) return false;
  if (m_fl2 != rhs.m_fl2) return false;
  if (!IsEqual(m_x1, rhs.m_x1, 1e-6)) return false;
  if (!IsEqual(m_x2, rhs.m_x2, 1e-6)) return false;
  if (!IsEqual(m_muf12, rhs.m_muf12, 1e-6)) return false;
  if (!IsEqual(m_muf22, rhs.m_muf22, 1e-6)) return false;
  if (!IsEqual(m_xf1, rhs.m_xf1, 1e-6)) return false;
  if (!IsEqual(m_xf2, rhs.m_xf2, 1e-6)) return false;
  return true;
}

bool PDF_Info::operator!=(const PDF_Info& rhs) const
{
  return !(*this == rhs);
}

bool Weight_Info::operator==(const Weight_Info& rhs) const
{
  // NOTE: Only compare nominals for now, since Weights_Map does not yet offer
  // a "=="-operator.
  if (!IsEqual(m_weightsmap.Nominal(), rhs.m_weightsmap.Nominal(), 1e-6))
    return false;
  if (!IsEqual(m_dxs, rhs.m_dxs, 1e-6)) return false;
  if (!IsEqual(m_ntrial, rhs.m_ntrial, 1e-6)) return false;
  if (m_pdf != rhs.m_pdf) return false;
  return true;
}

bool Weight_Info::operator!=(const Weight_Info& rhs) const
{
  return !(*this == rhs);
}

std::ostream & ATOOLS::operator<<(std::ostream & s,
                                 const ATOOLS::mewgttype::code & type)
{
  if (type==mewgttype::none) s<<"none";
  if (type&mewgttype::B)     s<<"B";
  if (type&mewgttype::VI)    s<<"VI";
  if (type&mewgttype::KP)    s<<"KP";
  if (type&mewgttype::DADS)  s<<"DADS";
  if (type&mewgttype::H)     s<<"H";
  if (type&mewgttype::RS)    s<<"RS";
  if (type&mewgttype::METS && type^mewgttype::METS)  s<<"|";
  if (type&mewgttype::METS)  s<<"METS";
  return s;
}

std::ostream & ATOOLS::operator<<(std::ostream & s,
                                  const ATOOLS::Cluster_Sequence_Info & csi)
{
  s<<"Cluster sequence: pdfwgt="<<csi.m_pdfwgt<<", flux="<<csi.m_flux
   <<", counter term="<<csi.m_ct;
  if (!csi.m_txfl.size()) s<<", no cluster steps";
  else                    s<<", steps: ";
  for (size_t i(0);i<csi.m_txfl.size();++i)
    s<<csi.m_txfl[i]<<" ";
  return s;
}

ME_Weight_Info &ME_Weight_Info::operator*=(const double &scal)
{
  m_B*=scal;
  m_VI*=scal;
  m_KP*=scal;
  if (m_type&mewgttype::VI)
    for (size_t i(0);i<m_wren.size();++i) m_wren[i]*=scal;
  if (m_type&mewgttype::KP)
    for (size_t i(0);i<m_wfac.size();++i) m_wfac[i]*=scal;
  for (size_t i(0);i<m_wass.size();++i) m_wass[i]*=scal;
  for (size_t i(0);i<m_dadsinfos.size();++i) m_dadsinfos[i].m_wgt*=scal;
  for (size_t i(0);i<m_rdainfos.size();++i) m_rdainfos[i].m_wgt*=scal;
  return *this;
}

void ME_Weight_Info::Reset()
{
  // undo DADS, METS, H settings as they are set event-by-event outside the
  // ME generators
  if (m_type&mewgttype::DADS) m_type^=mewgttype::DADS;
  if (m_type&mewgttype::METS) m_type^=mewgttype::METS;
  m_B=m_VI=m_KP=m_K=0.;
  m_dadsinfos.clear();
  m_rdainfos.clear();
  m_clusseqinfo=Cluster_Sequence_Info(1.,0.);
  m_x1=m_x2=m_y1=m_y2=1.;
  m_oqcd=m_oew=0;
  m_fl1=m_fl2=0;
  m_mur2=m_muf2=0.;
  if (m_type&mewgttype::VI) for (size_t i(0);i<m_wren.size();++i) m_wren[i]=0.;
  if (m_type&mewgttype::KP) for (size_t i(0);i<m_wfac.size();++i) m_wfac[i]=0.;
  for (size_t i(0);i<m_wass.size();++i) m_wass[i]=0.;
}

std::ostream & ATOOLS::operator<<(std::ostream & s,
                                  const ATOOLS::ME_Weight_Info & mwi)
{
  s<<"type="<<mwi.m_type<<", B="<<mwi.m_B<<", VI="<<mwi.m_VI<<", KP="<<mwi.m_KP
   <<", K="<<mwi.m_K<<std::endl;
  s<<"muR2="<<mwi.m_mur2<<", muF2="<<mwi.m_muf2
   <<", oqcd="<<mwi.m_oqcd<<", oew="<<mwi.m_oew
   <<", fl1="<<mwi.m_fl1<<", fl2="<<mwi.m_fl2
   <<", x1="<<mwi.m_x1<<", x2="<<mwi.m_x2
   <<", x1p="<<mwi.m_y1<<", x2p="<<mwi.m_y2<<std::endl;
  if (mwi.m_type&mewgttype::VI) s<<"wren="<<mwi.m_wren<<std::endl;
  if (mwi.m_type&mewgttype::KP) s<<"wfac="<<mwi.m_wfac<<std::endl;
  s<<"ass="<<mwi.m_wass<<std::endl;
  for (size_t i(0);i<mwi.m_dadsinfos.size();++i)
    s<<mwi.m_dadsinfos[i]<<std::endl;
  for (size_t i(0);i<mwi.m_rdainfos.size();++i)
    s<<mwi.m_rdainfos[i]<<std::endl;
  s<<mwi.m_clusseqinfo<<std::endl;
  return s;
}
