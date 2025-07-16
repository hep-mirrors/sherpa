#ifndef ATOOLS_Org_Integration_Info_inl_H
#define ATOOLS_Org_Integration_Info_inl_H

#include "ATOOLS/Org/Integration_Info.H"

namespace ATOOLS {

  inline Double_Container &Integration_Info::Doubles(const size_t valuekey)
  { return m_doubles[valuekey]; }

  inline Vector_Container &Integration_Info::Vectors(const size_t valuekey)
  { return m_vectors[valuekey]; }

  inline double &Integration_Info::Double(const size_t valuekey,const size_t i)
  { return m_doubles[valuekey][i]; }

  inline ATOOLS::Vec4D &Integration_Info::Vector(const size_t valuekey,const size_t i)
  { return m_vectors[valuekey][i]; }

  inline double Integration_Info::Double(const size_t valuekey,const size_t i) const
  { return m_doubles[valuekey][i]; }

  inline ATOOLS::Vec4D Integration_Info::Vector(const size_t valuekey,const size_t i) const
  { return m_vectors[valuekey][i]; }

  inline void Integration_Info::SetWeight(const size_t valuekey,const size_t weightkey,
					  const double weight)
  { m_weights[valuekey][weightkey]=weight; }

  inline double Integration_Info::Weight(const size_t valuekey,const size_t weightkey) const
  { return m_weights[valuekey][weightkey]; }

  inline si::code Integration_Info::Status(const size_t statuskey)
  { return m_status[statuskey]; }
  
  inline si::code Integration_Info::SetStatus(const size_t statuskey,const si::code status)
  { si::code old=m_status[statuskey]; m_status[statuskey]=status; return old; }

}// end of namespace ATOOLS

#endif
