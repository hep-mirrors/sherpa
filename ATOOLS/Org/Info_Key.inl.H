#ifndef ATOOLS_Org_Info_Key_inl_H
#define ATOOLS_Org_Info_Key_inl_H

#include "ATOOLS/Org/Info_Key.H"
#include "ATOOLS/Org/Integration_Info.H"

namespace ATOOLS {

  inline si::code Info_Key::SetStatus(const si::code status) const
  { return p_info->SetStatus(m_valuekey,status); }

  inline const std::string &Info_Key::Name() const
  { return m_name; }

  inline const std::string &Info_Key::Info() const
  { return m_info; }

  inline Double_Container &Info_Key::Doubles() const
  { return p_info->Doubles(m_valuekey); }

  inline Vector_Container &Info_Key::Vectors() const
  { return p_info->Vectors(m_valuekey); }

  inline si::code Info_Key::Status() const
  { return p_info->Status(m_valuekey); }

  inline double Info_Key::Weight() const
  { return p_info->Weight(m_valuekey,m_weightkey); }

  inline double Info_Key::Double(const size_t i) const
  { return p_info->Double(m_valuekey,i); }

  inline ATOOLS::Vec4D Info_Key::Vector(const size_t i) const
  { return p_info->Vector(m_valuekey,i); }

  inline bool Info_Key::Assigned() const
  { return p_info!=NULL; }

  inline double &Info_Key::operator[](const size_t i)
  { return p_info->Double(m_valuekey,i); }

  inline ATOOLS::Vec4D &Info_Key::operator()(const size_t i)
  { return p_info->Vector(m_valuekey,i); }

  inline void Info_Key::operator<<(const double weight)
  { p_info->SetWeight(m_valuekey,m_weightkey,weight); }

  inline Info_Key &Info_Key::operator=(const Info_Key &key)
  {
    for (size_t i=0;i<Doubles().size();++i) (*this)[i]=key.Double(i); 
    for (size_t j=0;j<Vectors().size();++j) (*this)(j)=key.Vector(j); 
    return *this;
  }

  inline bool Info_Key::operator==(const Info_Key &key)
  {
    if (key.m_valuekey==m_valuekey && key.m_weightkey==m_weightkey) return true; 
    return false;
  }

}// end of namespace ATOOLS

#endif
