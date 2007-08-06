#include "LEV_Group.H"

#include "Exception.H"

using namespace EXTRAXS;
using namespace ATOOLS;

LEV_Group::LEV_Group() {}

LEV_Group::~LEV_Group() 
{
  while (!m_levs.empty()) {
    delete m_levs.back();
    m_levs.pop_back();
  }
}

void LEV_Group::Add(LEV_Base *const lev) 
{ 
  m_levs.push_back(lev); 
  m_integrals.push_back(1.0); 
}

double LEV_Group::Value(const Vec4D &k1,const Vec4D &q1,
			const Vec4D &k2,const Vec4D &q2) const
{
  double value(0.0);
  for (size_t i(0);i<m_levs.size();++i) 
    value+=m_levs[i]->Value(k1,q1,k2,q2);
  return value;
}

double LEV_Group::MajorValue(const Vec4D &k1,const Vec4D &q1,
			     const Vec4D &k2,const Vec4D &q2) const
{
  double value(0.0);
  for (size_t i(0);i<m_levs.size();++i) 
    value+=m_levs[i]->MajorValue(k1,q1,k2,q2);
  return value;
}

double LEV_Group::MajorIntegral(const ATOOLS::Flavour &fl)
{
  double value(0.0);
  for (size_t i(0);i<m_levs.size();++i) {
    m_integrals[i]=value+=m_levs[i]->MajorIntegral(fl);
  }
  return value;
}

bool LEV_Group::SelectLEV(const double &rn)
{
  for (size_t i(0);i<m_levs.size();++i) {
    if (m_integrals[i]/m_integrals.back()>=rn) {
      p_selected=m_levs[i];
      return true;
    }
  }
  return false;
}

void LEV_Group::SetKT2Min(const double &kt2min) 
{
  for (size_t i(0);i<m_levs.size();++i) 
    m_levs[i]->SetKT2Min(kt2min);
  m_kt2min=kt2min; 
}

void LEV_Group::SetKT2Max(const double &kt2max) 
{
  for (size_t i(0);i<m_levs.size();++i) 
    m_levs[i]->SetKT2Max(kt2max);
  m_kt2max=kt2max; 
}

LEV_Base *LEV_Group::Selected() const
{
  return p_selected->Selected();
}

