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

double LEV_Group::Value(const double &y,const double &kt) const
{
  double value(0.0);
  for (size_t i(0);i<m_levs.size();++i) 
    value+=m_levs[i]->Value(y,kt);
  return value;
}

double LEV_Group::MajorValue(const double &y,const double &kt) const
{
  double value(0.0);
  for (size_t i(0);i<m_levs.size();++i) 
    value+=m_levs[i]->MajorValue(y,kt);
  return value;
}

double LEV_Group::MajorIntegral()
{
  double value(0.0);
  for (size_t i(0);i<m_levs.size();++i) {
    m_integrals[i]=value+=m_levs[i]->MajorIntegral();
  }
  return value;
}

bool LEV_Group::SelectSplitting(const double &rn)
{
  for (size_t i(0);i<m_levs.size();++i) {
    if (m_integrals[i]/m_integrals.back()>=rn) {
      p_selected=m_levs[i];
      return true;
    }
  }
  return false;
}

void LEV_Group::SetYMin(const double &ymin) 
{
  for (size_t i(0);i<m_levs.size();++i) 
    m_levs[i]->SetYMin(ymin);
  m_ymin=ymin; 
}

void LEV_Group::SetYMax(const double &ymax) 
{
  for (size_t i(0);i<m_levs.size();++i) 
    m_levs[i]->SetYMax(ymax);
  m_ymax=ymax; 
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
