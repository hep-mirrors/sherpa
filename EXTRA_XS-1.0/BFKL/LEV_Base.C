#include "LEV_Base.H"

#include "Exception.H"

using namespace EXTRAXS;
using namespace ATOOLS;

LEV_Base::~LEV_Base() {}

double LEV_Base::Value(const double &y,const double &kt) const
{
  THROW(fatal_error,"Virtual function called.");
  return 0.0;
}

double LEV_Base::MajorValue(const double &y,const double &kt) const
{
  THROW(fatal_error,"Virtual function called.");
  return 0.0;
}

double LEV_Base::MajorIntegral()
{
  THROW(fatal_error,"Virtual function called.");
  return 0.0;
}

bool LEV_Base::SelectSplitting(const double &rn)
{
  THROW(fatal_error,"Virtual function called.");
  return false;
}

void LEV_Base::SetYMin(const double &ymin) 
{
  m_ymin=ymin; 
}

void LEV_Base::SetYMax(const double &ymax) 
{
  m_ymax=ymax; 
}

void LEV_Base::SetKT2Min(const double &kt2min) 
{
  m_kt2min=kt2min; 
}

void LEV_Base::SetKT2Max(const double &kt2max) 
{
  m_kt2max=kt2max; 
}
