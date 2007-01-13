#include "LEV_Base.H"

#include "Exception.H"

using namespace EXTRAXS;
using namespace ATOOLS;

LEV_Base::~LEV_Base() {}

double LEV_Base::Value(const Vec4D &k1,const Vec4D &q1,
		       const Vec4D &k2,const Vec4D &q2) const
{
  THROW(fatal_error,"Virtual function called.");
  return 0.0;
}

double LEV_Base::MajorValue(const Vec4D &k1,const Vec4D &q1,
			    const Vec4D &k2,const Vec4D &q2) const
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

void LEV_Base::SetKT2Min(const double &kt2min) 
{
  m_kt2min=kt2min; 
}

void LEV_Base::SetKT2Max(const double &kt2max) 
{
  m_kt2max=kt2max; 
}

LEV_Base *LEV_Base::Selected() const
{
  return p_selected;
}

