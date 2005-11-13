#include "LO_LEV.H"

using namespace EXTRAXS;
using namespace ATOOLS;

G_GG::G_GG(): 
  LEV_Base(ATOOLS::kf::gluon,ATOOLS::kf::gluon,ATOOLS::kf::gluon) {}

double G_GG::Value(const double &y,const double &kt) const
{
  return 1.0;
}

double G_GG::MajorValue(const double &y,const double &kt) const
{
  return 1.0;
}

double G_GG::MajorIntegral()
{
  return 1.0;
}

