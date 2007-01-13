#include "LO_LEV.H"

#define N_C 3.0

using namespace EXTRAXS;
using namespace ATOOLS;

G_GG::G_GG(): 
  LEV_Base(ATOOLS::kf::gluon,ATOOLS::kf::gluon,ATOOLS::kf::gluon) {}

double G_GG::Value(const Vec4D &k1,const Vec4D &q1,
		   const Vec4D &k2,const Vec4D &q2) const
{
  return 2.0*N_C;
}

double G_GG::MajorValue(const Vec4D &k1,const Vec4D &q1,
			const Vec4D &k2,const Vec4D &q2) const
{
  return 2.0*N_C;
}

double G_GG::MajorIntegral()
{
  return 2.0*N_C;
}

