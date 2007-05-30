#include "LO_LEV.H"

#define N_C 3.0
#define C_A N_C
#define C_F (N_C*N_C-1.0)/(2.0*N_C)
#define T_R 0.5

using namespace EXTRAXS;
using namespace ATOOLS;

G_GG::G_GG(): 
  LEV_Base(ATOOLS::kf::gluon,ATOOLS::kf::gluon,ATOOLS::kf::gluon) {}

double G_GG::Value(const Vec4D &k1,const Vec4D &q1,
		   const Vec4D &k2,const Vec4D &q2) const
{
  return 2.0*C_A;
}

double G_GG::MajorValue(const Vec4D &k1,const Vec4D &q1,
			const Vec4D &k2,const Vec4D &q2) const
{
  return 2.0*C_A;
}

double G_GG::MajorIntegral(const ATOOLS::Flavour &fl)
{
  if (!fl.IsGluon()) return 0.0;
  return 2.0*C_A;
}

Q_GQ::Q_GQ(const ATOOLS::Flavour &q): 
  LEV_Base(q,ATOOLS::kf::gluon,q) {}

double Q_GQ::Value(const Vec4D &k1,const Vec4D &q1,
		   const Vec4D &k2,const Vec4D &q2) const
{
  return 2.0*C_F;
}

double Q_GQ::MajorValue(const Vec4D &k1,const Vec4D &q1,
			const Vec4D &k2,const Vec4D &q2) const
{
  return 2.0*C_F;
}

double Q_GQ::MajorIntegral(const ATOOLS::Flavour &fl)
{
  if (fl!=m_fla) return 0.0;
  return 2.0*C_F;
}

Q_QG::Q_QG(const ATOOLS::Flavour &q): 
  LEV_Base(q,ATOOLS::kf::gluon,q) {}

double Q_QG::Value(const Vec4D &k1,const Vec4D &q1,
		   const Vec4D &k2,const Vec4D &q2) const
{
  double z(q2.PPlus()/q1.PPlus());
  if (z>1) z=q2.PMinus()/q1.PMinus();
  return C_F*z;
}

double Q_QG::MajorValue(const Vec4D &k1,const Vec4D &q1,
			const Vec4D &k2,const Vec4D &q2) const
{
  return C_F;
}

double Q_QG::MajorIntegral(const ATOOLS::Flavour &fl)
{
  if (fl!=m_fla) return 0.0;
  return C_F;
}

G_QQ::G_QQ(const ATOOLS::Flavour &q): 
  LEV_Base(ATOOLS::kf::gluon,q,q.Bar()) {}

double G_QQ::Value(const Vec4D &k1,const Vec4D &q1,
		   const Vec4D &k2,const Vec4D &q2) const
{
  double z(q2.PPlus()/q1.PPlus());
  if (z>1) z=q2.PMinus()/q1.PMinus();
  return T_R*z;
}

double G_QQ::MajorValue(const Vec4D &k1,const Vec4D &q1,
			const Vec4D &k2,const Vec4D &q2) const
{
  return T_R;
}

double G_QQ::MajorIntegral(const ATOOLS::Flavour &fl)
{
  if (!fl.IsGluon()) return 0.0;
  return T_R;
}

