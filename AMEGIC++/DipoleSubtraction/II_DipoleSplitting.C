#include "AMEGIC++/DipoleSubtraction/II_DipoleSplitting.H"
#include "AMEGIC++/Main/ColorSC.H"

#include "ATOOLS/Org/My_Limits.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"

using namespace ATOOLS;
using namespace AMEGIC;
using namespace std;

void II_DipoleSplitting::SetMomenta(const Vec4D *mom)
{
  m_mom.clear();
  for(int i=0;i<=m_m;i++) m_mom.push_back(mom[i]);

  m_pi = mom[m_i];
  m_pj = mom[m_j];
  m_pk = mom[m_k];

  m_xijk = (m_pk*m_pi-m_pi*m_pj-m_pj*m_pk)/(m_pk*m_pi);

  m_ptk  = m_pk;
  m_ptij = m_xijk*m_pi;

  Vec4D K  = m_pi-m_pj+m_pk;
  Vec4D Kt = m_ptij+m_pk;
  Vec4D KKt = K+Kt;
  for(int i=2;i<=m_m;i++)
    m_mom[i]-=2.*(m_mom[i]*KKt/KKt.Abs2()*KKt-m_mom[i]*K/K.Abs2()*Kt);

  m_vi   = (m_pi*m_pj)/(m_pi*m_pk);
  m_a = m_vi;

  m_Q2 = (-m_pi+m_pj-m_pk).Abs2();
  m_kt2  = p_nlomc?p_nlomc->KT2(*p_subevt,m_xijk,m_vi,m_Q2):
    m_Q2*(1.-m_xijk-m_vi)/m_xijk*m_vi;

  double zijk(m_xijk);
  if (m_subtype==1) zijk=m_xijk+m_vi;
//   m_pt1  =    m_pj;
//   m_pt2  =-1.*m_vi*m_pk;
  m_pt1  =    m_pj-m_vi*m_pk;
  m_pt2  =    m_ptij;

  switch (m_ftype) {
  case spt::q2qg:
    m_sff = 2./(1.-m_xijk)-(1.+zijk);
    m_av  = m_sff;
    break;
  case spt::q2gq:
    m_sff = 1.-2.*zijk*(1.-zijk);
    m_av  = m_sff;
    break;
  case spt::g2qq:
    m_sff = zijk;
    m_av  = m_sff + 2.0*(1.0-m_xijk)/m_xijk;
    if (m_subtype==1) m_av += 2.0*(zijk/(sqr(zijk)+m_vi*(1.0-zijk))-1.0/m_xijk);
    if (m_subtype==2) m_av += 2.0*(1.0/(m_xijk+m_vi)-1.0/m_xijk);
    break;
  case spt::g2gg:
    m_sff = m_xijk/(1.-m_xijk)+zijk*(1.-zijk);
    m_av  = m_sff + (1.0-m_xijk)/m_xijk;
    if (m_subtype==1) m_av += zijk/(sqr(zijk)+m_vi*(1.0-zijk))-1.0/m_xijk;
    if (m_subtype==2) m_av += 1.0/(m_xijk+m_vi)-1.0/m_xijk;
    break;
  case spt::none:
    THROW(fatal_error, "Splitting type not set.");
  case spt::s2sg:
  case spt::s2gs:
  case spt::G2Gg:
  case spt::G2gG:
  case spt::V2Vg:
  case spt::V2gV:
    THROW(fatal_error, "DipoleSplitting can not handle splitting type "
        + ToString(m_ftype) + ".");
  }
  if (m_kt2<(p_nlomc?p_nlomc->KT2Min(1):0.0)) m_av=1.0;
}

double II_DipoleSplitting::GetValue()
{
  double h=1.0/(2.*m_pi*m_pj)/m_xijk;
  return h*m_fac*m_sff;
}

void II_DipoleSplitting::CalcDiPolarizations()
{
  double tc((1.-m_xijk)/m_xijk);
  if (m_subtype==1) tc+=(m_xijk+m_vi)/(sqr(m_xijk+m_vi)+m_vi*(1.0-m_xijk-m_vi))-1.0/m_xijk;
  if (m_subtype==2) tc+=1.0/(m_xijk+m_vi)-1.0/m_xijk;
  switch (m_ftype) {
  case spt::q2qg:
  case spt::q2gq:
    return;
  case spt::g2qq:
    CalcVectors(m_pt1,m_pt2,-m_sff/tc/4.);
    break;
  case spt::g2gg:
    CalcVectors(m_pt1,m_pt2,-m_sff/tc/2.);
    break;
  case spt::none:
    THROW(fatal_error, "Splitting type not set.");
  case spt::s2sg:
  case spt::s2gs:
  case spt::G2Gg:
  case spt::G2gG:
  case spt::V2Vg:
  case spt::V2gV:
    THROW(fatal_error, "DipoleSplitting can not handle splitting type "
        + ToString(m_ftype) + ".");
  }
}
