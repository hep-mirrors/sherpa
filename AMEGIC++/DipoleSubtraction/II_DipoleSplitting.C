#include "AMEGIC++/DipoleSubtraction/II_DipoleSplitting.H"
#include "AMEGIC++/Main/ColorSC.H"

#include "ATOOLS/Org/My_Limits.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

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
  for(int i=2;i<=m_m;i++) m_mom[i]-=2.*(m_mom[i]*KKt/KKt.Abs2()*KKt-m_mom[i]*K/K.Abs2()*Kt);

  m_vi   = (m_pi*m_pj)/(m_pi*m_pk);
  m_a = m_vi;

  m_Q2 = (-m_pi+m_pj-m_pk).Abs2();
  if (m_es==0) {
    m_kt2 = m_Q2*(1.-m_xijk)/m_xijk*m_vi;
  }
  else {
  m_kt2 = m_Q2/m_xijk*m_vi;
  switch (m_ft) {
  case 1:
    m_kt2*=(1.-m_xijk);
    break;
  case 4:
    m_kt2*=(1.-m_xijk);
    break;
  }
  }

//   m_pt1  =    m_pj;
//   m_pt2  =-1.*m_vi*m_pk;
  m_pt1  =    m_pj-m_vi*m_pk;
  m_pt2  =    m_ptij;

  switch (m_ft) {
  case 1:
    m_sff = 2./(1.-m_xijk)-(1.+m_xijk);
    m_av  = m_sff;
    break;
  case 2:
    m_sff = 1.-2.*m_xijk*(1.-m_xijk);
    m_av  = m_sff;
    break;
  case 3:
    m_sff = m_xijk;
    m_av  = m_sff + 2.0*(1.0-m_xijk)/m_xijk;
    break;
  case 4:
    m_sff = m_xijk/(1.-m_xijk)+m_xijk*(1.-m_xijk);
    m_av  = m_sff + (1.0-m_xijk)/m_xijk;
  }
  if (m_kt2<m_k0sqi) m_av=1.0;
}

double II_DipoleSplitting::GetValue()
{
  double h=1.0/(2.*m_pi*m_pj)/m_xijk;  
  switch (m_ft) {
  case 1:
    h*=m_sff;
    return h;
  case 2:
    h*=m_sff;   
    return h;
  case 3:
    return h*m_sff*CSC.TR/CSC.CA;
  case 4:
    h*=2.*m_sff;
    return h;
  }
  return 0.;
}

void II_DipoleSplitting::CalcDiPolarizations()
{
  switch (m_ft) {
  case 1:
  case 2:
    return;
  case 3:
    CalcVectors(m_pt1,m_pt2,-m_sff*m_xijk/(1.-m_xijk)/4.);
    break;
  case 4:
    CalcVectors(m_pt1,m_pt2,-m_sff*m_xijk/(1.-m_xijk)/2.);
    break;
  }
}

void II_MassiveDipoleSplitting::SetMomenta(const Vec4D * mom)
{
  // m_msub=1;
  //  if(m_alpha==1.) m_alpha =0.;
 
  Vec4D pa,pb,k;
  Vec4D Pab;
  double Q2;
 
  m_mom.clear();
  for(int i = 0; i <= m_m ; ++i) m_mom.push_back(mom[i]);
 
  // Translation:
  // i     ---- >    emitter   (a)
  // k     ---- >    spectator (b)
  // j     ---- >    emitted   (k)
 
  pa = mom[i_i];
  pb = mom[i_k];
  k  = mom[i_j];
 
  m_pi         = pa;
  m_pj         = k;
  m_pk         = pb;
 
  // Splitting kinematics
  Pab   = pa + pb - k ;  // momentum transferred
  m_xab = (pa*pb - pa*k - pb*k)/(pa*pb);
 
  double sab(0.), s(0.), lab(0.), vab(0.),
    slab(0.);
 
  sab            = 2.*pa*pb;
  s              = (pa+pb).Abs2();
  vab            = 2.*(k*pa)/sab;
  Q2             = Pab.Abs2();
  lab            = sqr(sab) - 4.*m_ma*m_mb;
  slab           = sqrt(lab);
  m_Q2           = Q2;
 
  double num     = sqrt(lambda(Q2,m_ma,m_mb));
 
  m_ptij         = num*(pa - pb*sab/2./m_mb)/slab +
    pb*(sab*m_xab)/(2.*m_mb);
  m_ptk          = pb;
  m_pt1          = k  - vab*pb;
  m_pt2          = m_ptij;
  m_vi           = vab;
  m_xijk         = m_xab;
  Vec4D Pabt, PPabt;
  Pabt           = m_ptij + m_ptk;
  PPabt          = Pabt + Pab;

  //////////////////
  // Check this //
  //////////////////

  if (m_es==0) {
    m_kt2 = m_Q2*(1.-m_xijk)/m_xijk*m_vi;
  }
  else {
    // This is the kt2 that gets used...doesn't include mass corrections...
    m_kt2 = m_Q2/m_xijk*m_vi;
  }
 
  for(size_t i(0); i <= m_m; ++i){
    m_mom[i]      = m_mom[i] - PPabt*(PPabt*m_mom[i])/(Q2 + Pab*Pabt)
      + 2.*(Pab*m_mom[i])*Pabt/(Q2);
  }
   
  m_a     = m_xab;//1.-m_xab;
  // m_xmin  = 2.*sqrt(m_ma*m_mb)/sab;
  
  // msg->SetPrecision(16);
  msg_Debugging() << "alpha:"   << m_alpha
		  << ", m_xab:" << m_xab
		  << ", Q2t:"   << m_Q2 - m_ma - m_mb
		  << ", sab:"   << sab
		  << ", s:"     << s << std::endl;
 
 
  // set the value
  if(m_ft == spt::qg){ // Pqq
    m_sff = 2./(1.-m_xab) - 1. - m_xab - m_xab*m_ma/(pa*k);
    m_av  = m_sff;
  }
  // else if(m_ft == spt::gq){ //Pqg
  //   if(m_xab >= m_xmin){
  // 		m_sff = 1.-2.*m_xijk*(1.-m_xijk) - m_xab*m_ma/(pa*k);
  // 		m_av  = m_sff;
  // 	}
  // 	else{
  // 		m_sff = m_av = 0.;
  // 	}
  // }
  // else if(m_ft == spt::qq){ // Pgq
  //   if(m_xab >= m_xmin){
  //     m_sff = sqr(m_xab)+sqr(1.-m_xab)+m_xab*m_mk/(pa*k);
  //     m_av  = m_sff;
  //   }
  // else{
  // 	m_sff = m_av = 0.;
  // }
  else
    THROW(not_implemented,"...Only available for q->qg and g->qq splittings");
}
 
double II_MassiveDipoleSplitting::GetValue()
{
  double h = 1./(2.*m_pi*m_pj*m_xab);
  if(m_test) return 1.;
  switch (m_ft) {
  case 1:
  case 2:
    h*=m_sff;
    return h;
  case 3:
    return h*m_sff*CSC.TR/CSC.CA;
  }
  return 0.;
}
  
void II_MassiveDipoleSplitting::CalcDiPolarizations()
{
  switch (m_ft) {
  case spt::qg:
  case spt::gq:
    return;
  case spt::qq:
    CalcVectors(m_pt1,m_pt2,-m_sff*m_xab/(1.-m_xab)/4.);
    break;
  case spt::gg:
    CalcVectors(m_pt1,m_pt2,-m_sff*m_xab/(1.-m_xab)/2.);
    break;
  }
 }
