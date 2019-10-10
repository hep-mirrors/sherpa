#include "AMEGIC++/DipoleSubtraction/FF_DipoleSplitting.H"
#include "AMEGIC++/Main/ColorSC.H"

#include "ATOOLS/Org/My_Limits.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/Poincare.H"
#include "AMEGIC++/DipoleSubtraction/DipoleSplitting_Base.H"
#include "MCATNLO/Showers/Kinematics_Base.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"

using namespace ATOOLS;
using namespace AMEGIC;
using namespace std;

double Lambda(const double& x,const double& y,const double& z)
{
  return sqr(x)+sqr(y)+sqr(z)-2.*(x*y+x*z+y*z);
}

double Vrel(const Vec4D& p1, const Vec4D& p2)
{
  double m21=(p1+p2).Abs2();
  double m1=p1.Abs2();
  double m2=p2.Abs2();
  return sqrt(Lambda(m21,m1,m2))/(m21-m1-m2);
}

void FF_DipoleSplitting::SetMomenta(const Vec4D* mom)
{
  m_mom.clear();
  for(int i=0;i<=m_m;i++) m_mom.push_back(mom[i]);

  m_pi = mom[m_i];
  m_pj = mom[m_j];
  m_pk = mom[m_k];

  switch(m_dipolecase){
    case IDa:
    case IDb:{
      if (m_dipolecase==IDa){
        if(m_j==5) m_n = mom[0]+mom[1]-mom[5]-mom[6]-mom[3] ; // b emitts
        if(m_j==6) m_n = mom[0]+mom[1]-mom[5]-mom[6]-mom[2] ; // bbar emitts
      }
      if (m_dipolecase==IDb){
        if(m_j==5) m_n = mom[0]+mom[1]-mom[5]-mom[2] ; // b emitts
        if(m_j==6) m_n = mom[0]+mom[1]-mom[6]-mom[3] ; // bbar emitts
      }
      m_vi_tilde = m_pi*m_pj/(m_pj*m_n);
      m_a = m_vi_tilde;                                     // alpha_min
      m_zain = m_pj*m_n / ((m_pj+m_pi)*m_n);
      m_viab = m_pj*m_pk / (m_pi*(m_pj+m_pk));
      m_xain = (m_pj-m_pi)*m_n / (m_pj*m_n);
    
      m_ptij = 1./m_zain * m_pj;
      /* trafo of spectator(s) */
      ATOOLS::Vec4D Ka = m_n - mom[m_i];
      ATOOLS::Vec4D Katilde = m_n - (1-m_xain)*mom[m_j];
      std::vector<size_t> kj;
      if(m_dipolecase==IDa){
        if(m_j==5) kj = {2};
        if(m_j==6) kj = {3};
      }
      if(m_dipolecase==IDb){
        if(m_j==5) kj = {3,6};
        if(m_j==6) kj = {2,5};
      }
      int i=0;
      for(auto j: kj){
        m_pIDspec[i] = mom[j] - 2.*(Ka+Katilde)*mom[j]/(Ka+Katilde).Abs2() * (Ka+Katilde)
                             + 2.*Ka*mom[j]/Ka.Abs2() * Katilde;
        i++;
      }

      /* correct v_iab to vtilde_iab in IDb case*/
      if (m_dipolecase==IDb){
        m_viab = m_pj*m_pIDspec[1] / (m_pi*(m_pj+m_pIDspec[1]));
      }
    
      if(m_es == 1){
        /* calculate kt2 (1st way) */
        double kt2 = ((m_pj*m_pi)*(m_pk*m_pi)) / (m_pj*m_pk);

        /* calculate kt2 (2nd way) */
        const ATOOLS::Vec4D pw = m_n - m_pi;
        m_Q2 = (pw+m_pj+m_pi).Abs2();
        const double mw2 = sqr(Flavour(24).Mass());
        const ATOOLS::Vec4D p_wtilde = m_pIDspec[0];
              ATOOLS::Vec4D p_minus  = p_wtilde - mw2/(m_Q2-mw2)*m_ptij;
        const double paipb           = m_ptij*m_pk;
        const double alpha           = m_pk*p_minus / (m_ptij*p_minus);
        const ATOOLS::Vec4D p        = m_ptij+p_minus;
        Poincare bst(p);
        bst.Boost(m_pi);
        bst.Boost(m_pk);
        bst.Boost(m_ptij);
        Vec3D pk_perp = Vec3D(m_pk) - Vec3D(m_pk)*Vec3D(m_ptij)/Vec3D(m_ptij).Sqr()*Vec3D(m_ptij);
        Vec3D pi_perp = Vec3D(m_pi) - Vec3D(m_pi)*Vec3D(m_ptij)/Vec3D(m_ptij).Sqr()*Vec3D(m_ptij);
        const double cos_phi_ib = pi_perp*pk_perp/(pi_perp.Abs()*pk_perp.Abs());
        DEBUG_VAR(acos(cos_phi_ib));
        bst.BoostBack(m_pi);
        bst.BoostBack(m_pk);
        bst.BoostBack(m_ptij);
        /* checked that m_kt2 and kt2 agree */
        m_kt2  = p_nlomc?p_nlomc->KT2(*p_subevt,m_zain,m_vi_tilde,m_Q2,paipb,alpha,cos_phi_ib):kt2;
      }
      if(m_es == 2){
        const ATOOLS::Vec4D pw = m_n - m_pi;
        m_Q2 = (pw+m_pj+m_pi).Abs2();
        double kt2 = m_pj*m_pi;
        /* checked that m_kt2 and kt2 agree */
        m_kt2 = p_nlomc?p_nlomc->KT2(*p_subevt,m_zain,m_vi_tilde,m_Q2):kt2;
      }

      switch (m_ftype) {
      case spt::q2qg:
        THROW(fatal_error, "g>qg should not appear, but g>gq");
        break;
      case spt::q2gq:
        m_sff = 2.*m_viab/m_zain - (1+m_zain);
        m_av  = m_sff;
        break;
      case spt::g2qq:
        m_sff = 1.;
        m_av  = m_sff;
        //THROW(fatal_error, "g>qq should not appear"); 
        //but does appear in construction, but is discarded right afterwards
        msg_Debugging() << "Must not appear (except in dipole construction)\n";
        break;
        break;
      case spt::g2gg:
        THROW(fatal_error, "g>gg in final-state not possible with psuedo-dipoles");
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
      if (m_kt2<(p_nlomc?p_nlomc->KT2Min(0):0.0)) m_av=1.0;
      break;
    }
    case CS:{ 
      m_yijk = m_pi*m_pj/(m_pi*m_pj+m_pj*m_pk+m_pk*m_pi);
      m_a = m_yijk;
    
      m_ptk  = 1./(1.-m_yijk)*m_pk;
      m_ptij = m_pi+m_pj-m_yijk/(1.-m_yijk)*m_pk;
    
      m_zi   = (m_pi*m_ptk)/(m_ptij*m_ptk);
      m_zj   = 1.-m_zi;
    
      m_Q2 = (m_pi+m_pj+m_pk).Abs2();
      m_kt2  = p_nlomc?p_nlomc->KT2(*p_subevt,m_zi,m_yijk,m_Q2):
        m_Q2*m_yijk*m_zi*m_zj;
    
      double zi(m_zi), zj(m_zj);
      if (m_subtype==1) {
        zi=1.0-(1.0-zi)*(1.0-m_yijk);
        zj=1.0-(1.0-zj)*(1.0-m_yijk);
      }
    //   m_pt1   =     m_zi*m_pi;
    //   m_pt2   = -1.*m_zj*m_pj;
      m_pt1   =     m_zi*m_pi-m_zj*m_pj;
      m_pt2   =     m_ptij;
    
      switch (m_ftype) {
      case spt::q2qg:
        m_sff = 2./(1.-m_zi*(1.-m_yijk))-(1.+zi);
        m_av  = m_sff;
        break;
      case spt::q2gq:
        m_sff = 2./(1.-m_zj*(1.-m_yijk))-(1.+zj);
        m_av  = m_sff;
        break;
      case spt::g2qq:
        m_sff = 1.;
        m_av  = m_sff - zi*(1.-zi) - zj*(1.-zj);
        break;
      case spt::g2gg:
        m_sff = 1./(1.-m_zi*(1.-m_yijk))+1./(1.-m_zj*(1.-m_yijk))-2.;
        m_av  = m_sff + ( zi*(1.-zi) + zj*(1.-zj) ) / 2.;
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
      if (m_kt2<(p_nlomc?p_nlomc->KT2Min(0):0.0)) m_av=1.0;
    }
  }
}

double FF_DipoleSplitting::GetValue()
{
  double h=1.0/(2.*m_pi*m_pj);
  return h*m_fac*m_sff;
}

void FF_DipoleSplitting::CalcDiPolarizations()
{
  switch (m_dipolecase){
  case CS: 
  case IDa: 
  case IDb: 
  case IDin: 
    double zi(m_zi), zj(m_zj);
    if (m_subtype==1) {
      zi=1.0-(1.0-zi)*(1.0-m_yijk);
      zj=1.0-(1.0-zj)*(1.0-m_yijk);
    }
    switch (m_ftype) {
    case spt::q2qg:
    case spt::q2gq:
      return;
    case spt::g2qq:
      CalcVectors(m_pt1,m_pt2,m_sff/(2.*(zi*(1.-zi)+zj*(1.-zj))));
      break;
    case spt::g2gg:
      CalcVectors(m_pt1,m_pt2,-m_sff/(zi*(1.-zi)+zj*(1.-zj)));
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
}


void FF_MassiveDipoleSplitting::SetMomenta(const Vec4D* mom)
{
  m_mom.clear();
  for(int i=0;i<=m_m;i++) m_mom.push_back(mom[i]);

  m_pi = mom[m_i];
  m_pj = mom[m_j];
  m_pk = mom[m_k];
  
  m_Q = m_pi+m_pj+m_pk;
  m_Q2 = m_Q.Abs2();

  m_ptk  = sqrt(Lambda(m_Q2,m_mij2,m_mk2)/Lambda(m_Q2,(m_pi+m_pj).Abs2(),m_mk2))
                *(m_pk-(m_Q*m_pk)/m_Q2*m_Q)
           +(m_Q2+m_mk2-m_mij2)/(2.*m_Q2)*m_Q;
  m_ptij = m_Q-m_ptk;

  m_yijk = m_pi*m_pj/(m_pi*m_pj+m_pj*m_pk+m_pk*m_pi);
  m_yp = 1. - 2.*(sqrt(m_mk2*m_Q2) - m_mk2)/(m_Q2 - m_mi2 - m_mj2 - m_mk2);
  m_a = m_yijk/m_yp;
  if ((m_mij2 && !IsEqual(m_ptij.Abs2(),m_mij2,1.0e-3)) ||
      (m_mk2 && !IsEqual(m_ptk.Abs2(),m_mk2,1.0e-3))) {
    msg_Tracking()<<METHOD<<"(): Kinematics unstable in {\n"
		  <<"  p_i = "<<m_pi<<" "<<sqrt(dabs(m_pi.Abs2()))<<"\n"
		  <<"  p_j = "<<m_pj<<" "<<sqrt(dabs(m_pj.Abs2()))<<"\n"
		  <<"  p_k = "<<m_pk<<" "<<sqrt(dabs(m_pk.Abs2()))<<"\n}"
		  <<std::endl;
    m_a=0.0;
  }

  m_zi   = (m_pi*m_pk)/(m_pi*m_pk+m_pj*m_pk);
  m_zj   = 1.-m_zi;

  m_kt2  = p_nlomc?p_nlomc->KT2(*p_subevt,m_zi,m_yijk,m_Q2):
    2.0*m_pi*m_pj*m_zi*m_zj-sqr(m_zi)*m_mj2-sqr(m_zj)*m_mi2;

  m_vijk = Vrel(m_pi+m_pj,m_pk);
  
  m_zim  = m_zi-0.5*(1.-m_vijk);
  m_zjm  = m_zj-0.5*(1.-m_vijk);
  m_zpm  = 0.;
  if (m_ftype==spt::g2qq || m_ftype==spt::g2gg) {
    m_zpm = sqr((2.*m_mi2+(m_Q2-m_mi2-m_mj2-m_mk2)*m_yijk)/
                (2.*(m_mi2+m_mj2+(m_Q2-m_mi2-m_mj2-m_mk2)*m_yijk)))
            *(1.-sqr(m_vijk*Vrel(m_pi+m_pj,m_pi)));
  }

  m_pt1   =     m_zim*m_pi-m_zjm*m_pj;
  m_pt2   =     m_ptij;

  double zi(m_zi), zj(m_zj);
  if (m_subtype==1) {
    zi=1.0-(1.0-zi)*(1.0-m_yijk);
    zj=1.0-(1.0-zj)*(1.0-m_yijk);
  }
  switch (m_ftype) {
  case spt::q2qg:
    m_sff = 2./(1.-m_zi*(1.-m_yijk))-Vrel(m_ptij,m_ptk)/m_vijk*(1.+zi+m_mij2/(m_pi*m_pj));
    m_av  = m_sff;
    break;
  case spt::q2gq:
    m_sff = 2./(1.-m_zj*(1.-m_yijk))-Vrel(m_ptij,m_ptk)/m_vijk*(1.+zj+m_mij2/(m_pi*m_pj));
    m_av  = m_sff;
    break;
  case spt::g2qq:
    m_sff = (1.-2.*m_kappa*(m_zpm-m_mi2/(m_pi+m_pj).Abs2()))/m_vijk;
    m_av  = m_sff - 2.0 * ( m_zi*m_zj - m_zpm )/m_vijk;
    if (m_subtype==1) m_av = m_sff - ( zi*(1.-zi) + zj*(1.-zj) - 2.*m_zpm )/m_vijk;
    break;
  case spt::g2gg:
    m_sff = 1./(1.-m_zi*(1.-m_yijk))+1./(1.-m_zj*(1.-m_yijk))
            -(2.-m_kappa*m_zpm)/m_vijk;
    m_av  = m_sff + ( m_zi*m_zj - m_zpm )/m_vijk;
    if (m_subtype==1) m_av = m_sff + ( zi*(1.-zi) + zj*(1.-zj) - 2.*m_zpm )/(2.*m_vijk);
    break;
  case spt::s2sg:
    m_sff = 2./(1.-m_zi*(1.-m_yijk))
            -Vrel(m_ptij,m_ptk)/m_vijk*(2.+m_mij2/(m_pi*m_pj));
    m_av  = m_sff;
    break;
  case spt::s2gs:
    m_sff = 2./(1.-m_zj*(1.-m_yijk))
            -Vrel(m_ptij,m_ptk)/m_vijk*(2.+m_mij2/(m_pi*m_pj));
    m_av  = m_sff;
    break;
  case spt::G2Gg:
    m_sff = 2./(1.-m_zi*(1.-m_yijk))
            -Vrel(m_ptij,m_ptk)/m_vijk*(1.+m_zi+m_mij2/(m_pi*m_pj));
    m_av  = m_sff;
    break;
  case spt::G2gG:
    m_sff = 2./(1.-m_zj*(1.-m_yijk))
            -Vrel(m_ptij,m_ptk)/m_vijk*(1.+m_zj+m_mij2/(m_pi*m_pj));
    m_av  = m_sff;
    break;
  case spt::V2Vg:
    if      (m_Vsubmode==0) m_sff = 2./(1.-m_zi*(1.-m_yijk))
                                    - Vrel(m_ptij,m_ptk)/m_vijk
                                      *(2.+m_mij2/(m_pi*m_pj));
    else if (m_Vsubmode==1) m_sff = 2./(1.-m_zi*(1.-m_yijk))
                                    - Vrel(m_ptij,m_ptk)/m_vijk
                                      *(1.+m_zi+m_mij2/(m_pi*m_pj));
    else if (m_Vsubmode==2) m_sff = (m_zi)/(1.-m_zi)
                                    -m_mij2/(m_pi*m_pj);
    m_av  = m_sff;
    break;
  case spt::V2gV:
    if      (m_Vsubmode==0) m_sff = 2./(1.-m_zj*(1.-m_yijk))
                                    - Vrel(m_ptij,m_ptk)/m_vijk
                                      *(2.+m_mij2/(m_pi*m_pj));
    else if (m_Vsubmode==1) m_sff = 2./(1.-m_zj*(1.-m_yijk))
                                    - Vrel(m_ptij,m_ptk)/m_vijk
                                      *(1.+m_zj+m_mij2/(m_pi*m_pj));
    else if (m_Vsubmode==2) m_sff = m_zj/(1.-m_zj)
                                    -m_mij2/(m_pi*m_pj);
    m_av  = m_sff;
    break;
  case spt::none:
    THROW(fatal_error, "Splitting type not set.");
  }
  if (m_kt2<(p_nlomc?p_nlomc->KT2Min(0):0.0)) m_av=1.0;
}

double FF_MassiveDipoleSplitting::GetValue()
{
  double h=1.0/((m_pi+m_pj).Abs2()-m_mij2);
  return h*m_fac*m_sff;
}

void FF_MassiveDipoleSplitting::CalcDiPolarizations()
{
  double zi(m_zi), zj(m_zj);
  if (m_subtype==1) {
    zi=1.0-(1.0-zi)*(1.0-m_yijk);
    zj=1.0-(1.0-zj)*(1.0-m_yijk);
  }
  switch (m_ftype) {
  case spt::q2qg:
  case spt::q2gq:
    return;
  case spt::g2qq:
    CalcVectors(m_pt1,m_pt2,m_sff*m_vijk/(2.*(zi*(1.-zi)+zj*(1.-zj)-2.0*m_zpm)));
    break;
  case spt::g2gg:
    CalcVectors(m_pt1,m_pt2,-m_sff*m_vijk/(zi*(1.-zi)+zj*(1.-zj)-2.0*m_zpm));
    break;
  case spt::s2sg:
  case spt::s2gs:
  case spt::G2Gg:
  case spt::G2gG:
  case spt::V2Vg:
  case spt::V2gV:
    return;
  case spt::none:
    THROW(fatal_error, "Splitting type not set.");
  }
}
