#include "AMEGIC++/DipoleSubtraction/FF_DipoleSplitting.H"
#include "AMEGIC++/Main/ColorSC.H"

#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/My_Limits.H"
#include "ATOOLS/Math/Poincare.H"
#include "PHASIC++/Channels/Antenna_Kinematics.H"
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
  if (m_subtype==subscheme::Alaric) // && m_ftype==spt::soft)
    return SetMomentaAlaric(mom);

  m_mom.clear();
  for(int i=0;i<=m_m;i++) m_mom.push_back(mom[i]);

  m_pi = mom[m_i];
  m_pj = mom[m_j];
  m_pk = mom[m_k];


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
  if (m_subtype==subscheme::Dire) {
    zi=1.0-(1.0-zi)*(1.0-m_yijk);
    zj=1.0-(1.0-zj)*(1.0-m_yijk);
  }
  m_pt1   =     m_zi*m_pi-m_zj*m_pj;
  m_pt2   =     m_ptij;

  switch (m_ftype) {
  case spt::soft:
    m_sff = 0;
    m_av = 0;
    break;
  case spt::q2qg:
    m_sff = 2./(1.-m_zi*(1.-m_yijk))-(1.+zi); //< CS eq. 5.7
    m_sff *= zi;
    if (m_subtype==subscheme::Alaric) m_sff = zi*(1.-zi);
    m_av  = m_sff;
    break;
  case spt::q2gq:
    m_sff = 2./(1.-m_zj*(1.-m_yijk))-(1.+zj); //< CS eq. 5.7
    m_sff *= zi;
    if (m_subtype==subscheme::Alaric) m_sff = zi*zi;
    m_av  = m_sff;
    break;
  case spt::g2qq:
    m_sff = 1.;
    m_sff *= zi;
    m_av  = m_sff - zi*(zi*(1.-zi) + zj*(1.-zj));
    break;
  case spt::g2gg:
    m_sff = 1./(1.-m_zi*(1.-m_yijk))+1./(1.-m_zj*(1.-m_yijk))-2.;
    m_sff *= zi;
    m_av  = m_sff + zi*( zi*(1.-zi) + zj*(1.-zj) ) / 2.; //< CS eq. 5.9
    if (m_subtype==subscheme::Alaric) {
      m_sff = zi;
      m_av = m_sff - zi*(zi*(1.-zi));
    }
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
  DEBUG_VAR(m_ftype<<" "<<m_sff<<" "<<m_av);
}

void FF_DipoleSplitting::SetMomentaAlaric(const ATOOLS::Vec4D* mom) {
  DEBUG_FUNC("");
  Vec4D n;
  if(m_ftype==spt::soft) {
    PHASIC::Ant_Args ff;
    Cluster_Amplitude* ampl = Cluster_Amplitude::New();
    ampl->SetNIn(2);
    for (int i = 0; i <= m_m; ++i) {
      ff.m_p.push_back(mom[i]);
      ampl->CreateLeg(i<2?-mom[i]:mom[i],i<2?p_subevt->p_real->p_fl[i].Bar():p_subevt->p_real->p_fl[i]);
    }

    ff.m_b=p_softrecoil->RecoilTags(ampl,m_i,m_j,m_k);
    PHASIC::ClusterAntenna(ff, m_i, m_j, m_k, 0.);

    m_pi = ff.m_pi;
    m_pj = ff.m_pj;
    m_pk = ff.m_pk;
    m_mom = ff.m_p;

    m_zi = ff.m_z;
    m_zj = 1.-ff.m_z;

    m_yijk = (m_pi*m_pk)/(m_pj*(m_pi+m_pk)); //< corresponds to v_{i,ab} in CS, eq. 5.169
    m_a = m_yijk;

    m_ptk  = ff.m_p[m_k];
    m_ptij = ff.m_pijt;

    m_Q2 = (m_pi+m_pj+m_pk).Abs2();
    // für Fixed_Order erstmal egal // wird mittels Getter gelöst, später
    m_kt2  = p_nlomc?p_nlomc->KT2(*p_subevt,m_zi,m_yijk,m_Q2):
      m_Q2*m_yijk*m_zi*m_zj;

    n = ff.m_n;
    m_pt1   =     (n*m_pi) / (m_pi*m_pj) * m_pj - n;
    m_pt2   =     m_ptij;

    ampl->Delete();
  }
  else {
    m_pi = mom[m_i];
    m_pj = mom[m_j];
    m_pk = mom[m_k];

    Cluster_Amplitude* ampl = Cluster_Amplitude::New();
    ampl->SetNIn(2);
    for (int i = 0; i <= m_m; ++i) {
      ampl->CreateLeg(i<2?-mom[i]:mom[i],i<2?p_subevt->p_real->p_fl[i].Bar():p_subevt->p_real->p_fl[i]);
    }

    Vec4D pij = mom[m_i]+mom[m_j];
    Vec4D K = p_collrecoil->Recoil(ampl,m_i,m_j,m_k);
    std::vector<int> tags = p_collrecoil->RecoilTags(ampl,m_i,m_j,m_k);
    int nk = std::count_if(tags.begin(),tags.end(),[](int t){return t&2;});

    double K2(K.Abs2());
    int mode = 0; // ?? what does it do?
    PHASIC::Kin_Args ff=PHASIC::ClusterFFDipole(0,0,0,K2,mom[m_i],mom[m_j],K,mode);
    if (ff.m_stat<0) {
      msg_Error()<<METHOD<<": Clustering failed in subtraction.\n";
    }

    m_mom.clear();
    if(nk>1) {
      Poincare oldcms(K), newcms(ff.m_pk);
      newcms.Invert();
      for(size_t i(0);i<ampl->Legs().size();++i) {
        if(tags[i]&2) {
          ampl->Leg(i)->SetMom(newcms*(oldcms*ampl->Leg(i)->Mom()));
        }
        m_mom.push_back((i<2?-1.:+1.)*ampl->Leg(i)->Mom());
      }
    }
    else {
      for(size_t i(0);i<ampl->Legs().size();++i) {
        if(tags[i]&2) {
          ampl->Leg(i)->SetMom(ff.m_pk);
        }
        m_mom.push_back((i<2?-1.:+1.)*ampl->Leg(i)->Mom());
      }
    }

    m_ptij = ff.m_pi;
    m_ptk = ampl->Leg(m_k)->Mom();

    n = ff.m_nb;

    m_zi = m_pi*n/((m_pi+m_pj)*n);
    m_zj = 1.-m_zi;

    /// TODO: some of these might need to become member variables
    ///       to be used in CalcDiPolarization
    double vijk = Vrel(pij,m_pk);
    double viji = Vrel(pij,m_pi);
    double vijj = Vrel(pij,m_pj);

    double zim = m_pi*pij/pij.Abs2() * (1+vijk*viji);
    double zjm = m_pj*pij/pij.Abs2() * (1+vijk*vijj);;
    double zpm = 0;

    m_pt1   =     zim*m_pi-zjm*m_pj;
    m_pt2   =     m_ptij;

    /// TODO: correct?
    m_yijk = m_pi*m_pj/(m_pi*m_pj+m_pj*m_pk+m_pk*m_pi);
    m_a = m_yijk;

    ampl->Delete();
  }


  double zi(m_zi), zj(m_zj);
  switch (m_ftype) {
  case spt::soft: {
    double sij(m_pi*m_pj), sik(m_pi*m_pk), skj(m_pj*m_pk);
    double D(sij*(m_pk*n)+skj*(m_pi*n));
    double A(2*sik/(sij*skj));
    A*=sij*skj*(m_pi*n)/D;
    m_sff = A;
    m_av  = m_sff;
    break;
  }
  case spt::q2qg:
    m_sff = zi*(1.-zi);
    m_av  = m_sff;
    break;
  case spt::q2gq:
    m_sff = zi*zi;
    m_av  = m_sff;
    break;
  case spt::g2qq:
    m_sff = 1.;
    m_sff *= zi;
    m_av  = m_sff - zi*(zi*(1.-zi) + zj*(1.-zj));
    break;
  case spt::g2gg:
      m_sff = zi;
      m_av = m_sff - zi*(zi*(1.-zi));
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
  DEBUG_VAR(m_ftype<<" "<<m_sff<<" "<<m_av);
}

double FF_DipoleSplitting::GetValue()
{
  double h=1.0/(2.*m_pi*m_pj);
  return h*m_fac*m_sff;
}

void FF_DipoleSplitting::CalcDiPolarizations()
{
  double zi(m_zi), zj(m_zj);
  if (m_subtype==subscheme::Dire) {
    zi=1.0-(1.0-zi)*(1.0-m_yijk);
    zj=1.0-(1.0-zj)*(1.0-m_yijk);
  }
  switch (m_ftype) {
  case spt::q2qg:
  case spt::q2gq:
    return;
  case spt::soft:
    CalcVectors(m_pt1,m_pt2);
    break;
  case spt::g2qq:
    CalcVectors(m_pt1,m_pt2,m_sff/zi/(2.*(zi*(1.-zi)+zj*(1.-zj))));
    break;
  case spt::g2gg: {
    if(m_subtype!=subscheme::Alaric) {
      CalcVectors(m_pt1,m_pt2,-m_sff/zi/(zi*(1.-zi)+zj*(1.-zj)));
    }
    else {
      double B = -m_sff/zi/(zi*(1.-zi)+zj*(1.-zj));
      CalcVectors(m_pt1,m_pt2,B);
      m_pfactors[0] = 0;
      m_pfactors[1] = -1/B;
    }
    break;
  }
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
  if (m_subtype==subscheme::Dire) {
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
    if (m_subtype==subscheme::Dire) m_av = m_sff - ( zi*(1.-zi) + zj*(1.-zj) - 2.*m_zpm )/m_vijk;
    break;
  case spt::g2gg:
    m_sff = 1./(1.-m_zi*(1.-m_yijk))+1./(1.-m_zj*(1.-m_yijk))
            -(2.-m_kappa*m_zpm)/m_vijk;
    m_av  = m_sff + ( m_zi*m_zj - m_zpm )/m_vijk;
    if (m_subtype==subscheme::Dire) m_av = m_sff + ( zi*(1.-zi) + zj*(1.-zj) - 2.*m_zpm )/(2.*m_vijk);
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
  if (m_subtype==subscheme::Dire) {
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
