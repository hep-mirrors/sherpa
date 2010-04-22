#include "HADRONS++/Current_Library/VA_P_V.H"

namespace HADRONS {
namespace VA_P_V_FFs {

class ISGW2 : public FormFactor_Base {
// hep-ph/
  double mb;
  double md;
  double betaB2;
  double mq;
  double betaX2;
  double mBbar;
  double Nf;
  double Cf;
  double mXbar;
  double Nfprime;

  bool m_excited;
  bool m_prime;

  double Getas( double massq, double massx );
public:
  ISGW2(struct GeneralModel model, double* masses, Flavour* flavs, int* i);
  void CalcFFs(ATOOLS::Vec4D p0, ATOOLS::Vec4D p1);
};

ISGW2::ISGW2(GeneralModel model, double* masses, Flavour* flavs, int* i) :
  FormFactor_Base(model, masses, flavs, i), m_excited(false), m_prime(false)
{
  switch(m_flavs[p_i[0]].Kfcode()) {
  case kf_B_plus:
  case kf_B:
    mb=5.2;
    md=0.33;
    betaB2=0.431*0.431;
    mBbar=0.75*5.325+0.25*5.279;
    Nf=4.0;
    break;
  case kf_D_plus:
  case kf_D:
    mb=1.82;
    md=0.33;
    betaB2=0.45*0.45;
    mBbar=1.963;
    Nf=3.0;
    break;
  case kf_D_s_plus:
    mb=1.82;
    md=0.55;
    betaB2=0.56*0.56;
    mBbar=1.968;
    Nf=3.0;
    break;
  case kf_B_s:
    mb=5.2;
    md=0.55;
    betaB2=0.54*0.54;
    mBbar=5.38;
    Nf=4.0;
    break;
  default:
    msg_Info()<<"Warning: Initializing ISGW2 form factor for "
              <<m_flavs[p_i[0]]<<" --> "<<m_flavs[p_i[1]]<<". "
              <<"The parameters have no defaults set for these, so make sure "
              <<"to have the correct parameters set in the decay channel file."
              <<std::endl;
    mb=1.0;
    md=1.0;
    betaB2=1.0;
    break;
  }

  switch(m_flavs[p_i[1]].Kfcode()) {
  case kf_rho_770:
  case kf_rho_770_plus:
  case kf_omega_782:
    if(m_flavs[p_i[0]].Kfcode()==kf_D_plus||m_flavs[p_i[0]].Kfcode()==kf_D)
      Cf=0.889;
    else
      Cf=0.905;
    mq=0.33;
    betaX2=0.299*0.299;
    mXbar=0.75*0.770+0.25*0.14;
    Nfprime=0.0;
    break;
  case kf_K_star_892:
  case kf_K_star_892_plus:
    if(m_flavs[p_i[0]].Kfcode()==kf_D_s_plus) {
      Cf=0.8731;
      mXbar=0.87;
    }
    else {
      Cf=0.928;
      mXbar=0.75*0.892+0.25*0.494;
    }
    mq=0.55;
    betaX2=0.33*0.33;
    Nfprime=2.0;
    break;
  case kf_D_star_2007:
  case kf_D_star_2010_plus:
    Cf=0.989;
    mq=1.82;
    betaX2=0.38*0.38;
    mXbar=0.75*2.01+0.25*1.87;
    Nfprime=3.0;
    break;
  case kf_D_s_star_plus:
    Cf=0.984;
    mq=1.82;
    betaX2=0.49*0.49;
    mXbar=0.75*2.11+0.25*1.97;
    Nfprime=3.0;
    break;
  case kf_h_1_1170:
  case kf_h_1_1380:
  case kf_b_1_1235:
  case kf_b_1_1235_plus:
    mq=0.33;
    betaX2=0.275*0.275;
    mXbar=(3.0*1.123+0.98+5.0*1.32+3.0*1.26)/12.0;
    Nfprime=0.0;
    m_excited=true;
    break;
  case kf_K_1_1270:
  case kf_K_1_1270_plus:
    mq=0.55;
    betaX2=0.30*0.30;
    mXbar=(3.0*1.27+1.43+5.0*1.43+3.0*1.4)/12.0;
    Nfprime=2.0;
    m_excited=true;
    break;
  case kf_phi_1020:
    Cf=0.911;
    mq=0.55;
    betaX2=0.37*0.37;
    mXbar=0.97;
    Nfprime=2.0;
    break;
  case kf_D_s1_2536_plus:
    mq=1.82;
    betaX2=0.41*0.41;
    mXbar=(5.0*2.61+3.0*2.54)/8.0;
    Nfprime=3.0;
    m_excited=true;
    break;
  case kf_D_1_2420:
  case kf_D_1_2420_plus:
    mq=1.82;
    betaX2=0.33*0.33;
    mXbar=(5.0*2.46+3.0*2.42)/8.0;
    Nfprime=3.0;
    m_excited=true;
    break;
  case kf_a_1_1260:
  case kf_a_1_1260_plus:
  case kf_f_1_1285:
  case kf_f_1_1420:
    mq=0.33;
    betaX2=0.275*0.275;
    mXbar=(3.0*1.23+0.98+5.0*1.32+3.0*1.26)/12.0;
    Nfprime=0.0;
    m_prime=true;
    break;
  case kf_K_1_1400:
  case kf_K_1_1400_plus:
    mq=0.55;
    betaX2=0.30*0.30;
    mXbar=(3.0*1.40+1.43+5.0*1.43+3.0*1.27)/12.0;
    Nfprime=2.0;
    m_prime=true;
    break;
  case kf_D_1_H:
  case kf_D_1_H_plus:
    mq=1.82;
    betaX2=0.33*0.33;
    mXbar=(3.0*2.49+2.40)/4.0;
    Nfprime=3.0;
    m_prime=true;
    break;
  case kf_D_s1_H:
    mq=1.82;
    betaX2=0.41*0.41;
    mXbar=(3.0*2.54+2.46)/4.0;
    Nfprime=3.0;
    m_prime=true;
    break;
  default:
    msg_Info()<<"Warning: Initializing ISGW2 form factor for particles "
              <<m_flavs[p_i[0]]<<" --> "<<m_flavs[p_i[1]]<<". "
              <<"The parameters have no defaults set for these, so make sure "
              <<"to have the correct parameters set in the decay channel file."
              <<std::endl;
    mq=1.0;
    betaX2=1.0;
    mXbar=1.0;
    Nfprime=0.0;
    Cf=1.0;
    break;
  }

  mb      = model("ISGW2_mb",mb);
  md      = model("ISGW2_md",md);
  betaB2  = model("ISGW2_betaB2",betaB2);
  mBbar   = model("ISGW2_mBbar",mBbar);
  Nf      = model("ISGW2_Nf",Nf);

  Cf      = model("ISGW2_Cf",Cf);
  mq      = model("ISGW2_mq",mq);
  betaX2  = model("ISGW2_betaX2",betaX2);
  mXbar   = model("ISGW2_mXbar",mXbar);
  Nfprime = model("ISGW2_Nfprime",Nfprime);
  
  m_excited   = bool(model("excited",m_excited?1.0:0.0));
  m_prime     = bool(model("prime",m_prime?1.0:0.0));
}

void ISGW2::CalcFFs( Vec4D p0, Vec4D p1 )
{
  double q2 = (p0-p1).Abs2();
  
  double mB  = mb+md;
  double mX  = mq+md;
  double muplus  = 1.0/(1.0/mq+1.0/mb);
  double muminus  = 1.0/(1.0/mq-1.0/mb);
  double betaBX2 = 0.5*(betaB2+betaX2);

  double tm   = (m_m0-m_m1)*(m_m0-m_m1);
//   if ( q2 > tm ) q2 = 0.99*tm;
  double wt = 1.0+(tm-q2)/(2.0*mBbar*mXbar);
  double muqm = 0.1;
  
  double r2 = 3.0/(4.0*mb*mq)+3.0*md*md/(2.0*mBbar*mXbar*betaBX2) +
    16.0/(mBbar*mXbar*(33.0-2.0*Nfprime))*log(Getas(muqm,muqm)/Getas(mq,mq));

  if(m_excited) {
    double F5 = sqrt(mX/mB)*pow(sqrt(betaX2*betaB2)/betaBX2,2.5) /
      (pow((1.0+r2*(tm-q2)/18.0),3.0));

    double F5v = F5*pow(mBbar/mB,-0.5)*pow(mXbar/mX,-0.5);
    double F5r = F5*pow(mBbar/mB,0.5)*pow(mXbar/mX,0.5);
    double F5sppsm = F5*pow(mBbar/mB,-1.5)*pow(mXbar/mX,0.5);
    double F5spmsm = F5*pow(mBbar/mB,-0.5)*pow(mXbar/mX,-0.5);
    
    double v,r,sp,sm,sppsm,spmsm;
    if(mq == md) {
      v = (mB*sqrt(betaB2)/(4.0*sqrt(2.0)*mb*mq*mX) +
          (wt-1.0)*md/(6.0*sqrt(2.0*betaB2)*mX))*F5v; // (118)
      r = mB*sqrt(betaB2/2.0)*(1.0/muplus+md*mX*sqr(wt-1.0)/(3.0*mq*betaB2))*F5r; // (119)
      sppsm = md/(sqrt(2.0*betaB2)*mB)*(1.0-md/mq+md*betaB2/(2.0*muplus*betaBX2))*F5sppsm; // (120)
      spmsm = md/(sqrt(2.0*betaB2)*mq)*
        ((4.0-wt)/3.0 - md*mq*betaB2/(2.0*mX*muplus*betaBX2))*F5spmsm; // (121)
    }
    else {
      v = -md*F5v/(2.0*sqrt(3.0*betaB2)*mX)*((wt+1.0)/2.0+betaB2*mB/(2.0*md*mq*mb)); // (134)
      r = -2.0*mB*sqrt(betaB2/3.0)*(1.0/mq + mX*md*(wt-1.0)/(2.0*betaB2)*
                      ((wt+1.0)/(2.0*mq)-md*betaB2/(2.0*muminus*mX*betaBX2)))*F5r; // (131)
      sppsm = -sqrt(3.0/betaB2)*md*F5sppsm/(2.0*mB)*
        (1.0 - md/(3.0*mq) - md*betaB2/(3.0*betaBX2)*(1.0/(2.0*muminus)-1.0/muplus)); // (132)
      spmsm = -md*F5spmsm/(2.0*sqrt(3.0*betaB2)*mX)*
        ((2.0-wt)*mX/mq + md*betaB2/betaBX2*(1.0/(2.0*muminus)-1.0/muplus)); // (133)
    }

    sp = (sppsm + spmsm)/2.0;
    sm = (sppsm - spmsm)/2.0;

    m_A1 = r/(m_m0+m_m1);
    m_A2 = -1.0*sp*(m_m0+m_m1);
    m_A3 = (m_m0+m_m1)/(2.0*m_m1)*m_A1
	    - (m_m0-m_m1)/(2.0*m_m1)*m_A2;
    m_A0 = m_A3 + q2*sm/(2.0*m_m1);
    m_V  = v*(m_m0+m_m1);
  }
  else if(m_prime) {
    double F5 = sqrt(mX/mB)*pow(sqrt(betaX2*betaB2)/betaBX2,2.5) /
      (pow((1.0+r2*(tm-q2)/18.0),3.0));

    double F5q = F5*pow(mBbar/mB,-0.5)*pow(mXbar/mX,-0.5);
    double F5l = F5*pow(( mBbar / mB ),0.5)*pow(mXbar/mX,0.5);
    double F5cppsm = F5*pow(mBbar/mB,-1.5)*pow(mXbar/mX,0.5);
    double F5cpmsm = F5*pow(mBbar/mB,-0.5)*pow(mXbar/mX,-0.5);

    double q,l,cp,cm,cppsm,cpmsm;
    if(mq==md) {
      q = -1.0*(md*(5.0+wt)*F5q/(2.0*mX*sqrt(betaB2)*6.0));
      l = -1.0*mB*sqrt(betaB2)*F5l*(1.0/muminus+ ( (md*mX*(wt-1.0)/betaB2)*
          ( (5.0+wt)/(6.0*mq)-(md*betaB2)/(2.0*muminus*mX*betaBX2))));
      cppsm = (-1.0*(md*mX*F5cppsm/(2.0*mq*mB*sqrt(betaB2)))*
              (1.0-(md*mq*betaB2)/(2.0*mX*muminus*betaBX2)));
      cpmsm = 1.0*(md*mX*F5cpmsm/(2.0*mq*mB*sqrt(betaB2)))*
              (((wt+2.0)/3.0)-(md*mq*betaB2)/(2.0*mX*muminus*betaBX2))
              *(mB/mX);
    }
    else {
      q = (1.0-betaB2*mB/(4.0*md*mq*mb))*md/(sqrt(6.0*betaB2)*mX)*F5q; // (138)
      l = sqrt(2.0/3.0)*mB*sqrt(betaB2)*(1.0/(2.0*mq) - 3.0/(2.0*mb) +
          md*mX*(wt-1.0)/betaB2*(1.0/mq-md*betaB2/(2.0*muminus*mX*betaBX2)))*F5l; // (135)
      cppsm = sqr(md)*betaX2*F5cppsm/(sqrt(6.0)*mB*mq*sqrt(betaB2)*betaBX2); // (136)
      cpmsm = -sqrt(2.0/3.0)*md/(mX*sqrt(betaB2))*(1.0+md*betaX2/(2.0*mq*betaBX2))*F5cpmsm;
      // (137)
    }

    cp = (cppsm + cpmsm)/2.0;
    cm = (cppsm - cpmsm)/2.0;

    m_A1 = l/(m_m0+m_m1);
    m_A2 = -1.0*cp*(m_m0+m_m1);
    m_A3 = (m_m0+m_m1)/(2.0*m_m1)*m_A1
	    - (m_m0-m_m1)/(2.0*m_m1)*m_A2;
    m_A0 = m_A3 + q2*cm/(2.0*m_m1);
    m_V  = q*(m_m0+m_m1);
  }
  else {
    double aI = -6.0/(33.0-2.0*Nf);
    double Cji = pow(Getas(mb,mb)/Getas(mq,mq),aI);
    double zji = mq/mb;

    double gammaji = 2.0*zji*log(1.0/zji)/(1.0-zji)-2.0;
    double chiji = -1.0-gammaji/(1.0-zji);
    
    double betaji_g = 2.0/3.0+gammaji;
    double betaji_f = -2.0/3.0+gammaji;
    double betaji_appam = -1.0-chiji+4.0/3.0/(1.0-zji)+
      2.0*(1.0+zji)*gammaji/(3.0*sqr(1.0-zji));
    double betaji_apmam = -4.0/(3.0*(1.0-zji))-chiji+1.0/3.0+
      (1.0-2.0*(1.0+zji)/(3.0*sqr(1.0-zji)))*gammaji;

    double R_g     = Cji*(1.0+betaji_g*Getas(mq,sqrt(m_m0*mq))/M_PI);
    double R_f     = Cji*(1.0+betaji_f*Getas(mq,sqrt(m_m0*mq))/M_PI);
    double R_apmam = Cji*(1.0+betaji_apmam*Getas(mq,sqrt(m_m0*mq))/M_PI);

    double F3=sqrt(mX/mB)*pow(sqrt(betaX2*betaB2)/betaBX2,1.5)/
      ((1.0+r2*(tm-q2)/12.0)*(1.0+r2*(tm-q2)/12.0));
    
    double F3f=sqrt(mBbar/mB)*sqrt(mXbar/mX)*F3;
    double F3g=sqrt(mB/mBbar)*sqrt(mX/mXbar)*F3;
    double F3appam=pow(mBbar/mB,-1.5)*sqrt(mXbar/mX)*F3;
    double F3apmam=sqrt(mB/mBbar)*sqrt(mX/mXbar)*F3;
    
    double f=Cf*mB*(1.0+wt+md*(wt-1.0)/(2.0*muplus))*F3f*R_f;
    double g=0.5*(1.0/mq-md*betaB2/(2.0*muminus*mX*betaBX2))*F3g*R_g;
    double appam=Cji*(md*betaX2*(1.0-md*betaX2/(2.0*mB*betaBX2))/
                      ((1.0+wt)*mq*mb*betaBX2)-
                      betaji_appam*Getas(mq,sqrt(mq*m_m0))/mB/M_PI)*F3appam;
    double apmam=-1.0/mX*(
      mB/mb -
      md*betaX2/(2.0*muplus*betaBX2) +
      wt*md*mB*betaX2*(1.0-md*betaX2/(2.0*mB*betaBX2))/((wt+1.0)*mq*mb*betaBX2)
      )*F3apmam*R_apmam;
    
    double ap=0.5*(appam+apmam);
    double am=0.5*(appam-apmam);  

    m_V  = g*(m_m0+m_m1);
    m_A1 = f/(m_m0+m_m1);
    m_A2 = -1.0*ap*(m_m0+m_m1);
    m_A3 = (m_m0+m_m1)/(2.0*m_m1)*m_A1 - (m_m0-m_m1)/(2.0*m_m1)*m_A2;
    m_A0 = m_A3 + ( (q2*am)/(2.0*m_m1));
  }
  
  m_calced = true;
}

double ISGW2::Getas( double massq, double massx )
{
  double lqcd2 = 0.04;
  double nflav = 4;
  double temp = 0.6;

  if ( massx > 0.6 ) {
    if ( massq < 1.85 ) {
      nflav = 3.0;}

    temp = 12.0*M_PI / ( 33.0 - 2.0*nflav) /
      log( massx*massx/lqcd2);
  }
  return temp;
}

}
} // namespace HADRONS
