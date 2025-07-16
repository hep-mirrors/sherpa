#include "HADRONS++/Current_Library/VA_P_P.H"

namespace HADRONS {
namespace VA_P_P_FFs {

class ISGW2 : public FormFactor_Base {
// hep-ph/
  double mb;
  double md;
  double betaB2;
  double mBbar;
  double Nf;
  double mXbar;
  double Nfprime;
  double mq;
  double betaX2;
  bool m_excited;
  
  double Getas( double massq, double massx );
public:
  ISGW2(GeneralModel model, double* masses, const Flavour_Vector& flavs,
        std::vector<int>& indices);
  void CalcFFs(ATOOLS::Vec4D p0, ATOOLS::Vec4D p1);
};

ISGW2::ISGW2(GeneralModel model, double* masses, const Flavour_Vector& flavs,
             std::vector<int>& indices) :
  FormFactor_Base(model, masses, flavs, indices), m_excited(false)
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
  case kf_pi:
  case kf_pi_plus:
  case kf_eta:
  case kf_eta_prime_958:
    mq=0.33;
    betaX2=0.406*0.406;
    mXbar=0.75*0.770+0.25*0.14;
    Nfprime=0.0;
    break;
  case kf_K:
  case kf_K_S:
  case kf_K_L:
  case kf_K_plus:
    mq=0.55;
    betaX2=0.44*0.44;
    mXbar=0.75*0.892+0.25*0.49767;
    Nfprime=2.0;
    break;
  case kf_D:
  case kf_D_plus:
    mq=1.82;
    betaX2=0.45*0.45;
    mXbar=0.75*2.01+0.25*1.87;
    Nfprime=3.0;
    break;
  case kf_a_0_980:
  case kf_a_0_980_plus:
  case kf_f_0_980:
  case kf_f_0_1370:
    mq=0.33;
    betaX2=0.275*0.275;
    mXbar=(3.0*1.23+0.98+5.0*1.32+3.0*1.26)/12.0;
    Nfprime=0.0;
    m_excited=true;
    break;
  case kf_K_0_star_1430:
  case kf_K_0_star_1430_plus:
    mq=0.55;
    betaX2=0.30*0.30;
    mXbar=(3.0*1.40+1.43+5.0*1.43+3.0*1.27)/12.0;
    Nfprime=2.0;
    m_excited=true;
    break;
  case kf_D_0_star:
  case kf_D_0_star_plus:
    mq=1.82;
    betaX2=0.33*0.33;
    mXbar=(3.0*2.49+2.40)/4.0;
    Nfprime=3.0;
    m_excited=true;
    break;
  case kf_D_s0_star:
    mq=1.82;
    betaX2=0.41*0.41;
    mXbar=(3.0*2.54+2.46)/4.0;
    Nfprime=3.0;
    m_excited=true;
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
    break;
  }

  mb      = model("ISGW2_mb",mb);
  md      = model("ISGW2_md",md);
  betaB2  = model("ISGW2_betaB2",betaB2);
  mBbar   = model("ISGW2_mBbar",mBbar);
  Nf      = model("ISGW2_Nf",Nf);

  mq      = model("ISGW2_mq",mq);
  betaX2  = model("ISGW2_betaX2",betaX2);
  mXbar   = model("ISGW2_mXbar",mXbar);
  Nfprime = model("ISGW2_Nfprime",Nfprime);

  m_excited = bool(model("excited",m_excited+0.5));
}

void ISGW2::CalcFFs( Vec4D p1, Vec4D p2 )
{
  Vec4D q = p1 - p2;
  double q2 = q*q;

  double mB = mb + md;
  double mX = mq + md;
  double muplus = 1.0/(1.0/mq+1.0/mb);
  double betaBX2= 0.5*(betaB2+betaX2);
  double tm  = (m_m0-m_m1)*(m_m0-m_m1);
//     if ( q2>tm ) { q2=0.99*tm; }

  double muqm = 0.1;
  double r2  = 3.0/(4.0*mb*mq)+3.0*md*md/(2.0*mBbar*mXbar*betaBX2) +
    16.0/(mBbar*mXbar*(33.0-2.0*Nfprime))*log(Getas(muqm,muqm)/Getas(mq,mq));

  double fppfm, fpmfm;
  if(m_excited) {
    double F5 = sqrt(mX/mB)*pow(sqrt(betaX2*betaB2)/betaBX2,2.5) /
	(pow((1.0+r2*(tm-q2)/18.0),3.0));

    double F5plus  = F5*pow(mBbar/mB,-0.5)*pow(mXbar/mX,0.5);
    double F5minus = F5*pow(mBbar/mB,0.5)*pow(mXbar/mX,-0.5);

    fppfm = -sqrt(2.0/3.0)*md/sqrt(betaB2)*F5plus;
    fpmfm = sqrt(2.0/3.0)*md*mB/sqrt(betaB2)/mX*F5minus;
  }
  else {
    double aI = -6.0/(33.0-2.0*Nf);
    double Cji = pow(Getas(mb,mb)/Getas(mq,mq),aI);
    double zji = mq/mb;
    double gammaji = 2.0*zji*log(1.0/zji)/(1.0-zji) - 2.0;
    double chiji = -1.0-gammaji/(1.0-zji);
    double betaji_plus  = gammaji-2.0/3.0*chiji;
    double betaji_minus = gammaji+2.0/3.0*chiji;
    double Rplus = Cji *(1.0 + betaji_plus*Getas( mq,sqrt(mb*mq) )/M_PI);
    double Rminus = Cji *(1.0 + betaji_minus*Getas( mq,sqrt(mb*mq) )/M_PI);
    double F3 = sqrt(mX/mB)*pow(sqrt(betaX2*betaB2)/betaBX2,1.5)/
      sqr(1.0+r2*(tm-q2)/12.0);
    double F3plus  = F3 * pow(mBbar/mB,-0.5) * pow(mXbar/mX,0.5);
    double F3minus = F3 * pow(mBbar/mB,0.5) * pow(mXbar/mX,-0.5);
    fppfm = (2.0-mX/mq*(1.0-md*mq*betaB2/(2.0*muplus*mX*betaBX2)))*F3plus*Rplus;
    fpmfm = mB/mq*(1.0-md*mq*betaB2/(2.0*muplus*mX*betaBX2))*F3minus*Rminus;
  }

  m_fplus = (fppfm + fpmfm)/2.0;
  double Fminus = (fppfm - fpmfm)/2.0;
  m_f0    = (Fminus/((m_m0*m_m0-m_m1*m_m1)/q2))+(m_fplus);

  m_calced = true;
}

double ISGW2::Getas( double massq, double massx )
{
  double pi = std::acos(-1.0);
  double lqcd2 = 0.04;
  double nflav = 4.0;
  double temp = 0.6;
  
  if ( massx > 0.6 ) {
    if ( massq < 1.85 ) nflav = 3.0;
    
    temp = 12.0*pi / ( 33.0 - 2.0*nflav) / log( massx*massx/lqcd2);
  }
  return temp;
}

}
} // namespace HADRONS
