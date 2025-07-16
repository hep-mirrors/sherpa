#include "HADRONS++/Current_Library/VA_P_T.H"

namespace HADRONS {
namespace VA_P_T_FFs {

class ISGW2 : public FormFactor_Base {
// hep-ph/
  double mb;
  double md;
  double betaB2;
  double mBbar;
  double Nf;
  double mq;
  double betaX2;
  double mXbar;
  double Nfprime;
  
  double Getas( double mass );
public:
  ISGW2(GeneralModel model, double* masses, const Flavour_Vector& flavs,
        const std::vector<int>& i);
  void CalcFFs(ATOOLS::Vec4D p0, ATOOLS::Vec4D p1);
};

ISGW2::ISGW2(GeneralModel model, double* masses, const Flavour_Vector& flavs,
             const std::vector<int>& i) :
  FormFactor_Base(model, masses, flavs, i)
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
  case kf_a_2_1320:
  case kf_a_2_1320_plus:
  case kf_f_2_1270:
  case kf_f_2_prime_1525:
    mq=0.33;
    betaX2=0.275*0.275;
    mXbar=(3.0*1.23+0.98+5.0*1.32+3.0*1.26)/12.0;
    Nfprime=0.0;
    break;
  case kf_K_2_star_1430:
  case kf_K_2_star_1430_plus:
    mq=0.55;
    betaX2=0.30*0.30;
    mXbar=(3.0*1.40+1.43+5.0*1.43+3.0*1.27)/12.0;
    Nfprime=2.0;
    break;
  case kf_D_2_star_2460:
  case kf_D_2_star_2460_plus:
    mq=1.82;
    betaX2=0.33*0.33;
    mXbar=(5.0*2.46+3.0*2.42)/8.0;
    Nfprime=3.0;
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

  mb = model("ISGW2_mb",mb);
  md = model("ISGW2_md",md);
  betaB2 = model("ISGW2_betaB2",betaB2);
  mBbar = model("ISGW2_mBbar",mBbar);
  Nf  = model("ISGW2_Nf",Nf);

  mq = model("ISGW2_mq",mq);
  betaX2 = model("ISGW2_betaX2",betaX2);
  mXbar = model("ISGW2_mXbar",mXbar);
  Nfprime = model("ISGW2_Nfprime",Nfprime);
}

void ISGW2::CalcFFs( Vec4D p0, Vec4D p1 )
{
  double q2 = (p0-p1).Abs2();
  double mB = mb + md;
  double mX = mq + md;
  
  double muplus = 1.0/(1.0/mq+1.0/mb);
  double muminus = 1.0/(1.0/mq-1.0/mb);
  double betaBX2= 0.5*(betaB2+betaX2);
  double tm  = (m_m0-m_m1)*(m_m0-m_m1);
//   if (q2>tm) q2 = 0.99*tm;
  double wt  = 1.0+(tm-q2)/(2.0*mBbar*mXbar);
  
  double muqm = 0.1;
  double r2 = 3.0/(4.0*mb*mq)+3.0*sqr(md)/(2.0*mBbar*mXbar*betaBX2)+
    16.0/(mBbar*mXbar*(33.0-2.0*Nfprime))*log(Getas(muqm)/Getas(mq));

  double F5 = sqrt(mX/mB)*pow(sqrt(betaX2*betaB2)/betaBX2,2.5)/
       (pow((1.0+r2*(tm-q2)/18.0),3.0));
  
  double F5h = F5*pow(mBbar/mB,-1.5)*pow(mXbar/mX,-0.5);
  double F5k = F5*pow(mBbar/mB,-0.5)*pow(mXbar/mX,0.5);
  double F5bppbm = F5*pow(mBbar/mB,-2.5)*pow(mXbar/mX,0.5);
  double F5bpmbm = F5*pow(mBbar/mB,-1.5)*pow(mXbar/mX,-0.5);
  
  m_h = md/(sqrt(8.0*betaB2)*mB)*(1.0/mq-md*betaB2/(2.0*muminus*mX*betaBX2))*F5h;
  m_k = md/sqrt(2.0*betaB2)*(1.0+wt)*F5k;
  double bppbm = sqr(md)*betaX2/(sqrt(32.0*betaB2)*mq*mb*mB*betaBX2)*
          (1.0-md*betaX2/(2.0*mB*betaBX2))*F5bppbm;
  double bpmbm = -md/(sqrt(2.0*betaB2)*mb*mX)*
    (1.0-md*mb*betaX2/(2.0*muplus*mB*betaBX2)+md*betaX2/(4.0*mq*betaBX2)*
     (1.0-md*betaX2/(2.0*mB*betaBX2)))*F5bpmbm;

  m_bplus = (bppbm + bpmbm)/2.0;
  m_bminus = (bppbm - bpmbm)/2.0;
  
  m_calced = true;
}

double ISGW2::Getas ( double mass )
{
  double lqcd2 = 0.04;
  double nflav = 4;
  double temp = 0.6;
  
  if ( mass > 0.6 ) {
    if ( mass < 1.85 ) {
      nflav = 3.0;}
    
    temp = 12.0*M_PI / ( 33.0 - 2.0*nflav) /
      log( mass*mass/lqcd2);
  }
  return temp;
}

} // namespace VA_P_T
} // namespace HADRONS
