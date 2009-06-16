#include "HADRONS++/Current_Library/VA_P_V.H"

namespace HADRONS {
namespace VA_P_V_FFs {

class ISGW : public FormFactor_Base {
  double mb;
  double md;
  double betaB2;
  double mq;
  double betaX2;
  double kapa2;

  bool m_excited;
  bool m_prime;
public:
  ISGW(struct GeneralModel model, double* masses, Flavour* flavs, int* i);
  void CalcFFs(ATOOLS::Vec4D p0, ATOOLS::Vec4D p1);
};

  ISGW::ISGW(GeneralModel model, double* masses, Flavour* flavs, int* i) :
  FormFactor_Base(model, masses, flavs,i), m_excited(false), m_prime(false)
{
  switch(m_flavs[p_i[0]].Kfcode()) {
  case kf_B_plus:
  case kf_B:
    mb=5.2;
    md=0.33;
    betaB2=0.41*0.41;
    break;
  default:
    msg_Info()<<"Warning: Initializing ISGW form factor for "
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
    mq=0.33;
    betaX2=0.31*0.31;
    break;
  case kf_D_star_2007:
  case kf_D_star_2010_plus:
    mq=1.82;
    betaX2=0.39*0.39;
    break;
  case kf_h_1_1170:
  case kf_h_1_1380:
  case kf_b_1_1235:
  case kf_b_1_1235_plus:
    mq=0.33;
    betaX2=0.27*0.27;
    m_excited=true;
    break;
  case kf_D_1_2420:
  case kf_D_1_2420_plus:
    mq=1.82;
    betaX2=0.34*0.34;
    m_excited=true;
    break;
  case kf_a_1_1260:
  case kf_a_1_1260_plus:
  case kf_f_1_1285:
  case kf_f_1_1420:
    mq=0.33;
    betaX2=0.27*0.27;
    m_prime=true;
    break;
  case kf_D_1_H:
  case kf_D_1_H_plus:
    mq=1.82;
    betaX2=0.34*0.34;
    m_prime=true;
    break;
  default:
    msg_Info()<<"Warning: Initializing ISGW form factor for particles "
              <<m_flavs[p_i[0]]<<" --> "<<m_flavs[p_i[1]]<<". "
              <<"The parameters have no defaults set for these, so make sure "
              <<"to have the correct parameters set in the decay channel file."
              <<std::endl;
    mq=1.0;
    betaX2=1.0;
    break;
  }

  mb     = model("ISGW_mb",mb);
  md     = model("ISGW_md",md);
  betaB2 = model("ISGW_betaB2",betaB2);

  mq     = model("ISGW_mq",mq);
  betaX2 = model("ISGW_betaX2",betaX2);
  kapa2  = model("ISGW_kapa2",0.7*0.7);

  m_excited   = bool(model("excited",m_excited?1.0:0.0));
  m_prime     = bool(model("prime",m_prime?1.0:0.0));
}

void ISGW::CalcFFs( Vec4D p0, Vec4D p1 )
{
  double q2 = (p0-p1).Abs2();
  
  double mB  = mb+md;
  double mX  = mq+md;
  double muminus  = 1.0/(1.0/mq-1.0/mb);
  double muplus  = 1.0/(1.0/mq+1.0/mb);
  double betaBX2 = 0.5*(betaB2+betaX2);

  double tm = sqr(m_m0-m_m1);

  if(m_excited) {
     // 1 ^1 P_1
    double F5 = sqrt(mX/mB)*pow(sqrt(betaX2*betaB2)/betaBX2,2.5)*
      exp(-1.0*((md*md*(tm-q2)/(4.0*mB*mX*kapa2*betaBX2))));
    double v = F5*mB*sqrt(betaB2)/(4.0*sqrt(2.0)*mb*mq*mX);
    double r = F5*mB*sqrt(betaB2)/(sqrt(2.0)*muplus);
    double splus = F5*md/(sqrt(2.0*betaB2)*mB)*(1.0+mb/(2.0*muminus)-
                     mb*mq*md*betaB2/(4.0*muplus*muminus*mX*betaBX2));
    double sminus = 0.0;

    m_A1 = r/(m_m0+m_m1);
    m_A2 = -1.0*splus*(m_m0+m_m1);
    m_A3 = (m_m0+m_m1)/(2.0*m_m1)*m_A1
	    - (m_m0-m_m1)/(2.0*m_m1)*m_A2;
    m_A0 = m_A3 + q2*sminus/(2.0*m_m1);
    m_V  = v*(m_m0+m_m1);
  }
  else if(m_prime) {
     // 1 ^3 P_1
    double F5 = sqrt(mX/mB)*pow(sqrt(betaX2*betaB2)/betaBX2,2.5)*
      exp(-1.0*((md*md*(tm-q2)/(4.0*mB*mX*kapa2*betaBX2))));
    double l = -F5*mB*sqrt(betaB2)*(1.0/muminus+md*(tm-q2)/(2.0*mB* kapa2*betaB2)*
                                    (1.0/mq-md*betaB2/(2.0*muminus*mX*betaBX2)));
    double q = F5*md/(2.0*mX*sqrt(betaB2));
    double cplus = F5*md*mb/(4.0*mB*sqrt(betaB2)*muminus)*
      (1.0-md*mq*betaB2/(2.0*mX*muminus*betaBX2));
    double cminus = 0.0;

    m_A1 = l/(m_m0+m_m1);
    m_A2 = -1.0*cplus*(m_m0+m_m1);
    m_A3 = (m_m0+m_m1)/(2.0*m_m1)*m_A1 - (m_m0-m_m1)/(2.0*m_m1)*m_A2;
    m_A0 = m_A3 + q2*cminus/(2.0*m_m1);
    m_V  = q*(m_m0+m_m1);
  }
  else {
     // 1 ^3 S_1
    double F3 = sqrt(mX/mB)*pow(sqrt(betaX2*betaB2)/betaBX2,1.5)*
      exp(-1.0*((md*md*(tm-q2)/(4.0*mB*mX*kapa2*betaBX2))));

    double f  = 2.0*mB*F3;
    double g  = 0.5*F3*(1.0/mq - md*betaB2/(2.0*muminus*mX*betaBX2));
    double aplus = -F3/(2.0*mX)*(1.0+md/mb*(betaB2-betaX2)/(betaB2+betaX2)-
                                 sqr(md)*sqr(betaX2)/(4.0*muminus*mB*sqr(betaBX2)));
    double aminus = 0.0;
    
    m_A1 = f/(m_m0+m_m1);
    m_A2 = -aplus*(m_m0+m_m1);
    m_A3 = (m_m0+m_m1)/(2.0*m_m1)*m_A1 - (m_m0-m_m1)/(2.0*m_m1)*m_A2;
    m_A0 = m_A3 + q2*aminus/(2.0*m_m1);
    m_V  = g*(m_m0+m_m1);
  }
  m_calced = true;
}

} // namespace VA_P_V
} // namespace HADRONS
