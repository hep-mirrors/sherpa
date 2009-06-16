#include "HADRONS++/Current_Library/VA_P_T.H"

namespace HADRONS {
namespace VA_P_T_FFs {

class ISGW : public FormFactor_Base {
  double mb;
  double md;
  double betaB2;
  double mq;
  double betaX2;
public:
  ISGW(struct GeneralModel model, double* masses, Flavour* flavs, int* i);
  void CalcFFs(ATOOLS::Vec4D p0, ATOOLS::Vec4D p1);
};

ISGW::ISGW(GeneralModel model, double* masses, Flavour* flavs, int* i) :
  FormFactor_Base(model, masses, flavs, i)
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
  case kf_a_2_1320:
  case kf_a_2_1320_plus:
  case kf_f_2_1270:
  case kf_f_2_prime_1525:
    mq=0.33;
    betaX2=0.27*0.27;
    break;
  case kf_D_2_star_2460:
  case kf_D_2_star_2460_plus:
    mq=1.82;
    betaX2=0.34*0.34;
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

  mb = model("ISGW_mb",mb);
  md = model("ISGW_md",md);
  betaB2 = model("ISGW_betaB2",betaB2);

  mq = model("ISGW_mq",mq);
  betaX2 = model("ISGW_betaX2",betaX2);
}

void ISGW::CalcFFs( Vec4D p0, Vec4D p1 )
{
  double q2 = (p0-p1).Abs2();
  double mB = mb + md;
  double mX = mq + md;
 
  double muplus = 1.0/(1.0/mq+1.0/mb);
  double muminus = 1.0/(1.0/mq-1.0/mb);
  double betaBX2= 0.5*(betaB2+betaX2);

  double tm=(m_m0-m_m1)*(m_m0-m_m1);
  double kapa = 0.7*0.7;

  double F5 = sqrt(mX/mB)*pow(sqrt(betaX2*betaB2)/betaBX2,2.5)*
       exp(-1.0*((md*md*(tm-q2)/(4.0*mB*mX*kapa*betaBX2))));

  m_h      = F5*md/(sqrt(8.0*betaB2)*mB)*(1.0/mq-md*betaB2/(2.0*mX*muminus*betaBX2));
  m_k      = sqrt(2.0/betaB2)*F5*md;
  m_bplus  = -F5*md/(sqrt(8.0*betaB2)*mX*mb)*
    (1.0-md*mb*betaX2/(2.0*muplus*mB*betaBX2)+
     md*mb*betaX2/(4.0*mB*muminus*betaBX2)*(1.0-md*betaX2/(2.0*mB*betaBX2)));
  m_bminus = 0.0;

  m_calced = true;
}

} // namespace VA_P_T
} // namespace HADRONS
