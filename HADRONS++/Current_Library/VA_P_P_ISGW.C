#include "HADRONS++/Current_Library/VA_P_P.H"

namespace HADRONS {
namespace VA_P_P_FFs {

class ISGW : public FormFactor_Base {
  double mb;
  double md;
  double mq;
  double betaB2;
  double betaX2;
  double kapa2;
  bool   m_excited;
public:
  ISGW(GeneralModel model, double* masses, const Flavour_Vector& flavs,
       std::vector<int>& indices);
  void CalcFFs(ATOOLS::Vec4D p0, ATOOLS::Vec4D p1);
};

  ISGW::ISGW(GeneralModel model, double* masses, const Flavour_Vector& flavs,
             std::vector<int>& indices) :
    FormFactor_Base(model, masses, flavs, indices), m_excited(false)
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
  case kf_pi:
  case kf_pi_plus:
  case kf_eta:
  case kf_eta_prime_958:
    mq=0.33;
    betaX2=0.31*0.31;
    break;
  case kf_D:
  case kf_D_plus:
    mq=1.82;
    betaX2=0.39*0.39;
    break;
  case kf_a_0_980:
  case kf_a_0_980_plus:
  case kf_f_0_980:
  case kf_f_0_1370:
    mq=0.33;
    betaX2=0.27*0.27;
    m_excited=true;
    break;
  case kf_D_0_star:
  case kf_D_0_star_plus:
    mq=1.82;
    betaX2=0.34*0.34;
    m_excited=true;
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

  m_excited = bool(model("excited",m_excited?1.0:0.0));
}

void ISGW::CalcFFs(Vec4D p0, Vec4D p1)
{
  Vec4D q = p0 - p1;
  double q2 = q*q;

  double mB = mb + md;
  double mX = mq + md;
  double muplus  = 1.0/(1.0/mq+1.0/mb);
  double muminus = 1.0/(1.0/mq-1.0/mb);
  double betaBX2 = 0.5*(betaB2+betaX2);
  double tm=sqr(m_m0-m_m1);

  if(m_excited) {
    double F5 = sqrt(mX/mB)*pow(sqrt(betaX2*betaB2)/betaBX2,2.5)*
                exp(-1.0*((md*md*(tm-q2)/(4.0*mB*mX*kapa2*betaBX2))));

    m_fplus = F5*md*mq*mb/(sqrt(6.0*betaB2)*mX*muminus);
    double Fminus = 0.0;
    m_f0 = Fminus/((m_m0*m_m0-m_m1*m_m1)/q2)+m_fplus;
  }

  else {
    double F3 = sqrt(mX/mB)*pow(sqrt(betaX2*betaB2)/betaBX2,1.5)*
                exp(-sqr(md)*(tm-q2)/(4.0*mB*mX*kapa2*betaBX2));
    m_fplus = F3*(1.0+(mb/(2.0*muminus))-(mb*mq*md*betaB2/(4.0*muplus*muminus*mX*betaBX2)));
    double Fminus = F3*(1.0-(mB+mX)*(0.5/mq-(md*betaB2/(4.0*muplus*mX*betaBX2))));
    m_f0 = Fminus/((m_m0*m_m0-m_m1*m_m1)/q2)+m_fplus;
  }

  m_calced = true;
}

} // namespace VA_P_V
} // namespace HADRONS
