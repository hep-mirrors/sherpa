#include "PHASIC++/Process/Process_Info.H"
#include "EXTRA_XS/Main/ME2_Base.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "PHASIC++/Process/Process_Info.H"
#include "PHASIC++/Process/External_ME_Args.H"
#include "MODEL/Main/Running_AlphaQED.H"

using namespace PHASIC;
using namespace ATOOLS;
using namespace MODEL;

namespace EXTRAXS {
  class PionPionGamma :  public ME2_Base {
  public:
    PionPionGamma(const External_ME_Args& args);
    ~PionPionGamma(){};

    double operator()(const ATOOLS::Vec4D_Vector& mom);
    double Initial();
    double Final();
    double IFI();
    double Full();
    inline double Pair(const Vec4D &p1,const Vec4D &p2) {return p1*p2;}
    inline double Den(const double& a, const double &b) {return (a==b?1:1./(a-b));}
    inline Vec4D k(const int &i) {return m_pmom[i-1];}
    double m_me, m_mpi, m_pij[5][5], m_s45;
    double m_betapi, m_s, m_sp;
    double ME2, MP2, m_alpha;
    double S34,T24,T14,S12;
    double S,U,T;
    Vec4D_Vector m_pmom;
  };
}

using namespace EXTRAXS;


PionPionGamma::PionPionGamma(const External_ME_Args& args) : ME2_Base(args)
{
  PRINT_INFO("initialised XS_PionPionGamma");
  Flavour_Vector outflavs = args.m_outflavs;
  // p_bornme = dynamic_cast<ME2_Base*>(PHASIC::Tree_ME2_Base::GetME2(bargs));
  // if (!p_bornme) THROW(fatal_error,"no born me found.");
  const Flavour_Vector fl = args.Flavours();
  double ME = fl[0].Mass();
  double MP = fl[3].Mass();
  ME2 = ME*ME;
  MP2 = MP*MP;
  m_oqcd = 1;
  m_oew  = 2;
  m_me  = outflavs[0].Mass();
  m_mpi = outflavs[2].Mass();
  m_s = sqr(rpa->gen.Ecms());
  m_betapi = sqrt(1-m_mpi*m_mpi/m_s);
}


double PionPionGamma::operator()
(const ATOOLS::Vec4D_Vector& p)
{
  m_pmom = p;
  m_sp = (p[2]+p[3]).Abs2();
  m_s45 = (p[3]+p[4]).Abs2();
  m_betapi = sqrt(1-m_mpi*m_mpi/m_sp);
  S = m_sp;
  S34 = (p[2]+p[3]).Abs2();
  T14 = (p[0]-p[3]).Abs2();
  T24 = (p[1]-p[3]).Abs2();
  T=(p[0]-p[2]).Abs2();;
  U=(p[0]-p[3]).Abs2();
  for (int i = 0; i < p.size(); ++i)
  {
    for (int j = i; j < p.size(); ++j)
    {
      m_pij[i][j]= p[i]*p[j];
      // if(IsZero(m_pij[i][j])){
      //   msg_Error()<<"Divide by Zero in "<<METHOD<<std::endl
      //              <<"p_"<<i<<j<<" = "<< m_pij[i][j] <<std::endl;
      // } 
    }
  }
  m_alpha = (*aqed)(m_sp);
  // PRINT_VAR(Full());
  double res;// = Initial()+Final()+IFI();
  double Q2 = S;
  double s15 = (p[0]+p[4]).Abs2(); 
  double s25 = (p[1]+p[4]).Abs2(); 
  double s35 = (p[2]+p[4]).Abs2();


  double ee2uug_ee = s15*s25*(-2*MP2*MP2*Q2 - MP2*(s15*s15 + s25*s25 + Q2*(-s15 + s25 + 2*s35 - 4*T))
      - Q2*((s15 - s35)*(Q2 + s15 + s25 - s35) + (2*Q2 + 3*s15 + s25 - 2*s35)*T + 2*T*T));
  ee2uug_ee = ee2uug_ee + ME2*(2*Q2*Q2*s15*s25 + 2*MP2*MP2*pow(s15 + s25,2)
    + 2*pow(s15*s15 + s25*T + s15*(s25 - s35 + T),2)
    - 4*MP2*(s15*s15*s15 + s25*s25*T + s15*s15*(2*s25 - s35 +T) + s15*s25*(Q2 + s25 - s35 + 2*T)) 
    + Q2*(2*s15*s15*s15 + 2*s25*s25*T + s15*s15*(7*s25 - 2*s35 + 2*T)
      + s15*s25*(3*s25 - 4*s35 + 8*T)));
  ee2uug_ee = ee2uug_ee + 2*ME2*ME2*ME2*pow(s15 + s25,2.);
  ee2uug_ee = ee2uug_ee + 2*ME2*ME2*(2*MP2*pow(s15 + s25,2) - Q2*(s15*s15 + 3*s15*s25 + s25*s25)
    - 2*(s15 + s25)*(s15*s15 + s25*T + s15*(s25 - s35 + T)));

  // ee2uug_ee = 256 * m_alpha * pow(M_PI*m_alpha,2) / Q2*Q2 / s15 / s15 / s25 / s25* ee2uug_ee;

  res = Full();
  if(IsBad(res)){
     for (int i = 1; i < p.size(); ++i)
    {
        for (int j = i; j < p.size(); ++j)
        {
            msg_Out()<<"pair(k("<<i<<"),("<<j<<")) = "<<Pair(k(i),k(j))<<std::endl;
        }
    }
  }
  return res;
}


double PionPionGamma::Full(){
  return -256*pow(M_PI*m_alpha,3)*(-(MP2*MP2* pow((Den(ME2 - 2*Pair(k(1),k(5)),ME2) - Den(ME2 - 2*Pair(k
 (2),k(5)),ME2))* Den(2*(MP2 + Pair(k(3),k(4))),0) + 2*Den(2*(ME2 + Pair(k(1),k(2))),0)* Den
 (MP2 + 2*Pair(k(3),k(5)),MP2),2)*(ME2 + Pair(k(1),k(2)))) + 2*((Den(ME2 - 2*Pair(k(1),k
 (5)),ME2) - Den(ME2 - 2*Pair(k(2),k(5)),ME2))* Den(2*(MP2 + Pair(k(3),k(4))),0) - 2*Den(2*
 (ME2 + Pair(k(1),k(2))),0)*Den(MP2 + 2*Pair(k(4),k(5)),MP2))*(((Den(ME2 - 2*Pair(k(1),k
 (5)),ME2) - Den(ME2 - 2*Pair(k(2),k(5)),ME2))* Den(2*(MP2 + Pair(k(3),k(4))),0) + 2*Den(2*
 (ME2 + Pair(k(1),k(2))),0)* Den(MP2 + 2*Pair(k(3),k(5)),MP2))*(Pair(k(1),k(4))*(MP2 - Pair(k(1),k
 (3)) - Pair(k(2),k(3)))* Pair(k(2),k(3)) - Pair(k(1),k(3))*(-MP2 + Pair(k(1),k(3)) + Pair(k(2),k
 (3)))*Pair(k(2),k(4)) + Pair(k(1),k(2))*(-MP2 + Pair(k(1),k(3)) + Pair(k(2),k(3)))* Pair(k(3),k
 (4)) + ME2*Pair(k(3),k(4))*(Pair(k(3),k(4)) + Pair(k(3),k(5)))) + Den(2*(MP2 + Pair(k(3),k
 (4))),0)*(-((4*Den(ME2 - 2*Pair(k(1),k(5)),ME2) + 5*Den(ME2 - 2*Pair(k(2),k(5)),ME2))*(MP2*Pair(k
 (1),k(2)) + Pair(k(1),k(3))*(ME2 - 2*Pair(k(2),k(3))) + ME2*Pair(k(2),k(3)))*Pair(k(2),k
 (5))) + ME2*(4*Den(ME2 - 2*Pair(k(1),k(5)),ME2) + 5*Den(ME2 - 2*Pair(k(2),k(5)),ME2))*Pair(k(2),k
 (5))*(Pair(k(3),k(4)) + Pair(k(3),k(5))) + MP2*(Den(ME2 - 2*Pair(k(1),k(5)),ME2) + Den
 (ME2 - 2*Pair(k(2),k(5)),ME2))*(-(Pair(k(1),k(5))*Pair(k(2),k(4))) + Pair(k(1),k(4))*Pair(k(2),k
 (5)) + ME2*(Pair(k(3),k(4)) + Pair(k(3),k(5)))) -(Den(ME2 - 2*Pair(k(1),k(5)),ME2) + Den
 (ME2 - 2*Pair(k(2),k(5)),ME2))*(4*MP2*(MP2*Pair(k(1),k(2)) + ME2*Pair(k(1),k(3))) - 6*Pair(k(1),k
 (3))*(MP2*Pair(k(1),k(2)) + ME2*Pair(k(1),k(3))) +(2*(2*MP2 - 3*Pair(k(1),k(3)))*(ME2 - 2*Pair(k
 (1),k(3))) -(4*ME2 + MP2 - 9*Pair(k(1),k(3)))*Pair(k(1),k(4)))* Pair(k(2),k(3)) + MP2*(MP2*Pair(k
 (1),k(2)) + Pair(k(1),k(3))*(ME2 - 2*Pair(k(2),k(3))) + ME2*Pair(k(2),k(3))) + Pair(k(1),k(4))*
 (-4*MP2*Pair(k(1),k(2)) - 4*ME2*Pair(k(1),k(3)) + pow(Pair(k(2),k(3)),2)) - Pair(k(1),k(3))*
 (-MP2 + Pair(k(1),k(3)) + Pair(k(2),k(3)))* Pair(k(2),k(4)) +(2*Pair(k(1),k(3))*(3*ME2 - 5*Pair(k
 (2),k(3))) + Pair(k(1),k(2))*(5*MP2 + Pair(k(1),k(3)) - Pair(k(2),k(3))) + 4*ME2*Pair(k(2),k
 (3)))*Pair(k(3),k(4)) - ME2*(4*MP2 - 6*Pair(k(1),k(3)) - 4*Pair(k(1),k(4)) + 5*Pair(k(3),k(4)))*
 (Pair(k(3),k(4)) + Pair(k(3),k(5)))))) + pow((Den(ME2 - 2*Pair(k(1),k(5)),ME2) - Den
 (ME2 - 2*Pair(k(2),k(5)),ME2))* Den(2*(MP2 + Pair(k(3),k(4))),0) - 2*Den(2*(ME2 + Pair(k(1),k
 (2))),0)*Den(MP2 + 2*Pair(k(4),k(5)),MP2), 2)*(-(ME2*(MP2*MP2)) + (MP2*Pair(k(1),k(2)) - 2*Pair(k
 (1),k(3))*Pair(k(2),k(3)))*(MP2 - 2*Pair(k(1),k(4)) - 2*Pair(k(2),k(4)) + 2*Pair(k(3),k
 (4))) - 2*ME2*MP2*Pair(k(4),k(5))) + 2*(MP2*pow((Den(ME2 - 2*Pair(k(1),k(5)),ME2) - Den
 (ME2 - 2*Pair(k(2),k(5)),ME2))* Den(2*(MP2 + Pair(k(3),k(4))),0) + 2*Den(2*(ME2 + Pair(k(1),k
 (2))),0)* Den(MP2 + 2*Pair(k(3),k(5)),MP2),2)*Pair(k(1),k(4))* Pair(k(2),k(4)) + Den(ME2 - 2*Pair
 (k(1),k(5)),ME2)* Den(2*(MP2 + Pair(k(3),k(4))),0)*((Den(ME2 - 2*Pair(k(1),k(5)),ME2) - Den
 (ME2 - 2*Pair(k(2),k(5)),ME2))* Den(2*(MP2 + Pair(k(3),k(4))),0) - 2*Den(2*(ME2 + Pair(k(1),k
 (2))),0)* Den(MP2 + 2*Pair(k(4),k(5)),MP2))*(Pair(k(1),k(5))*(MP2*Pair(k(1),k(2)) + Pair(k(1),k
 (3))*(ME2 - 2*Pair(k(2),k(3))) + ME2*Pair(k(2),k(3)) ) - 2*(MP2*Pair(k(1),k(2)) + Pair(k(1),k
 (3))*(ME2 - 2*Pair(k(2),k(3))) + ME2*Pair(k(2),k(3)))*Pair(k(3),k(5)) - ME2*(Pair(k(3),k
 (4)) + Pair(k(3),k(5)))*(-2*Pair(k(1),k(3)) - 2*Pair(k(1),k(4)) - Pair(k(1),k(5)) + 2*(MP2 + Pair
 (k(3),k(4)) + Pair(k(4),k(5)))))) + 2*Den(2*(MP2 + Pair(k(3),k(4))),0)*((Den(ME2 - 2*Pair(k(1),k
 (5)),ME2) - Den(ME2 - 2*Pair(k(2),k(5)),ME2))* Den(2*(MP2 + Pair(k(3),k(4))),0) + 2*Den(2*
 (ME2 + Pair(k(1),k(2))),0)*Den(MP2 + 2*Pair(k(3),k(5)),MP2))*(-(MP2*(Den(ME2 - 2*Pair(k(1),k
 (5)),ME2) + Den(ME2 - 2*Pair(k(2),k(5)),ME2))*(-(Pair(k(1),k(4))*Pair(k(2),k(3))) + Pair(k(1),k
 (2))*(-MP2 + Pair(k(1),k(3)) + Pair(k(2),k(3))) - Pair(k(1),k(3))*Pair(k(2),k(4)) + ME2*(Pair(k
 (3),k(4)) + Pair(k(3),k(5))))) +(Den(ME2 - 2*Pair(k(1),k(5)),ME2) + Den(ME2 - 2*Pair(k(2),k
 (5)),ME2))*((3*MP2 + 4*Pair(k(1),k(2)) - 5*Pair(k(1),k(3)) - 4*Pair(k(1),k(4)))*Pair(k(1),k
 (4))*Pair(k(2),k(3)) + Pair(k(1),k(4))*pow(Pair(k(2),k(3)),2) + 2*Pair(k(1),k(2))*((-MP2 + Pair
 (k(1),k(3)))*(-2*MP2 + 3*Pair(k(1),k(3)) + 2*Pair(k(1),k(4))) +(-2*MP2 + 3*Pair(k(1),k(3)))*Pair(k
 (2),k(3))) - Pair(k(1),k(3))*(-5*MP2 + 7*Pair(k(1),k(3)) + 4*Pair(k(1),k(4)) + Pair(k(2),k
 (3)))*Pair(k(2),k(4)) +(4*Pair(k(1),k(4))*Pair(k(2),k(3)) + Pair(k(1),k(2))*(5*MP2 - 5*(Pair(k
 (1),k(3)) + Pair(k(2),k(3)))) + 6*Pair(k(1),k(3))*Pair(k(2),k(4)))*Pair(k(3),k(4)) - ME2*
 (4*MP2 - 6*Pair(k(1),k(3)) - 4*Pair(k(1),k(4)) + 5*Pair(k(3),k(4)))*(Pair(k(3),k(4)) + Pair(k(3),k
 (5)))) -(4*Den(ME2 - 2*Pair(k(1),k(5)),ME2) + 5*Den(ME2 - 2*Pair(k(2),k(5)),ME2))*Pair(k(2),k(5))*
 (-(Pair(k(1),k(2))*(MP2 - 2*Pair(k(2),k(3)))) + Pair(k(1),k(3))*(-ME2 + Pair(k(2),k(5))) - Pair(k
 (2),k(3))*(Pair(k(1),k(5)) + 2*(Pair(k(2),k(3)) + Pair(k(2),k(5)) - Pair(k(3),k(5)))) + ME2*(Pair
 (k(2),k(3)) + Pair(k(3),k(4)) + Pair(k(3),k(5)))) + Den(ME2 - 2*Pair(k(1),k(5)),ME2)*(-((MP2*Pair
 (k(1),k(2)) -(ME2 + 2*Pair(k(1),k(2)))*Pair(k(1),k(3)) + 2*pow(Pair(k(1),k(3)),2))*Pair(k(1),k
 (5))) - 2*Pair(k(1),k(3))*pow(Pair(k(1),k(5)),2) + Pair(k(1),k(5))*(-ME2 + Pair(k(1),k(5)))*Pair
 (k(2),k(3)) - Pair(k(1),k(3))*Pair(k(1),k(5))*Pair(k(2),k(5)) + 2*(Pair(k(1),k(2))*(MP2 - 2*Pair(k
 (1),k(3))) + 2*pow(Pair(k(1),k(3)),2) +(ME2 - Pair(k(1),k(5)))*Pair(k(2),k(3)) + Pair(k(1),k
 (3))*(-ME2 + 3*Pair(k(1),k(5)) + Pair(k(2),k(5))))* Pair(k(3),k(5)) - 4*Pair(k(1),k(3))*pow(Pair
 (k(3),k(5)),2) + ME2*(Pair(k(3),k(4)) + Pair(k(3),k(5)))*(-2*Pair(k(1),k(3)) - 2*Pair(k(1),k
 (4)) - Pair(k(1),k(5)) + 2*(MP2 + Pair(k(3),k(4)) + Pair(k(4),k(5)))))) + 2*pow(Den(2*
 (MP2 + Pair(k(3),k(4))),0),2)*(pow(Den(ME2 - 2*Pair(k(1),k(5)),ME2) + Den(ME2 - 2*Pair(k(2),k
 (5)),ME2),2)*(-4*pow(Pair(k(1),k(4)),2)*(4*ME2 + Pair(k(2),k(3))) + pow(Pair(k(1),k(3)),2)*
 (-6*(7*ME2 + Pair(k(1),k(2))) + 7*Pair(k(2),k(4))) - MP2*((4*(ME2 + Pair(k(1),k(2))) + Pair(k(1),k
 (5)))* Pair(k(2),k(3)) + ME2*(25*MP2 + Pair(k(2),k(4)))) -(49*ME2*MP2 + 5*(ME2 + Pair(k(1),k
 (2)))*Pair(k(2),k(3)))* Pair(k(3),k(4)) - 24*ME2*pow(Pair(k(3),k(4)),2) + Pair(k(1),k(4))*
 (39*ME2*MP2 + pow(Pair(k(2),k(3)),2) + MP2*Pair(k(2),k(4)) + 40*ME2*Pair(k(3),k(4)) + 4*Pair(k
 (2),k(3))*(ME2 + MP2 + Pair(k(1),k(2)) + Pair(k(3),k(4)))) + Pair(k(1),k(3))*
 (64*ME2*MP2 + 6*ME2*Pair(k(2),k(3)) + MP2*Pair(k(2),k(3)) - 4*MP2*Pair(k(2),k(4)) - Pair(k(2),k
 (3))*Pair(k(2),k(4)) + Pair(k(1),k(4))*(-52*ME2 - 7*Pair(k(2),k(3)) + 4*Pair(k(2),k
 (4))) + MP2*Pair(k(2),k(5)) +(65*ME2 - 6*Pair(k(2),k(4)))*Pair(k(3),k(4)) + Pair(k(1),k(2))*
 (4*MP2 - 4*Pair(k(1),k(4)) + 6*Pair(k(2),k(3)) + 5*Pair(k(3),k(4))))) + ME2*(-(pow(4*Den
 (ME2 - 2*Pair(k(1),k(5)),ME2) + 5*Den(ME2 - 2*Pair(k(2),k(5)),ME2),2)* pow(Pair(k(2),k
 (5)),2)) - pow(Den(ME2 - 2*Pair(k(1),k(5)),ME2),2)* pow(Pair(k(1),k(5)) - 2*Pair(k(3),k
 (5)),2) - 2*Den(ME2 - 2*Pair(k(1),k(5)),ME2)*(4*Den(ME2 - 2*Pair(k(1),k(5)),ME2) + 5*Den
 (ME2 - 2*Pair(k(2),k(5)),ME2))*Pair(k(2),k(5))*(2*(ME2 + Pair(k(1),k(2))) - 3*Pair(k(1),k(5)) - 2*
 (Pair(k(1),k(3)) + Pair(k(2),k(3)) + Pair(k(2),k(5))) + 4*Pair(k(3),k(5))) + MP2*pow(Den
 (ME2 - 2*Pair(k(1),k(5)),ME2) + Den(ME2 - 2*Pair(k(2),k(5)),ME2),2)*Pair(k(4),k(5))) +(Den
 (ME2 - 2*Pair(k(1),k(5)),ME2) + Den(ME2 - 2*Pair(k(2),k(5)),ME2))*(Den(ME2 - 2*Pair(k(1),k
 (5)),ME2)*(Pair(k(1),k(5)) - 2*Pair(k(3),k(5)))*(2*pow(Pair(k(1),k(3)),2) - 8*ME2*Pair(k(1),k
 (4)) + Pair(k(1),k(5))*Pair(k(2),k(3)) + 10*ME2*(MP2 + Pair(k(3),k(4))) + Pair(k(1),k(3))*(-2*Pair
 (k(1),k(2)) + 2*Pair(k(1),k(5)) + 2*Pair(k(2),k(3)) + Pair(k(2),k(5)) - 2*(7*ME2 + Pair(k(3),k
 (5))))) +(4*Den(ME2 - 2*Pair(k(1),k(5)),ME2) + 5*Den(ME2 - 2*Pair(k(2),k(5)),ME2))*Pair(k(2),k
 (5))*(2*pow(Pair(k(2),k(3)),2) - 8*ME2*Pair(k(2),k(4)) + Pair(k(1),k(3))*Pair(k(2),k(5)) + Pair
 (k(2),k(3))*(-14*ME2 - 2*Pair(k(1),k(2)) + 2*Pair(k(1),k(3)) + Pair(k(1),k(5)) + 2*Pair(k(2),k
 (5)) - 2*Pair(k(3),k(5))) + 2*ME2*(5*MP2 + 5*Pair(k(3),k(4)) + 6*Pair(k(3),k(5)) + 4*Pair(k(4),k
 (5)))))));

}


double PionPionGamma::Initial(){
  double Aini;
  Aini = - 2 * pow(m_me, 2) * pow(m_s45, -2) * m_pij[0][1] * m_pij[3][4] * pow(m_pij[0][2], -2) 
           + 4 * pow(m_me, 2) * pow(m_s45, -2) * m_pij[0][1] * m_pij[3][4] * pow(m_pij[0][2], -1) * pow(m_pij[1][2], -1)
           - 2 * pow(m_me, 2) * pow(m_s45, -2) * m_pij[0][1] * m_pij[3][4] * pow(m_pij[1][2], -2) 
           + 2 * pow(m_me, 2) * pow(m_s45, -2) * m_pij[0][3] * m_pij[1][4] * pow(m_pij[0][2], -2) 
           + 2 * pow(m_me, 2) * pow(m_s45, -2) * m_pij[0][3] * m_pij[1][4] * pow(m_pij[1][2], -2) 
           - 2 * pow(m_me, 2) * pow(m_s45, -2) * m_pij[0][3] * m_pij[1][3] * pow(m_pij[0][2], -2) 
           - 2 * pow(m_me, 2) * pow(m_s45, -2) * m_pij[0][3] * m_pij[1][3] * pow(m_pij[1][2], -2) 
           - 2 * pow(m_me, 2) * pow(m_s45, -2) * m_pij[1][4] * m_pij[0][4] * pow(m_pij[0][2], -2) 
           - 2 * pow(m_me, 2) * pow(m_s45, -2) * m_pij[1][4] * m_pij[0][4] * pow(m_pij[1][2], -2) 
           + 2 * pow(m_me, 2) * pow(m_s45, -2) * m_pij[0][4] * m_pij[1][3] * pow(m_pij[0][2], -2) 
           + 2 * pow(m_me, 2) * pow(m_s45, -2) * m_pij[0][4] * m_pij[1][3] * pow(m_pij[1][2], -2) 
           - 0.5 * pow(m_me, 2) * pow(m_s45, -1) * pow(m_betapi, 2) * m_pij[0][1] * pow(m_pij[0][2], -2) 
           + pow(m_me, 2) * pow(m_s45, -1) * pow(m_betapi, 2) * m_pij[0][1] * pow(m_pij[0][2], -1) * pow(m_pij[1][2], -1) 
           - 0.5 * pow(m_me, 2) * pow(m_s45, -1) * pow(m_betapi, 2) * m_pij[0][1] * pow(m_pij[1][2], -2) 
           + 0.5 * pow(m_me, 2) * pow(m_s45, -1) * m_pij[0][1] * pow(m_pij[0][2], -2) 
           - pow(m_me, 2) * pow(m_s45, -1) * m_pij[0][1] * pow(m_pij[0][2], -1) * pow(m_pij[1][2], -1) 
           + 0.5 * pow(m_me, 2) * pow(m_s45, -1) * m_pij[0][1] * pow(m_pij[1][2], -2) 
           - 2 * pow(m_me, 4) * pow(m_s45, -2) * m_pij[3][4] * pow(m_pij[0][2], -2) 
           - 2 * pow(m_me, 4) * pow(m_s45, -2) * m_pij[3][4] * pow(m_pij[1][2], -2) 
           - 0.5 * pow(m_me, 4) * pow(m_s45, -1) * pow(m_betapi, 2) * pow(m_pij[0][2], -2) 
           - 0.5 * pow(m_me, 4) * pow(m_s45, -1) * pow(m_betapi, 2) * pow(m_pij[1][2], -2) 
           + 0.5 * pow(m_me, 4) * pow(m_s45, -1) * pow(m_pij[0][2], -2) 
           + 0.5 * pow(m_me, 4) * pow(m_s45, -1) * pow(m_pij[1][2], -2) 
           - 4 * pow(m_s45, -2) * m_pij[0][1] * m_pij[0][3] * m_pij[1][4] * pow(m_pij[0][2], -1) * pow(m_pij[1][2], -1) 
           + 4 * pow(m_s45, -2) * m_pij[0][1] * m_pij[0][3] * m_pij[1][3] * pow(m_pij[0][2], -1) * pow(m_pij[1][2], -1) 
           + 4 * pow(m_s45, -2) * m_pij[0][1] * m_pij[1][4] * m_pij[0][4] * pow(m_pij[0][2], -1) * pow(m_pij[1][2], -1) 
           - 4 * pow(m_s45, -2) * m_pij[0][1] * m_pij[0][4] * m_pij[1][3] * pow(m_pij[0][2], -1) * pow(m_pij[1][2], -1) 
           + 4 * pow(m_s45, -2) * pow(m_pij[0][1], 2) * m_pij[3][4] * pow(m_pij[0][2], -1) * pow(m_pij[1][2], -1) 
           + pow(m_s45, -1) * pow(m_betapi, 2) * pow(m_pij[0][1], 2) * pow(m_pij[0][2], -1) * pow(m_pij[1][2], -1) 
           - pow(m_s45, -1) * pow(m_pij[0][1], 2) * pow(m_pij[0][2], -1) * pow(m_pij[1][2], -1);
  if(IsBad(Aini)){
    msg_Error()<<"NaN in "<<METHOD<<std::endl;
  }
  return Aini;
}

double PionPionGamma::IFI(){
  double Aint;
  Aint =  - 2 * pow(m_me, 2) * pow(m_sp, -1) * pow(m_s45, -1) * m_pij[3][4] * m_pij[0][3] * pow(m_pij[0][2], -1) * pow(m_pij[2][3], -1)
            - 2 * pow(m_me, 2) * pow(m_sp, -1) * pow(m_s45, -1) * m_pij[3][4] * m_pij[1][4] * pow(m_pij[2][4], -1) * pow(m_pij[1][2], -1)
            + 2 * pow(m_me, 2) * pow(m_sp, -1) * pow(m_s45, -1) * m_pij[3][4] * m_pij[0][4] * pow(m_pij[0][2], -1) * pow(m_pij[2][4], -1)
            + 2 * pow(m_me, 2) * pow(m_sp, -1) * pow(m_s45, -1) * m_pij[3][4] * m_pij[1][3] * pow(m_pij[2][3], -1) * pow(m_pij[1][2], -1)
            - 0.5 * pow(m_me, 2) * pow(m_sp, -1) * pow(m_betapi, 2) * m_pij[0][3] * pow(m_pij[0][2], -1) * pow(m_pij[2][3], -1)
            - 0.5 * pow(m_me, 2) * pow(m_sp, -1) * pow(m_betapi, 2) * m_pij[1][4] * pow(m_pij[2][4], -1) * pow(m_pij[1][2], -1)
            + 0.5 * pow(m_me, 2) * pow(m_sp, -1) * pow(m_betapi, 2) * m_pij[0][4] * pow(m_pij[0][2], -1) * pow(m_pij[2][4], -1)
            + 0.5 * pow(m_me, 2) * pow(m_sp, -1) * pow(m_betapi, 2) * m_pij[1][3] * pow(m_pij[2][3], -1) * pow(m_pij[1][2], -1)
            + 0.5 * pow(m_me, 2) * pow(m_sp, -1) * m_pij[0][3] * pow(m_pij[0][2], -1) * pow(m_pij[2][3], -1)
            + 0.5 * pow(m_me, 2) * pow(m_sp, -1) * m_pij[1][4] * pow(m_pij[2][4], -1) * pow(m_pij[1][2], -1)
            - 0.5 * pow(m_me, 2) * pow(m_sp, -1) * m_pij[0][4] * pow(m_pij[0][2], -1) * pow(m_pij[2][4], -1)
            - 0.5 * pow(m_me, 2) * pow(m_sp, -1) * m_pij[1][3] * pow(m_pij[2][3], -1) * pow(m_pij[1][2], -1)
            - 2 * pow(m_sp, -1) * pow(m_s45, -1) * m_pij[0][1] * m_pij[3][4] * m_pij[0][3] * pow(m_pij[0][2], -1) * pow(m_pij[2][3], -1)
            - 2 * pow(m_sp, -1) * pow(m_s45, -1) * m_pij[0][1] * m_pij[3][4] * m_pij[1][4] * pow(m_pij[2][4], -1) * pow(m_pij[1][2], -1)
            + 2 * pow(m_sp, -1) * pow(m_s45, -1) * m_pij[0][1] * m_pij[3][4] * m_pij[0][4] * pow(m_pij[0][2], -1) * pow(m_pij[2][4], -1)
            + 2 * pow(m_sp, -1) * pow(m_s45, -1) * m_pij[0][1] * m_pij[3][4] * m_pij[1][3] * pow(m_pij[2][3], -1) * pow(m_pij[1][2], -1);

    Aint = Aint 
            - 2 * pow(m_sp, -1) * pow(m_s45, -1) * m_pij[0][3] * m_pij[1][4] * m_pij[0][4] * pow(m_pij[0][2], -1) * pow(m_pij[2][3], -1)
            - 2 * pow(m_sp, -1) * pow(m_s45, -1) * m_pij[0][3] * m_pij[1][4] * m_pij[0][4] * pow(m_pij[0][2], -1) * pow(m_pij[2][4], -1)
            - 2 * pow(m_sp, -1) * pow(m_s45, -1) * m_pij[0][3] * m_pij[1][4] * m_pij[1][3] * pow(m_pij[2][3], -1) * pow(m_pij[1][2], -1)
            - 2 * pow(m_sp, -1) * pow(m_s45, -1) * m_pij[0][3] * m_pij[1][4] * m_pij[1][3] * pow(m_pij[2][4], -1) * pow(m_pij[1][2], -1)
            + 2 * pow(m_sp, -1) * pow(m_s45, -1) * m_pij[0][3] * pow(m_pij[1][4], 2) * pow(m_pij[2][4], -1) * pow(m_pij[1][2], -1)
            + 2 * pow(m_sp, -1) * pow(m_s45, -1) * m_pij[0][3] * m_pij[0][4] * m_pij[1][3] * pow(m_pij[0][2], -1) * pow(m_pij[2][3], -1)
            + 2 * pow(m_sp, -1) * pow(m_s45, -1) * m_pij[0][3] * m_pij[0][4] * m_pij[1][3] * pow(m_pij[0][2], -1) * pow(m_pij[2][4], -1)
            + 2 * pow(m_sp, -1) * pow(m_s45, -1) * m_pij[0][3] * pow(m_pij[1][3], 2) * pow(m_pij[2][3], -1) * pow(m_pij[1][2], -1)
            + 2 * pow(m_sp, -1) * pow(m_s45, -1) * pow(m_pij[0][3], 2) * m_pij[1][4] * pow(m_pij[0][2], -1) * pow(m_pij[2][3], -1)
            - 2 * pow(m_sp, -1) * pow(m_s45, -1) * pow(m_pij[0][3], 2) * m_pij[1][3] * pow(m_pij[0][2], -1) * pow(m_pij[2][3], -1)
            + 2 * pow(m_sp, -1) * pow(m_s45, -1) * m_pij[1][4] * m_pij[0][4] * m_pij[1][3] * pow(m_pij[2][3], -1) * pow(m_pij[1][2], -1)
            + 2 * pow(m_sp, -1) * pow(m_s45, -1) * m_pij[1][4] * m_pij[0][4] * m_pij[1][3] * pow(m_pij[2][4], -1) * pow(m_pij[1][2], -1)
            + 2 * pow(m_sp, -1) * pow(m_s45, -1) * m_pij[1][4] * pow(m_pij[0][4], 2) * pow(m_pij[0][2], -1) * pow(m_pij[2][4], -1)
            - 2 * pow(m_sp, -1) * pow(m_s45, -1) * pow(m_pij[1][4], 2) * m_pij[0][4] * pow(m_pij[2][4], -1) * pow(m_pij[1][2], -1)
            - 2 * pow(m_sp, -1) * pow(m_s45, -1) * m_pij[0][4] * pow(m_pij[1][3], 2) * pow(m_pij[2][3], -1) * pow(m_pij[1][2], -1);

    Aint = Aint 
            - 2 * pow(m_sp, -1) * pow(m_s45, -1) * pow(m_pij[0][4], 2) * m_pij[1][3] * pow(m_pij[0][2], -1) * pow(m_pij[2][4], -1)
            - 0.5 * pow(m_sp, -1) * pow(m_betapi, 2) * m_pij[0][1] * m_pij[0][3] * pow(m_pij[0][2], -1) * pow(m_pij[2][3], -1)
            - 0.5 * pow(m_sp, -1) * pow(m_betapi, 2) * m_pij[0][1] * m_pij[1][4] * pow(m_pij[2][4], -1) * pow(m_pij[1][2], -1)
            + 0.5 * pow(m_sp, -1) * pow(m_betapi, 2) * m_pij[0][1] * m_pij[0][4] * pow(m_pij[0][2], -1) * pow(m_pij[2][4], -1)
            + 0.5 * pow(m_sp, -1) * pow(m_betapi, 2) * m_pij[0][1] * m_pij[1][3] * pow(m_pij[2][3], -1) * pow(m_pij[1][2], -1)
            + 0.5 * pow(m_sp, -1) * m_pij[0][1] * m_pij[0][3] * pow(m_pij[0][2], -1) * pow(m_pij[2][3], -1)
            + 0.5 * pow(m_sp, -1) * m_pij[0][1] * m_pij[1][4] * pow(m_pij[2][4], -1) * pow(m_pij[1][2], -1)
            - 0.5 * pow(m_sp, -1) * m_pij[0][1] * m_pij[0][4] * pow(m_pij[0][2], -1) * pow(m_pij[2][4], -1)
            - 0.5 * pow(m_sp, -1) * m_pij[0][1] * m_pij[1][3] * pow(m_pij[2][3], -1) * pow(m_pij[1][2], -1);
  if(IsBad(Aint)){
    msg_Error()<<"NaN in "<<METHOD<<std::endl;
  }
  return Aint;
}


double PionPionGamma::Final(){
  double Afin;
  Afin =  + 0.5 * pow(m_me, 2) * pow(m_sp, -2) * m_s45 * pow(m_betapi, 2) * m_pij[3][4] * pow(m_pij[2][3], -2)
          + pow(m_me, 2) * pow(m_sp, -2) * m_s45 * pow(m_betapi, 2) * m_pij[3][4] * pow(m_pij[2][3], -1) * pow(m_pij[2][4], -1)
          + 0.5 * pow(m_me, 2) * pow(m_sp, -2) * m_s45 * pow(m_betapi, 2) * m_pij[3][4] * pow(m_pij[2][4], -2)
          - 0.5 * pow(m_me, 2) * pow(m_sp, -2) * m_s45 * m_pij[3][4] * pow(m_pij[2][3], -2)
          - pow(m_me, 2) * pow(m_sp, -2) * m_s45 * m_pij[3][4] * pow(m_pij[2][3], -1) * pow(m_pij[2][4], -1)
          - 0.5 * pow(m_me, 2) * pow(m_sp, -2) * m_s45 * m_pij[3][4] * pow(m_pij[2][4], -2)
          - 0.25 * pow(m_me, 2) * pow(m_sp, -2) * pow(m_s45, 2) * pow(m_betapi, 2) * pow(m_pij[2][3], -2)
          - 0.25 * pow(m_me, 2) * pow(m_sp, -2) * pow(m_s45, 2) * pow(m_betapi, 2) * pow(m_pij[2][4], -2)
          + 0.125 * pow(m_me, 2) * pow(m_sp, -2) * pow(m_s45, 2) * pow(m_betapi, 4) * pow(m_pij[2][3], -2)
          + 0.125 * pow(m_me, 2) * pow(m_sp, -2) * pow(m_s45, 2) * pow(m_betapi, 4) * pow(m_pij[2][4], -2)
          + 0.125 * pow(m_me, 2) * pow(m_sp, -2) * pow(m_s45, 2) * pow(m_pij[2][3], -2)
          + 0.125 * pow(m_me, 2) * pow(m_sp, -2) * pow(m_s45, 2) * pow(m_pij[2][4], -2)
          + 4 * pow(m_me, 2) * pow(m_sp, -2) * pow(m_pij[3][4], 2) * pow(m_pij[2][3], -1) * pow(m_pij[2][4], -1)
          + 0.5 * pow(m_sp, -2) * m_s45 * pow(m_betapi, 2) * m_pij[0][1] * m_pij[3][4] * pow(m_pij[2][3], -2)
          + pow(m_sp, -2) * m_s45 * pow(m_betapi, 2) * m_pij[0][1] * m_pij[3][4] * pow(m_pij[2][3], -1) * pow(m_pij[2][4], -1)
          + 0.5 * pow(m_sp, -2) * m_s45 * pow(m_betapi, 2) * m_pij[0][1] * m_pij[3][4] * pow(m_pij[2][4], -2)
          - 0.5 * pow(m_sp, -2) * m_s45 * pow(m_betapi, 2) * m_pij[0][3] * m_pij[1][4] * pow(m_pij[2][3], -2)
          - 0.5 * pow(m_sp, -2) * m_s45 * pow(m_betapi, 2) * m_pij[0][3] * m_pij[1][4] * pow(m_pij[2][4], -2);

  Afin = Afin 
          + 0.5 * pow(m_sp, -2) * m_s45 * pow(m_betapi, 2) * m_pij[0][3] * m_pij[1][3] * pow(m_pij[2][3], -2)
          + 0.5 * pow(m_sp, -2) * m_s45 * pow(m_betapi, 2) * m_pij[0][3] * m_pij[1][3] * pow(m_pij[2][4], -2)
          + 0.5 * pow(m_sp, -2) * m_s45 * pow(m_betapi, 2) * m_pij[1][4] * m_pij[0][4] * pow(m_pij[2][3], -2)
          + 0.5 * pow(m_sp, -2) * m_s45 * pow(m_betapi, 2) * m_pij[1][4] * m_pij[0][4] * pow(m_pij[2][4], -2)
          - 0.5 * pow(m_sp, -2) * m_s45 * pow(m_betapi, 2) * m_pij[0][4] * m_pij[1][3] * pow(m_pij[2][3], -2)
          - 0.5 * pow(m_sp, -2) * m_s45 * pow(m_betapi, 2) * m_pij[0][4] * m_pij[1][3] * pow(m_pij[2][4], -2)
          - 0.5 * pow(m_sp, -2) * m_s45 * m_pij[0][1] * m_pij[3][4] * pow(m_pij[2][3], -2)
          - pow(m_sp, -2) * m_s45 * m_pij[0][1] * m_pij[3][4] * pow(m_pij[2][3], -1) * pow(m_pij[2][4], -1)
          - 0.5 * pow(m_sp, -2) * m_s45 * m_pij[0][1] * m_pij[3][4] * pow(m_pij[2][4], -2)
          + 0.5 * pow(m_sp, -2) * m_s45 * m_pij[0][3] * m_pij[1][4] * pow(m_pij[2][3], -2)
          + 0.5 * pow(m_sp, -2) * m_s45 * m_pij[0][3] * m_pij[1][4] * pow(m_pij[2][4], -2)
          - 0.5 * pow(m_sp, -2) * m_s45 * m_pij[0][3] * m_pij[1][3] * pow(m_pij[2][3], -2)
          - 0.5 * pow(m_sp, -2) * m_s45 * m_pij[0][3] * m_pij[1][3] * pow(m_pij[2][4], -2)
          - 0.5 * pow(m_sp, -2) * m_s45 * m_pij[1][4] * m_pij[0][4] * pow(m_pij[2][3], -2)
          - 0.5 * pow(m_sp, -2) * m_s45 * m_pij[1][4] * m_pij[0][4] * pow(m_pij[2][4], -2)
          + 0.5 * pow(m_sp, -2) * m_s45 * m_pij[0][4] * m_pij[1][3] * pow(m_pij[2][3], -2)
          + 0.5 * pow(m_sp, -2) * m_s45 * m_pij[0][4] * m_pij[1][3] * pow(m_pij[2][4], -2)
          - 0.25 * pow(m_sp, -2) * pow(m_s45, 2) * pow(m_betapi, 2) * m_pij[0][1] * pow(m_pij[2][3], -2)
          - 0.25 * pow(m_sp, -2) * pow(m_s45, 2) * pow(m_betapi, 2) * m_pij[0][1] * pow(m_pij[2][4], -2)
          + 0.125 * pow(m_sp, -2) * pow(m_s45, 2) * pow(m_betapi, 4) * m_pij[0][1] * pow(m_pij[2][3], -2)
          + 0.125 * pow(m_sp, -2) * pow(m_s45, 2) * pow(m_betapi, 4) * m_pij[0][1] * pow(m_pij[2][4], -2);

  Afin = Afin 
          + 0.125 * pow(m_sp, -2) * pow(m_s45, 2) * m_pij[0][1] * pow(m_pij[2][3], -2)
          + 0.125 * pow(m_sp, -2) * pow(m_s45, 2) * m_pij[0][1] * pow(m_pij[2][4], -2)
          + 4 * pow(m_sp, -2) * m_pij[0][1] * pow(m_pij[3][4], 2) * pow(m_pij[2][3], -1) * pow(m_pij[2][4], -1)
          - 4 * pow(m_sp, -2) * m_pij[3][4] * m_pij[0][3] * m_pij[1][4] * pow(m_pij[2][3], -1) * pow(m_pij[2][4], -1)
          + 4 * pow(m_sp, -2) * m_pij[3][4] * m_pij[0][3] * m_pij[1][3] * pow(m_pij[2][3], -1) * pow(m_pij[2][4], -1)
          + 4 * pow(m_sp, -2) * m_pij[3][4] * m_pij[1][4] * m_pij[0][4] * pow(m_pij[2][3], -1) * pow(m_pij[2][4], -1)
          - 4 * pow(m_sp, -2) * m_pij[3][4] * m_pij[0][4] * m_pij[1][3] * pow(m_pij[2][3], -1) * pow(m_pij[2][4], -1);
  if(IsBad(Afin)){
    msg_Error()<<"NaN in "<<METHOD<<std::endl;
  }
  return Afin;
}



DECLARE_TREEME2_GETTER(EXTRAXS::PionPionGamma,"PionPionGamma")
Tree_ME2_Base *ATOOLS::Getter
<PHASIC::Tree_ME2_Base,PHASIC::External_ME_Args,EXTRAXS::PionPionGamma>::
operator()(const External_ME_Args &args) const
{
  // if (pi.m_fi.m_nlotype==nlo_type::loop) {
  // Flavour_Vector fl(pi.ExtractFlavours());
  const Flavour_Vector fl = args.Flavours();
  if(fl.size()!=5) return NULL;
  if (fl[0]!=Flavour(kf_e) && fl[1]!=Flavour(kf_e)) return NULL; 
  int npion=0;
  int ngamma=0;
  for (int i = 2; i < 5; ++i)
  {
    if(fl[i].Kfcode()==kf_pi_plus || fl[i].Kfcode()==-kf_pi_plus) npion++;
    if(fl[i].IsPhoton()) ngamma++;
  }
  if(ngamma!=1 || npion!=2) return NULL;
  return new PionPionGamma(args);
  // if(fl[2])
  // }
  // return NULL;
}
