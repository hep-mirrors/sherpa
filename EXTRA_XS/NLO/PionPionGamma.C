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
  PRINT_INFO("Initialised XS_PionPionGamma");
  Flavour_Vector outflavs = args.m_outflavs;
  // p_bornme = dynamic_cast<ME2_Base*>(PHASIC::Tree_ME2_Base::GetME2(bargs));
  // if (!p_bornme) THROW(fatal_error,"no born me found.");
  const Flavour_Vector fl = args.Flavours();
  double ME = Flavour(kf_e).Mass();
  double MP = Flavour(kf_pi).Mass();
  ME2 = ME*ME;
  MP2 = MP*MP;
  m_oqcd = 1;
  m_oew  = 2;
  m_me  = ME;
  m_mpi = MP;
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
  S34 = (p[3]+p[4]).Abs2();
  T14 = (p[1]-p[4]).Abs2();
  T24 = (p[2]-p[4]).Abs2();
  T=(p[0]-p[2]).Abs2();
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
  return -128*pow(m_alpha,3)*pow(M_PI,3)*((2*
       (2*(ME2*ME2 + MP2*MP2) + 
         ME2*(4*MP2 - S - T - T14 - 
            2*(4*ME2 + 4*MP2 - S - S34 - T - T14 - U)) + (S + T + T14)*U - 
         MP2*(S + T + T14 + 2*U)))/(S*S34*(2*ME2 + 2*MP2 - S - T - T14)) + 
    (-4*(pow(ME2,3) + ME2*ME2*MP2) - 12*(ME2*(MP2*MP2) + pow(MP2,3)) + 
       6*(ME2*ME2)*(4*ME2 + 4*MP2 - S - S34 - T14 - U) + 
       (S34*T + 2*(T*T))*(4*ME2 + 4*MP2 - S - S34 - T - T14 - U) - 
       (S*(S34 + T) + S34*T14 - 
          (S + T + T14)*(4*ME2 + 4*MP2 - S - S34 - T - T14 - U))*U - 
       (S + T14)*(U*U) + 2*(MP2*MP2)*
        (4*T - T14 + 2*(4*ME2 + 4*MP2 - S - T - T14 - U) + U) - 
       T*(-pow(4*ME2 + 4*MP2 - S - S34 - T - T14 - U,2) + 
          S*(-4*ME2 - 4*MP2 + S + 2*S34 + T + T14 + U)) + 
       MP2*((2*S - S34)*T - 2*(T*T) - 
          (5*T - 3*T14)*(4*ME2 + 4*MP2 - S - S34 - T - T14 - U) - 
          pow(4*ME2 + 4*MP2 - S - S34 - T - T14 - U,2) + 
          (2*S + S34 - T + 3*T14 - 
             2*(4*ME2 + 4*MP2 - S - S34 - T - T14 - U))*U + U*U + 
          S34*(-4*ME2 - 4*MP2 + S + S34 + T + 2*T14 + U)) + 
       ME2*(-2*(T*T) + S34*(2*S - T + T14) - 
          (2*S + S34 + 7*T + T14)*(4*ME2 + 4*MP2 - S - S34 - T - T14 - U) - 
          pow(4*ME2 + 4*MP2 - S - S34 - T - T14 - U,2) + 
          (2*S + S34 - T + T14 - 
             2*(4*ME2 + 4*MP2 - S - S34 - T - T14 - U))*U + U*U - 
          2*MP2*(2*S - 5*T + 3*T14 - 
             3*(4*ME2 + 4*MP2 - S - S34 - T - T14 - U) + U)))/
     (S*S34*(2*ME2 + 2*MP2 - S - T - T14)*(2*ME2 + 2*MP2 - S34 - T - U)) - 
    (8*(ME2*(MP2*MP2) + pow(MP2,3)) - pow(T,3) - 
       T*(S34*S34 + 4*(ME2 + MP2)*T14 - T14*T14) - 
       MP2*(8*(ME2 + MP2)*S34 + S34*S34 + 4*(2*(ME2 + MP2) - S34)*T - 
          T*T - 2*(2*(ME2 + MP2) + S34 + 3*T)*T14 - T14*T14 + 
          2*(2*(ME2 + MP2) + S34)*U) + 
       2*(ME2*ME2*(S34 + T - T14) + MP2*MP2*(S34 - T - 5*T14 + 2*U)) + 
       ME2*(-4*MP2*(-2*(ME2 + MP2 + S34) + 3*T14 - U) - 
          2*(S34*S34 + (2*(ME2 + MP2) - S34)*T + 
             2*(ME2 + MP2)*(S34 - T14) + S34*U)) + 
       2*((2*(ME2 + MP2) - S34)*(T*T) + (ME2 + MP2)*(S34*S34 + S34*(T + U)))\
)/(S34*S34*pow(-2*ME2 - 2*MP2 + S34 + T + T14,2)) + 
    (-8*(ME2*(MP2*MP2) + pow(MP2,3)) + 
       (S34*S34 + 4*(ME2 + MP2)*(4*ME2 + 4*MP2 - S - S34 - T - T14 - U) - 
          pow(4*ME2 + 4*MP2 - S - S34 - T - T14 - U,2))*U + pow(U,3) + 
       2*ME2*(S34*S34 - 2*(ME2 + MP2)*
           (4*ME2 + 4*MP2 - S - S34 - T - T14 - 2*U) + 
          S34*(2*(ME2 + MP2) - T + U)) - 
       2*((ME2 + MP2)*(S34*S34 - S34*(T - 3*U)) + 
          (2*(ME2 + MP2) - S34)*(U*U) + 
          ME2*ME2*(-4*ME2 - 4*MP2 + S + 2*S34 + T + T14 + 2*U) + 
          MP2*MP2*(-4*ME2 - 4*MP2 + S + 2*S34 - T - 3*T14 + 4*U)) + 
       MP2*(S34*S34 - 2*(2*(ME2 + MP2) + S34)*T - 
          pow(4*ME2 + 4*MP2 - S - S34 - T - T14 - U,2) - U*U - 
          4*(ME2 + MP2)*(-4*ME2 - 4*MP2 + S + S34 + T + 3*T14 + U) - 
          2*ME2*(-4*ME2 - 4*MP2 + 2*S + S34 - 2*T14 + U + 5*(S34 + U)) + 
          2*((8*(ME2 + MP2) - 3*(4*ME2 + 4*MP2 - S - S34 - T - T14 - U))*
              U + S34*(-4*ME2 - 4*MP2 + 4*(ME2 + MP2) + S + S34 + T + T14 + 
                U))))/(S34*S34*pow(2*ME2 + 2*MP2 - S - T - T14,2)) - 
    (2*((2*(ME2*ME2 + 2*ME2*MP2 + MP2*MP2 + T*U - ME2*(T + U) - 
              MP2*(S + T + U)))/(S*S) + 
         (2*pow(ME2,3) - 18*pow(MP2,3) - S34*(S + 2*T)*U + T14*(U*U) + 
            MP2*MP2*(13*T14 + 12*(4*ME2 + 4*MP2 - S - T - T14 - U) + U) - 
            T*((4*ME2 + 4*MP2 - S - S34 - T - U)*U + U*U) - 
            ME2*ME2*(4*ME2 + 4*MP2 - S - S34 + 2*T - T14 - U + 
               2*(S34 + U)) + MP2*
             (-6*(ME2*ME2) - (S + 4*S34 + T)*T14 - 
               2*(S34*S34 + T14*T14 + 
                  pow(4*ME2 + 4*MP2 - S - S34 - T - T14 - U,2)) - 
               4*(S34 + T14)*(4*ME2 + 4*MP2 - S - S34 - T - T14 - U) + 
               U*(-4*ME2 - 4*MP2 + 4*S + S34 + 7*T - T14 + U)) + 
            ME2*(-26*(MP2*MP2) + S34*(S + 2*T) + 
               T*(4*ME2 + 4*MP2 - S - S34 - T - U) + 
               (4*ME2 + 4*MP2 - S + S34 + 3*T - 2*T14 - U)*U - 
               MP2*(2*S - 8*S34 + 5*T - 11*T14 - 
                  9*(4*ME2 + 4*MP2 - S - S34 - T - T14 - U) + 3*U)))/
          (S*S34*(2*ME2 + 2*MP2 - S - T - T14))))/
     (-2*ME2 - 2*MP2 + S + T + U) - 
    (-((32*pow(MP2,3) - S*(S34*S34) - pow(S34,3) + 
            S34*(2*S*T + T*T) - (S34*S34 - S34*T)*T14 - 
            (S34*S34 + S*T + 2*(T*T))*
             (4*ME2 + 4*MP2 - S - S34 - T - T14 - U) - 
            T*pow(4*ME2 + 4*MP2 - S - S34 - T - T14 - U,2) - 
            (2*(S34*S34) - 2*(S + S34)*T - T*T + 
               (S + 2*S34 - 2*T)*T14 + T14*T14 + 
               (S34 - 2*T)*(4*ME2 + 4*MP2 - S - S34 - T - T14 - U))*U - 
            (S34 - T + 2*T14)*(U*U) + 
            2*(ME2*ME2)*(-4*ME2 - 4*MP2 + S + S34 + 2*T + 2*U) + 
            MP2*(16*(ME2*ME2) + 10*(S34*S34) - (5*S + S34)*T - T*T + 
               S*(4*S34 - T14) + (7*S34 - 4*T)*T14 + T14*T14 + 
               (3*S + 11*S34 + 6*T + 4*T14)*
                (4*ME2 + 4*MP2 - S - S34 - T - T14 - U) + 
               5*pow(4*ME2 + 4*MP2 - S - S34 - T - T14 - U,2) - 
               (S - 11*S34 + 8*T - 6*T14 - 
                  4*(4*ME2 + 4*MP2 - S - S34 - T - T14 - U))*U + 3*(U*U)) \
+ MP2*MP2*(48*ME2 - 2*(16*S34 - T + 5*T14 + 
                  13*(4*ME2 + 4*MP2 - S - S34 - T - T14 - U) + 7*U)) + 
            ME2*(2*(16*ME2*MP2 + 8*(ME2*ME2 + MP2*MP2) - S*S34 + 
                  3*(S34*S34) - 2*(ME2 + MP2)*(S + S34) - 
                  (8*(ME2 + MP2) - S + S34)*T + 
                  2*(T*T + T*(4*ME2 + 4*MP2 - S - S34 - T - T14 - U)) - 
                  (S - S34 + 2*(4*ME2 + 4*MP2 - S - S34 - T14 - U))*U) + 
               4*MP2*(S - 9*S34 + 2*T - 
                  6*(4*ME2 + 4*MP2 - S - S34 - T - T14 - U) - 2*(T14 + U))\
))/(S34*S34*(2*ME2 + 2*MP2 - S - T - T14))) + 
       (2*(4*ME2*MP2 + 2*(ME2*ME2 + MP2*MP2) + 
             T*(4*ME2 + 4*MP2 - S34 - T - T14) - 
             MP2*(4*ME2 + 4*MP2 - S34 + T - T14) - 
             ME2*(4*ME2 + 4*MP2 - S34 - T + T14)) + 
          (-4*(pow(ME2,3) + ME2*ME2*MP2) - 
             12*(ME2*(MP2*MP2) + pow(MP2,3)) - 
             (T*T + T*(S34 - T14))*(4*ME2 + 4*MP2 - S34 - T - T14 - U) - 
             (S*(S34 + T - T14) - (S34 + T)*T14 - T14*T14)*U + 
             2*T14*(U*U) + 6*(ME2*ME2)*(T14 + U) + 
             2*(MP2*MP2)*(-4*ME2 - 4*MP2 + S + S34 + 2*T + T14 + 
                2*(S34 + T14) + 5*U) + 
             MP2*((2*S + S34)*T + T*T - (S34 + 2*T)*T14 - T14*T14 + 
                (S34 + 3*(T + T14))*
                 (4*ME2 + 4*MP2 - S - S34 - T - T14 - U) + 
                (2*S - S34 - T - 5*T14)*U - 2*(U*U)) - 
             ME2*(-(S34*T) - T*T - 2*S*(S34 + T) + 
                (S34 + 2*(S + T))*T14 + T14*T14 + 
                2*MP2*(2*S + T - 3*T14 + 
                   3*(4*ME2 + 4*MP2 - S - S34 - T - T14 - U) - 5*U) - 
                (S34 + T - T14)*
                 (4*ME2 + 4*MP2 - S - S34 - T - T14 - U) + 
                (S34 + T + 7*T14)*U + 2*(U*U)))/
           (2*ME2 + 2*MP2 - S34 - T - U) - 
          (2*(2*pow(ME2,3) - 18*pow(MP2,3) - S*S34*T + 
               MP2*MP2*(T + 12*(S34 + T14) + 
                  13*(4*ME2 + 4*MP2 - S - S34 - T - T14 - U)) + 
               T*T*(4*ME2 + 4*MP2 - S - S34 - T - T14 - U) - 
               (T*T + T*(4*ME2 + 4*MP2 - S + S34 - T - U))*U - 
               ME2*ME2*(2*(S34 + T) + T14 + 3*U) + 
               ME2*(-26*(MP2*MP2) + S34*(S + 2*T) + 
                  (4*ME2 + 4*MP2 - S + S34 + 3*T - U)*U + 
                  T*(-4*ME2 - 4*MP2 + S + S34 + T + 2*T14 + U)) + 
               MP2*(-6*(ME2*ME2) + 3*S*T - (4*S34 + T)*T14 - 
                  2*(S34*S34 + T14*T14 + 
                     pow(4*ME2 + 4*MP2 - S - S34 - T - T14 - U,2)) - 
                  (S + 2*T + 4*(S34 + T14))*
                   (4*ME2 + 4*MP2 - S - S34 - T - T14 - U) + 
                  U*(-4*ME2 - 4*MP2 + S + S34 + 7*T + T14 + U) - 
                  ME2*(2*S - 8*S34 + 3*T - 9*T14 - 
                     11*(4*ME2 + 4*MP2 - S - S34 - T - T14 - U) + 5*U))))/
           (-2*ME2 - 2*MP2 + S + T + U))/(S*S34))/
     (-2*ME2 - 2*MP2 + S34 + T + T14) - 
    (-2*((2*(-3*pow(MP2,3) + 
               MP2*(ME2*ME2 - ME2*(4*ME2 + 4*MP2 - S - S34 - T - U) + 
                  T14*(4*ME2 + 4*MP2 - S - S34 - T - T14 - U)) + 
               MP2*MP2*(-2*ME2 + S34 + T + U)))/
           pow(2*ME2 + 2*MP2 - S34 - T - U,2) + 
          (4*(pow(ME2,3) + ME2*(MP2*MP2)) + S*S34*T - 
             T*T*(4*ME2 + 4*MP2 - S - S34 - T - T14 - U) + 
             (S*S34 - T*(4*ME2 + 4*MP2 - S - S34 - T - U))*U - T14*(U*U) - 
             (2*(MP2*MP2) - MP2*(4*(ME2 + MP2) - 3*S - S34))*(T + U) + 
             ME2*(-2*S*S34 + T*T + 
                (4*ME2 + 4*MP2 - S - S34 - T + 2*T14 - U)*U + U*U + 
                T*(T14 + 3*(4*ME2 + 4*MP2 - S - S34 - T - T14 - U) + 
                   2*U) + 2*MP2*
                 (-4*ME2 - 4*MP2 + 3*S + S34 + T + U - 3*(T + U))) + 
             ME2*ME2*(8*MP2 - 2*
                 (4*ME2 + 4*MP2 - S - S34 - T - U + 2*(T + U))))/
           ((2*ME2 + 2*MP2 - S34 - T - U)*(-2*ME2 - 2*MP2 + S + T + U))) + 
       2*(2*ME2 - (2*(ME2*ME2 + MP2*MP2) - MP2*(4*ME2 + 4*MP2 - S - S34) + 
             ME2*(4*MP2 - 4*(ME2 + MP2) - S + S34) + 
             T*(4*ME2 + 4*MP2 - S34 - T - T14 - U) + (S + T14)*U)/
           (2*ME2 + 2*MP2 - S34 - T - U) - 
          (2*(2*pow(ME2,3) - ME2*ME2*MP2 - 12*ME2*(MP2*MP2) - 
               9*pow(MP2,3) - 
               MP2*(S34*S34 + T14*T14 + 
                  2*(S34*T14 + 
                     (S34 + T14)*(4*ME2 + 4*MP2 - S - S34 - T - T14 - U)) \
+ pow(4*ME2 + 4*MP2 - S - S34 - T - T14 - U,2) - 3*T*U) + 
               (4*ME2 + 4*MP2 - S - T - U)*(6*(MP2*MP2) - T*U) - 
               ME2*ME2*(4*ME2 + 4*MP2 - S - T - U + 2*(T + U)) + 
               ME2*(T*(4*ME2 + 4*MP2 - S - T - U) + 
                  (4*ME2 + 4*MP2 - S + T - U)*U - 
                  MP2*(-4*(4*ME2 + 4*MP2 - S - T - U) + 3*(T + U)))))/
           pow(-2*ME2 - 2*MP2 + S + T + U,2)))/(S*S));

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
  if(args.m_source!="Internal") return NULL;
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
