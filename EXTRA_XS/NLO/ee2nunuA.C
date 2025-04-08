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
  class ee2nunuA :  public ME2_Base {
  public:
    ee2nunuA(const External_ME_Args& args);
    ~ee2nunuA(){};

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
    double SW2, SW, CW, CW2, MZ2;
    Vec4D_Vector m_pmom;
  };
}

using namespace EXTRAXS;


ee2nunuA::ee2nunuA(const External_ME_Args& args) : ME2_Base(args)
{
  PRINT_INFO("initialised XS_ee2nunuA");
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
  SW2 = MODEL::s_model->ComplexConstant("csin2_thetaW").real();
  CW2 = 1.-SW*SW;
  MZ2 = Flavour(23).Mass()*Flavour(23).Mass();
  CW = sqrt(CW2);
  SW = sqrt(SW2);
}


double ee2nunuA::operator()
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


  res = Full();
  // if(IsBad(res)){
  //    for (int i = 1; i < p.size(); ++i)
  //   {
  //       for (int j = i; j < p.size(); ++j)
  //       {
  //           msg_Out()<<"pair(k("<<i<<"),("<<j<<")) = "<<Pair(k(i),k(j))<<std::endl;
  //       }
  //   }
  //   PRINT_VAR(ME2);
  //   PRINT_VAR(MP2);
  // }
  // PRINT_VAR(res);
  // PRINT_VAR(Initial());
  // PRINT_VAR(Final());
  // PRINT_VAR(IFI());
  // PRINT_VAR(res);
  // PRINT_VAR(p[4]);
  // return ee2uug_ee;
  return res;
}


double ee2nunuA::Full(){
  return pow(m_alpha*M_PI,3)*pow(Den(S34,MZ2),2)*
    (((16*(ME2*ME2)*(2*ME2 + 3*S34 + T + T14 - 2*T24) - 
         (4*(4*pow(ME2,3) + 
              2*(ME2*ME2*(10*ME2 - 2*S + 2*T - 3*T24 - U) + 
                 2*ME2*(2*(S + T) + T14)*U - (S + T + T14)*(U*U)) - 
              2*(ME2*ME2)*(4*ME2 - S34 + 2*(T + U)) + 
              ME2*(-2*S34*(2*ME2 - U) + 2*(3*S34 + 2*T14 + 4*T24)*U + 
                 8*(U*U) + 4*ME2*(T + U) - 
                 4*ME2*(2*(S34 + T) + T14 + 7*U))))/SW2 + 
         (4*pow(ME2,3) - ME2*ME2*
             (16*ME2 - 3*S - S34 - 7*T - 7*T14 - T24 - U) + 
            2*(2*ME2*(2*(S + T) + T14)*U - (S + T + T14)*(U*U)) + 
            ME2*((S34 + T24)*U + (4*T14 + 7*(S34 + T24))*U + 8*(U*U) + 
               4*ME2*(T + 3*U) - 4*ME2*(2*T + T14 + 9*U)))/(SW2*SW2) + 
         8*((-(S34*S34) + 2*ME2*(2*S34 - T))*T24 - S34*(T24*T24) + 
            (16*(ME2*ME2) - T24*T24 - 2*(ME2*(2*S34 + T14) + S34*T24))*U - 
            (8*ME2 - S34)*(U*U) + pow(U,3) + 
            ME2*(8*(ME2*ME2) - 6*ME2*S + S*S - T*T + T14*T14 - 
               2*(S34*(7*ME2 + T) + (4*ME2 - S)*T14 + (2*S + 3*T + T14)*U))))*pow(Den(3*ME2 - S - T - T14,ME2),2) - 
      (8*(-2*(S*S*S34 + S*(S34*S34)) + (S34*T + 2*(T*T))*T24 + 
            S*(T*T - (S34 - T)*T14 - (S34 - 2*T)*T24) + 
            T*(S*S - S*S34 + T24*T24) + 
            (S*S + T*T + (2*S + S34)*T14 + T14*T14 - 
               S*(S34 - 2*T - T24))*U + (S + T + 2*T14)*(U*U)) - 
         (8*(ME2*ME2*(8*ME2 - S - 5*S34 - 4*T14) + 
              (S*S34 - (8*ME2 - S - S34)*T + T*T)*T14 + 
              T*(4*(4*(ME2*ME2) - ME2*(S + S34)) + T14*T14) + 
              (T + T14)*(U*U) + 4*ME2*(-(S*S34) - T*T + (S - T)*U) + 
              ME2*(S34*S34 + 4*ME2*(2*S34 - T) - (S34 - 2*T)*T14 - 
                 (4*ME2 - S34 + T - T14)*T24 + S34*(T - U) - 
                 (16*ME2 - 3*T - 5*T14 - 4*(S34 + T24))*U + 2*(U*U))))/SW2+ 16*(ME2*(8*ME2*S34 + S34*S34 - 2*(5*ME2 - S34)*T + T*T - 
               (2*ME2 - 3*T)*T14 - (2*ME2 - 2*T - T14)*T24 - 
               (10*ME2 - 2*(S34 + T14) - 3*(T + T24))*U + U*U) + 
            ME2*ME2*(S - 3*S34 + 2*(T + U))) - 
         (16*pow(ME2,3) - ME2*ME2*
             (48*ME2 - 10*S - 18*S34 - 10*T14 - 2*T24 + 2*(T14 + T24)) - 
            ME2*(8*S*S34 + 4*(S34*S34) - 2*(T*T) + 2*S*(S34 + T) + 
               T*(S34 - T14) - T*(4*S - S34 + T14) - 
               2*((S - S34 + 5*T)*T24 + T24*T24) + 
               (2*S - 3*S34 - T24)*U - 
               (10*S - 3*S34 + 4*T + 2*T14 + 5*T24)*U - 4*(U*U)) - 
            2*((S*S34 - (8*ME2 - S - S34)*T + T*T)*T14 + 
               T*(4*(4*(ME2*ME2) - ME2*(S + S34)) + T14*T14) + 
               (T + T14)*(U*U) + 4*ME2*(-(S*S34) - T*T + (S - T)*U)))/
          (SW2*SW2))*Den(3*ME2 - S - T - T14,ME2)*
       Den(3*ME2 - S - T24 - U,ME2) + 
      (ME2*ME2*(-((3*(S34 + T) + T14 - 4*U)/(SW2*SW2)) + 
            16*(2*ME2 + 3*S34 - 2*T14 + T24 + U)) + 
         8*ME2*(8*(ME2*ME2) + S*S - 4*S*T + T24*T24 - U*U - 
            2*(ME2*(3*S + 7*S34) + (4*ME2 - S + T)*T24 + (S34 + 3*T)*U)) - 
         8*(-4*(4*(ME2*ME2) - ME2*S34)*T + (8*ME2 - S34)*(T*T) - 
            pow(T,3) + (S34*S34 - 2*S34*(2*ME2 - T))*T14 + 
            (S34 + T)*(T14*T14) + 2*ME2*(T*T24 + T14*U)) + 
         (4*(4*pow(ME2,3) + 
              2*((S34 + T)*(T14*T14) + 
                 T14*(S34*S34 + T*T - 
                    2*((2*ME2 - S34)*T + ME2*(2*S34 - U)))) + 
              ME2*(24*ME2*S34 - 2*(S34*S34) + 2*(4*ME2 - 3*S34)*T + 
                 2*(2*ME2 + S34)*T + 4*(ME2 - S34 - T)*T14 - 
                 2*(T*T + T14*T14) - 4*ME2*(S34 - U) - 4*(ME2 + T14)*U) - 
              2*(ME2*ME2)*(6*ME2 + 3*S34 - T - 3*T14 - 2*U + 2*(T + U))))/
          SW2 - (4*pow(ME2,3) + 
            2*((S34 + T)*(T14*T14) + 
               T14*(S34*S34 + T*T - 2*((2*ME2 - S34)*T + ME2*(2*S34 - U)))) - ME2*(-2*(S34*S34) - T*T + T*(4*ME2 - 3*S34 + T14) + 
               4*ME2*(2*S34 - U)) - 
            ME2*ME2*(12*ME2 + S34 + 5*T - 5*T14 - 4*U) + 
            ME2*(16*ME2*S34 + (16*ME2 - 7*S34)*T - 3*(T*T) + 
               (4*ME2 - 4*S34 - 3*T)*T14 - 2*(T14*T14) - 
               4*(S34*S34 + (ME2 + T14)*U)))/(SW2*SW2))*pow(Den(3*ME2 - S - T24 - U,ME2),2)))/(CW2*CW2);
}


DECLARE_TREEME2_GETTER(EXTRAXS::ee2nunuA,"ee2nunuA")
Tree_ME2_Base *ATOOLS::Getter
<PHASIC::Tree_ME2_Base,PHASIC::External_ME_Args,EXTRAXS::ee2nunuA>::
operator()(const External_ME_Args &args) const
{
  const Flavour_Vector fl = args.Flavours();
  PRINT_INFO(fl);
  if(fl.size()!=5) return NULL;
  return new ee2nunuA(args);
}
