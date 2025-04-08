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
    double m_me, m_mpi, m_pij[5][5], m_s45;
    double m_betapi, m_s, m_sp;
    double ME2, MP2, m_alpha;
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
  m_sp = (p[0]+p[1]).Abs2();
  m_s45 = (p[3]+p[4]).Abs2();
  m_betapi = sqrt(1-m_mpi*m_mpi/m_sp);
  // for (int i = 0; i < p.size(); ++i)
  // {
  //   for (int j = i; j < p.size(); ++j)
  //   {
  //     m_pij[i][j]= p[i]*p[j];
  //     // if(IsZero(m_pij[i][j])){
  //     //   msg_Error()<<"Divide by Zero in "<<METHOD<<std::endl
  //     //              <<"p_"<<i<<j<<" = "<< m_pij[i][j] <<std::endl;
  //     // } 
  //   }
  // }
  m_alpha = (*aqed)(m_sp);
  // PRINT_VAR(Full());
  double res;// = Initial()+Final()+IFI();
  res = Full();
  // PRINT_VAR(res);
  // PRINT_VAR(p[4]);
  // return 1;
  return Full();
}


double PionPionGamma::Full(){
  return -pow(M_PI*m_alpha,3)*
     (-512*ME2*
       pow(1/(2.*(ME2 + Pair(m_pmom[0],m_pmom[1]))) + 
         (MP2*(-1/2.*1/Pair(m_pmom[0],m_pmom[4]) - 1/(2.*Pair(m_pmom[1],m_pmom[4]))))/
          (2.*(MP2 + Pair(m_pmom[2],m_pmom[3]))),2) + 
      ((1/(2.*(ME2 + Pair(m_pmom[0],m_pmom[1]))) + 
           (MP2*(-1/2.*1/Pair(m_pmom[0],m_pmom[4]) - 1/(2.*Pair(m_pmom[1],m_pmom[4]))))/
            (2.*(MP2 + Pair(m_pmom[2],m_pmom[3]))))*
         ((-1/2.*1/Pair(m_pmom[0],m_pmom[4]) - 1/(2.*Pair(m_pmom[1],m_pmom[4])))*
            (512*(-(Pair(m_pmom[0],m_pmom[4])*Pair(m_pmom[1],m_pmom[2])) + 
                 Pair(m_pmom[0],m_pmom[2])*Pair(m_pmom[1],m_pmom[4])) - 
              1024*ME2*(4*MP2 - 6*Pair(m_pmom[0],m_pmom[2]) - 4*Pair(m_pmom[0],m_pmom[3]) + 
                 5*Pair(m_pmom[2],m_pmom[3]))) + 
           ME2*(-1024*(-2/Pair(m_pmom[0],m_pmom[4]) - 5/(2.*Pair(m_pmom[1],m_pmom[4])))*
               Pair(m_pmom[1],m_pmom[4]) - 
              (512*(Pair(m_pmom[0],m_pmom[4]) - 2*Pair(m_pmom[2],m_pmom[4])))/Pair(m_pmom[0],m_pmom[4]))))/
       (2.*(MP2 + Pair(m_pmom[2],m_pmom[3]))) - 
      256*(MP2*MP2)*(ME2 + Pair(m_pmom[0],m_pmom[1]))*
       pow((-1/2.*1/Pair(m_pmom[0],m_pmom[4]) + 1/(2.*Pair(m_pmom[1],m_pmom[4])))/
          (2.*(MP2 + Pair(m_pmom[2],m_pmom[3]))) + 
         1/(2.*(ME2 + Pair(m_pmom[0],m_pmom[1]))*Pair(m_pmom[2],m_pmom[4])),2) + 
      (-512*(MP2*Pair(m_pmom[0],m_pmom[1]) - 2*Pair(m_pmom[0],m_pmom[2])*Pair(m_pmom[1],m_pmom[2]) + 
            ME2*(Pair(m_pmom[0],m_pmom[2]) + Pair(m_pmom[1],m_pmom[2])))*
          (1/(2.*(ME2 + Pair(m_pmom[0],m_pmom[1]))) + 
            (MP2*(-1/2.*1/Pair(m_pmom[0],m_pmom[4]) - 1/(2.*Pair(m_pmom[1],m_pmom[4]))) + 
               (-2/Pair(m_pmom[0],m_pmom[4]) - 5/(2.*Pair(m_pmom[1],m_pmom[4])))*
                Pair(m_pmom[1],m_pmom[4]))/(2.*(MP2 + Pair(m_pmom[2],m_pmom[3])))) + 
         512*(ME2*(1/(2.*(ME2 + Pair(m_pmom[0],m_pmom[1]))) + 
               (MP2*(-1/2.*1/Pair(m_pmom[0],m_pmom[4]) - 1/(2.*Pair(m_pmom[1],m_pmom[4]))))/
                (2.*(MP2 + Pair(m_pmom[2],m_pmom[3]))))*
             (Pair(m_pmom[2],m_pmom[3]) + Pair(m_pmom[2],m_pmom[4])) + 
            ((-1/2.*1/Pair(m_pmom[0],m_pmom[4]) + 1/(2.*Pair(m_pmom[1],m_pmom[4])))/
                (2.*(MP2 + Pair(m_pmom[2],m_pmom[3]))) + 
               1/(2.*(ME2 + Pair(m_pmom[0],m_pmom[1]))*Pair(m_pmom[2],m_pmom[4])))*
             (Pair(m_pmom[0],m_pmom[3])*
                ((MP2 - Pair(m_pmom[0],m_pmom[2]))*Pair(m_pmom[1],m_pmom[2]) - 
                  pow(Pair(m_pmom[1],m_pmom[2]),2)) + 
               (-pow(Pair(m_pmom[0],m_pmom[2]),2) + 
                  Pair(m_pmom[0],m_pmom[2])*(MP2 - Pair(m_pmom[1],m_pmom[2])))*Pair(m_pmom[1],m_pmom[3]) 
+ Pair(m_pmom[0],m_pmom[1])*(-MP2 + Pair(m_pmom[0],m_pmom[2]) + Pair(m_pmom[1],m_pmom[2]))*
                Pair(m_pmom[2],m_pmom[3]) + 
               ME2*Pair(m_pmom[2],m_pmom[3])*(Pair(m_pmom[2],m_pmom[3]) + Pair(m_pmom[2],m_pmom[4])))) + 
         (512*(MP2*(-1/2.*1/Pair(m_pmom[0],m_pmom[4]) - 1/(2.*Pair(m_pmom[1],m_pmom[4])))*
                (-(Pair(m_pmom[0],m_pmom[4])*Pair(m_pmom[1],m_pmom[3])) + 
                  Pair(m_pmom[0],m_pmom[3])*Pair(m_pmom[1],m_pmom[4])) + 
               ME2*(-2/Pair(m_pmom[0],m_pmom[4]) - 5/(2.*Pair(m_pmom[1],m_pmom[4])))*
                Pair(m_pmom[1],m_pmom[4])*(Pair(m_pmom[2],m_pmom[3]) + Pair(m_pmom[2],m_pmom[4]))) - 
            512*(-1/2.*1/Pair(m_pmom[0],m_pmom[4]) - 1/(2.*Pair(m_pmom[1],m_pmom[4])))*
             (4*(MP2*MP2*Pair(m_pmom[0],m_pmom[1]) + ME2*MP2*Pair(m_pmom[0],m_pmom[2])) - 
               6*(MP2*Pair(m_pmom[0],m_pmom[1])*Pair(m_pmom[0],m_pmom[2]) + 
                  ME2*pow(Pair(m_pmom[0],m_pmom[2]),2)) + 
               (ME2*(4*MP2 - 6*Pair(m_pmom[0],m_pmom[2])) - 
                  8*MP2*Pair(m_pmom[0],m_pmom[2]) + 12*pow(Pair(m_pmom[0],m_pmom[2]),2) + 
                  (-4*ME2 - MP2 + 9*Pair(m_pmom[0],m_pmom[2]))*Pair(m_pmom[0],m_pmom[3]))*
                Pair(m_pmom[1],m_pmom[2]) + 
               Pair(m_pmom[0],m_pmom[3])*
                (-4*(MP2*Pair(m_pmom[0],m_pmom[1]) + ME2*Pair(m_pmom[0],m_pmom[2])) + 
                  pow(Pair(m_pmom[1],m_pmom[2]),2)) + 
               (-pow(Pair(m_pmom[0],m_pmom[2]),2) + 
                  Pair(m_pmom[0],m_pmom[2])*(MP2 - Pair(m_pmom[1],m_pmom[2])))*Pair(m_pmom[1],m_pmom[3]) 
+ (5*MP2*Pair(m_pmom[0],m_pmom[1]) + (6*ME2 + Pair(m_pmom[0],m_pmom[1]))*Pair(m_pmom[0],m_pmom[2]) + 
                  (4*ME2 - Pair(m_pmom[0],m_pmom[1]) - 10*Pair(m_pmom[0],m_pmom[2]))*
                   Pair(m_pmom[1],m_pmom[2]))*Pair(m_pmom[2],m_pmom[3]) + 
               ME2*(-((4*MP2 - 6*Pair(m_pmom[0],m_pmom[2]) - 4*Pair(m_pmom[0],m_pmom[3]))*
                     Pair(m_pmom[2],m_pmom[3])) - 5*pow(Pair(m_pmom[2],m_pmom[3]),2) - 
                  (4*MP2 - 6*Pair(m_pmom[0],m_pmom[2]) - 4*Pair(m_pmom[0],m_pmom[3]) + 
                     5*Pair(m_pmom[2],m_pmom[3]))*Pair(m_pmom[2],m_pmom[4]))))/
          (2.*(MP2 + Pair(m_pmom[2],m_pmom[3]))))*
       ((-1/2.*1/Pair(m_pmom[0],m_pmom[4]) + 1/(2.*Pair(m_pmom[1],m_pmom[4])))/
          (2.*(MP2 + Pair(m_pmom[2],m_pmom[3]))) - 
         1/(2.*(ME2 + Pair(m_pmom[0],m_pmom[1]))*Pair(m_pmom[3],m_pmom[4]))) + 
      pow((-1/2.*1/Pair(m_pmom[0],m_pmom[4]) + 1/(2.*Pair(m_pmom[1],m_pmom[4])))/
          (2.*(MP2 + Pair(m_pmom[2],m_pmom[3]))) - 
         1/(2.*(ME2 + Pair(m_pmom[0],m_pmom[1]))*Pair(m_pmom[3],m_pmom[4])),2)*
       (Pair(m_pmom[0],m_pmom[2])*Pair(m_pmom[1],m_pmom[2])*
          (512*MP2 - 1024*(MP2 - Pair(m_pmom[0],m_pmom[3]) - Pair(m_pmom[1],m_pmom[3]) + 
               Pair(m_pmom[2],m_pmom[3]))) + 
         (ME2 + Pair(m_pmom[0],m_pmom[1]))*(-256*(MP2*MP2) - 512*MP2*Pair(m_pmom[3],m_pmom[4]))) 
+ 512*(MP2*Pair(m_pmom[0],m_pmom[3])*Pair(m_pmom[1],m_pmom[3])*
          pow((-1/2.*1/Pair(m_pmom[0],m_pmom[4]) + 1/(2.*Pair(m_pmom[1],m_pmom[4])))/
             (2.*(MP2 + Pair(m_pmom[2],m_pmom[3]))) + 
            1/(2.*(ME2 + Pair(m_pmom[0],m_pmom[1]))*Pair(m_pmom[2],m_pmom[4])),2) - 
         (((-1/2.*1/Pair(m_pmom[0],m_pmom[4]) + 1/(2.*Pair(m_pmom[1],m_pmom[4])))/
               (2.*(MP2 + Pair(m_pmom[2],m_pmom[3]))) - 
              1/(2.*(ME2 + Pair(m_pmom[0],m_pmom[1]))*Pair(m_pmom[3],m_pmom[4])))*
            (Pair(m_pmom[0],m_pmom[4])*(MP2*Pair(m_pmom[0],m_pmom[1]) - 
                 2*Pair(m_pmom[0],m_pmom[2])*Pair(m_pmom[1],m_pmom[2]) + 
                 ME2*(Pair(m_pmom[0],m_pmom[2]) + Pair(m_pmom[1],m_pmom[2]))) + 
              (4*Pair(m_pmom[0],m_pmom[2])*Pair(m_pmom[1],m_pmom[2]) - 
                 2*(MP2*Pair(m_pmom[0],m_pmom[1]) + 
                    ME2*(Pair(m_pmom[0],m_pmom[2]) + Pair(m_pmom[1],m_pmom[2]))))*
               Pair(m_pmom[2],m_pmom[4]) + 
              ME2*(-((2*MP2 - 2*(Pair(m_pmom[0],m_pmom[2]) + Pair(m_pmom[0],m_pmom[3])) - 
                      Pair(m_pmom[0],m_pmom[4]))*Pair(m_pmom[2],m_pmom[3])) - 
                 (-2*(Pair(m_pmom[0],m_pmom[2]) + Pair(m_pmom[0],m_pmom[3])) - 
                    Pair(m_pmom[0],m_pmom[4]) + 2*Pair(m_pmom[2],m_pmom[3]))*Pair(m_pmom[2],m_pmom[4]) - 
                 2*(pow(Pair(m_pmom[2],m_pmom[3]),2) + MP2*Pair(m_pmom[2],m_pmom[4])) - 
                 2*(Pair(m_pmom[2],m_pmom[3]) + Pair(m_pmom[2],m_pmom[4]))*Pair(m_pmom[3],m_pmom[4]))))/
          (4.*Pair(m_pmom[0],m_pmom[4])*(MP2 + Pair(m_pmom[2],m_pmom[3])))) + 
      (ME2*(-512*pow(-2/Pair(m_pmom[0],m_pmom[4]) - 5/(2.*Pair(m_pmom[1],m_pmom[4])),2)*
             pow(Pair(m_pmom[1],m_pmom[4]),2) + 
            (128*(-pow(Pair(m_pmom[0],m_pmom[4]),2) + 
                 4*Pair(m_pmom[0],m_pmom[4])*Pair(m_pmom[2],m_pmom[4]) - 
                 4*pow(Pair(m_pmom[2],m_pmom[4]),2)))/pow(Pair(m_pmom[0],m_pmom[4]),2)) - 
         ((-1/2.*1/Pair(m_pmom[0],m_pmom[4]) - 1/(2.*Pair(m_pmom[1],m_pmom[4])))*
            (-1024*ME2*(-(Pair(m_pmom[0],m_pmom[4])*
                    (4*MP2 - 6*Pair(m_pmom[0],m_pmom[2]) - 4*Pair(m_pmom[0],m_pmom[3]) + 
                      5*Pair(m_pmom[2],m_pmom[3]))) - 
                 (-8*MP2 + 12*Pair(m_pmom[0],m_pmom[2]) + 8*Pair(m_pmom[0],m_pmom[3]) - 
                    10*Pair(m_pmom[2],m_pmom[3]))*Pair(m_pmom[2],m_pmom[4])) - 
              512*(-(pow(Pair(m_pmom[0],m_pmom[4]),2)*Pair(m_pmom[1],m_pmom[2])) - 
                 2*(pow(Pair(m_pmom[0],m_pmom[2]),2)*Pair(m_pmom[0],m_pmom[4]) + 
                    Pair(m_pmom[0],m_pmom[2])*
                     (pow(Pair(m_pmom[0],m_pmom[4]),2) + 
                       Pair(m_pmom[0],m_pmom[4])*Pair(m_pmom[1],m_pmom[2]))) + 
                 (4*pow(Pair(m_pmom[0],m_pmom[2]),2) + 
                    Pair(m_pmom[0],m_pmom[2])*
                     (-4*(ME2 + Pair(m_pmom[0],m_pmom[1])) + 6*Pair(m_pmom[0],m_pmom[4]) + 
                       4*Pair(m_pmom[1],m_pmom[2])) + 
                    2*(Pair(m_pmom[0],m_pmom[4])*Pair(m_pmom[1],m_pmom[2]) + 
                       Pair(m_pmom[0],m_pmom[2])*Pair(m_pmom[1],m_pmom[4])))*Pair(m_pmom[2],m_pmom[4]) 
+ Pair(m_pmom[0],m_pmom[2])*(Pair(m_pmom[0],m_pmom[4])*
                     (2*(ME2 + Pair(m_pmom[0],m_pmom[1])) - Pair(m_pmom[1],m_pmom[4])) - 
                    4*pow(Pair(m_pmom[2],m_pmom[4]),2)))))/(2.*Pair(m_pmom[0],m_pmom[4])) + 
         (-2/Pair(m_pmom[0],m_pmom[4]) - 5/(2.*Pair(m_pmom[1],m_pmom[4])))*Pair(m_pmom[1],m_pmom[4])*
          ((512*ME2*(2*(ME2 + Pair(m_pmom[0],m_pmom[1])) - 3*Pair(m_pmom[0],m_pmom[4]) - 
                 2*(Pair(m_pmom[0],m_pmom[2]) + Pair(m_pmom[1],m_pmom[2]) + 
                    Pair(m_pmom[1],m_pmom[4])) + 4*Pair(m_pmom[2],m_pmom[4])))/Pair(m_pmom[0],m_pmom[4]) 
+ (-1/2.*1/Pair(m_pmom[0],m_pmom[4]) - 1/(2.*Pair(m_pmom[1],m_pmom[4])))*
             (-512*(-2*pow(Pair(m_pmom[1],m_pmom[2]),2) + 
                  (-Pair(m_pmom[0],m_pmom[2]) - 2*Pair(m_pmom[1],m_pmom[2]))*
                   Pair(m_pmom[1],m_pmom[4]) + 
                  Pair(m_pmom[1],m_pmom[2])*
                   (-2*Pair(m_pmom[0],m_pmom[2]) - Pair(m_pmom[0],m_pmom[4]) + 
                     2*(ME2 + Pair(m_pmom[0],m_pmom[1]) + Pair(m_pmom[2],m_pmom[4])))) + 
               1024*ME2*(-6*Pair(m_pmom[1],m_pmom[2]) - 4*Pair(m_pmom[1],m_pmom[3]) + 
                  5*Pair(m_pmom[2],m_pmom[3]) + 6*(MP2 + Pair(m_pmom[2],m_pmom[4])) + 
                  4*Pair(m_pmom[3],m_pmom[4])))) + 
         pow(-1/2.*1/Pair(m_pmom[0],m_pmom[4]) - 1/(2.*Pair(m_pmom[1],m_pmom[4])),2)*
          (-512*(-(((6*ME2 + 2*MP2 + 6*Pair(m_pmom[0],m_pmom[1]))*
                     Pair(m_pmom[0],m_pmom[2]) + 
                    (4*ME2 + 3*MP2 + 4*Pair(m_pmom[0],m_pmom[1]) - 
                       7*Pair(m_pmom[0],m_pmom[2]))*Pair(m_pmom[0],m_pmom[3]) - 
                    4*pow(Pair(m_pmom[0],m_pmom[3]),2))*Pair(m_pmom[1],m_pmom[2])) - 
               Pair(m_pmom[0],m_pmom[3])*pow(Pair(m_pmom[1],m_pmom[2]),2) - 
               (ME2 + Pair(m_pmom[0],m_pmom[1]))*
                (4*MP2*Pair(m_pmom[0],m_pmom[2]) - 6*pow(Pair(m_pmom[0],m_pmom[2]),2) - 
                  4*(Pair(m_pmom[0],m_pmom[2])*Pair(m_pmom[0],m_pmom[3]) + 
                     MP2*Pair(m_pmom[1],m_pmom[2]))) - 
               (-5*MP2*Pair(m_pmom[0],m_pmom[2]) + 7*pow(Pair(m_pmom[0],m_pmom[2]),2) + 
                  Pair(m_pmom[0],m_pmom[2])*(4*Pair(m_pmom[0],m_pmom[3]) - Pair(m_pmom[1],m_pmom[2])))*
                Pair(m_pmom[1],m_pmom[3]) - 
               (MP2*Pair(m_pmom[0],m_pmom[1]) + 
                  (6*ME2 + 5*Pair(m_pmom[0],m_pmom[1]))*Pair(m_pmom[0],m_pmom[2]) + 
                  (-4*ME2 - 5*Pair(m_pmom[0],m_pmom[1]) + 4*Pair(m_pmom[0],m_pmom[3]))*
                   Pair(m_pmom[1],m_pmom[2]) - 6*Pair(m_pmom[0],m_pmom[2])*Pair(m_pmom[1],m_pmom[3]))*
                Pair(m_pmom[2],m_pmom[3]) + 
               ME2*(-48*MP2*Pair(m_pmom[0],m_pmom[2]) + 
                  36*pow(Pair(m_pmom[0],m_pmom[2]),2) + 
                  (-32*MP2 + 48*Pair(m_pmom[0],m_pmom[2]))*Pair(m_pmom[0],m_pmom[3]) + 
                  16*(MP2*MP2 + pow(Pair(m_pmom[0],m_pmom[3]),2)) + 
                  (39*MP2 - 59*Pair(m_pmom[0],m_pmom[2]) - 40*Pair(m_pmom[0],m_pmom[3]) + 
                     Pair(m_pmom[1],m_pmom[2]))*Pair(m_pmom[2],m_pmom[3]) + 
                  24*pow(Pair(m_pmom[2],m_pmom[3]),2))) + 
            512*(MP2*MP2*Pair(m_pmom[0],m_pmom[1]) + ME2*MP2*Pair(m_pmom[3],m_pmom[4])) + 
            MP2*(-512*(Pair(m_pmom[0],m_pmom[2])*Pair(m_pmom[1],m_pmom[2]) + 
                  Pair(m_pmom[0],m_pmom[3])*Pair(m_pmom[1],m_pmom[3])) + 
               512*(-(Pair(m_pmom[0],m_pmom[4])*Pair(m_pmom[1],m_pmom[3])) - 
                  Pair(m_pmom[0],m_pmom[3])*Pair(m_pmom[1],m_pmom[4]) + 
                  Pair(m_pmom[0],m_pmom[1])*Pair(m_pmom[3],m_pmom[4])))))/
       (4.*pow(MP2 + Pair(m_pmom[2],m_pmom[3]),2)) + 
      ((-1/2.*1/Pair(m_pmom[0],m_pmom[4]) + 1/(2.*Pair(m_pmom[1],m_pmom[4])))/
          (2.*(MP2 + Pair(m_pmom[2],m_pmom[3]))) + 
         1/(2.*(ME2 + Pair(m_pmom[0],m_pmom[1]))*Pair(m_pmom[2],m_pmom[4])))*
       (-512*((1/(2.*(ME2 + Pair(m_pmom[0],m_pmom[1]))) + 
               (MP2*(-1/2.*1/Pair(m_pmom[0],m_pmom[4]) - 1/(2.*Pair(m_pmom[1],m_pmom[4]))))/
                (2.*(MP2 + Pair(m_pmom[2],m_pmom[3]))))*
             (-(Pair(m_pmom[0],m_pmom[1])*
                  (MP2 - Pair(m_pmom[0],m_pmom[2]) - Pair(m_pmom[1],m_pmom[2]))) - 
               Pair(m_pmom[0],m_pmom[3])*Pair(m_pmom[1],m_pmom[2]) - 
               Pair(m_pmom[0],m_pmom[2])*Pair(m_pmom[1],m_pmom[3]) + 
               ME2*(Pair(m_pmom[2],m_pmom[3]) + Pair(m_pmom[2],m_pmom[4]))) + 
            ((-2/Pair(m_pmom[0],m_pmom[4]) - 5/(2.*Pair(m_pmom[1],m_pmom[4])))*
               (ME2*Pair(m_pmom[1],m_pmom[4])*(Pair(m_pmom[2],m_pmom[3]) + Pair(m_pmom[2],m_pmom[4])) + 
                 Pair(m_pmom[1],m_pmom[4])*
                  (-(MP2*Pair(m_pmom[0],m_pmom[1])) - 
                    ME2*(Pair(m_pmom[0],m_pmom[2]) - Pair(m_pmom[1],m_pmom[2])) - 
                    (-2*Pair(m_pmom[0],m_pmom[1]) + Pair(m_pmom[0],m_pmom[4]))*
                     Pair(m_pmom[1],m_pmom[2]) - 2*pow(Pair(m_pmom[1],m_pmom[2]),2) - 
                    (-Pair(m_pmom[0],m_pmom[2]) + 2*Pair(m_pmom[1],m_pmom[2]))*
                     Pair(m_pmom[1],m_pmom[4]) + 2*Pair(m_pmom[1],m_pmom[2])*Pair(m_pmom[2],m_pmom[4]))))/
             (2.*(MP2 + Pair(m_pmom[2],m_pmom[3])))) + 
         (512*(-1/2.*1/Pair(m_pmom[0],m_pmom[4]) - 1/(2.*Pair(m_pmom[1],m_pmom[4])))*
             (((3*MP2 + 4*Pair(m_pmom[0],m_pmom[1]) - 5*Pair(m_pmom[0],m_pmom[2]))*
                   Pair(m_pmom[0],m_pmom[3]) - 4*pow(Pair(m_pmom[0],m_pmom[3]),2))*
                Pair(m_pmom[1],m_pmom[2]) + 
               Pair(m_pmom[0],m_pmom[3])*pow(Pair(m_pmom[1],m_pmom[2]),2) + 
               Pair(m_pmom[0],m_pmom[1])*
                (4*(MP2*MP2) - 10*MP2*Pair(m_pmom[0],m_pmom[2]) + 
                  6*pow(Pair(m_pmom[0],m_pmom[2]),2) + 
                  (-4*MP2 + 4*Pair(m_pmom[0],m_pmom[2]))*Pair(m_pmom[0],m_pmom[3]) + 
                  (-4*MP2 + 6*Pair(m_pmom[0],m_pmom[2]))*Pair(m_pmom[1],m_pmom[2])) + 
               (-7*pow(Pair(m_pmom[0],m_pmom[2]),2) + 
                  Pair(m_pmom[0],m_pmom[2])*
                   (5*MP2 - 4*Pair(m_pmom[0],m_pmom[3]) - Pair(m_pmom[1],m_pmom[2])))*
                Pair(m_pmom[1],m_pmom[3]) + 
               (4*Pair(m_pmom[0],m_pmom[3])*Pair(m_pmom[1],m_pmom[2]) + 
                  Pair(m_pmom[0],m_pmom[1])*
                   (5*MP2 - 5*(Pair(m_pmom[0],m_pmom[2]) + Pair(m_pmom[1],m_pmom[2]))) + 
                  6*Pair(m_pmom[0],m_pmom[2])*Pair(m_pmom[1],m_pmom[3]))*Pair(m_pmom[2],m_pmom[3]) + 
               ME2*(-((4*MP2 - 6*Pair(m_pmom[0],m_pmom[2]) - 4*Pair(m_pmom[0],m_pmom[3]))*
                     Pair(m_pmom[2],m_pmom[3])) - 5*pow(Pair(m_pmom[2],m_pmom[3]),2) - 
                  (4*MP2 - 6*Pair(m_pmom[0],m_pmom[2]) - 4*Pair(m_pmom[0],m_pmom[3]) + 
                     5*Pair(m_pmom[2],m_pmom[3]))*Pair(m_pmom[2],m_pmom[4]))) - 
            (512*(-((MP2*Pair(m_pmom[0],m_pmom[1]) + 
                       (-ME2 - 2*Pair(m_pmom[0],m_pmom[1]))*Pair(m_pmom[0],m_pmom[2]) + 
                       2*pow(Pair(m_pmom[0],m_pmom[2]),2))*Pair(m_pmom[0],m_pmom[4])) - 
                  2*Pair(m_pmom[0],m_pmom[2])*pow(Pair(m_pmom[0],m_pmom[4]),2) - 
                  (ME2*Pair(m_pmom[0],m_pmom[4]) - pow(Pair(m_pmom[0],m_pmom[4]),2))*
                   Pair(m_pmom[1],m_pmom[2]) - 
                  Pair(m_pmom[0],m_pmom[2])*Pair(m_pmom[0],m_pmom[4])*Pair(m_pmom[1],m_pmom[4]) - 
                  (-2*MP2*Pair(m_pmom[0],m_pmom[1]) + 
                     (2*ME2 + 4*Pair(m_pmom[0],m_pmom[1]))*Pair(m_pmom[0],m_pmom[2]) - 
                     4*pow(Pair(m_pmom[0],m_pmom[2]),2) - 
                     6*Pair(m_pmom[0],m_pmom[2])*Pair(m_pmom[0],m_pmom[4]) + 
                     (-2*ME2 + 2*Pair(m_pmom[0],m_pmom[4]))*Pair(m_pmom[1],m_pmom[2]) - 
                     2*Pair(m_pmom[0],m_pmom[2])*Pair(m_pmom[1],m_pmom[4]))*Pair(m_pmom[2],m_pmom[4]) - 
                  4*Pair(m_pmom[0],m_pmom[2])*pow(Pair(m_pmom[2],m_pmom[4]),2)) - 
               512*ME2*(-((2*MP2 - 
                       2*(Pair(m_pmom[0],m_pmom[2]) + Pair(m_pmom[0],m_pmom[3])) - 
                       Pair(m_pmom[0],m_pmom[4]))*Pair(m_pmom[2],m_pmom[3])) - 
                  (-2*(Pair(m_pmom[0],m_pmom[2]) + Pair(m_pmom[0],m_pmom[3])) - 
                     Pair(m_pmom[0],m_pmom[4]) + 2*Pair(m_pmom[2],m_pmom[3]))*Pair(m_pmom[2],m_pmom[4]) - 
                  2*(pow(Pair(m_pmom[2],m_pmom[3]),2) + MP2*Pair(m_pmom[2],m_pmom[4])) - 
                  2*(Pair(m_pmom[2],m_pmom[3]) + Pair(m_pmom[2],m_pmom[4]))*Pair(m_pmom[3],m_pmom[4])))/
             (2.*Pair(m_pmom[0],m_pmom[4])))/(2.*(MP2 + Pair(m_pmom[2],m_pmom[3])))));
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
  PRINT_INFO(fl);
  if(fl.size()!=5) return NULL;
  return new PionPionGamma(args);
  // if(fl[2])
  // }
  // return NULL;
}
