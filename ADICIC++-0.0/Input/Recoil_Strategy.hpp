//bof
//Version: 2 ADICIC++-0.0/2004/08/09

//Possibility of having information at compile time.
//Globally defined parameter sets influencing Dipole_Handler.H.



#ifndef _Recoil_Strategy_hpp_
#define _Recoil_Strategy_hpp_ _Recoil_Strategy_hpp_





namespace ADICIC {



  namespace Recoil_Strategy {



    //Set the corresponding values:
    static const int Label_qqbar = 2;
    static const int Label_qg    = 3;
    static const int Label_gqbar = 1;
    static const int Label_gg    = 4;



    struct Unknown {};       //All other Labels.
    struct FixDir1 {};       //Label==1.
    struct Kleiss {};        //Label==2.
    struct FixDir3 {};       //Label==3.
    struct MinimizePt {};    //Label==4.
    struct Lonnblad {};      //Label==5.
    struct OldAdicic {};     //Label==6.
    struct Test {};          //Label==7.



    template<int label> struct Map {
      typedef Unknown Ret;
    };
    template<> struct Map<1> {
      typedef FixDir1 Ret;
    };
    template<> struct Map<2> {
      typedef Kleiss Ret;
    };
    template<> struct Map<3> {
      typedef FixDir3 Ret;
    };
    template<> struct Map<4> {
      typedef MinimizePt Ret;
    };
    template<> struct Map<5> {
      typedef Lonnblad Ret;
    };
    template<> struct Map<6> {
      typedef OldAdicic Ret;
    };
    template<> struct Map<7> {
      typedef Test Ret;
    };



    typedef Map<Label_qqbar>::Ret Ret_qqbar;
    typedef Map<Label_qg>::Ret    Ret_qg;
    typedef Map<Label_gqbar>::Ret Ret_gqbar;
    typedef Map<Label_gg>::Ret    Ret_gg;



  };    //eo namespace Recoil_Strategy





  struct TEMP {
    static bool CPTEST;
  };



}    //eo namespace ADICIC





#endif    //eo _Recoil_Strategy_hpp_



//eof
