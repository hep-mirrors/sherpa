//bof
//Version: 2 ADICIC++-0.0/2004/08/06

//Inline methods of Dipole_Parameter.H.





//#include <...>
//#include "..."





namespace ADICIC {



  //===========================================================================



  inline const bool Dipole_Parameter::IsAlphaSRunning() {    //Static.
    return s_isalphasrun;
  }
  inline const double Dipole_Parameter::AlphaSFix() {    //Static.
    return s_alphasfix;
  }
  inline const double Dipole_Parameter::MinOfK2t() {    //Static.
    return s_k2tmin;
  }
  inline const double Dipole_Parameter::MaxOfK2t() {    //Static.
    return s_k2tmax;
  }



  inline const int Dipole_Parameter::RecoilStrategyQQbar() {    //Static.
    return s_restratqqbar;
  }
  inline const int Dipole_Parameter::RecoilStrategyQG() {    //Static.
    return s_restratqg;
  }
  inline const int Dipole_Parameter::RecoilStrategyGQbar() {    //Static.
    return s_restratgqbar;
  }
  inline const int Dipole_Parameter::RecoilStrategyGG() {    //Static.
    return s_restratgg;
  }



  inline const int Dipole_Parameter::ChainEvolutionStrategy() {    //Static.
    return s_chevolstrat;
  }



  //===========================================================================



}    //eo namespace ADICIC





//eof
