//bof
//Version: 2 ADICIC++-0.0/2004/08/10    (surely preliminarily)

//A class is provided, which globally manages the Adicic parameters.



#ifndef _Dipole_Parameter_H_
#define _Dipole_Parameter_H_ _Dipole_Parameter_H_


//#include <...>
//#include "..."


#define DIPOLE_PARAMETER_OUTPUT DIPOLE_PARAMETER_OUTPUT
//#undef  DIPOLE_PARAMETER_OUTPUT





namespace ADICIC {



  class Dipole_Parameter_Init;





  class Dipole_Parameter {

    friend class Dipole_Parameter_Init;

  private:

    static const bool sf_start;

    static bool   s_isalphasrun;
    static double s_alphasfix;
    static double s_k2tmin;
    static double s_k2tmax;

    static int s_restratqqbar;
    static int s_restratqg;
    static int s_restratgqbar;
    static int s_restratgg;

    static int s_chevolstrat;

    static const bool Reset();

  public:

    Dipole_Parameter();
    ~Dipole_Parameter();

    static inline const bool   IsAlphaSRunning();
    static inline const double AlphaSFix();
    static inline const double MinOfK2t();
    static inline const double MaxOfK2t();

    static inline const int RecoilStrategyQQbar();
    static inline const int RecoilStrategyQG();
    static inline const int RecoilStrategyGQbar();
    static inline const int RecoilStrategyGG();

    static inline const int ChainEvolutionStrategy();

    static       void Show();
    static const bool Check(const double cmsenergy);
    static const bool ForceFirstInit();

  };



}    //eo namespace ADICIC





#include "Dipole_Parameter.inl.hh"


#endif    //eo _Dipole_Parameter_H_



//eof
