//bof
//Version: 4 ADICIC++-0.0/2006/06/30    (surely preliminarily)

//A class is provided, which globally manages the Adicic parameters.



#ifndef _Dipole_Parameter_H_
#define _Dipole_Parameter_H_ _Dipole_Parameter_H_


#include <vector>
#include "Sudakov_Utils.hpp"
#include "Recoil_Strategy.hpp"
#include "Recoil_Utils.hpp"
#include "Evolution_Strategy.hpp"


#define DIPOLE_PARAMETER_OUTPUT DIPOLE_PARAMETER_OUTPUT
//#undef  DIPOLE_PARAMETER_OUTPUT





namespace ADICIC {



  class Dipole_Parameter_Init;





  class Dipole_Parameter {

    friend class Dipole_Parameter_Init;

  private:

    static const bool sf_start;

    static const bool SetWithStatics();

  public:

    Dipole_Parameter();
    ~Dipole_Parameter();


    class Sud {
      friend class Dipole_Parameter;
      friend class Dipole_Parameter_Init;
    private:
      static bool            s_runalphas;
      static double          s_alphasfix;
      static unsigned short  s_nffix;    //No effect if AlphaS is running.
      static xbool           s_gsplit;
      static Radiation::Type s_radiatype;
      static double          s_k2tmin;    //FF dipoles.
      static double          s_k2tmax;
      static double          s_k2tiimin;    //II dipoles.
      static double          s_k2tiimax;
      static double          s_iieffexp;
      //-------------------------------------
      static double          s_k2tiifac;
      static double          s_k2tiifixscale;
      static double          s_k2tiivarscale;
    public:
      inline const bool            RunAlphaS() const;
      inline const double          AlphaSFix() const;
      inline const unsigned short  NfFix() const;
      inline const xbool           GsplitRule() const;
      inline const Radiation::Type RadiationType() const;
      inline const double          MinK2t() const;
      inline const double          MaxK2t() const;
      inline const double          MinIIK2t() const;
      inline const double          MaxIIK2t() const;
      inline const double          IIEffExp() const;
      //-------------------------------------------------
      inline void                  SetMaxIIScale(double d=-1.0);
    } sud;


    class Kin {
      friend class Dipole_Parameter;
      friend class Dipole_Parameter_Init;
    private:
      static int                                s_dsmode;
      static std::vector<Recoil_Strategy::Type> v_recostrat;
    public:
      inline const int                                 ShowerMode() const;
      inline const std::vector<Recoil_Strategy::Type>& RecoilStrategy() const;
    } kin;


    class Evo {
      friend class Dipole_Parameter;
      friend class Dipole_Parameter_Init;
    private:
      static fascat::code                                s_fascatype;
      static double                                      s_fmuf;
      static double                                      s_fmur;
      static std::vector<Chain_Evolution_Strategy::Type> v_chevostrat;
      //--------------------------------------------------------------
      static std::size_t                                 s_chpartlim;
      static std::size_t                                 s_chcorrlim;
      static double                                      s_scaoffset;
    public:
      inline const fascat::code FactScaleType() const;
      inline const double GetFactScaleFrom(const Multidouble&) const;
      inline const double GetFactScaleFrom(const Sudakov_Result&) const;
      inline const double FactScaleFactor() const;
      inline const double RenoScaleFactor() const;
      inline const std::vector<Chain_Evolution_Strategy::Type>&
      ChainEvolutionStrategy() const;
      inline const std::size_t ChainParticleLimit() const;
      inline const std::size_t ChainCorrelationLimit() const;
      //-------------------------------------------------------
      inline void SetChainParticleLimit(std::size_t n=7777777);
      inline void SetChainCorrelationLimit(std::size_t n=7777777);
      inline void SetFactScaleOffset(double d);
    } evo;


    static       void Show();
    static const bool Check(const double cmsenergy);

    inline const double MaxIISHat(double s) const;
    inline const double MaxIIInvScale(double s) const;
    inline const double HighestIIInvScale(double s) const;

  };

  extern const Dipole_Parameter dpa;
  extern       Dipole_Parameter dpv;



}    //eo namespace ADICIC





#include "Dipole_Parameter.inl.hh"


#endif    //eo _Dipole_Parameter_H_



//eof
