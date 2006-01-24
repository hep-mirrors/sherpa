//bof
//Version: 3 ADICIC++-0.0/2005/09/21

//Possibility of having information at compile time.
//Globally defined parameter sets influencing the calculation of recoils.



#ifndef _Recoil_Strategy_hpp_
#define _Recoil_Strategy_hpp_ _Recoil_Strategy_hpp_


#include <vector>
#include <enumextra>
#include "Vector.H"
#include "Poincare.H"





namespace ADICIC {



  namespace Recoil_Strategy {

    //The settings here must agree with the ones in ADICIC::MakeRecos(..),
    //see Recoil_Calculator.H

    enum Type {
      Unknown    = 0,
      Kleiss     = 2,
      FixDir1    = 1,
      FixDir3    = 3,
      MinimizePt = 4,
      Lonnblad   = 5,
      OldAdicic  = 6,
      Test       = 7,
      Ktii       = 8,
      stop       = -99
    };

    static const size_t NumberOfTypes = 10;

    static const Type List[NumberOfTypes]={
      Unknown,
      Kleiss,
      FixDir1,
      FixDir3,
      MinimizePt,
      Lonnblad,
      OldAdicic,
      Test,
      Ktii,
      stop    //Keep this at the bottom.
    };

  }    //eo namespace Recoil_Strategy





  //Collection of fourvectors.
  typedef std::vector<ATOOLS::Vec4D> Multi4vec;





  //Enhance readability! Do not change numbers.
  //Dipole shower mode dsm (j i.e. just).
  namespace dsm {
    enum code {
      off  = 0,
      jff  = 1,
      jfi  = 2,
      jif  = 4,
      jii  = 8,
      iiff = 9,
      all  = 15
    };
  }
  //Recoil labels rl (qbar i.e. a).
  namespace rl {
    enum code {
      qag   = 0,
      qgg   = 1,
      gag   = 2,
      ggg   = 3,
      iiaqg = 4,
      iiagg = 5,
      iigqg = 6,
      iiggg = 7,
      //-----
      qga   = 8,
      gaq   = 9,
      gga   = 10,
      ggq   = 11,
      iiaqa = 12,
      iiaqq = 13,
      iiagq = 14,
      iigqa = 15,
      //-----
      stop  = 16
    };
  }
  //Recoil results rr.
  namespace rr {
    enum code {
      p1   = 0,    //TopBranch.
      //----
      p2   = 1,
      axis = 1,
      //----
      p3   = 2,    //BotBranch.
      stop = 3
    };
  }
  //Recoil distribution type rdt.
  namespace rdt {
    enum code {
      local    = 0,
      iirecoil = 1,
      stop     = 2
    };
  }





  //Interfacing the results of the Recoil calculation.
  struct Recoil_Result {
    xbool     Poc;    //Partner of compensation.
    Multi4vec Vec;
    //------------------
    Recoil_Result();
    ~Recoil_Result() {}
    //------------------
    inline void Reset();
    void Print() const;
  };





  //Interfacing the results of the Recoil distribution.
  class Recoil_Tool {
  private:
    rdt::code         m_mode;
    ATOOLS::Poincare* p_fly;
    ATOOLS::Poincare* p_flyprime;
  public:
    Recoil_Tool(rdt::code c=rdt::local);
    ~Recoil_Tool();
    inline const rdt::code Mode() const;
    inline bool KeepThisBoost(ATOOLS::Poincare&);
    inline bool KeepThisBackBoost(ATOOLS::Poincare&);
    inline ATOOLS::Poincare& GetBoost() const;
    inline ATOOLS::Poincare& GetBackBoost() const;
  };





  struct TEMP {
    static bool CPTEST;
  };



}    //eo namespace ADICIC





#include "Recoil_Strategy.inl.hh"


#endif    //eo _Recoil_Strategy_hpp_



//eof
