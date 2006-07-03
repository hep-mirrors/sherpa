//bof
//Version: 4 ADICIC++-0.0/2006/05/26

//Encoded information, carriers and interface structures.



#ifndef _Recoil_Utils_hpp_
#define _Recoil_Utils_hpp_ _Recoil_Utils_hpp_


#include <vector>
#include <enumextra>
#include "Vector.H"
#include "Poincare.H"





namespace ADICIC {



  //Collection of fourvectors.
  typedef std::vector<ATOOLS::Vec4D> Multi4vec;





  //Enhance readability! Do not change numbers.
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





  //Interfacing the transformations appearing in the Recoil kinematics.
  class Multipoincare {
  private:
    typedef std::pair<int,ATOOLS::Poincare*> Set;
  private:    //Blocked!
    Multipoincare(const Multipoincare&);
    Multipoincare& operator=(const Multipoincare&);
  private:
    rdt::code         m_mode;
    std::vector<Set>  v_trafs;
  private:    //Methods.
    inline void Destruct();
  public:
    Multipoincare(rdt::code c=rdt::local);
    ~Multipoincare();
  public:    //Methods.
    inline const rdt::code Mode() const;
    inline const std::vector<Set>& GetItsVec() const;
    inline void Swap(Multipoincare&);
    inline void Keep(ATOOLS::Poincare&, bool, bool);
    const  bool Apply(ATOOLS::Vec4D&) const;
    void Clear();
    void Clear(rdt::code);
  };





  //Interfacing the results of the Recoil calculation.
  struct Recoil_Result {
    xbool         Poc;    //Partner of compensation.
    Multi4vec     Vec;
    Multipoincare Mup;
    //------------------
    Recoil_Result(rdt::code c=rdt::local);
    ~Recoil_Result() {}
    //------------------
    inline void Reset(rdt::code c=rdt::local);
    void Print() const;
  };





  struct TEMP {
    static bool CPTEST;
  };



}    //eo namespace ADICIC





#include "Recoil_Utils.inl.hh"


#endif    //eo _Recoil_Utils_hpp_



//eof
