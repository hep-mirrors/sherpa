//bof
//Version: 3 ADICIC++-0.0/2005/09/29

//Possibility of having information at compile time.
//Carriers or interface structures.
//Globally defined parameters influencing the Sudakov calculation.



#ifndef _Sudakov_Utils_hpp_
#define _Sudakov_Utils_hpp_ _Sudakov_Utils_hpp_


#include <string>
#include <enumextra>
#include <map>
#include <vector>
#include "Information.hpp"





namespace ADICIC {



  namespace Radiation {

    enum Group {
      gluon    = 10,    //F gluon emission.
      igluon   = 20,    //I gluon emission.
      qfront   = 21,
      qbarend  = 22,
      quark    = 50,    //Gluon splitting.
      qtop     = 51,    //FF.
      qbarbot  = 52,    //FF.
      qitop    = 53,    //II.
      qibot    = 54,    //II.
      incorrect=-9999
    };
    enum Type {
      d      = 1,
      du     = 2,
      dus    = 3,
      dusc   = 4,
      duscb  = 5,
      g      = 10,
      gd     = 11,
      gdu    = 12,
      gdus   = 13,
      gdusc  = 14,
      gduscb = 15
    };

  }





  //Collection of flavours and double numbers.
  typedef std::vector<ATOOLS::Flavour> Multiflavour;
  typedef std::vector<double>          Multidouble;





  //Enhance readability. Do not change numbers.
  namespace sf {
    enum code {
      plusini = 0,
      miusini = 1,
      plusfin = 2,
      miusfin = 3,
      stop    = 4
    };
  }
  namespace sr {
    enum code {
      xpini   = 0,
      xmini   = 1,
      mdip    = 2,
      shatmax = 3,
      expydip = 4,
      shat    = 5,
      kt      = 6,
      mt      = 7,
      expy    = 8,
      xpfin   = 9,
      xmfin   = 10,
      stop    = 11
    };
    class stringmap {
      std::map<code,std::string> m_map;
    public:
      stringmap();
      const std::string& operator[](code);
    };
    static stringmap name = stringmap();
  }





  //Coding the flavour of the radiation.
  struct Sudakov_Flavour {
    const Dipole_Gluon_Base*     Glu;
    const Dipole_Quark_Base*     Qua;
    const Dipole_Antiquark_Base* Aqu;
    //-------------------------------
    Sudakov_Flavour();
    ~Sudakov_Flavour() {}
    //-------------------------------
    inline void Nil();
    inline const bool Gluic() const;
    inline const bool Bqqic() const;
    inline const bool Quaic() const;
    inline const bool Aquic() const;
    void Print() const;
  };


  inline void Sudakov_Flavour::Nil() {
    Glu=NULL; Qua=NULL; Aqu=NULL;
  }
  inline const bool Sudakov_Flavour::Gluic() const {
    return Glu && !Qua && !Aqu;
  }
  inline const bool Sudakov_Flavour::Bqqic() const {
    if(!Glu && Qua && Aqu) return (*Qua)().Kfcode()==(*Aqu)().Kfcode();
    return false;
  }
  inline const bool Sudakov_Flavour::Quaic() const {
    return !Glu && Qua && !Aqu;
  }
  inline const bool Sudakov_Flavour::Aquic() const {
    return !Glu && !Qua && Aqu;
  }





  //Interfacing the results of the Sudakov calculation.
  struct Sudakov_Result {
    Radiation::Group Rad;
    Sudakov_Flavour  Sfc;
    double           P2t;
    double           Y;
    double           X1;
    double           X3;
    Multidouble      Isr;
    bool             Dir;   //Dipole direction: true..+to-, false..-to+ beam.
    //-------------------
    Sudakov_Result();
    ~Sudakov_Result() {}
    //-------------------
    inline void Reset();
    void Print() const;
  };


  inline void Sudakov_Result::Reset() {
    Rad=Radiation::incorrect;
    Sfc.Nil();
    P2t=0.0; Y=0.0; X1=1.0; X3=1.0;
    Isr.clear();
    Dir=true;
  }



}    //eo namespace ADICIC





#endif    //eo _Sudakov_Utils_hpp_



//eof
