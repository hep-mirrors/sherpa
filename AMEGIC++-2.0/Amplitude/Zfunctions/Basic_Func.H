#ifndef Basic_Func_h
#define Basic_Func_h

#include "Kabbala.H"
#include "Pfunc.H"
#include "Basic_Sfuncs.H"

namespace AMEGIC {

  class Virtual_String_Generator;

  class Spinor {
  public:
    enum code {None = 0,
               u    = 1,
	       ubar = 2,
               v    = 3,
               vbar = 4, 
               Unnown=99
    };
  };

  class Direction {
  public:
    enum code {None     =  0,
	       Incoming = -1,
               Outgoing =  1,
               Unknown  = 99
    };
  }; 

  class Argument {
  public:
    APHYTOOLS::kf::code        kfcode;
    Spinor::code    spinortype;
    Direction::code direction;
    int             numb;
    bool            maped;
    Argument() {
      spinortype = Spinor::None;
      direction  = Direction::Outgoing;
      numb       = -99;
      maped      = false;
    }   
  };

  inline bool operator==(const Argument& a, const Argument& b) 
  {
    if (a.spinortype!=b.spinortype) return false;
    if (a.direction!=b.direction) return false;
    if (a.numb!=b.numb) return false;

    return true;
  }

  inline bool operator!=(const Argument& a, const Argument& b) {return !(a==b);}


  class Basic_Func {
  public:
    Argument*                 arg;
    Complex*                  coupl;
    Argument*                 ps;
    std::list<Pfunc*>*        pl;
    Virtual_String_Generator* sgen;
    Basic_Sfuncs*             BS;

    Basic_Func(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : sgen(_sgen), BS(_BS) {}

    virtual Kabbala X(const int&,const int&,const int&,
 		      const int&,const int&,
 		      const Complex&,const Complex&) {
      std::cerr<<"Calling Basic_Func::X"<<std::endl;
      return Kabbala();
    }
    virtual Kabbala V(const int& a,const int &b)
    {
      std::cerr<<"Calling Basic_Func::V"<<std::endl;
      return Kabbala();
    }
    virtual Kabbala Vcplx(const int& a,const int &b)
    {
      std::cerr<<"Calling Basic_Func::Vcplx"<<std::endl;
      return Kabbala();
    }

    void SetArgCouplProp(int,Argument*,Complex*,int,Argument*,std::list<Pfunc*>*);
    void Map(int&,bool&);
  };

  class Basic_Yfunc : public virtual Basic_Func {
  public:
    Basic_Yfunc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : Basic_Func(_sgen,_BS)  {}    

    Kabbala Y(const int&);
    Complex Ycalc(const int&,const int&,const int&,const int&,const Complex&,const Complex&);

    template <int,int>
    Complex YT(const int&,const int&,const Complex&,const Complex&);
  };

  class Basic_Zfunc : public virtual Basic_Func {
  public:
    Basic_Zfunc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : Basic_Func(_sgen,_BS)  {}    
    Kabbala Z(const int&,const int&);
    Complex Zcalc(const int&,const int&,const int&,const int&,
		  const int&,const int&,const int&,const int&,
		  const Complex&,const Complex&,const Complex&,const Complex&);

    int Zmassless(const int&,const int&,const int&,const int&,
		  const int&,const int&,const int&,const int&,
		  const Complex&,const Complex&,const Complex&,const Complex&);
    template <int,int,int,int>
    inline Complex ZT(const int&,const int&,
		      const int&,const int&,
		      const Complex&,const Complex&,
		      const Complex&,const Complex&);

    template <int,int,int,int>
    inline Complex ZTM(const int&,const int&,
		       const int&,const int&,
		       const Complex&,const Complex&,
		       const Complex&,const Complex&);
  };

  class Basic_Xfunc : public virtual Basic_Func {
  public:
    Basic_Xfunc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : Basic_Func(_sgen,_BS)  {}    
    Kabbala X(const int&,const int&);
    Kabbala X(const int&,const int&,const int&,
	      const int&,const int&,
 	      const Complex&,const Complex&);
    Complex Xcalc(const int&,const int&,const int&,
		  const int&,const int&,
		  const Complex&,const Complex&);

    template <int,int>
    inline Complex XT(const int&,const int&,const int&,
		      const Complex&,const Complex&);
  };

  class Basic_Mfunc : public virtual Basic_Func {
  public:
    Basic_Mfunc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : Basic_Func(_sgen,_BS)  {}    
    Kabbala M(const int&);
  };

  class Basic_MassTermfunc : public virtual Basic_Func {
  public:
    Basic_MassTermfunc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : Basic_Func(_sgen,_BS)  {}    
    Kabbala MassTerm(int);
    Complex MassTermCalc(int,int);
    Complex MassTermCalc(int,APHYTOOLS::Flavour);
  };

  class Basic_Vfunc : public virtual Basic_Func {
  public:
    Basic_Vfunc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : Basic_Func(_sgen,_BS)  {}    
    Kabbala V(const int&,const int&);
    Complex Vcalc(const int&,const int&);
    Kabbala Vcplx(const int&,const int&);
    Complex Vcplxcalc(const int&,const int&);
  };

  class Basic_VNCfunc : public virtual Basic_Func {
  public:
    Basic_VNCfunc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : Basic_Func(_sgen,_BS)  {}    
    Kabbala VNC(const int&);
    Kabbala VNC(const int&,const int&);
    Complex VNCcalc(const int&,const int&);
  };
  
  class Basic_Pfunc : public virtual Basic_Func {
  public:
    Basic_Pfunc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : Basic_Func(_sgen,_BS)  {}    
    Kabbala P(const int&);
    Complex Propagator(double,APHYTOOLS::Flavour);
    Complex Pcalc(const int&,const int&);
    Complex Pcalc(const APHYTOOLS::Flavour&,const int&);
    Kabbala P(Pfunc*);
  };
}

#include "Basic_Func.icc"
//#include "String_Generator.H"


#endif
