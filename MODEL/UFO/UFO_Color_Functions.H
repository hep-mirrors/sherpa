#ifndef MODEL__UFO__UFO_Color_Functions_H
#define MODEL__UFO__UFO_Color_Functions_H

#include "MODEL/Main/Color_Function.H"

using namespace MODEL;

namespace UFO{

  class UFO_Color_Function : public Color_Function {
  public:
    UFO_Color_Function(cf::code _type, 
		      int _partarg0 = -1, int _partarg1 = -1, int _partarg2 = -1, 
		      Color_Function* n=NULL) :
      Color_Function(_type,
		     _partarg0,_partarg1,_partarg2, 
		     n) {};
    
    UFO_Color_Function operator *(const UFO_Color_Function & other) const;

    void SetParticleArg(int a, int b,int c=-1);
    void SetStringArg(char a, char b, char c='?');

    void SetNext(UFO_Color_Function* next);

  };

  class UFO_CF_1 : public UFO_Color_Function {
  public:
    UFO_CF_1() : UFO_Color_Function(cf::None) {
    };
  };

  class UFO_CF_Identity : public UFO_Color_Function {
  public:
    UFO_CF_Identity(const int& a, const int& b) : UFO_Color_Function(cf::D,a,b) {}
  };

  class UFO_CF_T : public UFO_Color_Function {
  public:
    UFO_CF_T(const int a, const int b, const int c) : UFO_Color_Function(cf::T,a,b,c) {}
  };

  class UFO_CF_f : public UFO_Color_Function {
  public:
    UFO_CF_f(const int a, const int b, const int c) : UFO_Color_Function(cf::F,a,b,c) {}
  };

}

#endif
