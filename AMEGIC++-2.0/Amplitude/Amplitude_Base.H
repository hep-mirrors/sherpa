#ifndef Amplitude_Base_h
#define Amplitude_Base_h

#include "MyComplex.H"
#include "Flavour.H"
#include "Zfunc.H"

namespace AMEGIC {

  class String_Handler;
  class Point;
  class Basic_Sfuncs;

  class Amplitude_Base {
  protected:
    int buildstring;
  public:
    int on;
    int     topnum, permnum;
    virtual ~Amplitude_Base() {} 
    virtual void    PrintGraph() = 0; 
    virtual Complex Zvalue(int,int*) = 0;
    virtual Complex Zvalue(String_Handler * sh,int ihel) =0;
    virtual void    ClearCalcList() = 0;
    virtual void    FillCoupling(String_Handler*) = 0;
    virtual void SetNumber(int&) = 0;
    virtual Amplitude_Base* GetAmplitude(const int) = 0;

    virtual void SetStringOn()  {buildstring = 1;}
    virtual void SetStringOff() {buildstring = 0;}
    virtual Point*  GetPointlist() 
      {std::cerr<<"Error: Virtual  Amplitude_Base::GetPointlist() called!"<<std::endl;return 0;}

    virtual void Add(Amplitude_Base*, int sign=1)
      {std::cerr<<"Error: Virtual  Amplitude_Base::Add() called!"<<std::endl;}
    virtual int Size()
      {std::cerr<<"Error: Virtual  Amplitude_Base::Size() called!"<<std::endl;return 0;}

    virtual Zfunc_List* GetZlist()
      {std::cerr<<"Error: Virtual  Amplitude_Base::GetZlist() called!"<<std::endl;return 0;}
    virtual Pfunc_List* GetPlist()
      {std::cerr<<"Error: Virtual  Amplitude_Base::GetPlist() called!"<<std::endl;return 0;}

    virtual int GetSign()
      {std::cerr<<"Error: Virtual  Amplitude_Base::GetSign() called!"<<std::endl;return 0;}

    virtual void SetSign(int)
      {std::cerr<<"Error: Virtual  Amplitude_Base::SetSign() called!"<<std::endl;}

    virtual void BuildGlobalString(int*,int,Basic_Sfuncs*,ATOOLS::Flavour*,String_Handler*)
      {std::cerr<<"Error: Virtual  Amplitude_Base::BuildGlobalString() called!"<<std::endl;}
  };
}
#endif












