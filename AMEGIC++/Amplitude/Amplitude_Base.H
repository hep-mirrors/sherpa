#ifndef AMEGIC_Amplitude_Amplitude_Base_H
#define AMEGIC_Amplitude_Amplitude_Base_H

#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Phys/Flavour.H"
#include "AMEGIC++/Amplitude/Zfunc.H"

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
    virtual Complex Zvalue(int) = 0;
    virtual Complex Zvalue(int,int*) = 0;
    virtual Complex Zvalue(String_Handler * sh,int ihel) =0;
    virtual void    ClearCalcList() = 0;
    virtual void    KillZList() = 0;
    virtual void    FillCoupling(String_Handler*) {};
    virtual void SetNumber(int&) = 0;
    virtual Amplitude_Base* operator[](const int i) { return this; }
    virtual Amplitude_Base* GetAmplitude(const int) = 0;
    virtual int  GetNumber() { return 0; }
    virtual bool IsGroup() { return false; }

    virtual void SetStringOn();  
    virtual void SetStringOff(); 
    virtual Point*  GetPointlist(); 
    virtual void Add(Amplitude_Base*, int sign=1);
    virtual int Size();

    virtual Zfunc_List* GetZlist();
    virtual Pfunc_List* GetPlist(); 

    virtual int GetSign();
    virtual void SetSign(int);

    virtual void BuildGlobalString(int*,int,Basic_Sfuncs*,ATOOLS::Flavour*,String_Handler*);

    virtual void DefineOrder(const std::vector<int> &o);
    virtual const std::vector<int> &GetOrder();
  };
}
#endif












