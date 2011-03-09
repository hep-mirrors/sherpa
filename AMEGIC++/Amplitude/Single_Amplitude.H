#ifndef AMEGIC_Amplitude_Single_Amplitude_H
#define AMEGIC_Amplitude_Single_Amplitude_H

#include "ATOOLS/Phys/Flavour.H"
#include "AMEGIC++/Main/Point.H"
#include "AMEGIC++/Main/Topology.H"
#include "AMEGIC++/Amplitude/Zfunctions/Basic_Sfuncs.H"
#include "AMEGIC++/Amplitude/Zfunc.H"
#include "AMEGIC++/String/String_Handler.H"
#include "MODEL/Interaction_Models/Color_Function.H"
#include "AMEGIC++/Amplitude/Single_Amplitude_Base.H"

namespace AMEGIC {
  class Process_Tags;
  class Single_Amplitude: public Single_Amplitude_Base {
    
    Point* Pointlist;
    MODEL::Color_Function* CFlist;
    MODEL::Color_Function* CCFlist;
    SpinorDirection* spind;
    int icoul;
  public:
    std::string CFColstring;
    std::string CFColstringC;
    void   Set_Coulomb() {icoul = 1;}
    int    Get_Coulomb() {return icoul;}  
    Single_Amplitude* Next;
    Single_Amplitude(Point*,int,int,int*,int,int,Topology*,Basic_Sfuncs*,ATOOLS::Flavour*,String_Handler*);
    Single_Amplitude(Point*,int*,int,int,Topology*,Basic_Sfuncs*,ATOOLS::Flavour*,String_Handler*);
    Single_Amplitude(int*,int,Process_Tags*,Single_Amplitude**,Basic_Sfuncs*,ATOOLS::Flavour*,String_Handler*);
    ~Single_Amplitude(); 
    void PrintGraph(); 
    void Zprojecting(ATOOLS::Flavour*,int,bool);
    void FillCoupling(String_Handler*);
    void MPolconvert(int,int);
    void Prop_Replace(ATOOLS::Flavour,int,int,int);
    Point*  GetPointlist();
    inline  MODEL::Color_Function* Get_CFlist()  {return CFlist;}
    inline  MODEL::Color_Function* Get_CCFlist() {return CCFlist;}

    void    AddSpinorDirection(const int&,const int&);    
    
    void    SetOrderQED();
    void    SetOrderQCD();

    int     GetOrderQED();
    int     GetOrderQCD();
    
  };
}
#endif












