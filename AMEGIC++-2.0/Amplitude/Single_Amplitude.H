#ifndef Single_Amplitude_h
#define Single_Amplitude_h

#include "Flavour.H"
#include "Point.H"
#include "Topology.H"
#include "Basic_Sfuncs.H"
#include "Zfunc.H"
#include "String_Handler.H"
#include "Color_Function.H"

#include "Single_Amplitude_Base.H"

namespace AMEGIC {
  class Single_Amplitude: public Single_Amplitude_Base {

    Point* Pointlist;
    Color_Function* CFlist;
    Color_Function* CCFlist;
    SpinorDirection* spind;
    int icoul;
  public:
    std::string Colstring;
    std::string ColstringC;
    std::string CFColstring;
    std::string CFColstringC;
    void   Set_Coulomb() {icoul = 1;}
    int    Get_Coulomb() {return icoul;}  
    Single_Amplitude* Next;
    Single_Amplitude(Point*,int,int,int*,int,int,Topology*,Basic_Sfuncs*,APHYTOOLS::Flavour*,String_Handler*);
    Single_Amplitude(Point*,int*,int,int,Topology*,Basic_Sfuncs*,APHYTOOLS::Flavour*,String_Handler*);
    ~Single_Amplitude(); 
    void PrintGraph(); 
    void Zprojecting(APHYTOOLS::Flavour*,int);
    void FillCoupling(String_Handler*);
    void MPolconvert(int,int);
    void Prop_Replace(APHYTOOLS::Flavour,int,int,int);
    Point*  GetPointlist() {return Pointlist;} 
    inline  Color_Function* Get_CFlist()  {return CFlist;}
    inline  Color_Function* Get_CCFlist() {return CCFlist;}

    void    AddSpinorDirection(const int&,const int&);    
  };
}
#endif












