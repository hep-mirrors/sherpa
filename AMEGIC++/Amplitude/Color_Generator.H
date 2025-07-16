#ifndef AMEGIC_Amplitude_Color_Generator_H
#define AMEGIC_Amplitude_Color_Generator_H

#include "AMEGIC++/Main/Point.H"

namespace AMEGIC {
  struct Cfunc {
    int type; // Number of args, multi gluons etc.
    int arg[3];
    Cfunc* Next;
  };
  
  class Color_Generator {
    MODEL::Color_Function* CFlist;
    MODEL::Color_Function* CCFlist;
    Cfunc* Clist;
    std::string stringchain;

    void   Replace(std::string&,std::string,std::string);
    std::string String(int);
    void   SearchnReplace(int,std::string&,char&,char&);   
    void   PropReplace(Point*,std::string&,char&,char&);
    void   Stringlist(Point*);
    
  public:
    Color_Generator() {Clist=NULL;CFlist = 0;CCFlist=0;stringchain=std::string("");}
    void Convert(Point*);
    void Kill(); 
    void C2string(int,std::string&,std::string&);
    Cfunc* Get() {return Clist;}
    void ColorString(int,Point*,std::string&,std::string&);

    void FillString(int,MODEL::Color_Function*,int&); 

    void CFKill();
    void CFConvert(int,int&,Point* p);
    void CFBuildString(int);
    MODEL::Color_Function* Get_CF() {return CFlist;}
    MODEL::Color_Function* Get_CCF(){return CCFlist;}
    //former private
    std::string CF2String(MODEL::Color_Function*);  
  };
}
#endif





