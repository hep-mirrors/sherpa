#ifndef Color_Generator_h
#define Color_Generator_h

#include "Point.H"

namespace AMEGIC {
  struct Cfunc {
    int type; // Number of args, multi gluons etc.
    int arg[3];
    Cfunc* Next;
  };
  
  class Color_Generator {
    Color_Function* CFlist;
    Color_Function* CCFlist;
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

    void FillString(int,Color_Function*,int&); 

    void CFKill();
    void CFConvert(int,int&,Point* p);
    void CFBuildString(int);
    Color_Function* Get_CF() {return CFlist;}
    Color_Function* Get_CCF(){return CCFlist;}
    //former private
    std::string CF2String(Color_Function*);  
  };
}
#endif





