#ifndef AMEGIC_String_String_Handler_H
#define AMEGIC_String_String_Handler_H

#include "ATOOLS/Math/Kabbala.H"
#include "ATOOLS/Math/MyComplex.H"
#include "AMEGIC++/String/String_Generator.H"
#include "AMEGIC++/String/String_Tree.H"
#include "AMEGIC++/String/Values.H"
#include "AMEGIC++/Main/Helicity.H"

namespace AMEGIC {
  class String_Handler {
    int                           gen_str,maxgraph,maxhel;
    Virtual_String_Generator    * sgen;
    int                           own_sgen,own_sk;
    int                           working;
    sknot                     *** sk;
    std::string                ** stringsk;
    String_Tree                   stree,sthelp;
    Values                      * val;
    std::string                   path;
    
    int m_mode;

    void                       Z_Kill();
    void                       KillStrings();
    sknot*                     MakeSknotAFill(std::string &);
  public:
    String_Handler(const int &,Basic_Sfuncs *,std::map<std::string,Complex>*);
    String_Handler(Virtual_String_Generator *); 
    String_Handler(Virtual_String_Generator *,sknot ***); 
    ~String_Handler();  
    bool                      SearchValues(const int,std::string &,Basic_Sfuncs *);
    void                      Initialize(const int &,const int &);
    void                      Set_String(const int &,const int &,const std::string &);
    void                      Calculate();
    void                      Complete(Helicity *);
    void                      Output(Helicity *);
    void                      Output(Helicity *,std::string);
    sknot ***                 GetSKnots()     { return sk; }
    Complex                   Zvalue(int,int); 
    Values *                  Set_Values(std::string &,Basic_Sfuncs *);
    
    Virtual_String_Generator* Get_Generator() { return sgen; }
    int                       Is_String()     { return working; }
    bool                      IsLibrary()     { return (val!=0); }

    inline int                NumberOfCouplings() { return sgen->NumberOfCouplings(); }
    inline int                NumberOfZfuncs()    { return sgen->NumberOfzxl(); }

    int                       LibMode() { return m_mode; }
  };
}

#endif
