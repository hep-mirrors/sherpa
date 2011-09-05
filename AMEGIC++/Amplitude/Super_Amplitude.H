#ifndef AMEGIC_Amplitude_Super_Amplitude_H
#define AMEGIC_Amplitude_Super_Amplitude_H

#include "AMEGIC++/Amplitude/Amplitude_Group.H"
#include "AMEGIC++/Amplitude/Single_Amplitude_Base.H"

namespace AMEGIC {
  class Super_Amplitude: public Amplitude_Group, public Single_Amplitude_Base {
    std::string str;

    int     FindNewNumber(int);
    void    ReduceZfuncs(std::string);
    int     NewSigns(std::vector<std::vector<int> > & zsignlists);
    void    SetZfuncSign();
  public:
    Super_Amplitude(int* _b,int _n,Basic_Sfuncs* _BS,
		    ATOOLS::Flavour* _fl,String_Handler* _shand); 
    ~Super_Amplitude();
    void Init(std::string);

    void PrintGraph();
    void FillCoupling(String_Handler * shand) {Amplitude_Group::FillCoupling(shand);}
    void Add(Amplitude_Base * ab)             {Amplitude_Group::Add(ab);}

    Complex Zvalue(String_Handler * sh,int ihel);
    Complex Zvalue(int ihel,int * signlist);
    Complex Zvalue(int ihel);
    void    ClearCalcList();   
    void    KillZList();         
    void    SetNumber(int & n);    
    Zfunc_List* GetZlist();    

    Amplitude_Base* GetAmplitude(const int n) {return (n==amplnumber) ? this : 0;}
    int GetNumber() { return amplnumber; }

    int GetOrderQED();
    int GetOrderQCD();
  };
}
#endif












