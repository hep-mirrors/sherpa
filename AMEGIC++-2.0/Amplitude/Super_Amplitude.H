#ifndef Super_Amplitude_h
#define Super_Amplitude_h

#include "Amplitude_Group.H"
#include "Single_Amplitude_Base.H"

namespace AMEGIC {
  class Super_Amplitude: public Amplitude_Group, public Single_Amplitude_Base {
    std::string str;

    int     FindNewNumber(int);
    void    ReduceZfuncs(std::string);
    int     NewSigns(std::vector<std::vector<int> > & zsignlists);
    void    SetZfuncSign();
  public:
    Super_Amplitude(int* _b,int _n,Basic_Sfuncs* _BS,
		    ATOOLS::Flavour* _fl,String_Handler* _shand) 
      : Single_Amplitude_Base(_b,_n,_BS,_fl,_shand) {groupname = std::string("Super-Amplitude");}
    ~Super_Amplitude() {}
    void Init(std::string);

    void PrintGraph();
    void FillCoupling(String_Handler * shand) {Amplitude_Group::FillCoupling(shand);}
    void Add(Amplitude_Base * ab)             {Amplitude_Group::Add(ab);}

    Complex Zvalue(String_Handler * sh,int ihel);
    Complex Zvalue(int ihel,int * signlist);
    Complex Zvalue(int ihel);
    void    ClearCalcList()                   {Single_Amplitude_Base::ClearCalcList();}
    void    KillZList()                       {Amplitude_Group::KillZList();}
    void    SetNumber(int & n)                {Single_Amplitude_Base::SetNumber(n);}
    Zfunc_List* GetZlist()                    {Single_Amplitude_Base::GetZlist();}

    Amplitude_Base* GetAmplitude(int n)       {return (n==amplnumber) ? this : 0;}
  };
}
#endif












