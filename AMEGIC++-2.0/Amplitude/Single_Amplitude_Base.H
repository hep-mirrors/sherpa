#ifndef Single_Amplitude_Base_h
#define Single_Amplitude_Base_h

#include "Amplitude_Base.H"

namespace AMEGIC {


  class Single_Amplitude_Base: public virtual Amplitude_Base {
  protected:
    int amplnumber;

    std::list<Zfunc*> zlist;
    std::list<Pfunc*> plist;
    int   N;
    int*  b;

    String_Handler* shand;
    Basic_Sfuncs* BS;
    APHYTOOLS::Flavour*      fl;

    virtual Kabbala Single_Zvalue(Zfunc*,vector<int>*,vector<int>*,int last=0);
    Kabbala Single_ZGroup_value(Zfunc* z, vector<int>* iz, vector<int>* iargs,int last=0);
    Kabbala Single_ZvalueTensor(Zfunc* z,vector<int>* iz, vector<int>* iargs,int);
    virtual int    Fill_Args(Zfunc* z,Argument* args,vector<int>*,vector<int>*);
    void    SetLoopVar(vector<int>&,vector<vector<int> >&);

    Kabbala Single_Mass_Terms(int iz,int iarg);
    APHYTOOLS::Flavour* Get_Pflav(int pn);
    Kabbala Get_Prop(int,Argument*,char);
    Complex Zvalue_sum(int,int*);
  public:
    int sign;
    Single_Amplitude_Base(int* _b,int _n, Basic_Sfuncs* _BS,
			  APHYTOOLS::Flavour* _fl,
			  String_Handler* _shand) 
      : b(_b), N(_n), shand(_shand), BS(_BS), fl(_fl) {}

    Complex Zvalue(int,int*);
    Complex Zvalue(String_Handler * sh,int ihel);

    inline  std::list<Zfunc*>* GetZlist() {return &zlist;}
    inline  std::list<Pfunc*>* GetPlist() {return &plist;}
    inline  int GetSign() {return sign;}
    inline  void SetSign(int s) {sign=s;}

    void ClearCalcList();
    void PrintGraph(); 
    void SetNumber(int& i) {amplnumber = i;i++;}
    Amplitude_Base* GetAmplitude(const int n) {return (n==amplnumber) ? this : 0;}
  };
}
#endif












