#ifndef Single_Amplitude_Base_h
#define Single_Amplitude_Base_h

#include "Amplitude_Base.H"

namespace AMEGIC {

  struct Block_Info {
    Zfunc* Z;
    int index;
    int l,r;
  };

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

    virtual Kabbala Single_Zvalue(Argument*,Zfunc*);
    virtual void    Fill_Args(Zfunc* z,Argument* args,int* signlist,
			      std::vector<int>* iz=0,std::vector<int>* ii=0,
			      std::vector<std::vector<int> >* iargs=0);
    void    SetLoopVar(std::vector<int>&,std::vector<std::vector<int> >&);
    Kabbala Mass_Terms(std::vector<int>& iz,std::vector<int>& ii,std::vector<std::vector<int> >& iargs);

    Kabbala Single_Mass_Terms_new(int iz,int iarg);
    Kabbala Generate_Block(std::vector<Block_Info>& BlockList,int bi,
			   std::vector<int>& iz,std::vector<int>& ii,
			   std::vector<std::vector<int> >& iargs,
			   int* signlist);
    Complex Zvalue_new_sum(int,int*);
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












