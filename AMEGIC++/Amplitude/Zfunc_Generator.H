#ifndef AMEGIC_Amplitude_Zfunc_Generator_H
#define AMEGIC_Amplitude_Zfunc_Generator_H

#include "AMEGIC++/Main/Point.H"
#include <vector>
#include "AMEGIC++/Amplitude/Zfunc.H"
#include "ATOOLS/Org/Smart_Pointer.H"

namespace AMEGIC {
  class Zfunc_Calc;
  struct SpinorDirection;
  class Basic_Sfuncs;
  class Virtual_String_Generator;

  typedef std::vector<SP(Zfunc_Calc)> ZF_Vector;

  class Zfunc_Generator {
    Zfunc_List zlist;
    static ZF_Vector zcalc;
    Basic_Sfuncs* BS;
    void  IsGaugeV(Point*,int&,int&);
    
    int   TryArgs(int,std::vector<MODEL::Lorentz_Function*> &,int,std::vector<MODEL::Lorentz_Function*> &);
    void  LFFill_Zfunc(Zfunc*,std::vector<MODEL::Lorentz_Function*> &,Point*,Point*,Point*);
  
    int   LFDetermine_Zfunc(Zfunc*,Point*,Point*,Point*);

    void  Lorentz_Sequence(Point*,std::vector<MODEL::Lorentz_Function*> &lflist);
    void  LFPrint(const std::vector<MODEL::Lorentz_Function*>&);
    std::string LFEff(const std::string &);

    int Compare(int,const std::vector<MODEL::Lorentz_Function*> &,int*,
		    const std::vector<MODEL::Lorentz_Function*> &,int*);
    void CopyOrder(std::vector<MODEL::Lorentz_Function*> &,std::vector<MODEL::Lorentz_Function*> &);
    void SearchNextProp(int,const std::vector<MODEL::Lorentz_Function*> &,int*,
			    const std::vector<MODEL::Lorentz_Function*> &,int*,
			int,int);
    void SetPropDirection(int,int,const std::vector<MODEL::Lorentz_Function*> &,int*,
			          const std::vector<MODEL::Lorentz_Function*> &,int*);
  public:
    Zfunc_Generator(Basic_Sfuncs* _BS) : BS(_BS) {}
    ~Zfunc_Generator();
    
    void  Set_In(Zfunc*,int,Point*,Point*,Point*);
    void  Set_Out(Zfunc*,int,Point*,Point*);
    void  Set_Tensor(Zfunc*,Point*);
    void  Set_FermionProp(Zfunc*,Point*,Point*);

    void SetArgs(Zfunc*,int*,int*,Point*,Point*,int&);
    void SetScalarArgs(Zfunc*,int&,Point*);

    void BuildZlist(Virtual_String_Generator*,Basic_Sfuncs*,int);
    void LorentzConvert(Point*);
    void MarkCut(Point*,int notcut,bool fromfermion=false,bool cutvectors=false);
    void Convert(Point*);

    void SetDirection(int,SpinorDirection*);
    void Get(Zfunc_List& _zlist) {
      for (Zfunc_Iterator zit=zlist.begin();zit!=zlist.end();++zit) 
	_zlist.push_back(*zit); 
    }
  };
}
#endif








