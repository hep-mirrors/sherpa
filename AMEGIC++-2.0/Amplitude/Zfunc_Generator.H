#ifndef Zfunc_Generator_h
#define Zfunc_Generator_h

#include "Point.H"
#include "Zfunc.H"
#include "Basic_Sfuncs.H"
#include "String_Generator.H"
#include "Zfunc_Calc.H"
#include <vector>

namespace AMEGIC {

  class Zfunc_Generator {
    std::list<Zfunc*> zlist;
    std::vector<Zfunc_Calc*> zcalc;
    Basic_Sfuncs* BS;
    void  Set_In(Zfunc*,int,Point*,Point*,Point*);
    void  Set_Out(Zfunc*,int,Point*,Point*);
    void  IsGaugeV(Point*,int&,int&);
    
    int   TryArgs(int,std::vector<Lorentz_Function> &,int,std::vector<Lorentz_Function> &);
    void  LFFill_Zfunc(Zfunc*,std::vector<Lorentz_Function> &,Point*,Point*,Point*);
  
    int   LFDetermine_Zfunc(Zfunc*,Point*,Point*,Point*);

    void  Lorentz_Sequence(Point*,std::vector<Lorentz_Function> &lflist);
    void  LFPrint(const std::vector<Lorentz_Function>&);
    void  LFPrint(const std::vector<Lorentz_Function*>&);
    lf::code LFEff(lf::code);

    int Compare(int,const std::vector<Lorentz_Function*> &,int*,
		    const std::vector<Lorentz_Function*> &,int*);
    void CopyOrder(std::vector<Lorentz_Function> &,std::vector<Lorentz_Function*> &);
    void SetArgs(Zfunc*,int*,int*,Point*,Point*,int&);
    void SearchNextProp(int,const std::vector<Lorentz_Function*> &,int*,
			    const std::vector<Lorentz_Function*> &,int*,
			int,int);
    void SetPropDirection(int,int,const std::vector<Lorentz_Function*> &,int*,
			          const std::vector<Lorentz_Function*> &,int*);
  public:
    Zfunc_Generator(Basic_Sfuncs* _BS) : BS(_BS) {}

    void BuildZlist(Virtual_String_Generator*,Basic_Sfuncs*);
    void LorentzConvert(Point*);
    void MarkCut(Point*,int notcut,bool fromfermion=false);
    void Convert(Point*);

    void SetDirection(int,SpinorDirection*);
    void Get(std::list<Zfunc*>& _zlist) {
      for (std::list<Zfunc*>::iterator zit=zlist.begin();zit!=zlist.end();++zit) 
	_zlist.push_back(*zit); 
    }
  };
}
#endif








