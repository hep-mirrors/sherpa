#ifndef Zfunc_h
#define Zfunc_h

namespace AMEGIC {
  struct Zfunc;
}

#include "Zfunc_Calc.H"

namespace AMEGIC {

  struct SpinorDirection {
    int from,to;
    SpinorDirection* Next;
  };

  class Zfunc {
  public:

    int       type; // Number of args, multi gluons etc.
    int       narg,ncoupl;
    int*      arg;
    Complex*  coupl;
    // massive propagators
    int       pn;
    Argument* psnew;
    //Lorentzfunction sign
    int       sign;
    std::string       str;
    std::vector<Kabbala>   value;
    std::vector<Argument*> calclist;
    Zfunc*            Equal;    
    Zfunc_Calc*       calculator;

    Zfunc() {Equal = this;narg = 0;ncoupl = 0;pn = 0;}
    Zfunc(const Zfunc& z) {
      narg = 0;ncoupl = 0;pn = 0;
      *this = z;
    }
    ~Zfunc() {
      if (narg!=0)   delete[] arg;
      if (ncoupl!=0) delete[] coupl;
      if (pn!=0)     delete[] psnew;
      for (int i=0;i<calclist.size();i++) delete[] calclist[i];
    }
    	
    Zfunc& operator=(const Zfunc& z) {
      if (this!=&z) {
	if (narg!=0)   delete[] arg;
	if (ncoupl!=0) delete[] coupl;
	if (pn!=0)     delete[] psnew;	
	
	type   = z.type;
	narg   = z.narg;
	ncoupl = z.ncoupl;
	pn     = z.pn;

	if (narg>0) {
	  arg = new int[narg];
	  for (int i=0;i<narg;i++) arg[i] = z.arg[i];
	}
	if (ncoupl>0) {
	  coupl = new Complex[ncoupl];
	  for (int i=0;i<ncoupl;i++) coupl[i] = z.coupl[i];
	}
	if (pn>0) {
	  psnew = new Argument[pn];
	  for (int i=0;i<pn;i++) psnew[i] = z.psnew[i];
	}
	
	sign  = z.sign;
	str   = z.str;
	value = z.value;
	
	for (int i=0;i<calclist.size();i++) delete[] calclist[i];
	calclist.clear();
	for (int i=0;i<z.calclist.size();i++) {
	  Argument* newargs = new Argument[2*narg];
	  for (int j=0;j<2*narg;j++) newargs[j] = (z.calclist[i])[j];
	  calclist.push_back(newargs);
	}

	Equal      = z.Equal;
	calculator = z.calculator;
      }
    }
    //virtual part
    virtual Zfunc* operator[](const int i) {return this;} 

    virtual void Print();
    virtual int  GetSize() {return 1;}

    virtual void SetSign(const int,const int) 
      {std::cerr<<"Called virtual function Zfunc::AddSign()!!!"<<std::endl;}
    virtual int  GetSign(const int i) {   
      std::cerr<<"Called virtual function Zfunc::GetSign()!!!"<<std::endl;
      return 1;
    }
  };

  class Zfunc_Group: public Zfunc {
  public:
    Zfunc_Group() : Zfunc() {}
    Zfunc_Group(const Zfunc& z) : Zfunc(z) {}
    
    std::vector<Zfunc*> zlist;
    Zfunc* operator[](const int i) {return zlist[i];} 
    void Print();
    int  GetSize() {return zlist.size();}
    
    std::vector<int>    zsign;
    void SetSign(const int pos,const int _sign) {zsign[pos] = _sign;}
    int  GetSign(const int i)                   {return zsign[i];}
  };
}

#endif
