#ifndef Pfunc_h
#define Pfunc_h

#include "Flavour.H"
#include <list>

namespace AMEGIC {
  class Pfunc {
  public:
    int*    arg;
    int     argnum;
    int     momnum;
    Complex value;
    int     on;
    int     haspol;
    int     zerowidth;
    ATOOLS::Flavour fl;


    Pfunc() {argnum = 0;}
    Pfunc(int a) {
      argnum = a;
      arg = new int[argnum];
      haspol = 0;
      zerowidth = 0;
    }
    Pfunc(const Pfunc& p) {
      argnum = 0;
      *this = p;
    }
    ~Pfunc() {
      if (argnum>0) delete[] arg;
    }
    
    Pfunc& operator=(const Pfunc& p) {
      if (this!=&p) {
	if (argnum>0) delete[] arg;
	argnum = p.argnum;
	if (argnum>0) {
	  arg = new int[argnum]; 
	  for (int i=0;i<argnum;i++) arg[i] = p.arg[i];
	} 

	momnum    = p.momnum;
	value     = p.value;
	on        = p.on;
	haspol    = p.haspol;
	zerowidth = p.zerowidth;
	fl        = p.fl;
      }
      return *this;
    }

    friend std::ostream& operator<<(std::ostream& os, Pfunc& pf) {
      os<<" mn:"<<pf.momnum<<"   argn:"<<pf.argnum<<"   on:"<<pf.on<<"   hp:"<<pf.haspol
	<<"   fl:"<<pf.fl<<std::endl;
      os<<" arg(";
      for (int i=0; i<pf.argnum; ++i) 
	os<<" "<<pf.arg[i];
      return os<<" )"<<std::endl;
    }

  };
  typedef std::vector<Pfunc*> Pfunc_List;
  typedef Pfunc_List::iterator Pfunc_Iterator;

}

#endif
