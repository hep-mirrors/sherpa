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
    ATOOLS::Flavour fl;


    Pfunc() {argnum = 0;}
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

	momnum = p.momnum;
	value  = p.value;
	on     = p.on;
	haspol = p.haspol;
	fl     = p.fl;
      }
      return *this;
    }

    std::ostream& operator<<(std::ostream& os) {
      os<<" mn:"<<momnum<<"   argn:"<<argnum<<"   on:"<<on<<"   hp:"<<haspol<<"   fl:"<<fl<<std::endl;
      os<<" arg(";
      for (int i=0; i<argnum; ++i) 
	os<<" "<<arg[i];
      os<<" )"<<std::endl;
      return os;
    }
  };
  typedef std::vector<Pfunc*> Pfunc_List;
  typedef Pfunc_List::iterator Pfunc_Iterator;

}

#endif
