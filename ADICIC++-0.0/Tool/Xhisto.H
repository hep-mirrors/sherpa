//bof
//Version: 3 ADICIC++-0.0/2005/07/28

//An additional histogram class.



#ifndef _Xhisto_H_
#define _Xhisto_H_ _Xhisto_H_


#include <list>
#include <string>
#include "Sudakov_Utils.hpp"


//#define XHISTO_OUTPUT XHISTO_OUTPUT
//#undef  XHISTO_OUTPUT





//namespace ... {



  class Xhisto {
  private:
    size_t m_i;
    size_t m_stop;
    std::list<ADICIC::Multidouble> l_sets;
  public:
    Xhisto(size_t stop=1000000000);
    ~Xhisto() {}
    void Reset();
    void Insert(const ADICIC::Multidouble&);
    void Output(const std::string&) const;
    const size_t Entries() const;
  };



//}    //eo namespace ...





#endif    //eo _Xhisto_H_



//eof
