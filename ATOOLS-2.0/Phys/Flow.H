#ifndef Flow_h
#define Flow_h

#include <iostream>
#include <map>
#include <set>

namespace ATOOLS {
  class Parton;
  class Flow {  
    friend std::ostream& operator<<( std::ostream& ostr, const Flow& f );
  private:
    Parton *            m_owner;
    std::map<int,int>   m_code;
    static long qcd_counter;
  public:
    Flow(Parton * _owner=NULL);
    ~Flow();
    const Parton *     Owner() const;
    int                Code(int _index=1);
    void               SetCode(int code_index=1,int code=-1);
    static void        ResetCounter() { qcd_counter = 600; }
  };  
}

#endif
