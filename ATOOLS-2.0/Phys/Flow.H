#ifndef Flow_h
#define Flow_h

#include <iostream>
#include <map>
#include <set>

namespace ATOOLS {
  class Particle;
  class Flow {  
    friend std::ostream& operator<<( std::ostream& ostr, const Flow& f );
  private:
    Particle *            m_owner;
    std::map<int,int>   m_code;
    static long qcd_counter;
  public:
    Flow(Particle * _owner=NULL);
    ~Flow();
    const Particle *     Owner() const;
    int                Code(int _index=1);
    void               SetCode(int code_index=1,int code=-1);
    static void        ResetCounter() { qcd_counter = 600; }
  };  
}

#endif
