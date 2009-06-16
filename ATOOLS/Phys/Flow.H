#ifndef ATOOLS_Phys_Flow_H
#define ATOOLS_Phys_Flow_H

#include <iostream>
#include <map>
#include <set>

namespace ATOOLS {

  class Particle;

  class Flow {  
  private:

    Particle *p_owner;

    std::map<unsigned int,unsigned int> m_code;

    static unsigned int s_qcd_counter;

    friend std::ostream& operator<<(std::ostream &ostr,const Flow &flow);

  public:

    // constructor
    Flow(Particle *const owner=NULL);
    Flow(const Flow &flow);

    // destructor
    ~Flow();
    
    // member functions
    void SetCode(const Flow &flow);
    void SetCode(const unsigned int code_index=1,const int code=-1);

    unsigned int Code(const unsigned int index=1) const;
    int          Index(const unsigned int code) const;

    // inline functions
    inline static void ResetCounter() { s_qcd_counter=600; }

    inline static unsigned int Counter() { return ++s_qcd_counter; }    

    inline const Particle *Owner() const { return p_owner; }

  };// end of class Flow  

}// end of namespace ATOOLS

#endif
