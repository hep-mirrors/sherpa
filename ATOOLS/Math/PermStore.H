#ifndef ATOOLS_Math_PermStore_H
#define ATOOLS_Math_PermStore_H

#include <vector>
#include "ATOOLS/Math/MyComplex.H"



namespace ATOOLS {


  class PermStore: public std::vector<PermStore*> {
  private:

    Complex amplitude;
    Complex colorstr;
   
    
  public:
    
    // constructor
    PermStore(size_t pnumber); 

    //destructor
    ~PermStore();

    void PutAmp(size_t* perm,Complex amp);
    void PutColor(size_t* perm,Complex color);
    Complex GetAmp(size_t* perm);
    Complex GetColor(size_t* perm);
    
  };

 


class PermStoreFast: public std::vector<PermStoreFast*> {
  private:

    Complex amplitude;
    Complex colorstr;
   
    
  public:
    
    // constructor
    PermStoreFast(size_t pnumber,size_t m_pnumber); 

    //destructor
    ~PermStoreFast();

    void PutAmp(size_t* perm,Complex amp,size_t pnumber);
    void PutColor(size_t* perm,Complex color,size_t pnumber);
    Complex GetAmp(size_t* perm,size_t pnumber);
    Complex GetColor(size_t* perm,size_t pnumber);
    
  };
 



  
}// end of namespace ATOOLS

#endif
