#ifndef Smart_Pointer_H
#define Smart_Pointer_H

namespace ATOOLS {

  template <class Class_Type> 
  class Smart_Pointer {
  public:

    typedef Class_Type Type;
    
  private:
    
    Type *p_this;

    bool m_owned;

    Smart_Pointer *p_copy;

    inline void New();   
    inline void Delete(); 

    inline void Reset(); 
    inline void Assign(const Smart_Pointer &pointer); 

  public:

    // constructors
    inline Smart_Pointer(): 
      p_this(new Type()), m_owned(true), m_copied(false) {}
    inline Smart_Pointer(const Smart_Pointer &pointer): 
      p_this(&pointer), m_owned(false), m_copied(false) {}

    // destructor
    inline ~Smart_Pointer() { Delete(); }

    // inline functions
    inline Type *const operator*() const  { return *p_this;  }
    inline Type *const operator->() const { return *p_this;  }
    inline Type *const operator-() const  { return p_this;   }
    
    inline Type *const operator=(const Smart_Pointer &ref) { Delete(); Assign(ref); }
    
    inline bool operator==(const Smart_Pointer &ref) const { return (-ref==p_this); }
    inline bool operator<(const Smart_Pointer &ref) const  { return (-ref<p_this);  }
    inline bool operator>(const Smart_Pointer &ref) const  { return (-ref>p_this);  }
    inline bool operator<=(const Smart_Pointer &ref) const { return !(p_this>ref);  }
    inline bool operator>=(const Smart_Pointer &ref) const { return !(p_this<ref);  }

  };// end of class Smart_Pointer

}// end of namespace ATOOLS

#include "Smart_Pointer.inl.H"

#endif
