#ifndef Smart_Pointer_H
#define Smart_Pointer_H

namespace ATOOLS {

  template <class Class_Type> 
  class Smart_Pointer {
  public:

    typedef Class_Type Type;
    
  private:
    
    Type *p_this;

    Smart_Pointer *p_owner, *p_copy;

    inline void Reset()
    {
      p_this=NULL;
      p_owner=this;
      p_copy=NULL;
    }

    inline void Assign(Smart_Pointer &ref)
    { 
      p_owner=&ref;
      p_this=ref.p_this; 
      p_copy=ref.p_copy;
      ref.p_copy=this; 
      if (p_copy!=NULL) p_copy->p_owner=this;
    }

  public:

    // constructors
    explicit inline Smart_Pointer(): 
      p_this(NULL), p_owner(this), p_copy(NULL) { }
    explicit inline Smart_Pointer(Type *const type): 
      p_this(type), p_owner(this), p_copy(NULL) { }
    explicit inline Smart_Pointer(Smart_Pointer &ref): 
      p_this(NULL), p_owner(this), p_copy(NULL) { Assign(ref); }

    // destructor
    inline ~Smart_Pointer() { Delete(); }

    // inline functions
    inline void New()
    { 
      Delete();
      p_this = new Type();
    } 

    inline void Delete()
    { 
      if (p_owner!=this) p_owner->p_copy=p_copy;
      else {
	if (p_copy==NULL && p_this!=NULL) delete p_this; 
	else p_copy->p_owner=p_copy; 
      } 
      Reset();
    } 

    inline Type *const operator*() const  { return *p_this; }
    inline Type *const operator->() const { return *p_this; }
    
    inline Type *const operator=(Type *const ref) 
    { Delete(); return p_this=ref; }
    inline Type *const operator=(Smart_Pointer &ref) 
    { Delete(); Assign(ref); return p_this; }
    
    inline bool operator==(const Smart_Pointer &ref) const 
    { return (ref.p_this==p_this); }
    inline bool operator!=(const Smart_Pointer &ref) const 
    { return !(*this==ref); }
    inline bool operator<(const Smart_Pointer &ref) const  
    { return (ref.p_this<p_this);  }
    inline bool operator<=(const Smart_Pointer &ref) const 
    { return !(*this>ref);  }
    inline bool operator>(const Smart_Pointer &ref) const  
    { return (ref.p_this>p_this);  }
    inline bool operator>=(const Smart_Pointer &ref) const 
    { return !(*this<ref);  }

  };// end of class Smart_Pointer

}// end of namespace ATOOLS

#define SP(Class_Type) \
  ATOOLS::Smart_Pointer<Class_Type> 

#endif
