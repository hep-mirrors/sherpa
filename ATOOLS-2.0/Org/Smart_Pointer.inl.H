#ifndef Smart_Pointer_inl_H
#define Smart_Pointer_inl_H

namespace ATOOLS {

  template <class Class_Type> 
  inline void Smart_Pointer<Class_Type>::Reset() 
  { 
    p_this=NULL; 
    p_copy=NULL; 
    p_owner=this;
  }
  
  template <class Class_Type> 
  inline void Smart_Pointer<Class_Type>::New()   
  { 
    Reset(); 
    p_this = new Type(); 
  }
  
  template <class Class_Type> 
  inline void Smart_Pointer<Class_Type>::Delete() 
  { 
    if (p_owner==this) { 
      if (p_copy==NULL) delete p_this; 
      else p_copy.p_owner=p_copy; 
    } 
    else {
      p_owner->p_copy=p_copy;
    }
    Reset(); 
  } 
  
  template <class Class_Type> 
  inline void Smart_Pointer<Class_Type>::Assign(const Smart_Pointer &pointer) 
  { 
    p_this=pointer.p_this; 
    p_copy=pointer.p_copy;
    if (pointer.p_copy!=NULL) pointer.p_copy.p_owner=this;
    pointer.p_copy=this; 
    m_owned=false; 
  }
  
}// end of namespace ATOOLS

#endif
