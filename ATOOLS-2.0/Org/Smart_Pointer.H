#ifndef Smart_Pointer_H
#define Smart_Pointer_H

#ifndef NULL
#define NULL __null
#endif

namespace ATOOLS {

  template <typename Class_Type> 
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
      if ((p_copy=ref.p_copy)!=NULL) p_copy->p_owner=this;
      ref.p_copy=this; 
    }

    inline void Remove()
    { 
      if (p_owner!=this) p_owner->p_copy=p_copy;
      else {
	if (p_copy!=NULL) p_copy->p_owner=p_copy; 
	else if (p_this!=NULL) delete p_this; 
      } 
      Reset();
    } 

    void RemoveBackward(); 
    void RemoveForward();  

    inline Smart_Pointer(const Smart_Pointer &ref) { }

  public:

    // constructors
    inline Smart_Pointer(): 
      p_this(NULL), p_owner(this), p_copy(NULL) { }
    inline Smart_Pointer(Type *const type): 
      p_this(type), p_owner(this), p_copy(NULL) { }
    inline Smart_Pointer(Smart_Pointer &ref): 
      p_this(NULL), p_owner(this), p_copy(NULL) { Assign(ref); }

    // destructor
    inline ~Smart_Pointer() { Remove(); }

    // inline functions
    inline void New()    { Remove(); p_this = new Type();     } 
    inline void Delete() { RemoveForward(); RemoveBackward(); } 

    inline Type &operator*() const        { return *p_this; }
    inline Type *const operator->() const { return p_this;  }
    
    inline Type *const operator=(Type *const ref) 
    { Remove(); return p_this=ref; }
    inline Type *const operator=(Smart_Pointer &ref) 
    { Remove(); Assign(ref); return p_this; }
    
    inline bool operator==(const Smart_Pointer &ref) const 
    { return ref.p_this==p_this; }
    inline bool operator!=(const Smart_Pointer &ref) const 
    { return !(*this==ref); }
    inline bool operator<(const Smart_Pointer &ref) const  
    { return ref.p_this<p_this;  }
    inline bool operator>=(const Smart_Pointer &ref) const 
    { return !(*this<ref);  }
    inline bool operator>(const Smart_Pointer &ref) const  
    { return ref.p_this>p_this;  }
    inline bool operator<=(const Smart_Pointer &ref) const 
    { return !(*this>ref);  }

    inline bool operator==(const Type *ref) const { return ref==p_this;    }
    inline bool operator!=(const Type *ref) const { return !(p_this==ref); }
    inline bool operator<(const Type *ref) const  { return ref<p_this;     }
    inline bool operator>=(const Type *ref) const { return !(p_this<ref);  }
    inline bool operator>(const Type *ref) const  { return ref>p_this;     }
    inline bool operator<=(const Type *ref) const { return !(p_this>ref);  }

  };// end of class Smart_Pointer

  /*!
    \class Smart_Pointer
    \brief Provides smart pointer handling

    This class provides smart pointer handling. Once a pointer to an 
    object of type Class_Type has been passed to a corresponding
    Smart_Pointer object it will automatically be deleted once the last 
    Smart_Pointer object assigned with it is destroyed. It will also be
    deleted if the Delete() method is called for one of the assigned
    Smart_Pointers. In this case all other Smart_Pointers are updated 
    accordingly, i.e. their references are set to the NULL pointer.

    The appropriate way to generate Smart_Pointers in the code 
    is as follows
    \code
       using namespace ATOOLS;
       ...
       void Example_Function()
       {
         ...
	 Example_Class *pointer = new Example_Class();
	 // call of constructor through #define directive 
	 SP(Example_Class) smart_pointer_1;
	 // pointer assignment
	 smart_pointer_1=pointer;
	 ...
	 // explicit call of constructor
	 Smart_Pointer<Example_Class> smart_pointer_2;
	 // smart pointer assignment
	 smart_pointer_2=smart_pointer_1;
         ...
       }
    \endcode

    To avoid inlining problems, each class which is desired to be
    used with the Delete() method must have the 
    INSTANTIATE_SMART_POINTER(Class_Name) statement in its C-file
    in order to instantiate the appropriate delete methods.
  */

#ifdef _CPP_IOSTREAM
  template <class Type> std::ostream &
  operator<<(std::ostream &str,const Smart_Pointer<Type> &p)
  { return str<<"SP("<<p<<")"; }
#endif

}// end of namespace ATOOLS

#define SP(Class_Type)				\
  ATOOLS::Smart_Pointer<Class_Type> 

#define INSTANTIATE_SMART_POINTER(Class_Type)			\
								\
  template <>							\
  void ATOOLS::Smart_Pointer<Class_Type>::RemoveBackward()	\
  {								\
    Remove();							\
    if (p_owner!=NULL && p_owner!=this)				\
      p_owner->RemoveBackward();				\
  }								\
								\
  template <>							\
  void ATOOLS::Smart_Pointer<Class_Type>::RemoveForward()	\
  {								\
    if (p_copy!=NULL) p_copy->RemoveForward();			\
    Remove();							\
  }								

#endif
