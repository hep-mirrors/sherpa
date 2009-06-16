#ifndef ATOOLS_Org_Smart_Pointer_H
#define ATOOLS_Org_Smart_Pointer_H

#ifndef NULL
#define NULL __null
#endif

#include <iostream>
#include <typeinfo>

namespace ATOOLS {

  template <typename Class_Type> 
  class Smart_Pointer {
  public:

    typedef Class_Type Type;
    /*!
      \typedef Type
      \brief Stores the object type
    */
    
  private:
    
    mutable Type *p_this;
    /*!
      \var p_this
      \brief Stores the actual object pointer.
    */

    mutable const Smart_Pointer *p_owner, *p_copy;
    /*!
      \var p_owner
      \brief Stores a pointer to the Smart_Pointer 
      where p_this has been taken from. Is this by default.
    *//*!
      \var p_copy
      \brief Stores a pointer to the first Smart_Pointer 
      having taken a copy of p_this. Is NULL by default.
    */

    void Connect(const Smart_Pointer &ref) const;
    /*!
      \fn inline void Connect(const Smart_Pointer &ref)
      \brief Connects to the chain of Smart_Pointers
      belonging to ref.
    */
    void Deconnect() const;
    /*!
      \fn inline void Deconnect()
      \brief Deconnects from connected Smart_Pointers
      and eventually deletes p_this.
    */

    const Smart_Pointer *FindOwner() const;
    /*!
      \fn const Smart_Pointer &FindOwner() const
      \brief Determines the first owner of p_this.
    */

  public:

    // constructors
    inline Smart_Pointer(): 
      p_this(NULL), p_owner(NULL), p_copy(NULL) { }
    inline Smart_Pointer(Type *const type): 
      p_this(type), p_owner(NULL), p_copy(NULL) { }
    inline Smart_Pointer(const Smart_Pointer &ref): 
      p_this(NULL), p_owner(NULL), p_copy(NULL) { Connect(ref); }

    // destructor
    inline ~Smart_Pointer() { Deconnect(); }
    /*!
      \fn inline ~Smart_Pointer()
      \brief Removes the Smart_Pointer from its current chain.
      Deletes p_this if no further Smart_Pointer holds it.
    */

    // inline functions
    inline Type &operator*() const { return *p_this; }
    /*!
      \fn inline Type &operator*() const
      \brief Returns a reference to p_this.
    */
    inline Type *const operator->() const { return p_this;  }
    /*!
      \fn inline Type *const operator->() const
      \brief Returns a pointer to p_this.
    */
    inline Type *const operator--() const { return p_this;  }
    /*!
      \fn inline Type *const operator--() const
      \brief Returns a pointer to p_this.
    */

    inline const Smart_Pointer &operator=(Type *const ref)
    {
      if (ref==p_this) return *this;
      Deconnect(); p_this=ref; return *this;
    }
    /*!
      \fn inline Type *const operator=(Type *const ref)
      \brief If not equal to ref, sets p_this to ref.
      Deconnects from connected Smart_Pointers 
      and eventually deletes p_this.
    */
    inline const Smart_Pointer &operator=(const Smart_Pointer &ref)
    {
      if (ref==p_this) return *this;
      Deconnect(); Connect(ref); return *this;
    }
    /*!
      \fn inline Type *const operator=(Smart_Pointer &ref) 
      \brief If not equal to ref.p_this, sets p_this to ref.p_this.
      Deconnects from connected Smart_Pointers 
      and eventually deletes p_this.
    */
    
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

    void PrintForward(std::ostream &str,const bool all=false) const;
    
  };// end of class Smart_Pointer

  /*!
    \class Smart_Pointer
    \brief Provides smart pointer handling

    This class provides smart pointer handling. Once a pointer to an 
    object of type Class_Type has been passed to a corresponding
    Smart_Pointer object it will automatically be deleted once the last 
    Smart_Pointer object assigned with it is destroyed. 

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

    To avoid inlining problems, each class which is to be
    used with Smart_Pointers must have the 
    template class SP(Class_Name) statement in its C-file.
  */

  template <class Class_Type> inline std::ostream &
  operator<<(std::ostream &str,const ATOOLS::Smart_Pointer<Class_Type> &p)
  { 
    p.PrintForward(str,true); 
    return str; 
  }

}// end of namespace ATOOLS

#define SP(Class_Type) ATOOLS::Smart_Pointer<Class_Type>

#endif
