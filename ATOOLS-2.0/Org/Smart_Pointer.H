#ifndef Smart_Pointer_H
#define Smart_Pointer_H

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

    void SetForward(Type *const ref) const;  
    /*!
      \fn void SetForward(Type *const ref) const
      \brief Sets p_this of all linked Smart_Pointers to ref.
    */

    const Smart_Pointer *FindOwner() const; 
    /*!
      \fn const Smart_Pointer &FindOwner() const
      \brief Determines the first owner of p_this.
    */
    const Smart_Pointer *FindCopy() const; 
    /*!
      \fn const Smart_Pointer &FindCopy() const
      \brief Determines the last copy of p_this.
    */

  public:

    // constructors
    inline Smart_Pointer(): 
      p_this(NULL), p_owner(NULL), p_copy(NULL) { }
    inline Smart_Pointer(Type *const type): 
      p_this(type), p_owner(NULL), p_copy(NULL) { }
    Smart_Pointer(const Smart_Pointer &ref); 

    // destructor
    ~Smart_Pointer();
    /*!
      \fn inline ~Smart_Pointer()
      \brief Removes the Smart_Pointer from its current chain.
      Deletes p_this if no further Smart_Pointer holds it.
    */

    // inline functions
    void Delete() const;
    /*!
      \fn inline void Delete()
      \brief Deletes p_this and updates all connected Smart_Pointers.
    */
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

    const Smart_Pointer &operator=(Type *const ref);
    /*!
      \fn inline Type *const operator=(Type *const ref)
      \brief If not equal to ref, deletes p_this and sets it to ref.
      Updates all connected Smart_Pointers.
    */
    const Smart_Pointer &operator=(const Smart_Pointer &ref);
    /*!
      \fn inline Type *const operator=(Smart_Pointer &ref) 
      \brief If not equal to ref.p_this, deletes p_this 
      and sets it to ref.p_this.
      Concatenates the two chains of connected Smart_Pointers.
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
    Idem applies to the output operator << and 
    the Print(..) method, respectively,
  */

  template <class Class_Type> inline std::ostream &
  operator<<(std::ostream &str,const ATOOLS::Smart_Pointer<Class_Type> &p)
  { 
    p.PrintForward(str,true); 
    return str; 
  }

}// end of namespace ATOOLS

#define SP(Class_Type)				\
  ATOOLS::Smart_Pointer<Class_Type>


#define INSTANTIATE_SMART_POINTER(Class_Type)				\
    template <> ATOOLS::Smart_Pointer<Class_Type>::              	\
    Smart_Pointer(const Smart_Pointer &ref):                            \
      p_this(NULL), p_owner(NULL), p_copy(NULL)                         \
    {                                                                   \
      p_owner=&ref;                                                     \
      p_this=ref.p_this;                                                \
      if ((p_copy=ref.p_copy)!=NULL) p_copy->p_owner=this;              \
      ref.p_copy=this;                                                  \
    }                                                                   \
    template <> ATOOLS::Smart_Pointer<Class_Type>::                     \
    ~Smart_Pointer()                                                    \
    {                                                                   \
      if (p_owner) {                                                    \
	if ((p_owner->p_copy=p_copy)!=NULL) p_copy->p_owner=p_owner;    \
      }                                                                 \
      else {                                                            \
	if (p_copy) p_copy->p_owner=NULL;                               \
	else if (p_this!=NULL) delete p_this;                           \
      }                                                                 \
    }                                                                   \
    template <> const ATOOLS::Smart_Pointer<Class_Type> *		\
    ATOOLS::Smart_Pointer<Class_Type>::FindOwner() const	        \
    {									\
      if (p_owner) return p_owner->FindOwner();			        \
      return this;				                        \
    }							                \
    template <> const ATOOLS::Smart_Pointer<Class_Type> *		\
    ATOOLS::Smart_Pointer<Class_Type>::FindCopy() const	       	        \
    {									\
      if (p_copy) return p_copy->FindCopy();			        \
      return this;				                        \
    }							                \
    template <> void ATOOLS::Smart_Pointer<Class_Type>::		\
    SetForward(Class_Type *const ref) const				\
    {									\
      if (p_copy!=NULL) p_copy->SetForward(ref);	        	\
      p_this=ref;             						\
    }									\
    template <> void ATOOLS::Smart_Pointer<Class_Type>::Delete() const  \
    {                                                                   \
      if (p_this) delete p_this;                                        \
      FindOwner()->SetForward(NULL);                                    \
    }                                                                   \
    template <> const Smart_Pointer<Class_Type> &ATOOLS::               \
    Smart_Pointer<Class_Type>::operator=(Class_Type *const ref)         \
    {                                                                   \
      if (ref==p_this) return *this;                                    \
      if (p_this) delete p_this;                                        \
      FindOwner()->SetForward(ref);                                     \
      return *this;                                                     \
    }                                                                   \
    template <> const Smart_Pointer<Class_Type> &ATOOLS::               \
    Smart_Pointer<Class_Type>::operator=(const Smart_Pointer &ref)      \
    {                                                                   \
      if (ref==p_this) return *this;                                    \
      if (p_this) delete p_this;                                        \
      const Smart_Pointer *own(FindOwner());                            \
      (own->p_copy=ref.FindCopy())->p_owner=own;                        \
      own->SetForward(ref.p_this);                                      \
      return *this;                                                     \
    }                                                                   \
    template <> void ATOOLS::Smart_Pointer<Class_Type>::		\
    PrintForward(std::ostream &str,const bool all) const		\
    {									\
      if (all) {							\
	str<<"("<<this<<")["<<typeid(p_this).name()<<"]: p_this = "	\
	   <<p_this<<" {\n";						\
	FindOwner()->PrintForward(str,false);				\
        str<<"}";							\
	return;								\
      }									\
      str<<"   ("<<this<<"): { p_owner = "<<p_owner			\
         <<", p_copy = "<<p_copy<<" }\n";				\
      if (p_copy!=NULL) p_copy->PrintForward(str);			\
    }                                                                   \


#endif
