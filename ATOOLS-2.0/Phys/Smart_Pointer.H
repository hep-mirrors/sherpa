#ifndef Smart_Pointer_H
#define Smart_Pointer_H

#ifndef NULL
#define NULL __null
#endif

#include <iostream>

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

    inline void Reset() const
    {
      p_this=NULL;
      p_owner=this;
      p_copy=NULL;
    }
    /*!
      \fn inline void Reset() const
      \brief Resets the Smart_Pointer
    */

    inline void Assign(const Smart_Pointer &ref) const
    { 
      p_owner=&ref;
      p_this=ref.p_this;
      if ((p_copy=ref.p_copy)!=NULL) p_copy->p_owner=this;
      ref.p_copy=this; 
    }
    /*!
      \fn void Assign(const Smart_Pointer &ref) const
      \brief Inserts this into the chain of Smart_Pointers 
      holding ref.p_this. Sets new links between this, ref
      and ref.p_copy.
    */

    inline void Remove() const
    { 
      if (p_owner!=this) {
	if ((p_owner->p_copy=p_copy)!=NULL) p_copy->p_owner=p_owner;
      }
      else {
	if (p_copy!=NULL) p_copy->p_owner=p_copy; 
	else if (p_this!=NULL) delete p_this; 
      } 
      Reset();
    } 
    /*!
      \fn inline void Remove() const
      \brief Removes this from the chain of Smart_Pointers
      holding p_this. Sets new links between p_owner and p_copy.
      Deletes p_this if p_owner is NULL and resets this.
    */

    void RemoveForward() const;  
    /*!
      \fn void RemoveForward() const
      \brief Removes all linked Smart_Pointers from the chain
      holding p_this, starting at p_copy.
    */

    const Smart_Pointer *FindOwner() const; 
    /*!
      \fn const Smart_Pointer &FindOwner() const
      \brief Determines the real owner of p_this.
    */

  public:

    // constructors
    inline Smart_Pointer(): 
      p_this(NULL), p_owner(this), p_copy(NULL) { }
    inline Smart_Pointer(Type *const type): 
      p_this(type), p_owner(this), p_copy(NULL) { }
    inline Smart_Pointer(const Smart_Pointer &ref): 
      p_this(NULL), p_owner(this), p_copy(NULL) { Assign(ref); }

    // destructor
    inline ~Smart_Pointer() { Remove(); }
    /*!
      \fn inline ~Smart_Pointer()
      \brief Removes the Smart_Pointer from its current chain.
      Deletes p_this if no further Smart_Pointer holds it.
    */

    // inline functions
    inline void New() { Remove(); p_this = new Type(); } 
    /*!
      \fn inline void New()
      \brief Removes the Smart_Pointer from its current chain.
      Deletes p_this if no further Smart_Pointer holds it and
      creates a new object using the default constructor.
    */
    inline void Delete() { FindOwner()->RemoveForward(); } 
    /*!
      \fn inline void Delete()
      \brief Removes all Smart_Pointers holding p_this 
      from the chain. Deletes p_this.
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

    inline Type *const operator=(Type *const ref) 
    { 
      if (ref==p_this) return p_this;
      Remove(); 
      return p_this=ref; 
    }
    /*!
      \fn inline Type *const operator=(Type *const ref)
      \brief Removes this from its current chain if ref 
      is not equal to p_this. 
      Assigns p_this to ref and returns p_this.
    */
    inline Type *const operator=(const Smart_Pointer &ref) 
    { 
      if (ref==p_this) return p_this;
      Remove(); 
      Assign(ref); 
      return p_this; 
    }
    /*!
      \fn inline Type *const operator=(Smart_Pointer &ref) 
      \brief Removes this from its current chain if ref.p_this 
      is not equal to p_this. 
      Assigns this to ref and returns p_this.
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
  namespace ATOOLS { \
									\
    template <> const Smart_Pointer<Class_Type> *			\
    ATOOLS::Smart_Pointer<Class_Type>::FindOwner() const			\
    {									\
      if (p_owner!=this) return p_owner->FindOwner();			\
				return this;							\
		}									\
	\
    template <> void Smart_Pointer<Class_Type>::			\
    RemoveForward() const							\
    {									\
			if (p_copy!=NULL) p_copy->RemoveForward();				\
			Remove();								\
		}									\
  									\
    template <> void Smart_Pointer<Class_Type>::			\
    PrintForward(std::ostream &str,const bool all) const			\
    {									\
			if (all) {								\
				str<<"("<<this<<")["<<typeid(p_this).name()<<"]: p_this = "	\
				   <<p_this<<" {\n";						\
				FindOwner()->PrintForward(str,false);				\
			  str<<"}";								\
				return;								\
			}									\
			str<<"   ("<<this<<"): { p_owner = "<<p_owner			\
			   <<", p_copy = "<<p_copy<<" }\n";					\
		  if (p_copy!=NULL) p_copy->PrintForward(str);			\
		} \
	\
	}

#endif

