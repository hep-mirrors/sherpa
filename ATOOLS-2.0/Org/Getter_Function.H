#ifndef Getter_Function_H
#define Getter_Function_H

#include <map>

namespace ATOOLS {

  template <class ObjectType,class ParameterType,class OwnerType>
  class Getter_Function {
  public:

    typedef ObjectType Object_Type;
    typedef ParameterType Parameter_Type;
    typedef OwnerType Owner_Type;
    typedef typename std::map<const std::string,
			      Getter_Function *const> String_Getter_Map;

  private:

    static String_Getter_Map *s_getters;

  protected:

    virtual Object_Type *const 
    operator()(const Parameter_Type &parameters,OwnerType *const owner) const;

    virtual void PrintInfo(std::ostream &str,const size_t width) const;

  public:
    
    // constructor
    Getter_Function(const std::string &name);

    // destructor
    virtual ~Getter_Function();

    // member functions
    static void PrintGetterInfo(std::ostream &str,const size_t width);

    static Object_Type *const GetObject(const std::string &name,
					const Parameter_Type &parameters,
					OwnerType *const owner);
    
  };// end of class Getter_Function

  template <class ObjectType,class ParameterType> class Primitive_Getter_Function:
    public Getter_Function<ObjectType,ParameterType,void> {
  public:

    typedef ObjectType Object_Type;
    typedef ParameterType Parameter_Type;
    typedef typename std::map<const std::string,
			      Primitive_Getter_Function *const> String_Getter_Map;
    
  public:

    // constructor
    Primitive_Getter_Function(const std::string &name);

    // member functions
    static Object_Type *const GetObject(const std::string &name,
					const Parameter_Type &parameters);
    
  };

}// end of namespace ATOOLS

#define DECLARE_GETTER(NAME,TAG,OBJECT,PARAMETER,OWNER)			\
  									\
  class NAME: public ATOOLS::Getter_Function<OBJECT,PARAMETER,OWNER> {	\
  private:								\
    									\
    static NAME s_initializer;						\
    									\
  protected:								\
    									\
    void PrintInfo(std::ostream &str,const size_t width) const;		\
    									\
    Object_Type *const							\
      operator()(const Parameter_Type &parameters,                      \
		 Owner_Type *const owner) const;	                \
    									\
  public:								\
    									\
    NAME(): ATOOLS::Getter_Function<OBJECT,PARAMETER,OWNER>(TAG) {}	\
      									\
  };									\
  NAME NAME::s_initializer				 

#define DECLARE_PRIMITIVE_GETTER(NAME,TAG,OBJECT,PARAMETER)		\
  									\
  class NAME:                                                           \
  public ATOOLS::Primitive_Getter_Function<OBJECT,PARAMETER> {	        \
  private:								\
    									\
    static NAME s_initializer;						\
    									\
  protected:								\
    									\
    void PrintInfo(std::ostream &str,const size_t width) const;		\
    									\
    Object_Type *const							\
      operator()(const Parameter_Type &parameters) const;	        \
    									\
  public:								\
    									\
    NAME(): ATOOLS::Primitive_Getter_Function<OBJECT,PARAMETER>(TAG) {}	\
      									\
  };									\
  NAME NAME::s_initializer				 

#endif
