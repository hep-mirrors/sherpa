#ifndef Getter_Function_H
#define Getter_Function_H

#include <map>

namespace ATOOLS {

  template <class ObjectType,class ParameterType>
  class Getter_Function {
  public:

    typedef ObjectType Object_Type;
    typedef ParameterType Parameter_Type;
    typedef typename std::map<const std::string,
			      Getter_Function *const> String_Getter_Map;

  private:

    static String_Getter_Map *s_getters;

  protected:

    virtual Object_Type *const 
    operator()(const Parameter_Type &parameters) const;

    virtual void PrintInfo(std::ostream &str) const;

  public:
    
    // constructor
    Getter_Function(const std::string &name);

    // destructor
    virtual ~Getter_Function();

    // member functions
    static void PrintGetterInfo(std::ostream &str);

    static Object_Type *const GetObject(const std::string &name,
					const Parameter_Type &parameters);
    
  };// end of class Getter_Function

}// end of namespace ATOOLS

#define DECLARE_GETTER(NAME,TAG,OBJECT,PARAMETER)	 \
  							 \
  class NAME: public Getter_Function<OBJECT,PARAMETER> { \
  private:						 \
    							 \
    static NAME s_initializer;				 \
    							 \
  public:						 \
    							 \
    NAME(): Getter_Function<OBJECT,PARAMETER>(TAG) {}	 \
    							 \
    Object_Type *const					 \
    operator()(const Parameter_Type &parameters) const;	 \
    							 \
  };							 \
  NAME NAME::s_initializer				 

#endif
