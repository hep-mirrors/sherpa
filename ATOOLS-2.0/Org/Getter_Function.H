#ifndef Getter_Function_H
#define Getter_Function_H

#include <map>
#include <string>

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

    static bool s_exactmatch;

  protected:

    virtual Object_Type * 
    operator()(const Parameter_Type &parameters) const;

    virtual void PrintInfo(std::ostream &str,const size_t width) const;

  public:
    
    // constructor
    Getter_Function(const std::string &name);

    // destructor
    virtual ~Getter_Function();

    // member functions
    static void PrintGetterInfo(std::ostream &str,const size_t width);

    static Object_Type *GetObject(const std::string &name,
				  const Parameter_Type &parameters);

    // inline functions
    inline static void SetExactMatch(const bool exactmatch)
    { s_exactmatch=exactmatch; }
    
    inline static bool ExactMatch() { return s_exactmatch; }
  
  };// end of class Getter_Function

}// end of namespace ATOOLS

#define DECLARE_GETTER(NAME,TAG,OBJECT,PARAMETER)			\
  									\
  class NAME: public ATOOLS::Getter_Function<OBJECT,PARAMETER> {	\
  private:								\
    									\
    static NAME s_initializer;						\
    									\
  protected:								\
    									\
    void PrintInfo(std::ostream &str,const size_t width) const;		\
    									\
    Object_Type *							\
    operator()(const Parameter_Type &parameters) const;	                \
    									\
  public:								\
    									\
    NAME(): ATOOLS::Getter_Function<OBJECT,PARAMETER>(TAG) {}		\
    									\
  };									\
  NAME NAME::s_initializer				 

#endif
