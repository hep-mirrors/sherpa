#ifndef Primitive_Interpreter_H
#define Primitive_Interpreter_H

#include "Node.H"

#include <string>
#include <set>
#include <map>

namespace ATOOLS {

  class Function {
  private:

    std::string m_tag;

  public:

    // constructor
    inline Function(const std::string &tag): m_tag(tag) {}

    // destructor
    virtual ~Function();

    // member functions
    virtual std::string Evaluate(const std::vector<std::string> &args) const;

    // inline functions
    inline const std::string &Tag() const { return m_tag; }

  };// end of class Function

  inline bool operator<(const Function &f1,const Function &f2)
  { return f1.Tag()<f2.Tag(); }

  class Operator: public Function {
  private:

    size_t m_priority;
    bool   m_binary;
    
  public:

    // constructor
    inline Operator(const std::string &tag,
		    const size_t &priority,const bool binary):
      Function(tag), m_priority(priority), m_binary(binary) {}

    // destructor
    ~Operator();

    // inline functions
    inline const size_t &Priority() const { return m_priority; }
    inline const bool    Binary() const   { return m_binary;   }

  };// end of class Operator

  class Primitive_Interpreter;

  class Interpreter_Function {
  protected:
    
    Primitive_Interpreter *p_interpreter;

  public:
    
    // constructor
    inline Interpreter_Function(Primitive_Interpreter *interpreter):
      p_interpreter(interpreter) {}

    // member functions
    virtual std::string 
    Interprete(const std::string &expr) = 0;

  };// end of class Interpreter_Function

  class Primitive_Interpreter {
  public:

    typedef std::map<std::string,std::string> String_Map;

    typedef std::set<Function*>             Function_Set;
    typedef std::pair<size_t,Operator*>     Operator_Pair;
    typedef std::multimap<size_t,Operator*> Operator_Map;
    typedef std::set<Interpreter_Function*> Interpreter_Set;

  private:

    Function_Set m_functions;
    Operator_Map m_operators;

    Interpreter_Set m_interpreters;

    String_Map m_tags;

    std::string &ReplaceTags(std::string &expr) const;

  public:

    // constructors
    Primitive_Interpreter(const bool standard=true);

    // destructor
    ~Primitive_Interpreter();

    // member functions
    std::string Interprete(const std::string &expr);
    std::string Iterate(const std::string &expr);

    void AddFunction(Function *const f);
    void AddOperator(Operator *const b);

    void AddTag(const std::string &tag,const std::string &value);
    void SetTags(const String_Map &tags);

    // inline functions
    inline const Function_Set &Functions() const { return m_functions; }
    inline const Operator_Map &Operators() const { return m_operators; }

  };// end of class Primitive_Interpreter

}// end of namespace ATOOLS

#define DEFINE_FUNCTION(NAME,TAG)					\
  class NAME: public Function {						\
  public:								\
    inline NAME(): Function(TAG) {}					\
    std::string Evaluate(const std::vector<std::string> &args) const;	\
  };									\
  std::string NAME::Evaluate(const std::vector<std::string> &args) const

#define DEFINE_UNARY_OPERATOR(NAME,TAG,PRIORITY)			\
  class NAME: public Operator {						\
  public:								\
    inline NAME(): Operator(TAG,PRIORITY,false) {}			\
    std::string Evaluate(const std::vector<std::string> &args) const;	\
  };									\
  std::string NAME::Evaluate(const std::vector<std::string> &args) const

#define DEFINE_BINARY_OPERATOR(NAME,TAG,PRIORITY)			\
  class NAME: public Operator {						\
  public:								\
    inline NAME(): Operator(TAG,PRIORITY,true) {}			\
    std::string Evaluate(const std::vector<std::string> &args) const;	\
  };									\
  std::string NAME::Evaluate(const std::vector<std::string> &args) const

#define DEFINE_INTERPRETER_FUNCTION(NAME)				\
  class NAME: public Interpreter_Function {				\
  public:								\
    inline NAME(Primitive_Interpreter *interpreter):			\
      Interpreter_Function(interpreter) {}				\
    std::string Interprete(const std::string &expr);			\
  };									\
  std::string NAME::Interprete(const std::string &expr)

#endif
