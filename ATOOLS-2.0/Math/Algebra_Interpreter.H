#ifndef Algebra_Interpreter_H
#define Algebra_Interpreter_H

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
    Function(const std::string &tag);

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

  class Algebra_Interpreter;

  class Interpreter_Function {
  protected:
    
    Algebra_Interpreter *p_interpreter;

  public:
    
    // constructor
    inline Interpreter_Function(Algebra_Interpreter *interpreter):
      p_interpreter(interpreter) {}

    // member functions
    virtual std::string 
    Interprete(const std::string &expr) = 0;

  };// end of class Interpreter_Function

  class Tag_Replacer {
  public:

    // destructor
    virtual ~Tag_Replacer();

    // member functions
    std::string &KillBlanks(std::string& expr) const;

    virtual std::string ReplaceTags(std::string &expr) const;    

  };// end of class Tag_Replacer

  class Algebra_Interpreter: public Tag_Replacer {
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

    Tag_Replacer *p_replacer;

    String_Map m_tags;

  public:

    // constructors
    Algebra_Interpreter(const bool standard=true);

    // destructor
    virtual ~Algebra_Interpreter();

    // member functions
    std::string Interprete(const std::string &expr);
    std::string Iterate(const std::string &expr);

    void AddFunction(Function *const f);
    void AddOperator(Operator *const b);

    std::string ReplaceTags(std::string &expr) const;    

    void AddTag(const std::string &tag,const std::string &value);
    void SetTags(const String_Map &tags);

    // inline functions
    inline void SetTagReplacer(Tag_Replacer *const replacer) 
    { p_replacer=replacer; }

    inline const Function_Set &Functions() const { return m_functions; }
    inline const Operator_Map &Operators() const { return m_operators; }

    inline Tag_Replacer *const TagReplacer() const { return p_replacer; }

  };// end of class Algebra_Interpreter

}// end of namespace ATOOLS

#define DEFINE_FUNCTION(NAME,TAG)					\
  class NAME: public Function {						\
  public:								\
    NAME();								\
    std::string Evaluate(const std::vector<std::string> &args) const;	\
  };									\
  NAME::NAME(): Function(TAG) {}					\
  std::string NAME::Evaluate(const std::vector<std::string> &args) const

#define DEFINE_UNARY_OPERATOR(NAME,TAG,PRIORITY)			\
  class NAME: public Operator {						\
  public:								\
    NAME();								\
    std::string Evaluate(const std::vector<std::string> &args) const;	\
  };									\
  NAME::NAME(): Operator(TAG,PRIORITY,false) {}				\
  std::string NAME::Evaluate(const std::vector<std::string> &args) const

#define DEFINE_BINARY_OPERATOR(NAME,TAG,PRIORITY)			\
  class NAME: public Operator {						\
  public:								\
    NAME();								\
    std::string Evaluate(const std::vector<std::string> &args) const;	\
  };									\
  NAME::NAME(): Operator(TAG,PRIORITY,true) {}				\
  std::string NAME::Evaluate(const std::vector<std::string> &args) const

#define DEFINE_INTERPRETER_FUNCTION(NAME)				\
  class NAME: public Interpreter_Function {				\
  public:								\
    inline NAME(Algebra_Interpreter *interpreter):			\
      Interpreter_Function(interpreter) {}				\
    std::string Interprete(const std::string &expr);			\
  };									\
  std::string NAME::Interprete(const std::string &expr)

#endif
