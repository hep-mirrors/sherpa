#ifndef ATOOLS_Math_Algebra_Interpreter_H
#define ATOOLS_Math_Algebra_Interpreter_H

#ifdef USING__Calc_only
#include "Tools.H"
#include "Node.H"
#include "Term.H"
#else
#include "ATOOLS/Org/STL_Tools.H"
#include "ATOOLS/Org/Node.H"
#include "ATOOLS/Math/Term.H"
#endif

#include <string>
#include <set>
#include <map>

namespace ATOOLS {

  class Algebra_Interpreter;

  class Function {
  protected:

    std::string m_tag;

    Algebra_Interpreter *p_interpreter;

  public:

    // constructor
    Function(const std::string &tag);

    // destructor
    virtual ~Function();

    // member functions
    virtual Term *Evaluate(const std::vector<Term*> &args) const;

    // inline functions
    inline const std::string &Tag() const { return m_tag; }

    inline void SetInterpreter(Algebra_Interpreter *const interpreter)
    { p_interpreter=interpreter; }

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

  class Interpreter_Function {
  protected:
    
    Algebra_Interpreter *p_interpreter;

  public:
    
    // constructor
    inline Interpreter_Function(Algebra_Interpreter *interpreter):
      p_interpreter(interpreter) {}

    // destructor
    virtual ~Interpreter_Function();
    
    // member functions
    virtual std::string Interprete(const std::string &expr) = 0;

  };// end of class Interpreter_Function

  class Tag_Replacer {
  public:

    // destructor
    virtual ~Tag_Replacer();

    // member functions
    std::string &KillBlanks(std::string& expr) const;

    virtual std::string  ReplaceTags(std::string &expr) const;    
    virtual Term        *ReplaceTags(Term *term) const;    

  };// end of class Tag_Replacer

  class Algebra_Interpreter: public Tag_Replacer {
  public:

    typedef std::map<std::string,std::string,String_Sort> String_Map;
    
    typedef std::pair<std::string,Function*> Function_Pair;
    typedef std::map<std::string,Function*>  Function_Map;
    typedef std::pair<size_t,Operator*>      Operator_Pair;
    typedef std::multimap<size_t,Operator*>  Operator_Map;

    typedef std::map<size_t,Interpreter_Function*> Interpreter_Map;

    typedef std::set<Term*> Term_Set;

  private:

    Function_Map m_functions, m_leafs;
    Operator_Map m_operators;

    Interpreter_Map       m_interpreters;
    Interpreter_Function *p_extractor;

    Tag_Replacer *p_replacer;

    String_Map m_tags;
    Term_Set   m_terms;

    Node<Function*> *p_root, *p_leaf;

    void PrintNode(Node<Function*> *const node) const;

  public:

    // constructors
    Algebra_Interpreter(const bool standard=true);

    // destructor
    virtual ~Algebra_Interpreter();

    // member functions
    std::string Interprete(const std::string &expr);
    std::string Iterate(const std::string &expr);

    Term *Calculate();
    Term *Iterate(Node<Function*> *const node);

    void AddFunction(Function *const f);
    void AddOperator(Operator *const b);

    void AddLeaf(Function *const f);
    void AddTerm(Term *const t);

    std::string  ReplaceTags(std::string &expr) const;    
    Term        *ReplaceTags(Term *expr) const;

    void AddTag(const std::string &tag,const std::string &value);
    void SetTags(const String_Map &tags);

    // inline functions
    inline void SetTagReplacer(Tag_Replacer *const replacer) 
    { p_replacer=replacer; }
    inline void SetLeaf(Node<Function*> *const leaf) 
    { p_leaf=leaf; }

    inline Interpreter_Function *const Extractor() const 
    { return p_extractor; }

    inline const Function_Map &Functions() const { return m_functions; }
    inline const Operator_Map &Operators() const { return m_operators; }

    inline Node<Function*> *const Leaf() const { return p_leaf; }

    inline Tag_Replacer *const TagReplacer() const { return p_replacer; }

    inline void PrintEquation() const { PrintNode(p_root); }

  };// end of class Algebra_Interpreter

}// end of namespace ATOOLS

#define DEFINE_FUNCTION(NAME,TAG)					\
  namespace ATOOLS { class NAME: public Function {                      \
  public:								\
    NAME();								\
    Term *Evaluate(const std::vector<Term*> &args) const;		\
  }; }									\
  NAME::NAME(): Function(TAG) {}

#define DEFINE_UNARY_OPERATOR(NAME,TAG,PRIORITY)			\
  namespace ATOOLS { class NAME: public Operator {			\
  public:								\
    NAME();								\
    Term *Evaluate(const std::vector<Term*> &args) const;		\
  }; }									\
  NAME::NAME(): Operator(TAG,PRIORITY,false) {}

#define DEFINE_BINARY_OPERATOR(NAME,TAG,PRIORITY)			\
  namespace ATOOLS { class NAME: public Operator {			\
  public:								\
    NAME();								\
    Term *Evaluate(const std::vector<Term*> &args) const;		\
  }; }									\
  NAME::NAME(): Operator(TAG,PRIORITY,true) {}

#define DEFINE_INTERPRETER_FUNCTION(NAME)				\
  namespace ATOOLS { class NAME: public Interpreter_Function {		\
  public:								\
    inline NAME(Algebra_Interpreter *interpreter):			\
      Interpreter_Function(interpreter) {}				\
    std::string Interprete(const std::string &expr);			\
  }; }									\
  std::string NAME::Interprete(const std::string &expr)

#endif
