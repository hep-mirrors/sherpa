#ifndef Color_H
#define Color_H

#include <iostream>
#include "Node.H"

namespace COLOR {

  const size_t NC(3);

  struct ctt {
    
    enum type {
      number      = 0,
      delta       = 1,
      fundamental = 2,
      adjoint     = 4
    };

  };// end of struct ctt

  class Expression;

  class Term {
  private:

    ctt::type m_type;

    Expression *p_owner;

  public:
    
    // constructor
    Term(const ctt::type &type,Expression *const owner);

    // destructor
    virtual ~Term();

    // member functions
    virtual bool Evaluate(Expression *const expression);
    virtual void Print() const;

    virtual Term *GetCopy(Expression *const expression) const;

    // inline functions
    inline ctt::type Type() const { return m_type; }

  };// end of class Term

  class Number: public Term {
  private:

    double m_n;

  public:
    
    // constructor
    Number(Expression *const owner,const double n);

    // member functions
    bool Evaluate(Expression *const expression);
    void Print() const;

    Term *GetCopy(Expression *const expression) const;

    // inline functions
    inline double operator()() { return m_n; }

  };// end of class Number

  class Delta: public Term {
  private:

    size_t m_i, m_j;

  public:
    
    // constructor
    Delta(Expression *const owner,
	  const size_t i,const size_t j);

    // member functions
    bool Evaluate(Expression *const expression);
    void Print() const;

    Term *GetCopy(Expression *const expression) const;

  };// end of class Delta

  class Fundamental: public Term {
  private:

    size_t m_a, m_i, m_j;

  public:
    
    // constructor
    Fundamental(Expression *const owner,
		const size_t a,const size_t i,const size_t j);

    // member functions
    bool Evaluate(Expression *const expression);
    void Print() const;

    Term *GetCopy(Expression *const expression) const;

  };// end of class Fundamental

  class Adjoint: public Term {
  private:

    size_t m_a, m_b, m_c;

  public:
    
    // constructor
    Adjoint(Expression *const owner,
	    const size_t a,const size_t b,const size_t c);

    // member functions
    bool Evaluate(Expression *const expression);
    void Print() const;

    Term *GetCopy(Expression *const expression) const;

  };// end of class Adjoint

  class Expression: public ATOOLS::Node<Term*> {
  public:

    typedef std::vector<Term*> Term_Vector;

    typedef std::vector<ATOOLS::Node<Term*>*> Expression_Vector;

  private:

    double m_result;

    size_t m_findex, m_aindex;
    size_t m_evaluated;

  public:
    
    // constructor
    Expression();
    Expression(const std::string &expression);

    // destructor
    ~Expression();

    // member functions
    bool Evaluate();
    void Print();

    size_t Add(Expression *const expression);

    size_t Size();
    size_t Evaluated();

    Expression *GetCopy() const;

    // inline functions
    inline double Result() const { return m_result; }

    inline size_t FIndex() { return ++m_findex; }
    inline size_t AIndex() { return ++m_aindex; }

  };// end of class Expression

}// end of namespace COLOR

#endif
