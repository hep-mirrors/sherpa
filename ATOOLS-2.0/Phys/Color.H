#ifndef Color_H
#define Color_H

#include <iostream>
#include "Node.H"
#include "MyComplex.H"

namespace ATOOLS {

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

    Complex m_n;

  public:
    
    // constructor
    Number(Expression *const owner,const Complex &n);

    // member functions
    bool Evaluate(Expression *const expression);
    void Print() const;

    Term *GetCopy(Expression *const expression) const;

    // inline functions
    inline Complex operator()() { return m_n; }

  };// end of class Number
  /*!
    \class Number
    \brief Represents a complex number.

    This class represents a complex number.
  */

  class Delta: public Term {
  private:

    size_t m_i, m_j;

  public:
    
    // constructor
    Delta(Expression *const owner,
	  const size_t &i,const size_t &j);

    // member functions
    bool Evaluate(Expression *const expression);
    /*!
      \fn bool Evaluate(Expression *const expression)
      \brief Contracts delta functions.

      This method contracts delta functions according to
      \f[
      \delta_{ij}\delta_{jk}\,=\;\delta_{ik}
      \f]
      and
      \f[
      \delta_{ii}\,=\;N_C\;.
      \f]
    */
    void Print() const;

    Term *GetCopy(Expression *const expression) const;

  };// end of class Delta
  /*!
    \class Delta
    \brief Represents a delta function.

    This class represents a delta function.
  */

  class Fundamental: public Term {
  private:

    size_t m_a, m_i, m_j;
    bool   m_fromf;
    /*!
      \var bool m_fromf
      \brief Tags that this objects stems from the expansion of 
      a matrix in the adjoint representation.
    */

  public:
    
    // constructor
    Fundamental(Expression *const owner,
		const size_t &a,const size_t &i,const size_t &j);

    // member functions
    bool Evaluate(Expression *const expression);
    /*!
      \fn bool Evaluate(Expression *const expression)
      \brief Expands the matrix.

      This method expands a contraction of two matrices 
      in the fundamental representation, \f$T^{a}_{ij}\f$ and 
      \f$T^{a}_{kl}\f$ in terms of delta functions, according to
      \f[
      T^{a}_{ij}T^{a}_{kl}\,=\;\frac{1}{2}\,\left(\,
      \delta_{il}\delta_{jk}-\frac{1}{N_C}\delta_{ij}\delta_{kl}
      \,\right)\;,
      \f]
      If the m_fromf flag is set, the \f$1/{N_C}\f$-term is dropped,
      since it does not contribute, see Adjoint.
    */
    void Print() const;

    Term *GetCopy(Expression *const expression) const;

    inline void SetFromF(const bool &fromf) { m_fromf=fromf; }

    inline const bool &FromF() const { return m_fromf; }

  };// end of class Fundamental
  /*!
    \class Fundamental
    \brief Represents a color matrix in the fundamental representation.

    This class represents a color matrix in the fundamental representation.
  */

  class Adjoint: public Term {
  private:

    size_t m_a, m_b, m_c;

    double SwapIndices(Adjoint *const adjoint);
    void   ReplacePair(Adjoint *const adjoint,
		       const size_t &i,const size_t &j,
		       const double &factor,
		       Expression *const expression);

  public:
    
    // constructor
    Adjoint(Expression *const owner,
	    const size_t &a,const size_t &b,const size_t &c);

    // member functions
    bool Evaluate(Expression *const expression);
    /*!
      \fn bool Evaluate(Expression *const expression)
      \brief Expands the matrix.

      This method expands a matrix of the adjoint representation, 
      \f$f^{abc}\f$ in terms of matrices in the fundamental representation.
      Since
      \f[
      f^{abc}\,=\;-2i\,{\rm tr}\left(\,T^a\left[T^b,T^c\right]\,\right)
      \,=\;-2i\,{\rm tr}\left(\,\left[T^a,T^b\right]T^c\,\right)
      \,=\;-2i\,{\rm tr}\left(\,\left[T^c,T^a\right]T^b\,\right)\;,
      \f]
      the generated Fundamental objects obtain a flag, which tags 
      that they stem from the expansion of an \f$f\f$ matrix 
      and therefore, in further evaluations terms of the form 
      \f$T^a_{ii}\f$ may be dropped.
    */
    void Print() const;

    Term *GetCopy(Expression *const expression) const;

  };// end of class Adjoint
  /*!
    \class Adjoint
    \brief Represents a color matrix in the adjoint representation.

    This class represents a color matrix in the adjoint representation.
  */

  class Expression: public ATOOLS::Node<Term*> {
  public:

    typedef std::vector<Term*> Term_Vector;

    typedef std::vector<ATOOLS::Node<Term*>*> Expression_Vector;

  private:

    Complex m_result;

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

    void PrintStatus(const bool endline=true,
		     const bool print=true);

    // inline functions
    inline Complex Result() const { return m_result; }

    inline size_t FIndex() { return ++m_findex; }
    inline size_t AIndex() { return ++m_aindex; }

  };// end of class Expression

}// end of namespace ATOOLS

#endif
