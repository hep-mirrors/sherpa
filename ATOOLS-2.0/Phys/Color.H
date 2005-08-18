#ifndef Color_H
#define Color_H
/*!
  \file Color.H
  \brief Declares the classes Color_Term,
  CNumber, Fundamental, Adjoint and Expression.
*/

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

  class Color_Term {
  private:

    ctt::type m_type;

  public:
    
    // constructor
    inline Color_Term(const ctt::type &type): 
      m_type(type) {}

    // destructor
    virtual ~Color_Term();

    // member functions
    virtual bool Evaluate(Expression *const expression);
    virtual void Print() const;

    virtual Color_Term *GetCopy() const;

    // inline functions
    inline ctt::type Type() const { return m_type; }

  };// end of class Color_Term
  /*!
    \class Color_Term
    \brief The base class of all objects in color space.

    This class is the base class of all classes
    representing objects in color space.
  */

  class CNumber: public Color_Term {
  private:

    Complex m_n;

  public:
    
    // constructor
    inline CNumber(const Complex &n):
      Color_Term(ctt::number), m_n(n) {}

    // member functions
    bool Evaluate(Expression *const expression);
    void Print() const;

    Color_Term *GetCopy() const;

    // inline functions
    inline Complex operator()() { return m_n; }

  };// end of class CNumber
  /*!
    \class CNumber
    \brief Represents a complex number.

    This class represents a complex number.
  */

  class Delta: public Color_Term {
  private:

    size_t m_i, m_j;

  public:
    
    // constructor
    inline Delta(const size_t &i,const size_t &j):
      Color_Term(ctt::delta), m_i(i), m_j(j) {}

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

    Color_Term *GetCopy() const;

  };// end of class Delta
  /*!
    \class Delta
    \brief Represents a delta function.

    This class represents a delta function.
  */

  class Fundamental: public Color_Term {
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
    inline Fundamental(const size_t &a,const size_t &i,const size_t &j,
		       const bool &fromf=false):
      Color_Term(ctt::fundamental), m_a(a), m_i(i), m_j(j), 
      m_fromf(fromf) {}

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

    Color_Term *GetCopy() const;

    inline const bool &FromF() const { return m_fromf; }

  };// end of class Fundamental
  /*!
    \class Fundamental
    \brief Represents a color matrix in the fundamental representation.

    This class represents a color matrix in the fundamental representation.
  */

  class Adjoint: public Color_Term {
  private:

    size_t m_a, m_b, m_c;

  public:
    
    // constructor
    inline Adjoint(const size_t &a,const size_t &b,const size_t &c):
      Color_Term(ctt::adjoint), m_a(a), m_b(b), m_c(c) {}

    // member functions
    bool Evaluate(Expression *const expression);
    /*!
      \fn bool Evaluate(Expression *const expression)
      \brief Expands the matrix.

      This method expands a matrix of the adjoint representation, 
      \f$f^{abc}\f$ in terms of matrices in the fundamental representation,
      according to
      \f[
      f^{abc}\,=\;-2i\,{\rm tr}\left(\,T^a\left[T^b,T^c\right]\,\right)
      \,=\;-2i\,{\rm tr}\left(\,T^c\left[T^a,T^b\right]\,\right)
      \,=\;-2i\,{\rm tr}\left(\,T^b\left[T^c,T^a\right]\,\right)\;.
      \f]
      The generated Fundamental objects, such as e.g. \f$T^a_{ij}\f$, 
      obtain a flag, which tags that they stem from the expansion 
      of an \f$f\f$ matrix and therefore, in the evaluation of 
      \f$T^a_{ij}T^a_{kl}\f$ the \f$1/N_C\f$-term may be dropped,
      see Fundamental.
    */
    void Print() const;

    Color_Term *GetCopy() const;

  };// end of class Adjoint
  /*!
    \class Adjoint
    \brief Represents a color matrix in the adjoint representation.

    This class represents a color matrix in the adjoint representation.
  */

  class Expression: public ATOOLS::Node<Color_Term*> {
  public:

    typedef std::vector<Color_Term*> Color_Term_Vector;

    typedef std::vector<ATOOLS::Node<Color_Term*>*> Expression_Vector;

  private:

    Complex m_result;

    size_t m_findex, m_aindex;
    size_t m_evaluated, m_cindex;

  public:
    
    // constructor
    inline Expression(): 
      Node<Color_Term*>(NULL,true), m_result(0.0,0.0),
      m_findex(0), m_aindex(0), m_evaluated(0), m_cindex(0) {}

    Expression(const std::string &expression);

    // destructor
    inline ~Expression()
    {
      for (Color_Term_Vector::iterator tit(begin());tit!=end();++tit) 
	delete *tit;
    }

    // member functions
    bool Evaluate();
    void Print();

    inline void Add(Expression *const expression)
    {
      (*this)().push_back(expression);
      (*expression)<<this;
    }

    size_t Size();
    size_t Evaluated();

    Expression *GetCopy() const;

    void PrintStatus(const bool endline=true,
		     const bool print=true);

    // inline functions
    inline Complex Result() const { return m_result; }

    inline size_t FIndex() { return ++m_findex; }
    inline size_t AIndex() { return ++m_aindex; }

    inline size_t CIndex() { return m_cindex; }

  };// end of class Expression
  /*!
    \class Expression
    \brief Represents a single color amplitude.
    
    This class represents a single color amplitude.
  */

}// end of namespace ATOOLS

#endif
