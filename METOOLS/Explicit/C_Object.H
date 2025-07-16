#ifndef METOOLS_Explicit_C_Object_H
#define METOOLS_Explicit_C_Object_H

#include "ATOOLS/Math/MyComplex.H"
#include <iostream>
#include <vector>

namespace METOOLS {

  /*!
    \class CObject
    \brief The base class of all currents.

    The functionality of this class is defined in its derived classes.
  */
  class CObject {
  protected:

    /*!
      \var int m_c[2]
      \brief The color and anti-colour indices of the current. The two elements are for color 
      and anti-color respectively. They take the values of +/-[1,2,3] for the SM.
      \var int m_h
      \brief The helicity state identifier of the current. This is later mapped to an entry 
      in a vector containing all possible helicity states, ie. confgurations of helicities 
      of each particle.
      \var int m_s
      \brief Whether the current is part of an NLO subtraction or not.
    */
    int m_c[2], m_h, m_s;

  public:

    // constructor
    inline CObject() {}

    // destructor
    virtual ~CObject();

    // member functions
    virtual CObject* Copy() const = 0;

    virtual void Delete() = 0;

    virtual void Add(const CObject *c) = 0;
    virtual void Divide(const double &d) = 0;
    virtual void Multiply(const Complex &c) = 0;
    virtual void Invert() = 0;

    virtual bool IsZero() const = 0;

    // inline functions
    /*!
      \fn inline int &operator()(const int i)
      \brief Returns the color or anticolor, depending on the index (0 and 1 respectively)
     */
    inline int &operator()(const int i) { return m_c[i]; }
    /*!
      \fn inline int operator()(const int i) const
      \brief Returns the color or anticolor, depending on the index (0 and 1 respectively)
     */
    inline int operator()(const int i) const { return m_c[i]; }

    inline void SetH(const int &h) { m_h=h; }
    inline void SetS(const int &s) { m_s=s; }

    inline const int &H() const { return m_h; }
    inline const int &S() const { return m_s; }

    /*!
      \fn inline bool operator==(const CObject &o) const
      \brief Compares the color and m_s of the two currents, not the helicity.
     */
    inline bool operator==(const CObject &o) const
    { return m_c[0]==o.m_c[0] && m_c[1]==o.m_c[1] && m_s==o.m_s; }

    /*!
      \fn template <typename CType> inline CType *Get()
      \brief Returns a pointer of CType to this object.
     */
    template <typename CType> inline CType *Get()
    { return static_cast<CType*>((void*)this); }

  };// end of class CObject

  struct CObject_Vector: public std::vector<CObject*> {

    CObject_Vector(const size_t &n=0): std::vector<CObject*>(n) {}

    template <typename CType> inline std::vector<CType*> *Get()
    { return static_cast<std::vector<CType*>*>((void*)this); }
    template <typename CType> inline const std::vector<CType*> *Get() const
    { return static_cast<const std::vector<CType*>*>((void*)this); }

  };// end of struct CObject_Vector

  struct CObject_Matrix: public std::vector<CObject_Vector> {

    template <typename CType> inline std::vector<std::vector<CType*> > *Get()
    { return static_cast<std::vector<std::vector<CType*> >*>((void*)this); }
    template <typename CType> inline const std::vector<std::vector<CType*> > *Get() const
    { return static_cast<const std::vector<std::vector<CType*> >*>((void*)this); }

  };// end of struct CObject_Matrix

  std::ostream &operator<<(std::ostream &str,const CObject &s);

}// end of namespace METOOLS

#endif
