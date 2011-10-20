#ifndef ATOOLS_Math_Term_H
#define ATOOLS_Math_Term_H

#include <string>

namespace ATOOLS {

  class Term {
  private:

    char m_type;

    std::string m_tag;
    size_t      m_id;

  public:

    inline Term(const char type): m_type(type) {}

    virtual ~Term();

    template <class TType> static Term *New(const TType &val);

    static Term *NewTerm(const std::string &tag);

    virtual void Delete() = 0;

    inline char Type() const { return m_type; }

    template <class TType> const TType &Get() const;
    template <class TType> void Set(const TType &val);

    void SetTerm(const std::string &val);

    void Print(std::ostream &ostr) const;

    Term *operator-() const;
    Term *operator!() const;

    Term *operator+(const Term &ref) const;
    Term *operator-(const Term &ref) const;
    Term *operator*(const Term &ref) const;
    Term *operator/(const Term &ref) const;

    Term *operator==(const Term &ref) const;
    Term *operator!=(const Term &ref) const;

    Term *operator<(const Term &ref) const;
    Term *operator>(const Term &ref) const;
    Term *operator<=(const Term &ref) const;
    Term *operator>=(const Term &ref) const;

    Term *operator%(const Term &ref) const;
    Term *operator<<(const Term &ref) const;
    Term *operator>>(const Term &ref) const;
    Term *operator&&(const Term &ref) const;
    Term *operator||(const Term &ref) const;
    Term *operator&(const Term &ref) const;
    Term *operator^(const Term &ref) const;
    Term *operator|(const Term &ref) const;

    Term *Real() const;
    Term *Imag() const;
    Term *Conj() const;
    Term *Comp(const Term &i) const;

    Term *Perp() const;
    Term *Plus() const;
    Term *Minus() const;

    Term *PPlus() const;
    Term *PMinus() const;
    Term *Abs2() const;
    Term *Mass() const;
    Term *PSpat() const;
    Term *PPerp() const;
    Term *PPerp2() const;
    Term *MPerp() const;
    Term *MPerp2() const;
    Term *Theta() const;
    Term *Eta() const;
    Term *Y() const;
    Term *Phi() const;

    Term *PPerp(const Term &ref) const;
    Term *Theta(const Term &ref) const;
    Term *DEta(const Term &ref) const;
    Term *DPhi(const Term &ref) const;

    inline void SetTag(const std::string &tag) { m_tag=tag; }
    inline void SetId(const size_t &id)        { m_id=id;   }

    inline const std::string &Tag() const { return m_tag; }
    inline const size_t      &Id() const  { return m_id;  }

  };// end of class Term

  std::ostream &operator<<(std::ostream &ostr,const Term &t);

  Term *TVec4D(const Term &t0,const Term &t1,
	       const Term &t2,const Term &t3);

  Term *TExp(const Term &t);
  Term *TLog(const Term &t);
  Term *TLog10(const Term &t);
  Term *TAbs(const Term &t);
  Term *TSgn(const Term &t);
  Term *TTheta(const Term &t1);
  Term *TSqr(const Term &t);
  Term *TSqrt(const Term &t);
  Term *TSin(const Term &t);
  Term *TCos(const Term &t);
  Term *TTan(const Term &t);
  Term *TSinh(const Term &t);
  Term *TCosh(const Term &t);
  Term *TTanh(const Term &t);
  Term *TASin(const Term &t);
  Term *TACos(const Term &t);
  Term *TATan(const Term &t);

  Term *TPow(const Term &t1,const Term &t2);
  Term *TMin(const Term &t1,const Term &t2);
  Term *TMax(const Term &t1,const Term &t2);

}// end of namespace ATOOLS

#endif
