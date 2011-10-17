#ifdef USING__Calc_only
#include "Term.H"
#include "Tools.H"
#include "Vector.H"
#else
#include "ATOOLS/Org/STL_Tools.H"
#include "ATOOLS/Math/Term.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#endif

namespace ATOOLS {

  inline bool IsNum(const char c) { return c>=48 && c<=57; }

  bool IsAlpha(const std::string &expr) 
  {
    bool ad(true), ae(true);
    for (size_t i=0;i<expr.length();++i) 
      if (!IsNum(expr[i])) {
	if (expr[i]=='.' && ad) {
	  ad=false;
	  continue;
	}
	if ((expr[i]=='e' || expr[i]=='E') && 
	    ae && i<expr.length()-1) {
	  if (expr[i+1]=='+' || expr[i+1]=='-') ++i;
	  ae=ad=false;
	  continue;
	}
	return true;
      }
    return false;
  }

  template <class _Type>
  class TermDelete_Vector:
    public std::vector<_Type*> {
  private:

#ifdef USING__Threading
    pthread_mutex_t m_mtx;
#endif

  public:

    TermDelete_Vector()
    {
#ifdef USING__Threading
      pthread_mutex_init(&m_mtx,NULL);
#endif
    }

    ~TermDelete_Vector()
    {
#ifdef USING__Threading
      pthread_mutex_destroy(&m_mtx);
#endif
      while (!this->empty()) {
	delete this->back();
	this->pop_back();
      }
    }

    inline void MtxLock()
    {
#ifdef USING__Threading
      pthread_mutex_lock(&m_mtx);
#endif
    }

    inline void MtxUnLock()
    {
#ifdef USING__Threading
      pthread_mutex_unlock(&m_mtx);
#endif
    }

  };// end of class TermDelete_Vector

  class DTerm: public Term {
  private:

    double m_this;
    friend class Term;

    static AutoDelete_Vector<DTerm> s_terms;

    inline DTerm(const double &val=0.0): 
      Term('D'), m_this(val) {}

  public:

    static DTerm *New(const double &val=0.0)
    {
      s_terms.MtxLock();
      if (s_terms.empty()) {
	s_terms.MtxUnLock();
	return new DTerm(val);
      }
      DTerm *term(s_terms.back());
      s_terms.pop_back();
      s_terms.MtxUnLock();
      term->m_this=val;
      return term;
    }

    void Delete()
    {
      s_terms.MtxLock();
      s_terms.push_back(this);
      s_terms.MtxUnLock();
    }

  };// end of class DTerm

  struct CTerm: public Term {
  private:

    Complex m_this;
    friend class Term;

    static AutoDelete_Vector<CTerm> s_terms;

    inline CTerm(const Complex &val=Complex(0.0,0.0)): 
      Term('C'), m_this(val) {}

  public:

    static CTerm *New(const Complex &val=Complex(0.0,0.0))
    {
      s_terms.MtxLock();
      if (s_terms.empty()) {
	s_terms.MtxUnLock();
	return new CTerm(val);
      }
      CTerm *term(s_terms.back());
      s_terms.pop_back();
      s_terms.MtxUnLock();
      term->m_this=val;
      return term;
    }

    void Delete()
    {
      s_terms.MtxLock();
      s_terms.push_back(this);
      s_terms.MtxUnLock();
    }

  };// end of class CTerm

  struct DV4Term: public Term {
  private:

    Vec4D m_this;
    friend class Term;

    static AutoDelete_Vector<DV4Term> s_terms;

    inline DV4Term(const Vec4D &val): 
      Term('V'), m_this(val) {}

  public:

    static DV4Term *New(const Vec4D &val)
    {
      s_terms.MtxLock();
      if (s_terms.empty()) {
	s_terms.MtxUnLock();
	return new DV4Term(val);
      }
      DV4Term *term(s_terms.back());
      s_terms.pop_back();
      s_terms.MtxUnLock();
      term->m_this=val;
      return term;
    }

    void Delete()
    {
      s_terms.MtxLock();
      s_terms.push_back(this);
      s_terms.MtxUnLock();
    }

  };// end of class DV4Term

  struct STerm: public Term {
  private:

    std::string m_this;
    friend class Term;

    static AutoDelete_Vector<STerm> s_terms;

    inline STerm(const std::string &val): 
      Term('S'), m_this(val) {}

  public:

    static STerm *New(const std::string &val)
    {
      s_terms.MtxLock();
      if (s_terms.empty()) {
	s_terms.MtxUnLock();
	return new STerm(val);
      }
      STerm *term(s_terms.back());
      s_terms.pop_back();
      s_terms.MtxUnLock();
      term->m_this=val;
      return term;
    }

    void Delete()
    {
      s_terms.MtxLock();
      s_terms.push_back(this);
      s_terms.MtxUnLock();
    }

  };// end of class STerm

  AutoDelete_Vector<DTerm> DTerm::s_terms;
  AutoDelete_Vector<CTerm> CTerm::s_terms;
  AutoDelete_Vector<DV4Term> DV4Term::s_terms;
  AutoDelete_Vector<STerm> STerm::s_terms;

  template <> const double  &Term::Get() const
  { return static_cast<const DTerm *>(this)->m_this; }
  template <> const Complex &Term::Get() const
  { return static_cast<const CTerm *>(this)->m_this; }
  template <> const Vec4D &Term::Get() const
  { return static_cast<const DV4Term *>(this)->m_this; }
  template <> const std::string &Term::Get() const
  { return static_cast<const STerm *>(this)->m_this; }

  Term::~Term() {}

  void Term::Print(std::ostream &ostr) const
  {
    if (m_type=='S') ostr<<Get<std::string>();
    else if (m_type=='V') ostr<<Get<Vec4D>();
    else if (m_type=='C') ostr<<Get<Complex>();
    else ostr<<Get<double>();
  }

  std::ostream &operator<<(std::ostream &ostr,const Term &t)
  {
    t.Print(ostr);
    return ostr;
  }

  template <> void Term::Set(const double &val)
  { static_cast<DTerm *>(this)->m_this=val; }
  template <> void Term::Set(const Complex &val)
  { static_cast<CTerm *>(this)->m_this=val; }
  template <> void Term::Set(const Vec4D &val)
  { static_cast<DV4Term *>(this)->m_this=val; }
  template <> void Term::Set(const std::string &val)
  { static_cast<STerm *>(this)->m_this=val; }

  void Term::SetTerm(const std::string &tag) 
  { 
    if (tag[0]!='(') {
      if (tag[0]=='"' && tag[tag.length()-1]=='"')
	static_cast<STerm*>(this)->m_this=tag.substr(1,tag.length()-2);
      else if (IsAlpha(tag)) static_cast<STerm*>(this)->m_this=(tag);
      else static_cast<DTerm*>(this)->m_this=ToType<double>(tag); 
    }
    else {
      size_t pos(tag.find(','));
      if (pos==std::string::npos) THROW(fatal_error,"Invalid syntax");
      if ((pos=tag.find(',',pos+1))!=std::string::npos)
	static_cast<DV4Term*>(this)->m_this=ToType<Vec4D>(tag);
      else static_cast<CTerm*>(this)->m_this=ToType<Complex>(tag);
    }
  }

  template <> Term *Term::New(const double &val)  { return DTerm::New(val); }
  template <> Term *Term::New(const Complex &val) { return CTerm::New(val); }
  template <> Term *Term::New(const Vec4D &val) { return DV4Term::New(val); }
  template <> Term *Term::New(const std::string &val) { return STerm::New(val); }

  Term *Term::NewTerm(const std::string &tag) 
  { 
    if (tag[0]!='(') {
      if (tag[0]=='"' && tag[tag.length()-1]=='"')
	return STerm::New(tag.substr(1,tag.length()-2)); 
      if (IsAlpha(tag)) return STerm::New(tag); 
      return DTerm::New(ToType<double>(tag)); 
    }
    else {
      size_t pos(tag.find(','));
      if (pos==std::string::npos) THROW(fatal_error,"Invalid syntax");
      if ((pos=tag.find(',',pos+1))!=std::string::npos)
	return DV4Term::New(ToType<Vec4D>(tag));
      return CTerm::New(ToType<Complex>(tag));
    }
  }

  Term *Term::operator-() const
  {
    if (m_type=='S') THROW(fatal_error,"Invalid syntax");
    if (m_type=='V') return DV4Term::New(-Get<Vec4D>());
    if (m_type=='C') return CTerm::New(-Get<Complex>());
    return DTerm::New(-Get<double>());
  }

  Term *Term::operator!() const
  {
    if (m_type=='C') return CTerm::New(!(int)(Get<Complex>().real()));
    if (m_type=='D') return DTerm::New(!(int)(Get<double>()));
    THROW(fatal_error,"Invalid syntax");
    return NULL;
  }

  Term *TVec4D(const Term &t0,const Term &t1,
	       const Term &t2,const Term &t3)
  {
    if (t0.Type()=='V' || t0.Type()=='C' || t0.Type()=='S' ||
	t1.Type()=='V' || t1.Type()=='C' || t1.Type()=='S' ||
	t2.Type()=='V' || t2.Type()=='C' || t2.Type()=='S' ||
	t3.Type()=='V' || t3.Type()=='C' || t3.Type()=='S')
      THROW(fatal_error,"Invalid syntax");
    return DV4Term::New(Vec4D(t0.Get<double>(),t1.Get<double>(),
			      t2.Get<double>(),t3.Get<double>()));
  }

#define DEFINE_BINARY_STERM_OPERATOR(OP)\
  Term *Term::operator OP(const Term &ref) const\
  {\
    if (m_type=='S' || ref.m_type=='S')\
      THROW(fatal_error,"Invalid syntax");\
    if (m_type=='V') {\
      if (ref.m_type=='V')\
	return DV4Term::New(Get<Vec4D>() OP ref.Get<Vec4D>());\
      THROW(fatal_error,"Invalid syntax");\
      return NULL;\
    }\
    if (m_type=='C') {\
      if (ref.m_type=='C')\
	return CTerm::New(Get<Complex>() OP ref.Get<Complex>());\
      if (ref.m_type=='D')\
        return CTerm::New(Get<Complex>() OP ref.Get<double>());\
      THROW(fatal_error,"Invalid syntax");\
      return NULL;\
    }\
    if (ref.m_type=='C')\
      return CTerm::New(Get<double>() OP ref.Get<Complex>());\
    return DTerm::New(Get<double>() OP ref.Get<double>());\
  }
 
  DEFINE_BINARY_STERM_OPERATOR(+)
  DEFINE_BINARY_STERM_OPERATOR(-)

  Term *Term::operator*(const Term &ref) const
  {
    if (m_type=='S' || ref.m_type=='S')
      THROW(fatal_error,"Invalid syntax");
    if (m_type=='V') {
      if (ref.m_type=='V')
	return DTerm::New(Get<Vec4D>()*ref.Get<Vec4D>());
      if (ref.m_type=='D')
	return DV4Term::New(ref.Get<double>()*Get<Vec4D>());
      THROW(fatal_error,"Invalid syntax");
      return NULL;
    }
    if (m_type=='C') {
      if (ref.m_type=='C')
	return CTerm::New(Get<Complex>()*ref.Get<Complex>());
      if (ref.m_type=='D')
        return CTerm::New(Get<Complex>()*ref.Get<double>());
      THROW(fatal_error,"Invalid syntax");
      return NULL;
    }
    if (ref.m_type=='V')
      return DV4Term::New(Get<double>()*ref.Get<Vec4D>());
    if (ref.m_type=='C')
      return CTerm::New(Get<double>()*ref.Get<Complex>());
    return DTerm::New(Get<double>()*ref.Get<double>());
  }
 
  Term *Term::operator/(const Term &ref) const
  {
    if (m_type=='S' || ref.m_type=='S')
      THROW(fatal_error,"Invalid syntax");
    if (m_type=='V') {
      if (ref.m_type=='D')
	return DV4Term::New(1.0/ref.Get<double>()*Get<Vec4D>());
      THROW(fatal_error,"Invalid syntax");
      return NULL;
    }
    if (m_type=='C') {
      if (ref.m_type=='C')
	return CTerm::New(Get<Complex>()/ref.Get<Complex>());
      if (ref.m_type=='D')
        return CTerm::New(Get<Complex>()/ref.Get<double>());
      THROW(fatal_error,"Invalid syntax");
      return NULL;
    }
    if (ref.m_type=='C')
      return CTerm::New(Get<double>()/ref.Get<Complex>());
    return DTerm::New(Get<double>()/ref.Get<double>());
  }
 
#define DEFINE_BINARY_BTERM_OPERATOR(OP)\
  Term *Term::operator OP(const Term &ref) const\
  {\
    if (m_type=='S' || ref.m_type=='S') {\
      if (m_type!='S' || ref.m_type!='S')\
	THROW(fatal_error,"Invalid syntax");\
      return DTerm::New(Get<std::string>() OP ref.Get<std::string>());\
    }\
    if (m_type=='V' || ref.m_type=='V')\
      THROW(fatal_error,"Invalid syntax");\
    if (m_type=='C') {\
      if (ref.m_type=='C')\
	return DTerm::New(Get<Complex>() OP ref.Get<Complex>());\
      return DTerm::New(Get<Complex>() OP ref.Get<double>());\
    }\
    if (ref.m_type=='C')\
      return DTerm::New(Get<double>() OP ref.Get<Complex>());\
    return DTerm::New(Get<double>() OP ref.Get<double>());\
  }

  bool operator<(const Complex &a,const Complex &b)
  { return a.real()<b.real() && a.imag()<b.imag(); }
  bool operator>(const Complex &a,const Complex &b)
  { return a.real()>b.real() && a.imag()>b.imag(); }

  bool operator<=(const Complex &a,const Complex &b) { return !(a>b); }
  bool operator>=(const Complex &a,const Complex &b) { return !(a<b); }

  DEFINE_BINARY_BTERM_OPERATOR(==)
  DEFINE_BINARY_BTERM_OPERATOR(!=)
  DEFINE_BINARY_BTERM_OPERATOR(<=)
  DEFINE_BINARY_BTERM_OPERATOR(>=)
  DEFINE_BINARY_BTERM_OPERATOR(<)
  DEFINE_BINARY_BTERM_OPERATOR(>)

#define DEFINE_BINARY_ITERM_OPERATOR(OP)\
  Term *Term::operator OP(const Term &ref) const\
  {\
    if (m_type=='S' || ref.m_type=='S' ||\
	m_type=='V' || ref.m_type=='V')\
      THROW(fatal_error,"Invalid syntax");\
    if (m_type=='C') {\
      if (ref.m_type=='C') \
	return DTerm::New((long int)(Get<Complex>().real()) OP\
			  (long int)(ref.Get<Complex>().real()));\
      return DTerm::New((long int)(Get<Complex>().real()) OP\
			(long int)(ref.Get<double>()));\
    }\
    if (ref.m_type=='C')\
      return DTerm::New((long int)(Get<double>()) OP\
			(long int)(ref.Get<Complex>().real()));\
    return DTerm::New((long int)(Get<double>()) OP\
		      (long int)(ref.Get<double>()));\
  }

  DEFINE_BINARY_ITERM_OPERATOR(%)
  DEFINE_BINARY_ITERM_OPERATOR(<<)
  DEFINE_BINARY_ITERM_OPERATOR(>>)
  DEFINE_BINARY_ITERM_OPERATOR(&&)
  DEFINE_BINARY_ITERM_OPERATOR(||)
  DEFINE_BINARY_ITERM_OPERATOR(&)
  DEFINE_BINARY_ITERM_OPERATOR(^)
  DEFINE_BINARY_ITERM_OPERATOR(|)

#define DEFINE_UNARY_DTERM_FUNCTION(NAME,FUNC)\
  Term *NAME(const Term &t)\
  {\
    if (t.Type()=='S' || t.Type()=='V')\
      THROW(fatal_error,"Invalid syntax");\
    if (t.Type()=='C') return NULL;\
    return DTerm::New(FUNC(t.Get<double>()));\
  }

  DEFINE_UNARY_DTERM_FUNCTION(TSgn,Sign)
  DEFINE_UNARY_DTERM_FUNCTION(TTheta,Theta)
  DEFINE_UNARY_DTERM_FUNCTION(TASin,asin)
  DEFINE_UNARY_DTERM_FUNCTION(TACos,acos)
  DEFINE_UNARY_DTERM_FUNCTION(TATan,atan)

#define DEFINE_UNARY_TERM_FUNCTION(NAME,FUNC)\
  Term *NAME(const Term &t)\
  {\
    if (t.Type()=='S' || t.Type()=='V')\
      THROW(fatal_error,"Invalid syntax");\
    if (t.Type()=='C')\
      return CTerm::New(FUNC(t.Get<Complex>()));\
    return DTerm::New(FUNC(t.Get<double>()));\
  }

  DEFINE_UNARY_TERM_FUNCTION(TExp,exp)
  DEFINE_UNARY_TERM_FUNCTION(TLog,log)
  DEFINE_UNARY_TERM_FUNCTION(TLog10,log10)
  DEFINE_UNARY_TERM_FUNCTION(TAbs,std::abs)
  DEFINE_UNARY_TERM_FUNCTION(TSqr,sqr)
  DEFINE_UNARY_TERM_FUNCTION(TSqrt,sqrt)
  DEFINE_UNARY_TERM_FUNCTION(TSin,sin)
  DEFINE_UNARY_TERM_FUNCTION(TCos,cos)
  DEFINE_UNARY_TERM_FUNCTION(TTan,tan)
  DEFINE_UNARY_TERM_FUNCTION(TSinh,sinh)
  DEFINE_UNARY_TERM_FUNCTION(TCosh,cosh)
  DEFINE_UNARY_TERM_FUNCTION(TTanh,tanh)

#define DEFINE_BINARY_DTERM_FUNCTION(NAME,FUNC)\
  Term *NAME(const Term &t1,const Term &t2)\
  {\
    if (t1.Type()!='D' || t2.Type()!='D')\
      THROW(fatal_error,"Invalid syntax");\
    return DTerm::New(FUNC(t1.Get<double>(),t2.Get<double>()));\
  }

  DEFINE_BINARY_DTERM_FUNCTION(TMin,Min)
  DEFINE_BINARY_DTERM_FUNCTION(TMax,Max)

#define DEFINE_BINARY_TERM_FUNCTION(NAME,FUNC)\
  Term *NAME(const Term &t1,const Term &t2)\
  {\
    if (t1.Type()=='S' || t2.Type()=='S' ||\
	t1.Type()=='V' || t2.Type()=='V')\
      THROW(fatal_error,"Invalid syntax");\
    if (t1.Type()=='C') {\
      if (t2.Type()=='C')\
        return CTerm::New(FUNC(t1.Get<Complex>(),t2.Get<Complex>()));\
      return CTerm::New(FUNC(t1.Get<Complex>(),t2.Get<double>()));\
    }\
    if (t2.Type()=='C')\
      return CTerm::New(FUNC(t1.Get<double>(),t2.Get<Complex>()));\
    return DTerm::New(FUNC(t1.Get<double>(),t2.Get<double>()));\
  }

  DEFINE_BINARY_TERM_FUNCTION(TPow,pow)

  Term *Term::Real() const
  {
    if (m_type=='S' || m_type=='V' || m_type=='D')
      THROW(fatal_error,"Invalid syntax");
    return DTerm::New(Get<Complex>().real());
  }

  Term *Term::Imag() const
  {
    if (m_type=='S' || m_type=='V' || m_type=='D')
      THROW(fatal_error,"Invalid syntax");
    return DTerm::New(Get<Complex>().imag());
  }

  Term *Term::Conj() const
  {
    if (m_type=='S' || m_type=='V' || m_type=='D')
      THROW(fatal_error,"Invalid syntax");
    return CTerm::New(std::conj(Get<Complex>()));
  }
  
  Term *Term::Comp(const Term &i) const
  {
    if (m_type=='V' && i.m_type=='D') return 
      DTerm::New(Get<Vec4D>()[(int)(i.Get<double>())]);
    THROW(fatal_error,"Invalid syntax");
    return NULL;
  }

#define DEFINE_UNARY_VVTERM_FUNCTION(FUNC)\
  Term *Term::FUNC() const\
  {\
    if (m_type=='V') return DV4Term::New(Get<Vec4D>().FUNC());\
    THROW(fatal_error,"Invalid syntax");\
    return NULL;\
  }

  DEFINE_UNARY_VVTERM_FUNCTION(Perp)
  DEFINE_UNARY_VVTERM_FUNCTION(Plus)
  DEFINE_UNARY_VVTERM_FUNCTION(Minus)

#define DEFINE_UNARY_VTERM_FUNCTION(FUNC)\
  Term *Term::FUNC() const\
  {\
    if (m_type=='V') return DTerm::New(Get<Vec4D>().FUNC());\
    THROW(fatal_error,"Invalid syntax");\
    return NULL;\
  }

  DEFINE_UNARY_VTERM_FUNCTION(PPlus)
  DEFINE_UNARY_VTERM_FUNCTION(PMinus)
  DEFINE_UNARY_VTERM_FUNCTION(Abs2)
  DEFINE_UNARY_VTERM_FUNCTION(Mass)
  DEFINE_UNARY_VTERM_FUNCTION(PPerp)
  DEFINE_UNARY_VTERM_FUNCTION(PPerp2)
  DEFINE_UNARY_VTERM_FUNCTION(MPerp)
  DEFINE_UNARY_VTERM_FUNCTION(MPerp2)
  DEFINE_UNARY_VTERM_FUNCTION(Theta)
  DEFINE_UNARY_VTERM_FUNCTION(Eta)
  DEFINE_UNARY_VTERM_FUNCTION(Y)
  DEFINE_UNARY_VTERM_FUNCTION(Phi)

#define DEFINE_BINARY_VTERM_FUNCTION(FUNC)\
  Term *Term::FUNC(const Term &ref) const\
  {\
    if (m_type=='V' && ref.m_type=='V')\
      return DTerm::New(Get<Vec4D>().FUNC(ref.Get<Vec4D>()));\
    THROW(fatal_error,"Invalid syntax");\
    return NULL;\
  }

  DEFINE_BINARY_VTERM_FUNCTION(PPerp)
  DEFINE_BINARY_VTERM_FUNCTION(Theta)
  DEFINE_BINARY_VTERM_FUNCTION(DEta)
  DEFINE_BINARY_VTERM_FUNCTION(DPhi)

}// end of namespace ATOOLS

