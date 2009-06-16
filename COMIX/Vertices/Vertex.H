#ifndef COMIX_Vertices_Vertex_H
#define COMIX_Vertices_Vertex_H

#include "COMIX/Amplitude/Vertex_Base.H"
#include "HELICITIES/Main/C_Scalar.H"
#include "HELICITIES/Main/C_Spinor.H"
#include "HELICITIES/Main/C_Vector.H"
#include "HELICITIES/Main/C_Tensor.H"

namespace COMIX {

  template <typename SType>
  class SSS_Vertex: public Vertex_Base {
  public:

    typedef ATOOLS::CScalar<SType>   CScalarType;
    typedef std::vector<CScalarType> CScalarType_Vector;

  protected:
    
    ATOOLS::CScalar<SType> Lorentz
    (ATOOLS::CScalar<SType> a,ATOOLS::CScalar<SType> b);
    virtual void Evaluate(const ATOOLS::CScalar<SType> &a,
			  const ATOOLS::CScalar<SType> &b) = 0;
    
  public:
    
    SSS_Vertex(const Vertex_Key &key,
	       const size_t &oew=0,const size_t &oqcd=0);
    
    void Evaluate();

  };// end of class SSS_Vertex

  template <typename SType>
  class FFS_Vertex: public Vertex_Base {
  public:

    typedef ATOOLS::CScalar<SType>   CScalarType;
    typedef std::vector<CScalarType> CScalarType_Vector;

    typedef ATOOLS::CSpinor<SType>   CSpinorType;
    typedef std::vector<CSpinorType> CSpinorType_Vector;

  protected:
    
    int m_dir;

    ATOOLS::CScalar<SType> 
    Lorentz(ATOOLS::CSpinor<SType> a,
	    ATOOLS::CSpinor<SType> b);
    ATOOLS::CSpinor<SType> 
    Lorentz(const ATOOLS::CSpinor<SType> &a,
	    const ATOOLS::CScalar<SType> &b);
    
    virtual void Evaluate(const ATOOLS::CSpinor<SType> &a,
			  const ATOOLS::CSpinor<SType> &b) = 0;
    virtual void Evaluate(const ATOOLS::CSpinor<SType> &a,
			  const ATOOLS::CScalar<SType> &b) = 0;
    
  public:
    
    FFS_Vertex(const Vertex_Key &key,
	       const size_t &oew=0,const size_t &oqcd=0);
    
    void Evaluate();

  };// end of class FFS_Vertex

  template <typename SType>
  class VVS_Vertex: public Vertex_Base {
  public:

    typedef ATOOLS::CVec4<SType>   CVec4Type;
    typedef std::vector<CVec4Type> CVec4Type_Vector;

    typedef ATOOLS::CScalar<SType>   CScalarType;
    typedef std::vector<CScalarType> CScalarType_Vector;

  protected:

    int m_dir;
    
    ATOOLS::CScalar<SType> 
    Lorentz(const ATOOLS::CVec4<SType> &a,
	    const ATOOLS::CVec4<SType> &b);
    ATOOLS::CVec4<SType> 
    Lorentz(const ATOOLS::CVec4<SType> &a,
	    const ATOOLS::CScalar<SType> &b);

    virtual void Evaluate(const ATOOLS::CVec4<SType> &a,
			  const ATOOLS::CVec4<SType> &b) = 0;
    virtual void Evaluate(const ATOOLS::CVec4<SType> &a,
			  const ATOOLS::CScalar<SType> &b) = 0;

  public:

    VVS_Vertex(const Vertex_Key &key,
	       const size_t &oew=0,const size_t &oqcd=0);

    void Evaluate();

  };// end of class VVS_Vertex

  template <typename SType>
  class FFV_Vertex: public Vertex_Base {
  public:

    typedef std::complex<SType> SComplex;

    typedef ATOOLS::CSpinor<SType>   CSpinorType;
    typedef std::vector<CSpinorType> CSpinorType_Vector;

    typedef ATOOLS::CVec4<SType>   CVec4Type;
    typedef std::vector<CVec4Type> CVec4Type_Vector;

  protected:

    int m_dir;

    ATOOLS::CVec4<SType> 
    LorentzLeft(ATOOLS::CSpinor<SType> a,
		ATOOLS::CSpinor<SType> b);
    ATOOLS::CVec4<SType> 
    LorentzRight(ATOOLS::CSpinor<SType> a,
		 ATOOLS::CSpinor<SType> b);

    inline bool CalcLeft(const ATOOLS::CSpinor<SType> &a,
			 const ATOOLS::CSpinor<SType> &b) 
    { return a.On()&2 && b.On()&1; }
    inline bool CalcRight(const ATOOLS::CSpinor<SType> &a,
			  const ATOOLS::CSpinor<SType> &b) 
    { return a.On()&1 && b.On()&2; }

    ATOOLS::CSpinor<SType> 
    LorentzLeft(const ATOOLS::CSpinor<SType> &a,
		const ATOOLS::CVec4<SType> &b);
    ATOOLS::CSpinor<SType> 
    LorentzRight(const ATOOLS::CSpinor<SType> &a,
		 const ATOOLS::CVec4<SType> &b);

    inline bool CalcLeft(const ATOOLS::CSpinor<SType> &a) 
    { return (a.B()>0 && a.On()&1) || (a.B()<0 && a.On()&2); }
    inline bool CalcRight(const ATOOLS::CSpinor<SType> &a) 
    { return (a.B()>0 && a.On()&2) || (a.B()<0 && a.On()&1); }

    virtual void Evaluate(const ATOOLS::CSpinor<SType> &a,
			  const ATOOLS::CSpinor<SType> &b) = 0;
    virtual void Evaluate(const ATOOLS::CSpinor<SType> &a,
			  const ATOOLS::CVec4<SType> &b) = 0;

  public:

    FFV_Vertex(const Vertex_Key &key,
	       const size_t &oew=0,const size_t &oqcd=0);

    void Evaluate();

  };// end of class FFV_Vertex

  template <typename SType>
  class VVV_Vertex: public Vertex_Base {
  public:

    typedef ATOOLS::CVec4<SType>   CVec4Type;
    typedef std::vector<CVec4Type> CVec4Type_Vector;

  protected:

    ATOOLS::CVec4<SType> 
    Lorentz(const ATOOLS::CVec4<SType> &a,
	    const ATOOLS::CVec4<SType> &b);

    virtual void Evaluate(const ATOOLS::CVec4<SType> &a,
			  const ATOOLS::CVec4<SType> &b) = 0;

  public:

    VVV_Vertex(const Vertex_Key &key,
	       const size_t &oew=0,const size_t &oqcd=0);

    void Evaluate();

  };// end of class VVV_Vertex

  template <typename SType>
  class VVT_Vertex: public Vertex_Base {
  public:

    typedef ATOOLS::CVec4<SType>   CVec4Type;
    typedef std::vector<CVec4Type> CVec4Type_Vector;

    typedef ATOOLS::CAsT4<SType>   CAsT4Type;
    typedef std::vector<CAsT4Type> CAsT4Type_Vector;

  protected:

    int m_dir;

    ATOOLS::CAsT4<SType> 
    Lorentz(const ATOOLS::CVec4<SType> &a,
	    const ATOOLS::CVec4<SType> &b);
    ATOOLS::CVec4<SType> 
    Lorentz(const ATOOLS::CVec4<SType> &a,
	    const ATOOLS::CAsT4<SType> &b);

    virtual void Evaluate(const ATOOLS::CVec4<SType> &a,
			  const ATOOLS::CVec4<SType> &b) = 0;
    virtual void Evaluate(const ATOOLS::CVec4<SType> &a,
			  const ATOOLS::CAsT4<SType> &b) = 0;

  public:

    VVT_Vertex(const Vertex_Key &key,
	       const size_t &oew=0,const size_t &oqcd=0);

    void Evaluate();

  };// end of class VVT_Vertex

}// end of namespace COMIX

#define DECLARE_VTX_GETTER(NAME)					\
  									\
  template <typename SType,char STag>					\
  class NAME: public ATOOLS::Getter_Function				\
  <Vertex_Base,Vertex_Key,std::less<std::string> > {			\
  protected:								\
    void PrintInfo(std::ostream &str,const size_t width) const;		\
    Object_Type *							\
      operator()(const Parameter_Type &parameters) const;		\
  public:								\
    NAME(const std::string &name):					\
      ATOOLS::Getter_Function<Vertex_Base,Vertex_Key,			\
      std::less<std::string> >(STag+name) {}				\
  }

#endif
