#include "METOOLS/Explicit/Current.H"
#include "METOOLS/Currents/C_Scalar.H"

namespace METOOLS {

  template <typename SType>
  class CS: public Current,
	    public Current_Contractor<SType> {
  public:

    typedef std::complex<SType>   SComplex;
    typedef std::vector<SComplex> SComplex_Vector;

    typedef CScalar<SType> CScalarType;
    typedef std::vector<CScalarType*> CScalarType_Vector;

  protected:

    bool m_pseudo;

    SComplex m_cmass2;

    std::string CLabel() const;

  public:

    CS(const Current_Key &key);

    void ConstructJ(const ATOOLS::Vec4D &p,const int ch,
		    const int cr,const int ca);
    void SetGauge(const ATOOLS::Vec4D &k);

    void AddPropagator();

    void SContract
    (const Current &c,const Int_Vector &pols,
     SComplex_Vector &ress,const size_t &offset) const;

    std::string Format(const CObject *c) const;

    char Type() const;    

  };// end of class CS

}// end of namespace METOOLS

#include "METOOLS/Explicit/Vertex.H"
#include "METOOLS/Explicit/Dipole_Kinematics.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/STL_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"

#define M_I SComplex(0.0,1.0)

using namespace METOOLS;
using namespace ATOOLS;

template <typename SType>
CS<SType>::CS(const Current_Key &key): 
  Current(key), m_pseudo(false)
{
  switch (key.m_fl.Kfcode()) {
  case kf_h0_qsc:
    m_pseudo=true;
    break;
  default:
    m_cmass2=SComplex(sqr(this->m_mass),-this->m_mass*this->m_width);
    break;
  }
}

template <typename SType>
void CS<SType>::ConstructJ(const ATOOLS::Vec4D &p,const int ch,
			   const int cr,const int ca)
{
  this->m_p=p;
  this->ResetJ();
  if (ch==0) {
    CScalarType *j(CScalarType::New(CScalarType(1.0,0,0)));
#ifdef DEBUG__BG
    msg_Debugging()<<METHOD<<"(): '+' "<<this->m_id<<" "<<*j
		   <<" "<<this->m_fl<<", m = "<<p.Mass()<<"\n";
#endif
    AddJ(j);
  }
}

template <typename SType>
void CS<SType>::SetGauge(const ATOOLS::Vec4D &k)
{
}

template <typename SType>
void CS<SType>::AddPropagator()
{
  // add propagator for off-shell leg
  SComplex prop(m_pseudo?M_I:M_I/(SType(this->m_p.Abs2())-m_cmass2));
  if (this->m_osd) prop=SComplex(M_I);
#ifdef DEBUG__BG
  msg_Debugging()<<"propagator: "<<prop<<" <- p^2 = "
		 <<this->m_p.Abs2()<<", m = "<<sqrt(m_cmass2)<<"\n";
#endif
  for (size_t i(0);i<m_j.size();++i) {
  CScalarType_Vector *j(m_j[i].Get<CScalarType>());
  for (typename CScalarType_Vector::iterator 
	 jit(j->begin());jit!=j->end();++jit) **jit*=prop;
  }
}

template <typename SType> void CS<SType>::SContract
(const Current &c,const Int_Vector &pols,
 SComplex_Vector &ress,const size_t &offset) const
{
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"(): {\n";
  msg_Indent();
#endif
  double phase(0.0);
  const std::vector<size_t> *pm(NULL);
  if (p_sub) {
    Vertex *v(p_sub->Sub()->In().front());
    if (v->Info()->Mode()==1) {
      phase=v->Kin()->Phase(offset==1?0:1);
      pm=&v->Kin()->PM();
    }
  }
  if (c.Type()!='S') THROW(fatal_error,"Invalid current type.");
  size_t i(0);
  for (typename CObject_Matrix::const_iterator 
	 ajit1(m_j.begin());ajit1!=m_j.end();++ajit1) {	
    const CScalarType_Vector *j(ajit1->Get<CScalarType>());
    for (typename CObject_Matrix::const_iterator 
	   ajit2(c.J().begin());ajit2!=c.J().end();++ajit2,++i) {
      // if (!pols[i]) continue;
      const CScalarType_Vector *cj(ajit2->Get<CScalarType>());
      for (typename CScalarType_Vector::const_iterator 
	     jit2(cj->begin());jit2!=cj->end();++jit2) 
	for (typename CScalarType_Vector::const_iterator 
	       jit1(j->begin());jit1!=j->end();++jit1)
	  if ((**jit1)(0)==(**jit2)(1) && (**jit1)(1)==(**jit2)(0) &&
	      (*jit1)->S()==offset && (*jit2)->S()==offset) {
#ifdef DEBUG__BG
	    msg_Debugging()<<"Add ("<<m_hm[i]<<")"
			   <<**jit1***jit2<<"["<<offset<<"]\n";
#endif
	    ress[m_hm[i]]+=**jit1***jit2;
	    if (offset && pm) {
#ifdef DEBUG__BG
	      msg_Debugging()<<"Add ("<<(*pm)[m_hm[i]]<<")"
			     <<**jit1***jit2*SType(phase)<<" ["
			     <<offset<<"] ( phase = "<<phase<<" )\n";
#endif
	      ress[(*pm)[m_hm[i]]]+=**jit1***jit2*SType(phase);
	    }
	  }
    }
  }
#ifdef DEBUG__BG
  msg_Debugging()<<"}\n";
#endif
}

template <typename SType>
std::string CS<SType>::Format(const CObject *c) const
{
  return ToString(*(CScalarType*)c,6);
}

template <typename SType>
char CS<SType>::Type() const
{
  return 'S';
}

template <typename SType>
std::string CS<SType>::CLabel() const
{
  return "dashes,label.side=right,label=$"+
    (this->m_out.empty()?this->m_fl.Bar():this->m_fl).TexName()+"$";
}

DECLARE_GETTER(CS<double>,"DS",Current,Current_Key);

Current *ATOOLS::Getter<Current,Current_Key,CS<double> >::
operator()(const Current_Key &key) const
{
  if (key.m_fl.IsScalar()) return new CS<double>(key);
  return NULL;
}

void ATOOLS::Getter<Current,Current_Key,CS<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"scalar current (double)";
}
