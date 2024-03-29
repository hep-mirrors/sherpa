#include "METOOLS/Explicit/Current.H"
#include "METOOLS/Currents/C_Vector.H"
#include "ATOOLS/Phys/Spinor.H"
/*!
  @file V_C.C
  @brief Implements the class CV.
*/

namespace METOOLS {

  template <typename SType>
  class CV: public Current,
	    public Current_Contractor<SType> {
  public:

    typedef std::complex<SType>   SComplex;
    typedef std::vector<SComplex> SComplex_Vector;

    typedef ATOOLS::Spinor<SType> SpinorType;
    typedef ATOOLS::Vec4<SType>   Vec4Type;

    typedef CVec4<SType> CVec4Type;
    typedef std::vector<CVec4Type*> CVec4Type_Vector;

  protected:

    /// @name (Squared) mass of the propagator
    /// @{
    SComplex m_cmass2, m_cmass;
    /// @}

    /// @brief returns the label for the graph output
    std::string CLabel() const;

  private:

    /// @brief gauge vector
    ATOOLS::Vec4D m_k;

    /// @brief gauge spinors, @brief gauge spinors
    SpinorType    m_kp, m_km;

    /// @brief constructs a vector from two spinors
    CVec4Type VT(const SpinorType &a,const SpinorType &b);

    /// @name massless minus and plus polarisation vectors
    /// @{
    CVec4Type EM(const ATOOLS::Vec4D &p,const int cr,const int ca);
    CVec4Type EP(const ATOOLS::Vec4D &p,const int cr,const int ca);
    /// @}

    /// @name massive minus, plus, and longitudinal polarisation vectors
    /// @{
    CVec4Type EMM(const ATOOLS::Vec4D &p,const int cr,const int ca);
    CVec4Type EMP(const ATOOLS::Vec4D &p,const int cr,const int ca);
    CVec4Type EML(const ATOOLS::Vec4D &p,const int cr,const int ca);
    /// @}

  public:

    /// @brief constructor using Current_Key key
    CV(const Current_Key &key);

    /// @brief constructs current
    void ConstructJ(const ATOOLS::Vec4D &p,const int ch,
		    const int cr,const int ca,const int mode);

    /// @brief sets the gauge vector
    void SetGauge(const ATOOLS::Vec4D &k);

    /// @brief adds the propagator to the current
    void AddPropagator();

    /// @brief contracts a current
    void SContract
    (const Current &c,const Int_Vector &pols,
     SComplex_Vector &ress,const size_t &offset) const;

    /// @brief formats the current into a string
    std::string Format(const CObject *c) const;

    /// @brief returns the current type
    char Type() const;    

  };// end of class CV
  /*!
    @class CV
    @brief This class represents a vector current.

    This class represents a vector current.
  */

}// end of namespace METOOLS

#include "METOOLS/Explicit/Vertex.H"
#include "METOOLS/Explicit/Color_Calculator.H"
#include "METOOLS/Explicit/Dipole_Kinematics.H"
#include "METOOLS/Explicit/Dipole_Color.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/STL_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"

#define M_I SComplex(0.0,1.0)

using namespace METOOLS;
using namespace ATOOLS;

template <typename SType>
CV<SType>::CV(const Current_Key &key): 
  Current(key), m_cmass2(0.0), m_cmass(0.0)
{
  m_cmass=sqrt(m_cmass2=SComplex(sqr(this->m_mass),-this->m_mass*this->m_width));
  if (key.m_n==1 && key.p_model->ScalarNumber("WidthScheme")!=1) 
    m_cmass=sqrt(m_cmass2=Complex(sqr(this->m_mass),0.0));
}

template <typename SType> CVec4<SType> 
CV<SType>::VT(const SpinorType &a,const SpinorType &b)
{
  CVec4Type e;
  e[0]=a.U1()*b.U1()+a.U2()*b.U2();
  e[SpinorType::R3()]=a.U1()*b.U1()-a.U2()*b.U2();
  e[SpinorType::R1()]=a.U1()*b.U2()+a.U2()*b.U1();
  e[SpinorType::R2()]=SComplex(0.0,1.0)*(a.U1()*b.U2()-a.U2()*b.U1());
  return e;
}

template <typename SType> CVec4<SType>
CV<SType>::EM(const Vec4D &p,const int cr,const int ca)
{
  SpinorType pp(1,p);
  CVec4Type e(VT(pp,m_km));
  e(0)=cr; e(1)=ca;
  e.SetH(1);
  static SType sqrttwo(sqrt(SType(2.0)));
  return e/(sqrttwo*std::conj(m_kp*pp));
}

template <typename SType> CVec4<SType>
CV<SType>::EP(const Vec4D &p,const int cr,const int ca)
{
  SpinorType pm(-1,p);
  CVec4Type e(VT(m_kp,pm));
  e(0)=cr; e(1)=ca;
  e.SetH(0);
  static SType sqrttwo(sqrt(SType(2.0)));
  return e/(sqrttwo*std::conj(m_km*pm));
}

template <typename SType> CVec4<SType>
CV<SType>::EMM(const Vec4D &p,const int cr,const int ca)
{
  return EM(p-p.Abs2()/(2.0*m_k*p)*m_k,cr,ca);
}

template <typename SType> CVec4<SType>
CV<SType>::EMP(const Vec4D &p,const int cr,const int ca)
{
  return EP(p-p.Abs2()/(2.0*m_k*p)*m_k,cr,ca);
}

template <typename SType> CVec4<SType>
CV<SType>::EML(const Vec4D &p,const int cr,const int ca)
{
  double p2(p.Abs2()), a(p2/(2.0*m_k*p));
  Vec4D b(p-a*m_k);
  SpinorType bm(-1,b), bp(1,b), am(-1,m_k), ap(1,m_k);
  CVec4Type e(VT(bp,bm)-SType(a)*VT(ap,am));
  e(0)=cr; e(1)=ca;
  e.SetH(2);
  return e/sqrt(SComplex(4.0*p2));
}

template <typename SType>
void CV<SType>::ConstructJ(const ATOOLS::Vec4D &p,const int ch,
			   const int cr,const int ca,const int mode)
{
  this->m_p=p;
  if (this->m_fl.Mass()==0.0 && p[1]==0.0 && p[2]==0.0)
    this->m_p[0]=this->m_p[0]<0.0?
      -std::abs(this->m_p[3]):std::abs(this->m_p[3]);
  this->ResetJ();
  if (ch>=0) {
    if (this->m_msv && (ch==0 || ch==3)) {
      CVec4Type j(EML(this->m_p,cr,ca));
      j=this->m_dir>0?-j:j.Conj();
      AddJ(CVec4Type::New(j));
#ifdef DEBUG__BG
      msg_Debugging()<<METHOD<<"(): "<<(this->m_dir>0?'I':'O')
		     <<"0 "<<this->m_id<<" "<<j
		     <<" "<<this->m_fl<<", m = "<<m_cmass<<"\n";
#endif
    }
    if (ch!=3) {
      CVec4Type j(this->m_msv?this->m_dir>0?
		  EMM(this->m_p,cr,ca):EMP(this->m_p,cr,ca):
		  this->m_dir>0?EM(this->m_p,cr,ca):EP(this->m_p,cr,ca));
      j=this->m_dir>0?j:j.Conj();
      CVec4Type *c(CVec4Type::New(j));
      AddJ(c);
#ifdef DEBUG__BG
      msg_Debugging()<<METHOD<<"(): "<<(this->m_dir>0?'I':'O')
		     <<"+ "<<this->m_id<<" "<<j
		     <<" "<<this->m_fl<<", m = "<<m_cmass<<"\n";
#endif
      if (p_sub) static_cast<Dipole_Color*>
	(p_sub->In().front()->Color().front())->AddJJK(c);
    }
  }
  if (ch<=0) {
    CVec4Type j(this->m_msv?this->m_dir>0?
		EMP(this->m_p,cr,ca):EMM(this->m_p,cr,ca):
		this->m_dir>0?EP(this->m_p,cr,ca):EM(this->m_p,cr,ca));
    j=this->m_dir>0?j:j.Conj();
    CVec4Type *c(CVec4Type::New(j));
    AddJ(c);
#ifdef DEBUG__BG
    msg_Debugging()<<METHOD<<"(): "<<(this->m_dir>0?'I':'O')
		   <<"- "<<this->m_id<<" "<<j
		   <<" "<<this->m_fl<<", m = "<<m_cmass<<"\n";
#endif
    if (p_sub) static_cast<Dipole_Color*>
      (p_sub->In().front()->Color().front())->AddJJK(c);
  }
#ifdef DEBUG__BG
  if (p_sub) Print();
#endif
}

template <typename SType>
void CV<SType>::SetGauge(const ATOOLS::Vec4D &k)
{
  m_k=k;
  m_kp=SpinorType(1,m_k);
  m_km=SpinorType(-1,m_k);
}

template <typename SType>
void CV<SType>::AddPropagator()
{
  // add propagator for off-shell leg
  SComplex p2(SType(this->m_p.Abs2())), prop(-M_I/(p2-m_cmass2));
  if (this->m_osd) prop=SComplex(M_I);
#ifdef DEBUG__BG
  msg_Debugging()<<"propagator: "<<prop<<"\n";
#endif
  for (size_t i(0);i<m_j.size();++i) {
  CVec4Type_Vector *j(m_j[i].template Get<CVec4Type>());
  if (!this->m_fl.IsGluon()) {
    if (!this->m_msv)
      for (typename CVec4Type_Vector::iterator 
	     jit(j->begin());jit!=j->end();++jit)
	**jit-=(**jit*Vec4Type(this->m_p))*CVec4Type(this->m_p)/p2;
    else
      for (typename CVec4Type_Vector::iterator 
	     jit(j->begin());jit!=j->end();++jit)
	**jit-=(**jit*Vec4Type(this->m_p))*CVec4Type(this->m_p)/m_cmass2;
  }
  for (typename CVec4Type_Vector::iterator 
	 jit(j->begin());jit!=j->end();++jit) **jit*=prop;
  }
}

template <typename SType> void CV<SType>::SContract
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
  if (c.Type()!='V') THROW(fatal_error,"Invalid current type.");
  size_t i(0);
  for (typename CObject_Matrix::const_iterator 
	 ajit1(m_j.begin());ajit1!=m_j.end();++ajit1) {	
    const CVec4Type_Vector *j(ajit1->Get<CVec4Type>());
    for (typename CObject_Matrix::const_iterator 
	   ajit2(c.J().begin());ajit2!=c.J().end();++ajit2,++i) {
      // if (!pols[i]) continue;
      const CVec4Type_Vector *cj(ajit2->Get<CVec4Type>());
      for (typename CVec4Type_Vector::const_iterator 
	     jit2(cj->begin());jit2!=cj->end();++jit2) 
	for (typename CVec4Type_Vector::const_iterator 
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
std::string CV<SType>::Format(const CObject *c) const
{
  return ToString(*(CVec4Type*)c,6);
}

template <typename SType>
char CV<SType>::Type() const
{
  return 'V';
}

template <typename SType>
std::string CV<SType>::CLabel() const
{
  switch (this->m_fl.Kfcode()) {
  case kf_gluon:
    return "gluon,label.side=right,label.dist=1.5curly_len,label=$g$";
  case kf_photon:
    return "photon,label.side=right,label.dist=1wiggly_len,label=$\\gamma$";
  case kf_Z:
    return "dots,label.side=right,label.dist=1wiggly_len,label=$Z^0$";
  case kf_Wplus:
    return "dots,label.side=right,label.dist=1wiggly_len,label=$"
      +(this->m_out.empty()?this->m_fl.Bar():this->m_fl).TexName()+"$";
  default: break;
  }
  return "wiggly,label.side=right,label.dist=1wiggly_len,label=$"
    +(this->m_out.empty()?this->m_fl.Bar():this->m_fl).TexName()+"$";
}

DECLARE_GETTER(CV<double>,"DV",Current,Current_Key);

Current *ATOOLS::Getter<Current,Current_Key,CV<double> >::
operator()(const Current_Key &key) const
{
  if (key.m_fl.IsVector()) return new CV<double>(key);
  return NULL;
}

void ATOOLS::Getter<Current,Current_Key,CV<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"vector current (double)";
}
