#include "COMIX/Vertices/Vertex.H"

#include "ATOOLS/Org/Exception.H"

using namespace COMIX;
using namespace ATOOLS;

template <typename SType>
SSS_Vertex<SType>::SSS_Vertex(const Vertex_Key &key,
			      const size_t &oew,const size_t &oqcd):
  Vertex_Base(key,oew,oqcd) {}

template <typename SType>
void SSS_Vertex<SType>::Evaluate()
{
  m_zero=true;
  if (p_a->Zero()||p_b->Zero()) return;
#ifdef DEBUG__BG
  msg_Debugging()<<*p_a<<"(+)"<<*p_b<<" SSS\n";
  msg_Indent();
#endif
  const CScalarType_Vector &ca(p_a->Current_Base::J<CScalarType>());
  const CScalarType_Vector &cb(p_b->Current_Base::J<CScalarType>());
  for (typename CScalarType_Vector::const_iterator 
	 ait(ca.begin());ait!=ca.end();++ait)
    for (typename CScalarType_Vector::const_iterator 
	   bit(cb.begin());bit!=cb.end();++bit)
      Evaluate(*ait,*bit);
}

template <typename SType> CScalar<SType> 
SSS_Vertex<SType>::Lorentz(CScalarType a,CScalarType b)
{
#ifdef DEBUG__BG
  msg_Debugging()<<"    "<<a<<"\n";
  msg_Debugging()<<"    "<<b<<"\n";
#endif
  CScalarType j(a*b);
  j.SetH(a.H(0)+b.H(0),0);
  j.SetH(a.H(1)+b.H(1),1);
  m_zero=false;
  return j;
}

template <typename SType>
FFS_Vertex<SType>::FFS_Vertex(const Vertex_Key &key,
			      const size_t &oew,const size_t &oqcd):
  Vertex_Base(key,oew,oqcd), 
  m_dir(key.p_b->Flav().IsFermion()?
	(key.p_a->Flav().IsFermion()?0:2):1) {}

template <typename SType>
void FFS_Vertex<SType>::Evaluate()
{
  m_zero=true;
  if (p_a->Zero()||p_b->Zero()) return;
#ifdef DEBUG__BG
  msg_Debugging()<<*p_a<<"(+)"<<*p_b<<" FFS("<<m_dir<<")\n";
  msg_Indent();
#endif
  switch (m_dir) {
  case 0: {
    const CSpinorType_Vector &ca(p_a->Current_Base::J<CSpinorType>());
    const CSpinorType_Vector &cb(p_b->Current_Base::J<CSpinorType>());
    for (typename CSpinorType_Vector::const_iterator 
	   ait(ca.begin());ait!=ca.end();++ait)
      for (typename CSpinorType_Vector::const_iterator 
	     bit(cb.begin());bit!=cb.end();++bit)
	Evaluate(*ait,*bit);
    break;
  }
  case 1: {
    const CSpinorType_Vector &ca(p_a->Current_Base::J<CSpinorType>());
    const CScalarType_Vector &cb(p_b->Current_Base::J<CScalarType>());
    for (typename CSpinorType_Vector::const_iterator 
	   ait(ca.begin());ait!=ca.end();++ait)
      for (typename CScalarType_Vector::const_iterator 
	     bit(cb.begin());bit!=cb.end();++bit)
	Evaluate(*ait,*bit);
    break;
  }
  case 2: {
    const CScalarType_Vector &ca(p_a->Current_Base::J<CScalarType>());
    const CSpinorType_Vector &cb(p_b->Current_Base::J<CSpinorType>());
    for (typename CScalarType_Vector::const_iterator 
	   ait(ca.begin());ait!=ca.end();++ait)
      for (typename CSpinorType_Vector::const_iterator 
	     bit(cb.begin());bit!=cb.end();++bit)
	Evaluate(*bit,*ait);
    break;
  }
  default:
    THROW(fatal_error,"Internal error");
  }
}

template <typename SType> CScalar<SType> 
FFS_Vertex<SType>::Lorentz(CSpinorType a,CSpinorType b)
{
#ifdef DEBUG__BG
  msg_Debugging()<<"<>  "<<a<<"\n";
  msg_Debugging()<<"    "<<b<<"\n";
#endif
  CScalarType j(a*b);
  j.SetH(a.H(0)+b.H(0),0);
  j.SetH(a.H(1)+b.H(1),1);
  m_zero=false;
  return j;
}

template <typename SType> CSpinor<SType>
FFS_Vertex<SType>::Lorentz(const CSpinorType &a,const CScalarType &b)
{
#ifdef DEBUG__BG
  msg_Debugging()<<"    "<<a<<"\n";
  msg_Debugging()<<"    "<<b<<"\n";
#endif
  CSpinorType j(a*b[0]);
  j.SetH(a.H(0)+b.H(0),0);
  j.SetH(a.H(1)+b.H(1),1);
  m_zero=false;
  return j;
}

template <typename SType>
VVS_Vertex<SType>::VVS_Vertex(const Vertex_Key &key,
			      const size_t &oew,const size_t &oqcd):
  Vertex_Base(key,oew,oqcd), 
  m_dir(key.p_b->Flav().IsVector()?
	(key.p_a->Flav().IsVector()?0:2):1) {}
 
template <typename SType>
void VVS_Vertex<SType>::Evaluate()
{
  m_zero=true;
  if (p_a->Zero()||p_b->Zero()) return;
#ifdef DEBUG__BG
  msg_Debugging()<<*p_a<<"(+)"<<*p_b<<" VVS("<<m_dir<<")\n";
  msg_Indent();
#endif
  switch(m_dir) {
  case 0: {
    const CVec4Type_Vector &ca(p_a->Current_Base::J<CVec4Type>());
    const CVec4Type_Vector &cb(p_b->Current_Base::J<CVec4Type>());
    for (typename CVec4Type_Vector::const_iterator 
	   ait(ca.begin());ait!=ca.end();++ait)
      for (typename CVec4Type_Vector::const_iterator 
	     bit(cb.begin());bit!=cb.end();++bit)
	Evaluate(*ait,*bit);
    break;
  }
  case 1: {
    const CVec4Type_Vector &ca(p_a->Current_Base::J<CVec4Type>());
    const CScalarType_Vector &cb(p_b->Current_Base::J<CScalarType>());
    for (typename CVec4Type_Vector::const_iterator 
	   ait(ca.begin());ait!=ca.end();++ait)
      for (typename CScalarType_Vector::const_iterator 
	     bit(cb.begin());bit!=cb.end();++bit)
	Evaluate(*ait,*bit);
    break;
  }
  case 2: {
    const CScalarType_Vector &ca(p_a->Current_Base::J<CScalarType>());
    const CVec4Type_Vector &cb(p_b->Current_Base::J<CVec4Type>());
    for (typename CScalarType_Vector::const_iterator 
	   ait(ca.begin());ait!=ca.end();++ait)
      for (typename CVec4Type_Vector::const_iterator 
	     bit(cb.begin());bit!=cb.end();++bit)
	Evaluate(*bit,*ait);
    break;
  }
  default:
    THROW(fatal_error,"Internal error");
  }
}
 
template <typename SType> CScalar<SType>
VVS_Vertex<SType>::Lorentz(const CVec4Type &a,
			   const CVec4Type &b)
{
#ifdef DEBUG__BG
  msg_Debugging()<<"    "<<a<<"\n";
  msg_Debugging()<<"    "<<b<<"\n";
#endif
  CScalarType j(a*b);
  j.SetH(a.H(0)+b.H(0),0);
  j.SetH(a.H(1)+b.H(1),1);
  m_zero=false;
  return j;
}

template <typename SType> CVec4<SType>
VVS_Vertex<SType>::Lorentz(const CVec4Type &a,
			   const CScalarType &b)
{
#ifdef DEBUG__BG
  msg_Debugging()<<"    "<<a<<"\n";
  msg_Debugging()<<"    "<<b<<"\n";
#endif
  CVec4Type j(a*b[0]);
  j.SetH(a.H(0)+b.H(0),0);
  j.SetH(a.H(1)+b.H(1),1);
  m_zero=false;
  return j;
}

template <typename SType>
FFV_Vertex<SType>::FFV_Vertex(const Vertex_Key &key,
			      const size_t &oew,const size_t &oqcd):
  Vertex_Base(key,oew,oqcd), 
  m_dir(key.p_b->Flav().IsFermion()?
	(key.p_a->Flav().IsFermion()?0:2):1) {}

template <typename SType>
void FFV_Vertex<SType>::Evaluate()
{
  m_zero=true;
  if (p_a->Zero()||p_b->Zero()) return;
#ifdef DEBUG__BG
  msg_Debugging()<<*p_a<<"(+)"<<*p_b<<" FFV("<<m_dir<<")\n";
  msg_Indent();
#endif
  switch (m_dir) {
  case 0: {
    const CSpinorType_Vector &ca(p_a->Current_Base::J<CSpinorType>());
    const CSpinorType_Vector &cb(p_b->Current_Base::J<CSpinorType>());
    for (typename CSpinorType_Vector::const_iterator 
	   ait(ca.begin());ait!=ca.end();++ait)
      for (typename CSpinorType_Vector::const_iterator 
	     bit(cb.begin());bit!=cb.end();++bit)
	Evaluate(*ait,*bit);
    break;
  }
  case 1: {
    const CSpinorType_Vector &ca(p_a->Current_Base::J<CSpinorType>());
    const CVec4Type_Vector &cb(p_b->Current_Base::J<CVec4Type>());
    for (typename CSpinorType_Vector::const_iterator 
	   ait(ca.begin());ait!=ca.end();++ait)
      for (typename CVec4Type_Vector::const_iterator 
	     bit(cb.begin());bit!=cb.end();++bit)
	Evaluate(*ait,*bit);
    break;
  }
  case 2: {
    const CVec4Type_Vector &ca(p_a->Current_Base::J<CVec4Type>());
    const CSpinorType_Vector &cb(p_b->Current_Base::J<CSpinorType>());
    for (typename CVec4Type_Vector::const_iterator 
	   ait(ca.begin());ait!=ca.end();++ait)
      for (typename CSpinorType_Vector::const_iterator 
	     bit(cb.begin());bit!=cb.end();++bit)
	Evaluate(*bit,*ait);
    break;
  }
  default:
    THROW(fatal_error,"Internal error");
  }
}

template <typename SType> CVec4<SType> 
FFV_Vertex<SType>::LorentzLeft(CSpinorType a,CSpinorType b)
{
  switch (a.B()) {
  case -1: break;
  case 1: std::swap<CSpinorType>(a,b); break;
  default:
    THROW(fatal_error,"Internal error");
  }
#ifdef DEBUG__BG
  msg_Debugging()<<"<>  "<<a<<"\n";
  msg_Debugging()<<"    "<<b<<"\n";
#endif
  SComplex j01(a[3]*b[1]), j02(a[2]*b[0]);
  SComplex j11(-a[2]*b[1]), j12(-a[3]*b[0]);
  CVec4Type j(j01+j02,j11+j12,SComplex(0.0,-1.0)*(j11-j12),j01-j02,
	      0,0,a.H(0)+b.H(0),a.H(1)+b.H(1));
  m_zero=false;
  return j;
}

template <typename SType> CVec4<SType>
FFV_Vertex<SType>::LorentzRight(CSpinorType a,CSpinorType b)
{
  switch (a.B()) {
  case -1: break;
  case 1: std::swap<CSpinorType>(a,b); break;
  default:
    THROW(fatal_error,"Internal error");
  }
#ifdef DEBUG__BG
  msg_Debugging()<<"<>  "<<a<<"\n";
  msg_Debugging()<<"    "<<b<<"\n";
#endif
  SComplex j01(a[0]*b[2]), j02(a[1]*b[3]);
  SComplex j11(a[0]*b[3]), j12(a[1]*b[2]);
  CVec4Type j(j01+j02,j11+j12,SComplex(0.0,-1.0)*(j11-j12),j01-j02,
	      0,0,a.H(0)+b.H(0),a.H(1)+b.H(1));
  m_zero=false;
  return j;
}

template <typename SType> CSpinor<SType>
FFV_Vertex<SType>::LorentzLeft(const CSpinorType &a,
			       const CVec4Type &b)
{
  switch (a.B()) {
  case -1: {
#ifdef DEBUG__BG
    msg_Debugging()<<"<|g "<<a<<"\n";
    msg_Debugging()<<"    "<<b<<"\n";
#endif
    CSpinorType j(a.R(),a.B(),a(),a.H(0)+b.H(0),a.H(1)+b.H(1),1);
    SComplex jp(a.PPlus(b)), jm(a.PMinus(b)), jt(a.PT(b)), jtc(a.PTC(b));
    j[0]=a[2]*jp+a[3]*jt;
    j[1]=a[2]*jtc+a[3]*jm;
    j[3]=j[2]=SComplex(0.0,0.0);
    m_zero=false;
    return j;
  }
  case 1: {
#ifdef DEBUG__BG
    msg_Debugging()<<"g|> "<<a<<"\n";
    msg_Debugging()<<"    "<<b<<"\n";
#endif
    CSpinorType j(a.R(),a.B(),a(),a.H(0)+b.H(0),a.H(1)+b.H(1),2);
    SComplex jp(a.PPlus(b)), jm(a.PMinus(b)), jt(a.PT(b)), jtc(a.PTC(b));
    j[1]=j[0]=SComplex(0.0,0.0);
    j[2]=a[0]*jp+a[1]*jtc;
    j[3]=a[0]*jt+a[1]*jm;
    m_zero=false;
    return j;
  }
  default:
    THROW(fatal_error,"Internal error");
  }
  return CSpinorType();
}

template <typename SType> CSpinor<SType>
FFV_Vertex<SType>::LorentzRight(const CSpinorType &a,
				const CVec4Type &b)
{
  switch (a.B()) {
  case -1: {
#ifdef DEBUG__BG
    msg_Debugging()<<"<|g "<<a<<"\n";
    msg_Debugging()<<"    "<<b<<"\n";
#endif
    CSpinorType j(a.R(),a.B(),a(),a.H(0)+b.H(0),a.H(1)+b.H(1),2);
    SComplex jp(a.PPlus(b)), jm(a.PMinus(b)), jt(a.PT(b)), jtc(a.PTC(b));
    j[1]=j[0]=SComplex(0.0,0.0);
    j[2]=a[0]*jm-a[1]*jt;
    j[3]=-a[0]*jtc+a[1]*jp;
    m_zero=false;
    return j;
  }
  case 1: {
#ifdef DEBUG__BG
    msg_Debugging()<<"g|> "<<a<<"\n";
    msg_Debugging()<<"    "<<b<<"\n";
#endif
    CSpinorType j(a.R(),a.B(),a(),a.H(0)+b.H(0),a.H(1)+b.H(1),1);
    SComplex jp(a.PPlus(b)), jm(a.PMinus(b)), jt(a.PT(b)), jtc(a.PTC(b));
    j[0]=a[2]*jm-a[3]*jtc;
    j[1]=-a[2]*jt+a[3]*jp;
    j[3]=j[2]=SComplex(0.0,0.0);
    m_zero=false;
    return j;
  }
  default:
    THROW(fatal_error,"Internal error");
  }
  return CSpinorType();
}

template <typename SType>
VVV_Vertex<SType>::VVV_Vertex(const Vertex_Key &key,
			      const size_t &oew,const size_t &oqcd): 
  Vertex_Base(key,oew,oqcd) {}

template <typename SType>
void VVV_Vertex<SType>::Evaluate()
{
  m_zero=true;
  if (p_a->Zero()||p_b->Zero()) return;
#ifdef DEBUG__BG
  msg_Debugging()<<*p_a<<"(+)"<<*p_b<<" VVV\n";
  msg_Indent();
#endif
  const CVec4Type_Vector &ca(p_a->Current_Base::J<CVec4Type>());
  const CVec4Type_Vector &cb(p_b->Current_Base::J<CVec4Type>());
  for (typename CVec4Type_Vector::const_iterator 
	 ait(ca.begin());ait!=ca.end();++ait)
    for (typename CVec4Type_Vector::const_iterator 
	   bit(cb.begin());bit!=cb.end();++bit)
      Evaluate(*ait,*bit);
}

template <typename SType> CVec4<SType>
VVV_Vertex<SType>::Lorentz(const CVec4Type &a,
			   const CVec4Type &b)
{
  CVec4Type j((a*b)*CVec4Type(p_a->P()-p_b->P())
	      +(a*CVec4Type(2.0*p_b->P()+p_a->P()))*b
	      -(b*CVec4Type(2.0*p_a->P()+p_b->P()))*a);
  j.SetH(a.H(0)+b.H(0),0);
  j.SetH(a.H(1)+b.H(1),1);
  m_zero=false;
  return j;
}

template <typename SType>
VVT_Vertex<SType>::VVT_Vertex(const Vertex_Key &key,
			      const size_t &oew,const size_t &oqcd):
  Vertex_Base(key,oew,oqcd), 
  m_dir(key.p_b->Flav().IsVector()?
	(key.p_a->Flav().IsVector()?0:2):1) {}

template <typename SType>
void VVT_Vertex<SType>::Evaluate()
{
  m_zero=true;
  if (p_a->Zero()||p_b->Zero()) return;
#ifdef DEBUG__BG
  msg_Debugging()<<*p_a<<"(+)"<<*p_b<<" VVT("<<m_dir<<")\n";
  msg_Indent();
#endif
  switch (m_dir) {
  case 0: {
    const CVec4Type_Vector &ca(p_a->Current_Base::J<CVec4Type>());
    const CVec4Type_Vector &cb(p_b->Current_Base::J<CVec4Type>());
    for (typename CVec4Type_Vector::const_iterator 
	   ait(ca.begin());ait!=ca.end();++ait) 
      for (typename CVec4Type_Vector::const_iterator 
	     bit(cb.begin());bit!=cb.end();++bit)
	Evaluate(*ait,*bit);
    break;
  }
  case 1: {
    const CVec4Type_Vector &ca(p_a->Current_Base::J<CVec4Type>());
    const CAsT4Type_Vector &cb(p_b->Current_Base::J<CAsT4Type>());
    for (typename CVec4Type_Vector::const_iterator 
	   ait(ca.begin());ait!=ca.end();++ait)
      for (typename CAsT4Type_Vector::const_iterator 
	     bit(cb.begin());bit!=cb.end();++bit)
	Evaluate(*ait,*bit);
    break;
  }
  case 2: {
    const CAsT4Type_Vector &ca(p_a->Current_Base::J<CAsT4Type>());
    const CVec4Type_Vector &cb(p_b->Current_Base::J<CVec4Type>());
    for (typename CAsT4Type_Vector::const_iterator 
	   ait(ca.begin());ait!=ca.end();++ait)
      for (typename CVec4Type_Vector::const_iterator 
	     bit(cb.begin());bit!=cb.end();++bit)
	Evaluate(*bit,*ait);
    break;
  }
  default:
    THROW(fatal_error,"Internal error");
  }
}

template <typename SType> CAsT4<SType>
VVT_Vertex<SType>::Lorentz(const CVec4Type &a,const CVec4Type &b)
{
  m_zero=false;
  return CAsT4Type(a,b);
}

template <typename SType> CVec4<SType>
VVT_Vertex<SType>::Lorentz(const CVec4Type &a,const CAsT4Type &b)
{
  CVec4Type j(a*b);
  j.SetH(a.H(0)+b.H(0),0);
  j.SetH(a.H(1)+b.H(1),1);
  m_zero=false;
  return j;
}

namespace COMIX {

  template class SSS_Vertex<double>;
  template class FFS_Vertex<double>;
  template class VVS_Vertex<double>;
  template class FFV_Vertex<double>;
  template class VVV_Vertex<double>;
  template class VVT_Vertex<double>;

  template class SSS_Vertex<long double>;
  template class FFS_Vertex<long double>;
  template class VVS_Vertex<long double>;
  template class FFV_Vertex<long double>;
  template class VVV_Vertex<long double>;
  template class VVT_Vertex<long double>;

}
