#include "COMIX/Currents/Current.H"
#include "ATOOLS/Phys/Spinor.H"

namespace COMIX {

  template <typename SType>
  class CS: public Current<ATOOLS::CScalar<SType> >,
	    public Current_Contractor_Base<SType> {
  public:

    typedef std::complex<SType>   SComplex;
    typedef std::vector<SComplex> SComplex_Vector;

    typedef ATOOLS::CScalar<SType>   CScalarType;
    typedef std::vector<CScalarType> CScalarType_Vector;

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

    SComplex SContract(const Current_Base &c,
		       const long unsigned int &hm,
		       const long unsigned int &hp) const;
    void SContract(const Current_Base &c,
		   const LongUIntMap_Matrix &hmp,
		   SComplex_Vector &ress) const;

    char Type() const;    

  };// end of class CS

  template <typename SType>
  class CF: public Current<ATOOLS::CSpinor<SType> >,
	    public Current_Contractor_Base<SType> {
  public:

    typedef std::complex<SType>   SComplex;
    typedef std::vector<SComplex> SComplex_Vector;

    typedef ATOOLS::CSpinor<SType>   CSpinorType;
    typedef std::vector<CSpinorType> CSpinorType_Vector;

  protected:

    SComplex m_cmass2, m_cmass;

    std::string CLabel() const;

  public:

    CF(const Current_Key &key);

    void ConstructJ(const ATOOLS::Vec4D &p,const int ch,
		    const int cr,const int ca);
    void SetGauge(const ATOOLS::Vec4D &k);

    void AddPropagator();

    SComplex SContract(const Current_Base &c,
		       const long unsigned int &hm,
		       const long unsigned int &hp) const;
    void SContract(const Current_Base &c,
		   const LongUIntMap_Matrix &hmp,
		   SComplex_Vector &ress) const;

    char Type() const;    

  };// end of class CF

  template <typename SType>
  class CV: public Current<ATOOLS::CVec4<SType> >,
	    public Current_Contractor_Base<SType> {
  public:

    typedef std::complex<SType>   SComplex;
    typedef std::vector<SComplex> SComplex_Vector;

    typedef ATOOLS::Spinor<SType> SpinorType;
    typedef ATOOLS::Vec4<SType>   Vec4Type;

    typedef ATOOLS::CVec4<SType>   CVec4Type;
    typedef std::vector<CVec4Type> CVec4Type_Vector;

  protected:

    SComplex m_cmass2, m_cmass;

    std::string CLabel() const;

  private:

    ATOOLS::Vec4D m_k;
    SpinorType    m_kp, m_km;

    CVec4Type VT(const SpinorType &a,const SpinorType &b);

    CVec4Type EM(const ATOOLS::Vec4D &p,const int cr,const int ca);
    CVec4Type EP(const ATOOLS::Vec4D &p,const int cr,const int ca);

    CVec4Type EMM(const ATOOLS::Vec4D &p,const int cr,const int ca);
    CVec4Type EMP(const ATOOLS::Vec4D &p,const int cr,const int ca);
    CVec4Type EML(const ATOOLS::Vec4D &p,const int cr,const int ca);

  public:

    CV(const Current_Key &key);

    void ConstructJ(const ATOOLS::Vec4D &p,const int ch,
		    const int cr,const int ca);
    void SetGauge(const ATOOLS::Vec4D &k);

    void AddPropagator();

    SComplex SContract(const Current_Base &c,
		       const long unsigned int &hm,
		       const long unsigned int &hp) const;
    void SContract(const Current_Base &c,
		   const LongUIntMap_Matrix &hmp,
		   SComplex_Vector &ress) const;

    char Type() const;    

  };// end of class CV

  template <typename SType>
  class CT: public Current<ATOOLS::CAsT4<SType> >,
	    public Current_Contractor_Base<SType> {
  public:

    typedef std::complex<SType>   SComplex;
    typedef std::vector<SComplex> SComplex_Vector;

    typedef ATOOLS::CAsT4<SType>   CAsT4Type;
    typedef std::vector<CAsT4Type> CAsT4Type_Vector;

  protected:

    SComplex m_prop;

    std::string CLabel() const;

  public:

    CT(const Current_Key &key);

    void ConstructJ(const ATOOLS::Vec4D &p,const int ch,
		    const int cr,const int ca);
    void SetGauge(const ATOOLS::Vec4D &k);

    void AddPropagator();

    SComplex SContract(const Current_Base &c,
		       const long unsigned int &hm,
		       const long unsigned int &hp) const;
    void SContract(const Current_Base &c,
		   const LongUIntMap_Matrix &hmp,
		   SComplex_Vector &ress) const;

    char Type() const;    

  };// end of class CT

}// end of namespace COMIX

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/STL_Tools.H"

#define M_I SComplex(0.0,1.0)

using namespace COMIX;
using namespace ATOOLS;

template <typename SType>
CS<SType>::CS(const Current_Key &key): 
  Current<CScalarType>(key), m_pseudo(false)
{
  switch (key.m_fl.Kfcode()) {
  case kf_h0:
    m_cmass2=SComplex(sqr(this->m_mass),-this->m_mass*this->m_width);
    break;
  case kf_h0_qsc:
    m_pseudo=true;
    break;
  default: break;
  }
}

template <typename SType>
void CS<SType>::ConstructJ(const ATOOLS::Vec4D &p,const int ch,
			   const int cr,const int ca)
{
  this->m_p=p;
  this->ResetJ();
  if (ch==0) {
    CScalarType j(CScalarType(1.0,0,0));
#ifdef DEBUG__BG
    msg_Debugging()<<METHOD<<"(): '+' "<<this->m_id<<" "<<j
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
#ifdef DEBUG__BG
  msg_Debugging()<<"propagator: "<<prop<<" <- p^2 = "
		 <<this->m_p.Abs2()<<", m = "<<sqrt(m_cmass2)<<"\n";
#endif
  for (typename CScalarType_Vector::iterator 
	 jit(this->m_j.begin());jit!=this->m_j.end();++jit) *jit*=prop;
}

template <typename SType> std::complex<SType> 
CS<SType>::SContract(const Current_Base &c,
		     const long unsigned int &hm,
		     const long unsigned int &hp) const
{
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"("<<hm<<","<<hp<<"): {\n";
  msg_Indent();
#endif
  if (c.Type()!='S') THROW(fatal_error,"Invalid current type.");
  Complex res(0.0);
  const CScalarType_Vector &cj(c.Current_Base::J<CScalarType>());
  for (typename CScalarType_Vector::const_iterator 
	 jit1(this->m_j.begin());jit1!=this->m_j.end();++jit1)
    for (typename CScalarType_Vector::const_iterator 
	   jit2(cj.begin());jit2!=cj.end();++jit2) 
      if (jit1->H(0)+jit2->H(0)==hm && jit1->H(1)+jit2->H(1)==hp) {
#ifdef DEBUG__BG
	msg_Debugging()<<"add "<<*jit1**jit2<<"\n";
	msg_Debugging()<<"    "<<*jit1<<"\n";
	msg_Debugging()<<"    "<<*jit2<<"\n";
#endif
	res+=*jit1**jit2;
      }
#ifdef DEBUG__BG
  msg_Debugging()<<"}\n";
#endif
  return res;
}

template <typename SType>
void CS<SType>::SContract(const Current_Base &c,
			  const LongUIntMap_Matrix &hmp,
			  SComplex_Vector &ress) const
{
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"(): {\n";
  msg_Indent();
#endif
  if (c.Type()!='S') THROW(fatal_error,"Invalid current type.");
  const CScalarType_Vector &cj(c.Current_Base::J<CScalarType>());
  for (typename CScalarType_Vector::const_iterator 
	 jit1(this->m_j.begin());jit1!=this->m_j.end();++jit1)
    for (typename CScalarType_Vector::const_iterator 
	   jit2(cj.begin());jit2!=cj.end();++jit2) 
      if ((*jit1)(0)==(*jit2)(1) && (*jit1)(1)==(*jit2)(0)) {
	LongUIntMap_Matrix::const_iterator 
	  hmit(hmp.find(jit1->H(0)+jit2->H(0)));
	if (hmit==hmp.end()) continue;
	LongUInt_Map::const_iterator 
	  hpit(hmit->second.find(jit1->H(1)+jit2->H(1)));
	if (hpit==hmit->second.end()) continue;
#ifdef DEBUG__BG
	msg_Debugging()<<"Add ("<<hmit->first<<","
		       <<hpit->first<<")"<<*jit1**jit2<<"\n";
#endif
	ress[hpit->second]+=*jit1**jit2;
      }
#ifdef DEBUG__BG
  msg_Debugging()<<"}\n";
#endif
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

DECLARE_TEMPLATE_GETTER(CS_Getter,"S",Current_Base,Current_Key);

template <typename SType,char STag> Current_Base *
CS_Getter<SType,STag>::operator()(const Current_Key &key) const
{
  if (key.m_fl.IsScalar()) return new CS<SType>(key);
  return NULL;
}

template <typename SType,char STag>
void CS_Getter<SType,STag>::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"scalar current "<<STag;
}

template <typename SType>
CF<SType>::CF(const Current_Key &key): 
  Current<CSpinorType>(key),
  m_cmass2(Complex(sqr(this->m_mass),-this->m_mass*this->m_width)), 
  m_cmass(sqrt(m_cmass2)) {}

template <typename SType>
void CF<SType>::ConstructJ(const ATOOLS::Vec4D &p,const int ch,
			   const int cr,const int ca)
{
  this->m_p=p;
  this->ResetJ();
  if (ch>=0) {
    CSpinorType j(this->m_fl.IsAnti()^(this->m_dir>0)?
		  CSpinorType(-1,-this->m_dir,1,p,
			      cr|ca,0,0,sqr(this->m_mass)):
		  CSpinorType(1,this->m_dir,1,p,
			      cr|ca,0,0,sqr(this->m_mass)));
    j.SetH(1<<this->m_key,this->m_dir>0?0:1);
#ifdef DEBUG__BG
    msg_Debugging()<<METHOD<<"(): "<<(this->m_dir>0?'I':'O')<<"+ "<<this->m_id
		   <<" "<<j<<" "<<(this->m_dir>0?this->m_fl.Bar():this->m_fl)
		   <<", m = "<<m_cmass<<" ("<<p.Mass()<<")\n";
#endif
    AddJ(j);
  }
  if (ch<=0) {
    CSpinorType j(this->m_fl.IsAnti()^(this->m_dir>0)?
		  CSpinorType(-1,-this->m_dir,-1,p,
			      cr|ca,0,0,sqr(this->m_mass)):
		  CSpinorType(1,this->m_dir,-1,p,
			      cr|ca,0,0,sqr(this->m_mass)));
    j.SetH(1<<this->m_key,this->m_dir>0?1:0);
#ifdef DEBUG__BG
    msg_Debugging()<<METHOD<<"(): "<<(this->m_dir>0?'I':'O')<<"- "<<this->m_id
		   <<" "<<j<<" "<<(this->m_dir>0?this->m_fl.Bar():this->m_fl)
		   <<", m = "<<m_cmass<<" ("<<p.Mass()<<")\n";
#endif
    AddJ(j);
  }
}

template <typename SType>
void CF<SType>::SetGauge(const ATOOLS::Vec4D &k)
{
}

template <typename SType>
void CF<SType>::AddPropagator()
{
  const CSpinorType hs;
  // add propagator for off-shell leg
  SComplex prop(M_I/(SType(this->m_p.Abs2())-m_cmass2));
  SComplex pp(hs.PPlus(this->m_p)), pm(hs.PMinus(this->m_p));
  SComplex pt(hs.PT(this->m_p)), ptc(hs.PTC(this->m_p));
#ifdef DEBUG__BG
  msg_Debugging()<<"propagator: "<<prop
		 <<" <- p^2 = "<<this->m_p.Abs2()<<", m = "<<m_cmass<<"\n";
  msg_Debugging()<<"pp = "<<pp<<", pm = "<<pm<<", pt = "<<pt<<"\n";
#endif
  for (typename CSpinorType_Vector::iterator 
	 jit(this->m_j.begin());jit!=this->m_j.end();++jit) {
    CSpinorType j(jit->R(),jit->B(),(*jit)(),jit->H(0),jit->H(1),
		  (jit->On()&1)<<1|(jit->On()&2)>>1);
    if (jit->B()>0) {
      j[0]=pm*(*jit)[2]-ptc*(*jit)[3];
      j[1]=-pt*(*jit)[2]+pp*(*jit)[3];
      j[2]=pp*(*jit)[0]+ptc*(*jit)[1];
      j[3]=pt*(*jit)[0]+pm*(*jit)[1];
    }
    else {
      j[0]=(*jit)[2]*pp+(*jit)[3]*pt;
      j[1]=(*jit)[2]*ptc+(*jit)[3]*pm;
      j[2]=(*jit)[0]*pm-(*jit)[1]*pt;
      j[3]=-(*jit)[0]*ptc+(*jit)[1]*pp;
    }
    if (jit->B()>0) *jit=(this->m_msv?j-*jit*m_cmass:j)*prop;
    else *jit=-(this->m_msv?j+*jit*m_cmass:j)*prop;
  }
}

template <typename SType> std::complex<SType> 
CF<SType>::SContract(const Current_Base &c,
		     const long unsigned int &hm,
		     const long unsigned int &hp) const
{
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"("<<hm<<","<<hp<<"): {\n";
  msg_Indent();
#endif
  if (c.Type()!='F') THROW(fatal_error,"Invalid current type.");
  Complex res(0.0);
  const CSpinorType_Vector &cj(c.Current_Base::J<CSpinorType>());
  for (typename CSpinorType_Vector::const_iterator 
	 jit1(this->m_j.begin());jit1!=this->m_j.end();++jit1)
    for (typename CSpinorType_Vector::const_iterator 
	   jit2(cj.begin());jit2!=cj.end();++jit2) 
      if ((*jit1)()==(*jit2)() &&
	  jit1->H(0)+jit2->H(0)==hm && jit1->H(1)+jit2->H(1)==hp) {
#ifdef DEBUG__BG
	msg_Debugging()<<"add "<<*jit1**jit2<<"\n";
	msg_Debugging()<<"    "<<*jit1<<"\n";
	msg_Debugging()<<"    "<<*jit2<<"\n";
#endif
	res+=*jit1**jit2;
      }
#ifdef DEBUG__BG
  msg_Debugging()<<"}\n";
#endif
  return res;
}

template <typename SType>
void CF<SType>::SContract(const Current_Base &c,
			  const LongUIntMap_Matrix &hmp,
			  SComplex_Vector &ress) const
{
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"(): {\n";
  msg_Indent();
#endif
  if (c.Type()!='F') THROW(fatal_error,"Invalid current type.");
  const CSpinorType_Vector &cj(c.Current_Base::J<CSpinorType>());
  for (typename CSpinorType_Vector::const_iterator 
	 jit1(this->m_j.begin());jit1!=this->m_j.end();++jit1)
    for (typename CSpinorType_Vector::const_iterator 
	   jit2(cj.begin());jit2!=cj.end();++jit2) 
      if ((*jit1)(0)==(*jit2)(1) && (*jit1)(1)==(*jit2)(0)) {
	LongUIntMap_Matrix::const_iterator 
	  hmit(hmp.find(jit1->H(0)+jit2->H(0)));
	if (hmit==hmp.end()) continue;
	LongUInt_Map::const_iterator 
	  hpit(hmit->second.find(jit1->H(1)+jit2->H(1)));
	if (hpit==hmit->second.end()) continue;
#ifdef DEBUG__BG
	msg_Debugging()<<"Add ("<<hmit->first<<","
		       <<hpit->first<<")"<<*jit1**jit2<<"\n";
#endif
	ress[hpit->second]+=*jit1**jit2;
      }
#ifdef DEBUG__BG
  msg_Debugging()<<"}\n";
#endif
}

template <typename SType>
char CF<SType>::Type() const
{
  return 'F';
}

template <typename SType>
std::string CF<SType>::CLabel() const
{
  return "fermion,label.side=right,label=$"+
    (this->m_out.empty()?this->m_fl.Bar():this->m_fl).TexName()+"$";
}

template <typename SType>
CV<SType>::CV(const Current_Key &key): 
  Current<CVec4Type>(key), m_cmass2(0.0), m_cmass(0.0)
{
  switch (key.m_fl.Kfcode()) {
  case kf_Wplus:
    m_cmass=sqrt(m_cmass2=key.p_model->Constant("m_{W,PS}^2"));
    break;
  case kf_Z:
    m_cmass=sqrt(m_cmass2=key.p_model->Constant("m_{Z,PS}^2"));
    break;
  default: break;
  }
}

DECLARE_TEMPLATE_GETTER(CF_Getter,"F",Current_Base,Current_Key);

template <typename SType,char STag> Current_Base *
CF_Getter<SType,STag>::operator()(const Current_Key &key) const
{
  if (key.m_fl.IsFermion()) return new CF<SType>(key);
  return NULL;
}

template <typename SType,char STag>
void CF_Getter<SType,STag>::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"fermion current "<<STag;
}

template <typename SType> CVec4<SType> 
CV<SType>::VT(const SpinorType &a,const SpinorType &b)
{
  CVec4Type e;
  e[0]=a.U1()*b.U1()+a.U2()*b.U2();
  e[3]=a.U1()*b.U1()-a.U2()*b.U2();
  e[1]=a.U1()*b.U2()+a.U2()*b.U1();
  e[2]=SComplex(0.0,1.0)*(a.U1()*b.U2()-a.U2()*b.U1());
  return e;
}

template <typename SType> CVec4<SType>
CV<SType>::EM(const Vec4D &p,const int cr,const int ca)
{
  SpinorType pp(1,p);
  CVec4Type e(VT(pp,m_km));
  e(0)=cr; e(1)=ca;
  e.SetH(1<<this->m_key,this->m_dir>0?1:0);
  static SType sqrttwo(sqrt(SType(2.0)));
  return e/(sqrttwo*std::conj(m_kp*pp));
}

template <typename SType> CVec4<SType>
CV<SType>::EP(const Vec4D &p,const int cr,const int ca)
{
  SpinorType pm(-1,p);
  CVec4Type e(VT(m_kp,pm));
  e(0)=cr; e(1)=ca;
  e.SetH(1<<this->m_key,this->m_dir>0?0:1);
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
  double a(p.Abs2()/(2.0*m_k*p));
  Vec4D b(p-a*m_k);
  SpinorType bm(-1,b), bp(1,b), am(-1,m_k), ap(1,m_k);
  CVec4Type e(VT(bp,bm)-SType(a)*VT(ap,am));
  e(0)=cr; e(1)=ca;
  e.SetH(1<<this->m_key,0);
  e.SetH(1<<this->m_key,1);
  return e/(SType(2.0)*m_cmass);
}

template <typename SType>
void CV<SType>::ConstructJ(const ATOOLS::Vec4D &p,const int ch,
			   const int cr,const int ca)
{
  this->m_p=p;
  if (this->m_fl.Mass()==0.0 && p[1]==0.0 && p[2]==0.0)
    this->m_p[0]=this->m_p[0]<0.0?
      -std::abs(this->m_p[3]):std::abs(this->m_p[3]);
  this->ResetJ();
  if (ch>=0) {
    if (this->m_msv && (ch==0 || ch==3)) {
      CVec4Type j(EML(this->m_p,cr,ca));
      AddJ(this->m_dir>0?j:j.Conj());
#ifdef DEBUG__BG
      msg_Debugging()<<METHOD<<"(): "<<(this->m_dir>0?'I':'O')
		     <<"0 "<<this->m_id<<" "<<this->m_j.back()
		     <<" "<<this->m_fl<<", m = "<<m_cmass<<"\n";
#endif
    }
    if (ch!=3) {
      CVec4Type j(this->m_msv?this->m_dir>0?
		  EMM(this->m_p,cr,ca):EMP(this->m_p,cr,ca):
		  this->m_dir>0?EM(this->m_p,cr,ca):EP(this->m_p,cr,ca));
      AddJ(this->m_dir>0?j:j.Conj());
#ifdef DEBUG__BG
      msg_Debugging()<<METHOD<<"(): "<<(this->m_dir>0?'I':'O')
		     <<"+ "<<this->m_id<<" "<<this->m_j.back()
		     <<this->m_fl<<", m = "<<m_cmass<<"\n";
#endif
    }
  }
  if (ch<=0) {
    CVec4Type j(this->m_msv?this->m_dir>0?
		EMP(this->m_p,cr,ca):EMM(this->m_p,cr,ca):
		this->m_dir>0?EP(this->m_p,cr,ca):EM(this->m_p,cr,ca));
    AddJ(this->m_dir>0?j:j.Conj());
#ifdef DEBUG__BG
    msg_Debugging()<<METHOD<<"(): "<<(this->m_dir>0?'I':'O')
		   <<"- "<<this->m_id<<" "<<this->m_j.back()
		   <<" "<<this->m_fl<<", m = "<<m_cmass<<"\n";
#endif
  }
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
  SComplex prop(-M_I/(SType(this->m_p.Abs2())-m_cmass2));
#ifdef DEBUG__BG
  msg_Debugging()<<"propagator: "<<prop<<"\n";
#endif
  if (!this->m_fl.Strong()) {
    if (!this->m_msv)
      for (typename CVec4Type_Vector::iterator 
	     jit(this->m_j.begin());jit!=this->m_j.end();++jit)
	*jit-=(*jit*Vec4Type(this->m_p))*CVec4Type(this->m_p)*prop;
    else
      for (typename CVec4Type_Vector::iterator 
	     jit(this->m_j.begin());jit!=this->m_j.end();++jit)
	*jit-=(*jit*Vec4Type(this->m_p))*CVec4Type(this->m_p)/m_cmass2;
  }
  for (typename CVec4Type_Vector::iterator 
	 jit(this->m_j.begin());jit!=this->m_j.end();++jit) *jit*=prop;
}

template <typename SType> std::complex<SType>
CV<SType>::SContract(const Current_Base &c,
		     const long unsigned int &hm,
		     const long unsigned int &hp) const
{
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"("<<hm<<","<<hp<<"): {\n";
  msg_Indent();
#endif
  if (c.Type()!='V') THROW(fatal_error,"Invalid current type.");
  SComplex res(0.0);
  const CVec4Type_Vector &cj(c.Current_Base::J<CVec4Type>());
  for (typename CVec4Type_Vector::const_iterator 
	 jit1(this->m_j.begin());jit1!=this->m_j.end();++jit1)
    for (typename CVec4Type_Vector::const_iterator 
	   jit2(cj.begin());jit2!=cj.end();++jit2) 
      if ((*jit1)(0)==(*jit2)(1) && (*jit1)(1)==(*jit2)(0) &&
	  jit1->H(0)+jit2->H(0)==hm && jit1->H(1)+jit2->H(1)==hp) {
#ifdef DEBUG__BG
	msg_Debugging()<<"add "<<*jit1**jit2<<"\n";
	msg_Debugging()<<"    "<<*jit1<<"\n";
	msg_Debugging()<<"    "<<*jit2<<"\n";
#endif
	res+=*jit1**jit2;
      }
#ifdef DEBUG__BG
  msg_Debugging()<<"}\n";
#endif
  return res;
}

template <typename SType>
void CV<SType>::SContract(const Current_Base &c,
			  const LongUIntMap_Matrix &hmp,
			  SComplex_Vector &ress) const
{
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"(): {\n";
  msg_Indent();
#endif
  if (c.Type()!='V') THROW(fatal_error,"Invalid current type.");
  const CVec4Type_Vector &cj(c.Current_Base::J<CVec4Type>());
  for (typename CVec4Type_Vector::const_iterator 
	 jit1(this->m_j.begin());jit1!=this->m_j.end();++jit1)
    for (typename CVec4Type_Vector::const_iterator 
	   jit2(cj.begin());jit2!=cj.end();++jit2) 
      if ((*jit1)(0)==(*jit2)(1) && (*jit1)(1)==(*jit2)(0)) {
	LongUIntMap_Matrix::const_iterator 
	  hmit(hmp.find(jit1->H(0)+jit2->H(0)));
	if (hmit==hmp.end()) continue;
	LongUInt_Map::const_iterator 
	  hpit(hmit->second.find(jit1->H(1)+jit2->H(1)));
	if (hpit==hmit->second.end()) continue;
#ifdef DEBUG__BG
	msg_Debugging()<<"Add ("<<hmit->first<<","
		       <<hpit->first<<")"<<*jit1**jit2<<"\n";
#endif
	ress[hpit->second]+=*jit1**jit2;
      }
#ifdef DEBUG__BG
  msg_Debugging()<<"}\n";
#endif
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

DECLARE_TEMPLATE_GETTER(CV_Getter,"V",Current_Base,Current_Key);

template <typename SType,char STag> Current_Base *
CV_Getter<SType,STag>::operator()(const Current_Key &key) const
{
  if (key.m_fl.IsVector()) return new CV<SType>(key);
  return NULL;
}

template <typename SType,char STag>
void CV_Getter<SType,STag>::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"vector current "<<STag;
}

template <typename SType>
CT<SType>::CT(const Current_Key &key): 
  Current<CAsT4Type>(key), m_prop(-M_I) 
{
}

template <typename SType>
void CT<SType>::ConstructJ(const ATOOLS::Vec4D &p,const int ch,
			   const int cr,const int ca)
{
}

template <typename SType>
void CT<SType>::SetGauge(const ATOOLS::Vec4D &k)
{
}

template <typename SType>
void CT<SType>::AddPropagator()
{
  // add propagator for off-shell leg
  for (typename CAsT4Type_Vector::iterator 
	 jit(this->m_j.begin());jit!=this->m_j.end();++jit) *jit*=m_prop;
}

template <typename SType> std::complex<SType>
CT<SType>::SContract(const Current_Base &c,
		     const long unsigned int &hm,
		     const long unsigned int &hp) const
{
  THROW(fatal_error,"Multiplication of tensor particles not allowed");
  return 0.0;
}

template <typename SType>
void CT<SType>::SContract(const Current_Base &c,
			  const LongUIntMap_Matrix &hmp,
			  SComplex_Vector &ress) const
{
  THROW(fatal_error,"Multiplication of tensor particles not allowed");
}

template <typename SType>
char CT<SType>::Type() const
{
  return 'T';
}

template <typename SType>
std::string CT<SType>::CLabel() const
{
  return "double,label.side=right,label=$"+
    (this->m_out.empty()?this->m_fl.Bar():this->m_fl).TexName()+"$";
}
#undef DECLARE_TEMPLATE_GETTER_BASE
#undef DECLARE_TEMPLATE_GETTER

#define DECLARE_TEMPLATE_GETTER_BASE(NAME,TAG,OBJECT,PARAMETER,SORT)	\
  									\
  template <typename TP,char TN>					\
  class NAME: public ATOOLS::Getter_Function<OBJECT,PARAMETER,SORT> {	\
  private:								\
    static NAME s_initializer;						\
  protected:								\
    void PrintInfo(std::ostream &str,const size_t width) const;		\
    Object_Type *							\
    operator()(const Parameter_Type &parameters) const;	                \
  public:								\
    NAME(const bool &d=true): ATOOLS::Getter_Function			\
      <OBJECT,PARAMETER,SORT>(std::string(1,TN)+TAG) { SetDisplay(d); }	\
  };									

#define DECLARE_TEMPLATE_GETTER(NAME,TAG,OBJECT,PARAMETER)		\
  DECLARE_TEMPLATE_GETTER_BASE(NAME,TAG,OBJECT,PARAMETER,		\
			       std::less<std::string>)			\
    template <typename TP,char TN>					\
  NAME<TP,TN> NAME<TP,TN>::s_initializer				 

DECLARE_TEMPLATE_GETTER(CT_Getter,"T",Current_Base,Current_Key);

template <typename SType,char STag> Current_Base *
CT_Getter<SType,STag>::operator()(const Current_Key &key) const
{
  if (key.m_fl.IsTensor()) return new CT<SType>(key);
  return NULL;
}

template <typename SType,char STag>
void CT_Getter<SType,STag>::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"tensor current "<<STag;
}

template class CS_Getter<double,'D'>;
template class CF_Getter<double,'D'>;
template class CV_Getter<double,'D'>;
template class CT_Getter<double,'D'>;

template class CS_Getter<long double,'Q'>;
template class CF_Getter<long double,'Q'>;
template class CV_Getter<long double,'Q'>;
template class CT_Getter<long double,'Q'>;
