#include "METOOLS/Explicit/Current.H"
#include "METOOLS/Currents/C_RaritaSchwinger.H"
#include "METOOLS/Currents/C_Vector.H"
#include "METOOLS/Currents/C_Spinor.H"
#include "ATOOLS/Phys/Spinor.H"
#include "METOOLS/Main/SpinFuncs.H"
//#define ZERO SComplex(0.0,0.0)

namespace METOOLS {

  template <typename SType>
  class CRS: public Current,
	    public Current_Contractor<SType> {
  public:

    typedef std::complex<SType>   SComplex;
    typedef std::vector<SComplex> SComplex_Vector;

    typedef ATOOLS::Spinor<SType> SpinorType;
    typedef ATOOLS::Vec4<SType>   Vec4Type;

    typedef CSpinor<SType> CSpinorType;
    typedef std::vector<CSpinorType*> CSpinorType_Vector;

    typedef CVec4<SType> CVec4Type;
    typedef std::vector<CVec4Type*> CVec4Type_Vector;

    typedef CRaritaSchwinger<SType> CRaScType;
    typedef std::vector<CRaScType*> CRaScType_Vector;

  protected:

    SComplex m_cmass2, m_cmass;

    //std::string CLabel() const;

  private:

    ATOOLS::Vec4D m_k;
    SpinorType    m_kp, m_km;

    // construction of polarization vector
    CVec4Type VT(const SpinorType &a, const SpinorType &b);

    CVec4Type EM(const ATOOLS::Vec4D &p, int cr, int ca);
    CVec4Type EP(const ATOOLS::Vec4D &p, int cr, int ca);

    CVec4Type EMM(const ATOOLS::Vec4D &p,const int cr,const int ca);
    CVec4Type EMP(const ATOOLS::Vec4D &p,const int cr,const int ca);
    CVec4Type EML(const ATOOLS::Vec4D &p,const int cr,const int ca);

    // different polarization states
    CRaScType RSMM(const ATOOLS::Vec4D &p, int r, int s, int b, int cr=0, int ca=0, int hh=0, int ms=0);
    CRaScType RSPP(const ATOOLS::Vec4D &p, int r, int s, int b, int cr=0, int ca=0, int hh=0, int ms=0);
    CRaScType RSP(const ATOOLS::Vec4D &p,const int cr,const int ca);
    CRaScType RSM(const ATOOLS::Vec4D &p,const int cr,const int ca);

    CRaScType SpinorVectorProduct(const CSpinorType spinor, const CVec4Type polvector, SType m2,
                                  int cr=0, int ca=0, int s=0);

  public:
    CRS(const Current_Key &key);
    void ConstructJ(const ATOOLS::Vec4D &p,const int ch,
		    const int cr,const int ca,const int mode);
    void SetGauge(const ATOOLS::Vec4D &k);

    void AddPropagator();

    void SContract
    (const Current &c,const Int_Vector &pols,
     SComplex_Vector &ress,const size_t &offset) const;

    std::string Format(const CObject *c) const;

    char Type() const;

    // Tests
    bool Test_WF_Properties(const ATOOLS::Vec4D &p);

    // TODO: Necessary or only for test purpose?
    inline CRaScType Get_RSPP(const ATOOLS::Vec4D &p, int r, int s, int b, int cr=0, int ca=0, int hh=0, int ms=0)
    { return RSPP(p, r, s, b);};
    inline CRaScType Get_RSMM(const ATOOLS::Vec4D &p, int r, int s, int b, int cr=0, int ca=0, int hh=0, int ms=0) const
    { return RSMM(p, r, s, b);};

  };

// end of class CV

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

// TODO: Bislang identisch zu Fermion- und Bosonkonstruktor: ist das richtig?
template <typename SType>
CRS<SType>::CRS(const Current_Key &key):
  Current(key), m_cmass2(0.0), m_cmass(0.0)
{
  m_cmass=sqrt(m_cmass2=SComplex(sqr(this->m_mass),-this->m_mass*this->m_width));
  if (key.m_n==1 && key.p_model->ScalarNumber("WidthScheme")!=1)
    m_cmass=sqrt(m_cmass2=Complex(sqr(this->m_mass),0.0));
}

// FUNKTIONS TO CALCULATE THE RARITA-SCHWINGER-WAVEFUNKTION
template <typename SType> CVec4<SType>
CRS<SType>::VT(const SpinorType &a,const SpinorType &b)
{
  CVec4Type e;
  e[0]=a.U1()*b.U1()+a.U2()*b.U2();
  e[SpinorType::R3()]=a.U1()*b.U1()-a.U2()*b.U2();
  e[SpinorType::R1()]=a.U1()*b.U2()+a.U2()*b.U1();
  e[SpinorType::R2()]=SComplex(0.0,1.0)*(a.U1()*b.U2()-a.U2()*b.U1());
  return e;
}
// TODO: PLUS AND MINUS EXCHANGED
template <typename SType> CVec4<SType>
CRS<SType>::EM(const Vec4D &p,const int cr,const int ca)
{
  SpinorType pp(1,p);
  CVec4Type e(VT(pp,m_km));
  e(0)=cr; e(1)=ca;
  e.SetH(1);
  static SType sqrttwo(sqrt(SType(2.0)));
  return e/(sqrttwo*std::conj(m_kp*pp));
}

template <typename SType> CVec4<SType>
CRS<SType>::EP(const Vec4D &p,const int cr,const int ca)
{
  SpinorType pm(-1,p);
  CVec4Type e(VT(m_kp,pm));
  e(0)=cr; e(1)=ca;
  e.SetH(0);
  static SType sqrttwo(sqrt(SType(2.0)));
  return e/(sqrttwo*std::conj(m_km*pm));
}
// TODO: p.Abs2() wirklich immer null, wenn wir masselose Teilchen wollen? Propagator sind extra, d.h. hier werden nur
//       onshell Teilchen erzeugt?
template <typename SType> CVec4<SType>
CRS<SType>::EMM(const Vec4D &p,const int cr,const int ca)
{
  return EM(p-p.Abs2()/(2.0*m_k*p)*m_k,cr,ca);
}

template <typename SType> CVec4<SType>
CRS<SType>::EMP(const Vec4D &p,const int cr,const int ca)
{
  return EP(p-p.Abs2()/(2.0*m_k*p)*m_k,cr,ca);
}

template <typename SType> CVec4<SType>
CRS<SType>::EML(const Vec4D &p,const int cr,const int ca)
{
  double p2(p.Abs2()), a(p2/(2.0*m_k*p));
  Vec4D b(p-a*m_k);
  SpinorType bm(-1,b), bp(1,b), am(-1,m_k), ap(1,m_k);
  CVec4Type e(VT(bp,bm)-SType(a)*VT(ap,am));
  e(0)=cr; e(1)=ca;
  e.SetH(2);
  return e/sqrt(SComplex(4.0*p2));
}

// TODO: PLUS AND MINUS EXCHANGED
// TODO: GENAUE DEFINITION BESCHREIBEN!!!
// TODO: Nach HELAS paper brauchen wir s, cr, ca, hh, ms nicht...
// TODO: Im HELAS-Paper hat die Wellenfunktion 18 Komponenten, wobei die letzten bei den Komponenten den Viererimpuls
//       entlang des Fermionzahlflusses enthalten -> brauchen wir das auch?
// - Füllungsreihenfolge wie in MadGraph /HELAS paper arXiv: 1010.4255
// - Form wie in S.F.Novaes & D.Spehler Nuclear Physics B 371 (1992), 618-636 Eq.(13) mit Phase theta = 0
template<typename SType> CRaritaSchwinger<SType>
CRS<SType>::RSPP(const ATOOLS::Vec4D &p, const int r, const int s, const int b, const int cr, const int ca,
                 const int hh, const int ms) {
  return
  b>0?METOOLS::CRS<SType>::SpinorVectorProduct(CSpinor(r, abs(b), 1, p, cr, ca, hh, s, p.Abs2(), ms),
                                                  EMM(p, cr, ca), p.Abs2(), cr, ca, s):
                                                  METOOLS::CRS<SType>::SpinorVectorProduct(CSpinor(r, abs(b), 1, p, cr,
                                                                                                       ca, hh, s,
                                                                                                       p.Abs2(), ms),
                                                                                               EMM(p, cr, ca), p.Abs2(),
                                                                                               cr, ca, s).Bar();
}

template<typename SType> CRaritaSchwinger<SType>
CRS<SType>::RSMM(const ATOOLS::Vec4D &p, const int r, const int s, const int b, const int cr, const int ca,
                 const int hh, const int ms) {
  return b>0?METOOLS::CRS<SType>::SpinorVectorProduct(CSpinor(r, abs(b), -1, p, cr, ca, hh, s, p.Abs2(), ms),
                                                  EMP(p, cr, ca), p.Abs2(), cr, ca, s):
                                                  METOOLS::CRS<SType>::SpinorVectorProduct(CSpinor(r, abs(b), -1, p, cr, ca,
                                                                                                   hh, s, p.Abs2(), ms),
                                                                                           EMP(p, cr, ca), p.Abs2(), cr,
                                                                                           ca, s).Bar();
}

template<typename SType>
CRaritaSchwinger<SType> CRS<SType>::SpinorVectorProduct(const CRS::CSpinorType spinor, const CRS::CVec4Type polvector,
                                                      const SType m2, const int cr, const int ca, const int s) {
  // set properties of new Rarita-Schwinger particle
  int vector_h = polvector.H();
  // convert in more approriate numbering, h=0 : longitudial, h=2: right, h=-2: left
  // TODO: ADJUST IF +- AND -SWITCH IS SOLVED!!!
  if (vector_h==2) vector_h=0;  // long
  else if (vector_h==1) vector_h=2; // right
  else if (vector_h==0) vector_h=-2; // left

  METOOLS::CRS<SType>::CRaScType RaSc(spinor.R(), spinor.B(), cr, ca, vector_h+spinor.H(), s);

  // TODO: Evtl. auch aus Spinor und Polarisationsvektor ablesbar -> Bedeutung herausfinden!!!!
  RaSc(0) = cr;
  RaSc(1) = ca;

  // TODO: Ist diese Version der Komponentenbefüllung richtig oder z.B. diese:
//  for (size_t i(0); i<8; ++i){
//    std::cout << "Fill components" << i << i/2 << std::endl;
//    RaSc[i] = spinor[i % 2] * polvector[i / 2];
//    RaSc[i+8] = spinor[i % 2 + 2] * polvector[i / 2];
// oder:     for (size_t i(0); i<16; ++i) {
//      RaSc[i] = spinor[i / 4] * polvector[i % 4];
//    }

  // Fill Rarita-Schwinger wave function; "ordering" of components according to HELAS subroutines for spin-3/2 particles
  // in 1010.4255 (four vector with Dirac spinors as components, i.e. first four components have Lorentz index 0 and
  // spinor indexes from 0 to 3, the second four components Lorentz index 1 ...), but using Dirac-spinors and
  // polarization vectors as already implemented in SHERPA
  if (IsZero(m2)){
    for (size_t i(0); i<16; ++i)
      RaSc[i] = spinor[i % 4] * polvector[i / 4];
  }
  bool on = RaSc.SetOn();
  return RaSc;
}


/*template <typename SType> CRaritaSchwinger<SType>
CRS<SType>::RSM(const Vec4D &p, const int r, const int cr,const int ca)
{p-p.Abs2()/(2.0*m_k*p)*m_k,cr,ca
  return SpinorVectorProduct(METOOLS::CSpinor(-1,1,const int &h,
  const ATOOLS::Vec4<Scalar> &p,
  const int cr=0,const int ca=0,
  const size_t &hh=0,const size_t &s=0,
  const Scalar &m2=-1.0,const int ms=1));
}

template <typename SType> CRaritaSchwinger<SType>
CRS<SType>::RSP(const Vec4D &p,const int cr,const int ca)
{
  return EP(p-p.Abs2()/(2.0*m_k*p)*m_k,cr,ca);
}

template <typename SType> CRaritaSchwinger<SType>
CRS<SType>::RSL(const Vec4D &p,const int cr,const int ca)
{
  double p2(p.Abs2()), a(p2/(2.0*m_k*p));
  Vec4D b(p-a*m_k);
  SpinorType bm(-1,b), bp(1,b), am(-1,m_k), ap(1,m_k);
  CVec4Type e(VT(bp,bm)-SType(a)*VT(ap,am));
  e(0)=cr; e(1)=ca;
  e.SetH(2);
  return e/sqrt(SComplex(4.0*p2));
}*/

template <typename SType>
void CRS<SType>::ConstructJ(const ATOOLS::Vec4D &p,const int ch,
			   const int cr,const int ca,const int mode)
{
  this->m_p=p;

  //TODO: Brauchen wir das alles?
  // TODO: Wozu ist das, sollte das nicht grundsätzlich gelten?
  if (this->m_fl.Mass()==0.0 && p[1]==0.0 && p[2]==0.0)
    this->m_p[0]=this->m_p[0]<0.0?
      -std::abs(this->m_p[3]):std::abs(this->m_p[3]);
  // TODO: Sollten wir Majorana-Spin-3/2 Teilchen unterstützen?
  bool anti(this->m_fl.IsAnti());
  if (this->m_fl.Majorana()) anti=(mode&1)?this->m_dir<0:this->m_dir>0;
  this->ResetJ();
  //TODO: Was bedeutet ch? Was für Werte kann ch für RaSc annehmen?
  // TODO: Wie wird dann am Ende h gesetzt bei Fermionen?
  // TODO: Brauchen wir noch Vec-Parameter?
  // TODO: Stimmt RSPP /RSMM? Haben die dann die richtigen helizitäten/Bars?
  if (ch>=0) {
    CRaScType j(anti^(this->m_dir>0)? RSPP(p, -1, 0, -this->m_dir, cr, ca): RSMM(p, 1, 0, this->m_dir, cr, ca));
    j.SetH(anti^(this->m_dir>0)?1:0);
#ifdef DEBUG__BG
    msg_Debugging()<<METHOD<<"(): "<<(this->m_dir>0?'I':'O')<<"+ "<<this->m_id
		   <<" "<<j<<" "<<(this->m_dir>0?this->m_fl.Bar():this->m_fl)
		   <<", m = "<<m_cmass<<" ("<<p.Mass()<<")\n";
#endif
    CRaScType *c(CRaScType::New(j));
    AddJ(c);
    if (p_sub) static_cast<Dipole_Color*>
      (p_sub->In().front()->Color().front())->AddJJK(c);
  }
  if (ch<=0) {
    CRaScType j(anti^(this->m_dir>0) ? RSPP(p, -1, 0, -this->m_dir, cr, ca) : RSMM(p, 1, 0, this->m_dir, cr, ca));
    j.SetH(anti^(this->m_dir>0)?0:1);
#ifdef DEBUG__BG
    msg_Debugging()<<METHOD<<"(): "<<(this->m_dir>0?'I':'O')<<"- "<<this->m_id
		   <<" "<<j<<" "<<(this->m_dir>0?this->m_fl.Bar():this->m_fl)
		   <<", m = "<<m_cmass<<" ("<<p.Mass()<<")\n";
#endif
    CRaScType *c(CRaScType::New(j));
    AddJ(c);
    if (p_sub) static_cast<Dipole_Color*>
      (p_sub->In().front()->Color().front())->AddJJK(c);
  }
#ifdef DEBUG__BG
  if (p_sub) Print();
#endif
}

template <typename SType>
void CRS<SType>::SetGauge(const ATOOLS::Vec4D &k)
{
  m_k=k;
  m_kp=SpinorType(1,m_k);
  m_km=SpinorType(-1,m_k);
}

template <typename SType>
void CRS<SType>::AddPropagator()
{}
/*
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
}*/

// For contracting the matrix element with its complex conjugate
template <typename SType> void CRS<SType>::SContract
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
  if (c.Type()!='R') THROW(fatal_error,"Invalid current type.");
  size_t i(0);
  for (typename CObject_Matrix::const_iterator 
	 ajit1(m_j.begin());ajit1!=m_j.end();++ajit1) {	
    const CRaScType_Vector *j(ajit1->Get<CRaScType>());
    for (typename CObject_Matrix::const_iterator 
	   ajit2(c.J().begin());ajit2!=c.J().end();++ajit2,++i) {
      // if (!pols[i]) continue;
      const CRaScType_Vector *cj(ajit2->Get<CRaScType>());
      for (typename CRaScType_Vector::const_iterator
	     jit2(cj->begin());jit2!=cj->end();++jit2) 
	for (typename CRaScType_Vector::const_iterator
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
/*
template <typename SType>
std::string CRS<SType>::CLabel() const
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
*/

template <typename SType>
std::string CRS<SType>::Format(const CObject *c) const
{
  return ToString(*(CRaScType *)c,6);
}

template <typename SType>
char CRS<SType>::Type() const
{
  return 'R';
}

// TODO: Wieso failt Normierungstest bei den einlaufenden Teilchen des Prozesses? (also mit deren Impulsen)
template<typename SType>
bool CRS<SType>::Test_WF_Properties(const ATOOLS::Vec4D &p) {
  CRaritaSchwinger<SType> rspp = RSPP(p, 1, 1, 1, 0, 0, 1);
  CRaritaSchwinger<SType> rsmm = RSMM(p, 1, 1, 1, 0, 0, -1);
  METOOLS::Gamma<SType> gammavec = Gamma<SType>();
  if (ATOOLS::IsZero(p.Abs2())){
    // normalization
    std::cout<<METHOD<<": Testing normalization of Rarita-Schwinger wave function..."<<std::endl;
    TCMatrix<SType> result1 = rspp.Contract4Index(rspp.Bar()) + rsmm.Contract4Index(rsmm.Bar()) + gammavec * p;

    for (size_t i(0); i<4; ++i){
      for (size_t j(0); j<4; ++j){
        if (std::abs(result1[i][j].real())> rspp.Accu()|| std::abs(result1[i][j].imag())>rspp.Accu()) {
          msg_Out()<<"Component " << i << j << " of resulting 4x4 matrix is " << result1[i][j] << ", instead of zero!"
          << std::endl;
          return false;
        }
      }
    }

    // equality between U++ and V-- / U-- and V++ for bar and non-bar

  }
  // normalization

  // completness relation
  return true;
}

DECLARE_GETTER(CRS<double>,"DR",Current,Current_Key);

Current *ATOOLS::Getter<Current,Current_Key,CRS<double> >::
operator()(const Current_Key &key) const
{
  if (key.m_fl.IsRaritaSchwinger()) return new CRS<double>(key);
  return NULL;
}

void ATOOLS::Getter<Current,Current_Key,CRS<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"rarita-schwinger current (double)";
}