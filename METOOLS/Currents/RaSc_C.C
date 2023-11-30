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

    CVec4Type EM(const ATOOLS::Vec4D &p, const int cr, const int ca, const ATOOLS::Vec4D m_k_mod);
    CVec4Type EP(const ATOOLS::Vec4D &p, int cr, int ca, const ATOOLS::Vec4D m_k_mod);

    CVec4Type EMM(const ATOOLS::Vec4D &p,const int cr,const int ca, const bool hel= true);
    CVec4Type EMP(const ATOOLS::Vec4D &p,const int cr,const int ca, const bool hel= true);
    CVec4Type EML(const ATOOLS::Vec4D &p,const int cr,const int ca, const bool hel= true);

    // different polarization states
    CRaScType RSMM(const ATOOLS::Vec4D &p, int r, int s, int b, int cr=0, int ca=0, int hh=0, int ms=0);
    CRaScType RSPP(const ATOOLS::Vec4D &p, int r, int s, int b, int cr=0, int ca=0, int hh=0, int ms=0);
    CRaScType RSP(const ATOOLS::Vec4D &p, int r, int s, int b, int cr=0, int ca=0, int hh=0, int ms=0);
    CRaScType RSM(const ATOOLS::Vec4D &p,int r, int s, int b, int cr=0, int ca=0, int hh=0, int ms=0);

    CRaScType SpinorVectorProduct(const CSpinorType spinor, const CVec4Type polvector, int spinor_h, int cr=0, int ca=0,
                                  int s=0);

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
    bool Test_WF_Properties(const ATOOLS::Vec4D &p, const bool &anti);

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
CRS<SType>::EM(const Vec4D &p, const int cr,const int ca, ATOOLS::Vec4D m_k_mod)
{
  SpinorType m_kp_mod=SpinorType(1,m_k_mod);
  SpinorType m_km_mod=SpinorType(-1,m_k_mod);
  SpinorType pp(1,p);
  CVec4Type e(VT(pp,m_km_mod));
  e(0)=cr; e(1)=ca;
  e.SetH(1);
  static SType sqrttwo(sqrt(SType(2.0)));
  return e/(sqrttwo*std::conj(m_kp_mod*pp));
}

template <typename SType> CVec4<SType>
CRS<SType>::EP(const Vec4D &p,const int cr,const int ca, ATOOLS::Vec4D m_k_mod)
{
  SpinorType m_kp_mod=SpinorType(1,m_k_mod);
  SpinorType m_km_mod=SpinorType(-1,m_k_mod);
  SpinorType pm(-1,p);
  CVec4Type e(VT(m_kp_mod,pm));
  e(0)=cr; e(1)=ca;
  e.SetH(0);
  static SType sqrttwo(sqrt(SType(2.0)));
  return e/(sqrttwo*std::conj(m_km_mod*pm));
}

template <typename SType> CVec4<SType>
CRS<SType>::EMM(const Vec4D &p,const int cr,const int ca, const bool hel)
{
  Vec4D m_k_mod;
  if (hel) m_k_mod=ATOOLS::Vec4D(1., -ATOOLS::Vec3D(p) / p.PSpat());
  else m_k_mod = m_k;
  return EM(p-p.Abs2()/(2.0*m_k_mod*p)*m_k_mod,cr,ca,m_k_mod);
}

template <typename SType> CVec4<SType>
CRS<SType>::EMP(const Vec4D &p,const int cr,const int ca, const bool hel)
{
  Vec4D m_k_mod;
  if (hel) m_k_mod=ATOOLS::Vec4D(1., -ATOOLS::Vec3D(p) / p.PSpat());
  else m_k_mod = m_k;
  return EP(p-p.Abs2()/(2.0*m_k_mod*p)*m_k_mod,cr,ca,m_k_mod);
}

template <typename SType> CVec4<SType>
CRS<SType>::EML(const Vec4D &p,const int cr,const int ca, const bool hel)
{
  Vec4D m_k_mod;
  if (hel) m_k_mod=ATOOLS::Vec4D(1., -ATOOLS::Vec3D(p) / p.PSpat());
  else m_k_mod = m_k;
  double p2(p.Abs2()), a(p2/(2.0*m_k_mod*p));
  Vec4D b(p-a*m_k_mod);
  SpinorType bm(-1,b), bp(1,b), am(-1,m_k_mod), ap(1,m_k_mod);
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
  CRaritaSchwinger<SType> wf(b>0?METOOLS::CRS<SType>::SpinorVectorProduct(CSpinor(r, abs(b), 1, p, cr, ca, hh, s,
                                                                                  this->m_msv?p.Abs2():0, ms),
                                                                          this->m_msv? EMM(p, cr, ca) : EM(p, cr, ca, m_k),
                                                                          1, cr, ca, s):
                                 METOOLS::CRS<SType>::SpinorVectorProduct(CSpinor(r, abs(b), 1, p, cr, ca, hh, s,
                                                                                  this->m_msv?p.Abs2():0, ms),
                                                                          this->m_msv? EMM(p, cr, ca) : EM(p, cr, ca, m_k),
                                                                          1, cr, ca, s).Bar());
  CRaritaSchwinger<SType> wf1 = METOOLS::CRS<SType>::SpinorVectorProduct(CSpinor(1, abs(b), 1, Vec4D(sqrt(p.Abs2()+1.0), 0, 0, 1.0), cr, ca, hh, s,
                                                                                 this->m_msv?p.Abs2():0, ms),
                                                                         this->m_msv? EMM(Vec4D(sqrt(p.Abs2()+1.0), 0, 0, 1.0), cr, ca) : EM(Vec4D(sqrt(p.Abs2()+1.0), 0, 0, 1.0), cr, ca, m_k),
  1, cr, ca, s);
  return wf;
}

template<typename SType> CRaritaSchwinger<SType>
CRS<SType>::RSP(const ATOOLS::Vec4D &p, int r, int s, int b, int cr, int ca, int hh, int ms){
  //SComplex exp_phi(SComplex(0, 0));
  /*if (p[1]>0) exp_phi = std::exp(SComplex(0, 1)*std::atan(p[2]/p[1]));
  else if (p[1]==0) exp_phi = std::exp(SComplex(0, 1)*(M_PI/2)*((p[2]>0)?SComplex(1,0):((p[2]==0)?0.0:SComplex(-1,0))));
  else if (p[1]<0 && p[2]>=0) exp_phi = std::exp(SComplex(0, 1)*(std::atan(p[2]/p[1])+M_PI));
  else if (p[1]<0 && p[2]>0) exp_phi = std::exp(SComplex(0, 1)*(std::atan(p[2]/p[1])-M_PI));*/
  if (!this->m_msv) THROW(fatal_error, "There is no massless Rarita-Schwinger-particle with Sz=1/2!")
  CRaritaSchwinger<SType> wf_p0(METOOLS::CRS<SType>::SpinorVectorProduct(CSpinor(r, abs(b), 1, p, cr, ca, hh, s, p.Abs2(),
                                                                              ms), EML(p, cr, ca), 1, cr, ca, s));
  CRaritaSchwinger<SType> wf_mp(METOOLS::CRS<SType>::SpinorVectorProduct(CSpinor(r, abs(b), -1, p, cr, ca, hh, s, p.Abs2(),
                                                                                 ms), EMM(p, cr, ca), 1, cr, ca, s));
  CRaritaSchwinger<SType> wf = b>0?(sqrt(2.0/3.0) * wf_p0 - sqrt(1.0/3.0) * wf_mp) :
                                   (sqrt(2.0/3.0) * wf_p0 - sqrt(1.0/3.0) * wf_mp).Bar();
  /*CRaritaSchwinger<SType> wf_p01(METOOLS::CRS<SType>::SpinorVectorProduct(CSpinor(1, abs(b), 1, Vec4D(sqrt(p.Abs2()+1.0), 0, 0, 1.0), cr, ca, hh, s, p.Abs2(),
                                                                                 ms), EML(Vec4D(sqrt(p.Abs2()+1.0), 0, 0, 1.0) , cr, ca), 1, cr, ca, s));
  CRaritaSchwinger<SType> wf_mp1(METOOLS::CRS<SType>::SpinorVectorProduct(CSpinor(1, abs(b), -1, Vec4D(sqrt(p.Abs2()+1.0), 0, 0, 1.0), cr, ca, hh, s, p.Abs2(),
                                                                                 ms), EMM(Vec4D(sqrt(p.Abs2()+1.0), 0, 0, 1.0), cr, ca), 1, cr, ca, s));
  SComplex exp_phi1(SComplex(0, 0));
  Vec4D p1 = Vec4D(sqrt(2), 0, 0, 1);
  if (p1[1]>0) exp_phi1 = std::exp(SComplex(0, 1)*std::atan(p1[2]/p1[1]));
  else if (p1[1]==0) exp_phi1 = std::exp(SComplex(0, 1)*(M_PI/2)*((p1[2]>0)?SComplex(1,0):((p1[2]==0)?0.0:SComplex(-1,0))));
  else if (p1[1]<0 && p1[2]>=0) exp_phi1 = std::exp(SComplex(0, 1)*(std::atan(p1[2]/p1[1])+M_PI));
  else if (p1[1]<0 && p1[2]>0) exp_phi1 = std::exp(SComplex(0, 1)*(std::atan(p1[2]/p1[1])-M_PI));
  std::cout << exp_phi1 << std::endl;
  CRaritaSchwinger<SType> wf1 = sqrt(2.0/3.0) * wf_p01 - sqrt(1.0/3.0) * wf_mp1;
  CRaritaSchwinger<SType> wf2 = sqrt(2.0/3.0) * wf_p01 + sqrt(1.0/3.0) * wf_mp1;
  CRaritaSchwinger<SType> wf3 = sqrt(2.0/3.0) * wf_p01 + exp_phi1*sqrt(1.0/3.0) * wf_mp1;
  std::cout << "Test-wf" << wf1 << std::endl;
  std::cout << "Test-wf" << wf2 << std::endl;
  std::cout << "Test-wf" << wf3 << std::endl;
  std::cout << "+" << wf << std::endl;*/
  return wf;
}

template<typename SType> CRaritaSchwinger<SType>
CRS<SType>::RSM(const ATOOLS::Vec4D &p, int r, int s, int b, int cr, int ca, int hh, int ms){
  if (!this->m_msv) THROW(fatal_error, "There is no massless Rarita-Schwinger-particle with Sz=1/2!")
  /*SComplex exp_phi(SComplex(0, 0));
  if (p[1]>0) exp_phi = std::exp(SComplex(0, 1)*std::atan(p[2]/p[1]));
  else if (p[1]==0) exp_phi = std::exp(SComplex(0, 1)*(M_PI/2)*((p[2]>0)?SComplex(1,0):((p[2]==0)?0.0:SComplex(-1,0))));
  else if (p[1]<0 && p[2]>=0) exp_phi = std::exp(SComplex(0, 1)*(std::atan(p[2]/p[1])+M_PI));
  else if (p[1]<0 && p[2]>0) exp_phi = std::exp(SComplex(0, 1)*(std::atan(p[2]/p[1])-M_PI));*/
  CRaritaSchwinger<SType> wf_pm(METOOLS::CRS<SType>::SpinorVectorProduct(CSpinor(r, abs(b), 1, p, cr, ca, hh, s, p.Abs2(),
                                                                                 ms), EMP(p, cr, ca), 1, cr, ca, s));
  CRaritaSchwinger<SType> wf_m0(METOOLS::CRS<SType>::SpinorVectorProduct(CSpinor(r, abs(b), -1, p, cr, ca, hh, s, p.Abs2(),
                                                                                 ms), EML(p, cr, ca), 1, cr, ca, s));
  CRaritaSchwinger<SType> wf = b>0?(sqrt(2.0/3.0) * wf_m0 + sqrt(1.0/3.0) * wf_pm) :
                               (sqrt(2.0/3.0) * wf_m0 + sqrt(1.0/3.0) * wf_pm).Bar();
  return wf;
}

template<typename SType> CRaritaSchwinger<SType>
CRS<SType>::RSMM(const ATOOLS::Vec4D &p, const int r, const int s, const int b, const int cr, const int ca,
                 const int hh, const int ms) {
  /*SComplex exp_phi(SComplex(0, 0));
  if (p[1]>0) exp_phi = std::exp(SComplex(0, 1)*std::atan(p[2]/p[1]));
  else if (p[1]==0) exp_phi = std::exp(SComplex(0, 1)*(M_PI/2)*((p[2]>0)?SComplex(1,0):((p[2]==0)?0.0:SComplex(-1,0))));
  else if (p[1]<0 && p[2]>=0) exp_phi = std::exp(SComplex(0, 1)*(std::atan(p[2]/p[1])+M_PI));
  else if (p[1]<0 && p[2]>0) exp_phi = std::exp(SComplex(0, 1)*(std::atan(p[2]/p[1])-M_PI));*/
  CRaritaSchwinger<SType> wf(b>0?METOOLS::CRS<SType>::SpinorVectorProduct(CSpinor(r, abs(b), -1, p, cr, ca, hh, s,
                                                                                  this->m_msv?p.Abs2():0, ms),
                                                                          this->m_msv? EMP(p, cr, ca) : EP(p, cr, ca, m_k), -1, cr, ca, s):
                             METOOLS::CRS<SType>::SpinorVectorProduct(CSpinor(r, abs(b), -1, p, cr, ca, hh, s, this->m_msv?p.Abs2():0,
                                                                              ms),
                                                                      this->m_msv? EMP(p, cr, ca) : EP(p, cr, ca, m_k), -1, cr,
                                                                      ca, s).Bar());
  return wf;
}

template<typename SType>
CRaritaSchwinger<SType> CRS<SType>::SpinorVectorProduct(const CRS::CSpinorType spinor, const CRS::CVec4Type polvector,
                                                        int spinor_h, const int cr, const int ca, const int s) {
  // set properties of new Rarita-Schwinger particle
  int vector_h = polvector.H();
  DEBUG_VAR(spinor);
  DEBUG_VAR(polvector);
  // convert in more approriate numbering, h=0 : longitudial, h=2: right, h=-2: left
  // TODO: ADJUST IF +- AND -SWITCH IS SOLVED!!!
  /*if (vector_h==2) vector_h=0;  // long
  else if (vector_h==1) vector_h=2; // right
  else if (vector_h==0) vector_h=-2; // left*/
  METOOLS::CRS<SType>::CRaScType RaSc(spinor.R(), spinor.B(), cr, ca, 0, s);
  RaSc(0) = cr;
  RaSc(1) = ca;

  // Fill Rarita-Schwinger wave function; "ordering" of components according to HELAS subroutines for spin-3/2 particles
  // in 1010.4255 (four vector with Dirac spinors as components, i.e. first four components have Lorentz index 0 and
  // spinor indexes from 0 to 3, the second four components Lorentz index 1 ...), but using Dirac-spinors and
  // polarization vectors as already implemented in SHERPA
  for (size_t i(0); i<16; ++i) {
    RaSc[i] = spinor[i % 4] * (spinor.R()>0?polvector[i / 4]:conj(polvector[i / 4]));
  }
  bool on = RaSc.SetOn();
  return RaSc;
}

template <typename SType>
void CRS<SType>::ConstructJ(const ATOOLS::Vec4D &p,const int ch,
			   const int cr,const int ca,const int mode)
{
  this->m_p=p;
  // TODO: Tests in RSPP, RSMM, und hier in DEBUG-Options sinnvoll einbauen
  if (this->m_fl.Mass()==0.0 && p[1]==0.0 && p[2]==0.0)
    this->m_p[0]=this->m_p[0]<0.0?
      -std::abs(this->m_p[3]):std::abs(this->m_p[3]);
  bool anti(this->m_fl.IsAnti());
  if (this->m_fl.Majorana()) anti=(mode&1)?this->m_dir<0:this->m_dir>0;
  this->ResetJ();
  int r(anti?-1:1);
  int b((anti^(this->m_dir>0))?1:-1);
  //TODO: Was bedeutet ch? Was für Werte kann ch für RaSc annehmen?
  // TODO: !!!Was ist nun der richtige Wert für SetH() 0 und 1 im masselosen Fall oder ganzzahlige Spinwerte (+-1, +-3)!!!
  if (ch>=0) {
    if (this->m_msv && (ch==0 || ch==3)) {
      //CRaScType j(anti^(this->m_dir>0)? RSP(p, -1, 0, -this->m_dir, cr, ca): RSP(p, 1, 0, this->m_dir, cr, ca));
      CRaScType j(RSP(p, r, 0, b, cr, ca));
      AddJ(CRaScType::New(j));
      // h=3 for bar vector-spinor
      //j.SetH(anti^(this->m_dir>0)?2:3);
      j.SetH(2);
      CRaScType *c(CRaScType::New(j));
      AddJ(c);
#ifdef DEBUG__BG
      msg_Debugging()<<METHOD<<"(): "<<(this->m_dir>0?'I':'O')
		     <<"0 "<<this->m_id<<" "<<j
		     <<" "<<this->m_fl<<", m = "<<m_cmass<<"\n";
#endif
      j.Test_Properties(p, j.R(), j.B());
    }
    if (ch!=3){
      //CRaScType j(anti^(this->m_dir>0)? RSPP(p, -1, 0, -this->m_dir, cr, ca): RSPP(p, 1, 0, this->m_dir, cr, ca));
      CRaScType j(RSPP(p, r, 0, b, cr, ca));
      //j.SetH(anti^(this->m_dir>0)?0:1);
      j.SetH(0);
#ifdef DEBUG__BG
      msg_Debugging()<<METHOD<<"(): "<<(this->m_dir>0?'I':'O')<<"+ "<<this->m_id
		   <<" "<<j<<" "<<(this->m_dir>0?this->m_fl.Bar():this->m_fl)
		   <<", m = "<<m_cmass<<" ("<<p.Mass()<<")\n";
#endif
      j.Test_Properties(p, j.R(), j.B());
      CRaScType *c(CRaScType::New(j));
      AddJ(c);
/*    if (p_sub) static_cast<Dipole_Color*>
      (p_sub->In().front()->Color().front())->AddJJK(c);*/
    }
  }
  if (ch<=0) {
    if (this->m_msv && (ch==0 || ch==-3)) {
      //CRaScType j(anti^(this->m_dir>0)? RSM(p, -1, 0, -this->m_dir, cr, ca): RSM(p, 1, 0, this->m_dir, cr, ca));
      CRaScType j(RSM(p, r, 0, b, cr, ca));
      AddJ(CRaScType::New(j));
      //j.SetH(anti^(this->m_dir>0)?3:2);
      j.SetH(3);
      CRaScType *c(CRaScType::New(j));
      AddJ(c);
#ifdef DEBUG__BG
      msg_Debugging()<<METHOD<<"(): "<<(this->m_dir>0?'I':'O')
		     <<"0 "<<this->m_id<<" "<<j
		     <<" "<<this->m_fl<<", m = "<<m_cmass<<"\n";
#endif
      j.Test_Properties(p, j.R(), j.B());
    }
    if (ch!=-3){
      //CRaScType j(anti^(this->m_dir>0) ? RSMM(p, -1, 0, -this->m_dir, cr, ca) : RSMM(p, 1, 0, this->m_dir, cr, ca));
      CRaScType j(RSMM(p, r, 0, b, cr, ca));
      //j.SetH(anti^(this->m_dir>0)?1:0);
      j.SetH(1);
#ifdef DEBUG__BG
      msg_Debugging()<<METHOD<<"(): "<<(this->m_dir>0?'I':'O')<<"- "<<this->m_id
		   <<" "<<j<<" "<<(this->m_dir>0?this->m_fl.Bar():this->m_fl)
		   <<", m = "<<m_cmass<<" ("<<p.Mass()<<")\n";
#endif
      j.Test_Properties(p, j.R(), j.B());
      CRaScType *c(CRaScType::New(j));
      AddJ(c);
/*    if (p_sub) static_cast<Dipole_Color*>
      (p_sub->In().front()->Color().front())->AddJJK(c);*/
    }
  }
  Test_WF_Properties(p, anti);
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
{
  for (size_t i(0);i<m_j.size();++i) {
    CRaScType_Vector *j(m_j[i].template Get<CRaScType>());
  }

}
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

// TODO: Not yet validated
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

// properties tested according to S.F.Novaes and D. Spehler, Nuclear Physics B 371 (1992) 618-636
template<typename SType>
bool CRS<SType>::Test_WF_Properties(const ATOOLS::Vec4D &p, const bool &anti) {
  CRaritaSchwinger<SType> rspp = anti?RSPP(p, -1, 1, 1, 0, 0, 1):RSPP(p, 1, 1, 1, 0, 0, 1);
  CRaritaSchwinger<SType> rsmm = anti?RSMM(p, -1, 1, 1, 0, 0, -1):RSMM(p, 1, 1, 1, 0, 0, -1);
  METOOLS::Gamma<SType> gammavec = Gamma<SType>();
  bool testresult(true);

  // properties of massless RS wave functions
  if (!this->m_msv){
    // normalization
    std::cout<<METHOD<<": Testing normalization of Rarita-Schwinger wave function..."<<std::endl;
    TCMatrix<SType> result1 = rspp.Contract4Index(rspp.Bar()) + rsmm.Contract4Index(rsmm.Bar()) + gammavec * p;
    for (size_t i(0); i<4; ++i){
      for (size_t j(0); j<4; ++j){
        if (std::abs(result1[i][j].real())> rspp.Accu()|| std::abs(result1[i][j].imag())>rspp.Accu()) {
          msg_Out()<<"Component " << i << j << " of resulting 4x4 matrix is " << result1[i][j] << ", instead of zero!"
          << std::endl;
          testresult = false;
        }
      }
    }
    if (testresult) msg_Out()<< "passed" << std::endl;

    // equality between U++ and V-- / U-- and V++ for bar and non-bar
    std::cout<<METHOD<<": Testing equality of U++ and V-- / U-- and V++ Rarita-Schwinger wave functions..."<<std::endl;
    CRaritaSchwinger<SType> rspp1 = anti?RSPP(p, 1, 1, 1, 0, 0, 1):RSPP(p, -1, 1, 1, 0, 0, 1);
    CRaritaSchwinger<SType> rsmm1 = anti?RSMM(p, 1, 1, 1, 0, 0, -1):RSMM(p, -1, 1, 1, 0, 0, -1);
    for (size_t i(0); i<16; ++i){
      if (std::abs((rspp1[i]-rsmm[i]).real()) > rsmm.Accu()|| std::abs((rspp1[i]-rsmm[i]).imag()) > rsmm.Accu()) {
        msg_Out()<<"Components " << i << " of the Rarita-Schwinger wave functions ";
        anti?(msg_Out()<<"U++ / V--" << rspp1[i] << "/" << rsmm[i]):
        (msg_Out()<<"U-- / V++" << rsmm[i] << "/" << rspp1[i]) << " are not equal!" << std::endl;
        testresult = false;
      }
      if (std::abs((rsmm1[i]-rspp[i]).real()) > rspp.Accu()|| std::abs((rsmm1[i]-rspp[i]).imag()) > rspp.Accu()) {
        msg_Out()<<"Components " << i << " of the Rarita-Schwinger wave functions ";
        anti?(msg_Out()<<"U-- / V++" << rsmm1[i] << "/" << rspp[i]):
        (msg_Out()<<"U++ / V--" << rspp[i] << "/" << rsmm1[i]) << " are not equal!" << std::endl;
        testresult = false;
      }
    }
    if (testresult) msg_Out()<< "passed" << std::endl;
  }
  else{
    // properties of massive RS wave functions
    // normalization
    std::cout<<METHOD<<": Testing normalization of Rarita-Schwinger wave function..."<<std::endl;
    CRaritaSchwinger<SType> rsp = anti?RSP(p, -1, 1, 1, 0, 0, 1):RSP(p, 1, 1, 1, 0, 0, 1);
    CRaritaSchwinger<SType> rsm = anti?RSM(p, -1, 1, 1, 0, 0, -1):RSM(p, 1, 1, 1, 0, 0, -1);
    std::vector<std::complex<SType>> result1(16);
    // result1 contains the result of contracting a bar wave function with a non-bar wave function, each entry is
    // another helicity combination of bar- and non-bar wave functions
    result1[0] = rspp.Bar()*rspp+2*rspp.R()*p.Mass();
    result1[1] = rspp.Bar()*rsp;
    result1[2] = rspp.Bar()*rsm;
    result1[3] = rspp.Bar()*rsmm;
    result1[4] = rsp.Bar()*rspp;
    result1[5] = rsp.Bar()*rsp+2*rspp.R()*p.Mass();
    result1[6] = rsp.Bar()*rsm;
    result1[7] = rsp.Bar()*rsmm;
    result1[8] = rsm.Bar()*rspp;
    result1[9] = rsm.Bar()*rsp;
    result1[10] = rsm.Bar()*rsm+2*rspp.R()*p.Mass();
    result1[11] = rsm.Bar()*rsmm;
    result1[12] = rsmm.Bar()*rspp;
    result1[13] = rsmm.Bar()*rsp;
    result1[14] = rsmm.Bar()*rsm;
    result1[15] = rsmm.Bar()*rsmm+2*rspp.R()*p.Mass();
    for (size_t i(0); i<16; ++i){
      if (std::abs(result1[i].real()) > rspp.Accu()|| std::abs(result1[i].imag()) > rspp.Accu()) {
        msg_Out()<<"Normalization of the Rarita-Schwinger wave functions" << i <<
        "do not fit! 0=++bar*++, 1=++bar*+, 2=++bar*-, 3=++bar*--, 4=+bar*++ ..." << std::endl;
        testresult = false;
      }
    }
    if (testresult) msg_Out()<< "passed" << std::endl;
    // TODO: WIP!
    // completeness relation according to Hagiwara et al. Eur. Phys. J. C (2011) 71: 1529
    SComplex propagator[4][4][4][4];
    ATOOLS::TCMatrix<SType> p_slash = gammavec*p;
    SComplex left_sum[4][4][4][4];

    // TODO: Is there a difference in the signs for anti-particles?
    ATOOLS::TCMatrix<SType> spropagator = ATOOLS::TCMatrix<SType>(p_slash + SComplex(sqrt(p.Abs2())) *
                                                                                      ATOOLS::TCMatrix<SType>(4, true));
    SComplex** lorentz_tensor = new SComplex*[4];
    for (int i=0;i<4;++i) lorentz_tensor[i] = new SComplex[4];
    //intermediate[0][0] = (*this)[0] * rs[0] + (*this)[1] * rs[1] + (*this)[2] * rs[2] + (*this)[3] * rs[3];
    //intermediate[0][1] = (*this)[0] * rs[4] + (*this)[1] * rs[5] + (*this)[2] * rs[6] + (*this)[3] * rs[7];
    for (size_t i(0); i<4; ++i){
      for (size_t j(0); j<4; ++j){
        for (size_t k(0); k<4; ++k){
          for (size_t l(0); l<4; ++l){
            lorentz_tensor[i][j] += (1.0 / 3.0)*gammavec[i][k][l]*gammavec[j][k][l];
            lorentz_tensor[i][j] -= (1.0 / (3.0*p.Abs2()))*(gammavec[i][k][l]*p[j]*p_slash[k][l]+p[i]*p_slash[k][l]*gammavec[j][k][l]);
            lorentz_tensor[i][j] += (1.0 / (3.0*p.Abs2()*p.Abs2()))*p[i]*p_slash[k][l]*p[j]*p_slash[k][l];
          }
        }
        // first summand: - metric tensor; second summand: momentum term
        lorentz_tensor[i][j] += (i == j ? (i > 0 ? 1.0 : -1.0) : 0.0) + (p[i] * p[j]) / p.Abs2();
      }
    }
    /*for (size_t i(0); i<4; ++i){
      for (size_t j(0); j<4; ++j){
        for (size_t k(0); k<4; ++k){
          for (size_t l(0); l<4; ++l){
            propagator[i][j][k][l]=spropagator[i][j]*lorentz_tensor[k][l];
            left_sum[i][j][k][l]=rspp[i][k]*rspp.Bar()[j][l] + rsp[i][k]*rsp.Bar()[j][l] + rsm[i][k]*rsm.Bar()[j][l]
                                 + rsmm[i][k]*rsmm.Bar()[j][l];
            if (std::abs((left_sum[i][j][k][l]-propagator[i][j][k][l]).real()) > rspp.Accu() ||
                std::abs((left_sum[i][j][k][l]-propagator[i][j][k][l]).imag()) > rspp.Accu()) {
              msg_Out() << "Completeness relation of the Rarita-Schwinger wave functions is not hold: "
                           "component " << i+k << j+l << "of the resulting 16 dimensional tensor is " <<
                           left_sum[i][j][k][l]-propagator[i][j][k][l] << " instead of zero!" << std::endl;
              testresult = false;

          }}
        }
      }
    }*/
    for (size_t i(0); i<4; ++i){
      for (size_t j(0); j<4; ++j){
        for (size_t k(0); k<4; ++k){
          for (size_t l(0); l<4; ++l){
            /*if (std::abs((left_sum[i+j][k]-spropagator[0][j]*lorentz_tensor[0][k]).real()) > rspp.Accu() ||
                std::abs((left_sum[i][j]-spropagator[k][l]*lorentz_tensor[i][j]).imag()) > rspp.Accu()) {
              msg_Out() << "Completeness relation of the Rarita-Schwinger wave functions is not hold: "
                           "component " << i+k << j+l << "of the resulting 16 dimensional tensor is " <<
                        left_sum[i+k][j+l]-spropagator[k][l]*lorentz_tensor[i][j] << " instead of zero!" << std::endl;
              testresult = false;
            }*/
          }
        }
      }
    }
  }

  // gauge invariance
  return testresult;
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