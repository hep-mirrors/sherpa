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
    typedef CSpinor<SType> CSpinorType;
    typedef CVec4<SType> CVec4Type;

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
    CRaScType RS(const ATOOLS::Vec4D &p, const int r, const int h, const int s, const int b, const int cr=0,
                 const int ca=0, const int hh=0, const int ms=0);
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
    bool Test_WF_Properties(const ATOOLS::Vec4D &p, const bool &anti, const int &dir);
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

// FUNCTIONS TO CALCULATE THE RARITA-SCHWINGER-WAVE-FUNCTION
// Polarisation vectors
// Representations are slightly changed compared to Spin-1 particles since calculation helicity basis is necessary
// to fit poalrisation vectors to Comix Dirac spinors (i.e. change reference vector m_k from constant default
// to -\vec{p}/|p|
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
  DEBUG_VAR(e/(sqrttwo*std::conj(m_kp_mod*pp)));
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
  DEBUG_VAR(e/(sqrttwo*std::conj(m_km_mod*pm)));
  return e/(sqrttwo*std::conj(m_km_mod*pm));
}

template <typename SType> CVec4<SType>
CRS<SType>::EMM(const Vec4D &p,const int cr,const int ca, const bool hel)
{
  Vec4D m_k_mod;
  if (hel)
    // explicit four vector form of polarisation vector for particles at rest (helicity basis not defined then)
    if (p.PSpat()==0) return CVec4Type(0., 1./ sqrt(2.), -std::complex<SType>(0., -1.)/ sqrt(2.), 0., cr, ca, 0);
    else m_k_mod=ATOOLS::Vec4D(1., -ATOOLS::Vec3D(p) / p.PSpat()); // helicity basis
  else m_k_mod = m_k; // standard setting for reference vector
  return EM(p-p.Abs2()/(2.0*m_k_mod*p)*m_k_mod,cr,ca,m_k_mod);
}

template <typename SType> CVec4<SType>
CRS<SType>::EMP(const Vec4D &p,const int cr,const int ca, const bool hel)
{
  Vec4D m_k_mod;
  // explicit four vector form of polarisation vector for particles at rest (helicity basis not defined then)
  if (hel)
    if (p.PSpat()==0) return CVec4Type(0., 1./ sqrt(2.), -std::complex<SType>(0., 1.)/ sqrt(2.), 0., cr, ca, 0);
    else m_k_mod=ATOOLS::Vec4D(1., -ATOOLS::Vec3D(p) / p.PSpat());  // helicity basis
  else m_k_mod = m_k; // standard setting for reference vector
  return EP(p-p.Abs2()/(2.0*m_k_mod*p)*m_k_mod,cr,ca,m_k_mod);
}

template <typename SType> CVec4<SType>
CRS<SType>::EML(const Vec4D &p,const int cr,const int ca, const bool hel)
{
  Vec4D m_k_mod;
  if (hel && p.PSpat()!=0) m_k_mod=ATOOLS::Vec4D(1., -ATOOLS::Vec3D(p) / p.PSpat()); // helicity basis
  else m_k_mod = m_k; // standard setting for reference vector
  double p2(p.Abs2()), a(p2/(2.0*m_k_mod*p));
  Vec4D b(p-a*m_k_mod);
  SpinorType bm(-1,b), bp(1,b), am(-1,m_k_mod), ap(1,m_k_mod);
  // if true, implement concrete four representation of Spin 1 particle at rest (helicity basis not defined then)
  CVec4Type e((hel && p.PSpat()==0.)?CVec4Type(0., 0., 0., (p[0]/fabs(p[0]))*sqrt(SComplex(4.0*p2))):VT(bp,bm)-SType(a)*VT(ap,am));
  e(0)=cr; e(1)=ca;
  e.SetH(2);
  DEBUG_VAR(e/sqrt(SComplex(4.0*p2)));
  return e/sqrt(SComplex(4.0*p2));
}

template<typename SType> CRaritaSchwinger<SType>
CRS<SType>::RS(const ATOOLS::Vec4D &p, const int r, const int h, const int s, const int b, const int cr, const int ca,
                 const int hh, const int ms) {
  CRaritaSchwinger<SType> wf;
  switch (h) {
    case 0: wf = RSPP(p, r, s, b, cr, ca, hh, ms); break;
    case 1: wf = RSMM(p, r, s, b, cr, ca, hh, ms); break;
    case 2: wf = RSP(p, r, s, b, cr, ca, hh, ms); break;
    case 3: wf = RSM(p, r, s, b, cr, ca, hh, ms); break;
    default:
      THROW(fatal_error,"Rarita-Schwinger particle has only four spin degrees of freedom!");
  }
  return wf;
}

// TODO: PLUS AND MINUS EXCHANGED
// TODO: Im HELAS-Paper hat die Wellenfunktion 18 Komponenten, wobei die letzten bei den Komponenten den Viererimpuls
//       entlang des Fermionzahlflusses enthalten -> brauchen wir das auch?
// RARITA-SCHWINGER WAVE FUNCTIONS
// - Order of indices same as in MadGraph / HELAS paper arXiv: 1010.4255 ("Four vector of Dirac spinors",
//   first four entries have same Lorentz index and ascending spinor index, fifth entry has spinor index 0 and Lorentz
//     index 1) ...
// - Reconstruction from polarisation vectors and spinors similar to
//   S.F.Novaes & D.Spehler Nuclear Physics B 371 (1992), 618-636 Eq.(13) mit Phase theta = 0
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
  // TODO: Check, whether m_r=0 for Majorana-Particles (should be inherited from spinor)
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
  if (this->m_fl.Mass()==0.0 && p[1]==0.0 && p[2]==0.0)
    this->m_p[0]=this->m_p[0]<0.0?-std::abs(this->m_p[3]):std::abs(this->m_p[3]);
  this->ResetJ();
  bool anti(this->m_fl.IsAnti());
  // mode describes whether Majorana particles are formally handled as Dirac particles or anti particles
  // should lead to the same result, mode only for consistency test purposes
  if (this->m_fl.Majorana()) anti=(mode&1)?this->m_dir<0:this->m_dir>0;
  // m_dir>0: originally incoming, m_dir<0: outgoing; for calculation all particles are handled as outgoing
  // hence: for incoming anti-particles, particle spinors are calculated and for incoming particles, anti-particle spinors
  // (anti^(this->m_dir>0))
  // TODO: !!!Richtige H-Werte?!!!
  if (ch>=0) {
    if (this->m_msv && (ch==0 || ch==3)) {
      CRaScType j(anti^(this->m_dir>0)?
                  RS(p, this->m_fl.Majorana()?-2:-1, this->m_fl.Majorana()?(mode?2:3):2, 0, -this->m_dir, cr, ca):
                  RS(p, this->m_fl.Majorana()?2:1, this->m_fl.Majorana()?(mode?3:2):2, 0, this->m_dir, cr, ca));
      AddJ(CRaScType::New(j));
      // h=3 for bar vector-spinor
      //j.SetH(anti^(this->m_dir>0)?2:3);
      j.SetH(2);
      CRaScType *c(CRaScType::New(j));
      AddJ(c);
#ifdef DEBUG__BG
      msg_Debugging()<<METHOD<<"(): "<<(this->m_dir>0?'I':'O')
		     <<"+ "<<this->m_id<<" "<<j
		     <<" "<<this->m_fl<<", m = "<<m_cmass<<"\n";
      j.Test_Properties(p, this->m_dir);
#endif
    }
    if (ch!=3){
      CRaScType j(anti^(this->m_dir>0)?
      RS(p, this->m_fl.Majorana()?-2:-1, this->m_fl.Majorana()?(mode?0:1):0, 0, -this->m_dir, cr, ca):
                  RS(p, this->m_fl.Majorana()?2:1, this->m_fl.Majorana()?(mode?1:0):0, 0, this->m_dir, cr, ca));
      //j.SetH(anti^(this->m_dir>0)?0:1);
      j.SetH(0);
#ifdef DEBUG__BG
      msg_Debugging()<<METHOD<<"(): "<<(this->m_dir>0?'I':'O')<<"++ "<<this->m_id
		   <<" "<<j<<" "<<(this->m_dir>0?this->m_fl.Bar():this->m_fl)
		   <<", m = "<<m_cmass<<" ("<<p.Mass()<<")\n";
      j.Test_Properties(p, this->m_dir);
#endif
      CRaScType *c(CRaScType::New(j));
      AddJ(c);
/*    if (p_sub) static_cast<Dipole_Color*>
      (p_sub->In().front()->Color().front())->AddJJK(c);*/
    }
  }
  if (ch<=0) {
    if (this->m_msv && (ch==0 || ch==-3)) {
      CRaScType j(anti^(this->m_dir>0)?
                  RS(p, this->m_fl.Majorana()?-2:-1, this->m_fl.Majorana()?(mode?3:2):3, 0, -this->m_dir, cr, ca):
                  RS(p, this->m_fl.Majorana()?2:1, this->m_fl.Majorana()?(mode?2:3):3, 0, this->m_dir, cr, ca));
      AddJ(CRaScType::New(j));
      //j.SetH(anti^(this->m_dir>0)?3:2);
      j.SetH(3);
      CRaScType *c(CRaScType::New(j));
      AddJ(c);
#ifdef DEBUG__BG
      msg_Debugging()<<METHOD<<"(): "<<(this->m_dir>0?'I':'O')
		     <<"- "<<this->m_id<<" "<<j
		     <<" "<<this->m_fl<<", m = "<<m_cmass<<"\n";
      j.Test_Properties(p, this->m_dir);
#endif
    }
    if (ch!=-3){
      CRaScType j(anti^(this->m_dir>0)?
                  RS(p, this->m_fl.Majorana()?-2:-1, this->m_fl.Majorana()?(mode?1:0):1, 0, -this->m_dir, cr, ca):
                  RS(p, this->m_fl.Majorana()?2:1, this->m_fl.Majorana()?(mode?0:1):1, 0, this->m_dir, cr, ca));
      //j.SetH(anti^(this->m_dir>0)?1:0);
      j.SetH(1);
#ifdef DEBUG__BG
      msg_Debugging()<<METHOD<<"(): "<<(this->m_dir>0?'I':'O')<<"-- "<<this->m_id
		   <<" "<<j<<" "<<(this->m_dir>0?this->m_fl.Bar():this->m_fl)
		   <<", m = "<<m_cmass<<" ("<<p.Mass()<<")\n";
      j.Test_Properties(p, this->m_dir);
#endif
      CRaScType *c(CRaScType::New(j));
      AddJ(c);
/*    if (p_sub) static_cast<Dipole_Color*>
      (p_sub->In().front()->Color().front())->AddJJK(c);*/
    }
  }
#ifdef DEBUG__BG
  Test_WF_Properties(p, anti, this->m_dir);
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

// TODO: Stimmt der Rest des Propagators so? Also Nenner, i im Zähler, Vorzeichen?
// TODO: Stimmt das mit unterschiedlichen Vorzeichen der Props bei pos/neg masse und B<0, B>0 Unterscheidung? (einfach
//       nur Contranktion von anderer Seite?) Ist noch etwas anderes zu beachten?
// TODO: Zählerkomponenten mit Mathematica vereinfachen, dafür mit Pythonskript (**jit) durch normalen String ohne *,
//       M_I durch i ersetzen und dann erst Mathematica übergeben; resultierender Ausdruck kann weiter vereinfacht
//       werden, minimaler Test: Vergleich mit Ergebnis des langen Ausdrucks für ein oder mehrere
//       Beispielwellenfunktionen
template <typename SType>
void CRS<SType>::AddPropagator()
{
  // add propagator for off-shell leg
  // denominator
  SComplex prop(M_I/(SType(this->m_p.Abs2())-m_cmass2));
  if (this->m_osd) prop=SComplex(M_I);
#ifdef DEBUG__BG
  msg_Debugging()<<"propagator: "<<prop
		 <<" <- p^2 = "<<this->m_p.Abs2()<<", m = "<<m_cmass<<"\n";
#endif
  Vec4D p = this->m_p;
  Complex m = m_cmass;
  Complex m2 = m_cmass2;
  for (size_t i(0);i<m_j.size();++i) {
    CRaScType_Vector *j(m_j[i].template Get<CRaScType>());
    for (typename CRaScType_Vector::iterator
           jit(j->begin());jit!=j->end();++jit) {
      CRaScType j((*jit)->R(),(*jit)->B(),(**jit)(0),(**jit)(1),
                    (*jit)->H(),(*jit)->S(),
                    ((*jit)->On()&1)<<1|((*jit)->On()&2)>>1);
      double r = m_fl.MassSign();
      // components of the propagator numerator from arXiv:1308.1668, Gl. (101) contracted with a RaSc wavefunction
      // from the right (B>0) / left (B<0), calculated by a python skript (see comment at the end of this file)
      if ((*jit)->B()>0) {// S(-p)
        j[0] = ((r*m)*(-1.0+4./(3.*m2)*p[0]*p[0]+1./3.-1./(3.*m2)*(p[0]+p[3])*p[0]-1./(3.*m2)*(p[0]-p[3])*p[0]))*(**jit)[0]
          - ((r*m)*(4./(3.*m2)*p[0]*p[1]-1./(3.*m2)*(p[0]+p[3])*p[1]-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[0]))*(**jit)[4]
          - ((r*m)*(4./(3.*m2)*p[0]*p[2]-1./(3.*m2)*(p[0]+p[3])*p[2]-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[0]))*(**jit)[8]
          - ((r*m)*(4./(3.*m2)*p[0]*p[3]-1./3.-1./(3.*m2)*(p[0]+p[3])*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[0]))*(**jit)[12]
          + ((r*m)*(-1./(3.*m2)*(p[1]-M_I*p[2])*p[0]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[0]))*(**jit)[1]
          - ((r*m)*(-1./3.-1./(3.*m2)*(p[1]-M_I*p[2])*p[1]-1./(3.*m2)*(-(p[0]-p[3]))*p[0]))*(**jit)[5]
          - ((r*m)*(M_I/3.-1./(3.*m2)*(p[1]-M_I*p[2])*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[0]))*(**jit)[9]
          - ((r*m)*(-1./(3.*m2)*(p[1]-M_I*p[2])*p[3]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[0]))*(**jit)[13]
          + ((p[0]-p[3])*(-1.0+4./(3.*m2)*p[0]*p[0]+1./3.-1./(3.*m2)*(p[0]-p[3])*p[0]-1./(3.*m2)*(p[0]+p[3])*p[0])
            + (-p[1]+M_I*p[2])*(-1./(3.*m2)*(-p[1]-M_I*p[2])*p[0]-1./(3.*m2)*(p[1]+M_I*p[2])*p[0]))*(**jit)[2]
          - ((p[0]-p[3])*(4./(3.*m2)*p[0]*p[1]-1./(3.*m2)*(p[0]-p[3])*p[1]-1./(3.*m2)*(p[1]-M_I*p[2])*p[0])
            + (-p[1]+M_I*p[2])*(1./3.-1./(3.*m2)*(-p[1]-M_I*p[2])*p[1]-1./(3.*m2)*(p[0]-p[3])*p[0]))*(**jit)[6]
          - ((p[0]-p[3])*(4./(3.*m2)*p[0]*p[2]-1./(3.*m2)*(p[0]-p[3])*p[2]-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[0])
            + (-p[1]+M_I*p[2])*(M_I/3.-1./(3.*m2)*(-p[1]-M_I*p[2])*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[0]))*(**jit)[10]
          - ((p[0]-p[3])*(4./(3.*m2)*p[0]*p[3]+1./3.-1./(3.*m2)*(p[0]-p[3])*p[3]-1./(3.*m2)*(p[0]+p[3])*p[0])
            + (-p[1]+M_I*p[2])*(-1./(3.*m2)*(-p[1]-M_I*p[2])*p[3]-1./(3.*m2)*(p[1]+M_I*p[2])*p[0]))*(**jit)[14]
          + ((p[0]-p[3])*(-1./(3.*m2)*(-p[1]+M_I*p[2])*p[0]-1./(3.*m2)*(p[1]-M_I*p[2])*p[0])
            + (-p[1]+M_I*p[2])*(-1.0+4./(3.*m2)*p[0]*p[0]+1./3.-1./(3.*m2)*(p[0]+p[3])*p[0]-1./(3.*m2)*(p[0]-p[3])*p[0]))*(**jit)[3]
          - ((p[0]-p[3])*(1./3.-1./(3.*m2)*(-p[1]+M_I*p[2])*p[1]-1./(3.*m2)*(p[0]+p[3])*p[0])
            + (-p[1]+M_I*p[2])*(4./(3.*m2)*p[0]*p[1]-1./(3.*m2)*(p[0]+p[3])*p[1]-1./(3.*m2)*(p[1]+M_I*p[2])*p[0]))*(**jit)[7]
          - ((p[0]-p[3])*(-M_I/3.-1./(3.*m2)*(-p[1]+M_I*p[2])*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[0])
            + (-p[1]+M_I*p[2])*(4./(3.*m2)*p[0]*p[2]-1./(3.*m2)*(p[0]+p[3])*p[2]-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[0]))*(**jit)[11]
          - ((p[0]-p[3])*(-1./(3.*m2)*(-p[1]+M_I*p[2])*p[3]-1./(3.*m2)*(-(p[1]-M_I*p[2]))*p[0])
            + (-p[1]+M_I*p[2])*(4./(3.*m2)*p[0]*p[3]-1./3.-1./(3.*m2)*(p[0]+p[3])*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[0]))*(**jit)[15];
        j[1] = ((r*m)*(-1./(3.*m2)*(p[1]+M_I*p[2])*p[0]-1./(3.*m2)*(-p[1]-M_I*p[2])*p[0]))*(**jit)[0]
          - ((r*m)*(-1./3.-1./(3.*m2)*(p[1]+M_I*p[2])*p[1]-1./(3.*m2)*(-(p[0]+p[3]))*p[0]))*(**jit)[4]
          - ((r*m)*(-M_I/3.-1./(3.*m2)*(p[1]+M_I*p[2])*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[0]))*(**jit)[8]
          - ((r*m)*(-1./(3.*m2)*(p[1]+M_I*p[2])*p[3]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[0]))*(**jit)[12]
          + ((r*m)*(-1.0+4./(3.*m2)*p[0]*p[0]+1./3.-1./(3.*m2)*(p[0]-p[3])*p[0]-1./(3.*m2)*(p[0]+p[3])*p[0]))*(**jit)[1]
          - ((r*m)*(4./(3.*m2)*p[0]*p[1]-1./(3.*m2)*(p[0]-p[3])*p[1]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[0]))*(**jit)[5]
          - ((r*m)*(4./(3.*m2)*p[0]*p[2]-1./(3.*m2)*(p[0]-p[3])*p[2]-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[0]))*(**jit)[9]
          - ((r*m)*(4./(3.*m2)*p[0]*p[3]+1./3.-1./(3.*m2)*(p[0]-p[3])*p[3]-1./(3.*m2)*(p[0]+p[3])*p[0]))*(**jit)[13]
          + ((-p[1]-M_I*p[2])*(-1.0+4./(3.*m2)*p[0]*p[0]+1./3.-1./(3.*m2)*(p[0]-p[3])*p[0]-1./(3.*m2)*(p[0]+p[3])*p[0])
            + (p[0]+p[3])*(-1./(3.*m2)*(-p[1]-M_I*p[2])*p[0]-1./(3.*m2)*(p[1]+M_I*p[2])*p[0]))*(**jit)[2]
          - ((-p[1]-M_I*p[2])*(4./(3.*m2)*p[0]*p[1]-1./(3.*m2)*(p[0]-p[3])*p[1]-1./(3.*m2)*(p[1]-M_I*p[2])*p[0])
            + (p[0]+p[3])*(1./3.-1./(3.*m2)*(-p[1]-M_I*p[2])*p[1]-1./(3.*m2)*(p[0]-p[3])*p[0]))*(**jit)[6]
          - ((-p[1]-M_I*p[2])*(4./(3.*m2)*p[0]*p[2]-1./(3.*m2)*(p[0]-p[3])*p[2]-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[0])
            + (p[0]+p[3])*(M_I/3.-1./(3.*m2)*(-p[1]-M_I*p[2])*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[0]))*(**jit)[10]
          - ((-p[1]-M_I*p[2])*(4./(3.*m2)*p[0]*p[3]+1./3.-1./(3.*m2)*(p[0]-p[3])*p[3]-1./(3.*m2)*(p[0]+p[3])*p[0])
            + (p[0]+p[3])*(-1./(3.*m2)*(-p[1]-M_I*p[2])*p[3]-1./(3.*m2)*(p[1]+M_I*p[2])*p[0]))*(**jit)[14]
          + ((-p[1]-M_I*p[2])*(-1./(3.*m2)*(-p[1]+M_I*p[2])*p[0]-1./(3.*m2)*(p[1]-M_I*p[2])*p[0])
            + (p[0]+p[3])*(-1.0+4./(3.*m2)*p[0]*p[0]+1./3.-1./(3.*m2)*(p[0]+p[3])*p[0]-1./(3.*m2)*(p[0]-p[3])*p[0]))*(**jit)[3]
          - ((-p[1]-M_I*p[2])*(1./3.-1./(3.*m2)*(-p[1]+M_I*p[2])*p[1]-1./(3.*m2)*(p[0]+p[3])*p[0])
            + (p[0]+p[3])*(4./(3.*m2)*p[0]*p[1]-1./(3.*m2)*(p[0]+p[3])*p[1]-1./(3.*m2)*(p[1]+M_I*p[2])*p[0]))*(**jit)[7]
          - ((-p[1]-M_I*p[2])*(-M_I/3.-1./(3.*m2)*(-p[1]+M_I*p[2])*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[0])
            + (p[0]+p[3])*(4./(3.*m2)*p[0]*p[2]-1./(3.*m2)*(p[0]+p[3])*p[2]-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[0]))*(**jit)[11]
          - ((-p[1]-M_I*p[2])*(-1./(3.*m2)*(-p[1]+M_I*p[2])*p[3]-1./(3.*m2)*(-(p[1]-M_I*p[2]))*p[0])
            + (p[0]+p[3])*(4./(3.*m2)*p[0]*p[3]-1./3.-1./(3.*m2)*(p[0]+p[3])*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[0]))*(**jit)[15];
        j[2] = ((p[0]+p[3])*(-1.0+4./(3.*m2)*p[0]*p[0]+1./3.-1./(3.*m2)*(p[0]+p[3])*p[0]-1./(3.*m2)*(p[0]-p[3])*p[0])
            + (p[1]-M_I*p[2])*(-1./(3.*m2)*(p[1]+M_I*p[2])*p[0]-1./(3.*m2)*(-p[1]-M_I*p[2])*p[0]))*(**jit)[0]
          - ((p[0]+p[3])*(4./(3.*m2)*p[0]*p[1]-1./(3.*m2)*(p[0]+p[3])*p[1]-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[0])
            + (p[1]-M_I*p[2])*(-1./3.-1./(3.*m2)*(p[1]+M_I*p[2])*p[1]-1./(3.*m2)*(-(p[0]+p[3]))*p[0]))*(**jit)[4]
          - ((p[0]+p[3])*(4./(3.*m2)*p[0]*p[2]-1./(3.*m2)*(p[0]+p[3])*p[2]-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[0])
            + (p[1]-M_I*p[2])*(-M_I/3.-1./(3.*m2)*(p[1]+M_I*p[2])*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[0]))*(**jit)[8]
          - ((p[0]+p[3])*(4./(3.*m2)*p[0]*p[3]-1./3.-1./(3.*m2)*(p[0]+p[3])*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[0])
            + (p[1]-M_I*p[2])*(-1./(3.*m2)*(p[1]+M_I*p[2])*p[3]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[0]))*(**jit)[12]
          + ((p[0]+p[3])*(-1./(3.*m2)*(p[1]-M_I*p[2])*p[0]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[0])
            + (p[1]-M_I*p[2])*(-1.0+4./(3.*m2)*p[0]*p[0]+1./3.-1./(3.*m2)*(p[0]-p[3])*p[0]-1./(3.*m2)*(p[0]+p[3])*p[0]))*(**jit)[1]
          - ((p[0]+p[3])*(-1./3.-1./(3.*m2)*(p[1]-M_I*p[2])*p[1]-1./(3.*m2)*(-(p[0]-p[3]))*p[0])
            + (p[1]-M_I*p[2])*(4./(3.*m2)*p[0]*p[1]-1./(3.*m2)*(p[0]-p[3])*p[1]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[0]))*(**jit)[5]
          - ((p[0]+p[3])*(M_I/3.-1./(3.*m2)*(p[1]-M_I*p[2])*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[0])
            + (p[1]-M_I*p[2])*(4./(3.*m2)*p[0]*p[2]-1./(3.*m2)*(p[0]-p[3])*p[2]-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[0]))*(**jit)[9]
          - ((p[0]+p[3])*(-1./(3.*m2)*(p[1]-M_I*p[2])*p[3]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[0])
            + (p[1]-M_I*p[2])*(4./(3.*m2)*p[0]*p[3]+1./3.-1./(3.*m2)*(p[0]-p[3])*p[3]-1./(3.*m2)*(p[0]+p[3])*p[0]))*(**jit)[13]
          + ((r*m)*(-1.0+4./(3.*m2)*p[0]*p[0]+1./3.-1./(3.*m2)*(p[0]-p[3])*p[0]-1./(3.*m2)*(p[0]+p[3])*p[0]))*(**jit)[2]
          - ((r*m)*(4./(3.*m2)*p[0]*p[1]-1./(3.*m2)*(p[0]-p[3])*p[1]-1./(3.*m2)*(p[1]-M_I*p[2])*p[0]))*(**jit)[6]
          - ((r*m)*(4./(3.*m2)*p[0]*p[2]-1./(3.*m2)*(p[0]-p[3])*p[2]-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[0]))*(**jit)[10]
          - ((r*m)*(4./(3.*m2)*p[0]*p[3]+1./3.-1./(3.*m2)*(p[0]-p[3])*p[3]-1./(3.*m2)*(p[0]+p[3])*p[0]))*(**jit)[14]
          + ((r*m)*(-1./(3.*m2)*(-p[1]+M_I*p[2])*p[0]-1./(3.*m2)*(p[1]-M_I*p[2])*p[0]))*(**jit)[3]
          - ((r*m)*(1./3.-1./(3.*m2)*(-p[1]+M_I*p[2])*p[1]-1./(3.*m2)*(p[0]+p[3])*p[0]))*(**jit)[7]
          - ((r*m)*(-M_I/3.-1./(3.*m2)*(-p[1]+M_I*p[2])*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[0]))*(**jit)[11]
          - ((r*m)*(-1./(3.*m2)*(-p[1]+M_I*p[2])*p[3]-1./(3.*m2)*(-(p[1]-M_I*p[2]))*p[0]))*(**jit)[15];
        j[3] = ((p[1]+M_I*p[2])*(-1.0+4./(3.*m2)*p[0]*p[0]+1./3.-1./(3.*m2)*(p[0]+p[3])*p[0]-1./(3.*m2)*(p[0]-p[3])*p[0])
            + (p[0]-p[3])*(-1./(3.*m2)*(p[1]+M_I*p[2])*p[0]-1./(3.*m2)*(-p[1]-M_I*p[2])*p[0]))*(**jit)[0]
          - ((p[1]+M_I*p[2])*(4./(3.*m2)*p[0]*p[1]-1./(3.*m2)*(p[0]+p[3])*p[1]-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[0])
            + (p[0]-p[3])*(-1./3.-1./(3.*m2)*(p[1]+M_I*p[2])*p[1]-1./(3.*m2)*(-(p[0]+p[3]))*p[0]))*(**jit)[4]
          - ((p[1]+M_I*p[2])*(4./(3.*m2)*p[0]*p[2]-1./(3.*m2)*(p[0]+p[3])*p[2]-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[0])
            + (p[0]-p[3])*(-M_I/3.-1./(3.*m2)*(p[1]+M_I*p[2])*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[0]))*(**jit)[8]
          - ((p[1]+M_I*p[2])*(4./(3.*m2)*p[0]*p[3]-1./3.-1./(3.*m2)*(p[0]+p[3])*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[0])
            + (p[0]-p[3])*(-1./(3.*m2)*(p[1]+M_I*p[2])*p[3]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[0]))*(**jit)[12]
          + ((p[1]+M_I*p[2])*(-1./(3.*m2)*(p[1]-M_I*p[2])*p[0]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[0])
            + (p[0]-p[3])*(-1.0+4./(3.*m2)*p[0]*p[0]+1./3.-1./(3.*m2)*(p[0]-p[3])*p[0]-1./(3.*m2)*(p[0]+p[3])*p[0]))*(**jit)[1]
          - ((p[1]+M_I*p[2])*(-1./3.-1./(3.*m2)*(p[1]-M_I*p[2])*p[1]-1./(3.*m2)*(-(p[0]-p[3]))*p[0])
            + (p[0]-p[3])*(4./(3.*m2)*p[0]*p[1]-1./(3.*m2)*(p[0]-p[3])*p[1]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[0]))*(**jit)[5]
          - ((p[1]+M_I*p[2])*(M_I/3.-1./(3.*m2)*(p[1]-M_I*p[2])*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[0])
            + (p[0]-p[3])*(4./(3.*m2)*p[0]*p[2]-1./(3.*m2)*(p[0]-p[3])*p[2]-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[0]))*(**jit)[9]
          - ((p[1]+M_I*p[2])*(-1./(3.*m2)*(p[1]-M_I*p[2])*p[3]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[0])
          + (p[0]-p[3])*(4./(3.*m2)*p[0]*p[3]+1./3.-1./(3.*m2)*(p[0]-p[3])*p[3]-1./(3.*m2)*(p[0]+p[3])*p[0]))*(**jit)[13]
          + ((r*m)*(-1./(3.*m2)*(-p[1]-M_I*p[2])*p[0]-1./(3.*m2)*(p[1]+M_I*p[2])*p[0]))*(**jit)[2]
          - ((r*m)*(1./3.-1./(3.*m2)*(-p[1]-M_I*p[2])*p[1]-1./(3.*m2)*(p[0]-p[3])*p[0]))*(**jit)[6]
          - ((r*m)*(M_I/3.-1./(3.*m2)*(-p[1]-M_I*p[2])*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[0]))*(**jit)[10]
          - ((r*m)*(-1./(3.*m2)*(-p[1]-M_I*p[2])*p[3]-1./(3.*m2)*(p[1]+M_I*p[2])*p[0]))*(**jit)[14]
          + ((r*m)*(-1.0+4./(3.*m2)*p[0]*p[0]+1./3.-1./(3.*m2)*(p[0]+p[3])*p[0]-1./(3.*m2)*(p[0]-p[3])*p[0]))*(**jit)[3]
          - ((r*m)*(4./(3.*m2)*p[0]*p[1]-1./(3.*m2)*(p[0]+p[3])*p[1]-1./(3.*m2)*(p[1]+M_I*p[2])*p[0]))*(**jit)[7]
          - ((r*m)*(4./(3.*m2)*p[0]*p[2]-1./(3.*m2)*(p[0]+p[3])*p[2]-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[0]))*(**jit)[11]
          - ((r*m)*(4./(3.*m2)*p[0]*p[3]-1./3.-1./(3.*m2)*(p[0]+p[3])*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[0]))*(**jit)[15];
        j[4] = ((r*m)*(4./(3.*m2)*p[1]*p[0]-1./(3.*m2)*(p[1]+M_I*p[2])*p[0]-1./(3.*m2)*(p[0]-p[3])*p[1]))*(**jit)[0] - ((r*m)*(1.0+4./(3.*m2)*p[1]*p[1]-1./3.-1./(3.*m2)*(p[1]+M_I*p[2])*p[1]-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[1]))*(**jit)[4] - ((r*m)*(4./(3.*m2)*p[1]*p[2]-M_I/3.-1./(3.*m2)*(p[1]+M_I*p[2])*p[2]-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[1]))*(**jit)[8] - ((r*m)*(4./(3.*m2)*p[1]*p[3]-1./(3.*m2)*(p[1]+M_I*p[2])*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[1]))*(**jit)[12] + ((r*m)*(1./3.-1./(3.*m2)*(p[0]-p[3])*p[0]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[1]))*(**jit)[1] - ((r*m)*(-1./(3.*m2)*(p[0]-p[3])*p[1]-1./(3.*m2)*(-(p[0]-p[3]))*p[1]))*(**jit)[5] - ((r*m)*(-1./(3.*m2)*(p[0]-p[3])*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[1]))*(**jit)[9] - ((r*m)*(1./3.-1./(3.*m2)*(p[0]-p[3])*p[3]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[1]))*(**jit)[13] + ((p[0]-p[3])*(4./(3.*m2)*p[1]*p[0]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[0]-1./(3.*m2)*(p[0]+p[3])*p[1]) + (-p[1]+M_I*p[2])*(-1./3.-1./(3.*m2)*(-(p[0]-p[3]))*p[0]-1./(3.*m2)*(p[1]+M_I*p[2])*p[1]))*(**jit)[2] - ((p[0]-p[3])*(1.0+4./(3.*m2)*p[1]*p[1]-1./3.-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[1]-1./(3.*m2)*(p[1]-M_I*p[2])*p[1]) + (-p[1]+M_I*p[2])*(-1./(3.*m2)*(-(p[0]-p[3]))*p[1]-1./(3.*m2)*(p[0]-p[3])*p[1]))*(**jit)[6] - ((p[0]-p[3])*(4./(3.*m2)*p[1]*p[2]-M_I/3.-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[2]-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[1]) + (-p[1]+M_I*p[2])*(-1./(3.*m2)*(-(p[0]-p[3]))*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[1]))*(**jit)[10] - ((p[0]-p[3])*(4./(3.*m2)*p[1]*p[3]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[3]-1./(3.*m2)*(p[0]+p[3])*p[1]) + (-p[1]+M_I*p[2])*(-1./3.-1./(3.*m2)*(-(p[0]-p[3]))*p[3]-1./(3.*m2)*(p[1]+M_I*p[2])*p[1]))*(**jit)[14] + ((p[0]-p[3])*(-1./3.-1./(3.*m2)*(-(p[0]+p[3]))*p[0]-1./(3.*m2)*(p[1]-M_I*p[2])*p[1]) + (-p[1]+M_I*p[2])*(4./(3.*m2)*p[1]*p[0]-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[0]-1./(3.*m2)*(p[0]-p[3])*p[1]))*(**jit)[3] - ((p[0]-p[3])*(-1./(3.*m2)*(-(p[0]+p[3]))*p[1]-1./(3.*m2)*(p[0]+p[3])*p[1]) + (-p[1]+M_I*p[2])*(1.0+4./(3.*m2)*p[1]*p[1]-1./3.-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[1]-1./(3.*m2)*(p[1]+M_I*p[2])*p[1]))*(**jit)[7] - ((p[0]-p[3])*(-1./(3.*m2)*(-(p[0]+p[3]))*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[1]) + (-p[1]+M_I*p[2])*(4./(3.*m2)*p[1]*p[2]+M_I/3.-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[2]-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[1]))*(**jit)[11] - ((p[0]-p[3])*(1./3.-1./(3.*m2)*(-(p[0]+p[3]))*p[3]-1./(3.*m2)*(-(p[1]-M_I*p[2]))*p[1]) + (-p[1]+M_I*p[2])*(4./(3.*m2)*p[1]*p[3]-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[1]))*(**jit)[15];
        j[5] = ((r*m)*(1./3.-1./(3.*m2)*(p[0]+p[3])*p[0]-1./(3.*m2)*(-p[1]-M_I*p[2])*p[1]))*(**jit)[0] - ((r*m)*(-1./(3.*m2)*(p[0]+p[3])*p[1]-1./(3.*m2)*(-(p[0]+p[3]))*p[1]))*(**jit)[4] - ((r*m)*(-1./(3.*m2)*(p[0]+p[3])*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[1]))*(**jit)[8] - ((r*m)*(-1./3.-1./(3.*m2)*(p[0]+p[3])*p[3]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[1]))*(**jit)[12] + ((r*m)*(4./(3.*m2)*p[1]*p[0]-1./(3.*m2)*(p[1]-M_I*p[2])*p[0]-1./(3.*m2)*(p[0]+p[3])*p[1]))*(**jit)[1] - ((r*m)*(1.0+4./(3.*m2)*p[1]*p[1]-1./3.-1./(3.*m2)*(p[1]-M_I*p[2])*p[1]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[1]))*(**jit)[5] - ((r*m)*(4./(3.*m2)*p[1]*p[2]+M_I/3.-1./(3.*m2)*(p[1]-M_I*p[2])*p[2]-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[1]))*(**jit)[9] - ((r*m)*(4./(3.*m2)*p[1]*p[3]-1./(3.*m2)*(p[1]-M_I*p[2])*p[3]-1./(3.*m2)*(p[0]+p[3])*p[1]))*(**jit)[13] + ((-p[1]-M_I*p[2])*(4./(3.*m2)*p[1]*p[0]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[0]-1./(3.*m2)*(p[0]+p[3])*p[1]) + (p[0]+p[3])*(-1./3.-1./(3.*m2)*(-(p[0]-p[3]))*p[0]-1./(3.*m2)*(p[1]+M_I*p[2])*p[1]))*(**jit)[2] - ((-p[1]-M_I*p[2])*(1.0+4./(3.*m2)*p[1]*p[1]-1./3.-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[1]-1./(3.*m2)*(p[1]-M_I*p[2])*p[1]) + (p[0]+p[3])*(-1./(3.*m2)*(-(p[0]-p[3]))*p[1]-1./(3.*m2)*(p[0]-p[3])*p[1]))*(**jit)[6] - ((-p[1]-M_I*p[2])*(4./(3.*m2)*p[1]*p[2]-M_I/3.-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[2]-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[1]) + (p[0]+p[3])*(-1./(3.*m2)*(-(p[0]-p[3]))*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[1]))*(**jit)[10] - ((-p[1]-M_I*p[2])*(4./(3.*m2)*p[1]*p[3]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[3]-1./(3.*m2)*(p[0]+p[3])*p[1]) + (p[0]+p[3])*(-1./3.-1./(3.*m2)*(-(p[0]-p[3]))*p[3]-1./(3.*m2)*(p[1]+M_I*p[2])*p[1]))*(**jit)[14] + ((-p[1]-M_I*p[2])*(-1./3.-1./(3.*m2)*(-(p[0]+p[3]))*p[0]-1./(3.*m2)*(p[1]-M_I*p[2])*p[1]) + (p[0]+p[3])*(4./(3.*m2)*p[1]*p[0]-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[0]-1./(3.*m2)*(p[0]-p[3])*p[1]))*(**jit)[3] - ((-p[1]-M_I*p[2])*(-1./(3.*m2)*(-(p[0]+p[3]))*p[1]-1./(3.*m2)*(p[0]+p[3])*p[1]) + (p[0]+p[3])*(1.0+4./(3.*m2)*p[1]*p[1]-1./3.-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[1]-1./(3.*m2)*(p[1]+M_I*p[2])*p[1]))*(**jit)[7] - ((-p[1]-M_I*p[2])*(-1./(3.*m2)*(-(p[0]+p[3]))*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[1]) + (p[0]+p[3])*(4./(3.*m2)*p[1]*p[2]+M_I/3.-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[2]-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[1]))*(**jit)[11] - ((-p[1]-M_I*p[2])*(1./3.-1./(3.*m2)*(-(p[0]+p[3]))*p[3]-1./(3.*m2)*(-(p[1]-M_I*p[2]))*p[1]) + (p[0]+p[3])*(4./(3.*m2)*p[1]*p[3]-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[1]))*(**jit)[15];
        j[6] = ((p[0]+p[3])*(4./(3.*m2)*p[1]*p[0]-1./(3.*m2)*(p[1]+M_I*p[2])*p[0]-1./(3.*m2)*(p[0]-p[3])*p[1]) + (p[1]-M_I*p[2])*(1./3.-1./(3.*m2)*(p[0]+p[3])*p[0]-1./(3.*m2)*(-p[1]-M_I*p[2])*p[1]))*(**jit)[0] - ((p[0]+p[3])*(1.0+4./(3.*m2)*p[1]*p[1]-1./3.-1./(3.*m2)*(p[1]+M_I*p[2])*p[1]-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[1]) + (p[1]-M_I*p[2])*(-1./(3.*m2)*(p[0]+p[3])*p[1]-1./(3.*m2)*(-(p[0]+p[3]))*p[1]))*(**jit)[4] - ((p[0]+p[3])*(4./(3.*m2)*p[1]*p[2]-M_I/3.-1./(3.*m2)*(p[1]+M_I*p[2])*p[2]-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[1]) + (p[1]-M_I*p[2])*(-1./(3.*m2)*(p[0]+p[3])*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[1]))*(**jit)[8] - ((p[0]+p[3])*(4./(3.*m2)*p[1]*p[3]-1./(3.*m2)*(p[1]+M_I*p[2])*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[1]) + (p[1]-M_I*p[2])*(-1./3.-1./(3.*m2)*(p[0]+p[3])*p[3]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[1]))*(**jit)[12] + ((p[0]+p[3])*(1./3.-1./(3.*m2)*(p[0]-p[3])*p[0]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[1]) + (p[1]-M_I*p[2])*(4./(3.*m2)*p[1]*p[0]-1./(3.*m2)*(p[1]-M_I*p[2])*p[0]-1./(3.*m2)*(p[0]+p[3])*p[1]))*(**jit)[1] - ((p[0]+p[3])*(-1./(3.*m2)*(p[0]-p[3])*p[1]-1./(3.*m2)*(-(p[0]-p[3]))*p[1]) + (p[1]-M_I*p[2])*(1.0+4./(3.*m2)*p[1]*p[1]-1./3.-1./(3.*m2)*(p[1]-M_I*p[2])*p[1]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[1]))*(**jit)[5] - ((p[0]+p[3])*(-1./(3.*m2)*(p[0]-p[3])*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[1]) + (p[1]-M_I*p[2])*(4./(3.*m2)*p[1]*p[2]+M_I/3.-1./(3.*m2)*(p[1]-M_I*p[2])*p[2]-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[1]))*(**jit)[9] - ((p[0]+p[3])*(1./3.-1./(3.*m2)*(p[0]-p[3])*p[3]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[1]) + (p[1]-M_I*p[2])*(4./(3.*m2)*p[1]*p[3]-1./(3.*m2)*(p[1]-M_I*p[2])*p[3]-1./(3.*m2)*(p[0]+p[3])*p[1]))*(**jit)[13] + ((r*m)*(4./(3.*m2)*p[1]*p[0]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[0]-1./(3.*m2)*(p[0]+p[3])*p[1]))*(**jit)[2] - ((r*m)*(1.0+4./(3.*m2)*p[1]*p[1]-1./3.-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[1]-1./(3.*m2)*(p[1]-M_I*p[2])*p[1]))*(**jit)[6] - ((r*m)*(4./(3.*m2)*p[1]*p[2]-M_I/3.-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[2]-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[1]))*(**jit)[10] - ((r*m)*(4./(3.*m2)*p[1]*p[3]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[3]-1./(3.*m2)*(p[0]+p[3])*p[1]))*(**jit)[14] + ((r*m)*(-1./3.-1./(3.*m2)*(-(p[0]+p[3]))*p[0]-1./(3.*m2)*(p[1]-M_I*p[2])*p[1]))*(**jit)[3] - ((r*m)*(-1./(3.*m2)*(-(p[0]+p[3]))*p[1]-1./(3.*m2)*(p[0]+p[3])*p[1]))*(**jit)[7] - ((r*m)*(-1./(3.*m2)*(-(p[0]+p[3]))*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[1]))*(**jit)[11] - ((r*m)*(1./3.-1./(3.*m2)*(-(p[0]+p[3]))*p[3]-1./(3.*m2)*(-(p[1]-M_I*p[2]))*p[1]))*(**jit)[15];
        j[7] = ((p[1]+M_I*p[2])*(4./(3.*m2)*p[1]*p[0]-1./(3.*m2)*(p[1]+M_I*p[2])*p[0]-1./(3.*m2)*(p[0]-p[3])*p[1]) + (p[0]-p[3])*(1./3.-1./(3.*m2)*(p[0]+p[3])*p[0]-1./(3.*m2)*(-p[1]-M_I*p[2])*p[1]))*(**jit)[0] - ((p[1]+M_I*p[2])*(1.0+4./(3.*m2)*p[1]*p[1]-1./3.-1./(3.*m2)*(p[1]+M_I*p[2])*p[1]-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[1]) + (p[0]-p[3])*(-1./(3.*m2)*(p[0]+p[3])*p[1]-1./(3.*m2)*(-(p[0]+p[3]))*p[1]))*(**jit)[4] - ((p[1]+M_I*p[2])*(4./(3.*m2)*p[1]*p[2]-M_I/3.-1./(3.*m2)*(p[1]+M_I*p[2])*p[2]-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[1]) + (p[0]-p[3])*(-1./(3.*m2)*(p[0]+p[3])*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[1]))*(**jit)[8] - ((p[1]+M_I*p[2])*(4./(3.*m2)*p[1]*p[3]-1./(3.*m2)*(p[1]+M_I*p[2])*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[1]) + (p[0]-p[3])*(-1./3.-1./(3.*m2)*(p[0]+p[3])*p[3]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[1]))*(**jit)[12] + ((p[1]+M_I*p[2])*(1./3.-1./(3.*m2)*(p[0]-p[3])*p[0]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[1]) + (p[0]-p[3])*(4./(3.*m2)*p[1]*p[0]-1./(3.*m2)*(p[1]-M_I*p[2])*p[0]-1./(3.*m2)*(p[0]+p[3])*p[1]))*(**jit)[1] - ((p[1]+M_I*p[2])*(-1./(3.*m2)*(p[0]-p[3])*p[1]-1./(3.*m2)*(-(p[0]-p[3]))*p[1]) + (p[0]-p[3])*(1.0+4./(3.*m2)*p[1]*p[1]-1./3.-1./(3.*m2)*(p[1]-M_I*p[2])*p[1]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[1]))*(**jit)[5] - ((p[1]+M_I*p[2])*(-1./(3.*m2)*(p[0]-p[3])*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[1]) + (p[0]-p[3])*(4./(3.*m2)*p[1]*p[2]+M_I/3.-1./(3.*m2)*(p[1]-M_I*p[2])*p[2]-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[1]))*(**jit)[9] - ((p[1]+M_I*p[2])*(1./3.-1./(3.*m2)*(p[0]-p[3])*p[3]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[1]) + (p[0]-p[3])*(4./(3.*m2)*p[1]*p[3]-1./(3.*m2)*(p[1]-M_I*p[2])*p[3]-1./(3.*m2)*(p[0]+p[3])*p[1]))*(**jit)[13] + ((r*m)*(-1./3.-1./(3.*m2)*(-(p[0]-p[3]))*p[0]-1./(3.*m2)*(p[1]+M_I*p[2])*p[1]))*(**jit)[2] - ((r*m)*(-1./(3.*m2)*(-(p[0]-p[3]))*p[1]-1./(3.*m2)*(p[0]-p[3])*p[1]))*(**jit)[6] - ((r*m)*(-1./(3.*m2)*(-(p[0]-p[3]))*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[1]))*(**jit)[10] - ((r*m)*(-1./3.-1./(3.*m2)*(-(p[0]-p[3]))*p[3]-1./(3.*m2)*(p[1]+M_I*p[2])*p[1]))*(**jit)[14] + ((r*m)*(4./(3.*m2)*p[1]*p[0]-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[0]-1./(3.*m2)*(p[0]-p[3])*p[1]))*(**jit)[3] - ((r*m)*(1.0+4./(3.*m2)*p[1]*p[1]-1./3.-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[1]-1./(3.*m2)*(p[1]+M_I*p[2])*p[1]))*(**jit)[7] - ((r*m)*(4./(3.*m2)*p[1]*p[2]+M_I/3.-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[2]-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[1]))*(**jit)[11] - ((r*m)*(4./(3.*m2)*p[1]*p[3]-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[1]))*(**jit)[15];
        j[8] = ((r*m)*(4./(3.*m2)*p[2]*p[0]-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[0]-1./(3.*m2)*(p[0]-p[3])*p[2]))*(**jit)[0] - ((r*m)*(4./(3.*m2)*p[2]*p[1]+M_I/3.-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[1]-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[2]))*(**jit)[4] - ((r*m)*(1.0+4./(3.*m2)*p[2]*p[2]-1./3.-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[2]-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[2]))*(**jit)[8] - ((r*m)*(4./(3.*m2)*p[2]*p[3]-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[2]))*(**jit)[12] + ((r*m)*(-M_I/3.-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[0]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[2]))*(**jit)[1] - ((r*m)*(-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[1]-1./(3.*m2)*(-(p[0]-p[3]))*p[2]))*(**jit)[5] - ((r*m)*(-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[2]))*(**jit)[9] - ((r*m)*(-M_I/3.-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[3]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[2]))*(**jit)[13] + ((p[0]-p[3])*(4./(3.*m2)*p[2]*p[0]-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[0]-1./(3.*m2)*(p[0]+p[3])*p[2]) + (-p[1]+M_I*p[2])*(-M_I/3.-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[0]-1./(3.*m2)*(p[1]+M_I*p[2])*p[2]))*(**jit)[2] - ((p[0]-p[3])*(4./(3.*m2)*p[2]*p[1]+M_I/3.-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[1]-1./(3.*m2)*(p[1]-M_I*p[2])*p[2]) + (-p[1]+M_I*p[2])*(-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[1]-1./(3.*m2)*(p[0]-p[3])*p[2]))*(**jit)[6] - ((p[0]-p[3])*(1.0+4./(3.*m2)*p[2]*p[2]-1./3.-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[2]-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[2]) + (-p[1]+M_I*p[2])*(-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[2]))*(**jit)[10] - ((p[0]-p[3])*(4./(3.*m2)*p[2]*p[3]-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[3]-1./(3.*m2)*(p[0]+p[3])*p[2]) + (-p[1]+M_I*p[2])*(-M_I/3.-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[3]-1./(3.*m2)*(p[1]+M_I*p[2])*p[2]))*(**jit)[14] + ((p[0]-p[3])*(M_I/3.-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[0]-1./(3.*m2)*(p[1]-M_I*p[2])*p[2]) + (-p[1]+M_I*p[2])*(4./(3.*m2)*p[2]*p[0]-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[0]-1./(3.*m2)*(p[0]-p[3])*p[2]))*(**jit)[3] - ((p[0]-p[3])*(-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[1]-1./(3.*m2)*(p[0]+p[3])*p[2]) + (-p[1]+M_I*p[2])*(4./(3.*m2)*p[2]*p[1]-M_I/3.-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[1]-1./(3.*m2)*(p[1]+M_I*p[2])*p[2]))*(**jit)[7] - ((p[0]-p[3])*(-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[2]) + (-p[1]+M_I*p[2])*(1.0+4./(3.*m2)*p[2]*p[2]-1./3.-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[2]-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[2]))*(**jit)[11] - ((p[0]-p[3])*(-M_I/3.-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[3]-1./(3.*m2)*(-(p[1]-M_I*p[2]))*p[2]) + (-p[1]+M_I*p[2])*(4./(3.*m2)*p[2]*p[3]-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[2]))*(**jit)[15];
        j[9] = ((r*m)*(M_I/3.-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[0]-1./(3.*m2)*(-p[1]-M_I*p[2])*p[2]))*(**jit)[0] - ((r*m)*(-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[1]-1./(3.*m2)*(-(p[0]+p[3]))*p[2]))*(**jit)[4] - ((r*m)*(-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[2]))*(**jit)[8] - ((r*m)*(-M_I/3.-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[3]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[2]))*(**jit)[12] + ((r*m)*(4./(3.*m2)*p[2]*p[0]-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[0]-1./(3.*m2)*(p[0]+p[3])*p[2]))*(**jit)[1] - ((r*m)*(4./(3.*m2)*p[2]*p[1]-M_I/3.-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[1]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[2]))*(**jit)[5] - ((r*m)*(1.0+4./(3.*m2)*p[2]*p[2]-1./3.-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[2]-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[2]))*(**jit)[9] - ((r*m)*(4./(3.*m2)*p[2]*p[3]-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[3]-1./(3.*m2)*(p[0]+p[3])*p[2]))*(**jit)[13] + ((-p[1]-M_I*p[2])*(4./(3.*m2)*p[2]*p[0]-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[0]-1./(3.*m2)*(p[0]+p[3])*p[2]) + (p[0]+p[3])*(-M_I/3.-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[0]-1./(3.*m2)*(p[1]+M_I*p[2])*p[2]))*(**jit)[2] - ((-p[1]-M_I*p[2])*(4./(3.*m2)*p[2]*p[1]+M_I/3.-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[1]-1./(3.*m2)*(p[1]-M_I*p[2])*p[2]) + (p[0]+p[3])*(-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[1]-1./(3.*m2)*(p[0]-p[3])*p[2]))*(**jit)[6] - ((-p[1]-M_I*p[2])*(1.0+4./(3.*m2)*p[2]*p[2]-1./3.-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[2]-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[2]) + (p[0]+p[3])*(-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[2]))*(**jit)[10] - ((-p[1]-M_I*p[2])*(4./(3.*m2)*p[2]*p[3]-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[3]-1./(3.*m2)*(p[0]+p[3])*p[2]) + (p[0]+p[3])*(-M_I/3.-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[3]-1./(3.*m2)*(p[1]+M_I*p[2])*p[2]))*(**jit)[14] + ((-p[1]-M_I*p[2])*(M_I/3.-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[0]-1./(3.*m2)*(p[1]-M_I*p[2])*p[2]) + (p[0]+p[3])*(4./(3.*m2)*p[2]*p[0]-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[0]-1./(3.*m2)*(p[0]-p[3])*p[2]))*(**jit)[3] - ((-p[1]-M_I*p[2])*(-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[1]-1./(3.*m2)*(p[0]+p[3])*p[2]) + (p[0]+p[3])*(4./(3.*m2)*p[2]*p[1]-M_I/3.-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[1]-1./(3.*m2)*(p[1]+M_I*p[2])*p[2]))*(**jit)[7] - ((-p[1]-M_I*p[2])*(-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[2]) + (p[0]+p[3])*(1.0+4./(3.*m2)*p[2]*p[2]-1./3.-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[2]-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[2]))*(**jit)[11] - ((-p[1]-M_I*p[2])*(-M_I/3.-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[3]-1./(3.*m2)*(-(p[1]-M_I*p[2]))*p[2]) + (p[0]+p[3])*(4./(3.*m2)*p[2]*p[3]-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[2]))*(**jit)[15];
        j[10] = ((p[0]+p[3])*(4./(3.*m2)*p[2]*p[0]-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[0]-1./(3.*m2)*(p[0]-p[3])*p[2]) + (p[1]-M_I*p[2])*(M_I/3.-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[0]-1./(3.*m2)*(-p[1]-M_I*p[2])*p[2]))*(**jit)[0] - ((p[0]+p[3])*(4./(3.*m2)*p[2]*p[1]+M_I/3.-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[1]-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[2]) + (p[1]-M_I*p[2])*(-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[1]-1./(3.*m2)*(-(p[0]+p[3]))*p[2]))*(**jit)[4] - ((p[0]+p[3])*(1.0+4./(3.*m2)*p[2]*p[2]-1./3.-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[2]-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[2]) + (p[1]-M_I*p[2])*(-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[2]))*(**jit)[8] - ((p[0]+p[3])*(4./(3.*m2)*p[2]*p[3]-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[2]) + (p[1]-M_I*p[2])*(-M_I/3.-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[3]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[2]))*(**jit)[12] + ((p[0]+p[3])*(-M_I/3.-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[0]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[2]) + (p[1]-M_I*p[2])*(4./(3.*m2)*p[2]*p[0]-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[0]-1./(3.*m2)*(p[0]+p[3])*p[2]))*(**jit)[1] - ((p[0]+p[3])*(-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[1]-1./(3.*m2)*(-(p[0]-p[3]))*p[2]) + (p[1]-M_I*p[2])*(4./(3.*m2)*p[2]*p[1]-M_I/3.-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[1]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[2]))*(**jit)[5] - ((p[0]+p[3])*(-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[2]) + (p[1]-M_I*p[2])*(1.0+4./(3.*m2)*p[2]*p[2]-1./3.-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[2]-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[2]))*(**jit)[9] - ((p[0]+p[3])*(-M_I/3.-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[3]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[2]) + (p[1]-M_I*p[2])*(4./(3.*m2)*p[2]*p[3]-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[3]-1./(3.*m2)*(p[0]+p[3])*p[2]))*(**jit)[13] + ((r*m)*(4./(3.*m2)*p[2]*p[0]-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[0]-1./(3.*m2)*(p[0]+p[3])*p[2]))*(**jit)[2] - ((r*m)*(4./(3.*m2)*p[2]*p[1]+M_I/3.-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[1]-1./(3.*m2)*(p[1]-M_I*p[2])*p[2]))*(**jit)[6] - ((r*m)*(1.0+4./(3.*m2)*p[2]*p[2]-1./3.-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[2]-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[2]))*(**jit)[10] - ((r*m)*(4./(3.*m2)*p[2]*p[3]-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[3]-1./(3.*m2)*(p[0]+p[3])*p[2]))*(**jit)[14] + ((r*m)*(M_I/3.-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[0]-1./(3.*m2)*(p[1]-M_I*p[2])*p[2]))*(**jit)[3] - ((r*m)*(-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[1]-1./(3.*m2)*(p[0]+p[3])*p[2]))*(**jit)[7] - ((r*m)*(-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[2]))*(**jit)[11] - ((r*m)*(-M_I/3.-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[3]-1./(3.*m2)*(-(p[1]-M_I*p[2]))*p[2]))*(**jit)[15];
        j[11] = ((p[1]+M_I*p[2])*(4./(3.*m2)*p[2]*p[0]-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[0]-1./(3.*m2)*(p[0]-p[3])*p[2]) + (p[0]-p[3])*(M_I/3.-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[0]-1./(3.*m2)*(-p[1]-M_I*p[2])*p[2]))*(**jit)[0] - ((p[1]+M_I*p[2])*(4./(3.*m2)*p[2]*p[1]+M_I/3.-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[1]-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[2]) + (p[0]-p[3])*(-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[1]-1./(3.*m2)*(-(p[0]+p[3]))*p[2]))*(**jit)[4] - ((p[1]+M_I*p[2])*(1.0+4./(3.*m2)*p[2]*p[2]-1./3.-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[2]-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[2]) + (p[0]-p[3])*(-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[2]))*(**jit)[8] - ((p[1]+M_I*p[2])*(4./(3.*m2)*p[2]*p[3]-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[2]) + (p[0]-p[3])*(-M_I/3.-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[3]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[2]))*(**jit)[12] + ((p[1]+M_I*p[2])*(-M_I/3.-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[0]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[2]) + (p[0]-p[3])*(4./(3.*m2)*p[2]*p[0]-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[0]-1./(3.*m2)*(p[0]+p[3])*p[2]))*(**jit)[1] - ((p[1]+M_I*p[2])*(-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[1]-1./(3.*m2)*(-(p[0]-p[3]))*p[2]) + (p[0]-p[3])*(4./(3.*m2)*p[2]*p[1]-M_I/3.-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[1]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[2]))*(**jit)[5] - ((p[1]+M_I*p[2])*(-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[2]) + (p[0]-p[3])*(1.0+4./(3.*m2)*p[2]*p[2]-1./3.-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[2]-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[2]))*(**jit)[9] - ((p[1]+M_I*p[2])*(-M_I/3.-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[3]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[2]) + (p[0]-p[3])*(4./(3.*m2)*p[2]*p[3]-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[3]-1./(3.*m2)*(p[0]+p[3])*p[2]))*(**jit)[13] + ((r*m)*(-M_I/3.-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[0]-1./(3.*m2)*(p[1]+M_I*p[2])*p[2]))*(**jit)[2] - ((r*m)*(-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[1]-1./(3.*m2)*(p[0]-p[3])*p[2]))*(**jit)[6] - ((r*m)*(-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[2]))*(**jit)[10] - ((r*m)*(-M_I/3.-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[3]-1./(3.*m2)*(p[1]+M_I*p[2])*p[2]))*(**jit)[14] + ((r*m)*(4./(3.*m2)*p[2]*p[0]-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[0]-1./(3.*m2)*(p[0]-p[3])*p[2]))*(**jit)[3] - ((r*m)*(4./(3.*m2)*p[2]*p[1]-M_I/3.-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[1]-1./(3.*m2)*(p[1]+M_I*p[2])*p[2]))*(**jit)[7] - ((r*m)*(1.0+4./(3.*m2)*p[2]*p[2]-1./3.-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[2]-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[2]))*(**jit)[11] - ((r*m)*(4./(3.*m2)*p[2]*p[3]-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[2]))*(**jit)[15];
        j[12] = ((r*m)*(4./(3.*m2)*p[3]*p[0]+1./3.-1./(3.*m2)*(p[0]+p[3])*p[0]-1./(3.*m2)*(p[0]-p[3])*p[3]))*(**jit)[0] - ((r*m)*(4./(3.*m2)*p[3]*p[1]-1./(3.*m2)*(p[0]+p[3])*p[1]-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[3]))*(**jit)[4] - ((r*m)*(4./(3.*m2)*p[3]*p[2]-1./(3.*m2)*(p[0]+p[3])*p[2]-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[3]))*(**jit)[8] - ((r*m)*(1.0+4./(3.*m2)*p[3]*p[3]-1./3.-1./(3.*m2)*(p[0]+p[3])*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[3]))*(**jit)[12] + ((r*m)*(-1./(3.*m2)*(p[1]-M_I*p[2])*p[0]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[3]))*(**jit)[1] - ((r*m)*(-1./3.-1./(3.*m2)*(p[1]-M_I*p[2])*p[1]-1./(3.*m2)*(-(p[0]-p[3]))*p[3]))*(**jit)[5] - ((r*m)*(M_I/3.-1./(3.*m2)*(p[1]-M_I*p[2])*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[3]))*(**jit)[9] - ((r*m)*(-1./(3.*m2)*(p[1]-M_I*p[2])*p[3]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[3]))*(**jit)[13] + ((p[0]-p[3])*(4./(3.*m2)*p[3]*p[0]-1./3.-1./(3.*m2)*(-(p[0]-p[3]))*p[0]-1./(3.*m2)*(p[0]+p[3])*p[3]) + (-p[1]+M_I*p[2])*(-1./(3.*m2)*(-p[1]-M_I*p[2])*p[0]-1./(3.*m2)*(p[1]+M_I*p[2])*p[3]))*(**jit)[2] - ((p[0]-p[3])*(4./(3.*m2)*p[3]*p[1]-1./(3.*m2)*(-(p[0]-p[3]))*p[1]-1./(3.*m2)*(p[1]-M_I*p[2])*p[3]) + (-p[1]+M_I*p[2])*(1./3.-1./(3.*m2)*(-p[1]-M_I*p[2])*p[1]-1./(3.*m2)*(p[0]-p[3])*p[3]))*(**jit)[6] - ((p[0]-p[3])*(4./(3.*m2)*p[3]*p[2]-1./(3.*m2)*(-(p[0]-p[3]))*p[2]-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[3]) + (-p[1]+M_I*p[2])*(M_I/3.-1./(3.*m2)*(-p[1]-M_I*p[2])*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[3]))*(**jit)[10] - ((p[0]-p[3])*(1.0+4./(3.*m2)*p[3]*p[3]-1./3.-1./(3.*m2)*(-(p[0]-p[3]))*p[3]-1./(3.*m2)*(p[0]+p[3])*p[3]) + (-p[1]+M_I*p[2])*(-1./(3.*m2)*(-p[1]-M_I*p[2])*p[3]-1./(3.*m2)*(p[1]+M_I*p[2])*p[3]))*(**jit)[14] + ((p[0]-p[3])*(-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[0]-1./(3.*m2)*(p[1]-M_I*p[2])*p[3]) + (-p[1]+M_I*p[2])*(4./(3.*m2)*p[3]*p[0]+1./3.-1./(3.*m2)*(p[0]+p[3])*p[0]-1./(3.*m2)*(p[0]-p[3])*p[3]))*(**jit)[3] - ((p[0]-p[3])*(-1./3.-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[1]-1./(3.*m2)*(p[0]+p[3])*p[3]) + (-p[1]+M_I*p[2])*(4./(3.*m2)*p[3]*p[1]-1./(3.*m2)*(p[0]+p[3])*p[1]-1./(3.*m2)*(p[1]+M_I*p[2])*p[3]))*(**jit)[7] - ((p[0]-p[3])*(M_I/3.-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[3]) + (-p[1]+M_I*p[2])*(4./(3.*m2)*p[3]*p[2]-1./(3.*m2)*(p[0]+p[3])*p[2]-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[3]))*(**jit)[11] - ((p[0]-p[3])*(-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[3]-1./(3.*m2)*(-(p[1]-M_I*p[2]))*p[3]) + (-p[1]+M_I*p[2])*(1.0+4./(3.*m2)*p[3]*p[3]-1./3.-1./(3.*m2)*(p[0]+p[3])*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[3]))*(**jit)[15];
        j[13] = ((r*m)*(-1./(3.*m2)*(-(p[1]+M_I*p[2]))*p[0]-1./(3.*m2)*(-p[1]-M_I*p[2])*p[3]))*(**jit)[0] - ((r*m)*(1./3.-1./(3.*m2)*(-(p[1]+M_I*p[2]))*p[1]-1./(3.*m2)*(-(p[0]+p[3]))*p[3]))*(**jit)[4] - ((r*m)*(M_I/3.-1./(3.*m2)*(-(p[1]+M_I*p[2]))*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[3]))*(**jit)[8] - ((r*m)*(-1./(3.*m2)*(-(p[1]+M_I*p[2]))*p[3]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[3]))*(**jit)[12] + ((r*m)*(4./(3.*m2)*p[3]*p[0]-1./3.-1./(3.*m2)*(-(p[0]-p[3]))*p[0]-1./(3.*m2)*(p[0]+p[3])*p[3]))*(**jit)[1] - ((r*m)*(4./(3.*m2)*p[3]*p[1]-1./(3.*m2)*(-(p[0]-p[3]))*p[1]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[3]))*(**jit)[5] - ((r*m)*(4./(3.*m2)*p[3]*p[2]-1./(3.*m2)*(-(p[0]-p[3]))*p[2]-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[3]))*(**jit)[9] - ((r*m)*(1.0+4./(3.*m2)*p[3]*p[3]-1./3.-1./(3.*m2)*(-(p[0]-p[3]))*p[3]-1./(3.*m2)*(p[0]+p[3])*p[3]))*(**jit)[13] + ((-p[1]-M_I*p[2])*(4./(3.*m2)*p[3]*p[0]-1./3.-1./(3.*m2)*(-(p[0]-p[3]))*p[0]-1./(3.*m2)*(p[0]+p[3])*p[3]) + (p[0]+p[3])*(-1./(3.*m2)*(-p[1]-M_I*p[2])*p[0]-1./(3.*m2)*(p[1]+M_I*p[2])*p[3]))*(**jit)[2] - ((-p[1]-M_I*p[2])*(4./(3.*m2)*p[3]*p[1]-1./(3.*m2)*(-(p[0]-p[3]))*p[1]-1./(3.*m2)*(p[1]-M_I*p[2])*p[3]) + (p[0]+p[3])*(1./3.-1./(3.*m2)*(-p[1]-M_I*p[2])*p[1]-1./(3.*m2)*(p[0]-p[3])*p[3]))*(**jit)[6] - ((-p[1]-M_I*p[2])*(4./(3.*m2)*p[3]*p[2]-1./(3.*m2)*(-(p[0]-p[3]))*p[2]-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[3]) + (p[0]+p[3])*(M_I/3.-1./(3.*m2)*(-p[1]-M_I*p[2])*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[3]))*(**jit)[10] - ((-p[1]-M_I*p[2])*(1.0+4./(3.*m2)*p[3]*p[3]-1./3.-1./(3.*m2)*(-(p[0]-p[3]))*p[3]-1./(3.*m2)*(p[0]+p[3])*p[3]) + (p[0]+p[3])*(-1./(3.*m2)*(-p[1]-M_I*p[2])*p[3]-1./(3.*m2)*(p[1]+M_I*p[2])*p[3]))*(**jit)[14] + ((-p[1]-M_I*p[2])*(-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[0]-1./(3.*m2)*(p[1]-M_I*p[2])*p[3]) + (p[0]+p[3])*(4./(3.*m2)*p[3]*p[0]+1./3.-1./(3.*m2)*(p[0]+p[3])*p[0]-1./(3.*m2)*(p[0]-p[3])*p[3]))*(**jit)[3] - ((-p[1]-M_I*p[2])*(-1./3.-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[1]-1./(3.*m2)*(p[0]+p[3])*p[3]) + (p[0]+p[3])*(4./(3.*m2)*p[3]*p[1]-1./(3.*m2)*(p[0]+p[3])*p[1]-1./(3.*m2)*(p[1]+M_I*p[2])*p[3]))*(**jit)[7] - ((-p[1]-M_I*p[2])*(M_I/3.-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[3]) + (p[0]+p[3])*(4./(3.*m2)*p[3]*p[2]-1./(3.*m2)*(p[0]+p[3])*p[2]-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[3]))*(**jit)[11] - ((-p[1]-M_I*p[2])*(-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[3]-1./(3.*m2)*(-(p[1]-M_I*p[2]))*p[3]) + (p[0]+p[3])*(1.0+4./(3.*m2)*p[3]*p[3]-1./3.-1./(3.*m2)*(p[0]+p[3])*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[3]))*(**jit)[15];
        j[14] = ((p[0]+p[3])*(4./(3.*m2)*p[3]*p[0]+1./3.-1./(3.*m2)*(p[0]+p[3])*p[0]-1./(3.*m2)*(p[0]-p[3])*p[3]) + (p[1]-M_I*p[2])*(-1./(3.*m2)*(-(p[1]+M_I*p[2]))*p[0]-1./(3.*m2)*(-p[1]-M_I*p[2])*p[3]))*(**jit)[0] - ((p[0]+p[3])*(4./(3.*m2)*p[3]*p[1]-1./(3.*m2)*(p[0]+p[3])*p[1]-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[3]) + (p[1]-M_I*p[2])*(1./3.-1./(3.*m2)*(-(p[1]+M_I*p[2]))*p[1]-1./(3.*m2)*(-(p[0]+p[3]))*p[3]))*(**jit)[4] - ((p[0]+p[3])*(4./(3.*m2)*p[3]*p[2]-1./(3.*m2)*(p[0]+p[3])*p[2]-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[3]) + (p[1]-M_I*p[2])*(M_I/3.-1./(3.*m2)*(-(p[1]+M_I*p[2]))*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[3]))*(**jit)[8] - ((p[0]+p[3])*(1.0+4./(3.*m2)*p[3]*p[3]-1./3.-1./(3.*m2)*(p[0]+p[3])*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[3]) + (p[1]-M_I*p[2])*(-1./(3.*m2)*(-(p[1]+M_I*p[2]))*p[3]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[3]))*(**jit)[12] + ((p[0]+p[3])*(-1./(3.*m2)*(p[1]-M_I*p[2])*p[0]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[3]) + (p[1]-M_I*p[2])*(4./(3.*m2)*p[3]*p[0]-1./3.-1./(3.*m2)*(-(p[0]-p[3]))*p[0]-1./(3.*m2)*(p[0]+p[3])*p[3]))*(**jit)[1] - ((p[0]+p[3])*(-1./3.-1./(3.*m2)*(p[1]-M_I*p[2])*p[1]-1./(3.*m2)*(-(p[0]-p[3]))*p[3]) + (p[1]-M_I*p[2])*(4./(3.*m2)*p[3]*p[1]-1./(3.*m2)*(-(p[0]-p[3]))*p[1]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[3]))*(**jit)[5] - ((p[0]+p[3])*(M_I/3.-1./(3.*m2)*(p[1]-M_I*p[2])*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[3]) + (p[1]-M_I*p[2])*(4./(3.*m2)*p[3]*p[2]-1./(3.*m2)*(-(p[0]-p[3]))*p[2]-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[3]))*(**jit)[9] - ((p[0]+p[3])*(-1./(3.*m2)*(p[1]-M_I*p[2])*p[3]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[3]) + (p[1]-M_I*p[2])*(1.0+4./(3.*m2)*p[3]*p[3]-1./3.-1./(3.*m2)*(-(p[0]-p[3]))*p[3]-1./(3.*m2)*(p[0]+p[3])*p[3]))*(**jit)[13] + ((r*m)*(4./(3.*m2)*p[3]*p[0]-1./3.-1./(3.*m2)*(-(p[0]-p[3]))*p[0]-1./(3.*m2)*(p[0]+p[3])*p[3]))*(**jit)[2] - ((r*m)*(4./(3.*m2)*p[3]*p[1]-1./(3.*m2)*(-(p[0]-p[3]))*p[1]-1./(3.*m2)*(p[1]-M_I*p[2])*p[3]))*(**jit)[6] - ((r*m)*(4./(3.*m2)*p[3]*p[2]-1./(3.*m2)*(-(p[0]-p[3]))*p[2]-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[3]))*(**jit)[10] - ((r*m)*(1.0+4./(3.*m2)*p[3]*p[3]-1./3.-1./(3.*m2)*(-(p[0]-p[3]))*p[3]-1./(3.*m2)*(p[0]+p[3])*p[3]))*(**jit)[14] + ((r*m)*(-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[0]-1./(3.*m2)*(p[1]-M_I*p[2])*p[3]))*(**jit)[3] - ((r*m)*(-1./3.-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[1]-1./(3.*m2)*(p[0]+p[3])*p[3]))*(**jit)[7] - ((r*m)*(M_I/3.-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[3]))*(**jit)[11] - ((r*m)*(-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[3]-1./(3.*m2)*(-(p[1]-M_I*p[2]))*p[3]))*(**jit)[15];
        j[15] = ((p[1]+M_I*p[2])*(4./(3.*m2)*p[3]*p[0]+1./3.-1./(3.*m2)*(p[0]+p[3])*p[0]-1./(3.*m2)*(p[0]-p[3])*p[3]) + (p[0]-p[3])*(-1./(3.*m2)*(-(p[1]+M_I*p[2]))*p[0]-1./(3.*m2)*(-p[1]-M_I*p[2])*p[3]))*(**jit)[0] - ((p[1]+M_I*p[2])*(4./(3.*m2)*p[3]*p[1]-1./(3.*m2)*(p[0]+p[3])*p[1]-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[3]) + (p[0]-p[3])*(1./3.-1./(3.*m2)*(-(p[1]+M_I*p[2]))*p[1]-1./(3.*m2)*(-(p[0]+p[3]))*p[3]))*(**jit)[4] - ((p[1]+M_I*p[2])*(4./(3.*m2)*p[3]*p[2]-1./(3.*m2)*(p[0]+p[3])*p[2]-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[3]) + (p[0]-p[3])*(M_I/3.-1./(3.*m2)*(-(p[1]+M_I*p[2]))*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[3]))*(**jit)[8] - ((p[1]+M_I*p[2])*(1.0+4./(3.*m2)*p[3]*p[3]-1./3.-1./(3.*m2)*(p[0]+p[3])*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[3]) + (p[0]-p[3])*(-1./(3.*m2)*(-(p[1]+M_I*p[2]))*p[3]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[3]))*(**jit)[12] + ((p[1]+M_I*p[2])*(-1./(3.*m2)*(p[1]-M_I*p[2])*p[0]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[3]) + (p[0]-p[3])*(4./(3.*m2)*p[3]*p[0]-1./3.-1./(3.*m2)*(-(p[0]-p[3]))*p[0]-1./(3.*m2)*(p[0]+p[3])*p[3]))*(**jit)[1] - ((p[1]+M_I*p[2])*(-1./3.-1./(3.*m2)*(p[1]-M_I*p[2])*p[1]-1./(3.*m2)*(-(p[0]-p[3]))*p[3]) + (p[0]-p[3])*(4./(3.*m2)*p[3]*p[1]-1./(3.*m2)*(-(p[0]-p[3]))*p[1]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[3]))*(**jit)[5] - ((p[1]+M_I*p[2])*(M_I/3.-1./(3.*m2)*(p[1]-M_I*p[2])*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[3]) + (p[0]-p[3])*(4./(3.*m2)*p[3]*p[2]-1./(3.*m2)*(-(p[0]-p[3]))*p[2]-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[3]))*(**jit)[9] - ((p[1]+M_I*p[2])*(-1./(3.*m2)*(p[1]-M_I*p[2])*p[3]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[3]) + (p[0]-p[3])*(1.0+4./(3.*m2)*p[3]*p[3]-1./3.-1./(3.*m2)*(-(p[0]-p[3]))*p[3]-1./(3.*m2)*(p[0]+p[3])*p[3]))*(**jit)[13] + ((r*m)*(-1./(3.*m2)*(-p[1]-M_I*p[2])*p[0]-1./(3.*m2)*(p[1]+M_I*p[2])*p[3]))*(**jit)[2] - ((r*m)*(1./3.-1./(3.*m2)*(-p[1]-M_I*p[2])*p[1]-1./(3.*m2)*(p[0]-p[3])*p[3]))*(**jit)[6] - ((r*m)*(M_I/3.-1./(3.*m2)*(-p[1]-M_I*p[2])*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[3]))*(**jit)[10] - ((r*m)*(-1./(3.*m2)*(-p[1]-M_I*p[2])*p[3]-1./(3.*m2)*(p[1]+M_I*p[2])*p[3]))*(**jit)[14] + ((r*m)*(4./(3.*m2)*p[3]*p[0]+1./3.-1./(3.*m2)*(p[0]+p[3])*p[0]-1./(3.*m2)*(p[0]-p[3])*p[3]))*(**jit)[3] - ((r*m)*(4./(3.*m2)*p[3]*p[1]-1./(3.*m2)*(p[0]+p[3])*p[1]-1./(3.*m2)*(p[1]+M_I*p[2])*p[3]))*(**jit)[7] - ((r*m)*(4./(3.*m2)*p[3]*p[2]-1./(3.*m2)*(p[0]+p[3])*p[2]-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[3]))*(**jit)[11] - ((r*m)*(1.0+4./(3.*m2)*p[3]*p[3]-1./3.-1./(3.*m2)*(p[0]+p[3])*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[3]))*(**jit)[15];
      }
      else {// S(p)
        j[0] = (**jit)[0]*((r*m)*(-1.0+4./(3.*m2)*p[0]*p[0]+1./3.-1./(3.*m2)*(p[0]+p[3])*p[0]-1./(3.*m2)*(p[0]-p[3])*p[0])) - (**jit)[4]*((r*m)*(4./(3.*m2)*p[1]*p[0]-1./(3.*m2)*(p[1]+M_I*p[2])*p[0]-1./(3.*m2)*(p[0]-p[3])*p[1])) - (**jit)[8]*((r*m)*(4./(3.*m2)*p[2]*p[0]-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[0]-1./(3.*m2)*(p[0]-p[3])*p[2])) - (**jit)[12]*((r*m)*(4./(3.*m2)*p[3]*p[0]+1./3.-1./(3.*m2)*(p[0]+p[3])*p[0]-1./(3.*m2)*(p[0]-p[3])*p[3])) + (**jit)[1]*((r*m)*(-1./(3.*m2)*(p[1]+M_I*p[2])*p[0]-1./(3.*m2)*(-p[1]-M_I*p[2])*p[0])) - (**jit)[5]*((r*m)*(1./3.-1./(3.*m2)*(p[0]+p[3])*p[0]-1./(3.*m2)*(-p[1]-M_I*p[2])*p[1])) - (**jit)[9]*((r*m)*(M_I/3.-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[0]-1./(3.*m2)*(-p[1]-M_I*p[2])*p[2])) - (**jit)[13]*((r*m)*(-1./(3.*m2)*(-(p[1]+M_I*p[2]))*p[0]-1./(3.*m2)*(-p[1]-M_I*p[2])*p[3])) + (**jit)[2]*((p[0]+p[3])*(-1.0+4./(3.*m2)*p[0]*p[0]+1./3.-1./(3.*m2)*(p[0]+p[3])*p[0]-1./(3.*m2)*(p[0]-p[3])*p[0]) + (p[1]-M_I*p[2])*(-1./(3.*m2)*(p[1]+M_I*p[2])*p[0]-1./(3.*m2)*(-p[1]-M_I*p[2])*p[0])) - (**jit)[6]*((p[0]+p[3])*(4./(3.*m2)*p[1]*p[0]-1./(3.*m2)*(p[1]+M_I*p[2])*p[0]-1./(3.*m2)*(p[0]-p[3])*p[1]) + (p[1]-M_I*p[2])*(1./3.-1./(3.*m2)*(p[0]+p[3])*p[0]-1./(3.*m2)*(-p[1]-M_I*p[2])*p[1])) - (**jit)[10]*((p[0]+p[3])*(4./(3.*m2)*p[2]*p[0]-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[0]-1./(3.*m2)*(p[0]-p[3])*p[2]) + (p[1]-M_I*p[2])*(M_I/3.-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[0]-1./(3.*m2)*(-p[1]-M_I*p[2])*p[2])) - (**jit)[14]*((p[0]+p[3])*(4./(3.*m2)*p[3]*p[0]+1./3.-1./(3.*m2)*(p[0]+p[3])*p[0]-1./(3.*m2)*(p[0]-p[3])*p[3]) + (p[1]-M_I*p[2])*(-1./(3.*m2)*(-(p[1]+M_I*p[2]))*p[0]-1./(3.*m2)*(-p[1]-M_I*p[2])*p[3])) + (**jit)[3]*((p[1]+M_I*p[2])*(-1.0+4./(3.*m2)*p[0]*p[0]+1./3.-1./(3.*m2)*(p[0]+p[3])*p[0]-1./(3.*m2)*(p[0]-p[3])*p[0]) + (p[0]-p[3])*(-1./(3.*m2)*(p[1]+M_I*p[2])*p[0]-1./(3.*m2)*(-p[1]-M_I*p[2])*p[0])) - (**jit)[7]*((p[1]+M_I*p[2])*(4./(3.*m2)*p[1]*p[0]-1./(3.*m2)*(p[1]+M_I*p[2])*p[0]-1./(3.*m2)*(p[0]-p[3])*p[1]) + (p[0]-p[3])*(1./3.-1./(3.*m2)*(p[0]+p[3])*p[0]-1./(3.*m2)*(-p[1]-M_I*p[2])*p[1])) - (**jit)[11]*((p[1]+M_I*p[2])*(4./(3.*m2)*p[2]*p[0]-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[0]-1./(3.*m2)*(p[0]-p[3])*p[2]) + (p[0]-p[3])*(M_I/3.-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[0]-1./(3.*m2)*(-p[1]-M_I*p[2])*p[2])) - (**jit)[15]*((p[1]+M_I*p[2])*(4./(3.*m2)*p[3]*p[0]+1./3.-1./(3.*m2)*(p[0]+p[3])*p[0]-1./(3.*m2)*(p[0]-p[3])*p[3]) + (p[0]-p[3])*(-1./(3.*m2)*(-(p[1]+M_I*p[2]))*p[0]-1./(3.*m2)*(-p[1]-M_I*p[2])*p[3]));
        j[1] = (**jit)[0]*((r*m)*(-1./(3.*m2)*(p[1]-M_I*p[2])*p[0]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[0])) - (**jit)[4]*((r*m)*(1./3.-1./(3.*m2)*(p[0]-p[3])*p[0]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[1])) - (**jit)[8]*((r*m)*(-M_I/3.-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[0]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[2])) - (**jit)[12]*((r*m)*(-1./(3.*m2)*(p[1]-M_I*p[2])*p[0]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[3])) + (**jit)[1]*((r*m)*(-1.0+4./(3.*m2)*p[0]*p[0]+1./3.-1./(3.*m2)*(p[0]-p[3])*p[0]-1./(3.*m2)*(p[0]+p[3])*p[0])) - (**jit)[5]*((r*m)*(4./(3.*m2)*p[1]*p[0]-1./(3.*m2)*(p[1]-M_I*p[2])*p[0]-1./(3.*m2)*(p[0]+p[3])*p[1])) - (**jit)[9]*((r*m)*(4./(3.*m2)*p[2]*p[0]-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[0]-1./(3.*m2)*(p[0]+p[3])*p[2])) - (**jit)[13]*((r*m)*(4./(3.*m2)*p[3]*p[0]-1./3.-1./(3.*m2)*(-(p[0]-p[3]))*p[0]-1./(3.*m2)*(p[0]+p[3])*p[3])) + (**jit)[2]*((p[0]+p[3])*(-1./(3.*m2)*(p[1]-M_I*p[2])*p[0]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[0]) + (p[1]-M_I*p[2])*(-1.0+4./(3.*m2)*p[0]*p[0]+1./3.-1./(3.*m2)*(p[0]-p[3])*p[0]-1./(3.*m2)*(p[0]+p[3])*p[0])) - (**jit)[6]*((p[0]+p[3])*(1./3.-1./(3.*m2)*(p[0]-p[3])*p[0]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[1]) + (p[1]-M_I*p[2])*(4./(3.*m2)*p[1]*p[0]-1./(3.*m2)*(p[1]-M_I*p[2])*p[0]-1./(3.*m2)*(p[0]+p[3])*p[1])) - (**jit)[10]*((p[0]+p[3])*(-M_I/3.-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[0]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[2]) + (p[1]-M_I*p[2])*(4./(3.*m2)*p[2]*p[0]-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[0]-1./(3.*m2)*(p[0]+p[3])*p[2])) - (**jit)[14]*((p[0]+p[3])*(-1./(3.*m2)*(p[1]-M_I*p[2])*p[0]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[3]) + (p[1]-M_I*p[2])*(4./(3.*m2)*p[3]*p[0]-1./3.-1./(3.*m2)*(-(p[0]-p[3]))*p[0]-1./(3.*m2)*(p[0]+p[3])*p[3])) + (**jit)[3]*((p[1]+M_I*p[2])*(-1./(3.*m2)*(p[1]-M_I*p[2])*p[0]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[0]) + (p[0]-p[3])*(-1.0+4./(3.*m2)*p[0]*p[0]+1./3.-1./(3.*m2)*(p[0]-p[3])*p[0]-1./(3.*m2)*(p[0]+p[3])*p[0])) - (**jit)[7]*((p[1]+M_I*p[2])*(1./3.-1./(3.*m2)*(p[0]-p[3])*p[0]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[1]) + (p[0]-p[3])*(4./(3.*m2)*p[1]*p[0]-1./(3.*m2)*(p[1]-M_I*p[2])*p[0]-1./(3.*m2)*(p[0]+p[3])*p[1])) - (**jit)[11]*((p[1]+M_I*p[2])*(-M_I/3.-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[0]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[2]) + (p[0]-p[3])*(4./(3.*m2)*p[2]*p[0]-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[0]-1./(3.*m2)*(p[0]+p[3])*p[2])) - (**jit)[15]*((p[1]+M_I*p[2])*(-1./(3.*m2)*(p[1]-M_I*p[2])*p[0]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[3]) + (p[0]-p[3])*(4./(3.*m2)*p[3]*p[0]-1./3.-1./(3.*m2)*(-(p[0]-p[3]))*p[0]-1./(3.*m2)*(p[0]+p[3])*p[3]));
        j[2] = (**jit)[0]*((p[0]-p[3])*(-1.0+4./(3.*m2)*p[0]*p[0]+1./3.-1./(3.*m2)*(p[0]-p[3])*p[0]-1./(3.*m2)*(p[0]+p[3])*p[0]) + (-p[1]+M_I*p[2])*(-1./(3.*m2)*(-p[1]-M_I*p[2])*p[0]-1./(3.*m2)*(p[1]+M_I*p[2])*p[0])) - (**jit)[4]*((p[0]-p[3])*(4./(3.*m2)*p[1]*p[0]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[0]-1./(3.*m2)*(p[0]+p[3])*p[1]) + (-p[1]+M_I*p[2])*(-1./3.-1./(3.*m2)*(-(p[0]-p[3]))*p[0]-1./(3.*m2)*(p[1]+M_I*p[2])*p[1])) - (**jit)[8]*((p[0]-p[3])*(4./(3.*m2)*p[2]*p[0]-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[0]-1./(3.*m2)*(p[0]+p[3])*p[2]) + (-p[1]+M_I*p[2])*(-M_I/3.-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[0]-1./(3.*m2)*(p[1]+M_I*p[2])*p[2])) - (**jit)[12]*((p[0]-p[3])*(4./(3.*m2)*p[3]*p[0]-1./3.-1./(3.*m2)*(-(p[0]-p[3]))*p[0]-1./(3.*m2)*(p[0]+p[3])*p[3]) + (-p[1]+M_I*p[2])*(-1./(3.*m2)*(-p[1]-M_I*p[2])*p[0]-1./(3.*m2)*(p[1]+M_I*p[2])*p[3])) + (**jit)[1]*((-p[1]-M_I*p[2])*(-1.0+4./(3.*m2)*p[0]*p[0]+1./3.-1./(3.*m2)*(p[0]-p[3])*p[0]-1./(3.*m2)*(p[0]+p[3])*p[0]) + (p[0]+p[3])*(-1./(3.*m2)*(-p[1]-M_I*p[2])*p[0]-1./(3.*m2)*(p[1]+M_I*p[2])*p[0])) - (**jit)[5]*((-p[1]-M_I*p[2])*(4./(3.*m2)*p[1]*p[0]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[0]-1./(3.*m2)*(p[0]+p[3])*p[1]) + (p[0]+p[3])*(-1./3.-1./(3.*m2)*(-(p[0]-p[3]))*p[0]-1./(3.*m2)*(p[1]+M_I*p[2])*p[1])) - (**jit)[9]*((-p[1]-M_I*p[2])*(4./(3.*m2)*p[2]*p[0]-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[0]-1./(3.*m2)*(p[0]+p[3])*p[2]) + (p[0]+p[3])*(-M_I/3.-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[0]-1./(3.*m2)*(p[1]+M_I*p[2])*p[2])) - (**jit)[13]*((-p[1]-M_I*p[2])*(4./(3.*m2)*p[3]*p[0]-1./3.-1./(3.*m2)*(-(p[0]-p[3]))*p[0]-1./(3.*m2)*(p[0]+p[3])*p[3]) + (p[0]+p[3])*(-1./(3.*m2)*(-p[1]-M_I*p[2])*p[0]-1./(3.*m2)*(p[1]+M_I*p[2])*p[3])) + (**jit)[2]*((r*m)*(-1.0+4./(3.*m2)*p[0]*p[0]+1./3.-1./(3.*m2)*(p[0]-p[3])*p[0]-1./(3.*m2)*(p[0]+p[3])*p[0])) - (**jit)[6]*((r*m)*(4./(3.*m2)*p[1]*p[0]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[0]-1./(3.*m2)*(p[0]+p[3])*p[1])) - (**jit)[10]*((r*m)*(4./(3.*m2)*p[2]*p[0]-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[0]-1./(3.*m2)*(p[0]+p[3])*p[2])) - (**jit)[14]*((r*m)*(4./(3.*m2)*p[3]*p[0]-1./3.-1./(3.*m2)*(-(p[0]-p[3]))*p[0]-1./(3.*m2)*(p[0]+p[3])*p[3])) + (**jit)[3]*((r*m)*(-1./(3.*m2)*(-p[1]-M_I*p[2])*p[0]-1./(3.*m2)*(p[1]+M_I*p[2])*p[0])) - (**jit)[7]*((r*m)*(-1./3.-1./(3.*m2)*(-(p[0]-p[3]))*p[0]-1./(3.*m2)*(p[1]+M_I*p[2])*p[1])) - (**jit)[11]*((r*m)*(-M_I/3.-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[0]-1./(3.*m2)*(p[1]+M_I*p[2])*p[2])) - (**jit)[15]*((r*m)*(-1./(3.*m2)*(-p[1]-M_I*p[2])*p[0]-1./(3.*m2)*(p[1]+M_I*p[2])*p[3]));
        j[3] = (**jit)[0]*((p[0]-p[3])*(-1./(3.*m2)*(-p[1]+M_I*p[2])*p[0]-1./(3.*m2)*(p[1]-M_I*p[2])*p[0]) + (-p[1]+M_I*p[2])*(-1.0+4./(3.*m2)*p[0]*p[0]+1./3.-1./(3.*m2)*(p[0]+p[3])*p[0]-1./(3.*m2)*(p[0]-p[3])*p[0])) - (**jit)[4]*((p[0]-p[3])*(-1./3.-1./(3.*m2)*(-(p[0]+p[3]))*p[0]-1./(3.*m2)*(p[1]-M_I*p[2])*p[1]) + (-p[1]+M_I*p[2])*(4./(3.*m2)*p[1]*p[0]-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[0]-1./(3.*m2)*(p[0]-p[3])*p[1])) - (**jit)[8]*((p[0]-p[3])*(M_I/3.-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[0]-1./(3.*m2)*(p[1]-M_I*p[2])*p[2]) + (-p[1]+M_I*p[2])*(4./(3.*m2)*p[2]*p[0]-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[0]-1./(3.*m2)*(p[0]-p[3])*p[2])) - (**jit)[12]*((p[0]-p[3])*(-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[0]-1./(3.*m2)*(p[1]-M_I*p[2])*p[3]) + (-p[1]+M_I*p[2])*(4./(3.*m2)*p[3]*p[0]+1./3.-1./(3.*m2)*(p[0]+p[3])*p[0]-1./(3.*m2)*(p[0]-p[3])*p[3])) + (**jit)[1]*((-p[1]-M_I*p[2])*(-1./(3.*m2)*(-p[1]+M_I*p[2])*p[0]-1./(3.*m2)*(p[1]-M_I*p[2])*p[0]) + (p[0]+p[3])*(-1.0+4./(3.*m2)*p[0]*p[0]+1./3.-1./(3.*m2)*(p[0]+p[3])*p[0]-1./(3.*m2)*(p[0]-p[3])*p[0])) - (**jit)[5]*((-p[1]-M_I*p[2])*(-1./3.-1./(3.*m2)*(-(p[0]+p[3]))*p[0]-1./(3.*m2)*(p[1]-M_I*p[2])*p[1]) + (p[0]+p[3])*(4./(3.*m2)*p[1]*p[0]-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[0]-1./(3.*m2)*(p[0]-p[3])*p[1])) - (**jit)[9]*((-p[1]-M_I*p[2])*(M_I/3.-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[0]-1./(3.*m2)*(p[1]-M_I*p[2])*p[2]) + (p[0]+p[3])*(4./(3.*m2)*p[2]*p[0]-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[0]-1./(3.*m2)*(p[0]-p[3])*p[2])) - (**jit)[13]*((-p[1]-M_I*p[2])*(-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[0]-1./(3.*m2)*(p[1]-M_I*p[2])*p[3]) + (p[0]+p[3])*(4./(3.*m2)*p[3]*p[0]+1./3.-1./(3.*m2)*(p[0]+p[3])*p[0]-1./(3.*m2)*(p[0]-p[3])*p[3])) + (**jit)[2]*((r*m)*(-1./(3.*m2)*(-p[1]+M_I*p[2])*p[0]-1./(3.*m2)*(p[1]-M_I*p[2])*p[0])) - (**jit)[6]*((r*m)*(-1./3.-1./(3.*m2)*(-(p[0]+p[3]))*p[0]-1./(3.*m2)*(p[1]-M_I*p[2])*p[1])) - (**jit)[10]*((r*m)*(M_I/3.-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[0]-1./(3.*m2)*(p[1]-M_I*p[2])*p[2])) - (**jit)[14]*((r*m)*(-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[0]-1./(3.*m2)*(p[1]-M_I*p[2])*p[3])) + (**jit)[3]*((r*m)*(-1.0+4./(3.*m2)*p[0]*p[0]+1./3.-1./(3.*m2)*(p[0]+p[3])*p[0]-1./(3.*m2)*(p[0]-p[3])*p[0])) - (**jit)[7]*((r*m)*(4./(3.*m2)*p[1]*p[0]-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[0]-1./(3.*m2)*(p[0]-p[3])*p[1])) - (**jit)[11]*((r*m)*(4./(3.*m2)*p[2]*p[0]-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[0]-1./(3.*m2)*(p[0]-p[3])*p[2])) - (**jit)[15]*((r*m)*(4./(3.*m2)*p[3]*p[0]+1./3.-1./(3.*m2)*(p[0]+p[3])*p[0]-1./(3.*m2)*(p[0]-p[3])*p[3]));
        j[4] = (**jit)[0]*((r*m)*(4./(3.*m2)*p[0]*p[1]-1./(3.*m2)*(p[0]+p[3])*p[1]-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[0])) - (**jit)[4]*((r*m)*(1.0+4./(3.*m2)*p[1]*p[1]-1./3.-1./(3.*m2)*(p[1]+M_I*p[2])*p[1]-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[1])) - (**jit)[8]*((r*m)*(4./(3.*m2)*p[2]*p[1]+M_I/3.-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[1]-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[2])) - (**jit)[12]*((r*m)*(4./(3.*m2)*p[3]*p[1]-1./(3.*m2)*(p[0]+p[3])*p[1]-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[3])) + (**jit)[1]*((r*m)*(-1./3.-1./(3.*m2)*(p[1]+M_I*p[2])*p[1]-1./(3.*m2)*(-(p[0]+p[3]))*p[0])) - (**jit)[5]*((r*m)*(-1./(3.*m2)*(p[0]+p[3])*p[1]-1./(3.*m2)*(-(p[0]+p[3]))*p[1])) - (**jit)[9]*((r*m)*(-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[1]-1./(3.*m2)*(-(p[0]+p[3]))*p[2])) - (**jit)[13]*((r*m)*(1./3.-1./(3.*m2)*(-(p[1]+M_I*p[2]))*p[1]-1./(3.*m2)*(-(p[0]+p[3]))*p[3])) + (**jit)[2]*((p[0]+p[3])*(4./(3.*m2)*p[0]*p[1]-1./(3.*m2)*(p[0]+p[3])*p[1]-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[0]) + (p[1]-M_I*p[2])*(-1./3.-1./(3.*m2)*(p[1]+M_I*p[2])*p[1]-1./(3.*m2)*(-(p[0]+p[3]))*p[0])) - (**jit)[6]*((p[0]+p[3])*(1.0+4./(3.*m2)*p[1]*p[1]-1./3.-1./(3.*m2)*(p[1]+M_I*p[2])*p[1]-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[1]) + (p[1]-M_I*p[2])*(-1./(3.*m2)*(p[0]+p[3])*p[1]-1./(3.*m2)*(-(p[0]+p[3]))*p[1])) - (**jit)[10]*((p[0]+p[3])*(4./(3.*m2)*p[2]*p[1]+M_I/3.-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[1]-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[2]) + (p[1]-M_I*p[2])*(-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[1]-1./(3.*m2)*(-(p[0]+p[3]))*p[2])) - (**jit)[14]*((p[0]+p[3])*(4./(3.*m2)*p[3]*p[1]-1./(3.*m2)*(p[0]+p[3])*p[1]-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[3]) + (p[1]-M_I*p[2])*(1./3.-1./(3.*m2)*(-(p[1]+M_I*p[2]))*p[1]-1./(3.*m2)*(-(p[0]+p[3]))*p[3])) + (**jit)[3]*((p[1]+M_I*p[2])*(4./(3.*m2)*p[0]*p[1]-1./(3.*m2)*(p[0]+p[3])*p[1]-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[0]) + (p[0]-p[3])*(-1./3.-1./(3.*m2)*(p[1]+M_I*p[2])*p[1]-1./(3.*m2)*(-(p[0]+p[3]))*p[0])) - (**jit)[7]*((p[1]+M_I*p[2])*(1.0+4./(3.*m2)*p[1]*p[1]-1./3.-1./(3.*m2)*(p[1]+M_I*p[2])*p[1]-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[1]) + (p[0]-p[3])*(-1./(3.*m2)*(p[0]+p[3])*p[1]-1./(3.*m2)*(-(p[0]+p[3]))*p[1])) - (**jit)[11]*((p[1]+M_I*p[2])*(4./(3.*m2)*p[2]*p[1]+M_I/3.-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[1]-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[2]) + (p[0]-p[3])*(-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[1]-1./(3.*m2)*(-(p[0]+p[3]))*p[2])) - (**jit)[15]*((p[1]+M_I*p[2])*(4./(3.*m2)*p[3]*p[1]-1./(3.*m2)*(p[0]+p[3])*p[1]-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[3]) + (p[0]-p[3])*(1./3.-1./(3.*m2)*(-(p[1]+M_I*p[2]))*p[1]-1./(3.*m2)*(-(p[0]+p[3]))*p[3]));
        j[5] = (**jit)[0]*((r*m)*(-1./3.-1./(3.*m2)*(p[1]-M_I*p[2])*p[1]-1./(3.*m2)*(-(p[0]-p[3]))*p[0])) - (**jit)[4]*((r*m)*(-1./(3.*m2)*(p[0]-p[3])*p[1]-1./(3.*m2)*(-(p[0]-p[3]))*p[1])) - (**jit)[8]*((r*m)*(-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[1]-1./(3.*m2)*(-(p[0]-p[3]))*p[2])) - (**jit)[12]*((r*m)*(-1./3.-1./(3.*m2)*(p[1]-M_I*p[2])*p[1]-1./(3.*m2)*(-(p[0]-p[3]))*p[3])) + (**jit)[1]*((r*m)*(4./(3.*m2)*p[0]*p[1]-1./(3.*m2)*(p[0]-p[3])*p[1]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[0])) - (**jit)[5]*((r*m)*(1.0+4./(3.*m2)*p[1]*p[1]-1./3.-1./(3.*m2)*(p[1]-M_I*p[2])*p[1]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[1])) - (**jit)[9]*((r*m)*(4./(3.*m2)*p[2]*p[1]-M_I/3.-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[1]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[2])) - (**jit)[13]*((r*m)*(4./(3.*m2)*p[3]*p[1]-1./(3.*m2)*(-(p[0]-p[3]))*p[1]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[3])) + (**jit)[2]*((p[0]+p[3])*(-1./3.-1./(3.*m2)*(p[1]-M_I*p[2])*p[1]-1./(3.*m2)*(-(p[0]-p[3]))*p[0]) + (p[1]-M_I*p[2])*(4./(3.*m2)*p[0]*p[1]-1./(3.*m2)*(p[0]-p[3])*p[1]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[0])) - (**jit)[6]*((p[0]+p[3])*(-1./(3.*m2)*(p[0]-p[3])*p[1]-1./(3.*m2)*(-(p[0]-p[3]))*p[1]) + (p[1]-M_I*p[2])*(1.0+4./(3.*m2)*p[1]*p[1]-1./3.-1./(3.*m2)*(p[1]-M_I*p[2])*p[1]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[1])) - (**jit)[10]*((p[0]+p[3])*(-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[1]-1./(3.*m2)*(-(p[0]-p[3]))*p[2]) + (p[1]-M_I*p[2])*(4./(3.*m2)*p[2]*p[1]-M_I/3.-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[1]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[2])) - (**jit)[14]*((p[0]+p[3])*(-1./3.-1./(3.*m2)*(p[1]-M_I*p[2])*p[1]-1./(3.*m2)*(-(p[0]-p[3]))*p[3]) + (p[1]-M_I*p[2])*(4./(3.*m2)*p[3]*p[1]-1./(3.*m2)*(-(p[0]-p[3]))*p[1]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[3])) + (**jit)[3]*((p[1]+M_I*p[2])*(-1./3.-1./(3.*m2)*(p[1]-M_I*p[2])*p[1]-1./(3.*m2)*(-(p[0]-p[3]))*p[0]) + (p[0]-p[3])*(4./(3.*m2)*p[0]*p[1]-1./(3.*m2)*(p[0]-p[3])*p[1]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[0])) - (**jit)[7]*((p[1]+M_I*p[2])*(-1./(3.*m2)*(p[0]-p[3])*p[1]-1./(3.*m2)*(-(p[0]-p[3]))*p[1]) + (p[0]-p[3])*(1.0+4./(3.*m2)*p[1]*p[1]-1./3.-1./(3.*m2)*(p[1]-M_I*p[2])*p[1]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[1])) - (**jit)[11]*((p[1]+M_I*p[2])*(-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[1]-1./(3.*m2)*(-(p[0]-p[3]))*p[2]) + (p[0]-p[3])*(4./(3.*m2)*p[2]*p[1]-M_I/3.-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[1]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[2])) - (**jit)[15]*((p[1]+M_I*p[2])*(-1./3.-1./(3.*m2)*(p[1]-M_I*p[2])*p[1]-1./(3.*m2)*(-(p[0]-p[3]))*p[3]) + (p[0]-p[3])*(4./(3.*m2)*p[3]*p[1]-1./(3.*m2)*(-(p[0]-p[3]))*p[1]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[3]));
        j[6] = (**jit)[0]*((p[0]-p[3])*(4./(3.*m2)*p[0]*p[1]-1./(3.*m2)*(p[0]-p[3])*p[1]-1./(3.*m2)*(p[1]-M_I*p[2])*p[0]) + (-p[1]+M_I*p[2])*(1./3.-1./(3.*m2)*(-p[1]-M_I*p[2])*p[1]-1./(3.*m2)*(p[0]-p[3])*p[0])) - (**jit)[4]*((p[0]-p[3])*(1.0+4./(3.*m2)*p[1]*p[1]-1./3.-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[1]-1./(3.*m2)*(p[1]-M_I*p[2])*p[1]) + (-p[1]+M_I*p[2])*(-1./(3.*m2)*(-(p[0]-p[3]))*p[1]-1./(3.*m2)*(p[0]-p[3])*p[1])) - (**jit)[8]*((p[0]-p[3])*(4./(3.*m2)*p[2]*p[1]+M_I/3.-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[1]-1./(3.*m2)*(p[1]-M_I*p[2])*p[2]) + (-p[1]+M_I*p[2])*(-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[1]-1./(3.*m2)*(p[0]-p[3])*p[2])) - (**jit)[12]*((p[0]-p[3])*(4./(3.*m2)*p[3]*p[1]-1./(3.*m2)*(-(p[0]-p[3]))*p[1]-1./(3.*m2)*(p[1]-M_I*p[2])*p[3]) + (-p[1]+M_I*p[2])*(1./3.-1./(3.*m2)*(-p[1]-M_I*p[2])*p[1]-1./(3.*m2)*(p[0]-p[3])*p[3])) + (**jit)[1]*((-p[1]-M_I*p[2])*(4./(3.*m2)*p[0]*p[1]-1./(3.*m2)*(p[0]-p[3])*p[1]-1./(3.*m2)*(p[1]-M_I*p[2])*p[0]) + (p[0]+p[3])*(1./3.-1./(3.*m2)*(-p[1]-M_I*p[2])*p[1]-1./(3.*m2)*(p[0]-p[3])*p[0])) - (**jit)[5]*((-p[1]-M_I*p[2])*(1.0+4./(3.*m2)*p[1]*p[1]-1./3.-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[1]-1./(3.*m2)*(p[1]-M_I*p[2])*p[1]) + (p[0]+p[3])*(-1./(3.*m2)*(-(p[0]-p[3]))*p[1]-1./(3.*m2)*(p[0]-p[3])*p[1])) - (**jit)[9]*((-p[1]-M_I*p[2])*(4./(3.*m2)*p[2]*p[1]+M_I/3.-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[1]-1./(3.*m2)*(p[1]-M_I*p[2])*p[2]) + (p[0]+p[3])*(-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[1]-1./(3.*m2)*(p[0]-p[3])*p[2])) - (**jit)[13]*((-p[1]-M_I*p[2])*(4./(3.*m2)*p[3]*p[1]-1./(3.*m2)*(-(p[0]-p[3]))*p[1]-1./(3.*m2)*(p[1]-M_I*p[2])*p[3]) + (p[0]+p[3])*(1./3.-1./(3.*m2)*(-p[1]-M_I*p[2])*p[1]-1./(3.*m2)*(p[0]-p[3])*p[3])) + (**jit)[2]*((r*m)*(4./(3.*m2)*p[0]*p[1]-1./(3.*m2)*(p[0]-p[3])*p[1]-1./(3.*m2)*(p[1]-M_I*p[2])*p[0])) - (**jit)[6]*((r*m)*(1.0+4./(3.*m2)*p[1]*p[1]-1./3.-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[1]-1./(3.*m2)*(p[1]-M_I*p[2])*p[1])) - (**jit)[10]*((r*m)*(4./(3.*m2)*p[2]*p[1]+M_I/3.-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[1]-1./(3.*m2)*(p[1]-M_I*p[2])*p[2])) - (**jit)[14]*((r*m)*(4./(3.*m2)*p[3]*p[1]-1./(3.*m2)*(-(p[0]-p[3]))*p[1]-1./(3.*m2)*(p[1]-M_I*p[2])*p[3])) + (**jit)[3]*((r*m)*(1./3.-1./(3.*m2)*(-p[1]-M_I*p[2])*p[1]-1./(3.*m2)*(p[0]-p[3])*p[0])) - (**jit)[7]*((r*m)*(-1./(3.*m2)*(-(p[0]-p[3]))*p[1]-1./(3.*m2)*(p[0]-p[3])*p[1])) - (**jit)[11]*((r*m)*(-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[1]-1./(3.*m2)*(p[0]-p[3])*p[2])) - (**jit)[15]*((r*m)*(1./3.-1./(3.*m2)*(-p[1]-M_I*p[2])*p[1]-1./(3.*m2)*(p[0]-p[3])*p[3]));
        j[7] = (**jit)[0]*((p[0]-p[3])*(1./3.-1./(3.*m2)*(-p[1]+M_I*p[2])*p[1]-1./(3.*m2)*(p[0]+p[3])*p[0]) + (-p[1]+M_I*p[2])*(4./(3.*m2)*p[0]*p[1]-1./(3.*m2)*(p[0]+p[3])*p[1]-1./(3.*m2)*(p[1]+M_I*p[2])*p[0])) - (**jit)[4]*((p[0]-p[3])*(-1./(3.*m2)*(-(p[0]+p[3]))*p[1]-1./(3.*m2)*(p[0]+p[3])*p[1]) + (-p[1]+M_I*p[2])*(1.0+4./(3.*m2)*p[1]*p[1]-1./3.-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[1]-1./(3.*m2)*(p[1]+M_I*p[2])*p[1])) - (**jit)[8]*((p[0]-p[3])*(-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[1]-1./(3.*m2)*(p[0]+p[3])*p[2]) + (-p[1]+M_I*p[2])*(4./(3.*m2)*p[2]*p[1]-M_I/3.-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[1]-1./(3.*m2)*(p[1]+M_I*p[2])*p[2])) - (**jit)[12]*((p[0]-p[3])*(-1./3.-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[1]-1./(3.*m2)*(p[0]+p[3])*p[3]) + (-p[1]+M_I*p[2])*(4./(3.*m2)*p[3]*p[1]-1./(3.*m2)*(p[0]+p[3])*p[1]-1./(3.*m2)*(p[1]+M_I*p[2])*p[3])) + (**jit)[1]*((-p[1]-M_I*p[2])*(1./3.-1./(3.*m2)*(-p[1]+M_I*p[2])*p[1]-1./(3.*m2)*(p[0]+p[3])*p[0]) + (p[0]+p[3])*(4./(3.*m2)*p[0]*p[1]-1./(3.*m2)*(p[0]+p[3])*p[1]-1./(3.*m2)*(p[1]+M_I*p[2])*p[0])) - (**jit)[5]*((-p[1]-M_I*p[2])*(-1./(3.*m2)*(-(p[0]+p[3]))*p[1]-1./(3.*m2)*(p[0]+p[3])*p[1]) + (p[0]+p[3])*(1.0+4./(3.*m2)*p[1]*p[1]-1./3.-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[1]-1./(3.*m2)*(p[1]+M_I*p[2])*p[1])) - (**jit)[9]*((-p[1]-M_I*p[2])*(-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[1]-1./(3.*m2)*(p[0]+p[3])*p[2]) + (p[0]+p[3])*(4./(3.*m2)*p[2]*p[1]-M_I/3.-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[1]-1./(3.*m2)*(p[1]+M_I*p[2])*p[2])) - (**jit)[13]*((-p[1]-M_I*p[2])*(-1./3.-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[1]-1./(3.*m2)*(p[0]+p[3])*p[3]) + (p[0]+p[3])*(4./(3.*m2)*p[3]*p[1]-1./(3.*m2)*(p[0]+p[3])*p[1]-1./(3.*m2)*(p[1]+M_I*p[2])*p[3])) + (**jit)[2]*((r*m)*(1./3.-1./(3.*m2)*(-p[1]+M_I*p[2])*p[1]-1./(3.*m2)*(p[0]+p[3])*p[0])) - (**jit)[6]*((r*m)*(-1./(3.*m2)*(-(p[0]+p[3]))*p[1]-1./(3.*m2)*(p[0]+p[3])*p[1])) - (**jit)[10]*((r*m)*(-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[1]-1./(3.*m2)*(p[0]+p[3])*p[2])) - (**jit)[14]*((r*m)*(-1./3.-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[1]-1./(3.*m2)*(p[0]+p[3])*p[3])) + (**jit)[3]*((r*m)*(4./(3.*m2)*p[0]*p[1]-1./(3.*m2)*(p[0]+p[3])*p[1]-1./(3.*m2)*(p[1]+M_I*p[2])*p[0])) - (**jit)[7]*((r*m)*(1.0+4./(3.*m2)*p[1]*p[1]-1./3.-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[1]-1./(3.*m2)*(p[1]+M_I*p[2])*p[1])) - (**jit)[11]*((r*m)*(4./(3.*m2)*p[2]*p[1]-M_I/3.-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[1]-1./(3.*m2)*(p[1]+M_I*p[2])*p[2])) - (**jit)[15]*((r*m)*(4./(3.*m2)*p[3]*p[1]-1./(3.*m2)*(p[0]+p[3])*p[1]-1./(3.*m2)*(p[1]+M_I*p[2])*p[3]));
        j[8] = (**jit)[0]*((r*m)*(4./(3.*m2)*p[0]*p[2]-1./(3.*m2)*(p[0]+p[3])*p[2]-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[0])) - (**jit)[4]*((r*m)*(4./(3.*m2)*p[1]*p[2]-M_I/3.-1./(3.*m2)*(p[1]+M_I*p[2])*p[2]-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[1])) - (**jit)[8]*((r*m)*(1.0+4./(3.*m2)*p[2]*p[2]-1./3.-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[2]-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[2])) - (**jit)[12]*((r*m)*(4./(3.*m2)*p[3]*p[2]-1./(3.*m2)*(p[0]+p[3])*p[2]-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[3])) + (**jit)[1]*((r*m)*(-M_I/3.-1./(3.*m2)*(p[1]+M_I*p[2])*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[0])) - (**jit)[5]*((r*m)*(-1./(3.*m2)*(p[0]+p[3])*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[1])) - (**jit)[9]*((r*m)*(-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[2])) - (**jit)[13]*((r*m)*(M_I/3.-1./(3.*m2)*(-(p[1]+M_I*p[2]))*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[3])) + (**jit)[2]*((p[0]+p[3])*(4./(3.*m2)*p[0]*p[2]-1./(3.*m2)*(p[0]+p[3])*p[2]-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[0]) + (p[1]-M_I*p[2])*(-M_I/3.-1./(3.*m2)*(p[1]+M_I*p[2])*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[0])) - (**jit)[6]*((p[0]+p[3])*(4./(3.*m2)*p[1]*p[2]-M_I/3.-1./(3.*m2)*(p[1]+M_I*p[2])*p[2]-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[1]) + (p[1]-M_I*p[2])*(-1./(3.*m2)*(p[0]+p[3])*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[1])) - (**jit)[10]*((p[0]+p[3])*(1.0+4./(3.*m2)*p[2]*p[2]-1./3.-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[2]-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[2]) + (p[1]-M_I*p[2])*(-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[2])) - (**jit)[14]*((p[0]+p[3])*(4./(3.*m2)*p[3]*p[2]-1./(3.*m2)*(p[0]+p[3])*p[2]-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[3]) + (p[1]-M_I*p[2])*(M_I/3.-1./(3.*m2)*(-(p[1]+M_I*p[2]))*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[3])) + (**jit)[3]*((p[1]+M_I*p[2])*(4./(3.*m2)*p[0]*p[2]-1./(3.*m2)*(p[0]+p[3])*p[2]-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[0]) + (p[0]-p[3])*(-M_I/3.-1./(3.*m2)*(p[1]+M_I*p[2])*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[0])) - (**jit)[7]*((p[1]+M_I*p[2])*(4./(3.*m2)*p[1]*p[2]-M_I/3.-1./(3.*m2)*(p[1]+M_I*p[2])*p[2]-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[1]) + (p[0]-p[3])*(-1./(3.*m2)*(p[0]+p[3])*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[1])) - (**jit)[11]*((p[1]+M_I*p[2])*(1.0+4./(3.*m2)*p[2]*p[2]-1./3.-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[2]-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[2]) + (p[0]-p[3])*(-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[2])) - (**jit)[15]*((p[1]+M_I*p[2])*(4./(3.*m2)*p[3]*p[2]-1./(3.*m2)*(p[0]+p[3])*p[2]-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[3]) + (p[0]-p[3])*(M_I/3.-1./(3.*m2)*(-(p[1]+M_I*p[2]))*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[3]));
        j[9] = (**jit)[0]*((r*m)*(M_I/3.-1./(3.*m2)*(p[1]-M_I*p[2])*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[0])) - (**jit)[4]*((r*m)*(-1./(3.*m2)*(p[0]-p[3])*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[1])) - (**jit)[8]*((r*m)*(-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[2])) - (**jit)[12]*((r*m)*(M_I/3.-1./(3.*m2)*(p[1]-M_I*p[2])*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[3])) + (**jit)[1]*((r*m)*(4./(3.*m2)*p[0]*p[2]-1./(3.*m2)*(p[0]-p[3])*p[2]-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[0])) - (**jit)[5]*((r*m)*(4./(3.*m2)*p[1]*p[2]+M_I/3.-1./(3.*m2)*(p[1]-M_I*p[2])*p[2]-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[1])) - (**jit)[9]*((r*m)*(1.0+4./(3.*m2)*p[2]*p[2]-1./3.-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[2]-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[2])) - (**jit)[13]*((r*m)*(4./(3.*m2)*p[3]*p[2]-1./(3.*m2)*(-(p[0]-p[3]))*p[2]-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[3])) + (**jit)[2]*((p[0]+p[3])*(M_I/3.-1./(3.*m2)*(p[1]-M_I*p[2])*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[0]) + (p[1]-M_I*p[2])*(4./(3.*m2)*p[0]*p[2]-1./(3.*m2)*(p[0]-p[3])*p[2]-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[0])) - (**jit)[6]*((p[0]+p[3])*(-1./(3.*m2)*(p[0]-p[3])*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[1]) + (p[1]-M_I*p[2])*(4./(3.*m2)*p[1]*p[2]+M_I/3.-1./(3.*m2)*(p[1]-M_I*p[2])*p[2]-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[1])) - (**jit)[10]*((p[0]+p[3])*(-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[2]) + (p[1]-M_I*p[2])*(1.0+4./(3.*m2)*p[2]*p[2]-1./3.-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[2]-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[2])) - (**jit)[14]*((p[0]+p[3])*(M_I/3.-1./(3.*m2)*(p[1]-M_I*p[2])*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[3]) + (p[1]-M_I*p[2])*(4./(3.*m2)*p[3]*p[2]-1./(3.*m2)*(-(p[0]-p[3]))*p[2]-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[3])) + (**jit)[3]*((p[1]+M_I*p[2])*(M_I/3.-1./(3.*m2)*(p[1]-M_I*p[2])*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[0]) + (p[0]-p[3])*(4./(3.*m2)*p[0]*p[2]-1./(3.*m2)*(p[0]-p[3])*p[2]-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[0])) - (**jit)[7]*((p[1]+M_I*p[2])*(-1./(3.*m2)*(p[0]-p[3])*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[1]) + (p[0]-p[3])*(4./(3.*m2)*p[1]*p[2]+M_I/3.-1./(3.*m2)*(p[1]-M_I*p[2])*p[2]-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[1])) - (**jit)[11]*((p[1]+M_I*p[2])*(-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[2]) + (p[0]-p[3])*(1.0+4./(3.*m2)*p[2]*p[2]-1./3.-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[2]-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[2])) - (**jit)[15]*((p[1]+M_I*p[2])*(M_I/3.-1./(3.*m2)*(p[1]-M_I*p[2])*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[3]) + (p[0]-p[3])*(4./(3.*m2)*p[3]*p[2]-1./(3.*m2)*(-(p[0]-p[3]))*p[2]-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[3]));
        j[10] = (**jit)[0]*((p[0]-p[3])*(4./(3.*m2)*p[0]*p[2]-1./(3.*m2)*(p[0]-p[3])*p[2]-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[0]) + (-p[1]+M_I*p[2])*(M_I/3.-1./(3.*m2)*(-p[1]-M_I*p[2])*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[0])) - (**jit)[4]*((p[0]-p[3])*(4./(3.*m2)*p[1]*p[2]-M_I/3.-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[2]-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[1]) + (-p[1]+M_I*p[2])*(-1./(3.*m2)*(-(p[0]-p[3]))*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[1])) - (**jit)[8]*((p[0]-p[3])*(1.0+4./(3.*m2)*p[2]*p[2]-1./3.-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[2]-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[2]) + (-p[1]+M_I*p[2])*(-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[2])) - (**jit)[12]*((p[0]-p[3])*(4./(3.*m2)*p[3]*p[2]-1./(3.*m2)*(-(p[0]-p[3]))*p[2]-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[3]) + (-p[1]+M_I*p[2])*(M_I/3.-1./(3.*m2)*(-p[1]-M_I*p[2])*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[3])) + (**jit)[1]*((-p[1]-M_I*p[2])*(4./(3.*m2)*p[0]*p[2]-1./(3.*m2)*(p[0]-p[3])*p[2]-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[0]) + (p[0]+p[3])*(M_I/3.-1./(3.*m2)*(-p[1]-M_I*p[2])*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[0])) - (**jit)[5]*((-p[1]-M_I*p[2])*(4./(3.*m2)*p[1]*p[2]-M_I/3.-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[2]-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[1]) + (p[0]+p[3])*(-1./(3.*m2)*(-(p[0]-p[3]))*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[1])) - (**jit)[9]*((-p[1]-M_I*p[2])*(1.0+4./(3.*m2)*p[2]*p[2]-1./3.-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[2]-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[2]) + (p[0]+p[3])*(-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[2])) - (**jit)[13]*((-p[1]-M_I*p[2])*(4./(3.*m2)*p[3]*p[2]-1./(3.*m2)*(-(p[0]-p[3]))*p[2]-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[3]) + (p[0]+p[3])*(M_I/3.-1./(3.*m2)*(-p[1]-M_I*p[2])*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[3])) + (**jit)[2]*((r*m)*(4./(3.*m2)*p[0]*p[2]-1./(3.*m2)*(p[0]-p[3])*p[2]-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[0])) - (**jit)[6]*((r*m)*(4./(3.*m2)*p[1]*p[2]-M_I/3.-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[2]-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[1])) - (**jit)[10]*((r*m)*(1.0+4./(3.*m2)*p[2]*p[2]-1./3.-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[2]-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[2])) - (**jit)[14]*((r*m)*(4./(3.*m2)*p[3]*p[2]-1./(3.*m2)*(-(p[0]-p[3]))*p[2]-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[3])) + (**jit)[3]*((r*m)*(M_I/3.-1./(3.*m2)*(-p[1]-M_I*p[2])*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[0])) - (**jit)[7]*((r*m)*(-1./(3.*m2)*(-(p[0]-p[3]))*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[1])) - (**jit)[11]*((r*m)*(-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[2])) - (**jit)[15]*((r*m)*(M_I/3.-1./(3.*m2)*(-p[1]-M_I*p[2])*p[2]-1./(3.*m2)*(M_I*(p[0]-p[3]))*p[3]));
        j[11] = (**jit)[0]*((p[0]-p[3])*(-M_I/3.-1./(3.*m2)*(-p[1]+M_I*p[2])*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[0]) + (-p[1]+M_I*p[2])*(4./(3.*m2)*p[0]*p[2]-1./(3.*m2)*(p[0]+p[3])*p[2]-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[0])) - (**jit)[4]*((p[0]-p[3])*(-1./(3.*m2)*(-(p[0]+p[3]))*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[1]) + (-p[1]+M_I*p[2])*(4./(3.*m2)*p[1]*p[2]+M_I/3.-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[2]-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[1])) - (**jit)[8]*((p[0]-p[3])*(-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[2]) + (-p[1]+M_I*p[2])*(1.0+4./(3.*m2)*p[2]*p[2]-1./3.-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[2]-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[2])) - (**jit)[12]*((p[0]-p[3])*(M_I/3.-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[3]) + (-p[1]+M_I*p[2])*(4./(3.*m2)*p[3]*p[2]-1./(3.*m2)*(p[0]+p[3])*p[2]-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[3])) + (**jit)[1]*((-p[1]-M_I*p[2])*(-M_I/3.-1./(3.*m2)*(-p[1]+M_I*p[2])*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[0]) + (p[0]+p[3])*(4./(3.*m2)*p[0]*p[2]-1./(3.*m2)*(p[0]+p[3])*p[2]-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[0])) - (**jit)[5]*((-p[1]-M_I*p[2])*(-1./(3.*m2)*(-(p[0]+p[3]))*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[1]) + (p[0]+p[3])*(4./(3.*m2)*p[1]*p[2]+M_I/3.-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[2]-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[1])) - (**jit)[9]*((-p[1]-M_I*p[2])*(-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[2]) + (p[0]+p[3])*(1.0+4./(3.*m2)*p[2]*p[2]-1./3.-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[2]-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[2])) - (**jit)[13]*((-p[1]-M_I*p[2])*(M_I/3.-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[3]) + (p[0]+p[3])*(4./(3.*m2)*p[3]*p[2]-1./(3.*m2)*(p[0]+p[3])*p[2]-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[3])) + (**jit)[2]*((r*m)*(-M_I/3.-1./(3.*m2)*(-p[1]+M_I*p[2])*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[0])) - (**jit)[6]*((r*m)*(-1./(3.*m2)*(-(p[0]+p[3]))*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[1])) - (**jit)[10]*((r*m)*(-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[2])) - (**jit)[14]*((r*m)*(M_I/3.-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[2]-1./(3.*m2)*(-M_I*(p[0]+p[3]))*p[3])) + (**jit)[3]*((r*m)*(4./(3.*m2)*p[0]*p[2]-1./(3.*m2)*(p[0]+p[3])*p[2]-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[0])) - (**jit)[7]*((r*m)*(4./(3.*m2)*p[1]*p[2]+M_I/3.-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[2]-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[1])) - (**jit)[11]*((r*m)*(1.0+4./(3.*m2)*p[2]*p[2]-1./3.-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[2]-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[2])) - (**jit)[15]*((r*m)*(4./(3.*m2)*p[3]*p[2]-1./(3.*m2)*(p[0]+p[3])*p[2]-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[3]));
        j[12] = (**jit)[0]*((r*m)*(4./(3.*m2)*p[0]*p[3]-1./3.-1./(3.*m2)*(p[0]+p[3])*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[0])) - (**jit)[4]*((r*m)*(4./(3.*m2)*p[1]*p[3]-1./(3.*m2)*(p[1]+M_I*p[2])*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[1])) - (**jit)[8]*((r*m)*(4./(3.*m2)*p[2]*p[3]-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[2])) - (**jit)[12]*((r*m)*(1.0+4./(3.*m2)*p[3]*p[3]-1./3.-1./(3.*m2)*(p[0]+p[3])*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[3])) + (**jit)[1]*((r*m)*(-1./(3.*m2)*(p[1]+M_I*p[2])*p[3]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[0])) - (**jit)[5]*((r*m)*(-1./3.-1./(3.*m2)*(p[0]+p[3])*p[3]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[1])) - (**jit)[9]*((r*m)*(-M_I/3.-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[3]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[2])) - (**jit)[13]*((r*m)*(-1./(3.*m2)*(-(p[1]+M_I*p[2]))*p[3]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[3])) + (**jit)[2]*((p[0]+p[3])*(4./(3.*m2)*p[0]*p[3]-1./3.-1./(3.*m2)*(p[0]+p[3])*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[0]) + (p[1]-M_I*p[2])*(-1./(3.*m2)*(p[1]+M_I*p[2])*p[3]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[0])) - (**jit)[6]*((p[0]+p[3])*(4./(3.*m2)*p[1]*p[3]-1./(3.*m2)*(p[1]+M_I*p[2])*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[1]) + (p[1]-M_I*p[2])*(-1./3.-1./(3.*m2)*(p[0]+p[3])*p[3]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[1])) - (**jit)[10]*((p[0]+p[3])*(4./(3.*m2)*p[2]*p[3]-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[2]) + (p[1]-M_I*p[2])*(-M_I/3.-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[3]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[2])) - (**jit)[14]*((p[0]+p[3])*(1.0+4./(3.*m2)*p[3]*p[3]-1./3.-1./(3.*m2)*(p[0]+p[3])*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[3]) + (p[1]-M_I*p[2])*(-1./(3.*m2)*(-(p[1]+M_I*p[2]))*p[3]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[3])) + (**jit)[3]*((p[1]+M_I*p[2])*(4./(3.*m2)*p[0]*p[3]-1./3.-1./(3.*m2)*(p[0]+p[3])*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[0]) + (p[0]-p[3])*(-1./(3.*m2)*(p[1]+M_I*p[2])*p[3]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[0])) - (**jit)[7]*((p[1]+M_I*p[2])*(4./(3.*m2)*p[1]*p[3]-1./(3.*m2)*(p[1]+M_I*p[2])*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[1]) + (p[0]-p[3])*(-1./3.-1./(3.*m2)*(p[0]+p[3])*p[3]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[1])) - (**jit)[11]*((p[1]+M_I*p[2])*(4./(3.*m2)*p[2]*p[3]-1./(3.*m2)*(-M_I*(p[1]+M_I*p[2]))*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[2]) + (p[0]-p[3])*(-M_I/3.-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[3]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[2])) - (**jit)[15]*((p[1]+M_I*p[2])*(1.0+4./(3.*m2)*p[3]*p[3]-1./3.-1./(3.*m2)*(p[0]+p[3])*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[3]) + (p[0]-p[3])*(-1./(3.*m2)*(-(p[1]+M_I*p[2]))*p[3]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[3]));
        j[13] = (**jit)[0]*((r*m)*(-1./(3.*m2)*(p[1]-M_I*p[2])*p[3]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[0])) - (**jit)[4]*((r*m)*(1./3.-1./(3.*m2)*(p[0]-p[3])*p[3]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[1])) - (**jit)[8]*((r*m)*(-M_I/3.-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[3]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[2])) - (**jit)[12]*((r*m)*(-1./(3.*m2)*(p[1]-M_I*p[2])*p[3]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[3])) + (**jit)[1]*((r*m)*(4./(3.*m2)*p[0]*p[3]+1./3.-1./(3.*m2)*(p[0]-p[3])*p[3]-1./(3.*m2)*(p[0]+p[3])*p[0])) - (**jit)[5]*((r*m)*(4./(3.*m2)*p[1]*p[3]-1./(3.*m2)*(p[1]-M_I*p[2])*p[3]-1./(3.*m2)*(p[0]+p[3])*p[1])) - (**jit)[9]*((r*m)*(4./(3.*m2)*p[2]*p[3]-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[3]-1./(3.*m2)*(p[0]+p[3])*p[2])) - (**jit)[13]*((r*m)*(1.0+4./(3.*m2)*p[3]*p[3]-1./3.-1./(3.*m2)*(-(p[0]-p[3]))*p[3]-1./(3.*m2)*(p[0]+p[3])*p[3])) + (**jit)[2]*((p[0]+p[3])*(-1./(3.*m2)*(p[1]-M_I*p[2])*p[3]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[0]) + (p[1]-M_I*p[2])*(4./(3.*m2)*p[0]*p[3]+1./3.-1./(3.*m2)*(p[0]-p[3])*p[3]-1./(3.*m2)*(p[0]+p[3])*p[0])) - (**jit)[6]*((p[0]+p[3])*(1./3.-1./(3.*m2)*(p[0]-p[3])*p[3]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[1]) + (p[1]-M_I*p[2])*(4./(3.*m2)*p[1]*p[3]-1./(3.*m2)*(p[1]-M_I*p[2])*p[3]-1./(3.*m2)*(p[0]+p[3])*p[1])) - (**jit)[10]*((p[0]+p[3])*(-M_I/3.-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[3]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[2]) + (p[1]-M_I*p[2])*(4./(3.*m2)*p[2]*p[3]-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[3]-1./(3.*m2)*(p[0]+p[3])*p[2])) - (**jit)[14]*((p[0]+p[3])*(-1./(3.*m2)*(p[1]-M_I*p[2])*p[3]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[3]) + (p[1]-M_I*p[2])*(1.0+4./(3.*m2)*p[3]*p[3]-1./3.-1./(3.*m2)*(-(p[0]-p[3]))*p[3]-1./(3.*m2)*(p[0]+p[3])*p[3])) + (**jit)[3]*((p[1]+M_I*p[2])*(-1./(3.*m2)*(p[1]-M_I*p[2])*p[3]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[0]) + (p[0]-p[3])*(4./(3.*m2)*p[0]*p[3]+1./3.-1./(3.*m2)*(p[0]-p[3])*p[3]-1./(3.*m2)*(p[0]+p[3])*p[0])) - (**jit)[7]*((p[1]+M_I*p[2])*(1./3.-1./(3.*m2)*(p[0]-p[3])*p[3]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[1]) + (p[0]-p[3])*(4./(3.*m2)*p[1]*p[3]-1./(3.*m2)*(p[1]-M_I*p[2])*p[3]-1./(3.*m2)*(p[0]+p[3])*p[1])) - (**jit)[11]*((p[1]+M_I*p[2])*(-M_I/3.-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[3]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[2]) + (p[0]-p[3])*(4./(3.*m2)*p[2]*p[3]-1./(3.*m2)*(M_I*(p[1]-M_I*p[2]))*p[3]-1./(3.*m2)*(p[0]+p[3])*p[2])) - (**jit)[15]*((p[1]+M_I*p[2])*(-1./(3.*m2)*(p[1]-M_I*p[2])*p[3]-1./(3.*m2)*(-p[1]+M_I*p[2])*p[3]) + (p[0]-p[3])*(1.0+4./(3.*m2)*p[3]*p[3]-1./3.-1./(3.*m2)*(-(p[0]-p[3]))*p[3]-1./(3.*m2)*(p[0]+p[3])*p[3]));
        j[14] = (**jit)[0]*((p[0]-p[3])*(4./(3.*m2)*p[0]*p[3]+1./3.-1./(3.*m2)*(p[0]-p[3])*p[3]-1./(3.*m2)*(p[0]+p[3])*p[0]) + (-p[1]+M_I*p[2])*(-1./(3.*m2)*(-p[1]-M_I*p[2])*p[3]-1./(3.*m2)*(p[1]+M_I*p[2])*p[0])) - (**jit)[4]*((p[0]-p[3])*(4./(3.*m2)*p[1]*p[3]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[3]-1./(3.*m2)*(p[0]+p[3])*p[1]) + (-p[1]+M_I*p[2])*(-1./3.-1./(3.*m2)*(-(p[0]-p[3]))*p[3]-1./(3.*m2)*(p[1]+M_I*p[2])*p[1])) - (**jit)[8]*((p[0]-p[3])*(4./(3.*m2)*p[2]*p[3]-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[3]-1./(3.*m2)*(p[0]+p[3])*p[2]) + (-p[1]+M_I*p[2])*(-M_I/3.-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[3]-1./(3.*m2)*(p[1]+M_I*p[2])*p[2])) - (**jit)[12]*((p[0]-p[3])*(1.0+4./(3.*m2)*p[3]*p[3]-1./3.-1./(3.*m2)*(-(p[0]-p[3]))*p[3]-1./(3.*m2)*(p[0]+p[3])*p[3]) + (-p[1]+M_I*p[2])*(-1./(3.*m2)*(-p[1]-M_I*p[2])*p[3]-1./(3.*m2)*(p[1]+M_I*p[2])*p[3])) + (**jit)[1]*((-p[1]-M_I*p[2])*(4./(3.*m2)*p[0]*p[3]+1./3.-1./(3.*m2)*(p[0]-p[3])*p[3]-1./(3.*m2)*(p[0]+p[3])*p[0]) + (p[0]+p[3])*(-1./(3.*m2)*(-p[1]-M_I*p[2])*p[3]-1./(3.*m2)*(p[1]+M_I*p[2])*p[0])) - (**jit)[5]*((-p[1]-M_I*p[2])*(4./(3.*m2)*p[1]*p[3]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[3]-1./(3.*m2)*(p[0]+p[3])*p[1]) + (p[0]+p[3])*(-1./3.-1./(3.*m2)*(-(p[0]-p[3]))*p[3]-1./(3.*m2)*(p[1]+M_I*p[2])*p[1])) - (**jit)[9]*((-p[1]-M_I*p[2])*(4./(3.*m2)*p[2]*p[3]-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[3]-1./(3.*m2)*(p[0]+p[3])*p[2]) + (p[0]+p[3])*(-M_I/3.-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[3]-1./(3.*m2)*(p[1]+M_I*p[2])*p[2])) - (**jit)[13]*((-p[1]-M_I*p[2])*(1.0+4./(3.*m2)*p[3]*p[3]-1./3.-1./(3.*m2)*(-(p[0]-p[3]))*p[3]-1./(3.*m2)*(p[0]+p[3])*p[3]) + (p[0]+p[3])*(-1./(3.*m2)*(-p[1]-M_I*p[2])*p[3]-1./(3.*m2)*(p[1]+M_I*p[2])*p[3])) + (**jit)[2]*((r*m)*(4./(3.*m2)*p[0]*p[3]+1./3.-1./(3.*m2)*(p[0]-p[3])*p[3]-1./(3.*m2)*(p[0]+p[3])*p[0])) - (**jit)[6]*((r*m)*(4./(3.*m2)*p[1]*p[3]-1./(3.*m2)*(-(-p[1]-M_I*p[2]))*p[3]-1./(3.*m2)*(p[0]+p[3])*p[1])) - (**jit)[10]*((r*m)*(4./(3.*m2)*p[2]*p[3]-1./(3.*m2)*(M_I*(-p[1]-M_I*p[2]))*p[3]-1./(3.*m2)*(p[0]+p[3])*p[2])) - (**jit)[14]*((r*m)*(1.0+4./(3.*m2)*p[3]*p[3]-1./3.-1./(3.*m2)*(-(p[0]-p[3]))*p[3]-1./(3.*m2)*(p[0]+p[3])*p[3])) + (**jit)[3]*((r*m)*(-1./(3.*m2)*(-p[1]-M_I*p[2])*p[3]-1./(3.*m2)*(p[1]+M_I*p[2])*p[0])) - (**jit)[7]*((r*m)*(-1./3.-1./(3.*m2)*(-(p[0]-p[3]))*p[3]-1./(3.*m2)*(p[1]+M_I*p[2])*p[1])) - (**jit)[11]*((r*m)*(-M_I/3.-1./(3.*m2)*(-M_I*(p[0]-p[3]))*p[3]-1./(3.*m2)*(p[1]+M_I*p[2])*p[2])) - (**jit)[15]*((r*m)*(-1./(3.*m2)*(-p[1]-M_I*p[2])*p[3]-1./(3.*m2)*(p[1]+M_I*p[2])*p[3]));
        j[15] = (**jit)[0]*((p[0]-p[3])*(-1./(3.*m2)*(-p[1]+M_I*p[2])*p[3]-1./(3.*m2)*(-(p[1]-M_I*p[2]))*p[0]) + (-p[1]+M_I*p[2])*(4./(3.*m2)*p[0]*p[3]-1./3.-1./(3.*m2)*(p[0]+p[3])*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[0])) - (**jit)[4]*((p[0]-p[3])*(1./3.-1./(3.*m2)*(-(p[0]+p[3]))*p[3]-1./(3.*m2)*(-(p[1]-M_I*p[2]))*p[1]) + (-p[1]+M_I*p[2])*(4./(3.*m2)*p[1]*p[3]-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[1])) - (**jit)[8]*((p[0]-p[3])*(-M_I/3.-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[3]-1./(3.*m2)*(-(p[1]-M_I*p[2]))*p[2]) + (-p[1]+M_I*p[2])*(4./(3.*m2)*p[2]*p[3]-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[2])) - (**jit)[12]*((p[0]-p[3])*(-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[3]-1./(3.*m2)*(-(p[1]-M_I*p[2]))*p[3]) + (-p[1]+M_I*p[2])*(1.0+4./(3.*m2)*p[3]*p[3]-1./3.-1./(3.*m2)*(p[0]+p[3])*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[3])) + (**jit)[1]*((-p[1]-M_I*p[2])*(-1./(3.*m2)*(-p[1]+M_I*p[2])*p[3]-1./(3.*m2)*(-(p[1]-M_I*p[2]))*p[0]) + (p[0]+p[3])*(4./(3.*m2)*p[0]*p[3]-1./3.-1./(3.*m2)*(p[0]+p[3])*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[0])) - (**jit)[5]*((-p[1]-M_I*p[2])*(1./3.-1./(3.*m2)*(-(p[0]+p[3]))*p[3]-1./(3.*m2)*(-(p[1]-M_I*p[2]))*p[1]) + (p[0]+p[3])*(4./(3.*m2)*p[1]*p[3]-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[1])) - (**jit)[9]*((-p[1]-M_I*p[2])*(-M_I/3.-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[3]-1./(3.*m2)*(-(p[1]-M_I*p[2]))*p[2]) + (p[0]+p[3])*(4./(3.*m2)*p[2]*p[3]-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[2])) - (**jit)[13]*((-p[1]-M_I*p[2])*(-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[3]-1./(3.*m2)*(-(p[1]-M_I*p[2]))*p[3]) + (p[0]+p[3])*(1.0+4./(3.*m2)*p[3]*p[3]-1./3.-1./(3.*m2)*(p[0]+p[3])*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[3])) + (**jit)[2]*((r*m)*(-1./(3.*m2)*(-p[1]+M_I*p[2])*p[3]-1./(3.*m2)*(-(p[1]-M_I*p[2]))*p[0])) - (**jit)[6]*((r*m)*(1./3.-1./(3.*m2)*(-(p[0]+p[3]))*p[3]-1./(3.*m2)*(-(p[1]-M_I*p[2]))*p[1])) - (**jit)[10]*((r*m)*(-M_I/3.-1./(3.*m2)*(M_I*(p[0]+p[3]))*p[3]-1./(3.*m2)*(-(p[1]-M_I*p[2]))*p[2])) - (**jit)[14]*((r*m)*(-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[3]-1./(3.*m2)*(-(p[1]-M_I*p[2]))*p[3])) + (**jit)[3]*((r*m)*(4./(3.*m2)*p[0]*p[3]-1./3.-1./(3.*m2)*(p[0]+p[3])*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[0])) - (**jit)[7]*((r*m)*(4./(3.*m2)*p[1]*p[3]-1./(3.*m2)*(-(-p[1]+M_I*p[2]))*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[1])) - (**jit)[11]*((r*m)*(4./(3.*m2)*p[2]*p[3]-1./(3.*m2)*(-M_I*(-p[1]+M_I*p[2]))*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[2])) - (**jit)[15]*((r*m)*(1.0+4./(3.*m2)*p[3]*p[3]-1./3.-1./(3.*m2)*(p[0]+p[3])*p[3]-1./(3.*m2)*(-(p[0]-p[3]))*p[3]));
      }
      **jit=j*prop;
    }
  }
}

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
      DEBUG_VAR(**jit1);
      DEBUG_VAR(**jit2);
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
bool CRS<SType>::Test_WF_Properties(const ATOOLS::Vec4D &p, const bool &anti, const int &dir) {
  CRaritaSchwinger<SType> rspp = anti?RSPP(p, -1, 1, 1, 0, 0, 1):RSPP(p, 1, 1, 1, 0, 0, 1);
  CRaritaSchwinger<SType> rsmm = anti?RSMM(p, -1, 1, 1, 0, 0, -1):RSMM(p, 1, 1, 1, 0, 0, -1);
  METOOLS::Gamma<SType> gammavec = Gamma<SType>();
  bool global_testresult(true);
  bool local_testresult(true);

  // properties of massless RS wave functions
  if (!this->m_msv){
    // normalization
    msg_Debugging()<<METHOD<<": Testing normalization of Rarita-Schwinger wave function..."<<"\n";
    TCMatrix<SType> result1 = rspp.Contract4Index(rspp.Bar()) + rsmm.Contract4Index(rsmm.Bar()) + gammavec * p;
    for (size_t i(0); i<4; ++i){
      for (size_t j(0); j<4; ++j){
        if (std::abs(result1[i][j].real())> rspp.Accu()|| std::abs(result1[i][j].imag())>rspp.Accu()) {
          msg_Debugging()<<"Component " << i << j << " of resulting 4x4 matrix is " << result1[i][j] << ", instead of zero!"
          << "\n";
          global_testresult = false;
        }
      }
    }
    if (global_testresult) msg_Debugging() << "passed" << "\n";

    // equality between U++ and V-- / U-- and V++ for bar and non-bar
    msg_Debugging()<<METHOD<<": Testing equality of U++ and V-- / U-- and V++ Rarita-Schwinger wave functions..."<<"\n";
    CRaritaSchwinger<SType> rspp1 = anti?RSPP(p, 1, 1, 1, 0, 0, 1):RSPP(p, -1, 1, 1, 0, 0, 1);
    CRaritaSchwinger<SType> rsmm1 = anti?RSMM(p, 1, 1, 1, 0, 0, -1):RSMM(p, -1, 1, 1, 0, 0, -1);
    for (size_t i(0); i<16; ++i){
      if (std::abs((rspp1[i]-rsmm[i]).real()) > rsmm.Accu()|| std::abs((rspp1[i]-rsmm[i]).imag()) > rsmm.Accu()) {
        if (anti){
          msg_Debugging()<<"Components " << i << " of the Rarita-Schwinger wave functions " << "U++ / V--" << rspp1[i] <<
                         "/" << rsmm[i] << " are not equal!" << "\n";}
        else{
          msg_Debugging()<<"Components " << i << " of the Rarita-Schwinger wave functions " << "U-- / V++" << rsmm[i] <<
          "/" << rspp1[i] << " are not equal!" << "\n";}
        global_testresult = false; local_testresult = false;
      }
      if (std::abs((rsmm1[i]-rspp[i]).real()) > rspp.Accu()|| std::abs((rsmm1[i]-rspp[i]).imag()) > rspp.Accu()) {
        if (anti){
          msg_Debugging()<<"Components " << i << " of the Rarita-Schwinger wave functions " << "U-- / V++" << rsmm1[i]
                         << "/" << rspp[i] << " are not equal!" << "\n";
        }
        else{
          msg_Debugging()<<"Components " << i << " of the Rarita-Schwinger wave functions " <<
                         "U++ / V--" << rspp[i] << "/" << rsmm1[i] << " are not equal!" << "\n";
        }
        global_testresult = false; local_testresult = false;
      }
    }
    if (local_testresult) msg_Debugging() << "passed" << "\n";
  }
  else{
    // properties of massive RS wave functions
    // normalization
    msg_Debugging()<<METHOD<<": Testing normalization of Rarita-Schwinger wave function..."<<"\n";
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
        msg_Debugging()<<"Normalization of the Rarita-Schwinger wave functions" << i <<
        "do not fit! 0=++bar*++, 1=++bar*+, 2=++bar*-, 3=++bar*--, 4=+bar*++ ..." << "\n";
        global_testresult = false;
      }
    }
    if (global_testresult) msg_Debugging() << "passed" << "\n";

    // completeness relation according to Hagiwara et al. Eur. Phys. J. C (2011) 71: 1529
    msg_Debugging()<<METHOD<<": Testing completeness of Rarita-Schwinger wave functions..."<<"\n";
    ATOOLS::TCMatrix<SType> p_slash = gammavec*p;
    //TODO: Stimmt das dir an dieser Stelle???
    ATOOLS::TCMatrix<SType> spropagator = ATOOLS::TCMatrix<SType>(p_slash + (-dir) * (anti?-1.0:1.0) * SComplex(p.Mass()) *
                                                                                      ATOOLS::TCMatrix<SType>(4, true));
    // --------------------------------------------------------------------------------------------------------------
    // Sum over Lorentz index
    /*TCMatrix<SType> left_sum_red =  rspp.Contract4Index(rspp.Bar()) + rsp.Contract4Index(rsp.Bar()) + rsm.Contract4Index(rsm.Bar())
    + rsmm.Contract4Index(rsmm.Bar());
    for (size_t k(0); k<4; ++k){
      for (size_t l(0); l<4; ++l){
        std::cout << left_sum_red[k][l] << "  " << (spropagator*TCMatrix<SType>(4, true)*2.)[k][l] << "\n";
        if (std::abs((left_sum_red[k][l]+(spropagator*TCMatrix<SType>(4, true)*2.)[k][l]).real()) > rspp.Accu() ||
          std::abs((left_sum_red[k][l]+(spropagator*TCMatrix<SType>(4, true)*2.)[k][l]).imag()) > rspp.Accu()) {
          msg_Debugging() << "Completeness relation of the Rarita-Schwinger wave functions is not hold: "
                       "component " << k << " " << l << "of the resulting 16 dimensional tensor is " <<
                    left_sum_red[k][l]+(spropagator*TCMatrix<SType>(4, true)*2.)[k][l] << " instead of zero!" << "\n";}}}*/
      // --------------------------------------------------------------------------------------------------------------
    // Full completeness relation
    for (size_t mu(0); mu < 4; ++mu){
      for (size_t nu(0); nu < 4; ++nu){
        for (size_t A(0); A < 4; ++A){
          for (size_t B(0); B < 4; ++B){
            // A,B spinor indexes; mu,nu lorentz indexes
            // helicity sum
            SComplex left_sum = rspp[A + 4 * mu] * rspp.Bar()[B + 4 * nu] + rsp[A + 4 * mu] * rsp.Bar()[B + 4 * nu]
              + rsm[A + 4 * mu] * rsm.Bar()[B + 4 * nu] + rsmm[A + 4 * mu] * rsmm.Bar()[B + 4 * nu];
            SComplex propagator(0., 0.);
            // components of the lorentz tensor part
            for (size_t C(0); C < 4; ++C){
              SComplex comp_lt_ten(0., 0.);
              comp_lt_ten += (mu == nu ? (mu > 0 ? 1. : -1.) : 0.) * (C == B ? 1.0 : 0.0); // metrischer Tensor -g_{\mu \nu} * 1_{CB}
              comp_lt_ten += (4. / (3. * p.Abs2())) * p[mu] * p[nu] * (C == B ? 1.0 : 0.0); // (4/(3 C^2))* p_{\mu} * p_{\D u} * 1_{CB}
              comp_lt_ten += (1. / 3.) * (gammavec[mu] * gammavec[nu])[C][B];  // (4/3) \gamma_{\mu, CE} \gamma_{\nu, EB}
              comp_lt_ten -= (1. / 3.) * (p[nu] / p.Abs2()) * (gammavec[mu] * p_slash)[C][B];
              comp_lt_ten -= (1. / 3.) * (p[mu] / p.Abs2()) * (p_slash * gammavec[nu])[C][B];
              // propagator is spinor propagator[A][C] * lorentz tensor[C][B] for fixed lorentz indizes mu, nu
              propagator += spropagator[A][C]*comp_lt_ten;
            }
            // TODO: Stimmt das mit dem dir hier?
            if (std::abs((left_sum - SComplex(-dir)*propagator).real()) > rspp.Accu() ||
            std::abs((left_sum - SComplex(-dir)*propagator).imag()) > rspp.Accu()) {
              msg_Debugging() << "Completeness relation of the Rarita-Schwinger wave functions is not fulfilled: "
                           "component " << A << " " << B << " " << mu << " " << nu << " of the resulting 16 dimensional tensor is " <<
                        left_sum - SComplex(-dir)*propagator << " instead of zero!" << "\n";
              global_testresult = false; local_testresult = false;
            }
          }
        }
      }
    }
    if (local_testresult) msg_Debugging() << "passed" << "\n";
  }

  // gauge invariance
  return global_testresult;
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

/* PYTHON SCRIPT TO DETERMINE THE COMPONENTS OF THE PROPAGATOR NUMERATOR
 * resulting components were first tested and were
 * subsequently further simplified by Mathematica; the result can be checked against the original expression from the
 * python script whether the same values are achieved
 *
 * import numpy as np

# gamma matrices in Sherpas Weyl basis
gamma0 = [[0., 0., 1.0, 0.], [0., 0., 0., 1.], [1., 0., 0., 0.], [0., 1., 0., 0.]]
gamma1 = [[0., 0., 0., 1.], [0., 0., 1., 0.], [0., -1., 0., 0.], [-1., 0., 0., 0.]]
gamma2 = [[0., 0., 0., -1j], [0., 0., 1j, 0.], [0., 1j, 0., 0.], [-1j, 0., 0., 0.]]
gamma3 = [[0., 0., 1.0, 0.], [0., 0., 0., -1.], [-1., 0., 0., 0.], [0., 1., 0., 0.]]
gamma = np.array([gamma0, gamma1, gamma2, gamma3], dtype=complex)

# RS wave function to contract with (written down in matrix form psi[A][mu] for a simpler calculation, components
# corresponds to the corresponding components in the vector notation of the wave function in Sherpa
# within the whole script, upper case latin letters denote spinor, greek letters Lorentz indices
psi = [["(**jit)[0]", "(**jit)[4]", "(**jit)[8]", "(**jit)[12]"],
 ["(**jit)[1]", "(**jit)[5]", "(**jit)[9]", "(**jit)[13]"], ["(**jit)[2]", "(**jit)[6]", "(**jit)[10]", "(**jit)[14]"],
 ["(**jit)[3]", "(**jit)[7]", "(**jit)[11]", "(**jit)[15]"]]
# variable name for the imaginary unit
imag_var = "M_I"

# gmunu-Term
g_term = np.zeros(shape=(4, 4, 4, 4))
for i in range(4):
    if i == 0:
        prop = -1
    else:
        prop = 1
    for j in range(4):
        g_term[j][j][i][i] = 1 * prop

# 1/3 gamma_mu gamma_nu term
gamma_prod = np.zeros(shape=(4, 4, 4, 4), dtype=complex)
for A in range(4):
    for mu in range(4):
        for B in range(4):
            for nu in range(4):
                for C in range(4):
                    gamma_prod[A][B][mu][nu] += gamma[mu][A][C] * gamma[nu][C][B]

# pmu pnu term
p = ["p[0]", "p[1]", "p[2]", "p[3]"]
p_term = [[[['' for _ in range(4)] for _ in range(4)] for _ in range(4)] for _ in range(4)]
for A in range(4):
    for mu in range(4):
        for B in range(4):
            for nu in range(4):
                if A==B:
                    p_term[A][B][mu][nu] = p[mu] + "*" + p[nu]

# gamma p terms (p dagger)
gamma_p = [['' for _ in range(4)] for _ in range(4)]
for A in range(4):
    for mu in range(4):
        for B in range(4):
            if mu == 0:
                prop = 1
            else:
                prop = -1
            if gamma[mu][A][B] != complex(0.0):
                if gamma[mu][A][B] * prop == complex(1., 0.):
                    if gamma_p[A][B] != '':
                        gamma_p[A][B] += "+"
                    gamma_p[A][B] += p[mu]
                elif gamma[mu][A][B] * prop == complex(-1., 0.):
                    gamma_p[A][B] += "-" + p[mu]
                elif gamma[mu][A][B] * prop == complex(0., 1.):
                    if gamma_p[A][B] != '':
                        gamma_p[A][B] += "+"
                    gamma_p[A][B] += imag_var + "*" + p[mu]
                elif gamma[mu][A][B] * prop == complex(0., -1.):
                    gamma_p[A][B] += "-" + imag_var + "*" + p[mu]
                else:
                    if gamma_p[A][B] != '':
                        gamma_p[A][B] += "+"
                    gamma_p[A][B] += str(gamma[mu][A][B] * prop) + "*" + p[mu]

# gamma_mu p_nu term
gammamu_pdagger = [[['' for _ in range(4)] for _ in range(4)] for _ in range(4)]
for A in range(4):
    for mu in range(4):
        for B in range(4):
            for C in range(4):
                if gamma[mu][A][C] != complex(0.) and gamma_p[C][B] != "":
                    if gamma[mu][A][C] == complex(1., 0.):
                        if gammamu_pdagger[A][B][mu] != '':
                            gammamu_pdagger[A][B][mu] += "+"
                        gammamu_pdagger[A][B][mu] += gamma_p[C][B]
                    elif gamma[mu][A][C] == complex(-1., 0.):
                        gammamu_pdagger[A][B][mu] += "-" + "(" + gamma_p[C][B] + ")"
                    elif gamma[mu][A][C] == complex(0., 1.) or gamma[mu][A][C] == complex(-0., 1.):
                        if gammamu_pdagger[A][B][mu] != '':
                            gammamu_pdagger[A][B][mu] += "+"
                        gammamu_pdagger[A][B][mu] += imag_var + "*" + "(" + gamma_p[C][B] + ")"
                    elif gamma[mu][A][C] == complex(0., -1.) or gamma[mu][A][C] == complex(-0., -1.):
                        gammamu_pdagger[A][B][mu] += "-" + imag_var + "*" + "(" + gamma_p[C][B] + ")"
                    else:
                        if gammamu_pdagger[A][B][mu] != '':
                            gammamu_pdagger[A][B][mu] += "+"
                        gammamu_pdagger[A][B][mu] += str(gamma[mu][A][C]) + "*" + "(" + gamma_p[C][B] + ")"

# gamma_mu p_nu term
pdagger_gammanu = [[['' for _ in range(4)] for _ in range(4)] for _ in range(4)]
for A in range(4):
    for nu in range(4):
        for B in range(4):
            for C in range(4):
                if gamma_p[A][C] != "" and gamma[nu][C][B] != complex(0., 0.):
                    if gamma[nu][C][B] == complex(1., 0.):
                        if pdagger_gammanu[A][B][nu] != "":
                            pdagger_gammanu[A][B][nu] += "+"
                        pdagger_gammanu[A][B][nu] += gamma_p[A][C]
                    elif gamma[nu][C][B] == complex(-1., 0.):
                        pdagger_gammanu[A][B][nu] += "-" + "(" + gamma_p[A][C] + ")"
                    elif gamma[nu][C][B] == complex(0., 1.):
                        if pdagger_gammanu[A][B][nu] != "":
                            pdagger_gammanu[A][B][nu] += "+"
                        pdagger_gammanu[A][B][nu] += imag_var + "*" + "(" + gamma_p[A][C] + ")"
                    elif gamma[nu][C][B] == complex(0., -1.):
                        pdagger_gammanu[A][B][nu] += "-" + imag_var + "*" + "(" + gamma_p[A][C] + ")"

# fermion propagator part
# r=1 for particles, r=-1 for anti-particles
spropagator = [['' for _ in range(4)] for _ in range(4)]
for A in range(4):
    for B in range(4):
        if gamma_p[A][B] !='':
            spropagator[A][B] += gamma_p[A][B]
        if (A==B):
            if spropagator[A][B] != "":
                spropagator[A][B] += "+"
            spropagator[A][B] += "r*m"

# sum all terms of the propagator beside the pure fermionic part
lorentz_prop = [[[['' for _ in range(4)] for _ in range(4)] for _ in range(4)] for _ in range(4)]
for A in range(4):
    for mu in range(4):
        for B in range(4):
            for nu in range(4):
                if g_term[A][B][mu][nu] != 0.0:
                    lorentz_prop[A][B][mu][nu] += str(g_term[A][B][mu][nu])
                if p_term[A][B][mu][nu] != "":
                    if lorentz_prop[A][B][mu][nu] != "":
                        lorentz_prop[A][B][mu][nu] += "+"
                    lorentz_prop[A][B][mu][nu] += "4./(3.*m2)" + "*" + p_term[A][B][mu][nu]
                if gamma_prod[A][B][mu][nu] != complex(0.0):
                    if gamma_prod[A][B][mu][nu] == complex(1., 0.):
                        if lorentz_prop[A][B][mu][nu] != "":
                            lorentz_prop[A][B][mu][nu] += "+"
                        lorentz_prop[A][B][mu][nu] += "1./3."
                    elif gamma_prod[A][B][mu][nu] == complex(-1., 0.):
                        lorentz_prop[A][B][mu][nu] += "-1./3."
                    elif gamma_prod[A][B][mu][nu] == complex(0., 1.):
                        if lorentz_prop[A][B][mu][nu] != "":
                            lorentz_prop[A][B][mu][nu] += "+"
                        lorentz_prop[A][B][mu][nu] += imag_var + "/3."
                    elif gamma_prod[A][B][mu][nu] == complex(0., -1.):
                        lorentz_prop[A][B][mu][nu] += "-" + imag_var + "/3."
                    else:
                        if lorentz_prop[A][B][mu][nu] != "":
                            lorentz_prop[A][B][mu][nu] += "+"
                        lorentz_prop[A][B][mu][nu] += "1./3." + "*" + gamma_prod[A][B][mu][nu]
                if gammamu_pdagger[A][B][mu] != "":
                    lorentz_prop[A][B][mu][nu] += "-1./(3.*m2)" + "*" + "(" + gammamu_pdagger[A][B][mu] + ")" + "*"
                                                  + p[nu]
                if pdagger_gammanu[A][B][nu] != "":
                    lorentz_prop[A][B][mu][nu] += "-1./(3.*m2)" + "*" + "(" + pdagger_gammanu[A][B][nu] + ")" + "*"
                                                  + p[mu]

# complete denominator of the propagator
result = [[[['' for _ in range(4)] for _ in range(4)] for _ in range(4)] for _ in range(4)]
for A in range(4):
    for mu in range(4):
        for B in range(4):
            for nu in range(4):
                for C in range(4):
                    if spropagator[A][C] != "" and lorentz_prop[C][B][mu][nu] != "":
                        if result[A][B][mu][nu] != "":
                            result[A][B][mu][nu] += " + "
                        result[A][B][mu][nu] += "(" + spropagator[A][C] + ")" + "*"
                                                    + "(" + lorentz_prop[C][B][mu][nu] + ")"

# contraction of the propagator with an RS wave function from the right
result_bpos = [['' for _ in range(4)] for _ in range(4)]
# B>0 contraction
for A in range(4):
    for mu in range(4):
        for B in range(4):
            for nu in range(4):
                    if result[A][B][mu][nu] != "":
                        if result_bpos[A][mu] != "" and nu == 0:
                            result_bpos[A][mu] += " + "
                        if nu == 0:
                            result_bpos[A][mu] += "(" + result[A][B][mu][nu] + ")" + "*" + psi[B][nu]
                        else:
                            result_bpos[A][mu] += " - " + "(" + result[A][B][mu][nu] + ")" + "*" + psi[B][nu]

# contraction of the propagator with an RS wave function from the left
result_bneg = [['' for _ in range(4)] for _ in range(4)]
# B<0 contraction
for B in range(4):
    for nu in range(4):
        for A in range(4):
            for mu in range(4):
                    if result[A][B][mu][nu] != "":
                        if result_bneg[B][nu] != "" and mu == 0:
                            result_bneg[B][nu] += " + "
                        if mu == 0:
                            result_bneg[B][nu] += psi[A][mu] + "*" + "(" +  result[A][B][mu][nu] + ")"
                        else:
                            result_bneg[B][nu] += " - " + psi[A][mu] + "*" + "(" + result[A][B][mu][nu] + ")"
# end of script
*/