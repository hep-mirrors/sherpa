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
    double m_A;
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
    // TODO: Make a proper implementation of this function
    std::vector<std::vector<Complex> > SGetCurrent() const { return {}; }

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

  Settings& s = Settings::GetMainSettings();
  m_A = s["COMIX_RS_OFFSHELL_PARAMETER"].SetDefault(-1).Get<double>();

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
  if (hel) {
    if (p.PSpat() == 0) return CVec4Type(0., 1. / sqrt(2.), -std::complex<SType>(0., 1.) / sqrt(2.), 0., cr, ca, 0);
    else m_k_mod = ATOOLS::Vec4D(1., -ATOOLS::Vec3D(p) / p.PSpat());  // helicity basis
  }
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
// TODO: Do RS anti-particles need a complex-conjugate polarisation vector as used in HELAS convention (s. 1308.1668)?
//       (W- in Sherpa are also described by non-conjugate vector)
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
  if (!this->m_msv) THROW(fatal_error, "There is no massless Rarita-Schwinger-particle with Sz=1/2!")

  CRaritaSchwinger<SType> wf_p0(METOOLS::CRS<SType>::SpinorVectorProduct(CSpinor(r, abs(b), 1, p, cr, ca, hh, s, p.Abs2(),
                                                                            ms), EML(p, cr, ca), 1, cr, ca, s));
  CRaritaSchwinger<SType> wf_mp(METOOLS::CRS<SType>::SpinorVectorProduct(CSpinor(r, abs(b), -1, p, cr, ca, hh, s, p.Abs2(),
                                                                               ms), EMM(p, cr, ca), 1, cr, ca, s));
  CRaritaSchwinger<SType> wf = b>0?(sqrt(2.0/3.0) * wf_p0 - sqrt(1.0/3.0) * wf_mp) :
                                   (sqrt(2.0/3.0) * wf_p0 - sqrt(1.0/3.0) * wf_mp).Bar();
  return wf;
}

template<typename SType> CRaritaSchwinger<SType>
CRS<SType>::RSM(const ATOOLS::Vec4D &p, int r, int s, int b, int cr, int ca, int hh, int ms){
  if (!this->m_msv) THROW(fatal_error, "There is no massless Rarita-Schwinger-particle with Sz=1/2!")

  CRaritaSchwinger<SType> wf_pm(METOOLS::CRS<SType>::SpinorVectorProduct(CSpinor(r, abs(b), 1, p, cr, ca, hh, s, p.Abs2(),
                                                                                 ms), EMP(p, cr, ca), 1, cr, ca, s));
  CRaritaSchwinger<SType> wf_m0(METOOLS::CRS<SType>::SpinorVectorProduct(CSpinor(r, abs(b), -1, p, cr, ca, hh, s, p.Abs2(),
                                                                                 ms), EML(p, cr, ca), 1, cr, ca, s));
  //std::complex<SType> ephi = p.PPerp()==0.0?1.0:(p[1]-std::complex<SType>(0.,1.)*p[2]) * (1./p.PPerp());
  CRaritaSchwinger<SType> wf = b>0?-(sqrt(2.0/3.0) * wf_m0 + sqrt(1.0/3.0) * wf_pm) :
                               (-(sqrt(2.0/3.0) * wf_m0 + sqrt(1.0/3.0) * wf_pm)).Bar();
  return wf;
}

template<typename SType> CRaritaSchwinger<SType>
CRS<SType>::RSMM(const ATOOLS::Vec4D &p, const int r, const int s, const int b, const int cr, const int ca,
                 const int hh, const int ms) {
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
      // h=2 for bar vector-spinor
      j.SetH(anti^(this->m_dir>0)?3:2);
      //j.SetH(2);
      CRaScType *c(CRaScType::New(j));
      AddJ(c);
      DEBUG_VAR(j);
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
      j.SetH(anti^(this->m_dir>0)?1:0);
      //j.SetH(0);
#ifdef DEBUG__BG
      msg_Debugging()<<METHOD<<"(): "<<(this->m_dir>0?'I':'O')<<"++ "<<this->m_id
		   <<" "<<j<<" "<<(this->m_dir>0?this->m_fl.Bar():this->m_fl)
		   <<", m = "<<m_cmass<<" ("<<p.Mass()<<")\n";
      j.Test_Properties(p, this->m_dir);
#endif
      CRaScType *c(CRaScType::New(j));
      AddJ(c);
      DEBUG_VAR(j);
/*    if (p_sub) static_cast<Dipole_Color*>
      (p_sub->In().front()->Color().front())->AddJJK(c);*/
    }
  }
  if (ch<=0) {
    if (this->m_msv && (ch==0 || ch==-3)) {
      CRaScType j(anti^(this->m_dir>0)?
                  RS(p, this->m_fl.Majorana()?-2:-1, this->m_fl.Majorana()?(mode?3:2):3, 0, -this->m_dir, cr, ca):
                  RS(p, this->m_fl.Majorana()?2:1, this->m_fl.Majorana()?(mode?2:3):3, 0, this->m_dir, cr, ca));
      j.SetH(anti^(this->m_dir>0)?2:3);
      //j.SetH(3);
      CRaScType *c(CRaScType::New(j));
      AddJ(c);
      DEBUG_VAR(j);
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
      j.SetH(anti^(this->m_dir>0)?0:1);
//      j.SetH(1);
#ifdef DEBUG__BG
      msg_Debugging()<<METHOD<<"(): "<<(this->m_dir>0?'I':'O')<<"-- "<<this->m_id
		   <<" "<<j<<" "<<(this->m_dir>0?this->m_fl.Bar():this->m_fl)
		   <<", m = "<<m_cmass<<" ("<<p.Mass()<<")\n";
      j.Test_Properties(p, this->m_dir);
#endif
      CRaScType *c(CRaScType::New(j));
      AddJ(c);
      DEBUG_VAR(j);
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
  // TODO: Factor i required only for the A-independent term or for both?
  // prop: denominator of the propagator, only A-independent terms are multiply with it
  SComplex prop(1./(SType(this->m_p.Abs2())-m_cmass2));
  if (this->m_osd) prop=SComplex(1.);
  // all terms are multiplied with prop1
  SComplex prop1(M_I);
#ifdef DEBUG__BG
  msg_Debugging()<<"propagator: "<<prop
		 <<" <- p^2 = "<<this->m_p.Abs2()<<", m = "<<m_cmass<<"\n";
#endif

  // command in for propagator test, spin-3/2 mass should not be to small
  // command out -p lines for B>0
  // this->m_p[0] = sqrt(this->m_p.PSpat2()+m_cmass2.real());
  Vec4D p = this->m_p;
  Complex m = m_cmass;
  Complex m2 = m_cmass2;

  SComplex pp(Spinor<SType>::PPlus(this->m_p));
  SComplex pm(Spinor<SType>::PMinus(this->m_p));
  SComplex pt(Spinor<SType>::PT(this->m_p));
  SComplex ptc(Spinor<SType>::PTC(this->m_p));

  for (size_t i(0);i<m_j.size();++i) {
    CRaScType_Vector *j(m_j[i].template Get<CRaScType>());
    for (typename CRaScType_Vector::iterator
           jit(j->begin());jit!=j->end();++jit) {
      CRaScType j((*jit)->R(),(*jit)->B(),(**jit)(0),(**jit)(1),
                    (*jit)->H(),(*jit)->S(),
                    ((*jit)->On()&1)<<1|((*jit)->On()&2)>>1);
      double r = m_fl.MassSign();

      // components of the propagator numerator from Lecture Notes in Physics
      // 830, Stefan Scherer, Matthias R. Schindler: A Primer for Chiral
      // Perturbation Theory, Springer Verlag 2012,
      // DOI 10.1007/978-3-642-19254-8, eq. (4.188) contracted with a RaSc
      // wavefunction from the right (B>0) / left (B<0), calculated by a python
      // script (see comment at the end of this file)
      // TODO: Wirklich -p gemacht hier? see F_C!!!
      if ((*jit)->B()>0) {// S(-p)
        // for -p
        p *= -1.;
        pp *= -1.;
        pm *= -1.;
        pt *= -1.;
        ptc *= -1.;

        // vector propagator
        j[0]= ((r * m) * (-prop + (2. / (3. * m2)) * prop * p[0] * p[0])) * (**jit)[0] - ((r * m) * ((2. / (3. * m2)) * prop * p[0] * p[1])) * (**jit)[4] - ((r * m) * ((2. / (3. * m2)) * prop * p[0] * p[2])) * (**jit)[8] - ((r * m) * ((2. / (3. * m2)) * prop * p[0] * p[3])) * (**jit)[12] + (pm * (-prop + (2. / (3. * m2)) * prop * p[0] * p[0])) * (**jit)[2] - (pm * ((2. / (3. * m2)) * prop * p[0] * p[1])) * (**jit)[6] - (pm * ((2. / (3. * m2)) * prop * p[0] * p[2])) * (**jit)[10] - (pm * ((2. / (3. * m2)) * prop * p[0] * p[3])) * (**jit)[14] + ((-ptc) * (-prop + (2. / (3. * m2)) * prop * p[0] * p[0])) * (**jit)[3] - ((-ptc) * ((2. / (3. * m2)) * prop * p[0] * p[1])) * (**jit)[7] - ((-ptc) * ((2. / (3. * m2)) * prop * p[0] * p[2])) * (**jit)[11] - ((-ptc) * ((2. / (3. * m2)) * prop * p[0] * p[3])) * (**jit)[15];
        j[1]= ((r * m) * (-prop + (2. / (3. * m2)) * prop * p[0] * p[0])) * (**jit)[1] - ((r * m) * ((2. / (3. * m2)) * prop * p[0] * p[1])) * (**jit)[5] - ((r * m) * ((2. / (3. * m2)) * prop * p[0] * p[2])) * (**jit)[9] - ((r * m) * ((2. / (3. * m2)) * prop * p[0] * p[3])) * (**jit)[13] + ((-pt) * (-prop + (2. / (3. * m2)) * prop * p[0] * p[0])) * (**jit)[2] - ((-pt) * ((2. / (3. * m2)) * prop * p[0] * p[1])) * (**jit)[6] - ((-pt) * ((2. / (3. * m2)) * prop * p[0] * p[2])) * (**jit)[10] - ((-pt) * ((2. / (3. * m2)) * prop * p[0] * p[3])) * (**jit)[14] + (pp * (-prop + (2. / (3. * m2)) * prop * p[0] * p[0])) * (**jit)[3] - (pp * ((2. / (3. * m2)) * prop * p[0] * p[1])) * (**jit)[7] - (pp * ((2. / (3. * m2)) * prop * p[0] * p[2])) * (**jit)[11] - (pp * ((2. / (3. * m2)) * prop * p[0] * p[3])) * (**jit)[15];
        j[2]= (pp*(-prop+ (2./(3. * m2)) * prop * p[0] * p[0])) * (**jit)[0] - (pp * ((2. / (3. * m2)) * prop * p[0] * p[1])) * (**jit)[4] - (pp * ((2. / (3. * m2)) * prop * p[0] * p[2])) * (**jit)[8] - (pp * ((2. / (3. * m2)) * prop * p[0] * p[3])) * (**jit)[12] + (ptc * (-prop + (2. / (3. * m2)) * prop * p[0] * p[0])) * (**jit)[1] - (ptc * ((2. / (3. * m2)) * prop * p[0] * p[1])) * (**jit)[5] - (ptc * ((2. / (3. * m2)) * prop * p[0] * p[2])) * (**jit)[9] - (ptc * ((2. / (3. * m2)) * prop * p[0] * p[3])) * (**jit)[13] + ((r * m) * (-prop + (2. / (3. * m2)) * prop * p[0] * p[0])) * (**jit)[2] - ((r * m) * ((2. / (3. * m2)) * prop * p[0] * p[1])) * (**jit)[6] - ((r * m) * ((2. / (3. * m2)) * prop * p[0] * p[2])) * (**jit)[10] - ((r * m) * ((2. / (3. * m2)) * prop * p[0] * p[3])) * (**jit)[14];
        j[3]= (pt*(-prop+ (2./(3. * m2)) * prop * p[0] * p[0])) * (**jit)[0] - (pt * ((2. / (3. * m2)) * prop * p[0] * p[1])) * (**jit)[4] - (pt * ((2. / (3. * m2)) * prop * p[0] * p[2])) * (**jit)[8] - (pt * ((2. / (3. * m2)) * prop * p[0] * p[3])) * (**jit)[12] + (pm * (-prop + (2. / (3. * m2)) * prop * p[0] * p[0])) * (**jit)[1] - (pm * ((2. / (3. * m2)) * prop * p[0] * p[1])) * (**jit)[5] - (pm * ((2. / (3. * m2)) * prop * p[0] * p[2])) * (**jit)[9] - (pm * ((2. / (3. * m2)) * prop * p[0] * p[3])) * (**jit)[13] + ((r * m) * (-prop + (2. / (3. * m2)) * prop * p[0] * p[0])) * (**jit)[3] - ((r * m) * ((2. / (3. * m2)) * prop * p[0] * p[1])) * (**jit)[7] - ((r * m) * ((2. / (3. * m2)) * prop * p[0] * p[2])) * (**jit)[11] - ((r * m) * ((2. / (3. * m2)) * prop * p[0] * p[3])) * (**jit)[15];
        j[4]= ((r * m) * ((2. / (3. * m2)) * prop * p[1] * p[0])) * (**jit)[0] - ((r * m) * (prop + (2. / (3. * m2)) * prop * p[1] * p[1])) * (**jit)[4] - ((r * m) * ((2. / (3. * m2)) * prop * p[1] * p[2])) * (**jit)[8] - ((r * m) * ((2. / (3. * m2)) * prop * p[1] * p[3])) * (**jit)[12] + (pm * ((2. / (3. * m2)) * prop * p[1] * p[0])) * (**jit)[2] - (pm * (prop + (2. / (3. * m2)) * prop * p[1] * p[1])) * (**jit)[6] - (pm * ((2. / (3. * m2)) * prop * p[1] * p[2])) * (**jit)[10] - (pm * ((2. / (3. * m2)) * prop * p[1] * p[3])) * (**jit)[14] + ((-ptc) * ((2. / (3. * m2)) * prop * p[1] * p[0])) * (**jit)[3] - ((-ptc) * (prop + (2. / (3. * m2)) * prop * p[1] * p[1])) * (**jit)[7] - ((-ptc) * ((2. / (3. * m2)) * prop * p[1] * p[2])) * (**jit)[11] - ((-ptc) * ((2. / (3. * m2)) * prop * p[1] * p[3])) * (**jit)[15];
        j[5]= ((r * m) * ((2. / (3. * m2)) * prop * p[1] * p[0])) * (**jit)[1] - ((r * m) * (prop + (2. / (3. * m2)) * prop * p[1] * p[1])) * (**jit)[5] - ((r * m) * ((2. / (3. * m2)) * prop * p[1] * p[2])) * (**jit)[9] - ((r * m) * ((2. / (3. * m2)) * prop * p[1] * p[3])) * (**jit)[13] + ((-pt) * ((2. / (3. * m2)) * prop * p[1] * p[0])) * (**jit)[2] - ((-pt) * (prop + (2. / (3. * m2)) * prop * p[1] * p[1])) * (**jit)[6] - ((-pt) * ((2. / (3. * m2)) * prop * p[1] * p[2])) * (**jit)[10] - ((-pt) * ((2. / (3. * m2)) * prop * p[1] * p[3])) * (**jit)[14] + (pp * ((2. / (3. * m2)) * prop * p[1] * p[0])) * (**jit)[3] - (pp * (prop + (2. / (3. * m2)) * prop * p[1] * p[1])) * (**jit)[7] - (pp * ((2. / (3. * m2)) * prop * p[1] * p[2])) * (**jit)[11] - (pp * ((2. / (3. * m2)) * prop * p[1] * p[3])) * (**jit)[15];
        j[6]= (pp*((2./(3. * m2)) * prop * p[1] * p[0])) * (**jit)[0] - (pp * (prop + (2. / (3. * m2)) * prop * p[1] * p[1])) * (**jit)[4] - (pp * ((2. / (3. * m2)) * prop * p[1] * p[2])) * (**jit)[8] - (pp * ((2. / (3. * m2)) * prop * p[1] * p[3])) * (**jit)[12] + (ptc * ((2. / (3. * m2)) * prop * p[1] * p[0])) * (**jit)[1] - (ptc * (prop + (2. / (3. * m2)) * prop * p[1] * p[1])) * (**jit)[5] - (ptc * ((2. / (3. * m2)) * prop * p[1] * p[2])) * (**jit)[9] - (ptc * ((2. / (3. * m2)) * prop * p[1] * p[3])) * (**jit)[13] + ((r * m) * ((2. / (3. * m2)) * prop * p[1] * p[0])) * (**jit)[2] - ((r * m) * (prop + (2. / (3. * m2)) * prop * p[1] * p[1])) * (**jit)[6] - ((r * m) * ((2. / (3. * m2)) * prop * p[1] * p[2])) * (**jit)[10] - ((r * m) * ((2. / (3. * m2)) * prop * p[1] * p[3])) * (**jit)[14];
        j[7]= (pt*((2./(3. * m2)) * prop * p[1] * p[0])) * (**jit)[0] - (pt * (prop + (2. / (3. * m2)) * prop * p[1] * p[1])) * (**jit)[4] - (pt * ((2. / (3. * m2)) * prop * p[1] * p[2])) * (**jit)[8] - (pt * ((2. / (3. * m2)) * prop * p[1] * p[3])) * (**jit)[12] + (pm * ((2. / (3. * m2)) * prop * p[1] * p[0])) * (**jit)[1] - (pm * (prop + (2. / (3. * m2)) * prop * p[1] * p[1])) * (**jit)[5] - (pm * ((2. / (3. * m2)) * prop * p[1] * p[2])) * (**jit)[9] - (pm * ((2. / (3. * m2)) * prop * p[1] * p[3])) * (**jit)[13] + ((r * m) * ((2. / (3. * m2)) * prop * p[1] * p[0])) * (**jit)[3] - ((r * m) * (prop + (2. / (3. * m2)) * prop * p[1] * p[1])) * (**jit)[7] - ((r * m) * ((2. / (3. * m2)) * prop * p[1] * p[2])) * (**jit)[11] - ((r * m) * ((2. / (3. * m2)) * prop * p[1] * p[3])) * (**jit)[15];
        j[8]= ((r * m) * ((2. / (3. * m2)) * prop * p[2] * p[0])) * (**jit)[0] - ((r * m) * ((2. / (3. * m2)) * prop * p[2] * p[1])) * (**jit)[4] - ((r * m) * (prop + (2. / (3. * m2)) * prop * p[2] * p[2])) * (**jit)[8] - ((r * m) * ((2. / (3. * m2)) * prop * p[2] * p[3])) * (**jit)[12] + (pm * ((2. / (3. * m2)) * prop * p[2] * p[0])) * (**jit)[2] - (pm * ((2. / (3. * m2)) * prop * p[2] * p[1])) * (**jit)[6] - (pm * (prop + (2. / (3. * m2)) * prop * p[2] * p[2])) * (**jit)[10] - (pm * ((2. / (3. * m2)) * prop * p[2] * p[3])) * (**jit)[14] + ((-ptc) * ((2. / (3. * m2)) * prop * p[2] * p[0])) * (**jit)[3] - ((-ptc) * ((2. / (3. * m2)) * prop * p[2] * p[1])) * (**jit)[7] - ((-ptc) * (prop + (2. / (3. * m2)) * prop * p[2] * p[2])) * (**jit)[11] - ((-ptc) * ((2. / (3. * m2)) * prop * p[2] * p[3])) * (**jit)[15];
        j[9]= ((r * m) * ((2. / (3. * m2)) * prop * p[2] * p[0])) * (**jit)[1] - ((r * m) * ((2. / (3. * m2)) * prop * p[2] * p[1])) * (**jit)[5] - ((r * m) * (prop + (2. / (3. * m2)) * prop * p[2] * p[2])) * (**jit)[9] - ((r * m) * ((2. / (3. * m2)) * prop * p[2] * p[3])) * (**jit)[13] + ((-pt) * ((2. / (3. * m2)) * prop * p[2] * p[0])) * (**jit)[2] - ((-pt) * ((2. / (3. * m2)) * prop * p[2] * p[1])) * (**jit)[6] - ((-pt) * (prop + (2. / (3. * m2)) * prop * p[2] * p[2])) * (**jit)[10] - ((-pt) * ((2. / (3. * m2)) * prop * p[2] * p[3])) * (**jit)[14] + (pp * ((2. / (3. * m2)) * prop * p[2] * p[0])) * (**jit)[3] - (pp * ((2. / (3. * m2)) * prop * p[2] * p[1])) * (**jit)[7] - (pp * (prop + (2. / (3. * m2)) * prop * p[2] * p[2])) * (**jit)[11] - (pp * ((2. / (3. * m2)) * prop * p[2] * p[3])) * (**jit)[15];
        j[10]= (pp*((2./(3. * m2)) * prop * p[2] * p[0])) * (**jit)[0] - (pp * ((2. / (3. * m2)) * prop * p[2] * p[1])) * (**jit)[4] - (pp * (prop + (2. / (3. * m2)) * prop * p[2] * p[2])) * (**jit)[8] - (pp * ((2. / (3. * m2)) * prop * p[2] * p[3])) * (**jit)[12] + (ptc * ((2. / (3. * m2)) * prop * p[2] * p[0])) * (**jit)[1] - (ptc * ((2. / (3. * m2)) * prop * p[2] * p[1])) * (**jit)[5] - (ptc * (prop + (2. / (3. * m2)) * prop * p[2] * p[2])) * (**jit)[9] - (ptc * ((2. / (3. * m2)) * prop * p[2] * p[3])) * (**jit)[13] + ((r * m) * ((2. / (3. * m2)) * prop * p[2] * p[0])) * (**jit)[2] - ((r * m) * ((2. / (3. * m2)) * prop * p[2] * p[1])) * (**jit)[6] - ((r * m) * (prop + (2. / (3. * m2)) * prop * p[2] * p[2])) * (**jit)[10] - ((r * m) * ((2. / (3. * m2)) * prop * p[2] * p[3])) * (**jit)[14];
        j[11]= (pt*((2./(3. * m2)) * prop * p[2] * p[0])) * (**jit)[0] - (pt * ((2. / (3. * m2)) * prop * p[2] * p[1])) * (**jit)[4] - (pt * (prop + (2. / (3. * m2)) * prop * p[2] * p[2])) * (**jit)[8] - (pt * ((2. / (3. * m2)) * prop * p[2] * p[3])) * (**jit)[12] + (pm * ((2. / (3. * m2)) * prop * p[2] * p[0])) * (**jit)[1] - (pm * ((2. / (3. * m2)) * prop * p[2] * p[1])) * (**jit)[5] - (pm * (prop + (2. / (3. * m2)) * prop * p[2] * p[2])) * (**jit)[9] - (pm * ((2. / (3. * m2)) * prop * p[2] * p[3])) * (**jit)[13] + ((r * m) * ((2. / (3. * m2)) * prop * p[2] * p[0])) * (**jit)[3] - ((r * m) * ((2. / (3. * m2)) * prop * p[2] * p[1])) * (**jit)[7] - ((r * m) * (prop + (2. / (3. * m2)) * prop * p[2] * p[2])) * (**jit)[11] - ((r * m) * ((2. / (3. * m2)) * prop * p[2] * p[3])) * (**jit)[15];
        j[12]= ((r * m) * ((2. / (3. * m2)) * prop * p[3] * p[0])) * (**jit)[0] - ((r * m) * ((2. / (3. * m2)) * prop * p[3] * p[1])) * (**jit)[4] - ((r * m) * ((2. / (3. * m2)) * prop * p[3] * p[2])) * (**jit)[8] - ((r * m) * (prop + (2. / (3. * m2)) * prop * p[3] * p[3])) * (**jit)[12] + (pm * ((2. / (3. * m2)) * prop * p[3] * p[0])) * (**jit)[2] - (pm * ((2. / (3. * m2)) * prop * p[3] * p[1])) * (**jit)[6] - (pm * ((2. / (3. * m2)) * prop * p[3] * p[2])) * (**jit)[10] - (pm * (prop + (2. / (3. * m2)) * prop * p[3] * p[3])) * (**jit)[14] + ((-ptc) * ((2. / (3. * m2)) * prop * p[3] * p[0])) * (**jit)[3] - ((-ptc) * ((2. / (3. * m2)) * prop * p[3] * p[1])) * (**jit)[7] - ((-ptc) * ((2. / (3. * m2)) * prop * p[3] * p[2])) * (**jit)[11] - ((-ptc) * (prop + (2. / (3. * m2)) * prop * p[3] * p[3])) * (**jit)[15];
        j[13]= ((r * m) * ((2. / (3. * m2)) * prop * p[3] * p[0])) * (**jit)[1] - ((r * m) * ((2. / (3. * m2)) * prop * p[3] * p[1])) * (**jit)[5] - ((r * m) * ((2. / (3. * m2)) * prop * p[3] * p[2])) * (**jit)[9] - ((r * m) * (prop + (2. / (3. * m2)) * prop * p[3] * p[3])) * (**jit)[13] + ((-pt) * ((2. / (3. * m2)) * prop * p[3] * p[0])) * (**jit)[2] - ((-pt) * ((2. / (3. * m2)) * prop * p[3] * p[1])) * (**jit)[6] - ((-pt) * ((2. / (3. * m2)) * prop * p[3] * p[2])) * (**jit)[10] - ((-pt) * (prop + (2. / (3. * m2)) * prop * p[3] * p[3])) * (**jit)[14] + (pp * ((2. / (3. * m2)) * prop * p[3] * p[0])) * (**jit)[3] - (pp * ((2. / (3. * m2)) * prop * p[3] * p[1])) * (**jit)[7] - (pp * ((2. / (3. * m2)) * prop * p[3] * p[2])) * (**jit)[11] - (pp * (prop + (2. / (3. * m2)) * prop * p[3] * p[3])) * (**jit)[15];
        j[14]= (pp*((2./(3. * m2)) * prop * p[3] * p[0])) * (**jit)[0] - (pp * ((2. / (3. * m2)) * prop * p[3] * p[1])) * (**jit)[4] - (pp * ((2. / (3. * m2)) * prop * p[3] * p[2])) * (**jit)[8] - (pp * (prop + (2. / (3. * m2)) * prop * p[3] * p[3])) * (**jit)[12] + (ptc * ((2. / (3. * m2)) * prop * p[3] * p[0])) * (**jit)[1] - (ptc * ((2. / (3. * m2)) * prop * p[3] * p[1])) * (**jit)[5] - (ptc * ((2. / (3. * m2)) * prop * p[3] * p[2])) * (**jit)[9] - (ptc * (prop + (2. / (3. * m2)) * prop * p[3] * p[3])) * (**jit)[13] + ((r * m) * ((2. / (3. * m2)) * prop * p[3] * p[0])) * (**jit)[2] - ((r * m) * ((2. / (3. * m2)) * prop * p[3] * p[1])) * (**jit)[6] - ((r * m) * ((2. / (3. * m2)) * prop * p[3] * p[2])) * (**jit)[10] - ((r * m) * (prop + (2. / (3. * m2)) * prop * p[3] * p[3])) * (**jit)[14];
        j[15]= (pt*((2./(3. * m2)) * prop * p[3] * p[0])) * (**jit)[0] - (pt * ((2. / (3. * m2)) * prop * p[3] * p[1])) * (**jit)[4] - (pt * ((2. / (3. * m2)) * prop * p[3] * p[2])) * (**jit)[8] - (pt * (prop + (2. / (3. * m2)) * prop * p[3] * p[3])) * (**jit)[12] + (pm * ((2. / (3. * m2)) * prop * p[3] * p[0])) * (**jit)[1] - (pm * ((2. / (3. * m2)) * prop * p[3] * p[1])) * (**jit)[5] - (pm * ((2. / (3. * m2)) * prop * p[3] * p[2])) * (**jit)[9] - (pm * (prop + (2. / (3. * m2)) * prop * p[3] * p[3])) * (**jit)[13] + ((r * m) * ((2. / (3. * m2)) * prop * p[3] * p[0])) * (**jit)[3] - ((r * m) * ((2. / (3. * m2)) * prop * p[3] * p[1])) * (**jit)[7] - ((r * m) * ((2. / (3. * m2)) * prop * p[3] * p[2])) * (**jit)[11] - ((r * m) * (prop + (2. / (3. * m2)) * prop * p[3] * p[3])) * (**jit)[15];

        // gammamu gammanu
        j[0]+= ((r * m) * ((1. / 3.) * prop)) * (**jit)[0] - ((r * m) * (-(1. / 3.) * prop)) * (**jit)[4 * SpinorType::R3()] - ((r * m) * (-(1. / 3.) * prop)) * (**jit)[4 * SpinorType::R1() + 1] - ((r * m) * ((M_I / 3.) * prop)) * (**jit)[4 * SpinorType::R2() + 1] + (pm * ((1. / 3.) * prop)) * (**jit)[2] - ((-ptc) * ((1. / 3.) * prop)) * (**jit)[4 * SpinorType::R1() + 2] - ((-ptc) * ((M_I / 3.) * prop)) * (**jit)[4 * SpinorType::R2() + 2] - (pm * ((1. / 3.) * prop)) * (**jit)[4 * SpinorType::R3() + 2] + ((-ptc) * ((1. / 3.) * prop)) * (**jit)[3] - (pm * ((1. / 3.) * prop)) * (**jit)[4 * SpinorType::R1() + 3] - (pm * (-(M_I / 3.) * prop)) * (**jit)[4 * SpinorType::R2() + 3] - ((-ptc) * (-(1. / 3.) * prop)) * (**jit)[4 * SpinorType::R3() + 3];
        j[1]+= - ((r * m) * (-(1. / 3.) * prop)) * (**jit)[4 * SpinorType::R1()] - ((r * m) * (-(M_I / 3.) * prop)) * (**jit)[4 * SpinorType::R2()] + ((r * m) * ((1. / 3.) * prop)) * (**jit)[1] - ((r * m) * ((1. / 3.) * prop)) * (**jit)[4 * SpinorType::R3() + 1] + ((-pt) * ((1. / 3.) * prop)) * (**jit)[2] - (pp * ((1. / 3.) * prop)) * (**jit)[4 * SpinorType::R1() + 2] - (pp * ((M_I / 3.) * prop)) * (**jit)[4 * SpinorType::R2() + 2] - ((-pt) * ((1. / 3.) * prop)) * (**jit)[4 * SpinorType::R3() + 2] + (pp * ((1. / 3.) * prop)) * (**jit)[3] - ((-pt) * ((1. / 3.) * prop)) * (**jit)[4 * SpinorType::R1() + 3] - ((-pt) * (-(M_I / 3.) * prop)) * (**jit)[4 * SpinorType::R2() + 3] - (pp * (-(1. / 3.) * prop)) * (**jit)[4 * SpinorType::R3() + 3];
        j[2]+=(pp*((1./3.)*prop))*(**jit)[0] - (ptc*(-(1./3.)*prop))*(**jit)[4*SpinorType::R1()] - (ptc*(-(M_I/3.)*prop))*(**jit)[4*SpinorType::R2()] - (pp*(-(1./3.)*prop))*(**jit)[4*SpinorType::R3()] + (ptc*((1./3.)*prop))*(**jit)[1] - (pp*(-(1./3.)*prop))*(**jit)[4*SpinorType::R1()+1] - (pp*((M_I/3.)*prop))*(**jit)[4*SpinorType::R2()+1] - (ptc*((1./3.)*prop))*(**jit)[4*SpinorType::R3()+1] + ((r * m) * ((1. / 3.) * prop)) * (**jit)[2] - ((r * m) * ((1. / 3.) * prop)) * (**jit)[4 * SpinorType::R3() + 2] - ((r * m) * ((1. / 3.) * prop)) * (**jit)[4 * SpinorType::R1() + 3] - ((r * m) * (-(M_I / 3.) * prop)) * (**jit)[4 * SpinorType::R2() + 3];
        j[3]+=(pt*((1./3.)*prop))*(**jit)[0] - (pm*(-(1./3.)*prop))*(**jit)[4*SpinorType::R1()] - (pm*(-(M_I/3.)*prop))*(**jit)[4*SpinorType::R2()] - (pt*(-(1./3.)*prop))*(**jit)[4*SpinorType::R3()] + (pm*((1./3.)*prop))*(**jit)[1] - (pt*(-(1./3.)*prop))*(**jit)[4*SpinorType::R1()+1] - (pt*((M_I/3.)*prop))*(**jit)[4*SpinorType::R2()+1] - (pm*((1./3.)*prop))*(**jit)[4*SpinorType::R3()+1] - ((r * m) * ((1. / 3.) * prop)) * (**jit)[4 * SpinorType::R1() + 2] - ((r * m) * ((M_I / 3.) * prop)) * (**jit)[4 * SpinorType::R2() + 2] + ((r * m) * ((1. / 3.) * prop)) * (**jit)[3] - ((r * m) * (-(1. / 3.) * prop)) * (**jit)[4 * SpinorType::R3() + 3];
        j[4*SpinorType::R1()]+= - ((r * m) * (-(1. / 3.) * prop)) * (**jit)[4 * SpinorType::R1()] - ((r * m) * (-(M_I / 3.) * prop)) * (**jit)[4 * SpinorType::R2()] + ((r * m) * ((1. / 3.) * prop)) * (**jit)[1] - ((r * m) * ((1. / 3.) * prop)) * (**jit)[4 * SpinorType::R3() + 1] + ((-ptc) * (-(1. / 3.) * prop)) * (**jit)[2] - (pm * (-(1. / 3.) * prop)) * (**jit)[4 * SpinorType::R1() + 2] - (pm * (-(M_I / 3.) * prop)) * (**jit)[4 * SpinorType::R2() + 2] - ((-ptc) * (-(1. / 3.) * prop)) * (**jit)[4 * SpinorType::R3() + 2] + (pm * (-(1. / 3.) * prop)) * (**jit)[3] - ((-ptc) * (-(1. / 3.) * prop)) * (**jit)[4 * SpinorType::R1() + 3] - ((-ptc) * ((M_I / 3.) * prop)) * (**jit)[4 * SpinorType::R2() + 3] - (pm * ((1. / 3.) * prop)) * (**jit)[4 * SpinorType::R3() + 3];
        j[4*SpinorType::R1()+1]+= ((r * m) * ((1. / 3.) * prop)) * (**jit)[0] - ((r * m) * (-(1. / 3.) * prop)) * (**jit)[4 * SpinorType::R3()] - ((r * m) * (-(1. / 3.) * prop)) * (**jit)[4 * SpinorType::R1() + 1] - ((r * m) * ((M_I / 3.) * prop)) * (**jit)[4 * SpinorType::R2() + 1] + (pp * (-(1. / 3.) * prop)) * (**jit)[2] - ((-pt) * (-(1. / 3.) * prop)) * (**jit)[4 * SpinorType::R1() + 2] - ((-pt) * (-(M_I / 3.) * prop)) * (**jit)[4 * SpinorType::R2() + 2] - (pp * (-(1. / 3.) * prop)) * (**jit)[4 * SpinorType::R3() + 2] + ((-pt) * (-(1. / 3.) * prop)) * (**jit)[3] - (pp * (-(1. / 3.) * prop)) * (**jit)[4 * SpinorType::R1() + 3] - (pp * ((M_I / 3.) * prop)) * (**jit)[4 * SpinorType::R2() + 3] - ((-pt) * ((1. / 3.) * prop)) * (**jit)[4 * SpinorType::R3() + 3];
        j[4*SpinorType::R1()+2]+=(ptc*((1./3.)*prop))*(**jit)[0] - (pp*(-(1./3.)*prop))*(**jit)[4*SpinorType::R1()] - (pp*(-(M_I/3.)*prop))*(**jit)[4*SpinorType::R2()] - (ptc*(-(1./3.)*prop))*(**jit)[4*SpinorType::R3()] + (pp*((1./3.)*prop))*(**jit)[1] - (ptc*(-(1./3.)*prop))*(**jit)[4*SpinorType::R1()+1] - (ptc*((M_I/3.)*prop))*(**jit)[4*SpinorType::R2()+1] - (pp*((1./3.)*prop))*(**jit)[4*SpinorType::R3()+1] - ((r * m) * (-(1. / 3.) * prop)) * (**jit)[4 * SpinorType::R1() + 2] - ((r * m) * (-(M_I / 3.) * prop)) * (**jit)[4 * SpinorType::R2() + 2] + ((r * m) * (-(1. / 3.) * prop)) * (**jit)[3] - ((r * m) * ((1. / 3.) * prop)) * (**jit)[4 * SpinorType::R3() + 3];
        j[4*SpinorType::R1()+3]+=(pm*((1./3.)*prop))*(**jit)[0] - (pt*(-(1./3.)*prop))*(**jit)[4*SpinorType::R1()] - (pt*(-(M_I/3.)*prop))*(**jit)[4*SpinorType::R2()] - (pm*(-(1./3.)*prop))*(**jit)[4*SpinorType::R3()] + (pt*((1./3.)*prop))*(**jit)[1] - (pm*(-(1./3.)*prop))*(**jit)[4*SpinorType::R1()+1] - (pm*((M_I/3.)*prop))*(**jit)[4*SpinorType::R2()+1] - (pt*((1./3.)*prop))*(**jit)[4*SpinorType::R3()+1] + ((r * m) * (-(1. / 3.) * prop)) * (**jit)[2] - ((r * m) * (-(1. / 3.) * prop)) * (**jit)[4 * SpinorType::R3() + 2] - ((r * m) * (-(1. / 3.) * prop)) * (**jit)[4 * SpinorType::R1() + 3] - ((r * m) * ((M_I / 3.) * prop)) * (**jit)[4 * SpinorType::R2() + 3];
        j[4*SpinorType::R2()]+= - ((r * m) * ((M_I / 3.) * prop)) * (**jit)[4 * SpinorType::R1()] - ((r * m) * (-(1. / 3.) * prop)) * (**jit)[4 * SpinorType::R2()] + ((r * m) * (-(M_I / 3.) * prop)) * (**jit)[1] - ((r * m) * (-(M_I / 3.) * prop)) * (**jit)[4 * SpinorType::R3() + 1] + ((-ptc) * (-(M_I / 3.) * prop)) * (**jit)[2] - (pm * ((M_I / 3.) * prop)) * (**jit)[4 * SpinorType::R1() + 2] - (pm * (-(1. / 3.) * prop)) * (**jit)[4 * SpinorType::R2() + 2] - ((-ptc) * (-(M_I / 3.) * prop)) * (**jit)[4 * SpinorType::R3() + 2] + (pm * ((M_I / 3.) * prop)) * (**jit)[3] - ((-ptc) * (-(M_I / 3.) * prop)) * (**jit)[4 * SpinorType::R1() + 3] - ((-ptc) * (-(1. / 3.) * prop)) * (**jit)[4 * SpinorType::R2() + 3] - (pm * (-(M_I / 3.) * prop)) * (**jit)[4 * SpinorType::R3() + 3];
        j[4*SpinorType::R2()+1]+= ((r * m) * ((M_I / 3.) * prop)) * (**jit)[0] - ((r * m) * (-(M_I / 3.) * prop)) * (**jit)[4 * SpinorType::R3()] - ((r * m) * (-(M_I / 3.) * prop)) * (**jit)[4 * SpinorType::R1() + 1] - ((r * m) * (-(1. / 3.) * prop)) * (**jit)[4 * SpinorType::R2() + 1] + (pp * (-(M_I / 3.) * prop)) * (**jit)[2] - ((-pt) * ((M_I / 3.) * prop)) * (**jit)[4 * SpinorType::R1() + 2] - ((-pt) * (-(1. / 3.) * prop)) * (**jit)[4 * SpinorType::R2() + 2] - (pp * (-(M_I / 3.) * prop)) * (**jit)[4 * SpinorType::R3() + 2] + ((-pt) * ((M_I / 3.) * prop)) * (**jit)[3] - (pp * (-(M_I / 3.) * prop)) * (**jit)[4 * SpinorType::R1() + 3] - (pp * (-(1. / 3.) * prop)) * (**jit)[4 * SpinorType::R2() + 3] - ((-pt) * (-(M_I / 3.) * prop)) * (**jit)[4 * SpinorType::R3() + 3];
        j[4*SpinorType::R2()+2]+=(ptc*((M_I/3.)*prop))*(**jit)[0] - (pp*((M_I/3.)*prop))*(**jit)[4*SpinorType::R1()] - (pp*(-(1./3.)*prop))*(**jit)[4*SpinorType::R2()] - (ptc*(-(M_I/3.)*prop))*(**jit)[4*SpinorType::R3()] + (pp*(-(M_I/3.)*prop))*(**jit)[1] - (ptc*(-(M_I/3.)*prop))*(**jit)[4*SpinorType::R1()+1] - (ptc*(-(1./3.)*prop))*(**jit)[4*SpinorType::R2()+1] - (pp*(-(M_I/3.)*prop))*(**jit)[4*SpinorType::R3()+1] - ((r * m) * ((M_I / 3.) * prop)) * (**jit)[4 * SpinorType::R1() + 2] - ((r * m) * (-(1. / 3.) * prop)) * (**jit)[4 * SpinorType::R2() + 2] + ((r * m) * ((M_I / 3.) * prop)) * (**jit)[3] - ((r * m) * (-(M_I / 3.) * prop)) * (**jit)[4 * SpinorType::R3() + 3];
        j[4*SpinorType::R2()+3]+=(pm*((M_I/3.)*prop))*(**jit)[0] - (pt*((M_I/3.)*prop))*(**jit)[4*SpinorType::R1()] - (pt*(-(1./3.)*prop))*(**jit)[4*SpinorType::R2()] - (pm*(-(M_I/3.)*prop))*(**jit)[4*SpinorType::R3()] + (pt*(-(M_I/3.)*prop))*(**jit)[1] - (pm*(-(M_I/3.)*prop))*(**jit)[4*SpinorType::R1()+1] - (pm*(-(1./3.)*prop))*(**jit)[4*SpinorType::R2()+1] - (pt*(-(M_I/3.)*prop))*(**jit)[4*SpinorType::R3()+1] + ((r * m) * (-(M_I / 3.) * prop)) * (**jit)[2] - ((r * m) * (-(M_I / 3.) * prop)) * (**jit)[4 * SpinorType::R3() + 2] - ((r * m) * (-(M_I / 3.) * prop)) * (**jit)[4 * SpinorType::R1() + 3] - ((r * m) * (-(1. / 3.) * prop)) * (**jit)[4 * SpinorType::R2() + 3];
        j[4*SpinorType::R3()]+= ((r * m) * ((1. / 3.) * prop)) * (**jit)[0] - ((r * m) * (-(1. / 3.) * prop)) * (**jit)[4 * SpinorType::R3()] - ((r * m) * (-(1. / 3.) * prop)) * (**jit)[4 * SpinorType::R1() + 1] - ((r * m) * ((M_I / 3.) * prop)) * (**jit)[4 * SpinorType::R2() + 1] + (pm * (-(1. / 3.) * prop)) * (**jit)[2] - ((-ptc) * ((1. / 3.) * prop)) * (**jit)[4 * SpinorType::R1() + 2] - ((-ptc) * ((M_I / 3.) * prop)) * (**jit)[4 * SpinorType::R2() + 2] - (pm * (-(1. / 3.) * prop)) * (**jit)[4 * SpinorType::R3() + 2] + ((-ptc) * ((1. / 3.) * prop)) * (**jit)[3] - (pm * (-(1. / 3.) * prop)) * (**jit)[4 * SpinorType::R1() + 3] - (pm * ((M_I / 3.) * prop)) * (**jit)[4 * SpinorType::R2() + 3] - ((-ptc) * (-(1. / 3.) * prop)) * (**jit)[4 * SpinorType::R3() + 3];
        j[4*SpinorType::R3()+1]+= - ((r * m) * ((1. / 3.) * prop)) * (**jit)[4 * SpinorType::R1()] - ((r * m) * ((M_I / 3.) * prop)) * (**jit)[4 * SpinorType::R2()] + ((r * m) * (-(1. / 3.) * prop)) * (**jit)[1] - ((r * m) * (-(1. / 3.) * prop)) * (**jit)[4 * SpinorType::R3() + 1] + ((-pt) * (-(1. / 3.) * prop)) * (**jit)[2] - (pp * ((1. / 3.) * prop)) * (**jit)[4 * SpinorType::R1() + 2] - (pp * ((M_I / 3.) * prop)) * (**jit)[4 * SpinorType::R2() + 2] - ((-pt) * (-(1. / 3.) * prop)) * (**jit)[4 * SpinorType::R3() + 2] + (pp * ((1. / 3.) * prop)) * (**jit)[3] - ((-pt) * (-(1. / 3.) * prop)) * (**jit)[4 * SpinorType::R1() + 3] - ((-pt) * ((M_I / 3.) * prop)) * (**jit)[4 * SpinorType::R2() + 3] - (pp * (-(1. / 3.) * prop)) * (**jit)[4 * SpinorType::R3() + 3];
        j[4*SpinorType::R3()+2]+=(pp*((1./3.)*prop))*(**jit)[0] - (ptc*((1./3.)*prop))*(**jit)[4*SpinorType::R1()] - (ptc*((M_I/3.)*prop))*(**jit)[4*SpinorType::R2()] - (pp*(-(1./3.)*prop))*(**jit)[4*SpinorType::R3()] + (ptc*(-(1./3.)*prop))*(**jit)[1] - (pp*(-(1./3.)*prop))*(**jit)[4*SpinorType::R1()+1] - (pp*((M_I/3.)*prop))*(**jit)[4*SpinorType::R2()+1] - (ptc*(-(1./3.)*prop))*(**jit)[4*SpinorType::R3()+1] + ((r * m) * (-(1. / 3.) * prop)) * (**jit)[2] - ((r * m) * (-(1. / 3.) * prop)) * (**jit)[4 * SpinorType::R3() + 2] - ((r * m) * (-(1. / 3.) * prop)) * (**jit)[4 * SpinorType::R1() + 3] - ((r * m) * ((M_I / 3.) * prop)) * (**jit)[4 * SpinorType::R2() + 3];
        j[4*SpinorType::R3()+3]+=(pt*((1./3.)*prop))*(**jit)[0] - (pm*((1./3.)*prop))*(**jit)[4*SpinorType::R1()] - (pm*((M_I/3.)*prop))*(**jit)[4*SpinorType::R2()] - (pt*(-(1./3.)*prop))*(**jit)[4*SpinorType::R3()] + (pm*(-(1./3.)*prop))*(**jit)[1] - (pt*(-(1./3.)*prop))*(**jit)[4*SpinorType::R1()+1] - (pt*((M_I/3.)*prop))*(**jit)[4*SpinorType::R2()+1] - (pm*(-(1./3.)*prop))*(**jit)[4*SpinorType::R3()+1] - ((r * m) * ((1. / 3.) * prop)) * (**jit)[4 * SpinorType::R1() + 2] - ((r * m) * ((M_I / 3.) * prop)) * (**jit)[4 * SpinorType::R2() + 2] + ((r * m) * ((1. / 3.) * prop)) * (**jit)[3] - ((r * m) * (-(1. / 3.) * prop)) * (**jit)[4 * SpinorType::R3() + 3];

        // pmu gammanu
        j[0]+= (pm*(-(1./(3. * m)) * prop * p[0])) * (**jit)[0] - ((-ptc) * (-(1. / (3. * m)) * prop * -p[0])) * (**jit)[4 * SpinorType::R1()] - ((-ptc) * (-(1. / (3. * m)) * prop * -M_I * p[0])) * (**jit)[4 * SpinorType::R2()] - (pm * (-(1. / (3. * m)) * prop * -p[0])) * (**jit)[4 * SpinorType::R3()] + ((-ptc) * (-(1. / (3. * m)) * prop * p[0])) * (**jit)[1] - (pm * (-(1. / (3. * m)) * prop * -p[0])) * (**jit)[4 * SpinorType::R1() + 1] - (pm * (-(1. / (3. * m)) * prop * M_I * p[0])) * (**jit)[4 * SpinorType::R2() + 1] - ((-ptc) * (-(1. / (3. * m)) * prop * p[0])) * (**jit)[4 * SpinorType::R3() + 1] + ((r * m) * (-(1. / (3. * m)) * prop * p[0])) * (**jit)[2] - ((r * m) * (-(1. / (3. * m)) * prop * p[0])) * (**jit)[4 * SpinorType::R3() + 2] - ((r * m) * (-(1. / (3. * m)) * prop * p[0])) * (**jit)[4 * SpinorType::R1() + 3] - ((r * m) * (-(1. / (3. * m)) * prop * -M_I * p[0])) * (**jit)[4 * SpinorType::R2() + 3];
        j[1]+= ((-pt)*(-(1./(3. * m)) * prop * p[0])) * (**jit)[0] - (pp * (-(1. / (3. * m)) * prop * -p[0])) * (**jit)[4 * SpinorType::R1()] - (pp * (-(1. / (3. * m)) * prop * -M_I * p[0])) * (**jit)[4 * SpinorType::R2()] - ((-pt) * (-(1. / (3. * m)) * prop * -p[0])) * (**jit)[4 * SpinorType::R3()] + (pp * (-(1. / (3. * m)) * prop * p[0])) * (**jit)[1] - ((-pt) * (-(1. / (3. * m)) * prop * -p[0])) * (**jit)[4 * SpinorType::R1() + 1] - ((-pt) * (-(1. / (3. * m)) * prop * M_I * p[0])) * (**jit)[4 * SpinorType::R2() + 1] - (pp * (-(1. / (3. * m)) * prop * p[0])) * (**jit)[4 * SpinorType::R3() + 1] - ((r * m) * (-(1. / (3. * m)) * prop * p[0])) * (**jit)[4 * SpinorType::R1() + 2] - ((r * m) * (-(1. / (3. * m)) * prop * M_I * p[0])) * (**jit)[4 * SpinorType::R2() + 2] + ((r * m) * (-(1. / (3. * m)) * prop * p[0])) * (**jit)[3] - ((r * m) * (-(1. / (3. * m)) * prop * -p[0])) * (**jit)[4 * SpinorType::R3() + 3];
        j[2]+= ((r * m) * (-(1. / (3. * m)) * prop * p[0])) * (**jit)[0] - ((r * m) * (-(1. / (3. * m)) * prop * -p[0])) * (**jit)[4 * SpinorType::R3()] - ((r * m) * (-(1. / (3. * m)) * prop * -p[0])) * (**jit)[4 * SpinorType::R1() + 1] - ((r * m) * (-(1. / (3. * m)) * prop * M_I * p[0])) * (**jit)[4 * SpinorType::R2() + 1] + (pp * (-(1. / (3. * m)) * prop * p[0])) * (**jit)[2] - (ptc * (-(1. / (3. * m)) * prop * p[0])) * (**jit)[4 * SpinorType::R1() + 2] - (ptc * (-(1. / (3. * m)) * prop * M_I * p[0])) * (**jit)[4 * SpinorType::R2() + 2] - (pp * (-(1. / (3. * m)) * prop * p[0])) * (**jit)[4 * SpinorType::R3() + 2] + (ptc * (-(1. / (3. * m)) * prop * p[0])) * (**jit)[3] - (pp * (-(1. / (3. * m)) * prop * p[0])) * (**jit)[4 * SpinorType::R1() + 3] - (pp * (-(1. / (3. * m)) * prop * -M_I * p[0])) * (**jit)[4 * SpinorType::R2() + 3] - (ptc * (-(1. / (3. * m)) * prop * -p[0])) * (**jit)[4 * SpinorType::R3() + 3];
        j[3]+= - ((r * m) * (-(1. / (3. * m)) * prop * -p[0])) * (**jit)[4 * SpinorType::R1()] - ((r * m) * (-(1. / (3. * m)) * prop * -M_I * p[0])) * (**jit)[4 * SpinorType::R2()] + ((r * m) * (-(1. / (3. * m)) * prop * p[0])) * (**jit)[1] - ((r * m) * (-(1. / (3. * m)) * prop * p[0])) * (**jit)[4 * SpinorType::R3() + 1] + (pt * (-(1. / (3. * m)) * prop * p[0])) * (**jit)[2] - (pm * (-(1. / (3. * m)) * prop * p[0])) * (**jit)[4 * SpinorType::R1() + 2] - (pm * (-(1. / (3. * m)) * prop * M_I * p[0])) * (**jit)[4 * SpinorType::R2() + 2] - (pt * (-(1. / (3. * m)) * prop * p[0])) * (**jit)[4 * SpinorType::R3() + 2] + (pm * (-(1. / (3. * m)) * prop * p[0])) * (**jit)[3] - (pt * (-(1. / (3. * m)) * prop * p[0])) * (**jit)[4 * SpinorType::R1() + 3] - (pt * (-(1. / (3. * m)) * prop * -M_I * p[0])) * (**jit)[4 * SpinorType::R2() + 3] - (pm * (-(1. / (3. * m)) * prop * -p[0])) * (**jit)[4 * SpinorType::R3() + 3];
        j[4]+= (pm*(-(1./(3. * m)) * prop * p[1])) * (**jit)[0] - ((-ptc) * (-(1. / (3. * m)) * prop * -p[1])) * (**jit)[4 * SpinorType::R1()] - ((-ptc) * (-(1. / (3. * m)) * prop * -M_I * p[1])) * (**jit)[4 * SpinorType::R2()] - (pm * (-(1. / (3. * m)) * prop * -p[1])) * (**jit)[4 * SpinorType::R3()] + ((-ptc) * (-(1. / (3. * m)) * prop * p[1])) * (**jit)[1] - (pm * (-(1. / (3. * m)) * prop * -p[1])) * (**jit)[4 * SpinorType::R1() + 1] - (pm * (-(1. / (3. * m)) * prop * M_I * p[1])) * (**jit)[4 * SpinorType::R2() + 1] - ((-ptc) * (-(1. / (3. * m)) * prop * p[1])) * (**jit)[4 * SpinorType::R3() + 1] + ((r * m) * (-(1. / (3. * m)) * prop * p[1])) * (**jit)[2] - ((r * m) * (-(1. / (3. * m)) * prop * p[1])) * (**jit)[4 * SpinorType::R3() + 2] - ((r * m) * (-(1. / (3. * m)) * prop * p[1])) * (**jit)[4 * SpinorType::R1() + 3] - ((r * m) * (-(1. / (3. * m)) * prop * -M_I * p[1])) * (**jit)[4 * SpinorType::R2() + 3];
        j[5]+= ((-pt)*(-(1./(3. * m)) * prop * p[1])) * (**jit)[0] - (pp * (-(1. / (3. * m)) * prop * -p[1])) * (**jit)[4 * SpinorType::R1()] - (pp * (-(1. / (3. * m)) * prop * -M_I * p[1])) * (**jit)[4 * SpinorType::R2()] - ((-pt) * (-(1. / (3. * m)) * prop * -p[1])) * (**jit)[4 * SpinorType::R3()] + (pp * (-(1. / (3. * m)) * prop * p[1])) * (**jit)[1] - ((-pt) * (-(1. / (3. * m)) * prop * -p[1])) * (**jit)[4 * SpinorType::R1() + 1] - ((-pt) * (-(1. / (3. * m)) * prop * M_I * p[1])) * (**jit)[4 * SpinorType::R2() + 1] - (pp * (-(1. / (3. * m)) * prop * p[1])) * (**jit)[4 * SpinorType::R3() + 1] - ((r * m) * (-(1. / (3. * m)) * prop * p[1])) * (**jit)[4 * SpinorType::R1() + 2] - ((r * m) * (-(1. / (3. * m)) * prop * M_I * p[1])) * (**jit)[4 * SpinorType::R2() + 2] + ((r * m) * (-(1. / (3. * m)) * prop * p[1])) * (**jit)[3] - ((r * m) * (-(1. / (3. * m)) * prop * -p[1])) * (**jit)[4 * SpinorType::R3() + 3];
        j[6]+= ((r * m) * (-(1. / (3. * m)) * prop * p[1])) * (**jit)[0] - ((r * m) * (-(1. / (3. * m)) * prop * -p[1])) * (**jit)[4 * SpinorType::R3()] - ((r * m) * (-(1. / (3. * m)) * prop * -p[1])) * (**jit)[4 * SpinorType::R1() + 1] - ((r * m) * (-(1. / (3. * m)) * prop * M_I * p[1])) * (**jit)[4 * SpinorType::R2() + 1] + (pp * (-(1. / (3. * m)) * prop * p[1])) * (**jit)[2] - (ptc * (-(1. / (3. * m)) * prop * p[1])) * (**jit)[4 * SpinorType::R1() + 2] - (ptc * (-(1. / (3. * m)) * prop * M_I * p[1])) * (**jit)[4 * SpinorType::R2() + 2] - (pp * (-(1. / (3. * m)) * prop * p[1])) * (**jit)[4 * SpinorType::R3() + 2] + (ptc * (-(1. / (3. * m)) * prop * p[1])) * (**jit)[3] - (pp * (-(1. / (3. * m)) * prop * p[1])) * (**jit)[4 * SpinorType::R1() + 3] - (pp * (-(1. / (3. * m)) * prop * -M_I * p[1])) * (**jit)[4 * SpinorType::R2() + 3] - (ptc * (-(1. / (3. * m)) * prop * -p[1])) * (**jit)[4 * SpinorType::R3() + 3];
        j[7]+= - ((r * m) * (-(1. / (3. * m)) * prop * -p[1])) * (**jit)[4 * SpinorType::R1()] - ((r * m) * (-(1. / (3. * m)) * prop * -M_I * p[1])) * (**jit)[4 * SpinorType::R2()] + ((r * m) * (-(1. / (3. * m)) * prop * p[1])) * (**jit)[1] - ((r * m) * (-(1. / (3. * m)) * prop * p[1])) * (**jit)[4 * SpinorType::R3() + 1] + (pt * (-(1. / (3. * m)) * prop * p[1])) * (**jit)[2] - (pm * (-(1. / (3. * m)) * prop * p[1])) * (**jit)[4 * SpinorType::R1() + 2] - (pm * (-(1. / (3. * m)) * prop * M_I * p[1])) * (**jit)[4 * SpinorType::R2() + 2] - (pt * (-(1. / (3. * m)) * prop * p[1])) * (**jit)[4 * SpinorType::R3() + 2] + (pm * (-(1. / (3. * m)) * prop * p[1])) * (**jit)[3] - (pt * (-(1. / (3. * m)) * prop * p[1])) * (**jit)[4 * SpinorType::R1() + 3] - (pt * (-(1. / (3. * m)) * prop * -M_I * p[1])) * (**jit)[4 * SpinorType::R2() + 3] - (pm * (-(1. / (3. * m)) * prop * -p[1])) * (**jit)[4 * SpinorType::R3() + 3];
        j[8]+= (pm*(-(1./(3. * m)) * prop * p[2])) * (**jit)[0] - ((-ptc) * (-(1. / (3. * m)) * prop * -p[2])) * (**jit)[4 * SpinorType::R1()] - ((-ptc) * (-(1. / (3. * m)) * prop * -M_I * p[2])) * (**jit)[4 * SpinorType::R2()] - (pm * (-(1. / (3. * m)) * prop * -p[2])) * (**jit)[4 * SpinorType::R3()] + ((-ptc) * (-(1. / (3. * m)) * prop * p[2])) * (**jit)[1] - (pm * (-(1. / (3. * m)) * prop * -p[2])) * (**jit)[4 * SpinorType::R1() + 1] - (pm * (-(1. / (3. * m)) * prop * M_I * p[2])) * (**jit)[4 * SpinorType::R2() + 1] - ((-ptc) * (-(1. / (3. * m)) * prop * p[2])) * (**jit)[4 * SpinorType::R3() + 1] + ((r * m) * (-(1. / (3. * m)) * prop * p[2])) * (**jit)[2] - ((r * m) * (-(1. / (3. * m)) * prop * p[2])) * (**jit)[4 * SpinorType::R3() + 2] - ((r * m) * (-(1. / (3. * m)) * prop * p[2])) * (**jit)[4 * SpinorType::R1() + 3] - ((r * m) * (-(1. / (3. * m)) * prop * -M_I * p[2])) * (**jit)[4 * SpinorType::R2() + 3];
        j[9]+= ((-pt)*(-(1./(3. * m)) * prop * p[2])) * (**jit)[0] - (pp * (-(1. / (3. * m)) * prop * -p[2])) * (**jit)[4 * SpinorType::R1()] - (pp * (-(1. / (3. * m)) * prop * -M_I * p[2])) * (**jit)[4 * SpinorType::R2()] - ((-pt) * (-(1. / (3. * m)) * prop * -p[2])) * (**jit)[4 * SpinorType::R3()] + (pp * (-(1. / (3. * m)) * prop * p[2])) * (**jit)[1] - ((-pt) * (-(1. / (3. * m)) * prop * -p[2])) * (**jit)[4 * SpinorType::R1() + 1] - ((-pt) * (-(1. / (3. * m)) * prop * M_I * p[2])) * (**jit)[4 * SpinorType::R2() + 1] - (pp * (-(1. / (3. * m)) * prop * p[2])) * (**jit)[4 * SpinorType::R3() + 1] - ((r * m) * (-(1. / (3. * m)) * prop * p[2])) * (**jit)[4 * SpinorType::R1() + 2] - ((r * m) * (-(1. / (3. * m)) * prop * M_I * p[2])) * (**jit)[4 * SpinorType::R2() + 2] + ((r * m) * (-(1. / (3. * m)) * prop * p[2])) * (**jit)[3] - ((r * m) * (-(1. / (3. * m)) * prop * -p[2])) * (**jit)[4 * SpinorType::R3() + 3];
        j[10]+= ((r * m) * (-(1. / (3. * m)) * prop * p[2])) * (**jit)[0] - ((r * m) * (-(1. / (3. * m)) * prop * -p[2])) * (**jit)[4 * SpinorType::R3()] - ((r * m) * (-(1. / (3. * m)) * prop * -p[2])) * (**jit)[4 * SpinorType::R1() + 1] - ((r * m) * (-(1. / (3. * m)) * prop * M_I * p[2])) * (**jit)[4 * SpinorType::R2() + 1] + (pp * (-(1. / (3. * m)) * prop * p[2])) * (**jit)[2] - (ptc * (-(1. / (3. * m)) * prop * p[2])) * (**jit)[4 * SpinorType::R1() + 2] - (ptc * (-(1. / (3. * m)) * prop * M_I * p[2])) * (**jit)[4 * SpinorType::R2() + 2] - (pp * (-(1. / (3. * m)) * prop * p[2])) * (**jit)[4 * SpinorType::R3() + 2] + (ptc * (-(1. / (3. * m)) * prop * p[2])) * (**jit)[3] - (pp * (-(1. / (3. * m)) * prop * p[2])) * (**jit)[4 * SpinorType::R1() + 3] - (pp * (-(1. / (3. * m)) * prop * -M_I * p[2])) * (**jit)[4 * SpinorType::R2() + 3] - (ptc * (-(1. / (3. * m)) * prop * -p[2])) * (**jit)[4 * SpinorType::R3() + 3];
        j[11]+= - ((r * m) * (-(1. / (3. * m)) * prop * -p[2])) * (**jit)[4 * SpinorType::R1()] - ((r * m) * (-(1. / (3. * m)) * prop * -M_I * p[2])) * (**jit)[4 * SpinorType::R2()] + ((r * m) * (-(1. / (3. * m)) * prop * p[2])) * (**jit)[1] - ((r * m) * (-(1. / (3. * m)) * prop * p[2])) * (**jit)[4 * SpinorType::R3() + 1] + (pt * (-(1. / (3. * m)) * prop * p[2])) * (**jit)[2] - (pm * (-(1. / (3. * m)) * prop * p[2])) * (**jit)[4 * SpinorType::R1() + 2] - (pm * (-(1. / (3. * m)) * prop * M_I * p[2])) * (**jit)[4 * SpinorType::R2() + 2] - (pt * (-(1. / (3. * m)) * prop * p[2])) * (**jit)[4 * SpinorType::R3() + 2] + (pm * (-(1. / (3. * m)) * prop * p[2])) * (**jit)[3] - (pt * (-(1. / (3. * m)) * prop * p[2])) * (**jit)[4 * SpinorType::R1() + 3] - (pt * (-(1. / (3. * m)) * prop * -M_I * p[2])) * (**jit)[4 * SpinorType::R2() + 3] - (pm * (-(1. / (3. * m)) * prop * -p[2])) * (**jit)[4 * SpinorType::R3() + 3];
        j[12]+= (pm*(-(1./(3. * m)) * prop * p[3])) * (**jit)[0] - ((-ptc) * (-(1. / (3. * m)) * prop * -p[3])) * (**jit)[4 * SpinorType::R1()] - ((-ptc) * (-(1. / (3. * m)) * prop * -M_I * p[3])) * (**jit)[4 * SpinorType::R2()] - (pm * (-(1. / (3. * m)) * prop * -p[3])) * (**jit)[4 * SpinorType::R3()] + ((-ptc) * (-(1. / (3. * m)) * prop * p[3])) * (**jit)[1] - (pm * (-(1. / (3. * m)) * prop * -p[3])) * (**jit)[4 * SpinorType::R1() + 1] - (pm * (-(1. / (3. * m)) * prop * M_I * p[3])) * (**jit)[4 * SpinorType::R2() + 1] - ((-ptc) * (-(1. / (3. * m)) * prop * p[3])) * (**jit)[4 * SpinorType::R3() + 1] + ((r * m) * (-(1. / (3. * m)) * prop * p[3])) * (**jit)[2] - ((r * m) * (-(1. / (3. * m)) * prop * p[3])) * (**jit)[4 * SpinorType::R3() + 2] - ((r * m) * (-(1. / (3. * m)) * prop * p[3])) * (**jit)[4 * SpinorType::R1() + 3] - ((r * m) * (-(1. / (3. * m)) * prop * -M_I * p[3])) * (**jit)[4 * SpinorType::R2() + 3];
        j[13]+= ((-pt)*(-(1./(3. * m)) * prop * p[3])) * (**jit)[0] - (pp * (-(1. / (3. * m)) * prop * -p[3])) * (**jit)[4 * SpinorType::R1()] - (pp * (-(1. / (3. * m)) * prop * -M_I * p[3])) * (**jit)[4 * SpinorType::R2()] - ((-pt) * (-(1. / (3. * m)) * prop * -p[3])) * (**jit)[4 * SpinorType::R3()] + (pp * (-(1. / (3. * m)) * prop * p[3])) * (**jit)[1] - ((-pt) * (-(1. / (3. * m)) * prop * -p[3])) * (**jit)[4 * SpinorType::R1() + 1] - ((-pt) * (-(1. / (3. * m)) * prop * M_I * p[3])) * (**jit)[4 * SpinorType::R2() + 1] - (pp * (-(1. / (3. * m)) * prop * p[3])) * (**jit)[4 * SpinorType::R3() + 1] - ((r * m) * (-(1. / (3. * m)) * prop * p[3])) * (**jit)[4 * SpinorType::R1() + 2] - ((r * m) * (-(1. / (3. * m)) * prop * M_I * p[3])) * (**jit)[4 * SpinorType::R2() + 2] + ((r * m) * (-(1. / (3. * m)) * prop * p[3])) * (**jit)[3] - ((r * m) * (-(1. / (3. * m)) * prop * -p[3])) * (**jit)[4 * SpinorType::R3() + 3];
        j[14]+= ((r * m) * (-(1. / (3. * m)) * prop * p[3])) * (**jit)[0] - ((r * m) * (-(1. / (3. * m)) * prop * -p[3])) * (**jit)[4 * SpinorType::R3()] - ((r * m) * (-(1. / (3. * m)) * prop * -p[3])) * (**jit)[4 * SpinorType::R1() + 1] - ((r * m) * (-(1. / (3. * m)) * prop * M_I * p[3])) * (**jit)[4 * SpinorType::R2() + 1] + (pp * (-(1. / (3. * m)) * prop * p[3])) * (**jit)[2] - (ptc * (-(1. / (3. * m)) * prop * p[3])) * (**jit)[4 * SpinorType::R1() + 2] - (ptc * (-(1. / (3. * m)) * prop * M_I * p[3])) * (**jit)[4 * SpinorType::R2() + 2] - (pp * (-(1. / (3. * m)) * prop * p[3])) * (**jit)[4 * SpinorType::R3() + 2] + (ptc * (-(1. / (3. * m)) * prop * p[3])) * (**jit)[3] - (pp * (-(1. / (3. * m)) * prop * p[3])) * (**jit)[4 * SpinorType::R1() + 3] - (pp * (-(1. / (3. * m)) * prop * -M_I * p[3])) * (**jit)[4 * SpinorType::R2() + 3] - (ptc * (-(1. / (3. * m)) * prop * -p[3])) * (**jit)[4 * SpinorType::R3() + 3];
        j[15]+= - ((r * m) * (-(1. / (3. * m)) * prop * -p[3])) * (**jit)[4 * SpinorType::R1()] - ((r * m) * (-(1. / (3. * m)) * prop * -M_I * p[3])) * (**jit)[4 * SpinorType::R2()] + ((r * m) * (-(1. / (3. * m)) * prop * p[3])) * (**jit)[1] - ((r * m) * (-(1. / (3. * m)) * prop * p[3])) * (**jit)[4 * SpinorType::R3() + 1] + (pt * (-(1. / (3. * m)) * prop * p[3])) * (**jit)[2] - (pm * (-(1. / (3. * m)) * prop * p[3])) * (**jit)[4 * SpinorType::R1() + 2] - (pm * (-(1. / (3. * m)) * prop * M_I * p[3])) * (**jit)[4 * SpinorType::R2() + 2] - (pt * (-(1. / (3. * m)) * prop * p[3])) * (**jit)[4 * SpinorType::R3() + 2] + (pm * (-(1. / (3. * m)) * prop * p[3])) * (**jit)[3] - (pt * (-(1. / (3. * m)) * prop * p[3])) * (**jit)[4 * SpinorType::R1() + 3] - (pt * (-(1. / (3. * m)) * prop * -M_I * p[3])) * (**jit)[4 * SpinorType::R2() + 3] - (pm * (-(1. / (3. * m)) * prop * -p[3])) * (**jit)[4 * SpinorType::R3() + 3];

        // pnu gammamu
        j[0]+= (pm*((1./(3. * m)) * prop * p[0])) * (**jit)[0] - (pm * ((1. / (3. * m)) * prop * p[1])) * (**jit)[4] - (pm * ((1. / (3. * m)) * prop * p[2])) * (**jit)[8] - (pm * ((1. / (3. * m)) * prop * p[3])) * (**jit)[12] + ((-ptc) * ((1. / (3. * m)) * prop * p[0])) * (**jit)[1] - ((-ptc) * ((1. / (3. * m)) * prop * p[1])) * (**jit)[5] - ((-ptc) * ((1. / (3. * m)) * prop * p[2])) * (**jit)[9] - ((-ptc) * ((1. / (3. * m)) * prop * p[3])) * (**jit)[13] + ((r * m) * ((1. / (3. * m)) * prop * p[0])) * (**jit)[2] - ((r * m) * ((1. / (3. * m)) * prop * p[1])) * (**jit)[6] - ((r * m) * ((1. / (3. * m)) * prop * p[2])) * (**jit)[10] - ((r * m) * ((1. / (3. * m)) * prop * p[3])) * (**jit)[14];
        j[1]+= ((-pt)*((1./(3. * m)) * prop * p[0])) * (**jit)[0] - ((-pt) * ((1. / (3. * m)) * prop * p[1])) * (**jit)[4] - ((-pt) * ((1. / (3. * m)) * prop * p[2])) * (**jit)[8] - ((-pt) * ((1. / (3. * m)) * prop * p[3])) * (**jit)[12] + (pp * ((1. / (3. * m)) * prop * p[0])) * (**jit)[1] - (pp * ((1. / (3. * m)) * prop * p[1])) * (**jit)[5] - (pp * ((1. / (3. * m)) * prop * p[2])) * (**jit)[9] - (pp * ((1. / (3. * m)) * prop * p[3])) * (**jit)[13] + ((r * m) * ((1. / (3. * m)) * prop * p[0])) * (**jit)[3] - ((r * m) * ((1. / (3. * m)) * prop * p[1])) * (**jit)[7] - ((r * m) * ((1. / (3. * m)) * prop * p[2])) * (**jit)[11] - ((r * m) * ((1. / (3. * m)) * prop * p[3])) * (**jit)[15];
        j[2]+= ((r * m) * ((1. / (3. * m)) * prop * p[0])) * (**jit)[0] - ((r * m) * ((1. / (3. * m)) * prop * p[1])) * (**jit)[4] - ((r * m) * ((1. / (3. * m)) * prop * p[2])) * (**jit)[8] - ((r * m) * ((1. / (3. * m)) * prop * p[3])) * (**jit)[12] + (pp * ((1. / (3. * m)) * prop * p[0])) * (**jit)[2] - (pp * ((1. / (3. * m)) * prop * p[1])) * (**jit)[6] - (pp * ((1. / (3. * m)) * prop * p[2])) * (**jit)[10] - (pp * ((1. / (3. * m)) * prop * p[3])) * (**jit)[14] + (ptc * ((1. / (3. * m)) * prop * p[0])) * (**jit)[3] - (ptc * ((1. / (3. * m)) * prop * p[1])) * (**jit)[7] - (ptc * ((1. / (3. * m)) * prop * p[2])) * (**jit)[11] - (ptc * ((1. / (3. * m)) * prop * p[3])) * (**jit)[15];
        j[3]+= ((r * m) * ((1. / (3. * m)) * prop * p[0])) * (**jit)[1] - ((r * m) * ((1. / (3. * m)) * prop * p[1])) * (**jit)[5] - ((r * m) * ((1. / (3. * m)) * prop * p[2])) * (**jit)[9] - ((r * m) * ((1. / (3. * m)) * prop * p[3])) * (**jit)[13] + (pt * ((1. / (3. * m)) * prop * p[0])) * (**jit)[2] - (pt * ((1. / (3. * m)) * prop * p[1])) * (**jit)[6] - (pt * ((1. / (3. * m)) * prop * p[2])) * (**jit)[10] - (pt * ((1. / (3. * m)) * prop * p[3])) * (**jit)[14] + (pm * ((1. / (3. * m)) * prop * p[0])) * (**jit)[3] - (pm * ((1. / (3. * m)) * prop * p[1])) * (**jit)[7] - (pm * ((1. / (3. * m)) * prop * p[2])) * (**jit)[11] - (pm * ((1. / (3. * m)) * prop * p[3])) * (**jit)[15];
        j[4*SpinorType::R1()]+= ((-ptc)*((1./(3. * m)) * prop * (-p[0]))) * (**jit)[0] - ((-ptc) * ((1. / (3. * m)) * prop * (-p[1]))) * (**jit)[4] - ((-ptc) * ((1. / (3. * m)) * prop * (-p[2]))) * (**jit)[8] - ((-ptc) * ((1. / (3. * m)) * prop * (-p[3]))) * (**jit)[12] + (pm * ((1. / (3. * m)) * prop * (-p[0]))) * (**jit)[1] - (pm * ((1. / (3. * m)) * prop * (-p[1]))) * (**jit)[5] - (pm * ((1. / (3. * m)) * prop * (-p[2]))) * (**jit)[9] - (pm * ((1. / (3. * m)) * prop * (-p[3]))) * (**jit)[13] + ((r * m) * ((1. / (3. * m)) * prop * p[0])) * (**jit)[3] - ((r * m) * ((1. / (3. * m)) * prop * p[1])) * (**jit)[7] - ((r * m) * ((1. / (3. * m)) * prop * p[2])) * (**jit)[11] - ((r * m) * ((1. / (3. * m)) * prop * p[3])) * (**jit)[15];
        j[4*SpinorType::R1()+1]+= (pp*((1./(3. * m)) * prop * (-p[0]))) * (**jit)[0] - (pp * ((1. / (3. * m)) * prop * (-p[1]))) * (**jit)[4] - (pp * ((1. / (3. * m)) * prop * (-p[2]))) * (**jit)[8] - (pp * ((1. / (3. * m)) * prop * (-p[3]))) * (**jit)[12] + ((-pt) * ((1. / (3. * m)) * prop * (-p[0]))) * (**jit)[1] - ((-pt) * ((1. / (3. * m)) * prop * (-p[1]))) * (**jit)[5] - ((-pt) * ((1. / (3. * m)) * prop * (-p[2]))) * (**jit)[9] - ((-pt) * ((1. / (3. * m)) * prop * (-p[3]))) * (**jit)[13] + ((r * m) * ((1. / (3. * m)) * prop * p[0])) * (**jit)[2] - ((r * m) * ((1. / (3. * m)) * prop * p[1])) * (**jit)[6] - ((r * m) * ((1. / (3. * m)) * prop * p[2])) * (**jit)[10] - ((r * m) * ((1. / (3. * m)) * prop * p[3])) * (**jit)[14];
        j[4*SpinorType::R1()+2]+= ((r * m) * ((1. / (3. * m)) * prop * (-p[0]))) * (**jit)[1] - ((r * m) * ((1. / (3. * m)) * prop * (-p[1]))) * (**jit)[5] - ((r * m) * ((1. / (3. * m)) * prop * (-p[2]))) * (**jit)[9] - ((r * m) * ((1. / (3. * m)) * prop * (-p[3]))) * (**jit)[13] + (ptc * ((1. / (3. * m)) * prop * p[0])) * (**jit)[2] - (ptc * ((1. / (3. * m)) * prop * p[1])) * (**jit)[6] - (ptc * ((1. / (3. * m)) * prop * p[2])) * (**jit)[10] - (ptc * ((1. / (3. * m)) * prop * p[3])) * (**jit)[14] + (pp * ((1. / (3. * m)) * prop * p[0])) * (**jit)[3] - (pp * ((1. / (3. * m)) * prop * p[1])) * (**jit)[7] - (pp * ((1. / (3. * m)) * prop * p[2])) * (**jit)[11] - (pp * ((1. / (3. * m)) * prop * p[3])) * (**jit)[15];
        j[4*SpinorType::R1()+3]+= ((r * m) * ((1. / (3. * m)) * prop * (-p[0]))) * (**jit)[0] - ((r * m) * ((1. / (3. * m)) * prop * (-p[1]))) * (**jit)[4] - ((r * m) * ((1. / (3. * m)) * prop * (-p[2]))) * (**jit)[8] - ((r * m) * ((1. / (3. * m)) * prop * (-p[3]))) * (**jit)[12] + (pm * ((1. / (3. * m)) * prop * p[0])) * (**jit)[2] - (pm * ((1. / (3. * m)) * prop * p[1])) * (**jit)[6] - (pm * ((1. / (3. * m)) * prop * p[2])) * (**jit)[10] - (pm * ((1. / (3. * m)) * prop * p[3])) * (**jit)[14] + (pt * ((1. / (3. * m)) * prop * p[0])) * (**jit)[3] - (pt * ((1. / (3. * m)) * prop * p[1])) * (**jit)[7] - (pt * ((1. / (3. * m)) * prop * p[2])) * (**jit)[11] - (pt * ((1. / (3. * m)) * prop * p[3])) * (**jit)[15];
        j[4*SpinorType::R2()]+= ((-ptc)*((1./(3. * m)) * prop * (-M_I * p[0]))) * (**jit)[0] - ((-ptc) * ((1. / (3. * m)) * prop * (-M_I * p[1]))) * (**jit)[4] - ((-ptc) * ((1. / (3. * m)) * prop * (-M_I * p[2]))) * (**jit)[8] - ((-ptc) * ((1. / (3. * m)) * prop * (-M_I * p[3]))) * (**jit)[12] + (pm * ((1. / (3. * m)) * prop * M_I * p[0])) * (**jit)[1] - (pm * ((1. / (3. * m)) * prop * M_I * p[1])) * (**jit)[5] - (pm * ((1. / (3. * m)) * prop * M_I * p[2])) * (**jit)[9] - (pm * ((1. / (3. * m)) * prop * M_I * p[3])) * (**jit)[13] + ((r * m) * ((1. / (3. * m)) * prop * (-M_I * p[0]))) * (**jit)[3] - ((r * m) * ((1. / (3. * m)) * prop * (-M_I * p[1]))) * (**jit)[7] - ((r * m) * ((1. / (3. * m)) * prop * (-M_I * p[2]))) * (**jit)[11] - ((r * m) * ((1. / (3. * m)) * prop * (-M_I * p[3]))) * (**jit)[15];
        j[4*SpinorType::R2()+1]+= (pp*((1./(3. * m)) * prop * (-M_I * p[0]))) * (**jit)[0] - (pp * ((1. / (3. * m)) * prop * (-M_I * p[1]))) * (**jit)[4] - (pp * ((1. / (3. * m)) * prop * (-M_I * p[2]))) * (**jit)[8] - (pp * ((1. / (3. * m)) * prop * (-M_I * p[3]))) * (**jit)[12] + ((-pt) * ((1. / (3. * m)) * prop * M_I * p[0])) * (**jit)[1] - ((-pt) * ((1. / (3. * m)) * prop * M_I * p[1])) * (**jit)[5] - ((-pt) * ((1. / (3. * m)) * prop * M_I * p[2])) * (**jit)[9] - ((-pt) * ((1. / (3. * m)) * prop * M_I * p[3])) * (**jit)[13] + ((r * m) * ((1. / (3. * m)) * prop * M_I * p[0])) * (**jit)[2] - ((r * m) * ((1. / (3. * m)) * prop * M_I * p[1])) * (**jit)[6] - ((r * m) * ((1. / (3. * m)) * prop * M_I * p[2])) * (**jit)[10] - ((r * m) * ((1. / (3. * m)) * prop * M_I * p[3])) * (**jit)[14];
        j[4*SpinorType::R2()+2]+= ((r * m) * ((1. / (3. * m)) * prop * M_I * p[0])) * (**jit)[1] - ((r * m) * ((1. / (3. * m)) * prop * M_I * p[1])) * (**jit)[5] - ((r * m) * ((1. / (3. * m)) * prop * M_I * p[2])) * (**jit)[9] - ((r * m) * ((1. / (3. * m)) * prop * M_I * p[3])) * (**jit)[13] + (ptc * ((1. / (3. * m)) * prop * M_I * p[0])) * (**jit)[2] - (ptc * ((1. / (3. * m)) * prop * M_I * p[1])) * (**jit)[6] - (ptc * ((1. / (3. * m)) * prop * M_I * p[2])) * (**jit)[10] - (ptc * ((1. / (3. * m)) * prop * M_I * p[3])) * (**jit)[14] + (pp * ((1. / (3. * m)) * prop * (-M_I * p[0]))) * (**jit)[3] - (pp * ((1. / (3. * m)) * prop * (-M_I * p[1]))) * (**jit)[7] - (pp * ((1. / (3. * m)) * prop * (-M_I * p[2]))) * (**jit)[11] - (pp * ((1. / (3. * m)) * prop * (-M_I * p[3]))) * (**jit)[15];
        j[4*SpinorType::R2()+3]+= ((r * m) * ((1. / (3. * m)) * prop * (-M_I * p[0]))) * (**jit)[0] - ((r * m) * ((1. / (3. * m)) * prop * (-M_I * p[1]))) * (**jit)[4] - ((r * m) * ((1. / (3. * m)) * prop * (-M_I * p[2]))) * (**jit)[8] - ((r * m) * ((1. / (3. * m)) * prop * (-M_I * p[3]))) * (**jit)[12] + (pm * ((1. / (3. * m)) * prop * M_I * p[0])) * (**jit)[2] - (pm * ((1. / (3. * m)) * prop * M_I * p[1])) * (**jit)[6] - (pm * ((1. / (3. * m)) * prop * M_I * p[2])) * (**jit)[10] - (pm * ((1. / (3. * m)) * prop * M_I * p[3])) * (**jit)[14] + (pt * ((1. / (3. * m)) * prop * (-M_I * p[0]))) * (**jit)[3] - (pt * ((1. / (3. * m)) * prop * (-M_I * p[1]))) * (**jit)[7] - (pt * ((1. / (3. * m)) * prop * (-M_I * p[2]))) * (**jit)[11] - (pt * ((1. / (3. * m)) * prop * (-M_I * p[3]))) * (**jit)[15];
        j[4*SpinorType::R3()]+= (pm*((1./(3. * m)) * prop * (-p[0]))) * (**jit)[0] - (pm * ((1. / (3. * m)) * prop * (-p[1]))) * (**jit)[4] - (pm * ((1. / (3. * m)) * prop * (-p[2]))) * (**jit)[8] - (pm * ((1. / (3. * m)) * prop * (-p[3]))) * (**jit)[12] + ((-ptc) * ((1. / (3. * m)) * prop * p[0])) * (**jit)[1] - ((-ptc) * ((1. / (3. * m)) * prop * p[1])) * (**jit)[5] - ((-ptc) * ((1. / (3. * m)) * prop * p[2])) * (**jit)[9] - ((-ptc) * ((1. / (3. * m)) * prop * p[3])) * (**jit)[13] + ((r * m) * ((1. / (3. * m)) * prop * p[0])) * (**jit)[2] - ((r * m) * ((1. / (3. * m)) * prop * p[1])) * (**jit)[6] - ((r * m) * ((1. / (3. * m)) * prop * p[2])) * (**jit)[10] - ((r * m) * ((1. / (3. * m)) * prop * p[3])) * (**jit)[14];
        j[4*SpinorType::R3()+1]+= ((-pt)*((1./(3. * m)) * prop * (-p[0]))) * (**jit)[0] - ((-pt) * ((1. / (3. * m)) * prop * (-p[1]))) * (**jit)[4] - ((-pt) * ((1. / (3. * m)) * prop * (-p[2]))) * (**jit)[8] - ((-pt) * ((1. / (3. * m)) * prop * (-p[3]))) * (**jit)[12] + (pp * ((1. / (3. * m)) * prop * p[0])) * (**jit)[1] - (pp * ((1. / (3. * m)) * prop * p[1])) * (**jit)[5] - (pp * ((1. / (3. * m)) * prop * p[2])) * (**jit)[9] - (pp * ((1. / (3. * m)) * prop * p[3])) * (**jit)[13] + ((r * m) * ((1. / (3. * m)) * prop * (-p[0]))) * (**jit)[3] - ((r * m) * ((1. / (3. * m)) * prop * (-p[1]))) * (**jit)[7] - ((r * m) * ((1. / (3. * m)) * prop * (-p[2]))) * (**jit)[11] - ((r * m) * ((1. / (3. * m)) * prop * (-p[3]))) * (**jit)[15];
        j[4*SpinorType::R3()+2]+= ((r * m) * ((1. / (3. * m)) * prop * (-p[0]))) * (**jit)[0] - ((r * m) * ((1. / (3. * m)) * prop * (-p[1]))) * (**jit)[4] - ((r * m) * ((1. / (3. * m)) * prop * (-p[2]))) * (**jit)[8] - ((r * m) * ((1. / (3. * m)) * prop * (-p[3]))) * (**jit)[12] + (pp * ((1. / (3. * m)) * prop * p[0])) * (**jit)[2] - (pp * ((1. / (3. * m)) * prop * p[1])) * (**jit)[6] - (pp * ((1. / (3. * m)) * prop * p[2])) * (**jit)[10] - (pp * ((1. / (3. * m)) * prop * p[3])) * (**jit)[14] + (ptc * ((1. / (3. * m)) * prop * (-p[0]))) * (**jit)[3] - (ptc * ((1. / (3. * m)) * prop * (-p[1]))) * (**jit)[7] - (ptc * ((1. / (3. * m)) * prop * (-p[2]))) * (**jit)[11] - (ptc * ((1. / (3. * m)) * prop * (-p[3]))) * (**jit)[15];
        j[4*SpinorType::R3()+3]+= ((r * m) * ((1. / (3. * m)) * prop * p[0])) * (**jit)[1] - ((r * m) * ((1. / (3. * m)) * prop * p[1])) * (**jit)[5] - ((r * m) * ((1. / (3. * m)) * prop * p[2])) * (**jit)[9] - ((r * m) * ((1. / (3. * m)) * prop * p[3])) * (**jit)[13] + (pt * ((1. / (3. * m)) * prop * p[0])) * (**jit)[2] - (pt * ((1. / (3. * m)) * prop * p[1])) * (**jit)[6] - (pt * ((1. / (3. * m)) * prop * p[2])) * (**jit)[10] - (pt * ((1. / (3. * m)) * prop * p[3])) * (**jit)[14] + (pm * ((1. / (3. * m)) * prop * (-p[0]))) * (**jit)[3] - (pm * ((1. / (3. * m)) * prop * (-p[1]))) * (**jit)[7] - (pm * ((1. / (3. * m)) * prop * (-p[2]))) * (**jit)[11] - (pm * ((1. / (3. * m)) * prop * (-p[3]))) * (**jit)[15];
      }
      else {// S(p)
        // vector propagator
        j[0]=(**jit)[0]*((r * m) * (-prop + (2. / (3. * m2)) * prop * p[0] * p[0])) - (**jit)[4] * ((r * m) * ((2. / (3. * m2)) * prop * p[1] * p[0])) - (**jit)[8] * ((r * m) * ((2. / (3. * m2)) * prop * p[2] * p[0])) - (**jit)[12] * ((r * m) * ((2. / (3. * m2)) * prop * p[3] * p[0])) + (**jit)[2] * (pp * (-prop + (2. / (3. * m2)) * prop * p[0] * p[0])) - (**jit)[6] * (pp * ((2. / (3. * m2)) * prop * p[1] * p[0])) - (**jit)[10] * (pp * ((2. / (3. * m2)) * prop * p[2] * p[0])) - (**jit)[14] * (pp * ((2. / (3. * m2)) * prop * p[3] * p[0])) + (**jit)[3] * (pt * (-prop + (2. / (3. * m2)) * prop * p[0] * p[0])) - (**jit)[7] * (pt * ((2. / (3. * m2)) * prop * p[1] * p[0])) - (**jit)[11] * (pt * ((2. / (3. * m2)) * prop * p[2] * p[0])) - (**jit)[15] * (pt * ((2. / (3. * m2)) * prop * p[3] * p[0]));
        j[1]=(**jit)[1]*((r * m) * (-prop + (2. / (3. * m2)) * prop * p[0] * p[0])) - (**jit)[5] * ((r * m) * ((2. / (3. * m2)) * prop * p[1] * p[0])) - (**jit)[9] * ((r * m) * ((2. / (3. * m2)) * prop * p[2] * p[0])) - (**jit)[13] * ((r * m) * ((2. / (3. * m2)) * prop * p[3] * p[0])) + (**jit)[2] * (ptc * (-prop + (2. / (3. * m2)) * prop * p[0] * p[0])) - (**jit)[6] * (ptc * ((2. / (3. * m2)) * prop * p[1] * p[0])) - (**jit)[10] * (ptc * ((2. / (3. * m2)) * prop * p[2] * p[0])) - (**jit)[14] * (ptc * ((2. / (3. * m2)) * prop * p[3] * p[0])) + (**jit)[3] * (pm * (-prop + (2. / (3. * m2)) * prop * p[0] * p[0])) - (**jit)[7] * (pm * ((2. / (3. * m2)) * prop * p[1] * p[0])) - (**jit)[11] * (pm * ((2. / (3. * m2)) * prop * p[2] * p[0])) - (**jit)[15] * (pm * ((2. / (3. * m2)) * prop * p[3] * p[0]));
        j[2]=(**jit)[0]*(pm*(-prop+ (2./(3. * m2)) * prop * p[0] * p[0])) - (**jit)[4] * (pm * ((2. / (3. * m2)) * prop * p[1] * p[0])) - (**jit)[8] * (pm * ((2. / (3. * m2)) * prop * p[2] * p[0])) - (**jit)[12] * (pm * ((2. / (3. * m2)) * prop * p[3] * p[0])) + (**jit)[1] * ((-pt) * (-prop + (2. / (3. * m2)) * prop * p[0] * p[0])) - (**jit)[5] * ((-pt) * ((2. / (3. * m2)) * prop * p[1] * p[0])) - (**jit)[9] * ((-pt) * ((2. / (3. * m2)) * prop * p[2] * p[0])) - (**jit)[13] * ((-pt) * ((2. / (3. * m2)) * prop * p[3] * p[0])) + (**jit)[2] * ((r * m) * (-prop + (2. / (3. * m2)) * prop * p[0] * p[0])) - (**jit)[6] * ((r * m) * ((2. / (3. * m2)) * prop * p[1] * p[0])) - (**jit)[10] * ((r * m) * ((2. / (3. * m2)) * prop * p[2] * p[0])) - (**jit)[14] * ((r * m) * ((2. / (3. * m2)) * prop * p[3] * p[0]));
        j[3]=(**jit)[0]*((-ptc)*(-prop+ (2./(3. * m2)) * prop * p[0] * p[0])) - (**jit)[4] * ((-ptc) * ((2. / (3. * m2)) * prop * p[1] * p[0])) - (**jit)[8] * ((-ptc) * ((2. / (3. * m2)) * prop * p[2] * p[0])) - (**jit)[12] * ((-ptc) * ((2. / (3. * m2)) * prop * p[3] * p[0])) + (**jit)[1] * (pp * (-prop + (2. / (3. * m2)) * prop * p[0] * p[0])) - (**jit)[5] * (pp * ((2. / (3. * m2)) * prop * p[1] * p[0])) - (**jit)[9] * (pp * ((2. / (3. * m2)) * prop * p[2] * p[0])) - (**jit)[13] * (pp * ((2. / (3. * m2)) * prop * p[3] * p[0])) + (**jit)[3] * ((r * m) * (-prop + (2. / (3. * m2)) * prop * p[0] * p[0])) - (**jit)[7] * ((r * m) * ((2. / (3. * m2)) * prop * p[1] * p[0])) - (**jit)[11] * ((r * m) * ((2. / (3. * m2)) * prop * p[2] * p[0])) - (**jit)[15] * ((r * m) * ((2. / (3. * m2)) * prop * p[3] * p[0]));
        j[4]=(**jit)[0]*((r * m) * ((2. / (3. * m2)) * prop * p[0] * p[1])) - (**jit)[4] * ((r * m) * (prop + (2. / (3. * m2)) * prop * p[1] * p[1])) - (**jit)[8] * ((r * m) * ((2. / (3. * m2)) * prop * p[2] * p[1])) - (**jit)[12] * ((r * m) * ((2. / (3. * m2)) * prop * p[3] * p[1])) + (**jit)[2] * (pp * ((2. / (3. * m2)) * prop * p[0] * p[1])) - (**jit)[6] * (pp * (prop + (2. / (3. * m2)) * prop * p[1] * p[1])) - (**jit)[10] * (pp * ((2. / (3. * m2)) * prop * p[2] * p[1])) - (**jit)[14] * (pp * ((2. / (3. * m2)) * prop * p[3] * p[1])) + (**jit)[3] * (pt * ((2. / (3. * m2)) * prop * p[0] * p[1])) - (**jit)[7] * (pt * (prop + (2. / (3. * m2)) * prop * p[1] * p[1])) - (**jit)[11] * (pt * ((2. / (3. * m2)) * prop * p[2] * p[1])) - (**jit)[15] * (pt * ((2. / (3. * m2)) * prop * p[3] * p[1]));
        j[5]=(**jit)[1]*((r * m) * ((2. / (3. * m2)) * prop * p[0] * p[1])) - (**jit)[5] * ((r * m) * (prop + (2. / (3. * m2)) * prop * p[1] * p[1])) - (**jit)[9] * ((r * m) * ((2. / (3. * m2)) * prop * p[2] * p[1])) - (**jit)[13] * ((r * m) * ((2. / (3. * m2)) * prop * p[3] * p[1])) + (**jit)[2] * (ptc * ((2. / (3. * m2)) * prop * p[0] * p[1])) - (**jit)[6] * (ptc * (prop + (2. / (3. * m2)) * prop * p[1] * p[1])) - (**jit)[10] * (ptc * ((2. / (3. * m2)) * prop * p[2] * p[1])) - (**jit)[14] * (ptc * ((2. / (3. * m2)) * prop * p[3] * p[1])) + (**jit)[3] * (pm * ((2. / (3. * m2)) * prop * p[0] * p[1])) - (**jit)[7] * (pm * (prop + (2. / (3. * m2)) * prop * p[1] * p[1])) - (**jit)[11] * (pm * ((2. / (3. * m2)) * prop * p[2] * p[1])) - (**jit)[15] * (pm * ((2. / (3. * m2)) * prop * p[3] * p[1]));
        j[6]=(**jit)[0]*(pm*((2./(3. * m2)) * prop * p[0] * p[1])) - (**jit)[4] * (pm * (prop + (2. / (3. * m2)) * prop * p[1] * p[1])) - (**jit)[8] * (pm * ((2. / (3. * m2)) * prop * p[2] * p[1])) - (**jit)[12] * (pm * ((2. / (3. * m2)) * prop * p[3] * p[1])) + (**jit)[1] * ((-pt) * ((2. / (3. * m2)) * prop * p[0] * p[1])) - (**jit)[5] * ((-pt) * (prop + (2. / (3. * m2)) * prop * p[1] * p[1])) - (**jit)[9] * ((-pt) * ((2. / (3. * m2)) * prop * p[2] * p[1])) - (**jit)[13] * ((-pt) * ((2. / (3. * m2)) * prop * p[3] * p[1])) + (**jit)[2] * ((r * m) * ((2. / (3. * m2)) * prop * p[0] * p[1])) - (**jit)[6] * ((r * m) * (prop + (2. / (3. * m2)) * prop * p[1] * p[1])) - (**jit)[10] * ((r * m) * ((2. / (3. * m2)) * prop * p[2] * p[1])) - (**jit)[14] * ((r * m) * ((2. / (3. * m2)) * prop * p[3] * p[1]));
        j[7]=(**jit)[0]*((-ptc)*((2./(3. * m2)) * prop * p[0] * p[1])) - (**jit)[4] * ((-ptc) * (prop + (2. / (3. * m2)) * prop * p[1] * p[1])) - (**jit)[8] * ((-ptc) * ((2. / (3. * m2)) * prop * p[2] * p[1])) - (**jit)[12] * ((-ptc) * ((2. / (3. * m2)) * prop * p[3] * p[1])) + (**jit)[1] * (pp * ((2. / (3. * m2)) * prop * p[0] * p[1])) - (**jit)[5] * (pp * (prop + (2. / (3. * m2)) * prop * p[1] * p[1])) - (**jit)[9] * (pp * ((2. / (3. * m2)) * prop * p[2] * p[1])) - (**jit)[13] * (pp * ((2. / (3. * m2)) * prop * p[3] * p[1])) + (**jit)[3] * ((r * m) * ((2. / (3. * m2)) * prop * p[0] * p[1])) - (**jit)[7] * ((r * m) * (prop + (2. / (3. * m2)) * prop * p[1] * p[1])) - (**jit)[11] * ((r * m) * ((2. / (3. * m2)) * prop * p[2] * p[1])) - (**jit)[15] * ((r * m) * ((2. / (3. * m2)) * prop * p[3] * p[1]));
        j[8]=(**jit)[0]*((r * m) * ((2. / (3. * m2)) * prop * p[0] * p[2])) - (**jit)[4] * ((r * m) * ((2. / (3. * m2)) * prop * p[1] * p[2])) - (**jit)[8] * ((r * m) * (prop + (2. / (3. * m2)) * prop * p[2] * p[2])) - (**jit)[12] * ((r * m) * ((2. / (3. * m2)) * prop * p[3] * p[2])) + (**jit)[2] * (pp * ((2. / (3. * m2)) * prop * p[0] * p[2])) - (**jit)[6] * (pp * ((2. / (3. * m2)) * prop * p[1] * p[2])) - (**jit)[10] * (pp * (prop + (2. / (3. * m2)) * prop * p[2] * p[2])) - (**jit)[14] * (pp * ((2. / (3. * m2)) * prop * p[3] * p[2])) + (**jit)[3] * (pt * ((2. / (3. * m2)) * prop * p[0] * p[2])) - (**jit)[7] * (pt * ((2. / (3. * m2)) * prop * p[1] * p[2])) - (**jit)[11] * (pt * (prop + (2. / (3. * m2)) * prop * p[2] * p[2])) - (**jit)[15] * (pt * ((2. / (3. * m2)) * prop * p[3] * p[2]));
        j[9]=(**jit)[1]*((r * m) * ((2. / (3. * m2)) * prop * p[0] * p[2])) - (**jit)[5] * ((r * m) * ((2. / (3. * m2)) * prop * p[1] * p[2])) - (**jit)[9] * ((r * m) * (prop + (2. / (3. * m2)) * prop * p[2] * p[2])) - (**jit)[13] * ((r * m) * ((2. / (3. * m2)) * prop * p[3] * p[2])) + (**jit)[2] * (ptc * ((2. / (3. * m2)) * prop * p[0] * p[2])) - (**jit)[6] * (ptc * ((2. / (3. * m2)) * prop * p[1] * p[2])) - (**jit)[10] * (ptc * (prop + (2. / (3. * m2)) * prop * p[2] * p[2])) - (**jit)[14] * (ptc * ((2. / (3. * m2)) * prop * p[3] * p[2])) + (**jit)[3] * (pm * ((2. / (3. * m2)) * prop * p[0] * p[2])) - (**jit)[7] * (pm * ((2. / (3. * m2)) * prop * p[1] * p[2])) - (**jit)[11] * (pm * (prop + (2. / (3. * m2)) * prop * p[2] * p[2])) - (**jit)[15] * (pm * ((2. / (3. * m2)) * prop * p[3] * p[2]));
        j[10]=(**jit)[0]*(pm*((2./(3. * m2)) * prop * p[0] * p[2])) - (**jit)[4] * (pm * ((2. / (3. * m2)) * prop * p[1] * p[2])) - (**jit)[8] * (pm * (prop + (2. / (3. * m2)) * prop * p[2] * p[2])) - (**jit)[12] * (pm * ((2. / (3. * m2)) * prop * p[3] * p[2])) + (**jit)[1] * ((-pt) * ((2. / (3. * m2)) * prop * p[0] * p[2])) - (**jit)[5] * ((-pt) * ((2. / (3. * m2)) * prop * p[1] * p[2])) - (**jit)[9] * ((-pt) * (prop + (2. / (3. * m2)) * prop * p[2] * p[2])) - (**jit)[13] * ((-pt) * ((2. / (3. * m2)) * prop * p[3] * p[2])) + (**jit)[2] * ((r * m) * ((2. / (3. * m2)) * prop * p[0] * p[2])) - (**jit)[6] * ((r * m) * ((2. / (3. * m2)) * prop * p[1] * p[2])) - (**jit)[10] * ((r * m) * (prop + (2. / (3. * m2)) * prop * p[2] * p[2])) - (**jit)[14] * ((r * m) * ((2. / (3. * m2)) * prop * p[3] * p[2]));
        j[11]=(**jit)[0]*((-ptc)*((2./(3. * m2)) * prop * p[0] * p[2])) - (**jit)[4] * ((-ptc) * ((2. / (3. * m2)) * prop * p[1] * p[2])) - (**jit)[8] * ((-ptc) * (prop + (2. / (3. * m2)) * prop * p[2] * p[2])) - (**jit)[12] * ((-ptc) * ((2. / (3. * m2)) * prop * p[3] * p[2])) + (**jit)[1] * (pp * ((2. / (3. * m2)) * prop * p[0] * p[2])) - (**jit)[5] * (pp * ((2. / (3. * m2)) * prop * p[1] * p[2])) - (**jit)[9] * (pp * (prop + (2. / (3. * m2)) * prop * p[2] * p[2])) - (**jit)[13] * (pp * ((2. / (3. * m2)) * prop * p[3] * p[2])) + (**jit)[3] * ((r * m) * ((2. / (3. * m2)) * prop * p[0] * p[2])) - (**jit)[7] * ((r * m) * ((2. / (3. * m2)) * prop * p[1] * p[2])) - (**jit)[11] * ((r * m) * (prop + (2. / (3. * m2)) * prop * p[2] * p[2])) - (**jit)[15] * ((r * m) * ((2. / (3. * m2)) * prop * p[3] * p[2]));
        j[12]=(**jit)[0]*((r * m) * ((2. / (3. * m2)) * prop * p[0] * p[3])) - (**jit)[4] * ((r * m) * ((2. / (3. * m2)) * prop * p[1] * p[3])) - (**jit)[8] * ((r * m) * ((2. / (3. * m2)) * prop * p[2] * p[3])) - (**jit)[12] * ((r * m) * (prop + (2. / (3. * m2)) * prop * p[3] * p[3])) + (**jit)[2] * (pp * ((2. / (3. * m2)) * prop * p[0] * p[3])) - (**jit)[6] * (pp * ((2. / (3. * m2)) * prop * p[1] * p[3])) - (**jit)[10] * (pp * ((2. / (3. * m2)) * prop * p[2] * p[3])) - (**jit)[14] * (pp * (prop + (2. / (3. * m2)) * prop * p[3] * p[3])) + (**jit)[3] * (pt * ((2. / (3. * m2)) * prop * p[0] * p[3])) - (**jit)[7] * (pt * ((2. / (3. * m2)) * prop * p[1] * p[3])) - (**jit)[11] * (pt * ((2. / (3. * m2)) * prop * p[2] * p[3])) - (**jit)[15] * (pt * (prop + (2. / (3. * m2)) * prop * p[3] * p[3]));
        j[13]=(**jit)[1]*((r * m) * ((2. / (3. * m2)) * prop * p[0] * p[3])) - (**jit)[5] * ((r * m) * ((2. / (3. * m2)) * prop * p[1] * p[3])) - (**jit)[9] * ((r * m) * ((2. / (3. * m2)) * prop * p[2] * p[3])) - (**jit)[13] * ((r * m) * (prop + (2. / (3. * m2)) * prop * p[3] * p[3])) + (**jit)[2] * (ptc * ((2. / (3. * m2)) * prop * p[0] * p[3])) - (**jit)[6] * (ptc * ((2. / (3. * m2)) * prop * p[1] * p[3])) - (**jit)[10] * (ptc * ((2. / (3. * m2)) * prop * p[2] * p[3])) - (**jit)[14] * (ptc * (prop + (2. / (3. * m2)) * prop * p[3] * p[3])) + (**jit)[3] * (pm * ((2. / (3. * m2)) * prop * p[0] * p[3])) - (**jit)[7] * (pm * ((2. / (3. * m2)) * prop * p[1] * p[3])) - (**jit)[11] * (pm * ((2. / (3. * m2)) * prop * p[2] * p[3])) - (**jit)[15] * (pm * (prop + (2. / (3. * m2)) * prop * p[3] * p[3]));
        j[14]=(**jit)[0]*(pm*((2./(3. * m2)) * prop * p[0] * p[3])) - (**jit)[4] * (pm * ((2. / (3. * m2)) * prop * p[1] * p[3])) - (**jit)[8] * (pm * ((2. / (3. * m2)) * prop * p[2] * p[3])) - (**jit)[12] * (pm * (prop + (2. / (3. * m2)) * prop * p[3] * p[3])) + (**jit)[1] * ((-pt) * ((2. / (3. * m2)) * prop * p[0] * p[3])) - (**jit)[5] * ((-pt) * ((2. / (3. * m2)) * prop * p[1] * p[3])) - (**jit)[9] * ((-pt) * ((2. / (3. * m2)) * prop * p[2] * p[3])) - (**jit)[13] * ((-pt) * (prop + (2. / (3. * m2)) * prop * p[3] * p[3])) + (**jit)[2] * ((r * m) * ((2. / (3. * m2)) * prop * p[0] * p[3])) - (**jit)[6] * ((r * m) * ((2. / (3. * m2)) * prop * p[1] * p[3])) - (**jit)[10] * ((r * m) * ((2. / (3. * m2)) * prop * p[2] * p[3])) - (**jit)[14] * ((r * m) * (prop + (2. / (3. * m2)) * prop * p[3] * p[3]));
        j[15]=(**jit)[0]*((-ptc)*((2./(3. * m2)) * prop * p[0] * p[3])) - (**jit)[4] * ((-ptc) * ((2. / (3. * m2)) * prop * p[1] * p[3])) - (**jit)[8] * ((-ptc) * ((2. / (3. * m2)) * prop * p[2] * p[3])) - (**jit)[12] * ((-ptc) * (prop + (2. / (3. * m2)) * prop * p[3] * p[3])) + (**jit)[1] * (pp * ((2. / (3. * m2)) * prop * p[0] * p[3])) - (**jit)[5] * (pp * ((2. / (3. * m2)) * prop * p[1] * p[3])) - (**jit)[9] * (pp * ((2. / (3. * m2)) * prop * p[2] * p[3])) - (**jit)[13] * (pp * (prop + (2. / (3. * m2)) * prop * p[3] * p[3])) + (**jit)[3] * ((r * m) * ((2. / (3. * m2)) * prop * p[0] * p[3])) - (**jit)[7] * ((r * m) * ((2. / (3. * m2)) * prop * p[1] * p[3])) - (**jit)[11] * ((r * m) * ((2. / (3. * m2)) * prop * p[2] * p[3])) - (**jit)[15] * ((r * m) * (prop + (2. / (3. * m2)) * prop * p[3] * p[3]));

        // gammamu gammanu
        j[0]+=(**jit)[0]*((r * m) * ((1. / 3.) * prop)) - (**jit)[4 * SpinorType::R3()] * ((r * m) * ((1. / 3.) * prop)) - (**jit)[4 * SpinorType::R1() + 1] * ((r * m) * ((1. / 3.) * prop)) - (**jit)[4 * SpinorType::R2() + 1] * ((r * m) * ((M_I / 3.) * prop)) + (**jit)[2] * (pp * ((1. / 3.) * prop)) - (**jit)[4 * SpinorType::R1() + 2] * (ptc * ((1. / 3.) * prop)) - (**jit)[4 * SpinorType::R2() + 2] * (ptc * ((M_I / 3.) * prop)) - (**jit)[4 * SpinorType::R3() + 2] * (pp * ((1. / 3.) * prop)) + (**jit)[3] * (pt * ((1. / 3.) * prop)) - (**jit)[4 * SpinorType::R1() + 3] * (pm * ((1. / 3.) * prop)) - (**jit)[4 * SpinorType::R2() + 3] * (pm * ((M_I / 3.) * prop)) - (**jit)[4 * SpinorType::R3() + 3] * (pt * ((1. / 3.) * prop));
        j[1]+= - (**jit)[4*SpinorType::R1()]*((r * m) * ((1. / 3.) * prop)) - (**jit)[4 * SpinorType::R2()] * ((r * m) * (-(M_I / 3.) * prop)) + (**jit)[1] * ((r * m) * ((1. / 3.) * prop)) - (**jit)[4 * SpinorType::R3() + 1] * ((r * m) * (-(1. / 3.) * prop)) + (**jit)[2] * (ptc * ((1. / 3.) * prop)) - (**jit)[4 * SpinorType::R1() + 2] * (pp * ((1. / 3.) * prop)) - (**jit)[4 * SpinorType::R2() + 2] * (pp * (-(M_I / 3.) * prop)) - (**jit)[4 * SpinorType::R3() + 2] * (ptc * (-(1. / 3.) * prop)) + (**jit)[3] * (pm * ((1. / 3.) * prop)) - (**jit)[4 * SpinorType::R1() + 3] * (pt * ((1. / 3.) * prop)) - (**jit)[4 * SpinorType::R2() + 3] * (pt * (-(M_I / 3.) * prop)) - (**jit)[4 * SpinorType::R3() + 3] * (pm * (-(1. / 3.) * prop));
        j[2]+=(**jit)[0]*(pm*((1./3.)*prop)) - (**jit)[4*SpinorType::R1()]*((-ptc)*(-(1./3.)*prop)) - (**jit)[4*SpinorType::R2()]*((-ptc)*(-(M_I/3.)*prop)) - (**jit)[4*SpinorType::R3()]*(pm*(-(1./3.)*prop)) + (**jit)[1]*((-pt)*((1./3.)*prop)) - (**jit)[4*SpinorType::R1()+1]*(pp*(-(1./3.)*prop)) - (**jit)[4*SpinorType::R2()+1]*(pp*(-(M_I/3.)*prop)) - (**jit)[4*SpinorType::R3()+1]*((-pt)*(-(1./3.)*prop)) + (**jit)[2]*((r * m) * ((1. / 3.) * prop)) - (**jit)[4 * SpinorType::R3() + 2] * ((r * m) * (-(1. / 3.) * prop)) - (**jit)[4 * SpinorType::R1() + 3] * ((r * m) * (-(1. / 3.) * prop)) - (**jit)[4 * SpinorType::R2() + 3] * ((r * m) * (-(M_I / 3.) * prop));
        j[3]+=(**jit)[0]*((-ptc)*((1./3.)*prop)) - (**jit)[4*SpinorType::R1()]*(pm*(-(1./3.)*prop)) - (**jit)[4*SpinorType::R2()]*(pm*((M_I/3.)*prop)) - (**jit)[4*SpinorType::R3()]*((-ptc)*((1./3.)*prop)) + (**jit)[1]*(pp*((1./3.)*prop)) - (**jit)[4*SpinorType::R1()+1]*((-pt)*(-(1./3.)*prop)) - (**jit)[4*SpinorType::R2()+1]*((-pt)*((M_I/3.)*prop)) - (**jit)[4*SpinorType::R3()+1]*(pp*((1./3.)*prop)) - (**jit)[4*SpinorType::R1()+2]*((r * m) * (-(1. / 3.) * prop)) - (**jit)[4 * SpinorType::R2() + 2] * ((r * m) * ((M_I / 3.) * prop)) + (**jit)[3] * ((r * m) * ((1. / 3.) * prop)) - (**jit)[4 * SpinorType::R3() + 3] * ((r * m) * ((1. / 3.) * prop));
        j[4*SpinorType::R1()]+= - (**jit)[4*SpinorType::R1()]*((r * m) * (-(1. / 3.) * prop)) - (**jit)[4 * SpinorType::R2()] * ((r * m) * ((M_I / 3.) * prop)) + (**jit)[1] * ((r * m) * (-(1. / 3.) * prop)) - (**jit)[4 * SpinorType::R3() + 1] * ((r * m) * ((1. / 3.) * prop)) + (**jit)[2] * (ptc * (-(1. / 3.) * prop)) - (**jit)[4 * SpinorType::R1() + 2] * (pp * (-(1. / 3.) * prop)) - (**jit)[4 * SpinorType::R2() + 2] * (pp * ((M_I / 3.) * prop)) - (**jit)[4 * SpinorType::R3() + 2] * (ptc * ((1. / 3.) * prop)) + (**jit)[3] * (pm * (-(1. / 3.) * prop)) - (**jit)[4 * SpinorType::R1() + 3] * (pt * (-(1. / 3.) * prop)) - (**jit)[4 * SpinorType::R2() + 3] * (pt * ((M_I / 3.) * prop)) - (**jit)[4 * SpinorType::R3() + 3] * (pm * ((1. / 3.) * prop));
        j[4*SpinorType::R1()+1]+=(**jit)[0]*((r * m) * (-(1. / 3.) * prop)) - (**jit)[4 * SpinorType::R3()] * ((r * m) * (-(1. / 3.) * prop)) - (**jit)[4 * SpinorType::R1() + 1] * ((r * m) * (-(1. / 3.) * prop)) - (**jit)[4 * SpinorType::R2() + 1] * ((r * m) * (-(M_I / 3.) * prop)) + (**jit)[2] * (pp * (-(1. / 3.) * prop)) - (**jit)[4 * SpinorType::R1() + 2] * (ptc * (-(1. / 3.) * prop)) - (**jit)[4 * SpinorType::R2() + 2] * (ptc * (-(M_I / 3.) * prop)) - (**jit)[4 * SpinorType::R3() + 2] * (pp * (-(1. / 3.) * prop)) + (**jit)[3] * (pt * (-(1. / 3.) * prop)) - (**jit)[4 * SpinorType::R1() + 3] * (pm * (-(1. / 3.) * prop)) - (**jit)[4 * SpinorType::R2() + 3] * (pm * (-(M_I / 3.) * prop)) - (**jit)[4 * SpinorType::R3() + 3] * (pt * (-(1. / 3.) * prop));
        j[4*SpinorType::R1()+2]+=(**jit)[0]*((-ptc)*((1./3.)*prop)) - (**jit)[4*SpinorType::R1()]*(pm*(-(1./3.)*prop)) - (**jit)[4*SpinorType::R2()]*(pm*((M_I/3.)*prop)) - (**jit)[4*SpinorType::R3()]*((-ptc)*((1./3.)*prop)) + (**jit)[1]*(pp*((1./3.)*prop)) - (**jit)[4*SpinorType::R1()+1]*((-pt)*(-(1./3.)*prop)) - (**jit)[4*SpinorType::R2()+1]*((-pt)*((M_I/3.)*prop)) - (**jit)[4*SpinorType::R3()+1]*(pp*((1./3.)*prop)) - (**jit)[4*SpinorType::R1()+2]*((r * m) * (-(1. / 3.) * prop)) - (**jit)[4 * SpinorType::R2() + 2] * ((r * m) * ((M_I / 3.) * prop)) + (**jit)[3] * ((r * m) * ((1. / 3.) * prop)) - (**jit)[4 * SpinorType::R3() + 3] * ((r * m) * ((1. / 3.) * prop));
        j[4*SpinorType::R1()+3]+=(**jit)[0]*(pm*((1./3.)*prop)) - (**jit)[4*SpinorType::R1()]*((-ptc)*(-(1./3.)*prop)) - (**jit)[4*SpinorType::R2()]*((-ptc)*(-(M_I/3.)*prop)) - (**jit)[4*SpinorType::R3()]*(pm*(-(1./3.)*prop)) + (**jit)[1]*((-pt)*((1./3.)*prop)) - (**jit)[4*SpinorType::R1()+1]*(pp*(-(1./3.)*prop)) - (**jit)[4*SpinorType::R2()+1]*(pp*(-(M_I/3.)*prop)) - (**jit)[4*SpinorType::R3()+1]*((-pt)*(-(1./3.)*prop)) + (**jit)[2]*((r * m) * ((1. / 3.) * prop)) - (**jit)[4 * SpinorType::R3() + 2] * ((r * m) * (-(1. / 3.) * prop)) - (**jit)[4 * SpinorType::R1() + 3] * ((r * m) * (-(1. / 3.) * prop)) - (**jit)[4 * SpinorType::R2() + 3] * ((r * m) * (-(M_I / 3.) * prop));
        j[4*SpinorType::R2()]+= - (**jit)[4*SpinorType::R1()]*((r * m) * (-(M_I / 3.) * prop)) - (**jit)[4 * SpinorType::R2()] * ((r * m) * (-(1. / 3.) * prop)) + (**jit)[1] * ((r * m) * (-(M_I / 3.) * prop)) - (**jit)[4 * SpinorType::R3() + 1] * ((r * m) * ((M_I / 3.) * prop)) + (**jit)[2] * (ptc * (-(M_I / 3.) * prop)) - (**jit)[4 * SpinorType::R1() + 2] * (pp * (-(M_I / 3.) * prop)) - (**jit)[4 * SpinorType::R2() + 2] * (pp * (-(1. / 3.) * prop)) - (**jit)[4 * SpinorType::R3() + 2] * (ptc * ((M_I / 3.) * prop)) + (**jit)[3] * (pm * (-(M_I / 3.) * prop)) - (**jit)[4 * SpinorType::R1() + 3] * (pt * (-(M_I / 3.) * prop)) - (**jit)[4 * SpinorType::R2() + 3] * (pt * (-(1. / 3.) * prop)) - (**jit)[4 * SpinorType::R3() + 3] * (pm * ((M_I / 3.) * prop));
        j[4*SpinorType::R2()+1]+=(**jit)[0]*((r * m) * ((M_I / 3.) * prop)) - (**jit)[4 * SpinorType::R3()] * ((r * m) * ((M_I / 3.) * prop)) - (**jit)[4 * SpinorType::R1() + 1] * ((r * m) * ((M_I / 3.) * prop)) - (**jit)[4 * SpinorType::R2() + 1] * ((r * m) * (-(1. / 3.) * prop)) + (**jit)[2] * (pp * ((M_I / 3.) * prop)) - (**jit)[4 * SpinorType::R1() + 2] * (ptc * ((M_I / 3.) * prop)) - (**jit)[4 * SpinorType::R2() + 2] * (ptc * (-(1. / 3.) * prop)) - (**jit)[4 * SpinorType::R3() + 2] * (pp * ((M_I / 3.) * prop)) + (**jit)[3] * (pt * ((M_I / 3.) * prop)) - (**jit)[4 * SpinorType::R1() + 3] * (pm * ((M_I / 3.) * prop)) - (**jit)[4 * SpinorType::R2() + 3] * (pm * (-(1. / 3.) * prop)) - (**jit)[4 * SpinorType::R3() + 3] * (pt * ((M_I / 3.) * prop));
        j[4*SpinorType::R2()+2]+=(**jit)[0]*((-ptc)*((M_I/3.)*prop)) - (**jit)[4*SpinorType::R1()]*(pm*(-(M_I/3.)*prop)) - (**jit)[4*SpinorType::R2()]*(pm*(-(1./3.)*prop)) - (**jit)[4*SpinorType::R3()]*((-ptc)*((M_I/3.)*prop)) + (**jit)[1]*(pp*((M_I/3.)*prop)) - (**jit)[4*SpinorType::R1()+1]*((-pt)*(-(M_I/3.)*prop)) - (**jit)[4*SpinorType::R2()+1]*((-pt)*(-(1./3.)*prop)) - (**jit)[4*SpinorType::R3()+1]*(pp*((M_I/3.)*prop)) - (**jit)[4*SpinorType::R1()+2]*((r * m) * (-(M_I / 3.) * prop)) - (**jit)[4 * SpinorType::R2() + 2] * ((r * m) * (-(1. / 3.) * prop)) + (**jit)[3] * ((r * m) * ((M_I / 3.) * prop)) - (**jit)[4 * SpinorType::R3() + 3] * ((r * m) * ((M_I / 3.) * prop));
        j[4*SpinorType::R2()+3]+=(**jit)[0]*(pm*(-(M_I/3.)*prop)) - (**jit)[4*SpinorType::R1()]*((-ptc)*((M_I/3.)*prop)) - (**jit)[4*SpinorType::R2()]*((-ptc)*(-(1./3.)*prop)) - (**jit)[4*SpinorType::R3()]*(pm*((M_I/3.)*prop)) + (**jit)[1]*((-pt)*(-(M_I/3.)*prop)) - (**jit)[4*SpinorType::R1()+1]*(pp*((M_I/3.)*prop)) - (**jit)[4*SpinorType::R2()+1]*(pp*(-(1./3.)*prop)) - (**jit)[4*SpinorType::R3()+1]*((-pt)*((M_I/3.)*prop)) + (**jit)[2]*((r * m) * (-(M_I / 3.) * prop)) - (**jit)[4 * SpinorType::R3() + 2] * ((r * m) * ((M_I / 3.) * prop)) - (**jit)[4 * SpinorType::R1() + 3] * ((r * m) * ((M_I / 3.) * prop)) - (**jit)[4 * SpinorType::R2() + 3] * ((r * m) * (-(1. / 3.) * prop));
        j[4*SpinorType::R3()]+=(**jit)[0]*((r * m) * (-(1. / 3.) * prop)) - (**jit)[4 * SpinorType::R3()] * ((r * m) * (-(1. / 3.) * prop)) - (**jit)[4 * SpinorType::R1() + 1] * ((r * m) * (-(1. / 3.) * prop)) - (**jit)[4 * SpinorType::R2() + 1] * ((r * m) * (-(M_I / 3.) * prop)) + (**jit)[2] * (pp * (-(1. / 3.) * prop)) - (**jit)[4 * SpinorType::R1() + 2] * (ptc * (-(1. / 3.) * prop)) - (**jit)[4 * SpinorType::R2() + 2] * (ptc * (-(M_I / 3.) * prop)) - (**jit)[4 * SpinorType::R3() + 2] * (pp * (-(1. / 3.) * prop)) + (**jit)[3] * (pt * (-(1. / 3.) * prop)) - (**jit)[4 * SpinorType::R1() + 3] * (pm * (-(1. / 3.) * prop)) - (**jit)[4 * SpinorType::R2() + 3] * (pm * (-(M_I / 3.) * prop)) - (**jit)[4 * SpinorType::R3() + 3] * (pt * (-(1. / 3.) * prop));
        j[4*SpinorType::R3()+1]+= - (**jit)[4*SpinorType::R1()]*((r * m) * ((1. / 3.) * prop)) - (**jit)[4 * SpinorType::R2()] * ((r * m) * (-(M_I / 3.) * prop)) + (**jit)[1] * ((r * m) * ((1. / 3.) * prop)) - (**jit)[4 * SpinorType::R3() + 1] * ((r * m) * (-(1. / 3.) * prop)) + (**jit)[2] * (ptc * ((1. / 3.) * prop)) - (**jit)[4 * SpinorType::R1() + 2] * (pp * ((1. / 3.) * prop)) - (**jit)[4 * SpinorType::R2() + 2] * (pp * (-(M_I / 3.) * prop)) - (**jit)[4 * SpinorType::R3() + 2] * (ptc * (-(1. / 3.) * prop)) + (**jit)[3] * (pm * ((1. / 3.) * prop)) - (**jit)[4 * SpinorType::R1() + 3] * (pt * ((1. / 3.) * prop)) - (**jit)[4 * SpinorType::R2() + 3] * (pt * (-(M_I / 3.) * prop)) - (**jit)[4 * SpinorType::R3() + 3] * (pm * (-(1. / 3.) * prop));
        j[4*SpinorType::R3()+2]+=(**jit)[0]*(pm*((1./3.)*prop)) - (**jit)[4*SpinorType::R1()]*((-ptc)*(-(1./3.)*prop)) - (**jit)[4*SpinorType::R2()]*((-ptc)*(-(M_I/3.)*prop)) - (**jit)[4*SpinorType::R3()]*(pm*(-(1./3.)*prop)) + (**jit)[1]*((-pt)*((1./3.)*prop)) - (**jit)[4*SpinorType::R1()+1]*(pp*(-(1./3.)*prop)) - (**jit)[4*SpinorType::R2()+1]*(pp*(-(M_I/3.)*prop)) - (**jit)[4*SpinorType::R3()+1]*((-pt)*(-(1./3.)*prop)) + (**jit)[2]*((r * m) * ((1. / 3.) * prop)) - (**jit)[4 * SpinorType::R3() + 2] * ((r * m) * (-(1. / 3.) * prop)) - (**jit)[4 * SpinorType::R1() + 3] * ((r * m) * (-(1. / 3.) * prop)) - (**jit)[4 * SpinorType::R2() + 3] * ((r * m) * (-(M_I / 3.) * prop));
        j[4*SpinorType::R3()+3]+=(**jit)[0]*((-ptc)*(-(1./3.)*prop)) - (**jit)[4*SpinorType::R1()]*(pm*((1./3.)*prop)) - (**jit)[4*SpinorType::R2()]*(pm*(-(M_I/3.)*prop)) - (**jit)[4*SpinorType::R3()]*((-ptc)*(-(1./3.)*prop)) + (**jit)[1]*(pp*(-(1./3.)*prop)) - (**jit)[4*SpinorType::R1()+1]*((-pt)*((1./3.)*prop)) - (**jit)[4*SpinorType::R2()+1]*((-pt)*(-(M_I/3.)*prop)) - (**jit)[4*SpinorType::R3()+1]*(pp*(-(1./3.)*prop)) - (**jit)[4*SpinorType::R1()+2]*((r * m) * ((1. / 3.) * prop)) - (**jit)[4 * SpinorType::R2() + 2] * ((r * m) * (-(M_I / 3.) * prop)) + (**jit)[3] * ((r * m) * (-(1. / 3.) * prop)) - (**jit)[4 * SpinorType::R3() + 3] * ((r * m) * (-(1. / 3.) * prop));

        // pmu gammanu
        j[0]+=(**jit)[0]*(pm*(-(1./(3. * m)) * prop * p[0])) - (**jit)[4] * (pm * (-(1. / (3. * m)) * prop * p[1])) - (**jit)[8] * (pm * (-(1. / (3. * m)) * prop * p[2])) - (**jit)[12] * (pm * (-(1. / (3. * m)) * prop * p[3])) + (**jit)[1] * ((-pt) * (-(1. / (3. * m)) * prop * p[0])) - (**jit)[5] * ((-pt) * (-(1. / (3. * m)) * prop * p[1])) - (**jit)[9] * ((-pt) * (-(1. / (3. * m)) * prop * p[2])) - (**jit)[13] * ((-pt) * (-(1. / (3. * m)) * prop * p[3])) + (**jit)[2] * ((r * m) * (-(1. / (3. * m)) * prop * p[0])) - (**jit)[6] * ((r * m) * (-(1. / (3. * m)) * prop * p[1])) - (**jit)[10] * ((r * m) * (-(1. / (3. * m)) * prop * p[2])) - (**jit)[14] * ((r * m) * (-(1. / (3. * m)) * prop * p[3]));
        j[1]+=(**jit)[0]*((-ptc)*(-(1./(3. * m)) * prop * p[0])) - (**jit)[4] * ((-ptc) * (-(1. / (3. * m)) * prop * p[1])) - (**jit)[8] * ((-ptc) * (-(1. / (3. * m)) * prop * p[2])) - (**jit)[12] * ((-ptc) * (-(1. / (3. * m)) * prop * p[3])) + (**jit)[1] * (pp * (-(1. / (3. * m)) * prop * p[0])) - (**jit)[5] * (pp * (-(1. / (3. * m)) * prop * p[1])) - (**jit)[9] * (pp * (-(1. / (3. * m)) * prop * p[2])) - (**jit)[13] * (pp * (-(1. / (3. * m)) * prop * p[3])) + (**jit)[3] * ((r * m) * (-(1. / (3. * m)) * prop * p[0])) - (**jit)[7] * ((r * m) * (-(1. / (3. * m)) * prop * p[1])) - (**jit)[11] * ((r * m) * (-(1. / (3. * m)) * prop * p[2])) - (**jit)[15] * ((r * m) * (-(1. / (3. * m)) * prop * p[3]));
        j[2]+=(**jit)[0]*((r * m) * (-(1. / (3. * m)) * prop * p[0])) - (**jit)[4] * ((r * m) * (-(1. / (3. * m)) * prop * p[1])) - (**jit)[8] * ((r * m) * (-(1. / (3. * m)) * prop * p[2])) - (**jit)[12] * ((r * m) * (-(1. / (3. * m)) * prop * p[3])) + (**jit)[2] * (pp * (-(1. / (3. * m)) * prop * p[0])) - (**jit)[6] * (pp * (-(1. / (3. * m)) * prop * p[1])) - (**jit)[10] * (pp * (-(1. / (3. * m)) * prop * p[2])) - (**jit)[14] * (pp * (-(1. / (3. * m)) * prop * p[3])) + (**jit)[3] * (pt * (-(1. / (3. * m)) * prop * p[0])) - (**jit)[7] * (pt * (-(1. / (3. * m)) * prop * p[1])) - (**jit)[11] * (pt * (-(1. / (3. * m)) * prop * p[2])) - (**jit)[15] * (pt * (-(1. / (3. * m)) * prop * p[3]));
        j[3]+=(**jit)[1]*((r * m) * (-(1. / (3. * m)) * prop * p[0])) - (**jit)[5] * ((r * m) * (-(1. / (3. * m)) * prop * p[1])) - (**jit)[9] * ((r * m) * (-(1. / (3. * m)) * prop * p[2])) - (**jit)[13] * ((r * m) * (-(1. / (3. * m)) * prop * p[3])) + (**jit)[2] * (ptc * (-(1. / (3. * m)) * prop * p[0])) - (**jit)[6] * (ptc * (-(1. / (3. * m)) * prop * p[1])) - (**jit)[10] * (ptc * (-(1. / (3. * m)) * prop * p[2])) - (**jit)[14] * (ptc * (-(1. / (3. * m)) * prop * p[3])) + (**jit)[3] * (pm * (-(1. / (3. * m)) * prop * p[0])) - (**jit)[7] * (pm * (-(1. / (3. * m)) * prop * p[1])) - (**jit)[11] * (pm * (-(1. / (3. * m)) * prop * p[2])) - (**jit)[15] * (pm * (-(1. / (3. * m)) * prop * p[3]));
        j[4*SpinorType::R1()]+=(**jit)[0]*((-ptc)*(-(1./(3. * m)) * prop * -p[0])) - (**jit)[4] * ((-ptc) * (-(1. / (3. * m)) * prop * -p[1])) - (**jit)[8] * ((-ptc) * (-(1. / (3. * m)) * prop * -p[2])) - (**jit)[12] * ((-ptc) * (-(1. / (3. * m)) * prop * -p[3])) + (**jit)[1] * (pp * (-(1. / (3. * m)) * prop * -p[0])) - (**jit)[5] * (pp * (-(1. / (3. * m)) * prop * -p[1])) - (**jit)[9] * (pp * (-(1. / (3. * m)) * prop * -p[2])) - (**jit)[13] * (pp * (-(1. / (3. * m)) * prop * -p[3])) + (**jit)[3] * ((r * m) * (-(1. / (3. * m)) * prop * -p[0])) - (**jit)[7] * ((r * m) * (-(1. / (3. * m)) * prop * -p[1])) - (**jit)[11] * ((r * m) * (-(1. / (3. * m)) * prop * -p[2])) - (**jit)[15] * ((r * m) * (-(1. / (3. * m)) * prop * -p[3]));
        j[4*SpinorType::R1()+1]+=(**jit)[0]*(pm*(-(1./(3. * m)) * prop * -p[0])) - (**jit)[4] * (pm * (-(1. / (3. * m)) * prop * -p[1])) - (**jit)[8] * (pm * (-(1. / (3. * m)) * prop * -p[2])) - (**jit)[12] * (pm * (-(1. / (3. * m)) * prop * -p[3])) + (**jit)[1] * ((-pt) * (-(1. / (3. * m)) * prop * -p[0])) - (**jit)[5] * ((-pt) * (-(1. / (3. * m)) * prop * -p[1])) - (**jit)[9] * ((-pt) * (-(1. / (3. * m)) * prop * -p[2])) - (**jit)[13] * ((-pt) * (-(1. / (3. * m)) * prop * -p[3])) + (**jit)[2] * ((r * m) * (-(1. / (3. * m)) * prop * -p[0])) - (**jit)[6] * ((r * m) * (-(1. / (3. * m)) * prop * -p[1])) - (**jit)[10] * ((r * m) * (-(1. / (3. * m)) * prop * -p[2])) - (**jit)[14] * ((r * m) * (-(1. / (3. * m)) * prop * -p[3]));
        j[4*SpinorType::R1()+2]+=(**jit)[1]*((r * m) * (-(1. / (3. * m)) * prop * p[0])) - (**jit)[5] * ((r * m) * (-(1. / (3. * m)) * prop * p[1])) - (**jit)[9] * ((r * m) * (-(1. / (3. * m)) * prop * p[2])) - (**jit)[13] * ((r * m) * (-(1. / (3. * m)) * prop * p[3])) + (**jit)[2] * (ptc * (-(1. / (3. * m)) * prop * p[0])) - (**jit)[6] * (ptc * (-(1. / (3. * m)) * prop * p[1])) - (**jit)[10] * (ptc * (-(1. / (3. * m)) * prop * p[2])) - (**jit)[14] * (ptc * (-(1. / (3. * m)) * prop * p[3])) + (**jit)[3] * (pm * (-(1. / (3. * m)) * prop * p[0])) - (**jit)[7] * (pm * (-(1. / (3. * m)) * prop * p[1])) - (**jit)[11] * (pm * (-(1. / (3. * m)) * prop * p[2])) - (**jit)[15] * (pm * (-(1. / (3. * m)) * prop * p[3]));
        j[4*SpinorType::R1()+3]+=(**jit)[0]*((r * m) * (-(1. / (3. * m)) * prop * p[0])) - (**jit)[4] * ((r * m) * (-(1. / (3. * m)) * prop * p[1])) - (**jit)[8] * ((r * m) * (-(1. / (3. * m)) * prop * p[2])) - (**jit)[12] * ((r * m) * (-(1. / (3. * m)) * prop * p[3])) + (**jit)[2] * (pp * (-(1. / (3. * m)) * prop * p[0])) - (**jit)[6] * (pp * (-(1. / (3. * m)) * prop * p[1])) - (**jit)[10] * (pp * (-(1. / (3. * m)) * prop * p[2])) - (**jit)[14] * (pp * (-(1. / (3. * m)) * prop * p[3])) + (**jit)[3] * (pt * (-(1. / (3. * m)) * prop * p[0])) - (**jit)[7] * (pt * (-(1. / (3. * m)) * prop * p[1])) - (**jit)[11] * (pt * (-(1. / (3. * m)) * prop * p[2])) - (**jit)[15] * (pt * (-(1. / (3. * m)) * prop * p[3]));
        j[4*SpinorType::R2()]+=(**jit)[0]*((-ptc)*(-(1./(3. * m)) * prop * -M_I * p[0])) - (**jit)[4] * ((-ptc) * (-(1. / (3. * m)) * prop * -M_I * p[1])) - (**jit)[8] * ((-ptc) * (-(1. / (3. * m)) * prop * -M_I * p[2])) - (**jit)[12] * ((-ptc) * (-(1. / (3. * m)) * prop * -M_I * p[3])) + (**jit)[1] * (pp * (-(1. / (3. * m)) * prop * -M_I * p[0])) - (**jit)[5] * (pp * (-(1. / (3. * m)) * prop * -M_I * p[1])) - (**jit)[9] * (pp * (-(1. / (3. * m)) * prop * -M_I * p[2])) - (**jit)[13] * (pp * (-(1. / (3. * m)) * prop * -M_I * p[3])) + (**jit)[3] * ((r * m) * (-(1. / (3. * m)) * prop * -M_I * p[0])) - (**jit)[7] * ((r * m) * (-(1. / (3. * m)) * prop * -M_I * p[1])) - (**jit)[11] * ((r * m) * (-(1. / (3. * m)) * prop * -M_I * p[2])) - (**jit)[15] * ((r * m) * (-(1. / (3. * m)) * prop * -M_I * p[3]));
        j[4*SpinorType::R2()+1]+=(**jit)[0]*(pm*(-(1./(3. * m)) * prop * M_I * p[0])) - (**jit)[4] * (pm * (-(1. / (3. * m)) * prop * M_I * p[1])) - (**jit)[8] * (pm * (-(1. / (3. * m)) * prop * M_I * p[2])) - (**jit)[12] * (pm * (-(1. / (3. * m)) * prop * M_I * p[3])) + (**jit)[1] * ((-pt) * (-(1. / (3. * m)) * prop * M_I * p[0])) - (**jit)[5] * ((-pt) * (-(1. / (3. * m)) * prop * M_I * p[1])) - (**jit)[9] * ((-pt) * (-(1. / (3. * m)) * prop * M_I * p[2])) - (**jit)[13] * ((-pt) * (-(1. / (3. * m)) * prop * M_I * p[3])) + (**jit)[2] * ((r * m) * (-(1. / (3. * m)) * prop * M_I * p[0])) - (**jit)[6] * ((r * m) * (-(1. / (3. * m)) * prop * M_I * p[1])) - (**jit)[10] * ((r * m) * (-(1. / (3. * m)) * prop * M_I * p[2])) - (**jit)[14] * ((r * m) * (-(1. / (3. * m)) * prop * M_I * p[3]));
        j[4*SpinorType::R2()+2]+=(**jit)[1]*((r * m) * (-(1. / (3. * m)) * prop * M_I * p[0])) - (**jit)[5] * ((r * m) * (-(1. / (3. * m)) * prop * M_I * p[1])) - (**jit)[9] * ((r * m) * (-(1. / (3. * m)) * prop * M_I * p[2])) - (**jit)[13] * ((r * m) * (-(1. / (3. * m)) * prop * M_I * p[3])) + (**jit)[2] * (ptc * (-(1. / (3. * m)) * prop * M_I * p[0])) - (**jit)[6] * (ptc * (-(1. / (3. * m)) * prop * M_I * p[1])) - (**jit)[10] * (ptc * (-(1. / (3. * m)) * prop * M_I * p[2])) - (**jit)[14] * (ptc * (-(1. / (3. * m)) * prop * M_I * p[3])) + (**jit)[3] * (pm * (-(1. / (3. * m)) * prop * M_I * p[0])) - (**jit)[7] * (pm * (-(1. / (3. * m)) * prop * M_I * p[1])) - (**jit)[11] * (pm * (-(1. / (3. * m)) * prop * M_I * p[2])) - (**jit)[15] * (pm * (-(1. / (3. * m)) * prop * M_I * p[3]));
        j[4*SpinorType::R2()+3]+=(**jit)[0]*((r * m) * (-(1. / (3. * m)) * prop * -M_I * p[0])) - (**jit)[4] * ((r * m) * (-(1. / (3. * m)) * prop * -M_I * p[1])) - (**jit)[8] * ((r * m) * (-(1. / (3. * m)) * prop * -M_I * p[2])) - (**jit)[12] * ((r * m) * (-(1. / (3. * m)) * prop * -M_I * p[3])) + (**jit)[2] * (pp * (-(1. / (3. * m)) * prop * -M_I * p[0])) - (**jit)[6] * (pp * (-(1. / (3. * m)) * prop * -M_I * p[1])) - (**jit)[10] * (pp * (-(1. / (3. * m)) * prop * -M_I * p[2])) - (**jit)[14] * (pp * (-(1. / (3. * m)) * prop * -M_I * p[3])) + (**jit)[3] * (pt * (-(1. / (3. * m)) * prop * -M_I * p[0])) - (**jit)[7] * (pt * (-(1. / (3. * m)) * prop * -M_I * p[1])) - (**jit)[11] * (pt * (-(1. / (3. * m)) * prop * -M_I * p[2])) - (**jit)[15] * (pt * (-(1. / (3. * m)) * prop * -M_I * p[3]));
        j[4*SpinorType::R3()]+=(**jit)[0]*(pm*(-(1./(3. * m)) * prop * -p[0])) - (**jit)[4] * (pm * (-(1. / (3. * m)) * prop * -p[1])) - (**jit)[8] * (pm * (-(1. / (3. * m)) * prop * -p[2])) - (**jit)[12] * (pm * (-(1. / (3. * m)) * prop * -p[3])) + (**jit)[1] * ((-pt) * (-(1. / (3. * m)) * prop * -p[0])) - (**jit)[5] * ((-pt) * (-(1. / (3. * m)) * prop * -p[1])) - (**jit)[9] * ((-pt) * (-(1. / (3. * m)) * prop * -p[2])) - (**jit)[13] * ((-pt) * (-(1. / (3. * m)) * prop * -p[3])) + (**jit)[2] * ((r * m) * (-(1. / (3. * m)) * prop * -p[0])) - (**jit)[6] * ((r * m) * (-(1. / (3. * m)) * prop * -p[1])) - (**jit)[10] * ((r * m) * (-(1. / (3. * m)) * prop * -p[2])) - (**jit)[14] * ((r * m) * (-(1. / (3. * m)) * prop * -p[3]));
        j[4*SpinorType::R3()+1]+=(**jit)[0]*((-ptc)*(-(1./(3. * m)) * prop * p[0])) - (**jit)[4] * ((-ptc) * (-(1. / (3. * m)) * prop * p[1])) - (**jit)[8] * ((-ptc) * (-(1. / (3. * m)) * prop * p[2])) - (**jit)[12] * ((-ptc) * (-(1. / (3. * m)) * prop * p[3])) + (**jit)[1] * (pp * (-(1. / (3. * m)) * prop * p[0])) - (**jit)[5] * (pp * (-(1. / (3. * m)) * prop * p[1])) - (**jit)[9] * (pp * (-(1. / (3. * m)) * prop * p[2])) - (**jit)[13] * (pp * (-(1. / (3. * m)) * prop * p[3])) + (**jit)[3] * ((r * m) * (-(1. / (3. * m)) * prop * p[0])) - (**jit)[7] * ((r * m) * (-(1. / (3. * m)) * prop * p[1])) - (**jit)[11] * ((r * m) * (-(1. / (3. * m)) * prop * p[2])) - (**jit)[15] * ((r * m) * (-(1. / (3. * m)) * prop * p[3]));
        j[4*SpinorType::R3()+2]+=(**jit)[0]*((r * m) * (-(1. / (3. * m)) * prop * p[0])) - (**jit)[4] * ((r * m) * (-(1. / (3. * m)) * prop * p[1])) - (**jit)[8] * ((r * m) * (-(1. / (3. * m)) * prop * p[2])) - (**jit)[12] * ((r * m) * (-(1. / (3. * m)) * prop * p[3])) + (**jit)[2] * (pp * (-(1. / (3. * m)) * prop * p[0])) - (**jit)[6] * (pp * (-(1. / (3. * m)) * prop * p[1])) - (**jit)[10] * (pp * (-(1. / (3. * m)) * prop * p[2])) - (**jit)[14] * (pp * (-(1. / (3. * m)) * prop * p[3])) + (**jit)[3] * (pt * (-(1. / (3. * m)) * prop * p[0])) - (**jit)[7] * (pt * (-(1. / (3. * m)) * prop * p[1])) - (**jit)[11] * (pt * (-(1. / (3. * m)) * prop * p[2])) - (**jit)[15] * (pt * (-(1. / (3. * m)) * prop * p[3]));
        j[4*SpinorType::R3()+3]+=(**jit)[1]*((r * m) * (-(1. / (3. * m)) * prop * -p[0])) - (**jit)[5] * ((r * m) * (-(1. / (3. * m)) * prop * -p[1])) - (**jit)[9] * ((r * m) * (-(1. / (3. * m)) * prop * -p[2])) - (**jit)[13] * ((r * m) * (-(1. / (3. * m)) * prop * -p[3])) + (**jit)[2] * (ptc * (-(1. / (3. * m)) * prop * -p[0])) - (**jit)[6] * (ptc * (-(1. / (3. * m)) * prop * -p[1])) - (**jit)[10] * (ptc * (-(1. / (3. * m)) * prop * -p[2])) - (**jit)[14] * (ptc * (-(1. / (3. * m)) * prop * -p[3])) + (**jit)[3] * (pm * (-(1. / (3. * m)) * prop * -p[0])) - (**jit)[7] * (pm * (-(1. / (3. * m)) * prop * -p[1])) - (**jit)[11] * (pm * (-(1. / (3. * m)) * prop * -p[2])) - (**jit)[15] * (pm * (-(1. / (3. * m)) * prop * -p[3]));

        // pnu gammamu
        j[0]+=(**jit)[0]*(pm*((1./(3. * m)) * prop * p[0])) - (**jit)[4 * SpinorType::R1()] * ((-ptc) * ((1. / (3. * m)) * prop * (-p[0]))) - (**jit)[4 * SpinorType::R2()] * ((-ptc) * ((1. / (3. * m)) * prop * (-M_I * p[0]))) - (**jit)[4 * SpinorType::R3()] * (pm * ((1. / (3. * m)) * prop * (-p[0]))) + (**jit)[1] * ((-pt) * ((1. / (3. * m)) * prop * p[0])) - (**jit)[4 * SpinorType::R1() + 1] * (pp * ((1. / (3. * m)) * prop * (-p[0]))) - (**jit)[4 * SpinorType::R2() + 1] * (pp * ((1. / (3. * m)) * prop * (-M_I * p[0]))) - (**jit)[4 * SpinorType::R3() + 1] * ((-pt) * ((1. / (3. * m)) * prop * (-p[0]))) + (**jit)[2] * ((r * m) * ((1. / (3. * m)) * prop * p[0])) - (**jit)[4 * SpinorType::R3() + 2] * ((r * m) * ((1. / (3. * m)) * prop * (-p[0]))) - (**jit)[4 * SpinorType::R1() + 3] * ((r * m) * ((1. / (3. * m)) * prop * (-p[0]))) - (**jit)[4 * SpinorType::R2() + 3] * ((r * m) * ((1. / (3. * m)) * prop * (-M_I * p[0])));
        j[1]+=(**jit)[0]*((-ptc)*((1./(3. * m)) * prop * p[0])) - (**jit)[4 * SpinorType::R1()] * (pm * ((1. / (3. * m)) * prop * (-p[0]))) - (**jit)[4 * SpinorType::R2()] * (pm * ((1. / (3. * m)) * prop * M_I * p[0])) - (**jit)[4 * SpinorType::R3()] * ((-ptc) * ((1. / (3. * m)) * prop * p[0])) + (**jit)[1] * (pp * ((1. / (3. * m)) * prop * p[0])) - (**jit)[4 * SpinorType::R1() + 1] * ((-pt) * ((1. / (3. * m)) * prop * (-p[0]))) - (**jit)[4 * SpinorType::R2() + 1] * ((-pt) * ((1. / (3. * m)) * prop * M_I * p[0])) - (**jit)[4 * SpinorType::R3() + 1] * (pp * ((1. / (3. * m)) * prop * p[0])) - (**jit)[4 * SpinorType::R1() + 2] * ((r * m) * ((1. / (3. * m)) * prop * (-p[0]))) - (**jit)[4 * SpinorType::R2() + 2] * ((r * m) * ((1. / (3. * m)) * prop * M_I * p[0])) + (**jit)[3] * ((r * m) * ((1. / (3. * m)) * prop * p[0])) - (**jit)[4 * SpinorType::R3() + 3] * ((r * m) * ((1. / (3. * m)) * prop * p[0]));
        j[2]+=(**jit)[0]*((r * m) * ((1. / (3. * m)) * prop * p[0])) - (**jit)[4 * SpinorType::R3()] * ((r * m) * ((1. / (3. * m)) * prop * p[0])) - (**jit)[4 * SpinorType::R1() + 1] * ((r * m) * ((1. / (3. * m)) * prop * p[0])) - (**jit)[4 * SpinorType::R2() + 1] * ((r * m) * ((1. / (3. * m)) * prop * M_I * p[0])) + (**jit)[2] * (pp * ((1. / (3. * m)) * prop * p[0])) - (**jit)[4 * SpinorType::R1() + 2] * (ptc * ((1. / (3. * m)) * prop * p[0])) - (**jit)[4 * SpinorType::R2() + 2] * (ptc * ((1. / (3. * m)) * prop * M_I * p[0])) - (**jit)[4 * SpinorType::R3() + 2] * (pp * ((1. / (3. * m)) * prop * p[0])) + (**jit)[3] * (pt * ((1. / (3. * m)) * prop * p[0])) - (**jit)[4 * SpinorType::R1() + 3] * (pm * ((1. / (3. * m)) * prop * p[0])) - (**jit)[4 * SpinorType::R2() + 3] * (pm * ((1. / (3. * m)) * prop * M_I * p[0])) - (**jit)[4 * SpinorType::R3() + 3] * (pt * ((1. / (3. * m)) * prop * p[0]));
        j[3]+= - (**jit)[4*SpinorType::R1()]*((r * m) * ((1. / (3. * m)) * prop * p[0])) - (**jit)[4 * SpinorType::R2()] * ((r * m) * ((1. / (3. * m)) * prop * (-M_I * p[0]))) + (**jit)[1] * ((r * m) * ((1. / (3. * m)) * prop * p[0])) - (**jit)[4 * SpinorType::R3() + 1] * ((r * m) * ((1. / (3. * m)) * prop * (-p[0]))) + (**jit)[2] * (ptc * ((1. / (3. * m)) * prop * p[0])) - (**jit)[4 * SpinorType::R1() + 2] * (pp * ((1. / (3. * m)) * prop * p[0])) - (**jit)[4 * SpinorType::R2() + 2] * (pp * ((1. / (3. * m)) * prop * (-M_I * p[0]))) - (**jit)[4 * SpinorType::R3() + 2] * (ptc * ((1. / (3. * m)) * prop * (-p[0]))) + (**jit)[3] * (pm * ((1. / (3. * m)) * prop * p[0])) - (**jit)[4 * SpinorType::R1() + 3] * (pt * ((1. / (3. * m)) * prop * p[0])) - (**jit)[4 * SpinorType::R2() + 3] * (pt * ((1. / (3. * m)) * prop * (-M_I * p[0]))) - (**jit)[4 * SpinorType::R3() + 3] * (pm * ((1. / (3. * m)) * prop * (-p[0])));
        j[4]+=(**jit)[0]*(pm*((1./(3. * m)) * prop * p[1])) - (**jit)[4 * SpinorType::R1()] * ((-ptc) * ((1. / (3. * m)) * prop * (-p[1]))) - (**jit)[4 * SpinorType::R2()] * ((-ptc) * ((1. / (3. * m)) * prop * (-M_I * p[1]))) - (**jit)[4 * SpinorType::R3()] * (pm * ((1. / (3. * m)) * prop * (-p[1]))) + (**jit)[1] * ((-pt) * ((1. / (3. * m)) * prop * p[1])) - (**jit)[4 * SpinorType::R1() + 1] * (pp * ((1. / (3. * m)) * prop * (-p[1]))) - (**jit)[4 * SpinorType::R2() + 1] * (pp * ((1. / (3. * m)) * prop * (-M_I * p[1]))) - (**jit)[4 * SpinorType::R3() + 1] * ((-pt) * ((1. / (3. * m)) * prop * (-p[1]))) + (**jit)[2] * ((r * m) * ((1. / (3. * m)) * prop * p[1])) - (**jit)[4 * SpinorType::R3() + 2] * ((r * m) * ((1. / (3. * m)) * prop * (-p[1]))) - (**jit)[4 * SpinorType::R1() + 3] * ((r * m) * ((1. / (3. * m)) * prop * (-p[1]))) - (**jit)[4 * SpinorType::R2() + 3] * ((r * m) * ((1. / (3. * m)) * prop * (-M_I * p[1])));
        j[5]+=(**jit)[0]*((-ptc)*((1./(3. * m)) * prop * p[1])) - (**jit)[4 * SpinorType::R1()] * (pm * ((1. / (3. * m)) * prop * (-p[1]))) - (**jit)[4 * SpinorType::R2()] * (pm * ((1. / (3. * m)) * prop * M_I * p[1])) - (**jit)[4 * SpinorType::R3()] * ((-ptc) * ((1. / (3. * m)) * prop * p[1])) + (**jit)[1] * (pp * ((1. / (3. * m)) * prop * p[1])) - (**jit)[4 * SpinorType::R1() + 1] * ((-pt) * ((1. / (3. * m)) * prop * (-p[1]))) - (**jit)[4 * SpinorType::R2() + 1] * ((-pt) * ((1. / (3. * m)) * prop * M_I * p[1])) - (**jit)[4 * SpinorType::R3() + 1] * (pp * ((1. / (3. * m)) * prop * p[1])) - (**jit)[4 * SpinorType::R1() + 2] * ((r * m) * ((1. / (3. * m)) * prop * (-p[1]))) - (**jit)[4 * SpinorType::R2() + 2] * ((r * m) * ((1. / (3. * m)) * prop * M_I * p[1])) + (**jit)[3] * ((r * m) * ((1. / (3. * m)) * prop * p[1])) - (**jit)[4 * SpinorType::R3() + 3] * ((r * m) * ((1. / (3. * m)) * prop * p[1]));
        j[6]+=(**jit)[0]*((r * m) * ((1. / (3. * m)) * prop * p[1])) - (**jit)[4 * SpinorType::R3()] * ((r * m) * ((1. / (3. * m)) * prop * p[1])) - (**jit)[4 * SpinorType::R1() + 1] * ((r * m) * ((1. / (3. * m)) * prop * p[1])) - (**jit)[4 * SpinorType::R2() + 1] * ((r * m) * ((1. / (3. * m)) * prop * M_I * p[1])) + (**jit)[2] * (pp * ((1. / (3. * m)) * prop * p[1])) - (**jit)[4 * SpinorType::R1() + 2] * (ptc * ((1. / (3. * m)) * prop * p[1])) - (**jit)[4 * SpinorType::R2() + 2] * (ptc * ((1. / (3. * m)) * prop * M_I * p[1])) - (**jit)[4 * SpinorType::R3() + 2] * (pp * ((1. / (3. * m)) * prop * p[1])) + (**jit)[3] * (pt * ((1. / (3. * m)) * prop * p[1])) - (**jit)[4 * SpinorType::R1() + 3] * (pm * ((1. / (3. * m)) * prop * p[1])) - (**jit)[4 * SpinorType::R2() + 3] * (pm * ((1. / (3. * m)) * prop * M_I * p[1])) - (**jit)[4 * SpinorType::R3() + 3] * (pt * ((1. / (3. * m)) * prop * p[1]));
        j[7]+= - (**jit)[4*SpinorType::R1()]*((r * m) * ((1. / (3. * m)) * prop * p[1])) - (**jit)[4 * SpinorType::R2()] * ((r * m) * ((1. / (3. * m)) * prop * (-M_I * p[1]))) + (**jit)[1] * ((r * m) * ((1. / (3. * m)) * prop * p[1])) - (**jit)[4 * SpinorType::R3() + 1] * ((r * m) * ((1. / (3. * m)) * prop * (-p[1]))) + (**jit)[2] * (ptc * ((1. / (3. * m)) * prop * p[1])) - (**jit)[4 * SpinorType::R1() + 2] * (pp * ((1. / (3. * m)) * prop * p[1])) - (**jit)[4 * SpinorType::R2() + 2] * (pp * ((1. / (3. * m)) * prop * (-M_I * p[1]))) - (**jit)[4 * SpinorType::R3() + 2] * (ptc * ((1. / (3. * m)) * prop * (-p[1]))) + (**jit)[3] * (pm * ((1. / (3. * m)) * prop * p[1])) - (**jit)[4 * SpinorType::R1() + 3] * (pt * ((1. / (3. * m)) * prop * p[1])) - (**jit)[4 * SpinorType::R2() + 3] * (pt * ((1. / (3. * m)) * prop * (-M_I * p[1]))) - (**jit)[4 * SpinorType::R3() + 3] * (pm * ((1. / (3. * m)) * prop * (-p[1])));
        j[8]+=(**jit)[0]*(pm*((1./(3. * m)) * prop * p[2])) - (**jit)[4 * SpinorType::R1()] * ((-ptc) * ((1. / (3. * m)) * prop * (-p[2]))) - (**jit)[4 * SpinorType::R2()] * ((-ptc) * ((1. / (3. * m)) * prop * (-M_I * p[2]))) - (**jit)[4 * SpinorType::R3()] * (pm * ((1. / (3. * m)) * prop * (-p[2]))) + (**jit)[1] * ((-pt) * ((1. / (3. * m)) * prop * p[2])) - (**jit)[4 * SpinorType::R1() + 1] * (pp * ((1. / (3. * m)) * prop * (-p[2]))) - (**jit)[4 * SpinorType::R2() + 1] * (pp * ((1. / (3. * m)) * prop * (-M_I * p[2]))) - (**jit)[4 * SpinorType::R3() + 1] * ((-pt) * ((1. / (3. * m)) * prop * (-p[2]))) + (**jit)[2] * ((r * m) * ((1. / (3. * m)) * prop * p[2])) - (**jit)[4 * SpinorType::R3() + 2] * ((r * m) * ((1. / (3. * m)) * prop * (-p[2]))) - (**jit)[4 * SpinorType::R1() + 3] * ((r * m) * ((1. / (3. * m)) * prop * (-p[2]))) - (**jit)[4 * SpinorType::R2() + 3] * ((r * m) * ((1. / (3. * m)) * prop * (-M_I * p[2])));
        j[9]+=(**jit)[0]*((-ptc)*((1./(3. * m)) * prop * p[2])) - (**jit)[4 * SpinorType::R1()] * (pm * ((1. / (3. * m)) * prop * (-p[2]))) - (**jit)[4 * SpinorType::R2()] * (pm * ((1. / (3. * m)) * prop * M_I * p[2])) - (**jit)[4 * SpinorType::R3()] * ((-ptc) * ((1. / (3. * m)) * prop * p[2])) + (**jit)[1] * (pp * ((1. / (3. * m)) * prop * p[2])) - (**jit)[4 * SpinorType::R1() + 1] * ((-pt) * ((1. / (3. * m)) * prop * (-p[2]))) - (**jit)[4 * SpinorType::R2() + 1] * ((-pt) * ((1. / (3. * m)) * prop * M_I * p[2])) - (**jit)[4 * SpinorType::R3() + 1] * (pp * ((1. / (3. * m)) * prop * p[2])) - (**jit)[4 * SpinorType::R1() + 2] * ((r * m) * ((1. / (3. * m)) * prop * (-p[2]))) - (**jit)[4 * SpinorType::R2() + 2] * ((r * m) * ((1. / (3. * m)) * prop * M_I * p[2])) + (**jit)[3] * ((r * m) * ((1. / (3. * m)) * prop * p[2])) - (**jit)[4 * SpinorType::R3() + 3] * ((r * m) * ((1. / (3. * m)) * prop * p[2]));
        j[10]+=(**jit)[0]*((r * m) * ((1. / (3. * m)) * prop * p[2])) - (**jit)[4 * SpinorType::R3()] * ((r * m) * ((1. / (3. * m)) * prop * p[2])) - (**jit)[4 * SpinorType::R1() + 1] * ((r * m) * ((1. / (3. * m)) * prop * p[2])) - (**jit)[4 * SpinorType::R2() + 1] * ((r * m) * ((1. / (3. * m)) * prop * M_I * p[2])) + (**jit)[2] * (pp * ((1. / (3. * m)) * prop * p[2])) - (**jit)[4 * SpinorType::R1() + 2] * (ptc * ((1. / (3. * m)) * prop * p[2])) - (**jit)[4 * SpinorType::R2() + 2] * (ptc * ((1. / (3. * m)) * prop * M_I * p[2])) - (**jit)[4 * SpinorType::R3() + 2] * (pp * ((1. / (3. * m)) * prop * p[2])) + (**jit)[3] * (pt * ((1. / (3. * m)) * prop * p[2])) - (**jit)[4 * SpinorType::R1() + 3] * (pm * ((1. / (3. * m)) * prop * p[2])) - (**jit)[4 * SpinorType::R2() + 3] * (pm * ((1. / (3. * m)) * prop * M_I * p[2])) - (**jit)[4 * SpinorType::R3() + 3] * (pt * ((1. / (3. * m)) * prop * p[2]));
        j[11]+= - (**jit)[4*SpinorType::R1()]*((r * m) * ((1. / (3. * m)) * prop * p[2])) - (**jit)[4 * SpinorType::R2()] * ((r * m) * ((1. / (3. * m)) * prop * (-M_I * p[2]))) + (**jit)[1] * ((r * m) * ((1. / (3. * m)) * prop * p[2])) - (**jit)[4 * SpinorType::R3() + 1] * ((r * m) * ((1. / (3. * m)) * prop * (-p[2]))) + (**jit)[2] * (ptc * ((1. / (3. * m)) * prop * p[2])) - (**jit)[4 * SpinorType::R1() + 2] * (pp * ((1. / (3. * m)) * prop * p[2])) - (**jit)[4 * SpinorType::R2() + 2] * (pp * ((1. / (3. * m)) * prop * (-M_I * p[2]))) - (**jit)[4 * SpinorType::R3() + 2] * (ptc * ((1. / (3. * m)) * prop * (-p[2]))) + (**jit)[3] * (pm * ((1. / (3. * m)) * prop * p[2])) - (**jit)[4 * SpinorType::R1() + 3] * (pt * ((1. / (3. * m)) * prop * p[2])) - (**jit)[4 * SpinorType::R2() + 3] * (pt * ((1. / (3. * m)) * prop * (-M_I * p[2]))) - (**jit)[4 * SpinorType::R3() + 3] * (pm * ((1. / (3. * m)) * prop * (-p[2])));
        j[12]+=(**jit)[0]*(pm*((1./(3. * m)) * prop * p[3])) - (**jit)[4 * SpinorType::R1()] * ((-ptc) * ((1. / (3. * m)) * prop * (-p[3]))) - (**jit)[4 * SpinorType::R2()] * ((-ptc) * ((1. / (3. * m)) * prop * (-M_I * p[3]))) - (**jit)[4 * SpinorType::R3()] * (pm * ((1. / (3. * m)) * prop * (-p[3]))) + (**jit)[1] * ((-pt) * ((1. / (3. * m)) * prop * p[3])) - (**jit)[4 * SpinorType::R1() + 1] * (pp * ((1. / (3. * m)) * prop * (-p[3]))) - (**jit)[4 * SpinorType::R2() + 1] * (pp * ((1. / (3. * m)) * prop * (-M_I * p[3]))) - (**jit)[4 * SpinorType::R3() + 1] * ((-pt) * ((1. / (3. * m)) * prop * (-p[3]))) + (**jit)[2] * ((r * m) * ((1. / (3. * m)) * prop * p[3])) - (**jit)[4 * SpinorType::R3() + 2] * ((r * m) * ((1. / (3. * m)) * prop * (-p[3]))) - (**jit)[4 * SpinorType::R1() + 3] * ((r * m) * ((1. / (3. * m)) * prop * (-p[3]))) - (**jit)[4 * SpinorType::R2() + 3] * ((r * m) * ((1. / (3. * m)) * prop * (-M_I * p[3])));
        j[13]+=(**jit)[0]*((-ptc)*((1./(3. * m)) * prop * p[3])) - (**jit)[4 * SpinorType::R1()] * (pm * ((1. / (3. * m)) * prop * (-p[3]))) - (**jit)[4 * SpinorType::R2()] * (pm * ((1. / (3. * m)) * prop * M_I * p[3])) - (**jit)[4 * SpinorType::R3()] * ((-ptc) * ((1. / (3. * m)) * prop * p[3])) + (**jit)[1] * (pp * ((1. / (3. * m)) * prop * p[3])) - (**jit)[4 * SpinorType::R1() + 1] * ((-pt) * ((1. / (3. * m)) * prop * (-p[3]))) - (**jit)[4 * SpinorType::R2() + 1] * ((-pt) * ((1. / (3. * m)) * prop * M_I * p[3])) - (**jit)[4 * SpinorType::R3() + 1] * (pp * ((1. / (3. * m)) * prop * p[3])) - (**jit)[4 * SpinorType::R1() + 2] * ((r * m) * ((1. / (3. * m)) * prop * (-p[3]))) - (**jit)[4 * SpinorType::R2() + 2] * ((r * m) * ((1. / (3. * m)) * prop * M_I * p[3])) + (**jit)[3] * ((r * m) * ((1. / (3. * m)) * prop * p[3])) - (**jit)[4 * SpinorType::R3() + 3] * ((r * m) * ((1. / (3. * m)) * prop * p[3]));
        j[14]+=(**jit)[0]*((r * m) * ((1. / (3. * m)) * prop * p[3])) - (**jit)[4 * SpinorType::R3()] * ((r * m) * ((1. / (3. * m)) * prop * p[3])) - (**jit)[4 * SpinorType::R1() + 1] * ((r * m) * ((1. / (3. * m)) * prop * p[3])) - (**jit)[4 * SpinorType::R2() + 1] * ((r * m) * ((1. / (3. * m)) * prop * M_I * p[3])) + (**jit)[2] * (pp * ((1. / (3. * m)) * prop * p[3])) - (**jit)[4 * SpinorType::R1() + 2] * (ptc * ((1. / (3. * m)) * prop * p[3])) - (**jit)[4 * SpinorType::R2() + 2] * (ptc * ((1. / (3. * m)) * prop * M_I * p[3])) - (**jit)[4 * SpinorType::R3() + 2] * (pp * ((1. / (3. * m)) * prop * p[3])) + (**jit)[3] * (pt * ((1. / (3. * m)) * prop * p[3])) - (**jit)[4 * SpinorType::R1() + 3] * (pm * ((1. / (3. * m)) * prop * p[3])) - (**jit)[4 * SpinorType::R2() + 3] * (pm * ((1. / (3. * m)) * prop * M_I * p[3])) - (**jit)[4 * SpinorType::R3() + 3] * (pt * ((1. / (3. * m)) * prop * p[3]));
        j[15]+= - (**jit)[4*SpinorType::R1()]*((r * m) * ((1. / (3. * m)) * prop * p[3])) - (**jit)[4 * SpinorType::R2()] * ((r * m) * ((1. / (3. * m)) * prop * (-M_I * p[3]))) + (**jit)[1] * ((r * m) * ((1. / (3. * m)) * prop * p[3])) - (**jit)[4 * SpinorType::R3() + 1] * ((r * m) * ((1. / (3. * m)) * prop * (-p[3]))) + (**jit)[2] * (ptc * ((1. / (3. * m)) * prop * p[3])) - (**jit)[4 * SpinorType::R1() + 2] * (pp * ((1. / (3. * m)) * prop * p[3])) - (**jit)[4 * SpinorType::R2() + 2] * (pp * ((1. / (3. * m)) * prop * (-M_I * p[3]))) - (**jit)[4 * SpinorType::R3() + 2] * (ptc * ((1. / (3. * m)) * prop * (-p[3]))) + (**jit)[3] * (pm * ((1. / (3. * m)) * prop * p[3])) - (**jit)[4 * SpinorType::R1() + 3] * (pt * ((1. / (3. * m)) * prop * p[3])) - (**jit)[4 * SpinorType::R2() + 3] * (pt * ((1. / (3. * m)) * prop * (-M_I * p[3]))) - (**jit)[4 * SpinorType::R3() + 3] * (pm * ((1. / (3. * m)) * prop * (-p[3])));
      }

#ifdef DEBUG__BG
      msg_Debugging() << METHOD << " Testing completeness relation of "
      "propagator: ";
      if (fabs(p.Mass()-m_cmass.real())<m_cmass.real()*1e-10){
        CRaritaSchwinger<SType> rspp = RSPP(this->m_p, 1, 1, 1);
        CRaritaSchwinger<SType> rsp = RSP(this->m_p, 1, 1, 1);
        CRaritaSchwinger<SType> rsm = RSM(this->m_p, 1, 1, 1);
        CRaritaSchwinger<SType> rsmm = RSMM(this->m_p, 1, 1, 1);
        std::vector<Complex> j_test(16);
        for (size_t B(0); B<16; ++B) {
          for (size_t A(0); A < 16; ++A) {
            SComplex fac = prop;
            if (A>3) fac = -prop;
            if ((*jit)->B()>0)
              j_test[B] += fac * (rspp[B] * rspp.Bar()[A] + rsp[B]
                * rsp.Bar()[A] + rsm[B] * rsm.Bar()[A] + rsmm[B]
                * rsmm.Bar()[A]) * (**jit)[A];
            else
              j_test[B] += fac * (**jit)[A] * (rspp[A] * rspp.Bar()[B] + rsp[A]
                * rsp.Bar()[B] + rsm[A] * rsm.Bar()[B] + rsmm[A]
                * rsmm.Bar()[B]);
            }
          if (std::abs((j_test[B] - j[B]).real()) >
          j.Accu()*std::abs(j[B].real()) || std::abs((j_test[B] - j[B]).imag())
          > j.Accu()*std::abs(j[B].imag()))
          msg_Out() << std::setprecision(12) << "Component " << B <<
          " is not the same between completeness relation "
            << j_test[B] << " and propagator " << j[B] << std::endl;
        }
        msg_Debugging() << " ... done: " << std::endl;
      }
      else
        msg_Debugging() << "Test only works for on-shell particles, consider "
                           "using DECAY_OS to set your intermediate RS "
                           "particle on-shell!" << std::endl;
#endif
      **jit=j*prop1;
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

/* #PYTHON SCRIPT TO DETERMINE THE COMPONENTS OF THE PROPAGATOR NUMERATOR
#resulting components were first tested and will
#subsequently be further simplified by Mathematica; the result can be checked against the original expression from the
#python script whether the same values are achieved

# two values can be adjusted in that script: gauge switched between the different COMIX_DEFAULT_GAUGES, b denotes,
# whether the current, the propagator should be contracted with, should be contracted from left (b<0) or right (b>0)

import numpy as np

# gamma matrices in Sherpas Weyl basis
gamma0 = [[0., 0., 1.0, 0.], [0., 0., 0., 1.], [1., 0., 0., 0.], [0., 1., 0., 0.]]
gamma1 = [[0., 0., 0., 1.], [0., 0., 1., 0.], [0., -1., 0., 0.], [-1., 0., 0., 0.]]
gamma2 = [[0., 0., 0., -1j], [0., 0., 1j, 0.], [0., 1j, 0., 0.], [-1j, 0., 0., 0.]]
gamma3 = [[0., 0., 1.0, 0.], [0., 0., 0., -1.], [-1., 0., 0., 0.], [0., 1., 0., 0.]]

# RS wave function to contract with (written down in matrix form psi[A][mu] for a simpler calculation, components
# corresponds to the corresponding components in the vector notation of the wave function in Sherpa
# within the whole script, upper case latin letters denote spinor, greek letters Lorentz indices
psi = [["(**jit)[0]", "(**jit)[4]", "(**jit)[8]", "(**jit)[12]"],
 ["(**jit)[1]", "(**jit)[5]", "(**jit)[9]", "(**jit)[13]"], ["(**jit)[2]", "(**jit)[6]", "(**jit)[10]", "(**jit)[14]"],
 ["(**jit)[3]", "(**jit)[7]", "(**jit)[11]", "(**jit)[15]"]]
# variable name for the imaginary unit
imag_var = "M_I"

# allows for the calculation of the propagator for different COMIX_DEFAULT_GAUGES
gauge = 0 #1 #2
if gauge == 0:
    gamma = np.array([gamma0, gamma1, gamma2, gamma3], dtype=complex)
    R1 = 1
    R2 = 2
    R3 = 3
elif gauge == 1:
    gamma = np.array([gamma0, gamma3, gamma1, gamma2], dtype=complex)
    R1 = 2
    R2 = 3
    R3 = 1
elif gauge == 2:
    gamma = np.array([gamma0, gamma2, gamma3, gamma1], dtype=complex)
    R1 = 3
    R2 = 1
    R3 = 2

# whether current should be contracted from the left (b<0) or right (b>0)
b = 1

# momentum of RS particle
p = ["p[0]", "p[1]", "p[2]", "p[3]"]

# gmunu-Term
g_term = np.zeros(shape=(4, 4, 4, 4))
for i in range(4):
    if i == 0:
        prop = -1
    else:
        prop = 1
    for j in range(4):
        g_term[j][j][i][i] = 1 * prop

# gamma_mu gamma_nu term
gamma_prod = np.zeros(shape=(4, 4, 4, 4), dtype=complex)
for A in range(4):
    for mu in range(4):
        for B in range(4):
            for nu in range(4):
                for C in range(4):
                    gamma_prod[A][B][mu][nu] += gamma[mu][A][C] * gamma[nu][C][B]

# p_mu gamma_nu term
pmu_gammanu = [[[['' for _ in range(4)] for _ in range(4)] for _ in range(4)] for _ in range(4)]
for A in range(4):
    for B in range(4):
        for mu in range(4):
            for nu in range(4):
                if (gamma[nu][A][B] != complex(0., 0.)):
                    if (gamma[nu][A][B] == complex(1., 0.)):
                        pmu_gammanu[A][B][mu][nu] = p[mu]
                    elif (gamma[nu][A][B] == complex(-1., 0.)):
                        pmu_gammanu[A][B][mu][nu] = "-" + p[mu]
                    elif (gamma[nu][A][B] == complex(0., 1.)):
                        pmu_gammanu[A][B][mu][nu] = imag_var + "*" + p[mu]
                    elif (gamma[nu][A][B] == complex(0., -1.)):
                        pmu_gammanu[A][B][mu][nu] = "-" + imag_var + "*" + p[mu]
                    else:
                        print("Error in pmu gammanu term: Gamma has different values than expected!")
                        exit(1)
for i in range(4):
    for j in range(4):
        print("mu nu", i, j, pmu_gammanu[i][j][0][1])
# gamma_mu p_nu term
pnu_gammamu = [[[['' for _ in range(4)] for _ in range(4)] for _ in range(4)] for _ in range(4)]
for A in range(4):
    for B in range(4):
        for mu in range(4):
            for nu in range(4):
                if (gamma[mu][A][B] != complex(0., 0.)):
                    if (gamma[mu][A][B] == complex(1., 0.)):
                        pnu_gammamu[A][B][mu][nu] = p[nu]
                    elif (gamma[mu][A][B] == complex(-1., 0.)):
                        pnu_gammamu[A][B][mu][nu] = "-" + p[nu]
                    elif (gamma[mu][A][B] == complex(0., 1.)):
                        pnu_gammamu[A][B][mu][nu] = imag_var + "*" + p[nu]
                    elif (gamma[mu][A][B] == complex(0., -1.)):
                        pnu_gammamu[A][B][mu][nu] = "-" + imag_var + "*" + p[nu]
                    else:
                        print("Error in pnu gammamu term: Gamma has different values than expected!")
                        exit(1)
for i in range(4):
    for j in range(4):
        print("nu mu", i, j, pnu_gammamu[i][j][0][1])
# pmu pnu term
p_term = [[[['' for _ in range(4)] for _ in range(4)] for _ in range(4)] for _ in range(4)]
for A in range(4):
    for mu in range(4):
        for B in range(4):
            for nu in range(4):
                if A==B:
                    p_term[A][B][mu][nu] = p[mu] + "*" + p[nu]

# pdagger gammamu gammanu
# gamma_rho gamma_mu gamma_nu term
gamma_prod3 = np.zeros(shape=(4, 4, 4, 4, 4), dtype=complex)
for A in range(4):
    for mu in range(4):
        for B in range(4):
            for nu in range(4):
                for rho in range(4):
                    for C in range(4):
                        gamma_prod3[A][B][mu][nu][rho] += gamma[rho][A][C] * gamma_prod[C][B][mu][nu]

pdagger_gammamu_gammanu = [[[['' for _ in range(4)] for _ in range(4)] for _ in range(4)] for _ in range(4)]
for A in range(4):
    for mu in range(4):
        for B in range(4):
            for nu in range(4):
                for rho in range(4):
                    if rho == 0:
                        prop = 1
                    else:
                        prop = -1
                    if gamma_prod3[A][B][mu][nu][rho] != complex(0.0):
                        if gamma_prod3[A][B][mu][nu][rho] * prop == complex(1., 0.):
                            if pdagger_gammamu_gammanu[A][B][mu][nu] != '':
                                pdagger_gammamu_gammanu[A][B][mu][nu] += "+"
                            pdagger_gammamu_gammanu[A][B][mu][nu] += p[rho]
                        elif gamma_prod3[A][B][mu][nu][rho] * prop == complex(-1., 0.):
                            pdagger_gammamu_gammanu[A][B][mu][nu] += "-" + p[rho]
                        elif gamma_prod3[A][B][mu][nu][rho] * prop == complex(0., 1.):
                            if pdagger_gammamu_gammanu[A][B][mu][nu] != '':
                                pdagger_gammamu_gammanu[A][B][mu][nu] += "+"
                            pdagger_gammamu_gammanu[A][B][mu][nu] += imag_var + "*" + p[rho]
                        elif gamma_prod3[A][B][mu][nu][rho] * prop == complex(0., -1.):
                            pdagger_gammamu_gammanu[A][B][mu][nu] += "-" + imag_var + "*" + p[rho]
                        else:
                            if pdagger_gammamu_gammanu[A][B][mu][nu] != '':
                                pdagger_gammamu_gammanu[A][B][mu][nu] += "+"
                            pdagger_gammamu_gammanu[A][B][mu][nu] += str(gamma_prod3[A][B][mu][nu][rho] * prop) + "*" \
                                                                     + p[rho]
# fermion propagator part
# r=1 for particles, r=-1 for anti-particles
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

spropagator = [['' for _ in range(4)] for _ in range(4)]
for A in range(4):
    for B in range(4):
        if gamma_p[A][B] !='':
            spropagator[A][B] += gamma_p[A][B]
        if (A==B):
            if spropagator[A][B] != "":
                spropagator[A][B] += "+"
            spropagator[A][B] += "r*mRS"
print(spropagator)
# sum all terms of the propagator later multiplied by the pure fermionic part
# prop contains the denominator of this propagator part
lorentz_prop = [[[['' for _ in range(4)] for _ in range(4)] for _ in range(4)] for _ in range(4)]
for A in range(4):
    for mu in range(4):
        for B in range(4):
            for nu in range(4):
                if g_term[A][B][mu][nu] != 0.0:
                    if g_term[A][B][mu][nu] == 1.0:
                        lorentz_prop[A][B][mu][nu] += "prop"
                    elif g_term[A][B][mu][nu] == -1.0:
                            lorentz_prop[A][B][mu][nu] += "-prop"
                    else:
                        print("Strange values of metric tensor")
                        exit(1)

                if gamma_prod[A][B][mu][nu] != complex(0.0):
                    if gamma_prod[A][B][mu][nu] == complex(1., 0.):
                        if lorentz_prop[A][B][mu][nu] != "":
                            lorentz_prop[A][B][mu][nu] += "+"
                        lorentz_prop[A][B][mu][nu] += "(1./3.)" + "*" + "prop"
                    elif gamma_prod[A][B][mu][nu] == complex(-1., 0.):
                        lorentz_prop[A][B][mu][nu] += "-(1./3.)" + "*" + "prop"
                    elif gamma_prod[A][B][mu][nu] == complex(0., 1.):
                        if lorentz_prop[A][B][mu][nu] != "":
                            lorentz_prop[A][B][mu][nu] += "+"
                        lorentz_prop[A][B][mu][nu] += "(" + imag_var + "/3.)" + "*" + "prop"
                    elif gamma_prod[A][B][mu][nu] == complex(0., -1.):
                        lorentz_prop[A][B][mu][nu] += "-" + "(" + imag_var + "/3.)" + "*" + "prop"
                    else:
                        if lorentz_prop[A][B][mu][nu] != "":
                            lorentz_prop[A][B][mu][nu] += "+"
                        lorentz_prop[A][B][mu][nu] += "(1./3.)" + "*" + "prop" + "*" + gamma_prod[A][B][mu][nu]

                if pmu_gammanu[A][B][mu][nu] != "":
                    lorentz_prop[A][B][mu][nu] += "-(1./(3.*mRS))" + "*" + "prop" + "*" + pmu_gammanu[A][B][mu][nu]

                if pnu_gammamu[A][B][mu][nu] != "":
                    if lorentz_prop[A][B][mu][nu] != "":
                        lorentz_prop[A][B][mu][nu] += "+"
                    lorentz_prop[A][B][mu][nu] += "(1./(3.*mRS))" + "*" + "prop" + "*" + pnu_gammamu[A][B][mu][nu]

                if p_term[A][B][mu][nu] != "":
                    if lorentz_prop[A][B][mu][nu] != "":
                        lorentz_prop[A][B][mu][nu] += "+"
                    lorentz_prop[A][B][mu][nu] += "(2./(3.*mRS2))" + "*" + "prop" + "*" + p_term[A][B][mu][nu]

# complete A-independent part of the propagator numerator
result = [[[['' for _ in range(4)] for _ in range(4)] for _ in range(4)] for _ in range(4)]
for A in range(4):
    for mu in range(4):
        for B in range(4):
            for nu in range(4):
                for C in range(4):
                    if spropagator[A][C] != "" and lorentz_prop[C][B][mu][nu] != "":
                        if result[A][B][mu][nu] != "":
                            result[A][B][mu][nu] += " + "
                        result[A][B][mu][nu] += "(" + spropagator[A][C] + ")" + "*" + "(" + lorentz_prop[C][B][mu][nu] \
                                                + ")"

# add A-dependent terms
prefactor = "(1./(3.*mRS2))" + "*" + "((1.+m_A)/(1.+2.*m_A))"
A_terms = [[[['' for _ in range(4)] for _ in range(4)] for _ in range(4)] for _ in range(4)]
for A in range(4):
    for mu in range(4):
        for B in range(4):
            for nu in range(4):
                if gamma_prod[A][B][mu][nu] != complex(0.0):
                    prefactor2 = "(m_A/(1.+2.*m_A))" + "*" + "mRS"
                    if gamma_prod[A][B][mu][nu] == complex(1., 0.):
                        if A_terms[A][B][mu][nu] != "":
                            A_terms[A][B][mu][nu] += "+"
                        A_terms[A][B][mu][nu] += prefactor + "*" + prefactor2
                    elif gamma_prod[A][B][mu][nu] == complex(-1., 0.):
                        A_terms[A][B][mu][nu] += "-" + prefactor + "*" + prefactor2
                    elif gamma_prod[A][B][mu][nu] == complex(0., 1.):
                        if A_terms[A][B][mu][nu] != "":
                            A_terms[A][B][mu][nu] += "+"
                        A_terms[A][B][mu][nu] += imag_var + "*" + prefactor + "*" + prefactor2
                    elif gamma_prod[A][B][mu][nu] == complex(0., -1.):
                        A_terms[A][B][mu][nu] += "-" + imag_var + "*" + prefactor + "*" + prefactor2
                    else:
                        if A_terms[A][B][mu][nu] != "":
                            A_terms[A][B][mu][nu] += "+"
                        A_terms[A][B][mu][nu] += prefactor + "*" + prefactor2 + "*" + gamma_prod[A][B][mu][nu]

                    prefactor3 = "((1.+m_A)/(2.+4.*m_A))"
                    if pdagger_gammamu_gammanu[A][B][mu][nu] != "":
                        A_terms[A][B][mu][nu] += "-" + prefactor + "*" + prefactor3 + "*" + \
                                                 pdagger_gammamu_gammanu[A][B][mu][nu]

                    if pnu_gammamu[A][B][mu][nu] != "":
                        A_terms[A][B][mu][nu] += "-" + prefactor + "*" + pnu_gammamu[A][B][mu][nu]

                    if pnu_gammamu[A][B][mu][nu] != "":
                        A_terms[A][B][mu][nu] += "-" + prefactor + "*" + prefactor2 + "*" + pmu_gammanu[A][B][mu][nu]

# Add A-dependent and A-independent terms
for A in range(4):
    for mu in range(4):
        for B in range(4):
            for nu in range(4):
                if A_terms[A][B][mu][nu] != "":
                    if result[A][B][mu][nu] != "":
                        result[A][B][mu][nu] += " + "
                    result[A][B][mu][nu] += A_terms[A][B][mu][nu]

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

# Write result into text file, can be directly copied to AddPropagator function
datei = open('textdatei.txt','w')

for mu in range(4):
    for A in range(4):
        # select, whether b>0 or b<0 should be printed out
        if b<0:
            result_print = result_bneg[A][mu]
        else:
            result_print = result_bpos[A][mu]
        # replace sums of momentum components
        result_print=result_print.replace("(p[0]+p[" + str(R3) + "])", "pp")
        result_print=result_print.replace("(p[0]-p[" + str(R3) + "])", "pm")
        result_print=result_print.replace("(p[" + str(R1) + "]+M_I*p[" + str(R2) + "])", "pt")
        result_print=result_print.replace("(M_I*p[" + str(R2) + "]+p[" + str(R1) + "])", "pt")
        result_print=result_print.replace("(p[" + str(R1) + "]-M_I*p[" + str(R2) + "])", "ptc")
        result_print=result_print.replace("(-M_I*p[" + str(R2) + "]+p[" + str(R1) + "])", "ptc")
        result_print=result_print.replace("(-p[" + str(R1) + "]-M_I*p[" + str(R2) + "])", "(-pt)")
        result_print=result_print.replace("(-M_I*p[" + str(R2) + "]-p[" + str(R1) + "])", "(-pt)")
        result_print=result_print.replace("(-p[" + str(R1) + "]+M_I*p[" + str(R2) + "])", "(-ptc)")
        result_print=result_print.replace("(M_I*p[" + str(R2) + "]-p[" + str(R1) + "])", "(-ptc)")
        # simplify expression further
        result_print=result_print.replace("(-(-ptc))", "ptc")
        result_print=result_print.replace("(-(-pt))", "pt")
        result_print=result_print.replace("(-M_I*(-ptc))", "M_I*ptc")
        result_print=result_print.replace("(-M_I*(-ptc))", "M_I*ptc")
        result_print=result_print.replace("-1./(3.*m2)*(-M_I", "+1./(3.*m2)*(M_I")
        result_print=result_print.replace("-1./(3.*m2)*(-pt)*", "+1./(3.*m2)*pt*")
        result_print=result_print.replace("-1./(3.*m2)*(-ptc)*", "+1./(3.*m2)*ptc*")
        result_print=result_print.replace("-1./(3.*m2)*(-pm)*", "+1./(3.*m2)*pm*")
        result_print=result_print.replace("-1./(3.*m2)*(-pp)*", "+1./(3.*m2)*pp*")

        datei.write("j[" + str(A+4*mu) + "]" +"=" + result_print + ";" + "\n")

# end of script
*/
