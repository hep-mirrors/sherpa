//bof
//Version: 4 ADICIC++-0.0/2006/07/28

//Inline methods of Sudakov_Calculator.H.





#include <cassert>
#include <cstdlib>
#include "Dipole_Parameter.H"





namespace ADICIC {



  //===========================================================================



  inline Sudakov_Calculator::Sudakov_Calculator()
    : p_dip(NULL), p_sur(NULL) {
    //There is only this constructor, so that the following is sufficient.
    static int ini=bool(AdjustEnvironment() || true);
    //assert(s_count);
    assert(ini==1);////////////////////////////////////////////////////////////
    s_count+=ini;
  }



  //---------------------------------------------------------------------------



  inline const bool Sudakov_Calculator::IsAlphaSRunning() {    //Static.
    //return (dpa.sud.RunAlphaS() && s_box.m_ras[0]);
    return bool(s_box.m_ras[0]);
  }
  inline const double Sudakov_Calculator::AlphaSApprox() {    //Static.
    return s_asapprox;
  }
  inline const double Sudakov_Calculator::AlphaSCorr(const double p2t) {
    //Static method.
    return GetAlphaSCorr(p2t);
  }
  inline const int Sudakov_Calculator::Nf(const double p2t) {
    //Static method.
    return GetNf(p2t);
  }
  inline const bool Sudakov_Calculator::ArePDFsInitialized() {    //Static.
    //return (bool(sf_pdf) && s_box.m_pdf[0] && s_box.m_pdf[1]);
    return bool(s_box.m_pdf[0] && s_box.m_pdf[1]);
  }
  inline const double
  Sudakov_Calculator::PlusPDFCorr(const Multiflavour& mufl,
				  const Sudakov_Result& sure) {
    //Static method.
    return GetPDFCorr(0,mufl,sure);
  }
  inline const double
  Sudakov_Calculator::MinusPDFCorr(const Multiflavour& mufl,
				   const Sudakov_Result& sure) {
    //Static method.
    return GetPDFCorr(1,mufl,sure);
  }





  inline const Dipole& Sudakov_Calculator::CurrentDipole() const {
    assert(p_dip); return *p_dip;
  }
  inline
  const Sudakov_Result& Sudakov_Calculator::CurrentSudakovResult() const {
    assert(p_sur); return *p_sur;
  }



  inline const double
  Sudakov_Calculator::PlusPDFCorr(const Multiflavour& mufl) const {
    assert(p_sur); assert(p_sur->Isr.size()>=sr::stop);
    return GetPDFCorr(0,mufl,*p_sur);
  }
  inline const double
  Sudakov_Calculator::MinusPDFCorr(const Multiflavour& mufl) const {
    assert(p_sur); assert(p_sur->Isr.size()>=sr::stop);
    return GetPDFCorr(1,mufl,*p_sur);
  }



  //---------------------------------------------------------------------------



  inline const double Sudakov_Calculator::FixAlphaSCorr(const double p2t) {
    //Static method.
    assert(dpa.sud.AlphaSFix()/s_asapprox==1.0);
    return 1.0;
  }
  inline const double Sudakov_Calculator::RunAlphaSCorr(const double p2t) {
    //Static method.
    //return (*s_box.m_ras[0])(p2t);    //Testing.
    //return (*s_box.m_ras[0])(p2t)/s_asapprox;
    assert(p2t>dpa.sud.MinK2t() || p2t>dpa.sud.MinIIK2t());
    double ret=(*s_box.m_ras[0])(dpa.evo.RenoScaleFactor()*p2t)/s_asapprox;
    assert(ret<=1.0);
    return ret;
  }
  inline const int Sudakov_Calculator::FixNf(const double p2t) {
    //Static method.
    return dpa.sud.NfFix();
  }
  inline const int Sudakov_Calculator::RunNf(const double p2t) {
    //Static method.
    int ret=static_cast<MODEL::Running_AlphaS*>(s_box.m_ras[0])
      ->Nf(dpa.evo.RenoScaleFactor()*p2t);
    assert(ret>=-1);
    return ret;
  }



  //===========================================================================



  inline Sudakov_Base::Sudakov_Base(const Sudakov_Flavour& c,
				    const Radiation::Group r)
    : m_code(c), m_split(between), m_rg(r), m_step(-1), m_mass(0.0) {
    ++s_count;
  }


  inline const Sudakov_Flavour& Sudakov_Base::RadCode() const {
    return m_code;
  }
  inline const xbool Sudakov_Base::RadPart() const {
    return m_split;
  }
  inline const Radiation::Group Sudakov_Base::RadGroup() const {
    return m_rg;
  }
  inline const double Sudakov_Base::RadMass() const {
    return m_mass;
  }


  inline void Sudakov_Base::InitRadParticle() {
    if(m_code.Gluic()) { m_mass=(*m_code.Glu)().Mass(); return;}
    if(m_code.Bqqic()) { m_mass=(*m_code.Qua)().Mass(); return;}
    assert(m_code.Quaic() || m_code.Aquic());
  }



  //===========================================================================



}    //eo namespace ADICIC





//eof
