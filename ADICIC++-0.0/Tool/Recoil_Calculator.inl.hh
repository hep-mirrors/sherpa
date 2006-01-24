//bof
//Version: 3 ADICIC++-0.0/2005/09/13

//Inline methods of Recoil_Calculator.H.





#include <cassert>
#include <cstdlib>
#include "Random.H"
#include "Poincare.H"





namespace ADICIC {



  //===========================================================================



  inline Recoil_Calculator::Recoil_Calculator()
    : p_dip(NULL), p_sur(NULL), p_rer(NULL) {
    ++s_count;
  }



  //===========================================================================



  template<Recoil_Strategy::Type ST>
  inline const Recoil_Strategy::Type Recoil<ST>::IsType() const {
    return ST;
  }





  template<Recoil_Strategy::Type ST> inline void Recoil<ST>::Which() const {
    std::cout<<"This is a non-specified Recoil_Calculator!\n";}

  template<> inline void Recoil<Recoil_Strategy::Kleiss>::Which() const {
    std::cout<<"Recoil_Calculator implementing the Kleiss prescription.\n";
  }
  template<> inline void Recoil<Recoil_Strategy::FixDir1>::Which() const {
    std::cout<<"Recoil_Calculator implementing the fixed-direction-1 idea.\n";
  }
  template<> inline void Recoil<Recoil_Strategy::FixDir3>::Which() const {
    std::cout<<"Recoil_Calculator implementing the fixed-direction-3 idea.\n";
  }
  template<> inline void Recoil<Recoil_Strategy::MinimizePt>::Which() const {
    std::cout<<"Recoil_Calculator implementing the original Pt minimization "
	     <<"strategy.\n";
  }
  template<> inline void Recoil<Recoil_Strategy::Lonnblad>::Which() const {
    std::cout<<"Recoil_Calculator implementing a Lonnblad idea.\n";
  }
  template<> inline void Recoil<Recoil_Strategy::OldAdicic>::Which() const {
    std::cout<<"Recoil_Calculator implementing the old Adicic approach.\n";
  }
  template<> inline void Recoil<Recoil_Strategy::Test>::Which() const {
    std::cout<<"Recoil_Calculator implementing a test strategy.\n";
  }
  template<> inline void Recoil<Recoil_Strategy::Ktii>::Which() const {
    std::cout<<"Recoil_Calculator implementing the hadronic lab frame gluon "
	     <<"kt strategy.\n";
  }





  template<Recoil_Strategy::Type ST>
  inline void Recoil<ST>::TestKey(Dipole_Handler::Key) const {}

  template<> inline void
  Recoil<Recoil_Strategy::Kleiss>::TestKey(Dipole_Handler::Key k) const {
    //std::cout<<k.first<<" "<<k.first/100<<"\n";
    assert(k.first/100==0);
  }
  template<> inline void
  Recoil<Recoil_Strategy::FixDir1>::TestKey(Dipole_Handler::Key k) const {
    assert(k.first/100==0);
  }
  template<> inline void
  Recoil<Recoil_Strategy::FixDir3>::TestKey(Dipole_Handler::Key k) const {
    assert(k.first/100==0);
  }
  template<> inline void
  Recoil<Recoil_Strategy::MinimizePt>::TestKey(Dipole_Handler::Key k) const {
    assert(k.first/100==0);
  }
  template<> inline void
  Recoil<Recoil_Strategy::Lonnblad>::TestKey(Dipole_Handler::Key k) const {
    assert(k.first/100==0);
  }
  template<> inline void
  Recoil<Recoil_Strategy::OldAdicic>::TestKey(Dipole_Handler::Key k) const {
    assert(k.first/100==0);
  }
  template<> inline void
  Recoil<Recoil_Strategy::Test>::TestKey(Dipole_Handler::Key k) const {
    assert(k.first/100==0);
  }
  template<> inline void
  Recoil<Recoil_Strategy::Ktii>::TestKey(Dipole_Handler::Key k) const {
    assert((k.first==Dipole::iiqbarq &&
	    (k.second==Radiation::gluon || k.second==Radiation::qfront ||
	     k.second==Radiation::qbarend))
	   ||
	   (k.first==Dipole::iiqbarg &&
	    (k.second==Radiation::gluon || k.second==Radiation::qfront))
	   ||
	   (k.first==Dipole::iigq &&
	    (k.second==Radiation::gluon || k.second==Radiation::qbarend))
	   ||
	   (k.first==Dipole::iigg && k.second==Radiation::gluon));
  }





  template<Recoil_Strategy::Type ST>
  inline const bool Recoil<ST>::GenerateMomenta(const Dipole& D,
						const Sudakov_Result& S,
						Recoil_Result& R) {
    p_dip=&D;
    p_sur=&S;
    p_rer=&R; assert(p_rer->Poc==both && p_rer->Vec.size()==rr::stop);
    m_e.clear();///////////////////////////////////////////////////////////////
    bool result=Generate();
    p_dip=NULL; p_sur=NULL; p_rer=NULL;
    return result;
  }



  //---------------------------------------------------------------------------



  template<Recoil_Strategy::Type ST>
  inline const bool Recoil<ST>::Generate() {
    cerr<<"\nMethod: "<<__PRETTY_FUNCTION__<<": "
	<<"Warning: Recoil strategy is unspecified!\n"<<endl;
    return false;
  }

  template<>
  inline const bool Recoil<Recoil_Strategy::Kleiss>::Generate() {
    return CmsMomenta();
  }
  template<>
  inline const bool Recoil<Recoil_Strategy::FixDir1>::Generate() {
    return CmsMomenta();
  }
  template<>
  inline const bool Recoil<Recoil_Strategy::FixDir3>::Generate() {
    return CmsMomenta();
  }
  template<>
  inline const bool Recoil<Recoil_Strategy::MinimizePt>::Generate() {
    return CmsMomenta();
  }
  template<>
  inline const bool Recoil<Recoil_Strategy::Lonnblad>::Generate() {
    return CmsMomenta();
  }
  template<>
  inline const bool Recoil<Recoil_Strategy::OldAdicic>::Generate() {
    return CmsMomenta();
  }
  template<>
  inline const bool Recoil<Recoil_Strategy::Test>::Generate() {
    return CmsMomenta();
  }
  template<>
  inline const bool Recoil<Recoil_Strategy::Ktii>::Generate() {
    return LabMomenta();
  }



  //===========================================================================



}    //eo namespace ADICIC





//eof
