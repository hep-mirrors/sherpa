//bof
//Version: 2 ADICIC++-0.0/2004/08/06

//Inline methods of Recoil_Calculator.H.





#include <cassert>
#include <cstdlib>
#include "Random.H"
#include "Poincare.H"
#include "Recoil_Strategy.hpp"





namespace ADICIC {



  //===========================================================================



  inline Recoil_Calculator::Recoil_Calculator()
    : f_recoil(Nil), m_p1(ATOOLS::Vec4D()), m_p3(ATOOLS::Vec4D()) {
    ++s_count;
  }


  Recoil_Calculator::~Recoil_Calculator() {
    --s_count;
  }



  //---------------------------------------------------------------------------



  inline void Recoil_Calculator::GetResult(Trio& recoilparticle,
					   ATOOLS::Vec4D& p1,
					   ATOOLS::Vec4D& p3) const {
    recoilparticle = f_recoil;
    p1             = m_p1;
    p3             = m_p3;
  }



  //===========================================================================



  template<class ST> inline void Recoil<ST>::Which() const {
    std::cout<<"This is a non-specified Recoil_Calculator!\n";
  }
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





  template<class ST>
  inline const bool Recoil<ST>::GenerateCmsMomenta(const Recoil_Setup& ini,
						   const ATOOLS::Vec4D& cax) {
    p_ini=&ini;
    m_cmsaxis=cax;

    if(Initialize()); else { p_ini=NULL; return false;}

    bool result=Calculate();

    p_ini=NULL;
    return result;

  }



  //===========================================================================



  inline const double Recoil_Setup::GetE1() const { return E1;}
  inline const double Recoil_Setup::GetE2() const { return E2;}
  inline const double Recoil_Setup::GetE3() const { return E3;}
  inline const double Recoil_Setup::GetQ1() const { return Q1;}
  inline const double Recoil_Setup::GetQ2() const { return Q2;}
  inline const double Recoil_Setup::GetQ3() const { return Q3;}



  //===========================================================================



}    //eo namespace ADICIC





//eof
