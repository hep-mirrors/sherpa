//bof
//Version: 1 ADICIC++-0.0/2004/05/10

//Inline methods of Recoil_Calculator.H.





#include <cassert>
#include <cstdlib>
#include "Random.H"
#include "Poincare.H"





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



  template<class ST> inline Recoil<ST>::Recoil()
    : Recoil_Calculator(),
      m_costheta(-1.0), m_sintheta(0.0), m_phi(0.0),
      p_ini(NULL), m_cmsaxis(Recoil_Calculator::ZAxis) {}



  //---------------------------------------------------------------------------



  template<class ST> inline void Recoil<ST>::Which() const {
    std::cout<<"This is a default Recoil_Calculator!\n";
  }
  template<> inline void Recoil<Kleiss_Strategy>::Which() const {
    std::cout<<"Recoil_Calculator implementing the Kleiss prescription.\n";
  }
  template<> inline void Recoil<FixDir1_Strategy>::Which() const {
    std::cout<<"Recoil_Calculator implementing the fixed-direction-1 idea.\n";
  }
  template<> inline void Recoil<FixDir3_Strategy>::Which() const {
    std::cout<<"Recoil_Calculator implementing the fixed-direction-3 idea.\n";
  }
  template<> inline void Recoil<MinimizePt_Strategy>::Which() const {
    std::cout<<"Recoil_Calculator implementing the original Pt minimization "
	     <<"strategy.\n";
  }
  template<> inline void Recoil<Lonnblad_Strategy>::Which() const {
    std::cout<<"Recoil_Calculator implementing a Lonnblad idea.\n";
  }
  template<> inline void Recoil<OldAdicic_Strategy>::Which() const {
    std::cout<<"Recoil_Calculator implementing the old Adicic approach.\n";
  }
  template<> inline void Recoil<Test_Strategy>::Which() const {
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



  //---------------------------------------------------------------------------



  template<class ST>
  inline const bool Recoil<ST>::Initialize() {
    m_costheta=( sqr(p_ini->GetE2())-sqr(p_ini->GetE1())-sqr(p_ini->GetE3()) )/
               ( 2.0*p_ini->GetE1()*p_ini->GetE3() );
    if(dabs(m_costheta)>1.0) {
      cerr<<"\nError: `m_costheta' value is out of range!\n";
      assert(dabs(m_costheta)<=1.0);
      return false;
    }
    m_sintheta=sqrt(1.0-sqr(m_costheta));
    m_phi=2.0*M_PI*ATOOLS::ran.Get();
    return true;
  }





  template<class ST>
  inline void Recoil<ST>::RotateOnto(const ATOOLS::Vec4D& axis) {
    ATOOLS::Poincare rot(Recoil_Calculator::ZAxis,axis);
    rot.Rotate(m_p1);
    rot.Rotate(m_p3);
  }
  template<class ST>
  inline void Recoil<ST>::Rotate(const ATOOLS::Vec4D& axis,
				 const ATOOLS::Vec4D& newaxis) {
    ATOOLS::Poincare rot(axis,newaxis);
    rot.Rotate(m_p1);
    rot.Rotate(m_p3);
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
