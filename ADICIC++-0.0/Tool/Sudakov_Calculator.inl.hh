//bof
//Version: 1 ADICIC++-0.0/2004/05/07

//Inline methods of Sudakov_Calculator.H.





#include <cassert>
#include <cstdlib>
#include <mathextra>
#include "Random.H"





namespace ADICIC {



  //===========================================================================



  inline Sudakov_Calculator::Sudakov_Calculator()
    : f_gsplit(false), m_p2t(0.0), m_x1(1.0), m_x3(1.0) {
    ++s_count;
  }


  Sudakov_Calculator::~Sudakov_Calculator() {
    --s_count;
  }



  //---------------------------------------------------------------------------



  inline const bool Sudakov_Calculator::IsAlphaSRunning() {
    return sf_alphasrun;
  }
  inline const double Sudakov_Calculator::MinOfK2t() {    //Static.
    return s_k2tmin;
  }
  inline const double Sudakov_Calculator::MaxOfK2t() {    //Static.
    return s_k2tmax;
  }
  inline const double Sudakov_Calculator::FixedAlphaS() {    //Static.
    return s_alphasfix;
  }



  //---------------------------------------------------------------------------



  inline void Sudakov_Calculator::GetResult(bool& gsplit, double& p2t,
					    double& x1, double& x3) const {
    gsplit = f_gsplit;
    p2t    = m_p2t;
    x1     = m_x1;
    x3     = m_x3;
  }


  void Sudakov_Calculator::Which() const {
    std::cout<<"Incomplete Sudakov_Calculator."<<std::endl;
  }



  //===========================================================================



  template<Dipole::Type DT, class AS> inline Sudakov<DT,AS>::Sudakov()
    : Sudakov_Calculator(),
      m_alphas(Sudakov_Calculator::FixedAlphaS()),
      m_s(Sudakov_Calculator::MaxOfK2t()),
      m_x2tmin(Sudakov_Calculator::MinOfK2t()/m_s), m_x2t(1.0),
      m_rap(0.0), m_corr(1.0) {}



  //---------------------------------------------------------------------------



  template<> inline void Sudakov<Dipole::qqbar,Alpha_S_Fix>::Which() const {
    std::cout<<"Sudakov_Calculator for a q-qbar dipole using fixed alpha_s.\n";
  }
  template<> inline void Sudakov<Dipole::qg,Alpha_S_Fix>::Which() const {
    std::cout<<"Sudakov_Calculator for a q-g dipole using fixed alpha_s.\n";
  }
  template<> inline void Sudakov<Dipole::gqbar,Alpha_S_Fix>::Which() const {
    std::cout<<"Sudakov_Calculator for a g-qbar dipole using fixed alpha_s.\n";
  }
  template<> inline void Sudakov<Dipole::gg,Alpha_S_Fix>::Which() const {
    std::cout<<"Sudakov_Calculator for a g-g dipole using fixed alpha_s.\n";
  }


  template<Dipole::Type DT, class AS>
  inline void Sudakov<DT,AS>::ShowSpecification() const {
    std::cout<<"Sudakov handling for Dipole::Type: "
	     <<Sudakov_Info<DT,AS>::Dipoletype
	     <<"\t ME correction powers: "
	     <<Sudakov_Info<DT,AS>::X1power<<","
	     <<Sudakov_Info<DT,AS>::X3power
	     <<"\t Colour factor: "<<Sudakov_Info<DT,AS>::Colourfactor
	     <<std::endl;
  }



  //---------------------------------------------------------------------------



  template<Dipole::Type DT, class AS>
  inline void Sudakov<DT,AS>::InitWith(const Dipole& dip) {

    //m_alphas is set with the alphas_fixed value.

    Dipole::Type test=Sudakov_Info<DT,AS>::Dipoletype;
    assert(dip.IsType()==test);

    m_s=dip.InvMass();
    m_x2tmin=Sudakov_Calculator::MinOfK2t()/m_s;
    m_x2t=Min(1.0,Sudakov_Calculator::MaxOfK2t()/m_s);
    m_rap=0.0;
    m_corr=1.0;

  }





  template<Dipole::Type DT, class AS>
  inline void Sudakov<DT,AS>::GenerateRap() {
    double ymax=-0.5*std::log(m_x2t);
    m_rap=ymax*(-1.0+2.0*ATOOLS::ran.Get());
  }





  template<Dipole::Type DT, class AS>
  inline void Sudakov<DT,AS>::GenerateCorr() {
    m_corr=0.5*( power<Sudakov_Info<DT,AS>::X1power >(m_x1) +
		 power<Sudakov_Info<DT,AS>::X3power >(m_x3));
  }





  template<Dipole::Type DT, class AS>
  inline const bool Sudakov<DT,AS>::TestEfracs() const {
    double sum=m_x1+m_x3;
    if(sum>1.0 && sum<2.0) return true;
    return false;
  }



  //===========================================================================



}    //eo namespace ADICIC





//eof
