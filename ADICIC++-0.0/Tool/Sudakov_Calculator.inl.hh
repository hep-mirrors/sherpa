//bof
//Version: 1 ADICIC++-0.0/2004/07/12

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


  inline Sudakov_Calculator::~Sudakov_Calculator() {
    --s_count;
  }



  //---------------------------------------------------------------------------



  inline const bool Sudakov_Calculator::IsAlphaSRunning() {    //Static.
    return (s_isalphasrun && s_pas);
  }
  inline const double Sudakov_Calculator::MinOfK2t() {    //Static.
    return s_k2tmin;
  }
  inline const double Sudakov_Calculator::MaxOfK2t() {    //Static.
    return s_k2tmax;
  }
  inline const double Sudakov_Calculator::AlphaSFix() {    //Static.
    return s_alphasfix;
  }
  inline const double Sudakov_Calculator::AlphaSApprox() {    //Static.
    return s_approx;
  }
  inline const double Sudakov_Calculator::AlphaSCorr(const double p2t) {
    //Static method.
    return GetAlphaSCorr(p2t);
  }



  //---------------------------------------------------------------------------



  inline const double Sudakov_Calculator::FixAlphaSCorr(const double p2t) {
    //Static method.
    return 1.0;
  }
  inline const double Sudakov_Calculator::RunAlphaSCorr(const double p2t) {
    //Static method.
    //return (*s_pas)(p2t);    //Testing.
    //return (*s_pas)(p2t)/s_approx;
    double ret=(*s_pas)(p2t)/s_approx; assert(ret<=1.0);
    return ret;
  }



  //---------------------------------------------------------------------------



  inline void Sudakov_Calculator::GetResult(bool& gsplit, double& p2t,
					    double& x1, double& x3) const {
    gsplit = f_gsplit;
    p2t    = m_p2t;
    x1     = m_x1;
    x3     = m_x3;
  }


  inline void Sudakov_Calculator::Which() const {
    std::cout<<"Incomplete Sudakov_Calculator."<<std::endl;
  }



  //===========================================================================



  template<> inline void Sudakov<Dipole::qqbar>::Which() const {
    std::cout<<"Sudakov_Calculator for a q-qbar dipole using ";
    if(IsAlphaSRunning()) std::cout<<"running alpha_s.\n";
    else std::cout<<"fixed alpha_s.\n";
  }
  template<> inline void Sudakov<Dipole::qg>::Which() const {
    std::cout<<"Sudakov_Calculator for a q-g dipole using ";
    if(IsAlphaSRunning()) std::cout<<"running alpha_s.\n";
    else std::cout<<"fixed alpha_s.\n";
  }
  template<> inline void Sudakov<Dipole::gqbar>::Which() const {
    std::cout<<"Sudakov_Calculator for a g-qbar dipole using ";
    if(IsAlphaSRunning()) std::cout<<"running alpha_s.\n";
    else std::cout<<"fixed alpha_s.\n";
  }
  template<> inline void Sudakov<Dipole::gg>::Which() const {
    std::cout<<"Sudakov_Calculator for a g-g dipole using ";
    if(IsAlphaSRunning()) std::cout<<"running alpha_s.\n";
    else std::cout<<"fixed alpha_s.\n";
  }


  template<Dipole::Type DT>
  inline void Sudakov<DT>::ShowSpecification() const {
    std::cout<<"Sudakov handling for Dipole::Type: "
	     <<Sudakov_Info<DT>::Dipoletype
	     <<"\t ME correction powers: "
	     <<Sudakov_Info<DT>::X1power<<","
	     <<Sudakov_Info<DT>::X3power
	     <<"\t Colour factor: "<<Sudakov_Info<DT>::Colourfactor
	     <<std::endl;
  }



  //---------------------------------------------------------------------------



  template<Dipole::Type DT>
  inline void Sudakov<DT>::InitWith(const Dipole& dip) {

    //m_alphas is set with the alphas_fixed value.

    assert(dip.IsType()==Sudakov_Info<DT>::Dipoletype);

    m_s=dip.InvMass();
    m_x2tmin=Sudakov_Calculator::MinOfK2t()/m_s;
    //m_x2t=Min(1.0,Sudakov_Calculator::MaxOfK2t()/m_s);
    //m_x2t=Min(dip.ProdScale()/m_s,Sudakov_Calculator::MaxOfK2t()/m_s);
    m_x2t=Min(dip.BootScale()/m_s,Sudakov_Calculator::MaxOfK2t()/m_s);
    m_rap=0.0;
    m_corr=1.0;

  }





  template<Dipole::Type DT>
  inline void Sudakov<DT>::GenerateRap() {
    double ymax=-0.5*std::log(m_x2t);
    m_rap=ymax*(-1.0+2.0*ATOOLS::ran.Get());
  }





  template<Dipole::Type DT>
  inline void Sudakov<DT>::GenerateCorr() {
    m_corr=AlphaSCorr(m_p2t)*0.5*
           ( power<Sudakov_Info<DT>::X1power >(m_x1) +
	     power<Sudakov_Info<DT>::X3power >(m_x3));
  }





  template<Dipole::Type DT>
  inline const bool Sudakov<DT>::TestEfracs() const {
    double sum=m_x1+m_x3;
    if(sum>1.0 && sum<2.0) return true;
    return false;
  }



  //===========================================================================



}    //eo namespace ADICIC





//eof
