//bof
//Version: 2 ADICIC++-0.0/2004/10/28

//Inline methods of Sudakov_Calculator.H.





#include <cassert>
#include <cstdlib>
#include <mathextra>
#include "Random.H"





namespace ADICIC {



  //===========================================================================



  inline Sudakov_Calculator::Sudakov_Calculator()
    : m_spl(between), m_rad(ATOOLS::kf::none),
      m_p2t(0.0), m_x1(1.0), m_x3(1.0) {
    ++s_count;
  }



  //---------------------------------------------------------------------------



  inline const bool Sudakov_Calculator::IsAlphaSRunning() {    //Static.
    return (s_isalphasrun && s_pas);
  }
  inline const int Sudakov_Calculator::NfFix() {    //Static.
    return s_nffix;
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
  inline const int Sudakov_Calculator::Nf(const double p2t) {
    //Static method.
    return GetNf(p2t);
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
  inline const int Sudakov_Calculator::FixNf(const double p2t) {
    //Static method.
    return s_nffix;
  }
  inline const int Sudakov_Calculator::RunNf(const double p2t) {
    //Static method.
    int ret=static_cast<MODEL::Running_AlphaS*>(s_pas)->Nf(p2t); assert(ret>0);
    return ret;
  }



  //---------------------------------------------------------------------------



  inline void Sudakov_Calculator::GetResult(Sudakov_Result& sudre) const {
    using namespace ATOOLS;
    switch(m_spl) {
    case between:
      sudre.Rad=Radiation::gluon; assert(m_rad==kf::gluon); break;
    case front  :
      sudre.Rad=Radiation::qtop;
      assert(m_rad==kf::d || m_rad==kf::u || m_rad==kf::s ||
	     m_rad==kf::c || m_rad==kf::b);
      break;
    case end    :
      sudre.Rad=Radiation::qbot;
      assert(m_rad==kf::d || m_rad==kf::u || m_rad==kf::s ||
	     m_rad==kf::c || m_rad==kf::b);
      break;
    default     :
      assert(0);
    }
    sudre.Kfc=m_rad;
    sudre.P2t=m_p2t;
    sudre.X1 =m_x1;
    sudre.X3 =m_x3;
  }



  //===========================================================================



  inline Sudakov_Base::Sudakov_Base(ATOOLS::kf::code c)
    : m_code(c), m_split(between), m_mass(0.0) {
    ++s_count;
  }


  inline const ATOOLS::kf::code Sudakov_Base::RadCode() const {
    return m_code;
  }
  inline const xbool Sudakov_Base::RadPart() const {
    return m_split;
  }
  inline const double Sudakov_Base::RadMass() const {
    return m_mass;
  }


  inline void Sudakov_Base::InitRadParticle() {
    m_mass=ATOOLS::Flavour(m_code).Mass();
  }



  //===========================================================================



  template<> inline void Sudakov_Group<Dipole::qqbar>::Which() const {
    std::cout<<"Sudakov_Calculator for a q-qbar dipole using ";
    if(IsAlphaSRunning()) std::cout<<"running alpha_s.\n";
    else std::cout<<"fixed alpha_s.\n";
  }
  template<> inline void Sudakov_Group<Dipole::qg>::Which() const {
    std::cout<<"Sudakov_Calculator for a q-g dipole using ";
    if(IsAlphaSRunning()) std::cout<<"running alpha_s.\n";
    else std::cout<<"fixed alpha_s.\n";
  }
  template<> inline void Sudakov_Group<Dipole::gqbar>::Which() const {
    std::cout<<"Sudakov_Calculator for a g-qbar dipole using ";
    if(IsAlphaSRunning()) std::cout<<"running alpha_s.\n";
    else std::cout<<"fixed alpha_s.\n";
  }
  template<> inline void Sudakov_Group<Dipole::gg>::Which() const {
    std::cout<<"Sudakov_Calculator for a g-g dipole using ";
    if(IsAlphaSRunning()) std::cout<<"running alpha_s.\n";
    else std::cout<<"fixed alpha_s.\n";
  }


  template<Dipole::Type DT>
  inline void Sudakov_Group<DT>::ShowSpecification() const {
    std::cout<<"Sudakov_Group for dipole and radiation type: "
	     <<Sudakov_Info<DT,Radiation::gluon>::Dipoletype<<" and "
	     <<m_radtype<<".\n";
    if(l_sud.empty()) { std::cout<<" No Sudakov's assigned.\n"; return;}
    short j=1;
    for(std::list<Sudakov_Base*>::const_iterator bit=l_sud.begin();
	bit!=l_sud.end(); ++bit) {
      std::cout<<" "<<j<<"."; (*bit)->Which();
      std::cout<<" "<<j<<"."; (*bit)->ShowSpecification();
      ++j;
    }
  }



  //===========================================================================



  template<>
  inline void Sudakov<Dipole::qqbar,Radiation::gluon>::Which() const {
    std::cout<<"Sudakov_Base for: q-qbar dipole producing g emission.\n";
  }
  template<>
  inline void Sudakov<Dipole::qg,Radiation::gluon>::Which() const {
    std::cout<<"Sudakov_Base for: q-g dipole producing g emission.\n";
  }
  template<>
  inline void Sudakov<Dipole::gqbar,Radiation::gluon>::Which() const {
    std::cout<<"Sudakov_Base for: g-qbar dipole producing g emission.\n";
  }
  template<>
  inline void Sudakov<Dipole::gg,Radiation::gluon>::Which() const {
    std::cout<<"Sudakov_Base for: g-g dipole producing g emission.\n";
  }


  template<>
  inline void Sudakov<Dipole::qg,Radiation::quark>::Which() const {
    std::cout<<"Sudakov_Base for: q-g dipole producing q-qbar splitting.\n";
  }
  template<>
  inline void Sudakov<Dipole::gqbar,Radiation::quark>::Which() const {
    std::cout<<"Sudakov_Base for: g-qbar dipole producing q-qbar splitting.\n";
  }
  template<>
  inline void Sudakov<Dipole::gg,Radiation::quark>::Which() const {
    std::cout<<"Sudakov_Base for: g-g dipole producing q-qbar splitting.\n";
  }





  template<Dipole::Type DT>
  inline void Sudakov<DT,Radiation::gluon>::ShowSpecification() const {
    std::cout<<"Sudakov for: dipole type="
	     <<Sudakov_Info<DT,Radiation::gluon>::Dipoletype
	     <<",   radiation group="
	     <<Sudakov_Info<DT,Radiation::gluon>::Radiationgroup
	     <<" and flavour code="<<m_code<<" ("<<m_mass<<" GeV)"
	     <<",   ME correction powers="
	     <<Sudakov_Info<DT,Radiation::gluon>::X1power<<","
	     <<Sudakov_Info<DT,Radiation::gluon>::X3power
	     <<",   colour factor="
	     <<Sudakov_Info<DT,Radiation::gluon>::Colourfactor
	     <<".\n";
  }


  template<Dipole::Type DT>
  inline void Sudakov<DT,Radiation::quark>::ShowSpecification() const {
    std::cout<<"Sudakov for: dipole type="
	     <<Sudakov_Info<DT,Radiation::quark>::Dipoletype
	     <<",   radiation group="
	     <<Sudakov_Info<DT,Radiation::quark>::Radiationgroup
	     <<" and flavour code="<<m_code<<" ("<<m_mass<<" GeV)"
	     <<",   colour factor="
	     <<Sudakov_Info<DT,Radiation::quark>::Colourfactor
	     <<".\n";
  }



  //===========================================================================



  template<Dipole::Type DT>
  inline const double Sudakov_Group<DT>::Sdip() const {
    return m_s;
  }

  template<Dipole::Type DT>
  inline const double Sudakov_Group<DT>::X2tmin() const {
    return m_x2tmin;
  }

  template<Dipole::Type DT>
  inline const double Sudakov_Group<DT>::X2t() const {
    return m_x2t;
  }

  template<Dipole::Type DT>
  inline const double Sudakov_Group<DT>::Ymax() const {
    return m_ymax;
  }





  template<Dipole::Type DT>
  inline const bool Sudakov_Group<DT>::TestEfracs(const double x1,
						  const double x3) const {
    double sum=x1+x3;
    if(sum>1.0 && sum<2.0) return true;
    return false;
  }

  template<Dipole::Type DT>
  inline void Sudakov_Group<DT>::Reset() {
    m_x2t=m_x2tmax;
    m_ymax=0.0;
    m_rap=0.0;
    m_corr=1.0;
  }



  //===========================================================================



  //===============
  //Gluon emission.
  //===============


  template<Dipole::Type DT>
  inline const double Sudakov<DT,Radiation::gluon>::CalculateRapLimit() {
    assert(m_step==0); ++m_step;
    //return -0.5*std::log(m_sgroup.X2t());    //Old.
    double ratio=0.25/m_sgroup.X2t();
    return std::log( sqrt(ratio)+sqrt(ratio-1.0) );
  }


  template<Dipole::Type DT>
  inline const double Sudakov<DT,Radiation::gluon>::GenerateRap() {
    assert(m_step==1); ++m_step;
    return m_sgroup.Ymax()*(-1.0+2.0*ATOOLS::ran.Get());
  }


  template<Dipole::Type DT>
  inline const double Sudakov<DT,Radiation::gluon>::
  GenerateCorr(const double xa, const double xb) {
    assert(m_step==2); ++m_step;
    double x2t=m_sgroup.X2t();
    return
      Sudakov_Calculator::AlphaSCorr(m_sgroup.Sdip()*x2t) *    //0.5 *
      ( -m_sgroup.Ymax()/std::log(x2t) ) *    //Old: exchange this by 0.5.
      ( power<Sudakov_Info<DT,Radiation::gluon>::X1power >(xa) +
	power<Sudakov_Info<DT,Radiation::gluon>::X3power >(xb) );
  }





  //==========================
  //Quark-antiquark splitting.
  //==========================


  template<Dipole::Type DT>
  inline const double Sudakov<DT,Radiation::quark>::CalculateRapLimit() {
    assert(m_step==0); ++m_step;
    double ratio=0.25/m_sgroup.X2t();
    return sqrt(ratio)+sqrt(ratio-1.0);
  }


  template<Dipole::Type DT>
  inline const double Sudakov<DT,Radiation::quark>::GenerateRap() {
    assert(m_step==1); ++m_step;
    double zmaxr=1.0/m_sgroup.Ymax();
    return m_split*std::log(zmaxr+(1.0/zmaxr-zmaxr)*ATOOLS::ran.Get());
  }


  template<Dipole::Type DT>
  inline const double Sudakov<DT,Radiation::quark>::
  GenerateCorr(const double xa, const double xb) {
    assert(m_step==2); ++m_step;
    double x2t=m_sgroup.X2t();
    double zmax=m_sgroup.Ymax();
    double corr=
      Sudakov_Calculator::AlphaSCorr(m_sgroup.Sdip()*x2t) *
      ( ATOOLS::sqr(xa+xb-1.0)+bool(m_split-1)*ATOOLS::sqr(1-xb)+bool(m_split+1)*ATOOLS::sqr(1-xa) ) *
      sqrt(x2t)*(zmax-1.0/zmax);
    if(Sudakov_Calculator::Ariadne) corr/=(2.0*sqrt(m_sgroup.Sdip()));
    return corr;
  }


  template<Dipole::Type DT>
  inline void Sudakov<DT,Radiation::quark>::SetSplitBranch() {
#ifdef SUDAKOV_CALCULATOR_OUTPUT
    std::cout<<"    xxx> "<<m_split<<", "<<bool(m_split+1)<<"\n";
#endif
  }
  template<>
  inline void Sudakov<Dipole::gg,Radiation::quark>::SetSplitBranch() {
    m_split=xbool(-1+2*bool(ATOOLS::ran.Get()>0.5));
    assert(m_split!=between);//////////////////////////////////////////////////
#ifdef SUDAKOV_CALCULATOR_OUTPUT
    std::cout<<"    xyx> "<<m_split<<", "<<bool(m_split+1)<<"\n";
#endif
  }



  //===========================================================================



}    //eo namespace ADICIC





//eof
