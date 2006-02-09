//bof
//Version: 4 ADICIC++-0.0/2006/02/09

//Inline methods of Sudakov_Group.H.





#include <cassert>
#include <cstdlib>
#include <mathextra>
#include "Random.H"





using ATOOLS::sqr;



namespace ADICIC {



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
	     <<DT<<" and "<<m_radtype<<".\n";
    if(l_sud.empty()) { std::cout<<" No Sudakov's assigned.\n"; return;}
    short j=1;
    for(std::list<Sudakov_Base*>::const_iterator bit=l_sud.begin();
	bit!=l_sud.end(); ++bit) {
      std::cout<<" "<<j<<"."; (*bit)->Which();
      std::cout<<" "<<j<<"."; (*bit)->ShowSpecification();
      ++j;
    }
  }



  //---------------------------------------------------------------------------



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



  //---------------------------------------------------------------------------



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
    std::cout<<"Sudakov for: dip.type="<<DT
	     <<",  rad.group="<<Radiation::gluon
	     <<" and flav.code="<<(*m_code.Glu)().Kfcode()
	     <<" ("<<m_mass<<" GeV)"
	     <<",  ME.powers="<<s_x1pow<<","<<s_x3pow
	     <<",  col.factor="<<s_colfac<<".\n";
  }


  template<Dipole::Type DT>
  inline void Sudakov<DT,Radiation::quark>::ShowSpecification() const {
    std::cout<<"Sudakov for: dip.type="<<DT
	     <<",  rad.group="<<Radiation::quark
	     <<" and flav.code="<<(*m_code.Qua)().Kfcode()
	     <<" ("<<m_mass<<" GeV)"
	     <<",  col.factor="<<s_colfac<<".\n";
  }



  //---------------------------------------------------------------------------



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
  GenerateCorr(const Multidouble& x) {
    assert(m_step==2); ++m_step;
    double x2t=m_sgroup.X2t();
    //Test the overestimation.
    //double mev=0.5*( power<s_x1pow>(x[0]) + power<s_x3pow>(x[1]) );
    //if(mev>0.97) std::cout<<"  -->  "<<mev<<std::endl;
    //assert(mev<=1.0);
    return
      Sudakov_Calculator::AlphaSCorr(m_sgroup.Sdip()*x2t) *    //0.5 *
      ( -m_sgroup.Ymax()/std::log(x2t) ) *    //Old: exchange this by 0.5.
      ( power(x[0],s_x1pow) + power(x[1],s_x3pow) );
      //( power<s_x1pow>(x[0]) + power<s_x3pow>(x[1]) );
  }





  //==========================
  //Quark-antiquark splitting.
  //==========================


  template<Dipole::Type DT>
  inline bool Sudakov<DT,Radiation::quark>::NfVeto(const double scl) const {
    assert(m_step==0);
    if((*m_code.Qua)().Kfcode() > Sudakov_Calculator::Nf(scl)) return true;
    return false;
  }


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
  GenerateCorr(const Multidouble& x) {
    assert(m_step==2); ++m_step;
    double x2t=m_sgroup.X2t();
    double zmax=m_sgroup.Ymax();
    double corr=
      Sudakov_Calculator::AlphaSCorr(m_sgroup.Sdip()*x2t) *
      ( sqr(x[0]+x[1]-1.0)+
	bool(m_split-1)*sqr(1-x[1])+
	bool(m_split+1)*sqr(1-x[0]) ) *
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
