//bof
//Version: 1 ADICIC++-0.0/2004/06/01

//Inline methods of Dipole_Handler.H.





#include <cassert>
#include <cstdlib>





namespace ADICIC {



  //===========================================================================



  inline const bool4::level Dipole_Handler::Status() const {
    if(!p_dix && !p_ban && !p_ati && !p_glu) return bool4::zero;
    if( p_dix && !p_ban && !p_ati &&  p_glu) return bool4::one;
    if(!p_dix &&  p_ban &&  p_ati && !p_glu) return bool4::two;
    return bool4::three;
  }


  inline const int Dipole_Handler::IsDocked() const {
    if(p_dip) return p_dip->Name; return 0;
  }


  inline const bool Dipole_Handler::IsDockedAt(const Dipole& dip) const {
    if(p_dip==&dip) return true; return false;
  }


  inline const bool Dipole_Handler::IsWaiting() const {
    return bool(f_gate);
  }



  //===========================================================================



  inline const bool Dipole_Handler::AttachDipole(Dipole* pD) {
    if( p_dix || p_ban || p_ati || p_glu) return false;
    if( p_dip==NULL && pD->IsHandledBy(*this) ) {
      p_dip=pD;
      p_tempa=s_map[p_dip->IsType()];
      p_sudakov=p_tempa->p_sud;
      p_recoil=p_tempa->p_rec;
      f_below=false;
      f_gate=0;
      return true;
    }
    return false;
  }


  inline const bool Dipole_Handler::DetachDipole(const Dipole* pD) {
    if(p_dip==pD && pD->IsHandled()==false) {
      f_gate=0;
      p_recoil=NULL;
      p_sudakov=NULL;
      p_tempa=NULL;
      p_dip=NULL;
      return true;
    }
    return false;
  }



  //===========================================================================



  inline void Dipole_Handler::ResetStatus() {
    f_below=false;
    if(p_dix) { delete p_dix; p_dix=NULL;}
    if(p_ban) { delete p_ban; p_ban=NULL;}
    if(p_ati) { delete p_ati; p_ati=NULL;}
    if(p_glu) { delete p_glu; p_glu=NULL;}
  }


  inline void Dipole_Handler::DecoupleNewDipole(Dipole*& pD, bool& below) {
    if(pD || !p_dix) return; pD=p_dix; p_dix=NULL; below=f_below;
  }


  inline void Dipole_Handler::DecoupleGlubranch(Dipole::Glubranch*& pG) {
    if(pG || !p_glu) return; pG=p_glu; p_glu=NULL;
  }



  //===========================================================================



  //===========================================================================



}    //eo namespace ADICIC





//eof
