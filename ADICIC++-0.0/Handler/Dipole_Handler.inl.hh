//bof
//Version: 2 ADICIC++-0.0/2004/09/03

//Inline methods of Dipole_Handler.H.





#include <cassert>
#include <cstdlib>





namespace ADICIC {



  //===========================================================================



  inline const bool4::level Dipole_Handler::Status() const {
    if(!p_dix && !p_ban && !p_ati && !p_glu) return bool4::zero;
    if( p_dix && !p_ban && !p_ati &&  p_glu) return bool4::one;
    if(!p_dix &&  p_ban &&  p_ati &&  p_glu) return bool4::two;
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
    if(p_dix || p_ban || p_ati || p_glu) return false;
    if(p_dip==NULL && pD->IsHandledBy(*this)) {
      p_dip=pD;
      m_key.first=p_dip->IsType();
      m_key.second=Radiation::incorrect;
      p_sudakov=s_sumap[m_key.first];
      p_recoil=s_remap[m_key];
      f_below=false;
      f_recoil=Nil;
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
      m_key=Key(Dipole::incorrect,Radiation::incorrect);
      //m_key=std::make_pair(Dipole::incorrect,Radiation::incorrect);
      p_dip=NULL;
      return true;
    }
    return false;
  }



  //===========================================================================



  inline void Dipole_Handler::DecoupleNew(Dipole*& pD,
					  Dipole::Glubranch*& pG,
					  Dipole::Antibranch*& pA,
					  Dipole::Branch*& pB,
					  bool& below,
					  Trio& recoil) {
    if(pD || pG || pA || pB) return;
    pD=p_dix; p_dix=NULL;
    pG=p_glu; p_glu=NULL;
    pA=p_ati; p_ati=NULL;
    pB=p_ban; p_ban=NULL;
    below=f_below; f_below=false;
    recoil=f_recoil; f_recoil=Nil;
  }



  //===========================================================================



}    //eo namespace ADICIC





//eof
