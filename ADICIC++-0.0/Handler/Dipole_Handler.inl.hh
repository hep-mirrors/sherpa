//bof
//Version: 4 ADICIC++-0.0/2006/05/24

//Inline methods of Dipole_Handler.H.





#include <cassert>
#include <cstdlib>





namespace ADICIC {



  //===========================================================================



  inline const char Dipole_Handler::Status() const {
    if(!p_dix && !p_ban && !p_ati && !p_glu) return 'n';
    if( p_dix && !p_ban && !p_ati &&  p_glu) return 'g';
    if(!p_dix &&  p_ban &&  p_ati &&  p_glu) return 's';
    if( p_dix &&  p_ban &&  p_ati &&  p_glu) return 'q';
    return '?';
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



  inline Dipole_Handler::Carrier::Carrier()
    : DipOrder(0), RecoComp(both),
      pQua(NULL), pAqu(NULL), pGlu(NULL), pDip(NULL),
      Mup() {}


  inline void Dipole_Handler::Carrier::Reset() {
    DipOrder=0; RecoComp=both;
    pQua=NULL; pAqu=NULL; pGlu=NULL; pDip=NULL;
    Mup.Clear();
  }



  //===========================================================================



}    //eo namespace ADICIC





//eof
