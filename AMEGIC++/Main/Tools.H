#ifndef AMEGIC__Main__Tools_H
#define AMEGIC__Main__Tools_H

#include "ATOOLS/Phys/Flavour.H"

namespace AMEGIC {

  inline bool IsSusy(const ATOOLS::Flavour &fl)
  { return 1000000<fl.Kfcode()&&fl.Kfcode()<3000000; }
  
  inline bool IsChargino(const ATOOLS::Flavour &fl)
  { return (fl.Kfcode()==kf_Chargino1||fl.Kfcode()==kf_Chargino2) && fl.IntSpin()==1; }
  inline bool IsNeutralino(const ATOOLS::Flavour &fl) 
  { return (fl.Kfcode()==kf_Neutralino1||fl.Kfcode()==kf_Neutralino2||
	    fl.Kfcode()==kf_Neutralino3||fl.Kfcode()==kf_Neutralino4) && fl.IntSpin()==1; }
  inline bool IsSlepton(const ATOOLS::Flavour &fl) 
  { return ((fl.Kfcode()>1000010&&fl.Kfcode()<1000017)||
	    (fl.Kfcode()>2000010&&fl.Kfcode()<2000017)) && fl.IntSpin()==0; }
  inline bool IsSneutrino(const ATOOLS::Flavour &fl) 
  { return fl.Kfcode()>1000010&&fl.Kfcode()<1000017&&fl.Kfcode()%2==0&&fl.IntSpin()==0; }

  inline bool IsSquark(const ATOOLS::Flavour &fl) 
  { return fl.Strong()&&(fl.StrongCharge()==3 || fl.StrongCharge()==-3)&&fl.IntSpin()==0&&!fl.Majorana(); }
  inline bool IsGluino(const ATOOLS::Flavour &fl) 
  { return fl.Kfcode()==kf_Gluino; }

  inline bool IsIno(const ATOOLS::Flavour &fl)    
  { return IsGluino(fl)||IsNeutralino(fl)||IsChargino(fl); }

  inline int LeptonNumber(const ATOOLS::Flavour &fl) 
  { if (fl.IsLepton()||IsSlepton(fl)||IsSneutrino(fl)) return fl.IsAnti()?-1:1; return 0; }
  inline double BaryonNumber(const ATOOLS::Flavour &fl) 
  { if (fl.IsQuark()||IsSquark(fl)) return fl.IsAnti()?-1./3.:1./3.; return 0.; 
    if (abs(fl.StrongCharge())==3) return 1./double(fl.StrongCharge()); }

}

#endif
